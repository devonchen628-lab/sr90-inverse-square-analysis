#!/usr/bin/env python3
"""Reproduce far-field analysis of Task16_2 and existing 18 mm Geant4 batches.

Dependencies: NumPy, pandas, SciPy, Matplotlib. No Geant4 run is performed.
The existing analysis modules are imported read-only from the sibling folder.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.stats import chi2

P=Path(__file__).resolve().parent
O=P.parent
sys.path.insert(0,str(O/'air_attenuation_comparison'))
from analysis import efficiency,read_data,verify_geometry
from analysis_staged import Model,TAU

DATA=Path('/Users/chenqiyue/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Radlab/Task16_2.csv')
G4=O/'geant4_aperture975_offset18'
OFFSET=.018
NBOOT=20000
SEED=6032026
AREA=.010**2
TOL=[.005,.01,.02]

def geometry_deviation(s):
    return np.abs(efficiency(np.atleast_1d(s))/(AREA/(4*np.pi*np.asarray(s)**2))-1)

def load_batches(kind,expected):
    files=sorted((G4/kind).glob('batch_*.csv'))
    assert len(files)==expected
    arrays=[]; seeds=[]; nsource=[]
    for f in files:
        frame=pd.read_csv(f); frame=frame[np.isclose(frame.threshold_keV,10)].sort_values('distance_m')
        assert len(frame)==32 and (frame.probability>0).all()
        if not arrays: ds=frame.distance_m.to_numpy(float)
        assert np.allclose(ds,frame.distance_m)
        if kind=='air':
            assert np.allclose(frame.active_distance_m-frame.distance_m,OFFSET)
            assert (frame.aperture_diameter_mm==9.75).all() and (frame.emitter_diameter_mm==22).all()
        else: assert np.allclose(frame.offset_m,OFFSET)
        assert frame.seed.nunique()==1 and frame.nsource.nunique()==1
        seeds.append(int(frame.seed.iloc[0]));nsource.append(int(frame.nsource.iloc[0]));arrays.append(frame.probability.to_numpy(float))
    assert len(set(seeds))==len(seeds)
    return ds,np.array(arrays),seeds,nsource

def reg_cov(samples):
    # Bootstrap covariance from 8/12 source libraries has limited rank. The
    # shrinkage is a predeclared numerical regularisation, not additional data.
    c=np.atleast_2d(np.cov(samples,rowvar=False,ddof=1)); diag=np.diag(np.diag(c))
    return .8*c+.2*diag

def anchored_mc(d,a,v,boot_a,boot_v,cut,ref):
    mask=d>=cut-1e-10; ind=np.flatnonzero(mask); ri=np.flatnonzero(np.isclose(d,ref))[0]
    assert mask[ri]
    air=a.mean(axis=0);vac=v.mean(axis=0); T=air/vac
    tb=boot_a/boot_v; nr=T/T[ri]; nb=tb/tb[:,ri,None]
    jj=ind[ind!=ri]; x=d[jj]-ref
    z=np.log(nr[jj]); zb=np.log(nb[:,jj]); C=reg_cov(zb)
    # Covariance is regularised because bootstrap rank is bounded by number of
    # independent batches. Coefficients and bootstrap bands are conditional.
    inv=np.linalg.inv(C); w=inv@x/(x@inv@x)
    mu=-float(w@z); mb=-(zb@w)
    residual=z+mu*x
    score=float(residual@inv@residual)
    # Null bootstrap: preserve correlations, centre each vector on its observed
    # value, refit anchored slope and calculate same quadratic discrepancy.
    dz=zb-z; fitnoise=dz-(dz@w)[:,None]*x
    q=np.einsum('bi,ij,bj->b',fitnoise,inv,fitnoise)
    p=(1+np.sum(q>=score))/(len(q)+1)
    return dict(mu=mu,mu_boot=mb,norm=nr,norm_boot=nb,T=T,Tboot=tb,
        score=score,p_boot=p,indices=ind,cov=C)

def vacuum_check(d,v,vb,geom,cut,ref):
    mask=d>=cut-1e-10; ri=np.flatnonzero(np.isclose(d,ref))[0]
    jj=np.flatnonzero(mask & ~np.isclose(d,ref))
    Q=v.mean(axis=0)/geom; q=Q/Q[ri]
    qb=(vb/geom); qb=qb/qb[:,ri,None]
    z=np.log(q[jj]);zb=np.log(qb[:,jj]);C=reg_cov(zb);inv=np.linalg.inv(C)
    score=float(z@inv@z);dz=zb-z
    qnull=np.einsum('bi,ij,bj->b',dz,inv,dz)
    p=(1+np.sum(qnull>=score))/(len(qnull)+1)
    return dict(q=q,qboot=qb,score=score,p_boot=p,idx=jj)

def exp_fit(df,offset=None,profile=False):
    f=Model(df,'M3',tau=TAU,offset=offset).run(); lo,hi=f.profile() if profile else (np.nan,np.nan)
    n=len(df); countdev=float(np.sum(f.opt.fun[:n]**2));bg=float(f.opt.fun[-1]**2)
    dpen=float(np.sum(f.opt.fun[n:2*n]**2)); openalty=((f.delta-.018)/.0015)**2 if offset is None else 0.
    assert np.isclose(countdev+bg+dpen+openalty,f.D,rtol=1e-8,atol=1e-8)
    # Existing code uses exp(-mu*d). Reparameterise C_phys=C_d exp(mu*delta)
    # for exp[-mu*(d+delta)]; predictions and likelihood are identical.
    Cphysical=f.C*np.exp(f.mu*f.delta)
    dd=df.d.to_numpy(float)
    lhs=f.C*efficiency(dd+f.delta)*np.exp(-f.mu*dd)
    rhs=Cphysical*efficiency(dd+f.delta)*np.exp(-f.mu*(dd+f.delta))
    assert np.allclose(lhs,rhs,rtol=1e-12,atol=0.)
    return f,dict(n=n,offset_mode='constrained' if offset is None else 'fixed',
        nominal_or_fixed_offset_mm=18. if offset is None else 1000*offset,
        fitted_offset_mm=1000*f.delta,C_rail_convention_per_s=f.C,C_physical_convention_per_s=Cphysical,
        mu_per_m=f.mu,mu_formal_se=f.muse,mu_profile_low=lo,mu_profile_high=hi,
        background_observed_per_s=f.beta,background_pre_deadtime_per_s=f.B,tau_us=f.tau*1e6,
        scan_poisson_deviance=countdev,background_control_deviance=bg,distance_penalty=dpen,
        offset_penalty=openalty,total_deviance=f.D,dof=n-2,reduced_total_deviance=f.D/(n-2),
        nominal_gof_p=chi2.sf(f.D,n-2),max_distance_pull_sigma=float(np.max(np.abs(f.z))))

def main():
    verify_geometry(); df,omitted=read_data(DATA);assert len(df)==32 and not omitted
    d,a,aseeds,an=load_batches('air',8);dv,v,vseeds,vn=load_batches('vacuum',12)
    assert np.allclose(d,dv) and np.allclose(d,df.d)
    assert set(aseeds).isdisjoint(vseeds), 'Air/vac seeds not independent'
    rng=np.random.default_rng(SEED)
    ab=a[rng.integers(0,len(a),(NBOOT,len(a)))].mean(axis=1)
    vb=v[rng.integers(0,len(v),(NBOOT,len(v)))].mean(axis=1)
    assert (ab>0).all() and (vb>0).all()
    g=efficiency(d+OFFSET); far=AREA/(4*np.pi*(d+OFFSET)**2); frac=abs(g/far-1)
    # Check monotonic decline densely before interpreting each root as a
    # threshold holding over the rest of the measured range.
    dense=np.linspace(d.min()+OFFSET,d.max()+OFFSET,1200); fdev=geometry_deviation(dense)
    assert (np.diff(fdev)<1e-10).all()
    thresholds=[]
    for tol in TOL:
        root=brentq(lambda s: geometry_deviation(s)[0]-tol,dense[0],dense[-1])
        rail=root-OFFSET; first=d[d>=rail-1e-12][0]
        thresholds.append(dict(tolerance=tol,physical_threshold_m=root,rail_threshold_m=rail,
            first_included_rail_m=first,first_included_physical_m=first+OFFSET,n_points=int(np.sum(d>=first-1e-10))))
    pd.DataFrame(thresholds).to_csv(P/'farfield_cutoff_analysis.csv',index=False)
    primary=next(r for r in thresholds if r['tolerance']==.01)
    cut=primary['first_included_rail_m']; ref=.30 if cut<=.30 else d[d>=cut][0]
    selected=df[df.d>=cut-1e-10].reset_index(drop=True)
    ex,exrow=exp_fit(selected,profile=True);exrow['cutoff_rail_m']=cut
    pd.DataFrame([exrow]).to_csv(P/'experimental_farfield_fit.csv',index=False)
    mc=anchored_mc(d,a,v,ab,vb,cut,ref)
    vac=vacuum_check(d,v,vb,g,cut,ref)
    bootmed,bootlo,boothi=np.percentile(mc['mu_boot'],[50,16,84])
    mcrow=dict(n_selected=len(selected),reference_rail_m=ref,cutoff_rail_m=cut,
       mu_best_per_m=mc['mu'],mu_bootstrap_median_per_m=bootmed,mu_bootstrap_p16_per_m=bootlo,
       mu_bootstrap_p84_per_m=boothi,bootstrap_samples=NBOOT,correlated_regularized_shape_statistic=mc['score'],
       shape_bootstrap_p_exploratory=mc['p_boot'],vacuum_shape_statistic=vac['score'],
       vacuum_bootstrap_p_exploratory=vac['p_boot'],relative_difference_mc_minus_exp=mc['mu']-ex.mu)
    pd.DataFrame([mcrow]).to_csv(P/'geant4_farfield_fit.csv',index=False)
    # Fixed windows were chosen before examining agreement; same observations
    # are used for experimental and simulation fits at each cutoff.
    sensitivity=[]
    for c in sorted(set([cut,.20,.30,.40])):
        subset=df[df.d>=c-1e-10].reset_index(drop=True);rf=.30 if c<=.30 else float(subset.d.iloc[0])
        ef,er=exp_fit(subset);mf=anchored_mc(d,a,v,ab,vb,c,rf)
        q16,q50,q84=np.percentile(mf['mu_boot'],[16,50,84])
        sensitivity.append(dict(cutoff_rail_m=c,first_physical_separation_m=c+OFFSET,n=len(subset),
            reference_rail_m=rf,mu_experimental=ef.mu,mu_experimental_formal_se=ef.muse,
            experimental_total_deviance=ef.D,experimental_df=len(subset)-2,mu_mc=mf['mu'],
            mu_mc_boot_median=q50,mu_mc_boot_p16=q16,mu_mc_boot_p84=q84))
    pd.DataFrame(sensitivity).to_csv(P/'cutoff_sensitivity.csv',index=False)
    offsetrows=[]
    for off in (.0165,.018,.0195):
        ef,er=exp_fit(selected,offset=off);offsetrows.append(er)
    pd.DataFrame(offsetrows).to_csv(P/'offset_sensitivity.csv',index=False)
    xx=d-ref; experimental=np.exp(-ex.mu*xx)
    elo=np.exp(-exrow['mu_profile_low']*xx);ehi=np.exp(-exrow['mu_profile_high']*xx)
    mlo,mhi=np.percentile(mc['norm_boot'],[16,84],axis=0)
    vlo,vhi=np.percentile(vac['qboot'],[16,84],axis=0)
    # Pointwise difference uncertainty, conditional, with statistically
    # independent MC and experimental mu. Avoid claiming a simultaneous band.
    mu_draw=rng.normal(ex.mu,ex.muse,NBOOT)
    dmuboot=mc['mu_boot']-mu_draw
    mcrow['difference_boot_median_per_m']=float(np.median(dmuboot))
    mcrow['difference_boot_p16_per_m']=float(np.percentile(dmuboot,16))
    mcrow['difference_boot_p84_per_m']=float(np.percentile(dmuboot,84))
    # Re-save after adding the conditional independent-MC/experiment difference.
    pd.DataFrame([mcrow]).to_csv(P/'geant4_farfield_fit.csv',index=False)
    diffboot=mc['norm_boot']-np.exp(-mu_draw[:,None]*xx)
    difflo,diffhi=np.percentile(diffboot,[16,84],axis=0)
    dsel=np.flatnonzero(d>=cut-1e-10)
    rate=selected.n.to_numpy(float)/selected.t.to_numpy(float);pred=ex.predict(selected.d)
    ep=efficiency(selected.d+ex.delta);epd=efficiency(selected.d+ex.delta,derivative=True)
    source=ex.C*ep*np.exp(-ex.mu*selected.d)
    derivative=source*(epd/ep-ex.mu)/(1+ex.tau*(source+ex.B))**2
    residual=(rate-pred)/np.sqrt(pred/selected.t.to_numpy(float)+(derivative*selected.d_unc.to_numpy(float))**2)
    trans=pd.DataFrame(dict(d_rail_m=d,s_physical_nominal_m=d+OFFSET,selected=d>=cut-1e-10,
        air_probability=a.mean(axis=0),vacuum_probability=v.mean(axis=0),T_mc=mc['T'],
        T_mc_normalized=mc['norm'],T_mc_p16=mlo,T_mc_p84=mhi,
        T_experimental_m3=experimental,T_experimental_p16=np.minimum(elo,ehi),T_experimental_p84=np.maximum(elo,ehi),
        difference_mc_minus_exp=mc['norm']-experimental,difference_boot_p16=difflo,difference_boot_p84=diffhi,
        vacuum_geometry_ratio=vac['q'],vacuum_geometry_p16=vlo,vacuum_geometry_p84=vhi))
    trans.to_csv(P/'transmission_comparison.csv',index=False)
    pd.DataFrame(dict(d_rail_m=selected.d,counts=selected.n,time_s=selected.t,predicted_rate_per_s=pred,
                      observed_rate_per_s=rate,standardized_residual=residual)).to_csv(P/'experimental_residuals.csv',index=False)
    intervals=[]
    for a1,b1 in [(cut,.30),(.30,.55),(.55,.85),(.30,.85)]:
        if not np.any(np.isclose(d,a1)) or not np.any(np.isclose(d,b1)): continue
        i=np.flatnonzero(np.isclose(d,a1))[0];j=np.flatnonzero(np.isclose(d,b1))[0]
        mb=-np.log(mc['Tboot'][:,j]/mc['Tboot'][:,i])/(b1-a1)
        intervals.append(dict(d1_m=a1,d2_m=b1,mu_effective=-np.log(mc['T'][j]/mc['T'][i])/(b1-a1),
                              median=np.median(mb),p16=np.percentile(mb,16),p84=np.percentile(mb,84)))
    pd.DataFrame(intervals).to_csv(P/'interval_slopes.csv',index=False)
    i30,i55,i85=[np.flatnonzero(np.isclose(d,x))[0] for x in (.30,.55,.85)]
    slope1=-np.log(mc['Tboot'][:,i55]/mc['Tboot'][:,i30])/.25
    slope2=-np.log(mc['Tboot'][:,i85]/mc['Tboot'][:,i55])/.30
    slope_delta=slope2-slope1
    mcrow['late_minus_middle_slope_boot_median']=float(np.median(slope_delta))
    mcrow['late_minus_middle_slope_boot_p16']=float(np.percentile(slope_delta,16))
    mcrow['late_minus_middle_slope_boot_p84']=float(np.percentile(slope_delta,84))
    pd.DataFrame([mcrow]).to_csv(P/'geant4_farfield_fit.csv',index=False)
    # Publication figures: physically selected window and common reference.
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,(ax,bx)=plt.subplots(2,1,figsize=(6.3,6.1),sharex=True,layout='constrained')
    ax.semilogy(d*100,frac*100,'o-',ms=3,color='#285b83');ax.axhline(1,color='#c04f30',ls='--',label='1% tolerance')
    ax.axvline(cut*100,color='#333333',ls=':',label=f'First selected point: {cut*100:g} cm')
    ax.set_ylabel('Geometry deviation (%)');ax.legend(frameon=False,fontsize=8)
    bx.fill_between(d*100,vlo,vhi,color='#285b83',alpha=.2,label='MC bootstrap 16–84%')
    bx.plot(d*100,vac['q'],'o',ms=3,color='#285b83');bx.axhline(1,color='#444444',ls='--')
    bx.axvline(cut*100,color='#333333',ls=':');bx.set_ylabel('Vacuum / geometry\nrelative to 30 cm');bx.set_xlabel('Rail separation (cm)');bx.legend(frameon=False,fontsize=8)
    for ext in ('png','pdf'):fig.savefig(P/f'figure1_farfield_selection.{ext}',dpi=300)
    plt.close(fig)
    fig,(ax,bx)=plt.subplots(2,1,figsize=(6.3,6),sharex=True,gridspec_kw={'height_ratios':[2,1]},layout='constrained')
    i=dsel;ax.fill_between(d[i]*100,mlo[i],mhi[i],color='#28638f',alpha=.20,label='Geant4 16–84%')
    ax.plot(d[i]*100,mc['norm'][i],'o',ms=3,color='#28638f',label='Geant4 air / vacuum')
    ax.fill_between(d[i]*100,np.minimum(elo,ehi)[i],np.maximum(elo,ehi)[i],color='#bc4e2d',alpha=.16,label='Experimental M3 formal interval')
    ax.plot(d[i]*100,experimental[i],color='#bc4e2d',lw=1.8,label='Experimental M3 shape')
    ax.axvline(ref*100,color='#777777',ls=':',lw=.8);ax.axvline(cut*100,color='#777777',ls='--',lw=.8)
    ax.text(.02,.97,f'Geometry cutoff: {cut*100:g} cm',transform=ax.transAxes,va='top',fontsize=7.5)
    ax.set_ylabel('Relative response\n(normalised at 30 cm)');ax.legend(frameon=False,fontsize=7.5)
    bx.fill_between(d[i]*100,difflo[i],diffhi[i],color='#506b80',alpha=.2,label='Pointwise conditional 16–84%')
    bx.plot(d[i]*100,(mc['norm']-experimental)[i],'o',ms=3,color='#285b83')
    bx.axhline(0,color='#333333',lw=.8);bx.axvline(ref*100,color='#777777',ls=':',lw=.8)
    bx.set_ylabel('Geant4 − M3');bx.set_xlabel('Rail separation (cm)');bx.legend(frameon=False,fontsize=7.5)
    for ext in ('png','pdf'):fig.savefig(P/f'figure2_transmission.{ext}',dpi=300)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.3,3.6),layout='constrained');i=dsel
    ax.fill_between(d[i]*100,np.log(mlo[i]),np.log(mhi[i]),color='#28638f',alpha=.2,label='Geant4 16–84%')
    ax.plot(d[i]*100,np.log(mc['norm'][i]),'o',ms=3,color='#28638f',label='Geant4 air/vacuum')
    ax.plot(d[i]*100,-mc['mu']*(d[i]-ref),color='#a44e2f',label=f'Anchored exponential, μ={mc["mu"]:.2f} m$^{{-1}}$')
    ax.axvline(ref*100,color='#777777',ls=':',lw=.8);ax.set_xlabel('Rail separation (cm)');ax.set_ylabel('ln(relative response)');ax.legend(frameon=False,fontsize=8)
    for ext in ('png','pdf'):fig.savefig(P/f'figure3_log_diagnostic.{ext}',dpi=300)
    plt.close(fig)
    full=pd.read_csv(O/'air_attenuation_comparison/results_final_table/M0_M3_comparison.csv');fullmu=float(full.loc[full.model=='M3','mu_per_m'].iloc[0])
    j=np.flatnonzero(np.isclose(d,.85))[0]
    summary=dict(primary=primary,experimental=exrow,geant4=mcrow,offset_sensitivity=offsetrows,
        cutoff_sensitivity=sensitivity,intervals=intervals,full_range_mu=fullmu,
        vacuum_relative_at_85=vac['q'][j],vacuum_85_boot_p16=vlo[j],vacuum_85_boot_p84=vhi[j],
        mc_T85=mc['norm'][j],mc_T85_p16=mlo[j],mc_T85_p84=mhi[j],exp_T85=experimental[j],
        reference_m=ref,input_file=str(DATA),input_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),
        air_seeds=aseeds,vacuum_seeds=vseeds,air_nsource=an,vacuum_nsource=vn,
        geometry_exponent_note='Existing code uses exp(-mu*d); C_physical=C_rail*exp(mu*offset) gives identical model exp(-mu*(d+offset)).',
        caveat='Gaussian constraints and Poisson counts conditional; nonparalyzable process is not exactly Poisson; MC covariance shrinkage and small batch bootstrap make shape tests exploratory.')
    (P/'summary.json').write_text(json.dumps(summary,indent=2,default=float)+'\n')
    print(json.dumps(summary,indent=2,default=float))

if __name__=='__main__':main()
