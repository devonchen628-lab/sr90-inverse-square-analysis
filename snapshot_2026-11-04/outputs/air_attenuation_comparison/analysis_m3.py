#!/usr/bin/env python3
"""Compare source models using raw Poisson counts and shared measurements.

Run: python analysis_m3.py data_used.csv --output-dir results_m3
Requires analysis.py in the same directory. No Geant4 is used.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, brentq
from scipy.special import xlogy
from scipy.stats import chi2
from analysis import efficiency, verify_geometry, read_data, DELTA0, DELTA_SD

TAU=1.97884e-6
TAU_SD=.02126e-6
NBG,TBG=8.,60.

def pres(n, lam):
    dev=np.maximum(2*(lam-n+xlogy(n,n/lam)),0.)
    res=np.sign(n-lam)*np.sqrt(dev)
    near=dev<1e-8
    der=np.empty_like(lam)
    der[~near]=np.sign(n[~near]-lam[~near])*(1-n[~near]/lam[~near])/np.sqrt(dev[~near])
    der[near]=-1/np.sqrt(lam[near])
    return res,der

class Model:
    def __init__(self,data,name,tau=TAU,offset=None,mu_fixed=None):
        self.data=data; self.name=name; self.tau=tau if name in ('M3','dead_only') else 0.
        self.with_bg=name in ('M3','background_only'); self.atten=name!='M1'
        self.free_mu=self.atten and mu_fixed is None; self.mu_fixed=mu_fixed
        self.offset=offset; self.joint=offset is None
        self.d=data.d.to_numpy(float); self.s=data.d_unc.to_numpy(float)
        self.n=data.n.to_numpy(float); self.t=data.t.to_numpy(float); self.N=len(data)
        # Background nuisance is common to all 32 points. For M1/M2/dead_only
        # it is irrelevant to the scan, so the independent control is saturated.
    def unpack(self,p):
        i=1; C=np.exp(p[0]); mu=0.
        if self.atten:
            if self.free_mu: mu=p[i]; i+=1
            else: mu=self.mu_fixed
        beta=NBG/TBG
        if self.with_bg: beta=np.exp(p[i]); i+=1
        delta=self.offset
        if self.joint: delta=DELTA0+DELTA_SD*p[i]; i+=1
        return C,mu,beta,delta,p[i:]
    def calc(self,p,jac=False):
        C,mu,beta,delta,z=self.unpack(p); x=self.d+self.s*z
        e=efficiency(x+delta); S=C*e*np.exp(-mu*x)
        B=beta/(1-self.tau*beta) if self.with_bg else 0.
        Q=S+B; den=1+self.tau*Q; R=Q/den; lam=self.t*R
        rc,dc=pres(self.n,lam)
        bres,bder=pres(np.array([NBG]),np.array([TBG*beta]))
        parts=[rc,z]
        if self.joint: parts.append(np.array([(delta-DELTA0)/DELTA_SD]))
        if self.with_bg: parts.append(bres)
        res=np.concatenate(parts)
        if not jac: return res
        J=np.zeros((len(res),len(p))); fac=dc*self.t/(den*den)
        de=efficiency(x+delta,derivative=True)/e
        J[:self.N,0]=fac*S; i=1
        if self.free_mu: J[:self.N,i]=-x*fac*S; i+=1
        if self.with_bg:
            J[:self.N,i]=fac*beta/(1-self.tau*beta)**2
            J[-1,i]=bder[0]*TBG*beta; i+=1
        if self.joint:
            J[:self.N,i]=fac*S*de*DELTA_SD
            J[2*self.N,i]=1.; i+=1
        J[np.arange(self.N),i+np.arange(self.N)]=fac*S*(de-mu)*self.s
        J[self.N+np.arange(self.N),i+np.arange(self.N)]=1.
        return J
    def run(self):
        delta=DELTA0 if self.joint else self.offset
        C0=self.n.sum()/np.sum(self.t*efficiency(self.d+delta))
        v=[np.log(C0)]; lo=[-20.]; hi=[35.]
        if self.free_mu: v+=[.4]; lo+=[0.]; hi+=[20.]
        if self.with_bg: v+=[np.log(NBG/TBG)]; lo+=[-20.]; hi+=[5.]
        if self.joint: v+=[0.]; lo+=[-8.]; hi+=[8.]
        v+=list(np.zeros(self.N)); lo+=list(np.maximum(-12.,-self.d/self.s+1e-5)); hi+=list(np.full(self.N,12.))
        opt=least_squares(self.calc,v,jac=lambda p:self.calc(p,True),bounds=(lo,hi),xtol=1e-10,ftol=1e-10,gtol=1e-9,max_nfev=700)
        if not opt.success: raise RuntimeError(opt.message)
        self.opt=opt; self.D=float(opt.fun@opt.fun)
        self.C,self.mu,self.beta,self.delta,self.z=self.unpack(opt.x)
        self.B=self.beta/(1-self.tau*self.beta) if self.with_bg else 0.
        self.cov=np.linalg.pinv(opt.jac.T@opt.jac,hermitian=True)
        self.Cse=self.C*np.sqrt(self.cov[0,0]); self.muse=np.sqrt(self.cov[1,1]) if self.free_mu else np.nan
        return self
    def predict(self,d):
        d=np.asarray(d); S=self.C*efficiency(d+self.delta)*np.exp(-self.mu*d)
        return (S+self.B)/(1+self.tau*(S+self.B))
    def profile(self):
        target=self.D+1
        def f(mu): return Model(self.data,self.name,self.tau,self.offset,mu).run().D-target
        a=0. if f(0.)<=0 else brentq(f,0.,self.mu,xtol=2e-5)
        b=self.mu+max(.1,2*self.muse)
        while b<20 and f(b)<0: b=min(20.,b*1.8+.1)
        return a,brentq(f,self.mu,b,xtol=2e-5) if f(b)>=0 else np.nan

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('input_csv',type=Path); p.add_argument('--output-dir',type=Path,default=Path(__file__).parent/'results_m3'); args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True); data,omitted=read_data(args.input_csv)
    if len(data)!=32: raise ValueError('Expected all 32 Task16_2 rows')
    checks=verify_geometry(); N=len(data)
    fits={name:Model(data,name).run() for name in ('M1','M2','M3','background_only','dead_only')}
    sensitivity=[]
    for offset in (.0165,.018,.0195,None):
        for tau in (TAU-TAU_SD,TAU,TAU+TAU_SD):
            f=Model(data,'M3',tau,offset).run()
            sensitivity.append(dict(offset_mode='constrained' if offset is None else 'fixed',fixed_offset_mm=offset*1000 if offset else np.nan,tau_us=tau*1e6,fitted_offset_mm=f.delta*1000,C_per_s=f.C,mu_per_m=f.mu,mu_formal_se_conditional_tau=f.muse,background_observed_per_s=f.beta,background_pre_deadtime_per_s=f.B,augmented_deviance=f.D))
    pd.DataFrame(sensitivity).to_csv(args.output_dir/'offset_deadtime_sensitivity.csv',index=False)
    intervals={name:fits[name].profile() for name in ('M2','M3')}
    rows=[]; base=fits['M1'].D
    for name in ('M3','M2','M1','background_only','dead_only'):
        f=fits[name]
        k=1 if name=='M1' else 2; dof=N-k
        rows.append(dict(model=name,n_distance_points=N,C_per_s=f.C,C_formal_se=f.Cse,mu_per_m=f.mu if name!='M1' else np.nan,mu_formal_se=f.muse,mu_profile_low=intervals.get(name,(np.nan,np.nan))[0],mu_profile_high=intervals.get(name,(np.nan,np.nan))[1],offset_fitted_mm=f.delta*1000,tau_fixed_us=f.tau*1e6,background_observed_rate_per_s=f.beta,background_pre_deadtime_rate_per_s=f.B if f.with_bg else np.nan,augmented_deviance=f.D,dof=dof,reduced_deviance=f.D/dof,nominal_gof_p=chi2.sf(f.D,dof),delta_deviance_vs_M1=f.D-base,delta_AIC_vs_M1=f.D-base+2*(k-1),delta_BIC_vs_M1=f.D-base+np.log(N)*(k-1)))
    pd.DataFrame(rows).to_csv(args.output_dir/'model_comparison.csv',index=False)
    data.to_csv(args.output_dir/'data_used.csv',index=False)
    d=data.d.to_numpy(float); sd=data.d_unc.to_numpy(float); t=data.t.to_numpy(float); n=data.n.to_numpy(float); R=n/t
    f=fits['M3']; S=f.C*efficiency(d+f.delta)*np.exp(-f.mu*d); B=f.B; tau=f.tau
    cases={'background_only':(S+B)/S-1,'dead_time_only':(S/(1+tau*S))/S-1,'background_plus_dead_time':((S+B)/(1+tau*(S+B)))/S-1}
    diagnostics=[]
    for key,val in cases.items():
        idx=np.argmax(np.abs(val)); diagnostics.append(dict(effect=key,max_absolute_relative_change_percent=abs(val[idx])*100,signed_change_percent=val[idx]*100,distance_cm=d[idx]*100,baseline='same fitted M3 source-only S; no refitting'))
    pd.DataFrame(diagnostics).to_csv(args.output_dir/'effect_sizes.csv',index=False)
    grid=np.linspace(.006,.875,500)
    fig,(ax,ar)=plt.subplots(2,1,figsize=(7.6,6.6),sharex=True,gridspec_kw={'height_ratios':[2.2,1],'hspace':.08})
    uerr=np.sqrt((np.sqrt(n)/t*d*d)**2+(2*R*d*sd)**2)
    ax.errorbar(d*100,R*d*d,yerr=uerr,xerr=sd*100,fmt='o',ms=3.6,capsize=2,lw=.7,color='black',label='Measurements (propagated 1σ)',zorder=4)
    colors={'M1':'#777777','M2':'#0877ae','M3':'#be4d23'}; styles={'M1':'--','M2':'-','M3':'-'}
    residuals={'d_m':d,'n':n,'t_s':t}
    for name in ('M3','M2','M1'):
        f=fits[name]; ax.plot(grid*100,f.predict(grid)*grid**2,color=colors[name],ls=styles[name],lw=2,zorder={'M1':1,'M2':2,'M3':3}[name],label={'M1':'M1: geometry','M2':'M2: geometry + attenuation','M3':'M3: M2 + background + dead time'}[name])
        pred=f.predict(d); source=f.C*efficiency(d+f.delta)*np.exp(-f.mu*d)
        derivative=source*(efficiency(d+f.delta,derivative=True)/efficiency(d+f.delta)-f.mu)/(1+f.tau*(source+f.B))**2
        res=(R-pred)/np.sqrt(pred/t+(derivative*sd)**2); residuals[name+'_standardized_residual']=res
        ar.scatter(d*100,res,s=17,color=colors[name],marker={'M1':'s','M2':'o','M3':'^'}[name],label=name,alpha=.85,zorder={'M1':1,'M2':2,'M3':3}[name])
    pd.DataFrame(residuals).to_csv(args.output_dir/'residuals.csv',index=False)
    ar.axhline(0,c='black',lw=.8)
    for y in (-2,2): ar.axhline(y,c='#bbbbbb',ls=':',lw=.8)
    ax.set_ylabel(r'$R d^2\ (\mathrm{m^2\,s^{-1}})$'); ar.set_ylabel('Standardized\nresidual'); ar.set_xlabel('Rail separation $d$ (cm)')
    ax.legend(frameon=False,fontsize=8.5,loc='best'); ar.legend(frameon=False,fontsize=8,ncol=3)
    ax.text(.02,.95,'(a)',transform=ax.transAxes,va='top',weight='bold'); ar.text(.02,.92,'(b)',transform=ar.transAxes,va='top',weight='bold')
    for a in (ax,ar): a.spines[['top','right']].set_visible(False)
    fig.subplots_adjust(left=.14,right=.97,top=.97,bottom=.1)
    for ext in ('png','pdf'): fig.savefig(args.output_dir/('comparison_m1_m2_m3.'+ext),dpi=300)
    plt.close(fig)
    meta=dict(input=str(args.input_csv.resolve()),sha256=hashlib.sha256(args.input_csv.read_bytes()).hexdigest(),N=N,background_observed_counts=8,background_time_s=60,tau_us=TAU*1e6,tau_formal_sigma_us=TAU_SD*1e6,geometry_checks=checks,omitted=omitted,notes='Counts and background control modelled as Poisson; shared background parameter in applicable fits; fixed tau sensitivity; shared Gaussian offset. M3 response is an unvalidated non-paralyzable approximation. Profile intervals include background and offset, not uncertainty in tau or unmodelled systematic effects. Delta AIC/BIC are heuristic due to per-point latent nuisance parameters and misspecification.')
    (args.output_dir/'metadata.json').write_text(json.dumps(meta,indent=2,default=str)+'\n')
    print(pd.DataFrame(rows).to_string(index=False)); print(pd.DataFrame(diagnostics).to_string(index=False)); print(json.dumps(checks))

if __name__=='__main__': main()
