"""Generate comparable M0–M2 fits, component deviances and sensitivity studies."""
from pathlib import Path
import argparse, hashlib, json, platform
import numpy as np
import pandas as pd
from scipy.optimize._numdiff import approx_derivative
from scipy.stats import chi2
from models import Model, TAU, TAU_SD, read_data, verify_geometry

CORE=('M0','M1','M2')

def summarize(f, previous=np.nan):
    N=f.N; r=f.opt.fun
    count=float(r[:N]@r[:N]); dist=float(r[N:2*N]@r[N:2*N])
    off=float(r[2*N]**2) if f.joint else 0.
    bg=float(r[-1]**2); total=count+dist+off+bg
    assert np.isclose(total,f.D,atol=1e-8,rtol=1e-10)
    # 32 scan counts + 32 distance readings + background control [+ offset reading]
    # minus 32 latent distances + normalization + background [+ offset] [+ mu].
    nobs=2*N+1+int(f.joint); npar=len(f.opt.x); dof=nobs-npar
    assert dof==N-1-int(f.free_mu)
    return dict(model=f.name,normalization=f.C,normalization_se=f.Cse,
        normalization_unit='m^2 s^-1' if f.name=='M0' else 's^-1',
        mu_per_m=f.mu if f.atten else np.nan,mu_formal_se_per_m=f.muse,
        offset_fitted_mm=f.delta*1000 if f.name!='M0' else np.nan,
        offset_formal_se_mm=f.deltase*1000,mu_offset_correlation=f.mu_delta_corr,
        background_observed_per_s=f.beta,background_observed_se_per_s=f.betase,
        background_physical_per_s=f.B,background_physical_se_per_s=f.betase/(1-f.tau*f.beta)**2,
        tau_fixed_us=f.tau*1e6,scan_poisson_deviance=count,background_control_deviance=bg,
        distance_constraint_penalty=dist,offset_constraint_penalty=off,total_deviance=total,
        n_likelihood_observations=nobs,n_fitted_parameters=npar,dof=dof,
        reduced_total_deviance=total/dof,nominal_gof_p=chi2.sf(total,dof),
        improvement_vs_previous=previous-total,improvement_vs_previous_percent=100*(previous-total)/previous,
        max_absolute_distance_pull=float(np.max(abs(f.z))),optimizer_active_bounds=int(np.count_nonzero(f.opt.active_mask)),
        uncertainty_condition=('local curvature; tau fixed; shared offset, background and distances profiled' if f.joint else 'local curvature; tau fixed; offset absent; background and distances profiled' if f.name=='M0' else 'local curvature; tau and offset fixed; background and distances profiled')+'; no model-error rescaling')

def analyze(input_csv, output):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    data,excluded=read_data(input_csv)
    if len(data)!=32 or excluded: raise ValueError('Expected original 32 scan rows')
    geometry=verify_geometry(); fits={m:Model(data,m).run() for m in CORE}
    rows=[];previous=np.nan
    for name in CORE:
        row=summarize(fits[name],previous); rows.append(row); previous=fits[name].D
    lo,hi=fits['M2'].profile(); rows[-1].update(mu_profile_low_per_m=lo,mu_profile_high_per_m=hi)
    table=pd.DataFrame(rows);table.to_csv(output/'model_comparison.csv',index=False)
    point=[summarize(Model(data,name).run()) for name in ('P1','P2')]
    pd.DataFrame(point).to_csv(output/'point_source_sensitivity.csv',index=False)
    sensitivity=[]
    for off in (.0165,.018,.0195):
        for name in ('M1','M2'):
            f=Model(data,name,offset=off).run();row=summarize(f)
            row.update(offset_mode='fixed scenario',fixed_offset_mm=off*1000)
            sensitivity.append(row)
    sdf=pd.DataFrame(sensitivity);sdf.to_csv(output/'offset_sensitivity.csv',index=False)
    tau_rows=[]
    for off in (.0165,.018,.0195,None):
        for tau in (TAU-TAU_SD,TAU,TAU+TAU_SD):
            f=Model(data,'M2',tau=tau,offset=off).run();row=summarize(f)
            row.update(offset_mode='shared constrained' if off is None else 'fixed scenario',
                fixed_offset_mm=off*1000 if off is not None else np.nan)
            tau_rows.append(row)
    pd.DataFrame(tau_rows).to_csv(output/'offset_deadtime_sensitivity.csv',index=False)
    d=data.d.to_numpy();s=data.d_unc.to_numpy();n=data.n.to_numpy();t=data.t.to_numpy();rate=n/t
    residuals=pd.DataFrame(dict(d_m=d,d_unc_m=s,n=n,t_s=t,observed_rate_per_s=rate))
    latent=[]; predictions=pd.DataFrame({'d_m':np.linspace(.006,.875,500)})
    multistart={}; jaccheck={}
    for name,f in fits.items():
        pred=f.predict(d); denom=np.sqrt(pred/t+(f.slope(d)*s)**2)
        residuals[name+'_predicted_rate_per_s']=pred
        residuals[name+'_residual_sd_per_s']=denom
        residuals[name+'_standardized_residual']=(rate-pred)/denom
        predictions[name+'_predicted_rate_per_s']=f.predict(predictions.d_m.to_numpy())
        for i in range(len(data)):
            latent.append(dict(model=name,d_m=d[i],latent_d_m=d[i]+s[i]*f.z[i],distance_pull=f.z[i],
                latent_d_formal_se_m=s[i]*np.sqrt(f.cov[len(f.opt.x)-f.N+i,len(f.opt.x)-f.N+i]),
                signed_scan_deviance_residual=f.opt.fun[i]))
        # Independent starting conditions test optimization stability, not new model freedom.
        other=Model(data,name).run(initial_mu=.8,initial_scale=1.5)
        multistart[name]=abs(other.D-f.D)
        assert np.isclose(other.D,f.D,rtol=1e-8,atol=1e-5)
        analytic=f.calc(f.opt.x,True);numeric=approx_derivative(f.calc,f.opt.x,method='3-point',rel_step=1e-4)
        relative=np.linalg.norm(analytic-numeric,axis=0)/np.maximum(np.linalg.norm(analytic,axis=0),1e-5)
        jaccheck[name]=float(relative.max());assert relative.max()<3e-4,(name,relative.max())
        assert not np.any(f.opt.active_mask), (name,'active optimization bound')
    residuals.to_csv(output/'residuals.csv',index=False)
    pd.DataFrame(latent).to_csv(output/'latent_distances.csv',index=False)
    predictions.to_csv(output/'model_predictions.csv',index=False)
    # Historical checks apply ONLY to unchanged finite-geometry and attenuation models.
    assert abs(fits['M1'].D-1022.7902136755895)<1e-4
    assert abs(fits['M2'].D-216.60539762144833)<1e-4
    assert abs(fits['M2'].mu-.44077770817129375)<1e-6
    reference=Path(__file__).resolve().parents[1]/'archive/2026-11-03_M0-M3/docs/reference_offset_sensitivity.csv'
    ref=pd.read_csv(reference);ref=ref[ref.model.eq('M3')]
    actual=sdf[sdf.model.eq('M2')]
    for col in ['total_deviance','normalization','mu_per_m']:
        assert np.allclose(actual[col],ref[col],rtol=2e-6,atol=1e-6)
    # The newly fitted M0 MUST include the control and nonzero fixed dead time.
    assert fits['M0'].with_bg and fits['M0'].tau==TAU and not fits['M0'].joint
    assert abs(fits['M0'].D-2894.5347284403992)>1
    meta=dict(passed=True,python=platform.python_version(),raw_sha256=hashlib.sha256(Path(input_csv).read_bytes()).hexdigest(),
        geometry_checks=geometry,n_scan_rows=len(data),background_counts=8,background_exposure_s=60,
        tau_fixed_us=TAU*1e6,tau_sensitivity_sigma_us=TAU_SD*1e6,
        multi_start_total_deviance_difference=multistart,jacobian_max_column_relative_error=jaccheck,
        no_active_optimizer_bounds=True,offset_reference_check_passed=True,
        M0_refitted_with_response=True,total_deviances=table.total_deviance.tolist(),
        M2_mu_per_m=fits['M2'].mu,M2_mu_profile_interval=[lo,hi],
        figure_width_mm=178,figure_height_mm=130,figure_png_dpi=300)
    (output/'verification.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(table.to_string(index=False));print(json.dumps(meta,indent=2))
    return fits,data,sdf

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_csv',type=Path);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();analyze(args.input_csv,args.output_dir)
