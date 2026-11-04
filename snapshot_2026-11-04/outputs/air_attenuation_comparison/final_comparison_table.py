#!/usr/bin/env python3
"""Create M0–M3 table and fixed-offset sensitivity from staged models."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from analysis_staged import Model, TAU, read_data, verify_geometry

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_csv',type=Path)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).parent/'results_final_table')
    args=parser.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    data,excluded=read_data(args.input_csv)
    if len(data)!=32 or excluded: raise ValueError('Expected 32 measurements without control rows')
    verify_geometry()
    models=('M0','M1','M2','M3')
    fits={m:Model(data,m).run() for m in models}
    lo,hi=fits['M3'].profile()
    rows=[]
    for i,name in enumerate(models):
        f=fits[name]; count_dev=float(np.sum(f.opt.fun[:len(data)]**2))
        bg_dev=float(f.opt.fun[-1]**2) if f.with_bg else 0.
        dist_dev=float(np.sum(f.opt.fun[len(data):2*len(data)]**2))
        offset_dev=((f.delta-.018)/.0015)**2
        total=count_dev+bg_dev+dist_dev+offset_dev
        assert np.isclose(total,f.D,rtol=1e-8,atol=1e-8)
        k=2 if name=='M3' else 1; df=len(data)-k
        previous=fits[models[i-1]].D if i else np.nan
        rows.append(dict(model=name,normalization=f.C,normalization_se=f.Cse,
            normalization_unit='m^2/s' if name=='M0' else '1/s',mu_per_m=f.mu if f.atten else np.nan,
            mu_formal_se=f.muse,mu_profile_low=lo if name=='M3' else np.nan,
            mu_profile_high=hi if name=='M3' else np.nan,offset_fitted_mm=f.delta*1000,
            background_observed_per_s=f.beta if f.with_bg else np.nan,
            tau_us=f.tau*1e6,poisson_scan_deviance=count_dev,
            background_control_deviance=bg_dev,distance_constraint_penalty=dist_dev,
            offset_constraint_penalty=offset_dev,total_deviance=total,dof=df,
            reduced_total_deviance=total/df,
            improvement_vs_previous_percent=100*(previous-total)/previous if i else np.nan,
            improvement_vs_M0_percent=100*(fits['M0'].D-total)/fits['M0'].D))
    df=pd.DataFrame(rows); df.to_csv(args.output_dir/'M0_M3_comparison.csv',index=False)
    sensitivity=[]
    for off in (.0165,.018,.0195):
        for name in models:
            f=Model(data,name,offset=off).run()
            sensitivity.append(dict(model=name,fixed_offset_mm=off*1000,normalization=f.C,
                normalization_unit='m^2/s' if name=='M0' else '1/s',mu_per_m=f.mu if f.atten else np.nan,
                mu_formal_se=f.muse,background_observed_per_s=f.beta if f.with_bg else np.nan,
                total_deviance=f.D,dof=len(data)-(2 if name=='M3' else 1)))
    sdf=pd.DataFrame(sensitivity); sdf.to_csv(args.output_dir/'offset_sensitivity_M0_M3.csv',index=False)
    print(df.to_string(index=False)); print(sdf.to_string(index=False))

if __name__=='__main__': main()
