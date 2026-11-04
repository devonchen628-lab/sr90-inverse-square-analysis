#!/usr/bin/env python3
"""Exploratory cross-check, existing weighted Geant4 batches versus M3.

Runs without Geant4. Full batch vectors are resampled together because source
histories were reused across distances. The diagonal-weight slope is only a
shape summary; the 32-point covariance has too few independent batches for an
ordinary full-covariance inverse.
"""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
ROOT=P.parent
sys.path.insert(0,str(ROOT/'air_attenuation_comparison'))
from analysis import efficiency

def load(medium,n):
    files=sorted((P/medium).glob('batch_*.csv'))
    assert len(files)==n
    values=[]
    for path in files:
        f=pd.read_csv(path); f=f[np.isclose(f.threshold_keV,10)].sort_values('distance_m')
        if medium=='air': assert np.allclose(f.active_distance_m-f.distance_m,.018)
        else: assert np.allclose(f.offset_m,.018)
        if not values: d=f.distance_m.to_numpy()
        assert np.allclose(d,f.distance_m) and np.all(f.probability>0)
        values.append(f.probability.to_numpy())
    return d,np.array(values)

d,a=load('air',8); dv,v=load('vacuum',12); assert np.allclose(d,dv)
data=pd.read_csv(ROOT/'air_attenuation_comparison/results_final_table/M0_M3_comparison.csv')
row=data[data.model=='M3'].iloc[0]; mu=row.mu_per_m; mulo=row.mu_profile_low; muhi=row.mu_profile_high
idx=np.flatnonzero(np.isclose(d,.35))[0]; ref=d[idx]
assert np.isclose(d[0],.01) and np.isclose(d[-1],.85)
pa,pv=a.mean(0),v.mean(0); ratio=pa/pv; norm=ratio/ratio[idx]
geom=efficiency(d+.018); gn=geom/geom[idx]; vn=pv/pv[idx]
rng=np.random.default_rng(639017); nboot=20000
ab=a[rng.integers(0,len(a),(nboot,len(a)))].mean(axis=1)
vb=v[rng.integers(0,len(v),(nboot,len(v)))].mean(axis=1)
assert np.all(ab>0) and np.all(vb>0)
vb_rel=vb/vb[:,idx,None]; rb=ab/vb; rb=rb/rb[:,idx,None]
err_log=np.std(np.log(ab/vb),axis=0,ddof=1)
weights=1/err_log**2; x=d; X=np.column_stack([np.ones(len(d)),x])
solver=np.linalg.solve(X.T@(weights[:,None]*X),X.T*weights)
coef=solver@np.log(ratio); bootstrap_coef=np.log(ab/vb)@solver.T
mu_mc=-coef[1]; mu_boot=-bootstrap_coef[:,1]
exp=np.exp(-mu*(d-ref)); exp_lo=np.exp(-mulo*(d-ref)); exp_hi=np.exp(-muhi*(d-ref))
vac_low,vac_high=np.percentile(vb_rel,[16,84],axis=0)
rat_low,rat_high=np.percentile(rb,[16,84],axis=0)
frame=pd.DataFrame(dict(d_m=d,air_probability=pa,vacuum_probability=pv,air_over_vacuum=ratio,
     geant_vacuum_relative=vn,geant_vacuum_rel_lo=vac_low,geant_vacuum_rel_hi=vac_high,
     geometry_relative_18mm=gn,air_vacuum_relative=norm,air_vacuum_rel_lo=rat_low,air_vacuum_rel_hi=rat_high,
     experimental_model_relative=exp))
frame.to_csv(P/'comparison_ref35.csv',index=False)
j=np.flatnonzero(np.isclose(d,.85))[0]
out={'reference_distance_m':ref,'distances_m':[float(d[0]),float(d[-1])],'points':len(d),
     'mu_exp_per_m':mu,'mu_exp_profile_interval_per_m':[mulo,muhi],
     'mu_MC_per_m':float(mu_mc),'mu_MC_bootstrap_16_84_per_m':np.percentile(mu_boot,[16,84]).tolist(),
     'mu_MC_bootstrap_std_per_m':float(mu_boot.std(ddof=1)),
     'MC_airvac_85_over_35':float(norm[j]),'MC_airvac_85_over_35_16_84':np.percentile(rb[:,j],[16,84]).tolist(),
     'exp_model_85_over_35':float(exp[j]),'vacuum_over_geo_normalized_85':float(vn[j]/gn[j]),
     'vacuum_over_geo_normalized_85_16_84':np.percentile(vb_rel[:,j]/gn[j],[16,84]).tolist(),
     'air_batches':len(a),'vacuum_batches':len(v),'bootstrap_draws':nboot,
     'slope_fit':'weighted straight-line fit to log air/vac ratio over all 32 distances; fixed diagonal bootstrap variance weights; entire source-library batches resampled to estimate uncertainty. Correlated distance estimates, not a full GLS fit; CI conditional on model and limited batches.',
     'caveats':['Geant4 emitter diameter 22 mm but analytical model emitter diameter 9.75 mm; aperture diameter 9.75 mm, not necessarily active diameter.','Geant4 and analytical geometry use 18 mm; experimental M3 profiles a common offset around 18 ± 1.5 mm and finds 19.04 mm.','Geant4 includes assumed 10 keV threshold, source housing/foil/acrylic and detector response, not present in analytical efficiency.','Independent air and vacuum batches; source histories reused across distances; no inverse of rank-limited 32x32 MC covariance.','MC intervals exclude geometry, materials, spectrum and threshold systematics; experimental profile interval excludes response model uncertainty.','Unnormalised Geant4 air/vac is a model-dependent response ratio; not measured transmission, nor proof of physical cause.']}
(P/'summary_ref35.json').write_text(json.dumps(out,indent=2)+'\n')
fig,(ax,bx)=plt.subplots(2,1,figsize=(8,7.2),sharex=True,layout='constrained')
ax.fill_between(d*100,vac_low/gn,vac_high/gn,color='#5481a6',alpha=.20,label='Geant4 vacuum: batch 16–84%')
ax.plot(d*100,vn/gn,'o',ms=3,color='#285b83',label='Vacuum / analytical geometry')
ax.axhline(1,ls='--',lw=1.5,color='#333333',label='Identical relative shape')
ax.set_ylabel('Relative-shape ratio\n(normalised at 35 cm)');ax.set_title('(a) Vacuum response compared with geometry',loc='left',fontsize=10);ax.legend(frameon=False,fontsize=8)
bx.fill_between(d*100,rat_low,rat_high,color='#5481a6',alpha=.20,label='Geant4: batch 16–84%')
bx.plot(d*100,norm,'o',ms=3,color='#285b83',label='Geant4 air / vacuum')
bx.plot(d*100,exp,lw=1.8,color='#ba4a32',label='Experimental M3 effective exponential')
bx.fill_between(d*100,np.minimum(exp_lo,exp_hi),np.maximum(exp_lo,exp_hi),color='#ba4a32',alpha=.16,label='M3 formal profile interval')
bx.axhline(1,color='#888888',lw=.7);bx.set_xlabel('Rail distance $d$ (cm)');bx.set_ylabel('Response / value at 35 cm');bx.set_title('(b) Air/vacuum ratio and experimental shape',loc='left',fontsize=10);bx.legend(frameon=False,fontsize=8)
for axis in (ax,bx): axis.spines[['top','right']].set_visible(False);axis.grid(alpha=.15)
for ext in ('png','pdf'):fig.savefig(P/f'crosscheck_ref35.{ext}',dpi=300)
print(json.dumps(out,indent=2))
