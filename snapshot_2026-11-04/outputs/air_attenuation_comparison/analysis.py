#!/usr/bin/env python3
"""Finite geometry versus effective attenuation; no Geant4 inputs.

Usage: python analysis.py INPUT.csv [--output-dir DIRECTORY]
Counts are modelled as Poisson, each distance measurement as an independent
normal constraint, and the common offset as one shared normal constraint.
Offset and distance uncertainties are interpreted as one standard deviation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, brentq
from scipy.special import xlogy
from scipy.stats import chi2
from numpy.polynomial.legendre import leggauss

RS = .004875
SIDE = .010
DELTA0, DELTA_SD = .018, .0015

def efficiency(r, radial=24, angular=96, derivative=False):
    r = np.atleast_1d(r).astype(float)
    if np.any(r <= 0):
        raise ValueError('Source–detector separation must be positive')
    z, w = leggauss(radial)
    rho = RS * np.sqrt((z + 1) / 2)
    phi = 2*np.pi*(np.arange(angular)+.5)/angular
    xs = rho[:, None]*np.cos(phi)[None, :]
    ys = rho[:, None]*np.sin(phi)[None, :]
    x1, x2 = -SIDE/2-xs, SIDE/2-xs
    y1, y2 = -SIDE/2-ys, SIDE/2-ys
    rr = r[:, None, None]
    def corner(a, b):
        s = np.sqrt(rr*rr+a*a+b*b)
        if derivative:
            return -a*b*(2*rr*rr+a*a+b*b)/(s*(rr*rr+a*a)*(rr*rr+b*b))
        return np.arctan2(a*b, rr*s)
    omega = corner(x2,y2)-corner(x2,y1)-corner(x1,y2)+corner(x1,y1)
    return np.sum(omega * (w[None,:,None]/2), axis=(1,2))/(angular*4*np.pi)

def verify_geometry():
    r = np.array([.0265, .05, .118, .818])
    e = efficiency(r)
    fine = efficiency(r, 48, 192)
    h = 1e-7
    deriv = efficiency(r, derivative=True)
    fd = (efficiency(r+h)-efficiency(r-h))/(2*h)
    # Independent quadrature over the detector plane as well as source disc.
    z, w = leggauss(24); phi = 2*np.pi*(np.arange(96)+.5)/96
    rho = RS*np.sqrt((z+1)/2)
    sx = (rho[:,None]*np.cos(phi)).ravel()
    sy = (rho[:,None]*np.sin(phi)).ravel()
    ws = np.repeat(w/2/96, 96)
    q, v = leggauss(16)
    dx, dy = np.meshgrid(q*SIDE/2, q*SIDE/2, indexing='ij')
    wd = np.outer(v,v).ravel()*(SIDE/2)**2
    dsq = (dx.ravel()[None,:]-sx[:,None])**2+(dy.ravel()[None,:]-sy[:,None])**2
    direct = np.array([np.sum(ws[:,None]*wd[None,:]*rr/(rr*rr+dsq)**1.5)/(4*np.pi) for rr in r])
    checks = {'quadrature_max_relative_difference': float(np.max(abs(e/fine-1))),
              'detector_integral_max_relative_difference': float(np.max(abs(e/direct-1))),
              'derivative_max_relative_difference': float(np.max(abs(deriv/fd-1)))}
    if max(checks.values()) > 2e-6: raise AssertionError(checks)
    return checks

def read_data(path):
    df = pd.read_csv(path)
    for col in ['d','d_unc','n','t']:
        if col not in df: raise ValueError('Missing required column '+col)
    background = df['d'].astype(str).str.strip().str.lower().eq('background')
    omitted = df.loc[background].to_dict('records')
    data = df.loc[~background].copy()
    for c in ['d','d_unc','n','t']:
        data[c] = pd.to_numeric(data[c], errors='raise')
    if not np.isfinite(data[['d','d_unc','n','t']].to_numpy()).all(): raise ValueError('Nonfinite input')
    if ((data.d <= 0) | (data.d > 2) | (data.d_unc <= 0) | (data.d_unc >= data.d) | (data.n < 0) | (data.t <= 0)).any():
        raise ValueError('Implausible SI values or uncertainties; inspect input units')
    if not np.allclose(data.n, np.round(data.n)): raise ValueError('n must be integer event counts')
    return data.reset_index(drop=True), omitted

class Fit:
    def __init__(self, data, model, offset=None, mu_fixed=None):
        self.x = data.d.to_numpy(float); self.s = data.d_unc.to_numpy(float)
        self.n = data.n.to_numpy(float); self.t = data.t.to_numpy(float)
        self.model, self.offset, self.mu_fixed = model, offset, mu_fixed
        self.free_mu = model == 'M2' and mu_fixed is None
        self.joint = offset is None
        self.p = 1+int(self.free_mu)+int(self.joint)
        self.nobs = len(data)
    def unpack(self, p):
        idx=1; mu=0.
        if self.model == 'M2':
            if self.free_mu: mu=p[idx]; idx+=1
            else: mu=self.mu_fixed
        delta=self.offset
        if self.joint: delta=DELTA0+DELTA_SD*p[idx]; idx+=1
        return np.exp(p[0]), mu, delta, p[idx:]
    def calc(self, p, jac=False):
        C, mu, delta, z = self.unpack(p)
        x=self.x+self.s*z
        eps=efficiency(x+delta)
        lam=self.t*C*eps*np.exp(-mu*x)
        dev=np.maximum(2*(lam-self.n+xlogy(self.n,self.n/lam)),0)
        resid=np.sign(self.n-lam)*np.sqrt(dev)
        allres=np.r_[resid, z, (delta-DELTA0)/DELTA_SD] if self.joint else np.r_[resid,z]
        if not jac: return allres
        de=efficiency(x+delta, derivative=True)/eps
        deriv=np.empty_like(lam); near=dev<1e-8
        deriv[~near]=np.sign(self.n[~near]-lam[~near])*(1-self.n[~near]/lam[~near])/np.sqrt(dev[~near])
        deriv[near]=-1/np.sqrt(lam[near])
        factor=deriv*lam
        J=np.zeros((len(allres),len(p)))
        J[:self.nobs,0]=factor
        idx=1
        if self.free_mu: J[:self.nobs,idx]=-x*factor; idx+=1
        if self.joint: J[:self.nobs,idx]=DELTA_SD*de*factor; J[-1,idx]=1; idx+=1
        J[np.arange(self.nobs),idx+np.arange(self.nobs)]=factor*self.s*(de-mu)
        J[self.nobs+np.arange(self.nobs),idx+np.arange(self.nobs)]=1
        return J
    def run(self):
        e=efficiency(self.x+(DELTA0 if self.joint else self.offset))
        c0=self.n.sum()/np.sum(self.t*e)
        vals=[np.log(c0)]; lo=[-20]; hi=[35]
        if self.free_mu: vals += [.3]; lo += [0]; hi += [20]
        if self.joint: vals += [0.]; lo += [-8]; hi += [8]
        vals+=list(np.zeros(self.nobs)); lo+=list(np.maximum(-12,-self.x/self.s+1e-5)); hi+=list(np.full(self.nobs,12))
        result=least_squares(self.calc, vals, jac=lambda p:self.calc(p,True), bounds=(lo,hi),
                             xtol=1e-10,ftol=1e-10,gtol=1e-9,max_nfev=600)
        if not result.success: raise RuntimeError(result.message)
        self.result=result; self.D=float(np.dot(result.fun,result.fun))
        self.C,self.mu,self.delta,self.z=self.unpack(result.x)
        self.cov=np.linalg.pinv(result.jac.T@result.jac, hermitian=True)
        self.C_se=self.C*np.sqrt(self.cov[0,0])
        self.mu_se=np.sqrt(self.cov[1,1]) if self.free_mu else np.nan
        return self
    def profile_interval(self):
        if not self.free_mu: return (np.nan,np.nan)
        target=self.D+1
        def f(mu): return Fit(data_global,'M2',self.offset,mu_fixed=mu).run().D-target
        left=0. if self.mu <= 1e-8 or f(0)<=0 else brentq(f,0,self.mu,xtol=2e-5)
        right=max(self.mu+.1, self.mu+2*self.mu_se)
        while right<20 and f(right)<0: right=min(20,right*1.8+.1)
        right=brentq(f,self.mu,right,xtol=2e-5) if f(right)>=0 else np.nan
        return left,right
    def predict(self, d):
        return self.C*efficiency(np.asarray(d)+self.delta)*np.exp(-self.mu*np.asarray(d))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_csv',type=Path); parser.add_argument('--output-dir',type=Path,default=Path(__file__).parent)
    args=parser.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    global data_global
    data_global, omitted=read_data(args.input_csv)
    checks=verify_geometry(); N=len(data_global)
    fits={}
    for offset in [.0165,.018,.0195,None]:
        for model in ['M1','M2']:
            fit=Fit(data_global,model,offset).run()
            fit.interval=fit.profile_interval() if model=='M2' and offset in [.018,None] else (np.nan,np.nan)
            fits[(offset,model)]=fit
    rows=[]
    for offset in [.0165,.018,.0195,None]:
        d0=fits[(offset,'M1')].D
        for model in ['M1','M2']:
            f=fits[(offset,model)]; k=1+(model=='M2'); dof=N-k
            rows.append(dict(model=model,offset_mode='constrained common offset' if offset is None else 'fixed',
                offset_nominal_mm=18 if offset is None else offset*1000,offset_fit_mm=f.delta*1000,
                offset_constraint_sigma_mm=1.5 if offset is None else 0,n_observations=N,C_per_s=f.C,
                C_formal_se_per_s=f.C_se,mu_per_m=f.mu if model=='M2' else np.nan,
                mu_formal_se_per_m=f.mu_se,mu_profile_low_per_m=f.interval[0],mu_profile_high_per_m=f.interval[1],
                augmented_deviance=f.D,degrees_freedom=dof,reduced_deviance=f.D/dof,
                nominal_gof_p=chi2.sf(f.D,dof),delta_deviance_vs_M1=f.D-d0,
                delta_AIC_vs_M1=f.D-d0+2*(k-1),delta_BIC_vs_M1=f.D-d0+np.log(N)*(k-1),
                offset_pull_sigma=(f.delta-DELTA0)/DELTA_SD if offset is None else np.nan))
    table=pd.DataFrame(rows); table.to_csv(args.output_dir/'model_comparison.csv',index=False)
    data_global.to_csv(args.output_dir/'data_used.csv',index=False)
    d=data_global.d.to_numpy(float); sd=data_global.d_unc.to_numpy(float)
    rate=data_global.n.to_numpy(float)/data_global.t.to_numpy(float)
    erate=np.sqrt(data_global.n.to_numpy(float))/data_global.t.to_numpy(float)
    grid=np.linspace(max(.001,d.min()*.6),d.max()*1.03,500)
    fig,(ax,axr)=plt.subplots(2,1,figsize=(7.6,6.6),sharex=True,gridspec_kw={'height_ratios':[2.2,1],'hspace':.08},layout=None)
    # Display propagation follows u=R d². These vertical and horizontal
    # uncertainties are correlated: the fit uses the original counts and d.
    uerr=np.sqrt((erate*d*d)**2+(2*rate*d*sd)**2)
    ax.errorbar(d*100,rate*d*d,yerr=uerr,xerr=sd*100,fmt='o',ms=3.8,lw=.8,capsize=2,color='black',label='Measurements (propagated 1σ)')
    for model,color,style in [('M1','#777777','--'),('M2','#075c95','-')]:
        f=fits[(.018,model)]; ax.plot(grid*100,f.predict(grid)*grid**2,style,color=color,lw=2,label={'M1':'M1: finite geometry','M2':'M2: geometry + effective attenuation'}[model])
    res={'d_m':d,'n':data_global.n,'t_s':data_global.t}
    for model,color,marker in [('M1','#777777','s'),('M2','#075c95','o')]:
        f=fits[(.018,model)]; pred=f.predict(d); slope=pred*(efficiency(d+.018,derivative=True)/efficiency(d+.018)-f.mu)
        denom=np.sqrt(pred/data_global.t.to_numpy(float)+(slope*sd)**2)
        standardized=(rate-pred)/denom
        res['standardized_'+model]=standardized
        axr.scatter(d*100,standardized,s=15,color=color,marker=marker,label=model,zorder=3)
    pd.DataFrame(res).to_csv(args.output_dir/'residuals.csv',index=False)
    axr.axhline(0,c='black',lw=.8); axr.axhline(2,c='#bbbbbb',lw=.7,ls=':'); axr.axhline(-2,c='#bbbbbb',lw=.7,ls=':')
    ax.set_ylabel(r'$R d^2$  ($\mathrm{m^2\,s^{-1}}$)'); axr.set_ylabel('Standardized\nresidual'); axr.set_xlabel('Rail separation $d$ (cm)')
    ax.legend(frameon=False,fontsize=9,loc='best'); axr.legend(frameon=False,fontsize=8,ncol=2)
    ax.text(.015,.95,'(a)',transform=ax.transAxes,va='top',weight='bold'); axr.text(.015,.92,'(b)',transform=axr.transAxes,va='top',weight='bold')
    for axis in (ax,axr): axis.spines[['top','right']].set_visible(False); axis.tick_params(direction='out')
    fig.subplots_adjust(left=.14,right=.97,top=.97,bottom=.1)
    for ext in ['png','pdf']: fig.savefig(args.output_dir/('model_comparison.'+ext),dpi=300)
    plt.close(fig)
    metadata={'input':str(args.input_csv.resolve()),'sha256':hashlib.sha256(args.input_csv.read_bytes()).hexdigest(),
        'rows_used':N,'distance_range_m':[float(d.min()),float(d.max())], 'omitted_background_rows':omitted,
        'geometry_checks':checks,'notes':'Poisson counts; independent Gaussian distance errors; shared Gaussian offset in constrained fits. No background, dead time or Geant4 correction. Figure uses fixed 18 mm offset. Residual denominator includes counting and independent distance errors, not shared offset or fitted-parameter covariance. Deviance p values and information criteria are approximate; systematic error not included.'}
    (args.output_dir/'metadata.json').write_text(json.dumps(metadata,indent=2,default=str)+'\n')
    print(table.to_string(index=False)); print(json.dumps(metadata,indent=2,default=str))

if __name__=='__main__': main()
