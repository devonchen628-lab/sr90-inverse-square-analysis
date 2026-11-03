#!/usr/bin/env python3
"""Common-response Poisson fits for ideal, finite-geometry and attenuation models."""
import numpy as np
from scipy.optimize import least_squares, brentq
from scipy.special import xlogy
from geometry import efficiency, verify_geometry, read_data, DELTA0, DELTA_SD

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
        if name not in ('M0', 'M1', 'M2', 'P1', 'P2'): raise ValueError(name)
        self.data=data; self.name=name; self.tau=tau
        self.with_bg=True; self.atten=name in ('M2','P2')
        self.free_mu=self.atten and mu_fixed is None; self.mu_fixed=mu_fixed
        self.offset=offset if name != 'M0' else 0.; self.joint=name != 'M0' and offset is None
        self.d=data.d.to_numpy(float); self.s=data.d_unc.to_numpy(float)
        self.n=data.n.to_numpy(float); self.t=data.t.to_numpy(float); self.N=len(data)
        # One background parameter is shared by the scan and control where enabled.
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
        e=(1/x**2) if self.name=='M0' else ((.010*.010/(4*np.pi*(x+delta)**2)) if self.name in ('P1','P2') else efficiency(x+delta))
        S=C*e*np.exp(-mu*x)
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
        de=np.zeros_like(x) if self.name=='M0' else ((-2/(x+delta)) if self.name in ('P1','P2') else efficiency(x+delta,derivative=True)/e)
        J[:self.N,0]=fac*S; i=1
        if self.free_mu: J[:self.N,i]=-x*fac*S; i+=1
        if self.with_bg:
            J[:self.N,i]=fac*beta/(1-self.tau*beta)**2
            J[-1,i]=bder[0]*TBG*beta; i+=1
        if self.joint:
            J[:self.N,i]=fac*S*de*DELTA_SD
            J[2*self.N,i]=1.; i+=1
        dlogdx=(-2/x) if self.name=='M0' else (de-mu)
        J[np.arange(self.N),i+np.arange(self.N)]=fac*S*dlogdx*self.s
        J[self.N+np.arange(self.N),i+np.arange(self.N)]=1.
        return J
    def run(self, initial_mu=.4, initial_scale=1.):
        delta=DELTA0 if self.joint else self.offset
        e0=(1/self.d**2) if self.name=='M0' else ((.010*.010/(4*np.pi*(self.d+delta)**2)) if self.name in ('P1','P2') else efficiency(self.d+delta))
        C0=self.n.sum()/np.sum(self.t*e0)
        v=[np.log(C0*initial_scale)]; lo=[-20.]; hi=[35.]
        if self.free_mu: v+=[initial_mu]; lo+=[0.]; hi+=[20.]
        if self.with_bg: v+=[np.log(NBG/TBG)]; lo+=[-20.]; hi+=[5.]
        if self.joint: v+=[0.]; lo+=[-8.]; hi+=[8.]
        # Same latent-distance domain for every model: physical x > 0, no sigma clipping.
        v+=list(np.zeros(self.N)); lo+=list(-self.d/self.s+1e-5); hi+=list(np.full(self.N,np.inf))
        opt=least_squares(self.calc,v,jac=lambda p:self.calc(p,True),bounds=(lo,hi),xtol=1e-10,ftol=1e-10,gtol=1e-9,max_nfev=700)
        if not opt.success: raise RuntimeError(opt.message)
        self.opt=opt; self.D=float(opt.fun@opt.fun)
        self.C,self.mu,self.beta,self.delta,self.z=self.unpack(opt.x)
        self.B=self.beta/(1-self.tau*self.beta) if self.with_bg else 0.
        self.cov=np.linalg.pinv(opt.jac.T@opt.jac,hermitian=True)
        self.Cse=self.C*np.sqrt(self.cov[0,0]); self.muse=np.sqrt(self.cov[1,1]) if self.free_mu else np.nan
        bi=1+int(self.free_mu)
        self.betase=self.beta*np.sqrt(self.cov[bi,bi])
        self.deltase=DELTA_SD*np.sqrt(self.cov[bi+1,bi+1]) if self.joint else np.nan
        self.mu_delta_corr=self.cov[1,bi+1]/np.sqrt(self.cov[1,1]*self.cov[bi+1,bi+1]) if self.joint and self.free_mu else np.nan
        return self
    def predict(self,d):
        d=np.asarray(d); e=(1/d**2) if self.name=='M0' else ((.010*.010/(4*np.pi*(d+self.delta)**2)) if self.name in ('P1','P2') else efficiency(d+self.delta))
        S=self.C*e*np.exp(-self.mu*d)
        return (S+self.B)/(1+self.tau*(S+self.B))
    def profile(self):
        target=self.D+1
        def f(mu): return Model(self.data,self.name,self.tau,self.offset,mu).run().D-target
        a=0. if f(0.)<=0 else brentq(f,0.,self.mu,xtol=2e-5)
        b=self.mu+max(.1,2*self.muse)
        while b<20 and f(b)<0: b=min(20.,b*1.8+.1)
        return a,brentq(f,self.mu,b,xtol=2e-5) if f(b)>=0 else np.nan

    def slope(self,d):
        d=np.asarray(d)
        if self.name=='M0': e=1/d**2; de=-2/d
        elif self.name in ('P1','P2'): e=.010*.010/(4*np.pi*(d+self.delta)**2); de=-2/(d+self.delta)
        else: e=efficiency(d+self.delta); de=efficiency(d+self.delta,derivative=True)/e
        S=self.C*e*np.exp(-self.mu*d)
        return S*(de-self.mu)/(1+self.tau*(S+self.B))**2
