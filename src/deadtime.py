from pathlib import Path
import os,csv,json,re
work=Path(__file__).parent
os.environ['MPLCONFIGDIR']=str(work/'mplcache')
import numpy as np
from scipy.optimize import minimize
from scipy.stats import exponnorm,chi2,norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

src=work.parent/'data/raw/deadtime.csv'
txt=src.read_text(encoding='utf-8-sig')
a=np.array(list(csv.reader(txt.split('BinIndex,Count',1)[1].strip().splitlines())),int)
n=a[:,1]
N=n.sum()
dt=31.25/500 # microseconds; settings confirmed by user
edges=np.arange(500)*dt # 499 ordinary bins; last bin is overflow
centers=(edges[:-1]+edges[1:])/2

def probs(p,model='soft'):
    mu,ls,ll=p
    sigma=np.exp(ls);lam=np.exp(ll)
    if model=='soft':
        cdf=exponnorm.cdf(edges,K=1/(lam*sigma),loc=mu,scale=sigma)
    else:
        cdf=-np.expm1(-lam*np.maximum(edges-mu,0))
    pp=np.r_[np.diff(cdf),1-cdf[-1]]
    return np.maximum(pp,1e-300)

def objective(p,model='soft'):
    pp=probs(p,model)
    return -np.sum(n*np.log(pp))

fit=minimize(objective,[1.95,np.log(.28),np.log(.005)],method='Nelder-Mead',options={'maxiter':5000,'xatol':1e-10,'fatol':1e-7})
p=fit.x
mu,sigma,lam=p[0],np.exp(p[1]),np.exp(p[2])
pp=probs(p);expected=N*pp
positive=n>0
dev=2*np.sum(n[positive]*np.log(n[positive]/expected[positive]))
# Local numerical curvature; statistical errors conditional on this descriptive model.
h=np.array([.0005,.001,.001])
H=np.zeros((3,3));f0=objective(p)
for i in range(3):
    ei=np.eye(3)[i]*h[i]
    H[i,i]=(objective(p+ei)-2*f0+objective(p-ei))/h[i]**2
    for j in range(i):
        ej=np.eye(3)[j]*h[j]
        H[i,j]=H[j,i]=(objective(p+ei+ej)-objective(p+ei-ej)-objective(p-ei+ej)+objective(p-ei-ej))/(4*h[i]*h[j])
cov=np.linalg.inv(H)
se=np.sqrt(np.diag(cov))

# Compare a genuinely sharp cutoff with the same total-count normalization.
def sharp_obj(q):return objective([q[0],np.log(.001),q[1]],'sharp')
sharp=minimize(sharp_obj,[1.15,np.log(.005)],method='Nelder-Mead',options={'maxiter':3000})
esharp=N*probs([sharp.x[0],np.log(.001),sharp.x[1]],'sharp')
sharpdev=2*np.sum(n[positive]*np.log(n[positive]/esharp[positive]))

# Shape sensitivity: permit an unconstrained tail amplitude in the displayed range.
def expected_free(q):
    m,ls,ll,la=q;s=np.exp(ls);l=np.exp(ll);amp=np.exp(la)
    cdf=exponnorm.cdf(edges,K=1/(l*s),loc=m,scale=s)
    return np.maximum(amp*np.diff(cdf),1e-300)
def free_obj(q):
    e=expected_free(q)
    return np.sum(e-n[:-1]*np.log(e))
free=minimize(free_obj,[*p,np.log(N)],method='Nelder-Mead',options={'maxiter':5000,'xatol':1e-9})

elapsed=int(re.search(r'AccumulatedElapsedMicroseconds:\s*(\d+)',txt).group(1))/1e6
summary={
 'source':'data/raw/deadtime.csv','confirmed_display_range_us':31.25,'bins':500,'bin_width_us':dt,
 'total_counts':int(N),'elapsed_s':elapsed,'mean_recorded_rate_cps':N/elapsed,
 'overflow_counts':int(n[-1]),'overflow_fraction':float(n[-1]/N),
 'first_nonzero_bin':int(np.flatnonzero(n)[0]),
 'first_nonzero_bin_interval_us':[float(edges[np.flatnonzero(n)[0]]),float(edges[np.flatnonzero(n)[0]+1])],
 'model':'Descriptive soft recovery: Gaussian-distributed delay plus exponential interval; includes overflow bin.',
 'mu_delay_us':float(mu),'mu_statistical_se_us':float(se[0]),
 'sigma_recovery_us':float(sigma),'sigma_statistical_se_us':float(sigma*se[1]),
 'exponential_rate_per_us':float(lam),'exponential_rate_statistical_se_per_us':float(lam*se[2]),
 'effective_10_50_90_percent_recovery_us':{str(q):float(mu+lam*sigma*sigma+norm.ppf(q)*sigma) for q in [.1,.5,.9,.95]},
 'multinomial_deviance':float(dev),'nominal_dof':497,
 'nominal_chi2_p_NOT_CALIBRATED_LOW_EXPECTED_EDGE_BINS':float(chi2.sf(dev,497)),
 'sharp_cutoff_us':float(sharp.x[0]),'sharp_cutoff_deviance':float(sharpdev),
 'mu_finite_range_free_amplitude_us':float(free.x[0]),
 'free_amplitude_tail_rate_per_us':float(np.exp(free.x[2])),
 'statistical_uncertainty_caveat':'Conditional on recovery model, known time-axis mapping and overflow convention; not a full physical dead-time uncertainty.',
}
(work.parent/'results/tables/deadtime_fit.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))

fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
for ax in axs:
    ax.errorbar(centers,n[:-1],yerr=np.sqrt(n[:-1]),fmt='.',ms=3,capsize=0,color='#285780',label='Measured intervals')
    ax.plot(centers,expected[:-1],color='#c16c29',lw=2,label='Soft-recovery fit')
    ax.grid(alpha=.18)
    ax.spines[['top','right']].set_visible(False)
    ax.set(xlabel='Interval between recorded events (µs)',ylabel='Counts per 62.5 ns bin')
axs[0].set(xlim=(0,4),ylim=(0,130),title='Recovery edge: not an abrupt cutoff')
axs[0].axvline(mu+lam*sigma*sigma,color='#6c3a79',ls='--',lw=1.2,label=f'50% recovery: {mu+lam*sigma*sigma:.3f} µs')
axs[0].axvspan(mu+lam*sigma*sigma+norm.ppf(.1)*sigma,mu+lam*sigma*sigma+norm.ppf(.9)*sigma,color='#6c3a79',alpha=.10)
axs[0].legend(frameon=False,fontsize=8,loc='upper left')
axs[1].set(xlim=(0,31.1875),ylim=(0,135),title='Finite interval range (overflow excluded from plot)')
axs[1].legend(frameon=False,fontsize=9,loc='lower right')
fig.suptitle('Dead-time diagnostic',fontsize=14)
fig.supxlabel('The overflow bin is included in the fit. Effective recovery time depends on the stated model and definition.',fontsize=9)
fig.savefig(work.parent/'results/figures/deadtime_recovery.png',dpi=190)
