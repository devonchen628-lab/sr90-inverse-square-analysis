#!/usr/bin/env python3
"""Draw the previously fitted full-range models, without refitting."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from report_style import apply,BLUE,RUST

P=Path(__file__).resolve().parent
ROOT=P.parent
sys.path.insert(0,str(ROOT/'air_attenuation_comparison'))
from analysis import efficiency

def main():
    apply()
    base=ROOT/'air_attenuation_comparison'
    data=pd.read_csv(base/'results_shared_background/data_used.csv')
    res=pd.read_csv(base/'results_shared_background/residuals.csv')
    models=pd.read_csv(base/'results_final_table/M0_M3_comparison.csv').set_index('model')
    assert len(data)==len(res)==32 and np.allclose(data.d,res.d_m)
    assert np.allclose(data.n,res.n) and np.allclose(data.t,res.t_s)
    d=data.d.to_numpy(float); sd=data.d_unc.to_numpy(float)
    n=data.n.to_numpy(float); t=data.t.to_numpy(float); rate=n/t
    u=rate*d**2
    uncertainty=np.sqrt((d*d*np.sqrt(n)/t)**2+(2*d*rate*sd)**2)
    # Previously saved Book2 propagation, independently checked.
    assert np.allclose(u,data.u,rtol=1e-9) and np.allclose(uncertainty,data.u_unc,rtol=1e-9)
    grid=np.linspace(.004,.875,600)
    def predict(name,x):
        r=models.loc[name]; delta=r.offset_fitted_mm/1000; tau=r.tau_us*1e-6
        mu=0. if name=='M2' else r.mu_per_m
        beta=r.background_observed_per_s; B=beta/(1-tau*beta)
        s=r.normalization*efficiency(x+delta)*np.exp(-mu*x)
        return (s+B)/(1+tau*(s+B))
    for name in ('M2','M3'):
        row=models.loc[name];delta=row.offset_fitted_mm/1000;tau=row.tau_us*1e-6
        mu=0. if name=='M2' else row.mu_per_m; beta=row.background_observed_per_s
        B=beta/(1-tau*beta);e=efficiency(d+delta)
        source=row.normalization*e*np.exp(-mu*d)
        slope=source*(efficiency(d+delta,derivative=True)/e-mu)/(1+tau*(source+B))**2
        pred=predict(name,d)
        standardized=(rate-pred)/np.sqrt(pred/t+(slope*sd)**2)
        assert np.allclose(standardized,res[f'{name}_standardized_residual'],rtol=2e-5,atol=2e-5)
    fig,(ax,bx)=plt.subplots(2,1,figsize=(3.5,4.0),sharex=True,
        gridspec_kw={'height_ratios':[2.15,1],'hspace':.06})
    ax.errorbar(d*100,u,yerr=uncertainty,fmt='o',ms=2.7,color='#111111',
                ecolor='#111111',elinewidth=.65,capsize=1.3,label='Data',zorder=5)
    ax.axhline(models.loc['M0','normalization'],color='#777777',ls=':',lw=1.25,label='M0: constant')
    ax.plot(grid*100,predict('M2',grid)*grid**2,color=BLUE,ls='--',lw=1.15,label='M2: finite geometry')
    ax.plot(grid*100,predict('M3',grid)*grid**2,color=RUST,lw=1.5,label='M3: geometry + attenuation')
    ax.set_ylabel(r'$R d^2$ ($\mathrm{m^2\,s^{-1}}$)')
    ax.set_ylim(0,14.5);ax.set_xlim(0,88)
    ax.text(.025,.94,'(a)',transform=ax.transAxes,va='top',weight='bold')
    handles,labels=ax.get_legend_handles_labels()
    by_label=dict(zip(labels,handles));order=['Data','M0: constant','M2: finite geometry','M3: geometry + attenuation']
    ax.legend([by_label[k] for k in order],order,frameon=False,loc='lower right',handlelength=2.1,labelspacing=.25,borderaxespad=.2)
    bx.axhline(0,color='#333333',lw=.7,zorder=0)
    for y in (-2,2): bx.axhline(y,color='#b0b0b0',ls=':',lw=.7,zorder=0)
    bx.plot(d*100,res.M2_standardized_residual,'s',ms=2.7,mfc='none',mec=BLUE,mew=.75,ls='none',label='M2')
    bx.plot(d*100,res.M3_standardized_residual,'o',ms=2.6,color=RUST,ls='none',label='M3')
    bx.set_ylabel('Standardized\nresidual');bx.set_xlabel('Rail separation $d$ (cm)')
    bx.set_ylim(-14.5,9);bx.set_yticks([-10,-5,0,5]);bx.legend(frameon=False,ncol=2,loc='upper right',handlelength=1)
    bx.text(.025,.92,'(b)',transform=bx.transAxes,va='top',weight='bold')
    fig.subplots_adjust(left=.18,right=.98,bottom=.11,top=.98)
    for ext in ('pdf','png'):fig.savefig(P/f'figure1_experimental_main.{ext}',dpi=400)
    plt.close(fig)
    caption=('Experimental model comparison. (a) Measured count rate multiplied by the square of the rail separation, '
      '$Rd^2$, with propagated counting and individual distance-reading uncertainties (one standard deviation). '
      'M0 is the fitted constant baseline; M2 includes finite source–detector geometry, background and a conditional '
      'non-paralyzable dead-time response; M3 additionally includes effective exponential attenuation. '
      '(b) Standardized residuals for M2 and M3, calculated using counting and independent distance-reading '
      'uncertainties, reveal remaining distance-dependent structure despite the improvement under M3. '
      'Shared offset, response-model uncertainty and fitted-parameter correlations are not included in these '
      'pointwise residual denominators. The fitted M2 offset is in tension with its independent constraint. '
      'The M0 constrained fit requires implausibly large distance pulls; '
      'its line is an exploratory reference rather than a validated fit.\n')
    (P/'figure1_caption.txt').write_text(caption)

if __name__=='__main__':main()
