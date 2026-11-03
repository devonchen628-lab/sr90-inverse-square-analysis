"""Three-panel report figure; designed at the final 178 mm two-column width."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

COLORS={'M0':'#666666','M1':'#0072B2','M2':'#D55E00'}

def make_figure(fits,data,sensitivity,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
        'legend.fontsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.65,
        'lines.linewidth':1.25,'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white',
        'axes.facecolor':'white','figure.facecolor':'white','mathtext.fontset':'dejavusans'})
    fig=plt.figure(figsize=(178/25.4,130/25.4))
    gs=fig.add_gridspec(2,2,height_ratios=[1.32,1],width_ratios=[1.7,1],
        left=.105,right=.985,bottom=.10,top=.985,hspace=.35,wspace=.38)
    a=fig.add_subplot(gs[0,:]);b=fig.add_subplot(gs[1,0]);c=fig.add_subplot(gs[1,1])
    d=data.d.to_numpy();sd=data.d_unc.to_numpy();n=data.n.to_numpy();t=data.t.to_numpy();R=n/t
    err=np.sqrt((np.sqrt(n)/t*d*d)**2+(2*R*d*sd)**2)
    grid=np.linspace(.006,.875,500)
    labels={'M0':'M0: ideal inverse square','M1':'M1: finite geometry','M2':'M2: geometry + attenuation'}
    for name in ('M0','M1','M2'):
        a.plot(grid*100,fits[name].predict(grid)*grid**2,color=COLORS[name],
            ls={'M0':':','M1':'--','M2':'-'}[name],label=labels[name],zorder=2)
    a.errorbar(d*100,R*d*d,xerr=sd*100,yerr=err,fmt='o',ms=2.7,mew=.5,
        mfc='white',mec='black',ecolor='black',color='black',capsize=1.4,elinewidth=.55,
        label=r'Data ($1\sigma$)',zorder=5)
    a.set(xlim=(0,88),ylim=(0,14.2),xlabel=r'Rail separation $d$ (cm)',
        ylabel=r'$Rd^2$ (m$^2$ s$^{-1}$)')
    a.legend(frameon=False,loc='lower right',ncol=2,columnspacing=1.6,handlelength=2.8)
    for name,marker in [('M1','s'),('M2','^')]:
        f=fits[name];pred=f.predict(d);res=(R-pred)/np.sqrt(pred/t+(f.slope(d)*sd)**2)
        b.scatter(d*100,res,s=12,marker=marker,facecolors='none',edgecolors=COLORS[name],
            linewidths=.8,label=name,zorder=3)
    b.axhline(0,color='black',lw=.7)
    for y in [-2,2]:b.axhline(y,color='.7',lw=.6,ls=':')
    b.set(xlim=(0,88),xlabel=r'Rail separation $d$ (cm)',ylabel='Standardised residual')
    b.legend(frameon=False,ncol=2,loc='lower left',handletextpad=.25)
    ss=sensitivity[sensitivity.model.eq('M2')].sort_values('fixed_offset_mm')
    c.errorbar(ss.fixed_offset_mm,ss.mu_per_m,yerr=ss.mu_formal_se_per_m,fmt='o',
        ms=3.7,color=COLORS['M2'],capsize=2,elinewidth=.8,label=r'Fixed $\delta$; conditional $1\sigma$')
    # No band or interpolated confidence interval: these are three refitted scenarios.
    c.set(xlim=(16.0,20.0),ylim=(.34,.48),xticks=[16.5,18,19.5],
        xlabel=r'Fixed common offset $\delta$ (mm)',ylabel=r'Fitted $\mu$ (m$^{-1}$)')
    c.text(.04,.96,'M2; fixed-offset fits\nBars: conditional statistical $1\\sigma$',transform=c.transAxes,va='top',fontsize=6.8,linespacing=1.6)
    for axis,label in [(a,'(a)'),(b,'(b)'),(c,'(c)')]:
        axis.spines[['top','right']].set_visible(False)
        axis.tick_params(direction='out',length=3,width=.6)
        axis.text(-.085 if axis is a else -.19,1.08 if axis is c else 1.0,label,transform=axis.transAxes,
            va='top',ha='left',fontsize=9,fontweight='bold')
    for ext in ('pdf','png'):fig.savefig(out/f'final_three_panel.{ext}',dpi=300)
    plt.close(fig)
