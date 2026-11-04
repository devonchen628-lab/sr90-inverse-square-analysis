from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
x=pd.read_csv(P/'data/pooled_measurements.csv'); S=json.loads((P/'data/analysis_summary.json').read_text()); z=np.load(P/'data/G2_model.npz')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
blue='#175e86'; orange='#b54b25'; gray='#59636e'
d=x.d_m.to_numpy(); c=x.raw_rate_cps.to_numpy(); sd=.0005/np.sqrt(3); rawerr=np.sqrt((d*d*np.sqrt(x.total_counts.to_numpy())/x.total_time_s.to_numpy())**2+(2*d*c*sd)**2)
fig,ax=plt.subplots(figsize=(3.4,2.45),layout='constrained'); ax.errorbar(d*100,x.u_recorded,yerr=rawerr,fmt='o',color=blue,markersize=3,elinewidth=.65,capsize=1,label='Measured product')
ref=float(x.loc[np.isclose(x.d_m,.30),'u_recorded'].iloc[0]);ax.axhline(ref,ls='--',color=gray,lw=.9,label='Constant, fixed at 30 cm')
ax.set(xlabel='Recorded rail distance $d$ (cm)',ylabel=r'$C_{\rm obs}d^2$ (m$^2$ s$^{-1}$)',xlim=(.85,100),ylim=(0,16));ax.set_xscale('log');ax.set_xticks([1,2.5,10,30,85],labels=['1','2.5','10','30','85']);ax.grid(alpha=.14);ax.legend(frameon=False,loc='lower right');
for e in ('png','pdf'):fig.savefig(P/'images'/f'observed_product.{e}',dpi=300)
plt.close(fig)
md=z['d_m']; prob=z['probability']; A=S['models']['G2']['A']; ix=np.array([np.argmin(abs(d-v)) for v in md]); r=d+.021; u=x.corrected_rate_cps.to_numpy()*r*r; uerr=np.sqrt((x.corrected_rate_se.to_numpy()*r*r)**2+(2*r*x.corrected_rate_cps.to_numpy()*sd)**2); mr=md+.021; model=A*prob*mr**2; merr=A*np.sqrt(np.diag(z['probability_covariance']))*mr**2; ratio=u[ix]/model; ratioerr=ratio*np.sqrt(np.maximum(0,np.diag(z['ratio_covariance_log'])))
fig,(ax,bx)=plt.subplots(2,1,figsize=(3.4,3.6),sharex=True,gridspec_kw={'height_ratios':[1.45,1]},layout='constrained')
ax.errorbar(d*100,u,yerr=uerr,fmt='o',markersize=2.6,color=blue,elinewidth=.6,label='Measured, corrected'); ax.errorbar(md*100,model,yerr=merr,fmt='s',markersize=4,mfc='white',mew=.9,color=orange,capsize=1.5,elinewidth=.75,label='Geant4 model');
for v in [.3,.5,.85]:ax.axvline(v*100,lw=.5,ls=':',color='#a1a7ad')
ax.set(ylabel='$C(d+21\\,\\mathrm{mm})^2$\n(m$^2$ s$^{-1}$)',ylim=(5,28));ax.legend(frameon=False,loc='upper right');ax.grid(alpha=.14)
bx.errorbar(md*100,ratio,yerr=ratioerr,fmt='o',ms=3,color=blue,capsize=1.5,elinewidth=.75); bx.axhline(1,ls='--',color=gray,lw=.9);bx.axvspan(1,10,color='#e7eff5',alpha=.55,lw=0);bx.text(1.2,.61,'Near-point test',color=blue,fontsize=7);bx.text(27,.61,'Scale points',color=gray,fontsize=7);bx.set_xscale('log');bx.set_xticks([1,2.5,10,30,50,85],labels=['1','2.5','10','30','50','85']);bx.set(xlabel='Recorded rail distance $d$ (cm)',ylabel='Measured / model',ylim=(.55,1.25),xlim=(.85,105));bx.grid(alpha=.14)
for e in ('png','pdf'):fig.savefig(P/'images'/f'model_test.{e}',dpi=300)
plt.close(fig)
