#!/usr/bin/env python3
"""Draw prior far-field results; does not rerun fits or simulations."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from report_style import apply,BLUE,RUST

P=Path(__file__).resolve().parent
F=P.parent/'geant4_farfield_validation'

def main():
    apply()
    all_points=pd.read_csv(F/'transmission_comparison.csv')
    fit=pd.read_csv(F/'experimental_farfield_fit.csv').iloc[0]
    mcfit=pd.read_csv(F/'geant4_farfield_fit.csv').iloc[0]
    cutoff=pd.read_csv(F/'farfield_cutoff_analysis.csv')
    tolerance=cutoff[np.isclose(cutoff.tolerance,.01)].iloc[0]
    d=all_points[all_points.selected].copy().sort_values('d_rail_m').reset_index(drop=True)
    x=d.d_rail_m.to_numpy(float);i=int(np.flatnonzero(np.isclose(x,.30))[0])
    assert len(d)==int(fit.n)==int(mcfit.n_selected)==30
    assert np.isclose(x[0],tolerance.first_included_rail_m) and np.isclose(x[0],.05)
    assert np.isclose(mcfit.reference_rail_m,.30)
    assert np.isclose(d.T_mc_normalized.iloc[i],1) and np.isclose(d.T_experimental_m3.iloc[i],1)
    assert np.allclose(d.T_experimental_m3,np.exp(-fit.mu_per_m*(x-.30)))
    assert np.allclose(d.difference_mc_minus_exp,d.T_mc_normalized-d.T_experimental_m3)
    fig,(ax,bx)=plt.subplots(2,1,figsize=(3.5,4.0),sharex=True,
        gridspec_kw={'height_ratios':[2.15,1],'hspace':.06})
    ax.fill_between(x*100,d.T_mc_p16,d.T_mc_p84,color=BLUE,alpha=.19,lw=0,label='Geant4 16–84%')
    ax.plot(x*100,d.T_mc_normalized,'o',ms=2.5,color=BLUE,ls='none',label='Geant4 air / vacuum',zorder=3)
    ax.fill_between(x*100,d.T_experimental_p16,d.T_experimental_p84,color=RUST,alpha=.15,lw=0,label='M3 formal interval')
    ax.plot(x*100,d.T_experimental_m3,color=RUST,lw=1.5,label='Experimental M3')
    ax.axvline(30,color='#999999',ls=':',lw=.7);ax.axvline(5,color='#777777',ls='--',lw=.7)
    ax.text(.045,.94,'(a)',transform=ax.transAxes,va='top',weight='bold')
    ax.text(.97,.94,'5 cm cutoff',transform=ax.transAxes,ha='right',va='top',fontsize=7)
    ax.set_xlim(3,88);ax.set_ylim(.36,1.25);ax.set_ylabel('Response / response\nat 30 cm')
    ax.legend(frameon=False,loc='lower left',handlelength=1.7,labelspacing=.22,borderaxespad=.2)
    bx.fill_between(x*100,d.difference_boot_p16,d.difference_boot_p84,color='#708fa4',alpha=.23,lw=0,label='Pointwise 16–84%')
    bx.plot(x*100,d.difference_mc_minus_exp,'o',ms=2.5,color=BLUE,ls='none')
    bx.axhline(0,color='#333333',lw=.7);bx.axvline(30,color='#999999',ls=':',lw=.7);bx.axvline(5,color='#777777',ls='--',lw=.7)
    bx.text(.045,.92,'(b)',transform=bx.transAxes,va='top',weight='bold')
    bx.set_ylabel('Geant4 − M3');bx.set_xlabel('Rail separation $d$ (cm)')
    bx.set_ylim(-.36,.35);bx.set_yticks([-.3,0,.3])
    bx.legend(frameon=False,loc='upper right',handlelength=1.4,borderaxespad=.2)
    fig.subplots_adjust(left=.19,right=.98,bottom=.11,top=.98)
    for ext in ('pdf','png'):fig.savefig(P/f'figure2_geant4_main.{ext}',dpi=400)
    plt.close(fig)
    caption=('Far-field comparison of the Geant4 air-to-vacuum response ratio and the effective attenuation '
      'inferred from experimental model M3. (a) Both responses are normalised at the 30 cm rail position; '
      'the far-field selection begins at 5 cm, from a 1% analytical finite-geometry criterion. '
      'Geant4 bands show pointwise 16–84% whole-batch bootstrap intervals; the narrower M3 band is a formal '
      'profile interval conditional on the experimental model. (b) Geant4 minus M3, with pointwise '
      'conditional 16–84% intervals combining the Geant4 bootstrap and formal experimental uncertainty. '
      'Both responses decrease overall, while Geant4 falls more strongly at the longest distances; this '
      'supports qualitative agreement but does not establish quantitative agreement. The intervals do not '
      'include uncertainties in source geometry, materials, threshold or response models; the vacuum and '
      'analytical-geometry responses also have different distance dependences.\n')
    (P/'figure2_caption.txt').write_text(caption)

if __name__=='__main__':main()
