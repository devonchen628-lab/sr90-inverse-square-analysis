#!/usr/bin/env python3
"""Single-panel far-field comparison; preserves every simulated position."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from report_style import apply, BLUE, RUST

P=Path(__file__).resolve().parent
F=P.parent/'geant4_farfield_validation'

def main():
    apply()
    plt.rcParams.update({'axes.labelsize':9.5, 'xtick.labelsize':8,
                         'ytick.labelsize':8, 'legend.fontsize':7.5})
    frame=pd.read_csv(F/'transmission_comparison.csv').sort_values('d_rail_m')
    fit=pd.read_csv(F/'experimental_farfield_fit.csv').iloc[0]
    mcfit=pd.read_csv(F/'geant4_farfield_fit.csv').iloc[0]
    cutoff=pd.read_csv(F/'farfield_cutoff_analysis.csv')
    cut=float(cutoff.loc[np.isclose(cutoff.tolerance,.01),'first_included_rail_m'].iloc[0])
    d=frame.d_rail_m.to_numpy(float)
    y=frame.T_mc_normalized.to_numpy(float)
    lo=frame.T_mc_p16.to_numpy(float);hi=frame.T_mc_p84.to_numpy(float)
    assert len(d)==32 and np.sum(frame.selected)==30 and np.isclose(cut,.05)
    assert np.array_equal(frame.selected.to_numpy(),d>=cut-1e-10)
    idx=np.flatnonzero(np.isclose(d,.30))[0]
    assert np.isclose(y[idx],1) and np.isclose(frame.T_experimental_m3.iloc[idx],1)
    assert np.isclose(mcfit.reference_rail_m,.30) and np.all(lo<=y) and np.all(hi>=y)
    assert np.allclose(frame.T_experimental_m3,np.exp(-fit.mu_per_m*(d-.30)))
    xline=np.linspace(cut,d.max(),500)
    model=np.exp(-fit.mu_per_m*(xline-.30))
    fig,ax=plt.subplots(figsize=(3.5,2.6))
    ax.axvspan(0,cut*100,color='#e8e8e8',alpha=.7,lw=0,zorder=0)
    ax.axvline(cut*100,color='#858585',ls=(0,(3,3)),lw=.65,zorder=1)
    ax.text(7.0,1.325,'Far-field cutoff: 5 cm',fontsize=7.5,color='#555555',va='center')
    eb=ax.errorbar(d*100,y,yerr=np.vstack([y-lo,hi-y]),fmt='o',ms=2.5,
        color=BLUE,ecolor=BLUE,elinewidth=.6,capsize=1.25,capthick=.6,
        label='Geant4 air/vacuum',zorder=3)
    line,=ax.plot(xline*100,model,color=RUST,lw=1.35,
        label='Experimental M3 (far-field fit)',zorder=2)
    mc_model=np.exp(-mcfit.mu_best_per_m*(xline-.30))
    mc_line,=ax.plot(xline*100,mc_model,color=BLUE,ls='--',lw=1.15,
        label='Geant4 exponential fit',zorder=2)
    ax.set_xlabel(r'Rail separation, $d$ (cm)')
    ax.set_ylabel('Relative detector response')
    ax.set_xlim(0,88);ax.set_ylim(.37,1.38)
    ax.legend([eb,mc_line,line],['Geant4 air/vacuum','Geant4 exponential fit','Experimental M3 (far-field fit)'],
              loc='lower left',frameon=False,handlelength=1.8,labelspacing=.25,borderaxespad=.3)
    fig.subplots_adjust(left=.21,right=.98,bottom=.19,top=.97)
    for ext in ('pdf','png'):fig.savefig(P/f'figure2_geant4_simplified.{ext}',dpi=400)
    plt.close(fig)
    caption=('Geant4 air-to-vacuum detector-response ratio and the experimental M3 fitted attenuation term, '
        'each normalised to unity at 30 cm rail separation. All 32 simulation positions are shown; '
        'vertical error bars show pointwise 16–84% batch-bootstrap intervals. The light-grey region '
        'marks distances excluded from the attenuation fits by the predefined 1% analytical geometry '
        'criterion. At large separations the simulated response declines more strongly than the '
        'experimental fitted term. The dashed blue line is the anchored far-field exponential summary '
        'of the simulated ratio, which is not required to follow an exact exponential law. '
        'The red curve represents a model fit rather than directly measured '
        'experimental transmission.\n')
    (P/'figure2_geant4_simplified_caption.txt').write_text(caption)

if __name__=='__main__':main()
