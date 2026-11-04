"""Consistent compact style for one-column journal figures."""
import matplotlib as mpl

def apply():
    mpl.rcParams.update({
        'font.family':'DejaVu Sans', 'font.size':8, 'axes.labelsize':8,
        'xtick.labelsize':7, 'ytick.labelsize':7, 'legend.fontsize':6.7,
        'axes.linewidth':.7, 'lines.linewidth':1.2,
        'xtick.major.size':3, 'ytick.major.size':3,
        'xtick.major.width':.7, 'ytick.major.width':.7,
        'xtick.direction':'out', 'ytick.direction':'out',
        'axes.spines.top':False, 'axes.spines.right':False,
        'pdf.fonttype':42, 'ps.fonttype':42, 'savefig.dpi':400,
    })

BLUE='#326b8d'
RUST='#b4422a'
