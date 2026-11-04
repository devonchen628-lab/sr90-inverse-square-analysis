"""Draw the conceptual rail-reference and active-plane distance diagram."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle

P=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(3.5,1.42))
ax.set_xlim(0,10);ax.set_ylim(0,4.3);ax.axis('off')
ax.add_patch(Rectangle((.55,.18),8.9,.21,fc='#333333'))
for x in (2.75,6.8):
    ax.add_patch(Rectangle((x-.38,.4),.76,.34,fc='#333333'))
    ax.add_patch(Rectangle((x-.1,.74),.2,1.35,fc='#888888'))
ax.add_patch(Rectangle((2.35,2.05),.38,1.1,fc='#333333'))
ax.add_patch(Rectangle((6.82,2.05),.48,1.1,fc='#386b85'))
ax.plot([2.73,2.73],[2.32,2.88],color='#bd4b32',lw=2)
ax.plot([6.82,6.82],[2.37,2.83],color='#bd4b32',lw=2)
ax.text(2.54,3.19,'Detector',ha='center',va='bottom')
ax.text(7.06,3.19,'Source',ha='center',va='bottom')
ax.annotate('',xy=(2.73,3.75),xytext=(6.82,3.75),arrowprops=dict(arrowstyle='<->',lw=.75,color='#bd4b32'))
ax.text(4.78,3.82,r'$r=d+\delta$',ha='center',va='bottom',color='#9b3a2b')
for x in (2.73,6.82):ax.plot([x,x],[3.22,3.72],color='#bd4b32',lw=.65,ls=':')
for x in (3.13,6.42):ax.plot([x,x],[.5,1.7],color='#386b85',lw=.65,ls=':')
ax.annotate('',xy=(3.13,1.18),xytext=(6.42,1.18),arrowprops=dict(arrowstyle='<->',lw=.75,color='#386b85'))
ax.text(4.77,1.26,r'$d$ (rail reference)',ha='center',va='bottom',color='#386b85',fontsize=7.5)
ax.text(5,.42,'rail',ha='center',va='bottom',color='#444444',fontsize=7)
fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.97)
for ext in ('pdf','png'):fig.savefig(P/'images'/('setup_schematic.'+ext),dpi=300)
