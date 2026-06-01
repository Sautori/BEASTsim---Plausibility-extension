"""Datasets run per method-variant, by platform."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B; PL=["MERFISH","Visium","Xenium"]
cov=C.coverage(); dl2pl={C.dlabel(p,v):p for p,v in C.DATASETS}
counts={me:{pl:0 for pl in PL} for me in C.METHODS}
for me in C.METHODS:
    for dl in C.DLABELS:
        if not np.isnan(cov.loc[me,dl]): counts[me][dl2pl[dl]]+=1
order=sorted(C.METHODS,key=lambda me:sum(counts[me].values()))
fig,ax=plt.subplots(figsize=(7.6,4.6)); ax.grid(axis="y",visible=False)
y=np.arange(len(order))
for i,me in enumerate(order):
    left=0
    for pl in PL:
        c=counts[me][pl]
        if c>0:
            ax.barh(i,c,left=left,height=0.62,color=B.PLATFORM[pl],edgecolor="white",linewidth=1.1)
            left+=c
    ax.text(left+0.08,i,str(int(left)),va="center",ha="left",fontsize=9,color="#33414F",fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels(order,fontsize=9); ax.set_xlim(0,5.5); ax.set_xticks(range(6))
ax.set_xlabel("datasets run (of 5)")
for s in ["top","right","left"]: ax.spines[s].set_visible(False)
ax.tick_params(length=0); ax.set_title("Evidence coverage")
B.legend_outside(ax,side="below",ncol=3,handles=[Line2D([0],[0],marker='s',ls='',mfc=B.PLATFORM[p],mec='white',ms=11,label=p) for p in PL])
C.save(fig,"fig_28_coverage_summary",None)
print("fig_28 coverage summary built")
