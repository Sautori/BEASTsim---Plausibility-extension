"""Correspondence-dependent ARI against label-based clustering quality."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, matplotlib.pyplot as plt
B=C.B
m=C.metric_long()
def mm(metric): return m[m.metric==metric].groupby("simulator").value.mean().reindex(C.METHODS)
x=mm("cluster_ARI_real_vs_sim").clip(lower=0)
ys=[("ARI vs cell-type label","cluster_ARI_real_vs_label"),("NMI vs cell-type label","cluster_NMI_real_vs_label")]
fig,axs=plt.subplots(1,2,figsize=(11.2,5.0),sharex=True,sharey=False)
for ax,(ylab,met) in zip(axs,ys):
    y=mm(met)
    ax.fill_between([0,1],[0,1],[1,1],color="#ECEEF1",zorder=0)          # upper-left: correspondence-limited
    ax.plot([0,1],[0,1],ls="--",color="#9AA6B2",lw=1.2,zorder=1)
    for me in C.METHODS:
        st=B.method_style(me)
        ax.scatter(x[me],y[me],s=165,facecolor=st["facecolor"],edgecolor=st["color"],
                   hatch=st["hatch"],linewidth=1.7,zorder=3)
    ax.set_xlim(-0.03,1.03); ax.set_ylim(0,0.72)
    ax.set_xlabel("clustering ARI vs simulated (needs cell correspondence)")
    B.ylabel_h(ax,ylab)
    ax.set_title(ylab.split(" vs ")[0],fontsize=11)
axs[0].text(0.03,0.66,"correspondence-limited\n(no cell match)",fontsize=8.5,color="#6C7B95",va="top")
B.legend_outside(axs[1],side="below",ncol=4,handles=B.variant_legend_handles())
fig.suptitle("Fair clustering: correspondence-free quality vs correspondence-dependent ARI",
             x=0.5,y=1.02,fontweight="bold",fontsize=13)
C.save(fig,"fig_30_clustering_fair",None)
print("fig_30 built")
