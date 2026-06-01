"""Cell-type spatial overlap (CT Dice) with real and simulated tissue examples."""
import sys, os; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
B=C.B; m=C.metric_long()
man=pd.read_csv(os.path.join(C.ROOT,"..","data","spatialsimbench_input_manifest.csv"))
ctd=m[(m.platform=="Visium")&(m.dataset_variant=="Spapros")&(m.metric=="CT_mean_dice_score")].groupby("simulator").value.mean().reindex(C.METHODS).dropna()
picks=[("SRTsim-rb-domain","high"),("SRTsim-rf-domain","partial"),("SpatialcoGCN-rb","low")]
realrel="data/"+os.path.basename(man[(man.platform=="Visium")&(man.dataset_variant=="Spapros")].iloc[0].real_h5ad)
xr,lr=C.load_tissue(realrel); cmap,order=B.celltype_cmap(lr)

fig=plt.figure(figsize=(9.8,9.6)); gs=GridSpec(3,3,height_ratios=[0.8,1.05,1.05],hspace=0.5,wspace=0.06)
# --- colmun plot (CT Dice), labels ON TOP of bars ---
axb=fig.add_subplot(gs[0,:]); xb=np.arange(len(ctd)); picknames={p[0]:p[1] for p in picks}
for i,me in enumerate(ctd.index):
    st=B.method_style(me)
    axb.bar(i,ctd[me],0.66,color=st["facecolor"],edgecolor=st["color"],hatch=st["hatch"],linewidth=1.4)
    if me in picknames: axb.text(i,ctd[me]+0.02,picknames[me],ha="center",fontsize=8,color="#22303C",fontweight="bold")
axb.set_xticks(xb); axb.set_xticklabels(list(ctd.index),rotation=20,ha="right",fontsize=8.3); axb.grid(axis="x",visible=False)
axb.set_ylim(0,1.12); B.ylabel_h(axb,"cell-type overlap (CT Dice)"); axb.set_title("Cell-type spatial overlap")
# --- example tissue (real / simulated), labels between bars and images ---
fig.text(0.06,0.37,"real",rotation=90,va="center",fontsize=11,color="#2B2B2B",fontweight="bold")
fig.text(0.06,0.155,"simulated",rotation=90,va="center",fontsize=11,color="#6C7B95",fontweight="bold")
for j,(me,lab_) in enumerate(picks):
    simrel=man[(man.platform=="Visium")&(man.dataset_variant=="Spapros")&(man.simulator==me)].iloc[0].simulated_h5ad
    xs,ls=C.load_tissue(simrel)
    axr=fig.add_subplot(gs[1,j]); axs=fig.add_subplot(gs[2,j])
    for ct in order[::-1]:
        mk=lr==ct
        if mk.any(): axr.scatter(xr[mk,0],xr[mk,1],s=7,color=cmap[ct],linewidths=0,rasterized=True)
        mk2=ls==ct
        if mk2.any(): axs.scatter(xs[mk2,0],xs[mk2,1],s=9,color=cmap.get(ct,"#CCC"),linewidths=0,rasterized=True)
    _px=(xr[:,0].max()-xr[:,0].min())*0.055; _py=(xr[:,1].max()-xr[:,1].min())*0.055
    for a in (axr,axs):
        a.set_aspect("equal"); a.axis("off"); a.invert_yaxis()
        a.set_xlim(xr[:,0].min()-_px,xr[:,0].max()+_px); a.set_ylim(xr[:,1].max()+_py,xr[:,1].min()-_py)
    st=B.method_style(me)
    axr.set_title(f"{me}\nCT Dice {ctd[me]:.2f}  ·  {lab_}",fontsize=9.2,color=st["color"],fontweight="bold",pad=3)
fig.savefig(os.path.join(C.OUT,"fig_13_spatial_overlap.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(C.OUT,"fig_13_spatial_overlap.pdf"),bbox_inches="tight",dpi=150); plt.close(fig)
print("Fig13 rebuilt")
