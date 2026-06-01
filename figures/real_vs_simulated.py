"""Real Visium section beside four simulated outputs."""
import sys, os; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B
man=pd.read_csv(os.path.join(C.ROOT,"..","data","spatialsimbench_input_manifest.csv"))
sel=man[(man.platform=="Visium")&(man.dataset_variant=="Spapros")&(man.simulator.isin(
        ["SRTsim-rb-domain","SRTsim-rf-domain","scCube-rfb","scDesign3-rb"]))].copy()
xr,lr=C.load_tissue("data/"+os.path.basename(sel.iloc[0].real_h5ad))
cmap,order=B.celltype_cmap(lr)
panels=[("real",xr,lr)]+[(r.simulator,*C.load_tissue(r.simulated_h5ad)) for _,r in sel.iterrows()]
fig,axs=plt.subplots(1,len(panels),figsize=(3.0*len(panels),3.3))
for ax,(name,xy,lab) in zip(axs,panels):
    for ct in order[::-1]:
        msk=lab==ct
        if msk.any(): ax.scatter(xy[msk,0],xy[msk,1],s=6,color=cmap.get(ct,"#CCCCCC"),linewidths=0,rasterized=True)
    ax.set_aspect("equal"); ax.axis("off"); ax.invert_yaxis()
    tc=B.SIM[B.family_of(name)] if name!="real" else "#2B2B2B"
    ax.set_title(name,fontsize=10,color=tc,fontweight="bold")
fig.suptitle("Real vs simulated tissue",x=0.5,y=1.06,fontweight="bold",fontsize=13)
leg_order=sorted(order,key=C.natkey)
hand=[Line2D([0],[0],marker='o',ls='',mfc=cmap[ct],mec='none',ms=9,label=f"cluster {ct}") for ct in leg_order]
fig.legend(handles=hand,loc="lower center",bbox_to_anchor=(0.5,-0.10),ncol=len(leg_order),frameon=False,fontsize=8.5,title="Visium cell-type clusters")
fig.savefig(os.path.join(C.OUT,"fig_04_real_vs_sim_tissue.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(C.OUT,"fig_04_real_vs_sim_tissue.pdf"),bbox_inches="tight",dpi=180); plt.close(fig)
print("Fig4 re-rendered | legend:",leg_order)
