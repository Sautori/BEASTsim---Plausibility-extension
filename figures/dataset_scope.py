"""Dataset scope table: cells, genes, cell types, runnable method-variants."""
import sys, os; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, anndata as ad, matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")
B=C.B
reals={("MERFISH","Spapros"):"real_st_spapros_merfish.h5ad",
       ("MERFISH","SpatialDE"):"real_st_spatialde_merfish.h5ad",
       ("Visium","Spapros"):"Visium_Mouse_Brain_SPAPROS_filtered_celltypes_annotated.h5ad",
       ("Visium","SpatialDE"):"Visium_Mouse_Brain_SpatialDE_filtered_celltypes_annotated.h5ad",
       ("Xenium","Xenium"):"Xenium-colon_annotated_preprocessed_0.5x_0.5y.h5ad"}
cov=C.coverage(); rows=[]
for (pl,dv),fn in reals.items():
    A=ad.read_h5ad(os.path.join(C.ROOT,"..","data",os.path.basename(fn)),backed="r")
    dl=C.dlabel(pl,dv)
    rows.append(dict(dataset=dl,platform=pl,cells=A.n_obs,genes=A.n_vars,
                     types=int(A.obs["cell_type"].nunique()),methods=int(cov[dl].notna().sum())))
    A.file.close()
sc=pd.DataFrame(rows)
cols=[("cells","cells"),("genes","genes"),("types","cell types"),("methods","method-variants")]
nrow=len(sc); ncol=len(cols)
fig,ax=plt.subplots(figsize=(7.6,3.0)); ax.axis("off"); ax.grid(False)
ax.set_xlim(0,ncol+1.25); ax.set_ylim(-0.4,nrow+1.1)
# header row (metric names)
for j,(_,lab) in enumerate(cols):
    ax.text(1.25+j+0.5,nrow+0.35,lab,ha="center",va="center",fontsize=9,color="#33414F",fontweight="bold")
ax.text(0.60,nrow+0.35,"dataset",ha="center",va="center",fontsize=9,color="#33414F",fontweight="bold")
ax.plot([0.05,ncol+1.2],[nrow+0.02,nrow+0.02],color="#C7CED7",lw=1.0)
# body rows
for i in range(nrow):
    y=nrow-1-i+0.5
    pl=sc.platform.iloc[i]
    # platform-coloured dataset chip
    ax.add_patch(plt.Rectangle((0.12,y-0.30),0.16,0.60,color=B.PLATFORM[pl],alpha=0.9))
    ax.text(0.36,y,sc.dataset.iloc[i],ha="left",va="center",fontsize=8.6,color="#1F2A36")
    for j,(col,_) in enumerate(cols):
        v=int(sc[col].iloc[i])
        ax.text(1.25+j+0.5,y,f"{v:,}",ha="center",va="center",fontsize=8.8,color="#1F2A36")
    if i<nrow-1: ax.plot([0.05,ncol+1.2],[y-0.5,y-0.5],color="#EDF0F4",lw=0.8)
ax.set_title("Dataset scope",loc="left",x=0.02,fontsize=12,fontweight="bold",pad=8)
from matplotlib.lines import Line2D
B.legend_outside(ax,side="below",ncol=3,handles=[Line2D([0],[0],marker='s',ls='',mfc=B.PLATFORM[p],mec='white',ms=10,label=p) for p in ["MERFISH","Visium","Xenium"]])
C.save(fig,"fig_03_dataset_scope",sc)
print("Fig3 -> table"); print(sc.to_string(index=False))
