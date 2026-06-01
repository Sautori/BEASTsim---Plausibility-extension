"""Cell-type neighbour co-occurrence matrices: real, a faithful sim, a distorting sim."""
import sys, os; sys.path.insert(0,".")
import common as C, numpy as np, matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
B=C.B
panels=[("real","data/Visium_Mouse_Brain_SPAPROS_filtered_celltypes_annotated.h5ad",None),
        ("SRTsim-rb  (proximity 1.00)","data/SRTsim_visium_spapros_rb_domain.h5ad","faithful"),
        ("SpatialcoGCN  (proximity 0.36)","data/spatialcoGCN_-_rb_-_spapros.h5ad","distorts")]
K=6; cats=[str(i) for i in range(6)]
def cooc(xy,lab):
    nn=NearestNeighbors(n_neighbors=K+1).fit(xy); _,idx=nn.kneighbors(xy)
    idx=idx[:,1:]                                   # drop self
    M=np.zeros((6,6))
    code={c:i for i,c in enumerate(cats)}
    lc=np.array([code.get(l,-1) for l in lab])
    for i in range(len(lab)):
        a=lc[i]
        if a<0: continue
        for j in idx[i]:
            b=lc[j]
            if b>=0: M[a,b]+=1
    rs=M.sum(1,keepdims=True); rs[rs==0]=1
    return M/rs
mats=[(t,cooc(*C.load_tissue(rel)),tag) for t,rel,tag in panels]
vmax=max(m.max() for _,m,_ in mats)
fig,axs=plt.subplots(1,3,figsize=(12.6,4.5))
for ax,(lab,Mx,tag) in zip(axs,mats):
    ax.grid(False); im=ax.imshow(Mx,cmap=B.SCORE_CMAP,vmin=0,vmax=vmax,aspect="auto")
    for i in range(6):
        for j in range(6):
            ax.text(j,i,f"{Mx[i,j]:.2f}",ha="center",va="center",fontsize=7.4,
                    color=B.heat_text_color(Mx[i,j],0,vmax,B.SCORE_CMAP))
    ax.set_xticks(range(6)); ax.set_xticklabels(["ct "+c for c in cats],fontsize=7.6)
    ax.set_yticks(range(6)); ax.set_yticklabels(["ct "+c for c in cats],fontsize=7.6)
    ax.set_xlabel("neighbour cell type",fontsize=9)
    if lab=="real": ax.set_ylabel("focal cell type",fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    col="#2B2B2B" if lab=="real" else (B.SIM["SRTsim"] if tag=="faithful" else B.ALERT)
    ax.set_title(lab,fontsize=10.5,color=col,pad=8)
cb=fig.colorbar(im,ax=axs,fraction=0.026,pad=0.03); cb.outline.set_visible(False)
cb.ax.set_title("P(neighbour\n| focal)",fontsize=7.5,color="#6C7B95",loc="left")
fig.suptitle("Cell-type neighbourhood co-occurrence  ·  Visium Spapros",x=0.5,y=1.04,fontweight="bold",fontsize=13)
fig.savefig(os.path.join(C.OUT,"fig_12_neighbourhood_matrices.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(C.OUT,"fig_12_neighbourhood_matrices.pdf"),bbox_inches="tight",dpi=180); plt.close(fig)
std,ok=C.greyscale_ok("fig_12_neighbourhood_matrices.png")
print(f"Fig12 rebuilt: real vs faithful vs distorting  [greyscale sigma={std:.1f} {'ok' if ok else 'LOW'}]")
