"""The three real tissues with their cell-type composition bars."""
import sys, os; sys.path.insert(0,".")
import common as C, numpy as np, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
B=C.B
sets=[("MERFISH","intestine","data/real_st_spapros_merfish.h5ad"),
      ("Visium","mouse brain","data/Visium_Mouse_Brain_SPAPROS_filtered_celltypes_annotated.h5ad"),
      ("Xenium","colon","data/Xenium-colon_annotated_preprocessed_0.5x_0.5y.h5ad")]
fig=plt.figure(figsize=(9.6,8.8)); gs=GridSpec(3,2,width_ratios=[3.7,1.7],hspace=0.40,wspace=0.10)
for r,(pl,organ,rel) in enumerate(sets):
    xy,lab=C.load_tissue(rel); cmap,order=B.celltype_cmap(lab)
    axt=fig.add_subplot(gs[r,0])
    for ct in order[::-1]:
        msk=lab==ct
        if msk.any(): axt.scatter(xy[msk,0],xy[msk,1],s=(3.0 if pl=="Xenium" else 9),color=cmap[ct],linewidths=0,rasterized=True)
    axt.set_aspect("equal"); axt.axis("off"); axt.invert_yaxis()
    axt.set_title(f"{pl}  ·  {organ}",loc="left",color=B.PLATFORM[pl],fontsize=13,fontweight="bold")
    # compositon bar only (no per-segment labels), caption centered over it
    axb=fig.add_subplot(gs[r,1])
    counts=np.array([(lab==ct).sum() for ct in order],float)
    si=np.argsort(-counts); order_b=[order[k] for k in si]; counts=counts[si]
    frac=counts/counts.sum(); left=0
    for ct,fr in zip(order_b,frac):
        axb.barh(0,fr,left=left,height=0.5,color=cmap[ct],edgecolor="white",linewidth=0.3); left+=fr
    axb.set_xlim(0,1); axb.set_ylim(-0.6,0.6); axb.axis("off")
    axb.text(0.5,-0.42,f"{len(order)} cell types · {len(lab):,} cells",ha="center",fontsize=8.4,color="#6C7B95")
fig.suptitle("The three tissues, and what they're made of",x=0.5,y=0.97,fontweight="bold",fontsize=14)
fig.savefig(os.path.join(C.OUT,"fig_02_dataset_gallery.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(C.OUT,"fig_02_dataset_gallery.pdf"),bbox_inches="tight",dpi=160); plt.close(fig)
print("Fig2 updated")
