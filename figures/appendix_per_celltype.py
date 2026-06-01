"""Per-cell-type spatial overlap on Visium Spapros, recomputed from the tissues."""
import sys, os; sys.path.insert(0,".")
import common as C, numpy as np, matplotlib.pyplot as plt
B=C.B
files={"SRTsim-rb-domain":"SRTsim_visium_spapros_rb_domain.h5ad",
       "SRTsim-rb-tissue":"SRTsim_visium_spapros_rb_tissue.h5ad",
       "SRTsim-rf-domain":"SRTsim_visium_spapros_rf_domain_grid_c1.h5ad",
       "SRTsim-rf-tissue":"SRTsim_visium_spapros_rf_tissue_grid_c1.h5ad",
       "scCube-rfb":"scCube_visium_spapros_rfb.h5ad",
       "scDesign3-rb":"scDesign3_visium_spapros_rb_p1_f1_b0.h5ad",
       "SpatialcoGCN-rb":"spatialcoGCN_-_rb_-_spapros.h5ad"}
real_file="Visium_Mouse_Brain_SPAPROS_filtered_celltypes_annotated.h5ad"
cats=[str(i) for i in range(6)]; G=24
def occ(rel):
    xy,lab=C.load_tissue("data/"+rel)
    mn=xy.min(0); rng=xy.max(0)-mn; rng[rng==0]=1
    u=(xy-mn)/rng                              # normalise to unit square (scale-invariant)
    gx=np.clip((u[:,0]*G).astype(int),0,G-1); gy=np.clip((u[:,1]*G).astype(int),0,G-1)
    out={}
    for c in cats:
        m=lab==c; out[c]=set(zip(gx[m],gy[m]))
    return out
R=occ(real_file)
def dice(a,b):
    if not a and not b: return np.nan
    inter=len(a&b); return 2*inter/(len(a)+len(b)) if (len(a)+len(b)) else np.nan
M=np.full((len(files),6),np.nan)
for i,(me,fn) in enumerate(files.items()):
    S=occ(fn)
    for j,c in enumerate(cats): M[i,j]=dice(R[c],S[c])
fig,ax=plt.subplots(figsize=(7.0,4.6)); ax.grid(False)
im=ax.imshow(M,cmap=B.SCORE_CMAP,aspect="auto",vmin=0,vmax=1)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v=M[i,j]
        if np.isnan(v): B.mark_not_run(ax,j,i)
        else: ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=8,color=B.heat_text_color(v,0,1,B.SCORE_CMAP))
ax.set_xticks(range(6)); ax.set_xticklabels(["ct "+c for c in cats],fontsize=8.5)
ax.set_yticks(range(len(files))); ax.set_yticklabels(list(files.keys()),fontsize=9)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.outline.set_visible(False)
cb.ax.set_title("occupancy\nDice",fontsize=8,color="#6C7B95",loc="left")
ax.set_title("Per-cell-type spatial overlap (Visium Spapros)")
C.save(fig,"fig_A2_per_celltype",None)
print("A2 built; row means:", np.round(np.nanmean(M,1),2))
