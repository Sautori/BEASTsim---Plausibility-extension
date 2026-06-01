"""Full 23-control perturbation ledger with expected direction."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
B=C.B; pr=C.perturb(); ev=C.expected().set_index("perturbation")
sums=["Balanced Realism","Spatial Tissue Fidelity","Novelty-Gated Plausibility"]
order=["copy_reference","copy_reference_jitter_low","copy_reference_jitter_medium","copy_reference_jitter_high",
 "gene_shuffle","label_shuffle","coordinate_shuffle","missing_label_control","invalid_coordinates_control",
 "dropout_inflation_low","dropout_inflation_medium","dropout_inflation_high",
 "library_size_resampling_low","library_size_resampling_medium","library_size_resampling_high",
 "spatial_jitter_low","spatial_jitter_medium","spatial_jitter_high",
 "spatial_smoothing_low","spatial_smoothing_medium","spatial_smoothing_high",
 "random_spatial_field_injection_low","random_spatial_field_injection_medium"]
def med(p,s):
    v=pr[(pr.perturbation==p)&(pr.summary==s)].delta; return float(v.median()) if len(v) else np.nan
rows=[]
for p in order:
    d={s:med(p,s) for s in sums}
    layer=ev.expected_damaged_layer.get(p,np.nan); layer="" if (isinstance(layer,float) and np.isnan(layer)) else str(layer)
    # expected: damage should degrade Balanced Realism (delta<0). gene_shuffle is the know Novelty-Gated blind spot.
    match = "ok" if d["Balanced Realism"]<-0.005 else "--"
    rows.append((p,layer,d["Balanced Realism"],d["Spatial Tissue Fidelity"],d["Novelty-Gated Plausibility"],match))
nrow=len(rows); fig,ax=plt.subplots(figsize=(10.4,9.2)); ax.axis("off")
cols=["control","expected layer","Δ Balanced\nRealism","Δ Spatial Tissue\nFidelity","Δ Novelty-\nGated","matches\nexpected"]
xpos=[0.005,0.30,0.52,0.66,0.80,0.93]; ax.set_xlim(0,1); ax.set_ylim(0,nrow+1.5)
for x,c in zip(xpos,cols): ax.text(x,nrow+0.6,c,fontsize=8.6,fontweight="bold",color="#33414F",va="center")
ax.plot([0,1],[nrow+0.05,nrow+0.05],color="#C7CED7",lw=1)
def cellcol(v):
    if np.isnan(v): return "#999"
    t=np.clip(v/0.4,-1,1); 
    return plt.get_cmap("RdBu_r")((t+1)/2)
for i,(p,layer,a,b,c,mk) in enumerate(rows):
    y=nrow-1-i+0.5
    ax.text(xpos[0],y,p.replace("_"," "),fontsize=7.8,va="center",color="#1F2A36")
    ax.text(xpos[1],y,layer.replace("_"," "),fontsize=7.6,va="center",color="#6C7B95")
    for x,v in zip(xpos[2:5],[a,b,c]):
        ax.add_patch(plt.Rectangle((x-0.055,y-0.34),0.11,0.66,color=cellcol(v),alpha=0.85,zorder=1))
        ax.text(x,y,("0.00" if (not np.isnan(v) and abs(v)<0.005) else (f"{v:+.2f}" if not np.isnan(v) else "·")),fontsize=7.4,ha="center",va="center",
                color=("white" if (not np.isnan(v) and abs(v)>0.25) else "#1F2A36"),zorder=2)
    ax.text(xpos[5],y,mk,fontsize=8.5,ha="center",va="center",color=(B.NEUTRAL_DARK if mk=="ok" else "#B79A3A"))
    if i<nrow-1: ax.plot([0,1],[y-0.5,y-0.5],color="#EEF1F4",lw=0.6)
ax.set_title("Perturbation ledger (all 23 controls)",loc="left",x=0.0,fontsize=12,fontweight="bold",pad=10)
ax.text(0.0,-0.4,"blue = degradation (expected under damage), red = increase. gene shuffle is the documented Novelty-Gated blind spot.",fontsize=7.6,color="#6C7B95")
C.save(fig,"fig_A3_ledger",None)
print("A3 built:",nrow,"controls")
