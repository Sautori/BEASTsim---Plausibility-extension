"""Are the fidelity axes independent evidence, or one quantity seen six times? Spearman
correlation across the eight variants among the six fidelity axes the radar shows, plus
copy-risk. Most of them co-vary and track copy-risk, so the radar is closer to one latent
fidelity dimension than to six independent wins. Cell-type proportion is the partial exception."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
m=C.load("metric_long.csv")
axes={"DE signal":("DE_lfc_correlation",1),"SVG recovery":("SVG_precision",1),
      "CT Dice":("CT_mean_dice_score",1),"clustering":("cluster_ARI_real_vs_sim",1),
      "CT proportion":("ct_proportion_JSD",-1),"spatial autocorr":("moran_profile_cosine",1)}
piv=m.pivot_table(index="simulator",columns="metric",values="value",aggfunc="mean")
F=pd.DataFrame({k:s*piv[v] for k,(v,s) in axes.items()}); F["copy-risk"]=piv["copy_similarity"]
R=F.corr(method="spearman"); labels=list(R.columns); n=len(labels)
fig,ax=plt.subplots(figsize=(6.6,5.6))
im=ax.imshow(R.values,cmap="plasma",vmin=0,vmax=1)
ax.set_xticks(range(n)); ax.set_yticks(range(n))
ax.set_xticklabels(labels,rotation=40,ha="right",fontsize=8.5); ax.set_yticklabels(labels,fontsize=8.5)
for i in range(n):
    for j in range(n):
        v=R.values[i,j]
        ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=7.6,color=("white" if v<0.6 else "#1b1b1b"))
ax.axhline(n-1.5,color="white",lw=2); ax.axvline(n-1.5,color="white",lw=2)
fa=[c for c in F.columns if c!="copy-risk"]
sub=R.loc[fa,fa].values; off=sub[np.triu_indices(len(fa),1)]
ax.set_title("Fidelity axes co-vary and track copy-risk")
ax.text(0.0,-0.34,f"mean pairwise corr among the six axes = {np.nanmean(off):.2f};  "
                  f"mean corr with copy-risk = {R.loc[fa,'copy-risk'].mean():.2f}",
        transform=ax.transAxes,fontsize=8.2,color="#6C7B95")
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.outline.set_visible(False)
cb.ax.set_title("Spearman",fontsize=7.5,color="#6C7B95",loc="left")
C.save(fig,"fig_33_fidelity_axis_correlation",None)
print("fig_33 built")
