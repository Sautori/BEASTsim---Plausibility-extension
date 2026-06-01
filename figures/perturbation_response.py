"""Median change in the candidate summaries under each perturbation control."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
B=C.B; pr=C.perturb()
sums=["Balanced Realism","Spatial Tissue Fidelity","Novelty-Gated Plausibility"]
g=pr[pr.summary.isin(sums)].groupby(["perturbation","summary"]).delta.median().unstack()[sums]
g=g.reindex(g["Balanced Realism"].sort_values().index)              # most-degrading at top
M=g.values; vmax=0.4
fig,ax=plt.subplots(figsize=(7.4,8.2)); ax.grid(False)
im=ax.imshow(M,cmap="RdBu_r",aspect="auto",vmin=-vmax,vmax=vmax)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v=M[i,j]
        if np.isnan(v): continue
        txt=("0.00" if abs(v)<0.005 else (f"{v:+.2f}" if abs(v)<1 else f"{v:+.1f}"))
        ax.text(j,i,txt,ha="center",va="center",fontsize=7,
                color=("white" if abs(v)>vmax*0.6 else "#222"))
ax.set_xticks(range(3)); ax.set_xticklabels(["Balanced\nRealism","Spatial Tissue\nFidelity","Novelty-Gated"],fontsize=9)
ax.set_yticks(range(M.shape[0])); ax.set_yticklabels([p.replace("_"," ") for p in g.index],fontsize=7.6)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=0.05,pad=0.03); cb.outline.set_visible(False)
cb.ax.set_title("median Δ\nfrom baseline",fontsize=7.5,color="#6C7B95",loc="left")
ax.set_title("Perturbation response",pad=30)
ax.text(0.0,1.010,"blue = degraded (expected) · red = increased",transform=ax.transAxes,fontsize=8,color="#6C7B95")
C.save(fig,"fig_21_expected_vs_observed",None)
print("controls:",M.shape[0],"| BR responds (neg):",int((g['Balanced Realism']<-0.02).sum()),"/",M.shape[0])
