"""Novelty-Gated Plausibility across a range of gate thresholds."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
B=C.B; nm=C.named()
bal=nm[nm.summary=="Balanced Realism"].groupby("simulator").score.mean().reindex(C.METHODS)
copy=nm[nm.summary=="copy_similarity"].groupby("simulator").score.mean().reindex(C.METHODS)
W=0.05; taus=[0.85,0.90,0.95,0.98,1.00]
def ng(b,c,t): return b*np.clip(1-max(c-t,0)/W,0,1)
M=np.array([[ng(bal[me],copy[me],t) for t in taus] for me in C.METHODS])
fig,ax=plt.subplots(figsize=(6.8,4.8)); ax.grid(False)
im=ax.imshow(M,cmap=B.SCORE_CMAP,aspect="auto",vmin=0,vmax=1)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j,i,f"{M[i,j]:.2f}",ha="center",va="center",fontsize=8,color=B.heat_text_color(M[i,j],0,1,B.SCORE_CMAP))
# highlight defualt 0.95 column
dj=taus.index(0.95)
ax.add_patch(Rectangle((dj-0.5,-0.5),1,M.shape[0],fill=False,edgecolor=B.ALERT,lw=2.4,zorder=4))
ax.text(dj,-0.85,"default",ha="center",color=B.ALERT,fontsize=8.5,fontweight="bold")
ax.set_xticks(range(len(taus))); ax.set_xticklabels([f"{t:.2f}" for t in taus])
ax.set_yticks(range(len(C.METHODS))); ax.set_yticklabels(C.METHODS,fontsize=9)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0); ax.set_xlabel("novelty-gate threshold τ")
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.outline.set_visible(False); cb.ax.set_title("Novelty-Gated\nPlausibility",fontsize=7.5,color="#6C7B95",loc="left")
ax.set_title("Gate calibration")
C.save(fig,"fig_26_gate_calibration",pd.DataFrame(M,index=C.METHODS,columns=[f"tau_{t}" for t in taus]).reset_index())
print("Fig26 heatmap done")
