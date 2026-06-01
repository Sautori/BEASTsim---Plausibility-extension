"""Simulated against real global Moran's I, per method-variant."""
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
B=C.B; m=C.metric_long()

mr=m[m.metric=="moran_global_real"].groupby("simulator").value.mean().reindex(C.METHODS)
ms=m[m.metric=="moran_global_sim"].groupby("simulator").value.mean().reindex(C.METHODS)
fig,ax=plt.subplots(figsize=(6.6,5.2))
lim=[min(mr.min(),ms.min())*0.9, max(mr.max(),ms.max())*1.05]
ax.plot(lim,lim,color=B.NEUTRAL_LIGHT,lw=1.4,ls="--",zorder=1)
for me in C.METHODS:
    st=B.method_style(me)
    ax.scatter(mr[me],ms[me],s=160,facecolor=st["facecolor"],edgecolor=st["color"],hatch=st["hatch"],linewidth=1.7,zorder=3)
ax.set_xlabel("Moran's I, real"); B.ylabel_h(ax,"Moran's I, simulated"); ax.set_title("Spatial autocorrelation")
ax.set_xlim(lim); ax.set_ylim(lim)
B.legend_outside(ax,side="below",ncol=4,handles=B.variant_legend_handles())
C.save(fig,"fig_07_spatial_autocorr",pd.DataFrame({"method":C.METHODS,"moran_real":mr.values,"moran_sim":ms.values}))
print("spatial autocorrelation saved")
