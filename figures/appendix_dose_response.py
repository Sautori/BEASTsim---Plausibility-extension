"""Dose-response of the candidate summaries to graded perturbation intensity."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B; pr=C.perturb()
BR,STF="#3A4654","#2A9D8F"   # neutral, non-family
fams=[("copy-reference jitter",["copy_reference_jitter_low","copy_reference_jitter_medium","copy_reference_jitter_high"]),
      ("dropout inflation",["dropout_inflation_low","dropout_inflation_medium","dropout_inflation_high"]),
      ("library-size resampling",["library_size_resampling_low","library_size_resampling_medium","library_size_resampling_high"]),
      ("spatial jitter",["spatial_jitter_low","spatial_jitter_medium","spatial_jitter_high"]),
      ("spatial smoothing",["spatial_smoothing_low","spatial_smoothing_medium","spatial_smoothing_high"]),
      ("random spatial field",["random_spatial_field_injection_low","random_spatial_field_injection_medium"])]
def med(p,s): 
    v=pr[(pr.perturbation==p)&(pr.summary==s)].delta
    return float(v.median()) if len(v) else np.nan
fig,axs=plt.subplots(2,3,figsize=(11.0,6.4)); axs=axs.ravel()
for ax,(name,ctrls) in zip(axs,fams):
    xs=np.arange(len(ctrls)); xlab=[c.rsplit("_",1)[1] for c in ctrls]
    for s,col in [("Balanced Realism",BR),("Spatial Tissue Fidelity",STF)]:
        y=[med(c,s) for c in ctrls]
        ax.plot(xs,y,marker="o",ms=7,color=col,lw=2.0,label=s)
        nonmono=[i for i in range(1,len(y)) if y[i]>y[i-1]+1e-6]   # damage should not lessen with intensity
        for i in nonmono: ax.scatter(xs[i],y[i],s=120,facecolor="none",edgecolor=B.ALERT,linewidth=1.8,zorder=4)
    ax.axhline(0,color="#C8D0DA",lw=1); ax.set_xticks(xs); ax.set_xticklabels(xlab,fontsize=8.5)
    ax.set_title(name,fontsize=10); ax.grid(axis="x",visible=False)
for ax in axs: B.ylabel_h(ax,"median Δ") if ax in (axs[0],axs[3]) else None
B.legend_outside(axs[-2],side="below",ncol=2,handles=[Line2D([0],[0],marker='o',color=BR,label='Balanced Realism'),
    Line2D([0],[0],marker='o',color=STF,label='Spatial Tissue Fidelity'),
    Line2D([0],[0],marker='o',ls='',mfc='none',mec=B.ALERT,mew=1.8,ms=11,label='non-monotone')])
fig.suptitle("Dose-response to graded perturbation intensity",x=0.5,y=1.01,fontweight="bold",fontsize=13)
fig.tight_layout()
C.save(fig,"fig_A1_dose_response",None)
print("A1 built")
