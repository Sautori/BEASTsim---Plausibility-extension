"""Candidate-summary panels: scores, criterion ranks, run coverage."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B; m=C.metric_long(); nm=C.named(); ce=C.components(); pr=C.perturb(); ev=C.expected()
RB,RF=B.REFMODE["rb"],B.REFMODE["rf"]
def mm(fam=None,metric=None):
    d=m.copy()
    if fam:d=d[d.family==fam]
    if metric:d=d[d.metric==metric]
    return d.groupby(["platform","dataset_variant","simulator"]).value.mean().reset_index()
def bym(df): return df.groupby("simulator").value.mean().reindex(C.METHODS)
SUMS=["Balanced Realism","Spatial Tissue Fidelity","Novelty-Gated Plausibility","Conservative Softmin Plausibility"]
def sum_by_method(s): return nm[nm.summary==s].groupby("simulator").score.mean().reindex(C.METHODS)

# ===== Fig19: BP criterion disagreement (bump chart) =====
R={s:sum_by_method(s).rank(ascending=False,method="min") for s in SUMS}
fig,ax=plt.subplots(figsize=(7.6,5.2)); xs=np.arange(len(SUMS))
for me in C.METHODS:
    st=B.method_style(me); ys=[R[s][me] for s in SUMS]
    ax.plot(xs,ys,color=st["color"],ls=st["ls"],lw=2.2,marker="o",ms=9,mfc=st["facecolor"],mec=st["color"],alpha=0.9,zorder=3)
    ax.text(xs[-1]+0.06,ys[-1],me,va="center",fontsize=7.6,color=st["color"])
ax.set_yticks(range(1,9)); ax.invert_yaxis(); ax.set_xticks(xs)
ax.set_xticklabels(["Balanced\nRealism","Spatial Tissue\nFidelity","Novelty-Gated","Conservative\nSoftmin"],fontsize=8.5)
ax.grid(axis="x",visible=False); B.ylabel_h(ax,"rank (1 = best)"); ax.set_title("Criterion disagreement"); ax.set_xlim(-0.2,len(SUMS)+0.9)
C.save(fig,"fig_19_criterion_disagreement",pd.DataFrame({s:R[s] for s in SUMS}))

# ===== Fig18: named-summary score heatmap (method x summary) =====
H=pd.DataFrame({s:sum_by_method(s) for s in SUMS}).reindex(C.METHODS)
fig,ax=plt.subplots(figsize=(6.8,4.7)); ax.grid(False)
im=ax.imshow(H.values,cmap=B.SCORE_CMAP,aspect="auto",vmin=0,vmax=1)
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        ax.text(j,i,f"{H.values[i,j]:.2f}",ha="center",va="center",fontsize=7.6,color=B.heat_text_color(H.values[i,j],0,1,B.SCORE_CMAP))
ax.set_xticks(range(4)); ax.set_xticklabels(["Balanced\nRealism","Spatial Tissue\nFidelity","Novelty-Gated","Conservative\nSoftmin"],fontsize=7.6)
ax.set_yticks(range(len(C.METHODS))); ax.set_yticklabels(C.METHODS,fontsize=9)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0); cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.outline.set_visible(False); cb.ax.set_title("score",fontsize=8,color="#6C7B95",loc="left")
ax.set_title("Candidate-summary scores")
C.save(fig,"fig_18_named_summaries",H.reset_index())

# ===== A15: coverage / gap map =====
cov=C.coverage()
fig,ax=plt.subplots(figsize=(6.2,4.7)); ax.grid(False)
for i,me in enumerate(C.METHODS):
    for j,dl in enumerate(C.DLABELS):
        if np.isnan(cov.loc[me,dl]): B.mark_not_run(ax,j,i)
        else:
            st=B.method_style(me)
            from matplotlib.patches import Rectangle
            ax.add_patch(Rectangle((j-0.5,i-0.5),1,1,facecolor=st["facecolor"] if st["facecolor"]!="white" else st["color"],
                                   edgecolor="white",lw=1.2,alpha=0.9,hatch=st["hatch"]))
ax.set_xlim(-0.5,len(C.DLABELS)-0.5); ax.set_ylim(len(C.METHODS)-0.5,-0.5)
ax.set_xticks(range(len(C.DLABELS))); ax.set_xticklabels(C.DLABELS); ax.set_yticks(range(len(C.METHODS))); ax.set_yticklabels(C.METHODS,fontsize=9)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0); ax.set_title("Run coverage")
B.legend_outside(ax,side="below",ncol=3,handles=[Line2D([0],[0],marker='s',ls='',mfc="#888",ms=11,label='run (29 pairs)'),
    Line2D([0],[0],marker='s',ls='',mfc=B.NOTRUN,ms=11,label='not run (gap)')])
C.save(fig,"fig_A15_coverage_map",cov.reset_index())

