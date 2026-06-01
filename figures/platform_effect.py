"""scDesign3 task metrics across the three platforms."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B; PL=["MERFISH","Visium","Xenium"]
m=C.metric_long(); d=m[m.simulator=="scDesign3-rb"]
# comparable, higher-is-better metric across platforms
specs=[("SVG\nJaccard","SVG_jaccard",lambda v:v),
       ("DE LFC\ncorrelation","DE_lfc_correlation",lambda v:v),
       ("CT domain\nDice","CT_mean_dice_score",lambda v:v),
       ("CT proportion\n(1 \u2212 JSD)","ct_proportion_JSD",lambda v:1-v),
       ("Moran profile\ncosine","moran_profile_cosine",lambda v:v),
       ("Clustering\nARI","cluster_ARI_real_vs_sim",lambda v:np.clip(v,0,1))]
vals={pl:[] for pl in PL}
for lab,met,fn in specs:
    for pl in PL:
        sub=d[(d.metric==met)&(d.platform==pl)].value
        vals[pl].append(fn(sub.mean()) if len(sub) else np.nan)
fig,ax=plt.subplots(figsize=(8.8,4.8)); ax.grid(axis="x",visible=False)
xb=np.arange(len(specs)); w=0.26
for k,pl in enumerate(PL):
    ax.bar(xb+(k-1)*w,vals[pl],w,color=B.PLATFORM[pl],edgecolor="white",linewidth=0.6,label=pl)
ax.set_xticks(xb); ax.set_xticklabels([s[0] for s in specs],fontsize=8.4)
ax.set_ylim(0,1.02); B.ylabel_h(ax,"score (higher = closer to real)")
ax.set_title("Platform effect (scDesign3)")
B.legend_outside(ax,side="below",ncol=3,handles=[Line2D([0],[0],marker='s',ls='',mfc=B.PLATFORM[p],mec='white',ms=11,label=p) for p in PL])
C.save(fig,"fig_29_platform_compare",None)
print("fig_29 platform comparison built")
