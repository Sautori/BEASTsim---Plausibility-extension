"""Further metric panels: SVG recovery, DE signal, components, family radar, sub-mode."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
B=C.B; m=C.metric_long(); nm=C.named(); ce=C.components()
RB,RF=B.REFMODE["rb"],B.REFMODE["rf"]
def mm(metric):
    d=m[m.metric==metric]
    return d.groupby("simulator").value.mean().reindex(C.METHODS)
def scatter_methods(ax,xv,yv):
    for me in C.METHODS:
        st=B.method_style(me)
        ax.scatter(xv[me],yv[me],s=160,facecolor=st["facecolor"],edgecolor=st["color"],hatch=st["hatch"],linewidth=1.7,zorder=3)

# Fig8: DE biological signal
x=mm("DE_lfc_correlation"); y=mm("DE_top_marker_overlap")
fig,ax=plt.subplots(figsize=(6.8,5.2)); scatter_methods(ax,x,y)
ax.set_xlabel("DE log-FC correlation"); B.ylabel_h(ax,"top-marker overlap"); ax.set_title("DE signal")
B.legend_outside(ax,side="below",ncol=4,handles=B.variant_legend_handles())
C.save(fig,"fig_08_de_signal",pd.DataFrame({"method":C.METHODS,"lfc_corr":x.values,"marker_overlap":y.values}))

# Fig6: SVG recovery precision-recall
x=mm("SVG_recall"); y=mm("SVG_precision"); jac=mm("SVG_jaccard")
fig,ax=plt.subplots(figsize=(6.8,5.2))
for me in C.METHODS:
    st=B.method_style(me)
    ax.scatter(x[me],y[me],s=80+700*float(np.nan_to_num(jac[me])),facecolor=st["facecolor"],edgecolor=st["color"],hatch=st["hatch"],linewidth=1.7,alpha=0.9,zorder=3)
ax.set_xlabel("SVG recall"); B.ylabel_h(ax,"SVG precision"); ax.set_title("SVG recovery")
sl=[0.1,0.3,0.5]
size_handles=[Line2D([0],[0],marker="o",linestyle="",markerfacecolor="#D5DBE3",markeredgecolor="#6C7B95",
              markersize=np.sqrt(80+700*j)/2.2,label=f"{j:.1f}") for j in sl]
leg_sz=ax.legend(handles=size_handles,title="Jaccard",loc="upper left",bbox_to_anchor=(1.01,1.0),frameon=False,
                 fontsize=8,title_fontsize=8,labelspacing=1.2,borderpad=0.8,handletextpad=1.1)
ax.add_artist(leg_sz)
B.legend_outside(ax,side="below",ncol=4,handles=B.variant_legend_handles())
C.save(fig,"fig_06_svg_recovery",pd.DataFrame({"method":C.METHODS,"recall":x.values,"precision":y.values,"jaccard":jac.values}))


# Fig17: component-evidence heatmap
comps=["expression","spatial_gene","layout","proximity_composition"]
H=pd.DataFrame({c:ce[ce.component==c].groupby("simulator").score.mean() for c in comps}).reindex(C.METHODS)
fig,ax=plt.subplots(figsize=(6.8,4.7)); ax.grid(False)
im=ax.imshow(H.values,cmap=B.SCORE_CMAP,aspect="auto",vmin=0,vmax=1)
for i in range(H.shape[0]):
    for j in range(H.shape[1]): ax.text(j,i,f"{H.values[i,j]:.2f}",ha="center",va="center",fontsize=7.6,color=B.heat_text_color(H.values[i,j],0,1,B.SCORE_CMAP))
ax.set_xticks(range(4)); ax.set_xticklabels(["expression","spatial\ngene","layout","proximity\ncomposition"],fontsize=8)
ax.set_yticks(range(len(C.METHODS))); ax.set_yticklabels(C.METHODS,fontsize=9)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0); cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.outline.set_visible(False); cb.ax.set_title("score",fontsize=8,color="#6C7B95",loc="left")
ax.set_title("Component evidence")
C.save(fig,"fig_17_component_evidence",H.reset_index())

# R2: SRTsim domain vs tissue x rb vs rf (similarity), grouped bar
sim=m[m.family=="BEASTsim_similarity_DSC"].groupby("simulator").value.mean()
srt=["SRTsim-rb-domain","SRTsim-rb-tissue","SRTsim-rf-domain","SRTsim-rf-tissue"]
fig,ax=plt.subplots(figsize=(6.8,4.8)); xb=np.arange(len(srt))
base=B.SIM["SRTsim"]; isrf=[False,False,True,True]
for i,me in enumerate(srt):
    c=B._lighten(base,0.50) if isrf[i] else base
    ax.bar(xb[i],sim[me],0.66,color=c,edgecolor=base,linewidth=1.4)
ax.set_xticks(xb); ax.set_xticklabels(["rb domain","rb tissue","rf domain","rf tissue"],fontsize=9); ax.grid(axis="x",visible=False)
B.ylabel_h(ax,"similarity to real (mean Dice)"); ax.set_title("SRTsim mode and sub-mode")
B.legend_outside(ax,side="below",ncol=2,handles=[Rectangle((0,0),1,1,facecolor=base,label='reference-based'),
    Rectangle((0,0),1,1,facecolor=B._lighten(base,0.50),edgecolor=base,label='reference-free')])
C.save(fig,"fig_R2_mode_submode",pd.DataFrame({"variant":srt,"similarity":[sim[v] for v in srt]}))

# Fig24: family-representative radar, the best variant of each families by mean similarity
fams={"SRTsim":["SRTsim-rb-domain","SRTsim-rb-tissue","SRTsim-rf-domain","SRTsim-rf-tissue"],
      "scCube":["scCube-rfb","scCube-rb"],"scDesign3":["scDesign3-rb"],"SpatialcoGCN":["SpatialcoGCN-rb"]}
rep={f:max(v,key=lambda me:sim.get(me,0)) for f,v in fams.items()}
axes_def=[("Data props",lambda d:1-d[d.family=="E1_data_properties"].value.mean()),
          ("Similarity",lambda d:d[d.family=="BEASTsim_similarity_DSC"].value.mean()),
          ("DE signal",lambda d:np.clip(d[d.metric=="DE_lfc_correlation"].value.mean(),0,1)),
          ("SVG",lambda d:d[d.metric=="SVG_jaccard"].value.mean()),
          ("Downstream",lambda d:np.clip(d[d.metric=="cluster_ARI_real_vs_sim"].value.mean(),0,1)),
          ("CT prop",lambda d:1-d[d.metric=="ct_proportion_JSD"].value.mean()),
          ("Novelty\n(trade-off)",lambda d:1-np.clip(d[d.metric=="copy_similarity"].value.mean(),0,1))]
labels=[a for a,_ in axes_def]; N=len(labels); ang=np.linspace(0,2*np.pi,N,endpoint=False).tolist(); ang+=ang[:1]
fig,ax=plt.subplots(figsize=(6.4,6.2),subplot_kw=dict(polar=True))
for f,me in rep.items():
    sub=m[m.simulator==me]; v=[fn(sub) for _,fn in axes_def]; v+=v[:1]
    ax.plot(ang,v,color=B.SIM[f],lw=2.2,label=f"{f} ({me.split('-',1)[1]})"); ax.fill(ang,v,color=B.SIM[f],alpha=0.07)
ax.set_xticks(ang[:-1]); ax.set_xticklabels(labels,fontsize=9); ax.set_ylim(0,1); ax.set_yticks([0.25,0.5,0.75,1.0]); ax.set_yticklabels(["","0.5","","1.0"],fontsize=7.5)
ax.set_title("Family fingerprints",pad=18)
B.legend_outside(ax,side="below",ncol=2)
C.save(fig,"fig_24_family_fingerprints",pd.DataFrame({rep[f]:[fn(m[m.simulator==rep[f]]) for _,fn in axes_def] for f in rep},index=labels).reset_index())
print("extras done")
