"""Core metric panels: data-property realism, similarity Dice, downstream, copy-risk."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B; m=C.metric_long(); nm=C.named()
def mean_metric(fam=None, metric=None):
    d=m.copy()
    if fam: d=d[d.family==fam]
    if metric: d=d[d.metric==metric]
    return d.groupby(["platform","dataset_variant","simulator"]).value.mean().reset_index()
def by_method(df,col="value"): return df.groupby("simulator")[col].mean().reindex(C.METHODS)
RB,RF=B.REFMODE["rb"],B.REFMODE["rf"]

# ---- R1: rb vs rf similarity dumbell (purple=rb, gold=rf) ----
sim=mean_metric(fam="BEASTsim_similarity_DSC").rename(columns={"value":"sim"})
srt=sim[sim.simulator.str.startswith("SRTsim")].copy(); srt["dl"]=[C.dlabel(p,d) for p,d in zip(srt.platform,srt.dataset_variant)]
ds=[d for d in C.DLABELS if d!="Xenium"]
rb=[srt[(srt.dl==d)&(~srt.simulator.apply(B.is_reference_free))].sim.mean() for d in ds]
rf=[srt[(srt.dl==d)&( srt.simulator.apply(B.is_reference_free))].sim.mean() for d in ds]
x=np.arange(len(ds))
fig,ax=plt.subplots(figsize=(7.0,4.7))
for i,(a,b) in enumerate(zip(rb,rf)):
    if "Spat" in ds[i]: ax.axvspan(x[i]-0.42,x[i]+0.42,color=B.NEUTRAL_LIGHT,alpha=0.55,zorder=0)
    ax.plot([x[i],x[i]],[a,b],color="#C2C8D0",lw=3.4,zorder=1,solid_capstyle="round")
    ax.scatter(x[i],a,s=210,color=RB,edgecolor="white",linewidth=1.4,zorder=3)
    ax.scatter(x[i],b,s=210,color=RF,edgecolor="white",linewidth=1.4,zorder=3)
ax.set_xticks(x); ax.set_xticklabels(ds); B.ylabel_h(ax,"similarity to real (mean Dice)")
ax.set_title("SRTsim: rb vs rf")
B.legend_outside(ax,side="below",ncol=3,handles=[
    Line2D([0],[0],marker='o',ls='',mfc=RB,mec='white',ms=12,label='reference-based'),
    Line2D([0],[0],marker='o',ls='',mfc=RF,mec='white',ms=12,label='reference-free'),
    Line2D([0],[0],marker='s',ls='',mfc=B.NEUTRAL_LIGHT,ms=13,label='SpatialDE-filtered')])
C.save(fig,"fig_R1_rb_vs_rf",pd.DataFrame({"dataset":ds,"rb":rb,"rf":rf}))

# ---- R1b: similarity vs copy-risk (family hue, cross-family) ----
cs=nm[nm.summary=="copy_similarity"][["platform","dataset_variant","simulator","score"]].rename(columns={"score":"copy"})
d=sim.merge(cs,on=["platform","dataset_variant","simulator"])
fig,ax=plt.subplots(figsize=(6.8,5.2))

# ---- Fig5: data-property realism heatmap (clean: solid gaps, horizontal cbar lable) ----
e1=mean_metric(fam="E1_data_properties"); e1["dl"]=[C.dlabel(p,v) for p,v in zip(e1.platform,e1.dataset_variant)]
P=e1.pivot_table(index="simulator",columns="dl",values="value").reindex(index=C.METHODS,columns=C.DLABELS)
fig,ax=plt.subplots(figsize=(6.9,4.7)); vmax=max(0.4,np.nanmax(P.values)); ax.grid(False)
im=ax.imshow(np.ma.masked_invalid(P.values),cmap=B.DIST_CMAP,aspect="auto",vmin=0,vmax=vmax)
for i in range(P.shape[0]):
    for j in range(P.shape[1]):
        v=P.values[i,j]
        if np.isnan(v): B.mark_not_run(ax,j,i)
        else: ax.text(j,i,f"{v:.02f}",ha="center",va="center",fontsize=7.6,color=B.heat_text_color(v,0,vmax,B.DIST_CMAP))
ax.set_xticks(range(P.shape[1])); ax.set_xticklabels(P.columns); ax.set_yticks(range(P.shape[0])); ax.set_yticklabels(P.index,fontsize=9)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.outline.set_visible(False)
cb.ax.set_title("mean KS\n(lower = closer)",fontsize=8,color="#6C7B95",loc="left")
ax.set_title("Data-property realism")
C.save(fig,"fig_05_data_property_realism",P.reset_index())

# ---- Fig11: similarity Dice components ----
comp={"CT":"CT_mean_dice_score","SVG":"SVG_mean_dice_score","Shape":"Shape_mean_dice_score"}
rows=[(me,lab,by_method(mean_metric(metric=met))[me]) for lab,met in comp.items() for me in C.METHODS]
dd=pd.DataFrame(rows,columns=["method","component","val"])
fig,ax=plt.subplots(figsize=(8.4,4.8)); xb=np.arange(len(C.METHODS)); w=0.26
for k,lab in enumerate(comp):
    vals=[dd[(dd.method==me)&(dd.component==lab)].val.values[0] for me in C.METHODS]
    cols=[B.method_style(me)["color"] for me in C.METHODS]
    ax.bar(xb+(k-1)*w,vals,w,color=cols,alpha=[1.0,0.62,0.34][k],edgecolor="white",linewidth=0.6,label=lab)
ax.set_xticks(xb); ax.set_xticklabels(C.METHODS,rotation=24,ha="right",fontsize=8.5); ax.grid(axis="x",visible=False)
B.ylabel_h(ax,"mean Dice (higher = closer)"); ax.set_title("Similarity components")
B.legend_outside(ax,side="below",ncol=3,title="Dice")
C.save(fig,"fig_11_similarity_dsc",dd)

# ---- Fig9: downstream preservation ----
ari=by_method(mean_metric(metric="cluster_ARI_real_vs_sim")); jsd=by_method(mean_metric(metric="ct_proportion_JSD"))
fig,ax=plt.subplots(figsize=(6.8,5.2))
for me in C.METHODS:
    st=B.method_style(me)
    ax.scatter(ari[me],1-jsd[me],s=160,facecolor=st["facecolor"],edgecolor=st["color"],hatch=st["hatch"],linewidth=1.7,zorder=3)
ax.set_xlabel("clustering ARI (real vs simulated)"); B.ylabel_h(ax,"cell-type proportion (1 - JSD)")
ax.set_title("Downstream preservation")
B.legend_outside(ax,side="below",ncol=4,handles=B.variant_legend_handles())
C.save(fig,"fig_09_downstream_tasks",pd.DataFrame({"method":C.METHODS,"ARI":ari.values,"one_minus_JSD":(1-jsd).values}))

# ---- Fig14: copy-risk calibration lollipop ----
cps=nm[nm.summary=="copy_similarity"].groupby("simulator").score.mean().reindex(C.METHODS); order=cps.sort_values().index.tolist()
fig,ax=plt.subplots(figsize=(7.0,4.5)); ax.grid(axis="x",visible=False)
ax.axvline(0.95,color="#B0182B",ls="--",lw=1.2,zorder=1)
ax.text(0.95,len(order)-0.4,"0.95",fontsize=8,color="#B0182B",ha="center",va="bottom")
for i,me in enumerate(order):
    st=B.method_style(me)
    ax.plot([0,cps[me]],[i,i],color=st["color"],lw=2.4,zorder=2,alpha=0.85)
    ax.scatter(cps[me],i,s=150,facecolor=st["facecolor"],edgecolor=st["color"],hatch=st["hatch"],linewidth=1.7,zorder=3)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order,fontsize=8.5); ax.set_xlabel("copy similarity (→ reference)")
ax.set_title("Copy-risk calibration")
C.save(fig,"fig_14_copy_risk",cps.reset_index())
print("core re-rendered (locked style)")
