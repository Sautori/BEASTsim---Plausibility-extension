"""Reference-based vs reference-free SRTsim similarity, by dataset."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, matplotlib.pyplot as plt
B=C.B
m=C.metric_long()
sim=m[m.family=="BEASTsim_similarity_DSC"]            # 3 Dice metrics, higher=better
sim=sim.groupby(["platform","dataset_variant","simulator"]).value.mean().reset_index()
srt=sim[sim.simulator.str.startswith("SRTsim")].copy()
srt["dl"]=[C.dlabel(p,d) for p,d in zip(srt.platform,srt.dataset_variant)]
# rb vs rf composite per dataset (mean of domain+tissue within mode). rf only exist non-Xenium
ds=[d for d in C.DLABELS if d!="Xenium"]
rb=[srt[(srt.dl==d)&(~srt.simulator.apply(B.is_reference_free))].value.mean() for d in ds]
rf=[srt[(srt.dl==d)&( srt.simulator.apply(B.is_reference_free))].value.mean() for d in ds]
col=B.SIM["SRTsim"]
fig,ax=plt.subplots(figsize=(7.4,4.6))
x=np.arange(len(ds))
for i,(a,b) in enumerate(zip(rb,rf)):
    win = b>a
    ax.plot([x[i],x[i]],[a,b],color=("#2A9D8F" if win else "#B0B0B0"),lw=3,zorder=1,solid_capstyle="round")
    ax.scatter(x[i],a,s=150,color=col,edgecolor=col,zorder=3)                 # rb solid
    ax.scatter(x[i],b,s=150,facecolor="white",edgecolor=col,linewidth=2.0,hatch="////",zorder=3) # rf open/hatched
    ax.annotate(f"{b-a:+.3f}",(x[i],max(a,b)),textcoords="offset points",xytext=(0,9),
                ha="center",fontsize=9,color=("#2A9D8F" if win else "#8A8A8A"),fontweight="bold")
# SpatialDE bracket (BEASTsim say rf recommended for SpatialDE-filtered)
for i,d in enumerate(ds):
    if "Spat" in d: ax.axvspan(x[i]-0.42,x[i]+0.42,color="#E7D7A8",alpha=0.35,zorder=0)
ax.set_ylim(0,1.05); ax.set_xticks(x); ax.set_xticklabels(ds)
ax.set_ylabel("similarity to real  (mean CT/SVG/Shape Dice →)")
ax.set_title("SRTsim: reference-based vs reference-free",loc="left")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=col,mec=col,ms=11,label='reference-based (rb)'),
                   Line2D([0],[0],marker='o',ls='',mfc='white',mec=col,mew=2,ms=11,label='reference-free (rf)'),
                   Line2D([0],[0],color='#2A9D8F',lw=3,label='rf closer to real'),
                   Line2D([0],[0],color='#E7D7A8',lw=6,alpha=0.9,label='SpatialDE-filtered (BEASTsim: rf rec.)')],
          loc="upper center",bbox_to_anchor=(0.5,-0.13),ncol=2,fontsize=8.5)
ax.text(0.0,-0.30,"Calibration comparison (similarity ≠ validated truth). rb = mean of domain+tissue; rf = mean of rf-domain+rf-tissue. Xenium omitted (rf not run).",
        transform=ax.transAxes,fontsize=7.2,color="#666",va="top")
src=__import__("pandas").DataFrame({"dataset":ds,"rb_similarity":rb,"rf_similarity":rf,"rf_minus_rb":np.array(rf)-np.array(rb)})
C.save(fig,"fig_R1_rb_vs_rf",src)
print(src.round(4).to_string(index=False))
