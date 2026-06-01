"""Like-for-like check: scDesign3 against reference-free SRTsim on the three datasets where
both actually run. The pooled Pareto point averages over non-overlapping dataset mixes, so this
restricts to the shared pairs and shows the realism and novelty gaps per dataset."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B
SCD=B.SIM["scDesign3"]; SRT=B.SIM["SRTsim"]

nm=C.named()
br=nm[nm.summary=="Balanced Realism"].copy()
cs=nm[nm.summary=="copy_similarity"].copy()
key=["platform","simulator","dataset_variant"]
m=br.merge(cs[key+["score"]],on=key,suffixes=("_realism","_copy"))
m["novelty"]=1-m["score_copy"]; m["pair"]=m.platform+"-"+m.dataset_variant

scd=m[m.simulator=="scDesign3-rb"].groupby("pair").agg(realism=("score_realism","mean"),novelty=("novelty","mean"))
rf=m[m.simulator.isin(["SRTsim-rf-domain","SRTsim-rf-tissue"])].groupby("pair").agg(realism=("score_realism","mean"),novelty=("novelty","mean"))
shared=[p for p in ["MERFISH-Spapros","MERFISH-SpatialDE","Visium-Spapros"] if p in scd.index and p in rf.index]
y=np.arange(len(shared))[::-1]

fig,axes=plt.subplots(1,2,figsize=(9.0,4.2),sharey=True)
for ax,col,title,xlim in [(axes[0],"realism","Balanced Realism",(0.84,0.99)),
                          (axes[1],"novelty","Novelty (1 - copy-risk)",(-0.02,0.10))]:
    for i,p in zip(y,shared):
        a=rf.loc[p,col]; b=scd.loc[p,col]
        ax.plot([a,b],[i,i],color="#C8D0DA",lw=2.4,zorder=1)
        ax.scatter(a,i,s=130,facecolor="white",edgecolor=SRT,linewidth=2.0,zorder=3)
        ax.scatter(b,i,s=130,facecolor=SCD,edgecolor=SCD,linewidth=2.0,zorder=3)
    ax.set_xlim(*xlim); ax.set_title(title,fontsize=12,fontweight="bold",pad=10)
    ax.set_yticks(y); ax.set_yticklabels(shared)
    ax.grid(axis="x",color="#ECECEC",lw=0.8); ax.set_axisbelow(True)
    ax.margins(y=0.18)
axes[0].set_xlabel("score (higher = more realistic)")
axes[1].set_xlabel("score (higher = more novel)")
h=[Line2D([0],[0],marker="o",ls="",mfc=SCD,mec=SCD,ms=10,label="scDesign3 (rb)"),
   Line2D([0],[0],marker="o",ls="",mfc="white",mec=SRT,mew=2,ms=10,label="SRTsim (reference-free)")]
fig.legend(handles=h,loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(0.5,0.005),fontsize=9.5)
fig.subplots_adjust(left=0.17,right=0.975,top=0.88,bottom=0.27,wspace=0.10)
C.save(fig,"fig_32_shared_dataset_domination",None)
print("fig_32 rebuilt:",shared)
