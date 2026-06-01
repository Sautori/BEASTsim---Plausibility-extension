"""Realism against novelty, with the Pareto front marked."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B; nm=C.named()
real=nm[nm.summary=="Balanced Realism"].groupby("simulator").score.mean().reindex(C.METHODS)
copy=nm[nm.summary=="copy_similarity"].groupby("simulator").score.mean().reindex(C.METHODS)
nov=1-copy
def dom(a,b): return (a[0]>=b[0] and a[1]>=b[1]) and (a[0]>b[0] or a[1]>b[1])
pts={me:(real[me],nov[me]) for me in C.METHODS}
rem=dict(pts); front={}; f=1
while rem:
    nd=[me for me in rem if not any(dom(rem[o],rem[me]) for o in rem if o!=me)]
    for me in nd: front[me]=f
    for me in nd: del rem[me]
    f+=1

order=sorted(C.METHODS,key=lambda me:(round(real[me],2),round(nov[me],2)))
groups=[]
for me in order:
    placed=False
    for g in groups:
        if abs(real[me]-np.mean([real[x] for x in g]))<0.02 and abs(nov[me]-np.mean([nov[x] for x in g]))<0.02:
            g.append(me); placed=True; break
    if not placed: groups.append([me])
pos={}
for g in groups:
    n=len(g)
    for k,me in enumerate(g):
        pos[me]=(real[me], nov[me]+(k-(n-1)/2)*0.022)

fig,ax=plt.subplots(figsize=(7.0,5.4))
f1=sorted([me for me in C.METHODS if front[me]==1],key=lambda me:real[me])
ax.plot([real[me] for me in f1],[nov[me] for me in f1],color=B.NEUTRAL_MID,lw=1.6,ls="--",zorder=1)
for me in C.METHODS:
    st=B.method_style(me); on=front[me]==1; x,y=pos[me]
    ax.scatter(x,y,s=230 if on else 140,facecolor=(st["facecolor"] if on else "none"),edgecolor=st["color"],
               hatch=(st["hatch"] if on else None),linewidth=(2.2 if on else 2.0),alpha=0.95,zorder=3 if on else 2)
# orientation annotations in clear space
ax.annotate("most novel",(real["SpatialcoGCN-rb"],nov["SpatialcoGCN-rb"]),xytext=(12,-2),textcoords="offset points",fontsize=8,color="#6C7B95")
ax.annotate("most realistic",(real["SRTsim-rb-domain"],pos["SRTsim-rb-domain"][1]),xytext=(-4,-16),textcoords="offset points",fontsize=8,color="#6C7B95",ha="right")
ax.set_xlabel("realism  (Balanced Realism →)"); B.ylabel_h(ax,"novelty  (1 - copy similarity →)")
ax.set_title("Realism-novelty trade-off"); ax.set_xlim(0.35,1.02); ax.set_ylim(-0.05,0.66)
H=[Line2D([0],[0],marker='o',ls='',mfc=B.SIM[f],mec=B.SIM[f],ms=11,label=f) for f in B.FAMILIES]
H+=[Line2D([0],[0],marker='o',ls='',mfc='#888',mec='#22303C',mew=2,ms=12,label='Pareto-optimal (filled)'),
    Line2D([0],[0],marker='o',ls='',mfc='none',mec='#888',mew=2,ms=11,label='dominated (hollow)')]
B.legend_outside(ax,side="below",ncol=3,handles=H)
C.save(fig,"fig_16_pareto_tradeoff",pd.DataFrame({"method":C.METHODS,"realism":real.values,"novelty":nov.values,"front":[front[me] for me in C.METHODS]}))
print("front1:",sorted(f1))
