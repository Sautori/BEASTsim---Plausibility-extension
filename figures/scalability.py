"""Compute time against peak memory, per method-variant."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=C.B; m=C.metric_long()
t=m[m.metric=="compute_seconds"].groupby("simulator").value.mean().reindex(C.METHODS)
mem=m[m.metric=="peak_memory_gb"].groupby("simulator").value.mean().reindex(C.METHODS)
pts=[(t[me],mem[me],me) for me in C.METHODS]
xr=t.max()-t.min(); yr=mem.max()-mem.min(); tolx,toly=0.05*xr,0.07*yr
clusters=[]
for x,y,me in pts:
    for cl in clusters:
        cx,cy=np.mean([p[0] for p in cl]),np.mean([p[1] for p in cl])
        if abs(x-cx)<=tolx and abs(y-cy)<=toly: cl.append((x,y,me)); break
    else: clusters.append([(x,y,me)])
fig,ax=plt.subplots(figsize=(7.4,5.0))
# moden: faint "ideal" corner (fast & light = lower-left), soft directional hint
ax.axhspan(mem.min()-0.2*yr, mem.min()+0.35*yr, xmin=0, xmax=0.34, color="#2A9D8F", alpha=0.05, zorder=0)
ax.annotate("faster · lighter", xy=(t.min(), mem.min()), xytext=(t.min()+0.04*xr, mem.min()+0.02*yr),
            fontsize=8.5, color="#2A9D8F", style="italic", va="bottom")
def lab_for(cl):
    if len(cl)==1: return cl[0][2]
    base=cl[0][2].rsplit("-",1)[0]; return base+" (domain/tissue)"
dy=0.06*yr
for cl in clusters:
    n=len(cl); cx=np.mean([p[0] for p in cl]); cy=np.mean([p[1] for p in cl])
    if n>1:
        ax.plot([cx,cx],[cy-(n-1)/2*dy,cy+(n-1)/2*dy],color=B.NEUTRAL_LIGHT,lw=1.4,zorder=1)
        for k,(x,y,me) in enumerate(sorted(cl,key=lambda p:p[2])):
            st=B.method_style(me); yy=cy+(k-(n-1)/2)*dy
            ax.scatter(cx,yy,s=170,facecolor=st["facecolor"],edgecolor=st["color"],hatch=st["hatch"],linewidth=1.8,zorder=3)
    else:
        x,y,me=cl[0]; st=B.method_style(me)
        ax.scatter(x,y,s=170,facecolor=st["facecolor"],edgecolor=st["color"],hatch=st["hatch"],linewidth=1.8,zorder=3)
    fam=B.SIM[B.family_of(cl[0][2])]
    ax.annotate(lab_for(cl),(cx,cy),xytext=(13,0),textcoords="offset points",va="center",
                fontsize=8.2,color=fam,fontweight="bold")
ax.set_xlim(t.min()-0.08*xr, t.max()+0.30*xr); ax.set_ylim(mem.min()-0.18*yr, mem.max()+0.16*yr)
ax.set_xlabel("compute seconds"); B.ylabel_h(ax,"peak memory (GB)"); ax.set_title("Scalability")
B.legend_outside(ax,side="below",ncol=2,handles=[
    Line2D([0],[0],marker='o',ls='',mfc="#888",mec="#888",ms=11,label='reference-based (solid)'),
    Line2D([0],[0],marker='o',ls='',mfc='white',mec="#888",mew=2,ms=11,label='reference-free (hatched/open)')])
C.save(fig,"fig_10_scalability",pd.DataFrame({"method":C.METHODS,"seconds":t.values,"mem_gb":mem.values}))
print("Fig10 modernized")
