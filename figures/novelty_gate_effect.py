"""Balanced Realism to Novelty-Gated Plausibility per method-variant."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
B=C.B; nm=C.named()
br=nm[nm.summary=="Balanced Realism"].groupby("simulator").score.mean().reindex(C.METHODS)
ng=nm[nm.summary=="Novelty-Gated Plausibility"].groupby("simulator").score.mean().reindex(C.METHODS)
df=pd.DataFrame({"br":br,"ng":ng}); df["drop"]=df.br-df.ng
df=df.sort_values("drop")            # smallest drop bottom, largest top
fig,ax=plt.subplots(figsize=(7.4,4.2))
yy=np.arange(len(df))
for i,(me,r) in enumerate(df.iterrows()):
    st=B.method_style(me); fam=st["color"]
    if r["drop"]>0.01:
        xs=np.linspace(r.ng,r.br,80); pts=np.array([xs,np.full_like(xs,i)]).T.reshape(-1,1,2)
        segs=np.concatenate([pts[:-1],pts[1:]],axis=1)
        cmap=LinearSegmentedColormap.from_list("g",[B.ALERT,fam])
        lc=LineCollection(segs,cmap=cmap,array=np.linspace(0,1,len(segs)),lw=6,capstyle="round",zorder=1)
        ax.add_collection(lc)
        ax.scatter(r.br,i,s=150,color=fam,edgecolor="white",linewidth=1.2,zorder=3)
        ax.scatter(r.ng,i,s=150,facecolor="white",edgecolor=B.ALERT,linewidth=2.2,zorder=3)
    else:   # gate did not fire: BR==NG, concentric marker
        ax.scatter(r.br,i,s=150,color=fam,edgecolor="white",linewidth=1.2,zorder=3)
        ax.scatter(r.br,i,s=360,facecolor="none",edgecolor=B.NEUTRAL_MID,linewidth=2.2,marker="s",zorder=4)
        ax.annotate("gate open",(r.br,i),xytext=(11,0),textcoords="offset points",va="center",fontsize=7.5,color=B.NEUTRAL_MID)
ax.set_yticks(yy); ax.set_yticklabels(df.index,fontsize=9); ax.grid(axis="y",visible=False)
ax.set_xlim(-0.03,1.06); ax.set_xlabel("plausibility score"); ax.set_title("Novelty-gate effect")
B.legend_outside(ax,side="below",ncol=3,handles=[
    Line2D([0],[0],marker='o',ls='',mfc='#3A4654',mec='white',ms=11,label='Balanced Realism (filled)'),
    Line2D([0],[0],marker='o',ls='',mfc='white',mec=B.ALERT,mew=2,ms=11,label='Novelty-Gated (open)'),
    Line2D([0],[0],marker='s',ls='',mfc='none',mec=B.NEUTRAL_MID,mew=2,ms=13,label='gate did not fire')])
C.save(fig,"fig_20_novelty_gate_effect",df.reset_index())
print("Fig20 redesigned")
