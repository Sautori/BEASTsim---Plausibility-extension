"""Top-ranked simulator as the gate threshold sweeps."""
import sys; sys.path.insert(0,".")
import common as C, numpy as np, matplotlib.pyplot as plt
B=C.B; nm=C.named()
bal=nm[nm.summary=="Balanced Realism"].groupby("simulator").score.mean()
copy=nm[nm.summary=="copy_similarity"].groupby("simulator").score.mean()
W=0.05; tau=np.linspace(0.90,1.0,400)
def ng(me): return bal[me]*np.clip(1-np.maximum(copy[me]-tau,0)/W,0,1)
cands=["scCube-rfb","scDesign3-rb","SRTsim-rb-domain"]
curves={me:ng(me) for me in cands}
winner=np.array([max(cands,key=lambda me:curves[me][i]) for i in range(len(tau))])
fig,ax=plt.subplots(figsize=(7.8,4.8))
# faint winner-regoin bands
for i in range(len(tau)-1):
    ax.axvspan(tau[i],tau[i+1],color=B.SIM[B.family_of(winner[i])],alpha=0.06,zorder=0,lw=0)
ax.axvline(0.95,color=B.ALERT,ls=":",lw=1.2,zorder=1); ax.text(0.95,1.04,"default 0.95",fontsize=8,ha="center",color=B.ALERT)
for me in cands:
    st=B.method_style(me)
    ax.plot(tau,curves[me],color=st["color"],lw=1.3,alpha=0.30,zorder=2)
    yv=np.where(winner==me,curves[me],np.nan)
    ax.plot(tau,yv,color=st["color"],lw=4.0,alpha=0.95,zorder=3,solid_capstyle="round")
    ax.text(1.002,curves[me][-1],me,color=st["color"],fontsize=8.6,va="center",fontweight="bold")
# crossover, labels placed low in clear space with leaders
flips=[i for i in range(1,len(winner)) if winner[i]!=winner[i-1]]
for k,i in enumerate(flips):
    xv,yv=tau[i],curves[winner[i]][i]
    ax.scatter(xv,yv,s=45,color="#22303C",zorder=4)
    ax.annotate(f"τ≈{xv:.2f}",(xv,yv),xytext=(9,-14),textcoords="offset points",fontsize=8,ha="left",color="#22303C",
                arrowprops=dict(arrowstyle="-",color="#C2C8D0",lw=0.8))
ax.set_xlabel("novelty-gate threshold τ"); B.ylabel_h(ax,"Novelty-Gated Plausibility (leader bold)")
ax.set_title("Recommendation vs gate limit"); ax.set_xlim(0.90,1.045); ax.set_ylim(0,1.08)
C.save(fig,"fig_27_reco_calibration",None)
print("Fig27 modernized | x starts 0.90 | flips:",[round(tau[i],3) for i in flips])
