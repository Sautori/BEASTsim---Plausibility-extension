"""Plotting palette and small style helpers. Simulator families use Okabe-Ito colours,
platforms a muted accent, and coverage gaps a hatched grey. matplotlib and seaborn only."""
import matplotlib as mpl, numpy as np

INK="#22303C"; TEAL="#2A9D8F"; GRID="#E6E6E6"; SLATE="#6C7B95"; CMAP="plasma"
# BEASTsim spectral gradient: dark-blue -> magenta/red -> orange -> yellow (tissue / cell-type ramp)
from matplotlib.colors import LinearSegmentedColormap as _LSC
BEAST_CMAP=_LSC.from_list("beast",["#15306B","#3B2A78","#7A1F57","#B81D34","#E8651B","#F39B12","#F7E018"])
# simulator FAMILY = primary colour (Okabe-Ito, colourblind-safe)
SIM={"real":"#2B2B2B","SRTsim":"#0072B2","scCube":"#D55E00","scDesign3":"#009E73","SpatialcoGCN":"#E69F00"}
# platform = secondary accent + marker shape (never one data marks)
PLATFORM={"MERFISH":"#5E81AC","Visium":"#BF616A","Xenium":"#8FBCBB"}
SHAPE={"MERFISH":"o","Visium":"^","Xenium":"s"}
# reference-mode palette (plasma-derived) for rb-vs-rf-centric figures: purple=rb, gold=rf
REFMODE={"rb":"#6A1B9A","rf":"#F4C20D"}
# RESERVED ROLE COLORS (one colour = one meaning, used consistantly across all figures):
ALERT="#B0182B"          # risk / threshold / problem: copy-risk zone, gate line, blind-spot increase
NEUTRAL_DARK="#3A4654"   # structural dark / "expected/good" direction
NEUTRAL_MID="#7E8A9C"    # structural mid
NEUTRAL_LIGHT="#C8D0DA"  # structural light / "flat" / connectors
SCORE_CMAP="plasma"      # higher-is-better score heatmaps
DIST_CMAP="plasma_r"     # lower-is-better distance heatmaps
# named-summary greyscale ladder (neutral, never clashes with family/refmode). The odd-one-out uses ALERT
SUMMARY_NEUTRAL={"Balanced Realism":"#B8C0CC","Spatial Tissue Fidelity":"#7E8A9C","Novelty-Gated Plausibility":ALERT}
# coverage gaps (method x dataset not run): light hatched grey, never imputed
NOTRUN="#DDDDDD"; NOTRUN_HATCH="////"; NOTRUN_EDGE="#BBBBBB"

# canonical method-variant order (for axes / legends)
METHOD_ORDER=["SRTsim-rb-domain","SRTsim-rb-tissue","SRTsim-rf-domain","SRTsim-rf-tissue",
              "scCube-rfb","scCube-rb","scDesign3-rb","SpatialcoGCN-rb"]
FAMILIES=["SRTsim","scCube","scDesign3","SpatialcoGCN"]

def family_of(method):
    m=str(method)
    for f in FAMILIES:
        if m.startswith(f): return f
    return "real" if m.lower().startswith("real") else m.split("-")[0]

def is_reference_free(method):
    """rf or rfb (reference-free / reference-free-bulk) -> True, rb -> False."""
    m=str(method).lower()
    return ("-rf" in m) or ("rfb" in m)

def submode(method):
    m=str(method).lower()
    if "domain" in m: return "domain"
    if "tissue" in m: return "tissue"
    return ""

def _lighten(hex_c, f=0.0):
    """blend hex toward white by fraction f in [0,1]."""
    h=hex_c.lstrip("#"); r,g,b=(int(h[i:i+2],16) for i in (0,2,4))
    r,g,b=(int(c+(255-c)*f) for c in (r,g,b))
    return "#%02x%02x%02x"%(r,g,b)

def method_style(method):
    """Return a dict of plot kwargs encoding family (hue) + reference mode (fill/line)
    + sub-mode (value-step), legible in greyscale. Use:
        st=method_style(m); ax.bar(...,color=st['facecolor'],edgecolor=st['color'],hatch=st['hatch'])
        ax.plot(...,color=st['color'],ls=st['ls'],marker=st['marker'],mfc=st['mfc'])"""
    fam=family_of(method); base=SIM.get(fam,"#888888")
    rf=is_reference_free(method); sub=submode(method)
    step={"domain":0.0,"tissue":0.28,"":0.0}[sub]          # tissue = lighter value-step of same hue
    color=_lighten(base,step)
    return dict(family=fam, color=color, base=base, ref="rf" if rf else "rb", submode=sub,
                facecolor=("white" if rf else color),     # rf = open/hatched, rb = solid
                hatch=("////" if rf else None),
                ls=("--" if rf else "-"),
                marker=("o"),
                mfc=("white" if rf else color),            # open marker for rf
                edgecolor=color, lw=1.6, alpha=0.92)

def variant_legend_handles():
    """Family-hue + rb/rf fill key, as a list of (label, Patch) for one shared legend."""
    from matplotlib.patches import Patch
    H=[Patch(facecolor=SIM[f], edgecolor=SIM[f], label=f) for f in FAMILIES]
    H+=[Patch(facecolor="#777777", edgecolor="#777777", label="reference-based (rb)"),
        Patch(facecolor="white", edgecolor="#777777", hatch="////", label="reference-free (rf/rfb)")]
    return H

def mark_not_run(ax, x, y, w=1.0, h=1.0, hatch=False):
    """Overlay the greyed not-run convention on a grid/heatmap cell centred at (x,y).
    Solid light grey by default (clean). Set hatch=True only where a legend needs it."""
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((x-w/2,y-h/2),w,h,facecolor=NOTRUN,edgecolor="white",
                           hatch=(NOTRUN_HATCH if hatch else None),lw=0.8,zorder=3))

def ylabel_h(ax, text, color="#6C7B95", fontsize=9.5):
    """Horizontal y-axis label placed above the axis (never rotated vertical)."""
    ax.set_ylabel("")
    ax.text(0.0,1.012,text,transform=ax.transAxes,ha="left",va="bottom",fontsize=fontsize,color=color)

def apply(font=11):
    import seaborn as sns
    sns.set_theme(style="whitegrid", context="notebook", font="DejaVu Sans")
    mpl.rcParams.update({
        "figure.facecolor":"white","axes.facecolor":"white","savefig.facecolor":"white",
        "font.size":font,"axes.titlesize":font+3,"axes.titleweight":"bold","axes.titlecolor":"#1F2A36",
        "axes.titlelocation":"left","axes.titlepad":24,
        "axes.labelcolor":"#3A4654","axes.labelsize":font,"text.color":"#3A4654",
        "xtick.color":"#6C7B95","ytick.color":"#6C7B95","xtick.labelsize":font-1.5,"ytick.labelsize":font-1.5,
        "axes.edgecolor":"#CBD2DA","axes.linewidth":1.0,
        "axes.spines.top":False,"axes.spines.right":False,
        "grid.color":"#EBEEF2","grid.linewidth":0.9,"axes.grid.axis":"y",
        "legend.frameon":False,"legend.fontsize":font-1.5,
        "pdf.fonttype":42,"savefig.bbox":"tight","savefig.dpi":300})

def legend_outside(ax, handles=None, labels=None, side="right", title=None, ncol=1, fontsize=9):
    """Place a frameless legend OUTSIDE the axes (right by default, or below)."""
    kw=dict(frameon=False,fontsize=fontsize,title=title,ncol=ncol)
    if side=="right":
        leg=ax.legend(handles=handles,labels=labels,loc="upper left",bbox_to_anchor=(1.02,1.0),**kw)
    else:
        leg=ax.legend(handles=handles,labels=labels,loc="upper center",bbox_to_anchor=(0.5,-0.16),**kw)
    if title: leg.get_title().set_fontsize(fontsize)
    return leg

def celltype_cmap(labels):
    """Abundance-ordered cell types sampled along the BEASTsim spectral gradient
    (dark-blue -> red/orange -> yellow). Identical mapping across all tissue panels."""
    order=[c for c,_ in sorted(zip(*np.unique(labels,return_counts=True)),key=lambda t:-t[1])]
    import matplotlib.pyplot as _plt
    _qual=list(_plt.get_cmap("tab20").colors)+list(_plt.get_cmap("tab20b").colors)
    return {c:_qual[i%len(_qual)] for i,c in enumerate(order)}, order

def heat_text_color(value, vmin, vmax, cmap="plasma"):
    """Readable annotation colour on a plasma(_r) cell."""
    import numpy as np
    t=0.0 if vmax==vmin else (value-vmin)/(vmax-vmin)
    if cmap.endswith("_r"): t=1-t
    return "#222222" if t>0.55 else "white"   # plasma: high end is light (yellow)
