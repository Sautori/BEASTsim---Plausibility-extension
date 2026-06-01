"""Shared helpers for the figure scripts: data loading, the colour palette, and a save wrapper."""
import os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
ROOT=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
import style as B
RUN=os.path.join(ROOT,"..","data")
OUT=os.path.join(ROOT,"..","outputs")
os.makedirs(OUT, exist_ok=True)
B.apply()

# canonical dataset order (platform, dataset_variant) and it's short labels
DATASETS=[("MERFISH","Spapros"),("MERFISH","SpatialDE"),("Visium","Spapros"),("Visium","SpatialDE"),("Xenium","Xenium")]
def dlabel(pl,dv): return "Xenium" if pl=="Xenium" else f"{pl[:3]}·{dv[:4]}"
DLABELS=[dlabel(*d) for d in DATASETS]
METHODS=B.METHOD_ORDER

def load(name): return pd.read_csv(os.path.join(RUN,name))
def metric_long(): return load("metric_long.csv")
def named():       return load("named_plausibility_baseline.csv")
def components():  return load("component_evidence.csv")
def perturb():     return load("perturbation_response.csv")
def expected():    return load("expected_vs_observed.csv")

def coverage():
    """which (dataset,method) cells were actually run, from metric_long pairs."""
    m=metric_long()[["platform","dataset_variant","simulator"]].drop_duplicates()
    run=set((r.platform,r.dataset_variant,r.simulator) for r in m.itertuples())
    grid=pd.DataFrame(index=METHODS, columns=DLABELS, dtype=float)
    for (pl,dv) in DATASETS:
        for me in METHODS:
            grid.loc[me, dlabel(pl,dv)] = 1.0 if (pl,dv,me) in run else np.nan
    return grid

def pivot_metric(metric_name, agg="mean"):
    """method x dataset matrix of one metric's value (NaN where not run)."""
    m=metric_long(); m=m[m.metric==metric_name].copy()
    m["dl"]=[dlabel(p,d) for p,d in zip(m.platform,m.dataset_variant)]
    p=m.pivot_table(index="simulator",columns="dl",values="value",aggfunc=agg)
    return p.reindex(index=METHODS, columns=DLABELS)

def greyscale_ok(path, min_std=6.0):
    """luminance spread sanity: ensures marks aren't indistinguishable in B/W."""
    from matplotlib import image as mpimg
    im=mpimg.imread(path)
    if im.ndim==3:
        lum=0.2126*im[...,0]+0.7152*im[...,1]+0.0722*im[...,2]
        std=float(np.std(lum*255))
        return std, std>=min_std
    return float("nan"), True

def save(fig, name, source_df=None):
    png=os.path.join(OUT,name+".png"); pdf=os.path.join(OUT,name+".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight"); fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    std,ok=greyscale_ok(png)
    print(f"  saved {name}  [greyscale σ={std:.1f} {'ok' if ok else 'LOW-CONTRAST'}]")
    return png


def load_tissue(rel):
    """Robustly load (xy float, cell_type labels) from a staged h5ad. Some sim files store
    obsm/spatial as byte-strings, coerced to float so layout is not mis-rendered as a grid."""
    import h5py, os, numpy as np
    p=os.path.join(ROOT,"..","data",os.path.basename(rel))
    with h5py.File(p,"r") as fh:
        sp=fh["obsm"]["spatial"][()]
        if sp.dtype==object or sp.dtype.kind in ("S","U"):
            xy=np.array([[float(a),float(b)] for a,b in sp])
        else:
            xy=np.asarray(sp,dtype=float)
        g=fh["obs"]["cell_type"]; codes=g["codes"][()]
        cats=[c.decode() if isinstance(c,bytes) else str(c) for c in g["categories"][()]]
        lab=np.array([cats[c] if c>=0 else "NA" for c in codes])
    return xy,lab

def natkey(x):
    try: return (0,float(x))
    except Exception: return (1,str(x))
