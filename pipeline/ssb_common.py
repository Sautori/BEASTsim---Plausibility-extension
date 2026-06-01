"""Shared I/O helpers for the SpatialSimBench metric job (SDU UCloud)."""
from __future__ import annotations
import numpy as np, anndata as ad

def read_h5ad(path):
    return ad.read_h5ad(path)

def dense(a):
    X = a.X
    return np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X, dtype=np.float64)

def common_genes(R, S):
    g = [x for x in R.var_names if x in set(S.var_names)]
    return R[:, g], S[:, g], g

def get_coords(a, key="obsm/spatial"):
    if key.startswith("obsm/"):
        k = key.split("/", 1)[1]
        if k in a.obsm:
            return np.asarray(a.obsm[k])[:, :2].astype(float)
    for k in ("spatial", "X_spatial"):
        if k in a.obsm:
            return np.asarray(a.obsm[k])[:, :2].astype(float)
    for cx, cy in (("X", "Y"), ("x", "y"), ("imagecol", "imagerow")):
        if cx in a.obs and cy in a.obs:
            return a.obs[[cx, cy]].values.astype(float)
    raise ValueError("no spatial coordinates found")

def get_labels(a, key="obs/Cell_Type"):
    cand = []
    if key and key.startswith("obs/"):
        cand.append(key.split("/", 1)[1])
    cand += ["Cell_Type", "cell_type", "Cell_type", "celltype", "domain", "annotation"]
    for c in cand:
        if c in a.obs:
            return a.obs[c].astype(str).values
    return None
