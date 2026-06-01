"""
SpatialSimBench metrics (Python port of Liang et al. 2025, SpatialSimBench/03_run_evaluation).

Implements:
  E1 Data properties  - distributional KS test (analog of ks::kde.test()$zstat) on
                        library size, effective library size (TMM), fraction-zero (gene/spot),
                        mean / variance / scaled variance, Pearson/Spearman correlation structure.
  E2 SVG detection    - per-gene Moran's I ranking -> precision / recall of real-vs-sim SVG overlap.
  E2 Moran's I        - global Moran's I (real, sim) + cosine similarity of per-gene profiles.
  E2 Spatial cluster  - PCA + KMeans (or Leiden if scanpy present), then ARI / NMI vs real labels.
  E2 CT proportion    - Jensen-Shannon divergence + RMSE of cell-type proportions.
Lower KS / JSD / RMSE = closer to real. Higher precision/recall/cosine/ARI/NMI = better.

Designed for parallel use: each pair is one worker, and gene-level ops are vectorised numpy.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import ks_2samp, spearmanr
from scipy.spatial.distance import jensenshannon
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def tmm_factors(counts: np.ndarray, ref_col: int | None = None,
                logratio_trim: float = 0.3, sum_trim: float = 0.05) -> np.ndarray:
    """edgeR-style TMM normalisation factors. counts: cells x genes."""
    libsize = counts.sum(1)
    cpm = counts / (libsize[:, None] + 1e-9)
    if ref_col is None:
        f75 = np.nanpercentile(np.where(counts > 0, counts, np.nan), 75, axis=1)
        if np.all(np.isnan(f75)):
            return np.ones(counts.shape[0])
        ref_col = int(np.nanargmin(np.abs(np.nan_to_num(f75 - np.nanmean(f75), nan=np.inf))))
    ref = counts[ref_col]; refsum = ref.sum() + 1e-9
    factors = np.ones(counts.shape[0])
    for i in range(counts.shape[0]):
        obs = counts[i]; nobs = obs.sum() + 1e-9
        mask = (obs > 0) & (ref > 0)
        if mask.sum() < 5:
            continue
        o, r = obs[mask], ref[mask]
        M = np.log2((o / nobs) / (r / refsum))
        A = 0.5 * np.log2((o / nobs) * (r / refsum))
        W = (nobs - o) / (nobs * o) + (refsum - r) / (refsum * r)
        lo, hi = np.quantile(M, [logratio_trim, 1 - logratio_trim])
        alo, ahi = np.quantile(A, [sum_trim, 1 - sum_trim])
        keep = (M >= lo) & (M <= hi) & (A >= alo) & (A <= ahi)
        if keep.sum() < 1:
            continue
        factors[i] = 2 ** (np.sum(W[keep] * M[keep]) / np.sum(W[keep]))
    return factors / np.exp(np.mean(np.log(factors + 1e-9)))


def ks_stat(real: np.ndarray, sim: np.ndarray) -> float:
    real = real[np.isfinite(real)]; sim = sim[np.isfinite(sim)]
    if len(real) < 5 or len(sim) < 5:
        return float("nan")
    return float(ks_2samp(real, sim).statistic)


def _corr_distribution(X: np.ndarray, n_sample: int, method: str, rng) -> np.ndarray:
    n = X.shape[0]; idx = rng.choice(n, min(n_sample, n), replace=False); M = X[idx]
    if method == "pearson":
        C = np.corrcoef(M)
    else:
        C = spearmanr(M, axis=1)[0]
        C = np.atleast_2d(C)
    iu = np.triu_indices(C.shape[0], 1)
    return C[iu]


def morans_i_per_gene(X: np.ndarray, xy: np.ndarray, k: int = 6) -> np.ndarray:
    n = X.shape[0]
    nn = NearestNeighbors(n_neighbors=min(k + 1, n)).fit(xy)
    _, idx = nn.kneighbors(xy); idx = idx[:, 1:]
    Xc = X - X.mean(0); denom = (Xc ** 2).sum(0) + 1e-12
    num = np.zeros(X.shape[1])
    for j in range(idx.shape[1]):
        num += (Xc * Xc[idx[:, j]]).sum(0)
    return (n / (n * idx.shape[1])) * (num / denom)


def data_property_metrics(XR, XS, n_corr_sample=1000, seed=0):
    rng = np.random.default_rng(seed); out = {}
    out["library_size"] = ks_stat(XR.sum(1), XS.sum(1))
    try:
        tr, ts = tmm_factors(XR), tmm_factors(XS)
        out["effective_library_size"] = ks_stat(XR.sum(1) * tr, XS.sum(1) * ts)
    except Exception:
        out["effective_library_size"] = ks_stat(XR.sum(1), XS.sum(1))  # TMM fallback
    out["fraction_zero_gene"] = ks_stat((XR == 0).mean(0), (XS == 0).mean(0))
    out["fraction_zero_spot"] = ks_stat((XR == 0).mean(1), (XS == 0).mean(1))
    out["mean_expression"] = ks_stat(XR.mean(0), XS.mean(0))
    out["variance_expression"] = ks_stat(XR.var(0), XS.var(0))
    libr = XR.sum(1, keepdims=True); libs = XS.sum(1, keepdims=True)
    out["scaled_variance"] = ks_stat((XR / (libr + 1e-9)).var(0), (XS / (libs + 1e-9)).var(0))
    out["pearson_cell_correlation"] = ks_stat(
        _corr_distribution(XR, n_corr_sample, "pearson", rng),
        _corr_distribution(XS, n_corr_sample, "pearson", rng))
    out["spearman_cell_correlation"] = ks_stat(
        _corr_distribution(XR, n_corr_sample, "spearman", rng),
        _corr_distribution(XS, n_corr_sample, "spearman", rng))
    return out


def svg_and_moran_metrics(XR, XS, xyR, xyS, top_frac=0.15, k=6):
    mr = morans_i_per_gene(XR, xyR, k); ms = morans_i_per_gene(XS, xyS, k)
    K = max(5, int(round(top_frac * XR.shape[1])))
    topR = set(np.argsort(-mr)[:K]); topS = set(np.argsort(-ms)[:K])
    ov = len(topR & topS)
    cos = float(np.dot(mr, ms) / (np.linalg.norm(mr) * np.linalg.norm(ms) + 1e-9))
    return {
        "SVG_precision": ov / K, "SVG_recall": ov / K,
        "SVG_jaccard": ov / (len(topR | topS) + 1e-9),
        "moran_global_real": float(np.nanmean(mr)),
        "moran_global_sim": float(np.nanmean(ms)),
        "moran_profile_cosine": cos,
    }


def clustering_metrics(XR, XS, labelsR, n_pc=20, seed=0):
    libr = XR.sum(1, keepdims=True); libs = XS.sum(1, keepdims=True)
    zr = PCA(min(n_pc, XR.shape[1] - 1), random_state=seed).fit_transform(np.log1p(XR / (libr + 1e-9) * 1e4))
    zs = PCA(min(n_pc, XS.shape[1] - 1), random_state=seed).fit_transform(np.log1p(XS / (libs + 1e-9) * 1e4))
    k = len(set(labelsR)) if labelsR is not None else 8
    clR = KMeans(k, n_init=10, random_state=seed).fit_predict(zr)
    clS = KMeans(k, n_init=10, random_state=seed).fit_predict(zs)
    out = {}
    if labelsR is not None:
        out["cluster_ARI_real_vs_label"] = adjusted_rand_score(labelsR, clR)
        out["cluster_NMI_real_vs_label"] = normalized_mutual_info_score(labelsR, clR)
    m = min(len(clR), len(clS))
    out["cluster_ARI_real_vs_sim"] = adjusted_rand_score(clR[:m], clS[:m])
    return out


def ct_proportion_metrics(labelsR, labelsS):
    if labelsR is None or labelsS is None:
        return {}
    cats = sorted(set(labelsR) | set(labelsS))
    pr = np.array([np.mean(labelsR == c) for c in cats])
    ps = np.array([np.mean(labelsS == c) for c in cats])
    return {"ct_proportion_JSD": float(jensenshannon(pr, ps)),
            "ct_proportion_RMSE": float(np.sqrt(np.mean((pr - ps) ** 2)))}


# ---- BEASTsim-side spatial components (feeds the biological-plausibility module) ----
def occupancy_grid_distance(xyR, xyS, grid=8):
    def occ(xy):
        if xy.std() < 1e-9: return np.ones((grid, grid)) / (grid * grid)
        u = (xy - xy.min(0)) / (np.ptp(xy, 0) + 1e-9)
        gx = np.clip((u[:, 0] * grid).astype(int), 0, grid - 1)
        gy = np.clip((u[:, 1] * grid).astype(int), 0, grid - 1)
        H = np.zeros((grid, grid))
        for a, b in zip(gx, gy): H[b, a] += 1
        return H / (H.sum() + 1e-9)
    pr, ps = occ(xyR), occ(xyS)
    return float(0.5 * np.abs(pr - ps).sum())  # total-variation distance


def _ct_cooccurrence(X, xy, labels, k=6):
    cats = sorted(set(labels)); ci = {c: i for i, c in enumerate(cats)}
    n = len(labels); nn = NearestNeighbors(n_neighbors=min(k + 1, n)).fit(xy)
    _, idx = nn.kneighbors(xy); idx = idx[:, 1:]
    M = np.zeros((len(cats), len(cats)))
    lab = np.array([ci[c] for c in labels])
    for i in range(n):
        for j in idx[i]:
            M[lab[i], lab[j]] += 1
    M = M / (M.sum(1, keepdims=True) + 1e-9)  # conditional prob
    return M, cats


def neighborhood_ctd_similarity(xyR, labelsR, xyS, labelsS, k=6):
    if labelsR is None or labelsS is None: return float("nan")
    cats = sorted(set(labelsR) | set(labelsS))
    def cm(xy, lab):
        ci = {c: i for i, c in enumerate(cats)}; n = len(lab)
        nn = NearestNeighbors(n_neighbors=min(k + 1, n)).fit(xy); _, idx = nn.kneighbors(xy); idx = idx[:, 1:]
        M = np.zeros((len(cats), len(cats))); la = np.array([ci[c] for c in lab])
        src = np.repeat(la, idx.shape[1]); tgt = la[idx.ravel()]
        np.add.at(M, (src, tgt), 1.0)
        return (M / (M.sum() + 1e-9)).ravel()
    a, b = cm(xyR, labelsR), cm(xyS, labelsS)
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-9))  # cosine similarity


def copy_risk_similarity(XR, XS):
    mr, ms = XR.mean(0), XS.mean(0)
    return float(np.dot(mr, ms) / (np.linalg.norm(mr) * np.linalg.norm(ms) + 1e-9))


def fast_components(XR, XS, xyR, xyS, labR, labS, gene_cap=800, rng=None):
    if gene_cap and XR.shape[1] > gene_cap:
        rng = rng or np.random.default_rng(0)
        gi = rng.choice(XR.shape[1], gene_cap, replace=False); XR = XR[:, gi]; XS = XS[:, gi]
    """Cheap component scores (used inside the perturbation loop). 0..1, higher = closer."""
    def ks(a, b): return ks_stat(a, b)
    expr = 1 - np.nanmean([ks(XR.sum(1), XS.sum(1)), ks((XR == 0).mean(0), (XS == 0).mean(0)),
                           ks(XR.mean(0), XS.mean(0)), ks(XR.var(0), XS.var(0))])
    sm = svg_and_moran_metrics(XR, XS, xyR, xyS)
    sgene = float(np.nanmean([sm["SVG_precision"], max(sm["moran_profile_cosine"], 0)]))
    layout = 1 - occupancy_grid_distance(xyR, xyS)
    prox = neighborhood_ctd_similarity(xyR, labR, xyS, labS)
    comp = {"expression": float(np.clip(expr, 0, 1)), "spatial_gene": float(np.clip(sgene, 0, 1)),
            "layout": float(np.clip(layout, 0, 1))}
    if np.isfinite(prox): comp["proximity_composition"] = float(np.clip(prox, 0, 1))
    return comp, copy_risk_similarity(XR, XS)


# ---- BEASTsim similarity benchmark: soft Dice / DSC (Fig. 2) ----
def soft_dice(p, q, pi=2):
    """BEASTsim _dice_score: 1 - sum|p-q|^pi / (sum p^pi + sum q^pi). 1 = identical."""
    p = np.asarray(p, float).ravel(); q = np.asarray(q, float).ravel()
    num = np.sum(np.abs(p - q) ** pi); den = np.sum(p ** pi) + np.sum(q ** pi) + 1e-12
    return float(1 - num / den)


def beastsim_similarity(XR, XS, xyR, xyS, labR, labS, grid=10):
    """CT / SVG / Shape mean Dice scores (DSC), BEASTsim similarity benchmark."""
    out = {}
    # Shape DSC: occupancy grid distribtions
    def occ(xy):
        if xy.std() < 1e-9: return np.ones(grid * grid) / (grid * grid)
        u = (xy - xy.min(0)) / (np.ptp(xy, 0) + 1e-9)
        gx = np.clip((u[:, 0] * grid).astype(int), 0, grid - 1)
        gy = np.clip((u[:, 1] * grid).astype(int), 0, grid - 1)
        H = np.zeros(grid * grid)
        for a, b in zip(gx, gy): H[b * grid + a] += 1
        return H / (H.sum() + 1e-9)
    out["Shape_mean_dice_score"] = soft_dice(occ(xyR), occ(xyS))
    # SVG DSC: per-gene Moran (spatial-gene signal) distributions
    mr = morans_i_per_gene(XR, xyR); ms = morans_i_per_gene(XS, xyS)
    out["SVG_mean_dice_score"] = soft_dice(np.clip(mr, 0, None), np.clip(ms, 0, None))
    # CT DSC: per-cell-type occupancy distribution similarity, averaged accross types
    if labR is not None and labS is not None:
        cats = sorted(set(labR) | set(labS)); ds = []
        for c in cats:
            mr_ = xyR[np.array(labR) == c]; ms_ = xyS[np.array(labS) == c]
            if len(mr_) >= 3 and len(ms_) >= 3:
                ds.append(soft_dice(occ(mr_), occ(ms_)))
        out["CT_mean_dice_score"] = float(np.mean(ds)) if ds else float("nan")
    return out


# ---- BEASTsim count-based biological signal: differential-expression preservation ----
def de_signal_preservation(XR, XS, labelsR, labelsS, top_k=50):
    """Per cell-type one-vs-rest log fold-change, then compare real vs sim DE structure.
    Returns LFC correlation (higher=better) and top-DE-gene overlap (higher=better)."""
    if labelsR is None or labelsS is None:
        return {}
    labelsR = np.asarray(labelsR); labelsS = np.asarray(labelsS)
    cats = [c for c in sorted(set(labelsR) & set(labelsS))]
    def lfc(X, lab, c):
        ing = X[lab == c]; out = X[lab != c]
        if len(ing) < 3 or len(out) < 3: return None
        return np.log2((ing.mean(0) + 1) / (out.mean(0) + 1))
    cors, ovs = [], []
    for c in cats:
        lr, ls = lfc(XR, labelsR, c), lfc(XS, labelsS, c)
        if lr is None or ls is None or np.std(lr) < 1e-9 or np.std(ls) < 1e-9:
            continue
        cors.append(float(np.corrcoef(lr, ls)[0, 1]))
        k = min(top_k, len(lr) // 3 + 1)
        tr, ts = set(np.argsort(-lr)[:k]), set(np.argsort(-ls)[:k])
        ovs.append(len(tr & ts) / k)
    out = {}
    if cors: out["DE_lfc_correlation"] = float(np.nanmean(cors))
    if ovs: out["DE_top_marker_overlap"] = float(np.mean(ovs))
    return out
