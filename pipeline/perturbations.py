"""Perturbation controls. Each returns (X', xy', labels') for the simulated object.
Severity is one of low, medium, high."""
from __future__ import annotations
import numpy as np
from sklearn.neighbors import NearestNeighbors

SEV = {"low": 0.25, "medium": 0.5, "high": 1.0}

BATTERY = [
    "baseline", "coordinate_shuffle", "label_shuffle", "gene_shuffle",
    "spatial_jitter_low", "spatial_jitter_medium", "spatial_jitter_high",
    "dropout_inflation_low", "dropout_inflation_medium", "dropout_inflation_high",
    "library_size_resampling_low", "library_size_resampling_medium", "library_size_resampling_high",
    "spatial_smoothing_low", "spatial_smoothing_medium", "spatial_smoothing_high",
    "copy_reference", "copy_reference_jitter_low", "copy_reference_jitter_medium", "copy_reference_jitter_high",
    "random_spatial_field_injection_low", "random_spatial_field_injection_medium",
    "invalid_coordinates_control", "missing_label_control",
]
# expected primary damaged layer (for the expected-vs-observed comparision)
EXPECTED = {
    "coordinate_shuffle": "spatial", "spatial_jitter": "spatial", "label_shuffle": "proximity",
    "gene_shuffle": "expression", "dropout_inflation": "expression", "library_size_resampling": "expression",
    "spatial_smoothing": "spatial_gene", "copy_reference": "similarity", "random_spatial_field_injection": "spatial_gene",
    "invalid_coordinates_control": "spatial", "missing_label_control": "proximity",
}


def _sev(name):
    for k, v in SEV.items():
        if name.endswith("_" + k):
            return v, name[: -len("_" + k)]
    return 1.0, name


def apply_perturbation(X, xy, labels, name, rng, real_X=None):
    s, base = _sev(name)
    X = X.astype(np.float64).copy(); xy = xy.astype(np.float64).copy()
    labels = None if labels is None else np.array(labels, dtype=object).copy()
    if base == "baseline":
        return X, xy, labels
    if base == "coordinate_shuffle":
        xy = xy[rng.permutation(xy.shape[0])]
    elif base == "spatial_jitter":
        sc = s * 0.1 * (xy.max(0) - xy.min(0) + 1e-9)
        xy = xy + rng.normal(0, sc, xy.shape)
    elif base == "label_shuffle" and labels is not None:
        labels = labels[rng.permutation(len(labels))]
    elif base == "gene_shuffle":
        X = X[:, rng.permutation(X.shape[1])]
    elif base == "dropout_inflation":
        mask = rng.random(X.shape) < (0.5 * s)
        X[mask] = 0.0
    elif base == "library_size_resampling":
        f = np.exp(rng.normal(0, 0.5 * s, X.shape[0]))
        X = X * f[:, None]
    elif base == "spatial_smoothing":
        k = 6; nn = NearestNeighbors(n_neighbors=k + 1).fit(xy); _, idx = nn.kneighbors(xy)
        sm = X[idx[:, 1:]].mean(1)
        X = (1 - s) * X + s * sm
    elif base == "copy_reference" and real_X is not None:
        m = min(X.shape[0], real_X.shape[0]); X[:m] = real_X[:m]
    elif base == "copy_reference_jitter" and real_X is not None:
        m = min(X.shape[0], real_X.shape[0]); X[:m] = real_X[:m] * np.exp(rng.normal(0, 0.2 * s, (m, X.shape[1])))
    elif base == "random_spatial_field_injection":
        grad = (xy[:, 0] - xy[:, 0].mean()) / (xy[:, 0].std() + 1e-9)
        X = X + s * np.abs(grad)[:, None] * X.mean()
    elif base == "invalid_coordinates_control":
        xy = np.zeros_like(xy)
    elif base == "missing_label_control" and labels is not None:
        labels = np.array(["NA"] * len(labels), dtype=object)
    return X, xy, labels
