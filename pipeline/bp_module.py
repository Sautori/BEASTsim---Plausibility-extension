"""Biological-plausibility components and the named candidate summaries.
Gate = 1 - max(copy_sim - 0.95, 0) / 0.05. Balanced is the component mean,
Spatial Tissue Fidelity the spatial-component mean, Novelty-Gated is Balanced times the gate,
and coverage Q is the fraction of available components."""
from __future__ import annotations
import numpy as np

SPATIAL_COMPONENTS = ["spatial_gene", "layout", "proximity_composition"]
ALL_COMPONENTS = ["expression", "spatial_gene", "layout", "proximity_composition"]


def novelty_gate(copy_similarity):
    if not np.isfinite(copy_similarity):
        return 1.0
    return float(np.clip(1.0 - max(copy_similarity - 0.95, 0.0) / 0.05, 0.0, 1.0))


def named_summaries(components: dict, copy_similarity: float) -> dict:
    def mean_of(keys):
        vals = [components[k] for k in keys if k in components and np.isfinite(components[k])]
        return float(np.mean(vals)) if vals else float("nan")
    balanced = mean_of(ALL_COMPONENTS)
    spatial = mean_of(SPATIAL_COMPONENTS)
    gate = novelty_gate(copy_similarity)
    novelty = balanced * gate if np.isfinite(balanced) else float("nan")
    avail = sum(1 for k in ALL_COMPONENTS if k in components and np.isfinite(components[k]))
    return {
        "Balanced Realism": balanced,
        "Spatial Tissue Fidelity": spatial,
        "Novelty-Gated Plausibility": novelty,
        "Conservative Softmin Plausibility": float(np.nanmin([components[k] for k in ALL_COMPONENTS if k in components])) if components else float("nan"),
        "coverage_Q": avail / len(ALL_COMPONENTS),
        "novelty_gate": gate,
        "copy_similarity": float(copy_similarity),
    }
