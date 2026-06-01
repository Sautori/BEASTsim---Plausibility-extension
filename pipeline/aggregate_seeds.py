"""Turn a multi-seed run into replicate statistics for the main comparison.

Reads named_plausibility_baseline.csv from an output folder (the one written by
thesis_pipeline.py, which carries a `seed` column) and reports:
  1. per-variant Balanced Realism across seeds: mean, sd, spread.
  2. the scDesign3 vs reference-free SRTsim comparison on the shared datasets, collapsed
     to one gap per dataset (mean over seeds) plus the across-seed sd, then a sign test on
     the dataset-level gaps.

Seeds vary the evaluation stochasticity (clustering, gene subsampling) on fixed simulated
outputs. Where a metric has no stochastic step for a given pair it is identical across seeds
(across-seed sd = 0), so the seeds are NOT extra independent samples of that quantity. The
test therefore runs on the dataset-level gaps, not on the seed-expanded rows, to avoid
pseudo-replication. It does not re-draw the simulators, which would need their generators.

Usage:
    python3 aggregate_seeds.py <output-folder>     (default: ../outputs)
"""
import sys, os
from math import comb
import pandas as pd, numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "outputs")
nm = pd.read_csv(os.path.join(OUT, "named_plausibility_baseline.csv"))
if "seed" not in nm.columns:
    nm["seed"] = 0
seeds = sorted(nm["seed"].unique())
print(f"[aggregate] {len(seeds)} seed(s): {seeds}")

br = nm[nm.summary == "Balanced Realism"].copy()
cs = nm[nm.summary == "copy_similarity"].copy()
key = ["platform", "simulator", "dataset_variant", "seed"]
m = br.merge(cs[key + ["score"]], on=key, suffixes=("_realism", "_copy"))
m["novelty"] = 1 - m["score_copy"]
m["pair"] = m.platform + "-" + m.dataset_variant

per_variant = m.groupby("simulator")["score_realism"].agg(["mean", "std", "min", "max"]).round(4)
per_variant.to_csv(os.path.join(OUT, "seed_stability_per_variant.csv"))
print("\n=== Balanced Realism across seeds (per variant) ===")
print(per_variant.to_string())

scd = m[m.simulator == "scDesign3-rb"].groupby(["pair", "seed"]).agg(
    realism=("score_realism", "mean"), novelty=("novelty", "mean")).reset_index()
rf = m[m.simulator.isin(["SRTsim-rf-domain", "SRTsim-rf-tissue"])].groupby(["pair", "seed"]).agg(
    realism=("score_realism", "mean"), novelty=("novelty", "mean")).reset_index()
shared = sorted(set(scd["pair"]) & set(rf["pair"]))
comp = scd.merge(rf, on=["pair", "seed"], suffixes=("_scd", "_rf"))
comp = comp[comp.pair.isin(shared)].copy()
comp["realism_gap"] = comp.realism_scd - comp.realism_rf
comp["novelty_gap"] = comp.novelty_scd - comp.novelty_rf
comp.to_csv(os.path.join(OUT, "seed_stability_shared_comparison.csv"), index=False)

# Collapse seeds to one gap per dataset, and records how much the seed moved it.
ds = comp.groupby("pair").agg(
    realism_gap=("realism_gap", "mean"),
    realism_gap_seed_sd=("realism_gap", "std"),
    novelty_gap=("novelty_gap", "mean")).reset_index()
ds["realism_gap_seed_sd"] = ds["realism_gap_seed_sd"].fillna(0.0)

print(f"\n=== scDesign3 vs reference-free SRTsim, shared datasets {shared} ===")
print(ds.round(5).to_string(index=False))
gaps = ds["realism_gap"].values
n = len(gaps); npos = int((gaps > 0).sum())
max_sd = float(ds["realism_gap_seed_sd"].max())
print(f"\nacross-seed sd of the realism gap: max {max_sd:.2e}  (0 means identical across seeds, no evaluation noise)")
print(f"realism gap > 0 on {npos}/{n} shared datasets; mean gap = {gaps.mean():+.4f}, range [{gaps.min():+.4f}, {gaps.max():+.4f}]")
if n >= 1:
    p_sign = sum(comb(n, k) for k in range(npos, n + 1)) / 2 ** n
    print(f"one-sided sign test (gap>0): p = {p_sign:.3f}  (n={n} independent datasets)")
    if n < 6:
        print(f"NOTE: {n} datasets is underpowered for significance. Report this as a consistent, "
              f"reproducible tendency, not a significant result.")
if len(seeds) > 1 and max_sd < 1e-9:
    print("NOTE: the gap is identical across all seeds on these pairs, so the seeds are not extra "
          "independent samples. A test over seed-expanded rows would be pseudo-replication.")
print(f"\nwrote seed_stability_per_variant.csv and seed_stability_shared_comparison.csv to {OUT}")
