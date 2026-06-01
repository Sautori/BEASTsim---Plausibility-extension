"""Regenerate every data-facing figure into ../outputs.

Run from this folder: ``python make_figures.py``. The CSV-driven figures work from the
summary tables in ../data. The tissue figures also need the per-pair .h5ad files from the
Zenodo deposit (see the README). They are skipped if those files are absent.
"""
import os, sys, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))

CSV_FIGURES = [
    "panels_core.py", "panels_extra.py", "panels_summary.py", "spatial_autocorrelation.py",
    "scalability.py", "perturbation_response.py", "novelty_gate_effect.py",
    "realism_novelty_pareto.py", "gate_calibration.py", "recommendation_sweep.py",
    "evidence_coverage.py", "platform_effect.py", "fair_clustering.py",
    "reference_mode_similarity.py", "shared_dataset_domination.py", "fidelity_axis_correlation.py",
    "appendix_dose_response.py", "appendix_perturbation_ledger.py",
]
TISSUE_FIGURES = [
    "tissue_gallery.py", "real_vs_simulated.py", "neighbourhood_cooccurrence.py",
    "celltype_spatial_overlap.py", "appendix_per_celltype.py", "dataset_scope.py",
]

def run(name):
    res = subprocess.run([sys.executable, name], cwd=HERE)
    if res.returncode != 0:
        print("  could not run %s (missing inputs?)" % name)

if __name__ == "__main__":
    for name in CSV_FIGURES:
        run(name)
    for name in TISSUE_FIGURES:
        run(name)
    print("figures written to ../outputs")
