"""
Build the run manifest over BEASTsim's simulator-variant set, treating each variant as a
distinct method so reference-based and reference-free can be compared.

Run it from the project root:
    python3 pipeline/build_beastsim_faithful_manifest.py

Deterministic, dataset-specific resolution: every (dataset, method) maps to ONE explicit
expected path under the simulated-data folder. If that path is absent the pair is skipped
and printed under NOT FOUND. It is never substituted with another dataset's file.
claim boundary: similarity = calibration, not validated truth.
"""
import os, re, csv, shutil
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIMROOT = os.path.join(REPO, "benchmarking", "data", "simulated")
DATA = os.path.join(REPO, "ucloud_thesis_run", "data")
MANI = os.path.join(REPO, "ucloud_thesis_run", "manifests", "spatialsimbench_input_manifest.csv")
os.makedirs(DATA, exist_ok=True)

# (platform, dataset_variant, real_h5ad_in_data, sim_subdir, sim_token)
DATASETS = [
    ("MERFISH", "Spapros",   "real_st_spapros_merfish.h5ad",                                  "spapros",   "merfish_spapros"),
    ("MERFISH", "SpatialDE", "real_st_spatialde_merfish.h5ad",                                "spatialde", "merfish_spatialde"),
    ("Visium",  "Spapros",   "Visium_Mouse_Brain_SPAPROS_filtered_celltypes_annotated.h5ad",  "spapros",   "visium_spapros"),
    ("Visium",  "SpatialDE", "Visium_Mouse_Brain_SpatialDE_filtered_celltypes_annotated.h5ad","spatialde", "visium_spatialde"),
    ("Xenium",  "Xenium",    "Xenium-colon_annotated_preprocessed_0.5x_0.5y.h5ad",            "xenium",    ""),
]

def methods_for(platform, subdir, token):
    """[(method_label, expected_relpath_under_SIMROOT)] using BEASTsim's canonical files.
    Deterministic: each label maps to exactly one expected path for THIS dataset."""
    if subdir == "xenium":
        # Xenium files carries no dataset token. reference-free SRTsim / scCube-rfb were not run for Xenium.
        return [
            ("SRTsim-rb-domain", "SRTsim/xenium/SRTsim_rb_domain.h5ad"),
            ("SRTsim-rb-tissue", "SRTsim/xenium/SRTsim_rb_tissue.h5ad"),
            ("scCube-rb",        "scCube/xenium/scCube_rb.h5ad"),          # only rb exists for Xenium (no rfb)
            ("scDesign3-rb",     "scDesign3/xenium/scDesign3_rb_p1_f1_b0.h5ad"),
        ]
    M = [
        ("SRTsim-rb-domain", "SRTsim/%s/SRTsim_%s_rb_domain.h5ad" % (subdir, token)),
        ("SRTsim-rb-tissue", "SRTsim/%s/SRTsim_%s_rb_tissue.h5ad" % (subdir, token)),
        ("SRTsim-rf-domain", "SRTsim/%s/SRTsim_%s_rf_domain_grid_c1.h5ad" % (subdir, token)),
        ("SRTsim-rf-tissue", "SRTsim/%s/SRTsim_%s_rf_tissue_grid_c1.h5ad" % (subdir, token)),
        ("scCube-rfb",       "scCube/%s/scCube_%s_rfb.h5ad" % (subdir, token)),
        ("scDesign3-rb",     "scDesign3/%s/scDesign3_%s_rb_p1_f1_b0.h5ad" % (subdir, token)),
    ]
    if platform == "Visium":
        M.append(("SpatialcoGCN-rb", "spatialcoGCN/spatialcoGCN - rb - %s.h5ad" % subdir))
    return M

def safe(p):
    return "data/" + re.sub(r"[^A-Za-z0-9._-]", "_", os.path.basename(p))

rows, missing = [], []
for platform, dv, real_name, subdir, token in DATASETS:
    real_rel = "data/" + real_name
    if not os.path.exists(os.path.join(DATA, real_name)):
        missing.append("REAL not staged: %s" % real_name)
    for label, rel in methods_for(platform, subdir, token):
        src = os.path.join(SIMROOT, rel)
        if not os.path.exists(src):
            missing.append("%s/%s/%s: expected source not found -> %s" % (platform, dv, label, rel))
            continue
        dst_rel = safe(src)
        dst_abs = os.path.join(REPO, "ucloud_thesis_run", dst_rel)
        if not (os.path.exists(dst_abs) and os.path.getsize(dst_abs) == os.path.getsize(src)):
            shutil.copy2(src, dst_abs)
        rows.append(dict(run_scope="beastsim_faithful", platform=platform, dataset_variant=dv, simulator=label,
            real_h5ad=real_rel, simulated_h5ad=dst_rel,
            coordinate_key_real="obsm/spatial", coordinate_key_sim="obsm/spatial",
            label_key_real="obs/cell_type", label_key_sim="obs/cell_type",
            notes="BEASTsim-faithful variant"))

cols = ["run_scope","platform","dataset_variant","simulator","real_h5ad","simulated_h5ad",
        "coordinate_key_real","coordinate_key_sim","label_key_real","label_key_sim","notes"]
with open(MANI, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

print("manifest written: %d pairs across %d datasets -> %s" % (len(rows), len(DATASETS), MANI))
print("methods:", dict(Counter(r["simulator"] for r in rows)))
print("staged unique sim files:", len({r["simulated_h5ad"] for r in rows}))

bind = {}
for r in rows:
    bind.setdefault(r["simulated_h5ad"], set()).add((r["platform"], r["dataset_variant"]))
clash = {k: v for k, v in bind.items() if len(v) > 1}
if clash:
    print("\n!! INTEGRITY FAIL (one sim file bound to >1 dataset):")
    for k, v in clash.items():
        print("  ", k, "->", sorted(v))
if missing:
    print("\n!! NOT FOUND / SKIPPED (genuinely absent in repo, NOT substituted):")
    for m in missing:
        print("  ", m)
if not clash and not missing:
    print("\nALL variants found, dataset-bound 1:1, and staged. Ready: attach THESIS_BATCH_FULL.sh and run.")
elif not clash:
    print("\nResolved cleanly except the explicit absences above (expected per dataset_inventory_audit). Ready to run.")
