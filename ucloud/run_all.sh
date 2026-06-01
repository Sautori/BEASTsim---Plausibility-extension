#!/usr/bin/env bash
# Full thesis evaluation run on UCloud.
# UCloud starts this script from /work, so it finds the job folder rather than assuming its
# own path. Machine used for the final run: cpu-amd-zen5, 128 cores, 384 GB.
set -uo pipefail

SEARCH_BASE="${SEARCH_BASE:-/work}"

# Bootstrap: if the scripts are still packed in an uploaded zip, unwrap it first.
if [ -z "$(find "${SEARCH_BASE}" -maxdepth 6 -name thesis_pipeline.py -print -quit 2>/dev/null)" ]; then
  BUNDLE="$(find "${SEARCH_BASE}" -maxdepth 4 -iname '*.zip' -print -quit 2>/dev/null || true)"
  if [ -n "${BUNDLE}" ]; then
    echo "[thesis] unwrapping ${BUNDLE} ..."
    unzip -oq "${BUNDLE}" -d "$(dirname "${BUNDLE}")"
  fi
fi

echo "[thesis] locating job folder under ${SEARCH_BASE} ..."
PIPELINE="$(find "${SEARCH_BASE}" -maxdepth 6 -name thesis_pipeline.py -print -quit 2>/dev/null || true)"
if [ -z "${PIPELINE}" ]; then
  echo "ERROR: scripts/thesis_pipeline.py not found under ${SEARCH_BASE}."
  echo "Upload the UPLOAD_TO_UCLOUD contents (scripts/ data/ manifests/ requirements.txt run_all.sh)."
  ls -la "${SEARCH_BASE}" 2>/dev/null
  exit 1
fi
JOB_ROOT="$(cd "$(dirname "${PIPELINE}")/.." && pwd)"
OUT="${JOB_ROOT}/outputs"; mkdir -p "${OUT}"
echo "[thesis] JOB_ROOT=${JOB_ROOT}"

# Diagnostics: show what actually got uploaded.
echo "[thesis] contents of JOB_ROOT:"; ls -la "${JOB_ROOT}"
NDATA="$(ls "${JOB_ROOT}/data"/*.h5ad 2>/dev/null | wc -l)"
NMAN="$(tail -n +2 "${JOB_ROOT}/manifests/spatialsimbench_input_manifest.csv" 2>/dev/null | wc -l)"
echo "[thesis] data files present: ${NDATA} | manifest rows: ${NMAN}"
if [ "${NDATA}" -eq 0 ]; then
  echo "ERROR: no .h5ad files under ${JOB_ROOT}/data. The datasets did not upload."
  exit 1
fi

export WORKERS="${WORKERS:-$(nproc)}"
# SEEDS: five seeds by default for the replication (sign test on the shared datasets).
# Set SEEDS="0" for the original single-pass run.
export SEEDS="${SEEDS:-0,1,2,3,4}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-2}"

echo "[thesis] installing requirements (workers=${WORKERS}, seeds=${SEEDS}) ..."
python3 -m pip install --user -r "${JOB_ROOT}/requirements.txt" \
  || python3 -m pip install --user numpy pandas scipy scikit-learn anndata h5py psutil \
  || echo "[thesis] WARNING: pip install reported an error; trying to run with whatever is present."

echo "[thesis] starting pipeline ..."
cd "${JOB_ROOT}/scripts"
python3 thesis_pipeline.py \
    --manifest    "${JOB_ROOT}/manifests/spatialsimbench_input_manifest.csv" \
    --data-root   "${JOB_ROOT}" \
    --output-root "${OUT}" \
    --workers     "${WORKERS}" \
    --corr-sample 5000 --max-obs 0 --pert-max-obs 5000 --replicates 5 --gene-cap 8000 --seeds "${SEEDS}" \
    2>&1 | tee "${OUT}/run_stdout.log"

# Replicate statistics for the main comparison (sign test on the shared datasets).
python3 aggregate_seeds.py "${OUT}" 2>&1 | tee "${OUT}/seed_stability.log" || true

echo "[thesis] DONE. Outputs in ${OUT}:"; ls -la "${OUT}"
