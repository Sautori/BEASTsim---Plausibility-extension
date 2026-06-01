#!/usr/bin/env bash
# Build a single self-contained UCloud bundle: scripts, manifest, requirements, run_all.sh,
# and the datasets, all in one zip. On UCloud you upload just this zip and attach run_all.sh
# as the batch script. It unwraps the bundle itself before running.
# Usage:
#   ucloud/make_upload.sh [DATA_DIR]
# DATA_DIR is an optional folder of .h5ad files (the Zenodo deposit). If given, the datasets
# are packed into the bundle. If not, a note is left and you add them before uploading.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
STAGE="${REPO}/ucloud_bundle"
ZIP="${REPO}/ucloud_bundle.zip"
DATA_DIR="${1:-}"

rm -rf "${STAGE}" "${ZIP}"
mkdir -p "${STAGE}/scripts" "${STAGE}/manifests" "${STAGE}/data"

cp "${REPO}/pipeline/"*.py                            "${STAGE}/scripts/"
cp "${REPO}/data/spatialsimbench_input_manifest.csv"  "${STAGE}/manifests/"
cp "${REPO}/requirements.txt"                         "${STAGE}/"
cp "${REPO}/ucloud/run_all.sh"                        "${STAGE}/"

if [ -n "${DATA_DIR}" ] && ls "${DATA_DIR}"/*.h5ad >/dev/null 2>&1; then
  cp "${DATA_DIR}"/*.h5ad "${STAGE}/data/"
  echo "packed $(ls "${STAGE}/data"/*.h5ad | wc -l) datasets into the bundle"
else
  printf 'Place the .h5ad datasets from the Zenodo deposit in this folder before zipping:\nhttps://zenodo.org/records/15775564\n' > "${STAGE}/data/PUT_DATASETS_HERE.txt"
  echo "no DATA_DIR given: add the Zenodo .h5ad files to ${STAGE}/data/ before zipping"
fi

# -0 (store) keeps it fast since the h5ad files are already compressed. Files sit at the top of the
# zip so it unravels cleanly.
( cd "${STAGE}" && zip -0 -rq "${ZIP}" scripts manifests data requirements.txt run_all.sh )
rm -rf "${STAGE}"
echo "built ${ZIP}"
echo "On UCloud: upload ucloud_bundle.zip, start terminal-ubuntu, attach run_all.sh as the batch script."
echo "run_all.sh unwraps the bundle, then runs. Set SEEDS=0,1,2,3,4 to replicate the main comparison."
