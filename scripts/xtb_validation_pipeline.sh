#!/usr/bin/env bash
# G1 — xtb / DFT validation pipeline.
#
# Compares DELFIN method outputs against xtb GFN2 relaxed reference geometries
# on a stratified sample drawn from the 2792332 voll-pool.
#
# Inputs (paths are fixed for this validation campaign):
#   - 2792332-aromatic-symmetry-VOLLPOOL  (UFF baseline)
#   - 065f6f4-fffree-FULLPOOL            (DELFIN fffree)
#   - D1v2-heal-VOLLPOOL                 (DELFIN heal-stack, optional per file)
#
# Outputs:
#   - paper_data/dft_validation_<N>_stratified_<DATE>.csv
#   - paper_data/dft_validation_summary.json
#   - paper_data/dft_validation_<N>_stratified_<DATE>.md
#
# Determinism: filename-sorted strata, every-Nth selection inside each stratum.
# xtb 6.7.1 (GFN2-xTB --opt loose --etemp 300, single-threaded per worker).
set -euo pipefail

REPO_ROOT=${REPO_ROOT:-/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-a4cb8d159d9524d5d}
PY=${PY:-/home/qmchem_max/micromamba/envs/delfin/bin/python}
XTB=${XTB:-/home/qmchem_max/micromamba/bin/xtb}
TOTAL=${TOTAL:-50}
WORKERS=${WORKERS:-16}
MAX_ATOMS=${MAX_ATOMS:-120}
MIN_ATOMS=${MIN_ATOMS:-6}
FRAMES_PER_FILE=${FRAMES_PER_FILE:-4}
PREFER_DIFFERING=${PREFER_DIFFERING:-1}

UFF_DIR=/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/2792332-aromatic-symmetry-VOLLPOOL
FFF_DIR=/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/065f6f4-fffree-FULLPOOL
HEAL_DIR=/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/D1v2-heal-VOLLPOOL

DATE=$(date +%Y_%m_%d)
OUT_CSV=${REPO_ROOT}/paper_data/dft_validation_${TOTAL}_stratified_${DATE}.csv
OUT_JSON=${REPO_ROOT}/paper_data/dft_validation_summary.json
OUT_MD=${REPO_ROOT}/paper_data/dft_validation_${TOTAL}_stratified_${DATE}.md
FALLBACK_FFF_JSON=${REPO_ROOT}/paper_data/dft_validation_fallback_fffree.json
FALLBACK_HEAL_JSON=${REPO_ROOT}/paper_data/dft_validation_fallback_heal.json
WORK_ROOT=/tmp/G1_xtb_validation/work

mkdir -p "${REPO_ROOT}/paper_data" "${WORK_ROOT}"

echo "[G1] xtb version:"
"${XTB}" --version 2>&1 | grep -E "xtb version" | head -1 || true
echo "[G1] config: TOTAL=${TOTAL} WORKERS=${WORKERS} MAX_ATOMS=${MAX_ATOMS}"

echo "[G1] scanning whole pool for UFF-fallback fraction (fffree) ..."
"${PY}" "${REPO_ROOT}/scripts/dft_rmsd_pool_fallback.py" \
    --uff-dir "${UFF_DIR}" --delfin-dir "${FFF_DIR}" \
    --workers "${WORKERS}" --label fffree --out-json "${FALLBACK_FFF_JSON}"

echo "[G1] scanning whole pool for UFF-fallback fraction (heal) ..."
"${PY}" "${REPO_ROOT}/scripts/dft_rmsd_pool_fallback.py" \
    --uff-dir "${UFF_DIR}" --delfin-dir "${HEAL_DIR}" \
    --workers "${WORKERS}" --label heal --out-json "${FALLBACK_HEAL_JSON}"

echo "[G1] running xtb RMSD comparator ..."
"${PY}" "${REPO_ROOT}/scripts/dft_rmsd_comparator.py" \
    --uff-dir   "${UFF_DIR}" \
    --fff-dir   "${FFF_DIR}" \
    --heal-dir  "${HEAL_DIR}" \
    --total     "${TOTAL}" \
    --workers   "${WORKERS}" \
    --max-atoms "${MAX_ATOMS}" \
    --min-atoms "${MIN_ATOMS}" \
    --frames-per-file "${FRAMES_PER_FILE}" \
    --prefer-differing "${PREFER_DIFFERING}" \
    --xtb-bin   "${XTB}" \
    --work-root "${WORK_ROOT}" \
    --out-csv   "${OUT_CSV}" \
    --out-json  "${OUT_JSON}"

echo "[G1] generating markdown report ..."
"${PY}" "${REPO_ROOT}/scripts/dft_rmsd_report.py" \
    --json "${OUT_JSON}" \
    --csv  "${OUT_CSV}" \
    --fallback-fff "${FALLBACK_FFF_JSON}" \
    --fallback-heal "${FALLBACK_HEAL_JSON}" \
    --md   "${OUT_MD}"

echo "[G1] done. CSV=${OUT_CSV}  JSON=${OUT_JSON}  MD=${OUT_MD}"
