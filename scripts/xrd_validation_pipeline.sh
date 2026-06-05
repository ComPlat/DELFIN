#!/bin/bash
# Mission G1' — XRD validation pipeline orchestrator.
#
# Phases:
#   1) Match-table:   identify SMILES with CCDC-resolved counterpart
#   2) RMSD compute:  per-builder pool vs CCDC ground-truth
#   3) Aggregate:     per-class table + summary
#   4) DELFIN build:  build NEW DELFIN F2-grip outputs for matched SMILES
#                     (only if --build is passed; off by default - slow)
#
# Usage:
#   bash scripts/xrd_validation_pipeline.sh                 # phases 1-3 only
#   bash scripts/xrd_validation_pipeline.sh --build         # adds phase 4
set -euo pipefail

cd /home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-aaf5fe6f6b5383c5e

PY=/home/qmchem_max/micromamba/envs/delfin/bin/python
ARCHIVE_ROOT=/home/qmchem_max/agent_workspace/quality_framework/xyz_archive
CCDC_XYZ=/home/qmchem_max/agent_workspace/quality_framework/reference/ccdc_cleaned_xyz
CCDC_INDEX=/home/qmchem_max/agent_workspace/quality_framework/reference/ccdc_tmc_index_cleaned.jsonl
OUTDIR=paper_data

BUILD_DELFIN=0
for arg in "$@"; do
  case "$arg" in
    --build) BUILD_DELFIN=1 ;;
  esac
done

export PYTHONHASHSEED=0

mkdir -p "$OUTDIR"

echo "================================================================"
echo "Mission G1' XRD validation pipeline"
echo "GUPPY HEAD: $(git log -1 --format='%h %s')"
echo "Archive root: $ARCHIVE_ROOT"
echo "CCDC ground-truth dir: $CCDC_XYZ"
echo "================================================================"

# ---------------------------------------------------------------------
# Phase 2 — Per-builder RMSD computation
# (Phase 1 is the match-table emitted by xrd_rmsd_aggregate.py from
#  the union of per-builder JSONL outputs.)
# ---------------------------------------------------------------------
echo
echo "## Phase 2: per-builder RMSD vs CCDC XRD"

for spec in \
    "iter27-uff-track1:uff_baseline" \
    "2792332-aromatic-symmetry-VOLLPOOL:delfin_aromatic" \
    "065f6f4-fffree-FULLPOOL:delfin_fffree"; do
  pool_dir="${spec%%:*}"
  label="${spec##*:}"
  out_jsonl="$OUTDIR/xrd_rmsd_${label}.jsonl"
  echo "  [$label] pool=$pool_dir"
  $PY scripts/xrd_rmsd_comparator.py \
      --pool "$ARCHIVE_ROOT/$pool_dir" \
      --label "$label" \
      --ccdc-xyz "$CCDC_XYZ" \
      --ccdc-index "$CCDC_INDEX" \
      --out "$out_jsonl" 2>&1 | tail -4
done

# Phase 4 — build NEW DELFIN F2/grip output for matched-refcode SMILES
if [ "$BUILD_DELFIN" = "1" ]; then
  echo
  echo "## Phase 4: building DELFIN F2-grip output for matched-refcode SMILES"
  if [ ! -f "$OUTDIR/xrd_subset_smiles.txt" ]; then
    echo "ERR: $OUTDIR/xrd_subset_smiles.txt not found"
    exit 2
  fi
  SUBSET=$OUTDIR/xrd_subset_smiles.txt
  LABEL="G1prime-XRD-subset"
  export DELFIN_FFFREE_BUILDER=1
  export DELFIN_FFFREE_GRIP=1
  export DELFIN_FFFREE_POST_GRIP_ALL=1
  export DELFIN_FFFREE_PURE_TRACK3=1
  export DELFIN_FFFREE_DD_RELAX=1
  export DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1
  export DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1
  export DELFIN_FFFREE_AROMATIC_BONDS=1
  export DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1
  export DELFIN_FFFREE_DONOR_VSEPR=1
  export DELFIN_FFFREE_RING_SCALE=1
  export DELFIN_FFFREE_GRIP_MAX_ITER=1000
  export DELFIN_FFFREE_GRIP_GTOL=1e-6
  export DELFIN_FFFREE_GRIP_RESTARTS=3
  export DELFIN_MAX_ISOMERS=300
  export DELFIN_FFFREE_FALLBACK_MODE=grip
  export DELFIN_INPROC=1
  export DELFIN_THREAD_WORKERS=1
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1

  echo "  smiles subset: $(wc -l < $SUBSET) entries"
  $PY /home/qmchem_max/agent_workspace/quality_framework/scripts/pool_evaluator.py \
      "$SUBSET" \
      --shadir /home/qmchem_max/ComPlat/DELFIN \
      --parallel 80 --timeout 600 \
      --no-retry-on-zero \
      --commit-label "$LABEL" \
      --xyz-archive "$ARCHIVE_ROOT" \
      --progress-every 5 2>&1 | tail -20 || true
  echo
  echo "  scoring new DELFIN build vs CCDC..."
  $PY scripts/xrd_rmsd_comparator.py \
      --pool "$ARCHIVE_ROOT/$LABEL" \
      --label "delfin_G1prime" \
      --ccdc-xyz "$CCDC_XYZ" \
      --ccdc-index "$CCDC_INDEX" \
      --out "$OUTDIR/xrd_rmsd_delfin_G1prime.jsonl" 2>&1 | tail -4
fi

# Phase 3 — Aggregate
echo
echo "## Phase 3: aggregate per-class + summary"
AGG_INPUTS=(
    --in "$OUTDIR/xrd_rmsd_uff_baseline.jsonl"
    --in "$OUTDIR/xrd_rmsd_delfin_aromatic.jsonl"
    --in "$OUTDIR/xrd_rmsd_delfin_fffree.jsonl"
)
if [ "$BUILD_DELFIN" = "1" ] && [ -s "$OUTDIR/xrd_rmsd_delfin_G1prime.jsonl" ]; then
  AGG_INPUTS+=(--in "$OUTDIR/xrd_rmsd_delfin_G1prime.jsonl")
fi
$PY scripts/xrd_rmsd_aggregate.py \
    "${AGG_INPUTS[@]}" \
    --out-dir "$OUTDIR"

echo
echo "## Done."
echo "Outputs:"
ls -la "$OUTDIR"/xrd_*.{csv,md,json,jsonl} 2>/dev/null
