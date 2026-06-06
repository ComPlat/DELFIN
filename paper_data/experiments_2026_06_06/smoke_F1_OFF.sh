#!/bin/bash
# Smoke 500 — ULTIMATE config MINUS F1_COVERAGE
# Hypothesis: F1's relaxed CshM gate (200/400/50 vs 20/60/5) is the root cause
# of structqual_heavy_collapse regression in ULTIMATE voll-pool
set -e
cd /home/qmchem_max/ComPlat/DELFIN
LABEL="SMOKE-F1-OFF-2026-06-06"

export PYTHONHASHSEED=0
export DELFIN_GRIP_LIB_PATH=/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz
export DELFIN_FFFREE_GRIP_LIB_COD_PATH=/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v6_cod.npz

# D1.v2 dispatch + ULTIMATE morgen-flags MINUS F1_COVERAGE
export DELFIN_FFFREE_BUILDER=1
export DELFIN_FFFREE_GRIP=1
export DELFIN_FFFREE_POST_GRIP_ALL=1
export DELFIN_FFFREE_PURE_TRACK3=1
# *** F1_COVERAGE OMITTED ***
export DELFIN_FFFREE_DD_RELAX=1
export DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1
export DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1
export DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY=1
export DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS=1
export DELFIN_FFFREE_BURNSIDE_CONFORMER=1
export DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING=1
export DELFIN_FFFREE_GRIP_HEALING_MODE=1
export DELFIN_FFFREE_HH_CLASH_INCLUDE=1
export DELFIN_FFFREE_SP3_H_UMBRELLA=1
export DELFIN_FFFREE_SP3_H_HEAL=1
export DELFIN_MOGUL_V3_TUNED=1
export DELFIN_GRIP_LOSS_WEIGHTS_TUNED=1
export DELFIN_FFFREE_AROMATIC_BONDS=1
export DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1
export DELFIN_FFFREE_SYMMETRY_EQUIVALENCE=1
export DELFIN_FFFREE_REALISM_SORT=1
export DELFIN_FFFREE_REALISM_TUNED_WEIGHTS=1
export DELFIN_FFFREE_REALISM_SOFT_GATES=1
export DELFIN_FFFREE_DONOR_VSEPR=1
export DELFIN_FFFREE_RING_SCALE=1

# D2 convergence + completeness (KEEP)
export DELFIN_FFFREE_GRIP_MAX_ITER=1000
export DELFIN_FFFREE_GRIP_GTOL=1e-6
export DELFIN_FFFREE_GRIP_RESTARTS=3
export DELFIN_MAX_ISOMERS=300
export DELFIN_FFFREE_GRIP_WEIGHT_BOND=2.5
export DELFIN_FFFREE_GRIP_WEIGHT_ANGLE=1.5

# F2 fallback mode (KEEP, mode=grip)
export DELFIN_FFFREE_FALLBACK_MODE=grip

# Forensik log
export DELFIN_FFFREE_FORENSIK_LOG=/tmp/${LABEL}_dispatch_forensik.tsv

# Fast pipeline
export DELFIN_INPROC=1
export DELFIN_THREAD_WORKERS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "=== SMOKE F1-OFF $LABEL ==="
date
echo "GUPPY HEAD: $(git log -1 --format=%h)"
echo "F1_COVERAGE: ${DELFIN_FFFREE_F1_COVERAGE:-UNSET}"
echo "FALLBACK_MODE: ${DELFIN_FFFREE_FALLBACK_MODE:-UNSET}"

/home/qmchem_max/micromamba/envs/delfin/bin/python \
    /home/qmchem_max/agent_workspace/quality_framework/scripts/pool_evaluator.py \
    /tmp/shared_smiles.txt \
    --shadir /home/qmchem_max/ComPlat/DELFIN \
    --parallel 80 --timeout 300 \
    --no-retry-on-zero \
    --commit-label "$LABEL" \
    --xyz-archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive \
    --progress-every 50
echo "=== Smoke DONE $(date) ==="
ARCH=/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/$LABEL
echo "XYZ count: $(ls $ARCH/*.xyz 2>/dev/null | wc -l)"
cd /home/qmchem_max/agent_workspace/quality_framework
echo "=== Detector battery ==="
PYTHONHASHSEED=0 DELFIN_DETECTOR_JOBS=8 \
  /home/qmchem_max/micromamba/envs/delfin/bin/python scripts/run_all_detectors.py \
  $ARCH --all --skip-xtb --cli-sample 50000 --parallel 16 2>&1 | tail -10
echo ""
echo "=== iter_gate vs ULTIMATE intersect (same SMILES test) ==="
PYTHONHASHSEED=0 /home/qmchem_max/micromamba/envs/delfin/bin/python scripts/iter_gate.py \
  $LABEL \
  --prev forensik_ULTIMATE_intersect 2>&1 | head -60
echo "=== Bisect DONE $(date) ==="
