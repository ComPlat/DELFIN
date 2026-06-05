#!/bin/bash
# Ultimate voll-pool: F1 (constructive expansion) + F2 (universal embed) + D2 (convergence)
# Target: 0% UFF, 100% coverage, FF-free throughout, deterministic
set -e
cd /home/qmchem_max/ComPlat/DELFIN
LABEL="${1:-ULTIMATE-VOLLPOOL}"

export PYTHONHASHSEED=0
export DELFIN_GRIP_LIB_PATH=/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz
export DELFIN_FFFREE_GRIP_LIB_COD_PATH=/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v6_cod.npz

# D1.v2 dispatch + all morgen-flags
export DELFIN_FFFREE_BUILDER=1
export DELFIN_FFFREE_GRIP=1
export DELFIN_FFFREE_POST_GRIP_ALL=1
export DELFIN_FFFREE_PURE_TRACK3=1
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

# D2 convergence + completeness
export DELFIN_FFFREE_GRIP_MAX_ITER=1000
export DELFIN_FFFREE_GRIP_GTOL=1e-6
export DELFIN_FFFREE_GRIP_RESTARTS=3
export DELFIN_MAX_ISOMERS=300
export DELFIN_FFFREE_GRIP_WEIGHT_BOND=2.5
export DELFIN_FFFREE_GRIP_WEIGHT_ANGLE=1.5

# F2 fallback mode (FF-free embed when constructive fails)
export DELFIN_FFFREE_FALLBACK_MODE=grip
export DELFIN_FFFREE_EMBED_GRIP_MAX_ITER=50

# Forensik log for paper-claim verification
export DELFIN_FFFREE_FORENSIK_LOG=/tmp/${LABEL}_dispatch_forensik.tsv

# Fast pipeline
export DELFIN_INPROC=1
export DELFIN_THREAD_WORKERS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "=== ULTIMATE voll-pool $LABEL ==="
date
echo "GUPPY HEAD: $(git log -1 --format=%h)"
env | grep -E "DELFIN_FFFREE|DELFIN_GRIP|DELFIN_MAX_ISOMERS|DELFIN_INPROC|DELFIN_MOGUL" | sort

/home/qmchem_max/micromamba/envs/delfin/bin/python \
    /home/qmchem_max/agent_workspace/quality_framework/scripts/pool_evaluator.py \
    /home/qmchem_max/agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt \
    --shadir /home/qmchem_max/ComPlat/DELFIN \
    --parallel 160 --timeout 300 \
    --no-retry-on-zero \
    --commit-label "$LABEL" \
    --xyz-archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive \
    --progress-every 200
echo "=== DONE $(date) ==="
