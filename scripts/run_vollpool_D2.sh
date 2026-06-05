#!/bin/bash
# run_vollpool_D2.sh — Mission D2 voll-pool launcher (2026-06-05)
#
# Mission D2: Beat UFF on UFF's own turf via tighter convergence + complete
# isomer/conformer enumeration.  Strategy:
#
#   * 5-7x deeper L-BFGS convergence (200 -> 1000 iter, 1e-4 -> 1e-6 gtol)
#   * 3-restart L-BFGS for robustness against shallow local minima
#   * GRIP loss weights bumped on bond (1.0 -> 2.5) + angle (0.5 -> 1.5) so
#     the Mahalanobis pull on internal geometry outranks the rest of the
#     loss budget.
#   * GRACE max_per_isomer 10 -> 100  (10x deeper conformer ensemble per isomer)
#   * smiles_converter max_isomers passed-through 50 -> 300  (gate-filtered)
#
# All new flags are default-OFF byte-identical with HEAD 8654d8f -- without
# them the build path is exactly the same as the C1 heal launcher.
#
# Determinism: PYTHONHASHSEED=0 + every restart uses a seeded perturbation,
# so 2-run output is byte-identical (verified by tests/test_d2_convergence.py).

set -e
# Mission D2 lives on branch feat-D2-completeness-convergence in this
# worktree; cd here so `sys.path[0]=''` resolves to the D2 delfin.
D2_WORKTREE=/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-a84fe85cb24392ed8
cd "$D2_WORKTREE"

LABEL="${1:-D2-VOLLPOOL}"
LOG=/tmp/${LABEL}.log

export PYTHONHASHSEED=0
# Force the D2 worktree onto PYTHONPATH so even if pool_evaluator runs
# from a different cwd, it imports the D2 delfin code.
export PYTHONPATH="$D2_WORKTREE:${PYTHONPATH:-}"
export DELFIN_GRIP_LIB_PATH=/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v3.npz

# ==== D1 heal-restore (2792332 race-stack) ====
export DELFIN_FFFREE_GRIP=1
export DELFIN_FFFREE_POST_GRIP_ALL=1
export DELFIN_FFFREE_PURE_TRACK3=1
export DELFIN_FFFREE_DD_RELAX=1
export DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1
export DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1
export DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY=1
export DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS=1
export DELFIN_FFFREE_BURNSIDE_CONFORMER=1

# ==== D1 verified-passing module flags (RESTORED) ====
export DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING=1
export DELFIN_FFFREE_GRIP_HEALING_MODE=1
export DELFIN_FFFREE_HH_CLASH_INCLUDE=1
export DELFIN_FFFREE_SP3_H_UMBRELLA=1
export DELFIN_FFFREE_SP3_H_HEAL=1
export DELFIN_MOGUL_V3_TUNED=1
export DELFIN_GRIP_LOSS_WEIGHTS_TUNED=1

# ==== Aromatic-fix + symmetry ====
export DELFIN_FFFREE_AROMATIC_BONDS=1
export DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1
export DELFIN_FFFREE_SYMMETRY_EQUIVALENCE=1

# ==== Realism sort ====
export DELFIN_FFFREE_REALISM_SORT=1
export DELFIN_FFFREE_REALISM_TUNED_WEIGHTS=1
export DELFIN_FFFREE_REALISM_SOFT_GATES=1

# ==== Construction defaults ====
export DELFIN_FFFREE_BUILDER=1
export DELFIN_FFFREE_DONOR_VSEPR=1
export DELFIN_FFFREE_RING_SCALE=1

# ==== Mission D2 new flags (TIGHT-CONVERGENCE STACK) ====
# Fix 1: deeper L-BFGS convergence (200 -> 1000 iter)
export DELFIN_FFFREE_GRIP_MAX_ITER=1000
# Fix 2: tighter gradient tolerance (1e-4 -> 1e-6, UFF-level precision)
export DELFIN_FFFREE_GRIP_GTOL=1e-6
# Fix 3: 3-restart L-BFGS for robustness
export DELFIN_FFFREE_GRIP_RESTARTS=3
export DELFIN_FFFREE_GRIP_RESTART_PERTURB=0.05
# Fix 6: per-class weight bumps (bond 1.0 -> 2.5; angle 0.5 -> 1.5)
export DELFIN_FFFREE_GRIP_WEIGHT_BOND=2.5
export DELFIN_FFFREE_GRIP_WEIGHT_ANGLE=1.5
# Fix 5: GRACE 10x deeper conformer ensemble per isomer
export DELFIN_FFFREE_GRACE_ENABLE=1
export DELFIN_FFFREE_GRACE_MAX_PER_ISOMER=100
export DELFIN_FFFREE_GRACE_MAX_ISOMERS=300

# ==== Fast pipeline ====
export DELFIN_INPROC=1
export DELFIN_THREAD_WORKERS=1

echo "=== D2 voll-pool $LABEL ===" | tee -a "$LOG"
date | tee -a "$LOG"
env | grep -E "DELFIN_FFFREE|DELFIN_GRIP|DELFIN_MOGUL|DELFIN_INPROC" | sort | tee -a "$LOG"
echo "GUPPY HEAD: $(git log -1 --format=%h)" | tee -a "$LOG"

/home/qmchem_max/micromamba/envs/delfin/bin/python \
    /home/qmchem_max/agent_workspace/quality_framework/scripts/pool_evaluator.py \
    /home/qmchem_max/agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt \
    --shadir /home/qmchem_max/ComPlat/DELFIN \
    --parallel 160 --timeout 300 \
    --no-retry-on-zero \
    --commit-label "$LABEL" \
    --xyz-archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive \
    --progress-every 100 \
    --max-isomers 300 \
    2>&1 | tee -a "$LOG"

echo "=== DONE $(date) ===" | tee -a "$LOG"
