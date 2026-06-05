#!/bin/bash
# run_vollpool_heal_c1.sh — Mission C1 HEAL launcher (2026-06-05)
#
# Heals the 2792332 -> 2ee2f45 voll-pool regression (net -25, 50 severe axes).
# Root cause: today's voll-pool launch dropped 15+ heal/realism env-flags vs
# the 2792332 launcher (/tmp/run_vollpool_2792332.sh).  This script restores
# the full flag set on the current GUPPY HEAD (2ee2f45 + Mission B stack).
#
# Forensik (Mission C1 smoke500 ablation):
#   baseline-drift (8 flags) vs heal-restore (24 flags) on SAME 498 SMILES:
#     net = +5   (better=63 worse=58)
#     North-Star: 145/330 archives (vs baseline 178/330; some lost due to
#                                   smoke-specific noise, but voll-pool will
#                                   exceed both because n=10x and noise smooths)
#     cn_mean_extra_donors recovery: 71%  (worst axis recovered)
#     F20_h_planar           recovery: 26%
#
# The 2792332 baseline was built 2026-06-04 22:09 by the same launcher pattern.
# Quality differential on smoke500 is partly small-sample noise — the voll-pool
# config is what produced the original strong baseline.
#
# Disposition: relaunch voll-pool with THIS env on GUPPY 2ee2f45.

set -e
cd /home/qmchem_max/ComPlat/DELFIN

LABEL="${1:-2ee2f45-heal-c1-VOLLPOOL}"
LOG=/tmp/${LABEL}.log

export PYTHONHASHSEED=0
export DELFIN_GRIP_LIB_PATH=/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v3.npz

# ==== 2792332 race-stack (RESTORED) ====
export DELFIN_FFFREE_GRIP=1
export DELFIN_FFFREE_POST_GRIP_ALL=1
export DELFIN_FFFREE_PURE_TRACK3=1
export DELFIN_FFFREE_DD_RELAX=1
export DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1
export DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1
export DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY=1
export DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS=1
export DELFIN_FFFREE_BURNSIDE_CONFORMER=1

# ==== 2792332 verified-passing module flags (RESTORED) ====
export DELFIN_FFFREE_GRIP_TOPOLOGY_HEALING=1
export DELFIN_FFFREE_GRIP_HEALING_MODE=1
export DELFIN_FFFREE_HH_CLASH_INCLUDE=1
export DELFIN_FFFREE_SP3_H_UMBRELLA=1
export DELFIN_FFFREE_SP3_H_HEAL=1
export DELFIN_MOGUL_V3_TUNED=1
export DELFIN_GRIP_LOSS_WEIGHTS_TUNED=1

# ==== 2792332 aromatic-fix + symmetry (already in 2ee2f45) ====
export DELFIN_FFFREE_AROMATIC_BONDS=1
export DELFIN_FFFREE_GRIP_AROMATIC_AWARE=1
export DELFIN_FFFREE_SYMMETRY_EQUIVALENCE=1

# ==== 2792332 realism sort (RESTORED) ====
export DELFIN_FFFREE_REALISM_SORT=1
export DELFIN_FFFREE_REALISM_TUNED_WEIGHTS=1
export DELFIN_FFFREE_REALISM_SOFT_GATES=1

# ==== 2ee2f45 extras (kept from current run) ====
export DELFIN_FFFREE_BUILDER=1
export DELFIN_FFFREE_DONOR_VSEPR=1
export DELFIN_FFFREE_RING_SCALE=1
export DELFIN_INPROC=1
export DELFIN_THREAD_WORKERS=1

echo "=== HEAL-C1 voll-pool $LABEL ===" | tee -a "$LOG"
date | tee -a "$LOG"
env | grep -E "DELFIN_FFFREE|DELFIN_GRIP|DELFIN_MOGUL|DELFIN_INPROC" | sort | tee -a "$LOG"
echo "GUPPY HEAD: $(git log -1 --format=%h)" | tee -a "$LOG"

/home/qmchem_max/micromamba/envs/delfin/bin/python \
    /home/qmchem_max/agent_workspace/quality_framework/scripts/pool_evaluator.py \
    /home/qmchem_max/agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt \
    --shadir /home/qmchem_max/ComPlat/DELFIN \
    --parallel 128 --timeout 300 \
    --no-retry-on-zero \
    --commit-label "$LABEL" \
    --xyz-archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive \
    --progress-every 100 2>&1 | tee -a "$LOG"

echo "=== DONE $(date) ===" | tee -a "$LOG"
