#!/bin/bash
# run_vollpool_D1.sh — Mission D1 voll-pool HEAL launcher (2026-06-05)
#
# Recovers the 2792332 -> 8654d8f voll-pool emission regression by decoupling
# DELFIN_FFFREE_PURE_TRACK3 from DELFIN_FFFREE_NO_FALLBACK auto-imply.
#
# Root cause + heal: see commit 24849e3 (or HEAD if rebased) and
# scripts/MISSION_D1_HEAL_FORENSIK.md.
#
# Pre-heal forensik (8654d8f + 22-flag stack):
#   silent zero-isomer:    9296 / 11452  (81.2 %)
#   zero-byte XYZ:         7196 / 9300   (77.4 %)
#   smoke10 emission:      1 / 10 non-empty
#
# Post-heal target:
#   silent zero-isomer:    < 20 %
#   zero-byte XYZ:         < 5 %
#   smoke10 emission:      >= 6 / 10 non-empty   (verified: 7/10)
#
# voll-pool acceptance:    net+ vs 2792332-aromatic-symmetry-VOLLPOOL
#                          AND < 10 severe regression axes

set -e
cd /home/qmchem_max/ComPlat/DELFIN

LABEL="${1:-D1-heal-VOLLPOOL}"
LOG=/tmp/${LABEL}.log

export PYTHONHASHSEED=0
export DELFIN_GRIP_LIB_PATH=/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v3.npz

# Exact 2792332 22-flag stack (proven production-quality)
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

# Fast-path: INPROC + THREAD_WORKERS bring voll-pool from ~3.5h to ~30 min.
# Per benchmark_runner docstring this is determinism-preserving (in-process
# vs subprocess only differs in fork startup tax, not in output).
export DELFIN_INPROC=1
export DELFIN_THREAD_WORKERS=1
# WITHHOLD the 2ee2f45 BUILDER + DONOR_VSEPR + RING_SCALE extras.
# 2792332 ran without them and won; keeping the proven minimum reduces
# risk of new flag-interaction regressions.

echo "=== D1-HEAL voll-pool $LABEL ===" | tee -a "$LOG"
date | tee -a "$LOG"
env | grep -E "DELFIN_FFFREE|DELFIN_GRIP|DELFIN_MOGUL" | sort | tee -a "$LOG"
echo "GUPPY HEAD: $(git log -1 --format=%h)" | tee -a "$LOG"
echo "Heal commit: $(git log --oneline -1 --grep='Mission D1' || git log -1 --format=%h)" | tee -a "$LOG"

/home/qmchem_max/micromamba/envs/delfin/bin/python \
    /home/qmchem_max/agent_workspace/quality_framework/scripts/pool_evaluator.py \
    /home/qmchem_max/agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt \
    --shadir /home/qmchem_max/ComPlat/DELFIN \
    --parallel 128 --timeout 300 \
    --commit-label "$LABEL" \
    --xyz-archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive \
    --progress-every 100 2>&1 | tee -a "$LOG"

echo "=== D1-HEAL DONE $(date) ===" | tee -a "$LOG"
