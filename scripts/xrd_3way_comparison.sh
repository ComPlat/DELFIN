#!/bin/bash
# DELFIN vs xtb-relaxed vs UFF — all compared to CCDC XRD ground truth
# Headline claim: DELFIN closer to XRD than xtb-relaxed
set -e
cd /home/qmchem_max/ComPlat/DELFIN

SAMPLE_SIZE="${1:-50}"

# Paths
XRD_REF=/home/qmchem_max/agent_workspace/quality_framework/reference/ccdc_tmc_index_cleaned.jsonl
UFF_POOL=/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/2792332-aromatic-symmetry-VOLLPOOL
DELFIN_POOL=/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/ULTIMATE-VOLLPOOL
XTB_POOL=/tmp/xtb_relaxed_stratified  # produced by G1
REPORT=/home/qmchem_max/agent_workspace/quality_framework/reports/xrd_3way_comparison.json

# Verify dependencies
[ -f "$XRD_REF" ] || { echo "ERROR: $XRD_REF missing"; exit 1; }
[ -d "$UFF_POOL" ] || { echo "ERROR: $UFF_POOL missing"; exit 1; }

if [ ! -d "$DELFIN_POOL" ]; then
  echo "INFO: DELFIN ULTIMATE-VOLLPOOL not yet built. This script needs to wait until it lands."
  exit 0
fi

if [ ! -d "$XTB_POOL" ]; then
  echo "INFO: xtb-relaxed structures not yet produced (G1 still running)."
  exit 0
fi

# Run comparison
PYTHONHASHSEED=0 /home/qmchem_max/micromamba/envs/delfin/bin/python << PYEOF
import json
from pathlib import Path

# Load CCDC XRD ground truth
xrd_index = {}
with open("$XRD_REF") as f:
    for line in f:
        e = json.loads(line)
        xrd_index[e["refcode"]] = e

print(f"Loaded {len(xrd_index)} CCDC XRD references")

# For each refcode found in both DELFIN and xtb outputs:
# - extract DELFIN output coords
# - extract xtb-relaxed coords  
# - extract UFF (2792332) coords
# - extract XRD ground truth coords
# - Kabsch-align all to XRD
# - Compute heavy-atom RMSD per method
# - Aggregate per-class

# Per-class RMSD table:
# Class | UFF | DELFIN | xtb | Winner
# Headline: DELFIN < xtb < UFF
# If true: paper-grade Knockout

print("3-way XRD-RMSD comparison would run here when DELFIN and xtb pools are ready")
PYEOF
