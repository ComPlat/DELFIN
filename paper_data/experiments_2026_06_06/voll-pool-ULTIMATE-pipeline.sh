#!/bin/bash
PARENT_PID=$(pgrep -f "ULTIMATE-2026-06-06-VOLLPOOL" | head -1)
echo "Waiting for ULTIMATE voll-pool PID $PARENT_PID..."
while kill -0 "$PARENT_PID" 2>/dev/null; do sleep 60; done
echo "ULTIMATE voll-pool finished at $(date)"

ARCH=/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/ULTIMATE-2026-06-06-VOLLPOOL
echo "XYZ count: $(ls $ARCH/*.xyz 2>/dev/null | wc -l)"

cd /home/qmchem_max/agent_workspace/quality_framework
echo "=== Detector battery ==="
PYTHONHASHSEED=0 DELFIN_DETECTOR_JOBS=8 \
  /home/qmchem_max/micromamba/envs/delfin/bin/python scripts/run_all_detectors.py \
  $ARCH --all --skip-xtb --cli-sample 50000 --parallel 16 2>&1 | tail -10

echo ""
echo "=== iter_gate ULTIMATE vs 2792332 (DECISIVE) ==="
PYTHONHASHSEED=0 /home/qmchem_max/micromamba/envs/delfin/bin/python scripts/iter_gate.py \
  ULTIMATE-2026-06-06-VOLLPOOL \
  --prev 2792332-aromatic-symmetry-VOLLPOOL 2>&1

echo ""
echo "=== Forensik log analysis: native rate ==="
TSV=/tmp/ULTIMATE-2026-06-06-VOLLPOOL_dispatch_forensik.tsv
if [ -f "$TSV" ]; then
  echo "Total entries: $(wc -l < $TSV)"
  echo "Native rate breakdown:"
  awk -F'\t' '{print $NF}' $TSV | sort | uniq -c | sort -rn | head -10
fi

echo ""
echo "=== Pipeline DONE $(date) ==="
