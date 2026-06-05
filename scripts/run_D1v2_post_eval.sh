#!/bin/bash
# Post-D1v2-voll-pool evaluation: detector battery + iter_gate vs 2792332
# To be run after `run_vollpool_D1.sh D1v2-heal-VOLLPOOL` completes.

set -e
LABEL="${1:-D1v2-heal-VOLLPOOL}"
PREV="${2:-2792332-aromatic-symmetry-VOLLPOOL}"
ARCHIVE_DIR="/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/$LABEL"
REPORTS_DIR="/home/qmchem_max/agent_workspace/quality_framework/reports"

echo "=== D1.v2 post-pool evaluation $(date) ==="
echo "Archive: $ARCHIVE_DIR"
echo "Comparing vs: $PREV"

# Sanity: archive must exist and be substantial
N=$(ls $ARCHIVE_DIR/ 2>/dev/null | wc -l)
NE=$(find $ARCHIVE_DIR/ -size +0 -name "*.xyz" 2>/dev/null | wc -l)
echo "Pool size: $N total, $NE non-empty"
if [ $N -lt 10000 ]; then
    echo "WARNING: pool is smaller than 10000 -- may be incomplete"
fi

# Detector battery (parallel CLI detectors)
export DELFIN_DETECTOR_JOBS=8
echo ""
echo "=== Detector battery ==="
/home/qmchem_max/micromamba/envs/delfin/bin/python \
  /home/qmchem_max/agent_workspace/quality_framework/scripts/run_all_detectors.py \
  $ARCHIVE_DIR \
  --output $REPORTS_DIR/all_metrics_$LABEL.json \
  --cli-sample 50000 \
  --all \
  --skip-xtb 2>&1 | tail -30

# iter_gate
echo ""
echo "=== iter_gate $LABEL vs $PREV ==="
/home/qmchem_max/micromamba/envs/delfin/bin/python \
  /home/qmchem_max/agent_workspace/quality_framework/scripts/iter_gate.py \
  $LABEL \
  --prev $PREV \
  --reports $REPORTS_DIR
GATE_EXIT=$?
echo "iter_gate exit code: $GATE_EXIT"

# compare_pools --all (North-Star)
echo ""
echo "=== compare_pools --all $LABEL (North-Star vs canonical archives) ==="
/home/qmchem_max/micromamba/envs/delfin/bin/python \
  /home/qmchem_max/agent_workspace/quality_framework/scripts/compare_pools.py \
  $LABEL --all --reports $REPORTS_DIR 2>&1 | tail -50

echo ""
echo "=== D1.v2 post-pool DONE $(date) ==="
echo "iter_gate exit: $GATE_EXIT"
echo "0 = PASS, 2 = HEAL-FIRST (net+ but >=1 severe), 1 = FAIL"
exit $GATE_EXIT
