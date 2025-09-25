#!/bin/bash
#BSUB -J delfin_run
#BSUB -W 8:00
#BSUB -n 12
#BSUB -R "rusage[mem=1300]"
#BSUB -q cpu
#BSUB -o delfin_%J.out
#BSUB -e delfin_%J.err

# Optional: load required modules
# module load python/3.11
# module load orca/6.1.0

# Activate Python environment
# source /path/to/venv/bin/activate

# Ensure DELFIN can place intermediates on the node-local scratch
export DELFIN_SCRATCH="$TMPDIR/delfin_${LSB_JOBID}"
mkdir -p "$DELFIN_SCRATCH"

# Copy CONTROL and related inputs to the compute node if necessary
# cp /project/path/CONTROL.txt $TMPDIR/
# cp /project/path/input.xyz $TMPDIR/
# cd $TMPDIR

# Run DELFIN programmatically; adapt CONTROL path as needed
python - <<'PY'
from delfin.api import run

run(
    control_file="CONTROL.txt",  # adjust path if CONTROL is elsewhere
    cleanup=True,
    recalc=False,
)
PY

# Optionally move results back to project storage
# cp -r $TMPDIR/*.out /project/path/results/