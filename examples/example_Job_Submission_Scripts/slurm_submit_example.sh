#!/bin/bash
#SBATCH --job-name=delfin_run
#SBATCH --time=08:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err

# Optional: load required modules
# module load python/3.11
# module load orca/6.1.0

# Activate Python environment
# source /path/to/venv/bin/activate

# Ensure DELFIN can place intermediates on the node-local scratch
export DELFIN_SCRATCH="$TMPDIR/delfin_${SLURM_JOB_ID}"
mkdir -p "$DELFIN_SCRATCH"

# Copy CONTROL and related inputs to the compute node if necessary
# cp /project/path/CONTROL.txt $SLURM_TMPDIR/
# cp /project/path/input.xyz $SLURM_TMPDIR/
# cd $SLURM_TMPDIR

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
# cp -r $SLURM_TMPDIR/*.out /project/path/results/
