#!/bin/bash
#SBATCH --job-name=delfin
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --partition=normal
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err

# ========================================================================
# DELFIN on SLURM - Optimal Configuration
# ========================================================================
#
# This script runs DELFIN efficiently on SLURM clusters.
# DELFIN will automatically detect SLURM and use all allocated cores!
#
# Key features:
# - Single SLURM job manages all OCCUPIER FoBs
# - Dynamic core allocation (92-98% efficiency)
# - No sub-job spawning (avoids queue overhead)
# - Optimal for 24-48 core nodes
# - Stage-in/stage-out: all I/O on node-local SSD ($TMPDIR)
#
# To disable local staging: export DELFIN_STAGE_IO=0 before sbatch
#
# ========================================================================

# Load required modules (adjust for your cluster)
module purge
module load orca/6.1.1
module load python/3.11
# module load conda  # if using conda environment

# Activate DELFIN environment
# source activate delfin
# OR
# source /path/to/delfin_venv/bin/activate

# Set OMP threads (ORCA uses OpenMP)
export OMP_NUM_THREADS=1  # DELFIN manages parallelism

# Print cluster information
echo "========================================="
echo "DELFIN SLURM Job Information"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: $SLURM_MEM_PER_NODE MB"
echo "Working directory: $(pwd)"
echo "TMPDIR: ${TMPDIR:-not set}"
echo "Started at: $(date)"
echo "========================================="
echo ""

# ------------------------------------------------------------------
# Stage-in: copy job data to node-local SSD
# ------------------------------------------------------------------
ORIGIN_DIR="$(pwd)"
STAGE_IO="${DELFIN_STAGE_IO:-1}"
WORK_DIR="$ORIGIN_DIR"

if [ "$STAGE_IO" = "1" ] && [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
  WORK_DIR="$TMPDIR/delfin_job"
  echo "[stage-in] Copying $ORIGIN_DIR -> $WORK_DIR"
  mkdir -p "$WORK_DIR"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --copy-links "$ORIGIN_DIR/" "$WORK_DIR/"
  else
    cp -rL "$ORIGIN_DIR/." "$WORK_DIR/"
  fi
  export DELFIN_SCRATCH="$TMPDIR"
  echo "[stage-in] Done. Working on local SSD."
  cd "$WORK_DIR"
fi

# Stage-out on any exit (success, failure, timeout, signal)
_DELFIN_PID=""

_stage_out() {
  set +e
  local rc=${1:-$?}
  if [ "$WORK_DIR" != "$ORIGIN_DIR" ]; then
    echo ""
    echo "[stage-out] Copying results $WORK_DIR -> $ORIGIN_DIR"
    if command -v rsync >/dev/null 2>&1; then
      rsync -a --update "$WORK_DIR/" "$ORIGIN_DIR/" 2>&1 || echo "[stage-out] WARNING: rsync errors"
    else
      cp -ru "$WORK_DIR/." "$ORIGIN_DIR/" 2>&1 || echo "[stage-out] WARNING: cp errors"
    fi
    echo "[stage-out] Done."
  fi
  [ -n "${ORCA_TMP:-}" ] && rm -rf "$ORCA_TMP" 2>/dev/null || true
  exit "$rc"
}

# Handle SLURM timeout (SIGTERM) — stop DELFIN, then stage-out results
_handle_signal() {
  echo ""
  echo "[signal] Received termination signal — stopping DELFIN and saving results..."
  if [ -n "$_DELFIN_PID" ]; then
    kill -TERM -- -"$_DELFIN_PID" 2>/dev/null || kill -TERM "$_DELFIN_PID" 2>/dev/null || true
    wait "$_DELFIN_PID" 2>/dev/null || true
  fi
  _stage_out 124
}

trap '_stage_out $?' EXIT
trap '_handle_signal' SIGTERM SIGINT SIGHUP

# ------------------------------------------------------------------
# Run DELFIN
# ------------------------------------------------------------------
if [ "$WORK_DIR" = "$ORIGIN_DIR" ]; then
  delfin
else
  delfin &
  _DELFIN_PID=$!
  wait "$_DELFIN_PID"
fi

EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job finished at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================="

exit $EXIT_CODE
