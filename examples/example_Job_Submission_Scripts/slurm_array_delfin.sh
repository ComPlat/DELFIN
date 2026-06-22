#!/bin/bash
#SBATCH --job-name=delfin_array
#SBATCH --array=1-10%3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --signal=B:USR1@300          # Send SIGUSR1 5min before walltime for preemptive sync
#SBATCH --partition=normal
#SBATCH --output=delfin_array_%A_%a.out
#SBATCH --error=delfin_array_%A_%a.err

# ========================================================================
# DELFIN on SLURM - Job Array for Multiple Systems
# ========================================================================
#
# This script runs DELFIN on multiple systems in parallel using SLURM
# job arrays. Perfect for high-throughput screening!
#
# Usage:
# 1. Create directories: system_1/, system_2/, ..., system_10/
# 2. Each directory contains CONTROL.txt and input files
# 3. Submit: sbatch slurm_array_delfin.sh
#
# Key features:
# - Up to 10 systems processed in parallel
# - Maximum 3 jobs running simultaneously (%3)
# - Each system gets full node (24 cores)
# - Optimal for parameter scans / screenings
# - Stage-in/stage-out: all I/O on node-local SSD ($TMPDIR)
#
# ========================================================================

# Load modules
module purge
module load orca/6.0.0
module load python/3.11

# Activate environment
# source activate delfin

# Set environment
export OMP_NUM_THREADS=1

# Array of system directories
SYSTEMS=(
    "system_1"
    "system_2"
    "system_3"
    "system_4"
    "system_5"
    "system_6"
    "system_7"
    "system_8"
    "system_9"
    "system_10"
)

# Get current system from array task ID
SYSTEM_DIR=${SYSTEMS[$SLURM_ARRAY_TASK_ID - 1]}

echo "========================================="
echo "DELFIN Job Array"
echo "========================================="
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "System: $SYSTEM_DIR"
echo "Node: $SLURM_JOB_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "TMPDIR: ${TMPDIR:-not set}"
echo "Started: $(date)"
echo "========================================="
echo ""

# Check if directory exists
if [ ! -d "$SYSTEM_DIR" ]; then
    echo "ERROR: Directory $SYSTEM_DIR not found!"
    exit 1
fi

# Go to system directory
cd $SYSTEM_DIR
ORIGIN_DIR="$(pwd)"

# ------------------------------------------------------------------
# Stage-in: copy job data to node-local SSD
# ------------------------------------------------------------------
STAGE_IO="${DELFIN_STAGE_IO:-1}"
SYNC_INTERVAL="${DELFIN_SYNC_INTERVAL:-900}"
WORK_DIR="$ORIGIN_DIR"

if [ "$STAGE_IO" = "1" ] && [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
  WORK_DIR="$TMPDIR/delfin_${SYSTEM_DIR}"
  echo "[stage-in] Copying $ORIGIN_DIR -> $WORK_DIR"
  mkdir -p "$WORK_DIR"
  _STAGE_OK=0
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --copy-links "$ORIGIN_DIR/" "$WORK_DIR/" && _STAGE_OK=1
  else
    cp -rL "$ORIGIN_DIR/." "$WORK_DIR/" && _STAGE_OK=1
  fi
  if [ "$_STAGE_OK" = "1" ]; then
    export DELFIN_SCRATCH="$TMPDIR"
    echo "[stage-in] Done. Working on local SSD."
    cd "$WORK_DIR"
  else
    echo "[stage-in] WARNING: copy failed (disk full?). Running directly on HOME."
    rm -rf "$WORK_DIR" 2>/dev/null || true
    WORK_DIR="$ORIGIN_DIR"
  fi
fi

# Periodic background sync (excludes temp files to minimize HOME I/O)
_SYNC_PID=""
_DELFIN_PID=""
RSYNC_EXCLUDES="--exclude=.orca_iso_* --exclude=*.tmp --exclude=__pycache__ --exclude=.exit_code_*"

_start_periodic_sync() {
  if [ "$WORK_DIR" = "$ORIGIN_DIR" ]; then return; fi
  if [ "$SYNC_INTERVAL" -le 0 ] 2>/dev/null; then return; fi
  (
    while true; do
      sleep "$SYNC_INTERVAL"
      rsync -a --update $RSYNC_EXCLUDES "$WORK_DIR/" "$ORIGIN_DIR/" 2>/dev/null || true
      echo "[sync] Periodic sync completed at $(date +%H:%M:%S)"
    done
  ) &
  _SYNC_PID=$!
  echo "[sync] Periodic background sync every ${SYNC_INTERVAL}s (PID $_SYNC_PID)"
}

_stop_periodic_sync() {
  if [ -n "$_SYNC_PID" ]; then
    kill "$_SYNC_PID" 2>/dev/null || true
    wait "$_SYNC_PID" 2>/dev/null || true
    _SYNC_PID=""
  fi
}

_stage_out() {
  set +e
  local rc=${1:-$?}
  _stop_periodic_sync
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
  exit "$rc"
}

_handle_signal() {
  echo ""
  echo "[signal] Received termination signal — stopping DELFIN and saving results..."
  if [ -n "$_DELFIN_PID" ]; then
    kill -TERM -- -"$_DELFIN_PID" 2>/dev/null || kill -TERM "$_DELFIN_PID" 2>/dev/null || true
    wait "$_DELFIN_PID" 2>/dev/null || true
  fi
  _stage_out 124
}

_handle_early_warning() {
  echo ""
  echo "[warning] Approaching walltime -- performing preemptive sync..."
  if [ "$WORK_DIR" != "$ORIGIN_DIR" ]; then
    rsync -a --update $RSYNC_EXCLUDES "$WORK_DIR/" "$ORIGIN_DIR/" 2>/dev/null || true
    echo "[warning] Preemptive sync done."
  fi
}

trap '_stage_out $?' EXIT
trap '_handle_signal' SIGTERM SIGINT SIGHUP
trap '_handle_early_warning' SIGUSR1

# ------------------------------------------------------------------
# Run DELFIN
# ------------------------------------------------------------------
if [ "$WORK_DIR" = "$ORIGIN_DIR" ]; then
  delfin
else
  _start_periodic_sync
  delfin &
  _DELFIN_PID=$!
  wait "$_DELFIN_PID"
fi

EXIT_CODE=$?

echo ""
echo "========================================="
echo "System $SYSTEM_DIR finished at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================="

exit $EXIT_CODE
