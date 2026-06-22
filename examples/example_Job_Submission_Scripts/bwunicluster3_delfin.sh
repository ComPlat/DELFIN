#!/bin/bash
#SBATCH --job-name=delfin
#SBATCH --partition=cpu              # Main partition (or dev_cpu for testing)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40           # 40 CPUs per node (some nodes have 28)
#SBATCH --threads-per-core=1         # Disable hyperthreading (important!)
#SBATCH --mem=120G                   # 120GB memory
#SBATCH --time=48:00:00              # Max runtime
#SBATCH --signal=B:USR1@300          # Send SIGUSR1 5min before walltime for preemptive sync
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --constraint=BEEOND          # Optional: fast local storage
#SBATCH --exclusive                  # Exclusive node access (required for BEEOND)

# ========================================================================
# DELFIN on BwUniCluster 3.0
# ========================================================================
#
# Optimized for BwUniCluster 3.0 hardware and SLURM configuration.
# Stage-in/stage-out moves all job I/O to node-local SSD ($TMPDIR)
# to avoid excessive HOME filesystem load.
#
# Hardware:
# - 40 CPU cores per node (some nodes have 28)
# - Hyperthreading available but disabled here
# - BeeOND for fast local storage
#
# Set DELFIN_STAGE_IO=0 to disable staging and run directly on HOME.
#
# ========================================================================

# Load modules (check names with 'module avail')
module purge
module load chem/orca/6.1.1         # Adjust ORCA version
module load devel/python/3.11       # Adjust Python version

# Activate DELFIN environment
# Copy venv to local SSD to minimize HOME I/O
VENV_DIR="$HOME/delfin_env"
VENV_TAR="$HOME/delfin_venv.tar"
if [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
    VENV_LOCAL="$TMPDIR/delfin_venv_${SLURM_JOB_ID}"
    if [ -f "$VENV_TAR" ]; then
        echo "Extracting venv to local SSD ($VENV_LOCAL)..."
        mkdir -p "$VENV_LOCAL"
        tar -xf "$VENV_TAR" --strip-components=1 -C "$VENV_LOCAL"
    else
        echo "WARNING: $VENV_TAR not found, using cp (more HOME I/O)."
        echo "         Run once: cd $HOME && tar -cf delfin_venv.tar delfin_env/"
        cp -a "$VENV_DIR" "$VENV_LOCAL"
    fi
    sed -i "s|$VENV_DIR|$VENV_LOCAL|g" "$VENV_LOCAL/bin/activate" 2>/dev/null || true
    source "$VENV_LOCAL/bin/activate"
    echo "venv loaded from local SSD."
else
    echo "WARNING: no TMPDIR available, loading venv directly from HOME."
    source "$VENV_DIR/bin/activate"
fi

# Further reduce Python HOME I/O
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

# Set environment variables
export OMP_NUM_THREADS=1            # DELFIN manages parallelism
export MKL_NUM_THREADS=1            # Intel MKL threading

# Scratch directory
# Option 1: BeeOND (fast, when --constraint=BEEOND is set)
if [ -n "$BEEOND_MOUNTPOINT" ]; then
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID}"
    echo "Using BeeOND storage: $BEEOND_MOUNTPOINT"
# Option 2: Standard scratch
else
    export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="/scratch/${USER}/orca_${SLURM_JOB_ID}"
fi

mkdir -p "$DELFIN_SCRATCH"
mkdir -p "$ORCA_TMP"

# ------------------------------------------------------------------
# Stage-in: copy job data to node-local SSD to avoid HOME filesystem I/O
# ------------------------------------------------------------------
ORIGIN_DIR="$(pwd)"
STAGE_IO="${DELFIN_STAGE_IO:-1}"
SYNC_INTERVAL="${DELFIN_SYNC_INTERVAL:-900}"
WORK_DIR="$ORIGIN_DIR"

if [ "$STAGE_IO" = "1" ] && [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
  WORK_DIR="$TMPDIR/delfin_job"
  echo "[stage-in] Copying job data $ORIGIN_DIR -> $WORK_DIR"
  mkdir -p "$WORK_DIR"
  _STAGE_OK=0
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --copy-links "$ORIGIN_DIR/" "$WORK_DIR/" && _STAGE_OK=1
  else
    cp -rL "$ORIGIN_DIR/." "$WORK_DIR/" && _STAGE_OK=1
  fi
  if [ "$_STAGE_OK" = "1" ]; then
    echo "[stage-in] Done. Working directory on local SSD."
    cd "$WORK_DIR"
  else
    echo "[stage-in] WARNING: copy failed (disk full?). Running directly on HOME."
    rm -rf "$WORK_DIR" 2>/dev/null || true
    WORK_DIR="$ORIGIN_DIR"
  fi
fi

# Periodic background sync: incrementally copy results every SYNC_INTERVAL
# seconds so intermediate results survive crashes/OOM kills.
# Excludes temp files to minimize HOME I/O.
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

# Stage-out: copy results back to HOME on any exit (success, failure, timeout)
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
  rm -rf "$DELFIN_SCRATCH" "$ORCA_TMP" "${VENV_LOCAL:-}" 2>/dev/null || true
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

# SIGUSR1: early warning before walltime (--signal=B:USR1@300).
# Sync results while DELFIN keeps running.
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

# Print job information
echo "========================================="
echo "DELFIN on BwUniCluster 3.0"
echo "========================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Node:            $SLURM_JOB_NODELIST"
echo "Partition:       $SLURM_JOB_PARTITION"
echo "CPUs:            $SLURM_CPUS_PER_TASK"
echo "Memory:          $SLURM_MEM_PER_NODE MB"
echo "Working Dir:     $(pwd)"
echo "Origin Dir:      $ORIGIN_DIR"
echo "Scratch:         $DELFIN_SCRATCH"
echo "ORCA Scratch:    $ORCA_TMP"
echo "Stage-IO:        $STAGE_IO"
echo "Started:         $(date)"
echo "========================================="
echo ""

# Run DELFIN
# DELFIN automatically detects SLURM_CPUS_PER_TASK and uses all available cores
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
echo "Job finished:    $(date)"
echo "Exit code:       $EXIT_CODE"
echo "========================================="

exit $EXIT_CODE
