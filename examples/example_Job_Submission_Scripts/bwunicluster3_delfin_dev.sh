#!/bin/bash
#SBATCH --job-name=delfin
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --threads-per-core=1
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --constraint=BEEOND
#SBATCH --exclusive

# ========================================================================
# DELFIN on BwUniCluster 3.0 - Development Installation
# ========================================================================
#
# Uses development installation (git clone + pip install -e .)
# instead of PyPI installation.
# Stage-in/stage-out moves all job I/O to node-local SSD.
#
# Set DELFIN_STAGE_IO=0 to disable staging.
#
# ========================================================================

# Load modules
module purge
module load chem/orca/6.1.1         # Adjust ORCA version
module load devel/python/3.11       # Adjust Python version

# Activate DELFIN development environment
# Copy venv to local SSD to minimize HOME I/O
DELFIN_DIR="$HOME/DELFIN"
VENV_DIR="$DELFIN_DIR/.venv"
VENV_TAR="$DELFIN_DIR/delfin_venv.tar"
if [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
    VENV_LOCAL="$TMPDIR/delfin_venv_${SLURM_JOB_ID}"
    if [ -f "$VENV_TAR" ]; then
        echo "Extracting venv to local SSD ($VENV_LOCAL)..."
        mkdir -p "$VENV_LOCAL"
        tar -xf "$VENV_TAR" --strip-components=1 -C "$VENV_LOCAL"
    else
        echo "WARNING: $VENV_TAR not found, using cp (more HOME I/O)."
        echo "         Run once: cd $DELFIN_DIR && tar -cf delfin_venv.tar .venv/"
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

# Environment variables
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Scratch directory
if [ -n "$BEEOND_MOUNTPOINT" ]; then
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID}"
    echo "Using BeeOND storage: $BEEOND_MOUNTPOINT"
else
    export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="/scratch/${USER}/orca_${SLURM_JOB_ID}"
fi

mkdir -p "$DELFIN_SCRATCH"
mkdir -p "$ORCA_TMP"

# Print job information
echo "========================================="
echo "DELFIN on BwUniCluster 3.0 (Dev)"
echo "========================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Node:            $SLURM_JOB_NODELIST"
echo "CPUs:            $SLURM_CPUS_PER_TASK"
echo "Memory:          $SLURM_MEM_PER_NODE MB"
echo "DELFIN Path:     $(which delfin)"
echo "DELFIN Version:  $(delfin --version)"
echo "Working Dir:     $(pwd)"
echo "Scratch:         $DELFIN_SCRATCH"
echo "Started:         $(date)"
echo "========================================="
echo ""

# Git commit info for reproducibility
cd ~/DELFIN
echo "Git Commit: $(git rev-parse --short HEAD)"
echo "Git Branch: $(git rev-parse --abbrev-ref HEAD)"
echo ""
cd - > /dev/null

# ------------------------------------------------------------------
# Stage-in: copy job data to node-local SSD to avoid HOME filesystem I/O
# ------------------------------------------------------------------
ORIGIN_DIR="$(pwd)"
STAGE_IO="${DELFIN_STAGE_IO:-1}"
WORK_DIR="$ORIGIN_DIR"

if [ "$STAGE_IO" = "1" ] && [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
  WORK_DIR="$TMPDIR/delfin_job"
  echo "[stage-in] Copying job data $ORIGIN_DIR -> $WORK_DIR"
  mkdir -p "$WORK_DIR"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --copy-links "$ORIGIN_DIR/" "$WORK_DIR/"
  else
    cp -rL "$ORIGIN_DIR/." "$WORK_DIR/"
  fi
  echo "[stage-in] Done. Working directory on local SSD."
  cd "$WORK_DIR"
fi

# Stage-out: copy results back to HOME on any exit (success, failure, timeout)
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

trap '_stage_out $?' EXIT
trap '_handle_signal' SIGTERM SIGINT SIGHUP

# Run DELFIN
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
echo "Job finished:    $(date)"
echo "Exit code:       $EXIT_CODE"
echo "========================================="

exit $EXIT_CODE
