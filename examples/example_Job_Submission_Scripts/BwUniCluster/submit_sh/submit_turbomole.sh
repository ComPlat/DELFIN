#!/bin/bash
#SBATCH --job-name=turbomole_job
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --threads-per-core=1
#SBATCH --mem=240G
#SBATCH --time=48:00:00
#SBATCH --output=turbomole_%j.out
#SBATCH --error=turbomole_%j.err
#SBATCH --signal=B:SIGTERM@120

# ========================================================================
# TURBOMOLE Central Submit Script for BwUniCluster
# ========================================================================
#
# ENVIRONMENT VARIABLES:
#   TM_JOB_NAME:     Job name for display
#   TM_MODULE:       TURBOMOLE module to run (ridft, dscf, ricc2, jobex, etc.)
#   TM_NPROCS:       Number of processors (default: from SLURM)
#   TM_PARA_ARCH:    Parallelization: SMP or MPI (default: SMP)
#   TM_RUN_DEFINE:   If "1", run define with TM_DEFINE_INPUT before calculation
#   TM_DEFINE_INPUT: Input for define (newline-separated commands)
#
# RESOURCE PARAMETERS (via sbatch command-line overrides):
#   --time, --cpus-per-task, --mem, --job-name
#
# ========================================================================

set -euo pipefail

# === Configuration from environment variables ===
DISPLAY_JOB_NAME="${TM_JOB_NAME:-${SLURM_JOB_NAME:-turbomole_job}}"
TM_MODULE="${TM_MODULE:-ridft}"
TM_PARA_ARCH="${TM_PARA_ARCH:-SMP}"
TM_RUN_DEFINE="${TM_RUN_DEFINE:-0}"
TM_DEFINE_INPUT="${TM_DEFINE_INPUT:-}"

# === Load modules ===
module purge
echo "Loading TURBOMOLE module..."
export TURBOMOLE_MODE="compute"

# Set parallelization before loading module
NPROCS="${TM_NPROCS:-${SLURM_CPUS_PER_TASK:-1}}"
if [ "$NPROCS" -gt 1 ]; then
    export PARA_ARCH="$TM_PARA_ARCH"
    export PARNODES="$NPROCS"
else
    unset PARA_ARCH
    unset PARNODES
fi

module load chem/turbomole/7.9

if [ -z "${TURBOMOLE_VERSION:-}" ]; then
    echo "ERROR: Could not load TURBOMOLE module"
    exit 1
fi

# === Environment setup ===
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
unset LANG
unset LC_CTYPE

# === Scratch directory setup ===
if [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
    runDIR="${TMPDIR}/turbomole_${SLURM_JOB_ID}"
    echo "Using local SSD (TMPDIR): $TMPDIR"
else
    runDIR="/scratch/${USER}/turbomole_${SLURM_JOB_ID}"
    echo "WARNING: Using /scratch (network filesystem)"
fi

# Set TURBOTMPDIR based on PARA_ARCH
if [ "${PARA_ARCH:-}" = "MPI" ]; then
    export TURBOTMPDIR="${runDIR}/run-dir"
else
    export TURBOTMPDIR="${runDIR}"
fi

mkdir -p "$runDIR"

# === Job info ===
echo "========================================"
echo "TURBOMOLE Job - $DISPLAY_JOB_NAME"
echo "========================================"
echo "Job ID:           $SLURM_JOB_ID"
echo "Node:             $SLURM_JOB_NODELIST"
echo "CPUs:             ${SLURM_CPUS_PER_TASK:-1}"
echo "Memory:           ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "Submit Dir:       $SLURM_SUBMIT_DIR"
echo "Run Dir:          $runDIR"
echo "TURBOMOLE Version: $TURBOMOLE_VERSION"
echo "TURBODIR:         $TURBODIR"
echo "PARA_ARCH:        ${PARA_ARCH:-sequential}"
echo "PARNODES:         ${PARNODES:-1}"
echo "TURBOTMPDIR:      $TURBOTMPDIR"
echo "Module to run:    $TM_MODULE"
echo "Started:          $(date)"
echo "========================================"
echo ""

# === Cleanup function ===
cleanup() {
    local signal_name="${1:-UNKNOWN}"
    echo ""
    echo "========================================"
    echo "Caught $signal_name signal, running cleanup..."
    echo "========================================"

    # Copy results back
    echo "Copying results back to $SLURM_SUBMIT_DIR..."
    if [ -d "$runDIR" ]; then
        # Remove scratch files before copying
        rm -f "$runDIR"/{dens,errvec,fock,oldfock,slave*} 2>/dev/null || true
        rsync -a "$runDIR"/ "$SLURM_SUBMIT_DIR"/ 2>/dev/null || true
        echo "Results copied successfully."
    fi

    rm -rf "$runDIR" 2>/dev/null || true
    echo "Cleanup completed at $(date)"
    exit 1
}

trap 'cleanup SIGTERM' SIGTERM
trap 'cleanup SIGINT' SIGINT
trap 'cleanup SIGHUP' SIGHUP

# === Copy input files to scratch ===
echo "Copying input files to scratch..."
cp -a "$SLURM_SUBMIT_DIR"/* "$runDIR"/ 2>/dev/null || true
rm -f "$runDIR"/turbomole_*.out "$runDIR"/turbomole_*.err 2>/dev/null || true

cd "$runDIR"

# === Run define if requested ===
if [ "$TM_RUN_DEFINE" = "1" ]; then
    # Try to read define input from file first, then from environment variable
    if [ -f "define_input.txt" ]; then
        echo "Running define with input from define_input.txt..."
        define < define_input.txt > define.out 2>&1
        DEFINE_EXIT=$?
    elif [ -n "$TM_DEFINE_INPUT" ]; then
        echo "Running define with input from environment..."
        echo -e "$TM_DEFINE_INPUT" | define > define.out 2>&1
        DEFINE_EXIT=$?
    else
        echo "WARNING: TM_RUN_DEFINE=1 but no define_input.txt or TM_DEFINE_INPUT found"
        DEFINE_EXIT=1
    fi

    if [ $DEFINE_EXIT -ne 0 ]; then
        echo "ERROR: define failed with exit code $DEFINE_EXIT"
        echo "=== define output ==="
        cat define.out
        echo "===================="
        cleanup "DEFINE_FAILED"
    fi
    echo "define completed successfully."
    echo ""
fi

# === Check for required files ===
if [ ! -f "control" ]; then
    echo "ERROR: No control file found in job directory."
    echo "Either run define first or provide a control file."
    cleanup "NO_CONTROL"
fi

if [ ! -f "coord" ]; then
    echo "ERROR: No coord file found in job directory."
    cleanup "NO_COORD"
fi

# === Run TURBOMOLE calculation ===
echo "Starting TURBOMOLE calculation: $TM_MODULE"
echo "========================================"
EXIT_CODE=0
set +e

case "$TM_MODULE" in
    ridft)
        time ridft > ridft.out 2>&1
        EXIT_CODE=$?
        ;;
    dscf)
        time dscf > dscf.out 2>&1
        EXIT_CODE=$?
        ;;
    ricc2)
        time ricc2 > ricc2.out 2>&1
        EXIT_CODE=$?
        ;;
    escf)
        time escf > escf.out 2>&1
        EXIT_CODE=$?
        ;;
    grad)
        time grad > grad.out 2>&1
        EXIT_CODE=$?
        ;;
    rdgrad)
        time rdgrad > rdgrad.out 2>&1
        EXIT_CODE=$?
        ;;
    rigrad)
        time rigrad > rigrad.out 2>&1
        EXIT_CODE=$?
        ;;
    jobex)
        # Geometry optimization
        if [ "${PARA_ARCH:-}" = "SMP" ] || [ -z "${PARA_ARCH:-}" ]; then
            time jobex -ri > jobex.out 2>&1
        else
            time jobex -ri -mpi > jobex.out 2>&1
        fi
        EXIT_CODE=$?
        ;;
    jobex-noopt)
        # Single point with jobex (no optimization)
        if [ "${PARA_ARCH:-}" = "SMP" ] || [ -z "${PARA_ARCH:-}" ]; then
            time jobex -ri -energy > jobex.out 2>&1
        else
            time jobex -ri -energy -mpi > jobex.out 2>&1
        fi
        EXIT_CODE=$?
        ;;
    aoforce)
        time aoforce > aoforce.out 2>&1
        EXIT_CODE=$?
        ;;
    NumForce)
        time NumForce -ri > numforce.out 2>&1
        EXIT_CODE=$?
        ;;
    freeh)
        time freeh > freeh.out 2>&1
        EXIT_CODE=$?
        ;;
    ccsdf12)
        time ccsdf12 > ccsdf12.out 2>&1
        EXIT_CODE=$?
        ;;
    pnoccsd)
        time pnoccsd > pnoccsd.out 2>&1
        EXIT_CODE=$?
        ;;
    custom)
        # Custom command sequence from TM_CUSTOM_CMD
        if [ -n "${TM_CUSTOM_CMD:-}" ]; then
            echo "Running custom command: $TM_CUSTOM_CMD"
            eval "$TM_CUSTOM_CMD" > custom.out 2>&1
            EXIT_CODE=$?
        else
            echo "ERROR: TM_MODULE=custom but TM_CUSTOM_CMD not set"
            EXIT_CODE=1
        fi
        ;;
    *)
        # Try to run the module directly
        echo "Running: $TM_MODULE"
        time $TM_MODULE > "${TM_MODULE}.out" 2>&1
        EXIT_CODE=$?
        ;;
esac

set -e

echo ""
echo "========================================"
echo "TURBOMOLE calculation finished: $(date)"
echo "Exit Code: $EXIT_CODE"
echo "========================================"

# === Cleanup scratch files ===
echo "Cleaning up scratch files..."
rm -f dens errvec fock oldfock slave* 2>/dev/null || true

# === Copy results back ===
echo "Copying results back to $SLURM_SUBMIT_DIR..."
rsync -a "$runDIR"/ "$SLURM_SUBMIT_DIR"/

# === Create result archive (optional) ===
tararchive="job_tm_${SLURM_JOB_ID}.tgz"
echo "Creating result archive: $tararchive"
cd "$SLURM_SUBMIT_DIR"
tar -czf "$tararchive" --exclude="$tararchive" --exclude="turbomole_*.out" --exclude="turbomole_*.err" . 2>/dev/null || true

# === Cleanup scratch ===
rm -rf "$runDIR" 2>/dev/null || true

echo ""
echo "========================================"
echo "Job completed successfully!"
echo "Results in: $SLURM_SUBMIT_DIR"
echo "Archive: $tararchive"
echo "========================================"

exit $EXIT_CODE
