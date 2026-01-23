#!/bin/bash
#SBATCH --job-name=delfin_job
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem=240G
#SBATCH --time=48:00:00
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --constraint=BEEOND
#SBATCH --exclusive
#SBATCH --signal=B:SIGTERM@120

# ========================================================================
# DELFIN Central Submit Script
# Supports all modes via environment variables:
#   DELFIN_MODE: delfin | delfin-recalc | orca | auto (default: auto)
#   DELFIN_INP_FILE: Specific .inp file for ORCA mode
#   DELFIN_JOB_NAME: Job name for display (optional)
#
# Resource parameters are passed via sbatch command-line overrides:
#   --time, --ntasks, --mem, --job-name
# ========================================================================

set -euo pipefail

# === Configuration from environment variables ===
MODE="${DELFIN_MODE:-auto}"
INP_FILE="${DELFIN_INP_FILE:-}"
DISPLAY_JOB_NAME="${DELFIN_JOB_NAME:-${SLURM_JOB_NAME:-delfin_job}}"

# Load modules
module purge
module load devel/python/3.11.7-gnu-14.2

# Initialize LD_LIBRARY_PATH if not set
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# Auto-detect base directory (search upwards for software/delfin)
if [ -z "${SLURM_SUBMIT_DIR:-}" ]; then
    echo "ERROR: SLURM_SUBMIT_DIR not set; cannot locate BASE_DIR"
    exit 1
fi
SEARCH_DIR="$SLURM_SUBMIT_DIR"
while [ "$SEARCH_DIR" != "/" ] && [ ! -d "$SEARCH_DIR/software/delfin" ]; do
    SEARCH_DIR="$(dirname "$SEARCH_DIR")"
done
if [ ! -d "$SEARCH_DIR/software/delfin" ]; then
    echo "ERROR: Could not locate software/delfin above $SLURM_SUBMIT_DIR"
    exit 1
fi
BASE_DIR="$SEARCH_DIR/software"
DELFIN_DIR="$BASE_DIR/delfin"

# Use custom OpenMPI 4.1.8 (compatible with ORCA)
if [ ! -d "$BASE_DIR/openmpi-4.1.8" ]; then
    echo "ERROR: OpenMPI 4.1.8 not found in $BASE_DIR/openmpi-4.1.8"
    echo "Please install it first. See installation instructions."
    exit 1
fi

echo "Using custom OpenMPI 4.1.8 from $BASE_DIR/openmpi-4.1.8"
export PATH="$BASE_DIR/openmpi-4.1.8/bin:$PATH"
export LD_LIBRARY_PATH="$BASE_DIR/openmpi-4.1.8/lib:$LD_LIBRARY_PATH"

# Set ORCA path
ORCA_BASE="$BASE_DIR/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2"
if [ ! -d "$ORCA_BASE" ]; then
    echo "ERROR: ORCA not found in $ORCA_BASE"
    echo "Please install ORCA 6.1.1 or update ORCA_BASE."
    exit 1
fi
export PATH="$ORCA_BASE:$PATH"
export LD_LIBRARY_PATH="$ORCA_BASE:$LD_LIBRARY_PATH"
export ORCA_DIR="$ORCA_BASE"
export ORCA_PLOT="$ORCA_BASE/orca_plot"

# =============================================================================
# MPI Configuration for ORCA on BwUniCluster
# Community-recommended stable settings (ORCA Input Library)
# =============================================================================

# Use ob1 PML with vader (shared memory) - sm was removed in OpenMPI 3.0
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp,vader

# OpenMPI settings optimized for ORCA stability
export OMPI_MCA_mpi_show_mca_params_file=0
export OMPI_MCA_mpi_yield_when_idle=1         # Reduce CPU spinning
export OMPI_MCA_coll_hcoll_enable=0           # Disable HCOLL (often problematic)

# Clear any inherited TCP interface include/exclude to avoid invalid if_inexclude warnings
unset OMPI_MCA_btl_tcp_if_include
unset OMPI_MCA_btl_tcp_if_exclude

# Process binding: let ORCA handle its own binding
export OMPI_MCA_hwloc_base_binding_policy=none
export OMPI_MCA_rmaps_base_mapping_policy=core
export OMPI_MCA_rmaps_base_oversubscribe=true # Allow ORCA's dynamic parallelism

# Activate DELFIN virtual environment
source "$DELFIN_DIR/.venv/bin/activate"

# Environment variables
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export DELFIN_ORCA_PROGRESS=0
export MPLBACKEND=Agg

# Scratch directory setup
if [ -n "${BEEOND_MOUNTPOINT:-}" ]; then
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID}"
    export ORCA_TMPDIR="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID}"
    echo "Using BeeOND: $BEEOND_MOUNTPOINT"
else
    export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
    export ORCA_TMPDIR="/scratch/${USER}/orca_${SLURM_JOB_ID}"
fi

RUN_DIR="$DELFIN_SCRATCH/run"
mkdir -p "$RUN_DIR" "$ORCA_TMPDIR"

# Cleanup function for trap (handles SIGTERM from timeout, SIGINT, etc.)
cleanup() {
    local signal_name="${1:-UNKNOWN}"
    echo ""
    echo "========================================"
    echo "Caught $signal_name signal, running cleanup..."
    echo "========================================"

    # Try DELFIN cleanup first (safe if not running DELFIN)
    cd "$RUN_DIR" 2>/dev/null && delfin --cleanup 2>/dev/null || true

    # CRITICAL: Copy ALL results back before cleanup
    echo "Copying results back to $SLURM_SUBMIT_DIR..."
    if [ -d "$RUN_DIR" ]; then
        # Remove useless .tmp files before copying (saves space and time)
        find "$RUN_DIR" -name "*.tmp" -delete 2>/dev/null || true
        cp -a "$RUN_DIR"/* "$SLURM_SUBMIT_DIR"/ 2>/dev/null || true
        echo "Results copied successfully."
    else
        echo "WARNING: RUN_DIR not found, nothing to copy."
    fi

    # Cleanup scratch (only after successful copy)
    rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR" 2>/dev/null || true
    echo "Cleanup completed at $(date)"
    exit 1
}

# Set trap for cleanup on termination signals (including SLURM timeout)
trap 'cleanup SIGTERM' SIGTERM
trap 'cleanup SIGINT' SIGINT
trap 'cleanup SIGHUP' SIGHUP

# Job info
echo "========================================"
echo "DELFIN Job - $DISPLAY_JOB_NAME"
echo "Mode: $MODE"
echo "========================================"
echo "Job ID:      $SLURM_JOB_ID"
echo "Node:        $SLURM_JOB_NODELIST"
echo "CPUs:        $SLURM_NTASKS"
echo "Memory:      ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "Time Limit:  ${SLURM_TIMELIMIT:-unknown}"
echo "Submit Dir:  $SLURM_SUBMIT_DIR"
echo "Scratch:     $DELFIN_SCRATCH"
echo "Started:     $(date)"
if [ -n "$INP_FILE" ]; then
    echo "Input File:  $INP_FILE"
fi
echo "========================================"
echo ""

# Check OpenMPI version
echo "OpenMPI: $(which mpirun)"
mpirun --version 2>&1 | head -1 || echo "MPI check failed"
echo ""

# Check ORCA version
orca_path=$(command -v orca || true)
echo "ORCA: ${orca_path:-not found}"
if [ -n "$orca_path" ]; then
    orca --version 2>&1 | head -5 || true
fi
echo ""

# Check DELFIN version
cd "$DELFIN_DIR"
echo "DELFIN Version: $(delfin --version 2>&1 || echo 'unknown')"
echo "Git Branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'N/A')"
echo "Git Commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'N/A')"
echo ""
cd - > /dev/null

# Copy ALL input files to scratch
echo "Copying input files to scratch..."
cp -a "$SLURM_SUBMIT_DIR"/* "$RUN_DIR"/ 2>/dev/null || true
# Remove output files from previous runs (if any)
rm -f "$RUN_DIR"/delfin_*.out "$RUN_DIR"/delfin_*.err 2>/dev/null || true

cd "$RUN_DIR"

# Auto-detect mode if set to "auto"
if [ "$MODE" = "auto" ]; then
    if [ -f "$SLURM_SUBMIT_DIR/CONTROL.txt" ] && [ -f "$SLURM_SUBMIT_DIR/input.txt" ]; then
        MODE="delfin"
        echo "Auto-detected mode: DELFIN (CONTROL.txt + input.txt found)"
    elif ls "$SLURM_SUBMIT_DIR"/*.inp 1>/dev/null 2>&1; then
        MODE="orca"
        echo "Auto-detected mode: ORCA (*.inp files found, no CONTROL.txt)"
    else
        echo "ERROR: No valid input files found in $SLURM_SUBMIT_DIR"
        echo "       Expected either: CONTROL.txt + input.txt (DELFIN mode)"
        echo "       Or: *.inp files (ORCA-only mode)"
        exit 1
    fi
fi

# Run appropriate mode
EXIT_CODE=0
set +e
case "$MODE" in
    delfin)
        echo "Starting DELFIN..."
        delfin
        EXIT_CODE=$?
        ;;
    delfin-recalc)
        echo "Starting DELFIN --recalc..."
        delfin --recalc
        EXIT_CODE=$?
        ;;
    orca)
        # Use specified input file or find the first one
        if [ -z "$INP_FILE" ]; then
            INP_FILE=$(ls *.inp 2>/dev/null | head -1)
        fi
        if [ -z "$INP_FILE" ]; then
            echo "ERROR: No .inp file found for ORCA mode"
            EXIT_CODE=1
        elif [ ! -f "$INP_FILE" ]; then
            echo "ERROR: Specified input file not found: $INP_FILE"
            EXIT_CODE=1
        else
            OUT_FILE="${INP_FILE%.inp}.out"
            echo "Starting ORCA: $INP_FILE -> $OUT_FILE"
            "$ORCA_BASE/orca" "$INP_FILE" > "$OUT_FILE" 2>&1
            EXIT_CODE=$?
        fi
        ;;
    *)
        echo "ERROR: Unknown mode: $MODE"
        echo "       Valid modes: delfin, delfin-recalc, orca, auto"
        EXIT_CODE=1
        ;;
esac
set -e

echo ""
echo "========================================"
echo "Job finished: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "========================================"

# Remove useless .tmp files before copying (saves space and time)
find "$RUN_DIR" -name "*.tmp" -delete 2>/dev/null || true

# Copy results back
cp -a "$RUN_DIR"/* "$SLURM_SUBMIT_DIR"/

# Cleanup scratch
rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR"

exit $EXIT_CODE
