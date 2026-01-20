#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --constraint=BEEOND
#SBATCH --exclusive

# ========================================================================
# DELFIN Short Job (24h)
# ========================================================================

set -euo pipefail

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

# MPI environment for BwUniCluster (Mellanox InfiniBand)
# UCX for optimal performance with mlx5 HCAs
export OMPI_MCA_pml=ucx
export OMPI_MCA_btl=^vader,tcp,openib
export UCX_NET_DEVICES=mlx5_0:1

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

# Cleanup function for trap
cleanup() {
    echo ""
    echo "========================================"
    echo "Caught signal, running cleanup..."
    echo "========================================"
    cd "$RUN_DIR" 2>/dev/null && delfin --cleanup 2>/dev/null || true
    # Copy results back before cleanup
    cp -a "$RUN_DIR"/* "$SLURM_SUBMIT_DIR"/ 2>/dev/null || true
    rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR" 2>/dev/null || true
    echo "Cleanup completed."
    exit 1
}

# Set trap for cleanup on termination signals
trap cleanup SIGTERM SIGINT SIGHUP

# Job info
echo "========================================"
echo "DELFIN Job - {job_name}"
echo "Job Type: Short (24h)"
echo "========================================"
echo "Job ID:      $SLURM_JOB_ID"
echo "Node:        $SLURM_JOB_NODELIST"
echo "CPUs:        $SLURM_NTASKS"
echo "Memory:      ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "Submit Dir:  $SLURM_SUBMIT_DIR"
echo "Scratch:     $DELFIN_SCRATCH"
echo "Started:     $(date)"
echo "========================================"
echo ""

# Check OpenMPI version
echo "OpenMPI: $(which mpirun)"
mpirun --version 2>&1 | head -1 || echo "MPI check failed"
echo ""

# Check ORCA version
echo "ORCA: $(which orca)"
orca --version 2>&1 | head -5 || echo "ORCA check failed"
echo ""

# Check DELFIN version
cd "$DELFIN_DIR"
echo "DELFIN Version: $(delfin --version 2>&1 || echo 'unknown')"
echo "Git Branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'N/A')"
echo "Git Commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'N/A')"
echo ""
cd - > /dev/null

# Verify required input files exist
if [ ! -f "$SLURM_SUBMIT_DIR/CONTROL.txt" ]; then
    echo "ERROR: CONTROL.txt not found in $SLURM_SUBMIT_DIR"
    exit 1
fi
if [ ! -f "$SLURM_SUBMIT_DIR/input.txt" ]; then
    echo "ERROR: input.txt not found in $SLURM_SUBMIT_DIR"
    exit 1
fi

# Copy inputs to scratch and work there
cp -a "$SLURM_SUBMIT_DIR/CONTROL.txt" "$SLURM_SUBMIT_DIR/input.txt" "$RUN_DIR"/

cd "$RUN_DIR"

# Run DELFIN
delfin

EXIT_CODE=$?

echo ""
echo "========================================"
echo "Job finished: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "========================================"

# Copy results back
cp -a "$RUN_DIR"/* "$SLURM_SUBMIT_DIR"/

# Cleanup scratch
rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR"

exit $EXIT_CODE
