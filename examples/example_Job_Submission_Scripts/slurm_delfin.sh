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
#
# ========================================================================

# Load required modules (adjust for your cluster)
module purge
module load orca/6.1.0
module load python/3.11
# module load conda  # if using conda environment

# Activate DELFIN environment
# source activate delfin
# OR
# source /path/to/delfin_venv/bin/activate

# Set OMP threads (ORCA uses OpenMP)
export OMP_NUM_THREADS=1  # DELFIN manages parallelism

# Optional: Set scratch directory for ORCA
export ORCA_TMP=/scratch/$SLURM_JOB_ID
mkdir -p $ORCA_TMP

# Print cluster information
echo "========================================="
echo "DELFIN SLURM Job Information"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: $SLURM_MEM_PER_NODE MB"
echo "Working directory: $(pwd)"
echo "Started at: $(date)"
echo "========================================="
echo ""

# Run DELFIN
# DELFIN automatically detects SLURM_CPUS_PER_TASK and uses all cores!
delfin

# Capture exit code
EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job finished at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================="

# Cleanup scratch
rm -rf $ORCA_TMP

exit $EXIT_CODE
