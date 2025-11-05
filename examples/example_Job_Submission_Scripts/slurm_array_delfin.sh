#!/bin/bash
#SBATCH --job-name=delfin_array
#SBATCH --array=1-10%3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=48:00:00
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

# Create scratch
export ORCA_TMP=/scratch/$SLURM_JOB_ID
mkdir -p $ORCA_TMP

# Run DELFIN
delfin

EXIT_CODE=$?

echo ""
echo "========================================="
echo "System $SYSTEM_DIR finished at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================="

# Cleanup
rm -rf $ORCA_TMP

exit $EXIT_CODE
