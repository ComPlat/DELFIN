#!/bin/bash
#SBATCH --job-name=delfin_highmem
#SBATCH --partition=highmem          # High-Memory Partition für große Systeme
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --threads-per-core=1
#SBATCH --mem=500G                   # High-Memory: bis zu 1TB möglich
#SBATCH --time=72:00:00              # Längere Laufzeit für große Systeme
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --constraint=BEEOND
#SBATCH --exclusive

# ========================================================================
# DELFIN auf BwUniCluster 3.0 - High-Memory Variante
# ========================================================================
#
# Für große Systeme mit hohem Speicherbedarf:
# - Große Moleküle (>200 Atome)
# - Korrelierte Methoden (CCSD(T), NEVPT2)
# - Basis-Sets mit vielen Funktionen
#
# ========================================================================

# Module laden
module purge
module load chem/orca/6.1.1
module load devel/python/3.11

# DELFIN Environment
source ~/delfin_env/bin/activate

# Umgebung
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Scratch
if [ -n "$BEEOND_MOUNTPOINT" ]; then
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID}"
else
    export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="/scratch/${USER}/orca_${SLURM_JOB_ID}"
fi

mkdir -p "$DELFIN_SCRATCH"
mkdir -p "$ORCA_TMP"

echo "========================================="
echo "DELFIN High-Memory Job"
echo "========================================="
echo "Job ID:      $SLURM_JOB_ID"
echo "Node:        $SLURM_JOB_NODELIST"
echo "CPUs:        $SLURM_CPUS_PER_TASK"
echo "Memory:      $SLURM_MEM_PER_NODE MB"
echo "Started:     $(date)"
echo "========================================="
echo ""

# DELFIN ausführen
delfin

EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job beendet: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "========================================="

# Cleanup
rm -rf "$DELFIN_SCRATCH"
rm -rf "$ORCA_TMP"

exit $EXIT_CODE
