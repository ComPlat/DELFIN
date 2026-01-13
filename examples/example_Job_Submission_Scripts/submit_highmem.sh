#!/bin/bash
#SBATCH --job-name=delfin_highmem
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --threads-per-core=1
#SBATCH --mem=500G
#SBATCH --time=72:00:00
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --constraint=BEEOND
#SBATCH --exclusive

# ========================================================================
# DELFIN High-Memory Job - BwUniCluster 3.0
# ========================================================================
#
# Für große Systeme mit hohem Speicherbedarf:
# - 40 CPUs, 500 GB RAM (bis 1 TB möglich), 72 Stunden
# - Große Moleküle (>200 Atome)
# - Korrelierte Methoden (CCSD(T), NEVPT2)
# - Große Basis-Sets
#
# Verwendung:
#   1. In dein Projektverzeichnis kopieren
#   2. Module-Namen anpassen (Zeile 29-30)
#   3. sbatch submit_highmem.sh
#
# ========================================================================

# Module laden - HIER ANPASSEN!
module purge
module load chem/orca/6.1.1         # Mit 'module avail orca' prüfen
module load devel/python/3.11       # Mit 'module avail python' prüfen

# DELFIN Environment aktivieren
source ~/DELFIN/.venv/bin/activate

# Umgebungsvariablen
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Scratch-Verzeichnis
if [ -n "$BEEOND_MOUNTPOINT" ]; then
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID}"
    echo "Using BeeOND: $BEEOND_MOUNTPOINT"
else
    export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="/scratch/${USER}/orca_${SLURM_JOB_ID}"
fi

mkdir -p "$DELFIN_SCRATCH"
mkdir -p "$ORCA_TMP"

# Job Info
echo "========================================"
echo "DELFIN High-Memory Job"
echo "========================================"
echo "Job ID:      $SLURM_JOB_ID"
echo "Node:        $SLURM_JOB_NODELIST"
echo "CPUs:        $SLURM_CPUS_PER_TASK"
echo "Memory:      $SLURM_MEM_PER_NODE MB"
echo "Working Dir: $(pwd)"
echo "Scratch:     $DELFIN_SCRATCH"
echo "Started:     $(date)"
echo "========================================"
echo ""

# DELFIN ausführen
delfin

EXIT_CODE=$?

echo ""
echo "========================================"
echo "Job beendet: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "========================================"

# Cleanup
rm -rf "$DELFIN_SCRATCH"
rm -rf "$ORCA_TMP"

exit $EXIT_CODE
