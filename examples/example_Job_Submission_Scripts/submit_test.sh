#!/bin/bash
#SBATCH --job-name=delfin_test
#SBATCH --partition=dev_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err

# ========================================================================
# DELFIN Test Job - BwUniCluster 3.0
# ========================================================================
#
# Für schnelle Tests und Debugging:
# - 8 CPUs, 16 GB RAM, 30 Minuten
# - Development Queue (schneller Durchsatz)
# - Ohne BeeOND (schneller Start)
#
# Verwendung:
#   1. In dein Projektverzeichnis kopieren
#   2. Module-Namen anpassen (Zeile 27-28)
#   3. sbatch submit_test.sh
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
export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
export ORCA_TMP="/scratch/${USER}/orca_${SLURM_JOB_ID}"
mkdir -p "$DELFIN_SCRATCH"
mkdir -p "$ORCA_TMP"

# Job Info
echo "========================================"
echo "DELFIN Test Job"
echo "========================================"
echo "Job ID:      $SLURM_JOB_ID"
echo "Node:        $SLURM_JOB_NODELIST"
echo "CPUs:        $SLURM_CPUS_PER_TASK"
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
