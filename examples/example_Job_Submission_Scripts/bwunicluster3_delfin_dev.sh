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
# DELFIN auf BwUniCluster 3.0 - Development Installation
# ========================================================================
#
# Nutzt Development-Installation (git clone + pip install -e .)
# statt PyPI-Installation
#
# ========================================================================

# Module laden
module purge
module load chem/orca/6.1.1         # ORCA Version anpassen
module load devel/python/3.11       # Python Version anpassen

# DELFIN Development Environment aktivieren
source ~/DELFIN/.venv/bin/activate  # Pfad zu deinem geklonten Repo

# Umgebungsvariablen
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Scratch-Verzeichnis
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

# Job Information
echo "========================================="
echo "DELFIN auf BwUniCluster 3.0 (Dev)"
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

# Git commit info für Reproduzierbarkeit
cd ~/DELFIN
echo "Git Commit: $(git rev-parse --short HEAD)"
echo "Git Branch: $(git rev-parse --abbrev-ref HEAD)"
echo ""
cd - > /dev/null

# DELFIN ausführen
delfin

EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job beendet:     $(date)"
echo "Exit Code:       $EXIT_CODE"
echo "========================================="

# Cleanup
rm -rf "$DELFIN_SCRATCH"
rm -rf "$ORCA_TMP"

exit $EXIT_CODE
