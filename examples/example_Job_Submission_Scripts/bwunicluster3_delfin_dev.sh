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
# venv auf lokale SSD kopieren um HOME I/O zu minimieren
DELFIN_DIR="$HOME/DELFIN"
VENV_DIR="$DELFIN_DIR/.venv"
VENV_TAR="$DELFIN_DIR/delfin_venv.tar"
if [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
    VENV_LOCAL="$TMPDIR/delfin_venv_${SLURM_JOB_ID}"
    if [ -f "$VENV_TAR" ]; then
        echo "Entpacke venv nach lokale SSD ($VENV_LOCAL)..."
        mkdir -p "$VENV_LOCAL"
        tar -xf "$VENV_TAR" --strip-components=1 -C "$VENV_LOCAL"
    else
        echo "WARNUNG: $VENV_TAR nicht gefunden, verwende cp (mehr HOME I/O)."
        echo "         Einmalig ausführen: cd $DELFIN_DIR && tar -cf delfin_venv.tar .venv/"
        cp -a "$VENV_DIR" "$VENV_LOCAL"
    fi
    sed -i "s|$VENV_DIR|$VENV_LOCAL|g" "$VENV_LOCAL/bin/activate" 2>/dev/null || true
    source "$VENV_LOCAL/bin/activate"
    echo "venv von lokaler SSD geladen."
else
    echo "WARNUNG: kein TMPDIR verfuegbar, lade venv direkt von HOME."
    source "$VENV_DIR/bin/activate"
fi

# Python HOME I/O weiter reduzieren
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

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
rm -rf "$DELFIN_SCRATCH" "$ORCA_TMP" "${VENV_LOCAL:-}"

exit $EXIT_CODE
