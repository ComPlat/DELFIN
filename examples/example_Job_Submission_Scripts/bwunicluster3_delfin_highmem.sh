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
# venv auf lokale SSD kopieren um HOME I/O zu minimieren
VENV_DIR="$HOME/delfin_env"
VENV_TAR="$HOME/delfin_venv.tar"
if [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
    VENV_LOCAL="$TMPDIR/delfin_venv_${SLURM_JOB_ID}"
    if [ -f "$VENV_TAR" ]; then
        echo "Entpacke venv nach lokale SSD ($VENV_LOCAL)..."
        mkdir -p "$VENV_LOCAL"
        tar -xf "$VENV_TAR" --strip-components=1 -C "$VENV_LOCAL"
    else
        echo "WARNUNG: $VENV_TAR nicht gefunden, verwende cp (mehr HOME I/O)."
        echo "         Einmalig ausführen: cd $HOME && tar -cf delfin_venv.tar delfin_env/"
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
rm -rf "$DELFIN_SCRATCH" "$ORCA_TMP" "${VENV_LOCAL:-}"

exit $EXIT_CODE
