#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem=240G
#SBATCH --time=120:00:00
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --constraint=BEEOND
#SBATCH --exclusive

# ========================================================================
# DELFIN Extra Long Job (120h / 5 days)
# ========================================================================

set -euo pipefail

# Module laden
module purge
module load devel/python/3.11.7-gnu-14.2

# Basisverzeichnis automatisch ermitteln (relativ zum Skript-Standort)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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

# Nutze selbst-installiertes OpenMPI 4.1.8 (kompatibel mit ORCA)
if [ ! -d "$BASE_DIR/openmpi-4.1.8" ]; then
    echo "ERROR: OpenMPI 4.1.8 not found in $BASE_DIR/openmpi-4.1.8"
    echo "Please install it first. See installation instructions."
    exit 1
fi

echo "Using custom OpenMPI 4.1.8 from $BASE_DIR/openmpi-4.1.8"
export PATH="$BASE_DIR/openmpi-4.1.8/bin:$PATH"
export LD_LIBRARY_PATH="$BASE_DIR/openmpi-4.1.8/lib:$LD_LIBRARY_PATH"

# ORCA Pfad setzen
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

# MPI Umgebungsvariablen fuer Cluster
export OMPI_MCA_btl_tcp_if_include=ib0
export OMPI_MCA_btl="^openib"

# DELFIN Environment aktivieren
source "$DELFIN_DIR/.venv/bin/activate"

# Umgebungsvariablen
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export DELFIN_ORCA_PROGRESS=0
export MPLBACKEND=Agg

# Scratch-Verzeichnis
if [ -n "${BEEOND_MOUNTPOINT:-}" ]; then
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID}"
    export ORCA_TMPDIR="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID}"
    echo "Using BeeOND: $BEEOND_MOUNTPOINT"
else
    export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
    export ORCA_TMPDIR="/scratch/${USER}/orca_${SLURM_JOB_ID}"
fi

mkdir -p "$DELFIN_SCRATCH" "$ORCA_TMPDIR"

# Job Info
echo "========================================"
echo "DELFIN Job - {job_name}"
echo "Job Type: Extra Long (120h / 5 days)"
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

# OpenMPI Version pruefen
echo "OpenMPI: $(which mpirun)"
mpirun --version 2>&1 | head -1 || echo "MPI check failed"
echo ""

# ORCA Version pruefen
echo "ORCA: $(which orca)"
orca --version 2>&1 | head -5 || echo "ORCA check failed"
echo ""

# DELFIN Version
cd "$DELFIN_DIR"
echo "DELFIN Version: $(delfin --version 2>&1 || echo 'unknown')"
echo "Git Branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'N/A')"
echo "Git Commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'N/A')"
echo ""
cd - > /dev/null

# Inputs ins Scratch kopieren und dort arbeiten
RUN_DIR="$DELFIN_SCRATCH/run"
mkdir -p "$RUN_DIR"

cp -a "$SLURM_SUBMIT_DIR"/{CONTROL.txt,input.txt} "$RUN_DIR"/

cd "$RUN_DIR"

# DELFIN ausfuehren
delfin

EXIT_CODE=$?

echo ""
echo "========================================"
echo "Job beendet: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "========================================"

# Ergebnisse zurueckkopieren
cp -a "$RUN_DIR"/* "$SLURM_SUBMIT_DIR"/

# Cleanup
rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR"

exit $EXIT_CODE
