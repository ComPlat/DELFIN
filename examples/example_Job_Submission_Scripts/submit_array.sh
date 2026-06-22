#!/bin/bash
#SBATCH --job-name=delfin_array
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # Weniger Cores pro Job
#SBATCH --threads-per-core=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-50%10              # 50 Jobs, max 10 gleichzeitig
#SBATCH --output=delfin_%A_%a.out
#SBATCH --error=delfin_%A_%a.err

# ========================================================================
# DELFIN Job Array - BwUniCluster 3.0
# ========================================================================
#
# Für High-Throughput Screening:
# - Mehrere Moleküle parallel verarbeiten
# - Automatische Workload-Verteilung
# - Array-Syntax: 1-50%10 = 50 Jobs, max 10 parallel
#
# Verwendung:
#   1. Erstelle Verzeichnisse: molecule_1/, molecule_2/, ...
#   2. Jedes mit CONTROL.txt und input.txt
#   3. sbatch submit_array.sh
#
# ========================================================================

# Module laden - HIER ANPASSEN!
module purge
module load chem/orca/6.1.1
module load devel/python/3.11

# DELFIN Environment
source ~/DELFIN/.venv/bin/activate

# Umgebung
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Array von Molekül-Verzeichnissen
MOLECULES=(
    "molecule_1"
    "molecule_2"
    "molecule_3"
    "molecule_4"
    "molecule_5"
    "molecule_6"
    "molecule_7"
    "molecule_8"
    "molecule_9"
    "molecule_10"
    # ... bis molecule_50
)

# Aktuelles Molekül aus Array-ID
MOLECULE_DIR=${MOLECULES[$SLURM_ARRAY_TASK_ID - 1]}

echo "========================================"
echo "DELFIN Job Array"
echo "========================================"
echo "Array Job ID:    $SLURM_ARRAY_JOB_ID"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Molecule:        $MOLECULE_DIR"
echo "Node:            $SLURM_JOB_NODELIST"
echo "CPUs:            $SLURM_CPUS_PER_TASK"
echo "Started:         $(date)"
echo "========================================"
echo ""

# Prüfe ob Verzeichnis existiert
if [ ! -d "$MOLECULE_DIR" ]; then
    echo "ERROR: Directory $MOLECULE_DIR not found!"
    exit 1
fi

# Wechsle ins Molekül-Verzeichnis
cd "$MOLECULE_DIR"

# Scratch
export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
export ORCA_TMP="/scratch/${USER}/orca_${SLURM_JOB_ID}"
mkdir -p "$DELFIN_SCRATCH"
mkdir -p "$ORCA_TMP"

# DELFIN ausführen
delfin

EXIT_CODE=$?

echo ""
echo "========================================"
echo "Molecule $MOLECULE_DIR finished: $(date)"
echo "Exit Code: $EXIT_CODE"
echo "========================================"

# Cleanup
rm -rf "$DELFIN_SCRATCH"
rm -rf "$ORCA_TMP"

exit $EXIT_CODE
