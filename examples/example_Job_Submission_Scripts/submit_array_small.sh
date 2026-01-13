#!/bin/bash
#SBATCH --job-name=delfin_array
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # 16 Cores pro Job
#SBATCH --threads-per-core=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-20%5               # 20 Jobs, max 5 gleichzeitig
#SBATCH --output=delfin_%A_%a.out
#SBATCH --error=delfin_%A_%a.err

# ========================================================================
# DELFIN Small Job Array - BwUniCluster 3.0
# ========================================================================
#
# Realistisch für normale Nutzer:
# - 20 Jobs total
# - Max 5 gleichzeitig (--array=1-20%5)
# - 16 Cores pro Job
# - Gesamt: 5 × 16 = 80 Cores gleichzeitig
#
# Struktur:
#   molecule_1/CONTROL.txt + input.txt
#   molecule_2/CONTROL.txt + input.txt
#   ...
#   molecule_20/CONTROL.txt + input.txt
#
# Verwendung:
#   sbatch submit_array_small.sh
#
# Das Array startet automatisch die nächsten Jobs wenn frühere fertig sind!
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

# Molekül-Verzeichnisse
MOLECULES=(
    "none"          # Index 0 (wird nicht genutzt, Array startet bei 1)
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
    "molecule_11"
    "molecule_12"
    "molecule_13"
    "molecule_14"
    "molecule_15"
    "molecule_16"
    "molecule_17"
    "molecule_18"
    "molecule_19"
    "molecule_20"
)

# Aktuelles Molekül
MOLECULE_DIR=${MOLECULES[$SLURM_ARRAY_TASK_ID]}

echo "========================================"
echo "DELFIN Job Array (Small)"
echo "========================================"
echo "Array Job ID:    $SLURM_ARRAY_JOB_ID"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_COUNT"
echo "Molecule:        $MOLECULE_DIR"
echo "Node:            $SLURM_JOB_NODELIST"
echo "CPUs:            $SLURM_CPUS_PER_TASK"
echo "Started:         $(date)"
echo "========================================"
echo ""

# Prüfe Verzeichnis
if [ ! -d "$MOLECULE_DIR" ]; then
    echo "ERROR: Directory $MOLECULE_DIR not found!"
    echo "Erstelle zuerst: mkdir $MOLECULE_DIR && cd $MOLECULE_DIR && delfin --define"
    exit 1
fi

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
