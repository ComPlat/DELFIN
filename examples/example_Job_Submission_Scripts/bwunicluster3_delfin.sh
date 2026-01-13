#!/bin/bash
#SBATCH --job-name=delfin
#SBATCH --partition=cpu              # Hauptpartition (oder dev_cpu für Tests)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40           # 40 CPUs pro Node (oder 28 je nach Node-Typ)
#SBATCH --threads-per-core=1         # Hyperthreading deaktivieren (wichtig!)
#SBATCH --mem=120G                   # 120GB Speicher
#SBATCH --time=48:00:00              # Max. Laufzeit
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --constraint=BEEOND          # Optional: Schnelles lokales Storage
#SBATCH --exclusive                  # Exklusiver Node-Zugriff (wenn BEEOND genutzt)

# ========================================================================
# DELFIN auf BwUniCluster 3.0
# ========================================================================
#
# Optimiert für BwUniCluster 3.0 Hardware und SLURM-Konfiguration
#
# Hardware:
# - 40 CPU Cores pro Node (einige Nodes haben 28)
# - Hyperthreading verfügbar aber hier deaktiviert
# - BeeOND für schnelles lokales Storage
#
# ========================================================================

# Module laden (Namen mit 'module avail' prüfen)
module purge
module load chem/orca/6.1.1         # ORCA Version anpassen
module load devel/python/3.11       # Python Version anpassen

# DELFIN Environment aktivieren
source ~/delfin_env/bin/activate

# Umgebungsvariablen setzen
export OMP_NUM_THREADS=1            # DELFIN managed Parallelisierung
export MKL_NUM_THREADS=1            # Intel MKL Threading

# Scratch-Verzeichnis
# Option 1: BEEOND (schnell, wenn --constraint=BEEOND gesetzt)
if [ -n "$BEEOND_MOUNTPOINT" ]; then
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID}"
    echo "Using BeeOND storage: $BEEOND_MOUNTPOINT"
# Option 2: Standard Scratch
else
    export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
    export ORCA_TMP="/scratch/${USER}/orca_${SLURM_JOB_ID}"
fi

mkdir -p "$DELFIN_SCRATCH"
mkdir -p "$ORCA_TMP"

# Job Information ausgeben
echo "========================================="
echo "DELFIN auf BwUniCluster 3.0"
echo "========================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Node:            $SLURM_JOB_NODELIST"
echo "Partition:       $SLURM_JOB_PARTITION"
echo "CPUs:            $SLURM_CPUS_PER_TASK"
echo "Memory:          $SLURM_MEM_PER_NODE MB"
echo "Working Dir:     $(pwd)"
echo "Scratch:         $DELFIN_SCRATCH"
echo "ORCA Scratch:    $ORCA_TMP"
echo "Started:         $(date)"
echo "========================================="
echo ""

# DELFIN ausführen
# DELFIN erkennt automatisch SLURM_CPUS_PER_TASK und nutzt alle verfügbaren Cores
delfin

# Exit-Code speichern
EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job beendet:     $(date)"
echo "Exit Code:       $EXIT_CODE"
echo "========================================="

# Scratch aufräumen
rm -rf "$DELFIN_SCRATCH"
rm -rf "$ORCA_TMP"

exit $EXIT_CODE
