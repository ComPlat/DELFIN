#!/bin/bash
#SBATCH --job-name=delfin_small
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # Nur 16 Cores (statt 40)
#SBATCH --threads-per-core=1
#SBATCH --mem=64G                    # Weniger RAM
#SBATCH --time=24:00:00              # Kürzere Zeit
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err

# ========================================================================
# DELFIN Small Job - BwUniCluster 3.0
# ========================================================================
#
# Für mehr Parallelität mit weniger Cores pro Job:
# - 16 CPUs, 64 GB RAM, 24 Stunden
# - Ermöglicht mehr gleichzeitige Jobs
# - Gut für kleinere/mittlere Moleküle
#
# Verwendung:
#   Wenn du mehr Jobs parallel laufen lassen willst:
#   10 Jobs × 16 Cores = 160 Cores (statt 5 Jobs × 40 Cores = 200 Cores)
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

# Scratch
export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
export ORCA_TMP="/scratch/${USER}/orca_${SLURM_JOB_ID}"
mkdir -p "$DELFIN_SCRATCH"
mkdir -p "$ORCA_TMP"

# Job Info
echo "========================================"
echo "DELFIN Small Job (16 Cores)"
echo "========================================"
echo "Job ID:      $SLURM_JOB_ID"
echo "Node:        $SLURM_JOB_NODELIST"
echo "CPUs:        $SLURM_CPUS_PER_TASK"
echo "Memory:      $SLURM_MEM_PER_NODE MB"
echo "Working Dir: $(pwd)"
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
