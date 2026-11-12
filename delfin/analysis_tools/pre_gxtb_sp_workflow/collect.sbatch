#!/usr/bin/env bash
#SBATCH --job-name=pre_SP_report
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00
#SBATCH --output=logs/collect-%j.out
#SBATCH --error=logs/collect-%j.err
set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Missing SLURM_SUBMIT_DIR}"
python3 "$SLURM_SUBMIT_DIR/prescreen.py" report --run-dir "$SLURM_SUBMIT_DIR"
