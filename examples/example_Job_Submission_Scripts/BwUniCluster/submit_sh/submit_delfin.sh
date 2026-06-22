#!/bin/bash
#SBATCH --job-name=delfin_job
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --threads-per-core=1
#SBATCH --mem=240G
#SBATCH --time=48:00:00
#SBATCH --output=delfin_%j.out
#SBATCH --error=delfin_%j.err
#SBATCH --signal=B:SIGTERM@120
# ========================================================================
# BwUniCluster convenience wrapper
# ========================================================================
#
# This script sets BwUniCluster-specific defaults and delegates to the
# central DELFIN SLURM submit script (delfin/submit_templates/submit_delfin.sh).
#
# All site-specific features are controlled via environment variables.
# See the central script for the full documentation.
# ========================================================================

# BwUniCluster defaults — override any of these before calling this script
export DELFIN_MODULES="${DELFIN_MODULES:-devel/python/3.11.7-gnu-14.2}"
export DELFIN_STAGE_ORCA="${DELFIN_STAGE_ORCA:-1}"
export DELFIN_STAGE_VENV="${DELFIN_STAGE_VENV:-1}"
export DELFIN_RUNTIME_CACHE="${DELFIN_RUNTIME_CACHE:-1}"
export DELFIN_AUTO_RESOURCES="${DELFIN_AUTO_RESOURCES:-1}"
export DELFIN_NODE_CORES="${DELFIN_NODE_CORES:-96}"
export DELFIN_NODE_MEM_MB="${DELFIN_NODE_MEM_MB:-$((384 * 1024))}"
export DELFIN_HIGHMEM_MB="${DELFIN_HIGHMEM_MB:-$((2304 * 1024))}"

# Locate the central submit script via the delfin package
_CENTRAL_SCRIPT="${DELFIN_SUBMIT_SCRIPT:-}"
if [ -z "$_CENTRAL_SCRIPT" ]; then
    _CENTRAL_SCRIPT="$(python3 -c '
from delfin.runtime_setup import get_packaged_submit_templates_dir
print(get_packaged_submit_templates_dir() / "submit_delfin.sh")
' 2>/dev/null || true)"
fi

if [ -z "$_CENTRAL_SCRIPT" ] || [ ! -f "$_CENTRAL_SCRIPT" ]; then
    # Fallback: look relative to this script's location
    _SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    _CENTRAL_SCRIPT="$(cd "$_SCRIPT_DIR/../../../../delfin/submit_templates" 2>/dev/null && pwd)/submit_delfin.sh"
fi

if [ ! -f "$_CENTRAL_SCRIPT" ]; then
    echo "ERROR: Cannot locate central DELFIN submit script."
    echo "       Expected at: $_CENTRAL_SCRIPT"
    echo "       Ensure delfin is installed or set DELFIN_SUBMIT_SCRIPT to the central script."
    exit 1
fi

exec bash "$_CENTRAL_SCRIPT" "$@"
