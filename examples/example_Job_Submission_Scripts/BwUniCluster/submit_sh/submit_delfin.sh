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
# Note: Add --constraint=BEEOND only for multi-node jobs that need shared scratch
# For single-node jobs, $TMPDIR (local SSD) is faster and automatically available

# ========================================================================
# DELFIN Central Submit Script for BwUniCluster
# ========================================================================
#
# MODES (set via DELFIN_MODE environment variable):
#   delfin | delfin-recalc | orca | build | auto (default: auto)
#   BUILD_MULTIPLICITY: Spin multiplicity for build mode (default: 1)
#   DELFIN_INP_FILE: Specific .inp file for ORCA mode
#   DELFIN_JOB_NAME: Job name for display (optional)
#
# RESOURCE PARAMETERS (via sbatch command-line overrides):
#   --time, --ntasks, --cpus-per-task, --mem, --job-name
#
# ORCA PARALLELISM NOTE (important):
#   Prefer: --ntasks=1 --cpus-per-task=<PAL>
#   Avoid:  --ntasks=<PAL> --cpus-per-task=1
#   ORCA manages its own parallel workers and is more stable when PAL is
#   expressed as cpus-per-task on a single task allocation.
#
# AUTO RESOURCES (optional):
#   DELFIN_AUTO_RESOURCES=1 ./submit_delfin.sh
#   Uses CONTROL.txt PAL/maxcore to set --ntasks/--mem and adds --exclusive
#   if the request is close to a full standard node.
#   DELFIN_FORCE_EXCLUSIVE=1 forces --exclusive even for smaller requests.
#
# SCRATCH DIRECTORY (automatic selection):
#   1) BeeOND (if --constraint=BEEOND): Parallel SSD across all nodes
#      - Best for multi-node jobs needing shared scratch
#      - Add: sbatch --constraint=BEEOND submit_delfin.sh
#   2) $TMPDIR (default): Local SSD on compute node
#      - Fastest option for single-node jobs
#      - Automatically cleaned after job
#   3) /scratch (fallback): Network filesystem
#      - Slower, only used if TMPDIR unavailable
#
# ========================================================================

set -euo pipefail

# === Optional auto-submit with resources derived from CONTROL.txt ===
if [ -z "${SLURM_JOB_ID:-}" ] && [ "${DELFIN_AUTO_RESOURCES:-0}" = "1" ]; then
    control_file="${DELFIN_CONTROL:-$PWD/CONTROL.txt}"
    pal_default=40
    maxcore_default=6000
    pal="$pal_default"
    maxcore="$maxcore_default"
    if [ -f "$control_file" ]; then
        pal_found="$(awk -F= '/^[[:space:]]*PAL[[:space:]]*=/ {gsub(/[^0-9]/,"",$2); print $2; exit}' "$control_file")"
        maxcore_found="$(awk -F= '/^[[:space:]]*maxcore[[:space:]]*=/ {gsub(/[^0-9]/,"",$2); print $2; exit}' "$control_file")"
        if [ -n "$pal_found" ]; then
            pal="$pal_found"
        fi
        if [ -n "$maxcore_found" ]; then
            maxcore="$maxcore_found"
        fi
    fi

    # ORCA is typically most stable with a single task and PAL expressed
    # as cpus-per-task.
    ntasks="1"
    cpus_per_task="$pal"
    mem_mb=$((pal * maxcore))

    node_cores=96
    node_mem_mb=$((384 * 1024))
    highmem_mb=$((2304 * 1024))
    partition="${DELFIN_PARTITION:-cpu}"
    exclusive=0
    if [ "$maxcore" -ge 9500 ] && [ "$partition" = "cpu" ]; then
        partition="highmem"
        node_mem_mb="$highmem_mb"
    elif [ "$mem_mb" -gt "$node_mem_mb" ] && [ "$partition" = "cpu" ]; then
        partition="highmem"
        node_mem_mb="$highmem_mb"
    fi
    if [ "$cpus_per_task" -ge $((node_cores * 9 / 10)) ] || [ "$mem_mb" -ge $((node_mem_mb * 9 / 10)) ]; then
        exclusive=1
    fi
    if [ "${DELFIN_FORCE_EXCLUSIVE:-0}" = "1" ]; then
        exclusive=1
    fi
    if [ "$mem_mb" -gt "$node_mem_mb" ]; then
        echo "ERROR: Requested mem ${mem_mb}MB exceeds partition capacity (${node_mem_mb}MB)."
        echo "Set DELFIN_PARTITION=highmem (or adjust PAL/maxcore) and submit again."
        exit 1
    fi

    sbatch_args=(
        --partition="${partition}"
        --nodes=1
        --ntasks="${ntasks}"
        --cpus-per-task="${cpus_per_task}"
        --threads-per-core=1
        --mem="${mem_mb}M"
    )
    if [ "$exclusive" -eq 1 ]; then
        sbatch_args+=(--exclusive)
    fi
    echo "Auto resources: PAL=${pal} maxcore=${maxcore} -> --ntasks=${ntasks} --cpus-per-task=${cpus_per_task} --mem=${mem_mb}M --partition=${partition}"
    if [ "$exclusive" -eq 1 ]; then
        echo "Auto resources: requesting --exclusive (near full node)"
    fi
    sbatch "${sbatch_args[@]}" "$0" "$@"
    exit 0
fi

# === Configuration from environment variables ===
MODE="${DELFIN_MODE:-auto}"
INP_FILE="${DELFIN_INP_FILE:-}"
OVERRIDE="${DELFIN_OVERRIDE:-}"
DISPLAY_JOB_NAME="${DELFIN_JOB_NAME:-${SLURM_JOB_NAME:-delfin_job}}"

# === PAL / SLURM sanity check (helps diagnose ORCA MPI/LINALG crashes) ===
PAL_FROM_CONTROL=""
MAXCORE_FROM_CONTROL=""
CONTROL_PATH="${SLURM_SUBMIT_DIR:-$PWD}/CONTROL.txt"
if [ -f "$CONTROL_PATH" ]; then
    PAL_FROM_CONTROL="$(awk -F= '/^[[:space:]]*PAL[[:space:]]*=/ {gsub(/[^0-9]/,"",$2); print $2; exit}' "$CONTROL_PATH")"
    MAXCORE_FROM_CONTROL="$(awk -F= '/^[[:space:]]*maxcore[[:space:]]*=/ {gsub(/[^0-9]/,"",$2); print $2; exit}' "$CONTROL_PATH")"
fi

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SLURM_TASKS="${SLURM_NTASKS:-1}"
    SLURM_CPT="${SLURM_CPUS_PER_TASK:-1}"
    if [ -n "$PAL_FROM_CONTROL" ]; then
        if [ "$SLURM_TASKS" -gt 1 ] && [ "$SLURM_CPT" -eq 1 ]; then
            echo "WARNING: SLURM allocation is --ntasks=${SLURM_TASKS} --cpus-per-task=${SLURM_CPT}."
            echo "         For ORCA stability, prefer --ntasks=1 --cpus-per-task=${PAL_FROM_CONTROL}."
            echo "         Mismatched task-style allocations can trigger ORCA linear algebra crashes"
            echo "         (e.g., BLAS_CholeskySolution: Wrong RHS dimension)."
        fi
        if [ "$SLURM_TASKS" -ne 1 ] && [ "$SLURM_CPT" -ne "$PAL_FROM_CONTROL" ]; then
            echo "WARNING: CONTROL.txt PAL=${PAL_FROM_CONTROL} but SLURM has ntasks=${SLURM_TASKS}, cpus-per-task=${SLURM_CPT}."
            echo "         Suggested submit command:"
            if [ -n "$MAXCORE_FROM_CONTROL" ]; then
                echo "         sbatch --ntasks=1 --cpus-per-task=${PAL_FROM_CONTROL} --mem=$((PAL_FROM_CONTROL * MAXCORE_FROM_CONTROL))M \"$0\""
            else
                echo "         sbatch --ntasks=1 --cpus-per-task=${PAL_FROM_CONTROL} \"$0\""
            fi
        fi
    fi
fi

# Load modules
module purge
module load devel/python/3.11.7-gnu-14.2

# Initialize LD_LIBRARY_PATH if not set
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# Auto-detect base directory (search upwards for software/delfin)
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

# Use custom OpenMPI 4.1.8 (compatible with ORCA)
if [ ! -d "$BASE_DIR/openmpi-4.1.8" ]; then
    echo "ERROR: OpenMPI 4.1.8 not found in $BASE_DIR/openmpi-4.1.8"
    echo "Please install it first. See installation instructions."
    exit 1
fi

echo "Using custom OpenMPI 4.1.8 from $BASE_DIR/openmpi-4.1.8"
export PATH="$BASE_DIR/openmpi-4.1.8/bin:$PATH"
export LD_LIBRARY_PATH="$BASE_DIR/openmpi-4.1.8/lib:$LD_LIBRARY_PATH"

# Set ORCA path (allow override from environment)
ORCA_BASE_DEFAULT="$BASE_DIR/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2"
ORCA_BASE="${DELFIN_ORCA_BASE:-$ORCA_BASE_DEFAULT}"
if [ ! -d "$ORCA_BASE" ]; then
    echo "ERROR: ORCA not found in $ORCA_BASE"
    echo "Please install ORCA or set DELFIN_ORCA_BASE."
    exit 1
fi
export PATH="$ORCA_BASE:$PATH"
export LD_LIBRARY_PATH="$ORCA_BASE:$LD_LIBRARY_PATH"
export ORCA_DIR="$ORCA_BASE"
export ORCA_PLOT="$ORCA_BASE/orca_plot"

# =============================================================================
# MPI Configuration for ORCA on BwUniCluster
# Community-recommended stable settings (ORCA Input Library)
# =============================================================================

# Use ob1 PML with vader (shared memory) - sm was removed in OpenMPI 3.0
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp,vader

# OpenMPI settings optimized for ORCA stability
export OMPI_MCA_mpi_show_mca_params_file=0
export OMPI_MCA_mpi_yield_when_idle=1         # Reduce CPU spinning
export OMPI_MCA_coll_hcoll_enable=0           # Disable HCOLL (often problematic)

# Clear any inherited TCP interface include/exclude to avoid invalid if_inexclude warnings
unset OMPI_MCA_btl_tcp_if_include
unset OMPI_MCA_btl_tcp_if_exclude

# Process binding: let ORCA handle its own binding
export OMPI_MCA_hwloc_base_binding_policy=none
export OMPI_MCA_rmaps_base_mapping_policy=core
export OMPI_MCA_rmaps_base_oversubscribe=true # Allow ORCA's dynamic parallelism

# Activate DELFIN virtual environment
source "$DELFIN_DIR/.venv/bin/activate"

# Environment variables
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export DELFIN_ORCA_PROGRESS=0
export MPLBACKEND=Agg

# Scratch directory setup
# Priority: 1) BeeOND (for multi-node), 2) $TMPDIR (local SSD), 3) /scratch (network, last resort)
# Normalize TMPDIR to SLURM_TMPDIR when available
if [ -z "${TMPDIR:-}" ] && [ -n "${SLURM_TMPDIR:-}" ]; then
    export TMPDIR="$SLURM_TMPDIR"
fi
if [ -n "${BEEOND_MOUNTPOINT:-}" ]; then
    # BeeOND: parallel filesystem across all job nodes (best for multi-node)
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID}"
    export ORCA_TMPDIR="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID}"
    echo "Using BeeOND (parallel SSD): $BEEOND_MOUNTPOINT (job-isolated: ${SLURM_JOB_ID})"
elif [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
    # $TMPDIR: local SSD on compute node (fastest for single-node jobs)
    # Note: Only visible on local node, auto-cleaned after job
    # Use SLURM_JOB_ID to guarantee isolation even if TMPDIR is shared
    export DELFIN_SCRATCH="$TMPDIR/delfin_${SLURM_JOB_ID}"
    export ORCA_TMPDIR="$TMPDIR/orca_${SLURM_JOB_ID}"
    echo "Using local SSD (TMPDIR): $TMPDIR (job-isolated: ${SLURM_JOB_ID})"
else
    # Fallback to /scratch (network filesystem - slower, avoid if possible)
    export DELFIN_SCRATCH="/scratch/${USER}/delfin_${SLURM_JOB_ID}"
    export ORCA_TMPDIR="/scratch/${USER}/orca_${SLURM_JOB_ID}"
    echo "WARNING: Using /scratch (network filesystem) - consider using TMPDIR or BeeOND for better I/O performance"
fi

# ORCA honors both ORCA_TMPDIR and ORCA_SCRDIR; set both for consistency.
export ORCA_SCRDIR="$ORCA_TMPDIR"

# Ensure delfin uses a unique run token for ORCA scratch isolation
export DELFIN_RUN_TOKEN="slurm_${SLURM_JOB_ID}"

RUN_DIR="$DELFIN_SCRATCH/run"
mkdir -p "$RUN_DIR" "$ORCA_TMPDIR"

# Write SLURM reason/timelimit info into run dir for debugging.
write_slurm_reason() {
    local out_file="$RUN_DIR/slurm_reason.txt"
    {
        echo "Timestamp: $(date)"
        echo "Job ID: ${SLURM_JOB_ID:-unknown}"
        echo "Node: ${SLURM_JOB_NODELIST:-unknown}"
        if command -v scontrol >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
            scontrol show job "$SLURM_JOB_ID" 2>/dev/null | tr ' ' '\n' | egrep 'JobId=|JobName=|State=|Reason=|TimeLimit=|RunTime=|TimeMin=|StartTime=|EndTime=' || true
        fi
        if command -v sacct >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
            sacct -j "$SLURM_JOB_ID" --format=JobID,JobName%30,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS,NodeList,Reason -P 2>/dev/null || true
        fi
    } > "$out_file" 2>/dev/null || true
}

write_slurm_reason

# Cleanup function for trap (handles SIGTERM from timeout, SIGINT, etc.)
cleanup() {
    local signal_name="${1:-UNKNOWN}"
    echo ""
    echo "========================================"
    echo "Caught $signal_name signal, running cleanup..."
    echo "========================================"

    # Try DELFIN cleanup first (safe if not running DELFIN)
    cd "$RUN_DIR" 2>/dev/null && delfin --cleanup 2>/dev/null || true

    # Capture scheduler reason at shutdown.
    write_slurm_reason

    # CRITICAL: Copy ALL results back before cleanup
    echo "Copying results back to $SLURM_SUBMIT_DIR..."
    if [ -d "$RUN_DIR" ]; then
        rsync -a --exclude='*.tmp' --exclude='.orca_iso*' "$RUN_DIR"/ "$SLURM_SUBMIT_DIR"/ 2>/dev/null || true
        echo "Results copied successfully."
    else
        echo "WARNING: RUN_DIR not found, nothing to copy."
    fi

    # Cleanup scratch (only after successful copy)
    rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR" 2>/dev/null || true
    echo "Cleanup completed at $(date)"
    exit 1
}

# Periodic copy-back to mitigate hard timeouts.
periodic_copy() {
    while true; do
        sleep 7200
        if [ -d "$RUN_DIR" ]; then
            rsync -a --exclude='*.tmp' --exclude='.orca_iso*' "$RUN_DIR"/ "$SLURM_SUBMIT_DIR"/ 2>/dev/null || true
        fi
    done
}

# Start periodic copy-back in background.
periodic_copy &
PERIODIC_COPY_PID=$!

# Get walltime in seconds from SLURM
get_walltime_seconds() {
    local timelimit
    timelimit=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | grep -oP 'TimeLimit=\K[^ ]+' || echo "")
    if [ -z "$timelimit" ] || [ "$timelimit" = "UNLIMITED" ]; then
        echo "0"
        return
    fi
    # Format: D-HH:MM:SS or HH:MM:SS or MM:SS
    echo "$timelimit" | awk -F'[-:]' '{
        if (NF==4) print $1*86400 + $2*3600 + $3*60 + $4
        else if (NF==3) print $1*3600 + $2*60 + $3
        else if (NF==2) print $1*60 + $2
        else print 0
    }'
}

# Schedule final backup 5 minutes before timeout
SAFETY_MARGIN=300  # 5 minutes in seconds
schedule_final_backup() {
    local walltime_sec
    walltime_sec=$(get_walltime_seconds)
    if [ "$walltime_sec" -le "$SAFETY_MARGIN" ]; then
        echo "Walltime too short for scheduled final backup, skipping."
        return
    fi
    local wait_time=$((walltime_sec - SAFETY_MARGIN))
    echo "Scheduled final backup in $((wait_time / 3600))h $((wait_time % 3600 / 60))m (5 min before timeout)"
    sleep "$wait_time"
    echo ""
    echo "========================================"
    echo "Final backup 5 min before timeout: $(date)"
    echo "========================================"
    if [ -d "$RUN_DIR" ]; then
        rsync -a --exclude='*.tmp' --exclude='.orca_iso*' "$RUN_DIR"/ "$SLURM_SUBMIT_DIR"/ 2>/dev/null || true
        echo "Final backup completed."
    fi
}

# Start scheduled final backup in background
schedule_final_backup &
FINAL_BACKUP_PID=$!

stop_periodic_copy() {
    if [ -n "${PERIODIC_COPY_PID:-}" ]; then
        kill "$PERIODIC_COPY_PID" 2>/dev/null || true
    fi
    if [ -n "${FINAL_BACKUP_PID:-}" ]; then
        kill "$FINAL_BACKUP_PID" 2>/dev/null || true
    fi
}

# Set trap for cleanup on termination signals (including SLURM timeout)
trap 'stop_periodic_copy' EXIT
trap 'cleanup SIGTERM' SIGTERM
trap 'cleanup SIGINT' SIGINT
trap 'cleanup SIGHUP' SIGHUP

# Job info
echo "========================================"
echo "DELFIN Job - $DISPLAY_JOB_NAME"
echo "Mode: $MODE"
echo "========================================"
echo "Job ID:      $SLURM_JOB_ID"
echo "Node:        $SLURM_JOB_NODELIST"
echo "CPUs:        $SLURM_NTASKS"
echo "Memory:      ${SLURM_MEM_PER_NODE:-unknown} MB"
# SLURM_TIMELIMIT is not always set; fall back to scontrol for visibility.
TIME_LIMIT="${SLURM_TIMELIMIT:-}"
if [ -z "$TIME_LIMIT" ] && [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
    TIME_LIMIT="$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | awk -F= '/TimeLimit=/{print $2; exit}')"
fi
echo "Time Limit:  ${TIME_LIMIT:-unknown}"
echo "Submit Dir:  $SLURM_SUBMIT_DIR"
echo "Scratch:     $DELFIN_SCRATCH"
echo "Started:     $(date)"
if [ -n "$INP_FILE" ]; then
    echo "Input File:  $INP_FILE"
fi
echo "========================================"
echo ""

# Check OpenMPI version
echo "OpenMPI: $(which mpirun)"
mpirun --version 2>&1 | head -1 || echo "MPI check failed"
echo ""

# Check ORCA version
orca_path=$(command -v orca || true)
echo "ORCA: ${orca_path:-not found}"
if [ -n "$orca_path" ]; then
    orca --version 2>&1 | head -5 || true
fi
echo ""

# Check DELFIN version
cd "$DELFIN_DIR"
echo "DELFIN Version: $(delfin --version 2>&1 || echo 'unknown')"
echo "Git Branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'N/A')"
echo "Git Commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'N/A')"
echo ""
cd - > /dev/null

# Copy ALL input files to scratch
echo "Copying input files to scratch..."
cp -a "$SLURM_SUBMIT_DIR"/* "$RUN_DIR"/ 2>/dev/null || true
# Remove output files from previous runs (if any)
rm -f "$RUN_DIR"/delfin_*.out "$RUN_DIR"/delfin_*.err 2>/dev/null || true

cd "$RUN_DIR"

# Auto-detect mode if set to "auto"
if [ "$MODE" = "auto" ]; then
    if [ -f "$SLURM_SUBMIT_DIR/CONTROL.txt" ] && [ -f "$SLURM_SUBMIT_DIR/input.txt" ]; then
        MODE="delfin"
        echo "Auto-detected mode: DELFIN (CONTROL.txt + input.txt found)"
    elif ls "$SLURM_SUBMIT_DIR"/*.inp 1>/dev/null 2>&1; then
        MODE="orca"
        echo "Auto-detected mode: ORCA (*.inp files found, no CONTROL.txt)"
    else
        echo "ERROR: No valid input files found in $SLURM_SUBMIT_DIR"
        echo "       Expected either: CONTROL.txt + input.txt (DELFIN mode)"
        echo "       Or: *.inp files (ORCA-only mode)"
        exit 1
    fi
fi

# Run appropriate mode
EXIT_CODE=0
set +e
case "$MODE" in
    delfin)
        echo "Starting DELFIN..."
        delfin
        EXIT_CODE=$?
        ;;
    delfin-recalc)
        echo "Starting DELFIN --recalc..."
        delfin --recalc
        EXIT_CODE=$?
        ;;
    delfin-recalc-override)
        if [ -z "$OVERRIDE" ]; then
            echo "ERROR: DELFIN_OVERRIDE not set for delfin-recalc-override mode"
            EXIT_CODE=1
        else
            echo "Starting DELFIN --recalc --occupier-override $OVERRIDE..."
            delfin --recalc --occupier-override "$OVERRIDE"
            EXIT_CODE=$?
        fi
        ;;
    orca)
        # Use specified input file or find the first one
        if [ -z "$INP_FILE" ]; then
            INP_FILE=$(ls *.inp 2>/dev/null | head -1)
        fi
        if [ -z "$INP_FILE" ]; then
            echo "ERROR: No .inp file found for ORCA mode"
            EXIT_CODE=1
        elif [ ! -f "$INP_FILE" ]; then
            echo "ERROR: Specified input file not found: $INP_FILE"
            EXIT_CODE=1
        else
            OUT_FILE="${INP_FILE%.inp}.out"
            echo "Starting ORCA: $INP_FILE -> $OUT_FILE"
            "$ORCA_BASE/orca" "$INP_FILE" > "$OUT_FILE" 2>&1
            EXIT_CODE=$?
        fi
        ;;
    build)
        # Build up metal complex step by step using XTB DOCKER
        BUILD_MULT="${BUILD_MULTIPLICITY:-1}"
        echo "Starting delfin-build (complex build-up)..."
        echo "  Multiplicity: $BUILD_MULT"
        python -m delfin.build_up_complex input.txt --goat --directory builder --multiplicity "$BUILD_MULT" --verbose
        EXIT_CODE=$?
        ;;
    *)
        echo "ERROR: Unknown mode: $MODE"
        echo "       Valid modes: delfin, delfin-recalc, delfin-recalc-override, orca, build, auto"
        EXIT_CODE=1
        ;;
esac
set -e

echo ""
echo "========================================"
echo "Job finished: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "========================================"

# Copy results back
rsync -a --exclude='*.tmp' --exclude='.orca_iso*' "$RUN_DIR"/ "$SLURM_SUBMIT_DIR"/

# Cleanup scratch
rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR"

exit $EXIT_CODE
