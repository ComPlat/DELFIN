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
#   delfin | delfin-recalc | delfin-recalc-classic | orca | build | guppy | hyperpol_xtb | tadf_xtb | delfin-co2-chain | auto (default: auto)
#   BUILD_MULTIPLICITY: Spin multiplicity for build mode (default: 1)
#   GUPPY_RUNS: Number of GUPPY sampling runs (default: 20)
#   GUPPY_PAL: Total PAL budget for GUPPY (default: DELFIN_PAL or SLURM_CPUS_PER_TASK)
#   GUPPY_MAXCORE: Maxcore per core in MB for GUPPY scheduler (default: DELFIN_MAXCORE or 6000)
#   GUPPY_PARALLEL_JOBS: Number of parallel GUPPY runs sharing resources (default: 4)
#   GUPPY_GOAT_TOPK: Number of top-ranked GUPPY candidates refined with GOAT (default: 3)
#   GUPPY_GOAT_PARALLEL_JOBS: Number of parallel GOAT refinement jobs (default: GUPPY_PARALLEL_JOBS)
#   DELFIN_CO2_SPECIES_DELTA: Redox species delta for delfin-co2-chain mode (default: 0)
#   DELFIN_INP_FILE: Specific .inp file for ORCA mode
#   DELFIN_JOB_NAME: Job name for display (optional)
#   DELFIN_XYZ_FILE: Specific .xyz file for browser-launched hyperpol_xtb/tadf_xtb modes
#   DELFIN_WORKFLOW_LABEL: Workflow label for browser-launched hyperpol_xtb/tadf_xtb modes
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

# Normalise SLURM scratch hints before touching the Python runtime.
if [ -z "${TMPDIR:-}" ] && [ -n "${SLURM_TMPDIR:-}" ]; then
    export TMPDIR="$SLURM_TMPDIR"
fi
STAGE_BASE="${TMPDIR:-${BEEOND_MOUNTPOINT:-}}"

# Activate DELFIN virtual environment from a single tarball to keep HOME I/O low.
VENV_TAR="$DELFIN_DIR/delfin_venv.tar"
VENV_LOCAL=""
if [ -n "${STAGE_BASE:-}" ] && [ -d "${STAGE_BASE}" ]; then
    VENV_LOCAL="$STAGE_BASE/delfin_venv_${SLURM_JOB_ID}"
    if [ ! -f "$VENV_TAR" ]; then
        echo "ERROR: Required runtime tarball not found: $VENV_TAR"
        echo "       Create it once after install: cd $DELFIN_DIR && tar -cf delfin_venv.tar .venv/"
        exit 1
    fi
    echo "Unpacking venv tarball to local SSD ($VENV_LOCAL) to minimise HOME I/O..."
    mkdir -p "$VENV_LOCAL"
    tar -xf "$VENV_TAR" --strip-components=1 -C "$VENV_LOCAL"
    # Rewrite shebangs/paths: activate script uses absolute paths.
    sed -i "s|$DELFIN_DIR/.venv|$VENV_LOCAL|g" "$VENV_LOCAL/bin/activate" 2>/dev/null || true
    source "$VENV_LOCAL/bin/activate"
    echo "venv loaded from local SSD."
else
    echo "WARNING: No local stage base detected; falling back to repository venv."
    echo "         Runtime overlay cache will be skipped in this job."
    source "$DELFIN_DIR/.venv/bin/activate"
fi

# Environment variables
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export DELFIN_ORCA_PROGRESS=0
export MPLBACKEND=Agg
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PIP_DISABLE_PIP_VERSION_CHECK=1
export PIP_NO_INPUT=1

if [ -n "${STAGE_BASE:-}" ] && [ -d "${STAGE_BASE}" ]; then
    export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$STAGE_BASE/.cache_${SLURM_JOB_ID}}"
    export MPLCONFIGDIR="${MPLCONFIGDIR:-$STAGE_BASE/.mplconfig_${SLURM_JOB_ID}}"
    export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$STAGE_BASE/.pip_${SLURM_JOB_ID}}"
    mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$PIP_CACHE_DIR"
fi

RUNTIME_CACHE_DIR="${DELFIN_RUNTIME_CACHE_DIR:-$DELFIN_DIR/.runtime_cache}"
RUNTIME_KEY=""
RUNTIME_BUILD_MODE=""
RUNTIME_WHEEL=""

compute_runtime_dirty_hash() {
    {
        git -C "$DELFIN_DIR" diff --binary HEAD -- delfin pyproject.toml README.md 2>/dev/null || true
        while IFS= read -r rel_path; do
            [ -f "$DELFIN_DIR/$rel_path" ] || continue
            printf '=== %s ===\n' "$rel_path"
            cat "$DELFIN_DIR/$rel_path"
            printf '\n'
        done < <(git -C "$DELFIN_DIR" ls-files --others --exclude-standard -- delfin pyproject.toml README.md 2>/dev/null)
    } | sha256sum | awk '{print $1}'
}

detect_runtime_key() {
    local head=""
    local tree_state=""
    local dirty_hash=""

    if [ -d "$DELFIN_DIR/.git" ] && command -v git >/dev/null 2>&1; then
        head="$(git -C "$DELFIN_DIR" rev-parse --verify HEAD 2>/dev/null || true)"
        if [ -n "$head" ]; then
            tree_state="$(git -C "$DELFIN_DIR" status --porcelain=v1 --untracked-files=normal -- delfin pyproject.toml README.md 2>/dev/null || true)"
            if [ -z "$tree_state" ]; then
                RUNTIME_KEY="git-${head}"
                RUNTIME_BUILD_MODE="git-archive"
                return 0
            fi
            dirty_hash="$(compute_runtime_dirty_hash)"
            if [ -n "$dirty_hash" ]; then
                RUNTIME_KEY="worktree-${head}-${dirty_hash}"
                RUNTIME_BUILD_MODE="live-repo"
                return 0
            fi
        fi
    fi

    RUNTIME_KEY="tree-default"
    RUNTIME_BUILD_MODE="live-repo"
}

build_runtime_context() {
    local build_dir="$1"
    mkdir -p "$build_dir"

    if [ "$RUNTIME_BUILD_MODE" = "git-archive" ]; then
        git -C "$DELFIN_DIR" archive --format=tar HEAD delfin pyproject.toml README.md | tar -xf - -C "$build_dir"
        return 0
    fi

    mkdir -p "$build_dir/delfin"
    cp -a "$DELFIN_DIR/delfin/." "$build_dir/delfin/"
    cp -a "$DELFIN_DIR/pyproject.toml" "$build_dir/"
    if [ -f "$DELFIN_DIR/README.md" ]; then
        cp -a "$DELFIN_DIR/README.md" "$build_dir/"
    fi
}

ensure_runtime_wheel() {
    local runtime_root=""
    local lock_fd_opened=0
    local stage_root=""
    local tmp_build_dir=""
    local tmp_wheel_dir=""
    local built_wheel=""

    if [ -n "$RUNTIME_WHEEL" ] && [ -f "$RUNTIME_WHEEL" ]; then
        return 0
    fi

    detect_runtime_key
    runtime_root="$RUNTIME_CACHE_DIR/$RUNTIME_KEY"
    mkdir -p "$runtime_root"

    if command -v flock >/dev/null 2>&1; then
        exec 9>"$runtime_root/.lock"
        flock 9
        lock_fd_opened=1
    fi

    RUNTIME_WHEEL="$(find "$runtime_root" -maxdepth 1 -type f -name 'delfin_complat-*.whl' | sort | tail -n 1)"
    if [ -n "$RUNTIME_WHEEL" ] && [ -f "$RUNTIME_WHEEL" ]; then
        if [ "$lock_fd_opened" -eq 1 ]; then
            flock -u 9
            exec 9>&-
        fi
        return 0
    fi

    echo "Building DELFIN runtime wheel cache for $RUNTIME_KEY..."
    stage_root="${STAGE_BASE:-$runtime_root}"
    tmp_build_dir="$stage_root/delfin_runtime_build_${SLURM_JOB_ID}_$$"
    tmp_wheel_dir="$stage_root/delfin_runtime_dist_${SLURM_JOB_ID}_$$"
    rm -rf "$tmp_build_dir" "$tmp_wheel_dir"
    mkdir -p "$tmp_build_dir" "$tmp_wheel_dir"

    build_runtime_context "$tmp_build_dir"
    python -m pip wheel \
        --quiet \
        --no-build-isolation \
        --no-deps \
        --wheel-dir "$tmp_wheel_dir" \
        "$tmp_build_dir"

    built_wheel="$(find "$tmp_wheel_dir" -maxdepth 1 -type f -name 'delfin_complat-*.whl' | sort | tail -n 1)"
    if [ -z "$built_wheel" ] || [ ! -f "$built_wheel" ]; then
        echo "ERROR: Failed to build runtime wheel for $RUNTIME_KEY"
        rm -rf "$tmp_build_dir" "$tmp_wheel_dir"
        if [ "$lock_fd_opened" -eq 1 ]; then
            flock -u 9
            exec 9>&-
        fi
        exit 1
    fi

    mv "$built_wheel" "$runtime_root/"
    printf 'runtime_key=%s\nbuild_mode=%s\ncreated=%s\n' \
        "$RUNTIME_KEY" "$RUNTIME_BUILD_MODE" "$(date -Is)" > "$runtime_root/runtime.meta"
    rm -rf "$tmp_build_dir" "$tmp_wheel_dir"

    RUNTIME_WHEEL="$(find "$runtime_root" -maxdepth 1 -type f -name 'delfin_complat-*.whl' | sort | tail -n 1)"
    if [ "$lock_fd_opened" -eq 1 ]; then
        flock -u 9
        exec 9>&-
    fi
}

install_cached_runtime_wheel() {
    if [ -z "${VENV_LOCAL:-}" ] || [ ! -x "$VENV_LOCAL/bin/python" ]; then
        export DELFIN_RUNTIME_KEY="editable-home-fallback"
        echo "WARNING: Skipping runtime wheel overlay because no local venv is staged."
        return 0
    fi

    ensure_runtime_wheel
    echo "Installing cached runtime wheel into local venv ($RUNTIME_KEY)..."
    "$VENV_LOCAL/bin/python" -m pip install \
        --quiet \
        --disable-pip-version-check \
        --no-deps \
        --force-reinstall \
        "$RUNTIME_WHEEL"
    export DELFIN_RUNTIME_KEY="$RUNTIME_KEY"
    echo "Using DELFIN runtime cache: $DELFIN_RUNTIME_KEY"
}

install_cached_runtime_wheel

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

# Auto-detect mode before copying files so fresh runs can avoid mirroring the
# complete HOME workspace into scratch.
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

copy_workspace_to_scratch() {
    local copy_profile="${1:-full}"
    local manifest_file="$DELFIN_SCRATCH/.fresh_inputs.txt"
    local staged_count=0
    echo "Copying files to scratch (profile: ${copy_profile})..."

    if [ "$copy_profile" = "fresh" ]; then
        : > "$manifest_file"
        if [ -f "$SLURM_SUBMIT_DIR/CONTROL.txt" ]; then
            printf '%s\n' "CONTROL.txt" >> "$manifest_file"
            input_entry="$(awk -F= '
                /^[[:space:]]*input_file[[:space:]]*=/ {
                    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
                    print $2
                    exit
                }
            ' "$SLURM_SUBMIT_DIR/CONTROL.txt")"
            input_entry="${input_entry:-input.txt}"
            if [[ "$input_entry" = /* ]]; then
                echo "Fresh staging fallback: CONTROL.txt references absolute input_file '$input_entry'."
                staged_count=0
            else
                printf '%s\n' "$input_entry" >> "$manifest_file"
            fi
        fi
        for optional_file in input.txt input.xyz start.txt co2.xyz; do
            [ -e "$SLURM_SUBMIT_DIR/$optional_file" ] && printf '%s\n' "$optional_file" >> "$manifest_file"
        done
        awk 'NF && !seen[$0]++' "$manifest_file" > "${manifest_file}.tmp" && mv "${manifest_file}.tmp" "$manifest_file"

        if [ -s "$manifest_file" ]; then
            while IFS= read -r rel_path; do
                [ -e "$SLURM_SUBMIT_DIR/$rel_path" ] || continue
                mkdir -p "$RUN_DIR/$(dirname "$rel_path")"
                cp -a "$SLURM_SUBMIT_DIR/$rel_path" "$RUN_DIR/$rel_path"
                staged_count=$((staged_count + 1))
            done < "$manifest_file"
        fi

        if [ -n "${DELFIN_XYZ_FILE:-}" ] && [ -f "$DELFIN_XYZ_FILE" ]; then
            cp -a "$DELFIN_XYZ_FILE" "$RUN_DIR/$(basename "$DELFIN_XYZ_FILE")"
        fi

        if [ "$staged_count" -gt 0 ]; then
            echo "Staged ${staged_count} explicit input paths to scratch."
        else
            echo "Falling back to rsync fresh profile because no safe explicit input set was found."
            rsync -a \
                --exclude='.git/' \
                --exclude='.venv/' \
                --exclude='__pycache__/' \
                --exclude='.ipynb_checkpoints/' \
                --exclude='.orca_iso*' \
                --exclude='delfin_*.out' \
                --exclude='delfin_*.err' \
                --exclude='*.out' \
                --exclude='*.gbw' \
                --exclude='*.hess' \
                --exclude='*.cube' \
                --exclude='*.densitiesinfo' \
                --exclude='*.property.txt' \
                --exclude='*.png' \
                --exclude='*.svg' \
                --exclude='*.pdf' \
                --exclude='*.docx' \
                --exclude='*.json' \
                --exclude='*.fprint' \
                --exclude='*.tmp' \
                --exclude='*_OCCUPIER/' \
                --exclude='ESD/' \
                --exclude='builder/' \
                "$SLURM_SUBMIT_DIR"/ "$RUN_DIR"/ 2>/dev/null || true
        fi
    else
        rsync -a \
            --exclude='.git/' \
            --exclude='.venv/' \
            --exclude='__pycache__/' \
            --exclude='.ipynb_checkpoints/' \
            --exclude='.orca_iso*' \
            "$SLURM_SUBMIT_DIR"/ "$RUN_DIR"/ 2>/dev/null || true
    fi
}

# Rescue files from any active .orca_iso* dirs into their parent directory.
# Non-.inp files are overwritten (latest ORCA state), .inp files are never
# overwritten (cp -n) so the original input is preserved.
rescue_iso_files() {
    [ -d "$RUN_DIR" ] || return 0
    while IFS= read -r -d '' iso_dir; do
        parent_dir="$(dirname "$iso_dir")"
        find "$iso_dir" -maxdepth 1 -type f | while IFS= read -r f; do
            fname="$(basename "$f")"
            if [[ "$fname" == *.inp ]]; then
                cp -n "$f" "$parent_dir/" 2>/dev/null || true   # no-clobber for .inp
            else
                cp -f "$f" "$parent_dir/" 2>/dev/null || true   # overwrite for everything else
            fi
        done
    done < <(find "$RUN_DIR" -name '.orca_iso*' -type d -print0 2>/dev/null)
}

# Delete GOAT per-conformer outputs that should never be copied back to HOME.
purge_xtb_goat_out_files() {
    [ -d "$RUN_DIR" ] || return 0
    find "$RUN_DIR" -type f -name 'XTB_GOAT.goat.*.out' -delete 2>/dev/null || true
}

SYNC_STAMP_MINIMAL="$DELFIN_SCRATCH/.last_result_sync.minimal"

collect_sync_results() {
    local list_file="$1"
    local profile="${2:-full}"
    local stamp_path=""

    case "$profile" in
        minimal)
            stamp_path="$SYNC_STAMP_MINIMAL"
            if [ -f "$stamp_path" ]; then
                find "$RUN_DIR" \
                    \( -type f -o -type l \) \
                    ! -path '*/.orca_iso*/*' \
                    ! -path '*/__pycache__/*' \
                    ! -path '*/.ipynb_checkpoints/*' \
                    ! -name '*.tmp' \
                    ! -name '*.pyc' \
                    ! -name '*.pyo' \
                    ! -name 'XTB_GOAT.goat.*.out' \
                    \( -name '*.out' -o -name '*.err' -o -name '*.log' -o -name '*.json' -o -name '*.txt' -o -name '*.xyz' -o -name '*.inp' -o -name '*.png' -o -name '*.svg' -o -name '*.pdf' -o -name '*.docx' -o -name '*.csv' -o -name '*.dat' \) \
                    -newer "$stamp_path" \
                    -printf '%P\n' | sort -u > "$list_file"
            else
                find "$RUN_DIR" \
                    \( -type f -o -type l \) \
                    ! -path '*/.orca_iso*/*' \
                    ! -path '*/__pycache__/*' \
                    ! -path '*/.ipynb_checkpoints/*' \
                    ! -name '*.tmp' \
                    ! -name '*.pyc' \
                    ! -name '*.pyo' \
                    ! -name 'XTB_GOAT.goat.*.out' \
                    \( -name '*.out' -o -name '*.err' -o -name '*.log' -o -name '*.json' -o -name '*.txt' -o -name '*.xyz' -o -name '*.inp' -o -name '*.png' -o -name '*.svg' -o -name '*.pdf' -o -name '*.docx' -o -name '*.csv' -o -name '*.dat' \) \
                    -printf '%P\n' | sort -u > "$list_file"
            fi
            ;;
        *)
            find "$RUN_DIR" \
                \( -type f -o -type l \) \
                ! -path '*/.orca_iso*/*' \
                ! -path '*/__pycache__/*' \
                ! -path '*/.ipynb_checkpoints/*' \
                ! -name '*.tmp' \
                ! -name '*.pyc' \
                ! -name '*.pyo' \
                ! -name 'XTB_GOAT.goat.*.out' \
                -printf '%P\n' | sort -u > "$list_file"
            ;;
    esac
}

sync_results_back() {
    local sync_label="${1:-sync}"
    local sync_profile="${2:-full}"
    local list_file="$DELFIN_SCRATCH/.sync_files.txt"
    local file_count

    if [ ! -d "$RUN_DIR" ]; then
        echo "WARNING: RUN_DIR not found, nothing to copy."
        return 0
    fi

    purge_xtb_goat_out_files
    if [ "$sync_profile" = "full" ]; then
        rescue_iso_files
    fi
    collect_sync_results "$list_file" "$sync_profile"

    if [ ! -s "$list_file" ]; then
        echo "No updated result files to copy (${sync_label})."
        if [ "$sync_profile" = "minimal" ]; then
            touch "$SYNC_STAMP_MINIMAL" 2>/dev/null || true
        fi
        return 0
    fi

    file_count="$(wc -l < "$list_file" 2>/dev/null || echo 0)"
    rsync -a --files-from="$list_file" "$RUN_DIR"/ "$SLURM_SUBMIT_DIR"/ 2>/dev/null || true
    if [ "$sync_profile" = "minimal" ]; then
        touch "$SYNC_STAMP_MINIMAL" 2>/dev/null || true
    fi
    echo "Copied ${file_count} result files (${sync_label}, profile=${sync_profile})."
}

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

    # CRITICAL: Rescue .orca_iso files first, then perform a full result sync.
    echo "Copying results back to $SLURM_SUBMIT_DIR..."
    sync_results_back "signal-${signal_name}" "full"

    # Cleanup scratch (only after successful copy)
    rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR" "${VENV_LOCAL:-}" 2>/dev/null || true
    echo "Cleanup completed at $(date)"
    exit 1
}

# Periodic copy-back to mitigate hard timeouts.
periodic_copy() {
    while true; do
        sleep 7200
        if [ -d "$RUN_DIR" ]; then
            sync_results_back "periodic" "minimal"
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

# Schedule final backup 10 minutes before timeout
SAFETY_MARGIN=600  # 10 minutes in seconds
schedule_final_backup() {
    local walltime_sec
    walltime_sec=$(get_walltime_seconds)
    if [ "$walltime_sec" -le "$SAFETY_MARGIN" ]; then
        echo "Walltime too short for scheduled final backup, skipping."
        return
    fi
    local wait_time=$((walltime_sec - SAFETY_MARGIN))
    echo "Scheduled final backup in $((wait_time / 3600))h $((wait_time % 3600 / 60))m (10 min before timeout)"
    sleep "$wait_time"
    echo ""
    echo "========================================"
    echo "Final backup 5 min before timeout: $(date)"
    echo "========================================"
    if [ -d "$RUN_DIR" ]; then
        sync_results_back "final-backup" "full"
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
echo "DELFIN Version: $(delfin --version 2>&1 || echo 'unknown')"
if [ "${DELFIN_SHOW_GIT_INFO:-0}" = "1" ] && [ -d "$DELFIN_DIR/.git" ]; then
    cd "$DELFIN_DIR"
    echo "Git Branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'N/A')"
    echo "Git Commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'N/A')"
    cd - > /dev/null
fi
echo ""

# Fresh DELFIN/BUILD/GUPPY/CO2 runs only need current inputs; recalc and raw ORCA
# runs keep the full workspace copy because they may depend on existing outputs.
COPY_PROFILE="full"
case "$MODE" in
    delfin|build|guppy|hyperpol_xtb|tadf_xtb|delfin-co2-chain)
        COPY_PROFILE="fresh"
        ;;
esac
copy_workspace_to_scratch "$COPY_PROFILE"
# Remove output files from previous runs (if any)
rm -f "$RUN_DIR"/delfin_*.out "$RUN_DIR"/delfin_*.err 2>/dev/null || true

cd "$RUN_DIR"
purge_xtb_goat_out_files
touch "$SYNC_STAMP_MINIMAL" 2>/dev/null || true

# Resolve the Python that belongs to the delfin installation
# (bare `python` may point to a different conda/mamba env)
DELFIN_PYTHON="$(head -1 "$(command -v delfin)" | sed 's/^#!//')"
if [ ! -x "$DELFIN_PYTHON" ]; then
    DELFIN_PYTHON="$(command -v python)"
fi

configure_local_xtb4stda_home() {
    local stage_base="${STAGE_BASE:-}"
    local short_link=""
    local target=""
    if [ -z "$stage_base" ] || [ ! -d "$stage_base" ]; then
        return 0
    fi
    short_link="$stage_base/xtb4stda_${SLURM_JOB_ID}"
    target="$("$DELFIN_PYTHON" - <<'PY' 2>/dev/null
from pathlib import Path
from delfin.qm_runtime import get_qm_tools_root
root = Path(get_qm_tools_root())
target = root / "share" / "xtb4stda"
print(target.resolve() if target.exists() else "")
PY
)"
    if [ -n "$target" ] && [ -d "$target" ]; then
        rm -f "$short_link"
        ln -s "$target" "$short_link" 2>/dev/null || true
        if [ -L "$short_link" ] || [ -d "$short_link" ]; then
            export XTB4STDAHOME="$short_link"
        fi
    fi
}

configure_local_xtb4stda_home

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
        export DELFIN_SMART_RECALC="${DELFIN_SMART_RECALC:-1}"
        delfin --recalc
        EXIT_CODE=$?
        ;;
    delfin-recalc-classic)
        echo "Starting DELFIN --recalc (classic marker mode)..."
        export DELFIN_SMART_RECALC=0
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
        "$DELFIN_PYTHON" -m delfin.build_up_complex input.txt --goat --directory builder --multiplicity "$BUILD_MULT" --verbose
        EXIT_CODE=$?
        ;;
    guppy)
        GUPPY_RUNS="${GUPPY_RUNS:-20}"
        GUPPY_PAL="${GUPPY_PAL:-${DELFIN_PAL:-${SLURM_CPUS_PER_TASK:-40}}}"
        GUPPY_MAXCORE="${GUPPY_MAXCORE:-${DELFIN_MAXCORE:-6000}}"
        GUPPY_PARALLEL_JOBS="${GUPPY_PARALLEL_JOBS:-4}"
        GUPPY_GOAT_TOPK="${GUPPY_GOAT_TOPK:-3}"
        GUPPY_GOAT_PARALLEL_JOBS="${GUPPY_GOAT_PARALLEL_JOBS:-$GUPPY_PARALLEL_JOBS}"
        echo "Starting GUPPY SMILES sampling..."
        echo "  Runs:         $GUPPY_RUNS"
        echo "  Charge:       auto (derived from SMILES)"
        echo "  Multiplicity: 1 (fixed closed-shell)"
        echo "  PAL total:    $GUPPY_PAL"
        echo "  Maxcore:      $GUPPY_MAXCORE"
        echo "  Parallel:     $GUPPY_PARALLEL_JOBS runs"
        echo "  GOAT top-k:   $GUPPY_GOAT_TOPK"
        echo "  GOAT parallel:$GUPPY_GOAT_PARALLEL_JOBS runs"
        GUPPY_CMD=(
            "$DELFIN_PYTHON" -m delfin.guppy_sampling input.txt
            --runs "$GUPPY_RUNS" \
            --pal "$GUPPY_PAL" \
            --maxcore "$GUPPY_MAXCORE" \
            --parallel-jobs "$GUPPY_PARALLEL_JOBS" \
            --goat-topk "$GUPPY_GOAT_TOPK" \
            --goat-parallel-jobs "$GUPPY_GOAT_PARALLEL_JOBS" \
            --output GUPPY_try.xyz
        )
        "${GUPPY_CMD[@]}"
        EXIT_CODE=$?
        ;;
    hyperpol_xtb)
        XYZ_FILE="${DELFIN_XYZ_FILE:-}"
        WORKFLOW_LABEL="${DELFIN_WORKFLOW_LABEL:-hyperpol_xtb}"
        TARGET_WORKDIR="$RUN_DIR"
        if [ -n "$XYZ_FILE" ] && [ -f "$RUN_DIR/$(basename "$XYZ_FILE")" ]; then
            XYZ_FILE="$RUN_DIR/$(basename "$XYZ_FILE")"
        fi
        if [ -z "$XYZ_FILE" ]; then
            echo "ERROR: DELFIN_XYZ_FILE not set for hyperpol_xtb mode"
            EXIT_CODE=1
        elif [ ! -f "$XYZ_FILE" ]; then
            echo "ERROR: XYZ file not found for hyperpol_xtb mode: $XYZ_FILE"
            EXIT_CODE=1
        else
            echo "Starting browser hyperpol_xtb workflow..."
            echo "  XYZ:          $XYZ_FILE"
            echo "  Label:        $WORKFLOW_LABEL"
            echo "  PAL:          ${DELFIN_PAL:-4}"
            echo "  Maxcore:      ${DELFIN_MAXCORE:-1000}"
            echo "  Target Dir:   $TARGET_WORKDIR"
            "$DELFIN_PYTHON" -m delfin.dashboard.browser_workflows hyperpol_xtb \
                --xyz-file "$XYZ_FILE" \
                --label "$WORKFLOW_LABEL" \
                --workdir "$TARGET_WORKDIR" \
                --engine std2 \
                --preopt none \
                --static-only \
                --energy-window 15 \
                --pal "${DELFIN_PAL:-4}" \
                --maxcore "${DELFIN_MAXCORE:-1000}" \
                --json-out "$TARGET_WORKDIR/hyperpol_xtb_summary.json" \
                | tee "$TARGET_WORKDIR/hyperpol_xtb.output"
            EXIT_CODE=${PIPESTATUS[0]}
        fi
        ;;
    tadf_xtb)
        XYZ_FILE="${DELFIN_XYZ_FILE:-}"
        WORKFLOW_LABEL="${DELFIN_WORKFLOW_LABEL:-tadf_xtb}"
        TARGET_WORKDIR="$RUN_DIR"
        if [ -n "$XYZ_FILE" ] && [ -f "$RUN_DIR/$(basename "$XYZ_FILE")" ]; then
            XYZ_FILE="$RUN_DIR/$(basename "$XYZ_FILE")"
        fi
        if [ -z "$XYZ_FILE" ]; then
            echo "ERROR: DELFIN_XYZ_FILE not set for tadf_xtb mode"
            EXIT_CODE=1
        elif [ ! -f "$XYZ_FILE" ]; then
            echo "ERROR: XYZ file not found for tadf_xtb mode: $XYZ_FILE"
            EXIT_CODE=1
        else
            echo "Starting browser tadf_xtb workflow..."
            echo "  XYZ:          $XYZ_FILE"
            echo "  Label:        $WORKFLOW_LABEL"
            echo "  PAL:          ${DELFIN_PAL:-4}"
            echo "  Maxcore:      ${DELFIN_MAXCORE:-1000}"
            echo "  Target Dir:   $TARGET_WORKDIR"
            "$DELFIN_PYTHON" -m delfin.dashboard.browser_workflows tadf_xtb \
                --xyz-file "$XYZ_FILE" \
                --label "$WORKFLOW_LABEL" \
                --workdir "$TARGET_WORKDIR" \
                --pal "${DELFIN_PAL:-4}" \
                --maxcore "${DELFIN_MAXCORE:-1000}" \
                --json-out "$TARGET_WORKDIR/tadf_xtb_summary.json" \
                | tee "$TARGET_WORKDIR/tadf_xtb.output"
            EXIT_CODE=${PIPESTATUS[0]}
        fi
        ;;
    delfin-co2-chain)
        CO2_SPECIES_DELTA="${DELFIN_CO2_SPECIES_DELTA:-0}"
        echo "Starting DELFIN + CO2 Coordinator chain..."
        echo "  CO2 Species Delta: $CO2_SPECIES_DELTA"

        # Step 1: DELFIN
        echo "=== Step 1/2: Running DELFIN ==="
        delfin
        DELFIN_EXIT=$?
        if [ "$DELFIN_EXIT" -ne 0 ]; then
            echo "ERROR: DELFIN failed (exit $DELFIN_EXIT). CO2 Coordinator skipped."
            EXIT_CODE=$DELFIN_EXIT
        else
            # Step 2: Chain Setup + CO2 Coordinator
            echo "=== Step 2/2: CO2 Coordinator ==="
            "$DELFIN_PYTHON" -m delfin.co2.chain_setup "$CO2_SPECIES_DELTA"
            SETUP_EXIT=$?
            if [ "$SETUP_EXIT" -ne 0 ]; then
                echo "ERROR: CO2 setup failed (exit $SETUP_EXIT)"
                EXIT_CODE=$SETUP_EXIT
            else
                cd CO2_coordination
                delfin co2
                EXIT_CODE=$?
                cd ..
            fi
        fi
        ;;
    *)
        echo "ERROR: Unknown mode: $MODE"
        echo "       Valid modes: delfin, delfin-recalc, delfin-recalc-classic, delfin-recalc-override, orca, build, guppy, hyperpol_xtb, tadf_xtb, delfin-co2-chain, auto"
        EXIT_CODE=1
        ;;
esac
set -e

echo ""
echo "========================================"
echo "Job finished: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "========================================"

# Copy results back (always rescue .orca_iso contents before the final full sync)
sync_results_back "job-end" "full"

# Cleanup scratch
rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR" "${VENV_LOCAL:-}"

exit $EXIT_CODE
