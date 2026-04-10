#!/usr/bin/env bash
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

set -euo pipefail

# ========================================================================
# DELFIN Central SLURM Submit Script
# ========================================================================
#
# Single submit script for all SLURM clusters. Site-specific features
# (module loading, ORCA staging, venv tarball, runtime wheel cache) are
# activated via environment variables so the same script works everywhere.
#
# MODES (set via DELFIN_MODE environment variable):
#   delfin | delfin-recalc | delfin-recalc-classic | orca | build | guppy
#   | hyperpol_xtb | tadf_xtb | censo_anmr | delfin-co2-chain | auto (default: auto)
#
# SITE-SPECIFIC FEATURES (all off by default):
#   DELFIN_AUTO_RESOURCES=1  Parse CONTROL.txt to derive sbatch args
#   DELFIN_MODULES="mod1 mod2"  Modules to load (module purge + load)
#   DELFIN_STAGE_ORCA=1      Stage ORCA + OpenMPI to node-local SSD
#   DELFIN_STAGE_VENV=1      Unpack venv tarball to SSD
#   DELFIN_RUNTIME_CACHE=1   Build/use runtime wheel cache
#
# RESOURCE PARAMETERS (via sbatch command-line overrides):
#   --time, --ntasks, --cpus-per-task, --mem, --job-name
#
# STAGE-IN / STAGE-OUT:
#   DELFIN_STAGE_IO=1 (default)  Copy workspace to local SSD
#   DELFIN_SYNC_INTERVAL=900     Periodic sync interval (seconds)
#
# AUTO RESOURCES (when DELFIN_AUTO_RESOURCES=1):
#   DELFIN_PARTITION=cpu         Default partition
#   DELFIN_HIGHMEM_PARTITION=highmem  Partition for high-memory jobs
#   DELFIN_NODE_CORES=96         Cores per node
#   DELFIN_NODE_MEM_MB=393216    RAM per node in MB
#   DELFIN_HIGHMEM_MB=2359296    RAM per highmem node in MB
#   DELFIN_FORCE_EXCLUSIVE=0     Force --exclusive
#
# ENVIRONMENT VARIABLES (passed through to job):
#   DELFIN_MODE, DELFIN_JOB_NAME, DELFIN_INP_FILE, DELFIN_XYZ_FILE,
#   DELFIN_WORKFLOW_LABEL, DELFIN_PAL, DELFIN_MAXCORE, DELFIN_ORCA_BASE,
#   BUILD_MULTIPLICITY, GUPPY_*, DELFIN_CO2_SPECIES_DELTA, etc.
#
# ========================================================================

# ======================================================================
# Section 0: Auto-resource detection (optional, pre-SLURM)
# ======================================================================
if [ -z "${SLURM_JOB_ID:-}" ] && [ "${DELFIN_AUTO_RESOURCES:-0}" = "1" ]; then
    control_file="${DELFIN_CONTROL:-$PWD/CONTROL.txt}"
    pal_default=40
    maxcore_default=6000
    pal="$pal_default"
    maxcore="$maxcore_default"
    if [ -f "$control_file" ]; then
        pal_found="$(awk -F= '/^[[:space:]]*PAL[[:space:]]*=/ {gsub(/[^0-9]/,"",$2); print $2; exit}' "$control_file")"
        maxcore_found="$(awk -F= '/^[[:space:]]*maxcore[[:space:]]*=/ {gsub(/[^0-9]/,"",$2); print $2; exit}' "$control_file")"
        [ -n "$pal_found" ] && pal="$pal_found"
        [ -n "$maxcore_found" ] && maxcore="$maxcore_found"
    fi

    ntasks="1"
    cpus_per_task="$pal"
    mem_mb=$((pal * maxcore))

    node_cores="${DELFIN_NODE_CORES:-96}"
    node_mem_mb="${DELFIN_NODE_MEM_MB:-$((384 * 1024))}"
    highmem_mb="${DELFIN_HIGHMEM_MB:-$((2304 * 1024))}"
    partition="${DELFIN_PARTITION:-cpu}"
    highmem_partition="${DELFIN_HIGHMEM_PARTITION:-highmem}"
    exclusive=0

    if [ "$maxcore" -ge 9500 ] && [ "$partition" != "$highmem_partition" ]; then
        partition="$highmem_partition"
        node_mem_mb="$highmem_mb"
    elif [ "$mem_mb" -gt "$node_mem_mb" ] && [ "$partition" != "$highmem_partition" ]; then
        partition="$highmem_partition"
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
        echo "Set DELFIN_PARTITION=${highmem_partition} (or adjust PAL/maxcore) and submit again."
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
    [ "$exclusive" -eq 1 ] && sbatch_args+=(--exclusive)

    echo "Auto resources: PAL=${pal} maxcore=${maxcore} -> --ntasks=${ntasks} --cpus-per-task=${cpus_per_task} --mem=${mem_mb}M --partition=${partition}"
    [ "$exclusive" -eq 1 ] && echo "Auto resources: requesting --exclusive (near full node)"
    sbatch "${sbatch_args[@]}" "$0" "$@"
    exit 0
fi

# ======================================================================
# Section 1: Configuration
# ======================================================================
MODE="${DELFIN_MODE:-auto}"
INP_FILE="${DELFIN_INP_FILE:-}"
OVERRIDE="${DELFIN_OVERRIDE:-}"
DISPLAY_JOB_NAME="${DELFIN_JOB_NAME:-${SLURM_JOB_NAME:-delfin_job}}"

export DELFIN_JOB_ID="${DELFIN_JOB_ID:-${SLURM_JOB_ID:-0}}"
export DELFIN_JOB_NAME="$DISPLAY_JOB_NAME"

STAGE_IO="${DELFIN_STAGE_IO:-1}"
SYNC_INTERVAL="${DELFIN_SYNC_INTERVAL:-900}"

# ======================================================================
# Section 2: PAL / SLURM sanity check
# ======================================================================
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
        fi
        if [ "$SLURM_TASKS" -ne 1 ] && [ "$SLURM_CPT" -ne "$PAL_FROM_CONTROL" ]; then
            echo "WARNING: CONTROL.txt PAL=${PAL_FROM_CONTROL} but SLURM has ntasks=${SLURM_TASKS}, cpus-per-task=${SLURM_CPT}."
            if [ -n "$MAXCORE_FROM_CONTROL" ]; then
                echo "         Suggested: sbatch --ntasks=1 --cpus-per-task=${PAL_FROM_CONTROL} --mem=$((PAL_FROM_CONTROL * MAXCORE_FROM_CONTROL))M"
            fi
        fi
    fi
fi

# ======================================================================
# Section 3: Module loading (optional)
# ======================================================================
if [ -n "${DELFIN_MODULES:-}" ]; then
    module purge 2>/dev/null || true
    for _mod in $DELFIN_MODULES; do
        module load "$_mod" 2>/dev/null || echo "WARNING: Failed to load module $_mod"
    done
fi

# Initialize LD_LIBRARY_PATH if not set
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ======================================================================
# Section 4: Software base auto-detection (optional)
# ======================================================================
# Locate the base directory that contains delfin, orca, openmpi, etc.
# Set DELFIN_SOFTWARE_BASE explicitly or let the script search upwards.
SOFTWARE_BASE="${DELFIN_SOFTWARE_BASE:-}"
DELFIN_DIR=""
if [ -z "$SOFTWARE_BASE" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    _search_dir="$SLURM_SUBMIT_DIR"
    while [ "$_search_dir" != "/" ] && [ ! -d "$_search_dir/software/delfin" ]; do
        _search_dir="$(dirname "$_search_dir")"
    done
    if [ -d "$_search_dir/software/delfin" ]; then
        SOFTWARE_BASE="$_search_dir/software"
    fi
fi

if [ -n "$SOFTWARE_BASE" ]; then
    DELFIN_DIR="$SOFTWARE_BASE/delfin"
fi

# Normalise SLURM scratch hints early
if [ -z "${TMPDIR:-}" ] && [ -n "${SLURM_TMPDIR:-}" ]; then
    export TMPDIR="$SLURM_TMPDIR"
fi
STAGE_BASE="${TMPDIR:-${BEEOND_MOUNTPOINT:-}}"

# ======================================================================
# Section 5: ORCA + OpenMPI staging to SSD (optional)
# ======================================================================
ORCA_LOCAL=""
OMPI_LOCAL=""
OMPI_HOME="${DELFIN_OMPI_HOME:-}"
ORCA_BASE="${DELFIN_ORCA_BASE:-}"

if [ -n "$SOFTWARE_BASE" ]; then
    [ -z "$OMPI_HOME" ] && [ -d "$SOFTWARE_BASE/openmpi-4.1.8" ] && OMPI_HOME="$SOFTWARE_BASE/openmpi-4.1.8"
    if [ -z "$ORCA_BASE" ]; then
        # Auto-detect ORCA directory
        for _candidate in "$SOFTWARE_BASE"/orca_*; do
            [ -d "$_candidate" ] && ORCA_BASE="$_candidate" && break
        done
    fi
fi

if [ "${DELFIN_STAGE_ORCA:-0}" = "1" ] && [ -n "${STAGE_BASE:-}" ] && [ -d "${STAGE_BASE}" ]; then
    if [ -n "$OMPI_HOME" ] && [ -d "$OMPI_HOME" ]; then
        OMPI_LOCAL="$STAGE_BASE/openmpi_${SLURM_JOB_ID}"
        echo "Staging OpenMPI to local SSD ($OMPI_LOCAL)..."
        cp -a "$OMPI_HOME" "$OMPI_LOCAL" && OMPI_HOME="$OMPI_LOCAL" \
            || echo "WARNING: OpenMPI staging failed, using HOME."
    fi
    if [ -n "$ORCA_BASE" ] && [ -d "$ORCA_BASE" ]; then
        ORCA_LOCAL="$STAGE_BASE/orca_${SLURM_JOB_ID}"
        echo "Staging ORCA to local SSD ($ORCA_LOCAL)..."
        cp -a "$ORCA_BASE" "$ORCA_LOCAL" && ORCA_BASE="$ORCA_LOCAL" \
            || echo "WARNING: ORCA staging failed, using HOME."
    fi
    [ -n "$OMPI_LOCAL" ] || [ -n "$ORCA_LOCAL" ] && echo "ORCA + OpenMPI staged to local SSD."
fi

if [ -n "$OMPI_HOME" ] && [ -d "$OMPI_HOME" ]; then
    echo "Using custom OpenMPI from $OMPI_HOME"
    export PATH="$OMPI_HOME/bin:$PATH"
    export LD_LIBRARY_PATH="$OMPI_HOME/lib:$LD_LIBRARY_PATH"
fi

if [ -n "$ORCA_BASE" ] && [ -d "$ORCA_BASE" ]; then
    export PATH="$ORCA_BASE:$PATH"
    export LD_LIBRARY_PATH="$ORCA_BASE:$LD_LIBRARY_PATH"
    export ORCA_DIR="$ORCA_BASE"
    export ORCA_PLOT="$ORCA_BASE/orca_plot"
    export DELFIN_ORCA_BASE="$ORCA_BASE"
fi

# ======================================================================
# Section 6: MPI configuration for ORCA stability
# ======================================================================
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp,vader
export OMPI_MCA_mpi_show_mca_params_file=0
export OMPI_MCA_mpi_yield_when_idle=1
export OMPI_MCA_coll_hcoll_enable=0
unset OMPI_MCA_btl_tcp_if_include 2>/dev/null || true
unset OMPI_MCA_btl_tcp_if_exclude 2>/dev/null || true
export OMPI_MCA_hwloc_base_binding_policy=none
export OMPI_MCA_rmaps_base_mapping_policy=core
export OMPI_MCA_rmaps_base_oversubscribe=true

# ======================================================================
# Section 7: Venv tarball staging (optional)
# ======================================================================
VENV_LOCAL=""
if [ "${DELFIN_STAGE_VENV:-0}" = "1" ] && [ -n "$DELFIN_DIR" ] && [ -n "${STAGE_BASE:-}" ] && [ -d "${STAGE_BASE}" ]; then
    VENV_TAR="${DELFIN_VENV_TAR:-$DELFIN_DIR/delfin_venv.tar}"
    if [ -f "$VENV_TAR" ] || [ -d "$DELFIN_DIR/.venv" ]; then
        VENV_LOCAL="$STAGE_BASE/delfin_venv_${SLURM_JOB_ID}"
        if [ ! -f "$VENV_TAR" ]; then
            echo "WARNING: $VENV_TAR not found, creating it now (one-time)..."
            tar -cf "$VENV_TAR" -C "$DELFIN_DIR" .venv/
            echo "Created $VENV_TAR"
        fi
        echo "Unpacking venv tarball to local SSD ($VENV_LOCAL) to minimise HOME I/O..."
        mkdir -p "$VENV_LOCAL"
        tar -xf "$VENV_TAR" --strip-components=1 -C "$VENV_LOCAL"
        sed -i "s|$DELFIN_DIR/.venv|$VENV_LOCAL|g" "$VENV_LOCAL/bin/activate" 2>/dev/null || true
        source "$VENV_LOCAL/bin/activate"
        echo "venv loaded from local SSD."
    fi
elif [ -n "$DELFIN_DIR" ] && [ -d "$DELFIN_DIR/.venv" ]; then
    echo "Using repository venv directly."
    source "$DELFIN_DIR/.venv/bin/activate"
fi

# ======================================================================
# Section 8: Runtime wheel cache (optional)
# ======================================================================
RUNTIME_CACHE_DIR="${DELFIN_RUNTIME_CACHE_DIR:-${DELFIN_DIR:+$DELFIN_DIR/.runtime_cache}}"
RUNTIME_KEY=""
RUNTIME_BUILD_MODE=""
RUNTIME_WHEEL=""

if [ "${DELFIN_RUNTIME_CACHE:-0}" = "1" ] && [ -n "$DELFIN_DIR" ]; then

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
    [ -f "$DELFIN_DIR/README.md" ] && cp -a "$DELFIN_DIR/README.md" "$build_dir/"
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
        [ "$lock_fd_opened" -eq 1 ] && { flock -u 9; exec 9>&-; }
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
        [ "$lock_fd_opened" -eq 1 ] && { flock -u 9; exec 9>&-; }
        exit 1
    fi

    mv "$built_wheel" "$runtime_root/"
    printf 'runtime_key=%s\nbuild_mode=%s\ncreated=%s\n' \
        "$RUNTIME_KEY" "$RUNTIME_BUILD_MODE" "$(date -Is)" > "$runtime_root/runtime.meta"
    rm -rf "$tmp_build_dir" "$tmp_wheel_dir"

    RUNTIME_WHEEL="$(find "$runtime_root" -maxdepth 1 -type f -name 'delfin_complat-*.whl' | sort | tail -n 1)"
    [ "$lock_fd_opened" -eq 1 ] && { flock -u 9; exec 9>&-; }
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

fi  # end DELFIN_RUNTIME_CACHE block

# ======================================================================
# Section 9: Environment variables
# ======================================================================
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export DELFIN_ORCA_PROGRESS=0
export MPLBACKEND=Agg
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PIP_DISABLE_PIP_VERSION_CHECK=1
export PIP_NO_INPUT=1

# Redirect cache directories to SSD when available
if [ -n "${STAGE_BASE:-}" ] && [ -d "${STAGE_BASE}" ]; then
    export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$STAGE_BASE/.cache_${SLURM_JOB_ID:-0}}"
    export MPLCONFIGDIR="${MPLCONFIGDIR:-$STAGE_BASE/.mplconfig_${SLURM_JOB_ID:-0}}"
    export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$STAGE_BASE/.pip_${SLURM_JOB_ID:-0}}"
    mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$PIP_CACHE_DIR"
fi

# ======================================================================
# Section 10: Python binary detection
# ======================================================================
PYTHON_BIN="${DELFIN_PYTHON_BIN:-}"
if [ -z "$PYTHON_BIN" ]; then
    # Try the delfin entry-point shebang first
    _delfin_bin="$(command -v delfin 2>/dev/null || true)"
    if [ -n "$_delfin_bin" ] && [ -x "$_delfin_bin" ]; then
        _shebang="$(head -1 "$_delfin_bin" 2>/dev/null || true)"
        if [ "${_shebang#\#!}" != "$_shebang" ]; then
            _candidate="${_shebang#\#!}"
            _candidate="${_candidate%% *}"
            if [ "$_candidate" = "/usr/bin/env" ]; then
                _env_py="$(printf '%s\n' "${_shebang#\#!}" | awk '{print $2}')"
                [ -n "$_env_py" ] && _candidate="$(command -v "$_env_py" 2>/dev/null || true)"
            fi
            [ -n "$_candidate" ] && [ -x "$_candidate" ] && PYTHON_BIN="$_candidate"
        fi
    fi
fi
if [ -z "$PYTHON_BIN" ]; then
    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN="$(command -v python)"
    elif command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="$(command -v python3)"
    else
        echo "ERROR: python/python3 not found in PATH." >&2
        exit 1
    fi
fi

# ======================================================================
# Section 11: Scratch directory setup
# ======================================================================
ORIGIN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

if [ -n "${BEEOND_MOUNTPOINT:-}" ]; then
    export DELFIN_SCRATCH="$BEEOND_MOUNTPOINT/delfin_${SLURM_JOB_ID:-0}"
    export ORCA_TMPDIR="$BEEOND_MOUNTPOINT/orca_${SLURM_JOB_ID:-0}"
    echo "Using BeeOND (parallel SSD): $BEEOND_MOUNTPOINT (job-isolated: ${SLURM_JOB_ID:-0})"
elif [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
    export DELFIN_SCRATCH="$TMPDIR/delfin_${SLURM_JOB_ID:-0}"
    export ORCA_TMPDIR="$TMPDIR/orca_${SLURM_JOB_ID:-0}"
    echo "Using local SSD (TMPDIR): $TMPDIR (job-isolated: ${SLURM_JOB_ID:-0})"
else
    export DELFIN_SCRATCH="/scratch/${USER:-$LOGNAME}/delfin_${SLURM_JOB_ID:-0}"
    export ORCA_TMPDIR="/scratch/${USER:-$LOGNAME}/orca_${SLURM_JOB_ID:-0}"
    echo "WARNING: Using /scratch (network filesystem) - consider using TMPDIR for better I/O"
fi

export ORCA_SCRDIR="$ORCA_TMPDIR"
export DELFIN_RUN_TOKEN="slurm_${SLURM_JOB_ID:-0}"

RUN_DIR="$DELFIN_SCRATCH/run"
mkdir -p "$RUN_DIR" "$ORCA_TMPDIR"

# ======================================================================
# Section 12: SLURM diagnostics
# ======================================================================
write_slurm_reason() {
    local out_file="$RUN_DIR/slurm_reason.txt"
    {
        echo "Timestamp: $(date)"
        echo "Job ID: ${SLURM_JOB_ID:-unknown}"
        echo "Node: ${SLURM_JOB_NODELIST:-unknown}"
        if command -v scontrol >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
            scontrol show job "$SLURM_JOB_ID" 2>/dev/null | tr ' ' '\n' | grep -E 'JobId=|JobName=|State=|Reason=|TimeLimit=|RunTime=|TimeMin=|StartTime=|EndTime=' || true
        fi
        if command -v sacct >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
            sacct -j "$SLURM_JOB_ID" --format=JobID,JobName%30,State,ExitCode,Elapsed,Timelimit,ReqMem,MaxRSS,NodeList,Reason -P 2>/dev/null || true
        fi
    } > "$out_file" 2>/dev/null || true
}

write_slurm_reason

# ======================================================================
# Section 13: Auto-detect mode
# ======================================================================
if [ "$MODE" = "auto" ]; then
    if [ -f "$ORIGIN_DIR/CONTROL.txt" ] && [ -f "$ORIGIN_DIR/input.txt" ]; then
        MODE="delfin"
        echo "Auto-detected mode: DELFIN (CONTROL.txt + input.txt found)"
    elif ls "$ORIGIN_DIR"/*.inp 1>/dev/null 2>&1; then
        MODE="orca"
        echo "Auto-detected mode: ORCA (*.inp files found)"
    else
        echo "ERROR: No valid input files found in $ORIGIN_DIR"
        echo "       Expected: CONTROL.txt + input.txt (DELFIN) or *.inp (ORCA)"
        exit 1
    fi
fi
export DELFIN_MODE="$MODE"

# ======================================================================
# Section 14: Workspace copy to scratch
# ======================================================================
copy_workspace_to_scratch() {
    local copy_profile="${1:-full}"
    local manifest_file="$DELFIN_SCRATCH/.fresh_inputs.txt"
    local staged_count=0
    echo "Copying files to scratch (profile: ${copy_profile})..."

    if [ "$copy_profile" = "fresh" ]; then
        : > "$manifest_file"
        if [ -f "$ORIGIN_DIR/CONTROL.txt" ]; then
            printf '%s\n' "CONTROL.txt" >> "$manifest_file"
            input_entry="$(awk -F= '
                /^[[:space:]]*input_file[[:space:]]*=/ {
                    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
                    print $2
                    exit
                }
            ' "$ORIGIN_DIR/CONTROL.txt")"
            input_entry="${input_entry:-input.txt}"
            if [[ "$input_entry" = /* ]]; then
                echo "Fresh staging fallback: CONTROL.txt references absolute input_file '$input_entry'."
                staged_count=0
            else
                printf '%s\n' "$input_entry" >> "$manifest_file"
            fi
        fi
        for optional_file in input.txt input.xyz start.txt co2.xyz; do
            [ -e "$ORIGIN_DIR/$optional_file" ] && printf '%s\n' "$optional_file" >> "$manifest_file"
        done
        awk 'NF && !seen[$0]++' "$manifest_file" > "${manifest_file}.tmp" && mv "${manifest_file}.tmp" "$manifest_file"

        if [ -s "$manifest_file" ]; then
            while IFS= read -r rel_path; do
                [ -e "$ORIGIN_DIR/$rel_path" ] || continue
                mkdir -p "$RUN_DIR/$(dirname "$rel_path")"
                cp -a "$ORIGIN_DIR/$rel_path" "$RUN_DIR/$rel_path"
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
                --exclude='.git/' --exclude='.venv/' --exclude='__pycache__/' \
                --exclude='.ipynb_checkpoints/' --exclude='.orca_iso*' \
                --exclude='delfin_*.out' --exclude='delfin_*.err' \
                --exclude='*.out' --exclude='*.gbw' --exclude='*.hess' \
                --exclude='*.cube' --exclude='*.densitiesinfo' \
                --exclude='*.property.txt' --exclude='*.png' --exclude='*.svg' \
                --exclude='*.pdf' --exclude='*.docx' --exclude='*.json' \
                --exclude='*.fprint' --exclude='*.tmp' \
                --exclude='*_OCCUPIER/' --exclude='ESD/' --exclude='builder/' \
                "$ORIGIN_DIR"/ "$RUN_DIR"/ 2>/dev/null || true
        fi
    else
        rsync -a \
            --exclude='.git/' --exclude='.venv/' --exclude='__pycache__/' \
            --exclude='.ipynb_checkpoints/' --exclude='.orca_iso*' \
            "$ORIGIN_DIR"/ "$RUN_DIR"/ 2>/dev/null || true
    fi
}

# Fresh DELFIN/BUILD/GUPPY/CO2 runs only need current inputs; recalc and raw ORCA
# runs keep the full workspace copy because they may depend on existing outputs.
COPY_PROFILE="full"
case "$MODE" in
    delfin|build|guppy|hyperpol_xtb|tadf_xtb|censo_anmr|delfin-co2-chain)
        COPY_PROFILE="fresh"
        ;;
esac

if [ "$STAGE_IO" = "1" ] && [ -d "$RUN_DIR" ]; then
    copy_workspace_to_scratch "$COPY_PROFILE"
    rm -f "$RUN_DIR"/delfin_*.out "$RUN_DIR"/delfin_*.err 2>/dev/null || true
else
    # No staging — run directly in origin dir
    RUN_DIR="$ORIGIN_DIR"
fi

cd "$RUN_DIR"

# ======================================================================
# Section 15: Result sync helpers
# ======================================================================
rescue_iso_files() {
    [ -d "$RUN_DIR" ] || return 0
    while IFS= read -r -d '' iso_dir; do
        parent_dir="$(dirname "$iso_dir")"
        find "$iso_dir" -maxdepth 1 -type f | while IFS= read -r f; do
            fname="$(basename "$f")"
            if [[ "$fname" == *.inp ]]; then
                cp -n "$f" "$parent_dir/" 2>/dev/null || true
            else
                cp -f "$f" "$parent_dir/" 2>/dev/null || true
            fi
        done
    done < <(find "$RUN_DIR" -name '.orca_iso*' -type d -print0 2>/dev/null)
}

purge_xtb_goat_out_files() {
    [ -d "$RUN_DIR" ] || return 0
    find "$RUN_DIR" -type f -name 'XTB_GOAT.goat.*.out' -delete 2>/dev/null || true
}

SYNC_STAMP_MINIMAL="$DELFIN_SCRATCH/.last_result_sync.minimal"

collect_sync_results() {
    local list_file="$1"
    local profile="${2:-full}"
    local stamp_path=""

    # Disable pipefail inside this function: find may encounter vanishing
    # files (ORCA MPI cleanup) and return non-zero, which with set -eo
    # pipefail would abort the entire script before the final sync runs.
    set +o pipefail 2>/dev/null || true

    case "$profile" in
        minimal)
            stamp_path="$SYNC_STAMP_MINIMAL"
            local find_filter=(
                \( -type f -o -type l \)
                ! -path '*/.orca_iso*/*'
                ! -path '*/__pycache__/*'
                ! -path '*/.ipynb_checkpoints/*'
                ! -name '*.tmp' ! -name '*.pyc' ! -name '*.pyo'
                ! -name 'XTB_GOAT.goat.*.out'
                \( -name '*.out' -o -name '*.err' -o -name '*.log' -o -name '*.json'
                   -o -name '*.txt' -o -name '*.xyz' -o -name '*.inp' -o -name '*.png'
                   -o -name '*.svg' -o -name '*.pdf' -o -name '*.docx' -o -name '*.csv'
                   -o -name '*.dat' \)
            )
            if [ -f "$stamp_path" ]; then
                find "$RUN_DIR" "${find_filter[@]}" -newer "$stamp_path" -printf '%P\n' 2>/dev/null | sort -u > "$list_file" || true
            else
                find "$RUN_DIR" "${find_filter[@]}" -printf '%P\n' 2>/dev/null | sort -u > "$list_file" || true
            fi
            ;;
        *)
            find "$RUN_DIR" \
                \( -type f -o -type l \) \
                ! -path '*/.orca_iso*/*' \
                ! -path '*/__pycache__/*' \
                ! -path '*/.ipynb_checkpoints/*' \
                ! -name '*.tmp' ! -name '*.pyc' ! -name '*.pyo' \
                ! -name 'XTB_GOAT.goat.*.out' \
                -printf '%P\n' 2>/dev/null | sort -u > "$list_file" || true
            ;;
    esac

    set -o pipefail 2>/dev/null || true
}

sync_results_back() {
    local sync_label="${1:-sync}"
    local sync_profile="${2:-full}"
    local list_file="$DELFIN_SCRATCH/.sync_files.txt"
    local file_count

    if [ "$RUN_DIR" = "$ORIGIN_DIR" ]; then return 0; fi
    if [ ! -d "$RUN_DIR" ]; then
        echo "WARNING: RUN_DIR not found, nothing to copy."
        return 0
    fi

    purge_xtb_goat_out_files
    [ "$sync_profile" = "full" ] && rescue_iso_files
    collect_sync_results "$list_file" "$sync_profile"

    if [ ! -s "$list_file" ]; then
        echo "No updated result files to copy (${sync_label})."
        [ "$sync_profile" = "minimal" ] && touch "$SYNC_STAMP_MINIMAL" 2>/dev/null || true
        return 0
    fi

    file_count="$(wc -l < "$list_file" 2>/dev/null || echo 0)"
    rsync -a --files-from="$list_file" "$RUN_DIR"/ "$ORIGIN_DIR"/ 2>/dev/null || true
    [ "$sync_profile" = "minimal" ] && touch "$SYNC_STAMP_MINIMAL" 2>/dev/null || true
    echo "Copied ${file_count} result files (${sync_label}, profile=${sync_profile})."
}

# ======================================================================
# Section 16: Periodic sync + signal handling
# ======================================================================
_DELFIN_PID=""
PERIODIC_COPY_PID=""
FINAL_BACKUP_PID=""

periodic_copy() {
    sleep 60
    [ -d "$RUN_DIR" ] && sync_results_back "initial-60s" "minimal"
    while true; do
        sleep "${SYNC_INTERVAL}"
        [ -d "$RUN_DIR" ] && sync_results_back "periodic" "minimal"
    done
}

periodic_copy &
PERIODIC_COPY_PID=$!

get_walltime_seconds() {
    local timelimit
    timelimit=$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | grep -oP 'TimeLimit=\K[^ ]+' || echo "")
    if [ -z "$timelimit" ] || [ "$timelimit" = "UNLIMITED" ]; then
        echo "0"
        return
    fi
    echo "$timelimit" | awk -F'[-:]' '{
        if (NF==4) print $1*86400 + $2*3600 + $3*60 + $4
        else if (NF==3) print $1*3600 + $2*60 + $3
        else if (NF==2) print $1*60 + $2
        else print 0
    }'
}

SAFETY_MARGIN=600
schedule_final_backup() {
    local walltime_sec
    walltime_sec=$(get_walltime_seconds)
    if [ "$walltime_sec" -le "$SAFETY_MARGIN" ]; then
        return
    fi
    local wait_time=$((walltime_sec - SAFETY_MARGIN))
    echo "Scheduled final backup in $((wait_time / 3600))h $((wait_time % 3600 / 60))m (10 min before timeout)"
    sleep "$wait_time"
    echo ""
    echo "========================================"
    echo "Final backup before timeout: $(date)"
    echo "========================================"
    [ -d "$RUN_DIR" ] && sync_results_back "final-backup" "full"
    echo "Final backup completed."
}

schedule_final_backup &
FINAL_BACKUP_PID=$!

stop_background_tasks() {
    [ -n "${PERIODIC_COPY_PID:-}" ] && kill "$PERIODIC_COPY_PID" 2>/dev/null || true
    [ -n "${FINAL_BACKUP_PID:-}" ] && kill "$FINAL_BACKUP_PID" 2>/dev/null || true
}

cleanup() {
    set +e
    local signal_name="${1:-UNKNOWN}"
    echo ""
    echo "========================================"
    echo "Caught $signal_name signal, running cleanup..."
    echo "========================================"

    # Stop DELFIN process if running
    if [ -n "$_DELFIN_PID" ]; then
        kill -TERM -- -"$_DELFIN_PID" 2>/dev/null || kill -TERM "$_DELFIN_PID" 2>/dev/null || true
        wait "$_DELFIN_PID" 2>/dev/null || true
    fi

    cd "$RUN_DIR" 2>/dev/null && delfin --cleanup 2>/dev/null || true
    write_slurm_reason

    echo "Copying results back to $ORIGIN_DIR..."
    sync_results_back "signal-${signal_name}" "full"

    rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR" "${VENV_LOCAL:-}" "${ORCA_LOCAL:-}" "${OMPI_LOCAL:-}" 2>/dev/null || true
    echo "Cleanup completed at $(date)"
    exit 1
}

handle_early_warning() {
    set +e
    echo ""
    echo "[warning] Approaching walltime -- performing preemptive sync..."
    [ -d "$RUN_DIR" ] && sync_results_back "preemptive-walltime" "full"
    echo "[warning] Preemptive sync done."
}

trap 'stop_background_tasks' EXIT
trap 'cleanup SIGTERM' SIGTERM
trap 'cleanup SIGINT' SIGINT
trap 'cleanup SIGHUP' SIGHUP
trap 'handle_early_warning' SIGUSR1

# ======================================================================
# Section 17: xtb4stda home configuration (optional)
# ======================================================================
configure_local_xtb4stda_home() {
    local stage_base="${STAGE_BASE:-}"
    if [ -z "$stage_base" ] || [ ! -d "$stage_base" ]; then return 0; fi
    local short_link="$stage_base/xtb4stda_${SLURM_JOB_ID:-0}"
    local target=""
    target="$("$PYTHON_BIN" - <<'PY' 2>/dev/null
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

# ======================================================================
# Section 18: NMR post-processing fallback
# ======================================================================
_postprocess_nmr_if_needed() {
    if [ "${DELFIN_MODE:-}" != "orca" ]; then return; fi
    local inp_file="${DELFIN_INP_FILE:-}"
    if [ -z "$inp_file" ]; then
        inp_file="$(find . -maxdepth 1 -type f -name '*.inp' | sort | head -n 1)"
        inp_file="${inp_file#./}"
    fi
    [ -z "$inp_file" ] || [ ! -f "$inp_file" ] && return
    grep -qi 'EPRNMR' "$inp_file" || return

    local stem="${inp_file%.inp}"
    local out_file="${stem}.out"
    local png_file="NMR_${stem}.png"
    [ ! -f "$out_file" ] && return
    [ -f "$png_file" ] && return

    echo "[nmr] Missing ${png_file}; running fallback NMR report generation"
    "$PYTHON_BIN" -m delfin.cli_nmr_report "$out_file" --table || \
        echo "[nmr] WARNING: fallback NMR report generation failed"
}

# ======================================================================
# Section 19: Job banner
# ======================================================================
echo "========================================"
echo "DELFIN Job - $DISPLAY_JOB_NAME"
echo "Mode: $MODE"
echo "========================================"
echo "Job ID:      ${SLURM_JOB_ID:-$DELFIN_JOB_ID}"
echo "Node:        ${SLURM_JOB_NODELIST:-$(hostname)}"
echo "CPUs:        ${SLURM_NTASKS:-$(nproc)}"
echo "Memory:      ${SLURM_MEM_PER_NODE:-unknown} MB"
TIME_LIMIT="${SLURM_TIMELIMIT:-}"
if [ -z "$TIME_LIMIT" ] && [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
    TIME_LIMIT="$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | grep -oP 'TimeLimit=\K[^ ]+' || echo "")"
fi
echo "Time Limit:  ${TIME_LIMIT:-unknown}"
echo "Submit Dir:  $ORIGIN_DIR"
echo "Scratch:     $DELFIN_SCRATCH"
echo "Started:     $(date)"
[ -n "$INP_FILE" ] && echo "Input File:  $INP_FILE"
echo "========================================"
echo ""

if command -v mpirun >/dev/null 2>&1; then
    echo "OpenMPI: $(which mpirun)"
    mpirun --version 2>&1 | head -1 || echo "MPI check failed"
    echo ""
fi

orca_path=$(command -v orca || true)
echo "ORCA: ${orca_path:-not found}"
[ -n "$orca_path" ] && { orca --version 2>&1 | head -5 || true; }
echo ""

_DELFIN_VER="$(delfin --version 2>&1 || echo 'unknown')"
[ -n "${DELFIN_RUNTIME_KEY:-}" ] && _DELFIN_VER="$_DELFIN_VER ($DELFIN_RUNTIME_KEY)"
echo "DELFIN Version: $_DELFIN_VER"
echo ""

purge_xtb_goat_out_files
touch "$SYNC_STAMP_MINIMAL" 2>/dev/null || true

# ======================================================================
# Section 20: Run DELFIN (always via local_runner.py)
# ======================================================================
_JOB_START_EPOCH=$(date +%s)

"$PYTHON_BIN" -m delfin.dashboard.local_runner &
_DELFIN_PID=$!
wait "$_DELFIN_PID"
EXIT_CODE=$?

_JOB_END_EPOCH=$(date +%s)
_JOB_ELAPSED=$(( _JOB_END_EPOCH - _JOB_START_EPOCH ))

echo ""
echo "========================================"
echo "Job finished: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "Elapsed:     $(printf '%02d:%02d:%02d' $((_JOB_ELAPSED/3600)) $(((_JOB_ELAPSED%3600)/60)) $((_JOB_ELAPSED%60)))"
echo "========================================"

if [ "$EXIT_CODE" -eq 0 ]; then
    _postprocess_nmr_if_needed
fi

# Final result sync
sync_results_back "job-end" "full"

# Cleanup scratch and staged binaries
rm -rf "$DELFIN_SCRATCH" "$ORCA_TMPDIR" "${VENV_LOCAL:-}" "${ORCA_LOCAL:-}" "${OMPI_LOCAL:-}" 2>/dev/null || true

exit $EXIT_CODE
