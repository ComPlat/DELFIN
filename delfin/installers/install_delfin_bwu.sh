#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# DELFIN install script for BwUniCluster 3.0
# Packaged copy used by the dashboard for PyPI/wheel installs.
# ============================================================================

DELFIN_REPO="${DELFIN_REPO:-$HOME/software/delfin}"
DELFIN_REPO_URL="${DELFIN_REPO_URL:-https://github.com/ComPlat/DELFIN.git}"
ORCA_TARBALL="${ORCA_TARBALL:-$HOME/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2.tar.xz}"
ORCA_DIR="${ORCA_DIR:-$HOME/software/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2}"
OPENMPI_PREFIX="${OPENMPI_PREFIX:-$HOME/software/openmpi-4.1.8}"
OPENMPI_URL="${OPENMPI_URL:-https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.8.tar.gz}"
MAKE_JOBS="${MAKE_JOBS:-8}"
DELFIN_AUTO_ACTIVATE_VENV="${DELFIN_AUTO_ACTIVATE_VENV:-1}"
REQUIRE_ORCA_PLOT="${REQUIRE_ORCA_PLOT:-1}"
DELFIN_CALC_DIR="${DELFIN_CALC_DIR:-$HOME/calc}"
DELFIN_ARCHIVE_DIR="${DELFIN_ARCHIVE_DIR:-$HOME/archive}"

log() { printf "[install] %s\n" "$*"; }
err() { printf "[install][error] %s\n" "$*" >&2; }

find_orca_dir() {
  if [ -x "$ORCA_DIR/orca" ]; then
    return 0
  fi
  local detected
  detected="$(
    find "$HOME/software" "$HOME/apps" "$HOME/local" "/opt" \
      -maxdepth 4 -type f -name orca 2>/dev/null | head -n1 || true
  )"
  if [ -n "$detected" ]; then
    ORCA_DIR="$(dirname "$detected")"
    log "Detected ORCA directory: $ORCA_DIR"
  fi
  [ -x "$ORCA_DIR/orca" ]
}

ensure_module_command() {
  if command -v module >/dev/null 2>&1; then
    return 0
  fi
  local init
  for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /etc/profile.d/lmod.sh; do
    if [ -f "$init" ]; then
      # shellcheck disable=SC1090
      source "$init"
      break
    fi
  done
  command -v module >/dev/null 2>&1
}

ensure_delfin_env_sourced() {
  local target="$1"
  local source_line='source "$HOME/.delfin_env.sh"'
  if [ ! -f "$target" ]; then
    touch "$target"
  fi
  if ! grep -Fq "$source_line" "$target" 2>/dev/null; then
    printf "\n# DELFIN environment\n%s\n" "$source_line" >> "$target"
    log "Added source \$HOME/.delfin_env.sh to $target"
  else
    log "$target already sources ~/.delfin_env.sh"
  fi
}

if ! ensure_module_command; then
  err "'module' command not found. Run this on BwUniCluster login node."
  exit 1
fi

log "Step 1/4: OpenMPI 4.1.8"
if [ -x "$OPENMPI_PREFIX/bin/mpirun" ]; then
  log "OpenMPI already installed at $OPENMPI_PREFIX"
else
  module purge
  module load compiler/gnu/14.2

  mkdir -p "$HOME/software"
  cd "$HOME/software"

  if [ ! -f openmpi-4.1.8.tar.gz ]; then
    log "Downloading OpenMPI 4.1.8"
    if command -v wget >/dev/null 2>&1; then
      wget "$OPENMPI_URL"
    else
      curl -L -O "$OPENMPI_URL"
    fi
  fi

  if [ ! -d openmpi-4.1.8 ]; then
    tar xzf openmpi-4.1.8.tar.gz
  fi

  cd openmpi-4.1.8
  ./configure --prefix="$OPENMPI_PREFIX" \
              --enable-mpi-cxx \
              --enable-mca-no-build=fs-gpfs \
              --disable-oshmem
  make -j "$MAKE_JOBS"
  make install

  log "OpenMPI installed: $OPENMPI_PREFIX"
  "$OPENMPI_PREFIX/bin/mpirun" --version | head -1
fi

log "Step 2/4: ORCA 6.1.1"
if [ -x "$ORCA_DIR/orca" ]; then
  log "ORCA already present at $ORCA_DIR"
else
  if [ -f "$ORCA_TARBALL" ]; then
    log "Extracting ORCA from $ORCA_TARBALL into $HOME/software"
    mkdir -p "$HOME/software"
    tar xf "$ORCA_TARBALL" -C "$HOME/software"
  else
    err "ORCA tarball not found. Please download it manually (license required)."
    err "Expected: $ORCA_TARBALL"
    err "After download, re-run this script."
    exit 2
  fi
fi

if ! find_orca_dir; then
  err "ORCA binary not found after extraction."
  err "Set ORCA_DIR correctly or provide a tarball containing the ORCA binaries."
  exit 3
fi

if [ "$REQUIRE_ORCA_PLOT" = "1" ] && [ ! -x "$ORCA_DIR/orca_plot" ]; then
  err "orca_plot not found at $ORCA_DIR/orca_plot"
  err "MO/ESP cube generation needs orca_plot. Install full ORCA package with orca_plot."
  exit 4
fi

log "Step 3/4: Python venv + DELFIN"
module purge
module load devel/python/3.11.7-gnu-14.2

mkdir -p "$HOME/software"
mkdir -p "$DELFIN_CALC_DIR"
mkdir -p "$DELFIN_ARCHIVE_DIR"

if [ ! -d "$DELFIN_REPO" ]; then
  log "Cloning DELFIN into $DELFIN_REPO"
  git clone "$DELFIN_REPO_URL" "$DELFIN_REPO"
else
  log "DELFIN repo already exists at $DELFIN_REPO"
fi

if [ ! -d "$DELFIN_REPO/.venv" ]; then
  log "Creating venv at $DELFIN_REPO/.venv"
  python3 -m venv "$DELFIN_REPO/.venv"
fi

# shellcheck disable=SC1090
source "$DELFIN_REPO/.venv/bin/activate"
python -m pip install --upgrade pip wheel

log "Installing DELFIN (pip install -e .)"
python -m pip install -e "$DELFIN_REPO"

log "DELFIN check"
delfin --version || true

log "Packaging local venv tarball for low-I/O job staging"
rm -f "$DELFIN_REPO/delfin_venv.tar"
tar -cf "$DELFIN_REPO/delfin_venv.tar" -C "$DELFIN_REPO" .venv
mkdir -p "$DELFIN_REPO/.runtime_cache"
log "Created $DELFIN_REPO/delfin_venv.tar"
log "Commit-pinned runtime wheels will be cached on demand in $DELFIN_REPO/.runtime_cache"

log "Writing user-local DELFIN runtime defaults"
python - <<EOF
from pathlib import Path
from delfin.user_settings import load_settings, save_settings

settings = load_settings()
runtime = settings.get("runtime", {}) or {}
runtime.setdefault("local", {})
runtime.setdefault("slurm", {})

if not runtime.get("backend") or runtime.get("backend") == "auto":
    runtime["backend"] = "slurm"
if not runtime.get("orca_base"):
    runtime["orca_base"] = str(Path("$ORCA_DIR").expanduser())
if not runtime.get("qm_tools_root"):
    runtime["qm_tools_root"] = str(Path("$DELFIN_REPO/delfin/qm_tools").expanduser())
if not runtime["local"].get("orca_base"):
    runtime["local"]["orca_base"] = str(Path("$ORCA_DIR").expanduser())
if not runtime["slurm"].get("orca_base"):
    runtime["slurm"]["orca_base"] = str(Path("$ORCA_DIR").expanduser())
if not runtime["slurm"].get("submit_templates_dir"):
    runtime["slurm"]["submit_templates_dir"] = str(
        Path("$DELFIN_REPO/examples/example_Job_Submission_Scripts/BwUniCluster/submit_sh").expanduser()
    )
if not runtime["slurm"].get("profile"):
    runtime["slurm"]["profile"] = "bwunicluster3"

settings["runtime"] = runtime
save_settings(settings)
print(f"Saved runtime defaults to {Path.home() / '.delfin_settings.json'}")
EOF

log "Step 4/4: Shell environment"
ENV_FILE="$HOME/.delfin_env.sh"
cat > "$ENV_FILE" <<EOF
# Auto-generated by install_delfin_bwu.sh
export DELFIN_REPO="$DELFIN_REPO"
export OPENMPI_PREFIX="$OPENMPI_PREFIX"
export ORCA_DIR="$ORCA_DIR"
export ORCA_BASE="$ORCA_DIR"
export ORCA_PLOT="$ORCA_DIR/orca_plot"
export PATH="$OPENMPI_PREFIX/bin:$ORCA_DIR:\$PATH"
export LD_LIBRARY_PATH="$OPENMPI_PREFIX/lib:$ORCA_DIR:\${LD_LIBRARY_PATH:-}"
export DELFIN_AUTO_ACTIVATE_VENV="\${DELFIN_AUTO_ACTIVATE_VENV:-$DELFIN_AUTO_ACTIVATE_VENV}"
if [ "\${DELFIN_AUTO_ACTIVATE_VENV}" = "1" ] \
   && [ -z "\${VIRTUAL_ENV:-}" ] \
   && [ -d "$DELFIN_REPO/.venv" ] \
   && [ -z "\${VSCODE_PID:-}" ] \
   && [ "\${TERM_PROGRAM:-}" != "vscode" ]; then
  source "$DELFIN_REPO/.venv/bin/activate"
fi
EOF

ensure_delfin_env_sourced "$HOME/.bashrc"
ensure_delfin_env_sourced "$HOME/.profile"

log "orca_plot check"
if [ -x "$ORCA_DIR/orca_plot" ]; then
  log "orca_plot available at $ORCA_DIR/orca_plot"
else
  log "orca_plot not found (check ORCA install)"
fi

# Apply env for current shell
# shellcheck disable=SC1090
source "$HOME/.delfin_env.sh"

log ""
log "============================================"
log "  DELFIN installation complete!"
log "============================================"
log ""
log "Directories created:"
log "  ~/calc       - working directory for calculations"
log "  ~/archive    - archive for completed calculations"
log ""
log "Next steps:"
log "  source ~/.bashrc        # reload shell environment"
log "  delfin --version        # verify installation"
log "  delfin-voila --port 9000  # launch dashboard (open URL in browser)"
log ""
log "Notes:"
log "  - ORCA download is manual (license required)."
log "  - Set DELFIN_AUTO_ACTIVATE_VENV=0 to disable auto-activation."
log "  - VS Code terminals are excluded from auto-activation."
