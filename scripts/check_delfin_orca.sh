#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Quick health check for DELFIN + ORCA installation
#
# Checks paths, versions, and Python imports without modifying anything.
# ============================================================================

DELFIN_REPO="${DELFIN_REPO:-$HOME/software/delfin}"
ORCA_DIR="${ORCA_DIR:-$HOME/software/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2}"
OPENMPI_PREFIX="${OPENMPI_PREFIX:-$HOME/software/openmpi-4.1.8}"
DELFIN_ENV_FILE="${DELFIN_ENV_FILE:-$HOME/.delfin_env.sh}"
EXPECTED_ORCA="$ORCA_DIR/orca"
EXPECTED_ORCA_PLOT="$ORCA_DIR/orca_plot"
EXPECTED_MPI="$OPENMPI_PREFIX/bin/mpirun"
EXPECTED_DELFIN_VENV="$DELFIN_REPO/.venv"

log() { printf "[check] %s\n" "$*"; }
warn() { printf "[check][warn] %s\n" "$*" >&2; }
real() { readlink -f "$1" 2>/dev/null || echo "$1"; }

log "Sourcing env"
if [ -f "$DELFIN_ENV_FILE" ]; then
  # shellcheck disable=SC1090
  source "$DELFIN_ENV_FILE"
else
  warn "Env file not found: $DELFIN_ENV_FILE"
fi

VENV_REAL="${VIRTUAL_ENV:-}"
EXPECTED_VENV_REAL="$(real "$EXPECTED_DELFIN_VENV")"

log "VIRTUAL_ENV: ${VIRTUAL_ENV:-<none>}"
log "Python: $(which python)"
log "Delfin: $(which delfin)"
log "OpenMPI: $(which mpirun)"
log "ORCA: $(which orca || echo 'not found')"
log "ORCA_PLOT: $(which orca_plot || echo 'not found')"
log "ORCA_PLOT env: ${ORCA_PLOT:-<unset>}"

if [ -x "$EXPECTED_MPI" ]; then
  if [ "$(which mpirun)" != "$EXPECTED_MPI" ]; then
    log "WARN: mpirun is not the expected one ($EXPECTED_MPI)"
  fi
fi

if [ -x "$EXPECTED_ORCA" ]; then
  if [ "$(which orca)" != "$EXPECTED_ORCA" ]; then
    log "WARN: orca is not the expected one ($EXPECTED_ORCA)"
  fi
else
  log "WARN: expected ORCA binary missing at $EXPECTED_ORCA"
fi

if [ -x "$EXPECTED_ORCA_PLOT" ]; then
  if [ "$(which orca_plot)" != "$EXPECTED_ORCA_PLOT" ]; then
    log "WARN: orca_plot is not the expected one ($EXPECTED_ORCA_PLOT)"
  fi
else
  log "WARN: expected orca_plot missing at $EXPECTED_ORCA_PLOT"
fi

if [ -z "${ORCA_PLOT:-}" ] && [ -x "$EXPECTED_ORCA_PLOT" ]; then
  log "WARN: ORCA_PLOT env is unset (run install/verify script to regenerate ~/.delfin_env.sh)"
fi

if [ -n "${ORCA_PLOT:-}" ] && [ ! -x "$ORCA_PLOT" ]; then
  log "WARN: ORCA_PLOT is set but not executable: $ORCA_PLOT"
fi

if [ -n "$VENV_REAL" ] && [ "$(real "$VENV_REAL")" != "$EXPECTED_VENV_REAL" ]; then
  log "WARN: active venv is not $EXPECTED_DELFIN_VENV"
fi

# Check working directories
log "Directories:"
for d in "$HOME/calc" "$HOME/archive"; do
  if [ -d "$d" ]; then
    log "  $d exists"
  else
    log "  WARN: $d missing"
  fi
done

log "Delfin version"
if command -v delfin >/dev/null 2>&1; then
  delfin --version 2>/dev/null || echo "delfin --version failed"
else
  echo "delfin not found"
fi

log "delfin-voila available"
if command -v delfin-voila >/dev/null 2>&1; then
  echo "yes"
else
  echo "no (run: pip install -e . in $DELFIN_REPO)"
fi

log "OpenMPI version"
if command -v mpirun >/dev/null 2>&1; then
  mpirun --version 2>/dev/null | head -5 || echo "mpirun --version failed"
else
  echo "mpirun not found"
fi

if command -v orca >/dev/null 2>&1; then
  log "ORCA version"
  if orca --version >/dev/null 2>&1; then
    orca --version | head -5
  elif orca -v >/dev/null 2>&1; then
    orca -v | head -5
  else
    orca 2>&1 | head -5 || true
  fi
fi

if [ -x "$EXPECTED_DELFIN_VENV/bin/python" ]; then
  log "Python imports (core packages)"
  "$EXPECTED_DELFIN_VENV/bin/python" - <<'PY'
mods = ("docx", "matplotlib", "pymol", "rdkit", "voila", "mammoth",
        "PIL", "pexpect", "sklearn", "ipywidgets", "py3Dmol")
missing = []
for m in mods:
    try:
        __import__(m)
    except Exception:
        missing.append(m)
if missing:
    print("missing:", ", ".join(missing))
else:
    print("all ok")
PY
else
  log "WARN: venv python not found at $EXPECTED_DELFIN_VENV/bin/python"
fi

log "DONE"
