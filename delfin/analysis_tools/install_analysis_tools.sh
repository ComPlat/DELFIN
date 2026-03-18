#!/usr/bin/env bash

set -euo pipefail

# ---------------------------------------------------------------------------
# Analysis tools installer for DELFIN.
#
# Installs:
#   - morfeus-ml      (pip, steric descriptors)
#   - CENSO           (conda/pip, conformer ensemble sorting)
#   - Multiwfn        (binary download, wavefunction analysis)
#
# Environment variables:
#   INSTALL_MORFEUS       Set to 1 to install morfeus-ml    (default: 1)
#   INSTALL_CENSO         Set to 1 to install CENSO         (default: 1)
#   INSTALL_MULTIWFN      Set to 1 to install Multiwfn      (default: 0, large binary)
#   FORCE_REINSTALL       Set to 1 to force reinstall       (default: 0)
#   DELFIN_ANALYSIS_TOOLS_ROOT  Override install root
# ---------------------------------------------------------------------------

ROOT="${DELFIN_ANALYSIS_TOOLS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="${ROOT}/logs"

INSTALL_MORFEUS="${INSTALL_MORFEUS:-1}"
INSTALL_CENSO="${INSTALL_CENSO:-1}"
INSTALL_MULTIWFN="${INSTALL_MULTIWFN:-0}"
FORCE_REINSTALL="${FORCE_REINSTALL:-0}"

log() {
  printf "[analysis_tools] %s\n" "$*"
}

warn() {
  printf "[analysis_tools] WARNING: %s\n" "$*" >&2
}

die() {
  printf "[analysis_tools] ERROR: %s\n" "$*" >&2
  exit 1
}

have() {
  command -v "$1" >/dev/null 2>&1
}

detect_python() {
  if have python; then
    command -v python
    return 0
  fi
  if have python3; then
    command -v python3
    return 0
  fi
  return 1
}

python_has_module() {
  local python_bin="$1"
  local module="$2"
  "${python_bin}" -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('${module}') else 1)" >/dev/null 2>&1
}

# ---------------------------------------------------------------------------
install_morfeus() {
  local python_bin="$1"

  if [ "${INSTALL_MORFEUS}" != "1" ]; then
    log "morfeus-ml: skipped (INSTALL_MORFEUS=0)"
    return 0
  fi

  if python_has_module "${python_bin}" "morfeus" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "morfeus-ml: already installed"
    return 0
  fi

  log "installing morfeus-ml..."
  "${python_bin}" -m pip install morfeus-ml 2>&1 | tee -a "${LOG_DIR}/morfeus_install.log"

  if python_has_module "${python_bin}" "morfeus"; then
    local version
    version="$("${python_bin}" -c "from importlib.metadata import version; print(version('morfeus-ml'))" 2>/dev/null || echo "?")"
    log "morfeus-ml v${version} installed successfully"
  else
    warn "morfeus-ml installation failed"
  fi
}

# ---------------------------------------------------------------------------
install_censo() {
  local python_bin="$1"

  if [ "${INSTALL_CENSO}" != "1" ]; then
    log "CENSO: skipped (INSTALL_CENSO=0)"
    return 0
  fi

  if python_has_module "${python_bin}" "censo" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "CENSO: already installed"
    return 0
  fi

  # Try conda first (preferred)
  if have conda || have micromamba || have mamba; then
    local conda_cmd
    if have micromamba; then
      conda_cmd="micromamba"
    elif have mamba; then
      conda_cmd="mamba"
    else
      conda_cmd="conda"
    fi

    log "installing CENSO via ${conda_cmd}..."
    "${conda_cmd}" install -y -c conda-forge censo 2>&1 | tee -a "${LOG_DIR}/censo_install.log" || true

    if python_has_module "${python_bin}" "censo"; then
      log "CENSO installed successfully via ${conda_cmd}"
      return 0
    fi
    log "conda install failed, trying pip from GitHub..."
  fi

  # Fallback: pip from GitHub
  log "installing CENSO via pip from GitHub..."
  "${python_bin}" -m pip install \
    git+https://github.com/grimme-lab/CENSO.git \
    2>&1 | tee -a "${LOG_DIR}/censo_install.log"

  if python_has_module "${python_bin}" "censo"; then
    log "CENSO installed successfully via pip"
  else
    warn "CENSO installation failed — check ${LOG_DIR}/censo_install.log"
  fi
}

# ---------------------------------------------------------------------------
install_multiwfn() {
  local python_bin="$1"

  if [ "${INSTALL_MULTIWFN}" != "1" ]; then
    log "Multiwfn: skipped (INSTALL_MULTIWFN=0, set INSTALL_MULTIWFN=1 to install)"
    log "  Download manually from: http://sobereva.com/multiwfn/"
    log "  Extract and add the directory containing 'Multiwfn' to PATH"
    return 0
  fi

  if have Multiwfn && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "Multiwfn: already available at $(command -v Multiwfn)"
    return 0
  fi

  log "Multiwfn is a standalone binary that must be downloaded manually."
  log "  1. Go to http://sobereva.com/multiwfn/"
  log "  2. Download the Linux version"
  log "  3. Extract and add to PATH:"
  log "     export PATH=/path/to/Multiwfn_directory:\$PATH"
  log "  4. For ORCA integration, edit settings.ini and set orca_2mkl path"
  warn "Automatic download is not supported due to licensing"
}

# ---------------------------------------------------------------------------
summary() {
  local python_bin

  log "============================================"
  log "  Analysis tools installation summary"
  log "============================================"

  if python_bin="$(detect_python)"; then
    for mod_label in "morfeus:morfeus-ml" "censo:CENSO"; do
      local mod="${mod_label%%:*}"
      local label="${mod_label##*:}"
      if python_has_module "${python_bin}" "${mod}"; then
        printf "  %-16s %s\n" "${label}" "installed"
      else
        printf "  %-16s %s\n" "${label}" "not installed"
      fi
    done

    if have Multiwfn; then
      printf "  %-16s %s\n" "Multiwfn" "installed ($(command -v Multiwfn))"
    else
      printf "  %-16s %s\n" "Multiwfn" "not installed (manual download required)"
    fi
  fi

  printf "  %-16s %s\n" "tools root" "${ROOT}"
  printf "\n"
  printf "Usage in Python:\n"
  printf "  from delfin.analysis_tools.morfeus_wrapper import full_steric_report\n"
  printf "  from delfin.analysis_tools.censo_wrapper import run_censo\n"
  printf "  from delfin.analysis_tools.multiwfn_wrapper import bond_order_analysis\n"
}

# ---------------------------------------------------------------------------
main() {
  local python_bin
  python_bin="$(detect_python)" || die "python/python3 not found"

  mkdir -p "${LOG_DIR}"
  install_morfeus "${python_bin}"
  install_censo "${python_bin}"
  install_multiwfn "${python_bin}"
  summary
}

main "$@"
