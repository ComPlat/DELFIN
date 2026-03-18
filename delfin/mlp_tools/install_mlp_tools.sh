#!/usr/bin/env bash

set -euo pipefail

# ---------------------------------------------------------------------------
# MLP (Machine Learning Potential) tools installer for DELFIN.
#
# Installs one or more ML potential backends:
#   - ANI-2x      (torchani)     — H,C,N,O,S,F,Cl
#   - AIMNet2     (aimnet2calc)  — 14 elements + charge support
#   - MACE-OFF    (mace-torch)   — broad element coverage
#   - CHGNet      (chgnet)       — universal potential (Materials Project)
#   - M3GNet      (matgl)        — materials graph neural network
#   - SchNetPack  (schnetpack)   — trainable equivariant GNNs
#   - NequIP      (nequip)       — E(3)-equivariant, data-efficient
#   - ALIGNN      (alignn)       — atomistic line graph NN
#
# Environment variables:
#   INSTALL_ANI2X       Set to 1 to install ANI-2x      (default: 1)
#   INSTALL_AIMNET2     Set to 1 to install AIMNet2      (default: 1)
#   INSTALL_MACE        Set to 1 to install MACE-OFF     (default: 0, large)
#   INSTALL_CHGNET      Set to 1 to install CHGNet       (default: 0)
#   INSTALL_M3GNET      Set to 1 to install M3GNet       (default: 0)
#   INSTALL_SCHNETPACK  Set to 1 to install SchNetPack   (default: 0)
#   INSTALL_NEQUIP      Set to 1 to install NequIP       (default: 0)
#   INSTALL_ALIGNN      Set to 1 to install ALIGNN       (default: 0)
#   FORCE_REINSTALL     Set to 1 to force reinstall      (default: 0)
#   DELFIN_MLP_TOOLS_ROOT  Override install root
# ---------------------------------------------------------------------------

ROOT="${DELFIN_MLP_TOOLS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="${ROOT}/logs"

INSTALL_ANI2X="${INSTALL_ANI2X:-1}"
INSTALL_AIMNET2="${INSTALL_AIMNET2:-1}"
INSTALL_MACE="${INSTALL_MACE:-0}"
INSTALL_CHGNET="${INSTALL_CHGNET:-0}"
INSTALL_M3GNET="${INSTALL_M3GNET:-0}"
INSTALL_SCHNETPACK="${INSTALL_SCHNETPACK:-0}"
INSTALL_NEQUIP="${INSTALL_NEQUIP:-0}"
INSTALL_ALIGNN="${INSTALL_ALIGNN:-0}"
FORCE_REINSTALL="${FORCE_REINSTALL:-0}"

log() {
  printf "[mlp_tools] %s\n" "$*"
}

warn() {
  printf "[mlp_tools] WARNING: %s\n" "$*" >&2
}

die() {
  printf "[mlp_tools] ERROR: %s\n" "$*" >&2
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
check_pytorch() {
  local python_bin="$1"

  if python_has_module "${python_bin}" "torch"; then
    local version
    version="$("${python_bin}" -c "import torch; print(torch.__version__)" 2>/dev/null || echo "unknown")"
    log "PyTorch ${version} found"

    local cuda_avail
    cuda_avail="$("${python_bin}" -c "import torch; print(torch.cuda.is_available())" 2>/dev/null || echo "False")"
    if [ "${cuda_avail}" = "True" ]; then
      log "CUDA available — GPU acceleration enabled"
    else
      log "CUDA not available — using CPU (still fast for small molecules)"
    fi
    return 0
  fi

  warn "PyTorch not found. Installing CPU version..."
  "${python_bin}" -m pip install --quiet torch --index-url https://download.pytorch.org/whl/cpu
}

# ---------------------------------------------------------------------------
install_ani2x() {
  local python_bin="$1"

  if [ "${INSTALL_ANI2X}" != "1" ]; then
    log "ANI-2x: skipped (INSTALL_ANI2X=0)"
    return 0
  fi

  if python_has_module "${python_bin}" "torchani" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "ANI-2x (torchani): already installed"
    return 0
  fi

  log "installing ANI-2x (torchani)..."
  "${python_bin}" -m pip install torchani 2>&1 | tee -a "${LOG_DIR}/ani2x_install.log"

  if python_has_module "${python_bin}" "torchani"; then
    local version
    version="$("${python_bin}" -c "import torchani; print(torchani.__version__)" 2>/dev/null || echo "?")"
    log "ANI-2x v${version} installed successfully"
  else
    warn "ANI-2x installation failed"
  fi
}

# ---------------------------------------------------------------------------
install_aimnet2() {
  local python_bin="$1"

  if [ "${INSTALL_AIMNET2}" != "1" ]; then
    log "AIMNet2: skipped (INSTALL_AIMNET2=0)"
    return 0
  fi

  if python_has_module "${python_bin}" "aimnet2calc" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "AIMNet2 (aimnet2calc): already installed"
    return 0
  fi

  # AIMNet2 requires torch-cluster which needs pre-built wheels
  log "installing torch-cluster (AIMNet2 dependency)..."
  local torch_version
  torch_version="$("${python_bin}" -c "import torch; v=torch.__version__.split('+')[0]; print(v)" 2>/dev/null || echo "")"
  if [ -n "${torch_version}" ]; then
    "${python_bin}" -m pip install torch-cluster \
      -f "https://data.pyg.org/whl/torch-${torch_version}+cpu.html" \
      2>&1 | tee -a "${LOG_DIR}/aimnet2_install.log" || true
  fi

  log "installing AIMNet2 (aimnet2calc from GitHub)..."
  "${python_bin}" -m pip install \
    git+https://github.com/zubatyuk/aimnet2calc.git \
    --no-build-isolation \
    2>&1 | tee -a "${LOG_DIR}/aimnet2_install.log"

  if python_has_module "${python_bin}" "aimnet2calc"; then
    log "AIMNet2 installed successfully"
  else
    warn "AIMNet2 installation failed — check ${LOG_DIR}/aimnet2_install.log"
  fi
}

# ---------------------------------------------------------------------------
install_mace() {
  local python_bin="$1"

  if [ "${INSTALL_MACE}" != "1" ]; then
    log "MACE-OFF: skipped (INSTALL_MACE=0, set INSTALL_MACE=1 to install)"
    return 0
  fi

  if python_has_module "${python_bin}" "mace" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "MACE (mace-torch): already installed"
    return 0
  fi

  log "installing MACE-OFF (mace-torch)... this may take a while"
  "${python_bin}" -m pip install mace-torch 2>&1 | tee -a "${LOG_DIR}/mace_install.log"

  if python_has_module "${python_bin}" "mace"; then
    log "MACE-OFF installed successfully"
  else
    warn "MACE-OFF installation failed"
  fi
}

# ---------------------------------------------------------------------------
install_chgnet() {
  local python_bin="$1"
  if [ "${INSTALL_CHGNET}" != "1" ]; then
    log "CHGNet: skipped (INSTALL_CHGNET=0)"
    return 0
  fi
  if python_has_module "${python_bin}" "chgnet" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "CHGNet: already installed"
    return 0
  fi
  log "installing CHGNet..."
  "${python_bin}" -m pip install chgnet 2>&1 | tee -a "${LOG_DIR}/chgnet_install.log"
  if python_has_module "${python_bin}" "chgnet"; then
    log "CHGNet installed successfully"
  else
    warn "CHGNet installation failed"
  fi
}

# ---------------------------------------------------------------------------
install_m3gnet() {
  local python_bin="$1"
  if [ "${INSTALL_M3GNET}" != "1" ]; then
    log "M3GNet (matgl): skipped (INSTALL_M3GNET=0)"
    return 0
  fi
  if python_has_module "${python_bin}" "matgl" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "M3GNet (matgl): already installed"
    return 0
  fi
  log "installing M3GNet (matgl)..."
  "${python_bin}" -m pip install matgl 2>&1 | tee -a "${LOG_DIR}/m3gnet_install.log"
  if python_has_module "${python_bin}" "matgl"; then
    log "M3GNet installed successfully"
  else
    warn "M3GNet installation failed"
  fi
}

# ---------------------------------------------------------------------------
install_schnetpack() {
  local python_bin="$1"
  if [ "${INSTALL_SCHNETPACK}" != "1" ]; then
    log "SchNetPack: skipped (INSTALL_SCHNETPACK=0)"
    return 0
  fi
  if python_has_module "${python_bin}" "schnetpack" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "SchNetPack: already installed"
    return 0
  fi
  log "installing SchNetPack..."
  "${python_bin}" -m pip install schnetpack 2>&1 | tee -a "${LOG_DIR}/schnetpack_install.log"
  if python_has_module "${python_bin}" "schnetpack"; then
    log "SchNetPack installed successfully"
  else
    warn "SchNetPack installation failed"
  fi
}

# ---------------------------------------------------------------------------
install_nequip() {
  local python_bin="$1"
  if [ "${INSTALL_NEQUIP}" != "1" ]; then
    log "NequIP: skipped (INSTALL_NEQUIP=0)"
    return 0
  fi
  if python_has_module "${python_bin}" "nequip" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "NequIP: already installed"
    return 0
  fi
  log "installing NequIP..."
  "${python_bin}" -m pip install nequip 2>&1 | tee -a "${LOG_DIR}/nequip_install.log"
  if python_has_module "${python_bin}" "nequip"; then
    log "NequIP installed successfully"
  else
    warn "NequIP installation failed"
  fi
}

# ---------------------------------------------------------------------------
install_alignn() {
  local python_bin="$1"
  if [ "${INSTALL_ALIGNN}" != "1" ]; then
    log "ALIGNN: skipped (INSTALL_ALIGNN=0)"
    return 0
  fi
  if python_has_module "${python_bin}" "alignn" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "ALIGNN: already installed"
    return 0
  fi
  log "installing ALIGNN..."
  "${python_bin}" -m pip install alignn 2>&1 | tee -a "${LOG_DIR}/alignn_install.log"
  if python_has_module "${python_bin}" "alignn"; then
    log "ALIGNN installed successfully"
  else
    warn "ALIGNN installation failed"
  fi
}

# ---------------------------------------------------------------------------
verify_backends() {
  local python_bin="$1"

  log "verifying MLP backends..."
  "${python_bin}" -c "
from delfin.mlp_tools import available_backends, get_torchani_version, get_aimnet2_version, get_mace_version
backends = available_backends()
print(f'  Available backends: {backends}')
v = get_torchani_version()
if v: print(f'  ANI-2x:   v{v}')
v = get_aimnet2_version()
if v: print(f'  AIMNet2:  v{v}')
v = get_mace_version()
if v: print(f'  MACE-OFF: v{v}')
if not backends:
    print('  WARNING: No MLP backends available!')
" 2>/dev/null || warn "backend verification failed (delfin.mlp_tools not importable from here)"
}

# ---------------------------------------------------------------------------
summary() {
  local python_bin

  log "============================================"
  log "  MLP tools installation summary"
  log "============================================"

  if python_bin="$(detect_python)"; then
    for mod_name in torchani:ANI-2x aimnet2calc:AIMNet2 mace:MACE-OFF chgnet:CHGNet matgl:M3GNet schnetpack:SchNetPack nequip:NequIP alignn:ALIGNN; do
      local mod="${mod_name%%:*}"
      local label="${mod_name##*:}"
      if python_has_module "${python_bin}" "${mod}"; then
        printf "  %-16s %s\n" "${label}" "installed"
      else
        printf "  %-16s %s\n" "${label}" "not installed"
      fi
    done

    local cuda_avail
    cuda_avail="$("${python_bin}" -c "import torch; print(torch.cuda.is_available())" 2>/dev/null || echo "unknown")"
    printf "  %-16s %s\n" "CUDA" "${cuda_avail}"
  fi

  printf "  %-16s %s\n" "mlp_tools root" "${ROOT}"
  printf "\n"
  printf "Usage in Python:\n"
  printf "  from delfin.mlp_tools.calculators import create_calculator\n"
  printf "  calc = create_calculator('ani2x')  # or 'aimnet2', 'mace_off'\n"
  printf "  atoms.calc = calc\n"
  printf "  energy = atoms.get_potential_energy()\n"
}

# ---------------------------------------------------------------------------
main() {
  local python_bin
  python_bin="$(detect_python)" || die "python/python3 not found"

  mkdir -p "${LOG_DIR}"
  check_pytorch "${python_bin}"
  install_ani2x "${python_bin}"
  install_aimnet2 "${python_bin}"
  install_mace "${python_bin}"
  install_chgnet "${python_bin}"
  install_m3gnet "${python_bin}"
  install_schnetpack "${python_bin}"
  install_nequip "${python_bin}"
  install_alignn "${python_bin}"
  verify_backends "${python_bin}"
  summary
}

main "$@"
