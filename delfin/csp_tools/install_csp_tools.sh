#!/usr/bin/env bash

set -euo pipefail

# ---------------------------------------------------------------------------
# CSP (Crystal Structure Prediction) tools installer for DELFIN.
#
# Installs Genarris (https://github.com/Yi5817/Genarris) from source.
#
# IMPORTANT lessons learned during integration:
#
# 1. MPI consistency: Genarris C extensions (cgenarris, rigid_press) must be
#    compiled against the SAME MPI implementation that mpi4py uses at runtime.
#    Mixing OpenMPI headers with Intel MPI libs (or vice versa) causes
#    "undefined symbol: ompi_mpi_char" at import time.
#
# 2. SWIG + clean build: The SWIG-generated C wrapper files and .so binaries
#    must be deleted before rebuilding, otherwise pip reuses cached artifacts
#    that were linked against the wrong libmpi.
#
# 3. numpy version: Genarris declares numpy>=2.0 but works fine with
#    numpy 1.26.x at runtime.  DELFIN requires numpy<2.  We install with
#    --no-build-isolation so that the build uses the environment's numpy,
#    and pass --no-deps to avoid pip pulling numpy 2.x as a dependency.
#    Pre-requisite packages (mpi4py, swig, Cython, setuptools) must already
#    be present in the target environment.
#
# 4. Conda Intel MPI: If conda/micromamba installed mpi4py, it often pulls
#    Intel MPI (impi_rt).  The system mpicc is usually OpenMPI.  Solution:
#    remove conda's impi_rt, reinstall mpi4py via pip (which picks up
#    system OpenMPI), then build Genarris.
#
# Environment variables:
#   GENARRIS_REPO          Git clone URL   (default: GitHub Yi5817/Genarris)
#   GENARRIS_BRANCH        Branch to clone (default: main)
#   FORCE_REINSTALL        Set to 1 to reinstall even if already present
#   MPICC                  MPI C compiler  (default: auto-detect)
#   DELFIN_CSP_TOOLS_ROOT  Override install root
# ---------------------------------------------------------------------------

ROOT="${DELFIN_CSP_TOOLS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
BIN_DIR="${ROOT}/bin"
BUILD_DIR="${ROOT}/.build"
LOG_DIR="${ROOT}/logs"
GENARRIS_REPO="${GENARRIS_REPO:-https://github.com/Yi5817/Genarris.git}"
GENARRIS_BRANCH="${GENARRIS_BRANCH:-main}"
FORCE_REINSTALL="${FORCE_REINSTALL:-0}"
MPICC="${MPICC:-}"

log() {
  printf "[csp_tools] %s\n" "$*"
}

warn() {
  printf "[csp_tools] WARNING: %s\n" "$*" >&2
}

die() {
  printf "[csp_tools] ERROR: %s\n" "$*" >&2
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
# Detect the right MPI C compiler.
# Prefers MPICC env var > system mpicc (OpenMPI).
# Warns if the detected mpicc is from conda (likely Intel MPI mismatch).
# ---------------------------------------------------------------------------
detect_mpicc() {
  if [ -n "${MPICC}" ] && have "${MPICC}"; then
    printf "%s\n" "$(command -v "${MPICC}")"
    return 0
  fi

  # Prefer system mpicc over conda mpicc
  local candidates=("/usr/bin/mpicc" "/usr/local/bin/mpicc")
  for candidate in "${candidates[@]}"; do
    if [ -x "${candidate}" ]; then
      printf "%s\n" "${candidate}"
      return 0
    fi
  done

  # Fall back to whatever is in PATH
  if have mpicc; then
    local found
    found="$(command -v mpicc)"
    # Warn if it looks like a conda environment mpicc (often Intel MPI)
    if [[ "${found}" == *"/envs/"* ]] || [[ "${found}" == *"/micromamba/"* ]] || [[ "${found}" == *"/conda/"* ]]; then
      warn "mpicc found at ${found} — this may be Intel MPI from conda."
      warn "If import fails with 'undefined symbol: ompi_mpi_char', set MPICC=/usr/bin/mpicc"
    fi
    printf "%s\n" "${found}"
    return 0
  fi

  return 1
}

# ---------------------------------------------------------------------------
# Check all system-level prerequisites.
# ---------------------------------------------------------------------------
check_system_deps() {
  log "checking system dependencies..."

  local missing=()
  local python_bin

  if ! have git; then
    missing+=("git")
  fi

  if ! detect_mpicc >/dev/null 2>&1; then
    warn "mpicc not found. Genarris needs an MPI C compiler."
    warn "Install via:  sudo apt install libopenmpi-dev"
    warn "         or:  conda install -c conda-forge openmpi"
    missing+=("mpicc")
  fi

  if ! have swig; then
    warn "SWIG not found. Genarris needs SWIG for C extensions."
    warn "Install via:  sudo apt install swig"
    warn "         or:  conda install -c conda-forge swig"
    missing+=("swig")
  fi

  python_bin="$(detect_python)" || { missing+=("python"); true; }
  if [ -n "${python_bin:-}" ]; then
    if ! python_has_module "${python_bin}" "mpi4py"; then
      warn "mpi4py not found in ${python_bin}."
      warn "Install via:  MPICC=/usr/bin/mpicc pip install mpi4py"
      warn "         or:  conda install -c conda-forge mpi4py"
      missing+=("mpi4py")
    fi
    if ! python_has_module "${python_bin}" "Cython"; then
      warn "Cython not found. Installing via pip..."
      "${python_bin}" -m pip install --quiet Cython || missing+=("Cython")
    fi
  fi

  if [ ${#missing[@]} -gt 0 ]; then
    die "missing dependencies: ${missing[*]}"
  fi

  log "system dependencies OK"
}

# ---------------------------------------------------------------------------
# Verify that mpi4py and the detected mpicc use the same MPI implementation.
# ---------------------------------------------------------------------------
verify_mpi_consistency() {
  local python_bin="$1"
  local mpicc_bin="$2"

  log "verifying MPI consistency..."

  # Check what libmpi mpi4py links against
  local mpi4py_libmpi
  mpi4py_libmpi="$("${python_bin}" -c "
import mpi4py.MPI, ctypes.util, pathlib, re, subprocess as sp
so = pathlib.Path(mpi4py.MPI.__file__)
out = sp.run(['ldd', str(so)], capture_output=True, text=True).stdout
for line in out.splitlines():
    if 'libmpi' in line:
        print(line.strip().split()[2] if '=>' in line else line.strip())
        break
" 2>/dev/null || true)"

  if [ -z "${mpi4py_libmpi}" ]; then
    warn "could not determine mpi4py's libmpi — skipping consistency check"
    return 0
  fi

  # Check if mpicc is OpenMPI
  local mpicc_flavor="unknown"
  if "${mpicc_bin}" -showme 2>/dev/null | grep -q "openmpi\|ompi"; then
    mpicc_flavor="openmpi"
  elif "${mpicc_bin}" --show 2>/dev/null | grep -q "impi\|intel"; then
    mpicc_flavor="intelmpi"
  fi

  # Check if mpi4py's libmpi is OpenMPI (has ompi_ symbols)
  local mpi4py_flavor="unknown"
  if nm -D "${mpi4py_libmpi}" 2>/dev/null | grep -q "ompi_mpi_"; then
    mpi4py_flavor="openmpi"
  elif nm -D "${mpi4py_libmpi}" 2>/dev/null | grep -q "PMPI_"; then
    # Could be either; if no ompi_ symbols, probably not OpenMPI
    if ! nm -D "${mpi4py_libmpi}" 2>/dev/null | grep -q "ompi_"; then
      mpi4py_flavor="non-openmpi"
    fi
  fi

  log "mpi4py links: ${mpi4py_libmpi}"
  log "mpicc flavor: ${mpicc_flavor}, mpi4py flavor: ${mpi4py_flavor}"

  if [ "${mpicc_flavor}" = "openmpi" ] && [ "${mpi4py_flavor}" = "non-openmpi" ]; then
    die "MPI mismatch: mpicc is OpenMPI but mpi4py links against non-OpenMPI (${mpi4py_libmpi}).
Remove conda Intel MPI and reinstall mpi4py:
  conda remove impi_rt mpi
  MPICC=${mpicc_bin} pip install --no-binary mpi4py mpi4py"
  fi

  log "MPI consistency OK"
}

# ---------------------------------------------------------------------------
# Clone (or update) and install Genarris.
# ---------------------------------------------------------------------------
install_genarris() {
  local python_bin mpicc_bin
  python_bin="$(detect_python)" || die "python/python3 not found"
  mpicc_bin="$(detect_mpicc)" || die "mpicc not found"

  if python_has_module "${python_bin}" "gnrs" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "Genarris (gnrs) already installed, skipping (set FORCE_REINSTALL=1 to force)"
    link_gnrs_cli "${python_bin}"
    return 0
  fi

  verify_mpi_consistency "${python_bin}" "${mpicc_bin}"

  # -- clone / update repo ------------------------------------------------
  log "cloning Genarris from ${GENARRIS_REPO} (branch: ${GENARRIS_BRANCH})..."
  mkdir -p "${BUILD_DIR}" "${LOG_DIR}"
  local genarris_dir="${BUILD_DIR}/Genarris"

  if [ -d "${genarris_dir}" ]; then
    log "updating existing clone..."
    cd "${genarris_dir}"
    git fetch origin
    git checkout "${GENARRIS_BRANCH}"
    git pull origin "${GENARRIS_BRANCH}"
  else
    git clone --branch "${GENARRIS_BRANCH}" "${GENARRIS_REPO}" "${genarris_dir}"
    cd "${genarris_dir}"
  fi

  log "initializing submodules (cgenarris)..."
  git submodule update --init --recursive

  # -- clean previous build artifacts (critical for MPI relinking) --------
  log "cleaning previous build artifacts..."
  rm -rf build/ dist/ *.egg-info
  find gnrs/ -name "*.so" -delete 2>/dev/null || true
  find gnrs/ -name "*_wrap.c" -delete 2>/dev/null || true

  # -- install with correct MPI ------------------------------------------
  log "installing Genarris with MPICC=${mpicc_bin}..."
  log "(using --no-build-isolation --no-deps to keep numpy<2 intact)"

  export CC="${mpicc_bin}"
  export MPICC="${mpicc_bin}"

  "${python_bin}" -m pip install . \
    --no-build-isolation \
    --no-deps \
    --no-cache-dir \
    --force-reinstall \
    2>&1 | tee "${LOG_DIR}/genarris_install.log"

  if python_has_module "${python_bin}" "gnrs"; then
    log "Genarris installed successfully"
  else
    die "Genarris installation failed — 'import gnrs' not possible. Check ${LOG_DIR}/genarris_install.log"
  fi

  # -- verify C extension actually loads ----------------------------------
  log "verifying C extension..."
  if "${python_bin}" -c "from gnrs.cgenarris import pygenarris_mpi" 2>/dev/null; then
    log "cgenarris C extension OK"
  else
    warn "cgenarris C extension failed to load!"
    warn "This usually means MPI mismatch. Check:"
    warn "  ldd \$(python -c 'import gnrs.cgenarris.src._pygenarris_mpi as m; print(m.__file__)')"
    warn "  ldd \$(python -c 'import mpi4py.MPI; print(mpi4py.MPI.__file__)')"
    warn "Both must link against the same libmpi.so version."
    die "C extension verification failed"
  fi

  link_gnrs_cli "${python_bin}"
}

# ---------------------------------------------------------------------------
# Link the gnrs CLI into csp_tools/bin.
# ---------------------------------------------------------------------------
link_gnrs_cli() {
  local python_bin="$1"
  local gnrs_bin

  gnrs_bin="$(command -v gnrs 2>/dev/null || true)"
  if [ -z "${gnrs_bin}" ]; then
    # Try to find it relative to the python binary
    local python_dir
    python_dir="$(dirname "${python_bin}")"
    if [ -x "${python_dir}/gnrs" ]; then
      gnrs_bin="${python_dir}/gnrs"
    fi
  fi

  if [ -n "${gnrs_bin}" ]; then
    mkdir -p "${BIN_DIR}"
    ln -sfn "${gnrs_bin}" "${BIN_DIR}/gnrs"
    log "linked gnrs -> ${gnrs_bin}"
  else
    warn "gnrs CLI not found in PATH after install"
  fi
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
summary() {
  local python_bin

  log "============================================"
  log "  CSP tools installation summary"
  log "============================================"

  if python_bin="$(detect_python)"; then
    if python_has_module "${python_bin}" "gnrs"; then
      local version
      version="$("${python_bin}" -c "import gnrs; print(gnrs.__version__)" 2>/dev/null || echo "unknown")"
      printf "  %-16s %s\n" "Genarris (gnrs)" "v${version}"
    else
      printf "  %-16s %s\n" "Genarris (gnrs)" "MISSING"
    fi

    if "${python_bin}" -c "from gnrs.cgenarris import pygenarris_mpi" 2>/dev/null; then
      printf "  %-16s %s\n" "cgenarris ext" "OK"
    else
      printf "  %-16s %s\n" "cgenarris ext" "FAILED (MPI mismatch?)"
    fi

    if python_has_module "${python_bin}" "mpi4py"; then
      printf "  %-16s %s\n" "mpi4py" "installed"
    else
      printf "  %-16s %s\n" "mpi4py" "MISSING"
    fi
  fi

  if [ -x "${BIN_DIR}/gnrs" ]; then
    printf "  %-16s %s\n" "gnrs CLI" "${BIN_DIR}/gnrs"
  else
    local gnrs_sys
    gnrs_sys="$(command -v gnrs 2>/dev/null || true)"
    if [ -n "${gnrs_sys}" ]; then
      printf "  %-16s %s\n" "gnrs CLI" "${gnrs_sys} (system)"
    else
      printf "  %-16s %s\n" "gnrs CLI" "not found"
    fi
  fi

  printf "  %-16s %s\n" "csp_tools root" "${ROOT}"
  printf "\n"
  printf "Activate with:\n"
  printf "  source %s/env.sh\n" "${ROOT}"
  printf "\n"
  printf "If the C extension fails to load, check MPI consistency:\n"
  printf "  python -c \"from gnrs.cgenarris import pygenarris_mpi; print('OK')\"\n"
}

main() {
  check_system_deps
  install_genarris
  summary
}

main "$@"
