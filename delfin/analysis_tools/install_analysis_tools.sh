#!/usr/bin/env bash

set -euo pipefail

# ---------------------------------------------------------------------------
# Analysis tools installer for DELFIN.
#
# Installs:
#   - morfeus-ml      (pip, steric descriptors)
#   - CENSO           (conda/pip, conformer ensemble sorting + c2anmr/nmrplot)
#   - ANMR            (ENSO release asset, ensemble NMR spectrum simulation)
#   - Multiwfn        (binary download, wavefunction analysis)
#   - cclib           (pip, log file parser)
#   - nglview         (pip, 3D molecular viewer)
#   - Packmol         (conda, packing tool for MD)
#
# Environment variables:
#   INSTALL_MORFEUS       Set to 1 to install morfeus-ml    (default: 1)
#   INSTALL_CENSO         Set to 1 to install CENSO         (default: 1)
#   INSTALL_ANMR          Set to 1 to install ANMR          (default: 1)
#   INSTALL_MULTIWFN      Set to 1 to install Multiwfn      (default: 0, large binary)
#   INSTALL_CCLIB         Set to 1 to install cclib          (default: 1)
#   INSTALL_NGLVIEW       Set to 1 to install nglview        (default: 1)
#   INSTALL_PACKMOL       Set to 1 to install Packmol        (default: 1)
#   FORCE_REINSTALL       Set to 1 to force reinstall       (default: 0)
#   DELFIN_ANALYSIS_TOOLS_ROOT  Override install root
# ---------------------------------------------------------------------------

ROOT="${DELFIN_ANALYSIS_TOOLS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="${ROOT}/logs"

INSTALL_MORFEUS="${INSTALL_MORFEUS:-1}"
INSTALL_CENSO="${INSTALL_CENSO:-1}"
INSTALL_ANMR="${INSTALL_ANMR:-1}"
INSTALL_MULTIWFN="${INSTALL_MULTIWFN:-0}"
INSTALL_CCLIB="${INSTALL_CCLIB:-1}"
INSTALL_NGLVIEW="${INSTALL_NGLVIEW:-1}"
INSTALL_PACKMOL="${INSTALL_PACKMOL:-1}"
FORCE_REINSTALL="${FORCE_REINSTALL:-0}"
CENSO_PREFER_LATEST="${CENSO_PREFER_LATEST:-1}"
# CENSO from commits, not a branch: pip checks the commit out and git
# verifies it by its id. v3.0.8 (the head of main on 2026-10-16, which is
# what "main" used to install) for Python >= 3.12, which it requires;
# v2.1.4 below that. Move a ref only after reading the new commit.
CENSO_V3_REF="${CENSO_V3_REF:-1a82c5ad5e0b46c17cbbc93ecd7890a3b8c8a1e3}"
CENSO_V2_REF="${CENSO_V2_REF:-db667bde403c916bd20e8a42544fa3b0bed8f7d7}"
# ANMR: the anmr binary of the ENSO release v.2.0.2 (2020-07-30, the newest;
# the asset the former "latest" lookup picked). GitHub records no digest
# for it; the SHA-256 was computed 2026-10-16 from the download. A mirror
# or another build needs both variables.
ENSO_ANMR_URL="${ENSO_ANMR_URL:-https://github.com/grimme-lab/enso/releases/download/v.2.0.2/anmr}"
ENSO_ANMR_SHA256="${ENSO_ANMR_SHA256:-9ad5fe7409eb383b86ed1787b7e8cb43f94437ad80a0311c7f6cafa4932c9b3a}"

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
  # The interpreter DELFIN is actually running in, when the caller says so.
  #
  # Taking the first `python` on the PATH was the default, and on a machine
  # where that is a different environment from the dashboard's -- which it
  # usually is -- every package installed here landed somewhere the dashboard
  # cannot import from. Measured: the installer chose Python 3.13 while the
  # dashboard ran 3.11, so cclib, nglview, censo, morfeus and torch were all
  # installed and all missing at the same time.
  if [[ -n "${DELFIN_PYTHON:-}" && -x "${DELFIN_PYTHON}" ]]; then
    printf "%s\n" "${DELFIN_PYTHON}"
    return 0
  fi
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

python_version_at_least() {
  local python_bin="$1"
  local major="$2"
  local minor="$3"
  "${python_bin}" - "$major" "$minor" <<'PY'
import sys
target_major = int(sys.argv[1])
target_minor = int(sys.argv[2])
sys.exit(0 if sys.version_info >= (target_major, target_minor) else 1)
PY
}

desired_censo_git_ref() {
  local python_bin="$1"
  printf '%s\n' "${CENSO_V3_REF}"
}

compatible_censo_git_ref() {
  local python_bin="$1"
  if python_version_at_least "${python_bin}" 3 12; then
    printf '%s\n' "${CENSO_V3_REF}"
  else
    printf '%s\n' "${CENSO_V2_REF}"
  fi
}

installed_censo_version() {
  local python_bin="$1"
  "${python_bin}" -m pip show censo 2>/dev/null | awk '/^Version: / {print $2; exit}'
}

version_at_least() {
  local actual="$1"
  local wanted="$2"
  ACTUAL="${actual}" WANTED="${wanted}" python - <<'PY'
import os
import sys

def parse(text: str):
    parts = [int(part) for part in str(text).strip().split(".") if part.strip()]
    while len(parts) < 3:
        parts.append(0)
    return tuple(parts[:3])

sys.exit(0 if parse(os.environ["ACTUAL"]) >= parse(os.environ["WANTED"]) else 1)
PY
}

minimum_supported_censo_version() {
  local python_bin="$1"
  if python_version_at_least "${python_bin}" 3 12; then
    printf '%s\n' "3.0.0"
  else
    printf '%s\n' "2.1.4"
  fi
}

describe_python_environment() {
  local python_bin="$1"
  local env_label=""
  local env_path=""

  if [ -n "${CONDA_PREFIX:-}" ]; then
    env_label="conda/micromamba"
    env_path="${CONDA_PREFIX}"
  elif [ -n "${VIRTUAL_ENV:-}" ]; then
    env_label="venv"
    env_path="${VIRTUAL_ENV}"
  else
    env_label="python"
    env_path="$("${python_bin}" - <<'PY'
import sys
print(sys.prefix)
PY
)"
  fi

  printf '%s (%s)' "${python_bin}" "${env_label}: ${env_path}"
}

python_has_module() {
  local python_bin="$1"
  local module="$2"
  "${python_bin}" -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('${module}') else 1)" >/dev/null 2>&1
}

# Where this environment's commands go: the venv's bin/, not the directory of
# the interpreter a venv's python links to. Resolved, that was /usr/bin for a
# venv made from the system python, and writing c2anmr, nmrplot and anmr there
# failed with "Permission denied" -- and stopped the run before cclib, nglview
# and Packmol.
python_bin_dir() {
  local python_bin="$1"
  "${python_bin}" - <<'PY'
import sysconfig
print(sysconfig.get_path("scripts"))
PY
}

have_nmrplot() {
  have nmrplot || have nmrplot.py
}

cleanup_tmp_dir() {
  local tmp_dir="${1:-}"
  if [ -n "${tmp_dir}" ] && [ -d "${tmp_dir}" ]; then
    rm -rf "${tmp_dir}"
  fi
}

download_to_file() {
  local url="$1"
  local target="$2"
  local python_bin=""
  if have curl; then
    curl -fsSL "$url" -o "$target"
    return $?
  fi
  if have wget; then
    wget -qO "$target" "$url"
    return $?
  fi
  python_bin="$(detect_python)" || return 1
  "${python_bin}" - "$url" "$target" <<'PY'
import sys
import urllib.request
url, target = sys.argv[1:3]
with urllib.request.urlopen(url) as response, open(target, "wb") as handle:
    handle.write(response.read())
PY
}

sha256_of() {
  # The digest alone on stdout; status 1 when no tool can compute one.
  if have sha256sum; then
    sha256sum "$1" | awk '{print $1}'
  elif have shasum; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif have python3; then
    python3 -I -c 'import hashlib, sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' "$1"
  else
    return 1
  fi
}

fetch_verified() {
  # download_to_file, then the SHA-256 check; status 1 and no file at
  # DEST when the download fails, no tool can compute the digest, or the
  # digest differs.
  local url="$1" dest="$2" want="$3" got
  if [ -z "${want}" ]; then
    warn "no SHA-256 is pinned for ${url}; it is not downloaded"
    return 1
  fi
  if ! download_to_file "${url}" "${dest}"; then
    rm -f "${dest}"
    warn "download failed: ${url}"
    return 1
  fi
  if ! got="$(sha256_of "${dest}")"; then
    rm -f "${dest}"
    warn "no sha256sum, shasum or python3 to verify $(basename "${dest}"); an unchecked download is not installed"
    return 1
  fi
  if [ "${got}" != "${want}" ]; then
    rm -f "${dest}"
    warn "checksum mismatch for $(basename "${dest}") from ${url}: expected ${want}, got ${got}. The file was removed and nothing was installed from it."
    return 1
  fi
  log "checksum verified: $(basename "${dest}")"
}

# The SHA-256 of a CENSO helper script (c2anmr, nmrplot) at a pinned
# commit, from raw.githubusercontent.com on 2026-10-16; the git blob ids
# of the downloads equal those in the commit's tree. The scripts are the
# same files in v2.1.3 and v2.1.4.
censo_helper_sha256() {
  case "$1/$2" in
    1a82c5ad5e0b46c17cbbc93ecd7890a3b8c8a1e3/c2anmr) printf '%s\n' f90e9e7a34056f6af4b9fae1e519919bae9f74923770cad8673950e8993e16a1 ;;
    1a82c5ad5e0b46c17cbbc93ecd7890a3b8c8a1e3/nmrplot) printf '%s\n' 6368c3dc83646701ac0058a5da2026ab0ca7a789f5f612c913c2a4ee1ee551c9 ;;
    db667bde403c916bd20e8a42544fa3b0bed8f7d7/c2anmr) printf '%s\n' fec66d2ce55143a95562d52b228b14f3f925719b7e45f03687876a6562f8a600 ;;
    db667bde403c916bd20e8a42544fa3b0bed8f7d7/nmrplot) printf '%s\n' 6276ea3d112b8c8b3ecb7fa4a9e148b001087b6d7c0916d619c834bb960fd384 ;;
    *) return 1 ;;
  esac
}

install_censo_helper_wrappers_from_repo() {
  local python_bin="$1"
  local censo_version=""
  local helper_store=""
  local bin_dir=""
  local helper_name=""
  local helper_target=""
  local wrapper_target=""
  local ref=""
  local digest=""
  local found_any=0
  local found_c2anmr=0
  local found_nmrplot=0
  local refs=()

  helper_store="${ROOT}/share/censo_scripts"
  bin_dir="$(python_bin_dir "${python_bin}")"
  mkdir -p "${helper_store}" "${bin_dir}"

  # The scripts of the pinned commit of the installed CENSO line first.
  censo_version="$("${python_bin}" -c "from importlib.metadata import version; print(version('censo'))" 2>/dev/null || true)"
  case "${censo_version}" in
    2.*) refs=("${CENSO_V2_REF}" "${CENSO_V3_REF}") ;;
    *)   refs=("${CENSO_V3_REF}" "${CENSO_V2_REF}") ;;
  esac

  for ref in "${refs[@]}"; do
    for helper_name in c2anmr nmrplot; do
      if [ "${helper_name}" = "c2anmr" ] && [ "${found_c2anmr}" = "1" ]; then
        continue
      fi
      if [ "${helper_name}" = "nmrplot" ] && [ "${found_nmrplot}" = "1" ]; then
        continue
      fi
      if ! digest="$(censo_helper_sha256 "${ref}" "${helper_name}")"; then
        warn "no SHA-256 is pinned for ${helper_name} at CENSO ${ref}; not downloaded"
        continue
      fi

      helper_target="${helper_store}/${helper_name}.py"
      log "downloading CENSO helper ${helper_name} at ${ref}..."
      if ! fetch_verified "https://raw.githubusercontent.com/grimme-lab/CENSO/${ref}/src/censo/scripts/${helper_name}.py" \
          "${helper_target}" "${digest}"; then
        continue
      fi
      chmod 755 "${helper_target}"

      wrapper_target="${bin_dir}/${helper_name}"
      cat > "${wrapper_target}" <<EOF
#!/usr/bin/env bash
exec "${python_bin}" "${helper_target}" "\$@"
EOF
      chmod 755 "${wrapper_target}"

      if [ "${helper_name}" = "nmrplot" ]; then
        cp -f "${wrapper_target}" "${bin_dir}/nmrplot.py"
        chmod 755 "${bin_dir}/nmrplot.py"
        found_nmrplot=1
      else
        found_c2anmr=1
      fi
      found_any=1
    done

    if [ "${found_c2anmr}" = "1" ] && [ "${found_nmrplot}" = "1" ]; then
      break
    fi
  done

  if [ "${found_c2anmr}" != "1" ]; then
    warn "Could not install c2anmr from the pinned CENSO sources"
  fi
  if [ "${found_nmrplot}" != "1" ]; then
    warn "Could not install nmrplot from the pinned CENSO sources"
  fi
  [ "${found_any}" = "1" ]
}

ensure_censo_scripts() {
  local python_bin="$1"
  if have c2anmr && have_nmrplot && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "CENSO helper scripts already available"
    return 0
  fi

  log "ensuring CENSO helper scripts (c2anmr, nmrplot)..."
  if ! "${python_bin}" -m pip install "censo[scripts]" 2>&1 | tee -a "${LOG_DIR}/censo_install.log"; then
    warn "censo[scripts] install failed; falling back to helper-script wrappers from the official CENSO source"
    install_censo_helper_wrappers_from_repo "${python_bin}" || true
  fi

  if ! have c2anmr || ! have_nmrplot; then
    install_censo_helper_wrappers_from_repo "${python_bin}" || true
  fi

  if have c2anmr && have_nmrplot; then
    log "CENSO helper scripts installed successfully"
  else
    warn "CENSO is installed but helper scripts c2anmr/nmrplot are still missing"
  fi
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
  local censo_ref=""
  local censo_git_url=""
  local compat_ref=""
  local min_version=""
  local current_version=""
  local install_ok=1

  if [ "${INSTALL_CENSO}" != "1" ]; then
    log "CENSO: skipped (INSTALL_CENSO=0)"
    return 0
  fi

  min_version="$(minimum_supported_censo_version "${python_bin}")"
  current_version="$(installed_censo_version "${python_bin}" 2>/dev/null || true)"

  if python_has_module "${python_bin}" "censo" && [ "${FORCE_REINSTALL}" != "1" ]; then
    if [ "${CENSO_PREFER_LATEST}" = "1" ] && [ -n "${current_version}" ] && ! version_at_least "${current_version}" "${min_version}"; then
      log "CENSO: installed version ${current_version} is below required ${min_version}; upgrading"
    else
      log "CENSO: already installed${current_version:+ (${current_version})}"
      ensure_censo_scripts "${python_bin}"
      return 0
    fi
  fi

  censo_ref="$(desired_censo_git_ref "${python_bin}")"
  compat_ref="$(compatible_censo_git_ref "${python_bin}")"

  for ref in "${censo_ref}" "${compat_ref}"; do
    [ -n "${ref}" ] || continue
    censo_git_url="git+https://github.com/grimme-lab/CENSO.git@${ref}"
    log "installing CENSO at commit ${ref} from GitHub..."
    if "${python_bin}" -m pip install --upgrade "${censo_git_url}" \
      2>&1 | tee -a "${LOG_DIR}/censo_install.log"; then
      install_ok=0
      current_version="$(installed_censo_version "${python_bin}" 2>/dev/null || true)"
      if [ -n "${current_version}" ] && version_at_least "${current_version}" "${min_version}"; then
        break
      fi
      if [ "${ref}" != "${compat_ref}" ]; then
        warn "Installed CENSO ${current_version:-unknown} is below preferred minimum ${min_version}; trying fallback."
      fi
    else
      warn "CENSO install attempt for ${ref} failed; trying next candidate if available."
    fi
  done

  if [ "${install_ok}" -eq 0 ] && python_has_module "${python_bin}" "censo"; then
    current_version="$(installed_censo_version "${python_bin}" 2>/dev/null || true)"
    log "CENSO installed successfully${current_version:+ (${current_version})}"
    ensure_censo_scripts "${python_bin}"
  else
    warn "CENSO installation failed — check ${LOG_DIR}/censo_install.log"
  fi
}

# ---------------------------------------------------------------------------
install_anmr() {
  local python_bin="$1"

  if [ "${INSTALL_ANMR}" != "1" ]; then
    log "ANMR: skipped (INSTALL_ANMR=0)"
    return 0
  fi

  if have anmr && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "ANMR: already available at $(command -v anmr)"
    return 0
  fi

  if [ -n "${ENSO_ANMR_RELEASE_TAG:-}" ]; then
    warn "ENSO_ANMR_RELEASE_TAG is no longer read; name the file with ENSO_ANMR_URL and its digest with ENSO_ANMR_SHA256"
  fi

  local bin_dir tmp_dir
  bin_dir="$(python_bin_dir "${python_bin}")"
  tmp_dir="$(mktemp -d)"

  log "downloading ANMR from ${ENSO_ANMR_URL}..."
  if ! fetch_verified "${ENSO_ANMR_URL}" "${tmp_dir}/anmr" "${ENSO_ANMR_SHA256}"; then
    warn "ANMR was not installed"
    cleanup_tmp_dir "${tmp_dir}"
    return 0
  fi

  mkdir -p "${bin_dir}"
  install -m 755 "${tmp_dir}/anmr" "${bin_dir}/anmr"
  log "ANMR installed successfully to ${bin_dir}/anmr"
  cleanup_tmp_dir "${tmp_dir}"
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
install_cclib() {
  local python_bin="$1"

  if [ "${INSTALL_CCLIB}" != "1" ]; then
    log "cclib: skipped (INSTALL_CCLIB=0)"
    return 0
  fi

  if python_has_module "${python_bin}" "cclib" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "cclib: already installed"
    return 0
  fi

  log "installing cclib..."
  "${python_bin}" -m pip install cclib 2>&1 | tee -a "${LOG_DIR}/cclib_install.log"

  if python_has_module "${python_bin}" "cclib"; then
    local version
    version="$("${python_bin}" -c "from importlib.metadata import version; print(version('cclib'))" 2>/dev/null || echo "?")"
    log "cclib v${version} installed successfully"
  else
    warn "cclib installation failed"
  fi
}

# ---------------------------------------------------------------------------
install_nglview() {
  local python_bin="$1"

  if [ "${INSTALL_NGLVIEW}" != "1" ]; then
    log "nglview: skipped (INSTALL_NGLVIEW=0)"
    return 0
  fi

  if python_has_module "${python_bin}" "nglview" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "nglview: already installed"
    return 0
  fi

  log "installing nglview..."
  "${python_bin}" -m pip install nglview 2>&1 | tee -a "${LOG_DIR}/nglview_install.log"

  if python_has_module "${python_bin}" "nglview"; then
    local version
    version="$("${python_bin}" -c "from importlib.metadata import version; print(version('nglview'))" 2>/dev/null || echo "?")"
    log "nglview v${version} installed successfully"
  else
    warn "nglview installation failed"
  fi
}

# ---------------------------------------------------------------------------
# micromamba wherever DELFIN or the user put it; fetched when there is none.
# The shared micromamba fetcher (delfin/installers/fetch_micromamba.sh): one
# pinned release per platform, held to its SHA-256. Beside this script when
# the installer was staged into a tool directory, one level up in a checkout
# or an installed package.
micromamba_fetcher() {
  local here candidate
  here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  for candidate in "${here}/fetch_micromamba.sh" "${here}/../installers/fetch_micromamba.sh"; do
    if [ -f "${candidate}" ]; then
      printf "%s\n" "${candidate}"
      return 0
    fi
  done
  return 1
}

find_micromamba() {
  local candidate
  for candidate in "${MAMBA_EXE:-}" "$(command -v micromamba 2>/dev/null || true)" \
      "${DELFIN_QM_TOOLS_ROOT:-${HOME}/.delfin/qm_tools}/bin/micromamba" \
      "${ROOT}/bin/micromamba" "${HOME}/micromamba/bin/micromamba" "${HOME}/.local/bin/micromamba"; do
    if [ -n "${candidate}" ] && [ -x "${candidate}" ]; then
      printf "%s\n" "${candidate}"
      return 0
    fi
  done
  local fetcher
  fetcher="$(micromamba_fetcher)" || return 1
  bash "${fetcher}" "${ROOT}/bin/micromamba" >&2 || return 1
  printf "%s\n" "${ROOT}/bin/micromamba"
}

# Packmol in an environment of its own, linked beside DELFIN's other commands.
#
# `micromamba install packmol` without a prefix put it into whatever base
# environment micromamba had, which is on nobody's PATH -- and a machine with
# DELFIN's own micromamba, but none on the PATH, was told Packmol "requires
# conda" and got nothing.
install_packmol() {
  local python_bin="$1"
  if [ "${INSTALL_PACKMOL}" != "1" ]; then
    log "Packmol: skipped (INSTALL_PACKMOL=0)"
    return 0
  fi

  if have packmol && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "Packmol: already installed at $(command -v packmol)"
    return 0
  fi

  local mamba env_dir bin_dir
  if ! mamba="$(find_micromamba)"; then
    warn "Packmol needs micromamba, and none was found or could be fetched."
    warn "  Install manually: conda install -c conda-forge packmol"
    return 0
  fi
  env_dir="${ROOT}/.mamba_env/packmol"
  bin_dir="$(python_bin_dir "${python_bin}")"
  if [ ! -x "${env_dir}/bin/packmol" ] || [ "${FORCE_REINSTALL}" = "1" ]; then
    log "installing Packmol into ${env_dir} with ${mamba}..."
    "${mamba}" create -y -p "${env_dir}" -c conda-forge packmol 2>&1 | tee -a "${LOG_DIR}/packmol_install.log" || true
  fi
  if [ -x "${env_dir}/bin/packmol" ]; then
    mkdir -p "${bin_dir}"
    ln -sfn "${env_dir}/bin/packmol" "${bin_dir}/packmol"
    log "Packmol installed: ${bin_dir}/packmol -> ${env_dir}/bin/packmol"
  else
    warn "Packmol could not be installed; see ${LOG_DIR}/packmol_install.log"
  fi
}

# ---------------------------------------------------------------------------
summary() {
  local python_bin

  log "============================================"
  log "  Analysis tools installation summary"
  log "============================================"

  if python_bin="$(detect_python)"; then
    for mod_label in "morfeus:morfeus-ml" "censo:CENSO" "cclib:cclib" "nglview:nglview"; do
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

    if have anmr; then
      printf "  %-16s %s\n" "ANMR" "installed ($(command -v anmr))"
    else
      printf "  %-16s %s\n" "ANMR" "not installed"
    fi

    if have c2anmr; then
      printf "  %-16s %s\n" "c2anmr" "installed ($(command -v c2anmr))"
    else
      printf "  %-16s %s\n" "c2anmr" "not installed"
    fi

    if have nmrplot || have nmrplot.py; then
      printf "  %-16s %s\n" "nmrplot" "installed ($(command -v nmrplot || command -v nmrplot.py))"
    else
      printf "  %-16s %s\n" "nmrplot" "not installed"
    fi

    if have packmol; then
      printf "  %-16s %s\n" "Packmol" "installed ($(command -v packmol))"
    else
      printf "  %-16s %s\n" "Packmol" "not installed"
    fi
  fi

  printf "  %-16s %s\n" "tools root" "${ROOT}"
  printf "\n"
  printf "Usage in Python:\n"
  printf "  from delfin.analysis_tools.morfeus_wrapper import full_steric_report\n"
  printf "  from delfin.analysis_tools.censo_wrapper import run_censo\n"
  printf "  from delfin.analysis_tools.multiwfn_wrapper import bond_order_analysis\n"
  printf "  from delfin.analysis_tools.cclib_wrapper import parse_output\n"
  printf "  from delfin.analysis_tools.packmol_wrapper import solvate\n"
}

# ---------------------------------------------------------------------------
main() {
  local python_bin
  python_bin="$(detect_python)" || die "python/python3 not found"

  mkdir -p "${LOG_DIR}"
  log "using Python environment: $(describe_python_environment "${python_bin}")"
  install_morfeus "${python_bin}"
  install_censo "${python_bin}"
  install_anmr "${python_bin}"
  install_multiwfn "${python_bin}"
  install_cclib "${python_bin}"
  install_nglview "${python_bin}"
  install_packmol "${python_bin}"
  summary
}

main "$@"
