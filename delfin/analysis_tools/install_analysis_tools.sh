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
  printf '%s\n' "main"
}

compatible_censo_git_ref() {
  local python_bin="$1"
  if python_version_at_least "${python_bin}" 3 12; then
    printf '%s\n' "main"
  else
    printf '%s\n' "v2.1.4"
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

python_bin_dir() {
  local python_bin="$1"
  "${python_bin}" - <<'PY'
from pathlib import Path
import sys
print(Path(sys.executable).resolve().parent)
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

extract_archive() {
  local archive_path="$1"
  local extract_dir="$2"
  local python_bin=""
  python_bin="$(detect_python)" || return 1
  "${python_bin}" - "$archive_path" "$extract_dir" <<'PY'
import shutil
import sys
import tarfile
import zipfile
from pathlib import Path

archive = Path(sys.argv[1])
target = Path(sys.argv[2])
target.mkdir(parents=True, exist_ok=True)
name = archive.name.lower()
if name.endswith(".zip"):
    with zipfile.ZipFile(archive) as zf:
        zf.extractall(target)
elif any(name.endswith(suffix) for suffix in (".tar.gz", ".tgz", ".tar.xz", ".txz", ".tar.bz2", ".tbz2", ".tar")):
    with tarfile.open(archive, "r:*") as tf:
        tf.extractall(target)
else:
    shutil.copy2(archive, target / archive.name)
PY
}

find_extracted_anmr() {
  local extract_dir="$1"
  local python_bin=""
  python_bin="$(detect_python)" || return 1
  "${python_bin}" - "$extract_dir" <<'PY'
import os
import sys
from pathlib import Path

root = Path(sys.argv[1])
candidates = []
for path in root.rglob("*"):
    if not path.is_file():
        continue
    name = path.name.lower()
    if name == "anmr" or name.startswith("anmr."):
        candidates.append(path)
for path in sorted(candidates, key=lambda item: (len(item.parts), item.name)):
    print(path)
    break
PY
}

select_enso_release_asset() {
  local python_bin="$1"
  local json_file="$2"
  "${python_bin}" - "$json_file" <<'PY'
import json
import sys
from pathlib import Path

data = json.loads(Path(sys.argv[1]).read_text())
assets = data.get("assets", [])


def score(name: str) -> tuple[int, int]:
    text = name.lower()
    value = 0
    if "linux" in text:
        value += 20
    if "x86_64" in text or "amd64" in text:
        value += 10
    if "windows" in text or "mac" in text or "darwin" in text:
        value -= 50
    if text.endswith(".zip") or ".tar" in text:
        value += 5
    if "anmr" in text:
        value += 3
    return value, -len(text)


ranked = sorted(
    (
        asset for asset in assets
        if asset.get("browser_download_url") and not asset.get("draft")
    ),
    key=lambda asset: score(asset.get("name", "")),
    reverse=True,
)
for asset in ranked:
    text = asset.get("name", "").lower()
    if "windows" in text or "darwin" in text or "mac" in text:
        continue
    print(asset["browser_download_url"])
    break
PY
}

find_censo_helper_source() {
  local extract_dir="$1"
  local helper_name="$2"
  local python_bin=""
  python_bin="$(detect_python)" || return 1
  "${python_bin}" - "$extract_dir" "$helper_name" <<'PY'
import sys
from pathlib import Path

root = Path(sys.argv[1])
helper = sys.argv[2]
names = []
if helper == "nmrplot":
    names = ["nmrplot.py", "nmrplot"]
elif helper == "c2anmr":
    names = ["c2anmr.py", "c2anmr"]
else:
    names = [helper]

for name in names:
    matches = sorted(root.rglob(name), key=lambda item: (len(item.parts), str(item)))
    if matches:
        print(matches[0])
        break
PY
}

install_censo_helper_wrappers_from_repo() {
  local python_bin="$1"
  local censo_version=""
  local helper_store=""
  local bin_dir=""
  local helper_name=""
  local helper_source=""
  local helper_target=""
  local wrapper_target=""
  local candidate_ref=""
  local found_any=0
  local found_c2anmr=0
  local found_nmrplot=0
  local tmp_dir=""
  local archive_path=""
  local extract_dir=""
  local repo_url=""

  helper_store="${ROOT}/share/censo_scripts"
  bin_dir="$(python_bin_dir "${python_bin}")"
  mkdir -p "${helper_store}" "${bin_dir}"

  censo_version="$("${python_bin}" -c "from importlib.metadata import version; print(version('censo'))" 2>/dev/null || true)"
  for candidate_ref in "${censo_version:+v${censo_version}}" "v2.1.3" "main"; do
    [ -n "${candidate_ref}" ] || continue
    tmp_dir="$(mktemp -d)"
    archive_path="${tmp_dir}/censo_repo.tar.gz"
    extract_dir="${tmp_dir}/extract"
    mkdir -p "${extract_dir}"
    if [ "${candidate_ref}" = "main" ]; then
      repo_url="https://codeload.github.com/grimme-lab/CENSO/tar.gz/refs/heads/main"
    else
      repo_url="https://codeload.github.com/grimme-lab/CENSO/tar.gz/refs/tags/${candidate_ref}"
    fi

    log "downloading CENSO helper scripts from ${repo_url}..."
    if ! download_to_file "${repo_url}" "${archive_path}" 2>&1 | tee -a "${LOG_DIR}/censo_install.log" >/dev/null; then
      cleanup_tmp_dir "${tmp_dir}"
      continue
    fi
    if ! extract_archive "${archive_path}" "${extract_dir}" 2>&1 | tee -a "${LOG_DIR}/censo_install.log" >/dev/null; then
      cleanup_tmp_dir "${tmp_dir}"
      continue
    fi

    for helper_name in c2anmr nmrplot; do
      if [ "${helper_name}" = "c2anmr" ] && [ "${found_c2anmr}" = "1" ]; then
        continue
      fi
      if [ "${helper_name}" = "nmrplot" ] && [ "${found_nmrplot}" = "1" ]; then
        continue
      fi

      helper_source="$(find_censo_helper_source "${extract_dir}" "${helper_name}")"
      if [ -z "${helper_source}" ] || [ ! -f "${helper_source}" ]; then
        continue
      fi

      helper_target="${helper_store}/${helper_name}.py"
      cp -f "${helper_source}" "${helper_target}"
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

    cleanup_tmp_dir "${tmp_dir}"
    if [ "${found_c2anmr}" = "1" ] && [ "${found_nmrplot}" = "1" ]; then
      break
    fi
  done

  if [ "${found_c2anmr}" != "1" ]; then
    warn "Could not locate c2anmr in the downloaded CENSO sources"
  fi
  if [ "${found_nmrplot}" != "1" ]; then
    warn "Could not locate nmrplot in the downloaded CENSO sources"
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
    if [ "${ref}" = "main" ]; then
      censo_git_url="git+https://github.com/grimme-lab/CENSO.git"
      log "installing latest CENSO from GitHub main..."
    else
      censo_git_url="git+https://github.com/grimme-lab/CENSO.git@${ref}"
      log "installing compatible CENSO fallback (${ref}) from GitHub..."
    fi
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

  local api_url release_selector bin_dir tmp_dir json_file archive_url archive_name archive_path extract_dir anmr_source
  release_selector="${ENSO_ANMR_RELEASE_TAG:-latest}"
  if [ "${release_selector}" = "latest" ]; then
    api_url="https://api.github.com/repos/grimme-lab/enso/releases/latest"
  else
    api_url="https://api.github.com/repos/grimme-lab/enso/releases/tags/${release_selector}"
  fi
  bin_dir="$(python_bin_dir "${python_bin}")"
  tmp_dir="$(mktemp -d)"
  json_file="${tmp_dir}/enso_release.json"
  archive_path="${tmp_dir}/enso_asset"
  extract_dir="${tmp_dir}/extract"
  mkdir -p "${extract_dir}"

  log "querying ENSO release metadata for ANMR..."
  download_to_file "${api_url}" "${json_file}" 2>&1 | tee -a "${LOG_DIR}/anmr_install.log" >/dev/null

  archive_url="$(select_enso_release_asset "${python_bin}" "${json_file}")"

  if [ -z "${archive_url}" ]; then
    warn "Could not determine a Linux ENSO release asset for ANMR"
    cleanup_tmp_dir "${tmp_dir}"
    return 0
  fi

  archive_name="$(basename "${archive_url}")"
  archive_path="${tmp_dir}/${archive_name}"
  log "downloading ENSO asset ${archive_name}..."
  download_to_file "${archive_url}" "${archive_path}" 2>&1 | tee -a "${LOG_DIR}/anmr_install.log" >/dev/null
  extract_archive "${archive_path}" "${extract_dir}" 2>&1 | tee -a "${LOG_DIR}/anmr_install.log" >/dev/null

  anmr_source="$(find_extracted_anmr "${extract_dir}")"
  if [ -z "${anmr_source}" ] || [ ! -f "${anmr_source}" ]; then
    warn "Downloaded ENSO asset did not contain an ANMR executable"
    cleanup_tmp_dir "${tmp_dir}"
    return 0
  fi

  mkdir -p "${bin_dir}"
  install -m 755 "${anmr_source}" "${bin_dir}/anmr"
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
install_packmol() {
  if [ "${INSTALL_PACKMOL}" != "1" ]; then
    log "Packmol: skipped (INSTALL_PACKMOL=0)"
    return 0
  fi

  if have packmol && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "Packmol: already installed at $(command -v packmol)"
    return 0
  fi

  if have micromamba || have mamba || have conda; then
    local conda_cmd
    if have micromamba; then
      conda_cmd="micromamba"
    elif have mamba; then
      conda_cmd="mamba"
    else
      conda_cmd="conda"
    fi

    log "installing Packmol via ${conda_cmd}..."
    "${conda_cmd}" install -y -c conda-forge packmol 2>&1 | tee -a "${LOG_DIR}/packmol_install.log" || true

    if have packmol; then
      log "Packmol installed successfully via ${conda_cmd}"
      return 0
    fi
  fi

  warn "Packmol installation requires conda/micromamba/mamba."
  warn "  Install manually: conda install -c conda-forge packmol"
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
  install_packmol
  summary
}

main "$@"
