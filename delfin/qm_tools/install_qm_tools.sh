#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${ROOT}/bin"
DOWNLOAD_DIR="${ROOT}/downloads"
BUILD_DIR="${ROOT}/build"
SHARE_DIR="${ROOT}/share"
XTB4STDA_SHARE_DIR="${SHARE_DIR}/xtb4stda"
MAMBA_ENV="${ROOT}/.mamba_env"
ORCA_LOCAL_DIR="${ROOT}/third_party/orca"

XTB4STDA_URL="${XTB4STDA_URL:-https://github.com/grimme-lab/xtb4stda/releases/download/v1.0/xtb4stda}"
STDA_URL="${STDA_URL:-https://github.com/grimme-lab/xtb4stda/releases/download/v1.0/stda_v1.6.1}"
XTB4STDA_RUNTIME_BASE_URL="${XTB4STDA_RUNTIME_BASE_URL:-https://raw.githubusercontent.com/grimme-lab/xtb4stda/master}"
STD2_TAG="${STD2_TAG:-v2.0.1}"
STD2_SRC_URL="${STD2_SRC_URL:-https://github.com/grimme-lab/std2/archive/refs/tags/${STD2_TAG}.tar.gz}"

USE_SYSTEM_TOOLS="${USE_SYSTEM_TOOLS:-1}"
INSTALL_STD2_FROM_SOURCE="${INSTALL_STD2_FROM_SOURCE:-0}"
FORCE_REDOWNLOAD="${FORCE_REDOWNLOAD:-0}"
STD2_BIN="${STD2_BIN:-}"

log() {
  printf "[qm_tools] %s\n" "$*"
}

warn() {
  printf "[qm_tools] WARNING: %s\n" "$*" >&2
}

die() {
  printf "[qm_tools] ERROR: %s\n" "$*" >&2
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

ensure_dirs() {
  mkdir -p "${BIN_DIR}" "${DOWNLOAD_DIR}" "${BUILD_DIR}" "${SHARE_DIR}" "${XTB4STDA_SHARE_DIR}" "${ROOT}/third_party" "${ROOT}/docs" "${ROOT}/logs"
}

download_file() {
  local url="$1"
  local dest="$2"

  if [[ -f "${dest}" && "${FORCE_REDOWNLOAD}" != "1" ]]; then
    log "reuse download ${dest}"
    return
  fi

  have curl || die "curl not found"
  log "download ${url}"
  curl -fL --retry 3 --retry-delay 2 -o "${dest}" "${url}"
}

link_into_bin() {
  local source_path="$1"
  local target_name="${2:-$(basename "${source_path}")}"
  local target_path="${BIN_DIR}/${target_name}"

  [[ -e "${source_path}" ]] || die "source path missing: ${source_path}"
  ln -sfn "${source_path}" "${target_path}"
  log "linked ${target_name} -> ${source_path}"
}

install_file_into_bin() {
  local source_path="$1"
  local target_name="${2:-$(basename "${source_path}")}"
  install -m 755 "${source_path}" "${BIN_DIR}/${target_name}"
  log "installed ${target_name}"
}

detect_existing_tool() {
  local prog="$1"
  local detected="" resolved=""
  if [[ "${USE_SYSTEM_TOOLS}" == "1" ]] && have "${prog}"; then
    detected="$(command -v "${prog}")"
    if [[ -n "${detected}" ]]; then
      if resolved="$(readlink -f "${detected}" 2>/dev/null)"; then
        if [[ -n "${resolved}" && -x "${resolved}" && "${resolved}" != "${BIN_DIR}/${prog}" ]]; then
          printf "%s\n" "${resolved}"
          return 0
        fi
      fi
      if [[ "${detected}" != "${BIN_DIR}/${prog}" && -x "${detected}" ]]; then
        printf "%s\n" "${detected}"
        return 0
      fi
    fi
  fi
  return 1
}

ensure_micromamba() {
  if have micromamba; then
    printf "%s\n" "$(command -v micromamba)"
    return 0
  fi
  if have conda; then
    printf "%s\n" "$(command -v conda)"
    return 0
  fi
  return 1
}

install_conda_stack() {
  local mamba

  mamba="$(ensure_micromamba)" || die "micromamba/conda not found; cannot auto-install xtb/crest/dftb+"

  if [[ ! -x "${MAMBA_ENV}/bin/xtb" || ! -x "${MAMBA_ENV}/bin/crest" || ! -x "${MAMBA_ENV}/bin/dftb+" ]]; then
    log "create/update local conda env at ${MAMBA_ENV}"
    if [[ "$(basename "${mamba}")" == "micromamba" ]]; then
      "${mamba}" create -y -p "${MAMBA_ENV}" -c conda-forge xtb crest dftbplus
    else
      "${mamba}" create -y -p "${MAMBA_ENV}" -c conda-forge xtb crest dftbplus
    fi
  fi

  link_into_bin "${MAMBA_ENV}/bin/xtb" xtb
  link_into_bin "${MAMBA_ENV}/bin/crest" crest
  link_into_bin "${MAMBA_ENV}/bin/dftb+" dftb+
}

install_xtb() {
  local path
  if path="$(detect_existing_tool xtb)"; then
    link_into_bin "${path}" xtb
    return
  fi
  install_conda_stack
}

install_crest() {
  local path
  if path="$(detect_existing_tool crest)"; then
    link_into_bin "${path}" crest
    return
  fi
  install_conda_stack
}

install_dftbplus() {
  local path
  if path="$(detect_existing_tool dftb+)"; then
    link_into_bin "${path}" dftb+
    return
  fi
  install_conda_stack
}

install_xtb4stda_bundle() {
  local xtb4stda_dl="${DOWNLOAD_DIR}/xtb4stda"
  local stda_dl="${DOWNLOAD_DIR}/stda_v1.6.1"
  local runtime_files=(
    ".xtb4stdarc"
    ".param_stda1.xtb"
    ".param_stda2.xtb"
    ".param_gbsa_acetone"
    ".param_gbsa_acetonitrile"
    ".param_gbsa_benzene"
    ".param_gbsa_ch2cl2"
    ".param_gbsa_chcl3"
    ".param_gbsa_cs2"
    ".param_gbsa_dmso"
    ".param_gbsa_ether"
    ".param_gbsa_h2o"
    ".param_gbsa_methanol"
    ".param_gbsa_thf"
    ".param_gbsa_toluene"
  )
  local runtime_file

  download_file "${XTB4STDA_URL}" "${xtb4stda_dl}"
  download_file "${STDA_URL}" "${stda_dl}"
  for runtime_file in "${runtime_files[@]}"; do
    download_file "${XTB4STDA_RUNTIME_BASE_URL}/${runtime_file}" "${XTB4STDA_SHARE_DIR}/${runtime_file}"
  done

  install_file_into_bin "${xtb4stda_dl}" xtb4stda
  install_file_into_bin "${stda_dl}" stda_v1.6.1
  link_into_bin "${BIN_DIR}/stda_v1.6.1" stda

  [[ -f "${XTB4STDA_SHARE_DIR}/.xtb4stdarc" ]] || die "xtb4stda runtime incomplete: missing .xtb4stdarc"
  [[ -f "${XTB4STDA_SHARE_DIR}/.param_stda1.xtb" ]] || die "xtb4stda runtime incomplete: missing .param_stda1.xtb"
  [[ -f "${XTB4STDA_SHARE_DIR}/.param_stda2.xtb" ]] || die "xtb4stda runtime incomplete: missing .param_stda2.xtb"
}

build_std2_from_source() {
  local std2_tar="${DOWNLOAD_DIR}/std2-${STD2_TAG}.tar.gz"
  local src_dir="${BUILD_DIR}/std2-${STD2_TAG}"
  local build_dir="${src_dir}/_build"
  local meson_bin ninja_bin fc_bin cc_bin std2_built

  download_file "${STD2_SRC_URL}" "${std2_tar}"
  rm -rf "${src_dir}"
  mkdir -p "${src_dir}"
  tar -xzf "${std2_tar}" -C "${src_dir}" --strip-components=1

  if have meson && have ninja && have gfortran && have gcc; then
    meson_bin="$(command -v meson)"
    ninja_bin="$(command -v ninja)"
    fc_bin="$(command -v gfortran)"
    cc_bin="$(command -v gcc)"
  else
    die "std2 source build requested, but meson/ninja/gfortran/gcc are not all available"
  fi

  rm -rf "${build_dir}"
  FC="${fc_bin}" CC="${cc_bin}" "${meson_bin}" setup "${build_dir}" "${src_dir}" -Dla_backend=openblas
  "${meson_bin}" compile -C "${build_dir}"

  std2_built="${build_dir}/std2"
  [[ -x "${std2_built}" ]] || die "std2 build finished without executable at ${std2_built}"
  link_into_bin "${std2_built}" std2
}

install_std2() {
  local path

  if [[ -n "${STD2_BIN}" ]]; then
    link_into_bin "${STD2_BIN}" std2
    return
  fi

  if path="$(detect_existing_tool std2)"; then
    link_into_bin "${path}" std2
    return
  fi

  if [[ "${INSTALL_STD2_FROM_SOURCE}" == "1" ]]; then
    build_std2_from_source
    return
  fi

  warn "std2 not found in PATH. Set STD2_BIN=/path/to/std2 or INSTALL_STD2_FROM_SOURCE=1"
}

summary() {
  local python_bin=""

  log "installation summary"
  for prog in xtb crest std2 stda xtb4stda dftb+; do
    if [[ -x "${BIN_DIR}/${prog}" ]]; then
      printf "  %-12s %s\n" "${prog}" "${BIN_DIR}/${prog}"
    else
      printf "  %-12s %s\n" "${prog}" "missing"
    fi
  done
  printf "  %-12s %s\n" "xtb4stda-data" "${XTB4STDA_SHARE_DIR}"
  printf "  %-12s %s\n" "orca" "not managed by this installer"
  printf "\n"
  printf "Activate with:\n"
  printf "  source %s/env.sh\n" "${ROOT}"
  printf "Verify with:\n"
  printf "  %s/check_qm_tools.sh\n" "${ROOT}"
  printf "Workflow preflight:\n"
  printf "  python -m delfin.cli tadf_xtb --check\n"

  if python_bin="$(detect_python)"; then
    if ! python_has_module "${python_bin}" "rdkit"; then
      warn "RDKit is missing in ${python_bin}; 'delfin tadf_xtb --smiles ...' will fail until project Python dependencies are installed"
    fi
  else
    warn "python/python3 not found in PATH; cannot validate RDKit for the tadf_xtb workflow"
  fi
}

main() {
  ensure_dirs
  install_xtb
  install_crest
  install_dftbplus
  install_xtb4stda_bundle
  install_std2
  summary
}

main "$@"
