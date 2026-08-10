#!/usr/bin/env bash

set -euo pipefail

# Where the tools go. Beside this script unless told otherwise -- and it was
# never told otherwise, which is the whole of the problem: the Settings tab
# ran a staged copy and filled the user's own directory, the Submit tab ran
# the packaged copy and filled the package, and only one of the two was being
# searched. A user pressed Install, the tools arrived, and they were reported
# missing. Both callers now name the directory the resolver reads.
ROOT="${DELFIN_QM_TOOLS_ROOT:-${DELFIN_QM_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
mkdir -p "${ROOT}"
ROOT="$(cd "${ROOT}" && pwd)"
BIN_DIR="${ROOT}/bin"
DOWNLOAD_DIR="${ROOT}/downloads"
BUILD_DIR="${ROOT}/build"
SHARE_DIR="${ROOT}/share"
XTB4STDA_SHARE_DIR="${SHARE_DIR}/xtb4stda"
MAMBA_ENV="${ROOT}/.mamba_env"
ORCA_LOCAL_DIR="${ROOT}/third_party/orca"

# g-xTB ships as a statically linked xtb 6.7.1 of its own -- the method is not
# in a released tblite yet, and the ordinary xtb accepts --gxtb and silently
# runs GFN2 -- so it is installed beside xtb rather than instead of it.
GXTB_TAG="${GXTB_TAG:-v2.0.1}"
GXTB_ASSET="${GXTB_ASSET:-xtb-6.7.1-gxtb-140526-linux-x86_64.tar.xz}"
GXTB_URL="${GXTB_URL:-https://github.com/grimme-lab/g-xtb/releases/download/${GXTB_TAG}/${GXTB_ASSET}}"
XTB4STDA_URL="${XTB4STDA_URL:-https://github.com/grimme-lab/xtb4stda/releases/download/v1.0/xtb4stda}"
STDA_URL="${STDA_URL:-https://github.com/grimme-lab/xtb4stda/releases/download/v1.0/stda_v1.6.1}"
XTB4STDA_RUNTIME_BASE_URL="${XTB4STDA_RUNTIME_BASE_URL:-https://raw.githubusercontent.com/grimme-lab/xtb4stda/master}"
STD2_TAG="${STD2_TAG:-v2.0.1}"
STD2_SRC_URL="${STD2_SRC_URL:-https://github.com/grimme-lab/std2/archive/refs/tags/${STD2_TAG}.tar.gz}"

# Whether a tool already on this system may stand in for one of ours at all.
USE_SYSTEM_TOOLS="${USE_SYSTEM_TOOLS:-1}"
# Whether it is *preferred* over building the managed environment. Off: an
# install should install. See install_xtb.
PREFER_SYSTEM_TOOLS="${PREFER_SYSTEM_TOOLS:-0}"
# The build we know answers. Unpinned, "install xtb" means whatever
# conda-forge published this morning, so two machines set up a week apart run
# two different programs -- and a fault inside xtb is a fault of one build, not
# of the structure. Measured on this one: water GFN2 -5.070325081194 Eh,
# GFN-FF -0.327252352397 Eh, and a doublet optimisation that finishes.
# Override to move deliberately, not by accident.
# The oldest xtb we are willing to run. 6.6.1 has a formatting fault inside
# its own optimiser -- a doublet optimisation dies with "Fortran runtime
# error: Missing comma between descriptors" at optimizer.f90:639 while 6.7.1
# finishes the same job normally. A floor rather than a pin, so newer builds
# are still taken.
XTB_MINIMUM="${XTB_MINIMUM:-6.7.1}"
INSTALL_STD2_FROM_SOURCE="${INSTALL_STD2_FROM_SOURCE:-0}"
FORCE_REDOWNLOAD="${FORCE_REDOWNLOAD:-0}"
FORCE_CONDA_UPDATE="${FORCE_CONDA_UPDATE:-0}"
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

# One environment per tool, not one environment for all of them.
#
# Solved together, conda-forge answers **xtb 6.6.1** -- the build whose
# optimiser dies with a Fortran format error -- because that is the newest xtb
# compatible with crest and dftbplus *at the same time*. Solved apart, each
# tool gets its own newest: xtb 6.7.1, crest 3.0.2, dftbplus 25.1. Measured
# with three dry runs of the solver, and it is why a user could press Install
# and end up with a broken xtb while the person beside them, whose xtb came
# from somewhere else, had no trouble at all.
#
# Pinning xtb inside the joint environment would have worked too, and it costs
# dftbplus 25.1 -> 21.1. Separate environments cost nothing but disk.
install_conda_tool() {
  local prog="$1" spec="$2" env_dir="${MAMBA_ENV}/$3" binary="${4:-$1}"
  local mamba

  mamba="$(ensure_micromamba)" || return 1

  if [[ ! -x "${env_dir}/bin/${binary}" ]]; then
    log "create ${prog} environment at ${env_dir}"
    "${mamba}" create -y -p "${env_dir}" -c conda-forge "${spec}" || return 1
  elif [[ "${FORCE_CONDA_UPDATE}" == "1" ]]; then
    log "update ${prog} in ${env_dir}"
    # 'update', not 'install': an unversioned 'install xtb' is a no-op when
    # xtb is already there -- the spec is satisfied -- so it never replaces an
    # outdated build.
    "${mamba}" update -y -p "${env_dir}" -c conda-forge "${spec}" \
      || log "WARNING: update of ${prog} did not complete"
  fi

  [[ -x "${env_dir}/bin/${binary}" ]] || return 1
  link_into_bin "${env_dir}/bin/${binary}" "${prog}"
  # Now that this tool has one of its own, the shared environment may have
  # nothing left in it worth keeping.
  retire_legacy_env
}

# The environment the tools used to share, once nothing points into it.
#
# Leaving it costs more than the disk. A binary that is still on disk can
# still be reached: an explicit path saved in the settings, a PATH somebody
# exported, a hard-coded fallback. So the xtb whose optimiser dies with
# "Missing comma between descriptors" can go on being the one that runs, long
# after a working one has been installed beside it.
#
# Removed only after the replacement is in place and only when no link in
# bin/ resolves into it any more -- a cleanup that runs first would leave a
# failed install with nothing at all.
retire_legacy_env() {
  local marker="${MAMBA_ENV}/conda-meta"
  [[ -d "${marker}" ]] || return 0

  local link resolved
  for link in "${BIN_DIR}"/*; do
    [[ -e "${link}" || -L "${link}" ]] || continue
    resolved="$(readlink -f "${link}" 2>/dev/null)" || continue
    case "${resolved}" in
      "${MAMBA_ENV}/bin/"*)
        log "keeping the shared environment: ${link##*/} still uses it"
        return 0
        ;;
    esac
  done

  local freed=""
  if have du; then
    freed="$(du -sh "${MAMBA_ENV}/bin" "${MAMBA_ENV}/lib" "${MAMBA_ENV}/share" \
             2>/dev/null | awk '{s+=$1} END {print ""}')"
  fi
  log "removing the shared environment nothing uses any more: ${MAMBA_ENV}"
  local part
  for part in bin lib share include conda-meta ssl etc libexec man sbin; do
    rm -rf "${MAMBA_ENV:?}/${part}"
  done
}

install_conda_stack() {
  # Kept for callers that want all three; each still gets its own environment.
  local failed=0
  install_conda_tool xtb "xtb>=${XTB_MINIMUM}" xtb || failed=1
  install_conda_tool crest crest crest || failed=1
  install_conda_tool dftb+ dftbplus dftbplus "dftb+" || failed=1
  return "${failed}"
}

# Install the build we know, and only adopt somebody else's if we cannot.
#
# Adopting first was the default, and it means "install" mostly did not
# install: it found whatever xtb stood first on the PATH and made a symlink to
# it. Two accounts on one cluster then run two different xtb builds from the
# same button -- one of them optimises and the other stops with a Fortran
# format error inside xtb's own optimiser, on the same structure. Pressing
# Install again cannot help, because it adopts the same broken build again.
#
# So the managed environment comes first and the system tool is the fallback,
# which keeps the case that made adoption the default in the first place: a
# cluster with no network and xtb behind a module. Whichever happened is said
# out loud, because "installed" and "adopted" are not the same claim.
install_managed_or_adopt() {
  local prog="$1" spec="$2" env_name="$3"
  local path binary="${prog}"
  if install_conda_tool "${prog}" "${spec}" "${env_name}" "${binary}"; then
    return 0
  fi
  log "WARNING: no managed environment could be built for ${prog}"
  if path="$(detect_existing_tool "${prog}")"; then
    log "adopting the ${prog} already on this system: ${path}"
    log "  (it was not installed by DELFIN, and DELFIN cannot vouch for it)"
    link_into_bin "${path}" "${prog}"
    return 0
  fi
  die "${prog} could not be installed and none was found on this system"
}

install_xtb() {
  if [[ "${PREFER_SYSTEM_TOOLS}" == "1" ]]; then
    local path
    if path="$(detect_existing_tool xtb)"; then
      log "using the xtb already on this system: ${path}"
      link_into_bin "${path}" xtb
      return
    fi
  fi
  install_managed_or_adopt xtb "xtb>=${XTB_MINIMUM}" xtb
}

install_gxtb() {
  local archive="${DOWNLOAD_DIR}/${GXTB_ASSET}"
  local sums="${archive}.sha256"
  local target="${ROOT}/third_party/gxtb"

  if [[ -x "${target}/bin/xtb" && "${FORCE_REDOWNLOAD}" != "1" ]]; then
    log "reuse g-xTB in ${target}"
    link_into_bin "${target}/bin/xtb" xtb-gxtb
    return
  fi

  download_file "${GXTB_URL}" "${archive}"
  # The release publishes a checksum beside the tarball.  A binary fetched
  # from the network and run on a user's structures is worth the one line it
  # costs to check that it is the one that was published.
  download_file "${GXTB_URL}.sha256" "${sums}"
  if have sha256sum; then
    local want got
    want="$(awk '{print $1}' "${sums}")"
    got="$(sha256sum "${archive}" | awk '{print $1}')"
    [[ "${want}" == "${got}" ]] || die "g-xTB checksum mismatch: expected ${want}, got ${got}"
    log "g-xTB checksum verified"
  else
    warn "sha256sum not found; the g-xTB download could not be verified"
  fi

  rm -rf "${target}"
  mkdir -p "${target}"
  # One directory deep in the tarball (xtb-6.7.1/), stripped so the layout
  # here does not depend on the version in its name.
  tar -xf "${archive}" -C "${target}" --strip-components=1
  [[ -x "${target}/bin/xtb" ]] || die "g-xTB archive did not contain bin/xtb"
  link_into_bin "${target}/bin/xtb" xtb-gxtb
}

install_crest() {
  if [[ "${PREFER_SYSTEM_TOOLS}" == "1" ]]; then
    local path
    if path="$(detect_existing_tool crest)"; then
      log "using the crest already on this system: ${path}"
      link_into_bin "${path}" crest
      return
    fi
  fi
  install_managed_or_adopt crest "crest" crest
}

install_dftbplus() {
  if [[ "${PREFER_SYSTEM_TOOLS}" == "1" ]]; then
    local path
    if path="$(detect_existing_tool dftb+)"; then
      log "using the dftb+ already on this system: ${path}"
      link_into_bin "${path}" dftb+
      return
    fi
  fi
  install_managed_or_adopt dftb+ "dftbplus" dftbplus
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
  for prog in xtb xtb-gxtb crest std2 stda xtb4stda dftb+; do
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

  if [[ $# -eq 0 ]]; then
    # No arguments: install everything (legacy behavior).
    install_xtb
    install_crest
    install_dftbplus
    install_xtb4stda_bundle
    install_std2
  else
    # Selective install: only the requested tools.
    for tool in "$@"; do
      case "${tool}" in
        xtb)          install_xtb ;;
        gxtb|g-xtb)   install_gxtb ;;
        crest)        install_crest ;;
        dftb+|dftbplus) install_dftbplus ;;
        xtb4stda|stda) install_xtb4stda_bundle ;;
        std2)         install_std2 ;;
        all)
          install_xtb
          install_crest
          install_dftbplus
          install_xtb4stda_bundle
          install_std2
          ;;
        *)
          die "Unknown tool: ${tool}. Available: xtb, gxtb, crest, dftb+, xtb4stda, stda, std2, all"
          ;;
      esac
    done
  fi

  summary
}

main "$@"
