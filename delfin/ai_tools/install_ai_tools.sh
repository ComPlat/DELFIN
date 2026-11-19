#!/usr/bin/env bash

set -euo pipefail

# ---------------------------------------------------------------------------
# AI/ML tools installer for DELFIN.
#
# Each tool is optional — set INSTALL_<TOOL>=1 to install.
# All default to 0 (not installed) to keep the environment lean.
#
# Environment variables:
#   INSTALL_MOLFORMER            (default: 0)  MoLFormer (HuggingFace transformers)
#   INSTALL_CHEMBERTA            (default: 0)  ChemBERTa (HuggingFace transformers)
#   INSTALL_UNIMOL               (default: 0)  Uni-Mol
#   INSTALL_REINVENT             (default: 0)  REINVENT4
#   INSTALL_SYNTHEMOL            (default: 0)  SyntheMol
#   INSTALL_GEOMOL               (default: 0)  GeoMol
#   INSTALL_TORSIONAL_DIFFUSION  (default: 0)  torsional-diffusion
#   INSTALL_MATTERGEN            (default: 0)  MatterGen
#   INSTALL_CDVAE                (default: 0)  CDVAE
#   INSTALL_AIZYNTHFINDER        (default: 0)  AiZynthFinder
#   INSTALL_LOCALRETRO           (default: 0)  LocalRetro
#   INSTALL_RXNMAPPER            (default: 0)  RXNMapper
#   INSTALL_DEEPCHEM             (default: 0)  DeepChem
#   INSTALL_ADMETLAB             (default: 0)  ADMETlab
#   INSTALL_MOLSIMPLIFY          (default: 0)  molSimplify
#   INSTALL_ARCHITECTOR          (default: 0)  architector
#   INSTALL_EPIC_MACE            (default: 0)  epic-MACE, in an environment of its own
#                                              (Python 3.7 + RDKit 2020.09, micromamba)
#   INSTALL_RXEMBED              (default: 0)  rxembed, own environment (Python 3.12, RDKit 2026.9)
#   INSTALL_MOLASSEMBLER         (default: 0)  Molassembler, own environment (Python 3.11)
#   INSTALL_METALLOGEN           (default: 0)  MetalloGen, own environment (Python 3.11; needs ORCA)
#   INSTALL_ORCA_MLIPS           (default: 0)  orca-mlips (ORCA ExtOpt MLIPs: AIMNet2, ORB), own
#                                              environment (Python 3.12, CUDA PyTorch)
#   INSTALL_ORCA_MLIPS_MACE      (default: 0)  orca-mlips MACE plugin (mace-torch) into that
#                                              environment; licensed, never part of INSTALL_ALL:
#                                              MACE weights are non-commercial (ASL)
#   INSTALL_ORCA_MLIPS_UMA       (default: 0)  orca-mlips UMA plugin (fairchem-core), own
#                                              environment; licensed, never part of INSTALL_ALL:
#                                              gated weights, needs your own HuggingFace login
#   INSTALL_ORCA_EXTERNAL_TOOLS  (default: 0)  orca-external-tools (ORCA ExtOpt server, AIMNet2),
#                                              own environment (Python 3.12, CUDA PyTorch)
#   INSTALL_PLOTLY               (default: 0)  plotly
#   INSTALL_ALL                  (default: 0)  Install everything
#   FORCE_REINSTALL              (default: 0)  Force reinstall
# ---------------------------------------------------------------------------

ROOT="${DELFIN_AI_TOOLS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="${ROOT}/logs"

INSTALL_ALL="${INSTALL_ALL:-0}"
INSTALL_MOLFORMER="${INSTALL_MOLFORMER:-${INSTALL_ALL}}"
INSTALL_CHEMBERTA="${INSTALL_CHEMBERTA:-${INSTALL_ALL}}"
INSTALL_UNIMOL="${INSTALL_UNIMOL:-${INSTALL_ALL}}"
INSTALL_REINVENT="${INSTALL_REINVENT:-${INSTALL_ALL}}"
INSTALL_SYNTHEMOL="${INSTALL_SYNTHEMOL:-${INSTALL_ALL}}"
INSTALL_GEOMOL="${INSTALL_GEOMOL:-${INSTALL_ALL}}"
INSTALL_TORSIONAL_DIFFUSION="${INSTALL_TORSIONAL_DIFFUSION:-${INSTALL_ALL}}"
INSTALL_MATTERGEN="${INSTALL_MATTERGEN:-${INSTALL_ALL}}"
INSTALL_CDVAE="${INSTALL_CDVAE:-${INSTALL_ALL}}"
INSTALL_AIZYNTHFINDER="${INSTALL_AIZYNTHFINDER:-${INSTALL_ALL}}"
INSTALL_LOCALRETRO="${INSTALL_LOCALRETRO:-${INSTALL_ALL}}"
INSTALL_RXNMAPPER="${INSTALL_RXNMAPPER:-${INSTALL_ALL}}"
INSTALL_DEEPCHEM="${INSTALL_DEEPCHEM:-${INSTALL_ALL}}"
INSTALL_ADMETLAB="${INSTALL_ADMETLAB:-${INSTALL_ALL}}"
INSTALL_MOLSIMPLIFY="${INSTALL_MOLSIMPLIFY:-${INSTALL_ALL}}"
INSTALL_ARCHITECTOR="${INSTALL_ARCHITECTOR:-${INSTALL_ALL}}"
INSTALL_EPIC_MACE="${INSTALL_EPIC_MACE:-${INSTALL_ALL}}"
INSTALL_RXEMBED="${INSTALL_RXEMBED:-${INSTALL_ALL}}"
INSTALL_MOLASSEMBLER="${INSTALL_MOLASSEMBLER:-${INSTALL_ALL}}"
INSTALL_METALLOGEN="${INSTALL_METALLOGEN:-${INSTALL_ALL}}"
INSTALL_ORCA_MLIPS="${INSTALL_ORCA_MLIPS:-${INSTALL_ALL}}"
# Licensed: installed only when named, never by INSTALL_ALL.
INSTALL_ORCA_MLIPS_MACE="${INSTALL_ORCA_MLIPS_MACE:-0}"
INSTALL_ORCA_MLIPS_UMA="${INSTALL_ORCA_MLIPS_UMA:-0}"
INSTALL_ORCA_EXTERNAL_TOOLS="${INSTALL_ORCA_EXTERNAL_TOOLS:-${INSTALL_ALL}}"
INSTALL_PLOTLY="${INSTALL_PLOTLY:-${INSTALL_ALL}}"
FORCE_REINSTALL="${FORCE_REINSTALL:-0}"

log()  { printf "[ai_tools] %s\n" "$*"; }
warn() { printf "[ai_tools] WARNING: %s\n" "$*" >&2; }
die()  { printf "[ai_tools] ERROR: %s\n" "$*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

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
  if have python; then command -v python; return 0; fi
  if have python3; then command -v python3; return 0; fi
  return 1
}

python_has_module() {
  local python_bin="$1" module="$2"
  "${python_bin}" -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('${module}') else 1)" >/dev/null 2>&1
}

pip_install() {
  local python_bin="$1" label="$2" module="$3"
  shift 3
  local packages=("$@")

  if [ "${!label:-0}" != "1" ]; then
    log "${label}: skipped (${label}=0)"
    return 0
  fi

  if python_has_module "${python_bin}" "${module}" && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "${label}: already installed"
    return 0
  fi

  local log_name
  log_name="$(echo "${label}" | tr '[:upper:]' '[:lower:]')_install.log"
  log "installing ${label}..."
  # One tool that pip cannot install is that tool's failure, not the run's.
  # Under pipefail it ended the whole script: REINVENT is not on PyPI, and
  # every tool after it was never tried.
  "${python_bin}" -m pip install "${packages[@]}" 2>&1 | tee -a "${LOG_DIR}/${log_name}" || true

  if python_has_module "${python_bin}" "${module}"; then
    log "${label} installed successfully"
  else
    warn "${label} installation failed — check ${LOG_DIR}/${log_name}"
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

# epic-MACE (Chernyshov & Pidko, JCTC 2024; GPL-3.0) in an environment of its own.
#
# It needs Python 3.7 and RDKit 2020.09, which DELFIN cannot run in, so it is
# never installed beside DELFIN: DELFIN starts it as an external program
# (delfin/common/external_builders.py) with the interpreter of this
# environment, which it finds at ${ROOT}/.mamba_env/epic_mace/bin/python;
# DELFIN_MACE_PYTHON overrides that. The package comes from a pinned GitHub
# commit: the PyPI release 0.5.0 has only the octahedron and the square, the
# commit adds the hapto ligands and TET/SPY/TBP/SAN.
EPIC_MACE_REF="${EPIC_MACE_REF:-efb5778e715ea461f80cf3bbc752929101bd0bb3}"
EPIC_MACE_URL="${EPIC_MACE_URL:-https://github.com/EPiCs-group/epic-mace/archive/${EPIC_MACE_REF}.tar.gz}"
EPIC_MACE_CONDA_SPECS="${EPIC_MACE_CONDA_SPECS:-python=3.7 rdkit=2020.09.5 numpy pyyaml pip}"

epic_mace_python() {
  # The environment's own interpreter, untouched by DELFIN's PYTHONPATH or a
  # user site of another Python version. Nor is its pip held to DELFIN's own
  # pins (PIP_CONSTRAINT, set for DELFIN's environment by delfin.installer):
  # an environment of its own exists because those pins cannot hold the tool.
  env -u PYTHONPATH -u PIP_CONSTRAINT PYTHONNOUSERSITE=1 "$@"
}

install_epic_mace() {
  if [ "${INSTALL_EPIC_MACE}" != "1" ]; then
    log "INSTALL_EPIC_MACE: skipped (INSTALL_EPIC_MACE=0)"
    return 0
  fi
  local env_dir="${ROOT}/.mamba_env/epic_mace"
  local py="${env_dir}/bin/python"
  local log_file="${LOG_DIR}/epic_mace_install.log"

  if [ -x "${py}" ] && epic_mace_python "${py}" -c "import mace" >/dev/null 2>&1 \
      && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "epic-MACE: already installed (${py})"
    return 0
  fi

  local mamba
  if ! mamba="$(find_micromamba)"; then
    warn "epic-MACE needs micromamba (or conda) for its Python 3.7 environment, and none was found or could be fetched."
    return 0
  fi
  if [ ! -x "${py}" ] || [ "${FORCE_REINSTALL}" = "1" ]; then
    log "creating the epic-MACE environment at ${env_dir} with ${mamba} (${EPIC_MACE_CONDA_SPECS})..."
    # shellcheck disable=SC2086
    "${mamba}" create -y -p "${env_dir}" -c conda-forge --override-channels ${EPIC_MACE_CONDA_SPECS} \
      2>&1 | tee -a "${log_file}" || true
  fi
  if [ ! -x "${py}" ]; then
    warn "epic-MACE: the environment could not be created; see ${log_file}"
    return 0
  fi
  log "installing epic-MACE ${EPIC_MACE_REF:0:8} into ${env_dir}..."
  # Its dependencies (numpy, pyyaml, rdkit) come from conda-forge above.
  epic_mace_python "${py}" -m pip install --no-deps --no-cache-dir --force-reinstall "${EPIC_MACE_URL}" \
    2>&1 | tee -a "${log_file}" || true
  if epic_mace_python "${py}" -c "import mace" >/dev/null 2>&1; then
    log "epic-MACE installed: ${py}"
    log "  DELFIN finds it there; DELFIN_MACE_PYTHON=<python> points it elsewhere."
  else
    warn "epic-MACE installation failed; see ${log_file}"
  fi
}

# rxembed, Molassembler and MetalloGen, each in an environment of its own.
#
# DELFIN starts them as external programs (delfin/common/external_builders.py,
# builds in delfin/cluster_bench/tool_workers/smiles_tools.py) with the
# interpreter of ${ROOT}/.mamba_env/<name>/bin/python; DELFIN_RXEMBED_PYTHON,
# DELFIN_MOLASSEMBLER_PYTHON and DELFIN_METALLOGEN_PYTHON override that.
# rxembed (MIT) needs Python >= 3.12 and RDKit >= 2026.9.1, newer than DELFIN's;
# it comes from a pinned GitHub commit (no PyPI release). Molassembler
# (scine-molassembler, BSD-3) is the PyPI wheel for CPython 3.11 on x86_64 Linux, held to its
# SHA-256 like the archives (no wheel for other platforms). MetalloGen (BSD-3) comes
# from a pinned GitHub commit; DELFIN runs it only with its ORCA calculator
# (GFN2-xTB through ORCA), so it needs ORCA, never Gaussian.
#
# Every version is exact: the tool itself is installed with --no-deps, its
# dependencies are named with == (rxembed: the paper benchmark's environment,
# KONKURRENZ_RXEMBED env_freeze; Molassembler: its benchmark environment), and
# every downloaded archive or wheel is held to its SHA-256 before pip sees it.
RXEMBED_REF="${RXEMBED_REF:-9c0ba5a6048bdbcba94450df5d6a68067bd250c8}"
RXEMBED_URL="${RXEMBED_URL:-https://github.com/aligfellow/rxembed/archive/${RXEMBED_REF}.tar.gz}"
RXEMBED_SHA256="${RXEMBED_SHA256:-58ad75e9f2fc652843f185540306560eba4258809ac991de10a34b297e2dfea5}"
RXEMBED_CONDA_SPECS="${RXEMBED_CONDA_SPECS:-python=3.12.16 pip}"
RXEMBED_PIP_DEPS="${RXEMBED_PIP_DEPS:-numpy==2.2.2 pillow==11.1.0 rdkit==2026.9.1 networkx==3.7 xyzgraph==1.7.0}"
MOLASSEMBLER_URL="${MOLASSEMBLER_URL:-https://files.pythonhosted.org/packages/f6/82/9c297a7368b0010dbde0a5b5e0dd5f4186cc59a7664f967b5d8e7580be3e/scine_molassembler-3.0.1-cp311-cp311-manylinux_2_28_x86_64.whl}"
MOLASSEMBLER_SHA256="${MOLASSEMBLER_SHA256:-b1035ce455ed2418d9b493fb19ebf9411eb7ac0d783ec5ed9081e838d54fec2b}"
MOLASSEMBLER_CONDA_SPECS="${MOLASSEMBLER_CONDA_SPECS:-python=3.11.18 pip}"
MOLASSEMBLER_PIP_DEPS="${MOLASSEMBLER_PIP_DEPS:-scine-utilities==10.1.0 numpy==2.4.6 scipy==1.17.1}"
METALLOGEN_REF="${METALLOGEN_REF:-fa26e2adb996ed2256e528e1c97fbc1c363c0808}"
METALLOGEN_URL="${METALLOGEN_URL:-https://github.com/kyunghoonlee777/MetalloGen/archive/${METALLOGEN_REF}.tar.gz}"
METALLOGEN_SHA256="${METALLOGEN_SHA256:-1a8a4f6a0223b7aa138b2eb0f57af4be6a84e056ed5c972bcbac2d6fdf3ab250}"
METALLOGEN_CONDA_SPECS="${METALLOGEN_CONDA_SPECS:-python=3.11.18 rdkit=2026.09.2 numpy=2.4.6 scipy=1.17.1 cclib=1.8.1 pulp=3.3.2 pip}"
METALLOGEN_PIP_DEPS="${METALLOGEN_PIP_DEPS:-}"

sha256_of() {
  if have sha256sum; then sha256sum "$1" | awk '{print $1}'
  elif have shasum; then shasum -a 256 "$1" | awk '{print $1}'
  elif have openssl; then openssl dgst -sha256 -r "$1" | awk '{print $1}'
  else return 1
  fi
}

fetch_checked() {
  # fetch_checked URL SHA256 DEST: download, and keep DEST only when its SHA-256 is SHA256.
  local url="$1" want="$2" dest="$3" got
  rm -f "${dest}"
  if have curl; then curl -fsSL -o "${dest}" "${url}" || return 1
  elif have wget; then wget -q -O "${dest}" "${url}" || return 1
  else warn "neither curl nor wget to fetch ${url}"; return 1
  fi
  if ! got="$(sha256_of "${dest}")" || [ "${got}" != "${want}" ]; then
    warn "${url}: SHA-256 ${got:-unknown} is not the pinned ${want}; not installed"
    rm -f "${dest}"
    return 1
  fi
}

install_own_env() {
  # install_own_env SWITCH LABEL ENV_NAME MODULE "CONDA_SPECS" "PIP_DEPS" SPEC_OR_URL SHA256
  # The tool goes in with --no-deps after its exactly pinned dependencies; with a SHA256 the
  # SPEC_OR_URL is an archive that is downloaded and checked first.
  local switch="$1" label="$2" name="$3" module="$4" specs="$5" deps="$6" pip_spec="$7" sha="$8"
  if [ "${!switch:-0}" != "1" ]; then
    log "${switch}: skipped (${switch}=0)"
    return 0
  fi
  local env_dir="${ROOT}/.mamba_env/${name}"
  local py="${env_dir}/bin/python"
  local log_file="${LOG_DIR}/${name}_install.log"

  if [ -x "${py}" ] && epic_mace_python "${py}" -c "import ${module}" >/dev/null 2>&1 \
      && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "${label}: already installed (${py})"
    return 0
  fi
  local mamba
  if ! mamba="$(find_micromamba)"; then
    warn "${label} needs micromamba (or conda) for its own environment, and none was found or could be fetched."
    return 0
  fi
  if [ ! -x "${py}" ] || [ "${FORCE_REINSTALL}" = "1" ]; then
    log "creating the ${label} environment at ${env_dir} with ${mamba} (${specs})..."
    # shellcheck disable=SC2086
    "${mamba}" create -y -p "${env_dir}" -c conda-forge --override-channels ${specs} \
      2>&1 | tee -a "${log_file}" || true
  fi
  if [ ! -x "${py}" ]; then
    warn "${label}: the environment could not be created; see ${log_file}"
    return 0
  fi
  if [ -n "${deps}" ]; then
    log "installing the pinned dependencies of ${label} (${deps})..."
    # shellcheck disable=SC2086
    epic_mace_python "${py}" -m pip install --no-cache-dir ${deps} \
      2>&1 | tee -a "${log_file}" || true
  fi
  local target="${pip_spec}"
  if [ -n "${sha}" ]; then
    target="${env_dir}/${name}-${sha:0:12}.tar.gz"
    # a wheel keeps its file name: pip reads the Python and platform tags from it
    case "${pip_spec}" in *.whl) target="${env_dir}/${pip_spec##*/}" ;; esac
    if ! fetch_checked "${pip_spec}" "${sha}" "${target}"; then
      warn "${label} installation failed: the archive could not be fetched and checked; see ${log_file}"
      return 0
    fi
  fi
  log "installing ${label} (${pip_spec}) into ${env_dir}..."
  epic_mace_python "${py}" -m pip install --no-cache-dir --no-deps "${target}" \
    2>&1 | tee -a "${log_file}" || true
  if epic_mace_python "${py}" -c "import ${module}" >/dev/null 2>&1; then
    log "${label} installed: ${py}"
  else
    warn "${label} installation failed; see ${log_file}"
  fi
}

install_smiles_tools() {
  install_own_env INSTALL_RXEMBED rxembed rxembed rxembed \
    "${RXEMBED_CONDA_SPECS}" "${RXEMBED_PIP_DEPS}" "${RXEMBED_URL}" "${RXEMBED_SHA256}"
  install_own_env INSTALL_MOLASSEMBLER Molassembler molassembler scine_molassembler \
    "${MOLASSEMBLER_CONDA_SPECS}" "${MOLASSEMBLER_PIP_DEPS}" "${MOLASSEMBLER_URL}" "${MOLASSEMBLER_SHA256}"
  # Its dependencies (rdkit, numpy, scipy, cclib, pulp) come pinned from conda-forge above.
  install_own_env INSTALL_METALLOGEN MetalloGen metallogen MetalloGen \
    "${METALLOGEN_CONDA_SPECS}" "${METALLOGEN_PIP_DEPS}" "${METALLOGEN_URL}" "${METALLOGEN_SHA256}"
}

# orca-mlips (MIT), the program ORCA's ExtOpt starts as ProgExt for an MLIP GOAT
# or optimisation (delfin/mlp_tools/orca_mlips.py), in an environment of its own:
# it needs a CUDA build of PyTorch, and its MLIP stacks pin numpy and friends
# beyond what DELFIN's own environment holds. Only the two plugins whose weights
# are licence-clean for any use go in: aimnet2 (MIT) and orb (Apache-2.0); the
# uma (gated) and mace (non-commercial) plugins are not installed here.
# ExtOpt names ${ROOT}/.mamba_env/orca_mlips/bin/orca-mlips-<backend> by its
# absolute path, so a compute node needs nothing on its PATH;
# DELFIN_ORCA_MLIPS_PYTHON points DELFIN at another environment.
# Every version is exact, as resolved by pip for orca-mlips[aimnet2,orb]==1.3.0:
# aimnet 0.2.6 needs nvalchemi-toolkit-ops >= 0.4.1, which orb-models 0.7 excludes,
# so orb-models is the newest release that installs beside it. The PyPI torch
# wheel is the CUDA build (it runs on the CPU where a job has no GPU).
ORCA_MLIPS_SPEC="${ORCA_MLIPS_SPEC:-orca-mlips==1.3.0}"
ORCA_MLIPS_CONDA_SPECS="${ORCA_MLIPS_CONDA_SPECS:-python=3.12.16 pip}"
ORCA_MLIPS_PIP_DEPS="${ORCA_MLIPS_PIP_DEPS:-torch==2.15.0 aimnet==0.2.6 orb-models==0.5.5 ase==3.29.0 numpy==2.4.6}"

# The two plugins whose weights are not licence-clean for any use, opt-in
# (INSTALL_ORCA_MLIPS_MACE / INSTALL_ORCA_MLIPS_UMA, never INSTALL_ALL), each
# saying its terms before anything is fetched.
# MACE goes into the environment above: mace-torch 0.3.17 needs torch>=2.5 and
# e3nn==0.4.4, and pip's resolution in that environment only adds packages
# (torch 2.15.0, numpy 2.4.6, aimnet and orb-models stay as they are).
ORCA_MLIPS_MACE_PIP_DEPS="${ORCA_MLIPS_MACE_PIP_DEPS:-mace-torch==0.3.17 e3nn==0.4.4 torch==2.15.0 numpy==2.4.6}"
# UMA has an environment of its own: fairchem-core 2.25.0 pins torch~=2.13.0 and
# e3nn>=0.5, which the torch above and MACE's e3nn==0.4.4 both exclude (pip:
# ResolutionImpossible). orca-mlips goes in with --no-deps after them, as above.
ORCA_MLIPS_UMA_CONDA_SPECS="${ORCA_MLIPS_UMA_CONDA_SPECS:-python=3.12.16 pip}"
ORCA_MLIPS_UMA_PIP_DEPS="${ORCA_MLIPS_UMA_PIP_DEPS:-torch==2.13.0 fairchem-core==2.25.0 ase==3.29.0 numpy==2.4.6}"

install_orca_mlips() {
  # The MACE plugin goes into this environment, so it is built first.
  if [ "${INSTALL_ORCA_MLIPS_MACE}" = "1" ]; then
    INSTALL_ORCA_MLIPS=1
  fi
  install_own_env INSTALL_ORCA_MLIPS orca-mlips orca_mlips orca_mlips_plugins \
    "${ORCA_MLIPS_CONDA_SPECS}" "${ORCA_MLIPS_PIP_DEPS}" "${ORCA_MLIPS_SPEC}" ""
  install_orca_mlips_mace
  install_orca_mlips_uma
}

install_orca_mlips_mace() {
  if [ "${INSTALL_ORCA_MLIPS_MACE}" != "1" ]; then
    log "INSTALL_ORCA_MLIPS_MACE: skipped (INSTALL_ORCA_MLIPS_MACE=0)"
    return 0
  fi
  log "LICENCE: MACE foundation-model weights (MACE-OMOL-0, MACE-OFF, MACE-MP) are under the"
  log "LICENCE: Academic Software License -- non-commercial (academic) use only."
  log "LICENCE: See the MACE documentation, section Foundation Models (mace-docs.readthedocs.io)."
  local env_dir="${ROOT}/.mamba_env/orca_mlips"
  local py="${env_dir}/bin/python"
  local log_file="${LOG_DIR}/orca_mlips_mace_install.log"
  if [ ! -x "${py}" ]; then
    warn "orca-mlips MACE plugin: the orca-mlips environment (${env_dir}) is not there; see ${LOG_DIR}/orca_mlips_install.log"
    return 0
  fi
  if epic_mace_python "${py}" -c "import mace, orca_mlips_plugins" >/dev/null 2>&1 \
      && [ "${FORCE_REINSTALL}" != "1" ]; then
    log "orca-mlips MACE plugin: already installed (${py})"
    return 0
  fi
  log "installing the MACE stack of orca-mlips (${ORCA_MLIPS_MACE_PIP_DEPS}) into ${env_dir}..."
  # shellcheck disable=SC2086
  epic_mace_python "${py}" -m pip install --no-cache-dir ${ORCA_MLIPS_MACE_PIP_DEPS} \
    2>&1 | tee -a "${log_file}" || true
  if epic_mace_python "${py}" -c "import mace, orca_mlips_plugins" >/dev/null 2>&1; then
    log "orca-mlips MACE plugin installed: ${env_dir}/bin/orca-mlips-mace"
  else
    warn "orca-mlips MACE plugin installation failed; see ${log_file}"
  fi
}

install_orca_mlips_uma() {
  if [ "${INSTALL_ORCA_MLIPS_UMA}" != "1" ]; then
    log "INSTALL_ORCA_MLIPS_UMA: skipped (INSTALL_ORCA_MLIPS_UMA=0)"
    return 0
  fi
  log "LICENCE: UMA weights are under the FAIR Chemistry licence and gated on HuggingFace."
  log "LICENCE: They download only with YOUR OWN HuggingFace login, after you accept the model"
  log "LICENCE: terms at https://huggingface.co/facebook/UMA. Log in once with"
  log "LICENCE:   ${ROOT}/.mamba_env/orca_mlips_uma/bin/hf auth login"
  log "LICENCE: (older huggingface_hub: huggingface-cli login). DELFIN never handles the token."
  install_own_env INSTALL_ORCA_MLIPS_UMA orca-mlips-uma orca_mlips_uma "orca_mlips_plugins, fairchem.core" \
    "${ORCA_MLIPS_UMA_CONDA_SPECS}" "${ORCA_MLIPS_UMA_PIP_DEPS}" "${ORCA_MLIPS_SPEC}" ""
}

# orca-external-tools (FACCTs, GPL-3.0; github.com/faccts/orca-external-tools),
# the ExtOpt provider with a server: `oet_server aimnet2` keeps the model
# resident for a whole job and `oet_client` (ORCA's ProgExt) only forwards each
# request (its README, "Server"; ORCA 6.1 tutorial "ORCA with External
# Methods"). Installed in an environment of its own like orca-mlips, from the
# pinned release v2.1.0 held to its SHA-256 (it is not on PyPI); DELFIN starts
# its commands and never ships its code. Only the aimnet2 extra (MIT weights):
# its mace (non-commercial weights) and uma (gated weights; conflicts with
# aimnet2) extras belong in installations of their own. Every version is exact,
# as pip resolved oet[aimnet2] 2.1.0 with torch 2.15.0 (the PyPI CUDA build).
# The archive has no git metadata, so its version is told to setuptools_scm.
OET_REF="${OET_REF:-v2.1.0}"
OET_URL="${OET_URL:-https://github.com/faccts/orca-external-tools/archive/refs/tags/${OET_REF}.tar.gz}"
OET_SHA256="${OET_SHA256:-c939a392232ce418fafc4007319326328109b65ce8ecb4d2a96cac3f4d1a7492}"
OET_CONDA_SPECS="${OET_CONDA_SPECS:-python=3.12.16 pip}"
OET_PIP_DEPS="${OET_PIP_DEPS:-torch==2.15.0 aimnet==0.2.6 numpy==2.5.4 flask==3.1.3 waitress==3.0.2 requests==2.34.2 psutil==7.2.2}"

install_oet() {
  export SETUPTOOLS_SCM_PRETEND_VERSION_FOR_OET="${OET_REF#v}"
  install_own_env INSTALL_ORCA_EXTERNAL_TOOLS orca-external-tools oet oet \
    "${OET_CONDA_SPECS}" "${OET_PIP_DEPS}" "${OET_URL}" "${OET_SHA256}"
  unset SETUPTOOLS_SCM_PRETEND_VERSION_FOR_OET
}

# ---------------------------------------------------------------------------
main() {
  local python_bin
  python_bin="$(detect_python)" || die "python/python3 not found"
  mkdir -p "${LOG_DIR}"

  # Foundation Models (MoLFormer and ChemBERTa both use transformers)
  pip_install "${python_bin}" INSTALL_MOLFORMER  "transformers" transformers torch
  pip_install "${python_bin}" INSTALL_CHEMBERTA  "transformers" transformers torch
  pip_install "${python_bin}" INSTALL_UNIMOL     "unimol_tools" unimol_tools

  # Generative
  pip_install "${python_bin}" INSTALL_REINVENT      "reinvent"       reinvent
  pip_install "${python_bin}" INSTALL_SYNTHEMOL     "synthemol"      synthemol

  # Conformers
  pip_install "${python_bin}" INSTALL_GEOMOL               "geomol"               geomol
  pip_install "${python_bin}" INSTALL_TORSIONAL_DIFFUSION  "torsional_diffusion"  torsional-diffusion

  # Crystal Generation
  pip_install "${python_bin}" INSTALL_MATTERGEN  "mattergen"  mattergen
  pip_install "${python_bin}" INSTALL_CDVAE      "cdvae"      cdvae

  # Retrosynthesis
  pip_install "${python_bin}" INSTALL_AIZYNTHFINDER "aizynthfinder"  aizynthfinder
  pip_install "${python_bin}" INSTALL_LOCALRETRO    "localretro"     localretro
  pip_install "${python_bin}" INSTALL_RXNMAPPER     "rxnmapper"      rxnmapper

  # Screening
  pip_install "${python_bin}" INSTALL_DEEPCHEM  "deepchem"   deepchem
  pip_install "${python_bin}" INSTALL_ADMETLAB  "admetlab3"  admetlab3

  # Metal Complex ML
  pip_install "${python_bin}" INSTALL_MOLSIMPLIFY "molSimplify" molSimplify
  pip_install "${python_bin}" INSTALL_ARCHITECTOR "architector"  architector
  install_epic_mace
  install_smiles_tools
  install_orca_mlips
  install_oet

  # Visualization
  pip_install "${python_bin}" INSTALL_PLOTLY "plotly" plotly

  # Summary
  log "============================================"
  log "  AI tools installation summary"
  log "============================================"
  for mod_label in \
    "transformers:MoLFormer/ChemBERTa" \
    "unimol_tools:Uni-Mol" \
    "reinvent:REINVENT4" \
    "synthemol:SyntheMol" \
    "geomol:GeoMol" \
    "torsional_diffusion:torsional-diffusion" \
    "mattergen:MatterGen" \
    "cdvae:CDVAE" \
    "aizynthfinder:AiZynthFinder" \
    "localretro:LocalRetro" \
    "rxnmapper:RXNMapper" \
    "deepchem:DeepChem" \
    "admetlab3:ADMETlab" \
    "molSimplify:molSimplify" \
    "architector:architector" \
    "plotly:plotly"; do
    local mod="${mod_label%%:*}"
    local label="${mod_label##*:}"
    if python_has_module "${python_bin}" "${mod}"; then
      printf "  %-24s %s\n" "${label}" "installed"
    else
      printf "  %-24s %s\n" "${label}" "not installed"
    fi
  done
  if [ -x "${ROOT}/.mamba_env/epic_mace/bin/python" ] \
      && epic_mace_python "${ROOT}/.mamba_env/epic_mace/bin/python" -c "import mace" >/dev/null 2>&1; then
    printf "  %-24s %s\n" "epic-MACE" "installed (own environment)"
  else
    printf "  %-24s %s\n" "epic-MACE" "not installed"
  fi
  for env_label in "rxembed:rxembed:rxembed" "molassembler:scine_molassembler:Molassembler" \
      "metallogen:MetalloGen:MetalloGen" "orca_mlips:orca_mlips_plugins:orca-mlips" \
      "orca_mlips:mace:orca-mlips MACE" "orca_mlips_uma:fairchem.core:orca-mlips UMA" \
      "oet:oet:orca-external-tools"; do
    local env_name="${env_label%%:*}" rest="${env_label#*:}"
    local env_module="${rest%%:*}" env_shown="${rest#*:}"
    if [ -x "${ROOT}/.mamba_env/${env_name}/bin/python" ] \
        && epic_mace_python "${ROOT}/.mamba_env/${env_name}/bin/python" -c "import ${env_module}" >/dev/null 2>&1; then
      printf "  %-24s %s\n" "${env_shown}" "installed (own environment)"
    else
      printf "  %-24s %s\n" "${env_shown}" "not installed"
    fi
  done
}

main "$@"
