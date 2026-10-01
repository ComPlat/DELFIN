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
MICROMAMBA_URL="${MICROMAMBA_URL:-https://micro.mamba.pm/api/micromamba/linux-64/latest}"

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
  have curl || return 1
  local work="${ROOT}/downloads/micromamba-$$"
  mkdir -p "${work}" "${ROOT}/bin"
  if ! curl -fsSL "${MICROMAMBA_URL}" | tar -xj -C "${work}" bin/micromamba 2>/dev/null; then
    rm -rf "${work}"
    return 1
  fi
  install -m 755 "${work}/bin/micromamba" "${ROOT}/bin/micromamba"
  rm -rf "${work}"
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
  # user site of another Python version.
  env -u PYTHONPATH PYTHONNOUSERSITE=1 "$@"
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
}

main "$@"
