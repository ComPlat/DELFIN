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
#   INSTALL_RXNMAPPER            (default: 0)  RXNMapper
#   INSTALL_DEEPCHEM             (default: 0)  DeepChem
#   INSTALL_MOLSIMPLIFY          (default: 0)  molSimplify
#   INSTALL_ARCHITECTOR          (default: 0)  architector
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
INSTALL_RXNMAPPER="${INSTALL_RXNMAPPER:-${INSTALL_ALL}}"
INSTALL_DEEPCHEM="${INSTALL_DEEPCHEM:-${INSTALL_ALL}}"
INSTALL_MOLSIMPLIFY="${INSTALL_MOLSIMPLIFY:-${INSTALL_ALL}}"
INSTALL_ARCHITECTOR="${INSTALL_ARCHITECTOR:-${INSTALL_ALL}}"
INSTALL_PLOTLY="${INSTALL_PLOTLY:-${INSTALL_ALL}}"
FORCE_REINSTALL="${FORCE_REINSTALL:-0}"

log()  { printf "[ai_tools] %s\n" "$*"; }
warn() { printf "[ai_tools] WARNING: %s\n" "$*" >&2; }
die()  { printf "[ai_tools] ERROR: %s\n" "$*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

detect_python() {
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
  "${python_bin}" -m pip install "${packages[@]}" 2>&1 | tee -a "${LOG_DIR}/${log_name}"

  if python_has_module "${python_bin}" "${module}"; then
    log "${label} installed successfully"
  else
    warn "${label} installation failed — check ${LOG_DIR}/${log_name}"
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
  pip_install "${python_bin}" INSTALL_RXNMAPPER     "rxnmapper"      rxnmapper

  # Screening
  pip_install "${python_bin}" INSTALL_DEEPCHEM "deepchem" deepchem

  # Metal Complex ML
  pip_install "${python_bin}" INSTALL_MOLSIMPLIFY "molSimplify" molSimplify
  pip_install "${python_bin}" INSTALL_ARCHITECTOR "architector"  architector

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
    "rxnmapper:RXNMapper" \
    "deepchem:DeepChem" \
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
}

main "$@"
