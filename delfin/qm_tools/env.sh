#!/usr/bin/env bash

# Local QM toolchain for DELFIN / ChemDarwin experiments.
# Usage:
#   source /path/to/DELFIN/delfin/qm_tools/env.sh

_this_file="${BASH_SOURCE[0]:-$0}"
_qm_root="$(cd "$(dirname "${_this_file}")" && pwd)"
_xtb4stda_default="${_qm_root}/share/xtb4stda"
_xtb4stda_home="${XTB4STDAHOME:-${_xtb4stda_default}}"

case ":${PATH}:" in
  *":${_qm_root}/bin:"*) ;;
  *) export PATH="${_qm_root}/bin:${PATH}" ;;
esac

if [ -z "${XTB4STDAHOME:-}" ] && [ -d "${_xtb4stda_default}" ] && [ ${#_xtb4stda_default} -gt 60 ]; then
  _xtb4stda_short="${HOME}/.delfin_xtb4stda"
  if [ -d "${_xtb4stda_short}" ] && [ -f "${_xtb4stda_short}/.param_stda1.xtb" ]; then
    _xtb4stda_home="${_xtb4stda_short}"
  elif [ ! -e "${_xtb4stda_short}" ] || [ -L "${_xtb4stda_short}" ]; then
    if ln -sfn "${_xtb4stda_default}" "${_xtb4stda_short}" 2>/dev/null; then
      _xtb4stda_home="${_xtb4stda_short}"
    fi
  fi
fi

export DELFIN_QM_ROOT="${_qm_root}"
export XTB4STDAHOME="${_xtb4stda_home}"
export STD2HOME="${STD2HOME:-${_qm_root}}"
export OMP_STACKSIZE="${OMP_STACKSIZE:-4G}"
