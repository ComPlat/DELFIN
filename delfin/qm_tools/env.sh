#!/usr/bin/env bash

# Local QM toolchain for DELFIN / ChemDarwin experiments.
# Usage:
#   source /path/to/DELFIN/delfin/qm_tools/env.sh

_this_file="${BASH_SOURCE[0]:-$0}"
_qm_root="$(cd "$(dirname "${_this_file}")" && pwd)"

case ":${PATH}:" in
  *":${_qm_root}/bin:"*) ;;
  *) export PATH="${_qm_root}/bin:${PATH}" ;;
esac

export DELFIN_QM_ROOT="${_qm_root}"
export XTB4STDAHOME="${XTB4STDAHOME:-${_qm_root}/share/xtb4stda}"
export STD2HOME="${STD2HOME:-${_qm_root}}"
export OMP_STACKSIZE="${OMP_STACKSIZE:-4G}"
