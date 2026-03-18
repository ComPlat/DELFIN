#!/usr/bin/env bash

# CSP toolchain for DELFIN (Crystal Structure Prediction).
# Usage:
#   source /path/to/DELFIN/delfin/csp_tools/env.sh

_this_file="${BASH_SOURCE[0]:-$0}"
_csp_root="$(cd "$(dirname "${_this_file}")" && pwd)"

case ":${PATH}:" in
  *":${_csp_root}/bin:"*) ;;
  *) export PATH="${_csp_root}/bin:${PATH}" ;;
esac

export DELFIN_CSP_TOOLS_ROOT="${_csp_root}"
