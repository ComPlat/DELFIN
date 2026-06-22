#!/usr/bin/env bash

# MLP toolchain for DELFIN (Machine Learning Potentials).
# Usage:
#   source /path/to/DELFIN/delfin/mlp_tools/env.sh

_this_file="${BASH_SOURCE[0]:-$0}"
_mlp_root="$(cd "$(dirname "${_this_file}")" && pwd)"

case ":${PATH}:" in
  *":${_mlp_root}/bin:"*) ;;
  *) export PATH="${_mlp_root}/bin:${PATH}" ;;
esac

export DELFIN_MLP_TOOLS_ROOT="${_mlp_root}"
