#!/usr/bin/env bash
# DELFIN installer. See delfin/installers/install_delfin.sh --help.
#
#   bash install.sh                 DELFIN, ORCA wiring, QM/analysis tools, Ketcher
#   bash install.sh --all           everything, ML stacks included
#   bash install.sh --only crest,gxtb
#   bash install.sh --dry-run       print the plan, change nothing
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DELFIN_REPO="${DELFIN_REPO:-$here}" exec bash "$here/delfin/installers/install_delfin.sh" "$@"
