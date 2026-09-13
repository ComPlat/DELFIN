#!/usr/bin/env bash
# The installer repairs what it finds broken and changes nothing that is in
# order, so verifying and repairing is running it again. For a check that
# changes nothing at all: install_delfin.sh --dry-run, or check_delfin_orca.sh.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$here/../delfin/installers/install_delfin.sh" --profile core "$@"
