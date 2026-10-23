#!/usr/bin/env bash
# Kept so existing instructions and older dashboards keep working. The installer
# is universal now and lives in delfin/installers/install_delfin.sh; this runs
# its core profile, which is what this script used to install. Any option given
# here is passed on, so --profile standard or --all still work.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$here/../delfin/installers/install_delfin.sh" --profile core "$@"
