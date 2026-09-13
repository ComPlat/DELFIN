#!/usr/bin/env bash
# Kept so existing instructions and older dashboards keep working. The installer
# is universal now: install_delfin.sh beside this file. This runs its core
# profile, which is what this script used to install; options are passed on.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$here/install_delfin.sh" --profile core "$@"
