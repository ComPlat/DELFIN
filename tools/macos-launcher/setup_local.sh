#!/bin/bash
# Sets DELFIN up to run on this Mac, for a person without Xcode, Homebrew or
# a Python of their own. Everything goes into one folder of this account:
#
#   ~/Library/Application Support/DELFIN Launcher/local/
#     bin/micromamba   one pinned release, held to its SHA-256
#                      (fetch_micromamba.sh, the file DELFIN's own
#                      installers use, copied into this package)
#     env/             Python 3.11 from conda-forge, and DELFIN in it
#
# DELFIN comes from PyPI (delfin-complat, the version named in VERSION beside
# this file). install.sh is not used here: it sets up a Linux login node
# (ORCA, OpenMPI, the QM tools) and needs git and a compiler, which a Mac has
# only after the Xcode command line tools. ORCA cannot be downloaded for
# anyone (licence); DELFIN on the Mac runs the dashboard and the agent, whose
# commands run in macOS's own Seatbelt sandbox.
#
# Run again to repair or update. Nothing outside the folder above is changed,
# except pip's download cache in ~/Library/Caches.
#
# DELFIN_LOCAL_SOURCE=PATH installs that DELFIN checkout instead of the PyPI
# release (the CI job tests the tree it runs on this way).
PATH=/usr/bin:/bin:/usr/sbin:/sbin
export PATH
set -e
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
STORE=${DELFIN_LAUNCHER_STORE:-"$HOME/Library/Application Support/DELFIN Launcher"}
PREFIX="$STORE/local"
say() { printf '[DELFIN setup] %s\n' "$*"; }
fail() { printf '[DELFIN setup] ERROR: %s\n' "$*" >&2; exit 1; }

version=$(tr -d '[:space:]' < "$HERE/VERSION")
[[ $version =~ ^[0-9]+\.[0-9]+\.[0-9]+([a-z0-9.]*)?$ ]] || fail 'VERSION is not a version number.'
if [ -n "${DELFIN_LOCAL_SOURCE:-}" ]; then
    [ -f "$DELFIN_LOCAL_SOURCE/pyproject.toml" ] || fail "DELFIN_LOCAL_SOURCE is no DELFIN checkout: $DELFIN_LOCAL_SOURCE"
    spec="$DELFIN_LOCAL_SOURCE[agent,docs]"
else
    spec="delfin-complat[agent,docs]==$version"
fi

say "DELFIN $version goes to $PREFIX."
say 'This downloads about 1 GB and takes 5 to 20 minutes. Nothing outside that folder changes.'
free_kb=$(df -Pk "$HOME" | awk 'NR == 2 {print $4}')
if [ -n "$free_kb" ] && [ "$free_kb" -lt 4000000 ]; then
    fail "Only $((free_kb / 1024 / 1024)) GB free; DELFIN on this Mac needs about 4 GB."
fi
mkdir -p "$PREFIX/bin"
chmod 700 "$STORE" "$PREFIX"

mamba="$PREFIX/bin/micromamba"
if [ ! -x "$mamba" ]; then
    bash "$HERE/fetch_micromamba.sh" "$mamba" || fail 'micromamba could not be downloaded and verified (see the line above).'
fi
export MAMBA_ROOT_PREFIX="$PREFIX/mamba"
if [ ! -x "$PREFIX/env/bin/python" ]; then
    say 'Creating Python 3.11 (conda-forge).'
    "$mamba" create -y -q -p "$PREFIX/env" --override-channels -c conda-forge 'python=3.11' pip \
        || fail 'Python could not be installed (no connection, or a proxy: export https_proxy=http://HOST:PORT and run again).'
fi
python="$PREFIX/env/bin/python"
say "Installing $spec (pip)."
"$python" -m pip install --disable-pip-version-check --upgrade pip >/dev/null
"$python" -m pip install --disable-pip-version-check --upgrade "$spec" \
    || fail 'DELFIN could not be installed (the pip output above says why).'
"$python" -I -c 'import delfin.cli_voila, delfin.agent.where, voila' \
    || fail 'DELFIN is installed but its dashboard does not import.'
[ -x "$PREFIX/env/bin/delfin-voila" ] || fail 'delfin-voila is missing from the environment.'
say 'DELFIN is set up on this Mac. In DELFIN choose "DELFIN on this Mac".'
