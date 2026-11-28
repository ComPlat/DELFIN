#!/usr/bin/env bash
# Linux half of Setup-Wsl.ps1: prepares the WSL distribution "DELFIN" that
# the launcher's local mode ("This PC (WSL)") runs in.
#
#   wsl -d DELFIN -u root -- bash -s -- system <user>   < wsl_setup.sh
#   wsl -d DELFIN -u <user> -- bash -s -- delfin <repo> <commit> < wsl_setup.sh
#
# It arrives on stdin. Everything sits in functions and runs from the last
# line, with stdin from /dev/null: bash has read the whole script before
# anything runs, so apt-get or git cannot swallow the rest of it.
#
# Both steps can run again: what is already there is kept.
set -euo pipefail

say() { printf 'DELFIN setup: %s\n' "$*"; }

valid_user() { [[ "$1" =~ ^[a-z_][a-z0-9_-]{0,31}$ ]] && [ "$1" != root ]; }

system_step() {
    local user=$1
    valid_user "$user" || { say "invalid user name"; exit 2; }
    [ "$(id -u)" = 0 ] || { say "the system step needs root"; exit 2; }
    export DEBIAN_FRONTEND=noninteractive
    command -v apt-get >/dev/null || { say "this distribution has no apt-get"; exit 3; }
    # gfortran, meson, ninja and OpenBLAS: std2 (standard profile) is built from source.
    say "installing tmux, git, curl and the compilers"
    apt-get update -qq
    apt-get install -y -qq tmux git curl ca-certificates bzip2 build-essential gfortran meson ninja-build pkg-config libopenblas-dev >/dev/null
    if ! id "$user" >/dev/null 2>&1; then
        say "creating the Linux user $user"
        # No password: nothing logs in with one. Root work goes through
        # "wsl -u root" from Windows, which needs no password.
        useradd -m -s /bin/bash "$user"
    fi
    # The user WSL starts as. Kept if the distribution already names one.
    local conf=${DELFIN_SETUP_WSL_CONF:-/etc/wsl.conf}
    if ! grep -qs '^\[user\]' "$conf"; then
        printf '[user]\ndefault=%s\n' "$user" >> "$conf"
    fi
    say "system ready"
}

delfin_step() {
    local repo=$1 commit=$2 head
    [[ "$repo" =~ ^https://[A-Za-z0-9./_-]+$ ]] || { say "invalid repository"; exit 2; }
    [[ "$commit" =~ ^[0-9a-f]{40}$ ]] || { say "invalid commit"; exit 2; }
    [ "$(id -u)" != 0 ] || { say "the DELFIN step runs as the normal user, not root"; exit 2; }
    cd "$HOME"
    if [ ! -d "$HOME/DELFIN/.git" ]; then
        # Pinned: the commit this launcher was released with. Its id is the
        # checksum of everything in it, so what is installed is exactly
        # that. main stays the tracked branch, so "install.sh --update"
        # moves on from there later.
        say "downloading DELFIN ${commit:0:12}"
        rm -rf "$HOME/DELFIN.partial"   # left by a run that was cut off
        git clone -q --branch main "$repo" "$HOME/DELFIN.partial"
        git -C "$HOME/DELFIN.partial" checkout -q -B main "$commit"
        head=$(git -C "$HOME/DELFIN.partial" rev-parse HEAD)
        [ "$head" = "$commit" ] || { say "the download is not commit $commit"; rm -rf "$HOME/DELFIN.partial"; exit 5; }
        mv "$HOME/DELFIN.partial" "$HOME/DELFIN"
    else
        # A checkout that is there stays at its revision (it may have been
        # updated on purpose).
        say "DELFIN checkout found at $(git -C "$HOME/DELFIN" rev-parse --short=12 HEAD 2>/dev/null || echo '?')"
    fi
    # Always run: install.sh keeps what is there and retries what failed.
    # standard = DELFIN, agent, xtb, g-xTB, CREST, MOPAC, DFTB+, sTDA.
    # No OpenMPI: only ORCA uses it, and ORCA cannot be fetched for the
    # person (licence). "install.sh --orca <tarball>" later adds both, with
    # the OpenMPI that ORCA build names.
    say "installing DELFIN (this takes a while)"
    bash "$HOME/DELFIN/install.sh" --profile standard --no-orca --no-openmpi
    [ -f "$HOME/.delfin_env.sh" ] || { say "install.sh did not write ~/.delfin_env.sh"; exit 4; }
    echo DELFIN_WSL_READY
}

main() {
    case "${1:-}" in
        system) system_step "${2:-}" ;;
        delfin) delfin_step "${2:-}" "${3:-}" ;;
        *) say "usage: system <user> | delfin <repository> <commit>"; exit 2 ;;
    esac
}

main "$@" </dev/null
