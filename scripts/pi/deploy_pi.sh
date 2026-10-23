#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Push the current DELFIN working tree to a Raspberry Pi agent test box.
# Runs on your PC (Git Bash / WSL / Linux / macOS). Needs: git, tar, ssh.
#
#   deploy_pi.sh --bootstrap pi@raspi    # once: set up the Pi (asks sudo pw)
#   deploy_pi.sh delfin@raspi            # deploy working tree (incl. uncommitted)
#   deploy_pi.sh --test delfin@raspi     # deploy + run pytest -m "not slow"
#   deploy_pi.sh --rollback delfin@raspi # back to the previous deploy
#   deploy_pi.sh --status delfin@raspi   # show deployed version
#
# Deploys exactly what git tracks plus untracked-but-not-ignored files, so
# .gitignore'd secrets, venvs and caches never leave this machine.
# ============================================================================

usage() { sed -n '7,12p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }

action=install
run_tests=""
case "${1:-}" in
    --bootstrap) action=bootstrap; shift ;;
    --rollback)  action=rollback;  shift ;;
    --status)    action=status;    shift ;;
    --test)      run_tests="--test"; shift ;;
    -h|--help|"") usage ;;
esac
HOST="${1:-}"
[[ -n "$HOST" ]] || usage

# Git Bash's ssh cannot resolve mDNS (*.local); Windows itself can.
if [[ "$HOST" == *.local && "$(uname -s)" == MINGW* ]]; then
    name="${HOST#*@}"
    ip="$(powershell.exe -NoProfile -Command \
        "(Resolve-DnsName -Type A '$name' -ErrorAction SilentlyContinue | Select-Object -First 1).IPAddress" | tr -d '\r')"
    [[ -n "$ip" ]] || { echo "error: cannot resolve $name" >&2; exit 1; }
    [[ "$HOST" == *@* ]] && HOST="${HOST%@*}@$ip" || HOST="$ip"
    echo "==> $name -> $ip"
fi

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(git -C "$HERE" rev-parse --show-toplevel)"

push_helper() {
    tr -d '\r' < "$HERE/pi_install.sh" \
        | ssh "$HOST" 'mkdir -p ~/delfin && cat > ~/delfin/pi_install.sh && chmod +x ~/delfin/pi_install.sh'
}

case "$action" in
bootstrap)
    tr -d '\r' < "$HERE/pi_bootstrap.sh" | ssh "$HOST" 'cat > /tmp/pi_bootstrap.sh'
    ssh -t "$HOST" 'sudo bash /tmp/pi_bootstrap.sh && rm -f /tmp/pi_bootstrap.sh'
    ;;
install)
    cd "$REPO"
    # Snapshot the working tree into a throwaway index -> tree, then archive it
    # with LF endings (a core.autocrlf=true checkout would otherwise ship CRLF
    # shell scripts that break on Linux).
    tmp_index="$(mktemp)"
    trap 'rm -f "$tmp_index"' EXIT
    cp "$(git rev-parse --git-path index)" "$tmp_index"
    GIT_INDEX_FILE="$tmp_index" git -c core.safecrlf=false add -A
    tree="$(GIT_INDEX_FILE="$tmp_index" git write-tree)"
    version="$(git rev-parse --short HEAD)"
    [[ "$tree" == "$(git rev-parse 'HEAD^{tree}')" ]] || version="$version-dirty"
    echo "==> deploying $version to $HOST"
    push_helper
    git -c core.autocrlf=false -c core.eol=lf archive --format=tar.gz "$tree" \
        | ssh "$HOST" "bash ~/delfin/pi_install.sh install '$version' $run_tests"
    ;;
rollback|status)
    push_helper
    ssh "$HOST" "bash ~/delfin/pi_install.sh $action"
    ;;
esac
