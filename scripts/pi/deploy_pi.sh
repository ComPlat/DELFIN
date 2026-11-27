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
#   deploy_pi.sh --landlock pi@raspi     # once, as the admin: Landlock kernel (~1 h, reboots)
#   deploy_pi.sh --shared pi@raspi       # one account per person: deploy the shared
#                                        # installation as the admin (asks sudo pw);
#                                        # --shared-test, --shared-rollback, --shared-status
#
# Deploys exactly what git tracks plus untracked-but-not-ignored files, so
# .gitignore'd secrets, venvs and caches never leave this machine.
# ============================================================================

usage() { sed -n '7,16p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }

action=install
run_tests=""
case "${1:-}" in
    --bootstrap) action=bootstrap; shift ;;
    --rollback)  action=rollback;  shift ;;
    --status)    action=status;    shift ;;
    --landlock)  action=landlock;  shift ;;
    --test)      run_tests="--test"; shift ;;
    --shared)          action=shared; shift ;;
    --shared-test)     action=shared; run_tests="--test"; shift ;;
    --shared-rollback) action=shared-rollback; shift ;;
    --shared-status)   action=shared-status; shift ;;
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

# DELFIN_PI_* set where this runs (DELFIN_PI_PYSCF=0, DELFIN_PI_SKIP_EXTRAS=...)
# reach pi_install.sh on the Pi, as `env NAME='value'` words. Only names of
# that form and values of letters, digits, space and ._,- pass; anything
# else is refused, so no value can carry shell syntax to the remote side.
pi_env() {
    local name value words=""
    for name in $(compgen -v DELFIN_PI_); do
        [[ "$name" =~ ^DELFIN_PI_[A-Z0-9_]+$ ]] || continue
        value="${!name}"
        if [[ ! "$value" =~ ^[A-Za-z0-9\ ._,-]*$ ]]; then
            echo "error: $name has characters a deploy does not forward: $value" >&2
            return 1
        fi
        words+=" $name='$value'"
    done
    [[ -z "$words" ]] || printf 'env%s ' "$words"
}

# The receiving script, staged in the agent account's own cache (not in
# ~/delfin: with one account per person that is root's shared installation,
# and pi_install.sh refuses to deploy it as the agent account).
push_helper() {
    tr -d '\r' < "$HERE/pi_install.sh" \
        | ssh "$HOST" 'mkdir -p ~/.cache/delfin-deploy && cat > ~/.cache/delfin-deploy/pi_install.sh && chmod +x ~/.cache/delfin-deploy/pi_install.sh'
}

# The working tree as a tree object (uncommitted changes included), and the
# version string a deploy records.
snapshot() {
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
}

case "$action" in
bootstrap)
    # A name from mktemp: created 0600 by the login user, and in the sticky
    # temp directory nobody else can put another file under it before sudo
    # runs it. A fixed name could be created or swapped by any local user.
    # The scripts beside it go along: one account per person
    # (DELFIN_PI_ACCOUNTS=1) installs root-owned copies of them.
    stage="$(mktemp -d)"
    trap 'rm -rf "$stage"' EXIT
    for f in pi_bootstrap.sh delfin_accounts.sh delfin_people.sh delfin_deploy.sh; do
        tr -d '\r' < "$HERE/$f" > "$stage/$f"
    done
    remote="$(tar -C "$stage" -cf - pi_bootstrap.sh delfin_accounts.sh delfin_people.sh delfin_deploy.sh \
        | ssh "$HOST" 'd="$(mktemp -d)" && tar -xf - -C "$d" && echo "$d"')"
    [[ "$remote" == /* ]] || { echo "error: upload of pi_bootstrap.sh failed" >&2; exit 1; }
    # DELFIN_PI_* reach the bootstrap too (sudo would drop them).
    ssh -t "$HOST" "sudo $(pi_env)bash '$remote/pi_bootstrap.sh'; rc=\$?; rm -rf '$remote'; exit \$rc"
    ;;
install)
    snapshot
    echo "==> deploying $version to $HOST"
    remote_env="$(pi_env)"
    push_helper
    git -c core.autocrlf=false -c core.eol=lf archive --format=tar.gz "$tree" \
        | ssh "$HOST" "${remote_env}bash ~/.cache/delfin-deploy/pi_install.sh install '$version' $run_tests"
    ;;
shared)
    # One account per person: the admin uploads the tarball into a mktemp
    # directory of its own and root's delfin-deploy builds it as the install
    # account. The agent account is not involved in the shared installation.
    snapshot
    echo "==> deploying $version to the shared installation on $HOST (as the admin)"
    remote_env="$(pi_env)"
    remote="$(git -c core.autocrlf=false -c core.eol=lf archive --format=tar.gz "$tree" \
        | ssh "$HOST" 'd="$(mktemp -d)" && cat > "$d/src.tgz" && echo "$d"')"
    [[ "$remote" == /* ]] || { echo "error: upload failed" >&2; exit 1; }
    ssh -t "$HOST" "sudo ${remote_env}/usr/local/sbin/delfin-deploy install '$version' '$remote/src.tgz' $run_tests; rc=\$?; rm -rf '$remote'; exit \$rc"
    ;;
shared-rollback|shared-status)
    ssh -t "$HOST" "sudo $(pi_env)/usr/local/sbin/delfin-deploy ${action#shared-}"
    ;;
landlock)
    # As the admin (the agent user has no root), from the agent's deployed
    # tree; tryboot falls back to the stock kernel by itself.
    k="~${AGENT_USER:-delfin}/delfin/src/scripts/pi/pi_kernel_landlock.sh"
    ssh -t "$HOST" "sudo bash $k build && sudo bash $k try"
    ;;
rollback|status)
    remote_env="$(pi_env)"
    push_helper
    ssh "$HOST" "${remote_env}bash ~/.cache/delfin-deploy/pi_install.sh $action"
    ;;
esac
