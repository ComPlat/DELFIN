#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# A box that already runs the one-account layout, moved to the per-person
# variant (and back). Called by pi_bootstrap.sh:
#
#   DELFIN_PI_ACCOUNTS=1 DELFIN_PI_MIGRATE=1 sudo bash pi_bootstrap.sh
#       the shared root is created closed (root, 0700); ~/delfin (a
#       directory) becomes ~/delfin.one-account, untouched; the agent's
#       workspace (work/) becomes ~/delfin-work (root's link $SHARED/work
#       leads there); agent.env becomes root's copy in the shared root; the
#       launcher's people mode in the agent account is switched off; what of
#       the persons is still in the agent account is warned about, loudly.
#       The bootstrap then sets up the shared root and the link as on a
#       fresh box; the next deploy_pi.sh builds the installation there and
#       restarts the running services onto it.
#   DELFIN_PI_MIGRATE=back sudo bash pi_bootstrap.sh
#       the same steps in reverse; the next deploy_pi.sh points the
#       starters at ~/delfin again.
#
# No data is deleted in either direction (back removes only the link and
# /etc/delfin/accounts.conf, which the bootstrap wrote). Who may log in --
# launcher keys and invitations in the agent account, the person accounts
# -- is listed with the command that changes it, never changed here.
#
#   pi_accounts_migrate.sh to|back AGENT_HOME SHARED_ROOT
#
# The installation moves by rename (one file system on a Pi: instant, and
# the open files of running services stay valid), so the services keep
# running on the old tree until the deploy restarts them onto the new one.
# Root touches the agent account's files only with that account's rights
# (runuser): a link it placed in its home never leads root anywhere.
# ============================================================================

direction="${1:-}"
home="${2:-}"
shared="${3:-}"
[[ "$direction" == to || "$direction" == back ]] && [[ "$home" == /* && "$shared" == /* ]] \
    || { echo "usage: pi_accounts_migrate.sh to|back AGENT_HOME SHARED_ROOT" >&2; exit 2; }
[[ "$shared" != / && "$home" != / ]] || { echo "error: refusing / as a path" >&2; exit 2; }

ACCOUNTS_CONF="${DELFIN_ACCOUNTS_CONF:-/etc/delfin/accounts.conf}"
current="$home/delfin"
kept="$home/delfin.one-account"
work="$home/delfin-work"
people_conf="$home/.config/delfin/launcher-people.conf"
people_off="$people_conf.off-accounts"
# The agent account: everything in its home is read and moved with its own
# rights where a link of its could lead root elsewhere.
owner="$(stat -c %U "$home")"
as_owner() {
    if [[ "$(id -u)" == 0 && "$owner" != root ]]; then
        runuser -u "$owner" -- "$@"
    else
        "$@"
    fi
}

# Who can still get into the agent account the way people mode did, and what
# of the persons still lies in it: named, loudly, with what changes it. An
# access decision is the admin's; this script removes nothing.
report_access() {
    local keys="$home/.ssh/authorized_keys" names homes="" dir
    names="$(as_owner grep -oE ' delfin-(launcher|invite)-[a-z0-9_-]+$' "$keys" 2>/dev/null | sed 's/^ //' | sort -u || true)"
    if [[ -d "$home/people" ]] && compgen -G "$home/people/*/" >/dev/null; then
        for dir in "$home"/people/*/; do homes+="${dir%/} "; done
    fi
    [[ -n "$names" || -n "$homes" ]] || return 0
    echo
    echo "################################################################################"
    echo "WARNING: the agent account still holds the persons of the old layout."
    if [[ -n "$names" ]]; then
        echo "  These keys still log in to the agent account ($owner):"
        printf '    %s\n' $names
        echo "  Give each person an account (sudo delfin-accounts add NAME), then remove"
        echo "  the old line, one per name:"
        echo "    sudo -u $owner sed -i '/ NAME\$/d' $keys"
    fi
    if [[ -n "$homes" ]]; then
        echo "  People-mode homes with the persons' data, still in the agent account:"
        printf '    %s\n' $homes
        echo "  Copy what a person keeps into their new account, then remove the old one:"
        echo "    sudo cp -a $home/people/OLD/. ~NAME/ && sudo chown -R NAME:NAME ~NAME"
        echo "  Until then the agent account can read and change those files, and their"
        echo "  dashboards keep running on the old installation until stopped."
    fi
    echo "################################################################################"
}

if [[ "$direction" == to ]]; then
    if [[ -L "$current" ]]; then
        echo "==> $current is already a link ($(readlink "$current")): nothing to move"
        report_access
        exit 0
    fi
    [[ -d "$current" ]] || { echo "==> no $current: nothing to move (a fresh box)"; exit 0; }
    [[ ! -e "$kept" && ! -L "$kept" ]] || { echo "error: $kept exists already: an earlier move? Nothing was changed." >&2; exit 1; }
    [[ ! -e "$work" && ! -L "$work" ]] || { echo "error: $work exists already. Nothing was changed." >&2; exit 1; }
    if [[ -e "$shared" || -L "$shared" ]] && [[ -n "$(ls -A "$shared" 2>/dev/null)" || -L "$shared" ]]; then
        echo "error: $shared exists and is not an empty directory. Nothing was changed." >&2
        exit 1
    fi
    # The shared root is created closed (root, 0700) before anything moves:
    # the bootstrap opens it to the readers' group once it is set up.
    install -d -m 700 "$shared"
    echo "==> moving the one-account installation aside: $kept"
    as_owner mv "$current" "$kept"
    if [[ -d "$kept/work" && ! -L "$kept/work" ]]; then
        # The agent's workspace is its data, not the installation: it stays in
        # its home (root's link $shared/work leads there), and goes back on `back`.
        as_owner mv "$kept/work" "$work"
        echo "==> agent workspace: $work"
    fi
    if [[ -f "$kept/agent.env" ]]; then
        # Read with the agent account's rights (a link of its reads as
        # nothing it could not read anyway), written as root's own file.
        if as_owner cat -- "$kept/agent.env" > "$shared/agent.env.tmp"; then
            chmod 644 "$shared/agent.env.tmp"
            mv -f "$shared/agent.env.tmp" "$shared/agent.env"
            echo "==> agent.env: root's copy in $shared"
        else
            rm -f "$shared/agent.env.tmp"
        fi
    fi
    if [[ -f "$people_conf" ]]; then
        as_owner mv "$people_conf" "$people_off"
        echo "==> launcher people mode in this account off (was: $people_conf)"
    fi
    printf '%s\n' "$shared" | as_owner tee "$kept/.moved-to-shared-root" >/dev/null
    report_access
    echo "    Back to the one-account layout: DELFIN_PI_MIGRATE=back sudo bash pi_bootstrap.sh"
    exit 0
fi

# back
[[ -d "$kept" && ! -L "$kept" ]] || { echo "error: no $kept to go back to. Nothing was changed." >&2; exit 1; }
if [[ -e "$current" && ! -L "$current" ]]; then
    echo "error: $current is a directory, not the shared root's link. Nothing was changed." >&2
    exit 1
fi
if [[ -L "$current" && "$(readlink -f "$current")" != "$(readlink -f "$shared")" ]]; then
    echo "error: $current points to $(readlink "$current"), not $shared. Nothing was changed." >&2
    exit 1
fi
if [[ -d "$work" && ! -L "$work" && ! -e "$kept/work" ]]; then
    as_owner mv "$work" "$kept/work"
    echo "==> agent workspace moved back: $kept/work"
fi
as_owner rm -f "$current"
as_owner rm -f "$kept/.moved-to-shared-root"
as_owner mv "$kept" "$current"
echo "==> one-account installation back in $current"
if [[ -f "$people_off" && ! -e "$people_conf" ]]; then
    as_owner mv "$people_off" "$people_conf"
    echo "==> launcher people mode on again"
fi
rm -f "$ACCOUNTS_CONF"
echo "    Kept: $shared (remove by hand once the box runs again)."
if command -v getent >/dev/null && [[ -n "$(getent group delfin-people 2>/dev/null | cut -d: -f4)" ]]; then
    echo "    Person accounts still exist: sudo delfin-accounts list / remove NAME --keep-home"
fi
echo "    Next, from your PC: scripts/pi/deploy_pi.sh <agent>@<pi> (starters back to $current)"
