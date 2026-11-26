#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# A box that already runs the one-account layout, moved to the per-person
# variant (and back). Called by pi_bootstrap.sh:
#
#   DELFIN_PI_ACCOUNTS=1 DELFIN_PI_MIGRATE=1 sudo bash pi_bootstrap.sh
#       ~/delfin (a directory) becomes ~/delfin.one-account, untouched;
#       the agent's workspace (work/) and agent.env go to the shared root;
#       the launcher's people mode in the agent account is switched off.
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
#   pi_accounts_migrate.sh to|back|status AGENT_HOME SHARED_ROOT
#
# status prints a loud warning (stderr) while people mode still holds homes
# or keys in the agent account; pi_bootstrap.sh runs it last.
#
# agent.env is copied into the shared root, which is the agent account's
# alone (0700 here, 0750 with the agent's own group after the bootstrap).
# Persons never read it: they source the copy `delfin-accounts publish`
# rebuilds from an allow-list of non-secret variables.
#
# The installation moves by rename (one file system on a Pi: instant, and
# the open files of running services stay valid), so the services keep
# running on the old tree until the deploy restarts them onto the new one.
# ============================================================================

direction="${1:-}"
home="${2:-}"
shared="${3:-}"
[[ "$direction" == to || "$direction" == back || "$direction" == status ]] && [[ "$home" == /* && "$shared" == /* ]] \
    || { echo "usage: pi_accounts_migrate.sh to|back|status AGENT_HOME SHARED_ROOT" >&2; exit 2; }
[[ "$shared" != / && "$home" != / ]] || { echo "error: refusing / as a path" >&2; exit 2; }

ACCOUNTS_CONF="${DELFIN_ACCOUNTS_CONF:-/etc/delfin/accounts.conf}"
current="$home/delfin"
kept="$home/delfin.one-account"
people_conf="$home/.config/delfin/launcher-people.conf"
people_off="$people_conf.off-accounts"

# Who can still get into the agent account the way people mode did: named,
# with what removes them. An access decision is the admin's.
report_access() {
    local keys="$home/.ssh/authorized_keys" names
    names="$(grep -oE ' delfin-(launcher|invite)-[a-z0-9_-]+$' "$keys" 2>/dev/null | sed 's/^ //' | sort -u || true)"
    if [[ -n "$names" ]]; then
        echo "    Keys of the launcher's people mode still log in to the agent account:"
        printf '      %s\n' $names
        echo "    Each person gets an account of their own: sudo delfin-accounts add NAME."
        echo "    Then remove the old line (one per name):"
        echo "      sudo -u $(stat -c %U "$home") sed -i '/ NAME\$/d' $keys"
    fi
    if [[ -d "$home/people" ]] && compgen -G "$home/people/*/" >/dev/null; then
        echo "    People-mode homes, kept as they are (copy what a person still needs):"
        for dir in "$home"/people/*/; do echo "      ${dir%/}"; done
        echo "    Their dashboards keep running on the old installation until stopped."
    fi
}

if [[ "$direction" == status ]]; then
    keys="$(grep -oE ' delfin-(launcher|invite)-[a-z0-9_-]+$' "$home/.ssh/authorized_keys" 2>/dev/null | sed 's/^ //' | sort -u || true)"
    homes=""
    [[ -d "$home/people" ]] && compgen -G "$home/people/*/" >/dev/null && homes=1
    if [[ -n "$keys" || -n "$homes" ]]; then
        {
            echo
            echo "################################################################"
            echo "WARNING: people mode still lives in the agent account ($home)."
            echo "  Its persons share ONE Unix user with the agent; the per-person"
            echo "  accounts protect nobody until this is moved and removed."
            report_access
            echo "################################################################"
        } >&2
    fi
    exit 0
fi

if [[ "$direction" == to ]]; then
    if [[ -L "$current" ]]; then
        echo "==> $current is already a link ($(readlink "$current")): nothing to move"
        exit 0
    fi
    [[ -d "$current" ]] || { echo "==> no $current: nothing to move (a fresh box)"; exit 0; }
    [[ ! -e "$kept" ]] || { echo "error: $kept exists already: an earlier move? Nothing was changed." >&2; exit 1; }
    if [[ -e "$shared" ]] && [[ -n "$(ls -A "$shared" 2>/dev/null)" ]]; then
        echo "error: $shared exists and is not empty. Nothing was changed." >&2
        exit 1
    fi
    echo "==> moving the one-account installation aside: $kept"
    # Closed from the first moment (a plain mkdir under root's umask 022
    # would leave the agent's workspace world-readable until the bootstrap
    # restricts it); the bootstrap then sets the intended group and mode.
    install -d -m 700 "$shared"
    mv "$current" "$kept"
    if [[ -d "$kept/work" ]]; then
        # The agent's workspace is data, not installation: it moves along
        # (a rename), and moves back on `back`.
        mv "$kept/work" "$shared/work"
        echo "==> agent workspace moved: $shared/work"
    fi
    if [[ -f "$kept/agent.env" ]]; then
        cp -p "$kept/agent.env" "$shared/agent.env"
        echo "==> agent.env copied"
    fi
    if [[ -f "$people_conf" ]]; then
        mv "$people_conf" "$people_off"
        echo "==> launcher people mode in this account off (was: $people_conf)"
    fi
    printf '%s\n' "$shared" > "$kept/.moved-to-shared-root"
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
if [[ -d "$shared/work" && ! -e "$kept/work" ]]; then
    mv "$shared/work" "$kept/work"
    echo "==> agent workspace moved back: $kept/work"
fi
rm -f "$current"
rm -f "$kept/.moved-to-shared-root"
mv "$kept" "$current"
echo "==> one-account installation back in $current"
if [[ -f "$people_off" && ! -e "$people_conf" ]]; then
    mv "$people_off" "$people_conf"
    echo "==> launcher people mode on again"
fi
rm -f "$ACCOUNTS_CONF"
echo "    Kept: $shared (remove by hand once the box runs again)."
if command -v getent >/dev/null && [[ -n "$(getent group delfin-people 2>/dev/null | cut -d: -f4)" ]]; then
    echo "    Person accounts still exist: sudo delfin-accounts list / remove NAME --keep-home"
fi
echo "    Next, from your PC: scripts/pi/deploy_pi.sh <agent>@<pi> (starters back to $current)"
