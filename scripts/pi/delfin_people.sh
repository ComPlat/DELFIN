#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# delfin-people: who works in this account (Windows launcher, people mode).
#
# Several people on one server account each get ~/people/<name> as their
# home (tools/windows-launcher/remote_bootstrap.sh). This lists them, read
# only: dashboard (port, running), tmux sessions, DELFIN processes. Installed
# by pi_install.sh as ~/.local/bin/delfin-people.
#
#   delfin-people                  overview
#   delfin-people --attach NAME    watch that person's dashboard tmux session
#   delfin-people --keys           the launcher keys that may log in
#   delfin-people --add-key NAME   let a new PC in (public key on stdin)
#
# They are one Unix user: this is an overview, not a boundary between them.
# ============================================================================

# -- letting a new person in -------------------------------------------------
# On a server that takes passwords, the launcher registers a PC by itself:
# the first login is with the password and remote_bootstrap.sh appends the
# key. Where PasswordAuthentication is off -- as on a hardened box -- that
# first login cannot happen, so a new person is locked out by a key that is
# only added once they are in. Measured on the lab Pi on 2026-11-19.
#
# --add-key is that one step, with the restriction the bootstrap applies.
# Appending by hand is what puts an unrestricted key in authorized_keys: a
# plain `ssh-ed25519 ...` line carries agent forwarding, X11 and ~/.ssh/rc,
# none of which the launcher uses.
add_key() {
    local name="$1" line body
    # The bootstrap's own rule, leading dash included: a name is used to
    # build a path and ends up in a key comment, and `-rf` is a name there
    # only until something passes it to a command.
    if [[ ! "$name" =~ ^[a-z0-9_][a-z0-9_-]{0,31}$ ]]; then
        echo "usage: delfin-people --add-key NAME < key.pub" >&2
        echo "NAME is the Windows user name: a-z, 0-9, _ and -" >&2
        exit 2
    fi
    if [[ -f "$HOME/.config/delfin/launcher-no-keys.conf" ]]; then
        echo "this account refuses launcher keys (launcher-no-keys.conf)" >&2
        exit 1
    fi
    IFS= read -r line || true
    line="${line%$'\r'}"
    # The bootstrap's own check, so one route cannot accept what the other
    # refuses. One line only: a second would be appended unchecked.
    if ! printf '%s\n' "$line" \
        | grep -Eq '^ssh-ed25519 [A-Za-z0-9+/]+={0,2}( [A-Za-z0-9@._-]{1,64})?$'; then
        echo "not one ed25519 public key line; nothing was added" >&2
        exit 1
    fi
    body="$(printf '%s' "$line" | cut -d' ' -f2)"
    if grep -qsF " $body" "$HOME/.ssh/authorized_keys"; then
        echo "that key may already log in; nothing was added"
        exit 0
    fi
    (umask 077 && mkdir -p "$HOME/.ssh" \
        && printf 'restrict,port-forwarding,pty %s delfin-launcher-%s\n' \
            "$(printf '%s' "$line" | cut -d' ' -f1-2)" "$name" \
            >> "$HOME/.ssh/authorized_keys")
    echo "added for $name. Revoke by deleting its delfin-launcher-$name line."
}

if [[ "${1:-}" == "--add-key" ]]; then add_key "${2:-}"; exit 0; fi

if [[ "${1:-}" == "--keys" ]]; then
    # Who may log in without a password, and under which restriction.
    if [[ ! -f "$HOME/.ssh/authorized_keys" ]]; then echo "no authorized_keys"; exit 0; fi
    printf '%-24s %s\n' WHO RESTRICTION
    while IFS= read -r line; do
        [[ "$line" == \#* || -z "$line" ]] && continue
        who="${line##* }"
        case "$line" in
            restrict,*) how="${line%% *}" ;;
            *) how="NONE - full access" ;;
        esac
        printf '%-24s %s\n' "$who" "$how"
    done < "$HOME/.ssh/authorized_keys"
    exit 0
fi

account_home="$(getent passwd "$(id -un)" | cut -d: -f6)"
base="$account_home/people"

sessions_of() {
    tmux ls -F '#S' 2>/dev/null | grep -E "^delfin-win-$1-" || true
}

if [[ "${1:-}" == "--attach" ]]; then
    name="${2:-}"
    [[ "$name" =~ ^[a-z0-9_-]{1,32}$ ]] || { echo "usage: delfin-people --attach NAME" >&2; exit 2; }
    session="$(sessions_of "$name" | head -1)"
    [[ -n "$session" ]] || { echo "no tmux session of $name" >&2; exit 1; }
    exec tmux attach -t "=$session"
fi

if [[ ! -d "$base" ]]; then
    echo "nobody yet: no $base (people mode: ~/.config/delfin/launcher-people.conf)"
    exit 0
fi

printf '%-16s %-10s %-6s %-6s %s\n' PERSON DASHBOARD PORT PROCS TMUX
for dir in "$base"/*/; do
    [[ -d "$dir" ]] || continue
    home="${dir%/}"
    name="$(basename "$home")"
    state="-"; port="-"
    record="$home/.delfin/dashboard_here.json"
    if [[ -f "$record" ]]; then
        port="$(jq -r '.port // "-"' "$record" 2>/dev/null || echo "?")"
        pid="$(jq -r '.pid // empty' "$record" 2>/dev/null || true)"
        if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then state="running"; else state="stopped"; fi
    fi
    procs=0
    for p in /proc/[0-9]*; do
        [[ -O "$p" ]] || continue
        # Caged processes (bubblewrap, no-new-privileges) keep their
        # environment unreadable even to their own user: skipped quietly.
        if { tr '\0' '\n' < "$p/environ"; } 2>/dev/null | grep -qx "HOME=$home" \
            && { tr '\0' ' ' < "$p/cmdline"; } 2>/dev/null | grep -q delfin; then
            procs=$((procs + 1))
        fi
    done
    tmux_sessions="$(sessions_of "$name" | tr '\n' ' ')"
    printf '%-16s %-10s %-6s %-6s %s\n' "$name" "$state" "$port" "$procs" "${tmux_sessions:--}"
done
