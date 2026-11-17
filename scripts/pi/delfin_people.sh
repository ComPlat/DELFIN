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
#
# They are one Unix user: this is an overview, not a boundary between them.
# ============================================================================

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
