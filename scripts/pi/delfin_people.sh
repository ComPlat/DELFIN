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
#   delfin-people --invite NAME    one file that lets a new person in by
#                                  themselves: ~/invites/DELFIN-NAME.zip
#                                  [--host HOST] [--hours N, default 48]
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

# -- an invitation: one file, and a new person is in --------------------------
# --add-key needs a round trip: the person presses Start, it fails, they send
# a file, the owner runs a command, they press Start again -- and on the
# first real login OpenSSH still asks them to compare a host fingerprint.
#
# --invite NAME builds DELFIN-NAME.zip: the Windows launcher package plus
#   invite.json  the server (user@host) and this machine's host keys, so the
#                first login asks nothing (the key is pinned, not trusted on
#                first use);
#   invite_key   a key made for this invitation only.
# Its line in authorized_keys can do exactly one thing: run --enroll, which
# reads the new PC's own launcher key from stdin, adds it with the launcher's
# restriction, and removes the invitation. restrict (no pty, no forwarding),
# a forced command, and an expiry time: unused, it stops working by itself.
# Whoever holds the file before the person uses it can enroll ONE PC; that is
# what an invitation is, so hand it over like a temporary password.
invite_hours=48
invite() {
    local name="$1" host="${2:-}" hours="${3:-$invite_hours}" package mac_package work key expiry line self out
    if [[ ! "$name" =~ ^[a-z0-9_][a-z0-9_-]{0,31}$ ]]; then
        echo "usage: delfin-people --invite NAME [--host HOST] [--hours N]" >&2
        echo "NAME: a-z, 0-9, _ and - (it labels the PC's key)" >&2
        exit 2
    fi
    [[ "$hours" =~ ^[0-9]{1,4}$ && "$hours" -ge 1 ]] || { echo "--hours: a whole number of hours" >&2; exit 2; }
    if [[ -f "$HOME/.config/delfin/launcher-no-keys.conf" ]]; then
        echo "this account refuses launcher keys (launcher-no-keys.conf)" >&2
        exit 1
    fi
    if [[ -z "$host" ]]; then
        host="$(hostname -I 2>/dev/null | tr ' ' '\n' | grep -E '^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$' | head -1 || true)"
    fi
    [[ "$host" =~ ^[A-Za-z0-9][A-Za-z0-9.-]{0,200}$ ]] || { echo "no usable address; pass --host NAME-OR-IP" >&2; exit 1; }
    package="${DELFIN_LAUNCHER_ZIP:-}"
    if [[ -z "$package" ]]; then
        for candidate in "$HOME/delfin/src/tools/delfin-windows-launcher.zip" \
                         "$(dirname "$(readlink -f "$0")")/../../tools/delfin-windows-launcher.zip"; do
            if [[ -f "$candidate" ]]; then package="$candidate"; break; fi
        done
    fi
    [[ -f "$package" ]] || { echo "launcher package not found; set DELFIN_LAUNCHER_ZIP" >&2; exit 1; }
    # The macOS package travels along when the server has it: the same
    # invitation then works on a PC (windows-launcher/) or a Mac
    # (macos-launcher/), whichever uses it first.
    mac_package="${DELFIN_MAC_LAUNCHER_ZIP:-}"
    if [[ -z "$mac_package" ]]; then
        mac_package="$(dirname "$package")/delfin-macos-launcher.zip"
    fi
    [[ -f "$mac_package" ]] || mac_package=""
    work="$(mktemp -d)"
    key="$work/invite_key"
    ssh-keygen -q -t ed25519 -N '' -C "delfin-invite-$name" -f "$key"
    expiry="$(date -u -d "+$hours hours" +%Y%m%d%H%MZ)"
    # The command is this very script: the one that made the invitation
    # takes it in, wherever it is installed.
    self="$(readlink -f "$0")"
    line="restrict,command=\"$self --enroll $name\",expiry-time=\"$expiry\" $(cut -d' ' -f1-2 "$key.pub") delfin-invite-$name"
    out="$HOME/invites"
    (umask 077 && mkdir -p "$out" "$HOME/.ssh")
    if ! DELFIN_INVITE_TARGET="$(id -un)@$host" DELFIN_INVITE_NAME="$name" DELFIN_INVITE_EXPIRY="$expiry" \
        DELFIN_HOST_KEYS="${DELFIN_HOST_KEYS:-/etc/ssh/ssh_host_*_key.pub}" \
        python3 - "$package" "$key" "$out/DELFIN-$name.zip" "$mac_package" <<'PY'
import glob, json, os, sys, zipfile
package, key, out, mac_package = sys.argv[1:5]
keys = []
for path in sorted(glob.glob(os.environ['DELFIN_HOST_KEYS'])):
    parts = open(path).read().split()
    if len(parts) >= 2:
        keys.append(parts[0] + ' ' + parts[1])
if not keys:
    sys.exit('no host keys found (' + os.environ['DELFIN_HOST_KEYS'] + ')')
invite = {'target': os.environ['DELFIN_INVITE_TARGET'], 'name': os.environ['DELFIN_INVITE_NAME'],
          'expires': os.environ['DELFIN_INVITE_EXPIRY'], 'host_keys': keys}
os.umask(0o077)
with zipfile.ZipFile(package) as source, zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as target:
    for item in source.infolist():
        target.writestr(item, source.read(item))
    target.writestr('windows-launcher/invite.json', json.dumps(invite, indent=2) + '\n')
    target.writestr('windows-launcher/invite_key', open(key).read())
    if mac_package:
        with zipfile.ZipFile(mac_package) as mac:
            for item in mac.infolist():
                target.writestr(item, mac.read(item))
        target.writestr('macos-launcher/invite.json', json.dumps(invite, indent=2) + '\n')
        target.writestr('macos-launcher/invite_key', open(key).read())
PY
    then
        rm -rf "$work"
        exit 1
    fi
    rm -rf "$work"
    (umask 077 && printf '%s\n' "$line" >> "$HOME/.ssh/authorized_keys")
    echo "Invitation for $name: $out/DELFIN-$name.zip"
    echo "  Server: $(id -un)@$host. Valid once, until $expiry (UTC)."
    echo "  Give the file to $name like a temporary password. They extract it,"
    echo "  double-click Install.cmd, and press Start in DELFIN. Nothing else."
    if [[ -n "$mac_package" ]]; then
        echo "  On a Mac: macos-launcher/Install.command, then choose the server in DELFIN."
    else
        echo "  (No macOS package here; put delfin-macos-launcher.zip beside the Windows one for Macs.)"
    fi
    echo "  Withdraw it: delete the delfin-invite-$name line (delfin-people --keys lists it)."
}

# The forced command of an invitation's key (above): never run by hand.
enroll() {
    local name="$1" line body file tmp
    [[ "$name" =~ ^[a-z0-9_][a-z0-9_-]{0,31}$ ]] || { echo "invalid invitation" >&2; exit 2; }
    file="$HOME/.ssh/authorized_keys"
    IFS= read -r line || true
    line="${line%$'\r'}"
    if ! printf '%s\n' "$line" \
        | grep -Eq '^ssh-ed25519 [A-Za-z0-9+/]+={0,2}( [A-Za-z0-9@._-]{1,64})?$'; then
        echo "not one ed25519 public key line; nothing was added" >&2
        exit 1
    fi
    body="$(printf '%s' "$line" | cut -d' ' -f2)"
    if ! grep -qsF " $body" "$file"; then
        (umask 077 && printf 'restrict,port-forwarding,pty %s delfin-launcher-%s\n' \
            "$(printf '%s' "$line" | cut -d' ' -f1-2)" "$name" >> "$file")
    fi
    # Used once: the invitation's own line goes.
    tmp="$(mktemp "$HOME/.ssh/authorized_keys.XXXXXX")"
    grep -v " delfin-invite-$name\$" "$file" > "$tmp" || true
    chmod 600 "$tmp"
    mv -f "$tmp" "$file"
    echo "DELFIN_ENROLLED"
}

if [[ "${1:-}" == "--invite" ]]; then
    name="${2:-}"; host=""; hours="$invite_hours"
    if [[ $# -ge 2 ]]; then shift 2; else shift $#; fi
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --host) host="${2:-}"; shift 2 ;;
            --hours) hours="${2:-}"; shift 2 ;;
            *) echo "unknown option: $1" >&2; exit 2 ;;
        esac
    done
    invite "$name" "$host" "$hours"
    exit 0
fi
if [[ "${1:-}" == "--enroll" ]]; then enroll "${2:-}"; exit 0; fi

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
