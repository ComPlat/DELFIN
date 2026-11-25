#!/bin/bash
# DELFIN for macOS: the Mac counterpart of tools/windows-launcher/DELFIN.ps1.
#
#   DELFIN.sh app                  the app's menu (DELFIN.app runs this)
#   DELFIN.sh connect ID           dashboard over SSH, in a Terminal window
#   DELFIN.sh local                the dashboard of DELFIN on this Mac
#   DELFIN.sh setup-local          set DELFIN up on this Mac (setup_local.sh)
#   DELFIN.sh install [options]    install the app from this package
#   DELFIN.sh uninstall [options]  remove it (saved connections on request)
#   DELFIN.sh validate             check the package against SHA256SUMS.txt
#
# Written for the /bin/bash every Mac has (3.2): no associative arrays, no
# ${var,,}, no mapfile, no `set -u` (an empty array is "unbound" there).
# Nothing here needs Python, Xcode or Homebrew on the Mac: ssh, scp,
# ssh-keygen, shasum, lsof, curl, open and osascript are part of macOS.
#
# The server half is the Windows launcher's, byte for byte
# (remote_bootstrap.sh, remote_launcher.py, copied in by build_package.py):
# the same check command, the same arguments, the same
# ~/.cache/delfin-launcher/<version>/ folder on the server, so a Windows PC
# and a Mac share one copy there.
#
# Security, as in the Windows launcher: StrictHostKeyChecking=ask,
# IdentitiesOnly with the launcher's own key, no agent or X11 forwarding, no
# forwarding but the one dashboard tunnel bound to 127.0.0.1, no stored
# passwords (ssh asks in the Terminal window itself; this script never sees
# them), and the browser gets the token only once the tunnel's listener is
# verified to be this connection's own ssh.

# Fixed system tools: a Homebrew or user ssh earlier on PATH is not taken
# silently (the Windows launcher pins System32\OpenSSH the same way).
PATH=/usr/bin:/bin:/usr/sbin:/sbin
export PATH
umask 077

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
# HERE is DELFIN.app/Contents/Resources; the package (Install.command,
# SHA256SUMS.txt, an invitation) is the folder that holds DELFIN.app.
BUNDLE=$(cd "$HERE/../.." && pwd -P)
PKG=$(dirname "$BUNDLE")
STORE=${DELFIN_LAUNCHER_STORE:-"$HOME/Library/Application Support/DELFIN Launcher"}
APPS=${DELFIN_LAUNCHER_APPS:-"$HOME/Applications"}
SSH_DIR=${DELFIN_LAUNCHER_SSH_DIR:-"$HOME/.ssh"}
DESKTOP=${DELFIN_LAUNCHER_DESKTOP:-"$HOME/Desktop"}
SSH=/usr/bin/ssh
SCP=/usr/bin/scp
KEYGEN=/usr/bin/ssh-keygen

say() { printf '%s\n' "$*"; }
warn() { printf '%s\n' "$*" >&2; }
die() { printf 'DELFIN: %s\n' "$*" >&2; exit 1; }

# -- values and checks (shared with the Windows launcher's rules) ------------

# POSIX single-quote escaping, as DELFIN.ps1's Shell-Quote.
shq() { printf "'%s'" "$(printf '%s' "$1" | sed "s/'/'\\\\''/g")"; }

is_uuid() { [[ $1 =~ ^[a-fA-F0-9]{8}-[a-fA-F0-9]{4}-[a-fA-F0-9]{4}-[a-fA-F0-9]{4}-[a-fA-F0-9]{12}$ ]]; }
new_uuid() { uuidgen 2>/dev/null | tr 'A-F' 'a-f' || cat /proc/sys/kernel/random/uuid; }

valid_target() { [[ $1 =~ ^[a-zA-Z0-9][a-zA-Z0-9_.@-]{0,200}$ ]]; }

# Working directory and DELFIN location are server paths: ~, ~/..., or plain
# path characters. Surrounding quotes, as people paste them, are dropped.
clean_path() {
    local path=$1
    path=${path#"${path%%[![:space:]]*}"}
    path=${path%"${path##*[![:space:]]}"}
    if [ ${#path} -ge 2 ]; then
        case $path in \"*\"|\'*\') path=${path:1:${#path}-2} ;; esac
    fi
    if [ -n "$path" ] && ! [[ $path =~ ^(~(/.*)?|[A-Za-z0-9_./\ -]+)$ ]]; then return 1; fi
    printf '%s' "$path"
}

valid_port() { [[ $1 =~ ^[0-9]{1,5}$ ]] && { [ "$1" -eq 0 ] || { [ "$1" -ge 1024 ] && [ "$1" -le 65535 ]; }; }; }

# Who sits at this Mac account; a server with several people on one account
# gives each a home of their own (remote_bootstrap.sh), others ignore it.
person_name() {
    local name
    name=$(id -un | tr 'A-Z' 'a-z' | sed 's/[^a-z0-9_-]/-/g; s/^-*//; s/-*$//')
    name=${name:0:32}
    name=$(printf '%s' "$name" | sed 's/-*$//')
    printf '%s' "$name"
}

sha256_of() {
    if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | awk '{print $1}'
    else sha256sum "$1" | awk '{print $1}'; fi
}

# -- profiles: one file per saved connection, key=value, read as data -------

profiles_dir() { printf '%s/profiles' "$STORE"; }
profile_file() { printf '%s/%s.conf' "$(profiles_dir)" "$1"; }
profile_get() { sed -n "s/^$2=//p" "$1" | head -n 1; }

# Profile ids, newest first (the last one saved or used leads the menu).
profile_ids() {
    local dir file
    dir=$(profiles_dir)
    [ -d "$dir" ] || return 0
    ls -t "$dir" 2>/dev/null | while IFS= read -r file; do
        case $file in *.conf) is_uuid "${file%.conf}" && printf '%s\n' "${file%.conf}" ;; esac
    done
}

# save_profile ID NAME TARGET DIRECTORY EXECUTABLE PORT KEEP WORKING REMEMBER
save_profile() {
    local id=$1 name=$2 target=$3 directory executable
    is_uuid "$id" || die 'Invalid connection ID.'
    valid_target "$target" || die 'Enter an SSH alias or user@hostname (no options).'
    directory=$(clean_path "$4") || die 'Working directory must be a server path: /path, ~/path or a relative path.'
    executable=$(clean_path "$5") || die 'DELFIN location must be a server path: /path, ~/path or a relative path.'
    valid_port "$6" || die 'Starting port must be 0 (automatic) or between 1024 and 65535.'
    case $name in *[![:print:]]*|'') name=$target ;; esac
    mkdir -p "$(profiles_dir)"
    {
        printf 'name=%s\n' "$name"
        printf 'target=%s\n' "$target"
        printf 'directory=%s\n' "$directory"
        printf 'executable=%s\n' "$executable"
        printf 'port=%s\n' "$6"
        printf 'keep=%s\n' "$7"
        printf 'working=%s\n' "$8"
        printf 'remember=%s\n' "$9"
    } > "$(profile_file "$id").tmp"
    mv -f "$(profile_file "$id").tmp" "$(profile_file "$id")"
}

# Loads P_* from a profile and checks it as DELFIN.ps1's Check-Profile does.
load_profile() {
    local file
    is_uuid "$1" || die 'Invalid connection ID.'
    file=$(profile_file "$1")
    [ -f "$file" ] || die 'Connection not found.'
    P_ID=$1
    P_NAME=$(profile_get "$file" name)
    P_TARGET=$(profile_get "$file" target)
    P_PORT=$(profile_get "$file" port); P_PORT=${P_PORT:-0}
    P_KEEP=$(profile_get "$file" keep); P_KEEP=${P_KEEP:-0}
    P_WORKING=$(profile_get "$file" working); P_WORKING=${P_WORKING:-0}
    P_REMEMBER=$(profile_get "$file" remember); P_REMEMBER=${P_REMEMBER:-1}
    valid_target "$P_TARGET" || die 'Enter an SSH alias or user@hostname (no options).'
    P_DIRECTORY=$(clean_path "$(profile_get "$file" directory)") || die 'Working directory must be a server path.'
    P_EXECUTABLE=$(clean_path "$(profile_get "$file" executable)") || die 'DELFIN location must be a server path.'
    valid_port "$P_PORT" || die 'Starting port must be 0 (automatic) or between 1024 and 65535.'
    [ "$P_KEEP" = 1 ] || P_KEEP=0
    [ "$P_WORKING" = 1 ] || P_WORKING=0
    [ "$P_REMEMBER" = 0 ] || P_REMEMBER=1
    [ -n "$P_NAME" ] || P_NAME=$P_TARGET
}

# -- the server command: identical to DELFIN.ps1's Get-RemoteCommand ---------

# Exit 0: both scripts present with these digests. Exit 10 (and
# DELFIN_SCRIPTS_MISSING): missing or different. Exit 11: no folder.
CHECK_SCRIPT='umask 077 && cd && mkdir -p $1 && cd $1 || exit 11; printf '"'"'%s  %s\n'"'"' $2 remote_bootstrap.sh $3 remote_launcher.py | sha256sum -c --status 2>/dev/null || { echo DELFIN_SCRIPTS_MISSING; exit 10; }'

# Sets V_FOLDER, V_BOOTSTRAP, V_LAUNCHER from the packaged server scripts.
script_version() {
    V_BOOTSTRAP=$(sha256_of "$HERE/remote_bootstrap.sh")
    V_LAUNCHER=$(sha256_of "$HERE/remote_launcher.py")
    V_FOLDER=".cache/delfin-launcher/${V_BOOTSTRAP:0:16}${V_LAUNCHER:0:16}"
    [[ $V_FOLDER =~ ^\.cache/delfin-launcher/[a-f0-9]{32}$ && $V_BOOTSTRAP =~ ^[a-f0-9]{64}$ && $V_LAUNCHER =~ ^[a-f0-9]{64}$ ]] \
        || die 'Invalid launcher script version.'
}

script_check_command() {
    printf '%s' "bash -c $(shq "$CHECK_SCRIPT") delfin $V_FOLDER $V_BOOTSTRAP $V_LAUNCHER"
}

# remote_command ROLE TUNNEL PUBLIC_KEY (profile in P_*, version in V_*)
remote_command() {
    local bootstrap=". ~/$V_FOLDER/remote_bootstrap.sh" launcher="~/$V_FOLDER/remote_launcher.py"
    printf '%s' "$(script_check_command) && bash -ilc $(shq "$bootstrap") delfin $P_PORT $(shq "$P_DIRECTORY") $(shq "$P_EXECUTABLE") $(shq "$launcher") $(shq "$1") $(shq "$2") $P_KEEP $P_WORKING $(shq "$(person_name)") $(shq "$3")"
}

# ssh_options KIND LOCAL REMOTE TUNNEL -> SSH_OPTS (DELFIN.ps1's Get-SSHOptions)
ssh_options() {
    local kind=$1 local_port=${2:-0} tunnel=${4:-} policy=ClearAllForwardings=yes private=0
    if [ "$kind" = Dashboard ] && [ "$local_port" -ge 1024 ] && [ "$local_port" -le 65535 ] \
        && [[ $tunnel =~ ^/tmp/delfin-win-tunnel-[a-f0-9-]{36}/dashboard\.sock$ ]]; then
        private=1
        policy=ClearAllForwardings=no
    fi
    SSH_OPTS=()
    [ "$kind" = Dashboard ] && SSH_OPTS=(-tt)
    SSH_OPTS=("${SSH_OPTS[@]}" -o StrictHostKeyChecking=ask -o ForwardAgent=no -o ForwardX11=no -o "$policy"
        -o PermitLocalCommand=no -o ServerAliveInterval=30 -o ServerAliveCountMax=3 -o ControlMaster=no -o ControlPath=none)
    if [ $private = 1 ]; then
        SSH_OPTS=("${SSH_OPTS[@]}" -o GatewayPorts=no -o ExitOnForwardFailure=yes -L "127.0.0.1:$local_port:$tunnel")
    fi
}

# The browser URL the server hands over: only the dashboard on 127.0.0.1,
# with a token, nothing else (DELFIN.ps1's Get-BrowserUrl).
browser_url() {
    local url
    [[ $1 =~ ^[A-Za-z0-9+/=]+$ ]] || return 1
    url=$(printf '%s' "$1" | base64 -d 2>/dev/null || printf '%s' "$1" | base64 -D 2>/dev/null) || return 1
    [[ $url =~ ^http://127\.0\.0\.1:([0-9]{4,5})(/|/voila/render/[A-Za-z0-9%_./~+-]*)\?token=[A-Za-z0-9_-]{32,128}(\&session=[A-Za-z0-9%_.~+-]+)?$ ]] || return 1
    [ "${BASH_REMATCH[1]}" -ge 1024 ] && [ "${BASH_REMATCH[1]}" -le 65535 ] || return 1
    case $url in *..*) return 1 ;; esac
    printf '%s' "$url"
}

# with_port URL PORT: the same URL on this Mac's end of the tunnel.
with_port() { printf '%s' "$1" | sed -E "s#^http://127\.0\.0\.1:[0-9]+#http://127.0.0.1:$2#"; }

port_in_use() { lsof -nP -iTCP:"$1" -sTCP:LISTEN -t >/dev/null 2>&1; }

# A free 127.0.0.1 port: the profile's starting port or the next free one,
# else a random one. ExitOnForwardFailure closes the remaining race.
free_local_port() {
    local port=${1:-0} last tries=0
    if [ "$port" -ne 0 ]; then
        last=$((port + 100)); [ $last -le 65535 ] || last=65535
        while [ "$port" -le $last ]; do port_in_use "$port" || { printf '%s' "$port"; return; }; port=$((port + 1)); done
    fi
    while [ $tries -lt 200 ]; do
        port=$((RANDOM % 30000 + 20000))
        port_in_use "$port" || { printf '%s' "$port"; return; }
        tries=$((tries + 1))
    done
    die 'No local port available.'
}

# Only this connection's own ssh may listen on the tunnel port before the
# token goes to the browser.
listener_is() { lsof -nP -a -p "$2" -iTCP@127.0.0.1:"$1" -sTCP:LISTEN -t >/dev/null 2>&1; }

# DELFIN_LAUNCHER_URL_FILE: the CI's stand-in for a browser (the URL goes
# into that file, 0600, instead of to `open`).
open_url() {
    if [ -n "${DELFIN_LAUNCHER_URL_FILE:-}" ]; then printf '%s\n' "$1" > "$DELFIN_LAUNCHER_URL_FILE"
    elif command -v open >/dev/null 2>&1; then open "$1"
    else xdg-open "$1" >/dev/null 2>&1; fi
}

# -- Remember login: one ed25519 key per Mac account --------------------------

# Sets LOGIN_KEY and LOGIN_PUBLIC, or leaves them empty (the login stays as is).
login_key() {
    LOGIN_KEY=''; LOGIN_PUBLIC=''
    local key="$SSH_DIR/delfin_launcher_ed25519" line
    if [ ! -f "$key" ]; then
        [ -x "$KEYGEN" ] || return 0
        mkdir -p "$SSH_DIR" && chmod 700 "$SSH_DIR"
        "$KEYGEN" -q -t ed25519 -N '' -C "delfin-launcher-$(person_name)" -f "$key" </dev/null || return 0
    fi
    [ -f "$key.pub" ] || return 0
    line=$(tr -d '\r\n' < "$key.pub")
    [[ $line =~ ^ssh-ed25519\ [A-Za-z0-9+/]+={0,2}(\ [A-Za-z0-9@._-]{1,64})?$ ]] || return 0
    LOGIN_KEY=$key; LOGIN_PUBLIC=$line
}

# -- the dashboard connection (runs in a Terminal window) --------------------

copy_scripts() {
    say 'Copying the DELFIN launcher scripts to the server (once per launcher version).'
    say 'Without a remembered login, ssh asks for the password and OTP again for this copy.'
    ssh_options Copy
    "$SCP" -q "${SSH_OPTS[@]}" "${IDENTITY[@]}" "$HERE/remote_bootstrap.sh" "$HERE/remote_launcher.py" "$P_TARGET:$V_FOLDER/" \
        || die 'Could not copy the launcher scripts to the server.'
}

enroll_invitation() {
    local invite="$STORE/invite-$P_ID.key" out
    [ -n "$LOGIN_KEY" ] && [ -f "$invite" ] || return 0
    say 'Registering this Mac with the server (once, no password needed)...'
    ssh_options Copy
    out=$(printf '%s\n' "$LOGIN_PUBLIC" | "$SSH" -T "${SSH_OPTS[@]}" -o BatchMode=yes -o IdentitiesOnly=yes -i "$invite" "$P_TARGET")
    case $out in
        *DELFIN_ENROLLED*) rm -f "$invite"; say 'This Mac is registered. From now on it logs in by itself.' ;;
        *) warn 'The invitation could not be used (expired, used already, or the server is not reachable).'
           warn 'Ask whoever runs the server for a new one. Trying the normal login instead.' ;;
    esac
}

# Waits for this connection's tunnel, then opens the dashboard.
browser_watch() {
    local url=$1 port=$2 owner=$3 deadline=$((SECONDS + 300)) code
    while [ $SECONDS -lt $deadline ]; do
        kill -0 "$owner" 2>/dev/null || return 0
        if listener_is "$port" "$owner"; then
            code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 2 "http://127.0.0.1:$port/login")
            if [ "$code" = 200 ]; then open_url "$url"; return 0; fi
        fi
        sleep 2
    done
    say 'Browser startup timed out after 5 minutes. Check this window and the SSH login.'
}

refused_note() {
    local note="$DESKTOP/DELFIN-key-for-the-server-owner.txt"
    {
        say "The DELFIN launcher could not log in to $P_TARGET."
        say ''
        say 'If you were given no password, this Mac is not known to the server yet.'
        say 'Send THIS FILE to whoever runs the server. It holds the public half of'
        say "this Mac's key - it is not secret, and it is not your password."
        say ''
        say 'They add it with one command:'
        say ''
        say "    delfin-people --add-key $(person_name) < this-file.txt"
        say ''
        say 'Then start DELFIN again. No password will be needed.'
        say ''
        say "$LOGIN_PUBLIC"
    } > "$note" 2>/dev/null || { say 'The login was refused. This Mac key could not be written to the desktop:'; say "$LOGIN_PUBLIC"; return; }
    say ''
    say '  The login was refused.'
    say '  If you were given no password, this Mac is not known to the server.'
    say '  A file is now on your desktop: DELFIN-key-for-the-server-owner.txt'
    say '  Send it to whoever runs the server, then start DELFIN again.'
}

connect_worker() {
    load_profile "$1"
    printf '\033]0;SSH Dashboard - %s\007' "$P_NAME"
    local local_port tunnel remote attempt=0 line url opened=0 missing=0 code=0 closed=0 marker
    # Global, not local: the EXIT trap runs after this function has returned.
    work=''; ssh_pid=''; watch_pid=''
    local_port=$(free_local_port "$P_PORT")
    tunnel="/tmp/delfin-win-tunnel-$(new_uuid)/dashboard.sock"
    IDENTITY=()
    LOGIN_KEY=''; LOGIN_PUBLIC=''
    [ "$P_REMEMBER" = 1 ] && login_key
    # IdentitiesOnly: offer THIS key and no other. Each key the agent or
    # ~/.ssh/id_* adds costs one of the server's MaxAuthTries, and a
    # hardened server then cuts the login before the password prompt.
    [ -n "$LOGIN_KEY" ] && IDENTITY=(-o IdentitiesOnly=yes -i "$LOGIN_KEY")
    enroll_invitation
    say ''
    say '  The first time, you may be asked for the password up to three times. That is normal:'
    say '  one question per connection, not a refusal.'
    [ -n "$LOGIN_KEY" ] && say '  After this once, this Mac is known and no password is asked again.'
    say ''
    say 'Password and OTP are handled by ssh itself; this app never sees them.'
    say 'On first login, compare the host fingerprint with the official server fingerprint.'
    say 'The browser signs in automatically.'
    if [ "$P_KEEP" = 1 ]; then
        say 'Keep session ON: close this window or detach with Ctrl+B, then D. Reconnect: start DELFIN again.'
    else
        say 'Keep session OFF: closing this window stops a dashboard newly started by this connection.'
    fi
    say 'Ctrl+C in the dashboard explicitly stops it, even with Keep session ON.'
    script_version
    remote=$(remote_command dashboard "$tunnel" "$LOGIN_PUBLIC")
    ssh_options Dashboard "$local_port" 0 "$tunnel"
    work=$(mktemp -d "${TMPDIR:-/tmp}/delfin-launcher.XXXXXX") || die 'No temporary folder.'
    trap '[ -n "$watch_pid" ] && kill "$watch_pid" 2>/dev/null; [ -n "$ssh_pid" ] && kill "$ssh_pid" 2>/dev/null; rm -rf "$work"' EXIT
    while [ $attempt -lt 2 ]; do
        attempt=$((attempt + 1)); missing=0; opened=0
        rm -f "$work/out"; mkfifo "$work/out"
        # The terminal stays ssh's input: passwords, tmux keys, Ctrl+C.
        "$SSH" "${SSH_OPTS[@]}" "${IDENTITY[@]}" "$P_TARGET" "$remote" < /dev/tty > "$work/out" &
        ssh_pid=$!
        exec 3< "$work/out"
        while IFS= read -r line <&3; do
            line=${line%$'\r'}
            if [ "$line" = DELFIN_SCRIPTS_MISSING ]; then missing=1; continue; fi
            if [[ $line =~ ^DELFIN_BROWSER:([A-Za-z0-9+/=]+)$ ]]; then
                url=$(browser_url "${BASH_REMATCH[1]}") || { warn 'Invalid browser URL from the server.'; kill "$ssh_pid" 2>/dev/null; break; }
                say "Dashboard on this Mac: 127.0.0.1:$local_port (tunnel to the server)"
                browser_watch "$(with_port "$url" "$local_port")" "$local_port" "$ssh_pid" &
                watch_pid=$!
                url=''
                opened=1
                break
            fi
            printf '%s\r\n' "$line"
        done
        if [ $opened = 1 ]; then
            marker="DELFIN_CONNECTION_CLOSED:$(basename "$(dirname "$tunnel")")"
            # tmux draws byte by byte: pass it through unchanged, and end the
            # connection once the server reports its cleanup complete.
            if [ -x /usr/bin/perl ]; then
                /usr/bin/perl -e '$|=1; my $m=shift; my $t=""; while (sysread(STDIN, my $b, 4096)) { syswrite(STDOUT, $b); $t .= $b; exit 3 if index($t, $m) >= 0; $t = substr($t, -512) }' "$marker" <&3
                [ $? = 3 ] && closed=1 && kill "$ssh_pid" 2>/dev/null
            else
                cat <&3
            fi
        fi
        exec 3<&-
        wait "$ssh_pid"; code=$?
        ssh_pid=''
        if [ $missing = 1 ] && [ $opened = 0 ] && [ $code = 10 ] && [ $attempt = 1 ]; then
            copy_scripts
            continue
        fi
        break
    done
    [ $closed = 1 ] && code=0
    say "SSH ended (exit code $code)."
    if [ $code != 0 ] && [ $opened = 0 ] && [ -n "$LOGIN_KEY" ]; then refused_note; fi
    if [ $code != 0 ]; then printf 'Press Enter to close. '; read -r _ < /dev/tty; fi
}

# -- DELFIN on this Mac -------------------------------------------------------

local_env() { printf '%s/local/env' "$STORE"; }

local_worker() {
    local voila port workdir
    printf '\033]0;DELFIN on this Mac\007'
    voila="$(local_env)/bin/delfin-voila"
    if [ ! -x "$voila" ]; then
        say 'DELFIN is not set up on this Mac yet. Choose "Set up DELFIN on this Mac" in DELFIN.'
        printf 'Press Enter to close. '; read -r _; exit 1
    fi
    workdir=${DELFIN_LOCAL_WORKDIR:-"$HOME/DELFIN"}
    mkdir -p "$workdir" && cd "$workdir" || die "Cannot use $workdir as the working folder."
    port=$(free_local_port 8866)
    # Jupyter's kernels talk over Unix sockets only when the runtime folder's
    # path is short enough (104 bytes on macOS); ~/Library/Jupyter/runtime
    # is not, and kernels then fall back to loopback TCP, which other
    # accounts of this Mac can reach. A short folder of this account's own,
    # taken only when it is really ours: no link, our owner, 0700.
    runtime="/tmp/delfin-jupyter-$(id -u)"
    mkdir -m 700 "$runtime" 2>/dev/null
    if [ -d "$runtime" ] && [ ! -L "$runtime" ] && [ -O "$runtime" ] \
        && [ "$(stat -f %Lp "$runtime" 2>/dev/null || stat -c %a "$runtime")" = 700 ]; then
        export JUPYTER_RUNTIME_DIR="$runtime"
    fi
    say "DELFIN on this Mac: dashboard on 127.0.0.1, working folder $workdir."
    say 'The browser opens with the dashboard and signs in with its token.'
    # Bound to 127.0.0.1 only; the agent's commands run in macOS's Seatbelt
    # sandbox (DELFIN_AGENT_SANDBOX=auto), without network, as on a server.
    # Started from a terminal, delfin-voila runs its server in the
    # background, one per account: starting again opens the same one, and
    # it ends itself 10 minutes after its last browser window closed
    # (DELFIN's own rule, as on a server). "Stop DELFIN on this Mac" in the
    # menu ends it at once.
    env DELFIN_AGENT_SANDBOX=auto DELFIN_AGENT_SANDBOX_NETWORK=0 \
        "$voila" --open-browser --ip 127.0.0.1 --port "$port"
    code=$?
    say ''
    if [ $code = 0 ]; then say 'You can close this window; the dashboard keeps running in the background.'
    else say "The dashboard did not start (exit $code); the lines above say why."; fi
    printf 'Press Enter to close. '; read -r _
    return $code
}

local_stop() {
    local voila out
    voila="$(local_env)/bin/delfin-voila"
    [ -x "$voila" ] || return 0
    out=$("$voila" stop 2>&1 | tail -n 3)
    tell_user "${out:-DELFIN on this Mac is stopped.}"
}

# -- the app: menus with the dialogs macOS has (osascript) --------------------

# ask_list PROMPT ITEM... -> the chosen item, or nothing on Cancel.
ask_list() {
    osascript - "$@" <<'EOF' 2>/dev/null
on run argv
    set promptText to item 1 of argv
    set choices to items 2 thru -1 of argv
    set picked to choose from list choices with title "DELFIN" with prompt promptText default items {item 1 of choices} OK button name "OK" cancel button name "Cancel"
    if picked is false then return ""
    return item 1 of picked
end run
EOF
}

# ask_text PROMPT DEFAULT -> the text, or nothing on Cancel.
ask_text() {
    osascript - "$1" "$2" <<'EOF' 2>/dev/null
on run argv
    try
        set answer to display dialog (item 1 of argv) default answer (item 2 of argv) with title "DELFIN" buttons {"Cancel", "OK"} default button "OK"
        return text returned of answer
    on error
        return ""
    end try
end run
EOF
}

tell_user() {
    osascript - "$1" <<'EOF' >/dev/null 2>&1
on run argv
    display dialog (item 1 of argv) with title "DELFIN" buttons {"OK"} default button "OK"
end run
EOF
}

# ask_flags -> up to three lines: keep, working, remember (chosen ones).
ask_flags() {
    osascript - "$@" <<'EOF' 2>/dev/null
on run argv
    set choices to {"Keep session on the server after disconnect", "Working shell (Ctrl+B, then N)", "Remember login (recommended)"}
    set picked to choose from list choices with title "DELFIN" with prompt "Options for this connection (Cmd-click to select several):" default items argv with multiple selections allowed and empty selection allowed
    if picked is false then return "CANCEL"
    set out to ""
    repeat with p in picked
        set out to out & p & linefeed
    end repeat
    return out
end run
EOF
}

# Opens a Terminal window running `DELFIN.sh SUBCOMMAND [ID]`. The .command
# file is written here (no download mark) and run by Terminal.
in_terminal() {
    local file="$STORE/start-$1.command"
    mkdir -p "$STORE"
    printf '#!/bin/bash\nrm -f %s\nexec /bin/bash %s %s %s\n' "$(shq "$file")" "$(shq "$HERE/DELFIN.sh")" "$1" "${2:+$(shq "$2")}" > "$file"
    chmod 700 "$file"
    open -a Terminal "$file"
}

edit_profile() {
    local id=$1 file name target directory='' executable='' port=0 keep=0 working=0 remember=1 defaults flags
    if [ -n "$id" ]; then
        load_profile "$id"
        name=$P_NAME target=$P_TARGET directory=$P_DIRECTORY executable=$P_EXECUTABLE port=$P_PORT keep=$P_KEEP working=$P_WORKING remember=$P_REMEMBER
    else
        id=$(new_uuid)
    fi
    target=$(ask_text 'SSH server: user@host, or an alias from ~/.ssh/config' "$target")
    [ -n "$target" ] || return 1
    if ! valid_target "$target"; then tell_user 'Enter an SSH alias or user@hostname (no options).'; return 1; fi
    name=$(ask_text 'Name of this connection' "${name:-$target}")
    [ -n "$name" ] || name=$target
    defaults=()
    [ "$keep" = 1 ] && defaults=("${defaults[@]}" 'Keep session on the server after disconnect')
    [ "$working" = 1 ] && defaults=("${defaults[@]}" 'Working shell (Ctrl+B, then N)')
    [ "$remember" = 1 ] && defaults=("${defaults[@]}" 'Remember login (recommended)')
    flags=$(ask_flags "${defaults[@]}")
    [ "$flags" = CANCEL ] && return 1
    case $flags in *Keep*) keep=1 ;; *) keep=0 ;; esac
    case $flags in *Working*) working=1 ;; *) working=0 ;; esac
    case $flags in *Remember*) remember=1 ;; *) remember=0 ;; esac
    if [ "$(ask_list 'Advanced settings (server paths, port)?' 'Leave automatic' 'Set them now')" = 'Set them now' ]; then
        directory=$(ask_text 'Working directory on the server (empty: automatic)' "$directory")
        executable=$(ask_text 'DELFIN location on the server (empty: automatic)' "$executable")
        port=$(ask_text 'Starting port (0: automatic)' "$port"); port=${port:-0}
    fi
    ( save_profile "$id" "$name" "$target" "$directory" "$executable" "$port" "$keep" "$working" "$remember" ) 2>"$STORE/.error" \
        || { tell_user "$(cat "$STORE/.error")"; rm -f "$STORE/.error"; return 1; }
    rm -f "$STORE/.error"
    printf '%s' "$id"
}

app_menu() {
    local ids=() items=() id choice index local_item
    mkdir -p "$STORE"
    while IFS= read -r id; do ids=("${ids[@]}" "$id"); done < <(profile_ids)
    index=1
    for id in "${ids[@]}"; do
        items=("${items[@]}" "$index. Connect: $(profile_get "$(profile_file "$id")" name)")
        index=$((index + 1))
    done
    if [ -x "$(local_env)/bin/delfin-voila" ]; then local_item='DELFIN on this Mac'; else local_item='Set up DELFIN on this Mac...'; fi
    items=("${items[@]}" "$local_item")
    [ "$local_item" = 'DELFIN on this Mac' ] && items=("${items[@]}" 'Stop DELFIN on this Mac')
    items=("${items[@]}" 'New server connection...')
    [ ${#ids[@]} -gt 0 ] && items=("${items[@]}" 'Edit a connection...' 'Remove a connection...')
    items=("${items[@]}" 'Uninstall DELFIN...')
    choice=$(ask_list 'Where do you want to work?' "${items[@]}")
    case $choice in
        '') return 0 ;;
        [0-9]*'. Connect: '*)
            index=${choice%%.*}
            id=${ids[$((index - 1))]}
            touch "$(profile_file "$id")"
            in_terminal connect "$id" ;;
        'DELFIN on this Mac') in_terminal local ;;
        'Stop DELFIN on this Mac') local_stop ;;
        'Set up DELFIN on this Mac...') in_terminal setup-local ;;
        'New server connection...') id=$(edit_profile '') && in_terminal connect "$id" ;;
        'Edit a connection...'|'Remove a connection...')
            items=(); index=1
            for id in "${ids[@]}"; do items=("${items[@]}" "$index. $(profile_get "$(profile_file "$id")" name)"); index=$((index + 1)); done
            local picked
            picked=$(ask_list "${choice%...}:" "${items[@]}")
            [ -n "$picked" ] || return 0
            id=${ids[$((${picked%%.*} - 1))]}
            if [ "$choice" = 'Edit a connection...' ]; then edit_profile "$id" >/dev/null
            else rm -f "$(profile_file "$id")" "$STORE/invite-$id.key"; fi ;;
        'Uninstall DELFIN...') in_terminal uninstall ;;
    esac
}

# -- package check, install, uninstall ----------------------------------------

# SHA256SUMS.txt lies at the package root, beside DELFIN.app and
# Install.command, and names every file of the package by its path there.
validate_package() {
    local dir=${1:-$PKG} name want got
    [ -f "$dir/SHA256SUMS.txt" ] || die 'SHA256SUMS.txt is missing: this is not a complete launcher package.'
    while read -r want name; do
        [ -n "$name" ] || continue
        case $name in /*|*..*) die "Invalid path in SHA256SUMS.txt: $name" ;; esac
        [ -f "$dir/$name" ] || die "Missing from the package: $name"
        got=$(sha256_of "$dir/$name")
        [ "$got" = "$want" ] || die "Changed since the package was built: $name"
    done < "$dir/SHA256SUMS.txt"
    for name in Install.command Uninstall.command DELFIN.app/Contents/MacOS/DELFIN         DELFIN.app/Contents/Resources/DELFIN.sh DELFIN.app/Contents/Resources/setup_local.sh         DELFIN.app/Contents/Resources/fetch_micromamba.sh DELFIN.app/Contents/Resources/remote_bootstrap.sh; do
        /bin/bash -n "$dir/$name" || die "Syntax error in $name"
    done
}

# The app's own processes: connections and local dashboards started from it.
running_app() { pgrep -f "$APPS/DELFIN.app/Contents/Resources/DELFIN.sh (connect|local)" 2>/dev/null; }

close_running() {
    local pids answer
    pids=$(running_app) || return 0
    [ -n "$pids" ] || return 0
    say 'DELFIN is open. It has to be closed first.'
    say 'Open connections end; a session kept on the server stays there.'
    if [ "$1" = yes ]; then answer=y; else printf 'Close DELFIN now and continue? (y/N) '; read -r answer; fi
    case $answer in y|Y) ;; *) die 'DELFIN is still open. Nothing was changed.' ;; esac
    for pid in $pids; do pkill -TERM -P "$pid" 2>/dev/null; kill "$pid" 2>/dev/null; done
    sleep 1
    [ -z "$(running_app)" ] || die 'DELFIN could not be closed. Close its Terminal windows and try again.'
}

# The app is built and checked with the package (and, for a signed
# release, signed and notarized there): it is copied as it is, never
# assembled or changed here, so a signature on it stays valid. ditto keeps
# the bundle's attributes and signature; a download mark the person chose to
# keep stays on the copy too, and macOS asks there in turn.
write_app() {
    local app="$APPS/DELFIN.app" new="$APPS/.DELFIN.app.new"
    [ -f "$BUNDLE/Contents/Info.plist" ] || die 'DELFIN.app is missing from the package.'
    mkdir -p "$APPS"
    rm -rf "$new"
    if command -v ditto >/dev/null 2>&1; then ditto "$BUNDLE" "$new"; else cp -Rp "$BUNDLE" "$new"; fi
    rm -rf "$app"
    mv "$new" "$app"
}

# An invitation (delfin-people --invite): invite.json and invite_key beside
# Install.command become the first saved connection, the server's host keys
# go into known_hosts (no fingerprint question), and the key is kept for the
# first Start, which registers this Mac's own key with it and deletes it.
import_invitation() {
    local json="$PKG/invite.json" key="$PKG/invite_key" target keys id existing server hostkey
    [ -f "$json" ] && [ -f "$key" ] || return 0
    target=$(grep -oE '"target"[[:space:]]*:[[:space:]]*"[^"]*"' "$json" | head -n 1 | sed -E 's/.*:[[:space:]]*"([^"]*)"$/\1/')
    [[ $target =~ ^[a-z_][a-z0-9_-]{0,31}@[A-Za-z0-9][A-Za-z0-9.-]{0,200}$ ]] || die 'The invitation names no valid server.'
    keys=$(grep -oE '"(ssh-ed25519|ecdsa-sha2-nistp(256|384|521)|ssh-rsa) [A-Za-z0-9+/]+={0,2}"' "$json" | tr -d '"')
    [ -n "$keys" ] || die 'The invitation carries no valid server key.'
    [ "$(head -n 1 "$key")" = '-----BEGIN OPENSSH PRIVATE'' KEY-----' ] || die 'The invitation key is not readable.'
    id=''
    for existing in $(profile_ids); do
        if [ "$(profile_get "$(profile_file "$existing")" target)" = "$target" ]; then id=$existing; break; fi
    done
    if [ -n "$id" ]; then
        load_profile "$id"
        save_profile "$id" "$P_NAME" "$target" "$P_DIRECTORY" "$P_EXECUTABLE" "$P_PORT" "$P_KEEP" "$P_WORKING" 1
    else
        id=$(new_uuid)
        save_profile "$id" "$target" "$target" '' '' 0 1 0 1
    fi
    mkdir -p "$SSH_DIR" && chmod 700 "$SSH_DIR"
    server=${target#*@}
    touch "$SSH_DIR/known_hosts"
    while IFS= read -r hostkey; do
        [ -n "$hostkey" ] || continue
        grep -qxF "$server $hostkey" "$SSH_DIR/known_hosts" || printf '%s %s\n' "$server" "$hostkey" >> "$SSH_DIR/known_hosts"
    done <<EOF
$keys
EOF
    # ssh ignores a private key others may read: owner only.
    rm -f "$STORE/invite-$id.key"
    cp "$key" "$STORE/invite-$id.key"
    chmod 600 "$STORE/invite-$id.key"
    INVITED=$target
}

install_app() {
    local yes=no
    [ "$1" = --close-running ] && yes=yes
    validate_package "$PKG"
    close_running $yes
    mkdir -p "$STORE" && chmod 700 "$STORE"
    rm -f "$STORE"/start-*.command
    write_app
    INVITED=''
    import_invitation
    say ''
    say 'DELFIN is installed.'
    [ -n "$INVITED" ] && say "  Your connection to $INVITED is set up. Start DELFIN and choose it: no password is needed."
    say "  Start it:  DELFIN in $APPS (Finder: Go > Home > Applications), or Spotlight (Cmd+Space, DELFIN)."
    say '  Dock:      while it runs, right-click its icon > Options > Keep in Dock.'
    say '  Remove:    Uninstall.command, or "Uninstall DELFIN..." in the app.'
    say '  Update:    run Install.command of a newer version. Saved connections are kept.'
}

# uninstall [--keep-profiles|--remove-profiles] [--keep-local|--remove-local] [--close-running]
uninstall_app() {
    local keep='' local_choice='' yes=no arg answer item
    for arg in "$@"; do
        case $arg in
            --keep-profiles) keep=1 ;; --remove-profiles) keep=0 ;;
            --keep-local) local_choice=keep ;; --remove-local) local_choice=remove ;;
            --close-running) yes=yes ;;
        esac
    done
    close_running $yes
    if [ -z "$keep" ]; then
        printf 'Keep your saved connections for a later install? (Y/n) '; read -r answer
        case $answer in n|N) keep=0 ;; *) keep=1 ;; esac
    fi
    if [ -d "$STORE/local" ] && [ -z "$local_choice" ]; then
        printf 'Also remove DELFIN installed on this Mac (%s)? Your work in ~/DELFIN stays. (y/N) ' "$STORE/local"; read -r answer
        case $answer in y|Y) local_choice=remove ;; *) local_choice=keep ;; esac
    fi
    [ -L "$APPS/DELFIN.app" ] && die 'The installed app is a link. Refusing recursive removal.'
    rm -rf "$APPS/DELFIN.app" "$APPS/.DELFIN.app.new"
    if [ -d "$STORE" ]; then
        [ -L "$STORE" ] && die 'The settings folder is a link. Refusing recursive removal.'
        for item in "$STORE"/* "$STORE"/.[!.]*; do
            [ -e "$item" ] || [ -L "$item" ] || continue
            case ${item##*/} in
                profiles|invite-*.key) [ "$keep" = 1 ] && continue ;;
                local) [ "$local_choice" = keep ] && continue ;;
            esac
            rm -rf "$item"
        done
        rmdir "$STORE" 2>/dev/null
    fi
    say ''
    say 'DELFIN is removed from this Mac.'
    say '  Nothing on any server was changed: your sessions and calculations there are untouched.'
    [ "$keep" = 1 ] && say "  Saved connections are kept in $STORE for a later install."
    [ "$local_choice" = keep ] && say "  DELFIN on this Mac is kept in $STORE/local."
    say "  The login key $SSH_DIR/delfin_launcher_ed25519 stays; delete it yourself if no server should know this Mac any more."
    say '  If you kept DELFIN in the Dock: right-click its icon > Options > Remove from Dock.'
}

if [ "${DELFIN_LAUNCHER_LIBRARY:-}" = 1 ]; then return 0 2>/dev/null; fi

case ${1:-app} in
    app) app_menu ;;
    connect) connect_worker "$2" ;;
    local) local_worker ;;
    setup-local)
        if /bin/bash "$HERE/setup_local.sh"; then
            printf 'Start DELFIN on this Mac now? (Y/n) '; read -r answer
            case $answer in n|N) ;; *) local_worker ;; esac
        else
            printf 'Press Enter to close. '; read -r _
        fi ;;
    install) shift; install_app "$@" ;;
    uninstall) shift; uninstall_app "$@" ;;
    validate) validate_package "${2:-$PKG}"; say 'Package OK.' ;;
    *) die "Unknown command: $1" ;;
esac
