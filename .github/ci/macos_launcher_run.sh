#!/bin/bash
# The macOS launcher for real, on macos-latest (.github/workflows/macos-launcher.yml):
#
#  1. DELFIN set up on the Mac by setup_local.sh (from this tree), and its
#     dashboard started by `DELFIN.sh local`: answering with the token,
#     refusing without it, listening on 127.0.0.1 only; the agent's shell
#     isolated by Seatbelt.
#  2. An invitation, imported by the installer, used by `DELFIN.sh connect`
#     against an sshd on this runner (127.0.0.1:2222, keys only): the
#     invitation key enrolls this Mac's own key and is deleted, the server
#     scripts are copied once, the dashboard comes up in tmux behind the
#     tunnel, and the browser URL (written to a file instead of opened)
#     reaches it through this connection's own ssh.
#
# Not covered: Finder, Gatekeeper's dialog, Terminal.app, a browser, and a
# password or OTP typed by a person.
set -eo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd -P)
work=$(mktemp -d /tmp/delfin-mac-ci.XXXXXX)
fail() { echo "FAIL: $*" >&2; exit 1; }
ditto -x -k "$repo/tools/delfin-macos-launcher.zip" "$work/x"
pkg="$work/x/macos-launcher"
export DELFIN_LAUNCHER_STORE="$work/store" DELFIN_LAUNCHER_APPS="$work/Applications"
cleanup() {
    pkill -f 'DELFIN.sh connect' 2>/dev/null || true
    [ -n "${local_pid:-}" ] && kill "$local_pid" 2>/dev/null || true
    [ -f "$work/sshd/pid" ] && kill "$(cat "$work/sshd/pid")" 2>/dev/null || true
    command -v tmux >/dev/null && tmux kill-server 2>/dev/null || true
}
trap cleanup EXIT

wait_for() { # SECONDS COMMAND...
    local deadline=$((SECONDS + $1)); shift
    until "$@"; do [ $SECONDS -lt $deadline ] || return 1; sleep 2; done
}

echo '== 1. DELFIN on this Mac'
DELFIN_LOCAL_SOURCE="$repo" /bin/bash "$pkg/setup_local.sh"
env="$work/store/local/env"
"$env/bin/python" -c 'from delfin.agent.isolation_posture import isolation_status as s; r = s(); print("isolation:", r); assert r.isolated and "eatbelt" in r.layer, r'
mkdir -p "$work/DELFIN"
DELFIN_LOCAL_WORKDIR="$work/DELFIN" BROWSER=true /bin/bash "$pkg/DELFIN.sh" local > "$work/local.log" 2>&1 < /dev/null &
local_pid=$!
wait_for 120 grep -q 'dashboard on 127.0.0.1:' "$work/local.log" || { cat "$work/local.log"; fail 'no port'; }
port=$(sed -n 's/.*dashboard on 127\.0\.0\.1:\([0-9]*\),.*/\1/p' "$work/local.log" | head -n 1)
record() { "$env/bin/python" -I -c 'import json, sys; from delfin.agent import where; r = where.dashboard(); print(r.get("token", "") if r.get("running") and int(r.get("port", 0)) == int(sys.argv[1]) else "")' "$port"; }
wait_for 180 test -n "$(record)" || { cat "$work/local.log"; fail 'the dashboard did not come up'; }
dashboard_key=$(record)
code=$(curl -s -o /dev/null -w '%{http_code}' -L --max-redirs 5 "http://127.0.0.1:$port/?token=$dashboard_key")
[ "$code" = 200 ] || fail "with the token: HTTP $code"
code=$(curl -s -o /dev/null -w '%{http_code}' "http://127.0.0.1:$port/voila/render/delfin_dashboard.ipynb")
case $code in 302|401|403) ;; *) fail "without the token: HTTP $code" ;; esac
lsof -nP -iTCP:"$port" -sTCP:LISTEN
if lsof -nP -iTCP:"$port" -sTCP:LISTEN | awk 'NR > 1 {print $9}' | grep -v '^127\.0\.0\.1:'; then fail 'listens beyond 127.0.0.1'; fi
for address in $(ifconfig | awk '/inet / && $2 != "127.0.0.1" {print $2}'); do
    if curl -s -o /dev/null --max-time 3 "http://$address:$port/"; then fail "reachable on $address"; fi
done
echo "local dashboard: token 200, no token $code, 127.0.0.1 only"
kill "$local_pid"; wait "$local_pid" 2>/dev/null || true; local_pid=''
wait_for 30 sh -c "! lsof -nP -iTCP:$port -sTCP:LISTEN >/dev/null" || fail 'the dashboard outlived its window'

echo '== 2. An invitation and a connection over SSH'
command -v tmux >/dev/null || brew install tmux
brew_bin=$(dirname "$(command -v tmux)")
# The check on the server is GNU sha256sum, as on a Linux server.
if command -v sha256sum >/dev/null; then sha_dir=$(dirname "$(command -v sha256sum)")
else brew install coreutils; sha_dir="$(brew --prefix coreutils)/libexec/gnubin"; fi
mkdir -p "$work/sshd"
ssh-keygen -q -t ed25519 -N '' -f "$work/sshd/host_key"
cat > "$work/sshd/config" <<EOF
Port 2222
ListenAddress 127.0.0.1
HostKey $work/sshd/host_key
PidFile $work/sshd/pid
PasswordAuthentication no
KbdInteractiveAuthentication no
UsePAM no
StrictModes no
SetEnv PATH=$brew_bin:$sha_dir:/usr/bin:/bin:/usr/sbin:/sbin
EOF
/usr/sbin/sshd -f "$work/sshd/config" -E "$work/sshd/log" || { cat "$work/sshd/log"; fail 'sshd'; }
user=$(id -un)
mkdir -p ~/.ssh && chmod 700 ~/.ssh
# The server's port lives in ~/.ssh/config, as for any server on another port.
printf '\nHost 127.0.0.1\n    Port 2222\n    HostKeyAlias 127.0.0.1\n' >> ~/.ssh/config
# What `delfin-people --invite ci` leaves (that script needs GNU date and is
# run on the Linux server; tests/test_the_macos_launcher_package_is_current.py
# runs it there): a key that can only enroll, and invite.json.
ssh-keygen -q -t ed25519 -N '' -C delfin-invite-ci -f "$work/invite_key"
expiry=$(date -u -v+1d +%Y%m%d%H%MZ)
printf 'restrict,command="/bin/bash %s --enroll ci",expiry-time="%s" %s delfin-invite-ci\n' \
    "$repo/scripts/pi/delfin_people.sh" "$expiry" "$(cut -d' ' -f1-2 "$work/invite_key.pub")" >> ~/.ssh/authorized_keys
chmod 600 ~/.ssh/authorized_keys
host_key=$(cut -d' ' -f1-2 "$work/sshd/host_key.pub")
printf '{\n  "target": "%s@127.0.0.1",\n  "name": "ci",\n  "expires": "%s",\n  "host_keys": [\n    "%s"\n  ]\n}\n' \
    "$user" "$expiry" "$host_key" > "$pkg/invite.json"
cp "$work/invite_key" "$pkg/invite_key"
/bin/bash "$pkg/DELFIN.sh" install < /dev/null
grep -qx "127.0.0.1 $host_key" ~/.ssh/known_hosts || fail 'the host key is not in known_hosts'
profile=$(ls "$work/store/profiles"/*.conf)
id=$(basename "$profile" .conf)
[ -f "$work/store/invite-$id.key" ] || fail 'the invitation key was not kept'
# The server's DELFIN is the one set up above (on a real server: its install).
sed -i '' "s#^executable=.*#executable=$env/bin/python#" "$profile"
app="$work/Applications/DELFIN.app/Contents/Resources"

connect() {
    rm -f "$work/url"
    # script(1) gives the connection the terminal Terminal.app would; the
    # sleep keeps its input open, as a window does.
    (sleep 900 | DELFIN_LAUNCHER_URL_FILE="$work/url" script -q "$work/connect-$1.log" /bin/bash "$app/DELFIN.sh" connect "$id") &
    wait_for 300 test -s "$work/url" || { cat "$work/connect-$1.log"; tail -n 30 "$work/sshd/log"; fail "no browser URL ($1)"; }
    url=$(cat "$work/url")
    echo "browser URL: ${url%%token=*}token=..."
    tunnel_port=$(printf '%s' "$url" | sed -E 's#^http://127\.0\.0\.1:([0-9]+)/.*#\1#')
    lsof -nP -iTCP:"$tunnel_port" -sTCP:LISTEN
    [ "$(lsof -nP -iTCP:"$tunnel_port" -sTCP:LISTEN | awk 'NR > 1 {print $1}' | sort -u)" = ssh ] || fail 'the tunnel is not ssh'
    code=$(curl -s -o /dev/null -w '%{http_code}' -L --max-redirs 5 "$url")
    [ "$code" = 200 ] || fail "through the tunnel with the token: HTTP $code"
    code=$(curl -s -o /dev/null -w '%{http_code}' "http://127.0.0.1:$tunnel_port/voila/render/delfin_dashboard.ipynb")
    case $code in 302|401|403) ;; *) fail "through the tunnel without the token: HTTP $code" ;; esac
}

connect first
[ ! -f "$work/store/invite-$id.key" ] || fail 'the invitation key was not used up'
grep -q 'delfin-launcher-' ~/.ssh/authorized_keys || fail "this Mac's key was not enrolled"
if grep -q 'delfin-invite-ci$' ~/.ssh/authorized_keys; then fail 'the invitation still works'; fi
ls ~/.cache/delfin-launcher/*/remote_bootstrap.sh || fail 'the scripts were not copied'
grep -q 'Copying the DELFIN launcher scripts' "$work/connect-first.log" || fail 'no copy on the first connection'

echo '== disconnect (window closed) and connect again: Keep session is on for an invitation'
pkill -HUP -f "$app/DELFIN.sh connect" || true
wait_for 30 sh -c "! pgrep -f '$app/DELFIN.sh connect' >/dev/null" || fail 'the connection did not end'
tmux ls || fail 'Keep session: the dashboard should stay on the server'
connect second
if grep -q 'Copying the DELFIN launcher scripts' "$work/connect-second.log"; then fail 'copied again'; fi
grep -q 'Registering this Mac' "$work/connect-second.log" && fail 'enrolled again'
echo 'OK: local dashboard and SSH connection on macOS'
