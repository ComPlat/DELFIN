#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# delfin-deploy: the narrow root step that deploys the SHARED installation
# (one Linux account per person, pi_bootstrap.sh with DELFIN_PI_ACCOUNTS=1).
# Installed root-owned as /usr/local/sbin/delfin-deploy by the bootstrap.
#
#   sudo delfin-deploy install VERSION TARBALL [--test]
#   sudo delfin-deploy rollback
#   sudo delfin-deploy status
#
# deploy_pi.sh --shared pi@<pi> runs it: the admin uploads the source
# tarball and types the sudo password. The agent account never deploys the
# shared installation, and nothing it can write is run or loaded by a
# person's account:
#   - root owns $DELFIN_SHARED_ROOT, swaps src/ in it, and writes agent.env,
#     settings-defaults.json (from /etc/delfin) and the launcher package that
#     invitations are built from (/usr/local/lib/delfin);
#   - the build itself (tar, uv, pip, the extras) runs as the install account
#     ($DELFIN_INSTALL_USER: no login, no group the agent or a person is in),
#     never as root;
#   - the agent account's own part (its settings, wrappers, services, lab
#     units) runs as the agent account and reads the installation only.
# A deploy ends with a check that the agent account owns nothing in it.
# ============================================================================

CONF="${DELFIN_ACCOUNTS_CONF:-/etc/delfin/accounts.conf}"
UV=/usr/local/lib/delfin/bin/uv
LIB=/usr/local/lib/delfin

die() { echo "delfin-deploy: $*" >&2; exit 1; }

[[ $EUID -eq 0 ]] || die "run with sudo"
[[ -f "$CONF" ]] || die "no $CONF: not a box with one account per person"
conf() { sed -n "s/^$1=//p" "$CONF" | head -n 1; }
ROOT="$(conf DELFIN_SHARED_ROOT)"
AGENT="$(conf DELFIN_DEPLOY_USER)"
INSTALLER="$(conf DELFIN_INSTALL_USER)"
[[ "$ROOT" =~ ^/[A-Za-z0-9._/-]+$ && "$ROOT" != / ]] || die "$CONF: no valid DELFIN_SHARED_ROOT"
[[ "$AGENT" =~ ^[a-z_][a-z0-9_-]{0,31}$ ]] || die "$CONF: no valid DELFIN_DEPLOY_USER"
[[ "$INSTALLER" =~ ^[a-z_][a-z0-9_-]{0,31}$ ]] || die "$CONF: no DELFIN_INSTALL_USER (run pi_bootstrap.sh with DELFIN_PI_ACCOUNTS=1 again)"
[[ "$INSTALLER" != "$AGENT" ]] || die "the install account must not be the agent account"
INSTALL_HOME="$(getent passwd "$INSTALLER" | cut -d: -f6)"
AGENT_HOME="$(getent passwd "$AGENT" | cut -d: -f6)"
AGENT_UID="$(id -u "$AGENT")"
[[ -x "$UV" ]] || die "$UV missing (pi_bootstrap.sh with DELFIN_PI_ACCOUNTS=1)"

# DELFIN_PI_* switches the admin passed (deploy_pi.sh forwards them), checked
# like deploy_pi.sh checks them: plain names, plain values.
pi_vars=()
while IFS='=' read -r name value; do
    [[ "$name" =~ ^DELFIN_PI_[A-Z0-9_]+$ && "$value" =~ ^[A-Za-z0-9\ ._,:=/+-]*$ ]] && pi_vars+=("$name=$value")
done < <(env)

# A child with nothing of root's environment.
as_installer() {
    runuser -u "$INSTALLER" -- env -i HOME="$INSTALL_HOME" USER="$INSTALLER" LOGNAME="$INSTALLER" \
        PATH=/usr/local/bin:/usr/bin:/bin LANG=C.UTF-8 \
        DELFIN_SHARED_ROOT="$ROOT" DELFIN_UV="$UV" UV_CACHE_DIR="$INSTALL_HOME/.cache/uv" \
        UV_PYTHON_INSTALL_DIR="$ROOT/python" "${pi_vars[@]}" bash -c 'cd "$HOME" && exec "$@"' _ "$@"
}
as_agent() {
    runuser -u "$AGENT" -- env -i HOME="$AGENT_HOME" USER="$AGENT" LOGNAME="$AGENT" \
        PATH="$AGENT_HOME/.local/bin:/usr/local/bin:/usr/bin:/bin" LANG=C.UTF-8 \
        XDG_RUNTIME_DIR="/run/user/$AGENT_UID" \
        DBUS_SESSION_BUS_ADDRESS="unix:path=/run/user/$AGENT_UID/bus" "${pi_vars[@]}" bash -c 'cd "$HOME" && exec "$@"' _ "$@"
}

# agent.env: root's file. The upgrades pi_install.sh applies on a one-account
# box (agent_env_upgrades), here as root.
root_agent_env() {
    local f="$ROOT/agent.env"
    [[ -f "$f" && ! -L "$f" ]] || die "$f missing or a link: run pi_bootstrap.sh with DELFIN_PI_ACCOUNTS=1 again"
    sed -i 's/^export DELFIN_AGENT_SANDBOX=off$/export DELFIN_AGENT_SANDBOX=auto/' "$f"
    grep -q '^export DELFIN_MCP_ISOLATION=' "$f" || printf '%s\n' \
        "# DELFIN's own MCP servers run in bwrap, seeing only their folders." \
        "export DELFIN_MCP_ISOLATION=builtin" >> "$f"
    grep -q '^export PLAYWRIGHT_BROWSERS_PATH=' "$f" || printf '%s\n' \
        "# One Chromium for every account (shared installation)." \
        "export PLAYWRIGHT_BROWSERS_PATH=$ROOT/ms-playwright" >> "$f"
    chown root:root "$f"
    chmod 644 "$f"
}

# Every person's defaults: from root's /etc/delfin/pi-agent.conf and
# pi-tabs.conf, never from the agent account's settings. The tab ids are read
# as text from the deployed dashboard (nothing of the tree is run as root).
root_settings_defaults() {
    python3 -I - /etc/delfin/pi-agent.conf /etc/delfin/pi-tabs.conf \
        "$ROOT/src/delfin/dashboard/__init__.py" "$ROOT/settings-defaults.json" <<'EOF'
import json, os, re, sys
agent_conf, tabs_conf, dashboard, out_path = sys.argv[1:5]
keys = ("mode", "provider", "model", "permission_profile", "effort")
agent = {}
if os.path.isfile(agent_conf):
    for line in open(agent_conf, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = (part.strip() for part in line.partition("="))
        if key in keys:
            agent[key] = value or None
        else:
            print(f"warning: {agent_conf}: unknown key {key!r}", file=sys.stderr)
out = {"agent": agent}
if os.path.isfile(tabs_conf) and os.path.isfile(dashboard):
    shown = [w for line in open(tabs_conf, encoding="utf-8")
             if not line.lstrip().startswith("#") for w in line.split()]
    ids = list(dict.fromkeys(re.findall(r"'id':\s*'([a-z_]+)'",
                                        open(dashboard, encoding="utf-8").read())))
    if ids:
        order = [t for t in shown if t in ids] + [t for t in ids if t not in shown and t != "settings"]
        hidden = [t for t in ids if t not in shown and t != "settings"]
        out["ui"] = {"tabs": {"order": order, "hidden": hidden}}
tmp = out_path + ".tmp"
fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC | os.O_NOFOLLOW, 0o644)
with os.fdopen(fd, "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2)
    fh.write("\n")
os.chmod(tmp, 0o644)
os.replace(tmp, out_path)
print("==> every person's defaults (root's): " + out_path)
EOF
    chown root:root "$ROOT/settings-defaults.json"
}

# What invitations are built from, and the root-owned tools: copied out of
# the tree the admin just deployed, into root's directory.
root_tools() {
    install -d -m 755 -o root -g root "$LIB"
    if [[ -f "$ROOT/src/tools/delfin-windows-launcher.zip" ]]; then
        install -m 644 -o root -g root "$ROOT/src/tools/delfin-windows-launcher.zip" "$LIB/delfin-windows-launcher.zip"
    fi
    install -m 755 -o root -g root "$ROOT/src/scripts/pi/delfin_people.sh" "$LIB/delfin-people"
    install -m 755 -o root -g root "$ROOT/src/scripts/pi/delfin_accounts.sh" /usr/local/sbin/delfin-accounts
}

# The boundary, measured after every deploy: the root directory and the two
# setting files are root's, and the agent account owns nothing in the tree.
check_boundary() {
    local bad
    [[ "$(stat -c '%U %a' "$ROOT")" == "root 750" ]] || die "$ROOT is not root's 750"
    for f in agent.env settings-defaults.json; do
        [[ "$(stat -c '%U %a' "$ROOT/$f")" == "root 644" ]] || die "$ROOT/$f is not root's 644"
    done
    bad="$(find "$ROOT" -xdev -path "$ROOT/work" -prune -o -user "$AGENT" -print -quit)"
    [[ -z "$bad" ]] || die "the agent account owns $bad inside $ROOT"
    echo "==> boundary: $ROOT root's; agent.env and settings-defaults.json root's; nothing in it is $AGENT's"
}

build_and_finish() {
    local run_tests="${1:-}"
    as_installer bash "$ROOT/src/scripts/pi/pi_install.sh" shared-build
    root_agent_env
    root_settings_defaults
    root_tools
    chown root:delfin-shared "$ROOT"
    chmod 750 "$ROOT"
    check_boundary
    as_agent bash "$ROOT/src/scripts/pi/pi_install.sh" account
    if [[ "$run_tests" == "--test" ]]; then
        echo "==> pytest (not slow), as $AGENT"
        as_agent bash -c "source '$ROOT/.venv/bin/activate' && cd '$ROOT/src' && python -m pytest -q -m 'not slow' -p no:cacheprovider"
    fi
}

cmd="${1:-}"
case "$cmd" in
install)
    version="${2:-}"; tarball="${3:-}"; run_tests="${4:-}"
    [[ "$version" =~ ^[A-Za-z0-9._+-]{1,80}$ ]] || die "VERSION: letters, digits, . _ + -"
    [[ -f "$tarball" && ! -L "$tarball" ]] || die "TARBALL: no such file"
    [[ -z "$run_tests" || "$run_tests" == "--test" ]] || die "unknown option: $run_tests"
    cd "$ROOT"
    # Unpacked by the install account, from a pipe root opened: tar never
    # runs as root, and nothing it unpacks can land outside src.new.
    rm -rf "$ROOT/src.new"
    install -d -m 755 -o "$INSTALLER" -g "$INSTALLER" "$ROOT/src.new"
    as_installer tar -xzf - -C "$ROOT/src.new" < "$tarball"
    as_installer bash -c 'printf "%s deployed %s\n" "$1" "$(date -Is)" > "$2/src.new/DEPLOYED_VERSION"' _ "$version" "$ROOT"
    if [[ -d "$ROOT/src" ]]; then
        as_installer bash "$ROOT/src.new/scripts/pi/pi_install.sh" carry-tui "$ROOT/src" "$ROOT/src.new"
    fi
    rm -rf "$ROOT/src.prev"
    [[ -d "$ROOT/src" ]] && mv "$ROOT/src" "$ROOT/src.prev"
    mv "$ROOT/src.new" "$ROOT/src"
    build_and_finish "$run_tests"
    echo "==> deployed $version (shared installation $ROOT)"
    ;;
rollback)
    cd "$ROOT"
    [[ -d "$ROOT/src.prev" ]] || die "no src.prev to roll back to"
    mv "$ROOT/src" "$ROOT/src.bad"
    mv "$ROOT/src.prev" "$ROOT/src"
    mv "$ROOT/src.bad" "$ROOT/src.prev"
    as_installer bash "$ROOT/src/scripts/pi/pi_install.sh" carry-tui "$ROOT/src.prev" "$ROOT/src"
    build_and_finish
    ;;
status)
    cat "$ROOT/src/DEPLOYED_VERSION" 2>/dev/null || echo "nothing deployed"
    ;;
*)
    sed -n '/^#   sudo delfin-deploy install/,/^#   sudo delfin-deploy status/p' "$0" | sed 's/^# \{0,1\}//' >&2
    exit 2
    ;;
esac
