#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# One-time setup of a Raspberry Pi 5 (aarch64) as a DELFIN agent test box.
#
# Runs ON the Pi as a sudo-capable login user (deploy_pi.sh --bootstrap does
# this for you). Idempotent: safe to run again.
#
# What this script does:
#   1. Install apt packages (git, build tools, X/GL libs for rdkit drawing)
#   2. Create a dedicated user (default: delfin) WITHOUT root: it may only
#      restart its dashboard service (AGENT_SUDO, below). The login user
#      (pi) administers the box.
#   3. Copy the invoking user's authorized_keys to that user
#   4. Install uv + Python 3.11 for that user and create ~/delfin/.venv
#   5. Write ~/delfin/agent.env and the ~/.local/bin/delfin-agent wrapper
#
# Code is NOT installed here - deploy_pi.sh pushes it to ~/delfin/src.
#
# Usage (on the Pi):
#   sudo bash pi_bootstrap.sh                   # agent may restart its service only
#   AGENT_SUDO=none sudo bash pi_bootstrap.sh   # agent user without sudo
#   AGENT_SUDO=full sudo bash pi_bootstrap.sh   # root without password (not advised)
#
# One Linux account per person (opt-in; the default is people mode, everyone
# in the agent user's account under ~/people/<name>):
#   DELFIN_PI_ACCOUNTS=1 sudo bash pi_bootstrap.sh
# puts the installation into $DELFIN_PI_SHARED_ROOT (default /opt/delfin,
# owned by the agent user, readable by the group delfin-people only, 0750;
# ~/delfin is a link to it),
# Python included, and installs delfin-accounts (delfin_accounts.sh, beside
# this script) for the admin: sudo delfin-accounts add NAME. On a box that
# already runs the one-account layout (~/delfin a directory), add
# DELFIN_PI_MIGRATE=1: pi_accounts_migrate.sh moves it aside first, and
# DELFIN_PI_MIGRATE=back sudo bash pi_bootstrap.sh is the way back.
# ============================================================================

AGENT_USER="${AGENT_USER:-delfin}"
# NO_SUDO=1 (the earlier switch) still means no sudo at all.
AGENT_SUDO="${AGENT_SUDO:-$([[ "${NO_SUDO:-0}" == "1" ]] && echo none || echo service)}"
PY_VERSION="${PY_VERSION:-3.11}"
ACCOUNTS="${DELFIN_PI_ACCOUNTS:-0}"
SHARED_ROOT="${DELFIN_PI_SHARED_ROOT:-/opt/delfin}"
MIGRATE="${DELFIN_PI_MIGRATE:-0}"
HERE="$(cd "$(dirname "$0")" && pwd)"

if [[ $EUID -ne 0 ]]; then
    echo "error: run with sudo" >&2
    exit 1
fi
[[ "$MIGRATE" == 0 || "$MIGRATE" == 1 || "$MIGRATE" == back ]] \
    || { echo "error: DELFIN_PI_MIGRATE must be 1 or back (got '$MIGRATE')" >&2; exit 1; }
[[ "$MIGRATE" != 1 || "$ACCOUNTS" == 1 ]] \
    || { echo "error: DELFIN_PI_MIGRATE=1 moves to the per-person variant: set DELFIN_PI_ACCOUNTS=1 too" >&2; exit 1; }
if [[ "$ACCOUNTS" == "1" || "$MIGRATE" == back ]]; then
    [[ "$SHARED_ROOT" =~ ^/[A-Za-z0-9._/-]+$ && "$SHARED_ROOT" != / ]]         || { echo "error: DELFIN_PI_SHARED_ROOT must be an absolute path (got '$SHARED_ROOT')" >&2; exit 1; }
    for f in delfin_accounts.sh delfin_people.sh delfin_deploy.sh pi_accounts_migrate.sh; do
        [[ -f "$HERE/$f" ]] || { echo "error: $HERE/$f missing: run this script from scripts/pi (deploy_pi.sh --bootstrap copies it)" >&2; exit 1; }
    done
fi
# The way back from the per-person variant to the one-account layout: only
# the moves, nothing installed (the next deploy_pi.sh repoints the starters).
if [[ "$MIGRATE" == back ]]; then
    back_home="$(getent passwd "$AGENT_USER" | cut -d: -f6)"
    [[ -n "$back_home" ]] || { echo "error: no user $AGENT_USER" >&2; exit 1; }
    bash "$HERE/pi_accounts_migrate.sh" back "$back_home" "$SHARED_ROOT"
    exit 0
fi
if [[ "$(uname -m)" != "aarch64" ]]; then
    echo "warning: expected aarch64, got $(uname -m)" >&2
fi

LOGIN_USER="${SUDO_USER:-}"

echo "==> [1/5] apt packages"
export DEBIAN_FRONTEND=noninteractive
apt-get update -q
apt-get install -y -q \
    git curl ca-certificates build-essential cmake pkg-config tmux htop rsync jq \
    libxrender1 libxext6 libsm6 libgl1 libglib2.0-0 libfreetype6 \
    openbabel openmpi-bin libopenmpi-dev bubblewrap nodejs npm gh \
    tesseract-ocr tesseract-ocr-eng tesseract-ocr-deu \
    libglew-dev libpng-dev libfreetype-dev libxml2-dev libglm-dev libnetcdf-dev

echo "==> [2/5] agent user '$AGENT_USER'"
if ! id "$AGENT_USER" >/dev/null 2>&1; then
    adduser --disabled-password --gecos "DELFIN agent" "$AGENT_USER"
fi
SUDOERS="/etc/sudoers.d/90-delfin-agent"
# What the agent user may do as root. Not "the Pi is the sandbox": what runs
# outside DELFIN's command cage (a CLI backend's shell, MCP servers, the
# agent process) runs as this user, and NOPASSWD:ALL made any approved
# `sudo ...` root. Administration is the login user's (pi).
#   service (default)  restart / query DELFIN's own services (dashboard, agent
#                      service, messenger bridge and gateway), nothing else
#   none               no sudo at all
#   full               NOPASSWD:ALL -- only when asked for, and said so
case "$AGENT_SUDO" in
    none)
        rm -f "$SUDOERS" ;;
    full)
        echo "warning: AGENT_SUDO=full: the agent user becomes root without a password" >&2
        printf '%s ALL=(ALL) NOPASSWD:ALL\n' "$AGENT_USER" > "$SUDOERS.tmp" ;;
    service)
        systemctl_bin="$(command -v systemctl)"
        # A deploy swaps the code; only a restart makes a running service use
        # it. With the dashboard alone allowed, the agent service and the
        # messenger kept yesterday's code after every deploy.
        rules=""
        for unit in delfin-dashboard delfin-agent-service delfin-messenger-bridge delfin-messenger-gateway; do
            rules+="${rules:+, }$systemctl_bin restart $unit, $systemctl_bin is-active $unit"
        done
        printf '%s ALL=(root) NOPASSWD: %s\n' "$AGENT_USER" "$rules" > "$SUDOERS.tmp" ;;
    *)
        echo "error: AGENT_SUDO must be service, none or full (got '$AGENT_SUDO')" >&2
        exit 1 ;;
esac
if [[ -f "$SUDOERS.tmp" ]]; then
    chmod 0440 "$SUDOERS.tmp"
    visudo -cf "$SUDOERS.tmp" >/dev/null
    mv "$SUDOERS.tmp" "$SUDOERS"
fi
AGENT_HOME="$(getent passwd "$AGENT_USER" | cut -d: -f6)"

# One pinned uv release, kept only when its SHA-256 matches (never curl | sh).
# install_uv DIR: uv and uvx into DIR (mode 755), owned by whoever runs it.
UV_VERSION=0.13.5
UV_SHA256_aarch64=4bb69cf95ed13f23713bb7bc09f241896ea0c75b7c557e615a20604b82f9c403
UV_SHA256_x86_64=f68d291c75783122dad58ab77def0e13f3233f358731db5d6d7e8715603b928f
install_uv() {
    local dir="$1" arch want_var want tmp got
    arch="$(uname -m)"
    want_var="UV_SHA256_$arch"
    want="${!want_var:-}"
    [[ -n "$want" ]] || { echo "error: no pinned uv digest for $arch" >&2; return 1; }
    tmp="$(mktemp -d)"
    curl -fsSL -o "$tmp/uv.tar.gz" \
        "https://github.com/astral-sh/uv/releases/download/$UV_VERSION/uv-$arch-unknown-linux-gnu.tar.gz"
    got="$(sha256sum "$tmp/uv.tar.gz" | cut -d' ' -f1)"
    if [[ "$got" != "$want" ]]; then
        echo "error: uv digest mismatch: got $got, want $want" >&2
        rm -rf "$tmp"; return 1
    fi
    tar -xzf "$tmp/uv.tar.gz" -C "$tmp"
    mkdir -p "$dir"
    install -m 755 "$tmp/uv-$arch-unknown-linux-gnu/uv" "$tmp/uv-$arch-unknown-linux-gnu/uvx" "$dir/"
    rm -rf "$tmp"
}

# One account per person: the trust boundary.
#   The agent account ($AGENT_USER) runs the agent; nothing it can write may
#   be executed or loaded as a setting by a person's account. So the shared
#   installation is NOT the agent account's:
#   - $SHARED_ROOT itself, agent.env and settings-defaults.json: root;
#   - what a deploy builds (src, .venv, python, ms-playwright, .install):
#     the system account $INSTALL_USER (no login, no password, no group a
#     person or the agent is in), driven by root through delfin-deploy;
#   - $SHARED_ROOT/work: root's link to the agent's own workspace in its
#     home, so every ~/delfin/work path keeps working.
#   Readable by the group delfin-shared (the persons and the agent account),
#   nothing for any other account. Deploys go through the admin
#   (deploy_pi.sh --shared pi@<pi> -> sudo delfin-deploy), never the agent.
INSTALL_USER="${DELFIN_PI_INSTALL_USER:-delfin-install}"
INSTALL_HOME="/var/lib/$INSTALL_USER"
UV_SHARED=/usr/local/lib/delfin/bin/uv
if [[ "$ACCOUNTS" == "1" ]]; then
    echo "==> one account per person: shared installation in $SHARED_ROOT (root, built by $INSTALL_USER)"
    if [[ -e "$AGENT_HOME/delfin" && ! -L "$AGENT_HOME/delfin" ]]; then
        if [[ "$MIGRATE" != 1 ]]; then
            echo "error: $AGENT_HOME/delfin is a directory: this box already has the one-account layout." >&2
            echo "       DELFIN_PI_MIGRATE=1 moves it aside first (pi_accounts_migrate.sh; nothing is" >&2
            echo "       deleted, DELFIN_PI_MIGRATE=back returns). Nothing was changed." >&2
            exit 1
        fi
        bash "$HERE/pi_accounts_migrate.sh" to "$AGENT_HOME" "$SHARED_ROOT"
    fi
    getent group delfin-people >/dev/null || groupadd --system delfin-people
    getent group delfin-shared >/dev/null || groupadd --system delfin-shared
    if ! id "$INSTALL_USER" >/dev/null 2>&1; then
        useradd --system --home-dir "$INSTALL_HOME" --no-create-home \
            --shell /usr/sbin/nologin --user-group --comment "DELFIN shared installation (builds only)" \
            "$INSTALL_USER"
    fi
    install -d -m 700 -o "$INSTALL_USER" -g "$INSTALL_USER" "$INSTALL_HOME"
    # Both read the shared root (root's, 0750): the agent account to run it,
    # the install account to reach the directories it builds in.
    agent_new_reader=0
    id -nG "$AGENT_USER" | tr ' ' '\n' | grep -qx delfin-shared || agent_new_reader=1
    usermod -aG delfin-shared "$AGENT_USER"
    usermod -aG delfin-shared "$INSTALL_USER"
    # Root's directory: the agent cannot add, rename or replace anything in it.
    install -d -m 750 -o root -g delfin-shared "$SHARED_ROOT"
    for d in src .venv python ms-playwright .install; do
        [[ -L "$SHARED_ROOT/$d" ]] && rm -f "$SHARED_ROOT/$d"
        install -d -m 755 -o "$INSTALL_USER" -g "$INSTALL_USER" "$SHARED_ROOT/$d"
    done
    # The agent's own workspace, in its home; the link is root's.
    install -d -m 700 -o "$AGENT_USER" -g "$AGENT_USER" "$AGENT_HOME/delfin-work"
    if [[ -e "$SHARED_ROOT/work" && ! -L "$SHARED_ROOT/work" ]]; then
        echo "error: $SHARED_ROOT/work is a directory; move it to $AGENT_HOME/delfin-work first" >&2
        exit 1
    fi
    ln -sfn "$AGENT_HOME/delfin-work" "$SHARED_ROOT/work"
    ln -sfn "$SHARED_ROOT" "$AGENT_HOME/delfin"
    chown -h "$AGENT_USER:$AGENT_USER" "$AGENT_HOME/delfin"
    install -d -m 755 /etc/delfin
    printf '%s\n' "# One Linux account per person (pi_bootstrap.sh, delfin-accounts, delfin-deploy)." \
        "DELFIN_SHARED_ROOT=$SHARED_ROOT" "DELFIN_DEPLOY_USER=$AGENT_USER" \
        "DELFIN_INSTALL_USER=$INSTALL_USER" > /etc/delfin/accounts.conf.tmp
    chmod 644 /etc/delfin/accounts.conf.tmp
    mv /etc/delfin/accounts.conf.tmp /etc/delfin/accounts.conf
    # The box's agent and tab defaults for every person: root's files, read
    # by delfin-deploy into settings-defaults.json. Never the agent account's
    # settings: it must not choose a person's permission profile.
    [[ -f /etc/delfin/pi-agent.conf ]] || printf '%s\n' \
        "# Every person's agent tab defaults (key=value); delfin-deploy writes them into" \
        "# $SHARED_ROOT/settings-defaults.json. mode: solo = Code; permission_profile:" \
        "# ask_all | repo_free | all_free (= Bypass). effort empty = the model's default." \
        "mode=solo" "provider=kit" "model=kit.glm-5.3-flash" \
        "permission_profile=ask_all" "effort=" > /etc/delfin/pi-agent.conf
    [[ -f /etc/delfin/pi-tabs.conf ]] || printf '%s\n' \
        "# Dashboard tabs every person sees (ids; Settings is always shown)." \
        "agent home archive office" > /etc/delfin/pi-tabs.conf
    chmod 644 /etc/delfin/pi-agent.conf /etc/delfin/pi-tabs.conf
    # agent.env of the shared installation: root's. An existing one (a move
    # from the one-account layout) is kept as root's copy.
    if [[ ! -f "$SHARED_ROOT/agent.env" ]]; then
        cat > "$SHARED_ROOT/agent.env.tmp" <<'EOF'
# Shared DELFIN installation (root's file; delfin-deploy keeps it current).
# auto: a command the dashboard runs for the agent is confined (bwrap, else
# Landlock); off would run it with the account's full rights.
export DELFIN_AGENT_SANDBOX=auto
export DELFIN_AGENT_SANDBOX_TIMEOUT=600
# DELFIN's own MCP servers run in bwrap, seeing only their folders.
export DELFIN_MCP_ISOLATION=builtin
EOF
        mv "$SHARED_ROOT/agent.env.tmp" "$SHARED_ROOT/agent.env"
    fi
    grep -q '^export PLAYWRIGHT_BROWSERS_PATH=' "$SHARED_ROOT/agent.env" || printf '%s\n' \
        "# One Chromium for every account (shared installation)." \
        "export PLAYWRIGHT_BROWSERS_PATH=$SHARED_ROOT/ms-playwright" >> "$SHARED_ROOT/agent.env"
    chown root:root "$SHARED_ROOT/agent.env"
    chmod 644 "$SHARED_ROOT/agent.env"
    # Root-owned tools: the admin runs delfin-accounts and delfin-deploy as
    # root, an invitation's forced command runs delfin-people as the person,
    # and the shared build runs this uv. A copy in the agent's tree could be
    # changed by the agent first.
    install -D -m 755 -o root -g root "$HERE/delfin_accounts.sh" /usr/local/sbin/delfin-accounts
    install -D -m 755 -o root -g root "$HERE/delfin_deploy.sh" /usr/local/sbin/delfin-deploy
    install -D -m 755 -o root -g root "$HERE/delfin_people.sh" /usr/local/lib/delfin/delfin-people
    ln -sfn /usr/local/lib/delfin/delfin-people /usr/local/bin/delfin-people
    if [[ ! -x "$UV_SHARED" ]] || ! "$UV_SHARED" --version 2>/dev/null | grep -q " $UV_VERSION"; then
        install_uv "$(dirname "$UV_SHARED")"
        chown root:root "$(dirname "$UV_SHARED")"/uv "$(dirname "$UV_SHARED")"/uvx
    fi
fi

# The administrator's floor (managed_policy), root-owned so neither the agent
# nor its user can loosen it: shell commands only ever in the cage (bwrap),
# and no session on a CLI backend, whose own shell runs outside the cage
# with the agent user's files. CLI_BACKENDS=allowed leaves those to the user.
install -d -m 755 /etc/delfin
if [[ "${CLI_BACKENDS:-refused}" == "allowed" ]]; then
    backends='[]'
else
    backends='["cli", "codex_cli"]'
fi
cat > /etc/delfin/managed_policy.json.tmp <<EOF
{
  "comment": "DELFIN agent box (pi_bootstrap.sh): isolation required; no CLI backend (its shell runs outside the cage).",
  "require_isolation": true,
  "container_backend": "bwrap",
  "disallowed_backends": $backends
}
EOF
chmod 644 /etc/delfin/managed_policy.json.tmp
mv /etc/delfin/managed_policy.json.tmp /etc/delfin/managed_policy.json

echo "==> [3/5] ssh keys"
install -d -m 700 -o "$AGENT_USER" -g "$AGENT_USER" "$AGENT_HOME/.ssh"
if [[ -n "$LOGIN_USER" ]]; then
    src="$(getent passwd "$LOGIN_USER" | cut -d: -f6)/.ssh/authorized_keys"
    if [[ -f "$src" ]]; then
        touch "$AGENT_HOME/.ssh/authorized_keys"
        sort -u "$src" "$AGENT_HOME/.ssh/authorized_keys" -o "$AGENT_HOME/.ssh/authorized_keys"
        chown "$AGENT_USER:$AGENT_USER" "$AGENT_HOME/.ssh/authorized_keys"
        chmod 600 "$AGENT_HOME/.ssh/authorized_keys"
    else
        echo "warning: $src not found - copy your key with ssh-copy-id first" >&2
    fi
fi

echo "==> [4/5] uv + Python $PY_VERSION + venv"
sudo -iu "$AGENT_USER" bash -s -- "$PY_VERSION" "$ACCOUNTS" <<'EOS'
set -euo pipefail
PY_VERSION="$1"
ACCOUNTS="$2"
export PATH="$HOME/.local/bin:$PATH"
# ~/delfin may be a link (one account per person): everything below that
# records a path (the venv, its Python) records the real one.
ROOT="$(readlink -f "$HOME/delfin" 2>/dev/null || echo "$HOME/delfin")"
# One pinned uv release, kept only when its SHA-256 matches (never curl | sh).
UV_VERSION=0.13.5
UV_SHA256_aarch64=4bb69cf95ed13f23713bb7bc09f241896ea0c75b7c557e615a20604b82f9c403
UV_SHA256_x86_64=f68d291c75783122dad58ab77def0e13f3233f358731db5d6d7e8715603b928f
if ! command -v uv >/dev/null; then
    arch="$(uname -m)"
    want_var="UV_SHA256_$arch"
    want="${!want_var:-}"
    [[ -n "$want" ]] || { echo "error: no pinned uv digest for $arch" >&2; exit 1; }
    tmp="$(mktemp -d)"
    curl -fsSL -o "$tmp/uv.tar.gz" \
        "https://github.com/astral-sh/uv/releases/download/$UV_VERSION/uv-$arch-unknown-linux-gnu.tar.gz"
    got="$(sha256sum "$tmp/uv.tar.gz" | cut -d' ' -f1)"
    if [[ "$got" != "$want" ]]; then
        echo "error: uv digest mismatch: got $got, want $want" >&2
        rm -rf "$tmp"; exit 1
    fi
    tar -xzf "$tmp/uv.tar.gz" -C "$tmp"
    mkdir -p "$HOME/.local/bin"
    install -m 755 "$tmp/uv-$arch-unknown-linux-gnu/uv" "$tmp/uv-$arch-unknown-linux-gnu/uvx" "$HOME/.local/bin/"
    rm -rf "$tmp"
fi
uv python install "$PY_VERSION"
# The shared installation (one account per person) is built by the install
# account below, never by this one.
if [[ "$ACCOUNTS" != "1" ]]; then
    mkdir -p "$ROOT" "$ROOT/work"
    [[ -x "$ROOT/.venv/bin/python" ]] || uv venv --python "$PY_VERSION" "$ROOT/.venv"
fi
# CLI backends of the agent (backend=cli). Login happens interactively later.
# Download first: `curl | bash </dev/null` would feed bash /dev/null, not the script.
if ! command -v claude >/dev/null && [[ ! -x "$HOME/.local/bin/claude" ]]; then
    # Into a directory of our own: a fixed name in /tmp could be created or
    # swapped by another local user between the download and the run.
    inst="$(mktemp -d)"
    curl -fsSL https://claude.ai/install.sh -o "$inst/install.sh" \
        && bash "$inst/install.sh" </dev/null \
        || echo "warning: Claude Code CLI install failed" >&2
    rm -rf "$inst"
fi
npm config set prefix "$HOME/.local" </dev/null
command -v codex >/dev/null || npm install -g @openai/codex </dev/null \
    || echo "warning: Codex CLI install failed" >&2
EOS

if [[ "$ACCOUNTS" == "1" ]]; then
    # Python and the venv of the shared installation, as the install account
    # with root's uv: its own Python under $SHARED_ROOT/python (readable by
    # the group, unlike a home), the venv's python a link to it.
    echo "==> shared Python $PY_VERSION + venv in $SHARED_ROOT (as $INSTALL_USER)"
    runuser -u "$INSTALL_USER" -- env -i HOME="$INSTALL_HOME" PATH=/usr/local/bin:/usr/bin:/bin \
        LANG=C.UTF-8 UV_PYTHON_INSTALL_DIR="$SHARED_ROOT/python" UV_CACHE_DIR="$INSTALL_HOME/.cache/uv" \
        bash -c 'cd "$HOME" && umask 022 && "$0" python install "$1" && { [ -x "$2/.venv/bin/python" ] || "$0" venv --allow-existing --python "$1" "$2/.venv"; }' \
        "$UV_SHARED" "$PY_VERSION" "$SHARED_ROOT"
fi

echo "==> [5/5] agent env + wrapper"
sudo -iu "$AGENT_USER" bash -s <<'EOS'
set -euo pipefail
ENV_FILE="$HOME/delfin/agent.env"
if [[ ! -f "$ENV_FILE" && ! -L "$HOME/delfin" ]]; then
    cat > "$ENV_FILE" <<'EOF'
# Sourced by ~/.local/bin/delfin-agent and the dashboard service. Edit freely.
# auto: a command the dashboard runs for the agent is confined (bwrap, else
# Landlock); off would run it with the account's full rights.
export DELFIN_AGENT_SANDBOX=auto
export DELFIN_AGENT_SANDBOX_TIMEOUT=600
# DELFIN's own MCP servers (delfin-tools, delfin-ops) run in bwrap, seeing
# only the folders they work in -- not ~/.ssh, not the credentials.
export DELFIN_MCP_ISOLATION=builtin
EOF
fi
mkdir -p "$HOME/.local/bin"
cat > "$HOME/.local/bin/delfin-agent" <<'EOF'
#!/usr/bin/env bash
# delfin-agent <args> inside the venv, always the deployed code: `python -m`
# would put the cwd first on sys.path and run a checkout's delfin instead
# (e.g. inside ~/delfin/work/DELFIN). Hence the entry point, else -P.
set -euo pipefail
source "$HOME/delfin/agent.env"
source "$HOME/delfin/.venv/bin/activate"
[[ -x "$HOME/delfin/.venv/bin/delfin-agent" ]] && exec "$HOME/delfin/.venv/bin/delfin-agent" "$@"
exec python -P -m delfin.agent.cli "$@"
EOF
chmod +x "$HOME/.local/bin/delfin-agent"
grep -q 'delfin/.venv/bin/activate' "$HOME/.bashrc" || cat >> "$HOME/.bashrc" <<'EOF'
export PATH="$HOME/.local/bin:$PATH"
source "$HOME/delfin/agent.env"
[[ -f "$HOME/delfin/.venv/bin/activate" ]] && source "$HOME/delfin/.venv/bin/activate"
EOF
EOS

echo "==> dashboard service (127.0.0.1:8866, starts once code is deployed)"
sudo -iu "$AGENT_USER" bash -s <<'EOS'
set -euo pipefail
umask 077
mkdir -p "$HOME/.delfin"
[[ -s "$HOME/.delfin/dashboard.token" ]] \
    || "$HOME/delfin/.venv/bin/python" -c "import secrets; print(secrets.token_urlsafe(32))" \
       > "$HOME/.delfin/dashboard.token"
umask 022
cat > "$HOME/.local/bin/delfin-dashboard" <<'EOF'
#!/usr/bin/env bash
# Run the dashboard on loopback with the persistent token (used by systemd).
set -euo pipefail
source "$HOME/delfin/agent.env"
source "$HOME/delfin/.venv/bin/activate"
# A service, not a launcher: without this the server ends itself after
# 10 min with no window open (resume_server.NEVER_CONNECTED_SECONDS).
export DELFIN_DASHBOARD_STAY_UP=1
cd "$HOME/delfin/work"
exec delfin-voila --ip 127.0.0.1 --port 8866 --no-browser --token "$(cat "$HOME/.delfin/dashboard.token")"
EOF
chmod +x "$HOME/.local/bin/delfin-dashboard"
EOS
# Enabled once, when the unit is new. Run again, the bootstrap leaves the
# state an admin chose (stopped and disabled until access is secured, say);
# `systemctl enable` on every run switched it back on for the next boot.
dashboard_unit_new=0
[[ -f /etc/systemd/system/delfin-dashboard.service ]] || dashboard_unit_new=1
cat > /etc/systemd/system/delfin-dashboard.service <<EOF
[Unit]
Description=DELFIN dashboard (loopback only)
After=network-online.target
ConditionPathExists=$AGENT_HOME/delfin/src/pyproject.toml

[Service]
User=$AGENT_USER
ExecStart=$AGENT_HOME/.local/bin/delfin-dashboard
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF
cat > /usr/local/bin/delfin-dashboard-url <<EOF
#!/usr/bin/env bash
# Print the dashboard URL incl. token (any sudo user may call this).
echo "http://127.0.0.1:8866/?token=\$(sudo cat $AGENT_HOME/.delfin/dashboard.token)"
EOF
chmod 755 /usr/local/bin/delfin-dashboard-url
cat > /usr/share/applications/delfin-dashboard.desktop <<'EOF'
[Desktop Entry]
Type=Application
Name=DELFIN Dashboard
Comment=DELFIN dashboard incl. agent (local)
Exec=sh -c 'x-www-browser "$(delfin-dashboard-url)"'
Icon=applications-science
Categories=Science;Education;
EOF
# Lab control: the e-stop owner and the dead-man are systemd user units of
# the agent user (pi_install.sh installs them); lingering lets them run
# without a login. GPIO (libgpiod, /dev/gpiochip*, group gpio via the
# Raspberry Pi OS udev rule) only when asked for: the group grants the whole
# chip, so it is an explicit choice (docs/design/lab-control.md section 9).
loginctl enable-linger "$AGENT_USER" || echo "warning: lingering not enabled for $AGENT_USER" >&2
# A user manager started before the agent account joined delfin-shared keeps
# its old groups, and so do the lab units it starts: they could not read the
# shared installation. It is restarted once, when the group is new (its
# sessions -- tmux, the lab units -- end and come back through lingering).
if [[ "$ACCOUNTS" == "1" && "${agent_new_reader:-0}" == "1" ]] \
    && systemctl is-active --quiet "user@$(id -u "$AGENT_USER").service"; then
    echo "==> $AGENT_USER joined delfin-shared: restarting its user manager (lab units, tmux sessions)"
    systemctl restart "user@$(id -u "$AGENT_USER").service" \
        || echo "warning: restart it yourself: sudo systemctl restart user@$(id -u "$AGENT_USER").service" >&2
fi
if [[ "${LAB_GPIO:-0}" == "1" ]] && getent group gpio >/dev/null; then
    usermod -aG gpio "$AGENT_USER"
fi
systemctl daemon-reload
if [[ "$dashboard_unit_new" == "1" ]]; then
    systemctl enable delfin-dashboard.service
else
    echo "    dashboard service left as it was: $(systemctl is-enabled delfin-dashboard.service 2>/dev/null)"
fi

echo
echo "Bootstrap done. Agent user: $AGENT_USER (sudo: $AGENT_SUDO)"
if [[ "$ACCOUNTS" == "1" ]]; then
    echo "One account per person: installation in $SHARED_ROOT (root's, built by $INSTALL_USER)."
    echo "Next (from your PC), as the admin -- the agent account cannot deploy it:"
    echo "  scripts/pi/deploy_pi.sh --shared ${LOGIN_USER:-pi}@<pi-host>"
    echo "then: sudo delfin-accounts add NAME    (or: delfin-people --add-person NAME)"
else
    echo "Next (from your PC): scripts/pi/deploy_pi.sh $AGENT_USER@<pi-host>"
fi
