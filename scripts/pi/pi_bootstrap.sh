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
#   2. Create a dedicated user (default: delfin) with passwordless sudo, so the
#      agent can administer the box. The Pi itself is the sandbox.
#   3. Copy the invoking user's authorized_keys to that user
#   4. Install uv + Python 3.11 for that user and create ~/delfin/.venv
#   5. Write ~/delfin/agent.env and the ~/.local/bin/delfin-agent wrapper
#
# Code is NOT installed here - deploy_pi.sh pushes it to ~/delfin/src.
#
# Usage (on the Pi):
#   sudo bash pi_bootstrap.sh            # full sudo for the agent user
#   NO_SUDO=1 sudo bash pi_bootstrap.sh  # agent user without sudo
# ============================================================================

AGENT_USER="${AGENT_USER:-delfin}"
NO_SUDO="${NO_SUDO:-0}"
PY_VERSION="${PY_VERSION:-3.11}"

if [[ $EUID -ne 0 ]]; then
    echo "error: run with sudo" >&2
    exit 1
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
if [[ "$NO_SUDO" == "1" ]]; then
    rm -f "$SUDOERS"
else
    printf '%s ALL=(ALL) NOPASSWD:ALL\n' "$AGENT_USER" > "$SUDOERS.tmp"
    chmod 0440 "$SUDOERS.tmp"
    visudo -cf "$SUDOERS.tmp" >/dev/null
    mv "$SUDOERS.tmp" "$SUDOERS"
fi
AGENT_HOME="$(getent passwd "$AGENT_USER" | cut -d: -f6)"

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
sudo -iu "$AGENT_USER" bash -s -- "$PY_VERSION" <<'EOS'
set -euo pipefail
PY_VERSION="$1"
export PATH="$HOME/.local/bin:$PATH"
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
mkdir -p "$HOME/delfin" "$HOME/delfin/work"
[[ -x "$HOME/delfin/.venv/bin/python" ]] || uv venv --python "$PY_VERSION" "$HOME/delfin/.venv"
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

echo "==> [5/5] agent env + wrapper"
sudo -iu "$AGENT_USER" bash -s <<'EOS'
set -euo pipefail
ENV_FILE="$HOME/delfin/agent.env"
if [[ ! -f "$ENV_FILE" ]]; then
    cat > "$ENV_FILE" <<'EOF'
# Sourced by ~/.local/bin/delfin-agent. Edit freely.
# The Pi is the sandbox: the agent's command runner needs no extra jail.
export DELFIN_AGENT_SANDBOX=off
export DELFIN_AGENT_SANDBOX_TIMEOUT=600
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
systemctl daemon-reload
systemctl enable delfin-dashboard.service

echo
echo "Bootstrap done. Agent user: $AGENT_USER ($([[ "$NO_SUDO" == "1" ]] && echo 'no sudo' || echo 'passwordless sudo'))"
echo "Next (from your PC): scripts/pi/deploy_pi.sh $AGENT_USER@<pi-host>"
