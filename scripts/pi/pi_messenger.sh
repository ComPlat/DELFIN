#!/usr/bin/env bash
# ============================================================================
# Provision the messenger gateway host: SimpleX daemon + DELFIN gateway.
#
#   sudo bash pi_messenger.sh [--whisper] [--pair] [--no-phone-sessions]
#                             [--no-nftables] [--uninstall]
#
# Works on any Linux with systemd and glibc >= 2.34 the pinned releases
# cover -- a Raspberry Pi (aarch64, Raspberry Pi OS or Ubuntu), a login
# node or a compute frontend (x86_64): Ubuntu 24.04 takes its own build,
# every other such system the Ubuntu 22.04 build (it needs glibc 2.34,
# OpenSSL 3, libgmp and zlib -- ldd/objdump of the pinned v7.0.3 binary).
#
# Does, idempotently:
#   - a system user `delfin-gw`, no sudo, no login shell;
#   - the pinned simplex-chat binary for this platform under
#     /opt/delfin-messenger/bin, owned by delfin-gw, 0755;
#   - bot state 0700 under /var/lib/delfin-messenger; the daemon creates
#     its bot profile "DELFIN" on first start and keeps it after;
#   - the gateway venv (/opt/delfin-messenger/bin/venv, Python 3.11 from
#     MESSENGER_PYTHON, uv or python3.11): websockets from
#     messenger-requirements.txt with --require-hashes, and a root-owned
#     copy of the delfin package next to this script (no build step);
#   - the agent service as a unit of its own (delfin-agent-service), run
#     with the service user's DELFIN venv and agent.env, so the bridge has
#     a service without an open terminal; sessions the phone starts
#     (neu NAME) get folders under ~<service user>/DELFIN-Handy
#     (--no-phone-sessions: none);
#   - with --pair: the pairing invitation for the phone, at the end;
#   - hardened systemd units simplex-chat-daemon.service and
#     delfin-messenger-gateway.service (the gateway unit only when the
#     DELFIN package is installed for the `delfin-gw` user -- on a server
#     whose DELFIN service runs as another user, the gateway must run in
#     that user's own service instead; see README "Messenger");
#   - an nftables rule letting only `meta skuid delfin-gw` reach the
#     daemon's local API port (skipped with --no-nftables, e.g. when the
#     host's firewall is managed elsewhere);
#   - with --whisper: local speech recognition for voice messages --
#     faster-whisper in a venv of its own (${PREFIX}/whisper-venv; it needs
#     numpy 2, DELFIN keeps numpy below 2), installed --require-hashes
#     from whisper-requirements-<arch>.txt, and the faster-whisper-small
#     model at a pinned revision, every file checked against its SHA-256,
#     under ${PREFIX}/models (root-owned, readable). Nothing is fetched at
#     run time; the gateway transcribes offline;
#   - nothing else. No fixed /tmp names: every temporary file comes from
#     mktemp and runs nothing from it.
#
# The binary is downloaded once, checked against its pinned SHA-256 below,
# and only installed when it matches. A mirror may be used with
# SIMPLEX_URL, but the pinned digest of this platform still decides.
# ============================================================================
set -euo pipefail

SIMPLEX_VERSION="v7.0.3"
# From the release notes of v7.0.3 (github.com/simplex-chat/simplex-chat,
# release page, "Desktop app" section), read 2026-10-25. One digest per
# (ubuntu version, machine); a machine with no digest here is refused.
SIMPLEX_SHA256_24_04_x86_64="895fb14cfaa662d1c0947f2f871141fc89680fe3e70bff64c366cbe9e59aa4f0"
SIMPLEX_SHA256_24_04_aarch64="2d2e62351f11bc51ae659618584722b38ea5a6796b806c9a388ee6e3124c90ac"
SIMPLEX_SHA256_22_04_x86_64="5afb1d25efe5ccf564a1ab124bc7f410a7a73171c974a4d8b7f4d8f2e3d62e77"
SIMPLEX_SHA256_22_04_aarch64="91f332edc813cab78a848da7ba75ed91bcb8d09abf1fc92eca74330dc29ec970"

GW_USER="delfin-gw"
# Where the systemd units live; overridable for tests and for systems whose
# unit dir lies elsewhere (a bind mount, a container image).
UNIT_DIR="${MESSENGER_UNIT_DIR:-/etc/systemd/system}"
PREFIX="${MESSENGER_PREFIX:-/opt/delfin-messenger}"
STATE_DIR="${MESSENGER_STATE_DIR:-/var/lib/delfin-messenger}"
SERVICE_USER="${MESSENGER_SERVICE_USER:-delfin}"
SERVICE_HOME="${MESSENGER_SERVICE_HOME:-$(getent passwd "${SERVICE_USER}" 2>/dev/null | cut -d: -f6)}"
# The bridge and the agent service run the service user's own DELFIN (the
# Pi's ~/delfin/.venv), so both speak the service's protocol; without it
# the bridge falls back to the gateway venv and no service unit is written.
USER_PYTHON="${SERVICE_HOME:+${SERVICE_HOME}/delfin/.venv/bin/python}"
if [ -n "${MESSENGER_BRIDGE_PYTHON:-}" ]; then
  BRIDGE_PYTHON="${MESSENGER_BRIDGE_PYTHON}"
elif [ -n "${USER_PYTHON}" ] && [ -x "${USER_PYTHON}" ]; then
  BRIDGE_PYTHON="${USER_PYTHON}"
else
  BRIDGE_PYTHON="${PREFIX}/bin/venv/bin/python"
fi
PHONE_ROOT="${SERVICE_HOME}/DELFIN-Handy"
BRIDGE_SOCKET="/run/delfin-messenger-bridge/messenger.sock"
API_PORT="5227"          # simplex-chat daemon local API (--chat-server-port)
NFT_TABLE="delfin_messenger"

# Local speech recognition (--whisper). The model is faster-whisper-small
# at one revision; each file is checked against its SHA-256 before install.
WHISPER_REPO="Systran/faster-whisper-small"
WHISPER_REVISION="536b0662742c02347bc0e980a01041f333bce120"
WHISPER_FILES="config.json:b55496ac7940a7ae47d2c01eab40edfd8701feec1229d9cce3b40014383fb828
model.bin:3e305921506d8872816023e4c273e75d2419fb89b24da97b4fe7bce14170d671
tokenizer.json:fb7b63191e9bb045082c79fd742a3106a12c99513ab30df4a0d47fa6cb6fd0ab
vocabulary.txt:34ce3fe1c5041027b3f8d42912270993f986dbc4bb34cf27f951e34a1e453913"
WHISPER_DIR="${PREFIX}/models/faster-whisper-small"
WHISPER_VENV="${PREFIX}/whisper-venv"
HERE="$(cd "$(dirname "$0")" && pwd)"

say() { printf '[messenger] %s\n' "$*" >&2; }
die() { printf '[messenger] ERROR: %s\n' "$*" >&2; exit 1; }

action=install
want_nftables=1
want_whisper=0
want_pair=0
want_phone_sessions=1
for arg in "$@"; do
  case "$arg" in
    --uninstall) action=uninstall ;;
    --no-nftables) want_nftables=0 ;;
    --whisper) want_whisper=1 ;;
    --pair) want_pair=1 ;;
    --no-phone-sessions) want_phone_sessions=0 ;;
    *) die "unknown argument: $arg" ;;
  esac
done

[ "$(id -u)" = "0" ] || die "run as root"

# --- platform ---------------------------------------------------------------
machine="$(uname -m 2>/dev/null)" || die "uname -m failed"
os_id="$(. "${MESSENGER_OS_RELEASE:-/etc/os-release}" 2>/dev/null; printf '%s' "${ID:-}")" || os_id=""
os_version="$(. "${MESSENGER_OS_RELEASE:-/etc/os-release}" 2>/dev/null; printf '%s' "${VERSION_ID:-}")" || os_version=""
glibc="$(getconf GNU_LIBC_VERSION 2>/dev/null | awk '{print $2}')" || glibc=""
glibc_major="${glibc%%.*}"
glibc_minor="${glibc#*.}"
glibc_minor="${glibc_minor%%.*}"
build=""
if [ "${os_id:-ubuntu}" = "ubuntu" ] && [ "${os_version%%.*}" = "24" ]; then
  build="24_04"
elif [ "${glibc_major}" = 2 ] && [ "${glibc_minor}" -ge 34 ] 2>/dev/null; then
  build="22_04"   # needs glibc 2.34, OpenSSL 3, libgmp, zlib
fi
digest_var="SIMPLEX_SHA256_${build}_${machine}"
SIMPLEX_SHA="${!digest_var:-}"
{ [ -n "${build}" ] && [ -n "${SIMPLEX_SHA}" ]; } || die "no pinned simplex-chat ${SIMPLEX_VERSION} digest for ${machine} / ${os_id:-unknown} ${os_version} (glibc ${glibc:-unknown}, needs >= 2.34); refusing to install anything"

asset="simplex-chat-ubuntu-${build}-${machine}"
case "${machine}" in
  x86_64|aarch64) : ;;
  *) die "unsupported machine: ${machine}" ;;
esac

# --- user -------------------------------------------------------------------
ensure_user() {
  if ! id "${GW_USER}" >/dev/null 2>&1; then
    useradd --system --home-dir "${STATE_DIR}" --create-home \
      --shell /usr/sbin/nologin "${GW_USER}"
    say "created user ${GW_USER}"
  fi
}

# --- pinned download --------------------------------------------------------
simplex_url() {
  printf 'https://github.com/simplex-chat/simplex-chat/releases/download/%s/%s\n' \
    "${SIMPLEX_VERSION}" "${asset}"
}

fetch_binary() {
  local dest="$1" url tmp downloaded
  url="${SIMPLEX_URL:-$(simplex_url)}"
  tmp="$(mktemp)"          # 0600, owner root; no fixed name in /tmp
  downloaded=0
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 2 -o "${tmp}" "${url}" && downloaded=1
  elif command -v wget >/dev/null 2>&1; then
    wget -q -O "${tmp}" "${url}" && downloaded=1
  else
    rm -f "${tmp}"
    die "need curl or wget to fetch ${asset}"
  fi
  if [ "${downloaded}" != 1 ] || [ ! -s "${tmp}" ]; then
    rm -f "${tmp}"; die "download failed: ${url}"
  fi
  local actual
  actual="$(sha256sum "${tmp}" | cut -d' ' -f1)"
  if [ "${actual}" != "${SIMPLEX_SHA}" ]; then
    rm -f "${tmp}"
    die "SHA-256 mismatch for ${asset}: got ${actual}, pinned ${SIMPLEX_SHA}"
  fi
  install -m 0755 -o "${GW_USER}" -g "${GW_USER}" "${tmp}" "${dest}"
  rm -f "${tmp}"
  say "installed ${dest} (SHA-256 verified)"
}

# --- systemd units ----------------------------------------------------------
daemon_unit() {
  cat <<EOF
[Unit]
Description=SimpleX chat daemon (local API for the DELFIN messenger gateway)
After=network.target

[Service]
User=${GW_USER}
Group=${GW_USER}
ExecStart=${PREFIX}/bin/simplex-chat -d ${STATE_DIR}/simplex.db --chat-server-port ${API_PORT} --files-folder ${STATE_DIR} --auto-accept-files 0 --create-bot-display-name DELFIN --create-bot-allow-files
Restart=on-failure
RestartSec=5

# hardening
NoNewPrivileges=yes
ProtectSystem=strict
ProtectHome=yes
PrivateTmp=yes
PrivateDevices=yes
ProtectKernelTunables=yes
ProtectKernelModules=yes
ProtectControlGroups=yes
RestrictAddressFamilies=AF_UNIX AF_INET AF_INET6 AF_NETLINK
CapabilityBoundingSet=
LockPersonality=yes
MemoryDenyWriteExecute=yes
RestrictRealtime=yes
RestrictSUIDSGID=yes
SystemCallArchitectures=native
ReadWritePaths=${STATE_DIR}

[Install]
WantedBy=multi-user.target
EOF
}

whisper_env() {
  # Only a model that install_whisper verified and an interpreter for it.
  if [ -f "${WHISPER_DIR}/model.bin" ] && [ -x "${WHISPER_VENV}/bin/python" ]; then
    printf 'Environment=DELFIN_WHISPER_MODEL=%s\n' "${WHISPER_DIR}"
    printf 'Environment=DELFIN_WHISPER_PYTHON=%s\n' "${WHISPER_VENV}/bin/python"
  fi
}

gateway_unit() {
  cat <<EOF
[Unit]
Description=DELFIN messenger gateway (phone <-> agent-service bridge)
After=simplex-chat-daemon.service delfin-messenger-bridge.service
Requires=simplex-chat-daemon.service delfin-messenger-bridge.service

[Service]
User=${GW_USER}
Group=${GW_USER}
Environment=DELFIN_MESSENGER_PORT=${API_PORT}
Environment=DELFIN_MESSENGER_SOCKET=${BRIDGE_SOCKET}
Environment=DELFIN_MESSENGER_DIR=${STATE_DIR}/pairing
$(whisper_env)
ExecStart=${PREFIX}/bin/venv/bin/python -m delfin.messenger run
Restart=on-failure
RestartSec=5

# hardening
NoNewPrivileges=yes
ProtectSystem=strict
ProtectHome=yes
PrivateTmp=yes
PrivateDevices=yes
ProtectKernelTunables=yes
ProtectKernelModules=yes
ProtectControlGroups=yes
RestrictAddressFamilies=AF_UNIX AF_INET AF_INET6 AF_NETLINK
CapabilityBoundingSet=
LockPersonality=yes
RestrictRealtime=yes
RestrictSUIDSGID=yes
SystemCallArchitectures=native
ReadWritePaths=${STATE_DIR}

[Install]
WantedBy=multi-user.target
EOF
}

agent_service_unit() {
  cat <<EOF
[Unit]
Description=DELFIN agent service (the sessions terminal, dashboard and phone share)
After=network-online.target

[Service]
User=${SERVICE_USER}
WorkingDirectory=${SERVICE_HOME}
# agent.env holds the sandbox settings the terminal wrapper sources too.
ExecStart=/bin/bash -c '[ -f ${SERVICE_HOME}/delfin/agent.env ] && . ${SERVICE_HOME}/delfin/agent.env; exec ${BRIDGE_PYTHON} -P -m delfin.agent.service.supervisor --idle-exit 315360000'
# A service a terminal started first holds the host: this one exits 0 then
# and takes over once that one ends.
Restart=always
RestartSec=15
# The service ends its sessions (they are saved) before the rest is killed.
KillMode=mixed
TimeoutStopSec=60

[Install]
WantedBy=multi-user.target
EOF
}

phone_env() {
  if [ "${want_phone_sessions}" = 1 ]; then
    printf 'Environment=DELFIN_MESSENGER_WORKSPACES=%s\n' "${PHONE_ROOT}"
  fi
}

bridge_unit() {
  cat <<EOF
[Unit]
Description=DELFIN agent service bridge for the messenger gateway
After=network.target delfin-agent-service.service
Wants=delfin-agent-service.service

[Service]
User=${SERVICE_USER}
Group=${GW_USER}
RuntimeDirectory=delfin-messenger-bridge
RuntimeDirectoryMode=0750
Environment=DELFIN_MESSENGER_PEER_UID=$(id -u "${GW_USER}")
Environment=DELFIN_MESSENGER_SOCKET=${BRIDGE_SOCKET}
$(phone_env)
ExecStart=${BRIDGE_PYTHON} -m delfin.agent.service.bridge
Restart=on-failure
RestartSec=5
NoNewPrivileges=yes
ProtectSystem=strict
ProtectHome=read-only
# Files from the phone land in a session's workspace, sessions it starts
# get a folder, /sync records live in ~/.delfin: the service home only.
ReadWritePaths=${SERVICE_HOME}
# The agent service socket is in the host's private per-user /tmp directory.
PrivateTmp=no
PrivateDevices=yes
CapabilityBoundingSet=
RestrictAddressFamilies=AF_UNIX

[Install]
WantedBy=multi-user.target
EOF
}

# The daemon must bind its API only on the loopback interface; the nftables
# rule below keeps even local users other than delfin-gw out.
daemon_conf() {
  cat <<EOF
# Written by scripts/messenger/provision_messenger.sh
server_port = ${API_PORT}
# The daemon listens on 127.0.0.1 only (started with the unit above); the
# socket is additionally guarded by nftables meta skuid below.
EOF
}

install_units() {
  install -d -m 0755 -o root -g root "${UNIT_DIR}"
  daemon_unit > "${UNIT_DIR}/simplex-chat-daemon.service"
  if [ -x "${PREFIX}/bin/venv/bin/python" ]; then
    id "${SERVICE_USER}" >/dev/null 2>&1 || die "agent service user does not exist: ${SERVICE_USER}"
    { [ -n "${SERVICE_HOME}" ] && [ -d "${SERVICE_HOME}" ]; } || die "no home directory for ${SERVICE_USER}"
    [ -x "${BRIDGE_PYTHON}" ] || die "bridge Python not executable: ${BRIDGE_PYTHON}"
    if [ "${BRIDGE_PYTHON}" = "${USER_PYTHON}" ]; then
      agent_service_unit > "${UNIT_DIR}/delfin-agent-service.service"
    else
      say "no DELFIN venv at ${SERVICE_HOME}/delfin/.venv -- no agent service unit; start the service as ${SERVICE_USER} yourself"
    fi
    if [ "${want_phone_sessions}" = 1 ]; then
      install -d -m 0700 -o "${SERVICE_USER}" -g "${SERVICE_USER}" "${PHONE_ROOT}"
    fi
    bridge_unit > "${UNIT_DIR}/delfin-messenger-bridge.service"
    gateway_unit > "${UNIT_DIR}/delfin-messenger-gateway.service"
  else
    say "no gateway venv under ${PREFIX}/bin/venv -- skipping the gateway unit (install DELFIN for ${GW_USER} first, then re-run)"
  fi
  systemctl daemon-reload
  systemctl enable --now simplex-chat-daemon.service
  if [ -f "${UNIT_DIR}/delfin-agent-service.service" ]; then
    systemctl enable --now delfin-agent-service.service
  fi
  if [ -f "${UNIT_DIR}/delfin-messenger-bridge.service" ]; then
    systemctl enable --now delfin-messenger-bridge.service
  fi
  [ -f "${UNIT_DIR}/delfin-messenger-gateway.service" ] \
    && systemctl enable --now delfin-messenger-gateway.service || true
  say "units installed and started"
}

# --- nftables ---------------------------------------------------------------
install_nftables() {
  command -v nft >/dev/null 2>&1 || { say "nft not found -- skipping the port rule"; return 0; }
  # One table, one chain, idempotent: recreated on each run.
  nft delete table inet ${NFT_TABLE} 2>/dev/null || true
  nft -f - <<EOF
table inet ${NFT_TABLE} {
  chain output {
    type filter hook output priority 0;
    ip daddr 127.0.0.1 tcp dport ${API_PORT} meta skuid "${GW_USER}" accept
    ip daddr 127.0.0.1 tcp dport ${API_PORT} drop
  }
}
EOF
  say "nftables: ${API_PORT} open only to ${GW_USER}"
}

# --- uninstall --------------------------------------------------------------
uninstall() {
  systemctl disable --now delfin-messenger-gateway.service 2>/dev/null || true
  systemctl disable --now delfin-messenger-bridge.service 2>/dev/null || true
  rm -f "${UNIT_DIR}/delfin-messenger-bridge.service"
  systemctl disable --now delfin-agent-service.service 2>/dev/null || true
  rm -f "${UNIT_DIR}/delfin-agent-service.service"
  systemctl disable --now simplex-chat-daemon.service 2>/dev/null || true
  rm -f "${UNIT_DIR}/delfin-messenger-gateway.service"
  rm -f "${UNIT_DIR}/simplex-chat-daemon.service"
  systemctl daemon-reload
  nft delete table inet ${NFT_TABLE} 2>/dev/null || true
  say "units, nftables table and state left in place: ${STATE_DIR} (delete it by hand if wanted)"
}

# --- the gateway venv -------------------------------------------------------
find_python311() {
  local uv=""
  if [ -n "${MESSENGER_PYTHON:-}" ]; then printf '%s\n' "${MESSENGER_PYTHON}"; return 0; fi
  # uv's Python carries ensurepip; Debian's python3.11 needs python3.11-venv.
  uv="$(command -v uv 2>/dev/null || true)"
  if [ -z "${uv}" ] && [ -n "${SERVICE_HOME}" ] && [ -x "${SERVICE_HOME}/.local/bin/uv" ]; then
    uv="${SERVICE_HOME}/.local/bin/uv"
  fi
  if [ -n "${uv}" ] && UV_PYTHON_INSTALL_DIR="${PREFIX}/python" "${uv}" python install 3.11 >&2; then
    UV_PYTHON_INSTALL_DIR="${PREFIX}/python" "${uv}" python find 3.11 && return 0
  fi
  command -v python3.11
}

install_gateway_venv() {
  local venv="${PREFIX}/bin/venv" lib="${PREFIX}/lib" src py site
  src="$(cd "${HERE}/../.." && pwd)"
  [ -f "${src}/delfin/messenger/gateway.py" ] || die "no DELFIN source next to this script (${src}/delfin)"
  if [ ! -x "${venv}/bin/python" ]; then
    py="$(find_python311)" || die "need Python 3.11: uv, python3.11, or MESSENGER_PYTHON=/path/to/python3.11"
    "${py}" -m venv "${venv}" \
      || die "could not create ${venv} (Debian/Ubuntu: apt install python3.11-venv)"
  fi
  "${venv}/bin/python" -m pip install --quiet --require-hashes --no-deps \
    -r "${HERE}/messenger-requirements.txt" || die "messenger packages did not install"
  # The deployed code, copied root-owned: the gateway user cannot change it,
  # and nothing is built or fetched for it.
  install -d -m 0755 -o root -g root "${lib}"
  rm -rf "${lib}/delfin-src.new"
  mkdir -p "${lib}/delfin-src.new"
  cp -R "${src}/delfin" "${lib}/delfin-src.new/delfin"
  find "${lib}/delfin-src.new" -name __pycache__ -type d -prune -exec rm -rf {} +
  chown -R root:root "${lib}/delfin-src.new"
  chmod -R u=rwX,go=rX "${lib}/delfin-src.new"
  rm -rf "${lib}/delfin-src"
  mv "${lib}/delfin-src.new" "${lib}/delfin-src"
  site="$("${venv}/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
  { [ -n "${site}" ] && [ -d "${site}" ]; } || die "no site-packages in ${venv}"
  printf '%s\n' "${lib}/delfin-src" > "${site}/delfin-messenger.pth"
  say "gateway venv ${venv}: websockets (hashed) + DELFIN from ${src}"
}

pair_phone() {
  # Prints the invitation and waits for the phone; refuses when paired.
  runuser -u "${GW_USER}" -- env DELFIN_MESSENGER_DIR="${STATE_DIR}/pairing" \
    DELFIN_MESSENGER_PORT="${API_PORT}" "${PREFIX}/bin/venv/bin/python" -m delfin.messenger pair \
    || say "not paired now (already paired? revoke first: docs/design/messenger-simplex.md)"
}

# --- local speech recognition (--whisper) -----------------------------------
install_whisper() {
  local base="${PREFIX}/bin/venv/bin/python" req tmpdir entry name want got
  [ -x "${base}" ] || die "--whisper needs the gateway venv (${PREFIX}/bin/venv) first"
  "${base}" -c 'import sys; sys.exit(sys.version_info[:2] != (3, 11))' \
    || die "--whisper: the pinned wheels are for CPython 3.11; the gateway venv is another version"
  req="${HERE}/whisper-requirements-${machine}.txt"
  [ -f "${req}" ] || die "no pinned whisper requirements for ${machine}"
  [ -x "${WHISPER_VENV}/bin/python" ] || "${base}" -m venv "${WHISPER_VENV}"
  "${WHISPER_VENV}/bin/python" -m pip install --quiet --require-hashes --no-deps -r "${req}" \
    || die "whisper packages did not install"
  tmpdir="$(mktemp -d)"
  for entry in ${WHISPER_FILES}; do
    name="${entry%%:*}"; want="${entry#*:}"
    if [ -f "${WHISPER_DIR}/${name}" ] \
       && [ "$(sha256sum "${WHISPER_DIR}/${name}" | cut -d' ' -f1)" = "${want}" ]; then
      cp "${WHISPER_DIR}/${name}" "${tmpdir}/${name}"
      continue
    fi
    curl -fL --retry 2 -o "${tmpdir}/${name}" \
      "https://huggingface.co/${WHISPER_REPO}/resolve/${WHISPER_REVISION}/${name}" \
      || { rm -rf "${tmpdir}"; die "download failed: ${name}"; }
    got="$(sha256sum "${tmpdir}/${name}" | cut -d' ' -f1)"
    [ "${got}" = "${want}" ] \
      || { rm -rf "${tmpdir}"; die "digest mismatch for ${name}: got ${got}, want ${want}"; }
  done
  install -d -m 0755 -o root -g root "${PREFIX}/models" "${WHISPER_DIR}"
  for entry in ${WHISPER_FILES}; do
    name="${entry%%:*}"
    install -m 0644 -o root -g root "${tmpdir}/${name}" "${WHISPER_DIR}/${name}"
  done
  rm -rf "${tmpdir}"
  say "speech recognition: ${WHISPER_DIR}, interpreter ${WHISPER_VENV}/bin/python"
}

# --- main -------------------------------------------------------------------
case "${action}" in
  uninstall) uninstall ;;
  install)
    ensure_user
    install -d -m 0755 -o root -g root "${PREFIX}" "${PREFIX}/bin"
    install -d -m 0700 -o "${GW_USER}" -g "${GW_USER}" "${STATE_DIR}"
    if [ -x "${PREFIX}/bin/simplex-chat" ]; then
      say "${PREFIX}/bin/simplex-chat already present -- keeping it (delete to re-download)"
    else
      fetch_binary "${PREFIX}/bin/simplex-chat"
    fi
    [ "${want_nftables}" = 1 ] && install_nftables
    install_gateway_venv
    [ "${want_whisper}" = 1 ] && install_whisper
    install_units
    if [ "${want_pair}" = 1 ]; then
      pair_phone
    else
      say "done. Pair the phone: sudo bash $0 --pair"
    fi
    ;;
esac
