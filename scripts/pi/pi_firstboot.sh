#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# First-boot provisioning for a freshly flashed Pi SD card / SSD.
#
# Lives on the boot partition as /boot/firmware/delfin/ next to:
#   pi_bootstrap.sh  pi_install.sh  delfin-src.tgz  VERSION
# and is started once by cloud-init (runcmd) as root. Needs network (LAN).
# Result: same state as `deploy_pi.sh --bootstrap` followed by a deploy.
# Log: /var/log/delfin-firstboot.log   Marker: /var/lib/delfin-firstboot.done
# ============================================================================

HERE="$(cd "$(dirname "$0")" && pwd)"
LOGIN_USER="${LOGIN_USER:-pi}"
AGENT_USER="${AGENT_USER:-delfin}"
MARKER=/var/lib/delfin-firstboot.done

exec >> /var/log/delfin-firstboot.log 2>&1
echo "=== delfin firstboot $(date -Is) ==="
[[ -f "$MARKER" ]] && { echo "already done"; exit 0; }

for _ in $(seq 60); do
    curl -fsS --max-time 5 https://astral.sh >/dev/null 2>&1 && break
    echo "waiting for network..."; sleep 5
done

SUDO_USER="$LOGIN_USER" AGENT_USER="$AGENT_USER" bash "$HERE/pi_bootstrap.sh"

agent_home="$(getent passwd "$AGENT_USER" | cut -d: -f6)"
install -o "$AGENT_USER" -g "$AGENT_USER" -m 755 "$HERE/pi_install.sh" "$agent_home/delfin/pi_install.sh"
sudo -iu "$AGENT_USER" bash "$agent_home/delfin/pi_install.sh" install "$(cat "$HERE/VERSION")" \
    < "$HERE/delfin-src.tgz"

# Chromium's system libraries are apt packages: root installs them here, once.
if "$agent_home/delfin/.venv/bin/python" -m playwright install-deps chromium </dev/null; then
    touch /var/lib/delfin-playwright-deps.done
fi

# Host lockdown: sshd key only, firewall, fail2ban, security updates.
LOGIN_USER="$LOGIN_USER" AGENT_USER="$AGENT_USER" \
    bash "$agent_home/delfin/src/scripts/pi/pi_harden_host.sh"

touch "$MARKER"
echo "=== delfin firstboot done $(date -Is) ==="

# The Raspberry Pi OS kernel again with Landlock (about an hour), then one
# tryboot into it; a kernel that does not come up falls back to the stock one.
# LANDLOCK_KERNEL=0 leaves the stock kernel alone.
if [[ "${LANDLOCK_KERNEL:-1}" == "1" ]]; then
    kernel_script="$agent_home/delfin/src/scripts/pi/pi_kernel_landlock.sh"
    if bash "$kernel_script" build; then
        bash "$kernel_script" try
    else
        echo "warning: Landlock kernel build failed; the stock kernel stays (retry: sudo bash $kernel_script build)"
    fi
fi
