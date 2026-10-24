#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Lock the Pi down at the host level (the agent's own cage is separate:
# Landlock, bubblewrap, the restricted agent user).
#
# Runs ON the Pi as root. Idempotent; first boot calls it, run it again any
# time. What it does:
#   1. unused network services off (rpcbind / NFS)
#   2. sshd: key only for every user, no root login, no X11 forwarding,
#      3 tries, only the login and the agent user may log in
#   3. nftables: inbound default deny; SSH only from private IPv4 nets and
#      the IPv6 prefixes in SSH_V6 (default: link-local only); mDNS, DHCP and
#      ICMP allowed; outbound open (agent egress is DELFIN's proxy's job)
#   4. fail2ban for sshd (journal backend, nftables bans)
#   5. unattended-upgrades for Debian (security) updates; Raspberry Pi
#      kernel/firmware stay manual, so a custom kernel is never replaced
#   6. sudo asks for the login user's password, once that user has one
#
# Over SSH the firewall is armed with a 120 s rollback: a successful NEW
# SSH login within that time keeps it (no marker file, no root needed).
#
# Environment: LOGIN_USER (pi), AGENT_USER (delfin),
#   SSH_V6="2001:db8::/32 ..." extra IPv6 prefixes allowed to reach SSH.
# ============================================================================

LOGIN_USER="${LOGIN_USER:-pi}"
AGENT_USER="${AGENT_USER:-delfin}"
SSH_V6="${SSH_V6:-}"
ROLLBACK_UNIT=delfin-fw-rollback

[[ $EUID -eq 0 ]] || { echo "error: run with sudo" >&2; exit 1; }

echo "==> [1/6] unused network services"
for unit in rpcbind.socket rpcbind.service nfs-blkmap.service; do
    systemctl disable --now "$unit" >/dev/null 2>&1 || true
done
systemctl mask rpcbind.socket rpcbind.service >/dev/null 2>&1 || true

echo "==> [2/6] sshd"
# Passwords go off below. A login user without a key would be locked out of
# SSH by that, so stop here rather than after the reload.
login_home="$(getent passwd "$LOGIN_USER" | cut -d: -f6)"
if ! grep -qE '^(ssh|ecdsa|sk-)' "$login_home/.ssh/authorized_keys" 2>/dev/null; then
    echo "error: $LOGIN_USER has no key in $login_home/.ssh/authorized_keys;" >&2
    echo "       SSH passwords would be turned off and lock $LOGIN_USER out." >&2
    echo "       Add a key first (ssh-copy-id $LOGIN_USER@<pi>), then run this again." >&2
    exit 1
fi
# Loaded before cloud-init's 50-* file; sshd keeps the first value it reads.
cat > /etc/ssh/sshd_config.d/40-delfin-hardening.conf <<EOF
# DELFIN agent box (pi_harden_host.sh): key only, no root, few tries.
PasswordAuthentication no
KbdInteractiveAuthentication no
PermitRootLogin no
X11Forwarding no
MaxAuthTries 3
LoginGraceTime 30
AllowUsers $LOGIN_USER $AGENT_USER
EOF
chmod 644 /etc/ssh/sshd_config.d/40-delfin-hardening.conf
# A password exception from earlier manual setups would undo the above.
rm -f /etc/ssh/sshd_config.d/60-delfin-password.conf
sshd -t
systemctl reload ssh

echo "==> [3/6] firewall (nftables)"
v6="fe80::/10"
for p in $SSH_V6; do v6="$v6, $p"; done
rules="$(mktemp)"
cat > "$rules" <<EOF
#!/usr/sbin/nft -f
# DELFIN agent box (pi_harden_host.sh): inbound default deny, outbound open.
flush ruleset
table inet filter {
    set ssh_v4 { type ipv4_addr; flags interval; elements = { 10.0.0.0/8, 172.16.0.0/12, 192.168.0.0/16 } }
    set ssh_v6 { type ipv6_addr; flags interval; elements = { $v6 } }
    chain input {
        type filter hook input priority filter; policy drop;
        ct state established,related accept
        ct state invalid drop
        iif "lo" accept
        meta l4proto icmp accept
        meta l4proto ipv6-icmp accept
        udp dport 68 accept comment "DHCPv4 client"
        udp dport 546 accept comment "DHCPv6 client"
        udp dport 5353 accept comment "mDNS"
        tcp dport 22 ip saddr @ssh_v4 accept
        tcp dport 22 ip6 saddr @ssh_v6 accept
    }
    chain forward {
        type filter hook forward priority filter; policy drop;
    }
    chain output {
        type filter hook output priority filter; policy accept;
    }
}
EOF
nft -c -f "$rules"
if [[ -n "${SSH_CONNECTION:-}" ]]; then
    # The proof that SSH still gets through is a login that succeeds after
    # the rules are in place; without one in 120 s the firewall goes again.
    # No marker file and no root needed to confirm: just log in anew.
    armed="$(date +%s)"
    systemctl stop "$ROLLBACK_UNIT.timer" 2>/dev/null || true
    systemd-run --unit "$ROLLBACK_UNIT" --on-active=120 /bin/sh -c \
        "journalctl -q --since @$armed -t sshd -t sshd-session --grep '^Accepted ' -n 1 | grep -q . \
         || /usr/sbin/nft flush ruleset" >/dev/null
    echo "    rollback armed: open a NEW SSH session within 120 s to keep the firewall;"
    echo "    without a successful new login it is removed again"
fi
nft -f "$rules"
install -m 755 "$rules" /etc/nftables.conf
rm -f "$rules"
systemctl enable nftables >/dev/null 2>&1

echo "==> [4/6] fail2ban"
export DEBIAN_FRONTEND=noninteractive
apt-get install -y -q fail2ban python3-systemd unattended-upgrades apt-listchanges >/dev/null
cat > /etc/fail2ban/jail.d/delfin-sshd.local <<'EOF'
# DELFIN agent box (pi_harden_host.sh)
[sshd]
enabled  = true
backend  = systemd
banaction = nftables-multiport
maxretry = 5
findtime = 10m
bantime  = 1h
bantime.increment = true
EOF
systemctl enable fail2ban >/dev/null 2>&1
systemctl restart fail2ban

echo "==> [5/6] unattended security upgrades"
cat > /etc/apt/apt.conf.d/20auto-upgrades <<'EOF'
APT::Periodic::Update-Package-Lists "1";
APT::Periodic::Unattended-Upgrade "1";
APT::Periodic::AutocleanInterval "7";
EOF
systemctl enable --now unattended-upgrades >/dev/null 2>&1

echo "==> [6/6] sudo password for $LOGIN_USER"
sudo_rule=/etc/sudoers.d/95-delfin-login-password
if [[ "$(passwd -S "$LOGIN_USER" | awk '{print $2}')" == "P" ]]; then
    # Later sudoers files win: this overrides cloud-init's NOPASSWD rule.
    printf '%s ALL=(ALL) ALL\n' "$LOGIN_USER" > "$sudo_rule.tmp"
    chmod 0440 "$sudo_rule.tmp"
    visudo -cf "$sudo_rule.tmp" >/dev/null
    mv "$sudo_rule.tmp" "$sudo_rule"
    echo "    sudo now asks for $LOGIN_USER's password"
else
    echo "    $LOGIN_USER has no password yet, sudo stays passwordless."
    echo "    Set one with 'sudo passwd $LOGIN_USER', then run this script again."
fi

echo
echo "Host hardening done."
