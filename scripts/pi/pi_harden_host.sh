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
#      3 tries, only the login and the agent user may log in (and, with
#      one account per person, the persons that are not locked)
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
# One account per person (delfin_accounts.sh): every member of this group
# whose account is not locked (expiry reached) may log in too. Without the
# group (the default setup) the list is the two accounts above, unchanged.
people_allowed() {
    local u expire today out=""
    today=$(( $(date +%s) / 86400 ))
    for u in $(getent group delfin-people | cut -d: -f4 | tr ',' ' '); do
        expire="$(getent shadow "$u" | cut -d: -f8)"
        [[ -n "$expire" && "$expire" -le "$today" ]] && continue
        out+=" $u"
    done
    printf '%s' "$out"
}
PEOPLE_ALLOW="$(people_allowed)"
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
#
# `Match all` first, and it is load-bearing. sshd reads the drop-ins in
# sorted order and CONCATENATES them: a `Match` in a file that sorts
# earlier (30-temp-password.conf, which a person-onboarding setup writes
# to let a new PC register its key with the account password) stays open
# across the file boundary, and every directive below would then apply to
# that one user instead of globally. Measured on a live box on
# 2026-11-18: the whole block below, AllowUsers included, was scoped to
# the agent account, so every other account on the machine kept the
# OpenSSH defaults -- six auth tries, a two-minute grace, and no
# allow-list at all. `Match all` closes any such block before these lines.
cat > /etc/ssh/sshd_config.d/40-delfin-hardening.conf <<EOF
# DELFIN agent box (pi_harden_host.sh): key only, no root, few tries.
Match all
PasswordAuthentication no
KbdInteractiveAuthentication no
PermitRootLogin no
X11Forwarding no
MaxAuthTries 3
LoginGraceTime 30
AllowUsers $LOGIN_USER $AGENT_USER$PEOPLE_ALLOW
EOF
chmod 644 /etc/ssh/sshd_config.d/40-delfin-hardening.conf
# A password exception from earlier manual setups would undo the above.
rm -f /etc/ssh/sshd_config.d/60-delfin-password.conf
# An exception that sorts BEFORE this file is not removed: it is how a new
# person's PC registers its key (the launcher sends the key on the first
# password login). It is named instead, because an exception nobody can
# see is an exception nobody withdraws.
for earlier in /etc/ssh/sshd_config.d/*.conf; do
    case "$(basename "$earlier")" in
        [0-3]*) : ;;
        *) continue ;;
    esac
    [ -e "$earlier" ] || continue
    if grep -qsE '^[[:space:]]*Match[[:space:]]' "$earlier"; then
        echo "note: $earlier opens a Match block and keeps a password" >&2
        echo "      exception for the accounts it names. The hardening" >&2
        echo "      above is global (Match all) and unaffected." >&2
        echo "      Remove that file once every PC has registered its key:" >&2
        echo "      sudo rm $earlier && sudo sshd -t && sudo systemctl reload ssh" >&2
    fi
done
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
# The messenger's port rule (pi_messenger.sh), kept across this flush; a
# glob that matches nothing is no error.
include "/etc/delfin-messenger/*.nft"
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
         || { /usr/sbin/nft flush ruleset; for f in /etc/delfin-messenger/*.nft; do [ -f \"\$f\" ] && /usr/sbin/nft -f \"\$f\"; done; }" >/dev/null
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
#
# maxretry counts SSH authentication failures, and one click in the Windows
# launcher makes up to THREE connections -- the script check, the copy, and
# the dashboard -- each asking for the password once. So a single attempt
# with a wrong password costs three failures, and at maxretry = 5 a person
# was banned on their SECOND click.
#
# What that looked like from the outside: the password was refused, then
# refused again, and then nothing worked at all, with the same message. A
# ban is indistinguishable from a wrong password to whoever is typing, so
# the person retries, which extends it. Measured on the lab Pi on
# 2026-11-19, where it cost an afternoon and 21 failed attempts before the
# jail was looked at.
#
# 10 allows three launcher attempts. 5m is enough to slow a brute-force
# attempt to nothing -- it needs millions of guesses, and bantime.increment
# lengthens every further ban of the same address -- while being a pause a
# person waits out rather than a wall they report as broken.
[sshd]
enabled  = true
backend  = systemd
banaction = nftables-multiport
maxretry = 10
findtime = 10m
bantime  = 5m
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
