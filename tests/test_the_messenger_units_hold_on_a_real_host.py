"""What the provisioned units must do on a real host, beyond starting.

Input: the provisioner, driven with the stand-ins of
test_the_messenger_provisioning_covers_every_pinned_server. Measured on a
Raspberry Pi after the first live setup:

- The bridge ran with the group ``delfin-gw``. The agent service reads a
  window's /proc/<pid>/environ to tell it from an agent's child, and Linux
  lets only a process with the same uid *and* gid read it; unreadable means
  "may watch, not act", so every phone command was refused (-32003). The
  bridge now runs with the service user's own group; only its socket
  directory carries ``delfin-gw`` (setgid, so the socket inherits it).
- The nftables rule that keeps every user but ``delfin-gw`` off the
  daemon's API port was loaded once with ``nft -f -`` and nowhere stored:
  a reboot or a reload of the host ruleset (``flush ruleset``) removed it,
  and the agent's user could then speak to the daemon directly. The rule
  is now a file the daemon unit loads before every start, and the daemon
  does not start without it.
- ``--lang=en`` sets the phone texts to English.
"""
from __future__ import annotations

from pathlib import Path

from test_the_messenger_provisioning_covers_every_pinned_server import run_script

HARDEN = Path(__file__).resolve().parents[1] / "scripts" / "pi" / "pi_harden_host.sh"


def _run(tmp_path, args=None):
    done = run_script(tmp_path, args=args,
                      env_extra={"MESSENGER_ETC_DIR": str(tmp_path / "etc")})
    assert done.returncode == 0, done.stderr
    return {p.name: p.read_text() for p in (tmp_path / "units").glob("*.service")}


def _lines(unit):
    return [line.strip() for line in unit.splitlines()]


def test_the_bridge_runs_with_the_service_users_own_group(tmp_path):
    bridge = _lines(_run(tmp_path)["delfin-messenger-bridge.service"])
    assert "Group=delfin-gw" not in bridge
    assert not any(line.startswith("Group=") and line != "Group=" for line in bridge)
    assert "RuntimeDirectoryMode=2750" in bridge
    pre = [line for line in bridge if line.startswith("ExecStartPre=+")]
    assert any("chgrp delfin-gw /run/delfin-messenger-bridge" in line for line in pre)
    assert any("chmod 2750 /run/delfin-messenger-bridge" in line for line in pre)


def test_the_port_rule_is_stored_and_loaded_before_every_daemon_start(tmp_path):
    daemon = _lines(_run(tmp_path)["simplex-chat-daemon.service"])
    rule = tmp_path / "etc" / "messenger.nft"
    text = rule.read_text()
    assert 'meta skuid "delfin-gw" accept' in text and "dport 5227 drop" in text
    assert any(line.startswith("ExecStartPre=+") and "-f " + str(rule) in line
               for line in daemon)
    assert any(line.startswith("After=") and "nftables.service" in line for line in daemon)
    assert "PartOf=nftables.service" in daemon


def test_without_nftables_the_daemon_loads_no_rule(tmp_path):
    daemon = _lines(_run(tmp_path, ["--no-nftables"])["simplex-chat-daemon.service"])
    assert not any(line.startswith("ExecStartPre=") for line in daemon)
    assert not (tmp_path / "etc" / "messenger.nft").exists()


def test_the_host_firewall_keeps_the_messenger_rule_when_it_reloads():
    text = HARDEN.read_text()
    assert "/etc/delfin-messenger/messenger.nft" in text


def test_english_phone_texts_are_one_flag_away(tmp_path):
    gateway = _run(tmp_path, ["--lang=en"])["delfin-messenger-gateway.service"]
    assert "Environment=DELFIN_MESSENGER_LANG=en" in gateway
