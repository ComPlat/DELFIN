"""One run of scripts/pi/pi_messenger.sh leaves a working messenger.

Input: the provisioner, driven with the stand-ins of
test_the_messenger_provisioning_covers_every_pinned_server (no network, no
root). Before, three steps were left to the person and one platform was
refused outright:

- Raspberry Pi OS (Debian) was refused: only Ubuntu was recognised. Any
  Linux with glibc >= 2.34 now takes the Ubuntu 22.04 build (it needs
  glibc 2.34, OpenSSL 3, libgmp, zlib; measured with ldd/objdump on the
  pinned v7.0.3 binary); Ubuntu 24.04 keeps its own build. Older glibc
  is refused before any download.
- The daemon asked for a display name on stdin and stopped there; it now
  creates its bot profile itself (``--create-bot-display-name``; a second
  start keeps the existing profile, both measured with the pinned binary).
- The gateway venv had to be made and filled by hand; the run now makes it
  from Python 3.11, installs ``messenger-requirements.txt`` with
  ``--require-hashes --no-deps`` and a root-owned copy of the deployed
  ``delfin`` package (no build step, nothing unhashed from PyPI).
- The bridge could not write the service user's home (ProtectHome=
  read-only), so files from the phone and new sessions failed; that home
  is now its one writable path.
- The agent service is a unit of its own, so the bridge has a service to
  talk to without a terminal session; phone-started sessions get their own
  folder root, which ``--no-phone-sessions`` leaves out.
"""
from __future__ import annotations

import re

from test_the_messenger_provisioning_covers_every_pinned_server import (
    SCRIPT, run_script)

DEBIAN_12 = 'ID=debian\nVERSION_ID="12"\n'
D22_AARCH64 = "91f332edc813cab78a848da7ba75ed91bcb8d09abf1fc92eca74330dc29ec970"


def _units(tmp_path):
    return {p.name: p.read_text() for p in (tmp_path / "units").glob("*.service")}


def test_raspberry_pi_os_takes_the_glibc_2_34_build(tmp_path):
    done = run_script(tmp_path, machine="aarch64", os_release=DEBIAN_12,
                      glibc="2.36", sha256=D22_AARCH64)
    assert done.returncode == 0, done.stderr
    calls = (tmp_path / "calls.log").read_text()
    assert "simplex-chat-ubuntu-22_04-aarch64" in calls


def test_an_older_glibc_is_refused_before_any_download(tmp_path):
    done = run_script(tmp_path, os_release=DEBIAN_12, glibc="2.31")
    assert done.returncode != 0 and "no pinned" in done.stderr
    assert "curl" not in (tmp_path / "calls.log").read_text()


def test_the_daemon_creates_its_own_bot_profile(tmp_path):
    assert run_script(tmp_path).returncode == 0
    daemon = _units(tmp_path)["simplex-chat-daemon.service"]
    assert "--create-bot-display-name DELFIN" in daemon
    assert "--create-bot-allow-files" in daemon


def test_the_gateway_venv_is_made_and_filled_from_hashes(tmp_path):
    done = run_script(tmp_path)
    assert done.returncode == 0, done.stderr
    calls = (tmp_path / "calls.log").read_text()
    assert re.search(r"-m venv \S+/prefix/bin/venv", calls)
    assert re.search(r"--require-hashes --no-deps -r \S+/messenger-requirements.txt", calls)
    package = tmp_path / "prefix/lib/delfin-src/delfin"
    assert (package / "messenger/gateway.py").is_file()
    assert not list(package.rglob("__pycache__"))
    pth = (tmp_path / "site/delfin-messenger.pth").read_text().strip()
    assert pth == str(tmp_path / "prefix/lib/delfin-src")
    units = _units(tmp_path)
    assert "delfin-messenger-gateway.service" in units
    assert "delfin-messenger-bridge.service" in units


def test_messenger_requirements_pin_websockets_with_hashes():
    text = (SCRIPT.parent / "messenger-requirements.txt").read_text()
    assert re.search(r"^websockets==\d", text, re.M)
    assert "--hash=sha256:" in text


def test_the_bridge_may_write_the_service_home_and_starts_sessions_there(tmp_path):
    assert run_script(tmp_path).returncode == 0
    home = tmp_path / "home-delfin"
    bridge = _units(tmp_path)["delfin-messenger-bridge.service"]
    assert "ProtectHome=read-only" in bridge
    assert f"ReadWritePaths={home}" in bridge
    assert f"Environment=DELFIN_MESSENGER_WORKSPACES={home}/DELFIN-Handy" in bridge
    assert (home / "DELFIN-Handy").is_dir()
    assert "delfin-agent-service.service" in bridge


def test_no_phone_sessions_leaves_the_folder_root_out(tmp_path):
    assert run_script(tmp_path, args=["--no-phone-sessions"]).returncode == 0
    bridge = _units(tmp_path)["delfin-messenger-bridge.service"]
    assert "DELFIN_MESSENGER_WORKSPACES" not in bridge


def test_the_agent_service_runs_as_its_own_unit(tmp_path):
    home = tmp_path / "home-delfin"
    venv = home / "delfin/.venv/bin/python"
    venv.parent.mkdir(parents=True)
    venv.write_text("#!/bin/sh\nexit 0\n")
    venv.chmod(0o755)
    assert run_script(tmp_path).returncode == 0
    units = _units(tmp_path)
    service = units["delfin-agent-service.service"]
    assert "User=delfin" in service
    assert f"{venv} -P -m delfin.agent.service.supervisor" in service
    assert f"{home}/delfin/agent.env" in service
    # The bridge talks with the same DELFIN the service runs.
    assert f"ExecStart={venv} -m delfin.agent.service.bridge" in units[
        "delfin-messenger-bridge.service"]


def test_uninstall_removes_the_new_unit_too():
    text = SCRIPT.read_text()
    assert "rm -f \"${UNIT_DIR}/delfin-agent-service.service\"" in text
