"""The macOS launcher ZIP holds what its sources say, and an invitation carries it.

tools/delfin-macos-launcher.zip is committed and published as it is. It
takes the server half from the Windows launcher and the micromamba fetch
from delfin/installers at build time; a change to any of them without a
rebuild would ship the old copy, and the Mac would then use a different
~/.cache/delfin-launcher/<version>/ folder on the server than a PC -- or
check its download against a checksum the sources no longer match.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import stat
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MAC = ROOT / "tools" / "macos-launcher"
ZIP = ROOT / "tools" / "delfin-macos-launcher.zip"
PEOPLE = ROOT / "scripts" / "pi" / "delfin_people.sh"


def _build():
    spec = importlib.util.spec_from_file_location("mac_build_package", MAC / "build_package.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_zip_holds_the_current_sources():
    build = _build()
    files = build.contents()
    with zipfile.ZipFile(ZIP) as archive:
        names = [n for n in archive.namelist() if not n.endswith("/")]
        assert sorted(names) == sorted("macos-launcher/" + n for n in build.NAMES + ["SHA256SUMS.txt"])
        for name in build.NAMES:
            assert archive.read("macos-launcher/" + name) == files[name], (
                f"{name} changed since the ZIP was built: python tools/macos-launcher/build_package.py "
                "(from an LF checkout)")
        assert archive.read("macos-launcher/SHA256SUMS.txt").decode() == build.checksums(files)
    assert (MAC / "SHA256SUMS.txt").read_text(encoding="ascii") == build.checksums(files)


def test_the_scripts_stay_executable_after_extraction():
    build = _build()
    with zipfile.ZipFile(ZIP) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            mode = info.external_attr >> 16
            assert info.create_system == 3 and stat.S_ISREG(mode), info.filename
            name = info.filename.split("/", 1)[1]
            assert bool(mode & 0o111) == (name in build.EXECUTABLE), info.filename


def test_the_server_half_is_the_windows_one_and_the_app_lists_the_package():
    build = _build()
    assert build.SHARED["remote_bootstrap.sh"] == ROOT / "tools" / "windows-launcher" / "remote_bootstrap.sh"
    assert build.SHARED["remote_launcher.py"] == ROOT / "tools" / "windows-launcher" / "remote_launcher.py"
    script = (MAC / "DELFIN.sh").read_text(encoding="utf-8")
    listed = script.split('PACKAGE_FILES="', 1)[1].split('"', 1)[0].split()
    assert sorted(listed) == sorted(build.NAMES), "the installer copies every file of the package"


@pytest.mark.skipif(shutil.which("ssh-keygen") is None or shutil.which("bash") is None
                    or sys.platform == "darwin",
                    reason="delfin-people runs on the Linux server (GNU date); needs ssh-keygen and bash")
def test_an_invitation_carries_the_mac_package_too(tmp_path):
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    packages = tmp_path / "tools"
    packages.mkdir()
    windows = packages / "delfin-windows-launcher.zip"
    with zipfile.ZipFile(windows, "w") as archive:
        archive.writestr("windows-launcher/Install.cmd", "@echo off\r\n")
    shutil.copy(ZIP, packages / "delfin-macos-launcher.zip")
    host_keys = tmp_path / "etc-ssh"
    host_keys.mkdir()
    subprocess.run(["ssh-keygen", "-q", "-t", "ed25519", "-N", "", "-f", str(host_keys / "ssh_host_ed25519_key")],
                   check=True)
    env = dict(os.environ, HOME=str(home), DELFIN_LAUNCHER_ZIP=str(windows),
               DELFIN_HOST_KEYS=str(host_keys / "ssh_host_*_key.pub"))
    result = subprocess.run(["bash", str(PEOPLE), "--invite", "anna", "--host", "192.0.2.10"], env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert "On a Mac: macos-launcher/Install.command" in result.stdout
    with zipfile.ZipFile(home / "invites" / "DELFIN-anna.zip") as archive:
        names = archive.namelist()
        assert "windows-launcher/Install.cmd" in names and "macos-launcher/Install.command" in names
        mac_invite = json.loads(archive.read("macos-launcher/invite.json"))
        assert mac_invite == json.loads(archive.read("windows-launcher/invite.json"))
        assert archive.read("macos-launcher/invite_key") == archive.read("windows-launcher/invite_key")
        mode = archive.getinfo("macos-launcher/Install.command").external_attr >> 16
        assert mode & 0o111, "the invitation keeps the Mac scripts executable"
        archive.extractall(tmp_path / "unpacked")

    # Step 2 on the Mac, run here: the installer imports the invitation.
    folder = tmp_path / "unpacked" / "macos-launcher"
    mac_home = tmp_path / "mac"
    mac_home.mkdir()
    install = subprocess.run(
        ["bash", str(folder / "DELFIN.sh"), "install"], capture_output=True, text=True, timeout=120,
        stdin=subprocess.DEVNULL,
        env={**os.environ, "HOME": str(mac_home), "DELFIN_LAUNCHER_STORE": str(mac_home / "store"),
             "DELFIN_LAUNCHER_APPS": str(mac_home / "Applications"),
             "DELFIN_LAUNCHER_SSH_DIR": str(mac_home / ".ssh")})
    assert install.returncode == 0, install.stderr + install.stdout
    user = subprocess.run(["id", "-un"], capture_output=True, text=True).stdout.strip()
    assert f"Your connection to {user}@192.0.2.10 is set up" in install.stdout
    host_key = " ".join((host_keys / "ssh_host_ed25519_key.pub").read_text().split()[:2])
    assert (mac_home / ".ssh" / "known_hosts").read_text() == f"192.0.2.10 {host_key}\n"
    assert len(list((mac_home / "store").glob("invite-*.key"))) == 1


def test_without_a_mac_package_the_invitation_is_as_before(tmp_path):
    """The existing test of the invitation covers this case in full; here
    only that the message says how to add the Mac package."""
    text = PEOPLE.read_text(encoding="utf-8")
    assert 'mac_package="${DELFIN_MAC_LAUNCHER_ZIP:-}"' in text
    assert "[[ -f \"$mac_package\" ]] || mac_package=\"\"" in text
