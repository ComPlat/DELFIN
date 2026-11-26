"""The macOS launcher ZIP holds what its sources say, and an invitation carries it.

tools/delfin-macos-launcher.zip is committed and published as it is. It
takes the server half from the Windows launcher and the micromamba fetch
from delfin/installers at build time; a change to any of them without a
rebuild would ship the old copy, and the Mac would then use a different
~/.cache/delfin-launcher/<version>/ folder on the server than a PC -- or
check its download against a checksum the sources no longer match.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import stat
import subprocess
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MAC = ROOT / "tools" / "macos-launcher"
ZIP = ROOT / "tools" / "delfin-macos-launcher.zip"
PEOPLE = ROOT / "scripts" / "pi" / "delfin_people.sh"
_KEY_HEADER = "-----BEGIN OPENSSH " + "PRIVATE KEY-----"


def _build():
    spec = importlib.util.spec_from_file_location("mac_build_package", MAC / "build_package.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


#: What signing adds (sign_and_notarize.sh): allowed in a signed ZIP, not listed.
SIGNATURE = ("macos-launcher/DELFIN.app/Contents/_CodeSignature/", "macos-launcher/DELFIN.app/Contents/CodeResources")


def test_the_zip_holds_the_current_sources():
    build = _build()
    files = build.contents()
    with zipfile.ZipFile(ZIP) as archive:
        names = [n for n in archive.namelist() if not n.endswith("/") and not n.startswith(SIGNATURE)]
        assert sorted(names) == sorted("macos-launcher/" + n for n in build.NAMES + ["SHA256SUMS.txt"])
        for name in build.NAMES:
            # Digests, not bytes: a byte diff of a mismatching icon or ZIP
            # entry is what pytest then renders, for minutes.
            got = hashlib.sha256(archive.read("macos-launcher/" + name)).hexdigest()
            assert got == hashlib.sha256(files[name]).hexdigest(), (
                f"{name} changed since the ZIP was built: python tools/macos-launcher/build_package.py "
                "(from an LF checkout)")
        assert archive.read("macos-launcher/SHA256SUMS.txt").decode() == build.checksums(files)
    assert (MAC / "SHA256SUMS.txt").read_text(encoding="ascii") == build.checksums(files)


def test_the_scripts_stay_executable_after_extraction():
    build = _build()
    with zipfile.ZipFile(ZIP) as archive:
        for name in build.NAMES:
            info = archive.getinfo("macos-launcher/" + name)
            mode = info.external_attr >> 16
            assert info.create_system == 3 and stat.S_ISREG(mode), name
            assert bool(mode & 0o111) == (name in build.EXECUTABLE), name


def test_the_package_is_a_finished_app_bundle():
    build = _build()
    files = build.contents()
    plist = files["DELFIN.app/Contents/Info.plist"].decode()
    assert f"<string>{build.version()}</string>" in plist and "@VERSION@" not in plist
    assert "<key>CFBundleExecutable</key><string>DELFIN</string>" in plist
    assert files[build.RESOURCES + "DELFIN.icns"][:4] == b"icns"
    windows = ROOT / "tools" / "windows-launcher"
    assert build.SOURCES[build.RESOURCES + "remote_bootstrap.sh"] == windows / "remote_bootstrap.sh"
    assert build.SOURCES[build.RESOURCES + "remote_launcher.py"] == windows / "remote_launcher.py"


def test_the_one_line_install_is_pinned_and_checked():
    """One command in Terminal: a fixed release asset, its SHA-256 checked
    before anything of it runs, never a script piped into a shell."""
    import hashlib
    build = _build()
    digest = hashlib.sha256(ZIP.read_bytes()).hexdigest()
    readme = (ROOT / "README.md").read_text(encoding="utf-8").replace("\r\n", "\n")
    assert build.readme_block(digest) in readme, (
        "README.md's one-line install does not name this ZIP: run build_package.py "
        "(or, for a signed ZIP, put sign_and_notarize.sh's digest there)")
    line = build.install_line(digest)
    assert f"/releases/download/v{build.version()}/delfin-macos-launcher.zip" in line
    assert line.index("shasum -a 256 -c") < line.index("ditto -x -k") < line.index("Install.command")
    assert "--proto =https" in line
    import re
    assert not re.search(r"\|\s*(sudo\s+)?(sh|bash|zsh)(\s|$)", line), "nothing is piped into a shell"


#: A fixed server host key: the test needs a key line, not a key pair.
_HOST_KEY = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIGJ1aWx0LWZvci10aGlzLXRlc3QtMDEyMzQ1Njc4OQ"


def _stand_ins(folder: Path) -> Path:
    """ssh-keygen and date for delfin-people, the same on every machine:
    the invitation's key pair is written as fixed text, and the expiry is a
    fixed time (delfin-people asks GNU date for it; a Mac has BSD date)."""
    folder.mkdir()
    keygen = folder / "ssh-keygen"
    keygen.write_text(
        "#!/bin/bash\n"
        "while [ $# -gt 0 ]; do case $1 in -f) out=$2; shift 2 ;; -C) comment=$2; shift 2 ;; *) shift ;; esac; done\n"
        f"printf '%s\\nc3RhbmQtaW4=\\n-----END OPENSSH PRIVATE KEY-----\\n' '{_KEY_HEADER}' > \"$out\"\n"
        "printf 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIHN0YW5kLWluLWZvci10aGUtaW52aXRhdGlvbjAw %s\\n'"
        " \"$comment\" > \"$out.pub\"\n",
        encoding="ascii")
    date = folder / "date"
    date.write_text("#!/bin/bash\necho 202612010000Z\n", encoding="ascii")
    for tool in (keygen, date):
        tool.chmod(0o755)
    return folder


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
    (host_keys / "ssh_host_ed25519_key.pub").write_text(_HOST_KEY + " root@server\n", encoding="ascii")
    stand_ins = _stand_ins(tmp_path / "stand-ins")
    env = dict(os.environ, HOME=str(home), DELFIN_LAUNCHER_ZIP=str(windows),
               DELFIN_HOST_KEYS=str(host_keys / "ssh_host_*_key.pub"),
               PATH=f"{stand_ins}{os.pathsep}{os.environ['PATH']}")
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
        ["bash", str(folder / "DELFIN.app" / "Contents" / "Resources" / "DELFIN.sh"), "install"], capture_output=True, text=True, timeout=120,
        stdin=subprocess.DEVNULL,
        env={**os.environ, "HOME": str(mac_home), "DELFIN_LAUNCHER_STORE": str(mac_home / "store"),
             "DELFIN_LAUNCHER_APPS": str(mac_home / "Applications"),
             "DELFIN_LAUNCHER_SSH_DIR": str(mac_home / ".ssh")})
    assert install.returncode == 0, install.stderr + install.stdout
    user = subprocess.run(["id", "-un"], capture_output=True, text=True).stdout.strip()
    assert f"Your connection to {user}@192.0.2.10 is set up" in install.stdout
    assert (mac_home / ".ssh" / "known_hosts").read_text() == f"192.0.2.10 {_HOST_KEY}\n"
    assert len(list((mac_home / "store").glob("invite-*.key"))) == 1


def test_without_a_mac_package_the_invitation_is_as_before(tmp_path):
    """The existing test of the invitation covers this case in full; here
    only that the message says how to add the Mac package."""
    text = PEOPLE.read_text(encoding="utf-8")
    assert 'mac_package="${DELFIN_MAC_LAUNCHER_ZIP:-}"' in text
    assert "[[ -f \"$mac_package\" ]] || mac_package=\"\"" in text
