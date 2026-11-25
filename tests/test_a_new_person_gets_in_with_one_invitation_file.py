"""A new person gets onto a hardened server with one file and no round trip.

Before: install, press Start, fail ("Permission denied"), send the key file
from the desktop to the owner, the owner runs `delfin-people --add-key`,
press Start again -- and then OpenSSH asks to compare a host fingerprint.

Now the owner runs `delfin-people --invite NAME` once. DELFIN-NAME.zip is the
launcher package plus the server's address and host keys and a key for this
invitation only, whose line in authorized_keys can do one thing: run
`--enroll`, which adds the new PC's own key (with the launcher's restriction)
and removes the invitation. It expires by itself.

The Windows half (Install.ps1 imports the invitation, DELFIN.ps1 enrolls on
the first Start) runs in Test-InstallRoundTrip.ps1 on Windows CI; here the
server half runs for real, and the Windows half is checked as text.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PEOPLE = ROOT / "scripts" / "pi" / "delfin_people.sh"
PKG = ROOT / "tools" / "windows-launcher"

#: The opening line of an OpenSSH private key, assembled so that no tracked
#: file carries it literally (tests/test_no_secrets_are_committed.py).
_KEY_HEADER = "-----BEGIN OPENSSH " + "PRIVATE KEY-----"


def _account(tmp_path: Path):
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    package = tmp_path / "delfin-windows-launcher.zip"
    with zipfile.ZipFile(package, "w") as archive:
        archive.writestr("windows-launcher/Install.cmd", "@echo off\r\n")
    host_keys = tmp_path / "etc-ssh"
    host_keys.mkdir()
    subprocess.run(["ssh-keygen", "-q", "-t", "ed25519", "-N", "", "-f", str(host_keys / "ssh_host_ed25519_key")], check=True)
    env = dict(os.environ, HOME=str(home), DELFIN_LAUNCHER_ZIP=str(package),
               DELFIN_HOST_KEYS=str(host_keys / "ssh_host_*_key.pub"))
    return home, env, host_keys


def _people(env, *args, stdin=""):
    return subprocess.run(["bash", str(PEOPLE), *args], env=env, input=stdin,
                          capture_output=True, text=True, timeout=60)


def _pc_key(tmp_path: Path) -> str:
    key = tmp_path / "pc_key"
    subprocess.run(["ssh-keygen", "-q", "-t", "ed25519", "-N", "", "-C", "delfin-launcher-anna", "-f", str(key)], check=True)
    return (tmp_path / "pc_key.pub").read_text().strip()


def test_an_invitation_is_one_file_with_server_host_keys_and_a_key_of_its_own(tmp_path):
    home, env, host_keys = _account(tmp_path)

    result = _people(env, "--invite", "anna", "--host", "192.0.2.10")

    assert result.returncode == 0, result.stderr
    out = home / "invites" / "DELFIN-anna.zip"
    with zipfile.ZipFile(out) as archive:
        names = archive.namelist()
        invite = json.loads(archive.read("windows-launcher/invite.json"))
        key = archive.read("windows-launcher/invite_key").decode()
    assert "windows-launcher/Install.cmd" in names, "the whole launcher package travels with it"
    assert invite["target"].endswith("@192.0.2.10") and invite["name"] == "anna"
    pinned = " ".join((host_keys / "ssh_host_ed25519_key.pub").read_text().split()[:2])
    assert invite["host_keys"] == [pinned]
    assert key.startswith(_KEY_HEADER)
    assert oct(out.stat().st_mode & 0o777) == "0o600", "an invitation is like a password"


def test_the_invitation_key_can_only_enroll_and_expires(tmp_path):
    home, env, _ = _account(tmp_path)
    _people(env, "--invite", "anna", "--host", "192.0.2.10", "--hours", "5")

    line = (home / ".ssh" / "authorized_keys").read_text().strip()
    options = line.split(" ssh-ed25519 ")[0]
    assert options.startswith("restrict,"), "no pty, no forwarding, no agent"
    assert re.search(r'command="[^"]*delfin_people\.sh --enroll anna"', options)
    assert re.search(r'expiry-time="\d{12}Z"', options)
    assert line.endswith(" delfin-invite-anna")


def test_enrolling_adds_the_pcs_key_restricted_and_uses_up_the_invitation(tmp_path):
    home, env, _ = _account(tmp_path)
    _people(env, "--invite", "anna", "--host", "192.0.2.10")
    pc = _pc_key(tmp_path)

    result = _people(env, "--enroll", "anna", stdin=pc + "\n")

    assert result.returncode == 0, result.stderr
    assert "DELFIN_ENROLLED" in result.stdout
    lines = (home / ".ssh" / "authorized_keys").read_text().splitlines()
    assert not [l for l in lines if "delfin-invite-anna" in l], "used once"
    added = [l for l in lines if l.endswith(" delfin-launcher-anna")]
    assert added == ["restrict,port-forwarding,pty " + " ".join(pc.split()[:2]) + " delfin-launcher-anna"]
    assert oct((home / ".ssh" / "authorized_keys").stat().st_mode & 0o777) == "0o600"


def test_enrolling_reads_one_line_and_never_a_second(tmp_path):
    home, env, _ = _account(tmp_path)
    _people(env, "--invite", "anna", "--host", "192.0.2.10")
    pc = _pc_key(tmp_path)

    result = _people(env, "--enroll", "anna", stdin=pc + "\nssh-ed25519 AAAAsecond attacker\n")

    assert result.returncode == 0, result.stderr
    text = (home / ".ssh" / "authorized_keys").read_text()
    assert "AAAAsecond" not in text and "attacker" not in text


@pytest.mark.parametrize("bad", ["", "ssh-rsa AAAAB3NzaC1yc2E x",
                                 'command="sh" ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAA x'])
def test_enrolling_refuses_anything_but_one_plain_ed25519_key(tmp_path, bad):
    home, env, _ = _account(tmp_path)
    _people(env, "--invite", "anna", "--host", "192.0.2.10")
    before = (home / ".ssh" / "authorized_keys").read_text()

    result = _people(env, "--enroll", "anna", stdin=bad + "\n")

    assert result.returncode != 0
    assert (home / ".ssh" / "authorized_keys").read_text() == before, "nothing added, invitation kept"


def test_a_name_that_is_not_a_plain_name_is_refused(tmp_path):
    home, env, _ = _account(tmp_path)
    for name in ("-rf", "Anna", "a b", "../x"):
        assert _people(env, "--invite", name, "--host", "192.0.2.10").returncode == 2
    assert not (home / ".ssh" / "authorized_keys").exists()


def test_the_windows_side_imports_and_uses_the_invitation():
    install = (PKG / "Install.ps1").read_text(encoding="utf-8")
    app = (PKG / "DELFIN.ps1").read_text(encoding="utf-8")
    # Install: first saved connection, pinned host keys, the key kept per connection.
    assert "Join-Path $PSScriptRoot 'invite.json'" in install
    assert "$profiles = @([pscustomobject]@{ Id=$profileId; Name=$target; Target=$target" in install
    assert "Join-Path $SshFolder 'known_hosts'" in install
    assert "('invite-' + $profileId + '.key')" in install
    # The invitation files are not copied into the app folder as they are.
    copied = re.search(r"foreach \(\$file in @\(([^)]*)\)\)", install).group(1)
    assert "invite" not in copied
    # First Start: the invitation key, offered alone, without a prompt; deleted when used.
    assert "('invite-' + $p.Id + '.key')" in app
    assert "@('-o','BatchMode=yes','-o','IdentitiesOnly=yes','-i',$inviteKeyFile,[string]$p.Target)" in app
    assert "$enrollProcess.StandardInput.WriteLine($loginKey.Public)" in app
    i_enrolled = app.index("if ($enrollOut.Result -match 'DELFIN_ENROLLED')")
    assert app.index("Remove-Item -LiteralPath $inviteKeyFile -Force", i_enrolled) - i_enrolled < 200
