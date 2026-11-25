"""The macOS launcher speaks to a server exactly as the Windows launcher does.

tools/macos-launcher/DELFIN.sh is the Mac's counterpart of
tools/windows-launcher/DELFIN.ps1. The server half is shared byte for byte,
so what matters on this side is that the Mac sends the same command, with
the same SSH options, accepts the same browser URL, and keeps the same rules
for profiles and invitations -- and that install and uninstall leave what
they promise. These tests run the script's own functions under bash; the
macOS-only parts (Gatekeeper's quarantine mark, the app icon, the Mac's
bash 3.2 and BSD tools, a real SSH round trip) run in
.github/workflows/macos-launcher.yml on macos-latest.
"""

from __future__ import annotations

import base64
import importlib.util
import os
import re
import shutil
import stat
import subprocess
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MAC = ROOT / "tools" / "macos-launcher"
SCRIPT = MAC / "DELFIN.sh"
WINDOWS_APP = (ROOT / "tools" / "windows-launcher" / "DELFIN.ps1").read_text(encoding="utf-8")
_KEY_HEADER = "-----BEGIN OPENSSH " + "PRIVATE KEY-----"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")


def _build():
    spec = importlib.util.spec_from_file_location("mac_build_package", MAC / "build_package.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _lib(snippet: str, env: dict | None = None) -> str:
    """Run SNIPPET with DELFIN.sh's functions loaded; return stdout."""
    program = f'DELFIN_LAUNCHER_LIBRARY=1 source "$1"\n{snippet}\n'
    result = subprocess.run(["bash", "-c", program, "delfin-test", str(SCRIPT)],
                            capture_output=True, text=True, timeout=60,
                            env={**os.environ, **(env or {})})
    assert result.returncode == 0, result.stderr
    return result.stdout


def _ps_check_script() -> str:
    found = re.findall(r"^\s*\$checkScript = '((?:[^']|'')*)'$", WINDOWS_APP, re.M)
    assert len(found) == 1
    return found[0].replace("''", "'")


def _quote(value: str) -> str:
    return "'" + value.replace("'", "'\\''") + "'"


def _ps_remote_command(folder, hb, hl, port, directory, executable, role, tunnel, keep, working, person, key):
    """Get-RemoteCommand of DELFIN.ps1, rebuilt from its own text."""
    body = re.search(r"^function Get-RemoteCommand\b.*?^\}", WINDOWS_APP, re.M | re.S).group(0)
    assert ("(Shell-Quote $launcher) + ' ' + (Shell-Quote $role) + ' ' + (Shell-Quote $tunnelPath)"
            in body), "the Windows argument order changed: change DELFIN.sh's remote_command with it"
    check = "bash -c " + _quote(_ps_check_script()) + f" delfin {folder} {hb} {hl}"
    return (check + " && bash -ilc " + _quote(". ~/" + folder + "/remote_bootstrap.sh") + " delfin "
            + f"{port} {_quote(directory)} {_quote(executable)} {_quote('~/' + folder + '/remote_launcher.py')} "
            + f"{_quote(role)} {_quote(tunnel)} {keep} {working} {_quote(person)} {_quote(key)}")


# -- the server command and the SSH options ----------------------------------

def test_the_check_script_is_the_windows_one():
    assert _lib('printf %s "$CHECK_SCRIPT"') == _ps_check_script()


@pytest.mark.parametrize("directory, executable, keep, working, key", [
    ("", "", 0, 0, ""),
    ("~/work it's", "/opt/delfin/bin/python", 1, 1, "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIAbc delfin-launcher-x"),
])
def test_the_server_command_is_the_windows_one(tmp_path, directory, executable, keep, working, key):
    for name in ("remote_bootstrap.sh", "remote_launcher.py"):
        shutil.copy(ROOT / "tools" / "windows-launcher" / name, tmp_path / name)
    tunnel = "/tmp/delfin-win-tunnel-0b8f1f7e-8b1c-4a4f-9f1e-3c5d7e9a1b2c/dashboard.sock"
    out = _lib(f'''HERE={_quote(str(tmp_path))}
P_PORT=8870 P_DIRECTORY={_quote(directory)} P_EXECUTABLE={_quote(executable)} P_KEEP={keep} P_WORKING={working}
script_version
printf '%s\\n' "$V_FOLDER" "$V_BOOTSTRAP" "$V_LAUNCHER" "$(person_name)"
remote_command dashboard {_quote(tunnel)} {_quote(key)}''')
    folder, hb, hl, person, command = out.split("\n", 4)
    import hashlib
    want_b = hashlib.sha256((tmp_path / "remote_bootstrap.sh").read_bytes()).hexdigest()
    want_l = hashlib.sha256((tmp_path / "remote_launcher.py").read_bytes()).hexdigest()
    assert (hb, hl) == (want_b, want_l)
    assert folder == f".cache/delfin-launcher/{want_b[:16]}{want_l[:16]}", "the server folder Windows uses"
    assert command == _ps_remote_command(folder, hb, hl, 8870, directory, executable, "dashboard",
                                         tunnel, keep, working, person, key)


def test_the_ssh_options_are_the_windows_ones():
    tunnel = "/tmp/delfin-win-tunnel-0b8f1f7e-8b1c-4a4f-9f1e-3c5d7e9a1b2c/dashboard.sock"
    dashboard = _lib(f"ssh_options Dashboard 40123 0 {tunnel}; printf '%s\\n' \"${{SSH_OPTS[@]}}\"").split()
    copy = _lib("ssh_options Copy; printf '%s\\n' \"${SSH_OPTS[@]}\"").split()
    windows = re.search(r"\$options = @\((.*?)\)\n", WINDOWS_APP).group(1)
    for option in re.findall(r"'([A-Za-z0-9]+=[A-Za-z0-9]+)'", windows):
        assert option in dashboard and option in copy, option
    assert dashboard[0] == "-tt" and "-tt" not in copy
    assert "ClearAllForwardings=no" in dashboard and "ClearAllForwardings=yes" in copy
    assert dashboard[-2:] == ["-L", f"127.0.0.1:40123:{tunnel}"], "the only forward, on loopback"
    assert "GatewayPorts=no" in dashboard and "ExitOnForwardFailure=yes" in dashboard
    for forbidden in ("ForwardAgent=yes", "ForwardX11=yes", "-A", "-X", "-Y", "-R", "-D"):
        assert forbidden not in dashboard + copy


def test_a_strange_tunnel_path_gets_no_forward():
    out = _lib("ssh_options Dashboard 40123 0 /tmp/elsewhere/dashboard.sock; printf '%s\\n' \"${SSH_OPTS[@]}\"").split()
    assert "-L" not in out and "ClearAllForwardings=yes" in out


def test_the_system_ssh_is_used_and_no_password_is_handled():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "\nSSH=/usr/bin/ssh\n" in text and "\nSCP=/usr/bin/scp\n" in text
    assert "\nPATH=/usr/bin:/bin:/usr/sbin:/sbin\n" in text
    assert "IdentitiesOnly=yes" in text and "StrictHostKeyChecking=ask" in text
    assert not re.search(r"sshpass|SSH_ASKPASS|security add-generic-password|-o\s+StrictHostKeyChecking=no", text)


# -- the browser URL -----------------------------------------------------------

TOKEN = "a" * 40


@pytest.mark.parametrize("url, ok", [
    (f"http://127.0.0.1:8866/?token={TOKEN}", True),
    (f"http://127.0.0.1:8866/voila/render/delfin_dashboard.ipynb?token={TOKEN}&session=abc-1", True),
    (f"http://localhost:8866/?token={TOKEN}", False),
    (f"http://192.0.2.1:8866/?token={TOKEN}", False),
    (f"https://127.0.0.1:8866/?token={TOKEN}", False),
    ("http://127.0.0.1:8866/?token=short", False),
    (f"http://127.0.0.1:8866/elsewhere?token={TOKEN}", False),
    (f"http://127.0.0.1:80/?token={TOKEN}", False),
    (f"http://u@127.0.0.1:8866/?token={TOKEN}", False),
    (f"http://127.0.0.1:8866/?token={TOKEN}#x", False),
    (f"http://127.0.0.1:8866/voila/render/../../x?token={TOKEN}", False),
    (f"http://127.0.0.1:8866/?token={TOKEN}\nhttp://evil", False),
])
def test_only_the_dashboard_url_is_opened(url, ok):
    encoded = base64.b64encode(url.encode()).decode()
    out = _lib(f"browser_url {encoded} && echo ACCEPTED; true")
    assert ("ACCEPTED" in out) is ok
    if ok:
        assert out.startswith(url)
        moved = _lib(f"with_port {_quote(url)} 40123")
        assert moved == url.replace(":8866/", ":40123/")


# -- profiles --------------------------------------------------------------------

@pytest.mark.parametrize("target, ok", [
    ("delfin@192.0.2.10", True), ("cluster-alias", True), ("a.b_c-d@host.example", True),
    ("-oProxyCommand=x", False), ("user@host -p 22", False), ("", False), ("user@host;id", False),
])
def test_a_target_is_an_alias_or_user_at_host(target, ok):
    assert ("YES" in _lib(f"valid_target {_quote(target)} && echo YES; true")) is ok


@pytest.mark.parametrize("path, want", [
    ("", ""), ("  '~/delfin'  ", "~/delfin"), ('"/opt/my env/bin/python"', "/opt/my env/bin/python"),
    ("relative/dir", "relative/dir"), ("C:\\Users\\x", None), ("$(id)", None), ("/x;rm", None),
])
def test_server_paths_are_paths(path, want):
    out = _lib(f"clean_path {_quote(path)} && echo; echo STATUS=$?")
    if want is None:
        assert "STATUS=1" in out
    else:
        assert out.startswith(want + "\n") and "STATUS=0" in out


def test_a_profile_is_saved_and_read_as_data(tmp_path):
    env = {"DELFIN_LAUNCHER_STORE": str(tmp_path / "store")}
    pid = "0b8f1f7e-8b1c-4a4f-9f1e-3c5d7e9a1b2c"
    _lib(f"save_profile {pid} 'Lab $(id)' delfin@192.0.2.10 '~/w' '' 0 1 0 1", env)
    file = tmp_path / "store" / "profiles" / f"{pid}.conf"
    assert stat.S_IMODE(file.stat().st_mode) == 0o600
    out = _lib(f"load_profile {pid}; printf '%s|' \"$P_NAME\" \"$P_TARGET\" \"$P_DIRECTORY\" \"$P_KEEP\" \"$P_REMEMBER\"", env)
    assert out == "Lab $(id)|delfin@192.0.2.10|~/w|1|1|"
    assert _lib("profile_ids", env).split() == [pid]


def test_the_person_name_follows_the_windows_rule():
    out = _lib("id() { echo 'Anna.Mueller-'; }; person_name")
    assert out == "anna-mueller"


# -- package, install, invitation, uninstall ----------------------------------------

def _package(tmp_path: Path, invite: dict | None = None) -> Path:
    build = _build()
    folder = tmp_path / "pkg" / "macos-launcher"
    folder.mkdir(parents=True)
    for name, data in build.contents().items():
        (folder / name).write_bytes(data)
    (folder / "SHA256SUMS.txt").write_text(build.checksums(build.contents()), encoding="ascii")
    if invite:
        import json
        (folder / "invite.json").write_text(json.dumps(invite, indent=2) + "\n", encoding="utf-8")
        (folder / "invite_key").write_text(_KEY_HEADER + "\nAAAA\n-----END OPENSSH PRIVATE KEY-----\n", encoding="ascii")
    return folder


def _run(folder: Path, home: Path, *args: str) -> subprocess.CompletedProcess:
    env = {**os.environ, "HOME": str(home),
           "DELFIN_LAUNCHER_STORE": str(home / "store"), "DELFIN_LAUNCHER_APPS": str(home / "Applications"),
           "DELFIN_LAUNCHER_SSH_DIR": str(home / ".ssh")}
    return subprocess.run(["bash", str(folder / "DELFIN.sh"), *args], capture_output=True, text=True,
                          env=env, timeout=120, stdin=subprocess.DEVNULL)


def test_a_changed_package_is_refused(tmp_path):
    folder = _package(tmp_path)
    assert _run(folder, tmp_path, "validate").returncode == 0
    (folder / "remote_launcher.py").write_text("print('changed')\n")
    result = _run(folder, tmp_path, "validate")
    assert result.returncode != 0 and "remote_launcher.py" in result.stderr


def test_install_with_an_invitation_then_uninstall_keeping_and_removing(tmp_path):
    host_key = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIHostKeyOfTheServer0000000000000000000"
    folder = _package(tmp_path, {"target": "delfin@192.0.2.10", "name": "anna", "expires": "202612010000Z",
                                 "host_keys": [host_key]})
    home = tmp_path / "home"
    home.mkdir()
    result = _run(folder, home, "install")
    assert result.returncode == 0, result.stderr + result.stdout
    assert "Your connection to delfin@192.0.2.10 is set up" in result.stdout
    app = home / "Applications" / "DELFIN.app" / "Contents"
    assert os.access(app / "MacOS" / "DELFIN", os.X_OK)
    assert "<key>CFBundleExecutable</key><string>DELFIN</string>" in (app / "Info.plist").read_text()
    for name in ("DELFIN.sh", "remote_bootstrap.sh", "remote_launcher.py", "setup_local.sh", "SHA256SUMS.txt"):
        assert (app / "Resources" / name).is_file(), name
    store = home / "store"
    assert stat.S_IMODE(store.stat().st_mode) == 0o700
    profiles = list((store / "profiles").glob("*.conf"))
    assert len(profiles) == 1
    text = profiles[0].read_text()
    assert "target=delfin@192.0.2.10\n" in text and "remember=1\n" in text and "keep=1\n" in text
    key = store / f"invite-{profiles[0].stem}.key"
    assert stat.S_IMODE(key.stat().st_mode) == 0o600, "ssh refuses a private key others can read"
    assert (home / ".ssh" / "known_hosts").read_text() == f"192.0.2.10 {host_key}\n"

    # A second install (an update) keeps the one connection and adds no second host line.
    assert _run(folder, home, "install").returncode == 0
    assert len(list((store / "profiles").glob("*.conf"))) == 1
    assert (home / ".ssh" / "known_hosts").read_text().count("\n") == 1

    (store / "local" / "env").mkdir(parents=True)
    result = _run(folder, home, "uninstall", "--keep-profiles", "--keep-local")
    assert result.returncode == 0, result.stderr
    assert not (home / "Applications" / "DELFIN.app").exists()
    assert profiles[0].is_file() and key.is_file() and (store / "local" / "env").is_dir()
    assert "Saved connections are kept" in result.stdout

    assert _run(folder, home, "install").returncode == 0
    result = _run(folder, home, "uninstall", "--remove-profiles", "--remove-local")
    assert result.returncode == 0, result.stderr
    assert not store.exists() and not (home / "Applications" / "DELFIN.app").exists()


@pytest.mark.parametrize("invite", [
    {"target": "root;id@host", "host_keys": ["ssh-ed25519 AAAA"]},
    {"target": "delfin@192.0.2.10", "host_keys": []},
])
def test_a_broken_invitation_is_refused(tmp_path, invite):
    folder = _package(tmp_path, invite)
    home = tmp_path / "home"
    home.mkdir()
    result = _run(folder, home, "install")
    assert result.returncode != 0
    assert "invitation" in result.stderr
    assert not (home / ".ssh" / "known_hosts").exists()


def test_the_script_runs_on_the_bash_macos_ships():
    """/bin/bash on a Mac is 3.2: these would fail there and pass here."""
    for path in (SCRIPT, MAC / "setup_local.sh", MAC / "Install.command", MAC / "Uninstall.command"):
        full = path.read_text(encoding="utf-8")
        text = "\n".join(line for line in full.splitlines() if not line.lstrip().startswith("#"))
        for construct in (r"declare -A", r"\$\{[A-Za-z_]+,,", r"\$\{[A-Za-z_]+\^\^", r"\bmapfile\b",
                          r"\breadarray\b", r"\|&", r"&>>", r"\bdeclare -n\b", r"\$\{[A-Za-z_]+@[QEPAa]\}",
                          r"^\s*set -u", r"\bcoproc\b"):
            assert not re.search(construct, text, re.M), f"{path.name}: {construct}"
        assert full.startswith("#!/bin/bash\n"), path.name
        assert b"\r" not in path.read_bytes()


def test_the_download_mark_is_offered_never_removed_silently():
    text = (MAC / "Install.command").read_text(encoding="utf-8")
    assert "xattr -d com.apple.quarantine" in text
    question = text.index("Remove the mark and continue? (y/N)")
    assert question < text.index("xattr -d com.apple.quarantine"), "asked before it is removed"
    assert "spctl" not in text and "--master-disable" not in text, "Gatekeeper itself is never touched"
    assert "cp -X" not in SCRIPT.read_text(encoding="utf-8"), "the install keeps the person's decision"
