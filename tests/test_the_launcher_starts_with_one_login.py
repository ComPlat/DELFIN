"""A start of the Windows launcher costs one SSH login, not two.

The launcher used to open a short SSH connection to check its two server
scripts, and then a second one for the dashboard. Each connection is a
password -- and on a cluster with one-time codes, a code: two codes for
every start. The check now runs on the dashboard connection itself, in
front of the login shell; only a missing or changed copy ends it with
DELFIN_SCRIPTS_MISSING and exit 10, after which the app copies the files
and connects once more.

The command is run here as a server runs it: the account's login shell gets
it as one string. A stand-in `bash` on PATH records the login shell start
instead of loading this machine's startup files.
"""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1] / "tools" / "windows-launcher"
APP = (ROOT / "DELFIN.ps1").read_text(encoding="utf-8")


def _check_script() -> str:
    found = re.findall(r"^\s*\$checkScript = '((?:[^']|'')*)'$", APP, re.M)
    assert len(found) == 1
    return found[0].replace("''", "'")


def _quote(value: str) -> str:
    """Shell-Quote of DELFIN.ps1."""
    return "'" + value.replace("'", "'\\''") + "'"


def _remote_command(folder: str, hb: str, hl: str) -> str:
    """What Get-RemoteCommand returns for a dashboard start (port 0, no
    directory, no location, Keep off, no working shell, no person, no key)."""
    body = re.search(r"^function Get-RemoteCommand\b.*?^\}", APP, re.M | re.S).group(0)
    assert "return $check + ' && ' + 'bash -ilc ' + (Shell-Quote $bootstrap)" in body
    check = "bash -c " + _quote(_check_script()) + f" delfin {folder} {hb} {hl}"
    bootstrap = ". ~/" + folder + "/remote_bootstrap.sh"
    launcher = "~/" + folder + "/remote_launcher.py"
    return (check + " && bash -ilc " + _quote(bootstrap) + " delfin 0 '' '' "
            + _quote(launcher) + " 'dashboard' '' 0 0 '' ''")


@pytest.fixture
def server(tmp_path):
    home = tmp_path / "account home"
    home.mkdir()
    stubs = tmp_path / "stubs"
    stubs.mkdir()
    real_bash = shutil.which("bash")
    record = tmp_path / "login-shell.txt"
    # `bash -c CHECK` is the real bash; `bash -ilc ...` is the login shell,
    # recorded with its arguments and the umask it was given.
    (stubs / "bash").write_text(
        "#!/bin/sh\n"
        f'if [ "$1" = -ilc ]; then umask > {_quote(str(record))}; '
        f'printf "%s\\n" "$@" >> {_quote(str(record))}; exit 0; fi\n'
        f'exec {_quote(real_bash)} "$@"\n')
    (stubs / "bash").chmod(0o755)
    hb = hashlib.sha256((ROOT / "remote_bootstrap.sh").read_bytes()).hexdigest()
    hl = hashlib.sha256((ROOT / "remote_launcher.py").read_bytes()).hexdigest()
    folder = ".cache/delfin-launcher/" + hb[:16] + hl[:16]
    env = {"PATH": f"{stubs}:/usr/bin:/bin", "HOME": str(home)}

    def start():
        # The account's login shell runs the remote command as one string;
        # the dashboard keeps the account's umask.
        return subprocess.run(["/bin/sh", "-c", "umask 022; " + _remote_command(folder, hb, hl)],
                              env=env, cwd=str(home), capture_output=True, text=True, timeout=30)

    return home, folder, record, start


def test_a_present_copy_starts_the_login_shell_on_the_same_connection(server):
    home, folder, record, start = server
    copy = home / folder
    copy.mkdir(parents=True)
    for name in ("remote_bootstrap.sh", "remote_launcher.py"):
        shutil.copyfile(ROOT / name, copy / name)

    result = start()

    assert result.returncode == 0, result.stderr
    assert "DELFIN_SCRIPTS_MISSING" not in result.stdout
    lines = record.read_text().splitlines()
    assert lines[0] == "0022", "the check's umask 077 must not reach the dashboard"
    assert lines[1:5] == ["-ilc", ". ~/" + folder + "/remote_bootstrap.sh", "delfin", "0"]
    assert "~/" + folder + "/remote_launcher.py" in lines
    assert "dashboard" in lines


def test_a_missing_copy_ends_the_connection_with_the_marker_and_no_login_shell(server):
    home, folder, record, start = server

    result = start()

    assert result.returncode == 10
    assert "DELFIN_SCRIPTS_MISSING" in result.stdout
    assert not record.exists(), "no dashboard may start on scripts that are not there"
    assert stat.S_IMODE((home / folder).stat().st_mode) == 0o700


def test_a_changed_copy_is_not_started(server):
    home, folder, record, start = server
    copy = home / folder
    copy.mkdir(parents=True)
    shutil.copyfile(ROOT / "remote_bootstrap.sh", copy / "remote_bootstrap.sh")
    (copy / "remote_launcher.py").write_text("print('someone else')\n")

    result = start()

    assert result.returncode == 10
    assert not record.exists()


def test_the_app_copies_and_connects_again_only_on_the_marker():
    dashboard = APP[APP.index("Get-RemoteCommand $p 'dashboard'"):]
    assert "if ($line -eq 'DELFIN_SCRIPTS_MISSING') { $scriptsMissing = $true; continue }" in dashboard
    retry = dashboard[dashboard.index("if ($scriptsMissing -and"):]
    retry = retry[:retry.index("break")]
    assert "-not $opened" in retry and "-not $scriptsCopied" in retry
    assert "$sshProcess.ExitCode -eq 10" in retry
    assert "Copy-ServerScripts ([string]$p.Target) $identity $version" in retry
    assert "$scriptsCopied = $true" in retry and "continue" in retry
    # No separate check connection is left before the dashboard.
    assert "Invoke-OpenSSH $ssh" not in APP
    assert "Sync-ServerScripts" not in APP


def test_the_command_carries_no_double_quotes_for_windows_argument_passing(server):
    # Quote-NativeArgument wraps each argument in double quotes for ssh.exe;
    # the remote command is kept free of them (Test-Launcher.ps1 checks the
    # PowerShell side, this the composed text).
    home, folder, record, start = server
    hb = hl = "0" * 64
    assert '"' not in _remote_command(folder, hb, hl)
