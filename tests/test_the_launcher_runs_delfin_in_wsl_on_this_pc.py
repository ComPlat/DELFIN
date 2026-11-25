"""The Windows launcher starts DELFIN in WSL on the same PC, not only on a server.

DELFIN does not run on Windows itself. Someone without a Linux server runs it
in a WSL distribution on their own PC. The launcher reaches it with wsl.exe
instead of ssh.exe: wsl.exe hands its command to the distribution's login
shell exactly as sshd does, so the same check-and-bootstrap command runs
unchanged. What differs:

* no tunnel: WSL forwards the dashboard's 127.0.0.1 port to Windows'
  127.0.0.1, so the server half starts no relay;
* the connection is still watched: closing the window (wsl.exe ends, the
  Linux side gets SIGHUP) stops a dashboard that is not kept, and a closing
  marker named for this connection tells Windows the cleanup is done;
* the two scripts are copied over stdin (no scp).

The DELFIN.ps1 side is checked by Test-Launcher.ps1 on Windows; here the
server half and the shared contract between the two.
"""

from __future__ import annotations

import os
import re
import shlex
import stat
import subprocess
import uuid
from types import SimpleNamespace

import pytest

from tests.test_windows_launcher import _ROOT, _bootstrap_env, _launcher_file, launcher

APP = (_ROOT / 'DELFIN.ps1').read_text(encoding='utf-8')
LOCAL_ID = 'delfin-win-local-' + str(uuid.UUID(int=0x1234))


@pytest.mark.parametrize('keep', [False, True])
def test_closing_the_wsl_window_stops_the_dashboard_unless_kept(monkeypatch, capsys, keep):
    from delfin.agent import where
    monkeypatch.setattr(where, 'dashboard', lambda: {})
    monkeypatch.setattr(launcher.shutil, 'which', lambda name: '/usr/bin/tmux')
    monkeypatch.setattr(launcher, 'select_port', lambda preferred: preferred)
    handlers = {}
    monkeypatch.setattr(launcher.signal, 'signal', lambda signum, handler: handlers.__setitem__(signum, handler))
    commands = []

    def run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(launcher.subprocess, 'run', run)

    def window_closed(*args, **kwargs):
        # What wsl.exe's end delivers: SIGHUP, turned into KeyboardInterrupt.
        handlers[launcher.signal.SIGHUP](launcher.signal.SIGHUP, None)
    monkeypatch.setattr(launcher, 'run_dashboard_connection', window_closed)
    monkeypatch.setattr(launcher.sys, 'argv', ['helper', '51234', '/home/me', '/env/bin/python', '',
                                               str(int(keep)), '0', LOCAL_ID])
    with pytest.raises(KeyboardInterrupt):
        launcher.main()
    assert launcher.signal.SIGHUP in handlers and launcher.signal.SIGTERM in handlers
    created = next(cmd for cmd in commands if cmd[1] == 'new-session')
    # The dashboard listens on 127.0.0.1 only, at the port Windows found free.
    assert '--ip 127.0.0.1 --strict-port --port 51234' in created[-1]
    assert any(cmd[1] == 'kill-session' for cmd in commands) is (not keep)
    assert 'DELFIN_CONNECTION_CLOSED:' + LOCAL_ID in capsys.readouterr().out


def test_a_wsl_connection_attaches_the_dashboard_and_stays_to_clean_up(monkeypatch, tmp_path, capsys):
    import base64
    from contextlib import nullcontext
    from delfin.agent import where
    from delfin.dashboard import session

    private_dir = tmp_path / 'private'
    private_dir.mkdir(mode=0o700)
    record_path = private_dir / 'dashboard.json'
    record_path.write_text('{}')
    record_path.chmod(0o600)
    token = 'b' * 43
    created = {}

    def dashboard():
        if not created:
            return {}
        return {'running': True, 'host': launcher.socket.gethostname(), 'tmux': created['name'],
                'port': 51234, 'token': token, 'url': '/voila/render/dashboard.ipynb'}
    monkeypatch.setattr(where, 'dashboard', dashboard)
    monkeypatch.setattr(where, 'record_path', lambda: record_path)
    monkeypatch.setattr(session, 'live_records', lambda: [])
    monkeypatch.setattr(launcher.shutil, 'which', lambda name: '/usr/bin/tmux')
    monkeypatch.setattr(launcher, 'select_port', lambda preferred: preferred)
    monkeypatch.setattr(launcher.socket, 'create_connection', lambda *a, **kw: nullcontext())
    monkeypatch.setattr(launcher.signal, 'signal', lambda *a: None)
    relays = []
    monkeypatch.setattr(launcher, 'DashboardRelay', lambda *a: relays.append(a))
    commands = []

    def run(command, **kwargs):
        commands.append(command)
        if command[1] == 'new-session':
            created['name'] = command[command.index('-s') + 1]
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(launcher.subprocess, 'run', run)

    def replaced(*args):
        raise AssertionError('a WSL connection must not exec tmux: it would skip the cleanup')
    monkeypatch.setattr(launcher.os, 'execv', replaced)
    monkeypatch.setattr(launcher.sys, 'argv', ['helper', '51234', '/home/me', '/env/bin/python', '', '0', '0', LOCAL_ID])
    launcher.main()
    assert relays == [], 'WSL forwards the port itself; no relay socket'
    verbs = [command[1] for command in commands]
    assert verbs == ['new-session', 'attach-session', 'kill-session']
    lines = capsys.readouterr().out.splitlines()
    url = base64.b64decode(lines[0].removeprefix('DELFIN_BROWSER:')).decode()
    assert url == f'http://127.0.0.1:51234/voila/render/dashboard.ipynb?token={token}'
    assert lines[-1] == 'DELFIN_CONNECTION_CLOSED:' + LOCAL_ID


@pytest.mark.parametrize('argv_tail', [
    ['', '0', '0', 'delfin-win-local-xyz'],
    ['', '0', '0', 'delfin-win-tunnel-' + str(uuid.UUID(int=1))],
    ['', '0', '0', LOCAL_ID + '; rm -rf ~'],
    ['/tmp/delfin-win-tunnel-' + str(uuid.UUID(int=1)) + '/dashboard.sock', '0', '0', LOCAL_ID],
])
def test_the_server_half_refuses_a_malformed_local_connection(monkeypatch, argv_tail):
    from delfin.agent import where
    monkeypatch.setattr(where, 'dashboard', lambda: {})
    monkeypatch.setattr(launcher.subprocess, 'run', lambda *a, **kw: pytest.fail('nothing may start'))
    monkeypatch.setattr(launcher.sys, 'argv', ['helper', '0', '/home/me', '/env/bin/python'] + argv_tail)
    with pytest.raises(RuntimeError, match='Invalid local connection'):
        launcher.main()


def _run_bootstrap_local(tmp_path, local_id):
    home = tmp_path / 'home'
    home.mkdir(exist_ok=True)
    env, python = _bootstrap_env(tmp_path, home)
    payload = _launcher_file(tmp_path, 'import sys; print("ARGS=" + "|".join(sys.argv[1:]))')
    result = subprocess.run(['bash', str(_ROOT / 'remote_bootstrap.sh'), '0', '', '', payload, 'dashboard',
                             '', '0', '0', '', '', local_id],
                            env=env, text=True, capture_output=True, timeout=30, cwd=str(home))
    return result, home, python


def test_the_bootstrap_hands_the_local_connection_to_the_launcher(tmp_path):
    result, home, python = _run_bootstrap_local(tmp_path, LOCAL_ID)
    assert result.returncode == 0, result.stderr
    assert f'ARGS=0|{home}|{python}||0|0|{LOCAL_ID}' in result.stdout


def test_an_ssh_start_passes_no_local_connection(tmp_path):
    result, home, python = _run_bootstrap_local(tmp_path, '')
    assert result.returncode == 0, result.stderr
    assert f'ARGS=0|{home}|{python}||0|0\n' in result.stdout


@pytest.mark.parametrize('bad', ['delfin-win-local-' + 'g' * 36, 'delfin-win-local-1234', 'x' + LOCAL_ID,
                                 LOCAL_ID + 'a', 'delfin-win-local-' + str(uuid.UUID(int=0xabcdef)).upper()])
def test_the_bootstrap_refuses_a_malformed_local_connection(tmp_path, bad):
    result, _, _ = _run_bootstrap_local(tmp_path, bad)
    assert result.returncode == 1
    assert 'Invalid local connection' in result.stderr
    assert 'ARGS=' not in result.stdout


def test_app_bootstrap_and_launcher_agree_on_the_local_connection():
    # DELFIN.ps1 makes the id from a GUID and appends it after the login key
    # (bootstrap $11); all three sides accept exactly that shape.
    body = re.search(r'^function Get-RemoteCommand\b.*?^\}', APP, re.M | re.S).group(0)
    assert "(Shell-Quote $publicKey) + $local" in body
    assert "$local = if ($localId) { ' ' + (Shell-Quote $localId) } else { '' }" in body
    assert "$localId = 'delfin-win-local-' + [guid]::NewGuid().ToString()" in APP
    pattern = re.search(r"\$localId -notmatch '\^(delfin-win-local-[^$']+)\$'", body).group(1)
    assert re.fullmatch(pattern, LOCAL_ID)
    assert launcher.LOCAL_ID.pattern == pattern
    server = (_ROOT / 'remote_bootstrap.sh').read_text(encoding='utf-8')
    assert 'local_id=${11:-}' in server
    assert '"$working" ${local_id:+"$local_id"}' in server


def _wsl_copy_script():
    """The sh script Copy-WslScripts runs in WSL, as sh receives it."""
    body = re.search(r'^function Copy-WslScripts\b.*?^\}', APP, re.M | re.S).group(0)
    found = re.findall(r"^\s*\$copyScript = '((?:[^']|'')*)'\r?$", body, re.M)
    assert len(found) == 1
    assert ("Get-WslArguments $distribution ('sh -c ' + (Shell-Quote $copyScript) + ' delfin ' + $version.Folder"
            " + ' ' + $name)") in body
    return found[0].replace("''", "'")


def test_the_scripts_reach_wsl_over_stdin_into_an_owner_only_folder(tmp_path):
    # What Copy-WslScripts does after the check reported exit 10: one
    # `wsl.exe --exec sh -c SCRIPT delfin FOLDER NAME` per file, the file on stdin.
    home = tmp_path / 'linux home'
    home.mkdir()
    folder = '.cache/delfin-launcher/' + 'c' * 32
    for name in ('remote_bootstrap.sh', 'remote_launcher.py'):
        subprocess.run(['sh', '-c', _wsl_copy_script(), 'delfin', folder, name],
                       input=(_ROOT / name).read_bytes(), env={'PATH': '/usr/bin:/bin', 'HOME': str(home)},
                       check=True, timeout=15, cwd=str(tmp_path))
        copy = home / folder / name
        assert copy.read_bytes() == (_ROOT / name).read_bytes()
        assert stat.S_IMODE(copy.stat().st_mode) == 0o600
    assert stat.S_IMODE((home / folder).stat().st_mode) == 0o700


def test_the_wsl_command_runs_through_the_login_shell_as_over_ssh(tmp_path):
    # wsl.exe -d X --cd ~ -- COMMAND hands COMMAND, verbatim from its raw
    # command line, to $SHELL -c, as sshd does: the remote command DELFIN.ps1
    # builds must survive one shell reading, with the local connection as
    # the last argument.
    body = re.search(r'^function Get-WslArguments\b.*?^\}', APP, re.M | re.S).group(0)
    assert "return '-d ' + $distribution + ' --cd ~ -- ' + $command" in body
    home = tmp_path / 'home'
    stubs = tmp_path / 'stubs'
    stubs.mkdir()
    home.mkdir()
    record = tmp_path / 'login-shell.txt'
    real_bash = subprocess.run(['sh', '-c', 'command -v bash'], capture_output=True, text=True).stdout.strip()
    (stubs / 'bash').write_text(
        '#!/bin/sh\n'
        f'if [ "$1" = -ilc ]; then printf "%s\\n" "$@" > {shlex.quote(str(record))}; exit 0; fi\n'
        f'exec {shlex.quote(real_bash)} "$@"\n')
    (stubs / 'bash').chmod(0o755)
    import hashlib
    hb = hashlib.sha256((_ROOT / 'remote_bootstrap.sh').read_bytes()).hexdigest()
    hl = hashlib.sha256((_ROOT / 'remote_launcher.py').read_bytes()).hexdigest()
    folder = '.cache/delfin-launcher/' + hb[:16] + hl[:16]
    (home / folder).mkdir(parents=True)
    for name in ('remote_bootstrap.sh', 'remote_launcher.py'):
        (home / folder / name).write_bytes((_ROOT / name).read_bytes())
    check = re.findall(r"^\s*\$checkScript = '((?:[^']|'')*)'$", APP, re.M)[0].replace("''", "'")
    q = shlex.quote
    command = (f'bash -c {q(check)} delfin {folder} {hb} {hl} && bash -ilc '
               + q('. ~/' + folder + '/remote_bootstrap.sh')
               + f" delfin 51234 '' '' {q('~/' + folder + '/remote_launcher.py')} 'dashboard' '' 0 0 'me' '' "
               + q(LOCAL_ID))
    result = subprocess.run(['/bin/sh', '-c', command], env={'PATH': f'{stubs}:/usr/bin:/bin', 'HOME': str(home)},
                            cwd=str(home), capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    args = record.read_text().splitlines()
    assert args[:4] == ['-ilc', '. ~/' + folder + '/remote_bootstrap.sh', 'delfin', '51234']
    assert args[-1] == LOCAL_ID


def test_the_token_goes_only_to_a_listener_the_launcher_trusts():
    # SSH: the connection's own ssh.exe. WSL: WSL's relay from its install
    # folder, or no Windows listener at all (mirrored networking).
    body = re.search(r'^function Test-TrustedListener\b.*?^\}', APP, re.M | re.S).group(0)
    assert "if ($transport -ne 'WSL') { return ($owners -contains $connectionProcess) }" in body
    assert "$process.ProcessName -ne 'wslrelay'" in body
    assert "(Join-Path $env:ProgramFiles 'WSL')" in body
    job = APP[APP.index('$browserJob = Start-Job'):]
    job = job[:job.index('Start-Process $url')]
    assert 'Test-TrustedListener $transport' in job
    assert "Invoke-WebRequest -Uri \"http://127.0.0.1:$local/login\"" in job


def test_wsl_profiles_take_no_ssh_options_login_key_or_scp():
    worker = APP[APP.index("$isWsl = (Get-Transport $p) -eq 'WSL'"):]
    assert "$options = if ($isWsl) { @() } else { @(Get-SSHOptions" in worker
    assert "if ($p.RememberLogin -and -not $isWsl)" in worker
    assert "if ($isWsl) { Copy-WslScripts ([string]$p.Target) $version } else { Copy-ServerScripts" in worker
    # No Windows quoting around the WSL command (wsl.exe would hand the
    # quotes to the shell); SSH keeps it.
    assert "$info.Arguments = Get-WslArguments ([string]$p.Target) $remote" in worker
