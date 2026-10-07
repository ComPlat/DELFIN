"""Checks for the portable launcher's security-sensitive return routing."""
import importlib.util
from pathlib import Path
import socket
from urllib.parse import parse_qs, urlsplit

import pytest

from delfin import cli_voila


spec = importlib.util.spec_from_file_location(
    'windows_remote_launcher',
    Path(__file__).resolve().parents[1] / 'tools/windows-launcher/remote_launcher.py',
)
launcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)


def test_return_uses_only_latest_kernel_of_authenticated_server():
    token = 'a' * 43
    record = {'port': 8866, 'token': token, 'url': '/voila/render/dashboard.ipynb',
              'resume_path': '/voila/render/return.ipynb'}
    host = socket.gethostname().split('.')[0]
    def session(name, stamp, port=8866, secret=token, node=host):
        return {'session_name': name, 'updated_at': stamp, 'host': node,
                'request_url': f'http://localhost:{port}/?token={secret}'}
    rows = [session('old', 1), session('new & kept', 2),
            session('other-port', 100, port=9000),
            session('other-server-token', 100, secret='b' * 43),
            session('other-host', 100, node='another-node'),
            session('bad-heartbeat', 'broken')]
    url = urlsplit(launcher.return_url(record, rows))
    assert url.hostname == '127.0.0.1'
    assert url.port == 8866
    assert url.path == '/voila/render/return.ipynb'
    assert parse_qs(url.query) == {'token': [token], 'session': ['new & kept']}


@pytest.mark.parametrize('resume_path', ['', '//attacker.invalid/', 'https://attacker.invalid/'])
def test_without_valid_resume_route_opens_local_authenticated_dashboard(resume_path):
    record = {'port': 8866, 'token': 'a' * 43, 'resume_path': resume_path}
    url = urlsplit(launcher.return_url(record, []))
    assert url.netloc == '127.0.0.1:8866'
    assert url.path == '/'
    assert parse_qs(url.query) == {'token': ['a' * 43]}


def test_strict_port_refuses_to_retarget_tunnel(monkeypatch, capsys):
    # Abort the preflight credential maintenance; this test never touches credentials.
    from delfin.agent import process_guard
    def no_credential_work(*args):
        raise RuntimeError('isolated test')
    monkeypatch.setattr(process_guard, 'protect', no_credential_work)
    monkeypatch.setattr(cli_voila, '_voila_is_available', lambda: True)
    monkeypatch.setattr(cli_voila, '_select_port', lambda port: port + 1)
    monkeypatch.setenv('DELFIN_LAUNCH_CWD', '/tmp')
    with pytest.raises(SystemExit) as exited:
        cli_voila.main(['--no-browser', '--strict-port', '--port', '8866'])
    assert exited.value.code == 1
    assert 'refusing to change the SSH tunnel target' in capsys.readouterr().err


@pytest.mark.parametrize("server_port", [8866, 9000])
@pytest.mark.parametrize("tmux_name", ["delfin-win-8866", "delfin", ""])
def test_reconnect_attaches_existing_server_without_starting_another(monkeypatch, tmp_path, capsys, server_port, tmux_name):
    import base64
    from contextlib import nullcontext
    from types import SimpleNamespace
    from delfin.agent import where
    from delfin.dashboard import session

    private_dir = tmp_path / 'private'
    private_dir.mkdir(mode=0o700)
    record_path = private_dir / 'dashboard.json'
    record_path.write_text('{}')
    record_path.chmod(0o600)
    token = 'a' * 43
    record = {'running': True, 'host': socket.gethostname(), 'tmux': tmux_name,
              'port': server_port, 'token': token, 'url': '/voila/render/dashboard.ipynb',
              'resume_path': '/voila/render/return.ipynb'}
    monkeypatch.setattr(where, 'dashboard', lambda: record)
    monkeypatch.setattr(where, 'record_path', lambda: record_path)
    monkeypatch.setattr(session, 'live_records', lambda: [])
    monkeypatch.setattr(launcher.shutil, 'which', lambda name: '/usr/bin/tmux')
    monkeypatch.setattr(launcher.socket, 'create_connection', lambda *a, **kw: nullcontext())
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(launcher.subprocess, 'run', run)
    class Attached(Exception):
        pass
    def attach(executable, command):
        commands.append(command)
        raise Attached()
    monkeypatch.setattr(launcher.os, 'execv', attach)
    if not tmux_name:
        def monitored(*args):
            raise Attached()
        monkeypatch.setattr(launcher.time, 'sleep', monitored)
    monkeypatch.setattr(launcher.sys, 'argv', ['helper', '8866', '/work', '/env/bin/delfin-voila'])
    with pytest.raises(Attached):
        launcher.main()
    if tmux_name:
        assert [command[1] for command in commands] == ['attach-session']
        assert commands[0][-1] == '=' + tmux_name
    else:
        assert commands == []
    output = capsys.readouterr().out.splitlines()[0]
    assert token not in output
    url = base64.b64decode(output.removeprefix('DELFIN_BROWSER:')).decode()
    assert url == f'http://127.0.0.1:{server_port}/voila/render/dashboard.ipynb?token={token}'


def test_remote_node_record_is_not_overwritten(monkeypatch):
    from delfin.agent import where
    monkeypatch.setattr(launcher.shutil, 'which', lambda name: '/usr/bin/tmux')
    monkeypatch.setattr(where, 'dashboard', lambda: {'host': 'other-login-node', 'running': None})
    monkeypatch.setattr(launcher.sys, 'argv', ['helper', '8866', '/work', '/env/bin/delfin-voila'])
    with pytest.raises(RuntimeError, match='another login node'):
        launcher.main()


def test_remote_port_search_skips_occupied_port():
    with socket.socket() as occupied:
        occupied.bind(('127.0.0.1', 0))
        port = occupied.getsockname()[1]
        selected = launcher.select_port(port)
        assert selected != port
        with socket.socket() as available:
            available.bind(('127.0.0.1', selected))


def test_remote_port_zero_uses_os_assignment():
    selected = launcher.select_port(0)
    assert 1024 <= selected <= 65535


@pytest.mark.parametrize('compatible', [False, True])
@pytest.mark.parametrize('source', ['repo-env', 'active-env', 'path', 'explicit-python', 'explicit-repo', 'home-software-dotvenv', 'home-software-venv'])
def test_bootstrap_discovers_environment_without_fixed_installation_path(tmp_path, source, compatible):
    import base64
    import os
    import subprocess
    import sys

    user_home = tmp_path / 'user home'
    user_home.mkdir()
    automatic = source.startswith('home-software-')
    repo = user_home / 'software/delfin' if automatic else tmp_path / 'repo with spaces'
    repo.mkdir(parents=True)
    for package in ['delfin', 'delfin/agent', 'delfin/dashboard']:
        path = repo / package
        path.mkdir(exist_ok=True)
        (path / '__init__.py').write_text('')
    for module in ['delfin/cli_voila.py', 'delfin/agent/where.py', 'delfin/dashboard/session.py']:
        (repo / module).write_text('')
    option = '--strict-port' if compatible else '--port'
    (repo / 'delfin/cli_voila.py').write_text(f'def main(argv):\n    print({option!r})\n')
    root = repo / 'env' if source in {'repo-env', 'explicit-repo'} else tmp_path / 'custom environment'
    if automatic: root = repo / ('.venv' if source.endswith('dotvenv') else 'venv')
    bin_dir = root / 'bin'
    bin_dir.mkdir(parents=True)
    python = bin_dir / 'python'
    python.symlink_to(sys.executable)
    (bin_dir / 'activate').write_text(f'export VIRTUAL_ENV={__import__("shlex").quote(str(root))}\n')
    env = dict(os.environ)
    env.pop('VIRTUAL_ENV', None)
    env.pop('CONDA_PREFIX', None)
    env['PATH'] = '/usr/bin:/bin'
    env['HOME'] = str(user_home)
    hint = ''
    if source == 'active-env':
        env['VIRTUAL_ENV'] = str(root)
    elif source == 'path':
        env['PATH'] = str(bin_dir) + ':' + env['PATH']
    elif source == 'explicit-python':
        hint = str(python)
    elif source == 'explicit-repo':
        hint = str(repo)
    payload = base64.b64encode(b'import os,sys; print("CHOSEN=" + sys.executable); print("ACTIVATED=" + os.environ.get("VIRTUAL_ENV", ""))').decode()
    script = Path(__file__).resolve().parents[1] / 'tools/windows-launcher/remote_bootstrap.sh'
    result = subprocess.run(['bash', str(script), '0', '' if automatic else str(repo), hint, payload, 'dashboard'],
                            env=env, text=True, capture_output=True, timeout=15)
    if not compatible:
        assert result.returncode == 1
        assert 'Server DELFIN is too old' in result.stderr
        assert 'CHOSEN=' not in result.stdout
        return
    assert result.returncode == 0, result.stderr
    assert 'CHOSEN=' + str(python) in result.stdout
    assert 'ACTIVATED=' + str(root) in result.stdout


@pytest.mark.parametrize('keep', [False, True])
def test_disconnect_cleanup_respects_keep_for_new_dashboard(monkeypatch, capsys, keep):
    from types import SimpleNamespace
    from delfin.agent import where
    monkeypatch.setattr(where, 'dashboard', lambda: {})
    monkeypatch.setattr(launcher.shutil, 'which', lambda name: '/usr/bin/tmux')
    monkeypatch.setattr(launcher, 'select_port', lambda preferred: 8866)
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=1 if command[1] == 'has-session' else 0)
    monkeypatch.setattr(launcher.subprocess, 'run', run)
    def disconnected(*args):
        raise RuntimeError('test disconnect')
    monkeypatch.setattr(launcher, 'run_dashboard_connection', disconnected)
    monkeypatch.setattr(launcher.sys, 'argv', ['helper','0','/work','/env/bin/python','/tmp/delfin-win-tunnel-test/dashboard.sock',str(int(keep)),'0'])
    with pytest.raises(RuntimeError, match='test disconnect'):
        launcher.main()
    assert any(cmd[1] == 'kill-session' for cmd in commands) is (not keep)
    created = next(cmd for cmd in commands if cmd[1] == 'new-session')
    assert ('--keep' in created[-1]) is keep
    assert 'DELFIN_CONNECTION_CLOSED:delfin-win-tunnel-test' in capsys.readouterr().out


def test_private_relay_forwards_and_flushes_half_close():
    import threading
    import uuid
    from pathlib import Path
    received = bytearray()
    payload = b'dashboard' * 30000
    with socket.socket() as server:
        server.bind(('127.0.0.1', 0))
        server.listen(1)
        def echo():
            peer, _ = server.accept()
            with peer:
                while True:
                    data = peer.recv(65536)
                    if not data:
                        break
                    received.extend(data)
                peer.sendall(received)
        worker = threading.Thread(target=echo)
        worker.start()
        path = Path('/tmp') / ('delfin-win-tunnel-' + str(uuid.uuid4())) / 'dashboard.sock'
        relay = launcher.DashboardRelay(str(path), server.getsockname()[1])
        try:
            assert path.parent.stat().st_mode & 0o777 == 0o700
            assert path.stat().st_mode & 0o777 == 0o600
            with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as client:
                client.settimeout(10)
                client.connect(str(path))
                client.sendall(payload)
                client.shutdown(socket.SHUT_WR)
                response = bytearray()
                while True:
                    data = client.recv(65536)
                    if not data:
                        break
                    response.extend(data)
            assert response == payload
        finally:
            relay.close()
            worker.join(timeout=5)
        assert not path.parent.exists()


def test_working_shell_uses_selected_environment_instead_of_tmux_server_environment(monkeypatch, tmp_path):
    import os
    import subprocess
    import sys
    repo = tmp_path / 'source repo'
    repo.mkdir()
    (repo / 'selected_module.py').write_text('VALUE = "selected-source"\n')
    environment = tmp_path / 'chosen environment'
    binary = environment / 'bin'
    binary.mkdir(parents=True)
    (binary / 'python').symlink_to(sys.executable)
    monkeypatch.setenv('PATH', str(binary) + ':/usr/bin:/bin')
    monkeypatch.setenv('VIRTUAL_ENV', str(environment))
    monkeypatch.setenv('PYTHONPATH', str(repo))
    monkeypatch.setenv('PS1', '(selected) $ ')
    command = launcher.working_shell_command(str(tmp_path))
    task = 'python -c \'import os,selected_module; print("ENV="+os.environ["VIRTUAL_ENV"]); print("SOURCE="+selected_module.VALUE)\'\nexit\n'
    result = subprocess.run(['/bin/bash', '-c', command], input=task, text=True,
                            env={'PATH':'/usr/bin:/bin'}, capture_output=True, timeout=15)
    assert result.returncode == 0, result.stderr
    assert 'ENV=' + str(environment) in result.stdout
    assert 'SOURCE=selected-source' in result.stdout
