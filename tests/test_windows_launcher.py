"""Checks for the portable launcher's security-sensitive return routing."""
import json
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

_ROOT = Path(__file__).resolve().parents[1] / 'tools/windows-launcher'


def _launcher_file(tmp_path, code):
    """A stand-in for the copied remote_launcher.py: the bootstrap gets its
    path as $4 and runs it as a plain script file."""
    folder = tmp_path / 'copied scripts'
    folder.mkdir(exist_ok=True)
    path = folder / 'remote_launcher.py'
    path.write_text(code)
    return str(path)


@pytest.mark.parametrize('browser_port', [8866, 58249])
def test_return_uses_only_latest_kernel_of_authenticated_server(browser_port):
    token = 'a' * 43
    record = {'port': 8866, 'token': token, 'url': '/voila/render/dashboard.ipynb',
              'resume_path': '/voila/render/return.ipynb'}
    host = socket.gethostname().split('.')[0]
    def session(name, stamp, port=browser_port, secret=token, node=host, browser_host='localhost'):
        return {'session_name': name, 'updated_at': stamp, 'host': node,
                'request_url': f'http://{browser_host}:{port}/?token={secret}'}
    rows = [session('old', 1), session('new & kept', 2),
            session('other-server', 100, port=9000, secret='b' * 43),
            session('foreign-url-host', 100, browser_host='attacker.invalid'),
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
    payload = _launcher_file(tmp_path, 'import os,sys; print("CHOSEN=" + sys.executable); print("ACTIVATED=" + os.environ.get("VIRTUAL_ENV", ""))')
    script = _ROOT / 'remote_bootstrap.sh'
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


def _installed_like_an_editable_finder(tmp_path):
    """A fake DELFIN reachable only through a meta-path finder appended to
    sys.meta_path (as setuptools' editable finder does), and a venv python.
    Returns (environment root, python, PYTHONPATH carrying the finder)."""
    import sys
    source = tmp_path / 'installed source'
    for package in ['delfin', 'delfin/agent', 'delfin/dashboard']:
        (source / package).mkdir(parents=True, exist_ok=True)
        (source / package / '__init__.py').write_text('')
    (source / 'delfin/__init__.py').write_text('__version__ = "test"\n')
    for module in ['delfin/agent/where.py', 'delfin/dashboard/session.py']:
        (source / module).write_text('')
    (source / 'delfin/cli_voila.py').write_text(
        'from delfin import __version__\n'
        'def main(argv):\n    print("--strict-port")\n')
    hook = tmp_path / 'finder hook'
    hook.mkdir()
    # The child must see ONLY this fake: on a machine with DELFIN installed
    # (CI, a developer's env) the real package -- a site-packages path or the
    # editable install's own finder -- would win and the test would neither
    # reproduce the shadowing nor prove the fix. sitecustomize runs after the
    # .pth files, so it can take those out; a folder WITHOUT __init__.py
    # (the ~/delfin install root under test) stays on sys.path.
    (hook / 'sitecustomize.py').write_text(
        'import os, sys\n'
        'sys.path[:] = [p for p in sys.path if not os.path.isfile('
        'os.path.join(p or ".", "delfin", "__init__.py"))]\n'
        'sys.meta_path[:] = [f for f in sys.meta_path '
        'if "editable" not in repr(f).lower()]\n'
        'from importlib.machinery import PathFinder\n'
        'class _Finder:\n'
        '    @classmethod\n'
        '    def find_spec(cls, name, path=None, target=None):\n'
        '        if name == "delfin":\n'
        f'            return PathFinder.find_spec(name, [{str(source)!r}])\n'
        '        return None\n'
        'sys.meta_path.append(_Finder)\n')
    root = tmp_path / 'venv'
    (root / 'bin').mkdir(parents=True)
    python = root / 'bin' / 'python'
    python.symlink_to(sys.executable)
    (root / 'bin' / 'activate').write_text(f'export VIRTUAL_ENV={__import__("shlex").quote(str(root))}\n')
    return root, python, str(hook)


def _run_bootstrap(env, home, payload_code, *, person=''):
    import subprocess
    payload = _launcher_file(home.parent, payload_code)
    script = _ROOT / 'remote_bootstrap.sh'
    return subprocess.run(['bash', str(script), '0', '', '', payload, 'dashboard', '', '0', '0', person],
                          env=env, text=True, capture_output=True, timeout=30, cwd=str(home))


def _bootstrap_env(tmp_path, home):
    import os
    root, python, hook = _installed_like_an_editable_finder(tmp_path)
    env = dict(os.environ)
    for key in ('CONDA_PREFIX', 'DELFIN_PERSON', 'DELFIN_RUNTIME_DIR', 'DELFIN_SETTINGS_DEFAULTS'):
        env.pop(key, None)
    env.update(PATH='/usr/bin:/bin', HOME=str(home), VIRTUAL_ENV=str(root), PYTHONPATH=hook)
    return env, python


def test_bootstrap_finds_delfin_although_home_holds_a_plain_delfin_folder(tmp_path):
    # The Raspberry Pi layout: ~/delfin is the install root (no __init__.py),
    # the package comes through the editable finder. From ~ a bare `python -c`
    # imported ~/delfin as an empty namespace package and found nothing.
    home = tmp_path / 'home'
    (home / 'delfin' / 'src').mkdir(parents=True)
    env, python = _bootstrap_env(tmp_path, home)
    result = _run_bootstrap(env, home, 'import delfin, sys; print("PKG=" + str(delfin.__version__)); print("CHOSEN=" + sys.executable)')
    assert result.returncode == 0, result.stderr
    assert 'PKG=test' in result.stdout
    assert 'CHOSEN=' + str(python) in result.stdout


def test_people_mode_gives_each_person_a_home_of_their_own(tmp_path):
    import stat
    home = tmp_path / 'home'
    (home / '.config/delfin').mkdir(parents=True)
    (home / '.config/delfin/launcher-people.conf').write_text('')
    (home / '.delfin').mkdir()
    (home / '.delfin/credentials.json').write_text('{}')
    (home / '.delfin_settings.json').write_text('{"ui": {"tabs": {"hidden": ["submit"]}}}')
    (home / '.delfin/agent_sessions').mkdir()
    env, _ = _bootstrap_env(tmp_path, home)
    probe = ('import os; print("HOME=" + os.environ["HOME"]); print("RUN=" + os.environ["DELFIN_RUNTIME_DIR"]);'
             ' print("WHO=" + os.environ["DELFIN_PERSON"]); print("CWD=" + os.getcwd());'
             ' print("DEFAULTS=" + os.environ["DELFIN_SETTINGS_DEFAULTS"])')
    results = {name: _run_bootstrap(env, home, probe, person=name) for name in ('anna', 'bert')}
    for name, result in results.items():
        assert result.returncode == 0, result.stderr
        person_home = home / 'people' / name
        assert f'HOME={person_home}' in result.stdout
        assert f'RUN={person_home}/.run' in result.stdout
        assert f'WHO={name}' in result.stdout
        assert f'CWD={person_home}/agent_workspace' in result.stdout
        assert stat.S_IMODE(person_home.stat().st_mode) == 0o700
        # Shared by link: the account's credentials; private: sessions. The
        # DELFIN settings (user_settings.SETTINGS_FILE_NAME) are layered: the
        # account's file is named as the defaults, the person has their own.
        from delfin.user_settings import SETTINGS_FILE_NAME, SETTINGS_DEFAULTS_ENV
        assert (person_home / '.delfin/credentials.json').resolve() == (home / '.delfin/credentials.json').resolve()
        own = person_home / SETTINGS_FILE_NAME
        assert own.is_file() and not own.is_symlink()
        assert json.loads(own.read_text()) == {}
        assert stat.S_IMODE(own.stat().st_mode) == 0o600
        assert SETTINGS_DEFAULTS_ENV == 'DELFIN_SETTINGS_DEFAULTS'
        assert f'DEFAULTS={home / SETTINGS_FILE_NAME}' in result.stdout
        assert not (person_home / '.delfin/agent_sessions').exists()
    assert json.loads((home / '.delfin_settings.json').read_text()) == {'ui': {'tabs': {'hidden': ['submit']}}}


def test_people_mode_is_off_unless_the_account_owner_enables_it(tmp_path):
    home = tmp_path / 'home'
    home.mkdir()
    env, _ = _bootstrap_env(tmp_path, home)
    result = _run_bootstrap(env, home, 'import os; print("HOME=" + os.environ["HOME"])', person='anna')
    assert result.returncode == 0, result.stderr
    assert f'HOME={home}' in result.stdout
    assert not (home / 'people').exists()


@pytest.mark.parametrize('person', ['../evil', 'a b', '-x', 'Anna'])
def test_people_mode_refuses_a_name_that_is_not_a_plain_name(tmp_path, person):
    home = tmp_path / 'home'
    (home / '.config/delfin').mkdir(parents=True)
    (home / '.config/delfin/launcher-people.conf').write_text('')
    env, _ = _bootstrap_env(tmp_path, home)
    result = _run_bootstrap(env, home, 'print("STARTED")', person=person)
    assert result.returncode == 1
    assert 'Invalid person name' in result.stderr
    assert 'STARTED' not in result.stdout
    assert not (home / 'people').exists()


def test_tmux_gets_home_and_delfin_settings_of_this_person(monkeypatch):
    monkeypatch.setenv('HOME', '/srv/account/people/anna')
    monkeypatch.setenv('DELFIN_PERSON', 'anna')
    monkeypatch.setenv('DELFIN_RUNTIME_DIR', '/srv/account/people/anna/.run')
    words = launcher.shell_environment().split()
    assert 'HOME=/srv/account/people/anna' in words
    assert 'DELFIN_PERSON=anna' in words
    assert 'DELFIN_RUNTIME_DIR=/srv/account/people/anna/.run' in words


def test_tmux_gets_the_settings_defaults_of_the_account(monkeypatch):
    # The dashboard and the agent service the person starts read the
    # account's settings under their own; the variable must reach them.
    monkeypatch.setenv('HOME', '/srv/account/people/anna')
    monkeypatch.setenv('DELFIN_SETTINGS_DEFAULTS', '/srv/account/.delfin_settings.json')
    words = launcher.shell_environment().split()
    assert 'DELFIN_SETTINGS_DEFAULTS=/srv/account/.delfin_settings.json' in words
    assert 'DELFIN_SETTINGS_DEFAULTS=/srv/account/.delfin_settings.json' in launcher.working_shell_command('/w')


def test_new_dashboard_tmux_session_is_named_for_the_person(monkeypatch):
    from types import SimpleNamespace
    from delfin.agent import where
    monkeypatch.setenv('DELFIN_PERSON', 'anna')
    monkeypatch.setattr(where, 'dashboard', lambda: {})
    monkeypatch.setattr(launcher.shutil, 'which', lambda name: '/usr/bin/tmux')
    monkeypatch.setattr(launcher, 'select_port', lambda preferred: 8866)
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(launcher.subprocess, 'run', run)
    def stop(*args):
        raise RuntimeError('stop here')
    monkeypatch.setattr(launcher, 'run_dashboard_connection', stop)
    monkeypatch.setattr(launcher.sys, 'argv', ['helper', '8866', '/work', '/env/bin/python', '', '1', '0'])
    with pytest.raises(RuntimeError, match='stop here'):
        launcher.main()
    created = next(cmd for cmd in commands if cmd[1] == 'new-session')
    assert created[created.index('-s') + 1].startswith('delfin-win-anna-8866-')
    assert 'sys.path[:]' in created[-1]


def test_a_server_can_make_keep_session_its_default(tmp_path):
    # The payload receives keep as argv[5] (remote_launcher.main). The app
    # sends 0 here; the account's launcher-keep.conf turns it on, and only
    # where that file exists.
    home = tmp_path / 'home'
    home.mkdir()
    env, _ = _bootstrap_env(tmp_path, home)
    probe = 'import sys; print("KEEP=" + sys.argv[5])'
    assert 'KEEP=0' in _run_bootstrap(env, home, probe).stdout
    (home / '.config/delfin').mkdir(parents=True)
    (home / '.config/delfin/launcher-keep.conf').write_text('')
    result = _run_bootstrap(env, home, probe)
    assert result.returncode == 0, result.stderr
    assert 'KEEP=1' in result.stdout


def _people_home(tmp_path):
    home = tmp_path / 'home'
    (home / '.config/delfin').mkdir(parents=True)
    (home / '.config/delfin/launcher-people.conf').write_text('')
    (home / '.delfin_settings.json').write_text('{"agent": {"permission_profile": "all_free"}}')
    person_home = home / 'people' / 'anna'
    person_home.mkdir(parents=True)
    return home, person_home


def test_people_mode_turns_a_settings_link_into_an_own_empty_file(tmp_path):
    # The earlier layout linked the person's settings to the account's. The
    # link's content was the account's, which is now the defaults layer: the
    # link goes, an empty personal file takes its place, the account's file
    # is untouched, and a second connection changes nothing more.
    home, person_home = _people_home(tmp_path)
    own = person_home / '.delfin_settings.json'
    own.symlink_to(home / '.delfin_settings.json')
    env, _ = _bootstrap_env(tmp_path, home)
    for _ in range(2):
        result = _run_bootstrap(env, home, 'print("OK")', person='anna')
        assert result.returncode == 0, result.stderr
    assert own.is_file() and not own.is_symlink()
    assert json.loads(own.read_text()) == {}
    assert json.loads((home / '.delfin_settings.json').read_text()) == {'agent': {'permission_profile': 'all_free'}}
    assert not list(person_home.glob('.delfin_settings.json.before-shared-*'))


def test_people_mode_keeps_an_own_settings_file_and_does_not_apply_a_kept_copy(tmp_path):
    # An own file is the person's choices: kept as it is, not moved aside.
    # A copy kept by the earlier layout (*.before-shared-*) stays where it is
    # and is not merged back.
    home, person_home = _people_home(tmp_path)
    own = person_home / '.delfin_settings.json'
    own.write_text('{"agent": {"permission_profile": "plan"}}')
    kept = person_home / '.delfin_settings.json.before-shared-20261101000000'
    kept.write_text('{"agent": {"model": "old"}}')
    env, _ = _bootstrap_env(tmp_path, home)
    result = _run_bootstrap(env, home, 'import os; print("DEFAULTS=" + os.environ["DELFIN_SETTINGS_DEFAULTS"])', person='anna')
    assert result.returncode == 0, result.stderr
    assert f'DEFAULTS={home}/.delfin_settings.json' in result.stdout
    assert json.loads(own.read_text()) == {'agent': {'permission_profile': 'plan'}}
    assert json.loads(kept.read_text()) == {'agent': {'model': 'old'}}
    assert list(person_home.glob('.delfin_settings.json.before-shared-*')) == [kept]


def test_people_mode_leaves_a_settings_link_elsewhere_alone(tmp_path):
    # Only the account's link is the earlier layout; a link the person made
    # to another file is their own setup, and a dangling one is not written
    # through.
    home, person_home = _people_home(tmp_path)
    own = person_home / '.delfin_settings.json'
    own.symlink_to(tmp_path / 'elsewhere.json')
    env, _ = _bootstrap_env(tmp_path, home)
    result = _run_bootstrap(env, home, 'print("OK")', person='anna')
    assert result.returncode == 0, result.stderr
    assert own.is_symlink() and not (tmp_path / 'elsewhere.json').exists()


def test_the_account_without_people_mode_gets_no_settings_layer(tmp_path):
    home = tmp_path / 'home'
    home.mkdir()
    env, _ = _bootstrap_env(tmp_path, home)
    result = _run_bootstrap(env, home, 'import os; print("DEFAULTS=" + os.environ.get("DELFIN_SETTINGS_DEFAULTS", "-"))', person='anna')
    assert result.returncode == 0, result.stderr
    assert 'DEFAULTS=-' in result.stdout


def _run_bootstrap_with_key(env, home, payload_code, key):
    import subprocess
    payload = _launcher_file(home.parent, payload_code)
    script = _ROOT / 'remote_bootstrap.sh'
    return subprocess.run(['bash', str(script), '0', '', '', payload, 'dashboard', '', '0', '0', '', key],
                          env=env, text=True, capture_output=True, timeout=30, cwd=str(home))


_KEY = 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIGxVZr1q3F0c0o1u2y3m4n5p6q7r8s9t0u1v2w3x4y5z delfin-launcher-anna'


def test_remember_login_registers_the_pcs_key_once_and_restricted(tmp_path):
    import stat
    home = tmp_path / 'home'
    home.mkdir()
    env, _ = _bootstrap_env(tmp_path, home)
    for _ in range(2):
        result = _run_bootstrap_with_key(env, home, 'print("OK")', _KEY)
        assert result.returncode == 0, result.stderr
    lines = (home / '.ssh/authorized_keys').read_text().splitlines()
    assert lines == ['restrict,port-forwarding,pty ' + _KEY]
    assert stat.S_IMODE((home / '.ssh').stat().st_mode) == 0o700
    assert stat.S_IMODE((home / '.ssh/authorized_keys').stat().st_mode) == 0o600


@pytest.mark.parametrize('bad', [
    'ssh-rsa AAAAB3NzaC1yc2E delfin',                       # only ed25519
    _KEY + '\nssh-ed25519 AAAAother attacker',              # a second line
    'command="sh" ' + _KEY,                                 # options smuggled in
    _KEY.replace('anna', 'anna; rm -rf ~'),                 # shell text in the comment
])
def test_remember_login_refuses_anything_but_one_plain_ed25519_key(tmp_path, bad):
    home = tmp_path / 'home'
    home.mkdir()
    env, _ = _bootstrap_env(tmp_path, home)
    result = _run_bootstrap_with_key(env, home, 'print("OK")', bad)
    assert result.returncode == 0, result.stderr
    assert not (home / '.ssh/authorized_keys').exists()


def test_a_server_can_refuse_remembered_logins(tmp_path):
    home = tmp_path / 'home'
    (home / '.config/delfin').mkdir(parents=True)
    (home / '.config/delfin/launcher-no-keys.conf').write_text('')
    env, _ = _bootstrap_env(tmp_path, home)
    assert _run_bootstrap_with_key(env, home, 'print("OK")', _KEY).returncode == 0
    assert not (home / '.ssh/authorized_keys').exists()


def test_the_app_and_the_server_agree_on_the_login_key_argument():
    # DELFIN.ps1 cannot run here; pin what the server side depends on: the
    # public key is the argument after the person name (bootstrap $10), it is
    # offered to ssh with -i, and the checkbox is saved with the connection.
    root = Path(__file__).resolve().parents[1] / 'tools/windows-launcher'
    app = (root / 'DELFIN.ps1').read_text(encoding='utf-8')
    server = (root / 'remote_bootstrap.sh').read_text(encoding='utf-8')
    assert "(Shell-Quote (Get-PersonName)) + ' ' + (Shell-Quote $publicKey)" in app
    assert 'person=${9:-}' in server and 'login_key=${10:-}' in server
    assert "$options += @('-i',$loginKey.Key)" in app
    assert 'RememberLogin=$rememberToggle.Checked' in app and 'RememberLogin=$true' in app


def _check_script():
    """The server-side check DELFIN.ps1 sends (Get-ScriptCheckCommand), as
    bash receives it: the PowerShell literal with '' read as '."""
    import re
    app = (_ROOT / 'DELFIN.ps1').read_text(encoding='utf-8')
    found = re.findall(r"^\s*\$checkScript = '((?:[^']|'')*)'$", app, re.M)
    assert len(found) == 1
    return found[0].replace("''", "'")


def _digest(path):
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_the_check_reports_a_missing_or_changed_copy_and_accepts_an_exact_one(tmp_path):
    # Exit 10 makes the app copy the two files (scp); exit 0 skips the copy.
    import shutil
    import stat
    import subprocess
    home = tmp_path / 'account home'
    home.mkdir()
    hb, hl = _digest(_ROOT / 'remote_bootstrap.sh'), _digest(_ROOT / 'remote_launcher.py')
    folder = '.cache/delfin-launcher/' + hb[:16] + hl[:16]
    def check():
        return subprocess.run(['bash', '-c', _check_script(), 'delfin', folder, hb, hl],
                              env={'PATH': '/usr/bin:/bin', 'HOME': str(home)},
                              capture_output=True, text=True, timeout=15).returncode
    assert check() == 10
    copy = home / folder
    assert copy.is_dir() and stat.S_IMODE(copy.stat().st_mode) == 0o700
    for name in ('remote_bootstrap.sh', 'remote_launcher.py'):
        shutil.copyfile(_ROOT / name, copy / name)
    assert check() == 0
    # An interrupted copy is not mistaken for a complete one.
    (copy / 'remote_launcher.py').write_bytes((_ROOT / 'remote_launcher.py').read_bytes()[:100])
    assert check() == 10


def test_the_login_shell_sources_the_copied_bootstrap_with_its_arguments(tmp_path):
    # What Get-RemoteCommand sends: bash -ilc '. ~/<folder>/remote_bootstrap.sh'
    # delfin ARGS, with $4 = '~/<folder>/remote_launcher.py'. -i and -l are
    # left out here so this machine's startup files stay out of the test.
    import shutil
    import subprocess
    home = tmp_path / 'home'
    home.mkdir()
    env, python = _bootstrap_env(tmp_path, home)
    folder = '.cache/delfin-launcher/' + 'a' * 32
    (home / folder).mkdir(parents=True)
    shutil.copyfile(_ROOT / 'remote_bootstrap.sh', home / folder / 'remote_bootstrap.sh')
    (home / folder / 'remote_launcher.py').write_text(
        'import sys; print("CHOSEN=" + sys.executable); print("ARGS=" + "|".join(sys.argv[1:]))')
    result = subprocess.run(['bash', '-c', '. ~/' + folder + '/remote_bootstrap.sh', 'delfin',
                             '0', '', '', '~/' + folder + '/remote_launcher.py', 'dashboard', '', '1', '0'],
                            env=env, cwd=str(home), capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert 'CHOSEN=' + str(python) in result.stdout
    assert f'ARGS=0|{home}|{python}||1|0' in result.stdout


def test_the_bootstrap_refuses_a_missing_launcher_file(tmp_path):
    home = tmp_path / 'home'
    home.mkdir()
    env, _ = _bootstrap_env(tmp_path, home)
    import subprocess
    result = subprocess.run(['bash', str(_ROOT / 'remote_bootstrap.sh'), '0', '', '',
                             str(tmp_path / 'absent.py'), 'dashboard'],
                            env=env, cwd=str(home), capture_output=True, text=True, timeout=30)
    assert result.returncode == 1
    assert 'launcher script is missing' in result.stderr


def test_the_app_builds_the_remote_command_from_plain_file_paths():
    import re
    app = (_ROOT / 'DELFIN.ps1').read_text(encoding='utf-8')
    body = re.search(r'^function Get-RemoteCommand\b.*?^\}', app, re.M | re.S).group(0)
    assert "$bootstrap = '. ~/' + $version.Folder + '/remote_bootstrap.sh'" in body
    assert "$launcher = '~/' + $version.Folder + '/remote_launcher.py'" in body
    assert "'bash -ilc ' + (Shell-Quote $bootstrap)" in body
    assert "(Shell-Quote $launcher)" in body
    for word in ('base64', 'Base64', 'gzip', 'GZip', 'eval', 'exec', 'ReadAllBytes'):
        assert word not in body, word
    # The two scripts reach the server as files (scp) after a short check.
    sync = re.search(r'^function Sync-ServerScripts\b.*?^\}', app, re.M | re.S).group(0)
    assert "'System32\\OpenSSH\\scp.exe'" in sync
    assert "@('remote_bootstrap.sh','remote_launcher.py',($target + ':' + $version.Folder + '/'))" in sync
    assert 'Get-ScriptCheckCommand $version' in sync
    # Both connections carry the same OpenSSH protections and the login key.
    assert "(Get-SSHOptions 'Check') + $identity" in sync and "(Get-SSHOptions 'Copy') + $identity" in sync
    server = (_ROOT / 'remote_bootstrap.sh').read_text()
    assert 'exec "$python" "$launcher_script"' in server
    assert 'base64' not in server and 'eval ' not in server


def _package_names():
    import importlib.util as util
    spec = util.spec_from_file_location('windows_launcher_build', _ROOT / 'build_package.py')
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NAMES


def test_the_package_holds_no_compiled_encoded_or_hidden_code():
    import re
    names = _package_names()
    assert not [n for n in names if n.endswith(('.cs', '.exe', '.dll', '.vbs', '.js'))]
    assert not [p.name for p in _ROOT.iterdir() if p.suffix in ('.cs', '.exe', '.dll')]
    forbidden = [
        (r'Add-Type\s+(-TypeDefinition|-MemberDefinition|-Path|-Language)', 'Add-Type with C# source'),
        (r'DllImport|GetDelegateForFunctionPointer|Reflection\.Emit|\[Reflection\.Assembly\]::Load', 'Win32 calls or loaded code'),
        (r'-EncodedCommand|-enc\b|-ec\b', 'encoded command'),
        (r'(?i)-WindowStyle\s+Hidden|-w\s+hidden', 'hidden window'),
        (r'(?i)csc\.exe|Invoke-Expression|\biex\b', 'compiler or expression evaluation'),
        (r'GZipStream|DeflateStream', 'packed data'),
        (r'(?i)DownloadString|DownloadFile|Invoke-WebRequest\s+-Uri\s+["\']?https?://(?!127\.0\.0\.1)', 'download'),
        (r'(?i)CurrentVersion\\Run\b|\\Startup\b|Set-MpPreference|ExecutionPolicy\s+(Unrestricted|Bypass)|Set-ExecutionPolicy', 'persistence or security change'),
    ]
    for name in names:
        if not name.endswith(('.ps1', '.cmd')):
            continue
        text = (_ROOT / name).read_text(encoding='utf-8')
        for pattern, what in forbidden:
            assert not re.search(pattern, text), f'{name}: {what}'
        # Base64 appears only for the browser URL the server reports
        # (Get-BrowserUrl decodes it) and the test URLs Test-Launcher encodes.
        if name != 'Test-Launcher.ps1':
            assert 'ToBase64String' not in text, name
    for name in ('remote_bootstrap.sh', 'remote_launcher.py'):
        text = (_ROOT / name).read_text()
        assert 'b64decode' not in text and 'exec(' not in text and 'eval ' not in text, name


def test_every_shipped_text_file_is_plain_text():
    for name in _package_names():
        if name.endswith(('.png', '.ico')):
            continue
        data = (_ROOT / name).read_bytes()
        text = data.decode('utf-8')
        assert '\x00' not in text, name
        # No long run of base64-like characters: nothing packed in a script.
        import re
        assert not re.search(r'[A-Za-z0-9+/=]{200,}', text), name


def test_install_creates_standard_shortcuts_to_windows_powershell():
    install = (_ROOT / 'Install.ps1').read_text(encoding='utf-8')
    assert '$shell = New-Object -ComObject WScript.Shell' in install
    assert "$shell.CreateShortcut((Join-Path ([Environment]::GetFolderPath('Desktop')) 'DELFIN.lnk'))" in install
    assert '$shortcut.TargetPath = "$env:SystemRoot\\System32\\WindowsPowerShell\\v1.0\\powershell.exe"' in install
    assert "$shortcut.Arguments = '-NoProfile -ExecutionPolicy RemoteSigned -File \"' + $app + '\"'" in install
    assert '$shortcut.WindowStyle = 7' in install
    assert "$shortcut.IconLocation = (Join-Path $destination 'DELFIN.ico') + ',0'" in install
    assert '$shortcut.Save()' in install
    assert "Join-Path $programs 'DELFIN.lnk'" in install
    copied = install.split('# Files of earlier versions')[0]
    assert 'Build-Starter' not in copied and 'Starter.cs' not in copied and 'Taskbar.cs' not in copied
    assert '& (Join-Path' not in install and 'DELFIN.exe' not in install.replace("@('DELFIN.exe',", '')
    # An update removes the starter program an earlier version compiled.
    assert "foreach ($old in @('DELFIN.exe','Starter.cs','Build-Starter.ps1','Taskbar.cs'))" in install
    # The app starts its dashboard console as a plain Windows PowerShell.
    app = (_ROOT / 'DELFIN.ps1').read_text(encoding='utf-8')
    assert "@('-NoProfile','-ExecutionPolicy','RemoteSigned','-File',$PSCommandPath,'-Mode','Dashboard'" in app
    assert 'Start-Process -FilePath $powershell -ArgumentList $arguments' in app


def test_the_checksums_match_the_shipped_files():
    import hashlib
    rows = (_ROOT / 'SHA256SUMS.txt').read_text(encoding='ascii').splitlines()
    listed = dict(reversed(row.split('  ', 1)) for row in rows)
    assert sorted(listed) == sorted(_package_names())
    for name, digest in listed.items():
        assert hashlib.sha256((_ROOT / name).read_bytes()).hexdigest() == digest, name
