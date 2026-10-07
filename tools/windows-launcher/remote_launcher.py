"""Executed in memory through the authenticated SSH connection; no upload/install."""
import base64
import os
import re
import signal
import shlex
import shutil
import socket
import stat
import subprocess
import sys
import time
from urllib.parse import parse_qs, urlencode, urlsplit


def return_url(record, sessions):
    """Return only this server's newest kept kernel, never a foreign record's URL."""
    port = int(record['port'])
    query = {'token': record['token']}
    eligible = []
    for item in sessions:
        try:
            url = urlsplit(item.get('request_url', ''))
            if (item.get('host') == socket.gethostname().split('.')[0]
                    and url.scheme in ('http', 'https') and url.hostname in ('localhost', '127.0.0.1', '::1')
                    and url.port is not None and 1024 <= url.port <= 65535
                    and parse_qs(url.query).get('token') == [record['token']]
                    and item.get('session_name')):
                eligible.append(item)
        except ValueError:
            continue
    # Go directly to the render route: root redirects can drop the token.
    render = urlsplit(record.get('url', ''))
    path = render.path if (not render.scheme and not render.netloc
                           and render.path.startswith('/voila/render/')) else '/'
    if eligible and record.get('resume_path', '').startswith('/voila/render/'):
        def updated(item):
            try:
                return float(item.get('updated_at', 0))
            except (TypeError, ValueError):
                return 0.0
        newest = max(eligible, key=updated)
        query['session'] = newest['session_name']
        path = record['resume_path']
    return f'http://127.0.0.1:{port}{path}?{urlencode(query)}'


def shell_environment():
    """Pass the selected environment explicitly; tmux may have an older server environment."""
    return ' '.join(shlex.quote(key + '=' + os.environ[key])
                    for key in ('PATH', 'VIRTUAL_ENV', 'CONDA_PREFIX', 'PYTHONPATH', 'PS1')
                    if key in os.environ)


def working_shell_command(directory):
    return 'cd ' + shlex.quote(directory) + ' && exec env ' + shell_environment() + ' bash --norc -i'


def select_port(preferred):
    """Choose a bindable loopback port; dashboard strict-port closes the later race."""
    candidates = range(preferred, min(preferred + 100, 65535) + 1) if preferred else [0]
    for candidate in candidates:
        with socket.socket() as probe:
            try:
                probe.bind(('127.0.0.1', candidate))
                return probe.getsockname()[1]
            except OSError:
                continue
    with socket.socket() as probe:
        probe.bind(('127.0.0.1', 0))
        return probe.getsockname()[1]


class DashboardRelay:
    """Private per-connection Unix socket reached by one SSH TCP->Unix forward."""
    def __init__(self, path, port):
        import threading
        import uuid
        from pathlib import Path
        candidate = Path(path)
        if candidate.name != 'dashboard.sock' or candidate.parent.parent != Path('/tmp'):
            raise RuntimeError('Invalid private tunnel path.')
        prefix = 'delfin-win-tunnel-'
        if not candidate.parent.name.startswith(prefix):
            raise RuntimeError('Invalid private tunnel directory.')
        uuid.UUID(candidate.parent.name[len(prefix):])
        self.path, self.port = candidate, port
        self.closed = False
        self.clients = set()
        self.lock = threading.Lock()
        os.mkdir(candidate.parent, 0o700)
        try:
            self.listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            self.listener.bind(str(candidate))
            os.chmod(candidate, 0o600)
            self.listener.listen(32)
            self.listener.settimeout(0.5)
        except BaseException:
            self.close()
            raise
        threading.Thread(target=self.serve, daemon=True).start()

    def serve(self):
        import threading
        while not self.closed:
            try:
                client, _ = self.listener.accept()
            except socket.timeout:
                continue
            except OSError:
                return
            threading.Thread(target=self.forward, args=(client,), daemon=True).start()

    def forward(self, client):
        import select
        upstream = None
        try:
            upstream = socket.create_connection(('127.0.0.1', self.port), timeout=10)
            client.setblocking(False)
            upstream.setblocking(False)
            with self.lock:
                if self.closed:
                    return
                self.clients.update((client, upstream))
            pending = {client: bytearray(), upstream: bytearray()}
            peers = {client: upstream, upstream: client}
            ended, shut = set(), set()
            while not self.closed:
                readers = [s for s in peers if s not in ended and len(pending[peers[s]]) < 262144]
                writers = [s for s in peers if pending[s]]
                if len(ended) == 2 and not writers:
                    return
                readable, writable, _ = select.select(readers, writers, [], 0.5)
                for src in readable:
                    chunk = src.recv(65536)
                    if not chunk:
                        ended.add(src)
                    else:
                        pending[peers[src]].extend(chunk)
                for dst in writable:
                    sent = dst.send(pending[dst])
                    del pending[dst][:sent]
                for src in ended:
                    dst = peers[src]
                    if not pending[dst] and dst not in shut:
                        dst.shutdown(socket.SHUT_WR)
                        shut.add(dst)
        except OSError:
            pass
        finally:
            with self.lock:
                self.clients.discard(client)
                if upstream:
                    self.clients.discard(upstream)
            client.close()
            if upstream:
                upstream.close()

    def close(self):
        self.closed = True
        if hasattr(self, 'listener'):
            self.listener.close()
        with self.lock:
            for client in self.clients:
                client.close()
        try:
            self.path.unlink(missing_ok=True)
            self.path.parent.rmdir()
        except OSError:
            pass


def main():
    from delfin.agent import where
    from delfin.dashboard import session as dashboard_session

    preferred, directory, python = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    if preferred != 0 and not 1024 <= preferred <= 65535:
        raise RuntimeError('Invalid starting port.')
    tunnel_path = sys.argv[4] if len(sys.argv) > 4 else ''
    keep = len(sys.argv) <= 5 or sys.argv[5] == '1'
    working = len(sys.argv) > 6 and sys.argv[6] == '1'
    name = f'delfin-win-{preferred}'
    port = preferred
    tmux = shutil.which('tmux')
    record = where.dashboard()
    if record and record.get('running') is None and record.get('host') != socket.gethostname():
        raise RuntimeError('Dashboard is on another login node. Connect to this exact node: ' + str(record.get('host')))
    existing = record.get('running') is True
    if existing:
        name = str(record.get('tmux') or '')
        if name and not re.fullmatch(r'[A-Za-z0-9_-]{1,128}', name):
            raise RuntimeError('Unsupported tmux session name in dashboard metadata.')
        port = int(record['port'])
    if not tmux and (not existing or name):
        raise RuntimeError('tmux is missing on the server. See README.')
    # Unverified leftover tmux sessions are not dashboards. Use a new unique name.
    exists = existing
    created_here = not exists
    if created_here:
        import uuid
        name += '-' + uuid.uuid4().hex
    def disconnected(signum, frame):
        raise KeyboardInterrupt()
    if tunnel_path:
        signal.signal(signal.SIGHUP, disconnected)
        signal.signal(signal.SIGTERM, disconnected)
    try:
        if created_here:
            port = select_port(preferred)
            command = ('cd ' + shlex.quote(directory)
                       + ' && exec env DELFIN_AGENT_SANDBOX=auto DELFIN_AGENT_SANDBOX_NETWORK=0 '
                       + ('DELFIN_KEEP_SESSIONS=1 ' if keep else 'DELFIN_KEEP_SESSIONS=0 ')
                       + shell_environment() + ' ' + shlex.quote(python)
                       + ' -c ' + shlex.quote('from delfin.cli_voila import main; raise SystemExit(main())')
                       + f' --no-browser --ip 127.0.0.1 --strict-port --port {port}'
                       + (' --keep --resume latest' if keep else ''))
            subprocess.run([tmux, 'new-session', '-d', '-s', name, '-x', '240', '-y', '50',
                            'bash -c ' + shlex.quote(command)], check=True)
        run_dashboard_connection(where, dashboard_session, tmux, name, tunnel_path, keep, working, directory, created_here)
    finally:
        if created_here and not keep:
            subprocess.run([tmux, 'kill-session', '-t', '=' + name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if tunnel_path:
            # Tell Windows cleanup is complete even if an unrelated descendant holds the PTY open.
            print('DELFIN_CONNECTION_CLOSED:' + os.path.basename(os.path.dirname(tunnel_path)), flush=True)


def run_dashboard_connection(where, dashboard_session, tmux, name, tunnel_path, keep, working, directory, created_here):
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        record = where.dashboard()
        if (record.get('running') is True and record.get('tmux') == name
                and 1024 <= int(record.get('port', 0)) <= 65535):
            port = int(record['port'])
            path = where.record_path()
            for private in (path, path.parent):
                info = private.stat()
                if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
                    raise RuntimeError('Dashboard metadata is not private (0600/0700 required).')
            try:
                with socket.create_connection(('127.0.0.1', port), timeout=1):
                    pass
            except OSError:
                time.sleep(0.5)
                continue
            token = record.get('token', '')
            if len(token) < 32 or any(c.isspace() for c in token):
                raise RuntimeError('No valid dashboard access token.')
            url = return_url(record, dashboard_session.live_records())
            # Private handoff over SSH, consumed in memory by the Windows launcher.
            relay = DashboardRelay(tunnel_path, port) if tunnel_path else None
            print('DELFIN_BROWSER:' + base64.b64encode(url.encode()).decode(), flush=True)
            try:
                if name:
                    if working:
                        result = subprocess.run([tmux, 'list-windows', '-t', '=' + name, '-F', '#{window_name}'], capture_output=True, text=True, check=True)
                        if 'DELFIN-work' not in result.stdout.splitlines():
                            shell = working_shell_command(directory)
                            subprocess.run([tmux, 'new-window', '-d', '-t', '=' + name, '-n', 'DELFIN-work', shell], check=True)
                        print('DELFIN: Working terminal is a tmux window. Switch with Ctrl+B, then N.', flush=True)
                    print('DELFIN: Keep session is ' + ('ON.' if keep else 'OFF.'), flush=True)
                    if not created_here:
                        print('DELFIN: Reusing a pre-existing dashboard; its lifetime is unchanged.', flush=True)
                    if not tunnel_path:
                        os.execv(tmux, [tmux, 'attach-session', '-t', '=' + name])
                    subprocess.run([tmux, 'attach-session', '-t', '=' + name], check=False)
                    return
                print('DELFIN: Existing dashboard has no tmux session. Monitoring it; keep its original terminal open.', flush=True)
                original_pid = record.get('pid')
                original_start = record.get('proc_start')
                while True:
                    time.sleep(2)
                    current = where.dashboard()
                    if (current.get('running') is not True or current.get('pid') != original_pid
                            or current.get('proc_start') != original_start):
                        return
            finally:
                if relay:
                    relay.close()
        if name and subprocess.run([tmux, 'has-session', '-t', '=' + name],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode:
            raise RuntimeError('Dashboard startup failed. delfin-voila --strict-port must be available; see README.')
        time.sleep(0.5)
    raise RuntimeError('Dashboard was not ready in time. The tmux session remains available for diagnosis.')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('DELFIN: Dashboard connection closed.', file=sys.stderr)
        sys.exit(0)
    except Exception as exc:
        print(f'DELFIN: {exc}', file=sys.stderr, flush=True)
        sys.exit(1)
