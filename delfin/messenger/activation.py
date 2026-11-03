"""Session opt-in owned by the DELFIN user, never by the phone gateway."""
from __future__ import annotations

import hashlib
import json
import os
import secrets
from pathlib import Path

_DIR = Path.home() / '.delfin' / 'messenger_sessions'


def _path(sid):
    return _DIR / (hashlib.sha256(sid.encode()).hexdigest() + '.json')


def record(sid):
    try:
        value = json.loads(_path(sid).read_text())
        return value if value.get('session_id') == sid else {}
    except (OSError, ValueError, AttributeError):
        return {}


def enabled(sid):
    return record(sid).get('enabled') is True


def records():
    result = []
    for path in _DIR.glob('*.json'):
        try:
            value = json.loads(path.read_text())
            if value.get('enabled') is True and value.get('session_id'):
                result.append(value)
        except (OSError, ValueError, AttributeError):
            continue
    return result


def configure(sid, args='', *, endpoint=None):
    if not sid:
        return 'Noch keine Sitzung. Zuerst eine Sitzung starten.'
    args = args.strip().lower()
    if args not in {'', 'on', 'off', 'status'}:
        return 'Verwendung: /sync on | off | status. Gilt nur für diese Sitzung.'
    if args == 'status':
        return 'SimpleX für diese Sitzung: ' + ('aktiviert.' if enabled(sid) else 'ausgeschaltet.')
    from delfin.agent.state_paths import write_text_atomic
    old = record(sid)
    on = args != 'off'
    same_endpoint = endpoint is None or (old.get('endpoint') or '') == endpoint
    value = {'session_id': sid, 'enabled': on,
             'generation': old.get('generation') if old.get('enabled') is on and same_endpoint else secrets.token_hex(12)}
    if endpoint:
        value['endpoint'] = endpoint
    elif endpoint is None and old.get('endpoint'):
        value['endpoint'] = old['endpoint']
    write_text_atomic(_path(sid), json.dumps(value))
    return ('SimpleX für diese Sitzung aktiviert. Der eingerichtete Messenger lädt dich in ihren Chat ein.'
            if on else 'SimpleX für diese Sitzung ausgeschaltet. Keine Nachrichten oder Handyfreigaben mehr.')


class LocalSession:
    """Dashboard adapter in the existing private per-user socket directory.

    Same peer check as the service: agent children cannot use this user door.
    The callback returns the same narrow pending/request shapes as the service.
    """
    def __init__(self, sid, callback):
        import socketserver
        import threading
        from delfin.agent.service import rundir
        from delfin.agent.service.supervisor import Supervisor, _peer_pid
        from types import SimpleNamespace
        self.sid = sid
        directory = rundir.ensure_runtime_dir()
        self.path = directory / ('msg-' + secrets.token_hex(8) + '.sock')
        parent = self
        class Handler(socketserver.StreamRequestHandler):
            def handle(self):
                if rundir.peer_uid(self.request) != os.getuid():
                    return
                identity = SimpleNamespace(sessions={}, starting=[SimpleNamespace(proc=SimpleNamespace(pid=os.getpid()))])
                if Supervisor._from_an_agent(identity, _peer_pid(self.request)):
                    return
                line = self.rfile.readline(16385)
                if len(line) > 16384:
                    return
                try:
                    command = json.loads(line)
                    if not enabled(sid) or command.get('sessionId') != sid:
                        reply = {'error': 'SimpleX ist für diese Sitzung ausgeschaltet.'}
                    else:
                        reply = callback(command)
                    self.wfile.write((json.dumps(reply) + '\n').encode())
                except Exception:
                    self.wfile.write(b'{"error":"Dashboard-Anfrage fehlgeschlagen."}\n')
        self.server = socketserver.ThreadingUnixStreamServer(str(self.path), Handler)
        self.server.daemon_threads = True
        os.chmod(self.path, 0o600)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def close(self):
        self.server.shutdown()
        self.server.server_close()
        self.path.unlink(missing_ok=True)


def ask_dashboard(sid, method, params):
    import socket
    path = record(sid).get('endpoint')
    if not path:
        raise ValueError('No dashboard endpoint')
    with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
        sock.settimeout(10)
        sock.connect(path)
        sock.sendall((json.dumps({'method': method, **params, 'sessionId': sid}) + '\n').encode())
        with sock.makefile('rb') as stream:
            data = stream.readline(65537)
        if len(data) > 65536:
            raise ValueError('Dashboard response too long')
        return json.loads(data)
