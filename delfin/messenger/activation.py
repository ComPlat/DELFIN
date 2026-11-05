"""Session opt-in owned by the DELFIN user, never by the phone gateway.

One record per session under ``~/.delfin/messenger_sessions`` (masked from
agent file tools and shells). ``directory`` names another person's records
(people mode: ``~/people/<name>/.delfin/messenger_sessions``); the bridge
passes it, everything run in a person's own session uses the default.

A record says whether the session is enabled, who switched it last
(``by``: ``session`` -- typed in the session -- or ``phone`` -- ``hol`` /
``weg`` on the person's verified phone) and the phone's label. A session
switched off in the session itself cannot be switched on from the phone.
The session's display name is not kept here: it is the agent service's
session title (session.setTitle), the one field every window and the phone
show; ``/sync NAME`` sets it where the session runs in the service.
"""
from __future__ import annotations

import hashlib
import json
import os
import secrets
from pathlib import Path

_DIR = Path.home() / '.delfin' / 'messenger_sessions'


def _base(directory=None):
    return Path(directory) if directory is not None else _DIR


def _path(sid, directory=None):
    return _base(directory) / (hashlib.sha256(sid.encode()).hexdigest() + '.json')


def record(sid, *, directory=None):
    try:
        value = json.loads(_path(sid, directory).read_text())
        return value if value.get('session_id') == sid else {}
    except (OSError, ValueError, AttributeError):
        return {}


def enabled(sid, *, directory=None):
    return record(sid, directory=directory).get('enabled') is True


def switched_off_in_session(sid, *, directory=None):
    """True when the session itself said /sync off: the phone may not undo it."""
    value = record(sid, directory=directory)
    return value.get('enabled') is False and value.get('by', 'session') != 'phone'


def _policy_path(directory=None):
    return _base(directory) / 'policy.json'


def pull_allowed(*, directory=None):
    """Whether this person's verified phones may pull a running session
    (``hol``): the person's own switch, ``/sync policy pull on|off`` at a
    window. Default on; an unreadable setting counts as off."""
    path = _policy_path(directory)
    if not path.exists():
        return True
    try:
        return json.loads(path.read_text()).get('pull') is True
    except (OSError, ValueError, AttributeError):
        return False


def set_pull(on, *, directory=None):
    from delfin.agent.state_paths import write_text_atomic
    write_text_atomic(_policy_path(directory), json.dumps({'pull': bool(on)}))


def records(*, directory=None):
    result = []
    for path in _base(directory).glob('*.json'):
        try:
            value = json.loads(path.read_text())
            if value.get('enabled') is True and value.get('session_id'):
                result.append(value)
        except (OSError, ValueError, AttributeError):
            continue
    return result


#: The longest session display name a phone chat shows.
NAME_MAX = 50


def clean_name(text):
    """A display name: one line, no control characters, trimmed, at most
    NAME_MAX characters ("" means the default: folder name and short id)."""
    try:
        from delfin.agent.service.protocol import clean_title
        return clean_title(text)[:NAME_MAX].strip()
    except Exception:  # noqa: BLE001 - a fallback with the same effect
        return ' '.join(''.join(c if c.isprintable() else ' ' for c in str(text or '')).split())[:NAME_MAX].strip()


#: What /sync understands (German matchers for what the user types).
SYNC_ON = frozenset({'', 'on', 'an', 'ein', 'enable'})
SYNC_OFF = frozenset({'off', 'aus', 'ausschalten', 'disable', 'stop'})
SYNC_USAGE = ('Usage: /sync [on|off|status] or /sync name NAME; '
              '/sync policy pull on|off for pulling from the phone. Anything else changes nothing.')


def parse_sync(args):
    """``/sync`` arguments as (action, name): action on | off | status |
    usage; a name only after ``name`` or in quotes. Anything else is
    ``usage`` -- it changes nothing (an unknown word never switches on)."""
    raw = str(args or '').strip()
    word = raw.lower()
    if word in SYNC_ON:
        return 'on', ''
    if word in SYNC_OFF:
        return 'off', ''
    if word == 'status':
        return 'status', ''
    head, _, rest = raw.partition(' ')
    if head.lower() == 'policy':
        setting = ' '.join(rest.lower().split())
        return ('policy', setting) if setting in {'pull on', 'pull off'} else ('usage', '')
    if head.lower() == 'name' and clean_name(rest):
        return 'on', clean_name(rest)
    if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in '"\'' and clean_name(raw[1:-1]):
        return 'on', clean_name(raw[1:-1])
    return 'usage', ''


def configure(sid, args='', *, endpoint=None, directory=None, by='session', phone=''):
    """``/sync`` [on | off | status | NAME]: switch or report. ``NAME``
    switches on like ``on``; naming the session is the caller's (the
    service worker sets its title). ``by``/``phone``: who switched (the
    bridge passes ``phone`` and the phone's label for ``hol`` and ``weg``)."""
    if not sid:
        return 'No session yet. Start a session first.'
    args, _name = parse_sync(args)
    if args == 'usage':
        return SYNC_USAGE
    if args == 'policy':
        on = _name == 'pull on'
        set_pull(on, directory=directory)
        return 'Pulling sessions from the phone (get): ' + ('on.' if on else 'off. /sync policy pull on allows it '
                                                                             'again.')
    if args == 'status':
        old = record(sid, directory=directory)
        text = 'SimpleX for this session: ' + ('on' if old.get('enabled') is True else 'off')
        if old.get('by') == 'phone' and old.get('phone'):
            text += ' (from the phone ' + str(old['phone'])[:40] + ')'
        return text + '.'
    from delfin.agent.state_paths import write_text_atomic
    old = record(sid, directory=directory)
    on = args != 'off'
    same_endpoint = endpoint is None or (old.get('endpoint') or '') == endpoint
    value = {'session_id': sid, 'enabled': on,
             'generation': old.get('generation') if old.get('enabled') is on and same_endpoint else secrets.token_hex(12),
             'by': 'phone' if by == 'phone' else 'session'}
    if by == 'phone' and phone:
        value['phone'] = str(phone)[:40]
    if endpoint:
        value['endpoint'] = endpoint
    elif endpoint is None and old.get('endpoint'):
        value['endpoint'] = old['endpoint']
    write_text_atomic(_path(sid, directory), json.dumps(value))
    if not on:
        return 'SimpleX switched off for this session. No more messages or phone approvals.'
    return 'SimpleX switched on for this session. The messenger invites you to its chat.'


#: A session's deliberate notes to the phone (push_notification): at most
#: one per NOTE_EVERY seconds, the last NOTE_KEEP kept.
NOTE_EVERY = 60.0
NOTE_KEEP = 5


def _notes_path(sid, directory=None):
    return _base(directory) / (hashlib.sha256(sid.encode()).hexdigest() + '.notes')


def notes(sid, *, directory=None):
    """The session's recent notes to the phone, oldest first:
    ``[{"id", "ts", "title", "body"}]``."""
    try:
        value = json.loads(_notes_path(sid, directory).read_text())
        return [n for n in value if isinstance(n, dict)] if isinstance(value, list) else []
    except (OSError, ValueError):
        return []


def post_note(sid, title, body='', *, directory=None, now=None):
    """Queue a note for the session's phone (the gateway sends it with its
    next look). Output: True, ``'rate_limited'`` (one per NOTE_EVERY s),
    or False when the session is not on the phone."""
    import time
    if not sid or not enabled(sid, directory=directory):
        return False
    now = time.time() if now is None else now
    kept = notes(sid, directory=directory)
    if kept and now - float(kept[-1].get('ts') or 0) < NOTE_EVERY:
        return 'rate_limited'
    from delfin.agent.state_paths import write_text_atomic
    kept.append({'id': secrets.token_hex(6), 'ts': now,
                 'title': ' '.join(str(title or '').split())[:80],
                 'body': str(body or '').strip()[:300]})
    write_text_atomic(_notes_path(sid, directory), json.dumps(kept[-NOTE_KEEP:]))
    return True


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
                        reply = {'error': 'SimpleX is switched off for this session.'}
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


def ask_dashboard(sid, method, params, *, directory=None, within=None):
    """One request to a live dashboard tab's adapter. ``within``: the
    runtime directory the endpoint must lie in (a person's own); an
    endpoint elsewhere is refused, never connected."""
    import socket
    path = record(sid, directory=directory).get('endpoint')
    if not path:
        raise ValueError('No dashboard endpoint')
    if within is not None:
        try:
            Path(path).resolve().relative_to(Path(within).resolve())
        except (OSError, ValueError):
            raise ValueError('Dashboard endpoint outside the person\'s runtime directory') from None
    with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
        sock.settimeout(10)
        sock.connect(path)
        sock.sendall((json.dumps({'method': method, **params, 'sessionId': sid}) + '\n').encode())
        with sock.makefile('rb') as stream:
            data = stream.readline(65537)
        if len(data) > 65536:
            raise ValueError('Dashboard response too long')
        return json.loads(data)
