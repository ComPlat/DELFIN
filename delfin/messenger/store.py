"""The pairing state: which contact is the user's phone.

Stored 0600 in a 0700 folder under $HOME/.delfin/messenger/pair.json:
the SimpleX contact id of the one allowed sender and, after ``pair``,
nothing else (no address, no name -- those live in the daemon's own
database, 0700 as well). ``revoke`` removes the file: the phone may no
longer command the agent.
"""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path


def state_dir() -> Path:
    return Path(os.environ.get("DELFIN_MESSENGER_DIR", "")
                or (Path.home() / ".delfin" / "messenger"))


def pair_file() -> Path:
    return state_dir() / "pair.json"


def load() -> dict:
    """The pairing: {"contact_id": "..."}; {} when not paired."""
    path = pair_file()
    try:
        data = json.loads(path.read_text())
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def save(contact_id: str) -> None:
    with _locked():
        data = load()
        if data.get("contact_id") != str(contact_id):
            data = {"contact_id": str(contact_id)}
        _write(data)


def revoke() -> bool:
    """Unpair. True when something was there to remove."""
    with _locked():
        try:
            pair_file().unlink()
            return True
        except FileNotFoundError:
            return False


def paired_contact() -> str | None:
    return str(load().get("contact_id") or "") or None


def _ensure_private(d: Path) -> None:
    try:
        mode = stat.S_IMODE(d.stat().st_mode)
        if mode & 0o077:
            os.chmod(d, 0o700)
    except OSError:
        pass


def rooms() -> dict[str, dict]:
    data = load().get('rooms') or {}
    return data if isinstance(data, dict) else {}


def put_room(chat: str, session_id: str | None, title: str, *, owner: str,
             ready: bool = False) -> None:
    """Bind a group permanently; a rename cannot redirect its conversation."""
    with _locked():
        data = load()
        if data.get('contact_id') != owner:
            raise ValueError('pairing changed')
        bindings = data.setdefault('rooms', {})
        old = bindings.get(chat)
        if old and old['session_id'] != session_id:
            raise ValueError('a chat cannot be rebound to another session')
        if any(r['session_id'] == session_id and key != chat for key, r in bindings.items()):
            raise ValueError('session already has a chat')
        bindings[chat] = {'session_id': session_id, 'title': title, 'ready': ready}
        _write(data)


def forget_room(chat: str) -> None:
    with _locked():
        data = load()
        data.setdefault('rooms', {}).pop(chat, None)
        _write(data)


def claim_message(chat: str, message_id: str) -> bool:
    """Durable at-most-once dispatch, including across gateway restarts."""
    import sqlite3
    directory = state_dir()
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    with sqlite3.connect(directory / 'received.sqlite3') as db:
        db.execute('CREATE TABLE IF NOT EXISTS received (chat TEXT, id TEXT, PRIMARY KEY(chat,id))')
        return db.execute('INSERT OR IGNORE INTO received VALUES (?,?)', (chat, message_id)).rowcount == 1


from contextlib import contextmanager


@contextmanager
def _locked():
    import fcntl
    directory = state_dir()
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    _ensure_private(directory)
    fd = os.open(directory / 'state.lock', os.O_CREAT | os.O_RDWR, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        os.close(fd)


def _write(data: dict) -> None:
    import tempfile
    fd, name = tempfile.mkstemp(dir=state_dir(), prefix='.pair-')
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(data, stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, pair_file())
    finally:
        if os.path.exists(name):
            os.unlink(name)
