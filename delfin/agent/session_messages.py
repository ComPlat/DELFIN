"""Messages between agent sessions.

A session can write to another open session: hand over a finding, ask it to
leave a file alone, say it has pushed. A message waits in the receiver's
inbox until that session takes it -- between the rounds of a turn that is
running, or as the next turn of one that is idle -- and always reads as
coming from another session, never from the user.

Inboxes are keyed like session_presence records, by the session's place in
its session list.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

_DIR = Path.home() / ".delfin" / "session_inbox"
_MAX_TEXT = 4000
# At most this many messages reach one prompt; older ones are counted, not read.
_MAX_TAKE = 20


def _inbox(key: str) -> Path:
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in key)[:80]
    return _DIR / f"{safe}.jsonl"


def _lock(fd: int) -> None:
    try:
        import fcntl
        fcntl.flock(fd, fcntl.LOCK_EX)
    except (ImportError, OSError):
        pass            # a filesystem without flock: best effort


def send(to_key: str, text: str, *, from_key: str = "",
         from_title: str = "") -> dict:
    """Leave ``text`` in session ``to_key``'s inbox. Returns the message."""
    message = {
        "to": str(to_key), "from": str(from_key or ""),
        "from_title": str(from_title or "")[:80],
        "text": str(text or "")[:_MAX_TEXT], "sent_at": time.time(),
    }
    line = (json.dumps(message, ensure_ascii=False) + "\n").encode("utf-8")
    from .state_paths import ensure_dir
    ensure_dir(_DIR)
    path = _inbox(str(to_key))
    for _attempt in range(5):
        fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o600)
        try:
            _lock(fd)
            # take() unlinks the inbox under the same lock. A sender that
            # opened it just before then holds a file nobody will read
            # again, so it starts over with the fresh one.
            try:
                current = os.stat(path).st_ino == os.fstat(fd).st_ino
            except FileNotFoundError:
                current = False
            if not current:
                continue
            os.write(fd, line)
            return message
        finally:
            os.close(fd)
    raise OSError(f"could not write to the inbox of {to_key}")


def take(key: str) -> list[dict]:
    """Every message waiting for session ``key``, removed from its inbox."""
    path = _inbox(str(key))
    try:
        fd = os.open(path, os.O_RDONLY)
    except FileNotFoundError:
        return []
    try:
        _lock(fd)
        try:
            os.unlink(path)
        except FileNotFoundError:
            return []
        chunks = []
        while True:
            chunk = os.read(fd, 65536)
            if not chunk:
                break
            chunks.append(chunk)
    finally:
        os.close(fd)
    out: list[dict] = []
    for raw in b"".join(chunks).decode("utf-8", "replace").splitlines():
        try:
            message = json.loads(raw)
        except ValueError:
            continue
        if isinstance(message, dict):
            out.append(message)
    # Everything waiting went into ONE prompt, uncapped: a sender in a loop
    # could fill the receiver's whole context with one take. The newest
    # messages are delivered; the receiver is told how many older ones
    # were not.
    if len(out) > _MAX_TAKE:
        dropped = len(out) - _MAX_TAKE
        out = out[-_MAX_TAKE:]
        out.insert(0, {"from": "", "from_title": "the session inbox",
                       "text": f"{dropped} earlier message(s) were not "
                               f"delivered: the inbox held more than "
                               f"{_MAX_TAKE}.", "sent_at": time.time()})
    return out


def render(message: dict) -> str:
    """A message as the receiving agent reads it."""
    sender = str(message.get("from") or "")
    who = str(message.get("from_title") or sender or "another session")
    reply = (f' Answer with session_message(to="{sender}", message=...) if '
             "it needs one." if sender else "")
    return (f'[Message from the session "{who}" — not from the user.{reply}]\n'
            f"{message.get('text') or ''}")
