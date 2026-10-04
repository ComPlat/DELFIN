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
import re
import os
import time
from pathlib import Path

_DIR = Path.home() / ".delfin" / "session_inbox"
_MAX_TEXT = 4000
# At most this many messages reach one prompt; older ones are counted, not read.
_MAX_TAKE = 20
# The operator is always reachable: a reserved mailbox, no presence needed. A
# message to it is queued until a session opens under this key (waves 11/12:
# it was refused with "no other open session 'operator'" unless a heartbeat
# process faked presence).
_RESERVED = ("operator",)


def _known_key(key: str) -> bool:
    """Whether ``key`` names a session that exists (or has existed).

    A presence record that is stale still names a real session -- one that has
    closed or is mid-restart -- and a recipient that has an inbox has received
    before. Either marks the address as known, so a message to a closed
    session queues for its next start instead of being refused.
    Order matters: a presence record may be reaped (deleted) by a check of
    open_sessions (its `_reap` drops same-host stale records), so look for the
    key's record file before asking whether it is open."""
    if key == "operator":
        return True
    if _inbox(key).exists():
        return True
    try:
        from . import session_presence as P
        # A record file (alive or stale) names a real session. Look before
        # open_sessions() runs `_reap`, which deletes same-host stale records.
        record = P._path(key)
        if record.exists():
            try:
                data = json.loads(record.read_text(encoding="utf-8"))
            except Exception:
                data = {}
            if isinstance(data, dict) and data.get("key") == key:
                return True
        for record in P.open_sessions():
            if record.get("key") == key:
                return True
    except Exception:
        return _inbox(key).exists()
    return False


def deliverable(to_key: str) -> bool:
    """Can a message be left for session ``to_key`` -- queued until it takes
    it, rather than refused as "no other open session"?

    The operator is always deliverable (a reserved mailbox). A session that is
    known -- announced, stale-but-recorded, or already holding an inbox -- is
    deliverable even while closed: the message waits and is picked up on its
    next start. An address that is neither reserved nor known is a typo and is
    refused."""
    to_key = str(to_key or "").strip()
    if not to_key:
        return False
    if to_key.lower() in _RESERVED:
        return True
    return _known_key(to_key)


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
    """Leave ``text`` in session ``to_key``'s inbox. Returns the message.

    Credentials are redacted first: the other session may run on another
    provider, and the text goes into its transcript and every later
    request it makes."""
    try:
        from .output_guard import scrub_secrets
        text = scrub_secrets(str(text or ""))
        from_title = scrub_secrets(str(from_title or ""))
    except Exception:
        pass
    message = {
        "to": str(to_key), "from": str(from_key or ""),
        "from_title": str(from_title or "")[:80],
        "text": str(text or "")[:_MAX_TEXT], "sent_at": time.time(),
    }
    line = (json.dumps(message, ensure_ascii=False) + "\n").encode("utf-8")
    from .state_paths import ensure_dir
    ensure_dir(_DIR)
    path = _inbox(str(to_key))
    from .bash_jobs import cross_process_lock
    # flock on the inbox is node-local on some shared filesystems; the
    # cross-process lock adds a lease that holds across login nodes.
    with cross_process_lock(path):
        return _append(path, line, message, to_key)


def _append(path, line: bytes, message: dict, to_key: str) -> dict:
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
    from .bash_jobs import cross_process_lock
    with cross_process_lock(path):
        return _take(path)


def _take(path) -> list[dict]:
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
    sender = re.sub(r"[^A-Za-z0-9_-]", "", str(message.get("from") or ""))[:40]
    # The title is another session's first message: it must not be able to
    # close the header's quotes or bracket and write text that reads as the
    # harness (security review 2026-09-16).
    who = re.sub(r'[\[\]"\n\r]', " ", str(message.get("from_title") or sender
                                          or "another session"))[:80].strip()
    # "if it needs one" read as an invitation: two sessions greeted each other
    # in a loop, each turn answering the last acknowledgement (driven
    # 2026-09-16). A reply is for a question or a request, never for thanks.
    reply = (f' Reply with session_message(to="{sender}", message=...) only if '
             "it asks you for something; a greeting, thanks or an "
             "acknowledgement gets no reply." if sender else "")
    return (f'[Message from the session "{who}" — not from the user.{reply}]\n'
            f"{message.get('text') or ''}")
