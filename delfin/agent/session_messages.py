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
import uuid
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
# How long an inbox file alone proves an address. One message to a typo (e.g.
# ``operatro``) writes an inbox that nobody ever takes; without a bound the
# key latches 'known' FOREVER, so every later message to the typo is silently
# accepted and dropped -- the waves 11/12 "sender never knows if read" bug
# re-created on the typo path (reviewer nacht-s17). QS decision (s26): an
# inbox OLDER than session_presence._STALE_S no longer marks an address
# deliverable; a real closed session that re-announces presence is open again,
# so legitimate delivery is preserved, while a stale typo reverts to being
# refused loudly.
# Delivery receipts live in one sidecar under _DIR, so ``status(message_id)``
# can still answer after take() unlinks the inbox. .receipts is not a *.jsonl
# inbox and is never scanned as one.
_RECEIPTS = ".receipts.json"


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
    if _inbox_known(key):
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
        return _inbox_known(key)
    return False


def _inbox_known(key: str) -> bool:
    """Whether an inbox alone proves ``key``: only while it is fresh.

    The inbox file's mtime refreshes on every append, so a live queue keeps
    its recipient known; a typo with one orphan message ages out (older than
    session_presence._STALE_S, the same window a dormant presence is reaped
    at) and a later message to it is refused instead of silently accepted
    forever."""
    inbox = _inbox(key)
    try:
        from . import session_presence as _presence
        window = _presence._STALE_S
    except Exception:
        window = 15 * 60.0
    try:
        return time.time() - inbox.stat().st_mtime < window
    except OSError:
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


def _reserved_canonical(key: str) -> str:
    """The one mailbox every case-variant of a reserved address shares.

    deliverable() is case-insensitive for a reserved mailbox
    (deliverable('OPERATOR') is True, lower() against _RESERVED), so a
    case-variant send must write the SAME file a read by its canonical key
    opens. Without this, send('OPERATOR') writes OPERATOR.jsonl while
    take('operator') reads operator.jsonl and the message is dropped."""
    low = str(key).lower()
    if low in _RESERVED:
        return low
    return str(key)


def _inbox(key: str) -> Path:
    key = _reserved_canonical(key)
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
        "id": uuid.uuid4().hex,
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
    # Delivery receipts: the messages handed to the receiver are now
    # delivered. Persist them (under the receipts lock) so status(message_id)
    # still answers after the inbox below was unlinked.
    _receipts_deliver([m for m in out if m.get("id")])
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


def _receipts_path() -> Path:
    return _DIR / _RECEIPTS


def _receipts() -> dict:
    try:
        data = json.loads(_receipts_path().read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _receipts_save(receipts: dict) -> None:
    from .state_paths import ensure_dir
    ensure_dir(_DIR)
    path = _receipts_path()
    from .bash_jobs import cross_process_lock
    with cross_process_lock(path):
        try:
            path.write_text(json.dumps(receipts, ensure_ascii=False),
                            encoding="utf-8")
        except OSError:
            pass


def _receipts_deliver(messages: list[dict]) -> None:
    """Record that each ``message`` was handed to its recipient (delivered)."""
    if not messages:
        return
    receipts = _receipts()
    now = time.time()
    for m in messages:
        mid = str(m.get("id") or "")
        if not mid:
            continue
        rec = receipts.get(mid) or {}
        rec.setdefault("to", m.get("to"))
        rec.setdefault("from", m.get("from"))
        rec.setdefault("from_title", m.get("from_title"))
        rec.setdefault("text", m.get("text"))
        rec.setdefault("sent_at", m.get("sent_at"))
        rec["taken_at"] = now
        rec["read_at"] = rec.get("read_at")
        receipts[mid] = rec
    _receipts_save(receipts)


def mark_read(message_id: str) -> None:
    """Mark a delivered message as read (the recipient turned to it)."""
    mid = str(message_id or "")
    if not mid:
        return
    receipts = _receipts()
    if mid not in receipts:
        return
    receipts[mid]["read_at"] = time.time()
    _receipts_save(receipts)


def _queued_inbox_of(message_id: str) -> "str | None":
    """The inbox key whose list holds ``message_id``, or None if it is not
    waiting anywhere (it was delivered/read or never existed)."""
    mid = str(message_id or "")
    try:
        inboxes = [p for p in _DIR.glob("*.jsonl") if p.is_file()]
    except OSError:
        inboxes = []
    for inbox in inboxes:
        try:
            for raw in inbox.read_text(encoding="utf-8",
                                       errors="replace").splitlines():
                line = json.loads(raw)
                if isinstance(line, dict) and line.get("id") == mid:
                    return inbox.stem
        except Exception:
            continue
    return None


def status(message_id: str) -> str:
    """What happened to a sent message: queued / delivered / read / unknown.

    ``queued`` - still in its recipient's inbox. ``delivered`` - taken by the
    recipient (receipt exists, not yet read). ``read`` - the recipient marked
    it read. ``unknown`` - no such message id exists."""
    mid = str(message_id or "")
    if not mid:
        return "unknown"
    rec = _receipts().get(mid)
    if rec is not None:
        if rec.get("read_at"):
            return "read"
        if rec.get("taken_at"):
            return "delivered"
    if _queued_inbox_of(mid) is not None:
        return "queued"
    return "unknown"


def ls() -> list[dict]:
    """Every known message (queued or with a receipt), for a `messages ls`
    listing. One row per message id: {id, to, from, status, sent_at, text}."""
    rows: dict = {}
    try:
        inboxes = [p for p in _DIR.glob("*.jsonl") if p.is_file()]
    except OSError:
        inboxes = []
    for inbox in inboxes:
        try:
            raw = inbox.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for line in raw.splitlines():
            try:
                m = json.loads(line)
            except ValueError:
                continue
            if isinstance(m, dict) and m.get("id"):
                rows[m["id"]] = {
                    "id": m["id"], "to": m.get("to"), "from": m.get("from"),
                    "status": "queued", "sent_at": m.get("sent_at"),
                    "text": m.get("text"),
                }
    for mid, rec in _receipts().items():
        st = "read" if rec.get("read_at") else "delivered"
        rows[mid] = {
            "id": mid, "to": rec.get("to"), "from": rec.get("from"),
            "status": st, "sent_at": rec.get("sent_at"),
            "text": rec.get("text"),
        }
    return sorted(rows.values(), key=lambda r: r.get("sent_at") or 0.0)
