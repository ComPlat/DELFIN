"""The gateway's pairing state: which contact is whose phone.

Stored 0600 in a 0700 folder ($HOME/.delfin/messenger or
DELFIN_MESSENGER_DIR), ``pair.json``, format 2::

    {"version": 2,
     "contacts": {"<contact id>": {"person": "", "state": "verified",
                                   "label": "...", "since": ts, ...}},
     "rooms": {"#<group>": {"session_id", "title", "ready", "person", ...}},
     "invitations": {"<conn id>": {"person", "request", "by", "expires", "tag"}},
     "retired": {"<conn id>": ts}, "issued": {"<person>": [ts, ...]},
     "revoked": [{"contact", "person", "label", "at", "by", "cleaned"}],
     "handled": {"<person>": {"<request id>": ts}},
     "schedules_room", "schedule_seq",
     "mac": "<HMAC-SHA256 of the rest>"}

``person`` is a person of people mode (delfin.messenger.people) or ``""``,
the account owner. A contact is bound to exactly one person; the binding is
never changed in place: :func:`bind` refuses a contact that is bound, and
a rebind is :func:`revoke_contact` plus a new pairing. A contact is
``unverified`` from the moment its phone connects until the person types the
code the phone received (:func:`check_code`); only ``verified`` contacts may
command anything.

Integrity: the file carries an HMAC under ``store.key`` (32 random bytes,
0600, beside it, created with the first write); the first write also
leaves a ``format`` marker. A file whose MAC does not match, a key file
that is unreadable or not 32 bytes, a format-1 file once the key or the
marker exists (a downgrade, also with the key deleted), a file that is not
JSON or a contact bound to a person name that is not valid refuses the
whole store: :func:`refused` is true, every reader sees no contact, every
writer raises :class:`StoreRefused`, and one ``messenger_store_refused``
audit event is recorded per process.

What this does and does not hold: it detects an edit, a damaged key and a
downgrade by anything that leaves the key or the marker in place. A process
of the gateway's own user can read the key, rewrite a valid MAC, or delete
key and marker together and write a format-1 file; that is not detected.
It is a consistency check, not a boundary against that user (the
gateway's state directory is 0700 of ``delfin-gw``).

Format 1 (``{"contact_id", "rooms", ...}``, one phone per installation) is
read and migrated on the first write without re-pairing: its contact becomes
the owner's verified phone, its rooms the owner's rooms.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
import stat
import threading
import time
from contextlib import contextmanager
from pathlib import Path

from . import pairing
from .people import OWNER, valid as valid_person

VERSION = 2
_REFUSED_SEEN: set[str] = set()
_REFUSED_LOCK = threading.Lock()


class StoreRefused(RuntimeError):
    """The pairing store failed its integrity or format check."""


class PairingRefused(ValueError):
    """A pairing step refused: ``reason`` is limit | rate | bound | person."""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


def state_dir() -> Path:
    return Path(os.environ.get("DELFIN_MESSENGER_DIR", "")
                or (Path.home() / ".delfin" / "messenger"))


def pair_file() -> Path:
    return state_dir() / "pair.json"


def key_file() -> Path:
    return state_dir() / "store.key"


def format_file() -> Path:
    """Written with the first format-2 write: from then on a format-1 file
    is a downgrade, whether or not the key is still there."""
    return state_dir() / "format"


# -- reading -------------------------------------------------------------------------

def _key(create: bool = False) -> bytes | None:
    path = key_file()
    try:
        fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        with os.fdopen(fd, "rb") as stream:
            data = stream.read(65)
        return data if len(data) == 32 else None
    except FileNotFoundError:
        if not create:
            return None
    except OSError:
        return None
    key = secrets.token_bytes(32)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o600)
    with os.fdopen(fd, "wb") as stream:
        stream.write(key)
        stream.flush()
        os.fsync(stream.fileno())
    return key


def _mac(key: bytes, data: dict) -> str:
    body = json.dumps({k: v for k, v in data.items() if k != "mac"}, sort_keys=True,
                      ensure_ascii=False).encode()
    return hmac.new(key, body, hashlib.sha256).hexdigest()


def _migrate(old: dict) -> dict:
    data = _empty()
    cid = str(old.get("contact_id") or "")
    if cid:
        data["contacts"][cid] = {"person": OWNER, "state": "verified", "label": "",
                                 "since": time.time(), "verified": time.time(),
                                 "request": "format-1", "legacy": True}
    for chat, room in (old.get("rooms") or {}).items():
        if isinstance(room, dict):
            data["rooms"][chat] = {**room, "person": OWNER, "members": [cid] if cid else []}
    for key in ("schedules_room", "schedule_seq"):
        if key in old:
            data[key] = old[key]
    return data


def _empty() -> dict:
    return {"version": VERSION, "contacts": {}, "rooms": {}, "invitations": {},
            "retired": {}, "issued": {}, "revoked": [], "handled": {}}


def _refuse(why: str) -> None:
    with _REFUSED_LOCK:
        first = why not in _REFUSED_SEEN
        _REFUSED_SEEN.add(why)
    if first:
        from . import audit
        audit.event("messenger_store_refused", f"pairing store refused: {why}", blocked=True)
    raise StoreRefused(why)


def _load_checked() -> dict:
    """The store as written, checked; raises StoreRefused."""
    path = pair_file()
    try:
        raw = path.read_bytes()
    except FileNotFoundError:
        return _empty()
    except OSError:
        _refuse("unreadable")
    try:
        data = json.loads(raw)
    except ValueError:
        _refuse("not JSON")
    if not isinstance(data, dict):
        _refuse("not an object")
    key = _key()
    if data.get("version") is None:
        if os.path.lexists(key_file()) or os.path.lexists(format_file()):
            _refuse("format 1 after format 2 (downgrade)")
        return _migrate(data)
    if data.get("version") != VERSION:
        _refuse("unknown format")
    if key is None or not hmac.compare_digest(str(data.get("mac") or ""), _mac(key, data)):
        _refuse("integrity check failed")
    contacts = data.get("contacts")
    if not isinstance(contacts, dict):
        _refuse("no contact table")
    for cid, rec in contacts.items():
        if (not isinstance(rec, dict) or not str(cid).isdigit()
                or not valid_person(rec.get("person"))
                or rec.get("state") not in {"verified", "unverified"}):
            _refuse("malformed contact")
    for chat, room in (data.get("rooms") or {}).items():
        if not isinstance(room, dict) or not valid_person(room.get("person", OWNER)):
            _refuse("malformed room")
    base = _empty()
    base.update({k: v for k, v in data.items() if k != "mac"})
    return base


def refused() -> bool:
    try:
        _load_checked()
        return False
    except StoreRefused:
        return True


def load() -> dict:
    """The store, or ``{}`` when it is refused (nobody is paired then)."""
    try:
        return _load_checked()
    except StoreRefused:
        return {}


def contacts() -> dict:
    return dict(load().get("contacts") or {})


def contact(cid: str) -> dict | None:
    rec = contacts().get(str(cid))
    return dict(rec) if isinstance(rec, dict) else None


def verified_person(cid: str) -> str | None:
    """The person of a verified contact; None for anyone else."""
    rec = contact(cid)
    return rec["person"] if rec and rec.get("state") == "verified" else None


def phones(person: str, *, verified_only: bool = False) -> list[tuple[str, dict]]:
    return [(cid, rec) for cid, rec in sorted(contacts().items())
            if rec.get("person") == person
            and (not verified_only or rec.get("state") == "verified")]


def persons() -> list[str]:
    """Every person with a verified phone, the owner first."""
    found = {rec["person"] for rec in contacts().values() if rec.get("state") == "verified"}
    return sorted(found, key=lambda p: (p != OWNER, p))


def paired_contact() -> str | None:
    """The owner's first verified phone (format 1 knew only this one)."""
    for cid, rec in phones(OWNER, verified_only=True):
        return cid
    return None


# -- writing -------------------------------------------------------------------------

def _ensure_private(d: Path) -> None:
    try:
        mode = stat.S_IMODE(d.stat().st_mode)
        if mode & 0o077:
            os.chmod(d, 0o700)
    except OSError:
        pass


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


@contextmanager
def _editing():
    """The store under the lock, written back when the block ends."""
    with _locked():
        data = _load_checked()
        yield data
        _write(data)


def _write(data: dict) -> None:
    import tempfile
    key = _key(create=True)
    if key is None:
        _refuse("key file unreadable or of the wrong length")
    if not os.path.lexists(format_file()):
        fd = os.open(format_file(), os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o600)
        with os.fdopen(fd, "w") as stream:
            stream.write("2\n")
    data = {**{k: v for k, v in data.items() if k != "mac"}, "version": VERSION}
    data["mac"] = _mac(key, data)
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


def save(contact_id: str, *, person: str = OWNER, label: str = "") -> None:
    """Bind ``contact_id`` to ``person`` as a VERIFIED phone. For the
    migration and for tests; a pairing goes through :func:`bind` and
    :func:`check_code`."""
    cid = str(contact_id)
    with _editing() as data:
        old = data["contacts"].get(cid)
        if old and old.get("person") != person:
            raise PairingRefused("bound")
        data["contacts"][cid] = {"person": person, "state": "verified", "label": label,
                                 "since": time.time(), "verified": time.time(),
                                 "request": "direct"}


def revoke() -> bool:
    """Unpair everybody (format 1's revoke). True when something was there."""
    with _locked():
        try:
            pair_file().unlink()
            return True
        except FileNotFoundError:
            return False


def _bump_issued(data: dict, person: str, now: float) -> None:
    times = [t for t in data["issued"].get(person, []) if now - t < pairing.INVITE_WINDOW_S]
    if len(times) >= pairing.INVITES_PER_WINDOW:
        raise PairingRefused("rate")
    data["issued"][person] = [*times, now]


def add_invitation(conn_id: str, *, person: str, request: str, by: str = "gateway",
                   tag: str = "", now: float | None = None) -> dict:
    """Record a new one-time invitation for ``person``; raises PairingRefused
    (limit: MAX_PHONES phones, rate: INVITES_PER_WINDOW in INVITE_WINDOW_S)."""
    now = time.time() if now is None else now
    if not valid_person(person) or not conn_id:
        raise PairingRefused("person")
    with _editing() as data:
        owned = [c for c, r in data["contacts"].items() if r.get("person") == person]
        if len(owned) >= pairing.MAX_PHONES:
            raise PairingRefused("limit")
        _bump_issued(data, person, now)
        entry = {"person": person, "request": str(request), "by": by, "created": now,
                 "expires": now + pairing.INVITATION_TTL_S, "tag": tag}
        data["invitations"][str(conn_id)] = entry
        return dict(entry)


def invitations() -> dict:
    return dict(load().get("invitations") or {})


def open_invitation(person: str, request: str) -> tuple[str, dict] | None:
    for conn, inv in invitations().items():
        if inv.get("person") == person and inv.get("request") == request:
            return conn, inv
    return None


def consume_invitation(conn_id: str, *, by: str = "gateway",
                       now: float | None = None) -> tuple[dict | None, str]:
    """Use an invitation once: (invitation, "ok"), or (None, why) with why
    ``replayed`` (used or expired before), ``expired``, ``other`` (another
    process's invitation) or ``unknown``."""
    now = time.time() if now is None else now
    conn_id = str(conn_id or "")
    with _editing() as data:
        if conn_id in data["retired"]:
            used = data["retired"][conn_id]
            # Used by the other process (the admin's pairing beside the
            # gateway): theirs to handle, not a replay.
            if isinstance(used, dict) and used.get("by") != by:
                return None, "other"
            # expired: dropped by the clock unused; anything else was used.
            return None, ("expired" if isinstance(used, dict) and used.get("why") == "expired"
                          else "replayed")
        inv = data["invitations"].get(conn_id)
        if inv is None:
            return None, "unknown"
        if inv.get("by") != by:
            return None, "other"
        del data["invitations"][conn_id]
        data["retired"][conn_id] = {"at": now, "by": by}
        if now >= float(inv.get("expires") or 0):
            return None, "expired"
        return dict(inv), "ok"


def cancel_invitations(person: str, *, request: str | None = None,
                       now: float | None = None) -> list[tuple[str, dict]]:
    now = time.time() if now is None else now
    gone = []
    with _editing() as data:
        for conn, inv in list(data["invitations"].items()):
            if inv.get("person") == person and (request is None or inv.get("request") == request):
                gone.append((conn, dict(inv)))
                del data["invitations"][conn]
                data["retired"][conn] = {"at": now, "by": inv.get("by", "gateway")}
    return gone


def bind(cid: str, *, person: str, label: str, request: str, code: str,
         now: float | None = None) -> dict:
    """Bind a contact that connected through ``person``'s invitation, as
    UNVERIFIED, with the code its phone is sent. A contact bound already
    (to anyone) is refused: a binding never changes in place."""
    now = time.time() if now is None else now
    cid = str(cid)
    salt = secrets.token_hex(8)
    with _editing() as data:
        if cid in data["contacts"]:
            raise PairingRefused("bound")
        owned = [c for c, r in data["contacts"].items() if r.get("person") == person]
        if len(owned) >= pairing.MAX_PHONES:
            raise PairingRefused("limit")
        rec = {"person": person, "state": "unverified", "label": str(label)[:40],
               "since": now, "verified": None, "request": str(request),
               "code": salt + ":" + hashlib.sha256((salt + pairing.digits(code)).encode()).hexdigest(),
               "attempts": 0}
        data["contacts"][cid] = rec
        return dict(rec)


def check_code(cid: str, code: str, *, request: str | None = None,
               now: float | None = None) -> tuple[str, dict | None]:
    """(``verified`` | ``wrong`` | ``failed`` | ``expired`` | ``none``, record).

    ``none``: no unverified contact (of that request). ``failed``: the
    MAX_ATTEMPTS-th wrong code; the contact is revoked. ``expired``: after
    VERIFY_TTL_S; the contact is revoked."""
    now = time.time() if now is None else now
    with _editing() as data:
        cid = str(cid)
        rec = data["contacts"].get(cid)
        if (rec is None or rec.get("state") != "unverified"
                or (request is not None and rec.get("request") != request)):
            return "none", None
        if now - float(rec.get("since") or 0) > pairing.VERIFY_TTL_S:
            _revoke_in(data, cid, by="expired", now=now)
            return "expired", dict(rec)
        salt, _, want = str(rec.get("code") or "").partition(":")
        got = hashlib.sha256((salt + pairing.digits(code)).encode()).hexdigest()
        if want and hmac.compare_digest(want, got):
            rec.update(state="verified", verified=now, attempts=0)
            rec.pop("code", None)
            return "verified", dict(rec)
        rec["attempts"] = int(rec.get("attempts") or 0) + 1
        if rec["attempts"] >= pairing.MAX_ATTEMPTS:
            _revoke_in(data, cid, by="wrong code", now=now)
            return "failed", dict(rec)
        return "wrong", dict(rec)


def unverified_of(person: str, request: str) -> str | None:
    for cid, rec in contacts().items():
        if rec.get("person") == person and rec.get("request") == request \
                and rec.get("state") == "unverified":
            return cid
    return None


def _revoke_in(data: dict, cid: str, *, by: str, now: float) -> dict | None:
    rec = data["contacts"].pop(str(cid), None)
    if rec is None:
        return None
    entry = {"contact": str(cid), "person": rec.get("person", OWNER),
             "label": rec.get("label", ""), "state": rec.get("state"),
             "at": now, "by": by, "cleaned": False}
    data["revoked"] = [*data.get("revoked", []), entry][-100:]
    for room in data["rooms"].values():
        if str(cid) in (room.get("members") or []):
            room["members"] = [m for m in room["members"] if m != str(cid)]
            room.setdefault("removed", []).append(str(cid))
    return entry


def revoke_contact(cid: str, *, by: str, now: float | None = None) -> dict | None:
    now = time.time() if now is None else now
    with _editing() as data:
        return _revoke_in(data, cid, by=by, now=now)


def revoke_person(person: str, *, by: str, label: str = "all",
                  now: float | None = None) -> list[dict]:
    """Revoke ``person``'s phones: all of them, or the ones whose label or
    contact id is ``label``."""
    now = time.time() if now is None else now
    gone = []
    with _editing() as data:
        for cid, rec in list(data["contacts"].items()):
            if rec.get("person") != person:
                continue
            if label not in {"all", "alle", "*"} and label not in {cid, rec.get("label")}:
                continue
            gone.append(_revoke_in(data, cid, by=by, now=now))
    return [g for g in gone if g]


def expire(now: float | None = None) -> dict:
    """Drop invitations past their time and phones not confirmed in time.
    Returns {"invitations": [(conn, inv)], "contacts": [revoked entries]}."""
    now = time.time() if now is None else now
    out = {"invitations": [], "contacts": []}
    with _editing() as data:
        for conn, inv in list(data["invitations"].items()):
            if now >= float(inv.get("expires") or 0) and inv.get("by") == "gateway":
                out["invitations"].append((conn, dict(inv)))
                del data["invitations"][conn]
                data["retired"][conn] = {"at": now, "by": "gateway", "why": "expired"}
        for cid, rec in list(data["contacts"].items()):
            if rec.get("state") == "unverified" and now - float(rec.get("since") or 0) > pairing.VERIFY_TTL_S:
                out["contacts"].append(_revoke_in(data, cid, by="not verified in time", now=now))
        data["retired"] = {c: u for c, u in data["retired"].items()
                           if now - float(u.get("at", 0) if isinstance(u, dict) else u) < 7 * 86400}
        for person, seen in list(data["handled"].items()):
            data["handled"][person] = {r: t for r, t in seen.items()
                                       if now - float(t) < 2 * pairing.REQUEST_TTL_S}
    return out


def uncleaned() -> list[dict]:
    return [dict(r) for r in load().get("revoked") or [] if not r.get("cleaned")]


def mark_cleaned(cid: str, at: float) -> None:
    with _editing() as data:
        for entry in data.get("revoked", []):
            if entry.get("contact") == str(cid) and entry.get("at") == at:
                entry["cleaned"] = True


def handled(person: str, rid: str) -> bool:
    return rid in (load().get("handled") or {}).get(person, {})


def mark_handled(person: str, rid: str, now: float | None = None) -> bool:
    """True the first time ``(person, rid)`` is marked: a replayed request
    id is never acted on twice."""
    now = time.time() if now is None else now
    with _editing() as data:
        seen = data["handled"].setdefault(person, {})
        if rid in seen:
            return False
        seen[rid] = now
        return True


# -- rooms ---------------------------------------------------------------------------

def rooms() -> dict[str, dict]:
    data = load().get('rooms') or {}
    return data if isinstance(data, dict) else {}


def _person_for(data: dict, owner: str | None, person: str | None) -> str:
    """The person a room write is for: the verified owner contact's, which
    must agree with ``person`` when both are given."""
    if owner is not None:
        rec = data["contacts"].get(str(owner))
        if not rec or rec.get("state") != "verified":
            raise ValueError('pairing changed')
        if person is not None and rec["person"] != person:
            raise ValueError('pairing changed')
        return rec["person"]
    if person is None or not any(r.get("person") == person and r.get("state") == "verified"
                                 for r in data["contacts"].values()):
        raise ValueError('pairing changed')
    return person


def put_room(chat: str, session_id: str | None, title: str, *, owner: str | None = None,
             person: str | None = None, ready: bool = False) -> None:
    """Bind a group permanently to one person's session (None: their forum);
    a rename cannot redirect its conversation, a group never changes person."""
    with _editing() as data:
        who = _person_for(data, owner, person)
        bindings = data['rooms']
        old = bindings.get(chat)
        if old and (old['session_id'] != session_id or old.get('person', OWNER) != who):
            raise ValueError('a chat cannot be rebound to another session')
        if any(r['session_id'] == session_id and r.get('person', OWNER) == who and key != chat
               for key, r in bindings.items()):
            raise ValueError('session already has a chat')
        # A group the gateway makes is bound before anyone is invited
        # (ready False): its member list starts empty and grows with each
        # invitation. A room without one predates member lists.
        fresh = {'members': []} if old is None and not ready else {}
        bindings[chat] = {**fresh, **(old or {}), 'session_id': session_id, 'title': title,
                          'ready': ready, 'person': who}


def add_member(chat: str, cid: str) -> None:
    with _editing() as data:
        room = data['rooms'].get(chat)
        rec = data['contacts'].get(str(cid))
        if not room or not rec or rec.get('person') != room.get('person', OWNER):
            raise ValueError('pairing or room changed')
        members = room.setdefault('members', [])
        if str(cid) not in members:
            members.append(str(cid))


def configure_room(chat: str, *, owner: str | None = None, person: str | None = None,
                   **settings) -> None:
    if set(settings) - {'forum', 'notifications', 'relay'}:
        raise ValueError('unknown room setting')
    with _editing() as data:
        who = _person_for(data, owner, person)
        room = data['rooms'].get(chat)
        if room is None or room.get('person', OWNER) != who:
            raise ValueError('pairing or room changed')
        room.update(settings)


def room_of_session(session_id: str | None, person: str = OWNER) -> str | None:
    """The chat bound to one person's session, or None."""
    return next((chat for chat, room in rooms().items()
                 if room.get('session_id') == session_id
                 and room.get('person', OWNER) == person), None)


def forum_of(person: str) -> str | None:
    return room_of_session(None, person)


def schedules_room() -> str | None:
    """The chat the owner's schedules are managed and reported in; None
    before it is opened."""
    return str(load().get('schedules_room') or '') or None


def set_schedules_room(chat: str, *, owner: str) -> None:
    with _editing() as data:
        if _person_for(data, owner, OWNER) != OWNER:
            raise ValueError('pairing changed')
        data['schedules_room'] = chat


def schedule_seq() -> int | None:
    value = load().get('schedule_seq')
    return value if isinstance(value, int) else None


def set_schedule_seq(seq: int, *, owner: str) -> None:
    with _editing() as data:
        _person_for(data, owner, OWNER)
        data['schedule_seq'] = int(seq)


def forget_room(chat: str) -> None:
    with _editing() as data:
        data['rooms'].pop(chat, None)


def claim_message(chat: str, message_id: str) -> bool:
    """Durable at-most-once dispatch, including across gateway restarts."""
    import sqlite3
    directory = state_dir()
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    with sqlite3.connect(directory / 'received.sqlite3') as db:
        db.execute('CREATE TABLE IF NOT EXISTS received (chat TEXT, id TEXT, PRIMARY KEY(chat,id))')
        return db.execute('INSERT OR IGNORE INTO received VALUES (?,?)', (chat, message_id)).rowcount == 1
