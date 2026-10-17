"""What a turn waits for, and when that wait is satisfied.

Package KG1 phase 2. A turn that WAITS — for a message from a session, a
job id, an approval id, a time — registers that here. Two effects:

1. turn_continuation.should_continue sends NO nudge while a wait is
   registered: "let me check every minute" must not be nudged forever;
   the wait is the mechanism that brings the session back.
2. The wait WAKES the session when its event arrives, once. due(now,
   events) matches incoming events against the stored waits, consumes
   what it fires, and the caller (the engine's idle path) delivers
   wait_note(w) as the one-line note naming the event.

The store is a small JSON file under the redirected user state
(state_paths.USER_STATE_RESOLVERS keeps every user-state module pointed
at a scratch root in tests; this module follows the same pattern with a
per-call resolver, never an import-time path). A time wait holds an epoch
value; the other kinds hold an exact string, matched exactly — a wait
names ONE sender/job/approval, and a differently-cased sender is a
different event (the mailbox folding in session_messages is about
DELIVERY, not about who the event came from).
"""

from __future__ import annotations

import json
import os
import time
import uuid
from pathlib import Path

#: The kinds a turn may wait for. Time waits fire when ``now >= value``;
#: every other kind fires when an event of the same kind and value arrives.
WAIT_KINDS = ("message_from", "job_id", "approval_id", "time")

_KH = ("message_from", "job_id", "approval_id")
_NOTE_BY_KIND = {
    "message_from": "A message from {value} arrived while the session waited.",
    "job_id": "Job {value} reported while the session waited.",
    "approval_id": "Approval {value} was answered while the session waited.",
    "time": "The wait until {value} is over.",
}


def _store_path(key: str = "") -> Path:
    """Resolved per call, so a redirected home moves it (wave-15a rule).

    One file per session key (the inbox pattern in session_messages): a
    wait belongs to the session that registered it. Reviewer nacht-s25,
    finding #4: keyed only by workspace, one session's wait suppressed
    the nudges of another. The key is the session's presence key -- the
    same address pause flags and mailboxes use.
    """
    from . import state_paths
    safe = "".join(c if c.isalnum() or c in "-_" else "_"
                   for c in str(key or ""))[:80] or "_default"
    return state_paths.workspace_state_dir("") / f"waits_{safe}.json"


def _load(key: str = "") -> list[dict]:
    try:
        raw = json.loads(_store_path(key).read_text(encoding="utf-8"))
    except (FileNotFoundError, NotADirectoryError, ValueError, OSError):
        return []
    return [w for w in raw if isinstance(w, dict)] if isinstance(raw, list) else []


def _save(waits: list[dict], key: str = "") -> None:
    from .state_paths import ensure_dir, write_text_atomic
    ensure_dir(_store_path(key).parent)
    write_text_atomic(_store_path(key), json.dumps(waits))


def register(kind: str, value: str, *, key: str = "") -> dict:
    """Register one wait for session ``key``. Returns the wait row.

    Raises ValueError for an unknown kind or an empty value: a wait with
    nothing to wait for can never fire and would suppress nudges forever.
    """
    if kind not in WAIT_KINDS:
        raise ValueError(f"wait_registry: unknown wait kind {kind!r}")
    value = str(value or "").strip()
    if not value:
        raise ValueError(f"wait_registry: empty value for {kind!r}")
    wait = {
        "id": uuid.uuid4().hex[:12],
        "kind": kind,
        "value": value if kind == "time" else value,
        "registered_at": time.time(),
    }
    waits = _load(key)
    waits.append(wait)
    _save(waits, key)
    return wait


def wait_note(wait: dict) -> str:
    """The one-line note a wake names its event with."""
    kind = str(wait.get("kind"))
    return _NOTE_BY_KIND.get(kind, "{kind}: {value}").format(
        kind=kind, value=wait.get("value"))


def due(now: float, events: list, *, key: str = "") -> list[dict]:
    """The waits of session ``key`` that ``now`` and ``events`` satisfy,
    removed from that session's store.

    ``events`` is an iterable of (kind, value) pairs the caller gathered
    this round (messages taken, job records seen, approvals answered).
    A time wait is satisfied when ``now >= value`` and needs no event.
    Time waits fire before event waits, so the idle prompt shows the
    deadline's note ahead of the mail that also arrived. What fires is
    consumed: the wake happens once.
    """
    now = float(now)
    waits = _load(key)
    hits: list[dict] = []
    rest: list[dict] = []
    for w in waits:
        kind, value = str(w.get("kind")), str(w.get("value"))
        if kind == "time":
            try:
                if now >= float(value):
                    hits.append(w)
                    continue
            except ValueError:
                pass
        elif any(k == kind and v == value for k, v in events):
            hits.append(w)
            continue
        rest.append(w)
    hits.sort(key=lambda w: (0 if w["kind"] == "time" else 1,
                             w.get("registered_at") or 0.0))
    if hits:
        _save(rest, key)
    return hits


class SessionWaits:
    """The person's view of the waits: list and clear (the status line and
    /state clear-wait read these). Keyed like the registry itself."""

    def __init__(self, key: str = ""):
        self.key = str(key or "")

    def list(self) -> list[dict]:
        return _load(self.key)

    def clear(self, wait_id: str) -> bool:
        waits = _load(self.key)
        keep = [w for w in waits if w.get("id") != str(wait_id)]
        if len(keep) == len(waits):
            return False
        _save(keep, self.key)
        return True
