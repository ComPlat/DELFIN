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
import time
import uuid
from pathlib import Path

#: The kinds a turn may wait for. Time waits fire when ``now >= value``;
#: every other kind fires when an event of the same kind and value arrives.
#: ``file`` fires when the file's CONTENT changes (a stat snapshot taken at
#: register time; a touch without a content change does not fire).
WAIT_KINDS = ("message_from", "job_id", "approval_id", "time", "file")

#: Who a registered wait may name. A free-form on_whom is an injection
#: surface (reviewer nacht-s19, wave 18): only these four exist.
ON_WHOM_VALUES = ("session", "operator", "person", "ci")

#: for_what is one line; a second line or markup has no business here.
FOR_WHAT_MAX = 200

_KH = ("message_from", "job_id", "approval_id")
_NOTE_BY_KIND = {
    "message_from": "A message from {value} arrived while the session waited.",
    "job_id": "Job {value} reported while the session waited.",
    "approval_id": "Approval {value} was answered while the session waited.",
    "time": "The wait until {value} is over.",
    "file": "The file {value} changed while the session waited.",
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


def _strip_markup(text: str) -> str:
    """Drop angle-bracket markup so a for_what like
    ``line one <script>alert(1)</script>`` cannot carry HTML into a
    renderer that does not escape. Content in the tags stays as plain
    text; the tags themselves go."""
    import re

    return re.sub(r"<[^>]*>", "", str(text or ""))


def _clean_oneline(value, limit: int = FOR_WHAT_MAX) -> str:
    """*value* as ONE bounded line: newlines and carriage returns become
    spaces, the rest is stripped, capped at *limit* chars. Markup is left
    as data — the renderer escapes (session_health._one_line and the
    cli printer both plain-text it); what must not happen is a second
    line or unbounded length reaching the display."""
    text = " ".join(str(value or "").replace("\r", " ").split("\n"))
    text = " ".join(text.split())
    return _strip_markup(text)[:limit]


def _clean_on_whom(value) -> str:
    """*value* if it is one of ON_WHOM_VALUES (case-folded), else ''.
    An unknown on_whom is rejected, not mapped: a wait that claims to wait
    on a name that is not one of the four is a lie about who can satisfy
    it."""
    text = str(value or "").strip().lower()
    return text if text in ON_WHOM_VALUES else ""


def _file_signature(path: Path) -> str:
    """A content-bound signature of *path*: size + sha256 of up to
    256 KiB. mtime is deliberately NOT in the signature — a touch without
    a content change must not fire the wait (target 4). Missing file ->
    'missing' (a later appearance IS a change)."""
    import hashlib

    try:
        if not path.is_file():
            return "missing"
        digest = hashlib.sha256()
        with path.open("rb") as fh:
            digest.update(fh.read(256 * 1024))
        return f"{path.stat().st_size}:{digest.hexdigest()}"
    except OSError:
        return "missing"


def _resolve_workspace_path(raw: str) -> Path:
    """A file-wait value as a workspace-relative path.

    The anchor is the SESSION'S OWN WORKSPACE — the process cwd. Each
    DELFIN session runs in its own workspace (worktree); the state dir
    (workspace_state_dir) is NOT a workspace, and anchoring there would
    silently watch a state file (reviewer nacht-s19, P2 finding).

    Reviewer constraint (wave 18): the value must stay inside that
    workspace — no absolute paths, no ``..`` escapes, no symlink escape.
    Absolute paths and escapes raise ValueError at register time; a
    symlink escape is caught by resolve() against the workspace root.
    """
    text = str(raw or "").strip()
    if not text:
        raise ValueError("wait_registry: empty file path")
    if text.startswith("/") or (len(text) > 1 and text[1] == ":"):
        raise ValueError(
            f"wait_registry: file wait must be workspace-relative, got {text!r}")
    candidate = Path(text)
    if any(part == ".." for part in candidate.parts):
        raise ValueError(
            f"wait_registry: file wait must not contain '..', got {text!r}")
    root = Path.cwd().resolve()
    resolved = (root / candidate).resolve()
    if resolved != root and root not in resolved.parents:
        raise ValueError(
            f"wait_registry: file wait escapes the workspace: {text!r}")
    return resolved


def register(kind: str, value: str, *, key: str = "",
             on_whom: str = "", for_what: str = "") -> dict:
    """Register one wait for session ``key``. Returns the wait row.

    Raises ValueError for an unknown kind or an empty value: a wait with
    nothing to wait for can never fire and would suppress nudges forever.

    ``on_whom`` names who can satisfy the wait — one of ON_WHOM_VALUES
    (a session name is carried in ``value``; the enum only says the KIND
    of satisfier). ``for_what`` is the one line the display shows. Both
    are optional so every existing caller keeps its signature.

    A ``file`` wait takes its content signature NOW (register time, not
    lazily — reviewer nacht-s19: a lazily-read snapshot is a race and a
    write amplification)."""
    if kind not in WAIT_KINDS:
        raise ValueError(f"wait_registry: unknown wait kind {kind!r}")
    value = str(value or "").strip()
    if not value:
        raise ValueError(f"wait_registry: empty value for {kind!r}")
    whom = _clean_on_whom(on_whom)
    if on_whom and not whom:
        raise ValueError(
            f"wait_registry: on_whom must be one of {ON_WHOM_VALUES}, "
            f"got {on_whom!r}")
    wait = {
        "id": uuid.uuid4().hex[:12],
        "kind": kind,
        "value": value if kind == "time" else value,
        "registered_at": time.time(),
        "on_whom": whom,
        "for_what": _clean_oneline(for_what),
    }
    if kind == "file":
        wait["file_sig"] = _file_signature(_resolve_workspace_path(value))
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
        elif kind == "file":
            # Evaluated HERE, inside due(), against the register-time
            # snapshot: the repl idle path passes events but a file
            # change has no event — the signature comparison IS the
            # wake route. No polling is added: due() runs at the idle
            # path's existing cadence (repl.py:4224).
            try:
                path = _resolve_workspace_path(value)
                current = _file_signature(path)
            except ValueError:
                current = None  # an escaping path can never fire
            if current is not None and current != str(w.get("file_sig")):
                hits.append(w)
                continue
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
