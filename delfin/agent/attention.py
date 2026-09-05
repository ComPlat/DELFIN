"""Durable attention/event layer for unattended runs.

When an agent run blocks on the user (confirm dialog, plan approval,
ask_user_question) or finishes/fails while nobody is watching the chat,
the moment must not be lost in scrollback. This module provides:

  - ``emit_attention(kind, ...)``
        Appends a durable event to ``~/.delfin/attention_inbox.jsonl``
        (atomic temp-file + ``os.replace``, the ``memory_store.
        _atomic_write`` pattern) AND best-effort delivers it out-of-band
        through the transports the codebase already has:

          * desktop notification (``notify.send_notification`` — the
            same helper the ``push_notification`` tool and the job
            monitor's ``announce()`` use),
          * the job monitor's webhook (``agent.job_monitor.webhook_url``
            setting via ``notify.send_remote_trigger`` with
            ``override_url`` — exactly like ``job_monitor.announce``),
          * an optional user shell hook ``agent.attention.
            notify_command`` — run detached (daemon thread, 5 s cap,
            argv list, never ``shell=True``) with the event title as
            the final argument, so users can wire desktop/matrix/email
            notifiers of their own.

  - The INBOX: ``list_pending()`` / ``resolve()`` / ``drain_resolved()``.
        Pending events survive process restarts (file-backed).
        ``resolve`` records the user's late answer (e.g. the reply to a
        question that timed out) and marks the event resolved.
        ``drain_resolved`` hands answered events to a session exactly
        once — an ``acknowledged`` flag persisted in the same file, the
        restart-safe pattern of ``bash_jobs.drain_finished_events``.

  - The USER SURFACE: ``render_inbox()`` / ``answer_item()`` /
        ``dismiss_item()`` / ``clear_all()``.
        Pure helpers for a dashboard slash command (``/attention``):
        compact grouped rendering with visible ids, answering a parked
        item (the engine later injects the answer via
        ``drain_resolved``), and dismissing items without replaying
        them to the model. All four never raise — they return an error
        string / ``{"ok": False, "error": ...}`` instead, because a
        broken inbox must never break the chat surface.

Delivery is best-effort and never raises; only the durable append (the
core contract) propagates errors to the caller.

TWO CLASSES OF EVENT, and nearly every rule below follows from the
difference:

  * BLOCKING (confirm/question/plan) — a turn is waiting for an answer
    only the user can give. These never expire, a blanket clear refuses
    to touch them without an explicit opt-in, and the answer is kept
    until a session has actually consumed it.
  * NOTICES (finished/failed/budget/security) — they report something
    that already happened. Nobody ever "answers" them, so nothing ever
    resolved them and they accumulated forever: a real inbox reached
    893 pending records, 15 days deep, where the one parked question
    was invisible among them. Notices therefore expire (per-kind TTL)
    and are additionally capped, which is what keeps the file small
    enough that the read-modify-write above stays cheap.
"""

from __future__ import annotations

import contextlib
import json
import os
import secrets
import shlex
import shutil
import subprocess
import tempfile
import threading
import time
from pathlib import Path
from typing import Any, Optional


ATTENTION_KINDS = (
    "confirm_pending",
    "question_pending",
    "plan_pending",
    "security_alert",
    "run_finished",
    "run_failed",
    "budget_warning",
)

# Kinds that mean "the agent is BLOCKED on the user right now" — a turn
# is waiting for an answer only the user can give.
_BLOCKING_KINDS = frozenset({
    "confirm_pending", "question_pending", "plan_pending",
})

# Kinds delivered with critical urgency. A failed run and a containment
# alert are loud, but they are not BLOCKING: nobody is waiting for an
# answer, so they follow the notice rules (expiry, blanket clear).
_URGENT_KINDS = _BLOCKING_KINDS | {"security_alert", "run_failed"}

# Notices and how long an unread one stays interesting. "Your run
# finished" is worthless after a couple of days; a failure or a
# containment alert deserves a week. Only pending notices expire —
# nothing the user still has to answer is ever dropped by age.
_NOTICE_TTL_S: dict[str, float] = {
    "run_finished": 2 * 86400,
    "run_failed": 7 * 86400,
    "budget_warning": 7 * 86400,
    "security_alert": 7 * 86400,
}
# Ceiling on pending notices regardless of age: a scheduler firing every
# few minutes can produce hundreds inside one TTL window, and every read
# of this file parses all of them.
_MAX_PENDING_NOTICES = 200

_INBOX_FILENAME = "attention_inbox.jsonl"
_FILE_LOCK = threading.Lock()
_NOTIFY_CMD_TIMEOUT_S = 5
# Retention: delivered events are dropped after this age; their number is
# additionally capped so the file cannot grow without bound. Blocking
# events and answers that no session has consumed yet are never pruned.
_ACKED_MAX_AGE_S = 14 * 86400
_MAX_NON_PENDING = 500

# Transport names recorded on an event / reported by transport_status().
TRANSPORTS = ("desktop", "webhook", "hook")

# Handle of the most recent notify_command worker thread (tests join it;
# production ignores it — the thread is a daemon and never blocks emit).
_LAST_NOTIFY_THREAD: Optional[threading.Thread] = None


# ---------------------------------------------------------------------------
# Storage (JSONL, atomic replace)
# ---------------------------------------------------------------------------

def _inbox_path() -> Path:
    # Resolved per call (not at import) so tests can monkeypatch Path.home.
    return Path.home() / ".delfin" / _INBOX_FILENAME


@contextlib.contextmanager
def _inbox_locked():
    """Serialize a whole load → mutate → write cycle on the inbox.

    ``_FILE_LOCK`` is a thread lock and this file is per USER: the
    dashboard, a headless run and a second CLI session all append to it
    from different PROCESSES. The atomic write kept them from reading a
    torn file and did nothing about the lost update — measured, four
    processes emitting forty events each left 48 of 160, and the events
    that vanished are exactly the ones telling the user the agent is
    blocked on them.
    """
    from .bash_jobs import cross_process_lock

    with _FILE_LOCK, cross_process_lock(_inbox_path()):
        yield


def _load_events() -> list[dict]:
    """Tolerant read: skip torn/foreign lines instead of failing the inbox."""
    try:
        text = _inbox_path().read_text(encoding="utf-8")
    except OSError:
        return []
    events: list[dict] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(rec, dict) and rec.get("id"):
            events.append(rec)
    return events


def _write_events(events: list[dict]) -> None:
    """Atomic rewrite (temp file + ``os.replace``, memory_store pattern):
    a dashboard and a headless daemon can share this file — a reader racing
    a plain truncate-write would observe an empty or torn inbox."""
    path = _inbox_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            for rec in events:
                fh.write(json.dumps(rec, ensure_ascii=False, default=str))
                fh.write("\n")
        os.replace(tmp, path)
        try:
            os.chmod(path, 0o600)
        except OSError:
            pass
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def _created_at(ev: dict) -> float:
    try:
        return float(ev.get("created_at") or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _is_notice(ev: dict) -> bool:
    """A record nobody is expected to answer (finished/failed/budget/
    security) — as opposed to a blocking confirm/question/plan."""
    return str(ev.get("kind", "")) in _NOTICE_TTL_S


def _notice_expired(ev: dict, now: float) -> bool:
    """A pending notice past its kind's TTL.

    Only pending ones: a notice the user answered is an answer in
    flight, and the age of the thing it reported on says nothing about
    whether the agent has received it.
    """
    if ev.get("status") != "pending":
        return False
    ttl = _NOTICE_TTL_S.get(str(ev.get("kind", "")))
    return ttl is not None and _created_at(ev) < now - float(ttl)


def _visible(events: list[dict], now: float) -> list[dict]:
    """Read-path view: expired notices are gone the moment they expire.

    Applied on read as well as on write so an inbox that nobody has
    emitted into lately still reports the truth — the file catches up at
    the next write (or via :func:`prune_inbox`).
    """
    return [ev for ev in events if not _notice_expired(ev, now)]


def _prune(events: list[dict], now: float) -> list[dict]:
    """Retention. Four rules, in order:

    1. delivered events (acknowledged) older than ``_ACKED_MAX_AGE_S`` go;
    2. pending NOTICES past their kind's TTL go;
    3. pending notices are capped at ``_MAX_PENDING_NOTICES``, oldest first;
    4. delivered events are capped at ``_MAX_NON_PENDING``, oldest first.

    What is never dropped: a pending BLOCKING event (the user's open
    to-do), and a resolved event no session has consumed yet. Rule 4
    used to count the latter, and since it drops the oldest first it
    deleted precisely the long-parked question the user had finally
    answered — on its way to the agent, one write before delivery.
    """
    cutoff = now - _ACKED_MAX_AGE_S
    kept = [
        ev for ev in events
        if not (ev.get("acknowledged")
                and float(ev.get("resolved_at") or ev.get("created_at") or 0.0)
                < cutoff)
    ]
    kept = _visible(kept, now)

    notices = [ev for ev in kept
               if ev.get("status") == "pending" and _is_notice(ev)]
    overflow = len(notices) - _MAX_PENDING_NOTICES
    if overflow > 0:
        notices.sort(key=_created_at)
        drop_ids = {ev["id"] for ev in notices[:overflow]}
        kept = [ev for ev in kept if ev["id"] not in drop_ids]

    delivered = [ev for ev in kept
                 if ev.get("status") != "pending" and ev.get("acknowledged")]
    overflow = len(delivered) - _MAX_NON_PENDING
    if overflow > 0:
        delivered.sort(key=_created_at)
        drop_ids = {ev["id"] for ev in delivered[:overflow]}
        kept = [ev for ev in kept if ev["id"] not in drop_ids]
    return kept


def prune_inbox() -> dict:
    """Apply retention now and rewrite the file. Never raises.

    Emitting applies the same rules, so this is only needed to compact an
    inbox that predates them (or one nothing emits into any more).
    Returns ``{"ok": True, "removed": n, "kept": m}``.
    """
    try:
        now = time.time()
        with _FILE_LOCK:
            events = _load_events()
            kept = _prune(events, now)
            if len(kept) != len(events):
                _write_events(kept)
        return {"ok": True, "removed": len(events) - len(kept),
                "kept": len(kept)}
    except Exception as exc:
        return {"ok": False, "error": str(exc), "removed": 0, "kept": 0}


# ---------------------------------------------------------------------------
# Settings + out-of-band delivery
# ---------------------------------------------------------------------------

def attention_settings(settings: dict | None = None) -> dict:
    """``agent.attention`` section of the user settings (job_monitor
    pattern). ``notify_command`` is a user-authored argv prefix; the event
    title is appended as the final argument."""
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    cfg = ((settings or {}).get("agent") or {}).get("attention") or {}
    return {
        "notify_command": str(cfg.get("notify_command", "") or ""),
        "desktop": bool(cfg.get("desktop", True)),
    }


def _run_notify_command(command: str, title: str) -> Optional[threading.Thread]:
    """Run the user's notify hook detached: argv list (no shell), output
    discarded, 5 s cap enforced in a daemon thread so emit never blocks
    and a hung hook never leaks. Never raises."""
    global _LAST_NOTIFY_THREAD
    try:
        argv = shlex.split(command)
    except ValueError:
        return None
    if not argv:
        return None
    argv.append(str(title))

    def _worker() -> None:
        try:
            subprocess.run(
                argv,
                timeout=_NOTIFY_CMD_TIMEOUT_S,
                check=False,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
        except Exception:
            pass

    t = threading.Thread(
        target=_worker, daemon=True, name="attention-notify-command")
    t.start()
    _LAST_NOTIFY_THREAD = t
    return t


def _webhook_url() -> str:
    """The job monitor's webhook URL — one setting for both surfaces."""
    try:
        from delfin.agent.job_monitor import monitor_settings
        return str(monitor_settings().get("webhook_url", "") or "")
    except Exception:
        return ""


def transport_status(settings: dict | None = None) -> dict:
    """Can this machine tell the user anything at all? Probe only — sends
    nothing.

    Returns ``{"usable": [names], "detail": {name: str}}``. The point of
    the report is the empty list: on a headless login node with no
    desktop notifier, no webhook and no hook command the whole fan-out
    below is a no-op, and every surface still reported healthy while the
    agent sat blocked on a question nobody could see.
    """
    cfg = {}
    try:
        cfg = attention_settings(settings)
    except Exception:
        cfg = {}
    detail: dict[str, str] = {}
    usable: list[str] = []

    if not cfg.get("desktop", True):
        detail["desktop"] = "disabled (agent.attention.desktop = false)"
    else:
        try:
            from delfin.agent import notify as _notify
            probe = _notify.notification_backend()
        except Exception as exc:
            probe = {"available": False, "detail": str(exc)}
        detail["desktop"] = str(probe.get("detail", ""))
        if probe.get("available"):
            usable.append("desktop")

    url = _webhook_url()
    if not url:
        detail["webhook"] = "agent.job_monitor.webhook_url not configured"
    elif not url.lower().startswith("https://"):
        detail["webhook"] = "configured but not https — send_remote_trigger " \
                            "refuses it"
    else:
        detail["webhook"] = "configured"
        usable.append("webhook")

    command = str(cfg.get("notify_command", "") or "")
    if not command:
        detail["hook"] = "agent.attention.notify_command not configured"
    else:
        try:
            argv = shlex.split(command)
        except ValueError:
            argv = []
        if not argv:
            detail["hook"] = f"unparseable command: {command!r}"
        elif not (shutil.which(argv[0]) or Path(argv[0]).exists()):
            detail["hook"] = f"{argv[0]} not found"
        else:
            detail["hook"] = f"{argv[0]}"
            usable.append("hook")
    return {"usable": usable, "detail": detail}


def _deliver(event: dict, *, skip: frozenset[str] = frozenset()) -> list[str]:
    """Best-effort out-of-band fan-out. Never raises, never blocks beyond
    the notify/webhook helpers' own short timeouts.

    Returns the transports that actually accepted the event. Each helper
    already reported success and every caller threw the answer away, so
    "the user was notified" was an assumption the code never once
    checked; the returned names are recorded on the event instead.

    ``skip`` names transports the caller has already driven itself (see
    ``job_monitor.announce``), so the user is not notified twice.
    """
    title = str(event.get("title", "")) or "delfin agent"
    detail = str(event.get("detail", ""))
    kind = str(event.get("kind", ""))
    urgency = "critical" if kind in _URGENT_KINDS else "normal"
    delivered: list[str] = []

    cfg = {}
    try:
        cfg = attention_settings()
    except Exception:
        cfg = {}

    # 1. Desktop notification — same helper as the push_notification tool
    #    and job_monitor.announce().
    if cfg.get("desktop", True) and "desktop" not in skip:
        try:
            from delfin.agent import notify as _notify
            if _notify.send_notification(
                    f"DELFIN: {title}", detail[:400], urgency=urgency):
                delivered.append("desktop")
        except Exception:
            pass

    # 2. Webhook — reuse the job monitor's URL setting + helper verbatim
    #    (agent.job_monitor.webhook_url via send_remote_trigger).
    webhook_url = "" if "webhook" in skip else _webhook_url()
    if webhook_url:
        try:
            from delfin.agent import notify as _notify
            result = _notify.send_remote_trigger(
                {
                    "event": f"attention.{kind}",
                    "id": event.get("id", ""),
                    "title": title,
                    "detail": detail[:800],
                    "session": event.get("session_id", ""),
                    "workspace": event.get("workspace", ""),
                    "options": event.get("options"),
                },
                override_url=webhook_url,
            )
            if getattr(result, "sent", False):
                delivered.append("webhook")
        except Exception:
            pass

    # 3. User shell hook (agent.attention.notify_command).
    notify_command = str(cfg.get("notify_command", "") or "")
    if notify_command and "hook" not in skip:
        try:
            # Detached by design: "handed to the hook" is the strongest
            # claim available without waiting for it.
            if _run_notify_command(notify_command, title) is not None:
                delivered.append("hook")
        except Exception:
            pass
    return delivered


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def emit_attention(
    kind: str,
    *,
    session_id: str = "",
    title: str,
    detail: str = "",
    options: Optional[list] = None,
    workspace: Optional[str | Path] = None,
    delivered_by_caller: Optional[dict[str, bool]] = None,
) -> str:
    """Record a durable attention event and fan it out out-of-band.

    Returns the event id. The append is the contract (raises on disk
    failure); every delivery transport is best-effort and never raises.

    ``delivered_by_caller`` maps transport name → whether the CALLER
    already sent this event that way (``job_monitor.announce`` keeps its
    own job-shaped webhook payload). Those transports are skipped here,
    and the ones that succeeded are recorded as delivered.

    The transports that accepted the event are written back onto the
    record (``delivered``), so "was the user actually told?" is a
    question the inbox and the doctor can answer instead of assume.
    """
    if kind not in ATTENTION_KINDS:
        raise ValueError(
            f"unknown attention kind {kind!r}; expected one of "
            f"{ATTENTION_KINDS}")
    now = time.time()
    by_caller = dict(delivered_by_caller or {})
    event: dict[str, Any] = {
        "id": f"att-{int(now * 1000):x}-{secrets.token_hex(3)}",
        "kind": kind,
        "session_id": str(session_id or ""),
        "title": str(title)[:200],
        "detail": str(detail)[:2000],
        "options": ([str(o) for o in options] if options else None),
        "workspace": str(workspace) if workspace else "",
        "created_at": now,
        "status": "pending",
        "answer": None,
        "resolved_at": None,
        "acknowledged": False,
        "delivered": [name for name, ok in by_caller.items() if ok],
    }
    with _inbox_locked():
        events = _load_events()
        events.append(event)
        _write_events(_prune(events, now))
    delivered = _deliver(event, skip=frozenset(by_caller))
    if delivered:
        _record_delivery(event["id"], delivered)
    return event["id"]


def _record_delivery(event_id: str, transports: list[str]) -> None:
    """Note on the stored event which transports accepted it. Best-effort:
    a delivery that could not be recorded must not undo the delivery."""
    try:
        with _FILE_LOCK:
            events = _load_events()
            for ev in events:
                if ev.get("id") == event_id:
                    names = list(ev.get("delivered") or [])
                    names += [t for t in transports if t not in names]
                    ev["delivered"] = names
                    _write_events(events)
                    return
    except Exception:
        pass


def list_pending(kind: Optional[str] = None) -> list[dict]:
    """Unresolved events (oldest first) — what the user still has to act
    on, across restarts. Optionally filtered by kind.

    Expired notices are filtered out here as well as pruned on write, so
    an inbox nobody has emitted into for a week does not report a week of
    stale "run finished" lines as open to-dos.
    """
    now = time.time()
    with _inbox_locked():
        events = _load_events()
    pending = [
        dict(ev) for ev in _visible(events, now)
        if ev.get("status") == "pending"
        and (kind is None or ev.get("kind") == kind)
    ]
    pending.sort(key=_created_at)
    return pending


def list_undelivered() -> list[dict]:
    """Answers the user has given that NO session has consumed yet.

    Neither pending nor delivered, so every surface that filtered on
    ``status == "pending"`` showed nothing: the user believed they had
    answered, the agent never saw it, and there was no place left where
    either of them could find it.
    """
    with _FILE_LOCK:
        events = _load_events()
    out = [dict(ev) for ev in events
           if ev.get("status") == "resolved" and not ev.get("acknowledged")]
    out.sort(key=_created_at)
    return out


def resolve(
    event_id: str,
    answer: Any = None,
    *,
    acknowledged: bool = False,
) -> bool:
    """Mark a pending event resolved, recording the user's ``answer``.

    ``acknowledged=True`` marks it consumed immediately — used when the
    answer was already delivered in-band (e.g. a live confirm click), so
    ``drain_resolved`` does not replay it to the model. Returns True when
    a pending event with this id was resolved."""
    with _inbox_locked():
        events = _load_events()
        for ev in events:
            if ev.get("id") == event_id and ev.get("status") == "pending":
                ev["status"] = "resolved"
                ev["answer"] = answer
                ev["resolved_at"] = time.time()
                ev["acknowledged"] = bool(acknowledged)
                _write_events(events)
                return True
    return False


def drain_resolved(
    session_id: str = "",
    *,
    workspace: Optional[str | Path] = None,
    kinds: Optional[frozenset[str] | set[str] | tuple[str, ...]] = None,
    include_unrouted: bool = True,   # accepted for call-site compatibility
) -> list[dict]:
    """Answered events for this agent — exactly once.

    Matching is by WORKSPACE and kind, never by session id. The session
    id an event was parked with is not stable: it is re-minted on a new
    cycle, cleared when the CLI backend restarts its process, and
    overwritten mid-turn by the id the backend reports. Draining by exact
    match therefore stranded the answer forever — the user had answered,
    the agent never saw it, and because the listing filtered on
    ``pending`` the record was invisible to all three of the inbox view,
    the user and the doctor. The id is kept on the record as context;
    ``session_id`` here is accepted for call-site compatibility and no
    longer filters anything.

    An event carrying no workspace matches any (the confirm broker has no
    workspace context); an explicit ``workspace`` argument limits the
    drain to events parked in that workspace.

    The ``acknowledged`` flag is persisted in the inbox file, so the
    exactly-once contract survives restarts (``bash_jobs.
    drain_finished_events`` pattern): a new process only sees what no
    previous process drained.
    """
    want_ws = str(workspace) if workspace else ""
    want_kinds = frozenset(kinds) if kinds else None
    out: list[dict] = []
    with _inbox_locked():
        events = _load_events()
        changed = False
        for ev in events:
            if ev.get("status") != "resolved" or ev.get("acknowledged"):
                continue
            if want_kinds is not None and str(ev.get("kind", "")) not in want_kinds:
                continue
            ev_ws = str(ev.get("workspace") or "")
            if want_ws and ev_ws and ev_ws != want_ws:
                continue
            ev["acknowledged"] = True
            changed = True
            out.append(dict(ev))
        if changed:
            _write_events(events)
    out.sort(key=_created_at)
    return out


# ---------------------------------------------------------------------------
# User surface (dashboard /attention command) — pure helpers, never raise
# ---------------------------------------------------------------------------

# Display order: kinds that block the agent first, informational last.
_KIND_ORDER = (
    "confirm_pending", "question_pending", "plan_pending",
    "security_alert", "run_failed", "run_finished", "budget_warning",
)
_KIND_LABELS = {
    "confirm_pending": "Confirmations waiting",
    "question_pending": "Questions waiting",
    "plan_pending": "Plans awaiting approval",
    "security_alert": "Security alerts",
    "run_failed": "Failed runs",
    "run_finished": "Finished runs",
    "budget_warning": "Budget warnings",
}


def _fmt_age(seconds: float) -> str:
    """Compact age: 45s / 12m / 3h / 5d. Negative/garbage clamps to 0s."""
    try:
        s = max(0.0, float(seconds))
    except (TypeError, ValueError):
        s = 0.0
    if s < 60:
        return f"{int(s)}s"
    if s < 3600:
        return f"{int(s // 60)}m"
    if s < 86400:
        return f"{int(s // 3600)}h"
    return f"{int(s // 86400)}d"


def _match_pending(item_id: str) -> tuple[Optional[dict], str]:
    """Find one pending event by full id or unique id prefix.

    Returns ``(event, "")`` on success, ``(None, error)`` otherwise —
    ids like ``att-19c8...-a1b2c3`` are painful to type in full, so a
    surface accepts any unambiguous prefix."""
    token = str(item_id or "").strip()
    if not token:
        return None, "missing item id"
    pending = list_pending()
    for ev in pending:
        if ev.get("id") == token:
            return ev, ""
    matches = [
        ev for ev in pending if str(ev.get("id", "")).startswith(token)]
    if len(matches) == 1:
        return matches[0], ""
    if not matches:
        return None, f"no pending item with id {token!r}"
    return None, (
        f"id prefix {token!r} is ambiguous "
        f"({len(matches)} pending items match)")


def _render_undelivered(kind: Optional[str], now: float) -> list[str]:
    """The "answered, not yet delivered" section.

    Answers waiting for a session to pick them up are neither pending nor
    done. Without a line of their own the user sees an inbox that looks
    finished while their answer has not reached the agent.
    """
    try:
        waiting = [ev for ev in list_undelivered()
                   if kind is None or ev.get("kind") == kind]
    except Exception:
        return []
    if not waiting:
        return []
    lines = [f"Answered, not yet delivered ({len(waiting)}) — the agent "
             "picks these up on its next turn:"]
    for ev in waiting[:10]:
        age = _fmt_age(now - float(ev.get("resolved_at") or now))
        title = str(ev.get("title", "")).strip() or "(no title)"
        lines.append(f"  {ev.get('id', '?')}  (answered {age} ago)  {title}")
    if len(waiting) > 10:
        lines.append(f"  ... and {len(waiting) - 10} more")
    return lines


def render_inbox(kind: Optional[str] = None, *, limit: int = 50) -> str:
    """Compact, grouped, id-visible view of the inbox.

    Groups by kind (blocking kinds first), oldest first inside each
    group, at most ``limit`` items total (a footer notes the overflow),
    followed by the answers that have not reached a session yet.
    Never raises — a broken inbox renders as an error line."""
    try:
        pending = list_pending(kind)
    except Exception as exc:
        return f"Attention inbox unavailable: {exc}"
    try:
        scope = f" of kind {kind!r}" if kind else ""
        now = time.time()
        if not pending:
            waiting = _render_undelivered(kind, now)
            head = f"Attention inbox: no pending items{scope}."
            return "\n".join([head] + waiting) if waiting else head
        try:
            cap = max(1, int(limit))
        except (TypeError, ValueError):
            cap = 50
        blocking = sum(
            1 for ev in pending if ev.get("kind") in _BLOCKING_KINDS)
        lines = [
            f"Attention inbox — {len(pending)} pending{scope}"
            + (f" ({blocking} blocking the agent)" if blocking else "")]
        shown = 0
        kinds_seen = [k for k in _KIND_ORDER
                      if any(ev.get("kind") == k for ev in pending)]
        kinds_seen += sorted(
            {str(ev.get("kind", "")) for ev in pending}
            - set(_KIND_ORDER))
        for k in kinds_seen:
            group = [ev for ev in pending if ev.get("kind") == k]
            if not group or shown >= cap:
                continue
            lines.append(f"{_KIND_LABELS.get(k, k)} ({len(group)}):")
            for ev in group:
                if shown >= cap:
                    break
                age = _fmt_age(now - float(ev.get("created_at") or now))
                title = str(ev.get("title", "")).strip() or "(no title)"
                lines.append(f"  {ev.get('id', '?')}  ({age})  {title}")
                opts = ev.get("options")
                if opts:
                    lines.append(
                        "      options: " + " | ".join(map(str, opts)))
                shown += 1
        if shown < len(pending):
            lines.append(f"  ... and {len(pending) - shown} more")
        lines.extend(_render_undelivered(kind, now))
        lines.append(
            "Answer: /attention answer <id> <text>   "
            "Dismiss: /attention dismiss <id|all>")
        if blocking:
            lines.append(
                "'dismiss all' clears notices only — dismiss a blocking "
                "item by its own id.")
        return "\n".join(lines)
    except Exception as exc:
        return f"Attention inbox unavailable: {exc}"


def answer_item(item_id: str, text: str) -> dict:
    """Record the user's late answer to a pending item (id or unique
    prefix). The engine injects it into the next turn via
    ``drain_resolved`` — so ``acknowledged`` stays False here. Never
    raises; returns ``{"ok": bool, ...}``."""
    try:
        answer = str(text or "").strip()
        if not answer:
            return {"ok": False, "error": "empty answer text"}
        target, err = _match_pending(item_id)
        if target is None:
            return {"ok": False, "error": err}
        if not resolve(target["id"], answer=answer):
            return {"ok": False,
                    "error": f"item {target['id']} is no longer pending"}
        return {"ok": True, "id": target["id"],
                "kind": str(target.get("kind", "")),
                "title": str(target.get("title", ""))}
    except Exception as exc:
        return {"ok": False, "error": str(exc)}


def dismiss_item(item_id: str) -> dict:
    """Clear a pending item without answering it (id or unique prefix).

    Resolved with ``acknowledged=True`` so ``drain_resolved`` never
    replays it to the model. Never raises."""
    try:
        target, err = _match_pending(item_id)
        if target is None:
            return {"ok": False, "error": err}
        if not resolve(target["id"], answer=None, acknowledged=True):
            return {"ok": False,
                    "error": f"item {target['id']} is no longer pending"}
        return {"ok": True, "id": target["id"],
                "kind": str(target.get("kind", "")),
                "title": str(target.get("title", ""))}
    except Exception as exc:
        return {"ok": False, "error": str(exc)}


def clear_all(kind: Optional[str] = None, *,
              include_blocking: bool = False) -> dict:
    """Dismiss pending items. Same ``acknowledged=True`` semantics as
    ``dismiss_item`` (nothing cleared here is ever replayed to the model).

    A BLANKET clear (no ``kind``) touches notices only. It used to walk
    every kind, so the one gesture the surface offered for an inbox full
    of run notices also destroyed the question or confirm the agent was
    blocked on — right after the timeout message had told the user they
    could still answer it there. Clearing one of those is an explicit
    act: name the kind, pass ``include_blocking=True``, or dismiss the
    item by its id.

    Never raises. Returns ``{"ok": True, "cleared": n, "kept": m}``,
    where ``kept`` counts the blocking items deliberately left alone.
    """
    try:
        pending = list_pending(kind)
        spare_blocking = kind is None and not include_blocking
        targets = [ev for ev in pending
                   if not (spare_blocking
                           and str(ev.get("kind", "")) in _BLOCKING_KINDS)]
        kept = len(pending) - len(targets)
        # One read-modify-write for the whole batch: resolving items one
        # by one re-read and rewrote the entire file per item, which on a
        # real 893-item inbox is 893 full passes over 352 kB.
        target_ids = {str(ev.get("id", "")) for ev in targets}
        cleared = 0
        now = time.time()
        with _FILE_LOCK:
            events = _load_events()
            for ev in events:
                if (ev.get("id") in target_ids
                        and ev.get("status") == "pending"):
                    ev["status"] = "resolved"
                    ev["answer"] = None
                    ev["resolved_at"] = now
                    ev["acknowledged"] = True
                    cleared += 1
            if cleared:
                _write_events(events)
        out = {"ok": True, "cleared": cleared, "kept": kept}
        if kept:
            out["hint"] = (f"{kept} blocking item(s) kept — dismiss one by "
                           f"its id, or clear_all(include_blocking=True)")
        return out
    except Exception as exc:
        return {"ok": False, "error": str(exc), "cleared": 0, "kept": 0}


__all__ = [
    "ATTENTION_KINDS",
    "TRANSPORTS",
    "answer_item",
    "attention_settings",
    "clear_all",
    "dismiss_item",
    "drain_resolved",
    "emit_attention",
    "list_pending",
    "list_undelivered",
    "prune_inbox",
    "render_inbox",
    "resolve",
    "transport_status",
]
