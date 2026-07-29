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
"""

from __future__ import annotations

import json
import os
import secrets
import shlex
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
    "run_finished",
    "run_failed",
    "budget_warning",
)

# Kinds that mean "the agent is blocked on the user right now".
_URGENT_KINDS = frozenset({
    "confirm_pending", "question_pending", "plan_pending", "run_failed",
})

_INBOX_FILENAME = "attention_inbox.jsonl"
_FILE_LOCK = threading.Lock()
_NOTIFY_CMD_TIMEOUT_S = 5
# Retention: acknowledged events are dropped after this age; the file is
# additionally capped so it cannot grow without bound. Pending events are
# never pruned — they are exactly the ones the user still has to act on.
_ACKED_MAX_AGE_S = 14 * 86400
_MAX_NON_PENDING = 500

# Handle of the most recent notify_command worker thread (tests join it;
# production ignores it — the thread is a daemon and never blocks emit).
_LAST_NOTIFY_THREAD: Optional[threading.Thread] = None


# ---------------------------------------------------------------------------
# Storage (JSONL, atomic replace)
# ---------------------------------------------------------------------------

def _inbox_path() -> Path:
    # Resolved per call (not at import) so tests can monkeypatch Path.home.
    return Path.home() / ".delfin" / _INBOX_FILENAME


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


def _prune(events: list[dict], now: float) -> list[dict]:
    """Drop old acknowledged events; cap non-pending records. Pending events
    are kept unconditionally — they are the user's open to-do items."""
    cutoff = now - _ACKED_MAX_AGE_S
    kept = [
        ev for ev in events
        if not (ev.get("acknowledged")
                and float(ev.get("resolved_at") or ev.get("created_at") or 0.0)
                < cutoff)
    ]
    non_pending = [ev for ev in kept if ev.get("status") != "pending"]
    overflow = len(non_pending) - _MAX_NON_PENDING
    if overflow > 0:
        non_pending.sort(
            key=lambda ev: float(ev.get("created_at") or 0.0))
        drop_ids = {ev["id"] for ev in non_pending[:overflow]}
        kept = [ev for ev in kept if ev["id"] not in drop_ids]
    return kept


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


def _deliver(event: dict) -> None:
    """Best-effort out-of-band fan-out. Never raises, never blocks beyond
    the notify/webhook helpers' own short timeouts."""
    title = str(event.get("title", "")) or "delfin agent"
    detail = str(event.get("detail", ""))
    kind = str(event.get("kind", ""))
    urgency = "critical" if kind in _URGENT_KINDS else "normal"

    cfg = {}
    try:
        cfg = attention_settings()
    except Exception:
        cfg = {}

    # 1. Desktop notification — same helper as the push_notification tool
    #    and job_monitor.announce().
    if cfg.get("desktop", True):
        try:
            from delfin.agent import notify as _notify
            _notify.send_notification(
                f"DELFIN: {title}", detail[:400], urgency=urgency)
        except Exception:
            pass

    # 2. Webhook — reuse the job monitor's URL setting + helper verbatim
    #    (agent.job_monitor.webhook_url via send_remote_trigger).
    try:
        from delfin.agent.job_monitor import monitor_settings
        webhook_url = monitor_settings().get("webhook_url", "")
    except Exception:
        webhook_url = ""
    if webhook_url:
        try:
            from delfin.agent import notify as _notify
            _notify.send_remote_trigger(
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
        except Exception:
            pass

    # 3. User shell hook (agent.attention.notify_command).
    notify_command = str(cfg.get("notify_command", "") or "")
    if notify_command:
        try:
            _run_notify_command(notify_command, title)
        except Exception:
            pass


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
) -> str:
    """Record a durable attention event and fan it out out-of-band.

    Returns the event id. The append is the contract (raises on disk
    failure); every delivery transport is best-effort and never raises.
    """
    if kind not in ATTENTION_KINDS:
        raise ValueError(
            f"unknown attention kind {kind!r}; expected one of "
            f"{ATTENTION_KINDS}")
    now = time.time()
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
    }
    with _FILE_LOCK:
        events = _load_events()
        events.append(event)
        _write_events(_prune(events, now))
    _deliver(event)
    return event["id"]


def list_pending(kind: Optional[str] = None) -> list[dict]:
    """Unresolved events (oldest first) — what the user still has to act
    on, across restarts. Optionally filtered by kind."""
    with _FILE_LOCK:
        events = _load_events()
    pending = [
        dict(ev) for ev in events
        if ev.get("status") == "pending"
        and (kind is None or ev.get("kind") == kind)
    ]
    pending.sort(key=lambda ev: float(ev.get("created_at") or 0.0))
    return pending


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
    with _FILE_LOCK:
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


def drain_resolved(session_id: str, *, include_unrouted: bool = True) -> list[dict]:
    """Answered events for ``session_id`` — exactly once.

    The ``acknowledged`` flag is persisted in the inbox file, so the
    exactly-once contract survives restarts (``bash_jobs.
    drain_finished_events`` pattern): a new process only sees what no
    previous process drained. Events emitted without a session id (e.g.
    the confirm broker has no session context) are included by default —
    routing them to the first session that drains beats losing them."""
    sid = str(session_id or "")
    out: list[dict] = []
    with _FILE_LOCK:
        events = _load_events()
        changed = False
        for ev in events:
            if ev.get("status") != "resolved" or ev.get("acknowledged"):
                continue
            ev_sid = str(ev.get("session_id") or "")
            if ev_sid != sid and not (include_unrouted and not ev_sid):
                continue
            ev["acknowledged"] = True
            changed = True
            out.append(dict(ev))
        if changed:
            _write_events(events)
    out.sort(key=lambda ev: float(ev.get("created_at") or 0.0))
    return out


# ---------------------------------------------------------------------------
# User surface (dashboard /attention command) — pure helpers, never raise
# ---------------------------------------------------------------------------

# Display order: kinds that block the agent first, informational last.
_KIND_ORDER = (
    "confirm_pending", "question_pending", "plan_pending",
    "run_failed", "run_finished", "budget_warning",
)
_KIND_LABELS = {
    "confirm_pending": "Confirmations waiting",
    "question_pending": "Questions waiting",
    "plan_pending": "Plans awaiting approval",
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


def render_inbox(kind: Optional[str] = None, *, limit: int = 50) -> str:
    """Compact, grouped, id-visible view of the pending inbox.

    Groups by kind (blocking kinds first), oldest first inside each
    group, at most ``limit`` items total (a footer notes the overflow).
    Never raises — a broken inbox renders as an error line."""
    try:
        pending = list_pending(kind)
    except Exception as exc:
        return f"Attention inbox unavailable: {exc}"
    try:
        scope = f" of kind {kind!r}" if kind else ""
        if not pending:
            return f"Attention inbox: no pending items{scope}."
        now = time.time()
        try:
            cap = max(1, int(limit))
        except (TypeError, ValueError):
            cap = 50
        urgent = sum(
            1 for ev in pending if ev.get("kind") in _URGENT_KINDS)
        lines = [
            f"Attention inbox — {len(pending)} pending{scope}"
            + (f" ({urgent} blocking the agent)" if urgent else "")]
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
        lines.append(
            "Answer: /attention answer <id> <text>   "
            "Dismiss: /attention dismiss <id|all>")
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


def clear_all(kind: Optional[str] = None) -> dict:
    """Dismiss every pending item (optionally one kind). Same
    ``acknowledged=True`` semantics as ``dismiss_item``. Never raises;
    returns ``{"ok": True, "cleared": n}`` on success."""
    try:
        cleared = 0
        for ev in list_pending(kind):
            if resolve(str(ev.get("id", "")), answer=None,
                       acknowledged=True):
                cleared += 1
        return {"ok": True, "cleared": cleared}
    except Exception as exc:
        return {"ok": False, "error": str(exc), "cleared": 0}


__all__ = [
    "ATTENTION_KINDS",
    "answer_item",
    "attention_settings",
    "clear_all",
    "dismiss_item",
    "drain_resolved",
    "emit_attention",
    "list_pending",
    "render_inbox",
    "resolve",
]
