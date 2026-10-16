"""What needs the user, from every DELFIN process, as the windows list it.

Two records already hold it, written by every DELFIN process of the user
(dashboards, the REPL, the job monitor daemon, the scheduler, sessions of
this service), on a cluster in the home every login node shares:

- the attention inbox (``delfin.agent.attention``): questions a session
  parked when nobody answered, confirmations, failed and finished runs,
  budget warnings, security alerts;
- the job monitor's findings (``job_monitor.load_findings``): a failed
  job's state, error signatures, the session its diagnosis was saved as,
  a summary and a prepared fix.

``snapshot`` turns them into the items the agent view shows under NEEDS
YOU (what blocks, what failed) and NOTICES (what happened). The inbox
decides what is open: a finding shows only while its job's ``run_failed``
entry is pending, and enriches it. An entry of a session this service
runs is left out -- that session's own rows say it. The files are read
through the modules' own readers and paths (tolerant of torn lines, and
redirected under scratch state as the rest of DELFIN's state is).

Nothing here writes, except ``answer`` and ``dismiss``, which go through
the inbox's own helpers (its lock, its rules for what may be cleared).

What goes out is masked as a session's events are (``redact.Redactor``):
the records hold what other processes wrote -- a confirmation's command
line, a job's output -- and every window, and every agent process that
watches the service, receives the items.
"""

from __future__ import annotations

import os
import re
import time
from pathlib import Path
from typing import Any

#: Inbox kinds -> (the window's kind, its group, whether it blocks someone).
KINDS: dict[str, tuple[str, str, bool]] = {
    "question_pending": ("question", "needs", True),
    "confirm_pending": ("confirm", "needs", True),
    "plan_pending": ("plan", "needs", True),
    "run_failed": ("failed", "needs", False),
    "security_alert": ("alert", "needs", False),
    "run_finished": ("finished", "notices", False),
    "budget_warning": ("budget", "notices", False),
}

#: Findings older than this have no pending inbox entry left (a failure
#: notice expires after seven days): not read.
FINDINGS_WINDOW_S = 8 * 86400.0

#: The most items sent: blocking ones and failures first, then the newest.
MAX_ITEMS = 300

_JOB_TITLE = re.compile(r"^Job\s+(\S+)\s+(failed|finished)\b", re.IGNORECASE)


def _inbox_path() -> Path:
    from .. import attention
    return attention._inbox_path()


def _findings_path() -> Path:
    from .. import job_monitor
    return Path(job_monitor._FINDINGS_PATH)


def fingerprint() -> tuple:
    """What changes when either record does: size and modification time."""
    out = []
    for path in (_inbox_path(), _findings_path()):
        try:
            st = path.stat()
            out.append((st.st_size, st.st_mtime_ns))
        except OSError:
            out.append(None)
    return tuple(out)


def _pending_events(now: float) -> list[dict]:
    """The inbox's open entries, expired notices left out (its own view)."""
    from .. import attention
    try:
        events = attention._visible(attention._load_events(), now)
    except Exception:
        return []
    return [e for e in events if e.get("status") == "pending"]


def _latest_findings(now: float) -> dict[str, dict]:
    """The newest finding per job id, of the last few days."""
    from .. import job_monitor
    try:
        found = job_monitor.load_findings(since=now - FINDINGS_WINDOW_S, path=_findings_path())
    except Exception:
        return {}
    out: dict[str, dict] = {}
    for rec in found:
        jid = str(rec.get("job_id") or "")
        if jid and float(rec.get("ts") or 0) >= float(out.get(jid, {}).get("ts") or 0):
            out[jid] = rec
    return out


def _first_line(text: str, limit: int = 160) -> str:
    line = next((l.strip() for l in str(text or "").splitlines() if l.strip()), "")
    return line if len(line) <= limit else line[: limit - 1] + "…"


class Titles:
    """What a session is called, by its id: an open session's presence
    record says it (and where it runs); a saved one's record in the store.
    A saved record is read again only when its file changed -- it is the
    whole conversation, and the feed is built every few seconds; one the
    store does not have yet (a dashboard's first turn still running) is
    asked again."""

    def __init__(self) -> None:
        self._saved: dict[str, tuple[Any, tuple[str, bool]]] = {}

    def presence(self) -> dict[str, dict]:
        try:
            from .. import session_presence
            return {str(r.get("session_id") or ""): r for r in session_presence.open_sessions()
                    if r.get("session_id")}
        except Exception:
            return {}

    def saved(self, sid: str) -> tuple[str, bool]:
        """(title, whether the store has the session)."""
        stamp = _record_stamp(sid)
        cached = self._saved.get(sid)
        if cached and stamp is not None and cached[0] == stamp:
            return cached[1]
        try:
            from .. import session_store
            data = session_store.load_session(sid)
        except Exception:
            data = None
        if not data:
            self._saved.pop(sid, None)
            return "", False
        found = (str(data.get("title") or data.get("session_title") or ""), True)
        if stamp is not None:
            self._saved[sid] = (stamp, found)
        return found


def _record_stamp(sid: str):
    """When a saved session's record last changed (renamed, saved again),
    or None when that cannot be told: then it is read each time."""
    try:
        from .. import session_store
        st = (Path(session_store._SESSIONS_DIR) / f"{sid}.json").stat()
        return (st.st_size, st.st_mtime_ns)
    except Exception:
        return None


def _held_elsewhere(sid: str) -> dict | None:
    """Who else writes the session, if anyone (the store's own lock)."""
    try:
        from .. import session_store
        return session_store.session_lock_holder(sid)
    except Exception:
        return None


def prepared_fix(finding: dict | None) -> dict:
    """The bounded fix a finding's assessment prepared (job_fix's payload
    has one under "fix" only when it could prepare one), or {}."""
    assessed = (finding or {}).get("fix")
    fix = assessed.get("fix") if isinstance(assessed, dict) else None
    return fix if isinstance(fix, dict) else {}


def _job_text(finding: dict | None, state: str) -> str:
    parts = [state] if state else []
    if finding:
        if finding.get("diagnosis_session"):
            parts.append("diagnosis ready")
        if prepared_fix(finding):
            parts.append("a fix is prepared")
    return " · ".join(parts) or "failed"


def item_of(event: dict, *, findings: dict[str, dict], presence: dict[str, dict],
            titles: Titles) -> dict:
    """One inbox entry as a row of the agent view."""
    kind, group, blocking = KINDS.get(str(event.get("kind") or ""), ("notice", "notices", False))
    title = str(event.get("title") or "").strip()
    detail = str(event.get("detail") or "").strip()
    sid = str(event.get("session_id") or "")
    workspace = str(event.get("workspace") or "")
    host = str(event.get("host") or "")
    item: dict[str, Any] = {
        "id": str(event.get("id") or ""), "kind": kind, "group": group, "blocking": blocking,
        "name": title, "text": _first_line(detail) or title, "detail": detail,
        "workspace": workspace, "host": host, "sessionId": sid,
        "createdAt": float(event.get("created_at") or 0.0),
        "options": [str(o) for o in (event.get("options") or [])],
    }
    job = _JOB_TITLE.match(title)
    if job and kind in ("failed", "finished"):
        jid = job.group(1)
        finding = findings.get(jid) if kind == "failed" else None
        state = ""
        m = re.search(r"\(([A-Z_]+)\)", title)
        if m:
            state = m.group(1)
        if finding:
            state = str(finding.get("state") or state)
            workspace = workspace or str(finding.get("folder") or "")
            item["job"] = {
                "id": jid, "state": state, "folder": str(finding.get("folder") or ""),
                "signatures": [str(s) for s in finding.get("signatures") or []],
                "diagnosisSession": str(finding.get("diagnosis_session") or ""),
                "summary": str(finding.get("summary") or ""),
                # job_fix's assessment: fix_class, recommendation, and under
                # "fix" the prepared fix, if any (summary_line, diff).
                "fix": finding.get("fix") if isinstance(finding.get("fix"), dict) else {},
            }
        else:
            item["job"] = {"id": jid, "state": state, "folder": workspace}
        item["name"] = f"job {jid} {job.group(2).lower()}"
        item["text"] = (_job_text(finding, state) if kind == "failed"
                        else _first_line(detail) or state or "finished")
        item["workspace"] = workspace
    elif kind == "question":
        question = detail or re.sub(r"^Question waiting:\s*", "", title)
        item["text"] = "question: " + _first_line(question, 140)
        item["detail"] = question
    elif kind == "confirm":
        item["text"] = re.sub(r"^Confirmation required:\s*", "confirm: ", title) or "a confirmation"
    # Whose it is, and where that runs: open in another process (its
    # presence record, or the store's lock) -- a window opens a copy -- and
    # whether the store has it at all (else it cannot be opened).
    if sid:
        saved_title, exists = titles.saved(sid)
        item["sessionSaved"] = exists
        open_record = presence.get(sid)
        if open_record:
            item["sessionTitle"] = str(open_record.get("title") or "") or saved_title
            item["host"] = item["host"] or str(open_record.get("host") or "")
            item["sessionOpen"] = True
        else:
            item["sessionTitle"] = saved_title
            holder = _held_elsewhere(sid) if exists else None
            item["sessionOpen"] = bool(holder)
            if holder and not item["host"]:
                item["host"] = str(holder.get("host") or "")
        if kind in ("question", "confirm", "plan") and item.get("sessionTitle"):
            item["name"] = item["sessionTitle"]
    if kind in ("question", "confirm", "plan") and not item.get("sessionTitle"):
        # Named for where it was asked: the inbox's title repeats the text.
        word = {"question": "a question", "confirm": "a confirmation", "plan": "a plan"}[kind]
        folder = Path(workspace).name if workspace else ""
        item["name"] = f"{word} in {folder}" if folder else word
    if not item["name"]:
        item["name"] = kind
    return item


def snapshot(own_sessions: set[str] | frozenset[str] = frozenset(), *, now: float | None = None,
             titles: Titles | None = None) -> list[dict]:
    """Every open item, the agent view's order: what blocks, then what
    failed, then notices; the newest first in each."""
    now = time.time() if now is None else now
    titles = titles or Titles()
    findings = _latest_findings(now)
    presence = titles.presence()
    items = []
    for event in _pending_events(now):
        if str(event.get("session_id") or "") in own_sessions:
            continue
        try:
            items.append(item_of(event, findings=findings, presence=presence, titles=titles))
        except Exception:
            continue
    rank = {True: 0, False: 1}
    items.sort(key=lambda it: (0 if it["group"] == "needs" else 1, rank[it["blocking"]],
                               -it["createdAt"]))
    redactor = _redactor()
    return [redactor.event(it, _KEPT) for it in items[:MAX_ITEMS]]


#: Fields of an item that are no text anyone wrote: never masked.
_KEPT = frozenset({"id", "kind", "group", "blocking", "sessionId", "createdAt",
                   "sessionOpen", "sessionSaved"})
_REDACTOR = None


def _redactor():
    global _REDACTOR
    if _REDACTOR is None:
        from .redact import Redactor
        _REDACTOR = Redactor()
    return _REDACTOR


def answer(item_id: str, text: str) -> dict:
    """The user's answer to a parked question: the next session in its
    workspace gets it (the inbox's own exactly-once hand-over)."""
    from .. import attention
    return attention.answer_item(item_id, text)


def dismiss(item_id: str) -> dict:
    """Cleared without an answer; never handed to a model."""
    from .. import attention
    return attention.dismiss_item(item_id)


def this_host() -> str:
    try:
        import socket
        return socket.gethostname()
    except Exception:
        return os.environ.get("HOSTNAME", "")
