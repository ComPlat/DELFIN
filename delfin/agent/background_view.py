"""What the agent has running in the background, in one view.

The agent starts work that outlives a tool call: shells, watched cluster
jobs and CI runs, scheduled wake-ups, sub-agents. Each kept its own records,
and nothing on screen said what was still out -- the subagent panel even
kept finished runs, while a background test suite or a pending CI watch did
not show at all. Claude Code shows one list of running tasks; this is that
list for the dashboard.

Reads local state only (job registries, watch files, the scheduler, the
sub-agent registry). It never asks SLURM or GitHub: the panel refreshes every
few seconds, and the watches already ask those, throttled, on their own.
"""

from __future__ import annotations

import html
import json
import time
from pathlib import Path
from typing import Any, Optional


def collect(workspace: Any, *, now: Optional[float] = None,
            session_id: Optional[str] = None) -> dict:
    """Everything still out for ``workspace``, grouped by kind.

    With ``session_id``, watches and wake-ups another session owns are left
    out: the panel belongs to one conversation, and several can share a
    workspace."""

    def _theirs(owner: Any) -> bool:
        # Strict: work nobody owns predates ownership and is no session's.
        return bool(session_id) and str(owner or "") != session_id

    now = time.time() if now is None else float(now)
    ws = str(workspace or "")
    view: dict[str, list[dict]] = {
        "shells": [], "watches": [], "agents": [], "wakeups": []}

    try:
        from . import bash_jobs as _bj
        live = {str(j.get("job_id")) for j in (_bj.live_jobs() or [])}
        records = (_bj._load_registry_file(ws).get("jobs") or {}) if ws else {}
        for jid, rec in records.items():
            if str(jid) not in live:
                continue
            rec = rec or {}
            if _theirs(rec.get("session_id")):
                continue
            view["shells"].append({
                "id": str(jid),
                "label": str(rec.get("description") or rec.get("command") or "")[:90],
                "since": float(rec.get("started_at") or now),
            })
    except Exception:
        pass

    try:
        from . import job_monitor as _jm
        watched = (_jm.load_watched(_jm._agent_watch_path(ws)).get("jobs") or {}) if ws else {}
        for jid, entry in watched.items():
            entry = entry or {}
            if entry.get("kind") == "bash" or _theirs(entry.get("session_id")):
                # A running bash job already has its shell row; a finished
                # one is reported to the agent, not waited for.
                continue
            view["watches"].append({
                "id": str(jid),
                "kind": str(entry.get("kind") or "job"),
                "label": str(entry.get("description") or "")[:90],
                "state": str(entry.get("last_state") or "waiting"),
                "since": float(entry.get("added_at") or now),
            })
    except Exception:
        pass

    try:
        from . import subagents as _sa
        for sa_id, entry in (_sa.read_running() or {}).items():
            entry = entry or {}
            if (not _sa._entry_owned_by_us(entry)
                    or _theirs(entry.get("owner_session"))):
                continue      # another session's delegate
            view["agents"].append({
                "id": str(sa_id),
                "label": (f"{entry.get('type') or 'agent'} · "
                          f"{str(entry.get('description') or '')[:70]}"),
                "since": float(entry.get("started_at") or now),
                "last": str(entry.get("last_action") or "")[:60],
            })
    except Exception:
        pass

    try:
        from . import scheduler as _sch
        _scheduler = _sch.get_scheduler()
        for entry in _scheduler.list_entries():
            if getattr(entry, "disabled", False):
                continue
            if session_id and _theirs(_scheduler.owner_of(entry.id)):
                continue
            view["wakeups"].append({
                "id": str(entry.id),
                "kind": str(getattr(entry, "kind", "once")),
                "label": str(entry.reason or entry.prompt or "")[:90],
                "at": float(getattr(entry, "next_fire_at", 0) or 0),
            })
    except Exception:
        pass

    for group in ("shells", "watches", "agents"):
        view[group].sort(key=lambda row: row["since"])
    view["wakeups"].sort(key=lambda row: row["at"])
    return view


def _duration(seconds: float) -> str:
    s = max(0, int(seconds))
    if s < 60:
        return f"{s}s"
    if s < 3600:
        return f"{s // 60}m{s % 60:02d}s"
    return f"{s // 3600}h{(s % 3600) // 60:02d}m"


# What the × on a row does, said where the user can read it first.
_STOP_TIPS = {
    "shells": "Stop this command",
    "watches": "Stop watching — the job itself keeps running",
    "agents": "Stop this background agent",
    "wakeups": "Cancel this wake-up",
}


def rows(view: dict, *, now: Optional[float] = None) -> list[dict]:
    """One row per item: ``kind``, ``label``, ``detail``; ``group`` and ``id``
    name it for :func:`cancel`, and ``tip`` says what stopping it does."""
    now = time.time() if now is None else float(now)
    out: list[dict] = []

    def _add(group: str, item_id: str, kind: str, label: str, detail: str):
        out.append({"group": group, "id": str(item_id), "kind": kind,
                    "label": label, "detail": detail,
                    "tip": _STOP_TIPS[group]})

    for row in view.get("shells", []):
        _add("shells", row["id"], "Shell", row["label"] or row["id"],
             f"running · {_duration(now - row['since'])}")
    for row in view.get("watches", []):
        kind = row["kind"].upper() if row["kind"] in ("ci", "slurm") else row["kind"]
        _add("watches", row["id"], f"Watch · {kind}", row["label"] or row["id"],
             f"{row['state'].lower()} · {_duration(now - row['since'])}")
    for row in view.get("agents", []):
        detail = f"running · {_duration(now - row['since'])}"
        if row.get("last"):
            detail += f" · {row['last']}"
        _add("agents", row["id"], "Agent", row["label"], detail)
    for row in view.get("wakeups", []):
        label = "Loop" if row["kind"] == "interval" else "Wake-up"
        when = (f"in {_duration(row['at'] - now)}" if row["at"] > now else "due")
        _add("wakeups", row["id"], label, row["label"] or row["id"], when)
    return out


def row_html(row: dict) -> str:
    return ("<div style='display:flex; gap:8px; align-items:baseline; "
            "font-size:11px; color:#546e7a;'>"
            f"<span style='min-width:7.5em; color:#37474f; font-weight:600;'>"
            f"{html.escape(row['kind'])}</span>"
            f"<span style='flex:1; min-width:0; overflow-wrap:anywhere;'>"
            f"{html.escape(row['label'])}</span>"
            f"<span style='color:#78909c; white-space:nowrap;'>"
            f"{html.escape(row['detail'])}</span></div>")


def header_html(count: int) -> str:
    return (f"<b style='font-size:11px; color:#546e7a;'>Background · "
            f"{count} running</b>") if count else ""


def render_html(view: dict, *, now: Optional[float] = None) -> str:
    """The Background panel; "" when nothing is out, so the panel hides."""
    items = rows(view, now=now)
    if not items:
        return ""
    return ("<div style='display:flex; flex-direction:column; gap:2px;'>"
            + header_html(len(items))
            + "".join(row_html(row) for row in items) + "</div>")


def cancel(workspace: Any, group: str, item_id: str) -> str:
    """Stop one listed item, and say what happened.

    A shell is terminated, a sub-agent is asked to stop, a wake-up is
    deleted. A watch only stops being watched: the cluster job or CI run it
    follows is not the dashboard's to end.
    """
    item_id = str(item_id or "")
    try:
        if group == "shells":
            from . import bash_jobs as _bj
            _ok, message = _bj.get_registry().kill(item_id)
            return f"Shell {item_id}: {message}"
        if group == "watches":
            from . import job_monitor as _jm
            if _jm.unwatch_agent_job(str(workspace or ""), item_id):
                return (f"Stopped watching {item_id}; the job itself keeps "
                        "running.")
            return f"{item_id} was not being watched."
        if group == "agents":
            from . import subagents as _sa
            if _sa.cancel_background(item_id):
                return f"Background agent {item_id} is stopping."
            return f"Background agent {item_id} is not running."
        if group == "wakeups":
            from . import scheduler as _sch
            if _sch.get_scheduler().delete(item_id):
                return f"Wake-up {item_id} cancelled."
            return f"Wake-up {item_id} was not scheduled."
    except Exception as exc:
        return f"Could not stop {item_id}: {exc}"
    return f"Nothing to stop for {group} {item_id}."


def finished_background_agents(already: set) -> list[dict]:
    """Background sub-agents of this session that finished and have not
    woken the agent yet, as wake entries. Adds each id to ``already``.

    Peeks at the pending-report markers without claiming them: the turn the
    wake starts drains the report itself, exactly once.
    """
    out: list[dict] = []
    try:
        from . import subagents as _sa
        running = set((_sa.read_running() or {}).keys())
        for f in sorted(Path(_sa._PENDING_DIR).glob("*.json")):
            try:
                rec = json.loads(f.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not isinstance(rec, dict) or not _sa._entry_owned_by_us(rec):
                continue
            sa_id = str(rec.get("sa_id") or f.stem)
            if sa_id in running or sa_id in already:
                continue
            already.add(sa_id)
            out.append({
                "kind": "agent", "job_id": sa_id, "state": "FINISHED",
                "description": (f"{rec.get('type') or 'agent'} · "
                                f"{str(rec.get('description') or '')[:70]}"),
            })
    except Exception:
        pass
    return out
