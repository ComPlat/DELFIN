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


def collect(workspace: Any, *, now: Optional[float] = None) -> dict:
    """Everything still out for ``workspace``, grouped by kind."""
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
            if entry.get("kind") == "bash":
                # A running one already has its shell row; a finished one is
                # reported to the agent, not waited for.
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
            if not _sa._entry_owned_by_us(entry):
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
        for entry in _sch.get_scheduler().list_entries():
            if getattr(entry, "disabled", False):
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


def render_html(view: dict, *, now: Optional[float] = None) -> str:
    """The Background panel; "" when nothing is out, so the panel hides."""
    now = time.time() if now is None else float(now)
    rows: list[tuple[str, str, str]] = []
    for row in view.get("shells", []):
        rows.append(("Shell", row["label"] or row["id"],
                     f"running · {_duration(now - row['since'])}"))
    for row in view.get("watches", []):
        kind = row["kind"].upper() if row["kind"] in ("ci", "slurm") else row["kind"]
        rows.append((f"Watch · {kind}", row["label"] or row["id"],
                     f"{row['state'].lower()} · {_duration(now - row['since'])}"))
    for row in view.get("agents", []):
        detail = f"running · {_duration(now - row['since'])}"
        if row.get("last"):
            detail += f" · {row['last']}"
        rows.append(("Agent", row["label"], detail))
    for row in view.get("wakeups", []):
        label = "Loop" if row["kind"] == "interval" else "Wake-up"
        when = (f"in {_duration(row['at'] - now)}" if row["at"] > now else "due")
        rows.append((label, row["label"] or row["id"], when))
    if not rows:
        return ""
    body = "".join(
        "<div style='display:flex; gap:8px; align-items:baseline;'>"
        f"<span style='min-width:7.5em; color:#37474f; font-weight:600;'>{html.escape(kind)}</span>"
        f"<span style='flex:1; min-width:0; overflow-wrap:anywhere;'>{html.escape(label)}</span>"
        f"<span style='color:#78909c; white-space:nowrap;'>{html.escape(detail)}</span>"
        "</div>"
        for kind, label, detail in rows)
    return ("<div style='font-size:11px; color:#546e7a; display:flex; "
            "flex-direction:column; gap:2px;'>"
            f"<b>Background · {len(rows)} running</b>{body}</div>")


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
