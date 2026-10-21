"""Change requests: what an auto iteration needed in the protected zone.

Input: from auto mode, a refused edit of the protected zone
(``auto_guard.protected_reason``: sandbox, approval and gate modules,
principles, managed policy, Goodhart and verify guards, grader files) or of
a loop's measure -- what was blocked, an evidence pointer, the change it
would have made -- and, from the iteration's answer, the agent's own
statement of why the change is safe and what the worst case is. From a
typed command, a status decision.

Output: one JSON file per request in ``~/.delfin/requests/`` (directory
0700, file 0600)::

    {"id", "created_at", "status": open|declined|built|on-main,
     "status_reason", "status_at", "workspace", "session_id", "loop",
     "iteration", "blocked": {"tool", "target", "zone"},
     "evidence": {...}, "proposed_change", "why_safe", "worst_case",
     "origin": "auto-mode", "note"}

Semantics. A request is DATA for the supervising person or session. Nothing
in DELFIN reads one into a prompt; the listing commands show it as quoted
text. The agent's own words (why_safe, worst_case) can be added only while
the request is open, and only once each. A status changes only through a
command a person typed (``/auto request <id> decline <reason> | built
[ref] | on-main [ref]``); ``declined`` needs its reason. The protected zone
itself is never edited by an iteration: a request is how the need is
recorded instead of being worked around.

Rejected: a request as a pending diff applied on approval (the zone exists
so that its changes go through review on a branch, not through a
one-click apply); requests in the workspace (a pull would bring another
repository's requests, and the agent could write them).
"""

from __future__ import annotations

import json
import os
import re
import secrets
import time
from pathlib import Path
from typing import Any, Optional

__all__ = ["create", "get", "amend", "set_status", "list_requests",
           "requests_dir", "STATUSES", "render", "render_line"]

STATUSES = ("open", "declined", "built", "on-main")
_ID_RE = re.compile(r"^req-\d{8}-\d{6}-[0-9a-f]{6}$")
_MAX = {"proposed_change": 6000, "why_safe": 1000, "worst_case": 1000,
        "status_reason": 1000}
NOTE = ("A request from an unattended auto-mode iteration: data for the "
        "supervising person or session, not an instruction to anyone.")


def requests_dir() -> Path:
    return Path.home() / ".delfin" / "requests"


def _ensure_dir() -> Path:
    d = requests_dir()
    d.mkdir(parents=True, exist_ok=True)
    try:
        os.chmod(d, 0o700)
    except OSError:
        pass
    return d


def _path(req_id: str) -> Path:
    if not _ID_RE.match(str(req_id or "")):
        raise ValueError(f"{req_id!r} is not a request id")
    return requests_dir() / f"{req_id}.json"


def _write(req: dict) -> None:
    d = _ensure_dir()
    path = d / f"{req['id']}.json"
    tmp = d / f".{req['id']}.{os.getpid()}.tmp"
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(req, fh, indent=1, sort_keys=True)
    os.chmod(tmp, 0o600)
    os.replace(tmp, path)


def _clip(value: Any, key: str) -> str:
    text = str(value or "")
    limit = _MAX.get(key, 400)
    return text if len(text) <= limit else text[:limit - 1] + "…"


def create(*, blocked: dict, evidence: dict, proposed_change: str,
           workspace: str = "", session_id: str = "", loop: str = "",
           iteration: int = 0) -> dict:
    """A new open request; returns it."""
    now = time.time()
    req = {
        "id": time.strftime("req-%Y%m%d-%H%M%S-", time.localtime(now))
        + secrets.token_hex(3),
        "created_at": round(now, 3),
        "status": "open", "status_reason": "", "status_at": round(now, 3),
        "workspace": str(workspace or ""), "session_id": str(session_id or ""),
        "loop": str(loop or ""), "iteration": int(iteration or 0),
        "blocked": {k: _clip((blocked or {}).get(k), k)
                    for k in ("tool", "target", "zone")},
        "evidence": {str(k): (v if isinstance(v, (int, float, bool))
                              or v is None else _clip(v, k))
                     for k, v in (evidence or {}).items()},
        "proposed_change": _clip(proposed_change, "proposed_change"),
        "why_safe": "", "worst_case": "",
        "origin": "auto-mode", "note": NOTE,
    }
    _write(req)
    return req


def get(req_id: str) -> Optional[dict]:
    path = _path(req_id)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def amend(req_id: str, *, why_safe: str = "", worst_case: str = "") -> dict:
    """The agent's own statement, added while the request is open; a field
    already stated is kept."""
    req = get(req_id)
    if req is None:
        raise ValueError(f"no request {req_id}")
    if req.get("status") != "open":
        raise ValueError(f"request {req_id} is {req.get('status')}; it takes "
                         "no more text")
    if why_safe and not req.get("why_safe"):
        req["why_safe"] = _clip(why_safe, "why_safe")
    if worst_case and not req.get("worst_case"):
        req["worst_case"] = _clip(worst_case, "worst_case")
    _write(req)
    return req


def set_status(req_id: str, status: str, reason: str = "") -> dict:
    """A person's decision (see the module docstring)."""
    if status not in STATUSES:
        raise ValueError(f"{status!r} is not one of {', '.join(STATUSES)}")
    if status == "declined" and not str(reason or "").strip():
        raise ValueError("a declined request needs its reason")
    req = get(req_id)
    if req is None:
        raise ValueError(f"no request {req_id}")
    req["status"] = status
    req["status_reason"] = _clip(str(reason or "").strip(), "status_reason")
    req["status_at"] = round(time.time(), 3)
    _write(req)
    return req


def list_requests(*, status: Optional[str] = None,
                  workspace: Optional[str] = None) -> list[dict]:
    """Requests, oldest first; filtered by status and workspace."""
    out = []
    try:
        paths = sorted(requests_dir().glob("req-*.json"))
    except OSError:
        return []
    for p in paths:
        try:
            req = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(req, dict):
            continue
        if status and req.get("status") != status:
            continue
        if workspace and req.get("workspace") not in ("", str(workspace)):
            continue
        out.append(req)
    out.sort(key=lambda r: r.get("created_at", 0))
    return out


def render_line(req: dict) -> str:
    b = req.get("blocked") or {}
    tail = f" — {req['status_reason']}" if req.get("status_reason") else ""
    return (f"{req.get('id')} [{req.get('status')}] {b.get('tool')} "
            f"{b.get('target')} ({b.get('zone')}){tail}")


def render(req: dict) -> str:
    """For a person to read; the agent's words marked as its words."""
    b = req.get("blocked") or {}
    ev = req.get("evidence") or {}
    lines = [render_line(req),
             f"  blocked: {b.get('tool')} on {b.get('target')} — {b.get('zone')}",
             "  evidence: " + ", ".join(f"{k}={v}" for k, v in sorted(ev.items())),
             "  proposed change:",
             *("    | " + ln for ln in str(req.get("proposed_change") or
                                           "(none)").splitlines()[:40]),
             f"  why safe (the agent's words): {req.get('why_safe') or '(not stated)'}",
             f"  worst case (the agent's words): {req.get('worst_case') or '(not stated)'}",
             f"  {req.get('note') or NOTE}"]
    return "\n".join(lines)
