"""Append-only retrospective log of agent failures.

Records `tool / command / error` triplets so the agent can:
- detect "I've broken this 3 times" patterns within or across sessions
- offer the user a suggestion to escalate, plan, or stop retrying
- surface a per-tool failure leaderboard via ``/failures``

The log lives at ``~/.delfin/failure_log.jsonl`` and is trimmed to the
most recent ``_MAX_LINES`` entries to keep it cheap to read on every
turn. Best-effort writes — silently skip if the home dir isn't
writable (the agent must never crash because logging failed).
"""

from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path


_LOG_PATH = Path.home() / ".delfin" / "failure_log.jsonl"
_MAX_LINES = 2000

_RECENT_WINDOW_S = 7 * 86_400  # 7 days = considered recent


def record_failure(
    tool: str,
    command: str,
    error: str,
    *,
    session_id: str = "",
    path: Path | None = None,
) -> None:
    """Append one failure record. Best-effort — never raises."""
    p = path or _LOG_PATH
    try:
        from .state_paths import ensure_dir, open_append, secure_file
        from .state_paths import write_text as _write_secure
        ensure_dir(p.parent)
        rec = {
            "ts": time.time(),
            "tool": (tool or "")[:80],
            "command": (command or "")[:400],
            "error": (error or "")[:400],
            "session_id": session_id or "",
        }
        # Created 0600 by ``open_append``: a chmod after the first write
        # leaves the file group-readable for the length of that write.
        with open_append(p) as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        secure_file(p)          # tighten a file that already existed
        # Trim if large; tolerant of trimming failure.
        try:
            if p.stat().st_size > 400_000:
                lines = p.read_text(encoding="utf-8").splitlines()
                if len(lines) > _MAX_LINES:
                    tail = lines[-_MAX_LINES:]
                    _write_secure(p, "\n".join(tail) + "\n")
        except OSError:
            pass
    except OSError:
        return


def read_failures(
    *,
    last_n: int | None = None,
    path: Path | None = None,
) -> list[dict]:
    """Load the failure log. Newest entries are at the end."""
    p = path or _LOG_PATH
    if not p.exists():
        return []
    out: list[dict] = []
    try:
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    except OSError:
        return []
    if last_n is not None and last_n > 0:
        return out[-last_n:]
    return out


def top_recurring(
    *,
    min_count: int = 3,
    within_s: float = _RECENT_WINDOW_S,
    path: Path | None = None,
) -> list[dict]:
    """Return the top recurring (tool, error-shape) pairs.

    Error-shape = first 40 chars of the error (low-cardinality clustering).
    Only pairs that occurred ``>= min_count`` times in the last
    ``within_s`` seconds are returned, newest signature first.
    """
    now = time.time()
    records = [
        r for r in read_failures(path=path)
        if (now - float(r.get("ts") or 0)) <= within_s
    ]
    if not records:
        return []
    shapes: Counter = Counter()
    latest_by_shape: dict[tuple[str, str], dict] = {}
    for r in records:
        shape = (r.get("tool") or "", (r.get("error") or "")[:40])
        shapes[shape] += 1
        latest_by_shape[shape] = r
    out: list[dict] = []
    for (tool, shape), n in shapes.most_common():
        if n < min_count:
            continue
        last = latest_by_shape[(tool, shape)]
        out.append({
            "tool": tool,
            "error_shape": shape,
            "count": n,
            "last_seen": last.get("ts"),
            "last_command": last.get("command", ""),
            "last_error": last.get("error", ""),
        })
    return out


def detect_repeat_for_current_task(
    tool: str,
    command: str,
    error: str,
    *,
    threshold: int = 3,
    within_s: float = 3600,
    path: Path | None = None,
) -> int:
    """Count how many times THIS specific (tool, command-shape, error-
    shape) has fired in the recent window. Used by the engine to short-
    circuit retry loops the user has already lost patience for."""
    now = time.time()
    cmd_shape = (command or "")[:80]
    err_shape = (error or "")[:80]
    count = 0
    for r in read_failures(path=path):
        if (now - float(r.get("ts") or 0)) > within_s:
            continue
        if r.get("tool") != tool:
            continue
        if (r.get("command") or "")[:80] != cmd_shape:
            continue
        if (r.get("error") or "")[:80] != err_shape:
            continue
        count += 1
    return count


__all__ = [
    "record_failure",
    "read_failures",
    "top_recurring",
    "detect_repeat_for_current_task",
]
