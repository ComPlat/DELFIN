"""Append-only audit log for every code-modifying KIT-agent action.

Every successful or denied call to write_file / edit_file / multi_edit /
notebook_edit / bash / bash_background / remember_permission(_bundle)
appends one JSON line to ``~/.delfin/audit.log``. Records carry:

* ``ts`` — ISO-8601 UTC timestamp
* ``session_id`` — engine session UUID, if any
* ``tool`` — tool name (lowercase, no MCP prefix)
* ``path`` / ``command`` — the target (file path / shell text)
* ``mode`` — KIT mode active at the time
* ``decision`` — ``auto`` (gate passed without a confirm) /
  ``confirmed`` (user clicked Erlauben) / ``denied`` (gate blocked or
  user clicked Ablehnen) / ``error`` (tool raised)
* ``persistence`` — for remember_* calls, what got written and where
* ``diff_lines_added`` / ``diff_lines_removed`` — for write/edit tools
* ``pid`` — the calling process pid (helps trace bash_background jobs)

Lines are written atomically (single ``write`` call with a final
``\\n``) so concurrent writes from multiple worker threads don't
interleave inside one record. Reads should treat each line
independently — corrupt lines (e.g. from a kill -9 mid-write) are
rare but possible; tooling should ``json.loads`` line-by-line and
skip what doesn't parse.

The log is ROTATED weekly: when the first line of a calendar week is
appended and the existing log is older than the week boundary, the
old log is renamed to ``audit-YYYY-Www.log`` and a fresh log starts.
This keeps the active file small enough to ``tail`` without pain
while preserving full history under predictable filenames.
"""

from __future__ import annotations

import json
import os
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


_LOG_LOCK = threading.Lock()


def _default_log_path() -> Path:
    return Path.home() / ".delfin" / "audit.log"


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _iso(ts: datetime) -> str:
    return ts.strftime("%Y-%m-%dT%H:%M:%SZ")


def _isocalendar_key(ts: datetime) -> str:
    """Return 'YYYY-Www' (ISO week) for rotation labelling."""
    iso = ts.isocalendar()
    return f"{iso.year}-W{iso.week:02d}"


def _rotate_if_needed(path: Path, now: datetime) -> None:
    """Move ``path`` aside when the active log spans across an ISO week.

    Called inside the lock just before each append. If the existing
    log's first-line timestamp belongs to an earlier ISO week than
    ``now``, we rename it to ``audit-YYYY-Www.log`` and start fresh.
    Cheap: read only the first line.
    """
    if not path.exists() or path.stat().st_size == 0:
        return
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            first_line = fh.readline()
        first = json.loads(first_line)
        first_ts = datetime.fromisoformat(first["ts"].rstrip("Z"))
        first_ts = first_ts.replace(tzinfo=timezone.utc)
    except Exception:
        return  # Corrupt header — leave it alone, append still works.
    if _isocalendar_key(first_ts) == _isocalendar_key(now):
        return
    rotated = path.with_name(f"audit-{_isocalendar_key(first_ts)}.log")
    try:
        path.rename(rotated)
    except OSError:
        pass


def append(record: dict, *, log_path: Optional[Path] = None) -> None:
    """Append one JSON record to the audit log.

    ``record`` is shallow-copied and decorated with the timestamp; the
    caller owns the rest of the schema. Failures are silent — the
    audit log must never block or break the tool path it observes.
    """
    path = log_path or _default_log_path()
    now = _now()
    enriched = {"ts": _iso(now), **record}
    line = json.dumps(enriched, ensure_ascii=False) + "\n"
    try:
        with _LOG_LOCK:
            path.parent.mkdir(parents=True, exist_ok=True)
            _rotate_if_needed(path, now)
            with path.open("a", encoding="utf-8") as fh:
                fh.write(line)
    except Exception:
        # Audit must not crash the agent. If the user removed the home
        # dir or ran out of disk, the action proceeds without logging.
        pass


def read_last_n(n: int = 100, *, log_path: Optional[Path] = None) -> list[dict]:
    """Return the last ``n`` parseable records (utility for the dashboard)."""
    path = log_path or _default_log_path()
    if not path.exists():
        return []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except Exception:
        return []
    out: list[dict] = []
    for line in lines[-(n * 2):]:        # over-read a little for parse errors
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except Exception:
            continue
    return out[-n:]


def _diff_line_counts(old: str, new: str) -> tuple[int, int]:
    """Quick added / removed line counts for the audit record.

    Not a real diff — just compares line sets to avoid bringing
    difflib into a hot path. Good enough for "approximately how big
    was this change" telemetry.
    """
    old_lines = old.splitlines() if old else []
    new_lines = new.splitlines() if new else []
    added = max(0, len(new_lines) - len(old_lines))
    removed = max(0, len(old_lines) - len(new_lines))
    # If the line counts match, count the actually-different lines.
    if added == 0 and removed == 0 and old_lines != new_lines:
        diffs = sum(1 for a, b in zip(old_lines, new_lines) if a != b)
        return diffs, diffs
    return added, removed


def make_record(
    *,
    tool: str,
    decision: str,
    mode: str = "",
    path: str = "",
    command: str = "",
    diff_lines_added: int = 0,
    diff_lines_removed: int = 0,
    persistence: Optional[dict] = None,
    session_id: str = "",
    extra: Optional[dict[str, Any]] = None,
) -> dict:
    """Construct a record dict with the standard fields."""
    rec: dict[str, Any] = {
        "session_id": session_id,
        "tool": tool,
        "decision": decision,
        "mode": mode,
        "pid": os.getpid(),
    }
    if path:
        rec["path"] = path
    if command:
        rec["command"] = command[:500]
    if diff_lines_added or diff_lines_removed:
        rec["diff_lines_added"] = diff_lines_added
        rec["diff_lines_removed"] = diff_lines_removed
    if persistence:
        rec["persistence"] = persistence
    if extra:
        rec.update(extra)
    return rec
