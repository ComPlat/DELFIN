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


def append(
    record: dict,
    *,
    log_path: Optional[Path] = None,
    session_id: str = "",
) -> None:
    """Append one JSON record to the audit log.

    ``record`` is shallow-copied and decorated with the timestamp; the
    caller owns the rest of the schema. ``session_id`` (optional) stamps
    the record with the engine session when the record itself carries
    none — callers that already put ``session_id`` into the record keep
    their value. Failures are silent — the audit log must never block
    or break the tool path it observes.
    """
    path = log_path or _default_log_path()
    now = _now()
    enriched = {"ts": _iso(now), **record}
    if session_id and not enriched.get("session_id"):
        enriched["session_id"] = session_id
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


# ---------------------------------------------------------------------------
# Changes report — the audit log's read side.
#
# "What did you change?" must be answered from the record, not from model
# memory. build_changes_report() aggregates the audit records of one
# session (or a time window) into a small dict; format_changes_report()
# renders that dict as compact markdown for chat display. Both are
# defensive end to end: torn/malformed lines are skipped and neither
# function ever raises.
# ---------------------------------------------------------------------------

_WRITE_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "notebook_edit",
    "apply_patch", "publish_report",
    # Document writes. The fillers record the blank form/template as
    # their ``path``; the file they produced is named in the result, so
    # the report shows what was worked on either way.
    "edit_sheet", "fill_pdf_form", "fill_docx_template", "create_docx",
    "fill_series", "merge_pdfs", "split_pdf", "create_pdf",
})
_COMMAND_TOOLS = frozenset({"bash", "bash_background"})
_PERSIST_TOOLS = frozenset({"remember_permission", "remember_permission_bundle"})
_DENIED_DECISIONS = frozenset({"denied", "error"})


def _read_all_records(log_path: Optional[Path] = None) -> list[dict]:
    """Every parseable record of the active log, oldest first.

    Corrupt / torn lines (kill -9 mid-write) are skipped silently. The
    weekly rotation keeps the active file small, so a full parse is
    cheap; rotated history files are intentionally not consulted.
    """
    path = log_path or _default_log_path()
    if not path.exists():
        return []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except Exception:
        return []
    out: list[dict] = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except Exception:
            continue
        if isinstance(rec, dict):
            out.append(rec)
    return out


def _under_workspace(target: str, workspace: str) -> bool:
    """True unless ``target`` is an absolute path outside ``workspace``.

    Relative paths and empty targets pass — audit records store
    workspace-relative paths for file tools, so only absolute strays
    (e.g. bash cwd escapes into a granted extra dir) are filtered.
    """
    if not target or not workspace:
        return True
    try:
        t = Path(target)
        if not t.is_absolute():
            return True
        t.relative_to(Path(workspace).resolve())
        return True
    except ValueError:
        return False
    except Exception:
        return True


def build_changes_report(
    session_id: Optional[str] = None,
    *,
    since_ts: Optional[str] = None,
    last_n: int = 200,
    workspace: Optional[str] = None,
    log_path: Optional[Path] = None,
) -> dict:
    """Aggregate audit records into a "what changed" report dict.

    Filtering: when ``session_id`` is given, only that session's records
    count (plus ``since_ts`` if also given); otherwise the time window
    ``since_ts`` applies (ISO strings compare lexicographically). At most
    the newest ``last_n`` matching records are aggregated. ``workspace``
    (optional) drops records whose absolute path/cwd lies outside it.

    Returns::

        {"files_written":        [{"path", "tool", "count"}],
         "commands":             [{"command", "cwd"}],
         "denied":               [{"tool", "target", "decision"}],
         "permissions_persisted":[{"tool", "persistence"}],
         "window": {"from_ts", "to_ts", "records"}}

    Never raises — on any failure the empty skeleton comes back.
    """
    empty: dict[str, Any] = {
        "files_written": [], "commands": [], "denied": [],
        "permissions_persisted": [],
        "window": {"from_ts": "", "to_ts": "", "records": 0},
    }
    try:
        records = _read_all_records(log_path)
        picked: list[dict] = []
        for rec in records:
            if session_id is not None and session_id != "":
                if rec.get("session_id", "") != session_id:
                    continue
            if since_ts:
                if str(rec.get("ts", "")) < str(since_ts):
                    continue
            if workspace:
                target = str(rec.get("path", "") or rec.get("cwd", ""))
                if not _under_workspace(target, str(workspace)):
                    continue
            picked.append(rec)
        try:
            cap = max(1, int(last_n))
        except (TypeError, ValueError):
            cap = 200
        picked = picked[-cap:]

        files: dict[str, dict[str, Any]] = {}
        file_tools: dict[str, set] = {}
        commands: list[dict[str, str]] = []
        denied: list[dict[str, str]] = []
        persisted: list[dict[str, Any]] = []
        for rec in picked:
            tool = str(rec.get("tool", ""))
            decision = str(rec.get("decision", ""))
            path = str(rec.get("path", ""))
            command = str(rec.get("command", ""))
            if decision in _DENIED_DECISIONS:
                denied.append({
                    "tool": tool,
                    "target": path or command,
                    "decision": decision,
                })
                continue
            if tool in _WRITE_TOOLS and path:
                entry = files.setdefault(
                    path, {"path": path, "tool": "", "count": 0})
                entry["count"] += 1
                file_tools.setdefault(path, set()).add(tool)
            elif tool in _COMMAND_TOOLS and command:
                commands.append({
                    "command": command,
                    "cwd": str(rec.get("cwd", "")),
                })
            if tool in _PERSIST_TOOLS or rec.get("persistence"):
                persisted.append({
                    "tool": tool,
                    "persistence": rec.get("persistence") or {},
                })
        for path, tools in file_tools.items():
            files[path]["tool"] = "+".join(sorted(tools))
        return {
            "files_written": list(files.values()),
            "commands": commands,
            "denied": denied,
            "permissions_persisted": persisted,
            "window": {
                "from_ts": str(picked[0].get("ts", "")) if picked else "",
                "to_ts": str(picked[-1].get("ts", "")) if picked else "",
                "records": len(picked),
            },
        }
    except Exception:
        return empty


def format_changes_report(report: dict) -> str:
    """Render a build_changes_report() dict as compact markdown.

    Chat-display companion of the report dict; tolerates missing keys
    and never raises.
    """
    try:
        files = report.get("files_written") or []
        commands = report.get("commands") or []
        denied = report.get("denied") or []
        persisted = report.get("permissions_persisted") or []
        window = report.get("window") or {}
        if not (files or commands or denied or persisted):
            return ("No recorded changes — the audit log has no matching "
                    "records for this session/window.")
        lines: list[str] = ["### Changes made"]
        if files:
            lines.append("Files written:")
            for f in files:
                lines.append(
                    f"- `{f.get('path', '?')}` — {f.get('tool', '?')}"
                    f" ×{f.get('count', 0)}")
        if commands:
            lines.append("Commands run:")
            for c in commands:
                cwd = c.get("cwd", "")
                suffix = f" (cwd: {cwd})" if cwd else ""
                lines.append(f"- `{c.get('command', '?')}`{suffix}")
        if denied:
            lines.append("Denied / failed:")
            for d in denied:
                lines.append(
                    f"- {d.get('tool', '?')}: {d.get('target', '?')}"
                    f" [{d.get('decision', 'denied')}]")
        if persisted:
            lines.append("Permissions persisted:")
            for p in persisted:
                lines.append(f"- {p.get('tool', '?')}: "
                             f"{json.dumps(p.get('persistence') or {}, ensure_ascii=False)}")
        n = window.get("records", 0)
        frm, to = window.get("from_ts", ""), window.get("to_ts", "")
        if n:
            span = f" ({frm} → {to})" if frm and to and frm != to else (
                f" ({frm})" if frm else "")
            lines.append(f"_{n} audit record(s){span}_")
        return "\n".join(lines)
    except Exception:
        return "No recorded changes — report unavailable."
