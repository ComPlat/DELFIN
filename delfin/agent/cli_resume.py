"""CLI resume: come back to a kept agent session from the terminal.

Reads ONLY the stores that already exist — the per-session JSON files in
``~/.delfin/agent_sessions/`` written by ``delfin.agent.session_store``
(and the chat history inside them). No new store, no new file format.

Three functions, no engine wiring (that lives in ``delfin/agent/cli.py``,
which this module deliberately does not touch):

- :func:`list_sessions`  -> one summary row per kept session, newest first
- :func:`resume_target`  -> resolve "" / id-prefix / index to ONE session
- :func:`render_sessions`-> plain-text table of :func:`list_sessions` rows

Resuming is refused — with a clear message, via :class:`ResumeError` —
when the session's workspace directory no longer exists on disk: a
resumed engine would otherwise silently run in whatever directory the
user happens to be standing in.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from .session_store import _SESSIONS_DIR  # location only; never mutated here


class ResumeError(Exception):
    """A resume request that must not proceed (unknown session, ambiguous
    id prefix, or a session whose workspace directory is gone)."""


def _iter_session_files(sessions_dir: Path) -> list[Path]:
    if not sessions_dir.is_dir():
        return []
    files = []
    for f in sessions_dir.glob("*.json"):
        if f.name.endswith(".turn.json"):
            continue  # mid-turn crash checkpoint, not a session
        files.append(f)
    return files


def _load_row(path: Path) -> dict[str, Any] | None:
    """One summary row from a session file, or None if unreadable.

    Mirrors session_store.list_sessions semantics (skip corrupt files,
    newest by ``updated_at`` first) but reads the directory we were given
    so tests can fabricate a store in ``tmp_path``.
    """
    try:
        data = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    chat = data.get("chat_messages") or []
    if not isinstance(chat, list):
        chat = []
    # Last user line = the last task the user actually typed. Tool results
    # and assistant prose carry no task; skip them.
    last_user = ""
    for msg in chat:
        if isinstance(msg, dict) and msg.get("role") == "user":
            content = msg.get("content")
            if isinstance(content, list):  # multimodal shape
                content = " ".join(
                    str(p.get("text", "")) for p in content
                    if isinstance(p, dict)
                )
            if isinstance(content, str) and content.strip():
                last_user = " ".join(content.split())
    return {
        "id": str(data.get("session_id") or path.stem),
        "started": float(data.get("created_at", 0) or 0),
        "updated_at": float(data.get("updated_at", 0) or 0),
        "model": str(data.get("model") or ""),
        "workspace": str(data.get("workspace") or ""),
        "turns": len(chat),
        "last_task_line": last_user,
    }


def list_sessions(*, limit: int = 20, sessions_dir: str | Path | None = None) -> list[dict]:
    """Summaries of kept sessions, newest first (by ``updated_at``).

    Rows: id, started, updated_at, model, workspace, turns,
    last_task_line. Unreadable/corrupt files are skipped, never raised.
    """
    d = Path(sessions_dir) if sessions_dir else _SESSIONS_DIR
    rows = []
    for f in _iter_session_files(d):
        row = _load_row(f)
        if row is not None:
            rows.append(row)
    rows.sort(key=lambda r: r.get("updated_at", 0), reverse=True)
    return rows[:limit]


def _workspace_exists(workspace: str) -> bool:
    """A session with no recorded workspace is resumable (it predates the
    field); one WITH a workspace requires that directory to still exist."""
    if not workspace:
        return True
    try:
        return Path(workspace).expanduser().is_dir()
    except (OSError, ValueError):
        return False


def resume_target(
    selector: str = "", *, sessions_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Resolve ``selector`` to exactly one resumable session.

    selector: "" -> the most recently updated session; a session-id
    prefix (must match exactly one session); or a 1-based index into the
    same newest-first order :func:`list_sessions` returns (also accepted
    as "#3").

    Raises :class:`ResumeError` — never returns None — when nothing
    matches, a prefix matches several sessions, the store is empty, or
    the session's workspace directory is gone. The error text is meant
    for the terminal and names the selector and what to do about it.
    """
    rows = list_sessions(limit=1000, sessions_dir=sessions_dir)
    if not rows:
        raise ResumeError("no kept sessions found — nothing to resume")

    sel = selector.strip().lstrip("#")

    row = None
    if not sel:
        row = rows[0]
    elif sel.isdigit():
        idx = int(sel)
        if not 1 <= idx <= len(rows):
            raise ResumeError(
                f"session index {idx} out of range "
                f"(1..{len(rows)}; run --sessions to list them)"
            )
        row = rows[idx - 1]
    else:
        matches = [r for r in rows if r["id"].startswith(sel)]
        if not matches:
            raise ResumeError(
                f"no kept session starts with {sel!r} "
                f"(run --sessions to list them)"
            )
        if len(matches) > 1:
            ids = ", ".join(m["id"] for m in matches[:5])
            raise ResumeError(
                f"session id prefix {sel!r} is ambiguous ({len(matches)} "
                f"matches: {ids}...) — give more characters"
            )
        row = matches[0]

    if not _workspace_exists(row["workspace"]):
        raise ResumeError(
            f"session {row['id']} worked in {row['workspace']!r}, "
            "which no longer exists — refusing to resume there. "
            "Recreate the directory or pick another session."
        )
    return row


def _fmt_age(ts: float) -> str:
    if not ts:
        return "?"
    delta = max(0.0, time.time() - ts)
    if delta < 3600:
        return f"{int(delta // 60)}m ago"
    if delta < 86400:
        return f"{int(delta // 3600)}h ago"
    return f"{int(delta // 86400)}d ago"


def _fmt_line(text: str, width: int = 48) -> str:
    line = " ".join(text.split())
    if len(line) > width:
        line = line[: width - 1] + "…"
    return line


def render_sessions(rows: list[dict]) -> str:
    """Render :func:`list_sessions` rows as a terminal table + hint lines.

    Never raises; empty input renders the 'no kept sessions' message.
    """
    if not rows:
        return "no kept sessions found"
    header = f"{'#':>3}  {'session id':<16}  {'age':>7}  {'model':<14}  task"
    lines = [header, "-" * len(header)]
    for i, r in enumerate(rows, 1):
        task = _fmt_line(r.get("last_task_line") or "(no user message kept)")
        model = r.get("model") or "?"
        lines.append(
            f"{i:>3}  {r['id'][:16]:<16}  {_fmt_age(r.get('updated_at', 0)):>7}  "
            f"{model[:14]:<14}  {task}"
        )
    lines.append("")
    lines.append("resume with: delfin-agent --resume <id-prefix | #index | ''>")
    return "\n".join(lines)
