"""Persistent session storage for the DELFIN Agent.

Saves and loads agent sessions so conversations survive dashboard restarts.
Each session is a JSON file in ``~/.delfin/agent_sessions/`` (per-user,
per-machine, 0600 permissions).  Sessions MUST NOT be committed to the repo.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any


def _chmod_user_only(path: Path) -> None:
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass

# Legacy mode name migration
_LEGACY_MODE_MAP = {
    "default": "quick",
    "high_risk": "reviewed",
    "runtime_cluster": "cluster",
    "release": "full",
}


def _migrate_mode(mode: str) -> str:
    return _LEGACY_MODE_MAP.get(mode, mode)


_SESSIONS_DIR = Path.home() / ".delfin" / "agent_sessions"


def _ensure_dir() -> Path:
    _SESSIONS_DIR.mkdir(parents=True, exist_ok=True)
    return _SESSIONS_DIR


def _transcript_archive_dir() -> Path:
    """Directory where full pre-compaction transcripts are archived so
    long sessions never truly lose information — the user can scroll
    back even after the in-memory ``messages`` list was compacted."""
    p = Path.home() / ".delfin" / "transcript_archive"
    p.mkdir(parents=True, exist_ok=True)
    return p


def archive_pre_compaction_transcript(
    session_id: str,
    messages: list[dict[str, Any]],
    *,
    info: dict[str, Any] | None = None,
) -> Path:
    """Persist the full pre-compaction conversation to JSONL so the user
    can browse it later (e.g. via /session archive ls). Returns the
    archive file path.

    The file is append-only: every compaction in a session adds a new
    record block at the end with ``info`` (msg count, tokens saved,
    timestamp). Cheap to read tail-first; safe to leave for months.
    """
    if not session_id:
        return Path(os.devnull)
    d = _transcript_archive_dir()
    p = d / f"{session_id}.jsonl"
    rec = {
        "compacted_at": time.time(),
        "n_messages": len(messages),
        "info": info or {},
        "messages": [
            {
                "role": m.get("role", ""),
                "content": (m.get("content", "") if isinstance(m.get("content"), str)
                            else json.dumps(m.get("content"), ensure_ascii=False)),
            }
            for m in messages
        ],
    }
    try:
        with p.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        _chmod_user_only(p)
    except OSError:
        return Path(os.devnull)
    return p


def list_transcript_archives() -> list[dict[str, Any]]:
    """Return one record per archived session: session_id, size, mtime,
    n_compaction_records (one line per compaction event)."""
    d = _transcript_archive_dir()
    out: list[dict[str, Any]] = []
    for p in d.glob("*.jsonl"):
        try:
            stat = p.stat()
            n_lines = sum(1 for _ in p.open("r", encoding="utf-8"))
        except OSError:
            continue
        out.append({
            "session_id": p.stem,
            "path": str(p),
            "bytes": stat.st_size,
            "mtime": stat.st_mtime,
            "compactions": n_lines,
        })
    out.sort(key=lambda r: r["mtime"], reverse=True)
    return out


def load_transcript_archive(session_id: str) -> list[dict[str, Any]]:
    """Return every compaction record in order (oldest → newest) for the
    given session. Each record contains the messages snapshot at that
    compaction. Empty list when no archive exists."""
    if not session_id:
        return []
    p = _transcript_archive_dir() / f"{session_id}.jsonl"
    if not p.exists():
        return []
    out: list[dict[str, Any]] = []
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
    return out


def save_session(
    session_id: str,
    *,
    mode: str = "quick",
    role_index: int = 0,
    route: list[str] | None = None,
    role_outputs: dict[str, str] | None = None,
    chat_messages: list[dict[str, Any]] | None = None,
    cycle_history: list[dict[str, Any]] | None = None,
    engine_messages: list[dict[str, Any]] | None = None,
    token_usage: dict[str, int] | None = None,
    cost_usd: float = 0.0,
    title: str = "",
    # Long-session extras (all optional, ignored if backend doesn't pass)
    perm_profile: str = "",
    provider: str = "",
    model: str = "",
    effort: str = "",
    active_gate: dict[str, Any] | None = None,
    last_compaction_info: dict[str, Any] | None = None,
    subagent_calls: list[dict[str, Any]] | None = None,
    pending_plan_body: str = "",
    todo_payload: list[dict[str, Any]] | None = None,
    transcript_archive_path: str = "",
) -> Path:
    """Save a session to disk.

    Parameters
    ----------
    session_id : str
        The CLI session ID (used as filename).
    mode : str
        Active DELFIN_AGENT_LITE mode.
    role_index : int
        Current role index in the route.
    route : list[str]
        The mode's role route.
    role_outputs : dict
        Collected outputs from completed roles.
    chat_messages : list
        The UI chat messages (user/assistant/system with role_label).
    cycle_history : list
        Recent gate / handoff / retry events for the Cycle Inspector.
    engine_messages : list
        The engine's conversation messages (role + content).
    token_usage : dict
        ``{"input": N, "output": N}``.
    cost_usd : float
        Accumulated cost.
    title : str
        Short description (auto-generated from first user message if empty).

    Returns
    -------
    Path
        Path to the saved session file.
    """
    if not session_id:
        return Path(os.devnull)

    d = _ensure_dir()

    # Auto-title from first user message
    if not title and chat_messages:
        for msg in chat_messages:
            if msg.get("role") == "user":
                text = msg.get("content", "")
                title = text[:80].replace("\n", " ").strip()
                if len(text) > 80:
                    title += "..."
                break

    data = {
        "session_id": session_id,
        "mode": mode,
        "role_index": role_index,
        "route": route or [],
        "role_outputs": role_outputs or {},
        "chat_messages": chat_messages or [],
        "cycle_history": cycle_history or [],
        "engine_messages": engine_messages or [],
        "token_usage": token_usage or {"input": 0, "output": 0},
        "cost_usd": cost_usd,
        "title": title,
        "updated_at": time.time(),
        # Long-session state (only persisted when callers pass it)
        "perm_profile": perm_profile or "",
        "provider": provider or "",
        "model": model or "",
        "effort": effort or "",
        "active_gate": active_gate or None,
        "last_compaction_info": last_compaction_info or None,
        "subagent_calls": subagent_calls or [],
        "pending_plan_body": pending_plan_body or "",
        "todo_payload": todo_payload or [],
        "transcript_archive_path": transcript_archive_path or "",
    }

    # Set created_at only on first save
    filepath = d / f"{session_id}.json"
    if filepath.exists():
        try:
            existing = json.loads(filepath.read_text())
            data["created_at"] = existing.get("created_at", data["updated_at"])
        except (json.JSONDecodeError, OSError):
            data["created_at"] = data["updated_at"]
    else:
        data["created_at"] = data["updated_at"]

    filepath.write_text(json.dumps(data, ensure_ascii=False, indent=1))
    _chmod_user_only(filepath)
    return filepath


def load_session(session_id: str) -> dict[str, Any] | None:
    """Load a session from disk.

    Returns None if the session file doesn't exist or is corrupt.
    """
    filepath = _SESSIONS_DIR / f"{session_id}.json"
    if not filepath.exists():
        return None
    try:
        return json.loads(filepath.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def list_sessions(limit: int = 50) -> list[dict[str, Any]]:
    """List saved sessions, newest first.

    Returns a list of lightweight session summaries (no full chat history).
    """
    d = _SESSIONS_DIR
    if not d.exists():
        return []

    sessions = []
    for f in d.glob("*.json"):
        try:
            data = json.loads(f.read_text())
            sessions.append({
                "session_id": data.get("session_id", f.stem),
                "title": data.get("title", "Untitled"),
                "mode": _migrate_mode(data.get("mode", "quick")),
                "role_index": data.get("role_index", 0),
                "route": data.get("route", []),
                "cost_usd": data.get("cost_usd", 0.0),
                "token_usage": data.get("token_usage", {}),
                "message_count": len(data.get("chat_messages", [])),
                "created_at": data.get("created_at", 0),
                "updated_at": data.get("updated_at", 0),
            })
        except (json.JSONDecodeError, OSError):
            continue

    sessions.sort(key=lambda s: s.get("updated_at", 0), reverse=True)
    return sessions[:limit]


def delete_session(session_id: str) -> bool:
    """Delete a session file. Returns True if deleted."""
    filepath = _SESSIONS_DIR / f"{session_id}.json"
    if filepath.exists():
        filepath.unlink()
        return True
    return False


def fork_session(
    source_id: str,
    new_id: str | None = None,
    *,
    title_suffix: str = " (fork)",
    sessions_dir: Path | None = None,
) -> str | None:
    """Duplicate an existing session under a new ID.

    The fork inherits all chat / engine / cycle history but gets a fresh
    ``session_id``, fresh timestamps, and an appended title suffix.  The
    ``cost_usd`` and ``token_usage`` counters stay attached so the fork's
    cumulative state is honest.

    Parameters
    ----------
    source_id : str
        The session ID to fork.  Must exist on disk.
    new_id : str, optional
        New session ID.  If omitted, a timestamp-based ID is generated.
    title_suffix : str
        Appended to the title to disambiguate from the source.
    sessions_dir : Path, optional
        Override the default ``~/.delfin/agent_sessions/`` path (for tests).

    Returns
    -------
    str | None
        The new session ID, or ``None`` if the source could not be loaded.
    """
    base_dir = sessions_dir or _SESSIONS_DIR
    src_path = base_dir / f"{source_id}.json"
    if not src_path.exists():
        return None
    try:
        data = json.loads(src_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

    if not new_id:
        new_id = f"fork-{int(time.time())}-{source_id[-6:] if source_id else 'src'}"

    # Mutate to make this a fresh standalone session
    data["session_id"] = new_id
    src_title = str(data.get("title") or "Untitled")
    if title_suffix and title_suffix not in src_title:
        data["title"] = f"{src_title}{title_suffix}"
    now = time.time()
    data["created_at"] = now
    data["updated_at"] = now
    data["forked_from"] = source_id

    out_path = base_dir / f"{new_id}.json"
    base_dir.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, ensure_ascii=False, indent=1))
    _chmod_user_only(out_path)
    return new_id


def session_lineage(session_id: str) -> list[dict[str, Any]]:
    """Walk the ``forked_from`` chain back to the root and return one
    summary record per ancestor, newest (the given session) first.

    Used by ``/session tree`` to show the branch lineage so the user can
    see where a fork came from. Stops at the first session without a
    ``forked_from`` field or when a parent is missing on disk.
    """
    out: list[dict[str, Any]] = []
    seen: set[str] = set()
    current = session_id
    while current and current not in seen:
        seen.add(current)
        data = load_session(current)
        if not data:
            break
        out.append({
            "session_id": current,
            "title": data.get("title", "Untitled"),
            "forked_from": data.get("forked_from", ""),
            "updated_at": data.get("updated_at", 0),
            "message_count": len(data.get("chat_messages") or []),
            "cost_usd": data.get("cost_usd", 0.0),
        })
        parent = data.get("forked_from") or ""
        current = parent
    return out


def session_children(session_id: str) -> list[dict[str, Any]]:
    """Return immediate child sessions (i.e. sessions whose
    ``forked_from`` points at ``session_id``). Newest first."""
    out: list[dict[str, Any]] = []
    for summary in list_sessions(limit=200):
        sid = summary.get("session_id", "")
        if not sid:
            continue
        data = load_session(sid)
        if not data:
            continue
        if data.get("forked_from") == session_id:
            out.append({
                "session_id": sid,
                "title": data.get("title", "Untitled"),
                "updated_at": data.get("updated_at", 0),
                "message_count": len(data.get("chat_messages") or []),
            })
    out.sort(key=lambda r: r.get("updated_at", 0), reverse=True)
    return out


def resume_latest(*, max_age_s: int | None = None) -> dict[str, Any] | None:
    """Return the most recently updated session, or None.

    Parameters
    ----------
    max_age_s : int, optional
        If set, ignore sessions whose ``updated_at`` is older than
        this many seconds. Useful for "resume only if I was just
        working on it" UX.
    """
    summaries = list_sessions(limit=1)
    if not summaries:
        return None
    head = summaries[0]
    if max_age_s is not None and head.get("updated_at"):
        if time.time() - float(head["updated_at"]) > max_age_s:
            return None
    return load_session(str(head["session_id"]))


def cleanup_old_sessions(*, max_age_days: int = 30) -> int:
    """Prune session files older than ``max_age_days``. Returns count."""
    if not _SESSIONS_DIR.exists():
        return 0
    cutoff = time.time() - max_age_days * 86_400
    removed = 0
    for f in _SESSIONS_DIR.glob("*.json"):
        try:
            mtime = f.stat().st_mtime
        except OSError:
            continue
        if mtime < cutoff:
            try:
                f.unlink()
                removed += 1
            except OSError:
                pass
    return removed


def find_sessions(query: str, *, limit: int = 20) -> list[dict[str, Any]]:
    """Substring search over session titles + first user message."""
    needle = query.strip().lower()
    if not needle:
        return list_sessions(limit=limit)
    out: list[dict[str, Any]] = []
    for s in list_sessions(limit=200):
        if needle in str(s.get("title", "")).lower():
            out.append(s)
            continue
        # Cheap dive into the file for first-message match
        full = load_session(str(s["session_id"]))
        if not full:
            continue
        msgs = full.get("chat_messages") or []
        first_user = ""
        for m in msgs:
            if m.get("role") == "user":
                first_user = str(m.get("content", ""))
                break
        if needle in first_user.lower():
            out.append(s)
        if len(out) >= limit:
            break
    return out
