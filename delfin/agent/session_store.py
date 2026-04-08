"""Persistent session storage for the DELFIN Agent.

Saves and loads agent sessions so conversations survive dashboard restarts.
Each session is a JSON file in ``~/.delfin/agent_sessions/``.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any

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
