"""Helpers used by the ``/hooks`` slash commands.

Wraps ``hooks.load_hooks`` for inspection and provides minimal
add/remove primitives that write to the user-global
``~/.delfin/settings.json`` so the user doesn't have to hand-edit JSON
to register a simple PreToolUse / PostToolUse hook.

All writers are best-effort and idempotent: failure to write surfaces
as a return value, never an exception.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any

from . import hooks as _hooks


_USER_SETTINGS = Path.home() / ".delfin" / "settings.json"
_VALID_EVENTS = ("PreToolUse", "PostToolUse", "UserPromptSubmit", "Stop")


def _read_settings(path: Path = _USER_SETTINGS) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {}


def _write_settings(data: dict[str, Any], path: Path = _USER_SETTINGS) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    tmp.replace(path)
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass


def list_hooks(workspace: Path | str | None = None) -> list[dict]:
    """Return one flat record per registered hook, across all settings
    layers, in event-order. Used by the ``/hooks`` slash command."""
    cfg = _hooks.load_hooks(workspace)
    out: list[dict] = []
    for event in _VALID_EVENTS:
        for i, cmd in enumerate(cfg.for_event(event)):
            out.append({
                "event": event,
                "index": i,
                "matcher": cmd.matcher,
                "command": cmd.command,
                "timeout_s": cmd.timeout_s,
                "type": cmd.type,
            })
    return out


def add_hook(
    event: str,
    matcher: str,
    command: str,
    *,
    timeout_s: float | None = None,
    settings_path: Path = _USER_SETTINGS,
) -> dict:
    """Append a hook to the user-global settings.json. Returns the
    persisted record so the caller can echo it back to the user."""
    if event not in _VALID_EVENTS:
        raise ValueError(
            f"unknown event {event!r}; valid: {list(_VALID_EVENTS)}"
        )
    if not (command or "").strip():
        raise ValueError("command must be non-empty")
    data = _read_settings(settings_path)
    hooks_obj = data.setdefault("hooks", {})
    bucket = hooks_obj.setdefault(event, [])
    entry = {
        "matcher": matcher or "",
        "hooks": [
            {
                "type": "command",
                "command": command,
                **({"timeout": timeout_s * 1000} if timeout_s else {}),
            }
        ],
    }
    bucket.append(entry)
    data.setdefault("_meta", {})["last_modified"] = int(time.time())
    _write_settings(data, settings_path)
    return {"event": event, "matcher": matcher, "command": command,
            "settings_path": str(settings_path)}


def remove_hook(
    event: str,
    index: int,
    *,
    settings_path: Path = _USER_SETTINGS,
) -> dict | None:
    """Remove the ``index``-th hook entry for ``event`` from the user
    settings file. Returns the removed entry or ``None`` if no such
    record existed."""
    if event not in _VALID_EVENTS:
        raise ValueError(
            f"unknown event {event!r}; valid: {list(_VALID_EVENTS)}"
        )
    data = _read_settings(settings_path)
    hooks_obj = data.get("hooks") or {}
    bucket = hooks_obj.get(event) or []
    if not (0 <= index < len(bucket)):
        return None
    removed = bucket.pop(index)
    if not bucket:
        hooks_obj.pop(event, None)
    data["hooks"] = hooks_obj
    _write_settings(data, settings_path)
    return removed


def dry_run_hook(
    event: str,
    *,
    tool_name: str = "",
    tool_input: dict | None = None,
    workspace: Path | str | None = None,
) -> list[dict]:
    """Fire ``event`` through ``hooks.run_hooks`` with synthesised args
    and return a serialisable summary of every command that actually
    ran. Always exits cleanly so the user can probe broken hooks."""
    cfg = _hooks.load_hooks(workspace)
    results = _hooks.run_hooks(
        event,
        cfg,
        tool_name=tool_name,
        arguments=tool_input or {},
        workspace=Path(workspace) if workspace else None,
    )
    out: list[dict] = []
    for r in results or []:
        out.append({
            "matched": bool(getattr(r, "matched", False)),
            "command": getattr(r, "command", ""),
            "decision": getattr(r, "decision", "") or "",
            "reason": getattr(r, "reason", "") or "",
            "exit_code": getattr(r, "exit_code", 0),
            "stderr_tail": (getattr(r, "stderr", "") or "")[-200:],
            "duration_s": round(getattr(r, "duration_s", 0.0), 3),
        })
    return out


__all__ = [
    "list_hooks", "add_hook", "remove_hook", "dry_run_hook",
    "_VALID_EVENTS",
]
