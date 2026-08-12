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

# Resolved on CALL, not bound as a default argument. A module-level default
# is evaluated once at definition time, so pointing `_USER_SETTINGS`
# elsewhere -- a caller with its own config root, or a test -- left every
# function still writing the original file. The user's real settings held
# allow_patterns ["^.*$"] and default_mode "bypassPermissions" for exactly
# this reason: a test persisted an approval rule that could not be
# redirected, and this file decides which commands run without asking.
_VALID_EVENTS = ("PreToolUse", "PostToolUse", "UserPromptSubmit", "Stop")


def _read_settings(path: Path | None = None) -> dict[str, Any]:
    path = _USER_SETTINGS if path is None else path
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {}


def _write_settings(data: dict[str, Any], path: Path | None = None) -> None:
    path = _USER_SETTINGS if path is None else path
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
    layers, in event-order. Used by the ``/hooks`` slash command.

    ``source`` is part of the record. Without it the listing showed
    event, matcher and command, so a hook the user wrote themselves and
    a hook a checked-out repository shipped were indistinguishable in
    the only place they are ever displayed.
    """
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
                "source": cmd.source,
            })
    return out


def hook_warnings(workspace: Path | str | None = None) -> list[str]:
    """Everything the load could not do: definitions withheld for lack of
    trust, entries dropped as malformed. Rendered by ``/hooks``.

    A guard the user cannot see is one they will work around, and a hook
    that silently stopped running looks exactly like a hook that found
    nothing to complain about."""
    return list(_hooks.load_hooks(workspace).warnings)


def trust_state(workspace: Path | str | None = None) -> str:
    """Human-readable trust state of *workspace* for hook definitions."""
    from . import workspace_trust as _trust
    return _trust.describe(workspace)


def trust_this_workspace(workspace: Path | str) -> dict:
    """Record that the USER trusts *workspace*'s hook definitions.

    Called only from the ``/hooks trust`` command -- i.e. only when the
    user typed it. The actor is passed explicitly so that no other
    caller can reach this by accident; ``workspace_trust`` refuses
    anything else.
    """
    from . import workspace_trust as _trust
    return _trust.trust_workspace(
        workspace, [_trust.KIND_HOOKS], actor=_trust.ACTOR_USER)


def untrust_this_workspace(workspace: Path | str) -> bool:
    """Withdraw hook trust for *workspace*. ``/hooks untrust``."""
    from . import workspace_trust as _trust
    return _trust.revoke_workspace(
        workspace, [_trust.KIND_HOOKS], actor=_trust.ACTOR_USER)


def add_hook(
    event: str,
    matcher: str,
    command: str,
    *,
    timeout_s: float | None = None,
    settings_path: Path | None = None,
) -> dict:
    """Append a hook to the user-global settings.json. Returns the
    persisted record so the caller can echo it back to the user."""
    if event not in _VALID_EVENTS:
        raise ValueError(
            f"unknown event {event!r}; valid: {list(_VALID_EVENTS)}"
        )
    if not (command or "").strip():
        raise ValueError("command must be non-empty")
    settings_path = (_USER_SETTINGS if settings_path is None
                     else settings_path)
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
    settings_path: Path | None = None,
) -> dict | None:
    """Remove the ``index``-th hook entry for ``event`` from the user
    settings file. Returns the removed entry or ``None`` if no such
    record existed."""
    if event not in _VALID_EVENTS:
        raise ValueError(
            f"unknown event {event!r}; valid: {list(_VALID_EVENTS)}"
        )
    settings_path = (_USER_SETTINGS if settings_path is None
                     else settings_path)
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
            "source": getattr(r, "source", ""),
            "decision": getattr(r, "decision", "") or "",
            "reason": getattr(r, "reason", "") or "",
            "exit_code": getattr(r, "exit_code", 0),
            "stderr_tail": (getattr(r, "stderr", "") or "")[-200:],
            "duration_s": round(getattr(r, "duration_s", 0.0), 3),
        })
    return out


__all__ = [
    "list_hooks", "add_hook", "remove_hook", "dry_run_hook",
    "hook_warnings", "trust_state", "trust_this_workspace",
    "untrust_this_workspace",
    "_VALID_EVENTS",
]
