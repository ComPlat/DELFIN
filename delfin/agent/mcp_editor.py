"""Tiny CRUD helpers for the user-global MCP server config so the
``/mcp`` slash command doesn't have to hand-edit JSON.

Writes to ``~/.delfin/mcp_servers.json`` (the user-global layer). Project-
scoped config remains hand-edited under ``<workspace>/.delfin/mcp_servers.json``
since project-level changes usually want code-review — and is honoured
only for a directory the user has trusted, which is what ``trust_state``
and ``trust_this_workspace`` here are for.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any


_USER_CONFIG = Path.home() / ".delfin" / "mcp_servers.json"


def _read(path: Path = _USER_CONFIG) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {"servers": {}}


def _write(data: dict[str, Any], path: Path = _USER_CONFIG) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    tmp.replace(path)
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass


def list_effective_mcp_servers(
    workspace: Path | str | None = None,
) -> list[dict[str, Any]]:
    """Every server that would actually load, each named by its source.

    ``list_mcp_servers`` reads one file. This reads what the client
    reads: the built-ins, the user's config, and a trusted workspace's
    config, with the file each entry came from. The ``/mcp`` listing used
    the single-file view, so a server a workspace contributed was spawned
    and never appeared — invisible in the one case where the source is
    the whole question.
    """
    from .mcp_client import effective_servers
    return effective_servers(Path(workspace) if workspace else None)


def mcp_trust_notice(workspace: Path | str | None = None) -> str:
    """What a lack of trust withheld for *workspace*, or ""."""
    from .mcp_client import trust_notice
    return trust_notice(Path(workspace) if workspace else None)


def trust_state(workspace: Path | str | None = None) -> str:
    """Human-readable trust state of *workspace* for server definitions."""
    from . import workspace_trust as _trust
    return _trust.describe(workspace)


def trust_this_workspace(workspace: Path | str) -> dict:
    """Record that the USER trusts *workspace*'s MCP server definitions.

    Called only from the ``/mcp trust`` command — i.e. only when the user
    typed it. The actor is passed explicitly so no other caller reaches
    this by accident; ``workspace_trust`` refuses anything else.
    """
    from . import workspace_trust as _trust
    return _trust.trust_workspace(
        workspace, [_trust.KIND_MCP_SERVERS], actor=_trust.ACTOR_USER)


def untrust_this_workspace(workspace: Path | str) -> bool:
    """Withdraw MCP-server trust for *workspace*. ``/mcp untrust``."""
    from . import workspace_trust as _trust
    return _trust.revoke_workspace(
        workspace, [_trust.KIND_MCP_SERVERS], actor=_trust.ACTOR_USER)


def list_mcp_servers(path: Path = _USER_CONFIG) -> list[dict[str, Any]]:
    """Return one record per configured server (name, command, args,
    enabled flag). Reads from the user-global config only."""
    data = _read(path)
    servers = data.get("servers") or {}
    out: list[dict[str, Any]] = []
    for name, cfg in servers.items():
        if not isinstance(cfg, dict):
            continue
        out.append({
            "name": name,
            "command": cfg.get("command", ""),
            "args": cfg.get("args") or [],
            "env": cfg.get("env") or {},
            "enabled": bool(cfg.get("enabled", True)),
        })
    out.sort(key=lambda r: r["name"])
    return out


def add_mcp_server(
    name: str,
    command: str,
    args: list[str] | None = None,
    env: dict[str, str] | None = None,
    *,
    enabled: bool = True,
    path: Path = _USER_CONFIG,
) -> dict[str, Any]:
    """Append or overwrite a server entry. Returns the persisted record."""
    if not name or not isinstance(name, str):
        raise ValueError("server name must be a non-empty string")
    if not command:
        raise ValueError("command must be a non-empty string")
    data = _read(path)
    servers = data.setdefault("servers", {})
    servers[name] = {
        "command": command,
        "args": list(args or []),
        "env": dict(env or {}),
        "enabled": bool(enabled),
    }
    data.setdefault("_meta", {})["last_modified"] = int(time.time())
    _write(data, path)
    return {"name": name, **servers[name]}


def remove_mcp_server(
    name: str, *, path: Path = _USER_CONFIG,
) -> dict[str, Any] | None:
    """Remove a server by name. Returns the deleted entry or None."""
    data = _read(path)
    servers = data.get("servers") or {}
    removed = servers.pop(name, None)
    if removed is None:
        return None
    data["servers"] = servers
    _write(data, path)
    return {"name": name, **removed}


def toggle_mcp_server(
    name: str, *, enabled: bool, path: Path = _USER_CONFIG,
) -> dict[str, Any] | None:
    """Flip ``enabled`` on an existing server. Returns the updated record
    or None if the server doesn't exist."""
    data = _read(path)
    servers = data.get("servers") or {}
    cfg = servers.get(name)
    if not isinstance(cfg, dict):
        return None
    cfg["enabled"] = bool(enabled)
    data["servers"] = servers
    _write(data, path)
    return {"name": name, **cfg}


__all__ = [
    "list_mcp_servers",
    "list_effective_mcp_servers",
    "mcp_trust_notice",
    "trust_state",
    "trust_this_workspace",
    "untrust_this_workspace",
    "add_mcp_server",
    "remove_mcp_server",
    "toggle_mcp_server",
]
