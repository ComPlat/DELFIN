"""Tiny CRUD helpers for the user-global MCP server config so the
``/mcp`` slash command doesn't have to hand-edit JSON.

Writes to ``~/.delfin/mcp_servers.json`` (the user-global layer). Project-
scoped config remains hand-edited under ``<workspace>/.delfin/mcp_servers.json``
since project-level changes usually want code-review.
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
    "add_mcp_server",
    "remove_mcp_server",
    "toggle_mcp_server",
]
