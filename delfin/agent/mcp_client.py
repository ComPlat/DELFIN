"""Minimal stdio-based MCP (Model Context Protocol) client.

Loads server definitions from a JSON config, spawns each server as
a long-lived subprocess speaking line-delimited JSON-RPC over stdin
/ stdout, discovers their tools via ``tools/list``, and surfaces
them to the agent with ``mcp__<server>__<tool>`` namespacing.

Supported only the subset DELFIN needs:

  - JSON-RPC 2.0 framing, line-delimited JSON (Content-Length is
    optional but most stdio servers ship LD-JSON).
  - ``initialize`` handshake (protocol version + minimal client info).
  - ``tools/list`` for discovery.
  - ``tools/call`` for invocation. Response ``result.content`` is
    flattened to a single string.

NOT supported (yet): HTTP/SSE transport, prompts, resources,
sampling, server-initiated requests, base64 image content.

Configuration shape (``~/.delfin/mcp_servers.json`` or per-project
``<workspace>/.delfin/mcp_servers.json``)::

    {
      "servers": {
        "fs": {
          "command": "npx",
          "args": ["-y", "@modelcontextprotocol/server-filesystem", "/tmp"],
          "env": {"FOO": "bar"},
          "enabled": true
        }
      }
    }

The registry is a singleton — first call to ``get_registry()``
loads all configured servers; further calls reuse the running
processes. Servers fail closed: a crash during discovery just
leaves that server's tool set empty.
"""

from __future__ import annotations

import itertools
import json
import os
import subprocess
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional


_DEFAULT_PROTOCOL_VERSION = "2025-06-18"
_RPC_TIMEOUT_S = 15.0
_NAMESPACE_PREFIX = "mcp__"


def _user_config_path() -> Path:
    return Path.home() / ".delfin" / "mcp_servers.json"


def _project_config_path(workspace: Path) -> Path:
    return Path(workspace) / ".delfin" / "mcp_servers.json"


def _load_configs(workspace: Path | None) -> dict[str, dict]:
    """Merge user-global and project-scoped MCP configs."""
    out: dict[str, dict] = {}
    for path in [_user_config_path()] + (
        [_project_config_path(workspace)] if workspace else []
    ):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        servers = data.get("servers") if isinstance(data, dict) else None
        if not isinstance(servers, dict):
            continue
        for name, cfg in servers.items():
            if isinstance(cfg, dict) and cfg.get("enabled", True):
                out[name] = cfg
    return out


@dataclass
class MCPTool:
    server: str
    name: str
    description: str
    schema: dict

    @property
    def namespaced_name(self) -> str:
        return f"{_NAMESPACE_PREFIX}{self.server}__{self.name}"


@dataclass
class MCPServer:
    name: str
    command: str
    args: list[str] = field(default_factory=list)
    env: dict[str, str] = field(default_factory=dict)
    proc: Optional[subprocess.Popen] = None
    _id_counter: itertools.count = field(default_factory=lambda: itertools.count(1))
    _lock: threading.Lock = field(default_factory=threading.Lock)
    _initialised: bool = False
    tools: list[MCPTool] = field(default_factory=list)
    last_error: str = ""

    def start(self) -> None:
        if self.proc is not None and self.proc.poll() is None:
            return
        env = dict(os.environ)
        env.update(self.env or {})
        try:
            self.proc = subprocess.Popen(
                [self.command, *self.args],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                env=env, text=True, bufsize=1,
            )
        except FileNotFoundError as exc:
            self.last_error = f"command not found: {self.command} ({exc})"
            self.proc = None
        except Exception as exc:    # pragma: no cover
            self.last_error = f"start failed: {exc}"
            self.proc = None

    def stop(self) -> None:
        if self.proc is None:
            return
        try:
            self.proc.terminate()
            self.proc.wait(timeout=2)
        except subprocess.TimeoutExpired:
            self.proc.kill()
        except Exception:    # pragma: no cover
            pass
        self.proc = None
        self._initialised = False

    def _send(self, method: str, params: dict | None = None) -> dict:
        if self.proc is None or self.proc.poll() is not None:
            return {"error": {"message": "server not running"}}
        rid = next(self._id_counter)
        req = {
            "jsonrpc": "2.0",
            "id": rid,
            "method": method,
            "params": params or {},
        }
        line = json.dumps(req) + "\n"
        with self._lock:
            try:
                assert self.proc.stdin is not None
                self.proc.stdin.write(line)
                self.proc.stdin.flush()
            except (BrokenPipeError, OSError) as exc:
                return {"error": {"message": f"write failed: {exc}"}}
            t_end = time.monotonic() + _RPC_TIMEOUT_S
            assert self.proc.stdout is not None
            while time.monotonic() < t_end:
                raw = self.proc.stdout.readline()
                if not raw:
                    return {"error": {"message": "EOF from server"}}
                raw = raw.strip()
                if not raw:
                    continue
                try:
                    msg = json.loads(raw)
                except json.JSONDecodeError:
                    continue
                if msg.get("id") == rid:
                    return msg
                # ignore notifications / unrelated messages
        return {"error": {"message": "rpc timeout"}}

    def initialize(self) -> bool:
        if self._initialised:
            return True
        if self.proc is None:
            self.start()
        if self.proc is None:
            return False
        resp = self._send("initialize", {
            "protocolVersion": _DEFAULT_PROTOCOL_VERSION,
            "capabilities": {},
            "clientInfo": {"name": "delfin-agent", "version": "0.1"},
        })
        if "error" in resp:
            self.last_error = json.dumps(resp["error"])[:240]
            return False
        # Required follow-up notification
        try:
            assert self.proc.stdin is not None
            self.proc.stdin.write(json.dumps({
                "jsonrpc": "2.0",
                "method": "notifications/initialized",
            }) + "\n")
            self.proc.stdin.flush()
        except Exception as exc:    # pragma: no cover
            self.last_error = f"initialized notify failed: {exc}"
        self._initialised = True
        return True

    def list_tools(self) -> list[MCPTool]:
        if not self.initialize():
            return []
        resp = self._send("tools/list")
        if "error" in resp:
            self.last_error = json.dumps(resp["error"])[:240]
            return []
        result = resp.get("result", {})
        items = result.get("tools", []) if isinstance(result, dict) else []
        out: list[MCPTool] = []
        for t in items:
            if not isinstance(t, dict):
                continue
            name = str(t.get("name", "")).strip()
            if not name:
                continue
            out.append(MCPTool(
                server=self.name,
                name=name,
                description=str(t.get("description", "")),
                schema=t.get("inputSchema") or {"type": "object"},
            ))
        self.tools = out
        return out

    def call_tool(self, tool_name: str, arguments: dict) -> str:
        if not self.initialize():
            return json.dumps({"error": self.last_error or "init failed"})
        resp = self._send("tools/call", {
            "name": tool_name,
            "arguments": arguments or {},
        })
        if "error" in resp:
            return json.dumps({"error": resp["error"]})
        result = resp.get("result", {})
        content = result.get("content", []) if isinstance(result, dict) else []
        # Flatten ``content: [{type:'text', text:'...'}]`` to a string.
        texts: list[str] = []
        for c in content:
            if isinstance(c, dict):
                if c.get("type") == "text":
                    texts.append(str(c.get("text", "")))
                else:
                    texts.append(json.dumps(c)[:1000])
        if texts:
            return "\n".join(texts)
        return json.dumps(result)


@dataclass
class MCPRegistry:
    servers: dict[str, MCPServer] = field(default_factory=dict)
    workspace: Optional[Path] = None
    loaded: bool = False

    def load(self, workspace: Path | None = None) -> None:
        configs = _load_configs(workspace)
        for name, cfg in configs.items():
            if name in self.servers:
                continue
            srv = MCPServer(
                name=name,
                command=str(cfg.get("command", "")),
                args=list(cfg.get("args", []) or []),
                env={k: str(v) for k, v in (cfg.get("env") or {}).items()},
            )
            self.servers[name] = srv
        self.workspace = Path(workspace) if workspace else self.workspace
        self.loaded = True

    def discover_all(self) -> list[MCPTool]:
        tools: list[MCPTool] = []
        for srv in self.servers.values():
            try:
                tools.extend(srv.list_tools())
            except Exception:   # pragma: no cover
                pass
        return tools

    def call(self, namespaced: str, arguments: dict) -> str:
        if not namespaced.startswith(_NAMESPACE_PREFIX):
            return json.dumps({"error": f"not an MCP tool: {namespaced!r}"})
        rest = namespaced[len(_NAMESPACE_PREFIX):]
        if "__" not in rest:
            return json.dumps({"error": f"malformed MCP name: {namespaced!r}"})
        server_name, _, tool_name = rest.partition("__")
        srv = self.servers.get(server_name)
        if srv is None:
            return json.dumps({"error": f"unknown MCP server: {server_name!r}"})
        return srv.call_tool(tool_name, arguments)

    def shutdown(self) -> None:
        for srv in self.servers.values():
            srv.stop()


_REGISTRY_LOCK = threading.Lock()
_REGISTRY: MCPRegistry | None = None


def get_registry(workspace: Path | None = None) -> MCPRegistry:
    """Return the process-wide MCP registry, lazily creating it."""
    global _REGISTRY
    with _REGISTRY_LOCK:
        if _REGISTRY is None:
            _REGISTRY = MCPRegistry()
        if not _REGISTRY.loaded:
            _REGISTRY.load(workspace)
    return _REGISTRY


def reset_registry() -> None:
    """Stop and clear the global registry. Used by tests."""
    global _REGISTRY
    with _REGISTRY_LOCK:
        if _REGISTRY is not None:
            _REGISTRY.shutdown()
        _REGISTRY = None


__all__ = [
    "MCPTool",
    "MCPServer",
    "MCPRegistry",
    "get_registry",
    "reset_registry",
]
