"""MCP (Model Context Protocol) client — stdio + Streamable HTTP.

Loads server definitions from a JSON config and reaches each one over
either transport, discovers their tools via ``tools/list``, and surfaces
them to the agent with ``mcp__<server>__<tool>`` namespacing.

Transports:

  - **stdio** (default): the server is spawned as a long-lived subprocess
    speaking line-delimited JSON-RPC 2.0 over stdin/stdout.
  - **http** (Streamable HTTP, MCP 2025-03-26+): the server is reached by
    HTTP ``POST`` to a single ``url``. The response is either a plain JSON
    body or an ``text/event-stream`` (SSE) carrying the JSON-RPC reply. The
    ``Mcp-Session-Id`` header returned on ``initialize`` is echoed on every
    subsequent request. Custom ``headers`` (e.g. bearer auth) are supported.

Protocol coverage: ``initialize`` handshake, ``tools/list`` + ``tools/call``,
``resources/list`` + ``resources/read``, ``prompts/list`` + ``prompts/get``.
NOT supported (yet): sampling, server-initiated requests, base64 image content.

Configuration shape (``~/.delfin/mcp_servers.json`` or per-project
``<workspace>/.delfin/mcp_servers.json``)::

    {
      "servers": {
        "fs": {
          "command": "npx",
          "args": ["-y", "@modelcontextprotocol/server-filesystem", "/tmp"],
          "env": {"FOO": "bar"},
          "enabled": true
        },
        "remote": {
          "type": "http",
          "url": "https://mcp.example.com/mcp",
          "headers": {"Authorization": "Bearer ..."},
          "enabled": true
        }
      }
    }

A server is treated as HTTP when ``type`` is ``http``/``sse``/
``streamable-http`` or a ``url`` is present; otherwise stdio.

The registry is a singleton — first call to ``get_registry()``
loads all configured servers; further calls reuse the running
processes / sessions. Servers fail closed: a crash during discovery
just leaves that server's tool set empty.
"""

from __future__ import annotations

import itertools
import json
import os
import re
import subprocess
import threading
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional


_DEFAULT_PROTOCOL_VERSION = "2025-06-18"
_RPC_TIMEOUT_S = 15.0
_NAMESPACE_PREFIX = "mcp__"
_HTTP_TYPES = {"http", "sse", "streamable-http", "streamablehttp"}


def _user_config_path() -> Path:
    return Path.home() / ".delfin" / "mcp_servers.json"


def _project_config_path(workspace: Path) -> Path:
    return Path(workspace) / ".delfin" / "mcp_servers.json"


# Built-in MCP servers that ship WITH DELFIN — always available so a mode like
# Pipeline works out-of-the-box without each user hand-registering them. The
# user config can override an entry's command/args or disable it entirely
# (set ``"enabled": false`` for that name).
_BUILTIN_SERVERS: dict[str, dict] = {
    "delfin-tools": {
        "command": "python",
        "args": ["-m", "delfin.tools.mcp_server"],
        "enabled": True,
    },
}


def _load_configs(workspace: Path | None) -> dict[str, dict]:
    """Merge built-in defaults + user-global + project-scoped MCP configs."""
    out: dict[str, dict] = {
        name: dict(cfg) for name, cfg in _BUILTIN_SERVERS.items()
    }
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
            if not isinstance(cfg, dict):
                continue
            if cfg.get("enabled", True):
                out[name] = cfg                 # add or override (incl. a builtin)
            else:
                out.pop(name, None)             # explicit disable (also disables a builtin)
    return out


def _flatten_content(content: Any) -> str:
    """Flatten an MCP content value (str | {type:text,text} | list) to text."""
    if content is None:
        return ""
    if isinstance(content, str):
        return content
    if isinstance(content, dict):
        if content.get("type") == "text" or content.get("text") is not None:
            return str(content.get("text", ""))
        return json.dumps(content)[:1000]
    if isinstance(content, list):
        out: list[str] = []
        for c in content:
            if isinstance(c, dict):
                if c.get("type") == "text" or c.get("text") is not None:
                    out.append(str(c.get("text", "")))
                else:
                    out.append(json.dumps(c)[:1000])
            elif isinstance(c, str):
                out.append(c)
        return "\n".join(out)
    return str(content)


def _extract_jsonrpc_from_sse(raw: str, rid: Any) -> dict:
    """Pull the JSON-RPC message matching ``rid`` out of an SSE body.

    SSE frames are blank-line-separated; ``data:`` lines carry the payload
    (concatenated with newlines per spec). Returns the frame whose ``id``
    matches, else the first JSON object seen, else an error dict.
    """
    best: dict | None = None
    for block in re.split(r"\n\n+", raw):
        datas = [ln[len("data:"):].lstrip()
                 for ln in block.splitlines() if ln.startswith("data:")]
        if not datas:
            continue
        try:
            msg = json.loads("\n".join(datas))
        except json.JSONDecodeError:
            continue
        if isinstance(msg, dict):
            if msg.get("id") == rid:
                return msg
            best = best or msg
    return best or {"error": {"message": "no JSON-RPC payload in SSE stream"}}


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
class MCPResource:
    server: str
    uri: str
    name: str = ""
    description: str = ""
    mime_type: str = ""


@dataclass
class MCPPrompt:
    server: str
    name: str
    description: str = ""
    arguments: list[dict] = field(default_factory=list)

    @property
    def namespaced_name(self) -> str:
        return f"{_NAMESPACE_PREFIX}{self.server}__{self.name}"


@dataclass
class MCPServer:
    name: str
    command: str
    args: list[str] = field(default_factory=list)
    env: dict[str, str] = field(default_factory=dict)
    # Transport. "stdio" (default) spawns a subprocess; "http" reaches the
    # server over Streamable HTTP at ``url`` with optional ``headers``.
    transport: str = "stdio"
    url: str = ""
    headers: dict[str, str] = field(default_factory=dict)
    session_id: str = ""        # Mcp-Session-Id, set from the initialize reply
    proc: Optional[subprocess.Popen] = None
    _id_counter: itertools.count = field(default_factory=lambda: itertools.count(1))
    _lock: threading.Lock = field(default_factory=threading.Lock)
    _initialised: bool = False
    tools: list[MCPTool] = field(default_factory=list)
    last_error: str = ""

    def start(self) -> None:
        if self.transport == "http":
            return  # HTTP is connectionless — nothing to spawn
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
        if self.transport == "http":
            self.session_id = ""
            self._initialised = False
            return
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
        if self.transport == "http":
            return self._send_http(method, params)
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

    def _send_http(
        self, method: str, params: dict | None = None,
        *, is_notification: bool = False,
    ) -> dict:
        """Streamable-HTTP JSON-RPC. Handles JSON and SSE responses and the
        Mcp-Session-Id lifecycle. Notifications return ``{}``."""
        body: dict[str, Any] = {"jsonrpc": "2.0", "method": method,
                                "params": params or {}}
        rid: Any = None
        if not is_notification:
            rid = next(self._id_counter)
            body["id"] = rid
        data = json.dumps(body).encode("utf-8")
        req_headers = {
            "Content-Type": "application/json",
            "Accept": "application/json, text/event-stream",
        }
        req_headers.update(self.headers or {})
        if self.session_id:
            req_headers["Mcp-Session-Id"] = self.session_id
        req = urllib.request.Request(
            self.url, data=data, headers=req_headers, method="POST")
        with self._lock:
            try:
                with urllib.request.urlopen(req, timeout=_RPC_TIMEOUT_S) as resp:
                    sid = resp.headers.get("Mcp-Session-Id")
                    if sid:
                        self.session_id = sid
                    ctype = (resp.headers.get("Content-Type") or "").lower()
                    raw = resp.read().decode("utf-8", "replace")
            except urllib.error.HTTPError as exc:
                return {"error": {"message": f"http {exc.code}: {exc.reason}"}}
            except Exception as exc:
                return {"error": {"message": f"http request failed: {exc}"}}
        if is_notification:
            return {}
        if "text/event-stream" in ctype:
            return _extract_jsonrpc_from_sse(raw, rid)
        if not raw.strip():
            return {"error": {"message": "empty response"}}
        try:
            msg = json.loads(raw)
        except json.JSONDecodeError:
            return {"error": {"message": "invalid JSON-RPC response"}}
        return msg if isinstance(msg, dict) else {
            "error": {"message": "unexpected response shape"}}

    def _notify(self, method: str, params: dict | None = None) -> None:
        """Fire a JSON-RPC notification (no id, no response expected)."""
        if self.transport == "http":
            self._send_http(method, params, is_notification=True)
            return
        try:
            assert self.proc is not None and self.proc.stdin is not None
            msg: dict[str, Any] = {"jsonrpc": "2.0", "method": method}
            if params:
                msg["params"] = params
            self.proc.stdin.write(json.dumps(msg) + "\n")
            self.proc.stdin.flush()
        except Exception as exc:    # pragma: no cover
            self.last_error = f"{method} notify failed: {exc}"

    def initialize(self) -> bool:
        if self._initialised:
            return True
        if self.transport != "http":
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
        # Required follow-up notification (over the active transport).
        self._notify("notifications/initialized")
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

    def list_resources(self) -> list[MCPResource]:
        if not self.initialize():
            return []
        resp = self._send("resources/list")
        if "error" in resp:
            self.last_error = json.dumps(resp["error"])[:240]
            return []
        result = resp.get("result", {})
        items = result.get("resources", []) if isinstance(result, dict) else []
        out: list[MCPResource] = []
        for r in items:
            if not isinstance(r, dict):
                continue
            uri = str(r.get("uri", "")).strip()
            if not uri:
                continue
            out.append(MCPResource(
                server=self.name, uri=uri,
                name=str(r.get("name", "")),
                description=str(r.get("description", "")),
                mime_type=str(r.get("mimeType", "")),
            ))
        return out

    def read_resource(self, uri: str) -> str:
        if not self.initialize():
            return ""
        resp = self._send("resources/read", {"uri": uri})
        if "error" in resp:
            self.last_error = json.dumps(resp["error"])[:240]
            return json.dumps({"error": resp["error"]})
        result = resp.get("result", {})
        contents = result.get("contents", []) if isinstance(result, dict) else []
        texts: list[str] = []
        for c in contents:
            if not isinstance(c, dict):
                continue
            if c.get("text") is not None:
                texts.append(str(c.get("text", "")))
            elif c.get("blob") is not None:
                texts.append(
                    f"[binary {c.get('mimeType', 'application/octet-stream')} "
                    f"resource {c.get('uri', uri)} — base64 omitted]"
                )
        return "\n".join(texts)

    def list_prompts(self) -> list[MCPPrompt]:
        if not self.initialize():
            return []
        resp = self._send("prompts/list")
        if "error" in resp:
            self.last_error = json.dumps(resp["error"])[:240]
            return []
        result = resp.get("result", {})
        items = result.get("prompts", []) if isinstance(result, dict) else []
        out: list[MCPPrompt] = []
        for p in items:
            if not isinstance(p, dict):
                continue
            name = str(p.get("name", "")).strip()
            if not name:
                continue
            out.append(MCPPrompt(
                server=self.name, name=name,
                description=str(p.get("description", "")),
                arguments=list(p.get("arguments", []) or []),
            ))
        return out

    def get_prompt(self, name: str, arguments: dict | None = None) -> str:
        if not self.initialize():
            return ""
        resp = self._send("prompts/get",
                          {"name": name, "arguments": arguments or {}})
        if "error" in resp:
            self.last_error = json.dumps(resp["error"])[:240]
            return json.dumps({"error": resp["error"]})
        result = resp.get("result", {})
        if not isinstance(result, dict):
            return ""
        parts: list[str] = []
        for m in result.get("messages", []) or []:
            if not isinstance(m, dict):
                continue
            role = str(m.get("role", "")).strip()
            text = _flatten_content(m.get("content"))
            if text:
                parts.append(f"{role}: {text}" if role else text)
        return "\n\n".join(parts)


def _server_from_config(name: str, cfg: dict) -> MCPServer:
    """Build an MCPServer from one config entry, picking the transport.

    HTTP when ``type`` is http/sse/streamable-http or a ``url`` is present;
    otherwise a stdio subprocess.
    """
    ttype = str(cfg.get("type", "")).strip().lower()
    if ttype in _HTTP_TYPES or cfg.get("url"):
        return MCPServer(
            name=name, command="", transport="http",
            url=str(cfg.get("url", "")),
            headers={k: str(v) for k, v in (cfg.get("headers") or {}).items()},
        )
    return MCPServer(
        name=name,
        command=str(cfg.get("command", "")),
        args=list(cfg.get("args", []) or []),
        env={k: str(v) for k, v in (cfg.get("env") or {}).items()},
    )


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
            self.servers[name] = _server_from_config(name, cfg)
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

    def discover_resources(self) -> list[MCPResource]:
        out: list[MCPResource] = []
        for srv in self.servers.values():
            try:
                out.extend(srv.list_resources())
            except Exception:   # pragma: no cover
                pass
        return out

    def discover_prompts(self) -> list[MCPPrompt]:
        out: list[MCPPrompt] = []
        for srv in self.servers.values():
            try:
                out.extend(srv.list_prompts())
            except Exception:   # pragma: no cover
                pass
        return out

    def read_resource(self, server: str, uri: str) -> str:
        srv = self.servers.get(server)
        if srv is None:
            return json.dumps({"error": f"unknown MCP server: {server!r}"})
        return srv.read_resource(uri)

    def get_prompt(self, namespaced: str,
                   arguments: dict | None = None) -> str:
        if not namespaced.startswith(_NAMESPACE_PREFIX):
            return json.dumps({"error": f"not an MCP prompt: {namespaced!r}"})
        rest = namespaced[len(_NAMESPACE_PREFIX):]
        if "__" not in rest:
            return json.dumps({"error": f"malformed MCP name: {namespaced!r}"})
        server_name, _, prompt_name = rest.partition("__")
        srv = self.servers.get(server_name)
        if srv is None:
            return json.dumps({"error": f"unknown MCP server: {server_name!r}"})
        return srv.get_prompt(prompt_name, arguments or {})

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
    "MCPResource",
    "MCPPrompt",
    "MCPServer",
    "MCPRegistry",
    "get_registry",
    "reset_registry",
]
