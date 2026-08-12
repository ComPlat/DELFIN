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
just leaves that server's tool set empty, and so does silence — every
stdio read is bounded by ``_RPC_TIMEOUT_S``, with the reason kept in
the server's ``last_error``.
"""

from __future__ import annotations

import itertools
import json
import os
import queue
import re
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional


_DEFAULT_PROTOCOL_VERSION = "2025-06-18"
_RPC_TIMEOUT_S = 15.0
# Lines a stdio server may run ahead of the client before it is made to wait.
_STDOUT_QUEUE_LINES = 1000
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
#
# ``delfin-ops`` is here because the role prompt MANDATES it: the first tool
# call for any ORCA question must be ``mcp__delfin-ops__extract_*`` or
# ``parse_orca_output``. Only ``delfin-tools`` was registered, and the config
# the dashboard auto-injected for the ops server was written to a different
# file (``~/.delfin/mcp_ops_config.json``) under a different top-level key
# (``mcpServers``) that reaches only the legacy CLI subprocess — this registry
# reads ``~/.delfin/mcp_servers.json`` / ``servers``. So the mandated first
# move answered "unknown MCP server: 'delfin-ops'" and the fallback the prompt
# forbids was the only thing that worked.
#
# ``delfin-docs`` is deliberately NOT here. Every one of its seven tools
# (search_docs, read_section, list_docs, list_sections, search_calcs,
# get_calc_info, calc_summary) is already advertised natively under exactly
# those names and served from the same index by ``_DocToolExecutor``.
# Registering it would pay for seven duplicate schemas on every request and
# give the model two names for one thing.
#
# ``sys.executable``, not ``"python"``: the interpreter running DELFIN is the
# one with DELFIN importable, and on a host that only has ``python3`` the
# literal spawns nothing — which, before the reachability fix below, was
# indistinguishable from a healthy server with no tools.
_BUILTIN_SERVERS: dict[str, dict] = {
    "delfin-tools": {
        "command": sys.executable or "python3",
        "args": ["-m", "delfin.tools.mcp_server"],
        "enabled": True,
    },
    "delfin-ops": {
        "command": sys.executable or "python3",
        "args": ["-m", "delfin.ops_server"],
        "enabled": True,
    },
}


def _load_configs(workspace: Path | None) -> dict[str, dict]:
    """Merge built-in defaults + user-global + project-scoped MCP configs."""
    out: dict[str, dict] = {
        name: dict(cfg) for name, cfg in _BUILTIN_SERVERS.items()
    }
    paths = [(_user_config_path(), True)]
    if workspace:
        paths.append((_project_config_path(workspace), False))
    for path, is_user in paths:
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
            if (not is_user and name in _BUILTIN_SERVERS
                    and cfg.get("enabled", True)):
                # A PROJECT config may DISABLE a builtin -- that is a
                # tightening, and a project that does not want DELFIN's
                # own tools is entitled to say so -- but it may not
                # REDEFINE one. `delfin-tools` is DELFIN's own server, so
                # an entry naming it silently redirects every
                # mcp__delfin-tools__* call the agent makes, and the /mcp
                # listing reads only the user-global file, so it cannot
                # even be seen. The person replacing it in their own
                # settings is the case the override exists for; a folder
                # doing it behind their back is not.
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
    # Lines drained off the subprocess by ``_reader_lines``; ``None`` is the
    # EOF sentinel. Bound to the process they came from, so a restart cannot
    # be served stale output from the previous one.
    _lines: Optional["queue.Queue"] = None
    _lines_proc: Optional[subprocess.Popen] = None

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
        # Drop the reader's queue with the process it belonged to, so a
        # restart cannot be handed lines the previous process left behind.
        self._lines = None
        self._lines_proc = None

    def _reader_lines(self) -> "queue.Queue":
        """Queue of stdout lines, filled by a daemon thread we may abandon.

        The deadline in ``_send`` has to be enforceable while NOTHING is
        arriving, and a blocking ``readline()`` cannot be interrupted: it
        returns on a newline or on EOF and on nothing else, so a clock checked
        around it bounds only servers that keep talking. The alternatives do
        not survive the pipe being opened in text mode (``text=True`` in
        ``start``): ``selectors`` on the raw fd can report "not ready" while
        the ``TextIOWrapper`` already holds a complete buffered line, and
        putting the fd in non-blocking mode makes that same wrapper raise on
        a short read instead of waiting. Moving the blocking read to a thread
        keeps the stream exactly as it was and lets the caller wait on
        ``Queue.get(timeout=...)`` -- on something it is allowed to give up
        on -- instead of on the kernel. The thread is a daemon and holds no
        lock, so a server that never speaks again costs one parked thread
        until the process is stopped rather than the turn it was blocking.
        """
        proc = self.proc
        if self._lines is not None and self._lines_proc is proc:
            return self._lines
        # Bounded on purpose. Reading only inside ``_send`` meant an idle
        # client left the server's output in the OS pipe buffer, which stops
        # a server that floods; draining into an unbounded queue would remove
        # that brake and let one chatty server grow the agent's memory instead
        # of its own. Full queue -> the pump blocks -> the pipe fills -> the
        # server waits, exactly as before.
        lines: "queue.Queue" = queue.Queue(maxsize=_STDOUT_QUEUE_LINES)
        self._lines = lines
        self._lines_proc = proc
        stdout = proc.stdout if proc is not None else None

        def _pump() -> None:
            try:
                while True:
                    raw = stdout.readline()     # type: ignore[union-attr]
                    if not raw:
                        break
                    lines.put(raw)
            except Exception:   # pragma: no cover - pipe torn down under us
                pass
            finally:
                lines.put(None)                 # EOF sentinel
        threading.Thread(
            target=_pump, name=f"mcp-reader-{self.name}", daemon=True,
        ).start()
        return lines

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
            lines = self._reader_lines()
            while True:
                remaining = t_end - time.monotonic()
                if remaining <= 0:
                    break
                try:
                    raw = lines.get(timeout=remaining)
                except queue.Empty:
                    break
                if raw is None:
                    # EOF is terminal for this process: put the sentinel back
                    # so every later call sees it too, rather than the first
                    # one consuming it and the rest waiting out the budget on
                    # a pipe that can no longer produce anything.
                    lines.put(None)
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
                # ignore notifications / unrelated messages -- including a
                # reply that a PREVIOUS call already gave up on, which the
                # id check discards here.
        # Recorded, not just returned: a server that accepts a request and
        # never answers looks identical to a healthy one from the outside,
        # and this field is the only place the reason survives the empty
        # tool list that discovery reports.
        self.last_error = (
            f"rpc timeout: no reply to {method} within {_RPC_TIMEOUT_S:g}s"
        )
        return {"error": {"message": self.last_error}}

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


# Discovery asks every configured server three questions (tools,
# resources, prompts) one after another, and each question can sit on the
# per-call deadline. A server that is dead still costs the full wait, and
# `initialize` sends its own request first, so one unreachable entry can
# hold the pass for twice the RPC deadline before anything moves on. Five
# of them, serially, is minutes -- and every failure was swallowed, so
# the only symptom was a turn that would not start.
#
# The per-call deadline cannot see this: it is doing its job each time.
# What was missing is a ceiling on the WHOLE pass. Once it is spent the
# remaining servers are not asked, and which ones they were is reported
# rather than dropped.
_DISCOVERY_BUDGET_S = 30.0


@dataclass
class _DiscoveryBudget:
    """Wall-clock left for one discovery pass, plus what it skipped."""

    total_s: float = _DISCOVERY_BUDGET_S
    started: float = field(default_factory=time.monotonic)
    skipped: list[str] = field(default_factory=list)
    failed: list[str] = field(default_factory=list)

    def spent(self) -> float:
        return time.monotonic() - self.started

    def exhausted(self) -> bool:
        return self.spent() >= self.total_s

    def note_skipped(self, name: str) -> None:
        if name not in self.skipped:
            self.skipped.append(name)

    def note_failed(self, name: str, exc: BaseException | str) -> None:
        entry = f"{name}: {exc}"[:200]
        if entry not in self.failed:
            self.failed.append(entry)

    def report(self) -> dict:
        return {
            "seconds": round(self.spent(), 2),
            "budget_s": self.total_s,
            "exhausted": self.exhausted(),
            "skipped": list(self.skipped),
            "failed": list(self.failed),
        }


@dataclass
class MCPRegistry:
    servers: dict[str, MCPServer] = field(default_factory=dict)
    workspace: Optional[Path] = None
    loaded: bool = False
    last_discovery: dict = field(default_factory=dict)

    def load(self, workspace: Path | None = None) -> None:
        configs = _load_configs(workspace)
        for name, cfg in configs.items():
            if name in self.servers:
                continue
            self.servers[name] = _server_from_config(name, cfg)
        self.workspace = Path(workspace) if workspace else self.workspace
        self.loaded = True

    def discover_all(self, budget: "_DiscoveryBudget | None" = None
                     ) -> list[MCPTool]:
        own = budget is None
        budget = budget or _DiscoveryBudget()
        tools: list[MCPTool] = []
        for name, srv in self.servers.items():
            if budget.exhausted():
                budget.note_skipped(name)
                continue
            try:
                tools.extend(srv.list_tools())
            except Exception as exc:
                budget.note_failed(name, exc)
            else:
                if getattr(srv, "last_error", ""):
                    budget.note_failed(name, srv.last_error)
        if own:
            self.last_discovery = budget.report()
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

    def discover_resources(self, budget: "_DiscoveryBudget | None" = None
                           ) -> list[MCPResource]:
        own = budget is None
        budget = budget or _DiscoveryBudget()
        out: list[MCPResource] = []
        for name, srv in self.servers.items():
            if budget.exhausted():
                budget.note_skipped(name)
                continue
            try:
                out.extend(srv.list_resources())
            except Exception as exc:
                budget.note_failed(name, exc)
        if own:
            self.last_discovery = budget.report()
        return out

    def discover_prompts(self, budget: "_DiscoveryBudget | None" = None
                         ) -> list[MCPPrompt]:
        own = budget is None
        budget = budget or _DiscoveryBudget()
        out: list[MCPPrompt] = []
        for name, srv in self.servers.items():
            if budget.exhausted():
                budget.note_skipped(name)
                continue
            try:
                out.extend(srv.list_prompts())
            except Exception as exc:
                budget.note_failed(name, exc)
        if own:
            self.last_discovery = budget.report()
        return out

    def discover_everything(self) -> tuple[list, list, list]:
        """One pass, one budget. Asking three times with three separate
        ceilings would let a broken configuration cost three times the
        ceiling, which is the thing the ceiling exists to prevent."""
        budget = _DiscoveryBudget()
        tools = self.discover_all(budget)
        resources = self.discover_resources(budget)
        prompts = self.discover_prompts(budget)
        self.last_discovery = budget.report()
        return tools, resources, prompts

    def discovery_notice(self) -> str:
        """One line for the user when a pass did not reach everything.

        Empty when it did. A discovery that quietly returned fewer tools
        than the config declares is indistinguishable from a config with
        fewer tools in it, and that is the state this reports."""
        rep = self.last_discovery or {}
        if not rep.get("skipped") and not rep.get("failed"):
            return ""
        parts = []
        if rep.get("skipped"):
            parts.append(
                f"{len(rep['skipped'])} server(s) were not asked "
                f"({', '.join(rep['skipped'])}) — the {rep.get('budget_s')}s "
                f"discovery budget was spent on the ones before them")
        if rep.get("failed"):
            parts.append(f"{len(rep['failed'])} did not answer: "
                         + "; ".join(rep["failed"]))
        return ("MCP discovery took "
                f"{rep.get('seconds')}s and did not complete: "
                + ". ".join(parts) + ".")

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


def count_live_tools(registry: "MCPRegistry") -> int:
    """How many tools the servers in ``registry`` actually offer.

    The dashboard used ``len(getattr(registry, "tools", {}) or {})``, and
    MCPRegistry has no ``tools`` attribute -- the default fired every
    time, so the number it printed beside the word "live" was zero for a
    healthy server exactly as much as for a dead one. Tools live on each
    server and appear only once someone asks.

    Asking is the point: it is the difference between reporting the
    config back and reporting what answered.
    """
    total = 0
    for server in (getattr(registry, "servers", None) or {}).values():
        try:
            total += len(server.list_tools() or ())
        except Exception:
            continue
    return total


def unreachable_servers(registry: "MCPRegistry") -> list[str]:
    """``name: reason`` for every server that could not be asked.

    The verdict comes from ``last_error`` AFTER the call, not from an
    exception. ``list_tools`` fails CLOSED by contract — a missing binary,
    a server that exits on startup, a refused HTTP request, a timeout: each
    records the cause in ``last_error`` and returns ``[]``. Deciding by
    ``except`` therefore never fired, and the list was empty for a dead
    server exactly as much as for a healthy one, while the doctor printed
    "configured + reachable" and the dashboard printed "0 tools, no
    problems". Measured on a server with a missing binary: ``list_tools ->
    []``, ``last_error -> "command not found: …"``, ``unreachable_servers
    -> []``.

    A server that answered WITH tools is reachable whatever an earlier call
    left behind, so a stale error is cleared rather than reported. A server
    that answered with an empty list and no error is not a failure — it is
    a server that offers nothing.
    """
    out: list[str] = []
    for name, server in (getattr(registry, "servers", None) or {}).items():
        tools: Any = None
        reason = ""
        try:
            tools = server.list_tools()
        except Exception as exc:
            reason = str(exc)
        if tools:
            continue
        try:
            recorded = str(getattr(server, "last_error", "") or "")
        except Exception:
            recorded = ""
        reason = recorded or reason
        if reason:
            out.append(f"{name}: {reason}")
    return out


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
