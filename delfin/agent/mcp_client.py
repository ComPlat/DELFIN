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

The user's own config is always read. A WORKSPACE's config is read only
when the user has explicitly trusted that directory for MCP servers —
see ``workspace_trust``. A server definition is the most powerful thing
a folder can ship: it is spawned with the parent environment while the
tool surface is being assembled, before the model emits a token, so a
repository the user merely checked out could run a command of its
choosing by containing one file. What a lack of trust withheld is
reported on the registry rather than dropped.

A registry is cached per resolved WORKSPACE — the first
``get_registry(ws)`` loads that workspace's configured servers and
further calls with the same workspace reuse the running processes /
sessions. Keyed by nothing, as it was, a second session in another repo
was silently handed the first repo's servers with the first repo's
environment, and its own project config was never read. Servers fail
closed: a crash during discovery
just leaves that server's tool set empty, and so does silence — every
stdio read is bounded by ``_RPC_TIMEOUT_S``, with the reason kept in
the server's ``last_error``. A stdio child's stderr is drained into a
bounded ring and, where the process is found dead, joined onto that
reason — discarding it meant a built-in that could not import its own
dependencies reached the caller as "server not running", which reads as
a server nobody configured rather than a broken install.
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
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional


_DEFAULT_PROTOCOL_VERSION = "2025-06-18"
_RPC_TIMEOUT_S = 15.0
# How much of a dead server's stderr is kept. A failing server explains
# itself there and the explanation has to survive, but a server is an
# arbitrary program and may write without end, so the account is held in a
# ring rather than a buffer: forty lines of four hundred characters is at
# most ~16 kB per server, and holds a Python traceback whole — the import
# error this was built for is eight lines.
_STDERR_KEEP_LINES = 40
_STDERR_LINE_CHARS = 400
# What is HELD and what is REPORTED are two different bounds. The ring above
# is sized so a traceback survives; this is what a doctor line or a /mcp
# listing is willing to print of it, so a server whose dying words fill the
# ring does not fill the terminal too.
_STDERR_REPORT_CHARS = 1200
# Only ever waited out on a process that has ALREADY exited, whose stderr
# is therefore closed and whose drain thread ends at EOF. It is a handover,
# not a wait on a live server.
_STDERR_JOIN_S = 0.5
# A child that dies closes stdout a moment BEFORE its exit status is
# reapable, and the reader sees the EOF first. Measured on a server that
# exits on an import error: reporting straight from the EOF gave "EOF from
# server" and lost the traceback that had already been captured. This is
# how long the reader waits for the exit that EOF announced before deciding
# the stream was merely closed by a server that is still running.
_EOF_EXIT_GRACE_S = 0.25
_NAMESPACE_PREFIX = "mcp__"
_HTTP_TYPES = {"http", "sse", "streamable-http", "streamablehttp"}


def _user_config_path() -> Path:
    return Path.home() / ".delfin" / "mcp_servers.json"


def _project_config_path(workspace: Path) -> Path:
    """Where a workspace's server config lives.

    Derived from the trust registry rather than spelled out here, so the
    file that SPAWNS processes and the file that needs a trust decision
    (and a protected glob) cannot drift apart.
    """
    from . import workspace_trust as _trust
    rel = _trust.get_kind(_trust.KIND_MCP_SERVERS).relative_paths[0]
    return Path(workspace) / rel


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


def _load_configs_with_sources(
    workspace: Path | None,
) -> tuple[dict[str, dict], dict[str, str], str]:
    """``(configs, source-per-server, trust notice)``.

    The trust gate lives HERE, in the loader, and not at the callers:
    the tool-surface assembly, ``/mcp reload``, the doctor and every
    future caller reach a workspace's config only through this function,
    so the question is asked once and cannot be forgotten by the next
    caller added.

    The source of each entry is carried out with it because ``/mcp``
    listed only the user's own file — a server a workspace added was
    spawned and never appeared in the listing at all.
    """
    out: dict[str, dict] = {
        name: dict(cfg) for name, cfg in _BUILTIN_SERVERS.items()
    }
    sources: dict[str, str] = {name: "built-in" for name in _BUILTIN_SERVERS}
    paths = [(_user_config_path(), True)]
    notice = ""
    if workspace:
        from . import workspace_trust as _trust
        decision = _trust.gate(_trust.KIND_MCP_SERVERS, workspace)
        notice = decision.notice
        paths.extend((p, False) for p in decision.paths)
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
                sources[name] = str(path)
            else:
                out.pop(name, None)             # explicit disable (also disables a builtin)
                sources.pop(name, None)
    return out, sources, notice


def _load_configs(workspace: Path | None) -> dict[str, dict]:
    """Merge built-in defaults + user-global + trusted project MCP configs."""
    return _load_configs_with_sources(workspace)[0]


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
    # Guards the process handle, the id draw and the stdin write — and
    # nothing else. It used to be held across the whole reply wait as
    # well, and the registry hands one server object to the parent turn
    # and to every background sub-agent, so a server that accepted and
    # never answered cost each waiter the full deadline one after the
    # other: four callers, four times fifteen seconds.
    _lock: threading.Lock = field(default_factory=threading.Lock)
    _initialised: bool = False
    tools: list[MCPTool] = field(default_factory=list)
    last_error: str = ""
    # Reply slots by request id, filled by the reader thread. Bound to the
    # process they belong to, so a restart cannot be served stale output
    # from the previous one.
    _waiters: dict = field(default_factory=dict)
    _waiters_lock: threading.Lock = field(default_factory=threading.Lock)
    _reader_proc: Optional[subprocess.Popen] = None
    _closed_reason: str = ""
    # The last of the child's stderr, kept so a server that died can say
    # what killed it. Bounded by the deque itself, not by trust in the
    # child: see ``_drain_stderr``.
    _stderr_tail_lines: Any = field(
        default_factory=lambda: deque(maxlen=_STDERR_KEEP_LINES))
    _stderr_thread: Optional[threading.Thread] = None

    def start(self) -> None:
        if self.transport == "http":
            return  # HTTP is connectionless — nothing to spawn
        with self._lock:
            if self.proc is not None and self.proc.poll() is None:
                return
            env = dict(os.environ)
            env.update(self.env or {})
            try:
                # stderr is PIPED, not discarded. It used to be DEVNULL,
                # which threw away the only account a dying server ever
                # gives: an unconditional ``mcp.server.fastmcp`` import in
                # the three built-ins met an SDK that had moved the module,
                # and every one of them exited 1 with a ModuleNotFoundError
                # on stderr — while the caller was told "server not
                # running", which is true, useless, and looks like a
                # configuration mistake rather than a broken install.
                self.proc = subprocess.Popen(
                    [self.command, *self.args],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env=env, text=True, bufsize=1,
                    # Its own process group. A stdio server is long-lived
                    # across turns, and under a terminal front-end it would
                    # otherwise sit in the foreground group and take every
                    # Ctrl+C the user aims at the agent — so interrupting
                    # one turn would tear down every configured server for
                    # the rest of the session. The dashboard never had a
                    # controlling terminal, so nothing surfaced this.
                    # stop() still terminates it explicitly.
                    start_new_session=True,
                )
                self._closed_reason = ""
                self._reader_proc = None
                # A restart must not be able to report the PREVIOUS
                # process's dying words as this one's. A fresh ring rather
                # than a cleared one, because the previous process's drain
                # thread may still be emptying its pipe: it holds the
                # object it was started with, so its late lines land in a
                # ring nobody reads instead of in this process's account.
                self._stderr_tail_lines = deque(maxlen=_STDERR_KEEP_LINES)
                self._drain_stderr(self.proc)
            except FileNotFoundError as exc:
                self.last_error = f"command not found: {self.command} ({exc})"
                self.proc = None
            except Exception as exc:    # pragma: no cover
                self.last_error = f"start failed: {exc}"
                self.proc = None

    def stop(self) -> None:
        """Shut the server down. Safe against a concurrent ``_send``.

        The handle is dropped under the lock — clearing it outside meant a
        caller could be dereferencing ``self.proc.stdin`` at the moment it
        became ``None``. Terminating is done after the lock is released:
        waiting two seconds for a child to die is exactly the blocking
        call this lock must not span."""
        if self.transport == "http":
            with self._lock:
                self.session_id = ""
                self._initialised = False
            return
        with self._lock:
            proc, self.proc = self.proc, None
            self._initialised = False
            self._reader_proc = None
            self._closed_reason = "server stopped"
        if proc is None:
            return
        try:
            proc.terminate()
            proc.wait(timeout=2)
        except subprocess.TimeoutExpired:
            proc.kill()
        except Exception:    # pragma: no cover
            pass
        # Anyone still waiting is waiting on a process that is gone.
        self._release_waiters("server stopped")

    def _release_waiters(self, reason: str) -> None:
        """Hand every outstanding caller the same ending. Never raises."""
        with self._waiters_lock:
            slots = list(self._waiters.values())
            self._waiters.clear()
        for slot in slots:
            try:
                slot.put_nowait({"error": {"message": reason}})
            except Exception:   # pragma: no cover - slot already answered
                pass

    def _drain_stderr(self, proc: subprocess.Popen) -> None:
        """Keep the child's stderr moving, and keep the last of it.

        Draining cannot be skipped once the stream is a pipe: a pipe nobody
        reads fills at the OS buffer — 64 kB on Linux — and the child then
        blocks forever on its next write. So the choice was never "pipe or
        not"; it was between discarding the stream and reading it
        continuously, and discarding it is what made a one-line import
        error indistinguishable from a server that was merely absent.

        Whatever arrives is kept in a bounded ring, so a server that logs
        without end costs a fixed amount of memory rather than growing one.
        Nothing here judges the content: a healthy server logs to stderr as
        a matter of course, and this thread runs for those too. The lines
        become an ERROR only where a failure has already been established
        by other means — see ``_exit_reason``.
        """
        stderr = proc.stderr
        if stderr is None:      # pragma: no cover - only if PIPE was refused
            return
        # Bound to the ring current at spawn, not to the attribute: see the
        # note in ``start`` about a restart overtaking a thread that is
        # still draining the process before it.
        ring = self._stderr_tail_lines

        def _pump() -> None:
            try:
                for raw in stderr:
                    line = raw.rstrip("\n")[:_STDERR_LINE_CHARS].strip()
                    if line:
                        ring.append(line)
            except Exception:   # pragma: no cover - pipe torn down under us
                pass
            finally:
                # The read end is this thread's to close; leaving it open
                # would cost a file descriptor per server started, where
                # DEVNULL cost none.
                try:
                    stderr.close()
                except Exception:   # pragma: no cover - already closed
                    pass
        thread = threading.Thread(
            target=_pump, name=f"mcp-stderr-{self.name}", daemon=True,
        )
        self._stderr_thread = thread
        thread.start()

    def _exit_reason(self, proc: Optional[subprocess.Popen]) -> str:
        """Why this server is not answering, in the child's own words.

        Called only from a path that has already established the process is
        gone. The drain thread is given a moment to finish first: the child
        has exited, so its stderr is closed and the thread is ending at EOF
        anyway, but without the handover the parent can look at the ring
        microseconds before the last line lands in it and report the same
        empty "not running" this exists to replace.

        The child's output is attached to an EXIT STATUS or not at all. With
        no handle there is no process to attribute the ring to — the handle
        is dropped by a deliberate ``stop``, not by a death — and a healthy
        server's startup banner offered as the cause of a failure would be
        worse than offering no cause, which is where this started.
        """
        if proc is None or proc.poll() is None:
            return "server not running"
        thread = self._stderr_thread
        if thread is not None and thread.is_alive():
            thread.join(timeout=_STDERR_JOIN_S)
        head = f"server exited with code {proc.poll()}"
        tail = " | ".join(self._stderr_tail_lines)
        if len(tail) > _STDERR_REPORT_CHARS:
            # The END is kept, not the beginning: a traceback's last line is
            # the one that names the error, and the frames above it are the
            # part a reader can do without.
            tail = "…" + tail[-_STDERR_REPORT_CHARS:]
        return f"{head}: {tail}" if tail else head

    def _record_rpc_error(self, err: Any) -> None:
        """Keep the fuller of the two accounts of a failed call.

        ``_send`` records the reason itself wherever it knows more than the
        JSON-RPC envelope can carry — a dead child's stderr, a deadline that
        expired. Re-encoding that here would JSON-escape it and cut it at
        240 characters, and for a traceback 240 characters is precisely the
        part BEFORE the one line that names the error. So a message ``_send``
        has already recorded is left as it stands. An error that came from a
        live SERVER is not something ``_send`` knows anything about, and is
        stored bounded, as it always was.
        """
        msg = err.get("message") if isinstance(err, dict) else None
        if msg and msg == self.last_error:
            return
        self.last_error = json.dumps(err)[:240]

    def _ensure_reader(self, proc: subprocess.Popen) -> None:
        """Start the daemon that reads this process's stdout, once.

        The deadline in ``_send`` has to be enforceable while NOTHING is
        arriving, and a blocking ``readline()`` cannot be interrupted: it
        returns on a newline or on EOF and on nothing else, so a clock checked
        around it bounds only servers that keep talking. The alternatives do
        not survive the pipe being opened in text mode (``text=True`` in
        ``start``): ``selectors`` on the raw fd can report "not ready" while
        the ``TextIOWrapper`` already holds a complete buffered line, and
        putting the fd in non-blocking mode makes that same wrapper raise on
        a short read instead of waiting. Moving the blocking read to a thread
        keeps the stream exactly as it was and lets each caller wait on a
        slot of its own -- on something it is allowed to give up on --
        instead of on the kernel. The thread is a daemon and holds no server
        lock, so a server that never speaks again costs one parked thread
        until the process is stopped rather than the turn it was blocking.

        This thread also does the DISPATCH: it matches each reply to the
        request that is waiting for it. Draining a shared queue inside
        ``_send`` meant the reply belonging to caller B could only be seen
        by whoever held the lock, which is why the wait had to hold it.
        Memory stays bounded without the old queue cap for a better
        reason than before: a line nobody is waiting for is dropped where
        it is read rather than buffered.
        """
        if self._reader_proc is proc:
            return
        self._reader_proc = proc
        stdout = proc.stdout

        def _pump() -> None:
            try:
                while True:
                    raw = stdout.readline()     # type: ignore[union-attr]
                    if not raw:
                        break
                    self._deliver(raw)
            except Exception:   # pragma: no cover - pipe torn down under us
                pass
            finally:
                # EOF is terminal for this process: recorded so every later
                # call sees it at once, rather than each waiting out the
                # budget on a pipe that can no longer produce anything.
                #
                # WHY it ended is worth more than the fact that it did. A
                # child that closed stdout because it DIED can still account
                # for itself on stderr; one that merely closed the stream and
                # kept running cannot, and must not be waited on — so the
                # exit is given a brief grace and the richer reason is used
                # only if it actually arrives.
                reason = "EOF from server"
                try:
                    proc.wait(timeout=_EOF_EXIT_GRACE_S)
                except Exception:
                    pass
                if proc.poll() is not None:
                    reason = self._exit_reason(proc)
                    # Recorded here as well, because this thread reaches the
                    # death before any caller does. Without it the FIRST
                    # discovery pass — the one that actually runs at startup
                    # — saw only the envelope this reason travels in, and
                    # re-encoded it into the 240-character bound meant for a
                    # live server's JSON-RPC errors.
                    self.last_error = reason
                if self._reader_proc is proc:
                    self._closed_reason = reason
                self._release_waiters(reason)
        threading.Thread(
            target=_pump, name=f"mcp-reader-{self.name}", daemon=True,
        ).start()

    def _deliver(self, raw: str) -> None:
        """Route one stdout line to the caller waiting on its id."""
        raw = raw.strip()
        if not raw:
            return
        try:
            msg = json.loads(raw)
        except json.JSONDecodeError:
            return
        if not isinstance(msg, dict):
            return
        rid = msg.get("id")
        if rid is None:
            return          # a notification: nobody is waiting for it
        with self._waiters_lock:
            slot = self._waiters.pop(rid, None)
        if slot is None:
            return          # a reply the caller already gave up on
        try:
            slot.put_nowait(msg)
        except Exception:   # pragma: no cover
            pass

    def _send(self, method: str, params: dict | None = None) -> dict:
        if self.transport == "http":
            return self._send_http(method, params)
        req_id: Any
        slot: "queue.Queue" = queue.Queue(maxsize=1)
        # The lock covers the write and the id registration. The WAIT is
        # deliberately outside it: it is this caller's wait, and holding
        # the lock across it made every other caller of the same server
        # queue up behind the full deadline.
        with self._lock:
            proc = self.proc
            if proc is None or proc.poll() is not None:
                # Recorded as well as returned: ``list_tools`` fails closed
                # and hands back ``[]``, so ``last_error`` is the only place
                # the reason survives to reach ``unreachable_servers``, the
                # doctor and the /mcp listing.
                self.last_error = self._exit_reason(proc)
                return {"error": {"message": self.last_error}}
            if self._closed_reason:
                return {"error": {"message": self._closed_reason}}
            self._ensure_reader(proc)
            req_id = next(self._id_counter)
            with self._waiters_lock:
                self._waiters[req_id] = slot
            line = json.dumps({
                "jsonrpc": "2.0",
                "id": req_id,
                "method": method,
                "params": params or {},
            }) + "\n"
            try:
                assert proc.stdin is not None
                proc.stdin.write(line)
                proc.stdin.flush()
            except (BrokenPipeError, OSError) as exc:
                with self._waiters_lock:
                    self._waiters.pop(req_id, None)
                return {"error": {"message": f"write failed: {exc}"}}
        try:
            return slot.get(timeout=_RPC_TIMEOUT_S)
        except queue.Empty:
            pass
        with self._waiters_lock:
            self._waiters.pop(req_id, None)
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
        Mcp-Session-Id lifecycle. Notifications return ``{}``.

        The lock covers the id draw and the session-id read/write, not the
        request: one server object is shared by the whole process, and a
        slow endpoint held every other caller for the full timeout each,
        one after another."""
        body: dict[str, Any] = {"jsonrpc": "2.0", "method": method,
                                "params": params or {}}
        rid: Any = None
        req_headers = {
            "Content-Type": "application/json",
            "Accept": "application/json, text/event-stream",
        }
        req_headers.update(self.headers or {})
        with self._lock:
            if not is_notification:
                rid = next(self._id_counter)
                body["id"] = rid
            if self.session_id:
                req_headers["Mcp-Session-Id"] = self.session_id
        data = json.dumps(body).encode("utf-8")
        req = urllib.request.Request(
            self.url, data=data, headers=req_headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=_RPC_TIMEOUT_S) as resp:
                sid = resp.headers.get("Mcp-Session-Id")
                if sid:
                    with self._lock:
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
            msg: dict[str, Any] = {"jsonrpc": "2.0", "method": method}
            if params:
                msg["params"] = params
            # Under the lock like every other write: interleaving with a
            # concurrent _send would put two half-lines on one pipe.
            with self._lock:
                proc = self.proc
                assert proc is not None and proc.stdin is not None
                proc.stdin.write(json.dumps(msg) + "\n")
                proc.stdin.flush()
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
            self._record_rpc_error(resp["error"])
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
            self._record_rpc_error(resp["error"])
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
            self._record_rpc_error(resp["error"])
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
            self._record_rpc_error(resp["error"])
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
            self._record_rpc_error(resp["error"])
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
            self._record_rpc_error(resp["error"])
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
    # Config file each server came from, by name. ``/mcp`` read only the
    # user's own file, so a server a workspace contributed was spawned
    # and never listed — invisible in exactly the case that matters.
    sources: dict[str, str] = field(default_factory=dict)
    # What a lack of trust withheld, in one sentence, or "".
    trust_notice: str = ""

    def load(self, workspace: Path | None = None) -> None:
        configs, sources, notice = _load_configs_with_sources(workspace)
        for name, cfg in configs.items():
            if name in self.servers:
                continue
            self.servers[name] = _server_from_config(name, cfg)
        self.sources.update(sources)
        self.trust_notice = notice
        self.workspace = Path(workspace) if workspace else self.workspace
        self.loaded = True

    def source_of(self, name: str) -> str:
        return self.sources.get(name, "")

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
        fewer tools in it, and that is the state this reports.

        The trust notice rides along: servers a workspace declared and
        the user has not trusted are the other reason the surface is
        smaller than the configuration reads."""
        rep = self.last_discovery or {}
        if not rep.get("skipped") and not rep.get("failed"):
            return self.trust_notice
        parts = []
        if rep.get("skipped"):
            parts.append(
                f"{len(rep['skipped'])} server(s) were not asked "
                f"({', '.join(rep['skipped'])}) — the {rep.get('budget_s')}s "
                f"discovery budget was spent on the ones before them")
        if rep.get("failed"):
            parts.append(f"{len(rep['failed'])} did not answer: "
                         + "; ".join(rep["failed"]))
        text = ("MCP discovery took "
                f"{rep.get('seconds')}s and did not complete: "
                + ". ".join(parts) + ".")
        return f"{text} {self.trust_notice}".strip()

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
# One registry PER WORKSPACE. It used to be a single global keyed by
# nothing: whoever called first fixed the configuration, and a second
# session in another repo silently got the first repo's servers, spawned
# with the first repo's environment — its own project config was never
# read and there was no symptom to see.
_REGISTRIES: dict[str, MCPRegistry] = {}


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


def registry_key(workspace: Path | None = None) -> str:
    """The identity a registry is cached under: the resolved workspace.

    Resolved, so ``/repo``, ``/repo/`` and a symlink to it are one entry
    rather than three sets of server subprocesses."""
    if not workspace:
        return ""
    try:
        return str(Path(workspace).expanduser().resolve())
    except (OSError, ValueError):
        return str(workspace)


def get_registry(workspace: Path | None = None) -> MCPRegistry:
    """Return the MCP registry for ``workspace``, lazily creating it."""
    key = registry_key(workspace)
    with _REGISTRY_LOCK:
        registry = _REGISTRIES.get(key)
        if registry is None:
            registry = _REGISTRIES[key] = MCPRegistry()
        if not registry.loaded:
            registry.load(workspace)
    return registry


def reset_registry(workspace: Path | None = None) -> None:
    """Stop and clear registries — one workspace's, or all of them.

    Used by ``/mcp reload`` and by tests. Shutdown happens after the lock
    is released: it terminates subprocesses and waits on them, which is
    not work to hold a lock across."""
    with _REGISTRY_LOCK:
        if workspace is None:
            stopping = list(_REGISTRIES.values())
            _REGISTRIES.clear()
        else:
            found = _REGISTRIES.pop(registry_key(workspace), None)
            stopping = [found] if found is not None else []
    for registry in stopping:
        try:
            registry.shutdown()
        except Exception:
            pass


def effective_servers(workspace: Path | None) -> list[dict]:
    """Every server that WOULD load, with the file it came from.

    The ``/mcp`` listing read the user-global config directly, so it
    could not show a builtin, could not show a workspace's entry, and
    could not say which file any line came from. This is what actually
    loads, named by source.
    """
    configs, sources, _notice = _load_configs_with_sources(workspace)
    out: list[dict] = []
    for name, cfg in configs.items():
        out.append({
            "name": name,
            "command": cfg.get("command", ""),
            "args": list(cfg.get("args") or []),
            "env": dict(cfg.get("env") or {}),
            "url": cfg.get("url", ""),
            "enabled": bool(cfg.get("enabled", True)),
            "source": sources.get(name, ""),
        })
    out.sort(key=lambda r: r["name"])
    return out


def trust_notice(workspace: Path | None) -> str:
    """What a lack of trust withheld for *workspace*, or ""."""
    return _load_configs_with_sources(workspace)[2]


__all__ = [
    "MCPTool",
    "MCPResource",
    "MCPPrompt",
    "MCPServer",
    "MCPRegistry",
    "effective_servers",
    "get_registry",
    "registry_key",
    "reset_registry",
    "trust_notice",
]
