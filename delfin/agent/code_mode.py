"""Code Mode (CM1): DELFIN's tools as a typed Python API for the model.

Instead of one model turn per tool call, the model writes ONE short program
against a typed API of DELFIN's read-only tools. The program runs in the
existing bash cage (``api_client._bash_isolation_argv``) and reaches the
tools only over a narrow local RPC that lands in the SAME central gate every
direct call takes (``KitToolExecutor._dispatch_under_rules``, api_client.py:18282).

Phase 1 = READ-ONLY tools only. A call to any write/submit/delete tool is
refused with a message naming the direct route. The allow-list is enforced
per call, twice: once in :class:`Connection` (this module) and again inside
the agent-side RPC server wrapper before the executor is touched.

The model runs a program like::

    python -m delfin.agent.code_mode run my_program.py

inside the cage. ``DELFIN_CODE_MODE_RPC`` (the socket path) is exported ONLY
to that caged child process; without it :func:`make_client` returns ``None``
and the stub functions raise — no program can fake a connection.
"""
from __future__ import annotations

import json
import os
import subprocess
import contextvars
from shlex import quote as _shlex_quote
import sys
from typing import Any, Callable, Iterable, Optional

#: Tools the model may call from a Code Mode program. Phase 1: read-only,
#: job status, output reading/parsing, listing, tables. Every name must exist
#: in ``_DOC_TOOLS_OPENAI`` (the real catalogue) — never an invented tool.
READ_ONLY_TOOLS: frozenset[str] = frozenset({
    "search_calcs", "get_calc_info", "search_docs", "read_section",
    "bash_status", "bash_output", "list_docs",
})

#: The tools a program must NOT call, and how to do it properly instead.
_DIRECT_ROUTE_HINT = (
    "Code Mode is read-only (phase 1). To use tool '%s' directly, the model "
    "calls it as a normal individual tool call — permissions, write scope "
    "and approvals apply there exactly as before."
)

#: Results are capped before they return to the model: a program that
#: tabulates a whole run tree must not flood the context.
MAX_RESULT_BYTES = 200_000

#: Default wall-clock limit for one program run.
DEFAULT_TIMEOUT_S = 120.0


class ProgramTimeout(RuntimeError):
    """The program ran longer than its wall-clock limit."""


class ProgramError(RuntimeError):
    """The program failed (non-zero exit); stderr is part of the message."""


def _check_allowed(name: str) -> Optional[str]:
    """Refusal text for a tool outside the read-only allow-list, else None."""
    if name in READ_ONLY_TOOLS:
        return None
    return _DIRECT_ROUTE_HINT % name


def _cap_result(result: str) -> str:
    if len(result) <= MAX_RESULT_BYTES:
        return result
    note = json.dumps({
        "error": (
            "Code Mode result truncated: %d bytes over the %d-byte cap "
            "(MAX_RESULT_BYTES). Narrow the query (e.g. max_results) or "
            "aggregate in the program before returning." % (
                len(result) - MAX_RESULT_BYTES, MAX_RESULT_BYTES)),
    }, ensure_ascii=False)
    return note


class Connection:
    """One program connection: enforces the allow-list per call."""

    def __init__(self, executor: Any, permissions: Any = None):
        self._executor = executor
        self._permissions = permissions

    def handle(self, request: dict) -> str:
        """One tool request from the program -> executor result or refusal."""
        name = str(request.get("tool", "") or "")
        arguments = request.get("arguments")
        if not isinstance(arguments, dict):
            arguments = {}
        refusal = _check_allowed(name)
        if refusal is not None:
            return json.dumps({"error": refusal}, ensure_ascii=False)
        # Public-entry parity (CM1 operator finding): the direct path enters
        # at KitToolExecutor.execute (api_client.py:17651), which carries
        # deny rules, plan-mode, sandbox, role and shape checks BEFORE
        # _dispatch_under_rules. The RPC must call the same single public
        # entry and nothing lower; rules are evaluated inside execute.
        return _cap_result(self._executor.execute(
            name, arguments, self._permissions))


def serve_request(request: dict, executor: Any,
                  permissions: Any = None) -> str:
    """One-shot request handling (no connection state kept)."""
    return Connection(executor, permissions).handle(request)


def runner_program_path(cmd: str) -> Optional[str]:
    """The program path if ``cmd`` is a Code Mode runner invocation, else None.

    Recognises ``... delfin.agent.code_mode run <program>`` in the bash
    command text. The env var this gates grants nothing by itself: the
    RPC server's own read-only allow-list is the boundary, so even a
    command that names the runner to pick up the socket can only reach
    the read-only tools through the gate.
    """
    import re
    match = re.search(
        r"delfin\.agent\.code_mode\s+run\s+([^\s;&|<>$`()\']+)", cmd or "")
    if match:
        return match.group(1)
    return None


def rpc_socket_path(directory: str, token: str) -> str:
    """Unique-per-invocation socket path under the given directory.

    The token (uuid4 hex) keeps two concurrent caged programs from
    cross-talking; under the directory the command runs in (the caged
    child's own write scope, never shared /tmp state).
    """
    import uuid
    base = os.path.join(directory or os.getcwd(), ".code_mode_rpc")
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, "rpc-%s.sock" % (token or uuid.uuid4().hex))


class RpcServer:
    """Agent-side RPC server wrapper — the boundary of Code Mode.

    The protected hook inside ``_dispatch_under_rules`` (api_client.py:18282,
    delivered as .gate/px copies) starts one of these per Code Mode bash
    command and hands each request payload here. The allow-list is checked
    HERE, before the executor — a fabricated payload naming a write tool is
    refused without touching the executor, whatever the client stub claims.
    """

    def __init__(self, executor: Any, permissions: Any = None):
        self._executor = executor
        self._permissions = permissions
        self._socket_path: Optional[str] = None
        # The hook serves requests from a daemon thread, but DELFIN's
        # executor reads contextvars (_RULE_CALL, ...) that a plain
        # threading.Thread does NOT inherit. Capture the constructing
        # (agent) thread's context and run every request under it.
        self._context = contextvars.copy_context()

    def handle_payload(self, payload: str) -> str:
        """One request payload -> reply string."""
        return self._context.run(self._handle_payload_in_context, payload)

    def _handle_payload_in_context(self, payload: str) -> str:
        try:
            request = json.loads(payload)
        except (ValueError, TypeError):
            return json.dumps({"error": "Code Mode: malformed request payload"})
        if not isinstance(request, dict):
            return json.dumps({"error": "Code Mode: request must be an object"})
        result = Connection(self._executor, self._permissions).handle(request)
        # Connection.handle already returns capped JSON; wrap as
        # {"result": ...} for the client's error/success split.
        try:
            data = json.loads(result)
        except ValueError:
            return json.dumps({"result": result})
        if isinstance(data, dict) and "error" in data:
            return json.dumps({"error": data["error"]}, ensure_ascii=False)
        return json.dumps({"result": data}, ensure_ascii=False)

    def bind(self, socket_path: str) -> str:
        """Bind the listening socket, clearing a stale file first.

        A timeout-killed previous run may leave a socket file; removing a
        file at THIS path before rebinding is safe (unique per invocation)
        and keeps the next invocation working.
        """
        try:
            os.unlink(socket_path)
        except FileNotFoundError:
            pass
        import socket as _socket
        self._socket_path = socket_path
        self._sock = _socket.socket(_socket.AF_UNIX, _socket.SOCK_STREAM)
        self._sock.bind(socket_path)
        self._sock.listen(8)
        return socket_path

    def serve_one(self, idle_timeout_s: Optional[float] = None) -> None:
        """Serve ONE program connection, then keep listening.

        The program's client opens a fresh socket per call, so the listener
        stays open. With ``idle_timeout_s`` a connection that never arrives
        within the window ends the server (the program died before its
        first call) — this is how the agent-side serve loop self-terminates
        instead of leaking a listener per invocation.
        """
        if idle_timeout_s is not None:
            self._sock.settimeout(idle_timeout_s)
        conn, _ = self._sock.accept()
        try:
            with conn:
                buf = conn.recv(1_000_000)
                conn.sendall(self.handle_payload(
                    buf.decode("utf-8", errors="replace")).encode("utf-8"))
        finally:
            # After a timeout accept() the socket timeout must be cleared for
            # the next serve_one call. OSError here only happens when close()
            # already ran, in which case there is nothing left to clear —
            # this except is a no-op guard, not an error swallow.
            try:
                self._sock.settimeout(None)
            except OSError:
                pass

    def close(self) -> None:
        """Close the listener and unlink the socket file (idempotent)."""
        try:
            self._sock.close()
        except OSError:
            pass
        path = self._socket_path
        if path:
            self._socket_path = None
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass


def make_client(rpc_path: str) -> Optional["RpcClient"]:
    """Client for the agent-side RPC server, or None without a socket path.

    The socket path reaches the program ONLY through the ``DELFIN_CODE_MODE_RPC``
    env var, which the runner sets for the caged child alone. An empty path
    (no env var) means no Code Mode connection — the stub raises instead of
    pretending to work.
    """
    if not rpc_path:
        return None
    return RpcClient(rpc_path)


class RpcClient:
    """Program-side client: one request per call over the local socket.

    Deliberately minimal — no retries, no fallback route. If the RPC is
    refused (unknown tool, gate refusal, transport error) the error is
    raised to the program, never swallowed.
    """

    def __init__(self, rpc_path: str):
        self.rpc_path = rpc_path

    def call(self, tool: str, arguments: Optional[dict] = None) -> Any:
        import socket
        payload = json.dumps({"tool": tool, "arguments": arguments or {}})
        # unix-family socket only; no network, no host resolution.
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as s:
            s.settimeout(60.0)
            s.connect(self.rpc_path)
            s.sendall(payload.encode("utf-8"))
            chunks = []
            while True:
                chunk = s.recv(65536)
                if not chunk:
                    break
                chunks.append(chunk)
        result = b"".join(chunks).decode("utf-8", errors="replace")
        data = json.loads(result)
        if "error" in data:
            raise RuntimeError("Code Mode tool call refused: %s"
                               % data["error"])
        return data.get("result", data)


def generate_stub() -> str:
    """The typed Python API the model programs against.

    Generated from the read-only allow-list so the model sees only tools
    Code Mode may execute. Every function goes over the RPC; nothing here
    bypasses the gate.
    """
    lines = [
        '"""Code Mode API — DELFIN\'s read-only tools as typed functions.',
        "",
        "Import this in your program (it is added to sys.path by the "
        "runner): from code_mode_api import delfin_tools",
        "",
        "Every call goes through the agent's central gate. Refusals raise.",
        '"""',
        "import json",
        "import os",
        "import sys",
        "",
        "_RPC_PATH = os.environ.get('DELFIN_CODE_MODE_RPC', '')",
        "if _RPC_PATH:",
        "    from delfin.agent.code_mode import make_client",
        "    _client = make_client(_RPC_PATH)",
        "else:",
        "    _client = None",
        "",
        "def _call(tool, arguments):",
        "    if _client is None:",
        "        raise RuntimeError(",
        "            'Code Mode is not connected: DELFIN_CODE_MODE_RPC is not "
        "set. Run programs with: python -m delfin.agent.code_mode run "
        "<program>')",
        "    return _client.call(tool, arguments)",
        "",
    ]
    for name in sorted(READ_ONLY_TOOLS):
        lines.append("def %s(**arguments):" % name)
        lines.append("    return _call(%r, arguments)" % name)
        lines.append("")
    lines.append("class delfin_tools:")
    for name in sorted(READ_ONLY_TOOLS):
        lines.append("    %s = staticmethod(%s)" % (name, name))
    lines.append("")
    lines.append("if __name__ == '__main__':")
    lines.append("    print(sorted(n for n in dir(delfin_tools)"
                 " if not n.startswith('_')))")
    return "\n".join(lines) + "\n"


def run_program(program_path: str, timeout_s: float = DEFAULT_TIMEOUT_S,
                rpc_path: str = "", workspace: str = "") -> str:
    """Run the model's program and return its stdout (capped).

    In the product the command runs through the bash tool, hence inside
    ``_bash_isolation_argv``'s cage. Here (agent side) we run it with the
    same wall-clock limit and cap so the process boundaries stay identical.
    The socket path is exported only to this child.
    """
    argv = ["/bin/bash", "-c",
            "%s -m delfin.agent.code_mode run %s"
            % (sys.executable, _shlex_quote(program_path))]
    env = dict(os.environ)
    if rpc_path:
        env["DELFIN_CODE_MODE_RPC"] = rpc_path
    else:
        env.pop("DELFIN_CODE_MODE_RPC", None)
    if workspace:
        env["PYTHONPATH"] = (workspace + os.pathsep + env.get("PYTHONPATH", "")
                             ).rstrip(os.pathsep)
    try:
        proc = subprocess.run(
            argv, capture_output=True, timeout=timeout_s, env=env,
            cwd=workspace or None)
    except subprocess.TimeoutExpired as exc:
        raise ProgramTimeout(
            "Code Mode program exceeded the %.0fs limit; batch smaller "
            "(fewer tool calls per run) or pass a higher timeout_s next run"
            % timeout_s
        ) from exc
    out = proc.stdout.decode("utf-8", errors="replace")
    err = proc.stderr.decode("utf-8", errors="replace")
    if proc.returncode != 0:
        raise ProgramError("exit %d: %s" % (proc.returncode, err.strip()))
    return _cap_result(out)


def main(argv: Optional[list[str]] = None) -> int:
    """``python -m delfin.agent.code_mode run <program>`` — the program side."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2 or argv[0] != "run":
        print("usage: python -m delfin.agent.code_mode run <program.py>",
              file=sys.stderr)
        return 2
    program = argv[1]
    sys.path.insert(0, os.path.dirname(os.path.abspath(program)))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    stub_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    # The stub is materialised next to the program by the caller (the bash
    # command writes it with the generated_stub() text); import by name.
    import runpy
    try:
        runpy.run_path(program, run_name="__main__")
    except SystemExit as exc:
        return int(exc.code or 0)
    return 0


if __name__ == "__main__":  # pragma: no cover
    import sys
    raise SystemExit(main())
