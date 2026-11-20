"""Code Mode (CM1): DELFIN's tools as a typed Python API for the model.

Instead of one model turn per tool call, the model writes ONE short program
against a typed API of DELFIN's read-only tools. The program runs in the
existing bash cage (``api_client._bash_isolation_argv``) and reaches the
tools only over a narrow local RPC that lands in the SAME central gate every
direct call takes (``KitToolExecutor.execute``).

There is exactly ONE launcher, and it is not in this module: the bash tool.
``api_client._code_mode_child_env`` recognises a runner command, starts the
:class:`RpcServer` and exports the socket to the caged child. A second
launcher here would build its own argv -- without the cage, the wall-clock
ceiling or the output truncation the bash tool owns -- and would be a way
round them.

Phase 1 = READ-ONLY tools only. A call to any write/submit/delete tool is
refused with a message naming the direct route. The allow-list is enforced
per call, twice: once in :class:`Connection` (this module) and again inside
the agent-side RPC server wrapper before the executor is touched.

The model obtains the API module and runs a program like::

    python -m delfin.agent.code_mode stub > tools_api.py
    python -m delfin.agent.code_mode run my_program.py

inside the cage. ``DELFIN_CODE_MODE_RPC`` (the socket path) is exported ONLY
to that caged child process; without it :func:`make_client` returns ``None``
and the stub functions raise — no program can fake a connection.
"""
from __future__ import annotations

import json
import os
import contextvars
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
        "Save this file as tools_api.py next to your program, then: "
        "import tools_api and call the functions below.",
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


def main(argv: Optional[list[str]] = None) -> int:
    """``python -m delfin.agent.code_mode <stub|run> ...``.

    Two subcommands, both documented in the Code Mode skill:

    * ``stub`` writes the typed API (:func:`generate_stub`) to stdout, so
      ``... stub > tools_api.py`` materialises the module a program imports.
      Without this branch the documented command printed the ``run`` usage and
      exited 2, so the API module could not be obtained at all.
    * ``run <program.py>`` executes the program: the program side of the RPC.

    There is no third subcommand and no in-module launcher. A program is
    started by the bash tool, where ``api_client._code_mode_child_env``
    recognises this command, starts the RPC server and exports the socket to
    the caged child; the bash tool owns the cage, the wall-clock ceiling and
    the output truncation for it.

    Returns 0 on success, 2 for a usage error.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv == ["stub"]:
        sys.stdout.write(generate_stub())
        return 0
    if len(argv) != 2 or argv[0] != "run":
        print("usage: python -m delfin.agent.code_mode stub | "
              "run <program.py>", file=sys.stderr)
        return 2
    program = argv[1]
    sys.path.insert(0, os.path.dirname(os.path.abspath(program)))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    # The stub is materialised next to the program by the caller; import by name.
    import runpy
    try:
        runpy.run_path(program, run_name="__main__")
    except SystemExit as exc:
        return int(exc.code or 0)
    return 0


if __name__ == "__main__":  # pragma: no cover
    import sys
    raise SystemExit(main())
