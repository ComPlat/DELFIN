"""The server thread must see the agent's contextvars snapshot.

RpcServer.handle_payload runs in the hook's daemon thread, but DELFIN's
executor reads contextvars (_RULE_CALL, ...) that a plain threading.Thread
does NOT inherit. The server therefore captures the constructing thread's
context (the agent's, at server construction) and runs every request under
it. Red on 49a6bcaea: handle_payload ran in whatever context the serving
thread had (empty contextvars), so downstream contextvar reads saw nothing.
"""
import contextvars
import threading

from delfin.agent import code_mode

_MARKER = contextvars.ContextVar("cm1_marker", default=None)


def _serve_once(server, payload):
    out = []
    t = threading.Thread(target=lambda: out.append(
        server.handle_payload(payload)))
    t.start()
    t.join(timeout=10)
    return out[0]


def test_handle_payload_runs_under_the_constructing_context(tmp_path):
    """A contextvar set in the agent thread is visible to a request handled
    in a DIFFERENT thread, because the server captured the context."""
    _MARKER.set("agent-snapshot")
    seen = []

    executor = code_mode.READ_ONLY_STUB if hasattr(
        code_mode, "READ_ONLY_STUB") else None

    class _Recording:
        def execute(self, name, arguments, permissions):
            seen.append(_MARKER.get())
            return "{}"

    server = code_mode.RpcServer(_Recording(), permissions="perms-obj")
    reply = _serve_once(server, '{"tool": "get_calc_info", "arguments": {}}')
    assert "error" not in reply
    assert seen == ["agent-snapshot"], "agent's contextvar not visible"


def test_without_the_captured_context_the_marker_is_lost():
    """The control: a plain thread does NOT inherit contextvars — this is
    the defect the capture fixes (documents the mechanism, not DELFIN)."""
    _MARKER.set("agent-snapshot")
    seen = []
    t = threading.Thread(target=lambda: seen.append(_MARKER.get()))
    t.start()
    t.join(timeout=10)
    assert seen == [None]
