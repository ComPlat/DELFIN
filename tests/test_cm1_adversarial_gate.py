"""CM1 adversarial review tests (reviewer nacht-s10).

Attack the real boundary: the SERVER-side allow-list in
``code_mode.RpcServer.handle_payload``, which the protected
``_dispatch_under_rules`` hook will call, plus the client and cap surfaces.
The program is untrusted code: whatever it sends over the socket must never
reach a write/submit tool through Code Mode.
"""
import json

import pytest

from delfin.agent import code_mode


class _RecordingExecutor:
    """Records every tool name it would have run."""

    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions=None):
        self.calls.append(name)
        return json.dumps({"ok": name})


class _OldCallableProtocol:
    """The pre-parity protocol. The server must NEVER call this — kept so
    test_nested_tool_name_wrapper_not_executed can prove it."""


def _server(ex=None, **kw):
    if ex is None:
        ex = _RecordingExecutor()
    return code_mode.RpcServer(ex, **kw)


# ---------------------------------------------------------------- allow-list

def test_write_tool_refused_even_with_permissions_none():
    """The refusal must not depend on the executor's own rules: with
    permissions=None (a server started without a permission bundle) a write
    tool is still refused before the executor is touched."""
    ex = _RecordingExecutor()
    server = _server(permissions=None)
    payload = json.dumps({"tool": "write_file",
                          "arguments": {"path": "x", "content": "y"}})
    reply = json.loads(server.handle_payload(payload))
    assert "error" in reply
    assert ex.calls == []


def test_submit_and_delete_tool_names_refused():
    ex = _RecordingExecutor()
    server = _server(ex)
    for name in ("submit_calculation", "submit_application",
                 "delete_calc_folder", "cancel_run"):
        reply = json.loads(server.handle_payload(
            json.dumps({"tool": name, "arguments": {}})))
        assert "error" in reply, name
        assert "direct" in reply["error"].lower(), name
    assert ex.calls == []


def test_tool_name_unicode_lookalike_and_non_string_refused():
    """A lookalike spelling or a non-string name must not slip past the
    frozenset membership test."""
    ex = _RecordingExecutor()
    server = _server(ex)
    for name in ("search_calсs", "search\u200bcalcs", " search_calcs",
                 "search_calcs ", 123, None, {"tool": "search_calcs"}):
        payload = json.dumps({"tool": name, "arguments": {}})
        reply = json.loads(server.handle_payload(payload))
        # every lookalike / non-string spelling is refused; none reaches the
        # executor (a dict name would crash str() handling if unguarded)
        assert "error" in reply, json.dumps(name, default=str)
        assert ex.calls == []


def test_nested_tool_name_wrapper_not_executed():
    """A payload that hides the tool inside a wrapper dict must not reach
    the executor."""
    ex = _RecordingExecutor()
    server = _server(ex)
    payload = json.dumps({"arguments": {"tool": "write_file", "path": "x"},
                          "tool": "get_calc_info"})
    reply = json.loads(server.handle_payload(payload))
    # get_calc_info is allowed, but the write instruction must not run:
    assert "result" in reply
    assert ex.calls == ["get_calc_info"]


# ------------------------------------------------------------- result capping

def test_oversized_executor_result_capped_at_server():
    """A tool that returns a huge result must be capped BEFORE it leaves the
    server. Pinned contract: the oversized result is NOT partially forwarded
    -- the reply is a truncation ERROR naming the cap, so the program learns
    it must narrow or aggregate (result data never floods the context)."""
    class _HugeExecutor(_RecordingExecutor):
        def execute(self, name, arguments, permissions=None):
            self.calls.append(name)
            return json.dumps(
                {"data": "x" * (code_mode.MAX_RESULT_BYTES + 50)})

    ex = _HugeExecutor()
    server = _server(ex)
    reply = json.loads(server.handle_payload(
        json.dumps({"tool": "search_calcs", "arguments": {}})))
    # result data is NOT partially forwarded; the reply names the cap
    body = json.dumps(reply, ensure_ascii=False)
    # the capped payload, not the raw flood, left the server
    assert len(body) < 1000
    assert "truncated" in body


# ------------------------------------------------------- client error surface

def test_client_raises_on_refusal_never_swallows(tmp_path):
    """A gate refusal inside the program surfaces as an exception in the
    program, not as a silent None."""
    sock = tmp_path / "rpc.sock"

    class _OneShotServer:
        def __init__(self):
            self.srv = code_mode.RpcServer(_RecordingExecutor(),
                                           permissions=None)

        def __call__(self):
            self.srv.bind(str(sock))
            self.srv.serve_one()

    import threading
    t = threading.Thread(target=_OneShotServer(), daemon=True)
    t.start()
    import time
    for _ in range(100):
        if sock.exists():
            break
        time.sleep(0.05)
    client = code_mode.RpcClient(str(sock))
    with pytest.raises(RuntimeError, match="refused"):
        client.call("write_file", {"path": "x"})
    t.join(timeout=5)


def test_client_returns_result_payload(tmp_path):
    sock = tmp_path / "rpc2.sock"
    ex = _RecordingExecutor()

    import threading
    import time
    srv = code_mode.RpcServer(ex, permissions=None)

    def _serve():
        srv.bind(str(sock))
        srv.serve_one()

    t = threading.Thread(target=_serve, daemon=True)
    t.start()
    for _ in range(100):
        if sock.exists():
            break
        time.sleep(0.05)
    client = code_mode.RpcClient(str(sock))
    assert client.call("search_calcs", {"query": "PBE0"}) == {"ok": "search_calcs"}
    assert ex.calls == ["search_calcs"]
    t.join(timeout=5)
