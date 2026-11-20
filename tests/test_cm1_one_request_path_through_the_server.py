"""CM1: the request-handling behaviour lives in one place, the server's.

``code_mode.serve_request`` was a one-line wrapper:

    return Connection(executor, permissions).handle(request)

``RpcServer._handle_payload_in_context`` makes the same call, plus the
``error``/``result`` split and the captured contextvar context. Two spellings
of one behaviour means a change to the gate has to be made twice, and the
second one has no transport to prove it on. The wrapper is gone; these tests
hold the capability to its real entry point, ``RpcServer.handle_payload``,
the method the product's serve loop calls through ``serve_one``.

Failure shapes are asserted against the SERVER's protocol, which differs from
what the wrapper returned: the server answers ``{"error": ...}`` or
``{"result": ...}``, so a refusal read from the server is one level up.
"""
import json

from delfin.agent import code_mode


class _FakeExecutor:
    """Stands in for KitToolExecutor.execute, the public entry the RPC calls."""

    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions=None):
        self.calls.append((name, arguments))
        if name == "get_calc_info":
            return json.dumps({"calc_id": arguments.get("calc_id", ""),
                               "energy": -1.0})
        return json.dumps({"error": "Unknown tool: %s" % name})


def _payload(request):
    return json.dumps(request)


def test_server_returns_the_executor_result_at_the_result_key():
    ex = _FakeExecutor()
    server = code_mode.RpcServer(ex, permissions=None)
    reply = json.loads(server.handle_payload(
        _payload({"tool": "get_calc_info", "arguments": {"calc_id": "c1"}})))
    assert ex.calls == [("get_calc_info", {"calc_id": "c1"})]
    assert reply["result"]["energy"] == -1.0
    assert "error" not in reply


def test_server_refuses_a_write_tool_without_touching_the_executor():
    ex = _FakeExecutor()
    server = code_mode.RpcServer(ex, permissions=None)
    reply = json.loads(server.handle_payload(_payload(
        {"tool": "write_file", "arguments": {"path": "x", "content": "y"}})))
    assert ex.calls == []
    assert "write_file" in reply["error"]


def test_server_refuses_a_submit_and_delete_tool_without_the_executor():
    ex = _FakeExecutor()
    server = code_mode.RpcServer(ex, permissions=None)
    for name in ("submit_calculation", "delete_calc"):
        reply = json.loads(server.handle_payload(
            _payload({"tool": name, "arguments": {}})))
        assert "error" in reply, name
        assert ex.calls == [], name


def test_server_surfaces_a_gate_refusal_as_an_error():
    ex = _FakeExecutor()
    server = code_mode.RpcServer(ex, permissions=None)
    reply = json.loads(server.handle_payload(
        _payload({"tool": "nonexistent_tool", "arguments": {}})))
    assert "nonexistent_tool" in reply["error"]


def test_server_caps_an_oversized_executor_result():
    class _BigExecutor:
        def execute(self, name, arguments, permissions=None):
            return json.dumps(
                {"blob": "x" * (code_mode.MAX_RESULT_BYTES + 1000)})

    server = code_mode.RpcServer(_BigExecutor(), permissions=None)
    reply = json.loads(server.handle_payload(
        _payload({"tool": "get_calc_info", "arguments": {}})))
    assert "truncated" in reply["error"]
    assert len(json.dumps(reply)) < code_mode.MAX_RESULT_BYTES
    assert "xxx" not in json.dumps(reply)


def test_server_refuses_a_malformed_payload():
    server = code_mode.RpcServer(_FakeExecutor(), permissions=None)
    assert "error" in json.loads(server.handle_payload("{not json"))
    assert "error" in json.loads(server.handle_payload("[1, 2, 3]"))


def test_the_duplicate_spelling_is_gone():
    """One request path, not two.

    ``serve_request`` returned ``Connection.handle``'s raw reply while the
    server wrapped it as ``{"error": ...}`` / ``{"result": ...}``; a caller
    that moved between them silently changed protocol. Red on the version
    that had both: the name existed there.
    """
    assert not hasattr(code_mode, "serve_request")
    assert hasattr(code_mode.RpcServer, "handle_payload")
    assert hasattr(code_mode.Connection, "handle")


def test_two_requests_on_one_connection_keep_the_per_call_allow_list():
    """A second, different call is judged on its own, never carried by the
    first one's success."""
    ex = _FakeExecutor()
    conn = code_mode.Connection(ex, permissions=None)
    ok = json.loads(conn.handle({"tool": "get_calc_info",
                                 "arguments": {"calc_id": "c"}}))
    assert "energy" in ok
    bad = json.loads(conn.handle({"tool": "bash",
                                  "arguments": {"command": "id"}}))
    assert "error" in bad
    assert len(ex.calls) == 1
