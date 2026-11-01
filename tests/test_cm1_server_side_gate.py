"""CM1 QS red tests (s11 comments): the server side is the real boundary.

These drive the SERVER wrapper (code_mode.RpcServer) that the protected
``_dispatch_under_rules`` hook will call. The client-side refusal is
cosmetic — the program is untrusted code; the server-side allow-list before
``_dispatch_under_rules`` is the boundary that must hold.
"""
import json
import os

from delfin.agent import code_mode


class _RecordingExecutor:
    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions=None):
        self.calls.append(name)
        return json.dumps({"ok": name})


def test_server_refuses_fabricated_write_call_even_if_client_allows():
    """The program may fabricate a direct socket payload with a WRITE tool —
    the server refuses it before the executor is touched, regardless of the
    client stub."""
    ex = _RecordingExecutor()
    server = code_mode.RpcServer(ex, permissions=None)
    payload = json.dumps({"tool": "write_file",
                          "arguments": {"path": "x", "content": "y"}})
    reply = json.loads(server.handle_payload(payload))
    assert "error" in reply
    assert "write_file" in reply["error"]
    assert ex.calls == []


def test_server_refuses_malformed_payload():
    ex = _RecordingExecutor()
    server = code_mode.RpcServer(ex, permissions=None)
    assert "error" in json.loads(server.handle_payload("{not json"))
    assert "error" in json.loads(server.handle_payload("[1,2,3]"))
    assert ex.calls == []


def test_server_passes_read_only_call_through():
    ex = _RecordingExecutor()
    server = code_mode.RpcServer(ex, permissions=None)
    reply = json.loads(server.handle_payload(
        json.dumps({"tool": "search_calcs", "arguments": {"query": "PBE0"}})))
    assert reply == {"result": {"ok": "search_calcs"}}
    assert ex.calls == ["search_calcs"]


def test_socket_path_is_unique_per_invocation(tmp_path):
    """Two concurrent caged programs cannot cross-talk: the socket path is
    unique per invocation (token in the name)."""
    a = code_mode.rpc_socket_path(directory=str(tmp_path), token="aaa")
    b = code_mode.rpc_socket_path(directory=str(tmp_path), token="bbb")
    assert a != b
    assert str(tmp_path) in a  # under the workspace, not /tmp shared state


def test_timeout_leaves_no_stale_socket(tmp_path):
    """After a killed run the next invocation must work: the socket file the
    timeout left behind is removed before binding a new one."""
    stale = os.path.join(str(tmp_path), "stale.sock")
    open(stale, "wb").close()
    assert code_mode.rpc_socket_path(directory=str(tmp_path),
                                     token="s").endswith(".sock")
    server = code_mode.RpcServer(_RecordingExecutor(), permissions=None)
    server.bind(stale)  # must not raise on the stale file
    assert os.path.exists(stale)
