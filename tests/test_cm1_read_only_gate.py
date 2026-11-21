"""CM1 Code Mode: the program's tool calls go back through the same gate.

Phase 1 is READ-ONLY tools only. Every rule tested here:
- a program call to a write/submit tool is refused, naming the direct route;
- a read-only call is executed through the PUBLIC executor entry (the same
  execute() a direct call takes; deny rules, plan-mode, sandbox live there)
  and its result is returned to the program;
- a gate refusal inside the program surfaces as an error to the program,
  not silently;
- the allow-list is enforced PER CALL: one connection may not smuggle a
  second, different call through after a first allowed one;
- results are size-capped before they reach the program.

Every case enters at ``RpcServer.handle_payload``, the request path the
product's serve loop (``serve_one``) calls. There is no second request entry
and no in-module program launcher: a program is started by the bash tool,
which owns the cage, the wall-clock ceiling and the output truncation.

The wall-clock limit and the stdout cap are therefore not tested here: the
launcher that owned them locally (``run_program``) is gone, and
``test_cm1_the_bash_route_runs_a_program.py`` holds the bash route to the
ceiling it inherits. Kept from the branch this file was merged with, because
a reader of this file needs to know where that coverage went.
"""
import json

from delfin.agent import code_mode


def _catalogue():
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    return {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}


def test_read_only_tools_exist_in_the_real_catalogue():
    """The allow-list names real, read-only tools — never an invented one."""
    assert code_mode.READ_ONLY_TOOLS
    assert code_mode.READ_ONLY_TOOLS <= _catalogue()


def test_read_only_tools_disjoint_from_write_tools():
    from delfin.agent.api_client import _WRITE_TOOL_NAMES
    assert not (code_mode.READ_ONLY_TOOLS & _WRITE_TOOL_NAMES)


class _FakeExecutor:
    """Stands in for KitToolExecutor._dispatch_under_rules: records every
    call and answers like it would."""

    def __init__(self):
        self.calls = []

    # CM1 public-entry parity: the server calls the PUBLIC entry
    # (executor.execute(name, arguments, permissions)), the same one a
    # direct call takes — rules, plan-mode, sandbox are evaluated inside.
    def execute(self, name, arguments, permissions=None):
        self.calls.append((name, arguments))
        if name == "get_calc_info":
            return json.dumps({"calc_id": arguments.get("calc_id", ""),
                               "energy": -1.0})
        return json.dumps({"error": "Unknown tool: %s" % name})


def _server(ex, permissions=None):
    """The one request path: what the product's serve loop calls."""
    return code_mode.RpcServer(ex, permissions=permissions)


def _payload(server, request):
    """Send one request through the server and decode its reply."""
    return json.loads(server.handle_payload(json.dumps(request)))


def test_read_only_call_reaches_the_executor_and_returns():
    ex = _FakeExecutor()
    reply = _payload(_server(ex), {"tool": "get_calc_info",
                                   "arguments": {"calc_id": "c1"}})
    assert ex.calls == [("get_calc_info", {"calc_id": "c1"})]
    assert reply["result"]["energy"] == -1.0


def test_write_tool_refused_and_names_the_direct_route():
    ex = _FakeExecutor()
    reply = _payload(_server(ex), {"tool": "write_file",
                                   "arguments": {"path": "x", "content": "y"}})
    assert "error" in reply
    assert ex.calls == []  # never dispatched
    assert "write_file" in reply["error"]


def test_submit_tool_refused_by_name():
    ex = _FakeExecutor()
    reply = _payload(_server(ex), {"tool": "submit_calculation",
                                   "arguments": {"folder": "f"}})
    assert "error" in reply
    assert ex.calls == []


def test_gate_refusal_surfaces_as_error_not_silence():
    ex = _FakeExecutor()
    reply = _payload(_server(ex), {"tool": "nonexistent_tool",
                                   "arguments": {}})
    assert "error" in reply


def test_allow_list_is_enforced_per_call():
    """A later disallowed call on the SAME connection is refused too."""
    ex = _FakeExecutor()
    conn = code_mode.Connection(ex, permissions=None)
    ok = json.loads(conn.handle({"tool": "get_calc_info",
                                 "arguments": {"calc_id": "c"}}))
    assert "energy" in ok
    bad = json.loads(conn.handle({"tool": "bash",
                                  "arguments": {"command": "id"}}))
    assert "error" in bad
    assert len(ex.calls) == 1


def test_result_is_size_capped(tmp_path):
    """An oversized result is replaced by a cap note, never the raw blob."""
    class _BigExecutor:
        def execute(self, name, arguments, permissions=None):
            return json.dumps(
                {"blob": "x" * (code_mode.MAX_RESULT_BYTES + 1000)})

    reply = _payload(_server(_BigExecutor()), {"tool": "get_calc_info",
                                               "arguments": {}})
    assert "error" in reply
    assert "truncated" in reply["error"]
    assert len(json.dumps(reply)) < code_mode.MAX_RESULT_BYTES
    assert "xxx" not in json.dumps(reply)


def test_program_without_rpc_env_gets_no_client():
    """Without the socket env var the stub must not pretend a connection."""
    client = code_mode.make_client(rpc_path="")
    assert client is None


def test_generated_stub_contains_only_read_only_tools():
    stub = code_mode.generate_stub()
    for name in code_mode.READ_ONLY_TOOLS:
        assert 'def %s(' % name in stub
    assert "write_file" not in stub
    assert "submit_calculation" not in stub
