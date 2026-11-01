"""CM1 reviewer tests for operator FINDING E: the RPC must enter at the
SINGLE public entry (executor.execute), not at _dispatch_under_rules.

The direct path enters KitToolExecutor.execute -> _execute_pinned
(api_client.py:17651/17672), which carries: managed/user DENY rules, plan-
mode deny-by-default, the no-sandbox refusal, the foreground-server and
apply_patch shape checks, the pre-tool hook. The hook's server entered at
_dispatch_under_rules (18282), which handles ONLY the ask rule — so a deny
rule on a read-only tool (e.g. read_file for a secrets path) and plan mode
are silently bypassed over the RPC. Red on f1ba3d461 until the server calls
executor.execute and the rule_evaluator plumbing is dropped.
"""
import json

from delfin.agent import code_mode


class _GateExecutor:
    """Proves the server calls execute() and nothing lower."""

    def __init__(self):
        self.via_execute = 0
        self.via_lower = 0

    def execute(self, name, arguments, permissions):
        self.via_execute += 1
        return json.dumps({"result": "ok"})

    def _dispatch_under_rules(self, name, arguments, permissions, rule_hit):
        self.via_lower += 1
        return json.dumps({"result": "ok"})

    def _dispatch(self, name, arguments, permissions):
        self.via_lower += 1
        return json.dumps({"result": "ok"})


class _DenyingExecutor:
    """Stands in for the real gate: execute() enforces a deny rule."""

    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions):
        self.calls.append(name)
        if name == "get_calc_info" and "denied" in json.dumps(arguments):
            return json.dumps({"error": (
                "get_calc_info is not allowed: path is denied by a managed "
                "deny rule")})
        return json.dumps({"result": "ok"})


def _server(executor, **kwargs):
    return code_mode.RpcServer(executor, **kwargs)


def test_server_calls_execute_not_the_lower_entries():
    """A stub proves the entry point: execute only, never
    _dispatch_under_rules or _dispatch."""
    ex = _GateExecutor()
    server = _server(ex)
    reply = json.loads(server.handle_payload(
        json.dumps({"tool": "search_calcs", "arguments": {}})))
    assert reply == {"result": {"result": "ok"}}
    assert ex.via_execute == 1
    assert ex.via_lower == 0, (
        "server bypassed the public entry: used a lower entry %d times"
        % ex.via_lower)


def test_deny_rule_on_read_only_tool_refuses_over_rpc_as_directly():
    """A managed deny on an ALLOW-LISTED read-only tool must refuse over
    the RPC exactly as a direct call — red while the server enters
    below execute(). read_file is not in the phase-1 allow-list, so it
    is refused earlier by the server; get_calc_info reaches execute()."""
    ex = _DenyingExecutor()
    server = _server(ex)
    reply = json.loads(server.handle_payload(json.dumps({
        "tool": "get_calc_info",
        "arguments": {"calc_id": "denied-id"},
    })))
    assert "error" in reply, "deny rule was bypassed over the RPC"
    assert "denied" in reply["error"]


def test_plan_mode_refusal_shape_survives_the_rpc():
    """Whatever execute() refuses (plan mode deny-by-default) surfaces as
    the program's error — the RPC adds no second, weaker entry."""
    class _PlanExecutor:
        def execute(self, name, arguments, permissions):
            return json.dumps({"error": (
                "plan mode is active: read-only until the plan is "
                "approved; make this call after approval")})
        # A lower entry that, if used, would NOT refuse:
        def _dispatch_under_rules(self, name, arguments, permissions,
                                  rule_hit):
            return json.dumps({"result": "leaked"})

    ex = _PlanExecutor()
    server = _server(ex)
    reply = json.loads(server.handle_payload(json.dumps({
        "tool": "write_file", "arguments": {"path": "x", "content": "y"},
    })))
    # write_file is refused by the server allow-list anyway; but the
    # shape check must come from execute()'s refusal, not silently pass.
    reply2 = json.loads(server.handle_payload(json.dumps({
        "tool": "search_docs", "arguments": {"query": "solvent"}})))
    # The server wrapper surfaces execute()'s refusal dict as the
    # top-level error payload (RpcServer.handle_payload error/result
    # split); the refusal must come from execute(), not "leaked".
    assert "plan mode is active" in str(reply2.get("error")), (
        "RPC did not surface execute()'s plan-mode refusal: %r" % reply2)
