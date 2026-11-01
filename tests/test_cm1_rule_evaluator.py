"""Pin for CM1 reviewer FINDING C (rule evaluation must be per call).

A DIRECT call evaluates permission rules per call (api_client.py:17772
_permission_rule_hit) and asks on an ask-rule hit (_dispatch_under_rules
api_client.py:18288-18295). The CM1 hook passed rule_hit=None into
RpcServer, so every Code Mode call skipped rule evaluation — a user's
ask rule on ANY tool name was silently bypassed over the RPC.

Contract now: RpcServer accepts a rule_evaluator callable
(perms, name, args) -> rule_hit (same shape _permission_rule_hit
returns); it evaluates per request and hands the FRESH result to the
executor. The red test: the executor stub must NOT always receive None.
Seen red on c1e81ce86 (executor stub received constant None).
"""
import json

from delfin.agent import code_mode


class _StubExecutor:
    def __init__(self):
        self.calls = []

    def __call__(self, name, arguments, permissions, rule_hit):
        self.calls.append((name, arguments, permissions, rule_hit))
        return json.dumps({"ok": name})


def test_rule_evaluator_runs_per_call_and_result_reaches_the_executor():
    ex = _StubExecutor()
    hits = iter(["hit1", "hit2"])
    server = code_mode.RpcServer(
        ex, permissions="perms",
        rule_evaluator=lambda perms, name, args: next(hits))
    server.handle_payload(json.dumps(
        {"tool": "get_calc_info", "arguments": {"folders": "a"}}))
    server.handle_payload(json.dumps(
        {"tool": "get_calc_info", "arguments": {"folders": "b"}}))
    assert ex.calls[0][3] == "hit1"
    assert ex.calls[1][3] == "hit2"


def test_no_evaluator_stays_none_for_tests_without_rules():
    ex = _StubExecutor()
    server = code_mode.RpcServer(ex, permissions="perms")
    server.handle_payload(json.dumps(
        {"tool": "get_calc_info", "arguments": {}}))
    assert ex.calls[0][3] is None


def test_evaluator_error_fails_closed_to_refusal():
    ex = _StubExecutor()

    def boom(perms, name, args):
        raise RuntimeError("rule store unreadable")

    server = code_mode.RpcServer(ex, permissions="perms",
                                 rule_evaluator=boom)
    reply = json.loads(server.handle_payload(json.dumps(
        {"tool": "get_calc_info", "arguments": {}})))
    assert "error" in reply
    assert ex.calls == []
