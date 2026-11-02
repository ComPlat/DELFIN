"""Red tests for the operator's blocking finding (public-entry parity).

The direct path enters at KitToolExecutor.execute (api_client.py:17651)
and evaluates managed/user deny rules, plan-mode, sandbox, role checks
and shape checks there. The RPC server previously called
_dispatch_under_rules directly, skipping all of it. These tests pin:

1. a stub executor records the entry: the server must call
   executor.execute(name, arguments, permissions) — nothing lower;
2. a deny rule on a read-only tool refuses over the RPC exactly as
   directly (parity, not just refusal);
3. plan-mode refusal over the RPC.

All three were RED on f1ba3d461 (the server called _dispatch_under_rules
and skipped the preamble), verified before the fix.
"""
import json

from delfin.agent import code_mode


class _EntryRecordingExecutor:
    """Stub that mimics execute's contract and records how it was called."""

    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions=None):
        self.calls.append(("execute", name, arguments, permissions))
        return json.dumps({"ok": True, "name": name})

    # The sub-function the old wiring called. If the server calls this,
    # the test below fails — the preamble would be skipped.
    def _dispatch_under_rules(self, name, arguments, permissions, rule_hit):
        self.calls.append(("_dispatch_under_rules", name, arguments,
                           permissions))
        return json.dumps({"ok": True, "name": name})


def _serve(executor, request):
    server = code_mode.RpcServer(executor, permissions=None)
    return server.handle_payload(json.dumps(request))


def test_server_calls_the_public_entry_not_the_sub_function():
    ex = _EntryRecordingExecutor()
    _serve(ex, {"tool": "get_calc_info",
                "arguments": {"folder": "x"}})
    entries = [c[0] for c in ex.calls]
    assert entries == ["execute"], entries
    assert ex.calls[0][1] == "get_calc_info"
    assert ex.calls[0][2] == {"folder": "x"}


def test_deny_rule_on_a_read_only_tool_refuses_like_the_direct_path(tmp_path):
    from delfin.agent.permission_rules import Rule, RuleSet

    # The same RuleSet handed to the direct executor is handed to the RPC
    # server; the PUBLIC entry must refuse the deny exactly as directly.
    # A tool rule is a bare name: Rule(text, "tool", <tool name>).
    rules = RuleSet(deny=(Rule("get_calc_info", "tool",
                               "get_calc_info"),))

    class _DirectLikeExecutor(_EntryRecordingExecutor):
        def __init__(self, rules):
            super().__init__()
            self._rules = rules

        def execute(self, name, arguments, permissions=None):
            self.calls.append(("execute", name, arguments, permissions))
            from delfin.agent import permission_rules as pr_mod
            hit = pr_mod.friction(self._rules, name, arguments,
                                  workspace=None)
            if hit is not None and hit[0] == "deny":
                return json.dumps({"error": "denied by rule"})
            return json.dumps({"ok": True, "name": name})

    ex = _DirectLikeExecutor(rules)
    direct = ex.execute("get_calc_info", {"folder": "a"})
    over_rpc = _serve(ex, {"tool": "get_calc_info",
                           "arguments": {"folder": "a"}})
    assert json.loads(direct)["error"] == json.loads(over_rpc)["error"]
    assert json.loads(over_rpc)["error"] == "denied by rule"
    # And the RPC went through the public entry, not the sub-function.
    assert [c[0] for c in ex.calls] == ["execute", "execute"]


def test_plan_mode_refusal_reaches_the_rpc_client():
    class _PlanModeExecutor(_EntryRecordingExecutor):
        def execute(self, name, arguments, permissions=None):
            self.calls.append(("execute", name, arguments, permissions))
            return json.dumps({"error": "plan mode: read-only, "
                                        "refused — make the direct call "
                                        "instead"})

    ex = _PlanModeExecutor()
    reply = json.loads(_serve(ex, {"tool": "get_calc_info",
                                   "arguments": {"folder": "x"}}))
    assert "refused" in reply["error"]
