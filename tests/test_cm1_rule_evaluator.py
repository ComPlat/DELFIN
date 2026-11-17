"""Rules live in the PUBLIC entry, not in the RPC layer.

Old decision this file pinned (P1d): the RPC layer evaluated permission
rules per call via a rule_evaluator callable. SUPERSEDED by the operator's
blocking finding (public-entry parity): the server must call
executor.execute(name, arguments, permissions) — the single public entry —
which evaluates deny rules, plan-mode, sandbox and shape checks itself.
The RPC layer must carry NO rule plumbing at all. These tests pin the new
contract.
"""
import json

from delfin.agent import code_mode


class _RecordingExecutor:
    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions=None):
        self.calls.append(("execute", name, arguments, permissions))
        return json.dumps({"ok": True, "name": name})

    def _dispatch_under_rules(self, name, arguments, permissions, rule_hit):
        self.calls.append(("_dispatch_under_rules", name, arguments,
                           permissions))
        return json.dumps({"ok": True, "name": name})

    # Any other rule plumbing the old wiring could have called.
    def __call__(self, name, arguments, permissions, rule_hit):
        self.calls.append(("callable", name, arguments, permissions))
        return json.dumps({"ok": True, "name": name})


def _serve(ex, tool, arguments):
    server = code_mode.RpcServer(ex, permissions=None)
    return json.loads(server.handle_payload(
        json.dumps({"tool": tool, "arguments": arguments})))


def test_rpc_layer_passes_no_rule_argument_to_the_entry():
    """The public entry receives exactly (name, arguments, permissions) —
    no rule_hit, no evaluator result."""
    ex = _RecordingExecutor()
    _serve(ex, "get_calc_info", {"folder": "x"})
    assert ex.calls == [("execute", "get_calc_info", {"folder": "x"}, None)]


def test_rpc_layer_never_touches_lower_entries():
    """Only execute() is called — never _dispatch_under_rules or the old
    callable protocol."""
    ex = _RecordingExecutor()
    _serve(ex, "get_calc_info", {"folder": "x"})
    assert [c[0] for c in ex.calls] == ["execute"]


def test_rpc_layer_carries_no_rule_evaluator_attribute():
    """The plumbing is gone: neither Connection nor RpcServer accepts or
    stores a rule evaluator."""
    import inspect
    for cls in (code_mode.Connection, code_mode.RpcServer):
        params = inspect.signature(cls.__init__).parameters
        assert not any("rule" in p for p in params), (cls, list(params))
