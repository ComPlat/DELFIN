"""A PreToolUse hook that blocks a call asks nobody first.

The sessions reach bash over MCP (kit-coding:bash), and on that path the
hooks and the gate were both evaluated before either result was looked
at: a hook that refused a call still put the gate's confirmation dialog
in front of the user, who decided for nothing before the hook's refusal
won anyway. The native path has always short-circuited.
"""

from __future__ import annotations

import json

from delfin.agent import api_client as A


def _patch(monkeypatch, hook_reason):
    calls = {"gate": 0}
    monkeypatch.setattr(A._doc_executor, "_run_pre_tool_hooks",
                        lambda base, args, perms: hook_reason)

    def gate(name, args, perms):
        calls["gate"] += 1
        return None
    monkeypatch.setattr(A._doc_executor, "_gate_mcp_tool", gate)
    return calls


def test_a_blocking_hook_answers_and_the_gate_never_asks(monkeypatch):
    calls = _patch(monkeypatch, "use the project's test runner")
    out = A._mcp_call_refusal("mcp__kit-coding__bash",
                              {"command": "python -m pytest"}, object())
    assert calls["gate"] == 0, "the gate (and its dialog) ran anyway"
    body = json.loads(out)
    assert body["error"] == "blocked_by_hook"
    assert body["reason"] == "use the project's test runner"


def test_without_a_hook_block_the_gate_decides_as_before(monkeypatch):
    calls = _patch(monkeypatch, None)
    assert A._mcp_call_refusal("mcp__kit-coding__bash",
                               {"command": "ls"}, object()) is None
    assert calls["gate"] == 1


def test_a_gate_refusal_still_refuses(monkeypatch):
    monkeypatch.setattr(A._doc_executor, "_run_pre_tool_hooks",
                        lambda base, args, perms: None)
    monkeypatch.setattr(A._doc_executor, "_gate_mcp_tool",
                        lambda name, args, perms: "user denied it")
    out = A._mcp_call_refusal("mcp__kit-coding__bash",
                              {"command": "rm -rf x"}, object())
    assert json.loads(out) == {"error": "user denied it"}


def test_the_hook_sees_the_base_name(monkeypatch):
    seen = {}
    monkeypatch.setattr(A._doc_executor, "_run_pre_tool_hooks",
                        lambda base, args, perms: seen.setdefault("b", base)
                        and None)
    monkeypatch.setattr(A._doc_executor, "_gate_mcp_tool",
                        lambda name, args, perms: None)
    A._mcp_call_refusal("mcp__kit-coding__bash", {"command": "ls"}, object())
    assert seen["b"] == "bash"
