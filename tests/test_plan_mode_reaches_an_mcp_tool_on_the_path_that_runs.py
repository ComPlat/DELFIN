"""Plan mode was enforced at one site the MCP calls never passed.

The deny-by-default plan check lives in ``_DocToolExecutor.execute()``. The
streaming loop routes every ``mcp__*`` name AROUND ``execute()``: it runs
the pre-tool hooks and ``_gate_mcp_tool``, then dispatches to the registry.
``_gate_mcp_tool`` re-implemented the role check and three families (shell
bases, network tools, the write map) and returned None for everything else
— with no plan-mode check anywhere in it. So in the profile labelled
read-only, ``mcp__github__create_pull_request``, ``mcp__jira__create_issue``
and ``mcp__db__delete_row`` all dispatched.

The test that was supposed to cover this called
``ex.execute("mcp__github__create_pull_request")`` — a path production
never takes, so it passed against the hole. These tests drive the branch
that really runs: the streaming loop with a stubbed model, asserting the
registry is never called, plus the gate the branch consults.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import mcp_client as M
from delfin.agent.api_client import KitToolPermissions, create_client, _doc_executor


# ---------------------------------------------------------------------------
# The gate the dispatch branch consults
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", [
    "mcp__github__create_pull_request",
    "mcp__jira__create_issue",
    "mcp__db__delete_row",
    "mcp__slack__post_message",
    "mcp__kit-coding__write_file",
    "mcp__kit-coding__bash",
])
def test_a_side_effecting_mcp_tool_is_refused_in_plan_mode(tmp_path, tool):
    perms = KitToolPermissions(workspace=tmp_path, mode="plan")
    blocked = _doc_executor._gate_mcp_tool(tool, {}, perms)
    assert blocked is not None and "plan mode" in blocked, tool


@pytest.mark.parametrize("tool", [
    "mcp__coding__read_file",
    "mcp__docs__search_docs",
    "mcp__coding__grep_file",
])
def test_investigating_through_a_server_still_works(tmp_path, tool):
    perms = KitToolPermissions(workspace=tmp_path, mode="plan")
    assert _doc_executor._gate_mcp_tool(tool, {}, perms) is None, tool


def test_a_dry_run_is_still_a_dry_run(tmp_path):
    """The same exemption the native check makes: check_only changes
    nothing, so plan mode has no reason to refuse it."""
    perms = KitToolPermissions(workspace=tmp_path, mode="plan")
    assert _doc_executor._gate_mcp_tool(
        "mcp__github__create_pull_request", {"check_only": True}, perms) is None


@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_no_other_mode_is_affected(tmp_path, mode):
    perms = KitToolPermissions(workspace=tmp_path, mode=mode)
    out = _doc_executor._gate_mcp_tool(
        "mcp__github__create_pull_request", {}, perms)
    assert out is None or "plan mode" not in out


def test_a_session_without_permissions_is_not_refused_for_the_mode(tmp_path):
    """No permissions means no mode, so the PLAN check cannot fire.

    This asserted ``is None`` — that a session without permissions passed
    the gate untouched — which was true only because the gate ended in a
    blanket ``return None``. With deny-by-default the call is still
    refused, for the reason it deserves: ``_permissions=None`` is the
    profile that "disables write/edit/bash", and reaching a server's own
    mutation from inside it is the same leak by another route. What the
    test is about — that plan mode is not blamed for it — is unchanged.
    """
    out = _doc_executor._gate_mcp_tool(
        "mcp__github__create_pull_request", {}, None)
    assert out is not None
    assert "plan mode" not in out
    assert "no permissions are configured" in out


# ---------------------------------------------------------------------------
# The branch that really runs
# ---------------------------------------------------------------------------

class _FnDelta:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TCDelta:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _FnDelta(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _Registry:
    """Records any dispatch that gets past the gate."""

    def __init__(self):
        self.calls: list[tuple] = []

    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []

    def call(self, name, args):
        self.calls.append((name, args))
        return json.dumps({"status": "created"})


def _run_tool_call(tmp_path, monkeypatch, tool_name, mode):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    registry = _Registry()
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: registry)
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))
    client._permissions = KitToolPermissions(workspace=tmp_path, mode=mode)

    streams = [
        _Stream([
            _Chunk([_Choice(_Delta(tool_calls=[_TCDelta(
                0, id="c1", name=tool_name,
                arguments='{"title": "x"}')]))]),
            _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
        ]),
        _Stream([_Chunk([_Choice(_Delta(content="done"), finish="stop")],
                        usage=_Usage())]),
    ]
    seen = {"n": 0}

    def _create(**kwargs):
        s = streams[min(seen["n"], len(streams) - 1)]
        seen["n"] += 1
        return s

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))
    results = [getattr(e, "tool_output", "") or ""
               for e in events if getattr(e, "type", "") == "tool_result"]
    return registry, results


def test_the_dispatch_branch_refuses_and_never_calls_the_server(
        tmp_path, monkeypatch):
    registry, results = _run_tool_call(
        tmp_path, monkeypatch, "mcp__github__create_pull_request", "plan")
    assert results, "no tool_result emitted"
    assert any("plan mode" in r for r in results), results
    assert registry.calls == [], "the MCP server was called in plan mode"


def test_the_same_branch_still_dispatches_a_read(tmp_path, monkeypatch):
    registry, results = _run_tool_call(
        tmp_path, monkeypatch, "mcp__coding__read_file", "plan")
    assert registry.calls, "a read-only MCP tool was refused in plan mode"
