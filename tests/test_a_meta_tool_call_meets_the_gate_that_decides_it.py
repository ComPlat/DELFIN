"""An MCP resource/prompt call is decided by the gate every other call meets.

The two meta-tools ``mcp_read_resource`` and ``mcp_get_prompt`` are the only
tool calls in a turn that reach a tool executor while skipping the decisions
the rest of the turn gets: the streaming loop dispatched them from their own
branch (``api_client.py``, ``elif is_mcp_meta:``), which asked the ROLE list
and then called the MCP registry directly.

Measured on cdadd9efb, driving ``OpenAIClient.stream_message`` with a scripted
stream (no network, no key):

* a ``deny ["mcp_read_resource"]`` in the project's ``.delfin/settings.json``
  -- a rule the SAME client refuses every other call for -- did not stop it,
  and the fake registry was called;
* ``--disallowed-tools mcp_read_resource`` (the session's own deny list) did
  not stop it either;
* a blocking PreToolUse hook on the name was never run.

What this test binds: the outcome for each case, through the product's own
entry -- ``stream_message``, the turn a model's tool call actually arrives at.
The control case (no rule, nobody refusing) asserts the behaviour that must
NOT change: the resource is read and its body reaches the model.

Rejected: calling ``_gate_mcp_tool`` (or ``_mcp_call_refusal``) directly. That
is the decider, not the path; a test that calls it proves the function exists
and says nothing about whether a turn reaches it -- which is the defect.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import mcp_client as M


# --- a fake MCP server that records what reached it -------------------------


class _FakeRegistry:
    def __init__(self):
        self.read_calls: list = []
        self.prompt_calls: list = []

    def discover_all(self):
        return []

    def discover_resources(self):
        return [M.MCPResource(server="remote", uri="file:///notes.md",
                              name="Notes")]

    def discover_prompts(self):
        return [M.MCPPrompt(server="remote", name="summarize",
                            description="Summarize text")]

    def read_resource(self, server, uri):
        self.read_calls.append((server, uri))
        return "RESOURCE-BODY"

    def get_prompt(self, name, arguments=None):
        self.prompt_calls.append((name, arguments))
        return "PROMPT-BODY"


# --- a scripted stream, so no model is needed -------------------------------


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


def _call_then_finish(name: str, args_json: str) -> _Stream:
    return _Stream([
        _Chunk([_Choice(_Delta(tool_calls=[
            _TCDelta(0, id="c1", name=name, arguments=args_json)]))]),
        _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
    ])


def _final_text() -> _Stream:
    return _Stream([
        _Chunk([_Choice(_Delta(content="done"), finish="stop")],
               usage=_Usage()),
    ])


# --- the turn ---------------------------------------------------------------


def _client(monkeypatch, tmp_path, *, settings=None):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    reg = _FakeRegistry()
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: reg)
    if settings is not None:
        conf = tmp_path / ".delfin"
        conf.mkdir(exist_ok=True)
        (conf / "settings.json").write_text(json.dumps(settings))
    from delfin.agent.api_client import create_client
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))
    return client, reg


def _turn_results(client, streams) -> list:
    """Run the turn; return the tool results the model was given."""
    calls = {"n": 0}

    def _create(**kwargs):
        s = streams[min(calls["n"], len(streams) - 1)]
        calls["n"] += 1
        return s

    client.client.chat.completions.create = _create
    out = []
    for ev in client.stream_message(
            "sys", [{"role": "user", "content": "go"}], max_tokens=100):
        if getattr(ev, "type", "") == "tool_result":
            out.append(str(getattr(ev, "tool_output", "")))
    return out


_RESOURCE_ARGS = json.dumps({"server": "remote", "uri": "file:///notes.md"})


def _one_tool_result(client, streams) -> str:
    results = _turn_results(client, streams)
    assert results, "the turn produced no tool result at all"
    return results[0]


# --- the cases --------------------------------------------------------------


def test_an_allowed_resource_read_still_runs(monkeypatch, tmp_path):
    """CONTROL: the behaviour that must not change.

    Passes on the base and must keep passing: a resource read with nobody
    refusing it reaches the server and its body comes back.
    """
    client, reg = _client(monkeypatch, tmp_path)
    result = _one_tool_result(
        client, [_call_then_finish("mcp_read_resource", _RESOURCE_ARGS),
                 _final_text()])
    assert reg.read_calls == [("remote", "file:///notes.md")], (
        f"the registry was not reached: {reg.read_calls}")
    assert "RESOURCE-BODY" in result, result


def test_a_deny_rule_on_a_resource_read_stops_it(monkeypatch, tmp_path):
    """The project's own deny rule, which every other call honours."""
    client, reg = _client(
        monkeypatch, tmp_path,
        settings={"permissions": {"deny": ["mcp_read_resource"]}})
    result = _one_tool_result(
        client, [_call_then_finish("mcp_read_resource", _RESOURCE_ARGS),
                 _final_text()])
    assert reg.read_calls == [], (
        "a deny rule on mcp_read_resource did not stop the call: the "
        f"server was reached with {reg.read_calls}")
    assert "mcp_read_resource" in result, result
    assert "RESOURCE-BODY" not in result, result


def test_the_sessions_own_deny_list_stops_a_resource_read(monkeypatch,
                                                          tmp_path):
    """``--disallowed-tools mcp_read_resource``: a name the user typed."""
    client, reg = _client(monkeypatch, tmp_path)
    client.narrow_tool_surface(denied=("mcp_read_resource",))
    result = _one_tool_result(
        client, [_call_then_finish("mcp_read_resource", _RESOURCE_ARGS),
                 _final_text()])
    assert reg.read_calls == [], (
        "the session's deny list did not stop mcp_read_resource: the server "
        f"was reached with {reg.read_calls}")
    assert "RESOURCE-BODY" not in result, result


def test_a_blocking_hook_stops_a_resource_read(monkeypatch, tmp_path):
    """A PreToolUse hook, the decision every other call in the turn gets."""
    from delfin.agent import api_client as ac

    seen: list = []

    def _blocks(name, args, perms):
        seen.append(str(name))
        return "hooks say no to that resource"

    client, reg = _client(monkeypatch, tmp_path)
    monkeypatch.setattr(ac._doc_executor, "_run_pre_tool_hooks", _blocks)
    result = _one_tool_result(
        client, [_call_then_finish("mcp_read_resource", _RESOURCE_ARGS),
                 _final_text()])
    assert seen, ("the PreToolUse hooks never ran for mcp_read_resource: "
                  "this call reaches the server without them")
    assert reg.read_calls == [], (
        f"a blocking hook did not stop the call: {reg.read_calls}")
    assert "RESOURCE-BODY" not in result, result


@pytest.mark.parametrize("tool,args,attr", [
    ("mcp_read_resource", _RESOURCE_ARGS, "read_calls"),
    ("mcp_get_prompt",
     json.dumps({"name": "mcp__remote__summarize",
                 "arguments": {"text": "X"}}), "prompt_calls"),
])
def test_neither_meta_tool_reaches_its_server_past_a_deny_rule(
        monkeypatch, tmp_path, tool, args, attr):
    """Both names, so a fix covering one spelling does not pass."""
    client, reg = _client(
        monkeypatch, tmp_path,
        settings={"permissions": {"deny": [tool]}})
    result = _one_tool_result(
        client, [_call_then_finish(tool, args), _final_text()])
    assert getattr(reg, attr) == [], (
        f"a deny rule on {tool} did not stop it: "
        f"{getattr(reg, attr)}")
    assert "BODY" not in result, result
