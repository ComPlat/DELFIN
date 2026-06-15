"""The MCP resource/prompt meta-tools (mcp_read_resource / mcp_get_prompt).

Verifies that when connected MCP servers expose resources/prompts, the agent
advertises the two meta-tools, and that a tool call routes through the MCP
registry. No live model/server: the OpenAI stream and the MCP registry are both
faked.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import mcp_client as M


# ---------------------------------------------------------------------------
# Fake MCP registry
# ---------------------------------------------------------------------------


class _FakeRegistry:
    def __init__(self):
        self.read_calls = []
        self.prompt_calls = []

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


# ---------------------------------------------------------------------------
# Fake OpenAI streaming
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


def _tool_call_then_finish(tool_name, args_json):
    """Round 1: emit one tool_call, then finish_reason=tool_calls."""
    return _Stream([
        _Chunk([_Choice(_Delta(tool_calls=[
            _TCDelta(0, id="c1", name=tool_name, arguments=args_json)]))]),
        _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
    ])


def _final_text():
    return _Stream([
        _Chunk([_Choice(_Delta(content="done"), finish="stop")],
               usage=_Usage()),
    ])


@pytest.fixture
def _client(monkeypatch, tmp_path):
    # Resolve caps as a normal tool-capable model (avoid network + no-tools gate)
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    reg = _FakeRegistry()
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: reg)
    from delfin.agent.api_client import create_client
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))
    return client, reg


def _run(client, streams):
    """Wire a stateful fake create() returning the given streams in order."""
    captured = {"tools_round1": None}
    calls = {"n": 0}

    def _create(**kwargs):
        if calls["n"] == 0:
            captured["tools_round1"] = kwargs.get("tools")
        s = streams[min(calls["n"], len(streams) - 1)]
        calls["n"] += 1
        return s

    client.client.chat.completions.create = _create
    list(client.stream_message("sys", [{"role": "user", "content": "go"}],
                               max_tokens=100))
    return captured


def _tool_names(tools):
    return {t["function"]["name"] for t in (tools or [])}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_meta_tools_advertised_when_resources_and_prompts_exist(_client):
    client, _ = _client
    captured = _run(client, [_final_text()])
    names = _tool_names(captured["tools_round1"])
    assert "mcp_read_resource" in names
    assert "mcp_get_prompt" in names


def test_read_resource_routes_through_registry(_client):
    client, reg = _client
    args = json.dumps({"server": "remote", "uri": "file:///notes.md"})
    _run(client, [_tool_call_then_finish("mcp_read_resource", args),
                  _final_text()])
    assert reg.read_calls == [("remote", "file:///notes.md")]


def test_get_prompt_routes_through_registry(_client):
    client, reg = _client
    args = json.dumps({"name": "mcp__remote__summarize",
                       "arguments": {"text": "X"}})
    _run(client, [_tool_call_then_finish("mcp_get_prompt", args),
                  _final_text()])
    assert reg.prompt_calls == [("mcp__remote__summarize", {"text": "X"})]


def test_meta_tools_absent_when_no_resources_or_prompts(monkeypatch, tmp_path):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)

    class _Empty(_FakeRegistry):
        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _Empty())
    from delfin.agent.api_client import create_client
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))
    captured = _run(client, [_final_text()])
    names = _tool_names(captured["tools_round1"])
    assert "mcp_read_resource" not in names
    assert "mcp_get_prompt" not in names
