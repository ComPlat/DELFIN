"""ACTION-protocol repair helpers: detection, extraction, repair results.

Covers the pure helpers in ``delfin/agent/action_protocol.py`` and — once
the repair branch is wired into the streaming tool loop — an end-to-end
check that an ACTION-style tool call no longer trips the
consecutive-identical-error abort and is registered on the text channel.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import action_protocol as ap


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [
    "ACTION", "action", "Action", "ACTION:",
    "mcp__delfin-docs__ACTION", "mcp__anything__action",
    "/tab", "/orca set functional B3LYP", "/recalc auto",
    "ACTION:/tab calc",
])
def test_action_style_names_are_detected(name):
    assert ap.is_action_style_call(name)


@pytest.mark.parametrize("name", [
    "search_docs", "read_file", "bash", "mcp__kit-coding__bash",
    "actions", "faction", "", "/", "//x",
    "/home/user/file.py",          # path lookalike, not a slash command
    "mcp_read_resource",
])
def test_ordinary_tool_names_are_not_detected(name):
    assert not ap.is_action_style_call(name)


def test_role_gate():
    assert ap.role_uses_action_protocol("dashboard_agent")
    assert not ap.role_uses_action_protocol("")
    assert not ap.role_uses_action_protocol("solo_agent")


def test_strip_tool_namespace():
    assert ap.strip_tool_namespace("mcp__delfin-docs__ACTION") == "ACTION"
    assert ap.strip_tool_namespace("ACTION") == "ACTION"
    assert ap.strip_tool_namespace("mcp__a__b__c") == "c"


# ---------------------------------------------------------------------------
# Extraction — liberal in what it accepts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name,args,expected", [
    # The confirmed field shape: tool named ACTION, command in arguments.
    ("ACTION", {"command": "/orca set functional B3LYP"},
     "/orca set functional B3LYP"),
    ("mcp__delfin-docs__ACTION", {"command": "/tab calc"}, "/tab calc"),
    # Alternative argument keys.
    ("ACTION", {"action": "/tab calc"}, "/tab calc"),
    ("ACTION", {"input": "/submit"}, "/submit"),
    ("ACTION", {"cmd": "/control set functional B3LYP"},
     "/control set functional B3LYP"),
    # First string value under an unknown key.
    ("ACTION", {"whatever": "/tab jobs"}, "/tab jobs"),
    # Redundant ACTION prefix / quoting / fencing inside the value.
    ("ACTION", {"command": "ACTION: /tab calc"}, "/tab calc"),
    ("ACTION", {"command": "`/tab calc`"}, "/tab calc"),
    ("ACTION", {"command": '"/tab calc"'}, "/tab calc"),
    # Missing slash under an imperative key is promoted.
    ("ACTION", {"command": "tab calc"}, "/tab calc"),
    # One nesting level.
    ("ACTION", {"input": {"command": "/tab calc"}}, "/tab calc"),
    # Arguments arriving as a JSON string or as plain text.
    ("ACTION", '{"command": "/tab calc"}', "/tab calc"),
    ("ACTION", "ACTION: /submit", "/submit"),
    ("ACTION", '"/tab calc"', "/tab calc"),
    # The tool NAME is the command.
    ("/orca set functional B3LYP", {}, "/orca set functional B3LYP"),
    ("/tab", {"command": "calc"}, "/tab calc"),
    ("/orca", {"args": "set functional B3LYP"},
     "/orca set functional B3LYP"),
    ("ACTION:/tab calc", {}, "/tab calc"),
    # Name and args both slash-shaped: the more complete one wins.
    ("/tab", {"command": "/tab calc"}, "/tab calc"),
])
def test_extract_slash_command(name, args, expected):
    assert ap.extract_slash_command(name, args) == expected


@pytest.mark.parametrize("name,args", [
    ("ACTION", {}),
    ("ACTION", None),
    ("ACTION", {"note": "please open the calculations tab"}),
    # Prose-carrying keys must not fabricate a command.
    ("ACTION", {"query": "how do I use /recalc here"}),
])
def test_extract_returns_empty_when_no_command(name, args):
    assert ap.extract_slash_command(name, args) == ""


def test_canonical_line_encodes_newlines():
    line = ap.canonical_action_line("/control set input\nline2")
    assert line == "ACTION: /control set input\\nline2"
    assert "\n" not in line


# ---------------------------------------------------------------------------
# Repair result — constructive, never error-shaped, never byte-identical
# ---------------------------------------------------------------------------


def test_repair_result_is_not_error_shaped():
    for cmd in ("/tab calc", ""):
        for attempt in (1, 2, 5):
            res = ap.build_repair_result(cmd, attempt=attempt)
            assert not res.lstrip().startswith('{"error"')
            payload = json.loads(res)
            assert next(iter(payload)) == "repaired"


def test_repair_result_states_interpretation_and_text_protocol():
    res = ap.build_repair_result(
        "/orca set functional B3LYP", role="dashboard_agent",
        registered=True, tool_name="ACTION")
    payload = json.loads(res)
    assert payload["repaired"] is True
    assert payload["understood_as"] == "ACTION: /orca set functional B3LYP"
    assert payload["registered"] is True
    assert "plain text" in payload["guidance"]
    assert "ACTION: /orca set functional B3LYP" in payload["guidance"]


def test_repair_result_varies_per_attempt():
    results = [
        ap.build_repair_result("/tab calc", attempt=n, registered=True)
        for n in (1, 2, 3)
    ]
    assert len(set(results)) == 3
    # Later occurrences must say the command is already registered.
    assert "already registered" in json.loads(results[1])["guidance"]


def test_repair_result_without_command_gives_example():
    res = ap.build_repair_result("", attempt=2)
    payload = json.loads(res)
    assert payload["repaired"] is False
    assert payload["understood_as"] is None
    assert payload["registered"] is False
    assert "ACTION: /tab calc" in payload["guidance"]
    # Attempt counter keeps repeats non-identical.
    assert res != ap.build_repair_result("", attempt=3)


# ---------------------------------------------------------------------------
# Loop integration — active once the repair branch is wired into
# api_client's streaming tool loop; skipped until then.
# ---------------------------------------------------------------------------


def _repair_wired() -> bool:
    from delfin.agent import api_client as ac
    return "action_protocol" in Path(ac.__file__).read_text(encoding="utf-8")


_REPAIR_WIRED = _repair_wired()


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


def _tool_round(tool_name, args_json, tc_id="c1"):
    return _Stream([
        _Chunk([_Choice(_Delta(tool_calls=[
            _TCDelta(0, id=tc_id, name=tool_name, arguments=args_json)]))]),
        _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
    ])


def _final():
    return _Stream([
        _Chunk([_Choice(_Delta(content="done"), finish="stop")],
               usage=_Usage()),
    ])


class _EmptyReg:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


@pytest.fixture
def dashboard_client(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    from delfin.agent import mcp_client as M
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    from delfin.agent.api_client import KitToolPermissions, create_client
    workspace = tmp_path / "ws"
    workspace.mkdir()
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(workspace))
    client.set_permissions(KitToolPermissions(
        workspace=workspace, agent_role="dashboard_agent"))
    return client


def _drive(client, streams):
    calls = {"n": 0}

    def _create(**kwargs):
        s = streams[min(calls["n"], len(streams) - 1)]
        calls["n"] += 1
        return s

    client.client.chat.completions.create = _create
    return list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))


@pytest.mark.skipif(
    not _REPAIR_WIRED,
    reason="ACTION repair branch not yet applied to api_client.py",
)
def test_repeated_action_tool_calls_do_not_abort_the_turn(dashboard_client):
    args = json.dumps({"command": "/orca set functional B3LYP"})
    events = _drive(dashboard_client, [
        _tool_round("ACTION", args, tc_id="c1"),
        _tool_round("ACTION", args, tc_id="c2"),
        _tool_round("ACTION", args, tc_id="c3"),
        _final(),
    ])
    stops = [e.stop_reason for e in events
             if getattr(e, "type", "") == "message_delta" and e.stop_reason]
    assert "consecutive_identical_errors" not in stops
    # The turn reached the final text instead of aborting.
    assert any(getattr(e, "type", "") == "text_delta"
               and "done" in (e.text or "") for e in events)
    # Tool results are constructive, non-error, and vary per occurrence.
    results = [e.tool_output for e in events
               if getattr(e, "type", "") == "tool_result"]
    assert len(results) == 3
    assert all(not r.lstrip().startswith('{"error"') for r in results)
    assert len(set(results)) == 3
    # The repaired command is registered on the text channel exactly once.
    text = "".join(e.text for e in events
                   if getattr(e, "type", "") == "text_delta" and e.text)
    assert text.count("ACTION: /orca set functional B3LYP") == 1


@pytest.mark.skipif(
    not _REPAIR_WIRED,
    reason="ACTION repair branch not yet applied to api_client.py",
)
def test_action_call_without_command_gets_guidance_not_error(dashboard_client):
    events = _drive(dashboard_client, [
        _tool_round("ACTION", "{}"),
        _final(),
    ])
    results = [e.tool_output for e in events
               if getattr(e, "type", "") == "tool_result"]
    assert len(results) == 1
    assert not results[0].lstrip().startswith('{"error"')
    assert "plain text" in results[0] or "ACTION: /tab calc" in results[0]
    # Nothing was registered on the text channel.
    text = "".join(e.text for e in events
                   if getattr(e, "type", "") == "text_delta" and e.text)
    assert "ACTION: /" not in text
