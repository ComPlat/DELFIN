"""Tests for delfin.agent.subagents (sub-agent delegation)."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from delfin.agent import subagents as SA
from delfin.agent.api_client import (
    KitToolPermissions, _DocToolExecutor,
)


def _fake_event(**kwargs):
    return SimpleNamespace(
        type=kwargs.get("type", "text_delta"),
        text=kwargs.get("text", ""),
        tool_name=kwargs.get("tool_name", ""),
        tool_input=kwargs.get("tool_input", ""),
        tool_output=kwargs.get("tool_output", ""),
        input_tokens=kwargs.get("input_tokens", 0),
        output_tokens=kwargs.get("output_tokens", 0),
    )


def _fake_client(text: str, tool_calls: int = 0):
    """Return a MagicMock client that streams a final message."""
    client = MagicMock()
    events = []
    for i in range(tool_calls):
        events.append(_fake_event(
            type="tool_use",
            tool_name=f"read_file_{i}",
            tool_input='{"path": "x"}',
        ))
    events.append(_fake_event(type="text_delta", text=text))
    events.append(_fake_event(
        type="message_delta", input_tokens=42, output_tokens=8,
    ))
    client.stream_message = MagicMock(return_value=iter(events))
    client._permissions = None
    client.set_permissions = MagicMock()
    return client


def test_list_subagents_returns_all_presets(monkeypatch, tmp_path):
    # Pin Path.home to a fresh tmp dir + force a reload so leaked user
    # fixtures from sibling test files don't poison this assertion.
    from pathlib import Path as _P
    monkeypatch.setattr(_P, "home", lambda: tmp_path)
    SA.reload_subagent_presets()
    out = SA.list_subagents()
    names = {p["name"] for p in out}
    # The four built-ins are always present; pack-shipped md files (e.g.
    # chemistry-reviewer) may register additional names, so use subset.
    assert {"explore", "plan", "code-reviewer", "general-purpose"} <= names


def test_run_subagent_unknown_type():
    res = SA.run_subagent(
        subagent_type="garbage",
        description="x",
        prompt="x",
        parent_client=_fake_client("hi"),
        parent_perms=None,
    )
    assert "unknown subagent_type" in res.error
    assert res.final_text == ""


def test_run_subagent_returns_text_and_tokens():
    client = _fake_client("found 3 places where X is used")
    res = SA.run_subagent(
        subagent_type="explore",
        description="find X",
        prompt="Where is X defined?",
        parent_client=client,
        parent_perms=None,
    )
    assert "found 3 places" in res.final_text
    assert res.input_tokens == 42
    assert res.output_tokens == 8
    assert res.error == ""


def test_run_subagent_records_tool_calls():
    client = _fake_client("done", tool_calls=3)
    res = SA.run_subagent(
        subagent_type="explore",
        description="t",
        prompt="explore stuff",
        parent_client=client,
        parent_perms=None,
    )
    assert len(res.tool_calls) == 3
    assert res.tool_calls[0]["name"] == "read_file_0"


def test_run_subagent_clones_perms_with_correct_mode():
    """An explore sub-agent's perms must be 'plan' even if parent is bypassPermissions."""
    seen_modes: list[str | None] = []

    def fake_set_perms(p):
        seen_modes.append(p.mode if p else None)

    with tempfile.TemporaryDirectory() as d:
        parent = KitToolPermissions(
            workspace=Path(d), mode="bypassPermissions",
        )
        client = _fake_client("ok")
        client.set_permissions.side_effect = fake_set_perms
        SA.run_subagent(
            subagent_type="explore",
            description="t",
            prompt="hi there test prompt",
            parent_client=client,
            parent_perms=parent,
        )
        # First call sets to sub-agent mode; second restores parent's saved mode.
        assert seen_modes[0] == "plan"
    # Parent unchanged
    assert parent.mode == "bypassPermissions"


def test_run_subagent_truncates_on_tool_call_budget():
    # 35 tool calls > the default 30 budget
    client = _fake_client("too many", tool_calls=35)
    res = SA.run_subagent(
        subagent_type="explore",
        description="t",
        prompt="hi there test prompt",
        parent_client=client,
        parent_perms=None,
        max_tool_calls=5,
    )
    assert res.truncated is True
    assert "tool-call budget" in res.error


def test_subagent_tool_dispatch_calls_runner():
    """The 'subagent' tool dispatches to perms.subagent_runner."""
    captured = {}

    def runner(*, subagent_type, description, prompt, isolation=""):
        captured["type"] = subagent_type
        captured["desc"] = description
        captured["prompt"] = prompt
        captured["isolation"] = isolation
        return {"result": "delegated", "subagent_type": subagent_type}

    with tempfile.TemporaryDirectory() as d:
        perms = KitToolPermissions(workspace=Path(d), mode="default")
        perms.subagent_runner = runner
        out = _DocToolExecutor().execute(
            "subagent",
            {
                "subagent_type": "explore",
                "description": "find audit-log uses",
                "prompt": "Search for callers of make_record() and report them.",
            },
            permissions=perms,
        )
        payload = json.loads(out)
    assert payload["result"] == "delegated"
    assert captured["type"] == "explore"
    assert captured["desc"] == "find audit-log uses"


def test_subagent_tool_rejects_unknown_type():
    with tempfile.TemporaryDirectory() as d:
        perms = KitToolPermissions(workspace=Path(d), mode="default")
        perms.subagent_runner = lambda **kw: {}
        out = _DocToolExecutor().execute(
            "subagent",
            {
                "subagent_type": "fancy",
                "description": "x",
                "prompt": "long enough briefing text here",
            },
            permissions=perms,
        )
        payload = json.loads(out)
    assert "error" in payload
    assert "available" in payload


def test_subagent_tool_rejects_short_prompt():
    with tempfile.TemporaryDirectory() as d:
        perms = KitToolPermissions(workspace=Path(d), mode="default")
        perms.subagent_runner = lambda **kw: {}
        out = _DocToolExecutor().execute(
            "subagent",
            {
                "subagent_type": "explore",
                "description": "x",
                "prompt": "tiny",
            },
            permissions=perms,
        )
        payload = json.loads(out)
    assert "error" in payload
    assert ">=20 chars" in payload["error"]


def test_subagent_tool_no_runner_attached():
    with tempfile.TemporaryDirectory() as d:
        perms = KitToolPermissions(workspace=Path(d), mode="default")
        # subagent_runner stays None
        out = _DocToolExecutor().execute(
            "subagent",
            {
                "subagent_type": "explore",
                "description": "x",
                "prompt": "long enough briefing text here",
            },
            permissions=perms,
        )
        payload = json.loads(out)
    assert "error" in payload
    assert "runner not attached" in payload["error"]


# ---------------------------------------------------------------------------
# Memory / instruction inheritance (inherit_context)
# ---------------------------------------------------------------------------


def _child_system_prompt(client) -> str:
    """The system prompt the child actually streamed with."""
    assert client.stream_message.call_args is not None
    return client.stream_message.call_args.kwargs["system"]


def _inheritance_workspace(tmp_path, monkeypatch):
    """A workspace with a DELFIN.md rule + two typed memories under a
    faked home, so both inheritance blocks have real content."""
    from pathlib import Path as _P
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(_P, "home", lambda: home)
    ws = tmp_path / "proj"
    ws.mkdir()
    (ws / "DELFIN.md").write_text(
        "Always run def2-TZVP single points after optimization.",
        encoding="utf-8")
    from delfin.agent import memory_store as MS
    MS.save_typed_memory(
        "Prefer ORCA 6.1 keywords for every DFT input deck.",
        repo_root=ws, memory_type="feedback", title="orca-keywords")
    MS.save_typed_memory(
        "Dashboard colours should stay colour-blind safe.",
        repo_root=ws, memory_type="user", title="dashboard-colours")
    return ws


def test_inherit_context_blocks_in_child_system_prompt(tmp_path, monkeypatch):
    ws = _inheritance_workspace(tmp_path, monkeypatch)
    perms = KitToolPermissions(workspace=ws, mode="default")
    client = _fake_client("done exploring")
    res = SA.run_subagent(
        subagent_type="explore",
        description="check dft conventions",
        prompt="Which ORCA keywords does this project use for DFT?",
        parent_client=client,
        parent_perms=perms,
    )
    assert res.error == ""
    sys_prompt = _child_system_prompt(client)
    assert "Project rules (inherited" in sys_prompt
    assert "def2-TZVP single points" in sys_prompt
    assert "Standing memory (inherited)" in sys_prompt
    assert "ORCA 6.1 keywords" in sys_prompt


def test_inherit_context_ranks_memories_by_prompt_relevance(
        tmp_path, monkeypatch):
    ws = _inheritance_workspace(tmp_path, monkeypatch)
    perms = KitToolPermissions(workspace=ws, mode="default")
    client = _fake_client("ok")
    SA.run_subagent(
        subagent_type="explore",
        description="dft keyword check",
        prompt="Which ORCA keywords does this project use for DFT?",
        parent_client=client,
        parent_perms=perms,
    )
    sys_prompt = _child_system_prompt(client)
    block = sys_prompt.split("Standing memory (inherited)", 1)[1]
    # BM25 puts the ORCA memory (prompt overlap) before the dashboard one.
    assert block.index("orca-keywords") < block.index("dashboard-colours")


def test_inherit_context_false_omits_blocks(tmp_path, monkeypatch):
    ws = _inheritance_workspace(tmp_path, monkeypatch)
    perms = KitToolPermissions(workspace=ws, mode="default")
    client = _fake_client("done")
    SA.run_subagent(
        subagent_type="explore",
        description="no context please",
        prompt="Which ORCA keywords does this project use for DFT?",
        parent_client=client,
        parent_perms=perms,
        inherit_context=False,
    )
    sys_prompt = _child_system_prompt(client)
    assert "Project rules (inherited" not in sys_prompt
    assert "Standing memory (inherited)" not in sys_prompt


def test_memory_failure_never_breaks_spawn(tmp_path, monkeypatch):
    ws = _inheritance_workspace(tmp_path, monkeypatch)
    perms = KitToolPermissions(workspace=ws, mode="default")

    def _boom(*a, **kw):
        raise RuntimeError("memory subsystem down")

    monkeypatch.setattr(
        "delfin.agent.project_memory.load_project_memory", _boom)
    monkeypatch.setattr(
        "delfin.agent.memory_store.list_typed_memories", _boom)
    client = _fake_client("survived without memories")
    res = SA.run_subagent(
        subagent_type="explore",
        description="resilience check",
        prompt="Long enough prompt for the child task here.",
        parent_client=client,
        parent_perms=perms,
    )
    assert res.error == ""
    assert "survived without memories" in res.final_text
    sys_prompt = _child_system_prompt(client)
    assert "Project rules (inherited" not in sys_prompt
    assert "Standing memory (inherited)" not in sys_prompt


def test_inherit_context_skipped_without_workspace():
    client = _fake_client("no workspace, still fine")
    res = SA.run_subagent(
        subagent_type="explore",
        description="perms-less run",
        prompt="Explore something without any permissions attached.",
        parent_client=client,
        parent_perms=None,
    )
    assert res.error == ""
    assert "Project rules (inherited" not in _child_system_prompt(client)


def test_runner_is_attached_by_openaiclient_set_permissions():
    """OpenAIClient.set_permissions should bind perms.subagent_runner."""
    from delfin.agent.api_client import OpenAIClient
    # We don't actually need a working API call, just instance + set_perms.
    with tempfile.TemporaryDirectory() as d:
        perms = KitToolPermissions(workspace=Path(d), mode="default")
        # Construct a client without making any HTTP calls
        client = OpenAIClient.__new__(OpenAIClient)
        client._permissions = None
        client.model = "x"
        client.set_permissions(perms)
        assert perms.subagent_runner is not None
        assert callable(perms.subagent_runner)
