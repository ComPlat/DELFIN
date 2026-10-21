"""Adversarial tests for RM1 phase 2 (reviewer s21).

Each test attacks a claim the builder's phase-2 commit makes, through a
real caller, not a source literal.
"""

from __future__ import annotations

import json

import pytest


# ---------------------------------------------------------------------------
# 1. The description says "no shell commands" — does the API path agree?
# ---------------------------------------------------------------------------

def test_research_role_denies_shell_on_the_api_path():
    """The dropdown description and the README both say research runs
    "no shell commands". On the CLI backend the _cli_tools list excludes
    Bash, but on the KIT/OpenAI backends the executor authorises by
    _tool_denied_for_role — and that denylist for research_agent names
    only the write tools. If bash passes it, the description is a promise
    the executor does not keep (api_client.py:7215-7229)."""
    from delfin.agent.api_client import (
        ToolSurfaceContext, advertisable_tools, _tool_denied_for_role)

    assert _tool_denied_for_role("research_agent", "bash"), (
        "research_agent may execute bash on the API path while the mode "
        "description promises 'no shell commands' — the description lies "
        "or the denylist must name bash")
    ctx = ToolSurfaceContext(role="research_agent")
    catalogue = [{"function": {"name": "bash", "description": "x",
                               "parameters": {}}}]
    assert advertisable_tools(catalogue, ctx) == [], (
        "bash is advertised to a research session; the user was told "
        "there are no shell commands")


# ---------------------------------------------------------------------------
# 2. A saved research session round-trips through the real store
# ---------------------------------------------------------------------------

def test_saved_research_session_round_trips_and_lists_by_mode(
        tmp_path, monkeypatch):
    """save_session writes mode='research'; load_session reads it back and
    list_sessions (which applies _LEGACY_MODE_MAP) still reports research,
    so the dropdown filter s.get('mode') == 'research' finds it."""
    from delfin.agent import session_store as store

    monkeypatch.setattr(store.Path, "home", lambda: tmp_path, raising=False)
    store.save_session("rm1-adv-1", mode="research",
                       chat_messages=[{"role": "user", "content": "hi"}])
    loaded = store.load_session("rm1-adv-1")
    assert loaded is not None
    assert loaded["mode"] == "research"
    listed = [s for s in store.list_sessions(limit=10)
              if s.get("mode") == "research"]
    assert any(s["session_id"] == "rm1-adv-1" for s in listed), (
        "a saved research session is not listed under mode research; the "
        "session list filters by mode for launch-independent modes")


# ---------------------------------------------------------------------------
# 3. The engine-rebuild decision on switching modes (workspace groups)
# ---------------------------------------------------------------------------

def test_every_research_switch_rebuilds_the_engine():
    """Research sits in its OWN workspace group (tab_agent.py:2519-2523):
    every switch into or out of research reports a workspace change, so
    the next send builds a fresh engine. That is stronger than sharing
    the dashboard group — research never inherits a Code engine (with
    its Write/Edit tool list) and never hands one over. The shared
    scratch folder would have allowed research<->dashboard to keep the
    engine, but the safety property the rebuild guarantees is the tool
    list and role stamp matching the mode, not the folder.

    (Original decision: this test first asserted research shares the
    dashboard group and must NOT rebuild; 0f5475e5 deliberately made
    research its own group so every research switch rebuilds. Rewritten
    along the intent — rebuilds are only wasted when the tool list and
    workspace would be identical, and dashboard's engine carries Write/
    Edit while research's does not, so rebuilding is required here.)"""
    from delfin.dashboard.tab_agent import _mode_workspace_differs
    for other in ("dashboard", "solo", "office"):
        assert _mode_workspace_differs("research", other), (
            f"research -> {other} did not report a workspace change; a "
            "kept engine would run the other mode on the research "
            "tool list or the research agent on a write-enabled list")
        assert _mode_workspace_differs(other, "research"), (
            f"{other} -> research did not report a workspace change; "
            "the research engine must be rebuilt read-only, not "
            "inherited from the previous mode")


def test_research_and_solo_differ_workspace_rebuild():
    """agent scratch folder vs the launch directory: switching MUST report
    a workspace change, so the next send builds an engine whose client
    matches the new folder and the new mode's tool list."""
    from delfin.dashboard.tab_agent import _mode_workspace_differs
    assert _mode_workspace_differs("research", "solo"), (
        "research binds the agent scratch folder, solo the launch "
        "directory; keeping the engine across this switch leaves the "
        "research client's read-only tool list in charge of a Code "
        "conversation")
    assert _mode_workspace_differs("solo", "research")
