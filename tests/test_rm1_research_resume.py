"""A research session saved elsewhere reopens as Research, read-only.

Phase-3 red test (RM1 wave 17). The dashboard restore path has three
parts the phase-2 test pins only in source: the mode dropdown keeps the
saved "research" value, the engine built on restore works in the fixed
agent scratch folder (ctx.agent_dir, launch-independent, not the
session's recorded workspace and not the launch directory), and the
tool list it builds carries no tool that can change a file or run a
shell command.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("ipywidgets")

_SID = "a1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6"
_MODEL = "kit.deepseek-v4-flash"

_WRITE_TOOLS = ("Write", "Edit", "MultiEdit", "NotebookEdit", "Bash")


@pytest.fixture
def env(tmp_path, monkeypatch):
    """A throwaway home, a saved research session whose recorded
    workspace names ANOTHER directory, and a dashboard tab whose launch
    context points at yet another agent_dir."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test")

    from delfin import user_settings
    from delfin.agent import scheduler as S
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "kit", "model": _MODEL,
        "permission_profile": "all_free", "probe_endpoint": False}}))
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))

    # The saved session was recorded in a workspace the launch context
    # does not know (simulating a resume from another launch directory).
    other_workspace = tmp_path / "somewhere-else"
    other_workspace.mkdir()

    from delfin.agent.session_store import save_session, release_session_lock
    save_session(
        _SID,
        mode="research",
        chat_messages=[{"role": "user", "content": "literature question"}],
        workspace=str(other_workspace),
        provider="kit",
        model=_MODEL,
    )
    release_session_lock(_SID)

    # The tab's scratch folder is its own place, not the launch cwd.
    agent_dir = tmp_path / "scratch"
    agent_dir.mkdir()

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    from delfin.dashboard.context import DashboardContext
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
        agent_dir=agent_dir,
    )
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "ws-launch")
    ctx.initial_session_id = _SID
    ctx.sessions_managed = True
    ctx.presence_key = "k1"

    from delfin.dashboard import tab_agent
    tab, refs = tab_agent.create_tab(ctx)
    try:
        yield SimpleNamespace(tab=tab, refs=refs, agent_dir=agent_dir,
                              other_workspace=other_workspace)
    finally:
        refs["shutdown"](save=False)
        release_session_lock(_SID)


def _dropdown(tab, description):
    import ipywidgets as widgets
    stack = [tab]
    while stack:
        node = stack.pop()
        if (isinstance(node, widgets.Dropdown)
                and node.description == description):
            return node
        stack.extend(getattr(node, "children", ()) or ())
    raise AssertionError(f"no {description!r} dropdown in the tab")


def _find_html(tab, prefix):
    import ipywidgets as widgets
    stack = [tab]
    while stack:
        node = stack.pop()
        if isinstance(node, widgets.HTML) and str(
                node.value or "").startswith(prefix):
            return node
        stack.extend(getattr(node, "children", ()) or ())
    return None


def test_research_session_restores_as_research(env):
    """The saved mode survives the restore: research stays research."""
    dd = _dropdown(env.tab, "Mode:")
    assert dd.value == "research", (
        f"a saved research session reopened as {dd.value!r}")


def test_restored_research_engine_runs_in_the_fixed_scratch_folder(env):
    """Like Dashboard, the restored research engine works in the fixed
    agent scratch folder — reachable from any launch, never pinned to
    the session's recorded workspace."""
    engine = env.refs["state"].get("engine")
    assert engine is not None, "restore built no engine"
    kp = getattr(engine, "kit_permissions", None)
    ws = str(getattr(kp, "workspace", "")
             or getattr(engine, "repo_dir", ""))
    assert ws == str(env.agent_dir), (
        f"research engine workspace {ws!r} is not the fixed scratch "
        f"folder {str(env.agent_dir)!r}")


def test_restored_research_engine_has_no_write_tool(env):
    """The restored research engine's tool list carries no tool that can
    change a file or run a shell command."""
    engine = env.refs["state"].get("engine")
    assert engine is not None, "restore built no engine"
    tools = getattr(engine, "cli_tools", None) or getattr(
        engine, "allowed_tools", None) or []
    names = {str(getattr(t, "name", t)) for t in tools}
    leaked = sorted(n for n in names
                    if any(n.startswith(w) or n == w for w in _WRITE_TOOLS))
    assert not leaked, (
        f"restored research engine offers write tools {leaked}")
    assert getattr(engine, "kit_permissions", None) is not None and (
        getattr(engine.kit_permissions, "agent_role", "")
        == "research_agent"), (
        "restored research engine is not role-locked to research_agent")
