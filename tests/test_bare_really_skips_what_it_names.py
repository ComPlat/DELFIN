"""``--bare`` skips what it names, and names what it still does not skip.

The flag promised to skip discovery and skipped one of four things. The
other three ran on regardless: hook definitions loaded at every tool call,
skill playbooks read while the tool surface was assembled, the memory
stores read while the system prompt was built. A flag that parses, prints
a reassuring line and changes almost nothing is the silent non-delivery
this area exists to close.

So every test here asserts against the thing that does the work — the MCP
registry, the hook config that reaches ``run_hooks``, the executor's own
answer about skills, the built system prompt — and never against the
namespace or the switch itself. A test that reads back the flag it set
would pass on the code as it was.

Everything the tests need is supplied here: the settings file through
``_user_settings_path``, the home directory through ``Path.home``, the
prompt pack as a temporary tree. Nothing is asserted about this machine.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from pathlib import Path

import pytest

from delfin.agent import api_client
from delfin.agent import cli as agent_cli
from delfin.agent import hooks as hooks_mod
from delfin.agent.prompt_loader import PromptLoader

PROBE_SKILL = "probe-skill"
PROBE_HOOK_COMMAND = "probe-hook-command"
PROBE_MEMO = "the probe memo line"


# ---------------------------------------------------------------------------
# What a session is given: a workspace, a settings file, a memory store and a
# prompt pack — all temporary, all supplied by the test.
# ---------------------------------------------------------------------------

def _args(**over) -> argparse.Namespace:
    base = dict(backend="", provider="", model="", mode="", cwd="",
                effort="", permission_mode="", settings_defaults=False)
    base.update(over)
    return argparse.Namespace(**base)


@pytest.fixture
def fake_home(tmp_path, monkeypatch):
    """A home directory of our own, so no real store can be read or hit."""
    home = tmp_path / "home"
    (home / ".delfin").mkdir(parents=True)
    monkeypatch.setattr(pathlib.Path, "home", staticmethod(lambda: home))
    return home


@pytest.fixture
def workspace(tmp_path):
    """A project directory carrying one project-scoped skill."""
    ws = tmp_path / "ws"
    skills_dir = ws / ".delfin" / "skills"
    skills_dir.mkdir(parents=True)
    (skills_dir / f"{PROBE_SKILL}.md").write_text(
        f"---\nname: {PROBE_SKILL}\ndescription: a probe playbook\n---\n\n"
        "# Probe\n\n1. step one of the probe playbook\n",
        encoding="utf-8")
    return ws


@pytest.fixture
def user_hooks(tmp_path, monkeypatch):
    """One PreToolUse hook in the user's own settings file.

    The user-level file is the case that matters: ``load_hooks`` reads it
    whatever workspace it is given, so switching hooks off cannot be done
    by narrowing the workspace.
    """
    settings = tmp_path / "user-settings.json"
    settings.write_text(json.dumps({"hooks": {"PreToolUse": [
        {"matcher": "*", "hooks": [
            {"type": "command", "command": PROBE_HOOK_COMMAND},
        ]},
    ]}}), encoding="utf-8")
    monkeypatch.setattr(hooks_mod, "_user_settings_path", lambda: settings)
    return settings


@pytest.fixture
def hook_spies(monkeypatch):
    """Record what the hook path loaded and fired, without running anything.

    ``run_hooks`` spawns shell commands; a test that let it would be
    asserting that this machine has a shell.
    """
    fired: list[tuple] = []
    loaded: list[tuple] = []
    real_load = hooks_mod.load_hooks

    def _run(event, cfg, **kwargs):
        fired.append((event, cfg))
        return []

    def _load(*a, **kw):
        loaded.append(a)
        return real_load(*a, **kw)

    monkeypatch.setattr(hooks_mod, "run_hooks", _run)
    monkeypatch.setattr(hooks_mod, "load_hooks", _load)
    return fired, loaded


@pytest.fixture
def memory_store(fake_home, workspace):
    """A project memory store where the loader really looks for one."""
    from delfin.agent.memory_store import _delfin_memory_dir

    mem = _delfin_memory_dir(workspace)
    mem.mkdir(parents=True, exist_ok=True)
    (mem / "MEMORY.md").write_text(
        f"# Project Memory\n- {PROBE_MEMO}.\n", encoding="utf-8")
    return mem


@pytest.fixture
def prompt_tree(tmp_path):
    """A minimal prompt pack, so the prompt build reads no shipped file."""
    tree = tmp_path / "pack_tree"
    (tree / "pack" / "agents").mkdir(parents=True)
    (tree / "pack" / "shared").mkdir(parents=True)
    (tree / "pack_lite" / "modes").mkdir(parents=True)
    (tree / "pack" / "agents" / "solo_agent.md").write_text(
        "# Solo Agent\nYou are the solo agent.", encoding="utf-8")
    return tree


class _StandInClient:
    """Carries the permissions object; that is all these paths need of it."""

    def __init__(self, permissions):
        self._permissions = permissions
        self.model = "kit.qwen3-coder-480b"


@pytest.fixture
def build_engine(monkeypatch, workspace, prompt_tree):
    """`_build_engine`, against a stand-in that carries the real objects.

    The switches live on the permissions object and on the prompt loader,
    so the stand-in holds real ones — a fake attribute bag would let the
    wiring pass while the loaders it is supposed to reach kept running.
    """
    import delfin.agent.engine as engine_mod

    class _StandInEngine:
        def __init__(self, **kwargs):
            self.client = _StandInClient(
                api_client.KitToolPermissions(workspace=workspace))
            self.loader = PromptLoader(prompt_tree)
            self.loader.workspace_root = workspace
            self.provider = kwargs.get("provider") or "kit"

        @property
        def kit_permissions(self):
            return getattr(self.client, "_permissions", None)

    monkeypatch.setattr(engine_mod, "AgentEngine", _StandInEngine)

    def _build(**over):
        return agent_cli._build_engine(_args(cwd=str(workspace), **over))

    return _build


@pytest.fixture
def clean_registries():
    """MCP registries are process-global and cached per workspace."""
    from delfin.agent import mcp_client

    mcp_client.reset_registry()
    yield mcp_client
    mcp_client.reset_registry()


def _prompt(loader: PromptLoader) -> str:
    return loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo", mode_description="solo",
        route=["solo_agent"], role_index=0)


def _pre_tool_hooks(perms) -> None:
    """Drive the hook point the executor runs before every tool call."""
    api_client._DocToolExecutor()._run_pre_tool_hooks(
        "bash", {"command": "true"}, perms)


# ---------------------------------------------------------------------------
# 1. MCP servers — the registry the tool surface reads
# ---------------------------------------------------------------------------

def test_bare_leaves_no_mcp_server_for_the_tool_surface(clean_registries,
                                                        workspace):
    """Asserted on the registry: no server configured, nothing discovered.

    This is the half that runs before the engine exists, because the
    registry is cached per workspace and loads its configuration once.
    """
    mcp = clean_registries
    assert mcp.get_registry(workspace).servers, (
        "the built-in servers should be configured before --bare runs")

    assert agent_cli._apply_bare(workspace) is True

    assert mcp.get_registry(workspace).servers == {}
    assert mcp.get_registry(workspace).discover_all() == []
    # The key a backend without a permissions object resolves to.
    assert mcp.get_registry(None).servers == {}


def test_without_bare_the_mcp_servers_stay_configured(clean_registries,
                                                      workspace):
    assert clean_registries.get_registry(workspace).servers != {}


# ---------------------------------------------------------------------------
# 2. Hooks — hooks.load_hooks, reached at every tool call
# ---------------------------------------------------------------------------

def test_bare_means_no_hook_definition_is_loaded_or_fired(
        build_engine, user_hooks, hook_spies):
    """Asserted on the loader and on what reached ``run_hooks``.

    The user's settings file holds a PreToolUse hook. After ``--bare`` the
    tool call must load nothing and fire nothing — not load it and decline
    to run it, because the file is executable configuration and reading it
    is the step that has to be skipped.
    """
    fired, loaded = hook_spies
    engine = build_engine(bare=True)

    _pre_tool_hooks(engine.kit_permissions)

    assert loaded == [], "load_hooks must not be reached under --bare"
    assert fired == [], "no hook may fire under --bare"


def test_without_bare_the_same_hook_is_loaded_and_fired(
        build_engine, user_hooks, hook_spies):
    """The control: the identical setup, minus the flag."""
    fired, loaded = hook_spies
    engine = build_engine()

    _pre_tool_hooks(engine.kit_permissions)

    assert loaded, "the hook loader is reached without --bare"
    assert [event for event, _ in fired] == ["PreToolUse"]
    cfg = fired[0][1]
    assert [c.command for c in cfg.for_event("PreToolUse")] == [
        PROBE_HOOK_COMMAND]


def test_a_sub_agent_inherits_the_hook_switch(build_engine, user_hooks,
                                              hook_spies):
    """A child must not be a route around what the parent was told to skip.

    ``_derive_perms`` clones with ``dataclasses.replace``, which copies
    field values — so the switch travels with the child by construction.
    """
    from delfin.agent import subagents

    fired, loaded = hook_spies
    engine = build_engine(bare=True)
    child = subagents._derive_perms(engine.kit_permissions, "default")

    _pre_tool_hooks(child)

    assert loaded == [] and fired == []


# ---------------------------------------------------------------------------
# 3. Skills — skills.discover_skills, during tool-surface assembly
# ---------------------------------------------------------------------------

def test_bare_means_no_skill_is_read_from_disk(build_engine, workspace,
                                               monkeypatch):
    """Asserted on the executor's answer and on the surface listing.

    A skill sits in the workspace. Under ``--bare`` neither the list built
    into the tool description nor the executor behind that tool may find
    it, and the playbook body must not be returned to the model.
    """
    from delfin.agent import skills as skills_mod

    reached: list = []
    real = skills_mod.discover_skills
    monkeypatch.setattr(
        skills_mod, "discover_skills",
        lambda *a, **kw: (reached.append(a) or real(*a, **kw)))

    engine = build_engine(bare=True)
    perms = engine.kit_permissions

    assert api_client._session_skills(perms) == []
    payload = json.loads(
        api_client._DocToolExecutor()._execute_skill({"name": PROBE_SKILL},
                                                     perms))
    assert payload.get("available") == []
    assert "step one of the probe playbook" not in json.dumps(payload)
    assert reached == [], "discover_skills must not be reached under --bare"


def test_without_bare_the_same_skill_is_found_and_returned(build_engine,
                                                           workspace):
    """The control: the identical workspace, minus the flag."""
    engine = build_engine()
    perms = engine.kit_permissions

    assert PROBE_SKILL in [s.name for s in api_client._session_skills(perms)]
    payload = json.loads(
        api_client._DocToolExecutor()._execute_skill({"name": PROBE_SKILL},
                                                     perms))
    assert PROBE_SKILL in json.dumps(payload)
    assert "step one of the probe playbook" in json.dumps(payload)


def test_a_sub_agent_inherits_the_skill_switch(build_engine, workspace):
    from delfin.agent import subagents

    engine = build_engine(bare=True)
    child = subagents._derive_perms(engine.kit_permissions, "default")

    assert api_client._session_skills(child) == []


# ---------------------------------------------------------------------------
# 4. Project memory — PromptLoader._load_external_memory_context
# ---------------------------------------------------------------------------

def test_bare_keeps_the_memory_stores_out_of_the_system_prompt(
        build_engine, memory_store):
    """Asserted on the built prompt, which is where the memories landed."""
    engine = build_engine(bare=True)

    prompt = _prompt(engine.loader)

    assert "External Memory" not in prompt
    assert PROBE_MEMO not in prompt
    assert engine.loader._load_external_memory_context() == ""


def test_without_bare_the_same_memory_reaches_the_prompt(build_engine,
                                                         memory_store):
    """The control: the identical store, minus the flag."""
    engine = build_engine()

    prompt = _prompt(engine.loader)

    assert "External Memory" in prompt
    assert PROBE_MEMO in prompt


# ---------------------------------------------------------------------------
# A bound on one session is not a change to every later one
# ---------------------------------------------------------------------------

def _snapshot(*roots: Path) -> dict[str, bytes]:
    out: dict[str, bytes] = {}
    for root in roots:
        for p in sorted(Path(root).rglob("*")):
            if p.is_file():
                out[str(p)] = p.read_bytes()
    return out


def test_bare_writes_nothing_to_any_settings_file(
        build_engine, clean_registries, workspace, fake_home, user_hooks,
        memory_store, prompt_tree):
    """Snapshot of the home directory and the workspace, before and after.

    ``--bare`` bounds THIS session. A settings file belongs to every later
    one, so the flag that switches discovery off must leave every file on
    disk exactly as it found it — including the settings file whose hooks
    it just declined to read.
    """
    roots = (fake_home, workspace, prompt_tree, user_hooks.parent)
    before = _snapshot(*roots)

    agent_cli._apply_bare(workspace)
    engine = build_engine(bare=True)
    _pre_tool_hooks(engine.kit_permissions)
    api_client._session_skills(engine.kit_permissions)
    _prompt(engine.loader)

    assert _snapshot(*roots) == before


# ---------------------------------------------------------------------------
# What the startup line and the help text may claim
# ---------------------------------------------------------------------------

def test_the_startup_line_names_each_switch_that_actually_took(
        build_engine, workspace):
    """Read back off the objects, so a name on the line means it arrived."""
    engine = build_engine(bare=True)
    skipped, unreached = agent_cli._bare_coverage(
        _args(bare=True, bare_mcp_skipped=True), engine)

    assert unreached == []
    assert skipped == list(agent_cli._BARE_SKIPS)

    line = "\n".join(agent_cli._bounding_notices(
        _args(bare=True, bare_mcp_skipped=True), engine))
    for name in agent_cli._BARE_SKIPS:
        assert name in line
    assert "skipped" in line
    # And the one thing it still does not cover is named, not left to be
    # discovered from a hook that fired anyway.
    assert agent_cli._BARE_NOT_SKIPPED in line


def test_a_backend_that_carries_none_of_the_objects_says_so(build_engine):
    """The can't-deliver path: nothing took, so nothing may be claimed."""
    engine = type("E", (), {"client": type("C", (), {"model": "m"})(),
                            "provider": "kit"})()
    notes = agent_cli._bounding_notices(
        _args(bare=True, bare_mcp_skipped=False), engine)

    assert any("REQUESTED but nothing was skipped" in n for n in notes)
    assert not any("skipped —" in n and "MCP servers skipped" in n
                   for n in notes)


def test_a_half_applied_bare_names_the_half_that_did_not_take(build_engine):
    """A backend with no permissions object still gets memory and MCP."""
    engine = build_engine(bare=True)
    engine.client._permissions = None

    line = "\n".join(agent_cli._bounding_notices(
        _args(bare=True, bare_mcp_skipped=True), engine))

    assert "MCP servers" in line and "project memory" in line
    assert "not reachable on this backend" in line
    assert "tool hooks" in line and "skills" in line


def test_without_bare_the_startup_line_says_nothing_about_it(build_engine):
    engine = build_engine()
    assert agent_cli._bounding_notices(_args(), engine) == []


def test_the_help_text_and_the_startup_line_come_from_one_source():
    """One source for both, so the promise cannot widen in one place."""
    import io
    import re

    buf = io.StringIO()
    parser = agent_cli.build_parser()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            action.choices["chat"].print_help(buf)
            break
    # argparse re-wraps help text, so the sentence is compared without the
    # line breaks it inserted.
    text = re.sub(r"\s+", " ", buf.getvalue())

    assert "--bare" in text
    for name in agent_cli._BARE_SKIPS:
        assert name in text
    assert agent_cli._BARE_NOT_SKIPPED in text
