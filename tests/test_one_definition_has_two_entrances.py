"""A specialist is reached two ways and defined once.

The generalist delegates to `workflow-operator`; a person who already
knows the task is a calculation chooses it directly with
`--mode workflow`. Two entry points, because they serve two different
people: a delegation costs an extra model roundtrip, and making it the
only route makes every direct user pay for a decision they had already
taken.

ONE DEFINITION. A mode used to mean a fifth hand-maintained role: a
prompt file in the pack, a model preference, a token budget, a steering
entry and an execution allow-list, in five separate tables. That shape
is how this repository's lists drifted before -- the subagent enum
offered 13 roles while the prompt's table named 7, and the six without a
cue were never chosen once in 79,056 tool calls.

So the mode's route names the PRESET, and everything about the
specialist keeps coming from that one file:

  * the prompt, because `load_role_prompt` falls back to the preset body
    where no `agents/<role>.md` exists;
  * the tool surface, because `_role_allowlist` derives the role's
    execution allow-list from the preset's `tools:` line.

AND THE SECURITY PROPERTY THIS FILE EXISTS FOR. A role absent from
`_ROLE_EXEC_ALLOWLIST` has NO role-level restriction -- `allow is None`
means "not refused" at both call sites. So without the derivation the
same specialist would be LESS restricted when chosen directly than when
delegated, and the second entry point would be a way around the first.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as ac
from delfin.agent import subagents as sa
from delfin.agent.prompt_loader import PromptLoader

ROLE = "workflow-operator"
MODE = "workflow"


@pytest.fixture(scope="module")
def loader():
    return PromptLoader()


@pytest.fixture(autouse=True)
def _fresh_registry():
    sa.reload_subagent_presets()
    yield
    sa.reload_subagent_presets()


# --- two entrances ------------------------------------------------------

def test_the_mode_exists_and_routes_to_the_preset(loader):
    modes = [m.get("id") for m in loader.load_lite_manifest().get("modes", [])]
    assert MODE in modes
    assert loader.load_mode(MODE)["route"] == [ROLE]


def test_the_mode_carries_a_description_a_person_can_read(loader):
    text = loader.load_mode(MODE)["description"]
    assert text.strip()
    assert "resolve_spec" in text


def test_the_specialist_is_still_delegable():
    """The other entrance must not have been replaced by this one."""
    assert ROLE in sa.subagent_type_names()


# --- one definition -----------------------------------------------------

def test_the_role_prompt_is_the_presets_own_body(loader):
    """Not a copy in the pack: the same file either way."""
    body = loader.load_role_prompt(ROLE)
    assert "resolve_spec" in body
    assert "ask_user_question" in body
    # the frontmatter is not part of a role prompt
    assert not body.lstrip().startswith("---")
    assert "write_globs:" not in body


def test_a_role_with_its_own_file_is_unaffected(loader):
    """The fallback is a fallback: a pack file still wins."""
    body = loader.load_role_prompt("solo_agent")
    assert body.strip()
    assert "resolve_spec" not in body, (
        "solo_agent must still get its own prompt, not a preset's")


def test_role_exists_answers_for_both_sources(loader):
    """One question, two sources. The directory stopped being the
    registry when a route could name a preset, and a caller still asking
    it refused a role the engine runs -- which is what
    `prompt-size --role workflow-operator` reported."""
    assert loader.role_exists("solo_agent")       # a pack file
    assert loader.role_exists(ROLE)               # a preset
    assert not loader.role_exists("no-such-role")
    assert not loader.role_exists("")


# --- the security property -----------------------------------------------

def test_the_direct_path_carries_the_presets_restriction():
    allow = ac._role_allowlist(ROLE)
    assert allow is not None, (
        "a preset-backed role with no restriction would make the direct "
        "entrance a way around the delegated one")
    expected = {str(t).rsplit("__", 1)[-1]
                for t in (sa.SUBAGENT_PRESETS[ROLE].tools or ())}
    assert allow == frozenset(expected)


@pytest.mark.parametrize("tool", [
    "resolve_spec", "run_application", "run_diagnostics",
    "ask_user_question", "mcp__delfin-tools__resolve_spec",
])
def test_what_the_specialist_needs_is_permitted(tool):
    assert ac._role_allowlist(ROLE)
    base = tool.rsplit("__", 1)[-1]
    assert base in ac._role_allowlist(ROLE)


@pytest.mark.parametrize("tool", ["git_commit", "apply_patch", "web_fetch"])
def test_what_it_has_no_business_with_is_refused(tool):
    """Not in the preset's tools line, so not in the role's either."""
    assert tool.rsplit("__", 1)[-1] not in ac._role_allowlist(ROLE)


def test_a_role_that_is_not_a_preset_keeps_no_restriction():
    """No behaviour change for the roles that had none: `allow is None`
    still means no role-level restriction."""
    assert ac._role_allowlist("solo_agent") is None
    assert ac._role_allowlist("") is None
    assert ac._role_allowlist("no-such-role") is None


def test_the_roles_that_had_a_restriction_still_have_theirs():
    for role in ("dashboard_agent", "office_agent"):
        allow = ac._role_allowlist(role)
        assert allow is not None and allow is ac._ROLE_EXEC_ALLOWLIST[role], (
            f"{role}: the dict must still answer first")


def test_a_preset_with_no_tools_line_adds_no_restriction():
    """SubagentPreset.tools empty means "no preset-level restriction",
    and it must not turn into "no tools" here -- the opposite default
    would let a preset without a frontmatter list silently become the
    most locked-down role in the system."""
    class _Bare:
        tools = ()

    saved = sa.SUBAGENT_PRESETS.get("_bare_probe")
    sa.SUBAGENT_PRESETS["_bare_probe"] = _Bare()
    try:
        assert ac._role_allowlist("_bare_probe") is None
    finally:
        if saved is None:
            sa.SUBAGENT_PRESETS.pop("_bare_probe", None)
        else:
            sa.SUBAGENT_PRESETS["_bare_probe"] = saved


# --- it composes ---------------------------------------------------------

def test_the_whole_prompt_composes_for_the_new_mode(loader):
    text = loader.build_system_prompt(ROLE, MODE)
    assert len(text) > 10_000, "the composed prompt looks empty"
    assert "resolve_spec" in text
