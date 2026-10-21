"""One preset was a specialist; the rest were shapes of "look at this".

A session could delegate research (explore), design (plan), a diff read
(code-reviewer), a chemistry audit (chemistry-reviewer) or unrestricted
construction (general-purpose). Everything else -- does this change widen
what the agent may do, does it break anything, which functional for this
property, what do these finished calculations actually say, where are the
turns going -- was general-purpose with a paragraph of instructions in
the prompt, which is a specialist written from scratch on every call and
forgotten afterwards.

Five more presets, as markdown in the pack rather than code: a preset is
data, and the discovery path (`*_subagent.md`) already reads it.

What is asserted here is what makes a preset a specialist rather than a
label: it narrows its tools, it carries a method, it says what it must
NOT claim, and -- because `_narrow_allowed_tools` only narrows -- it can
never reach past the session that spawned it.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import subagents as SA

_NEW = ("security-reviewer", "verifier", "method-researcher",
        "data-extractor", "friction-analyst")

_PACK = (pathlib.Path(SA.__file__).resolve().parent / "pack" / "agents")


@pytest.fixture(autouse=True)
def _fresh():
    SA.reload_subagent_presets()


def _preset(name):
    return SA.SUBAGENT_PRESETS[name]


# ---------------------------------------------------------------------------
# They exist, and they are reachable
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", _NEW)
def test_the_preset_is_discovered(name):
    assert name in SA.SUBAGENT_PRESETS


@pytest.mark.parametrize("name", _NEW)
def test_the_model_can_name_it(name):
    """The subagent tool's enum is built from this list, so a preset that
    is not in it cannot be asked for."""
    assert name in SA.subagent_type_names()


@pytest.mark.parametrize("name", _NEW)
def test_it_is_a_file_in_the_pack_not_code(name):
    assert (_PACK / f"{name}_subagent.md").is_file()


def test_the_older_presets_still_exist():
    for name in ("explore", "plan", "code-reviewer", "general-purpose",
                 "chemistry-reviewer"):
        assert name in SA.SUBAGENT_PRESETS


# ---------------------------------------------------------------------------
# What makes it a specialist
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", _NEW)
def test_it_narrows_its_tools(name):
    """A preset with no list is unrestricted, which is general-purpose's
    job and nobody else's."""
    tools = getattr(_preset(name), "tools", ()) or ()
    assert tools, f"{name} restricts nothing"
    assert len(tools) < 60


@pytest.mark.parametrize("name", _NEW)
def test_it_can_reach_the_session_running_it(name):
    tools = getattr(_preset(name), "tools", ()) or ()
    assert "subagent_message" in tools, name


@pytest.mark.parametrize("name", _NEW)
def test_it_carries_a_method_and_not_only_a_label(name):
    """The system prompt is what separates a specialist from a name: it
    says how to do the job. A one-line prompt is a label."""
    body = _preset(name).system_prompt
    assert len(body) > 400, f"{name} has no method, only a description"
    assert "You are a" in body


@pytest.mark.parametrize("name", _NEW)
def test_it_says_what_it_must_not_claim(name):
    """Every report-producing preset in this pack states a limit on its
    own conclusions. That is the half that keeps a delegate's report from
    reading as more than it measured."""
    body = _preset(name).system_prompt.lower()
    assert any(word in body for word in
               ("never", "do not", "not do", "rather than")), name


@pytest.mark.parametrize("name", _NEW)
def test_it_has_a_description_a_session_can_choose_by(name):
    desc = (_preset(name).description or "").strip()
    assert len(desc) > 40, name
    assert "use" in desc.lower() or "run" in desc.lower(), name


# ---------------------------------------------------------------------------
# Only one of them may act
# ---------------------------------------------------------------------------

_READ_ONLY_NEW = ("security-reviewer", "method-researcher",
                  "data-extractor", "friction-analyst")


@pytest.mark.parametrize("name", _READ_ONLY_NEW)
def test_a_reviewing_specialist_cannot_write(name):
    tools = set(getattr(_preset(name), "tools", ()) or ())
    for writer in ("write_file", "edit_file", "multi_edit", "bash",
                   "apply_patch", "notebook_edit"):
        assert writer not in tools, f"{name} may {writer}"
    assert _preset(name).mode == "plan", name


def test_the_verifier_may_run_commands_but_not_edit():
    """It exists to measure, and measuring a test suite means running it.
    Editing would let it make a failure go away instead of reporting it."""
    tools = set(getattr(_preset("verifier"), "tools", ()) or ())
    assert {"bash", "run_tests"} <= tools
    for writer in ("write_file", "edit_file", "multi_edit", "apply_patch"):
        assert writer not in tools


def test_only_the_verifier_among_the_new_ones_can_act():
    acting = [n for n in _NEW
              if "bash" in (getattr(_preset(n), "tools", ()) or ())]
    assert acting == ["verifier"]


# ---------------------------------------------------------------------------
# A preset cannot widen its session
# ---------------------------------------------------------------------------

def test_a_preset_cannot_hand_a_child_a_tool_the_session_denied():
    """The direction that matters: presets are files on disk, so one able
    to widen its parent would turn dropping a file into
    ~/.delfin/subagents into a privilege escalation."""
    from delfin.agent.api_client import KitToolPermissions

    import tempfile
    ws = pathlib.Path(tempfile.mkdtemp(prefix="preset-"))
    parent = KitToolPermissions(workspace=ws, mode="default",
                                session_allowed_tools=frozenset({"read_file"}))
    narrowed = SA._narrow_allowed_tools(
        parent, getattr(_preset("verifier"), "tools", ()) or ())
    assert narrowed is not None
    assert "bash" not in narrowed, (
        "a preset reached past the session's own allow list")
    assert "read_file" in narrowed


# ---------------------------------------------------------------------------
# The session is told they exist
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", _NEW)
def test_the_role_prompt_names_it(name):
    """A preset the prompt does not mention is one the model will not
    reach for. Asserted on the shipped file, which is what is injected."""
    text = (_PACK / "solo_agent.md").read_text(encoding="utf-8")
    assert f"`{name}`" in text, name
