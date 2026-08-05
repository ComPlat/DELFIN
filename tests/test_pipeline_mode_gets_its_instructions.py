"""Pipeline was never a mode; its instructions still have to arrive.

The file modes/pipeline.md said: call ``get_guide`` first,
``validate_spec`` before every run, work only through the
``delfin-tools`` server, never hand-roll pipeline logic. For a long time
none of it reached the model at all -- the gate that injects a mode
description sat after the branch that pipeline took.

Delivering it fixed the symptom and exposed the cause: pipeline routed
to ``solo_agent``, advertised the same tools, and differed from Code by
exactly one page of instructions. That is a procedure, which is what a
skill is. So the mode is retired and the page is the ``pipeline-build``
skill: reachable from Code without switching into a mode you can then be
stuck in, and paid for on the turn it is invoked rather than on every
turn of a session someone left the dropdown on.

These tests hold the guarantee across that move. The instructions must
still arrive, the old name must still open, and the retirement must not
have left a mode behind that no longer exists.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import skills as S
from delfin.agent.engine import _migrate_mode
from delfin.agent.prompt_loader import PromptLoader

_PACK = pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"


@pytest.fixture
def loader():
    return PromptLoader(_PACK)


# ---------------------------------------------------------------------------
# The instructions survived the move
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("marker", [
    "get_guide", "validate_spec", "delfin-tools", "save_application",
    "run_diagnostics",
])
def test_the_pipeline_instructions_still_exist(marker):
    sk = S.get_skill("pipeline-build")
    assert sk is not None, "the pipeline procedure was deleted, not moved"
    assert marker in sk.body, marker


def test_the_rule_the_page_exists_to_state_is_present():
    """"only through delfin-tools" is enforced by no mechanism, so stating
    it is the whole of its enforcement."""
    sk = S.get_skill("pipeline-build")
    assert "only" in sk.body and "delfin-tools" in sk.body


def test_it_is_offered_where_its_server_is():
    """The MCP server ships built-in, so the coding session can follow it;
    an office session cannot, and is not told about it."""
    assert "pipeline-build" in {
        s.name for s in S.discover_skills(domain="code")}
    assert "pipeline-build" not in {
        s.name for s in S.discover_skills(domain="office")}


# ---------------------------------------------------------------------------
# The mode is gone, and the old name still works
# ---------------------------------------------------------------------------

def test_the_mode_is_retired(loader):
    assert "pipeline" not in loader.available_modes()


def test_the_old_name_still_opens_a_session():
    """A saved session, a slash command, a script -- none may break."""
    assert _migrate_mode("pipeline") == "solo"


def test_the_dashboard_no_longer_offers_it():
    source = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
              / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    start = source.index("options=[(\"Dashboard\", \"dashboard\")")
    assert "\"pipeline\"" not in source[start:start + 200]


def test_nothing_still_keys_behaviour_on_the_retired_mode():
    """A leftover ``mode == "pipeline"`` branch is now unreachable code
    that reads as if the mode were alive."""
    for rel in ("agent/prompt_loader.py", "agent/engine.py"):
        path = _PACK.parent / rel
        offenders = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            # Comments explain the retirement; the migration table is the
            # one place the name legitimately survives, as a redirect; and
            # "pipeline" is also an ordinary word in the task-complexity
            # patterns, which have nothing to do with modes.
            if '"pipeline"' in line
            and "mode" in line
            and not line.lstrip().startswith("#")
            and '"pipeline":' not in line
        ]
        assert not offenders, f"{rel} still branches on the mode: {offenders}"


# ---------------------------------------------------------------------------
# What the retirement must not have broken
# ---------------------------------------------------------------------------

def test_the_coding_prompt_did_not_grow(loader):
    """The move must cost the Code mode nothing -- the skill body is paid
    for when it is called, not before."""
    prompt = loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo",
        mode_description=loader.load_mode("solo")["description"],
        task_text="build a pipeline", session_key="probe-solo")
    assert "get_guide" not in prompt


def test_the_pipeline_coordinators_still_get_their_mode_text(loader):
    """The multi-role routes always received theirs and must keep it."""
    prompt = loader.build_system_prompt(
        role_id="session_manager", mode_id="quick",
        mode_description=loader.load_mode("quick")["description"],
        session_key="probe-quick")
    assert "Mode: quick" in prompt
