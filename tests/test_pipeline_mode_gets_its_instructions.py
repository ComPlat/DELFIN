"""Pipeline mode's instructions were never sent to the model.

modes/pipeline.md says: call ``get_guide`` first, ``validate_spec``
before every run, work only through the ``delfin-tools`` server, never
hand-roll pipeline logic. None of it reached the model.

The gate that injects a mode description checks the ROLE, and it sits
after the ``solo_agent`` branch returns. Pipeline routes to solo_agent,
so the branch returned first and the gate was unreachable. Pipeline was
Code mode under a different label.

It was also a more expensive one. Lazy-module stripping and the context
distiller are both keyed on mode strings that did not include it, so a
pipeline turn shipped every optional prompt module and skipped
distillation: about 11,000 characters more prompt than solo, for the same
agent, with none of its own instructions.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent.prompt_loader import PromptLoader

_PACK = pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"


@pytest.fixture
def loader():
    return PromptLoader(_PACK)


def _prompt(loader, mode, role="solo_agent", task_text="", session_key=""):
    data = loader.load_mode(mode)
    return loader.build_system_prompt(
        role_id=role, mode_id=mode,
        mode_description=data["description"],
        task_text=task_text, session_key=session_key or f"probe-{mode}")


# ---------------------------------------------------------------------------
# The instructions arrive
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("marker", [
    "get_guide", "validate_spec", "delfin-tools", "Pipeline Builder",
])
def test_the_pipeline_instructions_reach_the_model(loader, marker):
    assert marker in _prompt(loader, "pipeline"), marker


def test_the_rule_the_prompt_exists_to_state_is_present(loader):
    """"only delfin-tools" is unenforced by any mechanism, so delivering it
    is the whole of its enforcement. Undelivered, it was not even a weak
    rule -- it was no rule."""
    text = _prompt(loader, "pipeline")
    assert "only" in text and "delfin-tools" in text


# ---------------------------------------------------------------------------
# ...without giving every other mode a second description of itself
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode,role", [
    ("solo", "solo_agent"),
    ("office", "office_agent"),
])
def test_a_mode_whose_file_only_describes_it_is_not_injected(loader, mode, role):
    """solo.md and dashboard.md tell a human what the dropdown entry means.
    Injecting those spends tokens restating what the role prompt says."""
    data = loader.load_mode(mode)
    first_line = next(
        (l.strip() for l in data["description"].splitlines()
         if len(l.strip()) > 30 and not l.strip().startswith("#")), "")
    assert first_line
    assert first_line[:40] not in _prompt(loader, mode, role)


def test_the_pipeline_coordinators_still_get_theirs(loader):
    """The two roles that always received it must not have lost it."""
    assert "Mode: quick" in _prompt(loader, "quick", "session_manager")


# ---------------------------------------------------------------------------
# And it costs less than before, not more
# ---------------------------------------------------------------------------

def test_pipeline_no_longer_ships_more_prompt_than_solo(loader):
    """It ran the same agent and paid ~11k characters more, because the
    stripping and the distiller were keyed on names that omitted it."""
    pipeline = len(_prompt(loader, "pipeline", task_text="build a pipeline",
                           session_key="s-pipeline"))
    solo = len(_prompt(loader, "solo", task_text="build a pipeline",
                       session_key="s-solo"))
    assert pipeline - solo < 4000, (
        f"pipeline is {pipeline - solo} chars above solo; the lazy-module "
        "stripping is off for it again")


def test_lazy_module_stripping_applies_to_pipeline():
    source = pathlib.Path(
        __import__("delfin.agent.prompt_loader", fromlist=["x"]).__file__
    ).read_text(encoding="utf-8")
    assert 'mode_id not in ("solo", "plan", "pipeline")' in source


def test_the_distiller_is_enabled_for_pipeline():
    source = pathlib.Path(
        __import__("delfin.agent.engine", fromlist=["x"]).__file__
    ).read_text(encoding="utf-8")
    # Anchor on the enable list itself: "ContextDistiller" appears
    # earlier as a type annotation, and indexing the first hit tested a
    # docstring rather than the decision.
    start = source.index("_enable = self.mode in (")
    assert '"pipeline"' in source[start:start + 200]
