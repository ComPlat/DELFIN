"""A navigation request was carrying the whole ORCA manual grounding.

``dashboard_agent.md`` had no module markers at all, so the entire role
prompt shipped on every turn — the manual-grounding rules, the Builder
field tables, the CONTROL.txt key reference and the calculation-failure
playbook went out with "wechsel zu Submit". Measured over the 28
dashboard utterances in the benchmark corpus, half trigger the chemistry
module and half trigger nothing.

Two things had to be true for the markers to do anything: the dashboard
mode had to be inside the trigger heuristic, and stripping had to run for
roles other than solo. Neither was, and a file with markers and no
stripping reads exactly like a working gate — the markers vanish either
way.
"""

from pathlib import Path

import pytest

from delfin.agent.prompt_loader import PromptLoader

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack"

_CHEMISTRY_ONLY = (
    "Ground every ORCA",
    "ORCA Builder capabilities",
    "CONTROL.txt — quick reference",
    "Calculation data search",
)
# What every dashboard turn needs, whatever the user asked for.
_ALWAYS = (
    "How `ACTION:` works",
    "Hard scope limits",
    "Safety rules",
    "Tools you may NOT use in dashboard mode",
    "Opening / reading files in calc folders",
)


def _built(task: str, key: str) -> str:
    return PromptLoader().build_system_prompt(
        role_id="dashboard_agent", mode_id="dashboard", task_text=task,
        session_key=key, model="kit.glm-5.3")


@pytest.mark.parametrize("section", _ALWAYS)
def test_the_dashboard_contract_ships_on_every_turn(section):
    assert section in _built("wechsel zu Submit", f"a-{section[:6]}")


@pytest.mark.parametrize("section", _CHEMISTRY_ONLY)
def test_a_navigation_request_carries_no_chemistry_reference(section):
    assert section not in _built("wechsel zu Submit", f"n-{section[:6]}")


@pytest.mark.parametrize("section", _CHEMISTRY_ONLY)
def test_a_chemistry_request_gets_all_of_it(section):
    assert section in _built("stell den orca functional auf b3lyp",
                             f"c-{section[:6]}")


def test_the_saving_is_real_and_measured():
    nav = PromptLoader().prompt_size_report(
        role_id="dashboard_agent", mode_id="dashboard",
        task_text="wechsel zu Submit", session_key="m1", model="kit.glm-5.3")
    chem = PromptLoader().prompt_size_report(
        role_id="dashboard_agent", mode_id="dashboard",
        task_text="stell den orca functional auf b3lyp", session_key="m2",
        model="kit.glm-5.3")
    saved = chem["stable_tokens"] - nav["stable_tokens"]
    assert saved > 1_200, saved


def test_the_dashboard_mode_is_inside_the_heuristic():
    """Without this the markers are stripped as comments and every block
    ships anyway — a gate that looks like a gate and gates nothing."""
    loader = PromptLoader()
    every = set(loader._MODULE_TRIGGERS)
    assert loader._detect_active_modules("wechsel zu Submit", "dashboard") != every
    # A mode that is NOT in the heuristic still gets everything.
    assert loader._detect_active_modules("wechsel zu Submit", "builder") == every


def test_stripping_runs_for_a_role_that_is_not_solo():
    """It used to run only inside the solo branch of compose_sections."""
    raw = (_PACK / "agents" / "dashboard_agent.md").read_text()
    assert "<!-- module:" in raw
    built = _built("wechsel zu Submit", "s1")
    assert len(built) < len(raw) + 20_000     # the addenda are the rest
    assert "Ground every ORCA" not in built


def test_a_role_file_without_markers_is_untouched():
    """The stripping call is new on this path; every other role must read
    exactly what it read before."""
    loader = PromptLoader()
    for role in ("builder_agent", "critic_agent", "office_agent"):
        raw = loader.load_role_prompt(role)
        if not raw or "<!-- module:" in raw:
            continue
        out = loader._strip_lazy_modules(
            raw, task_text="Hallo", mode_id="dashboard", model="kit.glm-5.3",
            session_key=f"u-{role}", role_id=role)
        assert out == raw, role
