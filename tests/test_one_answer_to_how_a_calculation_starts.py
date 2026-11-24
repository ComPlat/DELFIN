"""The prompt forbade what the gate admitted.

MEASURED, three mechanisms with three answers to "how do I start a
calculation":

  solo_agent.md      "Real calculations start only through the
                     sanctioned tools (/orca submit,
                     submit_calculation, pipeline_run) -- never by
                     launching ORCA, xTB or SLURM from bash."
  the shell gate     a checked `bash orca job.inp` is ADMITTED with no
                     dialog (calc_policy, `calculation_admitted`).
  the adapter route  was refused unattended for want of a dialog, until
                     judge_application_run.

So the prompt banned the route the gate rewarded, and named neither the
building blocks nor the specialist. The model took what worked: 20
hand-written inputs on one measured run.

ONE ANSWER NOW, and it matches the gate: the blocks or the CONTROL.txt
route, bash allowed and gate-checked, and what a by-hand run does not
carry. A prompt that forbids what the gate permits teaches that its
rules are decoration.

AND IT IS PAID FOR. The always-on line is 160 characters against the
165 of the prohibition it replaced, so the head's budget is unchanged;
the tool names and the specialist live in the lazy `chemistry` module,
which fires on "orca", "dft", "single point", "anregung" and the rest --
0 characters on a coding turn.
"""
from __future__ import annotations

import pytest

from delfin.agent.prompt_loader import PromptLoader

CODING = "please rename this variable and add a unit test"
CHEMISTRY = "compute a single point with ORCA for this complex"


@pytest.fixture(scope="module")
def loader():
    return PromptLoader()


def test_the_prohibition_the_gate_contradicts_is_gone(loader):
    for task in (CODING, CHEMISTRY):
        text = loader.build_system_prompt("solo_agent", "solo",
                                          task_text=task)
        assert "never by launching ORCA" not in text


def test_the_always_on_line_states_both_routes_and_the_cost(loader):
    text = loader.build_system_prompt("solo_agent", "solo", task_text=CODING)
    assert "building blocks or the CONTROL.txt route" in text
    assert "gate-checked" in text
    assert "declares no parameters" in text


@pytest.mark.parametrize("needle", ["resolve_spec", "run_application",
                                    "get_guide"])
def test_the_route_arrives_on_a_chemistry_turn(loader, needle):
    assert needle in loader.build_system_prompt(
        "solo_agent", "solo", task_text=CHEMISTRY)


@pytest.mark.parametrize("needle", ["resolve_spec", "run_application",
                                    "get_guide"])
def test_and_costs_a_coding_turn_nothing(loader, needle):
    assert needle not in loader.build_system_prompt(
        "solo_agent", "solo", task_text=CODING)


def test_the_decision_tree_has_a_row_for_starting_one(loader):
    """It had a row for reading every kind of result and none for
    running a calculation, which is the question the measured run got
    wrong."""
    text = loader.build_system_prompt("solo_agent", "solo",
                                      task_text=CHEMISTRY)
    assert "RUN / submit a calculation" in text
    assert "extract_imaginary_frequencies" in text, "the tree is still there"


def test_the_specialist_is_named_on_both(loader):
    """The subagent table is always on, so the specialist has a cue
    whatever the turn is about."""
    for task in (CODING, CHEMISTRY):
        assert "workflow-operator" in loader.build_system_prompt(
            "solo_agent", "solo", task_text=task)
