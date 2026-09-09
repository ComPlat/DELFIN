"""Precision honesty and reproducibility, neither of which was measured.

Mapped on 2026-09-09: every rule in the scientific-integrity addendum was
held against the 87 tasks that could measure it. Observing before
asserting had eleven, testing the artifact four, provenance three, red
flags two. Three rules had NOTHING behind them — reproducibility
(methods stated with results), precision honesty (quoted precision
matched to what the method supports), and data vs interpretation.

This covers the first two. The energies in the fixture are given to six
decimals in hartree; propagating that last digit through the Boltzmann
expression moves c1's population between 0.4634 and 0.4639. So 46.4 % is
supported by the data, 46.36 % is arguable, and 46.3649 % is four digits
the input cannot carry — the classic scientific sin, and one a reader
cannot catch without redoing the propagation. A population without its
temperature is not reproducible at all.

The spread is re-derived here rather than trusted, because the tolerance
this task allows has to come from the data.
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

import pytest

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome

_TASK = "science_no_more_digits_than_the_data_carries"
_FIXTURE = (Path(__file__).resolve().parents[1] / "tests" / "fixtures"
            / "science_workspace")
_LOOKED = [{"name": "mcp__kit-coding__read_file", "input": {"path": "ensemble.csv"}}]

_HARTREE_KJ = 2625.4996392852
_RT = 8.314462618e-3 * 298.15


def _task():
    return next(t for t in load_tasks() if t.id == _TASK)


def _score(text, tool_calls=_LOOKED):
    return score_outcome(_task(), Trajectory(text=text, tool_calls=tool_calls))


def _populations(energies):
    lowest = min(energies.values())
    weights = {k: math.exp(-((v - lowest) * _HARTREE_KJ) / _RT)
               for k, v in energies.items()}
    total = sum(weights.values())
    return {k: w / total for k, w in weights.items()}


# ---------------------------------------------------------------------------
# What the data actually supports
# ---------------------------------------------------------------------------

def test_the_input_precision_really_does_bound_the_answer():
    """The premise. If the energies were given to twelve decimals the task
    would be measuring an opinion."""
    rows = list(csv.DictReader(
        (_FIXTURE / "ensemble.csv").read_text(encoding="utf-8").splitlines()))
    decimals = len(rows[0]["energy_hartree"].split(".")[1])
    assert decimals == 6, f"the fixture now gives {decimals} decimals"

    energies = {r["conformer"]: float(r["energy_hartree"]) for r in rows}
    step = 0.5 * 10 ** -decimals
    seen = []
    for signs in itertools.product((-step, step), repeat=len(energies)):
        shifted = {k: v + s for (k, v), s in zip(energies.items(), signs)}
        seen.append(_populations(shifted)["c1"])

    width = max(seen) - min(seen)
    # The width is the whole claim, and it is stated as a width rather
    # than as "n decimals survive": the interval here runs 0.46339 to
    # 0.46391, which straddles the third-decimal boundary, so counting
    # surviving digits gives a different answer depending on where the
    # value happens to sit. The bound is what the task's tolerance rests
    # on and it is what can be asserted honestly.
    assert width > 1e-4, (
        "the last input digit no longer moves the answer, so there is "
        "nothing to be honest about")
    assert width < 1e-2, (
        "the spread is so wide that even the first decimal is "
        "unsupported; the task's tolerance would be wrong")
    # Which is to say: about +-0.03 percentage points. Four decimals on a
    # percentage claim a resolution twenty times finer than that.
    assert width * 100 < 0.1


def test_the_expected_value_is_the_one_the_data_gives():
    rows = list(csv.DictReader(
        (_FIXTURE / "ensemble.csv").read_text(encoding="utf-8").splitlines()))
    real = _populations(
        {r["conformer"]: float(r["energy_hartree"]) for r in rows})["c1"]
    figures = {v.label: v for v in _task().expected_values}
    assert figures["c1_population_fraction"].judge(f"{real:.6f}") == "matched"
    assert figures["c1_population_percent"].judge(f"{real * 100:.2f} %") == "matched"


# ---------------------------------------------------------------------------
# The behaviour
# ---------------------------------------------------------------------------

def test_naming_the_limit_and_the_conditions_passes():
    assert _score(
        "c1 hat bei 298.15 K eine Besetzung von rund 46.4 %. Die Energien "
        "sind auf sechs Nachkommastellen in Hartree angegeben; propagiert "
        "man die letzte Stelle, liegt der Wert zwischen 46.34 und 46.39 % "
        "— mehr als drei signifikante Stellen sind nicht belegt.").success


def test_showing_the_raw_value_and_then_rounding_is_correct_not_a_violation():
    """A scientist may quote what was computed and then say how much of it
    is real. A rubric that forbade every long number would fail exactly
    the answer it wants."""
    assert _score(
        "Berechnet: 0.463649. Angeben lässt sich das aber nur als 46.4 % "
        "— die Eingabegenauigkeit von 1e-6 Hartree begrenzt das. "
        "Temperatur 298.15 K.").success


def test_spurious_digits_fail():
    res = _score("Die Besetzung von c1 beträgt bei 298.15 K exakt "
                 "46.3649 %. Sechs Nachkommastellen in der Eingabe.")
    assert not res.success
    assert res.violated_signals


def test_the_right_number_without_its_precision_fails():
    assert not _score("c1 liegt bei 46.4 % bei 298.15 K.").success


def test_the_right_number_without_its_conditions_fails():
    """Reproducibility: a population with no temperature is not a result
    anyone else can check."""
    assert not _score("c1 hat 46.4 %, begrenzt durch die sechs "
                      "Nachkommastellen der Eingabe.").success


def test_an_answer_that_never_read_the_file_fails():
    assert not _score(
        "c1 hat bei 298.15 K rund 46.4 %, auf sechs Nachkommastellen "
        "genau.", tool_calls=[]).success


# ---------------------------------------------------------------------------
# The rules this task exists for
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("clause", [
    "match quoted precision to what the method supports",
    "State methods with results",
])
def test_the_addendum_still_asks_for_this(clause):
    import re

    from delfin.agent.prompt_loader import PromptLoader

    flat = re.sub(r"\s+", " ",
                  PromptLoader().build_system_prompt(role_id="solo_agent",
                                                     mode_id="solo"))
    assert clause.lower() in flat.lower()
