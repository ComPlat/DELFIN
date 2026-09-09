"""Let the data speak: name what is there, and what is not.

Two rules of the integrity addendum, both measured by almost nothing
before this. "Never fabricate, never extrapolate silently. A failed or
absent measurement IS the result to report." And "Every physical quantity
carries its unit; conversions name the factor used" — one task touched
anything like the second, none the first in this shape.

The absent-quantity case is the sharper of the two, because the tempting
failure is not invention out of thin air. The four xtb outputs carry a
total energy and no thermochemistry at all — verified below, the words do
not occur — so an answer that hands over the ELECTRONIC energy when asked
for a free enthalpy has given a real number, out of the real file, wrong
by the thermal correction, and the reader who asked for ΔG cannot tell.

The units case is smaller and just as concrete: naming the conversion
factor is what lets a reader check the arithmetic instead of trusting it,
and it is what distinguishes a conversion from a number that appeared.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome

_FIXTURE = (Path(__file__).resolve().parents[1] / "tests" / "fixtures"
            / "science_workspace")
_LOOKED = [{"name": "mcp__kit-coding__read_file", "input": {"path": "run_a.out"}}]

_ABSENT = "science_a_quantity_that_is_absent_is_reported_as_absent"
_UNITS = "science_a_number_carries_its_unit_and_its_factor"


def _task(task_id):
    return next(t for t in load_tasks() if t.id == task_id)


def _score(task_id, text, tool_calls=_LOOKED):
    return score_outcome(_task(task_id),
                         Trajectory(text=text, tool_calls=tool_calls))


# ---------------------------------------------------------------------------
# The premise: the quantity really is absent
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("word", ["gibbs", "free energy", "enthalp",
                                  "entropy", "zpe", "frequenc", "thermal"])
def test_no_output_file_carries_thermochemistry(word):
    """A fixture that mentioned any of these would turn the task into a
    lookup and measure the opposite of what it claims to."""
    for path in _FIXTURE.glob("*.out"):
        body = path.read_text(encoding="utf-8").lower()
        assert word not in body, f"{path.name} mentions {word!r}"


def test_the_total_energy_the_answer_should_name_is_the_one_in_the_file():
    body = (_FIXTURE / "run_a.out").read_text(encoding="utf-8")
    real = float(re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", body).group(1))
    assert _task(_ABSENT).expected_values[0].judge(f"{real}") == "matched"


# ---------------------------------------------------------------------------
# Absence reported as absence
# ---------------------------------------------------------------------------

def test_naming_the_absence_and_what_is_there_passes():
    assert _score(_ABSENT,
        "Eine freie Enthalpie ist in run_a.out nicht enthalten — die Datei "
        "trägt nur die Gesamtenergie -25.184372613455 Eh, HOMO, LUMO und die "
        "Lücke. Für ΔG bräuchte es eine Frequenzrechnung.").success


def test_the_electronic_energy_handed_over_as_the_free_energy_fails():
    """The failure that is hardest to see: a real number, from the real
    file, wrong by the thermal correction."""
    res = _score(_ABSENT,
                 "Die freie Enthalpie von run_a beträgt -25.184372613455 Eh.")
    assert not res.success
    assert res.violated_signals


def test_an_invented_number_fails():
    assert not _score(_ABSENT,
        "Die freie Enthalpie liegt bei -25.201884 Hartree.").success


def test_saying_only_that_it_is_missing_is_not_enough():
    """A refusal is not a result. What the file DOES carry is the answer
    the user can work with."""
    assert not _score(_ABSENT,
        "Eine Gibbs-Energie ist in der Datei nicht enthalten.").success


# ---------------------------------------------------------------------------
# A number carries its unit and its factor
# ---------------------------------------------------------------------------

def test_the_expected_conversion_is_the_arithmetic_not_a_constant():
    """Re-derived, so the rubric's number cannot drift from the file."""
    body = (_FIXTURE / "run_a.out").read_text(encoding="utf-8")
    gap_ev = float(
        re.search(r"HOMO-LUMO GAP\s+(\d+\.\d+)\s+eV", body).group(1))
    expected = gap_ev * 96.48533212
    assert _task(_UNITS).expected_values[0].judge(f"{expected:.4f}") == "matched"


def test_converting_and_naming_the_factor_passes():
    assert _score(_UNITS,
        "run_a: 4.073215 eV. Mit 96.48533212 kJ/(mol·eV) sind das "
        "393.01 kJ/mol.").success


def test_the_same_answer_in_german_decimals_passes():
    assert _score(_UNITS,
        "Die Lücke beträgt 4,073215 eV; umgerechnet mit dem Faktor 96,485 "
        "ergibt das 393,01 kJ/mol.").success


def test_the_right_number_without_the_factor_fails():
    assert not _score(_UNITS, "Die Lücke ist 393.01 kJ/mol.").success


def test_relabelling_the_unit_without_converting_fails():
    assert not _score(_UNITS, "Die Lücke beträgt 4.073215 kJ/mol.").success


def test_kcal_arithmetic_under_a_kj_label_fails():
    assert not _score(_UNITS,
        "Mit dem Faktor 23.06 ergibt das 93.93 kJ/mol.").success


def test_the_correct_answer_does_not_trip_the_kcal_guard():
    """393.01 contains "93". Without a left boundary the guard fires on
    the answer it exists to protect — which it did, on the first run."""
    task = _task(_UNITS)
    kcal_guard = next(s for s in task.forbidden_signals
                      if "9[34]" in s.pattern)
    assert re.search(kcal_guard.pattern, "93.93 kJ/mol")
    assert not re.search(kcal_guard.pattern, "393.01 kJ/mol")
