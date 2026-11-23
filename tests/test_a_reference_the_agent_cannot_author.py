"""An agent can fabricate a report; it cannot fabricate agreement.

The execution contract answers *did this run compute*. It cannot answer
*is the number right*, and a confident report over a wrong number costs
exactly as little to write as one over a right number. A value published
by someone else closes that gap, and a survey of 26 autonomous-research
systems found that under its own coding rule no system of the LLM era
has such an oracle inside the loop (arXiv:2608.05179).

What is asserted here, in order of what it protects:

  * **not_measured is never agreement**, and it does not raise coverage.
    A reference the run did not answer stays in the denominator. The
    opposite -- dropping the case it failed -- is how a coverage figure
    becomes a different measurement with the same name.
  * **an average over nothing is not zero.** mean_absolute_deviation is
    None with no comparable value, because 0.0 reads as perfect.
  * **the matching rule is named in every verdict.** Assigning a
    computed transition to a published one is the hard half: state
    ordering differs between methods and "closest in energy" rewards a
    wrong assignment. The rule is declared per reference and carried.
  * **a tolerance belongs to its reference**, not to this module. A
    coupled-cluster benchmark value and an experimental band maximum do
    not deserve the same window.
  * **the reference set is input.** Nothing is embedded, so nothing here
    inherits a data licence, and the comparator works with any set.

The fixtures below are ORCA output fragments written for this test, not
published data, and the reference values are invented numbers chosen to
exercise agreement and disagreement. Nothing in this file cites a real
measurement as though it were one.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import reference_oracle as ro

#: The block delfin.api.extract_excited_states reads. Two transitions:
#: a dark one at 3.00 eV and a bright one at 4.00 eV, so the matching
#: rules can be told apart -- "lowest" and "first bright" give different
#: answers, which is the point of declaring which one a reference uses.
_OUT = """\
-----------------------------------------------------------------------
         ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-----------------------------------------------------------------------
      Transition      Energy     Energy  Wavelength     fosc        D2
                       (eV)      (cm-1)      (nm)                (au**2)
-----------------------------------------------------------------------
  0-1A  ->  1-1A     3.000000   24196.0     413.3   0.000100000   0.00100
  0-1A  ->  2-1A     4.000000   32261.0     309.9   0.350000000   3.50000
"""


@pytest.fixture()
def workspace(tmp_path):
    """Two folders, as a run comparing two functionals leaves behind."""
    for name in ("pna_pbe0", "pna_wb97xd"):
        d = tmp_path / name
        d.mkdir()
        (d / f"{name}.out").write_text(_OUT, encoding="utf-8")
    return tmp_path


def _ref(**kw):
    base = dict(key="pna_pbe0", quantity="vertical_excitation_ev",
                value=4.00, tolerance=0.10,
                source="fixture value written for this test, not published",
                match=ro.MATCH_FIRST_BRIGHT)
    base.update(kw)
    return ro.Reference(**base)


def _by_key(report, key):
    return next(v for v in report.verdicts if v.key == key)


# --- the rule the oracle exists for -----------------------------------

def test_a_value_within_the_tolerance_agrees(workspace):
    rep = ro.compare(workspace, [_ref(value=4.00, tolerance=0.10)])
    v = _by_key(rep, "pna_pbe0")
    assert v.verdict == ro.AGREES
    assert v.computed == pytest.approx(4.0, abs=0.01)
    assert v.deviation == pytest.approx(0.0, abs=0.01)


def test_a_value_outside_the_tolerance_differs(workspace):
    rep = ro.compare(workspace, [_ref(value=3.20, tolerance=0.10)])
    v = _by_key(rep, "pna_pbe0")
    assert v.verdict == ro.DIFFERS
    assert v.deviation == pytest.approx(0.8, abs=0.02)


def test_a_reference_the_run_did_not_answer_is_not_agreement(workspace):
    """The whole point. A folder that is not there is not a pass."""
    rep = ro.compare(workspace, [_ref(key="never_computed")])
    v = _by_key(rep, "never_computed")
    assert v.verdict == ro.NOT_MEASURED
    assert v.computed is None
    assert "never_computed" in v.why


def test_a_missing_case_stays_in_the_denominator(workspace):
    """Coverage must not rise by dropping what failed. Two references,
    one answerable: coverage is a half, not a whole."""
    rep = ro.compare(workspace, [_ref(), _ref(key="never_computed")])
    assert rep.coverage() == pytest.approx(0.5)
    assert len(rep.verdicts) == 2


def test_an_average_over_nothing_is_not_zero(workspace):
    """0.0 eV mean deviation reads as perfect agreement. None says what
    is true: nothing was comparable."""
    rep = ro.compare(workspace, [_ref(key="never_computed")])
    assert rep.mean_absolute_deviation() is None
    assert rep.coverage() == 0.0


def test_the_average_ignores_the_cases_it_could_not_measure(workspace):
    """Averaging a deviation over references that produced none would
    mix two different quantities under one name."""
    rep = ro.compare(workspace, [
        _ref(value=4.00, tolerance=0.10),
        _ref(key="pna_wb97xd", value=3.60, tolerance=0.10),
        _ref(key="never_computed")])
    assert rep.mean_absolute_deviation() == pytest.approx(0.2, abs=0.02)


# --- the matching rule is the hard half, so it is declared ------------

def test_the_rules_give_different_answers_and_say_which_was_used(workspace):
    """A dark state sits below the bright one on purpose: "lowest" and
    "first bright" must not be interchangeable, and a verdict that did
    not say which it used would be unreadable."""
    bright = _by_key(ro.compare(workspace, [_ref(match=ro.MATCH_FIRST_BRIGHT)]),
                     "pna_pbe0")
    lowest = _by_key(ro.compare(workspace, [_ref(match=ro.MATCH_LOWEST)]),
                     "pna_pbe0")
    assert bright.computed == pytest.approx(4.0, abs=0.01)
    assert lowest.computed == pytest.approx(3.0, abs=0.01)
    assert bright.match == ro.MATCH_FIRST_BRIGHT
    assert lowest.match == ro.MATCH_LOWEST
    assert bright.why and lowest.why and bright.why != lowest.why


def test_an_index_past_the_table_is_not_measured(workspace):
    rep = ro.compare(workspace, [_ref(match=ro.MATCH_INDEX, index=9)])
    v = _by_key(rep, "pna_pbe0")
    assert v.verdict == ro.NOT_MEASURED
    assert "9" in v.why


def test_no_bright_transition_is_reported_not_substituted(tmp_path):
    """When nothing reaches the extractor's own bright threshold, the
    oracle says so instead of falling back to the lowest state -- a
    silent fallback would compare two different quantities."""
    d = tmp_path / "dark_only"
    d.mkdir()
    (d / "dark_only.out").write_text(
        _OUT.replace("0.350000000", "0.000000100"), encoding="utf-8")
    rep = ro.compare(tmp_path, [_ref(key="dark_only")])
    v = _by_key(rep, "dark_only")
    assert v.verdict == ro.NOT_MEASURED
    assert "bright" in v.why


# --- the reference set is input, and a bad row is reported ------------

def test_a_reference_file_is_read_and_its_bad_rows_are_named(tmp_path,
                                                             workspace):
    path = tmp_path / "refs.json"
    path.write_text(json.dumps({"references": [
        {"key": "pna_pbe0", "value": 4.0, "tolerance": 0.1,
         "source": "a fixture, not published data"},
        {"key": "broken", "value": "not a number", "tolerance": 0.1,
         "source": "x"},
        {"key": "unknown_rule", "value": 1.0, "tolerance": 0.1,
         "source": "x", "match": "closest"},
    ]}), encoding="utf-8")
    rep = ro.compare_file(workspace, path)
    assert len(rep.verdicts) == 1
    assert _by_key(rep, "pna_pbe0").verdict == ro.AGREES
    joined = " ".join(rep.notes)
    assert "row 1" in joined and "row 2" in joined
    assert "closest" in joined


def test_an_unreadable_reference_file_compares_nothing_and_says_so(tmp_path,
                                                                   workspace):
    bad = tmp_path / "nope.json"
    bad.write_text("{not json", encoding="utf-8")
    rep = ro.compare_file(workspace, bad)
    assert rep.verdicts == []
    assert rep.coverage() is None
    assert any("not the same as agreement" in n for n in rep.notes)


def test_every_verdict_carries_the_citation_it_was_judged_against():
    """A number with no traceable source is not a reference, and a
    reader of the verdict must be able to follow it."""
    rep = ro.compare(".", [_ref(source="J. Example Chem. 2026, 1, 1")])
    assert all(v.source for v in rep.verdicts)


def test_a_missing_workspace_is_not_measured_rather_than_a_difference(tmp_path):
    rep = ro.compare(tmp_path / "nowhere", [_ref()])
    assert _by_key(rep, "pna_pbe0").verdict == ro.NOT_MEASURED
    assert rep.coverage() == 0.0


def test_the_module_embeds_no_published_data():
    """The licence argument, as a test rather than a promise: a future
    edit that pastes a reference table in here would make this repository
    inherit that table's licence."""
    import pathlib
    src = pathlib.Path(ro.__file__).read_text(encoding="utf-8")
    assert "references" not in json.dumps(
        [ln for ln in src.splitlines() if ln.strip().startswith("{")])
    assert not [ln for ln in src.splitlines()
                if ln.strip().startswith('"value"')]


# --- what a person reads -----------------------------------------------

def test_the_rendering_puts_disagreement_first(workspace):
    rep = ro.compare(workspace, [
        _ref(value=4.00, tolerance=0.10),
        _ref(key="pna_wb97xd", value=1.00, tolerance=0.10),
        _ref(key="never_computed")])
    text = ro.render(rep)
    assert text.index(ro.DIFFERS) < text.index(ro.NOT_MEASURED)
    assert text.index(ro.NOT_MEASURED) < text.index("  " + ro.AGREES)
    assert "coverage 67%" in text


def test_the_report_is_json_serialisable(workspace):
    data = ro.compare(workspace, [_ref()]).as_dict()
    json.dumps(data)
    assert data["counts"][ro.AGREES] == 1
    assert data["coverage"] == 1.0
