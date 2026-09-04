"""The notes the user reads follow the session, and say each thing once.

From the recorded run, report 20260904-104748. The answer was English
throughout — the session-language pin held — and underneath it the
harness stapled this:

    [verify] Caveat: für 'elapsed_s' (nEOF) liegen zwei Werte vor,
    0.055 und 0.943 — 94.2% auseinander …

    > ⚠️ Diese Antwort nennt 27 components, 27 components, 27 components,
    aber die Ausgabe von … wurde in diesem Zug abgeschnitten …

Three defects in two lines, all in the guard rather than the model:
the caveats were hardcoded German under an English answer; `elapsed_s`
is the harness's own stopwatch and was compared as if it were a figure
about a record; and one claim was printed once per place it appeared.
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as vg


@pytest.fixture(autouse=True)
def _fresh():
    vg.reset_keyed_values()
    vg.set_caveat_language("de")
    yield
    vg.reset_keyed_values()
    vg.set_caveat_language("de")


# ------------------------------------------------- the stopwatch is not a figure


def test_a_tool_runtime_is_never_compared_as_a_figure():
    """0.055 s and 0.943 s are two shell calls, not a contradiction."""
    vg.record_keyed_values(
        '{"exit_code": 0, "elapsed_s": 0.055, "stdout": "x"}',
        source='{"path": "/w/archive/MOL_1/DELFIN_Data.json"}')
    vg.record_keyed_values(
        '{"exit_code": 0, "elapsed_s": 0.943, "stdout": "y"}',
        source='{"path": "/w/archive/MOL_1/summary.csv"}')
    assert vg.scan_for_conflicting_figures() == []


def test_a_real_figure_in_the_same_envelope_is_still_compared():
    """The exclusion is the envelope, not the whole result."""
    vg.record_keyed_values(
        '{"exit_code": 0, "elapsed_s": 0.055, "beta_HRS_au": 171232.01}',
        source='{"path": "/w/archive/MOL_1/DELFIN_Data.json"}')
    vg.record_keyed_values(
        '{"exit_code": 0, "elapsed_s": 0.943, "beta_HRS_au": 180721.43}',
        source='{"path": "/w/archive/MOL_1/summary.csv"}')
    flags = vg.scan_for_conflicting_figures()
    assert [f.field for f in flags] == ["beta_HRS_au"]


# ------------------------------------------------------- one claim, once


def test_a_claim_is_named_once_however_often_it_was_counted():
    caveat = vg.truncated_output_caveat(
        ["27 components", "27 components", "27 components"], ["bash"])
    assert caveat.count("27 components") == 1


def test_two_different_claims_are_both_named():
    caveat = vg.truncated_output_caveat(
        ["27 components", "404 folders", "27 components"], ["bash"])
    assert "27 components" in caveat and "404 folders" in caveat
    assert caveat.count("27 components") == 1


# --------------------------------------------------------- the language


def test_the_warning_is_german_in_a_german_session():
    vg.set_caveat_language("de")
    assert "abgeschnitten" in vg.truncated_output_caveat(["5 Zeilen"], ["bash"])


def test_the_warning_is_english_in_an_english_session():
    """The defect as filmed: an English answer with German warnings."""
    vg.set_caveat_language("en")
    caveat = vg.truncated_output_caveat(["27 components"], ["bash"])
    assert "was truncated in this turn" in caveat
    assert "abgeschnitten" not in caveat


def test_an_unknown_language_leaves_it_alone():
    """A half-translated caveat is worse than one consistent language."""
    vg.set_caveat_language("en")
    vg.set_caveat_language("fr")
    assert "was truncated" in vg.truncated_output_caveat(["3 rows"], ["bash"])
