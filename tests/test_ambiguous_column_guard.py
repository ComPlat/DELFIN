"""A figure over a column the reader could not read gets marked as such.

The reader already says a column is undecidable, and the office role
prompt already forbids computing on it. On the benchmark the model was
told exactly that, chose a reading anyway, and answered with a single
total plus its own note admitting the assumption. In a report handed to
someone else that figure is indistinguishable from a measured one — the
note travels with the answer, not with the number.

Neither the tool result nor the prompt bound the model, which is why
this is a scanner and not a sentence.
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as vg

_NOTE = (
    "column 'Anschaffungswert': values like '8.986' read as two different "
    "numbers depending on the convention, and nothing in the column decides "
    "it. Ask which reading is meant — do not compute on it as it is."
)


# ---------------------------------------------------------------------------
# Reading the reader's own note
# ---------------------------------------------------------------------------

def test_the_flagged_column_is_taken_from_the_note():
    """Read rather than re-derived, so the two cannot disagree about which
    column is in question."""
    assert vg.extract_ambiguous_columns(_NOTE) == ["Anschaffungswert"]


def test_the_example_value_does_not_break_the_match():
    """The gap carries "values like '8.986'" — the one character the note
    is about is a dot, and an earlier pattern excluded dots there."""
    assert "8.986" in _NOTE
    assert vg.extract_ambiguous_columns(_NOTE)


def test_several_columns_are_all_collected():
    both = _NOTE + "\nNOTE: column 'Restwert': values like '1.234' read as " \
                   "two different numbers depending on the convention"
    assert vg.extract_ambiguous_columns(both) == ["Anschaffungswert", "Restwert"]


def test_other_column_notes_are_not_mistaken_for_it():
    for other in ("column 'Betrag': numbers use a decimal comma",
                  "column 'Datum': dates are day-first",
                  "column 'Betrag': 1 of 4 value(s) are not numbers"):
        assert vg.extract_ambiguous_columns(other) == []


# ---------------------------------------------------------------------------
# What the scanner fires on, and what it must leave alone
# ---------------------------------------------------------------------------

def test_a_single_total_is_flagged():
    text = ("Gesamtwert: 25.136 EUR. Hinweis: Anschaffungswert im deutschen "
            "Zahlenformat interpretiert.")
    assert vg.scan_for_totals_over_ambiguous_columns(
        text, ["Anschaffungswert"]) == ["Anschaffungswert"]


def test_the_models_own_note_does_not_excuse_the_figure():
    """It travels with the answer, not with the number: whoever copies the
    figure into a report leaves the note behind."""
    text = "Gesamtwert Anschaffungswert: 25.136 EUR (Annahme: deutsches Format)"
    assert vg.scan_for_totals_over_ambiguous_columns(text, ["Anschaffungswert"])


def test_offering_both_readings_and_asking_passes():
    """This is the correct answer and must not be caveated."""
    text = ("Die Werte im Anschaffungswert sind mehrdeutig: 8.986 = 8986 oder "
            "8,986. Deutsche Lesart ergäbe Summe 25.136 EUR, englische "
            "25,136 EUR. Welche Lesart gilt für euer Inventar?")
    assert vg.scan_for_totals_over_ambiguous_columns(
        text, ["Anschaffungswert"]) == []


def test_a_question_alone_is_enough_to_pass():
    text = "Anschaffungswert: Summe wäre 25.136 EUR — welche Konvention gilt?"
    assert vg.scan_for_totals_over_ambiguous_columns(
        text, ["Anschaffungswert"]) == []


def test_an_answer_without_a_total_is_left_alone():
    text = "Die Datei enthält 3 Positionen in der Spalte Anschaffungswert."
    assert vg.scan_for_totals_over_ambiguous_columns(
        text, ["Anschaffungswert"]) == []


def test_a_total_over_a_different_column_is_left_alone():
    """A figure in a reply that never mentions the flagged column is some
    other figure."""
    text = "Gesamtbetrag der Buchungen: 1.986,40 EUR."
    assert vg.scan_for_totals_over_ambiguous_columns(
        text, ["Anschaffungswert"]) == []


def test_nothing_flagged_means_nothing_scanned():
    text = "Gesamtwert: 25.136 EUR."
    assert vg.scan_for_totals_over_ambiguous_columns(text, []) == []


def test_the_caveat_names_the_column_and_the_reason():
    caveat = vg.ambiguous_column_caveat(["Anschaffungswert"])
    assert "Anschaffungswert" in caveat
    assert "8.986" in caveat
    assert "nicht gemessen" in caveat


def test_no_columns_no_caveat():
    assert vg.ambiguous_column_caveat([]) == ""


# ---------------------------------------------------------------------------
# The engine collects the evidence and appends the caveat
# ---------------------------------------------------------------------------

class _Engine:
    """The two engine methods under test, without building a real engine."""

    def __init__(self):
        from delfin.agent.engine import AgentEngine

        self._ambiguous_columns_turn = []
        self.messages = [{"role": "assistant", "content": ""}]
        self._note_ambiguous_columns = AgentEngine._note_ambiguous_columns.__get__(self)
        self._scan_ambiguous_column_totals = (
            AgentEngine._scan_ambiguous_column_totals.__get__(self))
        self._append_ambiguous_column_caveat = (
            AgentEngine._append_ambiguous_column_caveat.__get__(self))


def test_the_engine_records_what_the_reader_flagged():
    eng = _Engine()
    eng._note_ambiguous_columns(_NOTE)
    assert eng._ambiguous_columns_turn == ["Anschaffungswert"]


def test_the_same_column_is_recorded_once():
    eng = _Engine()
    eng._note_ambiguous_columns(_NOTE)
    eng._note_ambiguous_columns(_NOTE)
    assert eng._ambiguous_columns_turn == ["Anschaffungswert"]


def test_a_tool_result_without_the_note_records_nothing():
    eng = _Engine()
    eng._note_ambiguous_columns("Betrag: number (decimal_comma)")
    assert eng._ambiguous_columns_turn == []


def test_the_caveat_lands_in_the_answer_and_the_transcript():
    eng = _Engine()
    eng._note_ambiguous_columns(_NOTE)
    # The answer as the model actually produced it on the benchmark.
    answer = ("Gesamtwert: 25.136 EUR. Hinweis: die Werte im "
              "Anschaffungswert sind im deutschen Zahlenformat.")
    eng.messages[-1]["content"] = answer
    flagged = eng._scan_ambiguous_column_totals(answer)
    out = eng._append_ambiguous_column_caveat(answer, flagged)
    assert "nicht gemessen" in out
    # Also in the message the next turn will read, so the claim does not
    # stand bare in the history.
    assert "nicht gemessen" in eng.messages[-1]["content"]


def test_a_good_answer_is_returned_untouched():
    eng = _Engine()
    eng._note_ambiguous_columns(_NOTE)
    answer = ("Anschaffungswert ist mehrdeutig — welche Lesart gilt? "
              "Summe wäre 25.136 EUR.")
    flagged = eng._scan_ambiguous_column_totals(answer)
    assert eng._append_ambiguous_column_caveat(answer, flagged) == answer


def test_the_ledger_is_cleared_for_each_user_turn():
    """Carrying it across turns would caveat a later, unrelated figure."""
    import inspect

    from delfin.agent.engine import AgentEngine

    source = inspect.getsource(AgentEngine.run_turn) if hasattr(
        AgentEngine, "run_turn") else inspect.getsource(AgentEngine)
    assert "self._ambiguous_columns_turn = []" in source


def test_a_total_that_never_names_the_column_is_a_known_gap():
    """Requiring the column name keeps an unrelated total from being
    caveated, and the price is a figure that names nothing. The observed
    failure named it; this records the limit rather than implying there
    is none."""
    assert vg.scan_for_totals_over_ambiguous_columns(
        "Gesamtwert: 25.136 EUR.", ["Anschaffungswert"]) == []
