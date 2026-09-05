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


# ---------------------------------------------------------------------------
# The ledger is cleared for each user turn -- driven, not read
# ---------------------------------------------------------------------------
# This was:
#
#     source = inspect.getsource(AgentEngine.run_turn) if hasattr(
#         AgentEngine, "run_turn") else inspect.getsource(AgentEngine)
#     assert "self._ambiguous_columns_turn = []" in source
#
# Three independent escape hatches, all three proven:
#
#   1. ``AgentEngine.run_turn`` does not exist, so the ``hasattr`` arm
#      never fired and the whole CLASS source was read instead -- silently.
#   2. The substring also occurs in ``_note_ambiguous_columns`` as
#      ``ledger = self._ambiguous_columns_turn = []``, a lazy init. So the
#      per-turn reset in ``stream_response`` can be deleted outright and
#      the assertion still holds: deleting it left all 20 tests in this
#      file green.
#   3. The ``_Engine`` stub above creates the ledger in its own
#      ``__init__``, so no behavioural test in this file could see the
#      reset missing either.
#
# Two real turns through ``stream_response`` is the shape that cannot be
# satisfied by a rename, a lazy init or a stub.

_ENGINE_NOTE = (
    "NOTE: column 'Anschaffungswert': values like '8.986' read as two "
    "different numbers depending on the convention, and nothing in the "
    "column decides it.")

# The answer the model actually produced on the benchmark: one total over
# the column the reader said it could not read.
_TOTAL_ANSWER = ("Gesamtwert Anschaffungswert: 25.136 EUR. Hinweis: im "
                 "deutschen Zahlenformat interpretiert.")


@pytest.fixture
def agent_tree(tmp_path):
    """The minimum pack an engine needs to build."""
    import textwrap

    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# quick mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def _two_turn_engine(agent_tree, first_turn_events):
    """An engine whose tool results arrive on the FIRST turn only."""
    from unittest.mock import MagicMock, patch

    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    turns = {"n": 0}

    def _stream(*a, **k):
        turns["n"] += 1
        if turns["n"] == 1:
            for event in first_turn_events:
                yield event
        yield StreamEvent(type="text_delta", text=_TOTAL_ANSWER)
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    client = MagicMock()
    client._observed_files_session = set()
    client.stream_message = MagicMock(side_effect=_stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


def _ambiguity_event():
    from delfin.agent.api_client import StreamEvent

    return StreamEvent(
        type="tool_result",
        tool_name="mcp__delfin-docs__read_document",
        tool_output="| Posten | Anschaffungswert |\n" * 40,
        output_truncated=True,
        output_chars=9000,
        output_notes=_ENGINE_NOTE,
    )


def test_a_flagged_column_caveats_the_answer_of_that_turn(agent_tree):
    """The turn that was told the column is undecidable."""
    engine = _two_turn_engine(agent_tree, (_ambiguity_event(),))
    answer = engine.stream_response("Was ist der Gesamtwert?")
    assert engine._ambiguous_columns_turn == ["Anschaffungswert"]
    assert "nicht gemessen" in answer


def test_the_next_turn_does_not_inherit_the_flag(agent_tree):
    """Carrying it across turns caveats a later, unrelated figure -- and a
    caveat on a figure that is fine teaches the reader to skip caveats."""
    engine = _two_turn_engine(agent_tree, (_ambiguity_event(),))
    first = engine.stream_response("Was ist der Gesamtwert?")
    assert "nicht gemessen" in first, "turn one never armed the ledger"

    second = engine.stream_response("Und wie viele Positionen sind es?")
    assert engine._ambiguous_columns_turn == [], (
        "the ledger survived into a turn that was told nothing")
    assert "nicht gemessen" not in second, (
        "an answer of a turn with no undecidable column was caveated anyway")


def test_a_turn_that_flags_nothing_caveats_nothing(agent_tree):
    """The control: without the reader's note there is no caveat at all,
    so the assertion above is about the reset and not about the text."""
    engine = _two_turn_engine(agent_tree, ())
    assert "nicht gemessen" not in engine.stream_response("Gesamtwert?")


def test_a_total_that_never_names_the_column_is_a_known_gap():
    """Requiring the column name keeps an unrelated total from being
    caveated, and the price is a figure that names nothing. The observed
    failure named it; this records the limit rather than implying there
    is none."""
    assert vg.scan_for_totals_over_ambiguous_columns(
        "Gesamtwert: 25.136 EUR.", ["Anschaffungswert"]) == []
