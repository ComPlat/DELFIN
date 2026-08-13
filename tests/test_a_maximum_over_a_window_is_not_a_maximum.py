"""The largest value among five rows of twenty-six is not the largest.

Every guarantee in the office layer used to be a sentence in a prompt --
"never extrapolate a total from the first page", "report coverage, not
just results" -- enforced by nothing. Prompt rules do not bind the model
this framework runs on; mechanisms do. The write side already had a real
one (edit_sheet refuses a file it has not read this session, and refuses
again if the mtime moved). The read side had none.

Two holes, both measured against the shipped scanner and the real
workbook before any of this was written:

    flags=0 | Die höchste Kostenstelle ist KST 4711 mit 128.430,55 €.
    flags=0 | Die höchste Buchung beträgt 1.234,50 EUR.

The first failed with an EMPTY ledger, so it was not about coverage: an
extreme was not a claim shape at all, only totals and derived figures
were. The second is the quieter one -- 1234,50 really was in the grid, so
the claim looked backed, and it was the maximum of the five rows read out
of twenty-six. Nothing carried the fact that the read was a window.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import office


BOOK = Path("tests/fixtures/office_workspace/Buchungen_2026.xlsx")


@pytest.fixture
def book() -> str:
    if not BOOK.exists():                       # generated before every run
        from delfin.agent.benchmark_fixtures import ensure_office_fixtures
        ensure_office_fixtures()
    return str(BOOK)


def _turn(label: str, result: dict) -> str:
    token = office.begin_figure_turn(label)
    office.reset_figure_ledger(token)
    office.record_tool_figures("read_document", result, token=token)
    return token


# ---------------------------------------------------------------------------
# The window, carried as data
# ---------------------------------------------------------------------------

def test_a_read_says_which_rows_it_actually_showed(book):
    part = office.read_sheet(book, max_rows=5)
    assert part["window"] == {"first_row": 1, "last_row": 5,
                              "total_rows": 26, "complete": False}
    whole = office.read_sheet(book, max_rows=500)
    assert whole["window"]["complete"] is True


def test_a_page_further_in_is_never_complete(book):
    """Starting at row 10 and reaching the end is still not the table: the
    rows before it were never seen."""
    rest = office.read_sheet(book, start_row=10, max_rows=500)
    assert rest["window"]["first_row"] == 10
    assert rest["window"]["complete"] is False


def test_the_window_costs_nothing_in_the_context(book):
    """The result dict is rendered into a report for the model, never
    handed over as it stands, which is why this field is free. Asserted on
    what the model is shown -- the grid and the notes -- rather than on
    the source text of the renderer: a test that reads code passes on a
    comment and breaks on a rename, and this file has one of those in its
    own history."""
    result = office.read_sheet(book, max_rows=5)
    shown = str(result["grid"]) + " ".join(result["notes"])
    assert "window" not in shown
    assert "first_row" not in shown and "total_rows" not in shown


def test_a_windowed_value_is_marked_and_a_count_is_not(book):
    token = _turn("marks", office.read_sheet(book, max_rows=5))
    ledger = office.figure_ledger(token)
    cells = [f for f in ledger if f.kind == "cell"]
    counts = [f for f in ledger if f.kind == "count"]
    assert cells and all(f.partial for f in cells)
    # `rows` is the file's own total, stated by the reader whatever it
    # showed. Marking it partial would caveat a correct row count.
    assert counts and not any(f.partial for f in counts)
    assert "von 26" in cells[0].window


def test_a_complete_read_marks_nothing(book):
    token = _turn("whole", office.read_sheet(book, max_rows=500))
    assert not any(f.partial for f in office.figure_ledger(token))


# ---------------------------------------------------------------------------
# The claim shape that was missing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("sentence", [
    "Die höchste Buchung beträgt 1.234,50 EUR.",
    "Der größte Posten liegt bei 1.234,50 EUR.",
    "Der teuerste Beleg kostet 1.234,50 EUR.",
    "Maximal wurden 1.234,50 EUR gebucht.",
    "Die niedrigste Buchung liegt bei 1.234,50 EUR.",
    "The largest booking is 1.234,50 EUR.",
])
def test_an_extreme_is_read_as_a_claim_about_the_whole_file(sentence):
    claims = office._claims_in_sentence(sentence, money_only=False)
    assert any(kind == "extremum" for _fig, kind in claims), claims


def test_a_reference_number_is_not_a_figure():
    """A guard that flags a cost-centre number as an invented amount is
    the kind of noise that gets guards switched off. The spelled-out word
    was already excluded; the abbreviation an office actually writes was
    not."""
    claims = office._claims_in_sentence(
        "Die höchste Kostenstelle ist KST 4711 mit 9.000,00 EUR.",
        money_only=False)
    assert "4711" not in [fig for fig, _kind in claims]


# ---------------------------------------------------------------------------
# ... and what the answer is told
# ---------------------------------------------------------------------------

def test_a_maximum_over_a_window_is_flagged_and_the_window_is_named(book):
    token = _turn("over", office.read_sheet(book, max_rows=5))
    caveat = office.figure_coverage_caveat(
        "Die höchste Buchung beträgt 1.234,50 EUR.", token=token)
    assert caveat, "the maximum of a page passed as the maximum of the file"
    assert "Höchst-" in caveat
    assert "Zeilen 1–5 von 26" in caveat
    # A different mistake needs a different repair: the number was really
    # there, so "where does this come from" would send the reader looking
    # for something that exists.
    assert "vollständig" in caveat


def test_the_same_maximum_is_accepted_from_a_complete_read(book):
    """A check that refuses honest work teaches the model to phrase tasks
    so nothing keys on them. This is the half that keeps that from
    happening."""
    token = _turn("honest", office.read_sheet(book, max_rows=500))
    assert office.figure_coverage_caveat(
        "Die höchste Buchung beträgt 1.234,50 EUR.", token=token) == ""


def test_an_invented_maximum_is_still_called_invented(book):
    token = _turn("invented", office.read_sheet(book, max_rows=500))
    caveat = office.figure_coverage_caveat(
        "Die höchste Buchung beträgt 99.999,99 EUR.", token=token)
    assert "stammt nicht aus einem Werkzeug-Ergebnis" in caveat


def test_quoting_a_value_that_was_shown_is_not_an_extreme(book):
    """Only the extreme claim is bounded by coverage. An answer naming a
    figure it read is quoting the tool, which is what the ledger is for."""
    token = _turn("quote", office.read_sheet(book, max_rows=5))
    assert office.figure_coverage_caveat(
        "Ein Beleg über 1.234,50 EUR ist enthalten.", token=token) == ""


def test_a_row_count_survives_a_windowed_read(book):
    token = _turn("count", office.read_sheet(book, max_rows=5))
    assert office.figure_coverage_caveat(
        "Die Datei hat 26 Zeilen.", token=token) == ""


def test_a_hedged_extreme_is_not_flagged(book):
    """"Vermutlich" is how German says not-measured, and caveating the
    honest answer is how a guard stops being read."""
    token = _turn("hedge", office.read_sheet(book, max_rows=5))
    assert office.figure_coverage_caveat(
        "Die höchste Buchung liegt vermutlich bei rund 1.234,50 EUR.",
        token=token) == ""
