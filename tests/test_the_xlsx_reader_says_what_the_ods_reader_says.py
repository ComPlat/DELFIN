"""The .ods reader warns about two things the .xlsx reader stays silent on.

Both were understood — and fixed only on the format nobody here uses.

**Hidden and filtered rows.** ``_read_ods`` counts them and says so:
"N row(s) of this sheet are hidden or filtered out. They are included in
the grid above, so the totals here can differ from what the file shows
on screen." The .xlsx branch never reads ``row_dimensions[...].hidden``
or ``auto_filter``. A finance workbook is normally saved with a filter
left on — the state Buchungen_2026.xlsx was in — and the agent's total
then silently includes the rows the user cannot see. The identical file
saved as .ods gets the sentence.

**Merged cells.** ``read_sheet`` reads through
``iter_rows(values_only=True)``, which yields None for every non-anchor
cell of a merged range; ``ws.merged_cells`` is never touched. The .ods
reader keeps ``covered-table-cell`` elements deliberately so nothing
shifts left. (``_fragile_features`` was reported as claiming to cover merges and not
doing so. Reading it: the docstring names them among the things openpyxl
CARRIES THROUGH a rewrite, which is why they are absent from that list.
The reported defect is not there — the reading problem is real and
belongs to ``read_sheet``.)

The German cost-centre layout is the case: column A holds the
Kostenstelle merged down each block, amounts per row in D. The grid
shows the label on the first row of the block and blank beneath, so
"which cost centre is highest" is answered from whichever row happened
to carry the label.

Neither is fixed by filling values in — that would invent data the file
does not contain. Both are fixed by saying what is there, which is what
the .ods side already does.
"""

from __future__ import annotations

import pytest

from delfin.agent import office as O

openpyxl = pytest.importorskip("openpyxl")


def _book(tmp_path, *, rows, hidden=(), merges=(), name="Buchungen.xlsx"):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "2026-02"
    for row in rows:
        ws.append(row)
    for r in hidden:
        ws.row_dimensions[r].hidden = True
    for span in merges:
        ws.merge_cells(span)
    p = tmp_path / name
    wb.save(p)
    return p


_ROWS = [
    ["Kostenstelle", "Beleg", "Betrag"],
    ["KST 4711", "R-001", 100],
    [None, "R-002", 200],
    [None, "R-003", 300],
    ["KST 4712", "R-004", 50],
]


def _notes(result) -> str:
    return " ".join(result.get("notes") or [])


# ---------------------------------------------------------------------------
# Hidden rows
# ---------------------------------------------------------------------------

def test_a_hidden_row_is_named(tmp_path):
    p = _book(tmp_path, rows=_ROWS, hidden=(3,))
    assert "hidden" in _notes(O.read_sheet(p)).lower()


def test_the_count_of_hidden_rows_is_given(tmp_path):
    p = _book(tmp_path, rows=_ROWS, hidden=(3, 4))
    assert "2" in _notes(O.read_sheet(p))


def test_it_says_the_totals_can_differ_from_the_screen(tmp_path):
    p = _book(tmp_path, rows=_ROWS, hidden=(3,))
    text = _notes(O.read_sheet(p)).lower()
    assert "screen" in text or "bildschirm" in text


def test_a_sheet_with_nothing_hidden_says_nothing_about_it(tmp_path):
    p = _book(tmp_path, rows=_ROWS)
    assert "hidden" not in _notes(O.read_sheet(p)).lower()


def test_the_hidden_rows_are_still_in_the_grid(tmp_path):
    """Same contract as the .ods side: they are included, and that is
    exactly why it has to be said."""
    p = _book(tmp_path, rows=_ROWS, hidden=(3,))
    assert "R-002" in O.read_sheet(p)["grid"]


# ---------------------------------------------------------------------------
# Merged cells
# ---------------------------------------------------------------------------

def test_a_merged_range_is_named(tmp_path):
    p = _book(tmp_path, rows=_ROWS, merges=("A2:A4",))
    assert "merge" in _notes(O.read_sheet(p)).lower()


def test_the_merged_range_is_identified(tmp_path):
    p = _book(tmp_path, rows=_ROWS, merges=("A2:A4",))
    assert "A2" in _notes(O.read_sheet(p))


def test_it_warns_that_the_blanks_below_are_not_empty(tmp_path):
    p = _book(tmp_path, rows=_ROWS, merges=("A2:A4",))
    text = _notes(O.read_sheet(p)).lower()
    assert "blank" in text or "empty" in text


def test_a_sheet_without_merges_says_nothing_about_them(tmp_path):
    p = _book(tmp_path, rows=_ROWS)
    assert "merge" not in _notes(O.read_sheet(p)).lower()


def test_the_values_are_not_invented(tmp_path):
    """Filling the merged range in would put data in the file's mouth."""
    p = _book(tmp_path, rows=_ROWS, merges=("A2:A4",))
    grid = O.read_sheet(p)["grid"]
    assert grid.count("KST 4711") == 1


# ---------------------------------------------------------------------------
# The docstring that claimed the coverage
# ---------------------------------------------------------------------------

def test_the_fragile_scan_is_about_rewriting_not_about_reading():
    """Read carefully rather than reported on: the docstring names merges
    among the things openpyxl CARRIES THROUGH a round-trip, which is why
    they are not in the fragile list. The reading problem they cause is a
    different one, and it belongs to read_sheet -- which is where the note
    above was added, not here."""
    doc = (O._fragile_features.__doc__ or "").lower()
    assert "round-trip" in doc
    import inspect
    assert "merged" not in inspect.getsource(O._fragile_features)
