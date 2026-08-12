"""An append reported itself verified without reading anything back.

Two defects that compound on the same call.

``_verify_cells`` reads back the cells in ``applied`` and reports
``verified``. Appended rows are never added to ``applied`` — the append
path only increments a counter — so an append-only ``edit_sheet``
returned ``verified: True`` having checked nothing at all. It is the
same shape as a delegate's report accepted without its evidence, and the
same shape as a cross-check that says "checked" after ruling on nothing.

And the row lands at ``max_row + 1``. A German bookkeeping sheet ends
with a total: ``Summe | =SUMME(D2:D250)``. The new booking is written to
row 251, below the total and outside the SUM range. The tool answers
``1 row(s) appended from row 251`` and the workbook's total is now
wrong — a number the user reads off their own file and believes.

Nothing detects a total row anywhere in the module.

The append is still performed: refusing to write below a total would
block ordinary work on any sheet that has one, and the caller may well
mean exactly that. What changes is that the row is read back like every
other write, and that a total row above the insertion point is named, so
the number that has just gone stale is the caller's decision rather than
their discovery.
"""

from __future__ import annotations

import pytest

from delfin.agent import office as O

openpyxl = pytest.importorskip("openpyxl")


def _book(tmp_path, rows, *, name="Buchungen.xlsx"):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "2026-02"
    for row in rows:
        ws.append(row)
    p = tmp_path / name
    wb.save(p)
    return p


_HEADER = ["Beleg", "Kreditor", "Betrag"]
_BOOKINGS = [
    _HEADER,
    ["R-001", "Meier GmbH", 100],
    ["R-002", "Schulze AG", 200],
]
_WITH_TOTAL = _BOOKINGS + [["Summe", None, "=SUM(C2:C3)"]]


def _notes(result) -> str:
    return " ".join(result.get("notes") or [])


# ---------------------------------------------------------------------------
# The append is read back
# ---------------------------------------------------------------------------

def test_an_append_only_edit_is_actually_verified(tmp_path):
    p = _book(tmp_path, _BOOKINGS)
    out = O.edit_sheet(p, append_records=[
        {"Beleg": "R-003", "Kreditor": "Institut", "Betrag": 300}])
    assert out["verified"] is True
    wb = openpyxl.load_workbook(p)
    assert wb.active.cell(row=4, column=1).value == "R-003"


def test_a_value_that_did_not_land_is_not_reported_as_verified(
    tmp_path, monkeypatch,
):
    """The flag has to come from a read-back, not from having tried."""
    p = _book(tmp_path, _BOOKINGS)

    # The real signature is (path, sheet_name, applied). Report every
    # appended cell as not holding its value, which is what a write that
    # silently did not take looks like from the read-back's side.
    monkeypatch.setattr(
        O, "_verify_cells",
        lambda path, sheet_name, applied: [e["cell"] for e in applied])
    out = O.edit_sheet(p, append_records=[
        {"Beleg": "R-003", "Kreditor": "Institut", "Betrag": 300}])
    assert out["verified"] is False


def test_appended_rows_by_position_are_verified_too(tmp_path):
    p = _book(tmp_path, _BOOKINGS)
    out = O.edit_sheet(p, append_rows=[["R-004", "Uni", 400]])
    assert out["verified"] is True


# ---------------------------------------------------------------------------
# A total above the insertion point is named
# ---------------------------------------------------------------------------

def test_a_total_row_above_the_append_is_reported(tmp_path):
    p = _book(tmp_path, _WITH_TOTAL)
    out = O.edit_sheet(p, append_records=[
        {"Beleg": "R-003", "Kreditor": "Institut", "Betrag": 300}])
    assert "summe" in _notes(out).lower() or "total" in _notes(out).lower()


def test_the_note_says_the_new_row_is_outside_it(tmp_path):
    p = _book(tmp_path, _WITH_TOTAL)
    out = O.edit_sheet(p, append_records=[
        {"Beleg": "R-003", "Kreditor": "Institut", "Betrag": 300}])
    text = _notes(out).lower()
    assert "outside" in text or "not included" in text


def test_the_row_is_still_appended(tmp_path):
    """Refusing would block ordinary work on any sheet that has a total."""
    p = _book(tmp_path, _WITH_TOTAL)
    O.edit_sheet(p, append_records=[
        {"Beleg": "R-003", "Kreditor": "Institut", "Betrag": 300}])
    wb = openpyxl.load_workbook(p)
    assert wb.active.cell(row=5, column=1).value == "R-003"


def test_a_sheet_without_a_total_says_nothing_about_one(tmp_path):
    p = _book(tmp_path, _BOOKINGS)
    out = O.edit_sheet(p, append_records=[
        {"Beleg": "R-003", "Kreditor": "Institut", "Betrag": 300}])
    text = _notes(out).lower()
    assert "summe" not in text and "total row" not in text


@pytest.mark.parametrize("label", ["Summe", "SUMME", "Gesamt", "Total",
                                   "Zwischensumme"])
def test_the_usual_german_and_english_labels_are_recognised(tmp_path, label):
    p = _book(tmp_path, _BOOKINGS + [[label, None, "=SUM(C2:C3)"]])
    out = O.edit_sheet(p, append_records=[
        {"Beleg": "R-003", "Kreditor": "Institut", "Betrag": 300}])
    assert _notes(out), label


def test_a_row_that_merely_says_summary_is_not_a_total(tmp_path):
    """A label alone is not enough — the row has to aggregate something."""
    p = _book(tmp_path, _BOOKINGS + [["Summe der Belege folgt", None, None]])
    out = O.edit_sheet(p, append_records=[
        {"Beleg": "R-003", "Kreditor": "Institut", "Betrag": 300}])
    assert "outside" not in _notes(out).lower()
