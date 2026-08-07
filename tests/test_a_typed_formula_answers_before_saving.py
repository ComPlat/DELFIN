"""Typing =SUM(A1:A5) shows the sum -- without saving, and without a wait.

Formula results came from an engine that models the whole workbook and reads
the *file*: about ten seconds for five hundred formulas, and nothing at all
until the file had been saved.  So the grid showed the formula text and the
user had to save to find out whether it was right.
"""

from __future__ import annotations

import json

import pytest

from delfin.dashboard import spreadsheet_view as sheet_view
from delfin.dashboard.context import DashboardContext

openpyxl = pytest.importorskip("openpyxl")


@pytest.fixture
def opened(tmp_path):
    """A workbook open in the browser tab, with the grid on screen."""
    from delfin.dashboard import tab_calculations_browser as browser

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    book = openpyxl.Workbook()
    ws = book.active
    ws.title = "Data"
    for row, value in enumerate([10, 20, 30, 40], start=1):
        ws.cell(row=row, column=1, value=value)
    path = tmp_path / "calc" / "numbers.xlsx"
    book.save(path)
    sheet_view.forget_sheet_cache()

    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = sent.append
    _widget, refs = browser.create_tab(ctx)
    refs["calc_list_directory"]()
    refs["calc_render_sheet"](path)
    sent.clear()
    return refs, sent, path


def _type(refs, sent, row, col, text):
    """One cell edit, as the browser reports it."""
    view = refs["sheet_state"]["sheet_view"]
    refs["calc_sheet_payload_input"].value = json.dumps({
        "action": "edit",
        "token": view["token"],
        "path": str(view["path"]),
        "sheet": view["sheet"],
        "kind": view["kind"],
        "ops": [{"op": "set", "row": row, "col": col, "text": text}],
    })
    refs["calc_sheet_action_btn"].click()
    return "\n".join(sent)


def test_a_typed_formula_is_worked_out_at_once(opened):
    refs, sent, _path = opened

    pushed = _type(refs, sent, 6, 1, "=SUM(A1:A4)")

    assert "__dsheetShow" in pushed, "no result was pushed to the grid"
    assert '"100"' in pushed or "[6, 1, \"100\"]" in pushed, (
        f"the sum was not 100 in what was pushed: {pushed[-400:]}"
    )


def test_the_file_is_not_touched_by_working_a_formula_out(opened):
    refs, sent, path = opened
    before = path.read_bytes()

    _type(refs, sent, 6, 1, "=SUM(A1:A4)")

    assert path.read_bytes() == before, "the file was written before Save"


def test_a_formula_sees_the_edits_that_are_not_saved_yet(opened):
    refs, sent, _path = opened
    _type(refs, sent, 6, 1, "=SUM(A1:A4)")
    sent.clear()

    pushed = _type(refs, sent, 1, 1, "70")   # 10 -> 70, so the sum is 160

    assert '"160"' in pushed, f"the sum ignored the unsaved edit: {pushed[-400:]}"


def test_an_error_is_reported_as_excel_reports_it(opened):
    refs, sent, _path = opened

    pushed = _type(refs, sent, 6, 1, "=1/0")

    assert "#DIV/0!" in pushed


def test_a_csv_keeps_its_text(tmp_path):
    """A csv holds text; Excel does not work formulas out in one either."""
    from delfin.dashboard import tab_calculations_browser as browser

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    path = tmp_path / "calc" / "table.csv"
    path.write_text("1,2\n3,4\n", encoding="utf-8")
    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = sent.append
    _widget, refs = browser.create_tab(ctx)
    refs["calc_list_directory"]()
    refs["calc_render_sheet"](path)
    sent.clear()

    pushed = _type(refs, sent, 3, 1, "=SUM(A1:A2)")

    assert "__dsheetShow" not in pushed
