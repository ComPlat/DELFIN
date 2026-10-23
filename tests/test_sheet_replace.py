"""Replace in a sheet, with the switches Excel's own dialog carries.

Excel's spreadsheet dialog has *match case* and *match entire cell contents* --
not Word's "whole word", which is a different dialog for a different thing.  A
replacement is an ordinary cell edit, so it is undoable, shows as unsaved, and
is written by the same code that writes a typed cell.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import spreadsheet_view as sheet_view


VALUES = [
    ['ligand', 'Ligand', 'ligands'],
    ['10', 'metal', ''],
    ['', '', '=SUM(A1:A2)'],
]
FORMULAS = [
    ['', '', ''],
    ['', '', ''],
    ['', '', '=SUM(A1:A2)'],
]


def _ops(**kwargs):
    return sheet_view.replace_ops(VALUES, FORMULAS, **kwargs)


def test_every_occurrence_is_replaced():
    ops = _ops(term='ligand', replacement='chelate')

    changed = {(op['row'], op['col']): op['text'] for op in ops}
    assert changed[(1, 1)] == 'chelate'
    assert changed[(1, 2)] == 'chelate'          # case-insensitive by default
    assert changed[(1, 3)] == 'chelates'         # inside a longer word


def test_matching_case_replaces_only_what_matches():
    ops = _ops(term='Ligand', replacement='Chelate', match_case=True)

    changed = {(op['row'], op['col']): op['text'] for op in ops}
    assert changed == {(1, 2): 'Chelate'}


def test_entire_cell_contents_is_a_different_question():
    ops = _ops(term='ligand', replacement='chelate', whole_cell=True)

    changed = {(op['row'], op['col']): op['text'] for op in ops}
    assert (1, 3) not in changed, "'ligands' is not the entire cell"
    assert changed[(1, 1)] == 'chelate'
    assert changed[(1, 2)] == 'chelate'


def test_a_formula_is_searched_in_its_formula_text():
    ops = sheet_view.replace_ops(VALUES, FORMULAS, term='SUM', replacement='MAX')

    changed = {(op['row'], op['col']): op['text'] for op in ops}
    assert changed[(3, 3)] == '=MAX(A1:A2)'


def test_the_result_is_an_ordinary_edit_journal():
    ops = _ops(term='ligand', replacement='chelate')

    # it has to survive the same validation a typed edit goes through
    assert sheet_view.validate_ops(ops) == ops
    assert all(op['op'] == 'set' for op in ops)


def test_blank_cells_and_misses_produce_nothing():
    assert _ops(term='tungsten', replacement='x') == []


def test_searching_for_nothing_is_refused():
    with pytest.raises(sheet_view.SpreadsheetError):
        _ops(term='', replacement='x')


def test_replacing_with_nothing_empties_the_cell():
    ops = _ops(term='metal', replacement='')

    changed = {(op['row'], op['col']): op['text'] for op in ops}
    assert changed[(2, 2)] == ''


def test_the_offered_number_formats_are_excels_own_codes():
    codes = dict(sheet_view.NUMBER_FORMATS)
    assert codes['General'] == 'General'
    assert codes['1,234.57'] == '#,##0.00'
    assert codes['Scientific'] == '0.00E+00'
    for _label, code in sheet_view.NUMBER_FORMATS:
        if code:
            assert sheet_view.check_number_format(code) == code


# ---------------------------------------------------------------------------
# through the tab, the way the browser asks
# ---------------------------------------------------------------------------
def test_replace_all_goes_through_the_tab(tmp_path):
    import json

    from delfin.dashboard import tab_calculations_browser as browser
    from delfin.dashboard.context import DashboardContext

    openpyxl = pytest.importorskip("openpyxl")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    book = openpyxl.Workbook()
    ws = book.active
    ws.title = "Data"
    ws.cell(row=1, column=1, value="ligand")
    ws.cell(row=2, column=1, value="Ligand field")
    path = tmp_path / "calc" / "words.xlsx"
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
    before = path.read_bytes()
    sent.clear()

    view = refs["sheet_state"]["sheet_view"]
    refs["calc_sheet_payload_input"].value = json.dumps({
        "action": "replace", "token": view["token"], "path": str(view["path"]),
        "sheet": view["sheet"], "kind": view["kind"], "ops": [],
        "find": "ligand", "repl": "chelate",
    })
    refs["calc_sheet_action_btn"].click()

    pushed = "\n".join(sent)
    assert "__dsheetApply" in pushed, "the grid was never told what changed"
    assert "chelate" in pushed
    assert path.read_bytes() == before, "replace must not write before Save"
    pending = refs["sheet_state"]["sheet_pending"][(str(path), "Data")]
    assert [op["text"] for op in pending] == ["chelate", "chelate field"]
