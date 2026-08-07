"""Fullscreen for a spreadsheet or a document: the document alone.

The tab's own fullscreen gives the tab the window and keeps the explorer beside
what is open, which is right when the point is a job's output next to the folder
it came from.  It is not right for a table.
"""

from __future__ import annotations

import json

import pytest

from delfin.dashboard import docx_view, spreadsheet_view as sheet_view
from delfin.dashboard.context import DashboardContext

openpyxl = pytest.importorskip("openpyxl")


@pytest.fixture
def opened(tmp_path):
    from delfin.dashboard import tab_calculations_browser as browser

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    book = openpyxl.Workbook()
    book.active.cell(row=1, column=1, value=1)
    path = tmp_path / "calc" / "table.xlsx"
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
    return refs, sent


def _ask_for_fullscreen(refs):
    view = refs["sheet_state"]["sheet_view"]
    refs["calc_sheet_payload_input"].value = json.dumps({
        "action": "zen", "token": view["token"], "path": str(view["path"]),
        "sheet": view["sheet"], "kind": view["kind"], "ops": [],
    })
    refs["calc_sheet_action_btn"].click()


def test_the_table_can_take_the_window(opened):
    refs, sent = opened

    _ask_for_fullscreen(refs)

    script = "\n".join(sent)
    assert "calc-zen-doc" in script
    assert "true" in script
    assert refs["sheet_state"]["zen_doc"] is True


def test_asking_again_gives_the_explorer_back(opened):
    refs, sent = opened
    _ask_for_fullscreen(refs)
    sent.clear()

    _ask_for_fullscreen(refs)

    assert refs["sheet_state"]["zen_doc"] is False
    assert "calc-zen-doc" in "\n".join(sent)


def test_the_explorer_is_what_gets_hidden(opened):
    refs, _sent = opened
    css = "\n".join(refs["sheet_state"].get("css_probe", []) or [])
    del css
    from delfin.dashboard import tab_calculations_browser as browser

    source = open(browser.__file__, encoding="utf-8").read()
    assert ".calc-tab.calc-zen-doc .calc-left { display:none !important; }" in source


def test_both_documents_offer_it_in_their_own_bar():
    grid = sheet_view.render_grid_html(
        sheet_view.SheetData(name="S", values=[["1"]], col_widths=[80]),
        token="t", kind="xlsx", editable=True,
    )
    assert "dsheet-zen" in grid
    assert "dsheet-zen" in sheet_view.grid_js("calc-scope-1", "t")

    assert "dw-zen" in docx_view.toolbar_html()
    assert "dw-zen" in docx_view.edit_js("calc-scope-1")


def test_escape_leaves_it():
    assert "Escape" in sheet_view.grid_js("calc-scope-1", "t")
    assert "calc-zen-doc" in sheet_view.grid_js("calc-scope-1", "t")
    assert "calc-zen-doc" in docx_view.edit_js("calc-scope-1")
