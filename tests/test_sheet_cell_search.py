"""Searching a spreadsheet finds cells, anywhere in the file.

The grid holds one window of rows at a time and the search box used to run
over the text of that window. So a term on row 900 of a 5000-row sheet was
reported as not found, and a term that was found could only be counted --
the offsets addressed a flattened copy, not a cell anyone could be sent to.

The search now walks the file and answers in cells. Going to a hit loads
the window it lives in, switches sheet if the hit is on another one, puts
the cursor on the cell and centres it.
"""

from __future__ import annotations

import json
import re

import pytest

openpyxl = pytest.importorskip('openpyxl')

from delfin.dashboard import spreadsheet_view as sheet
from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


# ---------------------------------------------------------------------------
# The search itself
# ---------------------------------------------------------------------------

@pytest.fixture
def workbook(tmp_path):
    book = openpyxl.Workbook()
    first = book.active
    first.title = 'Preise'
    for row in range(1, 1001):
        first.cell(row=row, column=1, value=f'Artikel {row}')
        first.cell(row=row, column=2, value=row * 2)
    first.cell(row=900, column=3, value='Sonderposten')
    second = book.create_sheet('Kunden')
    second.cell(row=4, column=2, value='Sonderposten Nord')
    path = tmp_path / 'katalog.xlsx'
    book.save(path)
    return path


def test_a_hit_far_past_the_loaded_window_is_found(workbook):
    hits, capped = sheet.search_cells(workbook, 'Sonderposten')
    assert not capped
    assert (hits[0].sheet, hits[0].row, hits[0].col) == ('Preise', 900, 3)


def test_the_search_crosses_sheets(workbook):
    hits, _capped = sheet.search_cells(workbook, 'Sonderposten')
    assert {hit.sheet for hit in hits} == {'Preise', 'Kunden'}


def test_hits_are_labelled_the_way_a_spreadsheet_labels_cells(workbook):
    hits, _capped = sheet.search_cells(workbook, 'Sonderposten')
    assert hits[0].label == 'C900'


def test_numbers_are_searched_as_they_read(workbook):
    hits, _capped = sheet.search_cells(workbook, '1998')
    assert any(hit.col == 2 and hit.row == 999 for hit in hits)


def test_the_search_is_case_insensitive(workbook):
    assert sheet.search_cells(workbook, 'SONDERPOSTEN')[0]


def test_an_empty_term_finds_nothing(workbook):
    assert sheet.search_cells(workbook, '') == ([], False)


def test_the_cap_is_reported_rather_than_hidden(workbook):
    hits, capped = sheet.search_cells(workbook, 'Artikel', limit=10)
    assert len(hits) == 10
    assert capped is True, 'a shortened list that says nothing reads as complete'


def test_csv_is_searched_the_same_way(tmp_path):
    path = tmp_path / 'liste.csv'
    path.write_text('a;b;c\n1;Treffer;3\n', encoding='utf-8')
    hits, _capped = sheet.search_cells(path, 'treffer')
    assert (hits[0].row, hits[0].col, hits[0].label) == (2, 2, 'B2')


# ---------------------------------------------------------------------------
# Windows
# ---------------------------------------------------------------------------

def test_the_window_for_a_row_lines_up_with_the_page_buttons():
    assert sheet.window_for_row(1, 400) == 0
    assert sheet.window_for_row(400, 400) == 0
    assert sheet.window_for_row(401, 400) == 400
    assert sheet.window_for_row(900, 400) == 800


def test_the_window_survives_nonsense():
    assert sheet.window_for_row(0, 400) == 0
    assert sheet.window_for_row(-5, 400) == 0
    assert sheet.window_for_row(10, 0) == 9


def test_the_grid_can_be_sent_to_a_cell():
    script = sheet.grid_js('calc-scope-1', 'tok')
    assert 'wrap.__dsheetGoto = function' in script
    body = script.split('wrap.__dsheetGoto = function')[1][:700]
    # The sheet's own row number, looked up rather than computed: the window
    # starts at an offset and a filter can hide rows in between.
    assert 'dataset.r' in body
    assert 'clientHeight / 2' in body, 'a hit at the edge reads as not found'
    assert 'dsheet-hit' in body


# ---------------------------------------------------------------------------
# In the browser
# ---------------------------------------------------------------------------

def _find(root, name):
    if name in (getattr(root, '_dom_classes', ()) or ()):
        return root
    for child in getattr(root, 'children', ()) or ():
        found = _find(child, name)
        if found is not None:
            return found
    return None


@pytest.fixture
def office(tmp_path, workbook):
    (tmp_path / 'archive').mkdir()
    scripts: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path, archive_dir=tmp_path / 'archive', office_dir=tmp_path)
    ctx.run_js = scripts.append
    widget, refs = browser.create_tab(ctx)
    refs['calc_list_directory']()
    file_list = refs['calc_file_list']
    match = [o for o in file_list.options if 'katalog.xlsx' in str(o)]
    value = match[0][1] if isinstance(match[0], tuple) else match[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value
    return widget, refs, scripts, ctx


def test_opening_a_sheet_switches_the_search_to_cells(office):
    _widget, refs, _scripts, _ctx = office
    assert refs['xyz_batch_state']['search_kind'] == 'sheet'


def test_searching_reports_the_cell_and_the_sheet(office):
    widget, refs, scripts, _ctx = office
    state = refs['xyz_batch_state']
    scripts.clear()
    refs['calc_search_input'].value = 'Sonderposten'

    assert len(state['sheet_hits']) == 2
    label = _find(widget, 'calc-content-area')  # noqa: F841 - keeps the tree alive
    assert 'C900' in refs['calc_search_result'].value
    assert 'Preise' in refs['calc_search_result'].value
    assert '__dsheetGoto(900, 3)' in '\n'.join(scripts)


def test_going_to_a_hit_loads_the_window_it_lives_in(office):
    _widget, refs, _scripts, _ctx = office
    state = refs['xyz_batch_state']
    refs['calc_search_input'].value = 'Sonderposten'
    assert state['sheet_view']['row_offset'] == 800, (
        'the hit is on row 900; the window that holds it has to be loaded')


def test_a_hit_on_another_sheet_switches_to_that_sheet(office):
    """The workbook opens on 'Preise'; the term is only on 'Kunden'."""
    _widget, refs, scripts, _ctx = office
    state = refs['xyz_batch_state']
    assert state['sheet_view']['sheet'] == 'Preise'

    scripts.clear()
    refs['calc_search_input'].value = 'Nord'

    assert state['sheet_hits'][0].sheet == 'Kunden'
    assert state['sheet_view']['sheet'] == 'Kunden', 'the sheet did not follow'
    assert '__dsheetGoto(4, 2)' in '\n'.join(scripts)
    assert 'Kunden · B4' in refs['calc_search_result'].value


def test_nothing_found_says_so(office):
    _widget, refs, _scripts, _ctx = office
    refs['calc_search_input'].value = 'gibtesnicht'
    assert '0 matches' in refs['calc_search_result'].value
    assert refs['xyz_batch_state']['sheet_hits'] == []
