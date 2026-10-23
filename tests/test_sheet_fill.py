"""Dragging the handle at the corner of a cell to fill down or across.

The difference between entering a formula once and entering it three
hundred times. The rules are worked out in the kernel rather than in the
browser: there is one of them either way, and here they can be read and
tested without a browser.

What the drag adds is written as one step, so one undo takes the whole
fill back rather than one cell of it.
"""

from __future__ import annotations

import json

import pytest

from delfin.dashboard import spreadsheet_view as sheet

openpyxl = pytest.importorskip('openpyxl')


# ---------------------------------------------------------------------------
# Moving a formula with the cell
# ---------------------------------------------------------------------------

def test_a_formula_follows_its_own_row():
    """=B2*C2 dragged down has to become =B3*C3, or every copy reads the
    row it came from."""
    assert sheet.shift_formula('=B2*C2', 1, 0) == '=B3*C3'
    assert sheet.shift_formula('=SUM(A1:A9)', 2, 0) == '=SUM(A3:A11)'


def test_a_locked_reference_stays_where_it_is():
    """Which is the whole reason $ exists."""
    assert sheet.shift_formula('=B$2*C2', 3, 0) == '=B$2*C5'
    assert sheet.shift_formula('=$B2+1', 0, 2) == '=$B2+1'
    assert sheet.shift_formula('=$B$2', 5, 5) == '=$B$2'


def test_columns_move_too():
    assert sheet.shift_formula('=A1', 0, 1) == '=B1'
    assert sheet.shift_formula('=Z1', 0, 1) == '=AA1'


def test_text_in_quotes_is_not_a_reference():
    assert sheet.shift_formula('="A1 bis A9"', 3, 0) == '="A1 bis A9"'
    assert sheet.shift_formula('=A1&" B2"', 1, 0) == '=A2&" B2"'


def test_a_reference_pushed_off_the_sheet_says_so():
    """Excel writes #REF! there. Clamping to A1 would read the wrong cell
    and look like it worked."""
    assert sheet.shift_formula('=A1', -5, 0) == '=#REF!'


def test_something_that_is_not_a_formula_is_left_alone():
    assert sheet.shift_formula('B2', 1, 0) == 'B2'
    assert sheet.shift_formula('', 1, 0) == ''


# ---------------------------------------------------------------------------
# Continuing a series
# ---------------------------------------------------------------------------

def test_one_number_is_copied_not_counted():
    """What a spreadsheet does, and what stops a column of identical
    prices turning into a count."""
    assert sheet.fill_series(['5'], 3, d_row=1) == ['5', '5', '5']


def test_two_numbers_set_the_step():
    assert sheet.fill_series(['1', '2'], 3, d_row=2) == ['3', '4', '5']
    assert sheet.fill_series(['10', '20'], 2, d_row=2) == ['30', '40']
    assert sheet.fill_series(['9', '7'], 2, d_row=2) == ['5', '3']


def test_a_continued_number_is_written_like_the_ones_before_it():
    assert sheet.fill_series(['1,5', '2,0'], 2, d_row=2) == ['2,5', '3,0']
    assert sheet.fill_series(['1.5', '2.0'], 1, d_row=2) == ['2.5']


def test_text_ending_in_digits_counts_on():
    assert sheet.fill_series(['Pos 7'], 3, d_row=1) == ['Pos 8', 'Pos 9', 'Pos 10']


def test_a_padded_number_keeps_its_width():
    assert sheet.fill_series(['A007'], 2, d_row=1) == ['A008', 'A009']


def test_a_single_date_steps_by_a_day_and_over_a_month_end():
    assert sheet.fill_series(['2026-01-30'], 3, d_row=1) == [
        '2026-01-31', '2026-02-01', '2026-02-02']


def test_several_dates_set_their_own_step():
    assert sheet.fill_series(['2026-01-01', '2026-01-08'], 2, d_row=2) == [
        '2026-01-15', '2026-01-22']


def test_anything_else_repeats_the_block_that_was_dragged():
    assert sheet.fill_series(['ja', 'nein'], 5, d_row=2) == [
        'ja', 'nein', 'ja', 'nein', 'ja']


def test_filling_nothing_returns_nothing():
    assert sheet.fill_series([], 3) == []
    assert sheet.fill_series(['x'], 0) == []


# ---------------------------------------------------------------------------
# Blocks
# ---------------------------------------------------------------------------

def test_each_column_of_a_block_continues_on_its_own():
    filled = sheet.fill_block([['=A1', '1'], ['=A2', '2']], rows=2, cols=0)
    assert filled == [['=A3', '3'], ['=A4', '4']]


def test_a_row_can_be_dragged_across():
    assert sheet.fill_block([['1', '2']], rows=0, cols=2) == [['3', '4']]


def test_a_fill_goes_one_way_at_a_time():
    with pytest.raises(sheet.SpreadsheetError):
        sheet.fill_block([['1']], rows=2, cols=2)


# ---------------------------------------------------------------------------
# The grid and the browser
# ---------------------------------------------------------------------------

def test_the_handle_is_only_offered_where_editing_is():
    script = sheet.grid_js('calc-scope-1', 'tok')
    assert 'dsheet-fill' in script
    assert 'if (editable)' in script.split('fillHandle = null')[1][:200]


def test_the_drag_picks_one_axis():
    """Like the handle in a spreadsheet: whichever way the pointer has
    travelled further."""
    script = sheet.grid_js('calc-scope-1', 'tok')
    body = script[script.index('var down = Math.max'):][:300]
    assert 'down >= across' in body


def test_the_whole_fill_is_one_undo_step():
    script = sheet.grid_js('calc-scope-1', 'tok')
    body = script[script.index('wrap.__dsheetApply'):][:600]
    assert 'writeCells(writes)' in body, (
        'one undo has to take the fill back, not one cell of it')


def test_the_grid_sends_what_it_has_on_screen():
    """Including edits that are not saved yet: the kernel would otherwise
    continue the series from the values on disk."""
    script = sheet.grid_js('calc-scope-1', 'tok')
    body = script[script.index("send('fill'"):][:400]
    assert 'block: block' in body
    assert 'at: [' in body


# ---------------------------------------------------------------------------
# End to end through the tab
# ---------------------------------------------------------------------------

@pytest.fixture
def office(tmp_path):
    from delfin.dashboard import tab_calculations_browser as browser
    from delfin.dashboard.context import DashboardContext

    (tmp_path / 'archive').mkdir()
    workbook = openpyxl.Workbook()
    ws = workbook.active
    ws['A1'] = 1
    ws['A2'] = 2
    ws['B1'] = '=A1*2'
    path = tmp_path / 'zahlen.xlsx'
    workbook.save(path)

    scripts: list[str] = []
    ctx = DashboardContext(calc_dir=tmp_path, archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path)
    ctx.run_js = scripts.append
    _widget, refs = browser.create_tab(ctx)
    refs['calc_list_directory']()
    file_list = refs['calc_file_list']
    match = [o for o in file_list.options if 'zahlen' in str(o)]
    value = match[0][1] if isinstance(match[0], tuple) else match[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value
    return path, refs, scripts


def _send(refs, action, **extra):
    state = refs['xyz_batch_state']
    message = {'action': action, 'token': state['sheet_view']['token'],
               'ops': [], 'cols': [], 'scroll': 0, 'cur': [0, 1]}
    message.update(extra)
    refs['calc_sheet_payload_input'].value = json.dumps(message)
    refs['calc_sheet_action_btn'].click()


def test_dragging_two_numbers_down_queues_the_continuation(office):
    path, refs, scripts = office
    scripts.clear()
    _send(refs, 'fill', block=[['1'], ['2']], rows=3, cols=0, at=[1, 1])

    state = refs['xyz_batch_state']
    pending = state['sheet_pending'][(str(path), state['sheet_view']['sheet'])]
    assert [(op['row'], op['text']) for op in pending] == [
        (3, '3'), (4, '4'), (5, '5')]
    assert '__dsheetApply' in '\n'.join(scripts)


def test_the_block_itself_is_not_written_again(office):
    _path, refs, _scripts = office
    _send(refs, 'fill', block=[['1'], ['2']], rows=1, cols=0, at=[1, 1])
    state = refs['xyz_batch_state']
    rows = [op['row'] for op in
            state['sheet_pending'][(str(_path), state['sheet_view']['sheet'])]]
    assert rows == [3], 'rows 1 and 2 already hold what was dragged'


def test_a_formula_dragged_down_moves_with_the_row(office):
    path, refs, _scripts = office
    _send(refs, 'fill', block=[['=A1*2']], rows=2, cols=0, at=[1, 2])
    state = refs['xyz_batch_state']
    pending = state['sheet_pending'][(str(path), state['sheet_view']['sheet'])]
    assert [op['text'] for op in pending] == ['=A2*2', '=A3*2']


def test_a_drag_that_went_nowhere_writes_nothing(office):
    path, refs, _scripts = office
    _send(refs, 'fill', block=[['1']], rows=0, cols=0, at=[1, 1])
    assert not refs['xyz_batch_state']['sheet_pending']


def test_junk_from_the_browser_is_ignored(office):
    _path, refs, _scripts = office
    _send(refs, 'fill', block='not a block', rows=2, cols=0, at=[1, 1])
    _send(refs, 'fill', block=[['1']], rows='many', cols=0, at=[1, 1])
    _send(refs, 'fill', block=[['1']], rows=2, cols=2, at=[1, 1])
    assert not refs['xyz_batch_state']['sheet_pending']


# ---------------------------------------------------------------------------
# Results on screen, formula in the cell
# ---------------------------------------------------------------------------

def test_the_grid_shows_a_result_without_calling_it_an_edit():
    """A worked-out value is not something the user typed: the cell keeps
    its formula, is not marked changed, and undo has nothing to take back."""
    script = sheet.grid_js('calc-scope-1', 'tok')
    body = script[script.index('wrap.__dsheetShow'):][:500]
    assert "getAttribute('data-f')" in body, 'only a formula cell shows a result'
    assert 'writeCells' not in body
    assert 'dsheet-dirty' not in body


def test_entering_a_cell_shows_the_formula_the_way_excel_does():
    """Double-click or F2. This is the behaviour people already have, and
    changing it would be the one thing worth avoiding."""
    script = sheet.grid_js('calc-scope-1', 'tok')
    assert "e.key === 'F2'" in script
    seed = script[script.index('function beginEdit'):][:400]
    assert "getAttribute('data-f')" in seed
    assert 'formula || disp' in seed


def test_results_are_worked_out_again_after_a_save(office):
    """The file changed, so its formulas come to something else now."""
    path, refs, scripts = office
    _send(refs, 'edit', ops=[{'op': 'set', 'row': 1, 'col': 1, 'text': '10'}])
    scripts.clear()
    _send(refs, 'save')

    sent = '\n'.join(scripts)
    assert '__dsheetSaved' in sent
    if pytest.importorskip('formulas', reason='needs the engine'):
        assert '__dsheetShow' in sent, 'the results on screen are now stale'


def test_the_results_are_worked_out_once_per_version_of_the_file(office):
    """Evaluating builds a graph of every cell, so doing it on each paging
    step would make moving through a workbook feel like waiting for one."""
    path, refs, _scripts = office
    state = refs['xyz_batch_state']
    first = state.get('formula_results')
    _send(refs, 'page', dir='next')
    assert state.get('formula_results') is first
