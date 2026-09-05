"""Undo and redo in the grid, the way a spreadsheet does it.

An undo is applied as a further change rather than by rewinding the
journal. The file may already have been saved, and undoing after a save is
allowed there and marks the file changed again -- which only works if the
step going back travels the same road as the step that went forward.

Each step therefore carries how to do itself and how to take itself back.
Deletions keep what they removed: an undo that put back empty rows would be
worse than no undo at all.
"""

from __future__ import annotations

import re

import pytest

from delfin.dashboard import spreadsheet_view as sheet


@pytest.fixture
def script():
    return sheet.grid_js('calc-scope-1', 'tok')


def _body(script, name):
    """The source of one function in the controller."""
    start = script.index(name)
    return script[start:start + 1400]


# ---------------------------------------------------------------------------
# The buttons
# ---------------------------------------------------------------------------

def _grid(**kwargs):
    data = sheet.SheetData(
        name='Sheet1', values=[['a', 'b'], ['c', 'd']],
        row_offset=0, total_rows=2, total_cols=2, col_widths=[110, 110],
    )
    return sheet.render_grid_html(data, token='tok', **kwargs)


def test_the_buttons_are_offered_where_the_office_feel_is():
    assert 'dsheet-undo' in _grid(office=True)
    assert 'dsheet-undo' not in _grid(office=False)


def test_they_start_disabled_because_nothing_has_happened_yet():
    markup = _grid(office=True)
    undo = markup[markup.index('dsheet-undo') - 40:markup.index('dsheet-undo') + 60]
    assert 'disabled' in undo


def test_they_say_which_keys_do_the_same(script):
    markup = _grid(office=True)
    assert 'Ctrl+Z' in markup
    assert "e.key === 'z'" in script
    assert 'e.shiftKey ? redo() : undo()' in script


# ---------------------------------------------------------------------------
# How a step is remembered
# ---------------------------------------------------------------------------

def test_a_step_carries_both_directions(script):
    assert 'function undo()' in script and 'function redo()' in script
    remember = _body(script, 'function remember(step)')
    assert 'history.push(step)' in remember
    assert 'future.length = 0' in remember, (
        'a new change has to end the redo branch, or redo would replay a '
        'future that no longer follows from the present')


def test_the_history_is_bounded(script):
    assert 'HISTORY_MAX' in script
    remember = _body(script, 'function remember(step)')
    assert 'history.shift()' in remember


def test_undoing_does_not_record_itself(script):
    """Otherwise the first undo would fill the history with its own
    reflection and the second one would undo the undo."""
    assert 'replaying = true' in script
    remember = _body(script, 'function remember(step)')
    assert 'if (replaying) return' in remember


def test_an_undone_step_still_reaches_python(script):
    """The file may already be saved. Undo has to be a change like any
    other, or the sheet on screen and the file on disk drift apart."""
    write = _body(script, 'function writeCells(entries)')
    assert 'push(' in write
    assert 'before' in write, 'nothing captured what was there'


# ---------------------------------------------------------------------------
# Structural steps
# ---------------------------------------------------------------------------

def test_inserting_rows_is_undone_by_deleting_them(script):
    body = _body(script, 'function insertRows(at, count)')
    assert 'undo: function(){ deleteRows(at, count); }' in body


def test_deleting_rows_keeps_what_it_removed(script):
    body = _body(script, 'function deleteRows(at, count)')
    assert 'removed.push(values)' in body, (
        'an undo that puts back empty rows is worse than no undo')
    assert 'insertRows(at, count)' in body
    assert 'writeCells(writes)' in body


def test_deleting_columns_keeps_what_it_removed(script):
    body = _body(script, 'function deleteCols(at, count)')
    assert 'removedCols.push(row)' in body
    assert 'insertCols(at, count)' in body


def test_a_cell_write_remembers_the_formula_not_the_shown_value(script):
    """A cell showing 42 may hold =6*7. Undo has to put the formula back."""
    body = _body(script, 'function cellText(td)')
    assert "getAttribute('data-f')" in body


# ---------------------------------------------------------------------------
# Every way of writing a cell goes through the history
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('caller', [
    'function commitEdit()',
    'function pasteTsv(text)',
    'function clearSelection()',
])
def test_every_cell_write_is_undoable(script, caller):
    body = _body(script, caller)
    assert 'writeCells(' in body, (
        f'{caller} writes cells without recording a step, so those edits '
        'cannot be undone')
