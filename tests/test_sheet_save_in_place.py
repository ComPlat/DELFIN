"""Saving a spreadsheet does not move the spreadsheet.

Saving used to re-render the whole grid. The scroll position was carried
across, but the rebuilt grid seeded its cursor at A1, so the selection was
gone and the view settled wherever the fresh markup measured. That is not
what pressing save in a spreadsheet does.

Now a save changes three things on screen -- the dirty marks, the two
buttons and the status line -- and touches nothing else. Discard still has
to rebuild, because the values on screen are the discarded ones and have to
come back from the file; it puts the user back on the cell they were on.
"""

from __future__ import annotations

import json
import re

import pytest

openpyxl = pytest.importorskip('openpyxl')

from delfin.dashboard import spreadsheet_view as sheet
from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


def _widget_with_class(root, name):
    if name in (getattr(root, '_dom_classes', ()) or ()):
        return root
    for child in getattr(root, 'children', ()) or ():
        found = _widget_with_class(child, name)
        if found is not None:
            return found
    return None


class Browser:
    """One built file browser, driven the way the page drives it."""

    def __init__(self, root, tmp_path):
        self.scripts: list[str] = []
        self.ctx = DashboardContext(
            calc_dir=root,
            archive_dir=tmp_path / 'archive',
            office_dir=root,
        )
        self.ctx.run_js = self.scripts.append
        self.widget, self.refs = browser.create_tab(self.ctx)
        self.content = _widget_with_class(self.widget, 'calc-content-area')
        self.payload = _widget_with_class(self.widget, 'calc-sheet-payload')
        self.action = _widget_with_class(self.widget, 'calc-sheet-action')
        assert self.content is not None
        assert self.payload is not None and self.action is not None

    def open(self, name):
        self.refs['calc_list_directory']()
        file_list = self.refs['calc_file_list']
        match = [o for o in file_list.options if name in str(o)]
        assert match, f'{name} not in {file_list.options}'
        value = match[0][1] if isinstance(match[0], tuple) else match[0]
        file_list.value = (value,) if isinstance(file_list.value, tuple) else value

    @property
    def token(self):
        found = re.search(r'data-token="([^"]+)"', self.content.value)
        assert found, 'no grid on screen'
        return found.group(1)

    def send(self, action, **extra):
        message = {
            'action': action, 'token': self.token,
            'ops': [], 'cols': [], 'scroll': 0, 'cur': [0, 1],
        }
        message.update(extra)
        self.payload.value = json.dumps(message)
        self.action.click()


@pytest.fixture
def book(tmp_path):
    (tmp_path / 'archive').mkdir()
    workbook = openpyxl.Workbook()
    ws = workbook.active
    for row in range(1, 40):
        ws.cell(row=row, column=1, value=f'row {row}')
        ws.cell(row=row, column=2, value=row)
    path = tmp_path / 'prices.xlsx'
    workbook.save(path)
    return path


def test_saving_leaves_the_grid_on_screen_untouched(book, tmp_path):
    page = Browser(book.parent, tmp_path)
    page.open('prices.xlsx')
    before = page.content.value

    page.send('edit', ops=[{'op': 'set', 'row': 3, 'col': 1, 'text': 'changed'}],
              cur=[12, 2], scroll=340)
    page.send('save', cur=[12, 2], scroll=340)

    assert page.content.value == before, (
        'the grid was re-rendered on save, which is what moved the view')
    assert openpyxl.load_workbook(book).active['A3'].value == 'changed'


def test_the_grid_is_told_it_was_saved(book, tmp_path):
    page = Browser(book.parent, tmp_path)
    page.open('prices.xlsx')
    token = page.token
    page.scripts.clear()

    page.send('edit', ops=[{'op': 'set', 'row': 2, 'col': 2, 'text': '7'}])
    page.send('save')

    sent = '\n'.join(page.scripts)
    assert '__dsheetSaved' in sent, 'nothing told the grid its edits are on disk'
    assert token in sent, 'the message went to a grid that is no longer on screen'


def test_the_saved_hook_clears_marks_and_buttons_without_rebuilding():
    controller = sheet.grid_js('calc-scope-1', 'tok')
    assert 'wrap.__dsheetSaved = function' in controller
    body = controller.split('wrap.__dsheetSaved = function')[1][:400]
    assert 'pending = 0' in body
    assert 'dsheet-dirty' in body
    assert 'reflectPending()' in body


def test_discard_rebuilds_but_puts_the_user_back(book, tmp_path):
    page = Browser(book.parent, tmp_path)
    page.open('prices.xlsx')
    before = page.content.value

    page.send('edit', ops=[{'op': 'set', 'row': 5, 'col': 1, 'text': 'oops'}])
    page.send('discard', cur=[9, 2], scroll=210)

    assert page.content.value != before, 'discard must bring the file values back'
    assert 'data-cursor="9,2"' in page.content.value
    assert 'data-scrolltop="210"' in page.content.value
    # And it really was discarded.
    assert openpyxl.load_workbook(book).active['A5'].value == 'row 5'


def test_a_first_render_has_no_cursor_to_restore(book, tmp_path):
    page = Browser(book.parent, tmp_path)
    page.open('prices.xlsx')
    assert 'data-cursor' not in page.content.value


def test_the_grid_seeds_the_cursor_from_the_markup():
    controller = sheet.grid_js('calc-scope-1', 'tok')
    assert "wrap.dataset.cursor" in controller
    # Seeded without revealing: the saved scroll position is applied right
    # after, and revealing first would drag the sheet back to the cursor.
    seed = controller.split('wrap.dataset.cursor')[1][:220]
    assert 'moveTo(' in seed and 'true)' in seed
