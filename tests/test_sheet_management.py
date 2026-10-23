"""Adding, renaming and removing sheets in a workbook.

A workbook of bookings gets a sheet per month, and until now the only way
to make one was to open the file somewhere else. The tab strip is where
that belongs, so it is shown even when there is only one sheet -- a strip
that appears once there are already two is a strip nobody finds.

Each of these writes the workbook immediately, so they are refused while
edits are pending, the same way paging is: the pending journal addresses
the sheet it was made on.
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
# Names
# ---------------------------------------------------------------------------

def test_a_name_has_to_be_something():
    with pytest.raises(sheet.SpreadsheetError):
        sheet.check_sheet_name('   ')


def test_a_name_excel_would_refuse_is_refused_here(recwarn):
    """Rewriting it silently would leave a sheet the user cannot find by
    the name they typed."""
    for bad in ('Januar/Februar', 'Buchungen[2026]', 'a:b', 'x?', 'y*'):
        with pytest.raises(sheet.SpreadsheetError) as excinfo:
            sheet.check_sheet_name(bad)
        assert 'cannot contain' in str(excinfo.value)


def test_a_name_that_is_too_long_says_the_limit():
    with pytest.raises(sheet.SpreadsheetError) as excinfo:
        sheet.check_sheet_name('x' * 32)
    assert '31' in str(excinfo.value)


def test_a_name_already_in_use_is_refused_case_insensitively():
    with pytest.raises(sheet.SpreadsheetError) as excinfo:
        sheet.check_sheet_name('januar', existing=['Januar'])
    assert 'already' in str(excinfo.value)


def test_a_usable_name_comes_back_trimmed():
    assert sheet.check_sheet_name('  Februar  ', existing=['Januar']) == 'Februar'


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

@pytest.fixture
def book(tmp_path):
    workbook = openpyxl.Workbook()
    workbook.active.title = 'Januar'
    workbook.active['A1'] = 'Betrag'
    workbook.active['A2'] = 12.5
    path = tmp_path / 'Buchungen_2026.xlsx'
    workbook.save(path)
    return path


def test_adding_a_sheet_keeps_the_one_that_was_there(book):
    _backup, name = sheet.add_sheet(book, 'Februar', backup=False)
    assert name == 'Februar'
    reopened = openpyxl.load_workbook(book)
    assert reopened.sheetnames == ['Januar', 'Februar']
    assert reopened['Januar']['A2'].value == 12.5


def test_renaming_keeps_the_contents(book):
    sheet.rename_sheet(book, 'Januar', 'Jaenner', backup=False)
    reopened = openpyxl.load_workbook(book)
    assert reopened.sheetnames == ['Jaenner']
    assert reopened['Jaenner']['A1'].value == 'Betrag'


def test_renaming_to_the_same_name_does_nothing(book):
    before = book.read_bytes()
    made, name = sheet.rename_sheet(book, 'Januar', 'Januar', backup=False)
    assert (made, name) == (None, 'Januar')
    assert book.read_bytes() == before


def test_deleting_says_which_sheet_to_show_next(book):
    sheet.add_sheet(book, 'Februar', backup=False)
    sheet.add_sheet(book, 'Maerz', backup=False)
    _backup, show = sheet.drop_sheet(book, 'Februar', backup=False)
    assert openpyxl.load_workbook(book).sheetnames == ['Januar', 'Maerz']
    assert show == 'Maerz', 'the tab that took its place'


def test_deleting_the_last_one_is_refused(book):
    with pytest.raises(sheet.SpreadsheetError) as excinfo:
        sheet.drop_sheet(book, 'Januar', backup=False)
    assert 'has to keep one' in str(excinfo.value)
    assert openpyxl.load_workbook(book).sheetnames == ['Januar']


def test_a_sheet_that_is_not_there(book):
    for call in (lambda: sheet.rename_sheet(book, 'Nope', 'X', backup=False),
                 lambda: sheet.drop_sheet(book, 'Nope', backup=False)):
        with pytest.raises(sheet.SpreadsheetError) as excinfo:
            call()
        assert 'no sheet called' in str(excinfo.value)


def test_a_copy_is_kept_before_the_workbook_is_changed(book, tmp_path):
    folder = tmp_path / 'Backups'
    made, _name = sheet.add_sheet(book, 'Februar', backup_dir=folder)
    assert made is not None and made.parent == folder
    assert openpyxl.load_workbook(made).sheetnames == ['Januar'], (
        'the copy is of the workbook before the sheet was added')


def test_a_failed_save_leaves_no_half_written_workbook(book, monkeypatch):
    before = book.read_bytes()

    def explode(self, target):
        raise OSError('disk full')

    monkeypatch.setattr(openpyxl.Workbook, 'save', explode)
    with pytest.raises(OSError):
        sheet.add_sheet(book, 'Februar', backup=False)

    assert book.read_bytes() == before
    assert not list(book.parent.glob('.dsheet-*'))


# ---------------------------------------------------------------------------
# The tab strip
# ---------------------------------------------------------------------------

def _grid(names, *, editable=True, kind='xlsx'):
    """The markup only. The stylesheet names every class it can style, so
    searching the whole document would find them all whatever is rendered."""
    data = sheet.SheetData(name=names[0], values=[['a']], total_rows=1,
                           total_cols=1, col_widths=[110])
    markup = sheet.render_grid_html(data, sheet_names=names, token='tok',
                                    kind=kind, editable=editable)
    return markup[markup.index('</style>'):]


def test_the_strip_is_there_for_a_single_sheet_too():
    """That is where the control for making the second one lives."""
    markup = _grid(['Januar'])
    assert 'dsheet-tabs' in markup
    assert 'dsheet-tab-add' in markup


def test_the_last_sheet_offers_no_delete():
    assert 'dsheet-tab-drop' not in _grid(['Januar'])
    assert 'dsheet-tab-drop' in _grid(['Januar', 'Februar'])


def test_a_csv_has_no_sheets():
    assert 'dsheet-tabs' not in _grid(['data.csv'], kind='csv')


def test_a_read_only_grid_offers_no_sheet_changes():
    markup = _grid(['Januar', 'Februar'], editable=False)
    assert 'dsheet-tab-add' not in markup
    assert 'dsheet-tab-drop' not in markup
    assert 'data-sheet="Februar"' in markup, 'switching sheets still works'


def test_deleting_asks_twice():
    """A sheet is not something to lose to a slip."""
    script = sheet.grid_js('calc-scope-1', 'tok')
    body = script[script.index('dropTab.addEventListener'):][:600]
    assert "dataset.armed" in body
    assert "'Delete?'" in body


def test_the_name_is_typed_where_the_name_will_appear():
    script = sheet.grid_js('calc-scope-1', 'tok')
    assert 'dsheet-name-input' in script
    assert "send(action, [], {name: name})" in script


# ---------------------------------------------------------------------------
# In the browser
# ---------------------------------------------------------------------------

@pytest.fixture
def office(tmp_path, book):
    (tmp_path / 'archive').mkdir()
    ctx = DashboardContext(calc_dir=tmp_path, archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path)
    ctx.run_js = lambda _s: None
    widget, refs = browser.create_tab(ctx)
    refs['calc_list_directory']()
    file_list = refs['calc_file_list']
    match = [o for o in file_list.options if 'Buchungen' in str(o)]
    value = match[0][1] if isinstance(match[0], tuple) else match[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value
    return book, refs


def _send(refs, action, **extra):
    state = refs['xyz_batch_state']
    message = {'action': action, 'token': state['sheet_view']['token'],
               'ops': [], 'cols': [], 'scroll': 0, 'cur': [0, 1]}
    message.update(extra)
    refs['calc_sheet_payload_input'].value = json.dumps(message)
    refs['calc_sheet_action_btn'].click()


def test_a_new_sheet_is_added_and_shown(office):
    path, refs = office
    _send(refs, 'new_sheet', name='Februar')
    assert openpyxl.load_workbook(path).sheetnames == ['Januar', 'Februar']
    assert refs['xyz_batch_state']['sheet_view']['sheet'] == 'Februar'


def test_a_rejected_name_says_why_and_changes_nothing(office):
    path, refs = office
    _send(refs, 'new_sheet', name='Januar')
    assert openpyxl.load_workbook(path).sheetnames == ['Januar']
    assert 'already' in refs['calc_file_info'].value


def test_sheets_cannot_be_changed_while_edits_are_pending(office):
    path, refs = office
    _send(refs, 'edit', ops=[{'op': 'set', 'row': 1, 'col': 1, 'text': 'x'}])
    _send(refs, 'new_sheet', name='Februar')
    assert openpyxl.load_workbook(path).sheetnames == ['Januar']
    assert 'Save or discard first' in refs['calc_file_info'].value


def test_the_backup_lands_in_the_office_folder(office):
    path, refs = office
    _send(refs, 'new_sheet', name='Februar')
    assert (path.parent / 'Backups' / 'Buchungen_2026.bak.xlsx').exists()
