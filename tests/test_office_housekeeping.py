"""What the Office file list shows, and where the copies go.

Saving keeps a copy of the original before overwriting it. Beside every
document that is a file list nobody can read -- one folder of spreadsheets
became .bak.xlsx, .bak.pdf and .bak.docx interleaved with the real files.
In Office the copies go into one folder that can be opened when something
has to come back, and DELFIN's own dot folder is not shown at all.

The calculations browser keeps writing copies beside the file: that is
where its users have always looked for them.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import spreadsheet_view as sheet
from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext

openpyxl = pytest.importorskip('openpyxl')


# ---------------------------------------------------------------------------
# make_backup
# ---------------------------------------------------------------------------

def test_a_copy_lands_beside_the_file_by_default(tmp_path):
    path = tmp_path / 'liste.xlsx'
    path.write_bytes(b'not really a workbook')
    made = sheet.make_backup(path)
    assert made == tmp_path / 'liste.bak.xlsx'
    assert made.exists()


def test_a_folder_takes_the_copy_instead(tmp_path):
    path = tmp_path / 'liste.xlsx'
    path.write_bytes(b'x')
    made = sheet.make_backup(path, folder=tmp_path / 'Backups')
    assert made == tmp_path / 'Backups' / 'liste.bak.xlsx'
    assert made.exists()
    assert not (tmp_path / 'liste.bak.xlsx').exists()


def test_the_copy_is_still_taken_only_once(tmp_path):
    path = tmp_path / 'liste.xlsx'
    path.write_bytes(b'first')
    folder = tmp_path / 'Backups'
    assert sheet.make_backup(path, folder=folder) is not None
    path.write_bytes(b'second')
    assert sheet.make_backup(path, folder=folder) is None
    assert (folder / 'liste.bak.xlsx').read_bytes() == b'first', (
        'the copy is of the original, not of the last save')


def test_an_unwritable_folder_falls_back_beside_the_file(tmp_path, monkeypatch):
    """Losing the copy would be worse than putting it somewhere untidy."""
    path = tmp_path / 'liste.xlsx'
    path.write_bytes(b'x')

    import pathlib
    real_mkdir = pathlib.Path.mkdir

    def refuse(self, *a, **k):
        if self.name == 'Backups':
            raise OSError('read-only')
        return real_mkdir(self, *a, **k)

    monkeypatch.setattr(pathlib.Path, 'mkdir', refuse)
    made = sheet.make_backup(path, folder=tmp_path / 'Backups')
    assert made == tmp_path / 'liste.bak.xlsx'
    assert made.exists()


# ---------------------------------------------------------------------------
# In the browser
# ---------------------------------------------------------------------------

def _tab(tmp_path, *, office):
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir(exist_ok=True)
    root = tmp_path / ('office' if office else 'calc')
    workbook = openpyxl.Workbook()
    workbook.active['A1'] = 'hello'
    workbook.save(root / 'liste.xlsx')
    (root / '.delfin').mkdir(exist_ok=True)
    (root / '.delfin' / 'state.json').write_text('{}', encoding='utf-8')

    ctx = DashboardContext(
        calc_dir=root, archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office')
    ctx.run_js = lambda _s: None
    _widget, refs = browser.create_tab(ctx)
    refs['calc_list_directory']()
    return root, refs


def _labels(refs):
    return [str(o[0] if isinstance(o, tuple) else o).replace(' ', ' ')
            for o in refs['calc_file_list'].options]


def test_office_does_not_list_the_dot_folder(tmp_path):
    _root, refs = _tab(tmp_path, office=True)
    assert not any('.delfin' in label for label in _labels(refs))
    assert any('liste.xlsx' in label for label in _labels(refs))


def test_calculations_still_lists_everything_it_did(tmp_path):
    _root, refs = _tab(tmp_path, office=False)
    assert any('.delfin' in label for label in _labels(refs))


def test_the_dot_folder_is_only_hidden_not_removed(tmp_path):
    root, _refs = _tab(tmp_path, office=True)
    assert (root / '.delfin' / 'state.json').exists()


def test_office_saves_put_the_copy_in_the_backup_folder(tmp_path):
    import json

    root, refs = _tab(tmp_path, office=True)
    file_list = refs['calc_file_list']
    match = [o for o in file_list.options if 'liste.xlsx' in str(o)]
    value = match[0][1] if isinstance(match[0], tuple) else match[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value

    state = refs['xyz_batch_state']
    token = state['sheet_view']['token']
    payload = refs['calc_sheet_payload_input']
    action = refs['calc_sheet_action_btn']
    for message in (
        {'action': 'edit', 'token': token, 'ops': [
            {'op': 'set', 'row': 1, 'col': 1, 'text': 'changed'}],
         'cols': [], 'scroll': 0, 'cur': [0, 1]},
        {'action': 'save', 'token': token, 'ops': [], 'cols': [],
         'scroll': 0, 'cur': [0, 1]},
    ):
        payload.value = json.dumps(message)
        action.click()

    assert (root / 'Backups' / 'liste.bak.xlsx').exists()
    assert not (root / 'liste.bak.xlsx').exists()
    assert openpyxl.load_workbook(root / 'liste.xlsx').active['A1'].value == 'changed'
