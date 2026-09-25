"""Reading a document is the same everywhere; the tabs differ in their job.

This file used to say the opposite -- that scrolling a PDF, filling the window,
editing a Word file and holding the viewport when a column is picked belonged
to the Office tab, and that the calculations browser was to work exactly as it
did. That was a decision about where a feature had landed, not about what a
folder is: none of those behaviours depends on whether the file underneath is a
calculation, a document or something in a home directory. They are the same in
all four explorers now.

What a tab *is* still differs, and that is what the second half holds: the
Office tab keeps a copy before it overwrites a document, the Archive tab moves
things back towards the calculations instead of into the archive, the
Calculations tab is the one that knows about jobs, and the Home tab starts
without the dot folders and does not walk a home directory looking for
workspaces.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext

openpyxl = pytest.importorskip('openpyxl')


def _find(root, name):
    if name in (getattr(root, '_dom_classes', ()) or ()):
        return root
    for child in getattr(root, 'children', ()) or ():
        found = _find(child, name)
        if found is not None:
            return found
    return None


class Tab:
    """One clone of the browser, on the folder the tab it stands for uses."""

    def __init__(self, tmp_path, *, kind='calculations'):
        for name in ('calc', 'archive', 'office'):
            (tmp_path / name).mkdir(exist_ok=True)
        roots = {'calculations': tmp_path / 'calc',
                 'office': tmp_path / 'office',
                 'archive': tmp_path / 'archive'}
        root = roots[kind]
        workbook = openpyxl.Workbook()
        for row in range(1, 20):
            workbook.active.cell(row=row, column=1, value=f'r{row}')
        workbook.save(root / 'liste.xlsx')

        self.kind = kind
        self.root = root
        self.scripts: list[str] = []
        self.ctx = DashboardContext(
            calc_dir=root,
            archive_dir=tmp_path / 'archive',
            office_dir=tmp_path / 'office',
        )
        self.ctx.run_js = self.scripts.append
        self.widget, self.refs = browser.create_tab(self.ctx)
        self.content = _find(self.widget, 'calc-content-area')

    def open(self, name):
        self.refs['calc_list_directory']()
        file_list = self.refs['calc_file_list']
        match = [o for o in file_list.options if name in str(o)]
        assert match, f'{name} not in {file_list.options}'
        value = match[0][1] if isinstance(match[0], tuple) else match[0]
        file_list.value = (value,) if isinstance(file_list.value, tuple) else value

    def send(self, action, **extra):
        token = re.search(r'data-token="([^"]+)"', self.content.value).group(1)
        payload = {'action': action, 'token': token, 'ops': [], 'cols': [],
                   'scroll': 0, 'cur': [0, 1]}
        payload.update(extra)
        _find(self.widget, 'calc-sheet-payload').value = json.dumps(payload)
        _find(self.widget, 'calc-sheet-action').click()


@pytest.fixture
def office(tmp_path):
    return Tab(tmp_path, kind='office')


@pytest.fixture
def calculations(tmp_path):
    return Tab(tmp_path, kind='calculations')


@pytest.fixture
def archive(tmp_path):
    return Tab(tmp_path, kind='archive')


# ---------------------------------------------------------------------------
# The same everywhere: reading a document
# ---------------------------------------------------------------------------

def test_the_viewing_switch_is_read_from_one_name():
    """Whatever is common is common through one name, not a hunt."""
    text = Path(browser.__file__).read_text(encoding='utf-8')
    assert text.count('_DOC_VIEW_FEEL = ') == 1
    assert text.count('_DOC_VIEW_FEEL') > 3, 'the switch is barely used'


def test_a_sheet_is_saved_in_place_in_every_tab(office, calculations):
    """Rebuilding the grid on save threw the user back to A1."""
    for tab in (office, calculations):
        tab.open('liste.xlsx')
        before = tab.content.value
        tab.send('edit', ops=[{'op': 'set', 'row': 2, 'col': 1, 'text': 'neu'}])
        tab.scripts.clear()
        tab.send('save')
        assert tab.content.value == before, f'{tab.kind} rebuilt the grid'
        assert '__dsheetSaved' in '\n'.join(tab.scripts)


def test_picking_a_column_holds_the_viewport_in_every_tab(office, calculations, archive):
    for tab in (office, calculations, archive):
        tab.open('liste.xlsx')
        assert 'data-office="1"' in tab.content.value, tab.kind


def test_the_grid_script_still_has_the_holding_behaviour(office):
    from delfin.dashboard import spreadsheet_view as sheet

    script = sheet.grid_js('calc-scope-1', 'tok')
    assert "var OFFICE = wrap.dataset.office === '1'" in script
    assert 'if (OFFICE)' in script


def test_a_document_can_fill_the_window_in_every_tab(office, calculations, archive):
    for tab in (office, calculations, archive):
        assert _find(tab.widget, 'calc-fullscreen-btn') is not None, tab.kind


def test_a_pdf_scrolls_through_in_every_tab(office, calculations, archive):
    from delfin.dashboard import pdf_view

    made = {}

    class Spy(pdf_view.PdfPanel):
        def __init__(self, **kwargs):
            made.update(kwargs)
            super().__init__(**kwargs)

    for tab in (office, calculations, archive):
        made.clear()
        original = pdf_view.PdfPanel
        pdf_view.PdfPanel = Spy
        try:
            tab.refs['xyz_batch_state']['pdf_panel'] = None
            (tab.root / 'x.pdf').write_bytes(b'%PDF-1.4 broken')
            tab.open('x.pdf')
        finally:
            pdf_view.PdfPanel = original
        assert made.get('continuous') is True, tab.kind          # page after page
        assert made.get('height_px') is None, tab.kind           # and the pane's height


def test_word_is_the_editable_view_in_every_tab(office, calculations, archive):
    docx = pytest.importorskip('docx')
    for tab in (office, calculations, archive):
        document = docx.Document()
        document.add_paragraph('Ein Absatz.')
        document.save(tab.root / 'brief.docx')
        tab.open('brief.docx')
        assert 'data-a="p:0"' in tab.content.value, (
            f'{tab.kind} lost the paragraph addresses an edit is written back to')
        assert tab.refs['xyz_batch_state']['search_kind'] == 'docx', tab.kind


# ---------------------------------------------------------------------------
# What each tab keeps for itself
# ---------------------------------------------------------------------------

def test_only_office_keeps_a_copy_before_it_overwrites(office, calculations):
    """A Backups folder belongs beside a document, not inside a calculation."""
    for tab in (office, calculations):
        tab.open('liste.xlsx')
        tab.send('edit', ops=[{'op': 'set', 'row': 3, 'col': 1, 'text': 'x'}])
        tab.send('save')

    assert (office.root / 'Backups').is_dir()
    assert not (calculations.root / 'Backups').exists()


def test_the_archive_tab_moves_the_other_way(archive, calculations):
    """Into the archive from the calculations; back out of it from the archive."""
    assert _find(calculations.widget, 'calc-move-archive-btn') is not None
    text = Path(browser.__file__).read_text(encoding='utf-8')
    assert "'archive' if _is_archive_tab else 'calculations'" in text
    assert text.count('_is_archive_tab') > 10, 'the archive tab lost its own behaviour'


def test_the_calculations_tab_is_the_one_that_knows_about_jobs():
    text = Path(browser.__file__).read_text(encoding='utf-8')
    for mark in ("folder_icon = '✅'",            # a finished run
                 "folder_icon = '🔵'",            # one that is running
                 'calc_report_btn',               # and can be reported on
                 'calc_collect_report_targets'):
        assert mark in text, mark


def test_the_home_tab_starts_quiet_and_does_not_walk_the_home(tmp_path, monkeypatch):
    from delfin.dashboard import tab_home

    fake_home = tmp_path / 'home'
    (fake_home / '.cache').mkdir(parents=True)
    (fake_home / 'messwerte').mkdir()
    monkeypatch.setattr(Path, 'home', staticmethod(lambda: fake_home))
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / 'calc', archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda script: None

    _widget, refs = tab_home.create_tab(ctx)
    refs['calc_list_directory']()
    shown = [str(o) for o in refs['calc_file_list'].options]

    assert refs['calc_dotfiles_btn'].value is False
    assert not any('.cache' in o for o in shown)
    assert Path(tab_home.__file__).read_text(encoding='utf-8').count(
        'browser_scan_is_bounded') == 1
