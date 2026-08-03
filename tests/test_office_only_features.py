"""The document behaviour is Office's; the calculations browser is untouched.

Selecting a column without travelling to the last row, saving without
rebuilding the view, scrolling a PDF instead of paging it, fullscreen and
editing Word all change how a file behaves. They belong to the Office tab.
The calculations browser is to work exactly as it did, so every one of them
is behind the same switch and these tests check both sides of it.

Two things are deliberately *not* behind it, because they are repairs
rather than changes: the startup script reaching every tab, and each tab
resolving its own DOM. The calculations browser was always the first
instance on the page, so both were already true for it.
"""

from __future__ import annotations

import json
import re

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
    def __init__(self, tmp_path, *, office):
        for name in ('calc', 'archive', 'office'):
            (tmp_path / name).mkdir(exist_ok=True)
        root = tmp_path / ('office' if office else 'calc')
        workbook = openpyxl.Workbook()
        for row in range(1, 20):
            workbook.active.cell(row=row, column=1, value=f'r{row}')
        workbook.save(root / 'liste.xlsx')

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
    return Tab(tmp_path, office=True)


@pytest.fixture
def calculations(tmp_path):
    return Tab(tmp_path, office=False)


# ---------------------------------------------------------------------------
# The switch itself
# ---------------------------------------------------------------------------

def test_the_switch_is_read_from_one_name():
    """Moving a feature across later should be that line, not a hunt."""
    source = browser.__file__
    text = open(source, encoding='utf-8').read()
    assert text.count('_OFFICE_DOC_FEEL = ') == 1
    assert text.count('_OFFICE_DOC_FEEL') > 3, 'the switch is barely used'


# ---------------------------------------------------------------------------
# Spreadsheets
# ---------------------------------------------------------------------------

def test_office_saves_in_place_and_calculations_re_renders(office, calculations):
    for tab, rebuilt in ((office, False), (calculations, True)):
        tab.open('liste.xlsx')
        before = tab.content.value
        tab.send('edit', ops=[{'op': 'set', 'row': 2, 'col': 1, 'text': 'neu'}])
        tab.scripts.clear()
        tab.send('save')
        changed = tab.content.value != before
        assert changed is rebuilt, (
            'Office must not rebuild the grid on save; Calculations must '
            'keep doing exactly what it did')
        if not rebuilt:
            assert '__dsheetSaved' in '\n'.join(tab.scripts)


def test_only_office_holds_the_viewport_when_a_column_is_picked(
        office, calculations):
    office.open('liste.xlsx')
    calculations.open('liste.xlsx')
    assert 'data-office="1"' in office.content.value
    assert 'data-office="0"' in calculations.content.value


def test_the_grid_script_keeps_both_behaviours(office):
    from delfin.dashboard import spreadsheet_view as sheet

    script = sheet.grid_js('calc-scope-1', 'tok')
    assert "var OFFICE = wrap.dataset.office === '1'" in script
    assert 'if (OFFICE)' in script


# ---------------------------------------------------------------------------
# Fullscreen
# ---------------------------------------------------------------------------

def test_fullscreen_is_offered_in_office_only(office, calculations):
    assert _find(office.widget, 'calc-fullscreen-btn') is not None
    assert _find(calculations.widget, 'calc-fullscreen-btn') is None


# ---------------------------------------------------------------------------
# PDF
# ---------------------------------------------------------------------------

def test_the_pdf_panel_scrolls_in_office_and_pages_elsewhere(office, calculations):
    from delfin.dashboard import pdf_view

    made = {}

    class Spy(pdf_view.PdfPanel):
        def __init__(self, **kwargs):
            made.update(kwargs)
            super().__init__(**kwargs)

    for tab, expected in ((office, True), (calculations, False)):
        made.clear()
        original = pdf_view.PdfPanel
        pdf_view.PdfPanel = Spy
        try:
            tab.refs['xyz_batch_state']['pdf_panel'] = None
            (tab.ctx.calc_dir / 'x.pdf').write_bytes(b'%PDF-1.4 broken')
            tab.open('x.pdf')
        finally:
            pdf_view.PdfPanel = original
        assert made.get('continuous') is expected
        assert (made.get('height_px') is None) is expected


# ---------------------------------------------------------------------------
# Word
# ---------------------------------------------------------------------------

def test_word_is_editable_in_office_and_a_preview_elsewhere(office, calculations):
    docx = pytest.importorskip('docx')
    for tab in (office, calculations):
        document = docx.Document()
        document.add_paragraph('Ein Absatz.')
        document.save(tab.ctx.calc_dir / 'brief.docx')
        tab.open('brief.docx')

    assert 'data-a="p:0"' in office.content.value, 'Office lost the addresses'
    assert 'data-a=' not in calculations.content.value, (
        'the calculations browser must keep the preview it had')


def test_the_search_only_changes_kind_where_the_view_did(office, calculations):
    docx = pytest.importorskip('docx')
    for tab in (office, calculations):
        document = docx.Document()
        document.add_paragraph('Rechnung')
        document.save(tab.ctx.calc_dir / 'brief.docx')
        tab.open('brief.docx')

    assert office.refs['xyz_batch_state']['search_kind'] == 'docx'
    assert calculations.refs['xyz_batch_state']['search_kind'] == 'text'
