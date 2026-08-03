"""Tests for the editable spreadsheet grid used by the file browser tab."""

import importlib.util
import sys
import zipfile
from pathlib import Path

import pytest

_MODULE_PATH = Path(__file__).resolve().parents[1] / 'delfin' / 'dashboard' / 'spreadsheet_view.py'
_SPEC = importlib.util.spec_from_file_location('_delfin_spreadsheet_view', _MODULE_PATH)
sv = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = sv
_SPEC.loader.exec_module(sv)


# ---------------------------------------------------------------------------
# Addressing
# ---------------------------------------------------------------------------

def test_col_letter_boundaries():
    assert sv.col_letter(0) == 'A'
    assert sv.col_letter(25) == 'Z'
    assert sv.col_letter(26) == 'AA'
    assert sv.col_letter(701) == 'ZZ'
    assert sv.col_letter(702) == 'AAA'


def test_col_letter_index_round_trip():
    for i in range(0, 800):
        assert sv.col_index(sv.col_letter(i)) == i


def test_col_index_rejects_junk():
    for bad in ('', '  ', '7B', 'A1', '-'):
        with pytest.raises(ValueError):
            sv.col_index(bad)


# ---------------------------------------------------------------------------
# Value coercion
# ---------------------------------------------------------------------------

def test_coerce_value_basic_types():
    assert sv.coerce_value('') is None
    assert sv.coerce_value(None) is None
    assert sv.coerce_value('42') == 42
    assert isinstance(sv.coerce_value('42'), int)
    assert sv.coerce_value('3.14') == pytest.approx(3.14)
    assert sv.coerce_value('1e3') == pytest.approx(1000.0)
    assert sv.coerce_value('-0.5') == pytest.approx(-0.5)
    assert sv.coerce_value('TRUE') is True
    assert sv.coerce_value('false') is False


def test_coerce_value_keeps_identifiers_as_text():
    # Leading zeros are sample IDs, not numbers.
    assert sv.coerce_value('007') == '007'
    assert sv.coerce_value('0031721') == '0031721'
    assert sv.coerce_value('X31721') == 'X31721'


def test_coerce_value_rejects_non_finite():
    # NaN and inf cannot be written to a valid workbook.
    for text in ('NaN', 'nan', 'inf', '-inf', 'Infinity'):
        assert isinstance(sv.coerce_value(text), str)


def test_coerce_value_keeps_formulas():
    assert sv.coerce_value('=SUM(A1:A2)') == '=SUM(A1:A2)'


def test_format_cell():
    import datetime as dt

    assert sv.format_cell(None) == ''
    assert sv.format_cell(True) == 'TRUE'
    assert sv.format_cell(5.0) == '5'
    assert sv.format_cell(5.25) == '5.25'
    assert sv.format_cell(dt.date(2026, 5, 5)) == '2026-05-05'
    assert sv.format_cell(dt.datetime(2026, 5, 5, 0, 0)) == '2026-05-05'


# ---------------------------------------------------------------------------
# Journal validation and replay
# ---------------------------------------------------------------------------

def test_validate_ops_accepts_known_ops():
    ops = sv.validate_ops([
        {'op': 'set', 'row': 2, 'col': 3, 'text': 'x'},
        {'op': 'insert_rows', 'at': 5, 'count': 2},
        {'op': 'delete_cols', 'at': 1, 'count': 1},
    ])
    assert [o['op'] for o in ops] == ['set', 'insert_rows', 'delete_cols']


@pytest.mark.parametrize('bad', [
    [{'op': 'drop_table'}],
    [{'op': 'set', 'row': 0, 'col': 1}],
    [{'op': 'set', 'row': 1, 'col': 0}],
    [{'op': 'insert_rows', 'at': 0, 'count': 1}],
    [{'op': 'insert_rows', 'at': 1, 'count': 99999}],
    ['not-a-dict'],
])
def test_validate_ops_rejects_bad_input(bad):
    with pytest.raises(ValueError):
        sv.validate_ops(bad)


def _grid(rows):
    return sv.SheetData(
        name='S',
        values=[list(r) for r in rows],
        formulas=[[''] * len(rows[0]) for _ in rows],
        col_widths=[100] * len(rows[0]),
        row_offset=0,
        total_rows=len(rows),
        total_cols=len(rows[0]),
    )


def test_replay_set_marks_dirty():
    sheet = _grid([['a', 'b'], ['c', 'd']])
    dirty = sv.replay_ops(sheet, sv.validate_ops([{'op': 'set', 'row': 2, 'col': 1, 'text': 'Z'}]))
    assert sheet.values == [['a', 'b'], ['Z', 'd']]
    assert dirty == {(1, 0)}


def test_replay_insert_and_delete_rows_shift_dirty_marks():
    sheet = _grid([['a'], ['b'], ['c']])
    ops = sv.validate_ops([
        {'op': 'set', 'row': 3, 'col': 1, 'text': 'C'},
        {'op': 'insert_rows', 'at': 1, 'count': 1},
    ])
    dirty = sv.replay_ops(sheet, ops)
    assert [r[0] for r in sheet.values] == ['', 'a', 'b', 'C']
    assert dirty == {(3, 0)}


def test_replay_delete_cols():
    sheet = _grid([['a', 'b', 'c'], ['d', 'e', 'f']])
    sv.replay_ops(sheet, sv.validate_ops([{'op': 'delete_cols', 'at': 2, 'count': 1}]))
    assert sheet.values == [['a', 'c'], ['d', 'f']]


def test_replay_respects_row_offset():
    sheet = _grid([['x'], ['y']])
    sheet.row_offset = 100
    sv.replay_ops(sheet, sv.validate_ops([{'op': 'set', 'row': 102, 'col': 1, 'text': 'hit'}]))
    assert sheet.values == [['x'], ['hit']]


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def test_render_escapes_cell_content():
    sheet = _grid([['<script>alert(1)</script>', 'a & b'], ['"quoted"', '']])
    html = sv.render_grid_html(sheet, sheet_names=['<b>S</b>'], token='tok', path='/tmp/x.xlsx')
    assert '<script>alert(1)' not in html
    assert '&lt;script&gt;' in html
    assert 'a &amp; b' in html
    assert 'data-token="tok"' in html


def test_render_keeps_cells_lean():
    # Per-cell styles or data attributes would blow up the widget comm payload.
    sheet = _grid([['a', 'b'], ['c', 'd']])
    html = sv.render_grid_html(sheet, token='t')
    assert html.count('<td') == 4
    assert 'data-r=' in html          # once per row
    assert '<td style=' not in html
    assert 'data-c=' not in html


def test_render_marks_dirty_cells():
    sheet = _grid([['a', 'b']])
    html = sv.render_grid_html(sheet, token='t', dirty={(0, 1)}, pending=1)
    assert 'dsheet-dirty' in html
    # Save is offered because there is something to save.
    assert 'dsheet-save"' in html


def test_render_readonly_hides_editing_controls():
    sheet = _grid([['a']])
    html = sv.render_grid_html(sheet, token='t', editable=False)
    assert 'dsheet-save' not in html
    assert 'data-act="insert_rows"' not in html


def test_grid_to_tsv():
    assert sv.grid_to_tsv(_grid([['a', 'b'], ['c', 'd']])) == 'a\tb\nc\td'


def test_grid_js_embeds_scope_and_token():
    js = sv.grid_js('calc-scope-42', 'tok-9')
    assert '"calc-scope-42"' in js
    assert '"tok-9"' in js


# ---------------------------------------------------------------------------
# Delimited files
# ---------------------------------------------------------------------------

def test_sniff_delimiter():
    assert sv.sniff_delimiter('a;b;c\n1;2;3\n') == ';'
    assert sv.sniff_delimiter('a,b,c\n1,2,3\n') == ','
    assert sv.sniff_delimiter('anything', '.tsv') == '\t'


def test_read_delimited_semicolon(tmp_path):
    path = tmp_path / 'mutations.csv'
    path.write_text('Label;SMILES\nA1;CCO\nA2;CCN\n', encoding='utf-8')
    sheet, delim = sv.read_delimited(path)
    assert delim == ';'
    assert sheet.values[0][:2] == ['Label', 'SMILES']
    assert sheet.values[2][:2] == ['A2', 'CCN']


def test_apply_ops_delimited_round_trip(tmp_path):
    path = tmp_path / 'data.csv'
    path.write_text('Label;SMILES\nA1;CCO\nA2;CCN\n', encoding='utf-8')
    sv.apply_ops_delimited(path, [{'op': 'set', 'row': 2, 'col': 2, 'text': 'OCC'}], ';')
    assert path.read_text(encoding='utf-8').splitlines() == ['Label;SMILES', 'A1;OCC', 'A2;CCN']


def test_apply_ops_delimited_preserves_quoting(tmp_path):
    path = tmp_path / 'data.csv'
    path.write_text('name,note\nx,"has, comma"\n', encoding='utf-8')
    sv.apply_ops_delimited(path, [{'op': 'set', 'row': 2, 'col': 1, 'text': 'y'}], ',')
    assert path.read_text(encoding='utf-8').splitlines() == ['name,note', 'y,"has, comma"']


def test_apply_ops_delimited_leaves_untouched_lines_byte_identical(tmp_path):
    # Re-encoding every line would drop quoting the user never asked us to change.
    path = tmp_path / 'data.csv'
    original = 'Label;SMILES;Kommentar\nA1;CCO;Ethanol\nA2;CCN;"hat, Komma"\nA3; spaced ;x\n'
    path.write_text(original, encoding='utf-8')
    sv.apply_ops_delimited(path, [{'op': 'set', 'row': 2, 'col': 2, 'text': 'OCC'}], ';')
    lines = path.read_text(encoding='utf-8').splitlines()
    assert lines[0] == 'Label;SMILES;Kommentar'
    assert lines[1] == 'A1;OCC;Ethanol'      # the one line we edited
    assert lines[2] == 'A2;CCN;"hat, Komma"'  # quoting kept
    assert lines[3] == 'A3; spaced ;x'        # spacing kept


def test_apply_ops_delimited_keeps_crlf(tmp_path):
    path = tmp_path / 'data.csv'
    path.write_bytes(b'a,b\r\n1,2\r\n')
    sv.apply_ops_delimited(path, [{'op': 'set', 'row': 2, 'col': 1, 'text': '9'}], ',')
    assert path.read_bytes() == b'a,b\r\n9,2\r\n'


def test_apply_ops_delimited_does_not_grow_on_blank_edits(tmp_path):
    path = tmp_path / 'data.csv'
    path.write_text('a,b\n1,2\n', encoding='utf-8')
    # Typing into a blank row past the end and clearing it again must not add lines.
    sv.apply_ops_delimited(path, [{'op': 'set', 'row': 2, 'col': 1, 'text': '1'}], ',')
    assert path.read_text(encoding='utf-8') == 'a,b\n1,2\n'


def test_apply_ops_delimited_structural(tmp_path):
    path = tmp_path / 'data.csv'
    path.write_text('a,b\n1,2\n3,4\n', encoding='utf-8')
    sv.apply_ops_delimited(path, [{'op': 'delete_rows', 'at': 2, 'count': 1}], ',')
    assert path.read_text(encoding='utf-8').splitlines() == ['a,b', '3,4']


# ---------------------------------------------------------------------------
# Backups
# ---------------------------------------------------------------------------

def test_backup_created_once(tmp_path):
    path = tmp_path / 'sheet.csv'
    path.write_text('a,b\n', encoding='utf-8')
    first = sv.make_backup(path)
    assert first is not None and first.name == 'sheet.bak.csv'
    original_backup = first.read_text(encoding='utf-8')

    path.write_text('changed,now\n', encoding='utf-8')
    second = sv.make_backup(path)
    assert second is None
    # The first-ever state is what a backup is for; it must not be overwritten.
    assert first.read_text(encoding='utf-8') == original_backup


# ---------------------------------------------------------------------------
# Workbook feature inspection (zipfile only, no openpyxl needed)
# ---------------------------------------------------------------------------

def _fake_xlsx(path, extra_parts=()):
    with zipfile.ZipFile(str(path), 'w') as zf:
        zf.writestr('[Content_Types].xml', '<Types/>')
        zf.writestr('xl/workbook.xml', '<workbook/>')
        for part in extra_parts:
            zf.writestr(part, 'x')


def test_inspect_workbook_features_clean(tmp_path):
    path = tmp_path / 'clean.xlsx'
    _fake_xlsx(path)
    assert sv.inspect_workbook_features(path) == []
    assert sv.describe_lossy_features([]) == ''


def test_inspect_workbook_features_detects_charts_and_media(tmp_path):
    path = tmp_path / 'rich.xlsx'
    _fake_xlsx(path, ['xl/charts/chart1.xml', 'xl/media/image1.png'])
    found = sv.inspect_workbook_features(path)
    assert 'charts' in found
    assert 'embedded images' in found
    note = sv.describe_lossy_features(found)
    assert 'charts' in note


def test_inspect_workbook_features_survives_garbage(tmp_path):
    path = tmp_path / 'broken.xlsx'
    path.write_bytes(b'not a zip')
    assert sv.inspect_workbook_features(path) == []


# ---------------------------------------------------------------------------
# Workbook round trip
# ---------------------------------------------------------------------------

openpyxl = pytest.importorskip('openpyxl')


def _make_workbook(path, rows, formulas=None, sheet_title='Tabelle1'):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet_title
    for row in rows:
        ws.append(list(row))
    for addr, formula in (formulas or {}).items():
        ws[addr] = formula
    wb.save(str(path))
    wb.close()


def test_read_xlsx_values_and_dimensions(tmp_path):
    path = tmp_path / 'book.xlsx'
    _make_workbook(path, [['Label', 'SMILES'], ['A1', 'CCO'], ['A2', 'CCN']])
    names, sheet = sv.read_xlsx(path)
    assert names == ['Tabelle1']
    assert sheet.name == 'Tabelle1'
    assert sheet.values[0][:2] == ['Label', 'SMILES']
    assert sheet.values[2][:2] == ['A2', 'CCN']
    # Blank rows are appended so the user can type past the end.
    assert sheet.n_rows > 3


def test_read_xlsx_shows_formula_when_no_cached_value(tmp_path):
    path = tmp_path / 'formula.xlsx'
    _make_workbook(path, [[1, 2]], formulas={'C1': '=A1+B1'})
    _, sheet = sv.read_xlsx(path)
    assert sheet.has_formulas
    assert sheet.formulas[0][2] == '=A1+B1'
    # openpyxl-written files carry no cached result, so the formula is displayed.
    assert sheet.values[0][2] == '=A1+B1'


def test_read_xlsx_row_window(tmp_path):
    path = tmp_path / 'big.xlsx'
    _make_workbook(path, [[f'r{i}'] for i in range(1, 60)])
    _, first = sv.read_xlsx(path, row_offset=0)
    assert first.row_offset == 0
    assert first.values[0][0] == 'r1'
    _, second = sv.read_xlsx(path, row_offset=10)
    assert second.row_offset == 10
    assert second.values[0][0] == 'r11'


def test_apply_ops_xlsx_round_trip(tmp_path):
    path = tmp_path / 'book.xlsx'
    _make_workbook(path, [['Label', 'SMILES'], ['A1', 'CCO']])
    backup = sv.apply_ops_xlsx(path, 'Tabelle1', [
        {'op': 'set', 'row': 2, 'col': 2, 'text': 'OCC'},
        {'op': 'set', 'row': 3, 'col': 1, 'text': '007'},
        {'op': 'set', 'row': 3, 'col': 2, 'text': '42'},
    ])
    assert backup is not None and backup.exists()

    wb = openpyxl.load_workbook(str(path))
    ws = wb['Tabelle1']
    assert ws['B2'].value == 'OCC'
    assert ws['A3'].value == '007'      # stayed text
    assert ws['B3'].value == 42         # became a number
    wb.close()


def test_apply_ops_xlsx_preserves_untouched_formulas(tmp_path):
    """The save pass must never use data_only -- that would bake in cached values."""
    path = tmp_path / 'formula.xlsx'
    _make_workbook(path, [[1, 2]], formulas={'C1': '=A1+B1'})
    sv.apply_ops_xlsx(path, 'Tabelle1', [{'op': 'set', 'row': 1, 'col': 1, 'text': '5'}])

    wb = openpyxl.load_workbook(str(path))
    ws = wb['Tabelle1']
    assert ws['A1'].value == 5
    assert ws['C1'].value == '=A1+B1'
    wb.close()


def test_apply_ops_xlsx_writes_formula(tmp_path):
    path = tmp_path / 'book.xlsx'
    _make_workbook(path, [[1, 2]])
    sv.apply_ops_xlsx(path, 'Tabelle1', [{'op': 'set', 'row': 1, 'col': 3, 'text': '=A1+B1'}])
    wb = openpyxl.load_workbook(str(path))
    assert wb['Tabelle1']['C1'].value == '=A1+B1'
    wb.close()


def test_apply_ops_xlsx_structural(tmp_path):
    path = tmp_path / 'book.xlsx'
    _make_workbook(path, [['a'], ['b'], ['c']])
    sv.apply_ops_xlsx(path, 'Tabelle1', [{'op': 'delete_rows', 'at': 2, 'count': 1}])
    wb = openpyxl.load_workbook(str(path))
    ws = wb['Tabelle1']
    assert [ws.cell(row=r, column=1).value for r in (1, 2)] == ['a', 'c']
    wb.close()


def test_apply_ops_xlsx_leaves_no_temp_files(tmp_path):
    path = tmp_path / 'book.xlsx'
    _make_workbook(path, [['a']])
    sv.apply_ops_xlsx(path, 'Tabelle1', [{'op': 'set', 'row': 1, 'col': 1, 'text': 'b'}])
    leftovers = [p.name for p in tmp_path.iterdir() if p.name.startswith('.dsheet-')]
    assert leftovers == []


def test_apply_ops_xlsx_empty_journal_is_a_noop(tmp_path):
    path = tmp_path / 'book.xlsx'
    _make_workbook(path, [['a']])
    assert sv.apply_ops_xlsx(path, 'Tabelle1', []) is None
    assert not sv.backup_path_for(path).exists()


# ---------------------------------------------------------------------------
# Text encoding: a German CSV is not utf-8
# ---------------------------------------------------------------------------

def test_a_cp1252_csv_reads_with_its_umlauts(tmp_path):
    """A CSV exported by a spreadsheet program on a German-language Windows
    system is cp1252. Read as utf-8 with errors="replace" it does not fail;
    it turns 'Grünwald' into 'Gr�nwald', and every name and street in
    a German data set carries one of those bytes."""
    path = tmp_path / 'kunden.csv'
    path.write_bytes(
        'Firma;Ort\nElektro Grünwald;Köln\nMüller & Söhne;Straße 5\n'
        .encode('cp1252'))

    sheet, _delimiter = sv.read_delimited(path)
    assert [row[:2] for row in sheet.values[:3]] == [
        ['Firma', 'Ort'],
        ['Elektro Grünwald', 'Köln'],
        ['Müller & Söhne', 'Straße 5'],
    ]


def test_a_utf8_csv_still_reads_as_utf8(tmp_path):
    """cp1252 maps every byte, so a detector that reaches it too eagerly
    would turn 'Grünwald' into 'GrÃ¼nwald' instead."""
    path = tmp_path / 'kunden.csv'
    path.write_text('Firma\nElektro Grünwald\n', encoding='utf-8')
    sheet, _delimiter = sv.read_delimited(path)
    assert sheet.values[1][0] == 'Elektro Grünwald'


def test_a_utf8_bom_is_not_shown_as_text(tmp_path):
    path = tmp_path / 'kunden.csv'
    path.write_text('Firma\nGrünwald\n', encoding='utf-8-sig')
    sheet, _delimiter = sv.read_delimited(path)
    assert sheet.values[0][0] == 'Firma'


def test_saving_keeps_the_file_in_the_encoding_it_had(tmp_path):
    """Writing utf-8 over a cp1252 file changes every umlaut in it,
    including on the rows the edit never touched."""
    path = tmp_path / 'kunden.csv'
    path.write_bytes(
        'Firma;Ort\r\nElektro Grünwald;Köln\r\nMüller;Straße 5\r\n'
        .encode('cp1252'))

    sv.apply_ops_delimited(
        path, [{'op': 'set', 'row': 2, 'col': 1, 'text': 'Elektro Grünwald GmbH'}],
        ';', backup=False)

    raw = path.read_bytes()
    assert b'\xef\xbf\xbd' not in raw, 'the umlauts were replaced in the file'
    assert raw.decode('cp1252').splitlines() == [
        'Firma;Ort', 'Elektro Grünwald GmbH;Köln', 'Müller;Straße 5']
    assert b'\r\n' in raw, 'the line endings changed too'


def test_an_untouched_row_of_a_cp1252_file_is_byte_identical(tmp_path):
    path = tmp_path / 'kunden.csv'
    path.write_bytes('a;b\nMüller;"x;y"\nSöhne;z\n'.encode('cp1252'))
    before = path.read_bytes().splitlines()[1]

    sv.apply_ops_delimited(
        path, [{'op': 'set', 'row': 3, 'col': 1, 'text': 'Andere'}], ';',
        backup=False)

    assert path.read_bytes().splitlines()[1] == before


def test_a_character_the_file_cannot_hold_is_refused_with_a_reason(tmp_path):
    """Losing the character, or silently rewriting the whole file in
    another encoding, are both worse than saying so."""
    path = tmp_path / 'kunden.csv'
    path.write_bytes('a\nMüller\n'.encode('cp1252'))
    before = path.read_bytes()

    with pytest.raises(sv.SpreadsheetError) as excinfo:
        sv.apply_ops_delimited(
            path, [{'op': 'set', 'row': 2, 'col': 1, 'text': 'Preis 5 ₹'}],
            ';', backup=False)

    assert 'cp1252' in str(excinfo.value)
    assert 'UTF-8' in str(excinfo.value)
    assert path.read_bytes() == before, 'the file was changed anyway'


def test_the_cell_search_reads_the_same_bytes_the_grid_does(tmp_path):
    path = tmp_path / 'kunden.csv'
    path.write_bytes('Firma\nElektro Grünwald\n'.encode('cp1252'))
    hits, _capped = sv.search_cells(path, 'grünwald')
    assert len(hits) == 1
    assert hits[0].text == 'Elektro Grünwald'


def test_xlsx_umlauts_were_never_the_problem(tmp_path):
    """A workbook stores its text as utf-8 XML; openpyxl reads it. Pinned
    so a future encoding change cannot quietly break the format that
    works."""
    workbook = openpyxl.Workbook()
    workbook.active['A1'] = 'Elektro Grünwald'
    workbook.active['B1'] = 'Straße'
    path = tmp_path / 'kunden.xlsx'
    workbook.save(path)

    _names, sheet = sv.read_xlsx(path, None)
    assert sheet.values[0][:2] == ['Elektro Grünwald', 'Straße']


def test_a_file_with_a_bom_keeps_it_and_one_without_does_not_gain_one(tmp_path):
    """utf-8-sig decodes a plain utf-8 file happily, but encoding with it
    writes a byte-order mark -- so a detector that reports it either way
    puts one at the top of every file that is saved."""
    plain = tmp_path / 'plain.csv'
    plain.write_text('a;b\nx;y\n', encoding='utf-8')
    marked = tmp_path / 'marked.csv'
    marked.write_text('a;b\nx;y\n', encoding='utf-8-sig')

    edit = [{'op': 'set', 'row': 2, 'col': 1, 'text': 'z'}]
    sv.apply_ops_delimited(plain, edit, ';', backup=False)
    sv.apply_ops_delimited(marked, list(edit), ';', backup=False)

    assert not plain.read_bytes().startswith(b'\xef\xbb\xbf')
    assert marked.read_bytes().startswith(b'\xef\xbb\xbf')
