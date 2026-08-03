"""Word documents: one renderer for showing, searching and editing.

The view used to be produced by mammoth. Good HTML, but nothing in it says
which paragraph of the document a line came from, so the search could count
matches and never jump to one, and there was nothing an edit could be
written back to.

Every block now carries its address, and an edit is spliced into the runs
that are already there rather than rebuilding the paragraph -- which is the
difference between editing a sentence and flattening every bold word in it.
"""

from __future__ import annotations

import pytest

docx = pytest.importorskip('docx')

from delfin.dashboard import docx_view as dv


@pytest.fixture
def letter(tmp_path):
    document = docx.Document()
    document.add_heading('Angebot 2026', level=1)
    para = document.add_paragraph('Sehr geehrte Damen, ')
    para.add_run('vielen Dank').bold = True
    para.add_run(' für Ihre Anfrage.')
    document.add_paragraph('Erster Punkt', style='List Bullet')
    table = document.add_table(rows=2, cols=2)
    table.cell(0, 0).text = 'Artikel'
    table.cell(0, 1).text = 'Preis'
    table.cell(1, 0).text = 'Schraube'
    table.cell(1, 1).text = '3,50'
    path = tmp_path / 'angebot.docx'
    document.save(path)
    return path


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def test_every_paragraph_carries_its_address(letter):
    read = dv.read_document(letter)
    addresses = [block.address for block in read.blocks]
    assert 'p:0' in addresses
    assert any(a.startswith('t:0:') for a in addresses), 'table cells are missing'
    assert len(set(addresses)) == len(addresses), 'two blocks share an address'


def test_headings_and_lists_are_recognised(letter):
    read = dv.read_document(letter)
    by_address = {block.address: block for block in read.blocks}
    assert by_address['p:0'].level == 1
    assert by_address['p:2'].listed is True
    assert by_address['p:1'].level == 0


def test_table_cells_know_which_cell_they_are(letter):
    read = dv.read_document(letter)
    cells = [b for b in read.blocks if b.in_table]
    assert cells and all(b.table is not None for b in cells)
    assert {b.text for b in cells} >= {'Artikel', 'Preis', 'Schraube'}


def test_a_doc_file_is_refused_with_a_reason(tmp_path):
    old = tmp_path / 'alt.doc'
    old.write_bytes(b'\xd0\xcf\x11\xe0legacy')
    with pytest.raises(dv.DocxError) as excinfo:
        dv.read_document(old)
    assert 'converted' in str(excinfo.value)


def test_an_oversized_document_is_named_not_truncated(letter, monkeypatch):
    monkeypatch.setattr(dv, 'MAX_FILE_BYTES', 10)
    with pytest.raises(dv.DocxError) as excinfo:
        dv.read_document(letter)
    assert 'MB' in str(excinfo.value)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def test_the_view_and_the_editor_render_the_same_blocks(letter):
    read = dv.read_document(letter)
    view = dv.render_html(read, editable=False)
    edit = dv.render_html(read, editable=True)
    assert view.count('data-a=') == edit.count('data-a=') == len(read.blocks)
    assert 'contenteditable' not in view
    assert edit.count('contenteditable="true"') == len(read.blocks)


def test_cell_text_is_escaped(tmp_path):
    document = docx.Document()
    document.add_paragraph('<script>alert(1)</script>')
    path = tmp_path / 'x.docx'
    document.save(path)
    html = dv.render_html(dv.read_document(path), editable=True)
    assert '<script>' not in html
    assert '&lt;script&gt;' in html


def test_a_truncated_document_is_shown_but_not_edited(letter, monkeypatch):
    monkeypatch.setattr(dv, 'MAX_BLOCKS', 2)
    read = dv.read_document(letter)
    assert len(read.blocks) == 2
    assert not dv.is_editable(read), (
        'editing a truncated view would address paragraphs by an index the '
        'view no longer agrees with')


# ---------------------------------------------------------------------------
# Searching
# ---------------------------------------------------------------------------

def test_search_finds_hits_in_body_and_tables(letter):
    read = dv.read_document(letter)
    assert [h.address for h in dv.search(read, 'preis')] == ['t:0:0:1:0']
    assert dv.search(read, 'dank')[0].address == 'p:1'


def test_search_offsets_point_into_the_block_not_the_document(letter):
    read = dv.read_document(letter)
    hit = dv.search(read, 'Anfrage')[0]
    block = next(b for b in read.blocks if b.address == hit.address)
    assert block.text[hit.start:hit.end] == 'Anfrage'


def test_search_is_case_insensitive_and_finds_repeats(tmp_path):
    document = docx.Document()
    document.add_paragraph('Rechnung und rechnung und RECHNUNG')
    path = tmp_path / 'r.docx'
    document.save(path)
    assert len(dv.search(dv.read_document(path), 'rechnung')) == 3


def test_an_empty_term_finds_nothing(letter):
    assert dv.search(dv.read_document(letter), '') == []


def test_the_jump_wraps_the_match_instead_of_rewriting_the_paragraph():
    """Rewriting innerHTML would flatten the formatting of the paragraph it
    landed in, and in a Word document that is the document."""
    script = dv.focus_js('calc-scope-1', 'p:3', 4, 9)
    assert 'surroundContents' in script
    assert 'innerHTML' not in script
    assert 'scrollIntoView' in script
    assert '"p:3"' in script


# ---------------------------------------------------------------------------
# Editing
# ---------------------------------------------------------------------------

def test_editing_a_word_keeps_the_bold_run_beside_it(letter):
    result = dv.apply_edits(
        letter, {'p:1': 'Sehr geehrte Damen, vielen Dank für Ihre Nachricht.'})
    dv.save(result['document'], letter)

    runs = docx.Document(str(letter)).paragraphs[1].runs
    assert [(r.text, r.bold) for r in runs if r.text] == [
        ('Sehr geehrte Damen, ', None),
        ('vielen Dank', True),
        (' für Ihre Nachricht.', None),
    ]


def test_the_changed_range_is_the_smallest_one():
    assert dv.changed_range('hello world', 'hello there') == (6, 11, 11)
    assert dv.changed_range('abc', 'abc') == (3, 3, 3)
    assert dv.changed_range('', 'new') == (0, 0, 3)
    assert dv.changed_range('gone', '') == (0, 4, 0)


def test_editing_a_table_cell_writes_to_that_cell(letter):
    result = dv.apply_edits(letter, {'t:0:1:1:0': '4,20'})
    dv.save(result['document'], letter)
    table = docx.Document(str(letter)).tables[0]
    assert table.cell(1, 1).text == '4,20'
    assert table.cell(1, 0).text == 'Schraube', 'a neighbouring cell changed'


def test_an_empty_paragraph_can_be_written_into(tmp_path):
    document = docx.Document()
    document.add_paragraph('')
    path = tmp_path / 'leer.docx'
    document.save(path)
    result = dv.apply_edits(path, {'p:0': 'jetzt mit Text'})
    dv.save(result['document'], path)
    assert docx.Document(str(path)).paragraphs[0].text == 'jetzt mit Text'


def test_an_address_that_does_not_resolve_is_an_error_not_a_silent_skip(letter):
    before = letter.read_bytes()
    with pytest.raises(dv.DocxError) as excinfo:
        dv.apply_edits(letter, {'p:999': 'nope'})
    assert 'left unchanged' in str(excinfo.value)
    assert letter.read_bytes() == before


def test_saving_nothing_is_refused(letter):
    with pytest.raises(dv.DocxError):
        dv.apply_edits(letter, {})


def test_a_round_trip_keeps_what_it_was_not_asked_to_change(letter):
    """python-docx writes back the package it read. The heading, the list
    style and the table survive an edit to an unrelated paragraph."""
    result = dv.apply_edits(letter, {'p:1': 'Kurz.'})
    dv.save(result['document'], letter)

    read = dv.read_document(letter)
    by_address = {b.address: b for b in read.blocks}
    assert by_address['p:0'].level == 1
    assert by_address['p:2'].listed is True
    assert read.tables == 1
    assert by_address['p:1'].text == 'Kurz.'


# ---------------------------------------------------------------------------
# The browser side
# ---------------------------------------------------------------------------

def test_a_block_is_reported_when_it_is_left_not_per_keystroke():
    script = dv.edit_js('calc-scope-1')
    assert "addEventListener('focusout'" in script
    for per_character in ("addEventListener('input'", "addEventListener('keyup'",
                          "addEventListener('keypress'"):
        assert per_character not in script, (
            f'{per_character} would put the kernel behind the typing')


def test_enter_does_not_create_a_paragraph_the_view_cannot_address():
    script = dv.edit_js('calc-scope-1')
    assert "e.key !== 'Enter'" in script
    assert 'preventDefault' in script


def test_marking_saved_clears_the_marks_without_rebuilding():
    script = dv.mark_saved_js('calc-scope-1')
    assert 'dw-dirty' in script
    assert 'innerHTML' not in script
