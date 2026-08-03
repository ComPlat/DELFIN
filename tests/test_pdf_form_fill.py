"""Filling a PDF form in the browser.

The pieces existed and only the agent could reach them: the field list and
the hardened writer in delfin.agent.office. The viewer now offers the same
thing to whoever is looking at the document.

Reading uses the library that draws the pages, because a field's rectangle
is what lets a click on it scroll to the place it fills in. Writing stays
with the agent's filler -- it is already hardened against the ways this
format fails quietly, and a second implementation would be a second set of
those failures.
"""

from __future__ import annotations

import importlib.util
import json

import pytest

if importlib.util.find_spec('fitz') is None:
    pytest.skip('PyMuPDF is not installed', allow_module_level=True)
pytest.importorskip('pypdf')
_canvas = pytest.importorskip('reportlab.pdfgen.canvas')
from reportlab.lib.pagesizes import A4

from delfin.dashboard import pdf_view as pv


def make_form_pdf(path, *, pages=1):
    """A PDF with a text field, a check box and a choice on each page."""
    c = _canvas.Canvas(str(path), pagesize=A4)
    for index in range(pages):
        c.setFont('Helvetica', 12)
        c.drawString(60, 760, f'Application, page {index + 1}')
        form = c.acroForm
        form.textfield(name=f'name{index}', x=60, y=700, width=200, height=20,
                       value='', borderStyle='inset')
        form.checkbox(name=f'agree{index}', x=60, y=660, size=16)
        form.choice(name=f'kind{index}', x=60, y=600, width=120, height=20,
                    options=[('One', 'One'), ('Two', 'Two')], value='One')
        c.showPage()
    c.save()
    return path


@pytest.fixture
def form(tmp_path):
    return make_form_pdf(tmp_path / 'antrag.pdf')


@pytest.fixture
def panel():
    made = pv.PdfPanel(run_js=lambda _s: None)
    try:
        yield made
    finally:
        made.close()


# ---------------------------------------------------------------------------
# Reading the fields
# ---------------------------------------------------------------------------

def test_fields_are_found_with_their_type_and_place(form):
    doc = pv.open_document(form)
    try:
        fields = pv.form_fields(doc)
    finally:
        doc.close()

    by_name = {f.name: f for f in fields}
    assert by_name['name0'].kind == 'text'
    assert by_name['agree0'].kind == 'check'
    assert by_name['kind0'].kind == 'choice'
    # The rectangle is the point: without it a click cannot go anywhere.
    assert by_name['name0'].rect[2] > by_name['name0'].rect[0]
    assert all(f.page == 0 for f in fields)


def test_a_document_without_fields_reports_none(tmp_path):
    plain = tmp_path / 'plain.pdf'
    c = _canvas.Canvas(str(plain), pagesize=A4)
    c.drawString(60, 700, 'no form here')
    c.save()
    doc = pv.open_document(plain)
    try:
        assert pv.form_fields(doc) == []
    finally:
        doc.close()


def test_fields_come_back_in_page_order(tmp_path):
    doc = pv.open_document(make_form_pdf(tmp_path / 'multi.pdf', pages=3))
    try:
        pages = [f.page for f in pv.form_fields(doc)]
    finally:
        doc.close()
    assert pages == sorted(pages)
    assert set(pages) == {0, 1, 2}


# ---------------------------------------------------------------------------
# The fields sit on the page
# ---------------------------------------------------------------------------

def _tell(panel, **message):
    message.setdefault('token', panel._token)
    panel._bridge_input.value = json.dumps(message)
    panel._bridge_btn.click()


def test_the_fields_are_drawn_over_the_page(panel, form):
    panel.open(form)
    markup = panel._overlays[0].value
    assert 'data-field="name0"' in markup
    assert 'type="checkbox"' in markup
    assert '<select' in markup
    assert 'position' not in markup, 'placement belongs in the stylesheet'
    assert 'left:' in markup and 'top:' in markup


def test_a_document_without_fields_draws_no_overlay(panel, tmp_path):
    plain = tmp_path / 'plain.pdf'
    c = _canvas.Canvas(str(plain), pagesize=A4)
    c.drawString(60, 700, 'nothing to fill')
    c.save()
    panel.open(plain)
    assert panel._overlays[0].value == ''
    assert panel.form_save_btn.layout.display == 'none', (
        'a Save button over a document with nothing to fill in can only '
        'disappoint')


def test_the_save_controls_appear_only_with_a_form(panel, form, tmp_path):
    panel.open(form)
    assert panel.form_save_btn.layout.display == ''
    assert panel.form_reset_btn.layout.display == ''
    assert 'field' in panel.form_status.value


def test_a_field_is_placed_where_the_page_puts_it(panel, form):
    """The same mapping the search highlights use, so a field and a hit on
    the same page agree about where the page is."""
    panel.open(form)
    doc = pv.open_document(form)
    try:
        fields = pv.form_fields(doc)
    finally:
        doc.close()
    zoom = pv.effective_dpi(*panel._sizes[0], panel.dpi()) / 72.0
    placed = pv.field_boxes(fields, 0, zoom)
    assert placed
    for found, (left, top, right, bottom) in placed:
        assert right > left and bottom > top
        assert f'left:{left}px' in panel._overlays[0].value


def test_the_fields_move_when_the_page_is_re_fitted(panel, form):
    panel.open(form)
    narrow = panel._overlays[0].value
    _tell(panel, action='width', px=1500)
    assert panel._overlays[0].value != narrow, (
        'the page was re-drawn at another size and the fields stayed put')


def test_a_read_only_field_is_shown_but_not_editable(panel, form):
    panel.open(form)
    panel._fields[0].readonly = True
    panel._paint_overlay(0)
    markup = panel._overlays[0].value
    assert 'disabled' in markup
    _tell(panel, action='fields', values={panel._fields[0].name: 'x'})
    assert panel._fields[0].name not in panel.form_values()


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def test_filling_and_saving_writes_the_values_and_keeps_a_copy(panel, form):
    from delfin.agent import office

    panel.open(form)
    _tell(panel, action='fields', values={'name0': 'Maxime Muster'})
    _tell(panel, action='fields', values={'agree0': True})
    _tell(panel, action='fields', values={'kind0': 'Two'})
    panel.save_form()

    assert 'saved' in panel.form_status.value
    written = {f['name']: f['value']
               for f in office.pdf_form_fields(form)['fields']}
    assert written['name0'] == 'Maxime Muster'
    assert written['kind0'] == 'Two'
    assert written['agree0'] not in ('', 'Off')
    assert list(form.parent.glob('antrag*.bak*')), 'no copy of the original'


def test_only_the_fields_that_were_touched_are_written(panel, form):
    """Writing back every field would rewrite ones nobody edited, and for a
    tick box that means deciding on the user's behalf what untouched meant."""
    panel.open(form)
    _tell(panel, action='fields', values={'name0': 'A'})
    assert set(panel.form_values()) == {'name0'}


def test_saving_nothing_says_so_rather_than_writing(panel, form):
    before = form.read_bytes()
    panel.open(form)
    panel.save_form()
    assert 'nothing filled in' in panel.form_status.value
    assert form.read_bytes() == before


def test_the_document_stays_usable_after_saving(panel, form):
    panel.open(form)
    _tell(panel, action='fields', values={'name0': 'A'})
    panel.save_form()
    assert panel.total_pages == 1
    assert panel.page_image(0) is not None
    assert 'data-field="name0"' in panel._overlays[0].value


def test_no_half_written_file_is_left_behind_on_failure(panel, form, monkeypatch):
    from delfin.agent import office

    panel.open(form)
    _tell(panel, action='fields', values={'name0': 'X'})
    monkeypatch.setattr(
        office, 'fill_pdf_form',
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError('nope')))
    panel.save_form()

    assert 'nope' in panel.form_status.value
    assert not list(form.parent.glob('*.filling.pdf'))
    assert form.exists()


def test_reset_clears_what_was_typed(panel, form):
    panel.open(form)
    _tell(panel, action='fields', values={'name0': 'typed something'})
    assert panel.form_values()
    panel.reset_form()
    assert panel.form_values() == {}
    assert 'typed something' not in panel._overlays[0].value


def test_the_one_page_viewer_offers_no_form(form):
    """The calculations browser keeps the viewer it had."""
    made = pv.PdfPanel(height_px=600, continuous=False)
    try:
        made.open(form)
        assert made._fields == []
        assert made.form_save_btn.layout.display == 'none'
    finally:
        made.close()


def make_paper_with_a_link_button(path):
    """A journal article: no form, but a pushbutton annotation on it."""
    import fitz

    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((60, 60), 'Physical Chemistry Chemical Physics')
    widget = fitz.Widget()
    widget.field_name = 'CrossMarkLinkButton'
    widget.field_type = fitz.PDF_WIDGET_TYPE_BUTTON
    widget.rect = fitz.Rect(60, 100, 120, 120)
    page.add_widget(widget)
    doc.save(str(path))
    doc.close()
    return path


def test_a_pushbutton_is_not_a_form_field(tmp_path):
    """Listing it turned a paper into a document with a one-field form
    that could not be saved."""
    path = make_paper_with_a_link_button(tmp_path / 'paper.pdf')
    doc = pv.open_document(path)
    try:
        assert pv.form_fields(doc) == []
    finally:
        doc.close()


def test_such_a_document_gets_no_form_panel(panel, tmp_path):
    panel.open(make_paper_with_a_link_button(tmp_path / 'paper.pdf'))
    assert panel._overlays[0].value == ''
    assert panel.form_save_btn.layout.display == 'none'


def test_several_fields_changed_at_once_all_arrive(panel, form):
    """Tabbing quickly through a form changes two fields in the same tick.
    Sent one at a time, each would write the single bridge field before the
    first had been read, and one of them would be lost."""
    panel.open(form)
    _tell(panel, action='fields',
          values={'name0': 'A', 'agree0': True, 'kind0': 'Two'})
    assert set(panel.form_values()) == {'name0', 'agree0', 'kind0'}


def test_the_browser_collects_them_before_sending():
    script = pv._pages_js('pdfv-1')
    assert 'pendingFields' in script
    assert "action: 'fields'" in script
    assert "action: 'field'," not in script, 'one message per field is the bug'


def test_a_junk_batch_is_ignored(panel, form):
    panel.open(form)
    _tell(panel, action='fields', values='not a dict')
    _tell(panel, action='fields', values={})
    assert panel.form_values() == {}


def test_a_field_covers_what_the_renderer_drew_underneath():
    """A filled field is drawn into the page image too, so a see-through box
    shows the value twice, half a pixel apart."""
    assert 'rgba' not in pv.FORM_CSS, 'a translucent field lets the page show through'
    assert 'background:#eaf1fd' in pv.FORM_CSS


# ---------------------------------------------------------------------------
# Fields that hold more than one line
# ---------------------------------------------------------------------------

def make_multiline_pdf(path):
    """A form with a one-line field and a box several lines tall."""
    import fitz

    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((60, 60), 'Antrag')

    single = fitz.Widget()
    single.field_name = 'name'
    single.field_type = fitz.PDF_WIDGET_TYPE_TEXT
    single.rect = fitz.Rect(60, 90, 300, 110)
    page.add_widget(single)

    wide = fitz.Widget()
    wide.field_name = 'begruendung'
    wide.field_type = fitz.PDF_WIDGET_TYPE_TEXT
    wide.field_flags = pv.TEXT_MULTILINE_FLAG
    wide.rect = fitz.Rect(60, 130, 400, 260)
    page.add_widget(wide)

    doc.save(str(path))
    doc.close()
    return path


def test_a_multiline_field_is_recognised(tmp_path):
    doc = pv.open_document(make_multiline_pdf(tmp_path / 'lang.pdf'))
    try:
        by_name = {f.name: f for f in pv.form_fields(doc)}
    finally:
        doc.close()
    assert by_name['begruendung'].multiline is True
    assert by_name['name'].multiline is False


def test_a_tall_box_counts_even_without_the_flag(tmp_path):
    """Plenty of forms are built by hand and never set the flag, and a
    two-centimetre box that scrolls one line sideways is unusable."""
    import fitz

    doc = fitz.open()
    page = doc.new_page()
    tall = fitz.Widget()
    tall.field_name = 'freitext'
    tall.field_type = fitz.PDF_WIDGET_TYPE_TEXT
    tall.rect = fitz.Rect(60, 100, 400, 200)     # 100 pt tall, no flag
    page.add_widget(tall)
    path = tmp_path / 'hoch.pdf'
    doc.save(str(path))
    doc.close()

    opened = pv.open_document(path)
    try:
        assert pv.form_fields(opened)[0].multiline is True
    finally:
        opened.close()


def test_a_multiline_field_is_drawn_as_a_wrapping_box(panel, tmp_path):
    panel.open(make_multiline_pdf(tmp_path / 'lang.pdf'))
    markup = panel._overlays[0].value
    assert '<textarea' in markup
    assert 'data-field="begruendung"' in markup
    # And the one-line field stays one line.
    assert 'data-field="name"' in markup
    assert markup.count('<textarea') == 1


def test_the_wrapping_box_wraps_rather_than_scrolling_sideways():
    assert 'white-space:pre-wrap' in pv.FORM_CSS
    assert 'textarea.pdfv-f' in pv.FORM_CSS


def test_a_multiline_field_is_not_given_letters_as_tall_as_its_box():
    """Sizing text to the height of the box is right for one line and
    absurd for a paragraph."""
    tall = pv.FormField(name='x', kind='text', multiline=True)
    one_line = pv.FormField(name='y', kind='text')
    assert pv.field_font_px(tall, 200, 1.0) < pv.field_font_px(one_line, 30, 1.0) * 2
    assert pv.field_font_px(tall, 200, 1.0) <= 14


def test_the_form_s_own_font_size_wins(tmp_path):
    """It is what the printed document will use."""
    stated = pv.FormField(name='x', kind='text', font_size=9.0)
    assert pv.field_font_px(stated, 40, 1.0) == 9
    assert pv.field_font_px(stated, 40, 2.0) == 18


def test_typing_into_a_wrapping_box_is_reported(panel, tmp_path):
    panel.open(make_multiline_pdf(tmp_path / 'lang.pdf'))
    _tell(panel, action='fields', values={'begruendung': 'Zeile eins\nZeile zwei'})
    assert panel.form_values()['begruendung'] == 'Zeile eins\nZeile zwei'


def test_the_browser_reports_a_text_area_the_same_way():
    script = pv._pages_js('pdfv-1')
    assert "el.tagName === 'TEXTAREA'" in script, (
        "a text area's type is not 'text', so it would never be reported")


def test_the_page_is_drawn_without_the_fields_when_the_overlay_has_them(panel, form):
    """A filled field is drawn into the page by the renderer, so an editable
    form showed its old value under the box being typed into -- and clearing
    a field left the previous text sitting there."""
    from delfin.agent import office

    # Give the form a value, so there is something that could show through.
    filled = form.with_name('gefuellt.pdf')
    office.fill_pdf_form(form, {'name0': 'Musterfrau'}, output=filled)

    panel.open(filled)
    with_fields = bytes(panel.page_image(0).value)

    doc = pv.open_document(filled)
    try:
        without = pv.render_page_png(doc, 0, panel.dpi(), annots=False)
        with_annots = pv.render_page_png(doc, 0, panel.dpi(), annots=True)
    finally:
        doc.close()

    assert with_annots != without, 'the value is drawn into the page'
    assert with_fields == without, 'the panel drew the value under its own box'


def test_a_document_without_a_form_still_shows_its_annotations(panel, tmp_path):
    """Leaving annotations out is for forms the overlay replaces, not for
    every stamp and comment in every document."""
    panel.open(make_paper_with_a_link_button(tmp_path / 'paper.pdf'))
    assert panel._fields == []
    doc = pv.open_document(tmp_path / 'paper.pdf')
    try:
        expected = pv.render_page_png(doc, 0, panel.dpi(), annots=True)
    finally:
        doc.close()
    assert bytes(panel.page_image(0).value) == expected


# ---------------------------------------------------------------------------
# Renaming a field
# ---------------------------------------------------------------------------

def make_badly_named_pdf(path):
    """What a form often arrives as: Text3, Text4, Check5."""
    import fitz

    doc = fitz.open()
    page = doc.new_page()
    for index, (name, top) in enumerate((('Text3', 100), ('Text4', 140))):
        widget = fitz.Widget()
        widget.field_name = name
        widget.field_type = fitz.PDF_WIDGET_TYPE_TEXT
        widget.rect = fitz.Rect(60, top, 300, top + 20)
        page.add_widget(widget)
    doc.save(str(path))
    doc.close()
    return path


def test_a_field_can_be_given_a_name_that_says_what_it_holds(tmp_path):
    path = make_badly_named_pdf(tmp_path / 'antrag.pdf')
    assert pv.rename_form_field(path, 'Text3', 'Nachname') == 'Nachname'

    doc = pv.open_document(path)
    try:
        names = [f.name for f in pv.form_fields(doc)]
    finally:
        doc.close()
    assert names == ['Nachname', 'Text4']


def test_the_new_name_is_what_the_writer_fills_by(tmp_path):
    """The name is the whole point: it is how the agent and anyone
    processing the document afterwards address the field."""
    from delfin.agent import office

    path = make_badly_named_pdf(tmp_path / 'antrag.pdf')
    pv.rename_form_field(path, 'Text3', 'Nachname')
    office.fill_pdf_form(path, {'Nachname': 'Musterfrau'},
                         output=tmp_path / 'out.pdf')
    written = {f['name']: f['value']
               for f in office.pdf_form_fields(tmp_path / 'out.pdf')['fields']}
    assert written['Nachname'] == 'Musterfrau'


def test_a_name_already_in_use_is_refused(tmp_path):
    path = make_badly_named_pdf(tmp_path / 'antrag.pdf')
    before = path.read_bytes()
    with pytest.raises(pv.PdfError) as excinfo:
        pv.rename_form_field(path, 'Text3', 'Text4')
    assert 'already' in str(excinfo.value)
    assert path.read_bytes() == before


def test_a_name_the_format_cannot_carry_is_refused(tmp_path):
    """A dot is the separator in the form's own hierarchy, so a field
    called a.b becomes a child of something that does not exist."""
    path = make_badly_named_pdf(tmp_path / 'antrag.pdf')
    for bad in ('a.b', 'x[1]', '   '):
        with pytest.raises(pv.PdfError):
            pv.rename_form_field(path, 'Text3', bad)


def test_renaming_a_field_that_is_not_there(tmp_path):
    path = make_badly_named_pdf(tmp_path / 'antrag.pdf')
    with pytest.raises(pv.PdfError) as excinfo:
        pv.rename_form_field(path, 'Nope', 'X')
    assert 'no field called' in str(excinfo.value)


def test_the_same_name_again_is_a_no_op(tmp_path):
    path = make_badly_named_pdf(tmp_path / 'antrag.pdf')
    before = path.read_bytes()
    assert pv.rename_form_field(path, 'Text3', 'Text3') == 'Text3'
    assert path.read_bytes() == before


def test_renaming_is_refused_while_entries_are_unsaved(panel, tmp_path):
    """The entries are held against the name they were typed into."""
    panel.open(make_badly_named_pdf(tmp_path / 'antrag.pdf'))
    _tell(panel, action='fields', values={'Text3': 'etwas'})
    _tell(panel, action='rename_field', old='Text3', name='Nachname')

    assert 'Save or reset' in panel.form_status.value
    doc = pv.open_document(tmp_path / 'antrag.pdf')
    try:
        assert [f.name for f in pv.form_fields(doc)] == ['Text3', 'Text4']
    finally:
        doc.close()


def test_renaming_through_the_panel_shows_the_new_name(panel, tmp_path):
    panel.open(make_badly_named_pdf(tmp_path / 'antrag.pdf'))
    _tell(panel, action='rename_field', old='Text3', name='Nachname')

    assert 'Nachname' in panel.form_status.value
    assert 'data-field="Nachname"' in panel._overlays[0].value


def test_the_rename_box_opens_on_the_field_itself():
    script = pv._pages_js('pdfv-1')
    assert "addEventListener('contextmenu'" in script
    body = script[script.index("addEventListener('contextmenu'"):][:900]
    assert 'pdfv-rename' in body
    assert "action: 'rename_field'" in body
