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
# The panel
# ---------------------------------------------------------------------------

def test_the_form_panel_appears_only_when_there_is_a_form(panel, form, tmp_path):
    panel.open(form)
    assert panel.form_box.layout.display == 'flex'

    plain = tmp_path / 'plain.pdf'
    c = _canvas.Canvas(str(plain), pagesize=A4)
    c.drawString(60, 700, 'nothing to fill')
    c.save()
    panel.open(plain)
    assert panel.form_box.layout.display == 'none', (
        'an empty form panel reads like a broken feature rather than one '
        'that does not apply here')


def test_each_field_gets_the_control_its_type_needs(panel, form):
    import ipywidgets as widgets

    panel.open(form)
    assert isinstance(panel._controls['name0'], widgets.Text)
    assert isinstance(panel._controls['agree0'], widgets.Checkbox)
    assert isinstance(panel._controls['kind0'], widgets.Dropdown)


def test_the_panel_sits_beside_the_page_not_under_it(panel, form):
    panel.open(form)
    assert panel.frame in panel.body.children
    assert panel.form_box in panel.body.children


def test_clicking_a_field_goes_to_its_page_and_marks_it(panel, tmp_path):
    panel.open(make_form_pdf(tmp_path / 'multi.pdf', pages=3))
    index = next(i for i, f in enumerate(panel._fields) if f.page == 2)
    panel.goto_field(index)
    assert panel.page == 2
    assert panel._focus_rect[0] == 2


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def test_filling_and_saving_writes_the_values_and_keeps_a_copy(panel, form):
    from delfin.agent import office

    panel.open(form)
    panel._controls['name0'].value = 'Maxime Muster'
    panel._controls['agree0'].value = True
    panel._controls['kind0'].value = 'Two'
    panel.save_form()

    assert 'saved' in panel.form_status.value
    written = {f['name']: f['value']
               for f in office.pdf_form_fields(form)['fields']}
    assert written['name0'] == 'Maxime Muster'
    assert written['kind0'] == 'Two'
    assert written['agree0'] not in ('', 'Off')

    backups = list(form.parent.glob('antrag*.bak*'))
    assert backups, 'the original was replaced without keeping a copy'


def test_the_document_stays_usable_after_saving(panel, form):
    panel.open(form)
    panel._controls['name0'].value = 'A'
    panel.save_form()
    assert panel.total_pages == 1
    assert panel.page_image(0) is not None


def test_no_half_written_file_is_left_behind_on_failure(panel, form, monkeypatch):
    from delfin.agent import office

    panel.open(form)
    panel._controls['name0'].value = 'X'
    monkeypatch.setattr(
        office, 'fill_pdf_form',
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError('nope')))
    panel.save_form()

    assert 'nope' in panel.form_status.value
    assert not list(form.parent.glob('*.filling.pdf'))
    assert form.exists()


def test_reset_puts_the_controls_back(panel, form):
    panel.open(form)
    panel._controls['name0'].value = 'typed something'
    panel.reset_form()
    assert panel._controls['name0'].value == ''


def test_a_read_only_field_is_neither_editable_nor_written(panel, form):
    panel.open(form)
    panel._fields[0].readonly = True
    panel._controls[panel._fields[0].name].disabled = True
    assert panel._fields[0].name not in panel.form_values()


def test_the_one_page_viewer_offers_no_form(form):
    """The calculations browser keeps the viewer it had."""
    made = pv.PdfPanel(height_px=600, continuous=False)
    try:
        made.open(form)
        assert made.form_box.layout.display == 'none'
        assert made._fields == []
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
    assert panel.form_box.layout.display == 'none'
