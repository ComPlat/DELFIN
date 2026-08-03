"""The PDF view scrolls like a document, not like a slide deck.

Pages used to be shown one at a time in a frame of a fixed 680 px, while
every other view in the browser fills the pane to the bottom. So the PDF
was both shorter than the file list beside it and impossible to read
across a page break.

Now every page has a widget sized from its own rectangle -- which is what
makes the scroll bar right before anything is rasterised -- and only the
pages near the viewport hold pixels. The width follows the pane: the
browser reports how wide the frame is, and fit-to-width picks the DPI from
that, so dragging the splitter or going fullscreen re-fits the page.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

# Checked without importing it: loading PyMuPDF pulls in the libstdc++ that
# ships beside ORCA, and rdkit cannot be imported afterwards. Whichever test
# module touches it first decides that for the whole session, so this one
# stays out of it until a test actually opens a document.
if importlib.util.find_spec('fitz') is None:
    pytest.skip('PyMuPDF is not installed', allow_module_level=True)

_canvas = pytest.importorskip('reportlab.pdfgen.canvas')
from reportlab.lib.pagesizes import A4

from delfin.dashboard import pdf_view as pv


def make_text_pdf(path, pages, *, term='Rechnung'):
    """A text PDF where page i is identifiable and searchable."""
    c = _canvas.Canvas(str(path), pagesize=A4)
    for index in range(pages):
        c.setFont('Helvetica', 12)
        c.drawString(60, 500, f'Seite {index + 1} Inhalt')
        c.drawString(60, 460, f'{term} Nummer {index + 1}')
        c.showPage()
    c.save()
    return Path(path)


@pytest.fixture
def panel():
    scripts: list[str] = []
    made = pv.PdfPanel(run_js=scripts.append)
    made.scripts = scripts
    try:
        yield made
    finally:
        made.close()


# ---------------------------------------------------------------------------
# Geometry, without a document
# ---------------------------------------------------------------------------

def test_fit_width_picks_the_dpi_that_fills_the_frame():
    # A4 is 595 pt across; at 72 dpi that is 595 px.
    assert pv.fit_width_dpi(595.0, 595) == 72
    assert pv.fit_width_dpi(595.0, 1190) == 144


def test_fit_width_is_bounded_at_both_ends():
    """A sliver of a pane must not ask for an unreadable page, and a very
    wide one must not ask for a rasterisation nobody can see."""
    assert pv.fit_width_dpi(595.0, 10) == pv.MIN_RENDER_DPI
    assert pv.fit_width_dpi(595.0, 100_000) == pv.FIT_WIDTH_MAX_DPI


def test_fit_width_survives_a_nonsense_page():
    assert pv.fit_width_dpi(0, 800) == pv.FIT_WIDTH_MAX_DPI
    assert pv.fit_width_dpi(595.0, 0) == pv.MIN_RENDER_DPI


def test_page_pixel_size_respects_the_pixel_budget():
    """An oversized sheet is drawn smaller rather than refused."""
    normal = pv.page_pixel_size(595, 842, 100)
    assert normal == (826, 1169)
    huge = pv.page_pixel_size(5000, 5000, 300)
    assert huge[0] * huge[1] <= pv.MAX_RENDER_PIXELS


# ---------------------------------------------------------------------------
# The page stack
# ---------------------------------------------------------------------------

def test_every_page_gets_a_placeholder_at_its_own_size(panel, tmp_path):
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 6))
    assert len(panel.pages_box.children) == 6
    for box in panel.pages_box.children:
        assert box.layout.height.endswith('px')
        assert int(box.layout.height[:-2]) > 0


def test_the_sizes_are_known_without_rendering(tmp_path, monkeypatch):
    """The scroll bar has to be right from the first moment, and it is not
    worth forty rasterisations to get it there."""
    calls = []
    real = pv.render_page_png
    monkeypatch.setattr(
        pv, 'render_page_png',
        lambda doc, i, *a, **k: (calls.append(i), real(doc, i, *a, **k))[1])
    made = pv.PdfPanel()
    try:
        made.open(make_text_pdf(tmp_path / 'doc.pdf', 30))
        assert len(made._sizes) == 30
        assert max(calls) <= pv.RENDER_WINDOW
    finally:
        made.close()


def test_opening_asks_the_browser_to_report_what_is_visible(panel, tmp_path):
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 3))
    sent = '\n'.join(panel.scripts)
    assert panel._token in sent
    assert 'ResizeObserver' in sent
    assert "action: 'pages'" in sent
    assert "action: 'width'" in sent


def test_going_to_a_page_scrolls_instead_of_swapping(panel, tmp_path):
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 8))
    panel.scripts.clear()
    panel.goto_page(5)
    sent = '\n'.join(panel.scripts)
    assert '__pdfvGoto(5)' in sent


# ---------------------------------------------------------------------------
# What the browser reports
# ---------------------------------------------------------------------------

def _tell(panel, **message):
    message.setdefault('token', panel._token)
    panel._bridge_input.value = json.dumps(message)
    panel._bridge_btn.click()


def test_reported_pages_are_drawn_and_the_others_released(panel, tmp_path):
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 20))
    _tell(panel, action='pages', want=[12, 13])

    assert panel.page_image(12).value != b''
    assert panel.page_image(13).value != b''
    assert panel.page_image(0).value == b'', 'the first screen was never let go'
    # The page counter follows what is on screen.
    assert panel.page == 12


def test_a_wider_frame_re_fits_the_pages(panel, tmp_path):
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 3))
    narrow = panel.pages_box.children[0].layout.width

    _tell(panel, action='width', px=1600)

    assert panel.dpi() > pv.MIN_RENDER_DPI
    assert panel.pages_box.children[0].layout.width != narrow


def test_a_nudge_of_a_few_pixels_is_not_a_resize(panel, tmp_path):
    """Re-rasterising the document every time a scroll bar appears would
    make dragging the splitter unusable."""
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 3))
    _tell(panel, action='width', px=900)
    width = panel.pages_box.children[0].layout.width
    _tell(panel, action='width', px=908)
    assert panel.pages_box.children[0].layout.width == width


def test_a_message_for_another_panel_is_ignored(panel, tmp_path):
    """Three file browsers can each hold a PDF panel on the same page."""
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 20))
    _tell(panel, token='pdfv-somebody-else', action='pages', want=[15])
    assert panel.page_image(15).value == b''


def test_a_report_left_over_from_a_longer_document_is_dropped(panel, tmp_path):
    """The observer stays bound to the panel across documents, so a report
    in flight can name a page the new document does not have."""
    panel.open(make_text_pdf(tmp_path / 'long.pdf', 30))
    panel.open(make_text_pdf(tmp_path / 'short.pdf', 3))
    _tell(panel, action='pages', want=[25, 26])
    assert panel.page == 0, 'the page counter followed a page that is not there'
    assert len(panel.pages_box.children) == 3


def test_junk_on_the_bridge_does_not_raise(panel, tmp_path):
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 2))
    panel._bridge_input.value = 'not json'
    panel._bridge_btn.click()
    _tell(panel, action='pages', want=['x'])
    _tell(panel, action='width', px='wide')
    _tell(panel, action='nonsense')


# ---------------------------------------------------------------------------
# Fitting into the pane
# ---------------------------------------------------------------------------

def test_the_frame_takes_the_space_it_is_given(panel):
    assert panel.frame.layout.height in (None, '')
    assert panel.frame.layout.flex == '1 1 0'
    assert panel.frame.layout.min_height == '0'


def test_a_caller_can_still_pin_a_height():
    made = pv.PdfPanel(height_px=500)
    assert made.frame.layout.height == '500px'


def test_fit_to_width_is_the_default_and_still_offers_the_steps(panel):
    assert panel.zoom_dd.value == pv.FIT_WIDTH
    labels = [label for label, _value in panel.zoom_dd.options]
    assert labels[0] == pv.FIT_WIDTH_LABEL
    for step in pv.ZOOM_LABELS:
        assert step in labels


def test_choosing_a_fixed_step_stops_following_the_frame(panel, tmp_path):
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 2))
    panel.set_zoom(len(pv.DPI_STEPS) - 1)
    assert panel.dpi() == pv.DPI_STEPS[-1]
    _tell(panel, action='width', px=1600)
    assert panel.dpi() == pv.DPI_STEPS[-1]
