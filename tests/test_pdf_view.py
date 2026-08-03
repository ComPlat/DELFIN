"""Tests for the PDF viewer used by the file browser tab.

Test documents are generated with reportlab so no binary fixture has to live
in the repository.
"""

import importlib.util
import io
import sys
from pathlib import Path

import numpy as np
import pytest

_MODULE_PATH = Path(__file__).resolve().parents[1] / 'delfin' / 'dashboard' / 'pdf_view.py'
_SPEC = importlib.util.spec_from_file_location('_delfin_pdf_view', _MODULE_PATH)
pv = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = pv
_SPEC.loader.exec_module(pv)

fitz = pytest.importorskip('fitz')
_canvas = pytest.importorskip('reportlab.pdfgen.canvas')
_pil_image = pytest.importorskip('PIL.Image')
from reportlab.lib.pagesizes import A4  # noqa: E402


def png_pixels(png_bytes):
    """Decode rendered PNG bytes into an (h, w, 3) uint8 array."""
    with _pil_image.open(io.BytesIO(bytes(png_bytes))) as img:
        return np.asarray(img.convert('RGB'))


# ---------------------------------------------------------------------------
# Document builders
# ---------------------------------------------------------------------------

def make_text_pdf(path, pages, *, term='Rechnung', per_page=1):
    """A text PDF where page i carries ``term`` ``per_page`` times."""
    c = _canvas.Canvas(str(path), pagesize=A4)
    for index in range(pages):
        c.setFont('Helvetica', 12)
        c.drawString(60, 500, f'Seite {index + 1} Inhalt')
        for k in range(per_page):
            c.drawString(60, 460 - 20 * k, f'{term} Nummer {index + 1}-{k + 1}')
        c.showPage()
    c.save()
    return Path(path)


def make_image_only_pdf(path, pages=2):
    """A PDF with drawn shapes and no text layer at all -- a stand-in for a scan."""
    c = _canvas.Canvas(str(path), pagesize=A4)
    for _ in range(pages):
        c.setFillColorRGB(0.2, 0.2, 0.2)
        c.rect(60, 400, 300, 120, fill=1, stroke=0)
        c.setFillColorRGB(0.6, 0.6, 0.6)
        c.circle(300, 600, 80, fill=1, stroke=0)
        c.showPage()
    c.save()
    return Path(path)


# ---------------------------------------------------------------------------
# Bounds (pure)
# ---------------------------------------------------------------------------

def test_clamp_page_keeps_index_in_range():
    assert pv.clamp_page(-5, 10) == 0
    assert pv.clamp_page(0, 10) == 0
    assert pv.clamp_page(9, 10) == 9
    assert pv.clamp_page(99, 10) == 9


def test_clamp_page_on_empty_document_is_zero():
    assert pv.clamp_page(3, 0) == 0
    assert pv.clamp_page(0, -2) == 0


def test_clamp_page_survives_junk():
    assert pv.clamp_page(None, 5) == 0
    assert pv.clamp_page('abc', 5) == 0
    assert pv.clamp_page('4', 5) == 4


def test_dpi_steps_are_ascending_and_index_clamps():
    assert list(pv.DPI_STEPS) == sorted(pv.DPI_STEPS)
    assert len(pv.DPI_STEPS) == len(pv.ZOOM_LABELS)
    assert len(pv.DPI_STEPS) >= 3
    assert pv.clamp_dpi_index(-1) == 0
    assert pv.clamp_dpi_index(99) == len(pv.DPI_STEPS) - 1
    assert pv.clamp_dpi_index(None) == pv.DEFAULT_DPI_INDEX
    assert pv.dpi_for_index(0) == pv.DPI_STEPS[0]
    assert pv.dpi_for_index(99) == pv.DPI_STEPS[-1]


def test_effective_dpi_keeps_normal_pages_untouched():
    # A4 at 300 dpi is roughly 8.7 megapixels and must not be downscaled.
    assert pv.effective_dpi(595, 842, 300) == 300
    assert pv.effective_dpi(595, 842, 100) == 100


def test_effective_dpi_downscales_oversized_pages():
    # A0 poster at 300 dpi would be well above the pixel budget.
    lowered = pv.effective_dpi(2384, 3370, 300)
    assert pv.MIN_RENDER_DPI <= lowered < 300
    pixels = (2384 / 72 * lowered) * (3370 / 72 * lowered)
    assert pixels <= pv.MAX_RENDER_PIXELS


# ---------------------------------------------------------------------------
# Geometry (pure)
# ---------------------------------------------------------------------------

def test_rect_to_pixels_scales_and_rounds_outwards():
    box = pv.rect_to_pixels((10.4, 20.6, 30.2, 40.1), 2.0)
    assert box == (20, 41, 61, 81)


def test_rect_to_pixels_subtracts_the_pixmap_origin():
    box = pv.rect_to_pixels((10, 20, 30, 40), 1.0, origin=(10, 20))
    assert box == (0, 0, 20, 20)


def test_rect_to_pixels_applies_a_rotation_matrix():
    # 90 degrees on a 600 pt wide page: (x, y) -> (H - y, x) in fitz notation.
    matrix = (0.0, 1.0, -1.0, 0.0, 600.0, 0.0)
    box = pv.rect_to_pixels((10, 20, 30, 40), 1.0, matrix=matrix)
    assert box == (560, 10, 580, 30)


def test_rect_to_pixels_clamps_to_image_bounds():
    box = pv.rect_to_pixels((-50, -50, 100, 100), 1.0, bounds=(80, 60))
    assert box == (0, 0, 80, 60)


def test_rect_to_pixels_returns_none_outside_the_image():
    assert pv.rect_to_pixels((500, 500, 520, 520), 1.0, bounds=(100, 100)) is None
    assert pv.rect_to_pixels((10, 10, 10, 40), 1.0) is None


# ---------------------------------------------------------------------------
# Highlight blending (pure)
# ---------------------------------------------------------------------------

def test_blend_boxes_only_touches_the_box():
    pixels = np.full((20, 20, 3), 255, dtype=np.uint8)
    pv.blend_boxes(pixels, [(5, 5, 10, 10)], (255, 0, 0), 0.5)
    assert (pixels[5:10, 5:10] != 255).any()
    assert (pixels[0:5, :] == 255).all()
    assert (pixels[10:, :] == 255).all()
    assert (pixels[:, 10:] == 255).all()


def test_blend_boxes_is_translucent_not_opaque():
    pixels = np.zeros((10, 10, 3), dtype=np.uint8)
    pv.blend_boxes(pixels, [(0, 0, 10, 10)], (255, 255, 255), 0.4)
    # Black under a white 40 % wash stays clearly distinguishable from white.
    assert 80 < int(pixels[0, 0, 0]) < 130


def test_blend_boxes_ignores_degenerate_and_offscreen_boxes():
    pixels = np.full((10, 10, 3), 255, dtype=np.uint8)
    pv.blend_boxes(pixels, [(5, 5, 5, 9), (50, 50, 60, 60), (-9, -9, -1, -1)])
    assert (pixels == 255).all()


def test_blend_boxes_with_no_boxes_is_a_noop():
    pixels = np.full((4, 4, 3), 7, dtype=np.uint8)
    assert (pv.blend_boxes(pixels, []) == 7).all()


# ---------------------------------------------------------------------------
# Opening documents
# ---------------------------------------------------------------------------

def test_open_document_reports_page_count(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 3)
    doc = pv.open_document(path)
    try:
        assert doc.page_count == 3
    finally:
        doc.close()


def test_open_document_rejects_a_damaged_file(tmp_path):
    path = tmp_path / 'broken.pdf'
    path.write_bytes(b'%PDF-1.7\nthis is not a pdf body')
    with pytest.raises(pv.PdfError):
        pv.open_document(path)


def test_open_document_rejects_a_non_pdf(tmp_path):
    path = tmp_path / 'note.pdf'
    path.write_text('just text pretending to be a pdf')
    with pytest.raises(pv.PdfError):
        pv.open_document(path)


def test_open_document_rejects_a_missing_file(tmp_path):
    with pytest.raises(pv.PdfError):
        pv.open_document(tmp_path / 'nope.pdf')


def test_open_document_rejects_a_password_protected_file(tmp_path):
    path = tmp_path / 'secret.pdf'
    c = _canvas.Canvas(str(path), pagesize=A4, encrypt='geheim')
    c.drawString(60, 500, 'vertraulich')
    c.save()
    with pytest.raises(pv.PdfError) as excinfo:
        pv.open_document(path)
    assert 'password protected' in str(excinfo.value).lower()


def test_open_document_rejects_an_oversized_file(tmp_path, monkeypatch):
    path = make_text_pdf(tmp_path / 'doc.pdf', 1)
    monkeypatch.setattr(pv, 'MAX_FILE_BYTES', 10)
    with pytest.raises(pv.PdfError):
        pv.open_document(path)


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def test_search_finds_hits_with_page_numbers(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 3, per_page=2)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, 'Rechnung')
    finally:
        doc.close()
    assert result.n_hits == 6
    assert result.pages_hit == [0, 1, 2]
    assert [h.page for h in result.hits_on(1)] == [1, 1]
    assert not result.capped
    assert result.has_text


def test_search_without_matches_reports_zero_not_a_scan(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 2)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, 'Zahlungserinnerung')
    finally:
        doc.close()
    assert result.n_hits == 0
    assert result.has_text is True
    assert '0 hits' in pv.search_summary_html(result)


def test_search_on_a_page_without_text_reports_a_scan(tmp_path):
    path = make_image_only_pdf(tmp_path / 'scan.pdf', pages=2)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, 'Rechnung')
        assert pv.probe_has_text(doc) == (False, 2)
    finally:
        doc.close()
    assert result.n_hits == 0
    assert result.has_text is False
    summary = pv.search_summary_html(result)
    assert 'scan' in summary
    assert '0 Treffer' not in summary


def test_search_is_capped_by_hit_count_and_says_so(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 4, per_page=3)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, 'Rechnung', max_hits=5)
    finally:
        doc.close()
    assert result.n_hits == 5
    assert result.hit_cap_reached
    assert result.capped
    assert 'stopped at' in pv.cap_note(result)
    assert '5' in pv.search_summary_html(result)


def test_search_is_capped_by_page_count_and_says_so(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 6)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, 'Rechnung', max_pages=2)
    finally:
        doc.close()
    assert result.total_pages == 6
    assert result.pages_searched == 2
    assert result.page_cap_reached
    assert result.pages_hit == [0, 1]
    note = pv.cap_note(result)
    assert '2' in note and '6' in note


def test_search_ignores_an_empty_term(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 2)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, '   ')
    finally:
        doc.close()
    assert result.n_hits == 0
    assert result.pages_searched == 0
    assert pv.search_summary_html(result) == ''


def test_hit_labels_carry_page_numbers(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 2)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, 'Rechnung')
    finally:
        doc.close()
    labels = pv.hit_labels(result)
    assert labels[0].endswith('page 1')
    assert labels[1].endswith('page 2')


def test_scan_hint_only_fires_without_text():
    assert pv.scan_hint_html(True, 8, 20) == ''
    hint = pv.scan_hint_html(False, 8, 20)
    assert 'first 8' in hint and 'scan' in hint
    assert 'All pages' in pv.scan_hint_html(False, 3, 3)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def test_render_page_returns_png_bytes(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 2)
    doc = pv.open_document(path)
    try:
        png = pv.render_page_png(doc, 0, 100)
    finally:
        doc.close()
    assert png.startswith(b'\x89PNG\r\n\x1a\n')


def test_render_page_zoom_steps_change_the_pixel_size(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 1)
    doc = pv.open_document(path)
    try:
        sizes = []
        for dpi in pv.DPI_STEPS:
            sizes.append(png_pixels(pv.render_page_png(doc, 0, dpi)).shape[1])
    finally:
        doc.close()
    assert sizes == sorted(sizes)
    assert sizes[0] < sizes[-1]


def test_render_page_clamps_an_out_of_range_page(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 2)
    doc = pv.open_document(path)
    try:
        last = pv.render_page_png(doc, 99, 100)
        real = pv.render_page_png(doc, 1, 100)
    finally:
        doc.close()
    assert last == real


def test_render_highlight_paints_the_hit_and_leaves_the_rest_alone(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 1)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, 'Rechnung')
        rect = result.hits[0].rect
        plain = png_pixels(pv.render_page_png(doc, 0, 100))
        marked = png_pixels(pv.render_page_png(doc, 0, 100, highlights=[rect]))
        zoom = 100 / 72.0
        box = pv.rect_to_pixels(rect, zoom, bounds=(plain.shape[1], plain.shape[0]))
    finally:
        doc.close()
    assert box is not None
    left, top, right, bottom = box
    assert not np.array_equal(plain[top:bottom, left:right],
                              marked[top:bottom, left:right])
    # Everything outside the hit rectangle stays byte-identical.
    outside_plain = plain.copy()
    outside_marked = marked.copy()
    outside_plain[top:bottom, left:right] = 0
    outside_marked[top:bottom, left:right] = 0
    assert np.array_equal(outside_plain, outside_marked)


def test_render_marks_the_active_hit_differently(tmp_path):
    path = make_text_pdf(tmp_path / 'doc.pdf', 1, per_page=2)
    doc = pv.open_document(path)
    try:
        result = pv.search_document(doc, 'Rechnung')
        first, second = result.hits[0].rect, result.hits[1].rect
        png = pv.render_page_png(doc, 0, 100, highlights=[second], active=first)
        pixels = png_pixels(png)
        zoom = 100 / 72.0
        bounds = (pixels.shape[1], pixels.shape[0])
        box_a = pv.rect_to_pixels(first, zoom, bounds=bounds)
        box_b = pv.rect_to_pixels(second, zoom, bounds=bounds)
    finally:
        doc.close()
    corner_a = pixels[box_a[1] + 1, box_a[0] + 1].tolist()
    corner_b = pixels[box_b[1] + 1, box_b[0] + 1].tolist()
    assert corner_a != corner_b


@pytest.mark.parametrize('rotation', [0, 90, 180, 270])
def test_render_highlight_lands_on_the_glyphs_of_a_rotated_page(tmp_path, rotation):
    # One word on an otherwise empty page: every dark pixel of the rendered
    # page belongs to the hit, so the highlight has to cover all of them.
    # Without the page rotation matrix it lands somewhere else entirely.
    path = tmp_path / 'rot.pdf'
    c = _canvas.Canvas(str(path), pagesize=A4)
    c.setFont('Helvetica', 14)
    c.drawString(72, 300, 'Rechnung')
    c.showPage()
    c.save()
    doc = pv.open_document(path)
    try:
        doc[0].set_rotation(rotation)
        assert doc[0].rotation == rotation
        result = pv.search_document(doc, 'Rechnung')
        assert result.hits, 'rotated page should still be searchable'
        plain = png_pixels(pv.render_page_png(doc, 0, 150))
        marked = png_pixels(
            pv.render_page_png(doc, 0, 150, highlights=[result.hits[0].rect]))
    finally:
        doc.close()
    changed = np.any(plain != marked, axis=2)
    assert changed.sum() > 0, 'highlight was drawn nowhere'
    dark = (plain < 128).any(axis=2)
    assert dark.sum() > 0
    assert int((dark & changed).sum()) == int(dark.sum())


# ---------------------------------------------------------------------------
# Panel behaviour (widgets, no browser)
# ---------------------------------------------------------------------------

def test_panel_opens_on_the_first_page(tmp_path):
    panel = pv.PdfPanel()
    try:
        panel.open(make_text_pdf(tmp_path / 'doc.pdf', 5))
        assert panel.page == 0
        assert panel.total_pages == 5
        assert panel.page_input.value == 1
        assert panel.page_input.max == 5
        assert panel.page_total.value == 'of 5'
        assert panel.prev_page_btn.disabled
        assert not panel.next_page_btn.disabled
        assert bytes(panel.page_image(0).value).startswith(b'\x89PNG')
    finally:
        panel.close()


def test_panel_page_navigation_clamps_at_both_ends(tmp_path):
    panel = pv.PdfPanel()
    try:
        panel.open(make_text_pdf(tmp_path / 'doc.pdf', 3))
        panel.step_page(1)
        assert panel.page == 1
        panel.step_page(-5)
        assert panel.page == 0
        assert panel.prev_page_btn.disabled
        panel.goto_page(99)
        assert panel.page == 2
        assert panel.next_page_btn.disabled
    finally:
        panel.close()


def test_panel_jumps_to_the_page_of_a_hit(tmp_path):
    panel = pv.PdfPanel()
    try:
        panel.open(make_text_pdf(tmp_path / 'doc.pdf', 4))
        panel.goto_page(3)
        result = panel.run_search('Rechnung')
        assert result.n_hits == 4
        assert panel.page == 0
        assert len(panel.hit_select.options) == 4
        panel.step_hit(1)
        assert panel.page == 1
        panel.step_hit(-1)
        assert panel.page == 0
        # Wrapping backwards from the first hit lands on the last one.
        panel.step_hit(-1)
        assert panel.page == 3
    finally:
        panel.close()


def test_panel_renders_only_around_where_the_reader_is(tmp_path, monkeypatch):
    """The pages scroll continuously, so every page has a widget -- but only
    the ones near the viewport hold pixels. A forty-page document must not
    cost forty rasterisations to open, and the pages left behind have to be
    let go of again, or the viewer is a memory leak with a scroll bar."""
    rendered = []
    real_render = pv.render_page_png

    def _spy(doc, page_index, dpi=pv.DPI_STEPS[pv.DEFAULT_DPI_INDEX], **kwargs):
        rendered.append(page_index)
        return real_render(doc, page_index, dpi, **kwargs)

    monkeypatch.setattr(pv, 'render_page_png', _spy)
    panel = pv.PdfPanel()
    window = 2 * pv.RENDER_WINDOW + 1
    try:
        panel.open(make_text_pdf(tmp_path / 'doc.pdf', 40))
        assert len(panel.pages_box.children) == 40, 'the scroll bar needs them all'
        assert rendered == [0, 1, 2], 'opening rasterised more than the first screen'

        rendered.clear()
        panel.run_search('Rechnung')
        # A search over 40 pages redraws the window, not the document.
        assert len(rendered) <= window

        rendered.clear()
        panel.goto_page(20)
        assert set(rendered) == {18, 19, 20, 21, 22}
        # Page 0 was left far behind and must not still be holding a PNG.
        assert panel.page_image(0).value == b''
    finally:
        panel.close()


def test_panel_zoom_changes_the_rendered_image(tmp_path):
    panel = pv.PdfPanel()
    try:
        panel.open(make_text_pdf(tmp_path / 'doc.pdf', 1))
        panel.set_zoom_index(0)
        small = panel.page_image(0).value
        panel.set_zoom_index(len(pv.DPI_STEPS) - 1)
        assert len(panel.page_image(0).value) != len(small)
        assert str(pv.DPI_STEPS[-1]) in panel.page_status.value
    finally:
        panel.close()


def test_panel_reports_a_scan_instead_of_zero_hits(tmp_path):
    panel = pv.PdfPanel()
    try:
        panel.open(make_image_only_pdf(tmp_path / 'scan.pdf', pages=2))
        assert 'scan' in panel.status.value
        panel.run_search('Rechnung')
        assert 'scan' in panel.status.value
        assert '0 Treffer' not in panel.status.value
    finally:
        panel.close()


def test_panel_close_releases_the_document(tmp_path):
    panel = pv.PdfPanel()
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 2))
    panel.close()
    assert panel.path is None
    assert panel.total_pages == 0
    assert panel.page_image(0) is None
    assert panel.pages_box.children == ()
    assert panel.widget.layout.display == 'none'
    # Navigation after closing must not raise.
    panel.step_page(1)
    panel.run_search('Rechnung')


def test_panel_show_error_replaces_the_page(tmp_path):
    panel = pv.PdfPanel()
    panel.open(make_text_pdf(tmp_path / 'doc.pdf', 1))
    panel.show_error('PDF konnte nicht geöffnet werden: kaputt')
    assert panel.widget.layout.display == 'flex'
    assert panel.toolbar.layout.display == 'none'
    assert 'kaputt' in panel.page_status.value
    assert panel.pages_box.children == ()


def test_panel_hit_stepping_after_manual_paging_enters_the_list_cleanly(tmp_path):
    panel = pv.PdfPanel()
    try:
        panel.open(make_text_pdf(tmp_path / 'doc.pdf', 4))
        panel.run_search('Rechnung')
        panel.goto_page(2)          # drops the hit cursor
        panel.step_hit(1)
        assert panel.page == 0
        panel.goto_page(2)
        panel.step_hit(-1)
        assert panel.page == 3
    finally:
        panel.close()


# ---------------------------------------------------------------------------
# Dispatch from the file browser tab
# ---------------------------------------------------------------------------

@pytest.fixture
def browser_tab(tmp_path):
    """A real browser tab rooted at an empty directory."""
    from delfin.dashboard.context import DashboardContext
    from delfin.dashboard import tab_calculations_browser as tb

    root = tmp_path / 'office'
    root.mkdir()
    ctx = DashboardContext(
        calc_dir=root, archive_dir=tmp_path / 'archive', office_dir=root,
    )
    _widget, refs = tb.create_tab(ctx)
    return root, refs


def _select(refs, filename):
    file_list = refs['calc_file_list']
    refs['calc_list_directory']()
    label = [o for o in file_list.options if o.endswith(filename)][0]
    file_list.value = ()
    file_list.value = (label,)


def test_browser_opens_a_pdf_in_the_page_viewer(browser_tab):
    root, refs = browser_tab
    make_text_pdf(root / 'Formular.pdf', 3)
    _select(refs, 'Formular.pdf')
    state = refs['xyz_batch_state']
    panel = state['pdf_panel']
    assert state['pdf_active'] is True
    assert panel.total_pages == 3
    assert bytes(panel.page_image(0).value).startswith(b'\x89PNG')


def test_browser_releases_the_pdf_when_another_file_is_selected(browser_tab):
    root, refs = browser_tab
    make_text_pdf(root / 'Formular.pdf', 2)
    (root / 'notiz.txt').write_text('nur Text\n')
    _select(refs, 'Formular.pdf')
    state = refs['xyz_batch_state']
    _select(refs, 'notiz.txt')
    assert state['pdf_active'] is False
    assert state['pdf_panel'].path is None


def test_browser_releases_the_pdf_when_the_folder_is_re_listed(browser_tab):
    root, refs = browser_tab
    make_text_pdf(root / 'Formular.pdf', 2)
    _select(refs, 'Formular.pdf')
    state = refs['xyz_batch_state']
    refs['calc_list_directory']()
    assert state['pdf_active'] is False
    assert state['pdf_panel'].path is None


def test_browser_shows_a_message_for_a_damaged_pdf(browser_tab):
    root, refs = browser_tab
    (root / 'kaputt.pdf').write_bytes(b'%PDF-1.7\nno body here')
    _select(refs, 'kaputt.pdf')
    panel = refs['xyz_batch_state']['pdf_panel']
    assert 'could not be opened' in panel.page_status.value
    assert panel.pages_box.children == ()
