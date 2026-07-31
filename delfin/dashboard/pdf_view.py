"""PDF page viewer for the dashboard file browser.

Renders one page at a time to a PNG and hands the bytes to an
``ipywidgets.Image``. Everything happens in the kernel: no PDF plugin, no
script tag, no CDN -- the same server-side contract the rest of the browser
tab follows, which is what makes it work under Voila.

Cost model: a document is opened once and kept open while it is on screen,
but only the page in view is rasterised. Full-text search walks pages lazily
and stops at :data:`MAX_SEARCH_HITS` / :data:`MAX_SEARCH_PAGES`; whenever it
stops early the result object says so, so the UI can report a partial answer
instead of presenting it as complete.

Search hits are drawn onto the rendered page by blending colour into the
pixel buffer. The document itself is never modified -- no annotations are
added, nothing is written back to disk.

The geometry, capping and scan-detection logic is kept in module-level
functions that take plain data, so it can be unit tested without a browser
and (except for the few functions that take a page) without PyMuPDF.
"""

from __future__ import annotations

import html as _html
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple

import ipywidgets as widgets
import numpy as np

PDF_SUFFIXES = ('.pdf',)

# Rendering is streamed page by page, so file size only bounds how much MuPDF
# has to map -- not how much is rasterised.
MAX_FILE_BYTES = 400 * 1024 * 1024

# Zoom steps in DPI. The first one is the default: an A4 page is then about
# 830 px wide, which fits the content pane without horizontal scrolling and
# keeps the PNG small enough to travel over the comm channel on every page
# turn. The higher steps exist for scans, where the default is often too
# coarse to read.
DPI_STEPS = (100, 150, 200, 300)
DEFAULT_DPI_INDEX = 0
ZOOM_LABELS = ('100 %', '150 %', '200 %', '300 %')

# A page image crosses the comm channel as PNG on every navigation step, and an
# oversized sheet (posters, plans) at 300 dpi would be tens of megapixels. The
# effective DPI is lowered for those pages instead of refusing to draw them.
MAX_RENDER_PIXELS = 12_000_000
MIN_RENDER_DPI = 36

# Search bounds. Both are reported when they bite -- a silently shortened hit
# list would read like a complete answer.
MAX_SEARCH_HITS = 400
MAX_SEARCH_PAGES = 300

# Pages sampled when a document is opened to tell a scan from a text PDF.
TEXT_PROBE_PAGES = 8

HIGHLIGHT_COLOR = (255, 214, 0)      # every hit on the page
ACTIVE_HIGHLIGHT_COLOR = (255, 109, 0)  # the hit the user is standing on
HIGHLIGHT_ALPHA = 0.42

IDENTITY_MATRIX = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)


class PdfError(RuntimeError):
    """Raised for anything the user should see as a message, not a traceback."""


# ---------------------------------------------------------------------------
# Data carriers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SearchHit:
    """One match, addressed by page and by rectangle in PDF points."""

    page: int                                   # 0-based
    rect: Tuple[float, float, float, float]


@dataclass
class SearchResult:
    """Outcome of one search run, including how far it actually got."""

    term: str = ''
    hits: List[SearchHit] = field(default_factory=list)
    total_pages: int = 0
    pages_searched: int = 0
    pages_with_text: int = 0
    hit_cap_reached: bool = False
    page_cap_reached: bool = False

    @property
    def n_hits(self) -> int:
        return len(self.hits)

    @property
    def capped(self) -> bool:
        return self.hit_cap_reached or self.page_cap_reached

    @property
    def has_text(self) -> bool:
        """True if any searched page carried extractable text."""
        return self.pages_with_text > 0

    @property
    def pages_hit(self) -> List[int]:
        """0-based page numbers carrying at least one hit, in order."""
        seen: List[int] = []
        for hit in self.hits:
            if hit.page not in seen:
                seen.append(hit.page)
        return seen

    def hits_on(self, page: int) -> List[SearchHit]:
        return [hit for hit in self.hits if hit.page == page]


# ---------------------------------------------------------------------------
# Bounds and geometry (pure)
# ---------------------------------------------------------------------------

def clamp_page(page: Any, total_pages: int) -> int:
    """Clamp a 0-based page index into ``[0, total_pages - 1]``.

    Returns 0 for an empty document so callers never have to special-case it.
    """
    try:
        index = int(page)
    except (TypeError, ValueError):
        index = 0
    last = int(total_pages) - 1
    if last < 0:
        return 0
    return max(0, min(index, last))


def clamp_dpi_index(index: Any) -> int:
    """Clamp a zoom step index into the range of :data:`DPI_STEPS`."""
    try:
        value = int(index)
    except (TypeError, ValueError):
        value = DEFAULT_DPI_INDEX
    return max(0, min(value, len(DPI_STEPS) - 1))


def dpi_for_index(index: Any) -> int:
    return DPI_STEPS[clamp_dpi_index(index)]


def effective_dpi(width_pt: float, height_pt: float, dpi: int,
                  max_pixels: int = MAX_RENDER_PIXELS) -> int:
    """Lower ``dpi`` until the page fits into the pixel budget.

    Page dimensions are in PDF points (1/72 inch).
    """
    dpi = max(int(MIN_RENDER_DPI), int(dpi))
    width_pt = max(1.0, float(width_pt))
    height_pt = max(1.0, float(height_pt))
    pixels = (width_pt / 72.0 * dpi) * (height_pt / 72.0 * dpi)
    if pixels <= max_pixels:
        return dpi
    scale = math.sqrt(max_pixels / pixels)
    return max(MIN_RENDER_DPI, int(dpi * scale))


def rect_to_pixels(rect: Sequence[float], zoom: float, *,
                   matrix: Sequence[float] = IDENTITY_MATRIX,
                   origin: Sequence[float] = (0.0, 0.0),
                   bounds: Optional[Tuple[int, int]] = None
                   ) -> Optional[Tuple[int, int, int, int]]:
    """Map a PDF rectangle onto the pixel grid of a rendered page.

    ``matrix`` is a PDF transform ``(a, b, c, d, e, f)`` -- the page rotation
    matrix, because search results are reported in unrotated page space while
    the rendered image is rotated. All four corners are transformed and the
    axis-aligned hull is taken, which is what makes the 90/270 cases land on
    the right pixels. Bounds are rounded outwards so a highlight never falls
    short of the glyphs it covers.

    Returns ``None`` when the rectangle is empty or lies outside ``bounds``.
    """
    x0, y0, x1, y1 = (float(v) for v in rect)
    a, b, c, d, e, f = (float(v) for v in matrix)
    zoom = float(zoom)
    xs: List[float] = []
    ys: List[float] = []
    for px, py in ((x0, y0), (x1, y0), (x1, y1), (x0, y1)):
        xs.append((a * px + c * py + e) * zoom - float(origin[0]))
        ys.append((b * px + d * py + f) * zoom - float(origin[1]))
    left = int(math.floor(min(xs)))
    top = int(math.floor(min(ys)))
    right = int(math.ceil(max(xs)))
    bottom = int(math.ceil(max(ys)))
    if bounds is not None:
        width, height = int(bounds[0]), int(bounds[1])
        left = max(0, min(left, width))
        right = max(0, min(right, width))
        top = max(0, min(top, height))
        bottom = max(0, min(bottom, height))
    if right <= left or bottom <= top:
        return None
    return left, top, right, bottom


def blend_boxes(pixels: np.ndarray, boxes: Sequence[Sequence[int]],
                color: Sequence[int] = HIGHLIGHT_COLOR,
                alpha: float = HIGHLIGHT_ALPHA) -> np.ndarray:
    """Blend ``color`` into ``pixels`` inside each box, in place.

    Translucent rather than opaque: the glyphs under a hit have to stay
    readable, which is the whole point of showing where the hit is.
    """
    if boxes is None or len(boxes) == 0:
        return pixels
    height, width = pixels.shape[0], pixels.shape[1]
    channels = pixels.shape[2] if pixels.ndim == 3 else 1
    alpha = max(0.0, min(1.0, float(alpha)))
    tint = np.array(list(color)[:channels], dtype=np.float32)
    if tint.size < channels:
        tint = np.resize(tint, channels)
    for box in boxes:
        left, top, right, bottom = (int(v) for v in box)
        left = max(0, min(left, width))
        right = max(0, min(right, width))
        top = max(0, min(top, height))
        bottom = max(0, min(bottom, height))
        if right <= left or bottom <= top:
            continue
        region = pixels[top:bottom, left:right].astype(np.float32)
        pixels[top:bottom, left:right] = (
            region * (1.0 - alpha) + tint * alpha
        ).astype(np.uint8)
    return pixels


# ---------------------------------------------------------------------------
# Document access
# ---------------------------------------------------------------------------

def _fitz():
    """Import PyMuPDF on demand so a missing binding is a message, not a crash."""
    try:
        import fitz  # noqa: PLC0415 -- optional at import time by design
    except Exception as exc:                     # pragma: no cover - env specific
        raise PdfError(
            'PDF-Anzeige nicht verfügbar: PyMuPDF ist nicht installiert '
            f'({exc}).'
        ) from exc
    return fitz


def open_document(path: Path):
    """Open a PDF for viewing.

    Every failure mode a user can actually hit -- missing file, oversized file,
    password protection, damaged container -- comes back as :class:`PdfError`
    with a sentence that says what happened.
    """
    fitz = _fitz()
    path = Path(path)
    try:
        size = path.stat().st_size
    except OSError as exc:
        raise PdfError(f'Datei nicht lesbar: {exc}') from exc
    if size > MAX_FILE_BYTES:
        raise PdfError(
            f'PDF zu groß für die Anzeige ({size / (1024 * 1024):.1f} MB).'
        )
    try:
        doc = fitz.open(str(path))
    except Exception as exc:
        raise PdfError(f'PDF konnte nicht geöffnet werden: {exc}') from exc
    if doc.needs_pass:
        # An empty owner password unlocks plenty of "protected" office output.
        unlocked = False
        try:
            unlocked = bool(doc.authenticate(''))
        except Exception:
            unlocked = False
        if not unlocked:
            doc.close()
            raise PdfError(
                'Das PDF ist passwortgeschützt und kann hier nicht '
                'angezeigt werden.'
            )
    try:
        page_count = int(doc.page_count)
    except Exception as exc:
        doc.close()
        raise PdfError(f'PDF konnte nicht gelesen werden: {exc}') from exc
    if page_count <= 0:
        doc.close()
        raise PdfError('Das PDF enthält keine Seiten.')
    return doc


def page_text(doc, page_index: int) -> str:
    """Extractable text of one page ('' for a scanned page)."""
    try:
        return doc[int(page_index)].get_text('text') or ''
    except Exception:
        return ''


def probe_has_text(doc, limit: int = TEXT_PROBE_PAGES) -> Tuple[bool, int]:
    """Sample the first pages and report whether any of them carries text.

    Returns ``(has_text, pages_probed)``. Sampling the head is enough to tell a
    scan from a text PDF and costs nothing on a 500-page document.
    """
    total = int(getattr(doc, 'page_count', 0) or 0)
    probed = max(0, min(int(limit), total))
    for index in range(probed):
        if page_text(doc, index).strip():
            return True, probed
    return False, probed


def search_document(doc, term: str, *,
                    max_hits: int = MAX_SEARCH_HITS,
                    max_pages: int = MAX_SEARCH_PAGES) -> SearchResult:
    """Find ``term`` across the document, page by page.

    Stops at ``max_hits`` or after ``max_pages`` pages and records which limit
    bit, so the caller can tell "that is all there is" from "that is all we
    looked at". Pages carrying no extractable text are counted separately --
    zero hits on a scan is a different answer from zero hits in a text PDF.
    """
    result = SearchResult(term=str(term or ''))
    result.total_pages = int(getattr(doc, 'page_count', 0) or 0)
    if not result.term.strip() or result.total_pages <= 0:
        return result

    limit_pages = min(result.total_pages, max(1, int(max_pages)))
    result.page_cap_reached = limit_pages < result.total_pages

    for index in range(limit_pages):
        result.pages_searched = index + 1
        try:
            page = doc[index]
        except Exception:
            continue
        try:
            if (page.get_text('text') or '').strip():
                result.pages_with_text += 1
        except Exception:
            pass
        try:
            rects = page.search_for(result.term) or []
        except Exception:
            rects = []
        for rect in rects:
            result.hits.append(SearchHit(
                page=index,
                rect=(float(rect.x0), float(rect.y0),
                      float(rect.x1), float(rect.y1)),
            ))
            if len(result.hits) >= max_hits:
                result.hit_cap_reached = True
                return result
    return result


def render_page_png(doc, page_index: int, dpi: int = DPI_STEPS[DEFAULT_DPI_INDEX],
                    *, highlights: Sequence[Sequence[float]] = (),
                    active: Optional[Sequence[float]] = None) -> bytes:
    """Rasterise one page to PNG bytes, with optional hit highlights.

    ``highlights`` and ``active`` are rectangles in PDF points as returned by
    the search; the active one is tinted differently so it stands out among
    the other hits on the same page.
    """
    fitz = _fitz()
    total = int(getattr(doc, 'page_count', 0) or 0)
    index = clamp_page(page_index, total)
    try:
        page = doc[index]
    except Exception as exc:
        raise PdfError(f'Seite {index + 1} konnte nicht geladen werden: {exc}') from exc

    rect = page.rect
    dpi = effective_dpi(rect.width, rect.height, dpi)
    zoom = dpi / 72.0
    try:
        pix = page.get_pixmap(
            matrix=fitz.Matrix(zoom, zoom), colorspace=fitz.csRGB, alpha=False)
    except Exception as exc:
        raise PdfError(f'Seite {index + 1} konnte nicht gezeichnet werden: {exc}') from exc

    boxes = list(highlights or ())
    if active is not None:
        boxes.append(active)
    if not boxes or pix.n != 3:
        return pix.tobytes('png')

    rot = page.rotation_matrix
    matrix = (rot.a, rot.b, rot.c, rot.d, rot.e, rot.f)
    origin = (float(pix.x), float(pix.y))
    bounds = (pix.width, pix.height)
    # np.frombuffer is read-only; the copy is also what keeps the blend from
    # reaching back into MuPDF's own buffer.
    pixels = np.frombuffer(pix.samples, dtype=np.uint8).reshape(
        pix.height, pix.width, pix.n).copy()

    plain = []
    for hit_rect in (highlights or ()):
        box = rect_to_pixels(hit_rect, zoom, matrix=matrix, origin=origin, bounds=bounds)
        if box is not None:
            plain.append(box)
    blend_boxes(pixels, plain, HIGHLIGHT_COLOR, HIGHLIGHT_ALPHA)
    if active is not None:
        box = rect_to_pixels(active, zoom, matrix=matrix, origin=origin, bounds=bounds)
        if box is not None:
            blend_boxes(pixels, [box], ACTIVE_HIGHLIGHT_COLOR, HIGHLIGHT_ALPHA)

    out = fitz.Pixmap(fitz.csRGB, pix.width, pix.height, pixels.tobytes(), False)
    return out.tobytes('png')


# ---------------------------------------------------------------------------
# Wording (pure, so the honesty rules are testable)
# ---------------------------------------------------------------------------

def hit_labels(result: SearchResult) -> List[str]:
    """One entry per hit for the jump list, carrying its page number."""
    return [
        f'Treffer {i + 1} · Seite {hit.page + 1}'
        for i, hit in enumerate(result.hits)
    ]


def cap_note(result: SearchResult) -> str:
    """Plain sentence about what the search did *not* look at."""
    notes = []
    if result.hit_cap_reached:
        notes.append(f'bei {result.n_hits} Treffern abgebrochen')
    if result.page_cap_reached:
        notes.append(
            f'nur Seite 1–{result.pages_searched} von {result.total_pages} durchsucht'
        )
    return ', '.join(notes)


def search_summary_html(result: SearchResult, current: int = -1) -> str:
    """Status line for a finished search.

    Distinguishes three outcomes that must never be collapsed into one: hits,
    no hits in text that exists, and no extractable text at all.
    """
    if not result.term.strip():
        return ''
    note = cap_note(result)
    note_html = f' <span style="color:#777;">({_html.escape(note)})</span>' if note else ''
    if not result.hits:
        if not result.has_text:
            scope = (
                f'den ersten {result.pages_searched} Seiten'
                if result.page_cap_reached else 'diesem Dokument'
            )
            return (
                '<span style="color:#b26a00;">Kein durchsuchbarer Text in '
                f'{_html.escape(scope)} – vermutlich ein Scan '
                '(Bildseiten ohne Textebene).</span>'
            )
        return f'<span style="color:#d32f2f;">0 Treffer</span>{note_html}'
    pages = len(result.pages_hit)
    page_word = 'Seite' if pages == 1 else 'Seiten'
    if current is None or current < 0:
        return (
            f'<span style="color:#2e7d32;">{result.n_hits} Treffer auf '
            f'{pages} {page_word}</span>{note_html}'
        )
    hit = result.hits[max(0, min(int(current), result.n_hits - 1))]
    return (
        f'<b>{int(current) + 1}/{result.n_hits}</b> '
        f'<span style="color:#555;">(Seite {hit.page + 1} von {result.total_pages})</span>'
        f'{note_html}'
    )


def scan_hint_html(has_text: bool, pages_probed: int, total_pages: int) -> str:
    """Note shown when a document is opened, before anyone searches."""
    if has_text or total_pages <= 0:
        return ''
    scope = 'Alle Seiten' if pages_probed >= total_pages else f'Die ersten {pages_probed} Seiten'
    return (
        f'<span style="color:#b26a00;">{scope} enthalten keinen Text – '
        'vermutlich ein Scan. Die Volltextsuche findet hier nichts.</span>'
    )


def page_status_html(page_index: int, total_pages: int, dpi: int) -> str:
    return (
        f'<span style="color:#555;">Seite {page_index + 1} von {total_pages} '
        f'· {dpi} dpi</span>'
    )


# ---------------------------------------------------------------------------
# Widget panel
# ---------------------------------------------------------------------------

class PdfPanel:
    """Page viewer widget: navigation, zoom and in-document search.

    Holds the open document while it is on screen. :meth:`close` must be
    called when another file is selected -- an open MuPDF document keeps the
    file mapped, and the last rendered PNG is worth freeing too.
    """

    def __init__(self, *, height_px: int = 700):
        self._doc = None
        self._path: Optional[Path] = None
        self._page = 0
        self._dpi_index = DEFAULT_DPI_INDEX
        self._result = SearchResult()
        self._hit = -1
        self._syncing = False

        def _step_button(symbol, tooltip, disabled=False):
            return widgets.Button(
                description=symbol, tooltip=tooltip, disabled=disabled,
                layout=widgets.Layout(width='40px', height='28px'),
            )

        self.prev_page_btn = _step_button('◀', 'Vorherige Seite')
        self.next_page_btn = _step_button('▶', 'Nächste Seite')
        self.page_input = widgets.BoundedIntText(
            value=1, min=1, max=1,
            layout=widgets.Layout(width='70px', height='28px'),
        )
        self.page_total = widgets.HTML('von 1')
        self.zoom_dd = widgets.Dropdown(
            options=[(label, i) for i, label in enumerate(ZOOM_LABELS)],
            value=DEFAULT_DPI_INDEX,
            layout=widgets.Layout(width='92px'),
        )
        # Searching a long document is not free, so the term is taken on Enter
        # or on leaving the field -- never on every keystroke.
        self.search_input = widgets.Text(
            value='', placeholder='Im Dokument suchen', continuous_update=False,
            layout=widgets.Layout(width='210px', height='28px'),
        )
        self.search_input.add_class('delfin-nospell')
        self.search_btn = widgets.Button(
            description='Suchen',
            layout=widgets.Layout(width='84px', height='28px'),
        )
        self.hit_prev_btn = _step_button('◀', 'Vorheriger Treffer', disabled=True)
        self.hit_next_btn = _step_button('▶', 'Nächster Treffer', disabled=True)
        self.hit_select = widgets.Dropdown(
            options=[], value=None,
            layout=widgets.Layout(width='190px', display='none'),
        )
        self.status = widgets.HTML('')
        self.page_status = widgets.HTML('')
        # No max-width: the image has to keep its rendered pixel size, or the
        # zoom steps would be scaled away again by the browser.
        self.image = widgets.Image(
            format='png',
            layout=widgets.Layout(margin='0'),
        )

        self.toolbar = widgets.HBox(
            [
                self.prev_page_btn, self.page_input, self.page_total, self.next_page_btn,
                widgets.HTML('&nbsp;│&nbsp;'), self.zoom_dd,
                widgets.HTML('&nbsp;│&nbsp;'),
                self.search_input, self.search_btn,
                self.hit_prev_btn, self.hit_next_btn, self.hit_select,
                self.status,
            ],
            layout=widgets.Layout(
                width='100%', margin='4px 0', gap='6px',
                flex_flow='row wrap', align_items='center',
            ),
        )
        # A zoomed page is meant to overflow and scroll; centring it would put
        # the top-left corner out of reach as soon as it does.
        self.frame = widgets.Box(
            [self.image],
            layout=widgets.Layout(
                width='100%', height=f'{int(height_px)}px',
                overflow='auto', border='1px solid #ddd', padding='6px',
                justify_content='flex-start', align_items='flex-start',
            ),
        )
        self.widget = widgets.VBox(
            [self.toolbar, self.page_status, self.frame],
            layout=widgets.Layout(width='100%', display='none'),
        )

        self.prev_page_btn.on_click(lambda _b: self.step_page(-1))
        self.next_page_btn.on_click(lambda _b: self.step_page(1))
        self.page_input.observe(self._on_page_input, names='value')
        self.zoom_dd.observe(self._on_zoom, names='value')
        self.search_btn.on_click(lambda _b: self.run_search())
        self.search_input.observe(self._on_search_input, names='value')
        self.hit_prev_btn.on_click(lambda _b: self.step_hit(-1))
        self.hit_next_btn.on_click(lambda _b: self.step_hit(1))
        self.hit_select.observe(self._on_hit_select, names='value')

    # -- state ---------------------------------------------------------------

    @property
    def path(self) -> Optional[Path]:
        return self._path

    @property
    def page(self) -> int:
        return self._page

    @property
    def total_pages(self) -> int:
        return int(getattr(self._doc, 'page_count', 0) or 0) if self._doc else 0

    @property
    def result(self) -> SearchResult:
        return self._result

    def open(self, path) -> None:
        """Show the first page of ``path``; raises :class:`PdfError` on failure."""
        self.close()
        doc = open_document(Path(path))
        self._doc = doc
        self._path = Path(path)
        self._page = 0
        self._result = SearchResult()
        self._hit = -1
        self._syncing = True
        try:
            self.search_input.value = ''
        finally:
            self._syncing = False
        self._sync_hit_controls()
        has_text, probed = probe_has_text(doc)
        self.status.value = scan_hint_html(has_text, probed, self.total_pages)
        self._sync_page_controls()
        self.toolbar.layout.display = 'flex'
        self.widget.layout.display = 'flex'
        self._render()

    def close(self) -> None:
        """Drop the document and the rendered page."""
        if self._doc is not None:
            try:
                self._doc.close()
            except Exception:
                pass
        self._doc = None
        self._path = None
        self._result = SearchResult()
        self._hit = -1
        self.image.value = b''
        self.status.value = ''
        self.page_status.value = ''
        self.widget.layout.display = 'none'

    def show_error(self, message: str) -> None:
        """Put a failure on screen in place of a page.

        The message goes below the (now pointless) toolbar, which is hidden:
        there is no document to page through, zoom or search.
        """
        self.close()
        self.widget.layout.display = 'flex'
        self.toolbar.layout.display = 'none'
        self.page_status.value = (
            f'<span style="color:#d32f2f;">{_html.escape(str(message))}</span>'
        )

    # -- navigation ----------------------------------------------------------

    def goto_page(self, page_index: int, *, keep_hit: bool = False) -> None:
        if self._doc is None:
            return
        target = clamp_page(page_index, self.total_pages)
        changed = target != self._page
        self._page = target
        if changed and not keep_hit:
            # The hit cursor addresses a specific page; hand-navigating away
            # from it must not leave the counter pointing somewhere else.
            self._hit = -1
            self._sync_hit_controls()
        self._sync_page_controls()
        self._render()

    def step_page(self, delta: int) -> None:
        self.goto_page(self._page + int(delta))

    def set_zoom_index(self, index: int) -> None:
        self._dpi_index = clamp_dpi_index(index)
        self._sync_page_controls()
        self._render()

    # -- search --------------------------------------------------------------

    def run_search(self, term: Optional[str] = None) -> SearchResult:
        if term is not None:
            self._syncing = True
            try:
                self.search_input.value = str(term)
            finally:
                self._syncing = False
        query = (self.search_input.value or '').strip()
        if self._doc is None:
            return SearchResult(term=query)
        if not query:
            self._result = SearchResult()
            self._hit = -1
            self._sync_hit_controls()
            self.status.value = ''
            self._render()
            return self._result
        self._result = search_document(self._doc, query)
        self._hit = 0 if self._result.hits else -1
        self._sync_hit_controls()
        if self._result.hits:
            self.goto_page(self._result.hits[0].page, keep_hit=True)
        else:
            self._render()
        self.status.value = search_summary_html(self._result, self._hit)
        return self._result

    def step_hit(self, delta: int) -> None:
        total = self._result.n_hits
        if total <= 0:
            return
        delta = int(delta)
        if self._hit < 0:
            # No hit selected (the user paged by hand): enter the list at the
            # end the step is heading for.
            self._hit = 0 if delta >= 0 else total - 1
        else:
            self._hit = (self._hit + delta) % total
        self._goto_hit(self._hit)

    def _goto_hit(self, index: int) -> None:
        total = self._result.n_hits
        if total <= 0:
            return
        self._hit = max(0, min(int(index), total - 1))
        self._syncing = True
        try:
            self.hit_select.value = self._hit
        finally:
            self._syncing = False
        self.goto_page(self._result.hits[self._hit].page, keep_hit=True)
        self.status.value = search_summary_html(self._result, self._hit)

    # -- rendering -----------------------------------------------------------

    def _render(self) -> None:
        if self._doc is None:
            return
        page_hits = [hit.rect for hit in self._result.hits if hit.page == self._page]
        active = None
        if 0 <= self._hit < self._result.n_hits:
            current = self._result.hits[self._hit]
            if current.page == self._page:
                active = current.rect
                page_hits = [r for r in page_hits if r != active]
        try:
            self.image.value = render_page_png(
                self._doc, self._page, dpi_for_index(self._dpi_index),
                highlights=page_hits, active=active,
            )
        except Exception as exc:
            # A page that cannot be drawn is one damaged page, not a broken
            # viewer -- the rest of the document stays reachable.
            self.image.value = b''
            self.page_status.value = (
                f'<span style="color:#d32f2f;">{_html.escape(str(exc))}</span>'
            )

    # -- widget sync ---------------------------------------------------------

    def _sync_page_controls(self) -> None:
        total = max(1, self.total_pages)
        self._syncing = True
        try:
            self.page_input.max = total
            self.page_input.value = self._page + 1
        finally:
            self._syncing = False
        self.page_total.value = f'von {total}'
        self.prev_page_btn.disabled = self._page <= 0
        self.next_page_btn.disabled = self._page >= total - 1
        self.page_status.value = page_status_html(
            self._page, total, dpi_for_index(self._dpi_index))

    def _sync_hit_controls(self) -> None:
        hits = self._result.n_hits
        self.hit_prev_btn.disabled = hits < 2
        self.hit_next_btn.disabled = hits < 2
        self._syncing = True
        try:
            self.hit_select.options = [
                (label, i) for i, label in enumerate(hit_labels(self._result))
            ]
            if hits:
                self.hit_select.value = max(0, self._hit)
                self.hit_select.layout.display = ''
            else:
                self.hit_select.layout.display = 'none'
        finally:
            self._syncing = False

    def _on_page_input(self, change) -> None:
        if self._syncing or self._doc is None:
            return
        self.goto_page(int(change.get('new') or 1) - 1)

    def _on_search_input(self, change) -> None:
        if self._syncing or self._doc is None:
            return
        self.run_search()

    def _on_zoom(self, change) -> None:
        if self._syncing:
            return
        self.set_zoom_index(change.get('new'))

    def _on_hit_select(self, change) -> None:
        if self._syncing or self._doc is None:
            return
        value = change.get('new')
        if value is None:
            return
        self._goto_hit(int(value))
