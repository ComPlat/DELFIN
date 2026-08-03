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
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

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

# Fit-to-width is the default: a document is read at the width it has room
# for, and that width changes when the splitter moves, the window resizes or
# the tab goes fullscreen. Kept out of DPI_STEPS because it is not a step --
# it is a rule for computing one.
FIT_WIDTH = -1
FIT_WIDTH_LABEL = 'Fit width'
FIT_WIDTH_MAX_DPI = 300
# Until the browser reports the real width. An A4 page at 100 dpi is ~830 px,
# which is what the pane used to be sized for.
ASSUMED_FRAME_WIDTH = 860

# Pages are stacked in one scroll container. Each page keeps its own widget,
# and only the ones near the viewport hold pixels.
PAGE_GAP_PX = 12
RENDER_WINDOW = 2

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


_PAGES_JS_TEMPLATE = r"""
(function(){
  var TOKEN = __TOKEN__;
  var tries = 0;

  function boot(){
    var root = document.querySelector('.' + TOKEN);
    var frame = root && root.querySelector('.pdfv-frame');
    var pages = root && root.querySelector('.pdfv-pages');
    if (!root || !frame || !pages || !pages.children.length) {
      if (++tries < 60) setTimeout(boot, 50);
      return;
    }
    if (root.dataset.bound === TOKEN) { report(); measure(); return; }
    root.dataset.bound = TOKEN;
    install(root, frame, pages);
  }

  function install(root, frame, pages){
    function send(payload){
      payload.token = TOKEN;
      var field = root.querySelector('.pdfv-bridge-input input');
      if (!field) return;
      var desc = Object.getOwnPropertyDescriptor(
        HTMLInputElement.prototype, 'value');
      if (desc && desc.set) desc.set.call(field, JSON.stringify(payload));
      else field.value = JSON.stringify(payload);
      field.dispatchEvent(new Event('input', {bubbles: true}));
      field.dispatchEvent(new Event('change', {bubbles: true}));
      var holder = root.querySelector('.pdfv-bridge-action');
      var btn = holder && (holder.tagName === 'BUTTON'
        ? holder : holder.querySelector('button'));
      if (btn) btn.click();
    }

    /* Which pages are on screen. Reported as one message after scrolling
       settles: asking per page would send a burst the kernel answers one
       render at a time, and the user would watch it catch up. */
    var pending = null;
    function report(){
      if (pending) clearTimeout(pending);
      pending = setTimeout(function(){
        pending = null;
        var top = frame.scrollTop, bottom = top + frame.clientHeight;
        var want = [];
        var kids = pages.children;
        for (var i = 0; i < kids.length; i++) {
          var el = kids[i];
          var elTop = el.offsetTop, elBottom = elTop + el.offsetHeight;
          if (elBottom > top - 40 && elTop < bottom + 40) want.push(i);
        }
        if (want.length) send({action: 'pages', want: want});
      }, 90);
    }
    root.__pdfvReport = report;

    var lastWidth = 0;
    function measure(){
      var width = Math.round(frame.clientWidth);
      if (!width || Math.abs(width - lastWidth) < 8) return;
      lastWidth = width;
      send({action: 'width', px: width});
    }
    root.__pdfvMeasure = measure;

    frame.addEventListener('scroll', report, {passive: true});
    if (window.ResizeObserver) {
      var timer = null;
      new ResizeObserver(function(){
        if (timer) clearTimeout(timer);
        timer = setTimeout(measure, 160);
      }).observe(frame);
    } else {
      window.addEventListener('resize', function(){ setTimeout(measure, 160); });
    }

    /* Scroll to a page, asked for by the kernel when the page number, a
       search hit or the step buttons move. */
    root.__pdfvGoto = function(index){
      var el = pages.children[index];
      if (!el) return;
      frame.scrollTop = Math.max(0, el.offsetTop - 6);
      report();
    };

    measure();
    report();
  }

  function report(){
    var root = document.querySelector('.' + TOKEN);
    if (root && root.__pdfvReport) root.__pdfvReport();
  }
  function measure(){
    var root = document.querySelector('.' + TOKEN);
    if (root && root.__pdfvMeasure) root.__pdfvMeasure();
  }

  boot();
})();
"""


def _pages_js(token: str) -> str:
    """Browser side of the page stack: report what is visible, measure width.

    Must go through the dashboard's JS channel: ipywidgets sets HTML content
    with ``innerHTML``, which never runs an inline script.
    """
    return _PAGES_JS_TEMPLATE.replace('__TOKEN__', json.dumps(str(token)))


def _goto_js(token: str, index: int) -> str:
    return (
        f"(function(){{var r=document.querySelector('.{token}');"
        f"if(r&&r.__pdfvGoto)r.__pdfvGoto({int(index)});}})();"
    )


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


def fit_width_dpi(width_pt: float, frame_px: float) -> int:
    """DPI at which a page of ``width_pt`` fills ``frame_px`` across.

    Bounded at both ends: below :data:`MIN_RENDER_DPI` the page is unreadable
    anyway, and above :data:`FIT_WIDTH_MAX_DPI` a wide pane would ask for a
    rasterisation nobody can see the benefit of but everybody waits for.
    """
    width_pt = max(1.0, float(width_pt))
    frame_px = max(1.0, float(frame_px))
    dpi = int(frame_px * 72.0 / width_pt)
    return max(MIN_RENDER_DPI, min(FIT_WIDTH_MAX_DPI, dpi))


def page_sizes(doc) -> List[Tuple[float, float]]:
    """Width and height in points for every page, without rendering any.

    This is what lets the scroll bar be right from the first moment: the
    placeholders are the size the pages will be.
    """
    sizes: List[Tuple[float, float]] = []
    for index in range(int(getattr(doc, 'page_count', 0) or 0)):
        try:
            rect = doc.load_page(index).rect
            sizes.append((float(rect.width), float(rect.height)))
        except Exception:
            # One unreadable page must not cost the whole document its
            # scroll bar; A4 is a better guess than nothing.
            sizes.append((595.0, 842.0))
    return sizes


def page_pixel_size(width_pt: float, height_pt: float, dpi: int
                    ) -> Tuple[int, int]:
    """On-screen size of a page at ``dpi``, after the pixel budget applies."""
    dpi = effective_dpi(width_pt, height_pt, dpi)
    scale = dpi / 72.0
    return max(1, int(width_pt * scale)), max(1, int(height_pt * scale))


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
            'PDF viewing is unavailable: PyMuPDF is not installed '
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
            f'PDF is too large to display ({size / (1024 * 1024):.1f} MB).'
        )
    try:
        doc = fitz.open(str(path))
    except Exception as exc:
        raise PdfError(f'PDF could not be opened: {exc}') from exc
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
                'The PDF is password protected and cannot be shown here.'
            )
    try:
        page_count = int(doc.page_count)
    except Exception as exc:
        doc.close()
        raise PdfError(f'PDF konnte nicht gelesen werden: {exc}') from exc
    if page_count <= 0:
        doc.close()
        raise PdfError('The PDF has no pages.')
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
        raise PdfError(f'Page {index + 1} could not be loaded: {exc}') from exc

    rect = page.rect
    dpi = effective_dpi(rect.width, rect.height, dpi)
    zoom = dpi / 72.0
    try:
        pix = page.get_pixmap(
            matrix=fitz.Matrix(zoom, zoom), colorspace=fitz.csRGB, alpha=False)
    except Exception as exc:
        raise PdfError(f'Page {index + 1} could not be drawn: {exc}') from exc

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
        f'Hit {i + 1} · page {hit.page + 1}'
        for i, hit in enumerate(result.hits)
    ]


def cap_note(result: SearchResult) -> str:
    """Plain sentence about what the search did *not* look at."""
    notes = []
    if result.hit_cap_reached:
        notes.append(f'stopped at {result.n_hits} hits')
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
                f'the first {result.pages_searched} pages'
                if result.page_cap_reached else 'this document'
            )
            return (
                '<span style="color:#b26a00;">No searchable text in '
                f'{_html.escape(scope)} - probably a scan '
                '(page images without a text layer).</span>'
            )
        return f'<span style="color:#d32f2f;">0 hits</span>{note_html}'
    pages = len(result.pages_hit)
    page_word = 'page' if pages == 1 else 'pages'
    if current is None or current < 0:
        return (
            f'<span style="color:#2e7d32;">{result.n_hits} hits on '
            f'{pages} {page_word}</span>{note_html}'
        )
    hit = result.hits[max(0, min(int(current), result.n_hits - 1))]
    return (
        f'<b>{int(current) + 1}/{result.n_hits}</b> '
        f'<span style="color:#555;">(page {hit.page + 1} of {result.total_pages})</span>'
        f'{note_html}'
    )


def scan_hint_html(has_text: bool, pages_probed: int, total_pages: int) -> str:
    """Note shown when a document is opened, before anyone searches."""
    if has_text or total_pages <= 0:
        return ''
    scope = 'All pages' if pages_probed >= total_pages else f'The first {pages_probed} pages'
    return (
        f'<span style="color:#b26a00;">{scope} hold no text - '
        'probably a scan. Full-text search will not find anything here.</span>'
    )


def page_status_html(page_index: int, total_pages: int, dpi: int) -> str:
    return (
        f'<span style="color:#555;">Page {page_index + 1} of {total_pages} '
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

    def __init__(self, *, height_px: Optional[int] = None, run_js=None,
                 continuous: bool = True):
        self._doc = None
        self._path: Optional[Path] = None
        self._page = 0
        self._continuous = bool(continuous)
        self._zoom = FIT_WIDTH if self._continuous else DEFAULT_DPI_INDEX
        self._result = SearchResult()
        self._hit = -1
        self._syncing = False
        self._run_js = run_js
        self._token = f'pdfv-{id(self):x}'
        self._sizes: List[Tuple[float, float]] = []
        self._page_images: Dict[int, widgets.Image] = {}
        self._page_boxes: Dict[int, widgets.Box] = {}
        self._filled: set = set()
        self._frame_px = ASSUMED_FRAME_WIDTH

        def _step_button(symbol, tooltip, disabled=False):
            return widgets.Button(
                description=symbol, tooltip=tooltip, disabled=disabled,
                layout=widgets.Layout(width='40px', height='28px'),
            )

        self.prev_page_btn = _step_button('◀', 'Previous page')
        self.next_page_btn = _step_button('▶', 'Next page')
        self.page_input = widgets.BoundedIntText(
            value=1, min=1, max=1,
            layout=widgets.Layout(width='70px', height='28px'),
        )
        self.page_total = widgets.HTML('of 1')
        # Fit-to-width only where the frame is measured; the one-page view
        # keeps the fixed steps it always had.
        _zoom_options = [(label, i) for i, label in enumerate(ZOOM_LABELS)]
        if self._continuous:
            _zoom_options.insert(0, (FIT_WIDTH_LABEL, FIT_WIDTH))
        self.zoom_dd = widgets.Dropdown(
            options=_zoom_options,
            value=self._zoom,
            layout=widgets.Layout(width='124px' if self._continuous else '92px'),
        )
        # Searching a long document is not free, so the term is taken on Enter
        # or on leaving the field -- never on every keystroke.
        self.search_input = widgets.Text(
            value='', placeholder='Search in document', continuous_update=False,
            layout=widgets.Layout(width='210px', height='28px'),
        )
        self.search_input.add_class('delfin-nospell')
        self.search_btn = widgets.Button(
            description='Search',
            layout=widgets.Layout(width='84px', height='28px'),
        )
        self.hit_prev_btn = _step_button('◀', 'Previous hit', disabled=True)
        self.hit_next_btn = _step_button('▶', 'Next hit', disabled=True)
        self.hit_select = widgets.Dropdown(
            options=[], value=None,
            layout=widgets.Layout(width='190px', display='none'),
        )
        self.status = widgets.HTML('')
        self.page_status = widgets.HTML('')

        # The pages, stacked. Built when a document is opened.
        self.pages_box = widgets.VBox(
            [],
            layout=widgets.Layout(width='100%', margin='0'),
        )
        self.pages_box.add_class('pdfv-pages')

        # Browser -> kernel. The browser reports which pages came into view
        # and how wide the frame is; the kernel answers by setting widget
        # values, which is a normal traitlet sync rather than another script.
        self._bridge_input = widgets.Text(
            value='', layout=widgets.Layout(width='1px', height='1px'))
        self._bridge_input.add_class('pdfv-bridge-input')
        self._bridge_btn = widgets.Button(
            description='', layout=widgets.Layout(width='1px', height='1px'))
        self._bridge_btn.add_class('pdfv-bridge-action')
        self._bridge_btn.on_click(self._on_bridge)
        self._bridge = widgets.Box(
            [self._bridge_input, self._bridge_btn],
            layout=widgets.Layout(display='none'),
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
        #
        # Height comes from flex, not from a pixel count: every other view in
        # the browser fills the pane down to the bottom, and a PDF frame with
        # a fixed height stopped short of it and stayed short in fullscreen.
        frame_layout = dict(
            width='100%', overflow='auto', border='1px solid #ddd',
            padding='6px', justify_content='flex-start',
            align_items='flex-start',
        )
        if height_px:
            frame_layout['height'] = f'{int(height_px)}px'
        else:
            frame_layout.update(flex='1 1 0', min_height='0')
        self.frame = widgets.Box(
            [self.pages_box], layout=widgets.Layout(**frame_layout))
        self.frame.add_class('pdfv-frame')

        self.widget = widgets.VBox(
            [self.toolbar, self.page_status, self.frame, self._bridge],
            layout=widgets.Layout(
                width='100%', display='none', flex='1 1 0', min_height='0'),
        )
        self.widget.add_class('pdfv-root')
        self.widget.add_class(self._token)

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
        """Show ``path``; raises :class:`PdfError` on failure."""
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
        self._sizes = page_sizes(doc)
        self._build_pages()
        self._sync_page_controls()
        self.toolbar.layout.display = 'flex'
        self.widget.layout.display = 'flex'
        self._render()
        self._install_observer()

    # -- the page stack ------------------------------------------------------

    def dpi(self) -> int:
        """The DPI pages are drawn at right now."""
        if self._zoom == FIT_WIDTH and self._continuous:
            widest = max((w for w, _h in self._sizes), default=595.0)
            # The frame's own padding and the scroll bar are not page.
            return fit_width_dpi(widest, max(120, self._frame_px - 28))
        return dpi_for_index(self._zoom)

    def _shown_pages(self) -> List[int]:
        """Which pages have a widget.

        All of them when the view scrolls; only the current one otherwise,
        which is the one-page-at-a-time view the calculations browser has
        always had and keeps.
        """
        if not self._continuous:
            return [self._page] if self._sizes else []
        return list(range(len(self._sizes)))

    def _build_pages(self) -> None:
        """A placeholder per shown page, at the size that page will be.

        In the scrolling view that makes the scroll bar right before
        anything has been rasterised, rather than growing under the reader.
        """
        dpi = self.dpi()
        self._page_images = {}
        self._page_boxes = {}
        self._filled = set()
        order = []
        for index in self._shown_pages():
            width_pt, height_pt = self._sizes[index]
            width_px, height_px = page_pixel_size(width_pt, height_pt, dpi)
            image = widgets.Image(
                format='png',
                layout=widgets.Layout(
                    width=f'{width_px}px', height=f'{height_px}px', margin='0'),
            )
            box = widgets.Box(
                [image],
                layout=widgets.Layout(
                    width=f'{width_px}px', height=f'{height_px}px',
                    margin=f'0 0 {PAGE_GAP_PX}px 0',
                    background_color='#f4f4f4',
                    border='1px solid #e0e0e0',
                ),
            )
            box.add_class('pdfv-page')
            self._page_images[index] = image
            self._page_boxes[index] = box
            order.append(box)
        self.pages_box.children = tuple(order)

    def _resize_pages(self) -> None:
        """Re-size the placeholders after a zoom or a frame-width change."""
        dpi = self.dpi()
        for index, box in self._page_boxes.items():
            width_pt, height_pt = self._sizes[index]
            width_px, height_px = page_pixel_size(width_pt, height_pt, dpi)
            for node in (box, self._page_images[index]):
                node.layout.width = f'{width_px}px'
                node.layout.height = f'{height_px}px'
        # Every rendered page is now the wrong size; drop and redraw.
        for index in list(self._filled):
            self._page_images[index].value = b''
        self._filled = set()

    def page_image(self, index: int) -> Optional[widgets.Image]:
        """The widget holding page ``index``, for tests and callers."""
        return self._page_images.get(int(index))

    def _fill_page(self, index: int) -> None:
        if self._doc is None or index not in self._page_images:
            return
        page_hits = [h.rect for h in self._result.hits if h.page == index]
        active = None
        if 0 <= self._hit < self._result.n_hits:
            current = self._result.hits[self._hit]
            if current.page == index:
                active = current.rect
                page_hits = [r for r in page_hits if r != active]
        try:
            self._page_images[index].value = render_page_png(
                self._doc, index, self.dpi(),
                highlights=page_hits, active=active,
            )
            self._filled.add(index)
        except Exception as exc:
            # A page that cannot be drawn is one damaged page, not a broken
            # viewer -- the rest of the document stays reachable.
            self._page_images[index].value = b''
            self._filled.discard(index)
            self.page_status.value = (
                f'<span style="color:#d32f2f;">{_html.escape(str(exc))}</span>'
            )

    def _free_outside(self, keep) -> None:
        """Let go of pages nobody is looking at.

        A rendered page is a PNG in the kernel and a copy in the browser; a
        few hundred of them is how a viewer turns into a memory leak.
        """
        keep = set(keep)
        for index in list(self._filled):
            if index not in keep:
                self._page_images[index].value = b''
                self._filled.discard(index)

    def _visible_window(self, centre: int) -> List[int]:
        if not self._continuous:
            return [self._page]
        first = max(0, centre - RENDER_WINDOW)
        last = min(len(self._sizes) - 1, centre + RENDER_WINDOW)
        return list(range(first, last + 1))

    # -- browser bridge ------------------------------------------------------

    def _emit(self, script: str) -> None:
        if callable(self._run_js) and script.strip():
            self._run_js(script)

    def _on_bridge(self, _button=None) -> None:
        """One message from the page stack in the browser."""
        raw = self._bridge_input.value or ''
        self._bridge_input.value = ''
        if not raw or self._doc is None:
            return
        try:
            message = json.loads(raw)
        except Exception:
            return
        if message.get('token') != self._token:
            return  # from a document that is no longer on screen

        action = message.get('action')
        if action == 'width':
            try:
                width = int(message.get('px') or 0)
            except (TypeError, ValueError):
                return
            if width <= 0 or abs(width - self._frame_px) < 24:
                return  # not a change worth re-rasterising the document for
            self._frame_px = width
            if self._zoom == FIT_WIDTH:
                self._resize_pages()
                self._render()
            return

        if action == 'pages':
            wanted = message.get('want') or []
            try:
                wanted = [int(p) for p in wanted]
            except (TypeError, ValueError):
                return
            wanted = [p for p in wanted if 0 <= p < len(self._page_images)]
            if not wanted:
                return
            # Which page the user is on follows what they are looking at, so
            # the page number in the toolbar keeps up while they scroll.
            self._page = clamp_page(wanted[0], self.total_pages)
            self._sync_page_controls()
            keep = set()
            for page in wanted:
                keep.update(self._visible_window(page))
            self._free_outside(keep)
            for page in wanted:
                if page not in self._filled:
                    self._fill_page(page)
            return

    def _install_observer(self) -> None:
        """Ask the browser to report what is on screen, and how wide it is.

        Only the scrolling view needs it: the one-page view shows what it
        was told to show and never has to ask.
        """
        if self._continuous:
            self._emit(_pages_js(self._token))

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
        self._sizes = []
        self._filled = set()
        self._page_images = {}
        self._page_boxes = {}
        self.pages_box.children = ()
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
        if not self._continuous:
            # One page on screen: going to another one replaces it.
            self._build_pages()
            self._render()
            return
        self._render()
        # The pages are all on screen already, so going to one means
        # scrolling to it rather than replacing what is shown.
        self._emit(_goto_js(self._token, target))

    def step_page(self, delta: int) -> None:
        self.goto_page(self._page + int(delta))

    def set_zoom(self, value: Any) -> None:
        """Fit to width, or one of the fixed steps."""
        try:
            value = int(value)
        except (TypeError, ValueError):
            value = FIT_WIDTH
        self._zoom = FIT_WIDTH if value == FIT_WIDTH else clamp_dpi_index(value)
        self._sync_page_controls()
        if self._doc is None:
            return
        self._resize_pages()
        self._render()
        if self._continuous:
            self._emit(_goto_js(self._token, self._page))

    def set_zoom_index(self, index: int) -> None:
        """Kept for callers that only know the fixed steps."""
        self.set_zoom(clamp_dpi_index(index))

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
        """Draw the pages around the current one.

        The browser refines this as soon as it reports what is actually on
        screen; this is what makes the first paint show something.
        """
        if self._doc is None:
            return
        window = self._visible_window(self._page)
        self._free_outside(window)
        for index in window:
            self._fill_page(index)

    # -- widget sync ---------------------------------------------------------

    def _sync_page_controls(self) -> None:
        total = max(1, self.total_pages)
        self._syncing = True
        try:
            self.page_input.max = total
            self.page_input.value = self._page + 1
        finally:
            self._syncing = False
        self.page_total.value = f'of {total}'
        self.prev_page_btn.disabled = self._page <= 0
        self.next_page_btn.disabled = self._page >= total - 1
        self.page_status.value = page_status_html(
            self._page, total, self.dpi())

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
        self.set_zoom(change.get('new'))

    def _on_hit_select(self, change) -> None:
        if self._syncing or self._doc is None:
            return
        value = change.get('new')
        if value is None:
            return
        self._goto_hit(int(value))
