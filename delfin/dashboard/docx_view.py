"""Word documents in the dashboard file browser: show, search, edit.

One renderer serves all three. That is the point of the module: the view
and the editor produce the same DOM, so what the search jumps to is what
the user sees, and what they type into is a block whose address in the
document is written on it.

The previous view converted the file with mammoth. Good HTML, but nothing
in it says which paragraph of the document a line came from, so it could
carry a view and nothing else. Here every block is rendered with its
address -- ``p:12`` for the twelfth body paragraph, ``t:1:3:2:0`` for the
first paragraph of a table cell -- and an edit is written straight back to
that paragraph.

Writing goes through :func:`delfin.agent.office._splice_runs`, the same
splice the template filler uses: the new text is placed into the runs that
are already there instead of rebuilding the paragraph, so bold, spacing and
font changes inside the line survive an edit that did not touch them.

What a round-trip preserves: python-docx writes back the package it read
and leaves untouched everything it was not asked to change -- headers and
footers, numbering, charts, text boxes, comments, tracked changes. That is
why editing text here is safe in a way that editing a workbook is not, and
why the warning a spreadsheet needs has no counterpart in this module.

Emphasis and paragraph style can be set. A heading is made by naming the
document's own style, never by turning text big and bold: the style is
what carries the look and what the table of contents reads.

What is deliberately not offered: adding or removing paragraphs, and
editing anything outside the body. Those shift the addresses every block
is written back through, and a control that half-does that would be worse
than none.
"""

from __future__ import annotations

import base64
import html as _html
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

WORD_SUFFIXES = ('.docx',)

# A document this far past the comm channel's comfort would arrive as a
# blank page after a long wait. Reported rather than truncated silently.
MAX_FILE_BYTES = 40 * 1024 * 1024
MAX_BLOCKS = 4000
MAX_IMAGE_BYTES = 4 * 1024 * 1024

# Paragraph styles that carry meaning worth keeping in the view. Anything
# else renders as body text -- the point is legibility, not a style engine.
_HEADING_RE = re.compile(r'^(?:heading|überschrift)\s*([1-6])$', re.IGNORECASE)
_LIST_RE = re.compile(r'(list|liste|aufz)', re.IGNORECASE)


class DocxError(RuntimeError):
    """A document could not be read or written, with a reason to show."""


@dataclass
class Block:
    """One editable paragraph, and where it lives in the document."""

    address: str
    text: str
    level: int = 0            # heading level, 0 for body text
    listed: bool = False
    table: Optional[Tuple[int, int, int]] = None   # table, row, cell
    # (text, bold, italic, underline) per run. What makes the view look like
    # the document rather than like a text dump of it; an edit still works on
    # the paragraph's plain text, and the splice on write puts the changed
    # span back into the run it belongs to.
    runs: List[Tuple[str, bool, bool, bool]] = field(default_factory=list)

    @property
    def in_table(self) -> bool:
        return self.table is not None


@dataclass
class DocxDocument:
    """What the browser needs to show and search a document."""

    path: Path
    blocks: List[Block] = field(default_factory=list)
    images: Dict[str, str] = field(default_factory=dict)
    tables: int = 0
    notes: List[str] = field(default_factory=list)

    @property
    def text(self) -> str:
        """The searchable buffer, matching the DOM's textContent exactly."""
        return '\n'.join(block.text for block in self.blocks)


def _docx():
    try:
        import docx  # python-docx
    except ImportError as exc:  # pragma: no cover - depends on the install
        raise DocxError(
            'Word support is missing. Install it with: pip install python-docx'
        ) from exc
    return docx


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def _runs_of(paragraph) -> List[Tuple[str, bool, bool, bool]]:
    """The paragraph's runs with the emphasis each one carries."""
    out: List[Tuple[str, bool, bool, bool]] = []
    for run in getattr(paragraph, 'runs', []) or []:
        text = run.text or ''
        if not text:
            continue
        out.append((text, bool(run.bold), bool(run.italic), bool(run.underline)))
    return out


def _style_of(paragraph) -> Tuple[int, bool]:
    """Heading level and whether the paragraph is a list item."""
    name = ''
    try:
        name = str(paragraph.style.name or '')
    except Exception:
        pass
    match = _HEADING_RE.match(name.strip())
    if match:
        return int(match.group(1)), False
    if name.strip().lower() in ('title', 'titel'):
        return 1, False
    return 0, bool(_LIST_RE.search(name))


def _cell_paragraphs(document) -> List[Tuple[str, Any, Tuple[int, int, int]]]:
    out = []
    for t_index, table in enumerate(getattr(document, 'tables', []) or []):
        for r_index, row in enumerate(table.rows):
            for c_index, cell in enumerate(row.cells):
                for p_index, para in enumerate(cell.paragraphs):
                    address = f't:{t_index}:{r_index}:{c_index}:{p_index}'
                    out.append((address, para, (t_index, r_index, c_index)))
    return out


def read_document(path) -> DocxDocument:
    """Read a .docx into addressed blocks.

    Body paragraphs first, then table cells, which is the order the view
    renders them in. Every block carries the address an edit is written
    back to, so the two can never drift apart.
    """
    path = Path(path)
    if path.suffix.lower() not in WORD_SUFFIXES:
        raise DocxError(
            f'{path.name}: only .docx is supported. The older .doc format '
            'has to be converted first.')
    try:
        size = path.stat().st_size
    except OSError as exc:
        raise DocxError(f'File cannot be read: {exc}') from exc
    if size > MAX_FILE_BYTES:
        raise DocxError(
            f'Document is {size / (1024 * 1024):.0f} MB; the viewer is '
            f'limited to {MAX_FILE_BYTES // (1024 * 1024)} MB.')

    docx = _docx()
    try:
        document = docx.Document(str(path))
    except Exception as exc:
        raise DocxError(f'Document could not be opened: {exc}') from exc

    result = DocxDocument(path=path)
    blocks: List[Block] = []
    for index, para in enumerate(document.paragraphs):
        level, listed = _style_of(para)
        blocks.append(Block(address=f'p:{index}', text=para.text,
                            level=level, listed=listed,
                            runs=_runs_of(para)))
    for address, para, where in _cell_paragraphs(document):
        level, listed = _style_of(para)
        blocks.append(Block(address=address, text=para.text,
                            level=level, listed=listed, table=where,
                            runs=_runs_of(para)))

    if len(blocks) > MAX_BLOCKS:
        result.notes.append(
            f'Document has {len(blocks)} paragraphs; the first {MAX_BLOCKS} '
            'are shown. Editing is therefore switched off here.')
        blocks = blocks[:MAX_BLOCKS]

    result.blocks = blocks
    result.tables = len(getattr(document, 'tables', []) or [])
    result.images = _inline_images(document)
    if not any(block.text.strip() for block in blocks):
        result.notes.append('The document holds no readable text.')
    return result


def _inline_images(document) -> Dict[str, str]:
    """Embedded images as data URIs, keyed by relationship id.

    Skipped rather than shrunk when oversized: the browser gets the page
    either way, and a hundred megabytes of base64 over the comm channel
    would cost the whole document.
    """
    images: Dict[str, str] = {}
    try:
        parts = document.part.related_parts
    except Exception:
        return images
    for rel_id, part in list(parts.items()):
        content_type = str(getattr(part, 'content_type', '') or '')
        if not content_type.startswith('image/'):
            continue
        try:
            blob = part.blob
        except Exception:
            continue
        if not blob or len(blob) > MAX_IMAGE_BYTES:
            continue
        encoded = base64.b64encode(blob).decode('ascii')
        images[str(rel_id)] = f'data:{content_type};base64,{encoded}'
    return images


def is_editable(document: DocxDocument) -> bool:
    """A truncated document is shown but not edited: an edit would address
    a paragraph by an index the view no longer agrees with."""
    return not any('Editing is therefore switched off' in note
                   for note in document.notes)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _run_html(runs, fallback: str) -> str:
    """The paragraph's text, with the emphasis it carries.

    Wrapped in spans rather than <b>/<i>, so that what the block's
    innerText reports back is exactly the paragraph text -- the addresses
    and the offsets an edit and a search jump use are both counted in that
    text.
    """
    if not runs:
        return _html.escape(fallback) or '&nbsp;'
    out = []
    for text, bold, italic, underline in runs:
        styles = []
        if bold:
            styles.append('font-weight:700')
        if italic:
            styles.append('font-style:italic')
        if underline:
            styles.append('text-decoration:underline')
        escaped = _html.escape(text)
        out.append(f'<span style="{";".join(styles)}">{escaped}</span>'
                   if styles else escaped)
    return ''.join(out) or '&nbsp;'


def _block_html(block: Block, editable: bool) -> str:
    text = _run_html(block.runs, block.text)
    classes = ['dw-b']
    if block.level:
        classes.append(f'dw-h{block.level}')
    if block.listed:
        classes.append('dw-li')
    if block.in_table:
        classes.append('dw-tc')
    attrs = (
        f' class="{" ".join(classes)}"'
        f' data-a="{_html.escape(block.address, quote=True)}"'
    )
    if editable:
        attrs += ' contenteditable="true" spellcheck="false"'
    return f'<div{attrs}>{text}</div>'


DOC_CSS = (
    '.dw-page { font-family: Calibri, Carlito, Arial, sans-serif;'
    ' font-size:15px; line-height:1.55; color:#1a1a1a; background:#fff;'
    ' padding:30px 34px; margin:0 auto; max-width:900px;'
    ' box-shadow:0 1px 4px rgba(0,0,0,0.14); box-sizing:border-box; }'
    '.dw-b { margin:0 0 9px 0; min-height:1em; white-space:pre-wrap;'
    ' overflow-wrap:anywhere; outline:none; border-radius:3px; }'
    '.dw-b[contenteditable="true"]:hover { background:#f4f8ff; }'
    '.dw-b[contenteditable="true"]:focus { background:#eef4ff;'
    ' box-shadow:0 0 0 2px #90caf9; }'
    '.dw-b.dw-dirty { background:#fff6cc; }'
    '.dw-h1 { font-size:1.7em; font-weight:600; color:#12447a;'
    ' margin:18px 0 9px 0; }'
    '.dw-h2 { font-size:1.4em; font-weight:600; color:#12447a;'
    ' margin:16px 0 8px 0; }'
    '.dw-h3 { font-size:1.18em; font-weight:600; margin:14px 0 6px 0; }'
    '.dw-h4, .dw-h5, .dw-h6 { font-size:1.05em; font-weight:600;'
    ' margin:12px 0 5px 0; }'
    '.dw-li { padding-left:22px; position:relative; }'
    '.dw-li:before { content:"•"; position:absolute; left:8px; color:#666; }'
    '.dw-table { border-collapse:collapse; margin:12px 0; width:100%; }'
    '.dw-table td { border:1px solid #d6d6d6; padding:6px 8px;'
    ' vertical-align:top; }'
    '.dw-img { max-width:100%; height:auto; margin:10px 0;'
    ' border:1px solid #e0e0e0; border-radius:3px; }'
    '.dw-match { background:#fff59d; }'
    '.dw-match.dw-current { background:#ffb74d; }'
    '.dw-bar { position:sticky; top:0; z-index:6; display:flex; gap:6px;'
    ' align-items:center; padding:7px 12px; margin:0 0 14px 0;'
    ' background:linear-gradient(#fdfdfe,#f1f4f7);'
    ' border:1px solid #d6dbe0; border-radius:5px;'
    ' box-shadow:0 1px 3px rgba(0,0,0,0.07); box-sizing:border-box; }'
    '.dw-btn { font-size:13px; min-width:30px; height:28px; cursor:pointer;'
    ' border:1px solid #c8ced4; border-radius:4px; background:#fff;'
    ' color:#1f2937; line-height:1; }'
    '.dw-btn:hover { border-color:#1565c0; color:#12447a; }'
    '.dw-btn:active { background:#e8f0fe; }'
    '.dw-bar-sep { width:1px; height:20px; background:#d6dbe0; margin:0 4px; }'
    '.dw-style { font-size:13px; height:28px; padding:0 6px;'
    ' border:1px solid #c8ced4; border-radius:4px; background:#fff;'
    ' color:#1f2937; min-width:150px; }'
    '.dw-style:hover { border-color:#1565c0; }'
)


def toolbar_html(current: str = 'Normal') -> str:
    """The formatting controls, shown while a document is being edited."""
    out = ['<div class="dw-bar">']
    out.append('<button class="dw-btn dw-b" title="Bold (Ctrl+B)"><b>B</b></button>')
    out.append('<button class="dw-btn dw-i" title="Italic (Ctrl+I)"><i>I</i></button>')
    out.append('<button class="dw-btn dw-u" title="Underline (Ctrl+U)">'
               '<u>U</u></button>')
    out.append('<span class="dw-bar-sep"></span>')
    out.append('<select class="dw-style" title="Paragraph style">')
    for label, code in PARAGRAPH_STYLES:
        chosen = ' selected' if code == current else ''
        out.append(f'<option value="{_html.escape(code, quote=True)}"{chosen}>'
                   f'{_html.escape(label)}</option>')
    out.append('</select>')
    out.append('</div>')
    return ''.join(out)


def render_html(document: DocxDocument, *, editable: bool = False) -> str:
    """The document as one block per paragraph, each carrying its address.

    Table cells are laid out as a table so the document still reads like
    one; the editable unit inside a cell is the same addressed block as
    everywhere else.
    """
    out: List[str] = []
    if editable:
        out.append(toolbar_html())
    out.append('<div class="dw-page">')
    for block in document.blocks:
        if block.in_table:
            continue
        out.append(_block_html(block, editable))

    by_table: Dict[int, Dict[Tuple[int, int], List[Block]]] = {}
    for block in document.blocks:
        if not block.in_table:
            continue
        t_index, r_index, c_index = block.table
        by_table.setdefault(t_index, {}).setdefault(
            (r_index, c_index), []).append(block)

    for t_index in sorted(by_table):
        cells = by_table[t_index]
        rows = sorted({r for r, _c in cells})
        out.append('<table class="dw-table"><tbody>')
        for r_index in rows:
            out.append('<tr>')
            for _r, c_index in sorted(k for k in cells if k[0] == r_index):
                out.append('<td>')
                for block in cells[(r_index, c_index)]:
                    out.append(_block_html(block, editable))
                out.append('</td>')
            out.append('</tr>')
        out.append('</tbody></table>')

    for uri in document.images.values():
        out.append(f'<img class="dw-img" src="{uri}" alt="">')

    out.append('</div>')
    return ''.join(out)


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

# The paragraph styles offered, by the name shown in the menu. Small on
# purpose: these are the ones a letter or a report is actually built from,
# and a document's own styles are what carry its look and its table of
# contents -- so a heading is made by naming the style, never by making the
# text big and bold.
PARAGRAPH_STYLES: Tuple[Tuple[str, str], ...] = (
    ('Body text', 'Normal'),
    ('Heading 1', 'Heading 1'),
    ('Heading 2', 'Heading 2'),
    ('Heading 3', 'Heading 3'),
    ('Bullet list', 'List Bullet'),
    ('Numbered list', 'List Number'),
)
_KNOWN_STYLES = {code for _label, code in PARAGRAPH_STYLES}


def check_style(name: Any) -> str:
    """A paragraph style this view offers."""
    text = str(name or '').strip()
    if text not in _KNOWN_STYLES:
        raise DocxError(f'{text!r} is not a paragraph style this view sets.')
    return text


def set_paragraph_style(document, address: str, style: str) -> None:
    """Give a paragraph one of the document's own styles.

    Named rather than imitated: a heading made by turning text big and bold
    is not a heading, and the table of contents will not have it.
    """
    paragraph = _paragraph_at(document, address)
    if paragraph is None:
        raise DocxError(f'There is no paragraph {address!r} in the document.')
    style = check_style(style)
    try:
        paragraph.style = document.styles[style]
    except KeyError as exc:
        raise DocxError(
            f'This document has no {style!r} style. It was made from a '
            'template that does not define one.') from exc


def _set_paragraph_runs(paragraph, runs: Sequence[Mapping[str, Any]]) -> None:
    """Rebuild a paragraph from runs that carry their own emphasis.

    Used only when the emphasis changed. A plain text edit goes through the
    splice instead, which keeps every run exactly as it was -- this path
    cannot, because the runs are being redrawn. What it can keep is the
    look the paragraph was written in, so the font, size and colour of the
    first run are carried onto all of them.
    """
    existing = list(paragraph.runs)
    base = existing[0] if existing else None
    for run in existing:
        run._element.getparent().remove(run._element)
    for spec in runs:
        text = str(spec.get('t') or '')
        if not text:
            continue
        run = paragraph.add_run(text)
        run.bold = bool(spec.get('b'))
        run.italic = bool(spec.get('i'))
        run.underline = bool(spec.get('u'))
        if base is not None:
            run.font.name = base.font.name
            run.font.size = base.font.size
            if base.font.color is not None and base.font.color.rgb is not None:
                run.font.color.rgb = base.font.color.rgb
    if not paragraph.runs:
        paragraph.add_run('')


def runs_differ(before: Sequence[Tuple[str, bool, bool, bool]],
                after: Sequence[Mapping[str, Any]]) -> bool:
    """Whether the emphasis changed, as opposed to only the text.

    Deciding this is what lets an ordinary edit keep the paragraph's runs
    untouched and only a formatting change redraw them.
    """
    was = [(bool(b), bool(i), bool(u)) for _t, b, i, u in before]
    now = [(bool(r.get('b')), bool(r.get('i')), bool(r.get('u')))
           for r in after if str(r.get('t') or '')]
    if len(was) != len(now):
        return True
    return was != now


def _paragraph_at(document, address: str):
    """Resolve a block address back to the paragraph it came from."""
    parts = str(address).split(':')
    try:
        if parts[0] == 'p' and len(parts) == 2:
            return document.paragraphs[int(parts[1])]
        if parts[0] == 't' and len(parts) == 5:
            _kind, t_index, r_index, c_index, p_index = parts
            table = document.tables[int(t_index)]
            cell = table.rows[int(r_index)].cells[int(c_index)]
            return cell.paragraphs[int(p_index)]
    except (IndexError, ValueError):
        return None
    return None


def changed_range(before: str, after: str) -> Tuple[int, int, int]:
    """The one span that differs, as ``(start, end_before, end_after)``.

    Splicing the whole paragraph would put all of its text into the first
    run and empty the others, which flattens every bold word and font
    change in the line -- including the ones the edit never touched. Only
    the span that actually changed is written, so the rest keeps the runs
    it already had.
    """
    start = 0
    limit = min(len(before), len(after))
    while start < limit and before[start] == after[start]:
        start += 1
    end_before, end_after = len(before), len(after)
    while (end_before > start and end_after > start
           and before[end_before - 1] == after[end_after - 1]):
        end_before -= 1
        end_after -= 1
    return start, end_before, end_after


def _set_paragraph_text(paragraph, text: str) -> None:
    """Put ``text`` into a paragraph without rebuilding it.

    A paragraph with no runs at all (an empty line) has nothing to splice
    into, so one run is added.
    """
    from delfin.agent.office import _splice_runs

    runs = list(paragraph.runs)
    if not runs:
        if text:
            paragraph.add_run(text)
        return
    current = ''.join(run.text for run in runs)
    if current == text:
        return
    if not current:
        # Runs that are all empty: there is no character to splice against,
        # and the insertion branch below would reach for current[0]. Word
        # leaves such runs behind routinely -- a deleted line, a stray
        # formatting mark -- so this is the ordinary case of typing into a
        # blank paragraph, not a broken document.
        runs[0].text = text
        return
    start, end_before, end_after = changed_range(current, text)
    if start == end_before and end_after > start:
        # A pure insertion is an empty range, and an empty range sits between
        # two runs rather than inside one: the splice skips every run and
        # writes nothing. Typing at the end of a paragraph is exactly that
        # case, so it silently saved the paragraph unchanged.
        #
        # Widen it onto one neighbouring character and carry that character
        # through. The text then lands in the run that character belongs to,
        # which is also where Word would put it.
        inserted = text[start:end_after]
        if start > 0:
            _splice_runs(runs, start - 1, start, current[start - 1] + inserted)
        else:
            _splice_runs(runs, 0, 1, inserted + current[0])
        return
    _splice_runs(runs, start, end_before, text[start:end_after])


def apply_edits(path, edits: Dict[str, str]) -> Dict[str, Any]:
    """Write block edits back into the document.

    Returns what was written and what was refused. An address that does
    not resolve is an error rather than a silent skip: it means the view
    and the file have drifted apart, and the user is owed that news
    before they close the tab believing the edit landed.
    """
    path = Path(path)
    if not isinstance(edits, dict) or not edits:
        raise DocxError('There is nothing to save.')

    docx = _docx()
    try:
        document = docx.Document(str(path))
    except Exception as exc:
        raise DocxError(f'Document could not be opened: {exc}') from exc

    written = 0
    unknown: List[str] = []
    for address, change in edits.items():
        paragraph = _paragraph_at(document, str(address))
        if paragraph is None:
            unknown.append(str(address))
            continue
        if isinstance(change, Mapping):
            style = change.get('style')
            if style:
                set_paragraph_style(document, str(address), style)
            runs = change.get('runs')
            if runs is not None:
                before = _runs_of(paragraph)
                if runs_differ(before, runs):
                    # The emphasis changed, so the runs have to be redrawn.
                    _set_paragraph_runs(paragraph, runs)
                else:
                    # Only the text did: the splice keeps every run as it is.
                    _set_paragraph_text(
                        paragraph, ''.join(str(r.get('t') or '') for r in runs))
            elif 'text' in change:
                _set_paragraph_text(paragraph, str(change.get('text') or ''))
        else:
            _set_paragraph_text(paragraph, '' if change is None else str(change))
        written += 1

    if unknown:
        raise DocxError(
            'These paragraphs are not in the document (any more): '
            + ', '.join(unknown[:5])
            + '. The file was left unchanged.')

    return {'document': document, 'written': written}


def save(document, path) -> None:
    """Write a python-docx document to ``path``."""
    try:
        document.save(str(path))
    except Exception as exc:
        raise DocxError(f'Saving failed: {exc}') from exc


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

@dataclass
class DocxHit:
    """One match, addressed so the view can jump to it."""

    block: int
    address: str
    start: int
    end: int


def search(document: DocxDocument, term: str, *, limit: int = 2000
           ) -> List[DocxHit]:
    """Every occurrence of ``term``, case-insensitively, block by block.

    Offsets are into the block's own text, not into a flattened document,
    so a jump lands on a character rather than near one.
    """
    term = str(term or '')
    if not term:
        return []
    needle = term.lower()
    hits: List[DocxHit] = []
    for index, block in enumerate(document.blocks):
        haystack = block.text.lower()
        start = haystack.find(needle)
        while start >= 0:
            hits.append(DocxHit(block=index, address=block.address,
                                start=start, end=start + len(needle)))
            if len(hits) >= limit:
                return hits
            start = haystack.find(needle, start + 1)
    return hits


def focus_js(scope_class: str, address: str, start: int, end: int) -> str:
    """Mark one hit and scroll it into view, without touching the rest.

    The match is wrapped with a Range rather than by rewriting innerHTML:
    a rewrite would flatten the formatting of the paragraph it landed in,
    and in a Word document that is the document.
    """
    return _FOCUS_JS_TEMPLATE.replace(
        '__SCOPE__', json.dumps(str(scope_class))
    ).replace(
        '__ADDR__', json.dumps(str(address))
    ).replace('__START__', str(int(start))).replace('__END__', str(int(end)))


def edit_js(scope_class: str) -> str:
    """Report a block to the kernel when the user leaves it.

    Reported on leaving rather than on every keystroke: a paragraph is a
    unit of text, and sending one message per character would put the
    kernel behind the typing.
    """
    return _EDIT_JS_TEMPLATE.replace('__SCOPE__', json.dumps(str(scope_class)))


def mark_saved_js(scope_class: str) -> str:
    """Clear the changed marks, without rebuilding anything.

    Rebuilding would take the caret out of the paragraph the user is
    still standing in, which is not what saving does.
    """
    return (
        "(function(){var r=document.querySelector('.' + "
        + json.dumps(str(scope_class))
        + ");if(!r)return;Array.prototype.forEach.call("
        "r.querySelectorAll('.dw-b.dw-dirty'),function(b){"
        "b.classList.remove('dw-dirty');});})();"
    )


_EDIT_JS_TEMPLATE = r"""
(function(){
  var SCOPE = __SCOPE__;
  var tries = 0;

  function boot(){
    var root = document.querySelector('.' + SCOPE);
    var page = root && root.querySelector('.dw-page');
    if (!page) { if (++tries < 40) setTimeout(boot, 50); return; }
    if (page.dataset.editBound === '1') return;
    page.dataset.editBound = '1';

    function send(payload){
      var field = root.querySelector('.calc-sheet-payload textarea')
               || root.querySelector('.calc-sheet-payload input');
      if (!field) return;
      var proto = field.tagName === 'TEXTAREA'
        ? HTMLTextAreaElement.prototype : HTMLInputElement.prototype;
      var desc = Object.getOwnPropertyDescriptor(proto, 'value');
      if (desc && desc.set) desc.set.call(field, JSON.stringify(payload));
      else field.value = JSON.stringify(payload);
      field.dispatchEvent(new Event('input', {bubbles: true}));
      field.dispatchEvent(new Event('change', {bubbles: true}));
      var holder = root.querySelector('.calc-sheet-action');
      var btn = holder && (holder.tagName === 'BUTTON'
        ? holder : holder.querySelector('button'));
      if (btn) btn.click();
    }

    /* The runs of a block, read back off the page. What is sent is the
       structure, not the text: an emphasis change has to arrive as one, or
       the paragraph would be written back with its old runs and the new
       words. */
    function runsOf(block){
      var runs = [];
      (function walk(node, bold, italic, under){
        for (var i = 0; i < node.childNodes.length; i++) {
          var child = node.childNodes[i];
          if (child.nodeType === 3) {
            if (child.nodeValue) {
              runs.push({t: child.nodeValue, b: bold ? 1 : 0,
                         i: italic ? 1 : 0, u: under ? 1 : 0});
            }
            continue;
          }
          if (child.nodeType !== 1) continue;
          var name = child.tagName;
          var style = child.style || {};
          var weight = style.fontWeight || '';
          walk(child,
               bold || name === 'B' || name === 'STRONG'
                 || weight === 'bold' || weight === '700' || weight === '600',
               italic || name === 'I' || name === 'EM'
                 || style.fontStyle === 'italic',
               under || name === 'U'
                 || (style.textDecoration || '').indexOf('underline') >= 0);
        }
      })(block, false, false, false);
      return runs;
    }

    function report(block, force){
      if (!block || block.getAttribute('contenteditable') !== 'true') return;
      if (!force && block.innerText === block.dataset.was) return;
      block.classList.add('dw-dirty');
      block.dataset.was = block.innerText;   /* sent; only resend on change */
      send({kind: 'docx', address: block.getAttribute('data-a'),
            runs: runsOf(block)});
    }

    /* The block being worked on, whether or not it still has focus: a
       toolbar button takes focus away the moment it is pressed. */
    var active = null;
    function currentBlock(){
      var node = document.activeElement;
      var block = node && node.closest ? node.closest('.dw-b') : null;
      return block || active;
    }

    function emphasise(command){
      var block = currentBlock();
      if (!block) return;
      block.focus();
      document.execCommand(command, false, null);
      report(block, true);
    }

    var bar = root.querySelector('.dw-bar');
    if (bar) {
      var press = {'.dw-b': 'bold', '.dw-i': 'italic', '.dw-u': 'underline'};
      Object.keys(press).forEach(function(sel){
        var button = bar.querySelector(sel);
        if (!button) return;
        /* mousedown, not click: by the time click fires, the selection in
           the paragraph is already gone. */
        button.addEventListener('mousedown', function(e){
          e.preventDefault();
          emphasise(press[sel]);
        });
      });
      var styleBox = bar.querySelector('.dw-style');
      if (styleBox) styleBox.addEventListener('change', function(){
        var block = currentBlock();
        if (!block) return;
        block.classList.add('dw-dirty');
        send({kind: 'docx', address: block.getAttribute('data-a'),
              style: styleBox.value});
        styleBox.blur();
      });
    }

    page.addEventListener('keydown', function(e){
      if (!(e.ctrlKey || e.metaKey)) return;
      var key = String(e.key || '').toLowerCase();
      var command = key === 'b' ? 'bold'
                  : key === 'i' ? 'italic'
                  : key === 'u' ? 'underline' : '';
      if (!command) return;
      e.preventDefault();
      emphasise(command);
    });

    page.addEventListener('focusin', function(e){
      var block = e.target.closest && e.target.closest('.dw-b');
      if (block) { block.dataset.was = block.innerText; active = block; }
    }, true);

    page.addEventListener('focusout', function(e){
      report(e.target.closest && e.target.closest('.dw-b'));
    }, true);

    /* Also while typing, once it pauses. Leaving the block is the precise
       moment, but it is not a reliable one: pressing Save is a click on a
       button, and whether that moves focus out of the paragraph first is up
       to the browser. When it does not, the edit was never sent and saving
       reported nothing to save while the text sat on screen. */
    var typing = null;
    page.addEventListener('input', function(e){
      var block = e.target.closest && e.target.closest('.dw-b');
      if (!block) return;
      if (typing) clearTimeout(typing);
      typing = setTimeout(function(){ typing = null; report(block); }, 350);
    });

    /* Enter would start a new paragraph, and a paragraph this view cannot
       address is a paragraph an edit cannot be written back to. */
    page.addEventListener('keydown', function(e){
      if (e.key !== 'Enter' || e.shiftKey) return;
      var block = e.target.closest && e.target.closest('.dw-b');
      if (!block) return;
      e.preventDefault();
      block.blur();
    });
  }

  boot();
})();
"""


_FOCUS_JS_TEMPLATE = r"""
(function(){
  var root = document.querySelector('.' + __SCOPE__);
  if (!root) return;
  var page = root.querySelector('.dw-page');
  if (!page) return;

  var old = page.querySelector('.dw-match.dw-current');
  if (old && old.parentNode) {
    var parent = old.parentNode;
    parent.replaceChild(document.createTextNode(old.textContent || ''), old);
    parent.normalize();
  }

  var block = page.querySelector('[data-a="' + __ADDR__ + '"]');
  if (!block) return;

  /* Walk to the character offset: a paragraph is several text nodes once
     it has any formatting, so an offset into its text is not an offset
     into its first node. */
  function locate(node, index){
    var walker = document.createTreeWalker(node, NodeFilter.SHOW_TEXT, null);
    var seen = 0, current = walker.nextNode();
    while (current) {
      var length = (current.nodeValue || '').length;
      if (seen + length >= index) return {node: current, offset: index - seen};
      seen += length;
      current = walker.nextNode();
    }
    return null;
  }

  var from = locate(block, __START__);
  var to = locate(block, __END__);
  if (!from || !to) return;

  var range = document.createRange();
  try {
    range.setStart(from.node, from.offset);
    range.setEnd(to.node, to.offset);
    var mark = document.createElement('mark');
    mark.className = 'dw-match dw-current';
    range.surroundContents(mark);
    mark.scrollIntoView({block: 'center'});
  } catch (_err) {
    /* A match spanning several runs cannot be wrapped in one element;
       showing the paragraph is still better than not moving at all. */
    block.scrollIntoView({block: 'center'});
  }
})();
"""
