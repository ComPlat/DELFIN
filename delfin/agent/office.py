"""Document primitives for office work: spreadsheets, PDFs, Word files.

Until now the agent could only reach these files through ``bash`` plus
generated Python. That works, but every failure mode has to be
rediscovered by the model on the spot, and several of them are silent:

* ``read_file`` decodes an ``.xlsx`` as UTF-8 and returns the raw ZIP
  container as mojibake — thousands of tokens of noise that look like
  data.
* A workbook written by a library carries formulas but **no cached
  values**. Reading it with ``data_only=True`` yields ``None`` for every
  computed cell, which reads as "the column is empty" rather than "the
  numbers were never evaluated".
* A PDF form filled without ``/NeedAppearances`` renders empty in
  viewers that do not build appearance streams themselves. The values
  are in the file; the printout is blank.
* Check boxes do not take ``true``. They take the on-state name declared
  by that particular field — ``/Yes``, ``/On``, ``/1`` — and silently
  stay off for anything else.
* XFA forms ignore AcroForm writes, so a "successful" fill changes
  nothing the user will see.
* A merge that opens its inputs while it writes produces a document
  that looks complete when the fourth attachment turns out to be
  encrypted. The pages missing from it are the ones nobody can name
  afterwards, so every input is opened before the first byte is
  written.
* An .ods cell repeats: LibreOffice writes one cell with
  ``number-columns-repeated="1013"`` where a naive reader sees one
  column, and every value after it lands under the wrong heading.
* A page with no text layer is a photograph of paper, not an empty
  page. Saying "no text" reports the file wrongly, and OCR of it is a
  reading of pixels rather than content of the document — so it is
  labelled as such and never overwrites a page that has real text.

Each of those is handled once, here, instead of in whatever code the
model happens to write that turn. The functions are pure in the sense
that matters for testing: they take paths and values, they touch only
the file named, and they raise only :class:`OfficeError`.

Where there is no reader, there is a refusal that names the way out
rather than a parser that copes. A legacy .xls read by an approximate
BIFF reader, or an HTML export that arrived under an .xls name, gives
back rows that look exactly like data — and an administrative file read
wrongly is worse than one not read at all, because nothing downstream
can tell. Those refusals name what the bytes actually are and the one
command that converts them.

Writing is deliberately narrow. ``edit_sheet`` changes cells and appends
rows; it does not restructure workbooks. Rewriting a spreadsheet through
``openpyxl`` reconstructs the whole container, so anything the library
does not model is at risk — the caller is told which fragile features a
workbook carries before it is rewritten, and the tool layer captures a
byte-exact pre-image so the change can be undone.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Optional

# Reading caps. A spreadsheet can be a million rows; the point of the
# read path is to let the model see the shape and the region it asked
# for, not to pour the file into the context window.
DEFAULT_MAX_ROWS = 200
DEFAULT_MAX_COLS = 40
MAX_TEXT_CHARS = 20000

SPREADSHEET_SUFFIXES = frozenset({".xlsx", ".xlsm", ".xltx", ".xltm"})
CSV_SUFFIXES = frozenset({".csv", ".tsv"})
PDF_SUFFIXES = frozenset({".pdf"})
WORD_SUFFIXES = frozenset({".docx"})
# OpenDocument. LibreOffice is the office suite of a great many public
# authorities and universities, so .ods and .odt arrive as ordinary
# attachments rather than as exotica. Reading only: see
# _CONVERSION_ROUTES for why writing them is refused instead of faked.
ODS_SUFFIXES = frozenset({".ods"})
ODT_SUFFIXES = frozenset({".odt"})

DOCUMENT_SUFFIXES = (
    SPREADSHEET_SUFFIXES | CSV_SUFFIXES | PDF_SUFFIXES | WORD_SUFFIXES
    | ODS_SUFFIXES | ODT_SUFFIXES
)

# Formats there is no reader for, mapped to the format a conversion has
# to produce. Deliberately NOT readers: every one of these has a parser
# somewhere that will return rows for most files and mangle the rest,
# and a wrong table out of an administrative file is worse than no
# table, because nothing about it looks wrong afterwards.
_CONVERSION_ROUTES: dict[str, tuple[str, str]] = {
    ".xls": ("legacy Excel workbook (BIFF, Excel 97-2003)", "xlsx"),
    ".xlt": ("legacy Excel template (BIFF, Excel 97-2003)", "xlsx"),
    ".doc": ("legacy Word document (Word 97-2003)", "docx"),
}

# Suffix -> (import name, distribution name). The distribution name is
# what the user has to install, and it differs from the import name for
# python-docx, which is the classic reason a "just pip install docx"
# hint sends someone to an unrelated package.
_BACKENDS: dict[str, tuple[str, str]] = {
    "spreadsheet": ("openpyxl", "openpyxl"),
    "pdf": ("pypdf", "pypdf"),
    "word": ("docx", "python-docx"),
    "csv": ("csv", ""),  # standard library
    # Reading a PDF and producing one are different libraries: pypdf
    # takes existing pages apart and puts them together, but it does not
    # lay out text. Writing a page needs reportlab.
    "pdf_write": ("reportlab", "reportlab"),
    # OpenDocument (.ods/.odt). The import name is 'odf' and the package
    # is 'odfpy' — the same trap python-docx sits in, which is why the
    # distribution name is carried separately.
    "opendocument": ("odf", "odfpy"),
}


class OfficeError(RuntimeError):
    """A document operation failed for a reason worth reporting verbatim."""


# ---------------------------------------------------------------------------
# Backends and file kinds
# ---------------------------------------------------------------------------

def document_kind(path: Any) -> Optional[str]:
    """The handling family of *path*, or None if there is no reader.

    ``spreadsheet`` / ``csv`` / ``pdf`` / ``word`` /
    ``opendocument_sheet`` / ``opendocument_text``. The OpenDocument
    kinds are separate from ``spreadsheet`` and ``word`` on purpose:
    they can be read but not written, and letting them share a kind
    would send an .ods into openpyxl, which fails with a message about
    a broken ZIP rather than about the format.
    """
    suffix = Path(str(path)).suffix.lower()
    if suffix in SPREADSHEET_SUFFIXES:
        return "spreadsheet"
    if suffix in CSV_SUFFIXES:
        return "csv"
    if suffix in PDF_SUFFIXES:
        return "pdf"
    if suffix in WORD_SUFFIXES:
        return "word"
    if suffix in ODS_SUFFIXES:
        return "opendocument_sheet"
    if suffix in ODT_SUFFIXES:
        return "opendocument_text"
    return None


def backend_available(kind: str) -> bool:
    """True if the library backing *kind* can be imported right now."""
    entry = _BACKENDS.get(kind)
    if entry is None:
        return False
    import importlib
    try:
        importlib.import_module(entry[0])
        return True
    except Exception:
        return False


def available_backends() -> dict[str, bool]:
    """Availability of every backend, for diagnostics and tool gating."""
    return {kind: backend_available(kind) for kind in _BACKENDS}


def have_office_support() -> bool:
    """True if at least one non-stdlib document backend is importable.

    The tool surface is gated on this: advertising a spreadsheet tool on
    a machine without ``openpyxl`` only buys a tool call that fails.
    """
    return any(
        backend_available(kind)
        for kind in ("spreadsheet", "pdf", "word", "opendocument")
    )


def _require(kind: str) -> Any:
    """Import the backend for *kind* or raise a fixable error."""
    entry = _BACKENDS.get(kind)
    if entry is None:
        raise OfficeError(f"unknown document kind {kind!r}")
    import importlib
    try:
        return importlib.import_module(entry[0])
    except Exception as exc:
        dist = entry[1] or entry[0]
        raise OfficeError(
            f"{kind} support needs the '{dist}' package, which is not "
            f"installed ({exc}). Install it with "
            f"'pip install delfin[office]' or 'pip install {dist}'."
        ) from exc


def _resolve(path: Any, *, must_exist: bool = True) -> Path:
    p = Path(str(path)).expanduser()
    if must_exist and not p.exists():
        raise OfficeError(f"file not found: {p}")
    if must_exist and p.is_dir():
        raise OfficeError(f"{p} is a directory, not a document")
    return p


def column_letter(index: int) -> str:
    """1-based column index to spreadsheet letters (1 -> A, 27 -> AA).

    Local rather than ``openpyxl.utils`` so the CSV path stays usable
    without the spreadsheet backend installed.
    """
    if index < 1:
        raise OfficeError(f"column index must be >= 1, got {index}")
    out = ""
    while index:
        index, rem = divmod(index - 1, 26)
        out = chr(ord("A") + rem) + out
    return out


_CELL_RE = re.compile(r"^\s*([A-Za-z]{1,3})\s*([1-9][0-9]{0,6})\s*$")


def parse_cell(ref: str) -> tuple[int, int]:
    """``"B7"`` -> ``(row, column)`` = ``(7, 2)``, both 1-based."""
    m = _CELL_RE.match(str(ref or ""))
    if not m:
        raise OfficeError(
            f"{ref!r} is not a cell reference — expected a form like 'B7'")
    letters, digits = m.group(1).upper(), m.group(2)
    col = 0
    for ch in letters:
        col = col * 26 + (ord(ch) - ord("A") + 1)
    return int(digits), col


# ---------------------------------------------------------------------------
# What a file really is, and what to do when it is not readable here
# ---------------------------------------------------------------------------
#
# A suffix is a claim, not evidence. Two cases turn up constantly in
# administrative mail and both end in a wrong table rather than an
# error: a report exported by a web system arrives as an HTML table
# named ``.xls``, and a workbook someone "saved as .ods" is an .xlsx
# underneath. Parsers meet those with a best effort — a row count, some
# plausible values — so the refusals below name what the bytes actually
# are instead of what the name says.

_CONTAINER_LABELS: dict[str, str] = {
    "ole2": "a legacy Microsoft Office document (OLE2/BIFF)",
    "xlsx": "an .xlsx workbook (Office Open XML)",
    "docx": "a .docx document (Office Open XML)",
    "ods": "an OpenDocument spreadsheet",
    "odt": "an OpenDocument text document",
    "odf": "an OpenDocument file",
    "zip": "a ZIP archive",
    "pdf": "a PDF",
    "html": "an HTML page or table",
    "xml": "an XML file",
    "text": "plain text",
}


def _zip_container_kind(path: Path) -> str:
    """Which OOXML/OpenDocument flavour a ZIP holds, or ``"zip"``."""
    import zipfile
    try:
        with zipfile.ZipFile(path) as archive:
            names = set(archive.namelist())
            if "mimetype" in names:
                try:
                    mime = archive.read("mimetype").decode("ascii", "replace")
                except Exception:
                    mime = ""
                if mime.endswith("spreadsheet"):
                    return "ods"
                if mime.endswith("text"):
                    return "odt"
                if "opendocument" in mime:
                    return "odf"
            if any(n.startswith("xl/") for n in names):
                return "xlsx"
            if any(n.startswith("word/") for n in names):
                return "docx"
    except Exception:
        return "zip"
    return "zip"


def container_kind(path: Any) -> str:
    """What *path* holds, judged by its bytes. ``""`` when undecidable.

    Cheap on purpose: the first few kilobytes decide, because this runs
    on the refusal path where the alternative is a wrong answer, not on
    a hot loop.
    """
    p = Path(str(path))
    try:
        with open(p, "rb") as fh:
            head = fh.read(4096)
    except OSError:
        return ""
    # Four bytes, not the full eight of the OLE2 signature: a file that
    # was truncated in transit still has to be named for what it was,
    # rather than falling through to "plain text".
    if head.startswith(b"\xd0\xcf\x11\xe0"):
        return "ole2"
    if head.startswith(b"%PDF"):
        return "pdf"
    if head.startswith(b"PK\x03\x04"):
        return _zip_container_kind(p)
    probe = head.lstrip()[:512].lower()
    if (probe.startswith(b"<!doctype html") or probe.startswith(b"<html")
            or b"<table" in probe):
        return "html"
    if probe.startswith(b"<?xml"):
        return "xml"
    if b"\x00" in head:
        return ""
    return "text" if head else ""


_LIBREOFFICE_HINT = (
    "libreoffice --headless --convert-to {target} --outdir <directory> "
    "\"{name}\"")


def _mismatch_note(path: Path, expected: tuple[str, ...]) -> str:
    """A sentence naming the real format when it is not the expected one.

    Empty when the bytes match the suffix or cannot be judged, so a
    caller can append it unconditionally.
    """
    found = container_kind(path)
    if not found or found in expected:
        return ""
    label = _CONTAINER_LABELS.get(found, f"a {found} file")
    return (
        f" Note that {path.name} does not hold what its name claims: the "
        f"bytes are {label}. Exports from web systems are routinely named "
        "for a format they are not, so convert or rename it to match what "
        "it really is before reading it.")


def _refuse_unreadable_format(path: Path) -> None:
    """Raise for a format with no reader, naming the conversion route.

    Reached before any parser sees the file. The message has to carry
    the whole way out, because the caller cannot look the command up:
    which format to convert to, and one command that does it.
    """
    suffix = path.suffix.lower()
    entry = _CONVERSION_ROUTES.get(suffix)
    if entry is None:
        return
    label, target = entry
    expected = ("ole2",) if suffix in (".xls", ".xlt", ".doc") else ()
    raise OfficeError(
        f"{path.name} is a {label}. There is no reader for it here, and "
        "nothing guesses at the bytes: a legacy Office file taken apart by "
        "the wrong parser gives back content that looks right and is not. "
        f"Convert it to .{target} once and read the result — "
        + _LIBREOFFICE_HINT.format(target=target, name=path.name)
        + f" — or open it and use 'Save as' with the .{target} format."
        + _mismatch_note(path, expected))


def _refuse_unwritable_format(path: Path, *, what: str) -> None:
    """Raise for a format this module reads but cannot write."""
    suffix = path.suffix.lower()
    if suffix in _CONVERSION_ROUTES:
        label, target = _CONVERSION_ROUTES[suffix]
        raise OfficeError(
            f"{path.name} is a {label}, which is neither read nor written "
            f"here. Work in .{target} instead: convert the original once — "
            + _LIBREOFFICE_HINT.format(target=target, name=path.name)
            + f" — and produce a .{target} as the result.")
    kind = document_kind(path)
    if kind == "opendocument_sheet":
        raise OfficeError(
            f"{path.name} is an OpenDocument spreadsheet. It can be read "
            f"here but not written, so {what} would either fail or produce "
            "a file that is no longer an .ods. Convert it to .xlsx first — "
            + _LIBREOFFICE_HINT.format(target="xlsx", name=path.name)
            + " — change that, and convert back if the .ods is what has to "
            "be delivered.")
    if kind == "opendocument_text":
        raise OfficeError(
            f"{path.name} is an OpenDocument text document. It can be read "
            f"here but not written, so {what} is not possible. Convert it "
            "to .docx first — "
            + _LIBREOFFICE_HINT.format(target="docx", name=path.name)
            + " — and work on that.")
    _refuse_unreadable_format(path)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def render_grid(
    rows: list[list[Any]],
    *,
    first_row: int = 1,
    first_col: int = 1,
    max_cell_chars: int = 40,
) -> str:
    """Render cells as an addressable text grid.

    The column letters and row numbers are part of the output on purpose:
    the model has to be able to name a cell (``Daten!B7``) in a follow-up
    edit, and it can only do that if it saw the coordinates.
    """
    if not rows:
        return "(empty)"
    width = max(len(r) for r in rows)
    cells = [
        [_fmt(r[i])[:max_cell_chars] if i < len(r) else "" for i in range(width)]
        for r in rows
    ]
    headers = [column_letter(first_col + i) for i in range(width)]
    row_label_w = max(len(str(first_row + len(rows) - 1)), 3)
    col_w = [
        max(len(headers[i]), *(len(c[i]) for c in cells)) for i in range(width)
    ]
    out = [
        " " * row_label_w + " | "
        + " | ".join(headers[i].ljust(col_w[i]) for i in range(width))
    ]
    out.append("-" * len(out[0]))
    for n, row in enumerate(cells):
        label = str(first_row + n).rjust(row_label_w)
        out.append(
            label + " | "
            + " | ".join(row[i].ljust(col_w[i]) for i in range(width))
        )
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Values: numbers and dates as they are actually written
# ---------------------------------------------------------------------------
#
# A cell reading "1.234,50" is a string. Handed to arithmetic as-is it
# either raises or, worse, parses as 1.234 — the thousands separator
# read as a decimal point, off by a factor of a thousand, with no error
# anywhere. The same trap sits in dates: 03.04.2026 and 04/03/2026 are
# the same day written under two conventions, and guessing per value
# gets it wrong half the time.
#
# The way out is to decide per COLUMN rather than per value. A single
# "1.234" is genuinely ambiguous; a column containing it alongside one
# value with both separators, or one with two dots, or one with a
# non-three-digit group, is not. So: gather the evidence over the whole
# column, and where the column stays ambiguous, say so instead of
# picking the reading that happens to look plausible.

# Noise a number may legitimately carry: whitespace (including the
# non-breaking and narrow kinds a spreadsheet emits), currency symbols and
# a currency code beside it. Nothing else — and that is the point. The
# earlier rule stripped every character that was not a digit or separator,
# which turned the document reference "R-001" into "-001" and read it as
# the number -1. A whole column of references profiled as numeric, and in
# a reconciliation "R-001" and "X-001" both became -1 and compared EQUAL:
# two different records reported as agreeing.
_CURRENCY_CODE_RE = re.compile(r"(?i)\b(?:eur|usd|chf|gbp|sfr)\b")
_NUM_NOISE_RE = re.compile(r"[\s\u00a0\u202f\u2009€$£¥']")
# Any run of digits, optionally broken by separators. Deliberately NOT
# "1-3 digits then groups": that shape is right for 1.234,50 and wrong
# for 1234.50, and rejecting the ungrouped form makes a plain numeric
# column read as text — which then compares by spelling rather than by
# value.
_NUMERIC_SHAPE_RE = re.compile(r"^[-+]?[0-9]+(?:[.,][0-9]+)*$")
_DATE_SEP_RE = re.compile(r"^\s*(\d{1,4})([./-])(\d{1,2})\2(\d{1,4})\s*$")


def _numeric_body(value: Any) -> Optional[str]:
    """The digits-and-separators core of a value, or None if it is not one.

    Returning None for anything carrying a letter is what keeps a column
    of document references out of the number machinery.
    """
    text = _CURRENCY_CODE_RE.sub("", str(value if value is not None else ""))
    text = _NUM_NOISE_RE.sub("", text).strip()
    if not text or not _NUMERIC_SHAPE_RE.match(text):
        return None
    return text

DECIMAL_COMMA = "decimal_comma"     # 1.234,50 — German / most of Europe
DECIMAL_POINT = "decimal_point"     # 1,234.50 — English
PLAIN_NUMBER = "plain"              # 1234.5 / 1234 — no grouping in play
AMBIGUOUS = "ambiguous"

DAY_FIRST = "day_first"             # 31.07.2026
MONTH_FIRST = "month_first"         # 07/31/2026
ISO_DATE = "iso"                    # 2026-07-31


def _number_evidence(text: str) -> Optional[str]:
    """What one value proves about its column's convention, if anything."""
    body = _numeric_body(text)
    if body is None:
        return None
    dots, commas = body.count("."), body.count(",")

    # Both separators present: the LAST one is the decimal separator.
    # This is decisive and needs no column context.
    if dots and commas:
        return DECIMAL_COMMA if body.rfind(",") > body.rfind(".") else DECIMAL_POINT
    # More than one of the same separator can only be grouping.
    if dots > 1:
        return DECIMAL_COMMA
    if commas > 1:
        return DECIMAL_POINT
    # Exactly one separator: a group that is not three digits cannot be a
    # thousands separator, so the separator is decimal.
    for sep, convention in ((".", DECIMAL_POINT), (",", DECIMAL_COMMA)):
        if body.count(sep) == 1:
            tail = body.split(sep)[1]
            if len(tail) != 3:
                return convention
            return None                      # 1.234 / 1,234 — undecidable
    return PLAIN_NUMBER                       # no separator at all


def detect_number_convention(values: list) -> tuple[str, str]:
    """Decide a column's number convention. Returns ``(convention, why)``.

    ``PLAIN_NUMBER`` means every value parses the same either way (no
    grouping separators in play). ``AMBIGUOUS`` means the column contains
    only values like "1.234" that read as 1234 or as 1.234 depending on
    the convention, and nothing in the column settles it — the caller
    must ask rather than pick.
    """
    votes: dict[str, int] = {}
    for value in values:
        evidence = _number_evidence(value)
        if evidence in (DECIMAL_COMMA, DECIMAL_POINT):
            votes[evidence] = votes.get(evidence, 0) + 1
    if len(votes) > 1:
        top = max(votes, key=lambda k: votes[k])
        other = min(votes, key=lambda k: votes[k])
        return AMBIGUOUS, (
            f"the column mixes conventions: {votes[top]} value(s) look like "
            f"{top}, {votes[other]} like {other}")
    if votes:
        convention = next(iter(votes))
        return convention, f"{votes[convention]} value(s) settle it"

    undecided = [
        v for v in values
        if _number_evidence(v) is None and _numeric_body(v) is not None
    ]
    if undecided:
        return AMBIGUOUS, (
            f"values like {str(undecided[0])!r} read as two different numbers "
            "depending on the convention, and nothing in the column decides it")
    return PLAIN_NUMBER, "no grouping separators in this column"


def parse_number(text: Any, convention: str = PLAIN_NUMBER) -> Optional[float]:
    """Parse one value under a known convention. None if it is not a number."""
    if isinstance(text, bool):
        return None
    if isinstance(text, (int, float)):
        return float(text)
    body = _numeric_body(text)
    if body is None:
        return None
    if convention == DECIMAL_COMMA:
        body = body.replace(".", "").replace(",", ".")
    elif convention == DECIMAL_POINT:
        body = body.replace(",", "")
    else:
        # No convention decided: a lone separator with a 3-digit tail is
        # left alone only when it cannot change the value.
        if body.count(",") == 1 and "." not in body:
            body = body.replace(",", ".")
        else:
            body = body.replace(",", "")
    try:
        return float(body)
    except ValueError:
        return None


def _date_parts(text: Any) -> Optional[tuple[int, int, int, str]]:
    """``(a, b, c, separator)`` of a date-shaped value, else None."""
    match = _DATE_SEP_RE.match(str(text or ""))
    if match is None:
        return None
    a, sep, b, c = match.group(1), match.group(2), match.group(3), match.group(4)
    return int(a), int(b), int(c), sep


def detect_date_convention(values: list) -> tuple[str, str]:
    """Decide a column's date convention. Returns ``(convention, why)``."""
    saw_iso = saw_dot = False
    day_first_proof = month_first_proof = 0
    for value in values:
        parts = _date_parts(value)
        if parts is None:
            continue
        a, b, _c, sep = parts
        if sep == "-" and a > 31:
            saw_iso = True
            continue
        if sep == ".":
            saw_dot = True
        # A first field above 12 can only be a day; a second field above
        # 12 can only be a day in second position.
        if a > 12:
            day_first_proof += 1
        elif b > 12:
            month_first_proof += 1
    if saw_iso and not (day_first_proof or month_first_proof):
        return ISO_DATE, "ISO dates (YYYY-MM-DD)"
    if day_first_proof and month_first_proof:
        return AMBIGUOUS, (
            f"{day_first_proof} value(s) can only be day-first and "
            f"{month_first_proof} can only be month-first — the column mixes "
            "two conventions")
    if day_first_proof:
        return DAY_FIRST, f"{day_first_proof} value(s) have a day above 12"
    if month_first_proof:
        return MONTH_FIRST, f"{month_first_proof} value(s) have a month above 12"
    if saw_dot:
        # Dotted dates are day-first everywhere they are written.
        return DAY_FIRST, "dot-separated dates are written day-first"
    if saw_iso:
        return ISO_DATE, "ISO dates (YYYY-MM-DD)"
    return AMBIGUOUS, (
        "every value could be read either way (no day above 12 anywhere)")


def parse_date(text: Any, convention: str = ISO_DATE) -> Optional[str]:
    """Parse one date under a known convention, returned as ISO ``YYYY-MM-DD``."""
    import datetime as _dt

    if isinstance(text, _dt.datetime):
        return text.date().isoformat()
    if isinstance(text, _dt.date):
        return text.isoformat()
    parts = _date_parts(text)
    if parts is None:
        return None
    a, b, c, sep = parts
    if convention == ISO_DATE or (sep == "-" and a > 31):
        year, month, day = a, b, c
    elif convention == MONTH_FIRST:
        month, day, year = a, b, c
    elif convention == DAY_FIRST:
        day, month, year = a, b, c
    else:
        return None
    if year < 100:                     # two-digit year, as spreadsheets do
        year += 2000 if year < 70 else 1900
    try:
        return _dt.date(year, month, day).isoformat()
    except ValueError:
        return None


def profile_column(values: list, *, name: str = "") -> dict:
    """Describe one column: what it holds, under which convention.

    ``kind`` is ``number`` / ``date`` / ``text`` / ``empty``. A column is
    only called numeric or date-like when most of its non-empty values
    parse that way — one stray "n/a" must not turn a money column into
    text, and one date-shaped code must not turn a text column into
    dates. Values that do NOT parse are listed, because those are the
    rows a total would silently omit.
    """
    filled = [v for v in values if str(v or "").strip() != ""]
    if not filled:
        return {"name": name, "kind": "empty", "convention": "",
                "values": 0, "parsed": 0, "unparsed": []}

    date_shaped = sum(1 for v in filled if _date_parts(v) is not None)
    numeric_shaped = sum(
        1 for v in filled
        if (isinstance(v, (int, float)) and not isinstance(v, bool))
        or _numeric_body(v) is not None
    )

    if date_shaped >= max(1, int(0.8 * len(filled))):
        convention, why = detect_date_convention(filled)
        parsed = [(v, parse_date(v, convention)) for v in filled]
        return {
            "name": name, "kind": "date", "convention": convention,
            "convention_reason": why, "values": len(filled),
            "parsed": sum(1 for _, p in parsed if p is not None),
            "unparsed": [str(v) for v, p in parsed if p is None][:10],
        }
    if numeric_shaped >= max(1, int(0.8 * len(filled))):
        convention, why = detect_number_convention(filled)
        parsed = [(v, parse_number(v, convention)) for v in filled]
        return {
            "name": name, "kind": "number", "convention": convention,
            "convention_reason": why, "values": len(filled),
            "parsed": sum(1 for _, p in parsed if p is not None),
            "unparsed": [str(v) for v, p in parsed if p is None][:10],
        }
    return {"name": name, "kind": "text", "convention": "",
            "values": len(filled), "parsed": len(filled), "unparsed": []}


def profile_table(rows: list[list], *, header: bool = True) -> list[dict]:
    """Profile every column of a table (first row taken as the header)."""
    if not rows:
        return []
    width = max(len(r) for r in rows)
    names = (
        [str(rows[0][i]) if i < len(rows[0]) else f"col{i + 1}"
         for i in range(width)]
        if header else [column_letter(i + 1) for i in range(width)]
    )
    body = rows[1:] if header else rows
    out = []
    for index in range(width):
        column = [r[index] if index < len(r) else "" for r in body]
        out.append(profile_column(column, name=names[index]))
    return out


def column_notes(profiles: list[dict]) -> list[str]:
    """The findings from a profile that a reader must not miss."""
    notes: list[str] = []
    for entry in profiles:
        name = entry.get("name") or "?"
        if entry.get("convention") == AMBIGUOUS:
            notes.append(
                f"column '{name}': {entry.get('convention_reason', '')}. "
                "Ask which reading is meant — do not compute on it as it is.")
        # The two findings are independent: a column can both use a
        # decimal comma AND hold values that are not numbers, and each
        # one on its own is enough to make a naive total wrong.
        if entry.get("kind") == "number" and entry.get("convention") == DECIMAL_COMMA:
            notes.append(
                f"column '{name}': numbers use a decimal comma "
                "(1.234,50 = 1234.5). Convert before computing — read as-is "
                "they are off by a factor of a thousand.")
        elif entry.get("kind") == "date" and entry.get("convention") == DAY_FIRST:
            notes.append(
                f"column '{name}': dates are day-first (31.07.2026 = "
                "2026-07-31).")
        unparsed = entry.get("unparsed") or []
        if unparsed:
            missing = entry["values"] - entry["parsed"]
            notes.append(
                f"column '{name}': {missing} of {entry['values']} value(s) "
                f"are not {entry['kind']}s (e.g. {', '.join(unparsed[:3])}). "
                "A total over this column would leave them out.")
    return notes


# ---------------------------------------------------------------------------
# Spreadsheets — reading
# ---------------------------------------------------------------------------

# Encodings tried in order when reading a text table. utf-8 first (it
# fails loudly on non-utf-8 bytes rather than producing plausible
# nonsense), then the two Windows encodings a spreadsheet program on a
# German-language system writes. cp1252 before latin-1 because it is a
# superset for the printable range and decodes the Euro sign and typographic
# quotes that latin-1 turns into control characters.
_TEXT_ENCODINGS: tuple[str, ...] = ("utf-8-sig", "cp1252", "latin-1")


def decode_text(raw: bytes) -> tuple[str, str]:
    """Decode table bytes, returning ``(text, encoding_used)``.

    A CSV exported from a spreadsheet program on a German-language
    Windows system is cp1252, not utf-8. Decoding it as utf-8 with
    ``errors="replace"`` does not fail — it turns 'Müller' into
    'M�ller' and 'Straße' into 'Stra�e', and the agent then
    writes those broken names into letters. Every name and street in a
    German data set carries at least one of those bytes, so guessing
    wrong here corrupts essentially the whole file, silently.
    """
    for encoding in _TEXT_ENCODINGS:
        try:
            return raw.decode(encoding), encoding
        except UnicodeDecodeError:
            continue
    # latin-1 maps every byte, so this is unreachable in practice; keep a
    # defined outcome rather than an exception if that ever changes.
    return raw.decode("utf-8", errors="replace"), "utf-8 (with damage)"


def sniff_delimiter(text: str, suffix: str = "") -> str:
    """Pick the separator that actually splits this file into columns.

    Counting occurrences is the obvious approach and the wrong one for
    German data: a decimal comma appears in every amount, so a
    semicolon-separated file looks comma-separated and collapses into one
    column. What distinguishes the real separator is that it produces the
    SAME number of fields in every row — the others land in the middle of
    values and vary from line to line.
    """
    if suffix.lower() == ".tsv":
        return "\t"
    import csv as _csv
    import io as _io

    sample = [ln for ln in text.splitlines()[:60] if ln.strip()][:40]
    if not sample:
        return ","
    best, best_score = ",", (-1.0, 0)
    for candidate in (";", ",", "\t", "|"):
        try:
            rows = list(_csv.reader(_io.StringIO("\n".join(sample)),
                                    delimiter=candidate))
        except _csv.Error:
            continue
        widths = [len(r) for r in rows if r]
        if not widths or max(widths) < 2:
            continue
        common = max(set(widths), key=widths.count)
        consistency = widths.count(common) / len(widths)
        score = (consistency, common)
        if score > best_score:
            best, best_score = candidate, score
    return best


def _read_csv(
    path: Path, *, max_rows: int, max_cols: int, start_row: int
) -> dict:
    import io

    # A .csv that is not text at all. Decoding a workbook or a PDF with
    # latin-1 never fails — every byte maps to some character — so this
    # would otherwise come back as a one-column table of mojibake that
    # has a row count and a header and means nothing.
    found = container_kind(path)
    if found in ("ole2", "zip", "xlsx", "docx", "ods", "odt", "odf", "pdf"):
        raise OfficeError(
            f"{path.name} is named like a text table but holds "
            f"{_CONTAINER_LABELS.get(found, found)}. Decoding it as text "
            "would return a table of noise rather than its contents — give "
            "it the right suffix and read it again, or convert it.")

    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise OfficeError(f"could not read {path.name}: {exc}") from exc
    text, encoding = decode_text(raw)

    delimiter = sniff_delimiter(text, path.suffix)
    all_rows = list(csv.reader(io.StringIO(text, newline=""),
                              delimiter=delimiter))

    total_rows = len(all_rows)
    total_cols = max((len(r) for r in all_rows), default=0)
    begin = max(0, start_row - 1)
    window = [r[:max_cols] for r in all_rows[begin:begin + max_rows]]
    notes = []
    if delimiter != ",":
        notes.append(f"delimiter detected: {delimiter!r}")
    if encoding != "utf-8-sig":
        notes.append(
            f"decoded as {encoding} (not utf-8) — this file came from a "
            "spreadsheet program on Windows. Anything you write back should "
            "keep the accented characters intact.")
    if total_rows > begin + len(window):
        notes.append(
            f"showing rows {begin + 1}-{begin + len(window)} of {total_rows} "
            f"— pass start_row to page further")
    if total_cols > max_cols:
        notes.append(f"showing {max_cols} of {total_cols} columns")
    # Rows that do not have the header's width. An unquoted separator
    # inside a value — "Digitalwaage 0,1 mg" in a comma-separated file —
    # shifts every field after it by one, so a personnel number column
    # quietly starts holding amounts. Nothing about the resulting table
    # looks wrong, which is why it has to be said.
    width = len(all_rows[0]) if all_rows else 0
    ragged = [i for i, r in enumerate(all_rows[1:], start=2)
              if len(r) != width and any(str(c).strip() for c in r)]
    if ragged and width:
        shown = ", ".join(str(i) for i in ragged[:8])
        notes.append(
            f"{len(ragged)} row(s) do not have the header's {width} column(s) "
            f"— rows {shown}{' …' if len(ragged) > 8 else ''}. A separator "
            "inside an unquoted value shifts every field after it, so the "
            "columns of those rows do not mean what their header says.")

    profiles = profile_table(all_rows[:2000])
    notes.extend(column_notes(profiles))
    return {
        "path": str(path),
        "kind": "csv",
        "sheets": [],
        "sheet": "",
        "rows": total_rows,
        "columns": total_cols,
        "grid": render_grid(window, first_row=begin + 1),
        "column_profile": profiles,
        "notes": notes,
    }


def _fragile_features(workbook: Any) -> list[str]:
    """Workbook contents that a rewrite puts at risk, named per sheet.

    ``openpyxl`` reconstructs the file on save. Recent versions carry
    charts, images, merges and number formats through a round-trip, but
    the library does not model every part of the format, and a workbook
    authored by Excel can hold constructs it has never seen. Naming what
    is present lets the caller decide, and lets the user hear about it
    before the file is rewritten rather than after.
    """
    found: list[str] = []
    for sheet in getattr(workbook, "worksheets", []) or []:
        title = getattr(sheet, "title", "?")
        for attr, label in (
            ("_charts", "chart"),
            ("_images", "image"),
            ("_pivots", "pivot table"),
            ("_tables", "table"),
        ):
            items = getattr(sheet, attr, None) or []
            try:
                count = len(items)
            except TypeError:
                count = 0
            if count:
                found.append(
                    f"{title}: {count} {label}{'s' if count != 1 else ''}")
    if getattr(workbook, "vba_archive", None) is not None:
        found.append("VBA macros")
    return found


# ---------------------------------------------------------------------------
# OpenDocument spreadsheets (.ods)
# ---------------------------------------------------------------------------
#
# Three properties of the format decide whether a reader reports the
# file or something adjacent to it:
#
# * Cells repeat. LibreOffice writes one cell with
#   ``table:number-columns-repeated="1013"`` rather than a thousand empty
#   ones, and one row with ``number-rows-repeated="1048571"`` for the
#   rest of the sheet. A reader that ignores the attribute shifts every
#   value after a repeated cell into the wrong column.
# * A cell carries both a stored value and a formatted display text:
#   ``office:value="1234.5"`` next to ``1.234,50 €``. The stored value is
#   what the file means; taking the display text would put the whole
#   decimal-comma problem back into a file that had already settled it.
# * A merged range leaves ``table:covered-table-cell`` elements behind.
#   They hold nothing but they occupy their position, so dropping them
#   moves everything to their right one column left.

MAX_ODS_ROWS = 100_000
# Nothing sensible has more columns than this; the number exists to stop
# a repeat count of a million from being materialised.
_ODS_MAX_COLS = 4096


def _ods_blank(value: Any) -> bool:
    return value is None or (isinstance(value, str) and not value.strip())


def _odf_load(p: Path) -> Any:
    """Open an OpenDocument file, or say what the file really is."""
    _require("opendocument")
    from odf.opendocument import load
    try:
        return load(str(p))
    except Exception as exc:
        raise OfficeError(
            f"could not open {p.name} as an OpenDocument file: "
            f"{type(exc).__name__}: {exc}."
            + (_mismatch_note(p, ("ods", "odt", "odf")) or
               " The file may be damaged or incomplete.")
        ) from exc


def _ods_cell_value(cell: Any) -> tuple[Any, bool]:
    """One cell as ``(value, is_uncached_formula)``.

    Typed values are returned in their own type — a float stays a float,
    a date comes back as a ``date`` — so the column profiler sees what
    the file stores rather than how a locale renders it.
    """
    import datetime as _dt

    from odf.namespaces import OFFICENS, TABLENS
    from odf import teletype

    attrs = dict(getattr(cell, "attributes", None) or {})
    value_type = attrs.get((OFFICENS, "value-type"))
    formula = attrs.get((TABLENS, "formula"))

    def _text() -> str:
        try:
            return teletype.extractText(cell)
        except Exception:
            return ""

    if value_type in ("float", "percentage", "currency"):
        try:
            return float(attrs.get((OFFICENS, "value"), "")), False
        except (TypeError, ValueError):
            return _text(), False
    if value_type == "date":
        raw = str(attrs.get((OFFICENS, "date-value"), ""))
        # A day is returned as a ``date``, not as a midnight ``datetime``:
        # rendered, the second one carries a 00:00:00 the file never
        # stated, and the column profiler then reads the whole column as
        # text rather than as dates.
        try:
            return (_dt.date.fromisoformat(raw) if len(raw) == 10
                    else _dt.datetime.fromisoformat(raw)), False
        except ValueError:
            return _text(), False
    if value_type == "boolean":
        raw = str(attrs.get((OFFICENS, "boolean-value"), "")).lower()
        if raw in ("true", "false"):
            return raw == "true", False
        return _text(), False
    if value_type == "string":
        stored = attrs.get((OFFICENS, "string-value"))
        return (stored if stored is not None else _text()), False
    if value_type == "time":
        # Stored as an ISO 8601 duration (PT10H30M00S), which is not what
        # anyone wrote in the cell. The display text is.
        return _text(), False

    text = _text()
    if formula and not text:
        # A formula with no stored result. Same case as a workbook
        # written by a library: the number does not exist in the file
        # yet, and showing an empty cell would report it as one.
        cleaned = str(formula)
        if cleaned.startswith("of:"):
            cleaned = cleaned[3:]
        return cleaned, True
    return text, False


def _ods_row_cells(row: Any) -> list:
    """The cells of one row, repeats expanded and trailing blanks cut."""
    from odf.namespaces import TABLENS

    cells: list[Any] = []
    for node in getattr(row, "childNodes", None) or []:
        qname = getattr(node, "qname", None)
        if not qname or qname[0] != TABLENS:
            continue
        if qname[1] not in ("table-cell", "covered-table-cell"):
            continue
        attrs = dict(getattr(node, "attributes", None) or {})
        try:
            repeat = int(attrs.get((TABLENS, "number-columns-repeated"), 1))
        except (TypeError, ValueError):
            repeat = 1
        repeat = max(1, min(repeat, _ODS_MAX_COLS - len(cells)))
        value = (None if qname[1] == "covered-table-cell"
                 else _ods_cell_value(node)[0])
        cells.extend([value] * repeat)
        if len(cells) >= _ODS_MAX_COLS:
            break
    while cells and _ods_blank(cells[-1]):
        cells.pop()
    return cells


def _ods_tables(document: Any) -> list:
    """Every sheet of an .ods, in document order."""
    from odf.namespaces import TABLENS

    body = getattr(document, "spreadsheet", None)
    if body is None:
        return []
    return [node for node in (getattr(body, "childNodes", None) or [])
            if getattr(node, "qname", None) == (TABLENS, "table")]


def _ods_sheet_names(document: Any) -> list[str]:
    from odf.namespaces import TABLENS
    names = []
    for index, table in enumerate(_ods_tables(document), start=1):
        attrs = dict(getattr(table, "attributes", None) or {})
        names.append(str(attrs.get((TABLENS, "name"), f"Tabelle{index}")))
    return names


def _ods_sheet(
    p: Path, sheet: Optional[str] = None
) -> tuple[list[str], str, list, int, int, int]:
    """``(sheet names, chosen, row entries, rows, columns, hidden rows)``.

    Row entries are ``(cells, repeat)`` pairs — unexpanded, because the
    repeat count at the end of a sheet routinely runs to a million and
    materialising it would turn a 4 kB file into a gigabyte of ``None``.
    """
    # Load first: it is the call that reports a missing odfpy as
    # something the reader can act on. Importing the namespaces above it
    # would raise a bare ImportError out of a document read instead.
    document = _odf_load(p)

    from odf.namespaces import TABLENS

    mimetype = str(getattr(document, "mimetype", "") or "")
    if mimetype and not mimetype.endswith(
            ("spreadsheet", "spreadsheet-template")):
        raise OfficeError(
            f"{p.name} is named like a spreadsheet but its OpenDocument "
            f"type is '{mimetype}'. Reading it as a table would invent a "
            "structure the file does not have.")

    names = _ods_sheet_names(document)
    tables = _ods_tables(document)
    if not tables:
        raise OfficeError(f"{p.name} holds no sheet")
    if sheet:
        if sheet not in names:
            raise OfficeError(
                f"no sheet named {sheet!r} — available: "
                f"{', '.join(names) or '(none)'}")
        table = tables[names.index(sheet)]
        chosen = sheet
    else:
        table, chosen = tables[0], names[0]

    entries: list[tuple[list, int]] = []
    concealed: list[bool] = []
    for node in getattr(table, "childNodes", None) or []:
        qname = getattr(node, "qname", None)
        if not qname or qname[0] != TABLENS:
            continue
        if qname[1] == "table-header-rows":
            # A repeated header is wrapped in its own element; its rows
            # are ordinary rows and belong at the top of the sheet.
            for inner in getattr(node, "childNodes", None) or []:
                if getattr(inner, "qname", None) == (TABLENS, "table-row"):
                    entries.append((_ods_row_cells(inner), 1))
                    concealed.append(False)
            continue
        if qname[1] != "table-row":
            continue
        attrs = dict(getattr(node, "attributes", None) or {})
        try:
            repeat = max(
                1, int(attrs.get((TABLENS, "number-rows-repeated"), 1)))
        except (TypeError, ValueError):
            repeat = 1
        entries.append((_ods_row_cells(node), repeat))
        concealed.append(
            str(attrs.get((TABLENS, "visibility"), "visible")) != "visible")

    while entries and not entries[-1][0]:
        entries.pop()
        concealed.pop()
    # Counted after the trim: the filler at the end of a sheet can carry
    # a visibility attribute too, and counting it would report a million
    # hidden rows for a file that has none.
    hidden = sum(repeat for (_cells, repeat), flag
                 in zip(entries, concealed) if flag)
    total_rows = sum(repeat for _cells, repeat in entries)
    total_cols = max((len(c) for c, _r in entries), default=0)
    return names, chosen, entries, total_rows, total_cols, hidden


def _ods_window(entries: list, start_row: int, max_rows: int) -> list[list]:
    """Materialise rows ``[start_row, start_row + max_rows)``, 1-based."""
    out: list[list] = []
    position = 1
    for cells, repeat in entries:
        if len(out) >= max_rows:
            break
        end = position + repeat
        if end > start_row:
            first = max(position, start_row)
            take = min(end, start_row + max_rows) - first
            out.extend([list(cells) for _ in range(max(0, take))])
        position = end
    return out


def _read_ods(
    p: Path, *, max_rows: int, max_cols: int, start_row: int,
    sheet: Optional[str] = None,
) -> dict:
    names, chosen, entries, total_rows, total_cols, hidden = _ods_sheet(
        p, sheet)
    window = [row[:max_cols] for row in
              _ods_window(entries, start_row, max_rows)]

    notes: list[str] = []
    uncached = sum(
        1 for row in window for value in row
        if isinstance(value, str) and value.startswith("=")
    )
    if uncached:
        notes.append(
            f"{uncached} formula cell(s) carry no stored result and are "
            "shown as their formula text. Nothing has evaluated this file, "
            "so those numbers do not exist in it yet — do not report them "
            "as computed.")
    if len(names) > 1:
        notes.append(f"sheets in this workbook: {', '.join(names)}")
    if total_rows > start_row + len(window) - 1:
        notes.append(
            f"showing rows {start_row}-{start_row + len(window) - 1} of "
            f"{total_rows} — pass start_row to page further")
    if total_cols > max_cols:
        notes.append(f"showing {max_cols} of {total_cols} columns")
    if hidden:
        notes.append(
            f"{hidden} row(s) of this sheet are hidden or filtered out. "
            "They are included in the grid above, so the totals here can "
            "differ from what the file shows on screen.")

    profiles = profile_table(window, header=(start_row == 1))
    notes.extend(column_notes(profiles))
    return {
        "path": str(p),
        # The content is a spreadsheet whatever the container is; the
        # format is carried separately so a caller can tell that writing
        # is not on the table.
        "kind": "spreadsheet",
        "format": "ods",
        "writable": False,
        "sheets": names,
        "sheet": chosen,
        "rows": total_rows,
        "columns": total_cols,
        "grid": render_grid(window, first_row=start_row),
        "column_profile": profiles,
        "notes": notes,
    }


def read_sheet(
    path: Any,
    *,
    sheet: Optional[str] = None,
    max_rows: int = DEFAULT_MAX_ROWS,
    max_cols: int = DEFAULT_MAX_COLS,
    start_row: int = 1,
) -> dict:
    """Read a window of a spreadsheet or CSV as an addressable grid.

    Formula cells are reported by their computed value where the file
    carries one. Where it does not — the usual case for a workbook
    written by a library and never opened in a spreadsheet program — the
    formula text is shown instead and a note says so, because silently
    printing an empty column would misreport the file.
    """
    p = _resolve(path)
    _refuse_unreadable_format(p)
    kind = document_kind(p)
    max_rows = max(1, min(int(max_rows), 2000))
    max_cols = max(1, min(int(max_cols), 200))
    start_row = max(1, int(start_row))

    if kind == "csv":
        return _read_csv(
            p, max_rows=max_rows, max_cols=max_cols, start_row=start_row)
    if kind == "opendocument_sheet":
        return _read_ods(p, max_rows=max_rows, max_cols=max_cols,
                         start_row=start_row, sheet=sheet)
    if kind != "spreadsheet":
        raise OfficeError(
            f"{p.name} is not a spreadsheet (suffix {p.suffix!r})")

    openpyxl = _require("spreadsheet")
    try:
        wb = openpyxl.load_workbook(p, data_only=False, read_only=False)
    except Exception as exc:
        raise OfficeError(f"could not open workbook: {exc}") from exc

    try:
        names = list(wb.sheetnames)
        if sheet:
            if sheet not in names:
                raise OfficeError(
                    f"no sheet named {sheet!r} — available: "
                    f"{', '.join(names) or '(none)'}")
            ws = wb[sheet]
        else:
            ws = wb.active

        total_rows = int(ws.max_row or 0)
        total_cols = int(ws.max_column or 0)
        window: list[list[Any]] = []
        has_formula = False
        for row in ws.iter_rows(
            min_row=start_row,
            max_row=min(total_rows, start_row + max_rows - 1),
            max_col=min(total_cols, max_cols),
            values_only=True,
        ):
            window.append(list(row))
            for value in row:
                if isinstance(value, str) and value.startswith("="):
                    has_formula = True

        notes: list[str] = []
        # Second pass only when it can change what is shown: swap in
        # cached values for the formula cells, and say which ones have
        # none rather than rendering them as blanks.
        if has_formula:
            uncached = 0
            try:
                wb_values = openpyxl.load_workbook(
                    p, data_only=True, read_only=True)
                ws_values = (
                    wb_values[ws.title] if ws.title in wb_values.sheetnames
                    else wb_values.active
                )
                for r_idx, row in enumerate(ws_values.iter_rows(
                    min_row=start_row,
                    max_row=min(total_rows, start_row + max_rows - 1),
                    max_col=min(total_cols, max_cols),
                    values_only=True,
                )):
                    if r_idx >= len(window):
                        break
                    for c_idx, cached in enumerate(row):
                        if c_idx >= len(window[r_idx]):
                            break
                        current = window[r_idx][c_idx]
                        if isinstance(current, str) and current.startswith("="):
                            if cached is None:
                                uncached += 1
                            else:
                                window[r_idx][c_idx] = cached
                wb_values.close()
            except Exception:
                notes.append(
                    "formula cells are shown as formulas — the cached "
                    "values could not be read")
                uncached = 0
            if uncached:
                notes.append(
                    f"{uncached} formula cell(s) carry no cached value and "
                    "are shown as their formula text. A spreadsheet program "
                    "has never evaluated this file, so those numbers do not "
                    "exist in it yet — do not report them as computed.")

        if len(names) > 1:
            notes.append(f"sheets in this workbook: {', '.join(names)}")
        if total_rows > start_row + len(window) - 1:
            notes.append(
                f"showing rows {start_row}-{start_row + len(window) - 1} of "
                f"{total_rows} — pass start_row to page further")
        if total_cols > max_cols:
            notes.append(f"showing {max_cols} of {total_cols} columns")
        fragile = _fragile_features(wb)
        if fragile:
            notes.append("fragile content present: " + "; ".join(fragile))

        # Profile the window that was read, treating its first row as the
        # header only when it really is one (paging into the middle of a
        # sheet must not promote a data row to a column name).
        profiles = profile_table(window, header=(start_row == 1))
        notes.extend(column_notes(profiles))

        return {
            "path": str(p),
            "kind": "spreadsheet",
            "sheets": names,
            "sheet": ws.title,
            "rows": total_rows,
            "columns": total_cols,
            "grid": render_grid(window, first_row=start_row),
            "column_profile": profiles,
            "notes": notes,
        }
    finally:
        try:
            wb.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Spreadsheets — writing
# ---------------------------------------------------------------------------

def plan_sheet_edits(
    path: Any,
    *,
    edits: Optional[list[dict]] = None,
    append_rows: Optional[list[list[Any]]] = None,
    sheet: Optional[str] = None,
) -> dict:
    """Resolve what :func:`edit_sheet` would change, without writing.

    Returns the per-cell before/after list, the fragile features the
    rewrite would put at risk, and a one-line summary. The tool layer
    shows this to the user as the approval preview: a binary diff is
    unreadable, whereas ``Daten!B7: 100 -> 140`` is exactly the question
    being asked.
    """
    p = _resolve(path)
    _refuse_unwritable_format(p, what="editing it here")
    if document_kind(p) != "spreadsheet":
        raise OfficeError(
            f"{p.name} is not a spreadsheet (suffix {p.suffix!r})")
    openpyxl = _require("spreadsheet")
    try:
        wb = openpyxl.load_workbook(p, data_only=False)
    except Exception as exc:
        raise OfficeError(f"could not open workbook: {exc}") from exc
    try:
        if sheet:
            if sheet not in wb.sheetnames:
                raise OfficeError(
                    f"no sheet named {sheet!r} — available: "
                    f"{', '.join(wb.sheetnames)}")
            ws = wb[sheet]
        else:
            ws = wb.active

        changes: list[dict] = []
        for edit in edits or []:
            if not isinstance(edit, dict):
                raise OfficeError(f"each edit must be an object, got {edit!r}")
            ref = str(edit.get("cell", "")).strip()
            row, col = parse_cell(ref)
            old = ws.cell(row=row, column=col).value
            changes.append({
                "cell": f"{column_letter(col)}{row}",
                "old": _fmt(old),
                "new": _fmt(edit.get("value")),
            })

        appended = len(append_rows or [])
        return {
            "path": str(p),
            "sheet": ws.title,
            "changes": changes,
            "appended_rows": appended,
            "first_appended_row": int(ws.max_row or 0) + 1 if appended else 0,
            "fragile": _fragile_features(wb),
        }
    finally:
        try:
            wb.close()
        except Exception:
            pass


def _verify_cells(path: Path, sheet_name: str, applied: list[dict]) -> list[str]:
    """Re-read the saved file and report cells that do not hold the new value.

    The form filler reads its result back before reporting success; a
    spreadsheet edit reported from the in-memory workbook would be
    claiming an outcome from the intent. Never raises — a verification
    that cannot run is reported as unverified, not as a failure.
    """
    if not applied:
        return []
    try:
        import openpyxl as _op
        wb = _op.load_workbook(path, data_only=False)
    except Exception:
        return []
    try:
        ws = wb[sheet_name] if sheet_name in wb.sheetnames else wb.active
        out: list[str] = []
        for change in applied:
            ref = str(change.get("cell", ""))
            try:
                row, col = parse_cell(ref)
            except OfficeError:
                continue
            if _fmt(ws.cell(row=row, column=col).value) != str(change.get("new", "")):
                out.append(ref)
        return out
    except Exception:
        return []
    finally:
        try:
            wb.close()
        except Exception:
            pass


def _sheet_header(ws: Any) -> list[str]:
    """The first row as column names (blank cells become empty strings)."""
    for row in ws.iter_rows(min_row=1, max_row=1, values_only=True):
        return [str(v).strip() if v is not None else "" for v in row]
    return []


def _locate_rows(ws: Any, key_column: str) -> tuple[int, dict, list[str]]:
    """``(column index, key -> row, duplicate keys)`` for a key column.

    Duplicates are collected rather than resolved: updating "the row of
    Müller" when there are two of them is not a decision a tool may make
    silently.
    """
    header = _sheet_header(ws)
    wanted = str(key_column).strip().lower()
    index = next((i for i, name in enumerate(header)
                  if name.strip().lower() == wanted), -1)
    if index < 0:
        raise OfficeError(
            f"no column {key_column!r} in this sheet. Columns: "
            + (", ".join(h for h in header if h) or "(none)"))
    by_key: dict[str, int] = {}
    duplicates: list[str] = []
    for number, row in enumerate(
            ws.iter_rows(min_row=2, values_only=True), start=2):
        raw = row[index] if index < len(row) else None
        token = str(raw).strip() if raw is not None else ""
        if not token:
            continue
        if token in by_key:
            if token not in duplicates:
                duplicates.append(token)
        else:
            by_key[token] = number
    return index, by_key, duplicates


def edit_sheet(
    path: Any,
    *,
    edits: Optional[list[dict]] = None,
    append_rows: Optional[list[list[Any]]] = None,
    key_column: Optional[str] = None,
    updates: Optional[list[dict]] = None,
    append_records: Optional[list[dict]] = None,
    sheet: Optional[str] = None,
) -> dict:
    """Set cells and append rows in an existing workbook, in place.

    Formulas elsewhere in the file are preserved (the workbook is loaded
    with ``data_only=False``; loading with cached values instead would
    replace every formula in the file with its last computed number).
    Macro-enabled workbooks keep their VBA project.

    A cell value beginning with ``=`` is stored as a formula, which
    means the file then carries no cached result for it until a
    spreadsheet program opens and recalculates the file. That is stated
    in the returned notes rather than left for the reader to discover.
    """
    p = _resolve(path)
    _refuse_unwritable_format(p, what="editing it here")
    if document_kind(p) != "spreadsheet":
        raise OfficeError(
            f"{p.name} is not a spreadsheet (suffix {p.suffix!r})")
    if not any((edits, append_rows, updates, append_records)):
        raise OfficeError(
            "nothing to do: pass edits/append_rows (by cell) or "
            "updates/append_records (by key and column name)")
    if updates and not key_column:
        raise OfficeError(
            "updates need key_column — the column whose value identifies a "
            "row (a document number, a case reference)")

    openpyxl = _require("spreadsheet")
    keep_vba = p.suffix.lower() in (".xlsm", ".xltm")
    try:
        wb = openpyxl.load_workbook(p, data_only=False, keep_vba=keep_vba)
    except Exception as exc:
        raise OfficeError(f"could not open workbook: {exc}") from exc

    try:
        if sheet:
            if sheet not in wb.sheetnames:
                raise OfficeError(
                    f"no sheet named {sheet!r} — available: "
                    f"{', '.join(wb.sheetnames)}")
            ws = wb[sheet]
        else:
            ws = wb.active

        applied: list[dict] = []
        wrote_formula = False

        # Record-oriented changes. Everything is validated BEFORE the first
        # cell is written: an unknown key, a duplicate key or an unknown
        # column stops the whole call. A half-applied batch on someone's
        # records is worse than a refused one, because nothing about the
        # file afterwards says which half went through.
        if updates or append_records:
            header = _sheet_header(ws)
            named = {h.strip().lower(): i for i, h in enumerate(header) if h}
            if not named:
                raise OfficeError(
                    "this sheet has no header row, so columns cannot be "
                    "addressed by name — use edits with cell references")
            wanted_columns: set[str] = set()
            for entry in list(updates or []) + list(append_records or []):
                if not isinstance(entry, dict):
                    raise OfficeError(f"expected an object, got {entry!r}")
                values = entry.get("set") if updates and "set" in entry else entry
                if not isinstance(values, dict):
                    raise OfficeError(
                        f"each update needs a 'set' object of column -> value, "
                        f"got {entry!r}")
                wanted_columns.update(str(c) for c in values
                                      if str(c).strip().lower() != "key")
            unknown = sorted(c for c in wanted_columns
                             if c.strip().lower() not in named)
            if unknown:
                raise OfficeError(
                    f"no column(s) {', '.join(unknown)} in this sheet. "
                    "Columns: " + ", ".join(h for h in header if h))

        if updates:
            _key_index, by_key, duplicates = _locate_rows(ws, key_column)
            header = _sheet_header(ws)
            named = {h.strip().lower(): i for i, h in enumerate(header) if h}
            wanted_keys = [str(u.get("key", "")).strip() for u in updates]
            missing = sorted({k for k in wanted_keys if k and k not in by_key})
            blank = [i for i, k in enumerate(wanted_keys, start=1) if not k]
            ambiguous = sorted(set(wanted_keys) & set(duplicates))
            if blank:
                raise OfficeError(
                    f"update {blank[0]} has no key — every update needs the "
                    f"value in '{key_column}' that identifies its row")
            if missing:
                raise OfficeError(
                    f"no row with {key_column} = "
                    + ", ".join(repr(k) for k in missing)
                    + f". The sheet holds {len(by_key)} distinct key(s).")
            if ambiguous:
                raise OfficeError(
                    f"{key_column} = " + ", ".join(repr(k) for k in ambiguous)
                    + " appears more than once — which row is meant is not "
                    "decidable, so nothing was changed.")
            for update in updates:
                row_number = by_key[str(update.get("key", "")).strip()]
                for column, value in (update.get("set") or {}).items():
                    col_index = named[str(column).strip().lower()] + 1
                    cell = ws.cell(row=row_number, column=col_index)
                    old_value = cell.value
                    cell.value = value
                    if isinstance(value, str) and value.startswith("="):
                        wrote_formula = True
                    applied.append({
                        "cell": f"{column_letter(col_index)}{row_number}",
                        "key": str(update.get("key", "")).strip(),
                        "column": str(column),
                        "old": _fmt(old_value), "new": _fmt(value),
                    })

        for edit in edits or []:
            if not isinstance(edit, dict):
                raise OfficeError(f"each edit must be an object, got {edit!r}")
            row, col = parse_cell(str(edit.get("cell", "")).strip())
            value = edit.get("value")
            cell = ws.cell(row=row, column=col)
            old = cell.value
            cell.value = value
            if isinstance(value, str) and value.startswith("="):
                wrote_formula = True
            applied.append({
                "cell": f"{column_letter(col)}{row}",
                "old": _fmt(old),
                "new": _fmt(value),
            })

        first_appended = 0
        appended = 0

        # Records are placed by column NAME. A positional list silently
        # lands in the wrong columns as soon as the sheet's order differs
        # from the caller's assumption, and the result reads perfectly.
        if append_records:
            header = _sheet_header(ws)
            named = {h.strip().lower(): i for i, h in enumerate(header) if h}
            for record in append_records:
                if not appended:
                    first_appended = int(ws.max_row or 0) + 1
                target_row = int(ws.max_row or 0) + 1
                for column, value in record.items():
                    col_index = named[str(column).strip().lower()] + 1
                    ws.cell(row=target_row, column=col_index).value = value
                    if isinstance(value, str) and value.startswith("="):
                        wrote_formula = True
                appended += 1

        for new_row in append_rows or []:
            if not isinstance(new_row, (list, tuple)):
                raise OfficeError(
                    f"each appended row must be a list, got {new_row!r}")
            if not appended:
                first_appended = int(ws.max_row or 0) + 1
            ws.append(list(new_row))
            appended += 1

        notes: list[str] = []
        fragile = _fragile_features(wb)
        if fragile:
            notes.append(
                "this workbook carries content a rewrite can affect ("
                + "; ".join(fragile)
                + ") — the file was rebuilt on save, so verify those parts "
                "survived before handing it over")
        if wrote_formula:
            notes.append(
                "a formula was written: the file now stores the formula but "
                "no computed result. Anything reading cached values will see "
                "an empty cell until a spreadsheet program recalculates it.")

        try:
            wb.save(p)
        except Exception as exc:
            raise OfficeError(f"could not save workbook: {exc}") from exc

        sheet_name = ws.title
        unverified = _verify_cells(p, sheet_name, applied)
        if unverified:
            notes.append(
                "these cells do not hold the new value in the saved file: "
                + ", ".join(unverified)
                + ". Do not report this edit as done.")

        return {
            "path": str(p),
            "sheet": sheet_name,
            "applied": applied,
            "appended_rows": appended,
            "first_appended_row": first_appended,
            "verified": not unverified,
            "notes": notes,
        }
    finally:
        try:
            wb.close()
        except Exception:
            pass


def create_sheet(
    path: Any,
    rows: list[list[Any]],
    *,
    sheet_name: str = "Sheet1",
    overwrite: bool = False,
) -> dict:
    """Write *rows* to a new workbook (or CSV, by suffix)."""
    p = _resolve(path, must_exist=False)
    _refuse_unwritable_format(p, what="writing it")
    if p.exists() and not overwrite:
        raise OfficeError(
            f"{p.name} already exists — pass overwrite=True to replace it, "
            "or use edit_sheet to change it in place")
    if not isinstance(rows, list) or not all(
            isinstance(r, (list, tuple)) for r in rows):
        raise OfficeError("rows must be a list of lists")

    kind = document_kind(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    if kind == "csv":
        delimiter = "\t" if p.suffix.lower() == ".tsv" else ","
        # utf-8 WITH a BOM: a spreadsheet program opening a plain utf-8 CSV
        # by double-click assumes the system code page and renders every
        # umlaut as mojibake. The BOM is what makes it detect utf-8. It is
        # invisible to every other reader that matters here (Python's csv
        # strips it via utf-8-sig, as does the reader above).
        with p.open("w", encoding="utf-8-sig", newline="") as fh:
            csv.writer(fh, delimiter=delimiter).writerows(
                [list(r) for r in rows])
    elif kind == "spreadsheet":
        openpyxl = _require("spreadsheet")
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = str(sheet_name or "Sheet1")[:31]
        for row in rows:
            ws.append(list(row))
        wb.save(p)
        wb.close()
    else:
        raise OfficeError(
            f"cannot create {p.suffix!r} as a spreadsheet — use .xlsx or .csv")
    return {
        "path": str(p),
        "kind": kind,
        "rows": len(rows),
        "columns": max((len(r) for r in rows), default=0),
    }


# ---------------------------------------------------------------------------
# Reconciling two tables
# ---------------------------------------------------------------------------

def _workbook_sheets(path: Any) -> list[str]:
    """Sheet names of a workbook, or an empty list for anything else."""
    p = Path(str(path))
    if document_kind(p) == "opendocument_sheet":
        try:
            return _ods_sheet_names(_odf_load(p))
        except Exception:
            return []
    if document_kind(p) != "spreadsheet":
        return []
    try:
        openpyxl = _require("spreadsheet")
        wb = openpyxl.load_workbook(p, read_only=True)
        try:
            return list(wb.sheetnames)
        finally:
            wb.close()
    except Exception:
        return []


def _table_rows(path: Any, sheet: Optional[str] = None) -> list[list]:
    """Every row of a table file, as lists (no window, no cap on width)."""
    p = _resolve(path)
    _refuse_unreadable_format(p)
    kind = document_kind(p)
    if kind == "csv":
        raw = p.read_bytes()
        text, _enc = decode_text(raw)
        import io
        return list(csv.reader(io.StringIO(text, newline=""),
                               delimiter=sniff_delimiter(text, p.suffix)))
    if kind == "opendocument_sheet":
        _names, _chosen, entries, total, _cols, _hidden = _ods_sheet(p, sheet)
        # Unlike the read path there is no window to bound the work here,
        # so the cap has to be a refusal: silently comparing the first
        # 100 000 rows of a longer sheet would report a reconciliation
        # that never saw the rest of the file.
        if total > MAX_ODS_ROWS:
            raise OfficeError(
                f"{p.name} has {total} rows, above the {MAX_ODS_ROWS} this "
                "reader will materialise at once. Narrow the file, or "
                "convert it to .xlsx and work on that.")
        return _ods_window(entries, 1, total)
    if kind != "spreadsheet":
        raise OfficeError(
            f"{p.name} is not a table (need .xlsx, .ods, .csv or .tsv)")
    openpyxl = _require("spreadsheet")
    try:
        wb = openpyxl.load_workbook(p, data_only=True, read_only=True)
    except Exception as exc:
        raise OfficeError(f"could not open workbook: {exc}") from exc
    try:
        if sheet:
            if sheet not in wb.sheetnames:
                raise OfficeError(
                    f"no sheet named {sheet!r} in {p.name} — available: "
                    f"{', '.join(wb.sheetnames)}")
            ws = wb[sheet]
        else:
            ws = wb.active
        return [list(r) for r in ws.iter_rows(values_only=True)]
    finally:
        try:
            wb.close()
        except Exception:
            pass


def _header_index(header: list, column: str, where: str) -> int:
    wanted = str(column).strip().lower()
    for index, name in enumerate(header):
        if str(name or "").strip().lower() == wanted:
            return index
    raise OfficeError(
        f"no column {column!r} in the {where} table. Columns: "
        + ", ".join(str(h) for h in header if str(h or "").strip()))


def _comparable(value: Any, profile: dict) -> Any:
    """Normalise one value for comparison under its column's convention."""
    kind = profile.get("kind")
    convention = profile.get("convention", "")
    if kind == "number":
        parsed = parse_number(value, convention)
        return parsed if parsed is not None else str(value or "").strip()
    if kind == "date":
        parsed = parse_date(value, convention)
        return parsed if parsed is not None else str(value or "").strip()
    return str(value or "").strip()


def compare_tables(
    left: Any,
    right: Any,
    *,
    key: str,
    right_key: Optional[str] = None,
    columns: Optional[list[str]] = None,
    left_sheet: Optional[str] = None,
    right_sheet: Optional[str] = None,
    max_report: int = 200,
) -> dict:
    """Reconcile two tables on a key column.

    Every row of both inputs lands in exactly one group — equal,
    differing, only-left, only-right, or not-comparable — and the counts
    add up to the inputs. That is the whole point: the failure mode of
    doing this by hand is not a wrong comparison, it is rows that quietly
    fall out of the report and are never noticed.

    Values are compared under the convention their column actually uses,
    so "1.234,50" and 1234.5 are equal and 31.07.2026 matches 2026-07-31.
    Duplicate keys are reported rather than joined: a duplicate key turns
    a comparison into a cross product, and the result looks plausible.
    """
    left_rows = _table_rows(left, left_sheet)
    right_rows = _table_rows(right, right_sheet)
    if len(left_rows) < 2 or len(right_rows) < 2:
        raise OfficeError(
            "both tables need a header row and at least one data row")

    left_header, right_header = left_rows[0], right_rows[0]
    # The two tables come from two systems, so the key is routinely spelled
    # differently on each side — "Beleg" against "Belegnummer". Requiring
    # one name would push the caller into renaming a column first, which
    # means writing to a file just to be able to read it.
    right_key_name = str(right_key).strip() if right_key else key
    left_key_index = _header_index(left_header, key, "left")
    right_key_index = _header_index(right_header, right_key_name, "right")

    # ``columns`` may be a list (same name on both sides) or a mapping
    # {left: right}. Two exports of the same facts routinely disagree on
    # every column name — "Betrag" against "Rechnungsbetrag" — and without
    # pairs those columns simply drop out of the comparison, which reports
    # agreement it never checked.
    pairs: dict[str, str] = {}
    if isinstance(columns, dict):
        pairs = {str(k): str(v) for k, v in columns.items()}
        wanted = list(pairs)
    elif columns:
        wanted = [str(c) for c in columns]
    else:
        # Default: every column both tables share, key excluded.
        right_names = {str(h or "").strip().lower() for h in right_header}
        wanted = [
            str(h) for h in left_header
            if str(h or "").strip()
            and str(h or "").strip().lower() in right_names
            and str(h or "").strip().lower() not in (
                str(key).strip().lower(), right_key_name.strip().lower())
        ]
    if not wanted:
        raise OfficeError(
            "the two tables share no comparable column besides the key — "
            "name the columns explicitly")

    left_idx = {c: _header_index(left_header, c, "left") for c in wanted}
    right_idx = {
        c: _header_index(right_header, pairs.get(c, c), "right")
        for c in wanted
    }

    # Each side is profiled on its own. The two tables routinely come from
    # different systems and are written under different conventions —
    # 1.234,50 out of a German export against 1234.50 out of a database.
    # Applying one side's convention to both would report every such pair
    # as a difference, which is exactly the noise that makes a
    # reconciliation useless.
    left_profiles = {
        c: profile_column([r[i] if i < len(r) else ""
                           for r in left_rows[1:]], name=c)
        for c, i in left_idx.items()
    }
    right_profiles = {
        c: profile_column([r[i] if i < len(r) else ""
                           for r in right_rows[1:]], name=pairs.get(c, c))
        for c, i in right_idx.items()
    }

    def _index(rows: list[list], key_at: int) -> tuple[dict, list, list]:
        by_key: dict[str, list] = {}
        blank: list[int] = []
        for offset, row in enumerate(rows[1:], start=2):
            raw = row[key_at] if key_at < len(row) else ""
            token = str(raw or "").strip()
            if not token:
                blank.append(offset)
                continue
            by_key.setdefault(token, []).append((offset, row))
        duplicates = [
            {"key": k, "rows": [o for o, _ in v]}
            for k, v in by_key.items() if len(v) > 1
        ]
        return by_key, duplicates, blank

    left_by_key, left_dupes, left_blank = _index(left_rows, left_key_index)
    right_by_key, right_dupes, right_blank = _index(right_rows, right_key_index)

    not_comparable: list[dict] = []
    for line in left_blank:
        not_comparable.append(
            {"side": "left", "row": line, "reason": "empty key"})
    for line in right_blank:
        not_comparable.append(
            {"side": "right", "row": line, "reason": "empty key"})
    for entry in left_dupes:
        not_comparable.append({
            "side": "left", "key": entry["key"], "rows": entry["rows"],
            "reason": "key appears more than once"})
    for entry in right_dupes:
        not_comparable.append({
            "side": "right", "key": entry["key"], "rows": entry["rows"],
            "reason": "key appears more than once"})

    duplicate_keys = ({e["key"] for e in left_dupes}
                      | {e["key"] for e in right_dupes})
    comparable_left = {k: v[0] for k, v in left_by_key.items()
                       if k not in duplicate_keys}
    comparable_right = {k: v[0] for k, v in right_by_key.items()
                        if k not in duplicate_keys}

    equal: list[str] = []
    differing: list[dict] = []
    for token, (line, row) in sorted(comparable_left.items()):
        if token not in comparable_right:
            continue
        other_line, other_row = comparable_right[token]
        diffs = []
        for column in wanted:
            a = row[left_idx[column]] if left_idx[column] < len(row) else ""
            b = (other_row[right_idx[column]]
                 if right_idx[column] < len(other_row) else "")
            left_p, right_p = left_profiles[column], right_profiles[column]
            if left_p.get("kind") == right_p.get("kind"):
                # Same kind on both sides: normalise each under its OWN
                # convention, so the comparison is of values, not spellings.
                same = _comparable(a, left_p) == _comparable(b, right_p)
            else:
                # One side is text where the other is a number or a date.
                # Comparing the normalised forms would be comparing two
                # different things, so fall back to the literal text.
                same = str(a or "").strip() == str(b or "").strip()
            if not same:
                diffs.append({"column": column,
                              "left": _fmt(a), "right": _fmt(b)})
        if diffs:
            differing.append({"key": token, "left_row": line,
                              "right_row": other_line, "differences": diffs})
        else:
            equal.append(token)

    only_left = sorted(set(comparable_left) - set(comparable_right))
    only_right = sorted(set(comparable_right) - set(comparable_left))

    notes: list[str] = []

    # Scope, before anything else. A workbook with one sheet per month
    # compared against a whole year reports hundreds of one-sided rows —
    # a result that reads like a catastrophe when the only thing wrong is
    # which sheet was taken. The active sheet is a silent default, so it
    # has to be said out loud.
    for side, source, chosen in (("left", left, left_sheet),
                                 ("right", right, right_sheet)):
        names = _workbook_sheets(source)
        if len(names) > 1 and not chosen:
            notes.append(
                f"the {side} workbook has {len(names)} sheets and none was "
                f"named, so its ACTIVE sheet was used: {names[0]!r} of "
                f"{', '.join(names)}. Check this is the intended scope."
            )

    for column in wanted:
        for side, profile in (("left", left_profiles[column]),
                              ("right", right_profiles[column])):
            if profile.get("convention") == AMBIGUOUS:
                notes.append(
                    f"column '{column}' ({side}): "
                    f"{profile.get('convention_reason', '')} — the comparison "
                    "treated its values as text.")
        if left_profiles[column].get("kind") != right_profiles[column].get("kind"):
            notes.append(
                f"column '{column}': the left side looks like "
                f"{left_profiles[column].get('kind')}s and the right side like "
                f"{right_profiles[column].get('kind')}s — compared as text, so "
                "a difference here may be a formatting difference only.")
    if duplicate_keys:
        notes.append(
            f"{len(duplicate_keys)} key(s) appear more than once and were NOT "
            "compared — a duplicate key would make the join a cross product. "
            "Resolve them or pick a different key.")
    if left_blank or right_blank:
        notes.append(
            f"{len(left_blank) + len(right_blank)} row(s) have an empty key "
            "and could not be matched.")

    # A gross size asymmetry is far more often the wrong scope than a real
    # gap of that size. Say which reading is the likelier one rather than
    # letting a count of hundreds stand as a finding.
    matched = len(equal) + len(differing)
    for side, rows, one_sided in (("left", len(left_rows) - 1, only_left),
                                  ("right", len(right_rows) - 1, only_right)):
        if rows and len(one_sided) > max(10, matched):
            notes.append(
                f"{len(one_sided)} of {rows} rows exist only on the {side} "
                "side, against " + str(matched) + " matched. A difference of "
                "that size is usually a scope mismatch — a different period, "
                "sheet or filter — rather than that many missing records. "
                "Confirm both tables cover the same range before reporting "
                "this.")

    # Every input row must be accounted for. If this ever fails, the
    # report is wrong in a way the caller could not see.
    accounted = (len(equal) + len(differing) + len(only_left)
                 + len(only_right) + len(left_blank) + len(right_blank)
                 + sum(len(e["rows"]) for e in left_dupes)
                 + sum(len(e["rows"]) for e in right_dupes))
    expected = (len(left_rows) - 1) + (len(right_rows) - 1)
    # Matched pairs consume one row from each side.
    balanced = accounted + len(equal) + len(differing) == expected

    return {
        "left": str(_resolve(left)),
        "right": str(_resolve(right)),
        "key": key,
        "right_key": right_key_name,
        "compared_columns": [
            c if pairs.get(c, c) == c else f"{c} / {pairs[c]}" for c in wanted],
        "left_rows": len(left_rows) - 1,
        "right_rows": len(right_rows) - 1,
        "equal": equal[:max_report],
        "equal_count": len(equal),
        "differing": differing[:max_report],
        "differing_count": len(differing),
        "only_left": only_left[:max_report],
        "only_left_count": len(only_left),
        "only_right": only_right[:max_report],
        "only_right_count": len(only_right),
        "not_comparable": not_comparable[:max_report],
        "not_comparable_count": len(not_comparable),
        "rows_accounted_for": balanced,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# Series: one document per row
# ---------------------------------------------------------------------------

_UNSAFE_NAME_RE = re.compile(r"[^\w.\- ()äöüÄÖÜß]+")

MAX_SERIES_ROWS = 500


def _safe_name(text: str) -> str:
    """A file name that cannot escape its directory or collide with shell."""
    cleaned = _UNSAFE_NAME_RE.sub("_", str(text or "").strip())
    # Collapse dot runs: separators are already gone, so "../etc" cannot
    # escape, but a name still carrying ".." invites doubt every time
    # someone reads the output listing.
    cleaned = re.sub(r"\.{2,}", ".", cleaned)
    cleaned = cleaned.strip(" .") or "unbenannt"
    return cleaned[:120]


def _series_name(pattern: str, record: dict, index: int, suffix: str) -> str:
    """Build one output file name from the row's own values."""
    out = pattern
    for column, value in record.items():
        out = out.replace("{" + column + "}", _safe_name(_fmt(value)))
    out = out.replace("{row}", str(index))
    out = _UNSAFE_NAME_RE.sub("_", out).strip(" .") or f"dokument_{index}"
    if not out.lower().endswith(suffix):
        out += suffix
    return out[:150]


def fill_series(
    table: Any,
    template: Any,
    *,
    output_dir: Any,
    mapping: Optional[dict] = None,
    constants: Optional[dict] = None,
    name_pattern: str = "",
    sheet: Optional[str] = None,
    limit: int = MAX_SERIES_ROWS,
) -> dict:
    """One filled document per row of *table*, from one form or template.

    The whole point is the report. A series that ends in "done" while six
    rows quietly failed is the expensive kind of mistake in
    administrative work, so every row gets its own outcome and the counts
    are what the caller reports.

    Everything that can be checked is checked BEFORE the first file is
    written: the template's fields, the table's columns, the mapping
    between them. A batch that fails halfway leaves a directory of
    documents nobody can tell apart from a complete one.

    Nothing is ever overwritten. A name that already exists — on disk or
    from an earlier row — fails that row and says which row took it,
    because two different records writing to one file is data loss that
    looks like success.
    """
    rows = _table_rows(table, sheet)
    if len(rows) < 2:
        raise OfficeError(
            "the table needs a header row and at least one data row")
    header = [str(h or "").strip() for h in rows[0]]
    if not any(header):
        raise OfficeError("the table's first row is empty — no column names")

    tpl = _resolve(template)
    _refuse_unwritable_format(tpl, what="filling it as a template")
    kind = document_kind(tpl)
    if kind == "pdf":
        available = {f["name"] for f in pdf_form_fields(tpl)["fields"]}
        if not available:
            raise OfficeError(
                f"{tpl.name} has no form fields to fill")
        suffix = ".pdf"
    elif kind == "word":
        available = {p["name"]
                     for p in docx_placeholders(tpl)["placeholders"]}
        if not available:
            raise OfficeError(
                f"{tpl.name} contains no {{{{placeholder}}}} markers")
        suffix = ".docx"
    else:
        raise OfficeError(
            f"{tpl.name} is neither a PDF form nor a .docx template")

    # Mapping field -> column. Without one, match by name: a template
    # built for this table usually already uses the column names.
    if mapping:
        resolved_map = {str(k): str(v) for k, v in mapping.items()}
    else:
        lowered = {h.lower(): h for h in header if h}
        resolved_map = {
            field: lowered[field.lower()]
            for field in available
            if field.lower() in lowered and field not in (constants or {})
        }
        if not resolved_map:
            raise OfficeError(
                "no field name matches a column name, so there is nothing to "
                f"map automatically. Template fields: {', '.join(sorted(available))}. "
                f"Columns: {', '.join(h for h in header if h)}. "
                "Pass mapping={field: column}.")

    # Fixed values: the same entry in every document — a file reference, a
    # date, a department. Without them a caller has to build a derived
    # table with the value copied into all 400 rows, which is what was
    # observed in the field: two constants cost a whole intermediate file.
    fixed = {str(k): v for k, v in (constants or {}).items()}
    both = sorted(set(fixed) & set(resolved_map))
    if both:
        raise OfficeError(
            f"field(s) {', '.join(both)} are given both a column and a fixed "
            "value — one of the two has to go")

    unknown_fields = sorted(
        f for f in list(resolved_map) + list(fixed) if f not in available)
    if unknown_fields:
        raise OfficeError(
            f"the template has no field(s) {', '.join(unknown_fields)}. "
            f"It has: {', '.join(sorted(available))}")
    unknown_columns = sorted(
        c for c in resolved_map.values() if c not in header)
    if unknown_columns:
        raise OfficeError(
            f"the table has no column(s) {', '.join(unknown_columns)}. "
            f"It has: {', '.join(h for h in header if h)}. If a value is the "
            "same for every row, pass it in constants instead of adding a "
            "column for it.")

    out_dir = _resolve(output_dir, must_exist=False)
    if out_dir.exists() and not out_dir.is_dir():
        raise OfficeError(f"{out_dir} exists and is not a directory")
    out_dir.mkdir(parents=True, exist_ok=True)

    pattern = name_pattern or (tpl.stem + "_{row}" + suffix)
    limit = max(1, min(int(limit or MAX_SERIES_ROWS), MAX_SERIES_ROWS))

    data_rows = rows[1:]
    truncated = len(data_rows) > limit
    used_names: dict[str, int] = {}
    results: list[dict] = []

    for offset, row in enumerate(data_rows[:limit], start=2):
        record = {
            header[i]: (row[i] if i < len(row) else "")
            for i in range(len(header)) if header[i]
        }
        values = {field: record.get(column, "")
                  for field, column in resolved_map.items()}
        values.update(fixed)
        empty = sorted(f for f, v in values.items()
                       if str(v or "").strip() == "")

        name = _series_name(pattern, record, offset, suffix)
        if name in used_names:
            results.append({
                "row": offset, "output": name, "status": "failed",
                "detail": (f"the file name is already used by row "
                           f"{used_names[name]} — two records would write "
                           "to one file"),
            })
            continue
        target = out_dir / name
        if target.exists():
            results.append({
                "row": offset, "output": name, "status": "failed",
                "detail": "a file of that name already exists; nothing was "
                          "overwritten",
            })
            continue

        try:
            if kind == "pdf":
                filled = fill_pdf_form(
                    tpl, {k: v for k, v in values.items()}, output=target)
                ok = filled["verified"]
                detail = "; ".join(filled.get("notes", []))
            else:
                filled = fill_docx_template(
                    tpl, {k: v for k, v in values.items()}, output=target,
                    strict=False)
                ok = filled["complete"]
                detail = "; ".join(filled.get("notes", []))
        except OfficeError as exc:
            results.append({"row": offset, "output": name,
                            "status": "failed", "detail": str(exc)})
            continue
        except Exception as exc:                      # pragma: no cover
            results.append({"row": offset, "output": name,
                            "status": "failed", "detail": str(exc)})
            continue

        used_names[name] = offset
        if empty:
            results.append({
                "row": offset, "output": name, "status": "incomplete",
                "detail": "no value for " + ", ".join(empty),
            })
        elif not ok:
            results.append({"row": offset, "output": name,
                            "status": "incomplete",
                            "detail": detail or "the result did not verify"})
        else:
            results.append({"row": offset, "output": name, "status": "ok"})

    counts = {
        "ok": sum(1 for r in results if r["status"] == "ok"),
        "incomplete": sum(1 for r in results if r["status"] == "incomplete"),
        "failed": sum(1 for r in results if r["status"] == "failed"),
    }
    notes: list[str] = []
    if truncated:
        notes.append(
            f"the table has {len(data_rows)} data rows and only the first "
            f"{limit} were processed. Raise limit, or narrow the table.")
    if counts["incomplete"]:
        notes.append(
            f"{counts['incomplete']} document(s) are incomplete — they exist "
            "but are missing values. They are not ready to hand over.")
    if counts["failed"]:
        notes.append(
            f"{counts['failed']} row(s) produced no document.")
    return {
        "table": str(_resolve(table)),
        "template": str(tpl),
        "output_dir": str(out_dir),
        "mapping": resolved_map,
        "constants": fixed,
        "rows": len(data_rows),
        "processed": len(results),
        "counts": counts,
        "results": results,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# PDF
# ---------------------------------------------------------------------------

def _pdf_reader(path: Path) -> Any:
    pypdf = _require("pdf")
    try:
        return pypdf.PdfReader(str(path))
    except Exception as exc:
        raise OfficeError(f"could not open PDF: {exc}") from exc


def _acroform(reader: Any) -> Optional[Any]:
    """The document's AcroForm dictionary, resolved, or None."""
    try:
        root = reader.trailer["/Root"]
        acro = root.get("/AcroForm")
        if acro is None:
            return None
        try:
            return acro.get_object()
        except Exception:
            return acro
    except Exception:
        return None


def _has_xfa(reader: Any) -> bool:
    """True if the document is an XFA form.

    XFA (LiveCycle) forms keep their data in an XML stream, not in the
    AcroForm field dictionaries. Writing AcroForm values on such a file
    reports success and changes nothing the user will see, so this has
    to be detected before a fill is claimed to have worked.
    """
    acro = _acroform(reader)
    try:
        return acro is not None and "/XFA" in acro
    except Exception:
        return False


def _is_dynamic_xfa(reader: Any) -> bool:
    """True for a *dynamic* XFA form, where the AcroForm is only a stub.

    The distinction decides whether a fill is worth attempting at all. A
    static (hybrid) XFA form carries real AcroForm fields and renders
    them; a dynamic one is laid out from the XML at open time and the
    catalogue says so with ``/NeedsRendering``. Writing AcroForm values
    into a dynamic form produces a file whose data is there and whose
    every page still shows blank.
    """
    try:
        root = reader.trailer["/Root"]
        flag = root.get("/NeedsRendering")
        try:
            flag = flag.get_object()
        except Exception:
            pass
        # Not ``bool(flag)``: a pypdf BooleanObject is an object, and
        # every object is truthy, so a /NeedsRendering false would read
        # as true.
        return str(getattr(flag, "value", flag)).strip().lower() == "true"
    except Exception:
        return False


# The XFA packets that hold something a reader can use. ``datasets`` is
# the filled-in data; ``template`` is the layout, which is far too large
# to be worth returning and holds no values.
_XFA_DATA_PACKETS = ("datasets", "xdp:xdp", "xdp")

MAX_XFA_VALUES = 300


def _xfa_packets(reader: Any) -> dict[str, bytes]:
    """The XFA packets of a form, by name. Empty for anything else.

    ``/XFA`` is either one stream holding the whole XDP or a flat array
    alternating names and streams. Both shapes occur in the wild, and a
    reader that assumes one of them silently returns nothing for the
    other.
    """
    acro = _acroform(reader)
    if acro is None:
        return {}
    try:
        xfa = acro.get("/XFA")
        xfa = xfa.get_object() if hasattr(xfa, "get_object") else xfa
    except Exception:
        return {}
    if xfa is None:
        return {}

    def _bytes(entry: Any) -> bytes:
        try:
            entry = entry.get_object()
        except Exception:
            pass
        for attr in ("get_data", "getData"):
            reader_fn = getattr(entry, attr, None)
            if callable(reader_fn):
                try:
                    data = reader_fn()
                    return (data if isinstance(data, bytes)
                            else str(data).encode())
                except Exception:
                    return b""
        if isinstance(entry, bytes):
            return entry
        return str(entry).encode("utf-8", "replace")

    packets: dict[str, bytes] = {}
    if isinstance(xfa, (list, tuple)):
        items = list(xfa)
        # Name/stream pairs; an odd tail is a malformed array and is
        # filed under its position rather than dropped silently.
        for index in range(0, len(items), 2):
            name = str(items[index])
            body = items[index + 1] if index + 1 < len(items) else None
            packets[name if body is not None else f"packet{index}"] = _bytes(
                body if body is not None else items[index])
    else:
        packets["xdp"] = _bytes(xfa)
    return packets


def _xfa_values(packets: dict[str, bytes]) -> list[dict]:
    """Field path -> value out of the XFA ``datasets`` packet.

    This is the only way to see what a dynamic XFA form actually holds:
    its values never reach the AcroForm dictionaries, so the field
    listing reports every field as empty while the form is full.
    """
    import xml.etree.ElementTree as ET

    candidates = [packets[name] for name in _XFA_DATA_PACKETS
                  if packets.get(name)]
    candidates += [body for name, body in packets.items()
                   if name not in _XFA_DATA_PACKETS and b"datasets" in body]
    out: list[dict] = []
    seen: set[str] = set()
    for body in candidates:
        try:
            root = ET.fromstring(body.decode("utf-8", "replace"))
        except Exception:
            continue
        data = root
        for element in root.iter():
            if element.tag.rsplit("}", 1)[-1] == "data":
                data = element
                break

        def walk(node: Any, path: str) -> None:
            children = list(node)
            if not children:
                text = (node.text or "").strip()
                if path and text and path not in seen:
                    seen.add(path)
                    out.append({"name": path, "value": text})
                return
            counts: dict[str, int] = {}
            for child in children:
                tag = child.tag.rsplit("}", 1)[-1]
                counts[tag] = counts.get(tag, 0) + 1
            used: dict[str, int] = {}
            for child in children:
                tag = child.tag.rsplit("}", 1)[-1]
                if counts[tag] > 1:
                    used[tag] = used.get(tag, 0) + 1
                    label = f"{tag}[{used[tag]}]"
                else:
                    label = tag
                walk(child, f"{path}.{label}" if path else label)

        for child in list(data):
            walk(child, child.tag.rsplit("}", 1)[-1])
        if out:
            break
    return out[:MAX_XFA_VALUES]


def pdf_xfa_data(path: Any) -> dict:
    """The values stored in an XFA form's dataset, so they can be read.

    Reading is all this offers. Writing the dataset back would mean
    rewriting an XML the form's own script layer validates, and a form
    that opens with the right numbers in the wrong places is exactly the
    outcome a refusal is worth more than.
    """
    p = _resolve(path)
    if document_kind(p) != "pdf":
        raise OfficeError(f"{p.name} is not a PDF")
    reader = _pdf_reader(p)
    if not _has_xfa(reader):
        raise OfficeError(
            f"{p.name} is not an XFA form — list its fields instead")
    packets = _xfa_packets(reader)
    values = _xfa_values(packets)
    notes: list[str] = []
    if not values:
        notes.append(
            "the XFA packets carry no filled-in dataset: either the form is "
            "empty or its data sits in a packet this reader does not "
            "understand (packets present: "
            + (", ".join(sorted(packets)) or "none") + ")")
    if len(values) == MAX_XFA_VALUES:
        notes.append(
            f"only the first {MAX_XFA_VALUES} values are listed")
    notes.append(
        "these values were read from the form's XML dataset, not from "
        "AcroForm fields. They cannot be written back here.")
    return {
        "path": str(p),
        "packets": sorted(packets),
        "dynamic": _is_dynamic_xfa(reader),
        "values": values,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# Scanned pages
# ---------------------------------------------------------------------------
#
# A PDF whose pages hold no text is not a broken PDF, it is a photograph
# of paper — the single most common attachment in an administrative
# in-tray. Two things matter here and neither is the OCR itself:
#
# * the diagnosis has to be specific. "OCR is not available" leaves the
#   reader with no way to change that; naming the component that is
#   missing and the command that installs it does.
# * text produced by OCR is a reading of pixels, not content of the
#   file. It is labelled as such everywhere it appears, and it never
#   replaces a page that has a real text layer.

OCR_MAX_PAGES = 3
OCR_DPI = 200
OCR_LANGUAGES = ("de", "en")
# Below this the recogniser is guessing at the shape of the glyphs. Such
# blocks are counted and reported rather than quietly dropped: a missing
# line is easier to notice than a wrong one.
OCR_LOW_CONFIDENCE = 0.5


def _easyocr_model_dir() -> Path:
    import os
    base = (os.environ.get("EASYOCR_MODULE_PATH")
            or os.environ.get("MODULE_PATH")
            or str(Path.home() / ".EasyOCR"))
    return Path(base).expanduser() / "model"


def ocr_availability() -> dict:
    """What OCR would need on this machine. Never raises.

    Returns ``{"available", "engine", "detail", "next_step"}``. The
    detail names the component that is missing rather than the
    capability, because "install tesseract-ocr" is actionable and "OCR
    unavailable" is not.
    """
    import importlib
    import shutil

    detail: list[str] = []
    renderer = False
    try:
        importlib.import_module("fitz")
        renderer = True
    except Exception:
        detail.append(
            "pymupdf is not installed, so a PDF page cannot be turned into "
            "an image at all (pip install pymupdf)")

    engine = ""
    try:
        importlib.import_module("pytesseract")
    except Exception:
        detail.append("pytesseract is not installed")
    else:
        if shutil.which("tesseract"):
            engine = engine or "pytesseract"
        else:
            detail.append(
                "pytesseract is installed but the 'tesseract' program it "
                "drives is not on PATH — pytesseract is only a wrapper, the "
                "recognition happens in that binary (Debian/Ubuntu: "
                "apt-get install tesseract-ocr tesseract-ocr-deu)")

    try:
        importlib.import_module("easyocr")
        importlib.import_module("torch")
    except Exception as exc:
        detail.append(f"easyocr is not usable ({type(exc).__name__})")
    else:
        weights = _easyocr_model_dir()
        try:
            present = len(list(weights.glob("*.pth")))
        except OSError:
            present = 0
        if present >= 2:
            engine = engine or "easyocr"
        else:
            detail.append(
                "easyocr is installed but its recognition models are not in "
                f"{weights} — the first call would download roughly 100 MB "
                "without asking, which a document read must not do. Fetch "
                "them once with network access: python -c \"import easyocr; "
                "easyocr.Reader(['de','en'])\"")

    available = bool(engine) and renderer
    if available:
        next_step = (
            f"pass ocr=true to read the scanned pages with {engine}; the "
            "result is a reading of the page image, not text out of the "
            "file")
    elif engine and not renderer:
        next_step = "install pymupdf, then pass ocr=true"
    else:
        next_step = (
            "install one of the two engines above, then pass ocr=true")
    return {
        "available": available,
        "engine": engine if available else "",
        "renderer": "pymupdf" if renderer else "",
        "detail": detail,
        "next_step": next_step,
    }


def _render_page(path: Path, number: int, dpi: int = OCR_DPI) -> bytes:
    """One PDF page as PNG bytes. Raises OfficeError if it cannot."""
    try:
        import fitz
    except Exception as exc:
        raise OfficeError(
            f"rendering a page needs pymupdf, which is not installed ({exc})"
        ) from exc
    try:
        document = fitz.open(str(path))
    except Exception as exc:
        raise OfficeError(f"could not open {path.name} for rendering: {exc}"
                          ) from exc
    try:
        return document[number - 1].get_pixmap(dpi=dpi).tobytes("png")
    except Exception as exc:
        raise OfficeError(
            f"page {number} of {path.name} could not be rendered: {exc}"
        ) from exc
    finally:
        try:
            document.close()
        except Exception:
            pass


_OCR_READER_CACHE: dict[tuple, Any] = {}


def _ocr_image(
    image: bytes, languages: tuple = OCR_LANGUAGES
) -> tuple[str, float, int]:
    """OCR one page image. Returns ``(text, mean confidence, weak blocks)``."""
    status = ocr_availability()
    if not status["available"] or status["engine"] != "easyocr":
        raise OfficeError(
            "OCR is not available here: "
            + ("; ".join(status["detail"]) or "no engine could be used")
            + ". " + status["next_step"])
    try:
        import easyocr
    except Exception as exc:                      # pragma: no cover
        raise OfficeError(f"easyocr could not be imported: {exc}") from exc

    key = tuple(languages)
    reader = _OCR_READER_CACHE.get(key)
    if reader is None:
        try:
            # Never on the GPU: this runs inside an interactive tool call
            # on machines that may be running something else on the card,
            # and a CUDA allocation failure would surface as an OCR
            # failure rather than as what it is.
            reader = easyocr.Reader(list(languages), gpu=False, verbose=False,
                                    download_enabled=False)
        except Exception as exc:
            raise OfficeError(
                f"the OCR engine could not be started: {exc}") from exc
        _OCR_READER_CACHE[key] = reader
    try:
        found = reader.readtext(image, detail=1)
    except Exception as exc:
        raise OfficeError(f"OCR failed on this page: {exc}") from exc

    lines = [str(text) for _box, text, _conf in found]
    scores = [float(conf) for _box, _text, conf in found]
    weak = sum(1 for s in scores if s < OCR_LOW_CONFIDENCE)
    mean = sum(scores) / len(scores) if scores else 0.0
    return "\n".join(lines), mean, weak


def _page_image_count(page: Any) -> int:
    """How many images a page draws, read from its resources.

    Deliberately not ``page.images``, which decodes every one of them: a
    300 dpi scan would be unpacked into memory only to be counted.
    """
    try:
        resources = page.get("/Resources")
        resources = (resources.get_object()
                     if hasattr(resources, "get_object") else resources)
        xobjects = resources.get("/XObject") if resources else None
        xobjects = (xobjects.get_object()
                    if hasattr(xobjects, "get_object") else xobjects)
        if not xobjects:
            return 0
        count = 0
        for key in xobjects:
            try:
                entry = xobjects[key].get_object()
            except Exception:
                continue
            if str(entry.get("/Subtype", "")) == "/Image":
                count += 1
        return count
    except Exception:
        return 0


def read_pdf(
    path: Any,
    *,
    pages: Optional[str] = None,
    max_chars: int = MAX_TEXT_CHARS,
    ocr: bool = False,
) -> dict:
    """Extract text from a PDF, page by page.

    ``pages`` accepts ``"3"``, ``"2-5"`` or ``"1,4,7"``; omitted means
    the whole document up to *max_chars*.

    ``ocr`` reads the pages that carry no text layer from their page
    image instead. It is off by default and it never touches a page that
    has real text: OCR output is a reading of pixels, and replacing text
    the file actually contains with a guess at it would trade a correct
    answer for a plausible one. Blocks it recognised are labelled as OCR
    in the output, with the confidence the engine reported.
    """
    p = _resolve(path)
    if document_kind(p) != "pdf":
        raise OfficeError(f"{p.name} is not a PDF")
    reader = _pdf_reader(p)
    total = len(reader.pages)
    wanted = _parse_pages(pages, total)

    extracted: dict[int, str] = {}
    empty: list[int] = []
    scanned: list[int] = []
    for number in wanted:
        try:
            page = reader.pages[number - 1]
            text = page.extract_text() or ""
        except Exception as exc:
            text = f"(page {number}: extraction failed: {exc})"
            page = None
        extracted[number] = text.strip()
        if not extracted[number]:
            empty.append(number)
            # A page with images and no text is a scan; a page with
            # neither is a blank sheet. OCR helps with the first and has
            # nothing to find on the second, and saying which is which
            # is the difference between a next step and a guess.
            if page is not None and _page_image_count(page):
                scanned.append(number)

    notes: list[str] = []
    ocr_pages: list[dict] = []
    if ocr and scanned:
        budget = scanned[:OCR_MAX_PAGES]
        for number in budget:
            try:
                image = _render_page(p, number)
                text, confidence, weak = _ocr_image(image)
            except OfficeError as exc:
                notes.append(f"page {number}: {exc}")
                break
            extracted[number] = text
            ocr_pages.append({"page": number,
                              "confidence": round(confidence, 3),
                              "low_confidence_blocks": weak})
        if ocr_pages:
            mean = sum(e["confidence"] for e in ocr_pages) / len(ocr_pages)
            weak_total = sum(e["low_confidence_blocks"] for e in ocr_pages)
            notes.append(
                "page(s) " + ", ".join(str(e["page"]) for e in ocr_pages)
                + f" were read by OCR ({ocr_availability()['engine']}, mean "
                f"confidence {mean:.2f}). That text is a reading of the page "
                "image, not content of the file: check every figure, name and "
                "reference number against the original before acting on it.")
            if weak_total:
                notes.append(
                    f"{weak_total} recognised block(s) scored below "
                    f"{OCR_LOW_CONFIDENCE} — treat those lines as unread "
                    "rather than as text.")
        if len(scanned) > len(budget):
            notes.append(
                f"{len(scanned)} page(s) have no text layer and only "
                f"{OCR_MAX_PAGES} were run through OCR (it takes tens of "
                "seconds per page). Pass pages to choose which.")
    elif ocr and not scanned:
        notes.append(
            "ocr was requested but no requested page is a scan — every page "
            "either has a text layer already or holds nothing at all.")

    chunks: list[str] = []
    used = 0
    truncated = False
    ocr_numbers = {e["page"]: e for e in ocr_pages}
    for number in wanted:
        entry = ocr_numbers.get(number)
        label = (f"--- page {number} (OCR, confidence "
                 f"{entry['confidence']:.2f}) ---" if entry
                 else f"--- page {number} ---")
        block = f"{label}\n{extracted[number]}"
        if used + len(block) > max_chars:
            chunks.append(block[: max(0, max_chars - used)])
            truncated = True
            break
        chunks.append(block)
        used += len(block)

    still_empty = [n for n in empty if n not in ocr_numbers]
    if truncated:
        notes.append(
            f"output truncated at {max_chars} characters — pass pages to "
            "read a specific range")
    if still_empty:
        status = ocr_availability()
        blank = [n for n in still_empty if n not in scanned]
        if len(still_empty) == len(wanted):
            head = ("no extractable text on any requested page"
                    if len(wanted) > 1 else "no extractable text on this page")
        else:
            head = f"{len(still_empty)} of {len(wanted)} pages held no text"
        if [n for n in still_empty if n in scanned]:
            head += (
                " — page(s) "
                + ", ".join(str(n) for n in still_empty if n in scanned)
                + " draw an image and carry no text layer, which is a scan")
            if status["available"]:
                head += ". " + status["next_step"]
            else:
                head += (". OCR is not possible on this machine: "
                         + ("; ".join(status["detail"]) or "no engine found")
                         + ". " + status["next_step"]
                         + ", or run the file through OCR elsewhere and read "
                         "the result.")
        if blank:
            head += (
                " — page(s) " + ", ".join(str(n) for n in blank)
                + " hold neither text nor an image and are simply empty; OCR "
                "would find nothing on them")
        notes.append(head)
    if _has_xfa(reader):
        dynamic = _is_dynamic_xfa(reader)
        found = len(_xfa_values(_xfa_packets(reader)))
        notes.append(
            ("this is a dynamic XFA form" if dynamic
             else "this is an XFA form")
            + " — its field data lives in an XML stream, not in the page text"
            + (f". {found} stored value(s) can be read: list the fields "
               "(fields=true) to see them" if found else
               ". No stored values were found in its dataset"))

    return {
        "path": str(p),
        "kind": "pdf",
        "pages": total,
        "pages_read": wanted,
        "pages_without_text": still_empty,
        "ocr_pages": ocr_pages,
        "text": "\n\n".join(chunks),
        "notes": notes,
    }


def _parse_pages(spec: Optional[str], total: int) -> list[int]:
    if not spec or not str(spec).strip():
        return list(range(1, total + 1))
    out: list[int] = []
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo_s, _, hi_s = part.partition("-")
            try:
                lo, hi = int(lo_s), int(hi_s)
            except ValueError as exc:
                raise OfficeError(f"bad page range {part!r}") from exc
            if lo > hi:
                lo, hi = hi, lo
            out.extend(range(lo, hi + 1))
        else:
            try:
                out.append(int(part))
            except ValueError as exc:
                raise OfficeError(f"bad page number {part!r}") from exc
    seen: set[int] = set()
    result = []
    for n in out:
        if 1 <= n <= total and n not in seen:
            seen.add(n)
            result.append(n)
    if not result:
        raise OfficeError(
            f"no valid page in {spec!r} — the document has {total} page(s)")
    return result


_FIELD_TYPES = {
    "/Tx": "text",
    "/Btn": "button/checkbox",
    "/Ch": "choice",
    "/Sig": "signature",
}


def _field_labels(path: Path) -> dict:
    """Pair each form field with the text printed next to it.

    A field name carries no meaning in many real forms — ``text1``,
    ``text2``, ``Kontrollkästchen3``. What identifies the field is the
    label printed beside it, and the two are connected only by their
    position on the page. Without that link a caller can do nothing but
    assume the field order matches the visual order, which in a form
    built with a designer it frequently does not: a value in the wrong
    field looks perfectly correct afterwards.

    Order of evidence: the field's own tooltip (``/TU``) when the form
    author set one, otherwise the nearest text to the LEFT on the same
    line, otherwise the text directly ABOVE. Returns
    ``{field: {"label", "source", "page"}}`` and never raises — a form
    whose geometry cannot be read simply yields no labels.
    """
    out: dict[str, dict] = {}
    try:
        import fitz
    except Exception:
        return out
    try:
        doc = fitz.open(str(path))
    except Exception:
        return out
    try:
        for number, page in enumerate(doc, start=1):
            try:
                widgets = list(page.widgets() or [])
            except Exception:
                continue
            if not widgets:
                continue
            try:
                words = page.get_text("words")     # x0, y0, x1, y1, word, ...
            except Exception:
                words = []
            for widget in widgets:
                name = str(getattr(widget, "field_name", "") or "")
                if not name or name in out:
                    continue
                tooltip = str(getattr(widget, "field_label", "") or "").strip()
                if tooltip:
                    out[name] = {"label": tooltip, "source": "tooltip",
                                 "page": number}
                    continue
                found = _nearest_label(widget.rect, words)
                if found:
                    out[name] = {"label": found[0], "source": found[1],
                                 "page": number}
    except Exception:
        return out
    finally:
        try:
            doc.close()
        except Exception:
            pass
    return out


def _nearest_label(rect: Any, words: list) -> Optional[tuple]:
    """The printed text that identifies a widget, plus where it sat.

    Left first: forms put the label on the same line, before the box.
    Above second, for stacked layouts. Both windows are deliberately
    narrow — a label further away than this is more likely to belong to
    a neighbouring field, and a wrong label is worse than none.
    """
    mid_y = (rect.y0 + rect.y1) / 2
    height = max(rect.y1 - rect.y0, 1.0)

    same_line = [
        w for w in words
        if w[2] <= rect.x0 + 2                      # ends left of the box
        and w[1] <= mid_y <= w[3] + height * 0.4    # vertically on the line
        and rect.x0 - w[2] < 320                    # not across the page
    ]
    if same_line:
        same_line.sort(key=lambda w: w[0])
        gap_limit = max(w[2] for w in same_line) - 260
        text = " ".join(w[4] for w in same_line if w[2] >= gap_limit)
        cleaned = text.strip(" :._-\t")
        if cleaned:
            return cleaned, "left"

    above = [
        w for w in words
        if w[3] <= rect.y0 + 2                      # ends above the box
        and rect.y0 - w[3] < height * 2.2           # directly above
        and w[2] > rect.x0 - 40 and w[0] < rect.x1  # horizontally overlapping
    ]
    if above:
        top = max(w[3] for w in above)
        text = " ".join(w[4] for w in sorted(above, key=lambda w: w[0])
                        if w[3] >= top - 2)
        cleaned = text.strip(" :._-\t")
        if cleaned:
            return cleaned, "above"
    return None


def pdf_form_fields(path: Any) -> dict:
    """List a PDF form's fields with their type, value and allowed states.

    The allowed states matter: a check box is not set with ``true`` but
    with the on-state name that field declares, and those differ between
    documents.
    """
    p = _resolve(path)
    if document_kind(p) != "pdf":
        raise OfficeError(f"{p.name} is not a PDF")
    reader = _pdf_reader(p)
    raw = {}
    try:
        raw = reader.get_fields() or {}
    except Exception as exc:
        raise OfficeError(f"could not read form fields: {exc}") from exc

    labels = _field_labels(p)
    fields = []
    for name, spec in raw.items():
        entry: dict[str, Any] = {
            "name": str(name),
            "type": _FIELD_TYPES.get(str(spec.get("/FT", "")), str(
                spec.get("/FT", "") or "unknown")),
            "value": _fmt(spec.get("/V")),
        }
        found = labels.get(str(name))
        if found:
            entry["label"] = found["label"]
            entry["label_source"] = found["source"]
            entry["page"] = found["page"]
        states = spec.get("/_States_")
        if states:
            entry["states"] = [str(s) for s in states]
        fields.append(entry)

    notes: list[str] = []
    unlabelled = [f["name"] for f in fields if not f.get("label")]
    if fields and unlabelled:
        notes.append(
            "no printed label could be located for: "
            + ", ".join(unlabelled[:10])
            + ". Do not infer them from the field order — in a form built "
            "with a designer it often does not follow the visual order.")
    xfa = _has_xfa(reader)
    dynamic = xfa and _is_dynamic_xfa(reader)
    xfa_values: list[dict] = []
    if xfa:
        try:
            xfa_values = _xfa_values(_xfa_packets(reader))
        except Exception:
            xfa_values = []
    if not fields and not xfa_values:
        notes.append(
            "this PDF has no form fields. Values cannot be filled in; "
            "producing a completed document means generating a new PDF "
            "or overlaying text.")
    if dynamic:
        notes.append(
            "dynamic XFA form (/NeedsRendering): the AcroForm entries above "
            "are a stub the viewer ignores — it builds the pages from the "
            "XML instead. Filling them is refused here rather than reported "
            "as done. To produce a completed form, fill it in a viewer that "
            "renders XFA, or rebuild it as a flat PDF.")
    elif xfa:
        notes.append(
            "XFA form: writing AcroForm field values will report success "
            "but will not change what a viewer displays. Confirm the "
            "result before handing this file over.")
    if xfa_values:
        # The values a dynamic form holds never reach the AcroForm
        # dictionaries, so the field listing above shows every one of
        # them as empty while the form is full. Reading them out of the
        # dataset is the only way to see what the document says.
        notes.append(
            f"{len(xfa_values)} value(s) were read out of the form's XML "
            "dataset (xfa_values). They are what this form actually holds; "
            "they cannot be written back here.")
    elif xfa:
        notes.append(
            "no values could be read from the XFA dataset — the form is "
            "either empty or stores its data in a packet this reader does "
            "not understand.")
    return {
        "path": str(p),
        "fields": fields,
        "field_count": len(fields),
        "xfa": xfa,
        "xfa_dynamic": dynamic,
        "xfa_values": xfa_values,
        "notes": notes,
    }


def fill_pdf_form(
    path: Any,
    values: dict,
    *,
    output: Any,
    flatten: bool = False,
) -> dict:
    """Fill an AcroForm and write the result to *output*.

    Behaviour that exists because the naive version silently fails:

    * ``/NeedAppearances`` is set on the form, so viewers that do not
      build appearance streams still show the values.
    * Boolean and loose check-box values are mapped to the on-state the
      individual field declares (``/Yes``, ``/On``, ``/1``, ...).
    * A field name that is not in the document is an error, not a
      silently ignored key — a typo would otherwise read as success.
    """
    src = _resolve(path)
    if document_kind(src) != "pdf":
        raise OfficeError(f"{src.name} is not a PDF")
    if not isinstance(values, dict) or not values:
        raise OfficeError("values must be a non-empty object of field -> value")
    out = _resolve(output, must_exist=False)
    if out.resolve() == src.resolve():
        raise OfficeError(
            "output must differ from the source — filling a form in place "
            "destroys the blank original")

    pypdf = _require("pdf")
    reader = _pdf_reader(src)
    # A dynamic XFA form is refused rather than filled. Its pages are
    # laid out from the XML at open time, so AcroForm values written
    # here are read by nothing: the call would report every field
    # written and the printed form would come out blank.
    if _has_xfa(reader) and _is_dynamic_xfa(reader):
        raise OfficeError(
            f"{src.name} is a dynamic XFA form (/NeedsRendering). Its pages "
            "are built from an XML layer that ignores the AcroForm fields, "
            "so filling those would report success and change nothing a "
            "viewer shows. Nothing was written. Fill it in a viewer that "
            "renders XFA, or ask for a flat PDF version of the form. Its "
            "current values can be read here (fields=true).")
    declared = {}
    try:
        declared = reader.get_fields() or {}
    except Exception:
        declared = {}
    if not declared:
        raise OfficeError(
            f"{src.name} has no form fields to fill — check the document "
            "with the field listing first")

    unknown = [k for k in values if str(k) not in declared]
    if unknown:
        raise OfficeError(
            f"no such field(s): {', '.join(sorted(unknown))}. Available: "
            f"{', '.join(sorted(str(k) for k in declared))}")

    resolved_values: dict[str, Any] = {}
    mapped: list[str] = []
    for key, value in values.items():
        spec = declared[str(key)]
        states = [str(s) for s in (spec.get("/_States_") or [])]
        if states:
            resolved = _checkbox_state(value, states)
            if resolved != value:
                mapped.append(f"{key}={value!r} -> {resolved}")
            resolved_values[str(key)] = resolved
        else:
            resolved_values[str(key)] = "" if value is None else str(value)

    try:
        writer = pypdf.PdfWriter(clone_from=str(src))
        _set_need_appearances(writer)
        _ensure_form_resources(writer)
        for page in writer.pages:
            writer.update_page_form_field_values(
                page, resolved_values, auto_regenerate=False)
        if flatten:
            _flatten_form(writer)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("wb") as fh:
            writer.write(fh)
    except OfficeError:
        raise
    except Exception as exc:
        raise OfficeError(f"could not fill the form: {exc}") from exc

    # Verify against the file that was actually written, rather than
    # reporting the intent as the outcome.
    written = {}
    try:
        written = {
            str(k): _fmt(v.get("/V"))
            for k, v in (pypdf.PdfReader(str(out)).get_fields() or {}).items()
        }
    except Exception:
        written = {}
    missing = [
        k for k, v in resolved_values.items()
        if written.get(k, "") != _fmt(v)
    ]

    notes: list[str] = []
    if mapped:
        notes.append("check-box values mapped to the field's own state: "
                     + "; ".join(mapped))
    if _has_xfa(reader):
        notes.append(
            "XFA form: the AcroForm values were written, but a viewer may "
            "render the XFA layer instead and show the form as empty. "
            "Open the result before treating this as done.")
    if missing:
        notes.append(
            "these fields did not hold the expected value when the written "
            "file was read back: " + ", ".join(sorted(missing)))
    if flatten:
        notes.append(
            "the fields were set to read-only, so the values can no longer "
            "be changed in a viewer. They remain form-field values and are "
            "not merged into the page content — text extraction will still "
            "not find them in the page text.")

    return {
        "source": str(src),
        "path": str(out),
        "filled": sorted(resolved_values),
        "verified": not missing,
        "notes": notes,
    }


def _checkbox_state(value: Any, states: list[str]) -> str:
    """Map a loose truthy/falsy value onto a field's declared states."""
    off = "/Off" if "/Off" in states else None
    on_states = [s for s in states if s != "/Off"]
    on = on_states[0] if on_states else "/Yes"
    if isinstance(value, str) and value in states:
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in ("true", "yes", "on", "x", "1", "ja"):
            return on
        if lowered in ("false", "no", "off", "", "0", "nein"):
            return off or "/Off"
        return value
    return on if value else (off or "/Off")


def _set_need_appearances(writer: Any) -> None:
    """Ask the viewer to build appearance streams for filled fields."""
    try:
        from pypdf.generic import BooleanObject, NameObject
        root = writer._root_object
        acro = root.get("/AcroForm")
        if acro is None:
            return
        try:
            acro = acro.get_object()
        except Exception:
            pass
        acro[NameObject("/NeedAppearances")] = BooleanObject(True)
    except Exception:
        # Cosmetic in viewers that render appearances themselves; never
        # worth failing an otherwise good fill.
        pass


def _ensure_form_resources(writer: Any) -> None:
    """Give the form the default-resources dictionary the writer expects.

    An AcroForm may legally omit ``/DR``, and some producers do. pypdf
    reads it while writing a field value and calls ``get_object()`` on the
    plain dictionary it substitutes, which fails with an AttributeError
    that names neither the form nor the field -- so a form that is simply
    missing an optional entry reads as a broken document.

    An empty but real ``/DR`` costs nothing here: ``/NeedAppearances`` is
    set, so the viewer builds the appearance itself.
    """
    try:
        from pypdf.generic import DictionaryObject, NameObject

        acro = writer._root_object.get("/AcroForm")
        if acro is None:
            return
        try:
            acro = acro.get_object()
        except Exception:
            pass

        def _resolved(container: Any, key: str) -> Any:
            value = container.get(key)
            try:
                return value.get_object() if value is not None else None
            except Exception:
                return None

        resources = _resolved(acro, "/DR")
        if not isinstance(resources, DictionaryObject):
            resources = DictionaryObject()
            acro[NameObject("/DR")] = resources
        if not isinstance(_resolved(resources, "/Font"), DictionaryObject):
            resources[NameObject("/Font")] = DictionaryObject()
    except Exception:
        # Best effort. Failing here would replace a confusing error with a
        # different confusing error.
        pass


def _flatten_form(writer: Any) -> None:
    """Set every form widget read-only.

    Deliberately NOT a true flatten: the values stay form-field values
    rather than being merged into the page content stream. Calling this
    "flattened" without saying so would overstate it — a reader
    extracting page text still will not find the entries.
    """
    try:
        from pypdf.generic import NameObject, NumberObject
        for page in writer.pages:
            for annot in page.get("/Annots", None) or []:
                try:
                    obj = annot.get_object()
                except Exception:
                    continue
                if obj.get("/Subtype") == "/Widget":
                    # Bit 1 (value 1) = read-only.
                    flags = int(obj.get("/Ff", 0) or 0)
                    obj[NameObject("/Ff")] = NumberObject(flags | 1)
    except Exception as exc:
        raise OfficeError(f"could not flatten the form: {exc}") from exc


# ---------------------------------------------------------------------------
# PDF assembly: merging, splitting, writing
# ---------------------------------------------------------------------------

# Caps on the batch sizes. A merge of a thousand files or a split of a
# thousand pages is far more likely to be a wrong argument than a wanted
# operation, and both would be discovered as a full directory.
MAX_MERGE_INPUTS = 100
MAX_SPLIT_PARTS = 200


def _pdf_pages(path: Path) -> tuple[Any, int]:
    """Open a PDF and count its pages, naming the file on every failure.

    An encrypted document opens without complaint and only refuses when a
    page is touched, so the page count is read here rather than during the
    write: a file that cannot be used has to be a refusal while nothing
    has been produced yet, not a half-finished result.
    """
    try:
        reader = _pdf_reader(path)
    except OfficeError as exc:
        raise OfficeError(f"{path.name}: {exc}") from exc
    try:
        return reader, len(reader.pages)
    except Exception as exc:
        if getattr(reader, "is_encrypted", False):
            raise OfficeError(
                f"{path.name} is encrypted and cannot be read without its "
                "password"
            ) from exc
        raise OfficeError(f"{path.name} could not be read: {exc}") from exc


def _pdf_has_fields(reader: Any) -> bool:
    try:
        return bool(reader.get_fields())
    except Exception:
        return False


def merge_pdfs(inputs: list, *, output: Any) -> dict:
    """Combine several PDFs into one new file, in the order given.

    Every input is opened and its pages are counted before anything is
    written. A merge that fails halfway would leave a document that looks
    finished while missing pages nobody can name afterwards, so an
    unreadable or encrypted input refuses the whole operation and says
    which file it was.

    An existing output is never replaced: in administrative work the file
    already sitting under that name is the one somebody else may have
    already sent.
    """
    if not isinstance(inputs, (list, tuple)) or len(inputs) < 2:
        raise OfficeError(
            "merging needs at least two PDFs, listed in the order they "
            "should appear in the result")
    if len(inputs) > MAX_MERGE_INPUTS:
        raise OfficeError(
            f"{len(inputs)} inputs exceeds the cap of {MAX_MERGE_INPUTS}")

    out = _resolve(output, must_exist=False)
    if document_kind(out) != "pdf":
        raise OfficeError(f"cannot write {out.suffix!r} as a PDF — use .pdf")
    if out.exists():
        raise OfficeError(
            f"{out.name} already exists — nothing was merged. Write to a "
            "new name.")

    pypdf = _require("pdf")
    sources: list[tuple[Path, Any, int]] = []
    for item in inputs:
        src = _resolve(item)
        if document_kind(src) != "pdf":
            raise OfficeError(
                f"{src.name} is not a PDF — nothing was merged")
        reader, count = _pdf_pages(src)
        sources.append((src, reader, count))

    expected = sum(count for _, _, count in sources)
    try:
        writer = pypdf.PdfWriter()
        for _src, reader, _count in sources:
            try:
                writer.append(reader)
            except AttributeError:      # pragma: no cover - older pypdf
                # No writer-level append: the pages still carry the
                # content, only the outline entries are lost.
                for page in reader.pages:
                    writer.add_page(page)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("wb") as fh:
            writer.write(fh)
    except OfficeError:
        raise
    except Exception as exc:
        # A truncated file under the requested name reads as a finished
        # merge to everyone who finds it later.
        try:
            out.unlink()
        except OSError:
            pass
        raise OfficeError(f"could not merge the documents: {exc}") from exc

    # Count the pages of the file that was written, not the pages that
    # were meant to go into it.
    written = 0
    try:
        written = len(pypdf.PdfReader(str(out)).pages)
    except Exception:
        written = 0

    notes: list[str] = []
    if written != expected:
        notes.append(
            f"the written file holds {written} page(s), not the {expected} "
            "the inputs add up to")
    with_fields = [src.name for src, reader, _c in sources
                   if _pdf_has_fields(reader)]
    if with_fields:
        notes.append(
            "form fields are present in " + ", ".join(with_fields)
            + " — fields carrying the same name collapse into one shared "
            "value when the documents are combined. Open the result before "
            "handing it over.")
    return {
        "path": str(out),
        "inputs": [{"path": str(src), "pages": count}
                   for src, _r, count in sources],
        "pages": written,
        "expected_pages": expected,
        "verified": written == expected and written > 0,
        "notes": notes,
    }


def _split_label(group: list[int]) -> str:
    """Name the page selection a part covers, for the output file name."""
    if len(group) == 1:
        return f"p{group[0]}"
    if group == list(range(group[0], group[-1] + 1)):
        return f"p{group[0]}-{group[-1]}"
    return "p" + "_".join(str(n) for n in group)


def split_pdf(
    path: Any, *, output_dir: Any, pages: Optional[str] = None
) -> dict:
    """Write page ranges of a PDF into separate files.

    ``pages`` takes the syntax :func:`read_pdf` uses, and every
    comma-separated part becomes one document: ``"2-5,9"`` writes a
    four-page file and a one-page file. Without it every page is written
    on its own.

    Each part is reported with its own outcome. A name already taken
    fails that part and leaves the existing file untouched, because the
    aggregate "12 files written" is exactly what hides the one page that
    silently replaced last week's version.
    """
    src = _resolve(path)
    if document_kind(src) != "pdf":
        raise OfficeError(f"{src.name} is not a PDF")

    pypdf = _require("pdf")
    reader, total = _pdf_pages(src)

    spec = str(pages or "").strip()
    if spec:
        # One output per comma-separated part, so "2-5" stays one
        # document. The parts themselves go through the reader's page
        # parser — there is one page syntax in this module, not two.
        groups = [_parse_pages(part, total)
                  for part in spec.split(",") if part.strip()]
        if not groups:
            raise OfficeError(f"no page selection in {spec!r}")
    else:
        groups = [[n] for n in _parse_pages(None, total)]
    if len(groups) > MAX_SPLIT_PARTS:
        raise OfficeError(
            f"{len(groups)} parts exceeds the cap of {MAX_SPLIT_PARTS} — "
            "pass a page selection")

    out_dir = _resolve(output_dir, must_exist=False)
    if out_dir.exists() and not out_dir.is_dir():
        raise OfficeError(f"{out_dir} exists and is not a directory")
    out_dir.mkdir(parents=True, exist_ok=True)

    stem = _safe_name(src.stem)
    used: dict[str, list[int]] = {}
    results: list[dict] = []
    for group in groups:
        name = f"{stem}_{_split_label(group)}.pdf"[:150]
        if name in used:
            results.append({
                "output": name, "pages": group, "status": "failed",
                "detail": (f"the name is already taken by pages "
                           f"{used[name]} of this call"),
            })
            continue
        target = out_dir / name
        if target.exists():
            results.append({
                "output": name, "pages": group, "status": "failed",
                "detail": "a file of that name already exists; nothing was "
                          "overwritten",
            })
            continue
        try:
            writer = pypdf.PdfWriter()
            for number in group:
                writer.add_page(reader.pages[number - 1])
            with target.open("wb") as fh:
                writer.write(fh)
        except Exception as exc:
            try:
                target.unlink()
            except OSError:
                pass
            results.append({"output": name, "pages": group,
                            "status": "failed", "detail": str(exc)})
            continue
        used[name] = group

        written = 0
        try:
            written = len(pypdf.PdfReader(str(target)).pages)
        except Exception:
            written = 0
        if written == len(group):
            results.append({"output": name, "pages": group, "status": "ok"})
        else:
            results.append({
                "output": name, "pages": group, "status": "incomplete",
                "detail": (f"the written file holds {written} page(s), not "
                           f"the {len(group)} requested"),
            })

    counts = {
        "ok": sum(1 for r in results if r["status"] == "ok"),
        "incomplete": sum(1 for r in results if r["status"] == "incomplete"),
        "failed": sum(1 for r in results if r["status"] == "failed"),
    }
    notes: list[str] = []
    if counts["failed"]:
        notes.append(f"{counts['failed']} part(s) produced no file.")
    if counts["incomplete"]:
        notes.append(
            f"{counts['incomplete']} file(s) do not hold the pages that were "
            "asked for. They are not ready to hand over.")
    if _pdf_has_fields(reader):
        notes.append(
            "the source is a form: the extracted pages keep the field "
            "widgets but not the document-level form dictionary, so filled "
            "values can render as empty. Check one result before handing "
            "it over.")
    return {
        "source": str(src),
        "output_dir": str(out_dir),
        "source_pages": total,
        "files": results,
        "counts": counts,
        "verified": counts["ok"] == len(groups),
        "notes": notes,
    }


def _pdf_markup(value: Any) -> str:
    """Escape a value for the paragraph markup of the PDF writer.

    Paragraph text is parsed as markup, so an unescaped ``&`` in a company
    name or a ``<`` in a note aborts the build or swallows the rest of the
    line. Line breaks are kept as breaks rather than collapsing into one
    run-on paragraph.
    """
    from xml.sax.saxutils import escape
    return escape(_fmt(value)).replace("\r\n", "\n").replace("\n", "<br/>")


_PROBE_RE = re.compile(r"[^\W\d_]{4,}", re.UNICODE)


def _text_probe(text: str) -> str:
    """The longest word of *text*, used to check the file reads back.

    A whole line cannot be searched for: the writer wraps lines wherever
    they run out of width. A single word is not hyphenated, so it survives
    extraction intact — and a word carrying an umlaut proves the encoding
    survived as well.
    """
    words = _PROBE_RE.findall(str(text or ""))
    if not words:
        return ""
    return max(words, key=len)[:24]


def _squash(text: str) -> str:
    return re.sub(r"\s+", "", str(text or ""))


def create_pdf(path: Any, blocks: list[dict]) -> dict:
    """Write a new PDF from the content blocks :func:`create_docx` takes.

    Each block is one of ``{"heading": str, "level": int}``,
    ``{"paragraph": str}``, ``{"table": [[...], ...], "header_row": bool}``
    or ``{"page_break": true}``.

    The result is read back before it is reported: the page count comes
    from the written file and the first block's text has to be findable in
    it, so a document that was produced but cannot be read does not come
    back as a success.
    """
    p = _resolve(path, must_exist=False)
    if document_kind(p) != "pdf":
        raise OfficeError(f"cannot write {p.suffix!r} as a PDF — use .pdf")
    if p.exists():
        raise OfficeError(
            f"{p.name} already exists — write to a new name; an existing "
            "document is never replaced")
    if not isinstance(blocks, list) or not blocks:
        raise OfficeError("blocks must be a non-empty list")

    _require("pdf_write")
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
        from reportlab.platypus import (
            PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
        )
    except Exception as exc:
        raise OfficeError(
            f"the 'reportlab' package is installed but incomplete: {exc}"
        ) from exc

    margin = 20 * mm
    doc = SimpleDocTemplate(
        str(p), pagesize=A4, title=p.stem,
        leftMargin=margin, rightMargin=margin,
        topMargin=margin, bottomMargin=margin)
    styles = getSampleStyleSheet()
    body = styles["BodyText"]
    cell_style = ParagraphStyle(
        "office_cell", parent=body, fontSize=9, leading=11)

    story: list[Any] = []
    counts = {"headings": 0, "paragraphs": 0, "tables": 0}
    probe = ""
    wide_tables = 0
    for index, block in enumerate(blocks, start=1):
        if not isinstance(block, dict):
            raise OfficeError(f"block {index} is not an object: {block!r}")
        if "heading" in block:
            level = max(0, min(int(block.get("level", 1) or 1), 6))
            style = (styles["Title"] if level == 0
                     else styles[f"Heading{level}"])
            text = str(block["heading"])
            story.append(Paragraph(_pdf_markup(text), style))
            counts["headings"] += 1
            probe = probe or _text_probe(text)
        elif "paragraph" in block:
            runs = _runs(block["paragraph"])
            if isinstance(block["paragraph"], (str, int, float)):
                for run in runs:
                    run["bold"] = bool(block.get("bold"))
                    run["italic"] = bool(block.get("italic"))
                    run["colour"] = parse_colour(
                        block.get("color", block.get("colour")))
            story.append(Paragraph(_runs_markup(runs), body))
            story.append(Spacer(1, 4))
            counts["paragraphs"] += 1
            probe = probe or _text_probe(
                "".join(r["text"] for r in runs))
        elif "table" in block:
            rows = block["table"]
            if not isinstance(rows, list) or not rows or not all(
                    isinstance(r, (list, tuple)) for r in rows):
                raise OfficeError(
                    f"block {index}: table must be a non-empty list of lists")
            width = max(len(r) for r in rows)
            # Equal column widths over the printable width: a table sized
            # from its content overflows the page silently, and a table
            # running off the paper is not a document anyone can file.
            col_width = doc.width / width
            header_row = bool(block.get("header_row"))

            def _cell(value, is_header):
                cell_runs = _runs(value)
                if is_header:
                    for run in cell_runs:
                        run["bold"] = True
                return Paragraph(_runs_markup(cell_runs), cell_style)

            data = [
                [_cell(r[c] if c < len(r) else "", header_row and i == 0)
                 for c in range(width)]
                for i, r in enumerate(rows)
            ]
            table = Table(data, colWidths=[col_width] * width,
                          repeatRows=1 if header_row else 0)
            commands = [
                ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
            if header_row:
                commands.append(
                    ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke))
            table.setStyle(TableStyle(commands))
            story.append(table)
            story.append(Spacer(1, 6))
            counts["tables"] += 1
            if width > 8:
                wide_tables += 1
            for row in rows:
                probe = probe or _text_probe(" ".join(_fmt(v) for v in row))
                if probe:
                    break
        elif block.get("page_break"):
            story.append(PageBreak())
        else:
            raise OfficeError(
                f"block {index}: expected one of heading / paragraph / "
                f"table / page_break, got keys {sorted(block)}")

    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        doc.build(story)
    except Exception as exc:
        try:
            p.unlink()
        except OSError:
            pass
        raise OfficeError(f"could not write the PDF: {exc}") from exc

    notes: list[str] = []
    pages = 0
    extracted = ""
    try:
        pypdf = _require("pdf")
        written = pypdf.PdfReader(str(p))
        pages = len(written.pages)
        extracted = " ".join(
            (page.extract_text() or "") for page in written.pages)
    except Exception as exc:
        notes.append(f"the written file could not be read back: {exc}")

    found = True
    if probe and pages:
        found = _squash(probe) in _squash(extracted)
    verified = pages > 0 and found
    if not pages:
        notes.append("the written file did not open as a PDF")
    elif not found:
        notes.append(
            f"the word {probe!r} was written but was not found when the text "
            "of the result was read back — check the document before "
            "handing it over")
    if wide_tables:
        notes.append(
            f"{wide_tables} table(s) have more than 8 columns; they were "
            "fitted into equal-width columns and may be hard to read on A4.")
    return {
        "path": str(p),
        "blocks": len(blocks),
        "pages": pages,
        "verified": verified,
        "notes": notes,
        **counts,
    }


# ---------------------------------------------------------------------------
# Word
# ---------------------------------------------------------------------------

_PLACEHOLDER_RE = re.compile(r"\{\{\s*([^{}]{1,80}?)\s*\}\}")


def _iter_paragraphs(container: Any) -> Any:
    """Every paragraph of a document part, tables and nesting included.

    Templates put placeholders in table cells and in the header and
    footer (the letterhead, the file reference, the page footer) at
    least as often as in the body. A filler that only walks
    ``doc.paragraphs`` silently leaves those in place, and the result
    looks complete until someone prints it.
    """
    for para in getattr(container, "paragraphs", []) or []:
        yield para
    for table in getattr(container, "tables", []) or []:
        for row in table.rows:
            for cell in row.cells:
                yield from _iter_paragraphs(cell)


def _document_parts(document: Any) -> Any:
    """The body plus every header/footer variant a section can carry."""
    yield document
    for section in getattr(document, "sections", []) or []:
        for attr in (
            "header", "footer",
            "first_page_header", "first_page_footer",
            "even_page_header", "even_page_footer",
        ):
            part = getattr(section, attr, None)
            if part is not None:
                yield part


def _splice_runs(runs: list, start: int, end: int, replacement: str) -> None:
    """Replace the character range [start, end) spread across *runs*.

    The replacement takes the formatting of the run the placeholder
    starts in; every other run keeps its own text and formatting. The
    naive alternative — writing the whole paragraph into its first run —
    would flatten bold, spacing and font changes across the line.
    """
    pos = 0
    first = True
    for run in runs:
        begin, stop = pos, pos + len(run.text)
        pos = stop
        if stop <= start or begin >= end:
            continue
        prefix = run.text[: max(0, start - begin)]
        suffix = run.text[max(0, end - begin):] if stop > end else ""
        if first:
            run.text = prefix + replacement + suffix
            first = False
        else:
            run.text = prefix + suffix


def _fill_paragraph(paragraph: Any, values: dict, unfilled: set) -> int:
    """Substitute placeholders in one paragraph. Returns the count.

    Works on the paragraph's joined text rather than run by run: Word
    splits a placeholder across runs whenever it feels like it (spell
    check, an edit, a formatting change), so ``{{anrede}}`` routinely
    arrives as ``{{`` + ``an`` + ``rede`` + ``}}``. Per-run replacement
    finds nothing and reports success.
    """
    filled = 0
    offset = 0
    for _ in range(500):                      # bound: pathological templates
        runs = list(paragraph.runs)
        if not runs:
            return filled
        text = "".join(r.text for r in runs)
        match = _PLACEHOLDER_RE.search(text, offset)
        if match is None:
            return filled
        name = match.group(1).strip()
        if name not in values:
            unfilled.add(name)
            offset = match.end()
            continue
        replacement = "" if values[name] is None else str(values[name])
        _splice_runs(runs, match.start(), match.end(), replacement)
        offset = match.start() + len(replacement)
        filled += 1
    return filled


def docx_placeholders(path: Any) -> dict:
    """List the ``{{name}}`` placeholders a template contains.

    Reported per location (body / table / header / footer) so a
    placeholder that only appears in the letterhead is visible before
    the fill, not after.
    """
    p = _resolve(path)
    if document_kind(p) != "word":
        raise OfficeError(f"{p.name} is not a .docx file")
    docx = _require("word")
    try:
        doc = docx.Document(str(p))
    except Exception as exc:
        raise OfficeError(f"could not open document: {exc}") from exc

    counts: dict[str, int] = {}
    for part in _document_parts(doc):
        for para in _iter_paragraphs(part):
            for match in _PLACEHOLDER_RE.finditer(para.text):
                name = match.group(1).strip()
                counts[name] = counts.get(name, 0) + 1

    notes: list[str] = []
    if not counts:
        notes.append(
            "no {{placeholder}} markers in this document. Word mail-merge "
            "fields (MERGEFIELD) are not text placeholders and are not "
            "handled — check what the template actually uses.")
    return {
        "path": str(p),
        "placeholders": [
            {"name": n, "occurrences": c} for n, c in sorted(counts.items())
        ],
        "notes": notes,
    }


def fill_docx_template(
    path: Any, values: dict, *, output: Any, strict: bool = True
) -> dict:
    """Substitute ``{{name}}`` placeholders and write a new document.

    The template is never written to, so it stays reusable. With
    *strict* a value whose placeholder does not exist is an error rather
    than a silent no-op — the same contract as the PDF form filler,
    because a typo would otherwise read as a completed document.

    Placeholders left without a value are reported and stay in the text
    verbatim: a half-filled letter that still shows ``{{betrag}}`` is
    recoverable, one that silently dropped it is not.
    """
    src = _resolve(path)
    if document_kind(src) != "word":
        raise OfficeError(f"{src.name} is not a .docx file")
    if not isinstance(values, dict) or not values:
        raise OfficeError("values must be a non-empty object of name -> value")
    out = _resolve(output, must_exist=False)
    if out.resolve() == src.resolve():
        raise OfficeError(
            "output must differ from the template — filling it in place "
            "destroys the blank original")

    known = {
        entry["name"] for entry in docx_placeholders(src)["placeholders"]
    }
    unknown = sorted(str(k) for k in values if str(k) not in known)
    if unknown and strict:
        raise OfficeError(
            f"no such placeholder(s): {', '.join(unknown)}. The template "
            f"contains: {', '.join(sorted(known)) or '(none)'}"
        )

    docx = _require("word")
    try:
        doc = docx.Document(str(src))
    except Exception as exc:
        raise OfficeError(f"could not open template: {exc}") from exc

    normalised = {str(k): v for k, v in values.items()}
    unfilled: set[str] = set()
    filled = 0
    for part in _document_parts(doc):
        for para in _iter_paragraphs(part):
            filled += _fill_paragraph(para, normalised, unfilled)

    try:
        out.parent.mkdir(parents=True, exist_ok=True)
        doc.save(str(out))
    except Exception as exc:
        raise OfficeError(f"could not save the document: {exc}") from exc

    # Verify against the written file rather than reporting the intent.
    remaining = {
        entry["name"]
        for entry in docx_placeholders(out)["placeholders"]
    }
    notes: list[str] = []
    if unknown:
        notes.append(
            "ignored, no such placeholder in the template: "
            + ", ".join(unknown))
    still_open = sorted(remaining)
    if still_open:
        notes.append(
            "these placeholders had no value and are still in the text: "
            + ", ".join(still_open)
            + ". The document is not ready to hand over as it stands.")
    leaked = sorted(n for n in normalised if n in remaining)
    if leaked:
        notes.append(
            "a value was given but the placeholder survived in the output: "
            + ", ".join(leaked) + " — treat this fill as failed.")
    return {
        "template": str(src),
        "path": str(out),
        "filled": filled,
        "unfilled": still_open,
        "complete": not still_open,
        "notes": notes,
    }


def create_docx(
    path: Any, blocks: list[dict], *, overwrite: bool = False
) -> dict:
    """Write a new .docx from a list of content blocks.

    Each block is one of ``{"heading": str, "level": int}``,
    ``{"paragraph": str}``, ``{"table": [[...], ...], "header_row":
    bool}`` or ``{"page_break": true}``.
    """
    p = _resolve(path, must_exist=False)
    if document_kind(p) != "word":
        raise OfficeError(
            f"cannot write {p.suffix!r} as a Word document — use .docx")
    if p.exists() and not overwrite:
        raise OfficeError(
            f"{p.name} already exists — pass overwrite=True to replace it")
    if not isinstance(blocks, list) or not blocks:
        raise OfficeError("blocks must be a non-empty list")

    docx = _require("word")
    doc = docx.Document()
    counts = {"headings": 0, "paragraphs": 0, "tables": 0}
    for index, block in enumerate(blocks, start=1):
        if not isinstance(block, dict):
            raise OfficeError(f"block {index} is not an object: {block!r}")
        if "heading" in block:
            level = max(0, min(int(block.get("level", 1) or 1), 9))
            doc.add_heading(str(block["heading"]), level=level)
            counts["headings"] += 1
        elif "paragraph" in block:
            para = doc.add_paragraph()
            runs = _runs(block["paragraph"])
            # Block-level bold/colour apply to a plain string; a list of
            # runs carries its own, so the label can stay plain while the
            # figure beside it is red.
            if isinstance(block["paragraph"], (str, int, float)):
                for run in runs:
                    run["bold"] = bool(block.get("bold"))
                    run["italic"] = bool(block.get("italic"))
                    run["colour"] = parse_colour(
                        block.get("color", block.get("colour")))
            _apply_runs(para, runs)
            counts["paragraphs"] += 1
        elif "table" in block:
            rows = block["table"]
            if not isinstance(rows, list) or not rows or not all(
                    isinstance(r, (list, tuple)) for r in rows):
                raise OfficeError(
                    f"block {index}: table must be a non-empty list of lists")
            width = max(len(r) for r in rows)
            table = doc.add_table(rows=len(rows), cols=width)
            if block.get("header_row"):
                try:
                    table.style = "Table Grid"
                except Exception:
                    pass
            for r_idx, row in enumerate(rows):
                for c_idx in range(width):
                    value = row[c_idx] if c_idx < len(row) else ""
                    cell = table.cell(r_idx, c_idx)
                    # A fresh cell already holds one empty run. Writing
                    # past it leaves the text in run 2, where anything
                    # reading the first run finds nothing.
                    target = cell.paragraphs[0]
                    for stale in list(target.runs):
                        stale._element.getparent().remove(stale._element)
                    cell_runs = _runs(value)
                    if r_idx == 0 and block.get("header_row"):
                        for run in cell_runs:
                            run["bold"] = True
                    _apply_runs(target, cell_runs)
            counts["tables"] += 1
        elif block.get("page_break"):
            doc.add_page_break()
        else:
            raise OfficeError(
                f"block {index}: expected one of heading / paragraph / "
                f"table / page_break, got keys {sorted(block)}")

    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        doc.save(str(p))
    except Exception as exc:
        raise OfficeError(f"could not save the document: {exc}") from exc
    return {"path": str(p), "blocks": len(blocks), **counts}


# Colours a caller may name. Deliberately a short list plus hex: an
# office document needs "this figure is red", not a palette, and a name
# that silently resolves to something else is worse than a refusal.
_COLOURS: dict[str, tuple[int, int, int]] = {
    "black": (0, 0, 0), "red": (192, 0, 0), "green": (0, 128, 0),
    "blue": (0, 0, 192), "orange": (224, 112, 0), "grey": (110, 110, 110),
    "gray": (110, 110, 110), "white": (255, 255, 255),
}

_HEX_RE = re.compile(r"^#?([0-9a-fA-F]{6})$")


def parse_colour(value: Any) -> Optional[tuple]:
    """A colour name or ``#RRGGBB`` as an RGB triple, else None."""
    if value is None or value == "":
        return None
    text = str(value).strip().lower()
    if text in _COLOURS:
        return _COLOURS[text]
    match = _HEX_RE.match(text)
    if match:
        digits = match.group(1)
        return tuple(int(digits[i:i + 2], 16) for i in (0, 2, 4))
    raise OfficeError(
        f"unknown colour {value!r}. Use one of "
        + ", ".join(sorted(_COLOURS)) + ", or #RRGGBB.")


def _runs(value: Any) -> list[dict]:
    """Normalise text into runs of ``{text, bold, italic, colour}``.

    A plain string is one unformatted run, so every existing caller keeps
    working. A list lets one paragraph carry an emphasised part without
    turning the whole line bold — which is what a report actually needs:
    the label plain, the figure red.
    """
    if value is None:
        return [{"text": "", "bold": False, "italic": False, "colour": None}]
    if isinstance(value, (str, int, float)):
        return [{"text": _fmt(value), "bold": False, "italic": False,
                 "colour": None}]
    if isinstance(value, dict):
        value = [value]
    if not isinstance(value, list):
        raise OfficeError(f"expected text or a list of runs, got {value!r}")
    out = []
    for part in value:
        if isinstance(part, (str, int, float)):
            out.append({"text": _fmt(part), "bold": False, "italic": False,
                        "colour": None})
            continue
        if not isinstance(part, dict):
            raise OfficeError(f"a run must be text or an object, got {part!r}")
        out.append({
            "text": _fmt(part.get("text", "")),
            "bold": bool(part.get("bold")),
            "italic": bool(part.get("italic")),
            "colour": parse_colour(part.get("color", part.get("colour"))),
        })
    return out or [{"text": "", "bold": False, "italic": False,
                    "colour": None}]


def _apply_runs(paragraph: Any, runs: list[dict]) -> None:
    """Write runs into a python-docx paragraph, keeping its font."""
    from docx.shared import RGBColor

    for run in runs:
        written = paragraph.add_run(run["text"])
        if run["bold"]:
            written.bold = True
        if run["italic"]:
            written.italic = True
        if run["colour"]:
            written.font.color.rgb = RGBColor(*run["colour"])


def _runs_markup(runs: list[dict]) -> str:
    """Runs as the inline markup reportlab understands, escaped."""
    from xml.sax.saxutils import escape

    out = []
    for run in runs:
        text = escape(run["text"])
        if run["colour"]:
            text = ('<font color="#%02x%02x%02x">%s</font>'
                    % (*run["colour"], text))
        if run["bold"]:
            text = f"<b>{text}</b>"
        if run["italic"]:
            text = f"<i>{text}</i>"
        out.append(text)
    return "".join(out)


def read_docx(path: Any, *, max_chars: int = MAX_TEXT_CHARS) -> dict:
    """Extract paragraphs and tables from a .docx file."""
    p = _resolve(path)
    _refuse_unreadable_format(p)
    if document_kind(p) == "opendocument_text":
        return read_odt(p, max_chars=max_chars)
    if document_kind(p) != "word":
        raise OfficeError(
            f"{p.name} is not a .docx file. The legacy .doc format is not "
            "supported — convert it first.")
    docx = _require("word")
    try:
        doc = docx.Document(str(p))
    except Exception as exc:
        raise OfficeError(f"could not open document: {exc}") from exc

    parts = [para.text for para in doc.paragraphs if para.text.strip()]
    tables = []
    for n, table in enumerate(doc.tables, start=1):
        rows = [[cell.text.strip() for cell in row.cells] for row in table.rows]
        tables.append(f"--- table {n} ---\n" + render_grid(rows))
    body = "\n".join(parts)
    if tables:
        body += "\n\n" + "\n\n".join(tables)
    notes: list[str] = []
    if len(body) > max_chars:
        body = body[:max_chars]
        notes.append(f"output truncated at {max_chars} characters")
    if not body.strip():
        notes.append("the document holds no extractable text")
    return {
        "path": str(p),
        "kind": "word",
        "paragraphs": len(parts),
        "tables": len(doc.tables),
        "text": body,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# OpenDocument text (.odt)
# ---------------------------------------------------------------------------

def _odt_cell_text(cell: Any) -> str:
    from odf import teletype
    try:
        return teletype.extractText(cell).strip()
    except Exception:
        return ""


def _odt_content(container: Any) -> tuple[list[str], list[list[list]]]:
    """``(paragraphs, tables)`` of an .odt body, in document order.

    Walks the tree rather than collecting paragraphs by type: in a
    letter the recipient block sits in a table and the enclosures sit in
    a list, and a reader that only takes top-level ``text:p`` elements
    returns a document with the addressee missing — which looks like a
    complete letter to anyone who did not have the original.
    """
    from odf.namespaces import TABLENS, TEXTNS
    from odf import teletype

    paragraphs: list[str] = []
    tables: list[list[list]] = []

    def walk(node: Any) -> None:
        for child in getattr(node, "childNodes", None) or []:
            qname = getattr(child, "qname", None)
            if not qname:
                continue
            namespace, name = qname
            if namespace == TEXTNS and name in ("p", "h"):
                try:
                    text = teletype.extractText(child).strip()
                except Exception:
                    text = ""
                if text:
                    paragraphs.append(text)
            elif namespace == TABLENS and name == "table":
                def _cells(row: Any) -> list:
                    return [
                        _odt_cell_text(cell)
                        for cell in (getattr(row, "childNodes", None) or [])
                        if (getattr(cell, "qname", None)
                            or ("", ""))[0] == TABLENS
                    ]

                rows = [
                    _cells(row)
                    for row in (getattr(child, "childNodes", None) or [])
                    if getattr(row, "qname", None) == (TABLENS, "table-row")
                ]
                if rows:
                    tables.append(rows)
            else:
                walk(child)

    walk(container)
    return paragraphs, tables


def read_odt(path: Any, *, max_chars: int = MAX_TEXT_CHARS) -> dict:
    """Extract paragraphs and tables from an OpenDocument text document."""
    p = _resolve(path)
    if document_kind(p) != "opendocument_text":
        raise OfficeError(f"{p.name} is not an .odt file")
    document = _odf_load(p)
    mimetype = str(getattr(document, "mimetype", "") or "")
    if mimetype and not mimetype.endswith(("text", "text-template")):
        raise OfficeError(
            f"{p.name} is named like a text document but its OpenDocument "
            f"type is '{mimetype}'. Read it under the name of what it is.")

    body = getattr(document, "text", None)
    if body is None:
        raise OfficeError(f"{p.name} holds no text body")
    paragraphs, tables = _odt_content(body)

    parts = list(paragraphs)
    rendered = [f"--- table {n} ---\n" + render_grid(rows)
                for n, rows in enumerate(tables, start=1)]
    text = "\n".join(parts)
    if rendered:
        text += ("\n\n" if text else "") + "\n\n".join(rendered)

    notes: list[str] = []
    if len(text) > max_chars:
        text = text[:max_chars]
        notes.append(f"output truncated at {max_chars} characters")
    if not text.strip():
        notes.append("the document holds no extractable text")
    notes.append(
        "read-only: an .odt cannot be written here. To produce a filled "
        "document, convert the template to .docx and fill that.")
    return {
        "path": str(p),
        "kind": "opendocument_text",
        "format": "odt",
        "writable": False,
        "paragraphs": len(parts),
        "tables": len(tables),
        "text": text,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

def read_document(path: Any, **kwargs: Any) -> dict:
    """Read any supported document, dispatching on its suffix."""
    p = _resolve(path)
    _refuse_unreadable_format(p)
    kind = document_kind(p)
    if kind == "opendocument_text":
        if kwargs.get("fields"):
            raise OfficeError(
                f"{p.name} is an OpenDocument text document; there is no "
                "template filler for it. Convert it to .docx — "
                + _LIBREOFFICE_HINT.format(target="docx", name=p.name)
                + " — and use the .docx as the template.")
        return read_odt(p)
    if kind in ("spreadsheet", "csv", "opendocument_sheet"):
        return read_sheet(
            p,
            sheet=kwargs.get("sheet"),
            max_rows=int(kwargs.get("max_rows") or DEFAULT_MAX_ROWS),
            max_cols=int(kwargs.get("max_cols") or DEFAULT_MAX_COLS),
            start_row=int(kwargs.get("start_row") or 1),
        )
    if kind == "pdf":
        if kwargs.get("fields"):
            return pdf_form_fields(p)
        return read_pdf(p, pages=kwargs.get("pages"),
                        ocr=bool(kwargs.get("ocr")))
    if kind == "word":
        # Same switch as for a PDF form: "show me what can be filled in"
        # is one question, whatever the format underneath.
        if kwargs.get("fields"):
            return docx_placeholders(p)
        return read_docx(p)
    raise OfficeError(
        f"{p.name}: no document reader for {p.suffix!r}. Supported: "
        + ", ".join(sorted(DOCUMENT_SUFFIXES))
    )
