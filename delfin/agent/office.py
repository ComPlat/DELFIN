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

Each of those is handled once, here, instead of in whatever code the
model happens to write that turn. The functions are pure in the sense
that matters for testing: they take paths and values, they touch only
the file named, and they raise only :class:`OfficeError`.

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

DOCUMENT_SUFFIXES = (
    SPREADSHEET_SUFFIXES | CSV_SUFFIXES | PDF_SUFFIXES | WORD_SUFFIXES
)

# Suffix -> (import name, distribution name). The distribution name is
# what the user has to install, and it differs from the import name for
# python-docx, which is the classic reason a "just pip install docx"
# hint sends someone to an unrelated package.
_BACKENDS: dict[str, tuple[str, str]] = {
    "spreadsheet": ("openpyxl", "openpyxl"),
    "pdf": ("pypdf", "pypdf"),
    "word": ("docx", "python-docx"),
    "csv": ("csv", ""),  # standard library
}


class OfficeError(RuntimeError):
    """A document operation failed for a reason worth reporting verbatim."""


# ---------------------------------------------------------------------------
# Backends and file kinds
# ---------------------------------------------------------------------------

def document_kind(path: Any) -> Optional[str]:
    """``"spreadsheet"`` / ``"csv"`` / ``"pdf"`` / ``"word"``, or None."""
    suffix = Path(str(path)).suffix.lower()
    if suffix in SPREADSHEET_SUFFIXES:
        return "spreadsheet"
    if suffix in CSV_SUFFIXES:
        return "csv"
    if suffix in PDF_SUFFIXES:
        return "pdf"
    if suffix in WORD_SUFFIXES:
        return "word"
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
        for kind in ("spreadsheet", "pdf", "word")
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
# Spreadsheets — reading
# ---------------------------------------------------------------------------

def _read_csv(
    path: Path, *, max_rows: int, max_cols: int, start_row: int
) -> dict:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        if path.suffix.lower() != ".tsv":
            try:
                delimiter = csv.Sniffer().sniff(
                    sample, delimiters=",;\t|").delimiter
            except csv.Error:
                delimiter = ","
        all_rows = list(csv.reader(fh, delimiter=delimiter))

    total_rows = len(all_rows)
    total_cols = max((len(r) for r in all_rows), default=0)
    begin = max(0, start_row - 1)
    window = [r[:max_cols] for r in all_rows[begin:begin + max_rows]]
    notes = []
    if delimiter != ",":
        notes.append(f"delimiter detected: {delimiter!r}")
    if total_rows > begin + len(window):
        notes.append(
            f"showing rows {begin + 1}-{begin + len(window)} of {total_rows} "
            f"— pass start_row to page further")
    if total_cols > max_cols:
        notes.append(f"showing {max_cols} of {total_cols} columns")
    return {
        "path": str(path),
        "kind": "csv",
        "sheets": [],
        "sheet": "",
        "rows": total_rows,
        "columns": total_cols,
        "grid": render_grid(window, first_row=begin + 1),
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
    kind = document_kind(p)
    max_rows = max(1, min(int(max_rows), 2000))
    max_cols = max(1, min(int(max_cols), 200))
    start_row = max(1, int(start_row))

    if kind == "csv":
        return _read_csv(
            p, max_rows=max_rows, max_cols=max_cols, start_row=start_row)
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

        return {
            "path": str(p),
            "kind": "spreadsheet",
            "sheets": names,
            "sheet": ws.title,
            "rows": total_rows,
            "columns": total_cols,
            "grid": render_grid(window, first_row=start_row),
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


def edit_sheet(
    path: Any,
    *,
    edits: Optional[list[dict]] = None,
    append_rows: Optional[list[list[Any]]] = None,
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
    if document_kind(p) != "spreadsheet":
        raise OfficeError(
            f"{p.name} is not a spreadsheet (suffix {p.suffix!r})")
    if not edits and not append_rows:
        raise OfficeError("nothing to do: pass edits and/or append_rows")

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

        return {
            "path": str(p),
            "sheet": ws.title,
            "applied": applied,
            "appended_rows": appended,
            "first_appended_row": first_appended,
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
        with p.open("w", encoding="utf-8", newline="") as fh:
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
# PDF
# ---------------------------------------------------------------------------

def _pdf_reader(path: Path) -> Any:
    pypdf = _require("pdf")
    try:
        return pypdf.PdfReader(str(path))
    except Exception as exc:
        raise OfficeError(f"could not open PDF: {exc}") from exc


def _has_xfa(reader: Any) -> bool:
    """True if the document is an XFA form.

    XFA (LiveCycle) forms keep their data in an XML stream, not in the
    AcroForm field dictionaries. Writing AcroForm values on such a file
    reports success and changes nothing the user will see, so this has
    to be detected before a fill is claimed to have worked.
    """
    try:
        root = reader.trailer["/Root"]
        acro = root.get("/AcroForm")
        if acro is None:
            return False
        try:
            acro = acro.get_object()
        except Exception:
            pass
        return "/XFA" in acro
    except Exception:
        return False


def read_pdf(
    path: Any,
    *,
    pages: Optional[str] = None,
    max_chars: int = MAX_TEXT_CHARS,
) -> dict:
    """Extract text from a PDF, page by page.

    ``pages`` accepts ``"3"``, ``"2-5"`` or ``"1,4,7"``; omitted means
    the whole document up to *max_chars*.
    """
    p = _resolve(path)
    if document_kind(p) != "pdf":
        raise OfficeError(f"{p.name} is not a PDF")
    reader = _pdf_reader(p)
    total = len(reader.pages)
    wanted = _parse_pages(pages, total)

    chunks: list[str] = []
    used = 0
    truncated = False
    empty_pages = 0
    for number in wanted:
        try:
            text = reader.pages[number - 1].extract_text() or ""
        except Exception as exc:
            text = f"(page {number}: extraction failed: {exc})"
        if not text.strip():
            empty_pages += 1
        block = f"--- page {number} ---\n{text.strip()}"
        if used + len(block) > max_chars:
            chunks.append(block[: max(0, max_chars - used)])
            truncated = True
            break
        chunks.append(block)
        used += len(block)

    notes: list[str] = []
    if truncated:
        notes.append(
            f"output truncated at {max_chars} characters — pass pages to "
            "read a specific range")
    if empty_pages == len(wanted) and wanted:
        notes.append(
            "no extractable text on any requested page. This is most likely "
            "a scanned document: the pages are images, and reading it needs "
            "OCR, which is not part of this tool.")
    elif empty_pages:
        notes.append(f"{empty_pages} of {len(wanted)} pages held no text")
    if _has_xfa(reader):
        notes.append(
            "this is an XFA form — its field data lives in an XML stream, "
            "not in the page text")

    return {
        "path": str(p),
        "kind": "pdf",
        "pages": total,
        "pages_read": wanted,
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

    fields = []
    for name, spec in raw.items():
        entry: dict[str, Any] = {
            "name": str(name),
            "type": _FIELD_TYPES.get(str(spec.get("/FT", "")), str(
                spec.get("/FT", "") or "unknown")),
            "value": _fmt(spec.get("/V")),
        }
        states = spec.get("/_States_")
        if states:
            entry["states"] = [str(s) for s in states]
        fields.append(entry)

    notes: list[str] = []
    if not fields:
        notes.append(
            "this PDF has no form fields. Values cannot be filled in; "
            "producing a completed document means generating a new PDF "
            "or overlaying text.")
    if _has_xfa(reader):
        notes.append(
            "XFA form: writing AcroForm field values will report success "
            "but will not change what a viewer displays. Confirm the "
            "result before handing this file over.")
    return {
        "path": str(p),
        "fields": fields,
        "field_count": len(fields),
        "xfa": _has_xfa(reader),
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
# Word
# ---------------------------------------------------------------------------

def read_docx(path: Any, *, max_chars: int = MAX_TEXT_CHARS) -> dict:
    """Extract paragraphs and tables from a .docx file."""
    p = _resolve(path)
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
# Dispatch
# ---------------------------------------------------------------------------

def read_document(path: Any, **kwargs: Any) -> dict:
    """Read any supported document, dispatching on its suffix."""
    p = _resolve(path)
    kind = document_kind(p)
    if kind in ("spreadsheet", "csv"):
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
        return read_pdf(p, pages=kwargs.get("pages"))
    if kind == "word":
        return read_docx(p)
    raise OfficeError(
        f"{p.name}: no document reader for {p.suffix!r}. Supported: "
        + ", ".join(sorted(DOCUMENT_SUFFIXES))
    )
