#!/usr/bin/env python3
"""Acceptance for the twelve everyday office tasks (package O1).

Every check reads the RESULT FILE the attempt wrote and compares it with
the expectation the setup module derives from the same fixture data
(``setup/o1_office_fixtures.py`` is imported, never copied: a graded
number typed twice drifts). Nothing here reads the answer text -- a
correct sentence over an untouched workbook is exactly the failure this
suite was written to catch.

WHAT A CHECK MUST REJECT, and why each rejection is one of these:

* **untouched** -- the result file does not exist. The measured failure
  (three sessions, then a correct script the shell could not run) leaves
  no file at all.
* **half-done** -- the file exists and is right about the easy half. The
  deliberately planted trap of each case is the half a plausible script
  gets wrong: 49 rows silently dropped, 18 sheet names cut to a shared
  prefix, 7 cells that name two people given one sheet between them, 60
  rows in another date format, 6 categories missing from a left join.
  ``_HALF_DONE`` below states, per case, the shape of that wrong answer,
  and ``tests/test_the_office_benchmark_grades_the_file.py`` builds each
  one and asserts this script REJECTS it.
* **wrong type** -- "1.234,56" as text where a number is expected.
  Office work is arithmetic afterwards; text that looks like a number is
  the failure that survives review.

A check that passes an empty workspace would measure nothing, so every
check starts by requiring the file.

argv[1] is the workspace. ``DELFIN_BENCH_TASK_ID`` names the task being
scored (the runner sets it); without one, ALL twelve are graded, which is
what the gate test wants.
"""
from __future__ import annotations

import importlib.util
import os
import re
import sys
from pathlib import Path

ROOT = "o1"
_SETUP_PATH = Path(__file__).resolve().parents[1] / "setup" / "o1_office_fixtures.py"

#: Per case, the plausible wrong answer that must NOT be accepted. Used by
#: the test that proves each check has teeth; stated here so the reason a
#: check is shaped the way it is lives beside the check.
_HALF_DONE = {
    "blaettern": "one sheet per RAW value (196 sheets, no merging of the "
                 "case/space variants, no shared sheet for the 7 two-person "
                 "cells, 49 rows dropped or put on a sheet of their own)",
    "zusammenfuehren": "an inner join: the orders without a creditor "
                       "disappear from the result instead of being listed",
    "dubletten": "only the second occurrence of a group marked, or the "
                 "nine same-Beleg-other-amount rows marked too",
    "monatssummen": "the 60 ISO-dated rows ignored (the totals still look "
                    "like totals) or the amounts read as text",
    "serienbrief": "one letter written for all rows, or the placeholder "
                   "left in the body",
    "formular": "one filled copy instead of one per application, or the "
                "fields filled from the wrong row",
    "pdf_tabelle": "the repeated header and the page footer taken into the "
                   "table, or two pages read as one",
    "bereinigen": "raw strings carried over (titles kept, ue/oe/ae not "
                  "resolved) or the split of 'Nachname, Vorname' missing",
    "pivot": "counts per category but only for the months that occur, or "
             "sums instead of counts",
    "csv_deutsch": "the amount or the date left as text",
    "regional_rollup": "a left join from the north: 'Auftragsfertigung' "
                       "and the two south-only categories missing",
    "rechnungs_dubletten": "a case- and space-sensitive comparison: the 14 "
                           "duplicate pairs read as 28 invoices",
}


def _setup_mod():
    spec = importlib.util.spec_from_file_location("o1_setup", _SETUP_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Reading what the attempt wrote
# ---------------------------------------------------------------------------


def _open_book(path: Path):
    from openpyxl import load_workbook
    return load_workbook(str(path), data_only=False)


def _sheet_rows(ws) -> list[list]:
    return [list(row) for row in ws.iter_rows(values_only=True)]


def _find(workspace: Path, *names: str) -> Path | None:
    """The first of ``names`` that exists, at the workspace root or one deep.

    One deep because a person who says "in den Ordner o1/blaettern" does
    not also say where the output belongs, and an attempt that writes it
    beside the input is right. A name is still required: the check is
    "this file exists", not "something exists".
    """
    for name in names:
        direct = workspace / name
        if direct.is_file():
            return direct
        for found in sorted(workspace.glob(f"*/{name}")):
            if found.is_file():
                return found
        for found in sorted(workspace.glob(f"*/*/{name}")):
            if found.is_file():
                return found
    return None


def _for_case(workspace: Path, case: str) -> Path:
    """The case's own directory, so a loose search cannot find a sibling's file."""
    return workspace / ROOT / case


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def _is_number(value) -> bool:
    """True for a real numeric cell.

    ``bool`` is not a number here: openpyxl reads a ticked box as True,
    and True == 1 would let a column of ticks pass as a column of counts.
    A string is never a number, however numeric it looks -- that
    distinction is the whole point of the "correct number types" case.
    """
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _as_number(value) -> float | None:
    """A cell as a float, accepting only what a spreadsheet would sum.

    Text is NOT accepted. A cell holding ``"7,5"`` is a string and the
    answer is "wrong type", which is reported rather than smoothed over.
    """
    if _is_number(value):
        return float(value)
    return None


def _norm(value) -> str:
    """Whitespace-collapsed text, for comparing a label."""
    return " ".join(str(value if value is not None else "").split())


# ---------------------------------------------------------------------------
# blaettern -- one sheet per value of a column
# ---------------------------------------------------------------------------


def _check_blaettern(case: Path, mod) -> str:
    path = _find(case, "blaettern.xlsx", "aufgeteilt.xlsx", "bestellungen_je_bearbeiter.xlsx")
    if path is None:
        return ("no result workbook: expected one of blaettern.xlsx / "
                "aufgeteilt.xlsx / bestellungen_je_bearbeiter.xlsx in "
                f"{case}")
    expected = mod.expected_sheets()
    wb = _open_book(path)
    got = list(wb.sheetnames)
    missing = sorted(set(expected) - set(got))
    extra = sorted(set(got) - set(expected))
    if missing or extra:
        detail = []
        if missing:
            detail.append(f"{len(missing)} sheet(s) missing, e.g. {missing[:4]}")
        if extra:
            detail.append(f"{len(extra)} sheet(s) that should not exist, "
                          f"e.g. {extra[:4]}")
        return ("the sheet map is not the expected one -- " + "; ".join(detail)
                + f" (expected {len(expected)} sheets, found {len(got)})")
    # The rows, not only the names: a workbook with the right 196 titles and
    # every sheet empty satisfies a name-only check.
    problems = []
    for title, belege in expected.items():
        rows = _sheet_rows(wb[title])
        if not rows:
            problems.append(f"{title}: empty sheet")
            continue
        header = [_norm(c) for c in rows[0]]
        if "Beleg" not in header:
            problems.append(f"{title}: header does not name the Beleg column")
            continue
        col = header.index("Beleg")
        got_belege = sorted(_norm(r[col]) for r in rows[1:]
                            if col < len(r) and _norm(r[col]))
        want = sorted(belege)
        if got_belege != want:
            problems.append(
                f"{title}: {len(got_belege)} row(s), expected {len(want)}")
        if len(problems) >= 4:
            break
    if problems:
        return ("the sheets exist but their contents are not the rows that "
                "belong to them: " + "; ".join(problems))
    total = sum(len(v) for v in expected.values())
    wb.close()
    return True, f"{len(expected)} sheets, {total} data rows distributed as expected"


# ---------------------------------------------------------------------------
# zusammenfuehren -- merge two tables on a key
# ---------------------------------------------------------------------------


def _check_zusammenfuehren(case: Path, mod) -> str:
    path = _find(case, "ergebnis_join.xlsx")
    if path is None:
        return f"no result workbook ergebnis_join.xlsx in {case}"
    orders, creditors, _ = mod.join_rows()
    master = {str(row[0]): (str(row[1]), str(row[2])) for row in creditors}
    unmatched = [str(r[1]) for r in orders[1:] if str(r[1]) not in master]
    wb = _open_book(path)
    rows = _sheet_rows(wb.active)
    if not rows:
        return "ergebnis_join.xlsx has no rows"
    header = [_norm(c) for c in rows[0]]
    for need in ("Beleg", "Name", "Ort", "Betrag"):
        if need not in header:
            return f"ergebnis_join.xlsx: header has no {need!r} column"
    data = [r for r in rows[1:] if any(_norm(c) for c in r)]
    # An inner join drops the 6 orders with no master record. Those are
    # the rows a plausible answer loses -- and the total still adds up.
    if len(data) < len(orders) - 1:
        return (f"ergebnis_join.xlsx has {len(data)} rows; the orders table "
                f"has {len(orders) - 1}, and only {len(unmatched)} of them "
                "have no master record, so rows were dropped")
    beleg_col = header.index("Beleg")
    got = {_norm(r[beleg_col]) for r in data}
    missing = [str(r[0]) for r in orders[1:] if str(r[0]) not in got]
    if len(missing) > len(unmatched):
        return (f"{len(missing)} orders missing from the result, but only "
                f"{len(unmatched)} lack a master record, e.g. {missing[:4]}")
    list_path = _find(case, "nicht_zugeordnet.txt", "unmatched.txt")
    if list_path is None:
        return "no nicht_zugeordnet.txt listing the orders without a master record"
    listed = {line.strip() for line in _text(list_path).splitlines()
              if line.strip()}
    want = {str(r[0]) for r in orders if str(r[1]) not in master}
    # Row 0 of orders is the header; its Beleg is the word "Beleg".
    want.discard("Beleg")
    if listed != want:
        return (f"nicht_zugeordnet.txt lists {len(listed)} order(s); expected "
                f"exactly the {len(want)} without a master record "
                f"(off by {sorted(listed ^ want)[:4]})")
    wb.close()
    return True, (f"{len(data)} rows kept, {len(want)} unmatched orders listed")


# ---------------------------------------------------------------------------
# dubletten -- find and mark duplicates
# ---------------------------------------------------------------------------


def _check_dubletten(case: Path, mod, *, name: str = "ergebnis_dubletten.xlsx",
                     column: str = "Dublette", yes: str = "ja") -> str:
    path = _find(case, name)
    if path is None:
        return f"no result workbook {name} in {case}"
    wb = _open_book(path)
    rows = _sheet_rows(wb.active)
    if not rows:
        return f"{name} has no rows"
    header = [_norm(c) for c in rows[0]]
    if column not in header:
        return f"{name}: no {column!r} column"
    col = header.index(column)
    marked = {idx for idx, row in enumerate(rows[1:], start=1)
              if col < len(row) and _norm(row[col]).lower() == yes.lower()}
    want = mod.duplicate_rows()[1]
    want_rows = {n - 1 for n in want}      # 1-based Excel row -> 0-based index
    if marked != want_rows:
        over = sorted(marked - want_rows)
        under = sorted(want_rows - marked)
        parts = []
        if under:
            parts.append(f"{len(under)} duplicate row(s) not marked, "
                         f"e.g. {[i + 2 for i in under[:4]]}")
        if over:
            parts.append(f"{len(over)} row(s) marked that are not duplicates, "
                         f"e.g. {[i + 2 for i in over[:4]]}")
        return "; ".join(parts)
    wb.close()
    return True, f"{len(marked)} rows marked as duplicates, exactly the expected set"


def _check_rechnungs_dubletten(case: Path, mod) -> str:
    path = _find(case, "register_dubletten.xlsx")
    if path is None:
        return f"no result workbook register_dubletten.xlsx in {case}"
    wb = _open_book(path)
    rows = _sheet_rows(wb.active)
    if not rows:
        return "register_dubletten.xlsx has no rows"
    header = [_norm(c) for c in rows[0]]
    if "Dublette" not in header:
        return "register_dubletten.xlsx: no 'Dublette' column"
    col = header.index("Dublette")
    marked = {idx for idx, row in enumerate(rows[1:], start=1)
              if col < len(row) and _norm(row[col]).lower() in ("ja", "yes")}
    want_rows = {n - 1 for n in mod.register_rows()[1]}
    if marked != want_rows:
        under = sorted(want_rows - marked)
        over = sorted(marked - want_rows)
        parts = []
        if under:
            parts.append(f"{len(under)} duplicate row(s) not marked "
                         f"(a case- or space-sensitive comparison misses these)")
        if over:
            parts.append(f"{len(over)} row(s) marked that are unique")
        return "; ".join(parts)
    wb.close()
    return True, f"{len(marked)} rows marked, exactly the expected set"


# ---------------------------------------------------------------------------
# monatssummen -- monthly totals from a date column with decimal commas
# ---------------------------------------------------------------------------


def _check_monatssummen(case: Path, mod) -> str:
    path = _find(case, "monatssummen.xlsx")
    if path is None:
        return f"no result workbook monatssummen.xlsx in {case}"
    expected = mod.expected_months()
    wb = _open_book(path)
    rows = _sheet_rows(wb.active)
    if not rows:
        return "monatssummen.xlsx has no rows"
    header = [_norm(c) for c in rows[0]]
    if not header or header[0].lower() not in ("monat", "month"):
        return f"monatssummen.xlsx: first column is {header[:1]}, not 'Monat'"
    got: dict[str, float] = {}
    for row in rows[1:]:
        if not row or not _norm(row[0]):
            continue
        value = _as_number(row[1] if len(row) > 1 else None)
        if value is None:
            return (f"monatssummen.xlsx: the sum for {_norm(row[0])!r} is not a "
                    f"number ({row[1]!r}) -- a total written as text cannot be "
                    "added up")
        got[_norm(row[0]).lower()] = value
    missing = [m for m in expected if m.lower() not in got]
    if missing:
        return f"monatssummen.xlsx: no row for {missing}"
    wrong = [(m, got[m.lower()], expected[m][0]) for m in expected
             if abs(got[m.lower()] - expected[m][0]) > 0.005]
    if wrong:
        month, was, want = wrong[0]
        return (f"{month}: {was:.2f} EUR, expected {want:.2f} EUR "
                f"({len(wrong)} month(s) off) -- a total over fewer rows "
                "still looks like a total")
    wb.close()
    n_rows = sum(n for _sum, n in expected.values())
    return True, f"{len(expected)} months, sums to the cent over all {n_rows} rows"


# ---------------------------------------------------------------------------
# serienbrief -- a mail merge into .docx letters
# ---------------------------------------------------------------------------


def _docx_paragraphs(path: Path) -> str:
    from docx import Document
    doc = Document(str(path))
    parts = [p.text for p in doc.paragraphs]
    for table in doc.tables:
        for row in table.rows:
            parts.extend(cell.text for cell in row.cells)
    return "\n".join(parts)


def _check_serienbrief(case: Path, mod) -> str:
    target = None
    for candidate in (case / "briefe", case / "letters", case):
        if candidate.is_dir() and any(candidate.glob("*.docx")):
            target = candidate
            break
    if target is None:
        return f"no letters: no .docx under {case}/briefe"
    letters = sorted(p for p in target.glob("*.docx"))
    rows = mod.mailmerge_rows()[1:]
    if len(letters) < len(rows):
        return (f"{len(letters)} letter(s) for {len(rows)} rows -- a merge "
                "that stops early leaves the rest of the list unposted")
    by_name: dict[str, Path] = {}
    for path in letters:
        for row in rows:
            surname, given = str(row[2]), str(row[1])
            if surname.lower() in path.name.lower() and \
                    given.lower() in path.name.lower():
                by_name[surname] = path
                break
    missing = [str(r[2]) for r in rows if str(r[2]) not in by_name]
    if missing:
        return (f"no letter named for {len(missing)} row(s), e.g. {missing[:3]} "
                "(briefe/<Nachname>_<Vorname>.docx)")
    problems = []
    for row in rows:
        path = by_name[str(row[2])]
        body = _docx_paragraphs(path)
        if "{{" in body:
            problems.append(f"{path.name}: a placeholder is still in the body")
            continue
        anrede, given, surname, betrag, ort = (str(v) for v in row)
        for need, label in ((surname, "Nachname"), (given, "Vorname"),
                            (betrag, "Betrag"), (ort, "Ort")):
            if need not in body:
                problems.append(f"{path.name}: {label} {need!r} not in the body")
        if anrede not in body:
            problems.append(f"{path.name}: Anrede {anrede!r} not in the body")
        if len(problems) >= 4:
            break
    if problems:
        return "; ".join(problems)
    return True, f"{len(rows)} letters, every placeholder replaced with its row's values"


# ---------------------------------------------------------------------------
# formular -- fill a PDF form from a table row
# ---------------------------------------------------------------------------


def _pdf_fields(path: Path) -> dict[str, str]:
    import pypdf
    fields = pypdf.PdfReader(str(path)).get_fields() or {}
    out = {}
    for name, spec in fields.items():
        value = spec.get("/V") if isinstance(spec, dict) else None
        out[str(name)] = "" if value is None else str(value)
    return out


def _check_formular(case: Path, mod) -> str:
    target = None
    for candidate in (case / "antraege_ausgefuellt", case / "ausgefuellt", case):
        if candidate.is_dir() and any(candidate.glob("*.pdf")):
            target = candidate
            break
    if target is None:
        return f"no filled forms: no .pdf under {case}/antraege_ausgefuellt"
    wanted = {str(row[0]): row for row in mod.application_rows()[1:]}
    forms = {}
    for path in sorted(target.glob("*.pdf")):
        number = path.stem.split("_")[0].split("-")[0]
        key = next((k for k in wanted if k in path.stem), None)
        if key is None:
            key = next((k for k in wanted
                        if k.replace("A-", "") in path.stem), None)
        if key is not None:
            forms[key] = path
        else:
            forms.setdefault(path.stem, path)
    missing = [k for k in wanted if k not in forms]
    if missing:
        return (f"no filled form for {missing} -- one filled copy was asked "
                "for per application")
    problems = []
    for number, row in wanted.items():
        fields = _pdf_fields(forms[number])
        want = {"Name": str(row[1]), "Kostenstelle": str(row[2]),
                "Betrag": str(row[3]), "Datum": str(row[4])}
        for field, value in want.items():
            got = _norm(fields.get(field, ""))
            if got != _norm(value):
                problems.append(f"{number}: field {field!r} is {got!r}, "
                                f"expected {_norm(value)!r}")
        if len(problems) >= 4:
            break
    if problems:
        return "; ".join(problems)
    return True, f"{len(wanted)} forms filled, every field from its own row"


# ---------------------------------------------------------------------------
# pdf_tabelle -- extract a table from a PDF into .xlsx
# ---------------------------------------------------------------------------


def _check_pdf_tabelle(case: Path, mod) -> str:
    path = _find(case, "stunden.xlsx")
    if path is None:
        return f"no result workbook stunden.xlsx in {case}"
    expected = mod.hours_rows()
    wb = _open_book(path)
    rows = [r for r in _sheet_rows(wb.active) if any(_norm(c) for c in r)]
    if not rows:
        return "stunden.xlsx has no rows"
    header = [_norm(c) for c in rows[0]]
    for need in ("Datum", "Gruppe", "Stunden", "Projekt"):
        if need not in header:
            return f"stunden.xlsx: header has no {need!r} column"
    body = rows[1:]
    # The two traps: the header repeated on page 2 becomes a body row, and
    # the footer becomes one too. Both are visible as a row that is not a
    # data row.
    footer = [r for r in body
              if any(re.match(r"(?i)^seite\s+\d+", _norm(c)) for c in r)]
    if footer:
        return (f"{len(footer)} page footer row(s) taken into the table, e.g. "
                f"{[_norm(c) for c in footer[0]][:2]} -- the footer is not data")
    repeats = [r for r in body if _norm(r[0]) == "Datum"]
    if repeats:
        return (f"the page-2 header row was taken as a data row "
                f"({len(repeats)} repeat(s))")
    if len(body) != len(expected):
        return (f"stunden.xlsx has {len(body)} data rows; the PDF holds "
                f"{len(expected)} over two pages")
    problems = []
    d_col, g_col = header.index("Datum"), header.index("Gruppe")
    s_col, p_col = header.index("Stunden"), header.index("Projekt")
    for got, want in zip(body, expected):
        want_datum, want_gruppe, want_hours, want_projekt = want
        if _norm(got[d_col]) != want_datum:
            problems.append(f"row {_norm(got[d_col])!r} != {want_datum!r}")
            continue
        hours = _as_number(got[s_col])
        if hours is None:
            problems.append(f"{want_datum}: Stunden {got[s_col]!r} is not a "
                            "number")
        elif abs(hours - float(want_hours.replace(",", "."))) > 0.001:
            problems.append(f"{want_datum}: Stunden {hours} != "
                            f"{want_hours}")
        if _norm(got[g_col]) != want_gruppe:
            problems.append(f"{want_datum}: Gruppe {_norm(got[g_col])!r}")
        if _norm(got[p_col]) != want_projekt:
            problems.append(f"{want_datum}: Projekt {_norm(got[p_col])!r}")
        if len(problems) >= 4:
            break
    if problems:
        return "; ".join(problems)
    wb.close()
    return True, f"{len(expected)} rows over two pages, footer and repeat excluded"


# ---------------------------------------------------------------------------
# bereinigen -- clean a column
# ---------------------------------------------------------------------------


def _check_bereinigen(case: Path, mod) -> str:
    path = _find(case, "personen_bereinigt.xlsx")
    if path is None:
        return f"no result workbook personen_bereinigt.xlsx in {case}"
    expected = mod.expected_clean()
    wb = _open_book(path)
    rows = _sheet_rows(wb.active)
    if not rows:
        return "personen_bereinigt.xlsx has no rows"
    header = [_norm(c).lower() for c in rows[0]]
    if "nachname" not in header or "vorname" not in header:
        return ("personen_bereinigt.xlsx: need the columns Nachname and "
                f"Vorname, found {header[:4]}")
    s_col, g_col = header.index("nachname"), header.index("vorname")
    body = rows[1:]
    if len(body) != len(expected):
        return (f"{len(body)} output rows for {len(expected)} input rows -- "
                "one output row per input row was asked for, order unchanged")
    problems = []
    for pos, (got, want) in enumerate(zip(body, expected), start=2):
        got_s = _norm(got[s_col]) if s_col < len(got) else ""
        got_g = _norm(got[g_col]) if g_col < len(got) else ""
        want_s, want_g = want
        if (got_s, got_g) != (_norm(want_s), _norm(want_g)):
            problems.append(f"row {pos}: ({got_s!r}, {got_g!r}) != "
                            f"({want_s!r}, {want_g!r})")
        if len(problems) >= 4:
            break
    if problems:
        return "; ".join(problems)
    wb.close()
    return True, f"{len(expected)} rows cleaned: titles dropped, umlauts restored, split"


# ---------------------------------------------------------------------------
# pivot -- counts per category
# ---------------------------------------------------------------------------


def _check_pivot(case: Path, mod) -> str:
    path = _find(case, "pivot.xlsx", "haeufigkeiten.xlsx")
    if path is None:
        return f"no result workbook pivot.xlsx in {case}"
    expected = mod.expected_pivot()
    months = list(mod.MONTHS)
    wb = _open_book(path)
    rows = _sheet_rows(wb.active)
    if not rows:
        return "pivot.xlsx has no rows"
    header = [_norm(c) for c in rows[0]]
    month_cols = {}
    for month in months:
        idx = next((i for i, h in enumerate(header)
                    if h.lower() == month.lower()), None)
        if idx is None:
            return f"pivot.xlsx: no column for {month}"
        month_cols[month] = idx
    got: dict[str, dict[str, int]] = {}
    for row in rows[1:]:
        if not row or not _norm(row[0]):
            continue
        label = _norm(row[0])
        got[label] = {}
        for month, idx in month_cols.items():
            value = row[idx] if idx < len(row) else None
            if not _is_number(value):
                return (f"pivot.xlsx {label}/{month}: {value!r} is not a "
                        "number -- a count written as text is not a count")
            got[label][month] = int(value)
    missing = [c for c in expected if c not in got]
    if missing:
        return f"pivot.xlsx: no row for {missing}"
    wrong = []
    for category, per_month in expected.items():
        for month, want in per_month.items():
            have = got[category].get(month)
            if have != want:
                wrong.append(f"{category}/{month}: {have} != {want}")
    if wrong:
        return ("counts are wrong: " + "; ".join(wrong[:4])
                + f" ({len(wrong)} cell(s)) -- sums instead of counts, or a "
                "month missing from the header")
    # A missing combination is a zero, and a zero must be written.
    empty = [(c, m) for c in expected for m in months
             if m not in expected[c] and got[c].get(m) not in (0, None)]
    if empty:
        return (f"{len(empty)} combination(s) with no sales carry a count "
                f"instead of 0, e.g. {empty[:3]}")
    wb.close()
    return True, f"{len(expected)} categories x {len(months)} months, counts as numbers"


# ---------------------------------------------------------------------------
# csv_deutsch -- a German-locale CSV into .xlsx with correct types
# ---------------------------------------------------------------------------


def _check_csv_deutsch(case: Path, mod) -> str:
    path = _find(case, "umsatz.xlsx")
    if path is None:
        return f"no result workbook umsatz.xlsx in {case}"
    expected = mod.expected_csv_deutsch()
    wb = _open_book(path)
    rows = _sheet_rows(wb.active)
    if not rows:
        return "umsatz.xlsx has no rows"
    header = [_norm(c) for c in rows[0]]
    for need in ("Kostenstelle", "Bezeichnung", "Datum", "Betrag"):
        if need not in header:
            return f"umsatz.xlsx: header has no {need!r} column"
    k_col = header.index("Kostenstelle")
    b_col = header.index("Bezeichnung")
    d_col = header.index("Datum")
    a_col = header.index("Betrag")
    body = [r for r in rows[1:] if any(_norm(c) for c in r)]
    if len(body) != len(expected):
        return f"umsatz.xlsx has {len(body)} rows, expected {len(expected)}"
    problems = []
    for got, want in zip(body, expected):
        want_k, want_b, want_d, want_a = want
        # The leading zeros are the trap: as a number, 004711 becomes 4711
        # and the row stops matching every other table that uses the key.
        if _norm(got[k_col]) != want_k:
            problems.append(f"Kostenstelle {got[k_col]!r} != {want_k!r} "
                            "(leading zeros lost)")
        if _norm(got[b_col]) != want_b:
            problems.append(f"Bezeichnung {got[b_col]!r} != {want_b!r}")
        amount = _as_number(got[a_col])
        if amount is None:
            problems.append(f"Betrag {got[a_col]!r} is not a number")
        elif abs(amount - want_a) > 0.005:
            problems.append(f"Betrag {amount} != {want_a}")
        datum = got[d_col]
        iso = None
        if hasattr(datum, "strftime"):
            iso = datum.strftime("%Y-%m-%d")
        else:
            text = _norm(datum)
            if re.fullmatch(r"\d{4}-\d{2}-\d{2}", text):
                iso = text
        if iso != want_d:
            problems.append(f"Datum {datum!r} is not the date {want_d} "
                            "(left as text?)")
        if len(problems) >= 4:
            break
    if problems:
        return "; ".join(problems)
    wb.close()
    return True, (f"{len(expected)} rows: amount a number, date a date, "
                 "Kostenstelle text with its leading zeros")


# ---------------------------------------------------------------------------
# regional_rollup -- (en) two regional sheets into one summary
# ---------------------------------------------------------------------------


def _check_regional_rollup(case: Path, mod) -> str:
    path = _find(case, "zusammenfassung.xlsx")
    if path is None:
        return f"no result workbook zusammenfassung.xlsx in {case}"
    expected = mod.expected_rollup()
    wb = _open_book(path)
    rows = [r for r in _sheet_rows(wb.active) if any(_norm(c) for c in r)]
    if not rows:
        return "zusammenfassung.xlsx has no rows"
    header = [_norm(c) for c in rows[0]]
    for need in ("Kategorie", "Umsatz"):
        if need not in header:
            return f"zusammenfassung.xlsx: header has no {need!r} column"
    c_col, a_col = header.index("Kategorie"), header.index("Umsatz")
    got: dict[str, float] = {}
    for row in rows[1:]:
        label = _norm(row[c_col])
        if not label:
            continue
        value = _as_number(row[a_col])
        if value is None:
            return (f"{label}: Umsatz {row[a_col]!r} is not a number -- "
                    "revenue written as text cannot be summed")
        got[label] = value
    missing = sorted(set(expected) - set(got))
    if missing:
        return (f"{len(missing)} categor(y/ies) from one of the two files "
                f"missing, e.g. {missing[:4]} -- an inner or one-sided join "
                "loses the categories that occur on one side only")
    wrong = [(c, got[c], expected[c]) for c in expected
             if abs(got[c] - expected[c]) > 0.005]
    if wrong:
        category, was, want = wrong[0]
        return (f"{category}: {was:.2f} != {want:.2f} "
                f"({len(wrong)} categor(y/ies) off)")
    extra = sorted(set(got) - set(expected))
    if extra:
        return f"categories that occur in neither file: {extra[:4]}"
    wb.close()
    return True, f"{len(expected)} categories, summed over both files as numbers"


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

_CHECKS = {
    "blaettern": _check_blaettern,
    "zusammenfuehren": _check_zusammenfuehren,
    "dubletten": _check_dubletten,
    "monatssummen": _check_monatssummen,
    "serienbrief": _check_serienbrief,
    "formular": _check_formular,
    "pdf_tabelle": _check_pdf_tabelle,
    "bereinigen": _check_bereinigen,
    "pivot": _check_pivot,
    "csv_deutsch": _check_csv_deutsch,
    "regional_rollup": _check_regional_rollup,
    "rechnungs_dubletten": _check_rechnungs_dubletten,
}

#: Case -> the task ids that grade it. Both modes exist for every case, so
#: a score in office mode and a score in solo mode ask the same question.
TASK_CASE = {}
for _case in _CHECKS:
    TASK_CASE[f"o1_office_{_case}"] = _case
    TASK_CASE[f"o1_solo_{_case}"] = _case


def check(workspace, task_id=None, **_ignored) -> tuple[bool, str]:
    """``(accepted, message)`` for ``workspace``.

    ``task_id=None`` grades ALL twelve cases (the gate test wants that); a
    single id grades only the case it names. An unknown id is refused
    rather than silently graded against everything -- a typo in a task
    file must not look like a passing task.
    """
    workspace = Path(workspace)
    mod = _setup_mod()
    if task_id is None:
        cases = list(_CHECKS)
    else:
        case = TASK_CASE.get(task_id)
        if case is None:
            return False, f"unknown task id {task_id!r}"
        cases = [case]
    done = []
    for case in cases:
        fn = _CHECKS[case]
        try:
            result = fn(_for_case(workspace, case), mod)
        except Exception as exc:                       # a broken file is a fail
            return False, f"{case}: cannot be graded ({type(exc).__name__}: {exc})"
        # The verdict is the TYPE, never a reading of the sentence. A check
        # returns ``(True, "<what it found>")`` when it accepted and a bare
        # string naming the defect when it did not. The first version
        # decided by matching the opening words of the message, which made a
        # refusal whose wording happened to start with "no " a pass -- the
        # assertion that cannot fail, one level up. A tuple is not a string,
        # so no sentence can be mistaken for a verdict.
        if isinstance(result, tuple) and len(result) == 2 and result[0] is True:
            done.append(f"{case}: {result[1]}")
        else:
            return False, f"{case}: {result}"
    return True, "; ".join(done)


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    tid = (sys.argv[2] if len(sys.argv) > 2
           else os.environ.get("DELFIN_BENCH_TASK_ID") or None)
    ok, msg = check(ws, tid)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
