#!/usr/bin/env python3
"""The reference solution for the twelve O1 office cases.

Not a fixture and not a grader: the answer a careful person's script gives,
written out. It exists for one reason -- ``accept/o1_office_accept.py`` has to
be shown to ACCEPT something, not only to reject. A check that rejects an
untouched workspace proves it has a floor; only a check that accepts a
correct result proves it has a ceiling, and a check with no ceiling rejects
everything and would silently make every task unpassable.

The companion test (``tests/test_the_office_benchmark_grades_the_file.py``)
runs this script in a scratch workspace, grades with the real acceptance
script, and asserts acceptance. It also runs it with ONE deliberate defect
per case (``--defect <case>``) and asserts rejection, so both sides are
measured: the reference is what "right" means, the defect is one plausible
wrong answer per case.

Usage::

    python3 o1_office_reference.py <workspace> [--defect CASE]

The workspace must already hold the fixtures
(``setup/o1_office_fixtures.py <workspace>``); this script writes its results
beside them, in ``<workspace>/o1/<case>/``, the same place the acceptance
script looks.
"""
from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from pathlib import Path

_SETUP_PATH = (Path(__file__).resolve().parents[1]
               / "setup" / "o1_office_fixtures.py")


def _fixtures():
    spec = importlib.util.spec_from_file_location("o1_setup", _SETUP_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# blaettern
# ---------------------------------------------------------------------------


def solve_blaettern(case: Path, fx, defect: str | None) -> str:
    """One sheet per collaborator, following HINWEISE.txt to the letter."""
    import openpyxl

    src = openpyxl.load_workbook(str(case / "bestellungen.xlsx"))
    sheet = src.active
    rows = [list(r) for r in sheet.iter_rows(values_only=True)]
    header, body = rows[0], rows[1:]
    col = header.index("Bearbeiter")

    if defect == "raw":
        # The plausible wrong answer: one sheet per raw cell value. Nobody
        # merges the case/space variants, the two-person cells get a sheet
        # of their own, and the 49 empty cells get one called "" (openpyxl
        # refuses it, so they are dropped).
        groups: dict[str, list[list]] = {}
        for row in body:
            value = str(row[col] or "")
            if not value.strip():
                continue
            groups.setdefault(value, []).append(row)
        wanted = [(v[:31].translate(str.maketrans('/\\?*[]:', '-------'))[:31], v)
                  for v in groups]
    else:
        # A sheet per PERSON, merged by case and spacing, split on the stated
        # separators, and one sheet for the rows with nobody in them.
        title_for: dict[str, str] = {}
        taken: dict[str, str] = {}
        groups: dict[str, list[list]] = {}
        for row in body:
            value = str(row[col] or "")
            if not value.strip():
                groups.setdefault("Ohne Bearbeiter", []).append(row)
                continue
            for person in fx._people_in(value):
                key = fx._merge_key(person)
                if key not in title_for:
                    base = fx._slug(person)
                    title, n = base, 2
                    while title in taken and taken[title] != key:
                        suffix = f"~{n}"
                        title = base[:fx.MAX_TITLE - len(suffix)] + suffix
                        n += 1
                    taken[title] = key
                    title_for[key] = title
                groups.setdefault(title_for[key], []).append(row)
        wanted = list(groups.items())

    out = openpyxl.Workbook()
    out.remove(out.active)
    for title, group in (wanted if defect != "raw" else
                         [(t, groups[v]) for t, v in wanted]):
        ws = out.create_sheet(title=title[:31])
        ws.append(header)
        for row in group:
            ws.append(row)
    out.save(str(case / "blaettern.xlsx"))
    out.close()
    src.close()
    return f"wrote {len(wanted)} sheets"


# ---------------------------------------------------------------------------
# zusammenfuehren
# ---------------------------------------------------------------------------


def solve_zusammenfuehren(case: Path, fx, defect: str | None) -> str:
    import openpyxl

    orders_wb = openpyxl.load_workbook(str(case / "bestellungen.xlsx"))
    master_wb = openpyxl.load_workbook(str(case / "kreditoren.xlsx"))
    orders = [list(r) for r in orders_wb.active.iter_rows(values_only=True)]
    master = [list(r) for r in master_wb.active.iter_rows(values_only=True)]
    lookup = {str(r[0]): (r[1], r[2]) for r in master[1:]}

    out = openpyxl.Workbook()
    ws = out.active
    ws.title = "Ergebnis"
    ws.append(["Beleg", "Name", "Ort", "Betrag"])
    unmatched = []
    for row in orders[1:]:
        key = str(row[1])
        if key not in lookup:
            unmatched.append(str(row[0]))
            if defect == "inner":
                # The plausible wrong answer: an inner join. The order with
                # no master record simply is not in the result, and the
                # total of the result still adds up.
                continue
            ws.append([row[0], "", "", row[2]])
        else:
            name, place = lookup[key]
            ws.append([row[0], name, place, row[2]])
    out.save(str(case / "ergebnis_join.xlsx"))
    out.close()
    orders_wb.close()
    master_wb.close()
    if defect != "inner":
        (case / "nicht_zugeordnet.txt").write_text(
            "\n".join(unmatched) + "\n", encoding="utf-8")
    return f"{len(orders) - 1} orders, {len(unmatched)} without a master record"


# ---------------------------------------------------------------------------
# dubletten
# ---------------------------------------------------------------------------


def solve_dubletten(case: Path, fx, defect: str | None) -> str:
    import openpyxl

    wb = openpyxl.load_workbook(str(case / "eingang.xlsx"))
    ws = wb.active
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    header = rows[0]
    seen: dict[tuple, list[int]] = {}
    for pos, row in enumerate(rows[1:], start=2):
        seen.setdefault((str(row[0]), str(row[1]), str(row[3])), []).append(pos)
    dup_positions = {p for positions in seen.values() if len(positions) > 1
                     for p in positions}
    if defect == "second-only":
        # The plausible wrong answer: keep the first copy, flag the rest.
        # Every group has the right number of rows and one row is "the
        # original", so the output looks reviewed.
        dup_positions = {p for positions in seen.values() if len(positions) > 1
                         for p in positions[1:]}
    ws.cell(row=1, column=len(header) + 1, value="Dublette")
    for pos in range(2, len(rows) + 1):
        ws.cell(row=pos, column=len(header) + 1,
                value="ja" if pos in dup_positions else None)
    wb.save(str(case / "ergebnis_dubletten.xlsx"))
    wb.close()
    return f"marked {len(dup_positions)} rows"


def solve_rechnungs_dubletten(case: Path, fx, defect: str | None) -> str:
    import openpyxl

    wb = openpyxl.load_workbook(str(case / "register.xlsx"))
    ws = wb.active
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    header = rows[0]

    def key(value):
        if defect == "strict":
            # The plausible wrong answer: an exact comparison. Two rows of
            # the same invoice typed by two people are two invoices.
            return str(value)
        return "".join(str(value).split()).upper()

    seen: dict[str, list[int]] = {}
    for pos, row in enumerate(rows[1:], start=2):
        seen.setdefault(key(row[0]), []).append(pos)
    dup = {p for positions in seen.values() if len(positions) > 1
           for p in positions}
    ws.cell(row=1, column=len(header) + 1, value="Dublette")
    for pos in range(2, len(rows) + 1):
        ws.cell(row=pos, column=len(header) + 1, value="ja" if pos in dup else None)
    wb.save(str(case / "register_dubletten.xlsx"))
    wb.close()
    return f"marked {len(dup)} rows"


# ---------------------------------------------------------------------------
# monatssummen
# ---------------------------------------------------------------------------


def _parse_de_amount(text: str) -> float:
    """'1.234,56' -> 1234.56, through DELFIN's own German parser.

    Not a second parser. ``delfin.agent.office.parse_number`` with
    ``DECIMAL_COMMA`` is the product's answer to this exact question, and
    a reference solution that hand-rolled its own would prove nothing
    about the product -- and would keep working after the product's
    parser broke. The assertion is the point: this line fails loudly if
    DELFIN can no longer read a German amount.
    """
    from delfin.agent.office import DECIMAL_COMMA, parse_number

    value = parse_number(text, DECIMAL_COMMA)
    if value is None:
        raise ValueError(f"not a German amount: {text!r}")
    return value


def _parse_date(text: str):
    """``(year, month, day)`` from either form the column holds.

    Both forms are in the fixture on purpose, so the convention cannot be
    guessed from the column: ISO rows are parsed as ISO and dotted rows
    day-first, stated per row. ``office.parse_date`` is the product's
    parser for both conventions, so this is again a use, not a copy.
    """
    from delfin.agent.office import DAY_FIRST, ISO_DATE, parse_date

    raw = str(text).strip()
    convention = ISO_DATE if re.fullmatch(r"\d{4}-\d{2}-\d{2}", raw) else DAY_FIRST
    iso = parse_date(raw, convention)
    if iso is None:
        raise ValueError(f"not a date in either form: {text!r}")
    year, month, day = (int(part) for part in iso.split("-"))
    return year, month, day


def solve_monatssummen(case: Path, fx, defect: str | None) -> str:
    import openpyxl

    rows = [line.split(";") for line in
            (case / "spenden.csv").read_text(encoding="utf-8").splitlines()[1:]
            if line.strip()]
    totals: dict[int, float] = {}
    for datum, _spender, amount in rows:
        if defect == "dotted-only" and "-" in datum:
            # The plausible wrong answer: one date format is understood and
            # the imported rows are silently skipped. Jan-June still come
            # out, and each monthly total still looks like a total.
            continue
        _year, month, _day = _parse_date(datum)
        totals[month] = round(totals.get(month, 0.0)
                              + _parse_de_amount(amount), 2)
    out = openpyxl.Workbook()
    ws = out.active
    ws.title = "Monatssummen"
    ws.append(["Monat", "Summe"])
    for idx, name in enumerate(fx.MONTHS, start=1):
        if defect == "text":
            # The plausible wrong answer: the German-formatted string
            # written straight into the cell. It reads correctly and is
            # not a number.
            ws.append([name, f"{totals.get(idx, 0.0):.2f}".replace(".", ",")])
        else:
            ws.append([name, round(totals.get(idx, 0.0), 2)])
    out.save(str(case / "monatssummen.xlsx"))
    out.close()
    return f"{len(totals)} months from {len(rows)} rows"


# ---------------------------------------------------------------------------
# serienbrief
# ---------------------------------------------------------------------------


def solve_serienbrief(case: Path, fx, defect: str | None) -> str:
    from docx import Document
    import openpyxl

    wb = openpyxl.load_workbook(str(case / "teilnehmer.xlsx"))
    rows = [list(r) for r in wb.active.iter_rows(values_only=True)]
    header, body = rows[0], rows[1:]
    wb.close()
    out_dir = case / "briefe"
    out_dir.mkdir(exist_ok=True)
    template = Document(str(case / "vorlage.docx"))

    written = []
    for row in body if defect != "first-only" else body[:1]:
        document = Document(str(case / "vorlage.docx"))
        values = dict(zip(header, row))
        for paragraph in document.paragraphs:
            text = paragraph.text
            for field, value in values.items():
                if defect == "keep-placeholder" and field == "Betrag":
                    continue
                text = text.replace("{{" + field + "}}", str(value))
            if paragraph.runs:
                paragraph.runs[0].text = text
                for extra in paragraph.runs[1:]:
                    extra.text = ""
        target = out_dir / f"{values['Nachname']}_{values['Vorname']}.docx"
        document.save(str(target))
        written.append(target.name)
    # Touch the template so an accidental edit shows up in a diff check.
    assert template.paragraphs, "template unreadable"
    return f"{len(written)} letters"


# ---------------------------------------------------------------------------
# formular
# ---------------------------------------------------------------------------


def solve_formular(case: Path, fx, defect: str | None) -> str:
    import pypdf

    rows = [line.split(";") for line in
            (case / "antraege.csv").read_text(encoding="utf-8").splitlines()]
    header, body = rows[0], rows[1:]
    out_dir = case / "antraege_ausgefuellt"
    out_dir.mkdir(exist_ok=True)
    written = []
    source = case / "antrag.pdf"
    for row in (body if defect != "one-only" else body[:1]):
        values = dict(zip(header, row))
        reader = pypdf.PdfReader(str(source))
        writer = pypdf.PdfWriter(clone_from=str(source))
        fields = {"Name": values["Name"], "Kostenstelle": values["Kostenstelle"],
                  "Betrag": values["Betrag"], "Datum": values["Datum"]}
        if defect == "wrong-row" and written:
            # The plausible wrong answer: the file is written per row but
            # filled with the FIRST row's values, so the count is right and
            # the content is not.
            fields = {"Name": body[0][1], "Kostenstelle": body[0][2],
                      "Betrag": body[0][3], "Datum": body[0][4]}
        writer.update_page_form_field_values(writer.pages[0], fields)
        target = out_dir / f"{values['Antragsnummer']}.pdf"
        with target.open("wb") as fh:
            writer.write(fh)
        reader.close()
        written.append(target.name)
    return f"{len(written)} forms"


# ---------------------------------------------------------------------------
# pdf_tabelle
# ---------------------------------------------------------------------------


def solve_pdf_tabelle(case: Path, fx, defect: str | None) -> str:
    import openpyxl
    import pypdf

    reader = pypdf.PdfReader(str(case / "stunden.pdf"))
    header = ["Datum", "Gruppe", "Stunden", "Projekt"]
    data: list[list] = []
    for page in reader.pages:
        text = page.extract_text(extraction_mode="layout") or ""
        for line in text.splitlines():
            if defect != "keep-footer" and re.match(r"(?i)^\s*Seite\s+\d+", line):
                continue
            cells = [c.strip() for c in re.split(r"\s{2,}", line.strip()) if c.strip()]
            if not cells:
                continue
            if cells[0] == "Datum" and header[:len(cells)] == cells[:len(header)]:
                if defect == "keep-footer":
                    # The plausible wrong answer: every line of every page
                    # is a row, so the page-2 header becomes data too.
                    data.append(cells[:4])
                continue
            if re.match(r"^\d{2}\.\d{2}\.\d{4}$", cells[0]) and len(cells) >= 4:
                datum, gruppe, stunden, projekt = cells[0], cells[1], cells[2], cells[3]
                if defect == "text-hours":
                    # The plausible wrong answer: the German-formatted hours
                    # carried over as text. It reads "7,5" exactly as the PDF
                    # does, and the column can no longer be added up.
                    hours = stunden
                else:
                    hours = float(stunden.replace(",", "."))
                data.append([datum, gruppe, hours, projekt])
    out = openpyxl.Workbook()
    ws = out.active
    ws.title = "Stunden"
    ws.append(header)
    for row in data:
        ws.append(row)
    out.save(str(case / "stunden.xlsx"))
    out.close()
    return f"{len(data)} rows"


# ---------------------------------------------------------------------------
# bereinigen
# ---------------------------------------------------------------------------


def solve_bereinigen(case: Path, fx, defect: str | None) -> str:
    import openpyxl

    wb = openpyxl.load_workbook(str(case / "personen.xlsx"))
    rows = [list(r) for r in wb.active.iter_rows(values_only=True)]
    wb.close()
    out = openpyxl.Workbook()
    ws = out.active
    ws.title = "Personal"
    ws.append(["Nachname", "Vorname"])
    for (cell,) in rows[1:]:
        text = " ".join(str(cell or "").split())
        if defect == "raw":
            # The plausible wrong answer: the cell split on the comma and
            # nothing else. Titles stay in the surname and the umlaut
            # digraphs stay as they were typed.
            surname, _, given = text.partition(",")
            ws.append([surname.strip(), given.strip()])
            continue
        for title in ("Dr. ", "Prof. "):
            if text.startswith(title):
                text = text[len(title):]
        if not text:
            ws.append(["", ""])
            continue
        surname, _, given = text.partition(",")
        ws.append([fx._to_umlaut(surname.strip()), fx._to_umlaut(given.strip())])
    out.save(str(case / "personen_bereinigt.xlsx"))
    out.close()
    return f"{len(rows) - 1} rows"


# ---------------------------------------------------------------------------
# pivot
# ---------------------------------------------------------------------------


def solve_pivot(case: Path, fx, defect: str | None) -> str:
    import openpyxl

    wb = openpyxl.load_workbook(str(case / "verkaeufe.xlsx"))
    rows = [list(r) for r in wb.active.iter_rows(values_only=True)]
    wb.close()
    months = list(fx.MONTHS)
    counts: dict[str, dict[str, int]] = {}
    sums: dict[str, dict[str, float]] = {}
    for category, month, amount in rows[1:]:
        counts.setdefault(str(category), {}).setdefault(str(month), 0)
        counts[str(category)][str(month)] += 1
        sums.setdefault(str(category), {}).setdefault(str(month), 0.0)
        sums[str(category)][str(month)] += float(amount)
    out = openpyxl.Workbook()
    ws = out.active
    ws.title = "Pivot"
    ws.append(["Kategorie"] + months)
    for category in fx.CATEGORIES:
        if defect == "sums":
            # The plausible wrong answer: the amount summed instead of the
            # rows counted. Every cell is a number and the table looks like
            # a pivot.
            ws.append([category] + [round(sums[category].get(m, 0.0), 2)
                                    for m in months])
        else:
            ws.append([category] + [counts.get(category, {}).get(m, 0)
                                    for m in months])
    out.save(str(case / "pivot.xlsx"))
    out.close()
    return f"{len(fx.CATEGORIES)} categories"


# ---------------------------------------------------------------------------
# csv_deutsch
# ---------------------------------------------------------------------------


def solve_csv_deutsch(case: Path, fx, defect: str | None) -> str:
    import datetime

    import openpyxl

    text = (case / "umsatz.csv").read_bytes().decode("cp1252")
    lines = [line.split(";") for line in text.splitlines() if line.strip()]
    header, body = lines[0], lines[1:]
    out = openpyxl.Workbook()
    ws = out.active
    ws.title = "Umsatz"
    ws.append(header)
    for number, label, datum, amount in body:
        day, month, year = (int(p) for p in datum.split("."))
        value = float(amount.replace(".", "").replace(",", "."))
        if defect == "text":
            # The plausible wrong answer: everything carried over as text.
            # It looks exactly like the CSV, and nothing can be summed.
            ws.append([number, label, datum, amount])
            continue
        key = number if defect == "number-key" else str(number)
        if defect == "number-key":
            ws.append([int(number), label, datetime.date(year, month, day), value])
        else:
            cell = ws.cell(row=ws.max_row + 1, column=1, value=key)
            cell.number_format = "@"
            ws.cell(row=ws.max_row, column=2, value=label)
            ws.cell(row=ws.max_row, column=3,
                    value=datetime.date(year, month, day))
            ws.cell(row=ws.max_row, column=4, value=value)
    out.save(str(case / "umsatz.xlsx"))
    out.close()
    return f"{len(body)} rows"


# ---------------------------------------------------------------------------
# regional_rollup
# ---------------------------------------------------------------------------


def solve_regional_rollup(case: Path, fx, defect: str | None) -> str:
    import openpyxl

    totals: dict[str, float] = {}
    for name in ("region_nord.xlsx", "region_sued.xlsx"):
        if defect == "north-only" and name == "region_sued.xlsx":
            # The plausible wrong answer: the north file treated as the
            # master and the south joined onto it, so the categories that
            # occur only in the south disappear.
            continue
        wb = openpyxl.load_workbook(str(case / name))
        rows = [list(r) for r in wb.active.iter_rows(values_only=True)]
        wb.close()
        for category, amount in rows[1:]:
            totals[str(category)] = round(
                totals.get(str(category), 0.0) + float(amount), 2)
    out = openpyxl.Workbook()
    ws = out.active
    ws.title = "Zusammenfassung"
    ws.append(["Kategorie", "Umsatz"])
    for category, amount in totals.items():
        ws.append([category, amount])
    out.save(str(case / "zusammenfassung.xlsx"))
    out.close()
    return f"{len(totals)} categories"


_SOLVERS = {
    "blaettern": solve_blaettern,
    "zusammenfuehren": solve_zusammenfuehren,
    "dubletten": solve_dubletten,
    "monatssummen": solve_monatssummen,
    "serienbrief": solve_serienbrief,
    "formular": solve_formular,
    "pdf_tabelle": solve_pdf_tabelle,
    "bereinigen": solve_bereinigen,
    "pivot": solve_pivot,
    "csv_deutsch": solve_csv_deutsch,
    "regional_rollup": solve_regional_rollup,
    "rechnungs_dubletten": solve_rechnungs_dubletten,
}


def solve(workspace: Path | str, defect: str | None = None,
          only: str | None = None) -> dict[str, str]:
    """Solve every case (or ``only``), with ``case:variant`` as a defect.

    ``defect`` is spelled ``<case>:<variant>`` -- one case's plausible wrong
    answer instead of the right one. The case is part of the string because
    every solver takes its variant by the same name it branches on, and the
    first version passed the CASE here instead: the variant then never
    matched anything, every "defect" run produced the correct answer, and
    the acceptance script accepted all seventeen of them. A defect mode that
    cannot produce a defect is not a test of the checker.

    EVERY case is solved on a defect run, not only the one named: the
    others must still be right, so a rejection can only be the named
    variant and not an accident of the workspace.
    """
    fx = _fixtures()
    ws = Path(workspace)
    target, _sep, variant = (defect or "").partition(":")
    if defect and (target not in _SOLVERS or not variant):
        raise SystemExit(
            f"--defect takes '<case>:<variant>' with a case from "
            f"{', '.join(sorted(_SOLVERS))}; got {defect!r}")
    out = {}
    for case, solver in _SOLVERS.items():
        if only and case != only:
            continue
        out[case] = solver(ws / fx.ROOT / case, fx,
                           variant if case == target else None)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("workspace")
    parser.add_argument("--defect", default="",
                        help="one case's plausible wrong answer instead")
    parser.add_argument("--only", default="")
    args = parser.parse_args()
    done = solve(args.workspace, args.defect or None, args.only or None)
    for case, msg in done.items():
        print(f"{case}: {msg}")


if __name__ == "__main__":
    main()
