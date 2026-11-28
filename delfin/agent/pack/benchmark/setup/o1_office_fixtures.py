#!/usr/bin/env python3
"""Setup for the twelve everyday office tasks (package O1).

Every fixture is BUILT here, from a spec that reads as a table, into
``<workspace>/o1/<case>/``. Nothing is committed as a binary, and no
name in this file belongs to a person: the collaborators are invented
combinations of a surname and one initial, the given names are a pool
used only by the "clean a column" case, and no combination refers to
anyone.

WHY THIS SUITE EXISTS. The office mode is the one a person does
administrative work in, and it was measured as the weakest: on the
installation where people actually do that work, a request to give every
collaborator a sheet of their own failed in three sessions running
(empty turns, nine repeated reads, and finally a correct script the shell
could not run), and 243 of 330 failed turns across those sessions were
empty turns. The pack's tasks were written for that mode but graded on
the ANSWER -- a sentence about the work -- while the user is graded on
the FILE. These tasks are graded on the file only: ``accept/`` opens what
the attempt wrote and compares it with the same expectation this module
derives.

WHAT IS SYNTHETIC, AND WHAT IS NOT. The data is invented. The SHAPE is
not: it is the shape measured on the real workbook, and ``self_check()``
returns the magnitudes so a test can assert them instead of trusting this
comment -- 737 data rows, 192 distinct values in the collaborator column,
49 empty cells, 47 values carrying a character a sheet name may not
carry, 18 longer than the 31-character limit, 4 that differ from another
only in case or spacing, and 7 that name several people at once.

DETERMINISM. ``random.Random(<fixed seed>)`` per generation, never the
global RNG, and no clock: two builds of one case are byte-identical in
content, so an acceptance run compares against the same numbers every
time.

argv[1] is the workspace. Exits non-zero when the fixtures cannot be
built -- a task whose question was never put is reported as unmeasured,
never as a model failure.
"""
from __future__ import annotations

import csv
import io
import random
import sys
from pathlib import Path

ROOT = "o1"

#: Every case this module builds, with the subdirectory it lives in.
CASES = (
    "blaettern",            # one sheet per value of a column (the measured case)
    "zusammenfuehren",      # merge two tables on a key
    "dubletten",            # find and mark duplicate rows
    "monatssummen",         # monthly totals from a date column with decimal commas
    "serienbrief",          # a mail merge into .docx letters
    "formular",             # fill a PDF form from a table row
    "pdf_tabelle",          # extract a table out of a PDF into .xlsx
    "bereinigen",           # trim, unify spelling, split "Nachname, Vorname"
    "pivot",                # count per category
    "csv_deutsch",          # a German-locale CSV into .xlsx with real types
    "regional_rollup",      # (en) two regional sheets into one summary
    "rechnungs_dubletten",  # (en) duplicate supplier invoices, case-insensitive
)

#: The two cases whose brief is written in English. A person working in
#: an English office types an English request; the FILE the task is graded
#: on is the same kind of file either way.
ENGLISH_CASES = frozenset({"regional_rollup", "rechnungs_dubletten"})

_SEED = 20260101

# ---------------------------------------------------------------------------
# Invented names
# ---------------------------------------------------------------------------

#: Surnames chosen for being common enough that no reader takes one for a
#: particular person, and paired below with a single initial so that no
#: full name occurs.
_SURNAMES = (
    "Albers", "Berger", "Cordes", "Dreyer", "Ehlers", "Faber", "Gerdes",
    "Hansen", "Ihde", "Jansen", "Kessler", "Lindner", "Mertens", "Nolte",
    "Ostermann", "Peters", "Quandt", "Reimers", "Schulte", "Tietjen",
    "Ulrich", "Voss", "Wendt", "Zimmer",
)
_INITIALS = "ABCDEFGHJKLMNPQRSTUVWXYZ"

#: Used by the "clean a column" case only, which has to split a real
#: given name out of "Nachname, Vorname".
_GIVEN = ("Jörg", "Katrin", "Ulrike", "Bastian", "Meike", "Sönke", "Anja",
          "Henning", "Ilse", "Rüdiger", "Gesine", "Torben", "Martina",
          "Kai", "Petra", "Lennart", "Birte", "Arne", "Heike", "Olaf")

_DEPTS_SHORT = ("NW", "SO", "Mitte", "Lab", "Technik")
#: Long enough to break the 31-character limit and free of BOTH the person
#: separators ("und", "&", ";") and the characters a sheet title may not
#: carry. Two fixture categories overlapping is how a count stops
#: describing the column that was measured: a department called
#: "... Infrastruktur und Beschaffung" made the multi-person rule fire on a
#: single person, and a "/" separator put the same value in the
#: forbidden-character group as well.
_DEPTS_LONG = (
    "Abteilung Technische Infrastruktur Bereich Beschaffung",
    "Referat Haushalt Innere Dienste Nord Bezirk West",
    "Sachgebiet Veranstaltungsorganisation Bereich Logistik",
)

#: The two surnames the "several people in one cell" values are built
#: from. Deliberately outside _SURNAMES: a person named in a shared cell
#: must not also be reachable through a single-name cell, or the expected
#: sheet map would depend on which one was seen first.
_SHARED = ("Ahrens", "Bohnhoff")

CATEGORIES = ("Verbrauchsmaterial", "Geräte", "Wartung", "Reisen",
               "Fortbildung", "Fremdleistungen")
MONTHS = ("Januar", "Februar", "März", "April", "Mai", "Juni")

# ---------------------------------------------------------------------------
# The collaborator column: the shape that was measured
# ---------------------------------------------------------------------------

#: Index ranges in the 192 raw values. The counts in the module docstring
#: are these ranges' lengths, and ``self_check()`` proves it against the
#: built data rather than against this comment.
_FORBIDDEN_UNTIL = 47        # 0..46 carry a character a sheet name may not
_LONG_UNTIL = _FORBIDDEN_UNTIL + 18   # 47..64 are longer than 31 characters
_SHARED_UNTIL = _LONG_UNTIL + 7       # 65..71 name several people
_VARIANTS_FROM = 188         # 188..191 repeat 0..3 in another case/spacing
_RAW_COUNT = 192

#: A character openpyxl refuses in a sheet title, and Excel refuses too.
#: Kept as one string so the rule and the fixture cannot drift apart.
FORBIDDEN_IN_TITLE = '/\\?*[]:'
MAX_TITLE = 31


def _plain_people() -> list[str]:
    return [f"{s}, {i}." for s in _SURNAMES for i in _INITIALS]


def _shared_people() -> list[tuple[str, str]]:
    """Seven cells, each naming two people, built from the shared pool.

    Fourteen distinct people, so a cell naming two of them never names the
    same person twice: two from each surname, seven cells.
    """
    people = [f"{s}, {i}." for s in _SHARED for i in "KMSVHUP"]
    assert len(people) == 14, people
    return [(people[2 * idx], people[2 * idx + 1]) for idx in range(7)]


def collaborator_values() -> list[str]:
    """The 192 distinct raw values of the column, in a fixed order.

    Ordering matters: the merge rule for values that differ only in case
    or spacing says the FIRST spelling seen wins, so the order is part of
    the contract, not an accident of generation.
    """
    pool = _plain_people()
    pairs = _shared_people()
    out: list[str] = []
    for idx in range(_RAW_COUNT):
        if idx < _FORBIDDEN_UNTIL:
            base = pool[idx]
            out.append(f"{base} /{_DEPTS_SHORT[idx % len(_DEPTS_SHORT)]}")
        elif idx < _LONG_UNTIL:
            dept = _DEPTS_LONG[(idx - _FORBIDDEN_UNTIL) % len(_DEPTS_LONG)]
            out.append(f"{pool[idx]} - {dept}")
        elif idx < _SHARED_UNTIL:
            left, right = pairs[idx - _LONG_UNTIL]
            out.append(f"{left} und {right}")
        elif idx < _VARIANTS_FROM:
            out.append(pool[idx])
        else:
            # Four plain values repeated in another case and with stray
            # spacing: the same person, typed twice. Cloned from the PLAIN
            # range on purpose -- cloning a value that already carries a
            # forbidden character would make one messy value stand for two
            # of the measured categories at once, and the counts in
            # self_check() would stop describing the column that was
            # measured.
            base = pool[idx - _VARIANTS_FROM + 100]
            out.append("  " + base.upper() + " ")
    return out


def _empty_row_indexes(n_rows: int) -> list[int]:
    """The 49 rows whose collaborator cell is empty.

    Every fifteenth row starting at seven: with 737 rows that is exactly
    49, and it is spread over the sheet rather than clustered, which is
    how an unfilled column really looks.
    """
    return [i for i in range(n_rows) if i % 15 == 7]


#: 737 data rows, the size of the workbook the request failed on.
N_ORDERS = 737


def order_rows() -> list[list]:
    """The shopping list: header row first, then the 737 data rows."""
    values = collaborator_values()
    empty = set(_empty_row_indexes(N_ORDERS))
    rng = random.Random(_SEED)
    creditors = ("Nordlicht GmbH", "Weser Bezug", "Hansa Labor",
                 "Elbe Technik", "Kanal Handel", "Marsch Papier")
    rows: list[list] = [["Beleg", "Datum", "Kostenstelle", "Kreditor",
                         "Bearbeiter", "Betrag"]]
    for idx in range(N_ORDERS):
        value = "" if idx in empty else values[idx % _RAW_COUNT]
        day = 1 + (idx * 3) % 28
        month = 1 + (idx // 28) % 6
        rows.append([
            f"B-{idx + 1:04d}",
            f"{day:02d}.{month:02d}.2026",
            4711 + idx % 3,
            creditors[rng.randrange(len(creditors))],
            value,
            round(20 + (idx * 37) % 4000 + ((idx * 7) % 100) / 100, 2),
        ])
    return rows


def _slug(value: str) -> str:
    """The sheet title the house rules ask for, or ``""`` for an empty value.

    The rules (``HINWEISE.txt`` in the fixture) in the order they apply:
    trim and collapse spacing, replace every character a title may not
    carry with ``-``, cut to 31 characters.
    """
    text = " ".join(str(value or "").split())
    if not text:
        return ""
    for ch in FORBIDDEN_IN_TITLE:
        text = text.replace(ch, "-")
    return text[:MAX_TITLE]


def _people_in(value: str) -> list[str]:
    """The people one cell names. Several are told apart by a separator.

    A cell naming two people counts for each of them -- which is the part
    a per-value split gets wrong: it gives the pair a sheet of their own
    and leaves both people without one.
    """
    text = " ".join(str(value or "").split())
    if not text:
        return []
    for sep in (" und ", " & ", "; "):
        if sep in text:
            return [part.strip() for part in text.split(sep) if part.strip()]
    return [text]


def _merge_key(value: str) -> str:
    """Case and spacing do not tell two people apart; nothing else does."""
    return "".join(str(value).split()).lower()


def expected_sheets() -> dict[str, list[str]]:
    """Sheet title -> the ``Beleg`` numbers that belong in it.

    The single source of truth for the ``blaettern`` case: the acceptance
    script compares a produced workbook against this map, and this map is
    computed from the same rows the setup writes. A graded number is never
    typed twice.
    """
    title_for: dict[str, str] = {}
    taken: dict[str, str] = {}
    members: dict[str, list[str]] = {}
    for row in order_rows()[1:]:
        beleg, value = str(row[0]), str(row[4])
        if not str(value).strip():
            members.setdefault("Ohne Bearbeiter", []).append(beleg)
            continue
        for person in _people_in(value):
            key = _merge_key(person)
            if key not in title_for:
                base = _slug(person)
                title = base
                # A truncation can make two different people collide; the
                # rules say the later one gets a suffix rather than losing
                # its rows into somebody else's sheet.
                n = 2
                while title in taken and taken[title] != key:
                    suffix = f"~{n}"
                    title = base[:MAX_TITLE - len(suffix)] + suffix
                    n += 1
                taken[title] = key
                title_for[key] = title
            members.setdefault(title_for[key], []).append(beleg)
    return members


def self_check() -> dict:
    """The magnitudes this fixture claims, measured on what it built.

    Returned rather than asserted so the test can state them: a fixture
    whose messiness silently disappeared would still produce a passing
    acceptance script, and would then measure a question nobody asked.
    """
    values = collaborator_values()
    rows = order_rows()[1:]
    titles = expected_sheets()
    return {
        "raw_values": len(set(values)),
        "data_rows": len(rows),
        "empty_cells": len(_empty_row_indexes(N_ORDERS)),
        "with_forbidden_char": sum(
            1 for v in values if any(c in v for c in FORBIDDEN_IN_TITLE)),
        "longer_than_31": sum(1 for v in values if len(v.strip()) > MAX_TITLE),
        "case_or_space_variants": sum(
            1 for i, v in enumerate(values)
            if _merge_key(v) in {_merge_key(x) for x in values[:i]}),
        "naming_several_people": sum(
            1 for v in values if len(_people_in(v)) > 1),
        "sheets": len(titles),
        "rows_without_collaborator": len(titles.get("Ohne Bearbeiter", [])),
    }


# ---------------------------------------------------------------------------
# Case fixtures
# ---------------------------------------------------------------------------


def _write_csv(path: Path, rows: list[list], *, delimiter: str = ";",
               encoding: str = "utf-8") -> None:
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter=delimiter, lineterminator="\r\n")
    writer.writerows(rows)
    path.write_bytes(buf.getvalue().encode(encoding))


def _write_xlsx(path: Path, rows: list[list], *, sheet: str = "Tabelle") -> None:
    import openpyxl

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet
    for row in rows:
        ws.append(row)
    wb.save(str(path))
    wb.close()


def _build_blaettern(target: Path) -> None:
    _write_xlsx(target / "bestellungen.xlsx", order_rows(), sheet="Bestellungen")
    (target / "HINWEISE.txt").write_text(
        "Blattnamen\n"
        "----------\n"
        "Jeder Bearbeiter bekommt ein eigenes Blatt. Der Blattname ist der\n"
        "Zellinhalt in Spalte E, dabei gilt:\n"
        "  * Zeichen, die ein Blattname nicht enthalten darf (/ \\ ? * [ ] :),\n"
        "    werden durch '-' ersetzt.\n"
        "  * Hoechstens 31 Zeichen; wird ein Name dadurch zweimal vergeben,\n"
        "    bekommt der spaetere '~2', '~3' ... angehaengt.\n"
        "  * Gross-/Kleinschreibung und ueberzaehlige Leerzeichen\n"
        "    unterscheiden keine Personen; die erste Schreibweise gewinnt.\n"
        "  * Eine Zelle mit mehreren Namen ('und', '&', ';') zaehlt fuer\n"
        "    jede der genannten Personen.\n"
        "  * Leere Zellen kommen auf ein Blatt 'Ohne Bearbeiter'.\n"
        "Jedes Blatt hat eine Kopfzeile wie im Ausgangsblatt.\n",
        encoding="utf-8")


def join_rows() -> tuple[list[list], list[list], list[str]]:
    """Orders, the creditor master, and the creditors nobody ordered from."""
    rng = random.Random(_SEED + 1)
    creditors = [[f"K-{n:03d}", f"Lieferant {n:02d}",
                  ("Nordstadt", "Suedhafen", "Westend")[n % 3]]
                 for n in range(1, 41)]
    known = [c[0] for c in creditors[:34]]
    rows: list[list] = [["Beleg", "Kreditor_Nr", "Betrag"]]
    for idx in range(120):
        key = known[idx % len(known)] if idx % 20 else f"K-{900 + idx:03d}"
        rows.append([f"P-{idx + 1:03d}", key,
                     round(50 + idx * 13.25, 2)])
    return rows, creditors, []


def _build_zusammenfuehren(target: Path) -> None:
    rows, creditors, _ = join_rows()
    _write_xlsx(target / "bestellungen.xlsx", rows, sheet="Bestellungen")
    _write_xlsx(target / "kreditoren.xlsx",
                [["Kreditor_Nr", "Name", "Ort"]] + creditors,
                sheet="Kreditoren")
    (target / "HINWEISE.txt").write_text(
        "Zuordnung\n"
        "---------\n"
        "Die Bestellungen tragen eine Kreditorennummer, das Stammblatt\n"
        "Kreditoren die Anschrift. Nicht zu jeder Bestellung gibt es einen\n"
        "Stammdatensatz. Ergebnis: 'ergebnis_join.xlsx', Spalten Beleg,\n"
        "Name, Ort, Betrag, und die Bestellungen ohne Stammdatensatz\n"
        "zusaetzlich als Liste 'nicht_zugeordnet.txt' (eine Belegnummer pro\n"
        "Zeile).\n",
        encoding="utf-8")


def duplicate_rows() -> tuple[list[list], list[int]]:
    """A register with 18 true duplicate pairs and 9 pairs that only look it.

    The row indexes returned are 1-based Excel rows of every row that
    belongs to a group with two identical members, first row included --
    the marking rule the fixture states.
    """
    rng = random.Random(_SEED + 2)
    base: list[list] = [["Beleg", "Datum", "Lieferant", "Betrag"]]
    for idx in range(264):
        base.append([f"E-{idx + 1:04d}",
                     f"{1 + idx % 28:02d}.{1 + idx // 28:02d}.2026",
                     f"Lieferant {rng.randrange(1, 20):02d}",
                     round(30 + idx * 7.5, 2)])
    dupes: list[list] = []
    # 18 exact copies, appended: same Beleg, same date, same amount.
    for idx in (3, 11, 19, 27, 33, 41, 52, 60, 68, 77, 85, 93, 101, 110,
                118, 126, 135, 143):
        dupes.append(list(base[idx]))
    # 9 near-misses: same Beleg and date, other amount. Not duplicates.
    for idx in (7, 25, 49, 71, 88, 104, 122, 139, 155):
        row = list(base[idx])
        row[3] = round(float(row[3]) + 0.01, 2)
        dupes.append(row)
    body = base[1:] + dupes
    rng.shuffle(body)
    rows = base[:1] + body
    marked: list[int] = []
    seen: dict[tuple, list[int]] = {}
    for pos, row in enumerate(rows[1:], start=2):
        key = (str(row[0]), str(row[1]), str(row[3]))
        seen.setdefault(key, []).append(pos)
    for key, positions in seen.items():
        if len(positions) > 1:
            marked.extend(positions)
    return rows, sorted(marked)


def _build_dubletten(target: Path) -> None:
    rows, _marked = duplicate_rows()
    _write_xlsx(target / "eingang.xlsx", rows, sheet="Eingang")
    (target / "HINWEISE.txt").write_text(
        "Dubletten\n"
        "---------\n"
        "Eine Zeile ist eine Dublette, wenn Beleg, Datum UND Betrag genau\n"
        "so schon vorkommen. Weicht der Betrag ab, ist es keine Dublette.\n"
        "Ergebnis: 'ergebnis_dubletten.xlsx' -- dasselbe Blatt wie\n"
        "'eingang.xlsx', zusaetzlich eine Spalte 'Dublette', die in JEDER\n"
        "Zeile einer solchen Gruppe 'ja' enthaelt (auch in der ersten) und\n"
        "sonst leer bleibt.\n",
        encoding="utf-8")


def donation_rows() -> list[list]:
    """240 rows whose dates and amounts are typed the way people type them."""
    rows: list[list] = [["Datum", "Spender", "Betrag"]]
    for idx in range(240):
        day = 1 + (idx * 5) % 28
        month = 1 + (idx // 40) % 6
        # Every fourth row was typed in ISO form by an import, the rest by
        # hand with dots. A reader that knows one form silently drops 60
        # rows -- and the total still looks like a total.
        datum = (f"2026-{month:02d}-{day:02d}" if idx % 4 == 0
                 else f"{day:02d}.{month:02d}.2026")
        cents = (idx * 3709) % 500000
        sign = "-" if idx % 23 == 0 else ""
        # German thousands separator AND decimal comma, because that is what
        # the column looks like: a naive float() reads "1.234,56" as 1.234
        # and the total is then wrong in the fourth digit only -- the kind of
        # error that survives a review. A leading minus is a repayment.
        euros, rest = cents // 100, cents % 100
        thousand = f"{euros:,}".replace(",", ".")
        rows.append([datum, f"Spender {idx + 1:03d}",
                     f"{sign}{thousand},{rest:02d}"])
    return rows


def expected_months() -> dict[str, tuple[float, int]]:
    """Month name -> (sum in EUR, number of usable rows).

    Derived from the same rows the setup writes, and deliberately by a
    route that does NOT reuse the acceptance script's parsing: the sum is
    accumulated here from the literal strings, so a mistake in one side
    shows up as a disagreement rather than cancelling out.
    """
    totals: dict[str, list] = {}
    for datum, _spender, amount in donation_rows()[1:]:
        if "-" in datum:
            month = int(datum.split("-")[1])
        else:
            month = int(datum.split(".")[1])
        cents = round(float(amount.replace(".", "").replace(",", ".")) * 100)
        entry = totals.setdefault(MONTHS[month - 1], [0, 0])
        entry[0] += cents
        entry[1] += 1
    return {name: (cents / 100.0, n) for name, (cents, n) in totals.items()}


def _build_monatssummen(target: Path) -> None:
    _write_csv(target / "spenden.csv", donation_rows())
    (target / "HINWEISE.txt").write_text(
        "Monatssummen\n"
        "------------\n"
        "Die Spalte Datum ist gemischt: von Hand mit Punkten geschrieben,\n"
        "per Import im ISO-Format. Beides ist derselbe Tag und muss gleich\n"
        "behandelt werden. Betraege sind Dezimalzahlen mit Komma; ein\n"
        "vorangestelltes Minus ist eine Rueckzahlung und wird abgezogen.\n"
        "Ergebnis: 'monatssummen.xlsx', Spalten Monat und Summe, eine Zeile\n"
        "je Monat von Januar bis Juni 2026, Summe als ZAHL (nicht als Text),\n"
        "auf den Cent genau.\n",
        encoding="utf-8")


def mailmerge_rows() -> list[list]:
    rows: list[list] = [["Anrede", "Vorname", "Nachname", "Betrag", "Ort"]]
    for idx in range(12):
        rows.append([
            "Frau" if idx % 2 == 0 else "Herr",
            _GIVEN[idx % len(_GIVEN)],
            _SURNAMES[idx],
            f"{120 + idx * 17},{(idx * 31) % 100:02d}",
            ("Nordstadt", "Suedhafen", "Westend", "Ostheim")[idx % 4],
        ])
    return rows


def _build_serienbrief(target: Path) -> None:
    from docx import Document

    rows = mailmerge_rows()
    _write_xlsx(target / "teilnehmer.xlsx", rows, sheet="Teilnehmer")
    doc = Document()
    doc.add_paragraph("Sehr geehrte {{Anrede}} {{Vorname}} {{Nachname}},")
    doc.add_paragraph("")
    doc.add_paragraph("wir bestaetigen Ihre Anmeldung fuer {{Ort}}.")
    doc.add_paragraph("Der Teilnahmebeitrag betraegt {{Betrag}} EUR.")
    doc.add_paragraph("")
    doc.add_paragraph("Mit freundlichen Gruessen")
    doc.add_paragraph("Team Anmeldung")
    doc.save(str(target / "vorlage.docx"))
    (target / "HINWEISE.txt").write_text(
        "Serienbrief\n"
        "----------\n"
        "Vorlage 'vorlage.docx' mit Platzhaltern der Form {{Feld}}. Fuer\n"
        "jede Zeile aus 'teilnehmer.xlsx' ein Brief als\n"
        "'briefe/<Nachname>_<Vorname>.docx', in dem alle Platzhalter durch\n"
        "die Werte der Zeile ersetzt sind -- im Brief darf kein '{{' mehr\n"
        "stehen.\n",
        encoding="utf-8")


def application_rows() -> list[list]:
    rows: list[list] = [["Antragsnummer", "Name", "Kostenstelle", "Betrag",
                         "Datum"]]
    for idx in range(3):
        rows.append([f"A-{idx + 1:03d}", _SURNAMES[idx],
                     4711 + idx,
                     f"{250 + idx * 100},{(idx * 17) % 100:02d}",
                     f"{1 + idx:02d}.03.2026"])
    return rows


def _build_formular(target: Path) -> None:
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    _write_csv(target / "antraege.csv", application_rows())
    path = target / "antrag.pdf"
    c = canvas.Canvas(str(path), pagesize=A4)
    c.setFont("Helvetica", 12)
    c.drawString(60, 800, "Antrag auf Kostenerstattung")
    for idx, name in enumerate(("Name", "Kostenstelle", "Betrag", "Datum")):
        c.setFont("Helvetica", 10)
        c.drawString(60, 742 - 34 * idx, name)
        c.acroForm.textfield(name=name, x=180, y=740 - 34 * idx,
                             width=260, height=18, fontSize=10)
    c.save()
    (target / "HINWEISE.txt").write_text(
        "Formular\n"
        "--------\n"
        "'antrag.pdf' ist ein Formular mit den Feldern Name, Kostenstelle,\n"
        "Betrag und Datum. 'antraege.csv' enthaelt je eine Zeile pro Antrag.\n"
        "Ergebnis: je Antrag eine ausgefuellte Kopie unter\n"
        "'antraege_ausgefuellt/<Antragsnummer>.pdf', in der die vier Felder\n"
        "den Werten der Zeile entsprechen (Betrag und Datum so, wie sie in\n"
        "der CSV stehen).\n",
        encoding="utf-8")


def hours_rows() -> list[list]:
    """24 data rows over two pages, each page repeating the header."""
    rows: list[list] = []
    for idx in range(24):
        day = 1 + idx
        rows.append([f"{day:02d}.03.2026",
                     f"Gruppe {'Alpha' if idx % 2 == 0 else 'Beta'}",
                     f"{7 + (idx % 3)},{'5' if idx % 4 else '0'}",
                     f"Projekt {1 + idx % 4}"])
    return rows


def _build_pdf_tabelle(target: Path) -> None:
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    rows = hours_rows()
    path = target / "stunden.pdf"
    c = canvas.Canvas(str(path), pagesize=A4)
    per_page = 12
    for page in range(2):
        c.setFont("Helvetica", 11)
        c.drawString(60, 800, "Stundenliste Maerz 2026")
        y = 770
        c.setFont("Helvetica-Bold", 10)
        for col, (x, head) in enumerate(((60, "Datum"), (180, "Gruppe"),
                                         (330, "Stunden"), (430, "Projekt"))):
            c.drawString(x, y, head)
        c.setFont("Helvetica", 10)
        y -= 18
        for row in rows[page * per_page:(page + 1) * per_page]:
            for x, cell in zip((60, 180, 330, 430), row):
                c.drawString(x, y, cell)
            y -= 18
        c.setFont("Helvetica-Oblique", 8)
        c.drawString(60, 60, f"Seite {page + 1} von 2")
        c.showPage()
    c.save()
    (target / "HINWEISE.txt").write_text(
        "Tabelle aus PDF\n"
        "---------------\n"
        "'stunden.pdf' hat zwei Seiten; auf jeder Seite steht die\n"
        "Kopfzeile erneut und unten eine Fusszeile 'Seite x von 2'.\n"
        "Ergebnis: 'stunden.xlsx' mit den 24 Datenzeilen und einer\n"
        "Kopfzeile -- die Kopfzeilenwiederholung und die Fusszeile gehoeren\n"
        "nicht in die Tabelle. Spalte Stunden als Zahl (7,5 -> 7.5).\n",
        encoding="utf-8")


def messy_people_rows() -> list[list]:
    """90 rows of what a column looks like after ten people typed into it."""
    rows: list[list] = [["Personal"]]
    for idx in range(90):
        name = _SURNAMES[idx % len(_SURNAMES)]
        given = _GIVEN[idx % len(_GIVEN)]
        form = idx % 6
        if form == 0:
            cell = f" {name} , {given} "
        elif form == 1:
            cell = f"{name.lower()}, {given.lower()}"
        elif form == 2:
            # ue/oe/ae typed by somebody whose keyboard had no umlaut.
            cell = (f"{name.replace('ö', 'oe').replace('ü', 'ue')}, "
                    f"{given.replace('ö', 'oe').replace('ü', 'ue')}")
        elif form == 3:
            cell = f"Dr. {name}, {given}"
        elif form == 4:
            cell = f"{name},   {given}  "
        else:
            cell = "" if idx % 30 == 5 else f"{name}, {given}"
        rows.append([cell])
    return rows


def expected_clean() -> list[tuple[str, str]]:
    """(Nachname, Vorname) per input row, in input order.

    The rule the fixture states: trim, drop a leading title, umlauts
    written as ue/oe/ae belong to the umlaut spelling, one output row per
    input row, order unchanged. An empty input row yields ("", "").
    """
    out: list[tuple[str, str]] = []
    for (cell,) in messy_people_rows()[1:]:
        text = " ".join(str(cell).split())
        for title in ("Dr. ", "Prof. "):
            if text.startswith(title):
                text = text[len(title):]
        if not text:
            out.append(("", ""))
            continue
        surname, _, given = text.partition(",")
        out.append((_to_umlaut(surname.strip()), _to_umlaut(given.strip())))
    return out


def _to_umlaut(text: str) -> str:
    """``ue`` -> ``ü`` and friends, then capitalise the first letter.

    Only the digraphs are rewritten -- a genuine ``ue`` inside a name is
    not distinguishable, and the fixture's names are chosen so it does not
    occur. The capitalisation restores what lower-casing the cell lost.
    """
    out = (text.replace("oe", "ö").replace("ue", "ü").replace("ae", "ä")
               .replace("Oe", "Ö").replace("Ue", "Ü").replace("Ae", "Ä"))
    return out[:1].upper() + out[1:]


def _build_bereinigen(target: Path) -> None:
    _write_xlsx(target / "personen.xlsx", messy_people_rows(), sheet="Personal")
    (target / "HINWEISE.txt").write_text(
        "Personenliste bereinigen\n"
        "------------------------\n"
        "Die Spalte Personal ist ueber Jahre gefuellt worden. Regeln:\n"
        "  * Leerzeichen am Rand und mehrfache Leerzeichen entfernen.\n"
        "  * Ein fuehrender Titel ('Dr. ', 'Prof. ') wird nicht uebernommen.\n"
        "  * 'ue', 'oe', 'ae' sind Umlaute, wenn der Name einen hat:\n"
        "    Mueller -> Müller.\n"
        "  * Aufteilung an der Form 'Nachname, Vorname'.\n"
        "  * Eine Ausgabezeile je Eingabezeile, Reihenfolge unveraendert;\n"
        "    eine leere Eingabezeile bleibt eine leere Ausgabezeile.\n"
        "Ergebnis: 'personen_bereinigt.xlsx', Spalten Nachname und Vorname.\n",
        encoding="utf-8")


def sales_rows() -> list[list]:
    rng = random.Random(_SEED + 3)
    rows: list[list] = [["Kategorie", "Monat", "Betrag"]]
    for idx in range(200):
        rows.append([CATEGORIES[idx % len(CATEGORIES)],
                     MONTHS[(idx // 6) % len(MONTHS)],
                     round(20 + rng.randrange(5000) / 10, 2)])
    return rows


def expected_pivot() -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for category, month, _amount in sales_rows()[1:]:
        counts.setdefault(category, {}).setdefault(month, 0)
        counts[category][month] += 1
    return counts


def _build_pivot(target: Path) -> None:
    _write_xlsx(target / "verkaeufe.xlsx", sales_rows(), sheet="Verkaeufe")
    (target / "HINWEISE.txt").write_text(
        "Haeufigkeiten\n"
        "-------------\n"
        "Wie viele Verkaeufe gibt es je Kategorie und Monat?\n"
        "Ergebnis: 'pivot.xlsx' -- eine Zeile je Kategorie, eine Spalte je\n"
        "Monat Januar bis Juni, in den Zellen die ANZAHL der Verkaeufe als\n"
        "Zahl. Die erste Spalte traegt die Kategorie, die erste Zeile die\n"
        "Monatsnamen.\n",
        encoding="utf-8")


def revenue_rows() -> list[list]:
    rows: list[list] = [["Kostenstelle", "Bezeichnung", "Datum", "Betrag"]]
    for idx in range(40):
        day = 1 + (idx * 2) % 28
        # Up to 9,000 EUR, so the thousands dot occurs as HINWEISE says.
        cents = (idx * 23417) % 900000
        euros, rest = cents // 100, cents % 100
        rows.append([
            # A typed number keeps its leading zeros only as text; as a
            # number 004711 becomes 4711, and the Kostenstelle stops
            # matching every other table.
            f"{4711 + idx % 2:06d}",
            ("Geräte & Zubehör" if idx % 3 == 0 else "Verbrauchsmaterial"),
            f"{day:02d}.{(idx % 6) + 1:02d}.2026",
            f"{euros:,}".replace(",", ".") + f",{rest:02d}",
        ])
    return rows


def expected_csv_deutsch() -> list[tuple[str, str, str, float]]:
    """(Kostenstelle, Bezeichnung, Datum ISO, Betrag) per data row."""
    out = []
    for number, label, datum, amount in revenue_rows()[1:]:
        day, month, year = datum.split(".")
        out.append((number, label, f"{year}-{month}-{day}",
                    float(amount.replace(".", "").replace(",", "."))))
    return out


def _build_csv_deutsch(target: Path) -> None:
    _write_csv(target / "umsatz.csv", revenue_rows(), encoding="cp1252")
    (target / "HINWEISE.txt").write_text(
        "CSV uebernehmen\n"
        "---------------\n"
        "'umsatz.csv' kommt aus einer deutschen Buchhaltung: Semikolon als\n"
        "Trennzeichen, Windows-Zeichensatz (cp1252), Betrag mit\n"
        "Tausendertrennzeichen und Dezimalkomma, Datum als TT.MM.JJJJ.\n"
        "Ergebnis: 'umsatz.xlsx' mit denselben vier Spalten, dabei\n"
        "  * Betrag als ZAHL (1234.56), nicht als Text,\n"
        "  * Datum als Datum (nicht als Text),\n"
        "  * Kostenstelle als TEXT, damit die fuehrenden Nullen bleiben.\n",
        encoding="utf-8")


def regional_rows() -> tuple[list[list], list[list]]:
    rng = random.Random(_SEED + 4)
    nord = [["Kategorie", "Umsatz"]]
    sued = [["Kategorie", "Umsatz"]]
    for idx, category in enumerate(CATEGORIES):
        nord.append([category, round(1000 + idx * 211.5, 2)])
        # Two categories only in the south: a left-join on that side
        # loses them, and the total still adds up.
        if idx < 4:
            sued.append([category, round(500 + rng.randrange(900), 2)])
    sued.append(["Auftragsfertigung", 812.4])
    return nord, sued


def expected_rollup() -> dict[str, float]:
    nord, sued = regional_rows()
    totals: dict[str, float] = {}
    for rows in (nord, sued):
        for category, amount in rows[1:]:
            totals[category] = round(totals.get(category, 0.0) + float(amount), 2)
    return totals


def _build_regional_rollup(target: Path) -> None:
    nord, sued = regional_rows()
    _write_xlsx(target / "region_nord.xlsx", nord, sheet="Nord")
    _write_xlsx(target / "region_sued.xlsx", sued, sheet="Sued")
    (target / "NOTES.txt").write_text(
        "Regional roll-up\n"
        "----------------\n"
        "The two regional files were exported from different systems, so\n"
        "they do not hold the same categories. Result: 'zusammenfassung.xlsx'\n"
        "with columns Kategorie and Umsatz, one row per category that\n"
        "occurs in EITHER file, Umsatz as a NUMBER and equal to the sum of\n"
        "both files.\n",
        encoding="utf-8")


def register_rows() -> tuple[list[list], list[int]]:
    """Supplier invoices, with 14 duplicates that differ in case/spacing."""
    rng = random.Random(_SEED + 5)
    rows: list[list] = [["Rechnungsnummer", "Lieferant", "Betrag"]]
    for idx in range(200):
        rows.append([f"RE-{1000 + idx}", f"Supplier {rng.randrange(1, 25):02d}",
                     round(60 + idx * 11.4, 2)])
    dupe_idx = (5, 14, 23, 37, 48, 59, 66, 78, 91, 103, 117, 129, 141, 158)
    for idx in dupe_idx:
        row = list(rows[idx + 1])
        # The same invoice, typed by a second person: other case, a
        # trailing space. A comparison that is not case- and
        # space-insensitive calls these two invoices.
        row[0] = (" " + str(row[0]).lower()) if idx % 2 else (str(row[0]) + " ")
        rows.append(row)
    marked: list[int] = []
    seen: dict[str, list[int]] = {}
    for pos, row in enumerate(rows[1:], start=2):
        key = "".join(str(row[0]).split()).upper()
        seen.setdefault(key, []).append(pos)
    for positions in seen.values():
        if len(positions) > 1:
            marked.extend(positions)
    return rows, sorted(marked)


def _build_rechnungs_dubletten(target: Path) -> None:
    rows, _marked = register_rows()
    _write_xlsx(target / "register.xlsx", rows, sheet="Register")
    (target / "NOTES.txt").write_text(
        "Duplicate invoices\n"
        "------------------\n"
        "Two rows are the same invoice when the Rechnungsnummer is the same,\n"
        "ignoring case and spaces. Result: 'register_dubletten.xlsx' -- the\n"
        "same sheet plus a column 'Dublette' that says 'ja' in EVERY row of\n"
        "such a group (the first one included) and stays empty otherwise.\n",
        encoding="utf-8")


_BUILDERS = {
    "blaettern": _build_blaettern,
    "zusammenfuehren": _build_zusammenfuehren,
    "dubletten": _build_dubletten,
    "monatssummen": _build_monatssummen,
    "serienbrief": _build_serienbrief,
    "formular": _build_formular,
    "pdf_tabelle": _build_pdf_tabelle,
    "bereinigen": _build_bereinigen,
    "pivot": _build_pivot,
    "csv_deutsch": _build_csv_deutsch,
    "regional_rollup": _build_regional_rollup,
    "rechnungs_dubletten": _build_rechnungs_dubletten,
}


def case_dir(workspace: Path | str, case: str) -> Path:
    return Path(workspace) / ROOT / case


def build(workspace: Path | str, cases=None) -> list[str]:
    """Build the fixtures for ``cases`` (all by default). Returns the names."""
    ws = Path(workspace)
    root = ws / ROOT
    if root.exists():
        raise SystemExit(f"refusing to overwrite {root}")
    written: list[str] = []
    for case in (cases or CASES):
        builder = _BUILDERS.get(case)
        if builder is None:
            raise SystemExit(f"unknown case {case!r}")
        target = root / case
        target.mkdir(parents=True)
        builder(target)
        written.append(case)
    return written


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: o1_office_fixtures.py <workspace> [case ...]",
              file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    cases = argv[2:] or None
    written = build(ws, cases)
    stats = self_check()
    print(f"prepared {ws / ROOT} -- {len(written)} cases "
          f"({', '.join(written)}); collaborator column: "
          f"{stats['raw_values']} distinct values over {stats['data_rows']} "
          f"rows, {stats['empty_cells']} empty, "
          f"{stats['with_forbidden_char']} with a character a sheet name may "
          f"not carry, {stats['longer_than_31']} over 31 characters, "
          f"{stats['case_or_space_variants']} case/space variants, "
          f"{stats['naming_several_people']} naming several people "
          f"-> {stats['sheets']} sheets expected")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
