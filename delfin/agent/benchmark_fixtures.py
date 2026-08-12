"""Workbook fixtures for the office benchmark, built from a readable spec.

The office suite scored a model against four small CSV files. Real
administrative work does not arrive as CSV: it arrives as ``.xlsx`` with
a merged block for a heading, a filter someone left switched on, a total
row at the bottom, a money column typed by hand, and more rows than any
reader shows at once. Every trap the office module learned to detect
this year lives in that format, and not one of them could fire against
the CSVs — so the suite measured a format the user does not work in.

Checking the workbooks into the repository was the obvious fix and the
wrong one: a binary fixture is a file nobody can review in a diff, and
the same objection is written into the workspace README. So the content
lives HERE, in a spec that reads as a table, and the runner materialises
it before the snapshot is taken. The diff stays reviewable; the format
under test becomes the one the user actually has on disk.

Three properties this has to keep:

* **Deterministic.** A benchmark baseline is only comparable against
  fixtures that are byte-stable, so nothing here is random and nothing
  is derived from the clock.
* **Rebuilt when the spec changes.** A stamp file carries a digest of
  the materialised rows. Stale workbooks from an older spec would score
  a model against a file the tasks no longer describe.
* **Loud when it cannot run.** Writing ``.xlsx`` needs the spreadsheet
  extra. Without it the workbook tasks cannot be measured at all, and
  the one thing that must not happen is that they quietly score as
  failures and drag a baseline down for a reason nobody can see.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any

OFFICE_WS_REL = Path("tests") / "fixtures" / "office_workspace"
STAMP_NAME = ".fixture-stamp"

# ---------------------------------------------------------------------------
# The spec
# ---------------------------------------------------------------------------
#
# Names, headings and labels are German because the files they stand in
# for are: the tasks are answered in the user's language and the reader's
# label matching is German-aware. The code around them is English.


def _kostenstellen_rows() -> list[list[Any]]:
    """A cost-centre overview whose grouping lives in merged cells.

    Rows 3-5 belong to 4711, rows 6-7 to 4712. Only the FIRST row of each
    block carries the number; the rest are empty because the cell is
    merged over them. A reader that does not report the merges shows
    those rows with a blank cost centre, and the obvious reading — that
    they are unassigned — is wrong.
    """
    return [
        ["Kostenstellenübersicht 2026", None, None, None],
        ["Kostenstelle", "Bezeichnung", "Position", "Budget"],
        [4711, "Anorganische Chemie", "Verbrauchsmaterial", 18500],
        [None, None, "Geräte", 42000],
        [None, None, "Wartung", 7300],
        [4712, "Beschaffung", "Verbrauchsmaterial", 5400],
        [None, None, "Fremdleistungen", 12750],
        [4713, "Verwaltung", "Bürobedarf", 3100],
    ]


_KOSTENSTELLEN_MERGES = ("A1:D1", "A3:A5", "B3:B5", "A6:A7", "B6:B7")


def _buchungen_rows() -> list[list[Any]]:
    """Bookings ending in a total row, with a filter left switched on.

    The total is a formula and the file has never been opened in a
    spreadsheet program, so it carries no cached value: the number the
    row promises does not exist in the file yet.
    """
    rows: list[list[Any]] = [
        ["Beleg", "Datum", "Kostenstelle", "Kreditor", "Betrag"],
    ]
    bookings = [
        ("R-001", "03.01.2026", 4711, "Meier GmbH", 1234.50),
        ("R-002", "07.01.2026", 4712, "Schulze AG", 289.90),
        ("R-003", "09.01.2026", 4711, "Laborbedarf Nord", 450.00),
        ("R-004", "14.01.2026", 4713, "Bürowelt", 96.40),
        ("R-005", "21.01.2026", 4711, "Meier GmbH", 780.00),
        ("R-006", "02.02.2026", 4712, "Chemie Direkt", 2145.75),
        ("R-007", "05.02.2026", 4711, "Laborbedarf Nord", 318.20),
        ("R-008", "11.02.2026", 4713, "Bürowelt", 142.00),
        ("R-009", "17.02.2026", 4712, "Schulze AG", 505.60),
        ("R-010", "24.02.2026", 4711, "Meier GmbH", 1670.00),
        ("R-011", "03.03.2026", 4711, "Glastechnik Süd", 894.30),
        ("R-012", "06.03.2026", 4712, "Chemie Direkt", 3320.00),
        ("R-013", "12.03.2026", 4713, "Reinigung Kern", 410.00),
        ("R-014", "19.03.2026", 4711, "Laborbedarf Nord", 265.85),
        ("R-015", "27.03.2026", 4712, "Schulze AG", 1180.40),
        ("R-016", "02.04.2026", 4711, "Meier GmbH", 940.00),
        ("R-017", "09.04.2026", 4713, "Bürowelt", 78.90),
        ("R-018", "15.04.2026", 4711, "Glastechnik Süd", 1512.60),
        ("R-019", "22.04.2026", 4712, "Chemie Direkt", 2760.00),
        ("R-020", "28.04.2026", 4711, "Laborbedarf Nord", 386.15),
        ("R-021", "05.05.2026", 4713, "Reinigung Kern", 410.00),
        ("R-022", "13.05.2026", 4711, "Meier GmbH", 1195.00),
        ("R-023", "20.05.2026", 4712, "Schulze AG", 642.30),
        ("R-024", "26.05.2026", 4711, "Glastechnik Süd", 2088.40),
    ]
    rows.extend([list(b) for b in bookings])
    last = 1 + len(bookings)
    rows.append(["Summe", None, None, None, f"=SUM(E2:E{last})"])
    return rows


# Someone filtered on cost centre 4711 and saved the file that way. These
# are the rows the filter hid — they are still in the file, still in the
# total, and invisible on screen.
_BUCHUNGEN_HIDDEN_ROWS = (3, 5, 10, 12, 16, 20, 24)


def _gutschriften_rows() -> list[list[Any]]:
    """A money column typed by a person, not by a program.

    Credit notes appear in two conventions in the same column — a
    trailing minus, and parentheses — and the rest are plain. The column
    is text throughout, so a reader that parses it naively either drops
    the credits or counts them as income, and the balance comes out
    positive when it is negative.
    """
    return [
        ["Beleg", "Kunde", "Datum", "Betrag"],
        ["A-2026-001", "Institut für Katalyse", "08.01.2026", "4.820,00"],
        ["A-2026-002", "Werkstoffprüfung Hamm", "15.01.2026", "1.240,50"],
        ["G-2026-001", "Institut für Katalyse", "22.01.2026", "1.500,00-"],
        ["A-2026-003", "Analytik Rhein", "04.02.2026", "962,40"],
        ["G-2026-002", "Werkstoffprüfung Hamm", "11.02.2026", "(340,00)"],
        ["A-2026-004", "Institut für Katalyse", "19.02.2026", "7.310,00"],
        ["G-2026-003", "Analytik Rhein", "26.02.2026", "962,40-"],
        ["A-2026-005", "Werkstoffprüfung Hamm", "05.03.2026", "2.155,80"],
        ["G-2026-004", "Institut für Katalyse", "12.03.2026", "(12.400,00)"],
        ["A-2026-006", "Analytik Rhein", "23.03.2026", "588,00"],
    ]


def _verbrauch_rows() -> list[list[Any]]:
    """More rows than a default read shows, with the maximum outside it.

    A reader shows the first 200 rows and says so. The largest value sits
    at row 243. Answering "which device consumes the most" from the first
    window gives a confident, wrong device — the note about the remaining
    rows is the only thing standing between the two.
    """
    rows: list[list[Any]] = [["Gerät", "Monat", "Verbrauch_kWh"]]
    devices = ("Trockenschrank", "Zentrifuge", "Abzug", "Kühlschrank",
               "Ofen", "Reinstwasser", "Pumpstand", "Schüttler")
    # A fixed, boring progression: reproducible, and nothing in the first
    # window comes near the outlier below.
    for i in range(260):
        device = devices[i % len(devices)]
        month = f"2026-{(i % 12) + 1:02d}"
        value = 120 + (i * 7) % 180
        rows.append([f"{device} {i // len(devices) + 1:02d}", month, value])
    # Row 243 in the sheet == index 242 in `rows` (row 1 is the header).
    rows[242][2] = 9840
    return rows


_WORKBOOKS: tuple[tuple[str, str, Any], ...] = (
    ("Kostenstellen.xlsx", "Übersicht", _kostenstellen_rows),
    ("Buchungen_2026.xlsx", "Buchungen", _buchungen_rows),
    ("Gutschriften.xlsx", "Offene Posten", _gutschriften_rows),
    ("Verbrauch_2026.xlsx", "Verbrauch", _verbrauch_rows),
)

WORKBOOK_NAMES: tuple[str, ...] = tuple(n for n, _s, _f in _WORKBOOKS)


# ---------------------------------------------------------------------------
# Digest and dependency
# ---------------------------------------------------------------------------


def spec_digest() -> str:
    """A digest of what the spec actually materialises.

    Hashing the rows rather than the source text means a comment change
    does not force a rebuild, and a changed number does.
    """
    payload = {
        "workbooks": [
            [name, sheet, rows()] for name, sheet, rows in _WORKBOOKS
        ],
        "merges": list(_KOSTENSTELLEN_MERGES),
        "hidden": list(_BUCHUNGEN_HIDDEN_ROWS),
    }
    blob = json.dumps(payload, sort_keys=True, default=str,
                      ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:16]


def missing_dependency_reason() -> str:
    """Empty when the workbooks can be written; otherwise why not."""
    try:
        import openpyxl  # noqa: F401
    except Exception as exc:
        return (f"openpyxl is not importable ({exc}) — install the "
                f"spreadsheet extra to measure the workbook tasks")
    return ""


# ---------------------------------------------------------------------------
# Materialising
# ---------------------------------------------------------------------------


def _write_workbook(target: Path, sheet_name: str,
                    rows: list[list[Any]], *,
                    merges: tuple[str, ...] = (),
                    hidden_rows: tuple[int, ...] = ()) -> None:
    """Write one workbook atomically.

    Atomically because two runs can be in flight at once, and a reader
    that opens a half-written workbook fails in a way that looks like a
    defect in the reader.
    """
    import openpyxl

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet_name
    for row in rows:
        ws.append(row)
    for ref in merges:
        ws.merge_cells(ref)
    for idx in hidden_rows:
        ws.row_dimensions[idx].hidden = True
    fd, tmp = tempfile.mkstemp(prefix=".wb-", suffix=".xlsx",
                               dir=str(target.parent))
    os.close(fd)
    try:
        wb.save(tmp)
        os.replace(tmp, target)
    finally:
        wb.close()
        if os.path.exists(tmp):
            os.unlink(tmp)


def fixtures_are_current(workspace: Path) -> bool:
    """True when every workbook exists and the stamp matches the spec."""
    stamp = workspace / STAMP_NAME
    try:
        if stamp.read_text(encoding="utf-8").strip() != spec_digest():
            return False
    except OSError:
        return False
    return all((workspace / name).is_file() for name in WORKBOOK_NAMES)


def ensure_office_fixtures(root: Path | str | None = None) -> tuple[list[str], str]:
    """Materialise the workbooks under ``root``'s office workspace.

    Returns ``(names_written, reason_not_written)``. Exactly one side is
    populated: either the workbooks are there and current, or the reason
    they are not is a sentence the caller can print.

    Must be called BEFORE the fixture snapshot is taken. The guard
    restores the workspace to whatever it held when the snapshot ran, so
    a workbook written afterwards is deleted again after the first
    attempt — and the task that needs it then fails for a reason that
    has nothing to do with the model.
    """
    base = Path(root) if root is not None else Path(os.getcwd())
    workspace = base / OFFICE_WS_REL
    if not workspace.is_dir():
        return [], f"no office workspace at {workspace}"
    if fixtures_are_current(workspace):
        return [], ""
    reason = missing_dependency_reason()
    if reason:
        return [], reason

    written: list[str] = []
    for name, sheet, rows in _WORKBOOKS:
        merges = (_KOSTENSTELLEN_MERGES if name == "Kostenstellen.xlsx"
                  else ())
        hidden = (_BUCHUNGEN_HIDDEN_ROWS if name == "Buchungen_2026.xlsx"
                  else ())
        try:
            _write_workbook(workspace / name, sheet, rows(),
                            merges=merges, hidden_rows=hidden)
        except Exception as exc:
            return written, f"could not write {name}: {exc}"
        written.append(name)
    try:
        (workspace / STAMP_NAME).write_text(
            spec_digest() + "\n", encoding="utf-8")
    except OSError as exc:
        return written, f"could not write the fixture stamp: {exc}"
    return written, ""


__all__ = [
    "OFFICE_WS_REL",
    "STAMP_NAME",
    "WORKBOOK_NAMES",
    "ensure_office_fixtures",
    "fixtures_are_current",
    "missing_dependency_reason",
    "spec_digest",
]
