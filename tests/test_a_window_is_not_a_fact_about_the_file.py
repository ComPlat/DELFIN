"""Three cuts, of which only one said what it was doing.

The reader shows a window of a table. Rows were signposted from the
start — "showing rows 1-200 of 40000 — pass start_row to page further" —
and that is the shape the other two never had.

**Cells** were cut at forty characters in silence. The grid is the only
representation of a cell the model ever sees, so

    Rechnung 2026-0001 Storno wegen Rücksend

was indistinguishable from a value that really ends there. The model
could quote it back as the cell's content, or copy it into an edit.

**Columns** said "showing 40 of 87 columns" and named nothing the caller
could do about it. The slice always began at column 1 and no reader took
a starting column at all, so columns 41 to 87 were unreachable through
this tool — and absent from the column profile, which therefore
described a subset of the file's columns without saying so.

**The profile** was computed from the window and rendered under a header
carrying the file's full row count, so

    40000 rows x 12 columns
    Betrag: number (decimal_comma), 3 of 200 not parseable

read as a property of the column. It is a property of rows 1-200. Worse,
the convention it decides there — decimal comma, day-first dates, what
counts as a number — is then applied to the whole file by the summing
and comparing tools.

A limit stated without a remedy is worse than no limit: it tells the
model something is missing and leaves it to guess, which in practice
means answering from the part it has.
"""

from __future__ import annotations

import pytest

from delfin.agent import office as O


def _notes(result) -> str:
    return " ".join(result.get("notes") or [])


# ---------------------------------------------------------------------------
# A cut cell says it was cut
# ---------------------------------------------------------------------------

def test_a_long_cell_ends_with_a_mark():
    rows = [["Beleg", "Zweck"],
            ["R-001", "Rechnung 2026-0001 Storno wegen Rücksendung Charge 7"]]
    grid = O.render_grid(rows)
    assert "…" in grid


def test_a_cell_that_fits_is_untouched():
    grid = O.render_grid([["Beleg"], ["R-001"]])
    assert "…" not in grid


def test_the_mark_does_not_widen_the_cell():
    long = "x" * 200
    grid = O.render_grid([[long]], max_cell_chars=40)
    cell = [ln for ln in grid.splitlines() if "x" in ln][0]
    assert cell.count("x") == 39 and "…" in cell


def test_the_reader_says_how_many_cells_were_cut(tmp_path):
    p = tmp_path / "belege.csv"
    p.write_text(
        "Beleg,Zweck\n"
        "R-001,Rechnung 2026-0001 Storno wegen Rücksendung Charge 7 Kunde X\n"
        "R-002,kurz\n", encoding="utf-8")
    text = _notes(O.read_sheet(p)).lower()
    assert "cell(s) are longer" in text
    assert "1 cell(s)" in text


def test_a_table_of_short_cells_says_nothing_about_clipping(tmp_path):
    p = tmp_path / "kurz.csv"
    p.write_text("a,b\n1,2\n", encoding="utf-8")
    assert "longer than" not in _notes(O.read_sheet(p))


def test_counting_clipped_cells_is_available_on_its_own():
    rows = [["x" * 60, "ok"], ["also fine", "y" * 41]]
    assert O.clipped_cell_count(rows, 40) == 2
    assert O.clipped_cell_count(rows, 100) == 0


# ---------------------------------------------------------------------------
# Columns can be paged across
# ---------------------------------------------------------------------------

@pytest.fixture
def wide(tmp_path):
    header = ",".join(f"S{i}" for i in range(1, 61))
    row = ",".join(str(i) for i in range(1, 61))
    p = tmp_path / "breit.csv"
    p.write_text(f"{header}\n{row}\n", encoding="utf-8")
    return p


def test_the_column_note_names_the_remedy(wide):
    text = _notes(O.read_sheet(wide, max_cols=40))
    assert "start_col" in text
    assert "1-40 of 60" in text


def test_paging_across_reaches_the_columns_past_the_cut(wide):
    r = O.read_sheet(wide, start_col=41, max_cols=40)
    assert "S41" in r["grid"]
    assert "S60" in r["grid"]


def test_the_column_letters_follow_the_page(wide):
    r = O.read_sheet(wide, start_col=41, max_cols=2)
    first = [ln for ln in r["grid"].splitlines() if "S41" in ln][0]
    header = r["grid"].splitlines()[0]
    assert "AO" in header, "column 41 is AO, not A"
    assert first


def test_a_narrow_table_needs_no_paging_note(tmp_path):
    p = tmp_path / "schmal.csv"
    p.write_text("a,b\n1,2\n", encoding="utf-8")
    assert "start_col" not in _notes(O.read_sheet(p))


def test_the_tool_offers_start_col():
    """A remedy the model cannot reach is not a remedy."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    i = src.index('"start_row": {')
    assert '"start_col"' in src[i:i + 1200]
    assert 'start_col=_as_int(arguments.get("start_col"), 1)' in src


# ---------------------------------------------------------------------------
# The profile says which rows it describes
# ---------------------------------------------------------------------------

def test_a_profile_over_part_of_the_file_says_so(tmp_path):
    """The CSV reader profiles up to 2000 rows, NOT the display window —
    so for a table shorter than that the profile really does describe the
    whole file, whatever the window shows. The claim only bites past that
    line. Asserted here at the real boundary rather than at the one the
    defect report assumed."""
    p = tmp_path / "lang.csv"
    lines = ["Beleg;Betrag"] + [f"R-{i:05d};{i},50" for i in range(1, 2500)]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    r = O.read_sheet(p, max_rows=100)
    assert r["profile_rows"] == 2000
    text = _notes(r)
    assert "rows 1-2000 of 2500" in text
    assert "not the whole file" in text


def test_a_window_shorter_than_the_profile_is_not_flagged(tmp_path):
    """A 400-row file read 100 rows at a time is still profiled in full;
    saying otherwise would be a caveat about nothing."""
    p = tmp_path / "mittel.csv"
    lines = ["Beleg;Betrag"] + [f"R-{i:04d};{i},50" for i in range(1, 400)]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    r = O.read_sheet(p, max_rows=100)
    assert r["profile_rows"] == r["rows"]
    assert "not the whole file" not in _notes(r)


def test_a_profile_over_the_whole_file_says_nothing(tmp_path):
    p = tmp_path / "kurz.csv"
    p.write_text("Beleg;Betrag\nR-1;1,50\nR-2;2,50\n", encoding="utf-8")
    r = O.read_sheet(p)
    assert r["profile_rows"] == r["rows"]
    assert "not the whole file" not in _notes(r)


def test_the_workbook_reader_reports_its_profile_scope():
    openpyxl = pytest.importorskip("openpyxl")
    from delfin.agent import benchmark_fixtures as BF
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[1]
    if not (root / BF.OFFICE_WS_REL / "Verbrauch_2026.xlsx").is_file():
        BF.ensure_office_fixtures(root)
    r = O.read_sheet(root / BF.OFFICE_WS_REL / "Verbrauch_2026.xlsx")
    assert r["profile_rows"] == 200
    assert r["rows"] > 200
    assert "not the whole file" in _notes(r)
