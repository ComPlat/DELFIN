"""The read side spoke German and the write side did not.

``detect_number_convention``, ``parse_number`` and the date readers all
handle decimal commas and day-first dates. Everything that puts a value
back into a document went through ``str(value)``.

Proven end to end before the fix: a row ``["Müller GmbH", 1234.5,
datetime(2026, 7, 31)]``, with the cells formatted ``#,##0.00 "€"`` and
``DD.MM.YYYY``, run through ``fill_series`` into a letter template
produced

    wir stellen Ihnen 1234.5 EUR in Rechnung, fällig am
    2026-07-31 00:00:00.
    counts: {'ok': 1, 'incomplete': 0, 'failed': 0}

— the en-US number, the ISO timestamp, and a run that reported success.
The cell's own ``number_format`` was available at every step and read at
none. This is a defect in the artefact that leaves the building and goes
to a customer.

The tests below fix the two values a German document is judged by
(``1.234,50`` and ``31.07.2026``) and, just as important, fix what must
NOT change: a string stays the string the user typed.
"""

from __future__ import annotations

import datetime

import pytest

from delfin.agent import german as G
from delfin.agent import office as O

openpyxl = pytest.importorskip("openpyxl")
docx = pytest.importorskip("docx")


# ---------------------------------------------------------------------------
# fixtures: one customer row, one letter template
# ---------------------------------------------------------------------------

@pytest.fixture
def kunden(tmp_path):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["kunde", "betrag", "faellig"])
    ws.append(["Müller GmbH", 1234.5, datetime.datetime(2026, 7, 31)])
    ws.cell(row=2, column=2).number_format = '#,##0.00 "€"'
    ws.cell(row=2, column=3).number_format = "DD.MM.YYYY"
    p = tmp_path / "kunden.xlsx"
    wb.save(p)
    return p


@pytest.fixture
def brief(tmp_path):
    doc = docx.Document()
    doc.add_paragraph("Sehr geehrte Damen und Herren von {{kunde}},")
    doc.add_paragraph(
        "wir stellen Ihnen {{betrag}} EUR in Rechnung, "
        "fällig am {{faellig}}.")
    p = tmp_path / "brief.docx"
    doc.save(p)
    return p


def _letter_text(path) -> str:
    return "\n".join(p.text for p in docx.Document(str(path)).paragraphs)


# ---------------------------------------------------------------------------
# the round trip that shipped wrong
# ---------------------------------------------------------------------------

def test_a_serienbrief_shows_the_amount_the_way_the_cell_does(
        tmp_path, kunden, brief):
    result = O.fill_series(kunden, brief, output_dir=tmp_path / "out",
                           name_pattern="brief_{row}.docx")
    assert result["counts"] == {"ok": 1, "incomplete": 0, "failed": 0}
    text = _letter_text(tmp_path / "out" / "brief_2.docx")
    assert "1.234,50 EUR" in text, text
    assert "1234.5" not in text


def test_a_serienbrief_shows_a_date_and_not_a_timestamp(
        tmp_path, kunden, brief):
    O.fill_series(kunden, brief, output_dir=tmp_path / "out",
                  name_pattern="brief_{row}.docx")
    text = _letter_text(tmp_path / "out" / "brief_2.docx")
    assert "31.07.2026" in text, text
    assert "2026-07-31" not in text
    assert "00:00:00" not in text


def test_the_umlaut_reaches_the_customer_unchanged(tmp_path, kunden, brief):
    O.fill_series(kunden, brief, output_dir=tmp_path / "out",
                  name_pattern="brief_{row}.docx")
    assert "Müller GmbH" in _letter_text(tmp_path / "out" / "brief_2.docx")


# ---------------------------------------------------------------------------
# the fillers themselves, called directly
# ---------------------------------------------------------------------------

def test_the_docx_filler_renders_a_value_it_is_handed_directly(
        tmp_path, brief):
    out = tmp_path / "direkt.docx"
    O.fill_docx_template(
        brief,
        {"kunde": "Özdemir KG", "betrag": 89.9,
         "faellig": datetime.date(2026, 8, 1)},
        output=out)
    text = _letter_text(out)
    assert "89,9 EUR" in text, text
    assert "01.08.2026" in text, text


def test_a_string_is_never_reformatted(tmp_path, brief):
    """A CSV already holds the text the user typed. Re-rendering it would
    be the same defect in reverse."""
    out = tmp_path / "text.docx"
    O.fill_docx_template(
        brief,
        {"kunde": "Meier", "betrag": "1.234,50", "faellig": "31.07.2026"},
        output=out)
    text = _letter_text(out)
    assert "1.234,50 EUR" in text
    assert "31.07.2026" in text


def test_a_csv_series_still_passes_its_text_through(tmp_path, brief):
    table = tmp_path / "kunden.csv"
    table.write_text("kunde;betrag;faellig\nMüller;1.234,50;31.07.2026\n",
                     encoding="utf-8")
    O.fill_series(table, brief, output_dir=tmp_path / "out",
                  name_pattern="csv_{row}.docx")
    text = _letter_text(tmp_path / "out" / "csv_2.docx")
    assert "1.234,50 EUR" in text
    assert "31.07.2026" in text


# ---------------------------------------------------------------------------
# the renderer, on its own
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value,fmt,expected", [
    (1234.5, '#,##0.00 "€"', "1.234,50"),
    (1234.5, "#,##0.00", "1.234,50"),
    (1234.5, "0.00", "1234,50"),
    (1234.5, "General", "1234,5"),
    (1234.5, "", "1234,5"),
    (-1234567.891, "#,##0.00", "-1.234.567,89"),
    (0.155, "0.0%", "15,5 %"),
    (1234, "#,##0", "1.234"),
    (42, "", "42"),
    (datetime.datetime(2026, 7, 31), "DD.MM.YYYY", "31.07.2026"),
    (datetime.datetime(2026, 7, 31), "", "31.07.2026"),
    (datetime.datetime(2026, 7, 31, 14, 5), "", "31.07.2026 14:05"),
    (datetime.date(2026, 7, 31), "", "31.07.2026"),
    ("1.234,50", '#,##0.00 "€"', "1.234,50"),
    ("R-001", "", "R-001"),
    (None, "", ""),
])
def test_the_renderer_reads_the_cells_own_format(value, fmt, expected):
    assert O._render_for_locale(value, fmt) == expected


def test_a_currency_literal_in_the_format_is_not_repeated():
    """The template already says EUR. Appending the format's own "€" would
    hand the customer "1.234,50 € EUR"."""
    assert "€" not in O._render_for_locale(1234.5, '#,##0.00 "€"')


def test_a_date_format_the_cell_states_wins_over_the_default():
    """A user who formatted the column ISO chose ISO."""
    assert O._render_for_locale(
        datetime.date(2026, 7, 31), "YYYY-MM-DD") == "2026-07-31"


def test_a_number_format_is_never_mistaken_for_a_date():
    assert G.excel_date_pattern('#,##0.00 "€"') == ""
    assert G.excel_date_pattern("General") == ""
    assert G.excel_date_pattern("DD.MM.YYYY") == "%d.%m.%Y"
