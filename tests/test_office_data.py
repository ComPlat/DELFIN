"""Numbers, dates and reconciliation — the arithmetic that fails quietly.

A cell reading "1.234,50" is a string. Handed to arithmetic it either
raises or, worse, parses as 1.234 — off by a factor of a thousand with
no error anywhere. Dates carry the same trap: 03.04.2026 and 04/03/2026
are the same day under two conventions and different days under one.

The decision is made per COLUMN. A single "1.234" is genuinely
ambiguous; a column holding it next to one value with both separators is
not. Where the column stays ambiguous the answer is to say so, not to
pick the reading that looks plausible.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import office
from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions

openpyxl = pytest.importorskip("openpyxl")


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    return d


def _perms(ws) -> KitToolPermissions:
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "office-data"
    return perms


# ---------------------------------------------------------------------------
# Number conventions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value,convention,expected", [
    ("1.234,50", office.DECIMAL_COMMA, 1234.5),
    ("1.234.567,89", office.DECIMAL_COMMA, 1234567.89),
    ("-1.234,50", office.DECIMAL_COMMA, -1234.5),
    ("0,5", office.DECIMAL_COMMA, 0.5),
    ("1,234.50", office.DECIMAL_POINT, 1234.5),
    ("1234.50", office.DECIMAL_POINT, 1234.5),
    ("12345", office.PLAIN_NUMBER, 12345.0),
    ("1.234,50 EUR", office.DECIMAL_COMMA, 1234.5),
    ("€ 89,90", office.DECIMAL_COMMA, 89.9),
])
def test_numbers_parse_under_their_convention(value, convention, expected):
    assert office.parse_number(value, convention) == pytest.approx(expected)


@pytest.mark.parametrize("value", ["n/a", "", "abc", "-", "1.2.3.x", None])
def test_non_numbers_return_none(value):
    assert office.parse_number(value, office.DECIMAL_COMMA) is None


def test_a_value_with_both_separators_settles_the_column():
    convention, why = office.detect_number_convention(
        ["1.234,50", "89,90", "12,00"])
    assert convention == office.DECIMAL_COMMA
    assert "settle" in why
    convention, _ = office.detect_number_convention(["1,234.50", "89.90"])
    assert convention == office.DECIMAL_POINT


def test_a_non_three_digit_group_settles_the_column():
    """"89,90" cannot be a thousands separator — the tail is two digits."""
    assert office.detect_number_convention(["89,90"])[0] == office.DECIMAL_COMMA
    assert office.detect_number_convention(["89.90"])[0] == office.DECIMAL_POINT


def test_repeated_separators_settle_the_column():
    assert office.detect_number_convention(
        ["1.234.567"])[0] == office.DECIMAL_COMMA
    assert office.detect_number_convention(
        ["1,234,567"])[0] == office.DECIMAL_POINT


def test_an_undecidable_column_is_reported_not_guessed():
    convention, why = office.detect_number_convention(["1.234", "2.500"])
    assert convention == office.AMBIGUOUS
    assert "two different numbers" in why


def test_a_column_mixing_conventions_is_ambiguous():
    convention, why = office.detect_number_convention(["1.234,50", "1,234.50"])
    assert convention == office.AMBIGUOUS
    assert "mixes" in why


def test_a_column_without_separators_needs_no_convention():
    assert office.detect_number_convention(["1", "22", "333"])[0] == \
        office.PLAIN_NUMBER


# ---------------------------------------------------------------------------
# Date conventions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value,convention,expected", [
    ("31.07.2026", office.DAY_FIRST, "2026-07-31"),
    ("1.8.2026", office.DAY_FIRST, "2026-08-01"),
    ("07/31/2026", office.MONTH_FIRST, "2026-07-31"),
    ("2026-07-31", office.ISO_DATE, "2026-07-31"),
    ("31.07.26", office.DAY_FIRST, "2026-07-31"),
])
def test_dates_parse_under_their_convention(value, convention, expected):
    assert office.parse_date(value, convention) == expected


def test_an_impossible_date_is_rejected():
    assert office.parse_date("31.02.2026", office.DAY_FIRST) is None
    assert office.parse_date("nicht am 1.", office.DAY_FIRST) is None


def test_a_day_above_twelve_settles_the_column():
    convention, why = office.detect_date_convention(
        ["01/02/2026", "31/07/2026"])
    assert convention == office.DAY_FIRST
    assert "above 12" in why
    assert office.detect_date_convention(
        ["01/02/2026", "07/31/2026"])[0] == office.MONTH_FIRST


def test_dotted_dates_are_day_first():
    assert office.detect_date_convention(["01.02.2026"])[0] == office.DAY_FIRST


def test_iso_dates_are_recognised():
    assert office.detect_date_convention(["2026-07-31"])[0] == office.ISO_DATE


def test_a_column_that_could_be_read_either_way_is_ambiguous():
    convention, why = office.detect_date_convention(
        ["01/02/2026", "03/04/2026"])
    assert convention == office.AMBIGUOUS
    assert "either way" in why


def test_a_column_mixing_date_conventions_is_ambiguous():
    assert office.detect_date_convention(
        ["31/07/2026", "07/31/2026"])[0] == office.AMBIGUOUS


# ---------------------------------------------------------------------------
# Column profiles reach the reader
# ---------------------------------------------------------------------------

_GERMAN = ("Name;Datum;Betrag;Anzahl\n"
           "Müller;31.07.2026;1.234,50;3\n"
           "Özdemir;01.08.2026;89,90;12\n"
           "Schmidt;12.09.2026;12.000,00;7\n"
           "Weber;02.10.2026;n/a;4\n")


def test_reading_a_table_reports_kinds_and_conventions(ws):
    p = ws / "kosten.csv"
    p.write_bytes(_GERMAN.encode("cp1252"))
    result = office.read_sheet(p)
    by_name = {c["name"]: c for c in result["column_profile"]}
    assert by_name["Name"]["kind"] == "text"
    assert by_name["Datum"]["kind"] == "date"
    assert by_name["Datum"]["convention"] == office.DAY_FIRST
    assert by_name["Betrag"]["kind"] == "number"
    assert by_name["Betrag"]["convention"] == office.DECIMAL_COMMA
    assert by_name["Anzahl"]["convention"] == office.PLAIN_NUMBER


def test_the_decimal_comma_warning_reaches_the_notes(ws):
    p = ws / "kosten.csv"
    p.write_bytes(_GERMAN.encode("cp1252"))
    notes = " ".join(office.read_sheet(p)["notes"])
    assert "decimal comma" in notes
    assert "factor of a thousand" in notes


def test_values_that_would_be_silently_dropped_are_named(ws):
    p = ws / "kosten.csv"
    p.write_bytes(_GERMAN.encode("cp1252"))
    notes = " ".join(office.read_sheet(p)["notes"])
    assert "n/a" in notes
    assert "leave them out" in notes


def test_both_findings_are_reported_for_the_same_column(ws):
    """A column can use a decimal comma AND hold unparseable values; each
    on its own is enough to make a naive total wrong."""
    p = ws / "kosten.csv"
    p.write_bytes(_GERMAN.encode("cp1252"))
    notes = office.read_sheet(p)["notes"]
    assert any("decimal comma" in n for n in notes)
    assert any("not numbers" in n for n in notes)


def test_an_ambiguous_column_tells_the_reader_to_ask(ws):
    p = ws / "amb.csv"
    p.write_text("Posten,Wert\nA,1.234\nB,2.500\n", encoding="utf-8")
    notes = " ".join(office.read_sheet(p)["notes"])
    assert "Ask which reading is meant" in notes


def test_the_profile_is_shown_next_to_the_grid(ws):
    p = ws / "kosten.csv"
    p.write_bytes(_GERMAN.encode("cp1252"))
    out = _DocToolExecutor()._dispatch(
        "read_document", {"path": str(p)}, _perms(ws))
    assert "Betrag: number (decimal_comma)" in out
    assert "Datum: date (day_first)" in out
    assert "not parseable" in out


def test_paging_does_not_promote_a_data_row_to_a_header(ws):
    wb = openpyxl.Workbook()
    wb.active.append(["Posten", "Betrag"])
    for n in range(1, 40):
        wb.active.append([f"P{n}", n])
    wb.save(ws / "lang.xlsx")
    wb.close()
    result = office.read_sheet(ws / "lang.xlsx", start_row=10, max_rows=5)
    names = [c["name"] for c in result["column_profile"]]
    assert names == ["A", "B"]        # positional, not the row that was read


# ---------------------------------------------------------------------------
# Reconciliation
# ---------------------------------------------------------------------------

_BOOKINGS = ("Beleg;Name;Betrag;Datum\n"
             "R-001;Müller;1.234,50;31.07.2026\n"
             "R-002;Özdemir;89,90;01.08.2026\n"
             "R-003;Schmidt;450,00;02.08.2026\n"
             "R-005;Weber;12,00;03.08.2026\n")
_INVOICES = ("Beleg,Name,Betrag,Datum\n"
             "R-001,Müller,1234.50,2026-07-31\n"
             "R-002,Özdemir,98.90,2026-08-01\n"
             "R-003,Schmidt,450.00,2026-08-02\n"
             "R-004,Neu,77.00,2026-08-06\n")


@pytest.fixture
def pair(ws):
    left = ws / "buchungen.csv"
    right = ws / "rechnungen.csv"
    left.write_bytes(_BOOKINGS.encode("cp1252"))
    right.write_text(_INVOICES, encoding="utf-8")
    return left, right


def test_the_same_value_written_two_ways_is_not_a_difference(pair):
    """The two tables come from different systems. 1.234,50 and 1234.50
    are the same amount, and 31.07.2026 is 2026-07-31."""
    left, right = pair
    result = office.compare_tables(left, right, key="Beleg")
    assert "R-001" in result["equal"]
    assert "R-003" in result["equal"]


def test_a_real_difference_is_caught(pair):
    left, right = pair
    result = office.compare_tables(left, right, key="Beleg")
    keys = {d["key"] for d in result["differing"]}
    assert keys == {"R-002"}
    diff = result["differing"][0]["differences"][0]
    assert diff["column"] == "Betrag"
    assert diff["left"] == "89,90" and diff["right"] == "98.90"


def test_one_sided_rows_are_reported_per_side(pair):
    left, right = pair
    result = office.compare_tables(left, right, key="Beleg")
    assert result["only_left"] == ["R-005"]
    assert result["only_right"] == ["R-004"]


def test_every_input_row_is_accounted_for(pair):
    left, right = pair
    result = office.compare_tables(left, right, key="Beleg")
    assert result["rows_accounted_for"] is True
    covered = (result["equal_count"] + result["differing_count"]
               + result["only_left_count"] + result["only_right_count"])
    assert covered == 5          # 3 matched + 1 left-only + 1 right-only


def test_a_duplicate_key_is_reported_instead_of_joined(ws):
    """A duplicate key turns the join into a cross product, and the result
    looks perfectly plausible."""
    left = ws / "l.csv"
    right = ws / "r.csv"
    left.write_text("K,V\nA,1\nA,2\nB,3\n", encoding="utf-8")
    right.write_text("K,V\nA,1\nB,3\n", encoding="utf-8")
    result = office.compare_tables(left, right, key="K")
    assert "A" not in result["equal"]
    reasons = [e["reason"] for e in result["not_comparable"]]
    assert any("more than once" in r for r in reasons)
    assert any("cross product" in n for n in result["notes"])
    assert result["equal"] == ["B"]


def test_an_empty_key_is_reported_not_dropped(ws):
    left = ws / "l.csv"
    right = ws / "r.csv"
    left.write_text("K,V\nA,1\n,2\n", encoding="utf-8")
    right.write_text("K,V\nA,1\n", encoding="utf-8")
    result = office.compare_tables(left, right, key="K")
    assert any(e["reason"] == "empty key" for e in result["not_comparable"])
    assert result["rows_accounted_for"] is True


def test_an_unknown_key_column_lists_the_real_ones(pair):
    left, right = pair
    with pytest.raises(office.OfficeError) as exc:
        office.compare_tables(left, right, key="Belegnummer")
    assert "Beleg" in str(exc.value)


def test_only_the_named_columns_are_compared(pair):
    left, right = pair
    result = office.compare_tables(
        left, right, key="Beleg", columns=["Name"])
    assert result["compared_columns"] == ["Name"]
    assert result["differing_count"] == 0     # the names all agree


def test_tables_can_be_compared_across_formats(ws, pair):
    left, _ = pair
    wb = openpyxl.Workbook()
    wb.active.append(["Beleg", "Betrag"])
    wb.active.append(["R-001", 1234.5])
    wb.active.append(["R-002", 89.9])
    xlsx = ws / "aus_system.xlsx"
    wb.save(xlsx)
    wb.close()
    result = office.compare_tables(left, xlsx, key="Beleg")
    assert result["compared_columns"] == ["Betrag"]
    assert set(result["equal"]) == {"R-001", "R-002"}


def test_comparison_through_the_tool_reports_every_group(ws, pair):
    left, right = pair
    out = _DocToolExecutor()._dispatch("compare_tables", {
        "left": str(left), "right": str(right), "key": "Beleg",
    }, _perms(ws))
    for label in ("equal:", "differing:", "only left:", "only right:",
                  "not comparable:"):
        assert label in out
    assert "R-002" in out
    assert "89,90 | 98.90" in out


def test_the_tool_requires_a_key(ws, pair):
    left, right = pair
    out = json.loads(_DocToolExecutor()._dispatch("compare_tables", {
        "left": str(left), "right": str(right),
    }, _perms(ws)))
    assert "key is required" in out["error"]


def test_comparing_something_that_is_not_a_table_is_refused(ws, pair):
    left, _ = pair
    doc = ws / "brief.docx"
    doc.write_bytes(b"PK\x03\x04not really")
    out = json.loads(_DocToolExecutor()._dispatch("compare_tables", {
        "left": str(left), "right": str(doc), "key": "Beleg",
    }, _perms(ws)))
    assert "not a table" in out["error"]
