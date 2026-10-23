"""Totalling a column says what it could not add.

The defect: nothing in the catalogue totalled anything. A sum was
arithmetic the model did in its head over a rendered grid, so the number
in the answer had no tool result behind it — nothing recorded which rows
it left out, and nothing downstream could tell a total that dropped a row
from one that did not. The fixture workbooks are built out of exactly
those rows: "1.500,00-" and "(340,00)" are credit notes, "n/a" is a
missing amount, seven rows are hidden behind a filter, and every one of
them reads as an ordinary part of the column.

What is pinned here: every row of a column lands in exactly one of
counted / skipped / blank, the three add up to the rows of the table, the
values that were skipped are named, an undecidable column is refused
rather than guessed, and the figures come back in a shape the answer-side
ledger consumes without re-reading the report.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import office
from delfin.agent.api_client import (
    _DOC_TOOLS_OPENAI,
    _OFFICE_TOOL_NAMES,
    _DocToolExecutor,
    KitToolPermissions,
    ToolSurfaceContext,
    advertisable_tools,
    tool_unavailable_reason,
)

openpyxl = pytest.importorskip("openpyxl")


@pytest.fixture(autouse=True)
def clean_ledger():
    office.reset_figure_ledger()
    yield
    office.reset_figure_ledger()


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    return d


@pytest.fixture
def posten(ws):
    """Open items the way a German ledger writes them: decimal comma, two
    credit-note conventions in one column, one row without an amount, and
    a filter that hides a row."""
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.title = "Offene Posten"
    for row in (
        ("Beleg", "Kunde", "Betrag"),
        ("A-001", "Institut", "4.820,00"),
        ("A-002", "Werkstoff", "1.240,50"),
        ("G-001", "Institut", "1.500,00-"),
        ("A-003", "Analytik", "962,40"),
        ("G-002", "Werkstoff", "(340,00)"),
        ("A-004", "Institut", ""),
    ):
        sheet.append(row)
    sheet.row_dimensions[3].hidden = True
    path = ws / "posten.xlsx"
    wb.save(path)
    wb.close()
    return path


def _perms(ws: Path) -> KitToolPermissions:
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "sum-column-test"
    return perms


# ---------------------------------------------------------------------------
# Coverage: every row lands somewhere, and the somewhere is named
# ---------------------------------------------------------------------------

def test_every_row_is_counted_skipped_or_empty(posten):
    result = office.sum_column(posten, "Betrag")
    assert (result["counted"] + len(result["skipped"]) + result["blank"]
            == result["rows"] == 6)


def test_the_values_it_could_not_add_are_named(posten):
    result = office.sum_column(posten, "Betrag")
    assert result["skipped"] == ["1.500,00-", "(340,00)"]
    note = " ".join(result["notes"])
    assert "NOT in the total" in note
    assert "1.500,00-" in note


def test_a_decimal_comma_column_totals_as_money(posten):
    """Read as-is these values are off by a factor of a thousand."""
    result = office.sum_column(posten, "Betrag")
    assert result["total"] == pytest.approx(4820.00 + 1240.50 + 962.40)
    assert result["convention"] == office.DECIMAL_COMMA


def test_hidden_rows_are_in_the_total_and_said_so(posten):
    result = office.sum_column(posten, "Betrag")
    assert any("hidden or filtered" in n for n in result["notes"])
    assert any(f["value"] == 1.0 and "hidden" in f["label"]
               for f in result["figures"])


def test_an_undecidable_column_is_refused_with_the_way_out(ws):
    table = ws / "inventar.csv"
    table.write_text("Gerät;Wert\nWaage;8.986\nOfen;1.250\n", encoding="utf-8")
    with pytest.raises(office.OfficeError) as excinfo:
        office.sum_column(table, "Wert")
    message = str(excinfo.value)
    assert "8.986" in message
    assert "decimal_comma" in message
    # And the way out works, with the source of the decision recorded:
    # under the German reading the dot groups thousands, so the same
    # column totals 10236 rather than 10,236.
    result = office.sum_column(table, "Wert", convention="decimal_comma")
    assert result["total"] == pytest.approx(10236.0)
    assert office.sum_column(
        table, "Wert", convention="decimal_point")["total"] == pytest.approx(
            10.236)
    assert any("given by the caller" in n for n in result["notes"])


def test_a_missing_column_names_the_columns_that_are_there(posten):
    with pytest.raises(office.OfficeError) as excinfo:
        office.sum_column(posten, "Summe")
    assert "Beleg" in str(excinfo.value)


def test_a_title_row_above_the_header_is_handled(ws):
    wb = openpyxl.Workbook()
    sheet = wb.active
    for row in (("Kostenstellenübersicht 2026",),
                ("Kostenstelle", "Budget"),
                ("4711", 18500), ("4712", 5400)):
        sheet.append(row)
    path = ws / "kostenstellen.xlsx"
    wb.save(path)
    wb.close()
    result = office.sum_column(path, "Budget", header_row=2)
    assert result["total"] == 23900
    assert result["counted"] == 2


def test_group_totals_add_up_to_the_total(posten):
    result = office.sum_column(posten, "Betrag", group_by="Kunde")
    assert sum(g["total"] for g in result["groups"]) == pytest.approx(
        result["total"])
    assert {g["key"] for g in result["groups"]} == {"Institut", "Werkstoff",
                                                    "Analytik"}


# ---------------------------------------------------------------------------
# The figures come back as figures, not as prose to be re-read
# ---------------------------------------------------------------------------

def test_the_total_is_returned_in_a_shape_the_ledger_can_consume(posten):
    result = office.sum_column(posten, "Betrag")
    figures = result["figures"]
    assert all({"kind", "value", "label"} <= set(f) for f in figures)
    assert any(f["kind"] == "sum" and f["value"] == result["total"]
               for f in figures)
    assert any(f["kind"] == "count" and f["value"] == result["counted"]
               for f in figures)
    # And the ledger takes them without touching the rendered report.
    recorded = office.record_tool_figures("sum_column", result)
    assert [f.value for f in recorded] == [f["value"] for f in figures]


def test_a_comparison_returns_its_counts_the_same_way(ws):
    left = ws / "a.csv"
    right = ws / "b.csv"
    left.write_text("Beleg;Betrag\nR-1;10,00\nR-2;20,00\n", encoding="utf-8")
    right.write_text("Beleg;Betrag\nR-1;10,00\nR-3;30,00\n", encoding="utf-8")
    result = office.compare_tables(left, right, key="Beleg")
    values = {f["label"]: f["value"] for f in result["figures"]}
    assert values[f"{result['equal_count']} equal row(s)"] == 1
    assert any("only in the left table" in label for label in values)


# ---------------------------------------------------------------------------
# The tool: dispatch, refusals, and what it puts in the ledger
# ---------------------------------------------------------------------------

def test_the_tool_reports_the_total_with_what_it_left_out(ws, posten):
    ex = _DocToolExecutor()
    out = ex._dispatch("sum_column",
                       {"path": str(posten), "column": "Betrag"},
                       _perms(ws))
    assert "7022.9" in out
    assert "1.500,00-" in out
    assert "counted 3 of 6" in out
    assert "hidden or filtered" in out


def test_the_tool_call_lands_in_the_figure_ledger(ws, posten):
    ex = _DocToolExecutor()
    ex._dispatch("sum_column", {"path": str(posten), "column": "Betrag"},
                 _perms(ws))
    assert any(f.kind == "sum" and f.value == pytest.approx(7022.9)
               for f in office.figure_ledger())
    # Which is what makes the answer that states it groundable.
    assert office.scan_answer_for_unledgered_figures(
        "Die Summe der offenen Posten beträgt 7.022,90 EUR.") == []


def test_a_call_without_a_column_says_what_is_missing(ws, posten):
    ex = _DocToolExecutor()
    out = json.loads(ex._dispatch("sum_column", {"path": str(posten)},
                                  _perms(ws)))
    assert "column is required" in out["error"]


def test_a_read_outside_the_workspace_is_refused(ws, tmp_path, posten):
    other = tmp_path / "outside"
    other.mkdir()
    ex = _DocToolExecutor()
    out = ex._dispatch("sum_column",
                       {"path": str(posten), "column": "Betrag"},
                       _perms(other))
    assert "error" in json.loads(out)


def test_a_reading_the_tool_refuses_is_reported_and_not_guessed(ws):
    table = ws / "inventar.csv"
    table.write_text("Gerät;Wert\nWaage;8.986\nOfen;1.250\n", encoding="utf-8")
    ex = _DocToolExecutor()
    out = json.loads(ex._dispatch(
        "sum_column", {"path": str(table), "column": "Wert"}, _perms(ws)))
    assert "8.986" in out["error"]
    assert office.figure_ledger() == []


# ---------------------------------------------------------------------------
# Surface registration
# ---------------------------------------------------------------------------

def test_the_tool_is_in_the_catalogue_and_the_office_family():
    names = {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}
    assert "sum_column" in names
    assert "sum_column" in _OFFICE_TOOL_NAMES


def test_it_is_dropped_without_a_spreadsheet_backend():
    assert tool_unavailable_reason(
        "sum_column", ToolSurfaceContext(has_office_libs=False)) is not None
    advertised = {t["function"]["name"] for t in advertisable_tools(
        _DOC_TOOLS_OPENAI, ToolSurfaceContext(has_office_libs=True))}
    assert "sum_column" in advertised


def test_it_is_read_only_and_reaches_the_office_role():
    from delfin.agent.api_client import _PLAN_READONLY_TOOLS, _tool_denied_for_role

    assert "sum_column" in _PLAN_READONLY_TOOLS
    assert not _tool_denied_for_role("office_agent", "sum_column")
