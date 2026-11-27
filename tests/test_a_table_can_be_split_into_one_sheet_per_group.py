"""A workbook can gain a sheet per group through edit_sheet.

Asked for "one sheet per collaborator in column E", the agent read the
workbook, found that edit_sheet sets cells and appends rows but cannot add
a sheet to an existing file, and ended its turns without a word. A sheet
name the workbook does not hold, sent with only append_rows, now adds that
sheet; every other shape keeps refusing a missing name, because there it
is a mistake rather than a request.

Driven through the tool entry (_DocToolExecutor.execute), the same route a
model's call takes.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

openpyxl = pytest.importorskip("openpyxl")


@pytest.fixture
def book(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", lambda: home)
    ws = tmp_path / "ws"
    ws.mkdir()
    p = ws / "platten.xlsx"
    wb = openpyxl.Workbook()
    sh = wb.active
    sh.title = "Screening"
    sh.append(["Platte", "Datum", "Ergebnis", "Notiz", "Collaborator"])
    sh.append(["P1", "2026-11-01", 0.4, "ok", "Anna"])
    sh.append(["P2", "2026-11-02", 0.7, "ok", "Ben"])
    sh.append(["P3", "2026-11-03", 0.5, "ok", "Anna"])
    wb.save(p)
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "sheet-per-group"
    return p, perms


def _call(perms, **args) -> dict:
    ex = _DocToolExecutor()
    ex.execute("read_document", {"path": args["path"]}, perms)
    out = ex.execute("edit_sheet", args, perms)
    start = out.find("{")
    return json.loads(out[start:]) if start >= 0 else {"raw": out}


def test_append_rows_under_a_new_name_adds_that_sheet(book):
    p, perms = book
    rows = [["Platte", "Datum", "Ergebnis", "Notiz", "Collaborator"],
            ["P1", "2026-11-01", 0.4, "ok", "Anna"],
            ["P3", "2026-11-03", 0.5, "ok", "Anna"]]
    out = _call(perms, path=str(p), sheet="Anna", append_rows=rows)
    assert "error" not in out, out
    wb = openpyxl.load_workbook(p)
    assert wb.sheetnames == ["Screening", "Anna"]
    assert [list(r) for r in wb["Anna"].values] == rows
    # The source sheet is untouched.
    assert wb["Screening"]["E3"].value == "Ben"
    assert wb["Screening"].max_row == 4


def test_a_missing_sheet_with_cell_edits_is_still_refused(book):
    p, perms = book
    out = _call(perms, path=str(p), sheet="Ben",
                edits=[{"cell": "A1", "value": "x"}])
    text = json.dumps(out)
    assert "no sheet named 'Ben'" in text
    assert "append_rows" in text          # says how a sheet is added
    assert openpyxl.load_workbook(p).sheetnames == ["Screening"]


def test_a_name_a_sheet_cannot_carry_is_refused(book):
    p, perms = book
    out = _call(perms, path=str(p), sheet="A/B", append_rows=[["x"]])
    assert "cannot be a sheet name" in json.dumps(out)
    assert openpyxl.load_workbook(p).sheetnames == ["Screening"]


def test_an_existing_sheet_still_appends_below_its_rows(book):
    p, perms = book
    out = _call(perms, path=str(p), sheet="Screening",
                append_rows=[["P4", "2026-11-04", 0.9, "ok", "Ben"]])
    assert "error" not in out, out
    sh = openpyxl.load_workbook(p)["Screening"]
    assert sh.max_row == 5 and sh["A5"].value == "P4"
