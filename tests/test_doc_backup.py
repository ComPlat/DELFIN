"""Version history for documents the agent changes.

The file browser copies a document aside before it saves an edit. The
agent has to write into the SAME history — same folder, same names —
because two mechanisms keeping two histories of one file is worse than
one: the user opens the folder they know and sees half the story.
"""

from __future__ import annotations

import json

import pathlib

import pytest

from delfin import doc_backup as bk
from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions

openpyxl = pytest.importorskip("openpyxl")


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "OFFICE"
    d.mkdir()
    return d


@pytest.fixture
def book(ws):
    wb = openpyxl.Workbook()
    wb.active.append(["Beleg", "Status"])
    wb.active.append(["R-001", "offen"])
    path = ws / "vorgaenge.xlsx"
    wb.save(path)
    wb.close()
    return path


def _perms(ws):
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "backup-test"
    return perms


# ---------------------------------------------------------------------------
# The naming rules, which belong to the file browser
# ---------------------------------------------------------------------------

def test_copies_go_into_a_folder_beside_the_document(book):
    made = bk.make_backup(book)
    assert made.parent == book.parent / "Backups"
    assert made.name == "vorgaenge.bak.xlsx"


def test_each_save_keeps_its_own_copy(book):
    """Without numbering the second save finds a backup already there and
    keeps nothing — a safety net for one edit, not a history."""
    names = []
    for _ in range(3):
        made = bk.make_backup(book)
        names.append(made.name)
        book.write_bytes(book.read_bytes() + b" ")
    assert names == ["vorgaenge.bak.xlsx", "vorgaenge.bak2.xlsx",
                     "vorgaenge.bak3.xlsx"]


def test_the_first_copy_carries_no_number(book):
    """Matching what the browser writes, so both fill one sequence
    instead of starting two."""
    assert bk.make_backup(book).name.endswith(".bak.xlsx")


def test_the_file_browser_writes_this_history_and_not_its_own(book):
    """The browser used to carry a copy of this numbering. It now asks for
    a versioned backup and gets one from here, so the two cannot drift --
    there is only one of them. What is pinned is that it still delegates:
    a local implementation reappearing here is what would split the folder
    the user opens into half a history each."""
    sheet_view = pytest.importorskip("delfin.dashboard.spreadsheet_view")
    source = pathlib.Path(sheet_view.__file__).read_text(encoding="utf-8")
    assert "def versioned_backup_path" not in source, (
        "the browser grew its own numbering again"
    )

    folder = book.parent / "Backups"
    made = sheet_view.make_backup(book, folder=folder, versioned=True)
    assert made == bk.versioned_backup_path(book, folder).with_name(made.name)
    assert made.name == "vorgaenge.bak.xlsx"

    book.write_bytes(b"changed")
    assert bk.make_backup(book, folder=folder).name == "vorgaenge.bak2.xlsx", (
        "the two started separate sequences"
    )


def test_the_browser_still_keeps_one_copy_beside_a_calculation(book):
    """Its own users have always found it there, and that half is not
    this module's business."""
    sheet_view = pytest.importorskip("delfin.dashboard.spreadsheet_view")
    assert sheet_view.backup_path_for(book).name == "vorgaenge.bak.xlsx"


def test_a_file_that_does_not_exist_yields_no_copy(ws):
    assert bk.make_backup(ws / "gibtsnicht.xlsx") is None


def test_backups_are_listed_oldest_first(book):
    for _ in range(3):
        bk.make_backup(book)
        book.write_bytes(book.read_bytes() + b" ")
    assert [p.name for p in bk.list_backups(book)] == [
        "vorgaenge.bak.xlsx", "vorgaenge.bak2.xlsx", "vorgaenge.bak3.xlsx"]


def test_another_documents_copies_are_not_listed(book):
    other = book.parent / "andere.xlsx"
    other.write_bytes(b"x")
    bk.make_backup(book)
    bk.make_backup(other)
    assert [p.name for p in bk.list_backups(book)] == ["vorgaenge.bak.xlsx"]


# ---------------------------------------------------------------------------
# The agent writes into that history
# ---------------------------------------------------------------------------

def test_an_edit_keeps_a_copy_and_names_it(ws, book):
    ex = _DocToolExecutor()
    perms = _perms(ws)
    ex.execute("read_document", {"path": str(book)}, perms)
    out = json.loads(ex._dispatch("edit_sheet", {
        "path": str(book), "key_column": "Beleg",
        "updates": [{"key": "R-001", "set": {"Status": "gebucht"}}],
    }, perms))
    assert out["status"] == "ok"
    assert out["backup"].endswith("Backups/vorgaenge.bak.xlsx")
    assert (ws / "Backups" / "vorgaenge.bak.xlsx").exists()


def test_repeated_edits_build_a_history(ws, book):
    ex = _DocToolExecutor()
    perms = _perms(ws)
    for value in ("gebucht", "storniert", "offen"):
        ex.execute("read_document", {"path": str(book)}, perms)
        ex._dispatch("edit_sheet", {
            "path": str(book), "key_column": "Beleg",
            "updates": [{"key": "R-001", "set": {"Status": value}}],
        }, perms)
    assert [p.name for p in bk.list_backups(book)] == [
        "vorgaenge.bak.xlsx", "vorgaenge.bak2.xlsx", "vorgaenge.bak3.xlsx"]


def test_the_copy_holds_the_state_from_before_the_edit(ws, book):
    ex = _DocToolExecutor()
    perms = _perms(ws)
    ex.execute("read_document", {"path": str(book)}, perms)
    ex._dispatch("edit_sheet", {
        "path": str(book), "key_column": "Beleg",
        "updates": [{"key": "R-001", "set": {"Status": "gebucht"}}],
    }, perms)
    wb = openpyxl.load_workbook(ws / "Backups" / "vorgaenge.bak.xlsx")
    assert wb.active["B2"].value == "offen"
    wb.close()
    wb = openpyxl.load_workbook(book)
    assert wb.active["B2"].value == "gebucht"
    wb.close()


def test_creating_a_new_file_needs_no_copy(ws):
    """There is no previous state to keep."""
    ex = _DocToolExecutor()
    out = json.loads(ex._dispatch("edit_sheet", {
        "path": str(ws / "neu.xlsx"), "create": True,
        "append_rows": [["A"], ["1"]],
    }, _perms(ws)))
    assert out["status"] == "ok"
    assert "backup" not in out
    assert not (ws / "Backups").exists()
