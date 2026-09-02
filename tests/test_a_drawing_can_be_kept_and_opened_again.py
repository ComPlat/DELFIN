"""A drawing used to live only in the frame it was drawn in.

The flow was one-way -- editor to Python, once, when TO SMILES was pressed --
and reloading the page lost it.  Nothing ever called ``setMolecule``, so there
was no way back in either: a structure that took ten minutes to draw could not
be put away and picked up again.

It is kept in ``<calc>/Ketcher`` now, beside the work it was drawn for rather
than in a store of its own that nothing else can see, and it is ``.ket`` by
default because that is the only one of these formats that keeps an arrow, a
text label and the layout.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.dashboard import ketcher
from delfin.dashboard.context import DashboardContext


# ---------------------------------------------------------------------------
# the folder
# ---------------------------------------------------------------------------
def test_it_is_kept_beside_the_calculations(tmp_path):
    assert ketcher.drawings_directory(tmp_path) == tmp_path / "Ketcher"
    assert ketcher.DRAWINGS_FOLDER == "Ketcher"
    assert ketcher.list_drawings(tmp_path) == [], "no folder is not an error"


def test_a_drawing_is_written_read_back_and_thrown_away(tmp_path):
    written = ketcher.save_drawing(tmp_path, "aspirin", '{"root":{}}')

    assert written["ok"] is True, written["status"]
    assert written["path"] == tmp_path / "Ketcher" / "aspirin.ket"
    assert [item.name for item in ketcher.list_drawings(tmp_path)] == ["aspirin.ket"]

    read = ketcher.read_drawing(written["path"])
    assert read["ok"] is True and read["text"] == '{"root":{}}'
    assert read["name"] == "aspirin.ket"

    assert ketcher.delete_drawing(written["path"])["ok"] is True
    assert ketcher.list_drawings(tmp_path) == []


def test_ket_is_the_default_because_it_is_the_one_that_keeps_everything(tmp_path):
    """A .mol cannot hold an arrow and a .rxn cannot hold a lone structure.
    Both are for handing something to a program that is not Ketcher."""
    assert ketcher.DRAWING_SUFFIXES[0] == ".ket"
    assert ketcher.save_drawing(tmp_path, "a", "x")["path"].suffix == ".ket"

    for suffix in (".mol", ".rxn", ".smi", ".cdxml"):
        made = ketcher.save_drawing(tmp_path, "a", "x", suffix)
        assert made["ok"] is True and made["path"].suffix == suffix

    refused = ketcher.save_drawing(tmp_path, "a", "x", ".exe")
    assert refused["ok"] is False and "not a format" in refused["status"]


def test_a_name_cannot_leave_the_folder(tmp_path):
    """The name is the user's and goes into a path.  ``safe_name`` is what
    every other named folder in DELFIN is made with, so it is what this is
    made with too -- and the resolved path is checked afterwards as well,
    because a rule that is only applied is a rule that can be got round."""
    made = ketcher.save_drawing(tmp_path, "../../etc/passwd", "x")

    assert made["ok"] is True
    assert made["path"].parent == tmp_path / "Ketcher"
    assert not (tmp_path.parent / "passwd.ket").exists()
    assert made["path"].name == "etc_passwd.ket"


def test_a_name_that_already_carries_the_suffix_does_not_get_it_twice(tmp_path):
    """Typing "acetone.ket" into the name box means the drawing is called
    acetone."""
    assert ketcher.save_drawing(
        tmp_path, "acetone.ket", "x", ".ket")["path"].name == "acetone.ket"


def test_an_empty_name_and_an_empty_drawing_are_both_said_plainly(tmp_path):
    assert "name" in ketcher.save_drawing(tmp_path, "  ", "x")["status"]
    assert "nothing drawn" in ketcher.save_drawing(tmp_path, "a", "  ")["status"]
    gone = ketcher.read_drawing(tmp_path / "Ketcher" / "never.ket")
    assert gone["ok"] is False and "not there" in gone["status"]


def test_only_what_the_editor_can_open_is_offered(tmp_path):
    folder = tmp_path / "Ketcher"
    folder.mkdir()
    for name in ("a.ket", "b.mol", "c.rxn", "notes.txt", "run.out"):
        (folder / name).write_text("x")

    assert [item.name for item in ketcher.list_drawings(tmp_path)] == [
        "a.ket", "b.mol", "c.rxn"]
    assert ketcher.is_drawing("x.KET") is True, "however it happens to be cased"
    assert ketcher.is_drawing("x.txt") is False


# ---------------------------------------------------------------------------
# putting it back
# ---------------------------------------------------------------------------
def test_a_kept_drawing_is_put_back_with_one_call_whatever_it_is():
    """Indigo reads what it is given, so a .ket, a .mol, a .rxn and a SMILES
    all go back in through the same line."""
    script = ketcher.load_js(".frame", '{"root":{}}')

    assert "setMolecule" in script
    assert '{\\"root\\":{}}' in script, "the drawing travels inside the script"
    # Ketcher is 30 MB and the global this reaches for does not exist for the
    # first second or two after the frame is put on the page.
    assert "setTimeout" in script and "tries" in script


# ---------------------------------------------------------------------------
# in the tab
# ---------------------------------------------------------------------------
@pytest.fixture
def editor(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit

    for name in ("calc", "archive", "office", "home"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    folder = tmp_path / ".delfin" / "ketcher"
    folder.mkdir(parents=True)
    (folder / "index.html").write_text("<html></html>")
    (folder / ".delfin-ketcher-version").write_text("3.17.0")
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    sent = []
    ctx.run_js = sent.append
    _widget, refs = tab_submit.create_tab(ctx)
    refs["submit_draw_open_btn"].value = True
    refs["sent_js"] = sent
    refs["tmp_path"] = tmp_path
    return refs


def test_saving_asks_the_editor_then_writes_what_comes_back(editor):
    """Two steps, because the drawing is in the browser and the disk is not."""
    editor["submit_draw_name"].value = "ester"
    editor["sent_js"].clear()

    editor["submit_draw_save_btn"].click()

    asked = "\n".join(editor["sent_js"])
    assert '"ket"' in asked, "the .ket the format box is showing"
    assert "submit-ketcher-file-sync" in asked, (
        "a channel of its own: the other one carries a structure on its way "
        "to the input box, and its shape is what pins that path"
    )

    editor["submit_draw_file_sync"].value = '1\nsave\n{"root":{}}'

    kept = editor["tmp_path"] / "calc" / "Ketcher" / "ester.ket"
    assert kept.read_text() == '{"root":{}}'
    assert editor["submit_draw_files_dd"].options == ("ester.ket",)


def test_a_drawing_with_no_name_is_not_asked_for_at_all(editor):
    editor["submit_draw_name"].value = "   "
    editor["sent_js"].clear()

    editor["submit_draw_save_btn"].click()

    assert editor["sent_js"] == [], "nothing is asked of the page"
    assert "name" in editor["mol_status"].value


def test_opening_one_shows_the_editor_and_hands_the_drawing_over(editor):
    (editor["tmp_path"] / "calc" / "Ketcher").mkdir()
    (editor["tmp_path"] / "calc" / "Ketcher" / "benzene.mol").write_text("drawn")
    editor["submit_draw_open_btn"].value = False
    editor["submit_draw_open_btn"].value = True   # refreshes the list
    editor["sent_js"].clear()

    editor["submit_draw_open_file_btn"].click()

    assert len(editor["sent_js"]) == 1, (
        "run_js clears its output before displaying the next script, so the "
        "focus binding travels with the drawing rather than after it"
    )
    script = editor["sent_js"][0]
    assert "setMolecule" in script and "drawn" in script
    assert "pointerenter" in script
    # The name box follows the file, so pressing SAVE writes it back where it
    # came from rather than under whatever was last typed.
    assert editor["submit_draw_name"].value == "benzene"
    assert editor["submit_draw_format_dd"].value == ".mol"


def test_a_drawing_handed_in_from_elsewhere_opens_the_editor_first(editor):
    """The Calculations browser comes in this way.  Handing a drawing to a
    folded-away editor would be putting it somewhere nobody can see it."""
    editor["submit_draw_open_btn"].value = False

    assert editor["open_drawing"]('{"root":{}}', "aspirin.ket") is True

    assert editor["submit_draw_open_btn"].value is True
    assert editor["submit_draw_frame"].layout.display == ""


def test_an_empty_file_is_said_rather_than_pushed(editor):
    assert editor["open_drawing"]("   ", "empty.ket") is False
    assert "empty" in editor["mol_status"].value
