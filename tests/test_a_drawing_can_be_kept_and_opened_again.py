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


def test_an_opened_drawing_is_fitted_to_the_frame_it_is_read_in():
    """``setMolecule`` fits it itself -- against the size the frame has at the
    moment it runs, which for a drawing opened into a tab that is not on
    screen yet is not the size it will be looked at.  That is how a structure
    ends up at 10% zoom somewhere off the side."""
    script = ketcher.load_js(".frame", '{"root":{}}')

    assert "zoomAccordingContent" in script
    assert "centerStruct" in script
    # and again once the frame has the size it is going to keep
    assert script.count("fit(api)") >= 3


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


def test_the_editors_own_save_button_writes_into_the_calc_folder(editor):
    """Ketcher already has the two buttons a chemist looks for, with a name
    field and a format list behind Save.  A second pair beside the frame was a
    second answer to a question the editor had already answered, so these are
    the ones that were made to work.

    The name comes back the way Ketcher wrote it -- "my aspirin.mol", stem and
    suffix together -- because it is the name typed into Ketcher's own dialog.
    """
    written = "\n  Ketcher  9 2\n\n  3  2  0  0\n"
    editor["submit_draw_file_sync"].value = (
        '1\nsave\nmy aspirin.mol\n' + written)

    kept = editor["tmp_path"] / "calc" / "Ketcher" / "my_aspirin.mol"
    assert kept.exists(), list((editor["tmp_path"] / "calc" / "Ketcher").glob("*"))
    assert kept.read_text() == written
    assert "Ketcher" in editor["mol_status"].value


def test_the_download_ketcher_would_have_made_never_happens(editor):
    """Ketcher writes a file the way every browser application does: a Blob, an
    object URL, a detached ``<a download>`` and a click dispatched at it.  None
    of that reaches the page's own listeners, so the two ends are taken
    instead -- and the click is swallowed rather than passed on, or the
    drawing would land in the browser's downloads folder as well."""
    wiring = ketcher.wire_js(".frame", ".sync")

    assert "createObjectURL" in wiring and "dispatchEvent" in wiring
    assert "hasAttribute('download')" in wiring
    assert "return true;" in wiring, "the click is swallowed, not forwarded"
    # Save as SVG or PNG is a picture for somewhere else, and still downloads.
    assert "KEEP.indexOf(ext)>=0" in wiring


def test_the_editors_own_open_button_is_answered_with_what_is_kept(editor):
    """A file picker shows the machine the *browser* is on, which down an SSH
    tunnel is not the machine the calculations are on.  So the button is
    answered with the drawings that are actually kept."""
    wiring = ketcher.wire_js(".frame", ".sync")
    assert "open-file-button" in wiring
    assert "ev.preventDefault(); ev.stopPropagation();" in wiring
    # and a genuinely local file is still one press away
    assert "__delfinKetcherPass" in wiring
    assert "From this computer" in wiring
    # Ketcher carries two Open buttons, one per mode, and the first in the
    # document is the hidden macromolecule editor's. Clicking that one does
    # nothing at all, which is what "From this computer" did.
    assert "offsetParent!==null" in wiring


def test_a_kept_drawing_asked_for_by_name_is_opened(editor):
    (editor["tmp_path"] / "calc" / "Ketcher").mkdir()
    (editor["tmp_path"] / "calc" / "Ketcher" / "benzene.mol").write_text("drawn")
    editor["submit_draw_open_btn"].value = False
    editor["submit_draw_open_btn"].value = True     # rescans the folder
    editor["sent_js"].clear()

    editor["submit_draw_file_sync"].value = "1\nopen\nbenzene.mol"

    script = "\n".join(editor["sent_js"])
    assert "setMolecule" in script and "drawn" in script


def test_a_name_that_is_gone_is_said_rather_than_pushed(editor):
    editor["submit_draw_file_sync"].value = "1\nopen\nnever_existed.ket"

    assert "not there any more" in editor["mol_status"].value


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


# ---------------------------------------------------------------------------
# the drawing goes with the job
# ---------------------------------------------------------------------------
def test_a_submitted_job_keeps_the_drawing_it_was_drawn_from(editor, tmp_path):
    """A CONTROL.txt says what was asked and an input.txt says of what.  A
    structure that was drawn rather than typed has a third thing to say, and
    until now it was said only inside a browser frame that the next page load
    empties."""
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(mol)
    editor["submit_draw_sync"].value = (
        "1\n" + Chem.MolToMolBlock(mol) + ketcher.KET_MARK + '{"root":{}}')
    assert editor["coords_widget"].value == "CCO", "the drawing was read"

    where = tmp_path / "calc" / "drawn_job"
    where.mkdir(parents=True)
    editor["keep_the_drawing"](where, "CCO")

    assert (where / "drawing.ket").read_text() == '{"root":{}}'


def test_a_job_that_was_typed_does_not_carry_someone_elses_drawing(editor,
                                                                   tmp_path):
    """The editor remembers the SMILES its drawing produced, so a job set up
    an hour later from a typed SMILES does not end up with an unrelated
    picture in its folder."""
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(mol)
    editor["submit_draw_sync"].value = (
        "1\n" + Chem.MolToMolBlock(mol) + ketcher.KET_MARK + '{"root":{}}')

    where = tmp_path / "calc" / "typed_job"
    where.mkdir(parents=True)
    editor["keep_the_drawing"](where, "c1ccccc1")

    assert not (where / "drawing.ket").exists()
    # and nothing drawn at all writes nothing
    editor["keep_the_drawing"](where, "")
    assert not (where / "drawing.ket").exists()


def test_the_drawing_is_written_where_control_and_input_are():
    """Beside them, in the same folder, by the same press."""
    from editor_source import TAB_SOURCE

    handler = TAB_SOURCE.split("def handle_submit(button)")[1]
    handler = handler.split("\n    def ")[0]
    order = [handler.index("CONTROL.txt"), handler.index("input.txt"),
             handler.index("_keep_the_drawing(")]
    assert order == sorted(order), "written after the two it stands beside"


# ---------------------------------------------------------------------------
# saved where it came from, in the format that keeps everything
# ---------------------------------------------------------------------------
def test_a_drawing_is_saved_into_whatever_folder_it_belongs_to(tmp_path):
    """The store is one particular folder, not a different kind of thing: a
    drawing opened out of a job's directory is saved back into that
    directory, the way a document is."""
    made = ketcher.save_into(tmp_path / "some_job", "aspirin", "x")

    assert made["ok"] is True
    assert made["path"] == tmp_path / "some_job" / "aspirin.ket"
    assert "some_job" in made["status"]
    assert [item.name for item in ketcher.list_in(tmp_path / "some_job")] == [
        "aspirin.ket"]
    # and the store is that, pointed at the Ketcher folder
    assert ketcher.save_drawing(tmp_path, "kept", "x")["path"] == (
        tmp_path / "Ketcher" / "kept.ket")


def test_the_save_dialog_is_opened_on_ket(tmp_path):
    """Of the formats Ketcher writes it is the only one that keeps an arrow, a
    text label and the layout, so a drawing saved and opened again is the
    drawing that was saved.

    The list is opened and the entry is pressed, which is what a person would
    do.  It is a Material select, and its own hidden input can be set without
    React noticing.
    """
    wiring = ketcher.wire_js(".frame", ".sync")

    assert 'save-dialog' in wiring
    assert "input.MuiSelect-nativeInput" in wiring
    assert "now.value==='ket'" in wiring, "already on Ket is left alone"
    assert "/^ket/i" in wiring
    assert "all[i].click()" in wiring, "pressed, not assigned"


# ---------------------------------------------------------------------------
# nothing is thrown away without asking
# ---------------------------------------------------------------------------
def test_replacing_unsaved_work_asks_first():
    """The clean mark is set when a drawing is opened and when one is saved,
    so a difference from it is work that opening another one would lose."""
    opening = ketcher.load_js(".frame", "x")

    assert "confirm(" in opening
    assert "__delfinKetcherClean" in opening
    assert "if(!go) return;" in opening, "saying no keeps the drawing"

    saving = ketcher.wire_js(".frame", ".sync")
    assert "__delfinKetcherClean" in saving, "saving makes it clean again"


def test_an_editor_nobody_has_opened_anything_into_has_nothing_to_lose():
    """Unknown means unknown, not dirty.  A first drawing must not be met with
    a question about losing something that was never there."""
    opening = ketcher.load_js(".frame", "x")

    assert "if(clean===undefined||clean===null){ go(true); return; }" in opening
