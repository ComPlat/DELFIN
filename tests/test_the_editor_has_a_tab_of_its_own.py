"""Ketcher at the size it was built to be used at.

The structure editor's Ketcher is a 560 px frame under a DRAW toggle, beside
the box it feeds.  That is the right size for "I would rather draw this than
type it" and the wrong size for drawing a reaction, keeping it, and opening it
again -- which is why the editor also has a tab.

The frame and the two channels are the same design the structure editor
established: a question out through ``run_js``, an answer back through a hidden
textarea, and a serial in front of it because a widget only reports a value
that changed.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.dashboard import ketcher, ketcher_panel
from delfin.dashboard.context import DashboardContext


def molblock(smiles):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    return Chem.MolToMolBlock(mol)


def rxnblock(reactants, products):
    text = "$RXN\n\n\n\n%3d%3d\n" % (len(reactants), len(products))
    for body in list(reactants) + list(products):
        text += "$MOL\n" + body
    return text


@pytest.fixture
def tab(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_ketcher

    for name in ("calc", "home"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    folder = tmp_path / ".delfin" / "ketcher"
    folder.mkdir(parents=True)
    (folder / "index.html").write_text("<html></html>")
    (folder / ".delfin-ketcher-version").write_text("3.17.0")

    ctx = DashboardContext(calc_dir=tmp_path / "calc")
    sent = []
    ctx.run_js = sent.append
    _widget, refs = tab_ketcher.create_tab(ctx)
    return {"ctx": ctx, "panel": refs["ketcher_panel"], "sent": sent,
            "tmp_path": tmp_path}


# ---------------------------------------------------------------------------
# what it says to the page
# ---------------------------------------------------------------------------
def test_nothing_is_sent_while_the_page_is_still_being_built(tab):
    """``create_dashboard`` collects every tab's startup script and sends them
    as one, after all of the tabs are built.  A run_js from inside a factory
    is a script against a page that has not been drawn yet -- and it would
    clear whatever the tab before it had sent."""
    assert tab["sent"] == []
    assert len(tab["ctx"].init_js_parts) == 1
    assert "pointerenter" in tab["ctx"].init_js_parts[0]


def test_every_press_sends_exactly_one_script(tab):
    """run_js clears its output before displaying the next, so two calls in a
    row can mean the first is thrown away before the browser has run it."""
    panel = tab["panel"]
    panel.name_box.value = "a"
    for press in (panel.smiles_btn, panel.save_btn, panel.smiles_copy_btn,
                  panel.rxn_copy_btn):
        tab["sent"].clear()
        press.click()
        assert len(tab["sent"]) == 1, press.description


def test_the_question_carries_the_focus_binding_with_it(tab):
    tab["sent"].clear()
    tab["panel"].smiles_btn.click()

    script = tab["sent"][0]
    assert "contentWindow.focus()" in script
    assert "containsReaction" in script, "an arrow means an RXN file"
    assert "getMolfile" in script, "and no arrow means a molfile"


def test_the_answer_says_which_question_it_answers(tab):
    """One hidden channel carries both -- what to put in the SMILES box and
    what to write to disk -- so the answer has to say which it is."""
    read = ketcher_panel.read_js("scope-1", "save", "ket")

    assert "'\\n'+kind+'\\n'+text" in read
    assert "Date.now()" in read, (
        "drawing the same thing twice would otherwise look like an answer "
        "that never came"
    )
    assert ".delfin-ketcher-sync.'+scope" in read, "never 'the' channel"


def test_the_structure_editor_asks_the_same_question_of_its_own_elements(tab):
    """It has a Ketcher of its own, under names of its own that predate this
    module.  The question is the same one, so it is asked from one place."""
    theirs = ketcher_panel.read_js(
        "scope-1", "save", "ket",
        frame_class="submit-ketcher-frame",
        sync_class="submit-ketcher-file-sync")

    assert ".submit-ketcher-file-sync.'+scope" in theirs
    assert ".submit-ketcher-frame.'+scope" in theirs
    assert ".delfin-ketcher-" not in theirs


# ---------------------------------------------------------------------------
# what comes back
# ---------------------------------------------------------------------------
def test_a_structure_and_a_reaction_do_not_share_a_box(tab):
    panel = tab["panel"]

    panel.sync.value = "1\nsmiles\n" + molblock("CCO")
    assert panel.smiles_out.value == "CCO"
    assert panel.rxn_out.value == ""

    panel.sync.value = "2\nsmiles\n" + rxnblock(
        [molblock("CCO")], [molblock("CC=O")])
    assert panel.rxn_out.value == "CCO>>CC=O"
    assert panel.smiles_out.value == "CCO", "the structure box is left alone"


def test_the_drawing_is_written_where_the_name_said_at_the_time(tab):
    """The name box belongs to the user and they may have moved on between
    pressing SAVE and the drawing coming back."""
    panel = tab["panel"]
    panel.name_box.value = "aspirin"
    panel.save_btn.click()
    panel.name_box.value = "something else entirely"

    panel.sync.value = '1\nsave\n{"root":{}}'

    kept = tab["tmp_path"] / "calc" / "Ketcher" / "aspirin.ket"
    assert kept.read_text() == '{"root":{}}'
    assert panel.files_dd.value == "aspirin.ket"


def test_the_structure_can_be_carried_to_the_tab_that_submits_it(tab):
    """Every field a calculation needs is in the Submit tab.  A second, staler
    copy of them here would be worse than a tab switch -- the same handover
    the reaction graph makes."""
    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    panel = tab["panel"]
    box = widgets.Textarea(value="")
    tab["ctx"].submit_refs = {"coords_widget": box}
    tab["ctx"].tabs_widget = type("Tabs", (), {"selected_index": None})()
    tab["ctx"].tab_indices = {"Submit Job": 0}
    panel.sync.value = "1\nsmiles\n" + molblock("CCO")

    panel.to_submit_btn.click()

    assert box.value == "CCO"
    assert tab["ctx"].tabs_widget.selected_index == 0


def test_nothing_drawn_is_not_carried_anywhere(tab):
    tab["ctx"].submit_refs = {}
    tab["panel"].to_submit_btn.click()

    assert "nothing to send" in tab["panel"].status.value


def test_an_editor_that_is_not_there_says_so_rather_than_nothing(tab):
    tab["panel"].sync.value = "1\nsmiles\n!no-editor"

    assert "not open yet" in tab["panel"].status.value


def test_a_half_written_answer_is_not_read_as_one(tab):
    """The channel is set once with all three parts.  A value with fewer is a
    value that is not an answer."""
    panel = tab["panel"]
    panel.smiles_out.value = "kept"

    panel.sync.value = "1"
    panel.sync.value = "1\nsmiles"

    assert panel.smiles_out.value == "kept"


# ---------------------------------------------------------------------------
# when it is not installed
# ---------------------------------------------------------------------------
def test_it_is_offered_rather_than_fetched(tmp_path, monkeypatch):
    """Thirty-odd megabytes over the network, and on a machine without one it
    is a wait that ends in nothing."""
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_ketcher

    (tmp_path / "calc").mkdir()
    (tmp_path / "home").mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    ctx = DashboardContext(calc_dir=tmp_path / "calc")
    ctx.run_js = lambda _script: None

    _widget, refs = tab_ketcher.create_tab(ctx)
    panel = refs["ketcher_panel"]

    assert panel.frame.value == "", "nothing to show"
    assert panel.install_btn.description == "FETCH KETCHER"
    assert panel.smiles_btn.disabled is True
    assert "32 MB" in panel.status.value

    building = pathlib.Path(ketcher.__file__).read_text(encoding="utf-8")
    assert "def install(" in building
    made = ketcher_panel.build.__doc__
    assert made, "the factory says what it makes"


def test_an_installed_editor_can_still_be_brought_up_to_date(tab):
    """An editor that is there is an editor that can be updated, and the tab
    is the place to do it from."""
    assert tab["panel"].install_btn.description == "UPDATE KETCHER"
    assert tab["panel"].install_btn.disabled is False


def test_a_drawing_cannot_be_opened_into_an_editor_that_is_not_there(
        tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_ketcher

    (tmp_path / "calc").mkdir()
    (tmp_path / "home").mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    ctx = DashboardContext(calc_dir=tmp_path / "calc")
    sent = []
    ctx.run_js = sent.append
    _widget, refs = tab_ketcher.create_tab(ctx)

    assert refs["open_drawing"]('{"root":{}}', "aspirin.ket") is False
    assert sent == [], "nothing is pushed into a frame that is not there"
    assert "FETCH KETCHER" in refs["ketcher_panel"].status.value
