"""An arrow on the canvas makes it a reaction, and a reaction is not a molfile.

Ketcher refuses to write one: asked for a molfile with an arrow on the canvas
it throws *"The structure cannot be saved as *.MOL due to reaction"*.  TO SMILES
asked for a molfile and nothing else, so a drawn reaction came back as that
error and the editor looked broken to anyone who drew one.

Asked for what it is, it answers -- and what comes back is the ordinary
reaction SMILES, ``reactants>agents>products``, which is what RDKit writes and
what the reaction SMARTS elsewhere in DELFIN already read.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.dashboard import ketcher
from delfin.dashboard.context import DashboardContext
from editor_source import SUBMIT_SOURCE


def molblock(smiles):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    AllChem.Compute2DCoords(mol)
    return Chem.MolToMolBlock(mol, kekulize=False)


def rxnblock(reactants, products):
    """An RXN file of the shape Ketcher hands back."""
    text = "$RXN\n\n\n\n%3d%3d\n" % (len(reactants), len(products))
    for body in list(reactants) + list(products):
        text += "$MOL\n" + body
    return text


def placed(smiles, x, y=0.0):
    """A molblock for *smiles*, drawn where it was put on the canvas.

    Coordinates are the whole of how a scheme is read: an RXN file keeps them
    and keeps nothing else about the layout.
    """
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    AllChem.Compute2DCoords(mol)
    frame = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        spot = frame.GetAtomPosition(i)
        frame.SetAtomPosition(i, (spot.x + x, spot.y + y, 0.0))
    return Chem.MolToMolBlock(mol, kekulize=False)


def canvas(arrows):
    """A KET holding just the arrows, which is all their positions need."""
    import json

    return json.dumps({'root': {'nodes': [
        {'type': 'arrow', 'data': {'mode': 'open-angle', 'pos': [
            {'x': x0, 'y': y, 'z': 0}, {'x': x1, 'y': y, 'z': 0}]}}
        for x0, x1, y in arrows]}})


# ---------------------------------------------------------------------------
# reading it
# ---------------------------------------------------------------------------
def test_two_things_drawn_into_one_come_back_joined_by_an_arrow():
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([molblock("CCO"), molblock("CC(=O)O")],
                 [molblock("CCOC(C)=O")]))

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "CC(=O)O.CCO>>CCOC(C)=O"
    assert ">>" in outcome["smiles"], "an arrow, not a dot"


def test_a_reaction_between_things_no_valence_table_likes():
    """The same tolerance the molecule path has.  A reaction whose reactant is
    a metal complex is still a reaction, and refusing it would refuse the
    reason for having a drawing editor."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([molblock("C[Pt](Cl)(Cl)N")], [molblock("CCO")]))

    assert outcome["ok"] is True, outcome["status"]
    assert "Pt" in outcome["smiles"] and ">>" in outcome["smiles"]


def test_an_arrow_with_nothing_after_it_is_said_plainly():
    """Half a reaction is not an error to raise, it is a thing to say."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([molblock("CCO")], []))

    assert outcome["ok"] is False
    assert "nothing after it" in outcome["status"]


def test_nothing_drawn_and_nonsense_are_both_answered_not_raised():
    assert ketcher.reaction_smiles_from_rxnfile("")["ok"] is False
    assert "Nothing was drawn" in ketcher.reaction_smiles_from_rxnfile("")["status"]
    assert ketcher.reaction_smiles_from_rxnfile("$RXN\nnot one")["ok"] is False


# ---------------------------------------------------------------------------
# telling one from the other
# ---------------------------------------------------------------------------
def test_the_first_line_says_which_it_is():
    """One call site for both.  An RXN file begins with ``$RXN`` and a molfile
    does not, so nothing has to be remembered about what was asked for."""
    drawn = ketcher.smiles_from_drawing(molblock("CCO"))
    assert drawn["reaction"] is False and drawn["smiles"] == "CCO"

    reaction = ketcher.smiles_from_drawing(
        rxnblock([molblock("CCO")], [molblock("CC=O")]))
    assert reaction["reaction"] is True
    assert reaction["smiles"] == "CCO>>CC=O"


def test_the_editor_is_asked_for_a_reaction_only_when_there_is_one():
    """``getRxn`` throws when there is no arrow, exactly as ``getMolfile``
    throws when there is one.  Which to ask for is decided in the page, where
    the canvas is, rather than guessed here."""
    handler = SUBMIT_SOURCE.split("def on_submit_draw_get")[1].split("\n    def ")[0]

    assert "containsReaction" in handler
    assert "api.getRxn()" in handler
    assert "api.getMolfile()" in handler
    # Still not the editor's own SMILES: everything downstream reads with
    # RDKit, and a SMILES RDKit wrote is one RDKit will certainly read back.
    assert "getSmiles" not in handler


# ---------------------------------------------------------------------------
# where it lands
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
    return refs


def test_a_reaction_lands_in_the_box_a_structure_lands_in(editor):
    """The same box.  It is where the rest of the tab looks, and a reaction
    SMILES asked for is a reaction SMILES wanted there -- a second field
    beside it would only hold a copy of what is already on screen."""
    editor["coords_widget"].value = ""

    editor["submit_draw_sync"].value = "1\n" + rxnblock(
        [molblock("CCO"), molblock("CC(=O)O")], [molblock("CCOC(C)=O")])

    assert editor["coords_widget"].value == "CC(=O)O.CCO>>CCOC(C)=O"


def test_it_says_that_convert_has_nothing_to_do_with_a_reaction(editor):
    """Convert builds one structure.  Said when the reaction arrives, rather
    than found out by pressing Convert and reading whatever it makes of
    "A.B>>C"."""
    editor["submit_draw_sync"].value = "1\n" + rxnblock(
        [molblock("CCO")], [molblock("CC=O")])

    said = editor["mol_status"].value
    assert "reaction" in said and "Convert" in said


def test_a_plain_drawing_still_lands_where_it_always_did(editor):
    """The reaction path must not have moved the ordinary one."""
    editor["submit_draw_sync"].value = "1\n" + molblock("CCO")

    assert editor["coords_widget"].value == "CCO"


def test_the_editor_carries_no_second_field_for_a_reaction(editor):
    for gone in ("submit_draw_rxn_out", "submit_draw_rxn_row",
                 "submit_draw_rxn_copy_btn"):
        assert gone not in editor, gone


# ---------------------------------------------------------------------------
# four places around an arrow, and only two of them survive an RXN file
# ---------------------------------------------------------------------------
def test_what_is_over_the_arrow_is_what_is_added():
    """Ketcher writes it as an agent already, and this keeps that."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCCCC1", 14), placed("CC", 8, 4)]),
        canvas([(6.0, 10.0, 0.0)]))

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "c1ccccc1>CC>C1CCCCC1"


def test_what_is_under_the_arrow_is_what_comes_off():
    """Indigo cannot tell the two apart -- measured, a cyclobutane over the
    arrow and a cyclopropane under it both came back as agents, as
    ``C1C=CC=CC=1>C1CCC1.C1CC1>C1CCCCC1``.  Under the arrow is what the step
    gives off, so it belongs with the products."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCCCC1", 14), placed("O", 8, -4)]),
        canvas([(6.0, 10.0, 0.0)]))

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "c1ccccc1>>C1CCCCC1.O"


def test_a_reactant_that_merely_reaches_across_the_line_is_not_a_reagent():
    """Reactants and products straddle the arrow's line -- measured at
    y -8.18..-6.17 against an arrow at -7.18 -- so being near it is not
    enough.  Over and under mean clear of it."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 8)], [placed("C1CCCCC1", 20)]),
        canvas([(14.0, 18.0, 0.0)]))

    assert outcome["smiles"] == "c1ccccc1>>C1CCCCC1"


def test_several_arrows_are_read_as_the_steps_they_were_drawn_as():
    """An RXN file holds one arrow, so Indigo flattens three steps into "the
    first thing, into everything else" and does not even keep the drawn order.
    The arrows survive in the KET, and every component keeps its coordinates
    in the same frame, so the steps come off the geometry."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCC1", 24), placed("C1CCCCC1", 12)]),
        canvas([(5.0, 8.0, 0.0), (17.0, 20.0, 0.0)]))

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "c1ccccc1>>C1CCCCC1>>C1CCC1"
    assert outcome["steps"] == 2


def test_a_reagent_belongs_to_the_step_it_is_drawn_over():
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCC1", 24), placed("C1CCCCC1", 12),
                  placed("CO", 18.5, 5)]),
        canvas([(5.0, 8.0, 0.0), (17.0, 20.0, 0.0)]))

    assert outcome["smiles"] == "c1ccccc1>>C1CCCCC1>CO>C1CCC1"


def test_without_a_canvas_the_rxn_files_own_split_is_used():
    """Which is the right answer for one arrow and no geometry to read."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("CCO", 0)], [placed("CC=O", 10)]))

    assert outcome["smiles"] == "CCO>>CC=O"
