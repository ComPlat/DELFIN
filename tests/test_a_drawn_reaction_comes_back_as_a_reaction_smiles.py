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


def sdf(*blocks):
    """The records Ketcher hands back for a canvas with an arrow on it."""
    return ''.join(block + "$$$$\n" for block in blocks)


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


def test_the_editor_is_asked_for_the_records_only_when_there_is_a_reaction():
    """``getMolfile`` throws when there is an arrow, so which to ask for is
    decided in the page, where the canvas is, rather than guessed here.

    And for a reaction it is the records, not a reaction file: an RXN drops
    components a scheme still needs.  Both tabs ask the same way -- the Submit
    tab used to ask for ``getRxn`` and send no canvas at all, which is why what
    it put in the input box was the ordinary three-field form with over and
    under run together."""
    handler = SUBMIT_SOURCE.split("def on_submit_draw_get")[1].split("\n    def ")[0]

    assert "containsReaction" in handler
    assert "api.getSdf()" in handler
    assert "api.getMolfile()" in handler
    assert "api.getKet()" in handler, "the arrows are only in there"
    assert "api.getRxn()" not in handler
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
# four places around an arrow, and four fields to put them in
# ---------------------------------------------------------------------------
#
# ``reactants>>in>out>>products``.  Ordinary reaction SMILES has three fields,
# and a by-product put in with the products is one nothing can pick out again:
# a dot there means "and this one too", which is what a second product drawn
# beside the first one means.  The fourth field keeps the two apart, so a dot
# only ever separates several things standing in the same place.
#
# The doubled marks are the arrow's own two ends, so steps run on one after
# another and each can still be seen: ``A>>in>out>>B>>in>out>>C``.
def test_what_is_over_the_arrow_is_what_is_added():
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCCCC1", 14), placed("CC", 8, 4)]),
        canvas([(6.0, 10.0, 0.0)]))

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "c1ccccc1>>CC>>>C1CCCCC1"


def test_what_is_under_the_arrow_is_what_comes_off():
    """Indigo cannot tell the two apart -- measured, a cyclobutane over the
    arrow and a cyclopropane under it both came back as agents, as
    ``C1C=CC=CC=1>C1CCC1.C1CC1>C1CCCCC1``.  Here they stand in different
    fields, so a dot in the products means a second product and nothing
    else."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCCCC1", 14), placed("O", 8, -4)]),
        canvas([(6.0, 10.0, 0.0)]))

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "c1ccccc1>>>O>>C1CCCCC1"


def test_both_at_once_read_as_one_in_and_one_out():
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("CCC", 0)],
                 [placed("CCC", 14), placed("C", 8, 4), placed("C", 8, -4)]),
        canvas([(6.0, 10.0, 0.0)]))

    assert outcome["smiles"] == "CCC>>C>C>>CCC"


def test_a_side_with_several_things_on_it_is_written_as_the_one_side_it_is():
    """``(A.B)`` rather than ``A.B``.  It is reaction SMILES' own grouping
    mark -- RDKit reads it as one template holding two fragments and writes it
    back the same way -- and here it says the same thing: this is one side of
    an arrow, and the dots inside it are the dots that were drawn."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("CCO", 0), placed("CC(=O)O", 0, 5)],
                 [placed("CCOC(C)=O", 14), placed("O", 14, 5)]),
        canvas([(6.0, 10.0, 0.0)]))

    assert outcome["smiles"] == "(CC(=O)O.CCO)>>>>>(CCOC(C)=O.O)"


def test_one_thing_on_a_side_needs_no_group():
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("CCC", 0)], [placed("CCC", 14)]),
        canvas([(6.0, 10.0, 0.0)]))

    assert outcome["smiles"] == "CCC>>>>>CCC"


def test_what_is_over_or_under_an_arrow_is_never_grouped():
    """Reaction SMILES allows the mark on the two sides only, and RDKit
    refuses it in the middle: "Problems constructing agent from SMARTS"."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCCCC1", 14), placed("CC", 7, 5),
                  placed("CO", 9, 5)]),
        canvas([(6.0, 10.0, 0.0)]))

    assert outcome["smiles"] == "c1ccccc1>>CC.CO>>>C1CCCCC1"
    assert "(CC.CO)" not in outcome["smiles"]


def test_a_reactant_that_merely_reaches_across_the_line_is_not_a_reagent():
    """Reactants and products straddle the arrow's line -- measured at
    y -8.18..-6.17 against an arrow at -7.18 -- so being near it is not
    enough.  Over and under mean clear of it."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 8)], [placed("C1CCCCC1", 20)]),
        canvas([(14.0, 18.0, 0.0)]))

    assert outcome["smiles"] == "c1ccccc1>>>>>C1CCCCC1"


def test_several_arrows_run_on_one_after_another():
    """An RXN file holds one arrow, so Indigo flattens three steps into "the
    first thing, into everything else" and does not even keep the drawn order.
    The arrows survive in the KET, and every component keeps its coordinates
    in the same frame, so the steps come off the geometry.

    What stands between two arrows is written once: it is the products of the
    step before it and the reactants of the step after it, which is what it is
    on the canvas."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCC1", 24), placed("C1CCCCC1", 12)]),
        canvas([(5.0, 8.0, 0.0), (17.0, 20.0, 0.0)]))

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "c1ccccc1>>>>>C1CCCCC1>>>>>C1CCC1"
    assert outcome["steps"] == 2


def test_a_reagent_belongs_to_the_step_it_is_drawn_over():
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCC1", 24), placed("C1CCCCC1", 12),
                  placed("CO", 18.5, 5)]),
        canvas([(5.0, 8.0, 0.0), (17.0, 20.0, 0.0)]))

    assert outcome["smiles"] == "c1ccccc1>>>>>C1CCCCC1>>CO>>>C1CCC1"


def test_what_a_step_gives_off_stays_in_that_steps_own_field():
    """It is not something the next step is made from, and it is not a second
    product either -- it has a field of its own."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCC1", 24), placed("C1CCCCC1", 12),
                  placed("O", 6.5, -5)]),
        canvas([(5.0, 8.0, 0.0), (17.0, 20.0, 0.0)]))

    assert outcome["smiles"] == "c1ccccc1>>>O>>C1CCCCC1>>>>>C1CCC1"
    first, rest = outcome["smiles"].split(">>C1CCCCC1", 1)
    assert first.endswith(">O"), "given off by the first step"
    assert "O" not in rest, "and nowhere in the second"


def test_without_a_canvas_the_rxn_files_own_split_is_used():
    """Which is the right answer for one arrow and no geometry to read: the
    ordinary three-field form, because that is all an RXN file can say."""
    outcome = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("CCO", 0)], [placed("CC=O", 10)]))

    assert outcome["smiles"] == "CCO>>CC=O"


def test_structures_with_no_arrow_between_them_are_one_dotted_smiles():
    """No arrow, no reaction: what is on the canvas is a set of structures,
    and a set of structures is a dotted SMILES -- which is exactly what goes
    into the input box and is read there as the reactants would be.

    Ketcher hands back one molfile holding every fragment, so this needs
    nothing of the geometry at all.
    """
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CCO.c1ccccc1.[Na+].[Cl-]")
    AllChem.Compute2DCoords(mol)

    outcome = ketcher.smiles_from_drawing(Chem.MolToMolBlock(mol))

    assert outcome["reaction"] is False
    assert outcome["smiles"] == "CCO.[Cl-].[Na+].c1ccccc1"
    assert ">" not in outcome["smiles"], "nothing to separate, so no arrow"


# ---------------------------------------------------------------------------
# and back again, as reactions RDKit can work with
# ---------------------------------------------------------------------------
def test_the_four_fields_are_read_back_apart():
    """The form is what the drawing says; this is how anything that speaks
    ordinary reaction SMILES gets to read it."""
    pytest.importorskip("rdkit")

    read = ketcher.parse_reaction_smiles("CCC>>C>C>>CCC")

    assert read["ok"] is True, read["status"]
    step, = read["steps"]
    assert step["reactants"] == ["CCC"]
    assert step["in"] == ["C"]
    assert step["out"] == ["C"]
    assert step["products"] == ["CCC"]


def test_each_step_comes_back_as_a_reaction_rdkit_can_hold():
    """A by-product is a product: putting it there is what keeps the step
    balanced, and it is the only place an ordinary reaction SMILES has."""
    pytest.importorskip("rdkit")

    step, = ketcher.parse_reaction_smiles("CCC>>C>C>>CCC")["steps"]

    assert step["smiles"] == "CCC>C>(CCC.C)"
    built = step["reaction"]
    assert built.GetNumReactantTemplates() == 1
    assert built.GetNumAgentTemplates() == 1
    # One template, because the products of a step are one side of one arrow;
    # the two fragments in it are the product and what came off.
    assert built.GetNumProductTemplates() == 1
    product, = built.GetProducts()
    from rdkit import Chem
    assert sorted(Chem.MolToSmiles(f) for f in
                  Chem.GetMolFrags(product, asMols=True)) == ["C", "CCC"]


def test_a_chain_comes_back_as_one_step_after_another():
    pytest.importorskip("rdkit")

    read = ketcher.parse_reaction_smiles(
        "c1ccccc1>>>O>>C1CCCCC1>>CO>>>C1CCC1")

    assert read["ok"] is True, read["status"]
    first, second = read["steps"]
    assert (first["reactants"], first["out"], first["products"]) == (
        ["c1ccccc1"], ["O"], ["C1CCCCC1"])
    assert (second["reactants"], second["in"], second["products"]) == (
        ["C1CCCCC1"], ["CO"], ["C1CCC1"])
    assert "O" not in second["reactants"], (
        "what the first step gave off is not what the second is made from"
    )


def test_what_the_drawing_writes_is_what_the_parser_reads():
    """The round trip is the point: the geometry is read off the canvas once,
    written down, and understood again without the canvas."""
    pytest.importorskip("rdkit")

    written = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("c1ccccc1", 0)],
                 [placed("C1CCCCC1", 14), placed("CC", 8, 4),
                  placed("O", 8, -4)]),
        canvas([(6.0, 10.0, 0.0)]))
    assert written["smiles"] == "c1ccccc1>>CC>O>>C1CCCCC1"

    read = ketcher.parse_reaction_smiles(written["smiles"])

    step, = read["steps"]
    assert step["in"] == ["CC"], "drawn over the arrow"
    assert step["out"] == ["O"], "drawn under it"
    assert step["products"] == ["C1CCCCC1"]


def test_the_ordinary_three_field_form_is_still_read():
    """An RXN file with no canvas beside it produces one, and so does every
    other tool."""
    pytest.importorskip("rdkit")

    step, = ketcher.parse_reaction_smiles("CCO>[Pd]>CC=O")["steps"]

    assert step["in"] == ["[Pd]"] and step["out"] == []
    assert step["smiles"] == "CCO>[Pd]>CC=O"


def test_something_that_is_not_a_reaction_is_said_rather_than_guessed():
    pytest.importorskip("rdkit")

    assert ketcher.parse_reaction_smiles("")["ok"] is False
    structures = ketcher.parse_reaction_smiles("CCO.CCC")
    assert structures["ok"] is False
    assert "no arrow" in structures["status"]
    crooked = ketcher.parse_reaction_smiles("CCO>>>CC=O")
    assert crooked["ok"] is False and "fields" in crooked["status"]


def test_several_things_in_one_place_are_separated_by_a_dot():
    """A dot only ever means "and this one too", in whichever of the four
    fields it stands -- and a step with nothing over and nothing under its
    arrow is five marks with nothing between them."""
    pytest.importorskip("rdkit")

    read = ketcher.parse_reaction_smiles(
        "CCO.CC(=O)O>>[Pd]>O>>CCOC(C)=O.CCC>>>>>C1CCC1.CO")

    assert read["ok"] is True, read["status"]
    first, second = read["steps"]
    assert first["reactants"] == ["CCO", "CC(=O)O"]
    assert first["products"] == ["CCOC(C)=O", "CCC"]
    assert second["in"] == [] and second["out"] == [], "nothing in, nothing out"
    assert second["products"] == ["C1CCC1", "CO"]


def test_a_follow_up_step_is_made_from_what_the_one_before_it_produced():
    """Written once, because on the canvas it is one set of structures: the
    products of the step to its left and the reactants of the step to its
    right."""
    pytest.importorskip("rdkit")

    written = ketcher.reaction_smiles_from_rxnfile(
        rxnblock([placed("CCO", 0), placed("CC(=O)O", 0, 6)],
                 [placed("CCOC(C)=O", 13), placed("CCC", 13, 6),
                  placed("C1CCC1", 27), placed("CO", 27, 6),
                  placed("[Pd]", 7, 5), placed("O", 7, -5)]),
        canvas([(5.0, 9.0, 0.0), (18.5, 22.5, 0.0)]))

    assert written["smiles"] == (
        "(CC(=O)O.CCO)>>[Pd]>O>>(CCC.CCOC(C)=O)>>>>>(C1CCC1.CO)")

    first, second = ketcher.parse_reaction_smiles(written["smiles"])["steps"]
    assert first["products"] == second["reactants"], (
        "the same set of structures, standing between the two arrows"
    )
    assert "O" not in second["reactants"], (
        "and what the first step gave off is not among them"
    )


# ---------------------------------------------------------------------------
# every component, which an RXN file does not carry
# ---------------------------------------------------------------------------
def test_a_scheme_is_read_from_the_records_that_keep_every_component():
    """An RXN file loses some.  Measured against the served build: benzene
    into cyclohexane into cyclobutane with a cyclopropane over the first
    arrow came back from ``getRxn`` as three ``$MOL`` blocks -- the one over
    the arrow simply was not there, and the scheme came out as
    ``c1ccccc1>>>>>C1CCCCC1>>>>>C1CCC1``, the reagent gone.  ``getSdf`` on the
    same canvas gave four records; ``getCml`` gave two.
    """
    outcome = ketcher.reaction_smiles_from_sdf(
        sdf(placed("c1ccccc1", 0), placed("C1CCCCC1", 14),
            placed("C1CCC1", 28), placed("CC", 8, 5)),
        canvas([(6.0, 10.0, 0.0), (20.0, 24.0, 0.0)]))

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "c1ccccc1>>CC>>>C1CCCCC1>>>>>C1CCC1"
    assert "CC" in outcome["smiles"], "the reagent an RXN file would have lost"


def test_the_editor_is_asked_for_the_records_and_not_for_a_reaction_file():
    """Which is the whole reason a reagent over an arrow survives."""
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import ketcher_panel

    reading = ketcher_panel.read_js("scope-1", "smiles", "auto")

    assert "arrow ? api.getSdf()" in reading, "the records, not a reaction file"
    assert "api.getKet()" in reading, "the arrows are only in there"
    # getRxn is still there, but only for saving a file the user asked for by
    # that name -- never for reading the canvas.
    assert "if(want==='rxn') return api.getRxn();" in reading


def test_records_with_no_arrow_are_a_set_of_structures():
    outcome = ketcher.reaction_smiles_from_sdf(
        sdf(placed("CCO", 0), placed("c1ccccc1", 10)), canvas([]))

    assert outcome["ok"] is True
    assert outcome["smiles"] == "CCO.c1ccccc1"
    assert outcome["steps"] == 0


def test_the_records_reach_the_reader_that_splits_them():
    """``smiles_from_drawing`` is the one call site, and it tells the three
    shapes apart by what they are."""
    body = sdf(placed("c1ccccc1", 0), placed("C1CCCCC1", 14))
    outcome = ketcher.smiles_from_drawing(
        body + ketcher.KET_MARK + canvas([(6.0, 10.0, 0.0)]))

    assert outcome["reaction"] is True
    assert outcome["smiles"] == "c1ccccc1>>>>>C1CCCCC1"


def test_the_grouping_mark_is_taken_back_off_when_it_is_read():
    """A SMILES carries parentheses of its own -- ``CC(=O)O`` -- but never at
    the front, because a branch has to hang off an atom.  So a leading one is
    always a group, and its partner is found by counting."""
    pytest.importorskip("rdkit")

    read = ketcher.parse_reaction_smiles(
        "(C1=CCC=C1.C1=CCC=C1)>>C1CCCCC1>O>>(C1CC1.C1CC1)")

    step, = read["steps"]
    assert step["reactants"] == ["C1=CCC=C1", "C1=CCC=C1"]
    assert step["in"] == ["C1CCCCC1"]
    assert step["out"] == ["O"]
    assert step["products"] == ["C1CC1", "C1CC1"]
    assert ketcher._ungrouped("CC(=O)O") == "CC(=O)O", "not a group"
    assert ketcher._ungrouped("(A.B)") == "A.B"


def test_the_ordinary_grouped_three_field_form_round_trips():
    """Which is what RDKit itself writes, character for character."""
    pytest.importorskip("rdkit")

    text = "(C1=CCC=C1.C1=CCC=C1)>C1CCCCC1.C1CCCCC1>(C1CC1.C1CC1)"
    step, = ketcher.parse_reaction_smiles(text)["steps"]

    assert step["reactants"] == ["C1=CCC=C1", "C1=CCC=C1"]
    assert step["in"] == ["C1CCCCC1", "C1CCCCC1"]
    assert step["out"] == []
    assert step["smiles"] == text
