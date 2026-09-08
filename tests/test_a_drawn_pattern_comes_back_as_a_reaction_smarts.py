"""What Ketcher draws has to arrive as a rule that does what was drawn.

The trap this suite exists for is on the product side.  ``getSmarts`` is
Indigo writing the drawing, and Indigo writes every atom as a query -- which
is right for the left of the arrow and wrong for the right, because a product
template is built rather than matched, and a query atom carries no element
aromaticity.  Taken at its word it turns benzene into tetrahydrofuran instead
of furan, and nothing errors on the way.

So most of what is asserted here is behavioural: the rule that comes out of a
drawing is compared against the hand-written rule it stands for by running
both and comparing what they make.
"""

from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")
Chem = pytest.importorskip("rdkit.Chem")
from rdkit.Chem import rdChemReactions                      # noqa: E402

from delfin.dashboard import chemdarwin_draw as draw        # noqa: E402
from delfin.dashboard import ketcher_smarts as ks           # noqa: E402


# The six rules the tab was built around, as they are written by hand ...
HAND = {
    'benzannulation': "[cH:1]:[cH:2]>>[c:1]1:[c:2]cccc1",
    'furan': "[c&H1:1]:[c&H1:2]:[c:3]:[c:4]:[c:5]:[c:6]"
             ">>[o:2]1[c:3][c:4][c:5][c:6]1",
    'thiophene': "[c&H1:1]:[c&H1:2]:[c:3]:[c:4]:[c:5]:[c:6]"
                 ">>[s:2]1[c:3][c:4][c:5][c:6]1",
    'pyridine': "[c&H1:1]:[c:2]:[c:3]:[c:4]:[c:5]:[c:6]"
                ">>[n:1]:[c:2]:[c:3]:[c:4]:[c:5]:[c:6]",
    'pyrimidine': "[cH:1]1:[c:2]:[c:3]:[cH:4]:[c:5]:[c:6]:1"
                  ">>[n:1]1:[c:2]:[c:3]:[n:4]:[c:5]:[c:6]:1",
    'pyranone': "[cH:1]1[c:2]:[c:3]:[cH:4]:[c:5]:[c:6]1"
                ">>[O:1]=[c:2]1[c:3][c:4][o][c:5][c:6]1",
}

# ... and the same six as Ketcher hands them over, all-query, atomic numbers.
DRAWN = {
    'benzannulation': "[#6;a;H1:1]:[#6;a;H1:2]"
                      ">>[#6;a:1]1:[#6;a:2]:[#6;a]:[#6;a]:[#6;a]:[#6;a]:1",
    'furan': "[#6;a;H1:1]:[#6;a;H1:2]:[#6;a:3]:[#6;a:4]:[#6;a:5]:[#6;a:6]"
             ">>[#8;a:2]1[#6;a:3][#6;a:4][#6;a:5][#6;a:6]1",
    'thiophene': "[#6;a;H1:1]:[#6;a;H1:2]:[#6;a:3]:[#6;a:4]:[#6;a:5]:[#6;a:6]"
                 ">>[#16;a:2]1[#6;a:3][#6;a:4][#6;a:5][#6;a:6]1",
    'pyridine': "[#6;a;H1:1]:[#6;a:2]:[#6;a:3]:[#6;a:4]:[#6;a:5]:[#6;a:6]"
                ">>[#7;a:1]:[#6;a:2]:[#6;a:3]:[#6;a:4]:[#6;a:5]:[#6;a:6]",
    'pyrimidine': "[#6;a;H1:1]1:[#6;a:2]:[#6;a:3]:[#6;a;H1:4]:[#6;a:5]:[#6;a:6]:1"
                  ">>[#7;a:1]1:[#6;a:2]:[#6;a:3]:[#7;a:4]:[#6;a:5]:[#6;a:6]:1",
}

PROBES = ("c1ccccc1", "c1ccc2ccccc2c1", "Cc1ccccc1", "Oc1ccccc1")


def _makes(smarts):
    """Everything this rule makes out of the probes, as canonical SMILES."""
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    assert rxn is not None, f"unreadable: {smarts}"
    made = set()
    for probe in PROBES:
        for group in rxn.RunReactants((Chem.MolFromSmiles(probe),)):
            mol = group[0]
            try:
                Chem.SanitizeMol(mol)
                made.add(Chem.MolToSmiles(mol))
            except Exception:                               # noqa: BLE001
                made.add("!unsanitizable")
    return made


# ---------------------------------------------------------------------------
# The product side
# ---------------------------------------------------------------------------

def test_indigos_own_spelling_would_have_saturated_the_ring():
    """The failure this module exists to prevent, pinned as it actually is."""
    verbatim = _makes(DRAWN['furan'])
    assert "C1CCOC1" in verbatim, "benzene should give the saturated ring here"
    assert "c1ccoc1" not in verbatim, (
        "if the verbatim Indigo spelling ever yields furan on its own, the "
        "product concretiser is no longer needed -- check before deleting it")
    assert "c1ccoc1" in _makes(HAND['furan'])


@pytest.mark.parametrize("name", sorted(DRAWN))
def test_a_drawn_rule_makes_what_the_hand_written_one_makes(name):
    outcome = ks.normalize_reaction_smarts(DRAWN[name])
    assert outcome['ok'], outcome['status']
    assert _makes(outcome['smarts']) == _makes(HAND[name])


@pytest.mark.parametrize("name", sorted(HAND))
def test_a_hand_written_rule_survives_being_normalised(name):
    outcome = ks.normalize_reaction_smarts(HAND[name])
    assert outcome['ok'], outcome['status']
    assert _makes(outcome['smarts']) == _makes(HAND[name])


def test_the_product_keeps_its_hydrogen_count_open():
    """``[c:3]`` and not ``[cH:3]``.

    Every mapped atom is written in brackets and a bracket atom states its
    hydrogens, so letting RDKit work them out pins the product at exactly the
    count the drawing happened to have -- and the rule stops firing at a fused
    position.  Naphthalene is the probe that tells the difference.
    """
    outcome = ks.normalize_reaction_smarts(DRAWN['benzannulation'])
    assert 'H' not in outcome['smarts'].split('>>')[1]
    assert "c1ccc2ccccc2c1" in _makes(outcome['smarts'])


# ---------------------------------------------------------------------------
# Spelling
# ---------------------------------------------------------------------------

def test_atomic_numbers_come_back_as_element_symbols():
    outcome = ks.normalize_reaction_smarts(DRAWN['furan'])
    assert outcome['smarts'] == HAND['furan']


def test_a_shortening_that_keeps_the_meaning_is_taken():
    """Only the spelling changes, and only on the side where it is safe."""
    verbose = "[#6&a&H1:1]:[#6&a&H1:2]>>[o:2]1[c:3][c:4][c:5][c:6]1"
    shorter = ks.prettify(verbose)
    assert shorter == "[c&H1:1]:[c&H1:2]>>[o:2]1[c:3][c:4][c:5][c:6]1"
    assert _makes(shorter) == _makes(verbose)


def test_a_shortening_that_would_change_the_meaning_is_dropped():
    """The guard, not the substitution: behaviour decides, not the spelling.

    ``[#6&a:1]`` on the *product* side is not ``[c:1]`` -- one is a query that
    builds a plain carbon, the other builds an aromatic one -- so the shorter
    spelling is refused here even though the same substitution is taken on the
    reactant side of the test above.
    """
    original = ("[#6&a&H1:1]:[#6&a&H1:2]:[#6&a:3]:[#6&a:4]:[#6&a:5]:[#6&a:6]"
                ">>[#8&a:2]1[#6&a:3][#6&a:4][#6&a:5][#6&a:6]1")
    assert "C1CCOC1" in _makes(original)
    assert "c1ccoc1" in _makes(HAND['furan'])
    assert ks.prettify(original) == original


# ---------------------------------------------------------------------------
# Atom mapping
# ---------------------------------------------------------------------------

def test_a_map_set_by_hand_is_never_moved():
    drawn = ("[#6;a;H1:7]:[#6;a;H1]:[#6;a]:[#6;a]:[#6;a]:[#6;a]"
             ">>[#7;a:7]:[#6;a]:[#6;a]:[#6;a]:[#6;a]:[#6;a]")
    outcome = ks.normalize_reaction_smarts(drawn)
    assert outcome['ok'], outcome['status']
    left, _, right = outcome['smarts'].partition('>>')
    # 7 was set on a carbon that becomes the nitrogen; it has to still join
    # those two and nothing else.
    assert '[c&H1:7]' in left
    assert '[n:7]' in right
    assert 7 not in outcome['auto_maps']


def test_a_drawing_with_no_maps_at_all_says_so_rather_than_guessing_quietly():
    outcome = ks.normalize_reaction_smarts("[#6;a:1]>>[#7;a:1]".replace(':1', ''))
    assert outcome['level'] == 'note'
    assert 'Nothing is mapped' in outcome['status']


def test_a_deleted_atom_is_reported_and_not_refused():
    """``Validate`` calls this an error; the furan rule does it on purpose."""
    outcome = ks.normalize_reaction_smarts(DRAWN['furan'])
    assert outcome['ok']
    assert outcome['deleted_maps'] == [1]
    assert 'is deleted' in outcome['status']


# ---------------------------------------------------------------------------
# Pattern boxes
# ---------------------------------------------------------------------------

def test_a_pattern_comes_back_without_an_arrow_or_a_map():
    outcome = ks.normalize_query_smarts("[#6;a;H1:3]:[#7;a]")
    assert outcome['ok']
    assert ':3' not in outcome['smarts']
    assert outcome['smarts'] == '[c&H1]:n'


def test_a_reaction_drawn_into_a_pattern_box_is_refused():
    outcome = ks.normalize_query_smarts("[c:1]>>[n:1]")
    assert not outcome['ok']
    assert 'arrow' in outcome['status']


def test_a_pattern_drawn_into_the_rule_box_is_refused():
    outcome = ks.normalize_reaction_smarts("[#6;a;H1]:[#6;a;H1]")
    assert not outcome['ok']
    assert 'Reaction Arrow' in outcome['status']


# ---------------------------------------------------------------------------
# inspect: judges, never rewrites
# ---------------------------------------------------------------------------

def test_inspect_leaves_the_line_alone():
    for line in HAND.values():
        assert ks.inspect(line)['ok'], line
    assert not ks.inspect("[c:1]:[c:2]")['ok']
    assert not ks.inspect("[c:1]:[c:2]>>[n:1")['ok']
    assert ks.inspect("[#8]~c~[#8]", reaction=False)['ok']
    assert not ks.inspect("nonsense (", reaction=False)['ok']


# ---------------------------------------------------------------------------
# The line arithmetic that binds a filter to its rule
# ---------------------------------------------------------------------------

def test_rule_lines_reads_the_box_the_way_the_engine_does():
    from delfin.dashboard import tab_chemdarwin as cd
    box = ("\n"
           "[cH:1]>>[c:1]Cl\n"
           "   \n"
           "a name that is not a rule\n"
           "[cH:1]>>[c:1]Br\n")
    assert draw.rule_lines(box) == ["[cH:1]>>[c:1]Cl", "[cH:1]>>[c:1]Br"]
    # and the engine agrees, which is the whole point of reading it twice
    assert len(cd.apply_custom_reaction_iter(
        "c1ccccc1", box, iterations=1, keep_rings=True)) > 0


def test_appending_a_rule_keeps_the_filter_lines_lined_up():
    forbidden = "[N+]\n-"
    assert draw.aligned(forbidden, 4).splitlines() == ["[N+]", "-", "-", "-"]
    # never cut: a filter typed ahead of the rule it is for stays
    assert draw.aligned("a\nb\nc", 1).splitlines() == ["a", "b", "c"]


def test_a_pattern_joins_the_line_of_the_rule_it_belongs_to():
    box = "-\n-\n-"
    box = draw.add_to_line(box, 1, "[N+]")
    assert box.splitlines() == ["-", "[N+]", "-"]
    box = draw.add_to_line(box, 1, "[#8]~c~[#8]")
    assert box.splitlines() == ["-", "[N+];[#8]~c~[#8]", "-"]
    # the same pattern twice is still one pattern
    box = draw.add_to_line(box, 1, "[N+]")
    assert box.splitlines()[1] == "[N+];[#8]~c~[#8]"


def test_a_placeholder_is_replaced_rather_than_appended_to():
    for nothing in ('-', 'none', 'N/A', '.', '*', ''):
        assert draw.add_to_line(nothing, 0, "[N+]") == "[N+]"


# ---------------------------------------------------------------------------
# The road back into the editor
# ---------------------------------------------------------------------------
# Measured against the Ketcher 3.17 the dashboard installs: ``setMolecule``
# reads a SMARTS and drops every atom map doing it, so a rule opened that way
# comes back mapped by guesswork.  An RXN keeps them.

def test_a_rule_goes_back_to_the_editor_as_an_rxn():
    block = ks.rxn_block_from_smarts(HAND['furan'])
    assert block and block.lstrip().startswith('$RXN')
    # the maps are what the whole detour is for
    assert 'V2000' in block
    back = ks.reaction_smarts_from_rxn_block(block)
    assert back and ':1' in back and ':2' in back


def test_the_product_ring_stays_aromatic_on_the_way_out():
    """A query bond has no order to write, and a single bond is what a molfile
    writer puts instead -- which the editor then draws as a saturated ring."""
    block = ks.rxn_block_from_smarts(HAND['furan'])
    back = ks.normalize_reaction_smarts(ks.reaction_smarts_from_rxn_block(block))
    assert back['ok'], back['status']
    assert "c1ccoc1" in _makes(back['smarts'])
    assert "C1CCOC1" not in _makes(back['smarts'])


def test_a_rule_that_would_come_back_different_is_flagged():
    # The benzannulation reactant is two atoms; losing H1 widens it to every
    # aromatic carbon pair, which is a different rule.
    assert not ks.survives_the_editor(HAND['benzannulation'])
    # The ring-swap rules are pinned by their six-atom skeleton either way.
    assert ks.survives_the_editor(HAND['furan'])


def test_the_map_ketcher_writes_in_the_wrong_place_is_put_right():
    """Ketcher writes ``[#6:3;v2]``; SMARTS wants the map last.

    Seen coming back out of the editor after a rule had been opened in it.
    RDKit refuses that spelling, so without this a good drawing reads as an
    unreadable one.
    """
    broken = "[#6:1]:[#6:2]>>[#8:2]1-[#6:3;v2]-[#6:4;v2]-[#6:5;v2]-[#6:6;v2]-1"
    assert Chem.MolFromSmarts(broken.split('>>')[1]) is None
    fixed = ks.repair_atom_maps(broken)
    assert '[#6;v2:3]' in fixed
    assert Chem.MolFromSmarts(fixed.split('>>')[1]) is not None
    # and it goes in on the way from the editor, not as a separate step
    assert ks.normalize_reaction_smarts(broken)['ok']


def test_a_recursive_query_is_left_alone():
    """The ``:`` in ``[$(c:c)]`` is a bond, not a map."""
    line = "[$([#6]:[#6]):1]"
    assert ks.repair_atom_maps(line) == line


# ---------------------------------------------------------------------------
# Kekule drawings against aromatic molecules
# ---------------------------------------------------------------------------
# Ketcher draws benzene with three alternating double bonds and Indigo writes
# that out as drawn; every molecule the rule will meet has been through
# RDKit's aromaticity perception. A '-' does not match an aromatic bond.

KEKULE_RING = ("[#6:1]1-[#6:2]=[#6:3]-[#6:4]=[#6:5]-[#6:6]=1"
               ">>[C:1]1[C:2][C:3][C:4][C:5][C:6]1")
KEKULE_PIECE = "[#6:1](~[#6:2])=[#6:3]~[#6:4]>>[C:1](~[C:2])=[N:3]~[C:4]"


def _runs_on(smarts, smiles):
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    return len(rxn.RunReactants((Chem.MolFromSmiles(smiles),))) if rxn else 0


def test_a_kekule_ring_would_have_matched_nothing():
    """The failure as reported: a right drawing and an empty result."""
    assert _runs_on(KEKULE_RING, "c1ccccc1") == 0


def test_a_kekule_ring_is_widened_and_then_matches():
    out = ks.normalize_reaction_smarts(KEKULE_RING)
    assert out['ok'], out['status']
    assert out['kekule'] == 6
    assert "C1CCCCC1" in _makes(out['smarts'])


def test_the_switch_is_the_whole_trade():
    """Cyclohexa-1,4-diene, drawn with the same two bond types benzene is.

    Nothing in the drawing says which one was meant, so the checkbox decides:
    on, it also matches the aromatic ring, which is what a ChemDarwin seed
    almost always is; off, it matches exactly the diene that was drawn.
    """
    diene = ("[#6:1]1-[#6:2]=[#6:3]-[#6:4]-[#6:5]=[#6:6]-1"
             ">>[C:1]1[C:2][C:3][C:4][C:5][C:6]1")
    off = ks.normalize_reaction_smarts(diene, aromatic=False)
    assert off['kekule'] == 0
    assert _runs_on(off['smarts'], "c1ccccc1") == 0
    assert _runs_on(off['smarts'], "C1=CCC=CC1") > 0

    on = ks.normalize_reaction_smarts(diene, aromatic=True)
    assert on['kekule'] == 6
    assert _runs_on(on['smarts'], "c1ccccc1") > 0
    # and it says so rather than widening quietly
    assert 'also match aromatic' in on['status']


def test_a_fragment_cut_out_of_a_ring_needs_the_switch():
    """An open fragment never perceives as aromatic, so the ring rule cannot
    reach it -- only the atom-by-atom pass can, and that is the checkbox."""
    off = ks.normalize_reaction_smarts(KEKULE_PIECE, aromatic=False)
    assert off['kekule'] == 0
    assert _runs_on(off['smarts'], "c1ccccc1") == 0
    on = ks.normalize_reaction_smarts(KEKULE_PIECE, aromatic=True)
    assert on['kekule'] == 1
    assert _runs_on(on['smarts'], "c1ccccc1") > 0


def test_an_explicitly_aliphatic_bond_is_never_widened():
    """[C] is a carbon drawn as aliphatic; [#6] is one drawn without saying."""
    out = ks.normalize_reaction_smarts("[C:1]-[C:2]>>[C:1]-[N:2]", aromatic=True)
    assert out['kekule'] == 0
    assert _runs_on(out['smarts'], "c1ccccc1") == 0
    assert _runs_on(out['smarts'], "CCCCCC") > 0


def test_a_rule_already_drawn_aromatic_is_untouched():
    out = ks.normalize_reaction_smarts(DRAWN['furan'])
    assert out['kekule'] == 0
    assert out['smarts'] == HAND['furan']


# ---------------------------------------------------------------------------
# The dry run against the seed in the box
# ---------------------------------------------------------------------------

def test_the_rule_is_tried_on_the_seed_before_run_is_pressed():
    good = ks.trial_on_seed(HAND['furan'], "c1ccccc1")
    assert good['level'] == 'ok' and 'product' in good['status']

    never = ks.trial_on_seed(KEKULE_RING, "c1ccccc1")
    assert never['level'] == 'note' and 'no match' in never['status']

    # matches, but the product is drawn aliphatic inside a ring that stays
    # aromatic, and RDKit will not build it
    broken = ks.trial_on_seed(
        "[#6:1](~[#6:2])=,:[#6]~[#6]>>[C:1](~[C:2])=N~C", "c1ccccc1")
    assert broken['level'] == 'note'
    assert 'no product survives' in broken['status']

    assert ks.trial_on_seed(HAND['furan'], '')['status'] == ''


def test_the_trial_names_unmapped_atoms_as_the_cause():
    """The rule matched sixteen times and built nothing, on anthracene.

    Two of its five atoms were left unmapped, so they were deleted and took
    the ring with them -- which is the thing to say, not just that RDKit
    refused what came out.
    """
    drawn = "[#6:1](~[#6:2])(~[#6:3])=,:[#6]~[#6]>>[C:1](~[C:2])(~[C:3])=N~C"
    tried = ks.trial_on_seed(drawn, "c1ccc2cc3ccccc3cc2c1")
    assert tried['level'] == 'note'
    assert 'matches 16x' in tried['status']
    assert 'unmapped and so deleted' in tried['status']

    # mapped through, and drawn aromatic, it builds
    mapped = ("[c:1](~[c:2])(~[c:3])=,:[c:4]~[c:5]"
              ">>[c:1](~[c:2])(~[c:3])~[n:4]~[c:5]")
    assert ks.trial_on_seed(mapped, "c1ccc2cc3ccccc3cc2c1")['level'] == 'ok'


# ---------------------------------------------------------------------------
# The connectivity constraint, drawn
# ---------------------------------------------------------------------------
# Measured through the running editor: a rule whose reacting fragment AND
# whose attachment to the rest of the molecule are both drawn -- the atoms
# carry Ketcher's "Substitution count" query property, which is how a
# ring-fusion carbon is said in a drawing -- comes back as this.

KETCHER_MESO = ("[#6:1;D2](:[#6:3;D3]):[#6:2;D3]"
                ">>[#7:1](:[#6:3;v2]):[#6:2;v2]")


def test_a_drawn_connectivity_constraint_arrives_intact():
    """Substitution count survives Ketcher; the map order does not."""
    # what the editor hands over is not readable as it stands
    assert Chem.MolFromSmarts(KETCHER_MESO.split('>>')[0]) is None
    fixed = ks.repair_atom_maps(KETCHER_MESO)
    assert '[#6;D2:1]' in fixed and '[#6;D3:3]' in fixed

    out = ks.normalize_reaction_smarts(KETCHER_MESO)
    assert out['ok'], out['status']
    left = out['smarts'].split('>>')[0]
    assert 'D2' in left and left.count('D3') == 2, out['smarts']


def test_the_drawn_meso_rule_makes_acridine_and_nothing_else():
    """D3 on both neighbours is what picks anthracene's two meso positions.

    X3 would not: it counts implicit hydrogens, so an aromatic CH is X3 as
    well and the rule would fire at every position.
    """
    out = ks.normalize_reaction_smarts(KETCHER_MESO)
    rxn = rdChemReactions.ReactionFromSmarts(out['smarts'])
    made = set()
    for group in rxn.RunReactants((Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1"),)):
        mol = group[0]
        try:
            Chem.SanitizeMol(mol)
            made.add(Chem.MolToSmiles(mol))
        except Exception:                                   # noqa: BLE001
            pass
    assert made == {"c1ccc2nc3ccccc3cc2c1"}                 # acridine

    loose = ks.normalize_reaction_smarts(
        KETCHER_MESO.replace('D3', 'X3'))
    others = set()
    for group in rdChemReactions.ReactionFromSmarts(
            loose['smarts']).RunReactants(
                (Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1"),)):
        mol = group[0]
        try:
            Chem.SanitizeMol(mol)
            others.add(Chem.MolToSmiles(mol))
        except Exception:                                   # noqa: BLE001
            pass
    assert len(others) > 1, "X3 does not single out the meso positions"


# ---------------------------------------------------------------------------
# A hydrogen drawn as an atom
# ---------------------------------------------------------------------------
# Putting an H on a free position is the natural thing to draw and the one
# thing that cannot work: [#1] matches an explicit hydrogen and a molecule
# RDKit has read carries its hydrogens implicitly.

ANTHRACENE = "c1ccc2cc3ccccc3cc2c1"
H_DRAWN = ("[#6:1](~[#6:2])(~[#6:3])=[#6:4]([#1])~[#6:5]"
           ">>[#6:1](~[#6:2])(~[#6:3])=[#7:4]~[#6:5]")


def test_a_drawn_hydrogen_would_have_matched_nothing():
    target = Chem.MolFromSmiles(ANTHRACENE)
    # the bond widened, so only the hydrogen is left to explain the difference
    with_h = "[#6:1](~[#6:2])(~[#6:3])=,:[#6:4]([#1])~[#6:5]"
    counted = "[#6:1](~[#6:2])(~[#6:3])=,:[#6;H1:4]~[#6:5]"
    assert len(target.GetSubstructMatches(Chem.MolFromSmarts(with_h))) == 0
    assert len(target.GetSubstructMatches(Chem.MolFromSmarts(counted))) == 8
    # it is the implicitness that does it: given real hydrogens it matches
    assert len(Chem.AddHs(target).GetSubstructMatches(
        Chem.MolFromSmarts(with_h))) == 8


def test_a_drawn_hydrogen_is_read_as_a_hydrogen_count():
    out = ks.normalize_reaction_smarts(H_DRAWN)
    assert out['ok'], out['status']
    assert out['drawn_h'] == 1
    assert 'H1' in out['smarts'].split('>>')[0]
    assert 'drawn hydrogen' in out['status']

    rxn = rdChemReactions.ReactionFromSmarts(out['smarts'])
    made = set()
    for group in rxn.RunReactants((Chem.MolFromSmiles(ANTHRACENE),)):
        mol = group[0]
        try:
            Chem.SanitizeMol(mol)
            made.add(Chem.MolToSmiles(mol))
        except Exception:                                   # noqa: BLE001
            pass
    assert "c1ccc2nc3ccccc3cc2c1" in made                   # acridine


def test_absorbing_a_hydrogen_keeps_the_rest_of_the_atom():
    out, count = ks.absorb_hydrogens(
        "[#6;D3:1](~[#6:2])(~[#6:3])=,:[#6;D2:4]([#1])~[#6:5]")
    assert count == 1
    assert '[#6&D2&H1:4]' in out                            # query and map kept
    assert '#1' not in out

    twice, count = ks.absorb_hydrogens("[#6:1]([#1])([#1])~[#6:2]")
    assert count == 2 and '[#6&H2:1]' in twice

    same, count = ks.absorb_hydrogens("[#6:1](~[#6:2])~[#6:3]")
    assert count == 0 and same == "[#6:1](~[#6:2])~[#6:3]"


def test_substitution_count_narrows_where_drawn_neighbours_cannot():
    """Three atoms with D<n> beat five atoms with Any bonds.

    Drawing a neighbour puts that atom *in* the pattern; saying how many
    neighbours an atom has constrains it without adding anything to match.
    """
    def products(smarts):
        out = ks.normalize_reaction_smarts(smarts)
        assert out['ok'], out['status']
        made = set()
        for group in rdChemReactions.ReactionFromSmarts(
                out['smarts']).RunReactants(
                    (Chem.MolFromSmiles(ANTHRACENE),)):
            mol = group[0]
            try:
                Chem.SanitizeMol(mol)
                made.add(Chem.MolToSmiles(mol))
            except Exception:                               # noqa: BLE001
                pass
        return made

    # five atoms, neighbours drawn: catches a peripheral position as well
    assert len(products(H_DRAWN)) == 2
    # three atoms, neighbours counted: the two meso positions and nothing else
    counted = "[#6;D3:1]~[#6;D2;H1:2]~[#6;D3:3]>>[#6:1]~[#7:2]~[#6:3]"
    assert products(counted) == {"c1ccc2nc3ccccc3cc2c1"}


def test_the_three_atom_rule_drawn_in_the_panel_makes_acridine():
    """Driven through the running dashboard, button by button.

    Three atoms joined by Any bonds, each carrying Ketcher's Substitution
    count -- 3 on the outer two, 2 on the middle -- and the middle one going
    to nitrogen. This is what getSmarts handed over, verbatim, malformed map
    order and product valence terms included.
    """
    drawn = "[#6:1;D3]~[#6:2;D2]~[#6:3;D3]>>[#6:1;v2]:[#7:2]:[#6:3;v2]"
    assert Chem.MolFromSmarts(drawn.split('>>')[0]) is None   # as it arrives

    out = ks.normalize_reaction_smarts(drawn)
    assert out['ok'], out['status']
    assert out['smarts'] == "[#6&D3:1]~[#6&D2:2]~[#6&D3:3]>>[c:1][n:2][c:3]"
    assert ks.trial_on_seed(out['smarts'], ANTHRACENE)['level'] == 'ok'

    # and through the tab's own engine, the way Run calls it
    from delfin.dashboard import tab_chemdarwin as cd
    made = cd.apply_custom_reaction_iter(
        ANTHRACENE, out['smarts'], iterations=1, keep_rings=False)
    assert sorted(s for s, _ in made) == ["c1ccc2nc3ccccc3cc2c1"]


def test_a_drawing_with_no_mapping_at_all_says_it_is_guessing():
    """The MCS fills the holes; with nothing drawn it is not filling holes."""
    nothing = "[#6]~[#6]=[#6]~[#6]>>[#6]~[#6](~[#6])=[#7]~[#6]"
    out = ks.normalize_reaction_smarts(nothing)
    assert out['drew_maps'] is False
    assert out['auto_maps']
    assert 'NOTHING was mapped' in out['status']
    assert out['level'] == 'note'

    # holes filled around maps that were drawn read differently
    some = "[#6:7]~[#6]=[#6]~[#6]>>[#6:7]~[#6]=[#7]~[#6]"
    filled = ks.normalize_reaction_smarts(some)
    assert filled['drew_maps'] is True
    assert 'NOTHING was mapped' not in filled['status']


def test_the_atoms_that_change_have_to_be_mapped_too():
    """Reported twice: mapped the context, left the changing atom out.

    An atom mapped on neither side is deleted from the reactant and built
    loose in the product, so the rule matches and makes nothing.
    """
    loose = "[#6:1](~[*:2])(~[*:3])=,:[#6]~*>>[C:1](~[*:2])(~[*:3])=N~*"
    assert ks.trial_on_seed(
        ks.normalize_reaction_smarts(loose)['smarts'], ANTHRACENE
    )['level'] == 'note'

    mapped = ("[#6:1](~[*:2])(~[*:3])=,:[#6:4]~[*:5]"
              ">>[C:1](~[*:2])(~[*:3])=[N:4]~[*:5]")
    out = ks.normalize_reaction_smarts(mapped)
    made = set()
    for group in rdChemReactions.ReactionFromSmarts(
            out['smarts']).RunReactants((Chem.MolFromSmiles(ANTHRACENE),)):
        mol = group[0]
        try:
            Chem.SanitizeMol(mol)
            made.add(Chem.MolToSmiles(mol))
        except Exception:                                   # noqa: BLE001
            pass
    assert "c1ccc2nc3ccccc3cc2c1" in made


# ---------------------------------------------------------------------------
# The map Ketcher drops on a generic atom
# ---------------------------------------------------------------------------
# Measured against Ketcher 3.17: A(:1)~C(:2)~A(:3) >> A(:1)~N(:2)~A(:3) drawn
# and mapped comes back from getSmarts with only the carbon's map. getRxn and
# getKet keep all six, so they are taken back off the drawing.

A_SMARTS = "[*;D3]~[#6:2;D2]~[*;D3]>>[*]~[#7:2]~[*]"


def _a_block():
    """The same rule as an RXN, which is what travels beside the SMARTS."""
    return ks.rxn_block_from_smarts(
        "[*:1]~[#6:2]~[*:3]>>[*:1]~[#7:2]~[*:3]")


def test_a_generic_atom_loses_its_map_in_the_smarts():
    fixed = ks.repair_atom_maps(A_SMARTS)
    left = Chem.MolFromSmarts(fixed.split('>>')[0])
    assert [a.GetAtomMapNum() for a in left.GetAtoms()] == [0, 2, 0]


def test_the_dropped_maps_are_taken_back_off_the_drawing():
    merged, recovered = ks.merge_maps(A_SMARTS, _a_block())
    assert recovered == 4
    left = Chem.MolFromSmarts(merged.split('>>')[0])
    assert [a.GetAtomMapNum() for a in left.GetAtoms()] == [1, 2, 3]

    out = ks.normalize_reaction_smarts(merged)
    assert out['ok'], out['status']
    made = set()
    for group in rdChemReactions.ReactionFromSmarts(
            out['smarts']).RunReactants((Chem.MolFromSmiles(ANTHRACENE),)):
        mol = group[0]
        try:
            Chem.SanitizeMol(mol)
            made.add(Chem.MolToSmiles(mol))
        except Exception:                                   # noqa: BLE001
            pass
    assert made == {"c1ccc2nc3ccccc3cc2c1"}                 # acridine


def test_maps_are_only_taken_from_a_drawing_that_matches():
    """Two readings of one canvas line up; two readings of two do not."""
    other = ks.rxn_block_from_smarts("[c:1]:[c:2]>>[n:1]:[c:2]")
    assert ks.merge_maps(A_SMARTS, other)[1] == 0
    assert ks.merge_maps(A_SMARTS, '')[1] == 0
    # an existing map is never overwritten
    merged, _ = ks.merge_maps(A_SMARTS, _a_block())
    assert ':2' in merged and merged.count(':2') == 2


# ---------------------------------------------------------------------------
# What Ketcher can express, and what survives the trip
# ---------------------------------------------------------------------------
# One drawing per construct, each loaded into the running editor and read back
# with getSmarts. These are the strings it produced -- every quirk in them is
# Ketcher's, not a paraphrase.

KETCHER_SAYS = {
    'generic A':        "[*]~[#6:2]~[*]>>[#6:1]~[#7:2]~[#6:3]",
    'generic AH':       "[*]~[#6:2]~[*]>>[#6:1]~[#7:2]~[#6:3]",
    'generic Q':        "[*;!#6:1]~[#6:2]~[*;!#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'generic X':        "[#6:1]~[F:2,Cl:2,Br:2,I:2,At:2]>>[#6:1]~[#8:2]",
    'generic M':        "[#6:1]~[!#6:2;!#7:2;!#8:2;!F:2;!#15:2;!#16:2;!Cl:2;"
                        "!#34:2;!Br:2;!I:2;!At:2;!He:2;!Ne:2;!Ar:2;!Kr:2;"
                        "!Xe:2;!Rn:2;*]>>[#6:1]~[#8:2]",
    'atom list':        "[#7,#8]~[#6:2]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'charge':           "[#6:1;+]~[#6:2]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'substitution':     "[#6:1;D3]~[#6:2;D2]~[#6:3;D3]>>[#6:1]~[#7:2]~[#6:3]",
    'ring bond count':  "[#6:1]~[#6:2;x2]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'unsaturated':      "[#6:1]~[#6:2;$([*,#1]=,#,:[*,#1])]~[#6:3]"
                        ">>[#6:1]~[#7:2]~[#6:3]",
    'hydrogen count':   "[#6:1]~[#6:2;H]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'bond single/dbl':  "[#6:1]!:;-,=[#6:2]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'bond single/arom': "[#6:1][#6:2]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'bond double/arom': "[#6:1]=,:[#6:2]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'bond triple':      "[#6:1]#[#6:2]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'bond aromatic':    "[#6:1]:[#6:2]:[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
    'bond any':         "[#6:1]~[#6:2]~[#6:3]>>[#6:1]~[#7:2]~[#6:3]",
}


@pytest.mark.parametrize("name", sorted(KETCHER_SAYS))
def test_every_construct_ketcher_offers_can_be_read(name):
    """None of these parse as they arrive; all of them do once repaired."""
    fixed = ks.repair_atom_maps(KETCHER_SAYS[name])
    for side in fixed.split('>>'):
        assert Chem.MolFromSmarts(side) is not None, f"{name}: {side}"


def test_a_map_repeated_on_every_term_collapses_to_one():
    """X and M are written as lists with the map on each branch."""
    halogen = ks.repair_atom_maps(KETCHER_SAYS['generic X'])
    assert '[F,Cl,Br,I,At:2]' in halogen
    metal = ks.repair_atom_maps(KETCHER_SAYS['generic M'])
    assert metal.count(':2') == 2                 # once per side, not per term

    # two different numbers in one atom is a real ambiguity, left to be seen
    assert ks.repair_atom_maps("[#6:1;#7:2]") == "[#6:1;#7:2]"


def test_a_bond_inside_a_recursive_query_is_not_mistaken_for_a_map():
    """Ketcher writes Unsaturated as $([*,#1]=,#,:[*,#1]) -- with a ':' in it."""
    fixed = ks.repair_atom_maps(KETCHER_SAYS['unsaturated'])
    assert '$([*,#1]=,#,:[*,#1])' in fixed         # the recursion is untouched
    assert Chem.MolFromSmarts(fixed.split('>>')[0]) is not None
    assert ks.repair_atom_maps("[$([#6]:[#6]):1]") == "[$([#6]:[#6]):1]"


def test_the_generic_atoms_mean_what_they_say():
    """X only halogens, M only metals -- neither reached anything before."""
    def hits(name, smiles):
        out = ks.normalize_reaction_smarts(ks.repair_atom_maps(KETCHER_SAYS[name]))
        assert out['ok'], out['status']
        rxn = rdChemReactions.ReactionFromSmarts(out['smarts'])
        made = set()
        for group in rxn.RunReactants((Chem.MolFromSmiles(smiles),)):
            mol = group[0]
            try:
                Chem.SanitizeMol(mol)
                made.add(Chem.MolToSmiles(mol))
            except Exception:                       # noqa: BLE001
                pass
        return made

    assert hits('generic X', "Clc1ccccc1") == {"Oc1ccccc1"}
    assert hits('generic X', "Cc1ccccc1") == set()
    assert hits('generic M', "[Fe]C") == {"CO"}
    assert hits('generic M', "Clc1ccccc1") == set()
