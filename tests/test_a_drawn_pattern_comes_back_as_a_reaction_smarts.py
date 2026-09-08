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
