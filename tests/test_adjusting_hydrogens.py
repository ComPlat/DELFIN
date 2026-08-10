"""Adjust H: whether an edit fills the valences it opened, or leaves them.

Filling them is right most of the time -- a carbon that has just lost a bond
should get a hydrogen back.  It is wrong when the open valence is the point:
a radical, a coordination site waiting for a ligand, a fragment about to be
joined to something else.  So it is a switch, and it stands beside the element
that a click draws with.
"""

from __future__ import annotations

import pytest

from delfin.dashboard.molecule_builder import (
    connect_atoms, delete_atoms, displaced_hydrogens, grow_from, place_atom,
    set_bond_order, set_element, structure_from_xyz,
)

ETHANE = (
    "8\nethane\n"
    "C 0 0 0\nC 1.53 0 0\n"
    "H -0.36 1.02 0\nH -0.36 -0.51 0.88\nH -0.36 -0.51 -0.88\n"
    "H 1.89 1.02 0\nH 1.89 -0.51 0.88\nH 1.89 -0.51 -0.88\n"
)

#: Two of them, far enough apart that nothing perceives a bond between the
#: carbons -- so drawing one is a real edit and not a correction.
TWO_METHANES = (
    "10\ntwo methanes\n"
    "C 0 0 0\n"
    "H 0.63 0.63 0.63\nH 0.63 -0.63 -0.63\n"
    "H -0.63 0.63 -0.63\nH -0.63 -0.63 0.63\n"
    "C 4.0 0 0\n"
    "H 4.63 0.63 0.63\nH 4.63 -0.63 -0.63\n"
    "H 3.37 0.63 -0.63\nH 3.37 -0.63 0.63\n"
)


def hydrogens(structure):
    return sum(1 for symbol in structure.symbols if symbol == "H")


def test_a_bond_that_becomes_double_gives_two_hydrogens_back():
    """Ethane to ethene is not a label change: two hydrogens have to go."""
    on = structure_from_xyz(ETHANE, {})
    set_bond_order(on, 0, 1, 2, adjust_h=True)
    assert hydrogens(on) == 4

    off = structure_from_xyz(ETHANE, {})
    set_bond_order(off, 0, 1, 2, adjust_h=False)
    assert hydrogens(off) == 6, "switched off, nothing is taken away"


def test_changing_the_element_re_satisfies_it_only_when_asked():
    on = structure_from_xyz(ETHANE, {})
    set_element(on, 0, "N", adjust_h=True)
    assert hydrogens(on) == 5, "nitrogen takes one fewer than carbon"

    off = structure_from_xyz(ETHANE, {})
    set_element(off, 0, "N", adjust_h=False)
    assert hydrogens(off) == 6


def test_a_drawn_bond_takes_the_place_of_a_hydrogen():
    """Two methanes bonded at the carbons are ethane, not C2H8."""
    on = structure_from_xyz(TWO_METHANES, {})
    assert on is not None and hydrogens(on) == 8

    mapping = connect_atoms(on, 0, 5, adjust_h=True)
    assert hydrogens(on) == 6, "one hydrogen from each end made room"
    assert len(on) == 8
    assert on.order(mapping.get(0, 0), mapping.get(5, 5)) == 1

    off = structure_from_xyz(TWO_METHANES, {})
    assert connect_atoms(off, 0, 5, adjust_h=False) == {}
    assert hydrogens(off) == 8, "switched off, the hydrogens stay"
    assert off.order(0, 5) == 1, "and the bond is drawn either way"


def test_the_hydrogen_that_goes_is_the_one_the_bond_replaces():
    """Not just any hydrogen: the one that was standing in the bond's way."""
    structure = structure_from_xyz(TWO_METHANES, {})
    doomed = displaced_hydrogens(structure, 0, 5)

    assert len(doomed) == 2
    # The first carbon sits at x=0 and its partner at x=4, so the hydrogen it
    # gives up is the one with the largest x; the other carbon's is the one
    # with the smallest.
    first = [j for j in doomed if structure.coords[j][0] < 2.0]
    second = [j for j in doomed if structure.coords[j][0] > 2.0]
    assert len(first) == len(second) == 1, "one from each end"
    assert structure.coords[first[0]][0] > 0, "it pointed at the partner"
    assert structure.coords[second[0]][0] < 4.0, "and so did that one"


def test_a_bond_to_a_metal_displaces_nothing():
    """Ammonia coordinating to platinum keeps all three of its hydrogens."""
    ammine = (
        "5\nammonia and a platinum\n"
        "N 0 0 0\nH 0.94 0 0.33\nH -0.47 0.81 0.33\nH -0.47 -0.81 0.33\n"
        "Pt 0 0 -2.4\n"
    )
    structure = structure_from_xyz(ammine, {})
    assert structure is not None

    assert displaced_hydrogens(structure, 0, 4) == []
    assert connect_atoms(structure, 0, 4, adjust_h=True) == {}
    assert hydrogens(structure) == 3
    assert structure.order(0, 4) == 1


def test_a_bond_that_is_already_there_takes_nothing():
    """Pressing Bond on a bond that exists is a correction, not an edit."""
    structure = structure_from_xyz(ETHANE, {})
    assert structure.order(0, 1) == 1

    assert connect_atoms(structure, 0, 1, adjust_h=True) == {}
    assert hydrogens(structure) == 6


def test_closing_a_ring_costs_both_ends_a_hydrogen():
    """Hexane's two ends joined is cyclohexane: C6H14 becomes C6H12."""
    chain = structure_from_xyz(ETHANE, {})
    # Grow it out to hexane, one carbon at a time.
    tail = 1
    for _ in range(4):
        tail = grow_from(chain, tail, "C")
    assert sum(1 for s in chain.symbols if s == "C") == 6
    assert hydrogens(chain) == 14

    head = 0
    connect_atoms(chain, min(head, tail), max(head, tail), adjust_h=True)
    assert hydrogens(chain) == 12


def test_every_edit_takes_the_switch():
    """One switch, and no edit that quietly ignores it."""
    import inspect

    for edit in (place_atom, grow_from, set_element, delete_atoms,
                 set_bond_order, connect_atoms):
        parameters = inspect.signature(edit).parameters
        assert "adjust_h" in parameters, f"{edit.__name__} ignores the switch"
        assert parameters["adjust_h"].default is True, (
            f"{edit.__name__} must fill valences unless it is told otherwise"
        )


def test_the_toolbar_hands_it_to_every_edit():
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "keep_h = not bool(submit_adjust_h_btn.value)" in handler
    assert handler.count("adjust_h=not keep_h") >= 6, (
        "place, grow, set element, delete and both bond edits"
    )


def test_the_bond_button_reads_the_switch_too():
    """Bond is the other way a bond gets drawn, and it has to make room."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def _edit_bond")[1].split("\n    def ")[0]
    assert "connect_atoms" in handler, "the Bond button ignores Adjust H"
    assert "if connect and bool(submit_adjust_h_btn.value):" in handler, (
        "only a bond being drawn makes room, and only with the switch on"
    )
    assert "_apply_structure(" in handler, (
        "losing an atom is a structural edit and has to be undoable"
    )


# ---------------------------------------------------------------------------
# the element it draws with
# ---------------------------------------------------------------------------
@pytest.fixture
def editor(tmp_path):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    sent = []
    ctx.run_js = sent.append
    _widget, refs = tab_submit.create_tab(ctx)
    refs = dict(refs)
    refs["_js"] = sent        # everything the page was told, in order
    return refs


def test_every_element_is_there(editor):
    choices = list(editor["submit_element_dd"].options)

    assert len(choices) == 118, f"{len(choices)} elements offered"
    for symbol in ("H", "C", "Fe", "Pt", "U", "Og", "Tc", "Xe", "La", "Am"):
        assert symbol in choices, symbol
    assert len(set(choices)) == len(choices), "no element twice"


def test_typing_a_symbol_reaches_the_element_it_names(editor):
    """A native select types ahead by label and matches in the order the
    options stand in.  Sorted by atomic number, pressing P would land on
    phosphorus but I would land on indium, which comes first.  The common ones
    are therefore at the top, and each of them is the first of its letter.
    """
    choices = list(editor["submit_element_dd"].options)

    for symbol in ("C", "H", "N", "O", "S", "P", "F", "B", "I"):
        first = next(c for c in choices if c.upper().startswith(symbol))
        assert first == symbol, (
            f"pressing {symbol} would land on {first}, not on {symbol}"
        )
    # and the rest are still reachable, by typing both letters
    assert "Pd" in choices and "In" in choices


def test_the_switch_stands_beside_the_element_it_draws_with(editor):
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    assert "submit_element_dd, submit_adjust_h_btn," in source, (
        "it belongs behind the element dropdown, where drawing happens"
    )
    assert editor["submit_adjust_h_btn"].value is True, (
        "filling the valences is what is wanted most of the time"
    )


# ---------------------------------------------------------------------------
# and it still has to be there after the next edit
# ---------------------------------------------------------------------------
def _three_methanes():
    rows = []
    for k in range(3):
        x = 4.0 * k
        rows += [f"C {x} 0 0",
                 f"H {x + 0.63} 0.63 0.63", f"H {x + 0.63} -0.63 -0.63",
                 f"H {x - 0.63} 0.63 -0.63", f"H {x - 0.63} -0.63 0.63"]
    return f"{len(rows)}\nthree methanes\n" + "\n".join(rows)


def test_the_first_drawn_bond_survives_the_second(editor):
    """Drawing the second bond rebuilds the picture, and a rebuild draws what
    the distances say -- so the first bond, drawn between two atoms that are
    nowhere near bonding distance, was simply gone from the view. It was still
    in force everywhere else, which is why pressing Optimise brought it back.
    """
    sent = editor["_js"]
    editor["coords_widget"].value = _three_methanes()

    editor["submit_cmd_sync"].value = "bondorder:1:0,5,1"
    rows = [line.split() for line in
            editor["coords_widget"].value.splitlines()[2:] if line.strip()]
    carbons = [i for i, row in enumerate(rows) if row[0] == "C"]
    assert len(rows) == 13, "two hydrogens made room for the first bond"

    sent.clear()
    editor["submit_cmd_sync"].value = f"bondorder:2:{carbons[1]},{carbons[2]},1"

    restored = [s for s in sent if "applyBondEdits" in s]
    assert restored, "the rebuild forgot the bond that was already drawn"
    # On the numbering of the structure that arrived, not the one left behind:
    # a removed hydrogen moves every atom after it.
    assert "[[0, 4, 1]]" in restored[0], restored[0]
    assert editor["editor_state"]["hand_bonds"] == {(0, 4): True, (4, 7): True}


def test_taking_the_edit_back_takes_its_bond_back_too(editor):
    sent = editor["_js"]
    editor["coords_widget"].value = _three_methanes()
    editor["submit_cmd_sync"].value = "bondorder:1:0,5,1"
    rows = [line.split() for line in
            editor["coords_widget"].value.splitlines()[2:] if line.strip()]
    carbons = [i for i, row in enumerate(rows) if row[0] == "C"]
    editor["submit_cmd_sync"].value = f"bondorder:2:{carbons[1]},{carbons[2]},1"

    sent.clear()
    editor["submit_cmd_sync"].value = "undo:3:"

    assert editor["editor_state"]["hand_bonds"] == {(0, 4): True}, (
        "the second bond went with the edit that drew it"
    )
    restored = [s for s in sent if "applyBondEdits" in s]
    assert restored and "[[0, 4, 1]]" in restored[0], (
        "and the first one came back with the structure it belongs to"
    )


def test_the_fullscreen_status_does_not_keep_a_finished_message(editor):
    """The overlay has its own status line, and it was only ever written to,
    never cleared: it kept saying "Quick convert (single structure)..." with
    the finished structure on screen beside it."""
    editor["coords_widget"].value = _three_methanes()
    editor["submit_cmd_sync"].value = "bondorder:1:0,99,1"   # says something
    assert editor["mol_status_fs"].value, "nothing to clear -- test says nothing"

    editor["coords_widget"].value = ETHANE                    # a new view
    assert editor["mol_status_fs"].value == "", (
        "the overlay is still showing the message the small view has dropped"
    )
    assert editor["mol_status"].value == ""


# ---------------------------------------------------------------------------
# what the switch means when you are drawing
# ---------------------------------------------------------------------------
def test_a_second_hydrogen_can_be_put_on_an_oxygen():
    """It could not be, and nothing said why.

    Growing a hydrogen ate the hydrogen already there -- three atoms went in,
    three came out, and the screen did not change. That is right while DELFIN
    is keeping the valences straight, and it is the whole problem when the
    user has switched that off and wants the atom they asked for.
    """
    water = '3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'

    on = structure_from_xyz(water, {})
    grow_from(on, 0, 'H', adjust_h=True)
    assert len(on) == 3, 'with the valence kept, the oxygen stays satisfied'

    off = structure_from_xyz(water, {})
    grown = grow_from(off, 0, 'H', adjust_h=False)
    assert len(off) == 4, 'the hydrogen that was asked for'
    assert grown is not None
    # And the two that were there are untouched.
    assert off.coords[1] == pytest.approx((0.96, 0.0, 0.0))
    assert off.coords[2] == pytest.approx((-0.24, 0.93, 0.0))


def test_it_lands_where_the_hand_let_go():
    """Only the direction of the drag used to be sent, so the atom always
    landed at the length its bond wanted. That is right when the chemistry is
    being kept straight and wrong when the user is placing by hand."""
    water = '3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'
    released = [0.0, -0.80, 0.55]

    off = structure_from_xyz(water, {})
    grown = grow_from(off, 0, 'H', adjust_h=False, at=released)
    assert off.coords[grown] == pytest.approx(tuple(released), abs=1e-6)

    # With the adjustment on it is the bond that decides, as before.
    on = structure_from_xyz(water, {})
    grown = grow_from(on, 0, 'H', adjust_h=True, at=released)
    reach = sum((on.coords[grown][i] - on.coords[0][i]) ** 2 for i in range(3)) ** 0.5
    assert reach == pytest.approx(0.97, abs=0.05), 'an O-H bond length'

    # A click that did not move is a click, not a placement.
    still = structure_from_xyz(water, {})
    grown = grow_from(still, 0, 'H', adjust_h=False, at=[0.02, 0.0, 0.0])
    reach = sum((still.coords[grown][i]) ** 2 for i in range(3)) ** 0.5
    assert reach > 0.3, 'it was put on top of the atom it grew from'


def test_the_page_sends_where_the_hand_let_go():
    from delfin.dashboard import tab_submit
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor = submit_manip_bootstrap_js()
    grow = editor.split("pushCommandToPython(scopeKey, 'grow'")[1][:400]
    assert 'to.x.toFixed' in grow and 'to.z.toFixed' in grow

    source = open(tab_submit.__file__, encoding='utf-8').read()
    handler = source.split('def on_submit_cmd')[1].split('\n    def ')[0]
    assert "verb == 'grow' and len(fields) in (6, 9)" in handler
    assert 'at=landed' in handler


def test_hydrogens_are_not_shuffled_by_an_edit_somewhere_else():
    """They were, on every edit that touched the centre.

    The ideal directions are worked out fresh each time and come out at
    whatever azimuth the construction gives, so an already-correct methyl
    group was matched against a rotated copy of itself and every hydrogen in
    it moved to a sibling's place -- for no reason a user could see, and
    undoing whatever they had positioned by hand.
    """
    from delfin.dashboard.molecule_builder import rearrange_hydrogens

    settled = structure_from_xyz(ETHANE, {})
    before = [tuple(c) for c in settled.coords]
    rearrange_hydrogens(settled, 0)
    assert [tuple(c) for c in settled.coords] == before, 'it shuffled them'

    # A hydrogen put somewhere by hand survives an edit elsewhere.
    dragged = structure_from_xyz(ETHANE, {})
    dragged.coords[2] = (-0.36, 1.30, 0.40)
    place_atom(dragged, 'C', [3.0, 0.0, 0.0])
    assert dragged.coords[2] == pytest.approx((-0.36, 1.30, 0.40))


def test_but_a_centre_that_changed_shape_is_still_rebuilt():
    """Which is what the rearrangement is for: ethene's carbons are trigonal
    planar, not two tetrahedra sharing an edge, and a double bond read back
    off a tetrahedral geometry comes out single."""
    import math

    structure = structure_from_xyz(ETHANE, {})
    set_bond_order(structure, 0, 1, 2)

    hydrogens = [i for i, s in enumerate(structure.symbols) if s == 'H']
    near = [h for h in hydrogens if structure.coords[h][0] < 0.7]
    assert len(near) == 2

    def angle(a, b, c):
        u = [structure.coords[a][i] - structure.coords[b][i] for i in range(3)]
        v = [structure.coords[c][i] - structure.coords[b][i] for i in range(3)]
        nu = sum(x * x for x in u) ** 0.5
        nv = sum(x * x for x in v) ** 0.5
        return math.degrees(math.acos(sum(u[i] * v[i] for i in range(3)) / (nu * nv)))

    assert angle(near[0], 0, near[1]) == pytest.approx(120.0, abs=2.0)
