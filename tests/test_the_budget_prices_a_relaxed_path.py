"""What a dragged geometry costs, and why it cannot be asked as it stands.

The thermal budget compares the energy of the structure under the hand against
the energy it started from, and refuses what the temperature cannot pay for.
Which energy is taken decides whether the editor is a chemistry tool or a
rule about aromatics:

  * The free relaxation undoes the drag.  xtb cannot be told to leave an atom
    alone, so it puts the hydrogen back and reports an intact benzene.  A
    proton could be pulled off with the line underneath saying all was well.

  * A single point on the geometry the mouse leaves overcharges everything
    that forms a bond, because a transition state consists of the rest of the
    molecule rearranging and a rigid geometry rearranges nothing.  Measured on
    butadiene and ethylene with GFN2 at 2.18 A: +124.3 kcal/mol against +3.6
    along the relaxed path.  A Diels-Alder that room temperature really does
    allow then needs 1498 K to get through -- which is how this was found.

  * Held and relaxed is neither.  It cannot undo the drag, because the drag is
    what is held, and it charges only for the distortion that survives once
    the structure has done what it can.  That is what a barrier is.
"""

from __future__ import annotations

import math
import shutil
import time

import pytest

import pathlib

from delfin.dashboard import gfn_optimize as gfn
from delfin.dashboard import structure_editor
from delfin.dashboard.structure_editor import thermal_temperature
from editor_source import EDITOR_SOURCE

_needs_xtb = pytest.mark.skipif(not shutil.which("xtb"), reason="xtb not installed")

#: Ethane, so that a drag has a bond to stretch and hydrogens to brush past.
_ETHANE = """8
ethane
C     0.000000    0.000000    0.000000
C     1.530000    0.000000    0.000000
H    -0.370000    1.020000    0.000000
H    -0.370000   -0.510000    0.880000
H    -0.370000   -0.510000   -0.880000
H     1.900000    0.510000    0.880000
H     1.900000    0.510000   -0.880000
H     1.900000   -1.020000    0.000000
"""


def _the_rise(said):
    """The number the verdict gives as the rise, however it qualifies it.

    It is qualified: a relaxed scan's top is the top of the best path that
    keeps to the coordinate somebody chose, so the saddle lies at or below it
    and the clause carries a "less than or equal".  Read as the word straight
    after a fixed phrase, a wording change turned into a ValueError -- which
    is a check breaking on the thing it was not about -- so it is read as a
    pattern with the qualifier optional, and the clause is now "rise <=4.9
    kcal/mol" rather than a sentence.
    """
    import re as _re

    found = _re.search(r'rise \u2264?(-?\d+(?:\.\d+)?)', said)
    assert found, said
    return float(found.group(1))

def test_a_drag_holds_the_contact_it_changed():
    """One atom under the hand is held to the neighbour it is leaving."""
    held = gfn.contacts_holding(_ETHANE, [2], most=1)
    assert len(held) == 1
    assert held[0]["kind"] == "distance"
    assert held[0]["mode"] == "drag"
    # Its own carbon, not the far one and not a hydrogen.
    assert set(held[0]["atoms"]) == {2, 0}
    assert held[0]["value"] == pytest.approx(1.09, abs=0.02)


def test_heavy_atoms_come_before_hydrogens():
    """A carbon approaching a carbon is a reaction coordinate; H...H is not.

    Holding the wrong one is not a smaller mistake, it is a different answer:
    the hydrogens are exactly what has to bend out of the way to make room for
    a bond, and a held H...H contact stops them.
    """
    # The whole of one methyl group is dragged.  Its nearest contacts to the
    # rest include hydrogen pairs, but the carbon-carbon bond is the one that
    # describes the drag.
    held = gfn.contacts_holding(_ETHANE, [1, 5, 6, 7], most=2)
    assert held, "a fragment leaving the rest has a contact"
    assert set(held[0]["atoms"]) == {1, 0}


def test_a_hydrogen_contact_is_not_held_alongside_a_heavy_one():
    """The second contact is there to stop a dragged group turning away.

    For that it has to be the same sort of approach as the first.  A dragged
    OH is held by its oxygen and its hydrogen left free; both bonds of a
    Diels-Alder are held, because both are carbon meeting carbon.
    """
    held = gfn.contacts_holding(_ETHANE, [1, 5, 6, 7], most=2)
    # Only one heavy-heavy contact exists here, so the second would have to be
    # a hydrogen one -- and is refused rather than taken.
    assert len(held) == 1


def test_two_bonds_forming_at_once_are_both_held():
    """That is a cycloaddition, and both halves of it are the coordinate."""
    apart = (
        "6\ntwo pairs approaching\n"
        "C  0.00 0.00 0.00\nC  1.40 0.00 0.00\nH -0.60 0.90 0.00\n"
        "C  0.00 0.00 2.40\nC  1.40 0.00 2.40\nH -0.60 0.90 2.40\n"
    )
    closer = apart.replace("2.40", "2.20")
    held = gfn.contacts_holding(closer, [3, 4], most=3, was=apart)
    pairs = {frozenset(one["atoms"]) for one in held
             if one["kind"] == "distance"}
    assert frozenset({3, 0}) in pairs
    assert frozenset({4, 1}) in pairs


def test_a_partner_with_something_in_the_way_is_not_a_contact():
    """What used to be an element rule is geometry, and says the same thing.\n\n    In an SN2 the arriving nucleophile closes on the carbon and on the leaving\n    group behind it at exactly the same rate, so both look like the coordinate\n    being driven.  One of them has a carbon in the way.  Held anyway the\n    bromide is nailed down -- C-Br stays at 1.94 to 1.99 A and the price runs\n    to +175 kcal/mol, against 2.01 -> 3.24 A and near zero when only the real\n    contact is held.\n    """
    line = (
        "3\nchloride, carbon, bromide in a row\n"
        "Cl  0.00 0.00 0.00\nC   0.00 0.00 2.50\nBr  0.00 0.00 4.44\n"
    )
    nearer = line.replace("Cl  0.00 0.00 0.00", "Cl  0.00 0.00 0.20")
    held = gfn.contacts_holding(nearer, [0], most=3, was=line)
    assert [frozenset(one["atoms"]) for one in held] == [frozenset({0, 1})]


def test_each_atom_is_held_once_and_to_a_partner_of_its_own():
    """Two contacts on the same atom describe one thing twice.

    They also over-hold: piling contacts on pins the two fragments rigidly to
    each other, which prices the same Diels-Alder geometry at +119 instead of
    +4.8 -- no better than not relaxing at all.
    """
    held = gfn.contacts_holding(_ETHANE, [2, 3], most=2)
    grabbed = [h["atoms"][0] for h in held]
    partners = [h["atoms"][1] for h in held]
    assert len(set(grabbed)) == len(grabbed)
    assert len(set(partners)) == len(partners)


def test_a_hand_grabs_what_hangs_off_what_it_grabs():
    """An atom dragged out of a molecule is not a physical motion.

    Left behind, the hydrogens of a dragged methyl end up at 0.92 and 2.12 A
    from their carbon, one of them 1.87 A from where it started.  On a ring it
    is worse: the axial hydrogen points the way the hand is pushing, so the
    whole push scores as a C-H stretch, that is what gets held, and the
    cyclohexane cannot be flipped.

    Terminal, not hydrogen -- an atom whose only bond is to something being
    dragged comes along whatever it is.
    """
    from delfin.atom_mapping import cov_radius

    rows = [line.split() for line in gfn.atom_lines(_BUTANE)]
    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    radius = [cov_radius(r[0]) for r in rows]
    # The far carbon brings its three hydrogens and nothing else.
    assert gfn.with_their_terminals(where, radius, [3]) == {3, 11, 12, 13}
    # An inner carbon brings its two, and not its carbon neighbours.
    assert gfn.with_their_terminals(where, radius, [1]) == {1, 7, 8}


def test_a_partner_across_the_room_is_a_tether_not_a_hold():
    """Two contacts only when the atom is really held *between* them.

    In an SN2 the two nearest heavy atoms to the arriving nucleophile are the
    carbon it attacks and the leaving group on the far side of it.  Holding
    both nails the leaving group down: measured on Cl- + CH3Br driven in from
    2.5 A, C-Br stays at 1.94 and the price runs to +175 kcal/mol, where with
    one contact C-Br opens to 3.24 and the price stays near zero -- which is
    the reaction.
    """
    sn2 = (
        "3\nchloride, carbon, bromide in a line\n"
        "Cl  0.00 0.00 0.00\n"
        "C   0.00 0.00 2.50\n"
        "Br  0.00 0.00 4.44\n"
    )
    held = gfn.contacts_holding(sn2, [0], most=2)
    assert len(held) == 1
    assert set(held[0]["atoms"]) == {0, 1}


def test_a_structure_moved_as_one_has_nothing_to_hold():
    """Sliding everything sideways changes no internal coordinate.

    It costs nothing and there is nothing to hold it against, so the caller is
    handed no values back and prices that drag some other way.
    """
    assert gfn.contacts_holding(_ETHANE, range(8), most=2) == []
    assert gfn.contacts_holding(_ETHANE, [], most=2) == []


def test_a_value_under_the_hand_is_held_softly_and_the_picture_is_still():
    """A fix rings: it is perturbed every frame and answers too hard.

    Counting how often an atom reverses direction between one answer and the
    next on a slow drag -- which is what a shaking molecule is -- a force
    constant of 20 turns 94 times and 1.0 turns not at all.  It costs almost
    nothing: the same benzene C-H asked to stand at 1.500 A comes back at
    1.499 and +30.2 kcal/mol under 20, and 1.478 and +29.2 under 1.0.  A
    constraint does not have to be stiffer than the bond it holds.
    """
    text = gfn.constraint_input(
        gfn.contacts_holding(_ETHANE, [2], most=1), atoms=8)["text"]
    assert "$constrain" in text
    assert f"force constant={gfn.DRAG_FORCE_CONSTANT}" in text
    assert gfn.DRAG_FORCE_CONSTANT < gfn.FIX_FORCE_CONSTANT
    # xtb counts from one.
    assert "distance: 3, 1," in text


def test_the_hand_sets_the_stiffness_for_the_whole_block():
    """xtb takes one force constant per block, so the drag decides it.

    A value the user is holding exactly goes two hundredths of an Angstrom
    looser while an atom is under the hand, and is exact again the moment it
    is let go -- which is better than a molecule that shakes.
    """
    mine = {"kind": "distance", "atoms": [0, 1], "value": 1.5, "mode": "fix"}
    alone = gfn.constraint_input([mine], atoms=8)["text"]
    assert f"force constant={gfn.FIX_FORCE_CONSTANT}" in alone
    during = gfn.constraint_input(
        [mine] + gfn.contacts_holding(_ETHANE, [2], most=1), atoms=8)["text"]
    assert f"force constant={gfn.DRAG_FORCE_CONSTANT}" in during


def test_letting_go_releases_the_leash_every_time():
    """Not only when the budget had run out.

    A leash is armed the moment an atom is picked up, so a drag that stayed
    well inside the budget still leaves marks and a reach on the page.  Only
    the walled case was released, so the next drag pulled against the last
    one\'s marks -- a molecule that will not move properly, with no line
    anywhere saying why.
    """
    body = EDITOR_SOURCE.split("def _clear_thermal_wall():")[1]
    body = body.split("def _end_gfn_follow")[0]
    assert "_push_thermal_wall(None)" in body
    for line in body.splitlines():
        if "_push_thermal_wall(None)" in line:
            # Not inside an "if it was walled" branch.
            assert line.startswith("        _push"), line


def test_switching_the_budget_off_releases_it_too():
    """Or the last drag\'s leash outlives the switch that made it."""
    off = EDITOR_SOURCE.split("def on_submit_thermal(change):")[1]
    off = off.split("def on_submit_thermal_anchor")[0]
    assert "_clear_thermal_wall()" in off


#: Butane, so that there is a torsion to turn and a bond to pull.
_BUTANE = """14
butane
C     -1.9000    0.3000    0.0000
C     -0.6000   -0.5000    0.0000
C      0.6000    0.4000    0.0000
C      1.9000   -0.4000    0.0000
H     -2.7800   -0.3400    0.0000
H     -1.9400    0.9400    0.8800
H     -1.9400    0.9400   -0.8800
H     -0.5700   -1.1500    0.8800
H     -0.5700   -1.1500   -0.8800
H      0.5700    1.0500    0.8800
H      0.5700    1.0500   -0.8800
H      2.7800    0.2400    0.0000
H      1.9400   -1.0400    0.8800
H      1.9400   -1.0400   -0.8800
"""


#: Cyclohexane, at its GFN2 chair: the case a naive fix would destroy.
#:
#: A ring flip is a torsion about a ring bond, so a rule that refused ring
#: torsions would refuse the most ordinary conformational change there is.
#: The six carbons come first, in ring order.
_CHAIR = """18
chair
C       1.0094     0.9471     0.4872
C      -0.1593     1.4088    -0.3785
C      -1.3684     0.5017    -0.1704
C      -1.0094    -0.9471    -0.4872
C       0.1593    -1.4088     0.3785
C       1.3684    -0.5017     0.1704
H       0.7325     1.0284     1.5410
H       1.8754     1.5911     0.3185
H       0.1368     1.3805    -1.4299
H      -0.4206     2.4404    -0.1324
H      -1.6989     0.5730     0.8687
H      -2.1937     0.8288    -0.8068
H      -0.7325    -1.0284    -1.5410
H      -1.8754    -1.5911    -0.3185
H      -0.1368    -1.3805     1.4299
H       0.4206    -2.4404     0.1324
H       1.6989    -0.5730    -0.8687
H       2.1937    -0.8288     0.8068
"""

#: Methylcyclohexane, at its GFN2 chair, for the flip driven as a gesture.
#:
#: The methyl carbon is first, the ring follows it in order.  A ring carbon is
#: held by two ring bonds and hardly moves when it is pushed at; a substituent
#: on one swings, and that is the grab that flips a chair.
_METHYL_CHAIR = """21
methylcyclohexane
C       2.2574    -0.8403     0.0669
C       0.8095    -0.4320     0.3213
C       0.6739     1.0984     0.3593
C      -0.7444     1.5597     0.0012
C      -1.7687     0.5217     0.4438
C      -1.5931    -0.7709    -0.3618
C      -0.1255    -0.9993    -0.7455
H       2.3650    -1.9221     0.1122
H       2.5814    -0.5053    -0.9166
H       2.9142    -0.3987     0.8136
H       0.5063    -0.8426     1.2903
H       0.9384     1.4543     1.3573
H       1.3833     1.5424    -0.3426
H      -0.8329     1.7017    -1.0776
H      -0.9509     2.5211     0.4751
H      -1.6319     0.3231     1.5084
H      -2.7811     0.9069     0.3106
H      -2.1965    -0.7243    -1.2707
H      -1.9584    -1.6141     0.2280
H       0.0591    -2.0674    -0.8755
H       0.0951    -0.5123    -1.6976
"""

#: The user's own structure: 57 atoms, manganese in an N4O2 salen-like
#: environment with four bromines, charge +1 and closed shell, at its GFN2
#: geometry.  Every ring at a heavy atom either runs through the metal or is
#: an aromatic one, which is what made every grab on it offer a torsion.
_MANGANESE = """57
manganese salen, +1
Br     -1.5068    -0.1914     7.0173
C      -1.0631    -0.6280     5.1985
C      -1.7209    -1.6898     4.5845
C      -1.3950    -1.9779     3.2821
Br     -2.2403    -3.4281     2.3800
C      -0.4475    -1.2529     2.5184
O      -0.2046    -1.5908     1.3111
Mn      0.5853    -0.4375    -0.1060
O       1.3728     0.7302    -1.5092
C       1.0813     1.9530    -1.7345
C       1.9930     2.7843    -2.4306
Br      3.6250     1.9543    -2.9605
C       1.7625     4.1050    -2.7280
C       0.5465     4.6566    -2.3352
Br      0.1911     6.5016    -2.7441
C      -0.4013     3.9074    -1.6767
C      -0.1587     2.5669    -1.3549
C      -1.1890     1.7732    -0.7648
N      -1.0395     0.5721    -0.2883
C      -2.2716    -0.1796    -0.1050
C      -2.5010    -1.1599    -1.2722
C      -1.4025    -2.2150    -1.4176
N      -0.1256    -1.5613    -1.7196
C       0.9692    -2.4498    -2.1414
C       1.7390    -2.9536    -0.9231
N       2.1481    -1.8224    -0.0757
C       3.4399    -1.2100    -0.3992
C       3.8549    -0.2576     0.7241
C       2.8899     0.9281     0.9232
N       1.5590     0.4642     1.2853
C       1.2793     0.5193     2.5543
C       0.2129    -0.1774     3.2004
C      -0.1059     0.1038     4.5341
H      -2.4615    -2.2669     5.1167
H       2.4973     4.6966    -3.2521
H      -1.3447     4.3544    -1.3993
H      -2.1890     2.2127    -0.7951
H      -2.2229    -0.7490     0.8218
H      -3.1326     0.5000    -0.0483
H      -3.4512    -1.6663    -1.0986
H      -2.5844    -0.5874    -2.2018
H      -1.6828    -2.9114    -2.2209
H      -1.3066    -2.7752    -0.4851
H      -0.3066    -0.9310    -2.5053
H       1.6460    -1.8675    -2.7710
H       0.5892    -3.3008    -2.7205
H       1.0808    -3.5926    -0.3293
H       2.6072    -3.5413    -1.2469
H       2.2370    -2.1906     0.8748
H       4.2162    -1.9795    -0.5172
H       3.3498    -0.6557    -1.3360
H       3.9298    -0.8118     1.6652
H       4.8388     0.1511     0.4910
H       3.2975     1.5790     1.7086
H       2.8437     1.4969    -0.0040
H       1.9375     1.0827     3.2203
H       0.4055     0.9092     5.0406
"""


def _moved(xyz, atom, by):
    """*xyz* with one atom shifted, the way a hand shifts it."""
    rows = [line.split() for line in gfn.atom_lines(xyz)]
    for n in range(3):
        rows[atom][n + 1] = f"{float(rows[atom][n + 1]) + by[n]:.6f}"
    return f"{len(rows)}\nmoved\n" + "\n".join(" ".join(r) for r in rows) + "\n"


def test_a_hand_that_turns_is_told_from_one_that_pulls():
    """A distance cannot hold a turn: it is the same all the way round.

    So the budget had no grip on a torsion at all, and neither did the drag.
    Dragging the carbon of a methyl round a chain bond, holding only its C-C,
    the hydrogens were left behind and the C-H distances went to 0.92 and
    2.12 A.  Holding the torsion instead, the same swing keeps C-H at 1.09
    and costs +4.7 kcal/mol over a full turn, which is a single bond.
    """
    # Across the C1-C2 bond is a turn; along it is a pull.
    across = _moved(_BUTANE, 3, (0.0, 0.0, 0.25))
    held = gfn.contacts_holding(across, [3], most=2, was=_BUTANE)
    assert any(one["kind"] == "dihedral" for one in held), held

    along = _moved(_BUTANE, 3, (0.25, 0.0, 0.0))
    held = gfn.contacts_holding(along, [3], most=2, was=_BUTANE)
    assert not any(one["kind"] == "dihedral" for one in held), held


def test_a_torsion_is_built_from_bonded_atoms_not_near_ones():
    """A torsion is a statement about connectivity.

    Taken by proximity the other two atoms came back as hydrogens on a
    neighbouring carbon, and that torsion moved 0.14 degrees while the chain
    swung round by ten.  Picking them by distance is not a rougher version of
    picking them by bonding; it is a different quantity.
    """
    across = _moved(_BUTANE, 3, (0.0, 0.0, 0.25))
    turn = next(one for one in gfn.contacts_holding(
        across, [3], most=2, was=_BUTANE) if one["kind"] == "dihedral")
    # The carbon chain, walked out: 4-3-2-1 in xyz numbering.
    assert turn["atoms"] == [3, 2, 1, 0]


def test_a_turn_is_held_by_the_angle_and_nothing_else():
    """A hand that turns did not ask for a bond length.

    The length it leaves behind -- an atom pushed rigidly through an arc, so
    the bond bent rather than swung -- is not worth defending, and defending
    it squeezed the ring bonds of a cyclohexane from 1.53 to 1.37 A.  Let go,
    they relax to what they should be and the flip costs +11.1 kcal/mol
    against a measured 10.7.
    """
    across = _moved(_BUTANE, 3, (0.0, 0.0, 0.25))
    held = gfn.contacts_holding(across, [3], most=3, was=_BUTANE)
    assert held and {one["kind"] for one in held} == {"dihedral"}


def test_a_bond_being_made_is_held_however_far_off_it_is():
    """The coordinate is the one the hand is changing, not the nearest one.

    Chosen by nearness, one end of a hexadiene pushed at the other for a Cope
    shift picks the carbon it is already bonded to -- whose length the push
    does not change -- so the drive read as a turn, nothing held the approach,
    and the two ends finished 4.19 A apart with the budget reporting +2.0 and
    no objection.  Ranked by what moved, the forming bond is held and the ends
    close to 1.96 A.
    """
    apart = (
        "6\ntwo fragments, one driven at the other\n"
        "C  0.00 0.00 0.00\nC  1.53 0.00 0.00\nH -0.60 0.90 0.00\n"
        "C  0.00 0.00 3.40\nC  1.53 0.00 3.40\nH -0.60 0.90 3.40\n"
    )
    closer = _moved(apart, 3, (0.0, 0.0, -0.20))
    held = gfn.contacts_holding(closer, [3], most=3, was=apart)
    assert held and held[0]["kind"] == "distance"
    # The carbon it is being driven at, not the carbon it is bonded to.
    assert set(held[0]["atoms"]) == {3, 0}


def test_the_first_answer_of_a_drag_turns_rather_than_stretches():
    """It has no previous geometry to compare against, so something has to
    stand in -- and the nearest contact is the wrong thing.

    For an atom in a chain the nearest contact is the bond it hangs on, which
    is the one coordinate a drag must not drive: a hand below what a bond
    holds can only stretch one by a tenth of an angstrom, and one above it
    tears the molecule.  Measured on a 2,4-hexadiene, a chain carbon dragged
    1.75 A: 0.09 A of movement with the pull, and three bonds broken with the
    rigid hand.

    What moves a group is the torsion about the bond one further in, with the
    grabbed atom as its outer end -- turning that swings it through an arc,
    where turning about its own attachment bond would spin the group and
    leave the atom where it was.  Which is what every conformer generator
    drives, and what GOAT freezes bonds in order to leave free.
    """
    held = gfn.contacts_holding(_BUTANE, [3], most=2)
    assert held and all(one["kind"] == "dihedral" for one in held), held
    assert held[0]["atoms"][0] == 3, held


def test_which_way_the_hand_pulls_decides_and_no_topology_can():
    """A torsion carries the grabbed atom one way and one way only.

    It swings it about the axis, so what it can express is the part of the
    hand's motion square to both the axis and the arm.  A hand pulling any
    other way is pushing on a coordinate that cannot say what it is doing:
    the drag moves almost nothing and spends a great deal not moving it.

    That this is the question, and not "is the axis inside a ring", is
    settled by one molecule.  The same end carbon of the same butane, dragged
    across the chain and dragged along it, wants a torsion in the first case
    and a bond in the second -- and the two grabs have identical topology, so
    no rule about rings, bonds or fragments can tell them apart.  The
    direction can: the swing about the C2-C3 axis is square to the plane the
    anti chain lies in, so a hand across the chain has all of it -- 1.42 A
    per radian along the pull against 0.16 for the bond -- and a hand along
    the chain has none of it, 0.00 against 0.89.

    Asked blind -- with the geometry the drag started from withheld -- both
    come back as the same torsion, which is the bug this is about.
    """
    across = _moved(_BUTANE, 3, (0.0, 0.0, 0.25))
    held = gfn.contacts_holding(across, [3], most=2, was=_BUTANE, opening=True)
    assert held[0]["kind"] == "dihedral", held
    assert held[0]["atoms"] == [3, 2, 1, 0], held

    along = _moved(_BUTANE, 3, (0.25, 0.0, 0.0))
    held = gfn.contacts_holding(along, [3], most=2, was=_BUTANE, opening=True)
    assert held[0]["kind"] == "distance", held
    assert set(held[0]["atoms"]) == {3, 2}, held

    # The same two grabs with nothing said about where the hand came from.
    assert gfn.contacts_holding(across, [3], most=2)[0]["kind"] == "dihedral"
    assert gfn.contacts_holding(along, [3], most=2)[0]["kind"] == "dihedral"


def test_a_ring_torsion_is_offered_by_the_gesture_and_not_by_the_ring():
    """The trap, and why refusing ring torsions is not the fix.

    A chelate's rings run through the metal, and turning about a bond in one
    of them deforms the ring rather than moving the atom -- which is what
    made a drag on the manganese complex feel dead and expensive at once.
    But a cyclohexane ring flip *is* a torsion about a ring bond, and it is
    the most ordinary conformational change there is.  A rule that refused
    ring torsions would break the case the editor is most wanted for.

    So the same ring carbon is grabbed twice.  Pushed out of the ring plane
    the ring torsions carry it -- 1.19 A per radian along the push, against
    0.16 for the bond it hangs on -- and one is offered, which is the pucker
    and the flip.  Pulled radially outward, in the plane, the same torsions
    carry 0.59 and the bond carries 0.62, so the bond is offered: pulling an
    atom out of a ring really is a bond stretch, and the price the wall then
    quotes is the price of one.
    """
    where = [tuple(float(one) for one in line.split()[1:4])
             for line in gfn.atom_lines(_CHAIR)]
    ring = list(range(6))
    middle = [sum(where[n][k] for n in ring) / 6 for k in range(3)]
    # The ring's own axis, summed over its edges: taken from any single pair
    # of atoms a chair gives a skewed one.
    axis = [0.0, 0.0, 0.0]
    for a in range(6):
        one = [where[ring[a]][k] - middle[k] for k in range(3)]
        two = [where[ring[(a + 1) % 6]][k] - middle[k] for k in range(3)]
        axis[0] += one[1] * two[2] - one[2] * two[1]
        axis[1] += one[2] * two[0] - one[0] * two[2]
        axis[2] += one[0] * two[1] - one[1] * two[0]
    span = math.sqrt(sum(v * v for v in axis))
    axis = [v / span for v in axis]

    out = [where[0][k] - middle[k] for k in range(3)]
    span = math.sqrt(sum(v * v for v in out))
    out = [v / span for v in out]

    def grabbed(unit, by=0.25):
        rows = [line.split() for line in gfn.atom_lines(_CHAIR)]
        for k in range(3):
            rows[0][k + 1] = f"{where[0][k] + by * unit[k]:.6f}"
        return (f"{len(rows)}\npushed\n"
                + "\n".join(" ".join(r) for r in rows) + "\n")

    held = gfn.contacts_holding(grabbed(axis), [0], most=3, was=_CHAIR,
                                opening=True)
    assert held[0]["kind"] == "dihedral", held
    # And it is a bond of the ring itself that is turned about, which is the
    # whole point: the ring is not being avoided, it is being puckered.
    assert set(held[0]["atoms"][1:3]) < set(ring), held

    held = gfn.contacts_holding(grabbed(out), [0], most=3, was=_CHAIR,
                                opening=True)
    assert held[0]["kind"] == "distance", held
    assert set(held[0]["atoms"]) < set(ring), held


def test_a_methyl_leaving_an_ethane_still_has_no_turn_to_offer():
    """A bond is only the wrong answer where a better one exists.

    There is no torsion between the two halves of an ethane, so nothing is
    offered and the C-C bond stands -- which is the coordinate that really
    does describe pulling one methyl off the other.  The direction the hand
    is pulling makes no difference to that, because there is nothing to
    choose between.
    """
    from delfin.atom_mapping import cov_radius

    rows = [line.split() for line in gfn.atom_lines(_ETHANE)]
    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    radius = [cov_radius(r[0]) for r in rows]
    held = gfn.with_their_terminals(where, radius, [1])
    assert gfn.turn_for(where, radius, [1], held) == []
    assert gfn.turn_for(where, radius, [1], held,
                        pulled=[1.0, 0.0, 0.0]) == []
    # Which is what the caller is left holding, whole methyl or single atom.
    assert set(gfn.contacts_holding(
        _ETHANE, [1, 5, 6, 7], most=2)[0]["atoms"]) == {1, 0}


def test_the_opening_decision_is_geometry_and_costs_nothing():
    """It runs on every grab, so it has to be free beside the calculation.

    Measured on the manganese complex, 57 atoms, every heavy atom grabbed in
    turn: 1.7 ms to decide what the hand is driving, against 2.97 s for the
    twenty relaxation cycles that answer the same frame under GFN2 -- about
    six hundredths of one per cent, and it is asked once per drag rather than
    once per answer.  Nothing in it evaluates an energy.
    """
    where = [tuple(float(one) for one in line.split()[1:4])
             for line in gfn.atom_lines(_MANGANESE)]
    heavy = [n for n, line in enumerate(gfn.atom_lines(_MANGANESE))
             if line.split()[0] != "H"]
    assert len(heavy) == 33
    began = time.perf_counter()
    for atom in heavy:
        rows = [line.split() for line in gfn.atom_lines(_MANGANESE)]
        rows[atom][1] = f"{where[atom][0] + 0.15:.6f}"
        put = (f"{len(rows)}\npushed\n"
               + "\n".join(" ".join(r) for r in rows) + "\n")
        gfn.contacts_holding(put, [atom], most=3, was=_MANGANESE,
                             opening=True)
    each = (time.perf_counter() - began) / len(heavy)
    # A wide margin: what is being defended is that this is nothing beside
    # one xtb answer, not a particular millisecond count.
    assert each < 0.05, f"{1000 * each:.1f} ms to decide one grab"


def test_the_complex_stops_driving_what_the_hand_cannot_express():
    """The finding, and what is left of it.

    Measured on the user's own structure -- 57 atoms, manganese in an N4O2
    salen-like environment with four bromines, charge +1 -- every heavy atom
    grabbed in turn and pulled three ways: outward from the middle of the
    molecule and along the two directions square to that.  Ninety-nine grabs.

    What the hand was made to drive, in Angstrom the grabbed atom moves per
    unit of the coordinate *along the pull*: a median of 0.542 before and
    0.955 after.  The grabs where what was driven could express almost none
    of the hand -- under 0.2 A -- fall from 19 to 6, and 87 torsions become
    69, the difference being grabs where the honest coordinate was the bond
    all along.  The chosen coordinate is never worse than the one that stood
    before, which is not luck: the same candidates are weighed and the best
    is taken, where before the first one the walk came to was.

    What that is worth when the drags are actually run -- the same 99 grabs,
    ten answers of 0.15 A each under GFN-FF, counting only the ground the
    hand kept while the structure was still inside a 298 K budget: the median
    goes from 0.034 to 0.100 of what the cursor asked for, the grabs that
    move at all (a tenth or better) from 34 to 51, and the dead ones (under a
    fiftieth) from 42 to 22.  Better on 37, unchanged on 51, worse on 11.
    Pulled outward alone, which is how the finding was first measured, 29 of
    33 drove a torsion and 5 moved; now 19 drive a torsion and 14 move.

    The complex is still a rigid complex: a tenth of the cursor is what a
    chelated salen gives, and the change is that the hand now spends its
    budget on a coordinate the drag can express instead of on one it cannot.
    """
    from delfin.atom_mapping import cov_radius

    lines = gfn.atom_lines(_MANGANESE)
    where = [tuple(float(one) for one in line.split()[1:4]) for line in lines]
    radius = [cov_radius(line.split()[0]) for line in lines]
    middle = [sum(one[k] for one in where) / len(where) for k in range(3)]

    def ways(atom):
        out = [where[atom][k] - middle[k] for k in range(3)]
        span = math.sqrt(sum(v * v for v in out)) or 1.0
        out = [v / span for v in out]
        other = [0.0, 0.0, 1.0] if abs(out[2]) < 0.9 else [1.0, 0.0, 0.0]
        one = [out[1] * other[2] - out[2] * other[1],
               out[2] * other[0] - out[0] * other[2],
               out[0] * other[1] - out[1] * other[0]]
        span = math.sqrt(sum(v * v for v in one)) or 1.0
        one = [v / span for v in one]
        two = [out[1] * one[2] - out[2] * one[1],
               out[2] * one[0] - out[0] * one[2],
               out[0] * one[1] - out[1] * one[0]]
        return out, one, two

    def carries(text, held, unit):
        put = [tuple(float(one) for one in line.split()[1:4])
               for line in gfn.atom_lines(text)]
        atoms = held[0]["atoms"]
        if held[0]["kind"] == "distance":
            i, j = atoms
            axis = [put[i][k] - put[j][k] for k in range(3)]
            span = math.sqrt(sum(v * v for v in axis)) or 1.0
            return gfn._carries([v / span for v in axis], unit)
        return gfn._carries(gfn._swing(put, *atoms[:3]), unit)

    heavy = [n for n, line in enumerate(lines) if line.split()[0] != "H"]
    got, turns = [], 0
    for atom in heavy:
        for unit in ways(atom):
            rows = [line.split() for line in lines]
            for k in range(3):
                rows[atom][k + 1] = f"{where[atom][k] + 0.15 * unit[k]:.6f}"
            put = (f"{len(rows)}\npushed\n"
                   + "\n".join(" ".join(r) for r in rows) + "\n")
            held = gfn.contacts_holding(put, [atom], most=3, was=_MANGANESE,
                                        opening=True)
            got.append(carries(put, held, unit))
            turns += held[0]["kind"] == "dihedral"
    assert len(got) == 99
    got.sort()
    # Margins, because a covalent radius or a cursor step could move a grab
    # from one side of the count to the other; the sizes are what is claimed.
    assert got[len(got) // 2] > 0.85, got[len(got) // 2]
    assert sum(1 for one in got if one < 0.2) <= 10
    assert 55 <= turns <= 80, turns


def _value_in(xyz, entry):
    """One held coordinate, read off a geometry -- the editor's own reader."""
    here = [tuple(float(one) for one in line.split()[1:4])
            for line in gfn.atom_lines(xyz)]
    at = [here[int(i)] for i in entry["atoms"]]
    if entry["kind"] == "distance":
        return math.dist(at[0], at[1])
    return gfn._dihedral(at, 0, 1, 2, 3)


def _gesture(xyz, atom, unit, frames=6, step=0.15, method="gfn2"):
    """A whole drag, driven the way the editor's follow loop drives one.

    The cursor moves a step between answers, each answer is what the next
    step starts from, and the coordinate the first answer chose sticks --
    which is thermal_was, thermal_turn and thermal_holding in the editor.
    """
    start = [tuple(float(one) for one in line.split()[1:4])
             for line in gfn.atom_lines(xyz)]
    rest = [n for n in range(len(start)) if n != atom]
    was, turning, holding, opening, settled = xyz, None, None, True, xyz
    trail = []
    for k in range(1, frames + 1):
        here = [tuple(float(one) for one in line.split()[1:4])
                for line in gfn.atom_lines(settled)]
        here[atom] = tuple(start[atom][n] + k * step * unit[n]
                           for n in range(3))
        rows = [line.split() for line in gfn.atom_lines(xyz)]
        for row, point in zip(rows, here):
            row[1:4] = [f"{one:.6f}" for one in point]
        current = (f"{len(rows)}\ndragged\n"
                   + "\n".join(" ".join(r) for r in rows) + "\n")
        contacts = gfn.contacts_holding(current, [atom], most=3, was=was,
                                        turning=turning, holding=holding)
        if contacts and opening:
            opening = False
            fresh = gfn.contacts_holding(current, [atom], most=3, was=was,
                                         opening=True)
            if fresh and str(fresh[0].get("kind")) == "dihedral":
                contacts = fresh
        pushes = gfn.as_pushes(contacts, was, 0.4 * gfn.A_BOND_HOLDS,
                               value_of=_value_in)
        answer = gfn.relax_steps(current, method=method, cycles=20,
                                 timeout=120.0, constraints=pushes)
        assert answer.get("ok"), answer.get("status")
        settled = gfn.settle_onto(answer["xyz"], current, rest)
        holding = [dict(one) for one in pushes]
        turned = [one for one in pushes if str(one.get("kind")) == "dihedral"]
        if turned:
            turning = list(turned[0]["atoms"])
        was = settled
        now = [tuple(float(one) for one in line.split()[1:4])
               for line in gfn.atom_lines(settled)]
        trail.append(sum((now[atom][n] - start[atom][n]) * unit[n]
                         for n in range(3)))
    return contacts, settled, trail


@_needs_xtb
def test_a_ring_flip_is_still_a_drag_that_works():
    """The case a rule about ring bonds would have destroyed, driven for real.

    A methyl on a chair, pushed back through the ring plane, six answers of
    0.15 A.  The coordinate the drag settles on is a torsion about a *ring*
    bond -- it has to be, that is what a flip is -- and the methyl carbon
    follows: 0.04, 0.10, 0.18, 0.35, 0.55, 0.66 A against 0.90 asked, so
    about three quarters of the last step's worth of ground, and the ring's
    torsions come out with three of the six signs reversed, which is the
    ring inverted rather than merely strained.  It costs +1.1 kcal/mol at its
    worst and +0.0 where it lands, against the 10.7 the flip really is.

    Identical before and after this change, to the third decimal: the carry
    ranks the two ring torsions at the grabbed methyl 1.33 and 1.11, and the
    one it takes is the one the walk used to reach first.  What the change
    does here is leave the case alone, which is the point of it.

    A ring *carbon* is a different grab and does not flip anything: it is
    held by two ring bonds and moves 0.01 A for the same push.  That is the
    chemistry and not the editor.
    """
    where = [tuple(float(one) for one in line.split()[1:4])
             for line in gfn.atom_lines(_METHYL_CHAIR)]
    ring = list(range(1, 7))
    middle = [sum(where[n][k] for n in ring) / 6 for k in range(3)]
    axis = [0.0, 0.0, 0.0]
    for a in range(6):
        one = [where[ring[a]][k] - middle[k] for k in range(3)]
        two = [where[ring[(a + 1) % 6]][k] - middle[k] for k in range(3)]
        axis[0] += one[1] * two[2] - one[2] * two[1]
        axis[1] += one[2] * two[0] - one[0] * two[2]
        axis[2] += one[0] * two[1] - one[1] * two[0]
    span = math.sqrt(sum(v * v for v in axis))
    axis = [v / span for v in axis]
    # Back through the plane, not further out of it: one way is the flip and
    # the other is a chair being made more of a chair, which is stiff.
    lean = sum((where[0][k] - middle[k]) * axis[k] for k in range(3))
    unit = [-math.copysign(1.0, lean) * v for v in axis]

    def torsions(text):
        put = [tuple(float(one) for one in line.split()[1:4])
               for line in gfn.atom_lines(text)]
        return [gfn._dihedral(put, ring[n], ring[(n + 1) % 6],
                              ring[(n + 2) % 6], ring[(n + 3) % 6])
                for n in range(6)]

    held, settled, trail = _gesture(_METHYL_CHAIR, 0, unit)
    assert held[0]["kind"] == "dihedral", held
    # About a bond of the ring itself.
    assert set(held[0]["atoms"][1:3]) < set(ring), held
    # The methyl goes where it is pushed, and goes on going.  Measured at
    # 0.04, 0.10, 0.18, 0.35, 0.55, 0.66 A, so both of these have room: xtb's
    # optimiser is not bit-reproducible under threading and what is claimed
    # here is the shape of the answer, not any one of its numbers.
    assert trail[-1] > 0.4, trail
    assert all(after > before - 0.05
               for before, after in zip(trail, trail[1:])), trail
    # And the ring is a different ring at the end of it.
    before, after = torsions(_METHYL_CHAIR), torsions(settled)
    flipped = sum(1 for one, two in zip(before, after) if one * two < 0)
    assert flipped >= 2, (before, after)


def test_the_atom_under_the_hand_stays_under_the_hand():
    """And when the hold is loose, the line says so instead of shaking.

    Drawn as it arrives, an answer puts the grabbed atom beside the cursor --
    0.102 A through the transition region of a backside attack -- the page
    draws it back on the next mouse move, and at seven answers a second that
    is a molecule that shakes.  Laid onto the whole structure it is no better:
    every other atom then moves the same 0.1 A instead, 0.210 against 0.108.

    So the atom goes back where the cursor has it, as before, and the
    divergence is reported rather than displayed.  On a palladium pushed at
    head on it reached 0.7 A and was the only sign anything was wrong.
    """
    assert "str(one.get('kind')) == 'dihedral'" in EDITOR_SOURCE
    assert "_gfn.hold_atoms_at(" in EDITOR_SOURCE
    assert "slipped > _SLIP_LOOSE" in EDITOR_SOURCE
    assert "The hold is loose here" in EDITOR_SOURCE
    # And the geometry a turn is measured against is the one this answer
    # handed back, not the one it was handed: against the latter the
    # difference holds the relaxation as well as the hand, and on a ring
    # being puckered the relaxation is the larger of the two.
    assert "state['thermal_was'] = settled" in EDITOR_SOURCE


def test_only_the_atoms_the_hand_took_hold_of_make_a_statement():
    """The travellers come along, but they do not each claim a coordinate.

    A methyl's hydrogens riding on their carbon each claimed one, and the two
    spurious distances they claimed pinned a cyclohexane so it could not be
    flipped -- +55.4 kcal/mol for a ring pucker that costs +6.9 when only the
    carbon speaks.  Three coordinates on one atom pin it outright: a Claisen
    whose new bond was still 4.3 A away then priced at +133.
    """
    across = _moved(_BUTANE, 3, (0.0, 0.0, 0.25))
    held = gfn.contacts_holding(across, [3], most=3, was=_BUTANE)
    owners = [one["atoms"][0] for one in held]
    assert owners == [3], owners

    # Two atoms under the hand make two statements, which is a cycloaddition.
    apart = (
        "6\ntwo pairs approaching\n"
        "C  0.00 0.00 0.00\nC  1.40 0.00 0.00\nH -0.60 0.90 0.00\n"
        "C  0.00 0.00 2.40\nC  1.40 0.00 2.40\nH -0.60 0.90 2.40\n"
    )
    closer = apart.replace("2.40", "2.20")
    held = gfn.contacts_holding(closer, [3, 4], most=3, was=apart)
    assert sorted(one["atoms"][0] for one in held) == [3, 4]


def test_the_follow_prices_the_relaxation_it_ran():
    """One calculation, not two: the follow and the price are the same run."""
    assert "if contacts and pull is None:\n                        priced = outcome" in EDITOR_SOURCE
    assert "_gfn.contacts_holding(" in EDITOR_SOURCE
    # And the held contacts go into the run that follows the drag, or the
    # relaxation would still be free to undo it.
    assert "constraints=keeping + contacts" in EDITOR_SOURCE


def test_a_held_follow_is_given_the_cycles_to_get_anywhere():
    """A relaxed scan cut short reads as a barrier that is not there.

    Measured on a Diels-Alder approach whose relaxed cost is +3.6 kcal/mol:
    three cycles say +40.0, six say +18.9, twenty say +8.3.  Each of the first
    two is a wall in the wrong place.
    """
    assert "_THERMAL_FOLLOW_CYCLES" in EDITOR_SOURCE
    line = next(one for one in EDITOR_SOURCE.splitlines()
                if one.strip().startswith("_THERMAL_FOLLOW_CYCLES ="))
    assert int(line.split("=")[1].strip()) >= 20


def test_more_cores_are_not_more_speed(monkeypatch):
    """The ladder is by size and it ends, because xtb's gain does.

    Measured at twenty GFN2 cycles: eighteen atoms are fastest on four cores,
    a hundred and eighty-five on sixteen, and sixty-four cores are three times
    worse than four.  A machine with 384 of them does not change that.

    Asked of a big enough machine, or the ladder is invisible behind the
    clamp: this first said 8 for 49 atoms on a host with four cores, where
    every rung reads as four.
    """
    monkeypatch.setattr(gfn.os, "cpu_count", lambda: 64)
    assert gfn.interactive_cores(18) == 4
    assert gfn.interactive_cores(49) == 8
    assert gfn.interactive_cores(101) == 16
    assert gfn.interactive_cores(5000) == 16
    # No size named: an interactive drag is mostly a small structure.
    assert gfn.interactive_cores() == gfn.INTERACTIVE_CORES


def test_a_small_machine_is_never_asked_for_more_than_it_has(monkeypatch):
    """The ladder is what xtb would like, not what the host can give."""
    monkeypatch.setattr(gfn.os, "cpu_count", lambda: 2)
    assert gfn.interactive_cores(5000) == 2
    assert gfn.interactive_cores(18) == 2


def test_the_cores_can_still_be_said_outright(monkeypatch):
    """Whoever sets the variable has said what they want."""
    monkeypatch.setattr(gfn.os, "cpu_count", lambda: 64)
    monkeypatch.setenv("DELFIN_GFN_CORES", "2")
    assert gfn.interactive_cores(5000) == 2


@_needs_xtb
def test_a_hydrogen_pulled_off_a_ring_stays_far_too_expensive():
    """The relaxation must not be able to undo the drag.

    Free, it puts the hydrogen back at 1.09 A and reports -1.2 kcal/mol
    whatever the hand does.  Held, there is nowhere for the ring to go and the
    price is the bond -- which no room temperature pays for.
    """
    from delfin.dashboard.structure_editor import thermal_ceiling

    benzene = (
        "12\nbenzene\n"
        "C   1.3970  0.0000  0.0000\nC   0.6985  1.2098  0.0000\n"
        "C  -0.6985  1.2098  0.0000\nC  -1.3970  0.0000  0.0000\n"
        "C  -0.6985 -1.2098  0.0000\nC   0.6985 -1.2098  0.0000\n"
        "H   2.4810  0.0000  0.0000\nH   1.2405  2.1487  0.0000\n"
        "H  -1.2405  2.1487  0.0000\nH  -2.4810  0.0000  0.0000\n"
        "H  -1.2405 -2.1487  0.0000\nH   1.2405 -2.1487  0.0000\n"
    )
    settled = gfn.optimize_with_gfn(benzene, "gfn2")
    assert settled.get("ok"), settled.get("status")
    rest = settled["xyz"]
    here = gfn.optimize_with_gfn(rest, "gfn2", optimise=False)["energy"]

    rows = [line.split() for line in gfn.atom_lines(rest)]
    rows[6][1] = f"{float(rows[6][1]) + 1.1:.6f}"        # the hydrogen, pulled
    pulled = f"{len(rows)}\npulled\n" + "\n".join(" ".join(r) for r in rows) + "\n"

    held = gfn.contacts_holding(pulled, [6], most=1)
    assert set(held[0]["atoms"]) == {6, 0}
    answer = gfn.relax_steps(pulled, method="gfn2", cycles=20,
                             timeout=120.0, constraints=held)
    assert answer.get("ok"), answer.get("status")
    spent = (answer["energy"] - here) * 627.5094740631
    assert spent > thermal_ceiling(298.15, 3600.0), (
        f"a hydrogen dragged 1.1 A off a ring priced at {spent:+.1f} kcal/mol")

    # And the free relaxation is what that has to be measured against: it says
    # the same drag costs nothing, because it has quietly undone it.
    free = gfn.relax_steps(pulled, method="gfn2", cycles=20, timeout=120.0)
    assert (free["energy"] - here) * 627.5094740631 < spent / 2


# --- The scan: a coordinate walked rather than a hand followed --------------


def test_the_ceiling_turns_around():
    """A drag asks what may happen at a temperature; a scan asks for the
    temperature.

    "+29 kcal/mol" means nothing until it is "and that wants 385 K within the
    hour", which is the difference between a number and an experiment.  The
    measured barriers land where the chemistry does: a cyclohexane flip at
    +11 wants 150 K, a [1,5]-hydrogen shift at +25.5 wants 340 K, a Claisen at
    +29 wants 385 K, a Cope at +34 wants 450 K.
    """
    from delfin.dashboard.structure_editor import (
        thermal_ceiling, thermal_temperature)

    for kcal in (5.0, 11.0, 22.3, 29.0, 34.0, 120.0):
        T = thermal_temperature(kcal, 3600.0)
        assert T is not None
        # It is the temperature whose ceiling is that barrier.
        assert thermal_ceiling(T, 3600.0) == pytest.approx(kcal, abs=0.05)

    assert thermal_temperature(22.3, 3600.0) == pytest.approx(298.15, abs=1.0)
    # Nothing under 5000 K crosses this, and saying so beats printing a number.
    assert thermal_temperature(2000.0, 3600.0) is None
    assert thermal_temperature("not a barrier", 3600.0) is None


def test_a_fast_crossing_is_said_in_units_it_fits():
    """An open path is the ordinary case, and it was reported as "about
    4.18e-06 s" -- a number in the wrong clothes."""
    body = EDITOR_SOURCE.split("def _thermal_wait(")[1].split("def ")[0]
    for unit in ("'ps'", "'ns'", "'us'", "'ms'"):
        assert unit in body, unit


def test_the_scan_walks_every_armed_leg_together():
    """One coordinate answers a rotation or a stretch; a concerted step needs
    two.

    Driven on the forming bond alone a Claisen prices at +74 kcal/mol; driven
    on both, it peaks at +29.0 where the literature says 30.
    """
    body = EDITOR_SOURCE.split("def on_submit_scan_run(")[1]
    body = body.split("def _scan_verdict(")[0]
    # Every leg gets the same fraction of its own span, so they move together.
    assert "share = n / float(steps)" in body
    assert "one['from'] + share * (one['to'] - one['from'])" in body
    # Held exactly: a scan point is a place, not a suggestion.
    assert "'mode': 'fix'" in body
    # And each point starts from the geometry the last one reached.
    assert "walked = outcome['xyz']" in body


def test_a_scan_point_is_relaxed_far_enough_to_mean_something():
    """A relaxed scan cut short reads as a barrier that is not there.

    Measured on a Diels-Alder approach whose relaxed cost is +3.6 kcal/mol:
    three cycles say +40.0, six +18.9, twenty +8.3.
    """
    line = next(one for one in EDITOR_SOURCE.splitlines()
                if one.strip().startswith("_SCAN_CYCLES ="))
    assert int(line.split("=")[1].strip()) >= 20


def test_the_scan_says_which_temperature_it_would_take():
    """Not only that it is closed."""
    body = EDITOR_SOURCE.split("def _scan_verdict(")[1].split("\n    def ")[0]
    assert "thermal_temperature(" in body
    assert "needs ~" in body
    # Wrapped across two lines in the source, so only the start of it.
    assert "K \\u2192 {stands}" in body


def test_a_scan_stops_at_the_next_minimum():
    """Past it, a scan stops describing a reaction and pushes into a structure.

    Measured on a tert-butyl bromide with a chloride, the carbon driven at the
    nucleophile past the product well: stopped at the minimum it walks 26 of 32
    points and reports the crossing at +25.2 kcal/mol, wanting 336 K.  Left to
    run it walks all 32, squeezes the new C-Cl to 1.20 A, and reports +202.9
    and 2566 K -- the barrier destroyed by the steps that came after it.

    Not detected by the path going flat: a scan keeps driving its coordinate,
    so past the well it climbs the far wall and never flattens.  That rule
    never fired once.
    """
    body = EDITOR_SOURCE.split("def _scan_arrived(")[1].split("def ")[0]
    # Over the top, down into something, and rising again.
    assert "_SCAN_OVER_THE_TOP" in body
    assert "_SCAN_CLIMBING" in body
    assert "b - a > _SCAN_UPHILL" in body
    # And it can be switched off for a deliberate overshoot.
    assert "submit_scan_whole" in EDITOR_SOURCE
    assert "not submit_scan_whole.value and _scan_arrived(path)" in EDITOR_SOURCE


def test_a_scan_comes_back_to_the_minimum_it_walked_through():
    """Stopping is not enough: the climb out is what recognises the minimum, so
    by then the walk stands two steps past it and the structure is squeezed
    that much.

    Measured on the same tert-butyl bromide: it now ends at C-Cl 1.87 A, which
    is the product bond, instead of two steps beyond it.
    """
    body = EDITOR_SOURCE.split("def on_submit_scan_run(")[1]
    body = body.split("def _scan_verdict(")[0]
    assert "summit, bottom = _descent(" in body
    assert "walked = bottom[1]" in body
    # And the geometry left in the box is that one.
    assert "def _done(final=walked):" in body
    assert "Scanned to the next minimum" in body


def test_the_temperature_is_named_whether_the_path_is_open_or_not():
    """The number a chemist came for was missing exactly when the news was
    good.

    "It needs 150 K and you have 298" is what makes an open path mean
    something; without it an open scan said only how long it took.
    """
    body = EDITOR_SOURCE.split("def _scan_verdict(")[1].split("\n    def ")[0]
    # One phrase, before the branch, so it is carried whichever way the
    # ceiling falls -- there is one line now rather than two nearly identical
    # ones, and "open" or "closed" is the whole of the difference.
    assert body.index("stands = 'open'") < body.index("if needs is None:")
    # And the word that says which it is, in the one clause both cases share.
    assert "stands = 'open' if rise <= ceiling else 'closed'" in body


def test_squeezing_is_refused_without_asking_the_energy():
    """A net that holds where the energy is weakest.

    Stretching is what a reaction does and the budget prices it; two atoms
    inside two thirds of a bond are not on any path at any temperature, and
    saying so needs no calculation -- which is why it also holds for a
    transition metal, where GFN2 is shaky.

    The floor sits between what was measured on paths that must stay open --
    a cyclohexane flipping came to 0.94 of a bond, an SN2 driven through its
    saddle and past its product 0.74, a [1,5]-hydrogen shift 0.78 -- and the
    one that must not: a ring carbon pushed across its own ring, 0.35.
    """
    assert 0.35 < gfn._TOO_CLOSE < 0.74

    ethane = (
        "8\nethane\n"
        "C  0.000  0.000  0.000\nC  1.530  0.000  0.000\n"
        "H -0.370  1.020  0.000\nH -0.370 -0.510  0.880\n"
        "H -0.370 -0.510 -0.880\nH  1.900  0.510  0.880\n"
        "H  1.900  0.510 -0.880\nH  1.900 -1.020  0.000\n"
    )
    share, i, j = gfn.closest_contact(ethane)
    assert share > gfn._TOO_CLOSE
    assert {i, j} == {0, 1} or share == pytest.approx(
        gfn.closest_contact(ethane)[0])

    squashed = ethane.replace("C  1.530  0.000  0.000", "C  0.700  0.000  0.000")
    share, _i, _j = gfn.closest_contact(squashed)
    assert share < gfn._TOO_CLOSE

    assert gfn.closest_contact("1\nlone\nC 0 0 0\n") == (None, None, None)


def test_a_loose_hold_stops_the_drag_rather_than_only_saying_so():
    """Advancing the mark would confirm a geometry nothing has priced.

    On a palladium pushed at head on the hold went 0.7 A loose and the drag
    carried on regardless, the picture showing a bromide 1.27 A from the metal
    while the price belonged to a structure where it had got out of the way.

    These two refuse at any temperature, so they are not the budget refusing
    and they are no longer said as though they were. Measured on a butane
    turned about its middle bond, the line read "+0.0 of 22.3 kcal/mol
    available at 298.15 K. Past the budget, so the last structure that was
    inside it is back." -- which sends a reader to the temperature box to fix
    something that is not there.
    """
    follow = EDITOR_SOURCE.split("def _thermal_wall(")[1].split("def ")[0]
    assert "def _thermal_wall" not in follow
    # On the highest price the drag has been at, not on the price of where it
    # is standing -- the ceiling is a barrier height.  The two are the same
    # number on the way up, which is the case this test is about.
    assert "if peak <= ceiling and not refuse:" in follow
    body = EDITOR_SOURCE.split("state['thermal_now'] = priced.get('energy')")[1]
    body = body.split("state['gfn_last_status'] = said")[0]
    assert "aside = (slipped > _SLIP_LOOSE) or crowded" in body
    assert "refuse=aside" in body
    assert "else current)[0]" in body
    # And the two are told apart where the refusal is reported.
    assert "'the last measured and allowed structure is '" in body
    assert "'past the budget -- the last structure inside '" in body


def test_the_line_says_how_steep_it_is():
    """It is already what shortens the leash, so a steep drag feels heavy.

    The number is what turns that into something to steer by: "+8 per A" and
    "+160 per A" are two different situations.
    """
    assert "state['thermal_slope'] = (float(spent)" in EDITOR_SOURCE
    assert "kcal/mol per A here." in EDITOR_SOURCE
    assert "'Climbing' if slope > 0 else 'Falling'" in EDITOR_SOURCE


def test_a_setting_is_not_taken_back():
    """The charge is read off the SMILES when the structure is built, once.

    Read afresh whenever a method is chosen, it put a hand-set charge back to
    zero -- and a charge that had come out of the SMILES correctly went to
    zero too, once the cached SMILES was no longer the one carrying it.  At
    zero those systems have an odd number of electrons and xtb refuses every
    step, so nothing runs, nothing moves, and the editor looks broken rather
    than misconfigured.
    """
    body = EDITOR_SOURCE.split("def _fill_charge_from_smiles(")[1]
    body = body.split("def on_submit_gfn_charge(")[0]
    # Once per structure ...
    assert "if state.get('charge_filled_for') == smiles:" in body
    assert "state['charge_filled_for'] = smiles" in body
    # ... and never over a number the user typed.
    assert "if state.get('charge_is_the_users'):" in body

    # The flag is raised by the user's own edits, not by the fill itself.
    said = EDITOR_SOURCE.split("def on_submit_gfn_charge(")[1].split("def ")[0]
    assert "state.get('charge_filling')" in said
    assert "state['charge_is_the_users'] = True" in said
    assert "submit_gfn_charge.observe(on_submit_gfn_charge" in EDITOR_SOURCE

    # A new structure from a new SMILES starts over: then the SMILES speaks.
    assert EDITOR_SOURCE.count("state['charge_is_the_users'] = False") >= 2


def test_a_scan_is_told_a_direction_rather_than_an_end():
    """Where it ends is the chemistry, not a number.

    A scan stops at the next minimum, so the end was the wrong thing to ask
    for -- and a number typed there is how two atoms came to be asked for
    0.60 A when the bond between them is 1.53.  The direction is the one thing
    that cannot be read off the selection, so that is what is asked; the
    number stays for the two cases that need it, following a figure from the
    literature and a coupled scan over two coordinates, where it is the ratio
    of the two ends that fixes the path.

    Driven on a real editor with nothing typed: further apart walks a C-C from
    1.52 to 4.02 A and wants 1607 K, closer walks it to 1.29 and wants 334 K.
    """
    assert "submit_scan_way" in EDITOR_SOURCE
    body = EDITOR_SOURCE.split("def _suggest_scan_target(")[1].split("def ")[0]
    # Outwards has two names now: the direction, and the verb "break" that
    # points the same way and carries a stopping rule with it.  Both walk to
    # the far end when nobody has said where to stop, which is what this is
    # about.
    assert "way in ('out', 'break')" in body
    assert "_SCAN_AS_FAR_AS" in body

    # The far end is only the brake for a walk with no next minimum -- two
    # fragments pulled apart never find one -- so it is generous.
    line = next(one for one in EDITOR_SOURCE.splitlines()
                if one.strip().startswith("_SCAN_AS_FAR_AS ="))
    assert "'distance': 2.5" in line

    # And it still gives way to the floor: nothing may be asked closer than
    # the bond it would make.
    arm = EDITOR_SOURCE.split("def on_submit_scan(")[1].split("def ")[0]
    assert "_suggest_scan_target(kind, here, submit_scan_way.value)" in arm
    assert "floor = _scan_floor_for(leg)" in arm
    assert arm.index("_suggest_scan_target(") < arm.index("floor = ")


def test_the_editor_writing_the_box_is_not_a_new_structure():
    """The host starts a structure over unless it is told it is the same one.

    Charge back to zero, multiplicity to one, and every held value, bond edit
    and history entry thrown away.  That is right for a structure someone
    loads and wrong for every write this editor makes -- and a scan makes one
    per point, so it lost the charge it was told to run at mid-scan, and the
    undo history with it.

    Measured with the host's own reset attached, a six-point scan on a charge
    of -1 and a multiplicity of 2: without the flag it came back at 0 and 1
    with no history at all, and the host had called it a new structure seven
    times; with it, -1 and 2, the history intact, and not once.
    """
    body = EDITOR_SOURCE.split("def _write_coords(")[1].split("\n    def ")[0]
    assert "state['structure_edit_inflight'] = True" in body
    assert "state['structure_edit_inflight'] = False" in body
    # Set around the write and cleared afterwards, whatever happens.
    assert body.index("= True") < body.index("coords_widget.value = text")
    assert "finally:" in body




def test_a_drag_is_not_a_new_structure_either():
    """The drag's own way back into the box was the one that killed the budget.

    Unguarded, the host started the structure over on every single drag --
    charge to zero, multiplicity to one, held values and history gone, and
    structure_changed dropping the thermal anchor with them.  Without an anchor
    there is no budget, so nothing could be refused and nothing taken back,
    which is why the budget appeared to do nothing at all.

    Measured with the host's reset attached, a ring hydrogen yanked to 1.90,
    2.60 and 3.20 A: before, the anchor was gone after the first yank and the
    viewer drew 1.96; after, the anchor is still set on all three, the host has
    called it a new structure not once, and both the viewer and the box hold
    1.08.
    """
    body = EDITOR_SOURCE.split("def on_submit_manip_sync(")[1]
    body = body.split("\n    def ")[0]
    assert "state['structure_edit_inflight'] = True" in body
    assert "coords_widget.value = payload" in body
    assert body.index("= True") < body.index("coords_widget.value = payload")
    assert "finally:" in body


def test_a_scan_streams_frames_instead_of_rewriting_the_box():
    """Every write to the box rebuilds the viewer from nothing.

    A thirty-point scan was thirty rebuilds: the picture crawled and the
    browser sometimes stopped answering altogether.  The follow has always
    streamed frames for exactly this reason.  Measured on a twenty-point scan:
    the box is written once, at the end, and twenty frames go down the
    channel.
    """
    body = EDITOR_SOURCE.split("def on_submit_scan_run(")[1]
    body = body.split("def _scan_verdict(")[0]
    assert "shown.append(_gfn.coordinates_of(walked))" in body
    assert "setattr(submit_gfn_frame, 'value', text)" in body
    # And the box only at the end, where the result is.
    assert "_write_coords" not in body.split("def _done(")[0]
    assert "_write_coords" in body


def test_the_rollback_is_drawn_down_the_channel_every_other_frame_uses():
    """After a crossing the budget refuses every step until the structure is
    back within 0.25 A of the last geometry it agreed to, so the picture goes
    between two geometries at the rate the answers arrive.  That is the wall
    working as designed.  What would not be working is a *second* writer
    drawing it -- a rollback that reached the viewer by a path of its own would
    be the thermal budget being a third mode rather than Dynamik Opt with a
    ceiling on it.

    There is one writer.  The follow writes the frame channel once, at the end
    of a step, whether the step was allowed or taken back: the geometry that
    survives is appended to the same trail and goes out under the same run
    number with the same pace on it.  The box is written once in the same
    place.  The wall has a field of its own and it carries marks for the
    browser's hand, never coordinates.

    Measured on a real editor at 298 K against a 22.3 kcal/mol ceiling with an
    energy costing 30 kcal/mol per angstrom of pull: four answers, four writes,
    all on run 1, all marked as a follow and all stamped at the slider's pace,
    the fourth carrying the third's geometry again -- and from there the follow
    stands still rather than shaking, which costs no xtb call and writes no
    frames at all.
    """
    follow = EDITOR_SOURCE.split("def _gfn_follow_step(", 1)[1].split(
        "_start_background(_work, 'The relaxation under the hand')", 1)[0]
    assert follow.count("submit_gfn_frame") == 1, (
        "the drag has a second frame writer in it")
    assert follow.count("_write_coords") == 1, (
        "one write to the box per step, allowed or refused")
    # Once it is clearly refusing it stands still: no xtb, no frames.
    assert "if _still_spent(current, holding):" in follow
    assert "continue" in follow.split(
        "if _still_spent(current, holding):", 1)[1].split("\n", 6)[5]

    wall = EDITOR_SOURCE.split("def _push_thermal_wall(", 1)[1].split(
        "\n    def ", 1)[0]
    assert "submit_gfn_wall" in wall
    assert "submit_gfn_frame" not in wall, "the wall must not draw geometries"
    assert "_write_coords" not in wall


def test_a_still_hand_holds_what_it_was_already_holding():
    """Nothing changed, so nothing about what is held should change either.

    Derived afresh, a hand that has stopped can be handed a *different* set --
    the nearest contact rather than the one the drag was driving -- and the
    structure then slides towards something nobody asked for while nothing is
    moving.
    """
    body = gfn.contacts_holding.__doc__ or ""
    assert "holding" in gfn.contacts_holding.__code__.co_varnames
    src = EDITOR_SOURCE.split("holding=state.get('thermal_holding')")
    assert len(src) > 1, "the follow has to hand the last set back in"
    assert "state['thermal_holding'] = [dict(one) for one in contacts]" \
        in EDITOR_SOURCE


# --- what an audit of the transitions turned up ----------------------------


def test_a_scan_honours_what_the_user_is_holding():
    """"Hold this while you walk that" is the whole point of a relaxed scan.

    The list was built from the scan's own legs and nothing else, so a held
    value never reached xtb.  Measured on an ethane: a C-H held at 1.60 came
    back at 1.080 after the scan and at 1.599 under Optimise -- the same hold,
    honoured by one and ignored by the other, with the list showing it
    throughout and the line saying nothing.
    """
    body = EDITOR_SOURCE.split("def on_submit_scan_run(")[1]
    body = body.split("def _scan_verdict(")[0]
    assert "walking = {tuple(one['atoms']) for one in legs}" in body
    assert "state.get('constraints')" in body
    # A coordinate that is both held and walked is walked: two values for one
    # thing cannot both be met, and the leg is what the user just asked for.
    assert "not in walking" in body


def test_a_scan_reads_its_start_from_the_structure_it_is_on():
    """Left as armed, a second press walked from a value the molecule no
    longer had.

    Measured: a C-C at 4.012 was compressed to 2.137 in one step and the run
    reported a 77 kcal/mol barrier with a temperature and a timescale
    attached, all of it invented.
    """
    body = EDITOR_SOURCE.split("def on_submit_scan_run(")[1]
    body = body.split("def _scan_verdict(")[0]
    assert "_value_of_constraint(one)" in body
    assert "dict(one, **{'from': now})" in body
    # And a leg naming atoms that are gone says so instead of throwing.
    assert "names atoms '" in body


def test_armed_legs_do_not_outlive_their_structure():
    """A scan armed on one molecule names atoms the next one may not have.

    Left armed it threw on the first click after loading a smaller structure,
    and run it walked a coordinate that was not there at all -- reported as a
    completed scan, because the run reads only whether xtb answered.
    """
    assert "state['scan_legs'] = []" in EDITOR_SOURCE
    from delfin.dashboard import tab_submit
    host = pathlib.Path(tab_submit.__file__).read_text(encoding="utf-8")
    assert "state['scan_legs'] = []" in host
    # And naming the atoms of a leg survives atoms that are gone.  One
    # function does it, because the profile puts the same pair on an axis,
    # the sentence names it, and the refusals about forming and breaking
    # name it too -- one place that names atoms, so one place that
    # survives their going.
    body = (EDITOR_SOURCE.split("def _leg_atoms_label(")[1]
            .split("\n    def ")[0])
    assert "if 0 <= index < len(known)" in body
    assert "_leg_atoms_label(leg)" in EDITOR_SOURCE.split(
        "def _describe_leg(")[1].split("\n    def ")[0]


def test_stop_comes_before_the_method_check():
    """Behind it, a scan started under GFN2 and then switched away could not be
    stopped at all: the button read Stop, was enabled, answered "a scan needs
    xtb", and the run it was meant to end walked to the last point and
    overwrote the box."""
    body = EDITOR_SOURCE.split("def on_submit_scan_run(")[1]
    body = body.split("def _scan_verdict(")[0]
    assert body.index("state.get('scan_run')") < body.index("A scan needs xtb")


def test_the_anchor_belongs_to_the_method_that_measured_it():
    """An energy of one method against energies of another is not a
    difference.

    Measured: an anchor taken under GFN2 and read against GFN-FF priced an
    untouched structure at +6384 kcal/mol against a 22.3 ceiling, so every
    drag sprang straight back.  The other way round is quieter and worse -- a
    C-C torn to 3 A reported "-6117.9 of 22.3 available" and the wall never
    fired.
    """
    body = EDITOR_SOURCE.split("def _thermal_budget(")[1].split("\n    def ")[0]
    assert "state.get('thermal_method') != str(submit_ff_dd.value)" in body
    assert "state['thermal_method'] = method" in EDITOR_SOURCE


def test_the_editor_does_not_forge_the_users_consent_on_the_charge():
    """reset_controls writes the charge, which fires the observer that
    remembers a number the user typed.

    So a reset the editor performed marked the charge as theirs, and from then
    on it was never read off a SMILES again -- measured, an acetate then ran at
    zero.
    """
    for name in ("def reset_controls(", "def _apply_controls("):
        body = EDITOR_SOURCE.split(name)[1].split("\n    def ")[0]
        assert "state['charge_filling'] = True" in body, name
        assert "finally:" in body, name


def test_the_scan_controls_belong_to_the_method_that_can_run_them():
    """Left visible under UFF or PM7 a whole scan could be armed, with the line
    saying "or press Run scan" -- an instruction that cannot work -- and the
    refusal arrived only on the press."""
    body = EDITOR_SOURCE.split("def _refresh_method_controls")[1]
    body = body.split("\n    def ")[0]
    assert "submit_scan_add_btn.layout.display = '' if xtb else 'none'" in body
    assert "submit_scan_run_btn" in body
    # The budget and the pull are decided elsewhere, because they answer to
    # the hand as well as to the method and one function has to hold both
    # rules or each refresh undoes the other.
    assert "_refresh_hand_controls()" in body
    assert "submit_thermal_btn.layout.display" not in body


def test_holding_a_value_is_a_step_of_its_own():
    """Judged on the picture alone a Hold was never one -- it changes no
    coordinate -- so Undo walked straight past it, wiped it on the way, and
    reported whatever it did land on.  Hold, then a scan, then Undo took back
    two actions on one press and named one of them."""
    # What an entry carries is gathered in one place now, and the holds are
    # part of it.
    marks = EDITOR_SOURCE.split("def _structure_marks(")[1].split("\n    def ")[0]
    assert "'constraints': [dict(one)" in marks
    assert "state.get('constraints') or []" in marks
    # And two entries are the same only when everything they carry is the
    # same -- not only the picture.  Compared on coordinates alone, a Hold
    # was never a step: it changes none.
    body = EDITOR_SOURCE.split("def _remember(what")[1].split("\n    def ")[0]
    assert "entry = dict(_structure_marks()," in body
    assert "last.get(key) == value for key, value in entry.items()" in body
    hold = EDITOR_SOURCE.split("def on_submit_hold(")[1].split("\n    def ")[0]
    assert "_remember(f'holding" in hold



def test_a_hold_that_names_nothing_is_not_priced():
    """The page names the atoms it is holding by their place in the structure
    it is looking at, and the kernel reads them against the structure it was
    sent.  When the two have moved apart -- the structure changed under the
    drag, or a frame arrived from before it did -- not one of those numbers
    names an atom that is there.

    :func:`gfn_optimize.contacts_holding` then correctly finds no coordinate
    to hold and returns nothing, which is indistinguishable from the honest
    empty case (the hand on everything, or on nothing) -- and the price fell
    straight through to the rigid single point.  That is the one that
    overcharges anything forming a bond twenty to thirty times over, and it
    did so with nothing said: the wall could take a drag back on a number
    computed for a structure neither side was looking at.
    """
    assert gfn.contacts_holding(_ETHANE, [99, 100]) == []
    source = EDITOR_SOURCE
    assert 'count = len(_gfn.atom_lines(current))' in source
    assert 'stale = bool(holding) and not any(' in source
    assert "pricing = _thermal_live() and not stale" in source
    # Nothing is charged, nothing is remembered, and the wall does not run.
    assert 'if not stale:\n' in source
    assert 'if pricing:' in source
    assert 'elif pricing:' in source
    assert 'else:\n                        priced = {}' in source
    # And it is said, because a step that is not priced looks exactly like a
    # step that cost nothing.
    assert 'does not have, so this step is not priced.' in source


def test_the_anchor_never_takes_the_drag_back():
    """Switching the budget on mid-drag asked for the drag to be measured, and
    answered by undoing it.

    The anchor relaxes first by default -- which is right when it is set on a
    structure at rest, and exactly wrong while a hand is on one: it optimised
    the geometry under the cursor and wrote the result into the coordinate
    box.  A hand on the molecule outranks the setting, and the structure the
    user is holding is measured as it stands.

    The same optimisation takes as long as it takes, and the editor is not
    frozen while it runs.  Its answer is only allowed to move the structure it
    was actually measuring; when the editor has moved on, the anchor belongs
    to the older structure, which is what the fingerprint already decides --
    so all that is left is to say so rather than overwrite what was done in
    the meantime.
    """
    source = EDITOR_SOURCE
    assert "held = bool(state.get('gfn_follow'))" in source
    assert 'if held:\n            wants_relax = False' in source
    assert "note = 'Measuring from the structure as it stands'" in source

    assert "moved_on = _current_xyz() != xyz" in source
    assert 'if wants_relax and not moved_on:' in source
    assert 'belongs to the older one' in source


def test_a_push_is_a_force_and_says_so_in_the_units_of_one():
    """xtb takes its distance force constants in Hartree per Bohr squared, so
    0.5 says nothing about how hard it pulls until it is converted -- it is
    1120 kcal/mol/A^2, half again as stiff as a C-C bond, which is why a
    "pull" drags an ethane's C-C from 1.52 to 2.46 A when asked for 2.50.

    A push is set in the units the ramp is thought in: the force it applies,
    in kcal/mol/A.  The restraint is held :data:`PUSH_REACH` ahead of wherever
    the coordinate now stands, so what the structure feels is k*reach and not
    a spring that pulls harder the further it has to go.
    """
    assert gfn.EH_PER_BOHR2_IN_KCAL == pytest.approx(2240.9, abs=0.5)
    # 20 kcal/mol/A held one angstrom ahead.  The two is xtb's restraint being
    # k*d^2 rather than the half of it a spring is usually written with --
    # measured, and left out it made every force here a claim about half of
    # what was really applied.
    assert gfn.push_constant(20.0) == pytest.approx(
        20.0 / (2 * gfn.EH_PER_BOHR2_IN_KCAL), rel=1e-12)
    # Half the reach, twice the stiffness, the same force.
    assert (gfn.push_constant(20.0, reach=0.5)
            == pytest.approx(2 * gfn.push_constant(20.0), rel=1e-9))
    # AFIR's own range, read through its collision parameter.
    assert 4.0 <= gfn.PUSH_FORCE_FROM < gfn.PUSH_FORCE_TO <= 400.0


def test_a_push_and_a_hold_cannot_share_one_force_constant():
    """xtb takes one force constant for a whole $constrain block -- a second
    ``force constant=`` line is read and ignored, measured -- and a push is
    three orders of magnitude softer than a hold.

    Run together the hold's stiffness wins and the push silently becomes an
    ordinary scan: the same picture, the same number, a different experiment.
    So the stiffer one wins, ``mixed`` says the two were asked for at once,
    and the scan refuses rather than running the wrong thing.
    """
    soft = gfn.push_constant(20.0)
    alone = gfn.constraint_input(
        [{'kind': 'distance', 'atoms': [0, 1], 'value': 2.0,
          'mode': 'push', 'force': soft}], atoms=8)
    assert alone['force'] == pytest.approx(soft)
    assert alone['mixed'] is False
    together = gfn.constraint_input(
        [{'kind': 'distance', 'atoms': [0, 1], 'value': 2.0,
          'mode': 'push', 'force': soft},
         {'kind': 'distance', 'atoms': [2, 3], 'value': 1.6, 'mode': 'fix'}],
        atoms=8)
    assert together['force'] == gfn.FIX_FORCE_CONSTANT
    assert together['mixed'] is True

    source = EDITOR_SOURCE
    assert "pushing = str(submit_scan_how.value) == 'push'" in source
    assert 'cannot share one force constant in ' in source
    # A force between two atoms drives a distance; an angle is what the two
    # bonds that make it already say.
    assert "A push is a force between two atoms, so it walks distances." in source


def test_a_push_is_priced_without_its_own_force_in_the_answer():
    """A push leaves a real restraint energy behind -- it is meant to, that is
    the force -- and xtb reports the biased total.

    Measured under GFN2: an ethane whose C-C was pulled towards 2.50 A comes
    back at 2.464 with a reported -7.201065 Eh, and a single point on that
    same geometry gives -7.203332.  The 1.42 kcal/mol between them is the
    restraint, and priced as it stands every point of a push would carry its
    own push in the barrier.
    """
    source = EDITOR_SOURCE
    assert 'def _unbiased(here, applied=()):' in source
    assert 'def _priced(got, applied):' in source
    assert "energy = (_priced(outcome, held) if pushing" in source
    # And the path records where the coordinate *got to*, because a push does
    # not dictate a value -- that is the whole point of it.  The geometry goes
    # down beside the price so that the finished profile can be handed to a
    # better method afterwards; for a push those are the only structures on
    # the way over the crossing.
    assert 'reached = (_value_in(walked, legs[0]) if pushing' in source
    assert 'path.append((reached, spent, walked))' in source


def test_a_push_ramps_geometrically_and_prices_what_it_falls_through():
    """The range is thirty to one, and where the reaction goes over is what
    matters: measured on butadiene and ethylene, equal steps of six put no
    point at all on the barrier.  Twenty geometric points from 4 to 120 put
    five between 9 and 20.

    Even then a crossing has a threshold -- a force either pays for it or does
    not -- so the last step falls through: 2.53 A at +4.4 kcal/mol to 1.54 A
    at -64.2, with the top never sampled.  The segment it fell through is then
    walked with the coordinate held, which is the one thing coordinate driving
    is good at now that both ends are known rather than guessed: that puts a
    point at 2.363 A and +6.3, which is the top.

    The push found the path and the walk priced it.  Neither could have done
    both -- which is why the editor does both.
    """
    source = EDITOR_SOURCE
    assert 'growth = (_gfn.PUSH_FORCE_TO / _gfn.PUSH_FORCE_FROM) ** (' in source
    assert 'force = force * growth if n > 1 else force' in source
    # Bisected from where the step started, so what is refined is the force
    # and not the path.
    assert 'force = (stood + force) / 2.0' in source
    assert 'while (tries < _PUSH_REFINE and outcome.get(\'ok\')' in source
    # And the crossing itself is priced, held exactly so no restraint energy
    # is left in the answer.
    assert 'def _across(here, was, now_at):' in source
    assert 'path.extend(_across(' in source
    assert "'mode': 'fix'," in source
    # The zero is the structure as it stands: a path that crosses on its first
    # step has to have something to have crossed from.
    assert 'base = _unbiased(walked)' in source
    # With the structure it was measured on, which is what lets the whole
    # profile be priced again later.
    assert "path.append((_value_in(walked, legs[0]), 0.0, walked))" in source


@_needs_xtb
def test_a_force_below_the_bond_settles_and_one_above_it_breaks():
    """Which is the whole reason a push is a force rather than a value.

    A restraint held :data:`PUSH_REACH` ahead of wherever the coordinate now
    stands applies at most ``k * reach``, and the structure keeps whatever of
    that it can answer.  Measured under GFN2 on an ethane at rest at 1.5212 A,
    the same push applied over and over:

        60 kcal/mol/A    1.6270, 1.6413, 1.6413, 1.6413, 1.6413, 1.6413
        120 kcal/mol/A   1.7613, 1.8815, 1.9736, 2.0686, 2.2524, 2.8766

    Sixty settles and stays there however long it is pushed -- the bond can
    answer it.  A hundred and twenty never settles: past the inflection of the
    curve the restoring force falls away, so once the push is stronger than
    the most the bond can pull back with, nothing balances it and the bond
    goes.

    That is the statement the ramp is built on.  Whether a deformation is
    possible stops being an accident of how far the mouse was dragged, or of
    what number was typed into a target box, and becomes a force the user set
    -- and the scan finds the reaction by raising it until something gives.
    """
    import math

    def bond(text):
        here = gfn.coordinates_of(text)
        return math.dist(here[0:3], here[3:6])

    start = gfn.optimize_with_gfn(_ETHANE, "gfn2", max_steps=200, timeout=300)
    assert start.get("ok"), start.get("status")

    def pushed(text, force):
        out = gfn.optimize_with_gfn(
            text, "gfn2", max_steps=60, timeout=300,
            constraints=[{"kind": "distance", "atoms": [0, 1],
                          "mode": "push", "force": gfn.push_constant(force),
                          "value": bond(text) + gfn.PUSH_REACH}])
        assert out.get("ok"), out.get("status")
        return out["xyz"]

    gentle, walk = start["xyz"], []
    for _ in range(4):
        gentle = pushed(gentle, 60.0)
        walk.append(bond(gentle))
    # It arrives somewhere and stays: the last three are the same structure.
    assert walk[-1] == pytest.approx(walk[-2], abs=2e-3), walk
    assert walk[-2] == pytest.approx(walk[-3], abs=2e-3), walk
    assert walk[-1] < bond(start["xyz"]) + 0.25, walk

    hard, tore = start["xyz"], []
    for _ in range(6):
        hard = pushed(hard, 120.0)
        tore.append(bond(hard))
    # And it never arrives: every step is further out than the last.
    assert tore == sorted(tore), tore
    assert tore[-1] > 2.5, tore

    # The reported energy carries the restraint, which is why a push is priced
    # by a single point on the geometry it reached.
    biased = gfn.optimize_with_gfn(
        start["xyz"], "gfn2", max_steps=60, timeout=300,
        constraints=[{"kind": "distance", "atoms": [0, 1], "mode": "push",
                      "force": gfn.push_constant(120.0),
                      "value": bond(start["xyz"]) + gfn.PUSH_REACH}])
    clean = gfn.optimize_with_gfn(biased["xyz"], "gfn2", timeout=300,
                                  optimise=False)
    assert biased["energy"] > clean["energy"], (biased["energy"],
                                                clean["energy"])


def test_the_scan_controls_say_what_they_are():
    """A bare 0, a bare 20 and a lone "closer" told nobody anything.

    The direction read as a setting with no subject and the two numbers as
    numbers with no units, so the row could only be used by someone who
    already knew what it did.  The words follow what is picked -- two atoms
    are closer together or further apart, three or four are narrower or wider
    -- and the numbers carry their own labels.

    And the end of the walk is gone unless it is asked for.  A scan walks to
    the next minimum, so where it ends is the chemistry; a field always on
    screen showing a zero read as something that had to be filled in, and a
    zero typed into a distance is how two atoms came to be asked for 0.60 A
    when the bond between them is 1.53.

    Asked for in the same box that asks the direction, because it is the same
    question: a value says which way the walk goes all by itself, and
    "further apart, to 2.40" was the same fact twice with nothing checking
    that its two halves agreed.  One box, three answers, and the field for the
    number under the third.
    """
    source = EDITOR_SOURCE
    assert ("options=[('closer together', 'in'), ('further apart', 'out'),"
            in source)
    assert "('to a value you give', 'to')]" in source
    assert "else [('narrower', 'in'), ('wider', 'out')," in source
    assert "kind = _CONSTRAINT_KINDS.get(picked)" in source
    assert "description='steps'" in source
    assert "description='to'" in source

    assert 'submit_scan_stop_at' not in source, (
        'the checkbox that only revealed the field is a control of its own '
        'again')
    assert "set_end = wanted == '' and str(submit_scan_way.value) == 'to'" in source
    assert "submit_scan_to.layout.display = '' if set_end else 'none'" in source
    # The direction is what is always used; a value only overrides it when one
    # was actually asked for.
    assert ('target = _suggest_scan_target(kind, here, submit_scan_way.value)'
            in source)
    assert "if str(submit_scan_way.value) == 'to':" in source
    # Opened on the value the coordinate has, not on a zero to be guessed at.
    assert 'submit_scan_to.value = float(submit_internal_value.value)' in source
    assert 'submit_scan_way.observe(on_submit_scan_way' in source


def test_a_drag_that_is_only_a_drag_stays_at_five_cycles():
    """Twenty cycles are the budget's, not the hold's.

    A price has to be for a properly relaxed path or the wall stands in the
    wrong place, which is what twenty buys.  Keyed on whether anything was
    being held -- and with the hand now a force there is always something --
    every step of every drag paid for an accuracy nothing was going to read:
    measured on a 102-atom complex, one xtb process is 0.06 s at one cycle,
    0.09 at five and 0.12 at ten, so a drag would have gone from about ten
    answers a second to two.

    Interactivity is the point of the whole mode, so it is keyed on the
    budget instead.
    """
    source = EDITOR_SOURCE
    assert 'cycles=(_THERMAL_FOLLOW_CYCLES if pricing' in source
    assert 'else _GFN_FOLLOW_CYCLES),' in source
    assert '_GFN_FOLLOW_CYCLES = 5' in source
    assert '_THERMAL_FOLLOW_CYCLES = 20' in source


@_needs_xtb
def test_the_restraints_own_energy_is_arithmetic_not_a_second_calculation():
    """xtb reports the total *including* what the restraint contributes, and a
    push is meant to leave a real one behind -- that residue is the force.

    Priced as it stands, every point of a push carries its own hand in the
    barrier.  The obvious answer is a second calculation on the same geometry
    with the constraints taken out, and that is a whole extra xtb process on
    every point of every scan and every step of every drag.

    It is not needed.  Measured under GFN2 by differencing the reported total
    against that single point: the ratio against ``0.5 * k * d^2`` is 2.00 for
    every force constant from 0.02 to 0.5, so the form is ``k * d^2``; and an
    angle's residue matches ``k * rad^2`` to four figures while ``k * deg^2``
    is out by three thousand.  Subtracted that way it agrees with the real
    thing to 0.0000 kcal/mol over a distance, an angle and a torsion at once.
    """
    import math

    def value_of(xyz, entry):
        here = gfn.coordinates_of(xyz)
        at = [(here[3 * i], here[3 * i + 1], here[3 * i + 2])
              for i in entry["atoms"]]
        if entry["kind"] == "distance":
            return math.dist(at[0], at[1])
        if entry["kind"] == "angle":
            one = [a - b for a, b in zip(at[0], at[1])]
            two = [a - b for a, b in zip(at[2], at[1])]
            na = math.sqrt(sum(v * v for v in one))
            nb = math.sqrt(sum(v * v for v in two))
            cosine = sum(a * b for a, b in zip(one, two)) / (na * nb)
            return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
        return gfn._dihedral(at, 0, 1, 2, 3)

    start = gfn.optimize_with_gfn(_ETHANE, "gfn2", max_steps=300, timeout=300)
    assert start.get("ok"), start.get("status")
    asked = [{"kind": "distance", "atoms": [0, 1], "mode": "push",
              "force": gfn.push_constant(80.0), "value": 2.5}]
    pushed = gfn.optimize_with_gfn(start["xyz"], "gfn2", max_steps=200,
                                   timeout=300, constraints=asked)
    assert pushed.get("ok"), pushed.get("status")
    clean = gfn.optimize_with_gfn(pushed["xyz"], "gfn2", timeout=300,
                                  optimise=False)

    bias = gfn.restraint_energy(pushed["xyz"], asked, value_of)
    # The residue is real and it is large enough to matter.
    assert bias * 627.5095 > 1.0, bias
    # And taking it off gives what the extra calculation gives.
    assert (pushed["energy"] - bias) == pytest.approx(clean["energy"],
                                                      abs=1e-6)

    # Nothing held is nothing to take off.
    assert gfn.restraint_energy(start["xyz"], [], value_of) == 0.0
    # And a coordinate that cannot be read is said so rather than guessed at,
    # so the caller falls back to asking for the calculation.
    assert gfn.restraint_energy(
        start["xyz"], asked, lambda _x, _e: None) is None


def test_the_price_is_taken_off_rather_than_calculated_again():
    """One xtb process per follow step, not two.

    Interactivity is the point of the mode: a step is one process at five
    cycles, measured at 0.09 s on a 102-atom complex, and asking for a second
    one to price it doubles that for a number the subtraction already has
    exactly.
    """
    source = EDITOR_SOURCE
    assert 'bias = _gfn.restraint_energy(' in source
    assert "dict(outcome, energy=float(outcome['energy']) - bias)" in source
    # And when it cannot be worked out, the calculation is asked for after all
    # rather than a wrong number being reported -- as it is when the step it
    # would be taken off has no energy, which is what a stopped one comes
    # back as.
    assert 'if bias is not None' in source
    assert "and outcome.get('energy') is not None else" in source

    # A held value and a pull cannot share one force constant, so the hold
    # stands and the hand goes back to placing -- said, not done quietly.
    assert 'held_too = bool(constraints) and pull is not None' in source
    assert 'force constant, so the hand places' in source


#: Butadiene and ethylene, 3.13 A apart, the pair of forming bonds C0-C10 and
#: C3-C11.  The canonical case for this file, and the one the whole thermal
#: budget was found by.
_DIELS_ALDER = """16
butadiene + ethylene
C    -1.513968    -0.052220    -0.000218
C    -0.727133     1.018345     0.000508
C     0.727133     1.018345    -0.000489
C     1.513968    -0.052221    -0.000515
H    -2.586669     0.037851    -0.000292
H    -1.123085    -1.055140     0.000253
H    -1.182283     2.001164     0.000336
H     1.182284     2.001164    -0.000174
H     2.586666     0.037852     0.000482
H     1.123087    -1.055140     0.000109
C    -0.658113    -0.352221     3.000000
C     0.658112    -0.352221     3.000000
H    -1.229821     0.561379     3.000000
H    -1.229821    -1.265820     3.000000
H     1.229820     0.561379     3.000000
H     1.229820    -1.265820     3.000000
"""


#: Anti-butane, relaxed under GFN2.  The C-C-C-C torsion driven out of it
#: climbs an eclipsed barrier and lands in gauche, which is about a kcal/mol
#: *above* anti -- so the whole path stays above where it started, which is
#: the shape the rollback got wrong.
_ANTI_BUTANE = """14
anti butane
C           -1.90036682238427       -0.40522890085261       -0.00000000065310
C           -0.62801675190148        0.43382104896482        0.00000000150501
C            0.62801675180062       -0.43382104895564        0.00000000084359
C            1.90036682257197        0.40522890071014       -0.00000000032059
H           -2.78101451548991        0.23301179233175       -0.00000000265027
H           -1.93775123301378       -1.04166818927924        0.88181350768924
H           -1.93775122992223       -1.04166818905944       -0.88181350801416
H           -0.61878667852378        1.07937568029498        0.88172390770606
H           -0.61878667758192        1.07937568065370       -0.88172390502591
H            0.61878667902179       -1.07937568022256        0.88172390743112
H            0.61878667712936       -1.07937568041878       -0.88172390533084
H            2.78101451539509       -0.23301179225031       -0.00000000153618
H            1.93775123262623        1.04166818899364        0.88181350716460
H            1.93775123027231        1.04166818908956       -0.88181350880857
"""


def _a_part(structure=_ETHANE):
    """One structure editor, built over a coordinate box of its own."""
    import tempfile

    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ("calc", "archive", "office"):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / "calc", archive_dir=room / "archive",
                           office_dir=room / "office")
    ctx.run_js = lambda _script: None
    return structure_editor.build(
        ctx, state={}, coords_widget=widgets.Textarea(value=structure),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None, get_smiles_charge=lambda *a, **k: None)


def _budgeted(structure, kelvin=298.15, hand="pull"):
    """One editor with a thermal budget on it, anchored on that structure.

    The anchor is an xtb run of its own and it is taken in a thread, so this
    waits for it: without one there is no budget at all and every assertion
    below would be about a ceiling that was never in force.
    """
    import time

    part = _a_part(structure)
    part.submit_ff_dd.value = "gfn2"
    part.submit_relax_btn.value = True
    part.submit_hand_dd.value = hand
    part.submit_temperature.value = kelvin
    part.submit_thermal_btn.value = True
    began = time.time()
    while part.state.get("thermal_e0") is None and time.time() - began < 300:
        time.sleep(0.05)
    assert part.state.get("thermal_e0") is not None, "the anchor never landed"
    return part


def _drag_message(xyz, note):
    """A geometry as the page sends one, with the reason on the comment line."""
    rows = gfn.atom_lines(xyz)
    return f"{len(rows)}\n{note}\n" + "\n".join(rows)


def _shifted(xyz, indices, dx):
    """*xyz* with those atoms shifted along x -- the mouse's own wish."""
    from delfin.dashboard.structure_editor import xyz_line

    here = gfn.coordinates_of(xyz)
    rows = []
    for i, line in enumerate(gfn.atom_lines(xyz)):
        at = [here[3 * i], here[3 * i + 1], here[3 * i + 2]]
        if i in indices:
            at[0] += dx
        rows.append(xyz_line(line.split()[0], *at))
    return f"{len(rows)}\nmoved\n" + "\n".join(rows)


def _quiet(state, seconds=300.0):
    """Wait until the follow has answered everything it was sent."""
    import time

    end = time.time() + seconds
    while time.time() < end:
        if (not state.get("gfn_follow_busy")
                and state.get("gfn_follow_xyz") is None):
            return
        time.sleep(0.02)
    raise AssertionError("the follow never came back")


def _turned_about_x(xyz, indices, degrees):
    """*xyz* with those atoms turned about the x axis.

    The chains here are laid along x, so this is the gesture of taking hold of
    one end and turning it -- the hand driving a torsion rather than a bond.
    """
    import math

    from delfin.dashboard.structure_editor import xyz_line

    here = gfn.coordinates_of(xyz)
    angle = math.radians(degrees)
    rows = []
    for i, line in enumerate(gfn.atom_lines(xyz)):
        x, y, z = here[3 * i], here[3 * i + 1], here[3 * i + 2]
        if i in indices:
            y, z = (y * math.cos(angle) - z * math.sin(angle),
                    y * math.sin(angle) + z * math.cos(angle))
        rows.append(xyz_line(line.split()[0], x, y, z))
    return f"{len(rows)}\nturned\n" + "\n".join(rows)


def _apart(xyz, i, j):
    import math

    here = gfn.coordinates_of(xyz)
    return math.dist(here[3 * i:3 * i + 3], here[3 * j:3 * j + 3])


def _costs(xyz, anchor):
    """What that geometry costs against the anchor, in kcal/mol."""
    from delfin.dashboard.structure_editor import _HARTREE_TO_KCAL

    got = gfn.optimize_with_gfn(xyz, "gfn2", timeout=300, optimise=False)
    assert got.get("energy") is not None, got.get("status")
    return (float(got["energy"]) - float(anchor)) * _HARTREE_TO_KCAL


def _dragged_apart(part, begin, far=2.0):
    """Pull the far methyl of an ethane out, one page message at a time.

    Absolute coordinates, the way the browser sends them: its model is where
    the cursor has taken the atom, not where the last answer put it.
    """
    part._begin_gfn_follow()
    part._arm_thermal_leash()
    methyl = {1, 5, 6, 7}
    step = far / 8.0
    for n in range(1, 9):
        part.submit_manip_sync.value = _drag_message(
            _shifted(begin, methyl, step * n), "DELFIN drag-follow held=1,5,6,7")
        _quiet(part.state)
    part.submit_manip_sync.value = _drag_message(
        _shifted(begin, methyl, far), "DELFIN drag-end")
    _quiet(part.state)
    return part.coords_widget.value


def _scanned(how, steps=20, seconds=600, structure=None, legs=None,
             energy='E'):
    """Run the editor's own Run scan on the Diels-Alder and hand back the lot.

    The real part, driven the way the button drives it -- this class of defect
    survives every reading of the source, because the source says exactly what
    it means to do.
    """
    import tempfile
    import time

    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ("calc", "archive", "office"):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / "calc", archive_dir=room / "archive",
                           office_dir=room / "office")
    ctx.run_js = lambda _script: None
    state = {}
    text = structure if structure is not None else _DIELS_ALDER
    box = widgets.Textarea(value=text)
    part = structure_editor.build(
        ctx, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None, get_smiles_charge=lambda *a, **k: None)
    state["current_xyz_for_copy"] = {"content": text}
    part.submit_ff_dd.value = "gfn2"
    part.submit_scan_how.value = how
    part.submit_scan_energy.value = energy
    state["scan_legs"] = legs or [
        {"kind": "distance", "atoms": [0, 10], "from": 3.134, "to": 0.7,
         "steps": steps, "structure": None},
        {"kind": "distance", "atoms": [3, 11], "from": 3.135, "to": 0.7,
         "steps": steps, "structure": None},
    ]
    began = time.time()
    part.on_submit_scan_run()
    while state.get("scan_run") and time.time() - began < seconds:
        time.sleep(0.05)
    assert not state.get("scan_run"), "the scan never finished"
    return part, state, box


def _forming(text):
    import math

    here = gfn.coordinates_of(text)
    return (math.dist(here[0:3], here[30:33]),
            math.dist(here[9:12], here[33:36]))


@_needs_xtb
def test_a_scan_hands_back_the_minimum_it_crossed_into_not_the_one_it_left():
    """"The scan walked 1 of 20 points. Highest +0.0 kcal/mol."

    It had walked all twenty, crossed a barrier and landed in a product -- and
    then handed back the structure it started from and described *that* as
    what it found.  A scan that works reads exactly like a scan that does
    nothing.

    The lowest point since the top was kept as the lowest point of the whole
    path.  A scan starts in a well and climbs out of it, so the lowest point
    of the whole path is almost always the first one; the walk then "came back
    to the minimum it walked through" by going back to the beginning, and the
    path was cut to that single point, which threw the barrier away and left a
    verdict describing the stump.

    Driven here on the real part.  Both ways of walking, because it was never
    about which: measured under GFN2, the value walked tops at +6.2 kcal/mol
    at 2.28 A and the force pushed at +6.3 at 2.36, both ending -64.2 in
    cyclohexene at about 1.55.
    """
    for how in ("hold", "push"):
        part, state, box = _scanned(how)
        said = part.mol_status.value

        assert state.get("scan_arrived"), (how, said)
        # It crossed something, and says how much.
        assert "rise " in said, (how, said)
        rise = _the_rise(said)
        assert 4.0 < rise < 9.0, (how, said)
        # And it is not the temperature of no barrier at all.
        wants = float(said.split("needs ~")[1].split()[0])
        assert 60 < wants < 120, (how, said)

        # What is left in the box is the product it walked into, not the
        # reactants it started from.
        close, other = _forming(box.value)
        assert close < 1.8 and other < 1.8, (how, close, other)
        assert "Scanned to the next minimum" in box.value

        # And the way back goes through the places the walk went through.
        #
        # This asserted one press, on the grounds that a scan is one action.
        # It still is -- one entry, recorded before it starts -- but a scan
        # also crosses two structures anyone would want again, and it works
        # them out anyway for the free energies: where it started and the
        # highest point it crossed.  Those are put into the history in the
        # order they were reached, so the presses walk back through them and
        # Redo walks forward again.  See
        # tests/test_a_scan_can_be_walked_back_through.py.
        #
        # What the one press was for is kept and is checked here: the walk
        # back reaches the structure from before the scan, and nothing in
        # between is a press that appears to do nothing.
        seen = [_forming(box.value)]
        for _ in range(3):
            part.on_submit_manip_undo()
            seen.append(_forming(box.value))
        assert len(set(seen)) == len(seen), (how, seen)
        back, also = seen[-1]
        assert back > 3.0 and also > 3.0, (how, seen)


def test_the_walk_comes_back_to_the_bottom_of_the_descent():
    """Not to the bottom of the path, which is usually where it started.

    A scan starts in a well and climbs out of it, so the lowest point of the
    whole path is almost always the first one.  Kept as "the minimum it walked
    through", every scan that crossed anything handed back the structure it
    began with -- and, because the path was then cut back to that one point,
    reported "The scan walked 1 of 20 points. Highest +0.0 kcal/mol.  It came
    back to the minimum it walked through."  Which is what it did.  It is also
    indistinguishable from a scan that did nothing.
    """
    descent = _a_part()._descent

    # Walking up: every point is a new highest, and the descent restarts at
    # each of them, so the bottom is always the point itself.
    summit, bottom = None, None
    for spent, where in ((0.0, 'a'), (3.0, 'b'), (8.0, 'c')):
        summit, bottom = descent(summit, bottom, spent, where, spent)
    # The top keeps its geometry too, because a free energy is asked for
    # there and the walk has moved on by the time it is.
    assert (summit, bottom) == ((8.0, 'c', 8.0), (8.0, 'c', 8.0))

    # Coming down again: the bottom follows the descent...
    summit, bottom = descent(summit, bottom, 4.0, 'd', 4.0)
    assert bottom == (4.0, 'd', 4.0)
    # ...and *not* back to the start, which is lower than any of it.
    summit, bottom = descent(summit, bottom, 5.0, 'e', 5.0)
    assert bottom == (4.0, 'd', 4.0), 'the descent is what is being tracked'
    assert summit == (8.0, 'c', 8.0)

    # A second, higher top starts the descent over.
    summit, bottom = descent(summit, bottom, 9.0, 'f', 9.0)
    assert (summit, bottom) == ((9.0, 'f', 9.0), (9.0, 'f', 9.0))

    source = EDITOR_SOURCE
    assert 'summit, bottom = _descent(' in source
    # The path is not cut back to it either: cut, the verdict describes the
    # stump instead of the walk, and the barrier that was just crossed is
    # thrown away.
    assert "state['scan_came_back'] = (bottom[2], bottom[0])" in source
    assert "came = state.get('scan_came_back')" in source
    assert 'ends = came[1] if came else path[-1][1]' in source
    assert "state['scan_came_back'] = None" in source


@_needs_xtb
def test_a_walk_that_never_settles_says_so_rather_than_pretending():
    """Anti-butane, driven a whole turn of the C-C-C-C torsion.

    Measured under GFN2, ten degrees a step: it climbs the eclipsed barrier to
    +2.40 kcal/mol at 120, falls into gauche at +0.44 near 70, and climbs the
    syn barrier to +4.78 at 0.  Never more than 1.96 below the highest point
    it has reached, and the rule wants two -- so the walk goes the whole way
    and does not claim to have arrived anywhere.

    Which is the honest answer.  A gauche well 1.96 deep is a real minimum and
    a real thing to want to stop at; what the scan may not do is *say* it
    stopped at one when it did not, and the two are distinguished here so that
    a change to the threshold cannot quietly turn one into the other.
    """
    part, state, box = _scanned(
        "hold", structure=_ANTI_BUTANE,
        legs=[{"kind": "dihedral", "atoms": [0, 1, 2, 3], "from": 180.0,
               "to": -60.0, "steps": 24, "structure": None}])
    said = part.mol_status.value
    assert not state.get("scan_arrived"), said
    assert "came back to the minimum" not in said, said
    # The steps that were asked for, and the narrowing counted apart from
    # them: the walk locates its own summit afterwards, and folding those
    # points into one total took the two numbers the user set out of the
    # answer.
    assert "24 of 24 points" in said, said
    assert "at the top" in said, said
    # And the barrier it reports is the syn one it went over.
    rise = _the_rise(said)
    assert 3.5 < rise < 6.5, said


@_needs_xtb
def test_a_torsion_restraint_is_periodic_and_not_a_spring():
    """Which is why the hand could not turn anything.

    A harmonic in an angle has a step in it where the angle wraps, so xtb does
    not use one for a torsion: it uses ``k * (1 - cos d)``.  Measured under
    GFN2 on anti-butane held 120 degrees from its target, the residue is
    4.7063 kcal/mol -- ``k * rad^2`` says 13.76 and ``k * (1 - cos)`` says
    4.7063.  An angle *is* a spring in radians, and a distance one in Bohr,
    which makes three shapes for the one force constant xtb takes per block.

    So the same number pulls a fifth as hard on a torsion as on a distance.
    Sized as a distance, the hand's largest torque was 6.2 kcal/mol per radian
    against a rotational barrier needing about 11 -- it stuck 21 degrees out
    and converged there at any number of cycles, which is a force balance and
    not a calculation that ran out of steps.
    """
    import math

    def torsion(text):
        here = gfn.coordinates_of(text)
        at = [(here[3 * i], here[3 * i + 1], here[3 * i + 2])
              for i in range(4)]
        return gfn._dihedral(at, 0, 1, 2, 3)

    start = gfn.optimize_with_gfn(_ANTI_BUTANE, "gfn2", max_steps=400,
                                  timeout=300)
    assert start.get("ok"), start.get("status")
    asked = [{"kind": "dihedral", "atoms": [0, 1, 2, 3], "mode": "push",
              "force": 0.005, "value": 60.0}]
    held = gfn.optimize_with_gfn(start["xyz"], "gfn2", max_steps=300,
                                 timeout=300, constraints=asked)
    assert held.get("ok"), held.get("status")
    clean = gfn.optimize_with_gfn(held["xyz"], "gfn2", timeout=300,
                                  optimise=False)
    residue = (held["energy"] - clean["energy"]) * 627.5095
    lag = 60.0 - torsion(held["xyz"])

    # It did not move at all: 0.005 is far too soft to turn a butane.
    assert abs(lag) > 100.0, lag
    # And what it cost is the periodic form, not the spring.
    periodic = 0.005 * (1.0 - math.cos(math.radians(lag))) * 627.5095
    spring = 0.005 * math.radians(lag) ** 2 * 627.5095
    assert residue == pytest.approx(periodic, rel=0.02), (residue, periodic)
    assert residue < spring / 2, (residue, spring)

    # Which is what restraint_energy has to subtract, or a pushed torsion is
    # priced with the hand in it.
    def value_of(xyz, entry):
        here = gfn.coordinates_of(xyz)
        at = [(here[3 * i], here[3 * i + 1], here[3 * i + 2])
              for i in entry["atoms"]]
        return gfn._dihedral(at, 0, 1, 2, 3)

    bias = gfn.restraint_energy(held["xyz"], asked, value_of)
    assert (held["energy"] - bias) == pytest.approx(clean["energy"], abs=1e-5)


def test_the_hardest_the_hand_can_pull_is_the_number_it_is_set_to():
    """In each coordinate's own units, which is the only way one number can
    mean one thing across three shapes.

    A distance stores ``k * d^2`` in Bohr, an angle the same in radians, and a
    torsion ``k * (1 - cos d)``.  So a force constant sized for one is the
    wrong strength for the others, and xtb takes one for the whole block --
    which is why the block's constant is sized for the coordinate the hand is
    actually driving, and a turn is a torsion.
    """
    import math

    hand = 44.6                       # kcal/mol/A at 298 K within the hour

    # A distance: hardest at full stretch, 2 k reach.
    k = gfn.push_constant(hand)
    hardest = 2 * k * (gfn.PUSH_REACH / gfn.BOHR_IN_ANGSTROM) * 627.5095 \
        / gfn.BOHR_IN_ANGSTROM
    assert hardest == pytest.approx(hand, rel=1e-9)

    # A torsion: k sin d, hardest at a right angle, so simply k.
    kd = gfn.push_constant(hand, kind="dihedral")
    assert kd * 627.5095 == pytest.approx(hand, rel=1e-9)
    # And it takes a much larger number to pull as hard.
    assert kd > 5 * k

    # The block follows the coordinate the hand is driving.
    turning = gfn.as_pushes(
        [{"kind": "dihedral", "atoms": [0, 1, 2, 3], "value": 60.0},
         {"kind": "distance", "atoms": [0, 3], "value": 3.0}],
        None, hand)
    assert all(one["force"] == pytest.approx(kd) for one in turning), turning
    stretching = gfn.as_pushes(
        [{"kind": "distance", "atoms": [0, 3], "value": 3.0}], None, hand)
    assert stretching[0]["force"] == pytest.approx(k)

    # The reach for an angle is the same reach, in the units xtb measures one
    # in: an angstrom is 1.8897 Bohr, so it is 1.8897 radians.
    assert gfn.PUSH_REACH_DEGREES == pytest.approx(
        math.degrees(gfn.PUSH_REACH / gfn.BOHR_IN_ANGSTROM))
    assert gfn.PUSH_REACH_DEGREES == pytest.approx(108.3, abs=0.1)


@_needs_xtb
def test_room_temperature_can_turn_a_molecule_into_its_own_conformers():
    """The thing the budget must never forbid.

    Anti-butane, the C-C-C-C torsion pushed by a hand sized for 298 K within
    the hour.  Measured: it turns 70 degrees in the first step and goes on
    turning, where a hand sized as a distance -- which is what it was --
    stopped 21 degrees out and stayed there however many cycles it was given.

    At 150 K it still turns, because a 4.8 kcal/mol rotational barrier needs
    about 11 kcal/mol per radian and the hand there is 22.  It is the wall
    that refuses what the temperature cannot pay for, not the hand.
    """
    def torsion(text):
        here = gfn.coordinates_of(text)
        at = [(here[3 * i], here[3 * i + 1], here[3 * i + 2])
              for i in range(4)]
        return gfn._dihedral(at, 0, 1, 2, 3)

    start = gfn.optimize_with_gfn(_ANTI_BUTANE, "gfn2", max_steps=400,
                                  timeout=300)
    assert start.get("ok"), start.get("status")
    was = torsion(start["xyz"])

    def turned(ceiling):
        force = gfn.push_force_for(ceiling)
        k = gfn.push_constant(force, kind="dihedral")
        out = gfn.optimize_with_gfn(
            start["xyz"], "gfn2", max_steps=60, timeout=300,
            constraints=[{"kind": "dihedral", "atoms": [0, 1, 2, 3],
                          "mode": "push", "force": k,
                          "value": was - gfn.PUSH_REACH_DEGREES}])
        assert out.get("ok"), out.get("status")
        return abs(was - torsion(out["xyz"]))

    # A hand for room temperature turns it a long way in a single step.
    assert turned(22.3) > 45.0, turned(22.3)
    # And one for 150 K still turns it past the eclipsed barrier at 120.
    assert turned(10.9) > 45.0, turned(10.9)
    # Sized as a distance -- 22.3 / reach, and then read as a torsion
    # constant -- it sticks, which is what this was.
    stuck = gfn.optimize_with_gfn(
        start["xyz"], "gfn2", max_steps=60, timeout=300,
        constraints=[{"kind": "dihedral", "atoms": [0, 1, 2, 3],
                      "mode": "push", "force": gfn.push_constant(22.3),
                      "value": was - gfn.PUSH_REACH_DEGREES}])
    assert abs(was - torsion(stuck["xyz"])) < 30.0, "it used to stick here"


#: Bromobenzene, relaxed under GFN2.  An aryl C-Br is the case that has to
#: work both ways: room temperature genuinely cannot break one within the
#: hour, and a hand with no temperature to answer to has to be able to.
_BROMOBENZENE = """12
bromobenzene
C   0.000  1.396  0.000
C   1.209  0.698  0.000
C   1.209 -0.698  0.000
C   0.000 -1.396  0.000
C  -1.209 -0.698  0.000
C  -1.209  0.698  0.000
Br  0.000  3.300  0.000
H   2.155  1.244  0.000
H   2.155 -1.244  0.000
H   0.000 -2.486  0.000
H  -2.155 -1.244  0.000
H  -2.155  1.244  0.000
"""


def test_without_a_temperature_the_hand_has_no_ceiling():
    """Drag further, pull harder.  Which is what pulling on something is like.

    The reach was applied always, so the hand could never do more than the
    slider was already set for however far it was dragged -- a ceiling on a
    setting that has no reason for one.  A ceiling is what a *temperature* is:
    with the budget on the reach comes back, and past the point where the hand
    is already pulling as hard as the temperature allows, dragging further
    buys nothing.
    """
    part = _a_part()
    part.submit_thermal_btn.value = False
    assert part._pull_most() is None, 'nothing may cap the hand'
    part.submit_thermal_btn.value = True
    part.submit_temperature.value = 298.15
    # Nothing caps the hand at any temperature: a temperature grants an
    # energy, and the wall is what enforces that on what was reached.
    assert part._pull_most() is None

    # What is *asked for* stays within the reach either way -- a target that
    # runs arbitrarily far ahead is overshot in a few cycles and comes back,
    # which on screen is a molecule that shakes.
    far = gfn.as_pushes(
        [{'kind': 'distance', 'atoms': [0, 6], 'value': 6.0}],
        'reference', 44.0, value_of=lambda _x, _e: 1.9)
    assert far[0]['value'] == pytest.approx(1.9 + gfn.PUSH_REACH)
    # The strength carries the excess instead: four reaches out, four times
    # as hard, which is exactly what an unclamped spring would apply there.
    near = gfn.as_pushes(
        [{'kind': 'distance', 'atoms': [0, 6], 'value': 2.9}],
        'reference', 44.0, value_of=lambda _x, _e: 1.9)
    assert far[0]['force'] == pytest.approx(4.1 * near[0]['force'], rel=0.02)
    # Unless something says otherwise, and only a temperature does.
    held = gfn.as_pushes(
        [{'kind': 'distance', 'atoms': [0, 6], 'value': 6.0}],
        'reference', 44.0, value_of=lambda _x, _e: 1.9, most=44.0)
    assert held[0]['force'] == pytest.approx(near[0]['force'])

    source = EDITOR_SOURCE
    assert 'most=_pull_most()' in source
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js
    js = submit_manip_bootstrap_js()
    assert 'function setPullStrength(scopeKey, share, most)' in js
    assert 'function pullOver(scopeKey, viewer)' in js
    assert 'var k = pullConstant(state) * pullOver(scopeKey, viewer);' in js
    assert 'if (state.pullMost > 0) {' in js


@_needs_xtb
def test_room_temperature_holds_an_aryl_bromide_and_nothing_else_does():
    """The two halves of the same setting, on the same bond.

    Measured under GFN2, the C-Br of a bromobenzene at rest at 1.909 A with
    the mouse dragged five angstroms out:

        thermal budget on at 298 K   held at 2.120
        budget off, hand at 0.4      torn to 4.367
        budget off, hand at 3.0      torn to 4.403

    An aryl C-Br really does not break at room temperature within the hour,
    and a hand with no temperature to answer to really should be able to break
    one -- so both of those are right, and they are the same code with one
    number changed.
    """
    import math

    def bond(text):
        here = gfn.coordinates_of(text)
        return math.dist(here[0:3], here[18:21])

    start = gfn.optimize_with_gfn(_BROMOBENZENE, "gfn2", max_steps=400,
                                  timeout=300)
    assert start.get("ok"), start.get("status")
    rest = bond(start["xyz"])

    def dragged(force, reach):
        here, wish = start["xyz"], rest
        for _ in range(8):
            wish += 0.5                       # the mouse keeps going
            target = wish if not reach else min(wish, bond(here) + reach)
            out = gfn.optimize_with_gfn(
                here, "gfn2", max_steps=20, timeout=300,
                constraints=[{"kind": "distance", "atoms": [0, 6],
                              "mode": "push",
                              "force": gfn.push_constant(force),
                              "value": target}])
            if not out.get("ok"):
                break
            here = out["xyz"]
        return bond(here)

    # The temperature holds it, however far the mouse goes.
    held = dragged(gfn.push_force_for(22.3), gfn.PUSH_REACH)
    assert held < rest + 0.6, held

    # And with no temperature to answer to, the hand takes it off -- at the
    # setting the editor opens on, not only at the top of the slider.
    torn = dragged(0.4 * gfn.A_BOND_HOLDS, 0.0)
    assert torn > 3.5, torn
    assert dragged(3.0 * gfn.A_BOND_HOLDS, 0.0) > 3.5


@_needs_xtb
def test_keeping_the_bonds_keeps_the_molecule():
    """A third thing the drag can be asked to respect.

    The budget says what a temperature can pay for; this says the bonds are
    the ones they were, whatever it costs -- which is what "move this group
    over there" means when the answer is not supposed to be a reaction.

    xtb cannot be told to keep a topology, so it is a wall like the budget's:
    the step runs, the bonding is read off what came back, and a step that
    made or broke one is replaced by the last that did not.  Measured on the
    real part, a bromide dragged off a tertiary carbon under GFN2 with the
    hand at its hardest: without it the C-Br goes from 2.055 to 7.204 A and
    the structure comes back with 21 bonds instead of 22; with it the bond
    holds at 2.066 and all 22 are there.
    """
    graph = gfn.bond_graph(_BROMOBENZENE)
    assert (0, 6) in graph, "the C-Br is a bond to begin with"
    assert len(graph) == 12, graph

    # Pulled apart, it is not.
    rows = [line.split() for line in gfn.atom_lines(_BROMOBENZENE)]
    rows[6][2] = f"{float(rows[6][2]) + 3.0:.6f}"
    apart = "12\ntorn\n" + "\n".join(" ".join(r) for r in rows) + "\n"
    after = gfn.bond_graph(apart)
    assert (0, 6) not in after

    # And what changed is said in a sentence a chemist can act on.
    symbols = [str(r[0]) for r in rows]
    assert gfn.graph_changed(graph, after, symbols) == "breaks C1-Br7"
    assert gfn.graph_changed(after, graph, symbols) == "makes C1-Br7"
    assert gfn.graph_changed(graph, graph, symbols) == ""


def test_the_bonds_are_kept_by_taking_the_step_back():
    """Not by asking xtb to hold them, which it cannot do."""
    source = EDITOR_SOURCE
    assert 'def _topology_wall(xyz):' in source
    assert 'submit_topology_btn = widgets.ToggleButton(' in source
    assert "description='Keep bonds'" in source
    # Judged on what the user would be left with, not on the raw answer.
    assert 'kept, changed = _topology_wall(settled)' in source
    assert 'bonding is being kept.' in source
    # The box as well as the picture, or letting go keeps the broken one.
    # One write at the end of the step now, whichever of the three things
    # happened to it, so that the box is never left holding a geometry no
    # answer produced.
    assert "why = 'Kept: the bonding would have changed'" in source
    # It belongs to this molecule and no other.
    assert "state.get('topology_for') != who" in source
    # And where it already holds, it says so instead of looking busy.
    assert 'already keeps its bonding, ' in source


def test_the_whole_profile_is_a_switch_that_says_what_it_does():
    """"all the way" said nothing about what it was all the way to.

    What it is for is a whole curve: a torsion turned right round has three
    minima and three barriers, and stopping at the first is exactly wrong when
    the profile is the point.  Off stays the default -- past the next minimum
    a reaction scan is pushing into a structure rather than following one.
    """
    source = EDITOR_SOURCE
    assert 'submit_scan_whole = widgets.ToggleButton(' in source
    assert "description='Whole profile'" in source
    assert 'submit_scan_whole.observe(on_submit_scan_whole' in source
    assert 'not submit_scan_whole.value and _scan_arrived(path)' in source


def test_the_topology_watch_is_the_one_a_conformer_search_uses():
    """A bond is anything inside 1.3 times the sum of the covalent radii.

    ORCA's GOAT uses that number and, by default, keeps the bonds a structure
    came with -- only its EXPLORE variant lets them break, and its REACT
    variant counts the topological difference as the sum of broken plus formed
    bonds.  The same question, and a conformer search is where it has been
    thought about hardest, so the same number.

    What a conformer search does not need and a drag does is hysteresis.  It
    compares finished optimisations; a drag asks ten times a second, so a bond
    resting on the threshold decides the answer differently from one answer to
    the next and the wall fires on a molecule that is not changing at all.  A
    bond that was there has to be clearly gone before it counts as broken and
    one that was not has to be clearly there: measured, a bond length moves
    two to five hundredths of an angstrom between answers, and for a C-Br the
    band is 2.548 to 2.744 A.
    """
    import math

    assert gfn.BOND_STARTS_AT == 1.3
    assert not hasattr(gfn, 'BOND_STOPS_AT'), 'one threshold, like GOAT'

    was = gfn.bond_graph(_BROMOBENZENE)
    rows = [line.split() for line in gfn.atom_lines(_BROMOBENZENE)]

    def pulled(out):
        moved = [list(one) for one in rows]
        moved[6][2] = f"{float(moved[6][2]) + out:.6f}"
        return "12\npulled\n" + "\n".join(" ".join(one) for one in moved) + "\n"

    def span(text):
        here = gfn.coordinates_of(text)
        return math.dist(here[0:3], here[18:21])

    # Inside the threshold, and held.
    inside = pulled(0.4)
    assert span(inside) < gfn.BOND_STARTS_AT * 1.96, span(inside)
    assert gfn.graph_holds(was, inside)[0], span(inside)

    # Past it, and gone -- said by name, and by the same measure that accepted
    # the one above.  Two measures and one wall is not a wall: accepting at
    # 1.4 and reporting at 1.3 let a bond stretch across several accepted
    # steps, and what came back was a geometry already broken by the stricter
    # one.  Measured on a 2,4-hexadiene with the rigid hand, five steps:
    # three bonds broken without the wall, and one broken *with* it.
    torn = pulled(1.4)
    holds, said = gfn.graph_holds(was, torn)
    assert not holds, (span(torn), said)
    assert said == "breaks C1-Br7", said

    # A bond that was never there has to be clearly there to count.
    apart = gfn.bond_graph(pulled(3.0))
    assert gfn.graph_holds(apart, _BROMOBENZENE)[1] == "makes C1-Br7"

    # And nothing at all is a molecule that did not change.
    assert gfn.graph_holds(was, _BROMOBENZENE) == (True, '')


@_needs_xtb
def test_a_free_energy_is_a_hessian_and_is_asked_for_never_assumed():
    """The ceiling is a free energy of activation and what is priced against
    it is an electronic one, which is an approximation and is meant as one.

    A Hessian is not free: measured under GFN2, 0.57 s against 0.29 for
    sixteen atoms and 3.72 against 0.76 for twenty-four.  A drag answers ten
    times a second, so there is no version of that which has a free energy in
    it -- and an entropy term roughly constant along a path largely cancels in
    a difference, which is what makes the approximation a reasonable one
    rather than merely a cheap one.

    A scan can afford it, and that is where it is offered.
    """
    warm = gfn.optimize_with_gfn(_ETHANE, "gfn2", timeout=300, optimise=False,
                                 free_energy=True, thermo_kelvin=298.15)
    assert warm.get("ok"), warm.get("status")
    assert warm.get("free_energy") is not None
    # A free energy is above the electronic one at room temperature -- the
    # zero point alone sees to that.
    assert warm["free_energy"] > warm["energy"]

    hot = gfn.optimize_with_gfn(_ETHANE, "gfn2", timeout=300, optimise=False,
                                free_energy=True, thermo_kelvin=800.0)
    assert hot.get("ok"), hot.get("status")
    # The same structure, so the same electronic energy; and G falls as the
    # temperature rises, because -TS does.
    assert hot["energy"] == pytest.approx(warm["energy"], abs=1e-8)
    assert hot["free_energy"] < warm["free_energy"]

    # Not asked for, not taken.
    plain = gfn.optimize_with_gfn(_ETHANE, "gfn2", timeout=300, optimise=False)
    assert plain.get("free_energy") is None


@_needs_xtb
def test_the_scan_can_be_priced_with_free_energies():
    """At the three places they are both affordable and meaningful: where the
    walk started, the highest point it crossed, and the minimum it came to.

    Not at every point -- twenty Hessians is a scan that takes minutes instead
    of seconds, and an RRHO free energy only means something at a stationary
    point anyway.

    Measured on the Diels-Alder, pushed, at 298.15 K: +6.3 kcal/mol to the top
    and -64.2 to the end as electronic energies, +3.3 and -58.5 as free ones.
    Both differences go the way they should -- the approach complex is held
    together by nothing much, so it pays entropy to reach the top, and two
    molecules becoming one pay more of it to stay there.
    """
    part, state, box = _scanned("push")
    assert state.get("scan_free") is None, "E unless G is asked for"
    plain = part.mol_status.value
    assert "free energ" not in plain, plain

    part, state, box = _scanned("push", energy="G")
    free = state.get("scan_free")
    assert free is not None, part.mol_status.value
    top, ends = free
    assert 1.0 < top < 6.0, free
    assert -70.0 < ends < -45.0, free

    said = part.mol_status.value
    assert "G, " in said and "E kcal/mol at 298.15 K" in said, said
    # And the temperature is worked out from the free energy, not from the
    # electronic one standing beside it: read the other way round, the line
    # quoted one number and reasoned from the other.
    wants = float(said.split("needs ~")[1].split()[0])
    assert wants == pytest.approx(
        thermal_temperature(top, 3600.0), abs=5.0), (wants, top)


def test_the_free_energy_is_taken_where_it_means_something():
    """Three Hessians, on geometries the walk keeps for the purpose."""
    source = EDITOR_SOURCE
    assert 'submit_scan_energy = widgets.Dropdown(' in source
    assert "('price with G', 'G')" in source
    assert 'def _free(here):' in source
    assert 'free_energy=True,' in source
    # Unconstrained, or the Hessian has the restraint's own curvature in it.
    assert 'with no restraint in it' in source
    # The start is the first point of the walk, which is the structure relaxed
    # at the value it already had.
    assert 'if began_at is None:' in source
    assert 'ends = bottom[1] if bottom is not None else walked' in source
    # And it is the barrier the temperature is worked out from.
    assert "free = state.get('scan_free')" in source
    assert 'rise = free[0]' in source


@_needs_xtb
def test_xtb_finds_its_own_way_between_the_two_ends_of_a_scan():
    """A scan drives a coordinate somebody chose; the path finder is given two
    structures and finds its own way between them.

    So it answers the question the scan can only approach -- and it needs the
    scan first, because the scan is what makes a product to aim at.  Measured
    on the butadiene and ethylene: 23 points in under four seconds, a forward
    barrier of 3 to 6 kcal/mol, 69 back, about -65 for the reaction, and an
    estimated transition state with the two forming bonds at 2.52 A.  The scan
    of the same reaction put its highest point at +6.3 and 2.36.  Two methods
    that share no machinery, agreeing.

    It is the same twice, which was written here the other way round and was
    wrong: three goes at one pair gave 5.755 kcal/mol to the digit, and four
    at another gave 43.7.  xtb seeds it fixed.  What moves the number is which
    two structures it is given -- the same reaction from a slightly different
    reactant complex gave 5.755 and 3.3 -- so a barrier from this is a
    statement about the two ends as much as about the reaction.
    """
    import math

    start = gfn.optimize_with_gfn(_DIELS_ALDER, "gfn2", max_steps=400,
                                  timeout=300)
    assert start.get("ok"), start.get("status")
    walked = start["xyz"]
    for target in (2.8, 2.4, 2.0, 1.7, 1.56):
        got = gfn.optimize_with_gfn(
            walked, "gfn2", max_steps=80, timeout=300,
            constraints=[{"kind": "distance", "atoms": [0, 10],
                          "mode": "fix", "value": target},
                         {"kind": "distance", "atoms": [3, 11],
                          "mode": "fix", "value": target}])
        assert got.get("ok"), got.get("status")
        walked = got["xyz"]
    product = gfn.optimize_with_gfn(walked, "gfn2", max_steps=400, timeout=300)
    assert product.get("ok"), product.get("status")

    found = gfn.walk_the_path(start["xyz"], product["xyz"], "gfn2")
    assert found.get("ok"), found.get("status")
    assert 1.0 < found["barrier"] < 12.0, found
    assert found["back"] > 40.0, found
    assert -75.0 < found["reaction"] < -50.0, found
    # How near it came to what it was aimed at, which is the one number that
    # says whether the answer is about the reaction that was asked for.
    assert found["rmsd"] is not None and found["rmsd"] < 0.3, found
    assert found["points"] and found["points"] > 5, found

    # And the transition state it estimates, with the two bonds half made.
    here = gfn.coordinates_of(found["ts"])
    assert len(here) == 48, len(here)
    for pair in ((0, 10), (3, 11)):
        span = math.dist(here[3 * pair[0]:3 * pair[0] + 3],
                         here[3 * pair[1]:3 * pair[1] + 3])
        assert 2.1 < span < 2.9, (pair, span)

    # Two structures of different molecules are not a reaction, and are said
    # so rather than walked between.
    assert not gfn.walk_the_path(_ETHANE, _BROMOBENZENE, "gfn2")["ok"]


def test_a_solvent_the_method_has_not_got_is_not_blamed_on_the_structures():
    """The walk's refusals must be about the walk, not about the chemistry.

    Optimise has refused a solvent under a method that has none for a long
    time; the path finder did not, and ran into it.  Driven: g-xTB handed
    ``--alpb water`` stops with "No ALPB/GBSA parameters found for the
    method/solvent", and what came back was "g-xTB found no path between the
    two structures" -- a sentence about the reaction, where the truth was
    about the build.  Nothing is run here: the answer arrives before the
    binary is asked for anything.
    """
    refused = gfn.walk_the_path(_ETHANE, _ETHANE, "gxtb", solvent="water")
    assert not refused["ok"]
    assert "no implicit solvation" in refused["status"], refused["status"]
    # A solvent nothing knows is still named as one, under a method that has
    # solvation at all.
    unknown = gfn.walk_the_path(_ETHANE, _ETHANE, "gfn2", solvent="unobtainium")
    assert not unknown["ok"]
    assert "not a solvent" in unknown["status"], unknown["status"]


@_needs_xtb
def test_a_chain_that_cannot_walk_never_climbs():
    """The chain is two halves, and the first one is allowed to refuse.

    A path walks atom 1 to atom 1, so two structures that are not the same
    molecule in the same order are not a path at all.  The refusal is the
    walk's and it arrives before any of ORCA's time is spent -- said as the
    walk's refusal, with the half it got to named, rather than as a saddle
    search that failed for reasons of its own.
    """
    from delfin.dashboard import saddle

    got = saddle.path_to_saddle(_ETHANE, _BROMOBENZENE, "gfn2")
    assert not got["ok"]
    assert got["stage"] == "path"
    assert got.get("xyz") is None
    assert "8" in got["status"] and "12" in got["status"], got["status"]


def test_the_path_is_offered_once_there_is_something_to_walk_between():
    """It cannot invent a product to aim at, so it appears when a scan has
    made one -- and what it finds is one step, which Undo takes back whole.

    The offer is a start in the box beside the press rather than two buttons
    of its own.  Both ways of walking between the two ends are still there --
    the path on its own, and the path carried on into a saddle -- and they are
    two answers to the one question the box next to it asks, so neither can
    appear without the other.
    """
    source = EDITOR_SOURCE
    assert 'submit_saddle_from = widgets.Dropdown(' in source
    assert '''("the scan's two ends", 'scan')''' in source
    assert "('the path only', 'walk')" in source
    assert "state['scan_ends'] = (" in source
    assert "if state.get('scan_ends'):" in source
    assert 'def _refresh_saddle_controls():' in source
    assert '_refresh_saddle_controls()' in source
    assert 'def _walk_the_path(ends, then_climb=False):' in source
    assert 'submit_saddle_btn.on_click(on_submit_saddle)' in source
    # Its own answer, kept as its own step.
    assert "_remember('the transition state the path finder estimated')" in source
    assert "'Estimated transition state, from the path'" in source
    # And how near it came, because the finder always reports a barrier and
    # only this says whether it is about the reaction that was asked for.
    assert 'RMSD from what it aimed ' in source


@_needs_xtb
def test_whether_an_estimated_transition_state_is_one():
    """"Estimated transition state" is a phrase; one imaginary frequency is a
    fact.

    A path finder returns an estimate whatever happened, and what makes a
    structure a transition state is one mode going the wrong way and no
    others.  That is a Hessian, and a Hessian is 0.6 s on sixteen atoms --
    nothing beside not knowing.  Measured on what the path finder handed back
    for the Diels-Alder: a single imaginary frequency at -131.4 cm-1, so it
    really is a first-order saddle at this level of theory and worth handing
    to ORCA's OptTS.  An ethane at rest has none.
    """
    rest = gfn.optimize_with_gfn(_ETHANE, "gfn2", max_steps=300, timeout=300)
    assert rest.get("ok"), rest.get("status")
    settled = gfn.optimize_with_gfn(rest["xyz"], "gfn2", timeout=300,
                                    optimise=False, free_energy=True)
    assert settled.get("imaginary") is not None
    assert settled["imaginary"]["count"] == 0, settled["imaginary"]

    # Not asked for, not counted -- it is a Hessian either way.
    plain = gfn.optimize_with_gfn(rest["xyz"], "gfn2", timeout=300,
                                  optimise=False)
    assert plain.get("imaginary") is None

    # And the same mode is not two of them: xtb prints its list more than
    # once, and counted twice a saddle point looks like a second-order one.
    modes = settled["imaginary"]["modes"]
    assert len(modes) == len(set(modes)), modes


def test_the_path_says_whether_what_it_found_is_a_transition_state():
    """And what to do about it either way.

    In the same words the saddle search says it in, and against the same
    thresholds, because it is the same question asked twice: a structure one
    of them names a second-order saddle must not be called a transition state
    by the other.  Both go through
    :func:`delfin.dashboard.saddle.verdict`, which is where the thresholds are
    kept and where the searches they were taken from are cited.

    Without its advice, though.  What to do about a converged second-order
    saddle is to displace along the second mode and climb again; an estimate
    is nobody's stationary point, and a climb from one with two modes going
    the wrong way may perfectly well end on a transition state.
    """
    source = EDITOR_SOURCE
    assert "state['path_shape'] = shape.get('imaginary')" in source
    assert "lines.extend(_said_modes(" in source
    assert "shape, 'The structure it estimates', advise=False))" in source
    # A minimum is not one, and neither is a second-order saddle -- both are
    # said rather than left for the user to notice.
    assert "if shape.get('count') == 0:" in source
    assert 'not sitting on' in source
    # And where it goes next, since xtb has no saddle optimiser at all: both
    # ways on from the estimate are the same press with the box beside it
    # moved, and the line names the two entries rather than two buttons.
    assert 'Set the box beside the press to "through ORCA" ' in source
    assert 'to sharpen it in seconds, or to "by hand" to ' in source


#: cis- and trans-2-butene, relaxed under GFN2.  The plainest reaction there
#: is, and the one a scan cannot do.
_CIS_BUTENE = """12
cis-2-butene
C            1.53845380780890        0.47156751190810        0.24815151049525
C            0.72259292978383       -0.62984247261506       -0.34631766873247
C           -0.59683088037401       -0.72716870418064       -0.39216635686013
C           -1.60437468580189        0.23974033726461        0.13894193453797
H            2.17838880916630        0.91117649333040       -0.51685933694657
H            2.18667179300093        0.07139161982068        1.02769435533968
H            0.92355533810252        1.25571935239631        0.67779792380508
H            1.31338637860703       -1.42624004246114       -0.78248673827834
H           -1.03299293753579       -1.59931845958150       -0.86402045115136
H           -2.24747192458536        0.58470382243324       -0.67065093780288
H           -2.23919720363920       -0.25507502633660        0.87390618327608
H           -1.14218042453326        1.10334556802157        0.60600958231770
"""


#: And trans, relaxed the same way: 1.5 kcal/mol below cis, which is where
#: a scan of the torsion never arrives.
_TRANS_BUTENE = """12
trans-2-butene, from cis by turning one end
C            1.52554923621953        0.45392706276874        0.23862932913643
C            0.72327853899832       -0.65382334827347       -0.35935996805432
C           -0.59751764029312       -0.70318862646274       -0.37912422524052
C           -1.39978867651523       -1.81093942967462       -0.97711418435257
H            2.16442245697151        0.90942661443140       -0.51801711143509
H            2.17270621886454        0.06940653401115        1.02697214268378
H            0.87928852393398        1.22065021589986        0.65896630885336
H            1.29895068609740       -1.45854709844859       -0.79997503795525
H           -1.17319015647323        0.10153484150752        0.06149079451936
H           -2.03866596473115       -2.26643618060440       -0.22046950012781
H           -2.04694150620205       -1.42641991623531       -1.76546095407231
H           -0.75352744687051       -2.57766461891954       -1.39744649395505
"""


@_needs_xtb
def test_a_scan_cannot_turn_a_double_bond_and_nothing_can_tell_it_so():
    """One dihedral does not pin four substituents, and the method cannot
    describe what is at the top.

    Driving C-C=C-C from cis to trans meets every value it is asked for and
    turns nothing: measured under GFN2 on a 2-butene, the constraint is met by
    pyramidalising the carbons -- the C=C goes from 1.327 to 1.477 A and the
    H-C=C-H dihedral to 113 degrees -- and the walk ends at 180 degrees and
    +64 kcal/mol above cis, where real trans lies 1.5 *below*.

    And it ends on a genuine minimum: released, that structure does not move
    at all, 0.00 A and 0.0 kcal/mol.  A twisted alkene is a biradical and a
    closed-shell single determinant cannot describe one, so GFN2 has a
    spurious minimum there -- which is a method failure and not something a
    scan can be taught to notice.  A scan has no target: it drives a
    coordinate and reports what it reaches, and nothing in it knows what the
    answer was supposed to be.

    Which is exactly what the path finder does have, and is why the answer to
    this reaction is to give it the two structures.  This test holds the
    measurement so the next person does not have to make it again.
    """
    part, state, box = _scanned(
        "hold", structure=_CIS_BUTENE, seconds=900,
        legs=[{"kind": "dihedral", "atoms": [0, 1, 2, 3], "from": 0.0,
               "to": 180.0, "steps": 18, "structure": None}])
    said = part.mol_status.value
    ends = float(said.split("\u00b7 end ")[1].split("\u00b7")[0].strip()
                 .lstrip("+"))
    assert ends > 40.0, said            # nowhere near trans, which is -1.5

    # And it really is a minimum of the method, which is why nothing catches
    # it: released, it stays.
    settled = gfn.optimize_with_gfn(box.value, "gfn2", max_steps=400,
                                    timeout=300)
    assert settled.get("ok"), settled.get("status")
    assert gfn.largest_shift(settled["xyz"], box.value) < 0.1, "it stays put"


@_needs_xtb
def test_two_structures_are_a_path_without_a_scan():
    """A great many questions arrive as two structures the user already has,
    and a cis/trans isomerisation is the plainest of them.

    Measured: given cis and trans, the path finder answers in under five
    seconds -- 51 to 55 kcal/mol forward, trans 1.5 below cis, within 0.01 A
    RMSD of what it was aiming at, and an estimated transition state near 85
    degrees, which is where a twisted alkene is.  Against the textbook, the
    isomerisation wants about 480 C and this says about 440.
    """
    cis = gfn.optimize_with_gfn(_CIS_BUTENE, "gfn2", max_steps=400,
                                timeout=300)
    assert cis.get("ok"), cis.get("status")
    trans = gfn.optimize_with_gfn(_TRANS_BUTENE, "gfn2", max_steps=400,
                                  timeout=300)
    assert trans.get("ok"), trans.get("status")
    # The real one, which is below cis -- a scan of the torsion never gets
    # here at all.
    assert (trans["energy"] - cis["energy"]) * 627.5095 < 0.0

    found = gfn.walk_the_path(cis["xyz"], trans["xyz"], "gfn2", points=40,
                              timeout=600)
    assert found.get("ok"), found.get("status")
    assert 40.0 < found["barrier"] < 75.0, found
    assert found["rmsd"] is not None and found["rmsd"] < 0.1, found

    here = gfn.coordinates_of(found["ts"])
    twist = abs(gfn._dihedral(
        [(here[3 * i], here[3 * i + 1], here[3 * i + 2]) for i in range(4)],
        0, 1, 2, 3))
    assert 60.0 < twist < 120.0, twist


def test_a_path_can_be_marked_one_structure_at_a_time():
    """The finder needs a start and an end and cannot invent either.  A scan
    leaves both; two structures the user has are marked one at a time.

    Marking stays a press of its own and does not become an entry in the box
    that chooses the start.  Marking happens once and the box is a setting
    that stands -- the same difference there is between Set and Hold -- and
    the two sources it produces are what that box is for.
    """
    source = EDITOR_SOURCE
    assert 'submit_path_from_btn = widgets.Button(' in source
    # It names which half of the gesture is next: pressed once it is the
    # beginning, and then it asks for the end. "Mark this end" said the same
    # thing before and after, and before there is a second structure nothing
    # visible happens -- which is what made it look like a press that does
    # nothing.
    assert "description='Mark the beginning'" in source
    assert "'Mark the end' if state.get('path_from')" in source
    assert "state['path_from'] = xyz" in source
    assert 'submit_path_from_btn.on_click(on_submit_path_from)' in source
    # Which of the two sources is used is asked rather than assumed: it used
    # to prefer the marked pair silently, so a scan walked after marking
    # something could not be walked between at all.
    # Two marks, not a mark and whatever is on screen: the end used to be
    # whatever happened to be showing, so the pair changed under the user
    # every time they looked at something else.
    assert "first = state.get('path_from')" in source
    assert "second = state.get('path_to')" in source
    assert "if first.strip() == second.strip():" in source
    assert "ends = state.get('scan_ends')" in source
    assert "if which == 'scan':" in source
    assert "('the end you marked', 'marked')" in source
    # And with neither, what to do is said rather than "run a scan first".
    assert 'Mark a beginning, build or load the ' in source


def test_the_method_says_when_it_has_stopped_being_able_to_answer():
    """A closed-shell single determinant describes two electrons in one
    orbital.  At the top of a bond-breaking barrier they are not in one, and
    the energy that comes back is about a state the method cannot represent.

    That cannot be fixed here -- it is what the method is.  It can be
    noticed, because the frontier gap says so: measured under GFN2 on a
    2-butene, 5.28 eV at the cis minimum and 2.42 at the twisted transition
    state, where GFN2 also invents a minimum 64 kcal/mol above cis that real
    trans lies 1.5 below.  A barrier quoted without that is a number somebody
    will use.

    The *fall* is the signal and not the value.  A great many perfectly
    ordinary systems simply have a small gap -- transition-metal complexes,
    long conjugated chains, anything with close-lying frontier orbitals -- and
    they are not in trouble for it.  Keyed on the value alone this would fire
    on every one of them, all the time, which is noise and teaches people to
    ignore it.
    """
    # A comfortable gap, and nothing to say.
    assert gfn.method_is_out_of_its_depth(5.28) == ''
    assert gfn.method_is_out_of_its_depth(4.0, 5.28) == ''
    # A system that simply has a small gap and keeps it: silent.
    assert gfn.method_is_out_of_its_depth(2.0, 2.1) == ''
    # And one that halves from very wide to still-wide: also silent, because
    # nothing has closed.
    assert gfn.method_is_out_of_its_depth(6.0, 12.0) == ''
    # But the same small-gap system when its gap closes further: said.
    assert 'frontier gap' in gfn.method_is_out_of_its_depth(0.8, 2.0)
    # Halved on the way up: said, with both numbers.
    fell = gfn.method_is_out_of_its_depth(2.42, 5.28)
    assert 'frontier gap is 2.4 eV' in fell, fell
    assert 'down from 5.3' in fell, fell
    assert 'worth checking open-shell' in fell, fell
    # Small on its own, with nothing to compare against.
    assert 'frontier gap is 1.2 eV' in gfn.method_is_out_of_its_depth(1.2)
    assert gfn.GAP_IS_SMALL < gfn.GAP_WORTH_SAYING
    # And nothing to read is nothing to say, rather than a warning about None.
    assert gfn.method_is_out_of_its_depth(None) == ''
    assert gfn.method_is_out_of_its_depth('') == ''


@_needs_xtb
def test_the_gap_is_read_from_the_run_that_was_going_to_happen_anyway():
    """No extra calculation for it: every xtb run prints it."""
    got = gfn.optimize_with_gfn(_ETHANE, "gfn2", timeout=300, optimise=False)
    assert got.get("ok"), got.get("status")
    assert got.get("gap") is not None
    # An ethane is nowhere near the edge of anything.
    assert got["gap"] > gfn.GAP_WORTH_SAYING, got["gap"]
    assert gfn.method_is_out_of_its_depth(got["gap"]) == ''


def test_the_gap_is_read_in_both_the_hands_xtb_writes_it_in():
    """xtb prints the gap twice, and g-xTB prints only one of the two.

    A summary block says ``:: HOMO-LUMO gap`` after every single point; a
    property block at the end says ``| HOMO-LUMO GAP``.  Only the second was
    being matched, and g-xTB has no property block -- measured on a propane,
    ``--sp``, ``--opt``, ``--ohess`` and ``--grad`` alike.  So the gap came
    back as nothing under the most accurate method the editor offers, and
    :func:`method_is_out_of_its_depth` -- which says nothing about a gap it
    cannot read -- said nothing about a barrier most likely to be quoted.

    The *last* one, because the summary is printed at every geometry the run
    passes through: on that propane under GFN2, ``--ohess`` prints 14.3656 for
    the structure that went in and 14.5837 for the one that came out, and it
    is the second that the answer is about.  Taking the last reproduces what
    was matched before, to every digit, for GFN2, GFN1 and GFN-FF across all
    four run types.
    """
    both = ('           :: HOMO-LUMO gap             14.365616087375 eV    ::\n'
            '  ...an optimisation happens here...\n'
            '           :: HOMO-LUMO gap             14.583665010763 eV    ::\n'
            '          | HOMO-LUMO GAP              14.583665065106 eV   |\n')
    assert [one.group(1) for one in gfn._GAP_RE.finditer(both)][-1] \
        == '14.583665065106'
    # g-xTB writes the summary line and nothing else, and it must be read.
    summary_only = (
        '           :: HOMO-LUMO gap             18.511149638525 eV    ::\n'
        '           :: HOMO-LUMO gap             18.656166584988 eV    ::\n')
    assert [one.group(1) for one in gfn._GAP_RE.finditer(summary_only)][-1] \
        == '18.656166584988'
    # GFN-FF has no orbitals and prints neither, which stays nothing to say.
    assert not gfn._GAP_RE.findall('| TOTAL ENERGY -5.070325081194 Eh |')


@pytest.mark.skipif(gfn.find_gxtb() is None, reason="g-xTB not installed")
def test_the_gap_comes_back_under_g_xtb_as_well_as_under_gfn2():
    """Driven rather than read: the same ethane through both binaries.

    g-xTB is the method somebody chooses when the number has to be worth
    quoting, so it is the one where a silent warning matters most.
    """
    for method in ('gfn2', 'gxtb'):
        got = gfn.optimize_with_gfn(_ETHANE, method, timeout=600,
                                    optimise=False)
        assert got.get('ok'), got.get('status')
        assert got.get('gap') is not None, method
        # An ethane is nowhere near the edge of anything under either.
        assert got['gap'] > gfn.GAP_WORTH_SAYING, (method, got['gap'])


def test_the_path_says_it_where_the_barrier_is_reported():
    """Beside the number it is about, not in a log."""
    source = EDITOR_SOURCE
    assert "state['path_depth'] = _gfn.method_is_out_of_its_depth(" in source
    assert "shape.get('gap'), start.get('gap'))" in source
    assert "if state.get('path_depth'):" in source
    assert "lines.append(state['path_depth'])" in source


@_needs_xtb
def test_a_scan_says_when_the_method_has_run_out_of_depth():
    """A walk runs into such a region without anyone choosing to.

    There are places on every surface where a closed-shell single determinant
    stops describing what is there -- a bond half broken, a ring opening, a
    double bond turned -- and the number that comes back from one is not wrong
    by a little.  The path finder said so already; a scan is where it matters
    more, because a scan has no target and reports whatever it reaches.

    Measured on the torsion of a 2-butene, which is the plainest example
    rather than the point: the frontier gap falls from 5.26 eV at the start to
    1.44 at the top, and the walk reports +84.6 kcal/mol -- a number somebody
    would otherwise use.

    Kept as the narrowest gap seen and compared against the first point,
    because what matters is where the walk went and not where it started.
    """
    part, state, box = _scanned(
        "hold", structure=_CIS_BUTENE, seconds=900,
        legs=[{"kind": "dihedral", "atoms": [0, 1, 2, 3], "from": 0.0,
               "to": 180.0, "steps": 18, "structure": None}])
    first, least = state.get("scan_gap_first"), state.get("scan_gap_least")
    assert first is not None and least is not None, (first, least)
    assert least < first, (first, least)
    said = part.mol_status.value
    assert "frontier gap" in said, said
    assert "worth checking open-shell" in said, said
    # And it names no chemistry: the rule is about closed-shell determinants,
    # and the system it was measured on belongs in a docstring rather than on
    # the screen of somebody working on a different one.
    for word in ("butene", "alkene", "cis", "trans", "Diels"):
        assert word not in said, (word, said)


@_needs_xtb
def test_the_geometry_that_survives_a_drag_is_one_the_temperature_can_reach():
    """The ceiling is about what is *kept*, and it was about nothing at all.

    Measured before this held, on an ethane whose far methyl is dragged out at
    298 K with 22.3 kcal/mol to spend under a pull: the box was left holding
    the page's own coordinates on every message of the drag, at +15.8, +45.3,
    +74.8, +99.1, +116.9, +129.0 and +136.7 kcal/mol, and at +141.2 once the
    hand let go -- while the line underneath read "+16.2 of 22.3 kcal/mol
    available at 298.15 K".

    The two numbers are about two different molecules.  The price belongs to
    the structure xtb relaxed under the hand, and under a pull that structure
    comes most of the way back -- the relaxed answers really did go +0.1,
    +1.1, +0.2, +4.0, +6.8, +12.2, +16.2 -- while the page's model is simply
    where the cursor is.  Priced one and kept the other, the ceiling had
    nothing to refuse: the wall wrote the box only when it said no, and the
    next message from the page wrote over it.

    Which is also the whole of "sometimes it works at room temperature".
    Under the rigid hand the wall fires early and the follow then stands still
    rather than shake, so from the first refusal onwards nothing was
    calculating at all while the mouse went on writing -- +141.2 kcal/mol
    again, with the status line frozen at "+16.8 of 22.3 available".  Whether
    a drag was stopped came down to whether the last message to land was the
    wall's or the page's, which is a race between xtb at a tenth of a second
    and a hand.
    """
    from delfin.dashboard.structure_editor import (
        _THERMAL_SECONDS, thermal_ceiling)

    start = gfn.optimize_with_gfn(_ETHANE, "gfn2", max_steps=400, timeout=300)
    assert start.get("ok"), start.get("status")
    begin = start["xyz"]
    ceiling = thermal_ceiling(298.15, _THERMAL_SECONDS)

    part = _budgeted(begin, kelvin=298.15, hand="pull")
    kept = _dragged_apart(part, begin, far=2.0)
    # The page asked for a C-C at 3.52 A.  What is left is a bond.
    assert _apart(kept, 0, 1) < 2.0, _apart(kept, 0, 1)
    spent = _costs(kept, part.state["thermal_e0"])
    assert spent <= ceiling + 2.0, f"{spent:+.1f} of {ceiling:.1f}"
    # And the box says which structure this is rather than wearing a claim
    # that was written for coordinates it no longer has.
    said = kept.splitlines()[1].lower()
    assert said.startswith(("within the budget", "past the budget",
                            "back to the last")), said


@_needs_xtb
def test_the_placing_hand_is_not_measured_and_does_not_pretend_to_be():
    """The same drag under the other hand, and it is not refused at all.

    That is the point of the change and it is a loss as well as a fix, so it
    is written down rather than left to be discovered.  Under a placing hand
    the atom goes where the cursor puts it -- that is what placing *is*, it is
    what building a structure needs, and it was asked for -- and the budget
    cannot make an exact statement about it, because what is kept afterwards
    is not the geometry that was priced: measured at 1.4 kcal/mol, +16.8
    priced against +18.2 kept, bounded by _SLIP_LOOSE at 0.25 A rather than by
    zero, and whether the wall landed before the next message came down to a
    race between xtb at 70-170 ms and a page reporting every 16-120 ms.

    So it does not act there and it is not on screen there.  A user who wants
    a placing hand that cannot pull the molecule apart has Keep bonds, which
    reads the bonding off what came back and takes the step back if it
    changed, and which is deliberately not tied to either hand.
    """
    from delfin.dashboard.structure_editor import (
        _THERMAL_SECONDS, thermal_ceiling)

    start = gfn.optimize_with_gfn(_ETHANE, "gfn2", max_steps=400, timeout=300)
    assert start.get("ok"), start.get("status")
    begin = start["xyz"]
    ceiling = thermal_ceiling(298.15, _THERMAL_SECONDS)

    part = _budgeted(begin, kelvin=298.15, hand="move")
    assert part._thermal_live() is False, "nothing reads the switch here"
    anchored = part.coords_widget.value
    _dragged_apart(part, begin, far=2.0)

    # Nothing was priced and nothing was refused.  The same drag under a pull
    # reported "+28.9 kcal/mol -- past the 22.3 this structure has at 298.15
    # K" and handed a C-C of 1.66 A back; here the line is the follow alone.
    said = part.state.get("gfn_last_status") or ""
    assert "follows the drag" in said, said
    for claim in ("kcal/mol", f"{ceiling:.1f}", "budget"):
        assert claim not in said, (claim, said)
    # And the wall never ran while the hand was down: the maximum it carries
    # is still the zero the anchor left it at, and the geometry it would hand
    # back is still the anchor's own.
    assert part.state.get("thermal_peak") == 0.0, part.state.get("thermal_peak")
    assert part.state.get("thermal_good") == anchored


@_needs_xtb
def test_the_same_drag_goes_through_at_a_temperature_that_can_pay_for_it():
    """A ceiling that refuses everything is not a ceiling either.

    The identical sequence of messages at 1500 K, where the hour buys 117.0
    kcal/mol instead of 22.3.  The step that was refused at room temperature
    -- the relaxed structure at +23.5 -- is kept here, and the bond is
    stretched further than room temperature would ever allow.
    """
    from delfin.dashboard.structure_editor import (
        _THERMAL_SECONDS, thermal_ceiling)

    start = gfn.optimize_with_gfn(_ETHANE, "gfn2", max_steps=400, timeout=300)
    begin = start["xyz"]
    cold = thermal_ceiling(298.15, _THERMAL_SECONDS)
    hot = thermal_ceiling(1500.0, _THERMAL_SECONDS)
    assert hot > 100.0 > cold

    part = _budgeted(begin, kelvin=1500.0, hand="pull")
    kept = _dragged_apart(part, begin, far=2.0)
    spent = _costs(kept, part.state["thermal_e0"])
    # Past what room temperature would have allowed, and inside what this
    # temperature does.
    assert spent > cold - 5.0, f"{spent:+.1f} kcal/mol"
    assert spent <= hot, f"{spent:+.1f} of {hot:.1f}"
    assert _apart(kept, 0, 1) > _apart(begin, 0, 1) + 0.2


@_needs_xtb
def test_a_turn_the_temperature_allows_is_not_refused_for_a_body_slide():
    """A held value is an internal coordinate, so the molecule may sit where
    it likes.

    xtb meets the dihedral it was given and is then free to place the whole
    structure anywhere that meets it, so the answer comes back turned and slid
    bodily away from the cursor.  That slide was being measured as a loose
    hold and refused at 0.25 A -- on a butane turned about its middle bond it
    reached 0.27 while the price stood at +0.0 of 22.3 kcal/mol, so the
    temperature was allowing the turn and the drag stopped anyway.  Worse, the
    line said "Past the budget", which sends a reader to the temperature box
    to fix something that is not there.

    Laid back on as a rigid body first -- which costs nothing, an energy does
    not depend on where a molecule is -- the residue is 0.005 A.
    """
    start = gfn.optimize_with_gfn(_ANTI_BUTANE, "gfn2", max_steps=400,
                                  timeout=300)
    assert start.get("ok"), start.get("status")
    begin = start["xyz"]
    end = (0, 4, 5, 6)
    # The gesture: that end of the chain turned about the axis the chain lies
    # along, which is what a hand on a terminal methyl does.
    turned = _turned_about_x(begin, end, 20.0)
    held = gfn.contacts_holding(turned, list(end), most=3, was=begin)
    assert held, "a turned methyl has something to hold"
    out = gfn.relax_steps(turned, method="gfn2", cycles=20, timeout=120,
                          constraints=held)
    assert out.get("ok"), out.get("status")

    raw = gfn.hold_atoms_at(out["xyz"], turned, end)
    slid = gfn.largest_shift(raw, out["xyz"])
    laid = gfn.settle_onto(out["xyz"], turned, end)
    residue = gfn.largest_shift(gfn.hold_atoms_at(laid, turned, end), laid)
    # The body's own freedom is the larger part by far, and it is not a
    # shortfall in the hold.
    assert slid > residue * 5, (slid, residue)
    assert residue < 0.05, residue
    # And it is the residue the editor measures.
    body = EDITOR_SOURCE.split("state['thermal_now'] = priced.get('energy')")[1]
    body = body.split("state['gfn_last_status'] = said")[0]
    assert "slipped = (_gfn.largest_shift(reached, laid)" in body


@_needs_xtb
def test_a_scan_hands_back_a_structure_the_temperature_can_reach():
    """The walk is reported whole; what is left in the box is not.

    A scan is asked to drive a coordinate wherever it is told, and saying how
    high the path went and what temperature would cross it is the answer it
    exists to give.  Leaving that geometry in the box is a different thing --
    it is the structure the user carries on from, and at this temperature they
    cannot get to it.  Otherwise the one place the ceiling is enforced can be
    walked round by pressing a different button.

    Measured on an ethane C-C walked from 1.52 to 3.0 A at 298 K: the walk
    ends at +114.0 kcal/mol against the anchor, and the box is left holding
    the last point that was inside the 22.3 available -- a C-C at 1.767 A
    costing +13.6.
    """
    import time

    from delfin.dashboard.structure_editor import (
        _THERMAL_SECONDS, thermal_ceiling)

    start = gfn.optimize_with_gfn(_ETHANE, "gfn2", max_steps=400, timeout=300)
    begin = start["xyz"]
    part = _budgeted(begin, kelvin=298.15)
    ceiling = thermal_ceiling(298.15, _THERMAL_SECONDS)
    part.state["scan_legs"] = [
        {"kind": "distance", "atoms": [0, 1], "from": 1.52, "to": 3.0,
         "steps": 12, "structure": None},
    ]
    part.submit_scan_how.value = "hold"
    part.on_submit_scan_run()
    began = time.time()
    while part.state.get("scan_run") and time.time() - began < 600:
        time.sleep(0.05)
    assert not part.state.get("scan_run"), "the scan never finished"

    # It went a long way past what the hour buys at this temperature.
    assert part.state.get("scan_walled") is not None
    assert part.state["scan_walled"] > 3 * ceiling, part.state["scan_walled"]
    kept = part.coords_widget.value
    assert _apart(kept, 0, 1) < 2.2, _apart(kept, 0, 1)
    assert _costs(kept, part.state["thermal_e0"]) <= ceiling
    assert kept.splitlines()[1].startswith("Scanned, and back to the last")


@_needs_xtb
def test_a_settle_that_goes_downhill_is_not_refused():
    """The release answers to the ceiling too, and must not block the ordinary
    case.

    A relaxation goes downhill, so this ordinarily has nothing to refuse --
    which is exactly why it is asked rather than assumed: a value the user is
    holding is restored on every step of one, and restoring one is uphill.
    Measured on an ethane anchored while strained, C-C 1.721 A: the settle
    takes it to 1.522 at -10.8 kcal/mol against the anchor, and goes into the
    box.
    """
    import time

    start = gfn.optimize_with_gfn(_ETHANE, "gfn2", max_steps=400, timeout=300)
    strained = _shifted(start["xyz"], {1, 5, 6, 7}, 0.20)
    part = _a_part(strained)
    part.submit_ff_dd.value = "gfn2"
    part.submit_thermal_relax.value = False      # anchor on the strain itself
    part.submit_temperature.value = 298.15
    part.submit_thermal_btn.value = True
    began = time.time()
    while part.state.get("thermal_e0") is None and time.time() - began < 300:
        time.sleep(0.05)
    assert part.state.get("thermal_e0") is not None

    part.submit_settle_btn.value = True
    part._gfn_settle_now()
    began = time.time()
    while part.state.get("gfn_settle_busy") and time.time() - began < 300:
        time.sleep(0.05)
    kept = part.coords_widget.value
    assert kept.splitlines()[1].startswith("Settled with")
    assert _apart(kept, 0, 1) < _apart(strained, 0, 1) - 0.1
    assert _costs(kept, part.state["thermal_e0"]) < 0.0


def test_every_claim_this_editor_writes_can_be_taken_back():
    """A comment the editor wrote is replaced when the coordinates change.

    Only the ones in ``_EDITOR_COMMENTS`` are; anything else in that line
    belongs to the user and is carried over whatever happens underneath it.
    So a claim this file writes and forgets to register becomes a sentence
    about a geometry that is no longer there, with nothing anywhere saying so:
    measured, "Past the budget: back to the last structure that was inside it"
    stood above a torn ethane at +141.2 kcal/mol, because the page's next
    message wrote its own coordinates under the comment it found.

    Read off the source rather than listed here, so that writing a new one and
    not registering it fails in this file rather than in a coordinate box.
    """
    import ast

    from delfin.dashboard import structure_editor
    from delfin.dashboard.structure_editor import _is_editor_comment

    def literals(node):
        """The comment a call writes, however that argument is built.

        An f-string keeps its fixed parts and stands a letter in for whatever
        is interpolated: "Settled with {label}" is registered as "settled with
        " and would not match its own leading fragment, which strips to
        "settled with" and loses the space that makes the prefix mean
        something.
        """
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return [node.value]
        if isinstance(node, ast.JoinedStr):
            return [''.join(
                one.value if isinstance(one, ast.Constant) else 'X'
                for one in node.values)]
        if isinstance(node, ast.IfExp):
            return literals(node.body) + literals(node.orelse)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            return literals(node.left)
        return []

    tree = ast.parse(open(structure_editor.__file__, encoding="utf-8").read())
    written = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if getattr(node.func, "id", "") != "xyz_document":
            continue
        if len(node.args) < 2:
            continue
        written.extend(literals(node.args[1]))
    assert len(written) >= 9, written
    for one in written:
        assert _is_editor_comment(one), one


def test_the_temperature_does_not_set_the_hand_and_no_longer_says_it_does():
    """It once did, and the docstring went on saying so long after it stopped.

    The hand was derived from the ceiling over a reach, which needs a length
    no temperature supplies: sized as a distance it was too weak to turn a
    torsion, so a molecule could not be put into its own conformers at exactly
    the temperature that certainly allows it.  What a temperature limits is
    the energy of what is kept, and the wall is what enforces that.
    """
    part = _a_part()
    part.submit_hand_dd.value = "pull"
    part.submit_pull_slider.value = 0.4
    for kelvin in (150.0, 298.15, 1500.0):
        part.submit_temperature.value = kelvin
        assert part._pull_force() == pytest.approx(0.4 * gfn.A_BOND_HOLDS)
    part.submit_pull_slider.value = 0.8
    assert part._pull_force() == pytest.approx(0.8 * gfn.A_BOND_HOLDS)
    # Nothing caps it at any temperature, which is the same statement.
    assert part._pull_most() is None

    body = EDITOR_SOURCE.split("def _pull_force():")[1]
    body = body.split("\n    def ")[0]
    assert "submit_temperature" not in body
    assert "the temperature sets it" not in body.lower()
    assert "the temperature does not come into it" in body.lower()


def test_the_page_does_not_write_its_own_coordinates_under_a_budget():
    """While something is pricing the drag, the box belongs to it.

    The page's model is where the cursor is, and it landed in the box on every
    message of every drag with nothing on that path ever asking what it cost.
    The wall ran a moment behind in a thread of its own and wrote the box only
    when it refused -- so the geometry that survived a drag was the one thing
    that had never been priced.
    """
    source = EDITOR_SOURCE
    assert "walled = ((dragging or released)" in source
    assert "and state.get('gfn_follow')" in source
    # An undo is not a drag: it hands back a geometry that was already there,
    # so held to the budget it would be answered with the very structure it is
    # undoing.
    assert "released = drag_ended" in source
    assert "and _thermal_budget()[0] is not None)" in source
    assert "if walled:" in source
    assert "_keep_the_priced_geometry()" in source
    # The follow writes what it priced, every step, rather than only when it
    # refuses -- otherwise nothing writes the box and the page fills the gap.
    assert "why = ('Within the budget at '" in source
    assert "_write_coords, xyz_document(rows, why), True)" in source
    # And the wall is asked about the geometry that will survive the step.
    assert "came_back = _thermal_wall(\n                            reached," in source


def test_a_budget_that_cannot_price_a_drag_does_not_stand_there_lit_up():
    """The page only reports a drag while the relaxation switch is down.

    Without those messages nothing runs, the ceiling has nothing to compare
    against, and it would sit there refusing nothing at all -- which is worse
    than being off, because it is off and says it is on.  So the two go on
    together and off together, and the line says which.
    """
    source = EDITOR_SOURCE
    on = source.split("def on_submit_thermal(change):")[1].split("\n    def ")[0]
    assert "if _server_method() and not submit_relax_btn.value:" in on
    assert "submit_relax_btn.value = True" in on
    off = source.split("def on_submit_relax_toggle(change):")[1]
    off = off.split("\n    def ")[0]
    assert "if submit_thermal_btn.value:" in off
    assert "submit_thermal_btn.value = False" in off
    assert "has gone off with it" in off


# ---------------------------------------------------------------------------
# The budget belongs to the hand that pulls
# ---------------------------------------------------------------------------

#: N-nitrosodimethylamine, relaxed with GFN2.  Its N-N rotation is a hindered
#: one with two equivalent minima and a real barrier between them, which is
#: what makes it the right shape for the tests further down; here it is only
#: a molecule with more than three atoms in it.
_NITROSAMINE = """11
N-nitrosodimethylamine, GFN2
C           -1.22971478919819       -0.31972219379423       -0.51910745579444
N            0.00850736080694        0.31228920818566       -0.15895433968563
C            0.90811938305689       -0.39664942551853        0.72445895670150
N            0.23961345405570        1.48707421451206       -0.66313769413187
O            1.29756115077536        2.04026454534117       -0.37046935361620
H           -1.82543404722282       -0.52752067104412        0.37221337786934
H           -1.77713038457114        0.35492985087427       -1.17340843733348
H           -1.04033128224138       -1.26021920368146       -1.04149064049211
H            0.42830215177442       -0.57816072484367        1.68695792258498
H            1.20393443451672       -1.34726557367726        0.27934610102296
H            1.78657156674635        0.23497897436250        0.86359255823555
"""


def _editor(text):
    """One structure editor over a coordinate box of its own.

    Built here rather than reached through a tab: what the wall does is the
    part's, and the tab hands out only the widgets it places.
    """
    import tempfile

    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ("calc", "archive", "office"):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / "calc", archive_dir=room / "archive",
                           office_dir=room / "office")
    ctx.run_js = lambda _script: None
    state = {}
    part = structure_editor.build(
        ctx, state=state, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None, get_smiles_charge=lambda *a, **k: None)
    return part, state


def _shown(widget):
    return (widget.layout.display or "") != "none"


def test_the_budget_goes_with_the_hand_that_pulls():
    """A budget prices what is kept, and a placing hand does not keep it.

    Under a pull the geometry that was priced is the geometry that is kept:
    the answer is laid back onto the atoms the hand never touched, which moves
    nothing chemical.  Under a placing hand the held atoms are put back where
    the cursor had them afterwards, and that is a second geometry -- measured
    at 1.4 kcal/mol, +16.8 priced against +18.2 kept, bounded by _SLIP_LOOSE
    at 0.25 A rather than by zero.  Whether a drag was stopped in time then
    came down to a race between xtb at 70-170 ms and a page reporting every
    16-120 ms.  A budget that cannot be made exact under one hand should not
    claim to be a budget there.

    So it disappears with that hand and comes back with the other, and while
    it is out of sight nothing reads it.  The switch keeps its value across
    the change: a hand is switched constantly while working, and measuring the
    anchor again on every switch is a calculation nobody asked for.
    """
    part, _state = _editor(_NITROSAMINE)
    part.submit_ff_dd.value = "gfn2"
    part.submit_hand_dd.value = "pull"
    part.submit_thermal_btn.value = True

    group = (part.submit_thermal_btn, part.submit_temperature,
             part.submit_thermal_relax, part.submit_thermal_anchor_btn)
    assert all(_shown(w) for w in group)
    assert _shown(part.submit_pull_slider)
    assert part._thermal_live() is True

    part.submit_hand_dd.value = "move"
    assert not any(_shown(w) for w in group), "it has nothing exact to say here"
    assert not _shown(part.submit_pull_slider), "a placement has no force"
    assert part._thermal_live() is False, "and nothing reads it while it is off"
    assert part.submit_thermal_btn.value is True, "the setting is not thrown away"

    # A refresh that belongs to the method must not put them back.  It used to:
    # the pull slider was shown unconditionally there, so choosing a method
    # while placing an atom brought back a control with nothing to act on.
    part._refresh_method_controls()
    assert not any(_shown(w) for w in group)
    assert not _shown(part.submit_pull_slider)

    part.submit_hand_dd.value = "pull"
    assert all(_shown(w) for w in group), "and back with the hand it belongs to"
    assert _shown(part.submit_pull_slider)
    assert part._thermal_live() is True


def test_the_placing_hand_keeps_everything_that_still_works():
    """A control wrongly hidden is as bad as one wrongly shown.

    Placing an atom exactly where it is meant to go is what building a
    structure is, and no force can do it -- the whole point of a force is that
    the chemistry gets a say.  So only what a force is needed for goes.
    """
    part, _state = _editor(_NITROSAMINE)
    part.submit_ff_dd.value = "gfn2"
    part.submit_hand_dd.value = "move"

    # How far the structure moves for how far the mouse moves is a property of
    # dragging, not of what dragging does.
    assert _shown(part.submit_sens_slider)
    # Keep bonds judges what a step did to the bonding; a rigid hand can put an
    # atom anywhere at all, so if anything it needs it more.
    assert _shown(part.submit_topology_btn)
    # A scan drives its own ramp of forces and never reads the hand's slider,
    # so it walks the same path whichever hand is chosen.
    #
    # With a selection, because arming answers to that now: + Add with nothing
    # picked can only say "pick 2, 3 or 4 atoms first", and a press that is
    # hidden for its own good reason is no witness to a hand taking things
    # away.
    _state['picked'] = [0, 1]
    part._refresh_scan()
    assert _shown(part.submit_scan_add_btn)
    assert _shown(part.submit_scan_run_btn)
    assert "submit_pull_slider" not in EDITOR_SOURCE.split(
        "def _push_target(")[1].split("\n        def ")[0]

    # Strength is how many steps the browser's own field takes per animation
    # frame -- the relaxation, not the hand -- and the rest of a structure
    # settles round a placed atom exactly as it settles round a pulled one.
    part.submit_ff_dd.value = "uff"
    part.submit_hand_dd.value = "move"
    assert _shown(part.submit_strength_slider)

    # And the rule is in one place, so the method's refresh and the hand's
    # cannot undo one another.
    body = EDITOR_SOURCE.split("def _refresh_hand_controls")[1].split(
        "\n    def ")[0]
    for spare in ("submit_sens_slider", "submit_topology_btn",
                  "submit_scan_add_btn", "submit_strength_slider"):
        assert spare not in body, spare


# ---------------------------------------------------------------------------
# And it belongs to the whole way the drag came, not to where it is standing
# ---------------------------------------------------------------------------

#: The four atoms of that nitrosamine's N-N rotation, in the order the
#: geometry above has them.
_TURN = (4, 3, 1, 0)


def _armed(part, state, text, energy, T=298.15, method="gfn2"):
    """The budget switched on with a known anchor, and a hand that pulls."""
    part.submit_ff_dd.value = method
    part.submit_hand_dd.value = "pull"
    part.submit_temperature.value = T
    part.coords_widget.value = text
    state["thermal_e0"] = float(energy)
    state["thermal_method"] = method
    state["thermal_for"] = part._structure_fingerprint(text)
    state["thermal_good"] = text
    state.pop("thermal_good_peak", None)
    state.pop("thermal_peak", None)


def _pulled(text, far):
    """The same molecule with one atom pulled *far* Angstrom out.

    A rigid shift would not do: the wall lays one geometry onto the other as a
    rigid body before comparing them, because an energy does not depend on
    where a molecule is and neither does whether it is the same structure.
    Only a change *inside* the molecule makes it a different place to be
    standing.
    """
    rows = [line.split() for line in gfn.atom_lines(text)]
    rows[0][1] = f"{float(rows[0][1]) + float(far):.6f}"
    return (f"{len(rows)}\npulled {far:.2f} A\n"
            + "\n".join(" ".join(one) for one in rows) + "\n")


def _stretched(part, text, pair, by):
    """The same molecule with one internal coordinate longer by *by*.

    A rigid shift of one atom along x is not that: which way it takes the
    coordinate the hand is driving depends on where the molecule happens to be
    pointing, and half the time it takes it the other way.  What the budget
    reads is the coordinate, so the coordinate is what is set here.
    """
    import math

    flat = gfn.coordinates_of(text)
    a, b = int(pair['atoms'][0]), int(pair['atoms'][1])
    away = [flat[3 * a + i] - flat[3 * b + i] for i in range(3)]
    span = math.sqrt(sum(v * v for v in away)) or 1.0
    rows = [line.split() for line in gfn.atom_lines(text)]
    for i in range(3):
        rows[a][1 + i] = f"{flat[3 * a + i] + away[i] / span * float(by):.6f}"
    return (f"{len(rows)}\nstretched {by:.2f} A\n"
            + "\n".join(" ".join(one) for one in rows) + "\n")


def _dihedral(text, a, b, c, d):
    """The torsion in degrees, for setting up a drive that means something."""
    import math

    flat = gfn.coordinates_of(text)
    p = [flat[i:i + 3] for i in range(0, len(flat), 3)]
    b0 = [p[b][i] - p[a][i] for i in range(3)]
    b1 = [p[c][i] - p[b][i] for i in range(3)]
    b2 = [p[d][i] - p[c][i] for i in range(3)]

    def cross(u, v):
        return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0]]

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    n1, n2 = cross(b0, b1), cross(b1, b2)
    return math.degrees(math.atan2(
        dot(cross(n1, b1), n2) / math.sqrt(dot(b1, b1)), dot(n1, n2)))


def test_the_wall_refuses_on_the_highest_point_the_drag_has_been_at():
    """A ceiling is a barrier height, so a point is not what to hold against it.

    A drag that crossed +32 and settled at 0 has still crossed +32, and at 298
    K -- where the hour buys 22.3 -- it could not have.  Priced on the point
    alone the far side was kept: the wall refused the step at the top, the next
    answer came back cheap, and that geometry became the one to fall back to.

    The numbers are the measured ones from the live run further down, so what
    this fixes is stated in the units it was found in.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -17.172582)

    def at(kcal):
        return -17.172582 + kcal / 627.5094740631

    def step(kcal, far):
        return part._thermal_wall(_pulled(_NITROSAMINE, far), at(kcal), [])

    # Up the near side, inside the budget the whole way, the hand travelling
    # as it goes -- a drag that never moves is not a drag.
    assert step(0.0, 0.0) is None
    assert step(13.4, 0.6) is None
    assert step(19.9, 1.2) is None
    # Over the top, and refused there as it always was.
    assert step(26.5, 1.8) is not None
    assert step(32.2, 2.4) is not None
    # And now the far side, which is cheap.  This is the one that was kept.
    assert step(13.7, 3.0) is not None
    assert step(0.0, 3.6) is not None, "it got here through what it cannot pay"
    # What the line quotes is the crossing, not where the structure stands.
    said = part._thermal_note(at(0.0))
    assert "went through +32.2" in said, said
    assert state["thermal_over"] == "path"
    # And the sentence about what was done says the same thing, so the two
    # halves of the line cannot disagree about why the structure sprang back.
    assert "past the budget on the way here -- the last " in EDITOR_SOURCE
    assert "state.get('thermal_over') == 'path'" in EDITOR_SOURCE


def test_an_excursion_that_is_taken_back_stops_refusing():
    """The failure mode in the other direction, which is the worse one.

    A maximum that never comes down makes the editor unusable: one jerk of the
    mouse past the ceiling and every later step of that grab is refused,
    including the ones walking straight back home.  What undoes a crossing is
    coming back over it, so the maximum drops to what was recorded for the last
    geometry the budget agreed to as soon as the structure is standing there
    again -- which, after a rollback, is exactly where it has just been put.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -17.172582)

    def at(kcal):
        return -17.172582 + kcal / 627.5094740631

    out = _pulled(_NITROSAMINE, 0.0)
    # Out past the 0.25 A that counts as standing in the same place, which is
    # the threshold the loose hold is already judged by.
    far = _pulled(_NITROSAMINE, 2.0)

    assert part._thermal_wall(out, at(5.0), []) is None
    # Over the ceiling, out where the hand has taken it: refused, and the last
    # affordable geometry comes back.
    assert part._thermal_wall(far, at(40.0), []) == out
    # Still out there and cheap: still refused, because it never came back.
    assert part._thermal_wall(far, at(1.0), []) == out
    # And back where the budget last agreed: allowed again, and the drag lives.
    assert part._thermal_wall(out, at(5.0), []) is None
    assert part._thermal_wall(out, at(12.0), []) is None


def test_letting_go_forgets_the_way_the_drag_came():
    """What a release leaves behind was inside the budget the whole way.

    Its own way here was affordable end to end -- that is what being kept
    means -- so the next grab starts from where it is standing rather than
    from the last one's worst moment.  Kept, one bad drag would have refused
    every drag after it until the anchor was set again.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -17.172582)
    state["thermal_peak"] = 40.0
    state["thermal_good_peak"] = 40.0
    part._clear_thermal_wall()
    assert "thermal_peak" not in state
    assert "thermal_good_peak" not in state
    # And it is the release that calls it, not only the walled case.
    off = EDITOR_SOURCE.split("def _end_gfn_follow()")[1].split("\n    def ")[0]
    assert "_clear_thermal_wall()" in off


def test_a_new_anchor_is_a_new_path():
    """Set here is pressed on the intermediate a drag has just reached.

    That is what it is for: something that has crossed a barrier thermalises,
    and the next elementary step is measured from where it landed.  Carrying
    the barrier it came over into the new budget would refuse everything from
    the moment it was pressed.
    """
    body = EDITOR_SOURCE.split("def _set_thermal_anchor(")[1].split(
        "\n    def ")[0]
    assert "state['thermal_peak'] = 0.0" in body
    assert "state['thermal_good_peak'] = 0.0" in body
    assert "state['thermal_good'] = _current_xyz() or xyz" in body


@_needs_xtb
def test_the_wall_prices_the_path_a_real_drag_took():
    """The whole of it, against xtb, on a rotation that really is hindered.

    N-nitrosodimethylamine, GFN2, the N-N torsion driven with everything else
    relaxed -- which is the question a drag under the budget asks.  Measured in
    15 degree steps from one minimum to the equivalent one:

        deg    -180   -150   -120    -90    -75    -60    -30      0
        kcal    0.0    3.6   13.4   26.5   32.2   13.7    3.5    0.0

    At 298.15 K the hour buys 22.3 kcal/mol.  The far minimum costs nothing at
    all, so a wall that prices the point keeps it -- and a structure that got
    there over a 32 kcal/mol barrier is one room temperature cannot make
    within the hour.  The decision is what is asserted and not the numbers:
    xtb's optimiser is not bit-reproducible under threading, and the same
    eight-step drag has been measured at +12.2 once and +28.9 another time.
    Every margin here is a third of the ceiling or more.
    """
    from delfin.dashboard.structure_editor import thermal_ceiling

    began = gfn.optimize_with_gfn(_NITROSAMINE, "gfn2", timeout=None)
    assert began.get("ok"), began.get("status")
    anchor, here = float(began["energy"]), began["xyz"]
    was = _dihedral(here, *_TURN)

    path, standing = [], here
    for n in range(13):
        got = gfn.optimize_with_gfn(
            standing, "gfn2", timeout=None,
            constraints=[{"kind": "dihedral", "atoms": list(_TURN),
                          "mode": "fix", "value": was + 15.0 * n}])
        assert got.get("ok"), got.get("status")
        standing = got["xyz"]
        path.append((standing, float(got["energy"])))

    kcal = [(e - anchor) * 627.5094740631 for _x, e in path]
    ceiling = thermal_ceiling(298.15, 3600.0)
    assert max(kcal) > ceiling * 1.3, kcal
    assert kcal[-1] < ceiling * 0.5, kcal[-1]

    # Room temperature: the far side is refused, because of how it was reached.
    part, state = _editor(here)
    _armed(part, state, here, anchor, T=298.15)
    verdicts = [part._thermal_wall(x, e, []) is None for x, e in path]
    assert verdicts[0] is True, "the anchor itself is affordable"
    assert verdicts[-1] is False, "and the far side was reached through the top"
    assert "went through" in part._thermal_note(path[-1][1])

    # The same drag where the temperature can pay for it.  A ceiling that
    # refused this too would not be a ceiling, it would be a rule.
    part, state = _editor(here)
    _armed(part, state, here, anchor, T=1500.0)
    assert thermal_ceiling(1500.0, 3600.0) > max(kcal) * 1.3
    assert all(part._thermal_wall(x, e, []) is None for x, e in path), (
        "1500 K buys the whole rotation")


@_needs_xtb
def test_a_drag_that_stays_cheap_is_never_taken_back():
    """The other half of the pair, and the one that would be easy to break.

    A rule that refuses a barrier it never crossed is worse than no rule.  Two
    drags room temperature certainly does allow, both under GFN2 with
    everything but the driven coordinate relaxed: an ethane C-C lengthened
    from 1.53 to 1.65 A, which costs +4.4 kcal/mol at the far end, and the same
    ethane turned from staggered to eclipsed, which costs +2.6 -- the whole of
    a torsional barrier that at room temperature turns over some 10^11 times a
    second.
    """
    from delfin.dashboard.structure_editor import thermal_ceiling

    began = gfn.optimize_with_gfn(_ETHANE, "gfn2", timeout=None)
    assert began.get("ok"), began.get("status")
    anchor, here = float(began["energy"]), began["xyz"]
    ceiling = thermal_ceiling(298.15, 3600.0)

    for leg, values in (
        ({"kind": "distance", "atoms": [0, 1]},
         [1.53, 1.56, 1.59, 1.62, 1.65]),
        ({"kind": "dihedral", "atoms": [2, 0, 1, 5]},
         [_dihedral(here, 2, 0, 1, 5) + 15.0 * n for n in range(5)]),
    ):
        path, standing = [], here
        for value in values:
            got = gfn.optimize_with_gfn(
                standing, "gfn2", timeout=None,
                constraints=[dict(leg, mode="fix", value=float(value))])
            assert got.get("ok"), got.get("status")
            standing = got["xyz"]
            path.append((standing, float(got["energy"])))
        kcal = [(e - anchor) * 627.5094740631 for _x, e in path]
        assert max(kcal) < ceiling * 0.5, (leg, kcal)

        part, state = _editor(here)
        _armed(part, state, here, anchor, T=298.15)
        assert all(part._thermal_wall(x, e, []) is None for x, e in path), (
            leg, kcal)


def test_the_drag_is_priced_with_an_electronic_energy_and_says_so():
    """The ceiling is a dG; what is held against it is a dE.

    They part company by T*dS.  While the drag leaves the structure in as many
    separate pieces as it found it that is small -- a torsion, an angle, a ring
    turning over -- and the two agree closely.  Where the number of pieces
    changes it is large and signed: taking something apart releases
    translation and rotation, so T*dS is of order ten kcal/mol at room
    temperature and the budget is too strict there; putting two things
    together is the same number the other way round and it is too lenient.

    Nothing corrects for it, and that is a decision.  A Hessian per drag step
    is not affordable -- measured under GFN2, 0.69 s on 21 atoms and 3.90 s on
    62, against 58 and 321 ms for the single point beside it, twelve times the
    cost at both sizes and one per answer of a control that reports several
    times a second -- and the
    only cheap guess at dS is how many pieces the structure is in, which is
    read off the same distance threshold already written down as flickering in
    a crowded coordination sphere, depends on a standard state nothing here
    knows, and is exactly the case where the method itself is least reliable.
    An invented number would be worse than the gap it papered over.

    So the gap is said.  And it is said in words that assume no chemistry:
    every kind of system is computed here.
    """
    from delfin.dashboard.structure_editor import _THERMAL_QUANTITY_SHORT as note

    assert "free energy" in note and "electronic" in note
    assert "scan" in note, "and where the free-energy answer can be had"
    # Universal: no class of compound, no reaction type, no element.
    for word in ("bond", "molecule", "amide", "ring", "metal", "alkene",
                 "organic", "protein", "carbon"):
        assert word not in note.lower(), word

    source = EDITOR_SOURCE
    # On the button, and on both lines that quote the ceiling.
    assert "+ _THERMAL_QUANTITY_SHORT" in source
    assert source.count("_THERMAL_QUANTITY_SHORT") >= 4
    # And the scan is where a free energy really is taken, at the three places
    # it is both affordable and meaningful.
    assert "is taking three Hessians for the free " in source



def _page_model(part):
    """The browser's own model of the structure, as far as this matters.

    Under a pull the frames the kernel sends are written to every atom, the
    dragged one included -- ``heldSerials`` hands back nothing when ``ffPull``
    is on, because how far the atom got is the answer to the question the drag
    asked.  So the page's model is the last answer, and the wish the cursor
    carries is laid on top of it.
    """
    import json

    held = {"rows": gfn.coordinates_of(part.coords_widget.value),
            "answers": 0}

    def frame(change):
        try:
            said = json.loads(str(change.get("new") or ""))
        except ValueError:
            return
        frames = said.get("frames") or []
        if frames and len(frames[-1]) == len(held["rows"]):
            held["rows"] = [float(one) for one in frames[-1]]
            held["answers"] += 1

    def box(change):
        rows = gfn.coordinates_of(str(change.get("new") or ""))
        if len(rows) == len(held["rows"]):
            held["rows"] = rows

    part.submit_gfn_frame.observe(frame, names="value")
    part.coords_widget.observe(box, names="value")
    return held


def _apart_in(rows, i, j):
    import math

    return math.dist(rows[3 * i:3 * i + 3], rows[3 * j:3 * j + 3])


def _drag_at(part, held, grabbed, toward, step, turns):
    """Grab, move the cursor *turns* times towards another atom, let go.

    Every message the page sends and nothing else: ``gfngrab``, one
    ``DELFIN drag-follow`` per mouse move, ``DELFIN drag-end``, ``gfnfree``.
    Hands back the contact after each answer.
    """
    import math

    symbols = [line.split()[0] for line in gfn.atom_lines(
        part.coords_widget.value)]

    def block(rows, note):
        body = "\n".join(
            f"{symbols[i]} {rows[3*i]:.8f} {rows[3*i+1]:.8f} {rows[3*i+2]:.8f}"
            for i in range(len(symbols)))
        return f"{len(symbols)}\n{note}\n{body}"

    part.submit_cmd_sync.value = f"gfngrab:{id(held) % 9973}:0"
    trail = []
    for turn in range(1, turns + 1):
        rows = list(held["rows"])
        here = rows[3 * grabbed:3 * grabbed + 3]
        there = rows[3 * toward:3 * toward + 3]
        arm = [b - a for a, b in zip(here, there)]
        far = math.sqrt(sum(one * one for one in arm)) or 1.0
        for n in range(3):
            rows[3 * grabbed + n] = here[n] + arm[n] / far * step
        part.submit_manip_sync.value = block(
            rows, f"DELFIN drag-follow held={grabbed} n={turn}")
        _quiet(part.state)
        if turn == 1:
            # What the first answer of this drag was told to hold, read
            # before the release clears it -- the wall forgets the way a
            # finished drag came, and this is the way it came.
            held["opened_on"] = [dict(one) for one
                                 in (part.state.get("thermal_holding") or ())]
        trail.append(_apart_in(held["rows"], grabbed, toward))
    part.submit_manip_sync.value = block(held["rows"], "DELFIN drag-end")
    part.submit_cmd_sync.value = f"gfnfree:{id(trail) % 9973}:"
    _quiet(part.state)
    return trail


@_needs_xtb
def test_a_second_grab_carries_on_from_where_the_structure_is():
    """With Auto off a release leaves the structure where the hand put it, and
    the next grab has to continue from there.

    It did not.  Reported from the viewer: drag, let go -- it stays, which is
    right -- grab the same atom again and it springs back towards the molecule
    before it starts following, and then goes on past the point it had been
    at.  The same gesture in two halves did not equal the gesture in one,
    which for an instrument meant to walk over a surface is not a cosmetic
    complaint.

    The cause is one line's worth, and it is not the pull.  Every answer of a
    drag sets ``thermal_was`` to the geometry it handed back, and the first
    answer had none -- and having none is not neutral.  With no geometry to
    compare against, :func:`gfn_optimize.contacts_holding` cannot see what the
    hand has changed and falls back to the nearest contact, which for a
    grabbed atom is usually a bond it is not driving; :func:`as_pushes` then
    asks for that bond at the length it already has, which is no force at all.
    So the first answer of every drag was a free relaxation with the
    coordinate the hand *is* driving left loose.

    Measured on butadiene and ethylene, an ethene carbon dragged at a diene
    terminus, six mouse moves each way.  Before: the first answer of each drag
    was told to hold the ethene's own C10=C11 at the length it already stood
    at -- 1.373 A on the first grab and 1.387 on the second.  The first drag
    closed 3.353 A to 2.772; the first answer of the second grab put it back
    to **3.033**, and three more answers went on
    recovering ground the hand had already covered -- six answers netted
    0.096 A where the first six netted 0.581.  After: 2.733 A after the first
    drag, then 2.695, 2.669, 2.650, 2.638, 2.629, 2.624 through the second,
    every one of them the way the hand is pulling.  Twelve answers in one drag
    land at 2.626 against 2.624 in two, so the gesture in two halves is the
    gesture in one to within a thousandth of an angstrom.
    """
    part = _a_part(_DIELS_ALDER)
    part.submit_ff_dd.value = "gfn2"
    part.submit_hand_dd.value = "pull"
    part.submit_relax_btn.value = True
    part.submit_auto_btn.value = False
    held = _page_model(part)
    began = _apart_in(held["rows"], 10, 0)

    first = _drag_at(part, held, 10, 0, 0.18, 6)
    # The coordinate the hand is driving is held from the very first answer,
    # rather than the bond the grabbed atom happens to hang on.
    opened = held.get("opened_on") or []
    assert opened and sorted(opened[0]["atoms"]) == [0, 10], opened
    assert first[0] < began - 0.05, (began, first)

    second = _drag_at(part, held, 10, 0, 0.18, 6)
    # Not one answer of the second drag moves away from where the hand is
    # pulling, and the first of them least of all.
    assert second[0] < first[-1], (first[-1], second[0])
    for before, after in zip([first[-1]] + second, second):
        assert after <= before + 1e-6, (first, second)


# --- what a refused drag costs, and what it says ---------------------------


def _drag_block(part, rows, note):
    """A geometry as the page sends one, with the reason on the comment line."""
    symbols = [line.split()[0] for line
               in gfn.atom_lines(part.coords_widget.value)]
    body = "\n".join(
        f"{symbols[i]} {rows[3*i]:.8f} {rows[3*i+1]:.8f} {rows[3*i+2]:.8f}"
        for i in range(len(symbols)))
    return f"{len(symbols)}\n{note}\n{body}"


def _pull_into_the_wall(part, held, grabbed, root, step=0.80, tries=14):
    """Drag *grabbed* straight out of *root* until the budget refuses.

    The cursor is put *step* ahead of where the atom actually is, which is what
    the browser sends: under a pull it writes the held atom at the point the
    hand is asking for, clamped to a reach ahead of the atom itself.
    """
    part.submit_cmd_sync.value = f"gfngrab:{id(held) % 9973}:0"
    for turn in range(tries):
        rows = list(held["rows"])
        here = rows[3 * grabbed:3 * grabbed + 3]
        stem = rows[3 * root:3 * root + 3]
        arm = [a - b for a, b in zip(here, stem)]
        far = math.sqrt(sum(one * one for one in arm)) or 1.0
        for n in range(3):
            rows[3 * grabbed + n] = here[n] + arm[n] / far * step
        part.submit_manip_sync.value = _drag_block(
            part, rows, f"DELFIN drag-follow held={grabbed} n={turn}")
        _quiet(part.state)
        if part.state.get("thermal_spent"):
            return True
    return False


@_needs_xtb
def test_a_refused_drag_starts_no_work_and_writes_no_line():
    """Standing still must not cost a worker and a redraw per mouse move.

    The rule that stops the kernel fighting a hand it has already refused --
    :func:`_still_spent` -- was asked inside the worker, so standing still
    meant beginning a thread, turning the ring on, finding there was nothing
    to compute and ending the thread again.  Once for every report the page
    sends, and the page reports at the rate a mouse moves.

    Measured on the user's manganese complex under GFN2 at charge +1, 298 K
    and a 22.3 kcal/mol ceiling, the phenolate oxygen pulled off the metal and
    then held against the wall for thirty seconds, driven through the page's
    own messages -- gfngrab, one drag-follow per mouse move, no release:

                               budget off   budget on, refusing
        reports from the page        1147                  1625
        workers started                 1                  1277
        status lines written            9                  2555
        frames drawn                    8                     1

    Forty-three threads a second, and the line that lies over the picture
    rewritten eighty-five times a second.  It alternated between exactly two
    values -- the same words with the ring and without it, 179 of each in
    eight seconds -- because nothing was being computed for the words to
    change.  So the picture stood perfectly still and the ring on top of it
    blinked twenty times a second, which is the shaking the budget adds.

    With the budget off the same hand costs one worker: xtb is genuinely busy,
    and every report after the first is folded into ``gfn_follow_xyz`` instead
    of starting anything.

    Asked before the worker is begun the rule costs no thread at all -- it is
    arithmetic over the atoms the hand is on.  Measured again over the same
    thirty seconds: 1698 reports, one worker, three status lines.
    """
    part = _budgeted(_ETHANE, kelvin=100.0)
    part.submit_pull_slider.value = 1.2
    held = _page_model(part)
    assert _pull_into_the_wall(part, held, 2, 0), "the wall never refused"

    writes = []
    part.mol_status.observe(lambda _c: writes.append(1), names="value")
    frames = []
    part.submit_gfn_frame.observe(lambda _c: frames.append(1), names="value")

    # The hand goes on pulling and the mouse goes on reporting, which is what
    # a hand that has met a wall does.  A real one never repeats a pixel, so
    # every message is a different geometry and none of them is dropped by
    # traitlets for being the same text.
    here = held["rows"][6:9]
    stem = held["rows"][0:3]
    arm = [a - b for a, b in zip(here, stem)]
    far = math.sqrt(sum(one * one for one in arm)) or 1.0
    for turn in range(30):
        rows = list(held["rows"])
        for n in range(3):
            rows[6 + n] = here[n] + arm[n] / far * (0.80 + 0.002 * turn)
        part.submit_manip_sync.value = _drag_block(
            part, rows, f"DELFIN drag-follow held=2 n={100 + turn}")
    _quiet(part.state)

    assert not writes, (
        f"a refused drag rewrote the status line {len(writes)} times")
    assert not frames, "a refused drag drew a frame"


@_needs_xtb
def test_a_refusal_reaches_the_page_as_a_stop_the_hand_can_feel():
    """The kernel stops answering, so the picture must stop promising.

    Nothing told the page a step had been refused.  The kernel stands still --
    there is nothing left to compute until the hand eases -- while the page
    went on running the wish out with the cursor: the band grew, the
    coordinates it reported went on changing, and the one thing that had
    actually happened, that the drag had reached its ceiling, was the one
    thing not on screen.

    So the held atom is marked where the budget last agreed it could stand and
    may go no further out than the hand had already run.  The page has had the
    rule for this all along -- ``thermalWallBlocks`` refuses a step only when
    it is both outside the reach and going further out -- and nothing ever
    armed it: ``_push_thermal_wall`` was called in one place in the whole
    editor, with ``None``.

    Coming back in is never blocked, and coming back in is exactly the event
    ``_still_spent`` reads as the hand easing off, so the two sides let the
    same gesture through instead of holding two opinions about it.
    """
    import json

    part = _budgeted(_ETHANE, kelvin=100.0)
    part.submit_pull_slider.value = 1.2
    held = _page_model(part)
    walls = []
    part.submit_gfn_wall.observe(
        lambda change: walls.append(json.loads(str(change["new"]))),
        names="value")
    assert _pull_into_the_wall(part, held, 2, 0), "the wall never refused"

    armed = [one for one in walls if one.get("wall")]
    assert armed, "the refusal never reached the page"
    stop = armed[-1]
    # The atom the hand is on, and nothing else: a mark on an atom nobody is
    # moving would hold back a hand that is not being refused.
    assert list(stop["wall"]) == ["2"], stop
    assert stop["reach"] > 0, stop
    # And the mark is where the budget agreed the atom could stand, which is
    # the structure the wall handed back.
    good = gfn.coordinates_of(part.state["thermal_good"])
    assert stop["wall"]["2"] == pytest.approx(good[6:9], abs=1e-9), stop
    # The reach is where the wish had run to, so the band stops rather than
    # snapping back to the atom: what the hand was asking for when it was
    # refused is a true thing to leave on the screen.
    assert stop["reach"] == pytest.approx(0.80, abs=0.05), stop


def test_the_standing_still_rule_is_asked_before_a_worker_is_spent():
    """And the worker keeps its own copy, for a report that beat the answer.

    Two places, one rule.  The page reports at the rate a mouse moves and xtb
    answers seconds later, so a report can arrive while a run is going -- that
    one is folded into ``gfn_follow_xyz`` and the worker meets it on its next
    turn, which is why the check inside the loop stays.  What must not happen
    is a *thread* per report, which is what asking only inside the loop meant
    once the answer was instant.
    """
    follow = EDITOR_SOURCE.split("def _gfn_follow_step(", 1)[1].split(
        "_start_background(_work, 'The relaxation under the hand')", 1)[0]
    head, inside = follow.split("def _work():", 1)
    assert "if _still_spent(xyz, holding):" in head
    # Before the report is taken in at all, so a refused drag does not even
    # queue a geometry for a run that is not going to happen.
    assert head.index("_still_spent") < head.index("state['gfn_follow_xyz']")
    assert "if _still_spent(current, holding):" in inside


def test_the_page_is_told_when_the_wall_refuses_and_when_it_stops():
    """One place decides, and it decides both ways.

    A stop that is armed and never taken down is worse than no stop: the next
    step the budget allows would be held back by the last one it refused, and
    nothing on screen would say why.  So the same call answers the allowed
    case, and it answers it by taking the stop away.
    """
    follow = EDITOR_SOURCE.split("def _gfn_follow_step(", 1)[1].split(
        "_start_background(_work, 'The relaxation under the hand')", 1)[0]
    assert "_stop_the_hand_at(" in follow
    assert "came_back if came_back is not None" in follow

    stop = EDITOR_SOURCE.split("def _stop_the_hand_at(", 1)[1].split(
        "\n    def ", 1)[0]
    # Allowed: the stop comes down, and only when one was up -- a widget
    # written on every answer is a message a second the page has to read.
    assert "if came_back is None:" in stop
    assert "_push_thermal_wall(None)" in stop
    # Refused: marks for the atoms the hand is on, and the reach the wish had
    # run to.  Never coordinates -- the wall is not a frame writer.
    assert "_push_thermal_wall(marks, reach)" in stop
    assert "submit_gfn_frame" not in stop
    assert "_write_coords" not in stop


def test_a_refusal_says_how_hot_and_not_only_how_long():
    """The ceiling is an inverted rate, so the same arithmetic run the other
    way turns a price already in hand into the temperature that pays for it.

    It costs nothing -- the energy is the expensive half of that sum and it
    has been paid for -- and it is a better answer than a refusal: "past the
    budget at 298 K" is where a user stops, "it wants about 406 K" is where
    they go next, and the second is the question they came with.  Measured
    against the hour: +11 kcal/mol wants 150 K, +30.6 wants 406, +45 wants
    591, +68 wants 883, and past about 402 nothing under 5000 K does it.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -17.172582)

    assert "~150 K (-123 C)" in part._thermal_wants(11.0)
    assert "~406 K (+133 C)" in part._thermal_wants(30.6)
    assert "~591 K (+318 C)" in part._thermal_wants(45.0)
    assert "~883 K (+610 C)" in part._thermal_wants(68.0)
    # And where there is no answer, that is said rather than a number printed.
    said = part._thermal_wants(500.0)
    assert "no temperature under 5000 K" in said, said
    # The window is named, because a temperature without one is not an answer:
    # 22.3 kcal/mol is 298 K for an hour and 250 K for a year.
    assert "within an hour" in part._thermal_wants(45.0)


def test_both_refusals_say_the_temperature_they_want():
    """The one about where the structure is standing, and the one about what
    the drag went through to put it there.

    Two different numbers refuse, so two different temperatures answer, and a
    line that quoted the wrong one would send the user to a temperature that
    does not open the way they came.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -17.172582)

    def at(kcal):
        return -17.172582 + kcal / 627.5094740631

    # Standing past the ceiling: the temperature is worked out from where the
    # structure stands, which is the number being refused.
    said = part._thermal_note(at(45.0))
    assert "+45.0 kcal/mol, past the 22.3" in said, said
    assert "~591 K (+318 C)" in said, said

    # Standing somewhere cheap, having crossed something that was not: the
    # crossing is what is refused, so the crossing is what is answered.
    state["thermal_peak"] = 25.0
    said = part._thermal_note(at(0.0))
    assert "went through +25.0" in said, said
    assert "~333 K (+60 C)" in said, said
    # And the waiting time stands beside it, because they are two questions:
    # how hot within the window, and how long at the temperature that is set.
    # +25 is four days at 298 K and inside the hour at 333, which is the pair
    # of answers a chemist chooses between.
    assert "~3.94 d" in said, said


def test_a_waiting_time_nobody_reads_is_said_in_words():
    """A refusal came out as "that is about 3.56e+29 years", which is a
    number in the wrong clothes.

    The same fault the fast end of this was fixed for, at the other end: the
    sentence it was written to produce is "longer than the age of the earth",
    and past about ten thousand million years every barrier reads the same.
    Measured at 298.15 K: +25 kcal/mol is about 4 days, +30 about 50 years,
    and +45 is past anything worth printing as a figure.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -17.172582)

    assert "3.94 d" in part._thermal_wait(25.0, 298.15)
    assert "49.9 years" in part._thermal_wait(30.0, 298.15)
    assert "longer than the universe" in part._thermal_wait(45.0, 298.15)
    # And the fast end is where it was: an open path answers in picoseconds.
    assert "ps" in part._thermal_wait(0.0, 298.15)


def test_a_refusal_says_it_priced_the_way_the_hand_went():
    """Which is narrower than a refusal reads, and cannot be answered here.

    What was priced is the path this hand took, one geometry at a time; the
    cheapest way from the anchor to where the hand was aiming is a different
    quantity, and a minimum over all paths is a search rather than a drag.  So
    it is said once, where the refusal lands, and it names what does search.
    """
    assert "this prices the way " in EDITOR_SOURCE
    assert "not the cheapest way " in EDITOR_SOURCE
    # Named, so the sentence is a next move rather than a disclaimer -- and
    # named as the path finder rather than as the press it sits on, which
    # reads "Find the path" or "To the saddle" depending on the box beside it.
    assert "Scan and the path finder " in EDITOR_SOURCE
    assert "look for that'" in EDITOR_SOURCE
    # Only where the temperature refused.  A hold that slipped, or two atoms
    # inside each other, is not a barrier -- and telling the user to look for
    # a cheaper way round would be advice about something else.
    body = EDITOR_SOURCE.split("state.get('thermal_over') == 'path' else")[1]
    body = body.split("# And whether this was a slope or a step.")[0]
    assert "if not aside:" in body, body


def test_a_refusal_says_when_the_drag_changed_how_many_pieces_there_are():
    """The one case where the budget is wrong by more than the method is.

    A ceiling is a free energy and a drag is priced with an electronic one.
    While the structure stays in as many pieces as it started in, the two
    agree to under 3 kcal/mol -- ethane turned to eclipsed is +2.592 against
    +2.568 -- and where a drag takes something apart they part company by
    about ten: a borazane pulled to 6 A is +22.5 electronic and +12.3 free,
    which at 298 K is the difference between past the 22.3 and comfortably
    inside it.  Nothing is corrected, because a number invented off a distance
    threshold would be worse than the gap; but the refusal can say which case
    it is in, and that is the difference between a verdict and a verdict the
    user can weigh.
    """
    part, state = _editor(_NITROSAMINE)

    assert part._pieces_in(_NITROSAMINE) == 1
    assert part._pieces_in(_DIELS_ALDER) == 2
    # One atom taken right out is a piece of its own, which is the thing being
    # detected: what a drag does when it finally lets go of something.  Five
    # rather than two, because the carbon that was pulled took nothing with
    # it and left its three hydrogens stranded one apiece -- which is the
    # honest answer about that geometry, and why the sentence quotes the count
    # rather than saying "it has come apart".
    assert part._pieces_in(_pulled(_NITROSAMINE, 5.0)) == 5
    # And nothing at all is nothing, rather than one.
    assert part._pieces_in("") == 0

    # Counted where the anchor is taken, because a count read at the moment of
    # a refusal has nothing to be a change from.
    assert "state['thermal_pieces'] = _pieces_in(" in EDITOR_SOURCE
    assert "began_in = state.get('thermal_pieces')" in EDITOR_SOURCE
    assert "if began_in and now_in > began_in:" in EDITOR_SOURCE
    # And what it says: which way the price is wrong, and what gives the
    # right one.
    assert "pieces here '" in EDITOR_SOURCE
    assert "price is strict " in EDITOR_SOURCE
    assert "by about ten kcal/mol" in EDITOR_SOURCE
    assert "with its energy set to G prices it " in EDITOR_SOURCE


def test_the_window_the_ceiling_is_quoted_over_stays_an_hour():
    """A control for it was weighed again and not built.

    Between a second and a year the ceiling moves ten kcal/mol while the
    distance between chemistry and nonsense is twenty against a hundred, and
    what a refusal now says covers the ground the knob would have: the waiting
    time is quoted at the temperature that is set, so how much longer it would
    take is already on the line, and the temperature it wants is quoted for
    the window.  Maeda's advice to run the permissive end is advice about
    searching, where missing a path finds nothing; this decides what stays in
    the box, where allowing what the temperature will not is the worse error.
    """
    from delfin.dashboard.structure_editor import _THERMAL_SECONDS

    assert _THERMAL_SECONDS == 3600.0
    assert "an hour" in EDITOR_SOURCE.split("def _timescale_label(")[1][:400]


@_needs_xtb
def test_the_free_energy_parts_company_where_the_pieces_change():
    """Why the budget prices an electronic energy against a free ceiling, in
    numbers rather than as an assurance.

    A borazane with its B-N pulled out is the case where the two disagree
    most: one molecule becoming two releases translation and rotation, so the
    electronic price is too strict by whatever T*dS is.  Measured under GFN2 at
    298.15 K against the relaxed adduct: at 3.5 A it costs +21.2 kcal/mol
    electronic and +14.2 free; once the two are apart at 6 A, +22.5 against
    +12.3.  Ten and a quarter kcal/mol at the end of it -- and it is the one
    place a verdict changes hands, since the hour at 298 K buys 22.3.

    The other way, a drag that leaves the structure in as many pieces as it
    found it, they agree: ethane turned to eclipsed is +2.592 electronic and
    +2.568 free, 0.02 apart on a barrier of 2.6.
    """
    from delfin.dashboard.structure_editor import (
        _HARTREE_TO_KCAL, thermal_ceiling)

    adduct = """8
borazane
B   0.000000  0.000000  0.830000
N   0.000000  0.000000 -0.830000
H   1.100000  0.000000  1.230000
H  -0.550000 -0.953000  1.230000
H  -0.550000  0.953000  1.230000
H  -0.950000  0.000000 -1.180000
H   0.475000  0.823000 -1.180000
H   0.475000 -0.823000 -1.180000
"""
    at_rest = gfn.optimize_with_gfn(adduct, "gfn2", timeout=600)
    assert at_rest.get("ok"), at_rest.get("status")

    def priced(xyz):
        got = gfn.optimize_with_gfn(
            xyz, "gfn2", timeout=600, optimise=False, free_energy=True,
            thermo_kelvin=298.15)
        assert got.get("energy") is not None, got.get("status")
        assert got.get("free_energy") is not None, got.get("status")
        return float(got["energy"]), float(got["free_energy"])

    e0, g0 = priced(at_rest["xyz"])
    pulled = gfn.optimize_with_gfn(
        at_rest["xyz"], "gfn2", timeout=600,
        constraints=[{"kind": "distance", "atoms": [0, 1], "value": 6.0,
                      "mode": "fix"}])
    assert pulled.get("ok"), pulled.get("status")
    e1, g1 = priced(pulled["xyz"])

    electronic = (e1 - e0) * _HARTREE_TO_KCAL
    free = (g1 - g0) * _HARTREE_TO_KCAL
    # xtb is not bit-reproducible under threading, so this is asserted with
    # room: what is claimed is "of order ten kcal/mol and in this direction",
    # not a decimal.
    assert 20.0 < electronic < 25.0, electronic
    assert 9.0 < free < 16.0, free
    assert electronic - free > 6.0, (electronic, free)
    # And the verdict really does change hands there, which is the whole
    # reason the size of the gap is written where the user can read it.
    ceiling = thermal_ceiling(298.15, 3600.0)
    assert electronic > ceiling > free, (electronic, ceiling, free)


def test_an_answer_that_arrives_after_the_release_marks_nothing():
    """"warum wird die force nicht mehr entfernt ich seh immer noch die
    orangenen kugeln im viewer?"

    An xtb round outlives the gesture that asked for it. From the journal that
    came with the report::

        1425.260s  wall  {"reach": 0.0, "wall": null}   the release clears
        1425.264s  page  DELFIN drag-end
        1425.610s        Optimise starts
        1425.977s  said  "following the drag: 3 step(s)..."   <- late

    The third answer of a drag arrived after the release had already cleared
    the marks. Drawing them again there leaves them on the page for good: the
    release has run, and only the *next* grab clears them.
    """
    from editor_source import EDITOR_SOURCE

    body = EDITOR_SOURCE.split('def _stop_the_hand_at(', 1)[1].split(
        '\n    def ', 1)[0]
    before = body.split('marks = {}', 1)[0]
    assert "if not state.get('gfn_follow'):" in before, before
    assert before.rstrip().endswith('return'), (
        'the late answer leaves before it marks anything')


def test_a_stopped_relaxation_step_is_not_priced():
    """Letting go is the way out of a drag frame now, and a way out has a
    shape: a stopped run comes back with a geometry and no energy.

    The pricing took the arithmetic of that energy.  It was safe while the
    only way a step could end early was a clock -- then it either answered or
    the whole frame failed -- and it stopped being safe the moment the step
    could be *asked* to stop.

    Reported from the field, three times in one session under the budget:
    "the relaxation under the hand stopped on an error: TypeError: float()
    argument must be a string or a real number, not 'NoneType'".
    """
    step = EDITOR_SOURCE.split('def _gfn_follow_step')[1].split('\n    def ')[0]

    # The restraint-energy shortcut, which is the one that took it.
    priced = step.split("dict(outcome, energy=float(outcome['energy']) - bias)")[1]
    assert "outcome.get('energy') is not None" in priced[:220], priced[:220]

    # And every other place the budget reads a price is guarded too, because
    # any of them can now be handed a step that was stopped.
    for reader in ("float(priced['energy'])",
                   "state['thermal_priced'] = float(priced['energy'])"):
        at = step.index(reader)
        near = step[max(0, at - 400):at + 400]
        assert "priced.get('energy') is not None" in near, reader


def test_the_hand_follows_more_slowly_when_the_answers_are_slow():
    """A fixed share is a fixed share *per answer*, and what it damps does not
    grow per answer -- it grows with the time between them.

    The cursor runs on while xtb thinks, so the longer an answer takes the
    further the structure is behind when it lands and the larger the force
    that is asked for.  Half of a much larger demand is still large enough to
    overshoot; then the hand goes slack, and then it overshoots the other way.

    From the field, on a slow system: answers 1.7 s apart, the pull sitting on
    its ceiling every single time and the ceiling itself swinging between 59
    and 88 kcal/mol/A, and the dragged atom going 1.3 A out and back once an
    answer -- nine such there-and-back triples in eleven seconds.  All of them
    written "Within the budget", so the wall had nothing to do with it.

    None of that is the molecule.  A structure under a steady pull does not
    swing 1.3 A at 0.6 Hz; what was swinging was the force, worked out from a
    geometry 1.7 seconds old.
    """
    part, state = _editor(_NITROSAMINE)
    follows = part._hand_follows_now

    # Nothing measured yet, and at the rate it was measured at: as before.
    state.pop('gfn_follow_took', None)
    assert follows() == pytest.approx(0.5)
    state['gfn_follow_took'] = 0.10
    assert follows() == pytest.approx(0.5)

    # Twice as slow, half as eager -- the step per answer stays the same size
    # however far the cursor has run on in between.
    state['gfn_follow_took'] = 0.20
    assert follows() == pytest.approx(0.25)

    # And the session that was reported gets a twentieth instead of a half.
    state['gfn_follow_took'] = 1.70
    assert follows() == pytest.approx(0.05)

    # A floor, because scaled all the way down is not shaking but is not
    # dragging either.  An overshoot needs a step larger than the distance to
    # the target, and a twentieth cannot be that.
    state['gfn_follow_took'] = 30.0
    assert follows() == pytest.approx(0.05)

    # The number it works from is the number the user is shown.
    assert '_hand_answered_in(began)' in EDITOR_SOURCE
    assert "state['gfn_follow_took'] = took" in EDITOR_SOURCE


def test_the_hand_is_not_reset_by_the_coordinate_changing_its_name():
    """The damping was being skipped on every other answer.

    The hand is on a pair.  Which coordinate that pair is expressed in is
    worked out afresh from the geometry every answer, and it flips -- a
    contact that is a distance on one answer is a torsion on the next.  Keyed
    on both, every flip was a key nobody had seen, and a new key starts at
    what this answer asks for so that beginning a drag is immediate.

    Simulated with the demand alternating 44 and 90 and the kind flipping with
    it, the hand applied 44, 90, 44, 90 -- which is exactly the ceiling seen
    swinging between 44 and 91 in the field, once an answer, while the atom
    went 1.3 A out and back.
    """
    part, state = _editor(_NITROSAMINE)

    def drag(demand, kinds, atoms, took=1.7):
        state['gfn_follow_run'] = 1
        state['gfn_follow_took'] = took
        for key in ('gfn_hand_force', 'gfn_hand_force_run',
                    'gfn_hand_force_most'):
            state.pop(key, None)
        out = []
        for want, kind, pair in zip(demand, kinds, atoms):
            got = part._steady_hand(
                [{'kind': kind, 'atoms': pair, 'force': float(want)}])
            out.append(got[0]['force'])
        return out

    # The reported case: same pair, the name of the coordinate flipping.
    applied = drag([44, 90] * 4, ['distance', 'dihedral'] * 4, [[0, 1]] * 8)
    assert max(applied) < 60, applied          # was 90 every other answer
    assert all(b >= a - 0.5 for a, b in zip(applied, applied[1:])), applied

    # A pair that really is new starts at what is asked, and that is not a
    # bug to be damped out: a fragment walked past a molecule does pick up a
    # different contact, and a hand that would not begin a drag until it had
    # been asked three times is not a hand.  There was a ceiling here on how
    # much one answer may multiply the force by; it bought this case and cost
    # six answers to reach a sustained pull, and the shaking it was aimed at
    # turned out to be the budget's -- see _within_the_budget.
    applied = drag([44, 90] * 4, ['distance', 'dihedral'] * 4,
                   [[0, 1], [2, 3]] * 4)
    assert applied[0] == pytest.approx(44.0), applied

    # And a fast system is untouched: half the demand per answer, as measured.
    applied = drag(list(range(10, 90, 10)), ['distance'] * 8, [[0, 1]] * 8,
                   took=0.10)
    assert applied[1] == pytest.approx(15.0)


def _reach(part, xyz, holding):
    return part._how_far_the_hand_has_asked(xyz, holding)


def _pulled_out(text, index, far):
    """The same molecule with one atom taken *far* Angstrom along +x.

    What the page sends while a hand is down: everything as it stands, with
    the held atom where the cursor is.
    """
    rows = [line.split() for line in gfn.atom_lines(text)]
    rows[index][1] = f"{float(rows[index][1]) + float(far):.10f}"
    return (f"{len(rows)}\ncursor\n"
            + "\n".join(" ".join(one) for one in rows) + "\n")


def test_the_budget_is_spent_forwards_so_nothing_has_to_be_taken_back():
    """The shaking is the rollback, and the rollback is avoidable.

    The wall is a refusal after the fact, and it has to be -- nothing can
    price a geometry before it has been relaxed.  But a refusal after the fact
    is a rollback, and a rollback while the cursor is still held out at the
    forbidden place is a loop: the hand asks, the structure is put back, the
    hand asks again.  On screen that is a molecule shaking once an answer.

    So the ceiling is spent forwards, and what is cut back is the wish itself:
    the held atoms are drawn back along the way they have come, and every
    thing after that -- which contact is perceived, what force it becomes,
    what xtb is asked -- follows from a wish that is already affordable.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -25.0)
    ceiling = part._thermal_budget()[1]
    assert ceiling == pytest.approx(22.3, abs=0.5)
    state['thermal_hand_from'] = _NITROSAMINE
    held = [0]

    # How far the hand has asked: one number for the whole gesture, and the
    # only one the budget needs.  It cannot flip, whatever the perception
    # makes of the geometry.
    assert _reach(part, _NITROSAMINE, held) == pytest.approx(0.0)
    assert _reach(part, _pulled_out(_NITROSAMINE, 0, 0.4), held) == \
        pytest.approx(0.4)

    # Nothing measured yet: the wall is the backstop it always was.
    state['thermal_peak'] = 0.0
    far = _pulled_out(_NITROSAMINE, 0, 2.0)
    assert part._within_the_budget(far, held) is far

    # Two answers of a real drag: the hand asked for a fifth of an Angstrom
    # more and the answer cost four kcal/mol.  Twenty per Angstrom asked.
    part._note_what_it_costs(0.0, 0.2)
    part._note_what_it_costs(4.0, 0.4)
    assert state['thermal_cost'] == pytest.approx(20.0)

    # And now the hand may ask for what is left of the ceiling and no more.
    state['thermal_peak'] = 4.0
    allowed = ((ceiling - 4.0) * structure_editor._BUDGET_SPENDS_AT_MOST
               / state['thermal_cost'])
    kept = part._within_the_budget(far, held)
    assert _reach(part, kept, held) == pytest.approx(0.4 + allowed, abs=1e-6)
    assert state.get('thermal_held_back') is True
    # Only the held atom was moved; the rest of the wish is untouched.
    assert gfn.atom_lines(kept)[3] == gfn.atom_lines(far)[3]

    # A cursor inside what is left is not touched at all.  The budget cuts a
    # wish back; it does not drive one.
    state.pop('thermal_held_back', None)
    near = _pulled_out(_NITROSAMINE, 0, 0.4 + allowed / 2)
    assert part._within_the_budget(near, held) is near
    assert state.get('thermal_held_back') is None

    # Coming back is never blocked.  A budget that could be entered and not
    # left would be a trap, and what tells the two apart is whether the wish
    # asks for more than the last one did.
    back = _pulled_out(_NITROSAMINE, 0, 0.1)
    assert part._within_the_budget(back, held) is back

    # Spent to the last kcal/mol, the hand is held where it already asked:
    # no further travel, and -- because a push is sized by how far its target
    # has run ahead -- no more force either.
    state['thermal_peak'] = ceiling
    assert _reach(part, part._within_the_budget(far, held), held) == \
        pytest.approx(0.4, abs=1e-6)


def test_what_a_step_costs_is_measured_on_the_gesture_and_not_the_contact():
    """The second half of the fault the hand's damping had.

    Which pair the perception names the hand by is worked out afresh every
    answer, and between two fragments it does not settle.  Measured on the
    field's butadiene and ethene, dragged together: the coordinate went
    distance[10,0], [10,5], [10,4], [10,5], [10,4], [10,0] over six answers.
    A rate kept per coordinate is dropped on every one of those changes, so
    the forward spend never engaged at all -- the user, on exactly this:
    "zappeln ist aber noch vor allem bei systemen die nicht verbunden sind
    sprich zwei molekuelen".

    So the rate is a price per Angstrom the *hand* has asked for.  That
    cannot flip, it is monotone with the gesture, and it is defined for a
    bonded coordinate and two unbonded fragments alike.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -25.0)

    # A hand that has barely asked for anything more says nothing about what
    # asking costs: the answer would be the noise of two relaxations divided
    # by about zero, and this number decides how far a hand may go.
    part._note_what_it_costs(0.0, 0.20)
    part._note_what_it_costs(0.3, 0.2005)
    assert state.get('thermal_cost') is None

    # Steeper is taken whole, on the answer it steepens: what this is for is
    # not overspending, and every approach to a wall stiffens.
    part._note_what_it_costs(0.0, 0.2)
    part._note_what_it_costs(2.0, 0.3)
    assert state['thermal_cost'] == pytest.approx(20.0)
    part._note_what_it_costs(8.0, 0.4)
    assert state['thermal_cost'] == pytest.approx(60.0)
    assert state['thermal_cost_grew'] == pytest.approx(2.0)

    # Flatter, half at a time: half of sixty and half of the thirty this one
    # measured.
    part._note_what_it_costs(11.0, 0.5)
    assert state['thermal_cost'] == pytest.approx(45.0)

    # Asking for more that costs nothing, or less than nothing, is the
    # structure settling rather than the surface, and there is nothing to
    # spend against.
    was = state['thermal_cost']
    part._note_what_it_costs(11.0, 0.6)
    assert state['thermal_cost'] == pytest.approx(was)
    part._note_what_it_costs(2.0, 0.7)
    assert state['thermal_cost'] == pytest.approx(was)

    # And a new grab starts without any of it.
    clears = EDITOR_SOURCE.split('def _clear_thermal_wall(')[1].split(
        '\n    def ')[0]
    for key in ('thermal_cost', 'thermal_cost_at', 'thermal_hand_from'):
        assert f"state.pop('{key}', None)" in clears, key
    begins = EDITOR_SOURCE.split('def _begin_gfn_follow(')[1].split(
        '\n    def ')[0]
    assert "state['thermal_hand_from'] = state['gfn_topology_source']" in begins


def test_a_hand_held_back_by_the_budget_says_so_and_shows_where():
    """A drag that stops following looks exactly like a drag that has broken.

    With the ceiling spent forwards the wall hardly ever fires, so the two
    things that used to tell the user a budget was refusing -- the sentence
    and the mark on the picture -- both went quiet, and what was left was a
    structure that had simply stopped moving under the cursor.  That is the
    same silence this whole wall exists to close.

    So the clause is on the line the drag is already writing, and the mark
    goes on the picture at the geometry that was reached, which is where the
    budget says the atom may stand.
    """
    body = EDITOR_SOURCE.split('def _gfn_follow_step(')[1].split(
        '\n    def ')[0]
    # The clause, on the same line rather than under it: the row stands above
    # the viewer and a second one moves the picture while an atom is aimed.
    assert "held at what the budget allows" in body
    assert "state.get('thermal_held_back')" in body
    # And the mark, through the same path a rollback marks by.
    assert 'else (reached if' in body
    assert '_stop_the_hand_at(' in body
    # Raised by the clamp and cleared before it, so it is about this answer.
    clamp = EDITOR_SOURCE.split('def _within_the_budget(')[1].split(
        '\n    def ')[0]
    assert "state['thermal_held_back'] = True" in clamp
    assert "state.pop('thermal_held_back', None)" in body


def test_a_hand_resting_on_the_ceiling_keeps_its_grip():
    """It rests against the ceiling; it does not let go of the structure.

    A wish drawn back to where the hand last asked is still a wish a long way
    ahead of the atoms, and it is that gap which is holding them out.  Drawn
    back to where the *atoms* are instead, the next answer would ask for the
    place they already occupy, which is no force at all: measured live, an
    ethane held at its ceiling relaxed the whole way home.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -25.0)
    ceiling = part._thermal_budget()[1]
    state['thermal_hand_from'] = _NITROSAMINE
    held = [0]

    part._note_what_it_costs(0.0, 0.8)
    part._note_what_it_costs(10.0, 1.0)
    state['thermal_peak'] = ceiling
    kept = part._within_the_budget(_pulled_out(_NITROSAMINE, 0, 3.0), held)
    # The reach it already had, not the nothing the atoms have moved.
    assert _reach(part, kept, held) == pytest.approx(1.0, abs=1e-6)
    assert _reach(part, kept, held) > 0.5

    # And the wall records it, so a rollback goes back to asking for what was
    # asked when the budget last agreed -- not for what overshot.
    wall = EDITOR_SOURCE.split('def _thermal_wall(')[1].split('\n    def ')[0]
    assert "state['thermal_good_ask'] = (" in wall
    assert "asked=float(held)" in wall


def test_the_allowance_is_worked_out_at_the_rate_the_next_answer_will_have():
    """A rate is a secant behind; the clamp spends it on the answer ahead.

    On any approach to a wall the surface is stiffening, so the secant behind
    is smaller than the tangent ahead and the allowance is a step too
    generous.  Measured live on an ethane under GFN2: the rate read 35.4
    kcal/mol per Angstrom asked, the answer realised 47, and the wall -- which
    is meant to be the backstop -- had to fire and take the step back.

    So the growth of the rate is carried and the allowance is divided by the
    rate the next answer is expected to have.  Bounded at a doubling: a growth
    factor read off two noisy secants can be anything, and no real stiffening
    doubles inside one answer of a drag.
    """
    part, state = _editor(_NITROSAMINE)
    _armed(part, state, _NITROSAMINE, -25.0)
    ceiling = part._thermal_budget()[1]
    state['thermal_hand_from'] = _NITROSAMINE
    held = [0]

    part._note_what_it_costs(0.0, 0.1)
    part._note_what_it_costs(2.0, 0.2)
    assert state['thermal_cost_grew'] == pytest.approx(1.0)
    part._note_what_it_costs(5.0, 0.3)
    assert state['thermal_cost'] == pytest.approx(30.0)
    assert state['thermal_cost_grew'] == pytest.approx(1.5)

    state['thermal_peak'] = 10.0
    left = (ceiling - 10.0) * structure_editor._BUDGET_SPENDS_AT_MOST
    kept = part._within_the_budget(_pulled_out(_NITROSAMINE, 0, 5.0), held)
    assert _reach(part, kept, held) == pytest.approx(
        0.3 + left / (30.0 * 1.5), abs=1e-6)


