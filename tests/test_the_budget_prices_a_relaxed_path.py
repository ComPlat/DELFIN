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

import shutil

import pytest

from delfin.dashboard import gfn_optimize as gfn
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


def test_two_contacts_are_held_when_both_are_the_same_kind():
    """Two bonds forming at once are both held: that is a cycloaddition."""
    # Two carbons of one fragment, each with its own carbon partner opposite.
    pair = (
        "4\ntwo approaching pairs\n"
        "C  0.00 0.00 0.00\n"
        "C  1.40 0.00 0.00\n"
        "C  0.00 0.00 2.20\n"
        "C  1.40 0.00 2.20\n"
    )
    held = gfn.contacts_holding(pair, [2, 3], most=2)
    assert len(held) == 2
    assert {frozenset(h["atoms"]) for h in held} == {
        frozenset({2, 0}), frozenset({3, 1})}
    for one in held:
        assert one["value"] == pytest.approx(2.20, abs=0.01)


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


def test_one_atom_alone_is_held_by_two_partners():
    """One distance does not hold an atom: it slides along the other bond.

    A ring carbon pulled out of a benzene, asked to stand at 1.72, 1.95 and
    2.20 A, came back at 1.41 every time however stiffly it was held -- the
    ring closed behind it and the price was that of a squeezed but whole
    benzene, +22 kcal/mol, which room temperature very nearly affords.  Held
    by two partners the same three pulls cost +36.6, +75.0 and +105.7.
    """
    ring = (
        "5\na carbon between two others\n"
        "C  0.00 0.00 0.00\n"
        "C  1.50 0.00 0.00\n"
        "C -0.75 1.30 0.00\n"
        "H  0.00 0.00 1.09\n"
        "H  2.00 0.90 0.00\n"
    )
    held = gfn.contacts_holding(ring, [0], most=2)
    assert len(held) == 2, "one contact lets it slide back"
    assert {h["atoms"][1] for h in held} == {1, 2}


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
    held = gfn.contacts_holding(across, [3], most=2, was=_BUTANE)
    assert [one["kind"] for one in held] == ["dihedral"]


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
    held = gfn.contacts_holding(closer, [3], most=2, was=apart)
    assert held and held[0]["kind"] == "distance"
    # The carbon it is being driven at, not the carbon it is bonded to.
    assert set(held[0]["atoms"]) == {3, 0}


def test_without_a_before_there_is_nothing_to_compare():
    """The first answer of a drag has no previous geometry, and says so."""
    held = gfn.contacts_holding(_BUTANE, [3], most=2)
    assert all(one["kind"] == "distance" for one in held)


def test_a_turn_is_not_put_back_where_the_cursor_had_it():
    """The torsion is what the hand asked for; a position would fight it.

    A methyl snapped back cannot spin about its own axis to get out of the
    way, and the same turn that costs +7.9 kcal/mol along the torsion costs
    +345 held rigid.
    """
    assert "str(one.get('kind')) == 'dihedral'" in EDITOR_SOURCE
    assert "settled = (outcome['xyz'] if turning else" in EDITOR_SOURCE
    # And the geometry a turn is measured against is the one this answer
    # handed back, not the one it was handed: against the latter the
    # difference holds the relaxation as well as the hand, and on a ring
    # being puckered the relaxation is the larger of the two.
    assert "state['thermal_was'] = settled" in EDITOR_SOURCE


def test_the_follow_prices_the_relaxation_it_ran():
    """One calculation, not two: the follow and the price are the same run."""
    assert "priced = outcome if contacts else (" in EDITOR_SOURCE
    assert "_gfn.contacts_holding(" in EDITOR_SOURCE
    # And the held contacts go into the run that follows the drag, or the
    # relaxation would still be free to undo it.
    assert "constraints=constraints + contacts" in EDITOR_SOURCE


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
