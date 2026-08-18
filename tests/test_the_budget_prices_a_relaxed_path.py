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


def test_without_a_before_there_is_nothing_to_compare():
    """The first answer of a drag has no previous geometry, and says so."""
    held = gfn.contacts_holding(_BUTANE, [3], most=2)
    assert all(one["kind"] == "distance" for one in held)


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
    assert "wants about" in body
    # Wrapped across two lines in the source, so only the start of it.
    assert "the whole path is " in body


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
    assert "if bottom is None or spent < bottom[0]:" in body
    assert "bottom = (spent, walked, held[0]['value'])" in body
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
    # One phrase, before the two branches, so both carry it.
    assert body.index("wants = (") < body.index("if top[1] <= ceiling:")
    assert "the whole path is open" in body.replace("'\n", "").replace(
        "                    f'", "")


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
    """
    follow = EDITOR_SOURCE.split("def _thermal_wall(")[1].split("def ")[0]
    assert "def _thermal_wall" not in follow
    assert "if spent <= ceiling and not refuse:" in follow
    body = EDITOR_SOURCE.split("state['thermal_now'] = priced.get('energy')")[1]
    body = body.split("state['gfn_last_status'] = said")[0]
    assert "refuse=(slipped > _SLIP_LOOSE) or crowded" in body
    assert "_gfn.closest_contact(current)" in body


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
    assert "way == 'out'" in body
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

