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

import pathlib

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
    assert "if contacts and pull is None:\n                        priced = outcome" in EDITOR_SOURCE
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
    assert "bottom = (spent, walked, reached)" in body
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
    assert body.index("wants = (") < body.index("if rise <= ceiling:")
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
    assert "else current)[0]" in body


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
    # And describing a leg survives atoms that are gone.
    body = EDITOR_SOURCE.split("def _describe_leg(")[1].split("\n    def ")[0]
    assert "if 0 <= index < len(known)" in body


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
    body = EDITOR_SOURCE.split("submit_thermal_btn.layout.display = '' if xtb")[1]
    body = body.split("\n    def ")[0]
    assert "submit_scan_btn.layout.display = '' if xtb else 'none'" in body
    assert "submit_scan_run_btn" in body


def test_holding_a_value_is_a_step_of_its_own():
    """Judged on the picture alone a Hold was never one -- it changes no
    coordinate -- so Undo walked straight past it, wiped it on the way, and
    reported whatever it did land on.  Hold, then a scan, then Undo took back
    two actions on one press and named one of them."""
    body = EDITOR_SOURCE.split("def _remember(what)")[1].split("\n    def ")[0]
    assert "history[-1].get('constraints')" in body
    assert "list(state.get('constraints') or [])" in body
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
    assert "pricing = bool(submit_thermal_btn.value) and not stale" in source
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
    # not dictate a value -- that is the whole point of it.
    assert 'reached = (_value_in(walked, legs[0]) if pushing' in source
    assert 'path.append((reached, spent))' in source


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
    assert "path.append((_value_in(walked, legs[0]), 0.0))" in source


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
    """
    source = EDITOR_SOURCE
    assert "options=[('closer together', 'in'), ('further apart', 'out')]" in source
    assert "else [('narrower', 'in'), ('wider', 'out')])" in source
    assert "kind = _CONSTRAINT_KINDS.get(picked)" in source
    assert "description='steps'" in source
    assert "description='to'" in source

    assert 'submit_scan_stop_at = widgets.Checkbox(' in source
    assert "description='to a set value'" in source
    assert "set_end = wanted == '' and bool(submit_scan_stop_at.value)" in source
    assert "submit_scan_to.layout.display = '' if set_end else 'none'" in source
    # The direction is what is always used; a value only overrides it when one
    # was actually asked for.
    assert ('target = _suggest_scan_target(kind, here, submit_scan_way.value)'
            in source)
    assert 'if submit_scan_stop_at.value:' in source
    # Opened on the value the coordinate has, not on a zero to be guessed at.
    assert 'submit_scan_to.value = float(submit_internal_value.value)' in source
    assert 'submit_scan_stop_at.observe(on_submit_scan_stop_at' in source


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
    # rather than a wrong number being reported.
    assert 'if bias is not None else' in source

    # A held value and a pull cannot share one force constant, so the hold
    # stands and the hand goes back to placing -- said, not done quietly.
    assert 'held_too = bool(constraints) and pull is not None' in source
    assert 'force constant in xtb, so the hand is ' in source
