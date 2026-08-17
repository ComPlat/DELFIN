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
    assert held[0]["mode"] == "fix"
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


def test_a_structure_moved_as_one_has_nothing_to_hold():
    """Sliding everything sideways changes no internal coordinate.

    It costs nothing and there is nothing to hold it against, so the caller is
    handed no values back and prices that drag some other way.
    """
    assert gfn.contacts_holding(_ETHANE, range(8), most=2) == []
    assert gfn.contacts_holding(_ETHANE, [], most=2) == []


def test_the_held_values_reach_xtb_as_exact_ones():
    """A held contact is not a suggestion: it is what the price is about."""
    text = gfn.constraint_input(
        gfn.contacts_holding(_ETHANE, [2], most=1), atoms=8)["text"]
    assert "$constrain" in text
    assert f"force constant={gfn.FIX_FORCE_CONSTANT}" in text
    # xtb counts from one.
    assert "distance: 3, 1," in text


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


def test_more_cores_are_not_more_speed():
    """The ladder is by size and it ends, because xtb's gain does.

    Measured on this machine at twenty GFN2 cycles: eighteen atoms are fastest
    on four cores, a hundred and eighty-five on sixteen, and sixty-four cores
    are three times worse than four.  A machine with 384 of them does not
    change that.
    """
    assert gfn.interactive_cores(18) == 4
    assert gfn.interactive_cores(49) == 8
    assert gfn.interactive_cores(101) == 16
    assert gfn.interactive_cores(5000) == 16
    # No size named: an interactive drag is mostly a small structure.
    assert gfn.interactive_cores() == gfn.INTERACTIVE_CORES


def test_the_cores_can_still_be_said_outright(monkeypatch):
    """Whoever sets the variable has said what they want."""
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
