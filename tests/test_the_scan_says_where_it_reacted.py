"""A scan is asked "how high", but what it is really asked is "where did it
react" -- and the two are not always the same place.

Measured in :func:`delfin.dashboard.gfn_optimize.paths_disagree` on the
Diels-Alder: the walk out puts its summit at 2.20 A with the forming bond -- the
one nobody is driving -- still at 2.92, and the walk back puts it at 2.90 with
that bond at 1.76.  Two maxima 0.7 A apart on the coordinate that was driven and
1.2 A apart on the one that was not.

So the walk carries its bond graph along.  It costs nothing: the geometries are
already in hand and the graph is covalent radii, the same test the picture draws
its lines with.  Where the bonding changes is where the molecule stopped being
the reactant, and whether the summit sits there decides whether the height above
it is a barrier or a deformation.

And the barrier is said as the bound it is.  A relaxed scan drives a coordinate
somebody chose, so its highest point is the top of the best path that keeps to
that coordinate; the saddle is the lowest such top over every path there is, and
lies at or below it.  Written as a plain number it reads as a measurement that a
saddle search would agree with.  It will not: it will be lower.
"""

from __future__ import annotations

import math

import pytest

from delfin.dashboard import gfn_optimize as gfn

from editor_source import EDITOR_SOURCE as SOURCE


def _verdict():
    return SOURCE.split('def _scan_verdict')[1].split('\n    def ')[0]


def test_the_barrier_is_said_as_a_bound():
    verdict = _verdict()
    assert "at_most = '\\u2264' if rise > 0 else ''" in verdict
    # On both sentences: the electronic one and the free-energy one.  Said on
    # one and not the other, a reader comparing two runs would think the mode
    # changed the meaning of the number.
    assert verdict.count('{at_most}') == 2
    # And not on a walk that found no rise, where "at most 0.0" is a bound on
    # nothing.
    assert 'if rise > 0' in verdict


def test_the_walk_says_where_the_bonding_changed():
    walk = SOURCE.split('def on_submit_scan_run')[1].split('\n    def ')[0]
    assert "state['scan_became'] = _where_the_bonding_changed(path)" in walk

    body = SOURCE.split('def _where_the_bonding_changed')[1] \
                 .split('\n            def ')[0]
    # The geometries the walk already has, and the graph the picture draws.
    assert '_gfn.bond_graph(' in body
    assert '_gfn.graph_changed(' in body
    # The first change only.  A walk that carries on past its reaction changes
    # the bonding again climbing out of the well, and a list of every change
    # is a list nobody reads.
    assert 'return {' in body
    assert body.count('return {') == 1
    # A point without a geometry is skipped rather than crashed on.
    assert 'if not geometry:' in body


def test_a_top_away_from_the_reaction_is_said_to_be_one():
    """The case worth the words.  A summit sitting on the bonding change is
    the ordinary one and needs no announcement; a summit a long way from it is
    a deformation being reported as a barrier, and the saddle search is what
    answers it."""
    verdict = _verdict()
    assert "became = state.get('scan_became')" in verdict
    assert 'at the top' in verdict
    assert 'a deformation, not this reaction' in verdict
    # Judged against the spacing of the walk rather than against a fixed
    # distance: a scan of a torsion moves in degrees and one of a bond in
    # Angstrom, and one number cannot be near in both.
    assert 'spacing' in verdict
    assert '1.5 * spacing' in verdict


def test_a_push_reaches_as_far_in_degrees_as_it_does_in_angstrom():
    """A reach in the coordinate's own units, and it was not.

    The scan's push holds its target a reach ahead of the structure and ramps
    the force, so what the molecule feels is a force and not a spring that
    pulls harder the further it has to go.  The reach was an Angstrom whatever
    the coordinate was -- so on an angle or a torsion, whose values are in
    degrees, the target sat one degree ahead where a reach is a hundred and
    eight of them.

    A push that asks for a hundredth of what it means to ask for does not
    refuse.  It ramps a force against a target the structure is already at,
    nothing happens, and the user meets it as "on an angle only walk the value
    works" -- which is what was reported.  gfn_optimize.PUSH_REACH_DEGREES
    existed the whole time and the drag had been converting with it; this had
    not.

    And the floor goes with it: a bond cannot be asked for at or below zero,
    but 0.2 degrees is not a floor on an angle, and a torsion is periodic and
    has no floor at all.
    """
    push = SOURCE.split('def _push_target(')[1].split('\n        def ')[0]
    assert '_gfn.PUSH_REACH_DEGREES' in push, push[:400]
    assert "kind == 'distance'" in push
    # The floor is a distance's, and is applied as one.
    assert "max(0.2, target) if kind == 'distance' else target" in push

    # And the two reaches really are the same reach, said twice.
    assert gfn.PUSH_REACH_DEGREES == pytest.approx(
        math.degrees(gfn.PUSH_REACH / gfn.BOHR_IN_ANGSTROM))
    assert gfn.PUSH_REACH_DEGREES > 100.0, 'a hundredth of this was the bug'
