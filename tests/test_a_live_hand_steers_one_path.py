"""The live hand walks one continuous path, which is the whole point of it.

The other drag engine pins the coordinate the hand has changed and re-minimises
the rest around it.  At a reaction top -- a bifurcation, a transition ridge --
that biased minimum has two solutions, and the minimiser flips between them
answer to answer: the two arrangements a molecule cannot be held on, reported
on a real cis/trans study as a bond that would not settle
(``0902-174939``).  More cycles flip more cleanly; it is inherent to
minimising a bifurcated biased surface.

:func:`climb.steer` does not minimise.  It springs each dragged atom towards
where the cursor has it, adds that to the field's own force, and moves every
atom down the total by a step capped per atom -- the quasi-static limit of
interactive molecular dynamics.  Because it carries the geometry forward from
the last answer it stays on one branch until the force pushes it over, so a
torsion driven by the hand advances one way rather than jumping back and forth.

Measured here through the same call the editor makes: a butane's central
torsion, driven by leading its terminal carbon tangentially, sweeps smoothly
and never reverses.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from delfin.dashboard import climb as C
from delfin.dashboard import gfn_optimize as gfn


_needs_xtb = pytest.mark.skipif(
    not C.have_fast_gradients() and gfn.find_xtb() is None,
    reason='no xtb to take gradients from')


# A relaxed n-butane (GFN2), the central C0-C1-C2-C3 torsion at anti.
_BUTANE = """14
butane, relaxed under GFN2
C  -1.54258393  -0.35536195   0.40104530
C  -0.70744036   0.40824399  -0.61910669
C   0.72604655  -0.11291589  -0.68783670
C   1.56472464   0.65918240  -1.69865454
H  -1.11190925  -0.25935801   1.39569002
H  -1.58190998  -1.41271644   0.14722368
H  -2.55983132   0.02872430   0.42912843
H  -0.69261017   1.46842827  -0.35515790
H  -1.16943788   0.31997548  -1.60535233
H   0.71098686  -1.17069122  -0.96132958
H   1.18541349  -0.03392440   0.30040534
H   1.61247655   1.71252240  -1.42993563
H   1.13116952   0.58052009  -2.69354809
H   2.57900527   0.26797097  -1.73457132
"""

_CARBONS = (0, 1, 2, 3)


def _coords(xyz):
    return np.array(gfn.coordinates_of(xyz), dtype=float).reshape(-1, 3)


def _dihedral(P, a, b, c, d):
    b0, b1, b2 = P[a] - P[b], P[c] - P[b], P[d] - P[c]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _wish(P, held, lead):
    """Where the cursor leads the held atom: tangent to its turn about C1-C2."""
    a, b, c, _d = _CARBONS
    axis = P[c] - P[b]
    axis = axis / np.linalg.norm(axis)
    arm = P[held] - P[c]
    arm = arm - np.dot(arm, axis) * axis
    tang = np.cross(axis, arm)
    tang = tang / (np.linalg.norm(tang) + 1e-9)
    out = P.copy()
    out[held] = P[held] + tang * lead
    return out


def _as_xyz(P):
    rows = [f'{s} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}'
            for s, r in zip('CCCCHHHHHHHHHH', P)]
    return f'{len(rows)}\nwish\n' + '\n'.join(rows) + '\n'


@_needs_xtb
def test_a_driven_torsion_advances_and_never_reverses():
    held = 3          # a terminal carbon, dragged around the central bond
    current = _BUTANE
    trail = [_dihedral(_coords(current), *_CARBONS)]
    for _ in range(18):
        P = _coords(current)
        wish = _as_xyz(_wish(P, held, 0.30))
        out = C.steer(current, wish, [held], method='gfn2', cores=4)
        assert out.get('ok'), out
        current = out['xyz']
        trail.append(_dihedral(_coords(current), *_CARBONS))

    # The hand moved the torsion a real distance -- this is a drag, not a
    # structure standing still.
    swept = trail[-1] - trail[0]
    assert abs(swept) > 12.0, trail

    # And it moved it one way.  Every step over the noise floor has the same
    # sign as the sweep: no answer walks back the way the one before it came,
    # which is exactly the alternation this engine exists to remove.
    steps = [trail[i] - trail[i - 1] for i in range(1, len(trail))]
    real = [d for d in steps if abs(d) > 0.5]
    assert real, trail
    forward = 1 if swept > 0 else -1
    assert all((1 if d > 0 else -1) == forward for d in real), trail


@_needs_xtb
def test_the_answer_carries_the_geometry_forward():
    """Continuous, not restarted: each answer is a short walk from the last,
    so no single atom jumps.  The cap is 0.12 A per atom over a few steps, so
    a whole answer moves any atom well under an angstrom."""
    held = 3
    P0 = _coords(_BUTANE)
    wish = _as_xyz(_wish(P0, held, 0.30))
    out = C.steer(_BUTANE, wish, [held], method='gfn2', cores=4)
    assert out.get('ok'), out
    moved = np.linalg.norm(_coords(out['xyz']) - P0, axis=1).max()
    assert moved < C.STEER_CAP * C.STEER_STEPS + 0.05, moved
    # And the returned structure is priced: an energy of the geometry handed
    # back, not of the one before its last move.
    assert out.get('energy') is not None


_ETHANE = """8
ethane, relaxed under GFN2
C   0.000000  0.000000  0.765000
C   0.000000  0.000000 -0.765000
H   0.000000  1.019000  1.163000
H  -0.882000 -0.510000  1.163000
H   0.882000 -0.510000  1.163000
H   0.000000 -1.019000 -1.163000
H   0.882000  0.510000 -1.163000
H  -0.882000  0.510000 -1.163000
"""


#: A relaxed cis-2-butene (GFN2), carbons 0-1-2-3, the C=C torsion near cis.
_CIS_BUTENE = """12
cis-2-butene, relaxed under GFN2
C   1.56133168  -0.30153550   0.34941920
C   0.58962547   0.81664449   0.15552991
C  -0.69758682   0.72956008  -0.14104061
C  -1.50476993  -0.50896586  -0.35699709
H   2.38089185  -0.20604371  -0.36277491
H   1.10146901  -1.27583353   0.21959173
H   1.98910506  -0.24943338   1.35057641
H   1.02420044   1.80112840   0.27987105
H  -1.26489964   1.64626528  -0.24753152
H  -1.93691854  -0.49814727  -1.35756982
H  -0.91382776  -1.41217448  -0.24475249
H  -2.32872082  -0.54156452   0.35577813
"""


def _relaxed(xyz):
    return gfn.optimize_with_gfn(xyz, 'gfn2', optimise=True, max_steps=200,
                                 timeout=120)['xyz']


def _drag_a_hydrogen_out(share, *, etemp=None, answers=26):
    """Drag H2 of an ethane straight out under a hand of *share* of a bond,
    and return the C0-H2 distance it reaches."""
    cap = share * gfn.A_BOND_HOLDS * C.FORCE_KCAL_PER_A_IN_AU
    current = gfn.optimize_with_gfn(_ETHANE, 'gfn2', optimise=True,
                                    max_steps=200, timeout=120)['xyz']
    for k in range(answers):
        P = _coords(current)
        wish = P.copy()
        wish[2] = P[2] + np.array([0.0, 0.0, 0.15 * (k + 1)])
        rows = [f'{s} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}'
                for s, r in zip('CCHHHHHH', wish)]
        wish_xyz = f'8\nw\n' + '\n'.join(rows) + '\n'
        out = C.steer(current, wish_xyz, [2], method='gfn2', cores=4,
                      max_force=cap, etemp=etemp)
        assert out.get('ok'), out
        current = out['xyz']
    P = _coords(current)
    return float(np.linalg.norm(P[0] - P[2]))


@_needs_xtb
def test_a_weak_hand_lags_and_a_strong_hand_tears():
    """The whole of the user's question: the live hand is a force with a
    ceiling, the pull's, so a gentle drag cannot tear a bond off however far
    the cursor runs ahead -- the atom lags instead -- and only a hand set as
    strong as a bond breaks one, deliberately.

    Measured here: an ethane hydrogen dragged 3.9 A out along its bond.  At a
    hand of 0.4 of a bond -- the default, what room temperature allows -- the
    C-H holds near its length; at 2.0 of a bond it goes.  An uncapped spring
    tears it at any strength, which is what this ceiling is here to stop.
    """
    held = _drag_a_hydrogen_out(0.4)
    assert held < 1.5, f'a gentle hand tore the C-H (reached {held:.2f} A)'

    # As strong as two bonds it breaks -- smeared, because a bond coming apart
    # closes the frontier gap and the closed-shell SCC needs the warmth, which
    # is exactly what the follow loop gives it.
    torn = _drag_a_hydrogen_out(2.0, etemp=1000.0)
    assert torn > 1.9, f'a hand as strong as two bonds did not tear it ({torn:.2f} A)'


@_needs_xtb
def test_the_drive_hand_forces_a_torsion_over_its_barrier():
    """What the live hand cannot do and this is for: force a chosen coordinate
    over its barrier.

    A Cartesian pull on an atom takes the softest way to move it, so a stiff
    torsion bends the whole molecule rather than twisting -- cis stays cis
    however hard it is dragged.  :func:`climb.steer_coordinate` restrains the
    coordinate itself, so ramping its target walks it across the top.  Measured
    on a cis-butene: the C=C torsion driven from cis towards trans crosses 90
    degrees -- the twisted top of a 65 kcal/mol barrier -- monotonically, and
    the C=C stays a bond.

    Nothing here is about butene.  The engine reads an energy and a gradient,
    which every SCC method gives, and drives whatever coordinate it is handed;
    the two tests below drive a distance and an angle with the same call.
    """
    cis = _CIS_BUTENE
    torsion = (0, 1, 2, 3)
    current = cis
    trail = [_dihedral(_coords(current), *torsion)]
    for target in range(10, 181, 10):
        out = C.steer_coordinate(current, 'dihedral', list(torsion),
                                 float(target), method='gfn2', cores=4,
                                 etemp=1000.0, steps=6)
        assert out.get('ok'), out
        current = out['xyz']
        trail.append(_dihedral(_coords(current), *torsion))

    crossed = max(abs(d) for d in trail)
    assert crossed > 100.0, f'the drive never crossed the barrier top: {trail}'
    # Monotone towards trans, no walking back -- the continuity the pin lacked.
    steps = [abs(trail[i]) - abs(trail[i - 1]) for i in range(1, len(trail))]
    assert sum(1 for s in steps if s < -3.0) == 0, f'it walked back: {trail}'
    P = _coords(current)
    assert np.linalg.norm(P[1] - P[2]) < 1.7, 'the C=C came apart'


@_needs_xtb
def test_the_drive_hand_is_universal_across_coordinates():
    """The same engine drives a distance and an angle, no per-system code."""
    # A bond, stretched.
    out = C.steer_coordinate(_ETHANE, 'distance', [0, 1], 2.2, method='gfn2',
                             cores=4, etemp=1000.0, steps=10)
    assert out.get('ok'), out
    P = _coords(out['xyz'])
    assert np.linalg.norm(P[0] - P[1]) > 1.9, 'the bond did not stretch'

    # An angle, opened.
    water = _relaxed('3\nwater\nO 0 0 0\nH 0 0.757 0.587\nH 0 -0.757 0.587\n')
    out = C.steer_coordinate(water, 'angle', [1, 0, 2], 135.0, method='gfn2',
                             cores=4, steps=10)
    assert out.get('ok'), out
    P = _coords(out['xyz'])
    u = P[1] - P[0]
    v = P[2] - P[0]
    ang = math.degrees(math.acos(
        float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))))
    assert ang > 120.0, f'the angle did not open: {ang:.0f}'


@_needs_xtb
def test_a_still_hand_barely_moves_the_structure():
    """The wish at the atom it already holds is no force at all, so the answer
    stays where it started -- a hand held still leaves a still molecule."""
    P0 = _coords(_BUTANE)
    wish = _as_xyz(P0)          # cursor exactly on the atom
    out = C.steer(_BUTANE, wish, [3], method='gfn2', cores=4)
    assert out.get('ok'), out
    moved = np.linalg.norm(_coords(out['xyz']) - P0, axis=1).max()
    assert moved < 0.05, moved
