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
