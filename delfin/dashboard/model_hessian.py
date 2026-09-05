"""A Hessian guessed from where the atoms are, costing no calculation at all.

A saddle search has to start from *some* second-derivative matrix, and
:meth:`delfin.dashboard.climb.Climb.numerical_hessian` buys one by differencing
the gradient: ``6N`` of them, which is 72 for butane and 330 for a manganese
complex -- thirteen minutes of waiting on this box before the climb takes its
first step.

Lindh's model (Lindh, Bernhardsson, Karlstrom, Malmqvist, *Chem. Phys. Lett.*
**241** (1995) 423) buys one from the geometry instead.  Every pair, triple and
quadruple of atoms contributes a stretch, a bend and a torsion, each weighted
by how close the atoms are::

    rho_ij = exp(alpha_ij (r_ij^ref^2 - r_ij^2))

and the matrix is the sum of ``k b b^T`` over those internals, ``b`` being an
internal's derivative with respect to the Cartesian coordinates.  Nothing has
to be decided about which atoms are bonded: the weight decides it, smoothly,
which is why no bond perception is involved and why a structure the hand has
pulled halfway through a reaction is not a special case.

**What it is right about, measured against differenced Hessians under GFN2
on twelve structures from 14 to 74 atoms.**  The soft *space*, on everything
rigid or moderately flexible: the model's softest mode sits inside the true
softest five at 0.96 to 0.99 for drags, a hydrogen-bonded dipeptide, a floppy
heptane, two molecules held apart, a 37-atom anion, a macrocycle, sucrose,
cholesterol, and a 57-atom manganese complex with two elements past the table
below.  Along a push it holds: displaced until the furthest atom has moved 1 A
and the energy is 20 kcal/mol up, the anion's model mode still sits 0.883 in
the true softest three.

**What it is wrong about.** The *magnitudes*: its curvatures come out four to
six times too stiff.  The *ordering* within the soft modes, which for
near-degenerate torsions is arbitrary -- heptane's model mode is 0.001 in the
true softest and 1.000 in the softest three, having picked a different member
of the same near-degenerate set.  *Two separated molecules*: the softest modes
there are the ones between the fragments, the weights decay with distance so
the model barely sees them, and it takes five modes to reach 0.975 where three
suffice elsewhere.

And -- the limit worth knowing before relying on any of the above -- **large
flexible molecules**, where it stops finding the soft space at all.  A
67-atom tetrapeptide puts the model's softest mode only 0.613 inside the true
softest five, and imatinib, 68 atoms with several rotatable bonds joining its
rings, 0.810.  Cholesterol at 74 atoms gives 0.977, so it is not size: it is
how many soft torsions there are to confuse with each other.  A guess is worth
what it is worth on the molecule in hand, and on that class it is not worth
much.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

#: Force constants, hartree per bohr^2 for a stretch and per rad^2 for the
#: other two.  Lindh's, unchanged.
STRETCH_HOLDS = 0.45
BEND_HOLDS = 0.15
TORSION_HOLDS = 0.005

#: An internal weighted below this contributes less than the noise in a
#: differenced Hessian and is left out.  It is also what keeps the count
#: linear in the number of atoms: at this threshold an atom has a handful of
#: neighbours whatever the size of the molecule, so the quadruples are
#: enumerated over neighbours rather than over every four atoms.
WORTH_KEEPING = 1e-4

#: A bend at 0 or 180 degrees and a torsion through one have no derivative;
#: the internal is dropped rather than divided by nothing.
STRAIGHT_WITHIN = 0.15

#: How many internals are assembled into the matrix at a time.  The working
#: array is ``chunk * 144`` doubles, so this is a memory bound and not a
#: correctness one.
AT_A_TIME = 4096

_STEP = 1.0e-5

_ALPHA: Dict[Tuple[int, int], float] = {
    (1, 1): 1.0000, (1, 2): 0.3949, (1, 3): 0.3949,
    (2, 2): 0.2800, (2, 3): 0.2800, (3, 3): 0.2800}
_REFERENCE: Dict[Tuple[int, int], float] = {
    (1, 1): 1.35, (1, 2): 2.10, (1, 3): 2.53,
    (2, 2): 2.87, (2, 3): 3.40, (3, 3): 3.40}


def _period(number: int) -> int:
    """Which of Lindh's three rows an element belongs to.

    The table stops at the third period.  Anything heavier is treated as the
    third, which is an extrapolation and not a lookup -- measured on a
    manganese complex with four bromines, whose softest model mode still came
    out 0.948 against the differenced one, so the extrapolation is reported
    here rather than refused.
    """
    if number <= 2:
        return 1
    return 2 if number <= 10 else 3


def _closeness(numbers: Sequence[int], bohr: np.ndarray) -> np.ndarray:
    """``rho`` for every pair: one where two atoms are bonded, zero far off."""
    rows = np.array([_period(int(one)) for one in numbers])
    gap = np.linalg.norm(bohr[:, None, :] - bohr[None, :, :], axis=2)
    alpha = np.zeros_like(gap)
    reference = np.zeros_like(gap)
    for first in (1, 2, 3):
        for second in (1, 2, 3):
            key = (min(first, second), max(first, second))
            here = (rows[:, None] == first) & (rows[None, :] == second)
            alpha[here] = _ALPHA[key]
            reference[here] = _REFERENCE[key]
    close = np.exp(alpha * (reference ** 2 - gap ** 2))
    np.fill_diagonal(close, 0.0)
    return close


def _stretch(spots: np.ndarray) -> np.ndarray:
    return np.linalg.norm(spots[:, 0] - spots[:, 1], axis=1)


def _bend(spots: np.ndarray) -> np.ndarray:
    """The angle at the middle atom of each triple."""
    one = spots[:, 0] - spots[:, 1]
    two = spots[:, 2] - spots[:, 1]
    one = one / np.linalg.norm(one, axis=1, keepdims=True)
    two = two / np.linalg.norm(two, axis=1, keepdims=True)
    return np.arccos(np.clip(np.einsum('ij,ij->i', one, two), -1.0, 1.0))


def _torsion(spots: np.ndarray) -> np.ndarray:
    """The dihedral of each quadruple, signed, in radians."""
    first = spots[:, 0] - spots[:, 1]
    middle = spots[:, 1] - spots[:, 2]
    last = spots[:, 3] - spots[:, 2]
    across = np.cross(first, middle)
    beyond = np.cross(last, middle)
    return np.arctan2(
        np.linalg.norm(middle, axis=1) * np.einsum('ij,ij->i', first, beyond),
        np.einsum('ij,ij->i', across, beyond))


def _across(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Two directions across a near-straight run, spanning the plane.

    *Which* two does not matter, and that is what makes this legitimate rather
    than a choice smuggled in: both are added with the same force constant, and
    a sum of outer products over an orthonormal pair spans their plane the same
    way whichever pair is picked.  So the result does not depend on how the
    molecule happens to be turned, even though these vectors do.
    """
    axis = axis / np.linalg.norm(axis, axis=1, keepdims=True)
    apart = np.zeros_like(axis)
    apart[np.arange(axis.shape[0]), np.abs(axis).argmin(axis=1)] = 1.0
    first = apart - axis * np.einsum('ij,ij->i', apart, axis)[:, None]
    first = first / np.linalg.norm(first, axis=1, keepdims=True)
    return first, np.cross(axis, first)


def _lean(spots: np.ndarray, way: np.ndarray) -> np.ndarray:
    """How far a near-straight run leans, measured across *way*.

    An angle cannot be used at 180 degrees: its derivative divides by the sine
    and there is none.  The sum of the two unit vectors is zero when the run is
    exactly straight and grows with the lean, so it differentiates where the
    angle will not.
    """
    one = spots[:, 0] - spots[:, 1]
    two = spots[:, 2] - spots[:, 1]
    one = one / np.linalg.norm(one, axis=1, keepdims=True)
    two = two / np.linalg.norm(two, axis=1, keepdims=True)
    return np.einsum('ij,ij->i', one + two, way)


def _slopes(what, spots: np.ndarray, circular: bool = False) -> np.ndarray:
    """Each internal's derivative with respect to the atoms it involves.

    Taken by differencing rather than from a formula, and differenced for
    every internal at once so it costs a couple of dozen array operations for
    the whole molecule.  No quantum chemistry hangs off these, so the
    difference is free; a dihedral's analytic derivative, by contrast, is
    exactly the kind of expression whose sign goes wrong without saying so.
    """
    count = spots.shape[1]
    out = np.zeros_like(spots)
    for atom in range(count):
        for axis in range(3):
            up = spots.copy()
            up[:, atom, axis] += _STEP
            down = spots.copy()
            down[:, atom, axis] -= _STEP
            rise = what(up) - what(down)
            if circular:
                # A dihedral runs round a circle, and either side of the seam
                # the plain difference of two angles is a whole turn out.
                rise = (rise + np.pi) % (2.0 * np.pi) - np.pi
            out[:, atom, axis] = rise / (2.0 * _STEP)
    return out


def _gather(built: np.ndarray, who: np.ndarray, slope: np.ndarray,
            holds: np.ndarray) -> None:
    """Add ``k b b^T`` for a family of internals into the whole matrix."""
    if who.shape[0] == 0:
        return
    wide = built.shape[0]
    spread = 3 * who.shape[1]
    for begin in range(0, who.shape[0], AT_A_TIME):
        end = begin + AT_A_TIME
        here = who[begin:end]
        arm = (slope[begin:end].reshape(here.shape[0], spread)
               * np.sqrt(holds[begin:end])[:, None])
        where = (here[:, :, None] * 3
                 + np.arange(3)[None, None, :]).reshape(here.shape[0], spread)
        block = arm[:, :, None] * arm[:, None, :]
        flat = (where[:, :, None] * wide + where[:, None, :]).reshape(-1)
        built += np.bincount(flat, weights=block.reshape(-1),
                             minlength=built.size).reshape(built.shape)


def _internals(close: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Which pairs, triples and quadruples are worth including.

    Enumerated through each atom's near neighbours rather than over all of
    them, which is what keeps this linear in the number of atoms: every four
    atoms would be 10.5 million combinations for a 57-atom complex, and its
    neighbour lists give 11210.
    """
    count = close.shape[0]
    near: List[np.ndarray] = [np.flatnonzero(close[i] > WORTH_KEEPING)
                              for i in range(count)]
    pairs = [(i, j) for i in range(count) for j in near[i] if i < j]
    triples = [(i, j, k)
               for j in range(count)
               for a, i in enumerate(near[j])
               for k in near[j][a + 1:]]
    # ``k > j`` orders the central bond, so a dihedral and its reverse --
    # i-j-k-l and l-k-j-i, which are the same coordinate -- are counted once
    # between them rather than twice.
    quads = [(i, j, k, l)
             for j in range(count)
             for k in near[j] if k > j
             for i in near[j] if i != k
             for l in near[k] if l != j and l != i]
    def shape(rows, width):
        return (np.asarray(rows, dtype=int) if rows
                else np.zeros((0, width), dtype=int))
    return shape(pairs, 2), shape(triples, 3), shape(quads, 4)


def model_hessian(numbers: Sequence[int], bohr: np.ndarray) -> np.ndarray:
    """Lindh's Hessian for this geometry, in hartree per bohr squared.

    *bohr* is ``(N, 3)`` and Cartesian, in bohr, as everything in
    :mod:`delfin.dashboard.climb` is.  What comes back is ``(3N, 3N)``,
    symmetric, and positive semi-definite by construction -- it is a sum of
    outer products with positive weights, so it can no more invent a negative
    curvature than a spring can pull.  That is a property worth naming,
    because it is also the model's limit: it cannot know that the structure it
    is looking at is a transition state.
    """
    spots = np.asarray(bohr, dtype=float).reshape(-1, 3)
    count = spots.shape[0]
    built = np.zeros((3 * count, 3 * count))
    if count < 2:
        return built
    close = _closeness(numbers, spots)
    pairs, triples, quads = _internals(close)

    if pairs.shape[0]:
        holds = STRETCH_HOLDS * close[pairs[:, 0], pairs[:, 1]]
        at = spots[pairs]
        _gather(built, pairs, _slopes(_stretch, at), holds)

    if triples.shape[0]:
        at = spots[triples]
        angle = _bend(at)
        straight = angle >= np.pi - STRAIGHT_WITHIN
        bent = (angle > STRAIGHT_WITHIN) & ~straight
        if bent.any():
            here, where = triples[bent], at[bent]
            holds = (BEND_HOLDS * close[here[:, 0], here[:, 1]]
                     * close[here[:, 1], here[:, 2]])
            _gather(built, here, _slopes(_bend, where), holds)
        if straight.any():
            # A run through 180 degrees still bends -- it is what a nitrile or
            # an alkyne does -- but it has no angle to differentiate, and
            # dropping it leaves the model with nothing holding that motion at
            # all.  Measured before this was here: carbon dioxide came out
            # flat in seven directions where a linear molecule has five, and
            # propargyl alcohol's softest model curvature was 0.653 against a
            # differenced 4.08.  Two leans across the run replace the one
            # angle, and carry the same force constant between them.
            here, where = triples[straight], at[straight]
            holds = (BEND_HOLDS * close[here[:, 0], here[:, 1]]
                     * close[here[:, 1], here[:, 2]])
            first, second = _across(where[:, 2] - where[:, 0])
            for way in (first, second):
                _gather(built, here,
                        _slopes(lambda spots, way=way: _lean(spots, way),
                                where), holds)

    if quads.shape[0]:
        at = spots[quads]
        # Both bends, each against both bounds.  A torsion is undefined when
        # EITHER of them is straight, and the smaller one tested against the
        # upper bound cannot see it: on an acetonitrile the H-C-C bend is
        # 109.4 degrees and the C-C-N is 180.0, so the minimum sails through
        # and the torsion is built anyway.  What it puts in the matrix is not
        # a small error -- the derivative diverges there.  Measured, the
        # largest eigenvalue of the guessed Hessian goes from 1.4 on an ethane
        # to 3.4e+09 on that acetonitrile and 3.1e+09 on a propyne, and the
        # six flat directions a molecule has to have collapse to one.
        one, two = _bend(at[:, :3]), _bend(at[:, 1:])
        open_enough = ((np.minimum(one, two) > STRAIGHT_WITHIN)
                       & (np.maximum(one, two) < np.pi - STRAIGHT_WITHIN))
        quads, at = quads[open_enough], at[open_enough]
        if quads.shape[0]:
            holds = (TORSION_HOLDS * close[quads[:, 0], quads[:, 1]]
                     * close[quads[:, 1], quads[:, 2]]
                     * close[quads[:, 2], quads[:, 3]])
            _gather(built, quads, _slopes(_torsion, at, circular=True), holds)

    return 0.5 * (built + built.T)
