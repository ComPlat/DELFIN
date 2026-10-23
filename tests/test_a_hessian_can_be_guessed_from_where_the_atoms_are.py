"""A saddle search has to start from a Hessian, and one of them is free.

:meth:`delfin.dashboard.climb.Climb.start` buys its first second-derivative
matrix by differencing the gradient: ``6N`` of them, which is 72 for butane,
210 for a 37-atom anion and 330 for a manganese complex -- thirteen minutes on
this box before the climb takes its first step, and the user watching a
"To the saddle" that has not moved.

:mod:`delfin.dashboard.model_hessian` guesses one from the geometry instead,
after Lindh (*Chem. Phys. Lett.* **241** (1995) 423).  What that guess is worth
was measured against differenced Hessians under GFN2 over twelve structures
from 14 to 74 atoms -- drags, a hydrogen-bonded dipeptide, a floppy chain, two
molecules held apart, an anion, a macrocycle, a sugar, a peptide, a drug, a
steroid -- and a manganese complex with four bromines.  The finding these
tests pin down is the *shape* of what it gets right: the soft subspace, not
the ordering inside it and not the magnitudes -- and not, on a large flexible
molecule, the subspace either.

What is checked here is what can be checked without a quantum chemistry
package: that the matrix is the kind of object it claims to be, that it does
not depend on where the molecule happens to sit or which way it faces, that
the fast enumeration agrees with the obvious slow one, and that the climb
still starts from the measured matrix unless it is asked not to.
"""

import numpy as np
import pytest

from delfin.dashboard import climb
from delfin.dashboard.model_hessian import (
    BEND_HOLDS, STRETCH_HOLDS, TORSION_HOLDS, WORTH_KEEPING,
    _bend, _closeness, _period, _slopes, _stretch, _torsion, model_hessian)

_needs_xtb = pytest.mark.skipif(
    not climb.have_fast_gradients() and climb._gfn.find_xtb() is None,
    reason='no xtb to take gradients from')

BUTANE = """14
butane
C   -1.9260   -0.2570    0.0000
C   -0.6420    0.5680    0.0000
C    0.6420   -0.5680    0.0000
C    1.9260    0.2570    0.0000
H   -2.8180    0.3760    0.0000
H   -1.9600   -0.8930    0.8880
H   -1.9600   -0.8930   -0.8880
H   -0.6320    1.2170    0.8830
H   -0.6320    1.2170   -0.8830
H    0.6320   -1.2170    0.8830
H    0.6320   -1.2170   -0.8830
H    2.8180   -0.3760    0.0000
H    1.9600    0.8930    0.8880
H    1.9600    0.8930   -0.8880
"""

WATER = """3
water
O    0.0000    0.0000    0.1173
H    0.0000    0.7572   -0.4692
H    0.0000   -0.7572   -0.4692
"""


def _at(text):
    found = climb._elements(text)
    return found['numbers'], found['angstrom'] / climb.BOHR


def _turned(spots, angle=0.7, shift=(3.0, -1.0, 2.0)):
    """The same molecule somewhere else, facing another way."""
    turn = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                     [np.sin(angle), np.cos(angle), 0.0],
                     [0.0, 0.0, 1.0]])
    return spots @ turn.T + np.asarray(shift)


def test_the_guess_is_a_hessian_and_not_merely_a_matrix():
    """Symmetric, never negative, and blind in exactly six directions.

    It is a sum of ``k b b^T`` with every ``k`` positive, so it can no more
    produce a negative curvature than a spring can pull -- which is worth
    stating as a test because it is also the model's limit: it cannot know
    that the structure it is looking at is a transition state, and a caller
    who reads an imaginary frequency out of it is reading its own assumption
    back.

    The six flat directions are the three translations and three rotations,
    which cost no energy and must cost none here either.
    """
    numbers, spots = _at(BUTANE)
    built = model_hessian(numbers, spots)

    assert built.shape == (3 * len(numbers), 3 * len(numbers))
    assert np.abs(built - built.T).max() == 0.0

    value = np.linalg.eigvalsh(built)
    assert value.min() > -1e-9, value[:3]
    assert int((np.abs(value) < 1e-9).sum()) == 6, value[:8]


def test_a_linear_molecule_is_flat_in_five_directions_and_not_six():
    """Falls out of the arithmetic rather than being special-cased.

    A rotation about a linear molecule's own axis moves nothing, so there is
    one fewer flat direction, and nothing here counts to six.
    """
    numbers, spots = _at("3\ncarbon dioxide\nO 0 0 -1.16\nC 0 0 0\n"
                         "O 0 0 1.16\n")
    value = np.linalg.eigvalsh(model_hessian(numbers, spots))
    assert int((np.abs(value) < 1e-9).sum()) == 5, value[:8]


def test_where_the_molecule_sits_and_which_way_it_faces_change_nothing():
    """Move it and turn it: the same Hessian arrives, turned with it.

    Not a formality.  The internals are built from distances and angles, which
    are invariant, but the derivatives are taken in Cartesian coordinates,
    which are not -- so an error in the assembly shows up here and nowhere
    else.
    """
    numbers, spots = _at(BUTANE)
    here = model_hessian(numbers, spots)
    angle = 0.7
    turn = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                     [np.sin(angle), np.cos(angle), 0.0],
                     [0.0, 0.0, 1.0]])
    there = model_hessian(numbers, _turned(spots, angle))

    whole = np.kron(np.eye(len(numbers)), turn)
    assert np.abs(whole @ here @ whole.T - there).max() < 1e-8
    assert np.abs(np.linalg.eigvalsh(here)
                  - np.linalg.eigvalsh(there)).max() < 1e-9


def test_the_fast_enumeration_finds_what_the_obvious_slow_one_finds():
    """Neighbour lists against every-combination, on a molecule small enough.

    The enumeration walks each atom's near neighbours rather than every pair,
    triple and quadruple, which is what keeps it linear in the number of atoms
    -- 10.5 million combinations for a 57-atom complex against 23426 that
    matter.  The saving is only worth having if the two agree, so here they
    are asked to.
    """
    numbers, spots = _at(WATER)
    count = len(numbers)
    close = _closeness(numbers, spots)

    slow = np.zeros((3 * count, 3 * count))
    for i in range(count):
        for j in range(count):
            if i < j and close[i, j] > WORTH_KEEPING:
                at = spots[[i, j]][None, :, :]
                arm = np.zeros(3 * count)
                arm[[3 * i, 3 * i + 1, 3 * i + 2,
                     3 * j, 3 * j + 1, 3 * j + 2]] = \
                    _slopes(_stretch, at)[0].reshape(-1)
                slow += STRETCH_HOLDS * close[i, j] * np.outer(arm, arm)
    for j in range(count):
        for i in range(count):
            for k in range(i + 1, count):
                if j in (i, k) or close[i, j] <= WORTH_KEEPING \
                        or close[j, k] <= WORTH_KEEPING:
                    continue
                at = spots[[i, j, k]][None, :, :]
                arm = np.zeros(3 * count)
                for slot, atom in enumerate((i, j, k)):
                    arm[3 * atom:3 * atom + 3] = _slopes(_bend, at)[0, slot]
                slow += (BEND_HOLDS * close[i, j] * close[j, k]
                         * np.outer(arm, arm))

    assert np.abs(model_hessian(numbers, spots) - slow).max() < 1e-8


def test_a_dihedral_is_one_coordinate_and_is_counted_once():
    """i-j-k-l and l-k-j-i are the same angle, and were once both counted.

    A prototype enumerated the central bond in both directions, which doubled
    every torsion and so quietly doubled :data:`TORSION_HOLDS` against the
    other two families.  Counting each once is not only correct but measurably
    better: butane's softest model mode went from 0.950 to 0.975 against the
    differenced one, and a heptane whose model mode had been 0.001 in the true
    softest -- it had picked a different member of a near-degenerate pair --
    came out at 0.999.
    """
    numbers, spots = _at(BUTANE)
    close = _closeness(numbers, spots)
    near = [np.flatnonzero(close[i] > WORTH_KEEPING) for i in range(len(numbers))]
    seen = set()
    for j in range(len(numbers)):
        for k in near[j]:
            if k <= j:
                continue
            for i in near[j]:
                if i == k:
                    continue
                for last in near[k]:
                    if last in (i, j):
                        continue
                    here = (i, j, int(k), int(last))
                    assert here[::-1] not in seen, here
                    seen.add(here)
    assert seen


def test_the_torsions_are_the_gentlest_of_the_three_families():
    """The ordering of Lindh's force constants, which is the whole physics.

    A bond resists far more than an angle and an angle far more than a
    twist -- ninety to one between the ends of it -- and that ordering is why
    the model's softest modes come out as torsions, which is what a search
    wants to be told about.
    """
    assert STRETCH_HOLDS > BEND_HOLDS > TORSION_HOLDS
    assert STRETCH_HOLDS / TORSION_HOLDS == 90.0


def test_elements_past_the_table_are_treated_as_the_last_row_it_has():
    """Lindh's parameters stop at the third period; heavier is extrapolated.

    Reported rather than refused, because it was measured: a manganese complex
    with four bromines -- both elements past the table -- still gave a softest
    model mode sitting 0.948 in the differenced one and 0.993 in its softest
    three.  Refusing it would have thrown away the case where the saving is
    largest, a Hessian that costs 330 gradients and thirteen minutes.
    """
    assert _period(1) == 1 and _period(2) == 1
    assert _period(6) == 2 and _period(8) == 2 and _period(10) == 2
    assert _period(11) == 3 and _period(17) == 3
    assert _period(25) == 3 and _period(35) == 3 and _period(92) == 3


def test_far_apart_atoms_are_left_out_rather_than_added_as_nothing():
    """The weights decay, and below the threshold the internal is dropped.

    This is what keeps the count linear, and it is also the model's one weak
    case: two molecules held apart have their softest modes *between* them,
    and those are exactly the internals this drops.  Measured on anthracene
    and maleic anhydride 2.3 A apart, the model's softest mode sat 0.779 in
    the true softest three where every connected structure gave 0.93 or
    better -- and 0.979 in the softest five, which is how a caller should use
    it there.
    """
    numbers, spots = _at(BUTANE)
    apart = np.array(spots)
    apart[7:] += np.array([0.0, 0.0, 40.0])
    close = _closeness(numbers, apart)
    crossing = close[:7, 7:]
    assert crossing.max() < WORTH_KEEPING, crossing.max()

    built = model_hessian(numbers, apart)
    assert np.abs(built[:21, 21:]).max() == 0.0


def test_a_climb_still_pays_for_its_hessian_unless_it_is_asked_not_to():
    """The default is the behaviour every climb has had, and a typo is refused.

    The switch exists so the swap can be measured against the shipped
    behaviour on structures that matter rather than argued about; until that
    is done the default stays where it was.  An unknown value is refused at
    construction because falling back to the expensive path is precisely the
    mistake a typo here would hide.
    """
    assert climb.FIRST_HESSIAN == 'measured'
    assert climb.HESSIANS_START_FROM == ('measured', 'model', 'corrected')
    with pytest.raises(ValueError, match='measured or model'):
        climb.Climb(BUTANE, method='gfn2', first_hessian='guess')


@_needs_xtb
def test_the_guessed_start_costs_no_gradients_and_the_measured_one_costs_6n():
    """What the swap is for, counted rather than asserted.

    Butane is 14 atoms, so a differenced Hessian is 84 gradients and the
    guess is none.  Both climbs then step identically -- the refreshes along
    the way stay measured, because by then the climb is somewhere the model
    has no claim to know about.
    """
    guessed = climb.Climb(BUTANE, method='gfn2', cores=2,
                          first_hessian='model')
    try:
        opened = guessed.start()
        assert opened['from'] == 'model'
        assert opened['gradients'] <= 2, opened['gradients']
        assert guessed.hessian is not None
        assert np.abs(guessed.hessian - guessed.hessian.T).max() < 1e-12
    finally:
        guessed.close()

    measured = climb.Climb(BUTANE, method='gfn2', cores=2)
    try:
        opened = measured.start()
        assert opened['from'] == 'measured'
        assert opened['gradients'] >= 6 * 14, opened['gradients']
    finally:
        measured.close()


@_needs_xtb
def test_the_correction_buys_back_the_one_thing_a_sum_of_springs_cannot_have():
    """A model Hessian has no negative curvature.  A saddle search needs one.

    :func:`model_hessian` is a sum of ``k b b^T`` with positive ``k``, so its
    lowest eigenvalue is zero and never less -- which is fine for starting an
    *optimisation* and is precisely wrong for starting a *climb*, because the
    climb's first act is to pick a mode to go up and there is nothing pointing
    that way.

    ``'corrected'`` measures the model's own softest modes -- two gradients
    each, twenty in all, against the 96 a differenced Hessian costs for
    sixteen atoms -- and puts the answer back in their place.  Measured on a
    Diels-Alder held at 1.10 A: the guess says +0.0037 where the surface says
    -0.0080, and the correction says -0.0071.
    """
    text = _DIELS_ALDER_HELD

    guessed = climb.Climb(text, method='gfn2', cores=2, first_hessian='model')
    corrected = climb.Climb(text, method='gfn2', cores=2,
                            first_hessian='corrected')
    try:
        was = corrected.engine.calls
        matrix = corrected.first_matrix()
        spent = corrected.engine.calls - was
        assert spent == 2 * climb.MEASURE_THE_SOFTEST, spent
        assert spent < 6 * len(corrected.symbols)

        edge = corrected._basis()
        soft = float(np.linalg.eigvalsh(edge.T @ matrix @ edge)[0])
        guess = float(np.linalg.eigvalsh(
            edge.T @ guessed.first_matrix() @ edge)[0])
        assert guess > 0.0, guess
        assert soft < 0.0, soft
    finally:
        guessed.close()
        corrected.close()


_DIELS_ALDER_HELD = """16
butadiene and ethene, pushed together by 1.10 A
C      -1.511770730818    -0.078246613170     0.455478938079
C      -0.727043823107     0.952677742712     0.160901830656
C       0.726619700240     0.953646170753     0.162987634387
C       1.511848933658    -0.076870712366     0.457589167929
H      -2.584826367129     0.006765122345     0.431252533856
H      -1.117792846415    -1.040460185517     0.733672017164
H      -1.183044319988     1.897378083656    -0.109168104333
H       1.182127851949     1.899269641514    -0.104743334819
H       2.585025855997     0.009342483188     0.435050712651
H       1.118707557773    -1.040095122972     0.731436720523
C      -0.658327621243    -0.317760213959     2.589418492198
C       0.658004039895    -0.316146847946     2.586413586994
H      -1.229591522504     0.583529699796     2.444725727827
H      -1.228837918327    -1.220460679314     2.736580781023
H       1.226373947517     0.586709172705     2.437742888343
H       1.231462082883    -1.217451713317     2.730574479147
"""


@_needs_xtb
def test_correcting_every_direction_gives_back_the_measured_hessian():
    """The correction's own control, and the one that makes the rest legible.

    ``'corrected'`` replaces the guess inside the space it measures and leaves
    the guess everywhere else.  Asked to measure *every* direction the climb
    works in, there is no "everywhere else" left, so what comes back must be
    the differenced Hessian -- the same matrix the shipped start pays 6N
    gradients for, arrived at through entirely different arithmetic.

    Without this, a sweep of how many directions are worth measuring could not
    be read: agreement failing at every count would be indistinguishable from
    the correction being wrong.

    The two do not agree exactly and cannot.  Both are differences of
    gradients, but along different displacements -- one along the Cartesian
    axes, the other along the modes -- so each carries its own truncation and
    its own residual from the self-consistent field.  Measured on water: 0.002
    where the matrix entries are around 3, which is the 0.2 percent the
    difference product was separately measured to be good to.  Demanding
    better than that is demanding better than a finite difference can be.
    """
    walk = climb.Climb(WATER, method='gfn2', cores=2, first_hessian='corrected')
    was = climb.MEASURE_THE_SOFTEST
    try:
        edge = walk._basis()
        climb.MEASURE_THE_SOFTEST = edge.shape[1]
        whole = walk.first_matrix()
        measured = walk.numerical_hessian()
        here = edge.T @ whole @ edge
        there = edge.T @ measured @ edge
        assert np.abs(here - there).max() < 0.01 * np.abs(there).max(), (
            np.abs(here - there).max(), np.abs(there).max())
    finally:
        climb.MEASURE_THE_SOFTEST = was
        walk.close()


def test_a_torsion_through_a_straight_bend_is_left_out():
    """Both bends decide, and the smaller one alone could not see it.

    A torsion i-j-k-l has no derivative when EITHER of its bends is straight,
    and the gate tested the smaller of the two against both bounds.  On an
    acetonitrile the H-C-C bend is 109.4 degrees and the C-C-N is 180.0: the
    minimum passes, the torsion is built, and its derivative diverges there.

    What that costs is not a small error.  The largest eigenvalue of the
    guessed Hessian went from 1.4 on an ethane to 3.4e+09 on that acetonitrile
    and 3.1e+09 on a propyne, and the six flat directions every molecule has --
    three translations and three rotations -- collapsed to one.  As the
    starting matrix for a climb that is not an inaccurate guess, it is noise.
    """
    bohr = 0.52917721092
    cases = {
        'acetonitrile': (
            np.array([6, 6, 7, 1, 1, 1]),
            np.array([[0.0, 0.0, 0.00], [0.0, 0.0, 1.46], [0.0, 0.0, 2.62],
                      [1.02, 0.0, -0.36], [-0.51, 0.884, -0.36],
                      [-0.51, -0.884, -0.36]])),
        'propyne': (
            np.array([6, 6, 6, 1, 1, 1, 1]),
            np.array([[0.0, 0.0, 0.00], [0.0, 0.0, 1.46], [0.0, 0.0, 2.66],
                      [0.0, 0.0, 3.72], [1.02, 0.0, -0.36],
                      [-0.51, 0.884, -0.36], [-0.51, -0.884, -0.36]])),
    }
    for name, (numbers, positions) in cases.items():
        guessed = model_hessian(numbers, positions / bohr)
        spectrum = np.linalg.eigvalsh(guessed)
        assert abs(spectrum).max() < 100.0, (
            f'{name}: a straight bend put {abs(spectrum).max():.3g} into the '
            f'matrix, which is a divided-by-nothing derivative and not a '
            f'force constant')
        flat = int((np.abs(spectrum) < 1e-9).sum())
        assert flat == 6, (
            f'{name}: {flat} flat directions where a molecule has six')
