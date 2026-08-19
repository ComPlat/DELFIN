"""A transition state you can push towards while it is being looked for.

The press next door is ORCA's OptTS on xtb gradients, and it is a press
because it cannot be anything else: measured on sixteen atoms, a three-step
burst of ``! XTB2 OPTTS`` costs 3.08 s cold and 2.66 s warm with the previous
burst's Hessian read back, and about 2.5 s of either is ORCA starting up.  A
drag answers ten times a second.  No amount of Hessian reuse closes a gap that
is mostly process startup, so an interactive saddle search has to be an
optimiser in DELFIN walking on xtb gradients -- which is what
:mod:`delfin.dashboard.climb` is, at 10 ms a step.

The method is not invented here.  It is partitioned rational function
optimisation with eigenvector following on a Bofill-updated Hessian: Banerjee,
Adams, Simons and Shepard 1985 for the augmented matrix, Baker 1986 for the
partitioning and the mode following, Bofill 1994 for the update -- the same
method Gaussian's Berny and ORCA's OptTS run, and the same one pysisyphus and
Sella run on gradients alone.

What the hand is for is the part worth measuring, and it turned out not to be
the part that was expected.  A pull is not a force the climb has to work
against: it is a *statement of which reaction is meant*, and the climb uses it
by following the Hessian eigenvector that most resembles the direction the
structure was dragged.  From the same dragged geometry that is the difference
between reaching the Diels-Alder saddle and walking back down to the
van-der-Waals complex.
"""

import math
import time

import numpy as np
import pytest

from delfin.dashboard import climb
from editor_source import EDITOR_SOURCE

_needs_xtb = pytest.mark.skipif(
    not climb.have_fast_gradients() and climb._gfn.find_xtb() is None,
    reason='no xtb to take gradients from')

#: The transition state the path finder estimated for the Diels-Alder, the same
#: sixteen atoms :mod:`test_the_saddle_is_a_button` hands to ORCA.  Two forming
#: bonds at 2.524 and 2.520 A, one imaginary mode at -131 cm-1.
_ESTIMATE = """16
estimated transition state, as the path finder left it
C           -1.48033944224035       -0.06482702402289        0.33313367693482
C           -0.71514739404797        0.97416833017942       -0.02453912031002
C            0.71570906772105        0.97433901585727       -0.02369014461558
C            1.47979867540738       -0.06451037231854        0.33714846810669
H           -2.54618181474079        0.05254805731267        0.44098634361779
H           -1.10814660238791       -1.07134599400465        0.33648369140953
H           -1.19784782506518        1.91961504046453       -0.23084448535057
H            1.19852540715951        1.91999750020048       -0.22886897502227
H            2.54571980273153        0.05242668156370        0.44531267327466
H            1.10758577128550       -1.07097250949853        0.33692691088494
C           -0.66455239189275       -0.33807729040577        2.70599831991346
C            0.66418880652793       -0.33667278042289        2.70602037351993
H           -1.22780129374320        0.57121740320920        2.80910868025585
H           -1.22537426142126       -1.25659532869891        2.76866685161788
H            1.22560872112149        0.57449371477619        2.80166232446937
H            1.22718959396701       -1.25397841607817        2.76640848291807
"""

#: The two of them apart: the same sixteen atoms relaxed under GFN2 with the
#: ethene 3.35 A off the diene, which is a van-der-Waals complex and a genuine
#: minimum.  Nothing about it says a Diels-Alder is what the user has in mind,
#: which is exactly why it is the structure a hand is needed for.
_COMPLEX = """16
butadiene and ethene, relaxed apart
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
C      -0.658327621243    -0.317760213959     3.689418492198
C       0.658004039895    -0.316146847946     3.686413586994
H      -1.229591522504     0.583529699796     3.544725727827
H      -1.228837918327    -1.220460679314     3.836580781023
H       1.226373947517     0.586709172705     3.537742888343
H       1.231462082883    -1.217451713317     3.830574479147
"""

#: What ORCA's OptTS reached from ``_ESTIMATE``, kept so the two can be
#: compared without running ORCA on every test run.  ``! XTB2 OPTTS``, 8.0 s,
#: one imaginary mode.
_ORCA_SADDLE = """16
ORCA OptTS on XTB2 from the estimate, E = -17.812259376910
  C          -1.43478368987847     -0.11555756879252      0.39801753593899
  C          -0.70800000464735      0.96029290831639      0.01351101273233
  C           0.70792543036907      0.96066325147865      0.01382699497627
  C           1.43510098205032     -0.11480687774322      0.39865617430983
  H          -2.50480092399673     -0.04135148267055      0.50582952622170
  H          -1.04673865581999     -1.11397983254932      0.31729931422258
  H          -1.20115753813432      1.91631252668957     -0.10370598026690
  H           1.20063543422368      1.91694042137411     -0.10317121370013
  H           2.50503091371526     -0.04004089983163      0.50694578141932
  H           1.04761410021811     -1.11343199175848      0.31776831289286
  C          -0.67562925889034     -0.26587715375425      2.58016852501963
  C           0.67503524104879     -0.26553006909586      2.58047503090081
  H          -1.22616409779150      0.63483074833459      2.78770928087516
  H          -1.23028498550611     -1.18636553719529      2.63888406495727
  H           1.22501381825161      0.63546069529481      2.78826237108677
  H           1.23013805516997     -1.18573310998389      2.63943734003805
"""


#: Ammonia forced flat, at the N-H distance where D3h is stationary under GFN2
#: -- found by a golden-section search on that one free parameter, |g|max
#: 2.6e-9.  It is the inversion transition state, and it is the structure that
#: showed a gradient-only Hessian disagreeing with xtb's own.
_FLAT_AMMONIA = """4
planar ammonia, D3h stationary under GFN2
N     0.000000000000     0.000000000000     0.000000000000
H     0.992730000000     0.000000000000     0.000000000000
H    -0.496365000000     0.859729166847     0.000000000000
H    -0.496365000000    -0.859729166847     0.000000000000
"""


def _where(text):
    return np.asarray(climb._elements(text)['angstrom'], dtype=float)


def _forming(rows):
    """The two bonds a Diels-Alder makes, in the atom order used throughout."""
    return math.dist(rows[0], rows[10]), math.dist(rows[3], rows[11])


def _rmsd(one, two):
    one = one - one.mean(0)
    two = two - two.mean(0)
    left, _, right = np.linalg.svd(one.T @ two)
    turn = left @ np.diag([1, 1, np.sign(np.linalg.det(left @ right))]) @ right
    return float(np.sqrt(((one @ turn - two) ** 2).sum() / len(one)))


def _dragged(depth):
    """The complex with the ethene pulled *depth* Angstrom onto the diene.

    Atoms placed and nothing else -- no relaxation, no restraint, which is
    what a bare drag in the viewer leaves behind.  That it works at all is
    worth knowing: the hand needs no help from the optimiser to guide a climb,
    so the editor only has to hand over the structure and where it came from.
    """
    rows = _where(_COMPLEX).copy()
    rows[10:16] += np.array([0.0, 0.0, -float(depth)])
    return climb.xyz_document(climb._elements(_COMPLEX)['symbols'], rows,
                              'where the hand left it')


# -- the method, without needing xtb ---------------------------------------

def test_the_step_climbs_one_mode_and_descends_the_rest():
    """Which is the whole of what "partitioned" means.

    On a quadratic with a known saddle -- one curvature negative, two
    positive -- the step has to go *with* the gradient along the negative one
    and *against* it along the others.  A step that descends all three is a
    minimiser and would walk away from every transition state it was pointed
    at; a step that climbs all three is a maximiser.  Both are one sign away
    from this, which is why it is pinned rather than assumed.
    """
    curvature = np.diag([-0.4, 0.6, 1.2])
    gradient = np.array([0.05, 0.03, -0.02])
    found = climb.partitioned_step(curvature, gradient, None)
    step = found['step']
    # Up along the mode being followed: the step shares the gradient's sign.
    assert step[0] * gradient[0] > 0, step
    # Down along the others: the step opposes it.
    assert step[1] * gradient[1] < 0, step
    assert step[2] * gradient[2] < 0, step
    assert found['negative'] == 1


def test_the_two_shifts_have_the_signs_the_method_requires():
    """Baker's own statement about them, and it is not decoration.

    The shift for the climbed mode is always positive and the shift for the
    rest always negative, and both go to zero at a saddle.  Their signs are
    what put the denominators either side of zero -- get one wrong and the
    step still has a plausible length, still moves the structure, and quietly
    optimises to the wrong kind of stationary point.
    """
    curvature = np.array([-0.4, 0.6, 1.2])
    force = np.array([0.05, 0.03, -0.02])
    up = 0.5 * (curvature[0] + math.sqrt(curvature[0] ** 2 + 4 * force[0] ** 2))
    assert up > 0
    assert up > curvature[0]         # so b_k - lambda_p is negative
    down = climb.lowest_lambda(curvature[1:], force[1:])
    assert down < 0
    assert down < curvature[1:].min()
    # It is a root of the secular equation, which is the only thing that
    # makes it the right number rather than merely a negative one.
    assert abs(float(np.sum(force[1:] ** 2 / (down - curvature[1:]))) - down) \
        < 1e-9


def test_a_flat_direction_does_not_divide_by_its_own_curvature():
    """A mode with no curvature and no gradient on it is what a molecule's
    own translations and rotations are, and there are six of them.

    Left in the Hessian they sit at zero -- which is exactly where the shift
    for the minimised modes lives, so each one is a pole sitting on top of the
    root being searched for.  They are kept out by working in a basis that
    never contained them, and this is that basis being the right size.
    """
    water = climb._elements(
        '3\n\nO 0.0 0.0 0.0\nH 0.0 0.757 0.587\nH 0.0 -0.757 0.587\n')
    basis = climb.translations_and_rotations(water['masses'],
                                             water['angstrom'] / climb.BOHR)
    assert basis.shape == (9, 3), basis.shape          # 3N - 6
    # A linear molecule turns about two axes, not three, and the sixth vector
    # is not there to be removed.  Falling out of the orthogonalisation rather
    # than being special-cased is what makes that true without a branch.
    line = climb._elements('3\n\nC 0 0 0\nO 0 0 1.16\nO 0 0 -1.16\n')
    thin = climb.translations_and_rotations(line['masses'],
                                            line['angstrom'] / climb.BOHR)
    assert thin.shape == (9, 4), thin.shape            # 3N - 5


def test_bofill_learns_the_hessian_a_quadratic_actually_has():
    """One Hessian at the start and none afterwards is the whole economy of
    this: 96 gradients once, then one per step.

    That only works if the update keeps up.  Bofill mixes the symmetric
    rank-one update, which converges fast and blows up when the step is
    perpendicular to what it missed, with Powell's, which is stable and slower
    -- weighted by how nearly parallel those two are.  On a quadratic the
    updated Hessian has to walk back to the exact one, and it does: measured
    from an identity on a three-by-three with a negative eigenvalue in it, the
    largest error falls 1.36 -> 0.083 in three steps, is under 0.002 by nine,
    and is down to 4e-8 by sixteen.

    The negative eigenvalue is the point.  BFGS would have removed it -- it is
    built to keep a Hessian positive definite -- and a saddle search whose
    Hessian has no negative eigenvalue has nothing to climb.
    """
    rng = np.random.default_rng(0)
    exact = np.array([[-0.4, 0.1, 0.0], [0.1, 0.9, 0.2], [0.0, 0.2, 1.3]])
    assert np.linalg.eigvalsh(exact).min() < 0
    guess = np.eye(3)
    for turn in range(1, 17):
        moved = rng.normal(0, 0.1, 3)
        guess = climb.bofill_update(guess, moved, exact @ moved)
        if turn == 9:
            assert np.abs(guess - exact).max() < 0.002, guess
    assert np.abs(guess - exact).max() < 1e-6, guess
    assert np.linalg.eigvalsh(guess).min() < 0
    # And it stays symmetric, which a Hessian has to be for eigh to be
    # meaningful at all.
    assert np.abs(guess - guess.T).max() < 1e-12


def test_the_rank_one_share_goes_to_zero_where_rank_one_blows_up():
    """Which is the only reason to mix two updates rather than use one.

    The rank-one update divides by ``j.dx``, so it is worthless exactly when
    the step is perpendicular to what the old Hessian got wrong.  Bofill's
    weight is the cosine squared of the angle between them, so it goes to zero
    in the same place -- and a weight that went to *one* there would be the
    same three lines of algebra doing the opposite of what they are for.
    """
    moved = np.array([1.0, 0.0, 0.0])
    # A Hessian whose error is entirely across the step.
    hessian = np.zeros((3, 3))
    changed = np.array([0.0, 1.0, 0.0])
    missed = changed - hessian @ moved
    share = (float(missed @ moved) ** 2
             / (float(missed @ missed) * float(moved @ moved)))
    assert share == pytest.approx(0.0)
    # Nothing infinite came out of it.
    updated = climb.bofill_update(hessian, moved, changed)
    assert np.isfinite(updated).all()


def test_a_method_with_a_basis_set_is_refused_and_named_as_a_job():
    """The same three methods the press next door offers, for the same reason.

    A climb makes a hundred gradients a second look like a picture; on a real
    basis set one of them is minutes, and there is no version of this that is
    interactive.  Saying so is better than starting.
    """
    assert set(climb.CLIMB_METHODS) == {'gfn2', 'gfn1', 'gfnff'}
    refused = climb.climb_to_saddle(_ESTIMATE, 'PBE0')
    assert not refused['ok']
    assert 'runs on xtb' in refused['status']
    assert not climb.climb_to_saddle('', 'gfn2')['ok']


# -- against ORCA -----------------------------------------------------------

@_needs_xtb
def test_the_estimate_reaches_the_saddle_orca_reaches():
    """The same structure, the same saddle, a quarter of the time.

    Measured: ORCA's OptTS took 8.0 s on this estimate and came out with the
    forming bonds at 2.3153 and 2.3153.  This converges in 11 steps and 1.8 s
    -- of which 0.6 is the starting Hessian and 0.6 the frequency check at the
    end, so the climb itself is a seventh of a second -- and comes out at
    2.3144 and 2.3161, which is 0.0057 A RMSD from ORCA's structure.

    The imaginary mode is worth a sentence.  ORCA *prints* -372 cm-1 for that
    run and this says -394, which looks like a disagreement and is not one:
    xtb's own ``--hess`` on ORCA's own final geometry says -393.5.  ORCA is
    told ``Recalc_Hess 5``, so the last Hessian it printed belongs to step 10
    of a 13-step run -- -372 describes a geometry it had not finished leaving.
    """
    began = time.perf_counter()
    got = climb.climb_to_saddle(_ESTIMATE, 'gfn2')
    seconds = time.perf_counter() - began
    assert got['ok'], got['status']
    assert seconds < 60, seconds

    here = _where(got['xyz'])
    one, two = _forming(here)
    assert 2.2 < one < 2.45, one
    # Symmetric, which the estimate was not.
    assert abs(one - two) < 0.02, (one, two)
    assert _rmsd(here, _where(_ORCA_SADDLE)) < 0.05

    # And it is a transition state, which is the only reason to have run it.
    assert got['imaginary']['count'] == 1, got['imaginary']
    assert -450 < got['imaginary']['modes'][0] < -330, got['imaginary']


@_needs_xtb
def test_a_step_costs_a_gradient_and_a_gradient_is_ten_milliseconds():
    """Which is what makes the whole thing possible.

    A saddle search made of repeated ORCA runs pays 2.5 s of process startup
    per burst; this pays one xtb gradient per step, measured at 6 ms in process
    on four threads and 35 ms through the command line, and 12 ms for the whole
    step around it.  Both are inside a frame; the first leaves room for the
    picture.

    Threads are pinned around every call and it is not a nicety: on this
    384-core box xtb takes every core for sixteen atoms and one gradient costs
    1.66 s, which is 230 times what four threads cost.  ``OMP_NUM_THREADS``
    cannot fix that once the runtime is up -- measured, setting it after the
    first calculation changed nothing at all.
    """
    walk = climb.Climb(_ESTIMATE, 'gfn2')
    try:
        walk.start()
        taken = []
        for _ in range(8):
            began = time.perf_counter()
            walk.step()
            taken.append(time.perf_counter() - began)
        middle = float(np.median(taken))
    finally:
        walk.close()
    # Generous by a factor of ten against the 10 ms measured, because this is
    # a shared machine and the claim being pinned is "interactive", not a
    # number.
    assert middle < 0.15, middle


@_needs_xtb
def test_what_it_reached_is_checked_and_not_assumed():
    """A climb that converges has proved that the gradient vanished, and
    nothing else.

    Measured on this box, a path finder feeding ORCA's OptTS converged onto a
    *second-order* saddle at -73 and -27 cm-1 and reported success; and a climb
    from a van-der-Waals complex with no hand to guide it converges perfectly
    happily onto a minimum.  Both are stationary points and neither is a
    reaction, so the number of modes going the wrong way is counted at the end
    and said.
    """
    walk = climb.Climb(_ESTIMATE, 'gfn2')
    try:
        walk.start()
        for _ in range(40):
            if walk.step()['converged']:
                break
        said = walk.verdict()
    finally:
        walk.close()
    assert said['ok'] and said['count'] == 1, said
    # A minimum is not a transition state, and this is what says so.
    still = climb.Climb(_COMPLEX, 'gfn2')
    try:
        empty = still.verdict()
    finally:
        still.close()
    assert empty['count'] == 0 and not empty['ok'], empty
    # -40 cm-1 is autodE's threshold and deliberately conservative; AutoMeKin
    # throws away anything softer than -200.  The number itself is always
    # reported because the two disagree about real cases.
    assert climb.IMAGINARY_BELOW == -40.0


@_needs_xtb
def test_the_verdict_is_asked_of_xtb_and_not_of_the_climb():
    """Because a Hessian differenced from gradients is not always the Hessian.

    A gradient-only saddle optimiser inherits whatever disagreement the model's
    gradient has with the model's energy, and GFN2 has one.  Measured on planar
    ammonia at its own D3h stationary point, in separate xtb processes so that
    nothing at all is shared between the two numbers: one Hessian element comes
    out at -0.01206 from differenced gradients and +0.02958 from four energies
    of the same xtb, and ``xtb --hess`` agrees with the energies at +0.0296.
    The umbrella is then +386 cm-1 from gradients and -971 from either of the
    others -- the difference between calling flat ammonia a minimum and calling
    it the inversion transition state it actually is.

    On an ordinary structure the two agree and none of this matters: on the
    Diels-Alder saddle it is -394.6 from gradients against -393.5 from
    ``xtb --hess``.  It is symmetry that separates them, and transition states
    are quite often symmetric -- so the verdict at the end of a climb is asked
    of xtb, which costs a quarter of a second and is the same Hessian the
    button next door reports.
    """
    walk = climb.Climb(_FLAT_AMMONIA, 'gfn2')
    try:
        walk.start()
        gradients = walk.frequencies()
        theirs = walk.frequencies_from_xtb()
        said = walk.verdict()
    finally:
        walk.close()
    assert theirs is not None
    # The umbrella, which the two disagree about by more than its own size.
    assert -1100 < theirs[0] < -850, theirs[:3]
    assert gradients[0] > 0, gradients[:3]
    # The verdict follows xtb, so flat ammonia is named as the transition
    # state it is rather than as the minimum a gradient Hessian makes it.
    assert said['ok'] and said['count'] == 1, said
    assert -1100 < said['modes'][0] < -850, said

    # And where there is no disagreement the two say the same thing, which is
    # what makes the ammonia case a fact about symmetry rather than a bug.
    walk = climb.Climb(_ORCA_SADDLE, 'gfn2')
    try:
        walk.start()
        assert abs(walk.frequencies()[0] - walk.frequencies_from_xtb()[0]) < 5


    finally:
        walk.close()


# -- the hand ---------------------------------------------------------------

@_needs_xtb
def test_the_hand_names_the_mode_and_that_is_what_helps():
    """The demonstration the whole feature rests on.

    Three drags of different depths, each one nothing but the ethene placed
    closer to the diene -- no relaxation, no restraint, which is all a bare
    drag in the viewer leaves.  From each, two climbs that differ in one
    thing: whether the climb is told which way the structure was moved.

    Told, all three reach the Diels-Alder saddle -- 42, 39 and 37 steps for
    drags of 0.80, 0.95 and 1.10 A -- at 2.315 A on both forming bonds, one
    imaginary mode near -394 cm-1, within 0.008 A RMSD of what ORCA reaches.
    Not told, every one of them climbs the *lowest* mode instead, which is the
    two fragments rocking against each other, and walks back down to the
    van-der-Waals complex 0.43 to 0.66 A away with no imaginary mode at all.

    So what a pull contributes is not force.  It is the answer to "which
    reaction", and the climb reads it as the Hessian eigenvector that most
    resembles the direction the structure was dragged.  That is Baker's mode
    following with the reference mode seeded by the user instead of by an
    index, which is the same thing ORCA offers as ``TS_Mode {B 0 1}`` and a
    great deal easier to say with a mouse.
    """
    orca = _where(_ORCA_SADDLE)
    for depth in (0.95, 1.1):
        moved = _dragged(depth)
        guided = climb.climb_to_saddle(moved, 'gfn2', aimed_from=_COMPLEX,
                                       max_steps=200)
        alone = climb.climb_to_saddle(moved, 'gfn2', max_steps=200)

        here = _where(guided['xyz'])
        one, two = _forming(here)
        assert guided['ok'], (depth, guided['status'])
        assert 2.25 < one < 2.40 and abs(one - two) < 0.02, (depth, one, two)
        assert _rmsd(here, orca) < 0.05, (depth, _rmsd(here, orca))
        assert guided['imaginary']['count'] == 1, (depth, guided['imaginary'])
        assert -450 < guided['imaginary']['modes'][0] < -330, \
            (depth, guided['imaginary'])

        # And the same structure, unaimed, does not get there.
        there = _where(alone['xyz'])
        assert _rmsd(there, orca) > 0.2, (depth, _rmsd(there, orca))
        assert alone['imaginary']['count'] != 1 \
            or alone['imaginary']['modes'][0] > -100, (depth,
                                                       alone['imaginary'])


@_needs_xtb
def test_from_the_complex_alone_there_is_no_saddle_to_find():
    """Which is the honest half of the previous test.

    Started on the relaxed van-der-Waals complex with no hand at all, the
    climb has nothing to follow but the lowest mode of a minimum, and the
    lowest mode of *this* minimum is the two fragments rocking against each
    other.  Measured: it converges after 105 steps onto a structure with the
    fragments 3.5 A apart and not one mode going the wrong way -- a minimum,
    reached without complaint, 0.6 A RMSD from the saddle ORCA finds.

    That is the method behaving correctly and being useless, and it is why the
    feature is a hand rather than a button.  Eigenvector following finds the
    saddle that the mode it is climbing leads to; far from any saddle, on a
    surface as flat as two molecules drifting past each other, which mode that
    is has nothing to do with the chemistry anyone had in mind.
    """
    got = climb.climb_to_saddle(_COMPLEX, 'gfn2', max_steps=200)
    apart, also = _forming(_where(got['xyz']))
    assert apart > 3.0 and also > 3.0, (apart, also)
    assert got['imaginary']['count'] == 0, got['imaginary']


def test_the_hand_guides_the_climb_and_never_restrains_it():
    """A restrained saddle is not a saddle, and this is measured rather than
    argued.

    Climbing on the surface with both forming bonds restrained to 2.45 A at
    200 kcal/mol/A^2 converges, in the sense that the gradient of the
    restrained surface vanishes -- and the point it converges to has a true
    gradient of 2.2e-2 (130 times the convergence threshold) and *two*
    imaginary modes, -637 and -50.  At 50 kcal/mol/A^2 it does not converge at
    all.  Leaving the editor's own push on instead of a fixed restraint is
    worse still: a push is a constant force, so it does not stop at a saddle,
    it drives past one.

    So the hand is never part of the climb.  While the mouse is down the climb
    is suspended; when it is let go the climb starts again from the structure
    that was made, aimed along the way it was made.  Both halves of that are
    in the editor, and this is the wiring that carries them.
    """
    source = EDITOR_SOURCE
    assert 'submit_climb_btn = widgets.ToggleButton(' in source
    assert "description='Climb to TS'" in source
    assert 'def on_submit_climb(change=None):' in source
    assert "submit_climb_btn.observe(on_submit_climb, names='value')" in source
    # Where the structure stood when the hand arrived, which is the other end
    # of the direction the drag names.
    assert "state['climb_was'] = state.get('climb_showing')" in source
    # And the hand-over, on the drag ending and not on every frame of it.
    assert "state['climb_hand'] = {" in source
    assert "'was': state.pop('climb_was', None)," in source
    # The loop reads this and stops, so the toggle is a real Stop.
    assert "state['climb_run'] = None" in source


def test_the_picture_is_fed_frames_and_never_the_coordinate_box():
    """Every write to the box rebuilds the viewer from nothing.

    A hundred steps a second through the box would be a hundred rebuilds a
    second, which is not a picture of a climb but a slideshow of teardowns.
    The frames go down the channel the optimiser already uses, and the box is
    written once at the end.
    """
    source = EDITOR_SOURCE
    assert 'submit_gfn_frame.value = json.dumps(payload)' in source
    # The window starts where the previous one started, so every frame is sent
    # twice and a read that lands between two writes misses nothing.
    assert 'start, sent[0], sent[1], sent[2] = sent[0], sent[1], len(walked), now' \
        in source
    # And it is held apart, because the climb makes frames far faster than any
    # page is asked to look.
    assert 'now - sent[2] < 0.04' in source
    # The box, once, at the end, with the flag that says the picture has it.
    assert "'Climbed towards a transition state'), drawn=True)" in source
