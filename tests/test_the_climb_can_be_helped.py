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


@_needs_xtb
def test_gfn_ff_is_not_allowed_to_write_into_the_directory_it_is_run_from(tmp_path):
    """A dashboard runs where the user launched it, which is often their project.

    Measured: an in-process GFN-FF drops ``gfnff_topo`` -- 142 kB of it -- into
    the process's working directory, and the library takes no working directory
    to be told about.  The command line has one of its own and left nothing at
    all behind in the same test, so GFN-FF goes out that way even where the
    library is installed.  The suite's own guard against writing into the
    checkout is what found this, which is the second reason to keep that guard.
    """
    import os

    was = os.getcwd()
    os.chdir(tmp_path)
    try:
        walk = climb.Climb(_COMPLEX, 'gfnff')
        try:
            assert isinstance(walk.engine, climb._CommandLine), type(walk.engine)
            walk.start()
            walk.step()
        finally:
            walk.close()
        left = sorted(p.name for p in tmp_path.iterdir())
    finally:
        os.chdir(was)
    assert left == [], left
    # GFN2 has no such file and keeps the faster engine.
    quick = climb.Climb(_COMPLEX, 'gfn2')
    try:
        assert isinstance(quick.engine, climb._InProcess) \
            or not climb.have_fast_gradients()
    finally:
        quick.close()


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


@_needs_xtb
def test_a_climb_already_running_takes_what_the_hand_made():
    """The loop the editor actually runs, without the editor.

    Start climbing from the bare complex, let it take a few steps going
    nowhere, then do what a mouse does -- put the ethene where the user wants
    it -- and hand that structure back with the geometry it came from.  The
    climb starts again from there, aimed along the way the hand moved it, and
    reaches the saddle.

    The Hessian is recomputed on the hand-over rather than carried across, and
    that is measured rather than careful: carried across, the same climb still
    gets there but in 62 steps against 15, because a Bofill update repairs a
    Hessian one step at a time and a hand moves further in one gesture than a
    climb does in twenty.  Six tenths of a second once, on the mouse being let
    go, is cheaper than 47 gradients and a picture that wanders on the way.
    """
    walk = climb.Climb(_COMPLEX, 'gfn2')
    try:
        walk.start()
        for _ in range(5):
            walk.step()
        # Nowhere near a saddle, and about to be told where to go.
        assert min(_forming(walk.angstrom)) > 3.0

        walk.took(_dragged(0.95), aimed_from=_COMPLEX)
        arrived = False
        for _ in range(200):
            if walk.step()['converged']:
                arrived = True
                break
        assert arrived, walk.steps
        one, two = _forming(walk.angstrom)
        assert 2.25 < one < 2.40 and abs(one - two) < 0.02, (one, two)
        assert _rmsd(walk.angstrom, _where(_ORCA_SADDLE)) < 0.05
        said = walk.verdict()
    finally:
        walk.close()
    assert said['ok'] and said['count'] == 1, said

    # A structure with a different number of atoms is a different molecule,
    # and the climb says so rather than reading one geometry as another.
    walk = climb.Climb(_COMPLEX, 'gfn2')
    try:
        walk.start()
        refused = walk.took('3\n\nO 0 0 0\nH 0 0.76 0.59\nH 0 -0.76 0.59\n')
    finally:
        walk.close()
    assert not refused['ok'] and 'changed' in refused['status']


def test_the_hand_guides_the_climb_and_never_restrains_it():
    """A restrained saddle is not a saddle, and this is measured rather than
    argued.

    Measured on the Diels-Alder, climbing on surfaces with both forming bonds
    restrained at 200 kcal/mol/A^2.  Held at 2.20 A the climb ends 0.53 A from
    the saddle with a true gradient of 4.2e-2 -- 138 times the threshold it
    thinks it met -- and *two* imaginary modes, -628 and -42.  Held at 2.60 A
    it converges in five steps onto a point 0.11 A away with no imaginary mode
    at all, which is a minimum.  Held at 2.45 A and only 50 kcal/mol/A^2 it
    converges onto a first-order saddle 0.66 A away that is a different saddle.
    Unrestrained, the same climb takes 11 steps and lands 0.006 A from ORCA.

    Leaving the editor's own push on instead of a fixed restraint is worse
    still: a push is a constant force, so it does not stop at a saddle, it
    drives past one.

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
    # of the direction the drag names.  The climb's own last frame if it has
    # made one, and what is in the box otherwise: a hand can arrive before the
    # first Hessian has finished.
    assert "state['climb_was'] = (state.get('climb_showing')" in source
    assert "or _current_xyz())" in source
    # And the hand-over, on the drag ending and not on every frame of it.
    assert "state['climb_hand'] = {" in source
    assert "'was': state.pop('climb_was', None)," in source
    # The loop reads this and stops, so the toggle is a real Stop.
    assert "state['climb_run'] = None" in source


@_needs_xtb
def test_the_button_climbs_and_the_drag_reaches_it(tmp_path):
    """The real part, driven the way a mouse drives it.

    Everything above tests the method; this tests the wiring, because the
    wiring is what a source-reading test cannot see.  A real editor is built,
    the toggle is pressed, and then the browser's own field is written the way
    the page writes it -- a ``DELFIN drag-follow`` frame while the mouse is
    down and a ``DELFIN drag-end`` when it is let go.  The climb has to stop
    for the hand, take the structure the hand made, remember where that
    structure came from, and get to the saddle.

    Measured: from the bare complex it goes nowhere; interrupted by the drag it
    converges to 2.315 A on both forming bonds with one imaginary mode near
    -394, and that is what is left in the coordinate box.
    """
    import time as clock

    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    room = DashboardContext(calc_dir=tmp_path / 'calc',
                            archive_dir=tmp_path / 'archive',
                            office_dir=tmp_path / 'office')
    room.run_js = lambda _script: None
    state = {}
    box = widgets.Textarea(value=_COMPLEX)
    part = structure_editor.build(
        room, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda call, *a, **k: call(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    state['current_xyz_for_copy'] = {'content': _COMPLEX}
    part.submit_ff_dd.value = 'gfn2'

    part.submit_climb_btn.value = True
    began = clock.time()
    while not state.get('climb_showing') and clock.time() - began < 60:
        clock.sleep(0.02)
    assert state.get('climb_run') is not None, 'the climb never started'

    # The hand: down, then let go somewhere else.  The page names the atoms it
    # is holding in the comment line, which is what tells the climb a hand has
    # arrived at all.
    moved = _dragged(0.95).splitlines()
    part.submit_manip_sync.value = '\n'.join(
        [moved[0], 'DELFIN drag-follow held=10,11'] + moved[2:])
    assert state.get('climb_was'), 'the climb did not remember where it stood'
    part.submit_manip_sync.value = '\n'.join(
        [moved[0], 'DELFIN drag-end'] + moved[2:])

    began = clock.time()
    while state.get('climb_run') is not None and clock.time() - began < 300:
        clock.sleep(0.05)
    part.submit_climb_btn.value = False
    assert state.get('climb_run') is None, 'the climb never finished'

    one, two = _forming(_where(box.value))
    assert 2.25 < one < 2.40 and abs(one - two) < 0.03, (one, two)
    assert _rmsd(_where(box.value), _where(_ORCA_SADDLE)) < 0.05
    said = part.mol_status.value
    assert 'transition state' in said, said
    assert 'cm-1' in said, said


def _an_editor(room):
    """A real editor over the van-der-Waals complex, driven from the kernel."""
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room.mkdir(parents=True, exist_ok=True)
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda _script: None
    state = {}
    box = widgets.Textarea(value=_COMPLEX)
    part = structure_editor.build(
        ctx, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda call, *a, **k: call(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    state['current_xyz_for_copy'] = {'content': _COMPLEX}
    part.submit_ff_dd.value = 'gfn2'
    return part, state, box


def _frames_seen(part):
    """Every write to the frame channel, as (run, how many frames)."""
    import json

    seen = []

    def note(change):
        if change.get('name') != 'value' or not change.get('new'):
            return
        try:
            said = json.loads(change['new'])
        except ValueError:
            return
        seen.append((said.get('run'), len(said.get('frames') or [])))

    part.submit_gfn_frame.observe(note, names='value')
    return seen


_SERIAL = [0]


def _command(part, verb, payload=''):
    """What the page sends when it cannot finish a gesture on its own."""
    _SERIAL[0] += 1
    part.submit_cmd_sync.value = f'{verb}:{_SERIAL[0]}:{payload}'


def _drag_towards_the_saddle(part, box, depth=0.95, steps=6, pause=0.25):
    """The page's own protocol, in the order the page sends it.

    ``gfngrab``, then a ``DELFIN drag-follow`` message per mouse move, then
    ``DELFIN drag-end`` carrying the geometry, then ``gfnfree``.  The last two
    are separate messages and arrive in that order because the page sends the
    geometry from its mouseup handler and says the hand has gone from the
    animation frame after it.
    """
    import time as clock

    from delfin.dashboard import climb as _c

    names = _c._elements(_COMPLEX)['symbols']
    start = _where(box.value)
    _command(part, 'gfngrab', '0')
    for turn in range(1, steps + 1):
        moved = start.copy()
        moved[10:16] += np.array([0.0, 0.0, -depth * turn / steps])
        part.submit_manip_sync.value = _c.xyz_document(
            names, moved, f'DELFIN drag-follow held=10,11 n={turn}')
        clock.sleep(pause)
    part.submit_manip_sync.value = _c.xyz_document(
        names, _where(box.value), 'DELFIN drag-end')
    _command(part, 'gfnfree', '')


def _settled(part, state, climb_on, seconds=120.0):
    """Wait for whatever owns the release to finish, then read the box."""
    import time as clock

    began = clock.time()
    while clock.time() - began < seconds:
        if climb_on and state.get('climb_run') is None:
            break
        if not climb_on and clock.time() - began > 25:
            break
        clock.sleep(0.1)
    part.submit_climb_btn.value = False
    clock.sleep(1.5)


@_needs_xtb
def test_a_drag_is_as_visible_with_the_climb_on_as_with_it_off(tmp_path):
    """The first thing the feature has to do is let you see what you are doing.

    It did not.  With Climb to TS on, the climb went on stepping while the
    hand was down -- on the geometry from *before* the drag, because that is
    what it had been given -- so it walked away from the hand while the hand
    was working.  On a structure a long way from any saddle it can converge in
    a second or two, and it then wrote what it had found into the coordinate
    box and switched its own toggle off in the middle of the gesture.
    Measured on this complex: the box came out of the drag at 3.53 A, further
    apart than it started, where the same drag with the climb off left it at
    2.45.

    So the climb now stands completely still while a hand is down: no step, no
    frame, no finishing.  Two things are pinned here, and neither is a frame
    count -- how many follow rounds fit inside a drag depends on what else the
    machine is doing, and a number that moves with the load is not a fact about
    this code.

    What is exact is *whose* frames those are.  During the drag every frame on
    the channel must belong to the run the hand being followed claimed, and
    none to the climb's; and the structure the drag leaves must be the one the
    hand put there.  Measured either way, with the climb off and on: 2.450 and
    2.445 A.
    """
    quiet, _, box = _an_editor(tmp_path / 'off')
    quiet.submit_relax_btn.value = True
    seen = _frames_seen(quiet)
    _drag_towards_the_saddle(quiet, box)
    without = [run for run, count in seen if count]
    plain = _forming(_where(box.value))
    assert without, 'Dynamik Opt drew nothing at all, so nothing is being read'

    busy, state, box = _an_editor(tmp_path / 'on')
    busy.submit_relax_btn.value = True
    busy.submit_climb_btn.value = True
    began = time.time()
    while not state.get('climb_showing') and time.time() - began < 120:
        time.sleep(0.02)
    assert state.get('climb_run') is not None, 'the climb never started'
    mine = state.get('climb_frame_run')
    seen = _frames_seen(busy)
    _drag_towards_the_saddle(busy, box)
    withal = [run for run, count in seen if count]
    climbed = _forming(_where(box.value))

    # The hand is still drawn, and the climb has not drawn over it.
    assert withal, 'the drag put nothing on the channel with the climb on'
    assert all(run != mine for run in withal), (mine, withal)
    # And the structure the drag leaves is the hand's, not the climb's.
    assert abs(climbed[0] - plain[0]) < 0.05, (plain, climbed)
    assert abs(climbed[1] - plain[1]) < 0.05, (plain, climbed)
    assert 2.3 < climbed[0] < 2.6, climbed
    _settled(busy, state, True)


@_needs_xtb
def test_letting_go_hands_over_to_the_climb_and_not_to_a_relaxation(tmp_path):
    """"Ich ziehe in TS, lasse los, und es faellt zurueck."

    Two things were armed by one release.  Dynamik Opt's own settle and Auto's
    minimisation walk downhill, the climb walks up along one mode, and downhill
    wins because it is a minimisation of the whole structure against a search
    that moves one direction at a time.  Nothing had ever told them about each
    other.

    Measured here, the same gesture with Dynamik Opt on either way: with the
    climb off the structure comes back to the van-der-Waals complex at 3.353
    and 3.348 A with no mode going the wrong way, and with the climb on it
    reaches 2.314 and 2.315 A with one mode at -394 cm-1.  Which is the
    difference between the feature working and the feature being the thing the
    user reported.
    """
    plain, state, box = _an_editor(tmp_path / 'off')
    plain.submit_relax_btn.value = True
    _drag_towards_the_saddle(plain, box)
    _settled(plain, state, False)
    fell = _forming(_where(box.value))

    part, state, box = _an_editor(tmp_path / 'on')
    part.submit_relax_btn.value = True
    part.submit_climb_btn.value = True
    began = time.time()
    while not state.get('climb_showing') and time.time() - began < 120:
        time.sleep(0.02)
    _drag_towards_the_saddle(part, box)
    _settled(part, state, True)
    reached = _forming(_where(box.value))

    # Without the climb, back to the complex it started from.
    assert fell[0] > 3.0 and fell[1] > 3.0, fell
    # With it, the saddle -- and it says so.
    assert 2.25 < reached[0] < 2.40, reached
    assert abs(reached[0] - reached[1]) < 0.03, reached
    assert _rmsd(_where(box.value), _where(_ORCA_SADDLE)) < 0.05
    said = part.mol_status.value
    assert 'transition state' in said, said
    assert 'cm-1' in said, said


@_needs_xtb
def test_the_climb_claims_a_run_again_every_time_it_starts_drawing(tmp_path):
    """"Ich sehe nicht das Update von Climb to TS im Viewer."

    The run number is how a writer says it is the current one, and a drag
    moves it: the hand being followed takes a number, and the settle after the
    release takes another.  The climb claimed one when its toggle went down
    and held it -- so by the time it had anything to draw it was two behind,
    it failed its own guard on every write, and it threw away every frame it
    made.  Silently, and correctly: it really was not the current run.

    Measured on this gesture: the climb claims run 1, the drag leaves run 2
    behind it, and the climb writes under run 3 afterwards -- so a climb
    holding the toggle-time number would have discarded all 94 frames of its
    own path.  It claims a run at each moment it starts producing frames
    instead, which is when it begins and after every hand-over.
    """
    part, state, box = _an_editor(tmp_path / 'run')
    part.submit_relax_btn.value = True
    part.submit_climb_btn.value = True
    began = time.time()
    while not state.get('climb_showing') and time.time() - began < 120:
        time.sleep(0.02)
    at_the_toggle = state.get('climb_frame_run')
    assert at_the_toggle, state.get('climb_frame_run')

    _drag_towards_the_saddle(part, box)
    after_the_drag = int(state.get('gfn_run', 0))
    seen = _frames_seen(part)
    _settled(part, state, True)

    drawn = [(run, count) for run, count in seen if count]
    assert drawn, 'the climb drew nothing after the release'
    # The drag moved the counter past the number the climb started with...
    assert after_the_drag > at_the_toggle, (at_the_toggle, after_the_drag)
    # ...and the climb took a fresh one rather than writing under a stale one.
    assert all(run > at_the_toggle for run, _ in drawn), (at_the_toggle, drawn)
    # Every frame of it would have been thrown away by the old guard.
    assert sum(count for run, count in drawn if run == at_the_toggle) == 0


def test_only_one_thing_writes_a_geometry_after_a_release():
    """A settle, an auto-minimisation and a climb all answer the same release.

    The invariant is put on each of the four things that would break it rather
    than on the places that call them, because the call sites are not all in
    one function and a fifth one added later would not know to ask.  Two of the
    four are the arming, which stops one being scheduled; the other two are the
    running, which stops one that already was.
    """
    source = EDITOR_SOURCE
    assert 'def _climb_owns_the_release():' in source
    # Armed, and about to run: both are guarded.
    for after in ('def _arm_gfn_settle(forced=False):',
                  'def _arm_gfn_minimise():',
                  'def _gfn_settle_now():',
                  'def _minimise_now():',
                  'def _arm_gfn_restart():'):
        body = source.split(after, 1)[1].split('\n    def ', 1)[0]
        assert '_climb_owns_the_release()' in body, after
    # And an optimisation interrupted by the grab is stood down rather than
    # left armed to come back over what the climb finds.
    assert "'The climb has this release." in source
    # Dynamik Opt itself is untouched: the follow under the hand is not gated
    # on the climb at all, only the downhill step at the release is.
    follow = source.split('def _gfn_follow_step(', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_climb_owns_the_release' not in follow


def test_the_climb_stands_still_while_a_hand_is_down():
    """And is told at the grab, not at the first drag-follow.

    A follow message is a tenth of a second away and the climb makes a step
    every ten milliseconds, so waiting for one is ten steps away from the hand
    before anything has told the climb a hand arrived.  The grab is the
    earliest word the page says.
    """
    source = EDITOR_SOURCE
    grab = source.split("if verb == 'gfngrab':", 1)[1].split(
        "if verb == 'gfnfree':", 1)[0]
    assert "state['climb_hand_down'] = True" in grab
    assert "state['climb_was'] = (state.get('climb_showing')" in grab
    # Cleared by the hand going, whatever else did or did not arrive.
    free = source.split("if verb == 'gfnfree':", 1)[1].split('\n        if ',
                                                             1)[0]
    assert "state.pop('climb_hand_down', None)" in free
    # The loop stands still rather than stepping.
    assert "if state.get('climb_hand_down'):" in source
    # And the structure is put down before the climb is let go, so the first
    # thing the climb does is take it rather than step once from where it was
    # standing.  The order is the whole of it, so the order is what is read.
    handover = source.split("state['climb_hand'] = {", 1)[1].split(
        '\n        elif ', 1)[0]
    assert "state.pop('climb_hand_down', None)" in handover
    # A hand that lands while the climb is reading what it reached -- xtb's
    # own Hessian is a third of a second, and a hand fits inside that -- keeps
    # what it made.  Writing the climb's geometry over an edit that came after
    # it is the one thing an editor may never do.
    assert "interrupted = state.pop('climb_hand', None) is not None" in source
    assert 'if rows and not interrupted:' in source
    assert 'You moved the structure while it was finishing' in source


def test_the_picture_is_fed_frames_and_never_the_coordinate_box():
    """Every write to the box rebuilds the viewer from nothing.

    A hundred steps a second through the box would be a hundred rebuilds a
    second, which is not a picture of a climb but a slideshow of teardowns.
    The frames go down the channel the optimiser already uses, and the box is
    written once at the end.
    """
    source = EDITOR_SOURCE
    # Through the same helper every other writer uses, so the climb is played
    # at the pace the slider is set to rather than at whatever the page last
    # heard -- and asked at the write, not at the answer, whether the run it
    # was prepared for is still the current one.
    assert 'submit_gfn_frame.value = _frame_payload(r, **fields)' in source
    assert 'if not _frame_run_is_current(r):' in source
    # The window starts where the previous one started, so every frame is sent
    # twice and a read that lands between two writes misses nothing.
    assert 'start, sent[0], sent[1], sent[2] = sent[0], sent[1], len(walked), now' \
        in source
    # And it is held apart, because the climb makes frames far faster than any
    # page is asked to look.
    assert 'now - sent[2] < 0.04' in source
    # The box, once, at the end, with the flag that says the picture has it.
    assert "'Climbed towards a transition state'), drawn=True)" in source
