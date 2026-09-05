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

And then the part that turned out to matter more.  Measured over sixteen hand
drags whose structures are whole, the climb aimed by the hand reaches the
reaction the hand meant 9 times, the same climb told to ignore the hand and
follow the softest mode reaches it 9 times, and ORCA's OptTS reaches it 7 --
and never the same nine, because each one sets off along a different
direction.  Tried cheapest first, with what each one reached checked against
the pair of atoms the hand was holding, they reach **13 of the 16**.  The
tests from ``test_the_held_stretch_...`` down are about that ladder: that the
check is a different question from counting imaginary modes, that every rung
is handed the structure the hand made rather than the last rung's answer, and
that a rung that was right is never followed by one that costs a minute.
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

#: The last rung.  Only the end-to-end demonstration needs it:
#: everything else about the hand-over is tested with a stand-in, because what
#: is being tested there is which structure is passed on and when, not what
#: ORCA does with it.
_needs_orca = pytest.mark.skipif(climb._saddle.find_orca() is None,
                                 reason='ORCA not installed')

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


def test_the_climbed_mode_is_looked_for_at_the_soft_end_and_nowhere_else():
    """Because a Bofill update can grow an eigenvalue no surface has.

    Measured while a climb was running, against a Hessian recomputed exactly
    at the same geometry: a carried curvature of -1379 where the exact one
    says +0.34, on a glycylglycine whose hydroxyl proton had been dragged
    towards the other oxygen.  A spurious mode like that is self-consistent --
    the overlap test reports 1.00 against it for ever, because Bofill keeps
    reinforcing what it invented -- so the guard cannot be the overlap.  It
    has to be the mode's *place in the ordering*, which is what every
    implementation of Baker's mode following uses: geomeTRIC hard-codes the
    climbed mode as the lowest and optking refuses anything else outright.

    Here the reference vector is deliberately the stiffest mode in the list,
    which is what a hand dragging a proton along an O-H bond looks like to an
    unrestricted search.  With the window it is not chosen.
    """
    curvature = np.diag([-0.4, 0.2, 0.6, 1.2, 3.0, 9.0, 40.0])
    gradient = np.array([0.05, 0.01, 0.03, -0.02, 0.01, 0.02, 0.04])
    stiffest = np.zeros(7)
    stiffest[6] = 1.0

    # Unrestricted -- which is what this did before it was measured -- the
    # stiffest mode is followed and the step goes uphill along it.
    wide = climb.partitioned_step(curvature, gradient, stiffest, within=0)
    assert wide['curvature'] == pytest.approx(40.0)
    assert wide['step'][6] * gradient[6] > 0, wide['step']

    # With the window, the search never reaches it: the best of the softest
    # five is taken, and the mode climbed has a curvature a molecule has.
    narrow = climb.partitioned_step(curvature, gradient, stiffest)
    assert narrow['curvature'] <= curvature[climb.AIM_WITHIN - 1,
                                            climb.AIM_WITHIN - 1]
    assert narrow['step'][6] * gradient[6] < 0, narrow['step']
    # And the window is the one the module says it is.
    assert climb.AIM_WITHIN == 5


@_needs_xtb
def test_the_hessian_is_computed_again_rather_than_updated_for_ever():
    """ORCA's ``Recalc_Hess``, which is what the button next door is given.

    :mod:`saddle` writes ``Recalc_Hess 5`` into every OptTS it runs, and
    pysisyphus says the same thing in its own documentation: "when the Hessian
    for the chosen computational method is reasonably cheap it is a good idea
    to recalculate it periodically; between recalculations it's updated using
    the Bofill update".  What it defends against was measured rather than
    taken on trust -- see :data:`climb.HESSIAN_EVERY` for the drift numbers,
    of which the shortest is that a Bofill Hessian carries one to two negative
    eigenvalues the exact Hessian at the same geometry does not have.

    This pins the bookkeeping rather than the chemistry: that a refresh
    happens on the interval, that it is not repeated for a step that was
    refused, and that a climb cannot spend Hessians for ever.  What the
    refresh buys is measured in
    :func:`test_the_hand_names_the_mode_and_that_is_what_helps`, whose three
    drags now converge in 23, 22 and 23 steps against 42, 39 and 37.
    """
    assert climb.HESSIAN_EVERY == 20
    assert climb.MOST_HESSIANS == 8

    walk = climb.Climb(_ESTIMATE, 'gfn2')
    try:
        # Counted rather than timed: a Hessian is 6N gradients, and on this
        # box that is 0.55 s at sixteen atoms against 9 ms for a step -- a
        # number that moves with the machine, where the count does not.
        walk.hessian = np.eye(3 * len(walk.symbols))
        walk.energy, walk.gradient = walk._measure(walk.bohr)
        walk.steps = climb.HESSIAN_EVERY
        walk.hessians = climb.MOST_HESSIANS
        before = int(walk.engine.calls)
        walk.step()
        # Out of Hessians: one gradient for the step and nothing more.
        assert int(walk.engine.calls) - before <= 2, walk.engine.calls

        walk.hessians = 0
        walk.steps = climb.HESSIAN_EVERY
        walk._exact_at = -1
        before = int(walk.engine.calls)
        walk.step()
        spent = int(walk.engine.calls) - before
        # Central differences: two gradients per Cartesian, so 6N of them,
        # and then the step's own -- unless the step was refused, which is
        # why this is a floor rather than an equality.
        assert spent >= 6 * len(walk.symbols), spent
        assert walk.hessians == 1
    finally:
        walk.close()


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

    Told, all three reach the Diels-Alder saddle -- 23, 22 and 23 steps for
    drags of 0.80, 0.95 and 1.10 A -- at 2.315 A on both forming bonds, one
    imaginary mode near -394 cm-1, within 0.002 A RMSD of what ORCA reaches.
    Not told, the shallower two climb the *lowest* mode instead, which is the
    two fragments rocking against each other, and walk back down to the
    van-der-Waals complex 0.43 A away with no imaginary mode at all.

    The deepest drag is the exception and it is worth keeping rather than
    hiding: at 1.10 A the ethene has already been placed at the separation the
    saddle has, so the lowest mode *is* the reaction mode and an unaimed climb
    reaches the saddle too -- in 41 steps against 23.  It did not before the
    exact Hessian was recomputed along the way, which is a second measurement
    of what that recomputation is worth.  So the pull is tested at the depths
    where the answer is in doubt, which is where a hand is actually needed.

    So what a pull contributes is not force.  It is the answer to "which
    reaction", and the climb reads it as the Hessian eigenvector that most
    resembles the direction the structure was dragged -- among the softest
    few, for the reason :data:`climb.AIM_WITHIN` gives.  That is Baker's mode
    following with the reference mode seeded by the user instead of by an
    index, which is the same thing ORCA offers as ``TS_Mode {B 0 1}`` and a
    great deal easier to say with a mouse.
    """
    orca = _where(_ORCA_SADDLE)
    for depth in (0.80, 0.95):
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

        # And the same structure, unaimed, does not get there: at these two
        # depths the lowest mode is the fragments rocking, not the reaction.
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

    So the hand is never part of the climb.  A grab ends it, through the same
    interrupt that ends a minimisation; letting go starts it again from the
    structure that was made, aimed along the way it was made.  Both halves of
    that are in the editor, and this is the wiring that carries them.
    """
    source = EDITOR_SOURCE
    assert 'submit_climb_btn = widgets.ToggleButton(' in source
    assert "description='Climb to TS'" in source
    assert 'def on_submit_climb(change=None):' in source
    assert "submit_climb_btn.observe(on_submit_climb, names='value')" in source
    # Where the structure stood when the hand arrived, which is the other end
    # of the direction the drag names.  The climb's own last frame while one
    # is really walking, and what is in the box otherwise: a hand can arrive
    # before the first Hessian has finished, and it can arrive when nothing is
    # climbing at all -- and that frame outlives the run that made it, so read
    # unconditionally it is the far end of a direction nobody pointed in.
    grab = source.split("if verb == 'gfngrab':", 1)[1].split(
        "if verb == 'gfnfree':", 1)[0]
    assert "(state.get('climb_showing') if climbing else None)" in grab
    assert 'or _current_xyz())' in grab
    # And it is read once, by the run that the release starts -- there is no
    # hand-over of its own any more, because the release path is the
    # minimisation's and it reads the structure out of the box like everything
    # else.
    climbing = source.split('def _climb_now():', 1)[1].split(
        '\n    def ', 1)[0]
    assert "aimed_from = state.pop('climb_was', None)" in climbing
    assert 'aimed_from=aimed_from,' in climbing
    assert 'constraints' not in climbing, 'the pull is never part of the climb'
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
    for the hand, remember where it stood when the hand arrived, start again
    from the structure the hand made, and get to the saddle.

    Measured: from the bare complex it goes nowhere; interrupted by the drag it
    converges to 2.315 A on both forming bonds with one imaginary mode near
    -394, and that is what is left in the coordinate box.
    """
    import time as clock

    pytest.importorskip('ipywidgets')

    part, state, box = _an_editor(tmp_path)

    part.submit_climb_btn.value = True
    began = clock.time()
    while not state.get('climb_showing') and clock.time() - began < 60:
        clock.sleep(0.02)
    assert state.get('climb_run') is not None, 'the climb never started'

    # The hand: down, then let go somewhere else.  The page names the atoms it
    # is holding in the comment line, which is what tells the editor a hand has
    # arrived at all when no grab came with it.
    moved = _dragged(0.95).splitlines()
    part.submit_manip_sync.value = '\n'.join(
        [moved[0], 'DELFIN drag-follow held=10,11'] + moved[2:])
    assert state.get('climb_was'), 'the climb did not remember where it stood'
    assert state.get('climb_run') is None, 'the hand did not stop the climb'
    part.submit_manip_sync.value = '\n'.join(
        [moved[0], 'DELFIN drag-end'] + moved[2:])

    # It comes back through the same restart a minimisation comes back through,
    # after the same wait, so it is not running the moment the hand lets go.
    _quiet(part, state)
    assert state.get('climb_run') is None, 'the climb never finished'

    one, two = _forming(_where(box.value))
    assert 2.25 < one < 2.40 and abs(one - two) < 0.03, (one, two)
    assert _rmsd(_where(box.value), _where(_ORCA_SADDLE)) < 0.05
    said = part.mol_status.value
    assert 'transition state' in said, said
    assert 'cm-1' in said, said
    part.submit_climb_btn.value = False
    _quiet(part, state)



def _an_editor(room):
    """A real editor over the van-der-Waals complex, driven from the kernel.

    The host is the Submit tab's, in the one respect these tests depend on:
    a write to the coordinate box that the viewer has already drawn still
    refreshes ``current_xyz_for_copy``, which is what ``_current_xyz`` answers
    with.  A host that skips it leaves every optimiser -- the minimisation as
    much as the climb -- starting from the structure the session opened on
    rather than from the one the hand just made, and both are then measured
    against a molecule nobody is looking at.
    """
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

    def update_view(*_a, **_k):
        raw = (box.value or '').strip()
        if not raw:
            return
        rows = [one for one in raw.split('\n') if one.strip()]
        body = rows[2:] if rows and rows[0].strip().isdigit() else rows
        state['current_xyz_for_copy'] = {
            'content': f'{len(body)}\nEdited in DELFIN viewer\n'
                       + '\n'.join(body)}
        state['manip_inflight'] = False

    part = structure_editor.build(
        ctx, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda call, *a, **k: call(*a, **k),
        update_view=update_view,
        get_smiles_charge=lambda *a, **k: None)
    box.observe(lambda _change: update_view(), names='value')
    update_view()
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


def _quiet(part, state, seconds=300.0):
    """Wait until nothing is armed and nothing is walking, whichever it was.

    One waiter for both columns, because there is one release path: a release
    arms a wait, the wait starts a walk, and the walk ends.  Asked for either
    optimiser by name it would be a different question for each, which is the
    shape the feature is being taken out of.
    """
    import time as clock

    began = clock.time()
    clock.sleep(2.0)                    # longer than the arming delay
    while clock.time() - began < seconds:
        busy = (state.get('climb_run') is not None
                or state.get('optimize_run') is not None
                or state.get('gfn_settle_busy')
                or state.get('gfn_restart_armed')
                or state.get('gfn_minimise_armed')
                or state.get('gfn_settle_armed')
                or part.submit_optimize_btn.value)
        if not busy:
            break
        clock.sleep(0.15)
    clock.sleep(1.0)
    return clock.time()


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

    A grab now ends the climb the way it ends a minimisation, through the one
    interrupt both go through, so there is nothing left to stand still.  Two
    things are pinned here, and neither is a frame count -- how many follow
    rounds fit inside a drag depends on what else the machine is doing, and a
    number that moves with the load is not a fact about this code.

    What is exact is *whose* frames those are.  During the drag every frame on
    the channel must belong to the run the hand being followed claimed, and
    none to the climb's; and the structure the drag leaves must be the one the
    hand put there.  Measured either way, with the climb off and on: 2.450 and
    2.450 A, from the same six follow answers at 3.201, 3.049, 2.898, 2.748,
    2.598, 2.450.
    """
    quiet, quiet_state, box = _an_editor(tmp_path / 'off')
    quiet.submit_relax_btn.value = True
    seen = _frames_seen(quiet)
    _drag_towards_the_saddle(quiet, box)
    without = [run for run, count in seen if count]
    plain = _forming(_where(box.value))
    assert without, 'Dynamik Opt drew nothing at all, so nothing is being read'
    assert all(run == quiet_state.get('gfn_follow_run') for run in without), (
        quiet_state.get('gfn_follow_run'), without)
    _quiet(quiet, quiet_state)

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
    assert all(run == state.get('gfn_follow_run') for run in withal), (
        state.get('gfn_follow_run'), withal)
    # And the structure the drag leaves is the hand's, not the climb's.
    assert abs(climbed[0] - plain[0]) < 0.05, (plain, climbed)
    assert abs(climbed[1] - plain[1]) < 0.05, (plain, climbed)
    assert 2.3 < climbed[0] < 2.6, climbed
    busy.submit_climb_btn.value = False
    _quiet(busy, state)


@_needs_xtb
def test_the_release_is_one_path_and_the_toggle_is_the_only_difference(tmp_path):
    """"Verhalten im Viewer von TS wie Dynamik Opt, nur dass man es beim
    Loslassen hin zu TS bringt, nicht Minimum."

    The same gesture, twice over, with Dynamik Opt and Auto on both times and
    Climb to TS the only thing that differs.  Everything the user can see
    before the release has to be the same and only where it ends may differ,
    because there is now one release path with the optimiser as its one
    parameter.

    Measured on this complex, and the drag half of it came out identical to
    three decimals: the same six follow answers at 3.201, 3.049, 2.898, 2.748,
    2.598, 2.450 A, the same geometry left behind at 2.450/2.445, the same
    0.903 A of hand.  Then the release parts them -- 3.304/3.302 A and a
    minimum with the toggle up, 2.315/2.315 A and one mode at -394 cm-1 with
    it down.

    Two things are deliberately not asserted.  Frame counts move with the
    machine's load, and how many steps a climb takes depends on where the drag
    happened to stop; what is pinned is that both columns drew, that the drag
    reached the same place in both, and that the ends differ in the one way
    they are meant to.
    """
    plain, plain_state, plain_box = _an_editor(tmp_path / 'off')
    plain.submit_relax_btn.value = True
    plain_seen = _frames_seen(plain)
    _drag_towards_the_saddle(plain, plain_box)
    plain_drag = _forming(_where(plain_box.value))
    plain_during = [run for run, count in plain_seen if count]
    _quiet(plain, plain_state)
    fell = _forming(_where(plain_box.value))

    part, state, box = _an_editor(tmp_path / 'on')
    part.submit_relax_btn.value = True
    part.submit_climb_btn.value = True
    began = time.time()
    while not state.get('climb_showing') and time.time() - began < 120:
        time.sleep(0.02)
    seen = _frames_seen(part)
    _drag_towards_the_saddle(part, box)
    climb_drag = _forming(_where(box.value))
    climb_during = [run for run, count in seen if count]
    _quiet(part, state)
    reached = _forming(_where(box.value))

    # Before the release the two are the same drag.
    assert plain_during and climb_during, (plain_during, climb_during)
    assert abs(plain_drag[0] - climb_drag[0]) < 0.05, (plain_drag, climb_drag)
    assert abs(plain_drag[1] - climb_drag[1]) < 0.05, (plain_drag, climb_drag)

    # After it they part, and only there.  Without the climb, back to the
    # complex it started from.
    assert fell[0] > 3.0 and fell[1] > 3.0, fell
    # With it, the saddle -- and it says so.
    assert 2.25 < reached[0] < 2.40, reached
    assert abs(reached[0] - reached[1]) < 0.03, reached
    assert _rmsd(_where(box.value), _where(_ORCA_SADDLE)) < 0.05
    said = part.mol_status.value
    assert 'transition state' in said, said
    assert 'cm-1' in said, said

    # The toggle is a mode, so it is still down and the next release still
    # goes up.  Lifted on convergence -- which it used to do -- the very next
    # drag walked the structure back to 3.353 A, which is "ich kann es nicht
    # beeinflussen": you can point it at a saddle exactly once.
    assert part.submit_climb_btn.value is True
    part.submit_climb_btn.value = False
    _quiet(part, state)


@_needs_xtb
def test_auto_carries_the_climb_on_with_nothing_pressed_in_between(tmp_path):
    """"Auch hier braucht es ja dann Auto Mode."

    Auto is what makes letting go of an atom finish the job, and it had one
    destination.  Here it is asked twice on one editor with no press of
    anything between the two gestures: drag, let go, arrive; drag again, let
    go again, arrive again.

    Measured: 2.315/2.315 A with one mode at -394 cm-1 after the first, and
    2.318/2.318 A at -390 after a second drag that pushed the two halves 0.30
    A further together.  The second climb is the whole of the test -- it was
    not possible at all before, because the toggle had lifted itself and the
    release had gone back to being a minimisation.
    """
    part, state, box = _an_editor(tmp_path / 'auto')
    part.submit_relax_btn.value = True
    assert part.submit_auto_btn.value is True, 'Auto is the default'
    part.submit_climb_btn.value = True
    began = time.time()
    while not state.get('climb_showing') and time.time() - began < 120:
        time.sleep(0.02)

    _drag_towards_the_saddle(part, box)
    _quiet(part, state)
    first = _forming(_where(box.value))
    assert 2.25 < first[0] < 2.40, first

    # Nothing is pressed here.  The switches stand exactly as the user left
    # them, which is the point.
    assert part.submit_climb_btn.value is True
    assert part.submit_auto_btn.value is True
    seen = _frames_seen(part)
    _drag_towards_the_saddle(part, box, depth=0.30)
    _quiet(part, state)
    second = _forming(_where(box.value))

    # It walked again, and it drew again while it walked.
    drawn = [(run, count) for run, count in seen if count]
    assert drawn, 'the second release drew nothing at all'
    assert any(run for run, _ in drawn
               if run == state.get('climb_frame_run')), (
        state.get('climb_frame_run'), drawn)
    # And it arrived at a transition state again rather than at a minimum.
    assert 2.25 < second[0] < 2.45, second
    assert abs(second[0] - second[1]) < 0.03, second
    said = part.mol_status.value
    assert 'transition state' in said, said
    part.submit_climb_btn.value = False
    _quiet(part, state)


@_needs_xtb
def test_the_climb_never_writes_under_a_run_number_a_drag_has_moved_past(
        tmp_path):
    """"Ich sehe nicht das Update von Climb to TS im Viewer."

    The run number is how a writer says it is the current one, and a drag
    moves it: the hand being followed takes a number, and whatever answers the
    release takes another.  The climb claimed one when its toggle went down
    and held it across all of that -- so by the time it had anything to draw
    it was two behind, it failed its own guard on every write, and it threw
    away every frame it made.  Silently, and correctly: it really was not the
    current run.

    It claims where the minimisation claims, through the same helper and at
    the same moment: when a run begins.  Measured on this gesture, the climb
    starts on run 1, the drag leaves run 3 behind it, and the climb draws
    under run 4 -- so a climb still holding the toggle-time number would have
    discarded every frame of its path.
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
    _quiet(part, state)

    drawn = [(run, count) for run, count in seen if count]
    assert drawn, 'the climb drew nothing after the release'
    # The drag moved the counter past the number the climb started with...
    assert after_the_drag > at_the_toggle, (at_the_toggle, after_the_drag)
    # ...and the climb took a fresh one rather than writing under a stale one.
    assert all(run > at_the_toggle for run, _ in drawn), (at_the_toggle, drawn)
    # Every frame of it would have been thrown away by the old guard.
    assert sum(count for run, count in drawn if run == at_the_toggle) == 0
    part.submit_climb_btn.value = False
    _quiet(part, state)


def test_the_release_has_one_path_and_the_optimiser_is_its_parameter():
    """Ein Pfad, zwei Optimierer.

    Every defect this button has had came from being a second implementation
    of something that already existed: it stepped on a stale geometry during
    the drag, it fought Dynamik Opt for the release, it claimed its run once
    and discarded its own frames, and it never got an auto mode because
    nothing on the auto path knew it was there.  Each was patched on its own.

    So what is read here is not that the two agree -- it is that there is only
    one of them to agree with.  A release is decided in _after_release and
    nowhere else, it is armed by one function, it lands on one function, and
    that one function is the only place in the file where the two optimisers
    are named apart.
    """
    source = EDITOR_SOURCE

    # One decision, and it does not name either optimiser.
    release = source.split('def _after_release():', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_arm_gfn_optimise()' in release
    assert '_arm_gfn_settle()' in release
    assert '_climb_now' not in release, 'the release picks a direction, not a run'
    assert 'submit_optimize_btn' not in release

    # One arming, and it does not name either optimiser either.
    arming = source.split('def _arm_gfn_optimise(asked=False):', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_climb' not in arming
    assert 'submit_optimize_btn' not in arming

    # And one landing, which is where -- and only where -- they part.
    landing = source.split('def _optimise_now():', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_climb_owns_the_release()' in landing
    assert '_climb_now()' in landing
    assert 'submit_optimize_btn.value = True' in landing

    # The toggle is the direction and nothing else: a mode, read as a mode.
    owns = source.split('def _climb_owns_the_release():', 1)[1].split(
        '\n    def ', 1)[0]
    assert 'return bool(submit_climb_btn.value)' in owns

    # Both optimisers are interrupted by one function, and both are brought
    # back by one function.  Named apart in either, a hand would once again be
    # able to stop one of them and not the other.
    stop = source.split('def _interrupt_gfn():', 1)[1].split('\n    def ', 1)[0]
    assert "state['climb_run'] = None" in stop
    assert "state['climb_interrupted'] = True" in stop
    assert "state['optimize_interrupted'] = token" in stop
    back = source.split('def _restart_gfn():', 1)[1].split('\n    def ', 1)[0]
    assert "state.pop('climb_interrupted', False)" in back
    assert '_climb_now()' in back
    assert "state.pop('optimize_interrupted', None)" in back

    # Dynamik Opt itself is untouched: the follow under the hand is not gated
    # on the climb at all, which is why the drag measures the same either way.
    # Cut at the last line of the function rather than at the next def: the
    # bounds the two walks are held to are written between the two, so a split
    # on "def" reads them as part of this one.
    follow = source.split('def _gfn_follow_step(', 1)[1].split(
        "_start_background(_work, 'The relaxation under the hand')", 1)[0]
    assert '_climb_owns_the_release' not in follow
    assert 'climb' not in follow


def test_the_stop_has_one_path_too_and_the_optimiser_is_its_parameter():
    """"Stop ist wie bei Dynamik Opt aus machen: das Frame, was man sieht."

    The release was unified and the Stop was not, so the climb still had a
    stop of its own -- and a stop of its own meant it did none of the three
    things a Stop does.  It never told the page, so the page walked out the
    rest of the trajectory: measured, a Stop at frame 13 of 117 at twelve
    frames a second drew 509 more times, frames 14 through 116.  It never
    cleared the halt mark when it started, so a halt could be swallowed by a
    grab that happened before it.  And it wrote what the walk had reached into
    the box, where the minimisation writes the frame the picture stopped on.

    Read here as one path rather than as two that agree: one function tells
    the page, one function cuts a path at the frame on screen, and both
    walkers call both.
    """
    source = EDITOR_SOURCE

    # One halt, and it names no optimiser.
    halt = source.split('def _halt_the_frames(run_id):', 1)[1].split(
        '\n    def ', 1)[0]
    assert "state.get('gfn_halt_sent')" in halt
    assert 'submit_gfn_frame' in halt
    assert '_frame_payload(run_id, halt=1, frames=[])' in halt, \
        'the halt is a write on the frame channel and is paced like one'
    # Neither optimiser is named in what it does -- only in what it says about
    # why it exists.  Read over the code alone, so a comment can go on telling
    # the story.
    doing = halt.split('"""', 2)[-1]
    assert 'climb' not in doing and 'optimize_run' not in doing

    # Both loops ask the same question and get the same answer written.
    optimising = source.split('def on_submit_optimize(change=None, '
                              'every_frame=False)', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_halt_the_frames(run_id)' in optimising
    climbing = source.split('def _climb_now():', 1)[1].split(
        '\n    def on_submit_', 1)[0]
    assert '_halt_the_frames(run)' in climbing
    # The loop is the module's now -- three searches rather than one, and one
    # Stop for all three -- so what the editor holds is the question, handed
    # over to be asked by every rung of it.
    assert 'should_stop=_stopped,' in climbing

    # And both claim their run the same way, which is what clears the halt
    # mark and the stale frame number for the run that is beginning.  The
    # frame number goes with the walk it was counted along, so both are
    # dropped together and there is one place that does it.
    assert "state['gfn_halt_sent'] = False" in climbing
    assert '_forget_the_shown_frame()' in climbing
    forget = EDITOR_SOURCE.split('def _forget_the_shown_frame')[1].split(
        '\n    def ')[0]
    assert "state.pop('gfn_shown_frame', None)" in forget
    assert "state.pop('gfn_shown_run', None)" in forget

    # One cutter, and the box is written by whichever of the worker and the
    # page's report arrives second -- they race, and a climb step is ten
    # milliseconds against an xtb round of seconds.
    assert 'def _frame_as_xyz(source, walked, shown, comment):' in source
    assert '_the_picture_stopped_here(' in climbing
    assert '_the_picture_stopped_here(' in optimising
    landing = source.split('def _land_the_stopped_frame():', 1)[1].split(
        '\n    def ', 1)[0]
    # The count the page reported, and only where it counts along the walk the
    # stopped path belongs to: the page names both, and a count that came from
    # another walk indexes this one's path perfectly well and wrongly.
    assert "_the_shown_frame_of(held.get('run'))" in landing
    assert '_write_coords(text, drawn=True)' in landing
    cmd = source.split('def on_submit_cmd(', 1)[1].split('\n    def ', 1)[0]
    assert '_land_the_stopped_frame()' in cmd

    # A Stop pays for no Hessian either -- a verdict is about where the walk
    # stands, and where the walk stands is not where the user stopped it.
    # Measured rather than read:
    # tests/test_gfn_methods_in_the_viewer.py drives a Stop against a climb
    # whose verdict() raises, and it passes.
    assert 'switched_off[0] = True' in climbing


def test_the_toggle_is_a_mode_and_never_lifts_itself():
    """Auto reads it, so a switch that turns itself off turns Auto off with it.

    Optimize lifts itself when it converges, and that is right for a run: one
    press, one minimisation.  Climb to TS is not a run.  It is what Auto
    consults to decide which way a release walks, and it lifted itself the
    moment a climb converged -- so the next drag was answered by a
    minimisation, which walked the structure straight back off the saddle it
    had just found.  Measured on the van-der-Waals complex: 2.316 A after the
    first release, 3.353 A after the second.

    Dynamik Opt does not turn itself off and neither does Auto.  Nor does this.
    """
    source = EDITOR_SOURCE
    climbing = source.split('def _climb_now():', 1)[1].split(
        '\n    def ', 1)[0]
    assert 'submit_climb_btn.value = False' not in climbing.split(
        'def _done():', 1)[1], 'the climb lifted its own toggle again'
    # It comes up for one reason: the method on screen cannot climb, which is
    # a press that could never have worked.
    handler = source.split('def on_submit_climb(change=None):', 1)[1].split(
        '\n    def ', 1)[0]
    assert handler.count('submit_climb_btn.value = False') == 1
    assert '_climb_can_run()' in handler
    # And standing a run down after a hand does not lift it either -- that
    # path lifts the Optimize switches, which are runs.
    stand = source.split('def _stand_down_after_interrupt(note=None):',
                         1)[1].split('\n    def ', 1)[0]
    assert 'submit_climb_btn' not in stand
    assert "state.pop('climb_interrupted', False)" in stand


def test_a_stop_keeps_the_frame_on_screen_and_a_hand_keeps_nothing():
    """The same two answers Optimize gives, told apart the same way.

    Pressing the toggle off is a Stop: what is kept is the frame you were
    looking at, and it lands in the box -- and the run number moves on first,
    so the frames the stopped walk still had in hand cannot play out over it
    afterwards.  A hand is not a Stop: the user has made a structure since,
    and writing what the climb reached over it is the one thing an editor may
    never do.  Which of the two it was is marked by token, on the walk itself,
    rather than read off a flag something else may have cleared in between.
    """
    source = EDITOR_SOURCE
    stop = source.split('def _interrupt_gfn():', 1)[1].split(
        '\n    def ', 1)[0]
    assert "state['climb_cut'] = state['climb_run']" in stop
    handler = source.split('def on_submit_climb(change=None):', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_claim_the_frame_run()' in handler, 'a Stop must refuse what is in flight'
    climbing = source.split('def _climb_now():', 1)[1].split(
        '\n    def on_submit_', 1)[0]
    assert "cut_by_a_hand = state.get('climb_cut') is token" in climbing
    assert "if state.get('climb_cut') is token:" in climbing
    # A Stop lands the frame the picture stopped on, which is the same thing
    # the minimisation's Stop lands and through the same function.
    assert '_the_picture_stopped_here(' in climbing
    # And a toggle switched off in the middle of a drag is a hand too, even
    # though no grab ever marked this walk as cut: by then it was already
    # stopped, so there was nothing left for the grab to cut.
    assert "holding = bool(state.get('gfn_follow'))" in climbing
    assert 'if rows and not holding:' in climbing
    assert 'You moved the structure while it was finishing' in climbing


def test_a_climb_that_is_going_nowhere_stops_and_says_so():
    """The minimisation has had a bound since it was written; this had none.

    Measured on the van-der-Waals complex with an aim that named nothing: 24616
    steps and still walking, a gradient every ten milliseconds for as long as
    the toggle stayed down.  Nothing ended it, because a climb converges or it
    does not and there was no third answer.

    There were then two third answers, and they disagreed: four hundred in the
    editor, picked as an order of magnitude above every climb that had
    finished at the time, and a hundred in the module, measured afterwards
    over twenty-one hand drags -- every climb that reached the reaction the
    hand pointed at arrived in 22 to 65 steps, and the ones that ran longer
    ran on to between 202 and 361 and were wrong at the end of all of them.
    The measured one is kept, and it is kept in one place: the module that
    walks the steps.  The other three hundred bought nothing but a longer wait
    before the next search was tried -- measured on a Diels-Alder drag whose
    first climb goes to a minimum, 211 steps and 83.8 s at four hundred
    against 61 steps and 21.8 s at a hundred, for the same wrong answer.
    """
    source = EDITOR_SOURCE
    assert climb.CLIMB_CEILING == 100
    assert '_CLIMB_STEPS' not in source, 'one ceiling, in the module'
    climbing = source.split('def _climb_now():', 1)[1].split(
        '\n    def on_submit_', 1)[0]
    assert 'max_steps=_climb.CLIMB_CEILING,' in climbing
    # And what it reached is still named, because a saddle search does not
    # fail -- it succeeds at arriving somewhere.  The sentence is the module's:
    # it is the one that knows which searches were tried and what each reached.
    assert "said = str(got.get('status') or '')" in climbing
    # With whatever the walk that handed this climb its start had to say in
    # front of it, so a barrier and the saddle it belongs to arrive together.
    assert 'lines = walked_said + [said]' in climbing


def test_the_picture_is_fed_frames_and_never_the_coordinate_box():
    """Every write to the box rebuilds the viewer from nothing.

    A hundred steps a second through the box would be a hundred rebuilds a
    second, which is not a picture of a climb but a slideshow of teardowns.
    The frames go down the channel the optimiser already uses -- through the
    same writer, not through a second one that happens to agree -- and the box
    is written once at the end.
    """
    source = EDITOR_SOURCE
    # One writer, and both walkers hand their path to it.
    assert 'def _stream_frames(run_id, frames, *, final=False, follow=False,' \
        in source
    climbing = source.split('def _climb_now():', 1)[1].split(
        '\n    def on_submit_', 1)[0]
    assert '_stream_frames(run, walked, final=final, least_apart=0.04)' \
        in climbing, 'a climb makes frames faster than any page is asked to look'
    # And not marked as a followed hand.  The mark was there to keep the page
    # from abandoning the climb's queue -- the page knew only one switch -- and
    # it cost the speed slider, because a followed hand is paced by how fast
    # xtb answers rather than by any setting.  The page reads the climb's own
    # toggle now, so the climb is streamed exactly as a minimisation is.
    assert 'follow=True' not in climbing, \
        'the climb is a walked path, not a hand being followed'
    optimising = source.split('def _push_frames(frames, final=False):',
                              1)[1].split('\n        #', 1)[0]
    assert '_stream_frames(run_id, frames, final=final)' in optimising

    writer = source.split('def _stream_frames(', 1)[1].split(
        '\n    def ', 1)[0]
    # The window starts where the previous one started, so every frame is sent
    # twice and a read that lands between two writes misses nothing.
    assert "state['gfn_push_start'] = int(state.get('gfn_push_end') or 0)" \
        in writer
    # Asked at the write and not at the answer, so a run that has been
    # replaced cannot draw over what replaced it.
    assert 'if not _frame_run_is_current(run_id):' in writer
    # The mark itself stays, for the writers that really are a followed hand:
    # the drag and the settle behind it have no switch to be read.
    assert "fields['follow'] = 1" in writer

    # The box, once, at the end -- and told whether the picture has it rather
    # than assuming so.  Assumed, a climb whose final frames were refused left
    # the viewer standing on a geometry the box no longer held.
    assert "'Climbed towards a transition '" in source
    assert 'drawn=_frame_run_is_current(run))' in source


# -- the ladder: one press, three tries -------------------------------------

#: A ring pulled open by hand: cyclohexene with C0-C1 dragged from 1.53 to
#: 2.39 A and everything else relaxed around it under GFN2, which is what the
#: editor's follow leaves behind.  It is kept because it is the case that
#: makes the ladder worth building -- the climb *converges* on it, in 84 steps,
#: onto a saddle at -534 cm-1 whose imaginary mode has nothing to do with the
#: bond that was pulled, and ORCA's OptTS from the same structure reaches the
#: retro-Diels-Alder at -394.
_RING_PULLED_OPEN = """16
cyclohexene with C0-C1 held at 2.4 A, the rest relaxed under GFN2
C      -1.488147806351    -0.537649585138    -0.047156919196
C      -0.524460978485     1.637328020803     0.123593624669
C       0.640367137115     1.383616704122    -0.627719776438
C       1.438612632379     0.322243533332    -0.371715260080
C       1.005644493407    -0.706951059580     0.582377782570
C      -0.308430050328    -1.407730372038     0.028238271626
H      -2.118269134393    -0.445955384136     0.819631273920
H      -1.958150864626    -0.366588567961    -0.997419539675
H      -1.263727705510     2.324616297746    -0.249981764156
H      -0.521621450610     1.483261424764     1.187825325833
H       0.842496407260     2.012178426677    -1.486681258414
H       2.293322249921     0.116025962908    -0.998008615585
H       1.764615633845    -1.476975358614     0.716192764571
H       0.747704844024    -0.288828344783     1.554887323158
H      -0.491468464467    -2.234289432240     0.718003717147
H      -0.058486943179    -1.814302265860    -0.952066949949
"""

#: And where that ring stood before the hand touched it, minimised under GFN2.
#: The climb is aimed by the difference between the two, which is what the
#: editor hands it -- and on this drag that matters more than anywhere else in
#: this file.  Aimed by the hand the climb goes to the wrong saddle; aimed at
#: nothing it follows the lowest mode, which on this structure happens to be
#: the retro-Diels-Alder already, and reaches it in 42 steps.  Leaving the
#: aim out here would test a case the editor never produces.
_CYCLOHEXENE = """16
cyclohexene, minimised under GFN2
C           -1.27487075798282       -0.27853481028821       -0.19519539689222
C           -0.89976507622432        1.19310896921734       -0.04260954294814
C            0.56710990115438        1.40677024953639       -0.26013402071706
C            1.45721400480502        0.42892655197901       -0.24664719797261
C            1.10436843117062       -1.00756919860520       -0.00861076010219
C           -0.29392253602956       -1.14407965287449        0.58793672785123
H           -2.29518149448465       -0.44312193976041        0.15693710051218
H           -1.23318757803729       -0.55612624086301       -1.25062596970557
H           -1.46707570112588        1.80332289934208       -0.75033844491857
H           -1.16217438037387        1.54164472951602        0.96248667462100
H            0.87790006082756        2.42982968963257       -0.42455774110749
H            2.50418105459721        0.63711041887047       -0.42213656974589
H            1.15672869739010       -1.54658216552424       -0.96112210924520
H            1.84485815128365       -1.46074100544348        0.65576344721829
H           -0.27812121693153       -0.81435948304502        1.62897619419220
H           -0.60806156003763       -2.18959901168981        0.56987760895903
"""


def test_the_held_stretch_is_the_one_direction_the_check_asks_about():
    """"Make this contact longer", in the metric the modes are written in.

    The check that decides a hand-over is a dot product, so the vector it dots
    against has to be a unit vector in the same space as the mode: the Wilson
    B row of the pair's distance, divided by the root masses because a B row
    differentiates with respect to Cartesians, and normalised.  Nothing here
    needs xtb -- it is the geometry of one direction.
    """
    walk = climb.Climb(_ORCA_SADDLE, 'gfn2')
    try:
        unit = walk.held_stretch((0, 10))
        assert unit is not None
        assert float(np.linalg.norm(unit)) == pytest.approx(1.0)
        rows = unit.reshape(-1, 3)
        # Only the two atoms named move, and they move against each other.
        moving = [i for i in range(len(rows))
                  if np.linalg.norm(rows[i]) > 1e-12]
        assert moving == [0, 10], moving
        assert float(rows[0] @ rows[10]) < 0

        # Mass weighting is not decoration: the same contact between a carbon
        # and a hydrogen puts most of the vector on the hydrogen, because that
        # is the atom that actually moves when the distance changes.
        light = walk.held_stretch((0, 4)).reshape(-1, 3)
        assert np.linalg.norm(light[4]) > 3.0 * np.linalg.norm(light[0])

        # A pair that is not a pair is refused rather than guessed at.
        assert walk.held_stretch((0, 0)) is None
        assert walk.held_stretch((0, 99)) is None
        assert walk.held_stretch(None) is None
    finally:
        walk.close()


@_needs_xtb
def test_one_imaginary_mode_is_not_the_same_question_as_the_right_one():
    """Which is why the verdict grew a second half rather than a better number.

    Counting imaginary modes is a test a wrong answer passes.  Measured over
    twenty-one hand drags, both this climb and ORCA's OptTS converge -- with
    exactly one mode going the wrong way, and reporting success -- onto methyl
    torsions at -48 cm-1 and onto two fragments rocking at -52.  So the
    structure is asked what its imaginary mode *does*, against the pair of
    atoms the hand was holding.

    The Diels-Alder saddle is the same structure either way; only the question
    changes.  Asked about a bond that is forming it scores 0.68, and asked
    about a C-H bond that is doing nothing it scores near zero -- and the gap
    between those is where :data:`climb.IS_THE_REACTION` sits.
    """
    walk = climb.Climb(_ORCA_SADDLE, 'gfn2')
    try:
        forming = walk.verdict(held=(0, 10))
        spectator = walk.verdict(held=(0, 4))
        quiet = walk.verdict()
    finally:
        walk.close()
    # One imaginary mode, and the old half of the verdict says so all three
    # times: it cannot tell these apart, which is the point.
    assert forming['ok'] and forming['count'] == 1, forming
    assert spectator['ok'] and spectator['count'] == 1, spectator
    # The new half separates them.  Margins rather than numbers: xtb under
    # threading is not bit-reproducible, and a threshold is what is being
    # tested anyway.
    assert forming['share'] > 0.5, forming
    assert spectator['share'] < 0.2, spectator
    assert forming['reaction'] is True, forming
    assert spectator['reaction'] is False, spectator
    # And with nothing held there is no second half to report, rather than a
    # guess presented as an answer.
    assert quiet['share'] is None and quiet['reaction'] is None, quiet


@_needs_xtb
def test_xtbs_own_hessian_carries_the_shapes_as_well_as_the_numbers():
    """The check needs a mode, not a frequency, and it is the same subprocess.

    ``xtb --hess`` prints its frequencies and writes the whole Cartesian
    Hessian to a file called ``hessian`` beside its input.  The module used to
    parse the printout, which gives numbers and no shapes; it reads the file,
    which gives both.  Measured over the twenty-one drags the two agree to
    0.11 cm-1 on every mode of every structure and cost the same -- 0.11 s at
    ten atoms, 0.4 at sixteen, 3.7 at thirty-three, 11.8 at fifty.
    """
    walk = climb.Climb(_ORCA_SADDLE, 'gfn2')
    try:
        got = walk.modes_from_xtb()
        assert got is not None
        # Translations and rotations are projected out, so 3N - 6 of them.
        assert got['cm'].shape[0] == 3 * len(walk.symbols) - 6
        assert got['shape'].shape[1] == got['cm'].shape[0]
        # Softest first, and the softest is the reaction coordinate.
        assert list(got['cm']) == sorted(got['cm'])
        assert -450 < got['cm'][0] < -340, got['cm'][:3]
        # The shapes are an orthonormal set, which is what makes an overlap
        # against one of them a number between nought and one.
        overlaps = got['shape'].T @ got['shape']
        assert np.abs(overlaps - np.eye(overlaps.shape[0])).max() < 1e-8
        # And the frequencies are the same object the older accessor returns.
        assert np.abs(np.asarray(walk.frequencies_from_xtb())
                      - got['cm']).max() < 1e-9
    finally:
        walk.close()


@_needs_xtb
def test_every_rung_is_handed_the_hands_structure_and_not_the_last_ones():
    """It looks like throwing work away, and it is measured to be right.

    Over the twelve drags of the twenty-one the aimed climb gets wrong, ORCA
    started from the hand's own guess reaches the reaction on **5** of them
    and ORCA started from where the climb stopped reaches it on **2**.  A
    search that is going wrong does not stop somewhere useful: handing on its
    endpoint loses the retro-Diels-Alder at 2.4 A, the Claisen at 2.0 and the
    salicylaldehyde proton transfer, all three of which the same optimiser
    reaches from the structure the hand made.

    That is also what makes a rung that misses cheap to allow: it costs its
    own wall time and changes nothing about what the next one is given.
    """
    handed = []

    def _pretend(xyz_text, method='gfn2', **kw):
        handed.append(xyz_text)
        return {'ok': False, 'xyz': None, 'seconds': 0.0,
                'status': 'not run in this test'}

    was_orca = climb._saddle.find_orca
    was_opt = climb._saddle.optimise_to_saddle
    climb._saddle.find_orca = lambda: '/nowhere/orca'
    climb._saddle.optimise_to_saddle = _pretend
    try:
        # The complex on its own has no saddle to find, so every check fails
        # and the ladder runs to the end.  Six steps a rung, because what is
        # being tested is which structure is passed on rather than where a
        # climb gets to.
        got = climb.reach_the_reaction(_COMPLEX, 'gfn2', held=(0, 10),
                                       aimed_from=_dragged(0.4), max_steps=6)
    finally:
        climb._saddle.find_orca = was_orca
        climb._saddle.optimise_to_saddle = was_opt

    assert len(handed) == 1, handed
    assert _where(handed[0]) == pytest.approx(_where(_COMPLEX)), \
        'the last rung is given the structure the hand made'
    assert not got['ok']


@_needs_xtb
def test_the_softest_mode_rung_is_skipped_when_it_would_be_the_same_climb():
    """With no gesture to aim along, the first rung already *is* the second.

    :meth:`aim` returns nothing without an *aimed_from*, and a climb with
    nothing to follow takes the lowest mode -- which is exactly what the
    second rung is for.  Running it would be the same climb twice and one more
    Hessian for a certainty.
    """
    calls = []
    was_climb = climb.climb_to_saddle

    def _counted(*args, **kw):
        calls.append(kw.get('aimed_from'))
        return was_climb(*args, **kw)

    was_orca = climb._saddle.find_orca
    climb.climb_to_saddle = _counted
    climb._saddle.find_orca = lambda: None
    try:
        climb.reach_the_reaction(_COMPLEX, 'gfn2', held=(0, 10), max_steps=6)
    finally:
        climb.climb_to_saddle = was_climb
        climb._saddle.find_orca = was_orca
    assert calls == [None], calls


@_needs_xtb
def test_a_climb_that_reached_the_reaction_never_pays_for_another_rung():
    """The common case is the fast one, which is the whole point of the order.

    Cheapest first and ORCA last reaches 15 of the twenty-one drags.  ORCA
    first, with the climb after it, reaches 9 -- started from where ORCA
    stopped the climb rescues 0 of the 12 ORCA gets wrong -- and started from
    the hand's structure instead it would reach the same 15 while paying
    ORCA's median 50.3 s on every press, including the nine the climb answers
    in 0.5 to 36.5 s.  So nothing beyond the first rung may run when the first
    rung was right, and this is what says so: both stand-ins raise if they are
    ever reached.
    """
    def _never(*_, **__):
        raise AssertionError('a later rung ran on a climb that was right')

    calls = []
    was_climb = climb.climb_to_saddle

    def _counted(*args, **kw):
        calls.append(kw.get('aimed_from'))
        return was_climb(*args, **kw)

    was_orca = climb._saddle.optimise_to_saddle
    climb.climb_to_saddle = _counted
    climb._saddle.optimise_to_saddle = _never
    try:
        got = climb.reach_the_reaction(_ESTIMATE, 'gfn2', held=(0, 10),
                                       aimed_from=_COMPLEX)
    finally:
        climb.climb_to_saddle = was_climb
        climb._saddle.optimise_to_saddle = was_orca
    assert len(calls) == 1, 'the softest-mode rung ran on a climb that arrived'
    assert got['ok'] and got['route'] == 'the climb', got.get('status')
    assert got['imaginary']['reaction'] is True, got['imaginary']
    # And the sentence names the contact the way the picture numbers it: the
    # viewer's ``#`` button draws ``String(i)`` over the atoms, counting from
    # nought, so a sentence that counted from one would name other atoms.
    assert 'C0-C10' in got['status'], got['status']
    assert 'reached the reaction you pointed at' in got['status']


@_needs_xtb
def test_with_no_contact_named_the_ladder_stops_at_the_first_rung():
    """Because every later rung would be judged by the test this one passed.

    Without a pair to check against, all that can be asked is how many modes
    go the wrong way -- and that is a test the wrong answer passes: both
    climbs and ORCA converge onto methyl torsions at -48 cm-1 with exactly one
    imaginary mode.  So a retry could only spend a minute to agree, and the
    sentence says which of the two questions was actually answered instead of
    letting a weaker check read like the stronger one.
    """
    def _never(*_, **__):
        raise AssertionError('a retry ran with nothing to judge it by')

    was = climb._saddle.optimise_to_saddle
    climb._saddle.optimise_to_saddle = _never
    try:
        got = climb.reach_the_reaction(_ESTIMATE, 'gfn2', aimed_from=_COMPLEX)
    finally:
        climb._saddle.optimise_to_saddle = was
    assert got['route'] == 'the climb', got.get('status')
    assert got['imaginary']['share'] is None, got['imaginary']
    assert 'Nothing said which contact you meant' in got['status'], \
        got['status']
    # Still called a transition state, because that much was established.
    assert got['ok'] and 'reached a transition state' in got['status']


def test_the_numbers_the_ladder_is_built_on_are_the_ones_that_were_measured():
    """Each of the three is a measurement, and each is written down where it is.

    * :data:`climb.IS_THE_REACTION` -- the structures that are the reaction
      score 0.40 to 0.95 and the torsions and rockings both routes wrongly
      converge onto score 0.00 to 0.07, so any threshold between 0.07 and 0.40
      gives the same answer on all twenty-one and a fifth is the middle of it.
    * :data:`climb.CLIMB_CEILING` -- every climb that reached the reaction
      arrived in 22 to 65 steps and none that ran longer ever did, so a
      hundred cannot cut a climb that is getting there.
    * :data:`climb.FALLBACK_STEPS` -- raising this was tried and dropped.
      Two of ORCA's nine reactions over the twenty-one needed 80 and 120
      optimiser cycles against the 60 it ships with, but both are among the
      five drags whose structure fell apart while the bench built it; over
      the sixteen whole ones every reaction ORCA reaches converges in 13 to
      49 cycles and 60 and 120 give it 7 of 16 either way.  Meanwhile the
      presses where ORCA's answer is no run to the ceiling by definition and
      take a median 63 s at 120 against 25 s at 60.
    """
    import inspect

    assert climb.IS_THE_REACTION == 0.20
    assert 0.07 < climb.IS_THE_REACTION < 0.40
    assert climb.CLIMB_CEILING == 100
    assert climb.CLIMB_CEILING > 65, 'the longest climb that arrived'
    assert climb.FALLBACK_STEPS == 60
    assert climb.FALLBACK_STEPS > 49, (
        'the most cycles any reaction ORCA reaches on a whole structure '
        'needed')
    default = inspect.signature(
        climb._saddle.optimise_to_saddle).parameters['max_steps'].default
    assert climb.FALLBACK_STEPS == default, (
        'pinned to what the optimiser already ships with, and pinned rather '
        'than inherited so this rung knows what it is asking for')


def test_the_ladder_refuses_what_no_engine_here_can_take():
    """The same two refusals the climb alone makes, before anything is started."""
    refused = climb.reach_the_reaction(_ESTIMATE, 'PBE0', held=(0, 10))
    assert not refused['ok'] and refused['route'] is None
    assert 'runs on xtb' in refused['status']
    empty = climb.reach_the_reaction('', 'gfn2', held=(0, 10))
    assert not empty['ok'] and empty['route'] is None


@_needs_xtb
def test_the_ladder_reaches_a_reaction_the_aimed_climb_walks_straight_past():
    """The demonstration the whole thing rests on, on a drag a hand could make.

    Cyclohexene with one ring bond pulled to 2.4 A and the rest relaxed around
    it.  The aimed climb *converges* on this one -- it is not a run that gives
    up -- onto a saddle at about -534 cm-1 whose imaginary mode is 0.00 of the
    bond that was pulled: a real stationary point, a real transition state,
    and the wrong reaction.  Counting imaginary modes calls that a success.

    Asking what the mode does calls it a miss, and the next rung is not ORCA
    but the same climb told to forget the gesture: the softest mode at this
    geometry is 0.94 of the C0-C1 stretch already, and following it reaches
    the retro-Diels-Alder at about -394 cm-1 in 42 steps.  Measured end to
    end, 6.4 s for the ladder, against 13.5 s for ORCA reaching the same
    saddle from the same structure and 3.2 s for a climb that is confidently
    wrong.  That rung is the reason this is three tries rather than two.
    """
    said = []
    was_orca = climb._saddle.find_orca
    # The third rung is taken away rather than allowed to rescue the test: what
    # is being proved is that the second one is enough here, and a run that
    # quietly fell through to ORCA would prove the opposite while passing.
    climb._saddle.find_orca = lambda: None
    try:
        got = climb.reach_the_reaction(_RING_PULLED_OPEN, 'gfn2', held=(0, 1),
                                       aimed_from=_CYCLOHEXENE,
                                       on_route=said.append)
    finally:
        climb._saddle.find_orca = was_orca
    assert got['route'] == 'the softest mode', got.get('status')
    assert got['ok'], got.get('status')
    shape = got['imaginary']
    assert shape['count'] == 1 and shape['share'] > 0.5, shape
    assert shape['modes'][0] < -300, shape
    # The bond the hand pulled is what came apart, rather than something else
    # that happened to be soft.
    assert math.dist(*_where(got['xyz'])[[0, 1]]) > 2.0

    # The user is told the route changed, while it changes rather than after,
    # and told why in the same sentence.
    assert len(said) == 1, said
    assert 'C0-C1' in said[0] and 'softest mode' in said[0], said[0]
    assert 'C0-C1' in got['status'], got['status']
    # And the sentence about the first rung is still in there, because "it
    # went somewhere else" is the reason the second one ran.
    assert 'different reaction from the one you pointed at' in got['status']


#: A proton half transferred: salicylaldehyde with the hydroxyl H dragged
#: across its own hydrogen bond to 1.10 A of the aldehyde oxygen.  This is the
#: one drag of the sixteen whole ones that only ORCA reaches, so it is what
#: the last rung is tested on.
_PROTON_HALF_MOVED = """15
salicylaldehyde with O0-H14 held at 1.10 A, the rest relaxed under GFN2
O       2.955825419586     0.673853493232     0.173427254023
C       1.998907772225    -0.294240634235     0.329980407639
C       0.664816508150    -0.198027998003     0.136885652934
C       0.022314964497     1.017033130549    -0.286793965667
C      -1.312646430845     1.074313535169    -0.466627571877
C      -2.118765421991    -0.081054418519    -0.236257919133
C      -1.579403981838    -1.250692889134     0.161533550317
C      -0.148540304432    -1.413708685350     0.380413263177
O       0.347859341899    -2.472693032314     0.737841158904
H       2.664663222262     1.597379672186    -0.126262849056
H       0.631982986848     1.896107541679    -0.460406733686
H      -1.790736168544     1.988971156377    -0.784660700128
H      -3.185729765118     0.005493157402    -0.390067805214
H      -2.189758488170    -2.124411362691     0.335935059237
H       3.039210345475    -0.418322666348     0.495061198530
"""

#: And where it stood before the hand touched it, minimised under GFN2.
_SALICYLALDEHYDE = """15
salicylaldehyde, minimised under GFN2
O            2.62872278490865        1.69083157857552       -0.13401011557226
C            2.14874100818517        0.60121369096412        0.05857344471600
C            0.72524663873148        0.28272125992550        0.00299308315211
C           -0.16768017606990        1.31329808144758       -0.32144963634084
C           -1.52187372037606        1.08967529935414       -0.39778161613399
C           -2.01740537195299       -0.18560109433387       -0.14824756717333
C           -1.16178375232469       -1.21678142720703        0.17474137862119
C            0.21442244508797       -1.00294037456158        0.25654048012243
O            0.97206085434591       -2.07127920130876        0.58408729152588
H            2.82108750314516       -0.25477342806656        0.29833449044618
H            0.25258994135455        2.29003272583950       -0.50924368461900
H           -2.19484023994761        1.89447265474598       -0.64881714668543
H           -3.07885997269744       -0.37376060890566       -0.20617568923744
H           -1.53013995234004       -2.21100772263288        0.37217230882333
H            1.90971200995189       -1.84610143383494        0.61828297835522
"""


@_needs_xtb
@_needs_orca
@pytest.mark.slow
def test_the_last_rung_is_orca_and_it_reaches_what_no_climb_here_does():
    """One of the sixteen whole drags is reached by ORCA's optimiser and nothing else.

    Salicylaldehyde with its hydroxyl proton dragged across its own hydrogen
    bond.  Both climbs settle from this structure onto something with no mode
    going the wrong way at all; ORCA, walking in redundant internal
    coordinates rather than in Cartesians, reaches the transfer at about
    -1174 cm-1 with 0.51 of its mode on the O0-H14 stretch.  Measured end to
    end: about 26 s, of which the two climbs are ten.

    Marked slow rather than trimmed.  It is the only test here that proves the
    third rung earns its place, and the whole of what it proves is that a
    minute of ORCA reaches something ten seconds of climbing does not -- there
    is no cheap version of that.
    """
    said = []
    got = climb.reach_the_reaction(_PROTON_HALF_MOVED, 'gfn2', held=(0, 14),
                                   aimed_from=_SALICYLALDEHYDE,
                                   on_route=said.append)
    assert got['route'] == 'ORCA', got.get('status')
    assert got['ok'], got.get('status')
    shape = got['imaginary']
    assert shape['count'] == 1 and shape['share'] > 0.4, shape
    assert shape['modes'][0] < -800, shape
    # Both retries were announced before they ran, not explained afterwards.
    assert len(said) == 2, said
    assert 'softest mode' in said[0] and 'ORCA' in said[1], said


def test_a_stop_does_not_pay_for_a_hessian_nobody_asked_for():
    """What a Stop is about is the frame on screen, not where the walk stood.

    The verdict is the expensive half of a climb -- xtb's own Hessian, 0.3 s
    at sixteen atoms and 11.8 at fifty -- and it used to be taken whatever
    ended the loop, including the press that ended it.  That is a third of a
    second and more spent at the exact moment somebody asked for the walk to
    stop, on a geometry they are not being left with: the editor keeps the
    frame the picture stopped on and throws the rest away.

    So a stopped climb says it was stopped and hands back no verdict at all,
    and the ladder reads that and stops too -- a rung is only tried when the
    one before it is known to have missed, and a rung that was stopped is not
    known to have missed anything.
    """
    got = climb.climb_to_saddle(_ESTIMATE, 'gfn2', should_stop=lambda: True)
    assert got['stopped'] is True
    assert got['imaginary'] is None
    assert 'stopped' in got['status']
    assert got['xyz'], 'where it stood is still handed back'


@_needs_xtb
def test_the_ladder_stops_where_the_press_was_stopped():
    """And says so, rather than reporting a minimum it never looked at."""
    got = climb.reach_the_reaction(_ESTIMATE, 'gfn2', held=(0, 10),
                                   should_stop=lambda: True)
    assert got['ok'] is False
    assert got['route'] == 'the climb', got.get('route')
    assert 'stopped' in got['status'], got['status']


def test_the_editor_hands_the_pair_over_and_never_guesses_it():
    """The ladder needs two atoms, and a wrong two are worse than none.

    Seven element-free rules were scored against twenty-one hand drags, each
    of them given the atom the page names and both geometries -- before the
    drag and after it.  The best is the contact the follow is actually
    holding, :func:`gfn_optimize.contacts_holding`, and it names the pair the
    gesture is about on **10 of the 21**; the other six score 1, 1, 1, 3, 7
    and 9.  They fail in three shapes: a drag at a carbon names the hydrogen
    on it, a drag that folds a chain names the torsion it turns rather than
    the contact it closes, and a proton walked between two oxygens names the
    one it started on.

    A wrong pair is not a smaller version of the right one -- it is the ladder
    throwing away an answer that was right.  Measured on the van-der-Waals
    complex with the pair guessed that way: the first climb reached the
    Diels-Alder saddle at 2.315 A and -394 cm-1, the check against C10-H4
    scored it at less than a fifth, and the press then spent 73 s on two more
    searches to report that it had not found what was already on the screen.

    So nothing is guessed.  The page names the atom under the cursor and the
    atoms that are picked, on messages it already sends and in the same
    numbering, and the pair is one of each: what you are dragging, and what
    you tapped to say what you are dragging it at.
    """
    source = EDITOR_SOURCE
    naming = source.split('def _name_the_pair_the_hand_means(holding):',
                          1)[1].split('\n    def ', 1)[0]
    assert "state.get('picked')" in naming
    assert "state['climb_held']" in naming
    # Both halves come off messages the page already sends, so nothing new is
    # recorded and a replayed journal names the same pair.
    assert '_gfn.contacts_holding' not in naming, \
        'the pair is named rather than guessed'
    climbing = source.split('def _climb_now():', 1)[1].split(
        '\n    def on_submit_', 1)[0]
    assert "held = state.pop('climb_held', None)" in climbing
    assert 'held=held,' in climbing
    # A pair belongs to the gesture that named it, so a new grab starts
    # without one rather than with the last drag's.
    grab = source.split("if verb == 'gfngrab':", 1)[1].split(
        "if verb == 'gfnfree':", 1)[0]
    assert "state.pop('climb_held', None)" in grab


@_needs_xtb
def test_the_pair_arrives_from_the_page_and_the_ladder_uses_it(tmp_path):
    """The gesture, through the page's own two messages and nothing else.

    An atom is tapped -- ``submit_pick_sync``, which is how the editor is told
    which atoms a user means everywhere else -- and another is dragged, which
    the page reports as ``DELFIN drag-follow held=10``.  Between them those
    two say which contact the gesture is about, and the release hands it to
    the climb.
    """
    import time as clock

    pytest.importorskip('ipywidgets')

    part, state, box = _an_editor(tmp_path)
    part.submit_relax_btn.value = True
    part.submit_pick_sync.value = '0'
    part.submit_climb_btn.value = True
    began = clock.time()
    while not state.get('climb_showing') and clock.time() - began < 60:
        clock.sleep(0.02)

    moved = _dragged(0.40).splitlines()
    part.submit_manip_sync.value = '\n'.join(
        [moved[0], 'DELFIN drag-follow held=10'] + moved[2:])
    assert state.get('climb_held') == (10, 0), state.get('climb_held')

    # And a gesture that names nothing leaves nothing behind, rather than the
    # last one's pair.
    part.submit_cmd_sync.value = 'gfngrab:991:0'
    assert 'climb_held' not in state
    part.submit_climb_btn.value = False
    _quiet(part, state)


def _pull_one_atom(part, box, grabbed, toward, reach, steps=6, pause=0.6):
    """One atom picked up and pulled at another, in the page's own messages.

    ``gfngrab``, a ``DELFIN drag-follow`` per mouse move naming the atom under
    the cursor, ``DELFIN drag-end`` and ``gfnfree`` -- in that order, because
    that is the order the page sends them in.  Each move is aimed straight at
    the other atom, which is what a hand does when it is dragging one thing at
    another.
    """
    import time as clock

    names = climb._elements(box.value)['symbols']
    start = _where(box.value)
    span = float(np.linalg.norm(start[toward] - start[grabbed]))
    _command(part, 'gfngrab', '0')
    for turn in range(1, steps + 1):
        rows = _where(box.value).copy()
        arm = rows[toward] - rows[grabbed]
        rows[grabbed] = rows[toward] - arm / float(np.linalg.norm(arm)) * (
            span - (span - reach) * turn / steps)
        part.submit_manip_sync.value = climb.xyz_document(
            names, rows, f'DELFIN drag-follow held={grabbed} n={turn}')
        clock.sleep(pause)
    part.submit_manip_sync.value = climb.xyz_document(
        names, _where(box.value), 'DELFIN drag-end')
    _command(part, 'gfnfree', '')


@_needs_xtb
@_needs_orca
@pytest.mark.slow
def test_the_gesture_reaches_a_saddle_that_one_climb_walks_away_from(tmp_path):
    """The whole of task one, as a gesture: tap, drag, let go.

    Dynamik Opt on so the molecule follows the hand, Climb to TS on so the
    release climbs, one ethene carbon dragged at a diene terminus until they
    stand 2.6 A apart, and the diene terminus tapped first so the editor knows
    which contact is meant.

    Measured on this drag, which is the point of the ladder rather than a
    detail of it: the climb aimed along the gesture converges -- in 61 steps
    and 21.8 s -- onto a structure with **no mode going the wrong way at
    all**, 2.693 A apart, and that is what the editor used to leave in the box
    after 211 steps and 83.8 s at its old ceiling of 400.  The same press with
    the three searches tried in turn ends on the Diels-Alder saddle: both
    forming bonds at 2.315 A, one imaginary mode at -393.5 cm-1, and 0.68 of
    that mode is the contact that was tapped.

    Slow, and there is no cheap version: the gesture is six xtb relaxations
    under the hand, then two climbs and an ORCA optimiser.
    """
    import time as clock

    pytest.importorskip('ipywidgets')

    part, state, box = _an_editor(tmp_path)
    part.submit_relax_btn.value = True          # Dynamik Opt
    part.submit_pick_sync.value = '0'           # the atom being aimed at
    part.submit_climb_btn.value = True          # Climb to TS
    began = clock.time()
    while not state.get('climb_showing') and clock.time() - began < 60:
        clock.sleep(0.05)

    _pull_one_atom(part, box, 10, 0, 2.6)
    assert state.get('climb_held') == (10, 0), state.get('climb_held')
    # What the hand left, and where the structure stood when it arrived, kept
    # so the same press can be run again with one search instead of three.
    released, was = box.value, state.get('climb_was')
    seen = _frames_seen(part)
    _quiet(part, state, seconds=600.0)

    # One press is one run, however many searches it takes.  The second climb
    # starts again from the structure the hand made, so the picture jumps back
    # to it and climbs away again -- and that must not be a new run: the
    # frames go into one list in the order they are drawn, and the count the
    # page reports back is an index into that list.
    runs = {run for run, _ in seen if run is not None}
    assert runs and max(runs) == state.get('climb_frame_run'), runs

    one, two = _forming(_where(box.value))
    assert 2.25 < one < 2.40 and abs(one - two) < 0.03, (one, two)
    said = part.mol_status.value
    assert 'reached it' in said, said
    assert 'C10-C0 stretch' in said, said

    # And the same gesture with only the first search, which is what this was
    # before: it converges, and onto a minimum.
    alone = climb.climb_to_saddle(released, 'gfn2', aimed_from=was,
                                  max_steps=climb.CLIMB_CEILING, held=(0, 10))
    assert alone['imaginary']['count'] == 0, alone['imaginary']
    assert not alone['imaginary']['reaction']
