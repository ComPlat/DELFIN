"""A transition state in seconds, on the machine the dashboard runs on.

xtb has no saddle-point optimiser -- 6.7.1 offers ``--hess``, ``--ohess``,
``--bhess`` and ``--path``, and nothing that climbs to a first-order saddle.
ORCA has one, and ORCA can be told to take its gradients from xtb.  That pair
is fast enough to be a button, which is what makes the interactive half of the
question answerable at all: pose a structure by hand or take what the path
finder estimated, press, and a few seconds later it is a converged saddle or
it is not, and which is said.

The other half stays a job.  A saddle search on a real basis set runs for
hours and belongs in the ORCA Builder, where OPTTS is now a job type.

Chained to the path finder it is a button that takes two structures and hands
back a converged saddle: twelve seconds on sixteen atoms, and the same saddle
a nudged elastic band reaches -- to 0.07 cm-1, for about twice the gradients.
(Seven minutes was written here for the band, and that was the serial number:
measured since, the same band is 272 s on one process and 39.4 s on eight,
because ORCA computes the images at once.  The band is the second entry in the
box beside this press because it costs more work, not because it takes longer
on a machine with cores.)  And it says what was reached,
which is not optional -- a saddle search does not fail when it goes wrong, it
succeeds at arriving somewhere, and a structure that is a maximum in two
directions at once is named here rather than reported as done.
"""

import pytest

from delfin.dashboard import saddle
from editor_source import EDITOR_SOURCE

_needs_orca = pytest.mark.skipif(saddle.find_orca() is None,
                                 reason="ORCA not installed")


def _has_xtb():
    from delfin.dashboard import gfn_optimize as gfn
    return gfn.find_binary("gfn2") is not None


_needs_xtb = pytest.mark.skipif(not _has_xtb(), reason="xtb not installed")

#: The transition state the path finder estimated for the Diels-Alder: the two
#: forming bonds at 2.524 and 2.520 A, one imaginary mode at -131.4 cm-1.  Not
#: symmetric, and not converged -- which is what makes it worth optimising.
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


def test_a_saddle_search_here_is_xtb_through_orca():
    """Because neither of them can do it alone.

    xtb has no saddle optimiser at all, and ORCA on a real basis set is a job
    rather than a press.  ORCA driving xtb is the pair that is both: measured
    on this estimate, sixteen atoms, under seven seconds.

    Three of the four are ORCA's own xtb, asked for by keyword.  g-xTB is not
    and cannot be -- it is a build of its own, and no ORCA keyword names it --
    so it is driven through ExtOpt, the interface ORCA publishes for programs
    it does not know; see :mod:`test_the_saddle_reaches_g_xtb`.
    """
    assert saddle.SADDLE_METHODS == {'gfn2': 'XTB2', 'gfn1': 'XTB1',
                                     'gfnff': 'XTBFF', 'gxtb': 'ExtOpt'}
    # Anything with a basis set is refused here and named as a job.
    refused = saddle.optimise_to_saddle(_ESTIMATE, 'PBE0')
    assert not refused['ok']
    assert 'runs on xtb through ORCA' in refused['status']
    # And nothing to work on is said rather than run.
    assert not saddle.optimise_to_saddle('', 'gfn2')['ok']


def test_what_it_reached_is_read_from_the_last_hessian():
    """A saddle search recalculates its Hessian as it goes, so the output holds
    several: the last is nearer the geometry it ended on than the earlier ones,
    which describe geometries it has left.

    Read the first -- which is what a naive parse does -- an OptTS that
    improved its structure reports the mode it started with, and the number on
    screen is about a structure that is no longer anywhere.

    The real modes come back too, and separately.  How many modes go the wrong
    way says what order of saddle a structure is; how soft the ones that go
    the right way are says how nearly it is one of a higher order, and that is
    a question published searches ask as well.
    """
    output = (
        ' projected vibrational frequencies (cm-1)\n'
        'eigval :      -0.00    -0.00    -0.00     0.00     0.00     0.00\n'
        'eigval :    -131.40    28.31   141.02\n'
        ':  # imaginary freq.    1  :\n'
        ' projected vibrational frequencies (cm-1)\n'
        'eigval :      -0.00    -0.00    -0.00     0.00     0.00     0.00\n'
        'eigval :    -372.26    86.79   167.99\n'
        ':  # imaginary freq.    1  :\n'
    )
    modes = saddle._last_modes(output)
    assert modes == {'count': 1, 'modes': [-372.26],
                     'real': [86.79, 167.99]}, modes
    # The six xtb projects out are printed as -0.00 and 0.00 exactly, and are
    # not vibrations; a genuinely soft mode of a few wavenumbers is, and is
    # the whole point of reading them.
    soft = saddle._last_modes(
        ' projected vibrational frequencies (cm-1)\n'
        'eigval :      -0.00    -0.00    -0.00     0.00     0.00     0.00\n'
        'eigval :     -18.20     3.41     6.02   140.11\n'
        ':  # imaginary freq.    1  :\n')
    assert soft['real'] == [3.41, 6.02, 140.11], soft
    # Nothing to read is None, not an empty answer that reads as "a minimum".
    assert saddle._last_modes('nothing here') is None


@_needs_orca
def test_the_estimate_becomes_a_converged_transition_state():
    """Measured: the estimate comes in with its two forming bonds at 2.524 and
    2.520 A and one imaginary mode at -131.4 cm-1, and comes out symmetric at
    2.315 and 2.315 with the mode at -393.5.

    Seven seconds, which is why this is a button.

    The -372 this said before came from the optimiser's own last Hessian
    block, which is about a geometry a few steps back: two runs of this same
    input reported -372 and -301.55 for the same converged structure, and a
    Hessian on the structure that came back says -393.5 both times.
    """
    import math

    from delfin.dashboard import gfn_optimize as gfn

    found = saddle.optimise_to_saddle(_ESTIMATE, 'gfn2', timeout=300)
    assert found['ok'], found['status']
    assert found['seconds'] < 120, found['seconds']

    here = gfn.coordinates_of(found['xyz'])
    one = math.dist(here[0:3], here[30:33])
    two = math.dist(here[9:12], here[33:36])
    assert 2.1 < one < 2.6, one
    # It converged to something symmetric, which the estimate was not.
    assert abs(one - two) < 0.02, (one, two)

    # And it is a transition state, which is the whole reason to run it.
    assert found['imaginary']['count'] == 1, found['imaginary']
    assert found['imaginary']['modes'][0] < -100.0, found['imaginary']


def test_it_runs_where_the_dashboard_runs_and_is_kept_small():
    """On a cluster that is the login node, which is a place with rules.

    Eight cores, two gigabytes and three minutes -- welcome exactly as long
    as it is small and short.  Eight because it was measured: on a
    fifty-one-atom molecule, 1, 4, 8, 16 and 32 cores took 34.5, 12.1, 7.9,
    6.3 and 5.9 s, so four leaves most of the gain on the table and sixteen
    buys a fifth of a speedup for four times the node.  On sixteen atoms it
    makes no difference at all -- 6.7 s on one core and 6.7 on sixteen --
    because the system is over before any of it can be shared out.

    And never more than the machine has, which matters on a small login node
    rather than on a big shared one.
    """
    import inspect

    source = inspect.getsource(saddle)
    assert 'cores: int = 8,' in source
    assert 'timeout: Optional[float] = 180.0' in source
    assert "'%maxcore 2000\\n'" in source
    # Eight processes where ORCA is doing the arithmetic; one where it is not,
    # because through ExtOpt the gradient belongs to another program and
    # ORCA's own parallel numerical Hessian came apart twice on it.
    assert 'ranks = 1 if own_program is not None else _share(cores)' in source
    assert 'nprocs {ranks}' in source
    assert 'min(want, os.cpu_count() or 1)' in source
    assert 'submit it as a job for more' in source
    # A saddle search needs a Hessian to know which mode to climb, and
    # recalculates it because the mode changes character as it moves.  The
    # cadence is a named constant now rather than a literal in the input,
    # because under g-xTB one Hessian is forty-odd separate processes and
    # dropping the recalculation is a tempting saving -- measured, it walks
    # the Diels-Alder estimate back down to the van-der-Waals complex.
    assert "'  Calc_Hess true\\n'" in source
    assert 'RECALC_HESS = 5' in source
    assert "f'  Recalc_Hess {RECALC_HESS}\\n'" in source


def test_the_button_works_on_whatever_is_in_the_box():
    """Not only on what the path finder left.

    A structure posed by hand into something that looks like a transition
    state is exactly the case this is for, and it is the interactive half of
    the question.
    """
    source = EDITOR_SOURCE
    assert 'submit_saddle_btn = widgets.Button(' in source
    assert "description='To the saddle'" in source
    assert 'def on_submit_saddle(_button=None):' in source
    assert 'submit_saddle_btn.on_click(on_submit_saddle)' in source
    # Whatever is on screen, and not a remembered pair of ends.
    assert 'xyz = _current_xyz()' in source
    # One step, which Undo takes back whole.
    assert "_remember('the saddle search')" in source
    # And the box is named by what was reached rather than by what was being
    # looked for: a comment that says "transition state" over a second-order
    # saddle is the one place the mistake would outlive the sentence.
    assert 'f\'Optimised to {said["name"]}\'' in source
    # What it reached is said, whichever it is -- in the words the module
    # that owns the thresholds says it in.
    assert 'def _said_modes(shape, what, advise=True):' in source
    assert "_saddle.verdict(shape, what, advise=advise)['lines']" in source
    assert "said = _saddle.verdict(found.get('imaginary')," in source
    # And where it goes next, when there is a next.
    assert "if said['first_order']:" in source
    assert 'Refine it with OPTTS in the ORCA Builder at ' in source


def test_only_whole_frames_are_read_from_a_trajectory_being_written():
    """ORCA appends to the trajectory while it climbs, so a read lands in the
    middle of a block as often as not.

    Half a geometry drawn is worse than a geometry drawn a tenth of a second
    later, so a partial block is left for the next look.
    """
    whole = ('2\nCoordinates from ORCA-job in E -17.813830573670\n'
             'C 0.0 0.0 0.0\nC 1.0 0.0 0.0\n')
    walked, spent = saddle.frames_in(whole * 2, 2)
    assert len(walked) == 2
    assert walked[0] == [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    # The energy rides on the comment line, which is where ORCA puts it.
    assert spent == [-17.813830573670, -17.813830573670]

    # A block cut off mid-write is not a frame yet.
    cut = whole + '2\nCoordinates from ORCA-job in E -17.8\nC 0.0 0.0 0.0\n'
    assert len(saddle.frames_in(cut, 2)[0]) == 1
    # And a file with nothing in it is no frames rather than an error.
    assert saddle.frames_in('', 2) == ([], [])


@_needs_orca
def test_the_climb_is_watched_while_it_climbs_and_can_be_stopped():
    """Measured on the sixteen-atom estimate: thirteen steps arriving between
    2.8 s and 6.8 s, about one every quarter second, each with its energy.

    Watching is worth more here than on a minimisation.  A minimisation that
    follows the wrong thing fails and says so; a climb does not -- it succeeds
    at arriving somewhere nobody wanted, and the only way to know that early
    is to see it happen.  So the same press stops it, and what it had reached
    is kept rather than thrown away.
    """
    arrived = []

    found = saddle.optimise_to_saddle(
        _ESTIMATE, 'gfn2', timeout=300,
        on_frame=lambda walked, spent: arrived.append((len(walked), spent[-1])))
    assert found['ok'], found['status']
    # The frames came in while it ran, not in one lump at the end.
    assert len(arrived) > 3, arrived
    assert [n for n, _ in arrived] == sorted(n for n, _ in arrived)
    assert len(found['path']) == arrived[-1][0]
    assert len(found['path'][0]) == 16 * 3
    # Every step carries the energy ORCA wrote beside it, so the climb can be
    # said as well as shown.
    assert all(one is not None for one in found['energies'])
    assert found['energies'][-1] > found['energies'][0], 'a climb goes up'

    # And the same press stops it, keeping what it had.
    got = {'n': 0}
    halted = saddle.optimise_to_saddle(
        _ESTIMATE, 'gfn2', timeout=300,
        on_frame=lambda walked, spent: got.update(n=len(walked)),
        should_stop=lambda: got['n'] >= 3)
    assert not halted['ok']
    assert halted['halted'] is True
    assert 'Stopped' in halted['status']
    assert halted['seconds'] < found['seconds'], 'it really did end early'
    # The structure it had reached is kept -- a climb that is stopped has
    # still been somewhere, and throwing that away wastes the whole run.
    assert len(_gfn_atoms(halted['xyz'])) == 16


def _gfn_atoms(text):
    from delfin.dashboard import gfn_optimize as gfn
    return gfn.atom_lines(text or '')


def test_the_button_shows_the_climb_and_the_same_press_stops_it():
    """The frame channel, not the coordinate box.

    Every write to the box rebuilds the viewer from nothing, so a climb of
    thirty steps would be thirty rebuilds -- which is what made the scan crawl
    before it was moved to the same channel.
    """
    source = EDITOR_SOURCE
    assert "state['saddle_frame_run'] = int(state.get('gfn_run', 0)) + 1" in source
    assert "state['gfn_run'] = state['saddle_frame_run']" in source
    assert 'on_frame=_watch,' in source
    assert "should_stop=lambda: bool(state.get('saddle_stop'))" in source
    # A run the page has been told to forget does not draw over the next one.
    assert "if state.get('gfn_run') != state.get('saddle_frame_run'):" in source
    # The same button, because there is only one thing to want while it runs.
    assert "submit_saddle_btn.description = 'Stop'" in source
    # Named where the press is named rather than in the run that stops:
    # what it goes back to is whatever the two boxes say at the time.
    assert "_name_the_saddle_press()" in source
    # Not disabled while it runs -- a button nobody can press cannot stop it.
    assert 'submit_saddle_btn.disabled = True' not in source


#: The two ends of the Diels-Alder, each relaxed under GFN2: s-cis butadiene
#: with an ethylene above it, and the cyclohexene the two of them make.  The
#: pair a scan leaves behind, and the input the chain is measured on.
_REACTANT = """16
butadiene and ethylene, relaxed
C           -1.51175969445403       -0.02647380860866       -0.06557436490527
C           -0.72674009436823        1.04324898539530       -0.13156348350727
C            0.72671613615743        1.04323688164824       -0.13200313736446
C            1.51174712621305       -0.02646817123412       -0.06592265675999
H           -2.58463144839309        0.06203170569435       -0.06842720691043
H           -1.11830994961107       -1.02635341307553       -0.00678112167498
H           -1.18275196117662        2.02412528656961       -0.19063949165969
H            1.18273106412604        2.02410195941484       -0.19113584108773
H            2.58463111371453        0.06202058682634       -0.06829703412017
H            1.11831903097964       -1.02629431053540       -0.00619815580858
C           -0.65812656382745       -0.39453753897557        3.15439731529169
C            0.65814270849405       -0.39445452845154        3.15457854734613
H           -1.22868654033087        0.51904761763282        3.14158888465065
H           -1.22983917433790       -1.30795394609327        3.16689639740950
H            1.22859345191931        0.51920014666317        3.14189155138364
H            1.22996179489533       -1.30780145287077        3.16718979771703
"""

_PRODUCT = """16
cyclohexene, relaxed
C           -1.36290194815897       -0.06313150970649        0.75107481477849
C           -0.66218260952155        1.09869132587048        0.12233880799005
C            0.66173055731454        1.09913795139967        0.12248552655879
C            1.36310557821274       -0.06226438766873        0.75126945964802
H           -2.43534601907205        0.11530592358156        0.83579350839840
H           -1.21914017994592       -0.94894767327948        0.12363419262958
H           -1.24659833279792        1.90123497817873       -0.30035668134291
H            1.24569990280881        1.90207650485640       -0.30007507370404
H            2.43541532174521        0.11689111711109        0.83621469395190
H            1.22007694226306       -0.94810437965944        0.12369377517092
C           -0.76948705831100       -0.34196841669303        2.14125370859433
C            0.76954534349931       -0.34170006236609        2.14127653272432
H           -1.13126031423455        0.42887199800881        2.82393080033168
H           -1.14957836569781       -1.29937167674247        2.50173435160042
H            1.13098517132138        0.42903787478979        2.82424151569521
H            1.14993301057964       -1.29908356768101        2.50149006697273
"""

#: What the chain converged onto when it was given a sloppy pair of ends: two
#: ethylenes face to face, held at arm's length, and the four-membered ring
#: they close into.  A stationary point -- ORCA converges on it and reports
#: one -- and not a transition state: two modes go the wrong way, at -54 and
#: -27 cm-1, because it is a maximum in two directions at once.
_SECOND_ORDER = """12
what the chain converged onto from a sloppy pair of ends
C            1.74870183904621       -0.80313936244350        0.01234407130288
C           -0.80897069838947        1.72403695528328       -0.01245004685042
C           -1.74814210715612        0.80256721072369        0.01475980854592
C            0.80840915289181       -1.72346647129475       -0.01466167881697
H            2.13488353013285       -0.40853974310695        0.93769394663060
H            2.17735181162523       -0.39821483871238       -0.88954509846445
H           -0.38622435805980        2.13518147547531        0.88940341366899
H           -0.41575507758763        2.11123600782113       -0.93794119413425
H           -2.14029055727280        0.41427600233891        0.94026074267805
H           -2.16984580770466        0.39034432396215       -0.88710971152635
H            0.37877639229415       -2.12729908157556        0.88723440214583
H            0.42110565797461       -2.11698286947635       -0.93998865504220
"""


def test_the_thresholds_belong_to_other_people_and_are_cited():
    """Not invented here, and not adjustable to taste.

    Whether a structure is a transition state is a question other searches
    have had to answer at scale, and the numbers they settled on are worth
    more than a number chosen to make this case come out well.  Each is cited
    where it is defined, and each was checked against its own source rather
    than against a paper about it:

    autodE keeps ``min_imag_freq = Frequency(-40, units="cm-1")`` in
    ``autode/config.py`` and refuses a structure whose most negative mode is
    above it, with the comment that "although most TSs have |v_imag| > 100
    cm-1 this threshold is designed to be conservative".

    AutoMeKin's ``imagmin`` is documented as "the minimum value for the
    imaginary frequency (in absolute value and cm-1) of the selected TS
    structures".  Its code default is 0 -- 200 is what its tutorial uses and
    what every input file it ships with sets, so it is said here as a
    recommendation and not as a default.  With ``tight_ts``, which it has on
    unless told otherwise, it also refuses a structure whose lowest real
    frequency is negative and one whose two lowest real frequencies sum to
    less than ten.

    YARP settles it in a line: ``is_valid_ts`` counts the modes below a
    numerical tolerance and returns ``n_imaginary == 1``.
    """
    import inspect

    assert saddle.MIN_IMAGINARY == -40.0
    assert saddle.CLEAR_IMAGINARY == 200.0
    assert saddle.FLAT_PAIR == 10.0
    assert saddle.FIRST_ORDER == 1
    source = inspect.getsource(saddle)
    for cited in ('autodE', 'AutoMeKin', 'YARP', 'min_imag_freq', 'imagmin',
                  'is_valid_ts', 'tight_ts'):
        assert cited in source, cited


def test_two_modes_going_the_wrong_way_are_named_not_reported_as_done():
    """A second-order saddle is not a worse transition state; it is not one.

    A first-order saddle is a point a path can pass through.  A structure
    that is a maximum in two directions at once is not a point anything
    passes through, and calling it a transition state with a caveat would
    leave the caveat to be read.  So it is named.

    The numbers are what a Hessian found on the structure a converged chain
    handed back, and nothing in that run had said anything was wrong.
    """
    said = saddle.verdict({'count': 2, 'modes': [-53.92, -27.21],
                           'real': [16.91, 34.41]})
    assert said['order'] == 2
    assert said['first_order'] is False
    assert said['name'] == 'a second-order saddle'
    assert 'second-order saddle' in said['lines'][0]
    assert 'not a transition state' in said['lines'][0]
    assert '-54, -27 cm-1' in said['lines'][0]
    # And what to do about it, since converging again from here finds this.
    assert 'along the second of those modes' in said['lines'][1]

    # Three is not two, and is not a special case either.
    assert saddle.verdict({'count': 3, 'modes': [-90.0, -40.0, -20.0]})[
        'name'] == 'a third-order saddle'
    assert saddle.verdict({'count': 7, 'modes': [-90.0]})['name'] == \
        'a saddle point of order 7'

    # The wording never assumes what is being computed.  Every kind of system
    # is, so a sentence that spoke of bonds forming or a reaction crossing
    # would be wrong about most of them.
    for shape in ({'count': 0, 'modes': []},
                  {'count': 1, 'modes': [-393.53], 'real': [119.3, 173.4]},
                  {'count': 2, 'modes': [-53.92, -27.21]}):
        for line in saddle.verdict(shape)['lines']:
            for assumed in ('bond', 'reaction', 'reactant', 'product',
                            'barrier', 'atom'):
                assert assumed not in line.lower(), (assumed, line)


def test_a_mode_too_shallow_to_believe_is_not_taken_for_a_saddle():
    """One imaginary mode is what a transition state has; a very small one is
    as easily the method's noise as the structure's shape.

    autodE draws the line at -40 cm-1 and says in its own comment that this is
    deliberately conservative, most transition states being past -100.
    Conservative is what is wanted here: this is a semiempirical answer meant
    to be handed on to a real one, so refusing something that was a saddle
    costs a press and accepting something that was not costs a job.
    """
    shallow = saddle.verdict({'count': 1, 'modes': [-22.0], 'real': [60.0]})
    assert shallow['first_order'] is False, shallow
    assert 'not one on this evidence' in shallow['lines'][0]
    assert '-40 cm-1' in shallow['lines'][0]
    # And not named one either.  The name is what goes into the coordinate
    # box, and a box that says "transition state" over a structure the
    # sentence has just refused is where the refusal would be lost.
    assert shallow['name'] == \
        'a stationary point with one shallow mode the wrong way'

    deep = saddle.verdict({'count': 1, 'modes': [-393.53],
                           'real': [119.3, 173.44]})
    assert deep['first_order'] is True
    assert deep['name'] == 'a transition state'
    assert deep['lines'] == [
        'What it reached is a transition state: one mode goes the wrong way, '
        'at -394 cm-1, and no others.'], deep['lines']

    # A soft one that is past the line is accepted and said to be soft: the
    # inputs AutoMeKin ships ask for 200 cm-1 before they will take a saddle.
    soft = saddle.verdict({'count': 1, 'modes': [-75.0], 'real': [40.0, 60.0]})
    assert soft['first_order'] is True
    assert '200 cm-1' in soft['lines'][1]

    # A structure flat in two of its remaining directions is close to being a
    # saddle of higher order, which is AutoMeKin's other refusal.
    flat = saddle.verdict({'count': 1, 'modes': [-393.53],
                           'real': [3.4, 6.0, 140.0]})
    assert flat['first_order'] is True
    assert 'summing to under the 10 cm-1' in flat['lines'][1]

    # And a mode on the wrong side of zero that the engine did not count is
    # said rather than dropped -- xtb counts against a cutoff of -20 cm-1, so
    # a mode at -11 is in the list and not in the count.
    under = saddle.verdict({'count': 1, 'modes': [-393.53, -11.25],
                            'real': [119.3]})
    assert 'still on the wrong side of zero (-11 cm-1)' in under['lines'][1]

    # Nothing to judge is said as nothing to judge, not as a minimum.
    assert saddle.verdict(None)['order'] is None
    assert 'could not be checked' in saddle.verdict(None)['lines'][0]


@_needs_orca
@_needs_xtb
def test_two_structures_become_a_converged_saddle_at_one_press():
    """The whole chain, measured end to end on this machine.

    Two ends of a Diels-Alder in, sixteen atoms, both relaxed under GFN2:
    xtb's path finder walks between them in 3.6 s and hands back an estimate
    with its forming bonds at 2.570 and 2.554 A and one imaginary mode at -79
    cm-1, and ORCA's OptTS on xtb gradients sharpens that in 8 s into a
    symmetric structure at 2.315 and 2.315 with one mode at -393.5.

    Twelve seconds for the pair.  The routes that do the same job inside ORCA
    alone were measured on the same system: ``ScanTS`` 98 s and -394.1,
    ``NEB-TS`` 416 s and -393.6.  All of them land on the same saddle to
    within a wavenumber, which is the whole argument for chaining the two
    engines rather than paying either of them to do both halves.

    The band's 416 s was a one-process wall time, and re-measured for
    ``neb_to_saddle`` the same band is 272 s on one process and 39.4 s on
    eight -- ORCA computes the images at once.  So the comparison that stands
    is the work: 203 gradients for the band, and rather fewer for this chain.
    The chain is offered first because it is cheaper, not because a band is
    slow on a machine with cores.
    """
    import math

    from delfin.dashboard import gfn_optimize as gfn

    seen = []
    stages = []
    found = saddle.path_to_saddle(
        _REACTANT, _PRODUCT, 'gfn2', timeout=300,
        on_frame=lambda walked, spent: seen.append(len(walked)),
        on_stage=stages.append)
    assert found['ok'], found['status']
    assert found['stage'] == 'saddle'
    assert found['seconds'] < 120, found['seconds']
    # Both halves' times, because which half a press was spent in is the
    # question anyone asks of a chain that felt slow.
    assert found['path_seconds'] > 0
    assert found['saddle_seconds'] > 0
    assert found['seconds'] >= found['path_seconds']

    # The walk's own answer is kept beside the climb's: a barrier forward and
    # back, and how near it came to the structure it was aimed at.
    assert 1.0 < found['barrier'] < 12.0, found['barrier']
    assert found['back'] > 40.0, found['back']
    assert found['rmsd'] is not None and found['rmsd'] < 0.3, found['rmsd']
    assert len(gfn.atom_lines(found['estimate'])) == 16

    # And the climb was watched while it climbed, one frame at a time, so the
    # viewer has something to draw before the answer arrives.
    assert len(seen) > 3, seen
    assert seen == sorted(seen)
    assert stages and 'sharpening' in stages[0]

    here = gfn.coordinates_of(found['xyz'])
    one = math.dist(here[0:3], here[30:33])
    two = math.dist(here[9:12], here[33:36])
    assert 2.1 < one < 2.6, one
    assert abs(one - two) < 0.02, (one, two)

    # And it is a transition state, said in those words.
    assert found['imaginary']['count'] == 1, found['imaginary']
    assert found['imaginary']['modes'][0] < -300.0, found['imaginary']
    assert found['verdict']['first_order'] is True
    assert found['verdict']['name'] == 'a transition state'


@_needs_orca
@_needs_xtb
def test_the_frequency_is_about_the_structure_that_came_back():
    """An OptTS recalculates its Hessian every few steps and converges
    somewhere between two of them, so the last block in its output belongs to
    a geometry it has since left.

    Measured on the sixteen-atom estimate: the run recalculated at cycles 1, 6
    and 11 and converged at 16, and its last block says -301.55 cm-1 while a
    Hessian on the structure it handed back says -393.53.  The same structure,
    two answers, and only the second is about anything on screen.  0.3 s to
    settle it, against a number that would otherwise be quoted for a geometry
    that exists nowhere.
    """
    told = saddle.optimise_to_saddle(_ESTIMATE, 'gfn2', timeout=300)
    assert told['ok'], told['status']
    assert told['confirmed'] is True
    assert told['imaginary']['modes'][0] < -380.0, told['imaginary']
    # The real modes come with it, which is what the flatness check needs.
    assert told['imaginary']['real'], told['imaginary']

    # Without the confirmation it is the optimiser's own last block, which is
    # about a geometry a few steps back -- kept as the fallback for when the
    # Hessian cannot be taken, and measurably not the same number.
    raw = saddle.optimise_to_saddle(_ESTIMATE, 'gfn2', timeout=300,
                                    confirm=False)
    assert raw['ok'], raw['status']
    assert raw['confirmed'] is False
    assert raw['imaginary']['modes'][0] > told['imaginary']['modes'][0]


@_needs_orca
@_needs_xtb
def test_a_converged_saddle_search_can_have_reached_something_else():
    """This is the reason the chain cannot be trusted blind.

    A minimisation that follows the wrong thing fails and says so.  A saddle
    search does not fail; it succeeds at arriving somewhere.  Measured here:
    given two ethylenes held face to face and the ring they close into -- a
    deliberately sloppy pair of ends -- the chain converged, reported a
    stationary point, and had arrived at a structure that is a maximum in two
    directions at once.  Nothing in the run said so.

    Restarting the search from that structure reaches it again in 4 s and six
    steps, which is what this test does: ORCA converges, reports a stationary
    point on XTB2, and the check names what was reached as a second-order
    saddle instead of reporting it as a transition state.  Two modes go the
    wrong way, at -53.92 and -27.21 cm-1, with a third at -14.7 under the
    cutoff xtb counts against.
    """
    found = saddle.optimise_to_saddle(_SECOND_ORDER, 'gfn2', timeout=300)
    # The run itself is a success, which is the point: nothing about it is
    # broken, and it has not found a transition state.
    assert found['ok'], found['status']
    assert 'reached a stationary point' in found['status']
    assert found['imaginary']['count'] == 2, found['imaginary']
    assert found['imaginary']['modes'][0] < -50.0, found['imaginary']
    assert found['imaginary']['modes'][1] < -20.0, found['imaginary']

    said = found['verdict']
    assert said['order'] == 2
    assert said['first_order'] is False
    assert said['name'] == 'a second-order saddle'
    assert 'is a second-order saddle, not a transition state' in said['lines'][0]
    # And the box would be labelled with that name rather than with the one
    # that was being looked for.
    assert 'transition state' not in f'Optimised to {said["name"]}'


def test_the_chain_is_one_press_and_the_same_press_stops_it():
    """One press for both halves, because a chain the user has to shepherd is
    two presses with a wait in the middle.

    A press of its own rather than something Find the path does on its way
    past.  The path's own answer -- a barrier forward, a barrier back, and how
    near it came to what it was aimed at -- is complete in a quarter of the
    time and needs no ORCA at all, and the climb can end somewhere nobody
    asked for.  A press that answers one question should not spend four times
    as long answering a second one and replace the first answer with a
    structure nobody has looked at yet.

    The climb streams into the viewer the way the saddle button's does: down
    the frame channel with a run id of its own, never into the coordinate box,
    because every write to the box rebuilds the viewer from nothing.
    """
    source = EDITOR_SOURCE
    assert 'def _path_then_orca(ends):' in source
    assert "('through ORCA', 'orca')" in source
    assert '_path_then_orca(ends)' in source
    assert 'submit_saddle_btn.on_click(on_submit_saddle)' in source
    assert '_saddle.path_to_saddle(' in source
    # Its own run id, so a run the page has been told to forget cannot draw
    # over this one.
    assert "state['chain_frame_run'] = int(state.get('gfn_run', 0)) + 1" in source
    assert "state['gfn_run'] = state['chain_frame_run']" in source
    assert "if state.get('gfn_run') != state.get('chain_frame_run'):" in source
    assert "'frames': [[round(float(v), 4)" in source
    # The same press stops it, and is not disabled while it runs -- a button
    # nobody can press cannot stop anything.
    assert "submit_saddle_btn.description = 'Stop'" in source
    assert "state['chain_stop'] = True" in source
    assert "should_stop=lambda: bool(state.get('chain_stop'))" in source
    assert 'submit_saddle_btn.disabled = True' not in source
    # One step for the whole chain, which Undo takes back whole.
    assert "_remember('the path and the climb')" in source
    # And two climbs never draw into one viewer.
    assert "if state.get('path_run') or state.get('chain_run'):" in source
    assert "if state.get('path_run') or state.get('saddle_run'):" in source


def _an_editor(text):
    """One structure editor over a coordinate box of its own.

    The real part, driven the way the button drives it: this class of defect
    survives every reading of the source, because the source says exactly what
    it means to do.
    """
    import pathlib
    import tempfile

    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda _script: None
    state = {}
    box = widgets.Textarea(value=text)
    part = structure_editor.build(
        ctx, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    return part, state, box


@_needs_orca
@_needs_xtb
def test_one_press_of_the_real_button_walks_climbs_and_draws_it_once():
    """The button, not the function underneath it.

    Two structures marked one at a time, one press, and 16 s later the box
    holds a converged saddle named by what it is.  Measured here: the walk 5.2
    s and the climb 10.5, sixteen frames drawn.

    Every frame arrives once and in order, down the frame channel.  That is
    the requirement the viewer lives by -- a write to the coordinate box
    rebuilds the whole viewer from nothing, so a sixteen-step climb written to
    the box would be sixteen rebuilds and the same trajectory shown over and
    over.  So the box is written once, at the end, with the answer.
    """
    import json
    import time

    from delfin.dashboard import gfn_optimize as gfn

    part, state, box = _an_editor(_PRODUCT)
    part.submit_ff_dd.value = 'gfn2'

    # Marked one at a time, which is how a pair of ends is given without
    # anything having to hold two structures at once.
    state['current_xyz_for_copy'] = {'content': _REACTANT}
    part.on_submit_path_from()
    state['current_xyz_for_copy'] = {'content': _PRODUCT}
    part.on_submit_path_from()
    # The pair is now a start the press can be given, and every way of
    # walking between the two ends is on the list beside it.
    part._refresh_saddle_controls()
    # One start is not a choice, so the box stays away and the press is what
    # appears -- which is the whole of what the user needs to see.
    assert part.submit_saddle_btn.layout.display == ''
    assert 'marked' in [value for _label, value
                        in part.submit_saddle_from.options]
    part.submit_saddle_from.value = 'marked'
    # A nudged elastic band arrived here later as a fourth way from the same
    # pair.  Second on the list, because the order is the recommendation: it
    # reaches the same saddle as the chain and spends about twice the
    # gradients doing it.
    assert [value for _label, value in part.submit_saddle_how.options] == \
        ['orca', 'neb', 'hand', 'walk']

    frames = []
    part.submit_gfn_frame.observe(
        lambda change: frames.append(json.loads(change['new'])), names='value')
    written = []
    box.observe(lambda change: written.append(change['new']), names='value')

    began = time.time()
    part.on_submit_saddle()
    # The same button stops it, so it says so while it runs.
    assert part.submit_saddle_btn.description == 'Stop'
    while state.get('chain_run') and time.time() - began < 600:
        time.sleep(0.05)
    assert not state.get('chain_run'), 'the chain never finished'
    assert part.submit_saddle_btn.description == 'To the saddle'

    # Every frame once, in order, and all of them this run's.
    assert len(frames) > 3, frames
    assert [one['from'] for one in frames] == list(range(len(frames)))
    assert {one['run'] for one in frames} == {state['chain_frame_run']}
    assert {len(one['frames'][0]) // 3 for one in frames} == {16}

    # And the box written once, at the end, with what was reached named in it.
    assert len(written) == 1, len(written)
    assert written[0].splitlines()[1] == \
        'From a path, optimised to a transition state'
    assert len(gfn.atom_lines(written[0])) == 16
    # Which Undo takes back whole: two structures went in and one came out,
    # and halfway back is not a place anyone asked for.
    history = state.get('history') or []
    assert len(history) == 1, history
    assert history[-1].get('what') == 'the path and the climb', history[-1]


def test_a_run_that_did_not_converge_is_still_told_what_it_reached():
    """Not converged is not the same as nothing reached.

    Reported from a real session on a C-Br bond breaking: OptTS pressed again
    and again, each time 'did not converge', and nothing said whether what it
    had reached was a saddle at all.  A dissociation saddle is flat, and ORCA
    can walk close to one and stop a hair short of its gradient thresholds in
    the cycles it has -- so the geometry it stopped at is often a first-order
    saddle whatever the gradient said.

    So the branch that used to say only 'did not converge' now takes the same
    Hessian the converged branch takes, on the geometry it stopped at, and
    branches on what a Hessian says it is: one imaginary mode of its own is a
    transition state to tighten as a job, and zero or several is named as the
    minimum or higher saddle it actually reached.
    """
    import inspect

    body = inspect.getsource(saddle.optimise_to_saddle)
    # The not-converged branch is the one guarded by the DONE marker being
    # absent; everything asserted here has to sit after that guard.
    after = body.split('if not _DONE_RE.search(output):', 1)[1]
    after = after.split('# What was reached, from a Hessian', 1)[0]
    # It takes a Hessian on what it reached ...
    assert 'optimize_with_gfn' in after, (
        'the not-converged branch reports without checking what it reached')
    assert 'free_energy=True' in after, 'the check has to be a Hessian'
    # ... reads the verdict from it ...
    assert 'verdict(' in after
    # ... and lets a first-order saddle come back as a reached transition
    # state rather than a bare failure.
    assert "read.get('first_order')" in after
    assert 'ok=True' in after and 'converged=False' in after
    # And where it is not a saddle, it names what it is.
    assert "read[\"name\"]" in after or "read['name']" in after


@_needs_orca
def test_a_flat_climb_that_stops_short_names_the_minimum_it_found():
    """Live, on a start OptTS cannot converge in the one cycle it is given: a
    pyramidal ammonia, whose planar saddle is several umbrella-mode cycles
    away.  The old branch said 'did not converge' and stopped; now a Hessian
    on the geometry it reached says which stationary point it is, so the user
    learns the climb went down to a minimum rather than up to the saddle."""
    pyramidal = ("4\nammonia pyramidal\n"
                 "N   0.0000   0.0000   0.35\n"
                 "H   0.9377   0.0000  -0.10\n"
                 "H  -0.4688   0.8121  -0.10\n"
                 "H  -0.4688  -0.8121  -0.10\n")
    got = saddle.optimise_to_saddle(pyramidal, 'gfn2', charge=0, uhf=0,
                                    max_steps=1, timeout=150, confirm=True)
    # It did not converge in one cycle from there.
    if got.get('ok'):
        # ORCA converged after all on this machine -- then the point of the
        # test does not arise, and a converged saddle is a fine outcome.
        assert got.get('imaginary') is not None
        return
    # Not converged, but now characterised: the verdict and the modes are
    # there, and the line names what it reached rather than only that it did
    # not converge.
    assert got.get('verdict') is not None, got.get('status')
    assert got.get('imaginary') is not None
    assert 'did not converge' in got['status']
    assert 'What it reached is' in got['status']
