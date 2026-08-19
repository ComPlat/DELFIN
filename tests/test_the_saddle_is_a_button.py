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
"""

import pytest

from delfin.dashboard import saddle
from editor_source import EDITOR_SOURCE

_needs_orca = pytest.mark.skipif(saddle.find_orca() is None,
                                 reason="ORCA not installed")

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
    """
    assert saddle.SADDLE_METHODS == {'gfn2': 'XTB2', 'gfn1': 'XTB1',
                                     'gfnff': 'XTBFF'}
    # Anything with a basis set is refused here and named as a job.
    refused = saddle.optimise_to_saddle(_ESTIMATE, 'PBE0')
    assert not refused['ok']
    assert 'runs on xtb through ORCA' in refused['status']
    # And nothing to work on is said rather than run.
    assert not saddle.optimise_to_saddle('', 'gfn2')['ok']


def test_what_it_reached_is_read_from_the_last_hessian():
    """A saddle search recalculates its Hessian as it goes, so the output holds
    several: the last describes the geometry it ended on and the earlier ones
    describe geometries it has left.

    Read the first -- which is what a naive parse does -- an OptTS that
    improved its structure reports the mode it started with, and the number on
    screen is about a structure that is no longer anywhere.
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
    assert modes == {'count': 1, 'modes': [-372.26]}, modes
    # Nothing to read is None, not an empty answer that reads as "a minimum".
    assert saddle._last_modes('nothing here') is None


@_needs_orca
def test_the_estimate_becomes_a_converged_transition_state():
    """Measured: the estimate comes in with its two forming bonds at 2.524 and
    2.520 A and one imaginary mode at -131.4 cm-1, and comes out symmetric at
    2.315 and 2.315 with the mode at -372.

    Seven seconds, which is why this is a button.
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

    One core, two gigabytes and three minutes -- welcome exactly as long as it
    is small and short.  A run that wants more than that is a job, and saying
    so is better than taking the node.
    """
    import inspect

    source = inspect.getsource(saddle)
    assert 'cores: int = 1,' in source
    assert 'timeout: Optional[float] = 180.0' in source
    assert "'%maxcore 2000\\n'" in source
    assert 'nprocs {max(1, int(cores))}' in source
    assert 'submit it as a job for more' in source
    # A saddle search needs a Hessian to know which mode to climb, and
    # recalculates it because the mode changes character as it moves.
    assert "'  Calc_Hess true\\n'" in source
    assert "'  Recalc_Hess 5\\n'" in source


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
    assert "'Optimised to a transition state'" in source
    # What it reached is said, whichever it is.
    assert 'def _said_modes(shape, what):' in source
    assert 'is a minimum, not a transition state' in source
    # And where it goes next.
    assert 'Refine it with OPTTS in the ORCA Builder at ' in source
