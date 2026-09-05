"""A transition state under the most accurate method the editor offers.

g-xTB approximates wB97M-V/def2-TZVPPD, and it is the method somebody chooses
when the transition state has to be worth quoting.  Under it there was no
route to a saddle at all: ORCA has keywords for ``XTB1``, ``XTB2`` and
``XTBFF`` and drives them through its own bundled xtb, and g-xTB is a
statically linked build of its own that no ORCA keyword names -- while an
ordinary xtb *accepts* ``--gxtb`` and silently runs GFN2, so there is nothing
to point an existing keyword at either.  Both buttons were hidden under it,
honestly, and the capability was missing.

ORCA publishes an interface for exactly this and it turns out to work:
**ExtOpt**.  ORCA writes a geometry and a request, calls whatever
``EXTOPTEXE`` names, reads an energy and a gradient back, and its own
optimiser does the rest -- the same shape GRRM uses to drive xtb through
``%link=non-supported``.  :mod:`delfin.dashboard.gxtb_engrad` is what answers
it with g-xTB.

Measured here on the sixteen-atom Diels-Alder estimate, ORCA 6.1.1:

* ``! ExtOpt OptTS`` with ``Calc_Hess true`` alone stops in PROPINT --
  "ERROR (SHARK): Failed to read input file" -- because ORCA's analytic
  Hessian has no basis set to work from.  ``NumHess true`` is ORCA's own
  answer to that and is what is asked for here.
* one numerical Hessian with Bofill updates after it is *not* enough: the two
  forming bonds walked from 2.524 A out to 3.03 and 3.05, the run hit its
  sixty-cycle bound, and what it left had no imaginary mode -- it had gone
  back down to the van-der-Waals complex.  So ``Recalc_Hess`` stays at five,
  the same cadence the xtb methods get.
* a g-xTB gradient is a whole process and that is where the cost is: on water,
  0.114 s for the process and 0.004 s of SCF inside it.  Which is also why
  this is a button and not the interactive climb next door -- see
  :func:`test_the_climb_is_still_the_three_xtb_methods_and_the_numbers_say_why`.
"""

import os
import subprocess
import tempfile
from pathlib import Path

import pytest

from delfin.dashboard import climb
from delfin.dashboard import gfn_optimize as gfn
from delfin.dashboard import gxtb_engrad
from delfin.dashboard import saddle

_needs_gxtb = pytest.mark.skipif(gfn.find_gxtb() is None,
                                 reason='g-xTB not installed')
_needs_orca = pytest.mark.skipif(saddle.find_orca() is None,
                                 reason='ORCA not installed')

#: The transition state the path finder estimated for the Diels-Alder, the
#: same sixteen atoms the rest of the saddle suite is built on.
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

#: And the saddle ORCA's OptTS on GFN2 converges to from it: the two forming
#: bonds at 2.3153 A, one imaginary mode at -393.5 cm-1 under GFN2 and -356.8
#: under g-xTB.  Written out here because it is the structure a g-xTB search
#: has somewhere to go from -- see
#: :func:`test_a_gfn2_estimate_is_not_where_a_g_xtb_saddle_is`.
_GFN2_SADDLE = """16
the Diels-Alder saddle, as ORCA on GFN2 converged it
  C          -1.43478368987752     -0.11555756879445      0.39801753593952
  C          -0.70800000464805      0.96029290831534      0.01351101273219
  C           0.70792543036838      0.96066325147959      0.01382699497615
  C           1.43510098205132     -0.11480687774120      0.39865617430935
  H          -2.50480092399590     -0.04135148267412      0.50582952622211
  H          -1.04673865581762     -1.11397983255077      0.31729931422325
  H          -1.20115753813639      1.91631252668784     -0.10370598026703
  H           1.20063543422156      1.91694042137579     -0.10317121370029
  H           2.50503091371612     -0.04004089982782      0.50694578141924
  H           1.04761410022072     -1.11343199175707      0.31776831289261
  C          -0.67562925889091     -0.26587715375125      2.58016852501954
  C           0.67503524104817     -0.26553006909885      2.58047503090094
  H          -1.22616409778789      0.63483074834021      2.78770928087501
  H          -1.23028498551070     -1.18636553718980      2.63888406495792
  H           1.22501381825542      0.63546069528911      2.78826237108677
  H           1.23013805516529     -1.18573310998941      2.63943734003727
"""

_WATER = """3
water
O   0.00000000  0.00000000  0.11779000
H   0.00000000  0.75545000 -0.47116000
H   0.00000000 -0.75545000 -0.47116000
"""


def _request(folder, xyz=_WATER, charge=0, multiplicity=1, cores=4,
             gradient=1, charges=''):
    """What ORCA writes before it calls the external program, as it writes it."""
    folder = Path(folder)
    (folder / 'in.xyz').write_text(xyz, encoding='utf-8')
    path = folder / 'in_EXT.extinp.tmp'
    path.write_text(
        f'in.xyz # xyz filename\n{charge} # charge\n'
        f'{multiplicity} # multiplicity\n{cores} # NCores\n'
        f'{gradient} # do gradient\n{charges} # point charge filename\n',
        encoding='utf-8')
    return path


# -- the table --------------------------------------------------------------

def test_orca_drives_g_xtb_through_the_interface_it_has_no_keyword_for():
    """The three it knows by name, and the one it is told how to ask for.

    ORCA runs GFN2, GFN1 and GFN-FF through its own bundled ``otool_xtb``.
    g-xTB is not in it and cannot be added to it: the method ships as its own
    xtb build, and the ordinary one accepts ``--gxtb`` and answers with GFN2
    -- measured elsewhere in this suite, to the last digit.  ExtOpt is the
    published way round that, and it is what this entry means.
    """
    assert saddle.SADDLE_METHODS['gfn2'] == 'XTB2'
    assert saddle.SADDLE_METHODS['gfn1'] == 'XTB1'
    assert saddle.SADDLE_METHODS['gfnff'] == 'XTBFF'
    assert saddle.SADDLE_METHODS['gxtb'] == 'ExtOpt'
    # And the ones ORCA drives itself are named, because they are the ones
    # that take no program of ours and no hook.
    assert 'gxtb' not in saddle._DRIVEN_BY_ORCAS_XTB
    # Anything with a basis set is still a job rather than a press.
    refused = saddle.optimise_to_saddle(_ESTIMATE, 'PBE0')
    assert not refused['ok']
    assert 'runs on xtb through ORCA' in refused['status']


def test_a_saddle_in_a_solvent_is_refused_under_g_xtb_rather_than_dropped(
        monkeypatch):
    """This build of g-xTB has no implicit solvation, and silence would lie.

    ALPB, GBSA and CPCM-X all stop it; ddCOSMO runs and only writes a file.
    Passing the solvent through would stop ORCA in the middle of a search;
    quietly dropping it would hand back a gas-phase saddle under a solvent
    nobody removed, which is the worse of the two.
    """
    monkeypatch.setattr(gfn, 'find_binary', lambda method=None: '/somewhere/xtb-gxtb')
    refused = saddle.optimise_to_saddle(_ESTIMATE, 'gxtb', solvent='water')
    assert not refused['ok']
    assert 'no implicit solvation' in refused['status']
    assert 'GFN2' in refused['status']          # and where to go instead


def test_a_missing_g_xtb_is_named_before_orca_is_started(monkeypatch):
    """A program that is not there is worth saying so about, not running into."""
    monkeypatch.setattr(gfn, 'find_binary', lambda method=None: None)
    refused = saddle.optimise_to_saddle(_ESTIMATE, 'gxtb')
    assert not refused['ok']
    assert 'was not found' in refused['status']


# -- the hook ---------------------------------------------------------------

def test_the_request_is_read_past_orcas_own_comments(tmp_path):
    """Six values, each with ORCA's explanation written after a hash."""
    asked = gxtb_engrad.read_request(
        _request(tmp_path, charge=-1, multiplicity=2, cores=7, gradient=0))
    assert asked == {'xyz': 'in.xyz', 'charge': -1, 'multiplicity': 2,
                     'cores': 7, 'gradient': False, 'point_charges': ''}


def test_the_answer_is_in_orcas_own_engrad_format():
    """Atom count, energy, gradient one component to a line, atoms in Bohr.

    ORCA checks the count against its own and stops when they disagree, so
    that is the one mistake this could make quietly -- and it cannot, because
    the count is written from the list rather than carried in beside it.
    """
    text = gxtb_engrad.engrad_document(
        ['O', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        -76.4374757521, [[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]])
    rows = [r for r in text.splitlines() if not r.startswith('#')]
    assert rows[0].strip() == '2'
    assert float(rows[1]) == pytest.approx(-76.4374757521)
    assert [float(r) for r in rows[2:8]] == [0.0, 0.0, -0.5, 0.0, 0.0, 0.5]
    first, second = rows[8].split(), rows[9].split()
    assert int(first[0]) == 8 and int(second[0]) == 1
    # One Angstrom apart in z, which is 1.889726 Bohr.
    assert float(second[3]) == pytest.approx(
        1.0 / gxtb_engrad.BOHR_IN_ANGSTROM, abs=1e-9)


def test_the_hook_runs_this_interpreter_and_not_whatever_is_on_the_path(tmp_path):
    """A dashboard is very often in an environment the PATH does not lead to.

    ``EXTOPTEXE`` has to name something executable, so a two-line shell script
    is written into the run's own directory -- this run's, not a shared one,
    so two searches at once cannot collide.
    """
    import sys

    hook = gxtb_engrad.write_hook(tmp_path)
    assert hook.parent == tmp_path
    assert os.access(str(hook), os.X_OK)
    said = hook.read_text()
    assert sys.executable in said
    assert 'delfin.dashboard.gxtb_engrad' in said


def test_what_it_cannot_do_is_said_rather_than_answered_anyway(tmp_path):
    """ORCA sees an exit code; the reason has to be written somewhere.

    Two ways this can be asked for something it has not got, and neither may
    come back as an energy.  External point charges are the QM/MM case -- g-xTB
    through here has no way to apply them, and an answer that ignored them
    would be an answer about a different system.  A missing binary is the
    other, and it is not the same as a broken one: g-xTB is a build of its
    own, and the xtb standing next to it would answer with GFN2.
    """
    with_charges = _request(tmp_path, charges='points.pc')
    assert gxtb_engrad.answer(with_charges) != 0
    assert not (tmp_path / 'in_EXT.engrad').is_file()

    plain = _request(tmp_path)
    assert gxtb_engrad.answer(plain, binary='/nowhere/xtb-gxtb') != 0
    assert not (tmp_path / 'in_EXT.engrad').is_file()

    # And a request that is not one at all.
    broken = tmp_path / 'nonsense_EXT.extinp.tmp'
    broken.write_text('in.xyz # xyz filename\n', encoding='utf-8')
    with pytest.raises(ValueError):
        gxtb_engrad.read_request(broken)


@_needs_gxtb
def test_the_hook_answers_with_the_numbers_g_xtb_itself_gives(tmp_path):
    """Driven, not read: the energy and gradient are compared with a run of its own.

    The energy comes out of the Turbomole ``gradient`` file rather than the
    printed one, which carries fewer digits -- so agreement here is to the
    twelfth decimal rather than to the sixth.
    """
    answered = gxtb_engrad.answer(_request(tmp_path))
    assert answered == 0
    written = (tmp_path / 'in_EXT.engrad').read_text()
    rows = [r for r in written.splitlines() if not r.startswith('#')]
    energy = float(rows[1])
    gradient = [float(r) for r in rows[2:11]]

    mine = Path(tempfile.mkdtemp())
    (mine / 'in.xyz').write_text(_WATER, encoding='utf-8')
    room = dict(os.environ, OMP_NUM_THREADS='4', MKL_NUM_THREADS='4')
    subprocess.run([gfn.find_gxtb(), 'in.xyz', '--gxtb', '--grad'],
                   cwd=str(mine), env=room, capture_output=True, text=True,
                   timeout=600)
    theirs, rows_of_it = gxtb_engrad._read_turbomole_gradient(
        mine / 'gradient', 3)
    assert energy == pytest.approx(theirs, abs=1e-9)
    flat = [one for row in rows_of_it for one in row]
    assert gradient == pytest.approx(flat, abs=1e-9)
    # And it is g-xTB rather than the GFN2 an ordinary xtb would answer with:
    # water is about -76.4 Hartree under g-xTB and about -5.07 under GFN2.
    assert energy < -50.0


@_needs_gxtb
def test_the_hook_runs_through_a_shell_the_way_orca_starts_it(tmp_path):
    """ORCA does not import it; it executes it, and passes one argument."""
    request = _request(tmp_path)
    hook = gxtb_engrad.write_hook(tmp_path)
    done = subprocess.run([str(hook), request.name], cwd=str(tmp_path),
                          capture_output=True, text=True, timeout=600)
    assert done.returncode == 0, done.stderr[-2000:]
    assert (tmp_path / 'in_EXT.engrad').is_file()


# -- the installer's own account of itself ---------------------------------

def test_the_installer_does_not_call_a_g_xtb_it_just_installed_absent(tmp_path):
    """A tool is not always called what its binary is called.

    g-xTB is asked for as "gxtb" and installed as "xtb-gxtb", deliberately:
    an ordinary xtb accepts --gxtb and silently runs GFN2, so the two builds
    must never be confused.  The closing summary looked the tool up under the
    name it was asked for, found nothing in bin/, and printed "gxtb absent"
    -- after downloading it, checking its sha256, unpacking it and linking it.
    Which is the one line a user reads.

    The function is run rather than the file read, on a bin/ made here.
    """
    script = gfn.install_script()
    assert script is not None
    lifted = subprocess.run(
        ['sed', '-n', '/^version_of() {/,/^}/p', str(script)],
        capture_output=True, text=True, check=True).stdout
    assert 'gxtb|g-xtb' in lifted, lifted

    binary = tmp_path / 'xtb-gxtb'
    binary.write_text('#!/bin/sh\necho "   * xtb version 6.7.1 (26dd68d)"\n',
                      encoding='utf-8')
    binary.chmod(0o755)
    said = subprocess.run(
        ['bash', '-c', f'{lifted}\nBIN_DIR="{tmp_path}"\nversion_of gxtb'],
        capture_output=True, text=True)
    assert said.stdout.strip() == '6.7.1', (said.stdout, said.stderr)
    # And one that really is absent still says so.
    missing = subprocess.run(
        ['bash', '-c', f'{lifted}\nBIN_DIR="{tmp_path / "empty"}"\n'
                       'version_of gxtb'],
        capture_output=True, text=True)
    assert missing.stdout.strip() == 'absent'


# -- what the climb was not given, and why ---------------------------------

def test_the_climb_is_still_the_three_xtb_methods_and_the_numbers_say_why():
    """The interactive climb is not extended to g-xTB, and this records why.

    The climb walks on gradients and needs them at the speed of a hand.  Its
    own budget, measured on this box under GFN2: a bare step is 9.1 ms at
    sixteen atoms and 116.5 ms at fifty, and one exact Hessian -- 6N gradients
    -- is 0.55 s at sixteen and 16.3 s at fifty, taken every twenty steps.

    g-xTB cannot be asked for a gradient that way at all.  It has no entry in
    the xtb Python module's ``Param`` -- the library the fast path uses is the
    ordinary xtb's -- so every gradient is a separate process, and a process
    is where the cost is: on water, 0.114 s for the whole run and 0.004 s of
    SCF inside it.  Measured on the same structures the climb was measured on,
    threads pinned to four, both methods through the command line and both
    keeping their directory so each has its own ``xtbrestart``:

        atoms   16      20      33      50
        GFN2    0.050   0.065   0.076   0.124  s
        g-xTB   0.186   0.199   0.300   0.543  s

    -- three to four times, process for process.  Against what the climb
    actually runs, which is GFN2 *in* process, one g-xTB gradient at sixteen
    atoms is 0.29 s against 6 ms: thirty times, or three frames a second where
    a drag wants ten.  And the Hessian settles it: 6N gradients at sixteen
    atoms is 96 processes, 17.9 s measured, against 0.55 -- so a climb that
    went nowhere would spend its eight Hessians on two and a half minutes at
    sixteen atoms and twenty at fifty.

    So the honest answer is that g-xTB is a press and not a drag, and that is
    what it was given: ORCA's OptTS through ExtOpt, above.  A hand that wants
    to push a structure towards a saddle can still do it under GFN2 and hand
    what it reaches to g-xTB, which is a single point of a fifth of a second.
    """
    assert set(climb.CLIMB_METHODS) == {'gfn2', 'gfn1', 'gfnff'}
    refused = climb.climb_to_saddle(_ESTIMATE, 'gxtb')
    assert not refused['ok']
    assert 'runs on xtb' in refused['status']


@_needs_gxtb
def test_a_gfn2_estimate_is_not_where_a_g_xtb_saddle_is():
    """Which is the reason a saddle under g-xTB is a second press and not a
    relabelling of the first.

    The path finder's estimate for this Diels-Alder is a saddle under GFN2 --
    xtb's own Hessian there has one mode at -131.4 cm-1 -- and under g-xTB the
    Hessian at the *same geometry* has no negative eigenvalue at all.  An
    optimiser started there under g-xTB does not fail; it correctly slides
    down, and measured, it converged in 268 s onto a structure with the two
    forming bonds at 3.48 and 3.51 A and no imaginary mode: the van-der-Waals
    complex.

    Where g-xTB does have one is at the saddle GFN2 itself converges to --
    2.3153 A on both forming bonds -- where its Hessian says -356.8 cm-1.  So
    the chain that works is GFN2 first and g-xTB after it, which is also the
    cheap way round: the GFN2 half is seconds.
    """
    at_the_estimate = gfn.optimize_with_gfn(
        _ESTIMATE, 'gxtb', optimise=False, free_energy=True, timeout=600)
    assert at_the_estimate.get('ok'), at_the_estimate.get('status')
    assert at_the_estimate['imaginary']['count'] == 0, \
        at_the_estimate['imaginary']

    at_the_saddle = gfn.optimize_with_gfn(
        _GFN2_SADDLE, 'gxtb', optimise=False, free_energy=True, timeout=600)
    assert at_the_saddle.get('ok'), at_the_saddle.get('status')
    assert at_the_saddle['imaginary']['count'] == 1, at_the_saddle['imaginary']
    # Around -357, and the digits move under threading, so it is bounded
    # rather than named.
    assert -500.0 < at_the_saddle['imaginary']['modes'][0] < -200.0


@pytest.mark.slow
@_needs_gxtb
@_needs_orca
def test_orca_climbs_to_a_saddle_on_g_xtb_gradients():
    """End to end, and what came back is asked of a Hessian rather than assumed.

    From the saddle GFN2 reached, because that is where g-xTB has one -- see
    the test above.  Measured: 101 s on sixteen atoms, against 97 s for the
    same optimiser on GFN2 taken on the same box within the hour, because
    every gradient is a separate process and every numerical Hessian is
    forty-odd of them.  A press for a molecule this size and a job above it,
    which is what the timeout is for.
    """
    got = saddle.optimise_to_saddle(_GFN2_SADDLE, 'gxtb', cores=8,
                                    timeout=1200.0)
    assert got.get('xyz'), got.get('status')
    assert got['ok'], got['status']
    # One mode going the wrong way, and it is g-xTB's own Hessian that says
    # so -- not the optimiser's last block, which belongs to a geometry it has
    # since left.
    assert got['confirmed'] is True
    assert got['imaginary']['count'] == 1, got['imaginary']
    assert got['imaginary']['modes'][0] < -40.0
    assert 'g-xTB' in got['status']
