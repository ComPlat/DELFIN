"""Optimise a structure to a transition state, fast enough to press a button.

xtb has no saddle-point optimiser -- 6.7.1 offers ``--hess``, ``--ohess``,
``--bhess`` and ``--path``, and nothing that climbs to a first-order saddle.
ORCA has one, and ORCA can be told to take its gradients from xtb.  That pair
is the whole of this file: ORCA's optimiser, xtb's speed.

Measured on the transition state the path finder estimated for a Diels-Alder,
sixteen atoms: ``! XTB2 OPTTS`` converged in 7 s.  The estimate came in with
its two forming bonds at 2.524 and 2.520 A and one imaginary mode at
-131.4 cm-1, and came out symmetric at 2.315 and 2.315 with the mode at
-393.5.  Seven seconds is a button; an OptTS on a real basis set is a job, and
belongs in the ORCA Builder where jobs are submitted.

The -372 written here before was the optimiser's own last Hessian block, and
two runs of the same input reported -372 and -301.55 for the same converged
geometry -- see :func:`_last_modes` for why, and why the frequency is now
taken on the structure that came back.

So the editor's chain runs end to end without leaving it: a scan walks the
reaction, the path finder finds a way between its two ends and estimates the
saddle, and this sharpens the estimate into a converged transition state --
with the one imaginary frequency that says it is one.  The last two of those
are one press in :func:`path_to_saddle`, which is twelve seconds and lands
within a wavenumber of a nudged elastic band that takes seven minutes.

And says when it is not one, which is the half that cannot be left out.  A
saddle search does not fail when it goes wrong; it succeeds at arriving
somewhere.  Measured here: the same chain, given a sloppy pair of ends,
converged and reported a stationary point that a Hessian found two modes
going the wrong way at.  :func:`verdict` names what was reached -- a minimum,
a transition state, a second-order saddle -- against thresholds taken from
searches that have had to make this judgement at scale, each cited where it is
defined.

The climb is watched while it happens rather than reported at the end.  ORCA
appends every accepted step to ``<base>_trj.xyz`` with its energy on the
comment line, so the path is readable as it is walked: the viewer can show
where the structure is going, and a climb that is heading somewhere the user
did not mean can be stopped instead of waited out.  A saddle search is the
one optimisation where watching matters most -- it is climbing, so a wrong
mode does not fail, it succeeds at reaching the wrong place.

It runs where the dashboard runs, which on a cluster is the login node, and
that is a place with rules.  So: a few cores, two gigabytes, and a timeout in
minutes rather than hours -- and only the methods that can finish inside
that.  Four is the most a login node is usually happy to lend; a saddle
search on a real basis set is a job, and the ORCA Builder is where jobs are
submitted.

How many cores is worth asking for was measured rather than assumed, and
the answer depends entirely on the system.  On the sixteen-atom estimate, 1,
2, 4, 8 and 16 cores took 6.7, 7.4, 6.7, 6.4 and 6.7 s -- flat, because a
system that small is over before any of it can be shared out.  On a
fifty-one-atom molecule, where the Hessian is the cost and grows with the
square of the atoms, the same counts took 34.5, 12.1, 7.9, 6.3 and 5.9 s.

So: eight.  Four is not wrong, but it leaves most of the gain on the table
on exactly the systems where the wait is long enough to mind -- and past
eight the curve has flattened, so the extra cores buy a fifth of a speedup
and cost a node other people are also using.  Never more than the machine
has, which matters on a small login node rather than on this one.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from . import gfn_optimize as _gfn

#: The methods ORCA can drive that are fast enough for this to be a button.
#: Anything with a basis set is a job, not a press.
SADDLE_METHODS = {
    'gfn2': 'XTB2',
    'gfn1': 'XTB1',
    'gfnff': 'XTBFF',
}

#: How often the trajectory on disk is looked at while the climb runs.
#:
#: A step of an xtb-driven OptTS takes a good deal longer than this even on a
#: small system, so nothing is missed by waiting; and the file is only read
#: when it has grown, so a look that finds nothing new costs a stat.
WATCH_EVERY = 0.12

#: ORCA says this when the optimiser is finished, and only then.
_DONE_RE = re.compile(r'\*\*\*\s*OPTIMIZATION RUN DONE\s*\*\*\*')
#: The frequencies come from the Hessian OptTS needs anyway, printed by the
#: engine underneath.  The last block is the one that describes what was
#: reached rather than what was started from.
_EIGVAL_RE = re.compile(r'eigval\s*:\s*(.+)')
_IMAGINARY_RE = re.compile(r'#\s*imaginary freq\.\s+(\d+)')
#: Every frame of the trajectory carries its energy, on the comment line:
#: ``Coordinates from ORCA-job in E -17.813830573670``.
_TRJ_ENERGY_RE = re.compile(r'E\s+(-?\d+\.\d+)')


def find_orca() -> Optional[str]:
    """Where ORCA is, through the same resolver the rest of DELFIN uses."""
    try:
        from ..qm_runtime import resolve_tool
    except Exception:
        return shutil.which('orca')
    try:
        found = resolve_tool('orca')
    except Exception:
        found = None
    return found.path if found else shutil.which('orca')


def _last_modes(output: str) -> Optional[Dict[str, Any]]:
    """The modes of the last Hessian a saddle search took, whichever that was.

    A saddle search recalculates its Hessian as it goes, so the output holds
    several; the last one is nearer the geometry the run ended on than the
    earlier ones, which are about geometries it has left.  Read the first --
    which is what a naive parse does -- an OptTS that improved its structure
    reports the mode it started with.

    The last one is still not the geometry that was reached, and the
    difference is not small.  Measured on the sixteen-atom Diels-Alder: the
    run recalculated at cycles 1, 6 and 11 and converged at 16, and its last
    block says -301.55 cm-1 while a Hessian taken on the structure it
    actually handed back says -393.53.  Same geometry, two answers, because
    the block belongs to a geometry five steps back.  So this is the fallback
    and :func:`optimise_to_saddle` takes a Hessian of its own -- 0.3 s on
    sixteen atoms, against a number that is about nothing on screen.
    """
    blocks = output.split('projected vibrational frequencies')
    if len(blocks) < 2:
        return None
    seen = []
    for row in _EIGVAL_RE.finditer(blocks[-1]):
        for word in row.group(1).split():
            try:
                seen.append(float(word))
            except ValueError:
                pass
        if len(seen) > 200:
            break
    counted = None
    for found in _IMAGINARY_RE.finditer(output):
        counted = int(found.group(1))
    modes = sorted({one for one in seen if one < -5.0})
    return {'count': counted if counted is not None else len(modes),
            'modes': modes[:4],
            'real': _gfn.lowest_real_modes(blocks[-1])}


#: The least an imaginary mode may be before a saddle is worth the name.
#:
#: autodE keeps ``min_imag_freq = Frequency(-40, units="cm-1")`` in
#: ``autode/config.py`` and refuses a structure whose most negative mode is
#: above it, with the comment: "Minimum magnitude of the imaginary frequency
#: (cm-1) to consider for a 'true' TS. For very shallow saddle points this may
#: need to be reduced to e.g. -10 cm-1. Although most TSs have |v_imag| > 100
#: cm-1 this threshold is designed to be conservative."
#:
#: Conservative is what is wanted here.  This is a semiempirical answer meant
#: to be handed on to a real one, and refusing something that was a saddle
#: costs a press, while accepting something that was not costs a job.
MIN_IMAGINARY = -40.0

#: And what a search that means to be strict asks for.
#:
#: AutoMeKin's ``imagmin`` is "the minimum value for the imaginary frequency
#: (in absolute value and cm-1) of the selected TS structures".  Its code
#: default is 0 -- the 200 is what every input file it ships with sets, and
#: what its own tutorial uses, so it is a recommendation rather than a
#: default and is said here as one.
CLEAR_IMAGINARY = 200.0

#: Two real modes that sum to less than this say the structure is nearly flat
#: in those directions, which is a saddle of higher order in all but name.
#:
#: AutoMeKin's check, with the ``tight_ts`` it has on by default: it refuses a
#: transition state outright when the lowest real frequency is negative
#: (``if(freq[2]<0) f=-1``), and refuses it as well when the two lowest real
#: ones sum to under ten (``else if(freq[2]+freq[3]<10) f=-2``), reported as
#: "Sum of 2 lowest real freqs < 10cm-1".
FLAT_PAIR = 10.0

#: How many modes going the wrong way is a transition state.  Exactly one, and
#: that is not a convention: a first-order saddle is a maximum in one
#: direction and a minimum in all the rest, and a structure that is a maximum
#: in two is not a point any single step passes through.  YARP says it in a
#: line -- ``is_valid_ts`` counts the modes below a numerical tolerance and
#: returns ``n_imaginary == 1``.
FIRST_ORDER = 1

#: What a saddle of each order is called, so a number is never printed where a
#: name would do.  Beyond these it is said as an order.
_ORDER_NAMES = {0: 'a minimum', 1: 'a transition state',
                2: 'a second-order saddle', 3: 'a third-order saddle'}


def _named(order: Optional[int]) -> str:
    if order is None:
        return 'a structure whose modes were not checked'
    return _ORDER_NAMES.get(order) or f'a saddle point of order {order}'


def verdict(shape: Optional[Dict[str, Any]],
            what: str = 'What it reached',
            advise: bool = True) -> Dict[str, Any]:
    """What a Hessian says a structure is, and whether that is a saddle.

    Returns ``{'order', 'first_order', 'name', 'lines'}``.  *order* is how
    many modes go the wrong way, *name* is what a structure with that many is
    called, and *lines* are the sentences to show.

    This exists because a saddle search cannot be trusted blind.  It does not
    fail when it goes wrong -- it succeeds at arriving somewhere, and a
    converged run that reports a stationary point may have converged onto a
    structure that is a maximum in two directions at once.  Measured on this
    machine: the chain from a path finder to ORCA's OptTS, given a deliberately
    sloppy pair of ends, converged in 14 s and reported a stationary point,
    and a Hessian on what it handed back found two modes going the wrong way
    at -52 and -26 cm-1.  Nothing in the run said so.  This does.

    A structure with two is named as a second-order saddle rather than
    reported as done, because the difference is not a detail of quality: a
    first-order saddle is a point a path can pass through, and a second-order
    one is not a point anything passes through.

    The thresholds are other people's, and each is cited where it is defined:
    see :data:`MIN_IMAGINARY`, :data:`CLEAR_IMAGINARY`, :data:`FLAT_PAIR` and
    :data:`FIRST_ORDER`.  They are stated as magnitudes and orders and never
    as anything about a particular kind of chemistry -- every sort of system
    is computed here, and a threshold that assumed one of them would be wrong
    about the rest.
    """
    if not shape:
        return {'order': None, 'first_order': False, 'name': _named(None),
                'lines': [f'{what} could not be checked for imaginary '
                          'modes.']}
    modes = [float(one) for one in (shape.get('modes') or [])]
    real = [float(one) for one in (shape.get('real') or [])]
    order = shape.get('count')
    order = len(modes) if order is None else int(order)
    # Only the ones the engine counted, because they are what the order is
    # being said from; whatever is below its cutoff is said separately and as
    # what it is.
    listed = ', '.join(f'{one:.0f}' for one in modes[:max(0, order)])
    lines: List[str] = []

    if order == FIRST_ORDER:
        deepest = modes[0] if modes else None
        if deepest is not None and deepest > MIN_IMAGINARY:
            # Shallow enough that the mode may be the method's noise rather
            # than the structure's shape.  Named as what it is instead of
            # counted as a saddle: autodE refuses exactly this.
            lines.append(
                f'{what} has one mode going the wrong way, at '
                f'{deepest:.0f} cm-1, and no others -- but that is shallower '
                f'than the {MIN_IMAGINARY:.0f} cm-1 autodE takes as the least '
                'it will call a transition state, so it is not one on this '
                'evidence.')
            first = False
        else:
            lines.append(
                f'{what} is a transition state: one mode goes the wrong way'
                + (f', at {deepest:.0f} cm-1' if deepest is not None else '')
                + ', and no others.')
            first = True
            if deepest is not None and abs(deepest) < CLEAR_IMAGINARY:
                lines.append(
                    'It is a soft one: the inputs AutoMeKin ships ask for '
                    f'{CLEAR_IMAGINARY:.0f} cm-1 before they will take a '
                    'saddle at all.')
    elif order == 0:
        lines.append(f'{what} is a minimum, not a transition state: no mode '
                     'goes the wrong way.')
        first = False
    else:
        lines.append(
            f'{what} is {_named(order)}, not a transition state: {order} '
            'modes go the wrong way'
            + (f' ({listed} cm-1)' if listed else '')
            + f'. A transition state has exactly {FIRST_ORDER}. A structure '
              'that is a maximum in more than one direction is not a point '
              'any single step passes through.')
        if advise:
            # Only where something has converged.  Told about an estimate --
            # which is nobody's stationary point -- it would be wrong twice
            # over: there is nothing to converge again from, and a climb from
            # an estimate with two modes going the wrong way may perfectly
            # well end on a transition state.
            lines.append(
                'Move it a little along the second of those modes and climb '
                'again, or start the search from a different structure -- '
                'converging again from here will only find this one.')
        first = False

    # The modes the engine did not count, and the directions that are nearly
    # not directions.  Both are the same warning: what is around the saddle
    # is as much of the answer as the mode being climbed.
    uncounted = modes[order:] if 0 <= order < len(modes) else []
    if uncounted:
        lines.append(
            'Below what the engine counts, '
            + ('one mode is' if len(uncounted) == 1
               else f'{len(uncounted)} more modes are')
            + ' still on the wrong side of zero ('
            + ', '.join(f'{one:.0f}' for one in uncounted)
            + ' cm-1), so this is within reach of a saddle one order higher.')
    if first and len(real) > 1 and (real[0] + real[1]) < FLAT_PAIR:
        lines.append(
            f'Its two softest real modes are {real[0]:.1f} and {real[1]:.1f} '
            f'cm-1, summing to under the {FLAT_PAIR:.0f} cm-1 AutoMeKin '
            'refuses -- flat enough in those directions to be a saddle of '
            'higher order in all but name.')
    return {'order': order, 'first_order': first, 'name': _named(order),
            'lines': lines}


def frames_in(text: str, atoms: int) -> Tuple[List[List[float]],
                                              List[Optional[float]]]:
    """The whole frames of a trajectory, and the energy each was written with.

    Read while ORCA is still appending, so a read lands in the middle of a
    block as often as not.  Only complete blocks are returned and the rest is
    left for the next look: half a geometry drawn is worse than a geometry
    drawn a tenth of a second later.
    """
    rows = text.splitlines()
    step = int(atoms) + 2
    if step < 3:
        return [], []
    walked: List[List[float]] = []
    spent: List[Optional[float]] = []
    for start in range(0, len(rows) - step + 1, step):
        head = rows[start].strip()
        if not head.isdigit() or int(head) != int(atoms):
            break
        flat: List[float] = []
        for line in rows[start + 2:start + step]:
            bits = line.split()
            if len(bits) < 4:
                break
            try:
                flat.extend(float(one) for one in bits[1:4])
            except ValueError:
                break
        if len(flat) != int(atoms) * 3:
            break
        walked.append(flat)
        found = _TRJ_ENERGY_RE.search(rows[start + 1])
        spent.append(float(found.group(1)) if found else None)
    return walked, spent


def _share(cores: Any) -> int:
    """How many cores to ask for, never more than there are.

    Asking for more than the machine has does not make it faster; it makes
    ORCA start processes that fight each other for the same cores, and on a
    small login node that is felt by everybody else on it.
    """
    try:
        want = int(cores)
    except (TypeError, ValueError):
        want = 1
    return max(1, min(want, os.cpu_count() or 1))


def _stop(running: subprocess.Popen) -> None:
    """Ask ORCA to stop, and insist if it does not.

    Terminate first, because a run that is given the chance closes its files
    and leaves the geometry it reached behind; kill only what will not go.
    """
    try:
        running.terminate()
    except OSError:
        return
    try:
        running.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:
            running.kill()
            running.wait(timeout=10)
        except (OSError, subprocess.TimeoutExpired):
            pass


def optimise_to_saddle(xyz_text: str, method: str = 'gfn2', *,
                       charge: int = 0, uhf: int = 0,
                       solvent: Optional[str] = None,
                       max_steps: int = 60,
                       cores: int = 8,
                       timeout: Optional[float] = 180.0,
                       confirm: bool = True,
                       on_frame: Optional[Callable[
                           [Sequence[Sequence[float]],
                            Sequence[Optional[float]]], None]] = None,
                       should_stop: Optional[Callable[[], bool]] = None,
                       ) -> Dict[str, Any]:
    """Climb to the nearest first-order saddle, and say whether one was found.

    Returns ``{'ok', 'xyz', 'imaginary', 'verdict', 'seconds', 'status',
    'path', 'energies'}``.  *imaginary* is what a Hessian says about the
    structure that was handed back: one mode going the wrong way and no others
    is a transition state, and anything else is not -- which is the whole
    reason to run this rather than trust an estimate.  *verdict* is
    :func:`verdict` on it, which names what was reached rather than reporting
    that something was.

    With *confirm* the Hessian is taken on the geometry the run ended on
    rather than read out of the optimiser's own output, and that is not
    pedantry.  An OptTS recalculates its Hessian every few steps and converges
    somewhere between two of them, so the last block in its output belongs to
    a geometry it has since left: measured on the sixteen-atom Diels-Alder,
    the last block says -301.55 cm-1 and a Hessian on what came back says
    -393.53 -- the same structure, and the second number is the one about it.
    0.3 s there, and it is the difference between a frequency about the
    structure on screen and a frequency about nothing.

    A saddle search needs a Hessian to know which mode to climb, so it asks
    for one and recalculates it every five steps: the mode being followed
    changes character as the structure moves, and one taken at the start stops
    describing it.

    *on_frame* is handed the path so far each time it grows, with the energy
    of every step, so the climb can be watched rather than waited for.
    *should_stop* is asked as often, and a true answer ends the run and keeps
    what it had reached -- a climb is the one optimisation where watching
    earns its keep, because a wrong mode does not fail, it succeeds at
    arriving somewhere nobody wanted.

    Eight cores and two gigabytes by default, and three minutes.  This runs
    where the dashboard runs, and on a cluster that is the login node -- a
    place where a calculation is welcome exactly as long as it is small and
    short.  Eight is measured to be where the gain stops being worth the
    node: on fifty-one atoms it is 4.4 times one core, and sixteen is only
    5.5.  A run that wants more than that is a job, and saying so is better
    than taking the node.
    """
    keyword = SADDLE_METHODS.get(str(method or '').lower())
    if keyword is None:
        return {'ok': False,
                'status': (f'A saddle search here runs on xtb through ORCA, '
                           f'and {method} is not one of '
                           f'{", ".join(sorted(SADDLE_METHODS))}.')}
    binary = find_orca()
    if binary is None:
        return {'ok': False,
                'status': ('A transition-state optimisation needs ORCA, which '
                           'was not found. xtb has no saddle-point optimiser '
                           'of its own.')}
    body = _gfn.atom_lines(xyz_text)
    if not body:
        return {'ok': False, 'status': 'There is no structure to work on.'}

    started = time.perf_counter()
    folder = Path(tempfile.mkdtemp(prefix='delfin-saddle-'))
    walked: List[List[float]] = []
    spent: List[Optional[float]] = []
    try:
        (folder / 'in.xyz').write_text(
            f'{len(body)}\nfrom the DELFIN viewer\n' + '\n'.join(body) + '\n',
            encoding='utf-8')
        wet = f' ALPB({solvent})' if solvent else ''
        (folder / 'in.inp').write_text(
            f'! {keyword} OPTTS{wet}\n'
            f'%pal\n  nprocs {_share(cores)}\nend\n'
            '%maxcore 2000\n'
            '%geom\n'
            '  Calc_Hess true\n'
            '  Recalc_Hess 5\n'
            f'  MaxIter {max(5, int(max_steps))}\n'
            'end\n'
            f'* xyzfile {int(charge)} {max(0, int(uhf)) + 1} in.xyz\n',
            encoding='utf-8')
        environment = dict(os.environ)
        # ORCA calls its own tools by name, so where it lives has to be on the
        # path for the run -- otherwise it finds its optimiser and not its
        # xtb interface, and stops with a message about neither.
        environment['PATH'] = (str(Path(binary).parent) + os.pathsep
                               + environment.get('PATH', ''))
        trail = folder / 'in_trj.xyz'
        log = folder / 'out.log'
        halted = False
        ran_out = False

        def _look() -> None:
            """Read the path so far, and hand on whatever is new."""
            nonlocal walked, spent
            if not trail.is_file():
                return
            try:
                text = trail.read_text(encoding='utf-8', errors='ignore')
            except OSError:
                return
            fresh, costs = frames_in(text, len(body))
            if len(fresh) <= len(walked):
                return
            walked, spent = fresh, costs
            if on_frame is not None:
                try:
                    on_frame(list(walked), list(spent))
                except Exception:
                    # A viewer that cannot draw is not a reason to stop
                    # climbing; the geometry at the end is the answer.
                    pass

        try:
            # Straight to a file rather than a pipe.  ORCA writes half a
            # megabyte here, and a pipe nobody is draining fills and stops the
            # process it belongs to -- which would look exactly like a run
            # that hangs, in the one place where hanging is hardest to tell
            # from working.
            with log.open('w', encoding='utf-8') as sink:
                running = subprocess.Popen(
                    [binary, 'in.inp'], cwd=str(folder), env=environment,
                    stdout=sink, stderr=subprocess.STDOUT, text=True)
                ends = None if not timeout else started + float(timeout)
                while running.poll() is None:
                    time.sleep(WATCH_EVERY)
                    _look()
                    if should_stop is not None and should_stop():
                        halted = True
                        _stop(running)
                        break
                    if ends is not None and time.perf_counter() > ends:
                        ran_out = True
                        _stop(running)
                        break
        except OSError as exc:
            return {'ok': False, 'status': f'ORCA did not run: {exc}'}
        _look()
        output = log.read_text(encoding='utf-8', errors='ignore') \
            if log.is_file() else ''
        seconds = time.perf_counter() - started
        reached = folder / 'in.xyz'
        # ORCA writes the optimised geometry back over the file it was given.
        text = reached.read_text(encoding='utf-8', errors='ignore') \
            if reached.is_file() else ''
        rest = {'seconds': seconds, 'path': walked, 'energies': spent,
                'xyz': text or None}
        if ran_out:
            return dict(rest, ok=False, halted=False, output=output[-4000:],
                        status=(f'The transition-state optimisation ran past '
                                f'{float(timeout):.0f} s and was stopped. This '
                                f'runs on the machine the dashboard is on, so '
                                f'it is kept short; '
                                f'submit it as a job for more.'))
        if halted:
            return dict(rest, ok=False, halted=True, output=output[-4000:],
                        status=('Stopped. The climb got '
                                f'{len(walked)} steps in, and the structure '
                                'it had reached is shown.'))
        if not _DONE_RE.search(output):
            return dict(rest, ok=False, halted=False, output=output[-4000:],
                        status=('The transition-state optimisation did not '
                                'converge; the structure it reached is '
                                'shown.'))
        # What was reached, from a Hessian on the geometry that was reached.
        shape = _last_modes(output)
        confirmed = False
        if confirm and text:
            checked = _gfn.optimize_with_gfn(
                text, method, charge=charge, uhf=uhf, solvent=solvent,
                solvation_model='alpb', optimise=False, free_energy=True,
                timeout=max(60.0, float(timeout or 180.0)))
            if checked.get('ok') and checked.get('imaginary') is not None:
                shape = checked['imaginary']
                confirmed = True
                # The Hessian is part of what the press cost, so it is in the
                # time -- but the climb is what ORCA did, and the sentence
                # says that.
                rest['seconds'] = time.perf_counter() - started
        return dict(rest, ok=True, halted=False,
                    imaginary=shape, confirmed=confirmed,
                    verdict=verdict(shape),
                    status=(f'ORCA reached a stationary point on {keyword} in '
                            f'{seconds:.1f} s.'))
    finally:
        shutil.rmtree(folder, ignore_errors=True)


def path_to_saddle(reactant: str, product: str, method: str = 'gfn2', *,
                   charge: int = 0, uhf: int = 0,
                   solvent: Optional[str] = None,
                   solvation_model: str = 'alpb',
                   points: int = _gfn.PATH_POINTS,
                   runs: int = _gfn.PATH_RUNS,
                   path_timeout: Optional[float] = 600.0,
                   max_steps: int = 60,
                   cores: int = 8,
                   timeout: Optional[float] = 180.0,
                   confirm: bool = True,
                   on_frame: Optional[Callable[
                       [Sequence[Sequence[float]],
                        Sequence[Optional[float]]], None]] = None,
                   should_stop: Optional[Callable[[], bool]] = None,
                   on_stage: Optional[Callable[[str], None]] = None,
                   ) -> Dict[str, Any]:
    """Two structures in, a converged saddle out: the path finder into OptTS.

    xtb's path finder walks between two ends and hands back the highest point
    it crossed, which is an estimate and is called one.  ORCA's OptTS sharpens
    an estimate into a stationary point but cannot find one.  Neither half is
    the answer and the pair is, so this is the pair.

    Measured on this machine, the sixteen-atom Diels-Alder, from the two ends
    a scan leaves:

    ==========================  =========  =========================
    route                       wall time  imaginary mode
    ==========================  =========  =========================
    the path finder alone           3.6 s  -79, and an estimate
    this chain                     12.0 s  -393.5
    ``! XTB2 ScanTS``                98 s  -394.1
    ``! XTB2 NEB-TS``               416 s  -393.6
    ==========================  =========  =========================

    Twelve seconds and seven minutes land on the same saddle to within a
    wavenumber.  That is the whole argument for the chain: the expensive
    routes are not buying a better answer here, they are buying the same one.

    A press of its own rather than something a path does on its way past, and
    the reason is in the same table.  The path's own answer -- a barrier
    forward, a barrier back, and how near it came to the structure it was
    aimed at -- is complete in 3.6 s and does not need ORCA at all, which
    matters on a machine that does not have it.  The climb is a different
    question, costs the rest of the time, and can end somewhere nobody asked
    for: measured on a deliberately sloppy pair of ends, it converged in 14 s
    onto a structure with two modes going the wrong way, at -52 and -26 cm-1.
    A press that answers one question should not spend four times as long
    answering a second one and replace the first answer with a structure
    nobody has looked at yet.  So: both halves under one press for whoever
    wants both, and the path alone still there for whoever wants that.

    *on_frame* is handed the climb as it climbs, so the picture keeps up.  The
    walk itself is not streamed: it ends at the product and the climb starts
    from the highest point of it, so playing both down one channel would run
    the picture forward to the end and then jump backwards -- which is the one
    thing a viewer should never do.  *should_stop* is asked between the two
    halves and all the way through the climb; xtb's path finder is a single
    call and cannot be interrupted part way, so a stop pressed during it takes
    effect when it returns.

    Returns the saddle search's answer -- ``{'ok', 'xyz', 'imaginary',
    'verdict', 'path', 'energies', 'halted', 'status'}`` -- with the walk's
    numbers beside it: ``'barrier'``, ``'back'``, ``'reaction'``, ``'rmsd'``,
    ``'estimate'``, ``'walk'``, and the two halves' times as
    ``'path_seconds'`` and ``'saddle_seconds'``.  ``'stage'`` says which half
    it got to.
    """
    started = time.perf_counter()
    walked = _gfn.walk_the_path(reactant, product, method, charge=charge,
                                uhf=uhf, solvent=solvent,
                                solvation_model=solvation_model,
                                points=points, runs=runs, timeout=path_timeout)
    path_seconds = time.perf_counter() - started
    numbers = {'stage': 'path', 'walk': walked, 'path_seconds': path_seconds,
               'seconds': path_seconds,
               'barrier': walked.get('barrier'), 'back': walked.get('back'),
               'reaction': walked.get('reaction'), 'rmsd': walked.get('rmsd'),
               'estimate': walked.get('ts')}
    if not walked.get('ok'):
        return dict(numbers, ok=False, halted=False,
                    status=str(walked.get('status')
                               or 'The path finder did not run.'))
    if not walked.get('ts'):
        return dict(numbers, ok=False, halted=False,
                    status=('The path finder walked a path but left no '
                            'structure for its highest point, so there is '
                            'nothing to climb from.'))
    if should_stop is not None and should_stop():
        return dict(numbers, ok=False, halted=True, xyz=walked.get('ts'),
                    status=(f'Stopped after the path, which took '
                            f'{path_seconds:.1f} s. The structure it '
                            'estimates for the saddle is shown.'))
    if on_stage is not None:
        try:
            on_stage(f'The path took {path_seconds:.1f} s; sharpening its '
                     'estimate into a stationary point...')
        except Exception:
            # A status line that cannot be written is not a reason to stop
            # halfway through a chain.
            pass
    climbed = optimise_to_saddle(
        walked['ts'], method, charge=charge, uhf=uhf, solvent=solvent,
        max_steps=max_steps, cores=cores, timeout=timeout, confirm=confirm,
        on_frame=on_frame, should_stop=should_stop)
    return dict(numbers, **dict(
        climbed, stage='saddle', walk=walked, estimate=walked.get('ts'),
        path_seconds=path_seconds, saddle_seconds=climbed.get('seconds'),
        seconds=time.perf_counter() - started,
        barrier=walked.get('barrier'), back=walked.get('back'),
        reaction=walked.get('reaction'), rmsd=walked.get('rmsd')))
