"""Optimise a structure to a transition state, fast enough to press a button.

xtb has no saddle-point optimiser -- 6.7.1 offers ``--hess``, ``--ohess``,
``--bhess`` and ``--path``, and nothing that climbs to a first-order saddle.
ORCA has one, and ORCA can be told to take its gradients from xtb.  That pair
is the whole of this file: ORCA's optimiser, xtb's speed.

Three of the four methods here are ORCA's own xtb, named by keyword.  g-xTB is
not: it is a build of its own that no ORCA keyword names, and an ordinary xtb
accepts ``--gxtb`` and silently runs GFN2, so there is nothing to point an
existing keyword at.  It is driven through ``! ExtOpt`` instead -- the
interface ORCA publishes for programs it does not know -- and
:mod:`delfin.dashboard.gxtb_engrad` is what answers it.  Measured on the
sixteen-atom Diels-Alder saddle, that is 101 s, against 97 s for the same
optimiser on GFN2 taken on the same box within the hour -- the box was carrying
a load average of 480 on 384 cores at the time, which is why the GFN2 number is
not the 7 s quoted above; the two are worth comparing with each other and not
with anything else.  A press at that size and a job above it, which is what
:data:`SECONDS_ALLOWED` is for.

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
    # And g-xTB, which ORCA has no keyword for and can still be driven.
    #
    # ORCA runs the three above through its own bundled otool_xtb, an xtb
    # 6.7.1 from 2024; g-xTB is a statically linked xtb 6.7.1 of its own from
    # 2026, and an ordinary xtb accepts --gxtb and silently runs GFN2.  So
    # there is no keyword to add and there never will be one.  What there is
    # is ExtOpt, the interface ORCA publishes for programs it does not know:
    # it writes a geometry and a request, calls whatever EXTOPTEXE names, and
    # reads back an energy and a gradient.  See
    # :mod:`delfin.dashboard.gxtb_engrad`, which is what answers it.
    #
    # It matters because g-xTB approximates wB97M-V/def2-TZVPPD and is the
    # most accurate method the editor offers -- the one somebody chooses when
    # the transition state has to be worth quoting -- and until this it was
    # the one method with no route to a saddle at all.
    'gxtb': 'ExtOpt',
}

#: Which of those are ORCA's own xtb rather than a program of ours.
_DRIVEN_BY_ORCAS_XTB = ('gfn2', 'gfn1', 'gfnff')

#: How often the optimiser computes the Hessian again rather than updating it.
#:
#: ORCA's own ``Recalc_Hess``, and five is what this file has always asked
#: for.  It is not a nicety, and the price of finding out was measured here on
#: the sixteen-atom Diels-Alder estimate under g-xTB, where one numerical
#: Hessian is forty-odd separate processes and dropping the recalculation is
#: tempting: with one Hessian and Bofill updates after it, the optimiser
#: walked the two
#: forming bonds from 2.524 A out to 3.03 and 3.05, ran to the sixty-cycle
#: bound without converging, and left a structure with no imaginary mode at
#: all -- it had gone back down to the van-der-Waals complex.  Which is the
#: same thing :mod:`delfin.dashboard.climb` measured about a carried Hessian
#: on twenty-one drags, and what Baker and Chan named in 1996: an update
#: cannot be relied on to preserve the eigenvalue structure a transition state
#: is defined by.
RECALC_HESS = 5

#: How long a press is given, by method, in seconds.
#:
#: Three minutes for the methods ORCA drives itself, which is what this file
#: has always allowed and what the login node it runs on is owed.
#:
#: g-xTB gets longer because every gradient is a separate process and every
#: numerical Hessian is forty-odd of them with an ORCA start-up each: measured
#: on the sixteen-atom Diels-Alder, ORCA's OptTS converged in 268 s from the
#: path finder's estimate over sixteen cycles and four Hessians, and in 101 s
#: from the saddle GFN2 reaches.  Three minutes would stop the first of those
#: a cycle or two from the answer, which is the worst of both -- the whole
#: cost and none of the result.  Ten minutes is what it needs at that size;
#: above it the run is a job, and the timeout says so and keeps what it had
#: reached.
SECONDS_ALLOWED: Dict[str, float] = {'gxtb': 600.0}
DEFAULT_SECONDS = 180.0


def seconds_for(method: Any) -> float:
    """What a press of this method is allowed, before it is a job instead."""
    return SECONDS_ALLOWED.get(str(method or '').lower(), DEFAULT_SECONDS)


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
    name = _named(order)

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
            # And not called one anywhere else either.  The name goes into
            # the coordinate box, and a box that says "transition state" over
            # a structure the sentence has just refused is where the refusal
            # would be lost.
            name = 'a stationary point with one shallow mode the wrong way'
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
    return {'order': order, 'first_order': first, 'name': name,
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
    key = str(method or '').lower()
    keyword = SADDLE_METHODS.get(key)
    label = (_gfn.GFN_METHODS.get(key) or {}).get('label') or str(method)
    if keyword is None:
        return {'ok': False,
                'status': (f'A saddle search here runs on xtb through ORCA, '
                           f'and {method} is not one of '
                           f'{", ".join(sorted(SADDLE_METHODS))}.')}
    own_binary = None
    if key not in _DRIVEN_BY_ORCAS_XTB:
        own_binary = _gfn.find_binary(key)
        if own_binary is None:
            return {'ok': False,
                    'status': (f'A saddle search on {label} is ORCA driving a '
                               'program of its own, and that program was not '
                               'found. Install it from Settings.')}
        if solvent:
            # Said rather than dropped.  Handed ALPB, this build of g-xTB
            # stops; run without it, the answer would be a gas-phase saddle
            # reported under a solvent nobody removed.
            return {'ok': False,
                    'status': (f'{label} has no implicit solvation in this '
                               'build, so a saddle in a solvent is not '
                               'something it can look for. Search in the gas '
                               'phase, or choose GFN2-xTB.')}
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
        # Whether ORCA is driving one of its own or one of ours.  For g-xTB it
        # is ours, and the hook that answers ORCA's requests is written into
        # this run's own directory so that two runs at once cannot share one.
        own_program = None
        if key not in _DRIVEN_BY_ORCAS_XTB:
            from . import gxtb_engrad

            own_program = gxtb_engrad.write_hook(folder)
        wet = f' ALPB({solvent})' if solvent and own_program is None else ''
        # A numerical Hessian, when the gradient is not ORCA's own to
        # differentiate.  ORCA's default for a method it drives itself is an
        # analytic one, and asked for that with ExtOpt it stops in PROPINT --
        # "ERROR (SHARK): Failed to read input file" -- because there is no
        # basis set for it to work with.  It is ORCA's own remedy: its message
        # elsewhere reads "for optimizations: %geom Calc_Hess true; NumHess
        # true end".
        numerical = '  NumHess true\n' if own_program is not None else ''
        # And ORCA runs on one process when it is not doing the arithmetic.
        #
        # Asked for eight, ORCA runs its numerical Hessian as an eight-process
        # job, and through ExtOpt that is where it comes apart: driven twice
        # here it failed twice and differently -- once with two of forty-six
        # displacements missing and "the Numerical calculation ISN'T
        # COMPLETE", once with "ERROR (ORCA_NUMCALC): Cannot open
        # in.hostnames" -- and both ended the whole search in nine seconds
        # having produced nothing.  There is nothing for those processes to do
        # anyway: with ExtOpt every gradient is somebody else's program, and
        # the cores are better handed to it, which is what
        # DELFIN_GXTB_CORES below does.
        ranks = 1 if own_program is not None else _share(cores)
        (folder / 'in.inp').write_text(
            f'! {keyword} OPTTS{wet}\n'
            f'%pal\n  nprocs {ranks}\nend\n'
            '%maxcore 2000\n'
            '%geom\n'
            '  Calc_Hess true\n'
            + numerical
            + f'  Recalc_Hess {RECALC_HESS}\n'
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
        if own_program is not None:
            # Where ORCA is to send its energy-and-gradient requests, and how
            # many threads whatever answers them may take.  The count is said
            # here rather than left to what ORCA writes into its request,
            # because a thread count is only read when the runtime starts:
            # measured on this box, an unpinned sixteen-atom gradient costs
            # 1.66 s against 6 ms on four threads, and setting the variable
            # afterwards changes nothing at all.
            environment['EXTOPTEXE'] = str(own_program)
            environment['DELFIN_GXTB_CORES'] = str(_share(cores))
            # And which binary, so the hook does not have to look: it is
            # started once per gradient and forty-odd times per Hessian, and
            # importing the resolver to find out costs 80 ms of that each
            # time -- measured, a bare interpreter starts in 37 ms and one
            # that imports gfn_optimize in 119.  This process asked the same
            # resolver a moment ago, before ORCA was started at all.
            environment['DELFIN_GXTB_BINARY'] = str(own_binary)
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
                    # The method, not the keyword.  ExtOpt is how ORCA was
                    # told to ask; g-xTB is what answered, and that is what
                    # the number is about.
                    status=(f'ORCA reached a stationary point on {label} in '
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


#: How many images the band is made of, not counting the two ends.
#:
#: ORCA's own default, and what everything below was measured with.  It is
#: also what says how many cores are worth asking for: the images are
#: independent gradients and ORCA runs them at once, so ``NProcs`` above
#: ``NImages`` has nothing left to hand a process.
NEB_IMAGES = 8

#: How long a nudged elastic band is given, in seconds.
#:
#: Longer than the three minutes :data:`DEFAULT_SECONDS` gives the chain,
#: because this is what is run when the chain's answer is not believed and a
#: band stopped one iteration short has bought nothing at all.  Ten minutes
#: leaves room at several times the size the measurement below was taken at.
NEB_SECONDS = 600.0

#: How many points of the straight line between the two ends are looked at
#: before anything is submitted.
#:
#: Enough to catch a fragment that comes away in the middle and few enough to
#: cost nothing: each point is a covalent-radius bond graph, which is
#: arithmetic on the coordinates and no calculation at all.
NEB_LOOK_AHEAD = 12

#: How long a band is given, by method, in seconds.
#:
#: g-xTB gets longer for the same reason it does under
#: :data:`SECONDS_ALLOWED`: every gradient is a separate process with an ORCA
#: start-up in front of it, and a band is two hundred of them.  Measured on
#: the sixteen-atom Diels-Alder, ``! ExtOpt NEB-TS`` converged in **716 s** --
#: past the ten minutes the methods ORCA drives itself get, which would have
#: stopped it a few iterations from the answer and bought nothing.
NEB_SECONDS_ALLOWED: Dict[str, float] = {'gxtb': 1800.0}


def neb_seconds_for(method: Any) -> float:
    """What a band on this method is allowed, before it is a job instead."""
    return NEB_SECONDS_ALLOWED.get(str(method or '').lower(), NEB_SECONDS)


#: ORCA says this when a band has climbed to a transition state, and only
#: then.  ``*** OPTIMIZATION RUN DONE ***`` is printed by the optimiser inside
#: the band as well -- measured, fourteen times in one converged run -- so
#: :data:`_DONE_RE` does not distinguish the two here.
_NEB_DONE_RE = re.compile(r'THE TS OPTIMIZATION HAS CONVERGED')
#: Where the band's own profile begins.  It is printed twice -- once when the
#: band converges and again with the sharpened saddle in it -- and the second
#: is the one that is about the answer.
_NEB_SUMMARY = 'PATH SUMMARY'
#: One row of it: an image number or ``TS``, and then floats.  How many floats
#: depends on which of the two tables it is (the first carries a distance
#: column the second does not), so the columns are counted from the right --
#: see :func:`_band_profile`.
_NEB_ROW_RE = re.compile(r'^\s*(TS|\d+)\s+((?:\s*-?\d+\.\d+){4,})\s*(?:<=.*)?$',
                         re.MULTILINE)
#: And what it cost, which is the number that says whether a band was worth
#: running rather than what it found.
_NEB_GRADIENTS_RE = re.compile(
    r'Number of SCF / gradient calculations:\s*\n'
    r'\s*NEB\s*\.\.\.\s*(\d+)[^\n]*\n'
    r'\s*TS optimization\s*\.\.\.\s*(\d+)')


def pieces(bonds: Any, atoms: int) -> int:
    """How many separate lumps *atoms* fall into under that bond graph.

    A plain union-find over the pairs :func:`gfn_optimize.bond_graph` hands
    back.  It answers one question -- did this geometry come apart -- and the
    answer has to be cheap enough to ask at every point of an interpolation
    before anything is submitted.
    """
    parent = list(range(max(0, int(atoms))))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for pair in (bonds or ()):
        a, b = int(pair[0]), int(pair[1])
        if 0 <= a < len(parent) and 0 <= b < len(parent):
            ra, rb = root(a), root(b)
            if ra != rb:
                parent[ra] = rb
    return len({root(i) for i in range(len(parent))})


def same_atoms(reactant: str, product: str) -> str:
    """Why these two structures are not two ends of one band, or ``''``.

    The same check :func:`gfn_optimize.walk_the_path` makes, for the same
    reason and in nearly the same words: a band joins atom 1 to atom 1, so
    two structures that are the same molecule in a different order have no
    band between them.  ORCA's own NEB tutorial names hand-prepared endpoints
    in the wrong order as the commonest way of getting nothing back.

    Worth making here even though the editor's usual way of producing a second
    end preserves the order by construction -- morphing one structure into
    another moves atoms and does not renumber them.  ``Adjust H`` is the
    exception: it fills and trims hydrogens around an edit, and a hydrogen
    taken out moves every atom after it up one.  A run that is going to spend
    minutes should find that out before it starts rather than after.

    What it does **not** catch, said plainly rather than left to be
    discovered: two atoms of the same element exchanged.  Measured -- swapping
    two carbons between the two ends of the Diels-Alder leaves the element
    column identical and this returns ``''``.  The band then runs and converges
    onto whatever that mapping describes, which is a real saddle of a
    different reaction.  There is no cheap test for it: any two orderings of
    the same molecule have the same elements in the same places, and telling
    them apart is the atom-mapping problem.  What defends against it is not
    building the second end separately, which is what the message says.
    """
    here = _gfn.atom_lines(reactant)
    there = _gfn.atom_lines(product)
    if not here or not there:
        return 'There are not two structures to put a band between.'
    if len(here) != len(there):
        return ('A band joins atom 1 to atom 1, so the two ends have to be '
                f'the same molecule -- these have {len(here)} and '
                f'{len(there)} atoms. Adjust H is the usual cause: it takes a '
                'hydrogen out of one end, and every atom after it moves up.')
    mine = [line.split()[0] for line in here]
    yours = [line.split()[0] for line in there]
    if mine != yours:
        first = next((n for n, (a, b) in enumerate(zip(mine, yours))
                      if a != b), 0)
        return ('A band joins atom 1 to atom 1, so the two ends have to be '
                f'the same molecule in the same order -- atom {first + 1} is '
                f'{mine[first]} in one and {yours[first]} in the other. Build '
                'the second end from the first rather than separately, and '
                'leave Adjust H off while you do.')
    return ''


def path_comes_apart(reactant: str, product: str,
                     look: int = NEB_LOOK_AHEAD) -> str:
    """Whether the way between the two ends tears a piece off, said in words.

    Measured rather than guessed at.  Given a proton four Angstrom from where
    it is going -- a hydrogen fluoride beside a fluoride at one end and the
    other way round at the other, so the straight line leaves the proton bonded
    to neither -- ``! XTB2 NEB-TS`` printed "Two fragments detected in initial
    structure. Preparing initial structure... Done", then failed three of its
    interior images with "It seems that the Numerical calculation ISN'T
    COMPLETE", and stopped on ``ERROR: GBW-File s_im1.gbw not found``.  Both
    endpoints computed perfectly.  Seventeen seconds there and a longer one on
    a bigger system, and either way nothing back -- and ORCA's own fragment
    handling did not save it, which is why this does not rely on it.

    The endpoints being fine is what makes it hard to see coming, and it is
    exactly what this looks past.

    ORCA does not interpolate in a straight line -- it uses the
    image-dependent pair potential, which exists to avoid this -- so a
    straight line is the pessimistic reading and this is deliberately a
    statement about the two ends rather than a simulation of what ORCA will
    do.  What it catches is the case IDPP cannot rescue either: two ends far
    enough apart that no short way between them keeps the molecule together.

    Returns why, or ``''``.
    """
    here = _gfn.coordinates_of(reactant)
    there = _gfn.coordinates_of(product)
    if not here or len(here) != len(there):
        return ''
    count = len(here) // 3
    rows = [line.split()[0] for line in _gfn.atom_lines(reactant)]
    ends = max(pieces(_gfn.bond_graph(reactant), count),
               pieces(_gfn.bond_graph(product), count))
    worst, worst_at, worst_graph = ends, 0.0, None
    look = max(2, int(look))
    for step in range(1, look):
        share = step / float(look)
        body = []
        for n in range(count):
            body.append('{} {:.8f} {:.8f} {:.8f}'.format(
                rows[n] if n < len(rows) else 'C',
                *[here[3 * n + k] + share * (there[3 * n + k]
                                             - here[3 * n + k])
                  for k in range(3)]))
        graph = _gfn.bond_graph(f'{count}\nbetween\n' + '\n'.join(body) + '\n')
        many = pieces(graph, count)
        if many > worst:
            worst, worst_at, worst_graph = many, share, graph
    if worst <= ends or worst_graph is None:
        return ''
    said = _gfn.graph_changed(_gfn.bond_graph(reactant), worst_graph, rows)
    return ('The way between these two ends comes apart in the middle: '
            f'{int(worst_at * 100)}% of the way across, the structure is in '
            f'{worst} pieces where both ends are in {ends}'
            + (f' ({said})' if said else '')
            + '. A band cannot be run through that -- both ends compute and '
              'every image between them fails -- so bring the two ends '
              'nearer one another first, or walk the reaction in steps.')


def _band_profile(output: str) -> Dict[str, Any]:
    """The barrier the band crossed and what it cost, out of ORCA's summary.

    The rise is read off the ``PATH SUMMARY`` table rather than worked out
    from the images, because the table already carries the transition state as
    a row of its own -- the climbing image after the optimiser has sharpened
    it.  Measured on the Diels-Alder, the highest image says 6.98 kcal/mol and
    the row marked ``TS`` says 6.74, and the second is the barrier.

    From the *last* such table, and its columns counted from the right, both
    for the same reason: ORCA prints the summary twice and the two are not the
    same shape.  The one it writes when the band converges carries a distance
    column that the one after the climb does not, so a fixed column index read
    a distance as an energy -- measured on the g-xTB run, which is where it
    was caught.  ``max(|Fp|)`` and ``RMS(Fp)`` are the last two either way,
    and the rise is the one in front of them.
    """
    text = output or ''
    cut = text.rfind(_NEB_SUMMARY)
    rows = _NEB_ROW_RE.findall(text[cut:] if cut >= 0 else text)
    out: Dict[str, Any] = {'barrier': None, 'reaction': None,
                           'gradients': None}
    read = []
    for label, numbers in rows:
        found = [float(one) for one in numbers.split()]
        if len(found) >= 4:
            read.append((label, found[-3]))
    if read:
        top = next((rise for label, rise in read if label == 'TS'), None)
        out['barrier'] = top if top is not None else max(r for _, r in read)
        out['reaction'] = read[-1][1]
    many = _NEB_GRADIENTS_RE.search(text)
    if many:
        out['gradients'] = int(many.group(1)) + int(many.group(2))
    return out


def neb_to_saddle(reactant: str, product: str, method: str = 'gfn2', *,
                  charge: int = 0, uhf: int = 0,
                  solvent: Optional[str] = None,
                  images: int = NEB_IMAGES,
                  cores: int = 8,
                  timeout: Optional[float] = NEB_SECONDS,
                  confirm: bool = True,
                  on_frame: Optional[Callable[
                      [Sequence[Sequence[float]],
                       Sequence[Optional[float]]], None]] = None,
                  should_stop: Optional[Callable[[], bool]] = None,
                  on_stage: Optional[Callable[[str], None]] = None,
                  ) -> Dict[str, Any]:
    """A nudged elastic band between two ends, and the climb off the top of it.

    ORCA's ``! NEB-TS``: a chain of images relaxed onto the way between two
    structures with energy-weighted springs, the highest of them turned into a
    climbing image, and that one handed to the transition-state optimiser.
    Asgeirsson and Jonsson, JCTC 17, 4929 (2021).

    **It is not the default and must not become one** -- but the reason is not
    the one this file used to give, and the difference was measured rather
    than assumed.

    What was measured, on the sixteen-atom Diels-Alder, from the two ends a
    scan leaves.  A band converged in 18 NEB iterations and 17 climbing steps,
    203 gradients in all, and reached an imaginary mode of **-393.60 cm-1**;
    :func:`path_to_saddle` on the same two ends within the same hour reached
    **-393.53**.  The same saddle to 0.07 cm-1 -- which is the whole argument
    for what this is *for*: a second opinion that shares no machinery with the
    first, worth having exactly when the first one is not believed.

    The wall times are a statement about the box and are given with the load
    they were taken at, on 384 cores:

    ==========================  ==========  =========================
    load average                this band   the chain (path + OptTS)
    ==========================  ==========  =========================
    740                             39.4 s  -
    909                            106.4 s  177.1 + 30.7 s
    912                            156.0 s  245.1 + 34.0 s
    ==========================  ==========  =========================

    So on a many-core machine a band is not the slow route in wall time; it is
    the *parallel* one.  The images are independent gradients and ORCA
    computes them at once, and how much that is worth was measured directly:
    the same band, same box, same hour, **272 s on one process and 39.4 s on
    eight** -- a factor of 6.9.  The 416 s this file used to quote for a band
    is the serial number.  xtb's path finder has no such axis: it is one
    metadynamics walk and it takes what it takes.

    What stays true is the work.  A band is 203 gradients whether or not there
    are cores to spread them over, and on a small login node -- which is where
    a dashboard on a cluster runs -- that is what will be paid.  So: the chain
    first because it is cheaper, and this beside it because it is different.
    ``NProcs`` above ``NImages`` is the one setting that cannot help, because
    there is nothing left to give a ninth process.

    Two things are checked before anything is submitted, because both fail in
    the way that costs the whole timeout and returns nothing:

    * :func:`same_atoms` -- a band joins atom 1 to atom 1.
    * :func:`path_comes_apart` -- two ends whose interpolation tears a
      fragment off compute perfectly at both ends and kill every image
      between them with ``abnormal termination of xtb``.

    It streams.  ORCA writes every band it accepts to ``<base>_MEP_ALL_trj``
    and the climb after it to ``<base>_trj``, both with the energy on the
    comment line, so *on_frame* is handed the band relaxing and then the
    climb, and *should_stop* ends it and keeps what it had reached.  Which
    matters more here than anywhere else in this file: this is the run that
    takes minutes, and a band heading somewhere nobody meant is the thing to
    be able to see it doing.

    Returns the shape :func:`optimise_to_saddle` returns -- ``{'ok', 'xyz',
    'imaginary', 'verdict', 'seconds', 'status', 'path', 'energies'}`` -- with
    ``'barrier'`` and ``'reaction'`` in kcal/mol off the band's own profile
    and ``'gradients'`` for what it cost.
    """
    key = str(method or '').lower()
    keyword = SADDLE_METHODS.get(key)
    label = (_gfn.GFN_METHODS.get(key) or {}).get('label') or str(method)
    if keyword is None:
        return {'ok': False,
                'status': (f'A band here is ORCA driving xtb, and {method} is '
                           f'not one of {", ".join(sorted(SADDLE_METHODS))}.')}
    wrong = same_atoms(reactant, product)
    if wrong:
        return {'ok': False, 'status': wrong}
    apart = path_comes_apart(reactant, product)
    if apart:
        return {'ok': False, 'status': apart}
    own_binary = None
    if key not in _DRIVEN_BY_ORCAS_XTB:
        own_binary = _gfn.find_binary(key)
        if own_binary is None:
            return {'ok': False,
                    'status': (f'A band on {label} is ORCA driving a program '
                               'of its own, and that program was not found. '
                               'Install it from Settings.')}
        if solvent:
            return {'ok': False,
                    'status': (f'{label} has no implicit solvation in this '
                               'build, so a band through a solvent is not '
                               'something it can relax. Run it in the gas '
                               'phase, or choose GFN2-xTB.')}
    binary = find_orca()
    if binary is None:
        return {'ok': False,
                'status': ('A nudged elastic band needs ORCA, which was not '
                           'found. xtb has no band of its own.')}
    body = _gfn.atom_lines(reactant)
    other = _gfn.atom_lines(product)

    started = time.perf_counter()
    folder = Path(tempfile.mkdtemp(prefix='delfin-neb-'))
    walked: List[List[float]] = []
    spent: List[Optional[float]] = []
    try:
        (folder / 'in.xyz').write_text(
            f'{len(body)}\nfrom the DELFIN viewer\n' + '\n'.join(body) + '\n',
            encoding='utf-8')
        (folder / 'to.xyz').write_text(
            f'{len(other)}\nthe other end\n' + '\n'.join(other) + '\n',
            encoding='utf-8')
        own_program = None
        if key not in _DRIVEN_BY_ORCAS_XTB:
            from . import gxtb_engrad

            own_program = gxtb_engrad.write_hook(folder)
        wet = f' ALPB({solvent})' if solvent and own_program is None else ''
        # One process when ORCA is not doing the arithmetic, for the reason
        # :func:`optimise_to_saddle` writes down: through ExtOpt every
        # gradient is a program of ours, and ORCA's own parallel driver
        # overwrites the request files mid-flight.  Otherwise as many as there
        # are images and never more -- past that there is nothing left to
        # compute in parallel.
        band = max(2, int(images))
        ranks = 1 if own_program is not None else min(_share(cores), band)
        # A numerical Hessian for the climb at the end, when the gradient is
        # not ORCA's own to differentiate.  Same remedy, same reason as in
        # :func:`optimise_to_saddle`.
        numerical = ('%geom\n  NumHess true\nend\n'
                     if own_program is not None else '')
        (folder / 'in.inp').write_text(
            f'! {keyword} NEB-TS{wet}\n'
            f'%pal\n  nprocs {ranks}\nend\n'
            '%maxcore 2000\n'
            '%neb\n'
            f'  NImages {band}\n'
            '  NEB_End_XYZFile "to.xyz"\n'
            'end\n'
            + numerical
            + f'* xyzfile {int(charge)} {max(0, int(uhf)) + 1} in.xyz\n',
            encoding='utf-8')
        environment = dict(os.environ)
        environment['PATH'] = (str(Path(binary).parent) + os.pathsep
                               + environment.get('PATH', ''))
        if own_program is not None:
            environment['EXTOPTEXE'] = str(own_program)
            environment['DELFIN_GXTB_CORES'] = str(_share(cores))
            environment['DELFIN_GXTB_BINARY'] = str(own_binary)
        # The band while it relaxes, and then the climb off the top of it.
        # Two files down one channel, in the order ORCA writes them, so what
        # the viewer plays is the whole run rather than half of it.
        trails = (folder / 'in_MEP_ALL_trj.xyz', folder / 'in_trj.xyz')
        log = folder / 'out.log'
        halted = False
        ran_out = False
        stage = ['band']

        def _look() -> None:
            """Read both trails, and hand on whatever is new."""
            nonlocal walked, spent
            fresh: List[List[float]] = []
            costs: List[Optional[float]] = []
            for trail in trails:
                if not trail.is_file():
                    continue
                try:
                    text = trail.read_text(encoding='utf-8', errors='ignore')
                except OSError:
                    continue
                more, paid = frames_in(text, len(body))
                fresh.extend(more)
                costs.extend(paid)
            if len(fresh) <= len(walked):
                return
            if stage[0] == 'band' and trails[1].is_file():
                stage[0] = 'climb'
                if on_stage is not None:
                    try:
                        on_stage('The band is relaxed; climbing from its '
                                 'highest image to the saddle...')
                    except Exception:
                        # A status line that cannot be written is not a
                        # reason to stop a run that is working.
                        pass
            walked, spent = fresh, costs
            if on_frame is not None:
                try:
                    on_frame(list(walked), list(spent))
                except Exception:
                    # A viewer that cannot draw is not a reason to stop; the
                    # geometry at the end is the answer.
                    pass

        try:
            # Straight to a file rather than a pipe, for the reason
            # :func:`optimise_to_saddle` gives: ORCA writes megabytes here and
            # a pipe nobody drains stops the process it belongs to.
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
        reached = folder / 'in_NEB-TS_converged.xyz'
        text = reached.read_text(encoding='utf-8', errors='ignore') \
            if reached.is_file() else ''
        profile = _band_profile(output)
        rest = {'seconds': seconds, 'path': walked, 'energies': spent,
                'xyz': text or None, 'images': band,
                'barrier': profile.get('barrier'),
                'reaction': profile.get('reaction'),
                'gradients': profile.get('gradients')}
        if ran_out:
            return dict(rest, ok=False, halted=False, output=output[-4000:],
                        status=(f'The band ran past {float(timeout):.0f} s '
                                'and was stopped. It runs on the machine the '
                                'dashboard is on, so it is kept short -- and '
                                'the path finder into OptTS reaches the same '
                                'saddle for a fraction of it.'))
        if halted:
            return dict(rest, ok=False, halted=True, output=output[-4000:],
                        status=(f'Stopped. The band got {len(walked)} frames '
                                'in, and what it had reached is shown.'))
        if not _NEB_DONE_RE.search(output):
            # The line that says why, when there is one.  :func:`why_it_stopped`
            # was written for xtb and reads ORCA well -- given the run that
            # died on a stranded fragment it returned "ERROR: GBW-File
            # s_im1.gbw not found." -- but its fallback names xtb, and this is
            # ORCA, so the fallback is not borrowed.
            why = (_gfn.why_it_stopped(output)
                   if 'error' in output.lower() else '')
            return dict(rest, ok=False, halted=False, output=output[-4000:],
                        status=('The band did not reach a transition state'
                                + (f': {why}' if why
                                   else '. What it reached is shown.')))
        # What was reached, from a Hessian on the geometry that was reached --
        # the same argument :func:`optimise_to_saddle` makes, and the same
        # 0.3 s.
        shape = _last_modes(output)
        confirmed = False
        if confirm and text:
            checked = _gfn.optimize_with_gfn(
                text, method, charge=charge, uhf=uhf, solvent=solvent,
                solvation_model='alpb', optimise=False, free_energy=True,
                timeout=max(60.0, float(timeout or NEB_SECONDS)))
            if checked.get('ok') and checked.get('imaginary') is not None:
                shape = checked['imaginary']
                confirmed = True
                rest['seconds'] = time.perf_counter() - started
        rise = profile.get('barrier')
        cost = profile.get('gradients')
        return dict(rest, ok=True, halted=False,
                    imaginary=shape, confirmed=confirmed,
                    verdict=verdict(shape),
                    status=(f'A {band}-image band on {label} reached a '
                            f'stationary point in {seconds:.1f} s'
                            + (f', over a barrier of {rise:.1f} kcal/mol'
                               if rise is not None else '')
                            + (f', and cost {cost} gradients'
                               if cost is not None else '') + '.'))
    finally:
        shutil.rmtree(folder, ignore_errors=True)
