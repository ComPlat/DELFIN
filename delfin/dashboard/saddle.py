"""Optimise a structure to a transition state, fast enough to press a button.

xtb has no saddle-point optimiser -- 6.7.1 offers ``--hess``, ``--ohess``,
``--bhess`` and ``--path``, and nothing that climbs to a first-order saddle.
ORCA has one, and ORCA can be told to take its gradients from xtb.  That pair
is the whole of this file: ORCA's optimiser, xtb's speed.

Measured on the transition state the path finder estimated for a Diels-Alder,
sixteen atoms: ``! XTB2 OPTTS`` converged in 7 s.  The estimate came in with
its two forming bonds at 2.524 and 2.520 A and one imaginary mode at
-131.4 cm-1, and came out symmetric at 2.315 and 2.315 with the mode at
-372.  Seven seconds is a button; an OptTS on a real basis set is a job, and
belongs in the ORCA Builder where jobs are submitted.

So the editor's chain runs end to end without leaving it: a scan walks the
reaction, the path finder finds a way between its two ends and estimates the
saddle, and this sharpens the estimate into a converged transition state --
with the one imaginary frequency that says it is one.

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
    """The modes of the structure that was reached, from the last Hessian.

    A saddle search recalculates its Hessian as it goes, so the output holds
    several; the last one is about the geometry the run ended on and the
    earlier ones are about geometries it has left.  Read the first, an OptTS
    that improved its structure reports the mode it started with.
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
            'modes': modes[:4]}


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
                       on_frame: Optional[Callable[
                           [Sequence[Sequence[float]],
                            Sequence[Optional[float]]], None]] = None,
                       should_stop: Optional[Callable[[], bool]] = None,
                       ) -> Dict[str, Any]:
    """Climb to the nearest first-order saddle, and say whether one was found.

    Returns ``{'ok', 'xyz', 'imaginary', 'seconds', 'status', 'path',
    'energies'}``.  *imaginary* is what the final Hessian says: one mode going
    the wrong way and no others is a transition state, and anything else is
    not -- which is the whole reason to run this rather than trust an
    estimate.

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
        return dict(rest, ok=True, halted=False,
                    imaginary=_last_modes(output),
                    status=(f'ORCA reached a stationary point on {keyword} in '
                            f'{seconds:.1f} s.'))
    finally:
        shutil.rmtree(folder, ignore_errors=True)
