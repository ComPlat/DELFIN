"""Optimise a structure to a transition state, fast enough to press a button.

xtb has no saddle-point optimiser -- 6.7.1 offers ``--hess``, ``--ohess``,
``--bhess`` and ``--path``, and nothing that climbs to a first-order saddle.
ORCA has one, and ORCA can be told to take its gradients from xtb.  That pair
is the whole of this file: ORCA's optimiser, xtb's speed.

Measured on the transition state the path finder estimated for a Diels-Alder,
sixteen atoms: ``! XTB2 OPTTS`` converged in 7 s.  The estimate came in with
its two forming bonds at 2.524 and 2.520 A and one imaginary mode at
-131.4 cm-1, and came out symmetric at 2.315 and 2.315 with the mode at
-243.9.  Seven seconds is a button; an OptTS on a real basis set is a job, and
belongs in the ORCA Builder where jobs are submitted.

So the editor's chain runs end to end without leaving it: a scan walks the
reaction, the path finder finds a way between its two ends and estimates the
saddle, and this sharpens the estimate into a converged transition state --
with the one imaginary frequency that says it is one.

It runs where the dashboard runs, which on a cluster is the login node, and
that is a place with rules.  So: one core, two gigabytes, and a timeout in
minutes rather than hours -- and only the methods that can finish inside
that.  A saddle search on a real basis set is a job, and the ORCA Builder is
where jobs are submitted; this is deliberately the other thing.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Optional

from . import gfn_optimize as _gfn

#: The methods ORCA can drive that are fast enough for this to be a button.
#: Anything with a basis set is a job, not a press.
SADDLE_METHODS = {
    'gfn2': 'XTB2',
    'gfn1': 'XTB1',
    'gfnff': 'XTBFF',
}

#: ORCA says this when the optimiser is finished, and only then.
_DONE_RE = re.compile(r'\*\*\*\s*OPTIMIZATION RUN DONE\s*\*\*\*')
#: The frequencies come from the Hessian OptTS needs anyway, printed by the
#: engine underneath.  The last block is the one that describes what was
#: reached rather than what was started from.
_EIGVAL_RE = re.compile(r'eigval\s*:\s*(.+)')
_IMAGINARY_RE = re.compile(r'#\s*imaginary freq\.\s+(\d+)')


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


def optimise_to_saddle(xyz_text: str, method: str = 'gfn2', *,
                       charge: int = 0, uhf: int = 0,
                       solvent: Optional[str] = None,
                       max_steps: int = 60,
                       cores: int = 1,
                       timeout: Optional[float] = 180.0) -> Dict[str, Any]:
    """Climb to the nearest first-order saddle, and say whether one was found.

    Returns ``{'ok', 'xyz', 'imaginary', 'seconds', 'status'}``.  *imaginary*
    is what the final Hessian says: one mode going the wrong way and no others
    is a transition state, and anything else is not -- which is the whole
    reason to run this rather than trust an estimate.

    A saddle search needs a Hessian to know which mode to climb, so it asks
    for one and recalculates it every five steps: the mode being followed
    changes character as the structure moves, and one taken at the start stops
    describing it.

    One core and two gigabytes by default, and three minutes.  This runs where
    the dashboard runs, and on a cluster that is the login node -- a place
    where a calculation is welcome exactly as long as it is small and short.
    A run that wants more than that is a job, and saying so is better than
    taking the node.
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
    try:
        (folder / 'in.xyz').write_text(
            f'{len(body)}\nfrom the DELFIN viewer\n' + '\n'.join(body) + '\n',
            encoding='utf-8')
        wet = f' ALPB({solvent})' if solvent else ''
        (folder / 'in.inp').write_text(
            f'! {keyword} OPTTS{wet}\n'
            f'%pal\n  nprocs {max(1, int(cores))}\nend\n'
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
        try:
            done = subprocess.run([binary, 'in.inp'], cwd=str(folder),
                                  env=environment, capture_output=True,
                                  text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return {'ok': False,
                    'status': (f'The transition-state optimisation ran past '
                               f'{timeout:.0f} s and was stopped. This runs on '
                               f'the machine the dashboard is on, so it is '
                               f'kept short; submit it as a job for more.')}
        except OSError as exc:
            return {'ok': False, 'status': f'ORCA did not run: {exc}'}
        output = (done.stdout or '') + (done.stderr or '')
        seconds = time.perf_counter() - started
        reached = folder / 'in.xyz'
        # ORCA writes the optimised geometry back over the file it was given.
        text = reached.read_text(encoding='utf-8') if reached.is_file() else ''
        if not _DONE_RE.search(output):
            return {'ok': False, 'seconds': seconds,
                    'xyz': text or None, 'output': output[-4000:],
                    'status': ('The transition-state optimisation did not '
                               'converge; the structure it reached is shown.')}
        return {'ok': True, 'seconds': seconds, 'xyz': text,
                'imaginary': _last_modes(output),
                'status': (f'ORCA reached a stationary point on {keyword} in '
                           f'{seconds:.1f} s.')}
    finally:
        shutil.rmtree(folder, ignore_errors=True)
