"""Optimise a structure with xtb's GFN Hamiltonians, from the dashboard.

UFF and MMFF94 run in the browser, which is what makes dragging an atom feel
live: a round trip per frame would cap it at about 13 Hz.  GFN cannot work that
way -- xtb is a compiled program on the server -- so it is offered where a round
trip is the natural shape of the thing: pressing Optimise.

Measured on a 102-atom Pd complex, one core:

    GFN-FF   single point 0.06 s   gradient 0.05 s   optimisation  0.4 s
    GFN2     single point 0.57 s   gradient 0.28 s   optimisation 13.5 s

So GFN-FF is a button press and GFN2 is a short wait, and neither belongs in a
drag.  Both are held to one core and a wall-clock limit: a login node is shared,
and a run that will not converge has to end by itself rather than by being
noticed.

Charge and spin are not optional here the way they are for a force field.  xtb
needs both, and a transition metal with the wrong number of unpaired electrons
gives a confidently wrong answer rather than an error -- so they are taken from
what the tab has on screen and named in the result.
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

__all__ = ['GFN_METHODS', 'is_gfn_method', 'xtb_available', 'optimize_with_gfn']

#: What the dropdown offers, and the flags each one means to xtb.
GFN_METHODS: Dict[str, Dict[str, Any]] = {
    'gfnff': {'label': 'GFN-FF', 'flags': ['--gfnff'], 'max_atoms': 600},
    'gfn2': {'label': 'GFN2-xTB', 'flags': ['--gfn', '2'], 'max_atoms': 250},
    'gfn1': {'label': 'GFN1-xTB', 'flags': ['--gfn', '1'], 'max_atoms': 250},
}

#: A run that has not finished by then is not going to be useful, and the
#: dashboard is not going to keep a login node busy waiting for it.
DEFAULT_TIMEOUT = 180.0

_ENERGY_RE = re.compile(r'TOTAL ENERGY\s+(-?\d+\.\d+)')


def is_gfn_method(method: Any) -> bool:
    return str(method or '').strip().lower() in GFN_METHODS


def xtb_available() -> bool:
    return bool(shutil.which('xtb'))


def _atom_count(xyz_text: str) -> int:
    lines = [line for line in str(xyz_text or '').splitlines() if line.strip()]
    if not lines:
        return 0
    try:
        return int(lines[0].split()[0])
    except (ValueError, IndexError):
        return max(0, len(lines) - 2)


def _read_optimised(folder: Path, fallback: str) -> Optional[str]:
    """xtb writes the relaxed geometry beside the input."""
    for name in ('xtbopt.xyz', 'xtblast.xyz'):
        candidate = folder / name
        if candidate.exists():
            text = candidate.read_text(encoding='utf-8', errors='replace')
            if text.strip():
                return text
    return None if fallback is None else None


def optimize_with_gfn(
    xyz_text: str,
    method: str = 'gfnff',
    *,
    charge: int = 0,
    uhf: int = 0,
    max_steps: Optional[int] = None,
    timeout: float = DEFAULT_TIMEOUT,
    max_atoms: Optional[int] = None,
) -> Dict[str, Any]:
    """Relax *xyz_text* with xtb and say what happened.

    Returns ``{'ok', 'xyz', 'energy', 'status', 'seconds', 'method'}``.  On any
    failure ``ok`` is False, ``xyz`` is the input unchanged and ``status`` says
    why in a sentence a chemist can act on -- a structure that silently comes
    back unrelaxed is worse than one that says it did not converge.
    """
    key = str(method or '').strip().lower()
    if key not in GFN_METHODS:
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'status': f'{method!r} is not a GFN method.'}
    spec = GFN_METHODS[key]
    label = spec['label']

    if not xtb_available():
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0,
                'status': f'{label} needs xtb on the PATH; it was not found.'}

    atoms = _atom_count(xyz_text)
    if atoms < 2:
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'status': 'There is nothing to optimise.'}
    ceiling = int(max_atoms if max_atoms is not None else spec['max_atoms'])
    if atoms > ceiling:
        return {
            'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
            'seconds': 0.0,
            'status': (f'{atoms} atoms is past the {label} limit of {ceiling} '
                       'for an interactive run; submit it as a job instead.'),
        }

    started = time.perf_counter()
    folder = Path(tempfile.mkdtemp(prefix='delfin-gfn-'))
    try:
        source = folder / 'input.xyz'
        source.write_text(xyz_text if xyz_text.endswith('\n') else xyz_text + '\n',
                          encoding='utf-8')
        command = ['xtb', source.name, *spec['flags'], '--opt',
                   '--chrg', str(int(charge)), '--uhf', str(max(0, int(uhf))),
                   '-P', '1']
        if max_steps:
            command += ['--cycles', str(int(max_steps))]
        # One core, and a scratch directory of its own: two of these must not
        # fight over xtbopt.xyz.
        environment = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
                           OMP_STACKSIZE='1G')
        try:
            finished = subprocess.run(
                command, cwd=str(folder), capture_output=True, text=True,
                timeout=timeout, env=environment,
            )
        except subprocess.TimeoutExpired:
            return {
                'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': time.perf_counter() - started,
                'status': (f'{label} did not finish within {int(timeout)} s '
                           'and was stopped.'),
            }
        output = (finished.stdout or '') + (finished.stderr or '')
        relaxed = _read_optimised(folder, xyz_text)
        energy = None
        found = _ENERGY_RE.search(output)
        if found:
            try:
                energy = float(found.group(1))
            except ValueError:
                energy = None
        seconds = time.perf_counter() - started

        if relaxed is None:
            reason = 'xtb wrote no optimised geometry'
            for line in reversed(output.splitlines()):
                if 'error' in line.lower() or 'abnormal' in line.lower():
                    reason = line.strip()
                    break
            return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                    'seconds': seconds, 'status': f'{label}: {reason}.'}

        converged = 'GEOMETRY OPTIMIZATION CONVERGED' in output
        if not converged:
            # The geometry is still better than the one that went in, so it is
            # handed back -- but not as though it were finished.
            return {
                'ok': True, 'xyz': relaxed, 'energy': energy, 'method': key,
                'seconds': seconds,
                'status': (f'{label} stopped before converging after '
                           f'{seconds:.1f} s; the geometry it reached is shown.'),
            }
        return {
            'ok': True, 'xyz': relaxed, 'energy': energy, 'method': key,
            'seconds': seconds,
            'status': f'{label} converged in {seconds:.1f} s.',
        }
    finally:
        shutil.rmtree(folder, ignore_errors=True)
