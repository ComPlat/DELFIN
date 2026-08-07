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
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Dict, Optional

__all__ = ['GFN_METHODS', 'atom_lines', 'find_xtb', 'is_gfn_method',
           'read_trajectory',
           'xtb_available', 'optimize_with_gfn', 'optimize_autospin',
           'electron_parity']

#: What the dropdown offers, and the flags each one means to xtb.
GFN_METHODS: Dict[str, Dict[str, Any]] = {
    'gfnff': {'label': 'GFN-FF', 'flags': ['--gfnff'], 'max_atoms': 600,
              'reports': 'GFN-FF'},
    'gfn2': {'label': 'GFN2-xTB', 'flags': ['--gfn', '2'], 'max_atoms': 250,
             'reports': 'GFN2-xTB'},
    'gfn1': {'label': 'GFN1-xTB', 'flags': ['--gfn', '1'], 'max_atoms': 250,
             'reports': 'GFN1-xTB'},
}

#: A run that has not finished by then is not going to be useful, and the
#: dashboard is not going to keep a login node busy waiting for it.
DEFAULT_TIMEOUT = 180.0

_ENERGY_RE = re.compile(r'TOTAL ENERGY\s+(-?\d+\.\d+)')
_VERSION_RE = re.compile(r'xtb version\s+([0-9.]+)')
# What the run says it did, taken from its own output rather than from the
# flags we passed: an xtb that ignored a flag would otherwise be indisting-
# uishable from one that honoured it.
_HAMILTONIAN_RE = re.compile(r'Hamiltonian\s+(GFN[0-9A-Za-z-]*)')
_GFNFF_BANNER = 'G F N - F F'


def is_gfn_method(method: Any) -> bool:
    return str(method or '').strip().lower() in GFN_METHODS


def find_xtb() -> Optional[str]:
    """Where xtb is, asked the way the rest of DELFIN asks.

    ``shutil.which`` alone is not enough and was the whole of the problem: the
    kernel a dashboard runs in does not inherit the login shell's PATH, so an
    xtb installed in the very environment the dashboard is running from was
    reported as missing.  DELFIN already has a resolver that knows about
    XTBHOME, its own tool directories and the cluster module paths; it is asked
    first, and the interpreter's own bin directory stands behind it -- an xtb
    beside the python that is running cannot sensibly be called absent.
    """
    try:
        from delfin.qm_runtime import find_tool_executable

        found = find_tool_executable('xtb')
        if found:
            return str(found)
    except Exception:
        pass
    found = shutil.which('xtb')
    if found:
        return found
    for prefix in (sys.prefix, getattr(sys, 'base_prefix', sys.prefix)):
        candidate = Path(prefix) / 'bin' / 'xtb'
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            return str(candidate)
    return None


def _where_it_looked() -> str:
    """The places tried, so a missing xtb is a fixable message, not a wall."""
    places = []
    try:
        from delfin.qm_runtime import get_qm_tools_bin_dir

        places.append(str(get_qm_tools_bin_dir()))
    except Exception:
        pass
    places += ['$DELFIN_XTB_BINARY', '$XTBHOME/$XTBPATH',
               str(Path(sys.executable).resolve().parent), 'PATH']
    return ', '.join(places)


def xtb_available() -> bool:
    return find_xtb() is not None


def atom_lines(xyz_text: str) -> list:
    """The coordinate lines themselves, header or no header.

    The count in the first line is not trusted: xtb reads that many atoms and
    silently ignores the rest, so a block whose header says 89 while carrying
    90 loses one -- and the loss is invisible, because everything else about
    the run succeeds.  The lines are counted instead, and the header written
    to match them.
    """
    raw = [line for line in str(xyz_text or '').splitlines() if line.strip()]
    if not raw:
        return []
    start = 0
    try:
        int(raw[0].split()[0])
        start = 2                     # a header and its comment line
    except (ValueError, IndexError):
        start = 0
    out = []
    for line in raw[start:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        out.append(f'{parts[0]} {parts[1]} {parts[2]} {parts[3]}')
    return out


def _atom_count(xyz_text: str) -> int:
    return len(atom_lines(xyz_text))


def read_trajectory(folder: Path) -> list:
    """Every cycle of the optimisation, as flat coordinate lists.

    xtb writes the whole path to xtbopt.log while it optimises, so the
    trajectory costs nothing extra -- no second run, no loop of processes.
    These are the geometries the optimiser really passed through.
    """
    log = folder / 'xtbopt.log'
    if not log.exists():
        return []
    frames: list = []
    current: list = []
    expected = 0
    for line in log.read_text(encoding='utf-8', errors='replace').splitlines():
        parts = line.split()
        if len(parts) == 1 and parts[0].isdigit():
            if current and len(current) == expected * 3:
                frames.append(current)
            expected = int(parts[0])
            current = []
            skip = True
            continue
        if len(parts) >= 4:
            try:
                current.extend(
                    (float(parts[1]), float(parts[2]), float(parts[3])))
            except ValueError:
                continue
    if current and expected and len(current) == expected * 3:
        frames.append(current)
    return frames


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
    should_stop: Optional[Callable[[], bool]] = None,
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
                'seconds': 0.0, 'frames': [], 'status': f'{method!r} is not a GFN method.'}
    spec = GFN_METHODS[key]
    label = spec['label']

    # What is wrong with the *structure* is worth saying whether or not xtb is
    # installed: too large is too large on any machine, and a caller that has
    # to hear about the missing program first cannot act on either.
    atoms = _atom_count(xyz_text)
    if atoms < 2:
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [], 'status': 'There is nothing to optimise.'}
    ceiling = int(max_atoms if max_atoms is not None else spec['max_atoms'])
    if atoms > ceiling:
        return {
            'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
            'seconds': 0.0, 'frames': [],
            'status': (f'{atoms} atoms is past the {label} limit of {ceiling} '
                       'for an interactive run; submit it as a job instead.'),
        }

    binary = find_xtb()
    if binary is None:
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [],
                'status': (f'{label} needs xtb, which was not found in '
                           f'{_where_it_looked()}.')}

    started = time.perf_counter()
    folder = Path(tempfile.mkdtemp(prefix='delfin-gfn-'))
    try:
        source = folder / 'input.xyz'
        # Written from the lines, with a header that counts them.
        body = atom_lines(xyz_text)
        source.write_text(f'{len(body)}\nfrom the DELFIN viewer\n'
                          + '\n'.join(body) + '\n', encoding='utf-8')
        command = [binary, source.name, *spec['flags'], '--opt',
                   '--chrg', str(int(charge)), '--uhf', str(max(0, int(uhf))),
                   '-P', '1']
        if max_steps:
            command += ['--cycles', str(int(max_steps))]
        # One core, and a scratch directory of its own: two of these must not
        # fight over xtbopt.xyz.
        environment = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
                           OMP_STACKSIZE='1G')
        try:
            if should_stop is None:
                finished = subprocess.run(
                    command, cwd=str(folder), capture_output=True, text=True,
                    timeout=timeout, env=environment,
                )
            else:
                # Started rather than run, so it can be stopped: an optimisation
                # a user has switched off has to end, not be waited out.
                running = subprocess.Popen(
                    command, cwd=str(folder), stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, text=True, env=environment,
                )
                waited = 0.0
                while running.poll() is None:
                    if should_stop():
                        running.terminate()
                        try:
                            running.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            running.kill()
                        # Keep what it reached.  Stopping is how a user
                        # says "that is far enough", not "undo it" -- xtb has
                        # already written the path and the last geometry on it.
                        walked = read_trajectory(folder)
                        # xtbopt.xyz is only written when the optimisation
                        # finishes, so a stopped run has its last geometry in
                        # the log and nowhere else.
                        reached = _read_optimised(folder, xyz_text)
                        if reached is None and walked:
                            symbols = [line.split()[0]
                                       for line in atom_lines(xyz_text)]
                            last = walked[-1]
                            if len(symbols) * 3 == len(last):
                                rows = [
                                    f'{symbols[i]} {last[3*i]:.8f} '
                                    f'{last[3*i+1]:.8f} {last[3*i+2]:.8f}'
                                    for i in range(len(symbols))
                                ]
                                reached = (f'{len(rows)}\nstopped after '
                                           f'{len(walked)} steps\n'
                                           + '\n'.join(rows) + '\n')
                        return {
                            'ok': bool(reached),
                            'xyz': reached or xyz_text, 'energy': None,
                            'method': key, 'frames': walked, 'engine': 'xtb',
                            'seconds': time.perf_counter() - started,
                            'status': (
                                f'{label} stopped after {len(walked)} step(s); '
                                'the geometry it had reached is kept.'
                                if reached else f'{label} was stopped.'),
                        }
                    if waited > timeout:
                        running.kill()
                        raise subprocess.TimeoutExpired(command, timeout)
                    time.sleep(0.05)
                    waited += 0.05
                out, err = running.communicate()
                finished = subprocess.CompletedProcess(
                    command, running.returncode, out, err)
        except subprocess.TimeoutExpired:
            return {
                'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': time.perf_counter() - started, 'frames': [],
                'status': (f'{label} did not finish within {int(timeout)} s '
                           'and was stopped.'),
            }
        output = (finished.stdout or '') + (finished.stderr or '')
        relaxed = _read_optimised(folder, xyz_text)
        frames = read_trajectory(folder)
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
                    'seconds': seconds, 'frames': [],
                    'status': (f'{label}: {reason}. The structure was left as '
                               'it was -- check the charge, the multiplicity '
                               'and whether any atoms overlap.')}

        # Which program, which version, which Hamiltonian -- read out of the
        # run.  Passing --gfn 2 and being given GFN2 are two different claims,
        # and only the second one is evidence.
        version = ''
        found_version = _VERSION_RE.search(output)
        if found_version:
            version = found_version.group(1)
        reported = ''
        found_hamiltonian = _HAMILTONIAN_RE.search(output)
        if found_hamiltonian:
            reported = found_hamiltonian.group(1).strip()
        elif _GFNFF_BANNER in output:
            reported = 'GFN-FF'
        wanted = str(spec.get('reports') or '')
        if reported and wanted and reported.upper() != wanted.upper():
            return {
                'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': time.perf_counter() - started, 'engine': 'xtb',
                'frames': [], 'version': version, 'hamiltonian': reported,
                'status': (f'{label} was asked for and xtb ran {reported}; '
                           'the result is not what it says on the button.'),
            }

        # xtb can leave a geometry behind and still have failed: it writes
        # xtbopt.log as it goes, so the files exist even when it stopped with
        # "something is totally wrong".  Reporting that as a partial success
        # hands over coordinates no one should use.
        if 'abnormal termination' in output or '[ERROR]' in output:
            said = ''
            for line in output.splitlines():
                stripped = line.strip()
                if stripped.startswith('-') and 'xtb' not in stripped[:4]:
                    said = stripped.lstrip('-0123456789 ').strip()
                    if said:
                        break
            return {
                'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': time.perf_counter() - started, 'engine': 'xtb',
                'version': version, 'hamiltonian': reported or wanted,
                'frames': [],
                'status': (f'{label} stopped with an error: '
                           f'{said or "abnormal termination"}. The structure '
                           'was left as it was -- check the charge, the '
                           'multiplicity and whether any atoms overlap.'),
            }

        converged = 'GEOMETRY OPTIMIZATION CONVERGED' in output
        if not converged:
            # The geometry is still better than the one that went in, so it is
            # handed back -- but not as though it were finished.
            return {
                'ok': True, 'xyz': relaxed, 'energy': energy, 'method': key,
                'seconds': seconds, 'engine': 'xtb', 'frames': frames, 'version': version,
                'hamiltonian': reported or wanted,
                'status': (f'{label} stopped before converging after '
                           f'{seconds:.1f} s; the geometry it reached is shown. '
                           f'(xtb {version}, {reported or wanted})'),
            }
        return {
            'ok': True, 'xyz': relaxed, 'energy': energy, 'method': key,
            'seconds': seconds, 'engine': 'xtb', 'frames': frames, 'version': version,
            'hamiltonian': reported or wanted,
            'status': (f'{label} converged in {seconds:.1f} s '
                       f'(xtb {version}, {reported or wanted}).'),
        }
    finally:
        shutil.rmtree(folder, ignore_errors=True)


def electron_parity(xyz_text: str, charge: int = 0) -> int:
    """Whether the molecule has an odd or an even number of electrons.

    An even count can only pair up to a singlet, triplet, quintet ...; an odd
    one to a doublet, quartet ... -- so the parity fixes which multiplicities
    are even possible, and scanning the others would be scanning nonsense.
    """
    from delfin.atom_mapping import _periodic_table

    table = _periodic_table()
    electrons = 0
    for line in str(xyz_text or '').splitlines()[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            electrons += int(table.GetAtomicNumber(parts[0].capitalize()))
        except Exception:
            continue
    return (electrons - int(charge)) % 2


def optimize_autospin(
    xyz_text: str,
    method: str = 'gfnff',
    *,
    charge: int = 0,
    spins: int = 3,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Optimise at several multiplicities and keep the lowest in energy.

    For an open-shell transition metal the spin is not a detail: a fixed guess
    gives a confidently wrong energy and, through it, a wrong geometry.  This
    tries the multiplicities the electron count allows -- the parity and the
    next two above it -- and returns the one that came out lowest, saying which
    it was.  It costs what it sounds like it costs: three runs instead of one.

    Falls back to the parity run if every attempt fails, so the caller always
    gets the same shape of answer.
    """
    parity = electron_parity(xyz_text, charge)
    attempts = []
    best = None
    for step in range(max(1, int(spins))):
        uhf = parity + 2 * step
        stopper = kwargs.get('should_stop')
        if stopper is not None and stopper():
            break                      # switched off between attempts
        result = optimize_with_gfn(xyz_text, method, charge=charge, uhf=uhf,
                                   **kwargs)
        attempts.append((uhf + 1, result.get('energy'), bool(result.get('ok'))))
        if not result.get('ok') or result.get('energy') is None:
            continue
        if best is None or result['energy'] < best['energy']:
            best = dict(result, uhf=uhf, multiplicity=uhf + 1)
    if best is None:
        fallback = optimize_with_gfn(xyz_text, method, charge=charge,
                                     uhf=parity, **kwargs)
        return dict(fallback, uhf=parity, multiplicity=parity + 1,
                    tried=attempts)
    tried = ', '.join(
        f'M={m}: {"failed" if e is None else f"{e:.6f}"}'
        for m, e, _ok in attempts
    )
    best['status'] = (
        f"{best['status']} Lowest of {len(attempts)} multiplicities "
        f"(M={best['multiplicity']}); {tried}."
    )
    best['tried'] = attempts
    return best


def relax_steps(
    xyz_text: str,
    *,
    charge: int = 0,
    uhf: int = 0,
    cycles: int = 5,
    timeout: float = 30.0,
) -> Dict[str, Any]:
    """A few optimisation cycles, for a loop that shows the structure settling.

    Measured on a 102-atom complex: one cycle 0.06 s, five 0.09 s, ten 0.12 s.
    So a loop of short runs moves at roughly ten steps a second and looks like
    a relaxation rather than a jump -- which a single 0.4 s call to
    :func:`optimize_with_gfn` does not.

    GFN-FF only.  GFN2 needs about a second for one cycle, and a loop of those
    is a slideshow that keeps a core busy for it.
    """
    result = optimize_with_gfn(
        xyz_text, 'gfnff', charge=charge, uhf=uhf,
        max_steps=max(1, int(cycles)), timeout=timeout,
    )
    result['converged'] = bool(
        result.get('ok') and 'converged in' in str(result.get('status') or '')
    )
    return result


def coordinates_of(xyz_text: str) -> list:
    """The flat [x, y, z, x, y, z, ...] the viewer writes positions from."""
    out: list = []
    for line in atom_lines(xyz_text):
        parts = line.split()
        out.extend((float(parts[1]), float(parts[2]), float(parts[3])))
    return out
