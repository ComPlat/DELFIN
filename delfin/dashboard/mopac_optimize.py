"""Optimise a structure with MOPAC's semiempirical Hamiltonians.

Beside xtb rather than instead of it, and for a measured reason.  On twelve
small organic molecules, comparing the characteristic bond against literature
values, GFN2 came out 11.3 mA off on average and PM6 5.0 -- because GFN2 draws
multiple and polar bonds short: C=C 1.316 against 1.339, C=O 1.192 against
1.208, C-O 1.406 against 1.428, while its plain C-C bonds are right to a
thousandth.  PM6 does not have that lean.

Non-covalent interaction is the other half of the picture and it reverses the
order twice over.  On four dimers against reference binding energies, plain
PM6 is the worst of everything tried -- the water dimer came apart, the two
molecules ending 12.3 A from each other with no interaction left at all -- and
PM7 is little better there (4.86 A, -0.34 kcal/mol against -5.0).  With the
dispersion and hydrogen-bond corrections, PM6-D3H4, the same dimer binds at
2.98 A and -4.59 kcal/mol, and the mean error over the four is 0.58 kcal/mol
where GFN2 gives 0.74.

So **PM6-D3H4** is the one worth having: better bonds than GFN2 and, on this
sample, better intermolecular behaviour.  PM6 and PM7 are offered because a
user may want to compare, and both are labelled with what they cannot do.

Conformers separate nobody: butane's gauche-anti gap is 0.60 (GFN2), 0.72
(PM6, PM7) and 0.82 (PM6-D3H4) against 0.67 kcal/mol.

Transition metals were the case where PM was expected to lose, and on five
closed-shell compounds started at their gas-phase geometries it does not:
mean deviations of 30.7 mA (PM6), 31.4 (PM6-D3H4), 51.8 (GFN2), 56.2 (GFN-FF)
and 60.5 (PM7).  GFN2 leans short here as it does with organics -- TiCl4
2.111 against 2.170, Cr(CO)6 1.881 against 1.918 -- and is 93 mA long on
ferrocene, where PM6 lands within 5.

The average hides the case that matters, though: on Fe(CO)5, PM6 is 103 mA
short where GFN2 is 48.  A method that is better on average and much worse
sometimes is a second opinion, not a replacement, and five neutral
closed-shell carbonyls and halides are not the open-shell charged complexes
most coordination chemistry is made of.  Which is why both are offered and
neither is the default.

Energies here are **heats of formation in kcal/mol**, which is what MOPAC
computes and not the same quantity as xtb's total energy in hartree.  Two
numbers from the two engines cannot be compared, and the unit is carried
through so that nobody tries.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from . import solvents as _solvents

__all__ = ['MOPAC_METHODS', 'find_mopac', 'is_mopac_method', 'mopac_available',
           'optimize_with_mopac', 'read_aux_charges', 'read_aux_frames',
           'frame_as_xyz']

#: What the dropdown offers and the keyword each one means to MOPAC.
#:
#: ``needs`` is what a user has to know before choosing it, in one sentence,
#: and every one of them is something that was measured rather than assumed.
MOPAC_METHODS: Dict[str, Dict[str, Any]] = {
    'pm6d3h4': {
        'label': 'PM6-D3H4', 'keywords': 'PM6-D3H4', 'max_atoms': 400,
        'reports': 'PM6-D3H4',
        'needs': ('PM6 with dispersion and a hydrogen-bond correction. The '
                  'one of these to use: on this machine it bound the water '
                  'dimer at 2.98 A and -4.59 kcal/mol against -5.0.'),
    },
    'pm6': {
        'label': 'PM6', 'keywords': 'PM6', 'max_atoms': 400, 'reports': 'PM6',
        'needs': ('No dispersion and no hydrogen-bond correction: the water '
                  'dimer came apart to 12.3 A. Bond lengths are good; '
                  'anything intermolecular is not.'),
    },
    'pm7': {
        'label': 'PM7', 'keywords': 'PM7', 'max_atoms': 400, 'reports': 'PM7',
        'needs': ('Dispersion is built in, but it made a poor job of the '
                  'water dimer here -- 4.86 A and -0.34 kcal/mol against '
                  '-5.0. Good for heats of formation.'),
    },
}

#: MOPAC names the spin state rather than counting unpaired electrons.
_SPIN_WORDS = {0: 'SINGLET', 1: 'DOUBLET', 2: 'TRIPLET', 3: 'QUARTET',
               4: 'QUINTET', 5: 'SEXTET', 6: 'SEPTET', 7: 'OCTET'}

DEFAULT_TIMEOUT = 300.0

_HEAT_RE = re.compile(r'FINAL HEAT OF FORMATION\s*=\s*(-?\d+\.\d+)')
_VERSION_RE = re.compile(r'MOPAC\s+v?([0-9][0-9.]*)')
#: One block per optimisation cycle in the AUX file, which is what makes a
#: MOPAC run watchable the way an xtb one is.
_FRAME_RE = re.compile(r'ATOM_X_UPDATED:ANGSTROMS\[\s*\d+\]=\s*\n((?:\s*-?\d.*\n)+)')
#: The partial charges, in the same file and written without being asked for.
#:
#: ``AUX`` is already on the keyword line -- it is what the trajectory is read
#: out of -- and MOPAC writes ``ATOM_CHARGES`` into it whether or not anybody
#: reads them.  Measured on a methane under PM7, the block is there after an
#: ordinary run with no extra keyword and no second calculation.  So the
#: editor can put charges on the atoms under PM6, PM6-D3H4 and PM7 for the
#: same nothing they cost under xtb.
#:
#: There is no bond order here.  MOPAC will print one, but only for the
#: ``BONDS`` keyword, and that is a different input -- so under MOPAC the
#: charges are offered and the bond orders are not, which is what the two
#: facts are.
_CHARGES_RE = re.compile(r'ATOM_CHARGES\[\s*\d+\]=\s*\n((?:\s*[-+]?\d.*\n)+)')

WATCH_INTERVAL = 0.01
FRAME_READ_INTERVAL = 0.05


def is_mopac_method(method: Any) -> bool:
    return str(method or '').strip().lower() in MOPAC_METHODS


def find_mopac() -> Optional[str]:
    """Where MOPAC is, asked the way the rest of DELFIN asks."""
    try:
        from delfin.qm_runtime import find_tool_executable

        found = find_tool_executable('mopac')
        if found:
            return str(found)
    except Exception:
        pass
    found = shutil.which('mopac')
    if found:
        return found
    for prefix in (os.sys.prefix, getattr(os.sys, 'base_prefix', os.sys.prefix)):
        candidate = Path(prefix) / 'bin' / 'mopac'
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            return str(candidate)
    return None


def mopac_available() -> bool:
    return find_mopac() is not None


def _atom_rows(xyz_text: str) -> list:
    rows = []
    for line in str(xyz_text or '').splitlines()[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            rows.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
        except ValueError:
            continue
    return rows


def read_aux_frames(path: Path, symbols: list) -> list:
    """Every geometry MOPAC has written so far, as flat coordinate lists.

    The AUX file grows a block per optimisation cycle, so a run can be watched
    while it walks -- the same as xtb's trajectory log, and the reason a MOPAC
    optimisation can be shown moving rather than only at the end.

    **Flat lists, not XYZ text**, because that is what a frame is everywhere
    else here: ``[x1, y1, z1, x2, y2, z2, ...]``, the shape the page's player
    hands straight to setPositions.  Returning text instead put strings where
    numbers belonged, and the viewer drew the molecule in pieces while the
    coordinates in the box were perfectly correct -- which is exactly how it
    looked to the user who reported it.
    """
    try:
        text = path.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return []
    wanted = 3 * len(symbols)
    frames = []
    for match in _FRAME_RE.finditer(text):
        numbers = match.group(1).split()
        if len(numbers) < wanted:
            continue
        try:
            frames.append([float(v) for v in numbers[:wanted]])
        except ValueError:
            continue
    return frames


def read_aux_charges(path: Path, atoms: int) -> Optional[list]:
    """The partial charges MOPAC wrote, or None if it wrote none for *atoms*.

    The last block in the file, because a run that optimised has one per
    energy evaluation and only the final one belongs to the geometry that
    comes back.  A block of the wrong length is refused rather than padded:
    charges that do not line up with the atoms would be drawn on the wrong
    ones, which is worse than not drawing them.
    """
    try:
        text = path.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return None
    values = None
    for match in _CHARGES_RE.finditer(text):
        try:
            found = [float(word) for word in match.group(1).split()]
        except ValueError:
            continue
        if len(found) == int(atoms):
            values = found
    return values


def frame_as_xyz(frame: list, symbols: list, comment: str = 'MOPAC') -> str:
    """One flat frame written out as a coordinate block."""
    lines = [f'{symbol} {frame[3*i]:.8f} {frame[3*i+1]:.8f} {frame[3*i+2]:.8f}'
             for i, symbol in enumerate(symbols)]
    return f'{len(symbols)}\n{comment}\n' + '\n'.join(lines) + '\n'


def _final_geometry(folder: Path, symbols: list) -> Optional[str]:
    """The optimised structure, from the archive MOPAC writes at the end."""
    archive = folder / 'in.arc'
    if not archive.is_file():
        return None
    rows = []
    for line in archive.read_text(encoding='utf-8', errors='replace').splitlines():
        parts = line.split()
        # symbol x flag y flag z flag -- the flags say which coordinates were
        # optimised, and are what tells a coordinate line from prose.
        if len(parts) == 7 and parts[2] in ('+1', '-1', '0', '1'):
            try:
                rows.append(f'{parts[0]} {float(parts[1]):.8f} '
                            f'{float(parts[3]):.8f} {float(parts[5]):.8f}')
            except ValueError:
                continue
    if len(rows) != len(symbols):
        return None
    return f'{len(rows)}\noptimised with MOPAC\n' + '\n'.join(rows) + '\n'


#: How many atoms each kind of held value names -- the editor's own list, so
#: the two engines are asked about the same things.
CONSTRAINT_ATOMS = {'distance': 2, 'angle': 3, 'dihedral': 4}


def freeze_flags(constraints: Any = (), atoms: Optional[int] = None) -> Dict[str, Any]:
    """Which atoms MOPAC is told not to move, and what that does not do.

    MOPAC's Cartesian input carries an optimisation flag per coordinate: 1 is
    optimise, 0 is leave it exactly where it is.  Measured on a propane with
    PM7, C1 and C3 started 3.000 A apart: free, they relax to 2.523 and each
    carbon moves 0.302 A; with their flags at 0 each moves 0.0000 A and the
    distance stays 3.000.  So the flag is a real constraint and not a hint.

    It is not xtb's constraint, and the difference matters twice.

    It fixes the atoms rather than the value between them, which is *more*
    than was asked: a held C-C also stops those two carbons turning and
    translating, and the rest of the molecule then relaxes around a frame that
    is more rigid than the one the user described.

    And it holds the value where it stands rather than moving to a target.
    xtb is handed the number and pulls the geometry to it; MOPAC is handed a
    geometry and told not to move parts of it.  In the editor's own order --
    Set puts the selection at the number, Hold then holds it -- those come to
    the same thing, because by the time it is held the geometry is already at
    the value.  Typing a new number and expecting a PM run to travel to it is
    the case where they part, and the caller is told so rather than left to
    compare a held 1.60 against a geometry that stayed at 1.55.

    A 'pull' cannot be expressed at all: there is one flag and it is on or
    off, so there is no negotiating with the chemistry.  Those are reported as
    not honoured rather than quietly frozen, which would hold as a fix
    something the user asked to have argued with.
    """
    frozen: set = set()
    held, pulls, dropped = 0, 0, []
    for entry in (constraints or ()):
        kind = str((entry or {}).get('kind') or '').strip().lower()
        wanted = CONSTRAINT_ATOMS.get(kind)
        indices = [int(i) for i in ((entry or {}).get('atoms') or ())]
        if wanted is None or len(indices) != wanted:
            dropped.append(entry)
            continue
        if any(i < 0 or (atoms is not None and i >= atoms) for i in indices):
            dropped.append(entry)
            continue
        if str((entry or {}).get('mode') or 'pull').lower() != 'fix':
            pulls += 1
            continue
        frozen.update(indices)
        held += 1
    return {'frozen': frozen, 'held': held, 'pulls': pulls, 'dropped': dropped}


def freeze_note(held: Dict[str, Any]) -> str:
    """What was held, said out loud, in MOPAC's own terms."""
    said = []
    if held['held']:
        said.append(
            f"{held['held']} held value(s) kept by fixing the "
            f"{len(held['frozen'])} atom(s) they name, where they stand -- "
            'MOPAC fixes atoms, not the value between them, so those atoms '
            'also stop turning and moving')
    if held['pulls']:
        said.append(f"{held['pulls']} pull(s) not honoured -- MOPAC's flag is "
                    'on or off, so it has no value to negotiate with; hold '
                    'them as fix, or optimise under a GFN method')
    if held['dropped']:
        said.append(f"{len(held['dropped'])} held value(s) dropped -- they "
                    'name atoms this structure does not have')
    return (' ' + '; '.join(said) + '.') if said else ''


def optimize_with_mopac(
    xyz_text: str,
    method: str = 'pm6d3h4',
    *,
    charge: int = 0,
    uhf: int = 0,
    max_steps: Optional[int] = None,
    timeout: Optional[float] = DEFAULT_TIMEOUT,
    max_atoms: Optional[int] = None,
    should_stop: Optional[Callable[[], bool]] = None,
    on_frames: Optional[Callable[[list], None]] = None,
    solvent: Optional[str] = None,
    constraints: Any = (),
) -> Dict[str, Any]:
    """Relax *xyz_text* with MOPAC and say what happened.

    The same answer shape as the xtb side, so a caller does not have to know
    which engine ran -- except for the energy, which is a heat of formation in
    kcal/mol here and a total energy in hartree there.  ``energy_unit`` says
    which, because the two are not the same quantity and adding a unit is
    cheaper than explaining a comparison nobody should have made.

    *solvent* is one of the names in :mod:`delfin.dashboard.solvents`, and is
    run with COSMO -- the only continuum MOPAC has.  It is switched on by being
    handed the dielectric constant of that solvent, so the PM run and the GFN
    run are asked about the same liquid.  Measured on a glycine in water: PM7
    -10.1 kcal/mol and PM6-D3H4 -9.5 against GFN2 with ALPB at -9.6.  It is not
    free but it is affordable -- 291 ms against 194 for PM7 -- which is what
    makes it usable while a structure is being dragged.
    """
    key = str(method or '').strip().lower()
    spec = MOPAC_METHODS.get(key)
    if spec is None:
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [],
                'status': f'{method!r} is not a MOPAC method.'}

    rows = _atom_rows(xyz_text)
    if not rows:
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [],
                'status': 'There are no coordinates to optimise.'}
    ceiling = max_atoms or spec.get('max_atoms')
    if ceiling and len(rows) > ceiling:
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [],
                'status': (f'{spec["label"]} is offered up to {ceiling} atoms '
                           f'and this has {len(rows)}.')}

    binary = find_mopac()
    if binary is None:
        # Needed and not there is a wait, not a wall -- the same rule the rest
        # of the tools follow. Telling a user that Settings can install it,
        # and then not installing it, is the worst of both.
        told = {}
        try:
            from delfin.qm_health import provide

            told = provide('mopac')
        except Exception:
            told = {}
        binary = find_mopac()
        if binary is None:
            reason = str(told.get('status') or '') if told else ''
            return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                    'seconds': 0.0, 'frames': [],
                    'status': (f'{spec["label"]} needs MOPAC. '
                               + (reason or 'It was not found and could not be '
                                  'installed; Settings shows the installer\'s '
                                  'own output.'))}

    wet = str(solvent or '').strip().lower()
    if wet:
        no = _solvents.refusal('cosmo', wet, key)
        if no:
            return {'ok': False, 'xyz': xyz_text, 'energy': None,
                    'method': key, 'seconds': 0.0, 'frames': [],
                    'status': f'{spec["label"]}: {no}'}

    words = [spec['keywords'], 'PRECISE', 'AUX']
    words += _solvents.mopac_words(wet)
    if int(charge):
        words.append(f'CHARGE={int(charge)}')
    unpaired = max(0, int(uhf))
    if unpaired:
        words.append(_SPIN_WORDS.get(unpaired, 'SINGLET'))
        words.append('UHF')
    if max_steps:
        words.append(f'CYCLES={int(max_steps)}')

    symbols = [row[0] for row in rows]
    started = time.perf_counter()
    folder = Path(tempfile.mkdtemp(prefix='delfin-mopac-'))
    try:
        # An atom a held value names is handed to MOPAC with its flags at
        # zero, which is the only constraint its Cartesian input has.
        holding = freeze_flags(constraints, atoms=len(rows))
        body = '\n'.join(
            '{0} {1:.6f} {4} {2:.6f} {4} {3:.6f} {4}'.format(
                s, x, y, z, 0 if i in holding['frozen'] else 1)
            for i, (s, x, y, z) in enumerate(rows))
        (folder / 'in.mop').write_text(
            ' '.join(words) + '\nfrom the DELFIN viewer\n\n' + body + '\n',
            encoding='utf-8')

        environment = dict(os.environ, OMP_NUM_THREADS='4', MKL_NUM_THREADS='4')
        running = subprocess.Popen(
            [binary, 'in.mop'], cwd=str(folder), env=environment,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            stdin=subprocess.DEVNULL)

        aux = folder / 'in.aux'
        waited, last_read, last_size, sent = 0.0, 0.0, -1, 0
        while running.poll() is None:
            if should_stop is not None and should_stop():
                running.kill()
                break
            if timeout and waited > timeout:
                running.kill()
                return {'ok': False, 'xyz': xyz_text, 'energy': None,
                        'method': key, 'seconds': time.perf_counter() - started,
                        'frames': [],
                        'status': (f'{spec["label"]} was still going after '
                                   f'{timeout:.0f} s and was stopped.')}
            # The same shape as the xtb watcher: asking whether the file grew
            # is free, reading it is not.
            try:
                size = aux.stat().st_size if aux.exists() else -1
            except OSError:
                size = -1
            if (on_frames is not None and size != last_size
                    and waited - last_read >= FRAME_READ_INTERVAL):
                last_read, last_size = waited, size
                walking = read_aux_frames(aux, symbols)
                if len(walking) > sent:
                    sent = len(walking)
                    try:
                        on_frames(walking)
                    except Exception:
                        pass
            time.sleep(WATCH_INTERVAL)
            waited += WATCH_INTERVAL
        running.wait(timeout=30)

        output = ''
        report = folder / 'in.out'
        if report.is_file():
            output = report.read_text(encoding='utf-8', errors='replace')
        frames = read_aux_frames(aux, symbols)
        if on_frames is not None and len(frames) > sent:
            try:
                on_frames(frames)
            except Exception:
                pass

        relaxed = _final_geometry(folder, symbols)
        # Out of the file the frames were just read from, before the folder
        # goes.  Nothing was asked of MOPAC for it -- see :data:`_CHARGES_RE`.
        charges = read_aux_charges(aux, len(symbols))
        found = _HEAT_RE.search(output)
        heat = float(found.group(1)) if found else None
        said = _VERSION_RE.search(output)
        version = said.group(1) if said else ''
        seconds = time.perf_counter() - started

        if relaxed is None and frames:
            # A run stopped at its cycle limit writes no archive -- and the
            # viewer's dynamic mode is nothing but step-limited runs, so
            # calling that a failure would have made every one of them one.
            # The geometry it reached is in the last frame, and it is better
            # than the one that went in even though it is not converged.
            relaxed = frame_as_xyz(frames[-1], symbols,
                                   f'stopped after {len(frames)} cycles')
            return {
                'ok': True, 'xyz': relaxed, 'energy': heat, 'held': holding,
                'energy_unit': 'kcal/mol (heat of formation)',
                'method': key, 'label': spec['label'], 'seconds': seconds,
                'charges': charges, 'bonds': None,
                'engine': 'mopac', 'version': version, 'frames': frames,
                'hamiltonian': spec['reports'], 'converged': False,
                'status': (f'{spec["label"]} stopped after {len(frames)} '
                           'cycles; the geometry it reached is shown.'
                           + freeze_note(holding)),
            }

        if relaxed is None:
            reason = 'MOPAC wrote no optimised geometry'
            for line in output.splitlines():
                stripped = line.strip().strip('*').strip()
                upper = stripped.upper()
                # The error block, not the citation request that stands beside
                # it -- picking the first starred line found "Please cite the
                # open-source release paper", which is true and useless.
                if not stripped or 'CITE' in upper or 'MOPAC' in upper[:6]:
                    continue
                if ('ERROR' in upper or 'EXCESS' in upper or 'FAILED' in upper
                        or 'IMPOSSIBLE' in upper or 'NOT ALLOWED' in upper):
                    reason = stripped
                    break
            return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                    'seconds': seconds, 'frames': [], 'engine': 'mopac',
                    'version': version, 'output': output[-4000:],
                    'status': (f'{spec["label"]}: {reason}. The structure was '
                               'left as it was.')}

        return {
            'ok': True, 'xyz': relaxed, 'energy': heat, 'held': holding,
            'energy_unit': 'kcal/mol (heat of formation)',
            'method': key, 'label': spec['label'], 'seconds': seconds,
            'charges': charges, 'bonds': None,
            'engine': 'mopac', 'version': version, 'frames': frames,
            'hamiltonian': spec['reports'], 'solvent': wet,
            'solvation_model': 'cosmo' if wet else '',
            # The brackets matter: written without them the conditional took
            # the whole concatenation as its first branch, so a MOPAC whose
            # version could not be read reported its result as ".".
            'status': (f'{spec["label"]} finished in {seconds:.1f} s'
                       + (f'; heat of formation {heat:.2f} kcal/mol'
                          if heat is not None else '')
                       + (f'. (MOPAC {version})' if version else '.')
                       + _solvents.note('cosmo', wet)
                       + freeze_note(holding)),
        }
    finally:
        shutil.rmtree(folder, ignore_errors=True)
