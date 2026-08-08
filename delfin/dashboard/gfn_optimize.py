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

__all__ = ['GFN_METHODS', 'atom_lines', 'constraint_input', 'find_xtb',
           'find_binary', 'find_gxtb',
           'held_note', 'hold_atoms_at', 'install_command', 'install_script',
           'install_xtb', 'is_gfn_method', 'read_trajectory', 'SOLVENTS',
           'solvent_note',
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
    # g-xTB approximates wB97M-V/def2-TZVPPD and is nearly as quick as GFN2:
    # measured here, a full optimisation is 0.29 s for 32 atoms, 0.68 s for 62
    # and 2.19 s for 92.  It needs a binary of its own -- see find_binary.
    # No solvent: this build takes none.  --alpb, --gbsa and --cpcmx all stop
    # it dead, and --cosmo runs but only writes a COSMO file -- it moved water
    # in water *up* by 1.6 kcal/mol, which is the wrong direction for a
    # solvation energy and means it is not one.
    'gxtb': {'label': 'g-xTB', 'flags': ['--gxtb'], 'max_atoms': 250,
             'reports': None, 'binary': 'gxtb', 'solvation': False},
}

#: The solvents this xtb accepts, asked of the binary rather than taken from a
#: manual: every name here was tried against xtb 6.7.1 and came back
#: parametrised, for GFN2 and for GFN-FF alike.  ALPB is the model -- it is
#: xtb's own recommendation, and the older GBSA knows only a subset of these
#: (no ethanol, no dioxane, no aniline, no ester, no alcohols beyond methanol).
#: Aliases xtb also takes are left out: h2o is water, n-hexane is hexane.
SOLVENTS: Dict[str, str] = {
    '': 'none (gas phase)',
    'acetone': 'acetone', 'acetonitrile': 'acetonitrile',
    'aniline': 'aniline', 'benzaldehyde': 'benzaldehyde',
    'benzene': 'benzene', 'ch2cl2': 'dichloromethane',
    'chcl3': 'chloroform', 'cs2': 'carbon disulfide',
    'dioxane': 'dioxane', 'dmf': 'DMF', 'dmso': 'DMSO',
    'ether': 'diethyl ether', 'ethanol': 'ethanol',
    'ethylacetate': 'ethyl acetate', 'furane': 'furan',
    'hexadecane': 'hexadecane', 'hexane': 'hexane',
    'methanol': 'methanol', 'nitromethane': 'nitromethane',
    'octanol': 'octanol', 'woctanol': 'octanol (wet)',
    'phenol': 'phenol', 'thf': 'THF', 'toluene': 'toluene',
    'water': 'water',
}

#: For callers that cannot be stopped by hand.  The dashboard passes None:
#: Optimise is a switch there, so the person watching decides when a run has
#: gone on long enough -- and a clock that stops a converging optimisation at
#: an arbitrary second is worse than one that never stops it.
DEFAULT_TIMEOUT = 180.0

_ENERGY_RE = re.compile(r'TOTAL ENERGY\s+(-?\d+\.\d+)')
_VERSION_RE = re.compile(r'xtb version\s+([0-9.]+)')
# What the run says it did, taken from its own output rather than from the
# flags we passed: an xtb that ignored a flag would otherwise be indisting-
# uishable from one that honoured it.
_HAMILTONIAN_RE = re.compile(r'Hamiltonian\s+(GFN[0-9A-Za-z-]*)')
_GFNFF_BANNER = 'G F N - F F'
#: What GFN-FF leaves behind that describes the molecule rather than the run.
_GFNFF_TOPOLOGY_FILES = ('gfnff_topo', 'gfnff_charges')


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


def find_gxtb() -> Optional[str]:
    """The xtb that has g-xTB in it, which is not the ordinary one.

    g-xTB is shipped as a statically linked xtb 6.7.1 of its own: the method is
    not in a released tblite yet, and an ordinary xtb **accepts** ``--gxtb``
    and silently runs GFN2 -- measured, the energy came back identical to the
    GFN2 one to the last digit.  So it is looked for under its own name and
    never confused with the xtb beside it.
    """
    override = os.environ.get('DELFIN_GXTB_BINARY')
    if override and Path(override).is_file():
        return override
    try:
        from delfin.qm_runtime import get_qm_tools_bin_dir

        candidate = Path(get_qm_tools_bin_dir()) / 'xtb-gxtb'
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            return str(candidate)
    except Exception:
        pass
    found = shutil.which('xtb-gxtb')
    return found or None


def find_binary(method: Any = None) -> Optional[str]:
    """Whichever program the chosen method needs."""
    spec = GFN_METHODS.get(str(method or '').strip().lower()) or {}
    if spec.get('binary') == 'gxtb':
        return find_gxtb()
    return find_xtb()


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


def install_script() -> Optional[Path]:
    """DELFIN's own installer for the QM tools, if it is where it should be."""
    candidate = (Path(__file__).resolve().parent.parent
                 / 'qm_tools' / 'install_qm_tools.sh')
    return candidate if candidate.is_file() else None


def install_command(tool: str = 'xtb') -> Optional[list]:
    """Exactly what installing *tool* would run, so it can be shown before it is.

    The script takes the tools to install as arguments; naming one keeps it to
    that one rather than fetching crest, dftb+ and the stda bundle behind it.
    """
    script = install_script()
    return ['bash', str(script), str(tool)] if script else None


def install_xtb(
    on_line: Optional[Callable[[str], None]] = None,
    timeout: float = 1800.0,
    tool: str = 'xtb',
) -> Dict[str, Any]:
    """Run DELFIN's installer for *tool* and say where it ended up.

    Never called by itself: fetching a few hundred megabytes through conda is
    not something a program should decide on a user's behalf, and on a cluster
    the answer is often "load the module instead".  The caller asks first.

    *on_line* is handed each line as it is printed, because the install takes
    minutes and a silent one is indistinguishable from a hung one.
    """
    command = install_command(tool)
    if command is None:
        return {'ok': False, 'binary': None, 'lines': [],
                'status': ('DELFIN\'s installer is not next to this copy of '
                           'the dashboard, so there is nothing to run.')}
    started = time.perf_counter()
    lines: list = []
    try:
        running = subprocess.Popen(
            command, cwd=str(Path(command[1]).parent),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            # Nothing to type into: a dashboard has no terminal to answer a
            # prompt on, and an installer waiting for one would look hung.
            stdin=subprocess.DEVNULL,
        )
    except OSError as exc:
        return {'ok': False, 'binary': None, 'lines': [],
                'status': f'The installer could not be started: {exc}'}
    try:
        # Read as it goes.  The output is small, but it is the only sign of
        # life for several minutes -- and an unread pipe would block it.
        for line in running.stdout:
            text = line.rstrip()
            lines.append(text)
            if on_line is not None:
                try:
                    on_line(text)
                except Exception:
                    pass
            if time.perf_counter() - started > timeout:
                running.kill()
                return {'ok': False, 'binary': None, 'lines': lines,
                        'status': (f'The install was still going after '
                                   f'{int(timeout / 60)} minutes and was '
                                   'stopped.')}
        running.wait(timeout=30)
    finally:
        try:
            running.stdout.close()
        except Exception:
            pass
    spent = time.perf_counter() - started
    binary = find_gxtb() if tool == 'gxtb' else find_xtb()
    name = 'g-xTB' if tool == 'gxtb' else 'xtb'
    if running.returncode == 0 and binary:
        return {'ok': True, 'binary': binary, 'lines': lines,
                'status': f'{name} installed at {binary} in {spent / 60:.1f} min.'}
    said = ''
    for text in reversed(lines):
        if 'ERROR' in text or 'error' in text:
            said = text
            break
    return {
        'ok': False, 'binary': binary, 'lines': lines,
        'status': (f'The install ended without a working {name}: {said}'
                   if said else
                   f'The install ended without a {name} the dashboard can '
                   'find.'),
    }


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


#: What the editor's two ways of holding a value mean to xtb, in its own
#: force-constant units.  Measured on a propane whose C1-C2 was asked for
#: 1.700 A and whose C1-C2-C3 angle was asked for 100.0 deg, GFN-FF:
#:
#:     0.25   1.6463 A   103.85 deg
#:     0.5    1.6704 A   102.23 deg
#:     1.0    1.6843 A   101.20 deg
#:     5.0    1.6967 A   100.26 deg
#:     20.0   1.6992 A   100.07 deg
#:
#: A pull is meant to negotiate with the chemistry and settle at a compromise,
#: which is what 0.5 does; a fix is meant to be met, and 20 meets it to a
#: thousandth of an Angstrom -- past the digits the value box shows.
PULL_FORCE_CONSTANT = 0.5
FIX_FORCE_CONSTANT = 20.0

#: How many atoms each kind of held value names.
CONSTRAINT_ATOMS = {'distance': 2, 'angle': 3, 'dihedral': 4}


def constraint_input(
    constraints: Any = (),
    atoms: Optional[int] = None,
) -> Dict[str, Any]:
    """The xtb input file for the values the editor is holding.

    Returns ``{'text', 'force', 'held', 'dropped', 'mixed'}``.  The text is
    empty when there is nothing to hold, and then no input file is written at
    all.

    Three things about xtb decide the shape of this, all of them measured on
    xtb 6.7.1 rather than taken from the manual.

    It numbers atoms from one, not from zero.

    It takes **one** force constant for the whole ``$constrain`` block: a
    second ``force constant=`` line, in the same block or in a second block, is
    read and ignored.  So a set with anything exact in it is held exact
    throughout, and the caller is handed ``mixed`` to say so rather than left
    to wonder why a pull stopped negotiating.

    And it holds *internal* coordinates only.  ``$constrain atoms:`` naming a
    single atom changes nothing at all -- the geometry comes back identical to
    the free one, digit for digit -- and ``$fix atoms:`` is worse than nothing
    in this build: three fixed carbons of a propane came back with their C-C
    at 1.4555 A instead of 1.5255 under GFN-FF, and at 0.4623 A under GFN2.
    Holding an atom where it is, is therefore not something xtb is asked to
    do; see :func:`optimize_with_gfn` for who does it instead.
    """
    kept, dropped = [], []
    exact = negotiated = False
    for entry in (constraints or ()):
        kind = str((entry or {}).get('kind') or '').strip().lower()
        wanted = CONSTRAINT_ATOMS.get(kind)
        indices = [int(i) for i in ((entry or {}).get('atoms') or ())]
        if wanted is None or len(indices) != wanted:
            dropped.append(entry)
            continue
        if any(i < 0 or (atoms is not None and i >= atoms) for i in indices):
            # A held value that names an atom this structure does not have is
            # dropped rather than passed on: xtb would stop with a parse error
            # about a line the user never typed.
            dropped.append(entry)
            continue
        value = (entry or {}).get('value')
        try:
            value = float(value)
        except (TypeError, ValueError):
            dropped.append(entry)
            continue
        mode = str((entry or {}).get('mode') or 'pull').lower()
        if mode == 'fix':
            exact = True
        else:
            negotiated = True
        kept.append((kind, indices, value))

    force = FIX_FORCE_CONSTANT if exact else PULL_FORCE_CONSTANT
    lines = []
    if kept:
        lines.append('$constrain')
        lines.append(f'  force constant={force}')
        for kind, indices, value in kept:
            numbers = ', '.join(str(i + 1) for i in indices)
            lines.append(f'  {kind}: {numbers}, {value:.6f}')
        lines.append('$end')
    return {
        'text': ('\n'.join(lines) + '\n') if lines else '',
        'force': force if kept else None,
        'held': len(kept),
        'dropped': dropped,
        # Both kinds in one set, and only one force constant to hold them with.
        'mixed': bool(exact and negotiated),
    }


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


def solvent_note(solvent: Any) -> str:
    """Which solvent a result is about, or that it is about none.

    A geometry optimised in the gas phase and one optimised in water are two
    different answers to two different questions, and a result that does not
    say which it is invites them to be compared.
    """
    wet = str(solvent or '').strip().lower()
    if not wet or wet not in SOLVENTS or wet == '':
        return ''
    return f' In {SOLVENTS[wet]} (ALPB).'


def held_note(held: Dict[str, Any]) -> str:
    """What was held, said out loud.

    A constraint that was quietly dropped, or a pull that was held as firmly as
    a fix because xtb has only one force constant to give, makes the result an
    answer to a different question than the one that was asked.
    """
    said = []
    if held['held']:
        said.append(f"{held['held']} held value(s) at force constant "
                    f"{held['force']:g}")
    if held['dropped']:
        said.append(f"{len(held['dropped'])} held value(s) dropped -- they "
                    'name atoms this structure does not have')
    if held['mixed']:
        said.append('xtb takes one force constant for the whole set, so the '
                    'pulls were held as firmly as the exact values')
    return (' ' + '; '.join(said) + '.') if said else ''


def optimize_with_gfn(
    xyz_text: str,
    method: str = 'gfnff',
    *,
    charge: int = 0,
    uhf: int = 0,
    max_steps: Optional[int] = None,
    timeout: Optional[float] = DEFAULT_TIMEOUT,
    max_atoms: Optional[int] = None,
    should_stop: Optional[Callable[[], bool]] = None,
    on_frames: Optional[Callable[[list], None]] = None,
    constraints: Any = (),
    topology: Optional[Path] = None,
    solvent: Optional[str] = None,
) -> Dict[str, Any]:
    """Relax *xyz_text* with xtb and say what happened.

    Returns ``{'ok', 'xyz', 'energy', 'status', 'seconds', 'method'}``.  On any
    failure ``ok`` is False, ``xyz`` is the input unchanged and ``status`` says
    why in a sentence a chemist can act on -- a structure that silently comes
    back unrelaxed is worse than one that says it did not converge.

    *constraints* are the values the editor is holding, in its own shape --
    ``{'kind', 'atoms', 'value', 'mode'}`` with atoms counted from zero; see
    :func:`constraint_input` for what each mode becomes.

    *topology* is a directory the caller owns, in which GFN-FF's perceived
    bonding is kept between runs.  GFN-FF works out its topology from whatever
    geometry it is handed, once, and then holds the molecule together with it.
    That is fine for a single run and wrong for a drag: measured on a propane,
    pulling a carbon out to a C-C of 1.87 A relaxes back to 1.49, and at 1.96 A
    the bond is no longer seen at all and the same relaxation pushes it to
    2.80.  A cliff, not a slope -- and a fresh perception on every step of a
    drag means the molecule can fall apart between one answer and the next.
    Given a directory, the perception made at the start is reused: at a C-C of
    2.33 A it still pulls back, to 1.51.

    *solvent* is one of :data:`SOLVENTS`, run with ALPB.  A geometry optimised
    in the gas phase and one optimised in water are different answers, and
    which was asked for belongs in the result rather than in the operator's
    memory -- so it is named in the status.
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

    # A multiplicity the electron count cannot produce is not refused by xtb:
    # water asked for a doublet came back with the singlet's energy to the last
    # digit, -5.070543980552 either way.  A confidently mislabelled answer is
    # the one failure this module exists to prevent, so it is refused here.
    try:
        parity = electron_parity(xyz_text, charge)
    except Exception:
        parity = None                 # cannot tell; not a reason to refuse
    if parity is not None and max(0, int(uhf)) % 2 != parity % 2:
        multiplicity = max(0, int(uhf)) + 1
        return {
            'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
            'seconds': 0.0, 'frames': [],
            'status': (
                f'{atoms} atoms at charge {charge} have an '
                f'{"odd" if parity else "even"} number of electrons, which '
                f'cannot make M = {multiplicity}. xtb does not refuse this -- '
                f'it answers with the nearest multiplicity it can and says '
                f'nothing -- so it is refused here. Try M = '
                f'{multiplicity - 1 if multiplicity > 1 else 2}.'),
        }

    wet = str(solvent or '').strip().lower()
    if wet and spec.get('solvation') is False:
        return {
            'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
            'seconds': 0.0, 'frames': [],
            'status': (f'{label} has no implicit solvation in this build: '
                       'ALPB, GBSA and CPCM-X all stop it, and COSMO only '
                       'writes a file. Optimise it in the gas phase, or '
                       'choose GFN2-xTB or GFN-FF for a solvent.'),
        }
    if wet and wet not in SOLVENTS:
        return {
            'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
            'seconds': 0.0, 'frames': [],
            'status': (f'{wet!r} is not a solvent this xtb is parametrised '
                       f'for. It knows: '
                       + ', '.join(n for n in SOLVENTS if n) + '.'),
        }

    binary = find_binary(key)
    if binary is None:
        wanted = ('its own xtb build, xtb-gxtb -- an ordinary xtb accepts '
                  '--gxtb and silently runs GFN2 instead'
                  if spec.get('binary') == 'gxtb' else 'xtb')
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [],
                'status': (f'{label} needs {wanted}, which was not found in '
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
        if wet:
            command += ['--alpb', wet]
        held = constraint_input(constraints, atoms=len(body))
        if held['text']:
            (folder / 'xtb.inp').write_text(held['text'], encoding='utf-8')
            command += ['--input', 'xtb.inp']
        # GFN-FF's perceived bonding, carried in from the caller's directory
        # if it has one from an earlier run on this molecule.  Only GFN-FF has
        # one; GFN2 works the bonding out afresh from the wavefunction every
        # time and cannot be told.
        keeping = topology is not None and key == 'gfnff'
        if keeping:
            for name in _GFNFF_TOPOLOGY_FILES:
                source_file = Path(topology) / name
                if source_file.is_file():
                    shutil.copy2(str(source_file), str(folder / name))
        # One core, and a scratch directory of its own: two of these must not
        # fight over xtbopt.xyz.
        environment = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
                           OMP_STACKSIZE='1G')
        # xtb's own output goes to a file rather than to a pipe.  A pipe that
        # nobody reads holds 64 KiB and then blocks the program writing to it,
        # and the loop below reads nothing until the process has ended: xtb
        # waits for the loop, the loop waits for xtb, and an optimisation that
        # takes half a second never finishes at all.  Measured on a decane --
        # 77 276 bytes of output from a GFN2 optimisation, which is about
        # thirty cycles, so every real molecule was past it.
        record = folder / 'xtb.out'
        sink = open(record, 'w', encoding='utf-8')
        sent = 0
        try:
            if should_stop is None and on_frames is None:
                subprocess.run(
                    command, cwd=str(folder), stdout=sink, env=environment,
                    stderr=subprocess.STDOUT, timeout=timeout,
                )
            else:
                # Started rather than run, so it can be stopped: an
                # optimisation a user has switched off has to end, not be
                # waited out.
                running = subprocess.Popen(
                    command, cwd=str(folder), stdout=sink,
                    stderr=subprocess.STDOUT, env=environment,
                )
                waited = 0.0
                last_read = 0.0
                last_size = -1
                log = folder / 'xtbopt.log'
                while running.poll() is None:
                    # Stopping comes first and every time round.  Reading the
                    # log is the expensive part -- it is parsed whole and it
                    # grows -- so on a long run it crowded the stop check out
                    # and the switch stopped working.
                    if should_stop is not None and should_stop():
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
                    if timeout and waited > timeout:
                        running.kill()
                        raise subprocess.TimeoutExpired(command, timeout)
                    # xtb writes the log as it optimises, so the path is handed
                    # over while it is still being walked -- but only when the
                    # file has actually grown, and no more than five times a
                    # second.  Parsing it on every pass is what starved the
                    # stop check.
                    if on_frames is not None and waited - last_read >= 0.2:
                        last_read = waited
                        try:
                            size = log.stat().st_size if log.exists() else -1
                        except OSError:
                            size = -1
                        if size != last_size:
                            last_size = size
                            try:
                                walking = read_trajectory(folder)
                            except Exception:
                                walking = []
                            if len(walking) > sent:
                                sent = len(walking)
                                try:
                                    on_frames(walking)
                                except Exception:
                                    pass
                    time.sleep(0.05)
                    waited += 0.05
                running.wait()
        except subprocess.TimeoutExpired:
            return {
                'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': time.perf_counter() - started, 'frames': [],
                'status': (f'{label} did not finish within {int(timeout)} s '
                           'and was stopped.'),
            }
        finally:
            sink.close()
        if keeping:
            # Kept for the next step of the same drag.  The atom count is
            # written beside it: a topology belongs to one molecule, and an
            # atom added or taken away makes it a different one.
            try:
                Path(topology).mkdir(parents=True, exist_ok=True)
                for name in _GFNFF_TOPOLOGY_FILES:
                    made = folder / name
                    if made.is_file():
                        shutil.copy2(str(made), str(Path(topology) / name))
            except OSError:
                pass
        output = record.read_text(encoding='utf-8', errors='replace')
        relaxed = _read_optimised(folder, xyz_text)
        frames = read_trajectory(folder)
        if on_frames is not None and len(frames) > sent:
            # The path, once, at the end.  The loop above reads it five times a
            # second at most, so a run shorter than that interval hands over
            # nothing at all -- and a settle of a small molecule is twenty
            # milliseconds.  The picture then never saw the relaxation it is a
            # picture of, kept the geometry the drag had left, and handed that
            # one to the next drag: the whole path was walked again, in front
            # of the user, as though the earlier drags were being redone.
            try:
                on_frames(frames)
            except Exception:
                pass
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
        if spec.get('reports') is None and reported:
            # g-xTB prints no Hamiltonian line of its own.  An ordinary xtb
            # given --gxtb prints one -- GFN2-xTB -- because it ignored the
            # flag and ran GFN2, and the energy is GFN2's to the last digit.
            # That is the one way this can go wrong silently, so it is the one
            # thing checked.
            return {
                'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': time.perf_counter() - started, 'engine': 'xtb',
                'frames': [], 'version': version, 'hamiltonian': reported,
                'status': (f'{label} was asked for and this xtb ran '
                           f'{reported} instead: it does not have g-xTB in it '
                           'and ignored the flag. Install the g-xTB build.'),
            }
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
                'hamiltonian': reported or wanted or label, 'held': held,
                'converged': False, 'solvent': wet,
                'status': (f'{label} stopped before converging after '
                           f'{seconds:.1f} s; the geometry it reached is shown. '
                           f'(xtb {version}, {reported or wanted or label})'
                           + solvent_note(wet) + held_note(held)),
            }
        return {
            'ok': True, 'xyz': relaxed, 'energy': energy, 'method': key,
            'seconds': seconds, 'engine': 'xtb', 'frames': frames, 'version': version,
            'hamiltonian': reported or wanted or label, 'held': held,
            'converged': True, 'solvent': wet,
            'status': (f'{label} converged in {seconds:.1f} s '
                       f'(xtb {version}, {reported or wanted or label}).'
                       + solvent_note(wet) + held_note(held)),
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
    method: str = 'gfnff',
    charge: int = 0,
    uhf: int = 0,
    cycles: int = 5,
    timeout: float = 30.0,
    constraints: Any = (),
    topology: Optional[Path] = None,
    solvent: Optional[str] = None,
) -> Dict[str, Any]:
    """A few optimisation cycles, for a loop that shows the structure settling.

    Measured on a 102-atom complex with GFN-FF: one cycle 0.06 s, five 0.09 s,
    ten 0.12 s.  So a loop of short runs moves at roughly ten steps a second
    and looks like a relaxation rather than a jump -- which a single 0.4 s call
    to :func:`optimize_with_gfn` does not.

    *method* is the caller's to choose, and it is the method the caller has on
    screen: a loop that quietly ran something else would be a picture of a
    calculation nobody asked for.  It is not free -- GFN2 needs about a second
    for a cycle where GFN-FF needs a twentieth of one, so a loop of GFN2 steps
    is a slideshow and keeps a core busy for it.  Whoever calls this is
    expected to say how long each step took, so the cost is visible rather
    than merely suffered.
    """
    result = optimize_with_gfn(
        xyz_text, method, charge=charge, uhf=uhf,
        max_steps=max(1, int(cycles)), timeout=timeout,
        constraints=constraints, topology=topology, solvent=solvent,
    )
    result['converged'] = bool(
        result.get('ok') and 'converged in' in str(result.get('status') or '')
    )
    return result


def hold_atoms_at(xyz_text: str, reference: str, indices: Any) -> str:
    """*xyz_text*, with the named atoms put back where *reference* has them.

    xtb cannot be asked to leave an atom alone -- it holds internal
    coordinates, not positions -- so an answer about a dragged structure has
    the dragged atom most of the way home again: 244 mA of a 250 mA drag, in
    five cycles.  Putting it back here rather than in the browser matters
    because the answer outlives the drag: applied after the release, when
    nothing is being held any more, it would take the atom with it.
    """
    wanted = sorted({int(i) for i in (indices or ())})
    if not wanted:
        return xyz_text
    rows = [line.split() for line in atom_lines(xyz_text)]
    source = [line.split() for line in atom_lines(reference)]
    if not rows or len(rows) != len(source):
        return xyz_text
    for index in wanted:
        if 0 <= index < len(rows):
            rows[index][1:4] = source[index][1:4]
    body = '\n'.join(' '.join(row) for row in rows)
    head = str(xyz_text or '').splitlines()
    comment = head[1] if len(head) > 1 else ''
    return f'{len(rows)}\n{comment}\n{body}\n'


def largest_shift(before: str, after: str) -> float:
    """How far the atom that moved most has moved, in Angstrom.

    Whether a relaxation is still getting anywhere is a different question from
    whether xtb calls it converged, and it is the one worth asking of a loop
    that runs until it is finished: held values that cannot all be met at once
    never converge, and a relaxation that will not end is a process per round
    for as long as the switch is down.
    """
    first = coordinates_of(before)
    second = coordinates_of(after)
    if not first or len(first) != len(second):
        return 0.0
    worst = 0.0
    for i in range(0, len(first), 3):
        moved = sum((first[i + n] - second[i + n]) ** 2 for n in range(3))
        worst = max(worst, moved)
    return worst ** 0.5


def coordinates_of(xyz_text: str) -> list:
    """The flat [x, y, z, x, y, z, ...] the viewer writes positions from."""
    out: list = []
    for line in atom_lines(xyz_text):
        parts = line.split()
        out.extend((float(parts[1]), float(parts[2]), float(parts[3])))
    return out
