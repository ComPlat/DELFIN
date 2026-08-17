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

import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from . import solvents as _solvents

__all__ = ['GFN_METHODS', 'atom_lines', 'constraint_input', 'contacts_holding',
           'find_xtb', 'find_binary', 'find_gxtb',
           'held_note', 'hold_atoms_at', 'install_command', 'install_root',
           'install_script',
           'install_xtb', 'is_gfn_method', 'read_trajectory', 'SOLVENTS',
           'solvent_note',
           'xtb_available', 'ensure_binary', 'auto_install_allowed',
           'optimize_with_gfn', 'optimize_autospin',
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
#: parametrised, for GFN2 and for GFN-FF alike.  Aliases xtb also takes are
#: left out: h2o is water, n-hexane is hexane.
#:
#: Which *model* wraps them -- ALPB, GBSA or ddCOSMO -- and which of the three
#: each method can be run with, lives in :mod:`delfin.dashboard.solvents`,
#: because MOPAC has to be told about the same liquids and cannot be told in
#: xtb's words.  This mapping stays for callers that only want the names.
SOLVENTS: Dict[str, str] = dict(
    [('', 'none (gas phase)')]
    + [(name, spec['label']) for name, spec in _solvents.SOLVENTS.items()])

#: For callers that cannot be stopped by hand.  The dashboard passes None:
#: Optimise is a switch there, so the person watching decides when a run has
#: gone on long enough -- and a clock that stops a converging optimisation at
#: an arbitrary second is worse than one that never stops it.
DEFAULT_TIMEOUT = 180.0

#: How often the watching loop looks at a running xtb.  It bounds two things:
#: how quickly a Stop is noticed, and how soon a finished frame can be shown.
#: Ten milliseconds of sleeping costs nothing and was fifty.
WATCH_INTERVAL = 0.01

#: And how often the trajectory may be read when it has grown.  Reading it
#: costs 1.73 ms for a 211 KiB log of 37 frames and the whole file is parsed
#: each time, so the work grows with the run; twenty times a second bounds
#: that at a few percent of one core.  It was five times a second, which held
#: a GFN-FF optimisation to five frames on screen while it was computing a
#: step every fourteen milliseconds.
FRAME_READ_INTERVAL = 0.05

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


#: The oldest xtb that can carry out an optimisation at all.  6.6.1 dies with
#: "Fortran runtime error: Missing comma between descriptors" at
#: optimizer.f90:639 on *every* optimisation -- GFN2 and GFN-FF alike, open
#: shell and closed -- while 6.7.1 finishes the same job normally.  Only a
#: single point without --opt survives on 6.6.1.  Measured on both builds.
XTB_MINIMUM_VERSION = (6, 7, 0)

#: What each binary turned out to be, so the question is asked once a session
#: rather than once a click.  Keyed by path and modification time, so a
#: reinstall at the same path is noticed.
_XTB_JUDGED: Dict[Any, Any] = {}


def _version_tuple(said: str) -> tuple:
    parts = []
    for piece in str(said or '').split('.'):
        digits = ''.join(ch for ch in piece if ch.isdigit())
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def judge_xtb(path: Any) -> Dict[str, Any]:
    """Whether this xtb can do the job, and what it is if it cannot.

    Asking costs one ``--version`` -- about ten milliseconds, once per binary
    per session.  Not asking cost a user a whole afternoon: the binary was
    there, it started, it reported itself, and every optimisation it was given
    died inside its own optimiser.
    """
    where = str(path or '')
    if not where:
        return {'ok': False, 'version': '', 'why': 'no xtb'}
    try:
        stamp = (where, Path(where).stat().st_mtime)
    except OSError:
        return {'ok': False, 'version': '', 'why': 'it is not there'}
    remembered = _XTB_JUDGED.get(stamp)
    if remembered is not None:
        return remembered

    said, why, ok = '', '', False
    try:
        told = subprocess.run([where, '--version'], capture_output=True,
                              text=True, timeout=30)
        text = (told.stdout or '') + (told.stderr or '')
        found = _VERSION_RE.search(text)
        said = found.group(1) if found else ''
        if 'error while loading shared libraries' in text:
            why = ('it cannot start: '
                   + text.split('shared libraries:')[-1].strip().split('\n')[0])
        elif not said:
            why = 'it does not say which version it is'
        elif _version_tuple(said) < XTB_MINIMUM_VERSION:
            wanted = '.'.join(str(part) for part in XTB_MINIMUM_VERSION)
            why = (f'this is xtb {said}, and below {wanted} an optimisation '
                   'dies inside xtb\'s own optimiser -- "Missing comma '
                   'between descriptors" -- whatever method it is given')
        else:
            ok = True
    except (OSError, subprocess.SubprocessError) as problem:
        why = f'it would not run: {problem}'

    verdict = {'ok': ok, 'version': said, 'why': why, 'path': where}
    _XTB_JUDGED[stamp] = verdict
    return verdict


def _xtb_candidates() -> list:
    """Every xtb worth considering, best-placed first."""
    found = []
    try:
        from delfin.qm_runtime import find_tool_executable

        resolved = find_tool_executable('xtb')
        if resolved:
            found.append(str(resolved))
    except Exception:
        pass
    on_path = shutil.which('xtb')
    if on_path:
        found.append(on_path)
    for prefix in (sys.prefix, getattr(sys, 'base_prefix', sys.prefix)):
        candidate = Path(prefix) / 'bin' / 'xtb'
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            found.append(str(candidate))
    ordered = []
    for item in found:
        if item not in ordered:
            ordered.append(item)
    return ordered


def find_xtb() -> Optional[str]:
    """Where xtb is, asked the way the rest of DELFIN asks -- and whether it works.

    ``shutil.which`` alone is not enough and was the first half of the problem:
    the kernel a dashboard runs in does not inherit the login shell's PATH, so
    an xtb installed in the very environment the dashboard is running from was
    reported as missing.  DELFIN already has a resolver that knows about
    XTBHOME, its own tool directories and the cluster module paths; it is asked
    first, and the interpreter's own bin directory stands behind it.

    The second half was taking the first answer on trust.  A tool directory is
    searched before the PATH, so one broken xtb left in it -- by an installer
    that used to solve xtb, crest and dftbplus together and got 6.6.1 -- beat
    the working one standing right behind it, and every optimisation failed
    for a user whose colleague on the same machine had no trouble at all.  So
    each candidate is asked what it is, and one that cannot optimise is
    stepped over rather than used.

    If none can, the best-placed one is still returned: a message about a real
    binary that cannot do the job is worth more than "no xtb found".
    """
    candidates = _xtb_candidates()
    for candidate in candidates:
        if judge_xtb(candidate)['ok']:
            return candidate
    return candidates[0] if candidates else None


def unusable_xtb_note() -> str:
    """Why the xtb that was found first is not the one being used, if so."""
    candidates = _xtb_candidates()
    chosen = find_xtb()
    for candidate in candidates:
        if candidate == chosen:
            return ''
        verdict = judge_xtb(candidate)
        if not verdict['ok']:
            return (f'{candidate} was found first and is not being used: '
                    f'{verdict["why"]}.')
    return ''


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


#: Lines that say a run ended, not why it ended.  Chosen as the reason, every
#: one of them leaves the user with nothing to act on -- which is exactly what
#: happened: an optimisation came back as
#: "GFN2-xTB: Error termination. Backtrace:." and the sentence above the
#: backtrace, the one naming the cause, was thrown away with the rest.
_TERMINATORS = (
    'error termination',
    'backtrace for this error',
    'error stop',
    'abnormal termination',
    'program stopped due to fatal error',
)

#: Where the reason is written, in the order worth looking.  Two error paths
#: run through xtb and they do not look alike: xtb's own handler prints an
#: ``[ERROR]`` block with numbered detail lines, while a Fortran runtime error
#: bypasses it entirely and prints ``At line N of file ...`` followed by
#: ``Fortran runtime error: ...``.  Only the second kind ends in a bare
#: backtrace, and that is the kind whose reason kept getting lost.
_REASON_MARKERS = (
    'fortran runtime error:',
    'program received signal',
)

#: How xtb words a complaint when it is not going through its own error block.
#: Every one of these was seen in a run that failed: a damaged parameter file
#: says "no basis found for atom", a solvent it does not know is "not
#: parametrized", an xyz whose count is wrong is a "missmatch", two atoms on
#: top of each other are a "*very* short distance".
_COMPLAINTS = (
    'no basis found', 'not parametrized', 'missmatch', 'mismatch',
    'short distance', 'cannot ', 'could not ', 'failed', 'is not ',
    'error', 'warning',
)


def why_it_stopped(output: Any) -> str:
    """The sentence in xtb's output that says what went wrong.

    A run that failed leaves thousands of lines and one of them is the reason.
    Picking the last line that contains the word "error" finds a terminator --
    ``Error termination. Backtrace:`` -- every single time a Fortran runtime
    error is what happened, because the numbered backtrace frames below it
    contain no words at all.  The reason sits two lines higher and went into
    the bin with the other thirteen kilobytes.

    Returns a sentence, ending in a full stop, or a plain statement that xtb
    said nothing usable -- never an empty string, and never a terminator.
    """
    lines = [line.rstrip() for line in str(output or '').splitlines()]
    stripped = [line.strip() for line in lines]

    def is_terminator(text: str) -> bool:
        low = text.lower()
        return any(mark in low for mark in _TERMINATORS)

    def is_decoration(text: str) -> bool:
        return not text.strip(' |-.:=*_#')

    # A Fortran runtime error, with the file and line it came from: that pair
    # is the whole diagnosis, so both are kept.
    for index, text in enumerate(stripped):
        low = text.lower()
        for marker in _REASON_MARKERS:
            if marker in low:
                said = text
                if index and stripped[index - 1].lower().startswith('at line '):
                    said = f'{stripped[index - 1]}: {text}'
                return said if said.endswith('.') else said + '.'

    # xtb's own error block: the numbered lines under it are the detail, and
    # all of them are worth having, not the last one.
    detail: list = []
    for index, text in enumerate(stripped):
        if 'program stopped due to fatal error' in text.lower():
            for follow in stripped[index + 1:]:
                if not follow:
                    if detail:
                        break
                    continue
                if is_terminator(follow) or is_decoration(follow):
                    break
                detail.append(follow.lstrip('-0123456789 ').strip())
            break
    if detail:
        said = '; '.join(part for part in detail if part)
        return said if said.endswith('.') else said + '.'

    # Otherwise: xtb's own words, wherever they ended up.
    #
    # Not "the line above the terminator", which is the tempting rule and the
    # wrong one. The two streams are merged into one file and they are not
    # buffered alike: the terminators go to stderr and arrive at once, while
    # everything xtb printed to stdout is flushed when it exits and lands
    # *after* them. Measured on a damaged parameter file -- "ERROR STOP" on
    # line 73 of the captured output and the reason,
    # "no basis found for atom 1 Z= 8", on line 102, the last line of all.
    # Reading upwards from the terminator gave the citation list.
    for text in reversed(stripped):
        low = text.lower()
        if not text or is_terminator(text) or is_decoration(text):
            continue
        if any(mark in low for mark in _COMPLAINTS):
            return text if text.endswith('.') else text + '.'

    return ('xtb stopped without writing a geometry and without saying why. '
            'Check that there is room on the disk for its scratch files, and '
            'that the charge and multiplicity fit the structure.')


def _said_version(output: Any) -> str:
    found = _VERSION_RE.search(str(output or ''))
    return found.group(1) if found else ''


def which_xtb_ran(binary: Any, output: Any = '') -> str:
    """Which program actually ran, named in the failure itself.

    Two accounts on one cluster do not necessarily run the same xtb: a module,
    a conda environment, a build in somebody's home directory. One of them
    optimises and the other stops with a Fortran format error inside xtb's own
    optimiser -- the same structure, the same DELFIN, a different binary. With
    the path and version in the message the two can be compared; without them
    there is nothing to compare.
    """
    said = _said_version(output)
    where = str(binary or '')
    if not where:
        return ''
    return (f'This ran {where}' + (f', xtb {said}' if said else '')
            + ' -- a failure inside xtb itself is a fault of that build, not '
              'of the structure, and another one may not have it.')


def parameter_home(binary: Any) -> Optional[str]:
    """The parameter directory that belongs to this binary, if it has one.

    A conda or system xtb keeps its parameters in ``share/xtb`` beside the
    ``bin`` it lives in.  g-xTB has them compiled in and ships no such
    directory -- running it with an empty home gives the identical energy --
    so there is nothing to point at and nothing is returned.
    """
    try:
        real = Path(str(binary)).resolve()
    except Exception:
        return None
    for share in (real.parent.parent / 'share' / 'xtb',
                  real.parent.parent / 'share'):
        try:
            if (share / 'param_gfn2-xtb.txt').is_file():
                return str(share)
        except OSError:
            continue
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


def install_root() -> Optional[Path]:
    """The tool directory an install should write into.

    The same one the resolver reads, which was not a given: this button used
    to run the packaged copy of the installer in place, while the Settings tab
    installed into the user's own copy -- so the two buttons filled two
    different directories and only one of them was being searched.
    """
    try:
        from delfin.qm_runtime import get_qm_tools_root

        return Path(get_qm_tools_root())
    except Exception:
        return None


def install_script() -> Optional[Path]:
    """DELFIN's own installer for the QM tools, if it is where it should be."""
    candidate = (Path(__file__).resolve().parent.parent
                 / 'qm_tools' / 'install_qm_tools.sh')
    return candidate if candidate.is_file() else None


def install_command(tool: str = 'xtb') -> Optional[list]:
    """Exactly what installing *tool* would run, so it can be shown before it is.

    The script takes the tools to install as arguments; naming one keeps it to
    that one rather than fetching crest, dftb+ and the stda bundle behind it.

    It is pointed at the directory the resolver reads.  Left to itself the
    script installs beside its own file, and this copy of it lives inside the
    package -- so this button filled one directory while the Settings tab
    filled another, and only one of the two was ever searched.  A user pressed
    Install, the tools arrived, and they were reported missing.
    """
    script = install_script()
    if not script:
        return None
    command = ['bash', str(script), str(tool)]
    root = install_root()
    if root is not None and root != script.parent:
        # env(1) rather than a shell assignment, so the command shown is the
        # command run.
        return ['env', f'DELFIN_QM_TOOLS_ROOT={root}',
                f'DELFIN_QM_ROOT={root}'] + command
    return command


#: Tools this session has already tried to put right by itself, so a machine
#: with no network spends the wait once rather than on every press.
_AUTO_TRIED: set = set()


def auto_install_allowed() -> bool:
    """Whether DELFIN may install a missing or broken tool by itself.

    On unless switched off. Fetching a few hundred megabytes is not nothing,
    but neither is the alternative that was measured on a real user: an xtb
    that could not optimise, a message that named no cause, and a day lost --
    with a working build one button press away that nobody knew to press.
    """
    said = str(os.environ.get('DELFIN_AUTO_INSTALL_QM_TOOLS', '1')).strip().lower()
    return said not in ('0', 'no', 'off', 'false')


def ensure_binary(method: Any = None,
                  on_line: Optional[Callable[[str], None]] = None,
                  timeout: float = 1800.0) -> Dict[str, Any]:
    """The program this method needs, installed if it is not there or cannot work.

    Three states are told apart, because they do not have the same answer:
    there is none, there is one that cannot start or cannot optimise, and
    there is one that is fine. Only the first two lead anywhere.

    Once per tool per session: a machine with no network must not spend the
    wait again on every press, and an installer that did not help the first
    time will not help the second.
    """
    tool = 'gxtb' if (GFN_METHODS.get(str(method or '').strip().lower()) or {}
                      ).get('binary') == 'gxtb' else 'xtb'
    found = find_binary(method)
    verdict = judge_xtb(found) if (found and tool == 'xtb') else {
        'ok': bool(found), 'why': '' if found else 'it is not installed'}
    if found and verdict['ok']:
        return {'path': found, 'installed': False, 'ok': True, 'status': ''}

    why = verdict.get('why') or 'it is not installed'
    if not auto_install_allowed():
        return {'path': found, 'installed': False, 'ok': False,
                'status': f'{tool}: {why}. Automatic installation is switched off.'}
    if tool in _AUTO_TRIED:
        return {'path': found, 'installed': False, 'ok': False,
                'status': (f'{tool}: {why}. It was already installed once this '
                           'session and is still not right -- install it from '
                           'Settings, where the installer\'s own output is shown.')}
    if install_script() is None:
        return {'path': found, 'installed': False, 'ok': False,
                'status': f'{tool}: {why}, and DELFIN\'s installer is not here.'}

    _AUTO_TRIED.add(tool)
    if on_line is not None:
        on_line(f'{tool}: {why}. Installing a working one -- a few minutes.')
    outcome = install_xtb(on_line=on_line, timeout=timeout, tool=tool)
    # Ask again from scratch: the answer kept about the old binary says
    # nothing about the one that has just replaced it.
    _XTB_JUDGED.clear()
    now = find_binary(method)
    fixed = judge_xtb(now) if (now and tool == 'xtb') else {'ok': bool(now), 'why': ''}
    if now and fixed['ok']:
        return {'path': now, 'installed': True, 'ok': True,
                'status': f'{tool} was not usable and has been installed: {now}.'}
    return {'path': now, 'installed': True, 'ok': False,
            'status': (f'{tool}: {why}. Installing one did not put it right: '
                       f'{fixed.get("why") or outcome.get("status") or "no reason given"}')}


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
        script = install_script()
        running = subprocess.Popen(
            # Where the script lives, not command[1] -- the command starts
            # with the directory it is being pointed at now, and reading a
            # working directory out of it gave a name that is not a path.
            command, cwd=str(script.parent) if script else None,
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
    raw = str(xyz_text or '').splitlines()
    while raw and not raw[0].strip():
        raw.pop(0)                    # a blank line before the header is none
    if not raw:
        return []
    start = 0
    try:
        int(raw[0].split()[0])
        # A header and its comment line -- and the comment line is allowed to
        # be empty, which is what a named block in the ORCA Builder writes.
        # Dropping blank lines before counting the two swallowed it, and with
        # it the first atom: a water came back from xtb as two hydrogens, and
        # everything about the run said it had succeeded.
        start = 2
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

#: And what a value under the hand is held with, which is neither.
#:
#: A fix rings.  The contact a drag holds is perturbed every frame -- the hand
#: moves, the answer moves the rest, and the next frame asks again -- and a
#: spring that stiff answers each perturbation with a step big enough to
#: overshoot the next one.  Measured on a slow drag across an ethylphenol,
#: counting how often an atom reverses direction between one answer and the
#: next, which is what a shaking molecule is:
#:
#:     k        20     10      5      2      1    0.5
#:     worst  0.131  0.102  0.071  0.044  0.045  0.040 A
#:     turns     94    118     23      7      0      1
#:
#: Most of that was an atom that was not properly held in the first place --
#: one distance lets it slide, and the optimiser was chasing it.  Held by two
#: (see contacts_holding) the same drag turns three times at 20.0, so the
#: stiffness matters far less than it looked.  A dragged *group* still gains
#: from softening, because each of its atoms has only one contact: 20.0 turns
#: three times and moves an atom 0.084 A in a step, 5.0 turns once, 1.0 not
#: at all.
#:
#: It costs almost nothing to soften, because the constraint does not have to
#: be stiffer than the bond it is holding.  A benzene C-H asked to stand at
#: 2.200 A comes back at 2.199 under 20.0 and 2.197 under 5.0, priced +89.1
#: either way.  So this sits at the still end of what is still exact.
#:
#: xtb takes one force constant for a whole block, so during a drag this is
#: also what any value the user is holding gets: thousandths of an Angstrom
#: looser while an atom is under the hand, exact again when it is let go.
DRAG_FORCE_CONSTANT = 5.0

#: How much further than the first a second partner may be, and still be
#: holding the same atom rather than tying it to something across the room.
#: Between the measured cases -- 1.00 for a ring carbon, 1.7 to 1.8 for an
#: SN2's chloride, 2.0 for a hydrogen on a ring -- with room on both sides.
_PIN_SPREAD = 1.4

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
    exact = negotiated = dragging = False
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
        elif mode == 'drag':
            dragging = True
        else:
            negotiated = True
        kept.append((kind, indices, value))

    # A value under the hand sets the block's stiffness whatever else is in
    # it: it is the one being perturbed every frame, so it is the one that
    # decides whether the picture is still.  See DRAG_FORCE_CONSTANT.
    force = (DRAG_FORCE_CONSTANT if dragging
             else FIX_FORCE_CONSTANT if exact else PULL_FORCE_CONSTANT)
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


def _dihedral(where: list, a: int, b: int, c: int, d: int) -> float:
    """The a-b-c-d torsion in degrees, from a list of (x, y, z)."""
    def minus(p, q):
        return (p[0] - q[0], p[1] - q[1], p[2] - q[2])

    def cross(p, q):
        return (p[1] * q[2] - p[2] * q[1],
                p[2] * q[0] - p[0] * q[2],
                p[0] * q[1] - p[1] * q[0])

    def dot(p, q):
        return sum(one * two for one, two in zip(p, q))

    b0, b1, b2 = minus(where[a], where[b]), minus(where[c], where[b]), \
        minus(where[d], where[c])
    length = math.sqrt(dot(b1, b1)) or 1.0
    b1 = tuple(one / length for one in b1)
    v = tuple(one - dot(b0, b1) * two for one, two in zip(b0, b1))
    w = tuple(one - dot(b2, b1) * two for one, two in zip(b2, b1))
    return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))


def _turned_by(now: float, before: float) -> float:
    """How far a torsion moved, in degrees, the short way round."""
    return abs((now - before + 180.0) % 360.0 - 180.0)


#: A turn worth holding, and how little a distance may move for the drag to
#: count as a turn at all.  A hand crosses about a tenth of an Angstrom
#: between one answer and the next, so a drive at something shows up as
#: tenths; a pure rotation changes a bonded distance by thousandths (0.1
#: squared over twice the bond length) and a torsion by degrees.
_TURN_DEGREES = 0.5

#: How far the hand must have changed a distance for that distance to be what
#: it is driving.  A hand crosses about a tenth of an Angstrom between one
#: answer and the next, so a drive at something shows up as tenths; a pure
#: rotation changes a bonded length by thousandths -- a tenth squared over
#: twice the bond -- and everything incidental sits below this.
_DRIVEN_ANGSTROM = 0.02


def _went_around(where, then, inside, outside, symbols) -> bool:
    """Whether the hand's step lay across its bond rather than along it."""
    moving = max(inside, key=lambda i: math.dist(where[i], then[i]))
    bonded = [j for j in outside if _is_a_bond(where, symbols, moving, j)]
    if not bonded:
        return False
    j = min(bonded, key=lambda n: math.dist(where[moving], where[n]))
    step = tuple(where[moving][n] - then[moving][n] for n in range(3))
    along = tuple(where[moving][n] - where[j][n] for n in range(3))
    span = math.sqrt(sum(one * one for one in along)) or 1.0
    along = tuple(one / span for one in along)
    radial = sum(one * two for one, two in zip(step, along))
    across = math.sqrt(max(0.0, sum(one * one for one in step) - radial * radial))
    return across > abs(radial)


def _is_a_bond(where, symbols, i, j) -> bool:
    """Whether these two are bonded, by covalent radii."""
    from delfin.atom_mapping import cov_radius
    return (math.dist(where[i], where[j])
            < 1.25 * (cov_radius(symbols[i]) + cov_radius(symbols[j])))


def _snapshot(was, count):
    """*was* as a list of points, or None if it does not fit this structure."""
    if was is None:
        return None
    before = coordinates_of(was) if isinstance(was, str) else list(was or ())
    if len(before) != 3 * count:
        return None
    return [(before[3 * n], before[3 * n + 1], before[3 * n + 2])
            for n in range(count)]


def contacts_holding(
    xyz_text: str,
    dragged: Any,
    most: int = 2,
    was: Any = None,
    turning: Any = None,
) -> list:
    """What a drag has changed, as values to hold while the rest relaxes.

    Returns constraint entries in the editor's own shape, so they go to xtb
    through :func:`constraint_input` like any other held value.

    A dragged geometry cannot be priced as it stands.  A transition state
    *consists* of everything rearranging -- the diene ends pyramidalise, bond
    lengths move -- and held rigid the approach is nothing but repulsion.
    Measured on butadiene and ethylene with GFN2, the two forming bonds at
    2.18 A: the geometry the mouse leaves costs +124.3 kcal/mol, the same
    separation along the relaxed path +3.6.  A budget fed the first number
    needs 1498 K to allow a Diels-Alder that room temperature really does
    allow, which is how this was found.

    Relaxing freely is no answer either: it undoes the drag.  Nor is holding
    the atoms where they are -- xtb cannot (``$fix atoms:`` is broken in 6.7.1
    and ``$constrain atoms:`` on one atom does nothing), and pinning positions
    would not help if it could.  With ``engine=inertial`` xtb does hold the
    atom, to 0.000 A, and hands back the energy of *intact benzene*: the ring
    simply walks the 1.5 A over and takes its hydrogen back.  One atom fixed in
    space holds no bond stretched.

    What is left is what a relaxed scan does: hold the internal coordinates the
    hand moved, free everything else.  Which coordinates matters as much as
    holding them.  Twelve contacts pin the two fragments to each other and
    price the same geometry at +119 -- nothing gained.  It has to be the few
    that carry the reaction:

        contacts held    1      2      3      4
        kcal/mol      +8.3  +19.8  +41.4  +75.2      (true: +3.6)

    So one per dragged atom, distinct partners, and heavy atoms first -- an
    H...H brush past is not a reaction coordinate and holding it stops the
    hydrogens bending out of the way.  For a single dragged atom that is its
    own neighbour, and a hydrogen pulled off a benzene is then priced at +89
    where it belongs, relaxed or not: there is nowhere for the ring to go.
    Breaking was always priced correctly; only forming was not.
    """
    held = {int(i) for i in (dragged or ())}
    rows = [line.split() for line in atom_lines(xyz_text)]
    # A hand that turns does not stop turning between two frames.  Decided
    # afresh every answer, the verdict flipped back and forth on a cyclohexane
    # -- turning, pulling, turning -- and the ring fell back towards the chair
    # each time the constraint set changed under it.  Once a drag is a turn it
    # stays one, and only the angle is measured again.
    if turning is not None and len(turning) == 4:
        spot = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
        if all(0 <= int(n) < len(spot) for n in turning):
            return [{'kind': 'dihedral', 'atoms': [int(n) for n in turning],
                     'mode': 'drag',
                     'value': _dihedral(spot, *[int(n) for n in turning])}]

    inside = sorted(i for i in held if 0 <= i < len(rows))
    outside = [j for j in range(len(rows)) if j not in held]
    if not inside or not outside:
        # The hand is on everything, or on nothing.  Either way no contact
        # between the two describes the drag, and a caller that gets no values
        # back is being told to price the geometry some other way.
        return []
    symbols = [str(r[0]) for r in rows]
    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    light = [one.strip()[:1].upper() == 'H' and len(one.strip()) == 1
             for one in symbols]
    then = _snapshot(was, len(rows))

    def moved(i, j):
        """How far the hand changed this distance, in Angstrom."""
        if then is None:
            return 0.0
        return abs(math.dist(where[i], where[j]) - math.dist(then[i], then[j]))

    # Three kinds of thing a hand can be doing, and they are asked for in the
    # order that tells them apart.  This is the whole of the rule.
    #
    #   1. Making or breaking a bond that is not there yet -- two heavy atoms
    #      being driven at each other.  That is a reaction coordinate whatever
    #      the distance happens to be, so it is looked for by *change* and not
    #      by nearness.  Chosen by nearness instead, one end of a hexadiene
    #      pushed at the other for a Cope shift picks the carbon it is already
    #      bonded to -- whose length the push does not change -- and the two
    #      ends finish 4.19 A apart with the budget reporting +2.0 and no
    #      objection.
    #
    #   2. Stretching a bond that is there.  Pulling a hydrogen off a ring,
    #      tearing a carbon out of one.
    #
    #   3. Neither, and then the hand went around something rather than along
    #      it: a torsion.  Turning a methyl, flipping a chair.
    #
    # Only heavy atoms open a new contact.  A hydrogen sweeping past during a
    # ring flip changes its distance to something by the whole of the step,
    # every step, and held it would stop the hydrogens bending out of the way
    # -- which is most of what makes room for anything.
    driven = sorted(
        ((-moved(i, j), math.dist(where[i], where[j]), i, j)
         for i in inside for j in outside
         if not light[i] and not light[j]
         and moved(i, j) >= _DRIVEN_ANGSTROM
         and not _is_a_bond(where, symbols, i, j)),
        )
    pulled = sorted(
        ((1 if (light[i] or light[j]) else 0),
         math.dist(where[i], where[j]), i, j)
        for i in inside for j in outside)
    pulled = [(0.0, span, i, j) for _rank, span, i, j in pulled]
    if then is not None and not driven and _went_around(where, then, inside,
                                                        outside, symbols):
        # Not at anything and not along anything: the hand went *around*.
        #
        # Which of the two it is cannot be read off a distance, because a
        # hand that turns bends the bond a little as well -- it drags in the
        # plane of the screen, not along an arc.  It can be read off the
        # step: how much of it lay along the bond being held, and how much
        # across it.  That needs no threshold in units that do not compare,
        # and it is what tells a methyl being swung from a hydrogen being
        # pulled off.
        pulled = []
    if then is None:
        # The first answer of a drag has nothing to compare against, so the
        # nearest contacts stand in -- heavy atoms first, as ever.
        driven = []
        pulled = sorted(
            ((1 if (light[i] or light[j]) else 0),
             math.dist(where[i], where[j]), i, j)
            for i in inside for j in outside)
        pulled = [(0.0, span, i, j) for _rank, span, i, j in pulled]

    pairs = driven or pulled
    kept: list = []
    taken_in: set = set()
    taken_out: set = set()
    # One contact per dragged atom -- except when there is only one, which
    # one distance does not hold at all.  A ring carbon pulled out along a
    # bond slides along the other one and the ring closes behind it: asked to
    # stand at 1.72, 1.95 and 2.20 A it came back at 1.41 every time, however
    # stiffly it was held, and the price was that of a squeezed but whole
    # benzene.  Two distances pin it, and the same three pulls then cost
    # +36.6, +75.0 and +105.7 kcal/mol, which is what tearing a ring open is.
    alone = len(inside) == 1
    first = None
    for _change, span, i, j in pairs:
        if (i in taken_in and not alone) or j in taken_out:
            continue
        # And only a partner it is really held *between*.  A second contact
        # much further off is not a hold, it is a tether: in an SN2 the two
        # nearest heavy atoms to the arriving chloride are the carbon it is
        # attacking and the bromide on the far side of it, and holding both
        # nails the leaving group down.  Measured on Cl- + CH3Br, chloride
        # driven in from 2.5 A: with one contact C-Br opens 2.01 -> 3.24 A and
        # the price stays near zero, which is the reaction; with two it stays
        # at 1.94 -> 1.99 and the price runs to +175, which is a bromide being
        # crushed into a carbon that already has four bonds.
        #
        # The ratio tells them apart cleanly.  A ring carbon sits between two
        # partners at 1.41 and 1.41 -- a ratio of one, and it needs both.  The
        # chloride's are 2.5 and 4.4, and a dragged hydrogen's own carbon and
        # the next ring carbon are 1.09 and 2.15: both far outside, and both
        # are held perfectly well by one.
        if first is not None and span > first * _PIN_SPREAD:
            break
        first = span if first is None else first
        kept.append({'kind': 'distance', 'atoms': [i, j],
                     'value': span, 'mode': 'drag'})
        taken_in.add(i)
        taken_out.add(j)
        if len(kept) >= max(1, int(most)):
            break

    if kept or then is None:
        return kept

    # Neither made nor stretched: the hand went around something.
    #
    # A distance cannot hold a turn -- rotating about a bond leaves that
    # bond's length exactly where it was -- so with only distances the budget
    # had no grip on a torsion at all, and neither did the drag.  Dragging the
    # carbon of a methyl round a chain bond, holding only its C-C: the
    # hydrogens were left behind and the C-H distances went to 0.92 and 2.12
    # A.  Holding the torsion instead, the same swing keeps C-H at 1.09 and
    # costs +4.7 kcal/mol over a full turn, which is a single bond.
    #
    # And the torsion alone: a hand that turns did not ask for a bond length,
    # and the one it leaves behind -- an atom pushed rigidly through an arc,
    # so the bond bent rather than swung -- is not worth defending.  Held
    # anyway it is defended, and a cyclohexane pushed through its own plane
    # came back with ring bonds squeezed from 1.53 to 1.37 A.  Let go, they
    # relax to what they should be, the flip costs +11.1 kcal/mol against a
    # measured 10.7, and the ring stays a ring.
    # Heavy first here too.  The nearest bond from a carbon is always one of
    # its own hydrogens, and a hydrogen has no second bond to walk out along,
    # so the torsion could not be built at all -- every turn fell through to
    # holding distances, and a cyclohexane could not be flipped.
    anchor = min(((1 if (light[i] or light[j]) else 0),
                  math.dist(where[i], where[j]), i, j)
                 for i in inside for j in outside
                 if _is_a_bond(where, symbols, i, j)) \
        if any(_is_a_bond(where, symbols, i, j)
               for i in inside for j in outside) else None
    if anchor is not None:
        _heavy, _span, i, j = anchor
        spun = _turn_holding(where, then, [i, j], held, symbols)
        if spun is not None:
            return [spun]
    # A hand that has stopped moving has changed nothing, and nothing is not
    # an answer: with no values held the price falls back to a single point on
    # the geometry as it stands, which is the number this whole thing exists
    # to stop using.  So a still hand holds what a fresh one would.
    return contacts_holding(xyz_text, dragged, most=most)



def _turn_holding(where, then, contact, dragged, symbols):
    """A dihedral to hold, if the hand turned this contact rather than pulled.

    ``None`` when it pulled, when there are not four atoms to make a torsion
    out of, or when nothing moved.

    The other two atoms have to be *bonded* ones, walked out along the chain,
    not merely the nearest in space.  Taken by distance they came back as two
    hydrogens on a neighbouring carbon, and the torsion through those barely
    moves while the chain swings right round: 0.14 degrees where the real one
    turned ten.  A torsion is a statement about connectivity, and picking its
    atoms by proximity is not a smaller version of that.
    """
    count = len(symbols)
    i, j = int(contact[0]), int(contact[1])

    from delfin.atom_mapping import cov_radius
    radius = [cov_radius(one) for one in symbols]

    def bonded_to(of, without):
        near = [n for n in range(count)
                if n != of and n not in without and n not in dragged
                and math.dist(where[of], where[n])
                < 1.25 * (radius[of] + radius[n])]
        # A torsion through the heavy atoms is the one that names a
        # conformer; through a hydrogen it names how a methyl happens to be
        # standing.
        heavy = [n for n in near
                 if str(symbols[n]).strip().capitalize() != 'H']
        pick = heavy or near
        return min(pick, key=lambda n: math.dist(where[of], where[n])) \
            if pick else None

    k = bonded_to(j, {i})
    if k is None:
        return None
    ell = bonded_to(k, {i, j})
    if ell is None:
        return None

    now = _dihedral(where, i, j, k, ell)
    if _turned_by(now, _dihedral(then, i, j, k, ell)) < _TURN_DEGREES:
        return None
    return {'kind': 'dihedral', 'atoms': [i, j, k, ell], 'mode': 'drag',
            'value': now}


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


def solvent_note(solvent: Any, model: Any = 'alpb') -> str:
    """Which solvent a result is about, and under which model.

    A geometry optimised in the gas phase and one optimised in water are two
    different answers to two different questions, and a result that does not
    say which it is invites them to be compared.  So is a geometry from ALPB
    and one from ddCOSMO: on a glycine in water the two disagreed by 2.8
    kcal/mol, which is more than most of what they are used to decide.
    """
    return _solvents.note(model, solvent)


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


#: How many cores one interactive run may use.
#:
#: One was the rule, on the grounds that a login node is shared -- and it is,
#: but a single Optimise press is seconds of one person's waiting rather than
#: a job. Measured on a 74-atom cholesterol: a GFN2 single point 133 ms on one
#: core and 81 on four, five optimisation steps 669 ms and 510. GFN-FF gains
#: nothing -- 23 ms against 19, and worse at eight -- because at that size 12
#: of those milliseconds are the program starting up.
#:
#: Four, not all of them: the same machine is running everybody's dashboard,
#: and beyond four the curve is flat anyway.
INTERACTIVE_CORES = 4


#: And what a bigger one gets.  More cores are not more speed: xtb stops
#: gaining well before a large machine runs out, and past that it loses.
#: Measured on this one, twenty GFN2 cycles, best of three:
#:
#:     atoms      1     2     4     8    16    32    64
#:        18   0.37  0.27  0.22  0.27  0.27  0.37
#:        49   1.47  1.37  1.28  0.87  0.87  1.48
#:        62   0.72  0.67  0.57  0.57  0.47  0.67
#:       101               3.44  2.68  2.28  2.59  7.31
#:       185              13.38 13.03 11.08 11.85 17.49
#:
#: Sixteen is the top of it even at 185 atoms, and sixty-four is three times
#: worse than four.  A machine with 384 cores does not change that -- the work
#: in one small SCC does not divide that far, and what is left is threads
#: waiting for each other.  So the ladder is by size, and it ends.
_CORE_LADDER = ((30, 4), (60, 8))
_CORE_CEILING = 16


def interactive_cores(atoms: Optional[int] = None) -> int:
    """Cores for one interactive run, never more than the machine has.

    *atoms* is how big the structure is; without it the small-system figure
    stands, which is what an interactive drag mostly is.
    """
    said = str(os.environ.get('DELFIN_GFN_CORES', '')).strip()
    if said.isdigit() and int(said) > 0:
        # Whoever sets this has said what they want, including more than the
        # measurement recommends.
        wanted = int(said)
    elif atoms is None:
        wanted = INTERACTIVE_CORES
    else:
        wanted = _CORE_CEILING
        for size, cores in _CORE_LADDER:
            if int(atoms) < size:
                wanted = cores
                break
    return max(1, min(wanted, os.cpu_count() or 1))


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
    solvation_model: str = 'alpb',
    optimise: bool = True,
) -> Dict[str, Any]:
    """Relax *xyz_text* with xtb and say what happened.

    With *optimise* false this is a single point: the energy of the geometry
    as it stands, and the geometry handed straight back.  Everything else --
    which binary, which parameters, how many cores, the solvent, the charge
    and the spin -- is the same, which is the reason it lives here rather than
    in a second function that would drift from this one.  A single point is
    what a thermal budget is anchored on: the energy of the structure the user
    switched the mode on with, before anything has been moved.

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

    *solvent* is one of :data:`SOLVENTS` and *solvation_model* one of ALPB,
    GBSA or ddCOSMO -- see :mod:`delfin.dashboard.solvents` for which method
    can be run with which, and what each was measured to cost.  A geometry
    optimised in the gas phase and one optimised in water are different
    answers, and so are two optimised under different models, so both are
    named in the status rather than left in the operator's memory.
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
    model = str(solvation_model or 'alpb').strip().lower()
    if wet and spec.get('solvation') is False:
        return {
            'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
            'seconds': 0.0, 'frames': [],
            'status': (f'{label} has no implicit solvation in this build: '
                       'ALPB, GBSA and CPCM-X all stop it, and COSMO only '
                       'writes a file. Optimise it in the gas phase, or '
                       'choose GFN2-xTB or GFN-FF for a solvent.'),
        }
    # Every way this combination can be wrong, said before the run.  Three of
    # the four produce no error from xtb itself: an unparametrised GBSA solvent
    # stops it with a message about a file, and ddCOSMO under GFN-FF does not
    # stop it at all -- it returns a destroyed structure.
    no = _solvents.refusal(model, wet, key)
    if no:
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [], 'status': f'{label}: {no}'}

    ready = ensure_binary(key)
    binary = ready['path'] if ready['ok'] else find_binary(key)
    if binary is not None and not ready['ok'] and ready.get('status'):
        # It is there and it cannot do the job -- say so before it is used,
        # rather than after it has failed in its own words.
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [], 'status': ready['status']}
    if binary is None:
        wanted = ('its own xtb build, xtb-gxtb -- an ordinary xtb accepts '
                  '--gxtb and silently runs GFN2 instead'
                  if spec.get('binary') == 'gxtb' else 'xtb')
        return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': 0.0, 'frames': [],
                'status': (f'{label} needs {wanted}, which was not found in '
                           f'{_where_it_looked()}.'
                           + (f' {ready["status"]}' if ready.get('status') else ''))}

    started = time.perf_counter()
    folder = Path(tempfile.mkdtemp(prefix='delfin-gfn-'))
    try:
        source = folder / 'input.xyz'
        # Written from the lines, with a header that counts them.
        body = atom_lines(xyz_text)
        source.write_text(f'{len(body)}\nfrom the DELFIN viewer\n'
                          + '\n'.join(body) + '\n', encoding='utf-8')
        cores = interactive_cores(len(body))
        command = [binary, source.name, *spec['flags'],
                   *(['--opt'] if optimise else []),
                   '--chrg', str(int(charge)), '--uhf', str(max(0, int(uhf))),
                   '-P', str(cores)]
        if max_steps:
            command += ['--cycles', str(int(max_steps))]
        command += _solvents.xtb_flags(model, wet)
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
        # A few cores, and a scratch directory of its own: two of these must
        # not fight over xtbopt.xyz.
        environment = dict(os.environ, OMP_NUM_THREADS=str(cores),
                           MKL_NUM_THREADS=str(cores), OMP_STACKSIZE='1G')
        # And the parameters that belong to this binary, not whatever is lying
        # around. xtb reads a `param_gfn2-xtb.txt` from XTBPATH or from the
        # home directory *instead of* the parameters compiled into it, and a
        # truncated one there turns every GFN2 run into
        #
        #     no basis found for atom 1 Z= 8
        #     ERROR STOP / Error termination. Backtrace:
        #
        # while GFN-FF, which does not read them, goes on working -- measured,
        # and the exact report a user sent in. Pointing XTBPATH at the share
        # directory beside the binary restored the right answer with the bad
        # file still in place, so that is what is done here.
        own_parameters = parameter_home(binary)
        if own_parameters:
            environment['XTBPATH'] = own_parameters
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
                    # file has actually grown.
                    #
                    # Asking whether it grew costs 0.002 ms, so it is asked
                    # every pass; reading it costs 1.73 ms for a 211 KiB log of
                    # 37 frames, so that is what the rate limit is for. It used
                    # to be five times a second, which held a GFN-FF
                    # optimisation to five frames on screen per second while it
                    # was computing a step every fourteen milliseconds -- the
                    # limit was the watching, not the calculating.
                    try:
                        size = log.stat().st_size if log.exists() else -1
                    except OSError:
                        size = -1
                    if (on_frames is not None and size != last_size
                            and waited - last_read >= FRAME_READ_INTERVAL):
                        last_read = waited
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
                    time.sleep(WATCH_INTERVAL)
                    waited += WATCH_INTERVAL
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
        # A single point writes no xtbopt.xyz and no path: the geometry that
        # went in is the geometry that comes back, and there is nothing to
        # play.  Reading them anyway would find the leftovers of whatever ran
        # in this directory before, which is nothing here, so it would hand
        # back None and read as a failure.
        relaxed = _read_optimised(folder, xyz_text) if optimise else xyz_text
        frames = read_trajectory(folder) if optimise else []
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
            return {'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                    'seconds': seconds, 'frames': [], 'output': output[-4000:],
                    'binary': str(binary), 'version': _said_version(output),
                    'status': (f'{label}: {why_it_stopped(output)} The '
                               'structure was left as it was. '
                               f'{which_xtb_ran(binary, output)}')}

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
            return {
                'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': time.perf_counter() - started, 'engine': 'xtb',
                'version': version, 'hamiltonian': reported or wanted,
                'frames': [],
                'output': output[-4000:], 'binary': str(binary),
                'status': (f'{label} stopped with an error: '
                           f'{why_it_stopped(output)} The structure was left '
                           f'as it was. {which_xtb_ran(binary, output)}'),
            }

        # A single point has no geometry to converge, and calling it
        # unconverged would put it among the runs that stopped short.
        converged = (not optimise) or 'GEOMETRY OPTIMIZATION CONVERGED' in output
        if not converged:
            # The geometry is still better than the one that went in, so it is
            # handed back -- but not as though it were finished.
            return {
                'ok': True, 'xyz': relaxed, 'energy': energy, 'method': key,
                'seconds': seconds, 'engine': 'xtb', 'frames': frames, 'version': version,
                'hamiltonian': reported or wanted or label, 'held': held,
                'converged': False, 'solvent': wet, 'solvation_model': model,
                'status': (f'{label} stopped before converging after '
                           f'{seconds:.1f} s; the geometry it reached is shown. '
                           f'(xtb {version}, {reported or wanted or label})'
                           + solvent_note(wet, model) + held_note(held)),
            }
        return {
            'ok': True, 'xyz': relaxed, 'energy': energy, 'method': key,
            'seconds': seconds, 'engine': 'xtb', 'frames': frames, 'version': version,
            'hamiltonian': reported or wanted or label, 'held': held,
            'converged': True, 'solvent': wet, 'solvation_model': model,
            'status': (f'{label} converged in {seconds:.1f} s '
                       f'(xtb {version}, {reported or wanted or label}).'
                       + solvent_note(wet, model) + held_note(held)),
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
    solvation_model: str = 'alpb',
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
        solvation_model=solvation_model,
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
