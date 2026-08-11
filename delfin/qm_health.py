"""Whether a tool is installed is not whether it works.

An installer that has run is not a tool that answers.  Everything that went
wrong for the users of this dashboard went wrong *after* a green banner:

* a binary copied out of the environment holding the libraries it is linked
  against -- present, executable, and unable to start;
* a ``param_gfn2-xtb.txt`` left in a home directory, which xtb reads instead
  of the parameters compiled into it, so every GFN2 run ends in
  ``Error termination. Backtrace:`` while GFN-FF goes on working;
* an ordinary xtb accepting ``--gxtb``, warning once, and quietly running
  GFN2 -- the answer differs by 71 Eh on water and nothing said so;
* an ``anmr`` whose every invocation ends in ``SIGSEGV`` being reported as
  ``status: ok`` with the crash message as its version.

So the question this module answers is not "is it there" but "does it give the
right answer", and when it does not, "what exactly is wrong and what would fix
it".  Every expected value below was measured against the build named beside
it on a machine where the tool was working -- none of it is from a manual.

The tools are found by :mod:`delfin.qm_runtime` and by nobody else here: a
second way to find a binary is how a repository ends up with nine of them.
"""

from __future__ import annotations

import os
import re
import resource
import shutil
import subprocess
import threading
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence

__all__ = [
    'ToolHealth', 'PROBES', 'known_tools', 'probe_environment',
    'check_tool', 'check_tools', 'repair_actions', 'repair_command',
    'repair_tool', 'ensure_tool', 'ensure_package', 'provide', 'history',
    'PACKAGES',
    'format_health',
]

#: Three atoms, so a probe costs milliseconds rather than seconds.
WATER = ('3\nwater\n'
         'O 0.00000000 0.00000000 0.00000000\n'
         'H 0.96000000 0.00000000 0.00000000\n'
         'H -0.24000000 0.93000000 0.00000000\n')

#: Two atoms and the smallest basis there is: enough to prove that ORCA finds
#: its basis-set library, which is the half of an ORCA install that breaks.
H2_INPUT = '! HF STO-3G SP\n\n*xyz 0 1\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n*\n'

_XTB_ENERGY = re.compile(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh')
_ORCA_ENERGY = re.compile(r'FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)')
_CREST_ENERGY = re.compile(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh')

#: What each tool is asked, and what a working one answers.
#:
#: ``energy`` is the number the probe must reproduce and ``tolerance`` how far
#: it may drift between builds -- 1e-4 Eh, which is far below any chemistry
#: and far above the last digits of a different compiler.
PROBES: Dict[str, Dict[str, Any]] = {
    'xtb': {
        'label': 'xtb',
        'version_args': ['--version'],
        'version_re': re.compile(r'xtb version\s+([0-9][0-9.]*)'),
        # GFN2 first, because GFN2 is what reads the parameter files and so
        # the only one of the two that a stray file can break.
        'answers': [
            {'args': ['in.xyz', '--gfn', '2', '--sp'], 'energy': -5.0703,
             'what': 'GFN2 single point on water'},
            {'args': ['in.xyz', '--gfnff', '--sp'], 'energy': -0.3273,
             'what': 'GFN-FF single point on water'},
        ],
        'files': {'in.xyz': WATER},
        'energy_re': _XTB_ENERGY,
        'tolerance': 1e-3,
        'measured_with': 'xtb 6.7.1',
    },
    'gxtb': {
        'label': 'g-xTB',
        'version_args': ['--version'],
        'version_re': re.compile(r'xtb version\s+([0-9][0-9.]*)'),
        'answers': [
            {'args': ['in.xyz', '--gxtb', '--sp'], 'energy': -76.4375,
             'what': 'g-xTB single point on water'},
        ],
        'files': {'in.xyz': WATER},
        'energy_re': _XTB_ENERGY,
        'tolerance': 1e-2,
        # An ordinary xtb takes --gxtb, warns once and runs GFN2 instead. The
        # energy is the only thing that says so, and it is 71 Eh out.
        'impostor_energy': -5.0703,
        'impostor_why': ('this is an ordinary xtb: it ignored --gxtb and ran '
                         'GFN2, which answers 71 Eh away from g-xTB'),
        'measured_with': 'g-xTB build 26dd68d',
    },
    'crest': {
        'label': 'CREST',
        'version_args': ['--version'],
        'version_re': re.compile(r'[Vv]ersion\s+([0-9][0-9.]*)'),
        'answers': [
            {'args': ['in.xyz', '--sp'], 'energy': -5.0703,
             'what': 'single point on water'},
        ],
        'files': {'in.xyz': WATER},
        'energy_re': _CREST_ENERGY,
        'tolerance': 1e-3,
        'measured_with': 'CREST 3.0.2',
    },
    'orca': {
        'label': 'ORCA',
        # It has no --version: given one it tries to open it as an input file
        # and exits 2. Reading the version out of the binary with `strings`
        # costs longer than running a job, so a job is run.
        'version_args': [],
        'version_re': re.compile(r'Program Version\s+([0-9][0-9.]*)'),
        'answers': [
            {'args': ['h2.inp'], 'energy': -1.1168,
             'what': 'HF/STO-3G on H2', 'marker': 'ORCA TERMINATED NORMALLY'},
        ],
        'files': {'h2.inp': H2_INPUT},
        'energy_re': _ORCA_ENERGY,
        'tolerance': 1e-3,
        'absolute_path': True,   # a parallel ORCA refuses a bare name
        'measured_with': 'ORCA 6.1.1',
    },
    'mopac': {
        'label': 'MOPAC',
        'version_args': [],
        'version_re': re.compile(r'MOPAC\s+v?([0-9][0-9.]*)'),
        # It wants an input file and writes its answer beside it, so the probe
        # is a whole tiny job -- which also proves it can write where it runs.
        'answers': [
            {'args': ['in.mop'], 'energy': -50.06,
             'what': 'PM6-D3H4 on water', 'marker': 'FINAL HEAT OF FORMATION'},
        ],
        'files': {'in.mop': ('PM6-D3H4 PRECISE\nprobe\n\n'
                             'O 0.00000 1 0.00000 1 0.00000 1\n'
                             'H 0.96000 1 0.00000 1 0.00000 1\n'
                             'H -0.24000 1 0.93000 1 0.00000 1\n')},
        'energy_re': re.compile(r'FINAL HEAT OF FORMATION\s*=\s*(-?\d+\.\d+)'),
        'reads': 'in.out',
        'tolerance': 5.0,
        'measured_with': 'MOPAC 23.2.5',
    },
    'anmr': {
        'label': 'ANMR',
        'runnable': False,
        # Measured: `anmr` and `anmr --help` both end in
        # `forrtl: severe (174): SIGSEGV`. Running it is not evidence either
        # way, so it is never run and never called broken for crashing.
        'why_no_probe': ('anmr has no probe mode -- every call without a full '
                         'ENSO directory ends in a segmentation fault, so it '
                         'is checked for being there and no further'),
    },
}

#: What the DELFIN installer can fetch, and therefore what a missing one can
#: be offered a repair for. ORCA, Turbomole and Multiwfn are not here on
#: purpose: they are licensed and are installed by hand.
INSTALLABLE = ('xtb', 'gxtb', 'crest', 'dftb+', 'xtb4stda', 'stda', 'std2',
                'mopac')

#: Tools that are found rather than probed: present or not, and nothing is
#: claimed about them beyond that.
PRESENCE_ONLY = ('dftb+', 'stda', 'xtb4stda', 'std2', 'packmol', 'censo',
                 'c2anmr', 'nmrplot', 'gnrs', 'multiwfn')


@dataclass(frozen=True)
class ToolHealth:
    """What is known about one tool, and what to do about it."""

    name: str
    label: str = ''
    present: bool = False
    path: str = ''
    source: str = ''
    runnable: bool = False
    version: str = ''
    healthy: bool = False
    level: str = 'absent'        # ok | warn | fail | absent
    why: str = ''
    fix: str = ''
    repair: str = ''
    seconds: float = 0.0
    evidence: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            'name': self.name, 'label': self.label or self.name,
            'present': self.present, 'path': self.path, 'source': self.source,
            'runnable': self.runnable, 'version': self.version,
            'healthy': self.healthy, 'level': self.level, 'why': self.why,
            'fix': self.fix, 'repair': self.repair,
            'seconds': round(self.seconds, 3), 'evidence': dict(self.evidence),
        }


def known_tools() -> tuple:
    """Every tool this module can say something about."""
    return tuple(PROBES) + PRESENCE_ONLY


# ---------------------------------------------------------------------------
# the environment the tools will run in
# ---------------------------------------------------------------------------

def _home_parameter_files(home: Optional[str] = None) -> List[str]:
    """Parameter files in the home directory, which xtb prefers to its own."""
    root = Path(home or os.path.expanduser('~'))
    found = []
    try:
        for entry in root.iterdir():
            name = entry.name
            if name.startswith('param_gfn') and name.endswith('.txt'):
                found.append(str(entry))
            elif name == '.xtbrc':
                found.append(str(entry))
    except OSError:
        return []
    return sorted(found)


def probe_environment(env: Optional[Mapping[str, str]] = None) -> Dict[str, Any]:
    """The conditions the tools will run under.  Runs nothing.

    Three of these were measured to matter and are not guesses:

    * **A parameter file in the home directory, or on XTBPATH, is used instead
      of the ones compiled into xtb.**  A truncated one turns every GFN2 run
      into ``no basis found for atom 1`` and ``Error termination. Backtrace:``
      while GFN-FF keeps working.
    * **The stack limit, not OMP_STACKSIZE, is what kills large jobs.**  A
      122-atom GFN2 Hessian died with signal 11 and no message at a 512 KiB
      soft limit and finished at 2 MiB.  A water single point passes at
      512 KiB, so no cheap probe finds this -- it is read, not tested.
    * **An unset thread count means one thread per hardware thread.**  On a
      384-core machine ``xtb --version`` alone spent 2.1 CPU seconds before
      printing a line.
    """
    source = dict(env if env is not None else os.environ)
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    xtbpath = source.get('XTBPATH', '')
    on_path = []
    for part in xtbpath.split(os.pathsep):
        if not part:
            continue
        candidate = Path(part) / 'param_gfn2-xtb.txt'
        if candidate.is_file():
            on_path.append(str(candidate))
    home_files = _home_parameter_files(source.get('HOME'))

    warnings: List[str] = []
    if home_files:
        warnings.append(
            'xtb reads ' + ', '.join(Path(f).name for f in home_files) +
            ' from the home directory instead of its own parameters; a '
            'damaged one there fails every GFN2 run with a backtrace while '
            'GFN-FF still works')
    if soft not in (resource.RLIM_INFINITY,) and soft < 2 * 1024 * 1024:
        warnings.append(
            f'the stack limit is {soft // 1024} KiB; a Hessian on a hundred '
            'atoms was measured to die with no message below about 2 MiB')
    if not source.get('OMP_NUM_THREADS'):
        warnings.append(
            'no thread count is set, so each run takes one thread per core')

    return {
        'threads': source.get('OMP_NUM_THREADS', ''),
        'stack_soft': soft, 'stack_hard': hard,
        'cores': os.cpu_count() or 0,
        'xtbpath': xtbpath,
        'xtbpath_parameters': on_path,
        'home_parameters': home_files,
        'warnings': warnings,
    }


# ---------------------------------------------------------------------------
# finding and running one tool
# ---------------------------------------------------------------------------

def find_tool(name: str) -> Dict[str, Any]:
    """Where a tool is, asked of the one resolver this repository has.

    A *broken* link is reported as such rather than as nothing.  A dangling
    ``qm_tools/bin/xtb`` is what the installer leaves behind when the
    environment it pointed into is removed, and ``is_file()`` follows the
    link, so today that reads as "xtb was not found" while ``ls`` plainly
    shows an xtb.  A user cannot act on that; "the link points at a path that
    is gone" they can.
    """
    path, source = '', ''
    if name == 'gxtb':
        # Asked for first, not as a fallback: a plain `gxtb` on the PATH is
        # the raw driver, while what DELFIN runs is an xtb build that takes
        # --gxtb and is resolved as xtb-gxtb. Finding the wrong one here would
        # report the wrong program healthy.
        try:
            from delfin.dashboard.gfn_optimize import find_gxtb

            found = find_gxtb()
            if found:
                return {'path': str(found), 'source': 'qm_tools', 'dangling': ''}
        except Exception:
            pass
    try:
        from delfin.qm_runtime import resolve_tool

        resolved = resolve_tool(name)
        if resolved is not None:
            path = str(getattr(resolved, 'path', '') or '')
            source = str(getattr(resolved, 'source', '') or '')
    except Exception:
        resolved = None
    if not path:
        found = shutil.which(name)
        if found:
            path, source = found, 'PATH'
    dangling = ''
    if not path:
        try:
            from delfin.qm_runtime import get_qm_tools_bin_dir

            link = Path(get_qm_tools_bin_dir()) / name
            if link.is_symlink() and not link.exists():
                dangling = os.readlink(link)
        except Exception:
            dangling = ''
    return {'path': path, 'source': source, 'dangling': dangling}


def _run(command: Sequence[str], *, cwd: Optional[str] = None,
         timeout: float = 60.0,
         env: Optional[Mapping[str, str]] = None) -> Dict[str, Any]:
    """One tool, one question, and never in the caller's directory.

    xtb leaves ``charges``, ``wbo``, ``xtbrestart`` and ``xtbtopo.mol`` behind,
    and GFN-FF adds ``gfnff_topo`` and ``gfnff_charges``.  A probe run where
    the user works drops those into their project.
    """
    started = time.perf_counter()
    settings = dict(env if env is not None else os.environ)
    # One thread: a probe that spreads over every core to answer in 20 ms is
    # a probe that costs more than it measures.
    settings.setdefault('OMP_NUM_THREADS', '1')
    settings.setdefault('MKL_NUM_THREADS', '1')
    try:
        done = subprocess.run(
            list(command), cwd=cwd, env=settings, timeout=timeout,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return {'exit': done.returncode, 'stdout': done.stdout,
                'stderr': done.stderr, 'seconds': time.perf_counter() - started}
    except subprocess.TimeoutExpired:
        return {'exit': None, 'stdout': '', 'stderr': f'no answer in {timeout:.0f} s',
                'timeout': True, 'seconds': time.perf_counter() - started}
    except OSError as problem:
        return {'exit': None, 'stdout': '', 'stderr': str(problem),
                'seconds': time.perf_counter() - started}


def _version_of(name: str, path: str, probe: Mapping[str, Any],
                timeout: float) -> Dict[str, Any]:
    args = list(probe.get('version_args') or [])
    if not args:
        return {'version': '', 'runnable': True, 'why': ''}
    answer = _run([path] + args, timeout=timeout)
    text = (answer.get('stdout') or '') + (answer.get('stderr') or '')
    pattern = probe.get('version_re')
    found = pattern.search(text) if pattern else None
    if answer.get('exit') is None:
        return {'version': '', 'runnable': False,
                'why': (answer.get('stderr') or 'it would not start').strip()}
    if 'error while loading shared libraries' in text:
        library = text.split('error while loading shared libraries:')[-1]
        return {'version': '', 'runnable': False,
                'why': ('it cannot start: ' + library.strip().split('\n')[0]),
                'repair': 'relink'}
    return {'version': found.group(1) if found else '', 'runnable': True,
            'why': ''}


def check_tool(name: str, *, depth: str = 'answer', timeout: float = 60.0,
               env: Optional[Mapping[str, str]] = None) -> ToolHealth:
    """Check one tool as deeply as *depth* asks.

    ``present`` locates it and runs nothing, ``runs`` starts it and reads a
    version, ``answer`` gives it a question whose answer is known.  The three
    differ by two orders of magnitude in cost, so the caller chooses.
    """
    started = time.perf_counter()
    probe = PROBES.get(name, {})
    label = str(probe.get('label') or name)
    where = find_tool(name)
    path, source = where['path'], where['source']

    if not path:
        if where['dangling']:
            return ToolHealth(
                name=name, label=label, present=False, source='qm_tools',
                level='fail', seconds=time.perf_counter() - started,
                why=(f'the link in the tool directory points at '
                     f'{where["dangling"]}, which is not there any more'),
                fix='point it at an installed binary, or install it again',
                repair='relink',
                evidence={'dangling': where['dangling']})
        return ToolHealth(name=name, label=label, level='absent',
                          seconds=time.perf_counter() - started,
                          why='not found', fix='install it',
                          repair='install' if name in INSTALLABLE else '')

    if depth == 'present' or probe.get('runnable') is False:
        return ToolHealth(
            name=name, label=label, present=True, path=path, source=source,
            runnable=bool(probe.get('runnable', True)),
            level='ok' if probe.get('runnable') is not False else 'warn',
            why=str(probe.get('why_no_probe') or ''),
            seconds=time.perf_counter() - started)

    told = _version_of(name, path, probe, timeout)
    if not told['runnable']:
        return ToolHealth(
            name=name, label=label, present=True, path=path, source=source,
            runnable=False, level='fail', why=told['why'],
            fix='install it again so it stands with its own libraries',
            repair=str(told.get('repair') or 'install'),
            seconds=time.perf_counter() - started)

    health = ToolHealth(
        name=name, label=label, present=True, path=path, source=source,
        runnable=True, version=told['version'], healthy=True, level='ok',
        seconds=time.perf_counter() - started)
    if depth != 'answer' or not probe.get('answers'):
        return health

    with tempfile.TemporaryDirectory(prefix='delfin-health-') as scratch:
        for filename, content in (probe.get('files') or {}).items():
            Path(scratch, filename).write_text(content, encoding='utf-8')
        for question in probe['answers']:
            command = [path] + list(question['args'])
            answer = _run(command, cwd=scratch, timeout=timeout, env=env)
            text = (answer.get('stdout') or '') + '\n' + (answer.get('stderr') or '')
            # Some of them answer in a file rather than on the terminal --
            # MOPAC writes its whole report beside the input and says nothing
            # on stdout, so reading only what was printed found nothing and
            # called a working install broken.
            written = question.get('reads') or probe.get('reads')
            if written:
                try:
                    text += '\n' + Path(scratch, written).read_text(
                        encoding='utf-8', errors='replace')
                except OSError:
                    pass
            evidence = {
                'command': ' '.join(command), 'exit': answer.get('exit'),
                'what': question.get('what', ''),
                'tail': '\n'.join(
                    line for line in text.splitlines() if line.strip())[-600:],
            }
            if not health.version and probe.get('version_re'):
                said = probe['version_re'].search(text)
                if said:
                    told['version'] = said.group(1)
            found = probe['energy_re'].search(text)
            value = float(found.group(1)) if found else None
            evidence['energy'] = value
            marker = question.get('marker')
            if marker and marker not in text:
                return _unhealthy(health, name, probe,
                                  f'{question["what"]} did not finish',
                                  evidence, started)
            if value is None:
                return _unhealthy(health, name, probe,
                                  f'{question["what"]} gave no energy',
                                  evidence, started)
            impostor = probe.get('impostor_energy')
            if impostor is not None and abs(value - impostor) < probe['tolerance']:
                return ToolHealth(
                    name=name, label=label, present=True, path=path,
                    source=source, runnable=True, version=told['version'],
                    healthy=False, level='fail',
                    why=str(probe.get('impostor_why') or 'the wrong binary'),
                    fix='install the g-xTB build; an ordinary xtb cannot do it',
                    repair='install',
                    seconds=time.perf_counter() - started, evidence=evidence)
            if abs(value - question['energy']) > probe['tolerance']:
                return _unhealthy(
                    health, name, probe,
                    (f'{question["what"]} answered {value:.4f} Eh where '
                     f'{question["energy"]:.4f} was expected'),
                    evidence, started)
            health = ToolHealth(
                name=name, label=label, present=True, path=path, source=source,
                runnable=True, version=told['version'], healthy=True,
                level='ok', seconds=time.perf_counter() - started,
                evidence=evidence)
    if name == 'xtb' and health.healthy:
        # It works -- but say if a broken one is standing in front of it,
        # because that is a trap set for the next tool, the next session, or
        # the colleague who resolves it differently.
        try:
            from delfin.dashboard.gfn_optimize import unusable_xtb_note

            note = unusable_xtb_note()
        except Exception:
            note = ''
        if note:
            return ToolHealth(
                name=name, label=health.label, present=True, path=health.path,
                source=health.source, runnable=True, version=health.version,
                healthy=True, level='warn', why=note,
                fix='remove it, or install over it, so nothing can find it',
                repair='install', seconds=health.seconds,
                evidence=health.evidence)
    return health


def _unhealthy(health: ToolHealth, name: str, probe: Mapping[str, Any],
               why: str, evidence: Dict[str, Any], started: float) -> ToolHealth:
    """A tool that started and then gave the wrong answer, and why that is."""
    fix = 'install it again'
    repair = 'install'
    tail = str(evidence.get('tail') or '')
    if name == 'xtb' and ('no basis found' in tail or 'Backtrace' in tail):
        # The one that has actually happened to a user: xtb reading a
        # parameter file that is not its own.
        fix = ("point xtb at its own parameters -- a param_gfn*-xtb.txt in "
               "the home directory or on XTBPATH is read instead of them")
        repair = 'xtbpath'
    return ToolHealth(
        name=name, label=health.label, present=True, path=health.path,
        source=health.source, runnable=True, version=health.version,
        healthy=False, level='fail', why=why, fix=fix, repair=repair,
        seconds=time.perf_counter() - started, evidence=evidence)


def check_tools(names: Optional[Sequence[str]] = None, *,
                depth: str = 'present', workers: int = 8,
                timeout: float = 60.0) -> List[ToolHealth]:
    """Several at once.

    Separate processes, so they overlap: eight probes were measured at 1.88 s
    one after another and 0.71 s together, after which the wait is whatever
    the slowest one costs on its own.
    """
    wanted = list(names) if names is not None else list(known_tools())
    if not wanted:
        return []
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as pool:
        return list(pool.map(
            lambda name: check_tool(name, depth=depth, timeout=timeout), wanted))


# ---------------------------------------------------------------------------
# putting it right
# ---------------------------------------------------------------------------

#: What each repair does, in the words a user reads before agreeing to it.
REPAIRS: Dict[str, str] = {
    'relink': 'point the tool directory at a binary that is actually there',
    'xtbpath': 'run xtb with the parameters that belong to it',
    'install': 'install the tool with the DELFIN installer',
}


def repair_actions(health: ToolHealth) -> List[str]:
    return [health.repair] if health.repair else []


def repair_command(name: str, action: str) -> Optional[List[str]]:
    """Exactly what a repair would run, so it can be read before it is run.

    The same contract as the installer preview the Submit tab already shows:
    nothing here happens without the user having seen the command first.
    """
    if action == 'install':
        try:
            from delfin.dashboard.gfn_optimize import install_command

            return install_command(tool=name)
        except Exception:
            return None
    if action == 'relink':
        where = find_tool(name)
        target = where['path'] or shutil.which(name)
        if not target:
            return None
        try:
            from delfin.qm_runtime import get_qm_tools_bin_dir

            return ['ln', '-sfn', target, str(Path(get_qm_tools_bin_dir()) / name)]
        except Exception:
            return None
    if action == 'xtbpath':
        try:
            from delfin.dashboard.gfn_optimize import find_xtb, parameter_home

            home = parameter_home(find_xtb())
            return ['export', f'XTBPATH={home}'] if home else None
        except Exception:
            return None
    return None


def repair_tool(name: str, action: str = '', *,
                on_line: Optional[Callable[[str], None]] = None,
                timeout: float = 1800.0) -> Dict[str, Any]:
    """Attempt one repair and say what happened.

    The tool is checked again afterwards, because a repair that is not
    re-proven is a claim.  ``health`` in the answer is that second check.
    """
    action = action or (check_tool(name, depth='answer').repair or '')
    if not action:
        return {'ok': False, 'action': '', 'status': 'nothing to repair',
                'lines': [], 'health': check_tool(name, depth='answer').as_dict()}
    command = repair_command(name, action)
    if not command:
        return {'ok': False, 'action': action, 'lines': [],
                'status': f'{action} cannot be done here',
                'health': check_tool(name, depth='answer').as_dict()}

    lines: List[str] = []
    if action == 'xtbpath':
        # Nothing is run and nothing in the home directory is touched: the
        # parameter file the user has may be there on purpose. DELFIN points
        # xtb at its own parameters when it starts it, which is where the
        # setting belongs.
        try:
            from delfin.dashboard.gfn_optimize import find_xtb, parameter_home

            home = parameter_home(find_xtb())
        except Exception:
            home = None
        if home:
            os.environ['XTBPATH'] = home
            lines.append(f'XTBPATH={home}')
    else:
        started = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1)
        try:
            for line in started.stdout or ():
                text = line.rstrip('\n')
                lines.append(text)
                if on_line is not None:
                    on_line(text)
            started.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            started.kill()
            lines.append(f'no answer in {timeout:.0f} s')
        finally:
            if started.stdout is not None:
                started.stdout.close()

    after = check_tool(name, depth='answer')
    return {
        'ok': after.healthy or after.level == 'ok',
        'action': action, 'lines': lines, 'health': after.as_dict(),
        'status': (f'{after.label} works now.' if after.healthy
                   else f'{after.label} still does not: {after.why}'),
    }


#: Tools this session has already tried to install by itself.
_TRIED: set = set()


def ensure_tool(name: str,
                on_line: Optional[Callable[[str], None]] = None,
                timeout: float = 1800.0) -> Dict[str, Any]:
    """The tool, installed if it is not there -- for anything the installer knows.

    The same rule as for xtb, applied to the rest: a tool that is needed and
    is not there is fetched once, with a word about why the wait is happening,
    and never a second time in one session.

    Licensed programs are not in :data:`INSTALLABLE` and are never fetched:
    ORCA, Turbomole, Gaussian and Multiwfn are downloaded by the person who
    holds the licence, and a program that goes and gets one on their behalf is
    not being helpful.
    """
    health = check_tool(name, depth='runs')
    if health.level == 'ok' and health.present:
        return {'ok': True, 'installed': False, 'health': health, 'status': ''}

    from delfin.dashboard.gfn_optimize import auto_install_allowed, install_xtb

    if name not in INSTALLABLE:
        return {'ok': False, 'installed': False, 'health': health,
                'status': (f'{health.label or name} is not there, and it is '
                           'not one DELFIN installs -- it is licensed and is '
                           'fetched by whoever holds the licence.')}
    if not auto_install_allowed():
        return {'ok': False, 'installed': False, 'health': health,
                'status': (f'{health.label or name}: {health.why or "not found"}. '
                           'Automatic installation is switched off.')}
    if name in _TRIED:
        return {'ok': False, 'installed': False, 'health': health,
                'status': (f'{health.label or name}: {health.why or "not found"}. '
                           'It was already installed once this session and is '
                           'still not right -- Settings shows the installer\'s '
                           'own output.')}

    _TRIED.add(name)
    if on_line is not None:
        on_line(f'{health.label or name}: {health.why or "not installed"}. '
                'Installing it -- a few minutes.')
    install_xtb(on_line=on_line, timeout=timeout, tool=name)
    after = check_tool(name, depth='runs')
    return {
        'ok': after.level == 'ok' and after.present, 'installed': True,
        'health': after,
        'status': (f'{after.label or name} was missing and has been installed.'
                   if after.present else
                   f'{after.label or name} could not be installed: {after.why}'),
    }


#: Python packages DELFIN can fetch, and what it takes to fetch each.
#:
#: ``family`` picks the installer, ``flag`` the switch inside it -- those
#: installers do every package they know unless told which one -- and ``size``
#: is said out loud before the wait, because "a few minutes" and "two
#: gigabytes over this network" are not the same sentence to a user on a
#: metered or slow connection.
PACKAGES: Dict[str, Dict[str, str]] = {
    'cclib':      {'family': 'analysis', 'flag': 'INSTALL_CCLIB',   'label': 'cclib'},
    'nglview':    {'family': 'analysis', 'flag': 'INSTALL_NGLVIEW', 'label': 'nglview'},
    'morfeus':    {'family': 'analysis', 'flag': 'INSTALL_MORFEUS', 'label': 'morfeus'},
    'torchani':   {'family': 'mlp', 'flag': 'INSTALL_TORCHANI', 'label': 'TorchANI',
                   'size': 'it brings PyTorch with it, which is a large download'},
    'aimnet2calc': {'family': 'mlp', 'flag': 'INSTALL_AIMNET2', 'label': 'AIMNet2',
                    'size': 'it brings PyTorch with it, which is a large download'},
    'mace':       {'family': 'mlp', 'flag': 'INSTALL_MACE', 'label': 'MACE',
                   'size': 'it brings PyTorch with it, which is a large download'},
    'chgnet':     {'family': 'mlp', 'flag': 'INSTALL_CHGNET', 'label': 'CHGNet',
                   'size': 'it brings PyTorch with it, which is a large download'},
}

_PACKAGES_TRIED: set = set()


def package_present(module: str) -> bool:
    import importlib.util

    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):
        return False


def ensure_package(module: str,
                   on_line: Optional[Callable[[str], None]] = None,
                   timeout: float = 3600.0) -> Dict[str, Any]:
    """A Python package DELFIN can use, installed if it is not there.

    Into the interpreter that will import it -- which is the whole point, and
    was the bug: the installers took the first python on the PATH, so packages
    were installed and missing at the same time.

    Never from a listing.  These are asked for where a package is about to be
    used, not where a status panel is drawn: a page that installs two
    gigabytes because somebody opened it is not a page anybody wants.
    """
    import importlib

    if package_present(module):
        return {'ok': True, 'installed': False, 'status': ''}

    spec = PACKAGES.get(module)
    label = (spec or {}).get('label', module)
    if spec is None:
        return {'ok': False, 'installed': False,
                'status': f'{label} is not one DELFIN installs.'}

    from delfin.dashboard.gfn_optimize import auto_install_allowed

    if not auto_install_allowed():
        return {'ok': False, 'installed': False,
                'status': f'{label} is not installed, and automatic '
                          'installation is switched off.'}
    if module in _PACKAGES_TRIED:
        return {'ok': False, 'installed': False,
                'status': (f'{label} was already installed once this session '
                           'and still cannot be imported -- Settings shows the '
                           'installer\'s own output.')}
    _PACKAGES_TRIED.add(module)

    size = spec.get('size')
    if on_line is not None:
        on_line(f'{label} is needed and not installed. Fetching it'
                + (f' -- {size}.' if size else ' -- a few minutes.'))

    from delfin.runtime_setup import (run_analysis_tools_installer,
                                      run_mlp_tools_installer)

    # Only the one that is wanted: these installers do everything they know
    # unless each switch is turned off by name.
    others = {other['flag']: '0' for other in PACKAGES.values()
              if other['family'] == spec['family']}
    others[spec['flag']] = '1'
    runner = (run_analysis_tools_installer if spec['family'] == 'analysis'
              else run_mlp_tools_installer)
    try:
        _target, result = runner(extra_env=others)
        output = getattr(result, 'stdout', '') or ''
    except Exception as problem:
        return {'ok': False, 'installed': True,
                'status': f'{label} could not be installed: {problem}'}

    # A package that has just arrived is not on the import path this process
    # already worked out.
    importlib.invalidate_caches()
    if package_present(module):
        return {'ok': True, 'installed': True,
                'status': f'{label} was missing and has been installed.'}
    return {'ok': False, 'installed': True,
            'status': (f'{label} could not be installed. '
                       + (output.strip().splitlines() or [''])[-1][:200])}


#: One lock per name, and a record of what was done under it.
#:
#: The dashboard answers on threads, and two of them wanting the same tool at
#: the same moment would have started two conda solves into one prefix -- which
#: does not end in two installs, it ends in a broken one. The second caller
#: waits for the first and then finds the work already done.
_LOCKS: Dict[str, Any] = {}
_LOCKS_GUARD = threading.Lock()
_HISTORY: List[Dict[str, Any]] = []


def _lock_for(name: str):
    with _LOCKS_GUARD:
        if name not in _LOCKS:
            _LOCKS[name] = threading.Lock()
        return _LOCKS[name]


def history() -> List[Dict[str, Any]]:
    """What was installed or repaired this session, in order.

    So the dashboard can say what it did on the user's behalf rather than
    leaving them to guess why something once took four minutes.
    """
    return list(_HISTORY)


def provide(what: str,
            on_line: Optional[Callable[[str], None]] = None,
            timeout: float = 3600.0) -> Dict[str, Any]:
    """Whatever is needed, made ready -- or a sentence saying why it cannot be.

    The one entry point.  A caller that is about to use something asks for it
    by name and gets back whether it can go ahead; what that took -- a binary
    fetched, a link repointed, a package installed into this interpreter, or
    nothing at all because it was already fine -- is this module's business
    and not the caller's.

    Three things can be wrong and they do not have the same answer, which is
    why asking "is it there" was never enough:

    * it is not there -- fetch it, if it is ours to fetch;
    * it is there and cannot do the job -- a binary that will not start, a
      link into an environment that has been deleted, an xtb too old to run an
      optimisation -- repair it, and prove the repair;
    * it is there and it works -- say so and get out of the way.

    Returns ``{'ok', 'action', 'status'}``.  ``action`` is one of ``''``,
    ``'installed'``, ``'repaired'`` -- what actually happened, so a caller can
    tell the user why it waited.
    """
    name = str(what or '').strip()
    if not name:
        return {'ok': False, 'action': '', 'status': 'nothing was asked for'}
    with _lock_for(name):
        answer = _provide_once(name, on_line, timeout)
    if answer.get('action'):
        _HISTORY.append({'what': name, 'action': answer['action'],
                         'status': answer.get('status', '')})
    return answer


def _provide_once(name: str, on_line, timeout: float) -> Dict[str, Any]:

    if name in PACKAGES or (name not in PROBES and name not in PRESENCE_ONLY
                            and package_present(name)):
        answer = ensure_package(name, on_line=on_line, timeout=timeout)
        return {'ok': answer['ok'],
                'action': 'installed' if answer.get('installed') else '',
                'status': answer.get('status', '')}

    health = check_tool(name, depth='runs')
    if health.level == 'ok' and health.present:
        return {'ok': True, 'action': '', 'status': ''}

    # There and unable: repair before reaching for the network. A repointed
    # link costs nothing and fixes the commonest of these outright.
    if health.present and health.repair and health.repair != 'install':
        if on_line is not None:
            on_line(f'{health.label or name}: {health.why} -- {health.fix}')
        done = repair_tool(name, health.repair, on_line=on_line)
        if done.get('ok'):
            return {'ok': True, 'action': 'repaired', 'status': done['status']}

    answer = ensure_tool(name, on_line=on_line, timeout=timeout)
    return {'ok': answer['ok'],
            'action': 'installed' if answer.get('installed') else '',
            'status': answer.get('status', '')}


def format_health(rows: Sequence[ToolHealth]) -> str:
    """The rows, and a count that says whether anything needs doing."""
    if not rows:
        return 'nothing checked'
    marks = {'ok': 'ok  ', 'warn': 'warn', 'fail': 'FAIL', 'absent': '--  '}
    width = max(len(r.label or r.name) for r in rows)
    out = []
    for row in rows:
        line = (f'{marks.get(row.level, "?   ")} '
                f'{(row.label or row.name):<{width}}  '
                f'{row.version or "":<8} {row.path or ""}')
        if row.why:
            line += f'\n     {row.why}'
            if row.fix:
                line += f'\n     -> {row.fix}'
        out.append(line.rstrip())
    counts = {}
    for row in rows:
        counts[row.level] = counts.get(row.level, 0) + 1
    out.append(', '.join(f'{counts[k]} {k}' for k in
                         ('ok', 'warn', 'fail', 'absent') if k in counts))
    return '\n'.join(out)
