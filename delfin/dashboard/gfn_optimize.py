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
from typing import Any, Callable, Dict, List, Optional

from . import solvents as _solvents

__all__ = ['GFN_METHODS', 'as_pushes', 'atom_lines', 'bond_graph',
           'bond_order_between', 'bond_order_note', 'read_bond_orders',
           'read_charges', 'not_a_stationary_point', 'rms_gradient',
           'constraint_input', 'contacts_holding', 'graph_changed',
           'bonds_to_freeze', 'graph_holds', 'method_is_out_of_its_depth',
           'a_rate_apart', 'paths_disagree', 'where_a_walk_jumped',
           'what_else_moved', 'pair_named',
           'gfnff_refusal', 'gfnff_pair_refusal', 'gfnff_would_form',
           'FOD_METHODS', 'can_measure_fod', 'fod_moved',
           'push_constant', 'turn_for',
           'restraint_energy', 'walk_the_path', 'lowest_real_modes',
           'find_xtb', 'find_binary', 'find_gxtb',
           'closest_contact', 'held_note', 'hold_atoms_at', 'settle_onto', 'install_command', 'install_root',
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
#: What xtb prints when it has been given a Hessian to work from.  The msRRHO
#: free energy at whatever temperature the $thermo block asked for.
_FREE_ENERGY_RE = re.compile(r'TOTAL FREE ENERGY\s+(-?\d+\.\d+)')
#: The rest of the same block, which was being thrown away.
#:
#: A Hessian is what a free energy costs, and once it has been paid for the
#: enthalpy and the zero-point energy are printed beside G in the same summary
#: -- reading three lines instead of one is free.  The entropy is not printed
#: as a total anywhere useful (xtb writes it per contribution, in cal/K/mol,
#: in a table whose columns move), so it is taken as ``T*S = H - G``: measured
#: on a methane at 298.15 K, H - G is 0.021114216337 Eh against the
#: 0.211142E-01 the table prints for T*S -- the same number to every digit
#: xtb gives.
_ENTHALPY_RE = re.compile(r'TOTAL ENTHALPY\s+(-?\d+\.\d+)')
_ZPE_RE = re.compile(r'zero point energy\s+(-?\d+\.\d+)')
#: How steep it is where the answer was taken, in Hartree per Bohr.
#:
#: This is what says whether a Hessian describes a stationary point or a point
#: on a slope.  Printed by every run, and never read until a Hessian could be
#: asked for on a structure somebody had dragged into shape -- which is
#: precisely the geometry where it is not safe to assume.
_GRADIENT_RE = re.compile(r'GRADIENT NORM\s+(-?\d+\.\d+)')
#: How many modes go the wrong way, which is what says whether a structure is
#: a minimum, a transition state, or neither.  xtb counts them itself, against
#: a cutoff of its own, so this is read rather than worked out again.
#:
#: The cutoff is printed beside the count and is -20 cm-1 in 6.7.1, so a mode
#: at -11 is listed among the frequencies and not among the imaginary ones.
#: That is why both are handed on: the count is the engine's own answer, and
#: the modes are what it was counting.
_IMAGINARY_RE = re.compile(r'#\s*imaginary freq\.\s+(\d+)')
#: How many electrons are not in a closed shell, which is the frontier gap's
#: question asked so that it has a number rather than a proxy.
#:
#: ``--fod`` runs one more SCF at 5000 K and adds up how far the occupations
#: are from two-and-zero; Grimme and Hansen, Angew. Chem. Int. Ed. 54 (2015)
#: 12308, https://doi.org/10.1002/anie.201501887.  Measured here under GFN2,
#: beside the gap the same run prints:
#:
#:     ethane C-C at 1.54 A      N_FOD 0.000   gap 15.08 eV
#:     ethane C-C at 3.50 A      N_FOD 1.726   gap  0.24 eV
#:     ethylene twisted 45 deg   N_FOD 0.075   gap  3.42 eV
#:     ethylene twisted 60 deg   N_FOD 0.251   gap  2.33 eV
#:     ethylene twisted 90 deg   N_FOD 2.000   gap  0.00 eV
#:
#: The gap and this agree about where the trouble is -- 60 degrees is where
#: :func:`method_is_out_of_its_depth` fires and where N_FOD leaves zero -- so
#: this is not a second opinion but the same one with a scale on it.
_NFOD_RE = re.compile(r'^\s*NFOD\s*:\s*(-?\d+\.\d+)', re.MULTILINE)
#: And where it sits.  xtb prints a Loewdin table under the total, one row an
#: atom, with the index and the element run together::
#:
#:      Loewdin FODpop      n(s)   n(p)   n(d)
#:          1C     0.8421   0.085  0.757  0.000
#:          8Mn    1.5533   0.001  0.002  1.550
#:
#: Which is the whole of "where is the trouble" for nothing: on the user's
#: 57-atom manganese complex, 1.553 of a total of 2.603 is on the manganese
#: and 1.550 of that is d.
#:
#: xtb writes a ``fod.cub`` beside it saying the same thing as a surface, and
#: that is read here and thrown away with the scratch directory.  It was
#: driven before it was dropped, so the next person need not: the cube is an
#: ordinary Gaussian one, and 3Dmol -- through
#: :func:`~delfin.dashboard.molecule_viewer.fukui_cube_isosurface_js`, which
#: already paints one -- parses and draws it.  Measured in chromium against
#: the 3Dmol this dashboard ships, at an isovalue of 0.005: an ethane pulled
#: to 3.50 A gives 384 vertices in 60 ms, an ethylene turned 90 degrees 536 in
#: 87 ms, and the manganese complex 572 in 301 ms.  Two things stand in the
#: way of it being a picture in the editor rather than a sentence.  The
#: manganese cube is 3.3 MB of text -- against 234 KB for the ethane -- and
#: every byte would have to cross the widget channel; and the editor's viewer
#: is the live one driven by ``window.__delfinSubmitManip``, which has no
#: channel for volumetric data at all, so it is a change to that script and
#: not a call to an existing function.  Both are worth doing and neither is
#: free.  (3Dmol's cube reader also takes one element too many off the end --
#: 15601 points for a 26x25x24 grid, the last of them NaN, from the file's
#: trailing newline.  Marching cubes indexes by grid position and never reads
#: it, so the picture is right; anything taking a minimum or a maximum over
#: ``VolumeData.data`` gets NaN.)
_FODPOP_HEADER = 'Loewdin FODpop'
_FODPOP_RE = re.compile(r'^\s*(\d+)([A-Za-z]{1,2})\s+(-?\d+\.\d+)')
#: How far apart the frontier orbitals are, which is how a single-determinant
#: method says whether it is still able to answer.
#:
#: It closes where two electrons stop being a pair.  Measured under GFN2 on a
#: 2-butene: 5.28 eV at the cis minimum and 2.42 at the twisted transition
#: state, where the pi bond is broken and the right description is two singly
#: occupied orbitals -- which a closed shell cannot be.  The number that comes
#: back there is not wrong by a little: GFN2 puts the twisted structure 49
#: kcal/mol up as a singlet, and there is a spurious minimum 64 above cis
#: where real trans lies 1.5 below.
#:
#: Nothing in this file can fix that; it is what the method is.  What it can
#: do is notice, because the gap says so.
#:
#: Both spellings, and the last of them.
#:
#: xtb writes the gap twice in two different hands: a summary block that says
#: ``:: HOMO-LUMO gap`` after every single point, and a property block that
#: says ``| HOMO-LUMO GAP`` at the end.  Only the second was matched, and
#: g-xTB does not print it -- measured on a propane, ``--sp``, ``--opt``,
#: ``--ohess`` and ``--grad`` alike, g-xTB prints the summary line and no
#: property block at all.  So the gap came back ``None`` under the one method
#: in the list a barrier is most likely to be quoted from, and
#: :func:`method_is_out_of_its_depth` -- which returns nothing at all for a
#: gap it cannot read -- said nothing about a structure it should have warned
#: about.
#:
#: The *last* one, because the summary block is printed at every geometry the
#: run passes through and the first of them describes the structure that went
#: in rather than the one that came out: on the same propane under GFN2,
#: ``--ohess`` prints 14.3656 for the input geometry and 14.5837 for the
#: relaxed one.  Taking the last reproduces what this matched before, to every
#: digit, for GFN2, GFN1 and GFN-FF across all four run types -- and gives
#: g-xTB the number it was missing.
_GAP_RE = re.compile(r'HOMO-LUMO\s+gap\s+(-?\d+\.\d+)', re.IGNORECASE)
_EIGVAL_RE = re.compile(r'eigval\s*:\s*(.+)')
_VERSION_RE = re.compile(r'xtb version\s+([0-9.]+)')
# What the run says it did, taken from its own output rather than from the
# flags we passed: an xtb that ignored a flag would otherwise be indisting-
# uishable from one that honoured it.
_HAMILTONIAN_RE = re.compile(r'Hamiltonian\s+(GFN[0-9A-Za-z-]*)')
_GFNFF_BANNER = 'G F N - F F'
#: What GFN-FF leaves behind that describes the molecule rather than the run.
_GFNFF_TOPOLOGY_FILES = ('gfnff_topo', 'gfnff_charges')


def _restart_named(topology, key, charge, uhf, solvent, model):
    """Where this run's wavefunction is kept between calls, or None.

    A drag is a sequence of xtb runs on geometries a fraction of an Angstrom
    apart, and every one of them starts its SCF from the extended Hueckel guess
    as though it had never seen the molecule.  xtb writes the converged
    wavefunction to ``xtbrestart`` and reads it back if it finds one, which is
    a better guess by far.

    How much it is worth depends on the shape of the run, and the difference
    is not small.  Measured on a 57-atom manganese complex under GFN2 at
    charge +1, cold and warm alternately so a machine carrying other people's
    work could not favour either, six of each:

        single point        0.42 s cold, 0.17 s warm   (59 % less)
        20-cycle relaxation 2.65 s cold, 2.46 s warm   ( 7 % less)

    A single point is one SCF, so the whole of it is the guess.  A twenty-cycle
    relaxation is twenty SCFs and xtb already restarts from the last one inside
    the run, so only the first of the twenty can be helped -- which is what a
    seventh of the cost is.  The drag's own answers are relaxations; its
    anchor, and the price of a placing hand, are single points.

    It changes no number.  Ten points along a Mn-O stretch priced both ways
    agree to 2.1e-5 kcal/mol at worst -- an SCF converged to its own threshold
    from two different starting guesses -- against a budget that refuses on
    tenths.  A six-answer walk priced three times cold and three times warm
    has the warm prices inside the cold spread, which is about 0.2 kcal/mol on
    that walk: xtb is not bit-reproducible under threading, and this sits
    under that.

    Keyed by everything that makes the wavefunction a different object.  A
    restart written for one charge, spin, method or solvent is not a guess for
    another, it is the wrong number of electrons in the wrong field; xtb would
    take it and the answer would be confidently wrong.  The directory is the
    caller's and already belongs to one molecule -- see the editor's
    ``_gfn_topology_dir``, which keys it on the element column -- so what is
    left to name here is the run.

    GFN-FF has no wavefunction and writes none, so it gets nothing.
    """
    if topology is None or str(key) == 'gfnff':
        return None
    wet = str(solvent or '').strip().lower() or 'gas'
    how = str(model or '').strip().lower() or 'alpb'
    return Path(topology) / (
        f'xtbrestart-{key}-q{int(charge)}-u{max(0, int(uhf))}-{how}-{wet}')


def _keep_restart(folder: Path, kept: Path, output: str) -> None:
    """Put this run's wavefunction where the next one will find it.

    Only from a run that finished.  A wavefunction from a run that stopped
    with an error is not a guess, and the next answer would start from it and
    have no way of knowing; so a bad run takes the stored one with it and the
    answer after that starts cold, which is how it worked before any of this.

    Written beside and then renamed, because the directory belongs to the
    molecule rather than to one run: an optimisation and a drag can be reading
    and writing it at the same moment, and a copy is not one operation.  Half a
    wavefunction handed to the next run is the one way this could produce a
    wrong answer instead of merely a slower one, and a rename is atomic.
    """
    made = folder / 'xtbrestart'
    finished = ('abnormal termination' not in output
                and '[ERROR]' not in output and made.is_file())
    try:
        if finished:
            kept.parent.mkdir(parents=True, exist_ok=True)
            beside = kept.with_name(kept.name + f'.{os.getpid()}.part')
            shutil.copy2(str(made), str(beside))
            os.replace(str(beside), str(kept))
        elif kept.is_file():
            kept.unlink()
    except OSError:
        pass


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


#: Below this a printed frequency is one of the six xtb projects out.
#:
#: The translations and the rotations come out of the Hessian as zero and are
#: printed as ``-0.00`` and ``0.00`` exactly, so a cut this fine separates
#: them from anything that is a vibration -- including a very soft one, which
#: is the case that matters here: a real mode of a few wavenumbers is what
#: says a stationary point is nearly flat in that direction, and rounding it
#: away with the trivial modes would throw out the very thing being looked
#: for.
_TRIVIAL_MODE = 0.005


def lowest_real_modes(output: Any, most: int = 4) -> List[float]:
    """The softest genuine vibrations of the last Hessian in *output*.

    The last one, because a run that took several -- a saddle search
    recalculates as it goes -- has left the geometries the earlier ones
    describe.

    Real rather than imaginary, because how flat a stationary point is in its
    remaining directions is a separate question from how many of them go the
    wrong way, and one that published searches ask: AutoMeKin, with its
    default ``tight_ts``, refuses a transition state when the two lowest real
    frequencies sum to less than 10 cm-1.  A structure that flat is on the
    edge of being a saddle point of higher order, and the number that says so
    is here rather than among the imaginary modes.
    """
    text = str(output or '')
    block = text.split('projected vibrational frequencies')[-1]
    seen: List[float] = []
    for row in _EIGVAL_RE.finditer(block):
        for word in row.group(1).split():
            try:
                seen.append(float(word))
            except ValueError:
                pass
        # The list is printed lowest first, so the soft end is at the front
        # and a big molecule's three hundred modes need not all be read.
        if len(seen) > 200:
            break
    return sorted(one for one in seen if one > _TRIVIAL_MODE)[:max(0, int(most))]


def _read_fod(output: Any) -> Optional[Dict[str, Any]]:
    """The fractional occupation this run printed, or None if it printed none.

    None rather than zero, and the difference matters: a molecule with no
    static correlation prints ``NFOD : 0.0000``, and a method that cannot be
    asked prints nothing at all.  Both would read as "there is nothing wrong
    here" if they came back the same, and one of them is a question that was
    never put -- see :data:`FOD_METHODS` for the two methods that do this.

    The per-atom breakdown is read out of the Loewdin table under the total,
    as a block bounded by its header and the blank line after it rather than
    by a pattern over the whole output: the row shape is ``8Mn 1.5533 ...``,
    which is loose enough to match other tables xtb prints further down.
    """
    text = str(output or '')
    total = _NFOD_RE.search(text)
    if not total:
        return None
    try:
        amount = float(total.group(1))
    except ValueError:
        return None
    on: List[tuple] = []
    lines = text.splitlines()
    for n, line in enumerate(lines):
        if _FODPOP_HEADER not in line:
            continue
        for row in lines[n + 1:]:
            if not row.strip():
                break
            found = _FODPOP_RE.match(row)
            if not found:
                break
            try:
                # xtb counts its atoms from one and this editor counts from
                # zero, and the two are within one line of each other
                # everywhere this is reported.
                on.append((found.group(2), int(found.group(1)) - 1,
                           float(found.group(3))))
            except ValueError:
                break
        break
    return {'total': amount, 'on': on}


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

#: What one xtb force constant is worth, in the units a chemist reads.
#:
#: xtb takes its distance force constants in Hartree per Bohr squared, so a
#: number like 0.5 says nothing about how hard it pulls until it is converted:
#: 0.5 is 1120 kcal/mol/A^2, half again as stiff as a C-C bond.  Which is why
#: a "pull" moves an ethane's C-C from 1.52 to 2.46 A when it is asked for
#: 2.50, measured under GFN2 -- it is not a gentle setting at all, it is
#: simply less than exact.
EH_PER_BOHR2_IN_KCAL = 627.5095 / (0.529177210903 ** 2)
BOHR_IN_ANGSTROM = 0.529177210903

#: What xtb's restraint actually costs, measured rather than assumed.
#:
#: ``E = k * d^2``, not the half of that a spring is usually written with,
#: with *d* in Bohr for a distance and in radians for an angle.  A torsion is
#: not that shape at all: it is ``k * (1 - cos d)``, which is what a periodic
#: coordinate has to be -- a harmonic in the angle would have a step in it
#: where the angle wraps.
#:
#: Measured under GFN2 on an ethane, a propane and a butane by differencing
#: the reported total against a single point on the same geometry.  The ratio
#: against ``0.5 * k * d^2`` came out at 2.00 for every one of four force
#: constants across two decades; an angle's residue matched ``k * rad^2`` to
#: four figures while ``k * deg^2`` was out by three thousand; and a butane
#: torsion held 120 degrees from its target cost 4.7063 kcal/mol where
#: ``k * rad^2`` says 13.76 and ``k * (1 - cos)`` says 4.7063.
#:
#: Two things follow.  The force a constant applies is ``2 * k * d``, so a
#: number offered to the user as kcal/mol/A has to be halved on the way in --
#: see :func:`push_constant`.  And the restraint's own energy can be
#: subtracted from the answer instead of asking for a second calculation
#: without it, which is one xtb process saved on every step of every drag.
RESTRAINT_IS_K_TIMES_D_SQUARED = True

#: How far ahead a push holds its target, in Angstrom.
#:
#: A harmonic restraint left at a fixed far target pulls harder the further it
#: is from being met, so it is a spring and not a force.  Moved up with the
#: structure each step -- always this far ahead of where the coordinate now
#: stands -- what the molecule feels is a constant k*reach, which is what an
#: artificial force is.  See :func:`push_constant`.
PUSH_REACH = 1.0

#: The ceiling force a push starts and ends with, in kcal/mol/A.
#:
#: The same xtb force constants throughout -- what changed when the factor of
#: two in the restraint was measured is what these numbers are *called*, not
#: what they do.
#:
#: A ceiling and not a fixed force: the restraint applies ``k * reach`` only
#: while the coordinate has not moved at all, and what the structure keeps is
#: whatever of it the structure can answer.  Which is the property the whole
#: thing rests on.  Measured under GFN2 on an ethane at rest at 1.5212 A, the
#: same push applied over and over:
#:
#:     30 kcal/mol/A   1.6270, 1.6413, 1.6413, 1.6413, 1.6413, 1.6413
#:     60 kcal/mol/A   1.7613, 1.8815, 1.9736, 2.0686, 2.2524, 2.8766
#:
#: Thirty settles and stays, however long it is pushed.  Sixty never settles:
#: past the inflection of the curve the restoring force falls away, so once
#: the push is stronger than the most the bond can pull back with, nothing
#: balances it.  So a force below a bond deforms and a force above it breaks
#: -- and which of the two happens is a number the user set.
#:
#: The range is AFIR's own, read through its collision parameter: gamma of 100
#: to 400 kJ/mol works out at roughly 25 to 100 kcal/mol/A of model force.  It
#: starts below anything that can break a C-C and ends above it, because
#: finding a reaction is raising the force until something gives.
PUSH_FORCE_FROM = 8.0
PUSH_FORCE_TO = 240.0

#: What a bond holds, as a force, in kcal/mol/A.
#:
#: When a frontier gap says a single-determinant answer has stopped being
#: one: how far it has to fall, and how small it has to get.
#:
#: The *fall* is the signal, not the value.  A great many perfectly ordinary
#: systems simply have a small gap -- transition-metal complexes, long
#: conjugated chains, anything with close-lying frontier orbitals -- and they
#: are not in trouble for it; a warning keyed on the value alone would fire on
#: every one of them all the time, which is noise and teaches people to
#: ignore it.  What means something is a gap that *closes on the way*: two
#: electrons that were a pair and are not one any more.
#:
#: So: fallen to about half, and low enough afterwards to matter.  Or very low
#: on its own, where there is nothing to compare against.  Measured under GFN2
#: driving a torsion through a double bond, 5.26 eV at the start and 1.44 at
#: the top; and on the transition state a path finder estimated, 5.28 and
#: 2.42.
GAP_IS_SMALL = 1.5
GAP_WORTH_SAYING = 3.0
GAP_HAS_FALLEN = 0.55


def method_is_out_of_its_depth(gap: Any, was: Any = None) -> str:
    """What to say when the frontier gap says the answer is not one.

    Empty when there is nothing to say.  A closed-shell single determinant
    describes two electrons in one orbital; where the gap closes they are not
    in one any more, and the energy that comes back is about a state the
    method cannot represent.  It cannot be fixed here -- it is what the method
    is -- but it can be said, and a barrier quoted without it is a number
    somebody will use.
    """
    try:
        gap = float(gap)
    except (TypeError, ValueError):
        return ''
    fell = False
    try:
        fell = (float(was) > 0 and gap < GAP_HAS_FALLEN * float(was)
                and gap < GAP_WORTH_SAYING)
    except (TypeError, ValueError):
        fell = False
    if gap >= GAP_IS_SMALL and not fell:
        return ''
    said = (f'The frontier gap is {gap:.1f} eV there'
            + (f', down from {float(was):.1f}' if fell else '')
            + '. A closed-shell method describes two electrons in one orbital, '
              'and where the gap closes they are not in one -- so this barrier '
              'is the method at the edge of what it can represent, and is '
              'worth checking open-shell.')
    return said


#: The methods that can be asked how much of the structure is not a closed
#: shell, which is not the same list as the methods that can be run.
#:
#: Both exclusions are measured, and they fail in opposite ways.  GFN-FF has
#: no wavefunction at all: given ``--fod`` it prints ``[WARNING] Runtime
#: exception occurred``, no NFOD and no cube, and then says ``normal
#: termination``.  g-xTB is worse -- it takes the flag, converges, terminates
#: normally, and prints no NFOD and writes no cube, so a caller that only
#: looked at the exit code would report a molecule with no static correlation
#: rather than a question that was never asked.  Silence that reads as good
#: news is the one failure this module exists to prevent, so the list is
#: positive: the two methods that were seen to answer.
FOD_METHODS = ('gfn1', 'gfn2')


def can_measure_fod(method: Any) -> bool:
    """Whether *method* answers ``--fod`` -- see :data:`FOD_METHODS`."""
    return str(method or '').strip().lower() in FOD_METHODS


#: How much of a change in N_FOD is worth a sentence.
#:
#: Set against the walk that calibrates it rather than against a threshold
#: from a paper.  Turning ethylene's double bond, GFN2: 0.008 flat, 0.075 at
#: 45 degrees, 0.251 at 60 -- which is where the frontier-gap rule already
#: fires -- and 2.000 at 90.  A tenth of an electron is the point where the
#: two rules start saying the same thing, so it is where this starts to speak.
FOD_HAS_MOVED = 0.10


def fod_moved(first: Any, later: Any, on: Any = ()) -> str:
    """What the fractional occupation did along a walk, or nothing to say.

    A *change*, never a number on its own, and that is not a stylistic choice.
    Measured against ORCA's own FOD at TPSS/def2-TZVP: benzene 0.032 against
    0.016 and ozone 0.426 against 0.439 -- three percent on the textbook
    diradicaloid -- but Ni(CO)4 0.335 against 0.085, a four-fold over-report
    on a closed-shell metal carbonyl.  A fixed cutoff would therefore tell
    somebody working on metal complexes that every one of their structures is
    a diradical.  What survives that is the difference between two points of
    one walk on one molecule, which is the same discipline
    :func:`method_is_out_of_its_depth` applies to the gap.

    Measured on a real one rather than only on a benchmark: a 57-atom
    manganese complex at charge +1, closed shell, reads 2.65 before anything
    has been moved -- a number an absolute rule would call a triradical, with
    1.55 of it on the metal and nearly all of that d.  Walking its Mn-Br bond
    out moves it to 3.42, and the +0.77 is the finding.  The 2.65 is what a
    threshold would have reported instead.

    *on* is the per-atom breakdown at the later point, as
    ``[(symbol, index, value), ...]``, and is used only to name where it went.
    """
    try:
        first = float(first)
        later = float(later)
    except (TypeError, ValueError):
        return ''
    if later - first < FOD_HAS_MOVED:
        return ''
    # Named only where the claim is true.  A bond broken symmetrically puts
    # half on each of two atoms -- measured on the ethane homolysis, 0.842 of
    # 1.726 on each carbon -- and "most of it on C0" about a 50/50 split is a
    # sentence that says something false about which end of the bond is the
    # trouble.  More than half is the bar, which is what tells that case apart
    # from the one where naming the atom *is* the finding: on the manganese
    # complex 1.553 of 2.603 sits on the metal, and 1.550 of that is d.
    busiest = ''
    ranked = sorted(
        [one for one in (on or ()) if len(one) >= 3],
        key=lambda one: -float(one[2]))
    if ranked and float(ranked[0][2]) > 0.5 * later:
        symbol, index = ranked[0][0], int(ranked[0][1])
        busiest = f', most of it on {symbol}{index}'
    return (
        f'Electrons in half-filled orbitals went from {first:.2f} to '
        f'{later:.2f} across the walk{busiest}. Near zero, one arrangement of '
        f'the electrons is the structure; near two it takes at least two, and '
        f'this method has one. Read the change, not the number: the same '
        f'measurement runs four times high on a closed-shell metal complex, '
        f'so no value of it means "bad" on its own.')


#: A force is a slope, so this is the whole meaning of the setting: the hand
#: climbs any part of the surface less steep than it is, and stops where the
#: surface gets steeper.  Which is why "as hard as a bond" is the right
#: yardstick -- a bond that holds at 110 kcal/mol/A is a wall 110 kcal/mol/A
#: steep, and a hand set below that walks up to it and no further.
#:
#: The yardstick the hand is set against, and it is worth having because it
#: turns out to be nearly one number.  Bisected under GFN2, the same push
#: applied over and over until the coordinate either settles or runs away:
#:
#:     C-C in ethane     holds 112, breaks by 116
#:     C-H in ethane      holds 98, breaks by 100
#:     C-O in methanol   holds 120, breaks by 122
#:
#: Three bonds of quite different length and strength, all of them giving way
#: between 98 and 122.  So "as hard as a bond" is a real quantity and not a
#: figure of speech, and a hand set at half of it deforms a molecule without
#: being able to take it apart -- whatever the molecule is made of.
A_BOND_HOLDS = 110.0

#: The same reach, in the units xtb measures an angle in.
#:
#: xtb takes one force constant for a whole block, so an angle in the same
#: block as a distance gets the distance's constant -- and its restraint is
#: ``k * radians^2`` where the distance's is ``k * Bohr^2``.  One angstrom is
#: 1.8897 Bohr, so the reach that gives an angle the same energy as a distance
#: is 1.8897 radians: 108.3 degrees.
#:
#: A round number was used instead, and twenty degrees of it.  That gave a
#: torsion 0.38 kcal/mol to spend where the distance had 11.15 -- a thirtieth
#: -- and a rotational barrier is about three.  So the hand could stretch
#: bonds and could not turn anything, and a molecule could not be put into its
#: own conformers at any temperature.  Which is the one thing room temperature
#: certainly does allow.
PUSH_REACH_DEGREES = math.degrees(PUSH_REACH / BOHR_IN_ANGSTROM)


def push_energy(force: float, reach: float = PUSH_REACH) -> float:
    """The most work a push of *force* can do over its reach, in kcal/mol.

    Half of force times reach, because a restraint is a spring: it reaches
    *force* only at full stretch and is weaker all the way there.  Which is
    why a hand set to the thermal ceiling divided by the reach could only ever
    spend half of what the temperature allowed.
    """
    return 0.5 * float(force) * (float(reach) or PUSH_REACH)


def push_force_for(energy: float, reach: float = PUSH_REACH) -> float:
    """The force whose push can spend *energy* kcal/mol over its reach.

    The inverse of :func:`push_energy`, and the one the thermal budget wants:
    what the temperature grants is an energy, and what a hand is set in is a
    force.

    The whole chain from a temperature to an xtb force constant, so that where
    it stops being derived is visible rather than buried:

    1. Eyring, inverted.  The largest barrier crossable at *T* within a time
       *tau* is ``R T ln(kB T tau / h)`` -- 22.3 kcal/mol at 298.15 K within
       the hour, 10.9 at 150 K, 117 at 1500 K.  See ``thermal_ceiling``.
       Nothing is chosen here.

       That ceiling is a free energy and what is priced against it is an
       electronic one, which is an approximation and is meant as one.  A free
       energy is a Hessian, and a Hessian is not free: measured under GFN2,
       0.57 s against 0.29 for sixteen atoms and 3.72 against 0.76 for
       twenty-four.  A drag answers ten times a second, so there is no version
       of this that has a free energy in it -- and an entropy term that is
       roughly constant along a path largely cancels in a difference, which is
       why the approximation is a reasonable one rather than merely a cheap
       one.
    2. That energy into a force.  A restraint is a spring: it reaches its
       force only at full stretch and is weaker all the way there, so what it
       can spend over a displacement is half of force times displacement.
       ``F = 2 * ceiling / reach``.  Nothing is chosen here either, *given* a
       reach.
    3. The force into a constant, per coordinate kind, from the shapes xtb
       actually uses -- measured, not assumed: ``k d^2`` in Bohr for a
       distance, the same in radians for an angle, ``k (1 - cos d)`` for a
       torsion.  See :func:`push_constant`.  Nothing is chosen here.

    The reach is chosen.  It is the length over which the barrier is taken to
    be climbed, and no temperature says what that should be -- it is a
    chemical scale: a bond being made or broken has its top something like
    half an angstrom to an angstrom from the minimum, so an angstrom is used,
    and the angular reach is that same length in the units xtb measures an
    angle in.  It is worth knowing how much rests on it: at four tenths of an
    angstrom -- the inflection of a Morse bond, which is another defensible
    reading -- room temperature would come out at 111 kcal/mol/A, which is
    exactly what a bond holds, and room temperature would break bonds.  An
    angstrom gives 45, four tenths of a bond, and that has been checked
    against real cases at three temperatures: it turns a molecule into its own
    conformers at 150 K and holds an aryl C-Br at 298 K.

    A hand that can *spend* a barrier's worth is not automatically a hand that
    can *drive* it -- a spring pulling a system over a hill has to be steeper
    than the hill at its steepest, which for a three-fold torsion of height B
    is 1.5 B per radian.  At 298 K this comes out at 45 kcal/mol per Angstrom
    and per radian against a 4.8 kcal/mol rotational barrier that needs about
    11, so there is room; at 150 K it is 22 against the same 11, and there
    still is.  What stops what the temperature forbids is the wall, which
    prices the structure that was actually reached -- not the hand, which only
    has to be able to ask.
    """
    return 2.0 * float(energy) / (float(reach) or PUSH_REACH)


def as_pushes(entries: Any, reference: Any, force: float,
              value_of: Any = None, reach: Any = PUSH_REACH,
              most: Any = None) -> list:
    """The same held coordinates, as forces instead of as values.

    *entries* are what :func:`contacts_holding` worked out the hand is
    changing, carrying the values the hand is *asking* for.  Each becomes a
    push of at most *force* kcal/mol/A towards that value, held no further
    than a reach from where the coordinate stands in *reference* -- so what
    the structure feels is a force with a ceiling, and how far it actually
    goes is the structure's answer.

    *value_of* reads one entry's coordinate off a geometry; the editor passes
    its own, which is the only place that arithmetic lives.

    The hand gets stronger the further it is dragged, which is what pulling
    on something is like -- but what is *asked for* stays within the reach.
    Those are two different things and conflating them cost both ways round.
    Capped, the hand could never do more than the slider was already set for
    however far it was dragged.  Uncapped, the target ran arbitrarily far
    ahead of the structure and a few cycles of a strong restraint overshot it
    and came back, which on screen is a molecule that shakes.

    So the target is held at the reach and the force constant carries the
    excess: at twice the reach the push is twice as strong, and the force at
    the target is exactly what an unclamped spring would have applied there.
    Same physics, nothing to overshoot.

    *most* is the largest force the hand may reach whatever the drag does,
    which is what a temperature is; ``None`` is no limit at all.
    """
    # One force constant for the whole block, sized for the coordinate the
    # hand is actually driving.  A turn is a torsion, and the distances that
    # come with it are anchors -- they are asked for where they already are,
    # so a constant that would be fierce for a drag costs them nothing and
    # holds the fragment together while the torsion negotiates.
    kinds = [str((one or {}).get('kind') or '') for one in (entries or ())]
    leads = ('dihedral' if 'dihedral' in kinds
             else 'angle' if 'angle' in kinds else 'distance')

    def _span(kind, far):
        return far if kind == 'distance' else math.degrees(far
                                                           / BOHR_IN_ANGSTROM)

    def _apart(kind, wanted, stands):
        """How far apart two values of this coordinate are.

        A torsion is periodic and 350 degrees away is ten degrees away.
        Subtracted plainly -- as this did -- a hand near the far side of a
        turn reads as being most of a circle out, the strength that carries
        the excess goes up by a factor of thirty, and what comes back is
        whatever an optimiser does with an absurd restraint.  Measured on a
        2,4-hexadiene, the same drag at three settings: 0.268 A at four
        tenths of a bond, 0.089 at one whole one and 0.111 at three.  A
        stronger hand moving an atom less is not physics.
        """
        gap = float(wanted) - float(stands)
        if kind == 'dihedral':
            gap = (gap + 180.0) % 360.0 - 180.0
        return gap

    # How far the hand has run ahead of the structure, in reaches.  Read off
    # the coordinate the hand is driving, since that is the one it is pulling.
    far = PUSH_REACH if reach is None else float(reach)
    over = 1.0
    if far > 0 and reference is not None and value_of is not None:
        for entry in (entries or ()):
            if str((entry or {}).get('kind') or '') != leads:
                continue
            try:
                stands = value_of(reference, entry)
                gap = abs(_apart(leads, entry.get('value'), stands))
            except (TypeError, ValueError, Exception):
                continue
            over = max(over, gap / _span(leads, far))
    pulling = float(force) * over
    if most is not None:
        pulling = min(pulling, float(most))
    constant = push_constant(pulling, kind=leads)
    out = []
    for entry in (entries or ()):
        kind = str((entry or {}).get('kind') or '')
        wanted = (entry or {}).get('value')
        try:
            wanted = float(wanted)
        except (TypeError, ValueError):
            continue
        span = _span(kind, far)
        now = None
        if reference is not None and value_of is not None and span > 0:
            try:
                now = value_of(reference, entry)
            except Exception:
                now = None
        if now is not None:
            gap = _apart(kind, wanted, now)
            if abs(gap) > span:
                gap = math.copysign(span, gap)
            wanted = now + gap
        out.append({'kind': kind,
                    'atoms': [int(i) for i in (entry.get('atoms') or ())],
                    'mode': 'push',
                    'force': constant,
                    'value': wanted})
    return out


def push_pulls_hardest(kind: str, reach: float = PUSH_REACH) -> float:
    """The hardest one unit of force constant can pull, over the reach.

    In kcal/mol per Angstrom for a distance and per radian for an angle or a
    torsion, which is the same statement in each coordinate's own units.  The
    three store energy in three different shapes -- ``k * d^2`` in Bohr, in
    radians, and ``k * (1 - cos d)`` -- so one number means three different
    strengths, and xtb takes one number for the whole block.  This is where
    the three are written down and nowhere else.
    """
    span = (float(reach) or PUSH_REACH) / BOHR_IN_ANGSTROM
    if str(kind) == 'dihedral':
        # k * (1 - cos d): the torque is k * sin d, hardest at a right angle.
        return 627.5095
    if str(kind) == 'angle':
        return 2.0 * 627.5095 * span
    # A distance is asked for in Angstrom and answered in Bohr.
    return 2.0 * EH_PER_BOHR2_IN_KCAL * (float(reach) or PUSH_REACH)


def push_pulls_now(kind: str, constant: float, gap: float) -> float:
    """What a push of this force constant is applying at this moment.

    The companion to :func:`push_pulls_hardest`, which says what the same
    constant would apply at full stretch.  A restraint is a spring: it reaches
    its ceiling only when the target is a whole reach away and is weaker all
    the way there, so the two numbers are different by however far the hand
    has run ahead of the structure -- which is the whole of what a mouse
    controls while a drag is running.

    They were the same number on screen, and that was a misreport rather than
    a rounding: the status line quoted the ceiling. Measured on the user's
    manganese complex, a phenolate oxygen dragged at 0.4 of a bond under
    GFN-FF, the hand ran about half a degree of torsion ahead of the structure
    per answer -- 0.4 kcal/mol per radian applied, reported as 44.

    *gap* is how far the target is from where the coordinate stands, in the
    coordinate's own units: Angstrom for a distance, degrees for an angle or a
    torsion.  The answer is in kcal/mol per Angstrom, or per radian.
    """
    far = abs(float(gap))
    if str(kind) == 'dihedral':
        # k * (1 - cos d), so the torque is k * sin d and falls away past a
        # right angle rather than growing.
        return 627.5095 * float(constant) * math.sin(math.radians(far))
    if str(kind) == 'angle':
        return 2.0 * 627.5095 * float(constant) * math.radians(far)
    return 2.0 * EH_PER_BOHR2_IN_KCAL * float(constant) * far


def push_constant(force: float, reach: float = PUSH_REACH,
                  kind: str = 'distance') -> float:
    """The xtb force constant whose ceiling is *force* kcal/mol/A.

    xtb's restraint is ``k * d^2``, so what it applies is ``2 * k * d`` -- the
    two is measured, not assumed, and leaving it out made every force in this
    file a claim about half of what was really being applied.

    Held *reach* from where the coordinate stands, the restraint applies at
    most ``2 * k * reach``, and less than that once the coordinate has moved
    towards it.  The caller puts the target *reach* ahead again at every step,
    so the ceiling is the same at every point of a path however far it has
    already come -- which is what makes a ramp of these a ramp of forces.

    *kind* is the coordinate the constant is meant for, and it matters a
    great deal: a torsion stores ``k * (1 - cos d)`` where a distance stores
    ``k * d^2``, so the same number pulls a fifth as hard on one as on the
    other.  Sized as a distance -- as it was -- the hand's largest torque fell
    just short of an ordinary rotational barrier, and a molecule could not be
    turned into its own conformers at any temperature.  Which is the one thing
    room temperature certainly does allow.

    Set here so that the hardest the hand can pull is *force*, in the
    coordinate's own units: kcal/mol per Angstrom, or per radian.

    Returned in xtb's units, which is the only place that conversion lives.
    """
    return float(force) / push_pulls_hardest(kind, reach)


def restraint_energy(xyz_text: str, constraints: Any, value_of: Any) -> Any:
    """What the restraints in *constraints* contribute to the reported energy.

    In Hartree, to be taken off the total xtb hands back.  ``None`` when it
    cannot be worked out, and then the caller has to ask for a calculation
    without them instead of guessing.

    A push is *meant* to leave a residue -- that residue is the force -- so an
    answer priced as it stands carries the hand as well as the structure.  The
    alternative was a second xtb process per step with the constraints taken
    out; this is the same number for none of the time.
    """
    prepared = constraint_input(constraints, atoms=_atom_count(xyz_text))
    k = prepared.get('force')
    if not k or not prepared.get('held'):
        return 0.0
    dropped = [id(one) for one in (prepared.get('dropped') or ())]
    total = 0.0
    for entry in (constraints or ()):
        if id(entry) in dropped:
            continue
        got = value_of(xyz_text, entry)
        if got is None:
            return None
        try:
            gap = float(entry.get('value')) - float(got)
        except (TypeError, ValueError):
            return None
        kind = str(entry.get('kind'))
        if kind == 'distance':
            total += k * (gap / BOHR_IN_ANGSTROM) ** 2
        elif kind == 'dihedral':
            # Periodic: 359 degrees away is one degree away, and the shape is
            # k * (1 - cos) rather than a harmonic -- measured.
            gap = (gap + 180.0) % 360.0 - 180.0
            total += k * (1.0 - math.cos(math.radians(gap)))
        else:
            total += k * math.radians(gap) ** 2
    return total


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
    exact = negotiated = dragging = pushing = False
    softest = None
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
        elif mode == 'push':
            # An artificial force rather than a value to be met: the entry
            # brings its own stiffness, because what the push is *for* is that
            # the structure decides where it ends up.
            pushing = True
            try:
                asked = float((entry or {}).get('force'))
            except (TypeError, ValueError):
                asked = push_constant(PUSH_FORCE_FROM)
            softest = asked if softest is None else min(softest, asked)
        else:
            negotiated = True
        kept.append((kind, indices, value))

    # A value under the hand sets the block's stiffness whatever else is in
    # it: it is the one being perturbed every frame, so it is the one that
    # decides whether the picture is still.  See DRAG_FORCE_CONSTANT.
    # A push is three orders of magnitude softer than a hold, and there is
    # one force constant for the block, so the two cannot both be had.  The
    # stiffer wins -- a value the user said to hold is held -- and ``mixed``
    # says so, which is what lets the caller refuse instead of quietly turning
    # a push into a scan.
    force = (DRAG_FORCE_CONSTANT if dragging
             else FIX_FORCE_CONSTANT if exact
             else PULL_FORCE_CONSTANT if negotiated
             else softest if pushing and softest is not None
             else PULL_FORCE_CONSTANT)
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
        'mixed': bool(exact and negotiated)
                 or bool(pushing and (exact or negotiated or dragging)),
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





#: How far a coordinate's change has to move atoms to count as one the hand is
#: driving, in Angstrom, and how small a share of the largest still counts.
#:
#: One currency for everything.  A bond, an angle and a torsion are scored by
#: how far their change moves atoms -- for a torsion that is the lever arm
#: times the angle in radians -- so they are directly comparable and there is
#: no threshold holding a degree against a length.  There is also no separate
#: branch for "the hand turned rather than pulled": a torsion is simply the
#: coordinate that scores highest when it did.
_MOVED_ANGSTROM = 0.02
_MOVED_SHARE = 0.25

#: How straight a contact has to be to the hand's own motion before it counts
#: as a bond being made.  A pair the hand is closing at nearly the full rate of
#: its own step is a reaction coordinate whatever the distance happens to be;
#: one it is sweeping past is not.
_ALONG_ENOUGH = 0.5


def _is_a_bond(where, radius, i, j, slack: float = 1.25) -> bool:
    """Whether these two are bonded, by covalent radii."""
    return math.dist(where[i], where[j]) < slack * (radius[i] + radius[j])


#: Where a bond starts and stops being one, as a share of the two radii.
#:
#: ORCA's GOAT calls a bond anything inside 1.3 times the sum of the covalent
#: radii, and by default keeps the bonds a structure came with -- only
#: ``GOAT-EXPLORE`` lets them break.  The same number is used here, because
#: this is the same question and a conformer search is where it has been
#: thought about hardest.
#:
#: One number, which is what the established procedures use.  GOAT keeps the
#: bonds a structure came with against that threshold; CREST takes the input
#: structure as its reference and discards what has changed against it.
#:
#: A second, looser one was tried here -- a bond that was there had to be
#: clearly gone before it counted as broken -- on the grounds that a bond
#: resting on the threshold would flicker from one answer to the next.  It
#: does not pay for itself: accepting at 1.4 and reporting at 1.3 lets a bond
#: stretch across several accepted steps and then hands back a geometry that
#: is already broken by the stricter measure.  Measured on a 2,4-hexadiene
#: with the rigid hand, five steps: three bonds broken without the wall, and
#: one broken *with* it.  Two measures and one wall is not a wall.
#:
#: Flicker is not a problem a threshold has to solve anyway.  A step whose
#: bonding differs is refused and the next one is tried; refusing one frame
#: of a drag that was not really changing costs a tenth of a second and
#: nothing else.
BOND_STARTS_AT = 1.3


def bond_graph(xyz_text: str, slack: float = BOND_STARTS_AT) -> frozenset:
    """Which atoms this geometry has bonded to which, by covalent radii.

    The same test the rest of this file contacts with, and the same one the
    viewer draws lines with, so what is compared is what is seen.  Returned as
    a frozen set of index pairs, which makes "did the molecule stay the same
    molecule" a set difference.

    Covalent radii and nothing else: a bond order would need a wavefunction,
    and this has to be answerable ten times a second on whatever the hand has
    just made -- including geometries no method would call a molecule.
    """
    rows = [line.split() for line in atom_lines(xyz_text)]
    if not rows:
        return frozenset()
    from delfin.atom_mapping import cov_radius

    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    radius = [cov_radius(str(r[0])) for r in rows]
    return frozenset(
        (i, j)
        for i in range(len(rows)) for j in range(i + 1, len(rows))
        if _is_a_bond(where, radius, i, j, slack=slack))


def graph_holds(before: Any, xyz_text: str) -> tuple:
    """Whether *xyz_text* still has the bonds *before* had, and what changed.

    Returns ``(holds, said)``.  One threshold, judged the same way for a bond
    that was there and one that was not -- see :data:`BOND_STARTS_AT` for why
    the two it had are worse than the one it has.
    """
    rows = [line.split() for line in atom_lines(xyz_text)]
    if not rows:
        return True, ''
    now = bond_graph(xyz_text)
    was = frozenset(tuple(sorted(one)) for one in (before or ()))
    if now == was:
        return True, ''
    return False, graph_changed(was, now, [str(r[0]) for r in rows])


def bonds_to_freeze(xyz_text: str) -> list:
    """Every bond, as a held value, so a relaxation cannot break one.

    GOAT's ``FREEZEBONDS``, which is on by default while it pushes a structure
    uphill: the way not to break a bond is not to let it move.  Refusing the
    step afterwards is the backstop and not the method -- measured on a
    2,4-hexadiene, refusing alone still let one bond go, because a bond can
    stretch across several accepted steps and the geometry handed back is
    then the last of those.

    Held at the length it has, which is what makes this cost nothing: the
    value is already met, so the restraint does no work until something tries
    to change it.
    """
    rows = [line.split() for line in atom_lines(xyz_text)]
    if not rows:
        return []
    from delfin.atom_mapping import cov_radius

    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    radius = [cov_radius(str(r[0])) for r in rows]
    held = []
    for i, j in sorted(bond_graph(xyz_text)):
        held.append({'kind': 'distance', 'atoms': [i, j], 'mode': 'fix',
                     'value': math.dist(where[i], where[j])})
    return held


def graph_changed(before: Any, after: Any, symbols: Any = None) -> str:
    """What happened between two bond graphs, said in words, or ''.

    Named by element and number rather than by index, because "C3-Br4 would
    break" is a sentence a chemist can act on and "(2, 3) left the set" is
    not.
    """
    made = sorted(set(after or ()) - set(before or ()))
    gone = sorted(set(before or ()) - set(after or ()))
    if not made and not gone:
        return ''

    def name(pair):
        if not symbols:
            return f'{pair[0] + 1}-{pair[1] + 1}'
        return '-'.join(f'{symbols[n]}{n + 1}' for n in pair
                        if 0 <= n < len(symbols))

    said = []
    if gone:
        said.append('breaks ' + ', '.join(name(one) for one in gone[:3]))
    if made:
        said.append('makes ' + ', '.join(name(one) for one in made[:3]))
    return ' and '.join(said)


def pair_named(pair: Any, symbols: Any = None) -> str:
    """Two atoms, as a chemist writes them: ``C3-Br4``.

    Counted from one, because that is what the picture and the atom list show
    and an index counted from zero in a sentence is a defect report waiting to
    be filed.  The same naming :func:`graph_changed` uses, out here where the
    rest of the file can reach it.
    """
    pair = tuple(int(n) for n in (pair or ()))
    if not symbols:
        return '-'.join(str(n + 1) for n in pair)
    return '-'.join(f'{symbols[n]}{n + 1}' for n in pair
                    if 0 <= n < len(symbols))


def _median(values: Any) -> float:
    """The middle of these numbers, or 0.0 when there are none."""
    kept = sorted(float(one) for one in values)
    if not kept:
        return 0.0
    middle = len(kept) // 2
    if len(kept) % 2:
        return kept[middle]
    return 0.5 * (kept[middle - 1] + kept[middle])


#: How far out of line a point has to be before a walk is said to have jumped.
#:
#: A driven scan holds one coordinate and relaxes everything else, and there
#: is no rule that the rest has to follow continuously.  When it does not --
#: some other coordinate slips over its own barrier between one point and the
#: next -- the profile has a step in it, and everything read off that profile
#: is about two different paths joined at a discontinuity.  Jonsson, Mills and
#: Jacobsen named it in 1998: "the path generated may be discontinuous ... some
#: atomic coordinates may slip near the saddle point region and the saddle
#: point configuration will then be missed."
#:
#: This is not the test for whether a scan can be quoted -- walking the second
#: leg and comparing is, and it separates cleanly (see :func:`paths_disagree`).
#: This is what *names the step*, so that the culprit can be named with it,
#: and its thresholds are set so that it never contradicts the second leg.
#:
#: Found on the second difference of the energy, which is what a step in a
#: curve is: the point before the fall and the point after it both disagree
#: with the straight line through their neighbours, and by the size of the
#: fall.  Measured over twenty-two legs -- eleven scans walked out and back
#: under GFN2, relaxed at every point -- with the two thresholds together,
#: large enough to be a different structure rather than a different geometry
#: *and* out of line with what this path itself does:
#:
#:     kept its bonding             largest second difference   ratio
#:       propane C-C-C angle                  0.3 kcal/mol       1.9x
#:       water dimer O-O                      0.6              20.6x
#:       glycol OCCO torsion                  1.1               1.8x
#:       ethane C-C stretch                   1.4               6.9x
#:       butanol, butane torsions             1.8, 2.1          2.3x
#:       SN2 Cl-/CH3Cl                        5.2              10.0x
#:     jumped
#:       Diels-Alder, backward               12.9              13.5x
#:       formate/water H transfer            19.0, 32.9         7.9x
#:       N-methylacetamide torsion           19.2              12.0x
#:       Diels-Alder, forward                33.5             130.7x
#:       ring opening                        60.6, 93.7        33.3x
#:
#: Eight legs of eight that jumped are found and none of the fourteen that did
#: not is called one.  Neither threshold does it alone, and that is the
#: finding: the ratios overlap -- a clean SN2 reaches 10.0x and a real
#: hydrogen transfer only 7.9x -- while the sizes do not, and a water dimer
#: whose two molecules turn over inside their hydrogen bond is 21 times its
#: own median at six tenths of a kcal/mol, which is not a jump by any reading.
#:
#: Nothing here is about a kind of coordinate.  An amide torsion jumps and an
#: SN2 does not; a distance driven through a Diels-Alder jumps and the same
#: kind of distance stretched in an ethane does not.  Which is why the second
#: leg is walked rather than a rule applied.
#:
#: The alternative that suggests itself, the RMSD between consecutive points,
#: is worse than either: measured on the Diels-Alder it is 0.536 A at the fall
#: and 0.685 A at a point where nothing happened, so it names the wrong step.
WALK_JUMPED_TIMES = 4.0
WALK_JUMPED_AT_LEAST = 10.0   # kcal/mol


def where_a_walk_jumped(spent: Any) -> Optional[Dict[str, Any]]:
    """Which step of a driven scan is a fall rather than a walk, or None.

    *spent* is the energy at each point in kcal/mol, in the order they were
    walked.  Returned as ``{'step', 'second', 'scale', 'times', 'fell'}``:
    *step* is the index of the point the walk landed on, so the discontinuity
    is between ``step - 1`` and ``step``; *fell* is how far the energy moved in
    that one step.

    Arithmetic on numbers the scan already has, and no geometry: a walk of
    forty points over a large structure cannot keep forty geometries, and it
    does not have to.
    """
    energies = []
    for one in (spent or ()):
        try:
            energies.append(float(one))
        except (TypeError, ValueError):
            return None
    if len(energies) < 4:
        # Three points make one second difference, which is its own median:
        # every path would then be exactly at its own scale and nothing could
        # ever stand out from it.
        return None
    second = [abs(energies[i + 1] - 2.0 * energies[i] + energies[i - 1])
              for i in range(1, len(energies) - 1)]
    scale = _median(second)
    worst = max(range(len(second)), key=lambda i: second[i])
    height = second[worst]
    if height < WALK_JUMPED_AT_LEAST:
        return None
    if scale > 0.0 and height < WALK_JUMPED_TIMES * scale:
        return None
    # A step shows in the second difference at both of its ends, so which of
    # the two steps around this point was the fall is decided by the plain
    # difference: the fall is the larger one.
    here = worst + 1                       # back into the energies' own index
    before = abs(energies[here] - energies[here - 1])
    after = (abs(energies[here + 1] - energies[here])
             if here + 1 < len(energies) else 0.0)
    step = here if before >= after else here + 1
    return {'step': step, 'second': height, 'scale': scale,
            'times': (height / scale) if scale > 0 else float('inf'),
            'fell': energies[step] - energies[step - 1]}


#: How near two atoms have to be, in either of two geometries, for the
#: distance between them to be worth watching.
#:
#: A jump is a bond made or broken, and neither happens at ten Angstrom.  The
#: bound is what keeps this arithmetic rather than a cost: every pair of atoms
#: is a square law, and at the 250-atom ceiling GFN2 is offered with that is
#: 31 000 pairs a step.  Six Angstrom is twice the longest contact the editor
#: draws and well past any bond.
_WORTH_WATCHING = 6.0


def what_else_moved(before: str, after: str,
                    driven: Any = ()) -> Optional[Dict[str, Any]]:
    """The internal coordinate that changed most between two points.

    Everything except the one the scan is driving, which is held and is
    supposed to change.  Returned as ``{'pair', 'was', 'now', 'moved'}``, or
    None when there is nothing to compare.

    This is what names the culprit.  A user told that a scan jumped can do
    nothing with that on its own; told *which* coordinate slipped, they can
    arm that one too and walk both together, which is what the editor's
    several-legs-at-once scan is for.  Measured at the fall of the
    Diels-Alder: the undriven forming C-C went from 2.915 to 1.558 A in one
    step, 1.357 A, against a median of 0.042 A over the rest of the path --
    32 times, on a path where nothing else came near it.

    Every pair of atoms is a square law, and it is still nothing beside the
    step it is describing: measured, 0.7 ms at 50 atoms, 12.7 at the 250 GFN2
    is offered with and 70.9 at GFN-FF's 600, against a relaxed scan point
    that is tenths of a second at the smallest of those.
    """
    was = coordinates_of(before or '')
    now = coordinates_of(after or '')
    if not was or len(was) != len(now):
        return None
    count = len(was) // 3
    keep = {tuple(sorted((int(a), int(b))))
            for a, b in _pairs_of(driven)}
    here = [(was[3 * i], was[3 * i + 1], was[3 * i + 2]) for i in range(count)]
    there = [(now[3 * i], now[3 * i + 1], now[3 * i + 2])
             for i in range(count)]
    best = None
    for i in range(count):
        for j in range(i + 1, count):
            if (i, j) in keep:
                continue
            one = math.dist(here[i], here[j])
            two = math.dist(there[i], there[j])
            if one > _WORTH_WATCHING and two > _WORTH_WATCHING:
                continue
            moved = abs(two - one)
            if best is None or moved > best['moved']:
                best = {'pair': (i, j), 'was': one, 'now': two, 'moved': moved}
    return best


def _pairs_of(driven: Any) -> List[tuple]:
    """Every pair of atoms the caller says is being driven.

    A leg of a scan is two atoms for a distance, three for an angle and four
    for a torsion, and in all three cases the atoms named are the ones the
    walk is dictating -- so every pair among them is a coordinate that is
    meant to change and is not evidence of anything.

    Given a list of legs or a single leg, because both readings of "the atoms
    being driven" are natural at the call site and getting the wrong one would
    quietly leave the driven coordinate in the comparison, where it is the
    largest change on every step of every scan.
    """
    legs = list(driven or ())
    if legs and all(isinstance(one, int) for one in legs):
        legs = [legs]
    out = []
    for leg in legs:
        try:
            numbers = [int(one) for one in leg]
        except (TypeError, ValueError):
            continue
        for a in range(len(numbers)):
            for b in range(a + 1, len(numbers)):
                out.append((numbers[a], numbers[b]))
    return out


#: How far two legs of the same walk may disagree and still be one path.
#:
#: A driven scan is a minimum-energy path only where the coordinate it is not
#: driving follows continuously; Bofill and Quapp (Mol. Phys. 2019) give the
#: condition exactly -- no turning point and no valley-ridge inflection.  Where
#: it fails the answer depends on which way the walk went, and the way to find
#: that out is to walk it back and compare.
#:
#: The number is what a difference in a barrier does to a rate.  At 298 K,
#: RT ln 10 is 1.36 kcal/mol: two barriers further apart than that are two
#: rates an order of magnitude apart, and two barriers nearer than it are the
#: same answer to any use a barrier is put to.  It is a temperature, so it is
#: worked out at the temperature the editor is set to rather than fixed here.
#:
#: Measured against it, on eleven scans run out and back under GFN2 -- torsions
#: of an alkane, an alcohol, a diol and an amide, a C-C stretch, a C-C-C angle,
#: a stretched hydrogen bond, an SN2, a ring opening, a hydrogen transfer and a
#: Diels-Alder:
#:
#:     ethane C-C stretch                     0.000 kcal/mol
#:     propane C-C-C angle                    0.001
#:     butanol, butane, glycol torsions       0.002, 0.004, 0.006
#:     SN2 Cl-/CH3Cl                          0.032
#:     water dimer O-O                        0.803
#:     ---------------------------------------------- RT ln 10 = 1.364
#:     N-methylacetamide torsion             14.2
#:     Diels-Alder                           23.8
#:     formate/water H transfer              60.8
#:     ring opening                         129.1
#:
#: A factor of eighteen of clear water between the two groups, so the
#: threshold is not a fitted number -- almost anything in the gap would do --
#: and the one with a meaning is the one worth using.  And nothing in the two
#: groups is about a kind of coordinate: an amide torsion jumps and an SN2
#: does not.
GAS_CONSTANT_KCAL = 1.987204259e-3


def a_rate_apart(kelvin: Any = 298.15) -> float:
    """How far apart two barriers have to be to be two different answers.

    A factor of ten in rate, which is where a difference stops being rounding
    and starts being chemistry.  Returned in kcal/mol at the temperature
    given, so a scan run hot is judged at the temperature it was run.
    """
    try:
        T = float(kelvin)
    except (TypeError, ValueError):
        T = 298.15
    return GAS_CONSTANT_KCAL * max(1.0, T) * math.log(10.0)


def paths_disagree(there: Any, back: Any) -> Optional[Dict[str, Any]]:
    """The largest gap between a walk and the same walk taken backwards.

    Both are ``[(value, energy)]`` on one zero, in the order they were walked;
    the second is the return leg, so its values run the other way.  Returned
    as ``{'at', 'gap', 'there', 'back', 'points'}`` over the coordinate values
    the two have in common, or None when they have none.

    Compared at the coordinate rather than at the top, because the top is
    exactly where the two legs are least likely to be about the same place.
    Measured on the Diels-Alder: forward puts its highest point at 2.20 A with
    the undriven forming bond at 2.92, backward at 2.90 A with the same bond
    at 1.76 -- two maxima 0.7 A apart on the driven coordinate and 1.2 A apart
    on the one nobody was driving.  Comparing the two heights would be
    comparing two different geometries; comparing them where the driven
    coordinate agrees is comparing the path with itself.
    """
    ours = [(float(v), float(e)) for v, e in (there or ())]
    theirs = [(float(v), float(e)) for v, e in (back or ())]
    if len(ours) < 2 or not theirs:
        return None
    spacing = _median([abs(b[0] - a[0]) for a, b in zip(ours, ours[1:])])
    near = 0.25 * spacing if spacing > 0 else 1e-6
    found = None
    for value, energy in ours:
        mate = min(theirs, key=lambda one: abs(one[0] - value))
        if abs(mate[0] - value) > near:
            continue
        gap = abs(energy - mate[1])
        if found is None or gap > found['gap']:
            found = {'at': value, 'gap': gap, 'there': energy,
                     'back': mate[1]}
    if found is None:
        return None
    found['points'] = sum(
        1 for value, _ in ours
        if abs(min(theirs, key=lambda one: abs(one[0] - value))[0]
               - value) <= near)
    return found


def gfnff_would_form(xyz_text: str, legs: Any) -> Optional[Dict[str, Any]]:
    """The leg of a scan that asks a fixed topology to grow a bond, or None.

    GFN-FF works its bonding out once, from the geometry it is first handed,
    and then holds the molecule together with it.  xtb's own documentation is
    blunt about what follows: "GFN-FF can only break bonds, dissociation
    reactions will therefore usually work fine, while association reactions
    are likely to fail."  ORCA's GOAT bars it from uphill steps for the same
    reason.

    So the boundary is a direction, not a method.  A distance that starts
    outside bonding range and is driven inside it is asking for a bond the
    force field has no term for; the same distance driven the other way is
    asking for one it has, and is exactly what a fast force field is good for.
    :data:`BOND_STARTS_AT` is where the editor already puts that line.

    Measured on butadiene and ethylene, one forming C-C driven from 3.40 A to
    1.60 in 0.1 A steps, everything relaxed at each point:

        GFN2       crosses at +7.3 kcal/mol at 2.20 A, ends -63.0 in the
                   product, and the *other* forming bond closes to 1.53 A
                   without being asked
        GFN-FF     climbs to +94.1 kcal/mol without crossing anything, and
                   the other forming bond ends at 3.39 A -- no reaction, and
                   87 kcal/mol of error in the one number a scan is for

    And the editor's topology cache makes that worse rather than better.  The
    cache exists so a drag cannot fall apart between one frame and the next,
    and it works: with the topology pinned the false profile is smooth and
    monotonic, which is what a wall of repulsion looks like when it is drawn
    carefully.  Unpinned, the same scan gives +108.6 with its maximum in a
    different place -- visibly wrong, and therefore less dangerous.
    """
    rows = [line.split() for line in atom_lines(xyz_text or '')]
    if not rows:
        return None
    from delfin.atom_mapping import cov_radius

    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    radius = [cov_radius(str(r[0])) for r in rows]
    for leg in (legs or ()):
        if str(leg.get('kind') or '') != 'distance':
            continue
        try:
            atoms = [int(one) for one in (leg.get('atoms') or ())]
            target = float(leg.get('to'))
        except (TypeError, ValueError):
            continue
        if len(atoms) != 2 or any(not (0 <= i < len(rows)) for i in atoms):
            continue
        i, j = atoms
        bonds_at = BOND_STARTS_AT * (radius[i] + radius[j])
        now = math.dist(where[i], where[j])
        if now >= bonds_at > target:
            return {'atoms': (i, j), 'now': now, 'to': target,
                    'bonds_at': bonds_at,
                    'symbols': [str(r[0]) for r in rows]}
    return None


def gfnff_pair_refusal(first: str, second: str) -> str:
    """Why GFN-FF cannot walk between these two structures, or ''.

    The same boundary as :func:`gfnff_would_form`, asked of a pair instead of
    a leg: when the second end has a bond the first has not, the path between
    them makes one, and a force field whose bonding was perceived at the first
    end has no term for it.

    Decidable here in a way it is not from a single structure, which is why
    this is refused and climbing from what is on screen is not: two ends say
    what the reaction is, and one geometry does not.

    Measured on the Diels-Alder, xtb's own path finder given the separated
    pair and cyclohexene:

        GFN2      a barrier of 20.8 kcal/mol and a reaction energy of -68.0
        GFN-FF    a barrier of 34.5 and a reaction energy of **+34.3**

    The sign is the whole of it.  GFN-FF never sees the product as a molecule,
    so it prices cyclohexene as a strained contact between the two things it
    still believes are there, and reports a reaction that goes uphill.
    """
    was, now = bond_graph(first or ''), bond_graph(second or '')
    made = sorted(set(now) - set(was))
    if not made:
        return ''
    rows = [line.split()[0] for line in atom_lines(second or '')]
    named = ', '.join(pair_named(one, rows) for one in made[:3])
    return (
        f'GFN-FF cannot walk between these two. The second end has bonds the '
        f'first has not -- {named} -- and GFN-FF perceives its bonding once '
        f'and then holds it, so it can break a bond and cannot make one. '
        f'Measured on a Diels-Alder: given the same two ends, GFN2 reports '
        f'the product 68 kcal/mol below the start and GFN-FF reports it 34 '
        f'above, because it never sees the product as a molecule. Choose '
        f'GFN2, GFN1 or g-xTB, or mark the two ends the other way round, '
        f'which is a reaction it can do.')


def gfnff_refusal(xyz_text: str, legs: Any) -> str:
    """Why GFN-FF cannot walk this scan, or '' when it can.

    Said before the run in the way an unparametrised solvent is (see
    :func:`delfin.dashboard.solvents.refusal`), and for the same reason: xtb
    does not refuse it.  It runs, it converges, it reports a number, and the
    number is a force field being squeezed.
    """
    found = gfnff_would_form(xyz_text, legs)
    if not found:
        return ''
    named = pair_named(found['atoms'], found['symbols'])
    return (
        f'GFN-FF cannot walk this one. It works its bonding out once and then '
        f'holds it, so it can break a bond and cannot make one -- xtb says so '
        f'itself: association reactions are likely to fail. This scan drives '
        f'{named} from {found["now"]:.2f} to {found["to"]:.2f} A, past the '
        f'{found["bonds_at"]:.2f} where those two would be bonded, and what '
        f'it would report is repulsion rather than a reaction. Measured on a '
        f'Diels-Alder: GFN2 crosses at +7.3 kcal/mol and lands 63 below the '
        f'start, GFN-FF climbs to +94 and crosses nothing. Choose GFN2, GFN1 '
        f'or g-xTB, or drive the bond that is breaking instead.')


def _occluded(where, radius, i, j) -> bool:
    """Whether a third atom sits between these two, so they are not in contact.

    This is what used to be "heavy atoms first", and it says the same thing
    without naming an element.  In an SN2 the arriving nucleophile closes on
    the carbon *and* on the leaving group behind it at exactly the same rate,
    so both look like the coordinate being driven; one of them has a carbon in
    the way.  Held anyway, the leaving group is nailed down -- measured on
    Cl- + CH3Br, C-Br stays at 1.94 to 1.99 A and the price runs to +175
    kcal/mol, against 2.01 -> 3.24 A and near zero when only the real contact
    is held.

    It is true of a bromide behind a carbon, a hydride behind a metal, and
    anything else behind anything else, which is the point.
    """
    axis = [where[j][n] - where[i][n] for n in range(3)]
    span = math.sqrt(sum(one * one for one in axis))
    if span <= 1e-9:
        return False
    axis = [one / span for one in axis]
    for k in range(len(where)):
        if k in (i, j):
            continue
        rel = [where[k][n] - where[i][n] for n in range(3)]
        along = sum(one * two for one, two in zip(rel, axis))
        if not (0.15 * span < along < 0.85 * span):
            continue
        off = math.sqrt(max(0.0, sum(one * one for one in rel) - along * along))
        if off < 0.75 * (radius[k] + min(radius[i], radius[j])):
            return True
    return False


def _swing(where, i, b, c) -> list:
    """Where a turn about the b-c axis takes atom *i*, per radian.

    A vector: which way the grabbed atom goes, and how fast.  Its length is
    the lever arm -- how far *i* stands from the axis -- so a torsion and a
    bond can be weighed against each other in one currency, and its direction
    is what says whether the turn goes the way the hand is pulling.
    """
    axis = [where[c][n] - where[b][n] for n in range(3)]
    span = math.sqrt(sum(one * one for one in axis)) or 1.0
    axis = [one / span for one in axis]
    rel = [where[i][n] - where[b][n] for n in range(3)]
    along = sum(one * two for one, two in zip(rel, axis))
    arm = [rel[n] - along * axis[n] for n in range(3)]
    return [axis[1] * arm[2] - axis[2] * arm[1],
            axis[2] * arm[0] - axis[0] * arm[2],
            axis[0] * arm[1] - axis[1] * arm[0]]


def _lever(where, i, b, c) -> float:
    """How far *i* stands from the b-c axis.

    What turns a torsion into a length, so an angle can be weighed against a
    bond in one currency.
    """
    return math.sqrt(sum(one * one for one in _swing(where, i, b, c)))


def _carries(vector, pulled) -> float:
    """How far this coordinate carries the grabbed atom the hand's own way.

    *vector* is where one unit of the coordinate takes the atom -- a radian of
    a torsion, an Angstrom of a distance -- and *pulled* is the unit vector
    the hand is moving it along.  What comes back is in Angstrom, for both
    kinds, which is what lets a turn and a stretch be compared at all.

    Without a direction it is the whole length, which is the older question --
    how far could this coordinate move the atom, in any direction -- and it is
    the right one to fall back to: a caller with no history does not know
    which way the hand went.
    """
    length = math.sqrt(sum(one * one for one in vector))
    if pulled is None or length <= 1e-12:
        return length
    return abs(sum(one * two for one, two in zip(vector, pulled)))


def turn_for(where, radius, grabbed, held, pulled=None) -> list:
    """The soft coordinate a hand on *grabbed* would really be moving.

    A bond is the wrong thing for a drag to drive.  It is the one coordinate
    that must not give, so a hand on it either does nothing -- a force below
    what a bond holds can only stretch one by a tenth of an angstrom -- or
    tears the molecule, if it is above.  Measured on a 2,4-hexadiene, a chain
    carbon dragged 1.75 A: 0.09 A of movement with the pull, and three bonds
    broken with the rigid hand.

    What moves is torsions.  That is what every conformer generator drives --
    RDKit's embedding, CREST's sampling, GOAT, which freezes bonds precisely
    in order to leave them free -- and it is what the user asked for in the
    first place: bond lengths hold, and the angles and torsions react.

    So: a torsion about the bond *one further in* than the one the grabbed
    atom hangs on, with the grabbed atom as its outer end -- turning that
    swings it through an arc.  Turning about its own attachment bond would
    spin the group and leave the atom where it was, which is why the bond it
    sits on is not the one to use.

    *Which* torsion is the question this used to skip.  There are usually
    several and it took the first the walk came to, in index order, without
    asking where it would take the atom -- and a torsion carries the grabbed
    atom in one direction only, at right angles to both the axis and the arm.
    A hand pulling any other way is pushing on a coordinate that cannot
    express it, which is a drag that does almost nothing and spends a great
    deal doing it.

    So each candidate is asked how far it would carry the grabbed atom *the
    way the hand is going*: :func:`_swing` says where a radian of it takes
    the atom, and :func:`_carries` takes the part of that which points along
    *pulled*.  The largest wins.  Without a *pulled* -- a caller with no
    history does not know which way the hand went -- it is the whole lever
    arm, which is the older question and the right thing to fall back to.

    That this is the question and not "is the axis in a ring" was measured on
    a butane.  The same end carbon of the same molecule, dragged across the
    chain and dragged along it, wants a torsion in the first case and a bond
    in the second; nothing topological can tell those two apart, because the
    topology is identical.  The direction can, and does: the swing about the
    C2-C3 axis is square to the chain's plane, so a hand across the chain has
    all of it and a hand along the chain has none.

    Nothing, when no such torsion exists, and then the caller keeps the
    contact it had.  Which is right: a whole methyl leaving the other half of
    an ethane has no torsion between the two, and the C-C bond really is the
    coordinate that describes that drag.  A bond is only the wrong answer
    where a better one exists -- and the caller weighs the two, in the one
    currency :func:`_carries` puts them both in.

    What this does not settle, said plainly because it was measured.  The
    carry says which coordinates *can* move the grabbed atom; it cannot say
    whether the molecule will let them, because that is an energy and there
    is none here.  A cyclohexane carbon and a benzene carbon pushed through
    their own ring planes are the same problem to every geometric measure --
    the same six-ring, the same offer, carries of 1.19 and 1.20 Angstrom per
    radian -- and the answers are 0.195 of the drag for +6.3 kcal/mol against
    0.042 for +13.9.  What separates them is conjugation, which no
    arrangement of coordinates can see.
    """
    if not grabbed:
        return []
    inside = set(held or ()) | set(grabbed)
    best: tuple = (0.0, None)
    for i in sorted(grabbed):
        near = [n for n in range(len(where))
                if n != i and _is_a_bond(where, radius, i, n)]
        # A terminal atom has no torsion that moves it: a hydrogen sits on
        # its bond, and turning about that bond leaves it exactly where it
        # was.  Pulled at, what changes is the bond and nothing else -- which
        # is the case the hold on "the neighbour it is leaving" was written
        # for, and it stays that.
        if len(near) < 2:
            continue
        for b in near:
            if b in grabbed:
                continue
            onward = [n for n in range(len(where))
                      if n not in (i, b) and n not in inside
                      and _is_a_bond(where, radius, b, n)]
            for c in onward:
                far = [n for n in range(len(where))
                       if n not in (i, b, c) and n not in inside
                       and _is_a_bond(where, radius, c, n)]
                if not far:
                    continue
                # Every d on this axis is the same turn and takes the grabbed
                # atom to the same place, so the axis is scored once and the
                # lowest-numbered d names it.
                carry = _carries(_swing(where, i, b, c), pulled)
                if carry > best[0]:
                    best = (carry, [i, b, c, far[0]])
    if best[1] is None:
        return []
    return [{'kind': 'dihedral', 'atoms': best[1], 'mode': 'drag',
             'value': _dihedral(where, *best[1])}]


def with_their_terminals(where, radius, dragged) -> set:
    """*dragged*, plus every atom that hangs off it and nothing else.

    A hand does not grab an atom out of a molecule; it grabs it with what is
    attached to it.  Left behind, the hydrogens of a dragged methyl end up at
    0.92 and 2.12 A from their carbon, one of them 1.87 A from where it
    started -- a picture nothing can be priced on.  On a ring it is worse than
    ugly: the axial hydrogen points the way the hand is pushing, so the whole
    push scores as a C-H stretch, that is the coordinate that gets held, and
    the torsion the hand meant never gets a look in.  A cyclohexane could not
    be flipped.

    Terminal, not hydrogen.  An atom whose only bond is to something being
    dragged comes along whatever element it is; one with bonds elsewhere stays
    and is negotiated with.  A methyl's hydrogens travel, a bridging atom does
    not, and nothing here knows what a hydrogen is.
    """
    held = {int(i) for i in (dragged or ())}
    out = set(held)
    for k in range(len(where)):
        if k in held:
            continue
        ties = [j for j in range(len(where))
                if j != k and _is_a_bond(where, radius, k, j)]
        if len(ties) == 1 and ties[0] in held:
            out.add(k)
    return out


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
    most: int = 3,
    was: Any = None,
    turning: Any = None,
    holding: Any = None,
    opening: bool = False,
    unchanged: bool = False,
) -> list:
    """The internal coordinates the hand moved, to hold while the rest relaxes.

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
    simply walks the 1.5 A over and takes its hydrogen back.  One atom fixed
    in space holds no bond stretched.

    What is left is what a relaxed scan does: hold the coordinates the hand
    moved, free everything else.  Which ones is decided by how far each one's
    change moves atoms -- one currency for bonds, contacts and torsions alike
    -- and by two geometric tests that replace what used to be a list of
    rules: :func:`_occluded` for a contact with something in the way, and the
    projection in ``_ALONG_ENOUGH`` for one the hand is merely sweeping past.

    Nothing here knows what a hydrogen or a metal is, which is the point:
    every rule that did know broke on the next kind of chemistry.  Two atoms
    closing on the same partner used to be forbidden, which is exactly what an
    oxidative addition is -- so the real Pd...Br contact went unheld, the
    bromide drove to 1.27 A of the metal, and the budget reported +13.2
    kcal/mol for a geometry worth several hundred.

    *opening* is the first answer of a drag, and it is a different question:
    what moved most is a poor guide when almost nothing has moved yet, and
    for a grabbed atom the coordinate that has changed most is usually the
    bond it hangs on rather than the one the hand is driving.  So the opening
    is decided by what each coordinate *can* do -- see :func:`turn_for` -- and
    the caller then makes that decision stick through *turning*.  It is still
    given *was*, because which way the hand went is exactly what tells a turn
    from a stretch, and a drag that has not been told cannot know.
    """
    rows = [line.split() for line in atom_lines(xyz_text)]
    if not rows:
        return []
    symbols = [str(r[0]) for r in rows]
    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    from delfin.atom_mapping import cov_radius
    radius = [cov_radius(one) for one in symbols]

    # A hand that turns does not stop turning between two frames.  Decided
    # afresh every answer the verdict flipped back and forth on a cyclohexane
    # -- turning, pulling, turning -- and the ring fell back towards the chair
    # each time the constraint set changed under it.
    if turning is not None and len(turning) == 4:
        want = [int(n) for n in turning]
        if all(0 <= n < len(where) for n in want):
            return [{'kind': 'dihedral', 'atoms': want, 'mode': 'drag',
                     'value': _dihedral(where, *want)}]

    # What came along, and what the hand actually took hold of.  The
    # travellers are part of the drag -- they must not read as the rest of the
    # molecule -- but they make no statement of their own: a methyl's
    # hydrogens riding on their carbon each claimed a coordinate, and the two
    # spurious distances they claimed pinned a cyclohexane so it could not be
    # flipped.
    grabbed = sorted({int(i) for i in (dragged or ())
                      if 0 <= int(i) < len(where)})
    held = with_their_terminals(where, radius, grabbed)
    inside = sorted(i for i in held if 0 <= i < len(where))
    outside = [j for j in range(len(where)) if j not in held]
    if not inside or not outside:
        # The hand is on everything, or on nothing.  Either way no coordinate
        # between the two describes the drag, and a caller that gets nothing
        # back is being told to price the geometry some other way.
        return []
    # The hand stood still, so the question has not changed and neither may
    # the answer.  Said by the caller, because only the caller knows where the
    # *hand* is: from in here the structure is all there is, and under a
    # budget the structure moves whether or not anybody is moving it.  Held
    # coordinates re-derived under that movement chase it -- see the editor's
    # own note where this is passed.
    if unchanged and holding:
        kept = [dict(one) for one in holding]
        for one in kept:
            atoms = [int(n) for n in one['atoms']]
            if any(not (0 <= n < len(where)) for n in atoms):
                break
            one['value'] = (math.dist(where[atoms[0]], where[atoms[1]])
                            if one['kind'] == 'distance'
                            else _dihedral(where, *atoms)
                            if one['kind'] == 'dihedral' else one['value'])
        else:
            return kept
    then = _snapshot(was, len(where))
    if then is None or opening:
        # The opening answer of a drag decides by what the coordinates *can*
        # do rather than by what they have already done -- the nearest contact
        # stands in, unless that contact is a *bond*.
        #
        # A bond is the one coordinate a drag must not drive: a hand below
        # what a bond holds can only stretch one by a tenth of an angstrom,
        # and one above it tears the molecule.  Measured on a 2,4-hexadiene, a
        # chain carbon dragged 1.75 A: 0.09 A of movement with the pull, three
        # bonds broken with the rigid hand.  What moves a group there is a
        # torsion; see :func:`turn_for`.
        #
        # A contact that is *not* a bond is a different matter and is exactly
        # right: two fragments closing on each other -- a Diels-Alder, an SN2
        # -- have no torsion between them, and the distance is the reaction.
        span, i, j = min(((math.dist(where[i], where[j]), i, j)
                          for i in grabbed for j in outside))
        # Which way the hand went, when the caller has said where the drag
        # came from.  A torsion carries the grabbed atom in one direction
        # only, so without this the choice between a turn and a stretch is
        # made blind -- and the blind choice is wrong the moment the hand
        # pulls the other way.
        pulled = None
        if then is not None:
            went = [where[i][n] - then[i][n] for n in range(3)]
            far = math.sqrt(sum(one * one for one in went))
            if far > 1e-9:
                pulled = [one / far for one in went]
        if _is_a_bond(where, radius, i, j):
            turn = turn_for(where, radius, grabbed, held, pulled=pulled)
            if turn:
                # Both in Angstrom the grabbed atom moves, so the two are
                # comparable and nothing has to be chosen: a radian of the
                # turn against an Angstrom of the stretch, which is the
                # currency the scored branch below already weighs them in.
                axis = [where[i][n] - where[j][n] for n in range(3)]
                length = math.sqrt(sum(one * one for one in axis)) or 1.0
                if (_carries(_swing(where, *turn[0]['atoms'][:3]), pulled)
                        > _carries([one / length for one in axis], pulled)):
                    return turn
        return [{'kind': 'distance', 'atoms': [i, j], 'value': span,
                 'mode': 'drag'}]

    best: dict = {}

    def offer(kind, atoms, value, score):
        key = (kind, tuple(atoms))
        if score > best.get(key, (0.0,))[0]:
            best[key] = (score, {'kind': kind, 'atoms': list(atoms),
                                 'value': value, 'mode': 'drag'})

    for i in grabbed:
        step = math.dist(where[i], then[i])
        if step <= 1e-9:
            continue
        for j in outside:
            bond = _is_a_bond(where, radius, i, j)
            if not bond:
                axis = [where[j][n] - where[i][n] for n in range(3)]
                span = math.sqrt(sum(one * one for one in axis)) or 1.0
                went = [where[i][n] - then[i][n] for n in range(3)]
                straight = abs(sum(one * two for one, two in zip(axis, went)))
                if straight / (span * step) < _ALONG_ENOUGH:
                    continue
            if _occluded(where, radius, i, j):
                continue
            offer('distance', (i, j), math.dist(where[i], where[j]),
                  abs(math.dist(where[i], where[j])
                      - math.dist(then[i], then[j])))
            if not bond:
                continue
            for k in range(len(where)):
                if k in (i, j) or k in held:
                    continue
                if not _is_a_bond(where, radius, j, k):
                    continue
                for ell in range(len(where)):
                    if ell in (i, j, k) or ell in held:
                        continue
                    if not _is_a_bond(where, radius, k, ell):
                        continue
                    now = _dihedral(where, i, j, k, ell)
                    offer('dihedral', (i, j, k, ell), now,
                          _lever(where, i, j, k)
                          * math.radians(_turned_by(
                              now, _dihedral(then, i, j, k, ell))))

    # One coordinate per dragged atom, best first.
    #
    # Not an element rule and not a guess: a hand puts each atom it holds
    # somewhere, and that is one statement about the structure each.  Held
    # three at a time on a single atom, the three pin it outright and the
    # relaxation has no freedom left -- a Claisen whose new bond was still
    # 4.3 A away priced at +133 kcal/mol, which is a strained cage rather
    # than a molecule on its way anywhere.  Two dragged atoms give two, which
    # is a cycloaddition; six give the two that carry it, because the rest
    # score below the share.
    ranked, spoken = [], set()
    for score, one in sorted(best.values(), key=lambda pair: -pair[0]):
        if score < _MOVED_ANGSTROM:
            continue
        owner = one['atoms'][0]
        if owner in spoken:
            continue
        spoken.add(owner)
        ranked.append((score, one))
    ranked = [one for _score, one in ranked]
    if not ranked:
        # A hand that has stopped moving has changed nothing, and nothing is
        # not an answer: with no values held the price falls back to a single
        # point on the geometry as it stands, which is the number this whole
        # thing exists to stop using.
        #
        # What it holds is what it was already holding.  Derived afresh, a
        # still hand can be handed a *different* set -- the nearest contact
        # rather than the one the drag was driving -- and the structure then
        # slides towards something nobody asked for while nothing is moving.
        # Nothing changed, so nothing about what is held should change either.
        if holding:
            kept = [dict(one) for one in holding]
            for one in kept:
                atoms = [int(n) for n in one['atoms']]
                if any(not (0 <= n < len(where)) for n in atoms):
                    return contacts_holding(xyz_text, dragged, most=most)
                one['value'] = (math.dist(where[atoms[0]], where[atoms[1]])
                                if one['kind'] == 'distance'
                                else _dihedral(where, *atoms)
                                if one['kind'] == 'dihedral' else one['value'])
            return kept
        return contacts_holding(xyz_text, dragged, most=most)
    top = sorted(best.values(), key=lambda one: -one[0])[0][0]
    kept = [one for one in ranked
            if best[(one['kind'], tuple(one['atoms']))][0] >= _MOVED_SHARE * top]
    return kept[:max(1, int(most))]

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


#: What xtb writes the atomic charges to, per Hamiltonian.
#:
#: Every run writes one of these into its own scratch directory and the
#: directory is then removed, so the numbers were computed and deleted on
#: every drag answer, every optimisation and every scan point.  The two names
#: are not two spellings of one thing: ``charges`` holds the Mulliken-like
#: partial charges of a wavefunction, and ``gfnff_charges`` holds the
#: electronegativity-equilibration charges GFN-FF is parametrised with, which
#: is what that method has instead of a wavefunction.  Both are per atom, in
#: the input order, one to a line -- measured on a methane, GFN2 gives the
#: carbon -0.153, GFN1 -0.130, g-xTB -0.359 and GFN-FF -0.092, which is the
#: spread of four different definitions of the same idea and the reason the
#: method is named wherever they are shown.
_CHARGE_FILES = ('charges', 'gfnff_charges')


def read_charges(folder: Path) -> Optional[List[float]]:
    """The partial charges of the run in *folder*, or None if it wrote none.

    Free: the file is already on disk when this is called and is one short
    line per atom -- 57 atoms is a 913-byte read.  Nothing here asks xtb for
    anything, and if it ever had to, this would be the wrong place for it.
    """
    for name in _CHARGE_FILES:
        found = folder / name
        if not found.is_file():
            continue
        values: List[float] = []
        try:
            text = found.read_text(encoding='utf-8', errors='replace')
        except OSError:
            return None
        for line in text.splitlines():
            word = line.strip()
            if not word:
                continue
            try:
                values.append(float(word))
            except ValueError:
                return None
        return values or None
    return None


def read_bond_orders(folder: Path) -> Optional[List[tuple]]:
    """Wiberg bond orders from the run in *folder*, as ``(i, j, order)``.

    Counted from zero, because everything else in this editor counts atoms
    from zero and xtb's ``wbo`` counts from one.

    None where the method has none.  Only the Hamiltonians have them: GFN1,
    GFN2 and g-xTB each write a ``wbo``, and GFN-FF writes none at all --
    measured, its scratch directory holds ``gfnff_charges`` and ``gfnff_topo``
    and nothing else.  Nothing is worked around there: GFN-FF's bonding is a
    topology perceived once and then held fixed, so it has no order of its own
    to report.

    These are a readout and nothing decides on them.  See
    :data:`BOND_WORTH_SAYING` for the measurement that settles why -- a
    closed-shell order does not fall when a bond breaks homolytically.

    Only pairs xtb actually printed are in the list; it prints nothing for
    pairs it found no bond order for, so an absent pair reads as zero.
    """
    found = folder / 'wbo'
    if not found.is_file():
        return None
    pairs: List[tuple] = []
    try:
        text = found.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return None
    for line in text.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            first, second, order = int(parts[0]), int(parts[1]), float(parts[2])
        except ValueError:
            continue
        if first < 1 or second < 1:
            continue
        pairs.append((first - 1, second - 1, order))
    return pairs or None


def bond_order_between(bonds: Any, i: Any, j: Any) -> Optional[float]:
    """What *bonds* says about this pair, or None when nothing was computed.

    Zero rather than None for a pair that is missing from a list that exists:
    xtb prints the pairs it found a bond order for, so a pair it left out has
    an order below what it prints.  Zero is the order, not a verdict about the
    bond -- see :data:`BOND_WORTH_SAYING` for why the two are not the same.
    """
    if bonds is None:
        return None
    try:
        first, second = int(i), int(j)
    except (TypeError, ValueError):
        return None
    for one, two, order in bonds:
        if (one, two) == (first, second) or (one, two) == (second, first):
            return float(order)
    return 0.0


#: How low a bond order has to be before it is worth putting in a sentence.
#:
#: A cut for *which* orders to mention, and nothing more than that.  It is
#: emphatically not a threshold for whether a bond exists, and the reason is
#: measured: a Wiberg order is a readout of the wavefunction and it is not a
#: test of whether two atoms are still bonded.
#:
#: Measured here under GFN2, the coordinate held and everything else relaxed
#: at each length:
#:
#:     ammonia borane, N-B          ethane, C-C
#:     d/A   order   gap/eV         d/A   order   gap/eV
#:     1.66  0.609    9.05          1.53  1.030   15.26
#:     1.86  0.511    7.69          1.73  1.010   10.39
#:     2.06  0.403    6.33          1.93  1.001    6.96
#:     2.36  0.244    5.11          2.03  0.999    5.73
#:     2.86  0.000    3.64          2.53  0.998    2.16
#:     3.36  0.000    3.08          3.03  1.000    0.75
#:
#: The right-hand column is the one to remember.  An ethane pulled to twice
#: its bond length still reads 1.00, because a single closed-shell determinant
#: holds that pair in one orbital however far the two carbons are taken; a
#: separate measurement on the same homolysis has the fractional occupation
#: density at 1.73 electrons -- a whole pair already come apart -- where the
#: order still reads 0.96.  So against the editor's geometric watch, which
#: calls the bond gone from 1.94 A, the order is not the more honest of the
#: two: it is the one that says everything is fine exactly where the method
#: has stopped describing the system.
#:
#: Which is why nothing here decides anything on a bond order.  The bond watch
#: stays :func:`_is_a_bond` and :func:`bond_graph`, on distances, and the order
#: is offered as what it is: a number a chemist can read, true and useful where
#: it is read as one -- 1.9 across a C=C, 0.61 for a dative N-B at its own
#: minimum, which is worth knowing and is not "a weak bond".
#:
#: The signal that does work for a bond coming apart is the fractional
#: occupation density, ``xtb --fod``: 0.000 at a C-C of 1.54 A, 1.727 at
#: 3.50 A, and exactly 2.000 on a ninety-degree-twisted ethylene.  It is not
#: built here -- it costs a further single point and has caveats of its own,
#: over-reporting by about fourfold on a closed-shell metal complex -- and it
#: is named so that nobody reaches for the bond order to do its job.
BOND_WORTH_SAYING = 0.7


def bond_order_note(order: Any, names: str = 'the pair',
                    gap: Any = None, was: Any = None) -> str:
    """A bond order as a readout, in a clause for the status line.

    Empty when there is no order to say anything about -- a force field, a
    method with no wavefunction, or an answer that has not arrived yet.

    It states the number and, where the number itself is a statement, what it
    is a statement of: an order well above one is multiple-bond character.  It
    never says a bond is there or gone, because a Wiberg order does not answer
    that question -- see :data:`BOND_WORTH_SAYING` for the two series that
    settle it.

    *gap* is the frontier gap of the same answer.  Where it has closed, the
    determinant the order was computed from has stopped describing the system,
    and the number should not be read at all; given one, this says so.
    """
    try:
        value = float(order)
    except (TypeError, ValueError):
        return ''
    if method_is_out_of_its_depth(gap, was):
        # Said with the gap in it rather than by calling that function's
        # sentence, which is about the energy and is said in its own right
        # wherever it belongs.  This one is about the number in this clause.
        try:
            said_gap = f'{float(gap):.1f} eV'
        except (TypeError, ValueError):
            said_gap = 'a gap this small'
        return (f'{names} reads {value:.2f}, which is not worth reading at a '
                f'frontier gap of {said_gap}: a closed-shell bond order holds '
                'a pair in one orbital however far the two are pulled, so a '
                'bond coming apart still reads about one -- measured on an '
                'ethane, 1.00 at a C-C of 3.03 A.')
    if value >= 1.4:
        return (f'{names} is at {value:.2f}, which is multiple-bond '
                'character.')
    # The bare number, because this clause is written onto the status line
    # several times a second while a drag is running, and the caveat that
    # belongs with it -- that an order is not a bond-existence test -- is a
    # sentence about the quantity rather than about this answer. It is said
    # once, where it can be read: on the control's own tooltip, and in the
    # line the "What is it?" press writes above the orders it lists.
    return f'{names} is at {value:.2f}.'


#: How steep a structure may be and still be a place rather than a slope, in
#: Hartree per Bohr per degree of freedom.
#:
#: The RMS gradient rather than the norm, so it means the same thing for five
#: atoms and for five hundred: xtb prints the norm, and the norm of a
#: converged structure grows with the square root of the number of atoms.
#:
#: The number is ORCA's ``TolRMSG`` and Gaussian's, which is also what
#: :mod:`delfin.dashboard.climb` converges its own walks on -- one number for
#: "this has stopped moving", wherever the question is asked here.  Ten times
#: it is where a Hessian stops being a statement about a stationary point:
#: measured on a hand-built methane with every C-H at 1.0897 A, the RMS
#: gradient is 2.6e-3 and xtb happily reports zero imaginary modes -- a
#: perfectly true statement about a geometry that is not a minimum.
GRADIENT_IS_STILL = 1.0e-4
GRADIENT_IS_A_SLOPE = 1.0e-3


def rms_gradient(norm: Any, atoms: Any) -> Optional[float]:
    """xtb's gradient norm as a per-coordinate RMS, or None."""
    try:
        count = int(atoms)
        return float(norm) / math.sqrt(3.0 * count) if count > 0 else None
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def not_a_stationary_point(norm: Any, atoms: Any) -> str:
    """What to say when a Hessian was taken somewhere nothing is resting.

    Empty when the structure is standing still, which is when the modes mean
    what a chemist reads them as meaning.

    A Hessian is the curvature at a point, and it is defined everywhere.  What
    is *not* defined everywhere is the reading: "one mode goes the wrong way,
    so this is a transition state" is a statement about a stationary point,
    and on a slope the same Hessian describes the side of a hill.  The
    frequencies are not wrong there -- they are simply not about a structure
    anything sits at, and a count of zero imaginary modes on a geometry with a
    large gradient says "downhill from here in every direction", which is true
    of most of the surface.

    So the honest answer for a structure that is still on a slope is that it
    is neither a minimum nor a saddle, and that is what this says.
    """
    rms = rms_gradient(norm, atoms)
    if rms is None or rms < GRADIENT_IS_A_SLOPE:
        return ''
    return (f'It is not standing still: the gradient is {rms:.1e} Hartree '
            f'per Bohr per coordinate, against the {GRADIENT_IS_STILL:.0e} an '
            'optimiser converges on. A Hessian here describes a point on a '
            'slope, so this structure is neither a minimum nor a saddle -- '
            'relax it first, and ask again about what it settles on.')


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
#:
#: Checked again on chemistry the table above has none of: a 57-atom manganese
#: complex with four bromines at charge +1, which is a very different SCC from
#: 57 atoms of hydrocarbon.  Twenty GFN2 cycles, best of three:
#:
#:     cores      1     2     4     8    16    32
#:     seconds 3.03  2.73  2.27  2.53  4.74  8.73
#:
#: The ladder hands this one eight and four is a tenth quicker, which is
#: inside the spread on a machine that is also running other people's work --
#: but sixteen is twice as slow and thirty-two nearly four times, so where the
#: curve turns over is the same place it was on the hydrocarbons.  Left alone:
#: the ladder is not costing anything worth a change, and the heavy elements
#: do not move where it should stop.
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


#: What xtb's path finder says when it has found one.
_PATH_FORWARD_RE = re.compile(r'forward\s+barrier \(kcal\)\s*:\s*(-?\d+\.\d+)')
_PATH_BACK_RE = re.compile(r'backward barrier \(kcal\)\s*:\s*(-?\d+\.\d+)')
_PATH_REACTION_RE = re.compile(r'reaction energy\s+\(kcal\)\s*:\s*(-?\d+\.\d+)')
_PATH_TAKEN_RE = re.compile(r'path\s+(\d+) taken with\s+(\d+) points')
_PATH_RMSD_RE = re.compile(
    r'run\s+(\d+)\s+barrier:\s*(-?\d+\.\d+)\s+dE:\s*(-?\d+\.\d+)'
    r'\s+product-end path RMSD:\s*(\d+\.\d+)')

#: How the path is asked for.  Two runs, because the first often finds a way
#: round rather than over: measured on a Diels-Alder, run 1 reported a barrier
#: of 51.49 kcal/mol and run 2 of 8.32, and xtb takes the better of them.
PATH_RUNS = 2
PATH_POINTS = 25


def walk_the_path(reactant: str, product: str, method: str = 'gfn2', *,
                  charge: int = 0, uhf: int = 0,
                  solvent: Optional[str] = None,
                  solvation_model: str = 'alpb',
                  points: int = PATH_POINTS, runs: int = PATH_RUNS,
                  timeout: Optional[float] = 600.0) -> Dict[str, Any]:
    """Find a way from *reactant* to *product*, and say what it costs.

    xtb's own path finder: metadynamics with a bias that pushes away from
    where it has been and pulls towards the product, run a few times with
    different strengths, and the best of them kept.  It answers the question a
    scan can only approach -- a scan drives a coordinate somebody chose, and
    this one is given two structures and finds its own way between them.

    Which is why it belongs after a scan rather than instead of one: the scan
    is what produces a product to aim at.  Measured on a butadiene and an
    ethylene, 2.6 s: forward barrier 5.754 kcal/mol, backward 69.586, reaction
    energy -63.831, and an estimated transition state with the two forming
    bonds at 2.524 and 2.520 A.  The scan of the same reaction put its highest
    point at +6.3 kcal/mol and 2.36 A -- two methods that share no machinery,
    agreeing.

    It is the same twice.  Metadynamics suggests otherwise and it was written
    here that way, wrongly: three goes at one pair of structures gave 5.755
    kcal/mol every time, to the digit, and four goes at another gave 43.7.
    xtb seeds it fixed.

    What does move it is which two structures it is given.  The same reaction
    from a slightly different reactant complex gave 5.755 and 3.3 -- so a
    barrier from this is a statement about the two ends as much as about the
    reaction, and the RMSD below is what says the second end was reached at
    all.  Within its runs xtb keeps the lowest barrier it found, which makes
    the number an upper bound over the paths it tried.

    Returns ``{'ok', 'barrier', 'back', 'reaction', 'ts', 'points', 'rmsd',
    'seconds', 'status'}``.  *rmsd* is how near the path actually came to the
    product it was aimed at, which is the one number that says whether the
    answer is about the reaction that was asked for.
    """
    key = str(method or 'gfn2').lower()
    spec = GFN_METHODS.get(key) or GFN_METHODS['gfn2']
    label = spec['label']
    # A solvent this method cannot be given, said before the walk rather than
    # blamed on the two structures afterwards.  :func:`optimize_with_gfn` has
    # refused this for a long time and this did not: driven, g-xTB handed
    # ``--alpb water`` stops with "No ALPB/GBSA parameters found for the
    # method/solvent" and the path finder reported "g-xTB found no path
    # between the two structures" -- which reads as "these two do not react",
    # and is a statement about the chemistry where the truth was about the
    # build.
    wet = str(solvent or '').strip().lower()
    if wet and spec.get('solvation') is False:
        return {'ok': False,
                'status': (f'{label} has no implicit solvation in this build, '
                           'so a path through a solvent is not something it '
                           'can walk. Walk it in the gas phase, or choose '
                           'GFN2-xTB or GFN-FF.')}
    no = _solvents.refusal(solvation_model, solvent, key)
    if no:
        return {'ok': False, 'status': f'{label}: {no}'}
    # And only then look for the program.  A solvent this method has not got
    # is true whether or not xtb is installed, so answering it first is both
    # the more useful sentence and the one that does not change with the
    # machine: asked on a box without xtb, the binary check spoke first and
    # the refusal that was actually about the request never came out.
    binary = find_binary(key)
    if binary is None:
        return {'ok': False, 'status': (f'A path needs xtb, which was not '
                                        f'found in {_where_it_looked()}.')}
    here = [line for line in atom_lines(reactant)]
    there = [line for line in atom_lines(product)]
    if not here or len(here) != len(there):
        return {'ok': False,
                'status': ('A path needs two structures of the same molecule; '
                           f'these have {len(here)} and {len(there)} atoms.')}
    # And in the same order.  A path is a walk from atom 1 to atom 1 and atom
    # 2 to atom 2; given the same molecule built twice, the two orderings
    # differ and there is no walk between them at all -- xtb answers "no path
    # found", which reads as "these do not react" and is not what happened.
    mine = [line.split()[0] for line in here]
    yours = [line.split()[0] for line in there]
    if mine != yours:
        first = next((n for n, (a, b) in enumerate(zip(mine, yours)) if a != b),
                     0)
        return {'ok': False,
                'status': ('A path walks atom 1 to atom 1, so the two '
                           'structures have to be the same molecule in the '
                           f'same order -- atom {first + 1} is {mine[first]} '
                           f'in one and {yours[first]} in the other. Build the '
                           'second one from the first rather than separately.')}

    started = time.perf_counter()
    folder = Path(tempfile.mkdtemp(prefix='delfin-path-'))
    try:
        (folder / 'from.xyz').write_text(
            f'{len(here)}\nfrom\n' + '\n'.join(here) + '\n', encoding='utf-8')
        (folder / 'to.xyz').write_text(
            f'{len(there)}\nto\n' + '\n'.join(there) + '\n', encoding='utf-8')
        (folder / 'path.inp').write_text(
            '$path\n'
            f'   nrun={max(1, int(runs))}\n'
            f'   npoint={max(5, int(points))}\n'
            '   anopt=10\n'
            '   kpush=0.003\n'
            '   kpull=-0.015\n'
            '   ppull=0.05\n'
            '$end\n', encoding='utf-8')
        cores = interactive_cores(len(here))
        command = [binary, 'from.xyz', *spec['flags'], '--path', 'to.xyz',
                   '--input', 'path.inp',
                   '--chrg', str(int(charge)), '--uhf', str(max(0, int(uhf))),
                   '-P', str(cores)]
        command += _solvents.xtb_flags(solvation_model, solvent)
        environment = dict(os.environ, OMP_NUM_THREADS=str(cores),
                           MKL_NUM_THREADS=str(cores), OMP_STACKSIZE='1G')
        own = parameter_home(binary)
        if own:
            environment['XTBPATH'] = own
        try:
            done = subprocess.run(command, cwd=str(folder), env=environment,
                                  capture_output=True, text=True,
                                  timeout=timeout)
        except subprocess.TimeoutExpired:
            return {'ok': False,
                    'status': (f'The path finder ran past {timeout:.0f} s and '
                               'was stopped.')}
        except OSError as exc:
            return {'ok': False, 'status': f'The path finder did not run: {exc}'}
        output = (done.stdout or '') + (done.stderr or '')
        forward = _PATH_FORWARD_RE.search(output)
        if not forward:
            # It ran and found nothing worth reporting.  Which is an answer:
            # these two structures are not joined by anything this method can
            # walk, and saying so beats a number that is not there.
            return {'ok': False, 'seconds': time.perf_counter() - started,
                    'output': output[-4000:],
                    'status': (f'{label} found no path between the two '
                               'structures.')}
        back = _PATH_BACK_RE.search(output)
        reaction = _PATH_REACTION_RE.search(output)
        taken = _PATH_TAKEN_RE.search(output)
        rmsd = None
        if taken:
            for one in _PATH_RMSD_RE.finditer(output):
                if one.group(1) == taken.group(1):
                    rmsd = float(one.group(4))
        guess = folder / 'xtbpath_ts.xyz'
        ts = guess.read_text(encoding='utf-8') if guess.is_file() else None
        seconds = time.perf_counter() - started
        return {
            'ok': True, 'seconds': seconds,
            'barrier': float(forward.group(1)),
            'back': float(back.group(1)) if back else None,
            'reaction': float(reaction.group(1)) if reaction else None,
            'points': int(taken.group(2)) if taken else None,
            'rmsd': rmsd, 'ts': ts, 'method': key,
            'status': (f'{label} walked a path in {seconds:.1f} s.'),
        }
    finally:
        shutil.rmtree(folder, ignore_errors=True)


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
    free_energy: bool = False,
    thermo_kelvin: float = 298.15,
    fod: bool = False,
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

    A successful run also carries what it computed on the way and used to
    throw away with its scratch directory: ``charges`` per atom, ``bonds`` as
    Wiberg orders for the pairs that have one, and ``gradient`` as the norm at
    the geometry it ended on.  None of the three costs a calculation -- see
    :func:`read_charges` and :func:`read_bond_orders` -- and with
    *free_energy* on, ``thermo`` carries H, T*S and the zero-point energy out
    of the same block the free energy was already being read from.

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

    The same directory carries the converged wavefunction between runs for the
    methods that have one, so a sequence of calculations on neighbouring
    geometries is not a sequence of unrelated ones -- 59 % off a single point
    and 7 % off a twenty-cycle relaxation on a manganese complex under GFN2.
    See :func:`_restart_named` for what it is keyed by, why the key has to be
    that, and for the measurement that it changes no number.

    *solvent* is one of :data:`SOLVENTS` and *solvation_model* one of ALPB,
    GBSA or ddCOSMO -- see :mod:`delfin.dashboard.solvents` for which method
    can be run with which, and what each was measured to cost.  A geometry
    optimised in the gas phase and one optimised in water are different
    answers, and so are two optimised under different models, so both are
    named in the status rather than left in the operator's memory.

    *fod* asks the same run how much of the structure is not a closed shell --
    see :data:`_NFOD_RE` for what the number is and :data:`FOD_METHODS` for
    which methods answer at all.  It is a single point by construction: the
    extra SCF runs at 5000 K, so optimising under it would be relaxing on a
    surface nobody wants, and the combination is refused rather than run.  The
    result gains ``'fod'`` as ``{'total', 'on'}`` with the per-atom breakdown
    in input order.  It is one more SCF and costs like one: measured as the
    minimum of seven runs against the same structure run ``--sp``, so that
    both estimates carry the same share of a shared machine, 0.96 s against
    0.49 for sixteen atoms and 10.5 against 3.5 for a 57-atom manganese
    complex -- two to three single points, and a small fraction of one
    optimisation step, which was 7.4 s and 35.6 s on the same two.
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

    # Asked of a method that cannot answer, ``--fod`` produces silence rather
    # than a refusal -- see :data:`FOD_METHODS` for what each of the two does
    # -- so the refusal is made here, where it can be a sentence.
    if fod and key not in FOD_METHODS:
        return {
            'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
            'seconds': 0.0, 'frames': [],
            'status': (
                f'{label} cannot be asked how much of this is not a closed '
                f'shell: it has no wavefunction to occupy fractionally, and '
                f'it answers the question with silence rather than a refusal. '
                f'Ask GFN2-xTB or GFN1-xTB.'),
        }
    if fod and optimise:
        return {
            'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
            'seconds': 0.0, 'frames': [],
            'status': ('Fractional occupation is measured at 5000 K, so a '
                       'geometry cannot be optimised under it. Ask for it as '
                       'a single point.'),
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
        # A free energy is a Hessian, and a Hessian is not free: measured
        # under GFN2, 0.57 s against 0.29 for sixteen atoms and 3.72 against
        # 0.76 for twenty-four.  So it is asked for and never assumed.
        command = [binary, source.name, *spec['flags'],
                   *(['--ohess'] if (optimise and free_energy)
                     else ['--opt'] if optimise
                     else ['--hess'] if free_energy else []),
                   *(['--fod'] if fod else []),
                   '--chrg', str(int(charge)), '--uhf', str(max(0, int(uhf))),
                   '-P', str(cores)]
        if max_steps:
            command += ['--cycles', str(int(max_steps))]
        command += _solvents.xtb_flags(model, wet)
        held = constraint_input(constraints, atoms=len(body))
        # One input file, whatever is in it: a second --input replaces the
        # first rather than adding to it.
        asked = held['text']
        if free_energy:
            asked += f'$thermo\n  temp={float(thermo_kelvin):g}\n$end\n'
        if asked:
            (folder / 'xtb.inp').write_text(asked, encoding='utf-8')
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
        # And the last wavefunction, for the methods that have one.  See
        # :func:`_restart_named` for what it is keyed by and why the key has
        # to be that and not less.
        #
        # Never for a fractional-occupation run.  That one converges at
        # 5000 K, and its wavefunction is not a guess at the ordinary one --
        # it is the smeared answer to a different question.  Written into the
        # store it would be handed to the next single point as its starting
        # guess, keyed as though it belonged to the same run, and nothing
        # downstream could tell.  So the check on this structure neither reads
        # nor writes what the walk around it is using.
        warm = (None if fod else
                _restart_named(topology, key, charge, uhf, wet, model))
        if warm is not None and warm.is_file():
            try:
                shutil.copy2(str(warm), str(folder / 'xtbrestart'))
            except OSError:
                pass
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
        if warm is not None:
            _keep_restart(folder, warm, output)
        # The two things this run computed and this directory was about to
        # take to the grave.  Read here, before the ``finally`` below removes
        # the folder, and read unconditionally: they are already written, they
        # are one line per atom and one line per bond, and nothing is asked of
        # xtb to get them.
        #
        # Measured: reading both is 0.053 ms for a methane and 0.158 ms for
        # the 57-atom manganese complex, against drag answers of 175 ms and
        # 2.28 s -- three hundredths of a percent and seven thousandths of
        # one.  The same drag timed with and without the two reads, alternated
        # so a shared machine could not favour either half, does not separate
        # them at all: 25 rounds on the methane gave 175.1 ms against 181.7,
        # and four rounds on the manganese complex gave 4.107 s against 4.053
        # -- the read is smaller than the noise it would have to be measured
        # in.  Which is the whole design: if reading these ever cost a
        # calculation, this would be the wrong place for it.
        charges = read_charges(folder)
        bonds = read_bond_orders(folder)
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
        gap = None
        # The last one the run printed, which is the geometry it ended on --
        # see :data:`_GAP_RE` for why there is more than one.
        for told_gap in _GAP_RE.finditer(output):
            try:
                gap = float(told_gap.group(1))
            except ValueError:
                gap = None
        # How steep it is where the run ended.  The last one for the same
        # reason the gap takes the last one: an optimisation prints the
        # summary at every geometry it passes through, and only the final one
        # is about the structure that comes back.
        gradient = None
        for told_gradient in _GRADIENT_RE.finditer(output):
            try:
                gradient = float(told_gradient.group(1))
            except ValueError:
                gradient = None
        spread = None
        if fod:
            spread = _read_fod(output)
        free = None
        wrong_way = None
        thermo = None
        if free_energy:
            told = _FREE_ENERGY_RE.search(output)
            if told:
                try:
                    free = float(told.group(1))
                except ValueError:
                    free = None
            # And the rest of what the Hessian paid for.  H and the
            # zero-point energy are printed in the same block as G and were
            # being read past; T*S is the difference of the two totals, which
            # is what it is by definition and what xtb's own table agrees
            # with to every digit it prints.
            told_enthalpy = _ENTHALPY_RE.search(output)
            told_zpe = _ZPE_RE.search(output)
            enthalpy = zpe = None
            try:
                enthalpy = float(told_enthalpy.group(1)) if told_enthalpy else None
            except ValueError:
                enthalpy = None
            try:
                zpe = float(told_zpe.group(1)) if told_zpe else None
            except ValueError:
                zpe = None
            if free is not None or enthalpy is not None or zpe is not None:
                warmth = float(thermo_kelvin) if thermo_kelvin else 0.0
                ts = (enthalpy - free
                      if (enthalpy is not None and free is not None) else None)
                thermo = {
                    'kelvin': warmth,
                    'free_energy': free,
                    'enthalpy': enthalpy,
                    'zpe': zpe,
                    'ts': ts,
                    # In Hartree per Kelvin, so every energy this module hands
                    # back is in one unit and whoever shows it decides what a
                    # chemist reads.
                    'entropy': (ts / warmth if (ts is not None and warmth > 0)
                                else None),
                }
            counted = _IMAGINARY_RE.search(output)
            if counted:
                # And the modes themselves, most negative first: "one
                # imaginary frequency" is what a transition state has, and
                # *which* one is what says it is the right transition state.
                seen = []
                for row in _EIGVAL_RE.finditer(output):
                    for word in row.group(1).split():
                        try:
                            seen.append(float(word))
                        except ValueError:
                            pass
                # Distinct values: xtb prints the list more than once, and
                # the same mode twice reads as two of them.  How many there
                # are is xtb's own count against its own cutoff, and is not
                # worked out again here.
                wrong_way = {
                    'count': int(counted.group(1)),
                    'modes': sorted({one for one in seen if one < -5.0})[:4],
                    # And the softest real ones, which say how flat what is
                    # left is.  One imaginary mode makes a transition state
                    # only if the directions around it are directions; a
                    # stationary point with a real mode of a few wavenumbers
                    # is on the edge of being a saddle of higher order, and
                    # published searches refuse it for that reason.
                    'real': lowest_real_modes(output),
                }
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

        # Asked for and not printed.  The list in :data:`FOD_METHODS` is the
        # two methods that were seen to answer, and it is a list rather than a
        # rule because the ways of not answering are silent: this is the same
        # refusal made again against what the run actually said, so a build
        # that stops printing it fails loudly instead of reporting a molecule
        # with no static correlation.
        if fod and spread is None:
            return {
                'ok': False, 'xyz': xyz_text, 'energy': None, 'method': key,
                'seconds': time.perf_counter() - started, 'engine': 'xtb',
                'frames': [], 'version': version,
                'status': (f'{label} was asked how much of this is not a '
                           f'closed shell and printed no answer. '
                           f'{which_xtb_ran(binary, output)}'),
            }

        # A single point has no geometry to converge, and calling it
        # unconverged would put it among the runs that stopped short.
        converged = (not optimise) or 'GEOMETRY OPTIMIZATION CONVERGED' in output
        if not converged:
            # The geometry is still better than the one that went in, so it is
            # handed back -- but not as though it were finished.
            return {
                'ok': True, 'xyz': relaxed, 'energy': energy, 'free_energy': free,
                'imaginary': wrong_way, 'gap': gap, 'method': key,
                'charges': charges, 'bonds': bonds, 'gradient': gradient,
                'thermo': thermo, 'fod': spread,
                'seconds': seconds, 'engine': 'xtb', 'frames': frames, 'version': version,
                'hamiltonian': reported or wanted or label, 'held': held,
                'converged': False, 'solvent': wet, 'solvation_model': model,
                'status': (f'{label} stopped before converging after '
                           f'{seconds:.1f} s; the geometry it reached is shown. '
                           f'(xtb {version}, {reported or wanted or label})'
                           + solvent_note(wet, model) + held_note(held)),
            }
        return {
            'ok': True, 'xyz': relaxed, 'energy': energy, 'free_energy': free,
            'imaginary': wrong_way, 'gap': gap, 'method': key,
            'charges': charges, 'bonds': bonds, 'gradient': gradient,
            'thermo': thermo, 'fod': spread,
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
    should_stop: Optional[Callable[[], bool]] = None,
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

    How much of a slideshow, measured on real structures at the twenty cycles
    a priced drag uses -- one whole follow answer, seconds:

        atoms   16     30     49     57    105
        GFN-FF   0.062  0.099  0.138  0.158  0.361
        GFN2     0.165  0.838  1.685  2.966 17.46

    A hundred atoms under GFN2 is a drag that answers once every seventeen
    seconds, and under GFN-FF the same drag answers three times a second: at
    that size the choice of method is the whole difference between an editor
    and a queue.  The four in the middle are transition-metal complexes, which
    is where the choice is hardest to make -- GFN-FF knows about the metal
    where UFF guesses, and GFN2 knows more still.
    """
    result = optimize_with_gfn(
        xyz_text, method, charge=charge, uhf=uhf,
        max_steps=max(1, int(cycles)), timeout=timeout,
        constraints=constraints, topology=topology, solvent=solvent,
        solvation_model=solvation_model, should_stop=should_stop,
    )
    result['converged'] = bool(
        result.get('ok') and 'converged in' in str(result.get('status') or '')
    )
    return result


#: How close two atoms may come, as a fraction of the bond they would make.
#:
#: A net that needs no energy, so it holds where the energy is weakest -- a
#: transition metal, an open shell, anything GFN2 is shaky about.  Squeezing is
#: never chemistry: stretching is what a reaction does and the budget already
#: prices it, but two atoms inside two thirds of a bond length are not on any
#: path at any temperature.
#:
#: Measured on the paths that must stay open, the closest any pair came:
#: a cyclohexane flipping 0.94, an SN2 driven through its saddle and past its
#: product 0.74, a [1,5]-hydrogen shift 0.78.  On the one that must not: a ring
#: carbon pushed across its own ring, 0.35.  And the diatomic cost at these
#: fractions, GFN2: 0.8 is between +10 and +45 kcal/mol, 0.6 between +153 and
#: +1083.  So this sits below everything real and above everything broken.
_TOO_CLOSE = 0.65


def closest_contact(xyz_text: str):
    """The tightest pair in the structure, as (fraction, i, j).

    A fraction of one is two atoms exactly a bond apart; below one they are
    closer than that.  ``(None, None, None)`` for a structure with fewer than
    two atoms.
    """
    from delfin.atom_mapping import cov_radius
    rows = [line.split() for line in atom_lines(xyz_text)]
    if len(rows) < 2:
        return None, None, None
    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    radius = [cov_radius(str(r[0])) for r in rows]
    worst = (None, None, None)
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            reach = radius[i] + radius[j]
            if reach <= 0:
                continue
            share = math.dist(where[i], where[j]) / reach
            if worst[0] is None or share < worst[0]:
                worst = (share, i, j)
    return worst


def settle_onto(xyz_text: str, reference: str, indices: Any) -> str:
    """*xyz_text* moved as a rigid body so the named atoms land where
    *reference* has them.

    A held value is an *internal* coordinate: xtb meets the distance it was
    given and is free to put the molecule anywhere that meets it.  So an
    answer comes back with the whole structure slid a little -- the geometry
    is exactly the one that was priced, and the atom under the cursor is no
    longer under the cursor.  Measured on a backside attack, chloride driven
    at a carbon: the answer moves the chloride 0.008 A away at the start and
    0.102 A through the transition region.  The page has already drawn it at
    the cursor, the answer draws it beside, the next mouse move draws it back
    -- about seven times a second, which is a molecule that shakes.

    Turning and shifting the whole structure fixes that and costs nothing:
    the energy of a molecule does not depend on where it is, so the geometry
    that was priced and the geometry that is drawn stay the same structure in
    every sense that matters.  This is what :func:`hold_atoms_at` should have
    been -- that one moves single atoms, which does change the structure, and
    it is kept only for the drags nothing is held for.

    Kabsch, on the atoms the hand is holding: one of them gives a shift, two
    a shift and a turn about their axis, three or more the whole rigid body.
    """
    wanted = sorted({int(i) for i in (indices or ())})
    rows = [line.split() for line in atom_lines(xyz_text)]
    source = [line.split() for line in atom_lines(reference)]
    if not wanted or not rows or len(rows) != len(source):
        return xyz_text
    if any(not (0 <= i < len(rows)) for i in wanted):
        return xyz_text
    here = [[float(r[1]), float(r[2]), float(r[3])] for r in rows]
    there = [[float(r[1]), float(r[2]), float(r[3])] for r in source]

    def centre(points, of):
        return [sum(points[i][n] for i in of) / len(of) for n in range(3)]

    mine, yours = centre(here, wanted), centre(there, wanted)
    moved = [[one[n] - mine[n] for n in range(3)] for one in here]
    if len(wanted) >= 2:
        # The 3x3 correlation of the two clouds, and its best rotation.
        cov = [[sum(moved[i][a] * (there[i][b] - yours[b]) for i in wanted)
                for b in range(3)] for a in range(3)]
        turn = _best_rotation(cov)
        if turn is not None:
            moved = [[sum(turn[a][b] * one[b] for b in range(3))
                      for a in range(3)] for one in moved]
    body = []
    for row, point in zip(rows, moved):
        body.append(f'{row[0]:<5}' + ''.join(
            f'{point[n] + yours[n]:>24.14f}' for n in range(3)))
    head = str(xyz_text or '').splitlines()
    comment = head[1] if len(head) > 1 else ''
    return f'{len(rows)}\n{comment}\n' + '\n'.join(body) + '\n'


def _best_rotation(cov):
    """The rotation matrix Kabsch gives for a 3x3 correlation, or None.

    Solved by Jacobi on the symmetric matrix rather than by pulling in a
    linear-algebra dependency for nine numbers.
    """
    try:
        import numpy as _np
    except Exception:
        return None
    try:
        matrix = _np.array(cov, dtype=float)
        u, _s, vt = _np.linalg.svd(matrix)
        sign = 1.0 if _np.linalg.det(u @ vt) > 0 else -1.0
        correct = _np.diag([1.0, 1.0, sign])
        return (u @ correct @ vt).T.tolist()
    except Exception:
        return None


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
