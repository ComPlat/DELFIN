"""Climb to a transition state one step at a time, so a hand can reach in.

:mod:`saddle` already presses a button and comes back with a converged saddle:
ORCA's OptTS driving xtb, seven to eight seconds for sixteen atoms.  That is a
press, and it is not something a hand can be in the middle of.  Measured on the
same sixteen atoms: a *three-step burst* of ``! XTB2 OPTTS`` costs 3.08 s cold
and 2.66 s warm with the previous burst's Hessian read back in, and about 2.5 s
of that is ORCA starting up.  A drag answers ten times a second.  Repeated ORCA
invocations are two orders of magnitude away from that and no amount of Hessian
reuse closes the gap, so an interactive climb cannot be ORCA restarted over and
over -- it has to be an optimiser here, walking on xtb gradients.

So this file is a saddle optimiser, and the method is the one production codes
use rather than one invented here: **partitioned rational function optimisation**
with **eigenvector following** (Banerjee, Adams, Simons and Shepard 1985;
Baker, *J. Comput. Chem.* **7**, 385, 1986) on a **Bofill-updated** Hessian
(*J. Comput. Chem.* **15**, 1, 1994).  One Hessian at the start, then one
gradient per step: the step maximises the energy along one chosen mode and
minimises it along every other, and the Hessian is updated from the step rather
than recomputed.  It is what Gaussian's Berny and ORCA's OptTS do, and the
gradient-only form of it is what pysisyphus and Sella do; their sources were
read for the step control and the mode tracking rather than derived here.

What it costs, measured on the sixteen-atom Diels-Alder estimate under GFN2:

* one gradient in process, four threads: **6 ms**; single threaded 8.9 ms
* one gradient through the xtb command line: **35 ms**, of which about 25 is
  the process starting -- still interactive, and the fallback when the xtb
  Python module is not installed
* one whole step -- gradient, step, Bofill update, projector: **10 ms**
* the starting Hessian, central differences, 6N gradients: **0.6 s** at sixteen
  atoms, and it grows as the cube of that -- roughly 17 s by sixty atoms, which
  is where a Hessian stops being something to press a button for
* a whole climb from the path finder's estimate: **11 steps, 1.7 s** end to
  end, of which 0.6 is the starting Hessian and 0.3 is xtb's own Hessian for
  the verdict at the end -- so the climb itself is a seventh of a second,
  against ORCA's 8.0 s for the same structure

and it lands in the same place: 2.3144 and 2.3161 A for the two forming bonds
against ORCA's 2.3153 and 2.3153, 0.0057 A RMSD, and the same one imaginary
mode at -394.0 cm-1.  (ORCA *prints* -372 for that run.  The difference is not
the geometry: xtb's own ``--hess`` on ORCA's own final structure says -393.5.
``Recalc_Hess 5`` means ORCA's last printed Hessian belongs to step 10 of a
13-step run, so -372 describes a geometry it had not finished leaving.)

Threads are pinned, and that is not a nicety.  On a 384-core node xtb's OpenMP
takes every core for a sixteen-atom molecule and one gradient costs **1.66 s**
-- 230 times what the same gradient costs on four threads.  Setting
``OMP_NUM_THREADS`` after the runtime has started does nothing (measured: 1.9 s
either way); ``omp_set_num_threads`` through the loaded OpenMP library works
(measured: 7.2 ms).  So that is what happens here, around each call and put
back afterwards.
"""

from __future__ import annotations

import ctypes
import ctypes.util
import math
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from . import gfn_optimize as _gfn

#: Angstrom in a Bohr, from the same place the rest of the editor takes it.
BOHR = _gfn.BOHR_IN_ANGSTROM

#: A square root of a Hartree per Bohr squared per electron mass, in cm-1.
HARTREE_IN_CM = 219474.6313702

#: Electron masses in an atomic mass unit.
ELECTRON_MASSES_PER_AMU = 1822.888486209

#: The methods this can climb on, and what xtb calls them.  The same three
#: :mod:`saddle` offers, for the same reason: anything with a basis set is a
#: job rather than a press, and a drag cannot wait for one.
CLIMB_METHODS = {
    'gfn2': 'GFN2xTB',
    'gfn1': 'GFN1xTB',
    'gfnff': 'GFNFF',
}

#: How far a single step may move any atom, in Bohr, before it is scaled back.
#: Baker's own step control: solve the step, and shorten it if it is longer
#: than the region the quadratic model is trusted in.
START_TRUST = 0.3
LARGEST_TRUST = 0.5
SMALLEST_TRUST = 0.02

#: ORCA's ``TolMaxG`` and ``TolRMSG``, in Hartree per Bohr.  Converging on the
#: same numbers as the button next door means the two can be compared at all.
#: 3e-4 on the largest gradient component is also Baker's own and Gaussian's,
#: so it is the one number the field agrees on.
GRADIENT_MAX = 3.0e-4
GRADIENT_RMS = 1.0e-4

#: How much of the mode it was climbing a step may leave behind before the step
#: is taken back rather than accepted.  Following a mode the structure has
#: moved out from under is how a climb converges onto a different saddle from
#: the one it set out for, and undoing the step and halving the trust radius is
#: what the literature does about it -- MOPAC's ``overlp``, whose ``omin`` for
#: a transition state is 0.8.
#:
#: Half, not four fifths, and the difference is measured.  MOPAC follows modes
#: in internal coordinates, where a mode keeps its shape while the molecule
#: turns; these are Cartesians, where a finite rotation changes every component
#: of an eigenvector without changing the mode at all.  At 0.8 the rule fires
#: on that rotation rather than on anything chemical: on eight test climbs it
#: refused 98, 97 and 97 steps out of 200 in three of them, and one that
#: reaches the Diels-Alder saddle at 0.5 stalled 0.126 A short of it at 0.8.
#: At 0.5 -- more than half of the mode still being the same mode -- it did not
#: fire on any of the eight, and the ratio test below did the useful work.  It
#: is kept as the guard against a genuine swap, which is what an overlap near
#: zero is.
LEAST_OVERLAP = 0.5

#: A step whose energy change is this far from what the quadratic model
#: predicted is not believed.  MOPAC's ``rmin``/``rmax`` for a transition
#: state, which are wider than a minimiser's because a saddle search is
#: allowed to go uphill.
WORST_RATIO = 0.0
BEST_RATIO = 4.0

#: Below this the disagreement is too small to be worth undoing a step for,
#: in Hartree -- MOPAC's own guard, so that a converging climb is not taken
#: back for a rounding error.
RATIO_MATTERS_ABOVE = 1.0e-5

#: Below this a mode counts as imaginary rather than as noise, in cm-1.
#:
#: A numerical Hessian puts the six modes that are nothing at all within a few
#: cm-1 of zero rather than at it, so some line has to be drawn or noise is
#: counted as chemistry.  autodE draws it at -40 and calls that deliberately
#: conservative; AutoMeKin ships ``imagmin 200`` and throws away anything
#: softer than that.  The conservative number is used, and the frequency itself
#: is always reported beside the count -- because the two choices disagree
#: about real cases, and a mode at -50 cm-1 is a first-order saddle by the
#: letter of the rule and is usually two fragments rocking against each other.
IMAGINARY_BELOW = -40.0

#: Element symbols by atomic number, and masses to match, which is all the
#: periodic table this needs: xtb takes atomic numbers and the mass weighting
#: takes masses.
_SYMBOLS = (
    '', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
    'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
    'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
    'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
    'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
    'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
)
_MASSES = (
    0.0, 1.008, 4.003, 6.94, 9.012, 10.81, 12.011, 14.007, 15.999, 18.998,
    20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948,
    39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933,
    58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904, 83.798,
    85.468, 87.62, 88.906, 91.224, 92.906, 95.95, 98.0, 101.07, 102.906,
    106.42, 107.868, 112.414, 114.818, 118.710, 121.76, 127.60, 126.904,
    131.293, 132.905, 137.327, 138.905, 140.116, 140.908, 144.242, 145.0,
    150.36, 151.964, 157.25, 158.925, 162.500, 164.930, 167.259, 168.934,
    173.045, 174.967, 178.49, 180.948, 183.84, 186.207, 190.23, 192.217,
    195.084, 196.967, 200.592, 204.38, 207.2, 208.980, 209.0, 210.0, 222.0,
    223.0, 226.0, 227.0, 232.038, 231.036, 238.029, 237.0, 244.0,
)
_NUMBER_OF = {name: z for z, name in enumerate(_SYMBOLS) if name}


def _elements(xyz_text: str) -> Optional[Dict[str, Any]]:
    """Symbols, atomic numbers, masses and coordinates, or ``None``.

    The atom count in the header is not trusted, for the reason
    :func:`gfn_optimize.atom_lines` gives: a block whose header says one thing
    and whose body says another loses atoms silently.
    """
    body = _gfn.atom_lines(xyz_text)
    if not body:
        return None
    names: List[str] = []
    where: List[List[float]] = []
    for line in body:
        parts = line.split()
        if len(parts) < 4:
            return None
        name = parts[0].strip().capitalize()
        if name not in _NUMBER_OF:
            return None
        try:
            where.append([float(parts[1]), float(parts[2]), float(parts[3])])
        except ValueError:
            return None
        names.append(name)
    numbers = np.array([_NUMBER_OF[one] for one in names], dtype=int)
    return {'symbols': names, 'numbers': numbers,
            'masses': np.array([_MASSES[z] for z in numbers], dtype=float),
            'angstrom': np.array(where, dtype=float)}


def xyz_document(symbols: Sequence[str], angstrom: np.ndarray,
                 comment: str = '') -> str:
    """The structure back as an XYZ block, with a header that counts the body."""
    rows = [f'{one:<3s} {r[0]:18.12f} {r[1]:18.12f} {r[2]:18.12f}'
            for one, r in zip(symbols, np.asarray(angstrom, dtype=float))]
    return f'{len(rows)}\n{comment}\n' + '\n'.join(rows) + '\n'


# --------------------------------------------------------------------------
# Gradients
# --------------------------------------------------------------------------

class _Threads:
    """Hold OpenMP to a few threads for the length of one call, then let go.

    On a 384-core login node xtb takes all of them for a sixteen-atom molecule
    and one gradient costs 1.66 s where four threads cost 5.5 ms -- 230 times.
    ``OMP_NUM_THREADS`` cannot fix it once the runtime is up: measured, setting
    it after the first calculation changed nothing at all (1.9 s either way),
    because OpenMP reads the variable when it starts and not afterwards.
    ``omp_set_num_threads`` does work (7.2 ms), so the library is asked
    directly and its previous setting is put back -- this process is not the
    only thing that may be using it.
    """

    def __init__(self, cores: int = 4):
        self.cores = max(1, int(cores))
        self.lib = None
        for name in ('gomp', 'omp', 'iomp5'):
            found = ctypes.util.find_library(name)
            if not found:
                continue
            try:
                lib = ctypes.CDLL(found)
                lib.omp_get_max_threads.restype = ctypes.c_int
                lib.omp_get_max_threads()
                self.lib = lib
                break
            except (OSError, AttributeError):
                continue
        self.was: Optional[int] = None

    def __enter__(self):
        if self.lib is not None:
            try:
                self.was = int(self.lib.omp_get_max_threads())
                self.lib.omp_set_num_threads(ctypes.c_int(self.cores))
            except (OSError, AttributeError, ValueError):
                self.was = None
        return self

    def __exit__(self, *_):
        if self.lib is not None and self.was:
            try:
                self.lib.omp_set_num_threads(ctypes.c_int(self.was))
            except (OSError, AttributeError):
                pass
        return False


def have_fast_gradients() -> bool:
    """Whether xtb can be called in process rather than started each time.

    Six times faster at sixteen atoms, and the difference is nearly all
    process startup, so it grows in importance as the molecule shrinks.  Not a
    requirement: the command line answers the same numbers to ten decimal
    places (measured: 7e-11 Hartree, 1.3e-7 Hartree/Bohr on the gradient) and
    is fast enough to drag against.
    """
    try:
        import xtb.interface        # noqa: F401
    except Exception:
        return False
    return True


class _InProcess:
    """xtb through its own library, which is where the speed is."""

    def __init__(self, numbers, bohr, method, charge, uhf, solvent, cores):
        from xtb.interface import Calculator, Param
        from xtb.libxtb import VERBOSITY_MUTED

        self.threads = _Threads(cores)
        name = CLIMB_METHODS.get(method) or 'GFN2xTB'
        with self.threads:
            self.calc = Calculator(getattr(Param, name), numbers,
                                   np.asarray(bohr, dtype=float),
                                   charge=float(charge), uhf=int(uhf))
            self.calc.set_verbosity(VERBOSITY_MUTED)
            # Nothing louder is asked for and nothing quieter is available:
            # ``set_output`` is the API for binding the Fortran side's own
            # output elsewhere, and bound to a file the very next ``update``
            # came back "Update of molecular structure failed" and every climb
            # stopped.  The one method that prints anyway is GFN-FF, and it
            # does not come through here -- see the engine choice in
            # :class:`Climb`.
            if solvent:
                try:
                    from xtb.utils import get_solvent
                    known = get_solvent(str(solvent))
                    if known is not None:
                        self.calc.set_solvent(known)
                except Exception:
                    pass
        self.last = None
        self.calls = 0

    def __call__(self, bohr):
        with self.threads:
            self.calc.update(np.asarray(bohr, dtype=float))
            self.last = self.calc.singlepoint(self.last)
            self.calls += 1
            return float(self.last.get_energy()), self.last.get_gradient().copy()


class _CommandLine:
    """xtb started once per gradient, for where the Python module is missing.

    35 ms for sixteen atoms against 5.5 in process, and about 25 of those 35
    are the process itself -- so this is the slow path and it is still under a
    frame.  Turbomole-format ``gradient`` is what ``--grad`` writes, and it is
    read rather than the printed energy, which carries fewer digits.
    """

    def __init__(self, numbers, bohr, method, charge, uhf, solvent, cores):
        self.binary = _gfn.find_xtb()
        self.symbols = [_SYMBOLS[int(z)] for z in numbers]
        self.method = method
        self.charge = int(charge)
        self.uhf = int(uhf)
        self.solvent = solvent
        self.cores = max(1, int(cores))
        self.folder = Path(tempfile.mkdtemp(prefix='delfin-climb-'))
        self.calls = 0

    def close(self):
        shutil.rmtree(self.folder, ignore_errors=True)

    def __call__(self, bohr):
        if self.binary is None:
            raise RuntimeError('no xtb')
        angstrom = np.asarray(bohr, dtype=float) * BOHR
        (self.folder / 'in.xyz').write_text(
            xyz_document(self.symbols, angstrom, 'from the DELFIN viewer'),
            encoding='utf-8')
        flag = {'gfn2': ['--gfn', '2'], 'gfn1': ['--gfn', '1'],
                'gfnff': ['--gfnff']}.get(self.method, ['--gfn', '2'])
        order = [self.binary, 'in.xyz', *flag, '--grad',
                 '--chrg', str(self.charge), '--uhf', str(self.uhf)]
        if self.solvent:
            order += ['--alpb', str(self.solvent)]
        room = dict(os.environ)
        room['OMP_NUM_THREADS'] = str(self.cores)
        room['MKL_NUM_THREADS'] = str(self.cores)
        subprocess.run(order, cwd=str(self.folder), env=room,
                       capture_output=True, text=True, timeout=120)
        self.calls += 1
        return _read_turbomole_gradient(self.folder / 'gradient',
                                        len(self.symbols))


def _read_turbomole_gradient(path: Path, count: int):
    """Energy and gradient out of what ``xtb --grad`` leaves behind."""
    rows = path.read_text(encoding='utf-8', errors='replace').splitlines()
    energy = None
    for row in rows:
        if 'energy =' in row:
            try:
                energy = float(row.split('energy =')[1].split()[0])
            except (IndexError, ValueError):
                pass
    body = [row for row in rows if row.strip() and not row.startswith('$')]
    wanted = body[1 + count:1 + 2 * count]
    grad = []
    for row in wanted:
        grad.append([float(word.replace('D', 'E').replace('d', 'E'))
                     for word in row.split()])
    if energy is None or len(grad) != count:
        raise RuntimeError('xtb wrote no gradient')
    return energy, np.array(grad, dtype=float)


# --------------------------------------------------------------------------
# The method
# --------------------------------------------------------------------------

def translations_and_rotations(masses: np.ndarray,
                               bohr: np.ndarray) -> np.ndarray:
    """An orthonormal basis of everything that is not a translation or a turn.

    Six vectors describe moving the whole molecule and turning it, and neither
    changes the energy.  The standard projector is built from them -- three
    mass-weighted translations and three infinitesimal rotations about the
    centre of mass -- but what is returned is a *basis* of what survives rather
    than the projector itself, and that is the point.

    ``P H P`` leaves six eigenvalues sitting at zero, and zero is exactly where
    the shift for the minimised modes lives: every one of them is a pole of
    ``sum F_i^2 / (lambda - b_i)`` sitting on top of the root being searched
    for.  Production codes get round it either by excluding those modes by
    index or by shifting them far up; working in a basis that never contained
    them does the first of those and leaves nothing to recognise afterwards.

    Linear molecules give five rather than six, which falls out of the
    orthogonalisation rather than being special-cased.
    """
    count = bohr.shape[0]
    root = np.sqrt(masses)
    centre = (masses[:, None] * bohr).sum(0) / masses.sum()
    arm = bohr - centre
    raw = []
    for axis in range(3):
        one = np.zeros((count, 3))
        one[:, axis] = root
        raw.append(one.reshape(-1))
    for axis in range(3):
        unit = np.zeros(3)
        unit[axis] = 1.0
        raw.append((np.cross(np.broadcast_to(unit, arm.shape), arm)
                    * root[:, None]).reshape(-1))
    kept: List[np.ndarray] = []
    for one in raw:
        for already in kept:
            one = one - already * float(already @ one)
        length = float(np.linalg.norm(one))
        if length > 1e-6:
            kept.append(one / length)
    whole = np.eye(3 * count)
    for one in kept:
        whole -= np.outer(one, one)
    value, vector = np.linalg.eigh(whole)
    return vector[:, value > 0.5]


def lowest_lambda(curvatures: np.ndarray, forces: np.ndarray) -> float:
    """The shift for the modes being minimised, from RFO's secular equation.

    ``lambda = sum F_i^2 / (lambda - b_i)`` with ``lambda`` below every
    ``b_i`` and below zero, which is the lowest root of the augmented matrix
    of Banerjee et al.  On that interval every term is negative and the
    left-hand side rises, so the function is strictly decreasing and bisection
    finds the root without ever needing a derivative -- and without the
    stepping past a pole that Newton does here.
    """
    if curvatures.size == 0:
        return 0.0
    edge = min(float(curvatures.min()), 0.0)

    def value(here: float) -> float:
        return float(np.sum(forces ** 2 / (here - curvatures))) - here

    high = np.nextafter(edge, -np.inf)
    low = edge - max(1e-8, abs(edge))
    for _ in range(400):
        if value(low) > 0:
            break
        low = edge - 2.0 * (edge - low)
    for _ in range(200):
        middle = 0.5 * (low + high)
        if value(middle) > 0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def partitioned_step(hessian: np.ndarray, gradient: np.ndarray,
                     following: Optional[np.ndarray]) -> Dict[str, Any]:
    """One P-RFO step: uphill along one mode, downhill along all the others.

    In the Hessian's own eigenbasis, with curvatures ``b`` and gradient
    components ``F``.  The mode being climbed gets its own two-by-two
    augmented matrix ``[[b_k, F_k], [F_k, 0]]`` whose larger root is
    ``(b_k + sqrt(b_k^2 + 4 F_k^2)) / 2`` and is never negative; everything
    else shares the augmented matrix whose lowest root is
    :func:`lowest_lambda` and is never positive.  Both give the same step,
    ``h_i = -F_i / (b_i - lambda)``, and the sign of the denominator is what
    makes one direction go up and the rest go down.

    *following* is the mode climbed last time.  Which mode to climb is chosen
    by overlap with it and not by its position in the list, because the list
    reorders as the structure moves: measured on a Diels-Alder, the reaction
    mode is the *second* lowest at the geometry a hand leaves, and a climb that
    takes the lowest one walks back down to the van-der-Waals complex.  With
    nothing to compare against, the lowest is taken, which for a structure that
    is already near a saddle is the imaginary one.
    """
    curvature, direction = np.linalg.eigh(hessian)
    force = direction.T @ gradient
    if following is None:
        chosen = 0
        overlap = None
    else:
        against = np.abs(direction.T @ following)
        chosen = int(np.argmax(against))
        overlap = float(against[chosen])
    up = 0.5 * (curvature[chosen]
                + math.sqrt(curvature[chosen] ** 2 + 4.0 * force[chosen] ** 2))
    others = np.array([i for i in range(curvature.size) if i != chosen],
                      dtype=int)
    down = (lowest_lambda(curvature[others], force[others])
            if others.size else 0.0)
    move = np.zeros_like(curvature)
    below = curvature[chosen] - up
    move[chosen] = -force[chosen] / below if abs(below) > 1e-300 else 0.0
    if others.size:
        move[others] = -force[others] / (curvature[others] - down)
    return {'step': direction @ move, 'mode': direction[:, chosen],
            'curvature': float(curvature[chosen]), 'overlap': overlap,
            'negative': int((curvature < -1e-10).sum())}


def bofill_update(hessian: np.ndarray, moved: np.ndarray,
                  changed: np.ndarray) -> np.ndarray:
    """Bofill's Hessian update: what makes one Hessian last a whole climb.

    A minimiser can use BFGS because it may assume the Hessian stays positive
    definite; a saddle search may not, and the update it needs has to be able
    to keep an eigenvalue below zero.  Bofill mixes the symmetric rank-one
    update of Murtagh and Sargent, which can, with the Powell-symmetric-Broyden
    update, which is stable where rank-one is not, in the proportion

        phi = (j.dx)^2 / ((j.j)(dx.dx))

    of the first to ``1 - phi`` of the second, where ``j = dg - H dx`` is what
    the old Hessian failed to predict.  Checked against pysisyphus's
    ``hessian_updates.py``, which is the same three expressions.
    """
    missed = changed - hessian @ moved
    square = float(moved @ moved)
    if square < 1e-20:
        return hessian
    along = float(missed @ moved)
    size = float(missed @ missed)
    powell = ((np.outer(moved, missed) + np.outer(missed, moved)) / square
              - along * np.outer(moved, moved) / square ** 2)
    if abs(along) < 1e-18 or size * square < 1e-40:
        return hessian + powell
    share = (along ** 2) / (size * square)
    rank_one = np.outer(missed, missed) / along
    return hessian + share * rank_one + (1.0 - share) * powell


class Climb:
    """A saddle search that takes one step at a time.

    Built to be driven from a loop that also has to answer a mouse, so nothing
    here blocks for longer than one gradient: :meth:`start` pays for the
    Hessian once, :meth:`step` costs one gradient, and :meth:`took` hands the
    climb a structure the user has just moved.

    The working coordinates are Cartesian, mass-weighted by ``m / mean(m)``.
    Both halves of that were measured rather than assumed.  Mass weighting is
    needed because the mode being climbed has to be the mode the frequency list
    names: unweighted, the ordering is by force constant, and on a Diels-Alder
    that put a hydrogen wag where the reaction coordinate should have been.
    *Relative* masses are used because RFO's step is not scale invariant -- its
    denominator is ``1 + h.h``, so the coordinates carry a length scale, and in
    electron masses a chemically ordinary step is a hundred and fifty units
    long and every step comes back over-damped.  Measured on the Diels-Alder
    estimate: 11 steps in relative masses, 77 in electron masses, for the same
    saddle.  Dividing by the mean leaves the eigenvectors and their order
    exactly as they were and puts the length scale back where RFO expects it.
    """

    def __init__(self, xyz_text: str, method: str = 'gfn2', *,
                 charge: int = 0, uhf: int = 0, solvent: Optional[str] = None,
                 cores: int = 4, trust: float = START_TRUST):
        found = _elements(xyz_text)
        if found is None:
            raise ValueError('there is no structure to climb from')
        self.symbols: List[str] = found['symbols']
        self.numbers = found['numbers']
        self.masses = found['masses']
        self.method = str(method or 'gfn2').lower()
        self.charge = int(charge)
        self.uhf = int(uhf)
        self.solvent = solvent
        self.cores = max(1, int(cores))
        self.bohr = found['angstrom'] / BOHR
        # The metric the optimiser works in, and the one frequencies want.
        self.scale = np.repeat(np.sqrt(self.masses / self.masses.mean()), 3)
        self.weight = np.repeat(
            np.sqrt(self.masses * ELECTRON_MASSES_PER_AMU), 3)
        self.trust = float(trust)
        self.hessian: Optional[np.ndarray] = None
        self.following: Optional[np.ndarray] = None
        self.energy: Optional[float] = None
        self.gradient: Optional[np.ndarray] = None
        self.steps = 0
        self.refused = 0
        self._before: Any = None
        # GFN-FF goes out through the command line even when the library is
        # there, and it is the file it writes that decides that: measured, an
        # in-process GFN-FF drops ``gfnff_topo`` -- 142 kB of it -- into
        # whatever directory the process happens to be in, which for a
        # dashboard is the directory the user launched it from.  The library
        # takes no working directory to be told about; the command line has
        # one of its own, and left nothing behind at all in the same test.
        # It also keeps the setup banner out of the kernel's log, since a
        # subprocess's output is captured rather than shared.  Three gradients
        # a second slower is worth not writing into somebody's project.
        fast = have_fast_gradients() and self.method != 'gfnff'
        self.engine = (_InProcess if fast else _CommandLine)(
            self.numbers, self.bohr, self.method, charge, uhf, solvent, cores)

    # -- what the caller reads -------------------------------------------

    @property
    def angstrom(self) -> np.ndarray:
        return self.bohr * BOHR

    def xyz(self, comment: str = 'Climbing to a transition state') -> str:
        return xyz_document(self.symbols, self.angstrom, comment)

    def frame(self) -> List[float]:
        """The flat coordinate list the viewer's frame channel is fed with."""
        return [round(float(one), 4) for one in self.angstrom.reshape(-1)]

    def close(self) -> None:
        closer = getattr(self.engine, 'close', None)
        if closer is not None:
            closer()

    # -- the calculation --------------------------------------------------

    def _basis(self) -> np.ndarray:
        whole = translations_and_rotations(self.masses, self.bohr)
        # The projector is built in the frequency metric; carried into the
        # working one it is no longer orthonormal, so it is made so again.
        carried = np.diag(self.weight / self.scale) @ whole
        made, _ = np.linalg.qr(carried)
        return made

    def _measure(self, bohr: np.ndarray):
        energy, gradient = self.engine(bohr)
        return float(energy), np.asarray(gradient, dtype=float).reshape(-1)

    def numerical_hessian(self, delta: float = 0.005) -> np.ndarray:
        """Central differences of the gradient, in the working metric.

        6N gradients and no cleverness: measured, 0.6 s for sixteen atoms in
        process and about 17 s by sixty, which is where this stops being
        something to press a button for.  Central rather than forward
        differences because the whole method rests on the sign of one
        eigenvalue, and forward differences get that wrong near a saddle where
        the curvature is small.  Checked against xtb's own ``--hess``: -394.76
        against -394.6 cm-1 on the same structure -- but see :meth:`verdict`
        for where the two part company, which is on symmetric structures and
        is a property of GFN2's gradient rather than of this differencing.
        """
        flat = self.bohr.reshape(-1)
        size = flat.size
        built = np.zeros((size, size))
        for axis in range(size):
            up = flat.copy()
            up[axis] += delta
            down = flat.copy()
            down[axis] -= delta
            _, forward = self._measure(up.reshape(-1, 3))
            _, backward = self._measure(down.reshape(-1, 3))
            built[axis] = (forward - backward) / (2.0 * delta)
        built = 0.5 * (built + built.T)
        return built / np.outer(self.scale, self.scale)

    def start(self, aimed_from: Any = None) -> Dict[str, Any]:
        """Pay for the Hessian, and choose the mode to climb.

        *aimed_from* is the structure as it stood before the user moved it.
        Given one, the displacement between then and now is the direction the
        hand asked for, and the mode climbed is whichever of the Hessian's
        eigenvectors most resembles it -- which is the whole of how a drag
        guides a saddle search.  See :meth:`took`.
        """
        began = time.perf_counter()
        self.hessian = self.numerical_hessian()
        self.energy, self.gradient = self._measure(self.bohr)
        self.steps = 0
        self.refused = 0
        self._before = None
        self.trust = START_TRUST
        self.following = None
        overlap = self.aim(aimed_from)
        modes = self.frequencies()
        return {'seconds': time.perf_counter() - began,
                'gradients': int(getattr(self.engine, 'calls', 0)),
                'aimed': overlap, 'modes': [float(one) for one in modes[:4]],
                'imaginary': int((modes < IMAGINARY_BELOW).sum())}

    def aim(self, aimed_from: Any) -> Optional[float]:
        """Point the climb along the way a hand has just moved the structure.

        Returns how much of that direction the chosen mode actually is, so a
        caller can say when the hand did not name anything in particular.  A
        pull along one bond in a molecule with a hundred modes will overlap
        weakly with all of them, and a weak overlap is worth reporting rather
        than acting on silently.
        """
        if aimed_from is None or self.hessian is None:
            return None
        before = _elements(aimed_from) if isinstance(aimed_from, str) \
            else {'angstrom': np.asarray(aimed_from, dtype=float)}
        if before is None:
            return None
        earlier = np.asarray(before['angstrom'], dtype=float)
        if earlier.shape != self.bohr.shape:
            return None
        moved = ((self.bohr - earlier / BOHR).reshape(-1)) * self.scale
        basis = self._basis()
        inside = basis.T @ moved
        length = float(np.linalg.norm(inside))
        if length < 1e-6:
            return None
        inside /= length
        curvature, direction = np.linalg.eigh(basis.T @ self.hessian @ basis)
        against = np.abs(direction.T @ inside)
        chosen = int(np.argmax(against))
        self.following = basis @ direction[:, chosen]
        self.following /= np.linalg.norm(self.following)
        return float(against[chosen])

    def took(self, xyz_text: str, *, aimed_from: Any = None) -> Dict[str, Any]:
        """Start again from a structure the user has just made.

        The Hessian is recomputed rather than carried across, and that is a
        measured choice rather than a careful one: carried across a drag it
        still reaches the same saddle, but in 62 steps against 15, because a
        Bofill update repairs a Hessian one step at a time and a hand moves
        further in one gesture than a climb does in twenty.  Six tenths of a
        second once, on the mouse being let go, is cheaper than 47 gradients
        and a picture that wanders on the way.
        """
        found = _elements(xyz_text)
        if found is None:
            return {'ok': False, 'status': 'There is no structure to climb.'}
        if len(found['symbols']) != len(self.symbols):
            return {'ok': False,
                    'status': ('The structure changed while the climb was '
                               'running, so it starts again from this one.')}
        self.bohr = found['angstrom'] / BOHR
        self.energy = None
        self.gradient = None
        began = self.start(aimed_from=aimed_from)
        began['ok'] = True
        return began

    def _keep(self) -> None:
        """Remember where the climb stands, in case the next step is refused."""
        self._before = (self.bohr.copy(), self.energy, self.gradient.copy(),
                        self.hessian.copy(), self.following.copy()
                        if self.following is not None else None)

    def _back_out(self) -> bool:
        """Put the last step back, and shorten the next one.

        Both reasons to do this come from the same place.  A step that lost the
        mode it was climbing and a step the quadratic model got badly wrong are
        both steps that went further than the model was good for, and what
        every implementation does about either is undo it and halve the trust
        radius rather than carry on from where it landed.
        """
        if getattr(self, '_before', None) is None:
            return False
        if self.trust <= SMALLEST_TRUST:
            return False               # already as careful as it can be
        (self.bohr, self.energy, self.gradient, self.hessian,
         self.following) = self._before
        self.bohr = self.bohr.copy()
        self.gradient = self.gradient.copy()
        self.hessian = self.hessian.copy()
        self.trust = max(SMALLEST_TRUST, 0.5 * self.trust)
        return True

    def step(self) -> Dict[str, Any]:
        """One gradient, one P-RFO step, one Bofill update.

        Measured at about 10 ms for sixteen atoms with the in-process engine,
        of which 6 are the gradient -- a hundred a second, so the picture and
        not the calculation is what limits how fast this can be shown.

        A step can also come back refused, having moved nothing: see
        :meth:`_back_out`.  It still cost a gradient, and the caller should
        simply ask for another.
        """
        if self.hessian is None:
            self.start()
        if self.energy is None or self.gradient is None:
            self.energy, self.gradient = self._measure(self.bohr)
        basis = self._basis()
        inside = basis.T @ self.hessian @ basis
        forces = basis.T @ (self.gradient / self.scale)
        aimed = None
        if self.following is not None:
            aimed = basis.T @ self.following
            length = float(np.linalg.norm(aimed))
            aimed = aimed / length if length > 1e-8 else None
        found = partitioned_step(inside, forces, aimed)
        lost = (found['overlap'] is not None
                and found['overlap'] < LEAST_OVERLAP)
        if lost and self._back_out():
            self.refused += 1
            return {'energy': self.energy, 'gmax': float(np.abs(
                        self.gradient).max()), 'grms': None, 'moved': 0.0,
                    'curvature': found['curvature'],
                    'overlap': found['overlap'],
                    'negative': found['negative'], 'trust': self.trust,
                    'ratio': None, 'steps': self.steps, 'refused': 'the mode',
                    'converged': False}
        self._keep()
        self.following = basis @ found['mode']
        self.following /= np.linalg.norm(self.following)
        weighted = basis @ found['step']
        move = weighted / self.scale
        reach = float(np.linalg.norm(move))
        shortened = reach > self.trust
        if shortened:
            move = move * (self.trust / reach)
            weighted = move * self.scale
        predicted = float(forces @ (basis.T @ weighted)
                          + 0.5 * (basis.T @ weighted) @ inside
                          @ (basis.T @ weighted))
        moved_to = self.bohr + move.reshape(-1, 3)
        energy, gradient = self._measure(moved_to)
        happened = energy - self.energy
        # Fletcher's ratio, and the trust radius it drives.  A saddle search
        # goes uphill along one mode, so neither the actual nor the predicted
        # change has a sign to expect: what is being asked is only whether the
        # quadratic model told the truth about this step.
        honest = happened / predicted if abs(predicted) > 1e-14 else 1.0
        if (abs(happened) > RATIO_MATTERS_ABOVE
                and not WORST_RATIO <= honest <= BEST_RATIO):
            if self._back_out():
                self.refused += 1
                return {'energy': self.energy, 'grms': None, 'moved': 0.0,
                        'gmax': float(np.abs(self.gradient).max()),
                        'curvature': found['curvature'],
                        'overlap': found['overlap'],
                        'negative': found['negative'], 'trust': self.trust,
                        'ratio': honest, 'steps': self.steps,
                        'refused': 'the energy', 'converged': False}
        if honest < 0.25:
            self.trust = max(SMALLEST_TRUST, 0.5 * min(reach, self.trust))
        elif honest > 0.75 and shortened:
            self.trust = min(LARGEST_TRUST, 2.0 * self.trust)
        self.hessian = bofill_update(
            self.hessian, weighted, (gradient - self.gradient) / self.scale)
        self.bohr, self.energy, self.gradient = moved_to, energy, gradient
        self.steps += 1
        after = self._basis().T @ gradient
        largest = float(np.abs(gradient).max())
        typical = float(np.sqrt(np.mean(after ** 2)))
        return {'energy': energy, 'gmax': largest, 'grms': typical,
                'moved': float(np.abs(move).max()) * BOHR,
                'curvature': found['curvature'], 'overlap': found['overlap'],
                'negative': found['negative'], 'trust': self.trust,
                'ratio': honest, 'steps': self.steps, 'refused': None,
                'converged': largest < GRADIENT_MAX and typical < GRADIENT_RMS}

    # -- what it reached --------------------------------------------------

    def frequencies(self, exact: bool = False) -> np.ndarray:
        """The modes in cm-1, translations and rotations taken out.

        *exact* recomputes the Hessian, which costs what the starting one cost.
        The updated Hessian is free and is right about how many modes go the
        wrong way; it is only roughly right about the number itself, so the
        verdict at the end of a climb is worth the recomputation and the ones
        along the way are not.
        """
        raw = (self.numerical_hessian() if exact or self.hessian is None
               else self.hessian)
        # Out of the working metric and into the frequency one.
        cartesian = raw * np.outer(self.scale, self.scale)
        weighted = cartesian / np.outer(self.weight, self.weight)
        basis = translations_and_rotations(self.masses, self.bohr)
        value = np.linalg.eigvalsh(basis.T @ weighted @ basis)
        return np.sort(np.sign(value) * np.sqrt(np.abs(value)) * HARTREE_IN_CM)

    def verdict(self, exact: bool = True) -> Dict[str, Any]:
        """Whether what was reached is a transition state, and what it is if not.

        Checked and not assumed, because a saddle search that converges has
        only proved that the gradient vanished.  Measured on this box, a path
        finder feeding ORCA's OptTS converged onto a *second-order* saddle at
        -73 and -27 cm-1 and reported success; and measured here, a climb from
        a van-der-Waals complex with no hand to guide it converges in 105 steps
        onto a structure with the fragments 3.5 A apart and not one mode going
        the wrong way.  Both are stationary points and neither is a reaction.

        xtb's own ``--hess`` is asked when there is an xtb to ask, and it is
        worth the second it costs.  A Hessian differenced from *gradients* --
        which is what this file otherwise has -- inherits whatever
        inconsistency the model's gradient has with the model's energy, and on
        flat, symmetric structures GFN2's is large enough to change the answer.
        Measured on planar ammonia at its own D3h stationary point, in separate
        xtb processes so nothing is shared: the gradient says one Hessian
        element is -0.01206 and four energies of the same xtb say +0.02958, and
        ``xtb --hess`` agrees with the energies (+0.0296).  The umbrella comes
        out at +386 cm-1 from gradients and -971 from either of the other two,
        which is the difference between calling flat ammonia a minimum and
        calling it the inversion transition state that it is.

        On an ordinary structure the two agree and the choice does not matter:
        on the Diels-Alder saddle, -394.6 cm-1 from gradients against -394.8
        from ``xtb --hess``.  It is symmetry that separates them, and a
        transition state is quite often symmetric.
        """
        modes = self.frequencies_from_xtb()
        if modes is None:
            modes = self.frequencies(exact=exact)
        wrong = [float(one) for one in modes if one < IMAGINARY_BELOW]
        return {'count': len(wrong), 'modes': wrong[:4],
                'lowest': float(modes[0]) if modes.size else None,
                'ok': len(wrong) == 1}

    def frequencies_from_xtb(self) -> Optional[np.ndarray]:
        """xtb's own Hessian for the structure as it now stands, or ``None``.

        A second at sixteen atoms, so this is for the end of a climb and not
        for the middle of one.  ``None`` when there is no xtb binary to run or
        it did not finish, and the caller falls back on the Hessian it already
        has -- an approximate verdict is worth more than none, as long as
        nothing pretends the two are the same thing.
        """
        binary = _gfn.find_xtb()
        if binary is None:
            return None
        flag = {'gfn2': ['--gfn', '2'], 'gfn1': ['--gfn', '1'],
                'gfnff': ['--gfnff']}.get(self.method)
        if flag is None:
            return None
        folder = Path(tempfile.mkdtemp(prefix='delfin-climb-hess-'))
        try:
            (folder / 'in.xyz').write_text(
                xyz_document(self.symbols, self.angstrom,
                             'from the DELFIN viewer'), encoding='utf-8')
            room = dict(os.environ)
            room['OMP_NUM_THREADS'] = room['MKL_NUM_THREADS'] = str(self.cores)
            done = subprocess.run(
                [binary, 'in.xyz', *flag, '--hess', '--chrg',
                 str(self.charge), '--uhf', str(self.uhf)]
                + (['--alpb', str(self.solvent)] if self.solvent else []),
                cwd=str(folder), env=room, capture_output=True, text=True,
                timeout=300)
        except (OSError, subprocess.SubprocessError):
            return None
        finally:
            shutil.rmtree(folder, ignore_errors=True)
        # The projected frequencies, printed once as a block of "eigval" rows.
        # The six that are nothing at all come out within a cm-1 of zero and
        # are dropped by count rather than by size: a real mode can be that
        # soft, and throwing away everything small would throw those away too.
        seen: List[float] = []
        for line in (done.stdout or '').splitlines():
            if line.strip().startswith('eigval'):
                for word in line.split(':', 1)[1].split():
                    try:
                        seen.append(float(word))
                    except ValueError:
                        pass
        wanted = 3 * len(self.symbols)
        if len(seen) < wanted:
            return None
        seen = sorted(seen[:wanted])
        # Translations and rotations are the ones nearest zero, not the lowest.
        spare = wanted - self._basis().shape[1]
        by_size = sorted(range(wanted), key=lambda i: abs(seen[i]))
        drop = set(by_size[:spare])
        return np.array([v for i, v in enumerate(seen) if i not in drop])


def climb_to_saddle(xyz_text: str, method: str = 'gfn2', *,
                    charge: int = 0, uhf: int = 0,
                    solvent: Optional[str] = None,
                    aimed_from: Any = None,
                    max_steps: int = 150, cores: int = 4,
                    on_frame: Any = None, should_stop: Any = None,
                    ) -> Dict[str, Any]:
    """Climb until it converges or runs out of steps, reporting as it goes.

    The whole climb without the loop, for a caller that only wants the answer:
    the editor drives :class:`Climb` a step at a time so it can answer a mouse
    between steps, and everything else can call this.  *on_frame* is handed
    the structure after every step, and *should_stop* is asked before each --
    the same two hooks :func:`saddle.optimise_to_saddle` streams ORCA through,
    so both saddle searches look the same from the outside.
    """
    if str(method or '').lower() not in CLIMB_METHODS:
        return {'ok': False,
                'status': (f'An interactive climb runs on xtb, and {method} '
                           f'is not one of '
                           f'{", ".join(sorted(CLIMB_METHODS))}.')}
    if _elements(xyz_text) is None:
        return {'ok': False, 'status': 'There is no structure to climb from.'}
    if not have_fast_gradients() and _gfn.find_xtb() is None:
        return {'ok': False,
                'status': ('An interactive climb needs xtb, which was not '
                           'found.')}
    began = time.perf_counter()
    walk = Climb(xyz_text, method, charge=charge, uhf=uhf, solvent=solvent,
                 cores=cores)
    try:
        opened = walk.start(aimed_from=aimed_from)
        arrived = False
        for _ in range(max(1, int(max_steps))):
            if should_stop is not None and should_stop():
                break
            outcome = walk.step()
            if on_frame is not None:
                on_frame(walk.frame(), outcome)
            if outcome.get('converged'):
                arrived = True
                break
        shape = walk.verdict()
        seconds = time.perf_counter() - began
        return {'ok': arrived, 'xyz': walk.xyz(
                    'Climbed to a transition state' if arrived
                    else 'Where the climb got to'),
                'seconds': seconds, 'steps': walk.steps,
                'started': opened, 'imaginary': shape,
                'gradients': int(getattr(walk.engine, 'calls', 0)),
                'status': (
                    f'The climb converged in {walk.steps} steps '
                    f'({seconds:.1f} s).' if arrived else
                    f'The climb did not converge in {walk.steps} steps; the '
                    f'structure it reached is shown.')}
    finally:
        walk.close()
