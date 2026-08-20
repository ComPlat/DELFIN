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
(*J. Comput. Chem.* **15**, 1, 1994), with the exact Hessian computed again
every twenty steps.  The step maximises the energy along one chosen mode and
minimises it along every other, and between recomputations the Hessian is
updated from the step.  It is what Gaussian's Berny and ORCA's OptTS do -- the
recomputation is ORCA's ``Recalc_Hess``, which :mod:`saddle` sets to 5 for the
button next door -- and the gradient-only form of it is what pysisyphus and
Sella do; their sources were read for the step control and the mode tracking
rather than derived here.

Two things about it were wrong until they were measured against twenty-one
drags a hand could actually make -- partly formed bonds, rings pulled open,
protons half transferred, on fifteen to fifty atoms -- rather than against the
one tidy Diels-Alder estimate the file was written for.  Both are recorded at
:data:`AIM_WITHIN` and :data:`HESSIAN_EVERY`, and both were the difference
between a climb that arrives and one that walks:

* the hand's direction was matched against the *whole* spectrum, so dragging a
  proton across a hydrogen bond aimed the climb at the O-H stretch at 1360
  cm-1 and P-RFO went uphill along it; and
* the Bofill update does not merely drift, it **invents curvature that is not
  there** -- one to two negative eigenvalues the exact Hessian at the same
  geometry does not have, and once at 1400 times its size.

With both put right, on those twenty-one drags: the climb reaches the reaction
the hand pointed at **9 times in 21, in a median of 5.7 s**, against **6 in 21**
before and against **9 in 21 in a median of 50.3 s** for ORCA's own OptTS
started from the same geometries.  "Reaches the reaction" is not "converged":
both routes also arrive, confidently, at methyl torsions at -48 cm-1 and at
fragments rocking at -52, so what is counted is whether xtb's own Hessian on
what came back has one imaginary mode *and that mode is the stretch of the pair
of atoms the hand was holding*.

And then the thing that mattered more than either number: **they are wrong
about different drags**.  Nine each, never the same nine, because the climb
sets off along the eigenvector nearest the gesture and ORCA sets off along the
lowest mode of its own Hessian -- one is listening to the hand and the other to
the surface.  A third search, this same climb told to forget the gesture, is a
third nine again.  :func:`reach_the_reaction` is what that is worth: the three
tried cheapest first, each one's answer checked against the contact the hand
held, **13 of 16** where any one of them alone is 7 to 9.  (Sixteen and not
twenty-one: five of the drags were built by a constrained relaxation that did
not finish and left atoms off on their own, so they are not gestures anybody
can make.  Both denominators are carried there.)  What is left after that is
not an optimiser's failure -- on two of the three it misses, the reaction
coordinate is not among the five softest modes of the structure the hand left
at all.  That is the honest state of a saddle search begun by hand, and it is
worth saying rather than implying otherwise.

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
* and with the exact Hessian recomputed every twenty steps, a step costs on
  average four times its bare price at sixteen atoms and eight times at fifty,
  because one Hessian is 6N gradients: 0.55 s at sixteen atoms, 1.10 at
  twenty, 4.92 at thirty-three and 16.3 at fifty, against 9.1, 14.9, 51.2 and
  116.5 ms for a step.  So the climb is interactive in the sense the word is
  used here -- tens of milliseconds a step -- up to about twenty atoms, and
  above that it is a short wait rather than an animation

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
from . import saddle as _saddle

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

#: How many of the softest modes the hand is allowed to have meant.
#:
#: A drag is answered by the Hessian eigenvector that most resembles it, and
#: for a long time that search ran over the whole spectrum.  Measured on
#: twenty-one drags, that is where most of them died: dragging a hydroxyl
#: proton across a hydrogen bond in a glycylglycine, the eigenvector nearest
#: the drag is the **O-H stretch at 1360 cm-1**, mode 31 of 45.  P-RFO then
#: maximises the energy along it, which is not a saddle search -- it is
#: pulling the proton off.  That climb passed +30 000 kcal/mol by its eleventh
#: step and never came back.  Eleven of the twenty-one aimed at a mode above
#: 170 cm-1, and not one of those reached the reaction the hand meant.
#:
#: So the hand may only name one of the softest few, which is where every
#: implementation of Baker's mode following looks: geomeTRIC hard-codes the
#: climbed mode as the lowest, optking refuses anything else outright, and
#: MOPAC's default for a transition state is mode 1.  The reason is physical
#: rather than numerical -- at a first-order saddle the reaction coordinate is
#: the mode with the negative eigenvalue, so on the way to one it is soft, and
#: a 1360 cm-1 stretch cannot become it without the molecule coming apart.
#:
#: Five, and the number is measured rather than chosen: over the same
#: twenty-one drags the count that reached the reaction the hand pointed at
#: was 7 at a window of one, 8 at three, **9 at five** and 8 at eight.  One is
#: too narrow because it throws the hand's information away -- it is the
#: lowest mode whatever was dragged, which is what the climb did before it was
#: told about the hand at all.  Wider than five and the hard modes start
#: coming back.  The overlap the hand achieved is reported either way, so a
#: gesture that named nothing in particular can be said out loud.
AIM_WITHIN = 5

#: How often the Hessian is computed again rather than updated, in steps.
#:
#: This is the expensive half of the fix and it is ORCA's own answer --
#: ``Recalc_Hess``, which :mod:`saddle` sets to 5 for the button next door.
#: pysisyphus says it plainly in its own documentation: "when the Hessian for
#: the chosen computational method is reasonably cheap it is a good idea to
#: recalculate it periodically; between recalculations it's updated using the
#: Bofill update", and "use as many exact Hessians as your computational
#: budget allows".
#:
#: What it is defending against was measured here rather than taken on trust.
#: A Bofill-updated Hessian does not merely drift: it **invents curvature that
#: is not there**.  Against a freshly computed Hessian at the same geometry,
#: every ten steps, on twenty-one climbs: 20 to 25 per cent relative Frobenius
#: distance within ten steps even on the climbs that succeed, the curvature of
#: the mode being climbed wrong by a factor of two, and one to two negative
#: eigenvalues that the exact Hessian does not have.  On a salicylaldehyde
#: proton transfer the exact Hessian has *no* negative eigenvalue from step 10
#: onwards and the carried one insists on one or two for three hundred steps,
#: so the climb spends them chasing a mode that does not exist.  On the
#: glycylglycine above it reached 1400 times the size of the real Hessian.
#: Baker and Chan named exactly this in 1996 (*J. Comput. Chem.* **17**, 888):
#: "for Cartesian coordinates, the commonly used Hessian update schemes are
#: unable to guarantee preservation of the necessary transition state
#: eigenvalue structure".
#:
#: Twenty steps, and what it costs is not hidden: one Hessian is 6N gradients,
#: measured on this box at 0.55 s for sixteen atoms, 1.10 for twenty, 4.92 for
#: thirty-three and 16.3 for fifty -- 61, 74, 96 and 140 steps' worth.  So a
#: refresh every twenty steps makes the climb about four times its bare cost
#: at sixteen atoms and eight times at fifty.  That is the price of the
#: eigenvalue structure, and it is the same price ORCA pays.
HESSIAN_EVERY = 20

#: And how many of those one climb may spend, so that a climb going nowhere
#: cannot spend for ever.  Measured: every climb that reached the reaction the
#: hand pointed at used four or fewer, the longest being 84 steps -- so this
#: bound never touches a climb that is arriving, only one that is not.
#:
#: It is also what bounds how the climb *looks*, which is the other half of
#: the price.  A refresh is a step that takes as long as a Hessian, and the
#: loop that draws the frames and watches for a Stop is between steps, so a
#: refresh is a pause in the animation and a Stop that waits.  A climb that
#: arrives now takes twenty to thirty steps, so it pauses once; one that is
#: lost pauses eight times and no more.
MOST_HESSIANS = 8

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

#: How much of the imaginary mode has to be the contact the hand held before
#: what was reached counts as the reaction the user pointed at.
#:
#: Counting imaginary modes is not enough and that is measured, not feared.
#: Over the twenty-one drags both this climb and ORCA's OptTS converge,
#: confidently and with exactly one mode going the wrong way, onto methyl
#: torsions at -48 cm-1 and onto two fragments rocking against each other at
#: -52.  Every one of those passes :func:`saddle.verdict` and none of them is
#: the reaction anybody asked for.  So what is asked is what the mode *does*:
#: the normalised stretch of the held pair, dotted into the normalised mode.
#:
#: A fifth, and the gap it sits in is wide rather than narrow.  Measured on
#: the twenty-one: the structures that are the reaction score 0.40 to 0.95 --
#: 0.68 on either forming bond of the Diels-Alder saddle, 0.63 on the
#: anthracene adduct, 0.51 on the salicylaldehyde proton transfer -- and the
#: torsions and rockings both routes wrongly converge onto score 0.00 to 0.07.
#: Nothing on either side of the twenty-one lands between 0.07 and 0.40, so
#: any threshold in that range gives the same answer and this one is the
#: middle of it.
IS_THE_REACTION = 0.20

#: How many steps a climb is given before it is taken to have gone wrong.
#:
#: Measured on the twenty-one drags: every climb that reached the reaction
#: the hand pointed at arrived in 22 to 65 steps, and no climb that ran
#: longer than that ever reached it -- the nine that did not converge ran on
#: to between 202 and 361 steps and were wrong at the end of all of them.
#: A hundred is half again as long as the longest arrival, so it cannot cut a
#: climb that is getting there, and it is what stops the other nine spending
#: three hundred steps to be wrong more slowly.
#:
#: It is a bound on the wait and not on the answer.  The fallback is handed
#: the structure the *hand* made rather than the one the climb reached -- see
#: :func:`reach_the_reaction` for the measurement that settled that -- so
#: where a doomed climb stops changes nothing about what comes back.
CLIMB_CEILING = 100

#: And how many optimiser cycles the last rung is allowed.
#:
#: Sixty, which is what :func:`saddle.optimise_to_saddle` already defaults to,
#: and it is written down here because raising it was tried and measured to be
#: a loss rather than because nobody looked.
#:
#: It looked like an obvious win.  Two of the nine reactions ORCA reaches over
#: the twenty-one drags needed 80 and 120 cycles, so at 60 they are lost --
#: and both of those two turn out to be among the five drags whose structure
#: fell apart while the bench was building it.  Over the sixteen whole ones
#: every reaction ORCA reaches converges in **13 to 49 cycles**, all of them
#: under 60, and 60 and 120 both give it 7 of 16.  The higher ceiling buys
#: nothing on a structure a hand can actually make.
#:
#: What it costs is not nothing, because an OptTS that is not going to arrive
#: runs to the ceiling by definition.  Measured over the same twenty-one, the
#: presses where ORCA's answer is no take a median of **63 s at 120 cycles and
#: 25 s at 60** -- two and a half times the wait, on exactly the presses where
#: there is nothing at the end of it.  So the ceiling stays where it was, and
#: it is pinned here rather than inherited so that the rung above knows what
#: it is asking for.
FALLBACK_STEPS = 60

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
                     following: Optional[np.ndarray],
                     within: int = AIM_WITHIN) -> Dict[str, Any]:
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

    *within* bounds that search to the softest few, for the reason
    :data:`AIM_WITHIN` gives -- and it is needed here and not only at the
    start, because a Bofill update can grow an eigenvalue no surface has.
    Measured on a glycylglycine, a carried curvature of -1379 where the exact
    Hessian at the same geometry says +0.34: once the followed vector has
    rotated into one of those, the overlap test reports 1.00 for ever and
    nothing brings the climb back.
    """
    curvature, direction = np.linalg.eigh(hessian)
    force = direction.T @ gradient
    if following is None:
        chosen = 0
        overlap = None
    else:
        against = np.abs(direction.T @ following)
        room = min(int(within), against.size) if within else against.size
        chosen = int(np.argmax(against[:room]))
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
        #: Exact Hessians spent since :meth:`start`, and which step the last
        #: one belongs to -- a refused step does not advance ``steps``, so
        #: without the second of these the same step could pay for two.
        self.hessians = 0
        self._exact_at = -1
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
        self.hessians = 0
        self._exact_at = -1
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

        Only the softest :data:`AIM_WITHIN` modes are candidates, and that one
        line is what most of the twenty-one measured drags turned on: over the
        whole spectrum the nearest eigenvector to a proton being dragged
        across a hydrogen bond is the O-H stretch, and climbing a stretch
        tears the molecule apart rather than finding a saddle.

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
        room = min(int(AIM_WITHIN), against.size) if AIM_WITHIN else against.size
        chosen = int(np.argmax(against[:room]))
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

        Every :data:`HESSIAN_EVERY` steps the Hessian is computed again
        instead of updated, and that step costs 6N gradients rather than one:
        0.55 s at sixteen atoms, 16.3 at fifty.  It is the expensive half of
        what makes the climb reach the reaction the hand pointed at rather
        than some other stationary point, and :data:`HESSIAN_EVERY` says what
        was measured about both halves of that trade.

        A step can also come back refused, having moved nothing: see
        :meth:`_back_out`.  It still cost a gradient, and the caller should
        simply ask for another.
        """
        if self.hessian is None:
            self.start()
        if (HESSIAN_EVERY and self.steps and self.hessians < MOST_HESSIANS
                and self.steps % int(HESSIAN_EVERY) == 0
                and self._exact_at != self.steps):
            # A refresh, not a restart: the trust radius, the mode being
            # followed and everything the climb has learned about where it is
            # going are kept.  Only the matrix that had stopped describing the
            # surface is replaced.
            self.hessian = self.numerical_hessian()
            self._exact_at = self.steps
            self.hessians += 1
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

    def held_stretch(self, held: Sequence[int]) -> Optional[np.ndarray]:
        """"Make this contact longer", as a unit vector a mode can be compared to.

        The Wilson B row of the distance between the two atoms the hand was
        holding, carried into the mass-weighted metric the modes live in and
        normalised, so that dotting it into a normalised mode gives a number
        between nought and one.  A pure two-atom stretch would be one; a mode
        in which that pair only comes along for the ride is a few hundredths.
        """
        try:
            one, two = (int(held[0]), int(held[1]))
        except (TypeError, ValueError, IndexError):
            return None
        count = self.bohr.shape[0]
        if not (0 <= one < count and 0 <= two < count) or one == two:
            return None
        arm = self.bohr[one] - self.bohr[two]
        length = float(np.linalg.norm(arm))
        if length < 1e-8:
            return None
        row = np.zeros((count, 3))
        row[one] = arm / length
        row[two] = -arm / length
        # A B row differentiates with respect to Cartesians, so it goes into
        # the mass-weighted metric divided by the root mass rather than
        # multiplied by it.
        got = row.reshape(-1) / self.weight
        return got / np.linalg.norm(got)

    def verdict(self, exact: bool = True,
                held: Optional[Sequence[int]] = None) -> Dict[str, Any]:
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

        *held* is the pair of atoms the hand was dragging, and given one this
        answers the harder question as well: not "is it a saddle" but "is it
        *the* saddle the hand pointed at".  ``share`` is how much of the
        imaginary mode is that contact stretching, and ``reaction`` is whether
        the structure passes both tests.  See :data:`IS_THE_REACTION` for what
        that separates and why counting imaginary modes alone does not.
        """
        got = self.modes_from_xtb()
        if got is None:
            modes, shapes = self.frequencies(exact=exact), None
        else:
            modes, shapes = got['cm'], got['shape']
        wrong = [float(one) for one in modes if one < IMAGINARY_BELOW]
        first = len(wrong) == 1
        share = None
        if held is not None and shapes is not None and modes.size:
            unit = self.held_stretch(held)
            if unit is not None:
                share = abs(float(unit @ shapes[:, 0]))
        return {'count': len(wrong), 'modes': wrong[:4],
                'lowest': float(modes[0]) if modes.size else None,
                'ok': first, 'share': share,
                'reaction': (None if share is None
                             else bool(first and share >= IS_THE_REACTION))}

    def _hessian_from_xtb(self) -> Optional[np.ndarray]:
        """xtb's own second derivatives for where the structure stands.

        Cartesian, unweighted, Hartree per Bohr squared, or ``None`` when
        there is no xtb binary to run or it did not finish.  ``xtb --hess``
        writes the matrix to a file called ``hessian`` beside its input, and
        that file is read rather than the frequencies it prints, because the
        matrix carries the *shapes* of the modes and the printed list carries
        only their numbers -- and which mode was reached is the question
        :meth:`verdict` has to answer.  It is the same one subprocess either
        way: measured over the twenty-one drags, 0.11 s at ten atoms, 0.4 at
        sixteen, 1.0 at twenty, 3.7 at thirty-three and 11.8 at fifty, and
        reading the file rather than the printout changed neither the cost
        nor the answer -- the two disagree by at most 0.11 cm-1 on any mode
        of any of the twenty-one.
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
            subprocess.run(
                [binary, 'in.xyz', *flag, '--hess', '--chrg',
                 str(self.charge), '--uhf', str(self.uhf)]
                + (['--alpb', str(self.solvent)] if self.solvent else []),
                cwd=str(folder), env=room, capture_output=True, text=True,
                timeout=300)
            written = folder / 'hessian'
            if not written.is_file():
                return None
            words: List[float] = []
            for line in written.read_text(
                    encoding='utf-8', errors='ignore').splitlines():
                if line.strip().startswith('$'):
                    continue
                for word in line.split():
                    try:
                        words.append(float(word))
                    except ValueError:
                        return None
        except (OSError, subprocess.SubprocessError):
            return None
        finally:
            shutil.rmtree(folder, ignore_errors=True)
        size = 3 * len(self.symbols)
        if len(words) < size * size:
            return None
        raw = np.array(words[:size * size]).reshape(size, size)
        return 0.5 * (raw + raw.T)

    def modes_from_xtb(self) -> Optional[Dict[str, Any]]:
        """The modes xtb's own Hessian has here: frequencies and shapes.

        ``{'cm': ..., 'shape': ...}``, the frequencies in cm-1 sorted from
        softest, and the shapes as the columns of an array in the
        mass-weighted Cartesian metric -- translations and rotations projected
        out, which is why there are ``3N - 6`` of them and not ``3N``.  A
        second at sixteen atoms, so this is for the end of a climb and not for
        the middle of one.  ``None`` when xtb could not be asked, and the
        caller falls back on the Hessian it already has -- an approximate
        verdict is worth more than none, as long as nothing pretends the two
        are the same thing.
        """
        raw = self._hessian_from_xtb()
        if raw is None:
            return None
        basis = translations_and_rotations(self.masses, self.bohr)
        weighted = raw / np.outer(self.weight, self.weight)
        value, vector = np.linalg.eigh(basis.T @ weighted @ basis)
        order = np.argsort(value)
        value = value[order]
        return {'cm': np.sign(value) * np.sqrt(np.abs(value)) * HARTREE_IN_CM,
                'shape': basis @ vector[:, order]}

    def frequencies_from_xtb(self) -> Optional[np.ndarray]:
        """Just the frequencies of :meth:`modes_from_xtb`, or ``None``."""
        got = self.modes_from_xtb()
        return None if got is None else got['cm']


def climb_to_saddle(xyz_text: str, method: str = 'gfn2', *,
                    charge: int = 0, uhf: int = 0,
                    solvent: Optional[str] = None,
                    aimed_from: Any = None,
                    max_steps: int = 150, cores: int = 4,
                    on_frame: Any = None, should_stop: Any = None,
                    held: Optional[Sequence[int]] = None,
                    ) -> Dict[str, Any]:
    """Climb until it converges or runs out of steps, reporting as it goes.

    The whole climb without the loop, for a caller that only wants the answer:
    the editor drives :class:`Climb` a step at a time so it can answer a mouse
    between steps, and everything else can call this.  *on_frame* is handed
    the structure after every step, and *should_stop* is asked before each.

    *held* is the pair of atoms the hand was dragging.  Given one, the verdict
    at the end says not only whether a saddle was reached but whether it is
    the reaction that pair names -- which is a different question, and the one
    :func:`reach_the_reaction` is built on.
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
        shape = walk.verdict(held=held)
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


def _contact(symbols: Sequence[str], held: Optional[Sequence[int]]) -> str:
    """The pair the hand held, written the way the viewer numbers atoms.

    Zero-based, because that is what the ``#`` button draws on the atoms --
    ``addLabel(String(i))`` over ``selectedAtoms``, counting from nought -- and
    a sentence that numbers them differently from the picture is worse than
    one that does not number them at all.
    """
    if held is None:
        return 'the contact you dragged'
    one, two = int(held[0]), int(held[1])
    return f'{symbols[one]}{one}-{symbols[two]}{two}'


def _named_landing(shape: Optional[Dict[str, Any]]) -> str:
    """What a Hessian found, in a phrase that can be put inside a sentence."""
    if not shape:
        return 'what it reached could not be checked'
    count = int(shape.get('count') or 0)
    modes = shape.get('modes') or []
    if count == 0:
        return 'it converged onto a minimum'
    if count == 1:
        return ('it converged onto a different saddle, at '
                f'{modes[0]:.0f} cm-1')
    return (f'it converged onto a saddle of order {count} ('
            + ', '.join(f'{one:.0f}' for one in modes[:3]) + ' cm-1)')


def _why_not(shape: Dict[str, Any], contact: str) -> str:
    """Why the climb's answer was not taken, said in one sentence."""
    count = int(shape.get('count') or 0)
    modes = shape.get('modes') or []
    if count == 0:
        return ('The climb settled somewhere with no mode going the wrong '
                'way, so it is a minimum and not a transition state.')
    if count > 1:
        return (f'The climb settled on a saddle of order {count} ('
                + ', '.join(f'{one:.0f}' for one in modes[:3])
                + ' cm-1), and a transition state has exactly one mode going '
                  'the wrong way.')
    if shape.get('share') is None:
        return (f'The climb reached a saddle at {modes[0]:.0f} cm-1, but '
                'nothing named the contact you were dragging, so there is no '
                'way to tell whether it is the reaction you were pointing at.')
    return (f'The climb reached a saddle at {modes[0]:.0f} cm-1, but its '
            f'imaginary mode is only {shape["share"]:.2f} of the {contact} '
            'stretch -- it is a different reaction from the one you pointed '
            'at.')


def _checked(xyz_text: Optional[str], method: str, *,
             held: Optional[Sequence[int]], charge: int, uhf: int,
             solvent: Optional[str], cores: int) -> Optional[Dict[str, Any]]:
    """xtb's own verdict on a structure some other optimiser handed back."""
    if not xyz_text or _elements(xyz_text) is None:
        return None
    walk = Climb(xyz_text, method, charge=charge, uhf=uhf, solvent=solvent,
                 cores=cores)
    try:
        return walk.verdict(held=held)
    finally:
        walk.close()


def reach_the_reaction(xyz_text: str, method: str = 'gfn2', *,
                       held: Optional[Sequence[int]] = None,
                       charge: int = 0, uhf: int = 0,
                       solvent: Optional[str] = None,
                       aimed_from: Any = None,
                       max_steps: int = CLIMB_CEILING,
                       cores: int = 4,
                       on_frame: Any = None, on_path: Any = None,
                       on_route: Any = None, should_stop: Any = None,
                       fallback: bool = True,
                       fallback_steps: int = FALLBACK_STEPS,
                       fallback_timeout: Optional[float] = 600.0,
                       ) -> Dict[str, Any]:
    """One press, three tries: climb, check what it reached, try again if wrong.

    A saddle search started by hand does not fail -- it succeeds at arriving
    somewhere, and usually somewhere else.  Measured over hand drags on
    fifteen to fifty atoms, three searches that all start from the structure
    the hand left reach the reaction the hand pointed at **about nine times
    each, and never the same nine**:

    ==================================  ========  ========  ===============
    where it starts climbing             of 16     of 21     told by
    ==================================  ========  ========  ===============
    the eigenvector nearest the drag        9         9      :meth:`aim`
    the surface's own softest mode          9         9      no *aimed_from*
    ORCA's OptTS, redundant internals       7         8      :mod:`saddle`
    **all three, in that order**          **13**    **14**
    ==================================  ========  ========  ===============

    Two columns because the bench has to be read twice.  Twenty-one drags were
    built by holding one contact and relaxing the rest under GFN2, and on five
    of them that relaxation did not finish and left atoms off on their own --
    an H2 floating beside a Claisen, four loose hydrogens on a Cope.  The
    editor's follow does not do that, so those five are not drags anybody can
    make and the sixteen are the honest denominator; the twenty-one are kept
    beside them because every earlier measurement in this file was taken on
    them.  Rebuilding the five and measuring again is the obvious next thing
    and has not been done.  It also matters more than a footnote: every one of
    the drags that made a higher optimiser ceiling look worth having was one
    of the torn five -- see :data:`FALLBACK_STEPS`.

    Either way the shape is the same and it is the whole reason this function
    exists: nine, nine and seven, and thirteen between them.  Closing that gap
    needs no new optimiser -- it needs the three that are already here, tried
    in order, with what each one reached actually checked.

    *Why* they differ is not a property of the molecules.  It is which
    direction each one sets off along, and each wins where *its* direction is
    already the reaction:

    * the aimed climb reached it on 9 of the 11 drags whose gesture is mostly
      one soft eigenvector (overlap 0.60 or better) and on **none of the 10**
      below that;
    * ORCA reached it on 8 of the 9 drags where the softest mode at the
      structure the hand left is already the held contact stretching, and on
      1 of the 12 where it is not.

    So on the Diels-Alder drags only the aimed climb gets, the softest mode is
    two fragments rocking -- 0.07 to 0.11 of the reaction -- and the gesture
    is 0.66 to 0.82 of it.  On the retro-Diels-Alder, Claisen and
    proton-transfer drags only ORCA gets, the softest mode is 0.80 to 0.94 of
    the reaction and the gesture is a muddle at 0.39 to 0.71.  One of them is
    listening to the hand and the other to the surface, and where those
    disagree is exactly where one beats the other.

    The middle rung is what that observation buys.  The same climb, told to
    forget the gesture and follow the softest mode, is ORCA's starting
    direction at the climb's price: it catches three drags the aimed climb
    misses -- the two retro-Diels-Alders and a Diels-Alder ORCA never reaches
    at all -- in 3.0, 4.9 and 3.9 s where ORCA needs 13.5, 20.6 and 76.2.  It
    is the same class, the same code and one more Hessian, and it is what
    takes the whole thing from 12 of the sixteen to 13.

    And the three of the sixteen that are never reached are not an optimiser's
    failure.  On two of them the reaction coordinate is not among the five
    softest modes at the structure the hand left at all (best overlap 0.00 and
    0.02), so nothing that follows a soft mode from there can find it; the
    third is a proton transfer whose drag geometry already has four imaginary
    modes.  Those are guesses that missed, and the sentence says so rather
    than implying that another minute of computing would help.

    Three things about the shape of this were measured rather than assumed.

    **The check is not optional.**  All three routes converge, confidently
    and with exactly one mode going the wrong way, onto methyl torsions at
    -48 cm-1 and onto fragments rocking at -52.  So what comes back is checked
    against the pair of atoms *held*, and only a first-order saddle whose
    imaginary mode is that contact stretching counts -- :data:`IS_THE_REACTION`
    for the numbers that separate.  Without *held* only the mode count can be
    checked; that is the weaker test, it cannot tell the rungs apart, and so
    the ladder stops at the first rung and the sentence says which test it was.

    **Every rung is handed the hand's own structure, never the last rung's.**
    That looks like throwing work away and it is measured to be right: over
    the twelve drags the aimed climb gets wrong, ORCA started from the hand's
    guess reaches the reaction on 5 and ORCA started from where the climb
    stopped reaches it on **2** -- handing the endpoint on loses the
    retro-Diels-Alder at 2.4 A, the Claisen at 2.0 and the salicylaldehyde
    proton transfer, all three of which the same optimiser reaches from the
    hand's structure.  A search that is going wrong does not stop somewhere
    useful.  So a rung that misses costs its wall time and nothing else, and
    :data:`CLIMB_CEILING` keeps that short.

    **This order and not the other.**  ORCA first, with the climb after it,
    reaches 9 of 21 rather than 14: started from where ORCA stopped the climb
    rescues **0 of the 12** ORCA gets wrong, because ORCA leaves it standing
    on a stationary point with no gradient to work with.  Started from the
    hand's structure instead it would reach the same fourteen -- it is the
    same union -- but it would pay ORCA's median 44 s on every press,
    including the nine the climb answers in 0.7 to 36.5 s.  Cheapest first is
    both the same answer and sooner.

    *on_frame* is handed the structure after every climb step, *on_path* the
    whole trajectory each time ORCA's grows -- the two hooks differ because
    the two optimisers report differently, and pretending otherwise would
    lose one of them.  *on_route* is called before each retry with the
    sentence explaining it, so a user is not left watching a still picture.
    *should_stop* is asked by all three.
    """
    if str(method or '').lower() not in CLIMB_METHODS:
        return {'ok': False, 'route': None,
                'status': (f'An interactive climb runs on xtb, and {method} '
                           f'is not one of '
                           f'{", ".join(sorted(CLIMB_METHODS))}.')}
    if _elements(xyz_text) is None:
        return {'ok': False, 'route': None,
                'status': 'There is no structure to climb from.'}
    contact = _contact(_elements(xyz_text)['symbols'], held)
    spent = 0.0
    tried: List[Dict[str, Any]] = []

    def _took(got: Dict[str, Any], route: str) -> Dict[str, Any]:
        """Book one attempt: name it, add its time to the press, keep it."""
        nonlocal spent
        spent += float(got.get('seconds') or 0.0)
        got['route'] = route
        got['seconds'] = spent
        tried.append(got)
        return got

    def _say(sentence: str) -> None:
        if on_route is None:
            return
        try:
            on_route(sentence)
        except Exception:                                   # noqa: BLE001
            # A page that cannot draw a sentence is not a reason to stop
            # searching; the structure at the end is the answer.
            pass

    def _arrived(got: Dict[str, Any], route: str, why: str) -> Dict[str, Any]:
        """The sentence for a rung that reached it."""
        shape = got.get('imaginary') or {}
        got['ok'] = True
        got['status'] = (
            why + f'{route} reached it in {got["seconds"]:.1f} s altogether: '
            f'one mode goes the wrong way, at '
            f'{(shape.get("modes") or [0.0])[0]:.0f} cm-1, and it is the '
            f'{contact} stretch.')
        return got

    # -- rung one: the eigenvector nearest the gesture ---------------------
    first = _took(climb_to_saddle(
        xyz_text, method, charge=charge, uhf=uhf, solvent=solvent,
        aimed_from=aimed_from, max_steps=max_steps, cores=cores,
        on_frame=on_frame, should_stop=should_stop, held=held), 'the climb')
    if not first.get('xyz'):
        return first
    reached = first.get('imaginary') or {}
    # Without a pair to check against, only the mode count can be asked, and
    # then the ladder is no use: every rung would be judged by the very test
    # this one has already passed, so the rest could only spend a minute to
    # agree.  The weaker test is used and the sentence says which one it was.
    told_which = reached.get('share') is not None
    if bool(reached.get('reaction') if told_which else reached.get('ok')):
        first['ok'] = True
        deepest = (reached.get('modes') or [0.0])[0]
        first['status'] = (
            (f'The climb reached the reaction you pointed at in '
             f'{first.get("steps", 0)} steps ({spent:.1f} s): one mode goes '
             f'the wrong way, at {deepest:.0f} cm-1, and it is the {contact} '
             f'stretch.')
            if told_which else
            (f'The climb reached a transition state in '
             f'{first.get("steps", 0)} steps ({spent:.1f} s): one mode goes '
             f'the wrong way, at {deepest:.0f} cm-1. Nothing said which '
             f'contact you meant, so whether it is the reaction you were '
             f'after is for you to read off the mode.'))
        return first
    why = _why_not(reached, contact)
    first['ok'] = False
    first['status'] = why
    if not told_which or not fallback:
        return first
    if should_stop is not None and should_stop():
        # Stopped between rungs.  ``ok`` is already false and the sentence is
        # already the one about where this rung landed, which is the honest
        # pair for a press the user ended: a climb that converged onto the
        # wrong saddle has not reached the reaction, whatever it converged on.
        return first

    # -- rung two: the surface's own softest mode --------------------------
    # Skipped when the first rung already was it -- with no gesture to aim
    # along, :meth:`aim` returns nothing and the climb follows the lowest
    # mode, so running it again would be the same climb twice.
    if aimed_from is not None:
        _say(why + ' Climbing again from your structure along the surface\'s '
                   'own softest mode instead of along your drag; that is '
                   'another few seconds.')
        again = _took(climb_to_saddle(
            xyz_text, method, charge=charge, uhf=uhf, solvent=solvent,
            aimed_from=None, max_steps=max_steps, cores=cores,
            on_frame=on_frame, should_stop=should_stop, held=held),
            'the softest mode')
        if (again.get('imaginary') or {}).get('reaction'):
            return _arrived(
                again, 'Climbing along the softest mode instead',
                why + ' ')
        again['ok'] = False
        again['status'] = (
            why + ' Climbing along the softest mode instead did not reach it '
            'either: ' + _named_landing(again.get('imaginary')) + '.')
        if should_stop is not None and should_stop():
            return again

    # -- rung three: ORCA's own optimiser ----------------------------------
    if not _saddle.find_orca():
        best = _closest(tried)
        best['ok'] = False
        best['status'] = (
            why + ' Climbing along the softest mode instead did not reach it '
            'either. The one optimiser left to try -- ORCA\'s OptTS, which '
            'walks in redundant internal coordinates rather than in '
            'Cartesians -- was not found, so this is where it stops.'
            if len(tried) > 1 else
            why + ' The optimiser that reaches some of these -- ORCA\'s '
            'OptTS -- was not found, so this is where it stops.')
        return best
    _say(why + ' Handing your structure to ORCA\'s optimiser, which walks in '
               'internal coordinates rather than Cartesians; that is a minute '
               'rather than a moment.')
    began = time.perf_counter()
    last = _saddle.optimise_to_saddle(
        xyz_text, method, charge=charge, uhf=uhf, solvent=solvent,
        max_steps=fallback_steps, cores=max(int(cores), 8),
        # ORCA's own confirming Hessian is turned off because the check below
        # is a better one and would otherwise be the second of two: the same
        # xtb Hessian, asked the extra question about the held contact.  0.4 s
        # at sixteen atoms and 11.8 at fifty, paid once instead of twice.
        timeout=fallback_timeout, confirm=False,
        on_frame=on_path, should_stop=should_stop)
    last['imaginary'] = _checked(last.get('xyz'), method, held=held,
                                 charge=charge, uhf=uhf, solvent=solvent,
                                 cores=cores)
    # And ORCA's own reading of what it reached goes with it.  It is taken
    # from the last Hessian block in the output, which belongs to a geometry
    # the run has since left -- measured on the sixteen-atom Diels-Alder, that
    # block says -301.55 cm-1 where a Hessian on the structure handed back
    # says -393.53.  Leaving both in the answer would be two verdicts on one
    # structure, and the wrong one is the one that reads like ORCA's.
    last.pop('verdict', None)
    # Timed here rather than taken from ORCA, because the check is part of
    # what the press cost and ORCA does not know about it.
    last['seconds'] = time.perf_counter() - began
    last = _took(last, 'ORCA')
    if (last.get('imaginary') or {}).get('reaction'):
        return _arrived(last, "ORCA's optimiser, given the structure you made,",
                        why + ' ')

    # Nothing reached it.  What comes back is whichever attempt got closest to
    # being a transition state at all, and the first keeps a tie because it is
    # the one the user watched.
    best = _closest(tried)
    best['ok'] = False
    best['status'] = (
        why + f' Climbing along the surface\'s own softest mode and then '
        f'ORCA\'s optimiser were given the same structure and did not reach '
        f'it either ({spent:.1f} s altogether). What is left on screen is the '
        f'nearest any of the three came: '
        + _named_landing(best.get('imaginary')) + '. All three start from '
        'where you let go, so it is the structure that has to change -- move '
        'the atoms nearer the arrangement the reaction passes through and '
        'press again.'
        if len(tried) > 2 else
        why + f' ORCA\'s optimiser was given the same structure and did not '
        f'reach it either ({spent:.1f} s altogether): '
        + _named_landing(best.get('imaginary')) + '. Both start from where '
        'you let go, so it is the structure that has to change -- move the '
        'atoms nearer the arrangement the reaction passes through and press '
        'again.')
    return best


def _closest(tried: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """Of the attempts that all missed, the one worth leaving on screen.

    A first-order saddle beats anything else even when it is the wrong
    reaction -- it is at least a structure the user can look at and move on
    from -- and the earliest attempt keeps a tie, because that is the one that
    was drawn while it walked.
    """
    def rank(got: Dict[str, Any]) -> int:
        return 1 if int((got.get('imaginary') or {}).get('count') or 0) == 1 \
            else 0

    best = tried[0]
    for got in tried[1:]:
        if got.get('xyz') and rank(got) > rank(best):
            best = got
    return best
