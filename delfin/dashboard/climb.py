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

With both put right, on the twenty-one drags as they then stood: the climb
reaches the reaction the hand pointed at **9 times in 21, in a median of
5.7 s**, against **6 in 21** before and against **9 in 21 in a median of
50.3 s** for ORCA's own OptTS started from the same geometries.  (Five of
those drags have since been rebuilt -- see below -- and on the bench as it now
is the same three routes reach 12, 12 and 10.)  "Reaches the reaction" is not
"converged": both routes also arrive, confidently, at methyl torsions at -48
cm-1 and at fragments rocking at -52, so what is counted is whether xtb's own
Hessian on what came back has one imaginary mode *and that mode is the stretch
of the pair of atoms the hand was holding*.

And then the thing that mattered more than either number: **they are wrong
about different drags**, because the climb sets off along the eigenvector
nearest the gesture and ORCA sets off along the lowest mode of its own Hessian
-- one is listening to the hand and the other to the surface.  A third search,
this same climb told to forget the gesture, is a third direction again.
:func:`reach_the_reaction` is what that is worth: the three tried cheapest
first, each one's answer checked against the contact the hand held, **16 of
21** where any one of them alone is 10 to 12.

Twenty-one, and it used to be sixteen.  Five of the drags had been built by a
constrained relaxation that tore atoms off -- an H2 floating beside a Claisen,
four loose hydrogens on a Cope -- so they were not gestures anybody could
make.  They are rebuilt the way the editor's follow really builds one: the
contact walked to where it ends a tenth of an angstrom at a time, five xtb
cycles between, everything else free.  All twenty-one are whole now by the
editor's own bond graph, and the three routes reproduce their old counts
exactly on the sixteen that never changed.  Two of the five needed the
molecule *posed* first as well, which is worth saying rather than hiding:
dragging one end of an extended 1,5-hexadiene onto the other is not a Cope but
a chain folding, and GFN2 closes two cyclopropanes instead at 3.1 A.  Posed
into the gauche conformer the Cope actually goes through, the same drag builds
it whole and all three searches reach the Cope saddle at -301 cm-1 in about
two seconds.

What is left after that is not an optimiser's failure, and the five drags none
of the three reaches say why in three different ways.  On two of them -- both
glycylglycine gestures -- the held contact is 0.00 and 0.02 of every one of
the five softest modes at the structure the hand left, so nothing that follows
a soft mode can find it from there.  One is a proton transfer whose drag
geometry already has four imaginary modes.  And the last two are the Claisen
and the Cope with their termini still 2.3 A apart, where the contact *is*
0.52 and 0.66 of a soft mode and no search converges onto it anyway -- while
the same two gestures taken 0.3 A further are reached by all three, in 1.4 and
2.0 s, at -493 and -301 cm-1.  How far the hand got is most of the answer, and
that is the honest state of a saddle search begun by hand.

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

And then what to do with the saddle once there is one, which is the last
section of this file rather than a module of its own for one reason: it is
:meth:`Climb.modes_from_xtb` that makes it possible, and that is here.  The
mode *shapes* are read out of xtb's ``hessian`` file, and a shape is what both
of the two questions about a transition state need -- what reaction is this,
which is the imaginary mode drawn (:func:`mode_pictures`), and does it join
what it was meant to join, which is that mode followed down both ways
(:func:`follow_the_mode_down`).  Measured against the rigorous form of the
second on the sixteen-atom Diels-Alder saddle: 1.0 s here against 207 s for
ORCA's own ``! XTB2 IRC``, and the same two ends.
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
from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np

from . import gfn_optimize as _gfn
from . import saddle as _saddle
from .model_hessian import model_hessian

#: Angstrom in a Bohr, from the same place the rest of the editor takes it.
BOHR = _gfn.BOHR_IN_ANGSTROM

#: A square root of a Hartree per Bohr squared per electron mass, in cm-1.
HARTREE_IN_CM = 219474.6313702

#: Kilocalories per mole in a Hartree, which is what an energy difference is
#: said in everywhere in this editor.
HARTREE_IN_KCAL = 627.5094740631

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

#: The interactive drag as overdamped steered dynamics -- see :func:`steer`.
#:
#: The hand is a spring from the dragged atom to where the cursor has it, and
#: every atom moves DOWN the total force (the field's own plus that spring) by
#: a step capped per atom.  This is a follow that carries the geometry forward
#: rather than re-minimising it, so it stays on one branch of the surface until
#: the force pushes it over -- which is the fix for the two arrangements a
#: constrained minimisation alternates between at a reaction top.
#:
#: The three numbers were measured on the prototype (a butane torsion swept
#: through a full turn with no reversal, and a butene C=C driven through its
#: twist).  ``STEER_SPRING`` is stiff enough to lead the atom and soft enough
#: that the chemistry answers; ``STEER_CAP`` is a trust radius, the most any
#: atom may move in one internal step; ``STEER_STEPS`` is how many internal
#: steps one follow answer takes, which sets how fast the structure keeps up
#: with the hand against how long an answer costs.
STEER_SPRING = 0.3          #: Hartree per Bohr squared, on the dragged atom(s)
STEER_CAP = 0.12            #: Angstrom, the most one atom moves in one step
STEER_STEPS = 4            #: gradient steps per follow answer

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

#: The displacement a differenced Hessian is taken with, in bohr.
HESSIAN_DELTA = 0.005

#: Where the *first* Hessian comes from: ``'measured'`` differences the
#: gradient for it, which is what every climb has done so far and costs ``6N``
#: gradients before the first step -- 72 for butane, 330 and thirteen minutes
#: for a manganese complex.  ``'model'`` guesses it from the geometry with
#: :func:`delfin.dashboard.model_hessian.model_hessian` and costs none.
#:
#: ``'corrected'`` is the middle course: the guess, with its own softest
#: :data:`MEASURE_THE_SOFTEST` modes measured after all, at two gradients each.
#:
#: **The default stays measured, and that is a measured decision.**  Over the
#: twenty-one drags this module was developed against -- each run twice, aimed
#: from the structure before the drag, everything else identical -- starting
#: from the guess *loses three real transition states*: a Diels-Alder at
#: -393 cm-1 replaced by one at -46 somewhere else, another at -721 lost
#: outright, and a proton transfer at -572 lost.  It gains one, at -82 cm-1
#: with the bond it is supposed to be breaking still at 1.544 A, which is a
#: conformational saddle and not the reaction.  Cheaper is not better when
#: cheaper changes the answer.
#:
#: Three attempts to make the swap safe were measured and none worked: the
#: correction above, adding the drag's own direction to what is measured
#: (0.804 against 0.804 -- it is dominated by the relaxation of everything
#: else), and scaling the model, which the surface is 4 to 6 times stiffer
#: than (1.00, 0.50 and 0.25 give the identical failure; the trust radius
#: absorbs the scale).  What goes wrong is two different things: on the
#: Diels-Alders the corrected start climbs the *same* mode -- 0.996 overlap --
#: and drifts anyway, taking 143 steps where the measured start takes 24; on a
#: proton transfer it climbs a different mode entirely, 0.002, because the
#: reaction there is a *stiff* direction and a sum of springs cannot hold one.
#:
#: Not the climb's fault, which was checked: two differenced Hessians a fifth
#: of a percent apart in step size reach the same saddle in the same 24 steps,
#: their energies equal to eight decimals.
FIRST_HESSIAN = 'measured'

#: What :data:`FIRST_HESSIAN` and the *first_hessian* argument may say.  An
#: unknown one is refused at construction rather than quietly falling back,
#: because falling back to the expensive path is exactly the mistake a typo
#: here would hide.
HESSIANS_START_FROM = ('measured', 'model', 'corrected')

#: How many of the model's own softest modes ``'corrected'`` measures rather
#: than guesses, at two gradients each.  Ten is 20 gradients against the 330 a
#: 57-atom complex costs outright.
#:
#: The number comes from where the guess is right and where it is not.
#: Measured against differenced Hessians, the model's softest mode lies inside
#: the true softest *five* to 0.96 or better on everything rigid or moderately
#: flexible, so the soft space is found there; but the *ordering* inside it is
#: not -- five model modes against five true ones overlap only 0.67 to 0.96 --
#: and it is the ordering that :meth:`Climb.aim` reads when it decides which
#: mode to climb.  Ten covers the five that matter with room for the mixing.
#:
#: On large flexible molecules the guess does not find the soft space at all
#: (0.613 for a 67-atom tetrapeptide, 0.810 for imatinib), so ten of its modes
#: are ten of the wrong ones.  That is a reason to leave the default where it
#: is rather than a number to tune.
MEASURE_THE_SOFTEST = 10

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
#: It looked like an obvious win.  Two of the nine reactions ORCA reached over
#: the bench as it then stood needed 80 and 120 cycles, so at 60 they are lost
#: -- and both of those two turned out to be among the five drags whose
#: structure had fallen apart while the bench was building it.  Over the
#: sixteen whole ones every reaction ORCA reaches converges in **13 to 49
#: cycles**, all of them under 60, and 60 and 120 both give it 7 of 16.  The
#: higher ceiling buys nothing on a structure a hand can actually make, and
#: the five have since been rebuilt so that a hand can: at 60 cycles ORCA
#: reaches **10 of the whole 21**, and the three of the five it reaches
#: converge in 4.5, 6.2 and 9.8 s.
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

    def __init__(self, numbers, bohr, method, charge, uhf, solvent, cores,
                 etemp=None):
        from xtb.interface import Calculator, Param
        from xtb.libxtb import VERBOSITY_MUTED

        self.threads = _Threads(cores)
        name = CLIMB_METHODS.get(method) or 'GFN2xTB'
        with self.threads:
            self.calc = Calculator(getattr(Param, name), numbers,
                                   np.asarray(bohr, dtype=float),
                                   charge=float(charge), uhf=int(uhf))
            self.calc.set_verbosity(VERBOSITY_MUTED)
            # Fermi smearing, the same rescue the follow loop reaches for at
            # a closed frontier gap: a closed-shell SCC fails where the two
            # frontier orbitals meet -- exactly where a bond being pulled
            # apart takes it -- and an electronic temperature holds it
            # together.  Set once here because a drag holds one gap for many
            # gradients, and it persists across ``update``.
            if etemp:
                try:
                    self.calc.set_electronic_temperature(float(etemp))
                except Exception:
                    pass
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

    def __init__(self, numbers, bohr, method, charge, uhf, solvent, cores,
                 etemp=None):
        self.binary = _gfn.find_xtb()
        self.symbols = [_SYMBOLS[int(z)] for z in numbers]
        self.method = method
        self.charge = int(charge)
        self.uhf = int(uhf)
        self.solvent = solvent
        self.cores = max(1, int(cores))
        self.etemp = float(etemp) if etemp else None
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
        if self.etemp:
            order += ['--etemp', str(self.etemp)]
        if self.solvent:
            order += ['--alpb', str(self.solvent)]
        room = dict(os.environ)
        room['OMP_NUM_THREADS'] = str(self.cores)
        room['MKL_NUM_THREADS'] = str(self.cores)
        # No clock.  A gradient is one xtb call and normally a fraction of a
        # second, but under g-xTB it is a process with an ORCA start-up in
        # front of it, and on a large system two minutes is a limit that stops
        # the climb rather than a hang that needed stopping.  The climb itself
        # is asked to stop between steps, which is the way out.
        subprocess.run(order, cwd=str(self.folder), env=room,
                       capture_output=True, text=True)
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


def _make_gradients(numbers, bohr, method, charge, uhf, solvent, cores,
                    etemp=None):
    """A callable ``(bohr) -> (energy, gradient)`` in atomic units.

    In process where the module is there and the method is not GFN-FF -- the
    same choice :class:`Climb` makes, and for the same reason: an in-process
    GFN-FF drops a topology file into the launch directory.
    """
    fast = have_fast_gradients() and str(method) != 'gfnff'
    return (_InProcess if fast else _CommandLine)(
        numbers, bohr, method, charge, uhf, solvent, cores, etemp=etemp)


def steer(start_xyz, wish_xyz, held, *, method='gfn2', charge=0, uhf=0,
          solvent=None, cores=1, etemp=None,
          k_spring=STEER_SPRING, cap_angstrom=STEER_CAP, steps=STEER_STEPS,
          should_stop=None):
    """One drag answer by overdamped steered dynamics -- the continuous follow.

    The current follow re-minimises every frame: it pins the dragged coordinate
    at the cursor and relaxes the rest.  At a reaction top that biased minimum
    has two solutions, and the minimiser flips between them answer to answer --
    the two arrangements a molecule cannot be held on.  This does not minimise.
    It puts a spring from each dragged atom to where the cursor has it, adds
    that to the field's own force, and moves every atom DOWN the total by a step
    no longer than ``cap_angstrom`` -- the quasi-static limit of interactive
    molecular dynamics.  Because it carries the geometry forward from the last
    answer, it stays on one branch until the force pushes it over, which is
    what a hand on a real molecule would feel and what lets every process --
    a rotation, a reaction, an isomerisation -- be walked rather than jumped.

    ``start_xyz`` is the geometry to advance (the last answer, not the cursor's
    wish); ``wish_xyz`` supplies only where the cursor has put the ``held``
    atoms, which is the spring's target.  The two carry the same atoms in the
    same order.  ``etemp`` smears the gradient the way the follow loop smears a
    relaxation at a closed gap.

    Returns a dict shaped like :func:`gfn_optimize.relax_steps`' -- ``ok``,
    ``xyz``, ``energy`` (of the geometry returned, Hartree), ``calls`` -- or
    ``ok`` ``False`` with a ``status`` where a gradient would not run.
    """
    start = _elements(start_xyz)
    wish = _elements(wish_xyz)
    if start is None or wish is None:
        return {'ok': False, 'xyz': start_xyz,
                'status': 'the structure could not be read'}
    numbers = start['numbers']
    symbols = start['symbols']
    if len(wish['angstrom']) != len(start['angstrom']):
        # A frame from before the structure changed under the drag -- there is
        # no target to spring towards, so nothing is driven.
        return {'ok': False, 'xyz': start_xyz,
                'status': 'the hand is on a structure that is no longer here'}
    held = [int(i) for i in (held or ())
            if 0 <= int(i) < len(numbers)]
    rb = np.asarray(start['angstrom'], dtype=float) / BOHR
    target = np.asarray(wish['angstrom'], dtype=float) / BOHR
    cap = float(cap_angstrom) / BOHR
    try:
        grad = _make_gradients(numbers, rb, method, charge, uhf, solvent,
                               cores, etemp=etemp)
    except Exception as oops:                        # noqa: BLE001
        return {'ok': False, 'xyz': start_xyz, 'status': str(oops)}
    energy = None
    try:
        for _ in range(max(1, int(steps))):
            if should_stop is not None and should_stop():
                break
            energy, gradient = grad(rb)
            force = -np.asarray(gradient, dtype=float)
            for i in held:
                # The spring only on the dragged atom(s): a force with a
                # ceiling that leads the atom towards the cursor and is weaker
                # the closer it gets, so the chemistry keeps a say.
                force[i] += k_spring * (target[i] - rb[i])
            norm = np.linalg.norm(force, axis=1, keepdims=True)
            # A trust radius per atom, not on the whole step: the dragged atom
            # is where the force is largest and the rest must not be dragged
            # bodily behind it.
            scale = np.minimum(1.0, cap / (norm + 1e-12))
            rb = rb + force * scale
        else:
            # The energy the loop leaves is the geometry before its last move,
            # so the returned structure is priced by one more single point --
            # skipped only when the hand let go mid-answer.
            if should_stop is None or not should_stop():
                energy, _ = grad(rb)
    except Exception as oops:                        # noqa: BLE001
        return {'ok': False, 'xyz': start_xyz, 'status': str(oops),
                'calls': int(getattr(grad, 'calls', 0))}
    finally:
        closer = getattr(grad, 'close', None)
        if closer is not None:
            closer()
    return {'ok': True,
            'xyz': xyz_document(symbols, rb * BOHR,
                                'Steered by the hand on the structure'),
            'energy': energy, 'calls': int(getattr(grad, 'calls', 0))}


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
                 cores: int = 4, trust: float = START_TRUST,
                 first_hessian: Optional[str] = None):
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
        self.first_hessian = str(first_hessian or FIRST_HESSIAN).lower()
        if self.first_hessian not in HESSIANS_START_FROM:
            raise ValueError(
                'a climb starts from '
                + ' or '.join(HESSIANS_START_FROM)
                + f', not from {first_hessian!r}')
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

    def numerical_hessian(self, delta: float = HESSIAN_DELTA) -> np.ndarray:
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

    def first_matrix(self, aimed_from: Any = None) -> np.ndarray:
        """The Hessian the climb sets off with, measured or guessed.

        Which of the two is :attr:`first_hessian`, and the difference is the
        whole of what it costs to begin: ``6N`` gradients against none.  The
        refreshes along the way are not affected -- those are always measured,
        because by then the climb is somewhere the model has no claim to know
        about and a guess would be replacing something better with something
        worse.
        """
        if self.first_hessian == 'measured':
            return self.numerical_hessian()
        guess = (model_hessian(self.numbers, self.bohr)
                 / np.outer(self.scale, self.scale))
        if self.first_hessian == 'model':
            return guess
        return self._corrected(guess, self._asked_for(aimed_from))

    def _asked_for(self, aimed_from: Any) -> Optional[np.ndarray]:
        """The way the hand moved the structure, as a working-metric unit."""
        if aimed_from is None:
            return None
        before = _elements(aimed_from) if isinstance(aimed_from, str) \
            else {'angstrom': np.asarray(aimed_from, dtype=float)}
        if before is None:
            return None
        earlier = np.asarray(before['angstrom'], dtype=float)
        if earlier.shape != self.bohr.shape:
            return None
        moved = ((self.bohr - earlier / BOHR).reshape(-1)) * self.scale
        size = float(np.linalg.norm(moved))
        return None if size < 1e-8 else moved / size

    def _product(self, way: np.ndarray) -> np.ndarray:
        """The real Hessian applied to one direction, for two gradients.

        A whole Hessian is 6N of these.  One is two, whatever the molecule's
        size, which is the only reason any of this is worth doing.
        """
        step = np.asarray(way, dtype=float) / self.scale
        reach = HESSIAN_DELTA / max(1e-12, float(np.linalg.norm(step)))
        flat = self.bohr.reshape(-1)
        _e, up = self._measure((flat + reach * step).reshape(-1, 3))
        _e, down = self._measure((flat - reach * step).reshape(-1, 3))
        return ((up - down) / (2.0 * reach)) / self.scale

    def _corrected(self, guess: np.ndarray,
                   asked: Optional[np.ndarray] = None) -> np.ndarray:
        """The model, with its own softest modes measured instead of guessed.

        The guess finds the soft *space* and gets the order inside it wrong,
        so the order is the only thing bought back: the model's softest
        :data:`MEASURE_THE_SOFTEST` directions are handed to the real surface,
        one gradient pair each, and that small block is replaced by what comes
        back.  Everything stiffer stays as the model had it, where being four
        to six times too stiff costs a climb nothing -- it is not going that
        way.

        *asked* is the way the hand moved the structure, and it is measured
        too, because the reaction is not always a soft mode.  Measured on a
        proton half transferred across a hydrogen bond: the true mode is a
        stretch at -0.387 in the working metric, nowhere near the soft end, and
        correcting only the softest ten brought the matrix to -0.072 -- the
        right sign and a fifth of the size.  A drag says where to look, and
        looking there costs one more gradient pair.

        The correction can also do the one thing the model never can.  A sum of
        springs cannot hold a negative curvature; a measured block can, and on
        a Diels-Alder held at 1.10 A it does -- the guess said +0.0037 where the
        surface says -0.0080, and the correction says -0.0071.  A climb that
        starts with nothing to climb is the failure this addresses.
        """
        basis = self._basis()
        small = basis.T @ guess @ basis
        _value, vector = np.linalg.eigh(0.5 * (small + small.T))
        want = int(min(MEASURE_THE_SOFTEST, vector.shape[1]))
        if want < 1:
            return guess
        ways = basis @ vector[:, :want]
        if asked is not None:
            # Orthogonal to the soft ones, or it is asking a second time for
            # something already bought.
            spare = asked - ways @ (ways.T @ asked)
            size = float(np.linalg.norm(spare))
            if size > 1e-3:
                ways = np.concatenate([ways, (spare / size)[:, None]], axis=1)
                want += 1
        image = np.stack([self._product(ways[:, i]) for i in range(want)],
                         axis=1)
        was = ways.T @ guess @ ways
        now = ways.T @ image
        now = 0.5 * (now + now.T)
        return guess + ways @ (now - was) @ ways.T

    def start(self, aimed_from: Any = None) -> Dict[str, Any]:
        """Pay for the Hessian, and choose the mode to climb.

        *aimed_from* is the structure as it stood before the user moved it.
        Given one, the displacement between then and now is the direction the
        hand asked for, and the mode climbed is whichever of the Hessian's
        eigenvectors most resembles it -- which is the whole of how a drag
        guides a saddle search.  See :meth:`took`.
        """
        began = time.perf_counter()
        self.hessian = self.first_matrix(aimed_from)
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
                'from': self.first_hessian,
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
        # And whether the Hessian was taken anywhere a Hessian means
        # anything.  Curvature is a statement about a stationary point; at a
        # point on a slope the same numbers are not the frequencies of
        # anything, and the module next door says so in as many words --
        # :func:`gfn_optimize.not_a_stationary_point`.  Measured on a
        # cyclohexene whose C0-C1 was driven to 1.90 A and then climbed: the
        # press reported "it converged onto a different saddle, at -205 cm-1"
        # while the gradient at that geometry was 3.3e-03 Hartree per Bohr
        # per coordinate, and the editor's own slope test, given the same
        # structure, answered "this structure is neither a minimum nor a
        # saddle".  Costs nothing: the gradient is the one the walk already
        # has.  ``None`` where there is none to read -- a verdict asked of a
        # structure this walk has not stepped on -- because unknown is not
        # the same as moving.
        still = gmax = grms = None
        if self.gradient is not None:
            gmax = float(np.abs(self.gradient).max())
            grms = float(np.sqrt(np.mean(
                (self._basis().T @ self.gradient) ** 2)))
            still = bool(gmax < GRADIENT_MAX and grms < GRADIENT_RMS)
        return {'count': len(wrong), 'modes': wrong[:4],
                'lowest': float(modes[0]) if modes.size else None,
                'ok': first, 'share': share, 'still': still,
                'gmax': gmax, 'grms': grms,
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
                # No clock.  This is the check that says whether what was
                # climbed to is a transition state at all, and cut off it says
                # nothing -- which reads as "no verdict" rather than as "the
                # verdict did not fit in five minutes".  A Hessian on fifty
                # atoms is already 11.8 s; on the sizes this editor is
                # reaching for now it is the answer worth waiting for.
                cwd=str(folder), env=room, capture_output=True, text=True)
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
        arrived = halted = False
        for _ in range(max(1, int(max_steps))):
            if should_stop is not None and should_stop():
                halted = True
                break
            outcome = walk.step()
            if on_frame is not None:
                on_frame(walk.frame(), outcome)
            if outcome.get('converged'):
                arrived = True
                break
        # A press the user ended is not a question about where the walk was
        # standing, and the verdict is the expensive part of this function:
        # xtb's own Hessian, 0.3 s at sixteen atoms and 11.8 at fifty, paid at
        # the exact moment somebody asked for it to stop.  What they are left
        # with is the frame on screen, which the caller knows about and this
        # does not.
        shape = None if halted else walk.verdict(held=held)
        seconds = time.perf_counter() - began
        # ``ok`` is the ladder's to rewrite -- a climb that converged onto
        # the wrong saddle has not reached the reaction -- so whether the
        # gradient actually vanished is kept under its own name as well.
        return {'ok': arrived, 'converged': arrived, 'xyz': walk.xyz(
                    'Climbed to a transition state' if arrived
                    else 'Where the climb got to'),
                'seconds': seconds, 'steps': walk.steps, 'stopped': halted,
                'started': opened, 'imaginary': shape,
                'gradients': int(getattr(walk.engine, 'calls', 0)),
                'status': (
                    f'The climb converged in {walk.steps} steps '
                    f'({seconds:.1f} s).' if arrived else
                    f'The climb was stopped after {walk.steps} steps.'
                    if halted else
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
    if shape.get('still') is False:
        return _on_a_slope(shape, short=True)
    count = int(shape.get('count') or 0)
    modes = shape.get('modes') or []
    if count == 0:
        return 'it converged onto a minimum'
    if count == 1:
        return ('it converged onto a different saddle, at '
                f'{modes[0]:.0f} cm-1')
    return (f'it converged onto a saddle of order {count} ('
            + ', '.join(f'{one:.0f}' for one in modes[:3]) + ' cm-1)')


def _on_a_slope(shape: Dict[str, Any], short: bool = False) -> str:
    """A Hessian taken where nothing was standing still, said plainly.

    The modes are still computed and still reported -- they are what the walk
    was following, and hiding them would leave the user with nothing to read.
    What is not said is that they are the modes of a saddle or of a minimum,
    because at a point on a slope they are the modes of neither.

    *short* is for the second place in a sentence that has already said this
    once: the ladder prefixes every later line with why the first rung missed,
    so the long form twice over is one line saying one thing twice.
    """
    largest = shape.get('gmax')
    count = int(shape.get('count') or 0)
    modes = shape.get('modes') or []
    deepest = f'{modes[0]:.0f} cm-1' if count else ''
    if short:
        said = 'it too ran out of tries before anything stood still'
        if deepest:
            said += f', at {deepest}'
        if largest is not None:
            said += f' with the gradient still {largest:.1e}'
        return said
    said = ('it stopped where it had run out of tries rather than where the '
            'structure had stopped moving')
    if largest is not None:
        said += (f' -- the gradient there is still {largest:.1e} Hartree per '
                 f'Bohr against the {GRADIENT_MAX:.0e} it converges on')
    if count == 1:
        said += f', so the {deepest} there is not a saddle\'s'
    elif count:
        said += (', so the ' + ', '.join(f'{one:.0f}' for one in modes[:3])
                 + " cm-1 there are not a saddle's")
    else:
        said += ', so nothing there says it is a minimum'
    return said


def _why_not(shape: Dict[str, Any], contact: str,
             advise: bool = True) -> str:
    """Why the climb's answer was not taken, said in one sentence.

    *advise* is off when the ladder has rungs left to try: this sentence is
    prefixed to every later one, and advice inside it would be repeated
    against whatever the last rung ends up saying.
    """
    if shape.get('still') is False:
        return ('The climb did not converge: ' + _on_a_slope(shape) + '.'
                + (' Pressing again starts a fresh climb from where it got '
                   'to, with its tries counted from nought.' if advise
                   else ''))
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
    somewhere, and usually somewhere else.  Measured over twenty-one hand
    drags on ten to fifty atoms, three searches that all start from the
    structure the hand left each reach the reaction the hand pointed at about
    half the time, **and not on the same drags**:

    ==================================  ========  ===============
    where it starts climbing              of 21    told by
    ==================================  ========  ===============
    the eigenvector nearest the drag        12     :meth:`aim`
    the surface's own softest mode          12     no *aimed_from*
    ORCA's OptTS, redundant internals       10     :mod:`saddle`
    **all three, in that order**          **16**
    ==================================  ========  ===============

    One column now, and it used to be two.  Twenty-one drags were built by
    holding one contact and relaxing the rest under GFN2, and on five of them
    that relaxation tore atoms off -- an H2 floating beside a Claisen, four
    loose hydrogens on a Cope -- so for a while the sixteen whole ones were
    the honest denominator and the counts were 9, 9, 7 and 13 of 16.  The five
    are rebuilt: the contact walked to where it ends a tenth of an angstrom at
    a time with five xtb cycles between, which is what the editor's follow
    does, and the bond graph checked after every answer.  All twenty-one are
    whole now, the three routes reproduce 9, 9 and 7 exactly on the sixteen
    that did not change, and of the five rebuilt drags three are reached and
    two are not.  Two of the five also had to be *posed* before they were
    dragged, which the module docstring says out loud rather than hides.

    Which is the useful half of rebuilding them.  A Claisen whose termini are
    2.0 A apart and a Cope at 2.0 A are reached by **all three** searches, in
    1.4 and 2.0 s, at -493 and -301 cm-1 -- so the ladder answers on its first
    rung where it used to have nothing to answer with.  The same two gestures
    stopped 0.3 A earlier are reached by none of the three.  How far the hand
    got is most of the answer.

    The shape is unchanged and it is the whole reason this function exists:
    twelve, twelve and ten, and sixteen between them.  Closing that gap needs
    no new optimiser -- it needs the three that are already here, tried in
    order, with what each one reached actually checked.

    *Why* they differ is not a property of the molecules.  It is which
    direction each one sets off along, and each wins where *its* direction is
    already the reaction:

    * the aimed climb reached it on **11 of the 15** drags whose gesture is
      mostly one soft eigenvector (overlap 0.60 or better) and on 1 of the 6
      below that;
    * ORCA reached it on **10 of the 11** drags where the softest mode at the
      structure the hand left is already the held contact stretching (0.60 or
      better), and on **none of the 10** where it is not.

    So on the Diels-Alder drags only the aimed climb gets, the softest mode is
    two fragments rocking -- 0.01 to 0.11 of the reaction -- and the gesture
    is 0.51 to 0.82 of it.  On the retro-Diels-Alder and the salicylaldehyde
    proton transfer only ORCA gets, the softest mode is 0.80 to 0.94 of the
    reaction.  One of them is listening to the hand and the other to the
    surface, and where those disagree is exactly where one beats the other.

    The middle rung is what that observation buys.  The same climb, told to
    forget the gesture and follow the softest mode, is ORCA's starting
    direction at the climb's price: it catches three drags the aimed climb
    misses -- the two retro-Diels-Alders and a Diels-Alder ORCA never reaches
    at all -- in 3.9, 2.4 and 2.9 s where ORCA needs 12.5, 8.7 and 22.6.  It
    is the same class, the same code and one more Hessian, and it is what
    takes the whole thing from 15 of the twenty-one to 16.

    And the five that are never reached are not an optimiser's failure.  On
    two of them -- both glycylglycine gestures -- the held contact is 0.00 and
    0.02 of every one of the five softest modes at the structure the hand
    left, so nothing that follows a soft mode from there can find it.  One is
    a proton transfer whose drag geometry already has four imaginary modes.
    And two are the Claisen and the Cope at 2.3 A, where the contact is 0.52
    and 0.66 of a soft mode and no search converges onto it anyway, while the
    same gestures 0.3 A further in are reached by all three.  Those are
    guesses that missed, and the sentence says so rather than implying that
    another minute of computing would help.

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
    That looks like throwing work away and it is measured to be right.  Over
    the twelve drags the aimed climb got wrong on the bench as it then stood,
    ORCA started from the hand's guess reaches the reaction on 5 and ORCA
    started from where the climb stopped reaches it on **2** -- handing the
    endpoint on loses the
    retro-Diels-Alder at 2.4 A, the Claisen at 2.0 and the salicylaldehyde
    proton transfer, all three of which the same optimiser reaches from the
    hand's structure.  A search that is going wrong does not stop somewhere
    useful.  So a rung that misses costs its wall time and nothing else, and
    :data:`CLIMB_CEILING` keeps that short.

    **This order and not the other.**  Measured on the bench as it then
    stood: ORCA first, with the climb after it, reaches 9 of 21 rather than
    14, because started from where ORCA stopped the climb rescues **0 of the
    12** ORCA gets wrong -- ORCA leaves it standing on a stationary point with
    no gradient to work with.  Started from the hand's structure instead it
    would reach the same union, but it would pay ORCA's median 44 s on every
    press, including the ones the climb answers in under two.  Cheapest first
    is both the same answer and sooner.

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
    if first.get('stopped'):
        # Ended inside the first rung.  There is nothing to check -- the
        # verdict was not asked for, because it costs a Hessian at the moment
        # somebody pressed Stop -- and nothing to try, because a rung is only
        # tried when the one before it is known to have missed.
        first['ok'] = False
        return first
    reached = first.get('imaginary') or {}
    # Without a pair to check against, only the mode count can be asked, and
    # then the ladder is no use: every rung would be judged by the very test
    # this one has already passed, so the rest could only spend a minute to
    # agree.  The weaker test is used and the sentence says which one it was.
    told_which = reached.get('share') is not None
    # A Hessian at a point on a slope is not a verdict about a saddle, so a
    # rung that ran out of tries has not reached anything however its modes
    # come out -- see :meth:`Climb.verdict`, where ``still`` is measured.  The
    # ladder carries on to the next rung, which is what a rung that missed has
    # always done.
    landed = reached.get('still') is not False
    if landed and bool(reached.get('reaction') if told_which
                       else reached.get('ok')):
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
    more = bool(told_which and fallback)
    why = _why_not(reached, contact, advise=not more)
    first['ok'] = False
    first['status'] = why
    if not more:
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
        if again.get('stopped'):
            again['ok'] = False
            again['status'] = (
                why + ' Climbing along the softest mode instead was stopped '
                'before it reached anything.')
            return again
        softest = again.get('imaginary') or {}
        if softest.get('still') is not False and softest.get('reaction'):
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
    orca = last.get('imaginary') or {}
    if orca.get('still') is not False and orca.get('reaction'):
        return _arrived(last, "ORCA's optimiser, given the structure you made,",
                        why + ' ')

    # Nothing reached it.  What comes back is whichever attempt got closest to
    # being a transition state at all, and the first keeps a tie because it is
    # the one the user watched.
    best = _closest(tried)
    best['ok'] = False
    # Whose fault it was.  "The structure has to change" is the right advice
    # for a search that ran to a stationary point somewhere else, and the
    # wrong one for a search that simply ran out of tries -- there the
    # structure may be perfectly good and the ceiling is what stopped it.
    advice = (
        ' Every one of them ran out of tries rather than arriving, so press '
        'again from what is on screen: each press starts a fresh climb from '
        'there with its tries counted from nought.'
        if (best.get('imaginary') or {}).get('still') is False else
        (' All three start from where you let go, so it is the structure '
         'that has to change -- move the atoms nearer the arrangement the '
         'reaction passes through and press again.'
         if len(tried) > 2 else
         ' Both start from where you let go, so it is the structure that has '
         'to change -- move the atoms nearer the arrangement the reaction '
         'passes through and press again.'))
    best['status'] = (
        why + f' Climbing along the surface\'s own softest mode and then '
        f'ORCA\'s optimiser were given the same structure and did not reach '
        f'it either ({spent:.1f} s altogether). What is left on screen is the '
        f'nearest any of the three came: '
        + _named_landing(best.get('imaginary')) + '.'
        if len(tried) > 2 else
        why + f' ORCA\'s optimiser was given the same structure and did not '
        f'reach it either ({spent:.1f} s altogether): '
        + _named_landing(best.get('imaginary')) + '.') + advice
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


# -- what to do with the saddle once you have one ----------------------
#
# A converged first-order saddle is a geometry and a number, and neither of
# those answers the two questions a chemist has about it.  What reaction is
# this -- the imaginary mode *is* the reaction coordinate, so drawing it says
# which bonds are forming and which are breaking, directly, without anybody
# having to read coordinates off a box.  And does it join what it was meant to
# join -- a different question, and the standard way to ask it is to push the
# structure a little way down that mode in each direction and see where it
# settles.
#
# Both start from :meth:`Climb.modes_from_xtb`, which is here already: xtb's
# own Hessian read out of the ``hessian`` file rather than off the printout,
# because the file carries the *shapes* of the modes and the printout carries
# only their numbers.  A shape is exactly what both of these need.  One xtb
# either way, measured at 0.41 s on sixteen atoms.


#: How far along a mode the picture is drawn, in Angstrom of the largest
#: single-atom displacement.
#:
#: It sits between two bounds, both measured on the same structure.  It has to
#: be big enough that the chemistry is unmistakable: on the sixteen-atom
#: Diels-Alder saddle, whose two forming bonds stand at 2.315 A, this
#: amplitude swings both of them between 1.63 A -- a formed C-C bond -- and
#: 3.01 A, which is no bond at all.  Nobody has to be told what that picture
#: is of.
#:
#: And it has to be small enough not to tear the molecule, which is the
#: question :data:`gfn_optimize._TOO_CLOSE` already answers for a hand: two
#: atoms inside two thirds of the bond they would make are not on any path at
#: any temperature.  On that same saddle the tightest contact anywhere in the
#: swing is 0.81 of a bond at this amplitude, and first drops under the 0.65
#: that counts as squeezed at an amplitude of 0.70.  So the default is half of
#: where this molecule first tears.
#:
#: Which is a fact about that molecule rather than about amplitudes, so it is
#: a starting point and not a rule: :func:`amplitude_that_fits` cuts it down
#: on a structure whose mode drives two atoms together sooner.  At the other
#: extreme, measured on planar ammonia at its inversion saddle -- where the
#: mode is three hydrogens swinging through the plane -- nothing is tight at
#: any amplitude up to 1.0 A, and this is nowhere near the limit.
MODE_AMPLITUDE = 0.35

#: Below this an amplitude is no longer a picture of anything, so the cutting
#: down stops rather than converging on nought.  A twentieth of an Angstrom is
#: under a twentieth of a bond, which no screen shows as motion.
LEAST_AMPLITUDE = 0.05

#: How many pictures make one swing: out one way, back through the structure,
#: out the other and back.
#:
#: Twenty-four of them at :data:`MODE_PACE_MS` is about a swing a second --
#: slow enough that the eye can follow one bond through it, fast enough that
#: it reads as motion rather than as a slide show.  The player interpolates
#: between pictures at the screen's own rate, so this is how many geometries
#: are computed and not how many are drawn.
MODE_PICTURES = 24

#: And how many swings one press draws before the picture is back where it
#: started.
#:
#: Six is about six seconds, which is long enough to watch the same motion
#: three or four times over and decide what it is.  Bounded on purpose: a
#: picture that runs until it is switched off is a second thing to remember to
#: switch off, and this one cannot be left running by accident.  Pressing
#: again is one press.
MODE_SWINGS = 6

#: How long one picture stands, in milliseconds.
#:
#: Named here rather than taken from the editor's playback slider, and the
#: difference between the two is the whole reason: the slider says how fast to
#: walk a path that was computed step by step, and at its top it means "do not
#: fall behind the calculation" -- which, for a path that is already finished,
#: is the whole animation inside one screen frame.  A mode is not a path being
#: walked.  It has a period, and the period is what makes it legible.
MODE_PACE_MS = 40

#: How far the structure is pushed down the mode before it is let go, in
#: Angstrom of the largest single-atom displacement.
#:
#: A saddle is a stationary point, so this only has to be enough to leave one;
#: the relaxation does the rest.  Measured on the Diels-Alder saddle at 0.1,
#: 0.2, 0.3 and 0.5 A, the two ends come back the same to two decimals every
#: time -- -70.64 kcal/mol with the ring closed at 1.54 A one way, -6.75 with
#: the two fragments 3.3 A apart the other -- so the answer does not depend on
#: it anywhere over that range.  ORCA's own IRC, asked for its default 2 mEh
#: drop on the same structure, chose 0.287 along its normalised mode, which is
#: the same size written in a different unit.
#:
#: 0.3 rather than 0.1 because a displacement that barely leaves a stationary
#: point leaves the optimiser to find its own way off one, and on a flat
#: saddle that is where it crawls or comes back.
DESCENT_STEP = 0.30


def modes_of(xyz_text: str, method: str = 'gfn2', *,
             charge: int = 0, uhf: int = 0, solvent: Optional[str] = None,
             cores: int = 4) -> Optional[Dict[str, Any]]:
    """The modes at this geometry, with shapes something can be drawn from.

    ``{'cm', 'ways', 'symbols', 'angstrom'}``.  *cm* are the frequencies in
    cm-1, sorted from softest; *ways* has one column per mode, in **Cartesian**
    Angstrom and scaled so that the atom moving furthest in it moves exactly
    one -- so an amplitude in Angstrom multiplies straight into a column and
    means what it says.  ``None`` when there is no xtb to ask.

    :meth:`Climb.modes_from_xtb` hands its shapes over mass-weighted, which is
    the metric the eigenvalues live in and not the one a picture is drawn in.
    Left that way a hydrogen moves a third of what it really does, and a proton
    transfer looks like the heavy atoms doing the work.  Dividing by the root
    mass puts it back.
    """
    walk = Climb(xyz_text, method, charge=charge, uhf=uhf, solvent=solvent,
                 cores=cores)
    try:
        got = walk.modes_from_xtb()
        if got is None:
            return None
        ways = np.asarray(got['shape'], dtype=float) / walk.weight[:, None]
        for column in range(ways.shape[1]):
            furthest = float(np.linalg.norm(
                ways[:, column].reshape(-1, 3), axis=1).max())
            if furthest > 1e-12:
                ways[:, column] /= furthest
        return {'cm': np.asarray(got['cm'], dtype=float), 'ways': ways,
                'symbols': list(walk.symbols),
                'angstrom': np.asarray(walk.angstrom, dtype=float)}
    finally:
        walk.close()


def imaginary_among(cm: Any) -> List[int]:
    """Which of these modes go the wrong way, softest first.

    Counted against the same :data:`IMAGINARY_BELOW` as everything else here,
    so a mode a verdict called imaginary is one that can be drawn and followed
    and no other is.
    """
    return [n for n, one in enumerate(np.asarray(cm, dtype=float))
            if float(one) < IMAGINARY_BELOW]


def displaced_along(angstrom: Any, way: Any, amplitude: float) -> np.ndarray:
    """This geometry moved *amplitude* Angstrom along *way*."""
    rows = np.asarray(angstrom, dtype=float)
    return rows + float(amplitude) * np.asarray(
        way, dtype=float).reshape(rows.shape)


def amplitude_that_fits(symbols: Sequence[str], angstrom: Any, way: Any,
                        wanted: float = MODE_AMPLITUDE) -> float:
    """The largest amplitude at or below *wanted* that tears nothing.

    Asked at both ends of the swing, because a mode is symmetric in the
    arithmetic and not in the molecule: one direction closes a contact and the
    other opens it, and only one of the two can be the one that squeezes.

    "Tears" is not a new idea here.  It is
    :func:`gfn_optimize.closest_contact` against
    :data:`gfn_optimize._TOO_CLOSE`, which is the floor a drag is already held
    to.  A bond graph would be the wrong test, and it is worth saying why: the
    whole point of this picture is that bonds appear and disappear in it, so a
    rule that refused a changed graph would refuse every mode worth looking
    at.  Measured on the Diels-Alder saddle, the graph changes at an amplitude
    of 0.2 A and nothing is tight until 0.70.

    Returns *wanted* unchanged where nothing is tight, which is the ordinary
    case; something smaller where it is; and never below
    :data:`LEAST_AMPLITUDE`, under which there is no picture left to protect.
    """
    rows = np.asarray(angstrom, dtype=float)
    names = list(symbols)
    amplitude = float(wanted)
    while amplitude > LEAST_AMPLITUDE:
        worst = None
        for sign in (1.0, -1.0):
            text = xyz_document(
                names, displaced_along(rows, way, sign * amplitude),
                'along a mode')
            tight = _gfn.closest_contact(text)[0]
            if tight is not None and (worst is None or tight < worst):
                worst = tight
        if worst is None or worst >= _gfn._TOO_CLOSE:
            return amplitude
        amplitude *= 0.8
    return LEAST_AMPLITUDE


def mode_pictures(angstrom: Any, way: Any, *, amplitude: float,
                  pictures: int = MODE_PICTURES,
                  swings: int = MODE_SWINGS) -> List[List[float]]:
    """The frames that show one mode, flat, the way the frame channel wants.

    ``x0 + A sin(phase)``: it begins at the structure, goes out one way, back
    through it, out the other, and **ends on the structure exactly**, because
    the last picture is a whole number of swings and the sine of that is
    nought.

    That last part is not a detail.  A geometry displaced along a mode is a
    picture and not a structure anybody chose, so the one thing this must
    never do is leave the viewer standing on one: whatever cuts the animation
    short, the frame it comes to rest on has to be the geometry the coordinate
    box holds.  Ending on it is how that is true when the animation runs out
    by itself, and the run's ``home`` frame is how it is true when a hand
    arrives in the middle of one.
    """
    rows = np.asarray(angstrom, dtype=float)
    step = np.asarray(way, dtype=float).reshape(rows.shape)
    out: List[List[float]] = []
    many = max(2, int(pictures))
    for n in range(many * max(1, int(swings)) + 1):
        far = float(amplitude) * math.sin(2.0 * math.pi * n / many)
        out.append([round(float(one), 4)
                    for one in (rows + far * step).reshape(-1)])
    return out


#: Two ends count as the same energy when they are within this, in kcal/mol.
#:
#: A tenth of a kcal/mol is under what a semiempirical method means by a
#: number, and well under what separates two conformers of anything.  It is
#: here to answer one question -- did both directions come back to the same
#: place -- and not to grade anything.
SAME_ENERGY = 0.1

#: And the same geometry when they are within this, in Angstrom of RMSD after
#: the best rotation.
#:
#: A twentieth of an Angstrom is a hundredth of a bond: two relaxations of the
#: same minimum from different starting points land inside it, and two
#: genuinely different arrangements never do.
SAME_PLACE = 0.05


def pieces_in(xyz_text: str) -> int:
    """How many separate structures this geometry is.

    By :func:`gfn_optimize.bond_graph`, so it is the same bonding the viewer
    draws lines with and the same one the topology wall compares.  It is the
    one thing about a relaxed structure that can be said without knowing any
    chemistry at all: two fragments 3 A apart are two pieces whatever they are
    made of.
    """
    rows = _gfn.atom_lines(xyz_text)
    if not rows:
        return 0
    home = list(range(len(rows)))

    def root(one: int) -> int:
        while home[one] != one:
            home[one] = home[home[one]]
            one = home[one]
        return one

    for left, right in _gfn.bond_graph(xyz_text):
        a, b = root(int(left)), root(int(right))
        if a != b:
            home[a] = b
    return len({root(one) for one in range(len(rows))})


def turned_onto(here: Any, there: Any) -> float:
    """RMSD between two geometries once one is turned onto the other.

    Kabsch, with the reflection taken back out: two structures that are mirror
    images are *not* the same structure, and an alignment that is allowed to
    invert would report them as one.  That case is not exotic -- an umbrella
    inversion has it at both ends of its own mode.
    """
    one = np.asarray(here, dtype=float)
    two = np.asarray(there, dtype=float)
    if one.shape != two.shape or one.size == 0:
        return float('nan')
    one = one - one.mean(0)
    two = two - two.mean(0)
    left, _, right = np.linalg.svd(one.T @ two)
    mirror = np.diag([1.0, 1.0, float(np.sign(np.linalg.det(left @ right)))])
    turn = left @ mirror @ right
    return float(np.sqrt(((one @ turn - two) ** 2).sum() / len(one)))


def _pair_named(symbols: Sequence[str], pair: Sequence[int],
                angstrom: Any) -> str:
    """One bond, named the way the viewer numbers its atoms and measured.

    Zero-based, for the reason :func:`_contact` gives: the ``#`` button draws
    the index on the atom counting from nought, and a sentence that numbers
    them differently from the picture is worse than one that does not number
    them at all.
    """
    one, two = int(pair[0]), int(pair[1])
    rows = np.asarray(angstrom, dtype=float)
    apart = float(np.linalg.norm(rows[one] - rows[two]))
    return f'{symbols[one]}{one}-{symbols[two]}{two} at {apart:.2f} A'


def _bonds_said(symbols: Sequence[str], changed: Sequence[Sequence[int]],
                angstrom: Any, most: int = 3) -> str:
    """A list of bonds, as a phrase, with the long tail counted rather than
    printed."""
    listed = [_pair_named(symbols, pair, angstrom)
              for pair in list(changed)[:most]]
    rest = len(list(changed)) - len(listed)
    return ', '.join(listed) + (f' and {rest} more' if rest > 0 else '')


#: Small counts written out, because a sentence with a digit in the middle of
#: it reads as a measurement and these are not measurements.
_SPELT = ('no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven',
          'eight', 'nine', 'ten')


def _many(count: int) -> str:
    """A small number in words, and a large one as itself."""
    return _SPELT[count] if 0 <= int(count) < len(_SPELT) else str(int(count))


def _pieces_said(count: int) -> str:
    """How many separate structures this is, as a phrase."""
    return 'one piece' if int(count) == 1 else f'{_many(count)} separate pieces'


def _end_said(end: Dict[str, Any], symbols: Sequence[str], which: str) -> str:
    """What one end of the descent is, in a sentence that assumes no chemistry.

    Everything in it is a measurement: how many separate pieces the structure
    is, which bonds it has that the saddle did not and which the saddle had
    that it does not, and what it costs against the saddle.  No word here
    names a kind of reaction or a role -- there is no reactant and no product,
    because which end of a mode is which is not a question about the saddle,
    and this editor is used on every sort of system.
    """
    if not end.get('ok'):
        return (f'{which} the relaxation did not finish: '
                + str(end.get('status') or 'it gave no answer.'))
    parts = []
    kcal = end.get('kcal')
    if kcal is not None:
        parts.append(f'{abs(kcal):.1f} kcal/mol '
                     + ('below' if kcal <= 0 else 'above') + ' the saddle')
    parts.append(_pieces_said(int(end.get('pieces') or 0)))
    if end.get('made'):
        parts.append('with ' + _bonds_said(symbols, end['made'], end['there'])
                     + ' that the saddle did not have')
    if end.get('broke'):
        parts.append('without ' + _bonds_said(symbols, end['broke'],
                                              end['here'])
                     + ' that the saddle had')
    if not end.get('made') and not end.get('broke'):
        parts.append('with the same bonds the saddle has')
    moved = end.get('moved')
    if moved is not None and moved == moved:
        parts.append(f'{moved:.2f} A RMSD from the saddle')
    return f'{which} it relaxed to a structure ' + ', '.join(parts) + '.'


def _modes_belong_to(modes: Dict[str, Any], xyz_text: str) -> bool:
    """Whether a kept Hessian is about this geometry.

    A tenth of a milliangstrom, which is under the four decimals a geometry is
    written out with: it separates "these coordinates, written and read back"
    from every real change.
    """
    found = _elements(xyz_text)
    here = np.asarray((modes or {}).get('angstrom'), dtype=float)
    if found is None or here.size == 0:
        return False
    there = np.asarray(found['angstrom'], dtype=float)
    if here.shape != there.shape:
        return False
    return float(np.abs(here - there).max()) <= 1e-4


def follow_the_mode_down(xyz_text: str, method: str = 'gfn2', *,
                         mode: int = 0, step: float = DESCENT_STEP,
                         charge: int = 0, uhf: int = 0,
                         solvent: Optional[str] = None, cores: int = 4,
                         timeout: Optional[float] = 120.0,
                         modes: Optional[Dict[str, Any]] = None,
                         on_stage: Optional[Callable[[str], None]] = None,
                         should_stop: Optional[Callable[[], bool]] = None,
                         ) -> Dict[str, Any]:
    """Push the structure down one mode both ways, and say where it lands.

    Returns ``{'ok', 'status', 'lines', 'ends', 'seconds', 'cm', 'order'}``.
    *ends* is the two directions, each with the structure it reached, what it
    costs against the saddle, how many pieces it is and which bonds it has
    that the saddle did not.

    This is the cheap confirmation and it is called one on purpose.  The
    rigorous version is a mass-weighted steepest descent -- an intrinsic
    reaction coordinate -- and ORCA has one; both were measured on the same
    sixteen-atom Diels-Alder saddle on this machine, and the numbers are the
    reason this is what the editor offers:

    * this, one Hessian and two relaxations: **1.0 s**, and it lands on the
      two minima themselves -- one piece with the ring closed at 1.54 A and
      -70.6 kcal/mol one way, two pieces 3.3 A apart at -6.7 the other.
    * ``! XTB2 IRC``, both directions, its own numerical Hessian: **207 s**,
      and it converges in the valley rather than at the bottom of it -- 2.83 A
      and -5.3 kcal/mol on the one side, 1.54 A and -70.0 on the other.  So it
      needs the same two relaxations afterwards to name what it reached, and
      then it agrees with this exactly.

    Two hundred seconds is not a press.  It is over the three minutes
    :func:`saddle.seconds_for` allows a saddle search itself, on the smallest
    case anybody would run, on a machine that is somebody's login node -- and
    it answers the same question.  An IRC is worth submitting as a job when
    the path itself is the result; it is not worth pressing to find out which
    two structures a saddle joins.

    It does work over ``ExtOpt``, which is how :data:`saddle.SADDLE_METHODS`
    drives g-xTB, and that was checked rather than assumed because there is a
    published interface between ORCA and a program of ours and nothing says an
    IRC uses it the way an optimiser does.  Measured on the same saddle,
    ``! ExtOpt IRC`` with ``InitHess calc_numfreq``: ORCA terminated normally
    and both directions converged, onto the same two ends -- 3.09 A and two
    pieces one way, 1.53 A and one piece the other.  It took **960 s**, 77 per
    cent of it in the external gradients, which makes the same point louder.

    *modes* is what :func:`modes_of` returned, for a caller that already paid
    for the Hessian; without it this pays for one.  They are checked against
    *xyz_text* before they are used and recomputed if they are about some
    other geometry -- a mode shape belongs to the structure it was taken at,
    and a kept Hessian outliving the structure it was taken on is the one way
    this could quietly answer about something else.
    """
    began = time.perf_counter()

    def _say(sentence: str) -> None:
        if on_stage is None:
            return
        try:
            on_stage(sentence)
        except Exception:                                   # noqa: BLE001
            # A page that cannot draw a sentence is not a reason to stop.
            pass

    if modes is not None and not _modes_belong_to(modes, xyz_text):
        modes = None
    if modes is None:
        _say('Taking a Hessian to find the mode...')
        modes = modes_of(xyz_text, method, charge=charge, uhf=uhf,
                         solvent=solvent, cores=cores)
    if not modes:
        return {'ok': False, 'lines': [], 'ends': [], 'cm': [], 'order': None,
                'seconds': time.perf_counter() - began,
                'status': ('The modes of this structure could not be '
                           'computed, so there is nothing to follow down. '
                           'It needs xtb.')}
    wrong = imaginary_among(modes['cm'])
    if not wrong:
        return {'ok': False, 'lines': [], 'ends': [], 'cm': list(modes['cm']),
                'order': 0, 'seconds': time.perf_counter() - began,
                'status': ('No mode of this structure goes the wrong way, so '
                           'there is no reaction coordinate to follow: it is '
                           'a minimum, and both directions would come '
                           'straight back to it.')}
    chosen = int(mode) if int(mode) in wrong else wrong[0]
    symbols = list(modes['symbols'])
    here = np.asarray(modes['angstrom'], dtype=float)
    way = np.asarray(modes['ways'], dtype=float)[:, chosen].reshape(here.shape)
    saddle_graph = _gfn.bond_graph(
        xyz_document(symbols, here, 'the structure the mode was taken at'))
    top = _gfn.optimize_with_gfn(
        xyz_document(symbols, here, 'the saddle'), method, charge=charge,
        uhf=uhf, solvent=solvent, optimise=False, timeout=timeout,
        should_stop=should_stop)
    summit = top.get('energy') if top.get('ok') else None

    ends: List[Dict[str, Any]] = []
    for sign, which in ((1.0, 'One way'), (-1.0, 'The other way')):
        # Between the two ways as well as inside each of them: the second is a
        # whole relaxation, and a press asking to stop should not have to wait
        # out a run it has already said it does not want.
        if should_stop is not None and should_stop():
            break
        _say(f'Pushing {which.lower()} down the mode, and relaxing...')
        pushed = xyz_document(
            symbols, displaced_along(here, way, sign * float(step)),
            'pushed down the imaginary mode')
        got = _gfn.optimize_with_gfn(
            pushed, method, charge=charge, uhf=uhf, solvent=solvent,
            optimise=True, timeout=timeout, should_stop=should_stop)
        end: Dict[str, Any] = {'which': which, 'ok': bool(got.get('ok')),
                               'status': got.get('status'),
                               'xyz': got.get('xyz'), 'here': here}
        if end['ok'] and got.get('xyz'):
            there = np.asarray(_elements(got['xyz'])['angstrom'], dtype=float)
            graph = _gfn.bond_graph(got['xyz'])
            end.update({
                'there': there,
                'graph': graph,
                'energy': got.get('energy'),
                'kcal': (None if (summit is None or got.get('energy') is None)
                         else (float(got['energy']) - float(summit))
                         * HARTREE_IN_KCAL),
                'pieces': pieces_in(got['xyz']),
                'made': sorted(graph - saddle_graph),
                'broke': sorted(saddle_graph - graph),
                'moved': turned_onto(here, there),
            })
        ends.append(end)

    lines = [_end_said(end, symbols, end['which']) for end in ends]
    both = [end for end in ends if end.get('ok')]
    same = None
    if len(both) == 2:
        first, second = both
        apart = turned_onto(first['there'], second['there'])
        cost = (None if None in (first.get('kcal'), second.get('kcal'))
                else abs(float(first['kcal']) - float(second['kcal'])))
        alike = first['graph'] == second['graph']
        same = bool(alike and apart == apart and apart < SAME_PLACE)
        if not alike:
            lines.append('The two ends do not have the same bonds, so this '
                         'mode goes between two different arrangements. '
                         'Which of them is which way round is not something '
                         'the saddle says.')
        elif same:
            lines.append(
                'Both ways came back to the same structure -- the same bonds '
                f'and {apart:.2f} A RMSD apart. Either nothing was crossed, '
                'or the push was too small to leave the saddle. Nothing here '
                'says a reaction was found.')
        elif cost is not None and cost < SAME_ENERGY:
            lines.append(
                'The two ends have the same bonds and the same energy to '
                f'within {cost:.2f} kcal/mol, and stand {apart:.2f} A RMSD '
                'apart: two arrangements of one structure rather than two '
                'different ones, which is what a passage between equivalent '
                'forms looks like.')
        else:
            lines.append(
                'The two ends have the same bonds but stand '
                f'{apart:.2f} A RMSD and '
                + ('' if cost is None else f'{cost:.1f} kcal/mol ')
                + 'apart, so the mode moves the structure without making or '
                  'breaking anything.')
    elif not both:
        lines.append('Neither direction relaxed to anything, so nothing is '
                     'claimed about what this saddle joins.')
    else:
        lines.append('Only one direction relaxed, so what this saddle joins '
                     'is half answered.')
    if len(wrong) > 1:
        lines.insert(0, (
            f'This structure has {_many(len(wrong))} modes going the wrong '
            'way, so it is not a transition state and the two ends below are '
            'only the two ways along one of them.'))
    seconds = time.perf_counter() - began
    return {'ok': bool(both), 'lines': lines, 'ends': ends, 'same': same,
            'cm': [float(one) for one in modes['cm']], 'order': len(wrong),
            'mode': chosen, 'seconds': seconds,
            'status': (f'Followed the mode at {float(modes["cm"][chosen]):.0f} '
                       f'cm-1 down both ways in {seconds:.1f} s.')}
