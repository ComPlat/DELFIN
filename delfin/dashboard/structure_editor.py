"""The structure editor, as a part rather than as a tab.

Everything the Submit tab's molecule preview can do -- picking, manipulating,
drawing, the browser force fields and the ones that run on the server, the
held internal coordinates, the polyhedra, the numbers on the atoms -- belongs
to a viewer, not to one tab. This module is where it lives so a second tab can
have the same editor over its own coordinates.

It is being assembled a piece at a time, and each piece has to leave the Submit
tab behaving exactly as it did. The first is the atom numbering, which was
reachable only from the ORCA Builder: its size control resolved the viewer as
``window._orcaBuildViewer``, one global for one tab. Numbering is a property of
a viewer, so it takes the viewer it is meant for and any tab can show it.

The numbers are a layer over the molecule and nothing else. They read the
atoms the viewer already holds, so they travel with the cores through a drag,
an optimisation or a dynamic run; switching them off removes sprites and draws
a frame, and leaves the model -- and with it the bonds -- untouched.
"""

from __future__ import annotations

import base64
import html
import importlib.resources
import json
import math
import os
import shutil
import threading
import time
from pathlib import Path

import ipywidgets as widgets

import py3Dmol

from . import climb as _climb
from . import editor_journal as _journal_mod
from . import gfn_optimize as _gfn
from . import ketcher as _ketcher
from . import mopac_optimize as _mopac
from . import saddle as _saddle
from . import scan_profile as _scan_profile
from . import separate_systems as _separate
from . import solvents as _solvents
from .input_processing import (
    append_hapto_previews_to_isomers, clean_input_data, smiles_to_xyz_isomers,
    smiles_to_xyz_quick_with_previews,
)
from .molecule_viewer import (
    apply_molecule_view_style, submit_manip_bootstrap_js, submit_manip_version,
)


#: How tall a digit comes out, in CSS pixels, per unit of the scale factor.
#: The label texture is 68 px tall and a digit fills about half of it; the
#: sprite is scaled against the drawing buffer and then divided by the device
#: pixel ratio again, so the two cancel and this holds whatever the viewer is
#: sized at. Measured in a browser at five settings: 9.5, 13, 17, 22 and 30 px
#: for 0.28, 0.38, 0.50, 0.66 and 0.86.
LABEL_PX_PER_SCALE = 34.0

#: What a size control offers: the height of a number in pixels, typed or
#: stepped. Five fixed rungs were not enough -- a crowded structure and a
#: three-atom one want different numbers, and neither wants the rung between.
LABEL_PX_DEFAULT = 6
LABEL_PX_MIN = 2
LABEL_PX_MAX = 48

#: On-screen size of the numbers: the factor the high-resolution label texture
#: is down-scaled by, so a larger number stays sharp instead of blurring.
LABEL_SCALE_DEFAULT = LABEL_PX_DEFAULT / LABEL_PX_PER_SCALE

#: How often the labels may be redrawn while something is running, in seconds.
#:
#: They travel by ``run_js``, which is the channel this editor keeps for
#: things that happen when a user presses something -- the frames and the
#: thermal wall each have a widget field of their own precisely because
#: run_js clears its output before displaying, so a message sent many times a
#: second is mostly overwritten before the page draws it.
#:
#: Labels can live with that where a frame cannot: each repaint is the whole
#: state rather than one step of a sequence, so the last to land is the right
#: one.  What is left is the cost of the script itself, and four times a
#: second is already faster than a charge visibly changes -- measured on an
#: ammonia borane pulled apart under GFN2, the nitrogen went -0.25, -0.26,
#: -0.27, -0.29, -0.33, -0.40 over six answers.
_LABEL_REPAINT_INTERVAL = 0.25


#: One layout for every coordinate this editor writes: the element in five
#: columns, then three fields of twenty-four with fourteen decimals.  It is
#: xtb's own, which is what most of these coordinates already are, and the
#: columns line up whether or not a number carries a minus sign.
#:
#: There were three layouts before, and which one the box ended up in said
#: where the coordinates had come from rather than anything about the
#: molecule: fourteen decimals from xtb, eight from the frame a stopped run
#: was showing, and six from the browser, whose model is serialised with
#: toFixed(6).  Two structures that differed only in that were two different
#: histories, and reading which was which off the decimal count is not
#: something a user should have to do.
_HARTREE_TO_KCAL = 627.5094740631

#: How long a structure is given to do something, in seconds.
#:
#: "Possible at this temperature" has two halves and this is the second: a
#: barrier is not possible or impossible, it is a waiting time.  At 298 K one
#: of 20 kcal/mol takes a minute and one of 25 takes four days, so the ceiling
#: depends on how long anyone is watching -- 17.5 kcal/mol for a second, 22.3
#: for an hour, 27.7 for a year.
#:
#: It was a control, and it did not earn the place.  Between a second and a
#: year the ceiling moves by ten kcal/mol, while the distance between
#: chemistry and nonsense is twenty against a hundred and more -- so the knob
#: changed the answer far less than it cost to understand, and it was asked
#: about twice.  An hour is a reaction set up in a flask, which is what this
#: is for.
#:
#: Asked again once a refusal started answering in both directions, and the
#: answer came out the same.  What a window control would be for is a user
#: saying "I would wait longer than that", and the two sentences a refusal now
#: carries answer that without anything being set: the waiting time is quoted
#: at the temperature that *is* set, so a barrier the hour refuses already
#: reads as "about 4 d", "about 50 years" or "longer than the universe has
#: existed", and the temperature it wants is quoted for the window.  A knob
#: that only moved the line between those two would be a third way of asking a
#: question now answered twice over without it.
#:
#: Maeda's advice for AFIR is to run the permissive end -- the gamma he
#: recommends is the ten-day 106.9 kJ/mol and not the one-hour 93.3 -- and it
#: is advice about *searching*: a search that misses a path has found nothing,
#: while one that finds too much is sorted out afterwards.  This is not a
#: search.  It decides what stays in the box, and there the two errors are not
#: the same size: a budget that allows what the temperature will not hands
#: back a structure with nothing anywhere saying it is impossible, and the
#: user goes on from it.  The strict end is right here for the same reason the
#: permissive end is right there.
_THERMAL_SECONDS = 3600.0
_XYZ_ELEMENT_COLUMNS = 5
_XYZ_NUMBER_COLUMNS = 24
_XYZ_DECIMALS = 14


def xyz_line(symbol, x, y, z):
    """One atom, in the layout every coordinate box in the editor uses."""
    return (f'{str(symbol):<{_XYZ_ELEMENT_COLUMNS}}'
            f'{float(x):>{_XYZ_NUMBER_COLUMNS}.{_XYZ_DECIMALS}f}'
            f'{float(y):>{_XYZ_NUMBER_COLUMNS}.{_XYZ_DECIMALS}f}'
            f'{float(z):>{_XYZ_NUMBER_COLUMNS}.{_XYZ_DECIMALS}f}')


def xyz_body(lines):
    """Coordinate lines re-laid-out, and left alone if they cannot be read.

    A line this cannot parse is passed through untouched rather than dropped:
    the box is the user's, and a comment or a fifth column someone is relying
    on is not this function's to throw away.
    """
    out = []
    for line in lines:
        if not str(line).strip():
            continue
        parts = str(line).split()
        if len(parts) < 4:
            out.append(str(line))
            continue
        try:
            out.append(xyz_line(parts[0], parts[1], parts[2], parts[3]))
        except (TypeError, ValueError):
            out.append(str(line))
    return out


def xyz_document(lines, comment):
    """A whole .xyz: the count, the comment, and the atoms in one layout."""
    body = xyz_body(lines)
    return f'{len(body)}\n{comment}\n' + '\n'.join(body)


#: Comments this editor writes about where a geometry came from, as opposed to
#: a name off a file or a note the user typed.  Only these may be replaced when
#: the coordinates underneath them are no longer the ones they were written
#: for; everything else in that line belongs to the user.
#:
#: A claim that is not in this list is treated as the user's and carried over
#: whatever happens to the coordinates beneath it, which turns it into a lie
#: nobody can see.  Measured: a drag that the budget rolled back left "Past
#: the budget: back to the last structure that was inside it" in the box, the
#: next message from the page wrote its own coordinates under that line, and
#: what the user was looking at was a torn ethane at +141.2 kcal/mol wearing
#: the sentence that says it had been taken back.  So every line this file
#: writes belongs here, and adding one to the writing means adding it here in
#: the same breath.
_EDITOR_COMMENTS = (
    'optimised in delfin viewer',
    'edited in delfin viewer',
    'settled with ',
    'stopped at the frame on screen',
    'stopped where you took hold',
    'relaxed, and the budget measured from here',
    'within the budget',
    'past the budget',
    'back to the last structure that was measured and allowed',
    'kept: the bonding would have changed',
    'scanned',
    'driven until the bonds were made and broken',
    'where the saddle search got to',
    'where the hand left it',
    'optimised to ',
    'climbed towards ',
    'where the chain got to',
    'where the band got to',
    'climbed to ',
    'from a path, optimised to ',
    'from a band, optimised to ',
    'optimised to a transition state',
    'estimated transition state, from the path',
    'delfin drag-end',
    'delfin drag-follow',
    'from the delfin viewer',
)


def _is_editor_comment(line):
    """Whether that comment line is this editor's own claim about a geometry."""
    text = str(line or '').strip().lower()
    return any(text.startswith(one) for one in _EDITOR_COMMENTS)


def fixed_atoms_note(held):
    """What a force field running here did with the held values, in its terms.

    The counterpart of :func:`gfn_optimize.held_note` and
    :func:`mopac_optimize.freeze_note`, and it reads the same shape MOPAC's
    does, because the two engines can do the same one thing: RDKit's UFF takes
    ``AddFixedPoint`` and Open Babel a list of fixed atoms, and neither has a
    restraint.  So a fix is met by holding the atoms that name it still --
    which is more than was asked, since those atoms also stop turning and
    moving -- and a pull cannot be said at all.

    Which of the two happened has to be said, because the alternative is what
    this branch of Optimise did before it was handed anything: measured on
    ethane with the bonding pinned, C0-H2 pulled out to 1.700 A and held
    exactly came back from one press at 1.1104 A, under a status line reading
    only "Optimised with UFF."  A held value that is silently given up makes
    the result an answer to a question nobody asked.
    """
    said = []
    if held['held']:
        said.append(
            f"{held['held']} held value(s) kept by fixing the "
            f"{len(held['frozen'])} atom(s) they name, where they stand -- "
            'the force field fixes atoms, not the value between them, so '
            'those atoms also stop turning and moving')
    if held.get('every_atom'):
        # A small molecule runs out of atoms quickly: an angle held on a water
        # names all three of them.  The press then cannot move anything, and
        # saying "Optimised" over a geometry nothing happened to would be the
        # same silence this function exists to end, one step further along.
        said.append('and that is every atom in the structure, so there was '
                    'nothing left to relax -- release one of them to let the '
                    'rest of the molecule move')
    if held['pulls']:
        said.append(f"{held['pulls']} pull(s) not honoured -- an atom here is "
                    'held or free, so there is no value to negotiate with; '
                    'hold them as fix, or optimise under a GFN method')
    if held['dropped']:
        said.append(f"{len(held['dropped'])} held value(s) dropped -- they "
                    'name atoms this structure does not have')
    return (' ' + '; '.join(said) + '.') if said else ''


#: Boltzmann, Planck and the gas constant, in the units this file speaks:
#: kcal/mol for energies, kelvin for temperature, seconds for time.
_BOLTZMANN_SI = 1.380649e-23          # J/K
_PLANCK_SI = 6.62607015e-34           # J s
_GAS_CONSTANT = 1.987204259e-3        # kcal/(mol K)
#: 13.8 thousand million years, in seconds.  The longest waiting time worth
#: printing as a number: past it every barrier reads the same, and a sentence
#: says what a mantissa and an exponent do not.
_UNIVERSE_SECONDS = 4.35e17


def thermal_ceiling(temperature, seconds):
    """The highest barrier a structure crosses at *temperature* within
    *seconds*, in kcal/mol.

    Eyring turned around.  A rate is k = (kB T / h) exp(-dG/RT), and a barrier
    is crossed about once in 1/k -- so asking how high it may be to be crossed
    within a given time gives

        dG = R T ln(kB T tau / h)

    which is the whole content of "possible at this temperature".  It is not
    one number: it depends on how long the user is prepared to wait, and the
    two make very different chemistry.  At 298 K a second buys 17.5 kcal/mol,
    an hour 22.3 and a year 27.7; at 773 K an hour buys 59.3.

    Measured against it on a benzene with GFN2, the ring bond stretched with
    everything else relaxed: 1.55 A costs 9.8 kcal/mol, 1.65 A costs 21.3 and
    1.75 A costs 34.7.  So at room temperature and an hour the bond can be
    pulled to about 1.66 A and no further -- the ring cannot be torn open, and
    that falls out of the energy rather than out of a rule about aromatics.
    At 773 K the same ceiling reaches 1.9 A, which is the temperature at which
    benzene really does come apart.
    """
    T = max(1.0, float(temperature))
    # A picosecond is the floor, and it is not arithmetic hygiene: below about
    # a tenth of that a molecule has not finished one vibration, so there is no
    # crossing to speak of and no barrier height that means anything.  It also
    # keeps the logarithm away from zero.
    tau = max(1e-12, float(seconds))
    inside = _BOLTZMANN_SI * T * tau / _PLANCK_SI
    if inside <= 1.0:
        return 0.0
    return _GAS_CONSTANT * T * math.log(inside)


def thermal_temperature(kcal, seconds=_THERMAL_SECONDS):
    """The temperature at which a barrier of *kcal* is crossed within *seconds*.

    :func:`thermal_ceiling` turned around again.  A drag has to be told what
    it may do at the temperature it was given, so there the question runs one
    way; a scan does not -- it walks the whole path and can then answer the
    question a chemist actually asks, which is *how hot*.  "+29 kcal/mol"
    means nothing until it is "and that wants 440 K to happen within the
    hour", which is the difference between a number and an experiment.

    The ceiling rises with temperature everywhere above a kelvin (T times the
    logarithm of T, both increasing), so bisection finds the answer and there
    is no need for the Lambert function this inverts to.  Returns ``None``
    when no temperature under 5000 K will do it -- past that a molecule is not
    a molecule, and saying so is better than printing a number.
    """
    try:
        wanted = float(kcal)
    except (TypeError, ValueError):
        return None
    if wanted <= thermal_ceiling(1.0, seconds):
        return 1.0
    low, high = 1.0, 5000.0
    if thermal_ceiling(high, seconds) < wanted:
        return None
    for _ in range(80):
        middle = 0.5 * (low + high)
        if thermal_ceiling(middle, seconds) < wanted:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


#: What the ceiling is, and what is held against it.
#:
#: The ceiling out of :func:`thermal_ceiling` is a free energy of activation:
#: the dG that Eyring inverts.  What a drag is priced with is an electronic
#: energy -- the difference between two relaxed points on the surface the
#: method describes -- and the two are not the same quantity.  They differ by
#: T*dS, and the size of that is not a detail that can be waved away:
#:
#:   * While the drag leaves the structure in as many separate pieces as it
#:     found it, dS is small.  A torsion, an angle, a ring turning over: what
#:     is lost from the vibrations at one end is found at the other, and the
#:     two numbers agree closely.  Measured here under GFN2 at 298.15 K, the
#:     same deformation priced both ways: ethane turned to eclipsed, +2.592
#:     electronic against +2.568 free -- 0.02 kcal/mol apart on a barrier of
#:     2.6; a benzene ring bond stretched to 1.62 A, +17.98 against +15.71;
#:     to 1.72 A, +30.58 against +27.68.  Under three kcal/mol on every one
#:     of them, which is well inside what the method itself is worth, and no
#:     verdict anywhere near the ceiling changes hands.
#:   * Where the drag changes how many pieces there are, dS is large and it
#:     has a sign.  Taking something apart releases translation and rotation,
#:     so T*dS is of order ten kcal/mol at room temperature and the electronic
#:     price is *too strict* -- the wall refuses something the temperature
#:     would in fact pay for.  Bringing two things together is the same number
#:     the other way round, and there it is too lenient.  Measured on a
#:     borazane with the B-N pulled out, GFN2 at 298.15 K: at 2.5 A +17.71
#:     electronic against +12.77 free, at 3.5 A +21.20 against +14.21, and
#:     once the two are apart at 6 A +22.51 against +12.26.  Ten and a quarter
#:     kcal/mol at the end of it, which is the "order ten" above arrived at by
#:     measurement, and it is the one place a verdict changes hands: against
#:     an hour at 298 K the electronic price is past the 22.3 and the free one
#:     is not, and where the electronic price asks for 301 K the free one asks
#:     for 167.
#:
#: Nothing here corrects for it, and that is a decision rather than an
#: oversight.  The exact answer is a Hessian, and a Hessian per drag step
#: cannot be afforded: measured here under GFN2, a numerical Hessian costs
#: 0.69 s on 21 atoms and 3.90 s on 62, against 58 and 321 ms for the single
#: point beside it -- twelve times the cost at both sizes, and it would have
#: to be taken on every answer of a control that reports several times a
#: second.  A cheap correction would have to guess dS, and the only guess
#: available is how many pieces the structure is in -- which is read off a
#: distance threshold, and that threshold is the one already written down as
#: flickering in a crowded coordination sphere.  It would also depend on a
#: standard state nothing in this editor knows, and it is exactly the case
#: where the method itself is least reliable.  A number invented there would
#: be worse than the gap it papered over.
#:
#: There is one published prescription for exactly this case, and measured, it
#: does not answer here.  xtb's ``--bhess`` takes a single-point Hessian for a
#: geometry that is not a stationary point -- Spicher and Grimme, J. Chem.
#: Theory Comput. 2021, 17, 1701 -- by optimising under an RMSD bias towards
#: the structure it was handed, so that the soft and imaginary modes of a
#: point nothing has relaxed do not wreck the entropy.  A dragged structure is
#: exactly a geometry that is not a stationary point, so it reads like the
#: answer to the paragraph above.  It is not, and the reason is that its bias
#: is sized in RMSD -- the target is 0.10 A -- while a drag is small in RMSD
#: and large in energy.  Measured under GFN2: a benzene ring bond stretched to
#: 1.72 A is +30.6 kcal/mol and 0.094 A of RMSD, xtb reads that as already
#: inside the target, prints ``final kpush: -0.000000``, optimises freely back
#: to the ring, and reports its TOTAL FREE ENERGY for relaxed benzene --
#: -0.0003 kcal/mol against the anchor.  Priced that way the wall would refuse
#: nothing whatever.  Held at the coordinate the hand is holding it does keep
#: the geometry, and then the restraint's own curvature is in the frequencies:
#: the same case comes back at +0.4 where an unbiased Hessian says -2.9, which
#: is the hold and not the chemistry.  What prices a dragged geometry is a
#: plain ``--hess`` on it, and that is what a scan set to G already runs.
#:
#: So the gap is said rather than filled, and a scan answers it properly: with
#: its energy set to G it takes three Hessians -- the start, the highest point
#: and the end -- and its verdict is then a free energy against a free energy.
_THERMAL_QUANTITY_SHORT = (
    'The ceiling is a free energy; a drag is priced with an electronic '
    'energy. While the structure stays in as many pieces as it started in the '
    'two agree to under 3 kcal/mol -- measured on a torsion and on stretched '
    'ring bonds. They part company where a drag changes how many separate '
    'pieces there are: about ten kcal/mol at room temperature, strict for '
    'taking something apart and lenient for putting it together. Run a scan '
    'with its energy set to G for the free-energy answer.'
)


def scale_for_px(px):
    """The scale factor that makes a digit *px* pixels tall."""
    try:
        wanted = float(px)
    except (TypeError, ValueError):
        return LABEL_SCALE_DEFAULT
    wanted = max(LABEL_PX_MIN, min(LABEL_PX_MAX, wanted))
    return wanted / LABEL_PX_PER_SCALE


#: What the editor knows about the structure it is working on, as opposed to
#: what it knows about itself. A tab that holds several structures and lets the
#: user step between them puts these aside for the one being left and hands
#: back the ones for the one being shown -- otherwise coming back to a
#: structure means perceiving it again from its coordinates, and an atom that
#: was pulled away from its neighbour comes back with the bond gone.
STRUCTURE_MEMORY_KEYS = (
    'perceived', 'perceived_for', 'bond_edits', 'hand_bonds', 'hyb_overrides',
    'constraints', 'poly_applied', 'poly_metal', 'poly_assignment',
    'poly_arrangements', 'poly_arrangement_index', 'history', 'structure_undo',
    'pristine_coords', 'gfn_topology', 'gfn_topology_source',
    # The charges belong to the structure they were computed on, and the tab
    # that steps between structures is exactly the one that would get this
    # wrong: a set of isomers has the same element column, so the fingerprint
    # that keeps one molecule's charges off another's atoms cannot tell two
    # isomers apart. Put aside with the structure being left, they are simply
    # not there for the one being shown until it has an answer of its own.
    'atom_charges', 'atom_charges_method', 'atom_charges_for',
    # What a saddle search found here, and the Hessian that was paid for to
    # draw and follow its mode. Both belong to one structure and to no other
    # -- they are checked against the geometry before they are believed --
    # and both are seconds of xtb, which is worth not spending twice because
    # somebody stepped to the next block and back.
    'saddle_found', 'saddle_modes',
)


def _atom_numbers_js():
    """The layer itself: ``window.__delfinAtomNumbers``.

    ``set(viewer, on, scale, texts)`` switches the labels on or off,
    ``refresh`` brings them back onto the atoms after those have moved,
    ``setScale`` resizes what is already there.

    *texts* is what each label says, one to an atom, and null means the atom's
    number.  The layer was written for the numbers and is not about them: what
    it really does is hold a sprite on an atom while that atom moves, hide it
    when something is in front of it, and keep it the size the toolbar asked
    for -- all of which is exactly as true of a partial charge as of an index.
    A second layer for charges would have been the same six hundred lines with
    a different string in one place, and would have drifted from this one.
    """
    return (
        "window.__delfinAtomNumbers=(function(){\n"
        "  var RAD=0.5, RAD2=RAD*RAD, EPS=0.05, DEPTH_SIGN=1;\n"
        "  var DEFAULT_SCALE=%.3f;\n"
        # Every model in the viewer is numbered, each from its own zero: the
        # ORCA Builder's overlay shows two structures at once and both want
        # the numbering of the file they came from.
        "  function modelsOf(v){\n"
        "    var out=[];\n"
        "    try{\n"
        "      var ms=v.models;\n"
        "      if(ms&&ms.length){for(var i=0;i<ms.length;i++){if(ms[i])out.push(ms[i]);}}\n"
        "      else{var m=v.getModel();if(m)out.push(m);}\n"
        "    }catch(e){}\n"
        "    return out;\n"
        "  }\n"
        "  function atomCount(v){\n"
        "    var ms=modelsOf(v),n=0;\n"
        "    for(var i=0;i<ms.length;i++){\n"
        "      try{n+=(ms[i].selectedAtoms({})||[]).length;}catch(e){}\n"
        "    }\n"
        "    return n;\n"
        "  }\n"
        # Taking the sprites out by hand rather than with removeLabel: that
        # one draws a frame per label, so a hundred numbers meant a hundred
        # renders. Nothing here reaches the model, which is the point -- the
        # numbers must not cost the user their bonds.
        "  function clear(v){\n"
        "    if(!v)return;\n"
        "    var L=v.__delfinLbls||[];\n"
        "    for(var i=0;i<L.length;i++){\n"
        "      var lab=L[i]&&L[i].l;if(!lab)continue;\n"
        "      try{if(lab.sprite&&v.modelGroup)v.modelGroup.remove(lab.sprite);}catch(e){}\n"
        "      try{var k=(v.labels||[]).indexOf(lab);if(k>=0)v.labels.splice(k,1);}catch(e){}\n"
        "      try{lab.dispose();}catch(e){}\n"
        "    }\n"
        "    v.__delfinLbls=[];v.__delfinProj=[];\n"
        "  }\n"
        "  function hook(v){\n"
        "    if(v.__delfinLabelHooked)return;\n"
        "    v.__delfinLabelHooked=true;\n"
        "    var settle=function(){\n"
        "      try{if(refresh(v)&&typeof v.render==='function')v.render();}catch(e){}\n"
        "    };\n"
        "    var hs=v.__delfinInteractionEndHandlers||(v.__delfinInteractionEndHandlers=[]);\n"
        "    hs.push(settle);\n"
        "    var before=v.__delfinCleanup;\n"
        "    v.__delfinCleanup=function(){\n"
        "      var pos=hs.indexOf(settle);if(pos>=0)hs.splice(pos,1);\n"
        "      v.__delfinLbls=[];v.__delfinProj=[];\n"
        "      if(typeof before==='function'){try{before();}catch(e){}}\n"
        "    };\n"
        "  }\n"
        "  function build(v,scale,texts){\n"
        "    if(!v||typeof v.addLabel!=='function')return 0;\n"
        "    clear(v);\n"
        "    if(scale!=null&&isFinite(+scale))v.__delfinLabelScale=+scale;\n"
        # Kept on the viewer, so a rebuild -- which happens by itself whenever
        # the atom count changes under a running drag -- says the same thing
        # it said before rather than falling back to the numbers.
        "    if(texts!==undefined)v.__delfinLabelTexts=texts;\n"
        "    var T=v.__delfinLabelTexts||null;\n"
        # A list of values that no longer matches the structure draws nothing
        # at all, rather than drawing the ones it still has and numbering the
        # rest. Half a set of charges and half a set of indices, in the same
        # typeface on the same atoms, is the one outcome worth writing code to
        # prevent: the next answer brings a fresh list and they come back.
        "    if(T&&T.length!==atomCount(v)){\n"
        "      v.__delfinLabelTexts=null;v.__delfinLbls=[];v.__delfinProj=[];\n"
        "      return 0;\n"
        "    }\n"
        # alignment:center anchors the text box on its centre, so the number
        # stays on the atom centre at every zoom (default corner-anchoring
        # drifts aside as atoms shrink). Fall back to the string form if the
        # enum is unavailable.
        "    var al=(window.$3Dmol&&$3Dmol.SpriteAlignment&&$3Dmol.SpriteAlignment.center)\n"
        "      ?$3Dmol.SpriteAlignment.center:'center';\n"
        "    var ms=modelsOf(v),L=[],n=0;\n"
        "    for(var mi=0;mi<ms.length;mi++){\n"
        "      var atoms=[];try{atoms=ms[mi].selectedAtoms({})||[];}catch(e){atoms=[];}\n"
        "      for(var i=0;i<atoms.length;i++,n++){\n"
        # fontSize 48 is a HIGH-RES texture kept sharp; the sprite is then
        # down-scaled in refresh() so the number appears small and crisp.
        # The fourth argument tells 3Dmol not to draw a frame per label.
        # The number is the atom's place in the model it belongs to, which is
        # what the ORCA Builder's two-structure overlay needs; a text handed
        # in is read off one list for the whole viewer, because whoever
        # computed it computed it for one structure.
        "        var a=atoms[i],lab=null;\n"
        "        var say=(T&&T[n]!=null)?String(T[n]):String(i);\n"
        "        try{lab=v.addLabel(say,{position:{x:a.x,y:a.y,z:a.z},\n"
        "          fontSize:48,fontColor:'black',alignment:al,\n"
        "          showBackground:false,inFront:true},undefined,true);}catch(e){lab=null;}\n"
        "        if(lab)L.push({a:a,l:lab});\n"
        "      }\n"
        "    }\n"
        "    v.__delfinLbls=L;v.__delfinProj=[];\n"
        "    if(L.length)hook(v);\n"
        "    return L.length;\n"
        "  }\n"
        # Called on every frame the molecule is drawn in, so a number is where
        # its atom is at that moment -- during a drag, an optimisation or a
        # dynamic run as much as after one. The atom object is held, not a
        # copy of its coordinates, which is why following it costs nothing.
        "  function refresh(v){\n"
        "    if(!v||v.__delfinDisposed)return false;\n"
        "    var L=v.__delfinLbls||[],m=L.length;\n"
        "    if(!m)return false;\n"
        "    if(atomCount(v)!==m&&!v.__delfinLabelRebuilding){\n"
        "      v.__delfinLabelRebuilding=true;\n"
        "      try{m=build(v,v.__delfinLabelScale);}finally{v.__delfinLabelRebuilding=false;}\n"
        "      L=v.__delfinLbls||[];\n"
        "      if(!m)return false;\n"
        "    }\n"
        "    var w;try{w=v.getView();}catch(e){return false;}\n"
        "    if(!w||w.length<8)return false;\n"
        "    var x=w[4],y=w[5],z=w[6],q=w[7];\n"
        "    var r11=1-2*(y*y+z*z), r12=2*(x*y-q*z), r13=2*(x*z+q*y);\n"
        "    var r21=2*(x*y+q*z), r22=1-2*(x*x+z*z), r23=2*(y*z-q*x);\n"
        "    var r31=2*(x*z-q*y), r32=2*(y*z+q*x), r33=1-2*(x*x+y*y);\n"
        "    var P=v.__delfinProj||(v.__delfinProj=[]);\n"
        "    var grid=Object.create(null);\n"
        "    for(var i=0;i<m;i++){\n"
        "      var a=L[i].a,sp=L[i].l&&L[i].l.sprite;\n"
        "      if(sp&&sp.position&&typeof sp.position.set==='function')\n"
        "        sp.position.set(a.x,a.y,a.z);\n"
        "      var p=P[i]||(P[i]=[0,0,0]);\n"
        "      p[0]=r11*a.x+r12*a.y+r13*a.z;\n"
        "      p[1]=r21*a.x+r22*a.y+r23*a.z;\n"
        "      p[2]=DEPTH_SIGN*(r31*a.x+r32*a.y+r33*a.z);\n"
        "      var gx=Math.floor(p[0]/RAD),gy=Math.floor(p[1]/RAD);\n"
        "      var key=gx+':'+gy,bucket=grid[key]||(grid[key]=[]);bucket.push(i);\n"
        "    }\n"
        # A number sits on top of its own atom (inFront) so it is never hidden
        # by it; one that is behind another atom is hidden here instead. The
        # grid keeps that to nearby candidates rather than every pair.
        "    for(var b0=0;b0<m;b0++){ var occ=false, pa=P[b0];\n"
        "      var agx=Math.floor(pa[0]/RAD),agy=Math.floor(pa[1]/RAD);\n"
        "      for(var ox=-1;ox<=1&&!occ;ox++){for(var oy=-1;oy<=1&&!occ;oy++){\n"
        "        var near=grid[(agx+ox)+':'+(agy+oy)]||[];\n"
        "        for(var n=0;n<near.length;n++){var b=near[n];if(b===b0)continue;var pb=P[b];\n"
        "          if(pb[2]>pa[2]+EPS){var ex=pa[0]-pb[0],ey=pa[1]-pb[1];\n"
        "            if(ex*ex+ey*ey<RAD2){occ=true;break;}}\n"
        "        }\n"
        "      }}\n"
        "      var s=L[b0].l&&L[b0].l.sprite; if(s) s.visible=!occ;\n"
        "    }\n"
        # The zoom factor f measures pixels-per-Angstrom in the SCREEN PLANE
        # at the rotation centre (world (0,0,rz), depth invariant under
        # rotation) along the screen-right model direction, so f tracks the
        # mouse-wheel zoom only, NOT rotation. DISP is the on-screen size the
        # toolbar's size control asked for.
        "    var DISP=(+v.__delfinLabelScale)||DEFAULT_SCALE, f=1;\n"
        "    if(typeof v.modelToScreen==='function'&&m>0){ try{\n"
        "      var pcx=-w[0],pcy=-w[1],pcz=-w[2];\n"
        "      var q1=v.modelToScreen({x:pcx,y:pcy,z:pcz});\n"
        "      var q2=v.modelToScreen({x:pcx+r11,y:pcy+r12,z:pcz+r13});\n"
        "      var ppa=Math.sqrt((q1.x-q2.x)*(q1.x-q2.x)+(q1.y-q2.y)*(q1.y-q2.y));\n"
        "      if(ppa>0){ var base=v.__delfinPPA0||(v.__delfinPPA0=ppa);\n"
        "        f=ppa/base; if(f<0.4)f=0.4; if(f>2.5)f=2.5; }\n"
        "    }catch(e){} }\n"
        "    for(var g=0;g<m;g++){ var sg=L[g].l&&L[g].l.sprite;\n"
        "      if(!sg||!sg.scale) continue;\n"
        "      if(!L[g].s0&&sg.scale.x>0) L[g].s0=[sg.scale.x,sg.scale.y];\n"
        "      if(L[g].s0){ sg.scale.x=L[g].s0[0]*DISP*f; sg.scale.y=L[g].s0[1]*DISP*f; }\n"
        "    }\n"
        "    return true;\n"
        "  }\n"
        "  function draw(v){\n"
        "    if(typeof v.render==='function'){try{v.render();}catch(e){}}\n"
        "  }\n"
        "  function setScale(v,scale){\n"
        "    if(!v||!isFinite(+scale))return false;\n"
        "    v.__delfinLabelScale=+scale;\n"
        "    if(!(v.__delfinLbls||[]).length)return false;\n"
        "    if(refresh(v))draw(v);\n"
        "    return true;\n"
        "  }\n"
        "  function set(v,on,scale,texts){\n"
        "    if(!v)return 0;\n"
        "    if(!on){\n"
        "      var had=(v.__delfinLbls||[]).length;\n"
        "      clear(v);v.__delfinLabelTexts=null;if(had)draw(v);\n"
        "      return 0;\n"
        "    }\n"
        "    var n=build(v,scale,texts);\n"
        "    if(n){refresh(v);draw(v);}\n"
        # A build that came to nothing still has to be drawn: what it cleared
        # is on screen until something renders over it.
        "    else draw(v);\n"
        "    return n;\n"
        "  }\n"
        "  return {set:set,build:build,clear:clear,refresh:refresh,setScale:setScale};\n"
        "})();"
    ) % LABEL_SCALE_DEFAULT


def atom_numbers_js():
    """The layer, installed once per page and only if it is not there yet."""
    return 'if(!window.__delfinAtomNumbers){\n' + _atom_numbers_js() + '\n}'


def show_atom_numbers_js(var='viewer', on=True, scale=None, texts=None):
    """Label the atoms of the viewer held in the JS variable *var*.

    The atoms are read off the model the viewer already has, so this is the
    same call whether a molecule was just rendered or has been on screen for
    an hour, and it never needs the coordinates handed to it a second time.

    With *texts* None the labels are the atom numbers, which is what this has
    always drawn.  Given a list -- one entry per atom, in the order of the
    coordinates -- each atom says what it was given instead.  Nothing here
    knows or cares what the values are; see :func:`atom_charge_texts` for the
    one thing that is drawn this way and for what it costs.
    """
    size = LABEL_SCALE_DEFAULT if scale is None else float(scale)
    said = 'null' if texts is None else json.dumps(
        [str(one) for one in texts])
    return (
        atom_numbers_js()
        + '\nwindow.__delfinAtomNumbers.set(%s,%s,%.3f,%s);'
        % (var, 'true' if on else 'false', size, said)
    )


#: How many decimals a partial charge is drawn to.
#:
#: Two.  A charge is not a measurement -- four methods give a methane carbon
#: -0.153, -0.130, -0.359 and -0.092, so the third decimal is a property of
#: the Hamiltonian rather than of the molecule -- and two is what fits on an
#: atom without the numbers running into each other.
CHARGE_DECIMALS = 2


def atom_charge_texts(charges, decimals=CHARGE_DECIMALS):
    """Partial charges as the strings that go on the atoms, or None.

    Signed always: the sign is the whole of what a chemist reads off a
    structure at a glance, and a charge of +0.15 written as 0.15 reads as a
    magnitude.
    """
    if not charges:
        return None
    said = []
    for one in charges:
        try:
            said.append(f'{float(one):+.{int(decimals)}f}')
        except (TypeError, ValueError):
            return None
    return said



from delfin.cli_manta import _CHAMPION_FLAGS as _MANTA_CHAMPION_FLAGS
_MANTA_OPT_TOPN = 10
_MANTA_OPT_WORKERS = 4
_MANTA_GIF_DATA_URI_CACHE = None

def _manta_gif_data_uri():
    """Return the packaged MANTA loading animation as a data URI (cached).

    Mirrors the dashboard logo loader: read the gif shipped under
    ``delfin/logo`` and inline it as a base64 ``data:`` URI so it renders
    inside the molecule viewer without needing a served static path.
    Returns ``''`` when the asset is unavailable.
    """
    global _MANTA_GIF_DATA_URI_CACHE
    if _MANTA_GIF_DATA_URI_CACHE is not None:
        return _MANTA_GIF_DATA_URI_CACHE
    uri = ''
    try:
        gif = importlib.resources.files('delfin').joinpath(
            'logo', 'MANTA_readme_demo.gif')
        if gif.is_file():
            data = gif.read_bytes()
            if data:
                uri = 'data:image/gif;base64,' + base64.b64encode(data).decode('ascii')
    except Exception:
        uri = ''
    _MANTA_GIF_DATA_URI_CACHE = uri
    return uri

def _manta_best_env(charge, construction="champion", method="gfn2", rank=True):
    """Env for the chosen construction config + GFN2 energy ranking, GFN2 charge
    from the SMILES. construction: 'champion' (full SHIP-31 rich + KAPPA4 reach,
    DEFAULT/maximum richness) | 'builder' (lean core + reach) | 'default' (legacy)."""
    env = {"DELFIN_GFNFF_CHARGE": str(int(charge))}
    if rank:
        env["DELFIN_FFFREE_GFNFF_RANK"] = "1"
        env["DELFIN_CONF_RANK_METHOD"] = method
    if construction != "default":
        env["DELFIN_FFFREE_BUILDER"] = "1"
        env["DELFIN_FRAME_RANK_FIX"] = "1"
        env["DELFIN_CHIRAL_ENUM"] = "1"        # Λ/Δ enantiomer enumeration (>=2 chelate pairs)
        if construction == "champion":
            for _f in _MANTA_CHAMPION_FLAGS:   # de-bloated set (KAPPA4 included; CONF_ENERGY_RANK dropped)
                env["DELFIN_FFFREE_" + _f] = "1"
        else:  # builder = lean core + reach
            env["DELFIN_FFFREE_KAPPA4"] = "1"
            env["DELFIN_FFFREE_SIGMA_ENSEMBLE"] = "1"
            env["DELFIN_FFFREE_CONF_ENERGY_RANK"] = "1"
    return env

def _manta_rank_only(isomers, charge, method="gfn2", spin="auto"):
    """RANK the manifold by xtb SINGLE-POINT energy: reorder best (lowest-energy) first WITHOUT
    changing any geometry.  Each item is ``(xyz_string, num_atoms, label)``; the emitted structures
    stay byte-identical to construction — only their ORDER changes.  spin='auto' -> parity-correct
    uhf per structure (even electrons=singlet, odd=doublet); a fixed multiplicity sets uhf=mult-1.
    Best-effort: any structure whose energy eval fails sinks to the end keeping its geometry.
    Returns the list unchanged if xtb is unavailable or there is nothing to reorder."""
    if not isomers or len(isomers) < 2:
        return isomers
    try:
        from delfin.manta import _gfnff_rank as _gff
    except Exception:
        return isomers
    if not _gff.available():
        return isomers
    import concurrent.futures as _cf

    def _uhf_for(xyz):
        if str(spin) != "auto":
            return max(0, int(spin) - 1)
        try:
            return _gff._n_electrons(xyz, int(charge)) % 2   # parity-correct ground-state multiplicity
        except Exception:
            return 0

    def _energy_one(item):
        xyz = item[0]
        try:
            return _gff.gfnff_energy(xyz, charge=int(charge), uhf=_uhf_for(xyz), method=method)
        except Exception:
            return None
    _max_workers = max(1, min(len(isomers), (os.cpu_count() or 4)))
    try:
        with _cf.ThreadPoolExecutor(max_workers=_max_workers) as ex:
            energies = list(ex.map(_energy_one, isomers))
    except Exception:
        return isomers
    # Ascending by energy; failed evals (None) sink to the end preserving their relative order.
    order = sorted(range(len(isomers)),
                   key=lambda i: (energies[i] is None, energies[i] if energies[i] is not None else 0.0, i))
    return [isomers[i] for i in order]

def _manta_opt_top(isomers, charge, topn=None, method="gfn2", spin="auto"):
    """Geometry-optimize the top-N ranked isomers in parallel (laptop-bounded),
    replace their geometry + label, re-sort the optimized head by opt energy. The
    opt ``method`` FOLLOWS the Rank selection (gfn2/gfnff/gfn1/gfn0) so one switch
    controls both. Each item is ``(xyz_string, num_atoms, label)``. Best-effort:
    any structure whose optimization fails keeps its unrelaxed geometry.
    ``topn`` (user-settable): None -> _MANTA_OPT_TOPN; 0 -> ALL structures (optimise
    the complete ranked manifold, slowest/best); N>0 -> top-N; N<0 -> none."""
    if not isomers:
        return isomers
    if topn is None:
        _n = _MANTA_OPT_TOPN
    elif int(topn) == 0:
        _n = len(isomers)                  # 0 = ALL (optimise everything)
    elif int(topn) < 0:
        return isomers                     # negative = none
    else:
        _n = int(topn)
    import concurrent.futures as _cf
    try:
        from delfin.manta import _gfnff_rank as _gff
    except Exception:
        return isomers
    if not _gff.available():
        return isomers
    head = list(isomers[:_n])
    tail = list(isomers[_n:])

    def _opt_one(item):
        xyz, _na, label = item
        try:
            if str(spin) == "auto":
                # auto-spin: scan multiplicity -> GFN2 ground state (parity-correct)
                r = _gff.gfnff_optimize_autospin(xyz, charge=int(charge), method=method)
            else:
                # fixed multiplicity chosen by the user: uhf = multiplicity - 1
                _uhf = max(0, int(spin) - 1)
                r = _gff.gfnff_optimize(xyz, charge=int(charge), uhf=_uhf, method=method)
        except Exception:
            r = None
        if r and r[0]:
            opt_xyz, e = r
            na = len([ln for ln in opt_xyz.splitlines() if ln.strip()])
            _m = method.upper()
            tag = (" [%s-opt %.1f kcal]" % (_m, e)) if e is not None else " [%s-opt]" % _m
            return ((opt_xyz, na, (label or "isomer") + tag), e)
        return (item, None)

    workers = max(1, min(_MANTA_OPT_WORKERS, len(head)))
    try:
        with _cf.ThreadPoolExecutor(max_workers=workers) as ex:
            opted = list(ex.map(_opt_one, head))
    except Exception:
        return isomers
    # optimized structures first, sorted by GFN2-opt energy (failed/None last)
    opted.sort(key=lambda t: (t[1] is None, t[1] if t[1] is not None else 0.0))
    return [it for (it, _e) in opted] + tail

class Editor:
    """One structure editor: what to place, and what the tab still calls.

    Everything is reachable under the name it had while this lived in the
    Submit tab, so a tab that holds one reads the same as it did.
    """

    def __init__(self, parts):
        self.__dict__.update(parts)

    @property
    def exported(self):
        """The widgets a tab hands out, under the names they have here."""
        keep = ('mol_output', 'mol_status', 'mol_status_fs', 'manta_button',
                'manta_settings_row', 'convert_smiles_button',
                'convert_smiles_quick_button', 'convert_smiles_uff_button',
                'isomer_nav_row', 'isomer_label', 'isomer_prev_btn',
                'isomer_next_btn', 'xyz_copy_btn', 'xyz_copy_status')
        return {name: value for name, value in self.__dict__.items()
                if name.startswith(('submit_', 'manta_')) or name in keep}


def build(ctx, *, state, coords_widget, viewer_height, schedule_ui_update,
          update_view, get_smiles_charge, set_buttons_disabled=None,
          offer_structures=None, read_input=None, write_input=None,
          list_structures=None, show_output=None):
    """Make one structure editor over the coordinates a tab keeps.

    *state* is the tab's own dictionary -- the editor keeps its history, its
    held internal coordinates and its bond edits in it, so a tab that reloads
    a structure can clear them the way it always has. *coords_widget* is where
    an edited structure is written; changing it is what makes the tab redraw.
    The callables are what the editor cannot know: when it is safe to touch
    the interface, how the tab shows a structure again, what charge the SMILES
    it came from implies, which of the tab's own buttons a long conversion has
    to hold shut, and -- for a tab that keeps more than one structure -- where
    a set of them goes.

    *offer_structures* is handed every structure a conversion produced, as
    ``[(xyz, atom_count, label), ...]``. A tab that has somewhere to put them
    all takes them and says so by returning true; the ORCA Builder writes them
    as named blocks, and its own stepper walks them. A tab that says nothing
    gets the editor's isomer stepper instead, which is what the Submit tab has
    always had.
    """
    if set_buttons_disabled is None:
        def set_buttons_disabled(*_args, **_kwargs):
            # A tab with no submit buttons of its own has nothing to hold shut
            # while a conversion runs.
            return None
    if read_input is None:
        # Where a SMILES is typed. In the Submit tab that is the same box the
        # editor writes a structure into; in the ORCA Builder it is the box
        # holding the named blocks, which is not the box the editor writes to
        # at all -- so a conversion has to be told where to look.
        def read_input():
            return coords_widget.value
    if show_output is None:
        # Where a picture or a line of text goes when the editor has something
        # to show instead of a structure -- the loader that runs while MANTA
        # builds, "please enter coordinates", a refusal. The editor's own
        # output widget by default; a tab that renders the viewer itself hands
        # over its own, or the animation runs where nobody can see it.
        def show_output(items):
            mol_output.outputs = tuple(items)
    if list_structures is None:
        # What "all" means. The Submit tab holds a set of isomers; the ORCA
        # Builder holds named blocks, and each of those is a frame too.
        def list_structures():
            return state.get('isomers') or []
    if write_input is None:
        # And where a drawing is put, which is the same box: a drawing comes
        # back as a SMILES, and the buttons that turn one into coordinates
        # read from there.
        def write_input(text):
            coords_widget.value = text
    # The handful of keys the editor reads without asking first. A tab that
    # has never shown a structure has not written them, and a second host has
    # not written any of them; seeding them here is what lets the editor stand
    # on a dictionary it has not been given a tour of.
    for _key, _value in (('manip_bootstrap_done', False), ('manip_inflight', False),
                         ('perceived', None), ('poly_applied', None),
                         ('poly_metal', None), ('gfn_generation', 0),
                         ('gfn_scanned_uhf', None), ('smiles_task_id', 0),
                         ('smiles_busy', False), ('isomers', []),
                         ('isomer_index', 0),
                         ('converted_xyz_cache', {'smiles': None, 'xyz': None}),
                         ('current_xyz_for_copy', {'content': None}),
                         # The tab's own busy flag, if it has one; a conversion
                         # asks after both before it says the page is idle.
                         ('batch_preview_busy', False),
                         # What is running in the background right now, as
                         # {token: what it is}. Non-empty means the spinner
                         # belongs on screen, whatever any status line says.
                         ('busy_jobs', {})):
        state.setdefault(_key, _value)

    # What this session did, kept in memory so that a defect can be walked
    # through again.  It records into a bounded ring and touches no disk until
    # the user presses send, and `record` is the whole of what the hot paths
    # below pay for it -- see :mod:`delfin.dashboard.editor_journal` for what
    # is kept, what is dropped first, and where a report lands.
    #
    # A host that builds two editors gets two journals, which is right: a bug
    # report is about the structure the user was looking at, and the other
    # tab's session is somebody else's evidence.
    journal = _journal_mod.Journal()
    state['editor_journal'] = journal
    record = journal.record

    # The status line lies ON the picture, along its bottom edge, rather than
    # in a row above it.
    #
    # Above it, every message of a different length moved the structure: a run
    # reports several times a second and the reports are not the same size --
    # one row while it follows the hand, two when it says what it reached --
    # so the atom the user was aiming stepped up and down under the cursor.
    # Giving that row a fixed height stopped the movement and spent the height
    # of two rows, empty most of the time, to do it.
    #
    # Anchored to the bottom of the viewer, it costs no layout at all: it grows
    # upwards into the picture when there is more to say and shrinks again,
    # and nothing outside it can move. The tab puts it and the viewer in a box
    # carrying `delfin-structure-viewer-stack`; the rules are in the shared
    # sheet, with the rest of what a molecule panel looks like.
    def _status_layout():
        return widgets.Layout(width='auto', margin='0')

    mol_status = widgets.HTML(value='', layout=_status_layout())
    mol_status.add_class('delfin-structure-status-over')
    # A second one, kept and written but no longer a member of anything.  It
    # existed because fullscreen relocated the status line by hand and
    # ipywidgets knows nothing about such a move, so the line borrowed for the
    # big view did not come back to the small one.  The line lives inside the
    # stack with the picture now, and the stack is what fullscreen takes -- it
    # goes and comes back with the thing it is about, so there is nothing left
    # for a stand-in to do.  A tab that places it still gets the same text.
    mol_status_fs = widgets.HTML(value='', layout=_status_layout())
    mol_status_fs.add_class('delfin-structure-status-over')
    # It lives in the ordinary layout so fullscreen can pick it up from there,
    # but it must not be seen next to the line it is a copy of -- the message
    # would simply be printed twice.  Hidden here, shown in the overlay.
    mol_status_fs.layout.display = 'none'
    mol_output = widgets.Output(layout=widgets.Layout(
        border='2px solid #1976d2', width='100%', height=f'{viewer_height}px',
        overflow='hidden', box_sizing='border-box',
    ))
    mol_output.add_class('submit-mol-output')

    submit_scope_id = f'submit-scope-{abs(id(coords_widget))}'

    submit_fullscreen_btn = widgets.Button(
        description='', icon='expand',
        tooltip='Toggle fullscreen (Esc to exit)',
        layout=widgets.Layout(width='40px', height='30px'),
        disabled=True,
    )
    # The one fullscreen there is.  The editor used to carry a button for a
    # machinery of its own and the Builder had to take that class off again and
    # put the shared one on, which is two overlays' worth of code answering one
    # click.  The button is the shared one wherever the editor is built.
    submit_fullscreen_btn.add_class('delfin-structure-fullscreen-btn')

    submit_select_btn = widgets.ToggleButton(
        value=False, description='Select', icon='crosshairs',
        button_style='',
        tooltip=(
            'Click atoms to pick or unpick them. Hold Shift and drag for a '
            'rectangle; add Ctrl to keep the previous selection.'
        ),
        layout=widgets.Layout(width='96px', height='30px'),
        disabled=True,
    )
    submit_manip_btn = widgets.ToggleButton(
        value=False, description='Manipulate', icon='arrows',
        button_style='',
        tooltip=(
            'Grab any atom and drag it; grabbing a selected atom moves the '
            'whole selection. Drag empty space to turn the view. Right-click '
            'an atom to set the pivot, right-drag to rotate the selection '
            'about it.'
        ),
        layout=widgets.Layout(width='112px', height='30px'),
        disabled=True,
    )
    from .molecule_builder import DRAW_ELEMENTS

    submit_draw_btn = widgets.ToggleButton(
        value=False, description='Draw', icon='pencil',
        button_style='',
        tooltip=(
            'Click empty space to place an atom of the chosen element, drag '
            'from an atom to grow a new one bonded to it, drag onto another '
            'atom to bond them at the chosen order, tap an atom to change its '
            'element, and press Delete to remove the selection. Free valences '
            'are filled with hydrogen. The right button still turns the view.'
        ),
        layout=widgets.Layout(width='88px', height='30px'),
        disabled=True,
    )
    # Every element, with the ones a molecule editor reaches for first at the
    # top.  A native select types ahead by label and matches in the order the
    # options stand in, so pressing P lands on phosphorus rather than on
    # palladium, and I on iodine rather than on indium -- which is why the
    # common ones come first rather than the list being sorted by number.
    # Typing two letters quickly still finds Pd or In wherever they are.
    # Each one-letter symbol before its two-letter relatives, or pressing B
    # lands on bromine and P on palladium.
    _COMMON_ELEMENTS = [
        'C', 'H', 'N', 'O', 'S', 'P', 'F', 'B', 'I',
        'Cl', 'Br', 'Si',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ru', 'Rh', 'Pd', 'Ag', 'Ir', 'Pt', 'Au',
    ]
    _ALL_ELEMENTS = [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
        'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
        'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
        'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
        'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
        'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
        'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
        'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
    ]
    _ELEMENT_CHOICES = _COMMON_ELEMENTS + [
        e for e in _ALL_ELEMENTS if e not in set(_COMMON_ELEMENTS)]
    submit_element_dd = widgets.Dropdown(
        options=_ELEMENT_CHOICES,
        value='C',
        tooltip=('What a click draws. Every element is here; the common ones '
                 'are at the top, and typing the symbol jumps to it.'),
        layout=widgets.Layout(width='72px', display='none'),
    )
    submit_adjust_h_btn = widgets.ToggleButton(
        value=True, description='Adjust H', icon='tint', button_style='info',
        tooltip=(
            'Fill or trim the hydrogens on the atoms an edit touches, so a '
            'carbon that has just lost a bond gets one back. Switch it off to '
            'leave them exactly as they are -- which is what a radical, an '
            'open coordination site, or a fragment about to be joined to '
            'something else needs.'
        ),
        layout=widgets.Layout(width='96px', height='30px'),
    )
    submit_manip_clear_btn = widgets.Button(
        description='Clear', button_style='warning',
        tooltip='Clear selection & pivot',
        layout=widgets.Layout(width='66px', height='30px'),
        disabled=True,
    )
    # When a structure has been dragged out of the picture -- or a relaxation
    # has carried it out -- this puts it back: the camera on the centre of
    # mass, and everything in view. It changes nothing about the molecule.
    submit_centre_btn = widgets.Button(
        description='Centre', button_style='', icon='crosshairs',
        tooltip=('Put the system back in the middle of the view. '
                 'Nothing about the structure changes.'),
        layout=widgets.Layout(width='90px', height='30px'),
        disabled=True,
    )
    #: Atom numbers, and how big they are drawn -- the same pair the ORCA
    #: Builder has, from the same code: numbering belongs to a viewer, so
    #: both tabs take it from the shared editor part.
    submit_labels_btn = widgets.ToggleButton(
        value=False, description='#',
        tooltip=('Number every atom, in the order the coordinates are in. '
                 'The numbers sit on the atom centres, and one behind another '
                 'atom is hidden, so a crowded structure stays readable.'),
        layout=widgets.Layout(width='46px', height='30px'),
        disabled=True,
    )
    submit_label_size = widgets.BoundedIntText(
        value=LABEL_PX_DEFAULT,
        min=LABEL_PX_MIN,
        max=LABEL_PX_MAX,
        step=1,
        tooltip=('How tall the labels are, in pixels. Type one or step it; '
                 'they resize as you go.'),
        layout=widgets.Layout(width='62px', height='30px', display='none'),
    )
    #: What the labels say, which is not a fourth button.
    #:
    #: Every GFN and PM answer this editor has ever made computed the partial
    #: charges and threw them away with its scratch directory -- xtb writes
    #: them to a file on every single point, every optimisation cycle set and
    #: every scan point, and MOPAC writes them into the AUX file the
    #: trajectory is already read out of. Showing them costs no calculation at
    #: all; see :func:`delfin.dashboard.gfn_optimize.read_charges`.
    #:
    #: So the question is only where they go, and the answer is: on the
    #: control that already draws things on atoms. The toolbar is the crowded
    #: row and the standing rule about it is that less is more, so this is not
    #: a switch of its own -- it appears beside the size box, when the labels
    #: are on, and goes away with them. The cost of having it at all is a
    #: 96 px box that is invisible until somebody has already said they want
    #: labels.
    #:
    #: The charge entry is only in the list under a method that computes one.
    #: A browser force field has no charges to show, and offering the word
    #: there would be the interface promising something it cannot do -- the
    #: visible set of controls is meant to be the answer to "what can I do
    #: now".
    submit_label_what = widgets.Dropdown(
        options=[('number', 'number')],
        value='number',
        tooltip=(
            'What the labels say. "number" is the atom\'s place in the '
            'coordinates. "charge" is the partial charge the last answer '
            'computed -- it costs nothing, because every GFN and PM run '
            'writes the charges out whether anybody reads them or not, and '
            'they follow the drag as fast as the answers arrive. They are the '
            'method\'s own definition of a charge and not a measured '
            'quantity: the same methane carbon is -0.15 under GFN2, -0.13 '
            'under GFN1, -0.36 under g-xTB and -0.09 under GFN-FF, so they '
            'are read against each other and not across methods.'
        ),
        layout=widgets.Layout(width='84px', display='none'),
    )
    submit_manip_undo_btn = widgets.Button(
        description='Undo', button_style='info', icon='undo',
        tooltip='Undo last move/rotate (Ctrl-Z)',
        layout=widgets.Layout(width='84px', height='30px'),
        disabled=True,
    )
    # Beside Undo, which is the one place anybody looks for it.  A second
    # button is a cost, and this is the pair where the convention is strong
    # enough to pay it: everywhere else Redo sits next to Undo, and a Redo
    # that exists but has no control cannot be found at all.  It is greyed
    # until there is a way forward, so it never claims to offer one.
    submit_manip_redo_btn = widgets.Button(
        description='Redo', button_style='info', icon='repeat',
        tooltip='Put back what Undo took away (Ctrl-Shift-Z, or Ctrl-Y)',
        layout=widgets.Layout(width='84px', height='30px'),
        disabled=True,
    )
    submit_relax_btn = widgets.ToggleButton(
        value=False, description='Dynamik Opt', icon='magic',
        button_style='',
        tooltip=(
            'Keep the force field running: the structure settles continuously, '
            'and an atom you grab drags the relaxing molecule along with it. '
            'Use the strength control to make it gentler. Toggle off to stop.'
        ),
        layout=widgets.Layout(width='132px', height='30px'),
        disabled=True,
    )
    submit_optimize_btn = widgets.ToggleButton(
        value=False,
        description='Optimize', button_style='success', icon='compress',
        tooltip=(
            'Minimise the frame on screen. A switch: press again to stop, and '
            'what is stopped is what you were looking at. Undo restores the '
            'geometry from before the run in one step.'
        ),
        layout=widgets.Layout(width='112px', height='30px'),
        disabled=True,
    )
    submit_optimize_all_btn = widgets.ToggleButton(
        value=False,
        description='all', button_style='success',
        tooltip=(
            'Minimise every frame that is loaded -- all isomers or batch '
            'entries, not only the one on screen. They run one after another, '
            'never side by side: a login node is shared, and a set of isomers '
            'would otherwise start one xtb per frame at once.'
        ),
        layout=widgets.Layout(width='52px', height='30px'),
        disabled=True,
    )
    submit_ff_dd = widgets.Dropdown(
        # GFN1 is in this list because it is in every other list: the climb
        # takes it, the saddle search takes it, the solvent table was measured
        # for it one solvent at a time, and three of this editor's own
        # refusals name it as one of the methods to choose -- while the box
        # they send the user to did not have it. A control that is missing an
        # option the refusals point at is the interface saying two things at
        # once, and the visible set of controls is meant to be the answer to
        # "what can I do now".
        #
        # In the ladder the four xtb methods really form, which is also the
        # order of what they cost: GFN-FF, GFN1, GFN2, g-xTB.
        options=[('UFF', 'uff'), ('MMFF94', 'mmff94'),
                 ('GFN-FF', 'gfnff'), ('GFN1-xTB', 'gfn1'),
                 ('GFN2-xTB', 'gfn2'), ('g-xTB', 'gxtb'),
                 ('PM6-D3H4', 'pm6d3h4'), ('PM6', 'pm6'), ('PM7', 'pm7')],
        value='uff',
        tooltip=(
            'What Optimise minimises with. UFF and MMFF94 run in the browser '
            'and also drive the live relaxation while you drag. GFN-FF, '
            'GFN1-xTB, GFN2-xTB and g-xTB run xtb on the server, and they '
            'know about the metal where UFF guesses. GFN1 is the older '
            'Hamiltonian and is kept because a structure GFN2 cannot converge '
            'sometimes converges under it; GFN2 is the one to reach for '
            'first. g-xTB approximates wB97M-V/def2-TZVPPD '
            'and needs a build of its own; Install g-xTB fetches it. '
            'PM6-D3H4, PM6 and PM7 run MOPAC on the server. Measured on '
            'twelve small organics, PM6 draws bonds closer to the literature '
            'than GFN2 (5.0 against 11.3 mA on average), because GFN2 pulls '
            'multiple and polar bonds short; on four dimers PM6-D3H4 is the '
            'best of all of them and plain PM6 the worst -- it has no '
            'dispersion, and the water dimer comes apart.'
        ),
        layout=widgets.Layout(width='128px'),
        disabled=True,
    )
    # xtb needs both, and there is no honest default for a metal complex: the
    # wrong number of unpaired electrons gives a confident wrong answer rather
    # than an error.  Shown only when a GFN method is chosen.
    submit_gfn_charge = widgets.IntText(
        value=0, description='q', step=1,
        style={'description_width': '12px'},
        layout=widgets.Layout(width='72px', display='none'),
    )
    # A multiplicity is M = 2S+1, so the smallest one there is is 1.  The box
    # took 0 and -3 as readily as 2, and neither reached xtb as itself: the
    # conversion to unpaired electrons floors at zero, so both were quietly
    # run as a singlet.  That is the confident wrong answer the charge box
    # above is warned about, arrived at from a number the user could see was
    # not what they typed.  The box refuses them now instead.
    submit_gfn_mult = widgets.BoundedIntText(
        value=1, min=1, max=20, description='M', step=1,
        style={'description_width': '14px'},
        layout=widgets.Layout(width='72px', display='none'),
    )
    #: The thermal budget: what this structure can do at this temperature.
    #:
    #: Switched on, it takes the energy of the structure as it stands as the
    #: anchor and reports every later geometry against it.  What that buys is
    #: a number that says whether a distortion is chemistry or nonsense: a
    #: benzene ring bond at 1.55 A is 9.8 kcal/mol and reachable at room
    #: temperature, at 1.75 A it is 34.7 and reachable at none.
    #:
    #: The energy costs nothing.  Every follow step already relaxes the whole
    #: structure except the atom the hand is holding and hands its energy
    #: back; nothing read it.  That is exactly the right quantity, too -- a
    #: point on the constrained path rather than the energy of a geometry
    #: nothing has relaxed around.
    submit_thermal_btn = widgets.ToggleButton(
        value=False, description='Thermal', icon='thermometer-half',
        tooltip=(
            'Measure every change against what this structure can actually do '
            'at the temperature below. The energy is read from the relaxation '
            'that already runs while you drag, so it costs nothing. It goes '
            'with the pulling hand: under a placing hand what is kept is not '
            'exactly what was priced, so it is not offered there. '
            + _THERMAL_QUANTITY_SHORT
        ),
        layout=widgets.Layout(width='96px', height='30px', display='none'),
        disabled=True,
    )
    submit_temperature = widgets.BoundedFloatText(
        value=298.15, min=1.0, max=5000.0, step=25.0, description='T/K',
        style={'description_width': '28px'},
        layout=widgets.Layout(width='104px', display='none'),
        disabled=True,
    )
    #: Whether switching on relaxes first.  A budget is measured from
    #: somewhere, and a hand-built structure is often well above its own
    #: minimum -- anchoring there gives a ceiling that is generous by however
    #: much strain was already in it.  Relaxing first anchors on the minimum
    #: and makes the number mean what it says; leaving it off keeps the
    #: geometry exactly as the user built it, which is sometimes the point.
    submit_thermal_relax = widgets.Checkbox(
        value=True, description='relax first', indent=False,
        tooltip=('Relax once when switching on, so the budget is measured '
                 'from a minimum. Off anchors on the structure as it stands, '
                 'strain and all.'),
        layout=widgets.Layout(width='104px', display='none'),
        disabled=True,
    )
    #: Re-anchor here.  Crossing a barrier lands in a new minimum, and that is
    #: an intermediate: a molecule that has reacted thermalises and starts
    #: again, so the next elementary step is measured from there.  Pressed by
    #: hand rather than done automatically, because only the user knows
    #: whether what they are looking at is an intermediate or a way station.
    submit_thermal_anchor_btn = widgets.Button(
        description='Set here', icon='map-pin', button_style='',
        tooltip=('Measure from this structure from now on -- the next step of '
                 'a mechanism starts from the intermediate it reached.'),
        layout=widgets.Layout(width='104px', height='30px', display='none'),
        disabled=True,
    )
    # The lines between the atoms are 3Dmol's own bond list, worked out once
    # when the model was built.  Moving atoms does not touch it, so a bond
    # pulled apart goes on being drawn -- the picture keeps the connectivity
    # the structure had rather than the one it has.  Off by default: it costs a
    # rebuild per frame, and in a crowded coordination sphere the perception is
    # at its limit and the lines flicker.
    submit_dyn_bonds_btn = widgets.ToggleButton(
        value=False, description='Dyn. bonds', icon='link',
        tooltip=(
            'Let the lines between the atoms follow the distances while the '
            'structure moves, instead of keeping the bonds it was drawn with. '
            'Only the picture changes -- what the calculation holds together '
            'is a separate question, and Bond and Unbond are how that is said.'
        ),
        layout=widgets.Layout(width='104px', height='30px'),
    )
    submit_gfn_solvent = widgets.Dropdown(
        options=[(label, name) for name, label in _gfn.SOLVENTS.items()],
        value='',
        tooltip=(
            'Optimise with an implicit solvent around the structure. A '
            'geometry optimised in the gas phase and one optimised in water '
            'are different answers; the status line says which you got. Which '
            'solvents are offered depends on the model beside this: ALPB '
            'covers all 25, GBSA fewer, and how many fewer depends on the '
            'method.'
        ),
        layout=widgets.Layout(width='140px', display='none'),
    )
    #: The continuum itself.  Its options are rebuilt whenever the method
    #: changes, because what is available is a property of the method: GFN-FF
    #: has ALPB and GBSA, GFN1 and GFN2 add ddCOSMO, and the PM methods have
    #: MOPAC's COSMO and nothing else.  A dropdown that offers a model the
    #: chosen method cannot run is a dropdown that can only produce a refusal.
    submit_gfn_solv_model = widgets.Dropdown(
        options=[('ALPB', 'alpb')], value='alpb',
        tooltip=(
            'Which continuum stands in for the solvent. ALPB costs nothing '
            'and covers every solvent; ddCOSMO is five to six times the price '
            'and is what COSMO-RS wants; GBSA is the older model, kept '
            'because published numbers were computed with it. PM methods get '
            'MOPAC\'s COSMO, given the dielectric constant of the same '
            'solvent, so the two engines are asked about the same liquid.'
        ),
        layout=widgets.Layout(width='108px', display='none'),
    )
    submit_gfn_autospin = widgets.Checkbox(
        value=False, description='auto M', indent=False,
        tooltip=(
            'Try the multiplicities the electron count allows and keep the '
            'one that comes out lowest. For an open-shell metal a fixed guess '
            'gives a confidently wrong energy and, through it, a wrong '
            'geometry -- but it costs three runs instead of one.'
        ),
        layout=widgets.Layout(width='86px', display='none'),
    )
    # Shown only when a GFN method is chosen and xtb cannot be found.  Two
    # presses, not one: the install fetches a conda environment of a few
    # hundred megabytes, and on a cluster the right answer is often "load the
    # module instead" -- so what it would run is shown before it runs.
    submit_xtb_install_btn = widgets.Button(
        description='Install xtb', icon='download',
        tooltip='Fetch xtb with DELFIN\'s own installer. It will say what it '
                'runs before it runs it.',
        layout=widgets.Layout(width='118px', height='30px', display='none'),
    )
    submit_xtb_confirm_btn = widgets.Button(
        description='Yes, install', button_style='warning',
        layout=widgets.Layout(width='108px', height='30px', display='none'),
    )
    submit_xtb_cancel_btn = widgets.Button(
        description='Cancel',
        layout=widgets.Layout(width='72px', height='30px', display='none'),
    )
    submit_internal_label = widgets.HTML(
        value=(
            '<span class="submit-internal-label" '
            'style="color:#888;font-size:0.9em;white-space:nowrap;">'
            'pick 2-4 atoms</span>'
        ),
        layout=widgets.Layout(margin='0 0 0 4px'),
    )
    submit_internal_value = widgets.FloatText(
        value=0.0, step=0.01,
        layout=widgets.Layout(width='92px', height='30px'),
        disabled=True,
    )
    submit_internal_value.add_class('submit-internal-value')
    submit_internal_btn = widgets.ToggleButton(
        value=False, description='Set', button_style='primary',
        tooltip=(
            'Turn the value by hand and watch it: while Set is on, the box '
            'drives the selection live -- the arrow keys step a bond by '
            '0.01 A, an angle by 0.1 and a dihedral by 0.5 degrees, and the '
            'fragment on the far side of the coordinate follows. Two atoms are a bond '
            'length, three an angle, four a dihedral. Hold is the other '
            'question: it keeps a value at its target while the field runs, '
            'with pull or fix.'
        ),
        layout=widgets.Layout(width='58px', height='30px'),
        disabled=True,
    )
    submit_manip_status = widgets.HTML(
        value='<span class="submit-manip-status" style="color:#888;font-size:0.9em;">— viewer empty —</span>',
        # Takes a share of the row when there is room -- fullscreen keeps the
        # toolbar on one line -- and wraps to its own line when there is not,
        # which is what happens inside the tab on a laptop.
        layout=widgets.Layout(
            flex='1 1 260px', min_width='0', overflow_x='hidden',
        ),
    )
    submit_manip_sync = widgets.Textarea(value='', layout=widgets.Layout(display='none'))
    submit_manip_sync.add_class('submit-manip-sync')
    # Coordinates from the kernel, one frame at a time.  Not through run_js:
    # that writes into a single Output and clears it first, so twenty scripts a
    # second overwrite each other before the page has rendered them -- the
    # relaxation ran, the last structure appeared, and nothing in between did.
    # A widget value is ordered, survives a background thread, and cannot be
    # clobbered by the next one.
    #: Where the hand may no longer go, polled by the page.
    #:
    #: Not run_js: that clears its output before displaying, so a second call
    #: wipes the first, and the wall is written about ten times a second while
    #: a drag is running.  It was sent that way at first and simply never
    #: arrived -- the budget said the structure was fifty-five kcal/mol past
    #: what it can spend and the ring came apart anyway.  A widget the page
    #: reads on its own clock cannot be overwritten before it is seen, which
    #: is why the trajectory travels this way too.
    #:
    #: A field of its own rather than sharing the frame field, because that is
    #: one slot as well: a wall written between two frames would be gone by
    #: the time the page looked.
    submit_gfn_wall = widgets.Text(value='', layout=widgets.Layout(display='none'))
    submit_gfn_wall.add_class('submit-gfn-wall')
    submit_gfn_frame = widgets.Textarea(value='', layout=widgets.Layout(display='none'))
    submit_gfn_frame.add_class('submit-gfn-frame')

    submit_strength_slider = widgets.IntSlider(
        value=20, min=1, max=200, step=1,
        description='Strength', continuous_update=False,
        readout=True, readout_format='d',
        style={'description_width': '58px'},
        layout=widgets.Layout(width='200px'),
        disabled=True,
    )
    #: Which kind of hand this is.
    #:
    #: Moving is the older behaviour and it is not a worse one: the atom goes
    #: exactly where the cursor puts it and the rest of the structure settles
    #: around it.  That is what placing something *is*, and it is the right
    #: tool for building -- putting a ligand where you want it, closing a ring,
    #: setting up a geometry to start from.  A force cannot do that, because
    #: the whole point of a force is that the chemistry gets a say.
    #:
    #: Pulling is the newer one, and it is the right tool for asking: drag,
    #: and how far the atom gets is the answer.  They are two questions, so
    #: they are two settings rather than one behaviour that replaced another.
    submit_hand_dd = widgets.Dropdown(
        options=[('pull with a force', 'pull'), ('move the atom', 'move')],
        value='pull',
        tooltip=('Pull: the atom follows as far as the chemistry allows, and '
                 'how far that is is the answer. Move: the atom goes where '
                 'the cursor puts it and the rest settles around it, which is '
                 'what building a structure wants.'),
        layout=widgets.Layout(width='158px'),
        disabled=True,
    )
    #: How hard the hand pulls, as a share of a bond.
    #:
    #: Dragging used to set the atom's coordinates: the hand won absolutely,
    #: which is why a molecule could be pulled into any shape at all.  A pull
    #: is a force instead -- the atom follows the mouse along the flattest way
    #: out of where it stands, and how far it gets is the field's answer, not
    #: the cursor's.
    #:
    #: The number is a share of what a bond holds, and both engines are set
    #: against that same yardstick even though they apply it very differently.
    #: In the browser it is a share of a C-H stretch, 662 kcal/mol/A^2, and a
    #: tenth of it drives torsions and angles while stretching a bond by less
    #: than a tenth of an angstrom.  Under xtb it is a share of 110 kcal/mol/A,
    #: which is where a C-C (112), a C-H (98) and a C-O (120) all give way.
    #:
    #: 1.0 is a hand as strong as the bond, which can break it.  0 is the old
    #: rigid hand, for placing an atom exactly where it is wanted.
    #:
    #: It opens at 0.4, which is about what room temperature allows: it turns
    #: a molecule into its own conformers and it does not break a bond.
    #:
    #: The temperature does not touch this, with the budget on or off.  It
    #: once did -- the hand was derived from the ceiling over a reach -- and
    #: that needs a length no temperature supplies; sized as a distance it was
    #: too weak to turn a torsion, so a molecule could not be put into its own
    #: conformers at exactly the temperature that certainly allows it.  What
    #: the temperature limits is the energy of what is *kept*, and the wall is
    #: what enforces that.  See :func:`_pull_force`.
    submit_pull_slider = widgets.FloatSlider(
        value=0.4, min=0.0, max=3.0, step=0.05,
        description='Pull', continuous_update=False,
        readout=True, readout_format='.2f',
        tooltip=('The hardest dragging may pull, as a share of a bond. 0.4 is '
                 'about what room temperature allows: it turns a molecule '
                 'into its own conformers and cannot break a bond. 1.0 is as '
                 'strong as the bond; 3.0 takes anything apart; 0 places the '
                 'atom outright. What is actually pulling is set by the '
                 'mouse: the hand is a spring between the cursor and the '
                 'atom, so the further you drag ahead of where the structure '
                 'has got to, the harder it pulls -- up to this. The line '
                 'under the viewer says both numbers while you drag. The '
                 'thermal budget does not change either of them: it limits '
                 'the energy of what you are left with.'),
        style={'description_width': '58px'},
        layout=widgets.Layout(width='200px'),
        disabled=True,
    )
    #: How far the structure moves for how far the mouse moves.
    #:
    #: One is the cursor and the atom staying together: let go and the atom is
    #: where the pointer is. Below one the hand travels further than the
    #: structure, which is what a crowded region wants; above one it travels
    #: less far, for reaching across a large system without running out of
    #: desk. It scales dragging and turning a selection only -- where a *new*
    #: atom is placed stays one to one, because an atom that appeared
    #: somewhere other than under the cursor would be a different kind of
    #: wrong.
    submit_sens_slider = widgets.FloatSlider(
        value=1.0, min=0.2, max=3.0, step=0.1,
        description='Mouse', continuous_update=False,
        readout=True, readout_format='.1f',
        tooltip=('How far the structure moves for how far the mouse moves. '
                 '1.0 keeps the atom under the cursor; lower is finer, higher '
                 'reaches further. Placing a new atom is always 1.0.'),
        style={'description_width': '58px'},
        layout=widgets.Layout(width='200px'),
        disabled=True,
    )
    #: How fast the trajectory is played back, in frames a second.
    #:
    #: It is a real choice now that the whole path reaches the page rather than
    #: a sample of it.  Slow is the useful end: the calculation runs on ahead
    #: while the picture walks the path, and grabbing an atom takes the frame
    #: on screen and throws the computed remainder away.  Fast keeps the
    #: picture at the calculation, and then there is little to throw away.
    #:
    #: The old pacing sped the playback up whenever the queue grew, on the
    #: grounds that a picture trailing its calculation was a fault.  That is
    #: the behaviour a slow setting exists to ask for, so the setting wins and
    #: the backlog rule is gone.
    #: Whole frames a second, from one up to live.
    #:
    #: The top is not a speed so much as "live", and what it means depends on
    #: the system: the frames only exist as fast as xtb makes them.  On
    #: something small the queue is full and sixty a second is a quick replay
    #: of a path already finished; on a hundred and fifty atoms the picture
    #: drains the queue and then waits for the next frame, which is the
    #: calculation watched as it happens.  So the top is the useful end for a
    #: large structure and the middle for a small one, and dialling down from
    #: live is the way to read it.
    #:
    #: Whole steps, because a tenth of a frame a second is a distinction
    #: nobody makes; one a second is already a second of looking at each
    #: geometry.
    #:
    #: It was briefly at the top of the range, on the theory that a viewer
    #: taking longer to arrive was the complaint.  It was not -- that was a
    #: flag left raised over a write that changed nothing, and it is fixed
    #: where it was.  At sixty a second a twenty-three-frame optimisation is
    #: over in 0.4 s, which is not a trajectory anyone can watch; it is the
    #: single jerk the playback was built to replace.  Twelve draws that same
    #: path in 1.9 s and a long one in proportion, and the slider is there for
    #: both ends.
    submit_play_speed = widgets.IntSlider(
        value=12, min=1, max=60, step=1,
        description='Speed', continuous_update=False,
        readout=True, readout_format='d',
        tooltip=('How many frames of the optimisation are drawn a second. '
                 'Slow lets the calculation run ahead of the picture: what '
                 'you see is where you are, and grabbing an atom there keeps '
                 'that frame and drops what was computed past it. At the top '
                 'the picture keeps up with the calculation, so on a large '
                 'structure that is the trajectory as it is computed; on a '
                 'small one it is a fast replay of a path already finished.'),
        style={'description_width': '48px'},
        layout=widgets.Layout(width='178px', display='none'),
        disabled=True,
    )
    submit_pick_sync = widgets.Text(value='', layout=widgets.Layout(display='none'))
    submit_pick_sync.add_class('submit-pick-sync')
    # Keyboard shortcuts for things Python owns. Unbond is not a picture edit:
    # it changes the topology the force field is built from, so the browser
    # cannot carry it out alone and has to ask through here.
    submit_cmd_sync = widgets.Text(value='', layout=widgets.Layout(display='none'))
    submit_cmd_sync.add_class('submit-cmd-sync')
    submit_poly_dd = widgets.Dropdown(
        options=[('— polyhedron —', '')], value='',
        layout=widgets.Layout(width='190px', display='none'),
        disabled=True,
    )
    submit_poly_turn_btn = widgets.Button(
        description='Turn', icon='refresh', button_style='',
        tooltip=(
            'Step to the next distinct arrangement on this polyhedron: which '
            'ligands take the axial or apical positions. Only shown where the '
            'vertices are not all alike -- an octahedron has nothing to turn, '
            'a trigonal bipyramid has ten arrangements.'
        ),
        layout=widgets.Layout(width='78px', height='30px', display='none'),
        disabled=True,
    )
    submit_hyb_dd = widgets.Dropdown(
        options=[('— hybridisation —', '')], value='',
        layout=widgets.Layout(width='190px', display='none'),
        disabled=True,
    )
    submit_hyb_auto_btn = widgets.Button(
        description='Types', icon='magic', button_style='',
        tooltip=(
            'Read each carbon\'s hybridisation off how many partners it is '
            'bonded to: four is tetrahedral, three trigonal planar, two '
            'linear. Stronger than perception, which goes through bond '
            'orders and misses a double bond it cannot see in the geometry. '
            'Applies to the selection, or to every carbon when nothing is '
            'selected.'
        ),
        layout=widgets.Layout(width='84px', height='30px'),
        disabled=True,
    )
    #: On to begin with, and only ever the browser's.  Leaving the strain of a
    #: drag in the structure is the surprising answer, not the useful one: a
    #: hand-placed atom measured 176 kcal/mol above a settled one on a real
    #: complex, and that is what would go to the queue.
    #:
    #: It was briefly off, because it was the one thing moving a structure
    #: while everything the user had touched said off -- but that was the
    #: server methods, where the switch is not on the toolbar to be found and
    #: what it started was the whole minimisation rather than a tidy-up.  That
    #: is fixed where it was broken: the toolbar takes it away under a server
    #: method and switches it off with it, and the gate on the release path no
    #: longer reads this widget at all.  So under the browser's field, where it
    #: means what its name says and costs nothing, it can go back to on.
    submit_settle_btn = widgets.ToggleButton(
        value=True, description='Settle', icon='level-down',
        button_style='info',
        tooltip=(
            'When you let go of an atom, let the structure relax around its '
            'new position instead of keeping the strain of the drag. Switch '
            'off to leave atoms exactly where you put them.'
        ),
        layout=widgets.Layout(width='92px', height='30px'),
        disabled=True,
    )
    #: What letting go of an atom leads to, while the molecule is following
    #: the hand.  It used to depend on something else entirely: a drag that
    #: interrupted a running Optimise brought that run back and the structure
    #: went down to a minimum, while the same drag with no run behind it left
    #: the structure where the hand put it.  One gesture, two outcomes, decided
    #: by a switch on the other side of the toolbar -- which reads as it
    #: working only sometimes.  It is asked for here instead.
    submit_auto_btn = widgets.ToggleButton(
        value=True, description='Auto', icon='angle-double-down',
        button_style='info',
        tooltip=(
            'While Dynamik Opt is on: when you let go of an atom, carry on -- '
            'down to a minimum, or up to a transition state while Climb to TS '
            'is down. Switch off to have it stop where your drag left it, so '
            'you can move something else first.'
        ),
        layout=widgets.Layout(width='78px', height='30px'),
        disabled=True,
    )
    submit_swap_btn = widgets.Button(
        description='Swap', button_style='', icon='exchange',
        tooltip=(
            'Exchange the two selected ligands on the polyhedron: they are '
            'pulled onto each other\'s vertex instead of back to their own.'
        ),
        layout=widgets.Layout(width='78px', height='30px', display='none'),
        disabled=True,
    )
    submit_bond_btn = widgets.Button(
        description='Bond', icon='link', button_style='',
        tooltip=(
            'Draw a bond between the two selected atoms. Distance-based '
            'perception is unreliable in a crowded coordination sphere, and '
            'the coordination number and the force field both follow from '
            'these bonds.'
        ),
        layout=widgets.Layout(width='74px', height='30px'),
        disabled=True,
    )
    submit_unbond_btn = widgets.Button(
        description='Unbond', icon='unlink', button_style='',
        tooltip='Remove the bond between the two selected atoms.',
        layout=widgets.Layout(width='90px', height='30px'),
        disabled=True,
    )
    submit_hold_mode = widgets.Dropdown(
        options=[('pull', 'pull'), ('fix', 'fix')], value='pull',
        layout=widgets.Layout(width='78px'),
        disabled=True,
    )
    submit_hold_btn = widgets.Button(
        description='Hold', button_style='warning', icon='thumb-tack',
        tooltip=(
            'Hold the value the selection describes while the field runs, '
            'instead of only setting it once. Held values appear in the list '
            'beside this button and can be dropped again there.'
        ),
        layout=widgets.Layout(width='72px', height='30px'),
        disabled=True,
    )
    #: The scan: the same selection as Hold, walked instead of held.
    #:
    #: Dragging answers "can this happen at this temperature" and answers it
    #: honestly, but the number beside it depends on how the hand came in --
    #: the same Diels-Alder measured +16.2 one way and +40.3 another.  A scan
    #: walks a coordinate the user has *named*, in equal steps, relaxing each
    #: point from the one before, which is the calculation that produced every
    #: reference figure in this file.  Grid spacing is not decoration either:
    #: sampled coarsely that approach reads +4.8 kcal/mol and finely +12.8, so
    #: the step count belongs to whoever is asking.
    #: Drag without letting the molecule become a different one.
    #:
    #: A third thing the drag can be asked to respect, beside the budget and
    #: the running field.  The budget says what a temperature can pay for; this
    #: says the bonds are the ones they were, whatever it costs -- which is
    #: what "move this group over there" means when the answer is not supposed
    #: to be a reaction.
    #:
    #: xtb cannot be told to keep a topology, so this is a wall like the
    #: budget's: the step runs, the graph is read off what came back, and a
    #: step that made or broke a bond is replaced by the last one that did
    #: not.  Under GFN-FF it is already true -- that method reads its bonding
    #: once and keeps it -- and under the browser's field it is true by
    #: construction, since the terms are assigned once.  Where it says
    #: something is GFN2 and its relatives, which decide the bonding afresh
    #: from the electrons at every step.
    submit_topology_btn = widgets.ToggleButton(
        value=False, description='Keep bonds', icon='link',
        tooltip=('Drag without making or breaking any bond. A step that '
                 'changes the bonding is taken back, so the molecule stays '
                 'the molecule it was.'),
        layout=widgets.Layout(width='118px', height='30px'),
        disabled=True,
    )
    submit_scan_btn = widgets.Button(
        description='Scan', button_style='info', icon='line-chart',
        tooltip=(
            'Walk the value the selection describes from where it is to where '
            'you say, relaxing everything else at every step. Several '
            'coordinates can be armed and are then walked together, which is '
            'what a concerted reaction needs.'
        ),
        layout=widgets.Layout(width='74px', height='30px'),
        disabled=True,
    )
    #: Where the walk goes, said in the words the coordinate is in.
    #:
    #: "closer" and "further" next to a bare 0 and a bare 20 said nothing
    #: about what any of the three were: the direction read as a setting with
    #: no subject, and the two numbers as numbers with no units.  The options
    #: are rewritten for the kind of coordinate that is picked -- two atoms
    #: are closer or further apart, three or four are narrower or wider -- and
    #: the numbers carry their own labels.
    #:
    #: The end of the walk is the third answer here and not a checkbox beside
    #: it.  Normally there is no end: a scan walks to the next minimum, so
    #: where it stops is the chemistry.  A number is wanted for two cases --
    #: following a figure from the literature, and fixing the ratio of two
    #: coordinates in a coupled scan -- and a value already says which way the
    #: walk goes, so "further apart, to 2.40" was the same fact twice with
    #: nothing checking that the two halves agreed.  One question, three
    #: answers, and the field for the number appears under the third.
    #:
    #: And two more for a pair of atoms, which are the same question answered
    #: as chemistry rather than as geometry: **form this bond** and **break
    #: this one**.  A direction says where the coordinate goes; a verb says
    #: what is supposed to have happened when the walk is over -- and that is
    #: the difference, because it is what the walk stops on.  Armed one each
    #: on two pairs and walked together, they are the instruction most
    #: reactions actually are: make this one while breaking that one.  A
    #: Diels-Alder is two forms, an SN2 is a form and a break.
    #:
    #: The shape is pyGSM's, whose driving-coordinate file is a list of atom
    #: picks and verbs -- ``ADD 4 12``, ``BREAK 1 11`` -- and SCINE's NT2 does
    #: the same thing with a force and stops on bond order.  Nothing new is
    #: being invented here: a form is an armed leg driven inwards, a break is
    #: one driven outwards, the existing ramp walks them together, and the
    #: only addition is what it stops on.  See :func:`_carried_out` for the
    #: stopping rule and why it is not a bond order.
    #:
    #: **When two instructions fight, the ramp is what settles it**, and that
    #: was measured rather than reasoned about.  Asked to form C1-C11 -- half
    #: a Diels-Alder, which is easy -- while breaking C1-C2, which is a
    #: butadiene double bond and one of the strongest things in the molecule,
    #: with both pulling on the same carbon under one shared force constant:
    #:
    #:     step  6   19.6 kcal/mol/A   the form is granted, -63.7 kcal/mol
    #:     step 12   57.3              C1-C2 has reached 1.61 A, -62.4
    #:     step 16  117.3              1.87 A, -35.1
    #:     step 17  140.3              2.13 A, -3.8, and both now hold
    #:
    #: So they do not deadlock: the cheap half is granted at a low force and
    #: the expensive half holds out until the force passes what that bond
    #: holds against -- :data:`gfn_optimize.A_BOND_HOLDS`, 110 kcal/mol/A, and
    #: it went at 140.  The price of the fight is the sixty kcal/mol the path
    #: climbs back through, which is on the profile and answers to the
    #: temperature like any other rise.  And when the ramp ends with a verb
    #: still unsatisfied, that is an answer too -- it ends above twice what a
    #: bond holds, so it is a statement about this method and this structure
    #: rather than a setting to turn up.
    #:
    #: They are offered for a pair and not for an angle or a torsion, because
    #: a bond is between two atoms and there is no third one to make or break.
    submit_scan_way = widgets.Dropdown(
        options=[('closer together', 'in'), ('further apart', 'out'),
                 ('form this bond', 'form'), ('break this bond', 'break'),
                 ('to a value you give', 'to')],
        value='in',
        tooltip=(
            'Where to walk. A scan stops at the next minimum, so where it '
            'ends is the chemistry rather than a number -- which way it goes '
            'is the one thing that cannot be read off the selection. Form or '
            'break says what is meant to have happened instead, and the walk '
            'stops when it has: arm one on each of two pairs to make one bond '
            'while breaking another. Give a value when the end is the point.'
        ),
        layout=widgets.Layout(width='178px', display='none'),
        disabled=True,
    )
    submit_scan_to = widgets.FloatText(
        value=0.0, step=0.1, description='to',
        style={'description_width': '18px'},
        layout=widgets.Layout(width='96px', display='none'),
        disabled=True,
    )
    submit_scan_steps = widgets.BoundedIntText(
        value=20, min=2, max=400, step=1, description='steps',
        style={'description_width': '38px'},
        layout=widgets.Layout(width='104px', display='none'),
        disabled=True,
    )
    #: Walk the whole thing, rather than stopping at the next minimum.
    #:
    #: "all the way" said nothing about what it was all the way *to*.  What it
    #: is for is a whole curve: a torsion turned right round has three minima
    #: and three barriers, and stopping at the first of them is exactly wrong
    #: when the profile is the point.  Off is the right default -- past the
    #: next minimum a reaction scan is pushing into a structure rather than
    #: following anything.
    submit_scan_whole = widgets.ToggleButton(
        value=False, description='Whole profile', icon='chart-line',
        tooltip=(
            'Keep walking past the next minimum, for the whole curve -- a '
            'torsion turned right round has three minima and three barriers. '
            'Off, the scan stops once it is over a barrier and has settled '
            'again, because past that it is pushing into a structure rather '
            'than following a reaction.'
        ),
        layout=widgets.Layout(width='146px', height='30px', display='none'),
        disabled=True,
    )
    #: Walk the coordinate, or push it with a force.
    #:
    #: Walking dictates the path: the value is told what to be at every step
    #: and everything else relaxes around it.  That is right when the
    #: coordinate really is the reaction -- a torsion, a ring flip -- and
    #: wrong when it is a guess, because the structure is then driven through
    #: whatever the guess implies.  It is also how atoms come to collapse into
    #: each other: a value is a value, and the scan meets it.
    #:
    #: A push is an artificial force between the atoms instead, ramped up
    #: until something happens -- Maeda and Morokuma's AFIR, which is how a
    #: reaction path is found rather than assumed.  Where the structure ends
    #: up at each force is the structure's answer, so it can slide sideways,
    #: take the flatter way round, and refuse to be squeezed: the repulsion
    #: wins as soon as the force cannot pay for it.  Measured on butadiene and
    #: ethylene 3.13 A apart under GFN2: at 10 kcal/mol/A they close to 2.47
    #: (+5.2 kcal/mol), and at 20 the Diels-Alder goes over the top and lands
    #: in cyclohexene at 1.53 (-64.1).  Nothing said which bonds to form.
    submit_scan_how = widgets.Dropdown(
        options=[('push with a force', 'push'), ('walk the value', 'hold')],
        value='push',
        tooltip=(
            'Walk: the coordinate is told what to be at every step. Push: an '
            'artificial force between the atoms, ramped up until the reaction '
            'happens -- the structure decides where it goes, so the path is '
            'found rather than assumed.'
        ),
        layout=widgets.Layout(width='150px', display='none'),
        disabled=True,
    )
    #: Price the walk with an electronic energy, or with a free one.
    #:
    #: E throughout is the default and is a fair approximation: an entropy
    #: term that is roughly constant along a path largely cancels in a
    #: difference.  Where it does not -- a reaction that ties two molecules
    #: into one, a rotor that stops turning at the top -- G is the number a
    #: chemist actually wants, and the ceiling it is compared against is a
    #: free energy anyway.
    #:
    #: Not at every point.  A free energy is a Hessian, and a Hessian is not
    #: free: measured under GFN2, 0.57 s against 0.29 for sixteen atoms and
    #: 3.72 against 0.76 for twenty-four, which over twenty points is a scan
    #: that takes minutes instead of seconds.  And an RRHO free energy only
    #: means something at a stationary point.  So it is asked for at the three
    #: places where it is both affordable and meaningful: where the walk
    #: started, the highest point it crossed, and the minimum it came to.
    submit_scan_energy = widgets.Dropdown(
        options=[('price with E', 'E'), ('price with G', 'G')], value='E',
        tooltip=('Electronic energies throughout, or a free energy at the '
                 'start, the top and the end -- three Hessians, at the '
                 'temperature above.'),
        layout=widgets.Layout(width='136px', display='none'),
        disabled=True,
    )
    #: Walk the same coordinate back again, and say whether the two agree.
    #:
    #: A driven scan holds one coordinate and relaxes everything else, and
    #: nothing makes the rest follow continuously.  Where it does not, the
    #: profile depends on which way the walk went, and neither direction is
    #: the path.  Jonsson, Mills and Jacobsen said it in 1998: "the path
    #: generated may be discontinuous and the procedure may depend on the
    #: direction of the drag ... some atomic coordinates may slip near the
    #: saddle point region and the saddle point configuration will then be
    #: missed."  Bofill and Quapp give the condition it holds under: no
    #: turning point and no valley-ridge inflection.
    #:
    #: There is no way to know which case a given scan is in except to walk it
    #: back.  Measured here on butadiene and ethylene under GFN2, one forming
    #: C-C driven from 3.40 A to 1.60 and back, 0.1 A at a time, everything
    #: relaxed at every point:
    #:
    #:     forward, apparent barrier            +7.3 kcal/mol at 2.20 A
    #:     backward, apparent barrier          +11.7 kcal/mol at 2.90 A
    #:     ORCA's converged saddle              +6.8
    #:     largest gap at the same coordinate   23.8 kcal/mol
    #:
    #: Both maxima carry exactly one imaginary frequency, so that test says
    #: nothing here.  A user reading the barrier off either leg alone is out
    #: by more than half, in opposite directions.
    #:
    #: And the null result is the reason this is on by default.  The same
    #: measurement over ten other scans -- torsions of an alkane, an alcohol,
    #: a diol and an amide, a C-C stretch, a C-C-C angle, a hydrogen bond, an
    #: SN2 and a ring opening -- gave gaps under 0.1 kcal/mol.  Those scans
    #: are worth quoting and the editor could not say so; now it can, and the
    #: price is one more leg of the walk that was just watched.  Off is one
    #: press for a large system where that price is felt.
    #:
    #: Walk mode only.  A push is a ramp of forces and not a grid of values,
    #: so there is no "the same coordinate, backwards" to walk; it finds its
    #: own crossing and prices it with :func:`_across` instead.
    submit_scan_back = widgets.ToggleButton(
        value=True, description='Walk it back', icon='rotate-left',
        button_style='info',
        tooltip=(
            'After the scan, walk the same coordinate back from where it '
            'ended and compare. A driven scan is only a path where nothing '
            'slips sideways, and the two legs disagreeing is how that shows. '
            'Costs a second leg of the same walk.'
        ),
        layout=widgets.Layout(width='138px', height='30px', display='none'),
        disabled=True,
    )
    #: Find the way between the two ends a scan has just produced.
    #:
    #: A scan drives a coordinate somebody chose; xtb's path finder is given
    #: two structures and finds its own way between them, with metadynamics
    #: pushing away from where it has been and pulling towards the product.
    #: So it answers the question the scan can only approach -- and it needs
    #: the scan first, because the scan is what makes a product to aim at.
    #:
    #: Measured on the butadiene and ethylene the scan had just walked: 3.6 s,
    #: forward barrier 5.755 kcal/mol, backward 69.586, reaction energy
    #: -63.831, and an estimated transition state with the two forming bonds
    #: at 2.524 and 2.520 A.  The scan of the same reaction put its highest
    #: point at +6.3 and 2.36.  Two methods that share no machinery, agreeing.
    #: Two structures, marked one at a time.
    #:
    #: The path finder needs a start and an end and cannot invent either.  A
    #: scan leaves both, which is where this began -- but a great many
    #: questions arrive as two structures the user already has, and a
    #: cis/trans isomerisation is the plainest of them: build cis, mark it,
    #: build trans, press.  Nothing about the reaction has to be guessed at,
    #: which is exactly the case a scan is worst at.
    #:
    #: Measured on 2-butene: a scan of the C-C=C-C torsion cannot do it at all
    #: -- one dihedral does not pin four substituents, so the constraint is
    #: met by pyramidalising the carbons rather than by twisting the bond, and
    #: the walk reports +95 kcal/mol at 150 degrees and a "trans" 64 above
    #: cis.  Given the two structures instead, the path finder answers in four
    #: seconds: 51.4 kcal/mol forward, trans 1.53 below cis, and its estimated
    #: transition state at 87.5 degrees, which is where a twisted alkene is.
    submit_path_from_btn = widgets.Button(
        description='Mark this end', icon='map-pin',
        tooltip=('Mark the structure on screen as one end of a path. Then '
                 'load or build the other one, and the press beside this can '
                 'start from the pair.'),
        layout=widgets.Layout(width='128px', height='30px'),
        disabled=True,
    )
    #: Where the search starts from.
    #:
    #: There were three buttons here and they were not three tools.  ``To the
    #: saddle`` climbed whatever was in the box, ``Path to saddle`` walked
    #: between two ends and climbed the estimate, and ``Climb to TS`` did the
    #: first of those by hand -- three of the six answers to two questions,
    #: with the difference between them living only in the tooltips.  This is
    #: the first question, and it is asked once.
    #:
    #: A start can come from three places and two of them appear on their own:
    #: what is on screen, the end marked beside this, and the two ends a
    #: finished scan leaves.  ``Path to saddle`` had both of the last two and
    #: no way to say which it would use -- it silently preferred the marked
    #: pair -- so a scan walked after marking something was unreachable.
    #:
    #: Shown only when there is more than one place to choose between, because
    #: a list of one is not a choice.  Nothing is hidden by that: with nothing
    #: marked and no scan walked, "the two ends" names a pair that does not
    #: exist, and the press starts where the only start is.
    submit_saddle_from = widgets.Dropdown(
        options=[('what is on screen', 'here')], value='here',
        tooltip=('Where the search starts. What is on screen, the end you '
                 'marked, or the two ends the last scan left.'),
        layout=widgets.Layout(width='158px', display='none'),
        disabled=True,
    )
    #: And how it gets there, which is the other question.
    #:
    #: Three answers, and each of them is offered exactly while it can run:
    #:
    #: * Through ORCA.  ORCA's OptTS on xtb gradients: measured on a
    #:   sixteen-atom Diels-Alder transition state, ``! XTB2 OPTTS`` converged
    #:   in under seven seconds from an estimate with its forming bonds at
    #:   2.524 and 2.520 A to a symmetric saddle at 2.315 and 2.315 with one
    #:   imaginary mode.  It needs a method ORCA has a keyword for, which is
    #:   not the whole GFN family.
    #: * By hand.  :mod:`climb` -- P-RFO with eigenvector following on a
    #:   Bofill-updated Hessian -- walking on xtb gradients at about 10 ms a
    #:   step, which is slow enough to watch and to interrupt and fast enough
    #:   to drag in the middle of.  It needs a method it knows how to ask for
    #:   a gradient from; see :data:`delfin.dashboard.climb.CLIMB_METHODS`.
    #: * The path only.  The walk between the two ends and nothing after it:
    #:   3.6 s against 12 for the pair on that same case, and a complete
    #:   answer -- a barrier forward, a barrier back, and how near the walk
    #:   came to the structure it was aimed at.  It was ``Find the path``, and
    #:   it is not a second tool: it is this one stopping a step earlier.
    #:
    #: The last two are offered from two ends and not from the structure on
    #: screen, where there is no walk to stop after and where climbing by hand
    #: is ``Climb to TS`` -- a switch, over with the other switches, and one
    #: press either way.  A second control that did the same thing under
    #: another name is what this whole row was suffering from.
    #:
    #: Optimising is the ordinary case and is the default.  The path alone is
    #: kept on the list rather than produced only when a climb fails, because
    #: it is the answer on a machine with no ORCA and the fast answer on one
    #: with it -- and something a user would go looking for must be findable
    #: before they know they need it.
    #: * A nudged elastic band.  ORCA's ``! NEB-TS``: a chain of images
    #:   relaxed onto the way between the two ends at once, the highest of
    #:   them climbed to the saddle.  The arbiter, and second on the list
    #:   rather than first: measured on the same sixteen-atom Diels-Alder from
    #:   the same two ends, it reaches the same saddle as the chain above to
    #:   0.07 cm-1 and spends 203 gradients doing it.  What it is for is the
    #:   case where the cheap answer is not believed, or where the two ends
    #:   are far enough apart that a cheap interpolation cannot bridge them --
    #:   it is a different method and not a longer run of the same one.
    #:
    #:   How long that is depends on the machine more than on the method, and
    #:   the number this editor used to quote -- seven minutes -- was the
    #:   serial one.  Measured here: the same band, same hour, 272 s on one
    #:   process and 39.4 s on eight, because the images are independent
    #:   gradients and ORCA computes them together.  On a box with cores it is
    #:   a press like the others; on a small login node it is 203 gradients
    #:   however they are arranged, which is what the timeout is for.
    submit_saddle_how = widgets.Dropdown(
        options=[('through ORCA', 'orca')], value='orca',
        tooltip=('How it gets there. Through ORCA converges it; NEB-TS '
                 'relaxes a whole chain of images between the two ends and '
                 'climbs the highest, which is slower and is what to reach '
                 'for when the fast answer is not believed; by hand you can '
                 'watch it, steer it and drag in the middle of it; the path '
                 'only stops at the structure the walk estimates.'),
        layout=widgets.Layout(width='142px', display='none'),
        disabled=True,
    )
    #: The one press, and what it does is whatever the two boxes beside it say.
    #:
    #: Neither engine answers alone from two ends -- xtb walks between them and
    #: estimates the top of what it crossed and has no saddle optimiser at
    #: all; ORCA has one and nothing that produces an estimate to give it --
    #: and the pair is quick.  Measured on the sixteen-atom Diels-Alder, from
    #: the two ends a scan leaves: 3.6 s for the walk, 12 s for the pair, and
    #: it arrives at -393.5 cm-1 where a nudged elastic band on the same two
    #: ends reaches -393.6 for about twice the gradients.  The same saddle to
    #: a tenth of a wavenumber, which is why the band is the second entry in
    #: the box beside this and not the first.
    #:
    #: From what is on screen it is the interactive half of the question: a
    #: structure posed by hand into something that looks like a transition
    #: state is exactly the case for it.  One core, two gigabytes and three
    #: minutes, because this runs where the dashboard runs and on a cluster
    #: that is the login node.  A saddle search on a real basis set is a job,
    #: and the ORCA Builder is where jobs are submitted.
    #:
    #: Its name says which of the two it is about to do, because "To the
    #: saddle" over a walk that stops at an estimate would be a promise the
    #: press does not keep.  While it runs it says ``Stop``, and stopping
    #: keeps what was reached.
    submit_saddle_btn = widgets.Button(
        description='To the saddle', icon='mountain', button_style='warning',
        tooltip=('Optimise to the nearest transition state, from wherever the '
                 'box beside this says. Says whether what it reached is one.'),
        layout=widgets.Layout(width='140px', height='30px'),
        disabled=True,
    )
    #: What the structure on screen actually is, asked of a Hessian.
    #:
    #: This is a press of its own, and the case for one is that nothing else
    #: on the row means it. The whole machinery -- the Hessian, the count of
    #: modes going the wrong way, the sentences that name a minimum, a
    #: transition state or a saddle of higher order -- was already built and
    #: was reachable only from inside the saddle search, the path walk and the
    #: scan. So somebody who dragged a structure into a shape, which is what
    #: this editor is for, could not ask whether the shape is anything. An
    #: absence on the toolbar is a statement about what can be done, and that
    #: one was not true.
    #:
    #: It stands beside "To the saddle" because it is the same question read
    #: the other way round: that press goes looking for a stationary point,
    #: and this one asks whether there is already one here. They share the
    #: wording, from :func:`delfin.dashboard.saddle.verdict`, so the two can
    #: never disagree about what a structure is called.
    #:
    #: Affordable, which is the other half of the case. Measured here through
    #: the editor's own path -- which is to say with the thread count set the
    #: way :func:`~delfin.dashboard.gfn_optimize.interactive_cores` sets it --
    #: under GFN2 on a shared machine carrying other people's work:
    #:
    #:     5 atoms, methane                    0.4 s  (4 cores)
    #:     23 atoms, heptane                   1.1 s  (4 cores)
    #:     57 atoms, a manganese complex, +1  23.7 s  (8 cores)
    #:
    #: So it is a press and a short wait, never something to put in a drag.
    #: The thread count is the whole difference at the top of that list: the
    #: same 57-atom Hessian run without it, letting OpenMP take what it liked
    #: on a box at a load average of 800, took 14 minutes 29 seconds of wall
    #: clock for 13 hours of CPU. A Hessian is the one thing here that will
    #: happily buy nothing with every core on the machine.
    #:
    #: The thermochemistry rides along free. A Hessian is what a free energy
    #: costs; H, T*S and the zero-point energy are printed in the same block
    #: as G and were being read past, and the temperature box the answer is
    #: computed at is already on this row.
    submit_shape_btn = widgets.Button(
        description='What is it?', icon='question', button_style='',
        tooltip=(
            'Take a Hessian on the structure exactly as it stands and say '
            'what it is: a minimum, a transition state, a saddle of higher '
            'order -- or none of those, which is the answer for a structure '
            'that is still on a slope. Nothing is optimised and nothing is '
            'moved. The free energy, enthalpy, entropy and zero-point energy '
            'come with it, at the temperature in the box on this row. '
            'Seconds for a small structure and minutes for a large one, so '
            'it is a press rather than something that happens by itself.'
        ),
        layout=widgets.Layout(width='106px', height='30px', display='none'),
        disabled=True,
    )
    #: Which mode the two presses after this are about, when there is a choice.
    #:
    #: There usually is not.  A transition state has exactly one mode going the
    #: wrong way, so on the structure the press next door is for there is
    #: nothing to choose between and this is not on screen -- the same rule the
    #: two boxes before it follow, that a list of one is not a choice.
    #:
    #: It appears where a search converged onto a saddle of higher order, and
    #: there it is the difference between reading the verdict and acting on it:
    #: the verdict's own advice is to move the structure along the *second* of
    #: those modes and climb again, and until now there was no way to look at
    #: the second of anything.
    #:
    #: The frequencies are the ones the verdict already reported, which is at
    #: most the four softest -- :func:`saddle._last_modes` and
    #: :meth:`climb.Climb.verdict` both cut their list there.  A structure with
    #: more than four modes going the wrong way is a long way from anything
    #: this row is for, and the four that can be named are the four offered.
    submit_mode_dd = widgets.Dropdown(
        options=[('the imaginary mode', 0)], value=0,
        tooltip=('Which of the modes going the wrong way to draw and to '
                 'follow down.'),
        layout=widgets.Layout(width='134px', display='none'),
        disabled=True,
    )
    #: Draw the imaginary mode, because the imaginary mode is the reaction.
    #:
    #: A saddle search reports "one mode goes the wrong way, at -394 cm-1" and
    #: that number is the whole of what the editor could say about it.  The
    #: mode itself is the reaction coordinate: the atoms it moves are the atoms
    #: the reaction moves, so watching it is what tells a chemist which bonds
    #: are forming and which are breaking, without reading a single coordinate.
    #:
    #: It is a picture and it never becomes a structure.  The frames go down
    #: the same channel a trajectory goes down and are drawn by the same
    #: player, and the coordinate box is not touched: what is in the box is
    #: what is being worked on, and a geometry displaced along an eigenvector
    #: is not something anybody chose.  It begins and ends on the structure,
    #: and a hand landing on it puts the picture back before the drag can push
    #: anything to the kernel -- see the ``home`` frame in the player.
    #:
    #: It costs one xtb Hessian, once: 0.41 s at sixteen atoms, 4.9 at
    #: thirty-three, 16 at fifty, and it is kept for as long as the structure
    #: in the box is the one it was taken at.
    submit_mode_btn = widgets.Button(
        description='Show the mode', icon='film', button_style='',
        tooltip=('Draw the imaginary mode -- the reaction coordinate itself, '
                 'so you can see which bonds are forming and which are '
                 'breaking. Six swings and it stops on the structure it '
                 'started from. It is a picture: the coordinates never '
                 'change.'),
        layout=widgets.Layout(width='142px', height='30px', display='none'),
        disabled=True,
    )
    #: And where the mode leads, which is the other question about a saddle.
    #:
    #: One imaginary mode makes a structure *a* transition state.  Whether it
    #: is *the* one for the reaction that was meant is a different question,
    #: and the standard way to ask it is to push the structure a little way
    #: down the mode each way, relax, and see which two structures come back.
    #:
    #: The rigorous form of that is a mass-weighted steepest descent -- an
    #: intrinsic reaction coordinate -- and ORCA has one.  Both were measured
    #: on the sixteen-atom Diels-Alder saddle on this machine: this route is
    #: 1.0 s and lands on the two minima, ``! XTB2 IRC`` is 207 s and stops in
    #: the valley above them, and the two then agree exactly.  Two hundred
    #: seconds is longer than :func:`saddle.seconds_for` allows the saddle
    #: search itself, on the smallest case anybody runs, on what is somebody's
    #: login node -- so the press is the cheap one and the IRC is a job.  See
    #: :func:`climb.follow_the_mode_down`.
    submit_ends_btn = widgets.Button(
        description='Follow it down', icon='code-fork', button_style='',
        tooltip=('Push the structure a little way down the imaginary mode '
                 'each way, relax both, and say which two structures this '
                 'saddle joins -- what they are and what they cost against '
                 'it. Seconds, and the box keeps the saddle.'),
        layout=widgets.Layout(width='146px', height='30px', display='none'),
        disabled=True,
    )
    #: The same climb, slowed down enough to put a hand in the middle of it.
    #:
    #: ORCA's OptTS is a press and cannot be anything else: measured, a
    #: three-step burst of it costs 2.7 to 3.1 s for sixteen atoms and 2.5 of
    #: those are ORCA starting up, so a saddle search made of repeated ORCA
    #: runs is two orders of magnitude away from a drag.  :mod:`climb` is the
    #: same method -- P-RFO with eigenvector following on a Bofill-updated
    #: Hessian -- walking on xtb gradients here, at 10 ms a step, and that is
    #: fast enough to watch and to interrupt.
    #:
    #: What the hand does to it is the point.  Grab an atom and the climb is
    #: interrupted exactly as a minimisation is -- it does not fight the hand
    #: and it does not walk on underneath it, because a saddle of a restrained
    #: surface is not a saddle of the real one.  Measured on the Diels-Alder,
    #: climbing with both forming bonds held: at 2.20 A it ends 0.53 A from
    #: the saddle with *two* imaginary modes and a true gradient 138 times the
    #: convergence threshold, and at 2.60 A it converges in five steps onto a
    #: point with no imaginary mode at all.  Unrestrained, the same climb takes
    #: 11 steps and lands 0.006 A from where ORCA lands.
    #:
    #: When the mouse is let go the climb starts again from the structure that
    #: was made -- and takes the direction of the drag as the mode to follow,
    #: which is what actually helps.  From the same dragged geometry, climbing
    #: the mode the hand pointed at reaches the Diels-Alder saddle in 39 steps
    #: at 2.315 A; climbing the lowest mode instead walks back down to the
    #: van-der-Waals complex 0.43 A away with no imaginary mode at all.
    #:
    #: So this is a *mode*, not a run: it says which way the release walks and
    #: it stays down across the walks it starts, the way Dynamik Opt and Auto
    #: do.  Everything else about a drag is theirs -- the follow under the
    #: hand, the interrupt at the grab, the restart after the release, the run
    #: number, the frame channel -- and there is one path through all of it
    #: with the optimiser as the only difference.  There used to be two, and
    #: every defect this button has had was the second one disagreeing with
    #: the first.
    submit_climb_btn = widgets.ToggleButton(
        value=False, description='Climb to TS', icon='hand-rock',
        tooltip=('Which way the optimiser walks: up to a transition state '
                 'rather than down to a minimum, a step at a time on xtb '
                 'gradients. Stays down like Dynamik Opt -- drag an atom to '
                 'point it at the reaction you mean, and Auto carries on '
                 'towards the saddle when you let go. Tap the atom you are '
                 'aiming at first and the search knows which contact to check '
                 'itself against, and tries two more ways when the first '
                 'misses.'),
        layout=widgets.Layout(width='140px', height='30px'),
        disabled=True,
    )
    submit_scan_run_btn = widgets.Button(
        description='Run scan', button_style='success', icon='play',
        tooltip='Walk every armed coordinate together, and say what it costs.',
        layout=widgets.Layout(width='104px', height='30px', display='none'),
        disabled=True,
    )
    #: The second opinion on a finished walk, and the moment it arrives is the
    #: whole of why it is a button.
    #:
    #: It appears when a scan has left a profile and goes away again when the
    #: next one starts, in the row that already comes and goes with the scan
    #: -- so its being there is the editor saying that something can be done
    #: now which could not be done a minute ago, the way the two ends were
    #: before f1be8954 folded them into a box.  There is nothing to set: the
    #: geometries are the ones just walked, and the method is the only one in
    #: the list better than what walked them, so a box beside it would have
    #: one entry and be a question with one answer.
    #:
    #: Not a setting chosen before the scan, which is where it would cost no
    #: pixels at all.  The question it answers -- "is that barrier really that
    #: small" -- is one nobody has until the barrier is on screen, and by then
    #: a setting is minutes of walking away while this is seconds: measured,
    #: 0.40 s a point against 0.9 s a point for the walk itself, on the same
    #: sixteen-atom complex.
    submit_scan_price_btn = widgets.Button(
        description='Price with g-xTB', button_style='info',
        icon='balance-scale',
        tooltip=(
            'Take the geometries this scan just walked and work out their '
            'energies again with g-xTB, which is much more accurate than the '
            'GFN methods and nearly as quick -- seconds for a whole profile. '
            'GFN2 systematically underestimates barriers where a bond is '
            'being broken, by about 10 kcal/mol on published test sets. What '
            'comes back is a screen, not a barrier to quote: the structures '
            'are still the ones GFN2 found.'
        ),
        layout=widgets.Layout(width='168px', height='30px', display='none'),
        disabled=True,
    )
    submit_scan_dd = widgets.Dropdown(
        options=[('nothing armed', '')], value='',
        layout=widgets.Layout(width='230px', display='none'),
        disabled=True,
    )
    submit_scan_del = widgets.Button(
        description='', icon='times', button_style='danger',
        tooltip='Drop this leg of the scan',
        layout=widgets.Layout(width='40px', height='30px', display='none'),
        disabled=True,
    )
    submit_constraint_dd = widgets.Dropdown(
        options=[('no constraints', '')], value='',
        layout=widgets.Layout(width='210px', display='none'),
        disabled=True,
    )
    submit_reset_btn = widgets.Button(
        description='Reset', icon='undo', button_style='danger',
        tooltip=(
            'Back to the structure as it was loaded, and drop everything set '
            'on it since: held values, polyhedron, hand-made bonds, '
            'hybridisation overrides and the edit history.'
        ),
        layout=widgets.Layout(width='84px', height='30px'),
        disabled=True,
    )
    submit_constraint_del = widgets.Button(
        description='', icon='times', button_style='danger',
        tooltip='Drop the selected constraint',
        layout=widgets.Layout(width='40px', height='30px', display='none'),
        disabled=True,
    )
    #: The internal-coordinate row: a row inside a row, and so a row that has
    #: to wrap and to give way like the one it sits in.
    #:
    #: It used to say ``row nowrap`` and ``flex: 0 0 auto``, and that is the
    #: whole of why the last buttons hung off the right of the screen.  The
    #: toolbar around it does wrap; a wrapping container breaks *between* its
    #: items and never inside one, and this group is one item.  Armed with a
    #: scan it is nineteen controls and about 1900 px of content, in a row 620
    #: px wide at a 1280 px window, laid out on a single line that had nowhere
    #: to go.
    #:
    #: Measured in chromium before this changed, with a scan set up: at 1280
    #: px, ten controls past the right edge of the toolbar embedded and eight
    #: in the fullscreen overlay, and nine and seven of those not the topmost
    #: element at their own centre -- unreachable, and not reachable by
    #: scrolling either, because every container between here and the document
    #: holds its overflow.  The three path buttons that used to sit at the end
    #: of this row were three of them, which is exactly what was reported.
    #: They have since become one press, and it stands with the switches that
    #: say what a walk does rather than at the end of this row.
    #:
    #: ``0 1 auto`` written out rather than left to the default, because the
    #: difference from ``0 0 auto`` is the point: when the group is wider than
    #: the line it is alone on it shrinks to the line and wraps inside itself
    #: instead of painting past it.  ``min_width: 0`` because a flex item's
    #: automatic minimum is its own content, and its content is the widest
    #: control in it.
    #:
    #: In the state the editor opens in -- the label, the value, Set, Hold and
    #: the mode, and nothing else on the row -- the group fits on its line
    #: either way and nothing moves: measured, the toolbar is the same height
    #: to the pixel at 1920, 1536, 1280 and 1024 px.
    submit_internal_group = widgets.HBox(
        [submit_internal_label, submit_internal_value,
         submit_internal_btn, submit_hold_btn, submit_hold_mode,
         submit_scan_btn, submit_scan_way, submit_scan_to,
         submit_scan_steps,
         submit_scan_dd, submit_scan_del, submit_scan_whole,
         submit_scan_how, submit_scan_energy, submit_scan_back,
         submit_scan_run_btn, submit_scan_price_btn],
        layout=widgets.Layout(
            gap='6px', align_items='center', flex_flow='row wrap',
            flex='0 1 auto', min_width='0', overflow='visible',
        ),
    )

    #: The labels and their two settings, as one item of the toolbar.
    #:
    #: A row inside a row, the way the internal-coordinate group is, and for
    #: the same reason: a flexbox breaks *between* its items and never inside
    #: one, so the number of direct children is what decides where the toolbar
    #: can wrap.  What the labels say is a setting of the labels, so it goes
    #: where the labels are; making the three of them one item is what keeps
    #: that from costing the row a place to wrap between.  The toolbar carries
    #: the same number of items it did before, plus the one new press.
    #:
    #: Where the group sits is a separate fact and is written where it is
    #: placed, in the toolbar's own list.
    #:
    #: ``flex 0 1 auto`` and ``min_width 0`` written out for the reason the
    #: internal group writes them out: a nested row that cannot shrink takes
    #: its whole content past the edge of the toolbar however narrow the
    #: window is.
    submit_label_group = widgets.HBox(
        [submit_labels_btn, submit_label_what, submit_label_size],
        layout=widgets.Layout(
            gap='6px', align_items='center', flex_flow='row wrap',
            flex='0 1 auto', min_width='0', overflow='visible',
        ),
    )

    #: A line break for a wrapping toolbar.  Flexbox has no "break here", so
    #: the break is an element: nothing wide, taking a whole line, which
    #: pushes everything after it onto the next row.  It is inert in the
    #: ordinary view -- where the toolbar sits beside the rest of the tab and
    #: wraps where it must -- and only takes effect inside the overlay, where
    #: the row is as wide as the screen and would otherwise put the two
    #: Optimise buttons at the far end of a very long first line.
    submit_fs_row_break = widgets.Box(
        [], layout=widgets.Layout(display='none'))
    submit_fs_row_break.add_class('submit-fs-row-break')

    submit_manip_toolbar = widgets.HBox(
        [
            submit_fullscreen_btn,
            submit_select_btn, submit_manip_btn, submit_draw_btn,
            submit_element_dd, submit_adjust_h_btn,
            submit_manip_clear_btn, submit_centre_btn,
            submit_manip_undo_btn, submit_manip_redo_btn, submit_reset_btn,
            # After Reset rather than before Undo, which is three buttons to
            # the right of where the numbering used to sit and is a placement
            # rather than a preference. Measured in chromium with a scan armed
            # and every control forced visible -- the widest state the toolbar
            # test builds -- at 1024 px embedded: an item added anywhere before
            # the charge box pushes that box to the end of a wrapped line, and
            # there its input paints 24 px past the 72 px container it was
            # given, two of them past the toolbar itself. Nothing about this
            # group causes that and nothing here can fix it; what this
            # placement does is not provoke it. After Reset the toolbar is
            # clear at 1024, 1280, 1536 and 1920 px, embedded and in the
            # overlay.
            submit_label_group,
            submit_ff_dd, submit_gfn_charge, submit_gfn_mult,
            submit_gfn_autospin, submit_gfn_solvent, submit_gfn_solv_model,
            submit_thermal_btn, submit_temperature,
            submit_thermal_relax, submit_thermal_anchor_btn,
            submit_topology_btn,
            submit_xtb_install_btn, submit_xtb_confirm_btn,
            submit_xtb_cancel_btn,
            submit_strength_slider, submit_hand_dd, submit_pull_slider,
            submit_sens_slider, submit_play_speed,
            submit_fs_row_break,
            # Climb to TS stands with the other switches that say what a walk
            # does, because that is what it is: the same release path with the
            # optimiser walking uphill instead of down (e3442010). It used to
            # sit beside the saddle press, where it read as a third way of
            # pressing for a transition state rather than as the sibling of
            # Dynamik Opt and Auto that it is.
            submit_optimize_btn, submit_optimize_all_btn,
            submit_relax_btn, submit_auto_btn, submit_climb_btn,
            submit_settle_btn,
            # And the one press for a transition state next to them, with the
            # two boxes that say where it starts and how it gets there. It is
            # the third thing that makes this structure move -- down, up under
            # a hand, or onto a saddle -- and it belongs with the other two
            # rather than at the far end of the internal-coordinate row, where
            # its three predecessors had accumulated.
            #
            # Measured in chromium with a scan armed, against the same row
            # before the three became one: the toolbar is shorter at seven of
            # the eight widths and places it is drawn in and the same at the
            # eighth -- 306/272 at 1920 embedded, 406/372 at 1536, 472/436 at
            # 1280, 604/568 at 1024, and 204/170, 236/204, 236/236, 304/272 in
            # the fullscreen overlay. Put after the internal group instead it
            # is 34 px taller at 1280 in the overlay; put inside it, taller
            # again at three widths.
            submit_path_from_btn, submit_saddle_from,
            submit_saddle_how, submit_saddle_btn, submit_shape_btn,
            # And what there is to do with a saddle once the press has found
            # one, immediately after it. They are absent until it has: the
            # visible controls are the answer to "what can I do now", so a
            # control for a structure nobody has yet is the row saying
            # something untrue. Their arriving is how the editor says the two
            # questions about a transition state can now be asked.
            submit_mode_dd, submit_mode_btn, submit_ends_btn,
            submit_poly_dd, submit_poly_turn_btn,
            submit_hyb_dd, submit_hyb_auto_btn,
            submit_internal_group,
            submit_bond_btn, submit_unbond_btn, submit_dyn_bonds_btn,
            submit_swap_btn, submit_constraint_dd, submit_constraint_del,
            submit_pick_sync, submit_cmd_sync,
            submit_manip_status, submit_manip_sync, submit_gfn_frame,
            # On the page, or the page cannot read it.  It was created and
            # written to and never placed, so the leash lived in the kernel's
            # memory alone: the line said the budget was spent and the atom
            # went on moving, because nothing on the other side had ever seen
            # a wall.  A widget that is not in the tree is not in the DOM.
            submit_gfn_wall,
        ],
        layout=widgets.Layout(
            display='none', gap='6px', align_items='center',
            width='100%', flex_flow='row wrap',
            margin='0 0 6px 0', overflow='visible',
        ),
    )
    # The toolbar goes into the overlay with the structure it acts on, wherever
    # this editor is built. Each host used to say so for itself, which is the
    # same fact written down twice and a third host free to forget it.
    submit_manip_toolbar.add_class('delfin-structure-fs-member')
    submit_manip_toolbar.add_class('delfin-structure-fs-toolbar')
    # And a name of its own, for the rules that hold wherever it is -- in the
    # Submit tab, in the ORCA Builder, and in the body-level overlay, which is
    # outside both tabs and so outside the tab sheets that used to be the only
    # thing keeping a control inside its column. The shared sheet is where
    # those rules live; see the toolbar section of
    # :data:`~delfin.dashboard.molecule_viewer.STRUCTURE_VIEWER_FULLSCREEN_CSS`.
    submit_manip_toolbar.add_class('delfin-structure-toolbar')

    #: Reporting what just went wrong, in the shape the DELFIN Agent tab
    #: already has it: a line to write in, and send.
    #:
    #: One button on the row until it is asked for, because the row is the
    #: crowded one and the standing rule about it is that less is more. The
    #: line and Send appear under the beetle and go away again, so the cost of
    #: having this at all is 44 px next to controls that are 200.
    #:
    #: Where it goes is said on the control itself rather than in a manual.
    #: The local copy is the user's own structures on his own machine and he
    #: can only judge that if he can see the directory it went to; the copy
    #: that goes on to his cluster is the same directory over the same
    #: transfer the dashboard uses for everything else, and the control names
    #: that destination before he presses rather than after.
    #:
    #: And the tooltip says which way round the button works. It was already
    #: recording -- into memory, from the moment the editor was built -- and
    #: the first thing the user asked when he saw it was whether it should
    #: not always record, which is the interface's answer to read, not his.
    #: A button that could equally be a start switch has to say it is not
    #: one, in the words that distinguish the two: what is here already, and
    #: handed over.
    submit_bug_btn = widgets.ToggleButton(
        value=False, description='🐞', tooltip=(
            'This viewer keeps what you do as you work, in memory. Nothing is '
            'written or sent until you press Send, and Send hands over what '
            'has already happened -- it does not start a recording.'),
        layout=widgets.Layout(width='44px', flex='0 0 auto'),
    )
    submit_bug_note = widgets.Text(
        value='', placeholder='What went wrong? (one line)',
        layout=widgets.Layout(width='260px', display='none'),
    )
    submit_bug_send = widgets.Button(
        description='Send', button_style='warning',
        tooltip='Hand over what has already happened',
        layout=widgets.Layout(width='72px', display='none'),
    )
    submit_bug_where = widgets.HTML(
        value='', layout=widgets.Layout(display='none'))
    #: A row inside a row has to wrap too, or its far end leaves the screen:
    #: a flexbox breaks between its items and never inside one, so a nested
    #: row saying nowrap puts all four of these on one line however narrow
    #: the toolbar is. That is how the internal-coordinate row came to hang
    #: off the right of the screen, and it is caught here by the test that
    #: was written for it rather than found again in a browser.
    submit_bug_group = widgets.HBox(
        [submit_bug_btn, submit_bug_note, submit_bug_send, submit_bug_where],
        layout=widgets.Layout(gap='6px', align_items='center',
                              flex_flow='row wrap',
                              flex='0 1 auto', min_width='0',
                              overflow='visible'),
    )
    # Appended to the row rather than written into the list that builds it.
    # The control then has one place of its own that stays put while the
    # contents of that list are rearranged around it, which is what a row
    # holding forty controls gets rearranged for.
    submit_manip_toolbar.children = (
        tuple(submit_manip_toolbar.children) + (submit_bug_group,))

    def _under_home(path):
        """A path as the user thinks of it, with his home written as ``~``."""
        text = str(path)
        home = str(Path.home())
        return '~' + text[len(home):] if text.startswith(home) else text

    #: A press that would file a report with no sentence in it, counted.
    #: Refused once, allowed the second time -- see :func:`on_submit_bug_send`.
    _bug_wordless = [0]

    def on_submit_bug_toggle(change=None):
        """Open the line to write in, or put it away again.

        What opens with it is the amount being held, in words, and the two
        places a press would put it. That the viewer remembers a session was
        true before this line existed and invisible, which is why the user
        asked for a feature he already had: a bounded ring filling up in
        memory has no appearance at all, and the button beside it could as
        easily have been a start switch. A count is the smallest thing that
        makes the difference visible -- it is a number that was already
        growing while he worked, so it cannot be read as an invitation to
        begin.

        Counted when the control is opened rather than kept live. A live
        counter would write to a widget several times a second on the same
        paths the journal is careful to cost nothing on, and it would be
        doing it to say a thing that does not change: that the session is
        being kept. The number is the evidence, not the instrument.
        """
        if change is not None and change.get('name') != 'value':
            return
        open_now = bool(submit_bug_btn.value)
        for widget in (submit_bug_note, submit_bug_send, submit_bug_where):
            widget.layout.display = 'flex' if open_now else 'none'
        if not open_now:
            # Putting the control away forgets that a wordless press was
            # already refused once, so the rule is per opening rather than
            # per session: somebody who comes back to this an hour later is
            # asked for a sentence again rather than filing without one on a
            # press he does not remember making.
            _bug_wordless[0] = 0
            return
        held = _journal_mod.how_much_is_held(journal.summary())
        here = _under_home(_journal_mod.resolve_archive_dir())
        target = _journal_mod.remote_target()
        goes = f'→ <code>{html.escape(here)}</code>'
        if target is None:
            goes += ', and no further -- no transfer host is configured'
            submit_bug_send.tooltip = (
                'Write the report. It stays on this machine: no transfer '
                'host is configured.')
        else:
            goes += (', then <code>'
                     + html.escape(target['shown']) + '</code>')
            submit_bug_send.tooltip = (
                'Write the report and send it to ' + target['shown']
                + '. The copy here is kept whatever the transfer does.')
        submit_bug_where.value = (
            '<span style="font-size:11px;opacity:.75">'
            f'Kept as you work: <b>{html.escape(held)}</b>. '
            'Send hands over what is already here.<br>'
            f'{goes}</span>')

    def on_submit_bug_send(_button=None):
        """Write everything this session did into a report, send it, say both.

        The sentence and the sequence go out together. Neither is much use
        alone: the sentence says what to look for and the sequence is what
        puts it back on the screen, and it is the pair that turns "es zappelt"
        into a test somebody can fail.

        A wordless press is refused once and let through on the second, which
        is neither of the two things it could have been. Filing silently
        without a sentence hands a maintainer a megabyte of coordinates and no
        idea what he is looking for -- of the four defects this feature was
        built for, the one that is still open is the one nobody could
        describe. Refusing outright would be worse: a user who cannot put a
        name to what he saw is exactly the user whose report is worth the
        most, and he must not be the one who cannot file. So the first press
        says what the sentence is for, and the second files it with the report
        marked so the maintainer knows to ask.

        The transfer runs on a thread. A third of a megabyte over SSH is
        nothing when the host answers and is a wedged connect when it does
        not, and the button is on the toolbar of a live editor: the local copy
        and its path are on the status line before the network is touched at
        all, and what the transfer did arrives after.
        """
        note = (submit_bug_note.value or '').strip()
        if not note:
            _bug_wordless[0] += 1
            if _bug_wordless[0] == 1:
                record('note', v='send refused: nothing was typed')
                _set_mol_status(
                    'Say in one line what went wrong, then press Send again.',
                    'The sequence is already kept and can be replayed by '
                    'anyone; your sentence is the only part of a report '
                    'nobody else can reconstruct. Send again to file it '
                    'without one.')
                return
        else:
            _bug_wordless[0] = 0

        record('note', v='the bug button was pressed')
        touched = journal.did_anything()
        try:
            report_dir = _journal_mod.write_report(
                journal,
                description=note,
                widgets=journal_watching,
                tab=str(state.get('editor_host') or ''),
            )
        except Exception as trouble:
            _set_mol_status(
                'The report could not be written: ' + str(trouble),
                'Nothing was lost -- what you did is still remembered, so '
                'Send can be pressed again.')
            return
        told = journal.summary()
        where = _under_home(report_dir)
        record('note', v=f'report written to {report_dir}')
        # The drop is said out loud rather than left in the file. A report
        # whose beginning the ring dropped replays from part-way through, and
        # the one person who can say whether that matters is standing here.
        went = (f' The session outran the buffer, so the oldest '
                f'{told["dropped"]} things you did are not in it.'
                if told['dropped'] else '')
        # Two presses in one session file two reports and the second is the
        # wider one, because a send does not empty the ring. Said here as well
        # as in the file: the person deciding whether to press again is
        # standing at the button, not reading the report.
        earlier = ('the first one' if journal.reports_written == 2
                   else 'the earlier ones')
        again = (f' This is report {journal.reports_written} from this '
                 f'session; it holds everything {earlier} did and what has '
                 f'happened since.'
                 if journal.reports_written > 1 else '')
        # And a press from a viewer nobody has touched. It is still a report
        # -- the structure that was loaded and every control as it stands are
        # in it, which is the whole of a complaint about what the viewer
        # opened on -- but it has no gesture to replay, and the difference
        # matters enough to say before he waits for somebody to reproduce it.
        empty = ('' if touched else
                 ' This viewer had sent nothing yet, so what is in it is the '
                 'structure it opened on and the controls as they stand, and '
                 'there is no sequence to replay.')
        wordless = ('' if note else
                    ' It went without a sentence, so it says to ask you what '
                    'went wrong.')
        # The count is not offered in that case. It would be a count of the
        # editor talking to itself -- the box it filled on load, the note this
        # press just wrote -- read under a phrase that says it is a count of
        # what the user did, which is the one number he cannot check.
        did = (f'{_journal_mod.how_much_is_held(told)} are in {where}'
               if touched else f'nothing you did is in it. It is in {where}')

        target = _journal_mod.remote_target()
        if target is None:
            _set_mol_status(
                f'Reported. {did}',
                'No transfer host is configured, so it stays on this '
                'machine. It can be played back into a viewer to see the '
                'same thing happen again.' + went + again + empty + wordless)
        else:
            _set_mol_status(
                f'Reported. {did}',
                f'Sending it to {target["shown"]}...'
                + went + again + empty + wordless)
            _start_background(
                lambda: _bug_report_travels(report_dir, where, target),
                'Sending the bug report')

        submit_bug_note.value = ''
        submit_bug_btn.value = False
        _bug_wordless[0] = 0

    def _bug_report_travels(report_dir, where, target):
        """Hand the written report to the configured transfer, from a thread.

        Both outcomes are said and neither is said quietly. A silent success
        leaves the user unable to tell a report that reached the cluster from
        one that never left, and a silent failure is worse than not offering
        the transfer at all: he would believe a maintainer has it.

        Whatever happens, the sentence ends on the local copy. That is not
        reassurance, it is the instruction -- the path is what he sends in a
        mail, or what somebody reads over his shoulder, when the cluster is
        the thing that is down.
        """
        # The destination it was told to go to is the one the status line
        # named a moment ago, handed on rather than resolved again.
        ok, said = _journal_mod.send_report(report_dir, target=target)
        record('note', v=('sent to ' + said) if ok
               else ('could not send: ' + said))
        schedule_ui_update(
            _set_mol_status,
            ('Sent to ' + said) if ok
            else ('It could not be sent: ' + said),
            (f'The copy in {where} is kept either way, and it is the same '
             f'report -- what went is the directory that is here.' if ok else
             f'The report is in {where} and nothing was lost -- that path is '
             f'what to pass on if the transfer stays down.'))

    submit_bug_btn.observe(on_submit_bug_toggle, names='value')
    submit_bug_send.on_click(on_submit_bug_send)

    def _set_mol_status(*lines, spinner=False):
        """Say what the editor is doing, in both copies of the status line.

        *spinner* says at the call site that this message is about work in
        flight. It no longer decides anything: whether the ring turns is
        decided by whether a worker is registered as running. Work that this
        editor does not start a thread for registers itself instead --
        _begin_round_trip, for a question put to the browser -- which is an
        override that says what is running rather than a flag that asserts it.

        It used to be the only thing that decided, and it defaulted to off --
        113 writes to it and 23 of them asking for it -- so a progress line
        written from inside a running worker ("step 12 of 40") turned the ring
        off, and the editor looked finished while xtb was still running.
        Someone who believes a calculation has stopped starts dragging into
        it, and then two things are moving the same structure.

        Giving it back its force as a floor was tried and measured to be
        wrong the other way round: the last frame of a drag arrives as a
        queued message asking for the ring, so a floor set by that message
        outlived the worker and the ring turned for good over a calculation
        that had finished. A ring nothing can stop is the same lie as one
        that never starts, so the registry is the only authority.
        """
        del spinner
        # Remembered, because the ring can come and go without anybody writing
        # a new line: a worker starting or ending re-renders what is already
        # there rather than inventing a message for it.
        state['mol_status_lines'] = tuple(lines)
        # In the plain lines it was given rather than the HTML they become: a
        # bug report is read by a person, and "es zappelt" is answered by what
        # the editor claimed at that second, not by its markup.
        record('status', lines=[str(one) for one in lines])
        _render_mol_status()

    def _render_mol_status():
        """Draw the remembered status line at the spinner state that holds now."""
        lines = state.get('mol_status_lines') or ()
        spinner = _busy_now()
        # Both copies always say the same thing; which one is on screen is the
        # overlay's business, not the caller's.
        rendered = [html.escape(str(line)) for line in lines if line not in (None, '')]
        # And what the page could not do with the frames, on the end of
        # whatever the row is saying -- see :func:`_trajectory_fault` for why
        # it is put on at the draw rather than written into the message.
        fault = _trajectory_fault()
        if fault:
            escaped = html.escape(fault)
            rendered = ([*rendered[:-1], f'{rendered[-1]} {escaped}']
                        if rendered else [escaped])
        spinner_html = (
            "<span class='delfin-busy' style='margin-right:6px; vertical-align:middle;' "
            "title='Working'></span>"
            if spinner else ''
        )
        text_html = '<br>'.join(rendered)
        if not spinner_html and not text_html:
            mol_status.value = ''
            mol_status_fs.value = ''
            return
        if spinner_html and text_html:
            first, *rest = rendered
            body = spinner_html + first
            if rest:
                body += '<br>' + '<br>'.join(rest)
        else:
            body = spinner_html + text_html
        rendered_html = (
            "<div style='font-family: monospace; white-space: pre-wrap; "
            "font-size: 13px; line-height: 1.35;'>"
            f"{body}</div>"
        )
        mol_status.value = rendered_html
        # The fullscreen copy is there to report work -- a spinner, a
        # trajectory, a result.  Asking for coordinates is not work, and in
        # fullscreen there is a structure on screen to look at, so the prompt
        # would be a permanent line saying nothing.
        #
        # The ring is not part of that: a prompt written while something is
        # running would otherwise take the ring off the overlay and leave it
        # in the inline line, and a spinner in one copy and none in the other
        # is the same bug seen from one side.
        prompt = any('enter XYZ' in str(line) or 'Load a structure' in str(line)
                     for line in lines)
        if not prompt:
            mol_status_fs.value = rendered_html
        elif spinner_html:
            mol_status_fs.value = (
                "<div style='font-family: monospace; white-space: pre-wrap; "
                "font-size: 13px; line-height: 1.35;'>"
                f"{spinner_html}</div>")
        else:
            mol_status_fs.value = ''

    def _clear_mol_status():
        # Both copies, the way _set_mol_status writes both. Clearing only the
        # small one left the overlay saying "Quick convert (single
        # structure)..." long after the structure was on screen: the finished
        # view comes through here, and in fullscreen that stale line was the
        # only thing the user had to go by.
        #
        # Through the same renderer, so that clearing the words does not clear
        # the ring: this is called from the middle of a conversion (the MANTA
        # loader takes the viewer over) and from every redraw, and a structure
        # arriving is not the same event as the work being over.
        _set_mol_status()

    # --- what is running, and therefore whether the ring is on screen -------
    #
    # The spinner used to be an argument of each status message: 113 writes to
    # the status line and 23 of them asked for it, so the other ninety turned
    # it off, including progress lines written from inside a worker that was
    # still running. That is not a list of missed call sites, it is the wrong
    # owner -- the ring is a fact about whether anything is running, and no
    # status message is in a position to know that. So the workers say, and
    # the renderer reads.
    _busy_serial = [0]

    def _busy_now():
        return bool(state.get('busy_jobs'))

    def _busy_begin(what):
        """Say that background work has started, and show it."""
        _busy_serial[0] += 1
        token = _busy_serial[0]
        state.setdefault('busy_jobs', {})[token] = str(what)
        _render_mol_status()
        return token

    def _busy_end(token):
        """Say that a piece of background work is over.

        Called from the worker's own thread, so the redraw is handed to the
        interface rather than done here.
        """
        jobs = state.get('busy_jobs') or {}
        if jobs.pop(token, None) is None:
            return
        # One turn of the interface behind, so that work handing over to more
        # work does not blink. An optimisation that has not converged presses
        # itself again from the callback that reports the round, and that
        # press arrives one turn after this one -- rendered immediately, the
        # ring went out and came back between two rounds of a press that never
        # stopped. Late is the safe direction: the ring stays up a fraction
        # longer than the work, never a fraction less.
        schedule_ui_update(lambda: schedule_ui_update(_render_mol_status))

    #: How long a question asked of the browser is given before the ring is
    #: taken down anyway. Reading a drawing takes milliseconds; a frame that
    #: has been folded away or reloaded never answers at all, and a ring that
    #: nothing can stop says a calculation is running when none is.
    _ROUND_TRIP_LEASH = 20.0

    def _end_round_trip(name):
        """The browser answered -- or something else made the question moot."""
        token = state.pop(f'{name}_round_trip', None)
        if token is not None:
            _busy_end(token)

    def _begin_round_trip(name, what, seconds=_ROUND_TRIP_LEASH):
        """A question asked of the page, which is work with no thread behind it.

        The drawing comes back through a hidden box the page writes into, so
        between the press and the answer there is no worker to count and the
        editor was, as far as the ring knew, idle -- while it was in fact
        waiting on something that may take a while or never come.
        """
        _end_round_trip(name)
        token = _busy_begin(what)
        state[f'{name}_round_trip'] = token

        def _leash():
            time.sleep(seconds)
            if state.get(f'{name}_round_trip') == token:
                state.pop(f'{name}_round_trip', None)
                _busy_end(token)

        # Not through _start_background: a leash is not work, and counting it
        # as work would hold the ring up for exactly as long as the leash.
        threading.Thread(target=_leash, daemon=True).start()

    def _start_background(work, what, *, guards=None, remember_in=None):
        """Run *work* on a thread, with the ring up for exactly its life.

        Every background worker in this editor goes through here, so the ring
        is on from before the thread starts until after the work returns --
        whatever any status line written in between happens to say.

        The count is lowered in a ``finally``. A worker that raises used to
        leave the ring turning for good: none of the fourteen threads was
        wrapped, and each of them cleared its own run flag inside a callback
        that an exception skips. *guards* are the state keys that finishing
        step would have cleared, given as {key: value}, so a failure leaves
        the editor pressable again rather than convinced a run is still going.

        *remember_in* is a state key the thread itself is put under before it
        starts, for the one caller that waits on the previous run to be gone
        before it doubles up on a shared login node.
        """
        token = _busy_begin(what)

        def _run():
            try:
                work()
            except Exception as exc:            # noqa: BLE001 - reported, not swallowed
                for key, value in (guards or {}).items():
                    state[key] = value
                schedule_ui_update(
                    _set_mol_status,
                    f'{what} stopped on an error: {exc.__class__.__name__}: '
                    f'{exc}')
            finally:
                _busy_end(token)

        thread = threading.Thread(target=_run, daemon=True)
        if remember_in is not None:
            # Before it starts, not after: the next press reads this key to
            # find the run it has to wait for.
            state[remember_in] = thread
        thread.start()
        return thread

    def _ensure_manip_bootstrap():
        # Once per page, not once per editor. The script is 159 KiB and does
        # nothing the second time -- it guards on its own version -- but with
        # two editors it was still being written out and shipped twice, and
        # each host kept its own "done" flag so neither could tell.
        if state.get('manip_bootstrap_done') or getattr(ctx, '_delfin_manip_sent', False):
            state['manip_bootstrap_done'] = True
            return
        try:
            ctx.run_js(submit_manip_bootstrap_js())
            state['manip_bootstrap_done'] = True
            try:
                ctx._delfin_manip_sent = True
            except Exception:
                pass
        except Exception:
            pass

    def _set_manip_toolbar_enabled(enabled):
        submit_fullscreen_btn.disabled = not enabled
        submit_select_btn.disabled = not enabled
        submit_manip_btn.disabled = not enabled
        submit_draw_btn.disabled = not enabled
        submit_element_dd.disabled = not enabled
        submit_adjust_h_btn.disabled = not enabled
        submit_manip_clear_btn.disabled = not enabled
        submit_relax_btn.disabled = not enabled
        submit_strength_slider.disabled = not enabled
        submit_pull_slider.disabled = not enabled
        submit_hand_dd.disabled = not enabled
        submit_play_speed.disabled = not enabled
        for widget in (submit_thermal_btn, submit_temperature,
                       submit_thermal_relax, submit_thermal_anchor_btn,
                       submit_topology_btn, submit_saddle_btn,
                       submit_saddle_from, submit_saddle_how,
                       submit_mode_dd, submit_mode_btn, submit_ends_btn,
                       submit_climb_btn, submit_shape_btn,
                       submit_path_from_btn):
            widget.disabled = not enabled
        submit_labels_btn.disabled = not enabled
        submit_sens_slider.disabled = not enabled
        submit_settle_btn.disabled = not enabled
        submit_auto_btn.disabled = not enabled
        submit_bond_btn.disabled = not enabled
        submit_unbond_btn.disabled = not enabled
        submit_hyb_auto_btn.disabled = not enabled
        submit_ff_dd.disabled = not enabled
        submit_optimize_btn.disabled = not enabled
        submit_optimize_all_btn.disabled = not enabled
        submit_internal_value.disabled = not enabled
        submit_internal_btn.disabled = not enabled
        submit_hold_btn.disabled = not enabled
        submit_hold_mode.disabled = not enabled
        submit_scan_btn.disabled = not enabled
        submit_manip_undo_btn.disabled = not enabled
        _refresh_undo_redo()
        submit_centre_btn.disabled = not enabled
        submit_reset_btn.disabled = not enabled
        submit_manip_toolbar.layout.display = 'flex' if enabled else 'none'
        if not enabled:
            if submit_select_btn.value:
                submit_select_btn.value = False
            if submit_manip_btn.value:
                submit_manip_btn.value = False

    # -- atom-selection / manipulation handlers -------------------------
    def _run_manip_js(code):
        try:
            ctx.run_js(code)
        except Exception:
            pass

    def _apply_manip_mode_js(mode):
        _ensure_manip_bootstrap()
        _run_manip_js(
            f'if(window.__delfinSubmitManip) '
            f'window.__delfinSubmitManip.setMode({json.dumps(submit_scope_id)}, '
            f'{json.dumps(mode)});'
        )

    def _mode_buttons_mutex(keep):
        """Only one mode at a time; the others switch themselves off."""
        for button in (submit_select_btn, submit_manip_btn, submit_draw_btn):
            if button is not keep and button.value:
                button.value = False

    def _refresh_draw_controls():
        drawing = bool(submit_draw_btn.value)
        submit_element_dd.layout.display = '' if drawing else 'none'

    def on_submit_select_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_select_btn.value)
        submit_select_btn.button_style = 'info' if active else ''
        if active:
            _mode_buttons_mutex(submit_select_btn)
        _apply_manip_mode_js('select' if active else 'off')

    def on_submit_manip_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_manip_btn.value)
        submit_manip_btn.button_style = 'info' if active else ''
        if active:
            _mode_buttons_mutex(submit_manip_btn)
        _apply_manip_mode_js('manipulate' if active else 'off')

    def on_submit_draw_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_draw_btn.value)
        submit_draw_btn.button_style = 'info' if active else ''
        if active:
            _mode_buttons_mutex(submit_draw_btn)
        _refresh_draw_controls()
        _apply_manip_mode_js('draw' if active else 'off')
        if active:
            on_submit_draw_choice(None)

    def on_submit_draw_choice(_change=None):
        """Hand over the element the browser draws with.

        Not the bond order: a drawn bond is single, and what it should be is
        decided afterwards by tapping the stick, where it can be seen. Having
        to choose in advance was a setting that was almost always wrong.
        """
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setDrawElement('
            f'{json.dumps(submit_scope_id)},{json.dumps(submit_element_dd.value)});'
        )

    def _set_ff_notes(notes):
        """Show what the force field had to approximate, under the viewer."""
        rendered = [html.escape(str(note)) for note in notes if str(note).strip()]
        if rendered and _server_method():
            # These come from the field that runs in the browser, which is UFF
            # whatever the method box says -- and read as though GFN had
            # produced them, which is how "GFN behaves exactly like UFF" gets
            # concluded from a panel that never claimed otherwise.
            label = _server_label(submit_ff_dd.value)
            rendered.insert(0, html.escape(
                f'These notes are about the live field, which is UFF: it runs '
                f'in the browser so that dragging follows the mouse. '
                f'{label} runs when Optimise is pressed, and says so in the '
                f'status line when it has.'))
        if not rendered:
            submit_ff_notes.value = ''
            return
        items = ''.join(
            f'<li style="margin:0 0 2px 0;">{note}</li>' for note in rendered
        )
        submit_ff_notes.value = (
            "<div style='font-size:12px; line-height:1.4; color:#5a6570; "
            "background:#f6f7f9; border:1px solid #e0e4e8; border-radius:4px; "
            "padding:6px 10px;'>"
            "<b style='color:#455a64;'>Force field notes</b>"
            f"<ul style='margin:4px 0 0 16px; padding:0;'>{items}</ul>"
            "</div>"
        )

    def _ensure_ff_bootstrap():
        # The same, for the 116 KiB force field.
        if state.get('ff_bootstrap_done') or getattr(ctx, '_delfin_ff_sent', False):
            state['ff_bootstrap_done'] = True
            return
        try:
            from .molecule_forcefield_js import molecule_ff_bootstrap_js
            ctx.run_js(molecule_ff_bootstrap_js())
            state['ff_bootstrap_done'] = True
            try:
                ctx._delfin_ff_sent = True
            except Exception:
                pass
        except Exception:
            pass

    def _structure_fingerprint(xyz):
        """Element column of an XYZ block -- what makes it the same molecule.

        A row counts only if what follows the symbol are three numbers. Four
        whitespace-separated words were enough before, and an XYZ comment line
        is free text: "Edited in DELFIN viewer" is four words, so it went into
        the fingerprint as though it were an atom called Edited. The molecule
        then looked like a different one every time a comment changed, the
        bonding was perceived again from the coordinates, and a bond stretched
        by dragging an atom away was gone.
        """
        out = []
        for line in (xyz or '').splitlines():
            row = line.split()
            if len(row) < 4:
                continue
            try:
                float(row[1]), float(row[2]), float(row[3])
            except ValueError:
                continue
            out.append(row[0])
        return tuple(out)

    def _perception_for(xyz):
        """Perceive the bonding once per structure and keep it.

        Bond orders are read from the geometry, and a twisted double bond stops
        looking like one: turning a C=C by 30 degrees is enough for perception
        to call it a single bond, which drops the barrier holding the two
        halves coplanar from 19.5 to 1.1 kcal/mol. Everything downstream then
        lets the double bond rotate freely, and nothing brings it back.

        Editing moves atoms; it must not be able to change what the molecule
        is. So the bonding is perceived from the structure as it arrived and
        reused until a genuinely different one is loaded.
        """
        from .molecule_forcefield import perceive_molecule

        fingerprint = _structure_fingerprint(xyz)
        cached = state.get('perceived')
        if cached and state.get('perceived_for') == fingerprint:
            return cached
        perceived = perceive_molecule(xyz)
        _apply_bond_edits(perceived)
        # After the bond edits, never before: rebuilding the typing molecule
        # sanitizes it, and sanitisation re-perceives hybridisation.
        _apply_hyb_overrides(perceived)
        state['perceived'] = perceived
        state['perceived_for'] = fingerprint
        return perceived

    def _apply_bond_edits(perceived):
        """Lay the user's bond corrections over what perception found.

        The correction has to reach the molecules the force-field parameters
        are read from, not only the bond list -- otherwise a drawn bond keeps
        the length it was drawn at instead of contracting.
        """
        from .molecule_forcefield import apply_bond_edits

        apply_bond_edits(perceived, state.get('bond_edits') or {})

    def _apply_hyb_overrides(perceived):
        """Force the hybridisations the user has chosen by hand.

        Bond orders are perceived from the geometry, and a double bond that is
        not seen leaves its carbon typed sp3: angles at 109.5 degrees, and a
        centre that puckers where it should stay flat.
        """
        from .molecule_forcefield import apply_hybridisation_overrides

        apply_hybridisation_overrides(perceived, state.get('hyb_overrides') or {})

    def _play_pace():
        """How fast the page walks the path, in milliseconds a frame.

        The slider counts frames a second and the player counts the delay
        between them.  The top of the slider is not a speed but "keep
        up", which is zero here and which the page answers by taking
        everything that has arrived in one go -- at 60 frames a second it
        was one frame per animation frame, so a burst of thirty answers
        put the picture half a second behind and it never caught up.
        """
        asked = int(submit_play_speed.value)
        pace = (0 if asked >= int(submit_play_speed.max)
                else max(1, int(round(1000.0 / max(1, asked)))))
        return pace

    def _frame_payload(run, **fields):
        """A write for the frame channel, named for its run and paced.

        The pace rides with the frames because the channel that carried
        it on its own -- run_js -- clears its output before it displays,
        and so races the start-up script that builds the player.  A page
        that lost that race played at the built-in 55 ms a frame however
        the slider was set, which is what "the top of the slider is not
        the fastest it goes" was.

        The slider is the *default* pace and not the only one.  A run may
        name its own in *fields*, and one does: a normal mode has a period
        rather than a rate, and the top of the slider -- "do not fall
        behind the calculation" -- is the whole animation inside one screen
        frame for a path that was finished before it was sent.  See
        :func:`_draw_the_mode`.
        """
        return json.dumps(dict(fields, run=run,
                               pace=fields.get('pace', _play_pace())))

    def _frame_run_is_current(run):
        """Whether a payload prepared for *run* may still be written.

        Every writer keeps the run number it started with, and three of
        them never asked again: a hand being followed, the settle after a
        release, and the scan.  Anything that moves the run number on
        while one of those has an answer in flight -- Optimise, a scan, a
        hand landing on the structure, an edit that abandons a run --
        left it writing frames about a structure that had been replaced.
        The player sees a run it does not know, resets, and walks the
        stale path over the new one: that is the trajectory shown twice,
        and the jump back to where it began.
        """
        return int(run or 0) == int(state.get('gfn_run', 0))

    def _note_the_run(run, by):
        """Write down that *by* has taken run number *run*.

        Said at every place that moves the counter, and said rather than
        worked out afterwards from the frames.  A run that claims a number
        and then never draws a frame is not an absence -- it is the shape two
        of this week's defects had, a walker failing its own guard on every
        write and throwing its whole path away silently and correctly.  Read
        back from the frames alone that case is indistinguishable from a
        number nobody ever took.

        And it is where the scan's profile goes.  The counter moves for two
        different reasons and only one of them ends a picture.  Something
        starts drawing over the structure -- the follow under a hand, the
        settle, an optimisation, a climb, a saddle search, a chain, the next
        scan -- and then the last walk's profile is no longer about what is on
        screen, and it goes at the moment the run starts rather than when it
        ends, because a picture beside a structure being pulled about has
        already stopped being true.  Or something is being *stopped* and the
        page's player has to be given a number it has never seen: that is
        ``press`` and ``abandoned``, claimed by Undo, Reset and the Stops, and
        it draws nothing at all.

        Measured on a twenty-four-point torsion, before that distinction was
        made here: the walk finished, the picture stood, and one press of Undo
        took it away over "Scanned: the highest point the walk crossed" --
        which is the walk the picture is of.  Undo goes back through the
        scan's own landmarks, so what arrives in the box is what decides that
        case, a few lines later in the same press (see
        :func:`_scan_plot_holds`).
        """
        record('run', v=int(run), by=by)
        if by not in ('press', 'abandoned'):
            _scan_plot_drop()
        return run

    def _claim_the_frame_run():
        """Take the next run number, for whatever is about to start drawing.

        The number is how a writer says its frames are about the structure on
        screen, and the protocol is that a run claims one when the run begins
        -- once, not at the moment the switch was pressed and not again while
        it walks.  The climb used to claim at the toggle and hold it, and the
        drag between the toggle and the first frame moved the counter twice:
        the hand being followed takes a number and the release takes another.
        The climb then failed its own guard on every write and threw away all
        94 frames of its own path, silently and correctly -- it really was not
        the current run.  It begins where the minimisation begins now, so the
        two cannot drift apart again.
        """
        run = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = run
        # Who took it, worked out here from what is walking rather than said
        # by the caller. The five callers of this are described by their exact
        # text in the tests that read this source -- a Stop must contain
        # ``_claim_the_frame_run()`` and be seen to -- and an argument they do
        # not need is not worth being the reason those tests are rewritten.
        # Where the distinction really carries something the site says it
        # itself: the follow, the settle, the scan, the saddle and the chain
        # each name themselves where they move the counter.
        _note_the_run(run,
                      'climb' if state.get('climb_run') is not None
                      else 'optimise' if state.get('optimize_run') is not None
                      else 'press')
        # A path put down by an earlier Stop and never claimed -- the page
        # never said where the picture stood -- belongs to a run that is over.
        # Left lying, it would be cut by the next Stop's frame number.
        state.pop('gfn_stopped_path', None)
        return run

    def _stream_frames(run_id, frames, *, final=False, follow=False,
                       least_apart=0.0):
        """Hand a path over while it is still being walked.

        One writer for both optimisers, because a path is a path: one walks
        down to a minimum and one walks up to a saddle, and neither of them
        needs its own way of being looked at.  Every frame, exactly once, and
        each one carried twice.

        The field is one slot, not a queue: a write that lands before the page
        has read the one before it replaces it, and those frames are gone.
        That is what an eight-frame tail was for -- it re-sent recent frames so
        a missed read still caught them -- but it was a *fixed* eight, and a
        walker makes frames far faster than the page is asked to look.  A
        benzene runs 23 cycles in a fraction of a second and a 149-atom chain
        260; everything between two reads beyond the last eight was never sent
        at all, and what reached the viewer was a sample of the path rather
        than the path.

        So the window starts where the *previous* window started rather than
        where it ended.  Every frame is therefore sent in two consecutive
        writes, which is the same insurance the tail gave, and nothing is
        skipped however fast the frames arrive.  It stays bounded -- a write
        carries at most two reads' worth -- and the coordinates are rounded to
        four decimals, which is a thousandth of a bond length and half the
        JSON.

        *final* is the write at the end of the run, and it goes out even when
        it carries nothing new.  Without it the last window is the one window
        sent only once -- the run ends before another write can repeat it --
        so a single missed read at exactly that moment leaves the picture
        short of the geometry the box holds.  That is the end of the path,
        which is the part that has to land.

        *follow* marks frames that no pressed Optimise stands behind.  The
        page abandons a queue when that switch is up, which is right for a
        minimisation the user stopped and fatal for anything else: a climb
        runs with Optimise up from beginning to end, so without this the page
        throws its path away as soon as it has queued any of it.

        *least_apart* holds two writes apart for a walker that makes frames
        faster than any page is asked to look: the climb makes a hundred a
        second where xtb's optimiser makes a few, and every write is a
        message.
        """
        walked = list(frames)
        if state.get('gfn_push_run') != run_id:
            state['gfn_push_run'] = run_id
            state['gfn_push_start'] = 0
            state['gfn_push_end'] = 0
            state['gfn_push_at'] = 0.0
        now = time.monotonic()
        if not final:
            if len(walked) <= int(state.get('gfn_push_end') or 0):
                return                  # nothing new since the last write
            if now - float(state.get('gfn_push_at') or 0.0) < least_apart:
                return
        start = int(state.get('gfn_push_start') or 0)
        state['gfn_push_start'] = int(state.get('gfn_push_end') or 0)
        state['gfn_push_end'] = len(walked)
        state['gfn_push_at'] = now
        fresh = [[round(float(one), 4) for one in frame]
                 for frame in walked[start:]]

        def _write(rows=fresh, first=start, last=bool(final)):
            # Asked at the write and not at the answer.  A run that has been
            # replaced does not draw: an interrupted one has frames in hand
            # when it is told to stop, and writing them afterwards played the
            # abandoned path over the structure the user had just made.
            if not _frame_run_is_current(run_id):
                return
            fields = {'from': first, 'frames': rows}
            if follow:
                fields['follow'] = 1
            if last:
                # Named, so this write differs from the one before it even
                # when it carries the same frames: the field is a widget
                # value and traitlets says nothing when a value is written
                # again unchanged, so the repeat that is the whole point of a
                # final write would never leave the kernel.
                fields['final'] = 1
            submit_gfn_frame.value = _frame_payload(run_id, **fields)

        schedule_ui_update(_write)

    def _halt_the_frames(run_id):
        """Tell the page the run it is playing has been switched off.

        One halt for both optimisers, because a Stop means one thing: the
        picture keeps the frame it is showing and drops what was computed past
        it.  The page cannot work that out for itself -- a run that ends and a
        run that is stopped look identical from there -- so it has to be told,
        and the climb never told it.  Measured on a climb stopped at frame 13
        of 117 with the slider at twelve a second: 509 more draws afterwards,
        frames 14 through 116, nine seconds of trajectory playing out over a
        structure nothing was calculating any more.  That is "ein Stop von
        Climb to TS laesst die ganze restliche Trajektorie noch nachspielen".

        Addressed to the run that is playing and not to whatever is current.
        The number has already moved on by the time this goes out -- moving it
        is what refuses the frames the stopped run still has in hand -- and a
        halt naming the newer run would clear *its* queue instead.

        Once per run.  A hand that interrupts marks it sent without sending
        anything, because an abandoned run is not a stopped one: its frames
        must not be landed on at all, and it says so with a run number the
        page has never seen.
        """
        if state.get('gfn_halt_sent'):
            return False
        state['gfn_halt_sent'] = True
        # Through _frame_payload like every other write on this channel, so
        # the pace rides with it too.  Written by hand it was the one payload
        # that carried none.
        schedule_ui_update(
            lambda: setattr(submit_gfn_frame, 'value',
                            _frame_payload(run_id, halt=1, frames=[])))
        return True

    def _install_gfn_frame_watcher():
        """Teach the page to play the trajectory, once.

        Two things were wrong with sending one frame at a time.  ipywidgets
        writes a new value into the DOM without firing an event, so there is
        nothing to listen for -- the value has to be read on a timer.  And the
        kernel produces frames faster than the reader can look, so all but the
        last of each burst were never seen: the structure jumped instead of
        moving.

        The field carries the whole trajectory instead, and the page keeps its
        place in it.  Nothing can be missed, because a frame that was written
        while the reader was busy is still there when it looks again.  Playback
        interpolates between frames -- twenty computed steps a second shown as
        continuous motion.  The positions in between are drawn, not calculated;
        every frame the structure actually passed through is one of the ends.
        """
        if state.get('gfn_watcher_installed'):
            return
        state['gfn_watcher_installed'] = True
        _emit_watcher_js(
            '(function(){\n'
            '  var scope=' + json.dumps(submit_scope_id) + ';\n'
            '  window.__delfinGfnPlay=window.__delfinGfnPlay||{};\n'
            '  if(window.__delfinGfnPlay[scope]) return;\n'
            '  var play={queue:[],at:0,started:0,last:null,seen:0,'
            'run:null,top:0,stopped:0,pace:' + json.dumps(_play_pace())
            + '};\n'
            '  window.__delfinGfnPlay[scope]=play;\n'
            '  var STEP_MS=55;\n'
            '  function stepMs(){\n'
            '    /* The pace the user asked for, in milliseconds a frame.  It\n'
            '       used to speed up whenever the queue grew -- 75 frames\n'
            '       arriving in 0.4 s would otherwise take 4 s to show, and the\n'
            '       picture trailed the calculation further and further.  But\n'
            '       trailing is what a slow setting is for: the whole path is\n'
            '       in the page now, so where the picture has got to is where\n'
            '       the user is, and grabbing an atom there keeps that frame\n'
            '       and drops what was computed past it.  A rule that quietly\n'
            '       sped the playback back up would take that away. */\n'
            '    var n=play.queue.length;\n'
            '    /* Zero is not a very small delay -- it is "as fast as this\n'
            '       machine can", and the loop below takes the whole queue in\n'
            '       one go for it.  The top of the slider means that, and so\n'
            '       does a hand on the structure: a drag is not a replay, and\n'
            '       a frame held back during one is a frame the picture is\n'
            '       behind the hand by.  Queued and then dropped at the\n'
            '       release, which is what it did, that is the springing back\n'
            '       the user sees -- and why it looks more alive once the\n'
            '       mouse is up.\n'
            '\n'
            '       Spreading a held frame over its interval was tried and is\n'
            '       not done, for a reason that is not about smoothness: the\n'
            '       drag pushes the *viewer\'s* coordinates back to the kernel\n'
            '       ten times a second as the structure the next xtb round\n'
            '       starts from, so a position halfway between two answers\n'
            '       would be handed to the optimiser as a geometry.  A drawn\n'
            '       in-between may be looked at; it may not be computed on. */\n'
            '    if(play.held) return 0;\n'
            '    /* One answer at a time -- a hand being followed -- is drawn\n'
            '       over exactly the time the next answer takes to come, or the\n'
            '       picture sits still between them.\n'
            '\n'
            '       This stands above the slider, and that is the whole of it:\n'
            '       the slider says how fast to walk a path that has already\n'
            '       been computed, and a followed hand has no such path -- it\n'
            '       is paced by how fast xtb answers, which no setting makes\n'
            '       faster.  Below the slider, as it was, the top of the slider\n'
            '       -- zero, "take everything that has arrived now" -- reached\n'
            '       it first and drew each raw answer the instant it landed.\n'
            '       Measured in node over one synthetic drag, single answers\n'
            '       70-170 ms apart: paced by the arrival gap the picture is\n'
            '       redrawn every 17.5 ms and moves 0.14 of a frame each time;\n'
            '       reached by the slider at its top instead it is redrawn\n'
            '       every 117 ms and moves a whole frame, standing still for a\n'
            '       sixth of a second between jumps.  That is the jitter.\n'
            '\n'
            '       Which of the two it is, is not a matter of counting the\n'
            '       queue: it is whichever of them is the *slower*.  The slider\n'
            '       asks not to be shown more than so many frames a second and\n'
            '       the arrivals cannot offer more than they offer, and obeying\n'
            '       both means taking the longer interval.  At the top of the\n'
            '       slider -- zero, no limit asked for -- that is the arrival\n'
            '       rate, which is what a followed hand wants; at Speed 12 over\n'
            '       a relaxation pouring out answers it is the slider, and the\n'
            '       queue grows, which is what asking to watch slowly means.\n'
            '\n'
            '       Three in the queue was the guard, and three is exactly the\n'
            '       size of a real message: the kernel reads the log twenty\n'
            '       times a second and hands over the two or three answers\n'
            '       walked since.  So the guard fired on every second message\n'
            '       and the two rules took turns -- creep at the follow pace,\n'
            '       then the whole queue emptied into one animation frame.\n'
            '       Measured over real GFN2 settles replayed at the arrival\n'
            '       times xtb produced, largest drawn step over median drawn\n'
            '       step: 16.2 with the count as the guard, 1.0 with the two\n'
            '       paces compared.  A queue that is running away is still\n'
            '       handed back -- a tab that stopped being drawn for a second\n'
            '       has fifty answers waiting and wants catching up, not\n'
            '       following -- but the bound stands well above a message\n'
            '       rather than at it. */\n'
            '    if(play.follow&&play.gap&&n<=12){\n'
            '      var asked=(play.pace!==undefined&&play.pace!==null)\n'
            '        ? play.pace : STEP_MS;\n'
            '      return asked?Math.max(asked,play.gap):play.gap;\n'
            '    }\n'
            '    if(play.pace!==undefined&&play.pace!==null) return play.pace;\n'
            '    if(n>60) return 8;\n'
            '    if(n>25) return 20;\n'
            '    if(n>10) return 35;\n'
            '    return STEP_MS;\n'
            '  }\n'
            '  function inScope(selector){\n'
            '    /* Every element carrying this editor\'s scope, not the first.\n'
            '       Fullscreen adds one, and the 2D drawing frame carries the\n'
            '       scope too and sits in the other column holding none of\n'
            '       these parts -- so taking the first found the drawing frame\n'
            '       and the report channel went nowhere at all. */\n'
            '    var roots=document.querySelectorAll("."+scope);\n'
            '    for(var i=0;i<roots.length;i++){\n'
            '      var hit=roots[i].querySelector(selector);\n'
            '      if(hit) return hit;\n'
            '    }\n'
            '    return null;\n'
            '  }\n'
            '  function inScopeAll(selector){\n'
            '    /* Every match, not the first.  A class worn by one element\n'
            '       can be read with inScope; a class worn by three cannot,\n'
            '       and the switch class is worn by three -- Optimize,\n'
            '       Optimize all and Climb to TS.  Taking the first found\n'
            '       asked only about Optimize, so a path started by either of\n'
            '       the other two read as a run whose switch was up. */\n'
            '    var out=[];\n'
            '    var roots=document.querySelectorAll("."+scope);\n'
            '    for(var i=0;i<roots.length;i++){\n'
            '      var hits=roots[i].querySelectorAll(selector);\n'
            '      for(var k=0;k<hits.length;k++) out.push(hits[k]);\n'
            '    }\n'
            '    var loose=document.querySelectorAll("."+scope+selector);\n'
            '    for(var j=0;j<loose.length;j++) out.push(loose[j]);\n'
            '    return out;\n'
            '  }\n'
            '  function setWall(w,reach){\n'
            '    var api=window.__delfinSubmitManip;\n'
            '    if(!api||!api.setThermalWall) return;\n'
            '    try{ api.setThermalWall(scope,w,reach); }catch(e){}\n'
            '  }\n'
            '  function readWall(){\n'
            '    /* Where the hand may no longer go.  Read on this clock and\n'
            '       not handed over by run_js: that clears its output before\n'
            '       displaying, so of two calls in quick succession only the\n'
            '       second survives -- and the kernel writes this about ten\n'
            '       times a second while a drag is running.  Sent that way it\n'
            '       never arrived, and a ring came apart under a budget that\n'
            '       said it could not. */\n'
            '    var field=inScope(".submit-gfn-wall input, '
            '.submit-gfn-wall textarea");\n'
            '    if(!field) return;\n'
            '    var text=field.value||"";\n'
            '    if(text===play.wallText) return;\n'
            '    play.wallText=text;\n'
            '    if(!text){ setWall(null); return; }\n'
            '    try{ var got=JSON.parse(text)||{};\n'
            '      setWall(got.wall||null, got.reach||0); }\n'
            '    catch(e){}\n'
            '  }\n'
            '  function read(arrivedAt){\n'
            '    /* Fullscreen moves the viewer into an overlay that carries\n'
            '       the same scope class, and the frame field is not one of the\n'
            '       things it takes -- so looking only inside the first element\n'
            '       with that class found the overlay and no field, and the\n'
            '       playback that worked in the small view showed nothing in\n'
            '       the big one.  Every element with the class is tried.\n'
            '       Nothing stands behind them any more: the document did, and\n'
            '       with a second editor on the page the first field in the\n'
            '       document belongs to whichever editor was written out first\n'
            '       -- so a Builder that lost sight of its own would have\n'
            '       played the Submit tab\'s trajectory into its viewer. Being\n'
            '       told there is none is the right answer. */\n'
            '    var field=inScope('
            '".submit-gfn-frame textarea, .submit-gfn-frame input");\n'
            '    if(!field) return;\n'
            '    var text=field.value||"";\n'
            '    /* The same text as last time is the same frames. Nothing\n'
            '       below has anything to do a second time -- frames past the\n'
            '       ones already queued is all it takes from here -- and\n'
            '       parsing it again ran 60 times a second on a payload that\n'
            '       had not changed, megabytes per second of work on the\n'
            '       thread that also has to draw. */\n'
            '    if(text===play.parsedText) return;\n'
            '    play.parsedText=text;\n'
            '    if(!text){ play.queue=[]; play.seen=0; play.last=null;'
            ' play.run=null; play.stopped=0; return; }\n'
            '    var data=null;\n'
            '    try{ data=JSON.parse(text); }catch(e){ return; }\n'
            '    /* The pace travels with the frames.  It reached the page\n'
            '       only through run_js before, which clears its output\n'
            '       before it displays -- so the one write that carries it\n'
            '       races the start-up script that creates this player, and\n'
            '       a page that lost that race ran at the built-in 55 ms a\n'
            '       frame while the slider said the top.  Carried on the\n'
            '       payload it cannot be lost: the writes that matter are\n'
            '       the ones the run is already making. */\n'
            '    if(data&&data.pace!==undefined&&data.pace!==null)\n'
            '      play.pace=data.pace;\n'
            '    var frames=(data&&data.frames)||[];\n'
            '    var run=(data&&data.run)||0;\n'
            '    /* A run number only ever goes up.  The kernel hands one\n'
            '       out for every optimisation, drag, settle and scan, and\n'
            '       each of those keeps the number it was given until it\n'
            '       stops answering -- which can be a whole xtb round after\n'
            '       the run it belongs to was replaced.  A write from one of\n'
            '       those arrives naming a run the page has already left,\n'
            '       and the check below reads it as a *new* run: it resets,\n'
            '       and walks the abandoned path over the structure the user\n'
            '       has just made.  That is the trajectory seen twice.\n'
            '       Refusing anything older than the newest run seen makes\n'
            '       that impossible however late the write is. */\n'
            '    if(play.top&&run&&run<play.top) return;\n'
            '    if(run>(play.top||0)) play.top=run;\n'
            '    if(data&&data.halt){\n'
            '      /* The run was switched off.  Playing out the queue after\n'
            '         that is the picture carrying on without the thing it is\n'
            '         a picture of.\n'
            '\n'
            '         The count of what has been shown stays.  Set from the\n'
            '         halt payload -- which carries no frames, so it set\n'
            '         zero -- the next write for the same run looked\n'
            '         entirely new, and its window began mid-stream:\n'
            '         measured, a Stop at frame 69 then drew 57, 59, 61, 63,\n'
            '         65, 67, 69.  Time running backwards and the tail\n'
            '         played twice, on every Stop.\n'
            '\n'
            '         And it is the halted run that stops, not whichever\n'
            '         one is playing.  A worker notices the switch went up\n'
            '         a whole xtb round late, by which time the next run can\n'
            '         be under way -- and a halt nobody addressed cleared\n'
            '         that one\'s queue and reported it stopped. */\n'
            '      if(run&&play.run!==null&&run!==play.run) return;\n'
            '      play.queue=[]; play.stopped=1;\n'
            '      if(!play.toldStop){ play.toldStop=1;\n'
            '        /* Which frame is on screen, and of which run.  Stopping\n'
            '           keeps that one: frames xtb had already computed but\n'
            '           nobody had seen are not what the user stopped at.  The\n'
            '           run rides with it for the same reason it rides with a\n'
            '           grab -- a count alone indexes any path at all. */\n'
            '        say("stopped at frame "+(play.shown||0)\n'
            '            +" of run "+(play.run||0)); }\n'
            '      return;\n'
            '    }\n'

            '    /* Whether these frames belong to a molecule following a hand\n'
            '       rather than to a minimisation.  The two are told apart\n'
            '       because Optimise is not pressed during a follow, and the\n'
            '       check that abandons a queue when that switch goes up would\n'
            '       otherwise throw every followed frame away. */\n'
            '    play.follow=(data&&data.follow)?1:0;\n'
            '    if(run!==play.run){\n'
            '      /* A new run. Without this the count of frames already\n'
            '         played carried over, so a shorter run than the one\n'
            '         before it played nothing at all -- which is what made\n'
            '         the playback look like it worked only sometimes. */\n'
            '      if(play.queue.length&&!data.abandoned){\n'
            '        /* Land the run that is ending on its last frame first.\n'
            '           Dropped there, the viewer keeps whatever it had drawn\n'
            '           while the kernel keeps the geometry it computed, and\n'
            '           the two drift apart: the next drag then hands over a\n'
            '           structure that is behind, and the relaxation nobody saw\n'
            '           is walked again -- which is every earlier drag being\n'
            '           made a second time.\n'
            '\n'
            '           Unless the run was abandoned, and then its queue is\n'
            '           exactly what must not be drawn: the kernel threw\n'
            '           those frames away because the user has just changed\n'
            '           the structure under them.  Landing on the newest of\n'
            '           them put the viewer on a geometry nobody has any\n'
            '           more, and nothing wrote over it afterwards. */\n'
            '        show(play.last,play.queue[play.queue.length-1],1);\n'
            '      }\n'
            '      play.run=run; play.seen=0; play.queue=[]; play.last=null;\n'
            '      play.shown=0; play.toldStop=0; play.complete=0;\n'
            '      play.stopped=0; play.home=null;\n'
            '    }\n'
            '    /* Where the picture has to be put back to if this run is cut\n'
            '       short.\n'
            '\n'
            '       Almost every run on this channel draws geometries somebody\n'
            '       is going to be left with: an optimiser walks to a minimum\n'
            '       and every frame of it is a structure that was computed, so\n'
            '       stopping on one of them is a choice the user made.  One run\n'
            '       is not like that.  A normal mode drawn out of a saddle is a\n'
            '       *picture* -- the structure displaced along an eigenvector,\n'
            '       which nobody chose and no method computed -- and the frame\n'
            '       it happens to be showing when a hand lands must never\n'
            '       become the structure.  A run that draws such frames sends\n'
            '       the geometry to return to with them, and the grab below\n'
            '       draws it before it lets go of the queue, so the coordinates\n'
            '       the page pushes back are the ones the box already holds.\n'
            '       Drawn here on the page and not asked of the kernel,\n'
            '       because a round trip is tens of milliseconds and the drag\n'
            '       pushes ten times a second: asked, the first push of every\n'
            '       drag would carry a displaced geometry. */\n'
            '    if(data&&data.home) play.home=data.home;\n'
            '    /* A run that was stopped stays stopped.  The kernel refuses\n'
            '       to write for it once the switch is up; this is the same\n'
            '       refusal on the page, for the write that was already in\n'
            '       flight when Stop was pressed.  Taken up, it walks the\n'
            '       very path the user stopped. */\n'
            '    if(play.stopped) return;\n'
            '    /* The whole path is in hand: play it out whatever the switch\n'
            '       does, because the switch going up is how a finished run\n'
            '       announces itself -- the kernel turns it off the moment the\n'
            '       last round lands, and the picture is still walking.  Set\n'
            '       above the run check instead, it was cleared by the very\n'
            '       write that carried it whenever a run was short enough to\n'
            '       finish in one message. */\n'
            '    if(data&&data.final) play.complete=1;\n'
            '    /* Where in the run these frames start.  A long run sends the\n'
            '       tail rather than the whole path -- every write is a message\n'
            '       and the whole path grows without end -- so counting from\n'
            '       the front of the message would replay frames already shown\n'
            '       and then stop showing new ones altogether. */\n'
            '    var from=(data&&data.from)||0;\n'
            '    if(play.follow&&play.seen===0&&frames.length>1){\n'
            '      /* A live run arriving from the start: only its newest frame\n'
            '         is worth anything.  The ones before it describe where the\n'
            '         structure was on the way here -- the drag that has just\n'
            '         happened, or a relaxation of it -- and playing those is\n'
            '         showing the user their own past.  What is wanted is the\n'
            '         frame that is current. */\n'
            '      from=from+frames.length-1;\n'
            '      frames=[frames[frames.length-1]];\n'
            '    }\n'
            '    if(from+frames.length>play.seen){\n'
            '      var start=Math.max(0,play.seen-from);\n'
            '      /* A frame arriving into an empty queue starts its own\n'
            '         interpolation.  Left over from the last one, the clock\n'
            '         reads far past a step and the frame is jumped to instead\n'
            '         of moved to -- which is every frame of a follow, where\n'
            '         they arrive further apart than a step. */\n'
            '      if(!play.queue.length) play.started=0;\n'
            '      /* How long the machine took over one answer, which is how\n'
            '         long one answer has to be drawn over.\n'
            '\n'
            '         Per answer and not per message.  The kernel reads xtb\'s\n'
            '         log twenty times a second and hands over everything\n'
            '         walked since, so one write regularly carries two or three\n'
            '         answers; paced as though it carried one, the picture drew\n'
            '         a third of what arrived, the queue filled, and at four in\n'
            '         hand the rule below took it back and emptied the whole\n'
            '         queue into a single animation frame.  What that\n'
            '         alternation looks like is not slowness but twitching, and\n'
            '         it lands hardest on the hydrogens, which move furthest\n'
            '         per answer: measured over a real GFN2 settle of an\n'
            '         ibuprofen -- 46 answers, 2.8 to a message, 53 ms apart --\n'
            '         the picture crept 0.018 A, 0.018 A, then jumped 0.24 A\n'
            '         and stood still, ten times a second, with the jump led by\n'
            '         a methyl hydrogen every time.  With the interval shared\n'
            '         out among the answers that came in it, and the guard\n'
            '         below asked of the clock instead of the queue, the same\n'
            '         settle is drawn in even steps of a frame each at the rate\n'
            '         the answers arrive.\n'
            '\n'
            '         One answer at a time -- the case this rule was written\n'
            '         for -- is untouched, because the divisor is then one.\n'
            '         The floor is a few milliseconds rather than 24: it is the\n'
            '         time one answer takes now, and answers that really do\n'
            '         arrive faster than the screen redraws must be allowed to,\n'
            '         or the queue grows for no reason other than the floor. */\n'
            '      if(play.arrived){\n'
            '        var many=Math.max(1,frames.length-start);\n'
            '        var measured=Math.min(600,Math.max(4,\n'
            '          (arrivedAt-play.arrived)/many));\n'
            '        /* Averaged, not taken raw.  Answers do not arrive evenly,\n'
            '           and pacing each frame by the last interval alone makes\n'
            '           the drawing speed jump about as much as the arrivals\n'
            '           do -- which is felt as jerkiness even when nothing is\n'
            '           being missed. */\n'
            '        play.gap=play.gap?(play.gap*0.6+measured*0.4):measured;\n'
            '      }\n'
            '      play.arrived=arrivedAt;\n'
            '      for(var i=start;i<frames.length;i++) play.queue.push(frames[i]);\n'
            '      play.seen=from+frames.length;\n'
            '      say("received "+play.seen+" frames");\n'
            '    }\n'
            '    /* The whole path is kept, however far behind it puts the\n'
            '       picture.  Twenty was the bound, on the grounds that xtb\n'
            '       makes frames faster than any frame rate shows them and a\n'
            '       queue that grows leaves the picture permanently behind the\n'
            '       calculation.  Measured in a browser against a sixty-frame\n'
            '       path arriving in bursts of twenty: thirty-five frames drawn\n'
            '       and twenty-five never shown at all -- the oldest of each\n'
            '       burst, thrown away to keep the bound.\n'
            '       Being behind is what the pace control asks for.  Where the\n'
            '       picture has got to is where the user is: taking hold of an\n'
            '       atom keeps the frame on screen and drops what was computed\n'
            '       past it, so the frames in front of it are the thing being\n'
            '       chosen between, not a backlog to be cleared.  A ceiling is\n'
            '       kept far above any real path, because a queue that grows\n'
            '       without one is a leak rather than a feature. */\n'
            '    if(play.queue.length>100000){\n'
            '      /* The oldest go, and the count of where the picture\n'
            '         has got to moves on by exactly as many.  It did not,\n'
            '         and that count is what a grab hands the kernel as the\n'
            '         frame to keep: measured over a cap of ten, the picture\n'
            '         stood at frame 70 and reported frame 1, so taking hold\n'
            '         of an atom would have cut the run 69 frames behind\n'
            '         what was on the screen -- the structure jumping back\n'
            '         to a geometry nobody had seen for a minute. */\n'
            '      play.shown=(play.shown||0)+(play.queue.length-100000);\n'
            '      play.queue=play.queue.slice(-100000);\n'
            '    }\n'
            '  }\n'
            '  function say(text){ send("gfnplay", text); }\n'
            '  function send(verb,text){\n'
            '    /* The page reports back through the command bridge the editor\n'
            '       already has, so a playback that does not appear says why by\n'
            '       itself instead of being read out of a console. */\n'
            '    try{\n'
            '      var wrap=inScope(".submit-cmd-sync");\n'
            '      var input=wrap&&wrap.querySelector("input, textarea");\n'
            '      if(!input) return;\n'
            '      play.serial=(play.serial||0)+1;\n'
            '      var proto=(input.tagName==="TEXTAREA")\n'
            '        ? window.HTMLTextAreaElement.prototype\n'
            '        : window.HTMLInputElement.prototype;\n'
            '      var setter=Object.getOwnPropertyDescriptor(proto,"value");\n'
            '      var line=verb+":"+play.serial+":"+text;\n'
            '      if(setter&&setter.set) setter.set.call(input,line);\n'
            '      else input.value=line;\n'
            '      input.dispatchEvent(new Event("input",{bubbles:true}));\n'
            '      input.dispatchEvent(new Event("change",{bubbles:true}));\n'
            '    }catch(e){}\n'
            '  }\n'
            '  function show(a,b,t){\n'
            '    if(!window.__delfinSubmitManip||'
            '!window.__delfinSubmitManip.setPositions){\n'
            '      if(!play.toldNoApi){ play.toldNoApi=1;'
            ' say("no setPositions on the page"); }\n'
            '      return;\n'
            '    }\n'
            '    var out=new Array(b.length);\n'
            '    if(!a||a.length!==b.length){ out=b; }\n'
            '    else { for(var i=0;i<b.length;i++) out[i]=a[i]+(b[i]-a[i])*t; }\n'
            '    var ok=false;\n'
            '    try{ ok=window.__delfinSubmitManip.setPositions('
            'scope,out,heldSerials()); }\n'
            '    catch(e){ ok=false; }\n'
            '    play.drawn=(play.drawn||0)+(ok?1:0);\n'
            '    if(!ok&&!play.toldNoDraw){ play.toldNoDraw=1;'
            ' say("setPositions did not draw"); }\n'
            '    if(ok&&!play.toldDrawing){ play.toldDrawing=1;'
            ' say("drawing"); }\n'
            '  }\n'
            '  function grabbed(){\n'
            '    /* Whether an atom is being moved right now.  A playback that\n'
            '       keeps writing positions during a drag puts the atom back\n'
            '       where xtb had it once per animation frame, so it cannot be\n'
            '       moved at all -- the drag has to own the picture while it\n'
            '       lasts.  A rectangle being pulled over the molecule selects\n'
            '       and moves nothing, so it is not a grab.\n'
            '\n'
            '       And neither does a tap.  A press that has not passed the\n'
            '       slop is a press that has moved nothing: translate and\n'
            '       rotate both do their work under movedEnough, so the\n'
            '       structure is untouched until it flips.  Counted as a grab\n'
            '       anyway, a tap sent gfngrab and gfnfree a tenth of a second\n'
            '       apart -- and a grab interrupts whatever is running and the\n'
            '       release arms it again, so a running relaxation or climb\n'
            '       stuttered and started over.  That was always true of a\n'
            '       click; it matters now because a tap is a gesture the user\n'
            '       makes on purpose, to name the atom a climb is aimed at.\n'
            '\n'
            '       Draw is not held to it: there a press that does not move\n'
            '       still places an atom, so it really is a change. */\n'
            '    var held=(window._submitManipStateByScope||{})[scope];\n'
            '    var drag=held&&held.drag;\n'
            '    if(!drag) return false;\n'
            '    if(drag.kind==="translate"||drag.kind==="rotate")\n'
            '      return !!drag.movedEnough;\n'
            '    return drag.kind==="draw";\n'
            '  }\n'
            '  function heldSerials(){\n'
            '    /* The atoms the hand is moving.  Coordinates that come back\n'
            '       describe where they were when they were sent, and the\n'
            '       cursor has moved on since -- so those are the ones the\n'
            '       playback must leave alone.\n'
            '\n'
            '       Unless the hand is a force, and then they are exactly the\n'
            '       ones it must draw: how far the atom got is the answer to\n'
            '       the question the drag asked, and keeping it under the\n'
            '       cursor is how a pull was made to look like a move. */\n'
            '    var st=(window._submitManipStateByScope||{})[scope];\n'
            '    if(st&&st.ffPull) return [];\n'
            '    return (st&&st.drag&&st.drag.targets)||[];\n'
            '  }\n'
            '  function followIsOn(){\n'
            '    /* Relax, read off the page rather than asked of the kernel --\n'
            '       the same reason as the switch below. */\n'
            '    var holder=inScope(".submit-gfn-follow")'
            '||(document.querySelector("."+scope+".submit-gfn-follow"));\n'
            '    if(!holder) return false;\n'
            '    var btn=(holder.tagName==="BUTTON")?holder'
            ':holder.querySelector("button");\n'
            '    if(!btn) return false;\n'
            '    return btn.classList.contains("mod-active");\n'
            '  }\n'
            '  function switchIsOn(){\n'
            '    /* Whether a walk of this structure still has a switch behind\n'
            '       it.  One question for all of them: Optimize walks down,\n'
            '       Optimize all walks down through a set and Climb to TS\n'
            '       walks up, and what the page does when the switch goes up\n'
            '       is the same in every case -- drop what was computed past\n'
            '       the picture and keep the frame on screen.  The climb was\n'
            '       not on this list, so it had to be marked as a followed\n'
            '       hand to survive here at all, and that mark is what took\n'
            '       its playback off the slider.\n'
            '\n'
            '       ipywidgets marks a pressed toggle with mod-active.  Reading\n'
            '       it here is instant; asking the kernel costs a round trip,\n'
            '       and the picture ran on for the length of it. */\n'
            '    var held=inScopeAll(".submit-optimize-switch");\n'
            '    if(!held.length) return true;\n'
            '    for(var i=0;i<held.length;i++){\n'
            '      var btn=(held[i].tagName==="BUTTON")?held[i]\n'
            '        :held[i].querySelector("button");\n'
            '      if(btn&&btn.classList.contains("mod-active")) return true;\n'
            '    }\n'
            '    return false;\n'
            '  }\n'
            '  function frame(now){\n'
            '    /* An atom picked up while xtb is running: the kernel is told\n'
            '       at the grab rather than at the release, because a GFN2 run\n'
            '       is thirteen seconds and every one of them would be spent\n'
            '       minimising a structure the user is in the middle of\n'
            '       changing.  The queue goes with it -- those frames belong to\n'
            '       the geometry that has just been altered. */\n'
            '    var held=grabbed();\n'
            '    if(held!==!!play.held){\n'
            '      play.held=held?1:0;\n'
            '      if(held){\n'
            '        /* The frames in front of the picture were computed for a\n'
            '           structure the hand is now changing, so they go.  Which\n'
            '           frame the picture stands on goes with the message: the\n'
            '           kernel keeps that one as the structure, and what xtb\n'
            '           had run on past it is thrown away.  Sent without it,\n'
            '           the kernel knew a hand had arrived and not where -- so\n'
            '           the geometry on screen lived only in the browser until\n'
            '           a drag happened to push it back. */\n'
            '        /* Unless this run was drawing a picture rather than a\n'
            '           path, and then the frame the hand landed on is exactly\n'
            '           what must not be kept: it is the structure displaced\n'
            '           along a mode, and the drag is about to push whatever\n'
            '           the viewer holds back to the kernel ten times a\n'
            '           second.  Put back first, then let go of the queue. */\n'
            '        if(play.home) show(play.last,play.home,1);\n'
            '        play.queue=[]; play.last=null;\n'
            '      }\n'
            '      else {\n'
            '        /* Let go of.  What is still queued describes the drag:\n'
            '           each of those frames carries the dragged atom where the\n'
            '           hand had it when the frame was computed, and while the\n'
            '           hand was down they were drawn around it.  Drawn now,\n'
            '           with nothing held any more, they walk that atom back\n'
            '           through the whole drag in front of the user.  Land on\n'
            '           the newest -- which is where the hand actually left it\n'
            '           -- and drop the rest. */\n'
            '        if(play.queue.length){\n'
            '          show(play.last,play.queue[play.queue.length-1],1);\n'
            '        }\n'
            '        play.queue=[]; play.last=null;\n'
            '        play.follow=0; play.pushed=0;\n'
            '      }\n'
            '      /* Which frame, and which run it is a frame of.  A count on\n'
            '         its own is a plausible index into any path: after a scan\n'
            '         the page is still standing on the scan\'s trajectory, and\n'
            '         a minimisation that starts a moment later is handed a\n'
            '         number that names a point of somebody else\'s walk.  The\n'
            '         kernel can only honour "the frame on screen" if it knows\n'
            '         the frame is one of the ones it is holding, so the run\n'
            '         goes with the number and it checks. */\n'
            '      send(held?"gfngrab":"gfnfree",\n'
            '           held?(String(play.shown||0)+","+(play.run||0)):"");\n'
            '    }\n'
            '    if(play.held&&!followIsOn()){\n'
            '      /* Still read.  The widget is one slot, so a write that\n'
            '         is not read is a write the next one erases: measured\n'
            '         with Relax off, a drag over six frames lost all six\n'
            '         and the picture jumped forward on release.  Reading\n'
            '         keeps the count honest -- and an abandoned run is\n'
            '         announced on this channel, so a player that stops\n'
            '         reading the moment a hand lands never hears that the\n'
            '         frames it is holding are about a structure that no\n'
            '         longer exists.  The wall is read with it. */\n'
            '      readWall(); read(now);\n'
            '      window.requestAnimationFrame(frame); return;\n'
            '    }\n'
            '    if(play.held){\n'
            '      /* Following: the geometry goes over while the mouse is\n'
            '         still down, and the pace is set by the machine rather\n'
            '         than by a clock.  The next one goes as soon as the last\n'
            '         answer has landed -- GFN-FF answers a small molecule in\n'
            '         under twenty milliseconds and a fixed fifth of a second\n'
            '         threw nine tenths of that away.  A floor keeps the\n'
            '         messages from becoming the bottleneck, and a ceiling\n'
            '         starts again if an answer never comes at all. */\n'
            '      /* The floor is the machine\'s own answer time, not a\n'
            '         constant: never ask more than twice as often as it can\n'
            '         answer, and never faster than an animation frame can show.\n'
            '         On a small molecule that is about twenty milliseconds; on\n'
            '         a hundred atoms it settles at sixty by itself. */\n'
            '      var floor=Math.max(16,Math.min(120,(play.gap||60)/2));\n'
            '      var since=now-(play.pushed||0);\n'
            '      var answered=play.seen>(play.pushedAt||0);\n'
            '      if(since>500||(answered&&since>floor)){\n'
            '        play.pushed=now; play.pushedAt=play.seen;\n'
            '        var api=window.__delfinSubmitManip;\n'
            '        if(!api||!api.pushXyz){\n'
            '          /* An editor from before the follow existed.  Swallowed,\n'
            '             this is a drag that does nothing and says nothing --\n'
            '             which is indistinguishable from a broken kernel. */\n'
            '          if(!play.toldNoPush){ play.toldNoPush=1;\n'
            '            say("this page has an editor without pushXyz; '
            'reload the page"); }\n'
            '        } else {\n'
            '          try{ api.pushXyz(scope,"drag-follow"); }\n'
            '          catch(e){ if(!play.toldNoPush){ play.toldNoPush=1;\n'
            '            say("pushXyz failed: "+e); } }\n'
            '        }\n'
            '      }\n'
            '    }\n'
            '    if(play.queue.length&&!play.follow&&!play.complete\n'
            '       &&!switchIsOn()){\n'
            '      /* Optimise going up abandons what it had computed but\n'
            '         nobody had seen.  A follow has no Optimise behind it --\n'
            '         checking that switch there threw away every frame the\n'
            '         drag produced, which is what "it does nothing" was.\n'
            '         And a run that finished is not a run that was stopped.\n'
            '         The switch goes up by itself the moment the last round\n'
            '         is in, so a path still being walked at a pace anyone can\n'
            '         watch was thrown away exactly then: the picture showed\n'
            '         the first few frames and stopped, and never arrived at\n'
            '         the minimum it had just computed.  The kernel says which\n'
            '         it is -- the last write of a finished run is marked, and\n'
            '         a stop arrives as a halt, which clears the queue above. */\n'
            '      play.queue=[];\n'
            '      if(!play.toldStop){ play.toldStop=1;\n'
            '        say("stopped at frame "+(play.shown||0)\n'
            '            +" of run "+(play.run||0)); }\n'
            '    }\n'
            '    readWall();\n'
            '    read(now);\n'
            '    if(play.queue.length){\n'
            '      if(!play.started) play.started=now;\n'
            '      var ms=stepMs();\n'
            '      if(ms<=0){\n'
            '        /* Everything that has arrived, this frame.  The picture\n'
            '           is then never behind what has been computed, which is\n'
            '           what the top of the slider promises and what a drag\n'
            '           needs whatever the slider says. */\n'
            '        var was=play.last;\n'
            '        play.shown=(play.shown||0)+play.queue.length;\n'
            '        play.last=play.queue[play.queue.length-1];\n'
            '        play.queue=[];\n'
            '        show(was,play.last,1);\n'
            '        play.started=now;\n'
            '        window.requestAnimationFrame(frame); return;\n'
            '      }\n'
            '      var t=(now-play.started)/ms;\n'
            '      if(t>=1){\n'
            '        /* However many steps are due, and the clock moves on by\n'
            '           exactly those -- not back to now.  This loop runs on\n'
            '           the animation frame, so one step per frame with the\n'
            '           clock reset each time snapped the rate to whole\n'
            '           divisors of the screen: 62.7 a second, then 31.3, and\n'
            '           nothing in between.  Every setting from 32 to 59 drew\n'
            '           31.3 -- measured -- so most of the slider did nothing.\n'
            '           Carrying the remainder makes the average the rate that\n'
            '           was actually asked for. */\n'
            '        var due=Math.floor(t);\n'
            '        var prev=play.last;\n'
            '        var left=due;\n'
            '        while(left-->0&&play.queue.length){\n'
            '          play.last=play.queue.shift();\n'
            '          play.shown=(play.shown||0)+1;\n'
            '        }\n'
            '        show(prev,play.last,1);\n'
            '        play.started+=due*ms;\n'
            '      } else if(play.last){\n'
            '        show(play.last,play.queue[0],t);\n'
            '      } else {\n'
            '        /* Counted, like every other frame that is drawn.  This\n'
            '           branch is the first frame of a run -- nothing to\n'
            '           interpolate from -- and it drew without counting, so\n'
            '           the number was one short for the whole of the run\n'
            '           after it.  Measured over ten frames written one at a\n'
            '           time: the picture stood on the tenth and the page\n'
            '           reported nine, and the kernel reads this as a count\n'
            '           (walked[shown-1]) -- so a grab kept the frame before\n'
            '           the one on screen, and a grab on the very first frame\n'
            '           reported zero, which the kernel reads as "no frame at\n'
            '           all" and answers with wherever the calculation had got\n'
            '           to instead. */\n'
            '        play.last=play.queue.shift(); show(null,play.last,1);\n'
            '        play.shown=(play.shown||0)+1;\n'
            '        play.started=now;\n'
            '      }\n'
            '    }\n'
            '    window.requestAnimationFrame(frame);\n'
            '  }\n'
            '  window.requestAnimationFrame(frame);\n'
            '})();'
        )

    def _emit_watcher_js(script):
        """Send the player with the start-up scripts, not through run_js.

        run_js writes into a single Output and clears it first, so a script
        sent at click time can be replaced before the page has run it -- which
        is how the player came to be missing while everything it depends on
        was working.  add_init_js is the channel the explorer's own JS arrives
        on, and that one has never been in doubt.
        """
        try:
            ctx.add_init_js(script)
            return
        except Exception:
            pass
        _run_manip_js(script)

    #: A follow step is a whole xtb process, so it is a few cycles rather than
    #: a minimisation.  Measured on a 102-atom complex: one cycle 0.06 s, five
    #: 0.09 s, ten 0.12 s -- five is about ten answers a second, which reads as
    #: the molecule following the hand rather than catching up with it.
    _GFN_FOLLOW_CYCLES = 5

    #: How far the answer may put a held atom from the cursor before the drag
    #: counts as under-determined.
    #:
    #: A held value is an internal coordinate: xtb meets it and may place the
    #: molecule anywhere that does.  Measured on a backside attack, chloride
    #: driven at a carbon, the answer sits 0.008 A off at the start and 0.102 A
    #: through the transition region -- ordinary freedom, and the price is
    #: about the structure on screen.  On a palladium pushed at head on it was
    #: 0.7 A, and there the price belonged to a geometry where the metal had
    #: got out of the way while the picture showed a bromide 1.27 A from it.
    #: That is worth saying rather than leaving to be discovered.
    _SLIP_LOOSE = 0.25

    #: And with the thermal budget on, more of them.
    #:
    #: The follow then holds the contacts the drag has changed and relaxes
    #: everything else, which is a relaxed scan and answers a harder question
    #: than "settle a little towards the hand": what is the lowest energy this
    #: structure can have with the atoms where the cursor has put them.
    #: Reaching it takes real reorganisation, and cut short the answer is
    #: simply too high -- measured on a Diels-Alder approach at 2.18 A, where
    #: the relaxed path costs +3.6 kcal/mol:
    #:
    #:     cycles       3      6     20
    #:     kcal/mol  +40.0  +18.9   +8.3
    #:
    #: Every one of those is a wall in the wrong place, so this is generous.
    #:
    #: Generous and not wasteful, which was worth checking on a real system
    #: rather than assuming either way.  ``--cycles`` is a *cap*: xtb stops
    #: when the relaxation converges, so a cap above what convergence needs
    #: costs nothing at all.  Measured on a 57-atom manganese complex under
    #: GFN2 at charge +1, one answer priced from the same geometry:
    #:
    #:     cycles      2      4      6      8     12   16   20   30   50
    #:     vs 200  +10.8   +2.3   +0.4  +0.17  +0.01  0.0  0.0  0.0  0.0
    #:
    #: Sixteen is already the converged answer to every digit, and twenty
    #: returns at convergence rather than at the cap.  So the four-fold cost
    #: over a plain drag is the price of relaxing properly, not of counting to
    #: twenty, and cutting the number would buy nothing on an easy answer and
    #: a wall in the wrong place on a hard one.
    #: It also does not start from nothing: the drag walks in small steps and
    #: each answer starts from the structure the last one relaxed, so the work
    #: is spread over the drag rather than repeated in it.  The cost is that
    #: an answer takes a second or two rather than a tenth, and the leash --
    #: which only ever lets the hand go as far as the last answer confirmed --
    #: turns that into a slower drag.  That is the honest price of walking a
    #: path instead of jumping along it, and it was asked for.
    _THERMAL_FOLLOW_CYCLES = 20

    def _begin_gfn_follow():
        """A drag has started and the molecule is to follow it."""
        if not (submit_relax_btn.value
                and _server_method()):
            return False
        # The method on screen is the method that runs.  It used to be GFN-FF
        # whatever the box said, which is a picture of a calculation nobody
        # asked for -- and indistinguishable, from the outside, from the right
        # one.
        state['gfn_follow_method'] = str(submit_ff_dd.value)
        if state.get('gfn_follow'):
            return True     # already following; it keeps the run it began
        # Whatever the last drag left standing, this one starts clean.
        #
        # The release clears it, and the release is one message from the page
        # -- sent from an animation frame when the pointer state changes.  A
        # drag that ends any other way never sends it: the tab is put in the
        # background and the frame loop stops, the pointer is cancelled, the
        # page is reloaded under a hand that is still down.  What was left
        # standing then was the *running maximum* of the last drag, and a
        # maximum is what the wall refuses on -- so the next grab was refused
        # for a crossing it had nothing to do with, and nothing anywhere said
        # so.  That is one drag in a session behaving differently from the
        # rest for no visible reason.
        #
        # Clearing it here as well makes it not depend on the message arriving:
        # a grab is the one moment there certainly is one.
        _clear_thermal_wall()
        # And the floor goes back, which the clearing takes away with the
        # rest.  Nothing has been crossed on the way to a drag that has not
        # started, so this drag's maximum is zero -- and zero rather than
        # absent, because absent lets the first answer set the maximum to
        # whatever it costs.  Below the anchor that is a negative maximum, and
        # a drag that dipped a little before it climbed would then be allowed
        # to climb that much further than the temperature pays for.  The
        # anchor is the zero; falling below it earns no credit.
        state['thermal_peak'] = 0.0
        # And the frontier gap this drag starts from, which is what says
        # later whether the bond order it is reporting can still be believed.
        # A drag that begins at a closed gap is a different statement from one
        # that closes it, so the baseline belongs to the grab.
        state.pop('gfn_follow_gap0', None)
        _gfn_new_generation()
        # What the molecule looked like before this drag: the bonding is read
        # from here, not from a frame that has already been pulled about.
        state['gfn_topology_source'] = _current_xyz()
        # And it is what the first answer measures the hand against.
        #
        # Every later answer sets this to the geometry it handed back, and the
        # first one used to have nothing -- but nothing is not neutral here.
        # With no geometry to compare against, contacts_holding cannot see
        # what the hand has changed and falls back to the nearest contact,
        # which for a grabbed atom is usually a bond it is not driving; and
        # as_pushes then asks for that bond at the length it already has,
        # which is no force at all.  So the first answer of every drag was a
        # free relaxation with the coordinate the hand *is* driving left
        # loose.
        #
        # Measured with Dynamik Opt on and Auto off, an ethene carbon dragged
        # at a diene terminus, six mouse moves each way: the first answer of
        # each drag was told to hold the ethene's own C=C at the length it
        # already stood at, 1.373 A and then 1.387.  On the first grab that costs a step -- 0.03 A where
        # the answers after it move 0.10 to 0.16.  On a second grab it is the
        # spring back the report is about: the hand had left the contact at
        # 2.772 A, the first answer relaxed it back to 3.033, and three more
        # answers went on recovering ground the user had already covered.
        # Same gesture in two halves did not equal the gesture in one.
        #
        # The structure on screen is where the hand arrived, so that is what
        # the hand is measured against, and the drag bites from its first
        # answer.
        state['thermal_was'] = state['gfn_topology_source']
        # The one thing that history takes away, put back below: with no
        # geometry to compare against, contacts_holding decides between a turn
        # and a stretch by a rule about bonds, and with one it decides by which
        # coordinate moved most.  Both are wanted -- see where this is read.
        state['gfn_follow_opening'] = True
        state['gfn_follow'] = True
        state['gfn_follow_steps'] = 0
        state['gfn_follow_frames'] = []
        # Whatever the page said about the last drag's playback is not about
        # this one.
        state['gfn_play_note'] = ''
        run = _note_the_run(int(state.get('gfn_run', 0)) + 1, 'follow')
        state['gfn_run'] = run
        state['gfn_follow_run'] = run
        return True

    def _clear_thermal_wall():
        """A drag that has ended has no hand to hold back.

        Sent every time, not only when the budget had actually run out.  A
        leash is armed the moment an atom is picked up and shortens as the
        budget is spent, so a drag that stayed well inside it still leaves
        marks and a reach on the page -- and only the walled case was being
        released.  The next drag then pulled against the last one's marks,
        which is a molecule that will not move properly and no line anywhere
        saying why.  Switching the budget off did not release it either, so
        it outlived the feature that made it.
        """
        state['thermal_safe'] = {}
        state.pop('thermal_walled', None)
        # The geometry a turn is measured against belongs to the drag that
        # is over.  Kept, the first answer of the next one would compare
        # against a structure from before whatever happened in between.
        state.pop('thermal_was', None)
        state.pop('thermal_turn', None)
        state.pop('thermal_holding', None)
        # The slope too.  Nothing cleared it, so the first answer of every
        # drag reported the last one's number -- on a different atom, and
        # after a rollback with the opposite sign: measured, an atom being
        # pulled strictly uphill was reported as "Falling 94 kcal/mol per A".
        state.pop('thermal_slope', None)
        state.pop('thermal_last', None)
        # And what the last answer of the last drag was priced at, which is
        # what "one answer cost this much" is measured against.
        state.pop('thermal_priced', None)
        # And the highest price this drag went through, with the one recorded
        # for the geometry it could fall back to.  A ceiling is a barrier
        # height, so the wall refuses on the maximum since the anchor -- and
        # that maximum belongs to the drag that is over.  What a release
        # leaves behind is a structure the budget agreed to, whose own way
        # here was affordable the whole length of it, so the next grab starts
        # from where it is standing rather than from the last one's worst
        # moment.  Kept, one bad drag would have refused every drag after it
        # until the anchor was set again.
        state.pop('thermal_peak', None)
        state.pop('thermal_good_peak', None)
        state.pop('thermal_over', None)
        # And what the hand was asking for when it ran out, which is the one
        # thing here that decides whether the *next* answer runs at all.  It
        # outlived the drag that recorded it, so a fresh grab could be met by
        # the last one's demand and stand still for it -- and now that the
        # question is asked before a worker is begun rather than inside one,
        # standing still means not starting.  The first report of a new drag
        # asks for the coordinate where it already stands, so this cleared
        # itself on the way past; it should not depend on that.
        state.pop('thermal_spent', None)
        _push_thermal_wall(None)

    def _end_gfn_follow():
        state['gfn_follow'] = False
        # The hand has gone, so there is nothing left to hold back.  Left
        # standing, the wall would meet the next drag with the positions of
        # the last one -- atoms that are not being held any more, marked at
        # places the structure has since moved away from.
        _clear_thermal_wall()

    #: What the page says about the playback that the user could not see for
    #: themselves.  "received 41 frames" and "drawing" are the playback
    #: working; a trajectory that is not drawing is not.
    _PLAY_NOTE_FAULTS = ('setPositions did not draw', 'no setPositions on the page')

    def _trajectory_fault():
        """What the page could not do with the frames, or ``''``.

        Three things report while a structure is being worked on: the kernel
        counts the follow steps it has run, the optimisation says what it
        reached, and the page counts the frames it has drawn.  They wrote the
        line in turn and in different shapes -- one row, then two, then one
        again, several times a second -- and the status line stands above the
        viewer in the column, so the picture stepped up and down with it while
        the user was trying to aim an atom.

        The page's count is not a row.  It was there to tell an invisible
        trajectory from a missing one, and "received 41 frames" is the case
        where nothing is wrong; a fault goes on the end of the line that is
        already there rather than under it.

        On the end of it at the moment it is drawn, which is the half that had
        to change.  This used to be composed into a message and written, and
        the message it was composed into was whatever the *optimiser* had last
        said -- so every report the page sent rewrote the row with a sentence
        about a run that might have finished minutes ago.  Measured: one press
        of Optimize, then a scan of the central C-C of butane under GFN2 in
        benzene, then the one message the browser sends when it takes the
        scan's frames -- and the row went from "The scan walked 8 of 8 points.
        Highest +52.8 kcal/mol at 2.2 ..." back to "Optimised with GFN2-xTB.
        E = -13.670432 Eh ...", which is the line the user was left looking at
        instead of the answer the scan had spent its minutes on.  Nothing was
        running: there was no optimisation behind that sentence at all.

        A decoration rather than a message, the way the spinner is: what the
        row says belongs to whoever wrote it, and the page adds to it without
        being able to take it away.  It cannot pile up either -- it is applied
        once per draw rather than folded into what is remembered.
        """
        note = str(state.get('gfn_play_note') or '')
        if note and any(fault in note for fault in _PLAY_NOTE_FAULTS):
            return f'Trajectory: {note}.'
        return ''

    # Whether there is a calculation behind what the line is reporting used to
    # be asked here, for the playback note to carry a spinner with it.  It is
    # asked of the register of running workers instead -- see _busy_now, which
    # is the only authority on the ring -- and the note no longer writes a
    # message that would need one.

    def _hand_pulls():
        """Whether the hand on the molecule is a force or a placement.

        One question, asked in one place, because several things in this
        toolbar mean nothing without a force behind the hand and every one of
        them used to ask for itself.  Both hands are wanted and neither is the
        lesser: placing an atom exactly where it is meant to go is what
        building a structure *is*, and no force can do it, because the whole
        point of a force is that the chemistry gets a say.  What differs is
        which controls have anything to act on -- see
        :func:`_refresh_hand_controls` for the list and for why each one is on
        it or is not.
        """
        return str(submit_hand_dd.value) == 'pull'

    def _pull_force():
        """How hard the hand may pull under a server method, in kcal/mol/A.

        ``None`` is the rigid hand: the coordinate is set and the calculation
        is told to meet it, which is what this always did.

        A share of what a bond holds, so the setting means the same thing here
        as it does in the browser's field even though the two apply it quite
        differently.  Measured under GFN2 by pushing until the coordinate
        either settles or runs away, a C-C gives at 56 kcal/mol/A, a C-H at
        49, a C-O at 60 -- near enough one number that "half a bond" is a
        statement about any molecule rather than about ethane.

        The temperature does not come into it, and this used to say that it
        did.  It was written when the budget derived the hand from the ceiling
        -- 22.3 kcal/mol over a reach was a hand of 45 kcal/mol per Angstrom
        -- and that chain was taken out because it needs a length no
        temperature supplies and it forbade what the temperature allows: sized
        as a distance the hand was too weak to turn a torsion, so a molecule
        could not be put into its own conformers at room temperature.  See
        :func:`_pull_most`, which is what is left of it.

        So the slider is the whole of the hand, with the budget on or off, and
        what the temperature limits is the *energy* of what is kept.  That is
        the wall's job: a hand strong enough to drive one barrier can always
        be dragged for longer than one barrier's worth, and no force ceiling
        anywhere can prevent that.  The two are different jobs and this is
        only the first of them.
        """
        if not _hand_pulls():
            # The rigid hand: the coordinate is set and whoever is calculating
            # is told to meet it.  Not a force at all, which is what None says.
            return None
        share = float(submit_pull_slider.value or 0.0)
        if share <= 0:
            return None
        return max(0.5, share * _gfn.A_BOND_HOLDS)

    def _still_spent(current, holding):
        """Whether the budget is spent and the hand has not come back in.

        Remembered from the step that spent it, as the coordinates the hand
        was asking for then.  The hand has come back when every one of them is
        closer to what the kept structure has than it was -- which is the hand
        easing off, and is the only thing that can make the drag affordable
        again.
        """
        stuck = state.get('thermal_spent')
        if not stuck:
            return False
        good = state.get('thermal_good')
        if not good or len(_gfn.atom_lines(good)) != len(
                _gfn.atom_lines(current)):
            state.pop('thermal_spent', None)
            return False
        for one in stuck:
            wanted = _value_in(current, one)
            standing = _value_in(good, one)
            if wanted is None or standing is None:
                continue
            if abs(wanted - standing) < abs(one['value'] - standing) - 1e-6:
                state.pop('thermal_spent', None)
                return False
        return True

    def _hand_share():
        """What the page should make of the slider, given which hand it is.

        Zero is the rigid hand there as it is here, so the two sides never
        disagree about which one is in the user's hand.
        """
        if not _hand_pulls():
            return 0.0
        return float(submit_pull_slider.value or 0.0)

    def _pull_most():
        """The hardest the hand may ever pull.  Nothing, now, ever.

        A temperature grants an *energy*, and that is the whole of what it
        says.  Steepness follows from it on its own: a steep part of the
        surface spends the budget in a short distance, so at a low temperature
        a steep path simply does not go far -- there is nothing left for a
        force ceiling to add.

        It was capping the force as well, derived from the ceiling over a
        chosen length.  That needs a length no temperature supplies, and it
        can forbid what the temperature allows: sized as a distance it left
        the hand too weak to turn a torsion, so a molecule could not be put
        into its own conformers at room temperature -- which is the one thing
        room temperature certainly does allow.

        What enforces the temperature is the wall, which prices the structure
        that was actually reached and hands back the last one that was inside
        the budget.  That is exact, it is coordinate-independent, and it needs
        no length at all.  The slider is the hand, the wall is the limit, and
        the two do not have to agree about anything.

        Kept as a function rather than deleted because the question "should
        anything cap this" is worth having an answer to in one place.
        """
        return None

    def _gfn_follow_step(xyz, holding=()):
        """Relax around the atom the hand is holding, and send that back.

        The dragged atom is not held by xtb: it cannot be.  ``$fix atoms:`` is
        broken in xtb 6.7.1 and ``$constrain atoms:`` naming one atom does
        nothing at all, both measured -- xtb holds internal coordinates, not
        positions.  The hold is the page's instead: the frames that come back
        are written to every atom *except* the ones being dragged, so the
        cursor keeps the one it is holding and xtb arranges the rest around it.

        Only one process at a time, and the newest geometry wins: a hand moves
        faster than xtb answers, and a queue of answers about where the atom
        used to be is worse than no answer at all.

        Where an answer's time goes, measured on a 57-atom manganese complex
        under GFN2 at charge +1 -- the system a user was dragging when he
        reported that reaching the thermal ceiling took ten minutes:

            xtb process        99.6 %
            everything else     0.4 %

        The "everything else" is all of this function's arithmetic together:
        the contact perception, turning the held values into pushes, the
        restraint energy, the rigid-body fit and the closest-contact scan.
        Under GFN-FF, where the process is a tenth of a second rather than
        three seconds, it is still 89 % xtb and about 8 % reading the
        trajectory back out of the log.  So there is nothing here worth
        optimising, and what an answer costs is what xtb costs.
        """
        if not state.get('gfn_follow') and not _begin_gfn_follow():
            return          # not following: Relax is up, or GFN is not chosen
        # A drag the budget is already refusing is not started again.
        #
        # The standing-still rule lived inside the worker, so standing still
        # cost a whole worker: the thread was begun, the ring was turned on,
        # the rule said there was nothing to compute, and the thread ended --
        # once per report from the page, and the page reports at the rate a
        # mouse moves.  Measured on the user's manganese complex at 298 K with
        # the phenolate oxygen held against the wall, thirty seconds of hand
        # on the structure:
        #
        #                          budget off   budget on, refusing
        #     reports from the page      1147                  1625
        #     workers started               1                  1277
        #     status lines written          9                  2555
        #     frames drawn                  8                     1
        #
        # Forty-three threads a second, and the line that lies over the
        # picture rewritten eighty-five times a second -- alternating on the
        # ring alone, since the words do not change while nothing is being
        # computed.  The picture is still and everything around it is not,
        # which is what the shaking is.  With the budget off the same hand
        # costs one worker, because xtb is genuinely busy and every report
        # after the first is folded into gfn_follow_xyz.
        #
        # Asked here it costs no thread at all: it is coordinate arithmetic
        # over the atoms the hand is on.  The worker keeps its own copy of the
        # question, for a report that arrived while xtb was running.
        if _still_spent(xyz, holding):
            return
        if state.get('optimize_run') is not None:
            # An optimisation owns the structure, so the follow does not get
            # to move it as well.  Picking an atom up interrupts a run before
            # the follow begins, so the ordinary way round never reaches this;
            # pressing Optimise while an atom is already held does, and two
            # xtb processes then arranged the same molecule around each other
            # -- the same collision the coordinate box was given an owner for,
            # one step earlier in the same path.  The press is the later of
            # the two things the user asked for, so it wins.
            return
        state['gfn_follow_xyz'] = (xyz, tuple(holding or ()))
        if state.get('gfn_follow_busy'):
            return
        state['gfn_follow_busy'] = True
        method = str(state.get('gfn_follow_method') or submit_ff_dd.value)
        label = _server_label(method)
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        constraints = list(state.get('constraints') or [])
        wet = str(submit_gfn_solvent.value or '') or None
        model = _solv_model()

        def _work():
            try:
                while state.get('gfn_follow'):
                    newest = state.pop('gfn_follow_xyz', None)
                    if newest is None:
                        return
                    current, holding = newest
                    began = time.perf_counter()
                    # Once the budget is spent, stand still rather than shake.
                    #
                    # The wall acts and then undoes: the step runs, the price
                    # arrives, and a geometry past the budget is replaced by
                    # the last one that was inside it.  Under a hand that goes
                    # on pulling, that is strain-and-snap several times a
                    # second -- the structure springs back to the same place
                    # over and over, which is what the shaking was.  It is
                    # spent; there is nothing left to compute.  So the answer
                    # is handed straight back until the hand comes in again,
                    # which costs no xtb at all and is perfectly still.
                    # And the same for the bonding, for the same reason.
                    #
                    # A rigid hand asking for a place that breaks a bond and a
                    # wall that will not have it are two answers to one
                    # question, and they were both being drawn: the page moved
                    # the atom, the answer came back and put it home, and the
                    # next mouse move did it again.  On screen that is a
                    # trajectory, a spring back, and another trajectory -- the
                    # motion twice over, several times a second.
                    #
                    # Once it is clearly refusing, it stands still: no xtb, no
                    # frames, and the line already says which two settings are
                    # disagreeing.
                    if (int(state.get('topology_refused') or 0) >= 3
                            and submit_topology_btn.value):
                        schedule_ui_update(
                            _set_mol_status,
                            state.get('gfn_last_status') or '',
                            spinner=True)
                        continue
                    if _still_spent(current, holding):
                        schedule_ui_update(
                            _set_mol_status,
                            state.get('gfn_last_status') or '',
                            spinner=True)
                        continue
                    # With the budget on, the drag is followed differently:
                    # the contacts the hand has changed are held and the rest
                    # is relaxed around them.  That one calculation is both
                    # the follow and the price -- see contacts_holding for why
                    # the geometry the mouse leaves cannot be priced as it
                    # stands, and why holding *these* values and no others is
                    # what makes the difference between a reaction and a
                    # structure being torn open.
                    #
                    # Empty when the hand is on every atom or on none: a whole
                    # structure moved as one has changed no internal
                    # coordinate and costs nothing, and there is nothing to
                    # hold it against.  The old single point stands in.
                    came_back = None
                    # The page names the atoms it is holding by their place in
                    # the structure it is looking at.  When not one of those
                    # numbers names an atom that is here -- the structure
                    # changed under the drag, or a frame arrived from before it
                    # did -- there is no coordinate to hold, and the price fell
                    # straight through to the rigid single point.  Which is the
                    # one that overcharges anything forming a bond twenty to
                    # thirty times over, and it did so with nothing said: the
                    # wall could roll a drag back on a number computed for a
                    # structure neither side was looking at.  So it is not
                    # priced at all, and the line says why.
                    count = len(_gfn.atom_lines(current))
                    stale = bool(holding) and not any(
                        0 <= int(i) < count for i in holding)
                    pricing = _thermal_live() and not stale
                    # The hand as a force rather than as a value.
                    #
                    # The coordinates the hand has changed were held exactly,
                    # so the calculation was told what the answer had to be and
                    # arranged the rest around it.  There was nothing the
                    # chemistry could refuse: a bond went wherever the mouse
                    # went, and the budget could only say afterwards that it
                    # should not have.  Pushed instead, the coordinate gets a
                    # force with a ceiling and settles where the two balance,
                    # which is what "the atom follows the mouse as far as it
                    # can" means.  See :func:`_pull_force` for the ceiling.
                    # MOPAC takes no held internals at all, so there is no
                    # coordinate for a force to act on there: the hand stays
                    # what it was, which is a placement.  Left in, a pull under
                    # PM7 would have had nothing to push and everything to
                    # undo -- the answer owns the atom, and without a push the
                    # answer is a free relaxation.
                    pull = (_pull_force()
                            if not stale and not _mopac.is_mopac_method(method)
                            else None)
                    # xtb takes one force constant for the whole block, and a
                    # push is three orders of magnitude softer than a hold.
                    # With both in it the hold's stiffness wins, so the pull
                    # would quietly go back to being a placement -- the same
                    # number, the same picture, a different experiment.  The
                    # user's hold is what they asked for, so it stands, and
                    # the line says what that costs.
                    held_too = bool(constraints) and pull is not None
                    if held_too:
                        pull = None
                    # Worked out whenever there is a hand to describe, and not
                    # only when there is a budget to charge.  Behind the
                    # budget, a pull with the budget off had no coordinate to
                    # act on at all: the follow relaxed freely and the answer
                    # -- which no longer puts the atom back under the cursor --
                    # undid the drag on every step.
                    contacts = (
                        _gfn.contacts_holding(
                            current, holding, most=3,
                            was=state.get('thermal_was'),
                            turning=state.get('thermal_turn'),
                            holding=state.get('thermal_holding'))
                        if ((pricing or pull is not None)
                            and not _mopac.is_mopac_method(method)) else [])
                    if contacts and state.pop('gfn_follow_opening', False):
                        # The first answer of a drag decides whether the hand
                        # is turning something or driving a contact, and the
                        # decision sticks -- see thermal_turn below, and the
                        # cyclohexane that fell back towards the chair every
                        # time it changed.  Asked by what has moved most the
                        # decision goes to the bond the grabbed atom hangs on,
                        # because on the first answer almost nothing else has
                        # moved at all; and driven, a bond is the one
                        # coordinate a hand must not have -- measured on a
                        # 2,4-hexadiene, a chain carbon dragged 1.75 A moved
                        # 0.09 A under the pull and tore three bonds under the
                        # rigid hand.
                        #
                        # So the question is put once more, as what each
                        # coordinate *can* do rather than what it has already
                        # done, and a turn that carries the hand wins.  It
                        # costs no calculation -- the same geometry, read
                        # again -- and it is asked once per drag.
                        #
                        # With the geometry the drag started from, which it
                        # used not to have.  Blind to which way the hand went,
                        # the turn it chose was whichever the walk came to
                        # first, and a torsion carries the grabbed atom in one
                        # direction only: on a chelate, where every torsion at
                        # a grabbed atom swings it out of a ring plane the
                        # hand is not pulling out of, that is a coordinate the
                        # drag cannot express.  A contact is still left to the
                        # scored answer, which is the half of this that
                        # history is better at.
                        opening = _gfn.contacts_holding(
                            current, holding, most=3,
                            was=state.get('thermal_was'), opening=True)
                        if opening and str(
                                opening[0].get('kind')) == 'dihedral':
                            contacts = opening
                    # Keep bonds, the way GOAT keeps them: frozen while the
                    # structure is being pushed, rather than the step refused
                    # afterwards.  The way not to break a bond is not to let
                    # it move -- measured on a 2,4-hexadiene, refusing alone
                    # still let one go, because a bond stretches across
                    # several accepted steps and what is handed back is the
                    # last of those.
                    #
                    # Only under the rigid hand.  A push needs a soft force
                    # constant and xtb takes one for the whole block, so with
                    # a pull these would make the hand a hold; a pull does not
                    # tear anything anyway, and the wall stays the backstop
                    # for it.
                    keeping = list(constraints)
                    if (submit_topology_btn.value and pull is None
                            and not _mopac.is_mopac_method(method)):
                        walking = {tuple(sorted(one['atoms']))
                                   for one in contacts
                                   if len(one.get('atoms') or ()) == 2}
                        keeping += [
                            one for one in _gfn.bonds_to_freeze(current)
                            if tuple(sorted(one['atoms'])) not in walking]
                    if pull is not None and contacts:
                        contacts = _gfn.as_pushes(
                            contacts, state.get('thermal_was') or current,
                            pull, value_of=_value_in, most=_pull_most())
                    if _mopac.is_mopac_method(method):
                        # MOPAC takes no held internals and no topology file,
                        # so it is given what it does take. A few cycles, the
                        # same as the xtb side: enough to move towards the
                        # hand, not so many that the answer is stale when it
                        # arrives. Measured on a benzophenone at 105 ms a run
                        # against GFN-FF's 36, and COSMO takes that to 155 for
                        # PM7 and 175 for PM6-D3H4 -- half again, and still
                        # inside the budget for a drag.
                        outcome = _mopac.optimize_with_mopac(
                            current, method, charge=charge, uhf=uhf,
                            max_steps=_GFN_FOLLOW_CYCLES, timeout=30.0,
                            constraints=constraints, solvent=wet)
                    else:
                        outcome = _gfn.relax_steps(
                            current, method=method, charge=charge, uhf=uhf,
                            # Twenty cycles are the *budget's*, not the
                            # hold's: a price has to be for a properly
                            # relaxed path or the wall stands in the wrong
                            # place.  Keyed on the contacts, a pull with no
                            # budget quadrupled the cost of every step of
                            # every drag -- 0.09 s became 0.36 on a hundred
                            # atoms -- to buy an accuracy nothing was going
                            # to read.  A drag that is only a drag stays at
                            # five, which is about ten answers a second.
                            cycles=(_THERMAL_FOLLOW_CYCLES if pricing
                                    else _GFN_FOLLOW_CYCLES),
                            timeout=30.0,
                            constraints=keeping + contacts, solvent=wet,
                            solvation_model=model,
                            topology=_gfn_topology_dir(current),
                        )
                    if not outcome.get('ok'):
                        note = str(outcome.get('status') or 'it did not run')
                        schedule_ui_update(
                            _set_mol_status,
                            f'The molecule stopped following: {note}')
                        return
                    # What this answer computed and used to throw away with
                    # its scratch directory. A list assignment and, with the
                    # labels off, nothing else at all -- see _repaint_labels.
                    _remember_charges(outcome)
                    if submit_labels_btn.value:
                        # Rate-limited: see _repaint_labels for why a drag
                        # must not put a script on the run_js channel once per
                        # answer, and why labels are the one thing that can be
                        # dropped there without being wrong.
                        schedule_ui_update(_repaint_labels, False)
                    # Say it out loud, every step.  A follow that is working
                    # and a follow that is not both look like a molecule that
                    # is not moving much, and the difference is the whole
                    # question when something in the chain is broken.
                    steps = int(state.get('gfn_follow_steps') or 0) + 1
                    state['gfn_follow_steps'] = steps
                    # How many atoms the hand is on.  Grabbing an atom that is
                    # part of a selection drags the whole selection, so a
                    # selection left over from earlier makes every drag move
                    # everything -- which reads as the molecule fighting
                    # itself, and is invisible unless it is counted here.
                    many = (f' holding {len(holding)} atoms,'
                            if len(holding) > 1 else '')
                    # What the hand is allowed to pull with, said in the units
                    # that decide what it can do.  A number nobody can see is
                    # a number nobody can act on, and this one is the whole
                    # difference between a deformation and a torn molecule.
                    hand = ''
                    if pull is not None and contacts:
                        # What is being applied, and the most this hand could
                        # apply -- two numbers, because they are two different
                        # things and only the second one was ever shown.
                        #
                        # A push is a spring: it reaches its ceiling when the
                        # target is a whole reach away and is weaker all the
                        # way there.  So what the hand really pulls with is
                        # set by how far the cursor has run ahead of the
                        # structure, which is the thing the mouse controls
                        # while the drag is running -- and that is exactly
                        # what the user asked whether he could do.  Measured
                        # on a manganese complex, a phenolate oxygen dragged
                        # at 0.4 of a bond: about half a degree of torsion of
                        # lead per answer, which is 0.4 kcal/mol per radian
                        # applied.  The line said 44.
                        driving = contacts[0]
                        kind = str(driving.get('kind') or 'distance')
                        constant = float(driving.get('force') or 0.0)
                        hardest = _gfn.push_pulls_hardest(kind) * constant
                        stands = _value_in(state.get('thermal_was') or current,
                                           driving)
                        applied = hardest
                        if stands is not None and driving.get('value') is not None:
                            applied = _gfn.push_pulls_now(
                                kind, constant,
                                float(driving['value']) - float(stands))
                        units = ('kcal/mol/A' if kind == 'distance'
                                 else 'kcal/mol per radian')
                        hand = (f' Pulling at {applied:.1f} of a possible '
                                f'{hardest:.0f} {units} -- drag further ahead '
                                f'of the atom to pull harder. The most this '
                                f'hand may reach is '
                                f'{hardest / _gfn.A_BOND_HOLDS:.2f} of what a '
                                f'bond holds.')
                    if held_too:
                        hand = (' A held value and a pull cannot share one '
                                'force constant in xtb, so the hand is '
                                'placing rather than pulling.')
                    said = (f'{label} is following the drag:{many} {steps} '
                            f'step(s), '
                            f'{(time.perf_counter() - began) * 1000:.0f} ms '
                            f'each.{hand}')
                    # The bond order of the pair the hand is driving, as a
                    # readout and as nothing else.
                    #
                    # It decides nothing. This was built to be the honest
                    # answer to "is the bond still there", on the grounds that
                    # the editor's own watch -- :func:`gfn_optimize._is_a_bond`,
                    # covalent radii with slack -- is a cliff at 1.94 A, and
                    # the measurement says the opposite: an ethane C-C held at
                    # 3.03 A with everything else relaxed still reads 1.000,
                    # because a single closed-shell determinant keeps that pair
                    # in one orbital however far the two carbons are taken. At
                    # the same geometry the frontier gap has fallen from 15.3
                    # to 0.75 eV and a separate measurement puts the fractional
                    # occupation density at 1.73 electrons. So as a
                    # bond-existence test the order is worse than the distance,
                    # not better: it says everything is fine exactly where the
                    # method has stopped working. The bond watch stays
                    # geometric, and this is a number to read.
                    #
                    # It is worth reading. Where the two fragments stay closed
                    # shells the order tracks the bond all the way -- an
                    # ammonia borane N-B runs 0.61, 0.51, 0.40, 0.24, 0.00 from
                    # 1.66 to 2.86 A -- and 1.9 across a C=C or 0.61 for a
                    # dative bond at its own minimum is something a chemist can
                    # use.
                    #
                    # On the pair the gesture named, not on the atoms: an
                    # order belongs to two atoms and nothing else, and summed
                    # onto one of them it would be a valence, which is a
                    # different quantity answering a different question. It
                    # goes in the line rather than into the picture for the
                    # same reason the price does -- the drag already has the
                    # user's eyes on the structure, and a second row of
                    # numbers over the atoms while an atom is being aimed is
                    # more than can be read.
                    #
                    # Nothing at all under GFN-FF, which computes no bond
                    # order.
                    driven = next((one for one in contacts
                                   if len(one.get('atoms') or ()) == 2), None)
                    if driven is not None:
                        pair = [int(i) for i in driven['atoms']]
                        order = _gfn.bond_order_between(
                            outcome.get('bonds'), pair[0], pair[1])
                        if order is not None:
                            rows = _gfn.atom_lines(outcome.get('xyz') or '')
                            names = 'The pair the hand is driving'
                            if all(0 <= i < len(rows) for i in pair):
                                names = (
                                    f'{rows[pair[0]].split()[0]}{pair[0]}-'
                                    f'{rows[pair[1]].split()[0]}{pair[1]}')
                            # With the gap of the same answer, and the gap
                            # this drag started from. Where the frontier gap
                            # has closed the order is being read off a
                            # description that has stopped working -- an
                            # ethane pulled to 3.03 A still reads 1.00 -- and
                            # the fall from where the drag began is what says
                            # so soonest.
                            note = _gfn.bond_order_note(
                                order, names, outcome.get('gap'),
                                state.get('gfn_follow_gap0'))
                            if note:
                                said = f'{said} {note}'
                    if state.get('gfn_follow_gap0') is None:
                        state['gfn_follow_gap0'] = outcome.get('gap')
                    # What the drag has cost so far: the energy of the
                    # relaxation that has just run, which held the drag and
                    # let everything else move.
                    #
                    # This was a fresh single point on the geometry the page
                    # sent, and before that the energy of the free follow.
                    # Both are wrong, in opposite directions, and the two
                    # mistakes are worth keeping written down because each
                    # looks right on its own.
                    #
                    # The free follow undoes the drag -- xtb cannot be asked
                    # to leave an atom alone, so it pulls it home and reports
                    # the repaired structure.  On a benzene with a hydrogen
                    # dragged out, five free cycles: C-H back at 1.08 A and
                    # -0.8 kcal/mol at 1.58, back at 1.09 and -1.2 at 2.58.
                    # A proton could be pulled off with the line underneath
                    # saying everything was fine.
                    #
                    # The single point overcharges anything that forms a
                    # bond, because a transition state consists of everything
                    # rearranging and a rigid geometry rearranges nothing.
                    # Butadiene and ethylene at 2.18 A: +124.3 against +3.6
                    # along the relaxed path.  So the wall said no to every
                    # reaction as firmly as it said no to tearing a ring
                    # open, and a Diels-Alder that room temperature allows
                    # needed 1498 K to get through.
                    #
                    # Held-and-relaxed is neither.  It cannot undo the drag,
                    # because the drag is what is held; and it charges for
                    # the distortion that is really left once the structure
                    # has done what it can, which is what a barrier is.
                    # A value held exactly leaves no restraint energy worth
                    # the name -- 20.0 meets it to a thousandth -- so that
                    # answer is its own price.  A push leaves a real one, it is
                    # meant to, and priced as it stands the budget would be
                    # charged for the hand as well as for the structure:
                    # measured on an ethane pulled to 2.46 A, 1.42 kcal/mol of
                    # it was the restraint.  So a pushed step is priced by a
                    # single point on the geometry it actually reached, which
                    # is also the geometry on screen.
                    if contacts and pull is None:
                        priced = outcome
                    elif contacts and pricing:
                        # The restraint's own energy is arithmetic, not a
                        # second calculation: xtb's is k*d^2 with d in Bohr,
                        # measured against a single point at 0.0000 kcal/mol
                        # over a distance, an angle and a torsion.  Asked for
                        # instead, it was a whole extra xtb process on every
                        # step of every drag.
                        bias = _gfn.restraint_energy(
                            outcome['xyz'], keeping + contacts, _value_in)
                        priced = (
                            dict(outcome, energy=float(outcome['energy']) - bias)
                            if bias is not None else
                            _gfn.optimize_with_gfn(
                                outcome['xyz'], method, charge=charge, uhf=uhf,
                                timeout=30.0, solvent=wet,
                                solvation_model=model,
                                topology=_gfn_topology_dir(current),
                                optimise=False))
                    elif pricing:
                        priced = _gfn.optimize_with_gfn(
                            current, method, charge=charge, uhf=uhf,
                            timeout=30.0, solvent=wet, solvation_model=model,
                            topology=_gfn_topology_dir(current),
                            optimise=False,
                        )
                    else:
                        priced = {}
                    # What this step would leave behind, worked out before the
                    # budget is asked about it rather than after.
                    #
                    # The price above belongs to the geometry xtb relaxed, and
                    # under a pull that is emphatically not the geometry the
                    # mouse asked for: the push is soft, the bond comes most of
                    # the way back, and the two structures part company further
                    # with every step of the drag.  Priced one and kept the
                    # other, the budget was answering about a molecule nobody
                    # was looking at -- measured on an ethane methyl pulled
                    # out at 298 K, a 22.3 ceiling: the relaxed answer went
                    # +0.1, +1.1, +0.2, +4.0, +6.8, +12.2, +16.2 kcal/mol and
                    # was allowed every time, while what was left in the box
                    # went +15.8, +45.3, +74.8, +99.1, +116.9, +129.0, +136.7.
                    #
                    # So the wall is asked about the structure that survives
                    # the step, and the two placements below are chosen to
                    # keep that honest: settle_onto is a rigid body, so the
                    # energy that was priced is exactly the energy of what is
                    # kept; hold_atoms_at moves single atoms, and how much
                    # that changes the structure is the slip measured further
                    # down, which refuses the step once it is more than a
                    # nudge.
                    #
                    # Under a rigid hand the held atoms go back where the
                    # cursor had them; that is what a rigid hand is.  Under a
                    # pull they emphatically do not -- how far the atom got is
                    # the answer to the question the drag asked, and putting it
                    # back under the cursor is how a force was made to look
                    # like a move.  The structure is laid onto the geometry it
                    # was handed by everything that is *not* held, so the
                    # molecule stays where it is on screen and only what was
                    # pulled has moved.
                    if pull is not None and contacts:
                        count = len(_gfn.atom_lines(current))
                        grabbed = {int(i) for i in (holding or ())}
                        rest = [i for i in range(count) if i not in grabbed]
                        laid = _gfn.settle_onto(
                            outcome['xyz'], current, rest or range(count))
                        reached = laid
                    else:
                        # Twice over, and the first of the two is what makes
                        # the second honest.  A held value is an internal
                        # coordinate: xtb meets it and is then free to put the
                        # whole molecule anywhere that does, so the answer
                        # comes back turned and slid bodily away from the
                        # cursor.  Laid back on as a rigid body that costs
                        # nothing -- an energy does not depend on where a
                        # molecule is -- and what is left over is the only
                        # part that is really a shortfall.
                        #
                        # Measured on a butane turned about its middle bond
                        # under a rigid hand: the body slid 0.27 A, which was
                        # counted as slip and refused the turn at 0.25 -- with
                        # the price standing at +0.0 of 22.3 kcal/mol, so the
                        # temperature was allowing it and the drag stopped
                        # anyway.  Laid on first, the residue is 0.005 A.
                        laid = _gfn.settle_onto(
                            outcome['xyz'], current, holding)
                        reached = _gfn.hold_atoms_at(laid, current, holding)
                    if not stale:
                        state['thermal_now'] = priced.get('energy')
                    if stale and _thermal_live():
                        said = (f'{said} The hand is on atoms this structure '
                                f'does not have, so this step is not priced.')
                    if pricing:
                        # How far the answer put the held atoms from where
                        # the cursor has them.  A held value is an internal
                        # coordinate, so xtb meets it and is free to place the
                        # molecule anywhere that does; while the values really
                        # determine the drag that freedom is small, and when it
                        # is not, the drag is under-determined and the price
                        # belongs to a different structure than the picture.
                        # On a palladium being pushed at, that was 0.7 A and
                        # the only sign that anything was wrong.
                        #
                        # Measured here, beside the geometry it is about, and
                        # said further down where the line is built.  Moving
                        # the whole of it down turned an UnboundLocalError into
                        # a NameError in a different function, which is what
                        # comes of chasing a use instead of the definition.
                        # Only meaningful for a hand that sets a value.  A
                        # push is *meant* to come up short -- that is the
                        # answer it asked for -- so measuring the shortfall
                        # against a rigid hand and refusing the drag for it
                        # would refuse every pull there is.
                        # Between the two placements above, so the molecule's
                        # freedom to sit anywhere is not counted as a
                        # shortfall.  Measured against the raw answer it was,
                        # and a butane turned about its middle bond was refused
                        # for a 0.27 A body slide while the price stood at
                        # +0.0 of 22.3 -- a refusal the temperature had nothing
                        # to do with.
                        slipped = (_gfn.largest_shift(reached, laid)
                                   if (contacts and pull is None) else 0.0)
                        # Two refusals that do not depend on the budget.
                        #
                        # A loose hold means the price is about a nearby
                        # structure rather than this one, so advancing the mark
                        # would confirm a geometry nothing has priced.  It was
                        # only reported before; on a palladium pushed at head
                        # on it reached 0.7 A and the drag went on regardless.
                        #
                        # And a squeezed contact is not chemistry at any
                        # temperature.  This one needs no energy at all, which
                        # is why it holds where GFN2 is weakest -- a transition
                        # metal, an open shell.  Stretching stays free: that is
                        # what a reaction does and the budget prices it.
                        # Read off where the structure actually is.  Under
                        # a push that is the answer, not the geometry that was
                        # sent: what was sent carries the hand's wish, and a
                        # wish an angstrom into somebody else's atom is not a
                        # squeezed contact -- it is a question, and the answer
                        # to it is exactly that the atoms did not go there.
                        tightest = _gfn.closest_contact(
                            outcome['xyz'] if (contacts and pull is not None)
                            else current)[0]
                        crowded = (tightest is not None
                                   and tightest < _gfn._TOO_CLOSE)
                        _thermal_slope(
                            (float(priced['energy'])
                             - float(_thermal_budget()[0] or 0.0))
                            * _HARTREE_TO_KCAL,
                            current, [int(i) for i in (holding or ())]
                        ) if priced.get('energy') is not None \
                            and _thermal_budget()[0] is not None else None
                        # Whether the temperature is what refused this, or one
                        # of the two things that are refused at any
                        # temperature.  Said apart, because "past the budget"
                        # over a step the budget was perfectly happy with sends
                        # the user to the temperature box to fix something that
                        # is not there: measured on a butane turned about its
                        # middle bond, "+0.0 of 22.3 kcal/mol available. Past
                        # the budget, so the last structure that was inside it
                        # is back."
                        aside = (slipped > _SLIP_LOOSE) or crowded
                        came_back = _thermal_wall(
                            reached, priced.get('energy'), holding,
                            refuse=aside)
                        # Said after the wall has run rather than before it.
                        # The wall is what works out the highest price this
                        # drag has been at, and the line has to quote the
                        # number it is being refused on; read the other way
                        # round it quoted the maximum as it stood one step
                        # ago, which is a statement about the previous answer.
                        spent = _thermal_note(priced.get('energy'))
                        if spent:
                            # On the end of the line that is already there,
                            # never under it: this row stands above the viewer
                            # and a second one moves the picture while the
                            # user is aiming an atom.
                            said = f'{said} {spent}'
                        # What the hand was asking for when it ran out, so the
                        # next steps can tell "still pulling" from "easing
                        # off" without running anything.
                        state['thermal_spent'] = (
                            [dict(one) for one in contacts]
                            if came_back is not None and contacts else None)
                        # And the page is told, so the picture stops promising
                        # travel the kernel is not going to deliver.
                        _stop_the_hand_at(came_back, current, holding)
                        if came_back is not None:
                            # Three cases rather than two.  Which number the
                            # budget refused on -- where the structure is
                            # standing, or what the drag went through to put
                            # it there -- is as much a difference as whether
                            # the budget refused at all: only the first can be
                            # answered by easing off where you are, so a line
                            # that did not tell them apart had the user pulling
                            # harder at a wall already behind them.
                            said = (f'{said} ' + (
                                'So the last structure that was measured and '
                                'allowed is back.' if aside else
                                'Past the budget on the way here, so the last '
                                'structure that was inside it is back.'
                                if state.get('thermal_over') == 'path' else
                                'Past the budget, so the last structure that '
                                'was inside it is back.'))
                            # And what the refusal is a refusal *of*, which is
                            # narrower than it reads.  What was priced is the
                            # way this hand went, one geometry at a time; the
                            # cheapest way from the anchor to where the hand
                            # was aiming is a different quantity and nothing
                            # here has looked for it.  A hand that takes a
                            # ligand off through the middle of a ring is
                            # refused on a barrier the reaction does not have.
                            # It cannot be answered by a drag at all -- a
                            # minimum over all paths is a search -- so it is
                            # said, once, at the moment the refusal lands, and
                            # it names what does search.  Kept to one clause:
                            # this row stands above the viewer and grows down
                            # the page, and a refusal already carries the two
                            # numbers, the retreat and the slope.
                            #
                            # The path finder rather than the press it is on,
                            # which is named "Find the path" or "To the
                            # saddle" depending on the box beside it -- see
                            # _name_the_saddle_press.  Naming the one it is
                            # not showing would send the user looking for a
                            # button that is not there.
                            if not aside:
                                said = (f'{said} This prices the way your '
                                        'hand went, not the cheapest way '
                                        'there -- Scan and the path finder '
                                        'look for that.')
                                # And whether this is the one case where the
                                # quantity itself is off.  A ceiling is a free
                                # energy and a drag is priced with an
                                # electronic one; that costs under 3 kcal/mol
                                # while the pieces stay as they were, and
                                # about ten once a drag has taken something
                                # apart -- with the electronic price the
                                # strict one, so a refusal here may be
                                # refusing what the temperature would pay for.
                                # Said rather than corrected: see
                                # _THERMAL_QUANTITY_SHORT for why a number
                                # invented off a distance threshold would be
                                # worse than the gap it filled.
                                began_in = state.get('thermal_pieces')
                                now_in = _pieces_in(reached)
                                if began_in and now_in > began_in:
                                    said = (
                                        f'{said} It is in {now_in} pieces '
                                        f'here where the budget was measured '
                                        f'on {began_in}, and an electronic '
                                        f'price is strict by about ten '
                                        f'kcal/mol there -- a scan with its '
                                        f'energy set to G prices it as a free '
                                        f'energy.')
                            # And whether this was a slope or a step.
                            #
                            # A refusal on a slope can be worked with: ease
                            # off, come at it more slowly, and the structure
                            # stops somewhere near the ceiling.  A refusal on
                            # a step cannot -- there is no geometry in between
                            # at any step size, so trying again more gently
                            # lands in exactly the same place, and the user
                            # who does not know that goes on trying.
                            #
                            # Measured on the user's manganese complex under
                            # GFN2, the phenolate oxygen pulled off the metal:
                            # the price rises smoothly to -1.6 kcal/mol and
                            # the next answer is +68.0, at a Mn-O of 3.19 A --
                            # and it is +68 whether the cursor advances a
                            # quarter of an Angstrom per answer or a fiftieth.
                            # Under GFN-FF the same drag reads -20.4, ... ,
                            # +2.5, +108.7 with the cursor moving 0.05 A at a
                            # time, and +108.9 at 0.02 A. There is no path up
                            # that face; what is on the other side is a
                            # different coordination sphere.
                            #
                            # Judged against the ceiling rather than against a
                            # fixed number of kcal/mol, so it means the same
                            # thing at 100 K and at 1500 K.
                            was = state.get('thermal_priced')
                            now = priced.get('energy')
                            jump = ((float(now) - float(was))
                                    * _HARTREE_TO_KCAL
                                    if was is not None and now is not None
                                    else None)
                            _, ceiling = _thermal_budget()
                            if jump is not None and jump > max(1.0, ceiling):
                                # What changed, named, because "the bonding is
                                # different on the far side" is the fact that
                                # makes the number make sense.  Perceived here
                                # and not every answer: it is a pass over every
                                # pair, and this is a refusal rather than the
                                # ordinary case.
                                bonding = ''
                                allowed = state.get('thermal_good')
                                if allowed:
                                    bonding = _gfn.graph_changed(
                                        _gfn.bond_graph(allowed),
                                        _gfn.bond_graph(reached),
                                        [line.split()[0] for line
                                         in _gfn.atom_lines(reached)])
                                said = (
                                    f'{said} That was a step and not a slope '
                                    f'-- one answer cost {jump:+.0f} kcal/mol '
                                    f'-- so coming at it more gently lands in '
                                    f'the same place.'
                                    + (f' On the far side the bonding '
                                       f'{bonding}.' if bonding else ''))
                        if crowded:
                            said = (f'{said} Two atoms are inside '
                                    f'{tightest:.2f} of a bond length, which '
                                    f'is no path at any temperature.')
                        # Both of these belong to the budget, so they are said
                        # where the budget is: read one indent further out,
                        # they are read on a drag that never set them, and the
                        # ordinary drag -- the one with no budget at all --
                        # dies on the first step.
                        slope = state.get('thermal_slope')
                        if slope is not None and abs(slope) > 1.0:
                            said = (f'{said} '
                                    + ('Climbing' if slope > 0 else 'Falling')
                                    + f' {abs(slope):.0f} kcal/mol per A here.')
                        if slipped > _SLIP_LOOSE:
                            said = (f'{said} The hold is loose here '
                                    f'({slipped:.2f} A), so the price is for '
                                    f'a nearby structure rather than this '
                                    f'one.')
                        # What this answer cost, so the next one can say
                        # whether the drag walked up to the wall or fell over
                        # a step onto the far side of it.  Written last, after
                        # everything that reads the previous value.
                        if priced.get('energy') is not None:
                            state['thermal_priced'] = float(priced['energy'])
                    state['gfn_last_status'] = said
                    schedule_ui_update(_set_mol_status, said, spinner=True)
                    # The atoms under the cursor go back where they were sent.
                    # xtb pulls them most of the way home in five cycles, and
                    # this answer outlives the drag -- applied after the
                    # release it would take them with it, which is the spring
                    # back that looked like the drag being undone.
                    # Unless a turn is being held, and then the torsion is
                    # what the hand asked for and a position would fight it.
                    # Put back where the cursor had them, a methyl turned
                    # about a chain bond cannot spin about its own axis to
                    # get out of the way, and the same turn that costs +7.9
                    # kcal/mol along the torsion costs +345 held rigid.
                    if contacts:
                        state['thermal_holding'] = [dict(one) for one in contacts]
                    turned = [one for one in contacts
                              if str(one.get('kind')) == 'dihedral']
                    if turned:
                        state['thermal_turn'] = list(turned[0]['atoms'])
                    # Price what you show, and show it where the hand left it.
                    #
                    # The atoms used to be put back where the cursor had them
                    # and the *relaxed* energy reported, so when the held
                    # values did not determine the drag the two came apart:
                    # the picture had a bromide 1.27 A from a palladium and
                    # the price belonged to a structure where the metal had
                    # got out of the way.  The relaxed geometry is the one
                    # that was priced, so it is the one that is drawn.
                    #
                    # But a held value is an internal coordinate, so xtb meets
                    # it and is free to put the molecule anywhere that does --
                    # and the answer comes back slid a little.  Drawn as it
                    # arrives, the atom leaves the cursor and the next mouse
                    # move brings it back, seven times a second: measured on a
                    # backside attack, 0.102 A each way through the transition
                    # region, which is a molecule that shakes.  Laid back onto
                    # the hand as a rigid body it costs nothing -- an energy
                    # does not depend on where a molecule is -- and the two are
                    # still the same structure.
                    #
                    # Which of the two placements this is was decided further
                    # up, where the price was taken, because the budget has to
                    # be asked about the geometry that is going to survive.
                    settled = came_back if came_back is not None else reached
                    # Whatever the budget said, the molecule stays the
                    # molecule it was if that was asked for.  Judged on what
                    # the user would be left with rather than on the raw
                    # answer, because that is what outlives the drag.
                    kept, changed = _topology_wall(settled)
                    if kept is not None:
                        settled = kept
                        said = (f'{said} Held: that step {changed}, and the '
                                f'bonding is being kept.')
                        # A rigid hand and Keep bonds are asking for different
                        # things, and when the cursor wants a place that
                        # breaks a bond one of them has to lose.  Keeping the
                        # bonding wins, and then nothing moves at all -- which
                        # on screen is indistinguishable from a drag that is
                        # not working.  Counted, so it is said rather than
                        # merely happening.
                        refused = int(state.get('topology_refused') or 0) + 1
                        state['topology_refused'] = refused
                        if refused >= 3 and pull is None:
                            said = (f'{said} Three steps running -- the hand '
                                    f'is set to move the atom, which is '
                                    f'asking for a place that breaks a bond. '
                                    f'Pull instead, or turn Keep bonds off.')
                    else:
                        state['topology_refused'] = 0
                    # What the box is to be left holding, and why.
                    #
                    # One write, at the end, of the geometry that survives the
                    # step -- rather than one write per refusal and nothing at
                    # all when the step was allowed.  The box is what outlives
                    # the drag: it is what Copy and Submit read and what the
                    # next calculation starts from, so whatever reaches it has
                    # to be something that has been priced.  Nothing wrote it
                    # on an allowed step, so the page's own model wrote it
                    # instead -- and the page's model is the mouse's wish, at
                    # +136.7 kcal/mol under a 22.3 ceiling while the line
                    # underneath read "+16.2 of 22.3 available".
                    why = ''
                    if kept is not None:
                        why = 'Kept: the bonding would have changed'
                    elif came_back is not None:
                        why = ('Back to the last structure that was measured '
                               'and allowed' if aside else
                               'Past the budget: back to the last structure '
                               'that was inside it')
                    elif pricing:
                        why = ('Within the budget at '
                               f'{float(submit_temperature.value):g} K')
                    if why:
                        rows = [line for line in settled.splitlines()[2:]
                                if line.strip()]
                        if rows:
                            schedule_ui_update(
                                _write_coords, xyz_document(rows, why), True)
                    # What the next answer measures the hand against: the
                    # geometry this one handed back, not the one it was
                    # handed.  Against the latter the difference holds the
                    # relaxation as well as the hand, and on a ring being
                    # puckered the relaxation is the larger of the two -- so
                    # every drag read as a drive at something, the torsion
                    # never fired, and a cyclohexane could not be flipped.
                    state['thermal_was'] = settled
                    frames = list(state.get('gfn_follow_frames') or [])
                    frames.append(_gfn.coordinates_of(settled))
                    state['gfn_follow_frames'] = frames
                    # The tail, not the whole drag: every write is a message,
                    # and a minute of dragging is three hundred frames of a
                    # hundred atoms in each of them.
                    trail = frames[-40:]
                    drag_run = state.get('gfn_follow_run')
                    payload = _frame_payload(
                        drag_run, **{'from': len(frames) - len(trail),
                                     'follow': 1, 'frames': trail})
                    # Asked at the write and not at the answer: an xtb
                    # round is long enough for the run under it to have
                    # been replaced, and the frame is then about a
                    # structure that is no longer anywhere.
                    schedule_ui_update(
                        lambda text=payload, r=drag_run: setattr(
                            submit_gfn_frame, 'value', text)
                        if _frame_run_is_current(r) else None)
            finally:
                state['gfn_follow_busy'] = False

        _start_background(_work, 'The relaxation under the hand')

    #: Letting go with Settle on.  More cycles than a follow step, because
    #: nothing is being held any more and the structure should reach somewhere
    #: it would stay rather than take one step towards it -- but still a bound,
    #: because this happens on every release.
    _GFN_SETTLE_CYCLES = 40

    #: How many rounds before it gives up on converging.  Held values that
    #: cannot all be met at once never converge -- measured on a propane with
    #: two fixed distances and an angle fighting each other, not converged in
    #: any round -- and a relaxation that will not end is a process per round
    #: for as long as the switch is down, and a structure that visibly jitters.
    #: How many xtb runs one press of Optimise may chain before it stops and
    #: says so.
    #:
    #: A run is one xtb ``--opt`` and xtb's optimiser has a cycle limit of its
    #: own: at the limit it hands back the geometry it reached and reports that
    #: it did not converge.  That was taken as the end, so the switch went up
    #: over a structure that was better than it had been and not at a minimum,
    #: and the user pressed again -- which worked, because pressing again is a
    #: new process from that geometry with a fresh cycle budget and a fresh
    #: optimiser history.  Measured on a pulled-about propane under GFN-FF at
    #: three cycles a run: three runs to converge, the largest shift falling
    #: 0.253, 0.066, 0.002 A.  The pressing is done here now.
    _OPTIMISE_ROUNDS = 12
    #: A structure that has stopped improving is finished as far as this can
    #: take it, whatever xtb says about convergence.  Same figure as a settle:
    #: below this, another round is another process for nothing.
    _OPTIMISE_STILL = 0.005
    #: How many of those in a row it takes to believe it.  One is not proof: a
    #: run that reached its cycle limit in a flat stretch moves almost nothing
    #: and its successor -- a new process, fresh budget, fresh optimiser
    #: history -- can still take a real step.  Two costs one xtb run at most.
    _OPTIMISE_STILL_ROUNDS = 2

    _GFN_SETTLE_ROUNDS = 12
    #: And a round that moved nothing has settled, whatever xtb calls it.
    _GFN_SETTLE_STILL = 0.005

    # The same bound for the walk that goes the other way -- how many
    # gradients one climb may spend before it stops and says so -- lives in
    # the module that walks it, as :data:`climb.CLIMB_CEILING`, and is not
    # repeated here.
    #
    # There were two of them and they disagreed: four hundred here, chosen as
    # an order of magnitude above every climb that had finished at the time,
    # and a hundred there, measured afterwards over twenty-one hand drags --
    # every climb that reached the reaction the hand pointed at arrived in 22
    # to 65 steps, and the ones that ran longer ran on to between 202 and 361
    # and were wrong at the end of all of them.  A hundred is the measured
    # one, so a hundred is what is kept; the other three hundred bought
    # nothing but a longer wait before the next rung was tried.

    def _write_coords(text, drawn=False, run=None):
        """Put a geometry in the box, and say whether the picture has it.

        *run* is the run number the geometry was computed under, for anything
        that answers from a thread. A run whose number has moved on is about a
        structure that is no longer on screen -- Undo or Reset put an older one
        there, or an edit did -- and its answer is refused rather than written.
        Without this, pressing Undo while xtb was minimising gave the geometry
        back for the few seconds the run had left to go and then lost it again,
        which from outside is a button that does nothing.

        *drawn* raises the flag the host's update_view consumes: the playback
        has already drawn this geometry, so redrawing it would rebuild the
        viewer and tear down what is still playing.

        The flag is consumed by update_view, and traitlets only calls that when
        the value actually changes.  Raised before a write that changes
        nothing, it stays raised -- and the next genuine change is the one that
        gets swallowed instead.  A structure already at its minimum written
        back by a second press of Optimise is exactly that write, and from
        there the picture stops following the box until something changes the
        value twice over, which is what cutting the coordinates out and pasting
        them back does.  So a write that changes nothing lowers the flag rather
        than raising it.
        """
        if run is not None and int(state.get('gfn_run', 0)) != int(run):
            return False
        if coords_widget.value == text:
            state['manip_inflight'] = False
            return False
        state['manip_inflight'] = bool(drawn)
        # The same molecule, changed -- never a different one.
        #
        # The host treats any write to the box as a structure it has not seen
        # unless it is told otherwise, and starts that one over: charge back
        # to zero, multiplicity to one, and every held value, bond edit and
        # history entry thrown away.  That is right for a structure someone
        # loads and wrong for every write this editor makes, and a scan makes
        # one per point -- so a scan lost the charge it was told to run at,
        # mid-scan, and the undo history with it.
        state['structure_edit_inflight'] = True
        try:
            coords_widget.value = text
        finally:
            state['structure_edit_inflight'] = False
        return True

    def _thermal_budget():
        """What this structure may spend, in kcal/mol, and from where.

        Returns ``(anchor_energy_Eh, ceiling_kcal)`` or ``(None, ceiling)``
        while there is no anchor yet.
        """
        ceiling = thermal_ceiling(submit_temperature.value, _THERMAL_SECONDS)
        # The anchor names the molecule it was measured on.  Asked of another
        # one it is not an anchor at all, and the difference between two
        # unrelated energies would read as a distortion of the second.  Keyed
        # on the element column the way the perception and the GFN-FF topology
        # are, so it guards itself rather than depending on a callback firing
        # -- structure_changed is called from one place, and a structure can
        # arrive without it.
        if state.get('thermal_for') != _structure_fingerprint(_current_xyz() or ''):
            return None, ceiling
        # And to the method that measured it.  An anchor is an energy, and an
        # energy of one method against energies of another is not a
        # difference at all: measured, an anchor taken under GFN2 and then
        # read against GFN-FF priced an untouched structure at +6384 kcal/mol
        # against a 22.3 ceiling, so every drag sprang straight back.  The
        # other way round is quieter and worse -- a C-C torn to 3 A reported
        # "-6117.9 of 22.3 available" and the wall never fired.
        if state.get('thermal_method') != str(submit_ff_dd.value):
            return None, ceiling
        return state.get('thermal_e0'), ceiling

    def _thermal_live():
        """Whether the budget is switched on *and* has a hand it can price.

        The switch alone is not enough.  A budget prices what is kept, and
        under a placing hand what is kept is not exactly what was priced --
        measured at 1.4 kcal/mol, +16.8 priced against +18.2 kept, bounded by
        _SLIP_LOOSE rather than by zero.  So the budget does not act there;
        see :func:`_refresh_hand_controls`, which is where it also stops being
        on screen.  The switch keeps its value across the change, because a
        hand is changed constantly and re-measuring the anchor each time is a
        calculation the user did not ask for.
        """
        return bool(submit_thermal_btn.value) and _hand_pulls()

    def _thermal_note(energy):
        """Where this geometry stands against the budget, said in one line.

        Two numbers, and the line says which of them it is talking about.
        Where the structure stands is one thing; what the drag went through to
        put it there is another, and only the second is a barrier.  "You are
        standing somewhere expensive" can be fixed by easing off where you
        are; "you got here through somewhere expensive" cannot be fixed at all
        except by coming back, and a line that showed only the first number
        would have the user pulling harder at a wall they had already gone
        past.
        """
        anchor, ceiling = _thermal_budget()
        if anchor is None or energy is None:
            return ''
        spent = (float(energy) - float(anchor)) * _HARTREE_TO_KCAL
        T = float(submit_temperature.value)
        was = state.get('thermal_peak')
        peak = spent if was is None else max(float(was), spent)
        if peak <= ceiling:
            return (f'{spent:+.1f} of {ceiling:.1f} kcal/mol available at '
                    f'{T:g} K.')
        # Said as a time rather than as a number: how long this would take is
        # the thing a chemist can judge, and 34.7 kcal/mol means nothing until
        # it is "longer than the age of the earth".
        if spent > ceiling:
            return (f'{spent:+.1f} kcal/mol -- past the {ceiling:.1f} this '
                    f'structure has at {T:g} K. {_thermal_wait(spent, T)} '
                    f'{_thermal_wants(spent)}')
        # Cheap here, and it got here over something that was not.  The
        # refusal is about the crossing, so the crossing is the number quoted
        # and the wait is worked out from it.
        #
        # And what to do about it, said in the same breath.  This is the one
        # refusal with an obvious next move and no way of guessing it: the
        # budget was spent getting here, the structure is standing somewhere
        # cheap, and whether that somewhere is a state the system would have
        # thermalised in is a chemist's judgement rather than a calculation.
        # Set here is how that judgement is entered, and a user who does not
        # know the button exists discovers instead that switching Dynamik Opt
        # off and on again makes the drag work -- which is the same thing
        # happening by accident.
        return (f'{spent:+.1f} of {ceiling:.1f} kcal/mol available at '
                f'{T:g} K, but this drag went through {peak:+.1f} to reach '
                f'it -- past the {ceiling:.1f} it has at {T:g} K. '
                f'{_thermal_wait(peak, T)} {_thermal_wants(peak)} '
                f'The budget is counted from where it was set; if the '
                f'structure has settled here, Set here measures from this '
                f'one instead.')

    def _thermal_wait(kcal, temperature):
        """How long a barrier of that height takes at that temperature.

        Down to picoseconds, not only down to seconds.  A scan that finds an
        open path is the ordinary case and it was reported as "about 4.18e-06
        s", which is a number in the wrong clothes: the answer wanted there is
        "4 microseconds", and below a picosecond there is no crossing to speak
        of anyway.

        And it stops at the top for the same reason it started at the bottom.
        A refusal on a barrier the temperature is nowhere near came out as
        "that is about 3.56e+29 years", which is a figure nobody reads as a
        quantity -- the sentence this was written to produce is "longer than
        the age of the earth", and past about ten billion years every number
        means the same thing.  Said in years up to there and in words past it.
        """
        T = max(1.0, float(temperature))
        rate = ((_BOLTZMANN_SI * T / _PLANCK_SI)
                * math.exp(-max(0.0, float(kcal)) / (_GAS_CONSTANT * T)))
        if rate <= 0:
            return 'It does not happen.'
        seconds = 1.0 / rate
        if seconds > _UNIVERSE_SECONDS:
            return 'That is longer than the universe has existed.'
        for limit, unit, name in ((1e-9, 1e-12, 'ps'), (1e-6, 1e-9, 'ns'),
                                  (1e-3, 1e-6, 'us'), (1.0, 1e-3, 'ms'),
                                  (60, 1, 's'), (3600, 60, 'min'),
                                  (86400, 3600, 'h'), (3.15576e7, 86400, 'd'),
                                  (float('inf'), 3.15576e7, 'years')):
            if seconds < limit:
                return f'That is about {seconds / unit:.3g} {name}.'
        return ''

    def _thermal_wants(kcal):
        """The temperature that price wants, said as an answer not a refusal.

        The ceiling is Eyring inverted, so running the same arithmetic the
        other way turns any price already in hand into the temperature at
        which it becomes possible.  It costs nothing -- the energy has been
        paid for, and this is two logarithms -- and it is a strictly better
        answer than "no": "refused at 298 K" is where a user stops, "it wants
        about 480 K" is where they go next, and the second is the question
        they came with.  The temperature box above is where the answer is
        used, which is why no control was added for it: the reverse question
        is answered in the sentence and typed into the box that was already
        there.

        Both halves are said, here and in :func:`_thermal_wait`, because they
        are two different questions with the same arithmetic behind them: how
        hot for the window, and how long at the temperature set.  A user who
        cannot heat it reads the first and stops; one who can wait reads the
        second and does not need a window control to find out.
        """
        needs = thermal_temperature(kcal, _THERMAL_SECONDS)
        if needs is None:
            return (f'No temperature under 5000 K crosses that within '
                    f'{_timescale_label()}.')
        return (f'It wants about {needs:.0f} K ({needs - 273.15:+.0f} C) '
                f'within {_timescale_label()}.')

    def _pieces_in(xyz):
        """How many separate pieces that geometry is in.

        Off the same bond graph the topology wall is judged by, so there is
        one perception in this file rather than two that can disagree, and
        asked only where a refusal has already landed -- it is a pass over
        every pair.

        What it is for is the one case where the budget is wrong by more than
        the method is.  The ceiling is a free energy and a drag is priced with
        an electronic one; while the pieces stay as they were the two agree to
        under 3 kcal/mol, and where a drag takes something apart they part
        company by about ten at room temperature, with the electronic price
        the strict one.  See _THERMAL_QUANTITY_SHORT.  Nothing is corrected
        here -- a number invented off a distance threshold would be worse than
        the gap -- but a refusal can say which case it is in, and that is the
        difference between a verdict and a verdict the user can weigh.
        """
        rows = _gfn.atom_lines(xyz or '')
        if not rows:
            return 0
        parent = list(range(len(rows)))

        def home(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for i, j in _gfn.bond_graph(xyz):
            one, two = home(i), home(j)
            if one != two:
                parent[one] = two
        return len({home(i) for i in range(len(rows))})

    def _thermal_slope(spent, xyz, serials):
        """What an angstrom costs here, from the last two answers.

        All that is left of the leash.  That held the hand to a fraction of an
        angstrom per answer so every geometry could be priced before the next
        was allowed, and it made the drag unusable -- and it was not even safe
        once it was long enough to feel right: lengthened on ground the last
        answers called flat it stood at 1.000 A going into a C-H bond 1.09 A
        long, which is a bond torn in a single frame.  The hand is free now and
        the price arrives behind it; see _thermal_wall.

        The number itself is worth keeping and saying.  "+8 per A" and "+160
        per A" are two different situations, and knowing which one you are in
        is the difference between working with the chemistry and against it.
        """
        here = _gfn.coordinates_of(xyz)
        last = state.get('thermal_last') or {}
        state['thermal_last'] = {
            'spent': float(spent),
            'at': {s: [here[3 * s], here[3 * s + 1], here[3 * s + 2]]
                   for s in serials if 3 * s + 2 < len(here)},
        }
        moved = 0.0
        for s in serials:
            was = (last.get('at') or {}).get(s)
            now = state['thermal_last']['at'].get(s)
            if was and now:
                moved = max(moved, math.dist(was, now))
        if moved <= 1e-6 or 'spent' not in last:
            return
        state['thermal_slope'] = (float(spent) - float(last['spent'])) / moved

    def _arm_thermal_leash():
        """Remember where the budget still agreed, so there is somewhere back to.

        There used to be a leash here: the hand was held to a fraction of an
        Angstrom per answer so every geometry on the way could be priced
        before the next one was allowed.  It is rigorous and it is unusable --
        the drag becomes a slideshow, and the one thing that made the editor
        worth using went with it.  Worse, it is not even safe once it is long
        enough to feel right: lengthened on flat ground it stood at 1.000 A
        going into a C-H bond 1.09 A long, which is a bond torn in one frame,
        before anything had been asked.

        So the hand is free and the answer arrives behind it.  What this keeps
        is the geometry to come back to.
        """
        if not _thermal_live():
            return
        anchor, _ceiling = _thermal_budget()
        if anchor is None:
            return
        here = _current_xyz() or ''
        if here:
            state['thermal_good'] = here
            # The pair travels together: a geometry to come back to, and what
            # the way to it cost.  Nothing has priced this one yet, so there
            # is no maximum to record against it.
            state.pop('thermal_good_peak', None)

    def _back_at(xyz, good):
        """Whether the structure is standing where the budget last agreed.

        Laid onto the other as a rigid body before the two are compared.  An
        energy does not depend on where a molecule is, and neither does
        whether it is the same structure; without the fit, a drag that had
        turned the whole molecule a little read as a different geometry at
        every atom, and the running maximum below would never have come down.

        The threshold is the one the loose hold is already judged by --
        _SLIP_LOOSE, 0.25 A on the atom that moved most -- so there is one
        number in this file for "near enough to be the same placement" rather
        than two that can drift apart.
        """
        count = len(_gfn.atom_lines(xyz or ''))
        if not good or not count or count != len(_gfn.atom_lines(good)):
            return False
        return _gfn.largest_shift(
            good, _gfn.settle_onto(xyz, good, range(count))) <= _SLIP_LOOSE

    def _thermal_wall(xyz, energy, holding, refuse=False):
        """Keep what the budget agreed to, and take back what it did not.

        *xyz* is the geometry that will be kept if this says yes, and *energy*
        is the energy of that same geometry.  Both halves of that sentence are
        load-bearing and neither used to hold: what was handed in was the
        page's own model -- the mouse's wish -- while the energy belonged to
        the structure xtb had relaxed around it.  Under a pull those two part
        company completely, because the push is soft and the bond comes most
        of the way back, so the budget was answering about one molecule and
        the user was left with another.  Measured on an ethane methyl pulled
        out at 298 K against a 22.3 ceiling, the relaxed answers were +0.1,
        +1.1, +0.2, +4.0, +6.8, +12.2, +16.2 kcal/mol and every one of them was
        allowed, while what stayed in the box went +15.8, +45.3, +74.8, +99.1,
        +116.9, +129.0, +136.7.  The geometry kept back as "the last one that
        was inside the budget" was the same wish, at +136.7.

        Act, then undo, rather than hold back.  The hand cannot be stopped at
        the right place in real time -- xtb is far too slow for that and always
        will be -- but nothing that was not allowed has to be *kept*.  So the
        drag runs free, the price arrives a moment behind it, and a geometry
        that turns out to be past the budget is replaced by the last one that
        was inside it.  The structure springs back, visibly, and what is left
        in the box has always been priced and allowed.

        What this gives up is that the path between two answers is no longer
        sampled without gaps.  For a drag whose price climbs steadily -- a bond
        stretched, an atom pulled out -- the endpoint catches it every time.
        Something that jumps a narrow barrier and is cheap again on the far
        side could slip through, and that is what the scan is for: it walks the
        coordinate rather than following a hand.

        What the ceiling is, though, is a barrier *height*, so what is held
        against it is the highest price this drag has been at since the anchor
        and not the price of wherever it happens to be standing.  A drag that
        crossed +32 and settled at 0 has still crossed +32, and at 298 K --
        where the hour buys 22.3 -- it could not have.  On the point alone the
        far side was kept: measured on the hindered N-N rotation of a
        nitrosamine under GFN2, the wall refused the two steps over the top at
        +26.5 and +32.2 and then allowed every step of the far side, +13.7,
        +7.9, +3.5, +0.8 and +0.0, so what the drag left behind was a
        structure reached over a barrier room temperature cannot pay for
        within the hour.  The maximum is carried, and it is what refuses.

        Carried, and not carried for ever.  A maximum that never comes down
        makes the editor unusable in the other direction: one jerk of the
        mouse past the ceiling and every later step of that grab is refused,
        including the ones walking straight back home.  What undoes a crossing
        is coming back over it, so the maximum drops to the one recorded for
        the last geometry the budget agreed to as soon as the structure is
        standing at that geometry again -- which, after a rollback, is exactly
        where it has just been put.  An overshoot then costs one refused step
        instead of the rest of the grab; a far side the structure is still
        standing on costs every step, because it never comes back.  Letting go
        clears it: what a release leaves behind was inside the budget the
        whole way, so the next grab starts on an affordable path.

        Where this still does not mean what it looks like: the sampling gap
        above is unchanged, so a top that falls between two answers is a top
        nothing has priced, and the retreat is judged on placement -- a
        crossing whose whole width is inside the 0.25 A that counts as
        standing in the same place is forgiven.  Both are what the scan is
        for, which walks the coordinate instead of following a hand.
        """
        anchor, ceiling = _thermal_budget()
        if anchor is None or energy is None:
            return None
        spent = (float(energy) - float(anchor)) * _HARTREE_TO_KCAL
        good = state.get('thermal_good')
        carried = state.get('thermal_peak')
        if _back_at(xyz, good):
            # The excursion has been retracted, so it is not on the way here
            # any more: what this geometry cost to reach is what it cost when
            # it was agreed to.  Read the other way -- while the structure is
            # still out where it went -- the maximum stands and every step is
            # refused, which is the whole point.
            carried = state.get('thermal_good_peak')
        peak = spent if carried is None else max(float(carried), spent)
        state['thermal_peak'] = peak
        # Which number is doing the refusing, for the line that says so.
        state['thermal_over'] = ('here' if spent > ceiling
                                 else 'path' if peak > ceiling else '')
        if peak <= ceiling and not refuse:
            state['thermal_good'] = xyz
            state['thermal_good_peak'] = peak
            return None
        # Past it.  Hand back the last geometry that was not, if there is one;
        # a drag that was already over budget when it started has nowhere to
        # go, and then the hand simply keeps what it has and the line says so.
        #
        # And it has to be *this* molecule.  Nothing cleared it, so a benzene
        # kept from an earlier drag was handed back during a drag on water:
        # measured, the coordinate box became twelve atoms of benzene and a
        # thirty-six float frame went to a three-atom viewer, with nothing
        # anywhere saying the molecule had been replaced.
        if good and len(_gfn.atom_lines(good)) != len(_gfn.atom_lines(xyz)):
            state.pop('thermal_good', None)
            state.pop('thermal_good_peak', None)
            return None
        return good

    def _settle_price(outcome, constraints):
        """What a relaxation reached, with any held value's own energy out.

        The same arithmetic the drag is priced by, for the same reason: a
        value held with a force leaves a real restraint energy in xtb's total,
        and charging the budget for the hand as well as for the structure
        prices something nobody is looking at.  A value met exactly leaves
        none worth the name, so the answer is then its own price.
        """
        if not outcome.get('ok') or outcome.get('energy') is None:
            return None
        if not constraints:
            return float(outcome['energy'])
        bias = _gfn.restraint_energy(outcome['xyz'], constraints, _value_in)
        return (float(outcome['energy']) - bias if bias is not None
                else float(outcome['energy']))

    def _keep_the_priced_geometry():
        """Leave the box holding the last geometry the budget agreed to.

        The follow writes every step it prices, so ordinarily the box has this
        already and nothing happens here.  What it is for is the release: the
        page sends its own model when the hand lets go, and that model is
        where the *cursor* was rather than where the chemistry allowed -- so
        the last word on a drag was the one geometry nothing had ever priced.
        """
        good = state.get('thermal_good')
        if not good:
            return
        rows = [line for line in good.splitlines()[2:] if line.strip()]
        here = coords_widget.value or ''
        if not rows or len(_gfn.atom_lines(here)) != len(rows):
            # A different molecule, so the kept geometry is not about this
            # one and handing it over would swap the structure underneath the
            # user.  The wall guards itself the same way, for the same reason.
            return
        if _gfn.coordinates_of(here) == _gfn.coordinates_of(good):
            return          # already holding it, and the comment says why
        _write_coords(
            xyz_document(rows, 'Within the budget at '
                               f'{float(submit_temperature.value):g} K'),
            True)

    def _topology_wall(xyz):
        """Keep the molecule the molecule it was, and take back what did not.

        The same shape as the budget's wall and for the same reason: xtb
        cannot be told to hold a topology, and a hand cannot be stopped at the
        right place in real time.  So the step runs, the bonding is read off
        what came back, and a step that made or broke one is replaced by the
        last one that did not.

        Returns the geometry to go back to, or None when nothing changed.  The
        second value is what changed, for the line.
        """
        if not submit_topology_btn.value or not xyz:
            state.pop('topology_good', None)
            return None, ''
        was = state.get('topology_graph')
        rows = [line.split() for line in _gfn.atom_lines(xyz)]
        who = _structure_fingerprint(xyz)
        if was is None or state.get('topology_for') != who:
            # First look, or a different molecule: this one is the one to
            # keep.  Carried over, a benzene's bonding would be held against
            # a water and every step would read as a change.
            state['topology_graph'] = _gfn.bond_graph(xyz)
            state['topology_for'] = who
            state['topology_good'] = xyz
            return None, ''
        # With hysteresis, or a bond resting on the threshold decides the
        # answer differently ten times a second and the drag sticks on a
        # molecule that is not changing at all.
        holds, said = _gfn.graph_holds(was, xyz)
        if holds:
            state['topology_good'] = xyz
            return None, ''
        return state.get('topology_good'), said

    def _push_thermal_wall(wall, reach=0.0):
        """Tell the page where the hand may no longer go, or that it may.

        Through a widget the page reads on its own clock, never through
        run_js: that clears its output before displaying, so of two calls in
        quick succession only the second survives -- and this is written about
        ten times a second while a drag is running.  Sent that way it never
        arrived at all.

        The serial makes two identical walls two different values, or
        traitlets would say nothing the second time and a wall that went down
        and up again in the same place would stay down.
        """
        state['thermal_wall_serial'] = int(
            state.get('thermal_wall_serial') or 0) + 1
        payload = {'n': state['thermal_wall_serial'],
                   'reach': float(reach),
                   'wall': ({str(k): v for k, v in wall.items()}
                            if wall else None)}
        schedule_ui_update(
            lambda text=json.dumps(payload): setattr(
                submit_gfn_wall, 'value', text))

    def _stop_the_hand_at(came_back, asked, holding):
        """A refused drag stops following the cursor; an allowed one follows.

        *came_back* is what the wall handed back, or ``None`` when it allowed
        the step.  *asked* is the geometry the page sent, which under a pull
        carries the wish rather than the atom: the browser writes the held
        atoms at the point the hand is asking for, clamped to a reach.

        Nothing told the page a step had been refused.  The kernel stops --
        _still_spent says there is nothing left to compute until the hand eases
        -- and the page went on running the wish out with the cursor, so the
        band grew, the coordinates the page reported went on changing, and the
        one thing that had actually happened, that the drag had reached its
        ceiling, was the one thing not on screen.

        So the wish is marked where the budget last agreed the atom could
        stand, and may go no further out than it had already run: an atom that
        has come up against something.  Coming back in is never blocked --
        thermalWallBlocks refuses only a step that is both outside the reach
        and going further out -- and coming back in is exactly the event
        _still_spent reads as the hand easing off.  The two sides then let the
        same gesture through instead of holding two opinions about it.

        The reach is where the wish had run to and not zero, so the band stops
        rather than snapping back to the atom: what the hand was asking for
        when it was refused is a true thing to leave on the screen.
        """
        if came_back is None:
            if state.pop('thermal_walled', None):
                _push_thermal_wall(None)
            return
        marks = {}
        reach = 0.0
        here = _gfn.coordinates_of(came_back)
        wished = _gfn.coordinates_of(asked or '')
        for one in (holding or ()):
            index = int(one)
            if 3 * index + 2 >= len(here) or 3 * index + 2 >= len(wished):
                continue
            mark = [here[3 * index + k] for k in range(3)]
            marks[index] = mark
            reach = max(reach, math.dist(
                mark, [wished[3 * index + k] for k in range(3)]))
        if not marks or reach <= 0:
            return
        state['thermal_walled'] = True
        _push_thermal_wall(marks, reach)

    def _set_thermal_anchor(relax=None, note='Measuring from here',
                            note_after=''):
        """Take the energy of the structure on screen as the budget's zero.

        A single point when the structure is to be kept as it is, one
        optimisation when it is not.  Either way what is stored is an energy
        of the *chosen method*, so the budget and the drag are the same
        calculation and their difference means something.

        *note_after* is anything the caller has to say about what switching
        the budget on has changed besides the anchor, said on the end of the
        line that reports the ceiling rather than in a line of its own -- this
        row stands above the viewer, and a second one moves the picture.
        """
        xyz = _current_xyz()
        method = str(submit_ff_dd.value)
        if not xyz or not _gfn.is_gfn_method(method):
            return
        wants_relax = (submit_thermal_relax.value if relax is None else relax)
        # A hand on the molecule outranks the setting.  Switching the budget on
        # mid-drag ran the relax-first anchor, which optimised the structure
        # and wrote the result into the coordinate box -- so the answer to
        # "measure what I am doing" was to undo it.  The structure the user is
        # holding is the one they mean; it is measured as it stands.
        held = bool(state.get('gfn_follow'))
        if held:
            wants_relax = False
            note = 'Measuring from the structure as it stands'
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        model = _solv_model()
        state['thermal_e0'] = None
        _set_mol_status(f'{note}: {_server_label(method)} is measuring the '
                        'energy of this structure...', spinner=True)

        def _work():
            outcome = _gfn.optimize_with_gfn(
                xyz, method, charge=charge, uhf=uhf, timeout=None,
                solvent=wet, solvation_model=model,
                topology=(_gfn_topology_dir(xyz)
                          if _gfn.is_gfn_method(method) else None),
                optimise=bool(wants_relax),
            )

            def _done():
                if not outcome.get('ok') or outcome.get('energy') is None:
                    _set_mol_status(
                        'The thermal budget has no anchor: '
                        + str(outcome.get('status') or 'the energy did not '
                              'come back.'))
                    return
                state['thermal_e0'] = float(outcome['energy'])
                state['thermal_method'] = method
                state['thermal_for'] = _structure_fingerprint(
                    outcome.get('xyz') or xyz)
                # A new zero is a new path, and nothing has been crossed on
                # the way to it.  Without this, Set here -- pressed on the
                # intermediate a drag has just reached, which is exactly what
                # it is for -- would be refused for ever by the barrier the
                # drag came over to get there.
                state['thermal_peak'] = 0.0
                state['thermal_good_peak'] = 0.0
                state.pop('thermal_over', None)
                # An optimisation takes as long as it takes, and the editor is
                # not frozen while it runs.  Writing its answer into the box
                # unconditionally overwrote whatever had been done since --
                # so the anchor is only allowed to move the structure it was
                # actually measuring.  When it has moved on, the anchor stands
                # for the structure it priced and the budget stays off until
                # one is taken here: which is what the fingerprint already
                # decides, so all that is needed is to say so.
                moved_on = _current_xyz() != xyz
                if wants_relax and not moved_on:
                    lines = [line for line in outcome['xyz'].splitlines()[2:]
                             if line.strip()]
                    if lines:
                        # Remembered here rather than at the press, the way the
                        # saddle search does it: the press only asks, and until
                        # the answer arrives there is nothing to take back.
                        # Without it, taking an anchor from a relaxed structure
                        # replaced the geometry and left no entry for it.
                        _remember('measuring the budget from a relaxed structure')
                        _write_coords(xyz_document(
                            lines, 'Relaxed, and the budget measured from here'))
                if moved_on:
                    _set_mol_status(
                        'The structure changed while its energy was being '
                        'measured, so that anchor belongs to the older one. '
                        'Press Set here to take one for this.')
                    return
                # The zero is also the last geometry the budget agreed
                # to, and its way here cost nothing by construction.  Left
                # pointing at an earlier drag's structure, the maximum would
                # be measured against a geometry this anchor knows nothing
                # about.
                state['thermal_good'] = _current_xyz() or xyz
                # And how many separate pieces it was measured on, so that a
                # refusal can say when the drag has changed that.  Taken here
                # and once: it belongs to the anchor the way the fingerprint
                # and the method do, and a count read at the moment of a
                # refusal would have nothing to be a change *from*.
                state['thermal_pieces'] = _pieces_in(state['thermal_good'])
                _, ceiling = _thermal_budget()
                _set_mol_status(
                    f'{note}. At {float(submit_temperature.value):g} K '
                    f'this structure has {ceiling:.1f} kcal/mol to spend '
                    f'within {_timescale_label()}.{note_after}',
                    _THERMAL_QUANTITY_SHORT)

            schedule_ui_update(_done)

        _start_background(_work, 'The energy this budget measures from')

    def _timescale_label():
        """The window the ceiling is quoted over, said in words."""
        return 'an hour'

    def _gfn_topology_dir(xyz):
        """Where GFN-FF's perceived bonding is kept while a structure is worked
        on.

        GFN-FF reads its topology out of the geometry it is handed, once.  On a
        drag that is a fresh perception per step, and it has a cliff in it:
        measured on a propane, a C-C pulled to 1.87 A relaxes back to 1.49,
        and at 1.96 A the bond is not seen at all and the same relaxation
        pushes it to 2.80.  So the perception is made once and kept, the way
        the browser's field assigns its parameters once when the switch goes
        down -- at 2.33 A it then still pulls back, to 1.51.

        It belongs to one molecule, and the element column is what says which:
        the count alone said a benzene and a cyclobutane were the same
        molecule, both being twelve atoms, so the second was optimised with the
        first one's bonding. It came back with a hydrogen 5.9 A from its carbon
        and the ring torn open -- at an energy that looked perfectly ordinary.
        """
        import tempfile

        who = _structure_fingerprint(xyz)
        atoms = len(who)
        kept = state.get('gfn_topology')
        if kept and kept.get('who') == who and atoms:
            return Path(kept['dir'])
        _drop_gfn_topology()
        folder = tempfile.mkdtemp(prefix='delfin-gfnff-topo-')
        state['gfn_topology'] = {'dir': folder, 'atoms': atoms, 'who': who}
        # Perceived here and now, from the structure as it stood before a hand
        # was laid on it.  Left to the first calculation that needs it, the
        # perception is made from a geometry that has already been pulled
        # about -- and if the drag has gone past where a bond is still
        # recognised, the topology that is then kept for the whole session is
        # one with the bond missing.  Measured: a propane whose C-C had been
        # pulled to 2.1 A perceived that way came back at 3.57 A.
        seed = (state.get('gfn_topology_source')
                or _current_xyz())
        if seed and len(_gfn.atom_lines(seed)) == atoms:
            try:
                _gfn.relax_steps(seed, cycles=1, timeout=30.0,
                                 topology=Path(folder))
            except Exception:
                pass
        return Path(folder)

    def _drop_gfn_topology():
        """Forget the bonding: the molecule is not the same one any more."""
        kept = state.pop('gfn_topology', None)
        if kept:
            shutil.rmtree(kept.get('dir') or '', ignore_errors=True)

    def _gfn_new_generation():
        """Everything in flight belongs to the structure it was started for.

        A drag, a Hold, a Set: each makes a new structure, and every timer and
        every worker started for the one before it is now about something that
        no longer exists.  They used to be left to finish -- writing their
        geometry over the new one, or holding a flag that made the next
        relaxation skip itself, which is why switching the toggle off and on
        again finished the job that letting go should have.  One counter
        settles all of it: whatever does not belong to the current generation
        does nothing at all.
        """
        state['gfn_generation'] = int(state.get('gfn_generation', 0)) + 1
        state['gfn_settle_forced'] = False
        state['gfn_settle_rounds'] = 0
        state['gfn_settle_again'] = False
        return state['gfn_generation']

    def _gfn_uhf_now():
        """How many unpaired electrons the live relaxation should assume.

        The box says a multiplicity and xtb counts unpaired electrons, so
        M = 2S+1 becomes M - 1.  With auto M on, the box is not the answer at
        all -- Optimise scans and keeps the lowest -- and a live relaxation
        running the box's fixed M while Optimise ran a scanned one is two
        answers about two different molecules, which is why pressing Optimise
        after it moved the structure again.  Scanning every round would cost
        three runs a round, so what the last scan settled on is used instead.
        """
        if submit_gfn_autospin.value and state.get('gfn_scanned_uhf') is not None:
            return int(state['gfn_scanned_uhf'])
        return max(0, int(submit_gfn_mult.value or 1) - 1)

    def _gfn_live_is_on():
        """Whether something on screen is meant to act on a change at once."""
        return (submit_relax_btn.value and _server_method()
                and _server_binary(submit_ff_dd.value) is not None)

    def _arm_gfn_takeup(note=''):
        """Take up a change to what is held, straight away.

        Setting a value, watching nothing happen and having to press Optimise
        for it is the switch claiming to be live and not being it.  Under the
        browser's field a held value acts the moment it is set, because the
        field is already running; here nothing is running between drags, so
        the change is what starts it.
        """
        if not _gfn_live_is_on():
            return
        _gfn_new_generation()
        if _climb_owns_the_release():
            # The same change, taken up the same way, walking the other way.
            # Routed to the settle it did nothing at all while Climb to TS was
            # down: the settle stands aside for the climb, and nothing else
            # was listening.  Asked for by a hand, so it happens whether or
            # not Auto is on -- exactly as the settle it replaces here does.
            _arm_gfn_optimise(asked=True)
            return
        state['gfn_settle_note'] = note
        _arm_gfn_settle(forced=True)

    def _climb_owns_the_release():
        """Which way the optimiser walks when a release hands it a structure.

        The toggle, and not whether a climb happens to be running.  That is
        the whole of what Climb to TS is: Dynamik Opt follows the hand and
        Auto carries on when the hand lets go, and this says which way "carry
        on" means.  It is a mode the way Dynamik Opt and Auto are modes, and
        it stays down across the runs it starts.

        Read as "a climb is running" it was a one-shot.  Measured on the
        van-der-Waals complex, dragging twice in a row with Dynamik Opt and
        Auto on: the first release climbed to 2.316 A and switched the toggle
        off behind it, and the second release -- the toggle now up, so downhill
        again -- walked the structure back to 3.353 A.  Which is "ich kann es
        nicht beeinflussen": you can point it at a saddle exactly once.

        Downhill and uphill cannot both answer one release.  A settle and an
        auto-minimisation walk down, a climb walks *up* along one mode while
        walking down along every other, and armed together downhill wins --
        which is what "I drag towards the transition state, let go, and it
        falls back" was.  Nothing is lost by standing the settle down: the
        tidy-up it performs is exactly what the climb's minimised subspace
        does on every one of its steps, and the climb does it without undoing
        the one direction the user pointed at.

        The drag itself is untouched either way.  Dynamik Opt goes on relaxing
        under the hand and streaming its frames exactly as it always did; only
        where the release ends up differs.
        """
        return bool(submit_climb_btn.value)

    def _arm_gfn_settle(forced=False):
        """Tidy the structure after a release, with the method on screen.

        The same promise Settle makes under the browser's field: what reaches
        the coordinate box is a structure that has relaxed around where the
        atom was put, rather than wherever the cursor happened to stop --
        measured at 176 kcal/mol above a settled one on a real complex.

        Waits a moment first, because the coordinates of the release arrive on
        a message of their own and this has to start from them.
        """
        if forced:
            # Asked for by a hand -- a value held, a value set.  That is not a
            # tidy-up to be skipped when something else is in the air; it is
            # the answer to something the user just did, and it has to happen.
            state['gfn_settle_forced'] = True
        if _climb_owns_the_release():
            return
        if not (_server_method()
                and _server_binary(submit_ff_dd.value) is not None):
            return
        # Only a hand may arm this: a value set, a value held.  The gate above
        # has already turned away everything but a server method, and on a
        # server method Settle is not on the toolbar -- so reading its widget
        # here, as this did, meant a switch the user could not see, left on
        # from an earlier method or restored along with a structure, ran a
        # full uncapped xtb optimisation on every release with Dynamik Opt and
        # Optimise both off.  What ran was the whole minimisation rather than a
        # tidy-up: this path stopped capping its cycles when it was made to
        # match Optimise.  Going to a minimum on release is Auto's, and Auto
        # is a switch that is on the toolbar saying so.
        if not state.get('gfn_settle_forced'):
            return
        if not forced and state.get('optimize_interrupted') is not None:
            # An optimisation was interrupted by this very drag and is coming
            # back.  A whole minimisation is more than a settle, not less.
            return
        # A round that follows another round waits only long enough not to
        # spin; a release waits for its coordinates to arrive first.
        state['gfn_settle_at'] = time.monotonic() + (
            0.05 if state.get('gfn_settle_again') else _GFN_RESTART_DELAY)
        if state.get('gfn_settle_armed'):
            return
        state['gfn_settle_armed'] = True
        generation = int(state.get('gfn_generation', 0))

        def _wait():
            while True:
                left = state.get('gfn_settle_at', 0.0) - time.monotonic()
                if left <= 0:
                    break
                time.sleep(min(left, 0.05))
            state['gfn_settle_armed'] = False
            if int(state.get('gfn_generation', 0)) != generation:
                return          # the structure it was armed for is history
            schedule_ui_update(_gfn_settle_now)

        _start_background(_wait, 'The settle after the release',
                          guards={'gfn_settle_armed': False})

    def _gfn_settle_now():
        if _climb_owns_the_release():
            # Armed before the climb was switched on, and firing after.  The
            # arming guards stop one being scheduled; this stops one that
            # already was, which is the difference between the invariant
            # holding and it holding most of the time.
            return
        if state.get('gfn_follow'):
            return                      # a new drag started in the meantime
        if state.get('gfn_settle_busy'):
            # One is already running.  What has just been asked for is a
            # different question from the one it is answering, so it is asked
            # again the moment that one is done -- dropped here, a value held
            # while something was in flight did nothing at all until the
            # toggle was switched off and on again.
            state['gfn_settle_pending'] = True
            return
        if state.get('optimize_run') is not None:
            return                      # Optimise is doing this already
        generation = int(state.get('gfn_generation', 0))
        method = str(submit_ff_dd.value)
        xyz = _current_xyz()
        # Either engine, not only xtb.  Settle was armed for the PM methods
        # along with the follow, but this gate let only GFN through -- so
        # letting go of an atom under PM7 with Settle on did nothing, and said
        # nothing about doing nothing.
        if not xyz or not _server_method(method):
            return
        label = _server_label(method)
        state['gfn_settle_busy'] = True
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        constraints = list(state.get('constraints') or [])
        wet = str(submit_gfn_solvent.value or '') or None
        model = _solv_model()
        rounds = int(state.get('gfn_settle_rounds') or 0) + 1
        state['gfn_settle_rounds'] = rounds
        # One run number for the whole relaxation, not one per round.  A new
        # run number resets the player, which drops what it had not drawn yet
        # and applies the next frame outright instead of moving to it: with a
        # round every few tenths of a second that is a twitch per round, and
        # what it looks like is the structure jittering.
        if rounds == 1:
            state['gfn_run'] = _note_the_run(
                int(state.get('gfn_run', 0)) + 1, 'settle')
            state['gfn_settle_offset'] = 0
        run = int(state.get('gfn_run', 0))
        offset = int(state.get('gfn_settle_offset') or 0)
        note = str(state.get('gfn_settle_note') or '')
        _set_mol_status(
            (f'{note}: {label} is moving the structure to it...' if note
             else f'{label} is settling the structure...')
            + (f' (round {rounds})' if rounds > 1 else ''), spinner=True)
        # Only GFN-FF has a topology to keep, and asking for one makes a
        # directory -- so a PM settle does not ask.
        perceived = (_gfn_topology_dir(xyz)
                     if _gfn.is_gfn_method(method) else None)
        # Auto M, and no scan has happened yet: this run does the scanning, so
        # that it and Optimise are asking about the same molecule.  Only xtb
        # can be asked to scan; under PM the multiplicity on screen is used.
        scanning = bool(submit_gfn_autospin.value
                        and _gfn.is_gfn_method(method)
                        and state.get('gfn_scanned_uhf') is None)

        def _settle_stopped():
            # A hand on an atom takes over from a relaxation of the whole
            # thing; Optimise is the same calculation asked for by hand and
            # supersedes it outright -- two xtb runs writing the same
            # coordinate box is how the first press came to look like it had
            # only worked out an energy.  And switching off means what it says,
            # whichever of the two switches was keeping this alive.
            #
            # A climb pressed while this is running is the same collision seen
            # from the other side, and worse: they walk opposite ways, so the
            # two would take the structure apart between them.  The press is
            # the later of the two things the user asked for.
            return bool(state.get('gfn_follow')
                        or state.get('optimize_run') is not None
                        or state.get('climb_run') is not None
                        or not (submit_settle_btn.value or _gfn_live_is_on()))

        def _push(frames):
            walked = list(frames)
            trail = walked[-8:]
            state['gfn_settle_walked'] = len(walked)
            text = _frame_payload(
                run, **{'from': offset + len(walked) - len(trail),
                        'follow': 1, 'frames': trail})
            # A settle outlives the release that started it, so it is the
            # writer most likely to still be answering after the user has
            # moved on to something else.
            schedule_ui_update(
                lambda t=text, r=run: setattr(
                    submit_gfn_frame, 'value', t)
                if _frame_run_is_current(r) else None)

        def _work():
            began = time.perf_counter()
            if _mopac.is_mopac_method(method):
                # No topology to carry, but the held values and the solvent
                # both are: a settle in water after a drag in water is the
                # same question the drag was asking, and a value held while
                # the drag ran is still held when the hand lets go.
                outcome = _mopac.optimize_with_mopac(
                    xyz, method, charge=charge, uhf=uhf,
                    max_steps=None, timeout=None, on_frames=_push,
                    constraints=constraints,
                    should_stop=_settle_stopped, solvent=wet)
            elif scanning:
                # Auto M with nothing scanned yet.  Optimise would scan and
                # keep the lowest; running the box's M here instead makes the
                # two answer about different molecules, and pressing Optimise
                # afterwards then moves the structure for no visible reason.
                # It costs three runs, once -- after that the answer is known.
                outcome = _gfn.optimize_autospin(
                    xyz, method, charge=charge, constraints=constraints,
                    on_frames=_push, topology=perceived, timeout=None,
                    solvent=wet, solvation_model=model,
                    should_stop=_settle_stopped,
                )
                if outcome.get('uhf') is not None:
                    state['gfn_scanned_uhf'] = int(outcome['uhf'])
            else:
                # No cycle cap and no clock: this is the ordinary
                # optimisation, run on the frame that is on screen now, on the
                # same terms Optimise runs on.  Chopping it into rounds bought
                # nothing -- the geometry between two rounds is not a place
                # anyone wants to stop at -- and cost a stutter at every
                # boundary and an early ending whenever a round moved little.
                outcome = _gfn.optimize_with_gfn(
                    xyz, method, charge=charge, uhf=uhf,
                    max_steps=None, timeout=None,
                    constraints=constraints, on_frames=_push, solvent=wet,
                    solvation_model=model,
                    topology=perceived, should_stop=_settle_stopped,
                )

            def _done():
                state['gfn_settle_busy'] = False
                if state.pop('gfn_settle_pending', False):
                    # Something was asked for while this was running.
                    schedule_ui_update(lambda: _arm_gfn_settle(forced=True))
                if (int(state.get('gfn_generation', 0)) != generation
                        or state.get('optimize_run') is not None):
                    # The structure moved on while this was running, or
                    # Optimise took over.  Either way this geometry is about a
                    # molecule that is no longer the current one, and writing
                    # it would be the past reaching into the present.
                    return
                state['gfn_settle_again'] = False
                # Where the next round's frames carry on from.
                state['gfn_settle_offset'] = offset + int(
                    state.pop('gfn_settle_walked', 0) or 0)
                if not outcome.get('ok'):
                    state['gfn_settle_forced'] = False
                    state['gfn_settle_rounds'] = 0
                    _set_mol_status(
                        f'{label} could not settle it: '
                        f'{outcome.get("status") or "it did not run"}')
                    return
                # The release answers to the ceiling as well.
                #
                # A settle is a relaxation and a relaxation goes downhill, so
                # this ordinarily has nothing to refuse -- which is exactly why
                # it is worth asking rather than assuming: a value the user is
                # holding is restored on every step of it, and restoring one
                # is uphill.  Anything a settle produces is what outlives the
                # drag, so it is held to the same ceiling the drag was.
                priced = _settle_price(outcome, constraints)
                anchor, ceiling = _thermal_budget()
                over = (None if not (submit_thermal_btn.value
                                     and anchor is not None
                                     and priced is not None)
                        else (float(priced) - float(anchor))
                        * _HARTREE_TO_KCAL)
                if over is not None and over > ceiling:
                    state['gfn_settle_forced'] = False
                    state['gfn_settle_rounds'] = 0
                    _set_mol_status(
                        f'{label} settled to a structure that costs '
                        f'{over:+.1f} kcal/mol, past the {ceiling:.1f} this '
                        f'one has at {float(submit_temperature.value):g} K, '
                        'so it has been left as it was. '
                        f'{_thermal_wait(over, submit_temperature.value)} '
                        f'{_thermal_wants(over)}')
                    return
                lines = [line for line in outcome['xyz'].splitlines()[2:]
                         if line.strip()]
                if lines:
                    # The playback has drawn this already; the box is what Copy
                    # and Submit read, and it has to be true whether or not a
                    # frame happened to land.
                    _write_coords(xyz_document(lines, f'Settled with {label}'),
                                  drawn=True, run=run)
                if over is not None:
                    state['thermal_good'] = outcome['xyz']
                # Not converged and the switch is still down: keep going.  That
                # is what makes this a relaxation rather than a single push.
                # It ends three ways -- converged, standing still, or out of
                # rounds -- because held values that cannot all be met at once
                # never converge, and a relaxation that will not end is a
                # process per round for as long as the switch is down.
                moved = _gfn.largest_shift(xyz, outcome['xyz'])
                if (not outcome.get('converged') and _gfn_live_is_on()
                        and not state.get('gfn_follow')
                        and rounds < _GFN_SETTLE_ROUNDS
                        and moved > _GFN_SETTLE_STILL):
                    state['gfn_settle_again'] = True
                    state['gfn_settle_forced'] = True
                    _arm_gfn_settle()
                    return
                state['gfn_settle_forced'] = False
                state['gfn_settle_note'] = ''
                took = time.perf_counter() - began
                if outcome.get('converged'):
                    said = (f'{label} relaxed the structure: converged after '
                            f'{rounds} round(s).')
                elif not _gfn_live_is_on() or state.get('gfn_follow'):
                    said = f'{label} settled the structure in {took:.1f} s.'
                else:
                    # Said out loud, because a structure that has stopped
                    # improving and one that is finished look identical, and
                    # only one of them is worth pressing Optimise on.
                    why = ('it is no longer moving' if moved <= _GFN_SETTLE_STILL
                           else f'{rounds} rounds is as far as this goes')
                    said = (
                        f'{label} stopped without converging: {why}. '
                        + ('Held values that cannot all be met at once are the '
                           'usual reason -- the list beside the toolbar is '
                           'what it is trying to satisfy.'
                           if state.get('constraints') else
                           'Optimise will run it without a round limit.'))
                state['gfn_settle_rounds'] = 0
                state['gfn_last_status'] = said
                _set_mol_status(said)

            schedule_ui_update(_done)

        _start_background(_work, 'The settle',
                          guards={'gfn_settle_busy': False})

    def _enable_live_forcefield():
        """Assign UFF parameters for the geometry now in the viewer.

        Runs once, when the toggle is switched on -- never during a drag. The
        browser relaxes from these parameters alone; a round trip per frame
        would cap the drag at about 13 Hz.

        Refused outright while a GFN method is chosen.  A dozen handlers call
        this -- Hold, a polyhedron, a hybridisation, a bond edit -- and any one
        of them installing UFF parameters put a UFF relaxation under a molecule
        whose method box said GFN: it settled on release, and the geometry that
        reached the coordinate box was UFF's. The method that is chosen is the
        method that acts, and nothing else may touch the structure.
        """
        if _server_method():
            _stop_browser_field()
            return
        xyz = _current_xyz()
        if not xyz:
            _set_mol_status('Load a structure before enabling Relax.')
            submit_relax_btn.value = False
            return
        try:
            from .molecule_forcefield import export_forcefield_terms
            polyhedron = None
            if state.get('poly_applied') and state.get('poly_metal') is not None:
                polyhedron = {
                    'metal': state['poly_metal'],
                    'geometry': state['poly_applied'],
                    # None means: work it out from where the ligands are now.
                    'assignment': state.get('poly_assignment'),
                }
            payload = export_forcefield_terms(
                xyz,
                perceived=_perception_for(xyz),
                # The live field is what the browser relaxes with on every
                # frame, so it is always one of the two that live there.
                method=_live_ff_method(),
                polyhedron=polyhedron,
                restraints=[
                    c for c in (state.get('constraints') or [])
                    if c.get('mode', 'pull') == 'pull'
                ],
            )
        except Exception as exc:
            _set_mol_status(f'Force field unavailable: {exc}')
            submit_relax_btn.value = False
            return
        if not payload.get('ok'):
            _set_mol_status('Force field could not be assigned for this structure.')
            submit_relax_btn.value = False
            return
        _ensure_manip_bootstrap()
        _ensure_ff_bootstrap()
        _push_bond_orders()
        # The resume flag is set here, in the same script that hands over the
        # parameters, and not earlier: setting it before the re-render meant
        # it could be consumed against the viewer that was going away, and the
        # relaxation came back stuck until the toggle was cycled by hand.
        resume = 'true' if submit_relax_btn.value else 'false'
        _run_manip_js(
            f'window.__delfinResumeAutoOpt = {resume};'
            'if(window.__delfinSubmitManip){'
            'window.__delfinSubmitManip.setForceField('
            f'{json.dumps(submit_scope_id)},{json.dumps(payload)});'
            'window.__delfinSubmitManip.setOptimizerStrength('
            f'{json.dumps(submit_scope_id)},{int(submit_strength_slider.value)});'
            'window.__delfinSubmitManip.setPullStrength('
            f'{json.dumps(submit_scope_id)},{_hand_share()},'
            f'{json.dumps(_pull_most())});'
            # Re-applied with the rest, so a reload keeps the feel the user set.
            'window.__delfinSubmitManip.setDragSensitivity('
            f'{json.dumps(submit_scope_id)},{float(submit_sens_slider.value)});'
            'window.__delfinSubmitManip.setSettleOnRelease('
            f'{json.dumps(submit_scope_id)},'
            f'{"true" if submit_settle_btn.value else "false"});'
            'window.__delfinSubmitManip.setFixedInternals('
            f'{json.dumps(submit_scope_id)},'
            + json.dumps([
                {'kind': c['kind'], 'atoms': c['atoms'], 'value': c['value']}
                for c in (state.get('constraints') or [])
                if c.get('mode') == 'fix'
            ])
            + ');'
            '}'
        )
        # Terms derived from the input geometry rather than real UFF typing --
        # the transition-metal case -- are worth saying out loud, and they
        # belong under the structure they describe rather than in the preview's
        # status line, which conversion messages keep overwriting.
        _set_ff_notes(payload.get('warnings') or [])

    def _stop_browser_field():
        """Take the browser's own field off the molecule.

        It has to be said to the page, not merely stopped being asked for: the
        relaxation runs on its own animation frames and keeps running until it
        is told, and the parameters stay installed until they are cleared --
        which is how a GFN method could be chosen while UFF went on relaxing
        underneath it.
        """
        _ensure_manip_bootstrap()
        _set_ff_notes([])
        _run_manip_js(
            'if(window.__delfinSubmitManip){'
            'window.__delfinSubmitManip.stopAutoOptimize('
            f'{json.dumps(submit_scope_id)});'
            'window.__delfinSubmitManip.setForceField('
            f'{json.dumps(submit_scope_id)},null);'
            '}'
        )

    def on_submit_relax_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_relax_btn.value)
        if active:
            # Where the continuous relaxation started, kept as one step. It
            # runs for as long as it is on and moves the structure the whole
            # time; one entry for the whole run is what a user means by "back
            # to before I switched it on".
            _remember('switching the continuous relaxation on')
        submit_relax_btn.button_style = 'info' if active else ''
        if _server_method():
            # There is no GFN engine in the browser to run per frame, so this
            # switch means something else here: while it is on, the molecule
            # follows the atom being dragged -- one short xtb run per push,
            # and nothing at all when nothing is being dragged.  A continuous
            # loop of processes on a shared machine is what this deliberately
            # is not.
            _end_gfn_follow()
            # Whether it goes on or off, the browser's field has no business
            # here.  Left running from a UFF session it went on relaxing under
            # a molecule whose method box said GFN.
            _stop_browser_field()
            label = _server_label(submit_ff_dd.value)
            if not active:
                state['gfn_settle_forced'] = False
                state['gfn_settle_rounds'] = 0
                # The budget goes with it, because it cannot price a drag
                # this switch is not answering.  The page only reports a drag
                # while it is down; without those messages nothing runs, the
                # ceiling has nothing to compare against, and it would sit
                # there lit up refusing nothing at all -- which is worse than
                # being off, because it is off and says it is on.
                if submit_thermal_btn.value:
                    submit_thermal_btn.value = False
                    _set_mol_status(
                        'The structure is no longer being relaxed, so the '
                        'thermal budget has nothing to measure a drag with '
                        'and has gone off with it. Switch the relaxation back '
                        'on to have changes priced again.')
                    return
                _set_mol_status('The structure is no longer being relaxed.')
                return
            if _server_binary(submit_ff_dd.value) is None:
                _set_mol_status(f'{label} needs a program that was not found.')
                submit_relax_btn.value = False
                return
            if not submit_manip_btn.value:
                submit_manip_btn.value = True  # dragging is what it is for
            _ensure_manip_bootstrap()
            _install_gfn_frame_watcher()
            said = [
                f'Drag an atom and the rest of the molecule follows it with '
                f'{label}. Letting go leaves it where you put it -- Optimise '
                'is what takes it downhill, and Settle asks for that on every '
                'release.'
            ]
            if submit_ff_dd.value != 'gfnff':
                said.append(f'{label} is the slow one -- if it drags heavily, '
                            'GFN-FF answers about twenty times faster.')
            # A solvent is free to drag in, except for the one that is not.
            # Measured on a benzophenone, one follow step: GFN2 167 ms in
            # vacuum, 117 with ALPB, 168 with GBSA -- and 1020 with ddCOSMO.
            # Said here rather than left to be discovered, because what six
            # times the cost per step looks like is a drag that is broken.
            wet_now = str(submit_gfn_solvent.value or '')
            if wet_now and _solv_model() == 'ddcosmo':
                said.append('ddCOSMO costs about six times what the other '
                            'models do per step -- a second each on 24 atoms '
                            '-- so this will move like a slideshow. ALPB is '
                            'the one to drag in.')
            elif wet_now:
                said.append(f'Following in {_solvents.label_of(wet_now)} '
                            f'({_solvents.model_label(_solv_model())}).')
            _set_mol_status(' '.join(said))
            return
        if not active:
            _ensure_manip_bootstrap()
            _set_ff_notes([])
            _run_manip_js(
                'if(window.__delfinSubmitManip){'
                'window.__delfinSubmitManip.stopAutoOptimize('
                f'{json.dumps(submit_scope_id)});'
                'window.__delfinSubmitManip.setForceField('
                f'{json.dumps(submit_scope_id)},null);'
                '}'
            )
            return
        if not submit_manip_btn.value:
            submit_manip_btn.value = True   # dragging is what it is there for
        _enable_live_forcefield()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.startAutoOptimize('
            f'{json.dumps(submit_scope_id)});'
        )

    def _frame_as_xyz(source, walked, shown, comment):
        """One frame of a walked path, as a structure, or ``None``.

        *shown* is the page's count of the frames it has drawn, so the frame
        the picture is standing on is the one before it.  The arithmetic is
        here and nowhere else: it was written out three times -- the
        optimisation's Stop, the optimisation cut off by a hand, and now the
        climb -- and three copies of an off-by-one is three chances to hand
        back a geometry nobody was looking at.
        """
        walked = list(walked or [])
        if not isinstance(shown, int) or not (0 < shown <= len(walked)):
            return None
        symbols = [line.split()[0] for line in _gfn.atom_lines(source or '')]
        frame = walked[shown - 1]
        if not symbols or len(symbols) * 3 != len(frame):
            return None
        return xyz_document(
            [xyz_line(symbols[i], frame[3 * i], frame[3 * i + 1],
                      frame[3 * i + 2]) for i in range(len(symbols))],
            comment)

    def _the_picture_stopped_here(run_id, source, walked, comment):
        """Put a stopped run's path down, to be cut where the picture stands.

        A Stop means the frame on screen, and the only thing that knows which
        frame that is is the page: it is the page that draws them.  It says so
        on the channel it says everything else on -- "stopped at frame 24" --
        and that message and the worker noticing it was stopped are two
        answers racing each other.  A climb takes a step every ten
        milliseconds and hears the switch inside one of them, so it finishes
        long before a round trip to the browser and back; an xtb round is
        seconds, so it usually loses the same race instead.  Neither order may
        decide what the user is left with, so the path is put down here and
        whichever of the two arrives second lands it.

        *run_id* is kept for nothing but the record: the number has moved on,
        which is exactly why the ordinary write at the end of a run is refused
        and this one is not.  What is written is the geometry the user is
        looking at, not one the run computed past it.
        """
        state['gfn_stopped_path'] = {
            'run': int(run_id or 0), 'source': source or '',
            'comment': comment,
            'frames': [list(one) for one in (walked or [])],
        }
        return _land_the_stopped_frame()

    def _land_the_stopped_frame():
        """Write the frame the picture stopped on, once both halves are in.

        Called from the worker that was stopped and from the page's report,
        and it does its work for whichever of the two is second.  Marked as
        already drawn, because it is: the picture is standing on that frame,
        and rebuilding the viewer around it would tear down a playback that
        has correctly stopped there.
        """
        held = state.get('gfn_stopped_path')
        if not held:
            return False
        text = _frame_as_xyz(held.get('source'), held.get('frames'),
                             _the_shown_frame_of(held.get('run')),
                             held.get('comment') or
                             'stopped at the frame on screen')
        if text is None:
            return False
        state.pop('gfn_stopped_path', None)
        return _write_coords(text, drawn=True)

    #: How long to wait after the last change before starting again.  Letting
    #: go of an atom arrives as a burst -- the release, then the coordinates,
    #: sometimes a settled version behind them -- and starting on the first of
    #: those would launch an xtb for each one.
    _GFN_RESTART_DELAY = 0.35

    def _keep_the_shown_frame(frame, walk):
        """Note which frame the picture stands on, and of whose walk.

        The page counts the frames it has drawn since the run it is playing
        began, so the number means nothing without the run it is counted
        along.  It used to travel alone, and after a scan that is a number
        indexing the scan's trajectory: a minimisation started a moment later
        was handed "frame 29" and cut its own path at its 29th point, which is
        a geometry nobody had seen and nowhere near the one on screen.
        """
        try:
            state['gfn_shown_frame'] = int(str(frame).strip())
            state['gfn_shown_run'] = int(str(walk).strip())
        except (TypeError, ValueError):
            # A page from before the run travelled with the number.  Refused
            # rather than guessed at: an index into the wrong path is worse
            # than no index, because the wrong one is acted on.
            state.pop('gfn_shown_frame', None)
            state.pop('gfn_shown_run', None)

    def _the_shown_frame_of(run):
        """The page's count, but only where it counts along *this* walk.

        ``None`` means the picture is not standing on any frame of this run --
        the number was reported for another one, or nothing of this one has
        been drawn -- and then the picture is whatever was on screen before
        the run started, which is what the coordinate box already holds.
        Answering with the run's own last geometry instead is how a structure
        came to jump to somewhere nobody had chosen.

        The arithmetic that turns a count into a frame is :func:`_frame_as_xyz`
        and stays there; this is the other question, and it is asked
        separately because the two go wrong in different ways -- one by a
        frame, the other by a whole trajectory.
        """
        if state.get('gfn_shown_run') != int(run or 0):
            return None
        return state.get('gfn_shown_frame')

    def _forget_the_shown_frame():
        """The number belongs to a run that is over, so it is not kept."""
        state.pop('gfn_shown_frame', None)
        state.pop('gfn_shown_run', None)

    def _interrupt_gfn():
        """End the running optimiser because the structure under it changed.

        Whichever of the two it is.  xtb is walking a geometry that stopped
        existing the moment an atom was moved, so the run is ended rather than
        raced.  It is not ended the way the switch ends it: nothing has been
        stopped from where the user is standing, and the frame they were shown
        is not a result to keep -- the walk is about to start again from what
        they have made.

        The climb used to be the one thing a hand could not interrupt.  It
        stood still instead, in a loop of its own, and was handed the released
        structure by a second path with its own hand-over rules -- which is
        where "es steppt auf einer alten Geometrie" and "es kaempft mit
        Dynamik Opt um das Loslassen" both came from.  Stopping it costs
        nothing that standing still did not already cost: what it resumes from
        is the structure the hand made either way, and the Hessian is
        recomputed either way, because a Bofill update repairs a Hessian one
        step at a time and a hand moves further in one gesture than a climb
        does in twenty -- carried across a drag it still reaches the same
        saddle, in 62 steps against 15.
        """
        token = state.get('optimize_run')
        climbing = state.get('climb_run') is not None
        if token is None and not climbing:
            return False
        if climbing:
            # The loop reads this and stops.  Marked as interrupted rather
            # than merely stopped, because that is what tells the release it
            # has something to bring back -- and marked by token as well,
            # which is what tells the run itself that it was cut off by a hand
            # rather than switched off.  A run switched off keeps what it
            # reached, the way a stopped Optimise does; one cut off by a hand
            # keeps nothing, because the user has made a structure since.
            state['climb_cut'] = state['climb_run']
            state['climb_run'] = None
            state['climb_interrupted'] = True
        if token is not None:
            state['optimize_run'] = None
            state['optimize_interrupted'] = token
        # No halt report: "stopped at frame 12" belongs to the switch.
        state['gfn_halt_sent'] = True
        # And no landing either.  A run a hand cut off has its own answer --
        # the geometry the user has just made -- and a path waiting to be cut
        # at a frame number would be written over it.
        state.pop('gfn_stopped_path', None)
        # A run number the page has never seen, carrying nothing.  It resets
        # the player, so the frames of the abandoned run cannot play out over
        # the geometry the user has just made.
        blank = _note_the_run(int(state.get('gfn_run', 0)) + 1, 'abandoned')
        state['gfn_run'] = blank
        # Said out loud, because a new run and an abandoned one look the
        # same to the page and want opposite things: a run that ended
        # cleanly should be landed on its last frame, and one the user has
        # just cut off must not be drawn at all -- those frames are about a
        # structure that no longer exists.
        submit_gfn_frame.value = _frame_payload(
            blank, **{'frames': [], 'abandoned': 1})
        return True

    def _stand_down_after_interrupt(note=None):
        """The run the change interrupted is not coming back by itself.

        The switch goes back up with it. Left lit over a structure nobody is
        minimising, it says a calculation is running that is not.

        Climb to TS is not one of those switches. It is a mode, the way
        Dynamik Opt and Auto are: it says which way the next release walks,
        and lifting it here would answer the next drag with a minimisation
        the user never asked for. So the climb's mark is cleared and its
        toggle is left exactly as the user set it -- which also means the
        sentence differs, because there is no switch left for the user to
        press: it is already down.

        *note* is what to say instead of the usual sentence, for the caller
        that is standing the run down because something else has taken the
        release rather than because nothing has.
        """
        state.pop('optimize_interrupted', None)
        climbing = bool(state.pop('climb_interrupted', False))
        for button in (submit_optimize_btn, submit_optimize_all_btn):
            if button.value:
                button.value = False
        _set_mol_status(note or (
            'Stopped where your change left it. '
            + ('Climb to TS is still down, so turning Auto on will climb from '
               'wherever you leave it.' if climbing else
               'Move what else you want to, then press Optimize to go down '
               'to a minimum.')))

    def _the_hand_interrupted():
        """What the hand cut off and is owed a restart, or ``''``.

        Two optimisers, one question. Asked in one place so that a path that
        knows about one of them cannot quietly forget the other, which is how
        the climb came to be the only run a drag could not interrupt.
        """
        if state.get('climb_interrupted'):
            return 'climb'
        if state.get('optimize_interrupted') is not None:
            return 'optimise'
        return ''

    def _arm_gfn_restart():
        """Start the walk again, once the changing has stopped."""
        cut = _the_hand_interrupted()
        if not cut:
            return
        if cut == 'optimise' and _climb_owns_the_release():
            # The release walks up and the run that was cut off walks down.
            # It is not left armed behind the climb -- armed, it would come
            # back the moment the climb ended and undo what it had just found.
            _stand_down_after_interrupt(
                'The climb has this release. Switch Climb to TS off to go '
                'down to a minimum instead.')
            return
        if not submit_auto_btn.value:
            # Auto is off, so nothing carries on by itself. The guard is here
            # rather than at the release, because a drag that moved something
            # sends its coordinates first and asks for the restart from there
            # -- gating only the release would have let that path through and
            # the switch would have worked for some drags and not others.
            _stand_down_after_interrupt()
            return
        state['gfn_restart_at'] = time.monotonic() + _GFN_RESTART_DELAY
        if state.get('gfn_restart_armed'):
            return                      # already waiting; it will wait longer
        state['gfn_restart_armed'] = True
        # Which structure this is waiting for.  A third of a second is long
        # enough to press Undo in, and the wait went on regardless: the
        # minimisation woke up and ran on the geometry the user had just put
        # back, so Undo looked as though it had done nothing.  Every timer in
        # here belongs to the structure it was started for, and the generation
        # counter is how the rest of them say so.
        generation = state.get('gfn_generation')

        def _wait():
            while True:
                left = state.get('gfn_restart_at', 0.0) - time.monotonic()
                if left <= 0:
                    break
                time.sleep(min(left, 0.05))
            state['gfn_restart_armed'] = False
            if state.get('gfn_generation') != generation:
                return
            schedule_ui_update(_restart_gfn)

        _start_background(_wait, 'The optimiser waiting to restart',
                          guards={'gfn_restart_armed': False})

    def _arm_gfn_optimise(asked=False):
        """Carry on to wherever the toggles say, once the changing has stopped.

        The other half of Auto: a drag with no run behind it has nothing to
        resume, so one is started. Which way it walks is Climb to TS and
        nothing else -- down to a minimum with it up, up to a saddle with it
        down. Same wait either way, and the wait is pushed out again by every
        release, so moving three atoms one after another is one walk at the
        end of it rather than three.

        *asked* is a hand -- a value set, a value held -- rather than a
        release. That is not something Auto may decline: it is the answer to
        something the user just did, the same way the settle it stands in for
        runs whether Auto is on or not.
        """
        if not (_server_method()
                and _server_binary(submit_ff_dd.value) is not None):
            return
        if asked:
            state['gfn_optimise_asked'] = True
        state['gfn_minimise_at'] = time.monotonic() + _GFN_RESTART_DELAY
        if state.get('gfn_minimise_armed'):
            return
        state['gfn_minimise_armed'] = True
        # The structure this wait belongs to; see _arm_gfn_restart for the
        # third of a second an Undo fits into.
        generation = state.get('gfn_generation')

        def _wait():
            while True:
                left = state.get('gfn_minimise_at', 0.0) - time.monotonic()
                if left <= 0:
                    break
                time.sleep(min(left, 0.05))
            state['gfn_minimise_armed'] = False
            if state.get('gfn_generation') != generation:
                state.pop('gfn_optimise_asked', None)
                return
            schedule_ui_update(_optimise_now)

        _start_background(_wait, 'The optimiser waiting to start',
                          guards={'gfn_minimise_armed': False})

    def _optimise_now():
        """Start whichever of the two the toggles chose."""
        asked = bool(state.pop('gfn_optimise_asked', False))
        if not (asked or submit_auto_btn.value):
            return                      # switched off while it was waiting
        if _climb_owns_the_release():
            if state.get('climb_run') is not None:
                return                  # one is already walking
            _climb_now()
            return
        if state.get('optimize_run') is not None:
            return                      # one is already running
        if submit_optimize_btn.value or submit_optimize_all_btn.value:
            return                      # a switch is already down
        # The switch, not the handler: it has to be seen to be on for as long
        # as it runs, and it is what the user presses to stop it again.
        submit_optimize_btn.value = True

    def _after_release():
        """What letting go of an atom leads to, and the switches that decide.

        Auto on: carry on, whether or not a run was going when the atom was
        picked up. That used to be the difference between the same gesture
        finishing the structure and leaving it strained -- a drag during a run
        interrupted it and the run came back, a drag with no run behind it got
        Settle's short tidy and nothing else.

        Auto off: it stops where the hand left it. Move something else, and
        press Optimize -- or Climb to TS -- when the structure is what you
        meant.

        Only while the molecule is following the hand. Dragging with Dynamik
        Opt off is placing an atom where you want it, and starting a walk on
        top of that would take it off the place you just put it -- which is
        what Settle is for, in the small.

        Climb to TS changes exactly one thing about all of that: which way the
        walk goes. Every line here is read the same either way, which is the
        point -- the two used to be separate paths that agreed about some of
        this and not the rest.
        """
        auto = bool(submit_auto_btn.value)
        if _the_hand_interrupted():
            # Comes back, or stands down; either way _arm_gfn_restart decides.
            _arm_gfn_restart()
            if auto:
                return                  # a whole walk is more than a settle
            _arm_gfn_settle()
            return
        if auto and submit_relax_btn.value:
            _arm_gfn_optimise()
            return
        _arm_gfn_settle()

    def _restart_gfn():
        if state.pop('climb_interrupted', False):
            # The hand cut a climb off, and this is it coming back: from the
            # structure the hand made, aimed along the way it was made. An
            # optimisation cut off by the same hand does not also come back --
            # they are two answers to one release, and the toggle chose.
            state.pop('optimize_interrupted', None)
            _climb_now()
            return
        if state.pop('optimize_interrupted', None) is None:
            return
        every_frame = bool(state.get('optimize_every_frame'))
        button = submit_optimize_all_btn if every_frame else submit_optimize_btn
        if not button.value:
            return                      # switched off while it was waiting
        state['gfn_restarting'] = True
        on_submit_optimize(None, every_frame=every_frame)

    def _optimise_carries_on(*, converged, moved, rounds, failed, every_frame,
                             still):
        """Whether one more xtb run belongs to the press that is running.

        A run is one xtb ``--opt`` and its optimiser has a cycle limit of its
        own; at the limit it hands back the geometry it reached and says it did
        not converge.  Another run from there is a fresh cycle budget and a
        fresh optimiser history, which is why pressing again always got
        further.  The pressing is done here instead.

        *failed* is a run that produced no geometry, and it ends the press.
        "Stopped before converging" is NOT that: it is a geometry, and it is
        the case these rounds exist for.  Counted among the failures -- as it
        was, because it is written into the same list the status line reads --
        this returned False for the only situation it was written for, and the
        user went on pressing.  That is why the two are counted apart, and why
        this is a function with a test that drives it rather than a condition
        buried in a closure.

        It ends the three ways a settle has always ended: converged, no longer
        moving, or out of rounds.  The last two are not decoration -- held
        values that cannot all be met at once never converge, and a run that
        will not end is one process per round for as long as the switch is down.

        *still* counts how many rounds in a row have moved next to nothing,
        this one included, and one of them is not enough to call it finished.
        The whole reason these rounds exist is that a new process starts with a
        fresh cycle budget AND a fresh optimiser history -- so a run that ran
        out of cycles in a flat stretch and moved four thousandths of an
        angstrom is exactly the run whose successor can still take a real step.
        Stopping on the first quiet one made "no longer moving" mean "did not
        move much just then", and the user pressed Optimise again and watched
        it move: the very thing these rounds were written to stop.  Two in a
        row is the proof, and it costs one xtb run at most.  Held values that
        cannot all be met at once oscillate, so they reach two quickly.
        """
        if failed or converged:
            return False
        if rounds >= _OPTIMISE_ROUNDS:
            return False
        if moved <= _OPTIMISE_STILL and still >= _OPTIMISE_STILL_ROUNDS:
            return False
        # And only while the switch the user pressed is still down.
        switch = submit_optimize_all_btn if every_frame else submit_optimize_btn
        return bool(switch.value)

    def on_submit_optimize_all(change=None):
        on_submit_optimize(change, every_frame=True)

    def on_submit_optimize(change=None, every_frame=False):
        """A switch, not a push: on starts it, off stops it, and it turns
        itself off when the optimisation has converged or failed.

        *every_frame* is the difference between the two buttons: Optimize takes
        the frame on screen, all takes the whole set -- one after another,
        because a login node is shared.
        """
        button = submit_optimize_all_btn if every_frame else submit_optimize_btn
        if isinstance(change, dict) and change.get('name') == 'value':
            if not button.value:
                running = state.get('optimize_run')
                state['optimize_run'] = None      # off: the run ends itself
                if running is not None:
                    # Pressed while it was still walking, so this is a Stop
                    # and not the switch coming up behind a run that
                    # finished.  The run number moves on, and every write
                    # the stopped run still has in hand is refused: it had
                    # frames computed and not yet sent when it was told to
                    # stop, and writing them afterwards put the player back
                    # where that window began.  Measured, a Stop at frame
                    # 69 then drew 57, 59, 61, 63, 65, 67, 69 -- the tail a
                    # second time, and time running backwards.
                    #
                    # Only for a real Stop.  A run that converged sets this
                    # to None itself and then lifts the switch, and moving
                    # the number there would refuse its own final write --
                    # which is the one frame that has to land.
                    #
                    # Through the one function that hands out run numbers, so
                    # a Stop here and a Stop on the climb move the counter the
                    # same way and clear the same things behind them.
                    _claim_the_frame_run()
                return
        elif not button.value:
            return
        if submit_optimize_btn.value and submit_optimize_all_btn.value:
            # One run at a time, whichever was asked for second stands down.
            other = submit_optimize_btn if every_frame else submit_optimize_all_btn
            other.value = False
        """Minimise every frame that is loaded, not just the one on screen.

        The Submit tab can hold a whole set at once -- generated isomers, or
        the frames of a batch -- and any of them can end up submitted, so
        optimising only the visible one would leave the rest untouched.

        The geometries from before the run are kept so Undo can put them back:
        the browser's own undo stack cannot, because the results arrive from
        Python and re-render the viewer.
        """
        frames = list(list_structures() or []) if every_frame else []
        single = _current_xyz()
        if not frames and not single:
            _set_mol_status('Load a structure before optimising.')
            return
        method = submit_ff_dd.value
        gfn = _gfn.is_gfn_method(method)
        pm = _mopac.is_mopac_method(method)
        label = _server_label(method)
        charge = int(submit_gfn_charge.value or 0)
        # xtb counts unpaired electrons, not multiplicity: M = 2S+1.
        uhf = max(0, int(submit_gfn_mult.value or 1) - 1)
        autospin = bool(submit_gfn_autospin.value)
        count = len(frames) or 1
        # Which switch is running, so a restart presses the same one.
        state['optimize_every_frame'] = bool(every_frame)
        again = bool(state.pop('gfn_restarting', False))
        # A round of the same press, rather than a press.  It keeps the undo
        # point and the round count of the press it belongs to, and it does not
        # announce itself: the line already says which round is running.
        carrying_on = bool(state.pop('optimize_carrying_on', False))
        if not carrying_on:
            state['optimize_rounds'] = 0
            state['optimize_still'] = 0
            _set_mol_status(
                (f'Moved while it ran; {label} starts again from the structure '
                 'you made...' if again else
                 f'Optimising {count} frame(s) with {label}...'), spinner=True,
            )
        if gfn:
            # One call, not two: run_js clears its output before displaying,
            # so a bootstrap followed immediately by a watcher is a bootstrap
            # the page may never have run.
            _ensure_manip_bootstrap()
            schedule_ui_update(_install_gfn_frame_watcher)
        played = [False]
        state['gfn_energy'] = None
        state['gfn_held'] = None
        state['gfn_halt_sent'] = False
        run_id = _claim_the_frame_run()
        # Which frame the picture stood on belongs to the run it was reported
        # for, and the page only reports it when a hand arrives or the switch
        # goes up.  Kept across runs, a number left over from an earlier grab
        # is a plausible index into this run's path -- so an edit that
        # interrupts this one would cut it at a frame nobody ever saw, from a
        # trajectory that no longer exists.  The run now travels with the
        # number and is checked, so this is belt as well as braces.
        _forget_the_shown_frame()

        def _push_frames(frames, final=False):
            """Hand the path over while xtb is still walking it.

            Through the one writer both optimisers use: see
            :func:`_stream_frames` for why the window is shaped the way it is.
            *played* is what the write at the end reads to decide whether the
            picture already has the geometry the box is about to be given.
            """
            played[0] = True
            _stream_frames(run_id, frames, final=final)

        # Where the optimisation started from, as one step: pressing Undo
        # after it comes back should return the geometry that was handed to
        # it, not one of the frames along the way.  A round that carries the
        # same press on does not take a second one: Undo after it should reach
        # the geometry the user pressed on, not the round before last.
        if not carrying_on:
            _remember('an optimisation')

        token = object()
        state['optimize_run'] = token

        def _stopped():
            """Whether this run is still the one, and the page told if not.

            The same two lines the climb's loop reads, through the same halt:
            one Stop, one thing it means.
            """
            halted = state.get('optimize_run') is not token
            if halted:
                _halt_the_frames(run_id)
            return halted

        # The run this one replaces, taken before it is overwritten below.
        earlier = state.get('optimize_thread')
        # What the user is holding.  Set and Hold mean the same to xtb as they
        # mean to the browser's field -- a pull negotiates, a fix is met -- so
        # a value held on screen is held in the optimisation too, rather than
        # being quietly given up the moment GFN is chosen.
        held = list(state.get('constraints') or [])
        if every_frame and held:
            # A held value names atoms of the structure it was set on, and
            # "all" walks a set of different molecules. Applied to the rest it
            # is not a constraint but a distortion: a C-C held at 1.700 A on a
            # cyclobutane was written into benzene's xtb input as
            # "distance: 1, 2, 1.700000" with force constant 20, pulling an
            # aromatic bond a third of an Angstrom out of the ring.
            state['held_set_aside'] = len(held)
            held = []
        else:
            state.pop('held_set_aside', None)
        wet = str(submit_gfn_solvent.value or '') or None
        model = _solv_model()

        def _work():
            from .molecule_forcefield import relax_xyz
            # One xtb at a time.  An interrupted run is still shutting its
            # process down when the replacement starts, and a login node is
            # shared -- so the new one waits for the old one to be gone rather
            # than briefly doubling up.
            if earlier is not None and earlier.is_alive():
                earlier.join(timeout=15)
            # Two lists, because they mean opposite things to the rounds below.
            # A failure is a run that did not produce a geometry, and it ends
            # the press. "Stopped before converging" IS a geometry, and it is
            # the whole case the rounds exist for -- counted as a failure, as
            # it was at first, it stopped the rounds before they ever ran.
            results, failures, unfinished = [], [], []
            # The path of the frame on screen, kept where _apply can reach it
            # however the loop ended: a run stopped before its first outcome
            # leaves that name unbound, and reading it there is a crash in the
            # one place that exists to handle a run being cut short.
            trail = [None]
            targets = frames or [(single, None, None)]
            for position, item in enumerate(targets):
                xyz = item[0]
                try:
                    if _stopped():
                        failures.append(f'frame {position + 1}: stopped')
                        results.append(item)
                        continue
                    # The bonding the editor has been working with, so that
                    # pressing Optimise on a structure that has been pulled
                    # about does not re-perceive it and shove the molecule
                    # apart -- the same cliff the drag had.  Only for the one
                    # frame on screen: a set of isomers is a set of different
                    # molecules, and one perception cannot describe them all.
                    perceived = (None if every_frame
                                 else _gfn_topology_dir(xyz))
                    if pm:
                        # MOPAC takes the spin state as a word and knows
                        # MOPAC takes the spin state as a word and knows
                        # nothing of xtb's topology files, so those are not
                        # quietly passed along as though they had been
                        # honoured.  The held values it does take, in its own
                        # terms: its Cartesian input carries an optimisation
                        # flag per coordinate, and an atom a held value names
                        # is handed over with those at zero -- measured on a
                        # propane, a frozen carbon moves 0.0000 A where a free
                        # one moves 0.302.  That fixes the atoms rather than
                        # the value between them, and a pull cannot be said at
                        # all, so the run says which of the two it did.  The
                        # solvent it takes as well: its COSMO is handed the
                        # dielectric constant of the same liquid the GFN side
                        # is given by name.
                        outcome = _mopac.optimize_with_mopac(
                            xyz, method, charge=charge, uhf=uhf,
                            should_stop=_stopped, timeout=None, solvent=wet,
                            constraints=held,
                            on_frames=_push_frames if position == 0 else None)
                    elif gfn and autospin:
                        outcome = _gfn.optimize_autospin(
                            xyz, method, charge=charge, should_stop=_stopped,
                            timeout=None, on_frames=_push_frames,
                            constraints=held, topology=perceived, solvent=wet,
                            solvation_model=model)
                    elif gfn:
                        outcome = _gfn.optimize_with_gfn(
                            xyz, method, charge=charge, uhf=uhf,
                            should_stop=_stopped, timeout=None,
                            constraints=held, topology=perceived, solvent=wet,
                            solvation_model=model,
                            on_frames=_push_frames if position == 0 else None)
                    else:
                        # The held values, in the only terms a force field
                        # here has.  RDKit's UFF takes AddFixedPoint and Open
                        # Babel a list of fixed atoms; neither can restrain a
                        # value.  So a fix is met by holding the atoms that
                        # name it still, and a pull cannot be said at all --
                        # the same reading MOPAC's flags get, from the same
                        # function, so that a value one engine drops is
                        # dropped by the other for the same reason.
                        #
                        # Handed nothing at all, as this branch was, the whole
                        # list went quietly: measured on ethane with the
                        # bonding pinned, C0-H2 pulled out to 1.700 A and held
                        # exactly came back from one press of Optimise at
                        # 1.1104 A, under a line that said only "Optimised
                        # with UFF."  Fixed, it comes back at 1.7000, and a
                        # Zn-N held at 2.600 through Open Babel likewise.
                        atoms = len(_gfn.atom_lines(xyz))
                        frozen = _mopac.freeze_flags(held,
                                                     atoms=atoms or None)
                        # Whether anything is left free.  Reaching for an
                        # angle on a water names all three atoms, and the
                        # minimisation then has nothing to do -- which is
                        # worth a sentence rather than an "Optimised" over a
                        # geometry that did not move.
                        frozen['every_atom'] = bool(
                            atoms and len(frozen['frozen']) >= atoms)
                        outcome = relax_xyz(
                            xyz,
                            fixed_indices=sorted(frozen['frozen']),
                            max_steps=500,
                            perceived=_perception_for(xyz),
                            method=method,
                        )
                        # Read by the status line below in this engine's terms,
                        # through the same key the other two report on.  It is
                        # cleared at the start of every press, so a reading
                        # made under one engine cannot be read under another.
                        outcome['held'] = frozen
                except Exception as exc:
                    failures.append(f'frame {position + 1}: {exc}')
                    results.append(item)
                    continue
                if outcome.get('ok'):
                    if outcome.get('energy') is not None and position == 0:
                        state['gfn_energy'] = float(outcome['energy'])
                        state['gfn_energy_unit'] = outcome.get('energy_unit')
                    if position == 0:
                        # The charges of the structure that is about to be on
                        # screen, out of the answer that is about to draw it.
                        # Only the frame being looked at: a batch of isomers
                        # optimises one after another, and the charges of the
                        # fourth of them belong to a structure nobody is
                        # showing.
                        _remember_charges(outcome)
                        # Whether this run finished the job, and how far it got
                        # if it did not.  An engine that does not report either
                        # counts as finished: a run that cannot say it is
                        # unconverged must not be run again for ever.
                        state['optimize_converged'] = bool(
                            outcome.get('converged', True))
                        try:
                            state['optimize_moved'] = _gfn.largest_shift(
                                xyz, outcome['xyz'])
                        except Exception:
                            state['optimize_moved'] = 0.0
                        state['gfn_held'] = outcome.get('held')
                        if outcome.get('multiplicity'):
                            # What the scan settled on, so the live relaxation
                            # asks the same question this one did rather than
                            # a different one with the box's fixed M.
                            state['gfn_scanned_uhf'] = int(outcome['uhf'])
                    kept = outcome['xyz']
                    if _stopped() and outcome.get('frames'):
                        # Stop means the frame that was on screen.  xtb runs
                        # ahead of the picture, and geometries nobody saw are
                        # not what the user stopped at.
                        kept = _frame_as_xyz(
                            xyz, outcome['frames'],
                            _the_shown_frame_of(run_id),
                            'stopped at the frame on screen') or kept
                    results.append((kept,) + tuple(item[1:]))
                    if gfn and outcome.get('frames') and position == 0:
                        trail[0] = outcome['frames']
                        played[0] = True
                        # xtb writes every cycle to xtbopt.log, so the path
                        # costs nothing extra -- one run, and the viewer plays
                        # what the optimiser really walked through.
                        _push_frames(outcome['frames'], final=True)
                    note = str(outcome.get('status') or '')
                    if 'before converging' in note:
                        # It came back with a geometry, but not a finished one.
                        unfinished.append(f'frame {position + 1}: {note}')
                else:
                    failures.append(
                        f"frame {position + 1}: {outcome.get('status') or 'failed'}"
                    )
                    results.append(item)

            def _apply():
                if state.get('optimize_interrupted') is token:
                    # The structure changed under this run.  What it reached is
                    # a minimum of a geometry that no longer exists, so the
                    # switch is not touched: the run that replaces it is
                    # already on its way.
                    #
                    # The coordinates are, and to the frame the picture stood
                    # on when the hand arrived -- not to where xtb had got to,
                    # which nobody saw.  Left alone, the box kept the geometry
                    # from before the run while the browser showed the frame
                    # the user had taken hold of, and the two only came back
                    # together if a drag happened to push the browser's model
                    # over.  Taking hold and letting go without moving left
                    # them apart.
                    text = _frame_as_xyz(single or '', trail[0],
                                         _the_shown_frame_of(run_id),
                                         'stopped where you took hold')
                    if text is not None:
                        _write_coords(text, drawn=True)
                    return
                if _stopped():
                    # A Stop keeps the frame the picture stopped on, and it
                    # keeps it in the box as well as on the screen.  The
                    # ordinary write at the end of a run is refused here --
                    # the run number moved on at the Stop, which is what
                    # refuses the frames the run still had in hand -- and
                    # refused with them was the one geometry that is not
                    # stale: the one the user is looking at.  Measured, a Stop
                    # left the picture standing at frame 14 and the box
                    # holding the structure the run had started from, so Copy,
                    # Submit, Undo and the next press all read a geometry that
                    # was not on the screen.
                    #
                    # Put down rather than written, because which frame it is
                    # comes from the page: see _the_picture_stopped_here.
                    _the_picture_stopped_here(run_id, single or '', trail[0],
                                              'stopped at the frame on screen')
                # A run is one xtb --opt, and xtb's optimiser has a cycle limit
                # of its own: when it reaches it, it hands back the geometry it
                # got to and reports that it did not converge.  That was taken
                # as the end -- the switch went up over a structure that was
                # better than it had been and not at a minimum, and the user
                # pressed again.  Pressing again is a NEW process from that
                # geometry, with a fresh cycle budget and a fresh optimiser
                # history, which is the whole reason it got further.
                #
                # So the pressing is done here.  It ends the three ways Settle
                # has always ended -- converged, no longer moving, or out of
                # rounds -- because a structure whose held values cannot all be
                # met at once never converges, and a run that will not end is a
                # process per round for as long as the switch is down.
                rounds = int(state.get('optimize_rounds', 0)) + 1
                moved_now = float(state.get('optimize_moved') or 0.0)
                # Rounds in a row that moved next to nothing, this one counted.
                still = (int(state.get('optimize_still', 0)) + 1
                         if moved_now <= _OPTIMISE_STILL else 0)
                carry_on = _optimise_carries_on(
                    converged=bool(state.get('optimize_converged', True)),
                    moved=moved_now,
                    rounds=rounds,
                    failed=bool(failures),
                    every_frame=every_frame,
                    still=still,
                )
                state['optimize_rounds'] = rounds if carry_on else 0
                state['optimize_still'] = still if carry_on else 0
                # Converged, failed or stopped -- the switch goes back up by
                # itself, so it never claims to be working when it is not.
                if state.get('optimize_run') is token:
                    state['optimize_run'] = None
                if not carry_on:
                    for switch in (submit_optimize_btn, submit_optimize_all_btn):
                        if switch.value:
                            switch.value = False
                        switch.disabled = False
                state['pre_optimize_frames'] = {
                    'isomers': frames,
                    'coords': coords_widget.value,
                }
                if frames:
                    # Back the way they came, whatever the picture is doing:
                    # a tab that keeps its structures as blocks gets them as
                    # blocks, and one that steps through isomers gets the one
                    # it was on. Only the showing is held back while a
                    # trajectory plays -- showing an isomer rebuilds the
                    # viewer, and that is what used to tear the playback down.
                    #
                    # Held back was the handing over as well, and that lost
                    # every optimised geometry in a tab whose structures are
                    # blocks: the status line said "Optimised 2 of 2 frames"
                    # over a coordinates box that had not changed a character.
                    _offer_isomers(results, show=not played[0])
                elif failures:
                    # Nothing came back.  A run that produced no geometry hands
                    # the input straight back, and writing that to the box
                    # under "Optimised in DELFIN viewer" labelled the geometry
                    # the user had made as the answer to a question that was
                    # never answered -- an unoptimised structure wearing the
                    # word for an optimised one, with only the failure line
                    # underneath to say otherwise.  The box is left alone.
                    pass
                else:
                    lines = [
                        line for line in results[0][0].splitlines()[2:] if line.strip()
                    ]
                    if played[0]:
                        # The trajectory is playing, and its last frame is this
                        # very geometry.  Writing the coordinates the ordinary
                        # way rebuilds the viewer, which tore the playback down
                        # milliseconds after it started -- so only the end of
                        # the optimisation was ever seen.  The box is updated
                        # for Copy and Submit; the picture is already right.
                        pass
                    _write_coords(
                        xyz_document(lines, 'Optimised in DELFIN viewer'),
                        drawn=played[0], run=run_id)
                done = count - len(failures)
                # "1 of 1 frame(s)" is a count of a thing there is one of, and
                # it cost the line the width that pushed it onto a second row.
                said = (f'{label} could not optimise it.' if not done else
                        f'Optimised with {label}.' if count == 1 else
                        f'Optimised {done} of {count} frame(s) with {label}.')
                # The energy, the way the force field shows one.  xtb reports
                # it in hartree; kcal/mol is what the rest of the tab speaks,
                # and both are given because a total energy is compared
                # against other totals and a difference against chemistry.
                energy = state.get('gfn_energy')
                if energy is not None:
                    # And in the unit the engine that produced it uses. MOPAC
                    # reports a heat of formation in kcal/mol, not a total
                    # energy in hartree: read as hartree and converted, its
                    # 15.35 came out as 9629.68 kcal/mol, which is not a
                    # number about this molecule at all.
                    unit = state.get('gfn_energy_unit')
                    if unit:
                        said += f' E = {energy:.4f} {unit}.'
                    else:
                        said += (f' E = {energy:.6f} Eh '
                                 f'({energy * 627.5094740631:.2f} kcal/mol).')
                # What was held while it ran, and what became of it.  A value
                # the user is holding on screen that the optimisation quietly
                # ignored would make the result an answer to a question nobody
                # asked.
                said += _solvents.note(_solv_model(), submit_gfn_solvent.value)
                # In the terms of whichever engine ran.  xtb restrains the
                # value with one force constant for the whole set; MOPAC and
                # the force fields here fix the atoms that name it and cannot
                # express a pull at all.  Read out with the wrong one, a MOPAC
                # result would claim a force constant that no MOPAC run has.
                kept = state.get('gfn_held')
                if pm:
                    said += _mopac.freeze_note(kept or {
                        'held': 0, 'pulls': 0, 'dropped': [], 'frozen': set()})
                elif gfn:
                    said += _gfn.held_note(kept or {
                        'held': 0, 'dropped': [], 'mixed': False,
                        'force': None})
                else:
                    said += fixed_atoms_note(kept or {
                        'held': 0, 'pulls': 0, 'dropped': [], 'frozen': set()})
                aside = int(state.pop('held_set_aside', 0) or 0)
                if aside:
                    # Said out loud rather than done quietly: a value held on
                    # one structure names atoms of that structure, and "all"
                    # walks a set of different molecules.
                    said += (f' The {aside} held value(s) were not applied -- '
                             'they name atoms of the structure they were set '
                             'on. Optimise that one on its own to keep them.')
                state['gfn_last_status'] = said
                if gfn:
                    if autospin:
                        picked = results[0][0] if results else None
                        del picked
                        said += f' charge {charge}, multiplicity scanned.'
                    else:
                        said += f' charge {charge}, multiplicity {uhf + 1}.'
                # The page's count of the frames it drew does not go under
                # this. It arrives while this message is being built, so the
                # two took turns -- one row, then two, then one -- and the
                # viewer stands under them and moved with every swap. What the
                # playback did wrong still gets said; that it worked does not
                # need a row of its own.
                state['gfn_last_status'] = said
                if carry_on:
                    # Where it has got to, rather than a result: the same press
                    # is still running and the coordinates above are what the
                    # next round starts from.
                    said = (f'{label} is still going: round {rounds}, '
                            f'{state.get("optimize_moved", 0.0):.3f} A moved.')
                    state['gfn_last_status'] = said
                    _set_mol_status(said, spinner=True)
                    state['optimize_carrying_on'] = True
                    schedule_ui_update(
                        lambda: on_submit_optimize(None, every_frame=every_frame))
                    return
                _set_mol_status(said, *(failures + unfinished)[:2])

            schedule_ui_update(_apply)

        _start_background(_work, 'The optimisation',
                          guards={'optimize_run': None},
                          remember_in='optimize_thread')

    def _clear_selection():
        """Drop the picks so the next constraint starts from a clean set."""
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.clearSelection('
            f'{json.dumps(submit_scope_id)});'
        )

    def _step_for_selection(indices):
        """How far one press of an arrow key moves the value.

        Three different quantities, three different steps.  A bond length is
        Angstrom, where a hundredth is fine.  An angle is degrees and a tenth
        is the useful step.  A dihedral is what gets turned through a whole
        rotation -- half a degree there, so holding the key sweeps it instead
        of creeping.
        """
        picked = len(indices or ())
        if picked == 2:
            submit_internal_value.step = 0.01     # bond length, in Angstrom
        elif picked == 4:
            submit_internal_value.step = 0.5      # dihedral, swept by hand
        else:
            submit_internal_value.step = 0.1      # angle

    def _apply_internal_now():
        """Put the selection at the value in the box, and leave it selected."""
        # The browser is about to move the atoms, so this is the last moment
        # the structure before the move exists anywhere. Set went into no
        # history at all: turning a dihedral by hand and pressing Undo took
        # back whatever had been done before the turn and left the turn
        # standing.
        #
        # One entry for the sweep, not one per press: with Set on, an arrow
        # key comes through here every time it repeats.
        picked = tuple(state.get('picked') or ())
        if len(picked) in _CONSTRAINT_KINDS:
            # Only where there is something to set. With nothing selected the
            # browser moves no atom, and an entry for it is a press of Undo
            # that appears to do nothing.
            _remember(f'setting {_describe_selection(picked)}',
                      gesture=('set', picked))
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setInternal('
            + json.dumps(submit_scope_id) + ','
            + repr(float(submit_internal_value.value)) + ');'
        )
        # The browser moves the atoms; under GFN the rest of the molecule has
        # to be given the chance to arrange itself around where they landed.
        _arm_gfn_takeup('Set')

    def on_submit_set_internal(change=None):
        """Set is a mode: while it is on, the box turns the selection by hand.

        Switching it on puts the selection at what the box says, and every
        further change of the box does the same -- so the arrow keys turn a
        dihedral a tenth of a degree at a time and the structure follows.  The
        picks are kept, because letting go of them after every step is what
        made turning something by hand impossible.
        """
        if not submit_internal_btn.value:
            return
        _apply_internal_now()

    def on_submit_internal_value(change):
        """The box changed.  Who owns it depends on what is selected."""
        if change.get('name') != 'value' or state.get('internal_quiet'):
            return
        if _selected_constraint()[1] is not None:
            return          # a held value is being retuned, not the geometry
        if submit_internal_btn.value:
            _apply_internal_now()

    def on_submit_pick_sync(change):
        """Offer the coordination polyhedra of a metal the moment it is picked.

        Which ones are possible follows from its coordination number, and the
        tables are MANTA's own -- the same ideal donor vectors it builds
        complexes with.
        """
        if change.get('name') != 'value':
            return
        raw = (submit_pick_sync.value or '').strip()
        record('pick', v=raw)
        indices = [int(part) for part in raw.split(',') if part.strip().isdigit()]
        state['picked'] = indices
        _step_for_selection(indices)
        _refresh_swap(indices)
        # The two scan fields belong to the selection: they appear when there
        # is a value to walk and go away again when there is not, so the row
        # does not carry two empty boxes about for the whole session.
        _refresh_scan()
        xyz = _current_xyz()
        options = None
        perceived = None
        if xyz and indices:
            # Perceive on demand. Waiting for the cache meant the offer only
            # appeared once the force field had been switched on at least
            # once -- tapping a metal before that did nothing at all, and said
            # nothing either.
            try:
                perceived = _perception_for(xyz)
            except Exception:
                perceived = None
        # A polyhedron is a set of restraints in the browser's own field, so
        # it is only offered where that field is the one running.  Under a
        # server method it was offered, accepted, and reported as "the donors
        # are pulled onto it" while the same press took the field off the page.
        on_the_server = _server_method()
        if (perceived is not None and len(indices) == 1
                and not on_the_server):
            try:
                from .molecule_forcefield import polyhedron_options
                options = polyhedron_options(perceived, indices[0])
            except Exception:
                options = None
        _refresh_hybridisation(indices, perceived)
        if not options:
            # Only the offer follows the selection. The applied polyhedron
            # stays: clearing it here meant that selecting three ligand atoms
            # to hold an angle silently threw the polyhedron away, and the very
            # next export went out without it.
            submit_poly_dd.layout.display = 'none'
            submit_poly_dd.disabled = True
            state['poly_offer_metal'] = None
            # Say why, when a single atom was picked and could have qualified.
            # Never under a server method: there the offer is absent because
            # the method has no use for it, and blaming the coordination
            # number would send the user looking for a table that is not the
            # reason.
            if (perceived is not None and len(indices) == 1
                    and not on_the_server):
                index = indices[0]
                symbol = perceived.symbols[index]
                if index in set(perceived.metal_indices or ()):
                    donors = sorted(
                        j for pair in perceived.bonds for j in pair
                        if index in pair and j != index
                    )
                    _set_mol_status(
                        f'{symbol}: coordination number {len(donors)} — no '
                        'polyhedron table for that (2 to 9 are covered).'
                    )
            return

        coordination, current, choices = options
        state['poly_offer_metal'] = indices[0]
        symbol = perceived.symbols[indices[0]]
        entries = [(f'{symbol} · CN {coordination} · free', '')]
        for code, label in choices:
            mark = ' (current)' if code == current else ''
            entries.append((f'{label}{mark}', code))
        state['poly_quiet'] = True
        try:
            submit_poly_dd.options = entries
            # Only a code this metal actually offers. A polyhedron held on one
            # metal has nothing to say about the next one picked, and assigning
            # a value the options do not contain raises.
            applied = state.get('poly_applied')
            offered = {code for code, _label in choices}
            submit_poly_dd.value = applied if applied in offered else ''
        finally:
            state['poly_quiet'] = False
        submit_poly_dd.disabled = False
        submit_poly_dd.layout.display = ''

    def _refresh_hybridisation(indices, perceived):
        """Offer to overrule the hybridisation of the picked atoms.

        Any number of them: a double bond that went unperceived usually cost
        both of its carbons their type, and retyping a ring one atom at a time
        is busywork. Metals are dropped from the selection rather than
        blocking it -- RDKit's UFF has no types for one at all, so its bonds
        and angles come from the geometry either way.

        Not offered at all under a server method: a type is what the browser's
        field builds its parameters from, and xtb and MOPAC have no such thing
        to be told. See :func:`_refresh_method_controls` for what was measured.
        """
        metals = set(perceived.metal_indices or ()) if perceived else set()
        # An index the structure no longer has: the browser pushes its picks
        # after a re-render, and an edit that deleted atoms renumbers them.
        # Asking the perception about one is an IndexError, and this handler
        # runs on every click in the viewer.
        chosen = [
            i for i in indices
            if i not in metals and 0 <= i < len(perceived.symbols)
        ] if perceived and not _server_method() else []
        if not chosen:
            submit_hyb_dd.layout.display = 'none'
            submit_hyb_dd.disabled = True
            state['hyb_offer_atoms'] = []
            return

        from .molecule_forcefield import (
            HYBRIDISATION_CHOICES, perceived_hybridisation_of,
        )

        overrides = state.get('hyb_overrides') or {}
        auto = {perceived_hybridisation_of(perceived, i) for i in chosen}
        if len(chosen) == 1:
            index = chosen[0]
            head = (f'{perceived.symbols[index]}{index} · '
                    f'{auto.pop() or "no type"} (automatic)')
        else:
            named = ''.join(sorted({perceived.symbols[i] for i in chosen}))
            found = auto.pop() if len(auto) == 1 else 'mixed'
            head = f'{len(chosen)} atoms ({named}) · {found or "no type"} (automatic)'
        entries = [(head, '')]
        for name in HYBRIDISATION_CHOICES:
            entries.append((f'force {name}', name))
        # Only a value they all already share can be shown as the current one.
        held = {overrides.get(i, '') for i in chosen}
        state['hyb_offer_atoms'] = chosen
        state['hyb_quiet'] = True
        try:
            submit_hyb_dd.options = entries
            submit_hyb_dd.value = held.pop() if len(held) == 1 else ''
        finally:
            state['hyb_quiet'] = False
        submit_hyb_dd.disabled = False
        submit_hyb_dd.layout.display = ''

    def on_submit_hyb_changed(change):
        if change.get('name') != 'value' or state.get('hyb_quiet'):
            return
        atoms = list(state.get('hyb_offer_atoms') or [])
        if not atoms:
            return
        index = atoms[0]
        chosen = submit_hyb_dd.value or ''
        # A type is what the field builds its angles from, so typing an atom
        # sp2 flattens the centre around it as surely as pulling on it would.
        _remember(f'typing {len(atoms)} atom(s) {chosen}' if chosen
                  else f'letting {len(atoms)} atom(s) back to the perceived type')
        overrides = dict(state.get('hyb_overrides') or {})
        for atom in atoms:
            if chosen:
                overrides[atom] = chosen
            else:
                overrides.pop(atom, None)
        state['hyb_overrides'] = overrides
        # Perception is cached by element sequence, which this does not
        # change, so the cache has to be dropped or the override never reaches
        # the force field.
        state['perceived'] = None
        state['perceived_for'] = None
        state['poly_assignment'] = None
        _enable_live_forcefield()
        perceived = state.get('perceived')
        if len(atoms) == 1:
            symbol = perceived.symbols[index] if perceived else '?'
            named = f'{symbol}{index}'
        else:
            named = f'{len(atoms)} atoms'
        if chosen:
            shape = {'sp': 'linear', 'sp2': 'trigonal planar',
                     'sp3': 'tetrahedral'}[chosen]
            _set_mol_status(f'{named} typed as {chosen}: {shape}.')
        else:
            _set_mol_status(f'{named} back to the perceived type.')
        _clear_selection()

    def _refresh_poly_turn():
        """Offer Turn only where the vertices are not all alike.

        An octahedron has nothing to turn -- every vertex is the same, and
        which ligand is trans to which is what Swap is for. A trigonal
        bipyramid has two kinds, so which pair is axial is a real choice.

        And never under a server method, where an arrangement is a rearranged
        set of restraints for a field that is not running.
        """
        geometry = state.get('poly_applied')
        metal = state.get('poly_metal')
        perceived = state.get('perceived')
        turnable = False
        if (geometry and metal is not None and perceived is not None
                and not _server_method()):
            try:
                from .molecule_forcefield import polyhedron_vertex_classes
                donors = len(perceived.neighbours()[int(metal)])
                grouped = polyhedron_vertex_classes(donors, geometry)
                turnable = bool(grouped) and len(set(grouped[0])) > 1
            except Exception:
                turnable = False
        submit_poly_turn_btn.layout.display = '' if turnable else 'none'
        submit_poly_turn_btn.disabled = not turnable
        # "Nothing to turn" and "not perceived yet" are different answers, and
        # only the first of them may throw the arrangements away.
        #
        # Undo and Redo both drop the perception before they call this -- the
        # structure has changed, so what was perceived is about another one --
        # and this then cleared the very arrangement they had just put back.
        # Measured on a trigonal-bipyramidal iron: Turn to arrangement 2, take
        # it back, put it forward again, and the index came back 0 with the
        # list empty, so the next press started from the beginning instead of
        # stepping on.  Undo had it too and it never showed, because what it
        # was restoring was empty anyway.
        undecided = (geometry and metal is not None and perceived is None)
        if not turnable and not undecided:
            state['poly_arrangements'] = []
            state['poly_arrangement_index'] = 0

    def on_submit_poly_turn(_button=None):
        """Step to the next way the ligands can sit on this polyhedron."""
        geometry = state.get('poly_applied')
        metal = state.get('poly_metal')
        xyz = _current_xyz()
        if not geometry or metal is None or not xyz:
            _set_mol_status('Choose a polyhedron for a metal first.')
            return
        try:
            from .molecule_forcefield import (
                describe_polyhedron_arrangement, parse_xyz,
                polyhedron_arrangements,
            )
            perceived = _perception_for(xyz)
            # The coordinates as they are now, not as they were perceived: a
            # ligand that has been dragged has to be scored where it sits.
            parsed = parse_xyz(xyz)
            coords = perceived.coords
            if parsed is not None and list(parsed[0]) == list(perceived.symbols):
                coords = parsed[1]
            arrangements = polyhedron_arrangements(
                perceived, int(metal), geometry, coords)
        except Exception as exc:
            _set_mol_status(f'Could not work out the arrangements: {exc}')
            return
        if len(arrangements) < 2:
            _set_mol_status(
                'Every vertex of this polyhedron is the same, so there is '
                'nothing to turn. Swap exchanges two ligands.'
            )
            return
        position = (int(state.get('poly_arrangement_index') or 0) + 1) % len(arrangements)
        # Before the arrangement changes. Turn moves no coordinate itself --
        # it says which vertex each ligand belongs on, and the field walks
        # them there over the next second -- so there was nothing for Undo to
        # find, and pressing it took back the step before instead while the
        # ligands went on moving into the arrangement it was meant to undo.
        _remember(f'turning to arrangement {position + 1}')
        state['poly_arrangements'] = arrangements
        state['poly_arrangement_index'] = position
        state['poly_assignment'] = arrangements[position]
        _enable_live_forcefield()
        described = describe_polyhedron_arrangement(
            perceived, geometry, arrangements[position])
        _set_mol_status(
            f'Arrangement {position + 1} of {len(arrangements)} — {described}.'
        )

    def on_submit_hyb_auto(_button=None):
        """Derive the carbon types from the connectivity and hold them."""
        xyz = _current_xyz()
        if not xyz:
            _set_mol_status('Load a structure first.')
            return
        try:
            from .molecule_forcefield import hybridisation_from_connectivity
            perceived = _perception_for(xyz)
            picked = list(state.get('picked') or [])
            derived = hybridisation_from_connectivity(perceived, picked or None)
        except Exception as exc:
            _set_mol_status(f'Could not read the types: {exc}')
            return
        if not derived:
            _set_mol_status(
                'No carbon in the selection — the count only fixes the shape '
                'for carbon, which has no lone pair.'
            )
            return
        changed = [
            i for i, name in derived.items()
            if (state.get('hyb_overrides') or {}).get(i) != name
        ]
        _remember(f'typing {len(derived)} carbon(s) from their partners')
        overrides = dict(state.get('hyb_overrides') or {})
        overrides.update(derived)
        state['hyb_overrides'] = overrides
        state['perceived'] = None
        state['perceived_for'] = None
        state['poly_assignment'] = None
        _enable_live_forcefield()
        where = f'{len(picked)} selected' if picked else 'the whole structure'
        counts = ', '.join(
            f'{sum(1 for v in derived.values() if v == name)}x {name}'
            for name in ('sp', 'sp2', 'sp3')
            if any(v == name for v in derived.values())
        )
        _set_mol_status(
            f'{len(derived)} carbons typed from their partners in {where} '
            f'({counts}); {len(changed)} changed.'
        )
        _clear_selection()

    #: How many structural edits can be taken back.  The browser keeps 50
    #: coordinate snapshots; there is no reason for this to be shorter.
    _STRUCTURE_UNDO_LIMIT = 50
    #: How many steps back a session keeps. Generous, because each is a
    #: coordinate block and a hundred of a 400-atom structure is a few
    #: megabytes -- against a user who cannot get back to where they started.
    _HISTORY_LIMIT = 200

    _CONSTRAINT_KINDS = {2: 'distance', 3: 'angle', 4: 'dihedral'}

    def _describe_selection(indices):
        """What a set of picked atoms describes, for a line about it."""
        kind = _CONSTRAINT_KINDS.get(len(indices or ()))
        return f'the {kind}' if kind else 'the selection'

    def _describe_constraint(entry):
        symbols = []
        perceived = state.get('perceived')
        for index in entry['atoms']:
            symbol = perceived.symbols[index] if perceived else '?'
            symbols.append(f'{symbol}{index}')
        unit = 'A' if entry['kind'] == 'distance' else 'deg'
        mode = entry.get('mode', 'pull')
        return f"{'-'.join(symbols)} = {entry['value']:.3g} {unit} ({mode})"

    def _carry_constraints_to(method):
        """Hand the held values over to the engine that has just been chosen.

        The list itself survives a change of method -- it describes the
        molecule, not the program -- but what an engine will do with it differs,
        and the three here differ completely:

        xtb holds internal coordinates, one force constant for the whole set.
        The browser's field holds the pulls and cannot fix anything, so a fix
        is met by pulling very hard rather than by the atom not moving.
        MOPAC is handed none of them: it takes no constraint block from this
        editor, so a value held on screen while PM7 runs is a value the run
        never hears about.

        That last one is why this says something instead of only re-installing
        the parameters.  A held bond that stays in the list and stops being
        held is the list describing a thing that is not happening, and the run
        that follows answers a question the user did not ask.
        """
        held = list(state.get('constraints') or [])
        # The browser's field carries its restraints in the parameters it was
        # given, so a new method means handing them over again.  A server
        # method reads the list when it runs and needs nothing installed.
        if not _server_method(method) and submit_relax_btn.value:
            _enable_live_forcefield()
        if not held:
            return ''
        if _mopac.is_mopac_method(method):
            reading = _mopac.freeze_flags(
                held, atoms=len(_gfn.atom_lines(_current_xyz() or '')) or None)
            return (f'{_server_label(method)} holds them its own way.'
                    + _mopac.freeze_note(reading))
        if _gfn.is_gfn_method(method):
            # Said from the same function the run uses, so what is promised
            # here is what will actually be written into xtb.inp.
            reading = _gfn.constraint_input(
                held, atoms=len(_gfn.atom_lines(_current_xyz() or '')) or None)
            if reading['dropped']:
                return (
                    f'{_server_label(method)} will hold '
                    f'{reading["held"]} of the {len(held)} value(s); '
                    f'{len(reading["dropped"])} name atoms this structure does '
                    'not have, or are not a distance, an angle or a dihedral.')
            if reading['mixed']:
                return (
                    f'{_server_label(method)} will hold all {len(held)} '
                    'value(s). xtb takes one force constant for the whole set, '
                    'so the pulls are held as firmly as the exact values.')
            return (f'{_server_label(method)} will hold all {len(held)} '
                    'held value(s).')
        fixed = sum(1 for entry in held if entry.get('mode') == 'fix')
        if fixed:
            return (
                f'{str(method).upper()} runs in the browser, which pulls '
                f'towards a held value rather than fixing it: the {fixed} '
                'exact one(s) are held very firmly, not held still. One press '
                'of Optimise meets them exactly, by fixing the atoms they '
                'name.')
        return ''

    def _selected_constraint():
        """The held entry the list is pointing at, or (None, None)."""
        key = submit_constraint_dd.value or ''
        if not key.startswith('c'):
            return None, None
        held = list(state.get('constraints') or [])
        position = int(key[1:])
        if not (0 <= position < len(held)):
            return None, None
        return position, held[position]

    def _sync_constraint_selection():
        """Mark the atoms of the selected entry and offer its value for editing.

        With nothing selected the picture keeps whatever the user picked and the
        value box goes back to serving the selection, which is what it does when
        no held value is being looked at.
        """
        _position, entry = _selected_constraint()
        if entry is None:
            return
        state['hold_mode_quiet'] = True
        try:
            submit_hold_mode.value = entry.get('mode', 'pull')
        finally:
            state['hold_mode_quiet'] = False
        kind = str(entry['kind'])
        unit = 'Å' if kind == 'distance' else '°'
        number = float(entry['value'])
        value = '{:.3f}'.format(number) if kind == 'distance' else '{:.1f}'.format(number)
        label = 'held <b>{}</b> ({})'.format(kind, unit)
        atoms = [int(i) for i in entry['atoms']]
        _run_manip_js(
            'if(window.__delfinSubmitManip&&window.__delfinSubmitManip.setPicks)'
            'window.__delfinSubmitManip.setPicks('
            + json.dumps(submit_scope_id) + ','
            + json.dumps(atoms) + ','
            + json.dumps(value) + ','
            + json.dumps(label) + ');'
        )

    def on_submit_constraint_selected(change):
        if change.get('name') != 'value' or state.get('constraint_quiet'):
            return
        if not (submit_constraint_dd.value or ''):
            _clear_selection()
            return
        _sync_constraint_selection()

    def on_submit_constraint_retune(change):
        """Editing the value while a held entry is selected retunes that entry."""
        if change.get('name') != 'value':
            return
        position, entry = _selected_constraint()
        if entry is None:
            return
        if abs(float(submit_internal_value.value) - float(entry['value'])) < 1e-9:
            return  # the box was filled with the held value, not edited
        _remember(f'retuning {_describe_constraint(entry)}')
        held = list(state.get('constraints') or [])
        held[position] = dict(entry, value=float(submit_internal_value.value))
        state['constraints'] = held
        _refresh_constraints()
        _enable_live_forcefield()
        # And under GFN, where nothing runs between drags, the change is what
        # has to start it.  This was the one of the three ways of altering a
        # held value that did not: Hold arms it, changing pull to fix arms it,
        # and typing a new number into the box did not -- so the list said one
        # thing, the structure went on standing at another, and pressing Hold
        # again was what made it happen.  The browser's field needs no such
        # push, which is why it only showed under a server method.
        _arm_gfn_takeup(f'Holding {_describe_constraint(held[position])}')
        _set_mol_status(f'Holding {_describe_constraint(held[position])}.')

    def _refresh_constraints():
        """Show what the field is currently being held to."""
        entries = []
        if state.get('poly_applied') and state.get('poly_metal') is not None:
            entries.append(('poly', f"polyhedron: {state['poly_applied']}"))
        for position, entry in enumerate(state.get('constraints') or []):
            entries.append((f'c{position}', _describe_constraint(entry)))
        visible = bool(entries)
        submit_constraint_dd.layout.display = '' if visible else 'none'
        submit_constraint_del.layout.display = '' if visible else 'none'
        submit_constraint_dd.disabled = not visible
        submit_constraint_del.disabled = not visible
        if not visible:
            # Leave nothing behind: a hidden dropdown that still holds the
            # last constraint shows it again the moment another one is set.
            state['constraint_quiet'] = True
            try:
                submit_constraint_dd.options = [('no constraints', '')]
                submit_constraint_dd.value = ''
            except Exception:
                pass
            finally:
                state['constraint_quiet'] = False
            return
        # Nothing is selected to begin with.  Selecting an entry marks the atoms
        # it holds, so a preselected one would mean the picture always shows a
        # marked set nobody asked for.
        previous = submit_constraint_dd.value
        state['constraint_quiet'] = True
        try:
            submit_constraint_dd.options = (
                [(f'{len(entries)} held · show which', '')]
                + [(label, key) for key, label in entries]
            )
            submit_constraint_dd.value = (
                previous if previous in dict(submit_constraint_dd.options).values() else ''
            )
        finally:
            state['constraint_quiet'] = False
        _sync_constraint_selection()

    def _refresh_swap(indices):
        """Offer an exchange whenever two donors of one metal are selected.

        It does not need a polyhedron: exchanging two ligands is a move across
        a barrier, which is useful with or without a target shape.
        """
        perceived = state.get('perceived')
        metals = set(getattr(perceived, 'metal_indices', ()) or ())
        donors = set()
        if perceived is not None and len(metals) == 1:
            metal = next(iter(metals))
            donors = {
                j for pair in perceived.bonds for j in pair
                if metal in pair and j != metal
            }
        ready = len(indices) == 2 and donors and all(i in donors for i in indices)
        submit_swap_btn.layout.display = '' if ready else 'none'
        submit_swap_btn.disabled = not ready

    def _edit_bond(connect):
        """Draw or remove a bond, and remember it.

        Bond perception reads distances, and in a crowded coordination sphere
        that is simply not reliable: on a real Pt complex it counted two ipso
        carbons of a phosphine's phenyls as donors, giving CN 6 for a
        four-coordinate metal, while the viewer's own perception invented a
        Pt-H bond instead. Neither is trustworthy, so the correction has to be
        remembered and re-applied -- otherwise the next perception, which runs
        from the geometry, would quietly undo it.
        """
        indices = list(state.get('picked') or [])
        if len(indices) != 2:
            _set_mol_status('Select exactly two atoms to change a bond.')
            return
        pair = (min(indices), max(indices))
        # Before anything is changed, and for every way through this: with
        # Adjust H on the structural path below records the step itself, and
        # with it off nothing did -- so drawing a bond, cutting one, or
        # correcting a perception went into no history at all and Undo walked
        # straight past to whatever had been done before it. Two presses of a
        # button that does not undo the thing that was just done are worse
        # than a button that says it has nothing to take back. The structural
        # path's own entry lands on this same state and is folded into it.
        _remember(f'{"drawing" if connect else "cutting"} the bond between '
                  f'{pair[0]} and {pair[1]}')
        # A bond drawn where there was none costs both ends a valence, and with
        # Adjust H on the hydrogen standing in its way goes: two methanes
        # bonded at the carbons are ethane, not C2H8. With it off nothing goes,
        # which is what correcting a perception needs -- the bond was always
        # there, so nothing has to make room for it. Removing an atom is a
        # structural edit and takes the structural path; the bond itself is
        # still only recorded and drawn, never relaxed into place.
        made_room = 0
        if connect and bool(submit_adjust_h_btn.value):
            from .molecule_builder import connect_atoms

            structure = _structure_now()
            if structure is not None:
                before = len(structure)
                ends = [structure.symbols[i] if i < len(structure) else '?'
                        for i in pair]
                mapping = connect_atoms(structure, pair[0], pair[1])
                if mapping:
                    made_room = before - len(structure)
                    _apply_structure(
                        structure,
                        f'Bonded {ends[0]}{pair[0]} and {ends[1]}{pair[1]}.')
                    pair = (mapping.get(pair[0], pair[0]),
                            mapping.get(pair[1], pair[1]))
        edits = {tuple(k): v for k, v in (state.get('bond_edits') or {}).items()}
        edits[pair] = bool(connect)
        state['bond_edits'] = edits
        # Separately from the bond list, because after a structural edit that
        # list is the whole structure and no longer says which bonds came from
        # a hand rather than from perception -- and only the ones from a hand
        # have to be put back into a picture that has been rebuilt.
        hand = {tuple(k): v for k, v in (state.get('hand_bonds') or {}).items()}
        hand[pair] = bool(connect)
        state['hand_bonds'] = hand
        # The perception is cached by element sequence, which a bond edit does
        # not change -- so the cache has to be dropped explicitly, or the
        # correction would never reach the force field at all.
        state['perceived'] = None
        state['perceived_for'] = None
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.editBond('
            f'{json.dumps(submit_scope_id)},{pair[0]},{pair[1]},'
            f'{"true" if connect else "false"});'
        )
        # Re-assign the parameters straight away: the topology decides the
        # bonds, angles and torsions the field works with, and until it is
        # rebuilt the relaxation is still holding the bond that was just cut.
        state['poly_assignment'] = None
        _enable_live_forcefield()
        verb = 'Bonded' if connect else 'Unbonded'
        room = (f' {made_room} hydrogen(s) made room for it.'
                if made_room else '')
        _set_mol_status(f'{verb} atoms {pair[0]} and {pair[1]}.{room}')
        _clear_selection()

    def _apply_structure(structure, note):
        """Put an edited structure back and let the tab rebuild around it.

        The coordinate box is the tab's single source of truth, so writing to
        it re-renders the viewer and re-perceives everything. What does *not*
        survive that is the bond orders: perception reads them off the
        geometry and does not get them back -- ethene built here came back as
        a single bond at 1.514 A. So the edited topology is re-seeded as the
        hand corrections it is, which is exactly what it is: the user built
        this, and their bonds outrank a distance table until a different
        structure is loaded.
        """
        from .molecule_builder import to_xyz

        _remember(note.rstrip('.') or 'an edit')

        # Before the coordinates are written, because writing them rebuilds the
        # picture and the rebuild is where the hand corrections are re-applied:
        # at that moment they have to name the atoms of the structure that is
        # arriving, not of the one being left behind. A removed hydrogen moves
        # every atom after it, so a bond re-applied on the old numbering would
        # be drawn between the wrong two atoms.
        moved = structure.renumbering()
        state['hand_bonds'] = {
            (min(moved[i], moved[j]), max(moved[i], moved[j])): connected
            for (i, j), connected in (state.get('hand_bonds') or {}).items()
            if i in moved and j in moved
        }

        # An atom added or taken away is a different molecule from the one xtb
        # has in hand, down to the number of coordinates -- so a run under it
        # is ended here and starts again on the molecule that now exists, and
        # the bonding GFN-FF had perceived goes with it.
        _drop_gfn_topology()
        _interrupt_gfn()
        _arm_gfn_restart()

        xyz = to_xyz(structure, note)
        lines = [line for line in xyz.splitlines()[2:] if line.strip()]
        # The write below re-renders through update_view, which
        # clears the history a new structure invalidates. This is not a new
        # structure, it is a step in the one being edited.
        state['structure_edit_inflight'] = True
        _mark_structure_edit()
        try:
            coords_widget.value = f'{len(lines)}\n{note}\n' + '\n'.join(lines)
        finally:
            state['structure_edit_inflight'] = False
        # After update_view, which clears them.
        state['bond_edits'] = {
            (int(i), int(j)): int(order)
            for (i, j), order in structure.bonds.items()
        }
        state['perceived'] = None
        state['perceived_for'] = None
        # Types, derived and held, after every build step -- which is what
        # pressing the button by hand after one was doing. Perception reads
        # bond orders off the geometry, and a structure that has just been
        # built is exactly where that geometry is least settled: a centre
        # comes back sp3 and its angles at 109.5 where they should be 120, so
        # the field pulls the new part into the wrong shape. The number of
        # partners says it outright and does not depend on the geometry at
        # all, which is why doing it by hand fixed things.
        # Everything the tab holds by index has to follow the renumbering an
        # edit causes -- a deleted hydrogen moves every atom after it. A held
        # value that quietly pointed at different atoms, or vanished, is worse
        # than one that is dropped and said so.
        renumber = structure.renumbering()

        def _follow(indices):
            moved = [renumber.get(int(i)) for i in indices]
            return None if any(x is None for x in moved) else moved

        kept, lost = [], 0
        for entry in (state.get('constraints') or []):
            moved = _follow(entry.get('atoms') or [])
            if moved is None:
                lost += 1
                continue
            kept.append(dict(entry, atoms=moved))
        state['constraints'] = kept

        state['hyb_overrides'] = {
            renumber[i]: name
            for i, name in (state.get('hyb_overrides') or {}).items()
            if i in renumber
        }
        metal = state.get('poly_metal')
        if metal is not None:
            if metal in renumber:
                state['poly_metal'] = renumber[metal]
                assignment = state.get('poly_assignment') or {}
                followed = {
                    renumber[d]: v for d, v in assignment.items() if d in renumber
                }
                state['poly_assignment'] = (
                    followed if len(followed) == len(assignment) else None)
            else:
                state['poly_applied'] = None
                state['poly_metal'] = None
                state['poly_assignment'] = None

        # Straight off the structure that was just built, not off a fresh
        # perception of it: perception can fail or come back empty at exactly
        # this moment, and then nothing was derived at all -- which is why the
        # button still had to be pressed by hand. The builder knows every bond
        # it made, so the count is certain.
        by_count = {2: 'sp', 3: 'sp2', 4: 'sp3'}
        derived = {}
        for index, symbol in enumerate(structure.symbols):
            if symbol != 'C':
                continue
            partners = structure.neighbours(index)
            # A side-on alkene is the one case where the metal does not count.
            side_on = [
                m for m in partners
                if structure.symbols[m] not in ('H', 'C', 'N', 'O', 'S', 'P')
                and any(o in structure.neighbours(m) for o in partners
                        if structure.symbols[o] == 'C')
            ]
            name = by_count.get(len(partners) - len(side_on))
            if name:
                derived[index] = name
        if derived:
            overrides = dict(state.get('hyb_overrides') or {})
            overrides.update(derived)
            state['hyb_overrides'] = overrides
            state['perceived'] = None
            state['perceived_for'] = None
        if lost:
            note = f'{note} {lost} held value(s) lost their atoms and were dropped.'
        _set_mol_status(note)
        _refresh_constraints()
        _push_bond_orders(structure.bonds)
        if submit_relax_btn.value or state.get('ff_bootstrap_done'):
            _enable_live_forcefield()

    def _structure_marks():
        """Everything about the structure the coordinates do not say.

        A held value, a bond drawn or cut by hand, an atom typed sp2, a
        polyhedron and the arrangement its ligands sit in: none of them is a
        coordinate, and every one of them decides what the field does with the
        coordinates next. Undo used to put back three of these and leave the
        rest standing, so taking back "type this carbon sp2" gave the geometry
        back and kept the typing -- and the field pulled the structure into the
        shape the typing asks for all over again.

        Whatever an action can change belongs in here, or Undo does not undo
        that action.
        """
        return {
            'bond_edits': dict(state.get('bond_edits') or {}),
            'hand_bonds': dict(state.get('hand_bonds') or {}),
            'constraints': [dict(one)
                            for one in (state.get('constraints') or [])],
            'hyb_overrides': dict(state.get('hyb_overrides') or {}),
            'poly_applied': state.get('poly_applied'),
            'poly_metal': state.get('poly_metal'),
            'poly_assignment': (dict(state['poly_assignment'])
                                if state.get('poly_assignment') else None),
            'poly_arrangements': [dict(one) for one
                                  in (state.get('poly_arrangements') or [])],
            'poly_arrangement_index': int(
                state.get('poly_arrangement_index') or 0),
        }

    def _remember_landmark(coords, what, comment):
        """Put a geometry that is *not* the one on screen into the history.

        :func:`_remember` records what is in the box at the moment it is
        called, which is the whole of what an ordinary action needs: it is
        called before the action, and before the action the box holds the
        state to come back to.  A scan cannot do that.  It walks through
        structures the user never chose and will want back -- where it
        started, the highest point it crossed -- and by the time it can name
        them it is standing at the far end of the walk.  So those go in from
        the geometry rather than from the box.

        Everything else about the state comes from where it stands now, which
        is right: a scan changes coordinates and nothing else.  A landmark
        that repeats the entry before it is not added -- a walk that only ever
        went downhill has its highest point at its first step, and two
        identical entries would be a press of Undo that appears to do nothing.

        *comment* is the line the box carries when the landmark comes back, so
        it says where the user is standing.  It has to be one this editor
        recognises as its own, or the next edit would keep it as though the
        user had typed it.
        """
        rows = [line for line in str(coords).splitlines()[2:] if line.strip()]
        if not rows:
            return False
        entry = dict(_structure_marks(),
                     coords=xyz_document(rows, comment),
                     what=str(what), gesture=None)
        history = list(state.get('history') or [])
        if history and history[-1].get('coords') == entry['coords']:
            return False
        history.append(entry)
        if len(history) > _HISTORY_LIMIT:
            history = history[:1] + history[-(_HISTORY_LIMIT - 1):]
        state['history'] = history
        return True

    def _remember(what, gesture=None):
        """Put the state as it is now into the history, under a name.

        One history for everything the viewer can do, and one entry per
        action. There used to be three -- a stack in the browser for drags, a
        stack here for structural edits, and a single slot for the geometry
        before an optimisation -- and none of them knew about the others. Undo
        walked whichever it found first, so it stopped part way back and there
        was no way to the beginning; and a re-render clears the browser's,
        which is to say every draw threw away every drag before it.

        Called *before* the thing it names happens, because what is worth
        keeping is the state that is about to be left.

        *gesture* names a continuous adjustment rather than a single act.
        Holding an arrow key with Set on turns a dihedral half a degree per
        press, and each press is a change of the structure: recorded one by
        one, a two-second sweep is two hundred entries and Undo creeps back
        through it half a degree at a time. Presses that carry the same
        gesture are therefore one step, ended by any other action -- which is
        what a run of typing is in a text editor, and a sweep is the same
        thing.
        """
        history = list(state.get('history') or [])
        entry = dict(_structure_marks(),
                     coords=coords_widget.value,
                     what=str(what),
                     gesture=gesture)
        last = history[-1] if history else None
        # The same picture *and* everything held with it.  Judged on the
        # picture alone, a Hold was never a step of its own -- it changes no
        # coordinate -- so Undo walked straight past it, wiped it on the way,
        # and reported whatever it did land on.  Hold, then a scan, then Undo
        # took back two actions on one press and named one of them.
        same = last is not None and all(
            last.get(key) == value for key, value in entry.items()
            if key not in ('what', 'gesture'))
        carrying_on = (gesture is not None and last is not None
                       and last.get('gesture') == gesture)
        if same or carrying_on:
            # Two actions from the same picture are one step back, and the
            # step is named after the first of them: going back to that state
            # undoes everything since, and the earliest is what the user would
            # look for. Overwriting the name here lost "the structure as it
            # was loaded" the first time anything happened at all.
            pass
        else:
            history.append(entry)
        # The first entry is the structure as it arrived and is never dropped:
        # it is what Reset goes back to, and a long session must not lose it.
        if len(history) > _HISTORY_LIMIT:
            history = history[:1] + history[-(_HISTORY_LIMIT - 1):]
        state['history'] = history
        # A new action makes the way forward unreachable, which is what every
        # editor does and what people rely on without knowing they do: undo
        # three things, do a fourth, and the three are gone rather than
        # waiting to be redone on top of a structure they were never part of.
        # Here rather than in each handler, because this is the one place
        # every structure-changing action already passes through.
        state['structure_undo'] = []
        _refresh_undo_redo()

    def _stop_what_is_running():
        """End everything in flight: it is about a structure that has gone.

        Undo and Reset are the two places where the structure on screen is
        replaced by one the user chose rather than one a calculation made, and
        both used to leave the calculation running. xtb went on minimising the
        geometry from before the press and wrote its answer into the box a few
        seconds later, over the structure that had just been put back -- so
        the press looked as though it had done nothing at all, and the only
        clue was that the coordinates were an optimisation's rather than the
        ones that had been there.

        Three things had to be told, and the run number tells all three: the
        page plays frames tagged with the run they belong to, the workers
        check it before writing an answer, and a fresh number belongs to no
        run at all. The generation counter goes with it, which is what an
        armed settle stands down on.

        Returns the sentence to add to the line saying so, or ``''``.
        """
        _gfn_new_generation()
        # Which frame the picture stood on belonged to a run that is over.
        # Left behind it is a plausible index into a path nobody is walking,
        # and the interrupted run would write that frame over the structure
        # being put back.
        _forget_the_shown_frame()
        running = (state.get('optimize_run') is not None
                   or state.get('climb_run') is not None)
        scanning = bool(state.get('scan_run'))
        if scanning:
            state['scan_stop'] = True
        if not _interrupt_gfn():
            # No optimisation to end, but a settle or a scan may still have
            # frames to push. A run number the page has never seen resets the
            # player and makes every one of them stale.
            #
            # Through the one function that hands them out, which also drops a
            # path a Stop put down and the page has not yet cut: an undo is
            # the one gesture where a geometry from before it must not arrive
            # afterwards, and that path would arrive on the page's next word.
            blank = _claim_the_frame_run()
            submit_gfn_frame.value = _frame_payload(blank, frames=[])
        else:
            # The switch goes back up with the run. It is not coming back by
            # itself -- nothing restarts it, because what replaced the
            # structure was the user asking for an older one.  That goes for
            # the climb's mark too; its toggle is a mode and stays where the
            # user set it, but an undo is not a hand pointing somewhere, so
            # there is nothing to resume and nothing to aim along.
            state.pop('optimize_interrupted', None)
            state.pop('climb_interrupted', None)
            state.pop('climb_was', None)
            state.pop('climb_held', None)
            for switch in (submit_optimize_btn, submit_optimize_all_btn):
                if switch.value:
                    switch.value = False
        if running or scanning:
            return (' The calculation that was running was about the '
                    'structure you took back, so it stopped.')
        return ''

    def _restore(entry, note):
        """Put a remembered state back on screen, and stop what is not."""
        aside = _stop_what_is_running()
        state['structure_edit_inflight'] = True
        _mark_structure_edit()
        # Before the write, because the write rebuilds the picture and the
        # hand corrections re-applied to it have to be the ones belonging to
        # the structure coming back.
        state['hand_bonds'] = dict(entry.get('hand_bonds') or {})
        try:
            coords_widget.value = entry['coords']
        finally:
            state['structure_edit_inflight'] = False
        state['bond_edits'] = dict(entry.get('bond_edits') or {})
        state['constraints'] = [dict(one)
                                for one in (entry.get('constraints') or [])]
        state['hyb_overrides'] = dict(entry.get('hyb_overrides') or {})
        state['poly_applied'] = entry.get('poly_applied')
        state['poly_metal'] = entry.get('poly_metal')
        # Copied out rather than handed over: what the history holds has to
        # survive being put back and changed again, or the second Undo returns
        # what the first one restored.
        assignment = entry.get('poly_assignment')
        state['poly_assignment'] = dict(assignment) if assignment else None
        state['poly_arrangements'] = [dict(one) for one
                                      in (entry.get('poly_arrangements') or [])]
        state['poly_arrangement_index'] = int(
            entry.get('poly_arrangement_index') or 0)
        # The dropdown as well, quietly: it is what says which polyhedron is
        # being held, and left showing the one that was just taken back it
        # would put it on again at the next thing the user touched.
        state['poly_quiet'] = True
        try:
            submit_poly_dd.value = entry.get('poly_applied') or ''
        except Exception:
            pass
        finally:
            state['poly_quiet'] = False
        state['perceived'] = None
        state['perceived_for'] = None
        _refresh_constraints()
        _refresh_poly_turn()
        _set_mol_status(note + aside)
        if submit_relax_btn.value or state.get('ff_bootstrap_done'):
            _enable_live_forcefield()

    def _keep_for_redo(what, back):
        """Put the state Undo is leaving on the way forward.

        *back* is the history entry Undo took off the stack, kept whole so
        that Redo can put the history back exactly as it was rather than
        rebuild something that looks like it.  It is ``None`` for the one
        Undo that takes nothing off -- the walk back to the structure as it
        was loaded, which leaves the first entry standing -- and then Redo
        puts nothing back either, or the history would grow an entry per
        round trip and Undo would start returning what it had just restored.
        """
        forward = list(state.get('structure_undo') or [])
        forward.append(dict(_structure_marks(),
                            coords=coords_widget.value,
                            what=str(what or 'the last step'),
                            back=back))
        if len(forward) > _HISTORY_LIMIT:
            forward = forward[-_HISTORY_LIMIT:]
        state['structure_undo'] = forward

    def _refresh_undo_redo():
        """Redo is offered only when there is something to go forward to."""
        submit_manip_redo_btn.disabled = not (state.get('structure_undo')
                                              and not submit_manip_undo_btn.disabled)

    def _undo_structure():
        """One step back, whatever kind of step it was.

        A drag, a drawn atom, a bond, an optimisation, switching the
        continuous relaxation on: each is one entry, so pressing Undo four
        times walks back through four actions rather than stopping at the
        first kind it does not recognise.
        """
        history = list(state.get('history') or [])
        if len(history) < 1:
            _set_mol_status('Nothing left to undo.')
            return
        # The last entry is the state before the most recent action, so going
        # back to it is undoing that action. The first entry stays: it is the
        # structure as it arrived, and Reset is what returns to it.
        at_start = len(history) == 1
        entry = history[0] if at_start else history.pop()
        state['history'] = history
        # What is being left, before it is left: Redo is Undo read backwards,
        # and the state a step is taken away from is the state that step comes
        # back to.
        _keep_for_redo(entry.get('what'), None if at_start else entry)
        left = 0 if at_start else len(history)
        _restore(entry, f'Took back: {entry.get("what") or "the last step"}.'
                        + (f' {left} more to go back through.' if left
                           else ' That is the structure as it was loaded.'))
        _refresh_undo_redo()

    def _redo_structure():
        """One step forward again, through what Undo took back.

        Excel and Word, exactly: Undo puts what it leaves on the way forward,
        Redo takes it off again and puts the step back on the way back, and
        any new action clears the way forward -- see :func:`_remember`, which
        is where that last one lives because every action passes through it.

        It stops what is running for the same reason Undo does: a run started
        on the structure Redo has just replaced would write its answer over
        the one that was put back a second or two later, which from outside is
        a button that does nothing.  :func:`_restore` is the one place that
        knows how, so both presses go through it.
        """
        forward = list(state.get('structure_undo') or [])
        if not forward:
            _set_mol_status('Nothing to redo. Undo something first, and this '
                            'puts it back.')
            return
        entry = forward.pop()
        state['structure_undo'] = forward
        back = entry.get('back')
        if back is not None:
            # The very entry Undo took off, back where it was, so a second
            # Undo takes this step away again rather than one before it.
            state['history'] = list(state.get('history') or []) + [back]
        left = len(forward)
        _restore(entry, f'Put back: {entry.get("what") or "the last step"}.'
                        + (f' {left} more to go forward through.' if left
                           else ''))
        _refresh_undo_redo()

    def _push_hand_bonds():
        """Put the bonds the user drew or cut back into a rebuilt picture.

        A rebuild draws what the distances say, and a bond drawn between two
        atoms that are not within bonding distance is not in them. So the first
        drawn bond vanished from the view the moment a second edit rebuilt it,
        while remaining in force everywhere else -- it was still there in the
        force field, and an optimisation pulled the two atoms together and made
        it visible again, which is exactly how it looked from outside: the bond
        was gone, and then it was back.

        Only what was corrected by hand, never the whole bond list: what
        perception found is the viewer's own business and is left alone.
        """
        hand = state.get('hand_bonds') or {}
        if not hand:
            return
        triples = [
            [int(i), int(j), 1 if connected else 0]
            for (i, j), connected in hand.items()
        ]
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.applyBondEdits('
            f'{json.dumps(submit_scope_id)},{json.dumps(triples)});'
        )

    def _push_bond_orders(bonds=None):
        """Let the picture show what the bonds are.

        A model read from an XYZ block has no orders in it -- the format
        carries none -- so every bond was drawn as one stick whatever it was.
        3Dmol draws a double as two cylinders and a triple as three once the
        model knows, so the orders are handed over after every render.
        """
        if bonds is None:
            xyz = _current_xyz()
            if not xyz:
                return
            try:
                perceived = _perception_for(xyz)
                from .molecule_forcefield import _orders_from_mol
                known = _orders_from_mol(perceived.typing_mol)
                bonds = {
                    pair: int(known.get(pair, 1)) for pair in perceived.bonds
                }
            except Exception:
                return
        triples = [
            [int(i), int(j), int(order)]
            for (i, j), order in dict(bonds).items() if int(order) > 1
        ]
        if not triples:
            return
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setBondOrders('
            f'{json.dumps(submit_scope_id)},{json.dumps(triples)});'
        )

    def _mark_structure_edit():
        """Tell the re-render that this is an edit, not a different molecule.

        Two things follow from it: the camera stays where the user put it, and
        the continuous relaxation picks up again with the atom that was just
        drawn in it -- which is the point of being able to draw while it runs.
        """
        _ensure_manip_bootstrap()
        _run_manip_js('window.__delfinStructureEdit = true;')

    def _current_xyz():
        """The structure the editor is working on, header and all.

        The Submit tab keeps it in ``state['current_xyz_for_copy']`` because
        its copy button wants it there too, and every force field here read it
        straight out of that key. A second tab has no reason to fill it, so in
        the ORCA Builder the editor believed there was no structure at all and
        Optimize, Dynamik Opt and Settle all answered "load a structure" over a
        molecule that was plainly on screen.

        So the coordinates box the editor was given is asked when the key is
        empty -- that box is the one thing every host has to keep pointed at
        what is being shown. Only when it holds a real XYZ header, since the
        Submit tab's box also holds SMILES, and a SMILES string read as
        coordinates is worse than no answer.
        """
        held = (state.get('current_xyz_for_copy') or {}).get('content')
        if held:
            return held
        lines = [line for line in (coords_widget.value or '').split('\n')
                 if line.strip()]
        if len(lines) < 3:
            return None
        head = lines[0].split()
        if len(head) != 1 or not head[0].isdigit():
            return None
        return coords_widget.value

    def _structure_now():
        from .molecule_builder import structure_from_xyz

        xyz = _current_xyz()
        if not xyz:
            return None
        return structure_from_xyz(xyz, state.get('bond_edits') or {})

    def on_submit_cmd(change):
        """Carry out a gesture the browser cannot finish on its own.

        The value is ``verb:serial:payload``; the serial only exists to make
        the same command twice in a row read as two changes. Placing an atom
        is cheap, but how many hydrogens it needs and where they go is decided
        here, where RDKit's valences and covalent radii are.
        """
        if change.get('name') != 'value':
            return
        # Recorded before it is acted on, and whatever it says: a message this
        # handler refuses is still a message the page sent, and a report that
        # showed only the ones the kernel understood would hide exactly the
        # kind of defect where the page says something the kernel does not.
        record('cmd', v=submit_cmd_sync.value or '')
        parts = (submit_cmd_sync.value or '').strip().split(':')
        if len(parts) != 3:
            return
        verb, payload = parts[0], parts[2]

        if verb == 'editor':
            # The page saying which editor it is running. Until it has, every
            # rendered structure carries a copy of the editor with it.
            state['manip_seen_version'] = str(payload)
            # And saying that a viewer has just been made ready, which is the
            # moment a running force field has to be worked out again. The
            # page clears the parameters it finds when a viewer appears -- they
            # were assigned for the geometry that just went away -- so anything
            # sent before this arrives is cleared with them, and Dynamik Opt
            # stood on with nothing running underneath it.
            if submit_relax_btn.value and not state.get('ff_reassigning'):
                state['ff_reassigning'] = True
                try:
                    _enable_live_forcefield()
                finally:
                    state['ff_reassigning'] = False
            return

        if verb == 'gfnplay':
            # What the playback is doing, said by the page.  Without this the
            # only way to tell an invisible trajectory from a missing one was
            # to read the browser's console.
            state['gfn_play_note'] = str(payload)
            if str(payload).startswith('stopped at frame '):
                # "stopped at frame 12 of run 4": the count and the walk it
                # counts along, because one without the other names nothing.
                words = str(payload).split()
                _keep_the_shown_frame(
                    words[3] if len(words) > 3 else None,
                    words[6] if len(words) > 6 else None)
                # The other half of a Stop.  The run that was stopped has put
                # its path down and cannot cut it, because only the page knows
                # where the picture got to; this is the page saying so.
                # Whichever of the two arrives second lands the geometry, so
                # the answer does not depend on a race between a browser and
                # an optimiser -- and the two run at opposite speeds, a climb
                # finishing inside ten milliseconds and an xtb round taking
                # seconds.
                _land_the_stopped_frame()
            # Drawn again rather than written.  The page reporting on frames
            # is not a result and has no sentence of its own to put in the
            # row: what it has is a fault, and a fault goes on the end of
            # whatever is already there -- see :func:`_trajectory_fault`.
            #
            # Written, it wrote the *optimiser's* last sentence, which is how
            # a scan's verdict came to be replaced by "Optimised with
            # GFN2-xTB..." from a run that had finished before the scan
            # started.  A message that carries no news must not be able to
            # take the row away from one that does.
            _render_mol_status()
            return

        if verb == 'grabbed':
            # A new grab is a new question: whatever the bonding wall
            # refused about the last one says nothing about this one.
            _forget_topology_refusals()
            # A hand on the structure. The step is recorded here rather than
            # when it is let go of, because by then the coordinate box already
            # holds what the drag made -- the relaxation pushes into it while
            # the mouse is still down -- and "before the drag" would not be
            # anywhere any more.
            _remember('a drag')
            return

        if verb == 'gfngrab':
            # An atom has been picked up while xtb was minimising.  Ending the
            # run at the grab rather than at the release is the whole point: a
            # GFN2 run is thirteen seconds, and all of them would otherwise be
            # spent on a structure the user is in the middle of changing.
            #
            # Which frame the picture stood on comes with the message, and it
            # is what the run is cut at: the geometries past it were computed
            # for a structure the hand is changing, and nobody has seen them.
            frame, _, walk = str(payload).strip().partition(',')
            _keep_the_shown_frame(frame, walk)
            # Whichever of the two was walking, and the same sentence for
            # both: a hand has arrived, so the run under it is about a
            # structure that has stopped existing.  Set at the grab rather
            # than at the first drag-follow, because a follow message is a
            # tenth of a second away and the climb makes a step every ten
            # milliseconds -- ten steps away from the hand before anything
            # told it a hand had arrived.
            climbing = state.get('climb_run') is not None
            if _interrupt_gfn():
                _set_mol_status(
                    'Moved while it ran; the '
                    + ('climb' if climbing else 'optimisation')
                    + ' stops there and starts again from what you make.',
                    spinner=True)
            # Where the structure stood when the hand arrived.  The difference
            # between this and where it is let go is the direction the user
            # asked for, and that direction -- not the pull -- is what guides
            # the climb.
            #
            # The climb's own last frame, but only while one was really
            # walking: it outlives the run that made it, and a hand can arrive
            # when no climb is running at all -- one converged a moment ago, or
            # a drag cut one off and Auto was down, and this gesture is what
            # points the next one.  Read unconditionally, that stale frame
            # would be the far end of a direction the user never pointed in.
            # Otherwise the box, which is where the picture is: a hand can also
            # arrive before the first Hessian has finished.
            if _climb_owns_the_release():
                state['climb_was'] = (
                    (state.get('climb_showing') if climbing else None)
                    or _current_xyz())
                # And the pair the last gesture was about is not this one's.
                # Kept, a climb after a drag that names no contact -- a turn,
                # or a hand that never moved -- would be checked against the
                # bond somebody dragged a minute ago, and a check against the
                # wrong pair is worse than no check: it is the ladder walking
                # past an answer that was right.
                state.pop('climb_held', None)
            _begin_gfn_follow()
            # The leash goes on before the hand has moved anything.  It used
            # to be set from the first follow answer, and the first answer is
            # a tenth of a second away -- by which time the mouse has taken
            # the atom an angstrom and the bond it was in is gone.  Waiting
            # for the calculation to say where the structure may stand meant
            # the one stretch that is never checked was the one that does the
            # damage.
            #
            # The structure on screen needs no checking: it is where the
            # budget was last agreed, by definition.  Every atom is marked,
            # not only the ones being dragged, because which those are is not
            # known until the first drag-follow says so -- and a mark on an
            # atom nobody is moving costs nothing.
            _arm_thermal_leash()
            return

        if verb == 'gfnfree':
            _end_gfn_follow()
            # Let go of.  Whether that carries on, and which way, is the Auto
            # and Climb to TS switches -- and it is asked in one place so the
            # answer cannot depend on whether a run happened to be going when
            # the atom was picked up, or on which of the two it was.
            _after_release()
            return

        if verb == 'undo':
            _undo_structure()
            return

        if verb == 'redo':
            # The keyboard's way in.  The browser keeps no way forward of its
            # own -- its snapshot stack is cleared by every re-render -- so
            # Redo is always this history's, whichever key asked for it.
            state.pop('pre_optimize_frames', None)
            _redo_structure()
            return

        if verb == 'unbond':
            indices = [int(p) for p in payload.split(',') if p.strip().isdigit()]
            if len(indices) == 2:
                state['picked'] = sorted(indices)
                _edit_bond(False)
            return

        from .molecule_builder import (
            delete_atoms, grow_from, normalise_element, place_atom,
            set_bond_order, set_element,
        )

        structure = _structure_now()
        if structure is None:
            _set_mol_status('Load a structure first.')
            return
        fields = payload.split(',')
        # Whether an edit fills or trims the hydrogens around it, or leaves
        # them exactly as they are.
        keep_h = not bool(submit_adjust_h_btn.value)
        try:
            if verb == 'addatom' and len(fields) == 4:
                element = normalise_element(fields[0])
                if element is None:
                    _set_mol_status(f'{fields[0]} is not an element.')
                    return
                place_atom(structure, element,
                           [float(v) for v in fields[1:4]],
                           adjust_h=not keep_h)
                _apply_structure(structure, f'Placed {element}.')
            elif verb == 'grow' and len(fields) in (6, 9):
                element = normalise_element(fields[1])
                if element is None:
                    _set_mol_status(f'{fields[1]} is not an element.')
                    return
                # Where the hand let go, when the page said. With the
                # hydrogens left alone the atom belongs there rather than at
                # the length its bond would prefer -- that is what switching
                # the adjustment off is for.
                landed = [float(v) for v in fields[6:9]] if len(fields) == 9 else None
                grown = grow_from(structure, int(fields[0]), element,
                                  order=int(fields[2]),
                                  direction=[float(v) for v in fields[3:6]],
                                  adjust_h=not keep_h, at=landed)
                if grown is None:
                    _set_mol_status(f'{element} could not be grown there.')
                    return
                _apply_structure(
                    structure,
                    f'Grew {element}.' if not keep_h else
                    f'Grew {element} where you let go.')
            elif verb == 'setelement' and len(fields) == 2:
                element = normalise_element(fields[1])
                if element is None:
                    _set_mol_status(f'{fields[1]} is not an element.')
                    return
                index = int(fields[0])
                was = structure.symbols[index] if index < len(structure) else '?'
                if set_element(structure, index, element, adjust_h=not keep_h):
                    _apply_structure(structure, f'{was}{index} is now {element}.')
            elif verb == 'bondcycle' and len(fields) == 2:
                first, second = (int(v) for v in fields)
                ends = [structure.symbols[i] for i in (first, second)]
                current = structure.order(first, second)
                stepped = None
                for step in (1, 2, 3):
                    candidate = (current - 1 + step) % 3 + 1
                    if candidate == current:
                        break
                    if set_bond_order(structure, first, second, candidate,
                                      adjust_h=not keep_h):
                        stepped = candidate
                        break
                if stepped is None:
                    _set_mol_status(
                        f'{ends[0]}{first}-{ends[1]}{second} can only be '
                        'single: neither end has valence for more.'
                    )
                else:
                    named = {1: 'single', 2: 'double', 3: 'triple'}[stepped]
                    _apply_structure(
                        structure,
                        f'{ends[0]}{first}-{ends[1]}{second} is now {named}.')
            elif verb == 'bondorder' and len(fields) == 3:
                first, second, order = (int(v) for v in fields)
                named = {1: 'single', 2: 'double', 3: 'triple'}.get(order, '')
                ends = [structure.symbols[i] for i in (first, second)]
                if not structure.order(first, second):
                    # A bond that is not there yet is made the way the Bond
                    # button makes one: the topology changes and nothing else.
                    # Moving fragments and re-placing hydrogens to go with it
                    # was more than was asked for, and it wrecked the molecule
                    # on the way -- while Bond, which does none of that, has
                    # always worked.
                    state['picked'] = [first, second]
                    _edit_bond(True)
                elif set_bond_order(structure, first, second, order,
                                    adjust_h=not keep_h):
                    _apply_structure(
                        structure,
                        f'{ends[0]}{first}-{ends[1]}{second} is now {named}.')
                else:
                    _set_mol_status(
                        f'{ends[0]}{first}-{ends[1]}{second} cannot be '
                        f'{named}: one of them has no valence left for it.'
                    )
            elif verb == 'delatoms':
                doomed = [int(p) for p in fields if p.strip().lstrip('-').isdigit()]
                gone = delete_atoms(structure, doomed, adjust_h=not keep_h)
                if gone:
                    _apply_structure(structure, f'Deleted {gone} atom(s).')
        except Exception as exc:
            _set_mol_status(f'That edit did not work: {exc}')

    def on_submit_centre(_button=None):
        """Put the system back in the middle of the picture.

        Nothing about the structure changes -- this is the camera and only the
        camera. It is a button rather than something that happens by itself
        because a view is the user's: moving in on one corner is a thing
        people do on purpose, and a picture that re-frames itself while they
        work is worse than one they have to bring back now and then.
        """
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            f'window.__delfinSubmitManip.recentreView({json.dumps(submit_scope_id)});'
        )
        _set_mol_status('Back in the middle. The structure is unchanged.')

    def on_submit_bond(_button=None):
        _edit_bond(True)

    def on_submit_unbond(_button=None):
        _edit_bond(False)

    def on_submit_swap(_button=None):
        """Exchange two ligands outright rather than dragging one at another.

        The two arrangements are separate minima and the relaxation only runs
        downhill, so it can never cross between them: a ligand dragged part of
        the way simply rolls back. The exchange is therefore performed in one
        step and the field is left to tidy up afterwards.
        """
        indices = list(state.get('picked') or [])
        metal = state.get('poly_metal')
        if metal is None:
            perceived = state.get('perceived')
            metals = list(getattr(perceived, 'metal_indices', ()) or ())
            metal = metals[0] if len(metals) == 1 else None
        if len(indices) != 2 or metal is None:
            _set_mol_status('Select two ligands of one metal to exchange them.')
            return
        # The line below promises that Undo puts them back, and until this was
        # here it did not: the exchange happens in the browser, arrives as a
        # geometry, and nothing recorded the arrangement it replaced -- so
        # Undo took back whatever had been done before the swap instead.
        _remember(f'exchanging the ligands at {indices[0]} and {indices[1]}')
        state['poly_assignment'] = None
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.exchangeLigands('
            f'{json.dumps(submit_scope_id)},{int(metal)},'
            f'{int(indices[0])},{int(indices[1])});'
        )
        _set_mol_status(
            'Exchanged the two ligands. The field is settling the result; '
            'Undo puts them back.'
        )

    #: How many rounds of relaxation each point of a scan may take.
    #:
    #: More than a follow step, because a scan is not trying to keep up with a
    #: hand: it is answering the question the drag can only estimate, and a
    #: relaxed scan cut short reads as a barrier that is not there.  Measured
    #: on a Diels-Alder approach whose relaxed cost is +3.6 kcal/mol: three
    #: cycles say +40.0, six +18.9, twenty +8.3.  Each point starts from the
    #: geometry the last one reached, so the work is spread along the path
    #: rather than repeated at every step -- which is why a dense scan is
    #: affordable at all, and a dense scan is the point.
    _SCAN_CYCLES = 60

    #: When a scan has arrived somewhere and should stop.
    #:
    #: A scan that runs past the next minimum stops describing a reaction and
    #: starts pushing into a structure.  Measured on a tert-butyl bromide with
    #: a chloride: the path crosses at +25.2 kcal/mol, falls to -4.7 in the
    #: product well, and the steps after that squeezed the new C-Cl to 1.28 A
    #: -- a number that means nothing, on a geometry the hold no longer
    #: determined.
    #:
    #: Where the next minimum lies cannot be known in advance, but it can be
    #: recognised.  Not by the path going flat: a scan keeps driving its
    #: coordinate, so past the well it climbs the far wall and never flattens
    #: at all -- that rule never fired.  By the climb itself.  Over the top,
    #: down into something, and then rising again for two steps running: the
    #: bottom of that is the minimum, and it is behind us.
    _SCAN_OVER_THE_TOP = 2.0      # kcal/mol below the highest point
    _SCAN_CLIMBING = 2            # steps of rising again that end it
    _SCAN_UPHILL = 0.2            # kcal/mol a step counts as rising

    def _climbed(spent):
        """The top of the path, and how far it stands above what came before.

        A barrier is a rise out of a minimum, and the first point of a path is
        not necessarily one.  Measured against the first point instead, a walk
        that starts strained reports a barrier smaller than it is -- and a
        push, whose first point relaxes what the hand posed, reports none at
        all: every later point sits below the start and the highest point on
        the path is the start itself, which came back as "it wants about 3 K",
        the temperature of no barrier.

        So the rise is measured from the lowest point *before* the top, which
        is the minimum the path actually left.
        """
        top = max(spent)
        where = spent.index(top)
        return top, top - min(spent[:where + 1])

    def _descent(summit, bottom, spent, geometry, value):
        """The highest point so far, and the lowest one since it.

        Returned as a pair, updated by one step at a time, because that is all
        a scan can afford to keep: the geometry of every point of a long walk
        over a large structure is a great deal of memory for a question with a
        two-value answer.

        The reset on a new highest point is the whole of it.  Kept over the
        whole path instead -- which is what this did -- the lowest point of a
        scan that starts in a well and climbs out of it is the *first* one, so
        every scan that found anything handed back the structure it began with
        and called it the minimum it had walked through.  On the page that
        read as "The scan walked 1 of 20 points. Highest +0.0 kcal/mol", which
        is indistinguishable from a scan that did nothing at all.
        """
        here = (spent, geometry, value)
        if summit is None or spent > summit[0]:
            return here, here
        if bottom is None or spent < bottom[0]:
            return summit, here
        return summit, bottom

    def _scan_arrived(path):
        """Whether the walk has crossed something and is climbing out again."""
        if len(path) < _SCAN_CLIMBING + 3:
            return False
        spent = [one[1] for one in path]
        top, rise = _climbed(spent)
        if rise < _SCAN_OVER_THE_TOP:
            return False              # nothing has been crossed yet
        floor = min(spent[spent.index(top):])
        if top - floor < _SCAN_OVER_THE_TOP:
            return False              # still up there
        last = spent[-(_SCAN_CLIMBING + 1):]
        return (all(b - a > _SCAN_UPHILL for a, b in zip(last, last[1:]))
                and spent[-1] - floor > _SCAN_UPHILL)

    def _value_of_constraint(entry):
        """What that coordinate measures on the structure now, or None.

        None when it names atoms this structure does not have -- which is what
        a scan armed on one molecule and run on another does, and what made the
        atom list throw on the next click rather than say so.
        """
        return _value_in(_current_xyz() or '', entry)

    def _value_in(xyz, entry):
        """The same, read off a geometry the caller names.

        A push has to read its coordinate back from what the last force
        actually did, and that geometry is in the middle of a run rather than
        in the box.
        """
        here = _gfn.coordinates_of(xyz or '')
        atoms = [int(i) for i in (entry.get('atoms') or ())]
        if not atoms or any(i < 0 or 3 * i + 2 >= len(here) for i in atoms):
            return None
        at = [(here[3 * i], here[3 * i + 1], here[3 * i + 2]) for i in atoms]
        if entry.get('kind') == 'distance' and len(at) == 2:
            return math.dist(at[0], at[1])
        if entry.get('kind') == 'angle' and len(at) == 3:
            first = [a - b for a, b in zip(at[0], at[1])]
            second = [a - b for a, b in zip(at[2], at[1])]
            one = math.sqrt(sum(v * v for v in first)) or 1.0
            two = math.sqrt(sum(v * v for v in second)) or 1.0
            cosine = sum(a * b for a, b in zip(first, second)) / (one * two)
            return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
        if entry.get('kind') == 'dihedral' and len(at) == 4:
            return _gfn._dihedral(at, 0, 1, 2, 3)
        return None

    def _scan_legs():
        return list(state.get('scan_legs') or [])

    def _leg_atoms_label(leg):
        """The atoms a leg drives, named the way the sentences name them.

        Its own function because the profile's axis wants exactly this and
        nothing else -- ``C0-C1`` -- while :func:`_describe_leg` wants it with
        the walk's two ends after it.  One place, so that the picture and the
        sentence cannot come to call the same pair of atoms different things.
        """
        symbols = []
        perceived = state.get('perceived')
        known = getattr(perceived, 'symbols', None) or ()
        if not known:
            # The perception comes back from the browser, and it has not
            # always been asked for by the time a leg is named -- an editor
            # driven from the outside may never have had a viewer at all.  The
            # element column of the structure in the box answers the same
            # question and is always there, so a leg is named "C0-C10" rather
            # than "?0-?10", on the picture and in the sentence alike.
            known = [line.split()[0] for line in
                     _gfn.atom_lines(_current_xyz() or '')]
        for index in leg['atoms']:
            symbol = known[index] if 0 <= index < len(known) else '?'
            symbols.append(f'{symbol}{index}')
        return '-'.join(symbols)

    def _describe_leg(leg):
        # A leg with a verb is read as the instruction it is rather than as
        # the pair of numbers underneath it: "form C1-C11" is what was asked
        # for, and "3.35 -> 1.53 A" is only how the force will be pointed.
        verb = str(leg.get('verb') or '')
        if verb in ('form', 'break'):
            return f'{verb} {_leg_atoms_label(leg)}'
        unit = 'A' if leg['kind'] == 'distance' else 'deg'
        return (f"{_leg_atoms_label(leg)} {leg['from']:.3g} -> {leg['to']:.3g} "
                f"{unit}")

    def _refresh_scan():
        """Show the armed legs, or nothing at all when there are none."""
        legs = _scan_legs()
        showing = '' if legs else 'none'
        for widget in (submit_scan_dd, submit_scan_del, submit_scan_whole,
                       submit_scan_how, submit_scan_energy, submit_scan_back,
                       submit_scan_run_btn):
            widget.layout.display = showing
            widget.disabled = not legs
        # Except the return leg, which belongs to walking a value and not to
        # pushing one: a push is a ramp of forces rather than a grid of
        # values, so there is no same-coordinate-backwards to walk.
        if legs and str(submit_scan_how.value) == 'push':
            submit_scan_back.layout.display = 'none'
        # The second opinion answers to a different question from the rest of
        # the row, and that is the point of it: everything above is here
        # because a scan is *armed*, and this is here because one has
        # *finished* and left a profile about the molecule on screen.  A press
        # that appeared with the arming would be offering to price a walk that
        # has not happened.
        showing_price = '' if _reprice_is_possible() else 'none'
        submit_scan_price_btn.layout.display = showing_price
        submit_scan_price_btn.disabled = not showing_price == ''
        options = [(_describe_leg(leg), str(n)) for n, leg in enumerate(legs)]
        submit_scan_dd.options = options or [('nothing armed', '')]
        if options:
            submit_scan_dd.value = options[0][1]
        # The two fields belong to the selection, not to the list, so they
        # follow whether something is selected rather than whether anything is
        # armed.
        picked = len(state.get('picked') or ())
        wanted = '' if picked in _CONSTRAINT_KINDS else 'none'
        for widget in (submit_scan_way, submit_scan_steps):
            widget.layout.display = wanted
            widget.disabled = not wanted == ''
        # Two atoms are closer or further apart; three or four are narrower or
        # wider.  The words follow what is picked, so the direction is never a
        # setting with no subject.
        #
        # And the two verbs are there for a pair and absent for an angle or a
        # torsion.  That absence is a statement, which is why it is made here
        # rather than by refusing the press afterwards: a bond is between two
        # atoms, so there is no bond for three of them to make or break.
        kind = _CONSTRAINT_KINDS.get(picked)
        options = (
            [('closer together', 'in'), ('further apart', 'out'),
             ('form this bond', 'form'), ('break this bond', 'break'),
             ('to a value you give', 'to')]
            if kind == 'distance'
            else [('narrower', 'in'), ('wider', 'out'),
                  ('to a value you give', 'to')])
        if list(submit_scan_way.options) != [tuple(one) for one in options]:
            # Rewritten around whatever was chosen.  A dropdown handed a new
            # list drops its value even where the new list still has it, and
            # the three values are the same three in either wording -- so
            # picking a torsion after setting an end silently turned the end
            # back into a direction, and the walk went the other way.
            want = submit_scan_way.value
            submit_scan_way.options = options
            if want in [value for _label, value in options]:
                submit_scan_way.value = want
        # The end of the walk only when one has been asked for.
        set_end = wanted == '' and str(submit_scan_way.value) == 'to'
        submit_scan_to.layout.display = '' if set_end else 'none'
        submit_scan_to.disabled = not set_end

    #: The shortest a scan may ask two atoms to be, as a share of the bond
    #: they would make.
    #:
    #: Not a calibration -- a statement.  A scan drives its coordinate wherever
    #: it is told, so a target under a bond length walks straight into one and
    #: the atoms collapse: asked for 0.60 A between two ring carbons whose bond
    #: is 1.53, it went there.  Slightly under one is real (a bond compresses
    #: in a transition state), well under one is not, and nothing has to be
    #: measured to know which side of that a number is on.
    _SCAN_NO_CLOSER = 0.85

    #: How a push ramps, and when it stops to look closer.
    #:
    #: Geometrically rather than in equal steps, because the range is thirty
    #: to one and what matters is where the reaction goes over: measured on
    #: butadiene and ethylene, 10 kcal/mol/A closes the pair to 2.47 A and 20
    #: is already in the product, so equal steps of six put no point at all on
    #: the barrier.  Twenty points from 4 to 120 put five between 9 and 20.
    #:
    #: And when a step still jumps -- the force finally pays for the crossing
    #: and the structure falls through to the other side -- the force between
    #: the last one and this one is bisected and the step taken again from
    #: where it started.  Without it the highest point on the path is whatever
    #: happened to be sampled before the fall, which reads as a low barrier
    #: rather than an unsampled one.
    #: And when even the finest force still falls through -- which it does,
    #: because a crossing has a threshold and a force either pays for it or
    #: does not -- the segment it fell through is walked afterwards with the
    #: coordinate held, which is the one thing coordinate driving is good at.
    #: Measured on the Diels-Alder: the push goes over between 2.53 A (+4.4
    #: kcal/mol) and 1.54 A (-64.2) in one step, so its own highest point is
    #: +4.4 and the top was never sampled.  Walking those two apart puts a
    #: point at 2.363 A and +6.3, which is the top, and -6.6 at 2.03 on the
    #: way down.  The push found the path and the walk priced it; neither
    #: could have done both.
    _PUSH_JUMP = 0.35             # A in one step that asks for a finer force
    _PUSH_REFINE = 3              # bisections before the step stands
    _PUSH_ACROSS = 6              # points the crossing itself is priced with

    def _scan_floor_for(leg):
        """The shortest this leg may be driven to, or None if it is an angle."""
        if leg['kind'] != 'distance' or len(leg['atoms']) != 2:
            return None
        rows = [line.split() for line in _gfn.atom_lines(_current_xyz() or '')]
        if any(not (0 <= i < len(rows)) for i in leg['atoms']):
            return None
        from delfin.atom_mapping import cov_radius
        reach = sum(cov_radius(str(rows[i][0])) for i in leg['atoms'])
        return _SCAN_NO_CLOSER * reach

    #: How far a scan may walk when nobody has said where to stop.
    #:
    #: It stops at the next minimum, so this is only the emergency brake for a
    #: walk that has no next minimum -- two fragments pulled apart never find
    #: one, they simply separate.  Generous enough that it is never what ends
    #: an ordinary scan.
    _SCAN_AS_FAR_AS = {'distance': 2.5, 'angle': 180.0, 'dihedral': 360.0}

    def _suggest_scan_target(kind, value, way='in'):
        """Where a scan walks when the user has given a direction and no end.

        The end used to be the thing the user typed, and it is the wrong thing
        to ask for: the scan stops at the next minimum anyway, so where it ends
        is the chemistry rather than a number -- and a number typed there is
        how two atoms came to be asked for 0.60 A when the bond between them is
        1.53.  The direction is the one thing that cannot be read off the
        selection, so that is what is asked, and the number stays for the two
        cases that need it: following a figure from the literature, and a
        coupled scan over two coordinates, where it is the ratio of the two
        ends that fixes the path.
        """
        far = _SCAN_AS_FAR_AS.get(kind, 360.0)
        return float(value) + (far if way in ('out', 'break') else -far)

    #: How far past the graph's own threshold a break is pointed.
    #:
    #: :data:`gfn_optimize.BOND_STARTS_AT` is where a bond stops being drawn,
    #: so a break aimed exactly there is aimed at the line it has to be over.
    #: Half an Angstrom past it, which for two carbons is 2.5 A rather than
    #: 2.0 -- and the number matters far less than it looks, because a push
    #: never meets its target: the target only says which way the force
    #: points, and :func:`_carried_out` is what decides the walk is finished.
    _BREAK_CLEAR_BY = 0.5

    def _verb_target(indices, verb):
        """Where a form or a break points its force, in Angstrom.

        A form is pointed at the bond those two atoms would make -- the sum of
        their covalent radii -- and a break at something safely clear of it.
        Neither is a value the walk is driven to: a push is a force, and where
        the structure ends up is the structure's answer.  What the target does
        is give the force a sign, and what it gives the user is a leg that
        reads as the chemistry it names.
        """
        rows = [line.split() for line in _gfn.atom_lines(_current_xyz() or '')]
        if any(not (0 <= i < len(rows)) for i in indices):
            return None
        from delfin.atom_mapping import cov_radius
        reach = sum(cov_radius(str(rows[i][0])) for i in indices)
        if verb == 'form':
            return reach
        return _gfn.BOND_STARTS_AT * reach + _BREAK_CLEAR_BY

    def _carried_out(xyz, legs):
        """Whether every instruction holds on this one geometry.

        The stopping rule for a form/break walk, and the whole of it: read the
        bond graph off what the last force actually produced and ask whether
        every armed verb is satisfied *at the same time*.  On the same
        geometry and not each at some point along the way, because a form can
        be satisfied at one step and undone two steps later by the break
        pulling the same fragment apart -- and a walk that stopped on the
        first would report a structure that no longer has the bond it says it
        made.

        It is the geometric test -- :func:`gfn_optimize.bond_graph`, covalent
        radii, the same one the viewer draws lines with and ``Keep bonds``
        judges against -- so what the rule reads is what the user sees.  It is
        deliberately **not** a bond order, and that was measured rather than
        assumed.

        SCINE's NT2 stops on Mayer bond orders: formed above 0.75, broken
        below 0.15.  Three measurements here, all under GFN2, reading xtb's
        own Wiberg orders.

        *Forming, and the two tests agree.*  Along the converged Diels-Alder
        band a forming C-C reads 0.000 at 2.64 A, 0.193 at 2.33, 0.524 at 2.09
        and 0.920 at 1.70, and the bond graph flips between those last two as
        well.  Same on the SN2's Cl-C, which is the harder case because the
        two radii differ: 0.000 at 2.58 A and 0.921 at 1.78, and the graph
        turns over at 2.31.

        *Breaking heterolytically, and they still agree.*  The same SN2's
        C-Br: 0.881 at 2.01 A, then 0.000 at 3.11.  The pair leaves with the
        bromide, which a restricted determinant describes perfectly well.

        *Breaking homolytically, and the order is simply wrong.*  An ethane's
        C-C stretched rigidly from its 1.5212 A equilibrium:

            1.52 A  1.030      3.00 A  0.958
            2.00 A  0.994      3.20 A  0.954
            2.50 A  0.973      3.50 A  0.913
                               4.00 A  0.264

        At 3.5 A -- two and a third times equilibrium, and a full Angstrom and
        a half past where the graph gave up on it -- the order still reads
        0.91, because a restricted single determinant cannot part an electron
        pair.  "Broken below 0.15" would not fire until about 4.5 A.

        So the rule is not "bond order is wrong"; it is that a bond order is
        right for a heterolytic break and wrong for a homolytic one -- and
        nothing here knows in advance which it has been asked for.  A test
        that is sound for one reaction and silently wrong for the next is
        worse than one that is neither, and the geometry is right in all three
        measurements above.  It also costs nothing: the graph is already
        computed for the topology wall, where an order would be another xtb
        file to read on every point.  One test for both halves, and it is the
        one the picture is drawn with.
        """
        graph = _gfn.bond_graph(xyz)
        for leg in legs:
            verb = str(leg.get('verb') or '')
            if verb not in ('form', 'break'):
                continue
            pair = tuple(sorted(int(i) for i in leg['atoms']))
            if (pair in graph) != (verb == 'form'):
                return False
        return True

    def on_submit_scan(_button=None):
        """Arm the value the selection describes as a leg of the scan."""
        indices = list(state.get('picked') or [])
        kind = _CONSTRAINT_KINDS.get(len(indices))
        if not kind:
            _set_mol_status('Pick 2, 3 or 4 atoms before arming a scan.')
            return
        here = float(submit_internal_value.value)
        # A value only when one was asked for.  Read out of the field either
        # way, an empty field was a target of zero, and a scan told to walk a
        # distance to zero walks it into a collision -- which is what the
        # floor below then had to catch.
        target = _suggest_scan_target(kind, here, submit_scan_way.value)
        if str(submit_scan_way.value) == 'to':
            asked = float(submit_scan_to.value)
            # A value says which way the walk goes all by itself, so there is
            # no direction left to fall back on when the value is the one the
            # coordinate already has.  Said rather than guessed: the guess
            # used to be inwards, and a scan that walks the opposite way from
            # the one the number implied is worse than no scan.
            if abs(asked - here) <= 1e-9:
                _set_mol_status(
                    f'{here:.3g} is where that coordinate already is, so '
                    'there is nothing to walk. Give another value, or choose '
                    'a direction and let the scan stop at the next minimum.')
                return
            target = asked
        # A verb points its own force, at the bond those two atoms would make
        # or at somewhere clear of it -- see :func:`_verb_target`.  It is not
        # where the walk is driven to; a push is a force, and the verb is what
        # says when the walk is finished.
        verb = str(submit_scan_way.value) if str(
            submit_scan_way.value) in ('form', 'break') else ''
        if verb:
            aimed = _verb_target(indices, verb)
            if aimed is None:
                _set_mol_status('That selection names atoms this structure '
                                'does not have.')
                return
            target = aimed
            # An instruction that is already carried out is not one, and it is
            # refused against the same test that will decide the walk is
            # finished -- :func:`_carried_out`.  Measured with a comparison of
            # lengths instead: an ethane's C-C is 1.521 A and the bond those
            # two carbons make is 1.520, so "form" on a bond slipped through
            # by a thousandth of an Angstrom and armed a walk that was over
            # before it started.  One rule, asked once, and there is no
            # thousandth to fall through.
            if _carried_out(_current_xyz() or '', [{'kind': 'distance',
                                                    'atoms': indices,
                                                    'verb': verb}]):
                named = _leg_atoms_label({'atoms': indices})
                _set_mol_status(
                    f'{named} is already '
                    + (f'bonded, at {here:.3g} A, so there is nothing to '
                       'form. Pick a pair that is not bonded yet.'
                       if verb == 'form' else
                       f'{here:.3g} A apart, which is not a bond, so there is '
                       'nothing to break.'))
                return
        # Same atoms in either order.  It matters here in a way it did not
        # before: with a verb, arming "break 11-1" over "form 1-11" would
        # otherwise leave both on the list, and a walk cannot make and break
        # the same bond at once.
        pair = sorted(indices)
        legs = [one for one in _scan_legs()
                if sorted(one['atoms']) != pair or len(one['atoms']) != len(
                    indices)]
        leg = {'kind': kind, 'atoms': indices, 'from': here,
               'to': target, 'steps': int(submit_scan_steps.value),
               'structure': _structure_fingerprint(_current_xyz() or '')}
        if verb:
            leg['verb'] = verb
        floor = _scan_floor_for(leg)
        clipped = ''
        if floor is not None and leg['to'] < floor:
            clipped = (f' Asked for {leg["to"]:.3g}, which is inside the bond '
                       f'those two would make, so it stops at {floor:.3g}.')
            leg['to'] = floor
        legs.append(leg)
        state['scan_legs'] = legs
        _refresh_scan()
        # An instruction is described as one, and what it is armed with is the
        # other half of the sentence: a form or a break is carried out by a
        # push, so a walk cannot do it and saying so here saves a refusal at
        # the press.
        if verb:
            others = [one for one in legs if one is not leg]
            _set_mol_status(
                f'Armed: {_describe_leg(leg)}. '
                + ('Arm the other half on a second pair -- one bond made '
                   'while another breaks is what most reactions are -- or '
                   'press Run scan.' if not others else
                   'Together: ' + ', '.join(_describe_leg(one) for one in legs)
                   + '.')
                + (' Set to "push with a force"; a walk drives a value and a '
                   'verb needs a force.'
                   if str(submit_scan_how.value) != 'push' else ''))
            _clear_selection()
            return
        _set_mol_status(
            f'Armed {_describe_leg(legs[-1])} in {legs[-1]["steps"]} steps. '
            + ('Arm another to walk them together, which is what a concerted '
               'step needs, or press Run scan.' if len(legs) == 1 else
               f'{len(legs)} legs, walked together.') + clipped)
        _clear_selection()

    def on_submit_scan_del(_button=None):
        legs = _scan_legs()
        try:
            which = int(submit_scan_dd.value)
        except (TypeError, ValueError):
            return
        if 0 <= which < len(legs):
            gone = legs.pop(which)
            state['scan_legs'] = legs
            _refresh_scan()
            _set_mol_status(f'Dropped {_describe_leg(gone)}.')

    def on_submit_scan_run(_button=None):
        """Walk every armed leg together and say what the path costs."""
        legs = _scan_legs()
        xyz = _current_xyz()
        method = str(submit_ff_dd.value)
        if not legs or not xyz:
            return
        # Stopping first.  Behind the method check, a scan started under GFN2
        # and then switched away could not be stopped at all: the button read
        # Stop, was enabled, and answered "a scan needs xtb" while the run it
        # was meant to end walked on to the last point and overwrote the box.
        if state.get('scan_run'):
            state['scan_stop'] = True
            _set_mol_status('Stopping the scan after this point.')
            return
        if not _gfn.is_gfn_method(method):
            _set_mol_status('A scan needs xtb: choose a GFN method.')
            return
        pushing = str(submit_scan_how.value) == 'push'
        # A verb is carried out by a force and cannot be carried out by a
        # walk.  A walk drives its coordinate to a value it was told, so
        # "break C2-Br3" as a walk is the editor pulling the bromide to 3.0 A
        # whatever the molecule thinks -- which is not the reaction, it is a
        # picture of one.  A push ramps until the structure gives, and what
        # gives is the structure's answer.
        instructed = [one for one in legs
                      if str(one.get('verb') or '') in ('form', 'break')]
        if instructed and not pushing:
            _set_mol_status(
                'Form and break are carried out by a force, so they need '
                '"push with a force". A walk drives the distance to a number '
                'instead, which draws the reaction rather than finding it.')
            return
        # xtb takes one force constant for the whole $constrain block, and a
        # push is three orders of magnitude softer than a hold.  Run together,
        # the hold's stiffness would win and the push would silently become an
        # ordinary scan -- the same number, the same picture, a different
        # experiment.  Said instead, with the two ways out of it.
        if pushing and (state.get('constraints') or []):
            _set_mol_status(
                'A push and a held value cannot share one force constant in '
                'xtb, and the hold would win. Drop what is held, or choose '
                '"walk the value".')
            return
        # And a push is a force between two atoms.  There is no force that
        # drives an angle towards a reaction: what it would push on is the
        # bond lengths that make the angle, which is what the two atoms at its
        # ends already say.
        if pushing and any(one['kind'] != 'distance' for one in legs):
            _set_mol_status(
                'A push is a force between two atoms, so it walks distances. '
                'Arm the distance that is forming or breaking, or choose '
                '"walk the value" for an angle or a torsion.')
            return
        # And the one question this method has no way of answering.
        #
        # GFN-FF perceives its bonding once and holds it, so a scan that
        # drives two atoms together across the line where they would be
        # bonded is asking a force field with no term for that bond what
        # making it costs.  Said here in the way an unparametrised solvent is
        # -- before the run, because xtb does not refuse it: it converges, and
        # reports repulsion as a barrier.  See :func:`_gfn.gfnff_refusal` for
        # the measurement and for why breaking is a different case.
        if str(method).strip().lower() == 'gfnff':
            no = _gfn.gfnff_refusal(xyz, legs)
            if no:
                _set_mol_status(no)
                return
        # From where the structure is *now*, not from where it was when the
        # leg was armed.  Left as armed, a second press walked from a value
        # the molecule no longer had: measured, a C-C at 4.012 was compressed
        # to 2.137 in one step and the run reported a 77 kcal/mol barrier with
        # a temperature and a timescale attached, all of it invented.
        fresh = []
        for one in legs:
            now = _value_of_constraint(one)
            if now is None:
                _set_mol_status(
                    f'The scan cannot start: {_describe_leg(one)} names atoms '
                    'this structure does not have. Drop it and arm it again.')
                return
            fresh.append(dict(one, **{'from': now}))
        legs = fresh
        state['scan_legs'] = legs
        _refresh_scan()
        steps = max(2, min(int(one['steps']) for one in legs))
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        model = _solv_model()
        # One entry, before anything moves.  A scan rewrites the structure
        # from end to end, and it went into no history at all: Undo stepped
        # back past it to whatever had been done before, and Reset -- which
        # returns to the structure as it was loaded -- had had that structure
        # replaced underneath it.  There was no way back to where the scan
        # started from.  Every step in this editor has to be one entry, and a
        # scan is a step.
        _remember(f'a scan of {_describe_leg(legs[0])}'
                  + (f' and {len(legs) - 1} more' if len(legs) > 1 else ''))
        state['scan_run'] = True
        state['scan_stop'] = False
        state['scan_arrived'] = False
        state['scan_came_back'] = None
        state['scan_free'] = None
        state['scan_ends'] = None
        state['scan_gap_first'] = None
        state['scan_gap_least'] = None
        state['scan_depth'] = ''
        state['scan_crowded'] = None
        state['scan_free_shaky'] = None
        state['scan_stopped_out'] = False
        # Whether a second leg was even on the table.  None for a push, which
        # has no grid of values to retrace; True or False for a walk, which
        # has.  It is what lets the verdict offer the return leg to the one
        # person who could have had it and did not.
        state['scan_back_wanted'] = (None if pushing
                                     else bool(submit_scan_back.value))
        # The two legs as they are walked, the return leg's verdict, and the
        # step the walk fell through if it fell through one.  Kept on the
        # state rather than only in the sentence, because a profile is a
        # picture before it is a paragraph and both legs belong on the same
        # axes -- see :func:`_scan_two_legs` for the shape and who reads it.
        state['scan_there'] = []
        state['scan_back'] = []
        state['scan_disagree'] = None
        state['scan_jumped'] = None
        # The walk that is starting replaces the one that finished, and its
        # second opinion with it.  Left standing, the press beside Run scan
        # would offer to re-price a profile that had just been walked over.
        state['scan_walk'] = None
        state['scan_repriced'] = None
        _refresh_scan()
        state['scan_carried_out'] = None
        state['scan_instructed'] = [_describe_leg(one) for one in instructed]
        state['scan_frame_run'] = _note_the_run(
            int(state.get('gfn_run', 0)) + 1, 'scan')
        state['gfn_run'] = state['scan_frame_run']
        _ensure_manip_bootstrap()
        schedule_ui_update(_install_gfn_frame_watcher)
        submit_scan_run_btn.description = 'Stop'
        submit_scan_run_btn.icon = 'stop'
        label = _server_label(method)

        def _push_target(here, leg):
            """Where a push holds its target: this far ahead of the atoms.

            Moved up with the structure at every step, so what the molecule
            feels is a force of ``k * reach`` and not a spring that pulls
            harder the further it has to go.  Which is what makes the ramp a
            ramp of forces, and what a force means the same thing at every
            point of the path.
            """
            now = _value_in(here, leg)
            if now is None:
                return None
            way = -1.0 if float(leg['to']) < float(leg['from']) else 1.0
            return max(0.2, now + way * _gfn.PUSH_REACH)

        # What the temperature will pay for, read once before the walk starts.
        #
        # A scan is a change to the structure like any other, and with the
        # budget on it answers to the same ceiling a drag does: the walk may
        # go wherever it likes -- that is what a scan is for, and the verdict
        # below reports the whole path however high it goes -- but what is
        # handed back into the box at the end has to be somewhere the
        # temperature can reach.  Otherwise the one place the ceiling is
        # enforced can be walked round by pressing a different button.
        scan_anchor, scan_ceiling = _thermal_budget()
        budgeted = bool(submit_thermal_btn.value) and scan_anchor is not None
        state['scan_walled'] = None

        def _work():
            # The path is (where the coordinate was, what it cost, the
            # geometry) at every point.  The geometry used to be dropped the
            # moment its frame had been drawn, which was enough for a picture
            # and not enough for anything else: a finished walk is a set of
            # structures somebody has just spent minutes computing, and
            # re-pricing them with a better method is seconds of work that
            # cannot be done at all once they are gone.  It costs what the
            # structures cost -- about 60 KB for nineteen points of a 57-atom
            # complex, against the 24 MB of coordinates the frame channel
            # already holds for a walk of four hundred.
            walked, path = xyz, []
            base = None
            bottom = None
            summit = None
            began_at = None
            standing = None
            shown = []
            # The values every leg was actually held at, point by point, so
            # the walk back is the same walk and not a fresh one over the same
            # range: a scan that stopped at the next minimum walked a part of
            # what was armed, and it is that part the return leg has to
            # retrace.
            drove = []
            # And what moved that nobody asked to move, one step at a time.
            #
            # Kept as a number and a pair of indices rather than a geometry:
            # the whole reason :func:`_descent` keeps two structures and not
            # forty is that a long walk over a large molecule cannot afford
            # forty, and this must not undo that.  Which step jumped is only
            # known once the walk is over, so the culprit for every step is
            # carried along and one of them is read out at the end.
            slipped = []
            # Where the walk was last inside the budget, and what every point
            # it kept costs against the anchor.  The structure it started from
            # is affordable by construction -- it is what the box was holding.
            affordable = xyz
            costs = {}
            # And the two of them the profile marks, kept as the pair the path
            # is made of rather than as geometries: where the walk started,
            # which is the zero the free energies are quoted against, and the
            # last point the budget could pay for, which is what the box is
            # given when the walk ends past the ceiling.  Both are numbers the
            # loop has already worked out; keeping them costs nothing and
            # saves the picture having to guess which point they were.
            began_point = None
            kept_point = None
            force = _gfn.PUSH_FORCE_FROM
            growth = (_gfn.PUSH_FORCE_TO / _gfn.PUSH_FORCE_FROM) ** (
                1.0 / max(1, steps - 1))

            def _push_once(here, hard):
                """One point of a push: the force applied, and what it did."""
                asked = [
                    {'kind': one['kind'], 'atoms': list(one['atoms']),
                     'mode': 'push', 'force': _gfn.push_constant(hard),
                     'value': _push_target(here, one)}
                    for one in legs
                ]
                return asked, _gfn.optimize_with_gfn(
                    here, method, charge=charge, uhf=uhf,
                    max_steps=_SCAN_CYCLES, timeout=None,
                    constraints=asked, solvent=wet, solvation_model=model,
                    topology=_gfn_topology_dir(here))

            def _across(here, was, now_at):
                """Price the crossing the push fell through in one step.

                The push has found where the reaction goes and landed on the
                other side of it; between those two geometries the coordinate
                is now known at both ends, so driving it is no longer a guess
                about which way the reaction goes -- it is the way of asking
                what the top of the segment costs.  Held exactly, so no
                restraint energy is left in the answer.
                """
                out, standing_here = [], here
                for m in range(1, _PUSH_ACROSS):
                    if state.get('scan_stop'):
                        break
                    share_here = m / float(_PUSH_ACROSS)
                    asked = [
                        {'kind': one['kind'], 'atoms': list(one['atoms']),
                         'mode': 'fix',
                         'value': was[k] + share_here * (now_at[k] - was[k])}
                        for k, one in enumerate(legs)
                    ]
                    got = _gfn.optimize_with_gfn(
                        standing_here, method, charge=charge, uhf=uhf,
                        max_steps=_SCAN_CYCLES, timeout=None,
                        constraints=asked, solvent=wet, solvation_model=model,
                        topology=_gfn_topology_dir(standing_here))
                    if not got.get('ok') or got.get('energy') is None:
                        break
                    standing_here = got['xyz']
                    # The geometry travels with the price.  For a push these
                    # are the barrier: the force fell through the crossing in
                    # one step and these are the only points on the way over
                    # it, so a walk that threw them away had nothing to hand
                    # to a re-pricing at the one place a re-pricing is for.
                    out.append((asked[0]['value'],
                                (float(got['energy']) - base)
                                * _HARTREE_TO_KCAL,
                                standing_here))
                return out

            def _free(here):
                """The free energy of this geometry, unconstrained.

                Its own calculation, with no restraint in it: a Hessian taken
                on a biased surface has the restraint's own curvature in its
                frequencies, and the whole point of asking for G is that the
                number means something.

                An RRHO free energy only means something at a stationary
                point, and the top of a scan is not one.  The published answer
                to exactly that is xtb's ``--bhess``, the single-point Hessian
                of Spicher and Grimme (JCTC 2021, 17, 1701), which biases the
                surface back towards the geometry it was handed.  It was tried
                here and it does not apply to a scan point, which is worth
                writing down so that it is not tried again:

                  * The bias is sized in RMSD against a target of 0.10 A, and
                    a scan point can be a long way up the surface without
                    being far in RMSD.  Measured on a benzene with one ring
                    C-C held at 1.72 A -- +30.4 kcal/mol -- the free
                    relaxation moves it only 0.094 A, which is inside the
                    target, so xtb settles on ``kpush = -0.000000``, applies
                    no bias at all, relaxes freely back to the ring and
                    reports the free energy of benzene.  The held bond came
                    back at 1.385 A.  It prices a different structure and says
                    nothing about having done so, which is worse than a plain
                    Hessian, because a plain Hessian at least does not move.
                  * Asking for a tighter target does not rescue it.  With
                    ``$metadyn rmsd=0.02`` the same structure still slips from
                    1.718 to 1.523 A, the restraint it settles on is
                    thirty times stronger, the thermostatistics move the other
                    way, and it costs 65 s against 1.25 for the plain
                    Hessian.
                  * Holding the coordinate during the Hessian does keep the
                    geometry -- 0.001 A -- and puts the hold's own curvature
                    into the frequencies: the same benzene gives G(RRHO)
                    0.0963 Eh held against 0.0677 free, which is 18 kcal/mol
                    of spring in the answer.

                So it stays a plain Hessian, and what the scan does instead is
                say so: where the Hessian at the top comes back with a mode
                that goes the wrong way, that is the surface saying this point
                is not a stationary point, and the verdict reports it rather
                than quoting a free energy as though it were one.
                """
                got = _gfn.optimize_with_gfn(
                    here, method, charge=charge, uhf=uhf, timeout=None,
                    solvent=wet, solvation_model=model,
                    topology=_gfn_topology_dir(here), optimise=False,
                    free_energy=True,
                    thermo_kelvin=float(submit_temperature.value or 298.15))
                shape = got.get('imaginary') or {}
                if int(shape.get('count') or 0) > 0:
                    # Kept as the worst of the three, because one point that
                    # is not stationary is enough to make the difference
                    # between them an estimate.
                    was = state.get('scan_free_shaky') or {}
                    if int(shape.get('count')) >= int(was.get('count') or 0):
                        state['scan_free_shaky'] = dict(shape)
                return got.get('free_energy')

            def _unbiased(here, applied=()):
                """The energy of this geometry with no force on it.

                A push leaves a real restraint energy behind -- it is meant to,
                that is the force -- and xtb reports the biased total.  Priced
                as it stands, a path would carry its own push in the barrier:
                measured on an ethane pulled to 2.46 A, +1.42 kcal/mol of the
                answer was the restraint.

                Taken off rather than calculated away.  xtb's restraint is
                k*d^2 with d in Bohr for a distance and in radians for an
                angle, measured, so the residue is arithmetic -- and it agrees
                with a real calculation without the constraints to 0.0000
                kcal/mol over a distance, an angle and a torsion.  Asking for
                that calculation instead doubled the cost of every point.
                """
                got = _gfn.optimize_with_gfn(
                    here, method, charge=charge, uhf=uhf, timeout=None,
                    solvent=wet, solvation_model=model,
                    topology=_gfn_topology_dir(here), optimise=False)
                return got.get('energy')

            def _priced(got, applied):
                """What the push just reached, with its own force taken out."""
                if got.get('energy') is None:
                    return None
                bias = _gfn.restraint_energy(got['xyz'], applied, _value_in)
                if bias is None:
                    return _unbiased(got['xyz'])
                return float(got['energy']) - bias

            try:
                if pushing:
                    # The zero is the structure as it stands: the first force
                    # is applied to *it*, so the first point already has a
                    # rise to show and the stop rule has something to compare
                    # against.  Taken from the first point instead, a path
                    # that crossed on its first step showed no crossing.
                    base = _unbiased(walked)
                    if base is None:
                        schedule_ui_update(
                            _set_mol_status,
                            'The push has no starting energy to measure from.')
                        return
                    standing = walked
                    path.append((_value_in(walked, legs[0]), 0.0, walked))
                for n in range(1, steps + 1):
                    if state.get('scan_stop'):
                        break
                    share = n / float(steps)
                    # What the user is holding goes with it.
                    #
                    # A scan built its list from its own legs and nothing else,
                    # so a held value -- the whole point of "hold this while
                    # you walk that" -- never reached xtb at all.  Measured: a
                    # C-H held at 1.60 came back at 1.080 after the scan, while
                    # the same hold under Optimise gave 1.599.  The list went
                    # on showing it throughout.
                    #
                    # The geometry this step starts from, kept for one step
                    # only: it is what the next one is compared against to see
                    # whether anything slipped, and it is the previous point
                    # rather than a copy of the path.
                    stood_at = walked
                    if pushing:
                        # The next force, and the same force again halfway
                        # down when the structure fell through the crossing in
                        # one step.  Bisected from where the step *started*,
                        # so what is refined is the force and not the path.
                        came_from = walked
                        was = [_value_in(walked, one) for one in legs]
                        was_at = was[0]
                        stood = force
                        force = force * growth if n > 1 else force
                        held, outcome = _push_once(walked, force)
                        tries = 0
                        while (tries < _PUSH_REFINE and outcome.get('ok')
                               and force - stood > 0.1):
                            went = _value_in(outcome['xyz'], legs[0])
                            if went is None or abs(went - was_at) <= _PUSH_JUMP:
                                break
                            force = (stood + force) / 2.0
                            held, outcome = _push_once(walked, force)
                            tries += 1
                    else:
                        held = [
                            {'kind': one['kind'], 'atoms': list(one['atoms']),
                             'mode': 'fix',
                             'value': one['from'] + share * (one['to'] - one['from'])}
                            for one in legs
                        ]
                        walking = {tuple(one['atoms']) for one in legs}
                        held += [dict(one) for one
                                 in (state.get('constraints') or [])
                                 if tuple(one.get('atoms') or ()) not in walking]
                        outcome = _gfn.optimize_with_gfn(
                            walked, method, charge=charge, uhf=uhf,
                            max_steps=_SCAN_CYCLES, timeout=None,
                            constraints=held, solvent=wet,
                            solvation_model=model,
                            topology=_gfn_topology_dir(walked))
                    if not outcome.get('ok') or outcome.get('energy') is None:
                        schedule_ui_update(
                            _set_mol_status,
                            'The scan stopped at step '
                            f'{n}: {outcome.get("status") or "it did not run"}')
                        return
                    walked = outcome['xyz']
                    energy = (_priced(outcome, held) if pushing
                              else outcome['energy'])
                    if energy is None:
                        schedule_ui_update(
                            _set_mol_status,
                            f'The push could not be priced at step {n}.')
                        return
                    if base is None:
                        base = float(energy)
                    # Whether the method can still answer here.
                    #
                    # A closed-shell single determinant describes two electrons
                    # in one orbital, and there are regions of every surface
                    # where they are not in one: a bond half broken, a ring
                    # opening, a double bond turned.  The frontier gap says so
                    # before the energy does, and a scan is where a walk runs
                    # into such a region without anyone choosing to.
                    #
                    # Kept as the narrowest seen and compared against the
                    # first point, because what matters is where the walk went
                    # and not where it started.
                    if outcome.get('gap') is not None:
                        if state.get('scan_gap_first') is None:
                            state['scan_gap_first'] = outcome['gap']
                        narrow = state.get('scan_gap_least')
                        if narrow is None or outcome['gap'] < narrow:
                            state['scan_gap_least'] = outcome['gap']
                    spent = (float(energy) - base) * _HARTREE_TO_KCAL
                    # The same floor the drag has.  A scan drives its
                    # coordinate wherever it was told, so a target set past the
                    # far side of a bond walks straight into one -- atoms
                    # collapsing into each other, which is no path at any
                    # temperature and no number worth reporting.  Two thirds of
                    # a bond, measured between what must stay open (an SN2
                    # through its saddle came to 0.74) and what must not (a ring
                    # carbon pushed across its ring, 0.35).
                    tightest = _gfn.closest_contact(walked)[0]
                    if tightest is not None and tightest < _gfn._TOO_CLOSE:
                        # And hand back the last point that was not crowded.
                        # Stopping is not enough on its own: the geometry that
                        # tripped the floor is the one that would be left in
                        # the box, and a scan that ends on two carbons 0.98 A
                        # apart has reported a collapse rather than refused
                        # one.
                        state['scan_crowded'] = tightest
                        walked = standing if standing is not None else walked
                        break
                    standing = walked
                    # Where the coordinate *is*, not where it was asked to be.
                    # A push does not dictate a value -- that is the whole
                    # point of it -- so the path is read off the structure.
                    #
                    # Read before the two blocks below rather than after them,
                    # because both of them name a point of the path and the
                    # coordinate is half of what a point is.
                    reached = (_value_in(walked, legs[0]) if pushing
                               else held[0]['value'])
                    if budgeted:
                        # Against the budget's own anchor, not against the
                        # point the walk happened to start from: the two are
                        # the same only when the scan starts where the budget
                        # was measured, and the question here is what the
                        # structure the user is left with costs.
                        costs[walked] = ((float(energy) - float(scan_anchor))
                                         * _HARTREE_TO_KCAL)
                        if costs[walked] <= scan_ceiling:
                            affordable = walked
                            kept_point = (reached, spent)
                    if began_at is None:
                        # The first point, which is the start relaxed at the
                        # value it already had -- a minimum in every direction
                        # but the one being walked, which is what a free
                        # energy wants.
                        began_at = walked
                        began_point = (reached, spent)
                    if (pushing and was_at is not None and reached is not None
                            and abs(reached - was_at) > _PUSH_JUMP):
                        # It fell through the crossing.  What the fall cost is
                        # the barrier, and it is between the two geometries
                        # rather than at either of them.
                        path.extend(_across(
                            came_from, was,
                            [_value_in(walked, one) for one in legs]))
                    path.append((reached, spent, walked))
                    if not pushing:
                        # The values this step really held, for the walk back,
                        # and the largest thing that moved without being
                        # asked, for naming what slipped if anything did.
                        drove.append([one['value'] for one in held[:len(legs)]])
                        slipped.append(_gfn.what_else_moved(
                            stood_at, walked,
                            [one['atoms'] for one in legs]))
                    # The lowest point *since the top*, kept with its
                    # geometry.
                    #
                    # The climb out is what says the minimum is behind us, so
                    # by the time it is recognised the walk stands two steps
                    # past it and the structure in the box is squeezed that
                    # much.  Stopping is not enough; it has to come back to
                    # the bottom it walked through, which is the geometry
                    # anyone would want to carry on from.
                    #
                    # Since the top, and not the lowest of the whole path:
                    # see :func:`_descent`.
                    summit, bottom = _descent(
                        summit, bottom, spent, walked, reached)
                    # The instruction, before the energy's own stop rule.
                    #
                    # They answer different questions and the instruction is
                    # the one that was asked: "over a barrier and settled
                    # again" is where a scan stops when nobody said what the
                    # reaction was, and here somebody did.  Measured on the
                    # Diels-Alder, both fire at the same step; measured on the
                    # SN2 -- Cl- and CH3Br, form Cl-C while breaking C-Br --
                    # the whole crossing is downhill from the complex and
                    # there is no barrier for the energy rule to notice at
                    # all, so the instruction is the only thing that could
                    # have stopped it.
                    #
                    # And it is not gated on Whole profile, which is a
                    # question about how much of a curve to draw.  An
                    # instruction that has been carried out is finished.
                    if instructed and _carried_out(walked, legs):
                        state['scan_carried_out'] = (
                            n, force,
                            _gfn.graph_changed(
                                _gfn.bond_graph(xyz), _gfn.bond_graph(walked),
                                [line.split()[0]
                                 for line in _gfn.atom_lines(walked)]))
                        break
                    if not submit_scan_whole.value and _scan_arrived(path):
                        state['scan_arrived'] = True
                        if bottom is not None:
                            walked = bottom[1]
                            # The path itself is *not* cut back to it.  Cut,
                            # the verdict describes the stump rather than the
                            # walk: the barrier that was just crossed is
                            # thrown away and what is reported is the one
                            # point that is left.
                            state['scan_came_back'] = (bottom[2], bottom[0])
                        break
                    schedule_ui_update(
                        _set_mol_status,
                        f'{label} is '
                        + (f'pushing at {force:.0f} kcal/mol/A: '
                           if pushing else 'walking the scan: ')
                        + f'step {n} of {steps}, '
                        f'{held[0]["kind"]} at {reached:.3g}, '
                        f'{spent:+.1f} kcal/mol so far.', spinner=True)
                    # Down the frame channel, not into the box.
                    #
                    # Every write to the box rebuilds the viewer from nothing,
                    # so a thirty-point scan was thirty rebuilds -- the picture
                    # crawled and the browser sometimes stopped answering
                    # altogether.  The follow has always streamed frames for
                    # exactly this reason; the scan was writing the box because
                    # it was easier, and thirty rebuilds is not a detail.
                    shown.append(_gfn.coordinates_of(walked))
                    scan_run = state.get('scan_frame_run')
                    schedule_ui_update(
                        lambda text=_frame_payload(
                            scan_run, **{'from': len(shown) - 1,
                                         'follow': 1,
                                         'frames': [shown[-1]]}),
                        r=scan_run: setattr(submit_gfn_frame, 'value', text)
                        if _frame_run_is_current(r) else None)
                # The walk as a pair of numbers per point, which is what
                # :func:`_scan_two_legs` promises and what the second leg is
                # compared against.  The geometry each point held travels with
                # the walk itself, not here: :func:`_keep_the_walk` is what
                # keeps the structures, for re-pricing and for marking a point,
                # and a leg is two axes.
                state['scan_there'] = [(one[0], one[1]) for one in path]
                # Whether it was the walk *out* that was interrupted.  A Stop
                # pressed during the return leg is a different thing and must
                # not turn the barrier into "where the walk was interrupted":
                # the walk out finished, and what it found stands.
                state['scan_stopped_out'] = bool(state.get('scan_stop'))
                # Whether any step of it was a fall rather than a step.
                #
                # Costs nothing: it is arithmetic on the energies the walk
                # already has and on one number per step that was worked out
                # while the step was taken.  Walking only -- a push means to
                # fall through its crossing and prices it afterwards with
                # :func:`_across`, so the same test on a push would fire on
                # the thing the push is for, and on a path whose points are
                # not evenly spaced besides.
                if not pushing:
                    fell = _gfn.where_a_walk_jumped(
                        [one[1] for one in path])
                    if fell is not None:
                        step = int(fell['step'])
                        fell['at'] = path[step][0]
                        fell['from'] = path[step - 1][0]
                        who = slipped[step] if step < len(slipped) else None
                        if who is not None:
                            rows = [line.split()[0] for line
                                    in _gfn.atom_lines(walked)]
                            fell['named'] = _gfn.pair_named(who['pair'], rows)
                            fell['moved'] = who['moved']
                            fell['was'] = who['was']
                            fell['now'] = who['now']
                        state['scan_jumped'] = fell
                # And the same coordinate walked back from where it ended.
                #
                # A driven scan is a minimum-energy path only where nothing
                # slips sideways, and there is no way to know which case a
                # given scan is in from one leg of it.  So the second leg is
                # walked over the values the first one really held -- which is
                # not the range that was armed, because a scan that stopped at
                # the next minimum walked a part of it.
                #
                # From where the walk ended and not from the minimum it came
                # back to: the question is whether retracing the same points
                # gives the same energies, and it has to start at the far end
                # of them.
                #
                # Not after a Stop, which is a user saying they have seen
                # enough, and not after a collapse, where the forward leg has
                # already reported that there is no path there.
                if (not pushing and bool(submit_scan_back.value)
                        and len(drove) > 2 and not state.get('scan_stop')
                        and state.get('scan_crowded') is None
                        and base is not None):
                    here = standing if standing is not None else walked
                    returned = [(drove[-1][0], path[-1][1])]
                    walking = {tuple(one['atoms']) for one in legs}
                    others = [dict(one) for one
                              in (state.get('constraints') or [])
                              if tuple(one.get('atoms') or ()) not in walking]
                    for back_n, values in enumerate(
                            reversed(drove[:-1]), start=1):
                        if state.get('scan_stop'):
                            break
                        asked = [
                            {'kind': one['kind'],
                             'atoms': list(one['atoms']), 'mode': 'fix',
                             'value': values[k]}
                            for k, one in enumerate(legs)
                        ] + others
                        got = _gfn.optimize_with_gfn(
                            here, method, charge=charge, uhf=uhf,
                            max_steps=_SCAN_CYCLES, timeout=None,
                            constraints=asked, solvent=wet,
                            solvation_model=model,
                            topology=_gfn_topology_dir(here))
                        if not got.get('ok') or got.get('energy') is None:
                            break
                        here = got['xyz']
                        returned.append(
                            (values[0], (float(got['energy']) - base)
                             * _HARTREE_TO_KCAL))
                        schedule_ui_update(
                            _set_mol_status,
                            f'{label} is walking it back: step {back_n} of '
                            f'{len(drove) - 1}, {legs[0]["kind"]} at '
                            f'{values[0]:.3g}, {returned[-1][1]:+.1f} '
                            'kcal/mol.', spinner=True)
                        # Shown as it happens, like the leg before it.  A
                        # return leg nobody can see is a number with no
                        # picture behind it, and the picture is what makes a
                        # jump obvious.
                        shown.append(_gfn.coordinates_of(here))
                        scan_run = state.get('scan_frame_run')
                        schedule_ui_update(
                            lambda text=_frame_payload(
                                scan_run, **{'from': len(shown) - 1,
                                             'follow': 1,
                                             'frames': [shown[-1]]}),
                            r=scan_run: setattr(submit_gfn_frame, 'value',
                                                text)
                            if _frame_run_is_current(r) else None)
                    state['scan_back'] = list(returned)
                    # Two points is not two legs.  A return leg that failed at
                    # its first step leaves only the geometry it started
                    # standing on, which agrees with itself by construction --
                    # and reporting that as agreement would be the check
                    # saying yes because it did not run.
                    if len(returned) > 2:
                        # The pairs, not the walk: both legs are (coordinate,
                        # energy) on one zero, and the geometries the walk out
                        # carries are no part of comparing two curves.
                        state['scan_disagree'] = _gfn.paths_disagree(
                            state['scan_there'], returned)
                # And the free energy, at the three places it is both
                # affordable and meaningful: where the walk started, the
                # highest point it crossed, and the minimum it came to.
                if (str(submit_scan_energy.value) == 'G' and began_at
                        and summit is not None and not state.get('scan_stop')):
                    schedule_ui_update(
                        _set_mol_status,
                        f'{label} is taking three Hessians for the free '
                        'energies...', spinner=True)
                    ends = bottom[1] if bottom is not None else walked
                    one, top, other = (_free(began_at), _free(summit[1]),
                                       _free(ends))
                    state['scan_free'] = (
                        None if None in (one, top, other) else
                        ((top - one) * _HARTREE_TO_KCAL,
                         (other - one) * _HARTREE_TO_KCAL))
                # Whether the method could describe what it walked through,
                # asked twice: from the frontier gap, which every point of the
                # walk reports for nothing, and by counting the electrons that
                # are not in a closed shell, which costs two more single
                # points.
                #
                # The count is skipped on a Stop, for the same reason the free
                # energies are: a press of Stop is somebody saying that is far
                # enough, and answering it by starting two more calculations
                # is the switch not working.  The gap costs nothing and is
                # said either way.
                #
                # Joined rather than concatenated: either half can be empty,
                # and the verdict puts a single space before whatever this is.
                state['scan_depth'] = ' '.join(part for part in (
                    _gfn.method_is_out_of_its_depth(
                        state.get('scan_gap_least'),
                        state.get('scan_gap_first')),
                    '' if state.get('scan_stop') else _fod_along_the_walk(
                        began_at, summit, method, charge, uhf, wet, model),
                ) if part)
                # The two ends, for whoever wants to walk between them.
                # A path finder is given two structures and finds its own way;
                # what it cannot do is invent a product to aim at, and this is
                # where one comes from.
                if began_at and path:
                    state['scan_ends'] = (
                        began_at, bottom[1] if bottom is not None else walked)
            finally:
                state['scan_run'] = False
                # The walk is reported whole; what is handed back is not.
                #
                # A path that climbs past the ceiling is exactly the answer a
                # scan is asked for, and the verdict says how high it went and
                # what temperature would cross it.  Leaving that geometry in
                # the box is a different thing: it is the structure the user
                # carries on from, and at this temperature they cannot get to
                # it.  So the box gets the last point the budget could pay for
                # and the line says which one that was.
                if budgeted and costs.get(walked, 0.0) > scan_ceiling:
                    state['scan_walled'] = costs[walked]
                    walked = affordable
                # And whether the instruction survives that.  It is the one
                # place the budget and a form/break instruction can disagree,
                # and they disagree by design: the ceiling prices exactly this
                # kind of deformation, and a bond driven apart by a force the
                # temperature cannot supply is a reaction that does not happen
                # at that temperature.  The verdict below says which of the
                # two the user is looking at, because the alternative is a
                # line reading "the bond broke" over a structure where it has
                # not.
                state['scan_carried_out_kept'] = bool(
                    instructed and state.get('scan_carried_out')
                    and _carried_out(walked, legs))

                def _done(final=walked):
                    submit_scan_run_btn.description = 'Run scan'
                    submit_scan_run_btn.icon = 'play'
                    # A run that walked nothing writes nothing and says only
                    # why.  It used to rewrite the box's comment to 'Scanned'
                    # over a geometry that had never moved, and then report
                    # 'The scan walked nothing' as though that were a result.
                    if not path:
                        return
                    # There are two ends now, and the press beside them is put
                    # on the thing this walk made possible -- see
                    # :func:`_scan_left_two_ends`, which is where the moment
                    # is answered rather than merely allowed for.
                    if state.get('scan_ends'):
                        _scan_left_two_ends()
                    # And the walk, kept whole with the charge, the spin and
                    # the solvent it ran under.  Kept here rather than where
                    # it finished, because a walk that was stopped or that
                    # failed part way is still a set of geometries somebody
                    # paid for: the press below prices whatever is in hand and
                    # says how much of the walk that was.
                    _keep_the_walk(path, method, charge, uhf, wet, model)
                    # And the press arrives with it, which is the row saying
                    # the profile can now be checked against something better
                    # -- the same news the closing sentence carries in words.
                    _refresh_scan()
                    # The places along the walk worth coming back to, put into
                    # the history in the order they were reached, so that Undo
                    # steps back through them and Redo forward again.
                    #
                    # Undo is where the user asked for this, and Undo is also
                    # where it belongs: a scan is one action and one entry
                    # already, and the three landmarks are that one action's
                    # own stations rather than three separate things that were
                    # done.  Nothing new to press, which is the point -- these
                    # are geometries the scan had already worked out for the
                    # three Hessians the free-energy mode takes, and they were
                    # being thrown away.
                    #
                    # The entry before them is the structure from before the
                    # scan, recorded when it started, so the walk back reads:
                    # where it came to, the highest point it crossed, where it
                    # started, and then the structure that was there before
                    # any of it.
                    if began_at:
                        _remember_landmark(
                            began_at, 'the scan, on from where it started',
                            'Scanned: where the walk started')
                    if summit is not None and summit[1]:
                        _remember_landmark(
                            summit[1],
                            'the scan, on from the highest point it crossed',
                            'Scanned: the highest point the walk crossed')
                    rows = [line for line in final.splitlines()[2:]
                            if line.strip()]
                    if rows:
                        _write_coords(xyz_document(
                            rows,
                            'Scanned, and back to the last point the '
                            'temperature can pay for'
                            if state.get('scan_walled') is not None else
                            'Driven until the bonds were made and broken'
                            if state.get('scan_carried_out') else
                            'Scanned to the next minimum'
                            if state.get('scan_arrived') else 'Scanned'),
                            run=state.get('scan_frame_run'))
                    _set_mol_status(*_scan_verdict(path, steps))

                schedule_ui_update(_done)
                # And the picture of the walk, after all of that: drawn here,
                # on the thread that did the walking, and handed to the
                # interface in a turn of its own.
                #
                # After rather than before, and measured: matplotlib's first
                # import is 0.3 s idle and over a second on a loaded machine,
                # and drawn before this line that whole delay fell between
                # the walk ending and the verdict reaching the screen -- six
                # of the scan tests caught it, finding the row still saying
                # "step 14 of 20" after the run had finished.  The answer
                # goes first; the picture follows it.
                #
                # Shown rather than drawn in the turn for the same reason:
                # a second spent in an interface callback is a second in
                # which the dashboard does not answer.
                schedule_ui_update(
                    _show_scan_profile,
                    _scan_profile_html(path, legs, pushing,
                                       began=began_point, kept=kept_point))

        _start_background(_work, 'The scan')

    def _fod_along_the_walk(began_at, summit, method, charge, uhf, wet, model):
        """How much of the walk was not a closed shell, said as a change.

        Two more single points on geometries the walk already produced: where
        it started, and the highest point it crossed.  That pair and not
        another, because the barrier is the number a scan exists to produce
        and the top is where it comes from -- so this is measured at the place
        the answer is quoted from, against the place it is quoted against.

        Beside the frontier gap rather than instead of it.  The gap is free --
        every point of the walk prints one -- and it is a proxy; this counts
        the electrons that are not paired and so has a scale.  Measured under
        GFN2 across an ethylene turned about its double bond, the two agree
        about where the trouble is: at 60 degrees, which is where
        :func:`~delfin.dashboard.gfn_optimize.method_is_out_of_its_depth`
        already fires, N_FOD has just left zero at 0.251, and by 90 degrees it
        is 2.000 with the gap shut.

        What it costs is two runs on a walk of twenty optimisations, and the
        ratio is what matters because the machine is shared.  Minimum of seven,
        against a plain single point and against one optimisation step of the
        same structure: 0.96 s, 0.49 s and 7.4 s for the sixteen-atom
        Diels-Alder complex, and 10.5 s, 3.5 s and 35.6 s for the 57-atom
        manganese complex.  So the pair of them is well under one step of a
        walk that takes twenty, which is why it is not a thing to switch on.

        Asked of GFN2 whatever walked the scan, which is the one decision in
        here that needs defending.  It is a probe rather than a result: the
        published one is TPSS/def2-TZVP used on everybody's geometries, and
        what it measures is a property of the structure -- how far from a
        closed shell it is -- not of the method that arrived at it.  Two of
        the four methods on this toolbar cannot be asked at all: GFN-FF has no
        electrons, and g-xTB takes ``--fod``, converges, terminates normally
        and prints nothing, which is silence that reads as good news.

        Under g-xTB that is not a corner case, it is the main one.  The
        frontier-gap rule this stands beside was calibrated on GFN2 and does
        not transfer: measured on an ethylene turned about its double bond,
        where 90 degrees is a textbook diradical that no closed shell
        describes,

            twist       0     30     60     90
            GFN2     5.47   4.37   2.32   0.00 eV
            g-xTB   12.17  10.55   7.94   5.12 eV

        so a g-xTB profile never trips
        :data:`~delfin.dashboard.gfn_optimize.GAP_IS_SMALL` and the warning is
        deaf under the most accurate method in the list -- which is the one a
        barrier is most likely to be quoted from.  Handing the two geometries
        to GFN2 gives every walk the same check for two single points, and the
        sentence names who was asked when it was not the method that walked.
        """
        if not began_at or summit is None or not summit[1]:
            return ''
        probe = (method if _gfn.can_measure_fod(method)
                 else _gfn.FOD_METHODS[-1])
        borrowed = str(probe).lower() != str(method).lower()
        schedule_ui_update(
            _set_mol_status,
            f'{_server_label(probe)} is measuring how much of the top is not '
            'a closed shell...', spinner=True)

        def _measure(here):
            got = _gfn.optimize_with_gfn(
                here, probe, charge=charge, uhf=uhf, timeout=None,
                # The solvent only where the probe is the method that walked.
                # GFN2 borrowed for a g-xTB walk has no solvent to inherit --
                # that build takes none -- and putting one on would make the
                # two points a different calculation from the walk they are
                # about.
                solvent=None if borrowed else wet,
                solvation_model=model, optimise=False, fod=True,
                # No topology, so no warm restart: the extra SCF converges at
                # 5000 K and its wavefunction is not a guess at the ordinary
                # one -- see :func:`optimize_with_gfn`, which refuses to store
                # it for that reason.
                topology=None)
            return got.get('fod') if got.get('ok') else None

        first, top = _measure(began_at), _measure(summit[1])
        if not first or not top:
            return ''
        said = _gfn.fod_moved(first['total'], top['total'], top['on'])
        if not said:
            return ''
        whose = (f' ({_server_label(probe)} was asked, on the geometries '
                 f'{_server_label(method)} walked.)' if borrowed else '')
        return f'{said}{whose}'

    def _keep_the_walk(path, method, charge, uhf, wet, model):
        """Put the finished walk where a second opinion can be taken on it.

        With the conditions it was run under, because a re-pricing that used
        a different charge or a different spin would not be a second opinion
        about the same thing -- and with the element column of what it walked
        on, so the press can be offered only while it is about the molecule on
        screen.  The two ends already answer to that rule
        (:func:`_scan_ends_here`) for the same reason: a set of geometries
        outlives the structure it was made from, and geometries belonging to
        another molecule are not this molecule's profile.
        """
        points = [one for one in path if len(one) > 2 and one[2]]
        if not points:
            state['scan_walk'] = None
            return
        state['scan_walk'] = {
            'points': [(one[0], one[1], one[2]) for one in points],
            'method': method, 'charge': int(charge), 'uhf': int(uhf),
            'solvent': wet or '', 'solvation_model': model,
            'structure': _structure_fingerprint(points[0][2]),
        }

    def _walk_here():
        """The last walk, while it is about the molecule on screen.

        Absent for that reason and no other, which is the rule every control
        in this editor that comes and goes is built on.
        """
        walk = state.get('scan_walk')
        if not walk or not walk.get('points'):
            return None
        if walk['structure'] != _structure_fingerprint(_current_xyz() or ''):
            return None
        return walk

    #: What a finished profile is re-priced with.
    #:
    #: One method and not a box of them, because there is one answer: g-xTB is
    #: the only thing in this editor's list that is better than what walked
    #: the scan, and offering a choice between GFN1 and GFN2 to somebody who
    #: has just walked with GFN2 is a box whose entries are all wrong.
    _REPRICE_WITH = 'gxtb'

    def _reprice_is_possible():
        """Whether there is a finished walk that something better can price."""
        walk = _walk_here()
        if not walk:
            return False
        # Nothing to improve on.  A walk made with the best method in the list
        # has no second opinion available here, and a press offering one would
        # be the toolbar claiming something it cannot do.
        return str(walk['method']).lower() != _REPRICE_WITH

    def on_submit_scan_price(_button=None):
        """Price the walk again with a better method, at its own geometries.

        Why it is worth a press.  GFN2 has no Fock exchange at all, which puts
        it at the far end of the self-interaction axis, and a transition state
        is where that hurts most: stretched bonds and near-degenerate weakly
        bound orbitals, which Bursch, Mewes, Hansen and Grimme (Angew. Chem.
        Int. Ed. 2022, https://doi.org/10.1002/anie.202205735) say "often
        leads to a systematic underestimation of their electronic energy and,
        in turn, barrier heights".  The error is structured rather than
        random: published mean absolute deviations of 10.22 kcal/mol on
        pericyclic barriers, 8.12 on BHDIV10 and 13.05 on BH9 aggregate,
        against 1.17 on BHROT27 -- rotations, where no bond is broken.  So a
        torsion profile is nearly right and a bond-breaking one is not.

        Measured here on a butadiene and an ethylene walked together under
        GFN2 in nineteen points, then re-priced: the barrier went from
        +6.31 kcal/mol to +21.67, a shift of +15.4 against a published
        pericyclic deviation of 10.22, and the reaction energy from -63.99 to
        -42.71.  One reaction is not a benchmark; the direction and the size
        are what the literature says to expect.

        Why it is a press on a finished scan rather than a setting before one.
        The geometries are already in hand, so this is seconds -- measured,
        0.40 s a point, 7.5 s for the whole nineteen -- while walking the scan
        again under another method is minutes.  And the question only arises
        after the barrier is on screen and looks too small to believe.

        It gets cheaper as the molecule gets bigger, which is the opposite of
        what a per-point cost usually does here.  Driven whole on a 57-atom
        manganese complex at charge +1, its Mn-Br bond walked out over six
        points on a loaded machine: the walk took 1097 s and pricing the same
        six geometries again took 35.8 s -- 5.96 s a point against 183 s a
        point, three percent of what it is checking.  A constrained
        optimisation grows with the molecule far faster than one single point
        does.
        """
        if state.get('reprice_run'):
            state['reprice_stop'] = True
            _set_mol_status('Stopping after this point.')
            return
        walk = _walk_here()
        if not walk:
            return
        points = walk['points']
        label = _server_label(_REPRICE_WITH)
        state['reprice_run'] = True
        state['reprice_stop'] = False
        submit_scan_price_btn.description = 'Stop'
        submit_scan_price_btn.icon = 'stop'
        # g-xTB takes no implicit solvation in this build, so a walk made in a
        # solvent is re-priced in the gas phase.  That is a different
        # experiment and it is said in the answer rather than refused: the
        # thing being checked is a barrier that may be ten kcal/mol too low,
        # and a solvation shift on a barrier is not that.
        dry = bool(walk.get('solvent'))

        def _work():
            priced, failed = [], ''
            try:
                for n, (value, _spent, here) in enumerate(points, start=1):
                    if state.get('reprice_stop'):
                        break
                    got = _gfn.optimize_with_gfn(
                        here, _REPRICE_WITH, charge=walk['charge'],
                        uhf=walk['uhf'], timeout=None, optimise=False)
                    if not got.get('ok') or got.get('energy') is None:
                        failed = got.get('status') or 'it did not run'
                        break
                    priced.append((value, float(got['energy'])))
                    schedule_ui_update(
                        _set_mol_status,
                        f'{label} is pricing the walk again: point {n} of '
                        f'{len(points)}.', spinner=True)
            finally:
                state['reprice_run'] = False

                def _done():
                    submit_scan_price_btn.description = (
                        f'Price with {label}')
                    submit_scan_price_btn.icon = 'balance-scale'
                    if not priced:
                        _set_mol_status(
                            f'{label} priced nothing: {failed}'
                            if failed else
                            f'{label} priced nothing before it was stopped.')
                        return
                    base = priced[0][1]
                    profile = [(value, (energy - base) * _HARTREE_TO_KCAL)
                               for value, energy in priced]
                    # Kept for whoever draws the profile.  Both curves are
                    # about the same geometries and are indexed by the same
                    # coordinate, so they go on one pair of axes against
                    # ``state['scan_walk']['points']``, which carries the
                    # first one as (coordinate, kcal/mol, geometry).
                    #
                    # The element column travels with it for the reason the
                    # press answers to it: a profile outlives the structure it
                    # was walked on, and a second curve drawn over a molecule
                    # it is not about is worse than no second curve.  Whoever
                    # draws this checks it the way :func:`_walk_here` does.
                    state['scan_repriced'] = {
                        'method': _REPRICE_WITH, 'label': label,
                        'points': profile,
                        'walked_with': _server_label(walk['method']),
                        'of': len(points),
                        'gas_phase': dry,
                        'structure': walk['structure'],
                    }
                    _set_mol_status(*_repriced_verdict(
                        walk, profile, label, dry, failed))

                schedule_ui_update(_done)

        _start_background(_work, 'The second opinion',
                          guards={'reprice_run': False})

    def _repriced_verdict(walk, profile, label, dry, failed):
        """Both profiles, and the three things that stop this being a barrier.

        The caveats are here and not in a docstring because they are the
        answer: a number that moved by fifteen kcal/mol is going to be quoted
        by somebody, and what it is worth has to travel with it.  Said in
        three sentences, in words that assume no chemistry.
        """
        walked = _server_label(walk['method'])
        first = [(one[0], one[1]) for one in walk['points']][:len(profile)]
        was_top = max(first, key=lambda one: one[1])
        now_top = max(profile, key=lambda one: one[1])
        moved = now_top[1] - was_top[1]
        many = ('' if len(profile) == len(walk['points']) else
                f' Only {len(profile)} of {len(walk["points"])} points were '
                f'priced'
                + (f': {failed}' if failed else ', because it was stopped.'))
        # The top can move as well as rise, and where it moved to is a fact
        # about the profile rather than a detail: measured on the
        # Diels-Alder, GFN2 put its highest point at 2.289 A and g-xTB put it
        # at 2.208 on the same nineteen geometries.
        elsewhere = ('' if abs(now_top[0] - was_top[0]) < 1e-9 else
                     f' -- and at {now_top[0]:.3g} rather than '
                     f'{was_top[0]:.3g}')
        return (
            f'{walked} walked it and {label} priced it again at the same '
            f'{len(profile)} geometries. The highest point goes from '
            f'{was_top[1]:+.1f} to {now_top[1]:+.1f} kcal/mol, a change of '
            f'{moved:+.1f}{elsewhere}, and the end from {first[-1][1]:+.1f} '
            f'to {profile[-1][1]:+.1f}.{many}',
            'Three things this is not. The energies are better and the '
            f'structures are not -- they are where {walked} put them, and it '
            'is out by about 0.2 A on a half-broken bond against 0.03 A on an '
            'ordinary one, which is the kind of bond a barrier is made of. So '
            'it is a screen: to report a barrier, optimise the top again with '
            f'the method you mean to report. {label} is a preprint method and '
            'this build says it is a development version differing from the '
            'paper. And the points are wherever the walk left them, which '
            'need not be anywhere in particular -- pricing changes the '
            'energies, never the geometries.'
            + (f' The walk ran in {walk["solvent"]} and {label} has no solvent '
               'in this build, so these are gas-phase energies at those '
               'geometries.' if dry else ''))

    def _said_modes(shape, what, advise=True):
        """What a Hessian says a structure is, in sentences.

        The wording and the thresholds it is said against live in
        :func:`delfin.dashboard.saddle.verdict`, with the searches each
        threshold was taken from cited where it is defined.  One place,
        because the path finder's estimate and the saddle search's answer are
        the same question asked twice, and a structure that is named a
        second-order saddle by one of them must not be called a transition
        state by the other.
        """
        return _saddle.verdict(shape, what, advise=advise)['lines']

    def _said_thermochemistry(thermo):
        """G, H, T*S and the zero-point energy, in one line, or nothing.

        In kcal/mol, which is what the rest of this editor prices things in,
        and the entropy as T*S rather than as S so that all four numbers are
        the same quantity and can be added and subtracted by eye.  The
        temperature is named because every one of them depends on it and the
        box that set it is at the other end of a crowded row.

        The zero-point energy is said apart from the rest: it is the only one
        of the four that has nothing to do with the temperature, and a
        structure with a soft mode has a small one for a reason worth
        noticing.
        """
        if not thermo:
            return []
        warmth = float(thermo.get('kelvin') or 0.0)
        free = thermo.get('free_energy')
        heat = thermo.get('enthalpy')
        ts = thermo.get('ts')
        zpe = thermo.get('zpe')
        if free is None and heat is None:
            return []
        parts = []
        if free is not None:
            parts.append(f'G = {float(free) * _HARTREE_TO_KCAL:.2f}')
        if heat is not None:
            parts.append(f'H = {float(heat) * _HARTREE_TO_KCAL:.2f}')
        if ts is not None:
            parts.append(f'T*S = {float(ts) * _HARTREE_TO_KCAL:.2f}')
        said = [f'At {warmth:g} K, in kcal/mol: ' + ', '.join(parts) + '.']
        if zpe is not None:
            said.append('The zero-point energy is '
                        f'{float(zpe) * _HARTREE_TO_KCAL:.2f} kcal/mol, which '
                        'no temperature takes away.')
        return said

    def _said_curvature(shape):
        """What the modes say where nothing is standing still.

        The naming in :func:`delfin.dashboard.saddle.verdict` is about
        stationary points -- "a minimum", "a transition state" -- and those
        words are not available here: a Hessian on a slope is the curvature of
        a hillside, and calling a hillside with no downward mode "a minimum"
        is exactly the false confidence this press exists to prevent.
        Measured on a hand-built methane with every C-H at 1.0897 A, which is
        nobody's stationary point: xtb counts zero imaginary modes there and
        is perfectly right to, and read through the naming it came out as "is
        a minimum".

        So on a slope the count is said as a count, and no structure is named.
        """
        if not shape:
            return ['Its modes could not be checked.']
        modes = [float(one) for one in (shape.get('modes') or [])]
        order = shape.get('count')
        order = len(modes) if order is None else int(order)
        if order == 0:
            return ['Every mode there curves upwards, which is what the '
                    'surface looks like on the way down into a well -- but '
                    'which well, and how far, is what Optimise answers.']
        listed = ', '.join(f'{one:.0f}' for one in modes[:max(0, order)])
        return [f'{order} mode(s) curve downwards there'
                + (f' ({listed} cm-1)' if listed else '')
                + '. On a slope that is the shape of the hillside and not a '
                  'saddle: a transition state is a stationary point, and this '
                  'is not one.']

    def _lowest_bond_orders_said(bonds, xyz, gap=None, most=2):
        """The lowest bond orders in this answer, named -- as a readout.

        Only the pairs xtb printed an order for, and only the ones under
        :data:`~delfin.dashboard.gfn_optimize.BOND_WORTH_SAYING`: a structure
        whose lowest order is a proper single bond has nothing to say here,
        and printing "the lowest is 0.98" on every press would train the
        reader to skip the line that matters.

        It says nothing about whether those bonds exist, and it is worth being
        explicit about why, because the opposite was built here first: a
        Wiberg order is a readout of the wavefunction, and under a closed-shell
        Hamiltonian a homolytically broken bond still reads about one -- an
        ethane at 3.03 A reads 1.000. The bond watch stays geometric. What
        this is for is the ordinary reading a chemist does: a long dative bond
        at 0.24, a partial bond in a structure being built.
        """
        if not bonds:
            return []
        rows = [line.split() for line in _gfn.atom_lines(xyz or '')]
        if not rows:
            return []
        lowest = sorted(
            (one for one in bonds if float(one[2]) < _gfn.BOND_WORTH_SAYING),
            key=lambda one: float(one[2]))[:int(most)]
        said = []
        for first, second, order in lowest:
            if not (0 <= first < len(rows) and 0 <= second < len(rows)):
                continue
            names = (f'{rows[first][0]}{first}-{rows[second][0]}{second}')
            said.append(_gfn.bond_order_note(order, names, gap))
        said = [one for one in said if one]
        if said:
            # Once, above them, and only where there are any to read. This is
            # the sentence that keeps the numbers from being read as a bond
            # watch, and it is the press that can afford it -- the drag line
            # quotes the bare number, several times a second.
            said.insert(0, 'The lowest bond orders here, as a readout -- an '
                            'order is not a test of whether a bond is there, '
                            'and under a closed shell a bond broken '
                            'homolytically still reads about one:')
        return said

    def on_submit_shape(_button=None):
        """One press: what is the structure on screen, and what does it cost.

        A Hessian on the geometry exactly as it stands -- nothing is
        optimised, nothing is moved, and the coordinate box is not written to.
        That is the point: the question is about the structure the user built,
        and a press that quietly relaxed it first would answer about a
        different one.

        What comes back is said in three parts, and the order matters.
        Whether it is standing still comes first, because everything after it
        is read differently if it is not: a Hessian at a geometry that is not
        a stationary point describes a point on a slope, and the honest
        sentence for a structure with a large gradient is that it is neither a
        minimum nor a saddle. Then what the modes say, in the same words the
        saddle search uses. Then what it costs at the temperature on the row.
        """
        if state.get('shape_run'):
            return
        xyz = _current_xyz()
        method = str(submit_ff_dd.value)
        if not xyz:
            _set_mol_status('There is no structure to ask about.')
            return
        if not _gfn.is_gfn_method(method):
            # It should not be reachable -- the press is not on the row under
            # anything else -- but a method can change while a press is in
            # flight, and a button that acts on the wrong engine is worse than
            # one that says so.
            _set_mol_status(
                'A Hessian here is xtb\'s, so choose GFN-FF, GFN1, GFN2 or '
                'g-xTB. Anything with a basis set is a job for the ORCA '
                'Builder.')
            return
        state['shape_run'] = True
        label = _server_label(method)
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        model = _solv_model()
        warmth = float(submit_temperature.value or 298.15)
        atoms = len(_gfn.atom_lines(xyz))
        _set_mol_status(
            f'{label} is taking a Hessian on the {atoms} atoms as they '
            'stand...', spinner=True)

        def _work():
            # No cap on the clock. A Hessian grows steeply with the structure
            # -- measured, 0.4 s for five atoms, 1.1 s for 23 and 23.7 s for
            # 57 -- and a press that gave up at an arbitrary second would
            # throw away exactly the answers that took long enough to be
            # worth asking for. The ring says it is running, and the press is
            # the user's to make.
            found = _gfn.optimize_with_gfn(
                xyz, method, charge=charge, uhf=uhf, timeout=None,
                solvent=wet, solvation_model=model,
                optimise=False, free_energy=True, thermo_kelvin=warmth)

            def _done():
                state['shape_run'] = False
                if not found.get('ok'):
                    _set_mol_status(str(found.get('status')
                                        or 'The Hessian did not run.'))
                    return
                _remember_charges(found)
                _repaint_labels()
                lines = [f'{label}: {found["seconds"]:.1f} s for the Hessian '
                         f'on {atoms} atoms, and the structure is untouched.']
                # First, because it decides how the rest may be worded: the
                # names a Hessian goes by -- a minimum, a transition state --
                # all mean "stationary point of this order", and none of them
                # is available for a structure that is still on its way
                # somewhere.
                slope = _gfn.not_a_stationary_point(found.get('gradient'),
                                                    atoms)
                if slope:
                    lines.append(slope)
                    lines.extend(_said_curvature(found.get('imaginary')))
                else:
                    lines.extend(_said_modes(
                        found.get('imaginary'), 'The structure as it stands'))
                lines.extend(_said_thermochemistry(found.get('thermo')))
                # And whether the method is still able to answer at all.
                # Before the bonds, because it is what decides whether their
                # orders may be read as evidence at all.
                depth = _gfn.method_is_out_of_its_depth(found.get('gap'))
                if depth:
                    lines.append(depth)
                lines.extend(_lowest_bond_orders_said(
                    found.get('bonds'), xyz, found.get('gap')))
                _set_mol_status(*lines)

            schedule_ui_update(_done)

        _start_background(_work, 'The Hessian',
                          guards={'shape_run': False})

    def on_submit_saddle(_button=None):
        """The one press for a transition state, and what the boxes say.

        There were three buttons and they were three of the six answers to two
        questions -- where the search starts, and how it gets there -- with
        nothing on screen saying so.  The questions are asked in the two boxes
        beside this, so all six are reachable and the press is one press.

        Stopping is the same press, whichever half is running, because while
        something is walking or climbing there is only one thing to want.  A
        minimisation that is following the wrong thing fails and says so; a
        climb does not -- it succeeds at arriving somewhere nobody wanted --
        and the only way to know that early is to watch it happen and be able
        to interrupt it.
        """
        if state.get('saddle_run'):
            # What it reached is kept.
            state['saddle_stop'] = True
            _set_mol_status('Stopping the climb...', spinner=True)
            return
        if state.get('chain_run'):
            state['chain_stop'] = True
            _set_mol_status('Stopping...', spinner=True)
            return
        if state.get('band_run'):
            # A band is the longest thing this press starts -- minutes, not
            # seconds -- so being able to end it matters more here than
            # anywhere else, and ORCA writes every band it accepts to disk, so
            # what it had reached is kept.
            state['band_stop'] = True
            _set_mol_status('Stopping the band...', spinner=True)
            return
        if state.get('path_run'):
            # The walk is a single call to xtb and cannot be interrupted part
            # way, so this says so rather than lighting a Stop that would do
            # nothing until the call returned anyway.
            _set_mol_status('The walk between the two ends is one call to '
                            'xtb and cannot be cut short; it is seconds.')
            return
        if state.get('climb_run') is not None:
            # The by-hand climb has a switch of its own and it is the same
            # run, so this is that switch going up: two walks over one
            # structure would each draw half a trajectory. What it reached is
            # kept and the climb's own ending says so.
            submit_climb_btn.value = False
            return
        start = str(submit_saddle_from.value or 'here')
        how = str(submit_saddle_how.value or 'orca')
        if start == 'here':
            if how == 'hand':
                _climb_from_here()
                return
            _saddle_from_here()
            return
        ends = _path_ends(start)
        if not ends:
            return
        # And the one pair this method cannot walk between.
        #
        # All three ways from two ends go through xtb's path finder first, so
        # all three are refused together.  Decidable here in a way it is not
        # from a single structure: two ends say what the reaction is, and that
        # is why climbing from what is on screen is left alone.  See
        # :func:`_gfn.gfnff_pair_refusal` for what GFN-FF answers instead.
        if str(submit_ff_dd.value).strip().lower() == 'gfnff':
            no = _gfn.gfnff_pair_refusal(ends[0], ends[1])
            if no:
                _set_mol_status(no)
                return
        if how == 'orca':
            _path_then_orca(ends)
            return
        if how == 'neb':
            _band_between(ends)
            return
        _walk_the_path(ends, then_climb=(how == 'hand'))

    def _climb_from_here():
        """Hand what is on screen to the climb, by putting its switch down.

        There is one climb and one switch for it.  The press does not start a
        second walk beside the switch -- that is the defect e3442010 was
        about, two paths that agreed about some of a drag and not the rest --
        it puts the switch down and lets the one path run.  Pressed while the
        switch is already down and nothing is walking, which is what a climb
        that has converged and stood still looks like, it starts again from
        where the structure is now.
        """
        if submit_climb_btn.value:
            _climb_now()
            return
        submit_climb_btn.value = True

    def _saddle_from_here():
        """Climb from the structure on screen to the nearest saddle point.

        ORCA's OptTS on xtb gradients: pose a transition state by hand or take
        the estimate a walk left, press, and a few seconds later it is a
        converged first-order saddle or it is not, and which is said.
        """
        xyz = _current_xyz()
        method = str(submit_ff_dd.value)
        if not xyz:
            _set_mol_status('There is no structure to climb from.')
            return
        if method.lower() not in _saddle.SADDLE_METHODS:
            _set_mol_status(
                'A saddle search here is ORCA\'s optimiser on a semiempirical '
                'gradient, so choose GFN2, GFN1, GFN-FF or g-xTB. Anything '
                'with a basis set is a job for the ORCA Builder.')
            return
        state['saddle_run'] = True
        state['saddle_stop'] = False
        # The picture belongs to this climb and to nothing else.  A number the
        # page has not seen yet, so whatever was queued from before cannot
        # play out over it.
        state['saddle_frame_run'] = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = state['saddle_frame_run']
        _note_the_run(state['saddle_frame_run'], 'saddle')
        _ensure_manip_bootstrap()
        schedule_ui_update(_install_gfn_frame_watcher)
        submit_saddle_btn.description = 'Stop'
        submit_saddle_btn.icon = 'stop'
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        _set_mol_status('Climbing to the nearest saddle point...', spinner=True)

        def _work():
            sent = [0]

            def _watch(walked, energies):
                """Every accepted step, as it is accepted.

                Down the frame channel and not into the box: a write to the
                box rebuilds the viewer from nothing, and a climb of thirty
                steps would be thirty rebuilds.  The energy comes with the
                geometry because ORCA puts it on the comment line, so the
                status can say where the climb is as well as show it.
                """
                if state.get('gfn_run') != state.get('saddle_frame_run'):
                    return
                climb_run = state.get('saddle_frame_run')
                for n in range(sent[0], len(walked)):
                    # Asked again at the write.  An ORCA step is long
                    # enough for the run to have been replaced between
                    # the answer and the moment this reaches the field,
                    # and a frame that arrives after the page has moved
                    # on is a climb drawn over whatever replaced it.
                    schedule_ui_update(
                        lambda text=_frame_payload(
                            climb_run, **{'from': n, 'follow': 1,
                                          'frames': [[round(float(v), 4)
                                                      for v in walked[n]]]}),
                        r=climb_run: setattr(submit_gfn_frame, 'value', text)
                        if _frame_run_is_current(r) else None)
                sent[0] = len(walked)
                climbed = None
                if len(energies) > 1 and None not in (energies[0],
                                                      energies[-1]):
                    climbed = (energies[-1] - energies[0]) * _HARTREE_TO_KCAL
                schedule_ui_update(
                    _set_mol_status,
                    f'Climbing: step {len(walked)}'
                    + (f', {climbed:+.1f} kcal/mol from where it started.'
                       if climbed is not None else '.'),
                    spinner=True)

            found = _saddle.optimise_to_saddle(
                xyz, method, charge=charge, uhf=uhf, solvent=wet,
                # What this method is allowed before it is a job instead.
                # Three minutes for the ones ORCA drives itself; longer for
                # g-xTB, where every gradient is a separate process -- see
                # saddle.SECONDS_ALLOWED for the measurement.
                timeout=_saddle.seconds_for(method),
                on_frame=_watch,
                should_stop=lambda: bool(state.get('saddle_stop')))

            def _done():
                state['saddle_run'] = False
                state['saddle_stop'] = False
                _name_the_saddle_press()
                rows = [line for line in (found.get('xyz') or '').splitlines()[2:]
                        if line.strip()]
                if not found.get('ok'):
                    if rows:
                        _remember('the saddle search')
                        _write_coords(xyz_document(
                            rows, 'Where the saddle search got to'))
                    _set_mol_status(str(found.get('status')
                                        or 'The saddle search did not run.'))
                    return
                said = _saddle.verdict(found.get('imaginary'),
                                       'What it reached')
                lines = [found['status'], *said['lines']]
                if rows:
                    _remember('the saddle search')
                    # Named in the box by what it is rather than by what was
                    # being looked for.  A comment that says "transition
                    # state" over a second-order saddle is the one place the
                    # mistake would outlive the sentence that said so.
                    kept = xyz_document(rows, f'Optimised to {said["name"]}')
                    _write_coords(kept)
                    _note_the_saddle(kept, found.get('imaginary'))
                    lines.append('It is in the box; Undo takes it back.')
                    lines.extend(_the_mode_is_offered(found.get('imaginary')))
                    if said['first_order']:
                        lines.append(
                            'Refine it with OPTTS in the ORCA Builder at a '
                            'level worth quoting.')
                _set_mol_status(*lines)

            schedule_ui_update(_done)

        _start_background(_work, 'The saddle search',
                          guards={'saddle_run': False})

    def _marked_pair():
        """The end that was marked and what is on screen, or nothing.

        Two structures, marked one at a time, so nothing has to hold two at
        once.  The same structure twice is not a pair.
        """
        marked = state.get('path_from')
        here = _current_xyz()
        if marked and here and marked.strip() != here.strip():
            return (marked, here)
        return None

    def _path_ends(which='marked'):
        """The two structures to walk between, or nothing and why.

        Asked for by name.  There are two ways a pair of ends comes about --
        a structure marked beside the press, and the two ends a finished scan
        leaves -- and this used to prefer the marked one silently whenever it
        had both, so a scan walked after marking something could not be walked
        between at all.  Which one is used is now the user's answer to the box
        that says so, and this is only the fetch.
        """
        if which == 'scan':
            # Asked of the one place that decides whether the pair exists, so
            # that the press cannot act on ends the box would not have offered
            # -- they name atoms, and the molecule they name may have been
            # replaced since.
            ends = _scan_ends_here()
            if ends:
                return ends
            _set_mol_status('The last scan left no two ends to walk between '
                            'on this structure; run one, or mark two '
                            'structures.')
            return None
        pair = _marked_pair()
        if pair:
            return pair
        _set_mol_status(
            'A path needs two structures, and the one marked is the one on '
            'screen. Press Mark this end on one of them and load or build '
            'the other -- or run a scan, which leaves both.')
        return None

    def _name_the_saddle_press():
        """Say on the press what the press is about to do.

        A button that walks between two ends and stops at the estimate is not
        going to the saddle, and a name that says it is would be a promise the
        press does not keep.  Two names, and the box beside it decides which.
        """
        if str(submit_saddle_how.value) == 'walk':
            submit_saddle_btn.description = 'Find the path'
            submit_saddle_btn.icon = 'route'
        else:
            submit_saddle_btn.description = 'To the saddle'
            submit_saddle_btn.icon = 'mountain'

    def _note_the_saddle(xyz, shape):
        """Write down that this structure has modes going the wrong way.

        Said by whichever search found it -- the press, the chain, the climb --
        because all three end the same way: a geometry in the box and a Hessian
        that says what it is.  What it buys is the two presses beside them,
        which are about a saddle and are meaningless without one.

        Kept against the *geometry* and not against a flag, so nothing has to
        remember to take it away: :func:`_saddle_here` compares what was found
        with what is in the box, and an edit, an Undo, a relaxation or a load
        makes them different structures and the presses go.  A flag would have
        to be cleared in every one of those places, and the one that was
        forgotten would leave a control offering to draw the mode of a
        structure nobody has any more.

        Said after the geometry has been written and not before, because what
        it records is checked against the structure the box holds -- said
        first, it would be a saddle nobody is standing on yet and the controls
        would stay hidden.

        And it asks for the row to be redrawn itself rather than leaving that
        to the write.  A write usually redraws, and the one that does not is
        exactly the one that matters most here: a geometry the playback has
        already drawn is written with *drawn* raised, and a host that sees
        that flag keeps its viewer and returns without rebuilding anything.
        That is how the interactive climb writes its answer -- so left to the
        write, the one route where the user watched the saddle arrive would be
        the one route where nothing appeared to do anything with it.  Which is
        the failure this row has already had once: a capability arrived and
        nothing on screen said so.
        """
        order = int((shape or {}).get('count') or 0)
        modes = [float(one) for one in ((shape or {}).get('modes') or [])]
        if not xyz or order < 1 or not modes:
            state.pop('saddle_found', None)
        else:
            state['saddle_found'] = {'xyz': xyz, 'order': order,
                                     'modes': modes}
        _refresh_saddle_controls()

    def _the_mode_is_offered(shape):
        """Say out loud that there is now something to do with this mode.

        The controls arriving is the editor's own statement that a capability
        has arrived -- the visible set of presses is the answer to "what can I
        do now" -- but two new presses appearing in a row of twenty is a thing
        somebody has to notice, and a saddle search ends by reading a
        sentence.  So the sentence names them, once, at the moment they become
        possible.  This row has already had the other failure: a capability
        arrived and nothing said so, and the user asked where the buttons had
        gone.

        Nothing at all where there is no mode going the wrong way, because
        then there are no presses either and a sentence about them would be
        the report describing a row that is not on screen.
        """
        order = int((shape or {}).get('count') or 0)
        if order < 1:
            return []
        if str(submit_ff_dd.value).lower() not in _climb.CLIMB_METHODS:
            # A saddle ORCA reached on g-xTB is a saddle, and the shapes of
            # its modes are not something xtb can be asked for: g-xTB is its
            # own build and takes no ``--hess``.  So the presses are not on
            # screen and this does not name them -- it names the method that
            # would put them there, which is the actionable half.
            return ['The shapes of the modes come from xtb\'s own Hessian, so '
                    'drawing this one and following it down are offered under '
                    'GFN2, GFN1 or GFN-FF rather than under this method.']
        if order == 1:
            return ['Show the mode draws it -- an imaginary mode is the '
                    'reaction coordinate, so it says which bonds are forming '
                    'and which are breaking -- and Follow it down says which '
                    'two structures it joins.']
        return ['Show the mode draws any of them and Follow it down follows '
                'one down both ways; the box beside them says which. That is '
                'how to look at the second mode before moving along it.']

    def _saddle_here():
        """What a search found, while the box still holds the structure it
        found it on.

        The comparison is on the coordinates rather than on the text, because
        the two are not the same thing: what goes into the box carries a
        comment line this editor wrote, and what comes back out of it has been
        through a viewer.  ``largest_shift`` is the same measure the rest of
        the editor asks "is this still that structure" with.

        Absent for that reason and no other, which is the rule
        :func:`_refresh_saddle_controls` is built on.
        """
        found = state.get('saddle_found')
        here = _current_xyz()
        if not found or not here or not found.get('xyz'):
            return None
        if len(_gfn.atom_lines(found['xyz'])) != len(_gfn.atom_lines(here)):
            return None
        if _gfn.largest_shift(found['xyz'], here) > _SAME_STRUCTURE:
            return None
        return found

    #: How far an atom may have moved before the box is holding a different
    #: structure, in Angstrom.  A ten-thousandth is below the four decimals a
    #: geometry is written with, so this separates "the same coordinates,
    #: written out and read back" from every real change.
    _SAME_STRUCTURE = 1e-4

    def _which_mode():
        """Which mode the two presses are about, as an index from the softest.

        The modes are sorted from the softest up wherever they are counted --
        :meth:`climb.Climb.modes_from_xtb` sorts them and
        :func:`climb.imaginary_among` reads the list in that order -- so the
        ones going the wrong way are the first few, and the box's value is an
        index straight into both.
        """
        found = _saddle_here()
        if not found:
            return 0
        wanted = int(submit_mode_dd.value or 0)
        return wanted if 0 <= wanted < int(found.get('order') or 1) else 0

    def _modes_for(xyz, method, charge, uhf, wet, on_wait=None):
        """xtb's modes and their shapes for this structure, paid for once.

        Both presses need the same thing and neither needs it until it is
        pressed, so nothing is spent by a search that finds a saddle nobody
        looks at.  Once spent it is kept for as long as the box holds the
        structure it was taken on and the method has not changed underneath
        it: the frequencies are a property of both, and answering a press
        about GFN2 with a Hessian taken under GFN-FF would be the one kind of
        wrong that never announces itself.

        Runs on a worker, so everything it needs is handed in rather than read
        off a widget.  *on_wait* is called only where the Hessian really is
        about to be taken, which is what keeps the line that says so off a
        press that answers instantly -- a status that appears and vanishes
        inside a frame is a flicker, and a user reads it as something having
        gone wrong.
        """
        under = (str(method).lower(), int(charge), int(uhf), str(wet or ''))
        kept = state.get('saddle_modes')
        if (kept and kept.get('under') == under and kept.get('xyz')
                and len(_gfn.atom_lines(kept['xyz'])) == len(
                    _gfn.atom_lines(xyz))
                and _gfn.largest_shift(kept['xyz'], xyz) <= _SAME_STRUCTURE):
            return kept['modes']
        if on_wait is not None:
            on_wait()
        # The cores every other climb this editor starts asks for, which is
        # the module's own default: this is one xtb of the same kind, and a
        # second opinion about how much of a login node to take would only be
        # a second number to keep in step with the first.
        got = _climb.modes_of(xyz, method, charge=charge, uhf=uhf, solvent=wet)
        if got is not None:
            state['saddle_modes'] = {'xyz': xyz, 'under': under, 'modes': got}
        return got

    def _mode_press_can_run():
        """Whether either press has a structure and a method to work with."""
        if not _saddle_here():
            return False
        if (state.get('optimize_run') is not None
                or state.get('climb_run') is not None
                or state.get('saddle_run') or state.get('chain_run')
                or state.get('path_run')):
            # The picture belongs to whatever is walking.  Both of these take
            # the frame channel over, and a run that loses it draws nothing
            # for the rest of its life -- so this is refused rather than
            # allowed to quietly blank somebody's trajectory.
            _set_mol_status('Something is walking this structure, and the '
                            'picture belongs to it until it stops.')
            return False
        if str(submit_ff_dd.value).lower() not in _climb.CLIMB_METHODS:
            _set_mol_status(
                'The modes come from xtb\'s own Hessian, so this needs GFN2, '
                'GFN1 or GFN-FF. A saddle reached under another method is '
                'still in the box -- switch to one of those to draw its mode.')
            return False
        return True

    def on_submit_show_mode(_button=None):
        """Draw the imaginary mode, and never let it into the box.

        The mode is the reaction coordinate, so this is the most informative
        thing there is to show about a saddle: which bonds close and which
        open, seen rather than read off a frequency.

        Everything about the playback is the one that already exists -- the
        frame channel, the run number, the player, the interpolation between
        frames.  A mode is frames like any other path; building a second
        playback for it would be a second set of the defects the first one has
        already had.  Two things are its own: the pace, because a mode has a
        period and the slider is about how fast to walk a computed path (see
        :func:`_frame_payload`), and the ``home`` frame, because these frames
        are a picture and not a structure.
        """
        if state.get('mode_run') or state.get('down_run'):
            return
        if not _mode_press_can_run():
            return
        xyz = _current_xyz()
        method = str(submit_ff_dd.value)
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        chosen = _which_mode()
        state['mode_run'] = True

        def _work():
            modes = _modes_for(
                xyz, method, charge, uhf, wet,
                on_wait=lambda: schedule_ui_update(
                    _set_mol_status, 'Taking a Hessian to find the mode...',
                    spinner=True))

            def _done():
                state['mode_run'] = False
                if not modes:
                    _set_mol_status(
                        'The modes of this structure could not be computed, '
                        'so there is nothing to draw. It needs xtb.')
                    return
                _draw_the_mode(modes, chosen)

            schedule_ui_update(_done)

        _start_background(_work, 'The mode',
                          guards={'mode_run': False})

    def _draw_the_mode(modes, chosen):
        """Put the mode on the frame channel, once, as a whole picture.

        One write and not a stream: the frames are arithmetic on a geometry
        that is already in hand, so there is nothing to wait for and nothing
        to watch arriving.  It is marked as the run's last write, which is
        what stops the player throwing the queue away when it notices that no
        Optimise switch stands behind it -- there is none, and there is
        nothing running for one to stand behind.

        The amplitude is cut to fit this molecule before anything is drawn.
        The default is a picture that reads clearly on the case it was
        measured on; a mode that drives two atoms together on some other
        structure would tear it, and a torn picture is worse than a small one
        because it looks like chemistry.
        """
        cm = list(modes['cm'])
        wrong = _climb.imaginary_among(cm)
        if not wrong:
            _set_mol_status(
                'No mode of the structure in the box goes the wrong way any '
                'more, so there is no reaction coordinate to draw.')
            # And the presses go with it: what the search reported and what a
            # Hessian says now disagree, and the Hessian is the later of the
            # two.
            _note_the_saddle(None, None)
            return
        picked = chosen if chosen in wrong else wrong[0]
        here = modes['angstrom']
        way = modes['ways'][:, picked].reshape(len(modes['symbols']), 3)
        amplitude = _climb.amplitude_that_fits(modes['symbols'], here, way)
        frames = _climb.mode_pictures(here, way, amplitude=amplitude)
        run = _claim_the_frame_run()
        _ensure_manip_bootstrap()
        schedule_ui_update(_install_gfn_frame_watcher)
        seconds = (len(frames) * _climb.MODE_PACE_MS) / 1000.0
        schedule_ui_update(
            lambda text=_frame_payload(
                run, pace=_climb.MODE_PACE_MS,
                **{'from': 0, 'frames': frames, 'final': 1,
                   'home': frames[-1]}),
            r=run: setattr(submit_gfn_frame, 'value', text)
            if _frame_run_is_current(r) else None)
        _set_mol_status(
            f'Drawing the mode at {float(cm[picked]):.0f} cm-1: '
            f'{_climb.MODE_SWINGS} swings, about {seconds:.0f} s, and it '
            'stops on the structure it started from.',
            'The atoms it moves are the atoms the reaction moves. It is a '
            'picture and nothing else: the coordinates in the box do not '
            'change, and taking hold of the structure puts it straight back.'
            + ('' if amplitude >= _climb.MODE_AMPLITUDE else
               f' Drawn at {amplitude:.2f} A rather than the usual '
               f'{_climb.MODE_AMPLITUDE:.2f} -- further than that and this '
               'molecule has two atoms inside each other.'))

    def on_submit_follow_down(_button=None):
        """Push the structure down the mode both ways and say where it lands.

        One imaginary mode makes a structure *a* transition state.  Whether it
        is the one for the reaction that was meant is the other question, and
        this is the standard way to ask it.  What comes back is described
        rather than named -- how many separate pieces each end is, which bonds
        it has that the saddle did not, what it costs against the saddle --
        because every sort of system is computed here and a sentence that
        assumed a kind of chemistry would be wrong about most of them.

        The saddle stays in the box.  Two structures came out and one box
        holds one structure, and the one being worked on is the saddle: these
        two are what it *joins*, which is a fact about it rather than a
        replacement for it.
        """
        if state.get('mode_run') or state.get('down_run'):
            return
        if not _mode_press_can_run():
            return
        xyz = _current_xyz()
        method = str(submit_ff_dd.value)
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        chosen = _which_mode()
        state['down_run'] = True
        _set_mol_status('Following the mode down both ways...', spinner=True)

        def _work():
            modes = _modes_for(
                xyz, method, charge, uhf, wet,
                on_wait=lambda: schedule_ui_update(
                    _set_mol_status, 'Taking a Hessian to find the mode...',
                    spinner=True))
            got = _climb.follow_the_mode_down(
                xyz, method, mode=chosen, charge=charge, uhf=uhf,
                solvent=wet, modes=modes,
                # The same allowance the saddle search itself has: this is two
                # relaxations of a structure that search has already converged
                # on, so anything it cannot do inside that is a job.
                timeout=_saddle.seconds_for(method),
                on_stage=lambda said: schedule_ui_update(
                    _set_mol_status, said, spinner=True))

            def _done():
                state['down_run'] = False
                lines = [str(got.get('status') or 'It did not run.')]
                lines.extend(str(one) for one in (got.get('lines') or ()))
                lines.append(
                    'The saddle is still in the box; these two are what it '
                    'joins, not a replacement for it.'
                    if got.get('ok') else
                    'Nothing was written to the box.')
                _set_mol_status(*lines)

            schedule_ui_update(_done)

        _start_background(_work, 'Following the mode down',
                          guards={'down_run': False})

    def _scan_left_two_ends():
        """A scan has finished: put the press on the pair it just made.

        Two buttons used to *appear* when a scan left two ends -- Find the
        path and Path to saddle -- and their arriving was the editor saying
        that something new could be done now.  Folded into the box that asks
        where a search starts (085cfd55), the entry appeared and nothing else
        did: the box went on standing at "what is on screen", so the press
        went on meaning climb what is in the box, and the one moment the user
        was being told about passed without anything happening on screen.
        The user, on that: "nach einem Scan konnten wir doch davor noch weiter
        den Path absuchen mit Buttons -- wo sind die hin".

        So the start moves.  A finished scan is the strongest statement there
        is about what the user is interested in now -- it is minutes of
        walking, and its two ends are what the walk was for -- and the press
        beside it then means the thing that walk made possible, one press
        away.  What changes on screen is the start naming the scan's ends and
        the second box arriving beside it with the ways that are open from a
        pair, which is the same news the two buttons used to carry.

        Nothing is taken away: "what is on screen" is still in the list, one
        selection away, and every way that could run from it can still run.

        Moved through the wish rather than by writing the box, because that is
        what the wish is for: a value that cannot be met is kept aside and put
        back the moment it can (see :func:`_keep_the_choice`).  Under a method
        with no route from a pair the start simply stays where it is and the
        move happens later, if the user picks a method that has one.
        """
        state['saddle_start_wish'] = 'scan'
        _refresh_saddle_controls()

    def _scan_ends_here():
        """The two ends the last scan left, while they are about this molecule.

        They outlive the structure they were walked on -- a scan leaves them
        and loading something else does not take them away -- and two ends
        belonging to another molecule are not a pair for this one: the walk
        between them would carry the previous structure into the box under the
        name of a path.  Told apart by the element column, the way the held
        values and the thermal anchor are, because that is what makes two
        geometries the same molecule.

        Absent for that reason and no other, which is the rule the whole of
        :func:`_refresh_saddle_controls` is built on.
        """
        ends = state.get('scan_ends')
        if not ends:
            return None
        here = _current_xyz() or ''
        if _structure_fingerprint(ends[0]) != _structure_fingerprint(here):
            return None
        return ends

    def _refresh_saddle_controls():
        """Offer the starts that exist and the ways that can run, and no more.

        Composed with :func:`_refresh_method_controls` rather than fighting
        it: that one decides what the chosen engine can do at all and calls
        this, and this decides what is left to choose between.  Nothing here
        is hidden for being advanced -- every option that is absent is absent
        because the thing it names does not exist at this moment: a pair of
        ends nobody has made, a walk with nothing to walk between, an engine
        this method has no keyword for.

        A box with one entry is not a choice, so it is not on screen.  The
        press alone then means the only thing it can mean, which is what the
        editor opens on: climb what is in the box.
        """
        if state.get('saddle_controls_quiet'):
            # Rewriting a box's entries can move its value, which comes back
            # through the observer that called this.  Once through is enough.
            return
        state['saddle_controls_quiet'] = True
        try:
            chosen = str(submit_ff_dd.value)
            starts = [('what is on screen', 'here')]
            if _marked_pair():
                starts.append(('the end you marked', 'marked'))
            if _scan_ends_here():
                starts.append(("the scan's two ends", 'scan'))

            def ways_from(start):
                """The ways to a saddle that can run from this start."""
                out = []
                if chosen.lower() in _saddle.SADDLE_METHODS:
                    out.append(('through ORCA', 'orca'))
                # By hand and the walk alone are answers only where there is
                # a walk.  From the structure on screen there is nothing
                # between two ends to stop after -- and climbing what is on
                # screen by hand is Climb to TS, which is a switch over with
                # the other switches and one press either way.  Offered here
                # as well it would be the same thing under two names in two
                # places, which is the whole complaint this is answering.
                if start != 'here':
                    # A band, which is another *how* from a two-ended start
                    # and not a button of its own: it answers the same
                    # question ORCA and by-hand answer, from the same pair, and
                    # a fourth press beside them would be the third time this
                    # row learnt that lesson.
                    #
                    # After ORCA rather than before it, because the order of
                    # the list is the recommendation.  Measured on the
                    # sixteen-atom Diels-Alder from the same two ends, the two
                    # reach the same saddle to 0.07 cm-1 and the band spends
                    # 203 gradients doing it -- so it is what to run when the
                    # cheap answer is not believed, or when the two ends are
                    # too far apart for a cheap interpolation to bridge, and
                    # it is not what to run first.
                    if chosen.lower() in _saddle.SADDLE_METHODS:
                        out.append(('through NEB-TS', 'neb'))
                    if chosen.lower() in _climb.CLIMB_METHODS:
                        out.append(('by hand', 'hand'))
                    if _gfn.is_gfn_method(chosen):
                        out.append(('the path only', 'walk'))
                return out

            # A start this engine can do nothing from is not offered, which is
            # what keeps the two boxes from arguing with each other.  Under a
            # method with neither saddle optimiser -- g-xTB is its own build
            # and ORCA cannot drive it -- the structure on screen has no way
            # up at all while two ends still have the walk between them, so
            # the list of starts is the pair and the press is the walk.  Left
            # to the boxes to sort out between them, the start stood at "what
            # is on screen", found no way, hid the press, and hid with it the
            # box that was the only way to choose the start that worked.
            starts = [one for one in starts if ways_from(one[1])]
            _keep_the_choice(submit_saddle_from, starts,
                             'saddle_start_wish')
            ways = ways_from(str(submit_saddle_from.value)) if starts else []
            _keep_the_choice(submit_saddle_how, ways,
                             'saddle_way_wish')
            # Nothing this method can do from anywhere, so nothing to press.
            # The mark stays: it describes two structures rather than a
            # program, and survives a change of method the way an armed scan
            # does.
            can_run = bool(starts and ways)
            submit_saddle_btn.layout.display = '' if can_run else 'none'
            submit_saddle_from.layout.display = (
                '' if can_run and len(starts) > 1 else 'none')
            submit_saddle_how.layout.display = (
                '' if can_run and len(ways) > 1 else 'none')
            submit_path_from_btn.layout.display = (
                '' if _gfn.is_gfn_method(chosen) else 'none')
            # And what there is to do with a saddle, which exists exactly when
            # a search has found one and the box still holds it.  Not "when a
            # transition state was found": a structure with two modes going
            # the wrong way has two reaction coordinates to look at and the
            # verdict's own advice is to move along the second of them, which
            # nobody can do without seeing it.  What is *not* offered is the
            # ordinary case, a minimum -- there is no mode going the wrong
            # way there and both presses would have nothing to say.
            found = _saddle_here()
            reads_modes = str(chosen).lower() in _climb.CLIMB_METHODS
            has_modes = bool(found) and reads_modes
            for button in (submit_mode_btn, submit_ends_btn):
                button.layout.display = '' if has_modes else 'none'
            # One mode is not a choice, which is the same rule the two boxes
            # above follow.  The frequencies are the verdict's own, so the box
            # names what the sentence named.
            order = int((found or {}).get('order') or 0)
            named = [float(one) for one in ((found or {}).get('modes') or ())]
            _keep_the_choice(
                submit_mode_dd,
                [(f'{one:.0f} cm-1', n) for n, one in enumerate(named)],
                'saddle_mode_wish')
            submit_mode_dd.layout.display = (
                '' if has_modes and order > 1 and len(named) > 1 else 'none')
            _name_the_saddle_press()
        finally:
            state['saddle_controls_quiet'] = False

    def _keep_the_choice(box, options, wish):
        """Rewrite a box's entries without throwing away what was chosen.

        A dropdown handed a list its value is not in loses that value, and
        losing it silently is how a setting somebody made comes back as the
        default one press later.  What is on these two lists is decided by the
        method and by which structures exist, and both of those change while
        the user is doing something else -- so the wish is kept aside under
        *wish* while it cannot be met and put back the moment it can.

        An empty list is left alone rather than written over with a placeholder
        entry: there is nothing on screen then anyway, and blanking it is how
        "through ORCA" became "by hand" on the way through a method that could
        run neither.
        """
        if not options:
            return
        want = state.get(wish) or box.value
        values = [value for _label, value in options]
        if list(box.options) != [tuple(one) for one in options]:
            # Only when they have actually changed.  This runs on every redraw
            # of the structure, and a scan redraws once a point: assigning the
            # same list each time is a widget message a point for a box nobody
            # is looking at.
            box.options = options
        if want in values:
            box.value = want
            state.pop(wish, None)
        else:
            state[wish] = want
            box.value = values[0]

    def _name_the_pair_the_hand_means(holding):
        """Which two atoms this gesture is about, when the user has said.

        The climb a release starts is judged on whether what it reached is
        *the reaction the hand pointed at*, and that is a question about two
        atoms -- see :func:`climb.reach_the_reaction`, where three searches
        reach 12, 12 and 10 of twenty-one hand drags and sixteen between
        them, and only because each one's answer can be checked against the
        pair.  Without a pair only the count of imaginary modes can be asked;
        every later rung would then be judged by the very test the first one
        has already passed, so the ladder stops at one rung and says which
        test it used.

        The page names the atom under the cursor and the atoms that are
        picked, both on messages it already sends and both counted the same
        way -- ``ffIndicesOf``, from nought.  So the pair is: the atom being
        dragged, and the atom the user tapped to say what they are dragging it
        at.  Picking two atoms and dragging them together names the pair
        outright and is taken as well.

        Nothing is guessed, and that is the measured part.  Seven element-free
        rules were scored against twenty-one hand drags, given the atom the
        page names and both geometries: the contact the follow is actually
        holding -- :func:`gfn_optimize.contacts_holding`, which is the best of
        them -- names the right pair on **10 of the 21**, and the other six on
        1, 1, 1, 3, 7 and 9.  They fail in three shapes: a drag at a carbon
        names the hydrogen on it, a drag that folds a chain names the torsion
        it turns rather than the contact it closes, and a proton walked
        between two oxygens names the one it started on.  A *wrong*
        pair is worse than none, because it is the ladder rejecting an answer
        that was right: measured on the van-der-Waals complex, the aimed climb
        reached the Diels-Alder saddle at 2.315 A and -394 cm-1, the check
        against C10-H4 scored it below a fifth, and the press then spent 73
        seconds on two more searches to say it had not found what was already
        on the screen.

        A pick the user has forgotten about would do the same, and that is the
        one case this cannot rule out.  It is a different kind of wrong,
        though: the selection is drawn on the structure, and the sentence at
        the end names the contact it checked against -- so it is visible
        before the press and said after it, which a rule nobody chose is not.
        """
        hand = [int(i) for i in (holding or ())]
        picked = [int(i) for i in (state.get('picked') or ())]
        if len(hand) == 2 and set(hand) == set(picked):
            # Both ends picked and dragged together.  It aims at nothing --
            # a rigid translation is projected out of the Hessian before the
            # gesture is matched against it -- but it names the contact, which
            # is the other half of what the ladder needs.
            state['climb_held'] = (hand[0], hand[1])
            return
        aimed = [i for i in picked if i not in hand]
        if len(hand) == 1 and len(aimed) == 1:
            state['climb_held'] = (hand[0], aimed[0])

    def _climb_can_run():
        """Whether an interactive climb is possible here, said out loud if not.

        Read before the toggle is believed and again before every run the
        release starts, because the method box can change between the two.
        """
        if not _current_xyz():
            _set_mol_status('There is no structure to climb from.')
            return False
        if str(submit_ff_dd.value).lower() not in _climb.CLIMB_METHODS:
            _set_mol_status(
                'An interactive climb runs on xtb, so choose GFN2, GFN1 or '
                'GFN-FF. Anything with a basis set is a job for the ORCA '
                'Builder.')
            return False
        return True

    def on_submit_climb(change=None):
        """Which way the optimiser walks, and a walk started at the press.

        A mode, not a run. Dynamik Opt follows the hand, Auto carries on when
        the hand lets go, and this says which way "carry on" means -- so it
        stays down across the runs it starts, exactly as those two do. It used
        to lift itself the moment a climb converged, which made it a one-shot:
        measured on the van-der-Waals complex, the drag after that walked
        straight back down to 3.353 A because the release was a minimisation
        again and nothing said so.

        Pressing it also starts a walk at once, the way pressing Optimize
        does, so the button is never a setting that does nothing until the
        next gesture. Pressing it again stops what is walking and puts the
        release back downhill.

        Everything else about it is Dynamik Opt's: the same follow under the
        hand, the same interrupt at the grab, the same restart after the
        release, the same run number and the same frame channel. See
        :func:`_after_release`, which is the one place the two are told apart.
        """
        if change is not None and change.get('name') != 'value':
            return
        if not submit_climb_btn.value:
            if state.get('climb_run') is not None:
                # Pressed while it was still walking, so this is a Stop, and
                # the same one Optimise's switch is: the run number moves on,
                # which refuses every write the stopped run still has in hand,
                # and the loop hears it and halts the page.  What is kept is
                # the frame the picture is standing on -- the steps it walked
                # past that one were never seen and were not chosen.
                _claim_the_frame_run()
            state['climb_run'] = None       # the loop reads this and stops
            state.pop('climb_interrupted', None)
            state.pop('climb_was', None)
            state.pop('climb_held', None)
            return
        if not _climb_can_run():
            submit_climb_btn.value = False
            return
        # A running optimisation and a climb are two answers to the same
        # question over one structure.  The press is the later of the two
        # things the user asked for, so it wins and the other is told to stop.
        if _interrupt_gfn():
            state.pop('climb_interrupted', None)
            _stand_down_after_interrupt(
                'The optimisation stopped: a climb and a minimisation cannot '
                'own one structure.')
        _climb_now()

    def _climb_now():
        """Walk to a transition state slowly enough to be interrupted.

        The same method ORCA's OptTS uses, run here on xtb gradients so that
        each step costs about ten milliseconds instead of the three seconds a
        fresh ORCA needs to start. That is the whole reason this exists beside
        the press: a press cannot be dragged in the middle of.

        Three searches rather than one, cheapest first, and each one's answer
        checked against the pair of atoms the hand held --
        :func:`climb.reach_the_reaction`, where the measurement is. A saddle
        search begun by hand does not fail; it succeeds at arriving somewhere,
        and usually somewhere else: all three of these converge, confidently
        and with exactly one mode going the wrong way, onto methyl torsions at
        -48 cm-1. So what is asked of what comes back is not "is it a saddle"
        but "is its imaginary mode the contact you were dragging", and a rung
        that misses hands the *hand's* structure to the next one rather than
        its own endpoint. The pair is what makes that question askable, and
        where it comes from is :func:`_name_the_pair_the_hand_means`; without
        one only the mode count can be checked, every rung would be judged by
        the test the first has already passed, and the ladder stops at one.

        Started from three places and identical from all three -- the toggle
        going down, a release with Auto on, and a release bringing back the
        climb a grab interrupted. What it starts from is whatever is on screen
        now, and what it aims along is ``climb_was``: where the structure stood
        when the hand arrived, so that the difference between then and now is
        the direction the user asked for. That direction, and not the pull, is
        what guides the climb -- from the same dragged geometry, climbing the
        mode the hand pointed at reaches the Diels-Alder saddle in 39 steps at
        2.315 A, and climbing the lowest mode instead walks back down to the
        van-der-Waals complex 0.43 A away with no imaginary mode at all.

        The pull itself is never part of it. A saddle of a restrained surface
        is not a saddle of the real one: measured on the Diels-Alder with both
        forming bonds held, at 2.20 A it ends 0.53 A from the saddle with two
        imaginary modes and a true gradient 138 times the convergence
        threshold, and at 2.60 A it converges in five steps onto a point with
        no imaginary mode at all.
        """
        if state.get('climb_run') is not None:
            return                      # one is already walking
        if not submit_climb_btn.value:
            return                      # switched off while this was waiting
        if not _climb_can_run():
            submit_climb_btn.value = False
            return
        xyz = _current_xyz()
        method = str(submit_ff_dd.value)
        # Where the structure stood when the hand arrived, if a hand is what
        # brought us here.  Taken once: a second run from the same mark would
        # aim along a gesture that has already been answered.
        aimed_from = state.pop('climb_was', None)
        # And which two atoms that gesture was about, written down while the
        # drag was running -- see :func:`_name_the_pair_the_hand_means`.
        # Taken once for the same reason: it belongs to the gesture, and the
        # gesture is answered here.
        held = state.pop('climb_held', None)
        token = object()
        state['climb_run'] = token
        state.pop('climb_interrupted', None)
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        _remember('the climb to a transition state')
        _ensure_manip_bootstrap()
        schedule_ui_update(_install_gfn_frame_watcher)
        # The run number, claimed where the minimisation claims its own: when
        # the run begins.  See _claim_the_frame_run for the 94 frames that
        # were thrown away when the climb claimed at the toggle instead and
        # then held the number through a drag.
        run = _claim_the_frame_run()
        state['climb_frame_run'] = run
        # The same two things the minimisation clears when it claims a run.
        # A halt marked as already sent -- by the very grab that brought this
        # climb back -- would swallow this run's own Stop, and a frame number
        # left over from an earlier grab is a plausible index into a path that
        # was walked by something else.
        state['gfn_halt_sent'] = False
        _forget_the_shown_frame()
        said = ('Carrying on from where you let go...' if aimed_from
                else 'Climbing to a transition state; drag an atom to point '
                     'it somewhere...')
        state['gfn_last_status'] = said
        _set_mol_status(said, spinner=True)
        walked: list = []

        def _push(final=False):
            """The same writer the minimisation uses, held further apart.

            A climb makes a hundred frames a second where xtb's optimiser
            makes a few, and every write is a message.

            Not marked as a follow.  It was, and for one reason only: the page
            abandons a queue when the Optimise switch is up, and nothing on
            the page knew that Climb to TS is a switch too, so the mark was
            what kept the climb's path from being thrown away.  It cost the
            slider.  A followed hand is paced by how fast xtb answers -- that
            is what the mark means to the player -- so at the top of the
            slider the climb was drawn at its arrival rate, one frame an
            animation frame, while the minimisation next to it drained
            everything that had arrived.  Measured over the same release,
            Speed 60: the minimisation drew 56 frames in 3 draws and never had
            anything in hand, the climb drew 113 in 71 and always did.  The
            page reads the climb's own switch now -- see switchIsOn -- so the
            mark is not needed and the two are paced by one rule.
            """
            _stream_frames(run, walked, final=final, least_apart=0.04)

        def _work():
            began = time.time()
            got = None
            rows = []
            steps, seconds = 0, 0.0
            # Whether the switch ended this, which is a different ending from
            # both of the other two: a run that arrived somewhere has a result
            # to name, one that failed has a reason, and one that was stopped
            # has neither -- it has the frame the user was looking at.
            switched_off = [False]
            # Frames drawn by the climbs, and how many of ORCA's have been
            # taken.  ORCA hands its whole trajectory over each time it grows
            # where a climb hands one step, so one of the two counts and the
            # other remembers where it had got to.
            counted = [0]
            from_orca = [0]

            def _stopped():
                """Whether this walk is still the one, and the page told if not.

                The same two lines the minimisation's loop reads, through the
                same halt: one Stop, one thing it means.  The climb had no
                such line at all -- it stopped walking and said nothing -- so
                the page went on drawing the path it already held.  Measured,
                a Stop at frame 13 of 117 with the slider at twelve a second
                drew 509 more times afterwards, frames 14 through 116.

                Asked by all three searches, so a Stop ends the press and not
                merely the rung that happens to be running.
                """
                halted = state.get('climb_run') is not token
                if halted:
                    _halt_the_frames(run)
                return halted

            def _drew(frame, _outcome):
                """One climb step, from whichever of the two climbs made it.

                Both rungs report here, so when the second one starts the
                picture jumps back to the structure the hand made and climbs
                away from it again.  That is what the second rung *is* -- every
                rung is handed the hand's own structure and never the last
                one's, because a search that is going wrong does not stop
                somewhere useful -- and the sentence that goes out before it
                says so.

                The run number is not troubled by the jump.  One press is one
                run: the frames go into one list in the order they are drawn,
                and the count the page reports back still indexes that list,
                which is what a Stop cuts the path at.

                What arrives is a flat coordinate list rather than xyz text,
                because that is what the frame channel speaks, so the geometry
                a grab reads off is rebuilt from it here.
                """
                walked.append([float(one) for one in frame])
                showing = _frame_as_xyz(xyz, walked, len(walked),
                                        'Climbing to a transition state')
                if showing:
                    state['climb_showing'] = showing
                _push()
                counted[0] += 1
                # Said while it walks, in the row every other calculation
                # writes.  Silent, the page's own playback notes rewrote
                # that row with whatever the drag had last said, so a
                # climb thirty seconds in was reported as a follow that
                # had finished long ago.
                if counted[0] % 8 == 0:
                    state['gfn_last_status'] = (
                        f'Climbing to a transition state: {counted[0]} '
                        f'step(s).')
                    schedule_ui_update(_set_mol_status,
                                       state['gfn_last_status'], spinner=True)

            def _walked(path, _energies):
                """ORCA's trajectory, which arrives whole each time it grows.

                The last rung is a different optimiser and reports differently
                -- every accepted step, all of them, each time -- so what is
                new is taken off the end and put on the one list the page is
                playing.  Pretending the two report the same way would lose
                one of them.
                """
                for n in range(from_orca[0], len(path)):
                    walked.append([float(one) for one in path[n]])
                from_orca[0] = len(path)
                _push()

            def _route(sentence):
                """A rung's own sentence, said before the rung it explains.

                The module writes these and they are shown as they come: a
                retry is seconds of a still picture otherwise, and a user
                watching one has no way of telling a search that has moved on
                to its next idea from one that has hung.
                """
                state['gfn_last_status'] = sentence
                schedule_ui_update(_set_mol_status, sentence, spinner=True)

            try:
                got = _climb.reach_the_reaction(
                    xyz, method, held=held, charge=charge, uhf=uhf,
                    solvent=wet, aimed_from=aimed_from,
                    max_steps=_climb.CLIMB_CEILING,
                    on_frame=_drew, on_path=_walked, on_route=_route,
                    should_stop=_stopped,
                    # What the last rung is allowed before it is a job
                    # instead, taken from the same place the press next door
                    # takes it.
                    fallback_timeout=_saddle.seconds_for(method))
                cut_by_a_hand = state.get('climb_cut') is token
                if cut_by_a_hand or _stopped():
                    # Nothing to name.  A climb a hand cut off is about a
                    # structure the user has since replaced, and one the
                    # switch stopped ended at the frame the picture was
                    # showing rather than where the walk had got to -- so a
                    # verdict on where the walk stands is a verdict on a
                    # geometry nobody is being left with.
                    #
                    # This is Optimise's Stop, said in the climb's words:
                    # neither of them reports a result, both keep the frame on
                    # screen, and the path is put down for the page's own
                    # report to cut.
                    if not cut_by_a_hand:
                        switched_off[0] = True
                        steps, seconds = counted[0], time.time() - began
                else:
                    _push(final=True)
                    steps, seconds = counted[0], time.time() - began
                    rows = [line for line in
                            str(got.get('xyz') or '').splitlines()[2:]
                            if line.strip()]
            except Exception as trouble:            # noqa: BLE001
                got, rows = None, []
                steps, seconds = 0, time.time() - began
                state['climb_error'] = str(trouble)

            def _done():
                if state.get('climb_cut') is token:
                    # A hand took the structure while this was walking or
                    # while it was finishing.  Whatever it reached is older
                    # than what the user has since made, and writing it over
                    # an edit that came afterwards is the one thing an editor
                    # may never do.  The path that interrupted it has already
                    # said what happens next.
                    state.pop('climb_error', None)
                    return
                if state.get('climb_run') is token:
                    state['climb_run'] = None
                # What a walk between two ends said, when this climb is the
                # second half of that press.  Its own opening line wrote over
                # the walk's report the moment it started, so the numbers are
                # kept and said again at the end: a barrier and the saddle it
                # belongs to are one answer, and one replacing the other on
                # the way past is how the barrier came to be lost.
                walked_said = list(state.pop('path_said', None) or ())
                if switched_off[0]:
                    # The switch, not an arrival and not a failure.  What is
                    # kept is the frame the picture stopped on and not where
                    # the walk had got to -- the same thing a Stop means to
                    # the minimisation, put down here for the page's own
                    # report to cut.
                    _the_picture_stopped_here(
                        run, xyz, walked, 'stopped at the frame on screen')
                    said = (f'The climb stopped at the frame you were looking '
                            f'at, {steps} step(s) in ({seconds:.1f} s). Press '
                            f'Climb to TS again to carry on from there.')
                    state['gfn_last_status'] = said
                    _set_mol_status(*walked_said, said)
                    return
                if got is None:
                    _set_mol_status('The climb could not run: '
                                    + str(state.pop('climb_error', '')))
                    return
                # Whether the picture already has this.  The frames go out
                # under the run number this walk claimed, so if that number is
                # still the current one they landed and a redraw would only
                # tear down what is playing; if something moved it on, they
                # were refused and the box has to draw itself.  Assumed
                # instead, a climb whose final frames were refused left the
                # viewer standing on a geometry the box no longer held.
                #
                # And not at all while a hand is on the structure.  Switching
                # the toggle off in the middle of a drag reaches here with a
                # geometry older than what the mouse is doing, and the grab
                # that would have marked it as cut never happened -- the climb
                # was already stopped by then.
                holding = bool(state.get('gfn_follow'))
                if rows and not holding:
                    kept = xyz_document(rows, 'Climbed towards a transition '
                                              'state')
                    _write_coords(kept, drawn=_frame_run_is_current(run))
                    # The climb's own Hessian on what it reached, which is the
                    # same question the press next door asks and is answered
                    # here already: the two presses beside them read it, so a
                    # saddle reached by hand offers exactly what a saddle
                    # reached by pressing offers.  This is also the write that
                    # a host does not redraw -- see :func:`_note_the_saddle`.
                    _note_the_saddle(kept, got.get('imaginary'))
                # The module's own sentence, whichever rung ended the press:
                # it is the one that knows which searches were tried, what
                # each reached and whether the mode it found is the contact
                # the hand held.  Written again here would be a second opinion
                # on one structure, and the wrong one is always the one that
                # is easier to reach.
                said = str(got.get('status') or '')
                lines = walked_said + [said]
                if held is None:
                    # The module has already said that nothing named the
                    # contact; what it cannot say is how to name one, because
                    # that is the editor's own gesture.  Said whether or not
                    # this press reached something: with no pair the other two
                    # searches are never tried at all, so it is the difference
                    # between one search and three either way.
                    lines.append(
                        'Tap the atom you are aiming at before you drag '
                        'another one at it, and the climb can check what it '
                        'reaches against that contact -- and try two more '
                        'searches when the first one misses.')
                if rows and not holding:
                    lines.extend(_the_mode_is_offered(got.get('imaginary')))
                lines.append(
                    'You moved the structure while it was finishing, so what '
                    'you made is what is in the box.' if holding else
                    'Drag an atom to point it at another one, or switch Climb '
                    'to TS off to go back down to minima. Undo takes the '
                    'whole climb back; refine it with OPTTS in the ORCA '
                    'Builder at a level worth quoting.')
                # The climb's own headline, not the walk's, because this is
                # what the row is asked to say again whenever the picture
                # redraws itself.
                state['gfn_last_status'] = said
                _set_mol_status(*lines)

            schedule_ui_update(_done)

        _start_background(_work, 'The climb to a transition state',
                          guards={'climb_run': None})

    def on_submit_path_from(_button=None):
        """Mark what is on screen as where a path starts.

        The other end is whatever is on screen when the press is made, so the
        two are marked one at a time and nothing has to hold two structures at
        once.  Marking is what puts "from the end you marked" in the box
        beside the press, which is where the pair is then chosen.
        """
        xyz = _current_xyz()
        if not xyz:
            return
        state['path_from'] = xyz
        _refresh_saddle_controls()
        _set_mol_status(
            f'Marked as one end of a path ({len(_gfn.atom_lines(xyz))} '
            'atoms). Load or build the other structure, and the box beside '
            'the press will offer to start from the pair.')

    def _path_then_orca(ends):
        """Two structures in, a converged saddle out, at one press.

        The path finder and the saddle optimiser are two halves of one
        question and neither engine holds both: xtb walks between two ends and
        estimates the top of what it crossed, and has no saddle optimiser at
        all; ORCA has one and nothing that produces an estimate to give it.
        Chained, the pair is twelve seconds on sixteen atoms -- and lands
        within a tenth of a wavenumber of a nudged elastic band on the same
        two ends, for about half the gradients.  Which is why the band is the
        other entry in the box beside this and not the default in it.

        Watched while it climbs and stopped by the same press, like the climb
        on its own.  The walk itself is not drawn: it ends at the product and
        the climb starts from the highest point of it, so playing both would
        run the picture forward and then jump back, and the picture is meant
        to be believable.

        What was reached is named rather than assumed, which is the half that
        cannot be left out -- a saddle search does not fail when it goes
        wrong, it succeeds at arriving somewhere.

        Stopping it after the walk and before the climb is what "stop at the
        estimate" asks for outright, and the two arrive at the same place:
        what the walk found is kept and said.
        """
        if state.get('path_run') or state.get('saddle_run'):
            _set_mol_status('Something is already running on this structure; '
                            'let it finish or stop it first.')
            return
        method = str(submit_ff_dd.value)
        if method.lower() not in _saddle.SADDLE_METHODS:
            _set_mol_status(
                'The path is xtb\'s and the climb is ORCA on a semiempirical '
                'gradient, so both halves want GFN2, GFN1, GFN-FF or g-xTB. '
                'Anything with a basis set is a job for the ORCA Builder.')
            return
        state['chain_run'] = True
        state['chain_stop'] = False
        # The picture belongs to this run and to nothing else: a number the
        # page has not seen yet, so whatever was queued from before cannot
        # play out over it.
        state['chain_frame_run'] = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = state['chain_frame_run']
        _note_the_run(state['chain_frame_run'], 'chain')
        _ensure_manip_bootstrap()
        schedule_ui_update(_install_gfn_frame_watcher)
        submit_saddle_btn.description = 'Stop'
        submit_saddle_btn.icon = 'stop'
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        model = _solv_model()
        label = _server_label(method)
        _set_mol_status(f'{label} is looking for a path between the two '
                        'ends...', spinner=True)

        def _work():
            sent = [0]

            def _watch(walked, energies):
                """Every accepted step of the climb, as it is accepted.

                Down the frame channel and not into the box: a write to the
                box rebuilds the viewer from nothing, and a climb of thirty
                steps would be thirty rebuilds.
                """
                if state.get('gfn_run') != state.get('chain_frame_run'):
                    return
                for n in range(sent[0], len(walked)):
                    schedule_ui_update(
                        lambda text=json.dumps({
                            'run': state.get('chain_frame_run'),
                            'from': n,
                            'follow': 1,
                            'frames': [[round(float(v), 4)
                                        for v in walked[n]]],
                        }): setattr(submit_gfn_frame, 'value', text))
                sent[0] = len(walked)
                climbed = None
                if len(energies) > 1 and None not in (energies[0],
                                                      energies[-1]):
                    climbed = (energies[-1] - energies[0]) * _HARTREE_TO_KCAL
                schedule_ui_update(
                    _set_mol_status,
                    f'Climbing from the estimate: step {len(walked)}'
                    + (f', {climbed:+.1f} kcal/mol from where it started.'
                       if climbed is not None else '.'),
                    spinner=True)

            found = _saddle.path_to_saddle(
                ends[0], ends[1], method, charge=charge, uhf=uhf,
                solvent=wet, solvation_model=model, on_frame=_watch,
                # The climbing half's allowance, which is the method's --
                # g-xTB takes a good deal longer than the xtb methods for the
                # reason saddle.SECONDS_ALLOWED gives.  The walking half's is
                # already generous enough for either: measured on sixteen
                # atoms, 13 s under GFN2 and 105 under g-xTB against 600.
                timeout=_saddle.seconds_for(method),
                should_stop=lambda: bool(state.get('chain_stop')),
                on_stage=lambda said: schedule_ui_update(
                    _set_mol_status, said, spinner=True))

            def _done():
                state['chain_run'] = False
                state['chain_stop'] = False
                _name_the_saddle_press()
                lines = []
                if found.get('barrier') is not None:
                    near = found.get('rmsd')
                    lines.append(
                        f'{label} walked a path in '
                        f'{found["path_seconds"]:.1f} s: '
                        f'{found["barrier"]:.1f} kcal/mol forward'
                        + (f', {found["back"]:.1f} back'
                           if found.get('back') is not None else '')
                        + '.'
                        + ('' if near is None else
                           f' It came within {near:.2f} A RMSD of the '
                           'structure it was aiming at.'))
                # What is shown is what the climb reached, or -- if it never
                # got that far -- the estimate the walk left, which is still
                # an answer and is named as the one it is.
                text = found.get('xyz') or found.get('estimate') or ''
                rows = [line for line in text.splitlines()[2:] if line.strip()]
                if not found.get('ok'):
                    lines.append(str(found.get('status')
                                     or 'The chain did not run.'))
                    if rows:
                        _remember('the path and the climb')
                        _write_coords(xyz_document(
                            rows, 'Where the chain got to'))
                    _set_mol_status(*lines)
                    return
                said = _saddle.verdict(found.get('imaginary'),
                                       'What it reached')
                lines.append(f'{found["status"]} '
                             f'{found["seconds"]:.1f} s for the pair.')
                lines.extend(said['lines'])
                if rows:
                    # One step for the whole chain, which Undo takes back
                    # whole: two structures went in and one came out, and
                    # halfway back is not a place anyone asked for.
                    _remember('the path and the climb')
                    kept = xyz_document(
                        rows, f'From a path, optimised to {said["name"]}')
                    _write_coords(kept)
                    _note_the_saddle(kept, found.get('imaginary'))
                    lines.append('It is in the box; Undo takes it back.')
                    lines.extend(_the_mode_is_offered(found.get('imaginary')))
                    if said['first_order']:
                        lines.append(
                            'Refine it with OPTTS in the ORCA Builder at a '
                            'level worth quoting.')
                _set_mol_status(*lines)

            schedule_ui_update(_done)

        _start_background(_work, 'The walk and the climb after it',
                          guards={'chain_run': False})

    def _band_between(ends):
        """A nudged elastic band between the two ends, and the climb off it.

        The arbiter.  Everything else this press does is fast, and fast is
        exactly what makes it worth having something slower to check against:
        the chain above walks its own way between two ends with metadynamics
        and climbs the highest point of it, and when the answer that comes
        back is not believable there has to be a second opinion that shares no
        machinery with the first.  A band is that -- eight images relaxed onto
        the way between the two ends at once, rather than one structure walked
        along it -- and :func:`delfin.dashboard.saddle.neb_to_saddle` is where
        the measurements live.

        Two things are checked before ORCA is started, and both of them were
        measured failing the slow way: a band between two ends whose
        interpolation pulls a fragment off computes both ends perfectly and
        kills every image in between, and a band between two ends whose atoms
        are in different orders has no path at all.  Either one costs the
        whole timeout and returns nothing, so both are refused here in the
        second it takes to read a bond graph.

        It streams, because it is the one press here that runs for minutes:
        ORCA writes every band it accepts and then every step of the climb,
        both with energies, and the same frame channel the climb uses plays
        them.  The same press stops it and keeps what it had.
        """
        if state.get('path_run') or state.get('saddle_run') or state.get(
                'chain_run'):
            _set_mol_status('Something is already running on this structure; '
                            'let it finish or stop it first.')
            return
        method = str(submit_ff_dd.value)
        if method.lower() not in _saddle.SADDLE_METHODS:
            _set_mol_status(
                'A band is ORCA relaxing a chain of images on a semiempirical '
                'gradient, so choose GFN2, GFN1, GFN-FF or g-xTB. Anything '
                'with a basis set is a job for the ORCA Builder.')
            return
        state['band_run'] = True
        state['band_stop'] = False
        state['band_frame_run'] = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = state['band_frame_run']
        _note_the_run(state['band_frame_run'], 'band')
        _ensure_manip_bootstrap()
        schedule_ui_update(_install_gfn_frame_watcher)
        submit_saddle_btn.description = 'Stop'
        submit_saddle_btn.icon = 'stop'
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        label = _server_label(method)
        # How long to say it will be is a question about the machine and not
        # about the method: the images are computed at once, so the same
        # sixteen-atom band measured 39 s on eight processes and 272 on one.
        # So the sentence says what it is doing and that the press ends it,
        # and does not promise a time it cannot know.
        _set_mol_status(
            f'{label}: relaxing a band of {_saddle.NEB_IMAGES} images between '
            'the two ends, then climbing from the highest of them. This is '
            'the thorough route rather than the quick one; the press stops '
            'it.', spinner=True)

        def _work():
            sent = [0]

            def _watch(walked, energies):
                """The band while it relaxes, then the climb, as they arrive.

                Down the frame channel and not into the box, for the reason
                every other watcher here gives: a write to the box rebuilds
                the viewer, and a band is a hundred and seventy frames.
                """
                if state.get('gfn_run') != state.get('band_frame_run'):
                    return
                for n in range(sent[0], len(walked)):
                    schedule_ui_update(
                        lambda text=json.dumps({
                            'run': state.get('band_frame_run'),
                            'from': n,
                            'follow': 1,
                            'frames': [[round(float(v), 4)
                                        for v in walked[n]]],
                        }): setattr(submit_gfn_frame, 'value', text))
                sent[0] = len(walked)
                schedule_ui_update(
                    _set_mol_status,
                    f'The band is relaxing: {len(walked)} images accepted so '
                    'far.', spinner=True)

            found = _saddle.neb_to_saddle(
                ends[0], ends[1], method, charge=charge, uhf=uhf,
                solvent=wet, on_frame=_watch,
                # The method's own allowance.  g-xTB is measured at 716 s on
                # the sixteen-atom Diels-Alder, where the methods ORCA drives
                # itself are 39 -- every gradient is a separate process --
                # so one number for both would stop one of them short.
                timeout=_saddle.neb_seconds_for(method),
                should_stop=lambda: bool(state.get('band_stop')),
                on_stage=lambda said: schedule_ui_update(
                    _set_mol_status, said, spinner=True))

            def _done():
                state['band_run'] = False
                state['band_stop'] = False
                _name_the_saddle_press()
                lines = []
                text = found.get('xyz') or ''
                rows = [line for line in text.splitlines()[2:] if line.strip()]
                if not found.get('ok'):
                    lines.append(str(found.get('status')
                                     or 'The band did not run.'))
                    if rows:
                        _remember('the band')
                        _write_coords(xyz_document(
                            rows, 'Where the band got to'))
                    _set_mol_status(*lines)
                    return
                said = _saddle.verdict(found.get('imaginary'),
                                       'What it reached')
                lines.append(found['status'])
                if found.get('reaction') is not None:
                    lines[-1] += (f' The two ends are '
                                  f'{found["reaction"]:+.1f} kcal/mol apart '
                                  'along it.')
                lines.extend(said['lines'])
                if rows:
                    # One step for the whole band, which Undo takes back
                    # whole: two structures went in and one came out.
                    _remember('the band')
                    _write_coords(xyz_document(
                        rows, f'From a band, optimised to {said["name"]}'))
                    lines.append('It is in the box; Undo takes it back.')
                    if said['first_order']:
                        lines.append(
                            'Refine it with OPTTS in the ORCA Builder at a '
                            'level worth quoting.')
                _set_mol_status(*lines)

            schedule_ui_update(_done)

        _start_background(_work, 'The band between the two ends',
                          guards={'band_run': False})

    def _walk_the_path(ends, then_climb=False):
        """Walk between two ends and keep what is found.

        Its own answer, not a refinement of a scan's: a scan drives a
        coordinate somebody chose and this finds its own way between two
        structures.  When the two agree, that is worth more than either of
        them alone; when they do not, the difference is the interesting part
        and both numbers are said.

        *then_climb* hands what it estimated to the climb, which is the
        by-hand half of the same question: the walk says where the top is and
        the climb walks onto it at ten milliseconds a step, so it can be
        watched and dragged in the middle of.  The estimate is in the box
        before the climb starts, which is what the climb reads, so a climb
        that goes somewhere nobody wanted still leaves the walk's own answer
        one Undo away.
        """
        method = str(submit_ff_dd.value)
        if not ends:
            return
        if not _gfn.is_gfn_method(method):
            _set_mol_status('A path needs xtb: choose a GFN method.')
            return
        if state.get('path_run') or state.get('chain_run'):
            return
        state['path_run'] = True
        # The press keeps its name through this half.  A walk is one call to
        # xtb and cannot be interrupted part way, so a "Stop" on it would be a
        # button that does nothing -- what says the walk is running is the
        # ring and the line under the picture, which say it for every other
        # calculation here too.
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        wet = str(submit_gfn_solvent.value or '') or None
        model = _solv_model()
        label = _server_label(method)
        _set_mol_status(f'{label} is looking for a path between the two '
                        'ends...', spinner=True)

        def _work():
            found = _gfn.walk_the_path(
                ends[0], ends[1], method, charge=charge, uhf=uhf,
                solvent=wet, solvation_model=model)
            # And a Hessian on what it estimates, because "estimated
            # transition state" is a phrase and one imaginary frequency is a
            # fact.  Cheap: 0.6 s on sixteen atoms.
            state['path_shape'] = None
            state['path_depth'] = ''
            if found.get('ok') and found.get('ts'):
                shape = _gfn.optimize_with_gfn(
                    found['ts'], method, charge=charge, uhf=uhf, timeout=None,
                    solvent=wet, solvation_model=model, optimise=False,
                    free_energy=True,
                    thermo_kelvin=float(submit_temperature.value or 298.15))
                state['path_shape'] = shape.get('imaginary')
                # And whether the method can still answer where the barrier
                # is.  A closed shell describes two electrons in one orbital;
                # at the top of a bond-breaking barrier they are not in one,
                # and the frontier gap says so before the energy does.
                # Measured on a 2-butene: 5.28 eV at cis and 2.42 at the
                # twisted top, where GFN2 also invents a minimum 64 kcal/mol
                # above cis that real trans lies 1.5 below.
                start = _gfn.optimize_with_gfn(
                    ends[0], method, charge=charge, uhf=uhf, timeout=None,
                    solvent=wet, solvation_model=model, optimise=False)
                state['path_depth'] = _gfn.method_is_out_of_its_depth(
                    shape.get('gap'), start.get('gap'))

            def _done():
                state['path_run'] = False
                _name_the_saddle_press()
                if not found.get('ok'):
                    _set_mol_status(str(found.get('status')
                                        or 'The path finder did not run.'))
                    return
                T = float(submit_temperature.value or 298.15)
                needs = thermal_temperature(found['barrier'], _THERMAL_SECONDS)
                # How near the path came to the structure it was aimed at.
                # Without it a barrier is a number about some other reaction:
                # the finder always reports one, and only this says whether it
                # is the one that was asked for.
                near = found.get('rmsd')
                aimed = ('' if near is None else
                         f' It came within {near:.2f} A RMSD of the structure '
                         f'it was aiming at.')
                lines = [
                    f'{label} walked {found.get("points") or "a"} points in '
                    f'{found["seconds"]:.1f} s: {found["barrier"]:.1f} '
                    f'kcal/mol forward'
                    + (f', {found["back"]:.1f} back' if found.get('back')
                       is not None else '')
                    + (f', {found["reaction"]:.1f} for the reaction'
                       if found.get('reaction') is not None else '')
                    + f'.{aimed}',
                ]
                lines.append(
                    'No temperature under 5000 K crosses that within '
                    f'{_timescale_label()}.' if needs is None else
                    f'It wants about {needs:.0f} K '
                    f'({needs - 273.15:+.0f} C) within {_timescale_label()}; '
                    f'you have {thermal_ceiling(T, _THERMAL_SECONDS):.1f} '
                    f'kcal/mol at {T:g} K.')
                # Whether it is one.
                #
                # A path finder returns an estimated transition state whatever
                # happened, and "estimated" is doing a great deal of work in
                # that phrase.  What makes a structure a transition state is
                # one mode going the wrong way and no others, and that is a
                # Hessian -- 0.6 s on sixteen atoms, which is nothing beside
                # not knowing.  Measured on this one: a single imaginary
                # frequency at -131.4 cm-1, so it really is a saddle point at
                # this level of theory and worth handing to ORCA's OptTS.
                if state.get('path_depth'):
                    lines.append(state['path_depth'])
                shape = state.get('path_shape')
                if shape is not None:
                    # Said in the same words the saddle search says it in, and
                    # without its advice: an estimate is nobody's stationary
                    # point, and what to do about one is the next line rather
                    # than a displacement along a mode.
                    lines.extend(_said_modes(
                        shape, 'The structure it estimates', advise=False))
                    if shape.get('count') == 0:
                        lines.append(
                            'The path went over something this estimate is '
                            'not sitting on.')
                    if not then_climb:
                        lines.append(
                            'Set the box beside the press to "through ORCA" '
                            'to sharpen it in seconds, or to "by hand" to '
                            'walk onto it yourself -- measured on a '
                            'sixteen-atom case, 3.6 s for the path and 12 s '
                            'for the pair.')
                rows = [line for line in (found.get('ts') or '').splitlines()[2:]
                        if line.strip()]
                if rows:
                    # One step, which Undo takes back whole -- the same as a
                    # scan, because this replaces the structure just as much.
                    _remember('the transition state the path finder estimated')
                    _write_coords(xyz_document(
                        rows, 'Estimated transition state, from the path'))
                    lines.append(
                        'The structure it estimates for the transition state '
                        'is in the box; Undo takes it back.')
                _set_mol_status(*lines)
                # And on into the climb, from the estimate that is now in the
                # box.  One turn of the interface later, because the climb
                # reads the box and the write above has to have landed in it
                # first -- and because the walk's own answer belongs on screen
                # before anything starts walking away from it.
                if then_climb and rows:
                    state['path_said'] = tuple(lines)
                    schedule_ui_update(_climb_from_here)

            schedule_ui_update(_done)

        _start_background(_work, 'The path search',
                          guards={'path_run': False})

    def _scan_two_legs():
        """The last scan's two legs and what they say about each other.

        One call rather than four state keys, because the profile wants to be
        drawn and a picture of it is somebody else's part: this is the shape
        that part is written against, and the keys behind it can move without
        the drawing having to.

        Returns ``{'there', 'back', 'disagree', 'jumped'}``.

        * ``there`` and ``back`` are ``[(coordinate, kcal/mol)]`` in the order
          each leg was walked -- so ``back`` runs the other way along the
          coordinate -- and both are against the same zero, which is the first
          point of the walk out.  They belong on one pair of axes.
        * ``back`` is empty when no second leg was walked: the toggle was up,
          the walk was stopped, or it was a push.
        * ``disagree`` is ``{'at', 'gap', 'there', 'back', 'points'}`` or
          None -- the coordinate value where the two legs are furthest apart,
          which is the one point a drawing should mark.
        * ``jumped`` is None or a dict whose ``step`` indexes ``there``: the
          point the walk landed on, so the discontinuity is the segment from
          ``step - 1`` to ``step``.  ``named``, ``was`` and ``now`` say which
          internal coordinate went with it, when one can be named.
        """
        return {'there': list(state.get('scan_there') or ()),
                'back': list(state.get('scan_back') or ()),
                'disagree': state.get('scan_disagree'),
                'jumped': state.get('scan_jumped')}

    def _scan_free_is_an_estimate():
        """What a Hessian at a scan point had to say about being one.

        An RRHO free energy is a sum over frequencies, and it means something
        at a stationary point.  Two of the three places a scan takes one are
        stationary; the top of a barrier is not, and the Hessian says so
        itself by coming back with modes that go the wrong way.  Measured
        under GFN2 at the top of a Diels-Alder scan: one mode at -128 cm-1,
        with the gradient 67 times the threshold that would call the geometry
        converged.

        There is a published way of taking a free energy at a geometry like
        that -- xtb's biased single-point Hessian -- and it is measured in
        :func:`_free` not to apply here: it relaxes off the point it was asked
        about.  So the mode is reported rather than removed.  A number that
        says what it is worth is a usable number; the same number without that
        is the one somebody quotes.
        """
        shape = state.get('scan_free_shaky')
        if not shape:
            return ''
        many = int(shape.get('count') or 0)
        modes = [one for one in (shape.get('modes') or ()) if one < 0]
        worst = f' ({modes[0]:.0f} cm-1)' if modes else ''
        return (f' One of those Hessians came back with '
                f'{"a mode" if many == 1 else f"{many} modes"} going the '
                f'wrong way{worst}, which is the surface saying that point is '
                f'not a stationary point -- a barrier top is not one -- so '
                f'the free energies are an estimate rather than the '
                f'thermodynamics of two minima.')

    def _scan_can_be_quoted(kelvin):
        """Whether the walk's own barrier is a number to quote, and why not.

        Three sentences at most, and often none: a push says nothing here
        because it has no second leg to compare against, and a walk that
        agreed with itself says so in one clause and stops.

        The order is the order a reader needs it in.  First whether the two
        legs of the walk are the same curve, which is the whole question.
        Then, if a step of it was a fall rather than a step, where that was
        and *what moved* -- because a user told only that their scan jumped
        can do nothing about it, and a user told which coordinate slipped can
        arm that one as well and walk both together, which is what this
        editor's several-legs-at-once scan is for.

        Wording that holds for whatever is being computed: a slipped
        coordinate is named by its two atoms and nothing is assumed about what
        kind of reaction it belongs to.
        """
        said = ''
        gap = state.get('scan_disagree')
        apart = _gfn.a_rate_apart(kelvin)
        if gap is not None:
            if gap['gap'] <= apart:
                said += (
                    f' Walked back over the same {gap["points"]} points, the '
                    f'two legs agree to {gap["gap"]:.2f} kcal/mol -- inside '
                    f'the {apart:.2f} that a factor of ten in rate is worth '
                    f'at {float(kelvin):g} K -- so this profile is the path '
                    f'and not the direction it was walked.')
            else:
                said += (
                    f' Walked back over the same {gap["points"]} points, the '
                    f'two legs disagree by {gap["gap"]:.1f} kcal/mol at '
                    f'{gap["at"]:.3g}: {gap["there"]:+.1f} on the way out '
                    f'against {gap["back"]:+.1f} on the way back. A driven '
                    f'scan is a path only while the coordinates nobody is '
                    f'driving follow it continuously; where one slips, each '
                    f'direction misses the crossing on its own side, so the '
                    f'height above is where this walk went and not the '
                    f'barrier. Two ends and a saddle search will answer what '
                    f'the walk cannot.')
        elif state.get('scan_back_wanted') is False:
            # Only where the second leg could have run and was not asked for.
            # A push has no grid to retrace, and a walk that was stopped or
            # that collapsed has already been told why it ended -- offering it
            # a second leg there would be answering a question nobody is in a
            # position to ask.
            said += (' Nothing walked it back, so whether this profile '
                     'depends on the direction it was walked is not known. '
                     '"Walk it back" beside Run scan answers that, for '
                     'another leg of the same walk.')
        fell = state.get('scan_jumped')
        if fell is not None:
            said += (
                f' It jumped between {fell["from"]:.3g} and {fell["at"]:.3g}: '
                f'{fell["fell"]:+.1f} kcal/mol in one step, where the rest of '
                f'the path bends by {fell["scale"]:.2f}.')
            if fell.get('named'):
                said += (f' What went with it was {fell["named"]}, '
                         f'{fell["was"]:.2f} to {fell["now"]:.2f} '
                         f'({fell["moved"]:.2f} A) -- arm that as well and '
                         f'walk both together.')
        return said

    def _scan_verdict(path, steps):
        """What the walk found, and the temperature it would take.

        Said as a temperature rather than only as a number: a drag has to be
        told what it may do at the temperature it was given, but a scan has
        walked the whole path and can answer the question a chemist actually
        asks.  "+29 kcal/mol" means nothing until it is "and that wants 385 K
        within the hour".
        """
        if not path:
            return ('The scan walked nothing.',)
        top = max(path, key=lambda one: one[1])
        # Where the walk was left, which after coming back to a minimum is not
        # the last point it took.
        came = state.get('scan_came_back')
        # The barrier is the rise out of the minimum the path left, not the
        # height above whatever the first point happened to be -- see
        # :func:`_climbed`.  Both are said: the rise is what a temperature is
        # worked out from, and the height is where the path stands.
        _, rise = _climbed([one[1] for one in path])
        ends = came[1] if came else path[-1][1]
        T = float(submit_temperature.value or 298.15)
        # The free energies, where they were asked for, and then they are the
        # barrier -- not a second opinion printed beside it.  Eyring inverts a
        # free energy of activation, so a ceiling compared against anything
        # else is an approximation standing in for one; when the real thing is
        # there it is what the temperature is worked out from.  Read after the
        # temperature was taken from the electronic rise, the line quoted one
        # number and reasoned from the other.
        free = state.get('scan_free')
        if free is not None:
            rise = free[0]
        ceiling = thermal_ceiling(T, _THERMAL_SECONDS)
        needs = thermal_temperature(rise, _THERMAL_SECONDS)
        arrived = (f' It came back to the minimum it walked through, at '
                   f'{came[0]:.3g}.'
                   if state.get('scan_arrived') and came else '')
        if state.get('scan_stopped_out') and not state.get('scan_arrived'):
            arrived = (' You stopped it there, so the highest point is where '
                       'the walk was interrupted rather than a barrier.')
        crowded = state.get('scan_crowded')
        if crowded:
            arrived = (f' It stopped: two atoms came inside {crowded:.2f} of a '
                       f'bond length, which is no path at any temperature. The '
                       f'target is past the far side of a bond.')
        # A push prices its crossing with points of its own, so the count
        # can pass the number of steps that were asked for -- said as it is
        # rather than as "27 of 20".
        many = (f'{len(path)} points' if len(path) > steps
                else f'{len(path)} of {steps} points')
        # Said as one quantity or as two, never as two wearing one name: a
        # highest point in E beside a rise in G reads as a contradiction, and
        # a reader has no way to tell which is which.
        if free is None:
            first = (f'The scan walked {many}. Highest '
                     f'{top[1]:+.1f} kcal/mol at {top[0]:.3g}, a rise of '
                     f'{rise:.1f} out of the minimum before it, ending '
                     f'{ends:+.1f}.{arrived}')
        else:
            first = (f'The scan walked {many}. Highest at {top[0]:.3g}: '
                     f'{free[0]:+.1f} kcal/mol as a free energy at {T:g} K '
                     f'({top[1]:+.1f} electronic), ending {free[1]:+.1f} '
                     f'({ends:+.1f}).{arrived} The free energies are from '
                     f'three Hessians -- where it started, the top, and where '
                     f'it came to -- and they are what the temperature below '
                     f'is worked out from.{_scan_free_is_an_estimate()}')
        # What was asked for, and whether it happened.  Said before the
        # temperature, because it is the question: a walk given a verb was not
        # asked how high the path was, it was asked to make a bond.
        #
        # And when it did not happen, said as a statement about the chemistry
        # rather than as a failure.  The ramp ends at
        # :data:`gfn_optimize.PUSH_FORCE_TO`, which is more than twice
        # :data:`gfn_optimize.A_BOND_HOLDS` -- the force a bond holds against
        # -- so a bond that has not formed under it is one this method will
        # not form from here, and that is an answer.
        instructed = list(state.get('scan_instructed') or ())
        done = state.get('scan_carried_out')
        if instructed and done:
            step, force, said = done[0], done[1], (done[2] if len(done) > 2
                                                   else '')
            first += (f' The instruction was carried out at step {step}, '
                      f'under {force:.0f} kcal/mol/A'
                      + (f': it {said}.' if said else '.'))
        elif instructed:
            first += (f' The instruction was not carried out: '
                      f'{", ".join(instructed)} still {"do" if len(instructed) > 1 else "does"} '
                      f'not hold at {_gfn.PUSH_FORCE_TO:.0f} kcal/mol/A, '
                      'which is twice what a bond holds against. From this '
                      'structure and on this method, that is the answer '
                      'rather than a setting to raise.')
        if state.get('scan_depth'):
            first += ' ' + state['scan_depth']
        # And whether this profile is one path or two joined at a fall.  Said
        # here, beside the barrier it is about, rather than on a line of its
        # own after the temperature: a caveat that arrives after the number
        # has been read is a caveat nobody applied.
        first += _scan_can_be_quoted(T)
        if free is None and str(submit_scan_energy.value) == 'G':
            first += (' The free energies could not be taken, so these are '
                      'electronic.')
        # The temperature is said either way.  Said only when the path was
        # closed, the number a chemist came for was missing exactly when the
        # answer was good news -- and "it needs 150 K and you have 298" is
        # what makes an open path mean something.
        if needs is None:
            wants = ('No temperature under 5000 K crosses that within '
                     f'{_timescale_label()}, which is another way of saying '
                     'it does not happen.')
        else:
            wants = (f'It wants about {needs:.0f} K ({needs - 273.15:+.0f} C) '
                     f'to be crossed within {_timescale_label()}.')
        # What was handed back, when it is not where the walk ended.  The walk
        # is reported whole either way -- that is what a scan is for -- but a
        # structure the temperature cannot reach is not one to carry on from,
        # so the box has the last point it could and this says so.
        walled = state.get('scan_walled')
        held_back = ('' if walled is None else
                     f' Where it ended costs {walled:+.1f} kcal/mol against '
                     f'this structure, which is past the {ceiling:.1f} '
                     f'available at {T:g} K, so the box has the last point '
                     f'that was inside it.')
        # Which for an instruction is the whole of it: the box now holds a
        # structure where the bonds were not made, and the line has to say so
        # or it is a claim about a geometry nobody has.
        if walled is not None and done and not state.get(
                'scan_carried_out_kept'):
            held_back += (' So the structure you have is from before the '
                          'bonds changed -- the budget priced the change and '
                          'this temperature cannot pay for it. Raise the '
                          'temperature to see what it would take, or switch '
                          'the budget off to keep what the force reached.')
        # And what the walk has made possible, said where the walk is being
        # reported.  The two ends are an entry in a box rather than the two
        # buttons that used to appear, and a box that has gained an entry is
        # quiet -- so the sentence says what the toolbar has just done, and
        # :func:`_scan_left_two_ends` is what makes the sentence true rather
        # than a description of something the user would have to find.
        left = ('' if not state.get('scan_ends') else
                ' It left two ends, and the press now starts from them: one '
                'press walks its own way between the two and climbs what it '
                'finds, without the coordinate you chose. The box beside it '
                'says how far to go, and the one before it goes back to the '
                'structure on screen.')
        if rise <= ceiling:
            return (first,
                    f'{wants} You have {ceiling:.1f} kcal/mol at {T:g} K, so '
                    f'the whole path is open. {_thermal_wait(rise, T)}'
                    + held_back + left)
        return (first,
                f'{wants} At {T:g} K only {ceiling:.1f} kcal/mol is '
                f'available, so the path is closed there. '
                f'{_thermal_wait(rise, T)}' + held_back + left)

    def _scan_plot_drop():
        """Take the profile off the page.

        A picture is a claim about what was measured, and it must not outlive
        the thing it describes: a profile standing beside a structure the walk
        never visited is the same lie as a control that offers something
        impossible.  So it goes at the two moments that end it -- a run that
        starts drawing over the structure (see :func:`_note_the_run`), and a
        geometry arriving in the box that is not this scan's own (see
        :func:`_scan_plot_holds`).

        Cheap enough to call from either: nothing is drawn most of the time,
        and then this is one dictionary lookup.
        """
        if not state.get('scan_plot'):
            return
        state['scan_plot'] = None
        submit_scan_plot.value = ''
        submit_scan_plot.layout.display = 'none'

    def _scan_plot_holds(text):
        """Whether the profile still describes the geometry now in the box.

        Two things have to be true, and they are the two the rest of the
        editor already tells structures apart by.  The comment line has to be
        one the scan wrote -- ``Scanned``, and the landmarks the walk left in
        the history, so stepping back through them with Undo keeps the picture
        that names them -- and the element column has to be the molecule the
        walk was made on, the way :func:`_scan_ends_here` and the thermal
        anchor decide the same question.

        Everything else takes it away: a drag ends with its own comment, an
        optimisation with its own, Reset and a newly loaded structure with
        theirs, and Undo past the scan lands on the entry from before it.
        """
        lines = str(text or '').splitlines()
        comment = lines[1].strip().lower() if len(lines) > 1 else ''
        if not comment.startswith('scanned'):
            return False
        walked_on = state.get('scan_plot_of')
        return bool(walked_on) and _structure_fingerprint(text) == walked_on

    def _on_box_for_scan_plot(change):
        """The box has changed: is the profile still about what is in it?"""
        if state.get('scan_plot') and not _scan_plot_holds(change.get('new')):
            _scan_plot_drop()

    def _scan_profile_html(path, legs, pushing, began=None, kept=None):
        """Draw the walk that has just finished, and hand back the picture.

        Once, here, and not a point at a time.  Measured with the interpreter
        the dashboard runs on: 0.17 s to build the figure and encode it, and
        40 to 66 kB of PNG depending on how many points there are.  The walks
        it was measured against cost 6 to 10 s a point under GFN2 -- a
        twenty-point Diels-Alder approach over sixteen atoms took 144 s, a
        twenty-four-point torsion 238 s -- so drawing at every point would
        cost a few per cent of the time and push a megabyte of pictures nobody
        waits for down the channel the frame stream needs.  At the end it is
        one write, after the last call to xtb has returned.

        Nothing is computed again.  Every number here was in hand when the
        verdict was written: the path, the free energies where they were
        asked for, and where the walk was left.

        Built here, on the scan's own worker thread, and shown by
        :func:`_show_scan_profile` when the interface's turn comes.  The two
        are separate because the first import of matplotlib costs 0.3 s on an
        idle machine and a second on a loaded one, and a second spent inside
        an interface callback is a second in which the dashboard does not
        answer.  Nothing here touches a widget.

        Wrapped, because a picture must not be able to take away the answer.
        A scan that has walked for minutes reports its verdict whatever
        matplotlib does with it.
        """
        # A point of the walk is where it was, what it cost, and the
        # structure it held there -- the third is what the marked points
        # and Undo are made of, and it is not drawn.  Taken by index
        # rather than unpacked, so a point that grows a fourth thing
        # later does not stop the picture; the pairs made here are the
        # same shape as *began*, *kept* and the point it came back to,
        # which everything below compares against.
        points = [(one[0], one[1]) for one in path
                  if one[0] is not None and one[1] is not None]
        # Two points are a line between two numbers the sentence already
        # gives.  A picture is worth its row when there is a shape to see.
        if len(points) < 3:
            return None
        try:
            unit = 'A' if legs[0]['kind'] == 'distance' else 'deg'
            symbols = _leg_atoms_label(legs[0])
            top = max(points, key=lambda one: one[1])
            came = state.get('scan_came_back')
            ended = came if came else points[-1]
            if state.get('scan_crowded'):
                ended_label = 'where it stopped: two atoms too close'
            elif state.get('scan_arrived') and came:
                ended_label = 'the minimum it came back to'
            elif state.get('scan_stop'):
                ended_label = 'where you stopped it'
            else:
                ended_label = 'where it ended'
            # The free energies at the three places they were taken, on the
            # same axis: both are kcal/mol above where the walk started, so
            # this is one scale and not two dressed as one.
            free = state.get('scan_free')
            drawn_free = ()
            if free is not None and began is not None:
                drawn_free = ((began[0], 0.0, 'where it started'),
                              (top[0], free[0], 'the highest point'),
                              (ended[0], free[1], 'where it came to'))
            # The ceiling where the verdict compares against it: over the
            # minimum the rise is measured out of, or -- when the free
            # energies are what the temperature is worked out from -- over
            # where the walk started, which is their zero.
            T = float(submit_temperature.value or 298.15)
            ceiling = thermal_ceiling(T, _THERMAL_SECONDS)
            _, rise = _climbed([one[1] for one in points])
            floor = 0.0 if free is not None else top[1] - rise
            said = []
            if len(legs) > 1:
                said.append('Walked together with '
                            + ', '.join(_describe_leg(one)
                                        for one in legs[1:]) + '.')
            said.append('Undo steps back through the marked points.')
            return _scan_profile.profile_html(
                points,
                x_label=f'{symbols} ({unit})',
                y_label='kcal/mol above the start',
                title=(f'{symbols}, pushed with a force' if pushing
                       else f'{symbols}, walked')
                + (f' and {len(legs) - 1} more together' if len(legs) > 1
                   else ''),
                started=began, top=top, ended=ended,
                ended_label=ended_label,
                free=drawn_free,
                ceiling=floor + ceiling,
                ceiling_label=f'{ceiling:.1f} kcal/mol at {T:g} K',
                kept=kept if state.get('scan_walled') is not None else None,
                note=' '.join(said))
        except Exception as exc:                  # a picture, not the answer
            record('note', v=f'the scan profile could not be drawn: {exc}')
            return None

    def _show_scan_profile(drawn):
        """Put a picture that has been drawn on the page, and say what it is of.

        The fingerprint is taken here rather than where the figure was built:
        the structure the picture is a claim about is the one now in the box,
        and the box is written a few lines before this.
        """
        if not drawn:
            return
        state['scan_plot_of'] = _structure_fingerprint(_current_xyz() or '')
        state['scan_plot'] = True
        submit_scan_plot.value = drawn
        submit_scan_plot.layout.display = ''

    def on_submit_hold(_button=None):
        """Hold the value the selection describes while the field runs."""
        indices = list(state.get('picked') or [])
        kind = _CONSTRAINT_KINDS.get(len(indices))
        if not kind:
            _set_mol_status('Pick 2, 3 or 4 atoms before holding a value.')
            return
        entry = {
            'kind': kind,
            'atoms': indices,
            'value': float(submit_internal_value.value),
            # 'pull' negotiates with the chemistry and settles at a compromise;
            # 'fix' is restored after every relaxation step, so the value is met
            # exactly and the rest of the molecule arranges itself around it.
            'mode': submit_hold_mode.value,
            # And which structure the numbers above belong to.  They are atom
            # indices and nothing more, so they mean something about every
            # structure with that many atoms -- which is how a C-C held at
            # 1.700 A on a cyclobutane reached an ORCA input for a benzene as
            # "{ B 0 1 1.7000 C }", both being twelve atoms.  The element
            # column is what tells the two apart, and it is written down here
            # while the structure it was set on is the one on screen.
            'structure': _structure_fingerprint(_current_xyz() or ''),
        }
        _remember(f'holding {_describe_constraint(entry)}')
        held = list(state.get('constraints') or [])
        held = [c for c in held if c['atoms'] != indices]
        held.append(entry)
        state['constraints'] = held
        _refresh_constraints()
        _set_mol_status(f'Holding {_describe_constraint(entry)}.')
        _enable_live_forcefield()
        # Under GFN nothing is running between drags, so the change is what
        # starts it: a value set while the switch is down used to sit there
        # until Optimise was pressed.
        _arm_gfn_takeup(f'Holding {_describe_constraint(entry)}')
        # A fresh set for the next one: several values can then be held at
        # once, which is the whole point of a list.
        _clear_selection()

    def on_submit_hold_mode(change):
        """Retune the selected constraint, so a mode can be changed without
        having to select the atoms and set the value again."""
        if change.get('name') != 'value':
            return
        if state.get('hold_mode_quiet'):
            return
        position, entry = _selected_constraint()
        if entry is None:
            return
        # A pull settles at a compromise and a fix is met exactly, so this
        # moves the structure -- and it armed the relaxation that moves it
        # without leaving anywhere to come back to.
        _remember(f'holding {_describe_constraint(entry)} as '
                  f'{submit_hold_mode.value}')
        held = list(state.get('constraints') or [])
        held[position] = dict(entry, mode=submit_hold_mode.value)
        state['constraints'] = held
        _refresh_constraints()
        _enable_live_forcefield()
        _arm_gfn_takeup(f'Holding {_describe_constraint(held[position])}')

    def on_submit_reset(_button=None):
        """Back to the structure that was loaded, with nothing set on it.

        Editing in the viewer is a one-way street otherwise: undo takes back
        one step at a time, and a structure that has been pulled apart over
        twenty of them has no way home short of pasting the coordinates again.

        And it is itself one more step, so Undo takes it back. It used to be
        the one action in the editor with no way back at all: an hour of work
        went at a press, and the history was cut to its first entry so that
        Undo -- which had been the way home from every other mistake -- landed
        on the structure Reset had just gone to and reported that there was
        nothing more to take back.
        """
        # The first thing in the history is the structure as it arrived. It
        # was a second remembered place before, set in the render path, and
        # the two could disagree -- which is how Reset came back to something
        # that was not the beginning.
        history = list(state.get('history') or [])
        pristine = (history[0].get('coords') if history
                    else state.get('pristine_coords'))
        if not pristine:
            _set_mol_status('Nothing to go back to yet.')
            return
        # Everything as it stands, before any of it is cleared. The write at
        # the bottom goes through the host's own "this is a structure I have
        # not seen" path -- which is what re-renders and re-perceives, and
        # what puts the toolbar back to its defaults -- and that path throws
        # the history away and seeds a new one. So it is taken back out
        # afterwards rather than being left to survive a write that is meant
        # to clear everything else.
        _remember('the reset')
        kept = list(state.get('history') or [])
        aside = _stop_what_is_running()
        state['constraints'] = []
        state['scan_legs'] = []
        state['bond_edits'] = {}
        state['hand_bonds'] = {}
        state['hyb_overrides'] = {}
        state['structure_undo'] = []
        state['poly_applied'] = None
        state['poly_metal'] = None
        state['poly_assignment'] = None
        state['poly_arrangements'] = []
        state['poly_arrangement_index'] = 0
        submit_poly_turn_btn.layout.display = 'none'
        submit_poly_turn_btn.disabled = True
        state['poly_quiet'] = True
        try:
            submit_poly_dd.value = ''
        except Exception:
            pass
        finally:
            state['poly_quiet'] = False
        state['hold_mode_quiet'] = True
        try:
            submit_hold_mode.value = 'pull'
        except Exception:
            pass
        finally:
            state['hold_mode_quiet'] = False
        submit_internal_value.value = 0.0
        _clear_selection()
        # Writing the coordinates is what re-renders and re-perceives; it also
        # clears everything above a second time, which is the point.
        if coords_widget.value == pristine:
            # Same text, so the write would be a no-op and nothing would be
            # redrawn -- the whole reason the viewer looks destroyed is that
            # the *coordinates* changed underneath it.
            update_view()
        else:
            coords_widget.value = pristine
        # The way back, put in after the write that cleared it. The first
        # entry of what is kept is the structure as it arrived, which is where
        # this has just gone; the last is everything as it stood a moment ago,
        # which is what Undo now returns to.
        state['history'] = kept
        state.pop('history_seed_pending', None)
        _refresh_constraints()
        _set_mol_status('Back to the structure as it was loaded. '
                        'Undo brings back what was here.' + aside)

    def on_submit_constraint_del(_button=None):
        key = submit_constraint_dd.value or ''
        if not key:
            return
        # Letting go of a held value lets the structure move again, and
        # setting one is a step -- so releasing one has to be, or Undo could
        # take back the holding and never the letting go.
        _remember('releasing the polyhedron' if key == 'poly'
                  else 'letting go of a held value')
        if key == 'poly':
            state['poly_applied'] = None
            state['poly_metal'] = None
            state['poly_assignment'] = None
            state['poly_arrangements'] = []
            state['poly_arrangement_index'] = 0
            submit_poly_turn_btn.layout.display = 'none'
            submit_poly_turn_btn.disabled = True
            state['poly_quiet'] = True
            try:
                submit_poly_dd.value = ''
            except Exception:
                pass
            finally:
                state['poly_quiet'] = False
        elif key.startswith('c'):
            held = list(state.get('constraints') or [])
            position = int(key[1:])
            if 0 <= position < len(held):
                held.pop(position)
            state['constraints'] = held
        _refresh_constraints()
        _enable_live_forcefield()

    def on_submit_poly_changed(change):
        if change.get('name') != 'value' or state.get('poly_quiet'):
            return
        # Putting a polyhedron on pulls the donors onto its vertices, which is
        # as much a change of the structure as dragging them there by hand.
        _remember(f'the {submit_poly_dd.label} polyhedron'
                  if submit_poly_dd.value else 'releasing the polyhedron')
        state['poly_applied'] = submit_poly_dd.value or None
        state['poly_assignment'] = None
        if state['poly_applied']:
            state['poly_metal'] = state.get('poly_offer_metal')
        else:
            state['poly_metal'] = None
        if state['poly_applied'] and state.get('poly_metal') is not None:
            try:
                from .molecule_forcefield import polyhedron_assignment
                state['poly_assignment'] = polyhedron_assignment(
                    state['perceived'], state['poly_metal'], state['poly_applied'],
                )
            except Exception:
                state['poly_assignment'] = None
        state['poly_arrangements'] = []
        state['poly_arrangement_index'] = 0
        _refresh_poly_turn()
        _refresh_constraints()
        # Re-assigning the parameters is what makes the pull start; with the
        # field running the complex visibly moves into the polyhedron.
        if submit_relax_btn.value or state.get('ff_bootstrap_done'):
            _enable_live_forcefield()
        # The metal has done its job once the polyhedron is on, exactly as
        # after Set or Hold: leaving it picked meant the next atom clicked
        # joined it, and the highlight sphere read as though something were
        # still waiting to be chosen. What is held stays in the list below.
        if state['poly_applied']:
            _set_mol_status(
                f'{submit_poly_dd.label}: the donors are pulled onto it.'
            )
        else:
            _set_mol_status('Polyhedron released.')
        _clear_selection()

    def on_submit_dyn_bonds(change):
        """Whether the drawn lines follow the distances."""
        if change.get('name') != 'value':
            return
        active = bool(submit_dyn_bonds_btn.value)
        submit_dyn_bonds_btn.button_style = 'info' if active else ''
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setDynamicBonds('
            f'{json.dumps(submit_scope_id)},{"true" if active else "false"});'
        )
        _set_mol_status(
            'The lines follow the distances now. Only the picture: what the '
            'calculation holds together is Bond and Unbond\'s business.'
            if active else
            'The lines keep the bonds the structure was drawn with.')

    def on_submit_auto_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_auto_btn.value)
        submit_auto_btn.button_style = 'info' if active else ''
        if not _server_method():
            # Nothing on the server to run: the browser's own field has no
            # minimisation to carry on with, so saying it would carry one on
            # would be a promise about a thing that is not there.
            return
        label = _server_label(submit_ff_dd.value)
        _set_mol_status(
            f'Letting go of an atom now runs {label} down to a minimum.'
            if active else
            'Letting go leaves the structure where you put it. Press '
            'Optimize when you want it taken down to a minimum.')

    def on_submit_settle_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_settle_btn.value)
        submit_settle_btn.button_style = 'info' if active else ''
        # Told to the browser either way, and told first.  Under a server
        # method this switch is gone and its answer is always no; returning
        # here before saying so left the page settling on release with a field
        # installed under some earlier method -- a structure relaxing when the
        # user had switched everything they could see to off.
        _ensure_manip_bootstrap()
        settling = active and not _server_method()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setSettleOnRelease('
            f'{json.dumps(submit_scope_id)},{"true" if settling else "false"});'
        )
        if _server_method():
            return
        _set_mol_status(
            'Letting go of an atom now lets the structure relax around it.'
            if active else
            'Atoms will stay exactly where you put them.')

    def _push_play_speed():
        """Tell the page how fast to walk the path.

        Frames a second on the slider, milliseconds a frame on the page: the
        user thinks in speed and the player counts in delay.
        """
        # The top of the slider is not a speed, it is "keep up": zero, which
        # the page reads as "take everything that has arrived, this frame".
        # This is the quick way and it reaches a page with nothing running;
        # the frames carry the same number, because this write races the
        # start-up script that builds the player and can lose.
        pace = _play_pace()
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinGfnPlay&&window.__delfinGfnPlay['
            + json.dumps(submit_scope_id) + ']){'
            'window.__delfinGfnPlay[' + json.dumps(submit_scope_id)
            + f'].pace={pace};' + '}'
        )

    def on_submit_play_speed(change):
        if change.get('name') != 'value':
            return
        _push_play_speed()
        asked = int(submit_play_speed.value)
        _set_mol_status(
            'The picture keeps up with the calculation: every frame is drawn '
            'as soon as it arrives.'
            if asked >= int(submit_play_speed.max) else
            f'The optimisation is drawn at {asked} frame(s) a second. Slower '
            'lets the calculation run ahead of the picture -- take hold of an '
            'atom and the frame you are looking at is the one that is kept. '
            'Dragging always keeps up, whatever this says.')

    def on_submit_thermal(change):
        """Switching the budget on anchors it; switching it off forgets."""
        if change.get('name') != 'value':
            return
        active = bool(submit_thermal_btn.value)
        submit_thermal_btn.button_style = 'info' if active else ''
        # The page's hand has a ceiling exactly while the budget does.
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setPullStrength('
            f'{json.dumps(submit_scope_id)},{_hand_share()},'
            f'{json.dumps(_pull_most() if active else None)});'
        )
        _refresh_hand_controls()
        if not active:
            state['thermal_e0'] = None
            state['thermal_for'] = None
            # And the page is told, or the last drag's leash outlives the
            # switch that made it: the marks stay, the reach stays, and the
            # next drag is held back by a budget that is no longer on.
            _clear_thermal_wall()
            _set_mol_status('The thermal budget is off. Drags are unmeasured '
                            'again.')
            return
        # A budget prices the relaxation that runs while the hand is down, so
        # that relaxation has to be running.  The page only reports a drag
        # while this switch is down -- without it nothing is calculated, the
        # ceiling has nothing to compare against, and every drag goes through
        # whatever the temperature says.  A switch that turns a limit off
        # without saying so is the defect this whole wall exists to close, so
        # the two go on together and the line says that they did.
        follow_too = ''
        if _server_method() and not submit_relax_btn.value:
            submit_relax_btn.value = True
            follow_too = (' The structure follows the hand while you drag, '
                          'because that relaxation is what the price is read '
                          'from.')
        _set_thermal_anchor(note_after=follow_too)

    def on_submit_scan_whole(change):
        """Only the light on the button; what it means is read where it is
        used."""
        if change.get('name') != 'value':
            return
        submit_scan_whole.button_style = (
            'info' if submit_scan_whole.value else '')

    def on_submit_scan_back(change):
        """The same, for the return leg."""
        if change.get('name') != 'value':
            return
        submit_scan_back.button_style = (
            'info' if submit_scan_back.value else '')

    def on_submit_scan_how(change):
        """A change of mode changes what is on the row.

        The return leg is walking's, and the row is redrawn rather than left
        showing a press that would do nothing under the mode now chosen.
        """
        if change.get('name') != 'value':
            return
        _refresh_scan()

    def _forget_topology_refusals(_change=None):
        """A new grab is a new question, so the count starts again."""
        state['topology_refused'] = 0

    def on_submit_topology(change):
        """Keep bonds, switched on or off.

        The graph is taken from the structure as it stands now, so what is
        kept is what the user is looking at when they ask for it -- not
        whatever was perceived when the molecule was loaded, which may be
        several edits ago.
        """
        if change.get('name') != 'value':
            return
        on = bool(submit_topology_btn.value)
        submit_topology_btn.button_style = 'info' if on else ''
        for key in ('topology_graph', 'topology_for', 'topology_good'):
            state.pop(key, None)
        _forget_topology_refusals()
        if not on:
            _set_mol_status('Bonds are free to make and break again.')
            return
        xyz = _current_xyz()
        if not xyz:
            return
        state['topology_graph'] = _gfn.bond_graph(xyz)
        state['topology_for'] = _structure_fingerprint(xyz)
        state['topology_good'] = xyz
        method = str(submit_ff_dd.value)
        # Where it already holds, say so rather than let it look like it is
        # doing something.  GFN-FF reads its bonding once and keeps it, so
        # under it the wall has nothing left to refuse.
        #
        # It used to be said for everything that was not GFN2 and its
        # relatives, which put it under the PM methods as well -- and that is
        # not true of them: MOPAC decides the bonding from the electrons like
        # any other semiempirical method, and the wall reads what came back
        # and takes the step away exactly as it does under GFN2.  The sentence
        # told the user the switch they had just pressed does nothing, while
        # it was working.
        already = ('' if method != 'gfnff'
                   else f' {_server_label(method)} already keeps its bonding, '
                        'so this changes nothing under it.')
        _set_mol_status(
            f'Keeping the {len(state["topology_graph"])} bonds this structure '
            f'has. A drag that would make or break one is taken back.'
            + already)

    def on_submit_thermal_anchor(_button=None):
        # Never relaxes: this is the user saying "from here", and the geometry
        # they are pointing at is the one they mean -- an intermediate they
        # have just reached, not the minimum near it.
        _set_thermal_anchor(relax=False, note='Measuring from this structure')

    def on_submit_thermal_terms(change):
        """T or the timescale moved: the anchor stands, the ceiling moves."""
        if change.get('name') != 'value' or not submit_thermal_btn.value:
            return
        _, ceiling = _thermal_budget()
        spent = _thermal_note(state.get('thermal_now'))
        # What the ceiling is worth as a slope, which is the part of it a
        # drag can be felt against: a barrier of that height climbed over an
        # angstrom is that steep, and the hand walks up anything shallower
        # than itself.  Said beside the hand's own setting rather than as a
        # cap on it -- the temperature limits the energy and nothing else.
        hand = _pull_force()
        steep = _gfn.push_force_for(float(ceiling)) if ceiling else None
        told = ('' if hand is None or steep is None else
                f' That is a slope of about {steep:.0f} kcal/mol/A over an '
                f'angstrom; the hand is set to {hand:.0f}, '
                f'{hand / _gfn.A_BOND_HOLDS:.2f} of what a bond holds.')
        _set_mol_status(
            f'At {float(submit_temperature.value):g} K this structure has '
            f'{ceiling:.1f} kcal/mol to spend within {_timescale_label()}.'
            + told + (f' {spent}' if spent else ''),
            # On its own row here, and only here: this is a press of the
            # temperature box, not a drag, so nothing is being aimed at while
            # the row appears.  The drag's own line never gains a second row
            # -- it stands above the viewer.
            _THERMAL_QUANTITY_SHORT)

    def on_submit_sens_changed(change):
        if change.get('name') != 'value':
            return
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setDragSensitivity('
            f'{json.dumps(submit_scope_id)},{float(submit_sens_slider.value)});'
        )

    def _submit_viewer_js():
        """The viewer this tab is showing, as an expression for the browser."""
        return '(window._submitMolViewerByScope||{})[%s]' % json.dumps(
            submit_scope_id)

    def _charges_here():
        """The charges of the last answer, if they still describe this one.

        Kept against the element column rather than against the coordinates:
        a charge is a property of the geometry and changes as it is dragged,
        so holding out for an exact match would show nothing at all during
        the one gesture the numbers are most worth watching. What it must
        never do is put one structure's charges on another structure's atoms,
        and the element column is what the rest of this editor tells those
        apart by -- the perception, the GFN-FF topology and the thermal anchor
        are all keyed on it, for the same reason.
        """
        charges = state.get('atom_charges')
        if not charges:
            return None
        if state.get('atom_charges_for') != _structure_fingerprint(
                _current_xyz() or ''):
            return None
        return charges

    def _label_texts_now():
        """What the labels are to say, or None for the atom numbers.

        An empty list where the charges are wanted and there are none yet,
        which draws nothing at all. Falling through to the numbers would have
        the box saying "charge" over a structure labelled with indices, and
        the two are told apart only by which of them happens to look like a
        charge -- atom 0 with a charge of 0 is the same sprite either way.
        Nothing on the atoms and a line underneath saying what would put
        something there is the honest state.
        """
        if str(submit_label_what.value) != 'charge':
            return None
        return atom_charge_texts(_charges_here()) or []

    def _repaint_labels(force=True):
        """Put the labels back on the atoms with whatever they say now.

        Nothing is recomputed and the molecule is not re-rendered: the sprites
        are rebuilt over the model the browser already holds.  A no-op with
        the labels off, so a drag under the ordinary settings pays nothing for
        this at all.

        *force* false is the drag's way in, and it is rate-limited on purpose.
        This travels by ``run_js``, which clears its output before displaying
        -- the reason the frames and the thermal wall each have a widget
        field of their own -- so a message sent twenty times a second is a
        message mostly overwritten before the page has drawn it.  Labels
        survive that where a frame would not, because each one is the whole
        state rather than one step of a sequence: the last to land is the
        right one and the ones lost in between were about geometries already
        gone.  What is left is the cost, and a 7 KB script at the rate GFN-FF
        answers a small molecule -- measured at 125 to 900 ms an answer on a
        borazane, and under 40 ms on the smallest -- is crowding a channel the
        drag has better uses for.

        Four times a second, which is faster than a charge visibly changes:
        measured on an ammonia borane pulled apart, the nitrogen went -0.25,
        -0.26, -0.27, -0.29, -0.33, -0.40 over six answers.
        """
        if not submit_labels_btn.value:
            return
        now = time.perf_counter()
        if not force:
            last = float(state.get('labels_drawn_at') or 0.0)
            if now - last < _LABEL_REPAINT_INTERVAL:
                return
        state['labels_drawn_at'] = now
        _run_manip_js(
            show_atom_numbers_js(
                var=_submit_viewer_js(), on=True,
                scale=scale_for_px(submit_label_size.value),
                texts=_label_texts_now())
        )

    def _remember_charges(outcome):
        """Keep what an answer computed, for the labels and for nothing else.

        Every server answer goes through here, including the ones a drag makes
        several times a second. It is a list assignment: the charges were read
        off a file the engine had already written, and this is where they
        stop being thrown away.
        """
        if not isinstance(outcome, dict):
            return
        charges = outcome.get('charges')
        if not charges:
            return
        # Named by the structure they were computed for and by the method
        # that computed them. Both matter: an answer that arrives after the
        # molecule has been replaced would otherwise draw one structure's
        # charges on another's atoms, and four methods give four different
        # numbers for the same atom.
        state['atom_charges'] = list(charges)
        state['atom_charges_method'] = str(outcome.get('method') or '')
        state['atom_charges_for'] = _structure_fingerprint(
            str(outcome.get('xyz') or ''))

    def _refresh_label_what():
        """Which labels this method can actually draw.

        The charge entry appears under the engines that compute one and is
        absent under the two that run in the browser, because there it would
        be a word for a quantity that does not exist. Absent rather than
        greyed: a control that is there and refuses is a question the user has
        to ask before finding out, and the answer never changes while the
        method does not.
        """
        offers = _server_method()
        entries = [('number', 'number')]
        if offers:
            entries.append(('charge', 'charge'))
        was = str(submit_label_what.value)
        state['label_what_quiet'] = True
        try:
            submit_label_what.options = entries
            submit_label_what.value = (
                was if was in {code for _name, code in entries} else 'number')
        finally:
            state['label_what_quiet'] = False
        submit_label_what.layout.display = (
            '' if submit_labels_btn.value else 'none')

    def on_submit_labels_toggle(change):
        """Labels on or off, and nothing else.

        The molecule is emphatically not rendered again. Rebuilding it from
        the coordinates is how a structure loses what only the browser knows:
        the bonds as they were perceived, and the ones made or broken by hand.
        Switching the numbers off used to do exactly that, so a molecule that
        had been optimised came back with bonds missing. The labels are a
        layer of sprites over the model -- they can be added and taken away
        without the model hearing about it.
        """
        if change.get('name') != 'value':
            return
        on = bool(submit_labels_btn.value)
        submit_label_size.layout.display = '' if on else 'none'
        submit_label_what.layout.display = '' if on else 'none'
        submit_labels_btn.button_style = 'info' if on else ''
        _run_manip_js(
            show_atom_numbers_js(
                var=_submit_viewer_js(), on=on,
                scale=scale_for_px(submit_label_size.value),
                texts=_label_texts_now())
        )
        if on:
            _say_what_the_labels_say()

    def _say_what_the_labels_say():
        """One line about charges, and only where it is worth a line.

        Nothing at all for the numbers, which say what they are. For the
        charges: which method's definition they are, because four methods give
        four different numbers for the same atom -- or, where no answer has
        been made yet, that there is nothing to draw and what would produce
        one. An empty structure with the word "charge" chosen and no
        explanation is the tool looking broken.
        """
        if str(submit_label_what.value) != 'charge':
            return
        charges = _charges_here()
        if not charges:
            _set_mol_status(
                'The labels will say the partial charges as soon as there is '
                'an answer to read them from -- press Optimise, or drag. '
                'Nothing extra is run for them: every server method writes '
                'them out already.')
            return
        method = str(state.get('atom_charges_method') or '')
        named = _server_label(method) if method else 'the last run'
        _set_mol_status(
            f'The labels are {named} partial charges, at the geometry of the '
            'last answer. Each method has its own definition of one, so they '
            'are read against each other within a structure and not across '
            'methods.')

    def on_submit_label_what(change):
        """Numbers or charges, on the labels that are already there."""
        if change.get('name') != 'value' or state.get('label_what_quiet'):
            return
        _repaint_labels()
        _say_what_the_labels_say()

    def on_submit_label_size(change):
        """Resize what is there in the viewer that is already there.

        Nothing is re-rendered: the browser rescales the label sprites it
        holds, so the size changes as the dropdown closes.
        """
        if change.get('name') != 'value':
            return
        _run_manip_js(
            atom_numbers_js()
            + 'window.__delfinAtomNumbers.setScale(%s,%.4f);'
            % (_submit_viewer_js(),
               scale_for_px(submit_label_size.value))
        )

    def on_submit_scan_way(change):
        """A direction, or an end the user gives, and the field that follows.

        Filled in with what the coordinate measures now the first time an end
        is asked for, so the field opens on a number that means something
        instead of on a zero that has to be guessed at.
        """
        if change.get('name') != 'value':
            return
        if str(submit_scan_way.value) == 'to' and not float(submit_scan_to.value):
            try:
                submit_scan_to.value = float(submit_internal_value.value)
            except (TypeError, ValueError):
                pass
        _refresh_scan()

    def _refresh_hand_controls():
        """Show only what the hand in use has something to act on.

        The hand's rule, kept in one place so that it and the method's rule --
        :func:`_refresh_method_controls`, which is about a different thing --
        compose rather than overwrite one another.  Both hands are wanted and
        neither is the lesser: placing an atom exactly where it is meant to go
        is what building a structure is, and no force can do it, because the
        whole point of a force is that the chemistry gets a say.

        What is on the list:

        * The pull slider.  It is the size of a force, and a placement has no
          force in it.  The page is told the same number through
          setPullStrength and its own field skips the pull entirely at a share
          of zero, which is what a placing hand sends.

        * The thermal budget, its temperature, its relax-first and its anchor.
          A budget has to price what is kept, and under a placing hand what is
          kept is not exactly what was priced: the answer comes back a little
          away from where the cursor asked and is laid back onto the hand,
          which is a second geometry.  Measured, +16.8 kcal/mol priced against
          +18.2 kept -- bounded by the same 0.25 A the loose hold is judged by
          and not by zero -- and whether a drag was stopped in time came down
          to a race between xtb at 70-170 ms and a page reporting every 16-120
          ms.  Under a pull that residual is exactly nothing: the geometry
          that was priced is the geometry that is kept, laid back only on the
          atoms the hand never touched.  A budget that cannot be made exact
          under one hand should not claim to be a budget there, so it goes
          away with that hand instead of sitting there being approximately
          true.

        What is deliberately not on the list, because taking away something
        that works is as bad as leaving something that does not:

        * Strength.  How many steps the browser's own field takes per
          animation frame -- the relaxation, not the hand.  The rest of the
          structure settles round a placed atom exactly as it settles round a
          pulled one.

        * Sensitivity.  How far the structure moves for how far the mouse
          moves, which is a property of dragging and not of what dragging
          does.

        * Keep bonds.  It reads the bonding off what came back and takes the
          step back if it changed; a placing hand is if anything the one that
          needs it more, because a rigid hand can put an atom anywhere at all.

        * The scan and the path finder.  A scan drives its own ramp of forces
          from PUSH_FORCE_FROM to PUSH_FORCE_TO and never reads the hand's
          slider, so it walks the same path whichever hand is chosen.  It does
          read the temperature box, for its own free energies and its verdict;
          that box is hidden here, and it was already hidden whenever the
          budget was switched off, so the scan is no worse off than it has
          always been and keeps whatever value it was set to.
        """
        pulling = _hand_pulls()
        submit_pull_slider.layout.display = '' if pulling else 'none'
        # Two conditions, about two different things.
        #
        # The method: the budget differences energies against an anchor, and
        # only xtb's follow hands back one that can be differenced -- MOPAC
        # reports a heat of formation, which is not the same quantity.  A
        # method that cannot price takes the switch down with it, because the
        # anchor it holds is an energy of the method that measured it and
        # means nothing under another.
        #
        # The hand: it does not take the switch down.  A hand is changed
        # constantly while working and losing the anchor each time would mean
        # measuring it again each time, so the setting survives out of sight
        # and comes back with the pull.  Nothing reads it in between -- see
        # the pricing in _gfn_follow_step, which asks _thermal_live.
        xtb = _gfn.is_gfn_method(submit_ff_dd.value)
        if not xtb and submit_thermal_btn.value:
            submit_thermal_btn.value = False
        shown = bool(xtb and pulling)
        submit_thermal_btn.layout.display = '' if shown else 'none'
        for widget in (submit_temperature, submit_thermal_relax,
                       submit_thermal_anchor_btn):
            widget.layout.display = ('' if (shown and submit_thermal_btn.value)
                                     else 'none')

    def on_submit_hand_changed(change):
        """Moving or pulling, and what each of them brings with it."""
        if change.get('name') != 'value':
            return
        pulling = _hand_pulls()
        _refresh_hand_controls()
        budget = bool(submit_thermal_btn.value)
        if state.get('hand_quiet'):
            # The method took the hand away or gave it back, which is the
            # toolbar following the choice that was just made.  The page still
            # has to be told below; the line belongs to that choice and says
            # so already, and a second sentence about the hand would land on
            # top of it.
            on_submit_pull_changed({'name': 'value'})
            return
        _set_mol_status(
            'Dragging pulls: the atom follows as far as the chemistry allows.'
            + (' The thermal budget measures it again.' if budget else '')
            if pulling else
            'Dragging moves: the atom goes where you put it and the rest '
            'settles around it.'
            + (' The thermal budget does not act on a placing hand -- what is '
               'kept there is not exactly what was priced -- so it is out of '
               'the way until you pull again.' if budget else ''))
        on_submit_pull_changed({'name': 'value'})

    def on_submit_pull_changed(change):
        if change.get('name') != 'value':
            return
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setPullStrength('
            f'{json.dumps(submit_scope_id)},{_hand_share()},'
            f'{json.dumps(_pull_most())});'
        )

    def on_submit_strength_changed(change):
        if change.get('name') != 'value':
            return
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setOptimizerStrength('
            f'{json.dumps(submit_scope_id)},{int(submit_strength_slider.value)});'
        )

    def _solv_model():
        """The continuum to run with, or '' when the method has none.

        Read from the dropdown, but never trusted past the method: the two are
        set separately and a stale model is how a run ends up asking for
        ddCOSMO under GFN-FF, which does not fail -- it returns a destroyed
        structure.
        """
        offered = _solvents.models_for(submit_ff_dd.value)
        if not offered:
            return ''
        chosen = str(submit_gfn_solv_model.value or '')
        return chosen if chosen in offered else offered[0]

    def _refresh_solvation_controls():
        """Put the two solvent controls in step with the chosen method.

        Both lists are properties of the method, not of the tab: GFN-FF has no
        ddCOSMO, GBSA has eleven fewer solvents than ALPB under GFN-FF and
        thirteen fewer under GFN1, and a PM method has COSMO alone.  Whatever
        was chosen is kept if the new method can still do it, and dropped with
        a word if it cannot -- silently keeping it would leave the control
        showing one thing and the run doing another.
        """
        offered = _solvents.models_for(submit_ff_dd.value)
        submit_gfn_solvent.layout.display = '' if offered else 'none'
        # One model is not a choice; it is a label, and the tooltip on the
        # solvent box already says which it is.
        submit_gfn_solv_model.layout.display = '' if len(offered) > 1 else 'none'
        state['solvent_quiet'] = True
        try:
            if not offered:
                submit_gfn_solvent.value = ''
                return
            model_options = [(_solvents.model_label(name), name)
                             for name in offered]
            if list(submit_gfn_solv_model.options) != model_options:
                keeping = str(submit_gfn_solv_model.value or '')
                submit_gfn_solv_model.options = model_options
                submit_gfn_solv_model.value = (
                    keeping if keeping in offered else offered[0])
            model = _solv_model()
            wet = str(submit_gfn_solvent.value or '')
            names = _solvents.solvents_for(model, submit_ff_dd.value)
            wanted = [('none (gas phase)', '')] + [
                (_solvents.label_of(name), name) for name in names]
            if list(submit_gfn_solvent.options) != wanted:
                submit_gfn_solvent.options = wanted
                submit_gfn_solvent.value = wet if wet in names else ''
                if wet and wet not in names:
                    _set_mol_status(
                        f'{_solvents.model_label(model)} is not parametrised '
                        f'for {_solvents.label_of(wet)} here, so the solvent '
                        'was cleared. ALPB covers every solvent in the list.')
        finally:
            state['solvent_quiet'] = False

    def on_submit_solv_model(change):
        if change.get('name') != 'value' or state.get('solvent_quiet'):
            return
        _refresh_solvation_controls()
        wet = str(submit_gfn_solvent.value or '')
        model = _solv_model()
        if not wet:
            return
        _set_mol_status(
            f'{_solvents.label_of(wet)} with '
            f'{_solvents.model_label(model)} from now on. The structure on '
            'screen was optimised under the previous model, and the two do '
            'not agree -- on a glycine in water ALPB and ddCOSMO differed by '
            '2.8 kcal/mol -- so it is worth optimising again.')

    def on_submit_solvent(change):
        if change.get('name') != 'value' or state.get('solvent_quiet'):
            return
        wet = str(submit_gfn_solvent.value or '')
        model = _solv_model()
        _set_mol_status(
            f'Optimising in {_solvents.label_of(wet)} from now on '
            f'({_solvents.model_label(model)}). The structure on screen was '
            'not, so it is worth optimising again.'
            if wet else
            'Optimising in the gas phase from now on.')

    def on_submit_autospin(change):
        if change.get('name') != 'value':
            return
        scanning = bool(submit_gfn_autospin.value)
        submit_gfn_mult.disabled = scanning
        if scanning:
            _set_mol_status(
                'The multiplicity is scanned: the ones the electron count '
                'allows are each optimised and the lowest is kept. Three runs '
                'instead of one.')

    def _fill_charge_from_smiles():
        """Take the charge off the SMILES the structure was built from.

        A SMILES states the formal charges outright, so asking the user for a
        number the input already carries is asking them to repeat themselves.
        Nothing else can be read that way -- a pasted XYZ says nothing about
        charge, and no input says anything about the spin -- so those stay the
        user's to set.  Returns the SMILES it read, or '' when there was none.

        Except once the user has set it.  A structure is often built from a
        neutral SMILES and then something charged is put next to it by hand --
        a chloride beside an alkyl bromide is exactly that -- and the SMILES
        still says zero.  Every change of method then quietly put it back to
        zero, and at zero those systems have an odd number of electrons and
        xtb refuses every step: nothing runs, nothing moves, and the editor
        looks broken rather than misconfigured.  A number the user has typed
        is theirs until they load something else.
        """
        smiles = str((state.get('converted_xyz_cache') or {}).get('smiles') or '')
        if not smiles:
            return ''
        if state.get('charge_is_the_users'):
            return smiles
        if state.get('charge_filled_for') == smiles:
            # Read once, when the structure was made, and not again.  Reading
            # it afresh every time a method is chosen is what put a hand-set
            # charge back to zero -- and worse, a charge that had come out of
            # the SMILES correctly went to zero too, because by then the
            # cached SMILES was no longer the one that carried it.  Deriving a
            # setting belongs to the moment the thing it describes is built.
            return smiles
        try:
            state['charge_filling'] = True
            submit_gfn_charge.value = int(get_smiles_charge(smiles))
            state['charge_filled_for'] = smiles
        except Exception:
            return ''
        finally:
            state['charge_filling'] = False
        return smiles

    def on_submit_gfn_charge(change):
        """Remember that the charge is the user's now."""
        if change.get('name') != 'value' or state.get('charge_filling'):
            return
        state['charge_is_the_users'] = True

    def _server_method(value=None):
        """Whether this method is computed on the server rather than the page.

        Two engines now: xtb for the GFN family, MOPAC for the PM one. What
        the rest of the tab needs to know about either is the same -- there is
        no force field in the browser for it, so the live relaxation and the
        drag behave differently -- so it asks one question instead of naming
        engines in a dozen places.
        """
        chosen = submit_ff_dd.value if value is None else value
        return _gfn.is_gfn_method(chosen) or _mopac.is_mopac_method(chosen)

    def _server_binary(value):
        """The program this method needs, whichever engine it belongs to."""
        if _mopac.is_mopac_method(value):
            return _mopac.find_mopac()
        return _gfn.find_binary(value)

    def _server_label(value):
        if _gfn.is_gfn_method(value):
            return _gfn.GFN_METHODS[value]['label']
        if _mopac.is_mopac_method(value):
            return _mopac.MOPAC_METHODS[value]['label']
        return str(value).upper()

    def _live_ff_method():
        """The method the browser relaxes with while an atom is dragged.

        GFN runs on the server; a round trip per frame would cap a drag at
        about 13 Hz, so dragging keeps the force field that lives in the
        browser.  Choosing GFN changes what Optimise does, not what a drag
        does, and the status line says so once rather than silently.
        """
        chosen = str(submit_ff_dd.value or 'uff')
        return 'uff' if _gfn.is_gfn_method(chosen) else chosen

    def _offer_xtb_install(offer):
        """Show the offer, or the two buttons that carry it out, or neither."""
        submit_xtb_install_btn.layout.display = '' if offer else 'none'
        submit_xtb_confirm_btn.layout.display = 'none'
        submit_xtb_cancel_btn.layout.display = 'none'

    def _missing_tool():
        """Which program the chosen method needs and has not got.

        g-xTB is not the xtb beside it: it ships as a build of its own, and an
        ordinary xtb accepts --gxtb and silently runs GFN2.  So what is missing
        is a question with two answers, and offering the wrong one would
        install something that changes nothing.
        """
        method = str(submit_ff_dd.value)
        if not _gfn.is_gfn_method(method):
            return None
        return None if _gfn.find_binary(method) else (
            'gxtb' if method == 'gxtb' else 'xtb')

    def _refresh_xtb_offer():
        """A missing program under a GFN method is a thing to fix, not a wall."""
        missing = _missing_tool()
        state['xtb_install_tool'] = missing
        submit_xtb_install_btn.description = (
            'Install g-xTB' if missing == 'gxtb' else 'Install xtb')
        _offer_xtb_install(
            missing is not None
            and _gfn.install_script() is not None
            and not state.get('xtb_installing'))

    def on_submit_xtb_install(_button=None):
        """The first press only says what the second one would do.

        Fetching a conda environment of a few hundred megabytes is not a thing
        to start on a single click, and on a cluster the right answer is often
        to load the module instead -- so the command is on screen before it
        runs, and cancelling is as easy as agreeing.
        """
        command = _gfn.install_command(state.get('xtb_install_tool') or 'xtb')
        if command is None:
            _set_mol_status('DELFIN\'s installer is not next to this copy of '
                            'the dashboard, so there is nothing to run.')
            return
        submit_xtb_install_btn.layout.display = 'none'
        submit_xtb_confirm_btn.layout.display = ''
        submit_xtb_cancel_btn.layout.display = ''
        _set_mol_status(
            'This will run:  ' + ' '.join(command),
            'It fetches xtb from conda-forge into DELFIN\'s own tool '
            'directory -- a few hundred megabytes, several minutes, and it '
            'needs the network. On a cluster, module load xtb may be the '
            'better answer.',
        )

    def on_submit_xtb_cancel(_button=None):
        _refresh_xtb_offer()
        _set_mol_status('Left alone. Optimise will say xtb is missing rather '
                        'than doing nothing.')

    def on_submit_xtb_confirm(_button=None):
        if state.get('xtb_installing'):
            return
        state['xtb_installing'] = True
        _offer_xtb_install(False)
        _set_mol_status(
            f'Installing {state.get("xtb_install_tool") or "xtb"}...',
            spinner=True)

        def _work():
            def _line(text):
                if text.strip():
                    schedule_ui_update(
                        _set_mol_status,
                        f'Installing {state.get("xtb_install_tool") or "xtb"}'
                        '...', text.strip(), spinner=True)

            outcome = _gfn.install_xtb(
                on_line=_line, tool=state.get('xtb_install_tool') or 'xtb')

            def _done():
                state['xtb_installing'] = False
                _refresh_xtb_offer()
                if outcome.get('ok'):
                    # Asked again rather than assumed: the resolver is what
                    # every run will use, and a binary it cannot see is not
                    # installed as far as the dashboard is concerned.
                    _set_mol_status(outcome['status'],
                                    'Optimise will use it from now on.')
                else:
                    _set_mol_status(outcome.get('status') or 'The install '
                                    'did not finish.',
                                    *outcome.get('lines', [])[-2:])

            schedule_ui_update(_done)

        _start_background(_work, 'The install',
                          guards={'xtb_installing': False})

    def _refresh_method_controls():
        """Show what the chosen method can do, and nothing else.

        Every control here belongs to one engine or the other, and a control
        that cannot act under the method on screen is worse than absent: it
        invites the press, does nothing, and says nothing about why.  Three
        engines share this toolbar -- the browser's own field, xtb and MOPAC --
        and what each of them has is different.

        Settle goes under both server methods, and that is a removal rather
        than a tidy-up.  It was a short relaxation on release, which is a thing
        the browser's field can do cheaply; on the server it became the
        ordinary optimisation run without a cycle cap, which is what Auto
        already does when an atom is let go.  Two switches for one behaviour is
        one too many, and the spare one was the one still doing it with Dynamik
        Opt and Optimise both switched off -- a structure relaxing on release
        with nothing on screen admitting to running it.
        """
        chosen = submit_ff_dd.value
        server = _server_method()
        xtb = _gfn.is_gfn_method(chosen)
        mopac = _mopac.is_mopac_method(chosen)
        # Charge and spin: the server engines are told both, the browser's
        # field has no notion of either.
        submit_gfn_charge.layout.display = '' if server else 'none'
        submit_gfn_mult.layout.display = '' if server else 'none'
        # Scanning the multiplicity is xtb's; MOPAC is given the one on screen.
        submit_gfn_autospin.layout.display = '' if xtb else 'none'
        if not xtb and submit_gfn_autospin.value:
            submit_gfn_autospin.value = False
        # What the labels may be asked to say, and whether the press that asks
        # what a structure is has anything to ask with. Both follow the
        # method, and both are absent rather than refusing where it cannot
        # answer -- see _refresh_label_what and submit_shape_btn.
        _refresh_label_what()
        submit_shape_btn.layout.display = '' if xtb else 'none'
        # An anchor and a set of charges belong to the method that measured
        # them. The budget already guards its own anchor that way; the charges
        # are simply dropped, because a label is drawn without anything
        # keyed on it and GFN-FF's charges on a GFN2 structure would be four
        # hundredths of an electron wrong with nothing on screen saying so.
        if state.get('atom_charges_method') not in ('', str(chosen), None):
            state.pop('atom_charges', None)
            state.pop('atom_charges_method', None)
            state.pop('atom_charges_for', None)
            _repaint_labels()
        # Strength is how many steps the browser's field takes per animation
        # frame, and that field does not run under a server method.
        submit_strength_slider.layout.display = 'none' if server else ''
        # Which hands there are.
        #
        # A pull is a force on an internal coordinate, so it needs an engine
        # that can be told to hold one: the browser's own field can, and xtb
        # can through its constrain block -- that is why the slider is no
        # longer hidden under a server method, which was what left dragging
        # under GFN setting coordinates outright.  MOPAC takes no held
        # internals from this editor at all, so there the pull has nothing to
        # act on and the follow step falls through to the rigid hand.
        #
        # Measured on an ethane with one hydrogen dragged 0.60 A along x, the
        # hand set to pull: GFN2 left the atom 0.44 A short of where the
        # cursor asked, which is the chemistry having its say, and PM7 put it
        # exactly there -- 0.0000 A from the cursor, the same geometry the
        # move hand gives to the last decimal.  Offering "pull with a force"
        # under MOPAC is offering the move hand under another name, and the
        # difference is invisible from the outside.  So under MOPAC there is
        # one hand, and the box says one.
        pulling = ('pull with a force', 'pull')
        moving = ('move the atom', 'move')
        state['hand_quiet'] = True
        try:
            if mopac:
                if str(submit_hand_dd.value) == 'pull':
                    # Given back on the way out, so a detour through PM7 does
                    # not quietly cost the hand the user was working with.
                    state['hand_was_a_pull'] = True
                submit_hand_dd.options = [moving]
            else:
                submit_hand_dd.options = [pulling, moving]
                if state.pop('hand_was_a_pull', False):
                    submit_hand_dd.value = 'pull'
        finally:
            state['hand_quiet'] = False
        # The slider is the pull's own setting, and so is the budget, so
        # both go where the pull goes -- and both also answer to the method,
        # which is why one function decides them and neither refresh can undo
        # the other.  Shown unconditionally the slider came back on every
        # change of method: pick the move hand under UFF, switch to GFN2, and
        # there it was again, offering to set the strength of a hand that is
        # not in use.
        _refresh_hand_controls()
        # And the scan is xtb's for the same reason.  Left visible under
        # UFF or PM7 a whole scan could be armed, with the line saying "or
        # press Run scan" -- an instruction that cannot work -- and the
        # refusal arrived only on the press.  Switching away with one armed
        # said nothing at all.
        submit_scan_btn.layout.display = '' if xtb else 'none'
        if not xtb:
            for widget in (submit_scan_way, submit_scan_to, submit_scan_steps,
                           submit_scan_dd, submit_scan_del, submit_scan_whole,
                           submit_scan_how, submit_scan_energy,
                           submit_scan_back, submit_scan_run_btn,
                           submit_scan_price_btn):
                widget.layout.display = 'none'
            if _scan_legs():
                _set_mol_status(
                    f'{_server_label(submit_ff_dd.value)} cannot walk a scan '
                    '-- that is xtb\'s. The armed legs are kept for when you '
                    'come back to a GFN method.')
        else:
            # And they come back with it.  The line above promises that the
            # armed legs are kept, and they were -- but nothing put their row
            # back on screen, so a detour through UFF and back to GFN2 left a
            # scan that was still armed with no list, no Run scan and no way
            # to reach either.  Kept and unreachable is the same as gone.
            _refresh_scan()
        # The transition-state press and the two boxes beside it: which starts
        # exist, and which of the ways to a saddle this engine can run.  Both
        # are read from the tables the runs themselves read -- ORCA's keywords
        # for the saddle optimiser, the climb's own list for the gradients it
        # knows how to ask for -- so a control and its refusal cannot drift
        # apart.  Neither is the whole GFN family, and they are no longer the
        # same three: ORCA has no keyword for g-xTB and never will, but it
        # publishes ExtOpt -- it writes an input, calls a program, reads an
        # energy and a gradient back -- so the saddle optimiser reaches it
        # where the climb cannot.
        #
        # The mark survives a change of method, the way an armed scan does: it
        # describes two structures, not a program.  Its button does not --
        # under UFF it used to take the press, say "marked as the start of a
        # path", and offer a walk that then answered "a path needs xtb", which
        # is two presses to be told the first could not have worked.
        _refresh_saddle_controls()
        # And the climb's own switch, which is not the press: it says which
        # way a release walks.  Left visible under a method it cannot ask for
        # a gradient from it refused only after the press, which is a button
        # that promises what it cannot do -- under the most accurate method in
        # the list, where a transition state is most worth having.
        submit_climb_btn.layout.display = (
            '' if str(chosen).lower() in _climb.CLIMB_METHODS else 'none')
        # Keep bonds works by watching what a follow step hands back and
        # taking back the ones that made or broke a bond -- so it needs a
        # follow step, and that is the kernel's, which runs for a server
        # method and for nothing else.  Under UFF the drag never leaves the
        # browser: _begin_gfn_follow answers no, the wall is never consulted,
        # and the switch was a switch with nothing behind it.
        #
        # Its value is left where it stands rather than switched off, for the
        # same reason as Auto: nothing reads it under a browser method, and
        # taking the toolbar away should not also take the setting.
        submit_topology_btn.layout.display = '' if server else 'none'
        # And the playback pace is the other way round: only a server engine
        # walks a path worth pacing.  The browser's field draws its own frames
        # as it computes them and there is nothing queued to play.
        submit_play_speed.layout.display = '' if server else 'none'
        if server:
            _push_play_speed()
        # Settle: the browser's alone, for the reason in the docstring.
        submit_settle_btn.layout.display = 'none' if server else ''
        if server and submit_settle_btn.value:
            submit_settle_btn.value = False
        # Auto is the other way round -- a server method's, and dead under a
        # browser one.  What it is for is going down to a minimum when an atom
        # is let go, and that is refused outright for a browser method: the
        # field is already running there, and Settle is the switch for what a
        # release does.  It kept one working half under UFF, resuming an
        # Optimise run a drag had interrupted, and that goes with it: the run
        # stands down and the switch comes back up, which is what the two
        # visible switches already say is happening.
        #
        # Hidden, and its value left alone -- unlike Settle, which is switched
        # off as well.  Nothing reads this under a browser method, so an Auto
        # left on there does nothing at all; switching it off would mean that
        # picking UFF for a moment and going back to GFN2 quietly cost the
        # user the switch they had set.  Settle is the opposite case: under a
        # browser method its value is what the page settles by, so it has to
        # be off in fact and not merely out of sight.
        submit_auto_btn.layout.display = '' if server else 'none'
        # Types, and the hybridisation box beside it.
        #
        # An atom type is a force-field idea and reaches the force field and
        # nothing else: every UFF parameter at a centre follows from its type,
        # so overruling one is how a double bond perception missed is put
        # back.  A server method has no such thing to overrule -- xtb and
        # MOPAC work the shape out from the electrons, every step -- and the
        # override never leaves the kernel.
        #
        # Driven on an ethene drawn long enough that the double bond is not
        # perceived, both carbons forced to sp2: under UFF and MMFF94 the
        # press hands a fresh set of parameters to the page, and under GFN-FF,
        # GFN2 and PM7 the same press takes the field *off* the page and hands
        # over nothing -- while the line said "2 carbons typed from their
        # partners ... 2 changed".  A report of a change that reached no
        # calculation is the worst of the three ways a control can fail here,
        # so the control goes where the change cannot land.
        submit_hyb_auto_btn.layout.display = 'none' if server else ''
        # The polyhedron is the same story with a different name: choosing one
        # installs a set of restraints pulling the donors onto its vertices,
        # and those restraints are terms in the browser's field.  Under a
        # server method the same choice said "the donors are pulled onto it"
        # and took the field off the page instead.  Turn goes with it -- it
        # steps between arrangements of a polyhedron that is not acting.
        #
        # Swap stays.  It is not a restraint but an edit: the page rotates the
        # two ligands onto each other's directions there and then, which is a
        # geometry every engine is handed afterwards.
        if server:
            for widget in (submit_hyb_dd, submit_poly_dd,
                           submit_poly_turn_btn):
                widget.layout.display = 'none'
                widget.disabled = True
        # A method without solvation gets no solvent box: a control that can
        # only produce a refusal is worse than no control.  Which models and
        # which solvents a method does have is the solvents module's answer,
        # and it differs for every one of them -- so both boxes are rebuilt
        # here rather than merely shown or hidden.
        _refresh_solvation_controls()

    def on_submit_ff_changed(change):
        if change.get('name') != 'value':
            return
        gfn = _server_method()
        # Which controls the new method has at all, before anything is said
        # about what it will do with them.
        _refresh_method_controls()
        # What the new engine will do with the values that are being held.
        # Said as a second line under whatever else the change has to report,
        # rather than in place of it: both are about the same choice.
        carried = _carry_constraints_to(submit_ff_dd.value)
        # Relax means the browser's own field running once per frame, and there
        # is no GFN engine in the browser to run.  Under GFN it means the other
        # half of the same idea: while an atom is being dragged, the rest of
        # the molecule follows it -- one short xtb run per push, and nothing at
        # all when nothing is being dragged.  It was switched off and hidden
        # here before there was anything for it to do.
        # A follow armed under the old method must not outlive the choice that
        # armed it: the toggle handler reads the method that is chosen *now*.
        _end_gfn_follow()
        submit_relax_btn.disabled = False
        # The page reads this switch itself, and only follows when the class is
        # on it: asking the kernel whether to follow would cost a round trip
        # per push, and a UFF drag would be pushing for no reason.
        if gfn:
            submit_relax_btn.add_class('submit-gfn-follow')
        else:
            submit_relax_btn.remove_class('submit-gfn-follow')
        # Dynamik Opt without an xtb behind it cannot do anything at all, so
        # it goes rather than sitting there being pressed to no effect.
        usable = not gfn or _gfn.find_binary(submit_ff_dd.value) is not None
        submit_relax_btn.layout.display = '' if usable else 'none'
        if not usable and submit_relax_btn.value:
            submit_relax_btn.value = False
        label = (_server_label(submit_ff_dd.value)
                 if gfn else str(submit_ff_dd.value).upper())
        submit_relax_btn.tooltip = (
            f'While this is on, dragging an atom lets the rest of the molecule '
            f'follow it -- {label} on the server, a few cycles per push, and '
            f'each step says how long it took. Nothing runs while nothing is '
            f'being dragged.'
            if gfn else
            'Relax the structure continuously while you work on it.'
        )
        submit_settle_btn.tooltip = (
            f'When you let go of an atom, let {label} tidy the structure '
            f'around its new position instead of leaving it where the cursor '
            f'stopped. One short run on the server per release.'
            if gfn else
            'When you let go of an atom, let the structure relax around its '
            'new position instead of keeping the strain of the drag. Switch '
            'off to leave atoms exactly where you put them.'
        )
        # The engine follows the method, and the switch keeps its position:
        # switching back and forth is something a user does constantly, and it
        # must not cost a press each time -- nor leave the previous engine
        # running underneath the new choice.
        if gfn:
            _stop_browser_field()
        elif submit_relax_btn.value:
            _enable_live_forcefield()
            _run_manip_js(
                'if(window.__delfinSubmitManip)'
                'window.__delfinSubmitManip.startAutoOptimize('
                f'{json.dumps(submit_scope_id)});'
            )
        _refresh_xtb_offer()
        said = ''
        if gfn:
            label = _server_label(submit_ff_dd.value)
            source = _fill_charge_from_smiles()
            if _server_binary(submit_ff_dd.value) is None:
                # Which program is missing depends on the method, and so does
                # whether there is an installer for it: MOPAC is not fetched
                # by the xtb button, and saying "needs xtb" under PM7 sent
                # people looking for the wrong thing.
                needs = ('MOPAC' if _mopac.is_mopac_method(submit_ff_dd.value)
                         else 'g-xTB' if submit_ff_dd.value == 'gxtb' else 'xtb')
                offered = (_missing_tool() is not None
                           and _gfn.install_script() is not None)
                said = (
                    f'{label} needs {needs}, which was not found. The button '
                    'fetches it with DELFIN\'s own installer, and says what '
                    'it will run before it runs it.'
                    if offered else
                    f'{label} needs {needs} on the PATH; it was not found. '
                    'Optimise will say so rather than doing nothing.')
            elif source:
                said = (
                    f'Optimise now uses {label}. Charge {submit_gfn_charge.value} '
                    f'read from the SMILES; the multiplicity is yours to set. '
                    'Switch Dynamik Opt on to have the molecule follow an atom '
                    'you drag.')
            else:
                said = (
                    f'Optimise now uses {label}. Set the charge (q) and the '
                    f'multiplicity (M): {label} needs both, and a wrong spin '
                    'on a metal gives a confident wrong answer rather than an '
                    'error. Switch Dynamik Opt on to have the molecule follow '
                    'an atom you drag.')
        if said or carried:
            _set_mol_status(said, carried)

    def on_submit_manip_clear(_button=None):
        _ensure_manip_bootstrap()
        _run_manip_js(
            f'if(window.__delfinSubmitManip) '
            f'window.__delfinSubmitManip.clear({json.dumps(submit_scope_id)});'
        )

    def on_submit_manip_undo(_button=None):
        """One step back through everything that has been done.

        Through the one history, never the browser's own stack: that one is
        cleared by every re-render, so a drawn atom used to throw away every
        drag before it, and Undo then stopped where the history it happened to
        be reading ran out.
        """
        state.pop('pre_optimize_frames', None)
        _undo_structure()

    def on_submit_manip_redo(_button=None):
        """One step forward again, through what Undo took back."""
        state.pop('pre_optimize_frames', None)
        _redo_structure()

    def on_submit_manip_sync(change):
        if change.get('name') != 'value':
            return
        new_xyz = submit_manip_sync.value
        if not new_xyz or not new_xyz.strip():
            return
        # The structure as the page sent it, comment line and all -- and the
        # comment line is the half that matters, because it is where the page
        # says whether this is a hand following, a hand letting go, an undo in
        # the browser or the force field reporting from under all three.
        record('manip', v=new_xyz)
        # Extract only the new coordinate lines; drop JS-side count + comment.
        new_lines = new_xyz.splitlines()
        if len(new_lines) >= 2:
            try:
                int(new_lines[0].strip())
                coord_lines = new_lines[2:]
            except ValueError:
                coord_lines = new_lines
        else:
            coord_lines = new_lines
        coord_body = '\n'.join(xyz_body(coord_lines))
        # The count belongs to the body, not to whatever was in the box
        # before.  Kept from the old header, three rows of water landed under
        # a twelve of benzene and the box held a malformed XYZ that everything
        # downstream went on reading.
        coord_rows = [one for one in coord_body.splitlines() if one.strip()]
        # Preserve the user's original header (atom count + comment line) if
        # present in the current coords_widget value.
        #
        # The user's -- a name off a file, a note they typed -- and not one of
        # this editor's own.  These coordinates come from the browser's model:
        # they are what the hand has just done to the structure, and carrying
        # the old comment over put "Optimised in DELFIN viewer" above a
        # geometry that had been dragged out of shape since.  The word then
        # says where the box was last written from rather than what is in it,
        # and a benzene with a hydrogen 2.66 A off its carbon reads as the
        # result of an optimisation.  A comment this editor wrote is replaced;
        # anything else is left exactly as it is.
        old_lines = coords_widget.value.splitlines()
        header = ''
        if len(old_lines) >= 2:
            try:
                int(old_lines[0].strip())
                kept = old_lines[1]
                if _is_editor_comment(kept):
                    kept = 'Edited in DELFIN viewer'
                # The comment is the user's and is kept; the count is the
                # body's.  Carried over from the old header, three rows of
                # water landed under a twelve of benzene and everything
                # downstream went on reading a malformed XYZ.
                header = f'{len(coord_rows)}\n{kept}\n'
            except ValueError:
                pass
        # A drag has just finished. If a polyhedron is being held, work out
        # again which donor is now nearest which vertex: dragging a ligand
        # towards another position and having the field haul it straight back
        # is not an exchange, it is a fight. Recomputing here means the
        # polyhedron accepts the ligand where it has been put and pulls it the
        # rest of the way onto the vertex it is now closest to.
        lines = new_xyz.splitlines()
        note = lines[1].strip() if len(lines) > 1 else ''
        drag_ended = note.startswith('DELFIN drag-end')
        dragging = note.startswith('DELFIN drag-follow')
        # The hand letting go, kept apart from drag_ended because an undo
        # joins that one further down and an undo is not a drag: it hands back
        # a geometry that was already there, so there is nothing to price and
        # nothing to refuse.  Held to the budget it would be answered with the
        # very structure it is undoing.
        released = drag_ended
        # Sent while the mouse is still down, so the molecule can follow the
        # atom rather than wait for it to be let go.  The comment line names
        # the atoms the hand is on, so the answer can keep them there.
        if dragging:
            holding = []
            for word in note.split():
                if word.startswith('held='):
                    holding = [int(n) for n in word[5:].split(',')
                               if n.strip().lstrip('-').isdigit()]
            _gfn_follow_step(new_xyz, holding)
            # A climb does not fight a hand: the grab stopped whatever was
            # walking, and this is the second way of hearing that a hand is
            # here, for a drag that reaches the kernel without one.  Only the
            # mark the climb aims along has to be caught up; a run that is
            # still going despite the grab is one this message can end.
            if _climb_owns_the_release():
                if not state.get('climb_was'):
                    state['climb_was'] = (
                        (state.get('climb_showing')
                         if state.get('climb_run') is not None else None)
                        or _current_xyz())
                # Read here, off the message that carries the hand, because
                # this is where the two halves of the pair are both in hand:
                # the atom being dragged comes with this message and the atom
                # it is being dragged at is the one already picked.
                _name_the_pair_the_hand_means(holding)
                _interrupt_gfn()
        if (drag_ended and state.get('poly_applied')
                and state.get('poly_metal') is not None):
            # Only a real end of a drag, not the twice-a-second heartbeat the
            # running optimiser sends: reassigning on every heartbeat reloaded
            # the whole field twice a second and never let a moved ligand
            # settle onto its new vertex.
            state['poly_assignment'] = None
            state['poly_recheck'] = True

        # Undoing a drag in the browser is a change the user made, the same as
        # making it was.  It arrived without a reason on it, so it took neither
        # of the two paths below: the optimisation went on running over a
        # geometry that had just been taken back, and whichever of the two
        # wrote the coordinate box last is what the user was left with.
        if note.startswith('DELFIN undo'):
            drag_ended = True

        # A drag the budget is pricing does not write its own coordinates.
        #
        # This is where the ceiling was being got round, and it needed no race
        # and no cleverness: the page's model landed in the box on every
        # message of every drag, unconditionally, and nothing on this path
        # ever asked what it cost.  The wall ran a moment behind in a thread
        # of its own and wrote the box only when it refused, so the geometry
        # that *survived* a drag had never been priced at all.  Measured, an
        # ethane methyl dragged out at 298 K under a 22.3 kcal/mol ceiling was
        # left at +141.2 -- and under the rigid hand it was worse, because the
        # first refusal sets thermal_spent and the follow then deliberately
        # stands still rather than shake, so from that moment nothing was
        # computing and every later mouse position went into the box untouched.
        #
        # That is also the whole of "sometimes it works at room temperature":
        # whether a drag was stopped came down to whether the wall happened to
        # write last, which is a race between xtb at a tenth of a second and a
        # mouse.  Let go a moment after a refusal and the ceiling looked real;
        # keep dragging and it was not there.
        #
        # So while something is pricing this drag, the box belongs to it: the
        # follow writes the geometry it priced, and the page's wish reaches the
        # picture and nothing else.  With the budget off, or with no anchor to
        # measure against, this is the drag it always was.
        #
        # The release lands while the follow is still on, and not by luck: the
        # page sends the geometry from its mouseup handler and says the hand
        # has gone from the animation frame after it, which the browser runs
        # once that handler has finished.  That order is what makes the last
        # write of a drag a priced one rather than a raw one, so it is worth
        # knowing that it is the event loop guaranteeing it.
        walled = ((dragging or released)
                  and bool(submit_thermal_btn.value)
                  and state.get('gfn_follow')
                  and _thermal_budget()[0] is not None)
        if drag_ended:
            # Set, Hold, a bond edit and a drag all arrive here.  Any of them
            # during an optimisation makes what xtb is doing about a structure
            # that is no longer on the screen, so it starts again from the one
            # that is -- and the geometry lands first, so it is the new one it
            # starts from.
            _interrupt_gfn()
            _arm_gfn_restart()
            # Whichever of the two was walking starts again from the structure
            # that was just made -- and a climb starts aimed along the way it
            # was made, which is the one thing about it that is not the
            # minimisation's.  The pull is never part of the climb, because a
            # restrained saddle is not a saddle: measured, a climb on a surface
            # with the forming bonds held converges onto points with two
            # imaginary modes or with none, half an angstrom from the real
            # saddle, while the same climb resumed after the hand let go
            # reaches it in 39 steps.  The geometry lands here first and the
            # restart reads the box, so there is nothing to hand over by hand.
        elif state.get('optimize_run') is not None:
            # Not an edit, and something is optimising.  The browser's own
            # field reports where it has got to twice a second and once more
            # when it is switched off; those are the field talking about a
            # structure the optimiser now owns, and writing them puts the
            # picture's idea of the geometry over the calculation's.  That is
            # the shape the "Optimised in DELFIN viewer" header was found in:
            # a box holding something no run had produced.  The optimisation
            # owns the box until it is done or a hand takes it back.
            return

        if walled:
            # The picture has the drag and the box has what was priced.
            #
            # At the release that is the last geometry the budget agreed to,
            # which the follow has usually written already -- so this is a
            # no-op whenever the answers kept up with the hand, and a spring
            # back to the last affordable structure when they did not.  Either
            # way what the user is left with is a structure that was measured
            # and allowed, which is the whole of what the temperature means.
            if released:
                _keep_the_priced_geometry()
            if state.pop('poly_recheck', False):
                schedule_ui_update(_enable_live_forcefield)
            return

        payload = header + coord_body
        # The guard is cleared by update_view, which traitlets only
        # calls when the value actually changes. Dragging an atom out and back,
        # or any edit that lands on the same coordinates, would otherwise leave
        # the flag set for the rest of the session and swallow the user's next
        # genuine edit of the coordinate box.
        if coords_widget.value == payload:
            state['manip_inflight'] = False
            return
        state['manip_inflight'] = True
        # The same molecule, moved -- not one the host has never seen.
        #
        # This is the drag's own way back into the box, and unguarded it made
        # the host start the structure over on every single drag: charge to
        # zero, multiplicity to one, held values and history gone, and
        # structure_changed dropping the thermal anchor with them.  Without an
        # anchor there is no budget, so nothing could be refused and nothing
        # taken back -- which is why the budget appeared to do nothing at all,
        # measured: the anchor was gone after the first yank and the charge
        # with it.
        state['structure_edit_inflight'] = True
        try:
            coords_widget.value = payload
        finally:
            state['structure_edit_inflight'] = False
        if state.pop('poly_recheck', False):
            # After the coordinates have landed, so the assignment is worked
            # out from where the ligands actually are now.
            schedule_ui_update(_enable_live_forcefield)

    # The journal listens before anything else does, and that ordering is the
    # whole of why it is written here rather than at the end of build().
    # Traitlets calls observers in the order they were registered, so a
    # recorder registered after the handlers would time-stamp a widget change
    # *after* the status line and the frames that change produced -- and the
    # readable half of a report would then say the editor answered a question
    # before it was asked.
    journal_watching = journal.watch(locals())
    submit_gfn_frame.observe(journal.on_frame, names='value')
    submit_gfn_wall.observe(journal.on_wall, names='value')
    # The coordinate box is the one channel the editor shares with its host,
    # and every write to it already carries a comment saying what produced it.
    # Watching it here means the geometry's history costs nothing to collect:
    # nobody has to remember to record it at each of the places that write.
    coords_widget.observe(journal.on_box, names='value')
    # The other half of what keeps the scan's profile honest.  A run that
    # starts takes it away at :func:`_note_the_run`; a geometry that arrives
    # in the box without one -- Undo, Reset, a structure someone loads, a drag
    # writing where it ended -- is caught here, and the picture stands only
    # while the box still holds the walk it is about.
    coords_widget.observe(_on_box_for_scan_plot, names='value')
    journal.opening(journal.snapshot(journal_watching), coords_widget.value)

    submit_select_btn.observe(on_submit_select_toggle, names='value')
    submit_manip_btn.observe(on_submit_manip_toggle, names='value')
    submit_manip_clear_btn.on_click(on_submit_manip_clear)
    submit_manip_undo_btn.on_click(on_submit_manip_undo)
    submit_manip_redo_btn.on_click(on_submit_manip_redo)
    submit_centre_btn.on_click(on_submit_centre)
    submit_relax_btn.observe(on_submit_relax_toggle, names='value')
    # On the page from the start: a player installed at click time races the
    # run it is meant to show.
    _install_gfn_frame_watcher()
    submit_ff_dd.observe(on_submit_ff_changed, names='value')
    # And once for the method the box already shows: a host that opens the
    # editor on a saved choice never fires the handler, so the toolbar would
    # stand there offering the previous method's controls.
    _refresh_method_controls()
    submit_xtb_install_btn.on_click(on_submit_xtb_install)
    submit_xtb_confirm_btn.on_click(on_submit_xtb_confirm)
    submit_xtb_cancel_btn.on_click(on_submit_xtb_cancel)
    submit_gfn_autospin.observe(on_submit_autospin, names='value')
    submit_gfn_solvent.observe(on_submit_solvent, names='value')
    submit_gfn_solv_model.observe(on_submit_solv_model, names='value')
    submit_labels_btn.observe(on_submit_labels_toggle, names='value')
    submit_label_what.observe(on_submit_label_what, names='value')
    submit_label_size.observe(on_submit_label_size, names='value')
    submit_strength_slider.observe(on_submit_strength_changed, names='value')
    submit_pull_slider.observe(on_submit_pull_changed, names='value')
    submit_hand_dd.observe(on_submit_hand_changed, names='value')
    submit_scan_way.observe(on_submit_scan_way, names='value')
    submit_topology_btn.observe(on_submit_topology, names='value')
    submit_scan_whole.observe(on_submit_scan_whole, names='value')
    submit_scan_back.observe(on_submit_scan_back, names='value')
    submit_scan_how.observe(on_submit_scan_how, names='value')
    submit_path_from_btn.on_click(on_submit_path_from)
    submit_saddle_btn.on_click(on_submit_saddle)
    submit_shape_btn.on_click(on_submit_shape)
    # The start decides which ways there are -- there is no walk to stop after
    # when the start is the structure on screen -- so it goes through the one
    # place that works both boxes out; the way only renames the press.
    submit_saddle_from.observe(lambda _c: _refresh_saddle_controls(),
                               names='value')
    submit_saddle_how.observe(lambda _c: _name_the_saddle_press(),
                              names='value')
    # Which mode the two presses are about is read when one of them is
    # pressed, so the box needs no observer of its own: nothing happens at the
    # moment it is changed, and a Hessian is not recomputed by picking a
    # different mode out of the one that was already taken.
    submit_mode_btn.on_click(on_submit_show_mode)
    submit_ends_btn.on_click(on_submit_follow_down)
    submit_climb_btn.observe(on_submit_climb, names='value')
    submit_sens_slider.observe(on_submit_sens_changed, names='value')
    submit_play_speed.observe(on_submit_play_speed, names='value')
    submit_thermal_btn.observe(on_submit_thermal, names='value')
    submit_temperature.observe(on_submit_thermal_terms, names='value')
    submit_thermal_anchor_btn.on_click(on_submit_thermal_anchor)
    submit_settle_btn.observe(on_submit_settle_toggle, names='value')
    submit_auto_btn.observe(on_submit_auto_toggle, names='value')
    submit_dyn_bonds_btn.observe(on_submit_dyn_bonds, names='value')
    submit_pick_sync.observe(on_submit_pick_sync, names='value')
    submit_poly_dd.observe(on_submit_poly_changed, names='value')
    submit_hyb_dd.observe(on_submit_hyb_changed, names='value')
    submit_hyb_auto_btn.on_click(on_submit_hyb_auto)
    submit_poly_turn_btn.on_click(on_submit_poly_turn)
    submit_cmd_sync.observe(on_submit_cmd, names='value')
    submit_draw_btn.observe(on_submit_draw_toggle, names='value')
    submit_element_dd.observe(on_submit_draw_choice, names='value')
    submit_hold_btn.on_click(on_submit_hold)
    submit_gfn_charge.observe(on_submit_gfn_charge, names='value')
    submit_scan_btn.on_click(on_submit_scan)
    submit_scan_del.on_click(on_submit_scan_del)
    submit_scan_run_btn.on_click(on_submit_scan_run)
    submit_scan_price_btn.on_click(on_submit_scan_price)
    submit_swap_btn.on_click(on_submit_swap)
    submit_hold_mode.observe(on_submit_hold_mode, names='value')
    submit_bond_btn.on_click(on_submit_bond)
    submit_unbond_btn.on_click(on_submit_unbond)
    submit_constraint_del.on_click(on_submit_constraint_del)
    submit_constraint_dd.observe(on_submit_constraint_selected, names='value')
    submit_internal_value.observe(on_submit_constraint_retune, names='value')
    submit_reset_btn.on_click(on_submit_reset)
    submit_internal_btn.observe(on_submit_set_internal, names='value')
    submit_internal_value.observe(on_submit_internal_value, names='value')
    # The page watches these buttons itself: waiting for the kernel to say the
    # switch went off costs a round trip, and the playback ran on for it.
    #
    # All three of them, because the page's question is "does a walk of this
    # structure still have a switch behind it" and three toggles can answer
    # yes.  Climb to TS was not one of them, so a climb read as a run whose
    # switch was up and its path was thrown away the moment it was queued --
    # which is why the climb had to pretend to be a followed hand, and why its
    # playback then ignored the speed slider.  See switchIsOn.
    submit_optimize_btn.add_class('submit-optimize-switch')
    submit_optimize_btn.observe(on_submit_optimize, names='value')
    submit_optimize_all_btn.add_class('submit-optimize-switch')
    submit_optimize_all_btn.observe(on_submit_optimize_all, names='value')
    submit_climb_btn.add_class('submit-optimize-switch')
    submit_manip_sync.observe(on_submit_manip_sync, names='value')

    # What the force field had to approximate belongs under the structure it
    # describes, not in the preview's status line where it competes with
    # conversion messages and scrolls away.
    submit_ff_notes = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', margin='4px 0 0 0'),
    )
    submit_ff_notes.add_class('submit-ff-notes')

    #: The picture of the path the last scan walked, under the structure it
    #: was walked on.
    #:
    #: Under it rather than on it.  The status line lies on the picture
    #: because it is written several times a second and a row above the viewer
    #: made the atom under the cursor step up and down; a profile is written
    #: once, at the end of a walk that took minutes, so it costs nothing to
    #: give it a row -- and a row of its own takes none of the pixels the
    #: structure is drawn in, which an overlay would.  It is ``display: none``
    #: until there is a walk to show, so an editor that has never scanned is
    #: laid out exactly as it was.
    #:
    #: It travels into fullscreen as a panel: bounded to 30vh and scrolling,
    #: which is the shared rule for results that belong to a picture (the RMSD
    #: pair, the Fukui numbers).  Fullscreen is still for the structure.
    submit_scan_plot = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', margin='4px 0 0 0',
                              display='none'),
    )
    submit_scan_plot.add_class('submit-scan-plot')
    submit_scan_plot.add_class('delfin-structure-fs-member')
    submit_scan_plot.add_class('delfin-structure-fs-panel')
    # ---- where a structure comes from ------------------------------

    convert_smiles_button = widgets.Button(
        description='CONVERT SMILES', button_style='info',
        layout=widgets.Layout(width='155px'),
        tooltip='Full isomer search',
    )
    convert_smiles_quick_button = widgets.Button(
        description='QUICK CONVERT SMILES', button_style='info',
        layout=widgets.Layout(width='185px'),
        tooltip='Fast single structure (no isomer search, no UFF)',
    )
    convert_smiles_uff_button = widgets.Button(
        description='CONVERT SMILES + UFF', button_style='info',
        layout=widgets.Layout(width='185px'),
    )


    manta_button = widgets.Button(
        description='MANTA', button_style='',
        layout=widgets.Layout(width='150px'),
        tooltip='MANTA: build the complete coordination-isomer manifold from the '
                'SMILES and rank it by GFN2 energy (best isomer/conformer first; '
                'needs xtb). Results show in the viewer with isomer navigation.',
    )
    # MANTA-logo teal/cyan.
    manta_button.style.button_color = '#1FA9C0'
    manta_button.style.font_weight = 'bold'

    # --- MANTA settings (the 5 keys a user actually needs; MANTA button sits BELOW) ---
    # Power-user knobs (construction / seeds / confs-per-isomer / rank-method /
    # merge-variants) are CLI-only on purpose: each has one sensible value for the
    # dashboard (construction is ALWAYS champion = best; the rest are redundant with
    # Quality), so exposing them only confuses users.  Pinned here.
    _MANTA_DASH_DEFAULTS = dict(construction="champion", num_confs=None,
                               collapse=False)
    # Quality preset -> conformer-seed count.  Selecting a preset auto-fills the
    # Seeds field (transparent: extreme = 60), but the field stays editable so a
    # user can dial a custom seed count on top of the preset's cap/templates.
    _MANTA_PROFILE_SEEDS = {'fast': 12, 'normal': 20, 'max': 40, 'extreme': 60}
    manta_quality = widgets.Dropdown(
        options=['fast', 'normal', 'max', 'extreme'], value='extreme',
        description='Quality:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='160px'),
        tooltip='Conformer depth (the main completeness<->speed switch). extreme guarantees the '
                'GFN2 global minimum is in the manifold (convergence-proven); fast/normal miss it '
                'by ~2.5 kcal/mol on multi-isomer systems.')
    manta_seeds = widgets.IntText(
        value=_MANTA_PROFILE_SEEDS['extreme'], description='Seeds:',
        style={'description_width': 'initial'}, layout=widgets.Layout(width='130px'),
        tooltip='ETKDG conformer seeds. Auto-fills from Quality (extreme=60) for transparency; '
                'edit it to use a custom seed count (more seeds = more conformers, slower).')

    def _sync_seeds_to_quality(change):
        try:
            manta_seeds.value = _MANTA_PROFILE_SEEDS.get(change['new'], manta_seeds.value)
        except Exception:
            pass
    manta_quality.observe(_sync_seeds_to_quality, names='value')
    manta_max_iso = widgets.IntText(
        value=0, description='Max isomers (0=ALL):', style={'description_width': 'initial'},
        layout=widgets.Layout(width='190px'),
        tooltip='Cap emitted isomers. 0 = COMPLETE manifold, never cut off (recommended).')
    # DEFAULT = No / No -> users get ONLY the manifold (no post-processing).  Rank and Opt are
    # SEPARATE opt-ins: Rank reorders by xtb single-point energy (geometry UNCHANGED); Opt xtb-
    # geometry-optimises structures (changes geometry).  You can Rank without Opt, or Opt without Rank.
    manta_rank = widgets.Dropdown(
        options=['No', 'gfn2', 'gfnff', 'gfn1', 'gfn0'], value='No',
        description='Rank:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='150px'),
        tooltip='Energy-RANK the manifold so the global-minimum isomer/conformer is first — '
                'single-point xtb energy, geometry UNCHANGED (pure reordering, no post-processing of '
                'the structures). No = keep construction order (already best-first by realism). '
                'gfn2 = standard; gfnff = fast force-field.')
    manta_opt = widgets.Dropdown(
        options=['No', 'Top 5', 'Top 10', 'Top 20', 'All'], value='No',
        description='Opt:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='140px'),
        tooltip='Geometry-OPTimise structures to a clean final geometry (xtb; method FOLLOWS Rank, or '
                'gfn2 if Rank=No). No = keep the construction geometry (pure manifold, DEFAULT). '
                'Top-N = optimise the N best; All = the whole manifold (slowest/best). '
                'Independent of Rank.')
    manta_spin = widgets.Dropdown(
        options=['auto', '1', '2', '3', '4', '5', '6', '7'], value='auto',
        description='Spin:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='130px'),
        tooltip='Spin multiplicity for the GFN2 energy/rank. auto = scan parity/+2/+4 and take the '
                'GFN2 ground state (parity-correct for open-shell TM). Or FIX it (1=singlet, 2=doublet, '
                '3=triplet, ...) when you know the state. NOTE: GFN2 spin energetics for TM are '
                'approximate -> set it explicitly for accuracy.')
    manta_det = widgets.Dropdown(
        options=['On', 'Off'], value='On',
        description='Determinism:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='160px'),
        tooltip='DETERMINISTIC construction (DEFAULT On) — byte-identical, reproducible manifold, '
                'IDENTICAL to what the CLI and the development loop build (ship = validate). '
                'Off = non-deterministic embedding (racy frame count; only for exploring extra '
                'ETKDG variety). Keep On for a stable, reproducible isomer x conformer manifold.')
    manta_settings_row = widgets.VBox([
        widgets.HTML("<b style='color:#1FA9C0'>MANTA settings</b> "
                     "<span style='color:#888;font-size:90%'>— complete coordination-isomer "
                     "&times; conformer manifold</span>"),
        widgets.HBox([manta_quality, manta_seeds, manta_max_iso, manta_rank, manta_opt, manta_spin, manta_det],
                     layout=widgets.Layout(gap='12px', flex_wrap='wrap', align_items='center')),
    ], layout=widgets.Layout(border='1px solid #d0e7ec', padding='8px', margin='4px 0'))


    xyz_copy_btn = widgets.Button(
        description='\U0001f4cb Copy Coordinates', button_style='success',
        layout=widgets.Layout(width='150px'), disabled=True,
    )
    xyz_copy_status = widgets.HTML(value='', layout=widgets.Layout(margin='0 0 0 6px'))

    isomer_prev_btn = widgets.Button(
        description='\u25c0', button_style='info',
        layout=widgets.Layout(width='35px'),
    )
    isomer_next_btn = widgets.Button(
        description='\u25b6', button_style='info',
        layout=widgets.Layout(width='35px'),
    )
    isomer_label = widgets.HTML(
        value='', layout=widgets.Layout(width='180px'),
    )
    isomer_nav_row = widgets.HBox(
        [isomer_prev_btn, isomer_label, isomer_next_btn],
        layout=widgets.Layout(gap='4px', align_items='center', display='none'),
    )


    def _set_smiles_conversion_busy(is_busy):
        state['smiles_busy'] = bool(is_busy)
        set_buttons_disabled(
            [
                convert_smiles_button,
                convert_smiles_quick_button,
                convert_smiles_uff_button,
                isomer_prev_btn,
                isomer_next_btn,
            ],
            is_busy,
        )
        ctx.set_busy(state['smiles_busy'] or state['batch_preview_busy'])


    def _build_mol_output_bundle(xyz_data):
        view = py3Dmol.view(width='100%', height=viewer_height)
        view.addModel(xyz_data, 'xyz')
        apply_molecule_view_style(view)
        scope_key_js = json.dumps(submit_scope_id)
        # The editor rides along only until the page has said it has it. It is
        # 136 KiB of the 159 a rendered structure weighs, and every conversion,
        # every edit and every optimisation result sent it again -- to a page
        # that had it already and threw the copy away on its version check.
        # Sending it at all, the first time, is the belt: the separate script
        # that installs it goes through an output widget whose content can be
        # replaced before the page has run it, and a structure without an
        # editor cannot be edited at all.
        carry_editor = (state.get('manip_seen_version')
                        != submit_manip_version())
        registration = (
            '\n'
            + (submit_manip_bootstrap_js() if carry_editor else '')
            + '\n(function(){\n'
            '  try {\n'
            '    window._submitMolViewerByScope = window._submitMolViewerByScope || {};\n'
            # Every render of this tab creates a fresh WebGL viewer. Without
            # releasing the previous one its context, observers and window-level
            # mouse listeners stay alive; browsers cap contexts and start
            # killing the oldest, which blacks out the viewers in other tabs.
            f'    var prev = window._submitMolViewerByScope[{scope_key_js}];\n'
            # Keep the camera across a re-render. Optimising or stepping to
            # another isomer rebuilds the viewer, and the view otherwise snapped
            # back to the default orientation every time.
            '    window._submitMolViewByScope = window._submitMolViewByScope || {};\n'
            '    if (prev && prev !== viewer_UNIQUEID) {\n'
            '      try {\n'
            '        var prevModel = prev.getModel ? prev.getModel() : null;\n'
            '        var prevCount = prevModel && prevModel.selectedAtoms\n'
            '          ? prevModel.selectedAtoms({}).length : -1;\n'
            f'        window._submitMolViewByScope[{scope_key_js}] =\n'
            '          {view: prev.getView(), atoms: prevCount};\n'
            '      } catch(e) {}\n'
            '      if (window.__delfinDisposeViewer) window.__delfinDisposeViewer(prev);\n'
            '    }\n'
            '    try {\n'
            f'      var saved = window._submitMolViewByScope[{scope_key_js}];\n'
            '      var model = viewer_UNIQUEID.getModel();\n'
            '      var count = model && model.selectedAtoms\n'
            '        ? model.selectedAtoms({}).length : -1;\n'
            # Same structure, just moved -- keep the view. A different one
            # deserves a fresh look at it. An edit that adds or removes an
            # atom is neither: the count changes but it is still the molecule
            # being worked on, so the camera stays put rather than snapping
            # back to the default every time something is drawn.
            '      var edited = !!window.__delfinStructureEdit;\n'
            '      window.__delfinStructureEdit = false;\n'
            '      if (saved && saved.view && (edited || saved.atoms === count)) {\n'
            '        viewer_UNIQUEID.setView(saved.view);\n'
            '        viewer_UNIQUEID.__delfinUserInteracted = true;\n'
            '        viewer_UNIQUEID.render();\n'
            '      }\n'
            '    } catch(e) {}\n'
            f'    window._submitMolViewerByScope[{scope_key_js}] = viewer_UNIQUEID;\n'
            '    var el = document.getElementById("3dmolviewer_UNIQUEID");\n'
            '    var fire = function(){\n'
            '      if (window.__delfinSubmitManip) {\n'
            f'        window.__delfinSubmitManip.onViewerReady({scope_key_js}, el);\n'
            '      }\n'
            '    };\n'
            '    setTimeout(fire, 80);\n'
            '  } catch(e) {}\n'
            '})();\n'
        )
        # The numbers, from the shared editor part -- the same code the ORCA
        # Builder draws them with.  py3Dmol names its viewer viewer_UNIQUEID
        # and rewrites that when it makes the HTML, so the labels can be
        # addressed to it by name like everything else here.
        labels = ''
        if submit_labels_btn.value:
            labels = show_atom_numbers_js(
                var='viewer_UNIQUEID',
                scale=scale_for_px(submit_label_size.value),
                texts=_label_texts_now())
        if hasattr(view, 'startjs'):
            view.startjs += registration
            if labels:
                view.startjs += '\n' + labels + '\n'
        html_payload = view._make_html()
        return ({
            'output_type': 'display_data',
            'data': {
                'application/3dmoljs_load.v0': html_payload,
                'text/html': html_payload,
            },
            'metadata': {},
        },)

    def _clear_mol_output():
        show_output(())

    def _replace_mol_output_text(*lines):
        _set_mol_status(*lines)
        _clear_mol_output()
        _set_manip_toolbar_enabled(False)


    def _replace_mol_output_view(xyz_data):
        # The structure as it arrived, once: the first thing in the history,
        # the thing Reset returns to, and the thing Undo stops at.
        if state.pop('history_seed_pending', False):
            state['history'] = []
            _remember('the structure as it was loaded')
        _clear_mol_status()
        _ensure_manip_bootstrap()
        mol_output.outputs = _build_mol_output_bundle(xyz_data)
        _set_manip_toolbar_enabled(True)
        # A marked end and the structure on screen are a pair only once they
        # are two different structures, so the start that names them appears
        # when the second one is drawn and not when the mark is made. Asked
        # here because this is the one place every host goes through when what
        # is on screen changes.
        _refresh_saddle_controls()
        # And the finished profile is about one molecule the same way. It
        # outlives the structure it was walked on -- loading another one does
        # not throw it away -- so the press that prices it again has to leave
        # when the molecule does, or it would offer a second opinion about
        # something that is no longer on screen.
        _refresh_scan()
        _push_hand_bonds()

    def _show_mol_busy(message):
        """Render the animated MANTA loader centered in the viewer.

        While the manifold build / GFN2 ranking / optimization runs there is
        no structure to show, so fill the 3D viewer box with the floating
        MANTA mark plus a status caption. ``_replace_mol_output_view`` later
        swaps it out for the finished structure (same as before).
        """
        _clear_mol_status()
        safe_msg = html.escape(str(message))
        gif_uri = _manta_gif_data_uri()
        if gif_uri:
            visual_html = (
                f"<img src='{gif_uri}' alt='MANTA working' "
                "style='width:92%; max-width:1050px; height:auto; "
                "object-fit:contain; image-rendering:auto;'/>"
            )
        else:
            visual_html = (
                "<span class='delfin-busy' style='width:42px; height:42px; "
                "border-width:4px;'></span>"
            )
        payload = (
            "<div class='manta-busy-stage' style='display:flex; "
            "flex-direction:column; align-items:center; justify-content:center; "
            f"width:100%; min-height:{viewer_height - 6}px; gap:20px; "
            "background:#fcfcfc;'>"
            f"{visual_html}"
            "<div class='manta-busy-caption' style='display:flex; "
            "align-items:center; justify-content:center; gap:8px; max-width:84%; "
            "font-family:monospace; font-size:13px; line-height:1.4; "
            "color:#1976d2; text-align:center;'>"
            "<span class='delfin-busy' style='vertical-align:middle;'></span>"
            f"<span>{safe_msg}</span></div></div>"
        )
        show_output(({
            'output_type': 'display_data',
            'data': {'text/html': payload},
            'metadata': {},
        },))
        _set_manip_toolbar_enabled(False)

    def _apply_smiles_conversion_result(task_id, *, quick, cleaned_data, result):
        if task_id != state['smiles_task_id']:
            return

        _set_smiles_conversion_busy(False)
        if quick:
            xyz_string = result.get('xyz_string')
            preview_items = result.get('preview_items') or []
            error = result.get('error')
            if error or not xyz_string:
                _replace_mol_output_text(f'Error: {error or "Conversion failed"}')
                return
            state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': xyz_string}
            state['charge_is_the_users'] = False
            state['charge_filled_for'] = None
            _offer_isomers(
                [(xyz_string, result['num_atoms'], 'quick')] + preview_items,
                quick=True)
            return

        isomers = result.get('isomers') or []
        error = result.get('error')
        if error or not isomers:
            _replace_mol_output_text(
                f'SMILES: {cleaned_data}',
                f'Error: {error or "No isomers generated"}',
            )
            state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
            state['isomers'] = []
            isomer_nav_row.layout.display = 'none'
            return

        state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': isomers[0][0]}
        state['charge_is_the_users'] = False
        state['charge_filled_for'] = None
        _offer_isomers(isomers)


    def on_xyz_copy(button):
        content = state['current_xyz_for_copy'].get('content')
        if not content:
            xyz_copy_status.value = '<span style="color:#d32f2f;">No XYZ to copy</span>'
            return
        text_payload = json.dumps(str(content))
        js_code = (
            "(function(){"
            f"const text={text_payload};"
            "function _manualPrompt(){"
            "try{window.prompt('Copy to clipboard (Cmd+C/Ctrl+C, Enter):', text);}catch(_e){}"
            "}"
            "function _legacyCopy(){"
            "try{"
            "const ta=document.createElement('textarea');"
            "ta.value=text;"
            "ta.setAttribute('readonly','readonly');"
            "ta.style.position='fixed';"
            "ta.style.top='-1000px';"
            "ta.style.left='-1000px';"
            "ta.style.opacity='0';"
            "document.body.appendChild(ta);"
            "ta.focus();"
            "ta.select();"
            "ta.setSelectionRange(0, ta.value.length);"
            "const ok=document.execCommand('copy');"
            "document.body.removeChild(ta);"
            "return !!ok;"
            "}catch(_e){return false;}"
            "}"
            "if(navigator.clipboard && navigator.clipboard.writeText){"
            "navigator.clipboard.writeText(text).catch(function(){"
            "if(!_legacyCopy()) _manualPrompt();"
            "});"
            "}else{"
            "if(!_legacyCopy()) _manualPrompt();"
            "}"
            "})();"
        )
        ctx.run_js(js_code)
        xyz_copy_status.value = '<span style="color:#388e3c;">Copied to clipboard</span>'

    #: How the switches read on an editor that has just been built, taken from
    #: the widgets themselves rather than written down a second time. A
    #: structure the editor has not seen starts here -- Settle on, the rest
    #: off, charge zero, multiplicity one. The charge is why this covers more
    #: than the switches: carrying a cation's charge into the next block is a
    #: wrong answer waiting to happen.
    def _structure_controls():
        """The switches that belong to the structure, not to the editor.

        Charge and multiplicity are the clearest case: a cation in one block
        and a neutral molecule in the next are two different calculations, and
        carrying one over to the other is a wrong answer waiting to happen. The
        same goes for the method, the solvent and whether a field is running --
        switching Dynamik Opt on for one structure says nothing about the next.

        Strength and Mouse are not in here. They are how the editor feels under
        the hand, and that does not change with the molecule.
        """
        return (submit_relax_btn, submit_settle_btn, submit_select_btn,
                submit_manip_btn, submit_draw_btn, submit_dyn_bonds_btn,
                submit_ff_dd, submit_gfn_charge, submit_gfn_mult,
                submit_gfn_autospin, submit_gfn_solvent, submit_gfn_solv_model,
                # The thermal budget is not in here.  T and the timescale are
                # how the user is working, like the method and the solvent,
                # and they do not stop being that because another structure
                # arrived.  The anchor is the part that belongs to one
                # structure, and it is an energy rather than a control -- it
                # is dropped in structure_changed, where the molecule it was
                # measured on stops being the one on screen.
                submit_hold_mode)

    def _apply_controls(values):
        # Same reason as reset_controls: giving a structure its switches back
        # is the editor's doing, not the user's.
        state['charge_filling'] = True
        try:
            for widget, value in zip(_structure_controls(), values or ()):
                try:
                    if widget.value != value:
                        widget.value = value
                except Exception:
                    # A list that no longer offers what it offered then.
                    pass
        finally:
            state['charge_filling'] = False
        # And then whatever the method on screen does not have, taken away
        # again.  A restore writes the switches one by one, and a memory made
        # under one method can hand a switch to another that has no such
        # thing: the method box may already hold the value it is being given,
        # so its handler never fires and nothing else would notice.
        _refresh_method_controls()

    def _controls_a_new_structure_resets():
        """What an unseen structure starts over with.

        Not everything the memory holds. The charge and the multiplicity are
        the calculation's own quantities and belong to one molecule: carrying a
        cation's charge into the next block is a wrong answer waiting to
        happen. The method, the solvent and its model are how the user is
        working, and they do not stop being that because another structure
        arrived -- picking GFN2 and having it fall back to UFF on every paste
        is its own kind of wrong.
        """
        return (submit_relax_btn, submit_settle_btn, submit_select_btn,
                submit_manip_btn, submit_draw_btn, submit_dyn_bonds_btn,
                submit_gfn_charge, submit_gfn_mult, submit_gfn_autospin,
                submit_hold_mode)

    _CONTROL_START = {id(w): w.value for w in _structure_controls()}

    def remember_structure():
        """What the editor knows about the structure it is on, to be given back.

        The bonding above all: it is perceived once from the structure as it
        arrived and kept, so that dragging an atom away from its neighbour does
        not decide the bond was never there. Step to another structure and back
        without this, and the coordinates are read afresh -- with the atom where
        it was dragged to, and no bond to it.
        """
        saved = {key: state.get(key) for key in STRUCTURE_MEMORY_KEYS}
        saved['controls'] = [w.value for w in _structure_controls()]
        return saved

    def restore_structure(saved):
        """Give back what was put aside, or start clean for an unseen one.

        The switches with it: a structure comes back the way it was left, and
        one the editor has not seen starts on the defaults.
        """
        for key in STRUCTURE_MEMORY_KEYS:
            if saved and key in saved:
                state[key] = saved[key]
            else:
                state.pop(key, None)
        # The shapes the rest of the editor reads without asking.
        for key, empty in (('perceived', None), ('poly_applied', None),
                           ('poly_metal', None), ('constraints', []),
                           ('scan_legs', []),
                           ('bond_edits', {}), ('hand_bonds', {}),
                           ('hyb_overrides', {}), ('history', []),
                           ('structure_undo', [])):
            state.setdefault(key, empty)
        if saved and saved.get('controls'):
            _apply_controls(saved['controls'])
        else:
            reset_controls()
        _refresh_constraints()

    def structure_changed():
        """A different structure is on screen now.

        The thermal anchor goes with it: an energy measured on one molecule
        says nothing about another, and a budget quietly carried over would
        report a difference between two unrelated numbers as though it were a
        distortion.

        A live force field is a set of parameters worked out for one molecule.
        Stepping from one block to another left the previous one's running
        under the new structure -- with the wrong number of atoms, even, so
        Dynamik Opt pulled at a molecule it had never been told about.

        The parameters are worked out again here, and again when the page says
        the new viewer is ready: a viewer appearing clears the parameters it
        finds, so whichever of the two comes second is the one that counts.
        """
        if state.get('thermal_e0') is not None:
            state['thermal_e0'] = None
            # And what the last structure's drag went through, which is not a
            # statement about this one.
            state.pop('thermal_peak', None)
            state.pop('thermal_good_peak', None)
            if submit_thermal_btn.value:
                _set_mol_status(
                    'A different structure: the thermal budget has no anchor '
                    'here. Press Set here to measure from it.')
        if submit_relax_btn.value:
            _enable_live_forcefield()

    def reset_controls():
        """Back to how the editor starts, for a structure it has not seen.

        A live force field belongs to the molecule its parameters were worked
        out for, and a mode belongs to the structure it was picked on. Leaving
        Dynamik Opt running across a change of structure means the next one is
        relaxed by the last one's parameters, and leaving Manipulate on means
        the first click lands on an atom nobody meant to move.

        To the defaults, not to off -- switching Settle off here took away
        something the editor is supposed to start with.
        """
        # Not the user's doing, so it must not be recorded as the user's.
        # This writes the charge, which fires the observer that remembers a
        # number the user typed -- so a reset the editor performed marked the
        # charge as theirs, and from then on it was never read off a SMILES
        # again.  Measured: an acetate then ran at zero.
        state['charge_filling'] = True
        try:
            for widget in _controls_a_new_structure_resets():
                start = _CONTROL_START.get(id(widget))
                try:
                    if widget.value != start:
                        widget.value = start
                except Exception:
                    pass
        finally:
            state['charge_filling'] = False

    def _offer_isomers(isomers, quick=False, show=True):
        """Every structure a conversion produced, to wherever they belong.

        A tab that keeps more than one -- the ORCA Builder, with its named
        blocks -- takes the whole set and steps through them itself, so the
        isomer stepper here stays out of the way. A tab that keeps one gets
        the first shown and the stepper to walk the rest, which is what the
        Submit tab has always done.

        *quick* says this came from the quick conversion, which is the one
        that answers with a structure rather than with a set to choose from.
        A tab may want plain coordinates for that one and blocks for the rest.
        """
        isomers = list(isomers or [])
        state['isomers'] = isomers
        if isomers and offer_structures is not None and offer_structures(
                isomers, quick):
            state['isomer_index'] = 0
            isomer_nav_row.layout.display = 'none'
            # Say the conversion is over. The other way out of here shows a
            # structure, and showing one clears the line on the way; this one
            # hands them all to the tab and never went past that, so the
            # "Converting SMILES..." spinner sat there for good.
            _clear_mol_status()
            return
        if show:
            _show_isomer_at_index(state.get('isomer_index', 0)
                                  if len(isomers) > 1 else 0)

    def _show_isomer_at_index(index):
        isomers = state['isomers']
        if not isomers:
            return
        index = index % len(isomers)
        state['isomer_index'] = index
        xyz_string, num_atoms, label = isomers[index]

        # Update navigation label and visibility
        if len(isomers) > 1:
            display_label = label or f'Isomer {index + 1}'
            isomer_label.value = (
                f'<span style="font-size:13px;">'
                f'{display_label} ({index + 1}/{len(isomers)})</span>'
            )
            isomer_nav_row.layout.display = ''
        else:
            isomer_nav_row.layout.display = 'none'

        xyz_data = f'{num_atoms}\nIsomer: {label}\n{xyz_string}'
        _replace_mol_output_view(xyz_data)

        # Update copy state
        state['current_xyz_for_copy'] = {'content': xyz_data}
        xyz_copy_btn.disabled = False
        xyz_copy_status.value = '<span style="color:#388e3c;">XYZ ready to copy</span>'

        # Keep converted_xyz_cache in sync for submit
        state['converted_xyz_cache']['xyz'] = xyz_string

        # Into the box without the tab drawing it again -- it is on screen
        # already. Unobserving the tab's own handler is not open to the editor:
        # it is handed a way to ask for a redraw, not the function the tab
        # registered, and taking the wrong one off the box throws. A flag both
        # hosts honour costs nothing and works for either.
        state['editor_quiet'] = True
        try:
            coords_widget.value = (
                f'{num_atoms}\nConverted from SMILES (isomer: {label})\n{xyz_string}')
        finally:
            state['editor_quiet'] = False
        # Another isomer is another structure: a running field belongs to the
        # one its parameters were worked out for.
        structure_changed()

    def _isomer_stepped_to(index):
        """Show another of the structures a conversion produced.

        And start its history: a different isomer is a different molecule, and
        the entries behind it belong to the one that was on screen a moment
        ago. Stepping is quiet -- it draws the structure itself rather than
        letting the write do it -- so the host's own "this is a structure I
        have not seen" path never runs, and Undo after a step put the previous
        isomer's geometry over this one. They need not even have the same
        number of atoms.

        Only the stepper does this. An optimisation of every frame hands its
        results back through the same shower, and there the entry from before
        the run is exactly what Undo has to reach.
        """
        _show_isomer_at_index(index)
        state['history'] = []
        _remember('the structure as it was loaded')
        # And what Reset goes back to, for the same reason.
        state['pristine_coords'] = coords_widget.value

    def handle_isomer_prev(button):
        if state['isomers']:
            _isomer_stepped_to(state['isomer_index'] - 1)

    def handle_isomer_next(button):
        if state['isomers']:
            _isomer_stepped_to(state['isomer_index'] + 1)

    def _start_smiles_conversion(*, apply_uff: bool, quick: bool, rank: bool = False,
                                 quality_mode=None, seeds_override=None,
                                 max_isomers=None, opt_topn=None, construction=None,
                                 method="gfn2", num_confs=None, collapse=None, spin="auto",
                                 deterministic: bool = True):
        # What is in the box wins whenever it is a SMILES.
        #
        # The quick conversion remembers the SMILES it last built, so that
        # pressing it again rolls another embedding of the same molecule --
        # by then the box holds coordinates and has nothing to offer. But it
        # took the remembered one even when the box held a different SMILES:
        # draw something new in Ketcher, hand it back with TO SMILES, press
        # convert, and the structure that came out was the one before it.
        typed = (read_input() or '').strip()
        cached_smiles = state['converted_xyz_cache'].get('smiles') if quick else None
        if typed and clean_input_data(typed)[1] == 'smiles':
            raw_input = typed
        else:
            raw_input = (cached_smiles or typed or '').strip()
        if not raw_input:
            _replace_mol_output_text('Please enter SMILES in the input box.')
            return

        cleaned_data, input_type = clean_input_data(raw_input)
        if input_type != 'smiles':
            _replace_mol_output_text('Please enter SMILES in the input box.')
            return

        state['smiles_task_id'] += 1
        task_id = state['smiles_task_id']
        _set_smiles_conversion_busy(True)

        # The MANTA logo animation must ALWAYS play while MANTA builds the manifold — WITH or WITHOUT
        # post-processing.  A MANTA submit is (not quick) AND has a construction preset (the convert
        # buttons pass construction=None); Rank=No / Opt<0 = manifold-only, but the logo still shows.
        if (not quick) and (construction is not None):
            _method_lbl = str(method).upper()
            _topn = -1 if opt_topn is None else int(opt_topn)
            _opt_on = _topn >= 0
            if not rank and not _opt_on:
                _caption = ('MANTA: building the complete coordination-isomer × conformer manifold '
                            '(no post-processing)...')
            else:
                _parts = []
                if rank:
                    _parts.append(f'{_method_lbl} energy-ranking')
                if _opt_on:
                    _parts.append('optimizing all' if _topn == 0 else f'optimizing top {_topn}')
                _caption = ('MANTA: building manifold + ' + ' + '.join(_parts) +
                            ' (needs xtb; takes a bit)...')
            _show_mol_busy(_caption)
        else:
            _clear_mol_output()
            if quick:
                _set_mol_status('Quick convert (single structure)...', spinner=True)
            elif apply_uff:
                _set_mol_status('Converting SMILES with UFF...', spinner=True)
            else:
                _set_mol_status('Converting SMILES (no UFF)...', spinner=True)

        def _worker():
            import os
            # MANTA "best version": derive the GFN2 charge from the SMILES, then
            # apply the SHIP-31 champion construction + GFN2-rank env for this
            # build only (snapshot + restore so the global env isn't polluted).
            _chg = 0
            if rank or construction:
                try:
                    from rdkit import Chem as _Chem
                    _m = _Chem.MolFromSmiles(cleaned_data, sanitize=False)
                    if _m is not None:
                        _chg = _Chem.GetFormalCharge(_m)
                except Exception:
                    _chg = 0
            # construction env applies ONLY for MANTA (construction set); the plain
            # convert/build-complex buttons pass construction=None -> unchanged behaviour.
            _best_env = (_manta_best_env(_chg, construction=construction, method=method,
                                         rank=rank) if construction else {})
            _saved_env = {k: os.environ.get(k) for k in _best_env}
            os.environ.update(_best_env)
            try:
                separate = _separate.has_separate_systems(cleaned_data)
                if quick and separate:
                    # A dot in a SMILES means two molecules that are not bonded
                    # to each other, and a converter handed both at once puts
                    # them in one another: measured on a
                    # hexaphenylbenzene.benzene, the benzene came out inside
                    # the other molecule, 0.877 A at the closest.  Built apart
                    # and set side by side they come out 5.1 A apart, which is
                    # a picture somebody can work in.
                    #
                    # The hapticity previews are made per part and travel with
                    # it.  They are the alternative ways a ligand can sit on
                    # its metal, so they belong to the part that has the metal;
                    # made from the whole string they would describe a molecule
                    # that is none of the frames.
                    per_part, error = [], None
                    for position, part in enumerate(
                            _separate.split_smiles(cleaned_data), start=1):
                        made, count, _m, previews, error = (
                            smiles_to_xyz_quick_with_previews(part))
                        if error or not made:
                            error = (f'part {position} could not be built: '
                                     f'{error or "nothing came back"}')
                            break
                        per_part.append([(made, count, 'quick')]
                                        + list(previews or []))
                    frames = ([] if error
                              else _separate.combine_isomers(per_part))
                    result = {
                        'error': error,
                        'xyz_string': frames[0][0] if frames else None,
                        'num_atoms': frames[0][1] if frames else 0,
                        'preview_items': frames[1:],
                        'separate_parts': len(per_part),
                    }
                elif quick:
                    xyz_string, num_atoms, _method, preview_items, error = (
                        smiles_to_xyz_quick_with_previews(cleaned_data)
                    )
                    result = {
                        'error': error,
                        'xyz_string': xyz_string,
                        'num_atoms': num_atoms,
                        'preview_items': preview_items,
                    }
                else:
                    # Interactive metal-complex conversion should prioritize
                    # isomer diversity over strict reproducibility.
                    _iso_kwargs = dict(
                        apply_uff=apply_uff,
                        collapse_label_variants=(bool(collapse) if collapse is not None else False),
                        include_binding_mode_isomers=True,
                        deterministic=deterministic,
                    )
                    # user-exposed completeness/speed switches (MANTA settings row);
                    # None -> library default. max_isomers None/0 -> COMPLETE (no cut).
                    if quality_mode:
                        _iso_kwargs["quality_mode"] = quality_mode
                    if seeds_override:
                        _iso_kwargs["seeds_override"] = int(seeds_override)
                    if max_isomers:
                        _iso_kwargs["max_isomers"] = int(max_isomers)
                    if num_confs:
                        _iso_kwargs["num_confs"] = int(num_confs)
                    if separate:
                        # Every part gets its own manifold, and the part with
                        # the most arrangements drives the navigation -- a
                        # counter-ion with one form does not multiply the
                        # complex's twelve into twelve of itself.
                        per_part, error = [], None
                        for position, part in enumerate(
                                _separate.split_smiles(cleaned_data), start=1):
                            made, error = smiles_to_xyz_isomers(
                                part, **_iso_kwargs)
                            if error or not made:
                                error = (f'part {position} could not be built: '
                                         f'{error or "nothing came back"}')
                                break
                            # Its own hapticity previews, from its own SMILES:
                            # the ways this ligand can sit on this metal, which
                            # is a question about this part and no other.
                            made = append_hapto_previews_to_isomers(
                                made, part, include_quick=apply_uff)
                            per_part.append(made)
                        isomers = ([] if error
                                   else _separate.combine_isomers(per_part))
                    else:
                        isomers, error = smiles_to_xyz_isomers(
                            cleaned_data, **_iso_kwargs)
                    # For a split input each part has had its own already,
                    # from its own SMILES.  Handed the whole string here they
                    # would describe a molecule that is none of the frames.
                    if not error and isomers and not separate:
                        isomers = append_hapto_previews_to_isomers(
                            isomers,
                            cleaned_data,
                            include_quick=apply_uff,
                        )
                    if rank and not error and isomers:
                        # RANK (opt-in): reorder best-first by xtb SINGLE-POINT energy.
                        # Geometry UNCHANGED — the emitted structures stay byte-identical to
                        # construction, only their order changes.
                        isomers = _manta_rank_only(isomers, _chg, method=method, spin=spin)
                    if (opt_topn is not None and int(opt_topn) >= 0
                            and not error and isomers):
                        # OPT (opt-in, independent of Rank): xtb geometry-optimise the top-N
                        # (0 = the whole manifold) for best-possible final geometry.
                        isomers = _manta_opt_top(isomers, _chg, topn=opt_topn, method=method, spin=spin)
                    result = {'error': error, 'isomers': isomers}
            except Exception as exc:
                result = {'error': str(exc)}
            finally:
                for _k, _v in _saved_env.items():
                    if _v is None:
                        os.environ.pop(_k, None)
                    else:
                        os.environ[_k] = _v

            schedule_ui_update(
                _apply_smiles_conversion_result,
                task_id,
                quick=quick,
                cleaned_data=cleaned_data,
                result=result,
            )

        _start_background(_worker, 'The conversion')

    def _convert_smiles(*, apply_uff: bool):
        _start_smiles_conversion(apply_uff=apply_uff, quick=False)

    def handle_convert_smiles(button):
        _convert_smiles(apply_uff=False)

    def handle_convert_smiles_quick(button):
        _start_smiles_conversion(apply_uff=False, quick=True)

    def handle_manta(button):
        # MANTA: full coordination-isomer manifold (with UFF cleanup) + GFN2
        # energy ranking, shown in the viewer with the existing isomer nav.
        # Read the exposed settings row (Quality/Seeds/Max-iso/Rank/Opt-top).
        _rank_sel = manta_rank.value
        # Opt dropdown -> top-N int: No = -1 (off, keep construction geometry); All = 0; Top-N = N.
        _opt_map = {'No': -1, 'Top 5': 5, 'Top 10': 10, 'Top 20': 20, 'All': 0}
        _start_smiles_conversion(
            apply_uff=True, quick=False,
            rank=(_rank_sel != 'No'),
            method=(_rank_sel if _rank_sel != 'No' else 'gfn2'),
            quality_mode=(manta_quality.value or None),
            # seeds field = preset value (transparent) unless the user edited it ->
            # custom seed count on top of the preset's cap/templates.
            seeds_override=(int(manta_seeds.value) or None),
            # 0 -> COMPLETE manifold (never cut off); else user cap
            max_isomers=(int(manta_max_iso.value) or 100000),
            opt_topn=_opt_map.get(manta_opt.value, -1),
            spin=str(manta_spin.value),     # 'auto' (scan) or fixed multiplicity (1/2/3/...)
            # Determinism toggle (default On) -> byte-identical, IDENTICAL to the CLI
            # and the development loop (ship = validate).  Off = non-deterministic embed.
            deterministic=(manta_det.value == 'On'),
            # construction always champion + power-user knobs CLI-only -> pinned here
            **_MANTA_DASH_DEFAULTS,
        )

    def handle_convert_smiles_uff(button):
        _convert_smiles(apply_uff=True)


    def _clean_xyz_block(raw_xyz):
        text = (raw_xyz or '').strip()
        if not text:
            return ''
        lines = text.splitlines()
        if len(lines) >= 3:
            try:
                int(lines[0].strip())
                return '\n'.join(lines[2:]).strip()
            except ValueError:
                pass
        return '\n'.join(lines).strip()

    xyz_copy_btn.on_click(on_xyz_copy)

    convert_smiles_button.on_click(handle_convert_smiles)
    convert_smiles_quick_button.on_click(handle_convert_smiles_quick)
    convert_smiles_uff_button.on_click(handle_convert_smiles_uff)

    manta_button.on_click(handle_manta)

    isomer_prev_btn.on_click(handle_isomer_prev)
    isomer_next_btn.on_click(handle_isomer_next)
    # ---- drawing one, in two dimensions -----------------------------
    # wanted rather than carried in the repository -- see dashboard/ketcher.py.
    submit_draw_open_btn = widgets.ToggleButton(
        value=False, description='DRAW', icon='pencil', button_style='info',
        tooltip=('Draw the structure in Ketcher and hand it back as a SMILES. '
                 'The editor is fetched the first time it is opened.'),
        layout=widgets.Layout(width='110px'),
    )
    submit_draw_get_btn = widgets.Button(
        description='TO SMILES', icon='arrow-down', button_style='success',
        tooltip='Put what is drawn into the input box above, as a SMILES.',
        layout=widgets.Layout(width='130px', display='none'),
    )
    submit_draw_update_btn = widgets.Button(
        description='Update', icon='refresh',
        tooltip='Fetch the newest published Ketcher.',
        layout=widgets.Layout(width='100px', display='none'),
    )
    submit_draw_frame = widgets.HTML(value='', layout=widgets.Layout(
        width='100%', display='none'))
    submit_draw_frame.add_class('submit-ketcher-frame')
    # Also stamped with this editor's scope further down, once it has one:
    # both tabs put their drawing frame outside the scope container, so the
    # class on the element is the only way to tell one editor's frame from
    # the other's -- and TO SMILES asked the page for "the" frame, found the
    # Submit tab's, was told there was no editor open in it, and answered into
    # that tab's channel. The Builder's "Reading the drawing..." never ended.
    # What the editor hands back, on the same kind of channel the viewer uses:
    # a widget value is ordered and survives a background thread, where a
    # script sent through run_js can be replaced before the page has run it.
    submit_draw_sync = widgets.Textarea(
        value='', layout=widgets.Layout(display='none'))
    submit_draw_sync.add_class('submit-ketcher-sync')
    # Which editor this drawing belongs to. Both tabs keep the frame outside
    # the scope container, so the class has to travel on the element itself.
    submit_draw_frame.add_class(submit_scope_id)
    submit_draw_sync.add_class(submit_scope_id)

    # -- drawing the structure ------------------------------------------
    def _draw_frame_html(url):
        """The editor itself, in a frame of its own.

        A frame rather than the page: Ketcher is a React application that owns
        its own globals, and the dashboard is another one.  Same origin, so the
        page may reach in and ask it for the drawing -- across origins there
        would be nothing to ask with, because Ketcher speaks no messages.
        """
        return (
            "<iframe src='" + html.escape(url, quote=True) + "' "
            "style='width:100%; height:560px; border:1px solid #d0d0d0; "
            "border-radius:6px; background:#fff;' "
            "title='Ketcher'></iframe>"
        )

    def _refresh_ketcher_controls():
        # Named for what it is.  This was called _refresh_draw_controls, which
        # is also the name of the one that shows the element dropdown for
        # drawing atoms in the viewer -- the later definition replaced the
        # earlier, so switching the viewer's Draw on stopped offering an
        # element to draw with.
        drawn = bool(submit_draw_open_btn.value)
        ready = _ketcher.is_installed()
        submit_draw_frame.layout.display = '' if (drawn and ready) else 'none'
        submit_draw_get_btn.layout.display = '' if (drawn and ready) else 'none'
        submit_draw_update_btn.layout.display = '' if (drawn and ready) else 'none'

    #: Keep the pane where it is while the editor loads.
    #:
    #: Opening DRAW does not move anything: the frame appears below the button
    #: and the button stays put.  Two to three seconds later Ketcher finishes
    #: loading inside the frame and takes the focus, and the browser scrolls
    #: the frame into view.  Measured: at +1.2 s the pane was still at 0 with
    #: the button 538 px down the viewport; at +3 s the pane was at 654 and the
    #: button at -116, off the top of the screen.
    #:
    #: So the first scroll within six seconds of opening is undone, and then
    #: the hold is gone -- one correction, not a clamp.  A clamp was tried and
    #: is worse than the jump: it also swallowed the user's own scrolling, and
    #: whether their wheel released it depended on what the pointer happened to
    #: be over. At most one deliberate scroll in those six seconds is undone,
    #: and the next one stands.
    _KETCHER_SCROLL_HOLD_JS = """
    (function() {
        var found = null;
        var buttons = document.querySelectorAll('button');
        for (var i = 0; i < buttons.length; i++) {
            if ((buttons[i].textContent || '').trim() === 'DRAW') {
                found = buttons[i];
                break;
            }
        }
        if (!found) return;
        var pane = found.parentElement;
        while (pane && pane !== document.body) {
            var how = window.getComputedStyle(pane);
            if (/(auto|scroll)/.test(how.overflowY)
                && pane.scrollHeight > pane.clientHeight + 4) break;
            pane = pane.parentElement;
        }
        if (!pane || pane === document.body) pane = document.scrollingElement;
        if (!pane) return;
        if (window.__delfinDrawScrollRelease) window.__delfinDrawScrollRelease();
        var keep = pane.scrollTop;
        var timer = window.setTimeout(release, 6000);
        function release() {
            window.clearTimeout(timer);
            pane.removeEventListener('scroll', onScroll);
            window.__delfinDrawScrollRelease = null;
        }
        function onScroll() {
            release();
            if (pane.scrollTop !== keep) pane.scrollTop = keep;
        }
        pane.addEventListener('scroll', onScroll);
        window.__delfinDrawScrollRelease = release;
    })();
        """

    def on_submit_draw_open(change):
        if change.get('name') != 'value':
            return
        if not submit_draw_open_btn.value:
            # Folded away, not thrown away.  Emptying the frame ends the
            # application inside it, and with it whatever had been drawn --
            # so folding it shut and open again lost the structure, which is
            # the opposite of what folding something away means.
            _refresh_ketcher_controls()
            return
        url = _ketcher.app_url()
        if url:
            version = _ketcher.installed_version()
            # Built once.  Setting the same frame again reloads it, and a
            # reloaded editor is an empty one.
            if url not in (submit_draw_frame.value or ''):
                submit_draw_frame.value = _draw_frame_html(url)
            _refresh_ketcher_controls()
            _run_manip_js(_KETCHER_SCROLL_HOLD_JS)
            _set_mol_status(
                f'Ketcher {version}: draw the structure, then press TO SMILES '
                'to put it in the input box.')
            return
        # Not there yet.  Offered rather than fetched: it is thirty-odd
        # megabytes, and on a machine without a network it is a wait that ends
        # in nothing.
        _refresh_ketcher_controls()
        if _ketcher.app_directory() is None:
            submit_draw_open_btn.value = False
            _set_mol_status(
                'The drawing editor needs a directory the browser can load it '
                'from, and this dashboard is not serving one.')
            return
        state['draw_installing'] = True
        _set_mol_status('Looking up the newest Ketcher...', spinner=True)

        def _work():
            newest = _ketcher.latest_release()

            def _ask():
                state['draw_installing'] = False
                if not newest['ok']:
                    submit_draw_open_btn.value = False
                    _set_mol_status(newest['status'])
                    return
                state['draw_offer'] = newest
                submit_draw_get_btn.layout.display = 'none'
                submit_draw_update_btn.description = 'Fetch it'
                submit_draw_update_btn.button_style = 'warning'
                submit_draw_update_btn.layout.display = ''
                _set_mol_status(
                    f'Ketcher {newest["version"]} is not here yet: '
                    f'{newest["size"] / 1e6:.0f} MB from '
                    'github.com/epam/ketcher, unpacked to about 32 MB beside '
                    'the dashboard. Press Fetch it to get it.')

            schedule_ui_update(_ask)

        _start_background(_work, 'The lookup of the newest Ketcher',
                          guards={'draw_installing': False})

    def on_submit_draw_update(_button=None):
        """Fetch the newest build -- the first time, or over an older one."""
        if state.get('draw_installing'):
            return
        state['draw_installing'] = True
        submit_draw_update_btn.layout.display = 'none'
        _set_mol_status('Fetching Ketcher...', spinner=True)

        def _work():
            def _line(text):
                schedule_ui_update(_set_mol_status, f'Ketcher: {text}',
                                    spinner=True)

            outcome = _ketcher.install(on_line=_line)

            def _done():
                state['draw_installing'] = False
                state['draw_offer'] = None
                submit_draw_update_btn.description = 'Update'
                submit_draw_update_btn.button_style = ''
                if not outcome['ok']:
                    submit_draw_open_btn.value = False
                    _refresh_ketcher_controls()
                    _set_mol_status(outcome['status'])
                    return
                url = _ketcher.app_url()
                submit_draw_frame.value = _draw_frame_html(url) if url else ''
                if not submit_draw_open_btn.value:
                    submit_draw_open_btn.value = True
                _refresh_ketcher_controls()
                _set_mol_status(outcome['status'],
                                'Draw the structure, then press TO SMILES.')

            schedule_ui_update(_done)

        _start_background(_work, 'The Ketcher download',
                          guards={'draw_installing': False})

    def on_submit_draw_get(_button=None):
        """Ask the editor for what has been drawn.

        The molfile, not the SMILES the editor could write itself: everything
        downstream here reads structures with RDKit, and a SMILES RDKit wrote
        is one RDKit will certainly read back.
        """
        _set_mol_status('Reading the drawing...')
        _begin_round_trip('draw_get', 'The drawing coming back')
        ctx.run_js(
            "(function(){\n"
            f"  var scope={json.dumps(submit_scope_id)};\n"
            "  var box=document.querySelector('.submit-ketcher-sync.'+scope);\n"
            "  var input=box&&box.querySelector('textarea, input');\n"
            "  function hand(text){\n"
            "    if(!input) return;\n"
            "    var proto=(input.tagName==='TEXTAREA')\n"
            "      ? window.HTMLTextAreaElement.prototype\n"
            "      : window.HTMLInputElement.prototype;\n"
            "    var setter=Object.getOwnPropertyDescriptor(proto,'value');\n"
            "    /* A serial in front, so drawing the same thing twice reads\n"
            "       as two answers rather than as one that never came. */\n"
            "    var line=(Date.now())+'\\n'+text;\n"
            "    if(setter&&setter.set) setter.set.call(input,line);\n"
            "    else input.value=line;\n"
            "    input.dispatchEvent(new Event('input',{bubbles:true}));\n"
            "    input.dispatchEvent(new Event('change',{bubbles:true}));\n"
            "  }\n"
            "  var host=document.querySelector('.submit-ketcher-frame.'+scope);\n"
            "  var frame=host&&host.querySelector('iframe');\n"
            "  var api=null;\n"
            "  try{ api=frame&&frame.contentWindow&&frame.contentWindow.ketcher; }\n"
            "  catch(e){ api=null; }\n"
            "  if(!api){ hand('!no-editor'); return; }\n"
            "  try{\n"
            "    Promise.resolve(api.getMolfile()).then(function(mol){\n"
            "      hand(mol||''); }, function(err){ hand('!'+err); });\n"
            "  }catch(e){ hand('!'+e); }\n"
            "})();"
        )

    def on_submit_draw_sync(change):
        """What the editor handed back."""
        if change.get('name') != 'value':
            return
        raw = submit_draw_sync.value or ''
        if '\n' not in raw:
            return
        _end_round_trip('draw_get')
        molfile = raw.split('\n', 1)[1]
        if molfile.startswith('!'):
            trouble = molfile[1:]
            _set_mol_status(
                'The editor is not open yet, so there is nothing to read.'
                if trouble == 'no-editor' else
                f'The drawing could not be read: {trouble}')
            return
        outcome = _ketcher.smiles_from_molfile(molfile)
        if not outcome['ok']:
            _set_mol_status(outcome['status'])
            return
        # Into the box the rest of the tab reads, so Convert, Build and every
        # other button downstream sees it exactly as a typed SMILES.
        write_input(outcome['smiles'])
        _set_mol_status(f'Drawn: {outcome["smiles"]}',
                        'It is in the input box -- Convert turns it into '
                        'coordinates.')

    submit_draw_open_btn.observe(on_submit_draw_open, names='value')
    submit_draw_get_btn.on_click(on_submit_draw_get)
    submit_draw_update_btn.on_click(on_submit_draw_update)
    submit_draw_sync.observe(on_submit_draw_sync, names='value')
    return Editor(locals())
