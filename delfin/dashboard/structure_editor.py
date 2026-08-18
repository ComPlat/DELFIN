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

from . import gfn_optimize as _gfn
from . import ketcher as _ketcher
from . import mopac_optimize as _mopac
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
_EDITOR_COMMENTS = (
    'optimised in delfin viewer',
    'edited in delfin viewer',
    'settled with ',
    'stopped at the frame on screen',
    'delfin drag-end',
    'delfin drag-follow',
    'from the delfin viewer',
)


def _is_editor_comment(line):
    """Whether that comment line is this editor's own claim about a geometry."""
    text = str(line or '').strip().lower()
    return any(text.startswith(one) for one in _EDITOR_COMMENTS)


#: Boltzmann, Planck and the gas constant, in the units this file speaks:
#: kcal/mol for energies, kelvin for temperature, seconds for time.
_BOLTZMANN_SI = 1.380649e-23          # J/K
_PLANCK_SI = 6.62607015e-34           # J s
_GAS_CONSTANT = 1.987204259e-3        # kcal/(mol K)


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
)


def _atom_numbers_js():
    """The layer itself: ``window.__delfinAtomNumbers``.

    ``set(viewer, on, scale)`` switches the numbers on or off, ``refresh``
    brings them back onto the atoms after those have moved, ``setScale``
    resizes what is already there.
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
        "  function build(v,scale){\n"
        "    if(!v||typeof v.addLabel!=='function')return 0;\n"
        "    clear(v);\n"
        "    if(scale!=null&&isFinite(+scale))v.__delfinLabelScale=+scale;\n"
        # alignment:center anchors the text box on its centre, so the number
        # stays on the atom centre at every zoom (default corner-anchoring
        # drifts aside as atoms shrink). Fall back to the string form if the
        # enum is unavailable.
        "    var al=(window.$3Dmol&&$3Dmol.SpriteAlignment&&$3Dmol.SpriteAlignment.center)\n"
        "      ?$3Dmol.SpriteAlignment.center:'center';\n"
        "    var ms=modelsOf(v),L=[];\n"
        "    for(var mi=0;mi<ms.length;mi++){\n"
        "      var atoms=[];try{atoms=ms[mi].selectedAtoms({})||[];}catch(e){atoms=[];}\n"
        "      for(var i=0;i<atoms.length;i++){\n"
        # fontSize 48 is a HIGH-RES texture kept sharp; the sprite is then
        # down-scaled in refresh() so the number appears small and crisp.
        # The fourth argument tells 3Dmol not to draw a frame per label.
        "        var a=atoms[i],lab=null;\n"
        "        try{lab=v.addLabel(String(i),{position:{x:a.x,y:a.y,z:a.z},\n"
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
        "  function set(v,on,scale){\n"
        "    if(!v)return 0;\n"
        "    if(!on){\n"
        "      var had=(v.__delfinLbls||[]).length;\n"
        "      clear(v);if(had)draw(v);\n"
        "      return 0;\n"
        "    }\n"
        "    var n=build(v,scale);\n"
        "    if(n){refresh(v);draw(v);}\n"
        "    return n;\n"
        "  }\n"
        "  return {set:set,build:build,clear:clear,refresh:refresh,setScale:setScale};\n"
        "})();"
    ) % LABEL_SCALE_DEFAULT


def atom_numbers_js():
    """The layer, installed once per page and only if it is not there yet."""
    return 'if(!window.__delfinAtomNumbers){\n' + _atom_numbers_js() + '\n}'


def show_atom_numbers_js(var='viewer', on=True, scale=None):
    """Number the atoms of the viewer held in the JS variable *var*.

    The numbers are read off the model the viewer already has, so this is the
    same call whether a molecule was just rendered or has been on screen for
    an hour, and it never needs the coordinates handed to it a second time.
    """
    size = LABEL_SCALE_DEFAULT if scale is None else float(scale)
    return (
        atom_numbers_js()
        + '\nwindow.__delfinAtomNumbers.set(%s,%s,%.3f);'
        % (var, 'true' if on else 'false', size)
    )



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
                         ('batch_preview_busy', False)):
        state.setdefault(_key, _value)

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
        tooltip=('How tall the numbers are, in pixels. Type one or step it; '
                 'the numbers resize as you go.'),
        layout=widgets.Layout(width='62px', height='30px', display='none'),
    )
    submit_manip_undo_btn = widgets.Button(
        description='Undo', button_style='info', icon='undo',
        tooltip='Undo last move/rotate (Ctrl-Z)',
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
        options=[('UFF', 'uff'), ('MMFF94', 'mmff94'),
                 ('GFN-FF', 'gfnff'), ('GFN2-xTB', 'gfn2'), ('g-xTB', 'gxtb'),
                 ('PM6-D3H4', 'pm6d3h4'), ('PM6', 'pm6'), ('PM7', 'pm7')],
        value='uff',
        tooltip=(
            'What Optimise minimises with. UFF and MMFF94 run in the browser '
            'and also drive the live relaxation while you drag. GFN-FF, '
            'GFN2-xTB and g-xTB run xtb on the server, and they know about the '
            'metal where UFF guesses. g-xTB approximates wB97M-V/def2-TZVPPD '
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
            'that already runs while you drag, so it costs nothing.'
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
            'While Dynamik Opt is on: when you let go of an atom, carry on '
            'down to a minimum. Switch off to have it stop where your drag '
            'left it, so you can move something else first -- then press '
            'Optimize when the structure is what you meant.'
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
    submit_scan_way = widgets.Dropdown(
        options=[('closer', 'in'), ('further', 'out')], value='in',
        tooltip=(
            'Which way to walk. A scan stops at the next minimum, so where it '
            'ends is the chemistry rather than a number -- the direction is '
            'the one thing that cannot be read off the selection.'
        ),
        layout=widgets.Layout(width='96px', display='none'),
        disabled=True,
    )
    submit_scan_to = widgets.FloatText(
        value=0.0, step=0.1, description='',
        layout=widgets.Layout(width='84px', display='none'),
        disabled=True,
    )
    submit_scan_steps = widgets.BoundedIntText(
        value=20, min=2, max=400, step=1, description='',
        layout=widgets.Layout(width='64px', display='none'),
        disabled=True,
    )
    submit_scan_whole = widgets.Checkbox(
        value=False, description='all the way', indent=False,
        tooltip=(
            'Keep walking past the next minimum. Off, the scan stops once it '
            'is over a barrier and has settled again, because past that it is '
            'pushing into a structure rather than following a reaction.'
        ),
        layout=widgets.Layout(width='118px', display='none'),
        disabled=True,
    )
    submit_scan_run_btn = widgets.Button(
        description='Run scan', button_style='success', icon='play',
        tooltip='Walk every armed coordinate together, and say what it costs.',
        layout=widgets.Layout(width='104px', height='30px', display='none'),
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
    submit_internal_group = widgets.HBox(
        [submit_internal_label, submit_internal_value,
         submit_internal_btn, submit_hold_btn, submit_hold_mode,
         submit_scan_btn, submit_scan_to, submit_scan_steps,
         submit_scan_way, submit_scan_dd, submit_scan_del, submit_scan_whole,
         submit_scan_run_btn],
        layout=widgets.Layout(
            gap='6px', align_items='center', flex_flow='row nowrap',
            flex='0 0 auto', overflow='visible',
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
            submit_labels_btn, submit_label_size,
            submit_manip_undo_btn, submit_reset_btn,
            submit_ff_dd, submit_gfn_charge, submit_gfn_mult,
            submit_gfn_autospin, submit_gfn_solvent, submit_gfn_solv_model,
            submit_thermal_btn, submit_temperature,
            submit_thermal_relax, submit_thermal_anchor_btn,
            submit_xtb_install_btn, submit_xtb_confirm_btn,
            submit_xtb_cancel_btn,
            submit_strength_slider, submit_sens_slider, submit_play_speed,
            submit_fs_row_break,
            submit_optimize_btn, submit_optimize_all_btn,
            submit_relax_btn, submit_auto_btn, submit_settle_btn,
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

    def _set_mol_status(*lines, spinner=False):
        # Both copies always say the same thing; which one is on screen is the
        # overlay's business, not the caller's.
        rendered = [html.escape(str(line)) for line in lines if line not in (None, '')]
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
        prompt = any('enter XYZ' in str(line) or 'Load a structure' in str(line)
                     for line in lines)
        mol_status_fs.value = '' if prompt else rendered_html

    def _clear_mol_status():
        # Both copies, the way _set_mol_status writes both. Clearing only the
        # small one left the overlay saying "Quick convert (single
        # structure)..." long after the structure was on screen: the finished
        # view comes through here, and in fullscreen that stale line was the
        # only thing the user had to go by.
        mol_status.value = ''
        mol_status_fs.value = ''

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
        submit_play_speed.disabled = not enabled
        for widget in (submit_thermal_btn, submit_temperature,
                       submit_thermal_relax, submit_thermal_anchor_btn):
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
            '  var play={queue:[],at:0,started:0,last:null,seen:0,run:null};\n'
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
            '    if(play.pace) return play.pace;\n'
            '    if(n>60) return 8;\n'
            '    if(n>25) return 20;\n'
            '    if(n>10) return 35;\n'
            '    /* One answer at a time -- a hand being followed -- is drawn\n'
            '       over exactly the time the next answer takes to come, or the\n'
            '       picture sits still between them.  A relaxation arrives in\n'
            '       bursts instead, and those are paced by the rules above: a\n'
            '       burst of thirty drawn at a follow pace of a third of a\n'
            '       second each would put the picture ten seconds behind a\n'
            '       calculation that took one. */\n'
            '    if(play.follow&&play.gap&&n<=3) return play.gap;\n'
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
            ' play.run=null; return; }\n'
            '    var data=null;\n'
            '    try{ data=JSON.parse(text); }catch(e){ return; }\n'
            '    if(data&&data.halt){\n'
            '      /* The run was switched off.  Playing out the queue after\n'
            '         that is the picture carrying on without the thing it is\n'
            '         a picture of. */\n'
            '      play.queue=[]; play.seen=(data.frames||[]).length;\n'
            '      if(!play.toldStop){ play.toldStop=1;\n'
            '        /* Which frame is on screen.  Stopping keeps that one:\n'
            '           frames xtb had already computed but nobody had seen\n'
            '           are not what the user stopped at. */\n'
            '        say("stopped at frame "+(play.shown||0)); }\n'
            '      return;\n'
            '    }\n'
            '    var frames=(data&&data.frames)||[];\n'
            '    var run=(data&&data.run)||0;\n'

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
            '      if(play.queue.length){\n'
            '        /* Land the run that is ending on its last frame first.\n'
            '           Dropped there, the viewer keeps whatever it had drawn\n'
            '           while the kernel keeps the geometry it computed, and\n'
            '           the two drift apart: the next drag then hands over a\n'
            '           structure that is behind, and the relaxation nobody saw\n'
            '           is walked again -- which is every earlier drag being\n'
            '           made a second time. */\n'
            '        show(play.last,play.queue[play.queue.length-1],1);\n'
            '      }\n'
            '      play.run=run; play.seen=0; play.queue=[]; play.last=null;\n'
            '      play.shown=0; play.toldStop=0; play.complete=0;\n'
            '    }\n'
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
            '      /* How long the machine took over the last answer, which is\n'
            '         how long this one has to be drawn over. */\n'
            '      if(play.arrived){\n'
            '        var measured=Math.min(600,Math.max(24,\n'
            '          arrivedAt-play.arrived));\n'
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
            '       and moves nothing, so it is not a grab. */\n'
            '    var held=(window._submitManipStateByScope||{})[scope];\n'
            '    var drag=held&&held.drag;\n'
            '    if(!drag) return false;\n'
            '    return drag.kind==="translate"||drag.kind==="rotate"'
            '||drag.kind==="draw";\n'
            '  }\n'
            '  function heldSerials(){\n'
            '    /* The atoms the hand is moving.  Coordinates that come back\n'
            '       describe where they were when they were sent, and the\n'
            '       cursor has moved on since -- so those are the ones the\n'
            '       playback must leave alone. */\n'
            '    var st=(window._submitManipStateByScope||{})[scope];\n'
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
            '    /* ipywidgets marks a pressed toggle with mod-active.  Reading\n'
            '       it here is instant; asking the kernel costs a round trip,\n'
            '       and the picture ran on for the length of it. */\n'
            '    var holder=inScope(".submit-optimize-switch")'
            '||(document.querySelector("."+scope+".submit-optimize-switch"));\n'
            '    if(!holder) return true;\n'
            '    var btn=(holder.tagName==="BUTTON")?holder'
            ':holder.querySelector("button");\n'
            '    if(!btn) return true;\n'
            '    return btn.classList.contains("mod-active");\n'
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
            '      send(held?"gfngrab":"gfnfree",\n'
            '           held?String(play.shown||0):"");\n'
            '    }\n'
            '    if(play.held&&!followIsOn()){\n'
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
            '        say("stopped at frame "+(play.shown||0)); }\n'
            '    }\n'
            '    readWall();\n'
            '    read(now);\n'
            '    if(play.queue.length){\n'
            '      if(!play.started) play.started=now;\n'
            '      var ms=stepMs();\n'
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
            '        play.last=play.queue.shift(); show(null,play.last,1);\n'
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
        _gfn_new_generation()
        # What the molecule looked like before this drag: the bonding is read
        # from here, not from a frame that has already been pulled about.
        state['gfn_topology_source'] = _current_xyz()
        state['gfn_follow'] = True
        state['gfn_follow_steps'] = 0
        state['gfn_follow_frames'] = []
        # Whatever the page said about the last drag's playback is not about
        # this one.
        state['gfn_play_note'] = ''
        run = int(state.get('gfn_run', 0)) + 1
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

    def _gfn_status_lines(said=None):
        """One row, whoever is writing it.

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
        """
        line = said if said is not None else (state.get('gfn_last_status') or '')
        note = str(state.get('gfn_play_note') or '')
        if note and any(fault in note for fault in _PLAY_NOTE_FAULTS):
            return ((f'{line} Trajectory: {note}.' if line
                     else f'Trajectory: {note}.'),)
        return (line,)

    def _gfn_is_working():
        """Whether there is a calculation behind what the line is reporting."""
        return bool(state.get('gfn_follow')
                    or state.get('optimize_run') is not None)

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
        """
        if not state.get('gfn_follow') and not _begin_gfn_follow():
            return          # not following: Relax is up, or GFN is not chosen
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
                    contacts = (
                        _gfn.contacts_holding(
                            current, holding, most=3,
                            was=state.get('thermal_was'),
                            turning=state.get('thermal_turn'),
                            holding=state.get('thermal_holding'))
                        if (submit_thermal_btn.value
                            and not _mopac.is_mopac_method(method)) else [])
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
                            cycles=(_THERMAL_FOLLOW_CYCLES if contacts
                                    else _GFN_FOLLOW_CYCLES),
                            timeout=30.0,
                            constraints=constraints + contacts, solvent=wet,
                            solvation_model=model,
                            topology=_gfn_topology_dir(current),
                        )
                    if not outcome.get('ok'):
                        note = str(outcome.get('status') or 'it did not run')
                        schedule_ui_update(
                            _set_mol_status,
                            f'The molecule stopped following: {note}')
                        return
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
                    said = (f'{label} is following the drag:{many} {steps} '
                            f'step(s), '
                            f'{(time.perf_counter() - began) * 1000:.0f} ms '
                            'each.')
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
                    priced = outcome if contacts else (
                        _gfn.optimize_with_gfn(
                            current, method, charge=charge, uhf=uhf,
                            timeout=30.0, solvent=wet, solvation_model=model,
                            topology=_gfn_topology_dir(current),
                            optimise=False,
                        ) if submit_thermal_btn.value else {})
                    state['thermal_now'] = priced.get('energy')
                    if submit_thermal_btn.value:
                        spent = _thermal_note(priced.get('energy'))
                        if spent:
                            # On the end of the line that is already there,
                            # never under it: this row stands above the viewer
                            # and a second one moves the picture while the
                            # user is aiming an atom.
                            said = f'{said} {spent}'
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
                        slipped = _gfn.largest_shift(
                            _gfn.hold_atoms_at(
                                outcome['xyz'], current, holding),
                            outcome['xyz']) if contacts else 0.0
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
                        tightest = _gfn.closest_contact(current)[0]
                        crowded = (tightest is not None
                                   and tightest < _gfn._TOO_CLOSE)
                        _thermal_slope(
                            (float(priced['energy'])
                             - float(_thermal_budget()[0] or 0.0))
                            * _HARTREE_TO_KCAL,
                            current, [int(i) for i in (holding or ())]
                        ) if priced.get('energy') is not None \
                            and _thermal_budget()[0] is not None else None
                        came_back = _thermal_wall(
                            current, priced.get('energy'), holding,
                            refuse=(slipped > _SLIP_LOOSE) or crowded)
                        if came_back is not None:
                            said = (f'{said} Past the budget, so the last '
                                    f'structure that was inside it is back.')
                            # The box as well, not only the picture.
                            #
                            # The frames go to the viewer, but the coordinate
                            # box is written from the page's own model when the
                            # hand lets go -- and that model still has the atom
                            # where the cursor left it.  So the drag sprang
                            # back on screen, the user let go, and what was
                            # kept was the torn structure after all: measured,
                            # the viewer drew 1.40 A while the box held 3.40.
                            # The rollback has to reach the thing that outlives
                            # the drag.
                            rows = [line for line
                                    in came_back.splitlines()[2:]
                                    if line.strip()]
                            if rows:
                                schedule_ui_update(
                                    _write_coords,
                                    xyz_document(
                                        rows,
                                        'Past the budget: back to the last '
                                        'structure that was inside it'),
                                    True)
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
                    state['gfn_last_status'] = said
                    schedule_ui_update(_set_mol_status,
                                       *_gfn_status_lines(said), spinner=True)
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
                    settled = (came_back if came_back is not None else
                               _gfn.hold_atoms_at(
                                   outcome['xyz'], current, holding))
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
                    payload = json.dumps({'run': state.get('gfn_follow_run'),
                                          'from': len(frames) - len(trail),
                                          'follow': 1, 'frames': trail})
                    schedule_ui_update(
                        lambda text=payload: setattr(
                            submit_gfn_frame, 'value', text))
            finally:
                state['gfn_follow_busy'] = False

        threading.Thread(target=_work, daemon=True).start()

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

    def _write_coords(text, drawn=False):
        """Put a geometry in the box, and say whether the picture has it.

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
        return state.get('thermal_e0'), ceiling

    def _thermal_note(energy):
        """Where this geometry stands against the budget, said in one line."""
        anchor, ceiling = _thermal_budget()
        if anchor is None or energy is None:
            return ''
        spent = (float(energy) - float(anchor)) * _HARTREE_TO_KCAL
        T = float(submit_temperature.value)
        if spent <= ceiling:
            return (f'{spent:+.1f} of {ceiling:.1f} kcal/mol available at '
                    f'{T:g} K.')
        # Said as a time rather than as a number: how long this would take is
        # the thing a chemist can judge, and 34.7 kcal/mol means nothing until
        # it is "longer than the age of the earth".
        return (f'{spent:+.1f} kcal/mol -- past the {ceiling:.1f} this '
                f'structure has at {T:g} K. {_thermal_wait(spent, T)}')

    def _thermal_wait(kcal, temperature):
        """How long a barrier of that height takes at that temperature.

        Down to picoseconds, not only down to seconds.  A scan that finds an
        open path is the ordinary case and it was reported as "about 4.18e-06
        s", which is a number in the wrong clothes: the answer wanted there is
        "4 microseconds", and below a picosecond there is no crossing to speak
        of anyway.
        """
        T = max(1.0, float(temperature))
        rate = ((_BOLTZMANN_SI * T / _PLANCK_SI)
                * math.exp(-max(0.0, float(kcal)) / (_GAS_CONSTANT * T)))
        if rate <= 0:
            return 'It does not happen.'
        seconds = 1.0 / rate
        for limit, unit, name in ((1e-9, 1e-12, 'ps'), (1e-6, 1e-9, 'ns'),
                                  (1e-3, 1e-6, 'us'), (1.0, 1e-3, 'ms'),
                                  (60, 1, 's'), (3600, 60, 'min'),
                                  (86400, 3600, 'h'), (3.15576e7, 86400, 'd'),
                                  (float('inf'), 3.15576e7, 'years')):
            if seconds < limit:
                return f'That is about {seconds / unit:.3g} {name}.'
        return ''

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
        if not submit_thermal_btn.value:
            return
        anchor, _ceiling = _thermal_budget()
        if anchor is None:
            return
        here = _current_xyz() or ''
        if here:
            state['thermal_good'] = here

    def _thermal_wall(xyz, energy, holding, refuse=False):
        """Keep what the budget agreed to, and take back what it did not.

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
        """
        anchor, ceiling = _thermal_budget()
        if anchor is None or energy is None:
            return None
        spent = (float(energy) - float(anchor)) * _HARTREE_TO_KCAL
        if spent <= ceiling and not refuse:
            state['thermal_good'] = xyz
            return None
        # Past it.  Hand back the last geometry that was not, if there is one;
        # a drag that was already over budget when it started has nowhere to
        # go, and then the hand simply keeps what it has and the line says so.
        return state.get('thermal_good')

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

    def _set_thermal_anchor(relax=None, note='Measuring from here'):
        """Take the energy of the structure on screen as the budget's zero.

        A single point when the structure is to be kept as it is, one
        optimisation when it is not.  Either way what is stored is an energy
        of the *chosen method*, so the budget and the drag are the same
        calculation and their difference means something.
        """
        xyz = _current_xyz()
        method = str(submit_ff_dd.value)
        if not xyz or not _gfn.is_gfn_method(method):
            return
        wants_relax = (submit_thermal_relax.value if relax is None else relax)
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
                state['thermal_for'] = _structure_fingerprint(
                    outcome.get('xyz') or xyz)
                if wants_relax:
                    lines = [line for line in outcome['xyz'].splitlines()[2:]
                             if line.strip()]
                    if lines:
                        _write_coords(xyz_document(
                            lines, 'Relaxed, and the budget measured from here'))
                _, ceiling = _thermal_budget()
                _set_mol_status(
                    f'Measuring from here. At {float(submit_temperature.value):g} K '
                    f'this structure has {ceiling:.1f} kcal/mol to spend '
                    f'within {_timescale_label()}.')

            schedule_ui_update(_done)

        threading.Thread(target=_work, daemon=True).start()

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
        state['gfn_settle_note'] = note
        _arm_gfn_settle(forced=True)

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

        threading.Thread(target=_wait, daemon=True).start()

    def _gfn_settle_now():
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
            state['gfn_run'] = int(state.get('gfn_run', 0)) + 1
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
            return bool(state.get('gfn_follow')
                        or state.get('optimize_run') is not None
                        or not (submit_settle_btn.value or _gfn_live_is_on()))

        def _push(frames):
            walked = list(frames)
            trail = walked[-8:]
            state['gfn_settle_walked'] = len(walked)
            schedule_ui_update(
                lambda t=trail, f=offset + len(walked) - len(trail): setattr(
                    submit_gfn_frame, 'value',
                    json.dumps({'run': run, 'from': f, 'follow': 1,
                                'frames': t})))

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
                lines = [line for line in outcome['xyz'].splitlines()[2:]
                         if line.strip()]
                if lines:
                    # The playback has drawn this already; the box is what Copy
                    # and Submit read, and it has to be true whether or not a
                    # frame happened to land.
                    _write_coords(xyz_document(lines, f'Settled with {label}'),
                                  drawn=True)
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

        threading.Thread(target=_work, daemon=True).start()

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

    #: How long to wait after the last change before starting again.  Letting
    #: go of an atom arrives as a burst -- the release, then the coordinates,
    #: sometimes a settled version behind them -- and starting on the first of
    #: those would launch an xtb for each one.
    _GFN_RESTART_DELAY = 0.35

    def _interrupt_gfn():
        """End the running optimisation because the structure under it changed.

        xtb is minimising a geometry that stopped existing the moment an atom
        was moved, so the run is ended rather than raced.  It is not ended the
        way the switch ends it: nothing has been stopped from where the user
        is standing, and the frame they were shown is not a result to keep --
        the optimisation is about to start again from what they have made.
        """
        token = state.get('optimize_run')
        if token is None:
            return False
        state['optimize_run'] = None
        state['optimize_interrupted'] = token
        # No halt report: "stopped at frame 12" belongs to the switch.
        state['gfn_halt_sent'] = True
        # A run number the page has never seen, carrying nothing.  It resets
        # the player, so the frames of the abandoned run cannot play out over
        # the geometry the user has just made.
        blank = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = blank
        submit_gfn_frame.value = json.dumps({'run': blank, 'frames': []})
        return True

    def _stand_down_after_interrupt():
        """The run the change interrupted is not coming back by itself.

        The switch goes back up with it. Left lit over a structure nobody is
        minimising, it says a calculation is running that is not.
        """
        state.pop('optimize_interrupted', None)
        for button in (submit_optimize_btn, submit_optimize_all_btn):
            if button.value:
                button.value = False
        _set_mol_status(
            'Stopped where your change left it. Move what else you want to, '
            'then press Optimize to go down to a minimum.')

    def _arm_gfn_restart():
        """Start the optimisation again, once the changing has stopped."""
        if state.get('optimize_interrupted') is None:
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

        def _wait():
            while True:
                left = state.get('gfn_restart_at', 0.0) - time.monotonic()
                if left <= 0:
                    break
                time.sleep(min(left, 0.05))
            state['gfn_restart_armed'] = False
            schedule_ui_update(_restart_gfn)

        threading.Thread(target=_wait, daemon=True).start()

    def _arm_gfn_minimise():
        """Go down to a minimum once the changing has stopped.

        The other half of Auto: a drag with no optimisation behind it has
        nothing to resume, so one is started.  Same wait as a restart, and the
        wait is pushed out again by every release -- moving three atoms one
        after another is one minimisation at the end of it, not three.
        """
        if not (_server_method()
                and _server_binary(submit_ff_dd.value) is not None):
            return
        state['gfn_minimise_at'] = time.monotonic() + _GFN_RESTART_DELAY
        if state.get('gfn_minimise_armed'):
            return
        state['gfn_minimise_armed'] = True

        def _wait():
            while True:
                left = state.get('gfn_minimise_at', 0.0) - time.monotonic()
                if left <= 0:
                    break
                time.sleep(min(left, 0.05))
            state['gfn_minimise_armed'] = False
            schedule_ui_update(_minimise_now)

        threading.Thread(target=_wait, daemon=True).start()

    def _minimise_now():
        if not submit_auto_btn.value:
            return                      # switched off while it was waiting
        if state.get('optimize_run') is not None:
            return                      # one is already running
        if submit_optimize_btn.value or submit_optimize_all_btn.value:
            return                      # a switch is already down
        # The switch, not the handler: it has to be seen to be on for as long
        # as it runs, and it is what the user presses to stop it again.
        submit_optimize_btn.value = True

    def _after_release():
        """What letting go of an atom leads to, and the switch that decides.

        Auto on: down to a minimum, whether or not an optimisation was running
        when the atom was picked up.  That used to be the difference between
        the same gesture finishing the structure and leaving it strained --
        a drag during a run interrupted it and the run came back, a drag with
        no run behind it got Settle's short tidy and nothing else.

        Auto off: it stops where the hand left it.  Move something else, and
        press Optimize when the structure is what you meant; that goes down to
        a minimum the way it always has.

        Only while the molecule is following the hand.  Dragging with Dynamik
        Opt off is placing an atom where you want it, and starting a
        minimisation on top of that would take it off the place you just put
        it -- which is what Settle is for, in the small.
        """
        auto = bool(submit_auto_btn.value)
        if state.get('optimize_interrupted') is not None:
            # Comes back, or stands down; either way _arm_gfn_restart decides.
            _arm_gfn_restart()
            if auto:
                return                  # a minimisation is more than a settle
            _arm_gfn_settle()
            return
        if auto and submit_relax_btn.value:
            _arm_gfn_minimise()
            return
        _arm_gfn_settle()

    def _restart_gfn():
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
                state['optimize_run'] = None      # off: the run ends itself
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
        run_id = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = run_id
        # Which frame the picture stood on belongs to the run it was reported
        # for, and the page only reports it when a hand arrives or the switch
        # goes up.  Kept across runs, a number left over from an earlier grab
        # is a plausible index into this run's path -- so an edit that
        # interrupts this one would cut it at a frame nobody ever saw, from a
        # trajectory that no longer exists.
        state.pop('gfn_shown_frame', None)

        def _push_frames(frames, final=False):
            """Hand the path over while xtb is still walking it.

            Every frame, exactly once, and each one carried twice.

            The field is one slot, not a queue: a write that lands before the
            page has read the one before it replaces it, and those frames are
            gone.  That is what the eight-frame tail was for -- it re-sent
            recent frames so a missed read still caught them -- but it was a
            *fixed* eight, and xtb makes frames far faster than the page is
            asked to look.  A benzene runs 23 cycles in a fraction of a second
            and a 149-atom chain 260; everything between two reads beyond the
            last eight was never sent at all, and what reached the viewer was
            a sample of the path rather than the path.  Measured, not argued:
            those two runs are in the round-loop measurements.

            So the window starts where the *previous* window started rather
            than where it ended.  Every frame is therefore sent in two
            consecutive writes, which is the same insurance the tail gave, and
            nothing is skipped however fast the frames arrive.  It stays
            bounded -- a write carries at most two reads' worth -- and the
            coordinates are rounded to four decimals, which is a thousandth of
            a bond length and half the JSON.

            *final* is the write at the end of the run, and it goes out even
            when it carries nothing new.  Without it the last window is the one
            window sent only once -- the run ends before another write can
            repeat it -- so a single missed read at exactly that moment leaves
            the picture short of the geometry the box holds.  That is the end
            of the path, which is the part that has to land.
            """
            played[0] = True
            walked = list(frames)
            if state.get('gfn_push_run') != run_id:
                state['gfn_push_run'] = run_id
                state['gfn_push_start'] = 0
                state['gfn_push_end'] = 0
            if not final and len(walked) <= int(state.get('gfn_push_end') or 0):
                return                      # nothing new since the last write
            start = int(state.get('gfn_push_start') or 0)
            state['gfn_push_start'] = int(state.get('gfn_push_end') or 0)
            state['gfn_push_end'] = len(walked)
            fresh = [[round(float(v), 4) for v in frame]
                     for frame in walked[start:]]

            def _write(t=fresh, first=start, last=bool(final)):
                # A run that has been replaced does not draw.  An interrupted
                # one has frames in hand when it is told to stop, and writing
                # them afterwards played the abandoned path over the structure
                # the user had just made.
                if state.get('gfn_run') != run_id:
                    return
                payload = {'run': run_id, 'from': first, 'frames': t}
                if last:
                    # Named, so this write differs from the one before it even
                    # when it carries the same frames: the field is a widget
                    # value and traitlets says nothing when a value is written
                    # again unchanged, so the repeat that is the whole point of
                    # a final write would never leave the kernel.
                    payload['final'] = 1
                submit_gfn_frame.value = json.dumps(payload)

            schedule_ui_update(_write)

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
            halted = state.get('optimize_run') is not token
            if halted and not state.get('gfn_halt_sent'):
                state['gfn_halt_sent'] = True
                schedule_ui_update(
                    lambda: setattr(submit_gfn_frame, 'value',
                                    json.dumps({'run': run_id, 'halt': 1,
                                                'frames': []})))
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
                        outcome = relax_xyz(
                            xyz,
                            max_steps=500,
                            perceived=_perception_for(xyz),
                            method=method,
                        )
                except Exception as exc:
                    failures.append(f'frame {position + 1}: {exc}')
                    results.append(item)
                    continue
                if outcome.get('ok'):
                    if outcome.get('energy') is not None and position == 0:
                        state['gfn_energy'] = float(outcome['energy'])
                        state['gfn_energy_unit'] = outcome.get('energy_unit')
                    if position == 0:
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
                        shown = state.get('gfn_shown_frame')
                        walked = outcome['frames']
                        if isinstance(shown, int) and 0 < shown <= len(walked):
                            symbols = [line.split()[0]
                                       for line in _gfn.atom_lines(xyz)]
                            frame = walked[shown - 1]
                            if len(symbols) * 3 == len(frame):
                                rows = [
                                    xyz_line(symbols[i], frame[3*i],
                                             frame[3*i+1], frame[3*i+2])
                                    for i in range(len(symbols))
                                ]
                                kept = (f'{len(rows)}\nstopped at the frame on '
                                        f'screen\n' + '\n'.join(rows) + '\n')
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
                    shown = state.get('gfn_shown_frame')
                    walked = trail[0]
                    if (isinstance(shown, int) and walked
                            and 0 < shown <= len(walked)):
                        symbols = [line.split()[0]
                                   for line in _gfn.atom_lines(single or '')]
                        frame = walked[shown - 1]
                        if symbols and len(symbols) * 3 == len(frame):
                            _write_coords(xyz_document(
                                [xyz_line(symbols[i], frame[3 * i],
                                          frame[3 * i + 1], frame[3 * i + 2])
                                 for i in range(len(symbols))],
                                'stopped where you took hold'), drawn=True)
                    return
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
                        drawn=played[0])
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
                # value with one force constant for the whole set; MOPAC fixes
                # the atoms that name it and cannot express a pull at all.
                # Read out with the wrong one, a MOPAC result would claim a
                # force constant that no MOPAC run has.
                kept = state.get('gfn_held')
                if pm:
                    said += _mopac.freeze_note(kept or {
                        'held': 0, 'pulls': 0, 'dropped': [], 'frozen': set()})
                else:
                    said += _gfn.held_note(kept or {
                        'held': 0, 'dropped': [], 'mixed': False,
                        'force': None})
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
                    _set_mol_status(*_gfn_status_lines(said), spinner=True)
                    state['optimize_carrying_on'] = True
                    schedule_ui_update(
                        lambda: on_submit_optimize(None, every_frame=every_frame))
                    return
                _set_mol_status(*_gfn_status_lines(said),
                                *(failures + unfinished)[:2])

            schedule_ui_update(_apply)

        worker = threading.Thread(target=_work, daemon=True)
        state['optimize_thread'] = worker
        worker.start()

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
        if perceived is not None and len(indices) == 1:
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
            if perceived is not None and len(indices) == 1:
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
        """
        metals = set(perceived.metal_indices or ()) if perceived else set()
        # An index the structure no longer has: the browser pushes its picks
        # after a re-render, and an edit that deleted atoms renumbers them.
        # Asking the perception about one is an IndexError, and this handler
        # runs on every click in the viewer.
        chosen = [
            i for i in indices
            if i not in metals and 0 <= i < len(perceived.symbols)
        ] if perceived else []
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
        overrides = dict(state.get('hyb_overrides') or {})
        chosen = submit_hyb_dd.value or ''
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
        """
        geometry = state.get('poly_applied')
        metal = state.get('poly_metal')
        perceived = state.get('perceived')
        turnable = False
        if geometry and metal is not None and perceived is not None:
            try:
                from .molecule_forcefield import polyhedron_vertex_classes
                donors = len(perceived.neighbours()[int(metal)])
                grouped = polyhedron_vertex_classes(donors, geometry)
                turnable = bool(grouped) and len(set(grouped[0])) > 1
            except Exception:
                turnable = False
        submit_poly_turn_btn.layout.display = '' if turnable else 'none'
        submit_poly_turn_btn.disabled = not turnable
        if not turnable:
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
                'exact one(s) are held very firmly, not held still. Optimise '
                'under a GFN method to have them met exactly.')
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

    def _remember(what):
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
        """
        history = list(state.get('history') or [])
        current = coords_widget.value
        if history and history[-1].get('coords') == current:
            # Two actions from the same picture are one step back, and the
            # step is named after the first of them: going back to that state
            # undoes everything since, and the earliest is what the user would
            # look for. Overwriting the name here lost "the structure as it
            # was loaded" the first time anything happened at all.
            pass
        else:
            history.append({
                'coords': current,
                'bond_edits': dict(state.get('bond_edits') or {}),
                'hand_bonds': dict(state.get('hand_bonds') or {}),
                'constraints': list(state.get('constraints') or []),
                'what': str(what),
            })
        # The first entry is the structure as it arrived and is never dropped:
        # it is what Reset goes back to, and a long session must not lose it.
        if len(history) > _HISTORY_LIMIT:
            history = history[:1] + history[-(_HISTORY_LIMIT - 1):]
        state['history'] = history

    def _restore(entry, note):
        """Put a remembered state back on screen."""
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
        state['constraints'] = list(entry.get('constraints') or [])
        state['perceived'] = None
        state['perceived_for'] = None
        _refresh_constraints()
        _set_mol_status(note)
        if submit_relax_btn.value or state.get('ff_bootstrap_done'):
            _enable_live_forcefield()

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
        left = 0 if at_start else len(history)
        _restore(entry, f'Took back: {entry.get("what") or "the last step"}.'
                        + (f' {left} more to go back through.' if left
                           else ' That is the structure as it was loaded.'))

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
                try:
                    state['gfn_shown_frame'] = int(str(payload).rsplit(' ', 1)[1])
                except ValueError:
                    pass
            # The same two rows the follow step writes, and the spinner with
            # them while there is a calculation behind it -- this line and that
            # one are the same message, not two taking turns.
            _set_mol_status(*_gfn_status_lines(), spinner=_gfn_is_working())
            return

        if verb == 'grabbed':
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
            try:
                state['gfn_shown_frame'] = int(str(payload).strip())
            except (TypeError, ValueError):
                pass
            if _interrupt_gfn():
                _set_mol_status('Moved while it ran; the optimisation stops '
                                'there and starts again from what you make.',
                                spinner=True)
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
            # Let go of.  Whether that goes down to a minimum, or leaves the
            # structure where the hand put it, is the Auto switch -- and it is
            # asked in one place so the answer cannot depend on whether a run
            # happened to be going when the atom was picked up.
            _after_release()
            return

        if verb == 'undo':
            _undo_structure()
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

    def _scan_arrived(path):
        """Whether the walk has crossed something and is climbing out again."""
        if len(path) < _SCAN_CLIMBING + 3:
            return False
        spent = [one[1] for one in path]
        top = max(spent)
        if top - spent[0] < _SCAN_OVER_THE_TOP:
            return False              # nothing has been crossed yet
        floor = min(spent[spent.index(top):])
        if top - floor < _SCAN_OVER_THE_TOP:
            return False              # still up there
        last = spent[-(_SCAN_CLIMBING + 1):]
        return (all(b - a > _SCAN_UPHILL for a, b in zip(last, last[1:]))
                and spent[-1] - floor > _SCAN_UPHILL)

    def _scan_legs():
        return list(state.get('scan_legs') or [])

    def _describe_leg(leg):
        symbols = []
        perceived = state.get('perceived')
        for index in leg['atoms']:
            symbol = perceived.symbols[index] if perceived else '?'
            symbols.append(f'{symbol}{index}')
        unit = 'A' if leg['kind'] == 'distance' else 'deg'
        return (f"{'-'.join(symbols)} {leg['from']:.3g} -> {leg['to']:.3g} "
                f"{unit}")

    def _refresh_scan():
        """Show the armed legs, or nothing at all when there are none."""
        legs = _scan_legs()
        showing = '' if legs else 'none'
        for widget in (submit_scan_dd, submit_scan_del, submit_scan_whole,
                       submit_scan_run_btn):
            widget.layout.display = showing
            widget.disabled = not legs
        options = [(_describe_leg(leg), str(n)) for n, leg in enumerate(legs)]
        submit_scan_dd.options = options or [('nothing armed', '')]
        if options:
            submit_scan_dd.value = options[0][1]
        # The two fields belong to the selection, not to the list, so they
        # follow whether something is selected rather than whether anything is
        # armed.
        picked = len(state.get('picked') or ())
        wanted = '' if picked in _CONSTRAINT_KINDS else 'none'
        for widget in (submit_scan_way, submit_scan_to, submit_scan_steps):
            widget.layout.display = wanted
            widget.disabled = not wanted == ''

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
        return float(value) + (far if way == 'out' else -far)

    def on_submit_scan(_button=None):
        """Arm the value the selection describes as a leg of the scan."""
        indices = list(state.get('picked') or [])
        kind = _CONSTRAINT_KINDS.get(len(indices))
        if not kind:
            _set_mol_status('Pick 2, 3 or 4 atoms before arming a scan.')
            return
        here = float(submit_internal_value.value)
        target = float(submit_scan_to.value)
        if abs(target - here) < 1e-9 or target == 0.0:
            target = _suggest_scan_target(kind, here, submit_scan_way.value)
        legs = [one for one in _scan_legs() if one['atoms'] != indices]
        leg = {'kind': kind, 'atoms': indices, 'from': here,
               'to': target, 'steps': int(submit_scan_steps.value),
               'structure': _structure_fingerprint(_current_xyz() or '')}
        floor = _scan_floor_for(leg)
        clipped = ''
        if floor is not None and leg['to'] < floor:
            clipped = (f' Asked for {leg["to"]:.3g}, which is inside the bond '
                       f'those two would make, so it stops at {floor:.3g}.')
            leg['to'] = floor
        legs.append(leg)
        state['scan_legs'] = legs
        _refresh_scan()
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
        if not _gfn.is_gfn_method(method):
            _set_mol_status('A scan needs xtb: choose a GFN method.')
            return
        if state.get('scan_run'):
            state['scan_stop'] = True
            _set_mol_status('Stopping the scan after this point.')
            return
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
        state['scan_crowded'] = None
        state['scan_frame_run'] = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = state['scan_frame_run']
        _ensure_manip_bootstrap()
        schedule_ui_update(_install_gfn_frame_watcher)
        submit_scan_run_btn.description = 'Stop'
        submit_scan_run_btn.icon = 'stop'
        label = _server_label(method)

        def _work():
            walked, path = xyz, []
            base = None
            bottom = None
            standing = None
            shown = []
            try:
                for n in range(1, steps + 1):
                    if state.get('scan_stop'):
                        break
                    share = n / float(steps)
                    held = [
                        {'kind': one['kind'], 'atoms': list(one['atoms']),
                         'mode': 'fix',
                         'value': one['from'] + share * (one['to'] - one['from'])}
                        for one in legs
                    ]
                    outcome = _gfn.optimize_with_gfn(
                        walked, method, charge=charge, uhf=uhf,
                        max_steps=_SCAN_CYCLES, timeout=None,
                        constraints=held, solvent=wet, solvation_model=model,
                        topology=_gfn_topology_dir(walked))
                    if not outcome.get('ok') or outcome.get('energy') is None:
                        schedule_ui_update(
                            _set_mol_status,
                            'The scan stopped at step '
                            f'{n}: {outcome.get("status") or "it did not run"}')
                        return
                    walked = outcome['xyz']
                    if base is None:
                        base = float(outcome['energy'])
                    spent = (float(outcome['energy']) - base) * _HARTREE_TO_KCAL
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
                    path.append((held[0]['value'], spent))
                    # The lowest point since the top, kept with its geometry.
                    #
                    # The climb out is what says the minimum is behind us, so
                    # by the time it is recognised the walk stands two steps
                    # past it and the structure in the box is squeezed that
                    # much.  Stopping is not enough; it has to come back to
                    # the bottom it walked through, which is the geometry
                    # anyone would want to carry on from.
                    if bottom is None or spent < bottom[0]:
                        bottom = (spent, walked, held[0]['value'])
                    if not submit_scan_whole.value and _scan_arrived(path):
                        state['scan_arrived'] = True
                        if bottom is not None:
                            walked = bottom[1]
                            path = path[:path.index(
                                (bottom[2], bottom[0])) + 1]
                        break
                    schedule_ui_update(
                        _set_mol_status,
                        f'{label} is walking the scan: step {n} of {steps}, '
                        f'{held[0]["kind"]} at {held[0]["value"]:.3g}, '
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
                    schedule_ui_update(
                        lambda text=json.dumps({
                            'run': state.get('scan_frame_run'),
                            'from': len(shown) - 1,
                            'follow': 1,
                            'frames': [shown[-1]],
                        }): setattr(submit_gfn_frame, 'value', text))
            finally:
                state['scan_run'] = False

                def _done(final=walked):
                    submit_scan_run_btn.description = 'Run scan'
                    submit_scan_run_btn.icon = 'play'
                    rows = [line for line in final.splitlines()[2:]
                            if line.strip()]
                    if rows:
                        _write_coords(xyz_document(
                            rows, 'Scanned to the next minimum'
                            if state.get('scan_arrived') else 'Scanned'))
                    _set_mol_status(*_scan_verdict(path, steps))

                schedule_ui_update(_done)

        threading.Thread(target=_work, daemon=True).start()

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
        ends = path[-1][1]
        T = float(submit_temperature.value or 298.15)
        ceiling = thermal_ceiling(T, _THERMAL_SECONDS)
        needs = thermal_temperature(top[1], _THERMAL_SECONDS)
        arrived = (f' It came back to the minimum it walked through, at '
                   f'{path[-1][0]:.3g}.'
                   if state.get('scan_arrived') else '')
        crowded = state.get('scan_crowded')
        if crowded:
            arrived = (f' It stopped: two atoms came inside {crowded:.2f} of a '
                       f'bond length, which is no path at any temperature. The '
                       f'target is past the far side of a bond.')
        first = (f'The scan walked {len(path)} of {steps} points. Highest '
                 f'{top[1]:+.1f} kcal/mol at {top[0]:.3g}, ending '
                 f'{ends:+.1f}.{arrived}')
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
        if top[1] <= ceiling:
            return (first,
                    f'{wants} You have {ceiling:.1f} kcal/mol at {T:g} K, so '
                    f'the whole path is open. {_thermal_wait(top[1], T)}')
        return (first,
                f'{wants} At {T:g} K only {ceiling:.1f} kcal/mol is '
                f'available, so the path is closed there. '
                f'{_thermal_wait(top[1], T)}')

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
        state['constraints'] = []
        state['bond_edits'] = {}
        state['hand_bonds'] = {}
        state['hyb_overrides'] = {}
        state['structure_undo'] = []
        # The history goes back to its first entry rather than being thrown
        # away: what Reset undid is one more thing that happened, and a user
        # who presses it by accident has to be able to come back.
        state['history'] = history[:1] if history else []
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
        _refresh_constraints()
        _set_mol_status('Back to the structure as it was loaded.')

    def on_submit_constraint_del(_button=None):
        key = submit_constraint_dd.value or ''
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
        pace = max(1, int(round(1000.0 / max(1, int(submit_play_speed.value)))))
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
        _set_mol_status(
            f'The optimisation is drawn at {int(submit_play_speed.value)} '
            'frame(s) a second. Slower lets the calculation run ahead of the '
            'picture -- take hold of an atom and the frame you are looking at '
            'is the one that is kept.')

    def on_submit_thermal(change):
        """Switching the budget on anchors it; switching it off forgets."""
        if change.get('name') != 'value':
            return
        active = bool(submit_thermal_btn.value)
        submit_thermal_btn.button_style = 'info' if active else ''
        for widget in (submit_temperature, submit_thermal_relax,
                       submit_thermal_anchor_btn):
            widget.layout.display = '' if active else 'none'
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
        _set_thermal_anchor()

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
        _set_mol_status(
            f'At {float(submit_temperature.value):g} K this structure has '
            f'{ceiling:.1f} kcal/mol to spend within {_timescale_label()}.'
            + (f' {spent}' if spent else ''))

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

    def on_submit_labels_toggle(change):
        """Numbers on or off, and nothing else.

        The molecule is emphatically not rendered again. Rebuilding it from
        the coordinates is how a structure loses what only the browser knows:
        the bonds as they were perceived, and the ones made or broken by hand.
        Switching the numbers off used to do exactly that, so a molecule that
        had been optimised came back with bonds missing. The numbers are a
        layer of sprites over the model -- they can be added and taken away
        without the model hearing about it.
        """
        if change.get('name') != 'value':
            return
        on = bool(submit_labels_btn.value)
        submit_label_size.layout.display = '' if on else 'none'
        submit_labels_btn.button_style = 'info' if on else ''
        _run_manip_js(
            show_atom_numbers_js(
                var=_submit_viewer_js(), on=on,
                scale=scale_for_px(submit_label_size.value))
        )

    def on_submit_label_size(change):
        """Resize them in the viewer that is already there.

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

        threading.Thread(target=_work, daemon=True).start()

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
        # Charge and spin: the server engines are told both, the browser's
        # field has no notion of either.
        submit_gfn_charge.layout.display = '' if server else 'none'
        submit_gfn_mult.layout.display = '' if server else 'none'
        # Scanning the multiplicity is xtb's; MOPAC is given the one on screen.
        submit_gfn_autospin.layout.display = '' if xtb else 'none'
        if not xtb and submit_gfn_autospin.value:
            submit_gfn_autospin.value = False
        # Strength is how many steps the browser's field takes per animation
        # frame, and that field does not run under a server method.
        submit_strength_slider.layout.display = 'none' if server else ''
        # The thermal budget needs an energy per follow step, and that is
        # xtb's: MOPAC's follow reports a heat of formation, which is not the
        # same quantity and cannot be differenced against an xtb anchor.
        submit_thermal_btn.layout.display = '' if xtb else 'none'
        if not xtb and submit_thermal_btn.value:
            submit_thermal_btn.value = False
        for widget in (submit_temperature, submit_thermal_relax,
                       submit_thermal_anchor_btn):
            widget.layout.display = ('' if (xtb and submit_thermal_btn.value)
                                     else 'none')
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
                    'multiplicity (M): xtb needs both, and a wrong spin on a '
                    'metal gives a confident wrong answer rather than an error. '
                    'Switch Dynamik Opt on to have the molecule follow an atom '
                    'you drag.')
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

    def on_submit_manip_sync(change):
        if change.get('name') != 'value':
            return
        new_xyz = submit_manip_sync.value
        if not new_xyz or not new_xyz.strip():
            return
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
                header = f'{old_lines[0]}\n{kept}\n'
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
        # Sent while the mouse is still down, so the molecule can follow the
        # atom rather than wait for it to be let go.  The comment line names
        # the atoms the hand is on, so the answer can keep them there.
        if note.startswith('DELFIN drag-follow'):
            holding = []
            for word in note.split():
                if word.startswith('held='):
                    holding = [int(n) for n in word[5:].split(',')
                               if n.strip().lstrip('-').isdigit()]
            _gfn_follow_step(new_xyz, holding)
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

        if drag_ended:
            # Set, Hold, a bond edit and a drag all arrive here.  Any of them
            # during an optimisation makes what xtb is doing about a structure
            # that is no longer on the screen, so it starts again from the one
            # that is -- and the geometry lands first, so it is the new one it
            # starts from.
            _interrupt_gfn()
            _arm_gfn_restart()
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

    submit_select_btn.observe(on_submit_select_toggle, names='value')
    submit_manip_btn.observe(on_submit_manip_toggle, names='value')
    submit_manip_clear_btn.on_click(on_submit_manip_clear)
    submit_manip_undo_btn.on_click(on_submit_manip_undo)
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
    submit_label_size.observe(on_submit_label_size, names='value')
    submit_strength_slider.observe(on_submit_strength_changed, names='value')
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
    # The page watches this button itself: waiting for the kernel to say the
    # switch went off costs a round trip, and the playback ran on for it.
    submit_optimize_btn.add_class('submit-optimize-switch')
    submit_optimize_btn.observe(on_submit_optimize, names='value')
    submit_optimize_all_btn.add_class('submit-optimize-switch')
    submit_optimize_all_btn.observe(on_submit_optimize_all, names='value')
    submit_manip_sync.observe(on_submit_manip_sync, names='value')

    # What the force field had to approximate belongs under the structure it
    # describes, not in the preview's status line where it competes with
    # conversion messages and scrolls away.
    submit_ff_notes = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', margin='4px 0 0 0'),
    )
    submit_ff_notes.add_class('submit-ff-notes')
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
                scale=scale_for_px(submit_label_size.value))
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
        for widget, value in zip(_structure_controls(), values or ()):
            try:
                if widget.value != value:
                    widget.value = value
            except Exception:
                # A list that no longer offers what it offered then.
                pass
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
        for widget in _controls_a_new_structure_resets():
            start = _CONTROL_START.get(id(widget))
            try:
                if widget.value != start:
                    widget.value = start
            except Exception:
                pass

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

    def handle_isomer_prev(button):
        if state['isomers']:
            _show_isomer_at_index(state['isomer_index'] - 1)

    def handle_isomer_next(button):
        if state['isomers']:
            _show_isomer_at_index(state['isomer_index'] + 1)

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

        threading.Thread(target=_worker, daemon=True).start()

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

        threading.Thread(target=_work, daemon=True).start()

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

        threading.Thread(target=_work, daemon=True).start()

    def on_submit_draw_get(_button=None):
        """Ask the editor for what has been drawn.

        The molfile, not the SMILES the editor could write itself: everything
        downstream here reads structures with RDKit, and a SMILES RDKit wrote
        is one RDKit will certainly read back.
        """
        _set_mol_status('Reading the drawing...', spinner=True)
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
