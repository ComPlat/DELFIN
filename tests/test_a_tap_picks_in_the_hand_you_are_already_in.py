"""Naming an atom in Manipulate mode, driven as a mouse drives it.

The editor had two hands and one of them could not point.  In Manipulate mode
every press on an atom was a drag, so selecting the atom a gesture is *about*
meant switching to Select, tapping, and switching back before the drag could
start: four steps for one gesture.

That pair is not decoration.  The interactive climb is checked against the two
atoms the gesture is about -- the one being dragged and the one it is being
driven towards -- and deriving that pair from the geometry was measured and
rejected: the best rule reaches 10 of 21 hand drags, and a wrong pair made the
climb throw away a saddle it had already found and spend 73 s more looking.  So
the pair has to be named by hand, and naming it has to cost one tap.

The rule is that a press which does not move is a tap and a press which moves
is a drag.  Everything here is about where that line sits and what it costs,
because a threshold that ate the first millimetre of every drag would be a far
worse bargain than the friction it removes.

There is no browser on this box, so the viewer's script -- which is a Python
string -- is handed to node against a page stubbed down to what the pointer
handling touches: an overlay that records its listeners, a canvas with a known
rectangle, and a viewer whose projection and whose ``getPixelToWorld`` are
given the *same* scale, so a pixel is worth an exactly known number of
Angstrom and what a gesture asked for can be read off the coordinates.

What was measured before the change, on a straight 20 px drag at 0.03 A per
pixel: the atom travelled 0.51 A of the 0.60 A the hand asked for.  The
threshold's pixels were spent and then thrown away -- ``lastX`` was advanced on
every move, including the moves that were still deciding -- so the atom lagged
the cursor by that 0.09 A for the rest of the gesture, and by 0.90 A on a view
ten times further out.  After: 0.60 A of 0.60 A, at every zoom, to within
5e-14.  The threshold is now free, which is what makes it safe to widen.
"""

from __future__ import annotations

import json
import pathlib
import shutil
import subprocess
import tempfile

import pytest

from delfin.dashboard import climb as _climb
from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js


_NODE = shutil.which("node")

#: A water dimer, laid flat in the plane of the screen: H1 is the hydrogen
#: that bridges, 0.96 A from its own oxygen and 1.90 A from the acceptor O3.
#: "Drag H1 at O3" is the kind of contact the interactive climb is checked
#: against, and it is the gesture this file drives.
#:
#: A model rather than a measured structure -- what is under test is the
#: wiring, not the chemistry -- but a closed-shell neutral one, so a climb can
#: really be started over it, and flat, with every atom further from its
#: neighbours than their drawn discs are wide, so a synthetic click on an
#: atom's centre can only ever be that atom.
_DIMER = """6
a water dimer, flat in the plane of the screen
O   0.00  0.00  0.00
H   0.96  0.00  0.00
H  -0.24  0.93  0.00
O   2.86  0.00  0.00
H   3.24  0.90  0.00
H   3.24 -0.90  0.00
"""

#: Which atom the gesture names, and which it drives at it.
_NAMED, _DRAGGED = 3, 1


def _atoms_of(xyz):
    """The molecule as the viewer holds it: serials from nought, in order.

    Derived from the document rather than written out beside it, so the two
    cannot drift apart -- an atom moved in one and not the other would move
    every synthetic click with it.
    """
    rows = [line.split() for line in xyz.splitlines()[2:] if line.strip()]
    return json.dumps([
        {"serial": i, "elem": row[0], "x": float(row[1]),
         "y": float(row[2]), "z": float(row[3])}
        for i, row in enumerate(rows)])


_ATOMS = _atoms_of(_DIMER)


_DRIVER = r"""
// --- a page, as much of one as the pointer handling touches ---------------
var W = 600, H = 600;

function el(tag) {
  return {
    tagName: (tag || 'DIV').toUpperCase(),
    style: {}, className: '', children: [], _h: {},
    classList: {add: function(){}, remove: function(){},
                contains: function(){ return false; }},
    addEventListener: function(t, fn){
      (this._h[t] = this._h[t] || []).push(fn);
    },
    removeEventListener: function(t, fn){
      var a = this._h[t] || []; var i = a.indexOf(fn); if (i >= 0) a.splice(i, 1);
    },
    appendChild: function(c){ this.children.push(c); return c; },
    removeChild: function(c){
      var i = this.children.indexOf(c); if (i >= 0) this.children.splice(i, 1);
    },
    querySelector: function(){ return null; },
    querySelectorAll: function(){ return []; },
    getBoundingClientRect: function(){
      return {left: 0, top: 0, right: W, bottom: H, width: W, height: H};
    },
    contains: function(){ return true; },
    dispatchEvent: function(){ return true; },
    fire: function(t, ev){
      (this._h[t] || []).slice().forEach(function(fn){ fn(ev); });
    }
  };
}

var canvas = el('canvas');
canvas.clientWidth = W; canvas.clientHeight = H;
var viewerEl = el('div');
viewerEl.querySelector = function(sel){ return sel === 'canvas' ? canvas : null; };

// The two channels the kernel reads a gesture through.
var pickInput = el('input'); pickInput.value = '';
var pickSync = el('div'); pickSync.querySelector = function(){ return pickInput; };
var xyzInput = el('textarea'); xyzInput.value = '';
var manipSync = el('div'); manipSync.querySelector = function(){ return xyzInput; };
var root = el('div');
root.querySelector = function(sel){
  if (sel === '.submit-pick-sync') return pickSync;
  if (sel === '.submit-manip-sync') return manipSync;
  return null;
};

globalThis.document = {
  querySelector: function(sel){ return sel === '.s1' ? root : null; },
  querySelectorAll: function(sel){ return sel === '.s1' ? [root] : []; },
  createElement: function(tag){ return el(tag); },
  addEventListener: function(){},
  body: {appendChild: function(){}, removeChild: function(){}},
  activeElement: null
};
var winH = {};
globalThis.window = {
  document: document,
  addEventListener: function(t, fn){ (winH[t] = winH[t] || []).push(fn); },
  removeEventListener: function(){},
  requestAnimationFrame: function(fn){ fn(); return 1; },
  cancelAnimationFrame: function(){},
  setTimeout: function(fn){ fn(); return 1; },
  clearTimeout: function(){},
  getComputedStyle: function(){ return {position: 'relative'}; },
  HTMLInputElement: {prototype: {}}, HTMLTextAreaElement: {prototype: {}},
  Event: function(){},
  performance: {now: function(){ return Date.now(); }}
};
globalThis.setTimeout = window.setTimeout;

__BOOTSTRAP__

// --- a molecule, and a camera whose scale is known exactly -----------------
var MODEL = __ATOMS__;
var PER_PX = 0.03;                       // Angstrom of world per screen pixel

function makeViewer(perPx) {
  var atoms = MODEL.map(function(a){
    return {serial: a.serial, elem: a.elem, x: a.x, y: a.y, z: a.z};
  });
  var model = {atoms: atoms, selectedAtoms: function(){ return atoms; }};
  var camZ = 150;
  // One Angstrom is exactly 1/perPx pixels in the projection *and* in
  // getPixelToWorld.  Letting those two disagree would mean picking and
  // dragging were measuring different screens.
  var s = 2 / (perPx * W);
  var fov = 2 * Math.atan(perPx * H / (2 * camZ)) * 180 / Math.PI;
  var shapes = [];
  return {
    atoms: atoms,
    getModel: function(){ return model; },
    selectedAtoms: function(){ return atoms; },
    addSphere: function(o){ var sh = {o: o}; shapes.push(sh); return sh; },
    addLine: function(o){ var sh = {o: o}; shapes.push(sh); return sh; },
    addLabel: function(o){ return {o: o}; },
    removeShape: function(sh){
      var i = shapes.indexOf(sh); if (i >= 0) shapes.splice(i, 1);
    },
    removeAllLabels: function(){}, render: function(){},
    setSlab: function(){}, getSlab: function(){ return {near: -50, far: 50}; },
    setClickable: function(){}, setStyle: function(){},
    rotationGroup: {position: {z: camZ},
                    matrix: {elements: [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]}},
    modelGroup: {matrixWorld: {elements: [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]},
                 updateMatrixWorld: function(){}},
    camera: {
      fov: fov,
      matrixWorldInverse: {elements: [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,-camZ,1]},
      projectionMatrix: {elements: [s,0,0,0, 0,s,0,0, 0,0,-1,0, 0,0,0,1]},
      updateMatrixWorld: function(){}
    }
  };
}

var SCOPE = 's1';
var LIVE = null;

function fresh(perPx, tweak) {
  var scale = (perPx === undefined) ? PER_PX : perPx;
  window._submitMolViewerByScope = {};
  window._submitManipStateByScope = {};
  winH.mousemove = []; winH.mouseup = []; winH.mousedown = [];
  viewerEl._h = {}; canvas._h = {};
  pickInput.value = ''; xyzInput.value = '';
  var viewer = makeViewer(scale);
  window._submitMolViewerByScope[SCOPE] = viewer;
  var api = window.__delfinSubmitManip;
  api.onViewerReady(SCOPE, viewerEl);
  var st = window._submitManipStateByScope[SCOPE];
  // The rigid hand: the atom goes where the cursor puts it, so what the
  // gesture asked for can be read straight off the coordinates.  Under the
  // soft hand the same delta becomes a force and the field answers it, which
  // is a different measurement.
  st.pullShare = 0;
  st.dragSensitivity = 1;
  if (tweak) tweak(st);
  if (st.ffActive) {
    window.__delfinFF = {
      load: function(){ return true; }, grab: function(){ return true; },
      step: function(){ return {moved: 0}; }, release: function(){},
      energy: function(){ return 0; }, stats: function(){ return null; }
    };
    window._delfinFFByScope = {}; window._delfinFFByScope[SCOPE] = {n: MODEL.length};
  } else {
    window.__delfinFF = null; window._delfinFFByScope = {};
  }
  api.setMode(SCOPE, 'manipulate');
  LIVE = {viewer: viewer, state: st, overlay: st.overlay, perPx: scale};
  return LIVE;
}

function screenOf(serial) {
  var a = LIVE.viewer.atoms[serial];
  return {x: 300 + a.x / LIVE.perPx, y: 300 - a.y / LIVE.perPx};
}
function ev(x, y, extra) {
  var e = {clientX: x, clientY: y, button: 0, shiftKey: false, ctrlKey: false,
           metaKey: false, preventDefault: function(){},
           stopPropagation: function(){}, stopImmediatePropagation: function(){},
           target: null};
  for (var k in (extra || {})) e[k] = extra[k];
  return e;
}
function down(x, y, extra){ LIVE.overlay.fire('mousedown', ev(x, y, extra)); }
function move(x, y, extra){
  (winH.mousemove || []).forEach(function(fn){ fn(ev(x, y, extra)); });
}
function up(x, y, extra){
  (winH.mouseup || []).forEach(function(fn){ fn(ev(x, y, extra)); });
}
function click(x, y, extra){ canvas.fire('click', ev(x, y, extra)); }

function picks(){ return LIVE.state.picks.map(function(p){ return p.serial; }); }
function at(serial){ var a = LIVE.viewer.atoms[serial]; return [a.x, a.y, a.z]; }
function reported(){ return pickInput.value; }

function tap(serial, extra) {
  var p = screenOf(serial);
  down(p.x, p.y, extra); up(p.x, p.y, extra);
}
function drag(serial, dx, dy, steps) {
  var p = screenOf(serial);
  steps = steps || Math.max(Math.abs(dx), Math.abs(dy));
  down(p.x, p.y);
  for (var i = 1; i <= steps; i++) {
    move(p.x + dx * i / steps, p.y + dy * i / steps);
  }
  up(p.x + dx, p.y + dy);
}

// The largest travel, in pixels, that is still a tap here: bisected on real
// gestures rather than read off a constant, because the constant is only half
// of the answer -- the zoom is the other half.
function tapSlopMeasured(perPx, tweak) {
  var lo = 0, hi = 600;
  for (var k = 0; k < 40; k++) {
    var mid = (lo + hi) / 2;
    fresh(perPx, tweak);
    var p = screenOf(1);
    down(p.x, p.y);
    move(p.x + mid, p.y);
    var isDrag = !!(LIVE.state.drag && LIVE.state.drag.movedEnough);
    up(p.x + mid, p.y);
    if (isDrag) hi = mid; else lo = mid;
  }
  return lo;
}

var NAMED = 3, DRAGGED = 1, HOME = MODEL[DRAGGED].x;
var out = {};
"""

_PROBE = r"""
// --- a clean tap picks; a press that moves does not ------------------------
fresh();
tap(NAMED);
out.tap = {picks: picks(), where: at(NAMED), reported: reported()};

fresh();
drag(DRAGGED, 20, 0);
out.drag = {picks: picks(), where: at(DRAGGED), asked: 20 * PER_PX};

// --- and the drag keeps every pixel it was given --------------------------
out.lost = {};
[0.3, 0.03, 0.003].forEach(function(p) {
  fresh(p);
  drag(DRAGGED, 60, 0);
  out.lost[String(p)] = Math.abs((60 * p) - (at(DRAGGED)[0] - HOME));
});

// --- a drag that creeps before it moves -----------------------------------
fresh();
var start = screenOf(DRAGGED);
down(start.x, start.y);
move(start.x + 1, start.y);
move(start.x + 2, start.y);
out.creepEarly = !!LIVE.state.drag.movedEnough;
for (var i = 3; i <= 30; i++) move(start.x + i, start.y);
up(start.x + 30, start.y);
out.creep = {picks: picks(), asked: 30 * PER_PX,
             lost: Math.abs(30 * PER_PX - (at(DRAGGED)[0] - HOME))};

// --- a tap after a drag, on the atom that was dragged ---------------------
fresh();
drag(DRAGGED, 20, 0);
out.afterDrag = picks();
var moved = screenOf(DRAGGED);
down(moved.x, moved.y); up(moved.x, moved.y);
out.tapAfterDrag = picks();

// --- a second tap, and a tap with a modifier ------------------------------
fresh();
tap(NAMED); var once = picks();
tap(NAMED); var twice = picks();
out.repeat = [once, twice];
fresh();
tap(0);
tap(NAMED, {shiftKey: true}); var added = picks();
tap(NAMED, {shiftKey: true}); var again = picks();
out.additive = [added, again];

// --- empty space, in both modes -------------------------------------------
fresh();
tap(NAMED);
down(30, 560); up(30, 560);
out.emptyManipulate = {picks: picks(), drag: !!LIVE.state.drag};
fresh();
window.__delfinSubmitManip.setMode(SCOPE, 'select');
var p2 = screenOf(NAMED);
click(p2.x, p2.y);
var picked = picks();
click(30, 560);
out.select = {tap: picked, empty: picks()};
click(p2.x, p2.y);
out.selectRepeat = picks();

// --- what the threshold actually is, at five zooms ------------------------
out.slop = {};
[0.3, 0.03, 0.01, 0.003, 0.0003].forEach(function(p) {
  var px = tapSlopMeasured(p);
  out.slop[String(p)] = {px: Math.round(px * 100) / 100,
                         angstrom: Math.round(px * p * 1e5) / 1e5};
});
// Sensitivity is part of what the hand asks for, so it is part of the answer.
out.slopSensitive = Math.round(
  tapSlopMeasured(0.003, function(st){ st.dragSensitivity = 5; }) * 100) / 100;

// --- a tap keeps nothing where it was put ---------------------------------
function withField(tweak) {
  return function(st){
    st.ffActive = true; st.settleOnRelease = false;
    if (tweak) tweak(st);
  };
}
fresh(undefined, withField());
tap(DRAGGED);
out.tapPinned = {pinned: (LIVE.state.pinned || []).slice(), picks: picks()};
fresh(undefined, withField());
drag(DRAGGED, 20, 0);
out.dragPinned = {pinned: (LIVE.state.pinned || []).slice(), picks: picks()};

// --- the pair, named and driven in one mode -------------------------------
fresh();
tap(NAMED);                               // the atom the gesture is about
out.pairNamed = {picks: picks(), reported: reported()};
var h = screenOf(DRAGGED);
down(h.x, h.y);                           // and the atom that is driven at it
for (i = 1; i <= 20; i++) move(h.x + i, h.y);
// What the page sends the kernel while the mouse is still down: the whole
// message, kept as it stands so the kernel can be handed exactly it.
window.__delfinSubmitManip.pushXyz(SCOPE, 'drag-follow');
out.pairFollowMessage = xyzInput.value;
out.pairFollow = xyzInput.value.split('\n')[1];
out.pairHolding = LIVE.state.drag.targets.slice();
out.pairPicksDuring = picks();
up(h.x + 20, h.y);
out.pairAfter = {picks: picks(), reported: reported(),
                 closed: at(DRAGGED)[0] - HOME};

console.log(JSON.stringify(out));
"""


def _drive():
    """Run every gesture in one node process and hand back what it measured."""
    script = (_DRIVER.replace("__BOOTSTRAP__", submit_manip_bootstrap_js())
                     .replace("__ATOMS__", _ATOMS) + _PROBE)
    folder = pathlib.Path(tempfile.mkdtemp())
    path = folder / "gestures.js"
    path.write_text(script)
    done = subprocess.run([_NODE, str(path)], capture_output=True, text=True,
                          timeout=300)
    assert done.returncode == 0, done.stderr[-2000:]
    return json.loads(done.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def gestures():
    if _NODE is None:
        pytest.skip("node is not installed; the viewer cannot be driven")
    return _drive()


def test_a_tap_on_an_atom_picks_it_without_leaving_the_hand(gestures):
    """The whole point: one tap, in the mode the hand is already in.

    Before, the pointer handler turned every press on an atom into a drag and
    nothing else, so this tap selected nothing at all -- ``picks`` came back
    empty and the pick channel the kernel reads stayed blank.  Now the atom is
    selected and reported, and it has not moved: the tap is a name, not an
    edit.
    """
    assert gestures["tap"]["picks"] == [_NAMED]
    assert gestures["tap"]["where"] == [2.86, 0.0, 0.0], "a tap moved the atom"
    # As a model index, which is how the pick channel and the drag's own
    # ``held=`` both number atoms -- one numbering for both halves of a pair.
    assert gestures["tap"]["reported"] == str(_NAMED)


def test_a_press_that_moves_is_still_a_drag_and_picks_nothing(gestures):
    """The gesture people use constantly, unchanged.

    Twenty pixels at 0.03 A per pixel is 0.60 A of asked-for movement; the
    atom is there, and the selection has not been touched on the way.
    """
    assert gestures["drag"]["picks"] == []
    assert gestures["drag"]["where"][0] == pytest.approx(0.96 + 0.60,
                                                        abs=1e-9)


def test_the_drag_keeps_every_pixel_it_was_asked_for(gestures):
    """Widening the threshold had to be free, and this is what made it free.

    ``lastX`` used to advance on every move, including the moves that were
    still deciding whether this was a drag at all -- so the frame that crossed
    the threshold moved the atom by that one frame's pixels and the travel
    before it was gone.  What was lost is the threshold's own three pixels,
    whatever the drag: measured on a straight 60 px drag before the change,
    1.71 A of the 1.80 A asked at 0.03 A per pixel, 17.10 of 18.00 at 0.3 and
    0.171 of 0.180 at 0.003 -- and the atom lagged the cursor by that gap for
    the whole of the gesture, not only its start.

    Held at the start until the press is a drag, the deciding frame carries
    the whole travel.  Nothing is lost at any zoom, which is why the threshold
    may now be as wide as telling a tap from a drag actually needs.
    """
    for zoom, lost in gestures["lost"].items():
        assert lost < 1e-12, (zoom, lost)


def test_a_drag_that_creeps_before_it_moves_keeps_its_first_frames(gestures):
    """Two sub-threshold pixels and then a real movement is one drag.

    Not a tap that turned into a drag and not a drag that lost its opening:
    the two creeping pixels are still there at the end, and nothing was
    selected on the way.
    """
    assert gestures["creepEarly"] is False, "two pixels started a drag"
    assert gestures["creep"]["picks"] == []
    assert gestures["creep"]["lost"] < 1e-12, gestures["creep"]


def test_a_tap_after_a_drag_picks_the_atom_that_was_dragged(gestures):
    """The two gestures do not contaminate each other.

    A drag leaves the selection alone, and a tap on the same atom afterwards
    -- at the place the drag left it -- names it.  That order is the one the
    climb wants when a hand changes its mind about which contact it meant.
    """
    assert gestures["afterDrag"] == []
    assert gestures["tapAfterDrag"] == [_DRAGGED]


def test_a_second_tap_takes_the_pick_back_and_a_modifier_does_not(gestures):
    """One rule for both modes rather than one rule each.

    Manipulate goes through the same ``togglePick`` a tap in Select mode goes
    through, so a plain tap toggles, a tap with Shift, Ctrl or Cmd adds and
    never subtracts, and neither mode has to be learnt twice.
    """
    assert gestures["repeat"] == [[_NAMED], []]
    assert gestures["additive"] == [[0, _NAMED], [0, _NAMED]]


def test_a_tap_on_empty_space_leaves_the_selection_where_it_was(gestures):
    """In both modes, and for the same reason: nothing was named.

    In Manipulate a press on empty space is handed to the viewer so the camera
    still turns, and it never becomes a drag of this editor's -- so there is
    nothing at the release to select and nothing to clear.  Select mode has
    always behaved that way, and clearing on a stray click would be the more
    annoying of the two possible rules: a selection is expensive to build and
    a click beside the molecule is cheap to make.
    """
    assert gestures["emptyManipulate"] == {"picks": [_NAMED], "drag": False}
    assert gestures["select"] == {"tap": [_NAMED], "empty": [_NAMED]}


def test_a_tap_in_select_mode_still_picks_the_way_it_did(gestures):
    """The mode that could already do this was not touched.

    It resolves through 3Dmol's own click path rather than through the
    overlay, and it still toggles: one tap selects, the next takes it back.
    """
    assert gestures["select"]["tap"] == [_NAMED]
    assert gestures["selectRepeat"] == []


def test_the_threshold_is_angstrom_held_between_a_floor_and_a_ceiling(gestures):
    """The unit a threshold in this editor has to be stated in.

    The drag turns pixels into world movement through the camera, so a
    threshold in pixels means a different thing at every zoom: the three
    pixels this editor used were 0.90 A on a wide view and 0.0009 A on a close
    one.  On the close one that is nothing at all -- a tap that wobbles five
    pixels asks the structure to move by a thousandth of a bond -- and it
    counted as a drag and named nothing, which is exactly the case tapping was
    added for.  0.05 A is the size below which the picture does not move: it
    is under the width of the stick an atom is drawn on.

    The floor and the ceiling are where the Angstrom stops being honest.
    Zoomed far out a single pixel is worth more than the threshold and every
    tremor would be a drag, so three pixels -- what this editor has always
    used -- is the least it may ask.  Zoomed far in the threshold would be
    hundreds of pixels and the atom would stand still while the hand swept the
    screen, so twelve pixels is the most: still a tap by the touch slop the
    platforms use, 8 dp on Android and about 10 px on iOS.

    Measured across five zooms, in pixels and in the Angstrom they are worth:
    0.3 A/px -> 3 px (0.9 A, floor), 0.03 -> 3 px (0.09 A, floor), 0.01 -> 5 px
    (0.05 A, the rule itself), 0.003 -> 12 px (0.036 A, ceiling), 0.0003 ->
    12 px (0.0036 A, ceiling).  No time bound was needed with it: a press that
    is held still does nothing in this mode whether it is held for a moment or
    a minute, so there is no second gesture for a clock to tell apart -- and a
    clock would only invent a class of presses that look like taps and are
    not.
    """
    slop = gestures["slop"]
    assert slop["0.3"]["px"] == pytest.approx(3, abs=0.01)
    assert slop["0.03"]["px"] == pytest.approx(3, abs=0.01)
    assert slop["0.01"]["angstrom"] == pytest.approx(0.05, abs=1e-4)
    assert slop["0.003"]["px"] == pytest.approx(12, abs=0.01)
    assert slop["0.0003"]["px"] == pytest.approx(12, abs=0.01)
    # Never below the floor and never above the ceiling, whatever the zoom.
    for zoom, found in slop.items():
        assert 3 - 0.01 <= found["px"] <= 12 + 0.01, (zoom, found)
    # Sensitivity is part of what the hand asks for, so it is part of the
    # answer.  Told to move the structure five times as far per pixel, the
    # same 0.05 A is reached in 3.33 px rather than the 16.7 px it would take
    # at one -- which is above the ceiling, and was cut to 12 there.
    assert gestures["slopSensitive"] == pytest.approx(0.05 / (0.003 * 5),
                                                     abs=0.02)


def test_a_tap_keeps_nothing_where_it_was_put(gestures):
    """With settling on release switched off, a release pins what it held.

    That is what the switch is for -- an atom placed by hand stays where it
    was put while the rest of the molecule settles around it -- but a press
    that never moved has put nothing anywhere.  Measured before the change: a
    click on an atom in Manipulate mode with the switch off pinned it and
    selected nothing, so the gesture this file is about would have quietly
    frozen every atom it named.  The release is handed an empty list when the
    press was a tap; a real drag still pins.
    """
    assert gestures["tapPinned"] == {"pinned": [], "picks": [1]}
    assert gestures["dragPinned"] == {"pinned": [1], "picks": []}


def test_one_gesture_in_one_mode_names_the_pair_the_climb_is_checked_against(
        gestures):
    """The point of the exercise, end to end on the page.

    Tap the oxygen the proton is being driven towards, then drag the proton at
    it -- both in Manipulate mode, with nothing switched in between.  The two
    halves of the pair leave the page on their own channels and in the same
    numbering: the pick channel carries the atom that was named, and the
    geometry the drag pushes while the mouse is still down carries ``held=``
    with the atom that is being moved.  The drag does not disturb the name it
    is aimed at, which is what makes the two readable together.

    Before, the tap did neither: the pick channel stayed empty and the pair
    had one half.
    """
    assert gestures["pairNamed"] == {"picks": [_NAMED],
                                     "reported": str(_NAMED)}
    assert gestures["pairHolding"] == [_DRAGGED]
    assert gestures["pairPicksDuring"] == [_NAMED], "the drag lost the name"
    assert "DELFIN drag-follow" in gestures["pairFollow"]
    assert f"held={_DRAGGED}" in gestures["pairFollow"], gestures["pairFollow"]
    assert gestures["pairAfter"]["picks"] == [_NAMED]
    assert gestures["pairAfter"]["reported"] == str(_NAMED)
    # And the drag did what it looked like: 20 px at 0.03 A/px took the
    # bridging hydrogen 0.60 A of the 1.90 A towards the oxygen.
    assert gestures["pairAfter"]["closed"] == pytest.approx(0.60, abs=1e-9)


def test_the_pair_the_gesture_names_is_the_contact_the_climb_asks_about():
    """And the pair is a legal answer to the question the ladder asks.

    ``held`` is not a label the climb carries around: it is turned into the
    one direction the verdict dots the imaginary mode against, the mass
    weighted stretch of that contact.  The pair this gesture produces --
    ``(dragged, named)`` = (1, 3) -- is exactly that, on the very structure the
    gesture was made on: only the bridging hydrogen and the oxygen it was
    driven at move, they move against each other, and the hydrogen carries
    most of it because it is the atom that actually moves when that distance
    changes.

    Nothing here needs xtb; it is the geometry of one direction.
    """
    import numpy as np

    walk = _climb.Climb(_DIMER, "gfn2")
    try:
        unit = walk.held_stretch((_DRAGGED, _NAMED))
        assert unit is not None
        assert float(np.linalg.norm(unit)) == pytest.approx(1.0)
        rows = unit.reshape(-1, 3)
        moving = [i for i in range(len(rows))
                  if float(np.linalg.norm(rows[i])) > 1e-12]
        assert moving == [_DRAGGED, _NAMED], moving
        assert float(rows[_DRAGGED] @ rows[_NAMED]) < 0
        assert (np.linalg.norm(rows[_DRAGGED])
                > 3.0 * np.linalg.norm(rows[_NAMED]))
    finally:
        walk.close()


def _an_editor(room):
    """A real structure editor over the dimer, driven from the kernel side.

    The host does the one thing these tests depend on: a write to the
    coordinate box refreshes what ``_current_xyz`` answers with, so an
    optimiser started here is started from the structure on screen rather than
    from the one the session opened on.
    """
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office"):
        (room / name).mkdir(parents=True, exist_ok=True)
    ctx = DashboardContext(calc_dir=room / "calc", archive_dir=room / "archive",
                           office_dir=room / "office")
    ctx.run_js = lambda _script: None
    state: dict = {}
    box = widgets.Textarea(value=_DIMER)

    def update_view(*_a, **_k):
        rows = [one for one in (box.value or "").split("\n") if one.strip()]
        body = rows[2:] if rows and rows[0].strip().isdigit() else rows
        state["current_xyz_for_copy"] = {
            "content": f"{len(body)}\nEdited in DELFIN viewer\n" + "\n".join(body)}

    part = structure_editor.build(
        ctx, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda call, *a, **k: call(*a, **k),
        update_view=update_view, get_smiles_charge=lambda *a, **k: None)
    box.observe(lambda _change: update_view(), names="value")
    update_view()
    return part, state, box


def test_the_atom_the_tap_named_reaches_the_kernel(tmp_path):
    """The last step of the first half, in a real editor.

    The page writes the selection into the hidden pick field as model indices;
    this is the kernel reading it.  Driven on a real ``structure_editor``
    rather than asserted about its source, because what is being checked is
    that the value written by a tap lands in the state a pair is built from --
    and that is an observer, which source cannot show working.
    """
    pytest.importorskip("ipywidgets")
    part, state, _box = _an_editor(tmp_path)

    part.submit_pick_sync.value = str(_NAMED)
    assert state.get("picked") == [_NAMED]
    # And a second tap, which clears the field, clears it there too -- so the
    # kernel never holds a name the picture is no longer showing.
    part.submit_pick_sync.value = ""
    assert state.get("picked") == []


@pytest.mark.skipif(not _climb.have_fast_gradients()
                    and _climb._gfn.find_xtb() is None,
                    reason="no xtb to take gradients from")
def test_the_gesture_hands_the_climb_the_pair_it_is_checked_against(
        gestures, tmp_path):
    """The point of the exercise, from the mouse to ``climb_held``.

    The two messages a real editor is given here are not written by hand: they
    are the bytes the node run above got out of the viewer for one gesture in
    Manipulate mode -- the pick channel after the tap, and the whole
    ``DELFIN drag-follow`` document the drag pushes while the mouse is still
    down, ``held=`` and all.  With Climb to TS on, the editor reads one atom
    off each and writes the pair down.

    That is the whole reason for the tap.  The climb's answer is judged on
    whether what it reached is the reaction the hand pointed at, which is a
    question about two atoms; with no pair only the count of imaginary modes
    can be asked, and then every later rung is judged by the test the first
    one already passed -- so the ladder stops at one rung.  Before this
    change, both halves could not be produced without leaving the mode: the
    tap in Manipulate named nothing, so ``picked`` stayed empty and
    ``_name_the_pair_the_hand_means`` had one atom and wrote nothing down.
    """
    pytest.importorskip("ipywidgets")

    import time as clock

    part, state, _box = _an_editor(tmp_path)
    part.submit_ff_dd.value = "gfn2"

    # The tap, exactly as the page reported it.
    part.submit_pick_sync.value = gestures["pairNamed"]["reported"]
    assert state.get("picked") == [_NAMED]

    # Climb to TS on.  Dynamik Opt is left alone here: what the pair is read
    # from is the drag-follow message, and arming the live relaxation as well
    # would only add an xtb round between the hand and the assertion.
    part.submit_climb_btn.value = True
    assert part.submit_climb_btn.value, "the method offered no climb"

    # And the drag, as the page pushed it mid-gesture.
    part.submit_manip_sync.value = gestures["pairFollowMessage"]

    # The atom being dragged first, the atom it is being driven at second --
    # which is the order the verdict's ``held_stretch`` reads it in; that the
    # pair is a legal direction on this structure is measured just above.
    assert state.get("climb_held") == (_DRAGGED, _NAMED), state.get("climb_held")

    part.submit_climb_btn.value = False
    began = clock.time()
    while state.get("climb_run") is not None and clock.time() - began < 120:
        clock.sleep(0.05)
    assert state.get("climb_run") is None, "a climb was left walking"


def test_a_tap_does_not_grab_and_so_cuts_nothing_that_is_running():
    """The other half of "a tap is not a drag", on the kernel's side.

    The page's frame loop watches one predicate -- ``grabbed()`` -- and sends
    ``gfngrab`` when it turns true and ``gfnfree`` when it turns false.  Those
    two messages are not observations: the grab interrupts whatever is running
    and the release arms it again, so a press and a release a tenth of a second
    apart cut a relaxation or a climb and started it over.

    Any click has always done that.  What changed is that a tap is now a
    gesture the user makes on purpose -- it is how the atom a climb is aimed at
    gets named -- so it goes from an accident to something he does while a walk
    is running, and the walk stutters for no reason he can connect to anything.

    A press that has not passed the slop has moved nothing: ``translate`` and
    ``rotate`` both do their work under ``movedEnough``, which is asserted by
    the node runs above.  So the predicate reads it too, and a tap sends
    neither message.  Draw is not held to it, because there a press that does
    not move still places an atom.
    """
    from editor_source import EDITOR_SOURCE

    watcher = EDITOR_SOURCE.split("function grabbed()", 1)[1].split(
        "function heldSerials()", 1)[0]
    assert 'drag.kind==="translate"||drag.kind==="rotate")' in watcher
    assert "return !!drag.movedEnough;" in watcher
    # And draw still counts whether or not it moved.
    assert 'return drag.kind==="draw";' in watcher
    # The one predicate is what both messages hang off, so reading it is
    # reading both.
    loop = EDITOR_SOURCE.split("var held=grabbed();", 1)[1].split(
        "function heldSerials", 1)[0]
    assert 'send(held?"gfngrab":"gfnfree"' in loop


def test_what_a_tap_would_have_reset_had_it_been_a_grab():
    """Why this is worth more than a stutter: the budget hangs off the same two
    messages.

    ``gfngrab`` begins a follow, and beginning one clears what the last drag
    went through and takes the geometry the budget may fall back to from the
    structure as it stands; ``gfnfree`` ends it, and ending one clears the
    running maximum again.  For a tap all of that is a no-op in effect -- the
    structure did not change between the press and the release, so the
    fallback geometry is the same one and the maximum was already cleared by
    the previous release -- but it is a no-op by luck rather than by design,
    and the interrupt beside it is not one at all.

    Read here so that the claim is on the record rather than assumed: these are
    the things that hang off the messages a tap no longer sends.
    """
    from editor_source import EDITOR_SOURCE

    grab = EDITOR_SOURCE.split("if verb == 'gfngrab'", 1)[1].split(
        "if verb == 'gfnfree'", 1)[0]
    assert "_interrupt_gfn()" in grab
    assert "_begin_gfn_follow()" in grab
    assert "_arm_thermal_leash()" in grab

    free = EDITOR_SOURCE.split("if verb == 'gfnfree'", 1)[1].split(
        "if verb == 'undo'", 1)[0]
    assert "_end_gfn_follow()" in free
    assert "_after_release()" in free

    begins = EDITOR_SOURCE.split("def _begin_gfn_follow", 1)[1].split(
        "def _clear_thermal_wall", 1)[0]
    assert "_clear_thermal_wall()" in begins
    ends = EDITOR_SOURCE.split("def _end_gfn_follow", 1)[1].split(
        "\n    #:", 1)[0]
    assert "_clear_thermal_wall()" in ends
