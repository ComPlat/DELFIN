"""A torque on the right button, driven the way a mouse drives one.

Asked for while the field was running: "kann man im Dynamik opt vielleicht
auch mit rechter maus taste eine drehkraft ansetzen an einem atom ... um
beispielsweise einfach bindungen zu rotieren", turning left or right by
dragging up or down.

xtb has no torque to give.  What it has is a torsion, and a turn *is* a
torsion -- which is what the left hand has driven all along: ``turn_for``
picks the soft coordinate a drag really moves, because a bond is the one
coordinate a hand must not have.  So this is the same force on the same kind
of coordinate.  What is new is that the user names the bond instead of the
editor inferring it from the direction they happened to drag.

The button was free.  While the field runs, pivot picking is off and the
right button fell through to the scene pan; empty space still pans, so
nothing has been taken away.

There is no browser on this box, so the viewer's script -- a Python string --
is handed to node against a page stubbed down to what the pointer handling
touches.  The driver is the one the tap tests use, so the two cannot measure
different screens.
"""

from __future__ import annotations

import json
import pathlib
import shutil
import subprocess
import tempfile

import pytest

from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js
from test_a_tap_picks_in_the_hand_you_are_already_in import _DRIVER, _atoms_of


_NODE = shutil.which("node")

#: Ethane, opened out so that no atom hides another in projection: a synthetic
#: click on an atom's centre can then only ever be that atom.  The central
#: bond is the one a hand takes hold of, and both halves carry three
#: hydrogens, so the torsion is a real one whichever way it is read.
_ETHANE = """8
ethane, spread for the screen
C   0.00  0.00  0.00
C   1.54  0.00  0.00
H  -0.55  1.00  0.00
H  -0.55 -1.00  0.00
H  -0.55  0.00  1.00
H   2.09  1.00  0.00
H   2.09 -1.00  0.00
H   2.09  0.00 -1.00
"""

#: Cyclopropane: every carbon-carbon bond is a ring bond, and a torsion about
#: a ring bond turns nothing -- both ends stay joined the other way round.
_RING = """9
cyclopropane, flat in the plane of the screen
C   0.00  0.87  0.00
C  -0.75 -0.43  0.00
C   0.75 -0.43  0.00
H   0.00  1.95  0.90
H   0.00  1.95 -0.90
H  -1.69 -0.98  0.90
H  -1.69 -0.98 -0.90
H   1.69 -0.98  0.90
H   1.69 -0.98 -0.90
"""


def _drive(xyz, probe):
    """Run one page's worth of gestures in node and hand back what it saw."""
    script = (_DRIVER.replace("__BOOTSTRAP__", submit_manip_bootstrap_js())
                     .replace("__ATOMS__", _atoms_of(xyz)) + probe)
    folder = pathlib.Path(tempfile.mkdtemp())
    path = folder / "turning.js"
    path.write_text(script)
    done = subprocess.run([_NODE, str(path)], capture_output=True, text=True,
                          timeout=300)
    assert done.returncode == 0, done.stderr[-2000:]
    return json.loads(done.stdout.strip().splitlines()[-1])


#: Everything the ethane page is asked, in one process.
#:
#: ``fresh`` is given the field switch, because that is the whole condition:
#: with it off the right button picks a pivot as it always has, and with it on
#: the same press is a torque.
_ETHANE_PROBE = r"""
function dihedral(i, j, k, l) {
  var A = LIVE.viewer.atoms;
  function sub(a,b){return {x:a.x-b.x,y:a.y-b.y,z:a.z-b.z};}
  function cross(a,b){return {x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,
                              z:a.x*b.y-a.y*b.x};}
  function dot(a,b){return a.x*b.x+a.y*b.y+a.z*b.z;}
  function len(a){return Math.sqrt(dot(a,a));}
  var b1 = sub(A[j],A[i]), b2 = sub(A[k],A[j]), b3 = sub(A[l],A[k]);
  var n1 = cross(b1,b2), n2 = cross(b2,b3), m = cross(n1, {
    x: b2.x/len(b2), y: b2.y/len(b2), z: b2.z/len(b2)});
  return Math.atan2(dot(m,n2), dot(n1,n2)) * 180 / Math.PI;
}
function span(i, j) {
  var A = LIVE.viewer.atoms;
  var dx = A[i].x-A[j].x, dy = A[i].y-A[j].y, dz = A[i].z-A[j].z;
  return Math.sqrt(dx*dx + dy*dy + dz*dz);
}
function rightDrag(serial, dy) {
  var p = screenOf(serial);
  down(p.x, p.y, {button: 2});
  for (var i = 1; i <= Math.max(1, Math.abs(dy)); i++) {
    move(p.x, p.y + dy * i / Math.max(1, Math.abs(dy)), {button: 2});
  }
  return LIVE.state.drag;
}
// Dynamik Opt, through the door the editor uses.  Setting autoOpt here
// instead would test a mode the editor switches *off* under every server
// method -- which is how the torque came to be unreachable in the first
// place.
function field(st){ window.__delfinSubmitManip.setLiveHand(SCOPE, true); }

var out = {};

// --- the press names a torsion, and only while the field runs -------------
fresh(undefined, field);
var p0 = screenOf(0);
down(p0.x, p0.y, {button: 2});
out.started = LIVE.state.drag ? {kind: LIVE.state.drag.kind,
                                 idx: LIVE.state.drag.idx.slice(),
                                 held: LIVE.state.drag.targets.slice()} : null;
up(p0.x, p0.y, {button: 2});

fresh();                                  // the field off: pivot, as before
down(p0.x, p0.y, {button: 2});
out.withoutField = {kind: LIVE.state.drag && LIVE.state.drag.kind,
                    pivot: LIVE.state.pivot && LIVE.state.pivot.serial};
up(p0.x, p0.y, {button: 2});

// --- up and down turn opposite ways ---------------------------------------
fresh(undefined, field);
var idx = null;
var before = null, upTurn = null;
var d = rightDrag(0, -100);               // up the screen
idx = d.idx.slice();
upTurn = dihedral(idx[0], idx[1], idx[2], idx[3]);
up(0, 0, {button: 2});

fresh(undefined, field);
before = (function(){ var q = turnStart(); return q; })();
function turnStart() {
  var p = screenOf(0);
  down(p.x, p.y, {button: 2});
  var q = LIVE.state.drag.idx.slice();
  var was = dihedral(q[0], q[1], q[2], q[3]);
  up(p.x, p.y, {button: 2});
  return {idx: q, was: was};
}
fresh(undefined, field);
rightDrag(0, 100);                        // down the screen
var downTurn = dihedral(idx[0], idx[1], idx[2], idx[3]);
up(0, 0, {button: 2});
out.turn = {was: before.was, up: upTurn, down: downTurn};

// --- and nothing is stretched by it ---------------------------------------
fresh(undefined, field);
var spans = [[0,1],[0,2],[0,3],[0,4],[1,5],[1,6],[1,7]];
var wasSpans = spans.map(function(s){ return span(s[0], s[1]); });
rightDrag(0, -80);
out.stretched = spans.map(function(s, n){
  return Math.abs(span(s[0], s[1]) - wasSpans[n]); });
up(0, 0, {button: 2});

// --- the torsion rides with the coordinates -------------------------------
fresh(undefined, field);
rightDrag(0, -60);
// The push the player makes while a hand is down, made here by hand: this
// page carries the manipulator and not the player, and what is under test is
// the comment the manipulator writes.
window.__delfinSubmitManip.pushXyz(SCOPE, 'drag-follow');
out.note = xyzInput.value.split('\n')[1];
up(0, 0, {button: 2});

// --- a press that does not move is not a turn -----------------------------
fresh(undefined, field);
var pt = screenOf(0);
down(pt.x, pt.y, {button: 2});
move(pt.x, pt.y + 1, {button: 2});
out.tapIsNoTurn = {moved: !!LIVE.state.drag.movedEnough,
                   dihedral: dihedral(LIVE.state.drag.idx[0],
                                      LIVE.state.drag.idx[1],
                                      LIVE.state.drag.idx[2],
                                      LIVE.state.drag.idx[3])};
up(pt.x, pt.y + 1, {button: 2});

// --- empty space still belongs to the viewer ------------------------------
fresh(undefined, field);
down(5, 5, {button: 2});
out.emptySpace = LIVE.state.drag ? LIVE.state.drag.kind : null;

console.log(JSON.stringify(out));
"""


_RING_PROBE = r"""
function field(st){ window.__delfinSubmitManip.setLiveHand(SCOPE, true); }
var out = {};
fresh(undefined, field);
var p = screenOf(0);
down(p.x, p.y, {button: 2});
out.drag = LIVE.state.drag ? LIVE.state.drag.kind : null;
console.log(JSON.stringify(out));
"""


@pytest.fixture(scope="module")
def turning():
    if _NODE is None:
        pytest.skip("node is not installed; the viewer cannot be driven")
    return _drive(_ETHANE, _ETHANE_PROBE)


@pytest.fixture(scope="module")
def ring():
    if _NODE is None:
        pytest.skip("node is not installed; the viewer cannot be driven")
    return _drive(_RING, _RING_PROBE)


def test_a_right_press_under_the_field_names_a_torsion(turning):
    """Four atoms, and the bond the hand took hold of in the middle of them.

    The axis is the bond the grabbed atom hangs on -- take hold of an atom
    and turn what is beyond it -- which is a different question from the one
    the left hand asks, where the axis is a bond one further in so that the
    grabbed atom itself swings through an arc.
    """
    started = turning['started']
    assert started is not None, 'the right button did nothing'
    assert started['kind'] == 'turn'
    assert len(started['idx']) == 4
    assert started['idx'][1:3] == [0, 1], (
        'the axis is the bond the grabbed atom hangs on', started['idx'])
    assert started['held'], 'the atoms that turn are named for the kernel'


def test_with_the_field_off_the_right_button_still_picks_a_pivot(turning):
    """Nothing is taken away. The torque lives in the space the right button
    was not using: while the field runs, pivot picking is off."""
    assert turning['withoutField']['kind'] != 'turn'
    assert turning['withoutField']['pivot'] == 0


def test_dragging_up_and_down_turns_opposite_ways(turning):
    """A torsion has one number and one sign, so up is one way round and down
    is the other, and sideways is left to mean nothing."""
    was, upward, downward = (turning['turn']['was'], turning['turn']['up'],
                             turning['turn']['down'])

    def moved(a, b):
        return ((b - a + 540) % 360) - 180

    assert abs(moved(was, upward)) > 20, (was, upward)
    assert abs(moved(was, downward)) > 20, (was, downward)
    assert moved(was, upward) * moved(was, downward) < 0, (
        'up and down turned the same way', was, upward, downward)


def test_a_turn_stretches_nothing(turning):
    """The whole promise of driving a torsion rather than a position: bond
    lengths hold and the torsion reacts. A hand that moved atoms in space
    would be stretching the bonds it turns about."""
    assert max(turning['stretched']) < 1e-6, turning['stretched']


def test_the_torsion_rides_with_the_coordinates(turning):
    """Named on the way to the kernel rather than worked out again there.
    Which coordinate a hand is driving is a decision, and one the user made
    by choosing a bond must not be re-taken from the geometry a frame
    later."""
    note = turning['note']
    assert 'drag-follow' in note, note
    assert ' turn=' in note, note
    named = [int(n) for n in note.split('turn=')[1].split()[0].split(',')]
    assert len(named) == 4 and named[1:3] == [0, 1], note


def test_a_press_that_does_not_move_turns_nothing(turning):
    """The same rule the other two gestures answer to: a press that has not
    passed the slop has moved nothing. It matters here because a tap on an
    atom is a gesture people make on purpose."""
    assert turning['tapIsNoTurn']['moved'] is False


def test_empty_space_still_pans(turning):
    assert turning['emptySpace'] is None


def test_a_ring_bond_is_refused_before_the_gesture_starts(ring):
    """A torsion about a ring bond turns nothing -- both ends stay joined the
    other way round -- and driving one fights the ring. Refused at the press,
    so the refusal costs no gesture."""
    assert ring['drag'] is None


# ---------------------------------------------------------------------------
# And the other end of it
# ---------------------------------------------------------------------------


_ETHANE_XYZ = '\n'.join(
    ['8', 'DELFIN drag-follow held=0,2,3,4 turn=2,0,1,5']
    + [line for line in _ETHANE.splitlines()[2:] if line.strip()])


def _an_editor(text):
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda _script: None
    state = {}
    part = structure_editor.build(
        ctx, state=state, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    return part, state


def test_the_torsion_the_hand_named_is_the_one_that_is_held():
    """It goes into the same slot a left hand's own decision lands in, so
    :func:`gfn_optimize.contacts_holding` holds this torsion and nothing else
    and :func:`gfn_optimize.as_pushes` makes a force of it -- under the
    slider and the thermal ceiling, like every other hand."""
    part, state = _an_editor(_ETHANE)
    part.submit_manip_sync.value = _ETHANE_XYZ
    assert state.get('thermal_turn') == [2, 0, 1, 5]


def test_a_named_turn_is_not_asked_about_again():
    """The opening answer of a drag asks what each coordinate *can* do,
    because on the first answer almost nothing has moved. There is nothing
    left to ask once the hand has said which coordinate it is on."""
    part, state = _an_editor(_ETHANE)
    state['gfn_follow_opening'] = True
    part.submit_manip_sync.value = _ETHANE_XYZ
    assert not state.get('gfn_follow_opening')


def test_a_drag_that_names_no_turn_leaves_the_decision_alone():
    """The left hand still decides for itself, and a message without a turn
    in it must not clear one that a turn put there."""
    part, state = _an_editor(_ETHANE)
    state['thermal_turn'] = [2, 0, 1, 5]
    part.submit_manip_sync.value = _ETHANE_XYZ.replace(
        ' turn=2,0,1,5', '')
    assert state.get('thermal_turn') == [2, 0, 1, 5]


def test_the_four_atoms_are_held_as_a_dihedral_and_nothing_else():
    """The end of the chain, in the module that speaks to xtb: given a turn,
    it returns that torsion and stops looking."""
    from delfin.dashboard import gfn_optimize as gfn

    held = gfn.contacts_holding(_ETHANE, [0, 2, 3, 4], turning=[2, 0, 1, 5])
    assert len(held) == 1, held
    assert held[0]['kind'] == 'dihedral'
    assert held[0]['atoms'] == [2, 0, 1, 5]


def _an_editor_watching_js(text):
    """An editor whose page script is captured rather than run."""
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    said = []
    ctx.run_js = said.append
    state = {}
    part = structure_editor.build(
        ctx, state=state, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    return part, state, said


def test_the_page_is_told_when_dynamik_opt_goes_on():
    """The half of this that was missing, and the whole of the report.

    The page has a switch it can see for itself -- ``autoOpt`` -- and it is
    not this question: it is the *browser's* own field, which Dynamik Opt
    switches off under every server method so that the molecule can follow the
    hand through xtb instead. A right button reading it was told "no" in
    exactly the mode most of this editor is used in, fell through to the pivot,
    and a pivot rotate turns the selection: "drehmoment geht doch aber nur wenn
    alles ausgewaehlt ist ... ich kann nicht nur ein atom nehmen".
    """
    part, _state, said = _an_editor_watching_js(_ETHANE)
    part.submit_relax_btn.value = True
    on = [one for one in said if 'setLiveHand' in one]
    assert on, 'the page was never told the hand is live'
    assert 'true' in on[-1], on[-1]

    said.clear()
    part.submit_relax_btn.value = False
    off = [one for one in said if 'setLiveHand' in one]
    assert off and 'false' in off[-1], off
