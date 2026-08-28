"""Behavioural invariants of the Submit-tab molecule editor.

These lock in the interaction model that was verified by driving the real
editor in a headless browser: one toggle per click, a rubber band that projects
instead of dispatching synthetic events, direct atom dragging, and a camera
that can still be turned while manipulating.
"""
import re

import pytest

from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js
from editor_source import FULLSCREEN_CSS
from editor_source import EDITOR_SOURCE as _EDITOR_PY
from editor_source import SUBMIT_SOURCE
from editor_source import SUBMIT_SOURCE as _TAB_AND_EDITOR

EDITOR = submit_manip_bootstrap_js()

#: Planar benzene, atom order C0..C5 then the hydrogens.  C0 and C3 sit across
#: the ring from each other, 2.795 A apart -- far enough that a bond drawn
#: between them is unmistakably not at its equilibrium length.
_BENZENE = """12
benzene
C  1.3970  0.0000  0.0000
C  0.6985  1.2098  0.0000
C -0.6985  1.2098  0.0000
C -1.3970  0.0000  0.0000
C -0.6985 -1.2098  0.0000
C  0.6985 -1.2098  0.0000
H  2.4810  0.0000  0.0000
H  1.2405  2.1486  0.0000
H -1.2405  2.1486  0.0000
H -2.4810  0.0000  0.0000
H -1.2405 -2.1486  0.0000
H  1.2405 -2.1486  0.0000
"""


def _body(name):
    """Return the source of one JS function in the editor bootstrap."""
    start = EDITOR.index(f'function {name}(')
    depth, i = 0, EDITOR.index('{', start)
    for j in range(i, len(EDITOR)):
        if EDITOR[j] == '{':
            depth += 1
        elif EDITOR[j] == '}':
            depth -= 1
            if depth == 0:
                return EDITOR[start:j + 1]
    raise AssertionError(f'unbalanced braces in {name}')


def test_a_click_toggles_a_pick_exactly_once():
    """3Dmol's picker and the projection picker used to be two independent
    click handlers that each toggled, so clicking an atom selected and
    immediately deselected it."""
    clickable = _body('attachClickable')
    assert 'togglePick' not in clickable, (
        'the native picker must only record its hit, not toggle')
    assert '_nativeHit' in clickable

    fallback = _body('installCanvasClickFallback')
    assert fallback.count('togglePick') == 1
    # Deferred, so 3Dmol's own picker has had its chance to record a hit first.
    assert 'setTimeout' in fallback
    assert 'raycastAtom' in fallback


def test_rubber_band_projects_instead_of_dispatching_synthetic_events():
    """The band used to dispatch a mousedown/mouseup/click triple every 7 px,
    which froze the tab and scrambled the selection if Shift came up early."""
    finish = _body('finishRect')
    assert 'dispatchEvent' not in finish
    assert 'new MouseEvent' not in finish
    assert 'projectAllAtoms' in finish
    # A non-additive band replaces the selection instead of only ever growing it.
    assert 'if (!additive) state.picks = [];' in finish


def test_picking_prefers_the_atom_nearest_the_camera():
    """Screen-distance-only picking selects atoms hidden behind the structure
    in fused rings and crowded coordination spheres."""
    raycast = _body('raycastAtom')
    assert 'depth' in raycast
    assert 'bestDepth' in raycast
    project = _body('projectWithDepth')
    assert 'projectionMatrix' in project
    assert 'matrixWorldInverse' in project


def test_manipulate_mode_supports_direct_atom_drag():
    """Avogadro's gesture: press on an atom and drag it, with no prior
    selection round-trip through Select mode."""
    assert "kind: 'translate'" in EDITOR
    assert 'targets:' in EDITOR
    assert 'grabbed:' in EDITOR
    # Grabbing a selected atom moves the whole selection.
    assert 'inSelection' in EDITOR
    # And the drag acts on the drag's own target set, not blindly on picks.
    assert 'applyTranslate(scopeKey, delta, d.targets)' in EDITOR


def test_manipulate_mode_still_lets_the_camera_be_turned():
    """3Dmol binds its mouse handlers on the canvas and the overlay is the
    canvas's sibling, so a press on empty space has to be handed over
    explicitly — bubbling cannot reach it."""
    assert '_handleMouseDown' in EDITOR
    assert re.search(r'v\._handleMouseDown\(e\)', EDITOR)


def test_editor_ships_no_debug_logging():
    """One of these fired on every mousemove of a rotate drag."""
    assert 'console.log' not in EDITOR


def test_dragging_an_atom_relaxes_the_rest_through_the_force_field():
    """Avogadro's manipulation: the grabbed atom follows the cursor and the
    molecule settles around it. The relaxation runs in the browser because a
    per-frame round trip to the kernel costs 45 ms.

    Exactly, when the hand is set to rigid; otherwise it is pulled and the
    field decides how far it gets -- see the pull test below."""
    # Frozen atoms are a gradient mask on the engine side, set once per drag.
    assert 'ffBeginDrag(scopeKey, state.drag.targets)' in EDITOR
    assert 'ffApplyFrozen(scopeKey, ffIndicesOf(viewer, targets))' in EDITOR
    # The relaxation runs after the grabbed atoms are placed, not before.
    grab_then_relax = EDITOR.index('applyTranslate(scopeKey, delta, d.targets)')
    assert EDITOR.index('ffRelaxAsync(scopeKey', grab_then_relax) > grab_then_relax
    # Releasing is now expressed as re-applying the frozen set with nothing
    # extra held, so a drag and a pinned atom go through one code path.
    assert 'ffApplyFrozen(scopeKey, [])' in _body('ffEndDrag')

    # The engine adapts its batch to a full-frame budget, so the time it is
    # given has to include the renderer's share of the frame.
    relax = _body('ffRelaxFrame')
    assert 'state.ffFrameMs = nowMs() - t0;' in relax

    # Parameters are assigned for one geometry; a re-render invalidates them.
    ready = _body('onViewerReady')
    assert 'state.ffActive = false;' in ready


def test_force_field_is_switched_on_from_python_not_polled():
    from delfin.dashboard import tab_submit
    source = _EDITOR_PY
    assert 'def _enable_live_forcefield' in source
    assert 'export_forcefield_terms' in source
    assert 'setForceField' in source
    # Assigning parameters is a one-off; nothing may talk to Python mid-drag.
    assert 'Runs once, when the toggle is switched on' in source


def test_internal_coordinates_move_the_far_fragment():
    """Setting a bond, angle or dihedral has to move a whole fragment, not one
    atom: which one follows from the model's own connectivity, cutting the bond
    the coordinate turns about."""
    assert 'function setInternal' in EDITOR
    assert 'function readInternal' in EDITOR
    frag = _body('fragmentFrom')
    assert 'cutA' in frag and 'cutB' in frag

    body = _body('setInternal')
    # Every edit is undoable and reaches the coordinate box; the geometry work
    # itself lives in applyInternalValue, shared with a fixed constraint.
    assert 'snapshotForUndo(scopeKey)' in body
    assert "pushXyzToPython(scopeKey, 'drag-end')" in body
    assert 'applyInternalValue(' in body
    body = _body('applyInternalValue')
    # A ring keeps both halves attached, so the caller is told rather than
    # silently tearing the ring open.
    assert 'ring: only the second atom moved' in body
    assert 'that dihedral turns about a ring bond' in body
    # The sign of a rotation is verified against the value actually reached.
    assert 'if (Math.abs(angleV(' in body


def test_force_field_choice_is_honest_about_what_it_changes():
    """MMFF94's bond, angle, torsion and van der Waals forms are all different
    from the ones the browser engine evaluates, so selecting it must not just
    relabel UFF terms."""
    from delfin.dashboard import molecule_forcefield as mff

    assert mff.normalise_method('MMFF94') == 'mmff94'
    assert mff.normalise_method('mmff-94') == 'mmff94'
    assert mff.normalise_method(None) == 'uff'
    assert mff.normalise_method('nonsense') == 'uff'

    source = open(mff.__file__, encoding='utf-8').read()
    assert 'Silently relabelling' in source
    assert 'MMFFGetMoleculeForceField' in source
    assert 'MMFF94 has no transition-metal parameters' in source


def test_optimisation_covers_every_frame_and_is_undoable():
    """The Submit tab can hold a whole set at once -- generated isomers or the
    frames of a batch -- and any of them can end up submitted, so optimising
    only the visible one would leave the rest untouched. The browser's undo
    stack cannot cover this: results arrive from Python and re-render."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    body = source.split('def on_submit_optimize(change=None, every_frame=False)')[1].split('\n    def ')[0]
    assert "frames = list(list_structures() or []) if every_frame else []" in body, (
        "Optimize takes the frame on screen; all takes the set"
    )
    assert "state['pre_optimize_frames']" in body
    assert 'relax_xyz(' in body and 'max_steps=500' in body
    # A 500-step minimisation per frame takes seconds; it must not block the UI
    # -- and it goes on a thread through the one place that counts the work, so
    # the spinner is up for the whole of it and comes down even if it raises.
    assert '_start_background(_work,' in body
    # The thread is kept where the next press can find it: an interrupted run
    # is still shutting its process down when the replacement starts, and a
    # login node is shared.
    assert "remember_in='optimize_thread'" in body
    # One bad frame must not lose the others.
    assert 'results.append(item)' in body

    undo = source.split('def on_submit_manip_undo')[1].split('\n    def ')[0]
    # One history now: the geometry before an optimisation is an entry in it
    # like any other step, rather than a slot of its own that Undo had to look
    # in first and that nothing could get past.
    assert "state.pop('pre_optimize_frames', None)" in undo
    assert '_undo_structure()' in undo
    assert "_remember('an optimisation')" in source


def test_optimize_toggle_runs_the_field_continuously():
    """Avogadro's auto optimisation keeps the force field running on a timer,
    so the structure settles whether or not the mouse is down and a grabbed
    atom drags the relaxing molecule along."""
    assert 'function startAutoOptimize' in EDITOR
    assert 'function stopAutoOptimize' in EDITOR
    tick = _body('autoOptimizeTick')
    # A drag already relaxes in its own handler; a second batch here would
    # double the step budget for that frame.
    assert 'if (state.drag) {' in tick
    assert tick.index('if (state.drag) {') < tick.index('ffRelaxAsync')
    assert 'window.requestAnimationFrame' in tick
    # The coordinate box follows at a readable rate -- each push is a round trip.
    assert "state.autoPushed" in tick

    stop = _body('stopAutoOptimize')
    assert 'cancelAnimationFrame' in stop
    # Named like the heartbeat, and for the same reason: this is the last word
    # of the field that has just been switched off, not an edit -- and a push
    # that cannot say which it is cannot be kept off a running optimisation.
    assert "pushXyzToPython(scopeKey, 'field')" in stop

    # Starting is undoable, and a re-render must not leave a loop spinning on
    # a viewer that no longer exists.
    assert 'snapshotForUndo(scopeKey);\n        state.autoOpt = true;' in EDITOR
    ready = _body('onViewerReady')
    assert 'state.autoOpt = false;' in ready


def test_ctrl_z_belongs_to_whatever_is_being_typed_in():
    """The editor took Ctrl-Z globally, so undoing a typo while editing the
    coordinate box silently moved atoms instead. Delete, added later for
    Unbond, would have cut a bond on a backspace in the same box."""
    guard = _body('typingInAField')
    assert 'document.activeElement' in guard
    assert "tag === 'INPUT' || tag === 'TEXTAREA'" in guard
    assert 'focused.isContentEditable' in guard

    handler = EDITOR.split("on(window, 'keydown'")[1].split('}, true);')[0]
    # Every shortcut asks first, before it looks at any molecule.
    for shortcut in ("key === 'z'", "key === 'Delete'"):
        after = handler[handler.index(shortcut):]
        assert after.index('typingInAField()') < after.index('_submitManipStateByScope')


def test_value_box_says_what_it_sets_and_shows_the_current_value():
    """It was an unlabelled number field: nothing said it wanted a bond length,
    and it did not show what the selection currently measures."""
    readout = _body('updateInternalReadout')
    assert '.submit-internal-label' in readout
    assert '.submit-internal-value input' in readout
    assert 'pick 2-4 atoms' in readout
    # Angstrom for a bond, degrees for an angle or dihedral.
    assert "info.unit === 'A'" in readout
    # Typing must not be overwritten under the user's fingers.
    assert 'document.activeElement !== box' in readout
    # And it has to refresh whenever the selection changes.
    status = _body('updateStatus')
    assert 'updateInternalReadout(scopeKey)' in status
    # The status line carries state, not instructions: gesture hints belong in
    # the tooltips, and the value box already labels itself.
    assert 'turns the view' not in status
    assert 'press <b>Set</b>' not in status


def test_optimise_starts_a_second_row_in_fullscreen():
    """A screen-wide toolbar put Optimise at the far end of one long line.

    Flexbox cannot be told to break, so the break is an element that takes a
    whole line and no height, sitting in front of Optimise.  Measured in a
    browser: at 1440 wide inside the overlay the toolbar is two rows, the
    second beginning with Optimise; at 760 wide outside it the element is
    display:none and the toolbar wraps where it always did -- three rows, not
    four, so the break costs nothing where it is not wanted.
    """
    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    children = source.split('submit_manip_toolbar = widgets.HBox')[1].split(']')[0]
    # Anchored on the method box rather than on the strength slider: the
    # sliders have gone onto the picture, where the eye already is, and what
    # this asks is only that the break comes before Optimise.
    order = ['submit_ff_dd', 'submit_fs_row_break',
             'submit_optimize_btn', 'submit_optimize_all_btn']
    positions = [children.index(name) for name in order]
    assert positions == sorted(positions), children

    made = source.split('submit_fs_row_break = widgets.Box')[1].split(')\n')[0]
    assert "display='none'" in made, 'it must be inert outside the overlay'
    # The rule is in the shared sheet: the break has to work in the ORCA
    # Builder's overlay too, which holds the same toolbar.
    rule = FULLSCREEN_CSS.split(
        '.delfin-structure-fs-overlay .submit-fs-row-break {')[1]
    rule = rule.split('}')[0]
    assert 'display: block' in rule
    assert 'flex: 1 0 100%' in rule
    assert 'height: 0' in rule


def test_toolbar_wraps_instead_of_clipping_its_own_controls():
    """On a laptop the row was wider than the panel, and nowrap plus
    overflow hidden simply cut off whatever did not fit."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    toolbar = source.split('submit_manip_toolbar = widgets.HBox')[1].split(')\n')[0]
    assert "flex_flow='row wrap'" in toolbar
    assert "overflow='hidden'" not in toolbar
    # The status line is not on this row at all any more: what is picked is a
    # fact about the structure, so it lies in the picture's top-left corner
    # with the run's own sentence along the bottom and the view controls in
    # the other corner. On the toolbar it took a share of a row it had nothing
    # to do with, and on a laptop it wrapped to a line of its own.
    status = source.split('submit_manip_status = widgets.HTML')[1].split(')\n')[0]
    assert "flex='1 1 260px'" not in status
    assert "delfin-structure-picks-over" in source

    # The label, the value box and Set are still one group, in that order: a
    # label that lands on a different row from its field explains nothing.
    #
    # This used to be enforced by giving the group ``flex_flow='row nowrap'``,
    # and that is the line the toolbar's overflow came out of: the toolbar
    # around it wraps, a wrapping container breaks between its items and never
    # inside one, and the group is one item -- so with a scan armed its
    # nineteen controls were laid out on a single line about 1900 px wide,
    # inside a row 620 px wide at a 1280 px window, and the last of them (Path
    # from here, Find the path, Path to saddle) were off the screen with no way
    # to reach them. The declaration bought nothing either: measured in
    # chromium with the group wrapping, at six widths from 1920 to 800 px, in
    # three states and in both the tab and the fullscreen overlay, the label,
    # its value box and Set are on the same row in all seventy-two -- they lead
    # the group, so they lead a line. Where they end up is measured in
    # tests/test_the_toolbar_stays_on_the_screen.py; the order is what is
    # asked here.
    group = source.split('submit_internal_group = widgets.HBox')[1].split(')\n')[0]
    for name in ('submit_internal_label', 'submit_internal_value',
                 'submit_internal_btn', 'submit_hold_btn'):
        assert name in group, name
    assert "flex_flow='row nowrap'" not in group

    # Force field, then its strength, then the one-shot run, then the
    # continuous one, then the internal-coordinate group.
    children = source.split('submit_manip_toolbar = widgets.HBox')[1].split(']')[0]
    # The strength slider is no longer on this row -- it acts on the feel of
    # the hand rather than on the structure, and it lives on the picture with
    # the others of its kind (see submit_view_panel).
    order = [
        'submit_ff_dd',
        'submit_optimize_btn', 'submit_relax_btn', 'submit_internal_group',
    ]
    positions = [children.index(name) for name in order]
    assert positions == sorted(positions), children


def test_relaxation_strength_is_adjustable():
    """At full strength the field moves the structure so far per frame that a
    dragged atom has to be fought rather than led."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'submit_strength_slider' in source
    assert 'setOptimizerStrength' in source

    strength = _body('setOptimizerStrength')
    assert 'maxChunk' in strength
    # Re-applied whenever parameters are assigned, so it survives a reload.
    assert 'setOptimizerStrength(scopeKey, state.ffStrength)' in EDITOR

    from delfin.dashboard import molecule_forcefield_js as ffjs
    engine = ffjs.MOLECULE_FF_BOOTSTRAP_JS
    assert 'st.maxChunk' in engine
    assert 'opts.maxChunk' in engine


def test_camera_survives_a_rebuild_of_the_same_structure():
    """Optimising or stepping to another isomer rebuilds the viewer, and the
    view snapped back to the default orientation every time."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    assert '_submitMolViewByScope' in source
    assert 'prev.getView()' in source
    assert 'viewer_UNIQUEID.setView(saved.view)' in source
    # Only for the same structure: a different molecule deserves a fresh look.
    assert 'saved.atoms === count' in source


def test_panning_calls_the_function_that_actually_moves_the_scene():
    """viewer.translate() exists in this 3Dmol build but does nothing;
    translateScene() is the call that moves it. Probing translate first meant
    the working call was never reached, so right-drag never panned at all."""
    from delfin.dashboard import molecule_viewer as mv

    body = mv.RIGHT_MOUSE_TRANSLATE_PATCH_JS.split('var translateNow')[1].split('};')[0]
    assert body.index('translateScene') < body.index('viewer.translate(')


def test_editing_gets_room_so_atoms_are_not_clipped_or_fogged_away():
    """The default slab is about the size of the loaded molecule, so an atom
    pulled towards the viewer crosses the near plane: it washes out and then
    disappears. Measured brightness went from 210.5 (nearly gone) back to
    147.4, against 145.7 at rest, once the slab was widened."""
    widen = _body('widenSlabForEditing')
    assert 'viewer.setSlab(-EDIT_SLAB, EDIT_SLAB)' in widen
    # The viewing slab is restored when editing ends.
    assert 'state.slabSaved' in widen
    assert 'viewer.setSlab(state.slabSaved.near, state.slabSaved.far)' in widen
    mode = _body('setMode')
    assert "widenSlabForEditing(scopeKey, state.mode !== 'off')" in mode


def test_the_wheel_still_zooms_while_a_mode_is_active():
    """The overlay covers the canvas while editing and 3Dmol listens for the
    wheel on the canvas, so scrolling stopped zooming the moment a mode was
    switched on."""
    assert "ov.addEventListener('wheel'" in EDITOR
    assert 'viewer._handleMouseScroll(e)' in EDITOR
    assert '{passive: false}' in EDITOR


def test_right_drag_pans_while_the_field_is_running():
    """With the structure moving under the cursor, shifting the view is what
    the user reaches for -- not a pivot rotation."""
    assert 'if (s.autoOpt) continue;' in EDITOR
    assert 'if (state.autoOpt) return;' in EDITOR
    # The context menu stays suppressed either way, or panning pops a menu.
    assert "if (s.mode === 'manipulate' || s.autoOpt) {" in EDITOR


def test_hydrogens_can_be_grabbed_in_their_own_right():
    """Picking took the atom nearest the camera within a flat 14 px radius, so
    a hydrogen was routinely stolen by the fatter carbon bonded to it — the
    cursor sat on the hydrogen and the carbon won anyway.

    Each atom is now only a candidate where its own drawn disc is under the
    cursor, and depth decides between the candidates. Measured on cholesterol:
    clicking precisely on each of its 66 hydrogens grabs 60 of them, and all
    six exceptions are genuinely covered by an atom in front."""
    body = _body('raycastAtom')
    assert 'elementRadius(q.atom)' in body
    assert 'DEFAULT_ATOM_SCALE' in body
    assert 'getPixelToWorld' in body
    # Depth still shields a covered atom, but only among real candidates.
    assert 'q.depth < bestDepth' in body
    # A hydrogen must be much smaller than a carbon or the distinction is moot.
    assert EDITOR.index('H: 1.10') > 0
    radii = _body('elementRadius')
    assert 'DEFAULT_VDW' in radii


def test_toolbar_parts_are_found_even_when_fullscreen_moves_them():
    """Fullscreen lifts the toolbar into a floating overlay, outside the tab's
    scope container. A scope-only lookup then found nothing: the value box
    stayed empty and, worse, edits never reached Python at all, because the
    sync input had left the scope too.

    The overlay carries the scope's class as well, so the answer is to look in
    every element that has it. It used to fall back to the whole page instead,
    on the grounds that there is one Submit tab per dashboard -- which stopped
    being true when the ORCA Builder got an editor of its own.
    """
    finder = _body('findInScope')
    assert "document.querySelectorAll('.' + scopeKey)" in finder
    assert 'roots[i].querySelector(selector)' in finder
    assert 'document.querySelector(selector)' not in finder
    for name in ('getSyncInput', 'getStatusEl', 'updateInternalReadout'):
        body = _body(name)
        assert 'findInScope(' in body, name
        assert 'getRoot(scopeKey)' not in body, name


def test_empty_space_belongs_to_the_viewer_for_both_buttons():
    """In Manipulate mode the buttons split the same way on empty space as on
    an atom: left turns the view, right pans it. The editor only takes the
    right button where it lands on an atom, to set the pivot.

    Measured in a browser: left on empty space rotates and does not pan, right
    on empty space pans and does not rotate, right on an atom sets the pivot
    and moves neither."""
    window_steal = EDITOR.split("on(window, 'mousedown'")[1].split('}, true);')[0]
    # The probe has to happen before the event is taken, not after.
    assert window_steal.index('probeClickAtom') < window_steal.index('e.preventDefault()')
    assert 'if (!picked) continue;' in window_steal

    overlay = EDITOR.split("ov.addEventListener('mousedown'")[1]
    # The manipulate branch specifically: draw mode has its own right button
    # now, and it means something else there.
    manipulate = overlay.split("if (state.mode === 'manipulate') {")[1]
    right = manipulate.split('if (e.button === 2) {')[1][:3000]
    assert right.index('probeClickAtom') < right.index('e.preventDefault();')
    # Two branches now: with the field running the right button is a torque on
    # a bond, and without it the pivot it always was. Both probe for an atom
    # first and both hand the press back when there is none, because empty
    # space belongs to the viewer whatever the button does on atoms.
    assert 'if (!picked) return;' in right
    assert 'if (!turning) return;' in right
    # Asked of the switch that means "a hand here will be answered", which
    # is not the browser field's own flag: under every server method Dynamik
    # Opt switches that one off, so a right button reading it fell through to
    # the pivot in the mode the editor is mostly used in.
    turning = right.split('if (handIsLive(state)) {')[1]
    assert turning.index('probeClickAtom') < turning.index('e.preventDefault()')


def test_bonding_is_perceived_once_and_not_re_read_from_a_dragged_geometry():
    """Bond orders are read from the geometry, and a twisted double bond stops
    looking like one. Measured on stilbene: turning the central C=C from 180
    to 150 degrees is already enough for perception to report zero double
    bonds, which turns the torsion from two-fold with a 19.5 kcal/mol barrier
    into three-fold with 1.1 — a free single bond. Nothing then pulls the
    double bond back to planar, which is exactly what a user dragging a
    stilbene sees.

    Reusing the perception taken from the structure as loaded keeps n = 2 and
    19.5 kcal/mol at any twist."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'def _perception_for' in source
    assert "state.get('perceived_for') == fingerprint" in source
    # Both the live parameters and the release-time relaxation use it.
    assert 'perceived=_perception_for(xyz)' in source
    assert source.count('perceived=_perception_for(xyz)') >= 2
    # A genuinely different molecule must be perceived afresh.
    assert 'def _structure_fingerprint' in source


def test_remembered_bonding_is_dropped_when_the_structure_really_changes():
    """The element sequence cannot tell two constitutional isomers apart, so
    the cache must not rest on it alone. Every genuine change — a paste, a
    conversion, an isomer step, an optimisation result — comes through
    update_molecule_view, and a drag does not: it takes the manip_inflight
    branch, which is what keeps a dragged double bond a double bond."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    view = source.split('def update_molecule_view')[1].split('\n    def ')[0]
    assert "state['perceived'] = None" in view
    # It must sit after the inflight early return, or a drag would clear it.
    assert view.index("state.get('manip_inflight')") < view.index("state['perceived'] = None")


def test_letting_go_settles_the_structure():
    """The dragged atom was left wherever the cursor stopped, and the strain
    the drag put in reached the coordinate box and from there the calculation.
    Measured on cholesterol at the shipped strength: 176 kcal/mol above what a
    relaxation gives, and 1447 kcal/mol above after a harder pull."""
    end = _body('ffEndDrag')
    assert 'settleAfterDrag(scopeKey)' in end

    settle = _body('settleAfterDrag')
    # The continuous optimiser is already doing this; do not run both.
    assert 'if (state.autoOpt) return;' in settle
    # Stops on convergence, on a stall, or on a bounded number of frames.
    assert 'stats.converged || stats.stalled' in settle
    assert 'SETTLE_MAX_FRAMES' in settle
    # And the settled geometry is what reaches Python.
    assert "pushXyzToPython(scopeKey, 'drag-end')" in settle
    # A new drag takes over from a settle in progress.
    assert 'stopSettling(scopeKey);' in EDITOR


def test_the_viewer_shows_the_energy():
    badge = _body('updateEnergyBadge')
    assert 'kcal/mol' in badge
    # Shown from the moment the field is loaded, before anything is relaxed.
    assert '__delfinFF.energy(scopeKey, ffReadPositions(viewer))' in badge
    # And it says which parameters produced it.
    assert 'state.ffInfo' in badge
    ensure = _body('ensureEnergyBadge')
    assert "badge.className = 'submit-energy-badge'" in ensure
    assert 'pointer-events:none' in ensure
    # A re-render replaces the viewer element, so the badge must be forgotten.
    assert 'state.energyBadge = null;' in _body('onViewerReady')


def test_settling_on_release_can_be_switched_off():
    """Placing an atom somewhere and having it stay there is sometimes exactly
    what is wanted, strain and all. Measured on cholesterol: with the switch
    on, 1539 kcal/mol at release settles to 126.5; with it off, 1598.6 stays
    1598.6."""
    end = _body('ffEndDrag')
    assert 'state.settleOnRelease === false' in end
    # Switching off still pushes the geometry the user placed.
    assert end.index('settleOnRelease === false') < end.index('settleAfterDrag')
    assert "pushXyzToPython(scopeKey, 'drag-end');" in end

    setter = _body('setSettleOnRelease')
    assert 'stopSettling(scopeKey)' in setter

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'submit_settle_btn' in source
    assert 'setSettleOnRelease' in source
    # On by default: it is what keeps a submitted geometry sane.  Only under
    # the browser's field, though -- the toolbar takes it away under a server
    # method, where what it ran was the whole minimisation and not a tidy-up.
    settle = source.split('submit_settle_btn = widgets.ToggleButton')[1].split(')\n')[0]
    assert 'value=True' in settle
    # And the choice survives a re-assignment of the parameters.
    assert source.count('setSettleOnRelease') >= 2


def test_placed_atoms_are_held_against_the_running_field():
    """With settling off, releasing the atom was not enough: the continuous
    field simply relaxed it back, so the switch did nothing whenever the field
    was on — the combination a user is most likely to have. Placed atoms stay
    frozen in the gradient mask instead, and the rest of the molecule goes on
    settling around them.

    Measured on cholesterol with the field running: the held atom moves
    0.0000 A while the rest moves up to 5.41 A. With settling on, the same
    gesture moves it 3.28 A."""
    end = _body('ffEndDrag')
    assert 'state.settleOnRelease === false' in end
    # The mouseup handler clears state.drag before ending the drag, so the
    # held atoms have to be handed over rather than read back -- and only
    # when the press actually moved something.  A press that did not move
    # is a tap, which selects the atom; pinning it as well would be one
    # gesture doing two things, and it would freeze every atom the user
    # named.
    assert 'ffEndDrag(scopeKey, d.movedEnough ? d.targets : [])' in EDITOR
    assert 'heldSerials' in end
    assert 'state.pinned = pinned;' in end

    apply_frozen = _body('ffApplyFrozen')
    assert 'state.pinned' in apply_frozen
    assert '__delfinFF.grab(scopeKey, list)' in apply_frozen
    # A drag freezes what it holds on top of what is already held.
    assert 'ffApplyFrozen(scopeKey, ffIndicesOf(viewer, targets))' in EDITOR

    # Three ways to let go again: Clear, switching settling back on, a re-render.
    assert 'unpinAll(scopeKey);' in _body('clearPicks')
    assert 'unpinAll(scopeKey);' in _body('setSettleOnRelease')
    assert 'state.pinned = [];' in _body('onViewerReady')
    # And the count is visible, or the mode is invisible.
    assert "'</b> held'" in _body('updateStatus')


def test_tapping_a_metal_offers_its_coordination_polyhedra():
    """The editor reports the selection to Python as model indices, so a metal
    tapped in Select mode can be offered the polyhedra its coordination number
    allows -- from MANTA's own ideal-donor tables."""
    push = _body('pushPicksToPython')
    assert '.submit-pick-sync' in push
    # Model indices, not 3Dmol serials: the payload is indexed the same way.
    assert 'ffIndicesOf(viewer, serials)' in push
    # Sent through the native setter, the way ipywidgets notices a change.
    assert "dispatchEvent(new Event('input', {bubbles: true}))" in push
    assert 'pushPicksToPython(scopeKey);' in _body('updateStatus')

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'submit_pick_sync' in source
    assert 'polyhedron_options(perceived, indices[0])' in source
    # Only for a single selected atom, and only when it is a metal.
    assert 'len(indices) == 1' in source
    # Choosing one re-assigns the parameters, which is what starts the pull.
    changed = source.split('def on_submit_poly_changed')[1].split('\n    def ')[0]
    assert '_enable_live_forcefield()' in changed
    # Repopulating the list must not look like a user choice.
    assert "state.get('poly_quiet')" in changed


def test_polyhedra_are_offered_without_the_force_field_being_on_first():
    """The offer hung on the cached perception, which only exists once the
    force field has been switched on. Tapping a metal before that did nothing
    at all — and said nothing either."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    handler = source.split('def on_submit_pick_sync')[1].split('\n    def ')[0]
    assert '_perception_for(xyz)' in handler
    assert "state.get('perceived')" not in handler
    # And when nothing can be offered for a metal, the reason is shown.
    assert 'no ' in handler and 'polyhedron table' in handler
    assert 'coordination number' in handler


def test_values_can_be_held_and_dropped_again():
    """Set was a one-off edit. Holding a value means the field works against it
    for as long as it is listed, alongside a polyhedron if one is set.

    Measured on a Ni complex: square planar alone reaches 3.7 degrees of
    deviation, a ligand angle held alone pulls 126.6 to 111.9 degrees, and both
    together give 3.3 degrees and 111.8 -- with the worst ligand bond changing
    0.1377 A against 0.1344 for a plain relaxation."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'submit_hold_btn' in source
    hold = source.split('def on_submit_hold')[1].split('\n    def ')[0]
    assert "_CONSTRAINT_KINDS.get(len(indices))" in hold
    # Re-holding the same atoms replaces rather than stacks.
    assert "held = [c for c in held if c['atoms'] != indices]" in hold
    assert '_enable_live_forcefield()' in hold

    # Everything the field is being held to is listed, the polyhedron included.
    refresh = source.split('def _refresh_constraints')[1].split('\n    def ')[0]
    assert "state.get('poly_applied')" in refresh
    assert "state.get('constraints')" in refresh
    # And each entry can be dropped again.
    drop = source.split('def on_submit_constraint_del')[1].split('\n    def ')[0]
    assert "state['poly_applied'] = None" in drop
    assert 'held.pop(position)' in drop
    assert 'restraints=[' in source


def test_an_exchange_lets_the_polyhedron_hand_over_the_vertices():
    """The vertex bookkeeping was swapped by hand once, which the field then
    fought: it had no way to cross into the other arrangement. The exchange
    moves the ligands instead, and the polyhedron works its vertices out
    afresh from where they have landed."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    handler = source.split('def on_submit_swap')[1].split('\n    def ')[0]
    # No hand-written assignment any more; it is re-derived.
    assert "state['poly_assignment'] = None" in handler
    assert 'assignment[first], assignment[second]' not in handler
    assert 'exchangeLigands(' in handler


def test_a_held_polyhedron_survives_selecting_something_else():
    """The offer follows the selection; what has been applied must not.
    Clearing the metal whenever the selection was not a single metal meant
    that picking three ligand atoms to hold an angle silently threw the
    polyhedron away, and the export that followed went out without it."""
    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    handler = source.split('def on_submit_pick_sync')[1].split('\n    def ')[0]
    assert "state['poly_offer_metal'] = None" in handler
    assert "state['poly_metal'] = None" not in handler
    # Applying a geometry is what binds it to a metal, and removing it unbinds.
    changed = source.split('def on_submit_poly_changed')[1].split('\n    def ')[0]
    assert "state['poly_metal'] = state.get('poly_offer_metal')" in changed
    drop = source.split('def on_submit_constraint_del')[1].split('\n    def ')[0]
    assert "state['poly_metal'] = None" in drop
    # A different structure makes every stored index meaningless.
    view = source.split('def update_molecule_view')[1].split('\n    def ')[0]
    assert "state['constraints'] = []" in view
    assert "state['poly_applied'] = None" in view


def test_dragging_a_ligand_hands_it_the_vertex_it_lands_on():
    """Dragging a ligand towards another position and having the field haul it
    straight back is not an exchange, it is a fight. The assignment is worked
    out again when the drag ends, so the polyhedron accepts the ligand where it
    was put.

    Measured on a Ni complex: dragging two ligands onto each other's side
    exchanges their vertices (2 to 0 and 0 to 2), and after relaxing the moved
    ligand sits 0.38 A from where it was dropped rather than the 2.60 A back at
    its old place."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    sync = source.split('def on_submit_manip_sync')[1].split('\n    def ')[0]
    assert "state['poly_assignment'] = None" in sync
    assert 'schedule_ui_update(_enable_live_forcefield)' in sync


def test_a_typed_value_is_not_overwritten_by_the_running_field():
    """The value box refreshed on every status update, and the continuous
    optimiser triggers one per frame. So a value typed in was replaced by the
    current measured one within a frame, and Set or Hold then acted on the
    geometry instead of on what was asked for — which looked exactly like the
    buttons doing nothing at all in dynamic mode.

    Verified in a browser: 180 typed into the box is still 180 after 2.5
    seconds of relaxation, and changing the selection gives a fresh reading."""
    readout = _body('updateInternalReadout')
    assert 'state.readoutFor' in readout
    assert 'selectionChanged' in readout
    assert 'selectionChanged && document.activeElement !== box' in readout
    # A re-render must not make the next selection look unchanged.
    assert 'state.readoutFor = null;' in _body('onViewerReady')


def test_the_polyhedron_reconsiders_once_per_drag_not_twice_a_second():
    """Reassigning donors to vertices was hung on the coordinate push, and the
    running optimiser pushes twice a second as a heartbeat. So the whole field
    was re-exported and reloaded every 500 ms, and a ligand dragged towards
    another vertex never got the chance to settle on it.

    Verified in a browser: three seconds of relaxation produce six pushes and
    no drag-end marker at all, and one drag produces exactly one."""
    push = _body('pushXyzToPython')
    assert "var note = reason ? ('DELFIN ' + reason) : null;" in push
    assert "serializeXyz(viewer, note,\n" in push
    # The heartbeat inside the relaxation loop is the field saying where it
    # has got to, and it is named that.  It carried no reason at all once,
    # which made it indistinguishable on the kernel side from anything else
    # that was not a drag -- so it could not be told apart from an edit, and
    # it wrote the coordinate box over a running optimisation.
    tick = _body('autoOptimizeTick')
    assert "pushXyzToPython(scopeKey, 'field')" in tick
    assert "drag-end" not in tick
    # A drag that ends does, whether it settles or not.
    assert "pushXyzToPython(scopeKey, 'drag-end')" in _body('ffEndDrag')
    assert "pushXyzToPython(scopeKey, 'drag-end')" in _body('settleAfterDrag')

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    sync = source.split('def on_submit_manip_sync')[1].split('\n    def ')[0]
    # The comment line now carries the atoms the hand is on as well, so the
    # marker is read as a prefix rather than as the whole line.
    assert "note.startswith('DELFIN drag-end')" in sync
    # And it runs after the coordinates have landed, or the assignment would be
    # worked out from where the ligands used to be.
    assert sync.index("state['poly_recheck'] = True") < sync.index('coords_widget.value = payload')
    assert "state.pop('poly_recheck', False)" in sync


def test_force_field_notes_sit_under_the_structure_they_describe():
    """What the field had to approximate — a metal with no RDKit parameters, a
    polyhedron being forced — was written into the preview's status line, where
    conversion messages overwrite it and it scrolls away from the structure it
    is about."""
    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    assert 'submit_ff_notes' in source
    assert "submit_ff_notes.add_class('submit-ff-notes')" in source
    # Below the copy row, in the panel the viewer lives in. Read off the box
    # itself rather than the first line that happens to name those three: the
    # order in the overlay is the order of these children, so this is the
    # thing worth pinning.
    children = source.split('submit_right = widgets.VBox([')[1].split('])')[0]
    # mol_viewer_stack is the picture: the viewer with the status line lying
    # along its bottom edge rather than standing in a row above it.
    order = ['mol_viewer_stack', 'isomer_nav_row', 'xyz_copy_row',
             'submit_ff_notes']
    positions = [children.index(name) for name in order]
    assert positions == sorted(positions), children

    enable = source.split('def _enable_live_forcefield')[1].split('\n    def ')[0]
    assert "_set_ff_notes(payload.get('warnings') or [])" in enable
    assert '_set_mol_status(*warnings' not in enable
    # Escaped, because these carry element symbols and free text.
    notes = source.split('def _set_ff_notes')[1].split('\n    def ')[0]
    assert 'html.escape' in notes
    # Switching the field off clears them.
    assert '_set_ff_notes([])' in source


def test_a_held_value_can_be_a_pull_or_an_exact_fix():
    """A pull negotiates with the chemistry and settles at a compromise; a fix
    is restored after every relaxation step, so the value is met exactly and
    the rest of the molecule arranges itself around it.

    Measured on cholesterol, the same angle asked for 180 degrees: pull settles
    at 137.4 with the field at 206 kcal/mol, fix reaches 180.0 at 2091 -- the
    strain an exact value costs, which the energy readout shows."""
    assert 'function setFixedInternals' in EDITOR
    assert 'function applyFixedInternals' in EDITOR
    # The exact value is restored after the field has had its say, not before.
    relax = _body('ffRelaxFrame')
    assert relax.index('ffWritePositions(viewer, out)') < relax.index(
        'applyFixedInternals(scopeKey)'
    )
    # Set and a fix share one geometry routine, which does no bookkeeping.
    core = _body('applyInternalValue')
    assert 'snapshotForUndo' not in core
    assert 'pushXyzToPython' not in core

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'submit_hold_mode' in source
    # Only pulls become force-field terms; fixes are re-imposed in the browser.
    assert "if c.get('mode', 'pull') == 'pull'" in source
    assert "if c.get('mode') == 'fix'" in source
    # A mode can be changed without setting the constraint again.
    mode = source.split('def on_submit_hold_mode')[1].split('\n    def ')[0]
    assert 'mode=submit_hold_mode.value' in mode


def test_holding_a_value_is_a_constraint_and_not_a_fight():
    """A held value used to be enforced by *setting* the coordinate after each
    step and then superimposing the molecule back to cancel the drift that
    caused -- an overdamped deform-and-restore cycle swims, 35.5 degrees of
    rotation and 1.6 A of drift in four seconds on a real Re complex.

    The superposition cured the swimming but not the cause: the field pulled
    one way and the enforcement snapped back every frame, so the structure
    settled into a cycle rather than a minimum. On cholesterol with one angle
    held at 95 degrees, the rest was still moving 17 to 33 milliangstrom per
    frame after three seconds.

    The correction is built from the coordinate's own gradient now. That
    gradient is orthogonal to every rigid translation and rotation -- moving
    along it changes the value and nothing else -- so no rigid-body motion is
    injected and the superposition is gone with it. Sweeping until each
    constraint is satisfied is SHAKE, which is also what stops several of them
    undoing one another. Measured again: 0.11 to 1.38 milliangstrom per frame,
    the same as the molecule relaxing freely, with the angle still at 95.00."""
    body = _body('applyFixedInternals')
    assert 'superimposeOnto' not in body
    assert 'entry.value - current' in body            # the error, not the value
    assert 'lambda * grad[3 * d]' in body             # moved along the gradient
    assert 'CONSTRAINT_SWEEPS' in body                # and swept until met
    # A dihedral is periodic; the short way round.
    assert 'while (error > 180) error -= 360;' in body
def test_two_ligands_are_exchanged_rather_than_dragged_past_each_other():
    """Two arrangements of the same ligand set are separate minima, and a
    steepest-descent relaxation only runs downhill: it can never cross between
    them. A ligand dragged part of the way rolls back into the basin it came
    from -- for an octahedron the saddle is a Bailar or Ray-Dutt twist, a long
    way from either end. So the exchange is performed in one step.

    An animated swap would be worse than a jump, not better: both ligands
    travelling the same arc in opposite senses meet in the middle, whereas a
    jump has no middle. The landing is the real risk, and each ligand is spun
    about its own new metal-donor axis to find the roomiest one.

    Measured on a real Re complex, exchanging the nitrosyl with a bromide: the
    N lands 0.35 A from where the Br sat and the Br 0.35 A from where the N
    sat, closest contact 2.28 A, and after four seconds of settling the N is
    still 2.79 A from its own old place -- it stayed swapped."""
    swap = _body('exchangeLigands')
    # Rotation about the metal, so each ligand keeps its own bond length.
    assert 'rotateAtomsAbout(atoms, fragA.atoms, centre' in swap
    assert 'rotateAtomsAbout(atoms, fragB.atoms, centre' in swap
    # Two arms of one chelate cannot trade places.
    assert 'same chelate' in swap
    # A trans pair has no unique rotation axis; any perpendicular one serves.
    turn = _body('rotationBetween')
    assert 'Math.PI' in turn and 'seed' in turn
    # The landing is chosen for clearance and reported.
    assert 'spinForClearance' in swap
    assert 'contact: contact' in swap
    # Undoable, and the polyhedron works its vertices out afresh afterwards.
    assert 'snapshotForUndo(scopeKey)' in swap
    assert "pushXyzToPython(scopeKey, 'drag-end')" in swap

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    handler = source.split('def on_submit_swap')[1].split('\n    def ')[0]
    assert 'exchangeLigands(' in handler
    # Offered without a polyhedron too: an exchange is useful either way.
    refresh = source.split('def _refresh_swap')[1].split('\n    def ')[0]
    assert 'metal_indices' in refresh
    assert "state.get('poly_applied')" not in refresh


def test_bonds_can_be_drawn_and_removed_by_hand():
    """Bond perception reads distances, and in a crowded coordination sphere
    that is not reliable. On a real Pt complex the two ipso carbons of a
    phosphine's phenyls sit 2.32 and 2.53 A from the metal, in among the true
    donors at 2.25 to 2.40, and were counted as donors: CN 6 for a
    four-coordinate square-planar Pt. The viewer's own perception was no
    better, inventing a Pt-H bond at 1.88 A instead. Neither can be trusted,
    so the user has to be able to say what is bonded.

    Removing the two spurious Pt-C bonds takes the complex from CN 6 offering
    octahedral and trigonal prismatic to CN 4 offering square planar, which is
    what a Pt(II) complex is."""
    edit = _body('editBond')
    # It edits the model's own bond list, which is what the sticks are drawn
    # from and what every geometry operation here already reads. The two halves
    # of a bond are linked and unlinked by their own functions, because putting
    # a correction back after a rebuild does exactly the same thing.
    assert 'connectOne(i, j); connectOne(j, i);' in edit
    assert 'disconnectOne(i, j); disconnectOne(j, i);' in edit
    assert 'atoms[a].bonds.push(b)' in _body('linkOne')
    assert 'list.splice(at, 1)' in _body('unlinkOne')
    assert 'atoms[a].bondOrder' in _body('linkOne')
    assert 'snapshotForUndo(scopeKey)' in edit
    assert "pushXyzToPython(scopeKey, 'drag-end')" in edit

    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    assert 'submit_bond_btn' in source and 'submit_unbond_btn' in source
    # The correction is remembered and laid over perception, or the next
    # perception -- which runs from the geometry -- would quietly undo it.
    apply_edits = source.split('def _apply_bond_edits')[1].split('\n    def ')[0]
    assert 'apply_bond_edits(perceived, state.get(\'bond_edits\') or {})' in apply_edits
    assert '_apply_bond_edits(perceived)' in source
    # And it is dropped when a different structure arrives.
    view = source.split('def update_molecule_view')[1].split('\n    def ')[0]
    assert "state['bond_edits'] = {}" in view

    pytest.importorskip('rdkit')
    from delfin.dashboard.molecule_forcefield import apply_bond_edits, perceive_molecule

    perceived = perceive_molecule(_BENZENE)
    assert (0, 3) not in perceived.bonds and (0, 1) in perceived.bonds
    assert apply_bond_edits(perceived, {(0, 3): True, (0, 1): False}) is True
    assert (0, 3) in perceived.bonds and (0, 1) not in perceived.bonds
    # Out-of-range pairs are ignored rather than raising.
    assert apply_bond_edits(perceived, {(0, 99): True}) is False


def test_a_drawn_bond_is_parameterised_as_a_bond_not_as_a_distance():
    """A hand-drawn bond used to keep the length it was drawn at.

    Only the bond *list* was corrected; the molecule the UFF parameters are
    read from still had the original connectivity, so RDKit had no entry for
    the new bond and the exporter fell back to its geometric estimate --
    equilibrium = the distance it happened to be drawn at. Joining two carbons
    across a benzene ring gave r0 = 2.798 A with k = 111 instead of 1.514 A
    with k = 700, so the bond never contracted and neither carbon changed
    hybridisation.

    The typing molecule is rebuilt from the corrected bonds now, so RDKit
    re-perceives the chemistry: both carbons go sp3, and the two neighbours
    they took a double bond from pair up again -- Dewar benzene, whose
    remaining C=C are real double bonds rather than adjacent radicals."""
    pytest.importorskip('rdkit')
    from delfin.dashboard.molecule_forcefield import (
        apply_bond_edits, export_forcefield_terms, perceive_molecule,
    )

    perceived = perceive_molecule(_BENZENE)
    apply_bond_edits(perceived, {(0, 3): True})
    payload = export_forcefield_terms(_BENZENE, perceived=perceived, method='uff')
    terms = {(b['i'], b['j']): b for b in payload['bonds']}

    drawn = terms[(0, 3)]
    assert 1.45 < drawn['r0'] < 1.58, drawn      # an sp3 C-C bond, not 2.798
    assert drawn['k'] > 600, drawn               # and a real force constant
    # The carbons it joined are 2.795 A apart in the ring, so the field now
    # pulls them together instead of holding them where they were drawn.
    assert drawn['r0'] < 2.0

    # What the edit did not touch keeps its own chemistry: C1=C2 is the double
    # bond of Dewar benzene, well short of the 1.514 A an all-single fallback
    # would have given every bond in the molecule.
    assert terms[(1, 2)]['r0'] < 1.40, terms[(1, 2)]


def test_a_bond_edit_actually_reaches_the_force_field():
    """Cutting a bond changed the picture and the coordination number but not
    what the relaxation was doing: the perception is cached by element
    sequence, which a bond edit does not alter, so the cached one came back
    unchanged and no re-export was triggered at all. The field went on holding
    the bond that had just been cut.

    On the Pt complex, removing the two spurious Pt-C bonds takes the terms at
    the metal from 6 bonds and 15 angles to 4 and 6."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    edit = source.split('def _edit_bond')[1].split('\n    def ')[0]
    # The cache has to be dropped, or the correction never leaves the picture.
    assert "state['perceived'] = None" in edit
    assert "state['perceived_for'] = None" in edit
    # And the parameters are rebuilt at once.
    assert '_enable_live_forcefield()' in edit
    assert edit.index("state['perceived'] = None") < edit.index('_enable_live_forcefield()')


def test_the_selection_is_released_once_a_value_has_been_set():
    """Picks accumulate, so leaving them standing after Hold or a bond edit
    meant the next atom clicked joined them: three atoms became four and the
    next constraint was built from the wrong set. That is why several values
    could not be held at once. The highlight spheres stayed on the molecule
    too, reading as though something were still selected.

    Set is the exception, and deliberately so: it is a mode for turning a
    value by hand, and letting go of the picks after every tenth of a degree
    is what made turning something by hand impossible.

    Verified in a browser: after clearing, picks and highlight spheres both go
    to zero while the pivot and any held atoms stay."""
    body = _body('clearSelection')
    assert 'state.picks = [];' in body
    # It asks for a drawing -- marks only, because dropping a selection
    # moves nothing.
    assert 'redrawHighlights(scopeKey, true)' in body
    # Narrower than Clear: the pivot and the held atoms are not touched.
    assert 'unpinAll' not in body
    assert 'state.pivot' not in body

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'def _clear_selection' in source
    # Choosing a polyhedron is the same kind of act: the metal has done its
    # job, and leaving it picked meant the next atom clicked joined it.
    for handler in ('on_submit_hold', '_edit_bond', 'on_submit_poly_changed'):
        body = source.split(f'def {handler}')[1].split('\n    def ')[0]
        assert '_clear_selection()' in body, handler
    # ... and Set, which is a mode now, must NOT let go of them.
    body = source.split('def on_submit_set_internal')[1].split('\n    def ')[0]
    assert '_clear_selection()' not in body


def test_undo_answers_for_operations_not_for_relaxation_frames():
    """The dynamic optimiser used to push an undo snapshot every two seconds.

    The stack holds 50, so a minute and forty seconds of it evicted every real
    operation: Undo then stepped back one relaxation frame at a time instead of
    taking back the angle that had just been set. Measured in a browser --
    9 s of Dynamik Opt after setting an angle left the stack 7 deep, and one
    Undo landed 0.0022 A from where the optimiser had it, which is nothing at
    all. Without the periodic snapshot the stack is 2 deep and the same Undo
    lands exactly on the geometry the Set produced, 0.093 A back.

    Switching the field on still snapshots: that is an operation too, and it is
    what lets Undo return to the geometry from before the optimiser ran."""
    body = _body('autoOptimizeTick')
    assert 'snapshotForUndo' not in body
    assert 'pushXyzToPython' in body          # the coordinate box still follows
    assert 'snapshotForUndo(scopeKey);' in _body('startAutoOptimize')
    # The operations that do answer to Undo.
    for name in ('setInternal', 'editBond', 'exchangeLigands'):
        assert 'snapshotForUndo(scopeKey);' in _body(name), name


def test_the_hybridisation_of_a_picked_atom_can_be_overruled():
    """Perception reads bond orders off the geometry, and carbons that should
    be sp2 come back sp3 -- so their angles are typed at 109.5 degrees and the
    centre puckers instead of staying trigonal planar.

    The offer follows a single picked atom, like the polyhedron does, and is
    withheld for metals: RDKit's UFF has no types for one at all, so its bonds
    and angles come from the geometry either way."""
    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    assert 'submit_hyb_dd' in source
    assert 'submit_hyb_dd.observe(on_submit_hyb_changed' in source

    offer = source.split('def _refresh_hybridisation')[1].split('\n    def ')[0]
    assert 'perceived.metal_indices' in offer          # not for a metal
    assert 'perceived_hybridisation_of' in offer       # names what automatic means

    handler = source.split('def on_submit_hyb_changed')[1].split('\n    def ')[0]
    # The cache is keyed by element sequence, which this does not change.
    assert "state['perceived'] = None" in handler
    assert '_enable_live_forcefield()' in handler
    assert '_clear_selection()' in handler

    # Applied after the bond edits, because rebuilding the typing molecule
    # sanitizes it and sanitisation re-perceives hybridisation.
    perception = source.split('def _perception_for')[1].split('\n    def ')[0]
    assert perception.index('_apply_bond_edits') < perception.index('_apply_hyb_overrides')
    # And dropped when a different structure arrives.
    view = source.split('def update_molecule_view')[1].split('\n    def ')[0]
    assert "state['hyb_overrides'] = {}" in view


def test_a_bond_can_be_clicked_instead_of_both_of_its_atoms():
    """Removing a bond meant picking its two atoms first, which is two clicks
    for something drawn as one object.

    Clicking a stick selects the two atoms it joins, so everything that reads
    the selection -- Unbond, the value box, Set, Hold -- keeps working
    unchanged. Only the middle of the stick counts: verified in a browser on
    cholesterol by walking along one bond, the hit is the atom at 0.00 and
    0.15, the bond at 0.35, 0.50 and 0.65, and the far atom at 0.85 and 1.00.
    Clicking the same stick twice clears the selection again."""
    body = _body('raycastBond')
    # Point-to-segment in the same projection the atom picker uses.
    assert 'projectWithDepth' in body
    # Only the middle of the stick belongs to the bond; the ends belong to
    # the atoms. Without that, a tap aimed at an atom whose drawn disc is
    # small was answered with the bond, so a three-atom selection for an angle
    # silently became two and Hold had nothing sensible to hold.
    assert 'if (t < BOND_PICK_MARGIN || t > 1.0 - BOND_PICK_MARGIN) continue;' in body
    assert 'var BOND_PICK_MARGIN = 0.3;' in EDITOR
    assert 'if (b <= a) continue;' in body          # every bond exactly once
    assert 'depth < bestDepth' in body              # the near stick shields the far one

    pick = _body('pickBond')
    assert 'state.picks =' in pick
    assert 'redrawHighlights(scopeKey)' in pick

    click = _body('installCanvasClickFallback')
    # An atom squarely under the cursor wins; only then does the stick get its
    # turn; 3Dmol's own picker and the slack pass both answer with an atom for
    # a click in the middle of a bond, so both have to wait.
    assert 'raycastAtom(scopeKey, cx, cy, true)' in click
    assert click.index('raycastAtom(scopeKey, cx, cy, true)') < click.index('raycastBond')
    assert click.index('raycastBond') < click.index('getAtomBySerial')
    assert 'if (exactOnly) return null;' in _body('raycastAtom')


def test_delete_removes_the_selected_bond_and_nothing_else():
    """Unbond is not a picture edit -- it changes the topology the force field
    is built from -- so the browser cannot carry it out alone and asks through
    a hidden widget.

    Verified in a browser, all in select mode: one atom picked sends
    ``delatoms:1:0``, three send ``delatoms:2:0,5,9``, two bonded atoms picked
    one at a time send ``delatoms:3:0,1`` -- and the stick between those same
    two, tapped, sends ``unbond:4:0,1``. Delete while the caret is in a text
    field sends nothing at all."""
    editor = EDITOR
    assert "key === 'Delete' || key === 'Backspace'" in editor
    guard = editor.split("key === 'Delete' || key === 'Backspace'")[1][:2200]
    assert 'typingInAField()' in guard
    assert 'bondsOf(' in guard                      # only a bond that is there
    assert "pushCommandToPython(names[s], 'unbond'" in guard
    # What the key means follows from how the selection was made, not from the
    # mode: a stick that was tapped comes off as a bond, atoms picked one at a
    # time are deleted. So it reads the same wherever the user is.
    assert 'scope.pickedAsBond' in guard
    assert "pushCommandToPython(names[s], 'delatoms'" in guard
    assert 'state.pickedAsBond = false;' in _body('togglePick')
    assert 'state.pickedAsBond = !isExactly;' in _body('pickBond')

    push = _body('pushCommandToPython')
    assert '.submit-cmd-sync' in push
    # The counter is what makes the value change: the same command twice in a
    # row is two commands, and a widget only reports a change.
    assert 'commandSerial += 1;' in push

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert "submit_cmd_sync.add_class('submit-cmd-sync')" in source
    assert 'submit_cmd_sync.observe(on_submit_cmd' in source
    handler = source.split('def on_submit_cmd')[1].split('\n    def ')[0]
    assert "if verb == 'unbond':" in handler
    assert '_edit_bond(False)' in handler


def test_a_hybridisation_can_be_forced_on_a_whole_selection():
    """One unperceived double bond usually costs both of its carbons their
    type, and retyping a ring one atom at a time is busywork.

    Driven through the real tab: picking C0, C1 and C2 of benzene and choosing
    sp3 sends a force field whose angles at all three are 109.5 degrees while
    C4, which was not picked, keeps 120.0 -- and choosing the automatic entry
    again puts C0 back to 120.0."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    offer = source.split('def _refresh_hybridisation')[1].split('\n    def ')[0]
    # Every picked atom, metals dropped rather than blocking the offer -- and
    # an index the structure no longer has dropped too, because the browser
    # pushes its picks after a re-render and an edit renumbers the atoms.
    assert 'i not in metals' in offer
    assert 'len(perceived.symbols)' in offer
    assert "state['hyb_offer_atoms'] = chosen" in offer
    # A value can only be shown as current when they all share it.
    assert 'held.pop() if len(held) == 1' in offer

    handler = source.split('def on_submit_hyb_changed')[1].split('\n    def ')[0]
    assert 'for atom in atoms:' in handler

    # The perception has to be available for any pick count, not just one --
    # it used to be fetched only inside the single-atom polyhedron branch.
    sync = source.split('def on_submit_pick_sync')[1].split('\n    def ')[0]
    assert 'if xyz and indices:' in sync
    assert sync.index('_perception_for(xyz)') < sync.index('polyhedron_options(')


def test_a_second_rubber_band_adds_to_the_selection():
    """Drawing a new box threw away everything the last one had picked, so a
    molecule could only ever be selected in one rectangle.

    The band already needs Shift, and a plain click on an atom has always
    accumulated -- the band was the one gesture that cleared. It used to add
    only when Ctrl was held on top of Shift, which is two modifiers for what
    the gesture already means."""
    editor = EDITOR
    band = editor.split("kind: 'maybe-rect',")[1][:400]
    assert 'additive: true,' in band
    assert 'ctrlKey' not in band


def test_carbon_types_can_be_read_off_the_connectivity():
    """Perception goes through bond orders, and those are read from the
    geometry -- a double bond twisted out of plane, or at an unusual length,
    is simply not seen and its carbon comes back sp3. The number of partners
    says it outright, because carbon has no lone pair."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'submit_hyb_auto_btn' in source
    assert 'submit_hyb_auto_btn.on_click(on_submit_hyb_auto)' in source
    handler = source.split('def on_submit_hyb_auto')[1].split('\n    def ')[0]
    assert 'hybridisation_from_connectivity' in handler
    # The selection when there is one, the whole structure when there is not.
    assert 'picked or None' in handler
    assert "state['perceived'] = None" in handler
    assert '_enable_live_forcefield()' in handler


def test_a_polyhedron_held_on_one_metal_is_not_offered_for_the_next():
    """The dropdown is rebuilt from whichever metal is picked, and its value
    was then set to whatever polyhedron happened to be applied -- even when
    that metal does not offer it. Assigning a value the options do not contain
    raises, so picking a four-coordinate metal while an octahedron was held on
    another one took the whole handler down.

    Only reachable once the selection is released after applying a polyhedron:
    before that the second metal made a two-atom pick, which offers nothing."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    body = source.split('def on_submit_pick_sync')[1].split('\n    def ')[0]
    assert 'offered = {code for code, _label in choices}' in body
    assert 'submit_poly_dd.value = applied if applied in offered else \'\'' in body


def test_turn_is_offered_only_where_the_vertices_differ():
    """Which ligands take the axial positions of a trigonal bipyramid is a
    real choice; an octahedron has nothing to turn, because every vertex is
    the same and exchanging two ligands there is what Swap is for.

    Driven through the real tab: choosing a trigonal bipyramid on a
    five-coordinate iron shows Turn and ten clicks walk through ten distinct
    axial pairs before coming back; choosing an octahedron or a trigonal
    prism on the Re complex leaves it hidden."""
    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    assert 'submit_poly_turn_btn' in source
    assert 'submit_poly_turn_btn.on_click(on_submit_poly_turn)' in source

    offer = source.split('def _refresh_poly_turn')[1].split('\n    def ')[0]
    assert 'polyhedron_vertex_classes' in offer
    assert 'len(set(grouped[0])) > 1' in offer

    handler = source.split('def on_submit_poly_turn')[1].split('\n    def ')[0]
    assert 'polyhedron_arrangements' in handler
    # The coordinates as they are now: a ligand that has been dragged has to
    # be scored where it sits, not where it was perceived.
    assert 'parse_xyz' in handler
    assert "% len(arrangements)" in handler
    assert '_enable_live_forcefield()' in handler
    # Choosing a different polyhedron starts the cycle over.
    changed = source.split('def on_submit_poly_changed')[1].split('\n    def ')[0]
    assert "state['poly_arrangement_index'] = 0" in changed
    assert '_refresh_poly_turn()' in changed


def test_draw_mode_hands_every_gesture_to_python():
    """What the browser contributes is where the user pointed. How many
    hydrogens an atom needs and where they go is decided in Python, where
    RDKit's valences and covalent radii are -- so each gesture ends as one
    command and the structure comes back rendered.

    Driven through the real tab from methane: growing a carbon gives ethane,
    retyping a hydrogen to nitrogen gives C2NH7, placing an oxygen adds a
    water, raising the C-C bond to double drops two hydrogens, and deleting
    the oxygen takes its own two with it. A nonsense index or element changes
    nothing."""
    assert "mode === 'draw'" in EDITOR
    finish = _body('finishDraw')
    # Tap an atom: retype it. Drag onto another: bond them. Drag into space:
    # grow. Nothing under the press at all: place one where the cursor is.
    for verb in ('setelement', 'bondorder', 'grow', 'addatom'):
        assert f"'{verb}'" in finish, verb
    assert 'drag.movedEnough' in finish

    # Depth cannot be read back from a click, so it is borrowed from something
    # already in the scene -- otherwise a placed atom lands behind the molecule.
    world = _body('screenToWorld')
    assert 'getCameraBasis' in world and 'getPixelToWorld' in world
    assert 'projectWithDepth' in world

    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    assert 'submit_draw_btn' in source
    # No order to choose in advance: a drawn bond is single, and what it
    # should be is decided afterwards by tapping the stick, where it can be
    # seen. Having to set it beforehand was a control that was nearly always
    # on the wrong value.
    assert 'submit_element_dd' in source
    assert 'submit_order_dd' not in source
    handler = source.split('def on_submit_cmd')[1].split('\n    def ')[0]
    for verb in ('addatom', 'grow', 'setelement', 'bondorder', 'delatoms'):
        assert f"'{verb}'" in handler, verb
    # The edited topology is re-seeded, because perception reads bond orders
    # off the geometry and does not get them back: ethene built by hand came
    # back as a single bond at 1.514 A.
    apply_ = source.split('def _apply_structure')[1].split('\n    def ')[0]
    assert "state['bond_edits'] = {" in apply_
    assert 'coords_widget.value' in apply_


def test_undo_reaches_structural_edits_too():
    """A snapshot of coordinates cannot bring back an atom that was deleted or
    take away one that was placed, so structural edits keep their own stack on
    the Python side.

    The order stays right without either side keeping a clock: every
    structural edit re-renders, and a re-render clears the browser's stack --
    so an empty stack there means the next thing to undo is a structural edit.
    Verified through the real tab: two grows and then Undo walks back
    C3NH8O to C3NH7 to C2NH5, and going on undoes the whole session's edits
    back to the methane it started from. In a browser, Undo with an empty
    stack sends ``undo:6:structure`` rather than doing nothing."""
    body = _body('undo')
    assert 'if (!state.undo.length) {' in body
    assert "pushCommandToPython(scopeKey, 'undo', 'structure')" in body

    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    handler = source.split('def on_submit_cmd')[1].split('\n    def ')[0]
    assert "if verb == 'undo':" in handler

    step = source.split('def _undo_structure')[1].split('\n    def ')[0]
    assert "state['history']" in step
    # The putting-back is its own function now, because Reset does it too.
    assert '_restore(entry' in step
    restore = source.split('def _restore')[1].split('\n    def ')[0]
    assert "coords_widget.value = entry['coords']" in restore

    apply_ = source.split('def _apply_structure')[1].split('\n    def ')[0]
    # The previous structure is remembered before the new one is written, and
    # the write must not be mistaken for a different molecule arriving.
    assert '_remember(' in apply_, 'the edit is not in the history'
    assert "state['structure_edit_inflight'] = True" in apply_
    view = source.split('def update_molecule_view')[1].split('\n    def ')[0]
    assert "if not state.get('structure_edit_inflight'):" in view


def test_drawing_does_not_reset_the_camera_or_stop_the_field():
    """Two things used to happen on every atom placed: the view snapped back
    to the default, and the continuous relaxation stopped.

    Both follow from the re-render. The camera was kept only when the atom
    count matched, which a structural edit by definition breaks -- but an edit
    is not a different molecule, it is the one being worked on. And
    onViewerReady clears the running loop, so drawing while Dynamik Opt ran
    silently switched it off; the parameters are re-assigned for the structure
    that includes the new atom, which is the moment to pick it up again."""
    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    assert 'var edited = !!window.__delfinStructureEdit;' in source
    assert '(edited || saved.atoms === count)' in source
    mark = source.split('def _mark_structure_edit')[1].split('\n    def ')[0]
    assert '__delfinStructureEdit = true' in mark
    # The resume flag is set with the parameters, not before the re-render:
    # setting it earlier meant it could be consumed against the viewer that
    # was going away, and the relaxation came back stuck until the toggle was
    # cycled by hand.
    enable = source.split('def _enable_live_forcefield')[1].split('\n    def ')[0]
    assert '__delfinResumeAutoOpt = {resume}' in enable
    assert 'submit_relax_btn.value' in enable
    assert enable.index('__delfinResumeAutoOpt') < enable.index('setForceField')

    field = _body('setForceField')
    assert 'window.__delfinResumeAutoOpt' in field
    assert 'startAutoOptimize(scopeKey)' in field


def test_tapping_a_bond_in_draw_mode_retypes_it():
    """A tap steps it on: single, double, triple, single. The hydrogens and
    the length follow, and the shape with them -- ethane becomes ethene with
    its carbons trigonal planar, then ethyne, linear."""
    finish = _body('finishDraw')
    assert 'drag.bond' in finish
    assert "'bondcycle'" in finish
    assert 'drag.bond[0]' in finish and 'drag.bond[1]' in finish
    # Only a tap, never the end of a drag -- that gesture already means grow.
    assert '!drag.movedEnough' in finish
    assert 'raycastBond(scopeKey, e.clientX, e.clientY)' in EDITOR


def test_a_double_bond_is_drawn_as_two_sticks():
    """A model read from an XYZ block has no bond orders in it -- the format
    carries none -- so every bond was drawn as one stick whatever it was.

    3Dmol draws a double as two cylinders and a triple as three once the model
    knows, so the orders are handed over after every render that changes them.
    Measured in a browser on ethene: the C-C order goes 1, 2, 3 and the
    geometry built for it goes 36, 40, 44 vertices."""
    body = _body('setBondOrders')
    assert 'atoms[i].bondOrder[at] = order;' in body
    assert 'atoms[j].bondOrder[back] = order;' in body   # both ends
    assert 'invalidateGeometry(viewer)' in body
    assert 'if (order < 1 || order > 3) continue;' in body

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    push = source.split('def _push_bond_orders')[1].split('\n    def ')[0]
    assert 'setBondOrders(' in push
    # Only what differs from a plain stick needs sending.
    assert 'if int(order) > 1' in push
    apply_ = source.split('def _apply_structure')[1].split('\n    def ')[0]
    assert '_push_bond_orders(structure.bonds)' in apply_


def test_a_placed_atom_lands_under_the_cursor():
    """It landed a quarter of the way further out than the click, and further
    the further from the centre it was -- a clean scale error of 1.254.

    getPixelToWorld derives its scale from the field of view and the camera
    distance, and the model group carries a scale of its own on top of that.
    So the transform is measured instead of worked out: projecting one world
    unit along each screen axis says what it actually is, whatever else is in
    the chain. Measured in a browser at five points across the canvas, each
    now lands 0.0 px from where it was clicked."""
    world = _body('screenToWorld')
    assert 'getPixelToWorld(' not in world      # named only in the comment
    assert 'var alongRight = probe(basis.right);' in world
    assert 'var alongUp = probe(basis.up);' in world
    # A 2x2 solve, guarded against a degenerate projection.
    assert 'Math.abs(det) < 1e-9' in world


def test_the_editor_version_is_its_own_content():
    """It was a number to bump by hand, and it was not bumped once across a
    day of changes -- so an open dashboard kept running the editor it had
    loaded and every fix shipped in it was invisible, which is the very thing
    the version was added to prevent."""
    import re

    from delfin.dashboard.molecule_viewer import SUBMIT_MANIP_BOOTSTRAP_JS

    assert "var MANIP_VERSION = '__DELFIN_MANIP_VERSION__';" in SUBMIT_MANIP_BOOTSTRAP_JS
    stamped = submit_manip_bootstrap_js()
    assert '__DELFIN_MANIP_VERSION__' not in stamped
    match = re.search(r"var MANIP_VERSION = '([0-9a-f]{12})';", stamped)
    assert match, 'the version is not stamped in'
    # And it moves when the script does -- including the program the worker
    # runs, which is part of the editor and goes into the same hash.
    import hashlib
    import json

    from delfin.dashboard.molecule_viewer import FF_WORKER_LOOP_JS

    full = SUBMIT_MANIP_BOOTSTRAP_JS.replace(
        '__DELFIN_FF_WORKER_LOOP__', json.dumps(FF_WORKER_LOOP_JS))
    assert match.group(1) == hashlib.sha256(full.encode('utf-8')).hexdigest()[:12]


def test_the_right_button_takes_things_away_in_draw_mode():
    """An atom under it goes, a stick under it loses its bond. On empty space
    it still belongs to the viewer, so the scene can be turned and panned
    without leaving the mode -- which is the only way to see what was built."""
    overlay = EDITOR.split("ov.addEventListener('mousedown'")[1]
    draw = overlay.split("if (state.mode === 'draw') {")[1].split(
        "if (state.mode === 'manipulate')")[0]
    right = draw.split('if (e.button === 2) {')[1]
    assert "pushCommandToPython(scopeKey, 'delatoms'" in right
    assert "pushCommandToPython(scopeKey, 'unbond'" in right
    assert 'raycastBond(scopeKey, e.clientX, e.clientY)' in right
    # Nothing under the cursor: the press is left alone.
    assert 'if (!stick) return;' in right


def test_a_held_value_follows_the_atoms_through_an_edit():
    """Constraints name atoms by index, and drawing renumbers: a deleted
    hydrogen moves every atom after it. A held value that quietly pointed at
    different atoms afterwards, or vanished without a word, is worse than one
    that is dropped and said so.

    The builder records where every atom came from, so everything the tab
    holds by index -- constraints, forced types, the metal a polyhedron sits
    on -- is carried across. Driven through the real tab: an angle held on
    three carbons survives saturating one of them, and is dropped with a note
    when the atom it names is deleted."""
    from delfin.dashboard import tab_submit
    from delfin.dashboard.molecule_builder import Structure

    structure = Structure(['C', 'C', 'C'],
                          [(0, 0, 0), (1.54, 0, 0), (3.08, 0, 0)], {})
    assert structure.renumbering() == {0: 0, 1: 1, 2: 2}
    structure.remove_atoms([0])
    assert structure.renumbering() == {1: 0, 2: 1}
    structure.add_atom('H', (0.0, 0.0, 1.0))
    assert structure.renumbering() == {1: 0, 2: 1}   # a new atom came from nowhere

    source = _TAB_AND_EDITOR
    apply_ = source.split('def _apply_structure')[1].split('\n    def ')[0]
    assert 'renumber = structure.renumbering()' in apply_
    assert "state['constraints'] = kept" in apply_
    assert "state['hyb_overrides'] = {" in apply_
    assert "state['poly_metal'] = renumber[metal]" in apply_
    assert 'lost' in apply_                       # and it says so
    # An edit is not a different molecule, so none of this is thrown away
    # underneath it.
    view = source.split('def update_molecule_view')[1].split('\n    def ')[0]
    assert "if not state.get('structure_edit_inflight'):" in view
    assert view.index("structure_edit_inflight") < view.index("state['constraints'] = []")


# ---------------------------------------------------------------------------
# a bond drawn by hand has to survive the next rebuild of the picture
# ---------------------------------------------------------------------------
def test_hand_drawn_bonds_are_put_back_into_a_rebuilt_picture():
    """A rebuild draws what the distances say.

    A bond drawn between two atoms that are not within bonding distance is not
    in the coordinates, so perception never finds it again: the first drawn
    bond disappeared from the view the moment a second edit rebuilt it, while
    still being in force everywhere else -- which is why an optimisation, which
    pulled the two atoms together, made it visible again.
    """
    body = _body('applyBondEdits')

    # It puts one back and takes one away, both directions of the correction.
    assert 'linkOne(atoms, i, j); linkOne(atoms, j, i);' in body
    assert 'unlinkOne(atoms, i, j); unlinkOne(atoms, j, i);' in body
    # Nothing is touched that already agrees, so repeating it is free.
    assert 'if (connect === linked) continue;' in body
    # A newer call wins over an older one's pending attempts, so a correction
    # the user has since taken back is not re-asserted.
    assert 'window.__delfinBondEditGeneration !== generation' in body
    assert body.count('setTimeout(once,') == 2

    assert 'applyBondEdits: applyBondEdits,' in EDITOR


def test_only_the_corrections_from_a_hand_are_put_back():
    """Never the whole bond list: what perception found is the viewer's own
    business, and after a structural edit the bond list is the whole structure
    and no longer says which bonds came from a hand."""
    from delfin.dashboard import tab_submit

    source = _TAB_AND_EDITOR
    push = source.split('def _push_hand_bonds')[1].split('\n    def ')[0]
    assert "state.get('hand_bonds')" in push
    assert 'bond_edits' not in push

    # Recorded where a bond is drawn or cut, and nowhere else.
    edit = source.split('def _edit_bond')[1].split('\n    def ')[0]
    assert "hand[pair] = bool(connect)" in edit

    # Carried across the renumbering an edit causes, before the coordinates
    # are written -- writing them is what rebuilds the picture.
    apply_ = source.split('def _apply_structure')[1].split('\n    def ')[0]
    assert "state['hand_bonds'] = {" in apply_
    assert apply_.index("state['hand_bonds'] = {") < apply_.index('coords_widget.value =')

    # And put back by the rebuild itself, so every path that rebuilds is covered.
    replace = source.split('def _replace_mol_output_view')[1].split('\n    def ')[0]
    assert '_push_hand_bonds()' in replace


def test_a_different_molecule_keeps_none_of_the_corrections():
    """They name atoms by index, which says nothing about another structure."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    view = source.split('def update_molecule_view')[1].split('\n    def ')[0]
    assert "state['hand_bonds'] = {}" in view
    assert view.index('structure_edit_inflight') < view.index("state['hand_bonds'] = {}")


def test_the_fullscreen_status_is_cleared_with_the_small_one():
    """Both copies say the same thing, or the overlay lies.

    Clearing only the small one left fullscreen showing "Quick convert (single
    structure)..." long after the structure was on screen -- and in fullscreen
    that stale line is the only thing there is to go by.

    Driven rather than read, because clearing goes through the one renderer
    now: it has to leave the ring alone while something is still running, and
    only the renderer knows whether anything is.
    """
    part, _state, _pump = _a_part()
    part._set_mol_status('Quick convert (single structure)...')
    assert part.mol_status.value
    assert part.mol_status_fs.value

    part._clear_mol_status()
    assert part.mol_status.value == ''
    assert part.mol_status_fs.value == ''


# ---------------------------------------------------------------------------
# the force field beside the page rather than in front of it
# ---------------------------------------------------------------------------
def test_the_relaxation_does_not_own_the_frame_the_page_needs():
    """A batch aims at a whole frame, and a whole frame spent on the physics
    is a whole frame the page does not have.

    Measured in a browser with the shipped engine, one second of Dynamik Opt on
    a 100-atom peptide: 411 ms of main-thread JavaScript, and a task that
    wanted to run "now" waited 44 ms at the 95th percentile. Everything else --
    a click, a widget update, a message from the kernel -- was behind that,
    which is what made the whole dashboard drag while the optimisation ran.
    Computed in a Worker instead: 2 ms and 2 ms, with the relaxation itself no
    slower (977 -> 1012 steps/s, same energy after four seconds).
    """
    make = _body('getFFWorker')
    assert 'new Worker(url)' in make
    assert 'window.__delfinFFSource' in make
    # An old browser, or a policy that forbids blob: -- everything then goes on
    # being computed here, exactly as before.
    assert 'ffWorkerRefused' in make
    assert 'done(ffRelaxFrame(scopeKey));' in _body('ffRelaxAsync')
    assert 'function ffRelaxFrame' in EDITOR, 'the fallback has to still exist'
    # A source that has not arrived yet is not a browser that cannot do it:
    # the engine's script and this one arrive separately.
    assert 'if (!source) return null;' in make

    relax = _body('ffRelaxAsync')
    # One batch out at a time: a second would only queue work the picture has
    # already moved past.
    assert 'state.ffBusy' in relax
    # And a batch that never comes back must not hold the loop for ever.
    assert 'ffBusySince' in relax
    # What the *rigid* hand holds stays where the hand put it: the answer
    # describes where it was when it was asked for, which under the cursor is
    # the past.  Under a pull the field owns those atoms too and its answer is
    # written back in full -- that is the whole point of pulling.
    assert ('var held = (!state.ffPull && state.drag && state.drag.targets)'
            in relax)
    assert 'skipSerials' in _body('ffWritePositions')


def test_both_copies_of_the_field_are_told_the_same_thing():
    """The page keeps an engine for what cannot wait -- whether a field is
    loaded at all, the energy of a geometry nobody has relaxed yet -- and the
    worker does the work. Anything that is not a batch reaches both, or they
    would relax different molecules."""
    for command in ('load', 'configure', 'grab', 'dispose'):
        assert f"cmd: '{command}'" in EDITOR, f'the worker never hears about {command}'
    # Frozen atoms above all: a worker that has not been told would relax the
    # atom the user is holding.
    frozen = _body('ffApplyFrozen')
    assert 'window.__delfinFF.grab(scopeKey, list)' in frozen
    assert "ffTellWorker({cmd: 'grab', scope: scopeKey, list: list})" in frozen
    # Statistics come back with every batch, so the energy badge and the
    # settle's convergence test read one place.
    assert 'function ffStatsOf' in EDITOR
    assert 'state.ffStats' in _body('ffRelaxAsync')


def test_the_worker_answers_with_the_positions_and_the_statistics():
    """Both are read every frame; a second round trip for the second would
    cost more than computing either."""
    from delfin.dashboard.molecule_viewer import FF_WORKER_LOOP_JS

    assert 'var window = self;' in EDITOR, 'the engine speaks of window'
    assert 'FF.step(message.scope, message.positions || null' in FF_WORKER_LOOP_JS
    assert 'stats: FF.stats(message.scope)' in FF_WORKER_LOOP_JS
    # The buffer stops being ours the moment it is transferred, and the engine
    # goes on using its own array.
    assert 'new Float64Array(out)' in FF_WORKER_LOOP_JS
    assert '[answer.buffer]' in FF_WORKER_LOOP_JS


def test_the_editor_is_not_sent_again_to_a_page_that_has_it():
    """It is 136 KiB of the 159 a rendered structure weighs.

    Every conversion, every edit and every optimisation result carried a whole
    copy of the editor to a page that already had one and threw it away on its
    version check -- because nothing could tell the kernel otherwise. Now the
    page says which editor it is running, and the copy stops.

    Sending it the first time is the belt: the separate script that installs it
    goes through an output widget whose content can be replaced before the page
    has run it, and a structure without an editor cannot be edited at all.
    """
    from delfin.dashboard import tab_submit
    from delfin.dashboard.molecule_viewer import submit_manip_version

    ready = _body('onViewerReady')
    assert "pushCommandToPython(scopeKey, 'editor', MANIP_VERSION)" in ready

    source = _TAB_AND_EDITOR
    bundle = source.split('def _build_mol_output_bundle')[1].split('\n    def ')[0]
    assert "carry_editor = (state.get('manip_seen_version')" in bundle
    assert "submit_manip_bootstrap_js() if carry_editor else ''" in bundle

    handler = source.split('def on_submit_cmd')[1].split('\n    def ')[0]
    assert "state['manip_seen_version'] = str(payload)" in handler

    # The version is the editor's own content, so a page running an older one
    # is not taken for a page that is up to date.
    assert len(submit_manip_version()) == 12
    assert submit_manip_version() in submit_manip_bootstrap_js()


def test_a_confirmed_page_gets_a_much_smaller_picture(tmp_path):
    """Measured through the tab itself, not asserted from the source."""
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext
    from delfin.dashboard.molecule_viewer import submit_manip_version

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda _script: None
    _widget, refs = tab_submit.create_tab(ctx)
    water = '3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'

    def weight():
        outputs = refs['mol_output'].outputs
        return len(outputs[0]['data']['text/html']) if outputs else 0

    refs['coords_widget'].value = water
    before = weight()

    refs['submit_cmd_sync'].value = f'editor:1:{submit_manip_version()}'
    refs['coords_widget'].value = water.replace('0.96', '0.97')
    after = weight()

    assert after < before / 4, f'{before} -> {after} bytes'

    # An editor from an older build is not one the kernel may rely on.
    refs['submit_cmd_sync'].value = 'editor:2:000000000000'
    refs['coords_widget'].value = water.replace('0.96', '0.98')
    assert weight() > before / 2


def test_the_picture_is_drawn_once_a_frame_however_often_it_is_asked_for():
    """A drag asks twice per mouse event -- once when the cursor moves the
    atom and once when the relaxation answers -- and a mouse reports far more
    often than a screen refreshes.

    Measured over a 60-step drag in a browser: 61 events produced 116 full
    geometry rebuilds and 181 ms of drawing, of which the display could show
    at most a third. Collapsed to one drawing per frame: 65 rebuilds and
    107 ms, for the same drag.
    """
    ask = _body('redrawHighlights')
    assert 'redrawPending' in ask
    assert 'requestAnimationFrame' in ask
    assert 'drawHighlightsNow' in ask
    # The drawing itself is still there, under its own name.
    assert 'function drawHighlightsNow' in EDITOR

    # A picture that is being replaced must not have a frame pending on it.
    ready = _body('onViewerReady')
    assert 'cancelAnimationFrame' in ready
    assert 'state.redrawPending = false;' in ready


def test_marking_atoms_does_not_rebuild_the_whole_model():
    """Picking three atoms moves nothing.

    Only the translucent markers over them change, and rebuilding a
    four-hundred-atom model to draw three spheres costs four times what
    drawing the frame costs -- measured at 405 atoms, 2.3 ms against 0.6, and
    the drawing of a selection fell from 3-5 ms to 1.

    The default is the safe one: a caller that says nothing gets the rebuild,
    so a path that does move atoms cannot go wrong by forgetting to ask.
    """
    ask = _body('redrawHighlights')
    assert 'marksOnly' in ask
    assert 'if (!marksOnly) state.highlightsOnly = false;' in ask

    draw = _body('drawHighlightsNow')
    assert 'if (!state.highlightsOnly || state.dynamicBonds) {' in draw
    assert 'state.highlightsOnly = false;' in draw

    # Only the selection functions claim it.
    for marks_only in ('setPicks', 'clearSelection'):
        assert 'redrawHighlights(scopeKey, true)' in _body(marks_only), marks_only
    for moves in ('applyTranslate', 'setPositions'):
        body = _body(moves)
        assert 'redrawHighlights(scopeKey, true)' not in body, moves


def test_a_drag_does_not_queue_relaxations_it_cannot_use():
    """A mouse reports faster than a batch comes back.

    Every request that found the worker busy hung another retry on the end,
    and that retry hung another: measured at 405 atoms, a thirty-step drag
    asked for 263 animation frames to draw 30, and the pile grew with the
    length of the drag -- which is what made dragging under Dynamik Opt feel
    like wading. One waiting, never a queue: 263 became 58.

    Nothing is lost by dropping the extra requests. The atom is already where
    the cursor put it, and the next batch reads the positions as they are
    then.
    """
    relax = _body('ffRelaxAsync')
    assert 'if (state.ffWaiting) { done(false); return; }' in relax
    assert 'state.ffWaiting = true;' in relax
    assert 'state.ffWaiting = false;' in relax

    # And a rebuilt picture starts with nothing pending on the old one.
    ready = _body('onViewerReady')
    assert 'state.ffWaiting = false;' in ready


def test_no_gutter_is_kept_for_an_output_number_that_will_never_come():
    """A notebook reserves a column on the left for "Out[7]:".

    There is no Out[7] here and never will be, but the column is reserved all
    the same -- and in fullscreen it is a white band down the left of the
    picture, inside the blue frame. The viewer should start where its frame
    starts.
    """
    from delfin.dashboard import tab_submit

    # The ordinary view is this tab's own; the overlay is everyone's, so the
    # gutter is taken away for the Builder and the browser tabs in the same
    # breath -- they show a structure in an Output widget as well.
    for scope, source in (('.submit-mol-output', _TAB_AND_EDITOR),
                          ('.delfin-structure-fs-overlay', FULLSCREEN_CSS)):
        assert f'{scope} .jp-OutputPrompt' in source, scope
        assert f'{scope} .jp-OutputArea-prompt' in source, scope
    # And the child it sat in keeps no padding where it was.
    assert '.delfin-structure-fs-overlay .jp-OutputArea-child' in FULLSCREEN_CSS


def test_the_camera_turns_about_the_system_not_about_where_it_was_loaded():
    """The camera orbits a point, and that point was set once, at load.

    Drag an atom a long way, or let a relaxation carry the molecule, and the
    point is no longer in the molecule -- so turning the view swings the whole
    thing off the screen and into the white. Measured: with the structure 18 A
    off the old centre, none of its atoms were on screen after half a turn;
    re-centred first, all of them were.
    """
    centre = _body('systemCentre')
    # Mass-weighted, as a chemist means it -- not the middle of the box.
    assert 'ATOMIC_MASS[atoms[i].elem]' in centre
    assert 'mx / total' in centre
    assert 'centre.radius' in centre, 'it needs the reach to judge the drift'

    keep = _body('centreOnSystem')
    # The first three of getView are the negated point the camera turns about.
    assert 'view[0] = -centre.x;' in keep
    assert 'viewer.setView(view)' in keep
    # Only when it has drifted: re-centring on every press would undo a
    # deliberate pan.
    assert 'if (!force && drift <= allowed) return false;' in keep
    assert '0.25 * Math.max(centre.radius, 1.0)' in keep

    # And it happens as a turn begins, never during one.
    assert 'centreOnSystem(scopeKey, false);' in EDITOR
    assert 'centreOnSystem: centreOnSystem,' in EDITOR
    turn = EDITOR[EDITOR.index('centreOnSystem(scopeKey, false);'):][:400]
    assert '_handleMouseDown' in turn, 'it must run before the camera is given the press'


def test_a_button_brings_the_system_back_into_view():
    """A view is the user's. Moving in on one corner is something people do on
    purpose, so the picture must not re-frame itself while they work -- but
    when a structure has been dragged out of it, or a relaxation has carried
    it out, there has to be a way back.

    Measured: a molecule turned and pushed 40 A off centre had 0% of its atoms
    on screen at 74 atoms and 25% at 405; after Centre, 100% in both, in a
    millisecond, and not one coordinate changed.
    """
    back = _body('recentreView')
    # Two halves of "I cannot see it": fit what is there, then put the point
    # the camera turns about on the centre of mass.
    assert 'viewer.zoomTo()' in back
    assert 'centreOnSystem(scopeKey, true)' in back
    # It is the camera and only the camera.
    assert 'atoms[' not in back and 'setPositions' not in back
    assert 'recentreView: recentreView,' in EDITOR

    from delfin.dashboard import tab_submit

    source = _EDITOR_PY
    handler = source.split('def on_submit_centre')[1].split('\n    def ')[0]
    assert 'recentreView' in handler
    assert 'structure is unchanged' in handler
    assert 'submit_centre_btn.on_click(on_submit_centre)' in source
    # It lives with the other viewer controls -- which now lie on the picture
    # itself, because that is what they act on: the camera and nothing else.
    panel = source.split('submit_view_body = widgets.VBox')[1].split(')\n')[0]
    assert 'submit_centre_btn' in panel, panel
    assert 'submit_centre_btn.disabled = not enabled' in source


# ---------------------------------------------------------------------------
# one history, and a way back to the beginning
# ---------------------------------------------------------------------------
@pytest.fixture
def tab(tmp_path):
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda _script: None
    _widget, refs = tab_submit.create_tab(ctx)
    return refs


def _atoms(refs):
    return [line for line in refs['coords_widget'].value.splitlines()[2:]
            if line.strip()]


def test_undo_walks_back_through_every_kind_of_step(tab):
    """It could not reach the beginning, and it is plain why.

    There were three histories: a stack in the browser for drags, one here for
    structural edits, and a single slot for the geometry before an
    optimisation. None knew about the others, so Undo walked whichever it
    found first and stopped where that one ran out -- and a re-render clears
    the browser's, which is to say every drawn atom threw away every drag
    before it.
    """
    refs = tab
    water = '3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'
    refs['coords_widget'].value = water

    # a drag, a drawn atom, and the continuous relaxation switched on
    refs['submit_cmd_sync'].value = 'grabbed:1:'
    refs['submit_manip_sync'].value = (
        '3\nDELFIN drag-end\nO 0 0 0\nH 1.40 0 0\nH -0.24 0.93 0\n')
    refs['submit_cmd_sync'].value = 'addatom:2:C,3.0,0.0,0.0'
    refs['submit_relax_btn'].value = True

    named = [entry['what'] for entry in refs['editor_state']['history']]
    assert named[0] == 'the structure as it was loaded'
    assert 'Placed C' in named
    assert any('continuous relaxation' in name for name in named)

    assert len(_atoms(refs)) == 8

    refs['submit_cmd_sync'].value = 'undo:10:'
    assert len(_atoms(refs)) == 8, 'the relaxation switch, not the atom'

    refs['submit_cmd_sync'].value = 'undo:11:'
    assert len(_atoms(refs)) == 3, 'the drawn atom is gone'
    assert '1.40' in _atoms(refs)[1], 'and the drag is still there'

    refs['submit_cmd_sync'].value = 'undo:12:'
    assert '0.96' in _atoms(refs)[1], 'now the drag too'

    # And it stays there rather than falling off the end.
    refs['submit_cmd_sync'].value = 'undo:13:'
    assert '0.96' in _atoms(refs)[1]
    assert 'as it was loaded' in refs['mol_status'].value


def test_reset_goes_to_the_structure_that_was_loaded(tab):
    """It jumped to a second remembered place, set in the render path, and the
    two could disagree -- which is how Reset came back to something that was
    not the beginning."""
    refs = tab
    refs['coords_widget'].value = '3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'
    refs['submit_cmd_sync'].value = 'addatom:1:C,3.0,0.0,0.0'
    refs['submit_cmd_sync'].value = 'addatom:2:N,5.0,0.0,0.0'
    built = len(_atoms(refs))
    assert built > 3

    refs['submit_reset_btn'].click()
    assert len(_atoms(refs)) == 3
    assert refs['editor_state']['history'][0]['what'] == (
        'the structure as it was loaded')
    # And what Reset undid is one more thing that happened, so it is the last
    # entry: somebody who pressed it by accident presses Undo and their work
    # is back.  The history used to be cut to its first entry here, which is
    # the state Reset had just gone to -- so Undo landed on it and said there
    # was nothing more to take back.
    assert refs['editor_state']['history'][-1]['what'] == 'the reset'
    refs['submit_manip_undo_btn'].click()
    assert len(_atoms(refs)) == built, 'Reset threw the structure away for good'


def test_a_step_is_recorded_before_it_happens_not_after(tab):
    """A drag is recorded when the atom is picked up. By the time it is let
    go of, the coordinate box already holds what the drag made -- the
    relaxation pushes into it while the mouse is still down -- and "before the
    drag" would not be anywhere any more."""
    from delfin.dashboard import tab_submit
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor = submit_manip_bootstrap_js()
    assert editor.count("pushCommandToPython(scopeKey, 'grabbed', '')") == 2, (
        'both ways of starting a drag have to say so'
    )

    source = _EDITOR_PY
    handler = source.split('def on_submit_cmd')[1].split('\n    def ')[0]
    assert "if verb == 'grabbed':" in handler
    assert "_remember('a drag')" in handler
    # And an optimisation is remembered before it is started.
    assert "_remember('an optimisation')" in source
    optimise = source.split("_remember('an optimisation')")[1][:200]
    assert 'state[\'optimize_run\'] = token' in optimise


def test_the_hand_pulls_the_atom_instead_of_placing_it():
    """A drag is a force on the atom, not a new coordinate for it.

    Setting the coordinate is why a molecule could be pulled into any shape at
    all: the hand won absolutely and the field only got to tidy up afterwards,
    so there was nothing that could refuse.  Nobody who does this well works
    that way -- VMD/NAMD's interactive MD turns a drag into a force, Narupa
    biases the dynamics through a field whose strength the user sets, and
    Reiher's real-time quantum chemistry hands the response back as force
    feedback.  Pull, and the bond pulls back; the atom follows exactly as far
    as the balance allows, which is the flattest way out of where it stands.

    So the grabbed atom is emphatically *not* frozen -- a frozen atom has
    nothing for a force to act on -- and nothing is moved by the mouse handler
    at all.  The target moves; the atoms arrive with the relaxation.
    """
    begin = _body('ffBeginDrag')
    # A pull takes hold of atoms without freezing them.  Only what the user
    # pinned stays frozen, which is what the empty extra list means.
    assert 'ffApplyFrozen(scopeKey, [])' in begin
    assert 'ffPullApply(scopeKey);' in begin
    # Rigid is still reachable, and is also what happens when there is nothing
    # to take hold of -- a grab on an atom that is not there falls back to the
    # old path rather than silently doing nothing.
    assert 'state.pullShare > 0' in begin
    assert 'state.ffPull = want.length ? {want: want} : null;' in begin
    assert 'ffApplyFrozen(scopeKey, ffIndicesOf(viewer, targets));' in begin
    # The pull is set up before the engine is asked about, because under a
    # server method there is no engine on this page at all and the hand still
    # has to be a force -- the kernel answers it instead.
    assert begin.index('state.ffPull =') < begin.index('if (!ffEnabled(state))')

    # The mouse moves the wanted point and nothing else.
    assert 'if (state2.ffPull) {' in EDITOR
    assert 'ffPullWant(scopeKey, delta);' in EDITOR
    want = _body('ffPullWant')
    assert 'w.x += deltaWorld.x;' in want
    # The thermal leash bounds what may be *asked* for, not where the atom is:
    # under a pull the atom never goes anywhere the field did not take it.
    assert 'thermalWallBlocks(scopeKey, w, w.atom, deltaWorld)' in want

    # Letting go clears the pull before anything else, and whether or not the
    # field is still there -- one left behind would go on dragging the atom
    # towards a point the mouse abandoned the next time the engine ran.
    end = _body('ffEndDrag')
    assert end.index('ffLetGo(scopeKey);') < end.index('if (!ffEnabled(state)) return;')
    assert 'ffTether(scopeKey, []);' in _body('ffLetGo')

    # Under a pull the field owns the held atoms too, so its answer is written
    # back in full.  Under the rigid hand they are the page's and are skipped.
    assert 'var held = (!state.ffPull && state.drag && state.drag.targets)' in EDITOR

    # The engine runs in a worker as often as not, and it has to hear about
    # the pull as well -- posted before the step, so it is in effect for it.
    from delfin.dashboard.molecule_viewer import FF_WORKER_LOOP_JS
    assert "message.cmd === 'tether'" in FF_WORKER_LOOP_JS
    assert 'FF.tether(message.scope, message.list);' in FF_WORKER_LOOP_JS
    assert "ffTellWorker({cmd: 'tether'" in EDITOR


def test_how_hard_the_hand_pulls_is_the_users_to_set():
    """And it is set in units of the thing it pulls against.

    A force constant means nothing on its own; a share of a bond means
    something immediately.  Measured on the browser field: at a tenth of a
    bond a dragged hydrogen travels 1.97 A and the C-H it hangs on gives
    0.05 A -- the deformation goes into the angles and torsions, which are an
    order of magnitude softer.  At one whole bond the same drag stretches it
    to 1.68 A, and at three to 2.63 A.

    Zero is the old rigid hand, kept deliberately: placing an atom exactly
    where it is wanted is sometimes the point, and it is one click away.
    """
    setter = _body('setPullStrength')
    # A share of a bond, and the bond is the one the measurements used.
    assert 'var PULL_LIKE_A_BOND = 662.0;' in EDITOR
    assert 'state.pullShare = Math.min(3, asked);' in setter
    assert "if (!isFinite(asked) || asked < 0) asked = 0;" in setter
    # Changed mid-drag it takes effect without letting go.
    assert 'if (state.ffPull) ffPullApply(scopeKey);' in setter
    assert 'setPullStrength: setPullStrength,' in EDITOR

    # The ceiling on the force is the engine's to enforce, so it travels with
    # every pull rather than being clamped once per mouse event.
    assert 'var PULL_REACH = 1.0;' in EDITOR
    assert 'reach: PULL_REACH' in _body('ffPullApply')

    # On the toolbar, wired both ways, and re-applied when the field is
    # handed its parameters so a reload keeps the feel the user set.
    source = SUBMIT_SOURCE
    assert 'submit_pull_slider = widgets.FloatSlider(' in source
    assert "description='Pull'" in source
    assert 'submit_pull_slider.observe(on_submit_pull_changed' in source
    # Handed over with the parameters, moved by hand, and re-sent when the
    # budget goes on or off -- the page's hand has a ceiling exactly while
    # the budget does.
    assert source.count('setPullStrength(') == 3, source.count('setPullStrength(')
    # Shown under both engines.  Hidden under a server method -- as it was at
    # first, on the grounds that the pull belonged to the browser's field --
    # the one place where the budget and the scan live was the one place the
    # hand stayed absolute, which is exactly how it read.
    assert "submit_pull_slider.layout.display = ''" in source
    # Enabled and disabled with the rest of the manipulation toolbar.
    assert 'submit_pull_slider.disabled = not enabled' in source


#: The band is viewer JavaScript, so it is driven the way the force field is:
#: node runs exactly the code the browser will run, against a stand-in viewer
#: that records what was drawn.  Skipped cleanly where node is not installed.
_NODE = __import__("shutil").which("node")

_BAND_DRIVER = r"""
var shapes = [], drawn = [];
globalThis.document = {
  querySelector: function(){ return null; }, querySelectorAll: function(){ return []; },
  createElement: function(){ return {style:{}, classList:{add:function(){},remove:function(){}}}; },
  addEventListener: function(){}, body: {appendChild: function(){}}
};
globalThis.window = {
  document: document, addEventListener: function(){},
  requestAnimationFrame: function(fn){ fn(); return 1; },
  cancelAnimationFrame: function(){}, getComputedStyle: function(){ return {}; },
  HTMLInputElement: {prototype:{}}, HTMLTextAreaElement: {prototype:{}},
  Event: function(){}
};
__BOOTSTRAP__
var scope = 's1';
var atoms = [{serial:0, elem:'C', x:0, y:0, z:0}, {serial:1, elem:'C', x:1.53, y:0, z:0}];
var model = {atoms: atoms, selectedAtoms: function(){ return atoms; }};
var viewer = {
  getModel: function(){ return model; }, selectedAtoms: function(){ return atoms; },
  addLine: function(o){ var s={kind:'line',o:o}; shapes.push(s); drawn.push(s); return s; },
  addSphere: function(o){ var s={kind:'sphere',o:o}; shapes.push(s); drawn.push(s); return s; },
  removeShape: function(s){ var i=shapes.indexOf(s); if(i>=0) shapes.splice(i,1); },
  render: function(){}
};
window._submitMolViewerByScope[scope] = viewer;
var st = window._submitManipStateByScope[scope] = {
  mode:'manip', scopeKey:scope, ffActive:false, settleOnRelease:true, pullShare:0.1,
  pinned:[], ffFrameMs:16, picks:[], pivot:null, shapes:[], pivotShape:null,
  undo:[], overlay:null, viewerEl:null, canvas:null, rect:null, drag:null
};
function draw(){ drawn = []; window.__delfinSubmitManip.setPositions(scope, [0,0,0, 1.53,0,0]); }
var out = {};
draw(); out.quiet = drawn.length;
st.ffPull = {want:[{atom:0, serial:0, x:-0.4, y:0, z:0}]};
draw();
out.short = drawn.map(function(s){ return s.kind; });
out.shortEnd = drawn[drawn.length-1].o.center;
out.color = drawn[0].o.color;
st.ffPull = {want:[{atom:0, serial:0, x:-3.0, y:0, z:0}]};
draw(); out.farEnd = drawn[drawn.length-1].o.center;
out.liveWhilePulling = shapes.length;
st.ffPull = {want:[{atom:0, serial:0, x:-0.001, y:0, z:0}]};
draw(); out.aClick = drawn.length;
st.ffPull = null;
draw(); out.afterLetGo = drawn.length; out.liveAfter = shapes.length;
console.log(JSON.stringify(out));
"""


@pytest.mark.skipif(_NODE is None, reason="node not installed")
def test_the_hand_is_drawn_because_the_answer_is_slow():
    """An atom being pulled and not giving way looks like an atom that is not
    listening.

    Under a server method every frame the picture gets comes from the kernel,
    and the kernel answers about ten times a second.  Nothing else on screen
    moves at the rate of the hand, so a drag that is doing exactly what it
    should -- applying a force the chemistry is refusing -- reads as a drag
    that is broken.  Interactive molecular dynamics has drawn a band between
    the cursor and the atom since VMD, for this reason.

    It is drawn to the *clamped* point rather than to the cursor, so its
    length is the force: short while the atom is keeping up, and at full
    stretch when the hand is pulling as hard as it is allowed to.  Drawn to
    the cursor it would promise something the structure is not being asked
    for, since past the reach the pull stops growing.

    Driven in node against a stand-in viewer that records what was drawn.
    """
    import json
    import subprocess
    import tempfile

    script = _BAND_DRIVER.replace('__BOOTSTRAP__', EDITOR)
    with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as fh:
        fh.write(script)
        path = fh.name
    done = subprocess.run([_NODE, path], capture_output=True, text=True,
                          timeout=120)
    assert done.returncode == 0, done.stderr[-800:]
    out = json.loads(done.stdout.strip().splitlines()[-1])

    # Nothing held, nothing drawn.
    assert out['quiet'] == 0, out
    # Held: dashes and a mark where the hand is.  Dashes rather than a line,
    # because a solid segment between two points in a molecule reads as a
    # bond, and this is the one thing on screen that is not one.
    assert out['short'] == ['line', 'line', 'line', 'sphere'], out
    assert out['color'] == '#ff6d00', out
    # Inside the reach the mark is the wish itself.
    assert out['shortEnd'] == {'x': -0.4, 'y': 0, 'z': 0}, out
    # Past it the band stops growing, because so does the pull: three
    # angstroms of mouse and one angstrom of band.
    assert out['farEnd'] == {'x': -1, 'y': 0, 'z': 0}, out
    # A press that turns out to be a click draws nothing at all.
    assert out['aClick'] == 0, out
    # And letting go takes it away, leaving nothing behind on the viewer.
    assert out['afterLetGo'] == 0, out
    assert out['liveWhilePulling'] == 4, out
    assert out['liveAfter'] == 0, out


def test_the_band_moves_at_the_rate_of_the_hand():
    """Not at the rate of the answers, which is the whole point of it.

    Redrawn from the mouse handler as well as from every path that moves
    atoms, so it follows the cursor between one kernel frame and the next --
    and rendered there only when nothing else is about to, because with the
    browser's field running the relaxation draws the same frame a moment
    later and a second render is a second frame's worth of work.
    """
    move = EDITOR.split("if (state2.ffPull) {")[1].split('} else {')[0]
    assert 'ffPullWant(scopeKey, delta);' in move
    assert 'drawPull(scopeKey);' in move
    assert 'if (!ffEnabled(state2)) {' in move

    # Redrawn with everything else that moves atoms, so a frame from the
    # kernel takes the band with it.
    assert 'drawPull(scopeKey);' in _body('drawHighlightsNow')
    # Letting go clears it, and that path does render -- nothing else will.
    assert 'drawPull(scopeKey);' in _body('ffLetGo')
    # Nothing is rendered inside the drawing itself: every caller renders once
    # for its own reasons.
    assert 'viewer.render()' not in _body('drawPull')
    # And a new viewer starts with no band left over from the old one.
    assert 'state.pullShapes = [];' in _body('onViewerReady')


def test_the_hand_lets_go_over_a_few_frames_rather_than_at_once():
    """A structure held out of shape and then released in one frame springs.

    The field answers the whole strain in a single step, and what the user
    sees is a snap -- which is the same thing they see when the budget's wall
    fires, and the same complaint.  It ends in the same place either way,
    because the place is the nearest minimum; the difference is whether the
    way there can be watched.

    A fifth off each frame is gone in about twenty of them, a third of a
    second: slow enough to read, quick enough not to feel like the hand is
    still attached.
    """
    assert 'var PULL_FADES_BY = 0.8;' in EDITOR
    assert 'var PULL_IS_OVER = 0.02;' in EDITOR

    # Letting go hands the pull to the fade rather than dropping it, and what
    # it holds is what the hand held -- only the strength comes down.
    go = _body('ffLetGo')
    assert 'state.ffFading = {want: state.ffPull.want, share: 1.0};' in go
    assert 'state.ffPull = null;' in go
    # Unless there is no field here to ease off against, which is a server
    # method: there the kernel owns the frames and nothing on the page runs.
    assert 'if (!ffEnabled(state)) {' in go
    assert 'ffTether(scopeKey, []);' in go

    fade = _body('ffFadeStep')
    assert 'fading.share *= PULL_FADES_BY;' in fade
    assert 'if (fading.share < PULL_IS_OVER) {' in fade
    # The same clamp and the same over-factor as a live pull, so the hand on
    # the way out is the hand that was there, weaker.
    assert 'pullConstant(state) * pullOver(scopeKey, viewer) * fading.share' in fade
    assert 'reach: PULL_REACH' in fade

    # Stepped from the relaxation itself, so the structure is answering a hand
    # that is getting weaker rather than one that has vanished.
    assert 'if (state.ffFading) ffFadeStep(scopeKey);' in _body('ffRelaxFrame')
    assert 'if (state.ffFading) ffFadeStep(scopeKey);' in _body('ffRelaxAsync')


# --- the ring, and whether it says the truth about what is running ---------
#
# Measured by driving the real handlers with a queued schedule_ui_update -- the
# way Voila's io_loop queues one -- and sampling both status labels from the
# interface thread while the worker really runs on its own.  Reading the source
# cannot answer this class of question: the source says exactly what it means
# to do at every call site, and the defect was in which call sites there are.

_RING = "class='delfin-busy'"


def _a_part(xyz=_BENZENE):
    """One structure editor, and the queue its interface updates go through.

    Queued rather than called straight through, because a worker that runs to
    completion before anything is sampled proves nothing about what is on
    screen while it runs.
    """
    import collections
    import pathlib
    import tempfile

    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ("calc", "archive", "office"):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / "calc", archive_dir=room / "archive",
                           office_dir=room / "office")
    ctx.run_js = lambda _script: None
    ctx.add_init_js = lambda _script: None
    queue = collections.deque()
    state = {}
    box = widgets.Textarea(value=xyz)
    part = structure_editor.build(
        ctx, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: queue.append((func, a, k)),
        update_view=lambda *a, **k: None, get_smiles_charge=lambda *a, **k: None)
    state["current_xyz_for_copy"] = {"content": xyz}
    part.submit_ff_dd.value = "gfn2"

    def pump():
        while queue:
            func, a, k = queue.popleft()
            func(*a, **k)

    return part, state, pump


def _watch(part, state, pump, press, seconds=20.0):
    """Press something and sample both labels until the work is over.

    Each sample is (work in flight, ring in the inline label, ring in the
    fullscreen copy, the text).
    """
    import time

    samples = []
    press()
    began = time.perf_counter()
    while True:
        pump()
        inline, full = part.mol_status.value, part.mol_status_fs.value
        running = bool(state.get("busy_jobs"))
        samples.append((running, _RING in inline, _RING in full,
                        re.sub("<[^>]+>", "", inline)))
        if not running or time.perf_counter() - began > seconds:
            break
        time.sleep(0.01)
    for _ in range(3):
        pump()
        inline, full = part.mol_status.value, part.mol_status_fs.value
        samples.append((bool(state.get("busy_jobs")), _RING in inline,
                        _RING in full, re.sub("<[^>]+>", "", inline)))
        time.sleep(0.01)
    return samples


def _slow(seconds, answer):
    import time

    def _fake(*_a, **_k):
        time.sleep(seconds)
        return answer
    return _fake


def test_the_ring_turns_for_the_whole_of_a_run_and_stops_with_it():
    """A settle driven for real, sampled while xtb is being waited on.

    The settle was one of the workers that already announced itself properly,
    and this holds it to that now that the announcement no longer decides
    anything.  What the ring is decided by was the defect: sampled the same
    way, a scan showed no ring at all for its first third of a second -- the
    line on screen was the note left over from choosing GFN2, written without
    one -- and the 0.35 s each release waits before a settle or a minimisation
    starts showed none either.  Someone who reads that as finished picks an
    atom up, and then two calculations are moving the same structure.
    """
    from delfin.dashboard import gfn_optimize as _gfn

    part, state, pump = _a_part()
    real = _gfn.optimize_with_gfn
    _gfn.optimize_with_gfn = _slow(0.3, {"ok": True, "xyz": _BENZENE,
                                         "converged": True, "energy": -7.0})
    try:
        samples = _watch(part, state, pump, part._gfn_settle_now)
    finally:
        _gfn.optimize_with_gfn = real

    running = [s for s in samples if s[0]]
    assert len(running) > 5, "the settle was over before it could be sampled"
    assert all(s[1] and s[2] for s in running), (
        "the ring went out while the settle was still running: "
        f"{[s[3] for s in running if not (s[1] and s[2])][:1]}")
    assert not any(s[1] or s[2] for s in samples if not s[0]), (
        "the ring was still turning after the settle had finished")


def test_a_progress_line_cannot_take_the_ring_off_a_running_job():
    """The defect, at its smallest: 127 status writes, 23 asking for the ring.

    The other 104 turned it off, and several of them are written from inside a
    worker that is still running -- a scan saying which step it is on, an
    optimisation saying how far the last round moved.  The ring belongs to
    whether anything is running, which no single status message is in a
    position to know, so it is read from the register of running work instead.
    """
    part, state, _pump = _a_part()
    token = part._busy_begin("something long")
    try:
        part._set_mol_status("step 12 of 40")
        assert _RING in part.mol_status.value
        assert _RING in part.mol_status_fs.value
        # Including the one that clears the words: a structure arriving is not
        # the same event as the work being over.
        part._clear_mol_status()
        assert _RING in part.mol_status.value
        assert _RING in part.mol_status_fs.value
    finally:
        part._busy_end(token)
    assert state.get("busy_jobs") == {}


def test_a_worker_that_dies_stops_the_ring_instead_of_turning_for_ever():
    """The opposite failure, and just as bad: a ring that cannot be stopped.

    None of the fourteen background workers was wrapped, and each cleared its
    own run flag inside a callback an exception skips.  Measured before this:
    a saddle search whose engine raised left ``saddle_run`` true and the ring
    turning over a calculation that had stopped -- for the rest of the
    session, with the button still offering to stop it.
    """
    import time

    from delfin.dashboard import saddle as _saddle

    part, state, pump = _a_part()
    real = _saddle.optimise_to_saddle

    def _boom(*_a, **_k):
        time.sleep(0.2)
        raise RuntimeError("xtb went missing halfway through")

    _saddle.optimise_to_saddle = _boom
    try:
        samples = _watch(part, state, pump, part.on_submit_saddle)
    finally:
        _saddle.optimise_to_saddle = real

    assert any(s[0] for s in samples), "the climb never started"
    assert not state.get("busy_jobs"), "the register still says work is running"
    assert not any(s[1] or s[2] for s in samples if not s[0]), (
        "the ring outlived the worker that died")
    assert state.get("saddle_run") is False, (
        "the run flag stayed set, so the button still says Stop")
    assert "stopped on an error" in re.sub("<[^>]+>", "", part.mol_status.value)


def test_a_question_asked_of_the_page_counts_as_work_and_is_on_a_leash():
    """Reading a drawing has no thread behind it, and is still a wait.

    The molfile comes back through a hidden box the page writes into, so
    between the press and the answer there was no worker to count and the ring
    was decided by a status message alone.  A frame that has been folded away
    never answers at all -- which is why the wait is on a leash rather than
    trusted to end.
    """
    part, state, _pump = _a_part()
    part.on_submit_draw_get()
    assert state.get("busy_jobs"), "waiting on the page did not count as work"
    assert _RING in part.mol_status.value
    assert _RING in part.mol_status_fs.value

    # The page answers, and the wait is over.
    part.submit_draw_sync.value = "1\n!no-editor"
    assert not state.get("busy_jobs")
    assert _RING not in part.mol_status.value
    assert _RING not in part.mol_status_fs.value


def test_both_status_labels_always_agree_about_the_ring():
    """A ring in the inline line and none in the overlay is the same bug.

    Which of the two is on screen is the overlay's business, so a reader in
    fullscreen and a reader outside it have to be told the same thing.  The
    one line the overlay drops -- the prompt asking for coordinates, which
    would sit there for ever over a structure that is plainly on screen --
    keeps the ring even so.
    """
    part, state, _pump = _a_part()
    token = part._busy_begin("something long")
    try:
        for line in ("Optimising 1 frame(s) with GFN2-xTB...",
                     "step 12 of 40",
                     "Load a structure before optimising.",
                     "Paste or enter XYZ coordinates."):
            part._set_mol_status(line)
            assert (_RING in part.mol_status.value) is True, line
            assert (_RING in part.mol_status_fs.value) is True, line
    finally:
        part._busy_end(token)


def test_every_background_worker_is_started_through_the_one_place():
    """So that the next one added cannot forget the ring.

    Fixing the call sites one at a time leaves the next status message anyone
    writes re-introducing the defect, and the next thread anyone starts
    leaving the register untouched.  There are two raw threads in the editor
    and they are the mechanism itself: the wrapper that counts the work, and
    the leash on a question the page may never answer.
    """
    started = re.findall(r'threading\.Thread\(target=(\w+)', _EDITOR_PY)
    assert sorted(started) == ["_leash", "_run"], (
        f"a background worker bypasses _start_background: {started}")
    assert "_start_background(" in _EDITOR_PY


# ---------------------------------------------------------------------------
# the toolbar shows what the chosen method can do, and nothing else
# ---------------------------------------------------------------------------
#: What the method box offers, in the order it offers it: two force fields
#: that run in the browser, three that run xtb on the server, three that run
#: MOPAC there.
_METHODS = ('uff', 'mmff94', 'gfnff', 'gfn2', 'gxtb', 'pm6d3h4', 'pm6', 'pm7')
_BROWSER = ('uff', 'mmff94')
_XTB = ('gfnff', 'gfn2', 'gxtb')
_MOPAC = ('pm6d3h4', 'pm6', 'pm7')

_ETHANE = """8
ethane
C  0.000  0.000  0.765
C  0.000  0.000 -0.765
H  1.019  0.000  1.163
H -0.510 -0.883  1.163
H -0.510  0.883  1.163
H  1.019  0.000 -1.163
H -0.510 -0.883 -1.163
H -0.510  0.883 -1.163
"""

#: Ethene with its double bond drawn long, so perception reads it as single
#: and both carbons come back sp3 -- the case Types exists for.
_ETHENE = """6
ethene, drawn long enough that the double bond is not perceived
C  0.000  0.000  0.780
C  0.000  0.000 -0.780
H  0.930  0.000  1.320
H -0.930  0.000  1.320
H  0.930  0.000 -1.320
H -0.930  0.000 -1.320
"""

#: A six-coordinate centre, with the donors far enough out that perception
#: finds all six of them.
_COMPLEX = """7
a six-coordinate centre
Fe  0.000  0.000  0.000
F   2.000  0.000  0.000
F  -2.000  0.000  0.000
F   0.000  2.000  0.000
F   0.000 -2.000  0.000
F   0.000  0.000  2.000
F   0.000  0.000 -2.000
"""


@pytest.fixture
def editor(tmp_path):
    """One structure editor over a coordinate box of its own.

    Built rather than read: what a control does under a method is a property
    of the running editor, and the source says only what it means to do.
    """
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    sent = []
    ctx.run_js = sent.append
    state = {}

    def build(text=_ETHANE):
        box = widgets.Textarea(value=text)
        part = structure_editor.build(
            ctx, state=state, coords_widget=box, viewer_height=560,
            schedule_ui_update=lambda func, *a, **k: func(*a, **k),
            update_view=lambda *a, **k: None,
            get_smiles_charge=lambda *a, **k: None)
        state['current_xyz_for_copy'] = {'content': text}
        # What a loaded structure does to the toolbar, so the test drives the
        # editor a user would be looking at rather than an empty one.
        part._set_manip_toolbar_enabled(True)
        part.sent = sent
        return part

    return build


def _visible(widget):
    return str(getattr(widget.layout, 'display', '') or '') != 'none'


def _said(part):
    return ' '.join(re.sub('<[^>]+>', ' ', part.mol_status.value or '').split())


def test_a_type_is_a_force_fields_idea_and_is_offered_where_one_runs(editor):
    """Types, and the hybridisation box beside it, under all eight methods.

    An atom type is what a force field builds its parameters from: forcing a
    carbon to sp2 types it C_2 and its angles come back at 120 degrees, which
    is how a double bond that perception missed is put back. The methods that
    run on the server have no such thing to be told -- xtb and MOPAC work the
    shape out from the electrons at every step -- so the override has nowhere
    to land.

    Driven, on an ethene drawn long enough that the double bond is not
    perceived, both carbons typed from their partners: under UFF and MMFF94
    the press hands a fresh parameter set to the page, and under GFN-FF, GFN2
    and PM7 the same press hands over ``setForceField(scope, null)`` -- the
    field taken off -- while the status line said "2 carbons typed from their
    partners ... 2 changed". A report of a change that reached no calculation
    is the worst of the ways a control can fail, and it is why this is hidden
    rather than left to refuse on the press. Where the press lands is the test
    below; this one is only about which methods offer it at all.
    """
    part = editor(_ETHENE)
    for method in _METHODS:
        part.submit_ff_dd.value = method
        wanted = method in _BROWSER
        assert _visible(part.submit_hyb_auto_btn) is wanted, method
        # And the box does not come back with a selection, which is what
        # showed it: picking two carbons under GFN2 offered the override again.
        part.submit_pick_sync.value = '0,1'
        assert _visible(part.submit_hyb_dd) is wanted, method
        part.submit_pick_sync.value = ''


def test_a_type_forced_under_a_server_method_reached_nothing(editor):
    """The measurement the rule above rests on.

    Kept as a test of its own because it is the fact rather than the
    consequence: if a server method ever does start reading the override, this
    is what will say so.

    On the ethene above, so the press has a real correction to make rather
    than a type to confirm.
    """
    pytest.importorskip('rdkit')
    call = re.compile(r'window\.__delfinSubmitManip\.setForceField\('
                      r'"[^"]*",(.{0,8})', re.S)

    def where_it_landed(method):
        part = editor(_ETHENE)
        part.submit_ff_dd.value = method
        part.sent.clear()
        part.state['hyb_overrides'] = {}
        part.on_submit_hyb_auto()
        return [m.group(1).startswith('null')
                for script in part.sent for m in call.finditer(script)]

    for method in _BROWSER:
        assert where_it_landed(method) == [False], (
            f'{method} must be handed the parameters the types produced')
    for method in ('gfnff', 'gfn2') + _MOPAC:
        assert where_it_landed(method) == [True], (
            f'{method} takes the field off instead, so the types went nowhere')


def test_the_polyhedron_is_offered_where_its_restraints_can_act(editor):
    """A polyhedron is a set of restraints, and they are terms in one field.

    Choosing one installs pulls that draw the donors onto its vertices, and
    those pulls live in the browser's own field. Driven on a six-coordinate
    centre under GFN2 and PM7, choosing an octahedron said "the donors are
    pulled onto it" and handed the page ``setForceField(scope, null)`` in the
    same breath -- the field taken off, the promise unkept.

    Turn goes with it: it steps between arrangements of a polyhedron that is
    not acting. Swap stays, because it is not a restraint but an edit -- the
    page rotates the two ligands onto each other's directions there and then,
    and every engine is handed the geometry that results.
    """
    part = editor(_COMPLEX)
    for method in _METHODS:
        part.submit_ff_dd.value = method
        part.submit_pick_sync.value = '0'
        wanted = method in _BROWSER
        assert _visible(part.submit_poly_dd) is wanted, method
        if not wanted:
            # And no explanation blaming the coordination number, which would
            # send the user looking for a table that is not the reason.
            assert 'polyhedron table' not in _said(part), method
        part.submit_pick_sync.value = ''


def test_the_saddle_search_is_offered_where_orca_can_drive_the_method(editor):
    """It is ORCA's optimiser on somebody's gradients, and the button is
    offered wherever ORCA can be told how to ask for them.

    That was once the ORCA keywords alone -- XTB1, XTB2, XTBFF -- and under
    the other five methods the button answered "a saddle search here runs on
    xtb through ORCA, so choose GFN2, GFN1 or GFN-FF", which is a button that
    promises what it cannot do.

    g-xTB has no keyword and never will: it is a build of its own, and an
    ordinary xtb accepts --gxtb and silently runs GFN2. It is driven through
    ExtOpt instead, ORCA's own interface for a program it does not know, and
    the button is offered under it too. Under UFF, MMFF94 and the three MOPAC
    methods there is still nothing to drive, and the button is not there.

    Read from the table the run itself reads, so the button and the refusal
    cannot drift apart.
    """
    from delfin.dashboard import saddle as _saddle

    part = editor()
    for method in _METHODS:
        part.submit_ff_dd.value = method
        assert _visible(part.submit_saddle_btn) is (
            method in _saddle.SADDLE_METHODS), method
    # Which, for the methods this box offers, is these three and no others.
    assert [m for m in _METHODS if m in _saddle.SADDLE_METHODS] == [
        'gfnff', 'gfn2', 'gxtb']


def test_path_from_here_no_longer_promises_a_path_nothing_can_walk(editor):
    """Both ends of a path are walked by xtb's own path finder.

    Find the path always refused anything else; Path from here asked nothing
    at all. Driven under UFF: the first press was accepted -- "marked as the
    start of a path (8 atoms) ... press Find the path" -- and put Find the
    path on the toolbar, and the second press answered "a path needs xtb:
    choose a GFN method". Two presses to be told the first could not have
    worked.

    The mark itself survives a change of method, the way an armed scan does:
    it describes two structures, not a program.
    """
    part = editor()
    for method in _METHODS:
        part.submit_ff_dd.value = method
        assert _visible(part.submit_path_from_btn) is (method in _XTB), method


def test_keep_bonds_is_offered_where_a_step_can_be_taken_back(editor):
    """The switch works by watching what a follow step hands back.

    xtb cannot be told to hold a topology, so a step that made or broke a bond
    is replaced by the last one that did not -- which needs a follow step, and
    that is the kernel's. It runs for a server method and for nothing else:
    under UFF the drag never leaves the browser, ``_begin_gfn_follow`` answers
    no, and the wall is never consulted. Measured that way, with Dynamik Opt
    on and the switch down: False under UFF, True under GFN2 and PM7.
    """
    part = editor()
    for method in _METHODS:
        part.submit_ff_dd.value = method
        assert _visible(part.submit_topology_btn) is (
            method not in _BROWSER), method
    # And the value is left where it stands rather than switched off, so a
    # detour through UFF does not cost the setting.
    part.submit_ff_dd.value = 'gfn2'
    part.submit_topology_btn.value = True
    part.submit_ff_dd.value = 'uff'
    assert part.submit_topology_btn.value is True
    part.submit_ff_dd.value = 'gfn2'
    assert _visible(part.submit_topology_btn)


def test_keep_bonds_does_not_tell_a_pm_user_it_changes_nothing(editor):
    """The sentence was said for everything that was not GFN2 and its
    relatives, which put it under the PM methods too.

    It is not true of them. MOPAC decides the bonding from the electrons like
    any other semiempirical method, and the wall reads what came back and
    takes the step away exactly as it does under GFN2. The line told the user
    that the switch they had just pressed does nothing, while it was working.

    GFN-FF is the one that really does keep its bonding -- it reads its
    topology once -- and it is the one that still says so.
    """
    for method, keeps in (('gfnff', True), ('gfn2', False), ('pm7', False)):
        part = editor()
        part.submit_ff_dd.value = method
        part.submit_topology_btn.value = True
        assert ('already keeps its bonding' in _said(part)) is keeps, method


def test_the_pull_hand_is_not_offered_where_it_is_a_placement(editor):
    """Two hands, and under MOPAC only one of them exists.

    A pull is a force on an internal coordinate, so it needs an engine that
    can be told to hold one: the browser's field can, and xtb can through its
    constrain block. MOPAC takes no held internals from this editor at all, so
    the follow step falls through to the rigid hand.

    Measured on an ethane with one hydrogen dragged 0.60 A along x, the hand
    set to pull: GFN2 left the atom 0.4399 A short of where the cursor asked,
    which is the chemistry having its say, and PM7 put it 0.0000 A from the
    cursor -- the same geometry the move hand gives, to the last decimal.
    "Pull with a force" under MOPAC was the move hand under another name, and
    the difference was invisible from the outside.
    """
    part = editor()
    for method in _METHODS:
        part.submit_ff_dd.value = method
        offered = [value for _label, value in part.submit_hand_dd.options]
        assert offered == (['move'] if method in _MOPAC
                           else ['pull', 'move']), method
        if method in _MOPAC:
            assert part.submit_hand_dd.value == 'move'
            assert not _visible(part.submit_pull_slider)


def test_a_detour_through_a_pm_method_gives_the_hand_back(editor):
    """Switching methods is something a user does constantly, and it must not
    cost a setting each time.

    The pull is taken away under MOPAC because it cannot act there; it is
    handed back on the way out, which is a different thing from it never
    having been taken.
    """
    part = editor()
    part.submit_ff_dd.value = 'gfn2'
    part.submit_hand_dd.value = 'pull'
    part.submit_ff_dd.value = 'pm7'
    assert part.submit_hand_dd.value == 'move'
    part.submit_ff_dd.value = 'gfn2'
    assert part.submit_hand_dd.value == 'pull'
    assert _visible(part.submit_pull_slider)


def test_the_pull_slider_belongs_to_the_hand_and_not_to_the_method(editor):
    """It was shown unconditionally, so it came back on every change of
    method: pick the move hand under UFF, switch to GFN2, and there it was
    again -- a control for setting the strength of a hand that is not in use.
    """
    part = editor()
    part.submit_ff_dd.value = 'uff'
    part.submit_hand_dd.value = 'move'
    assert not _visible(part.submit_pull_slider)
    part.submit_ff_dd.value = 'gfn2'
    assert not _visible(part.submit_pull_slider)
    assert part.submit_hand_dd.value == 'move'
    part.submit_hand_dd.value = 'pull'
    assert _visible(part.submit_pull_slider)


def test_a_hidden_control_stays_hidden_when_a_structure_loads(editor):
    """Two owners for one attribute is how a control comes back under a method
    that cannot run it.

    ``_set_manip_toolbar_enabled`` writes ``disabled`` over the whole toolbar
    every time a structure arrives, so a method-based *disable* would be
    undone by the next load. The method rule uses ``display`` for that reason,
    and this is the check that the two do not fight.
    """
    part = editor()
    part.submit_ff_dd.value = 'pm7'
    hidden = [part.submit_saddle_btn, part.submit_path_from_btn,
              part.submit_hyb_auto_btn]
    assert not any(_visible(w) for w in hidden)
    part._set_manip_toolbar_enabled(False)
    part._set_manip_toolbar_enabled(True)
    assert not any(_visible(w) for w in hidden), (
        'a structure loading must not bring back what the method cannot use')


_STOP_DRIVER = r"""
// The rule as it ships, with a stand-in for the one thing it reads.
var STATE = {s1: {thermalWall: null, thermalReach: 0}};
function getState(key) { return STATE[key]; }
__WALL__
__SET__
var scope = 's1';
var out = {free: [], stopped: [], back: null, released: null};
var atom = {x: 0, y: 0, z: 0};
var step = {x: 0.2, y: 0, z: 0};
var backwards = {x: -0.2, y: 0, z: 0};
function move(delta) {
  if (thermalWallBlocks(scope, atom, 0, delta)) return;
  atom.x += delta.x; atom.y += delta.y; atom.z += delta.z;
}
function where() { return Math.round(atom.x * 1e6) / 1e6; }

// Nothing marked: the hand takes the atom wherever it goes.
for (var i = 0; i < 6; i++) { move(step); out.free.push(where()); }

// The budget refuses.  The atom is marked where it was last allowed to
// stand, with the reach the hand had already run past it.
atom.x = 0;
setThermalWall(scope, {0: [0, 0, 0]}, 0.5);
for (var j = 0; j < 6; j++) { move(step); out.stopped.push(where()); }

// Easing back is never refused, however far out it starts.
move(backwards);
out.back = where();

// And with the stop down the hand has the structure again.
setThermalWall(scope, null, 0);
move(step);
out.released = where();
console.log(JSON.stringify(out));
"""


@pytest.mark.skipif(_NODE is None, reason="node not installed")
def test_a_refused_drag_comes_up_against_a_stop_it_can_push_at():
    """What a refusal has to do to the picture, driven in node.

    The kernel stops answering a drag it has refused -- there is nothing left
    to compute until the hand eases -- and nothing told the page.  So the page
    went on running the wish out with the cursor: the band grew, the
    coordinates it reported went on changing, and the one thing that had
    happened, that the drag had reached its ceiling, was the one thing not on
    screen.

    The rule was already here and nothing ever armed it: ``_push_thermal_wall``
    was called in one place in the whole editor, with ``None``.  A marked atom
    may stand at most a reach from where the budget last agreed, and a step is
    refused only when it is both outside that reach and going further out.  So
    the hand comes up against something it can push at, and pushing harder
    does nothing -- which is the truth about where it is.

    Coming back in is never refused, however far out it starts.  That matters
    twice over: it is the one thing a user must be able to do from a structure
    the budget will not pay for, and it is exactly the gesture the kernel's own
    ``_still_spent`` reads as the hand easing off.  The two sides then let the
    same move through instead of holding two opinions about it.

    Driven in node over the shipped source of the two functions, with the mark
    at the origin, half an angstrom of reach and a hand moving a fifth of an
    angstrom a step.
    """
    import json
    import subprocess
    import tempfile

    script = (_STOP_DRIVER
              .replace('__WALL__', _body('thermalWallBlocks'))
              .replace('__SET__', _body('setThermalWall')))
    with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as fh:
        fh.write(script)
        path = fh.name
    done = subprocess.run([_NODE, path], capture_output=True, text=True,
                          timeout=120)
    assert done.returncode == 0, done.stderr[-800:]
    out = json.loads(done.stdout.strip().splitlines()[-1])

    # Nothing marked: six steps of a fifth of an angstrom go where they are put.
    assert out['free'] == [0.2, 0.4, 0.6, 0.8, 1.0, 1.2], out
    # Marked: the atom walks out to the reach and stops there.  The last step
    # inside it is taken whole -- a stop is a place, not a rounding.
    assert out['stopped'] == [0.2, 0.4, 0.4, 0.4, 0.4, 0.4], out
    # Pushing at it does not accumulate: four more refused steps leave the atom
    # exactly where the first refusal found it.
    assert out['stopped'][-1] == out['stopped'][1], out
    # Easing back is applied.
    assert out['back'] == pytest.approx(0.2, abs=1e-9), out
    # With the stop down the hand has the structure again.
    assert out['released'] == pytest.approx(0.4, abs=1e-9), out
