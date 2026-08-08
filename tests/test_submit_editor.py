"""Behavioural invariants of the Submit-tab molecule editor.

These lock in the interaction model that was verified by driving the real
editor in a headless browser: one toggle per click, a rubber band that projects
instead of dispatching synthetic events, direct atom dragging, and a camera
that can still be turned while manipulating.
"""
import re

import pytest

from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

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
    """Avogadro's manipulation: the grabbed atom follows the cursor exactly and
    the molecule settles around it. The relaxation runs in the browser because
    a per-frame round trip to the kernel costs 45 ms."""
    # Frozen atoms are a gradient mask on the engine side, set once per drag.
    assert 'ffBeginDrag(scopeKey, state.drag.targets)' in EDITOR
    assert 'ffApplyFrozen(scopeKey, ffIndicesOf(viewer, targets))' in EDITOR
    # The relaxation runs after the grabbed atoms are placed, not before.
    grab_then_relax = EDITOR.index('applyTranslate(scopeKey, delta, d.targets)')
    assert EDITOR.index('ffRelaxFrame(scopeKey)', grab_then_relax) > grab_then_relax
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
    body = source.split('def on_submit_optimize(change=None, every_frame=False)')[1].split('\n    def ')[0]
    assert "frames = list(state.get('isomers') or []) if every_frame else []" in body, (
        "Optimize takes the frame on screen; all takes the set"
    )
    assert "state['pre_optimize_frames']" in body
    assert 'relax_xyz(' in body and 'max_steps=500' in body
    # A 500-step minimisation per frame takes seconds; it must not block the UI.
    assert 'threading.Thread(target=_work, daemon=True)' in body
    assert 'worker.start()' in body
    # One bad frame must not lose the others.
    assert 'results.append(item)' in body

    undo = source.split('def on_submit_manip_undo')[1].split('\n    def ')[0]
    assert "state.pop('pre_optimize_frames', None)" in undo
    assert "state['isomers'] = snapshot['isomers']" in undo


def test_optimize_toggle_runs_the_field_continuously():
    """Avogadro's auto optimisation keeps the force field running on a timer,
    so the structure settles whether or not the mouse is down and a grabbed
    atom drags the relaxing molecule along."""
    assert 'function startAutoOptimize' in EDITOR
    assert 'function stopAutoOptimize' in EDITOR
    tick = _body('autoOptimizeTick')
    # A drag already relaxes in its own handler; a second batch here would
    # double the step budget for that frame.
    assert 'if (!state.drag) {' in tick
    assert 'window.requestAnimationFrame' in tick
    # The coordinate box follows at a readable rate -- each push is a round trip.
    assert "state.autoPushed" in tick

    stop = _body('stopAutoOptimize')
    assert 'cancelAnimationFrame' in stop
    assert 'pushXyzToPython(scopeKey)' in stop

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


def test_toolbar_wraps_instead_of_clipping_its_own_controls():
    """On a laptop the row was wider than the panel, and nowrap plus
    overflow hidden simply cut off whatever did not fit."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
    toolbar = source.split('submit_manip_toolbar = widgets.HBox')[1].split(')\n')[0]
    assert "flex_flow='row wrap'" in toolbar
    assert "overflow='hidden'" not in toolbar
    # The status line takes a share of the row when there is room -- so
    # fullscreen stays one line -- and wraps when there is not.
    status = source.split('submit_manip_status = widgets.HTML')[1].split(')\n')[0]
    assert "flex='1 1 260px'" in status

    # The label, the value box and Set are one wrapping unit: a label that
    # lands on a different row from its field explains nothing.
    group = source.split('submit_internal_group = widgets.HBox')[1].split(')\n')[0]
    for name in ('submit_internal_label', 'submit_internal_value',
                 'submit_internal_btn', 'submit_hold_btn'):
        assert name in group, name
    assert "flex_flow='row nowrap'" in group

    # Force field, then its strength, then the one-shot run, then the
    # continuous one, then the internal-coordinate group.
    children = source.split('submit_manip_toolbar = widgets.HBox')[1].split(']')[0]
    order = [
        'submit_ff_dd', 'submit_strength_slider',
        'submit_optimize_btn', 'submit_relax_btn', 'submit_internal_group',
    ]
    positions = [children.index(name) for name in order]
    assert positions == sorted(positions), children


def test_relaxation_strength_is_adjustable():
    """At full strength the field moves the structure so far per frame that a
    dragged atom has to be fought rather than led."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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
    sync input had left the scope too."""
    finder = _body('findInScope')
    assert 'document.querySelector(selector)' in finder
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
    right = manipulate.split('if (e.button === 2) {')[1][:900]
    assert right.index('probeClickAtom') < right.index('e.preventDefault();')
    assert 'if (!picked) return;' in right


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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
    assert 'submit_settle_btn' in source
    assert 'setSettleOnRelease' in source
    # On by default: it is what keeps a submitted geometry sane.
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
    # held atoms have to be handed over rather than read back.
    assert 'ffEndDrag(scopeKey, d.targets)' in EDITOR
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
    sync = source.split('def on_submit_manip_sync')[1].split('\n    def ')[0]
    assert "state['poly_assignment'] = None" in sync
    assert '_schedule_ui_update(_enable_live_forcefield)' in sync


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
    assert "serializeXyz(viewer, reason ? ('DELFIN ' + reason) : null)" in push
    # The heartbeat inside the relaxation loop carries no reason.
    tick = _body('autoOptimizeTick')
    assert 'pushXyzToPython(scopeKey);' in tick
    assert "drag-end" not in tick
    # A drag that ends does, whether it settles or not.
    assert "pushXyzToPython(scopeKey, 'drag-end')" in _body('ffEndDrag')
    assert "pushXyzToPython(scopeKey, 'drag-end')" in _body('settleAfterDrag')

    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
    sync = source.split('def on_submit_manip_sync')[1].split('\n    def ')[0]
    assert "lines[1].strip() == 'DELFIN drag-end'" in sync
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
    assert 'submit_ff_notes' in source
    assert "submit_ff_notes.add_class('submit-ff-notes')" in source
    # Below the copy row, in the panel the viewer lives in.
    children = source.split('mol_output, isomer_nav_row, xyz_copy_row')[1][:120]
    assert 'submit_ff_notes' in children

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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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
    # from and what every geometry operation here already reads.
    assert 'atoms[a].bonds.push(b)' in edit
    assert 'list.splice(at, 1)' in edit
    assert 'atoms[a].bondOrder' in edit
    assert 'snapshotForUndo(scopeKey)' in edit
    assert "pushXyzToPython(scopeKey, 'drag-end')" in edit

    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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
    assert 'redrawHighlights(scopeKey)' in body
    # Narrower than Clear: the pivot and the held atoms are not touched.
    assert 'unpinAll' not in body
    assert 'state.pivot' not in body

    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
    handler = source.split('def on_submit_cmd')[1].split('\n    def ')[0]
    assert "if verb == 'undo':" in handler

    step = source.split('def _undo_structure')[1].split('\n    def ')[0]
    assert "state['structure_undo']" in step
    assert 'coords_widget.value = snapshot' in step

    apply_ = source.split('def _apply_structure')[1].split('\n    def ')[0]
    # The previous structure is remembered before the new one is written, and
    # the write must not be mistaken for a different molecule arriving.
    assert "state['structure_undo'] = history[-_STRUCTURE_UNDO_LIMIT:]" in apply_
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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
    # And it moves when the script does.
    import hashlib
    assert match.group(1) == hashlib.sha256(
        SUBMIT_MANIP_BOOTSTRAP_JS.encode('utf-8')).hexdigest()[:12]


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

    source = open(tab_submit.__file__, encoding='utf-8').read()
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
