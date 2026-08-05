"""Behavioural invariants of the Submit-tab molecule editor.

These lock in the interaction model that was verified by driving the real
editor in a headless browser: one toggle per click, a rubber band that projects
instead of dispatching synthetic events, direct atom dragging, and a camera
that can still be turned while manipulating.
"""
import re

from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

EDITOR = submit_manip_bootstrap_js()


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
    assert 'window.__delfinFF.grab(scopeKey, ffIndicesOf(viewer, targets))' in EDITOR
    # The relaxation runs after the grabbed atoms are placed, not before.
    grab_then_relax = EDITOR.index('applyTranslate(scopeKey, delta, d.targets)')
    assert EDITOR.index('ffRelaxFrame(scopeKey)', grab_then_relax) > grab_then_relax
    assert 'window.__delfinFF.release(scopeKey)' in EDITOR

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
    # Every edit is undoable and reaches the coordinate box.
    assert 'snapshotForUndo(scopeKey)' in body
    assert 'pushXyzToPython(scopeKey)' in body
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
