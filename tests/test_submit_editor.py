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
