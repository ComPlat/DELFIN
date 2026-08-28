"""How far the structure moves for how far the mouse moves.

One is the cursor and the atom staying together: let go and the atom is where
the pointer is, which is what makes placing one by hand trustworthy. Below one
the hand travels further than the structure, which is what a crowded region
wants; above one it travels less far, for reaching across a large system
without running out of desk.

It scales dragging and turning a selection, and nothing else. Where a *new*
atom is placed stays one to one -- an atom that appeared somewhere other than
under the cursor would be a different kind of wrong, and that is the one thing
about this editor that was hardest to get right.

Verified in a browser against the real module: at 1.0 a 40 px right, 25 px down
drag with pixelToWorld 0.05 gives {x: 2, y: -1.25}, byte for byte what the code
produced before the factor existed; 0.5 halves it and 2.0 doubles it. Driven
through the running viewer, setDragSensitivity returned and stored 1, 0.5 and
2.5 unchanged, and turned 0, 'abc' and 99 into 1, 1 and 5.
"""

from delfin.dashboard import molecule_viewer, tab_submit
from editor_source import SUBMIT_SOURCE


VIEWER = open(molecule_viewer.__file__, encoding='utf-8').read()
SUBMIT = SUBMIT_SOURCE


def _setter():
    return VIEWER.split('function setDragSensitivity(scopeKey, value) {')[1] \
                 .split('\n    }')[0]


def test_anything_unusable_falls_back_to_one():
    """Zero, negative and nonsense all mean 'do not scale', not 'freeze'."""
    body = _setter()

    assert 'if (!isFinite(asked) || asked <= 0) asked = 1;' in body
    assert 'Math.max(0.1, Math.min(5, asked))' in body


def test_it_scales_the_drag():
    move = VIEWER.split("if (d.kind === 'translate' && d.movedEnough)")[1] \
                 .split('} else if')[0]

    assert '* (state2.dragSensitivity || 1)' in move
    # The fallback matters: the state exists before the slider has ever been
    # applied, and an undefined factor must mean one, not zero.
    assert 'state2.dragSensitivity || 1' in move


def test_it_scales_turning_a_selection_too():
    turn = VIEWER.split("if (d.kind === 'rotate' && d.movedEnough)")[1] \
                 .split('} else if')[0]

    assert 'state2.dragSensitivity || 1' in turn
    assert 'applyRotate(scopeKey, dx * turn, dy * turn)' in turn


def test_placing_an_atom_is_left_alone():
    """The draw path must not be scaled: a new atom belongs under the cursor."""
    draw = VIEWER.split("d.kind === 'draw'")[1][:4000]

    assert 'dragSensitivity' not in draw


def test_the_viewer_offers_it_and_keeps_it_across_a_reload():
    assert 'setDragSensitivity: setDragSensitivity,' in VIEWER
    # Re-applied where setOptimizerStrength is, so a rebuilt page keeps the
    # feel the user set.
    reapply = VIEWER.split('setOptimizerStrength(scopeKey, state.ffStrength);')[1][:200]
    assert 'setDragSensitivity(scopeKey, state.dragSensitivity)' in reapply


def test_the_tab_puts_it_behind_strength():
    # On the picture now, with the other controls that act on the view and on
    # the feel of the hand rather than on the structure -- and still behind
    # the strength, which is the order this asks about.
    panel = SUBMIT.split('submit_view_body = widgets.VBox')[1].split(')\n')[0]

    assert 'submit_sens_slider' in panel
    assert panel.index('submit_strength_slider') < panel.index(
        'submit_sens_slider')

    made = SUBMIT.split('submit_sens_slider = widgets.FloatSlider')[1].split(')\n')[0]
    assert 'value=1.0' in made, 'the default has to be what it did before'
    assert 'min=0.2' in made and 'max=3.0' in made


def test_moving_it_reaches_the_viewer():
    handler = SUBMIT.split('def on_submit_sens_changed')[1].split('\n    def ')[0]

    assert 'setDragSensitivity(' in handler
    assert 'float(submit_sens_slider.value)' in handler
    assert 'submit_sens_slider.observe(on_submit_sens_changed' in SUBMIT
