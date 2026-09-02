"""The mode field offers what the armed legs can be walked with, and no more.

A push is a force between two atoms, so it walks distances.  There is no force
that drives an *angle* towards a reaction -- what one would push on is the two
bonds that make it, which is what the atoms at its ends already are.  The run
said so, and said it after the leg was armed and the button pressed, which is
the worst moment for it: the work is done and then comes the correction.

So it comes off the field instead, exactly the way ``form`` and ``break``
already come off the direction field for an angle.  When the legs go back to
being distances the three are offered again.
"""
from __future__ import annotations

import pathlib

_EDITOR = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()


def _refresh():
    return _EDITOR.split('def _refresh_scan(')[1].split('\n    def ')[0]


def test_a_leg_that_is_not_a_distance_takes_the_push_off_the_field():
    body = _refresh()
    assert "turning = any(str(one.get('kind')) != 'distance' for one in legs)" \
        in body
    # Walking is what is left, and it is what the field then holds.
    assert "[('walk the value', 'hold')] if turning else" in body
    assert "if submit_scan_how.value not in [v for _l, v in offered]:" in body
    assert "submit_scan_how.value = offered[0][1]" in body


def test_and_the_three_come_back_when_the_legs_are_distances_again():
    body = _refresh()
    after = body.split('else')[1] if 'else' in body else body
    for one in ('push with a force', 'walk the value', 'pull along the arrows'):
        assert one in body, one


def test_the_field_narrowing_itself_does_not_call_the_refresh_again():
    """The refresh sets the field; the field's own handler calls the refresh."""
    body = _refresh()
    assert "state['scan_how_quiet'] = True" in body
    handler = _EDITOR.split('def on_submit_scan_how(')[1].split('\n    def ')[0]
    assert "if state.get('scan_how_quiet'):" in handler
    assert handler.index("scan_how_quiet") < handler.index('_refresh_scan()')


def test_the_run_still_refuses_it_if_it_ever_gets_there():
    """The field is the way it is prevented; the refusal is the floor under it.

    A field can be got round -- a scan armed under one mode and the legs
    changed under another -- and the reason it cannot work does not depend on
    which control was used to ask for it.
    """
    run = _EDITOR.split('def on_submit_scan_run(')[1].split('\n    def ')[0]
    assert "if pushing and any(one['kind'] != 'distance' for one in legs):" \
        in run
