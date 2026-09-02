"""Reset restores the structure.  It must not put the tools away with it.

The write Reset makes goes through the host's "a structure I have not seen"
path, and that path switches Manipulate, Select, Draw and Load back off.  For
a structure that has just arrived that is right, and the host's own note
argues it: the first click would otherwise land on an atom nobody meant to
move.  Reset is not that.  It goes back to the same molecule the user is
working on, and taking away the mode they are working in leaves a picture that
does not answer a click.

Reported twice inside an hour -- "der viewer ist nach reset einfach
eingefroren" and "warum reagiert nach reset der viewer nicht mehr" -- from a
session whose journal shows the press, then ``submit_manip_btn = False``, and
then nothing the hand did reaching anything at all.

The same reasoning the host already applies to the method and the solvent:
they are how the user is working, and a structure arriving does not stop them
being that.
"""
from __future__ import annotations

import pathlib

_EDITOR = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()


def _reset():
    return _EDITOR.split('def on_submit_reset(')[1].split('\n    def ')[0]


def test_the_four_mode_buttons_are_put_back_after_the_write():
    body = _reset()
    assert 'modes_were = [(one, one.value) for one in' in body
    for one in ('submit_select_btn', 'submit_manip_btn',
                'submit_draw_btn', 'submit_load_btn'):
        assert one in body.split('modes_were')[1][:260], one
    # Read before the write and put back after it, or the write would take
    # the values with it.
    assert body.index('modes_were = [') < body.index('_clear_selection()')
    assert body.index('for widget, was in modes_were:') > \
        body.index('_clear_selection()')


def test_and_what_is_drawn_on_the_structure_goes_with_it():
    """An arrow survives a redraw on purpose -- it is about atoms by number,
    and the structure that arrives has the same ones -- so nothing in the
    write takes one away.  After a Reset the pull vectors were still there
    over a structure that had never been pulled."""
    body = _reset()
    assert "state['loads'] = []" in body
    assert '_tell_the_page_the_loads()' in body
    assert "state['pivot'] = None" in body
