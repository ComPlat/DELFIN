"""The lowest structure a hand-search has been through, kept and gone back to.

Asked for: "ich brauch ein button mit dem man immer die struktur mit der
niedrigsten energie sammeln kann ... wenn ich eigenstaendig den
konformationsraum abtaste ... und das ich die struktur dann auch rein holen
kann wieder in den editor".

Searching a conformer space by hand is a walk through minima, and the one
that matters is usually not the last one reached: you drag, let go, watch it
settle, and go on -- and the good one is three gestures back, in a history
that also holds every gesture in between.

Nothing here is computed.  Every number was in hand when the answer was
written, so a search that never presses the button pays nothing for it.  What
it costs is the button, and only once there is somewhere to go back to.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest


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


def _geo(d, comment='Settled with GFN2-xTB'):
    return (f'2\n{comment}\n'
            f'C 0.000000 0.000000 0.000000\nO {d:.6f} 0.000000 0.000000\n')


def _shown(widget):
    return (widget.layout.display or '') != 'none'


def _a_search(part, walk=((1.20, -100.0000), (1.15, -100.0050),
                          (1.30, -100.0010))):
    """A hand-search: three releases, three minima, the best in the middle."""
    for where, energy in walk:
        part._keep_the_best(_geo(where), energy, 'settle')
        part.coords_widget.value = _geo(where)
    part._refresh_the_best()
    return walk


def test_nothing_is_offered_until_a_search_has_found_something():
    part, _state = _an_editor(_geo(1.20, 'start'))
    assert not _shown(part.submit_best_btn)


def test_the_lowest_is_the_one_kept_and_not_the_last():
    """Which is the whole point. The last minimum reached is where the walk
    happens to have stopped; the lowest is what the walk was for."""
    part, state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    assert _shown(part.submit_best_btn)
    assert state['best_kept']['energy'] == pytest.approx(-100.0050)
    assert state['best_kept']['seen'] == 3


def test_the_press_says_how_far_below_you_it_is():
    """The gap is against the geometry in the box and not against the worst
    thing ever seen: what a person wants to know is whether going back is
    worth it from where they are."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    said = part.submit_best_btn.description
    assert 'Best of 3' in said, said
    assert '-2.5 kcal/mol' in said, said


def test_standing_on_the_best_is_said_rather_than_offered():
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part, walk=((1.30, -100.0010), (1.15, -100.0050)))
    assert 'you are on it' in part.submit_best_btn.description
    assert part.submit_best_btn.disabled


def test_going_back_is_a_move_like_any_other():
    """Into the box and onto the screen together, with one step in the
    history: the structure that was being worked on is not lost by looking at
    a better one."""
    part, state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    left = part.coords_widget.value

    part.on_submit_best()
    assert (part._geometry_key(part.coords_widget.value)
            == part._geometry_key(_geo(1.15)))
    said = ' '.join(state.get('mol_status_lines') or ())
    assert '3 minima' in said, said

    part.on_submit_manip_undo()
    assert part.coords_widget.value == left


def test_a_different_question_is_not_this_search_s_record():
    """A total energy is an answer to a question -- the method, the charge,
    the spin, the solvent, and any value being held. Change one and the
    lowest number so far stops meaning "the best arrangement" and starts
    meaning "the cheapest question asked"."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    assert _shown(part.submit_best_btn)

    part.submit_gfn_charge.value = 1
    part._refresh_the_best()
    assert not _shown(part.submit_best_btn)

    part.submit_gfn_charge.value = 0
    part._refresh_the_best()
    assert _shown(part.submit_best_btn), 'and the record is still there'


def test_a_held_value_is_part_of_the_question():
    """A structure relaxed with a bond pinned at 2.5 A is a minimum of a
    different surface from the same molecule left free."""
    part, state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    state['constraints'] = [{'kind': 'distance', 'atoms': [0, 1],
                             'mode': 'fix', 'value': 2.5}]
    part._refresh_the_best()
    assert not _shown(part.submit_best_btn)


def test_the_record_belongs_to_the_molecule_it_was_made_on():
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    part.coords_widget.value = '3\nother\nO 0 0 0\nH 1 0 0\nH 0 1 0\n'
    part._refresh_the_best()
    assert not _shown(part.submit_best_btn)


def test_only_a_free_relaxation_is_collected():
    """A geometry standing under a hand or a push carries the restraint's own
    energy in its total: it is a number about the hand as much as about the
    structure, and kept as a best it would claim something it is not.

    Read off the two places that collect, which are the two a free relaxation
    lands in -- a release settling and Optimise reaching a minimum -- and off
    the arithmetic they price with, which takes a held value's own energy back
    out.
    """
    from editor_source import EDITOR_SOURCE

    assert EDITOR_SOURCE.count('_keep_the_best(') == 3     # two calls, one def
    settle = EDITOR_SOURCE.split("f'Settled with {label}'", 1)[1][:600]
    assert '_keep_the_best(' in settle, 'a release does not reach the search'
    press = EDITOR_SOURCE.split("state['gfn_energy_unit']", 1)[1][:800]
    assert '_keep_the_best(' in press and '_settle_price(' in press
    assert 'if not _stopped():' in press, (
        'a stopped run keeps the frame on screen, and the energy in hand is '
        'about the whole run rather than that frame')
