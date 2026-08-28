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
        part._keep_a_point(_geo(where), energy, 'minimum', 'settle')
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
    assert part._best_here()['energy'] == pytest.approx(-100.0050)
    assert len(part._points_here()) == 3, 'three distinct minima'
    assert state


def test_the_press_says_how_far_below_you_it_is():
    """The gap is against the geometry in the box and not against the worst
    thing ever seen: what a person wants to know is whether going back is
    worth it from where they are."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    # The press is named, not measured: a count and a difference in the label
    # make the button a readout that changes width every time an answer lands.
    # The numbers are on the list beside it, where they can be compared, and
    # in the tooltip of the press itself.
    assert part.submit_best_btn.description == 'Best of'
    said = part.submit_best_btn.tooltip
    assert 'Best of 3' in said, said
    assert '-2.5 kcal/mol' in said, said


def test_standing_on_the_best_is_said_rather_than_offered():
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part, walk=((1.30, -100.0010), (1.15, -100.0050)))
    assert part.submit_best_btn.description == 'Best of'
    assert 'You are on it' in part.submit_best_btn.tooltip
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

    assert EDITOR_SOURCE.count('_keep_a_point(') == 4   # three calls, one def
    # Sliced to the next thing that happens rather than to a count of
    # characters: comments grow, and a window measured in bytes turns a rule
    # about *where* the collecting happens into a rule about how much is
    # written around it.
    settle = EDITOR_SOURCE.split("f'Settled with {label}'", 1)[1].split(
        'def ', 1)[0]
    assert '_keep_a_point(' in settle, 'a release does not reach the search'
    press = EDITOR_SOURCE.split("state['gfn_energy_unit']", 1)[1].split(
        'results.append', 1)[0]
    assert '_keep_a_point(' in press and '_settle_price(' in press
    assert 'if not _stopped():' in press, (
        'a stopped run keeps the frame on screen, and the energy in hand is '
        'about the whole run rather than that frame')


# ---------------------------------------------------------------------------
# All of them, not only the lowest
# ---------------------------------------------------------------------------
#
# "vielleicht kann man nicht nur best so far gehen sondern wir merken uns alle
# und man hat da eine liste ... sind dort quasi sortiert", and then the line
# that makes it a set rather than a recording: "also nur minimum TS nicht
# zwischendinger also nur extrema auf der PES".


def test_the_same_minimum_reached_again_is_one_entry():
    """A search by hand falls back into the same minimum again and again --
    that is what searching is -- and a list where eighteen of twenty entries
    are three conformers is worse than no list."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part, walk=((1.20, -100.0000), (1.15, -100.0050),
                          (1.20, -100.0000), (1.30, -100.0010),
                          (1.20, -100.0000)))
    points = part._points_here()
    assert len(points) == 3, [one['energy'] for one in points]
    again = [one for one in points if one['seen'] > 1]
    assert len(again) == 1 and again[0]['seen'] == 3


def test_the_list_is_sorted_and_says_how_far_up_each_one_is():
    """Within one question the energies are one comparison, and the
    difference is the whole of what a conformer set is for."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    assert _shown(part.submit_best_dd)
    labels = [label for label, _value in part.submit_best_dd.options][1:]
    assert labels[0].startswith('minimum, +0.00'), labels
    costs = [float(one.split(', ')[1].split()[0]) for one in labels]
    assert costs == sorted(costs), costs
    assert costs[-1] > 0, costs


def test_a_transition_state_joins_the_list_and_carries_no_number():
    """It is an extremum and belongs beside the minima. Without an energy:
    none is in hand where a search reports what it reached, and pricing it
    there would be a single point nobody asked for -- so it is listed as what
    it is rather than borrowing a number about something else."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    part._keep_a_point(_geo(1.55, 'Optimised to a transition state'), None,
                       'transition state', 'a saddle search')
    part._refresh_the_best()
    labels = [label for label, _value in part.submit_best_dd.options]
    saddles = [one for one in labels if one.startswith('transition state')]
    assert len(saddles) == 1, labels
    assert 'no energy in hand' in saddles[0]
    assert labels[-1] == saddles[0], 'and it sorts after the priced minima'


def test_a_saddle_search_is_what_puts_one_there():
    """One mode going the wrong way, reported by a search, and not a guess:
    a second-order saddle is not a transition state and does not go in."""
    from editor_source import EDITOR_SOURCE

    noted = EDITOR_SOURCE.split('def _note_the_saddle', 1)[1].split(
        'def ', 1)[0]
    assert 'if order == 1:' in noted, noted[-400:]
    assert "'transition state'" in noted


def test_another_question_keeps_its_own_list_beside_this_one():
    """Change the method or the charge and what was found does not become
    wrong -- it becomes the answer to a question you are no longer asking, and
    the neutral's best arrangement is still the neutral's best arrangement
    once you have gone on to the cation."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    part.submit_gfn_charge.value = 1
    part._keep_a_point(_geo(1.25), -99.5000, 'minimum', 'settle')
    part.coords_widget.value = _geo(1.25)
    part._refresh_the_best()
    assert len(part._points_here()) == 1, 'the cation has its own list'

    part.submit_gfn_charge.value = 0
    part.coords_widget.value = _geo(1.30)
    part._refresh_the_best()
    assert len(part._points_here()) == 3, 'and the neutral still has its own'
    labels = [label for label, _value in part.submit_best_dd.options]
    other = [one for one in labels if 'its own best' in one]
    assert len(other) == 1 and 'q +1' in other[0], labels


def test_the_other_questions_carry_no_number():
    """Two totals from two different questions are not comparable, and
    putting them one under the other is an invitation to subtract them."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    part.submit_gfn_charge.value = 1
    part._keep_a_point(_geo(1.25), -99.5000, 'minimum', 'settle')
    part.submit_gfn_charge.value = 0
    part.coords_widget.value = _geo(1.30)
    part._refresh_the_best()
    other = [label for label, _v in part.submit_best_dd.options
             if 'its own best' in label][0]
    assert 'kcal/mol' not in other, other


def test_the_multiplicity_is_part_of_the_question():
    """Per charge the spin state is its own question, which is where the user
    put it: "pro charge ist M interessant, das zaehlt da am besten mit rein"."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    part.submit_gfn_mult.value = 3
    assert part._points_here() == [], 'a triplet is not the singlet"s search'
    part.submit_gfn_mult.value = 1
    assert len(part._points_here()) == 3


def test_going_to_another_question_changes_nothing_above():
    """Bringing a structure back is one thing; switching the method, the
    charge or the spin under the user is another, and doing the second
    silently in order to do the first would answer a question nobody asked."""
    part, _state = _an_editor(_geo(1.20, 'start'))
    _a_search(part)
    part.submit_gfn_charge.value = 1
    part._keep_a_point(_geo(1.25), -99.5000, 'minimum', 'settle')
    part.submit_gfn_charge.value = 0
    part.coords_widget.value = _geo(1.30)
    part._refresh_the_best()

    which = [value for label, value in part.submit_best_dd.options
             if 'its own best' in label][0]
    part.submit_best_dd.value = which
    assert part.submit_gfn_charge.value == 0, 'the charge was changed under us'
    said = ' '.join(part.state.get('mol_status_lines') or ())
    assert 'not comparable' in said, said
    assert (part._geometry_key(part.coords_widget.value)
            == part._geometry_key(_geo(1.25)))


def test_the_list_is_bounded_and_says_when_it_drops_one():
    """A session is not bounded and this list must be. What matters more than
    the number is that dropping one is recorded rather than done quietly."""
    part, state = _an_editor(_geo(1.20, 'start'))
    for n in range(part._POINTS_KEPT + 5):
        part._keep_a_point(_geo(1.10 + 0.05 * n), -100.0 + 0.01 * n,
                           'minimum', 'settle')
    part.coords_widget.value = _geo(1.10)
    part._refresh_the_best()
    assert len(part._points_here()) == part._POINTS_KEPT
    book = state['points_by'][part._best_conditions()]
    assert book['dropped'] == 5, book['dropped']
