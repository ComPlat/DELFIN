"""Which way the live relaxation walks, said once instead of twice.

Dynamik Opt and Climb to TS were two switches, and two switches are four
combinations of which two mean nothing: climbing with the field off is not a
thing that happens.  It is one question with two answers -- and the toolbar's
own comment had already said so, that Climb to TS "is the same release path
with the optimiser walking uphill instead of down".

Saying it as a direction is also what makes the other press explainable.  With
two switches, "Climb to TS" and "To the saddle" read as two ways of pressing
for a transition state and nothing said which to use.  As a direction and a
press they are plainly different things: one is where your hand walks, the
other is DELFIN going to look on its own.

What is not changed is the machinery.  Everything about a climb hangs off the
switch -- the release path, the refusals, the page's own reading of which
switches are pressed -- and the box drives that switch rather than replacing
it.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest


def _an_editor():
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
    part = structure_editor.build(
        ctx, state={},
        coords_widget=widgets.Textarea(
            value='3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    part._set_manip_toolbar_enabled(True)
    return part


def _method(part, name):
    part.submit_ff_dd.value = [value for _label, value in
                               part.submit_ff_dd.options
                               if name in str(value).lower()][0]


def _shown(widget):
    return (widget.layout.display or '') != 'none'


def test_it_starts_downhill():
    """A relaxation goes downhill. Climbing is the thing you ask for."""
    part = _an_editor()
    assert part.submit_climb_way.value == 'down'
    assert part.submit_climb_btn.value is False


def test_the_box_is_what_is_seen_and_the_switch_is_what_acts():
    part = _an_editor()
    _method(part, 'gfn2')
    assert _shown(part.submit_climb_way)
    assert not _shown(part.submit_climb_btn), (
        'two controls for one question is what this removes')

    part.submit_climb_way.value = 'up'
    assert part.submit_climb_btn.value is True
    part.submit_climb_way.value = 'down'
    assert part.submit_climb_btn.value is False


def test_the_box_follows_whatever_moved_the_switch():
    """A refusal that puts the switch back down has to be able to say so, or
    the row reads "walks uphill" over a release that goes down."""
    part = _an_editor()
    _method(part, 'gfn2')
    part.submit_climb_btn.value = True
    assert part.submit_climb_way.value == 'up'
    part.submit_climb_btn.value = False
    assert part.submit_climb_way.value == 'down'


def test_a_method_that_cannot_climb_offers_no_direction():
    """Left visible under a method it cannot ask for a gradient from, it
    refused only after the press -- a control that promises what it cannot
    do, under the most accurate method in the list."""
    part = _an_editor()
    _method(part, 'gfn2')
    part.submit_climb_way.value = 'up'
    assert _shown(part.submit_climb_way)

    _method(part, 'uff')
    assert not _shown(part.submit_climb_way)
    assert part.submit_climb_btn.value is False, (
        'a direction that survived the change would be a release walking a '
        'way the method has no answer for')
    assert part.submit_climb_way.value == 'down'


def test_it_stands_with_the_switches_that_say_what_a_walk_does():
    part = _an_editor()
    row = list(part.submit_manip_toolbar.children)
    assert row.index(part.submit_climb_way) == row.index(
        part.submit_auto_btn) + 1
    assert row.index(part.submit_relax_btn) < row.index(part.submit_climb_way)


# ---------------------------------------------------------------------------
# The hand the page thinks it has
# ---------------------------------------------------------------------------


def _watching():
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
    part = structure_editor.build(
        ctx, state={},
        coords_widget=widgets.Textarea(
            value='3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    part._set_manip_toolbar_enabled(True)
    return part, said


def _shares(said):
    import re
    return [m.group(1) for one in said
            for m in re.finditer(r'setPullStrength\("[^"]+",([^,]+),', one)]


def test_the_page_is_told_what_the_hand_is_worth_before_the_first_drag():
    """The page decides between a force and a placement on one number, and at
    a share of nothing it skips the pull entirely and sets the coordinate --
    which is the placing hand.  Nothing is what it holds until it is told.

    It was told only when a control moved, and under a browser method also
    when the field was handed its parameters.  Under a server method neither
    happens on the way in: choosing GFN2 sent nothing, switching Dynamik Opt
    on sent nothing, and the first drag was a placement however the box read
    -- "es ist jetzt immer move the atom egal was ich hier rein mache".
    """
    part, said = _watching()
    gfn = [value for _label, value in part.submit_ff_dd.options
           if 'gfn2' in str(value).lower()][0]

    said.clear()
    part.submit_ff_dd.value = gfn
    assert _shares(said), 'the method was chosen and the page was told nothing'
    assert _shares(said)[-1].startswith('0.4'), _shares(said)

    said.clear()
    part.submit_relax_btn.value = True
    assert _shares(said), 'dragging became possible and the page was told nothing'


def test_the_rule_that_says_it_tells_it():
    """`_refresh_hand_controls` describes itself as telling the page the same
    number through setPullStrength. It now does."""
    part, said = _watching()
    said.clear()
    part._refresh_hand_controls()
    assert _shares(said), 'the rule says it tells the page, and it must'


def test_the_bug_button_stays_at_the_right_hand_end():
    """It was pushed there by the status line taking the space between. That
    line lies on the picture now, so the push has to be its own."""
    part, _said = _watching()
    assert part.submit_bug_group.layout.margin == '0 0 0 auto'
    assert list(part.submit_manip_toolbar.children)[-1] is part.submit_bug_group
