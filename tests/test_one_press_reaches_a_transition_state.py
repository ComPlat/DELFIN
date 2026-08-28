"""Three buttons that reached a transition state, made into one press.

The user, on the row this editor's toolbar had grown into: "die Devise ist oft
weniger ist mehr bei Knoepfen ... wir wollen dem User die gesamte
Funktionalitaet geben, aber ihn nicht ueberfordern" -- and then "Knoepfe
zusammenlegen, semantisch richtig, und dass es bedienerfreundlich und dynamisch
ist".

``To the saddle`` climbed whatever was in the box through ORCA, ``Path to
saddle`` walked between two ends and climbed the estimate, and ``Climb to TS``
did the first of those by hand.  They are not three tools: they are three of
the six answers to two questions -- where the search starts, and how it gets
there -- and the difference between them lived only in the tooltips.  ``Find
the path`` was a fourth button differing from ``Path to saddle`` by one step,
which is a setting.

What is asserted here is the thing the user asked for and nothing narrower:
every capability that existed before is still reachable, the combinations that
no button covered are reachable too, and no box is on screen while it names
only one thing.  Driven on the real editor rather than read out of the source,
because which control is on screen when is a fact about the widgets and the
state they answer to.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest


_ETHANE = """8
ethane
C            0.00000000000000        0.00000000000000        0.00000000000000
C            1.53000000000000        0.00000000000000        0.00000000000000
H           -0.36000000000000        1.02000000000000        0.00000000000000
H           -0.36000000000000       -0.51000000000000        0.88400000000000
H           -0.36000000000000       -0.51000000000000       -0.88400000000000
H            1.89000000000000        0.51000000000000        0.88400000000000
H            1.89000000000000        0.51000000000000       -0.88400000000000
H            1.89000000000000       -1.02000000000000        0.00000000000000
"""
_STRETCHED = _ETHANE.replace('1.53000000000000', '2.53000000000000')
_FURTHER = _ETHANE.replace('1.53000000000000', '3.53000000000000')


def _an_editor(text=_ETHANE):
    """One structure editor over a coordinate box of its own.

    The real part, driven the way a user drives it.  This class of question --
    which control is on screen, and what the press does with the boxes beside
    it set that way -- cannot be answered by reading the source, because the
    source says what it means to do.
    """
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
    part.submit_manip_toolbar.layout.display = 'flex'
    return part, state


def _method(part, prefix):
    for _label, value in part.submit_ff_dd.options:
        if str(value).lower().startswith(prefix):
            return value
    return None


def _shown(widget):
    return (widget.layout.display or '') != 'none'


def _values(box):
    return [value for _label, value in box.options]


def test_the_two_questions_are_asked_once_each_and_reach_every_combination():
    """Where the search starts, and how it gets there: two boxes, one press.

    The three buttons covered four of the six combinations and could not say
    which of two pairs of ends ``Path to saddle`` would use -- it preferred
    the marked pair silently, so a scan walked after marking something could
    not be walked between at all.  Asked as two questions, every combination
    is reachable and each is reachable by saying what you want rather than by
    knowing which button means which.

    Nine of them now rather than the six the three buttons could not cover:
    three starts, with one way up from the structure on screen and four from a
    pair of ends.  The point is that the count is free -- a way added to the
    second box is reachable from every start it works from, without a button.
    """
    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')

    # Nothing marked and no scan walked: no pair, so no press. Searching
    # between two ends needs two ends, and what is on screen is the other
    # press's business -- offered where the structure is on its way up rather
    # than from any minimum at all.
    assert not _shown(part.submit_saddle_btn)
    assert not _shown(part.submit_optts_btn)
    part.submit_climb_btn.value = True
    part._refresh_saddle_controls()
    assert _shown(part.submit_optts_btn)
    assert not _shown(part.submit_saddle_btn), 'there is still no pair'
    part.submit_climb_btn.value = False
    part._refresh_saddle_controls()

    # A scan leaves two ends, and the second start appears with them.
    state['scan_ends'] = (_ETHANE, _STRETCHED)
    part._refresh_saddle_controls()
    assert _values(part.submit_saddle_from) == ['scan']
    assert _shown(part.submit_saddle_btn), (
        'one start is not a choice, so what appears is the press')

    # And a marked end is the third, named apart from the scan's -- which is
    # the combination no button reached.
    state['path_from'] = _FURTHER
    part._refresh_saddle_controls()
    assert _values(part.submit_saddle_from) == ['marked', 'scan']

    # What is on screen is not on this list at all: it is the other press,
    # which has one way and needs no box to say so.

    # From two ends there are four, and the walk alone is one of them.  The
    # order of the list is the recommendation: the interpolated search first,
    # the band after it for when that answer is not believed, then the two
    # that run here.
    for start in ('marked', 'scan'):
        part.submit_saddle_from.value = start
        assert _values(part.submit_saddle_how) == ['orca', 'neb', 'hand',
                                                   'walk']
        assert _shown(part.submit_saddle_how)


def test_the_marked_pair_appears_when_the_second_structure_is_drawn():
    """A mark and the structure on screen are a pair only once they differ.

    Pressing the mark on what is already on screen marks one structure, not
    two, so the start that names the pair cannot appear then -- it has to
    appear when the other one is loaded or built.  Asked where every host goes
    through when what is on screen changes, or the option would wait for the
    next change of method.
    """
    from editor_source import EDITOR_SOURCE

    drawing = EDITOR_SOURCE.split('def _replace_mol_output_view(', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_refresh_saddle_controls()' in drawing

    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')
    state['current_xyz_for_copy'] = {'content': _ETHANE}
    part.on_submit_path_from()
    # The same structure twice is not a pair, so nothing is offered yet.
    part._refresh_saddle_controls()
    assert _values(part.submit_saddle_from) == ['here']
    state['current_xyz_for_copy'] = {'content': _STRETCHED}
    part._refresh_saddle_controls()
    assert _values(part.submit_saddle_from) == ['marked']


def test_the_press_says_which_of_the_two_things_it_is_about_to_do():
    """"To the saddle" over a walk that stops at an estimate is a promise the
    press does not keep, so the name follows the box beside it."""
    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')
    state['scan_ends'] = (_ETHANE, _STRETCHED)
    part._refresh_saddle_controls()
    part.submit_saddle_from.value = 'scan'

    part.submit_saddle_how.value = 'orca'
    assert part.submit_saddle_btn.description == 'To the saddle'
    part.submit_saddle_how.value = 'hand'
    assert part.submit_saddle_btn.description == 'To the saddle'
    part.submit_saddle_how.value = 'walk'
    assert part.submit_saddle_btn.description == 'Find the path'


def test_the_press_runs_what_the_two_boxes_say(monkeypatch):
    """Every combination, driven through the press, checked at the engine.

    The engines are stood in for -- what is being measured is which of them
    the press reaches and with which pair of structures, and that is decided
    before any of them is called.
    """
    part, state = _an_editor()
    from delfin.dashboard import gfn_optimize as gfn_mod
    from delfin.dashboard import saddle as saddle_mod

    called = []

    def _no_saddle(xyz_text, method='gfn2', **kw):
        called.append(('orca-here', xyz_text))
        return {'ok': False, 'status': 'stood in for'}

    def _no_chain(reactant, product, method='gfn2', **kw):
        called.append(('orca-chain', reactant, product))
        return {'ok': False, 'status': 'stood in for'}

    def _no_walk(reactant, product, method='gfn2', **kw):
        called.append(('walk', reactant, product))
        return {'ok': False, 'status': 'stood in for'}

    monkeypatch.setattr(saddle_mod, 'optimise_to_saddle', _no_saddle)
    monkeypatch.setattr(saddle_mod, 'path_to_saddle', _no_chain)
    monkeypatch.setattr(gfn_mod, 'walk_the_path', _no_walk)

    part.submit_ff_dd.value = _method(part, 'gfn2')
    state['current_xyz_for_copy'] = {'content': _ETHANE}
    state['scan_ends'] = (_ETHANE, _STRETCHED)
    state['path_from'] = _FURTHER
    part._refresh_saddle_controls()

    wanted = [
        ('marked', 'orca', ('orca-chain', _FURTHER, _ETHANE)),
        ('scan', 'orca', ('orca-chain', _ETHANE, _STRETCHED)),
        ('marked', 'walk', ('walk', _FURTHER, _ETHANE)),
        ('scan', 'walk', ('walk', _ETHANE, _STRETCHED)),
        # By hand from two ends is the walk first and then the climb on what
        # it estimated: the combination no button covered.
        ('marked', 'hand', ('walk', _FURTHER, _ETHANE)),
    ]
    for start, how, expected in wanted:
        called.clear()
        part.submit_saddle_from.value = start
        part.submit_saddle_how.value = how
        part.on_submit_saddle()
        _wait_for(called)
        assert called and called[0] == expected, (start, how, called)
        # Whatever the stand-in left behind, so the next press is not refused
        # for a run that is still going.
        for key in ('saddle_run', 'chain_run', 'path_run'):
            state.pop(key, None)

    # And the press that takes no pair, which runs what start-from-here ran:
    # one press per operation, and the machinery underneath unchanged.
    called.clear()
    part.on_submit_optts()
    _wait_for(called)
    assert called and called[0] == ('orca-here', _ETHANE), called


def _wait_for(called, seconds=10.0):
    """The runs go on a thread; the stand-ins return at once."""
    import time
    began = time.time()
    while not called and time.time() - began < seconds:
        time.sleep(0.01)


def test_a_method_reached_through_another_program_offers_both_starts():
    """What a box hides must be what cannot run, never what is merely deeper.

    g-xTB is a build of its own and ORCA has no keyword for it, so for a while
    there was no way up from the structure on screen at all -- and left to the
    two boxes to sort out between them, the start stood at "what is on screen",
    found no way, hid the press, and hid with it the only control that could
    have chosen the start that worked.

    ORCA publishes ExtOpt, which is its interface to a program it does not
    know, so the saddle optimiser reaches g-xTB after all and both starts are
    real again.  The interesting half of this test is therefore not the answer
    but where it comes from: the box is filled from the same table the run
    reads, so a route that appears in one appears in the other without anybody
    remembering to add it twice.
    """
    part, state = _an_editor()
    gxtb = _method(part, 'gxtb')
    if gxtb is None:
        pytest.skip('this build offers no g-xTB')
    state['scan_ends'] = (_ETHANE, _STRETCHED)
    part.submit_ff_dd.value = gxtb

    from delfin.dashboard import climb as _climb
    from delfin.dashboard import saddle as _saddle
    assert 'gxtb' in _saddle.SADDLE_METHODS      # through ExtOpt
    assert 'gxtb' not in _climb.CLIMB_METHODS    # a process per gradient
    assert _values(part.submit_saddle_from) == ['scan']
    assert _shown(part.submit_saddle_btn)
    # By hand is the one way that is still missing, and it is missing for a
    # measured reason rather than an oversight: a g-xTB gradient is a whole
    # process, 0.29 s at sixteen atoms against 6 ms, and one exact Hessian is
    # 17.9 s against 0.55.  So the switch that walks a release uphill stays
    # off rather than refusing after it has been pressed.
    assert 'hand' not in _values(part.submit_saddle_how)
    assert not _shown(part.submit_climb_btn)
    # And the mark stays, because it describes two structures, not a program.
    assert _shown(part.submit_path_from_btn)


def test_a_choice_that_cannot_be_met_is_kept_rather_than_overwritten():
    """A method that offers less must not silently reset what was chosen.

    Both lists are decided by the method and by which structures exist, and
    both of those change while the user is doing something else.  Handed a
    list its value is not in, a dropdown loses that value -- which is how
    "the path only" came back as "through ORCA" after a detour through a
    method that could run neither.
    """
    part, state = _an_editor()
    state['scan_ends'] = (_ETHANE, _STRETCHED)
    part.submit_ff_dd.value = _method(part, 'gfn2')
    part.submit_saddle_from.value = 'scan'
    part.submit_saddle_how.value = 'hand'

    part.submit_ff_dd.value = _method(part, 'uff')
    assert not _shown(part.submit_saddle_btn)

    part.submit_ff_dd.value = _method(part, 'gfn2')
    assert part.submit_saddle_from.value == 'scan'
    assert part.submit_saddle_how.value == 'hand'


def test_the_climb_stands_with_the_switches_and_not_with_the_press():
    """It is a mode, not a third way of pressing for a transition state.

    Climb to TS says which way the optimiser walks when an atom is let go --
    the same release path as Dynamik Opt and Auto with the optimiser going
    uphill instead of down (e3442010).  Beside the saddle press it read as one
    of three buttons that all reached a transition state, which is exactly the
    confusion this row was reported for.
    """
    part, _state = _an_editor()
    children = list(part.submit_manip_toolbar.children)
    switches = {part.submit_optimize_btn, part.submit_optimize_all_btn,
                part.submit_relax_btn, part.submit_auto_btn,
                part.submit_settle_btn}

    # It is said as a direction now rather than as a switch of its own: two
    # switches are four combinations of which two mean nothing, and climbing
    # with the field off is not a thing that happens. The box stands where
    # the switch stood, among the ones that say what a walk does.
    where = children.index(part.submit_climb_way)
    assert children[where - 1] in switches, type(children[where - 1]).__name__
    assert children[where + 1] is part.submit_climb_btn

    # And the switch itself is never on the row: one question, one control.
    assert part.submit_climb_btn.layout.display == 'none'

    # It is still a switch underneath, not a press: it stays down across the
    # walks it starts, which is what makes it a mode -- and everything about a
    # climb still hangs off it, which is why the box drives it rather than
    # replacing it.
    import ipywidgets as widgets
    assert isinstance(part.submit_climb_btn, widgets.ToggleButton)


def test_an_armed_scan_comes_back_with_the_method_that_can_walk_it():
    """"The armed legs are kept" has to mean they can be reached again.

    A method that cannot walk a scan takes the row off the toolbar and says
    the legs are kept for when a GFN method comes back.  They were -- and
    nothing put the row back, so a detour through UFF left a scan still armed
    with no list, no Run scan and no way to reach either.  Kept and
    unreachable is the same as gone.
    """
    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')
    state['scan_legs'] = [{'kind': 'distance', 'atoms': [0, 1], 'from': 1.53,
                           'to': 3.0, 'steps': 20, 'structure': 'x'}]
    part._refresh_scan()
    assert _shown(part.submit_scan_run_btn)

    part.submit_ff_dd.value = _method(part, 'uff')
    assert not _shown(part.submit_scan_run_btn)

    part.submit_ff_dd.value = _method(part, 'gfn2')
    assert _shown(part.submit_scan_run_btn)
    assert _shown(part.submit_scan_dd)
    assert _shown(part.submit_scan_whole)


def test_where_a_scan_walks_is_one_question_and_what_is_picked_answers_it():
    """A direction and a checkbox that only revealed a field were two controls
    for one question.

    A value says which way the walk goes all by itself, so "further apart, to
    2.40" was the same fact twice with nothing checking that the two halves
    agreed.  The field for the number appears under the third answer and
    nowhere else, and it opens on what the coordinate measures rather than on
    a zero that has to be guessed at.

    How many answers there are is what is picked: a pair of atoms can be asked
    to make or break the bond between them, and three or four atoms cannot --
    there is no bond between three atoms to form.  That absence is the list
    saying so, rather than a press refused afterwards.
    """
    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')
    assert not hasattr(part, 'submit_scan_stop_at')

    state['picked'] = [0, 1]
    part.submit_internal_value.value = 1.53
    part._refresh_scan()
    assert _values(part.submit_scan_way) == ['in', 'out', 'form', 'break',
                                             'to']
    assert not _shown(part.submit_scan_to)

    part.submit_scan_way.value = 'to'
    assert _shown(part.submit_scan_to)
    assert part.submit_scan_to.value == pytest.approx(1.53)

    # And the words follow what is picked -- three answers for a torsion, the
    # two verbs gone with the bond they were about.
    state['picked'] = [0, 1, 2, 3]
    part._refresh_scan()
    assert _values(part.submit_scan_way) == ['in', 'out', 'to']
    assert part.submit_scan_way.value == 'to', 'the answer given was lost'
    assert [label for label, _value in part.submit_scan_way.options][:2] == \
        ['narrower', 'wider']


def test_a_walk_to_where_the_coordinate_already_is_is_refused_rather_than_guessed():
    """A value implies its own direction, and there is none left to fall back
    on when the value is the one the coordinate has.

    The fallback used to be inwards, silently, because the direction was a
    separate box that still held an answer.  A scan that walks the opposite
    way from the one the number implied is worse than no scan.
    """
    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')
    state['picked'] = [0, 1]
    part.submit_internal_value.value = 1.53
    part._refresh_scan()
    part.submit_scan_way.value = 'to'
    part.submit_scan_to.value = 1.53

    part.on_submit_scan()
    assert not state.get('scan_legs')
    assert 'already is' in ' '.join(state.get('mol_status_lines') or ())
