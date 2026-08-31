"""A finished scan announces what it made possible, and keeps its answer.

Two things the user found in one sitting, and they are the same moment seen
from two sides -- the end of a scan.

**The announcement.**  Two buttons used to appear when a scan left two ends:
``Find the path`` and ``Path to saddle``.  Their arriving was the editor
saying that something new could be done now.  085cfd55 folded them into one
press and two boxes, which is the shape that was asked for -- "weniger ist
mehr bei Knoepfen" -- but the moment went with them: the box gained an entry
and went on standing at "what is on screen", so the press went on meaning
climb what is in the box.  The user: "nach einem Scan konnten wir doch davor
noch weiter den Path absuchen mit Buttons -- wo sind die hin, geht das nicht
mehr?"

Measured on the real editor before the change, GFN2 with a scan's two ends in
hand::

    starts ['here', 'scan'] standing at here   ways ['orca'] (box hidden)
    one press: climb what is on screen, through ORCA

and after::

    starts ['here', 'scan'] standing at scan   ways ['orca', 'neb', 'hand',
                                                       'walk']
    one press: walk the scan's two ends, then climb through ORCA

Nothing left the matrix: the combinations the two boxes can be put into
between them are the same ones before and after -- nine under GFN-FF and GFN2
with a marked end and a scan's ends in hand, seven under g-xTB, which has no
by-hand climb -- and "what is on screen" is one selection away.

(Seven and five when this was written.  A nudged elastic band arrived later
as a fourth *how* from a two-ended start, which is what the shape of this
control is for: another way of answering the second question is an entry in
the box, not a fourth press beside it.)

**The answer.**  The user, in the same message: "und wenn ich Scan mache, sehe
ich am Ende nicht den Text, er verschwindet und ich sehe ``Optimised with
GFN2-xTB. E = -23.177979 Eh ...``".  The scan's verdict -- the barrier, where
it ended, the temperature it wants -- was written and then replaced by a
sentence about an optimisation.  No optimisation was running.  The page
reports on the frames it draws, a scan streams a frame a point, and that
report rewrote the row with whatever the *optimiser* had last said, however
long ago that was.  Reproduced on the real tab: one press of Optimize, a scan
of butane's central C-C under GFN2 in benzene, then the single message the
browser sends when it takes the scan's frames::

    after the scan   The scan walked 8 of 8 points. Highest +52.8 kcal/mol ...
    after the page   Optimised with GFN2-xTB. E = -13.670432 Eh ...

A result that took minutes to compute must not be erased by a message that
carries no news, so what the page has to say is put on the end of the row at
the moment it is drawn -- the way the spinner is -- and can no longer take the
row away from whoever wrote it.
"""

from __future__ import annotations

import pathlib
import shutil
import tempfile
import time

import pytest


_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'),
                                reason='xtb not installed')

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

#: Water: another molecule, for the ends that must stop being offered.
_WATER = """3
water
O            0.00000000000000        0.00000000000000        0.00000000000000
H            0.75800000000000        0.58600000000000        0.00000000000000
H           -0.75800000000000        0.58600000000000        0.00000000000000
"""

#: Butane, walked apart at the central C-C.  A distance rather than a torsion
#: because the leg starts from the value read off the structure at the press,
#: so the number in the box is only the target and the test does not depend on
#: the browser having measured anything.
_BUTANE = """14
butane
C -1.9026 -0.2415  0.0000
C -0.6178  0.5776  0.0000
C  0.6178 -0.5776  0.0000
C  1.9026  0.2415  0.0000
H -1.9741 -0.8776  0.8834
H -1.9741 -0.8776 -0.8834
H -2.7756  0.4128  0.0000
H -0.5892  1.2266  0.8794
H -0.5892  1.2266 -0.8794
H  0.5892 -1.2266  0.8794
H  0.5892 -1.2266 -0.8794
H  1.9741  0.8776  0.8834
H  1.9741  0.8776 -0.8834
H  2.7756 -0.4128  0.0000
"""


def _an_editor(text=_ETHANE):
    """One structure editor over a coordinate box of its own.

    The real part, driven the way a user drives it: which control is on screen
    and what the press means with the boxes standing where they stand are
    facts about the widgets, not about what the source says it intends.
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


@pytest.fixture
def tab(tmp_path):
    """A whole Submit tab, which is where the scan's own path runs."""
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda _script: None
    _widget, refs = tab_submit.create_tab(ctx)
    refs['coords_widget'].value = _BUTANE
    return refs


def _shown(widget):
    return (widget.layout.display or '') != 'none'


def _values(box):
    return [str(value) for _label, value in box.options]


def _method(part, prefix):
    for _label, value in part.submit_ff_dd.options:
        if str(value).lower().startswith(prefix):
            return value
    return None


def _reachable(part):
    """Every (start, way) the two boxes can be put into between them."""
    was = str(part.submit_saddle_from.value)
    out = set()
    for start in _values(part.submit_saddle_from):
        part.submit_saddle_from.value = start
        for way in _values(part.submit_saddle_how):
            out.add((start, way))
    # The box may no longer offer what it held: "what is on screen" left this
    # press for the one beside it, and putting it back would be asking for a
    # start that does not exist.
    if was in _values(part.submit_saddle_from):
        part.submit_saddle_from.value = was
    return out


def _run_the_scan(refs, steps=8):
    """Butane's central C-C, walked from where it is out to 2.20 A."""
    refs['submit_ff_dd'].value = 'gfnff'
    refs['submit_scan_how'].value = 'hold'
    refs['submit_pick_sync'].value = '1,2'
    refs['submit_internal_value'].value = 1.53
    refs['submit_scan_way'].value = 'to'
    refs['submit_scan_to'].value = 2.20
    refs['submit_scan_steps'].value = steps
    refs['submit_scan_whole'].value = True
    refs['submit_scan_add_btn'].click()
    refs['submit_scan_run_btn'].click()
    began = time.time()
    while refs['editor_state'].get('scan_run') and time.time() - began < 240:
        time.sleep(0.05)
    time.sleep(0.3)


@_needs_xtb
def test_a_finished_scan_leaves_the_next_thing_one_press_away(tab):
    """After a real scan, the press means the thing the scan made possible.

    The acceptance the user stated: the next thing he would want is one press
    away and he can see that it is.  What he can see is the start naming the
    scan's two ends and the second box arriving beside it -- which is the same
    news the two buttons used to carry by appearing.
    """
    refs = tab
    assert _values(refs['submit_saddle_from']) == ['here']
    assert not _shown(refs['submit_saddle_from'])
    assert not _shown(refs['submit_saddle_how'])

    _run_the_scan(refs)
    assert refs['editor_state'].get('scan_ends'), 'the scan left no two ends'

    # The start stands on the pair, and the press starts from it.
    assert refs['submit_saddle_from'].value == 'scan'
    # What is on screen has left this press for the one beside it, so the
    # pair is the whole list -- and a list of one is not shown.
    assert _values(refs['submit_saddle_from']) == ['scan']
    assert _shown(refs['submit_saddle_btn'])
    # And the second box has arrived: from a pair there are three ways, and
    # from the structure on screen there was only ever one, so this box was
    # not on the toolbar at all a moment ago.
    # A band arrived here later, and is another way from the same pair --
    # so the box that was not on the toolbar a moment ago now holds four.
    assert _values(refs['submit_saddle_how']) == ['orca', 'neb', 'hand',
                                                  'walk']
    assert _shown(refs['submit_saddle_how'])

    # What is on screen is not lost either -- it is the press beside this
    # one, which converges a guess rather than searching between two ends and
    # so has one way and needs no box to say so.
    assert 'here' not in _values(refs['submit_saddle_from'])


@_needs_xtb
def test_the_verdict_says_what_the_toolbar_has_just_done(tab):
    """The closing sentence describes the toolbar as it now stands.

    It said the box "will now start from them" while the box stood on the
    structure on screen, which is a description of something the user would
    have had to go and find.  Now it is true when it is read.
    """
    refs = tab
    _run_the_scan(refs)
    said = ' '.join(refs['editor_state'].get('mol_status_lines') or ())
    assert 'The scan walked' in said, said
    assert 'It left two ends' in said, said
    assert 'saddle press now starts from them' in said, said
    assert refs['submit_saddle_from'].value == 'scan'


@_needs_xtb
def test_a_real_scans_verdict_is_still_there_after_the_page_reports(tab):
    """The user's own sequence, end to end and on the real tab.

    A scan, and then the one message the browser sends when it takes the
    frames the scan streamed.  The sentence it used to be replaced by is put
    in place first, because that is what a press of Optimize leaves behind and
    the whole point is that a run which finished before the scan could reach
    forward and take its answer.
    """
    refs = tab
    refs['editor_state']['gfn_last_status'] = (
        'Optimised with GFN2-xTB. E = -23.177979 Eh (-14544.40 kcal/mol). '
        'In benzene (ALPB). charge 0, multiplicity 1.')
    _run_the_scan(refs)
    verdict = tuple(refs['editor_state'].get('mol_status_lines') or ())
    assert verdict and verdict[0].startswith('The scan walked'), verdict

    refs['submit_cmd_sync'].value = 'gfnplay:9:received 8 frames'
    assert tuple(refs['editor_state'].get('mol_status_lines') or ()) == verdict
    assert 'Optimised with GFN2-xTB' not in refs['mol_status'].value


def test_the_move_costs_the_matrix_nothing():
    """Every combination that could be reached before can still be reached.

    Six answers to two questions, and the change moves where the boxes rest
    rather than what they hold.  Driven over each engine that can reach a
    saddle at all, with a marked end and a scan's ends in hand::

        method  starts                     ways                      resting
        gfnff   here, marked, scan         orca, neb, hand, walk     scan
        gfn2    here, marked, scan         orca, neb, hand, walk     scan
        gxtb    here, marked, scan         orca, neb, walk           scan

    g-xTB has no by-hand climb -- there is no gradient the climb knows how to
    ask it for -- which is why its row is one shorter, and that was true
    before the move as well.
    """
    for prefix, ways in (('gfnff', {'orca', 'neb', 'hand', 'walk'}),
                         ('gfn2', {'orca', 'neb', 'hand', 'walk'}),
                         ('gxtb', {'orca', 'neb', 'walk'})):
        part, state = _an_editor()
        method = _method(part, prefix)
        if method is None:
            continue
        part.submit_ff_dd.value = method
        # A pair is two marks: the beginning and the end are separate
        # slots, so marking only one offers nothing to search between.
        state['path_from'] = _ETHANE
        state['path_to'] = _FURTHER
        state['scan_ends'] = (_ETHANE, _STRETCHED)
        part._refresh_saddle_controls()
        before = _reachable(part)

        part._scan_left_two_ends()
        after = _reachable(part)

        assert after == before, prefix
        assert before == {(start, way)
                          for start in ('marked', 'scan')
                          for way in ways}, prefix
        # And the resting place is the pair the scan left.
        assert part.submit_saddle_from.value == 'scan', prefix


def test_two_ends_from_another_molecule_are_not_offered():
    """They name atoms, and atoms belong to the structure they were walked on.

    The ends outlive the structure: a scan leaves them and loading something
    else does not take them away.  Offered against a different molecule, the
    walk between them would put the previous structure into the box under the
    name of a path -- and the move above would make that the resting answer
    rather than one the user had to choose.
    """
    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')
    state['scan_ends'] = (_ETHANE, _STRETCHED)
    part._refresh_saddle_controls()
    assert 'scan' in _values(part.submit_saddle_from)

    # Another molecule in the box, the ends still in hand.
    part.coords_widget.value = _WATER
    part._refresh_saddle_controls()
    assert _values(part.submit_saddle_from) == ['here']
    # And the press refuses them rather than walking them.
    assert part._path_ends('scan') is None
    assert 'on this structure' in ' '.join(state.get('mol_status_lines') or ())


def test_the_page_reporting_on_frames_cannot_take_the_row():
    """A message that carries no news must not replace one that does.

    The page reports on the frames it draws, a scan streams a frame a point,
    and that report used to rewrite the row with whatever the optimiser had
    last said.  Nothing had to be running for it to happen: the sentence is
    remembered, so a run that ended before the scan started could take the
    scan's verdict away minutes later.
    """
    part, state = _an_editor()
    # What an optimisation left behind, exactly as the user quoted it.
    state['gfn_last_status'] = (
        'Optimised with GFN2-xTB. E = -23.177979 Eh (-14544.40 kcal/mol). '
        'In benzene (ALPB). charge 0, multiplicity 1.')
    verdict = ('The scan walked 8 of 8 points. Highest +52.8 kcal/mol at 2.2.',
               'It wants about 691 K to be crossed within an hour.')
    part._set_mol_status(*verdict)
    was = part.mol_status.value

    part.submit_cmd_sync.value = 'gfnplay:9:received 8 frames'
    assert tuple(state.get('mol_status_lines') or ()) == verdict
    assert part.mol_status.value == was
    assert 'Optimised with GFN2-xTB' not in part.mol_status.value


def test_a_fault_still_reaches_the_row_and_does_not_pile_up():
    """The one thing the page has to say goes on the end of what is there.

    "received 41 frames" is the case where nothing is wrong and it was never
    a row; a playback that cannot draw is the case the reporting was added
    for.  Put on at the draw rather than written into the message, so it
    decorates whatever the row says and cannot be said twice.
    """
    part, state = _an_editor()
    part._set_mol_status('The scan walked 8 of 8 points.')

    part.submit_cmd_sync.value = 'gfnplay:1:setPositions did not draw'
    once = part.mol_status.value
    assert 'setPositions did not draw' in once
    assert once.count('<br>') == 0, 'it goes on the end, not under it'
    assert once.count('Trajectory:') == 1

    part.submit_cmd_sync.value = 'gfnplay:2:setPositions did not draw'
    assert part.mol_status.value == once
    # And what is remembered is the message, not the message plus the note --
    # otherwise the next draw says it again.
    assert tuple(state.get('mol_status_lines') or ()) == (
        'The scan walked 8 of 8 points.',)


def test_marking_a_pair_wins_over_the_ends_a_scan_left():
    """A scan moves the start onto its own two ends, which is right: minutes
    of walking is a strong statement about what the user is interested in.

    Marking is a stronger one.  It is two deliberate presses, one end at a
    time, and it has to win over a pair that was left standing there by
    something else -- otherwise the marks appear in the list, the box goes on
    saying "the scan's two ends", and the press searches between the scan's
    while the line says it searches between the marks.

    Reported: a scan, then To the saddle, then a beginning and an end marked
    by hand -- "hat einfach aus dem scan die enden genommen".
    """
    part, state = _an_editor()
    method = _method(part, 'gfn2')
    if method is None:
        pytest.skip('no method that can reach a saddle')
    part.submit_ff_dd.value = method

    # A finished scan, which moves the start onto its ends.
    state['scan_ends'] = (_ETHANE, _STRETCHED)
    part._scan_left_two_ends()
    assert part.submit_saddle_from.value == 'scan'

    # Now two marks by hand.
    state['path_from'] = _ETHANE
    part.coords_widget.value = _FURTHER
    part.on_submit_path_from(None)

    assert state.get('path_to')
    assert part.submit_saddle_from.value == 'marked', (
        'the press would still search between the ends the scan left')

    # And the scan's pair is still there, one selection away -- nothing was
    # taken from the user, the newest statement simply won.
    assert 'scan' in [value for _label, value in part.submit_saddle_from.options]
