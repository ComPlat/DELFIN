"""What is on the toolbar after a scan is an account of what the scan left.

The user, having restarted the dashboard on f1be8954 and run one: "Habe ich
jetzt nach so einem Scan nicht die Moeglichkeit, dafuer einen TS zu optimieren
oder noch andere Methoden den Path zu untersuchen? Davor hatten wir das doch,
oder?"  He had ``To the saddle`` and the ``Climb to TS`` switch and neither
box beside them -- which is the toolbar saying there is nothing further to do
with this structure.  It was not true.

The rule he stated for the whole toolbar: "Funktionen -- also Buttons --
kommen nur dann, wenn sie etwas machen und bringen, damit man immer weiss, was
man machen kann."  So an absence is itself a statement, and a wrong absence is
the editor claiming something about the world that is false.

**What was measured.**  f1be8954 puts the press on the pair at the end of a
scan, and every scan that finishes its walk does leave one -- driven on the
real editor with real xtb, the whole profile, the walk that stops at the next
minimum, the push with a force, two legs walked together and a walk cut short
with the Stop button all end with ``['here', 'scan']`` standing at ``scan``
and ``['orca', 'neb', 'hand', 'walk']`` beside it.  A walk that ends any
*other* way
did not::

    GFN2, butane's central C-C out to 2.20 A in 8 steps, with an xtb that
    stops converging at step 4 (a real one, driven into an SCF it cannot
    solve):

        from box  ['here']  HIDDEN        how box  ['orca']  HIDDEN
        press     To the saddle
        row       "The scan walked 3 of 8 points. Highest +11.8 kcal/mol at
                   1.88 ... You have 22.3 kcal/mol at 298.15 K, so the whole
                   path is open."

Three points had been walked and their two ends were in hand.  Nothing offered
them, and the sentence saying why the walk had stopped had been written and
then replaced by that verdict.  The cause is one line of control flow: the
pair was assigned at the end of the walk, and the four early exits leave
through the ``finally`` instead, so a walk that gave up never reached it.  The
pair is now assigned in the ``finally``, where every way out goes.

And where there really is no pair -- a walk that priced one point has its
start and its end at the same geometry -- the verdict says so, because silence
is what the user was complaining about::

        "It left no two ends to walk between, so there is no path to
         investigate from it and the box beside the press has no pair to
         offer -- what is on screen is all there is to climb."

**And back.**  The user again: "Undo muss dann natuerlich in Teilen auch
wieder Funktionen nehmen."  A scan makes the pair, so Undo back past the scan
takes it away again -- the entry, the ways from it, and the box's own value --
and Redo puts all three back.  Measured over the three entries a scan leaves,
so the toolbar is true at every station of the walk back and not only at its
ends.
"""

from __future__ import annotations

import pathlib
import shutil
import tempfile
import time

import pytest


_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'),
                                reason='xtb not installed')

#: Butane, walked at its central C-C.  A distance rather than a torsion
#: because the leg starts from the value read off the structure at the press.
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

#: Water: a different molecule, for the marks that must stop being offered.
_WATER = """3
water
O            0.00000000000000        0.00000000000000        0.00000000000000
H            0.75800000000000        0.58600000000000        0.00000000000000
H           -0.75800000000000        0.58600000000000        0.00000000000000
"""


@pytest.fixture
def tab(tmp_path):
    """A whole Submit tab, which is where a scan's own path runs."""
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


def _an_editor(text=_ETHANE):
    """One structure editor over a coordinate box of its own."""
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


def _shown(widget):
    return (widget.layout.display or '') != 'none'


def _values(box):
    return [str(value) for _label, value in box.options]


def _method(part, prefix):
    for _label, value in part.submit_ff_dd.options:
        if str(value).lower().startswith(prefix):
            return value
    return None


def _said(refs):
    return ' '.join(refs['editor_state'].get('mol_status_lines') or ())


def _arm(refs, atoms, value, way='to', to=2.20, steps=8):
    refs['submit_pick_sync'].value = atoms
    refs['submit_internal_value'].value = value
    refs['submit_scan_way'].value = way
    if way == 'to':
        refs['submit_scan_to'].value = to
    refs['submit_scan_steps'].value = steps
    refs['submit_scan_btn'].click()


def _wait(refs, budget=300):
    began = time.time()
    while refs['editor_state'].get('scan_run') and time.time() - began < budget:
        time.sleep(0.02)
    time.sleep(0.2)


def _run_a_scan(refs, method='gfnff', how='hold', whole=True, way='to',
                to=2.20, steps=8, legs=(('1,2', 1.53),)):
    refs['submit_ff_dd'].value = method
    refs['submit_scan_how'].value = how
    refs['submit_scan_whole'].value = whole
    for atoms, value in legs:
        _arm(refs, atoms, value, way=way, to=to, steps=steps)
    refs['submit_scan_run_btn'].click()
    _wait(refs)


def _the_path_work_is_on_screen(refs):
    """Both boxes, the start on the pair, and the three ways from it."""
    # The pair is the only start there is now -- what is on screen belongs to
    # the press beside this one -- so the box holds one entry and follows the
    # rule every box here follows: a list of one is not a choice and is not
    # shown. What is on screen is the press, and the four ways from a pair.
    return (_values(refs['submit_saddle_from']) == ['scan']
            and refs['submit_saddle_from'].value == 'scan'
            and _shown(refs['submit_saddle_btn'])
            and _values(refs['submit_saddle_how']) == ['orca', 'neb',
                                                            'hand', 'walk']
            and _shown(refs['submit_saddle_how']))


@_needs_xtb
@pytest.mark.parametrize('shape,how', [
    ('the whole profile', dict(whole=True, way='to', to=2.20, steps=8)),
    ('walking to the next minimum', dict(whole=False, way='out', steps=8)),
    ('pushing with a force', dict(how='push', whole=False, way='in', steps=6)),
    ('two legs walked together',
     dict(whole=True, way='to', to=2.00, steps=5,
          legs=(('1,2', 1.53), ('2,3', 1.53)))),
])
def test_every_shape_of_scan_leaves_the_path_work_on_screen(tab, shape, how):
    """Not only the happy path: four real walks, four real toolbars.

    The one the existing tests drive is the whole profile of a distance, and
    it was the one that worked.  These are the others the user has in front of
    him -- the walk that stops at the next minimum, which is what the toolbar
    opens on; the push with a force, which is what the ``how`` box opens on;
    and a concerted step with two coordinates armed together.  Each ends with
    the pair in hand, the start standing on it and the three ways from a pair
    beside it, because what a walk left has nothing to do with which kind of
    walk it was.
    """
    refs = tab
    assert _values(refs['submit_saddle_from']) == ['here']
    _run_a_scan(refs, **how)

    assert refs['editor_state'].get('scan_ends'), shape
    assert _the_path_work_is_on_screen(refs), (
        shape, _values(refs['submit_saddle_from']),
        _values(refs['submit_saddle_how']))
    assert 'It left two ends' in _said(refs), shape


@_needs_xtb
def test_a_walk_cut_short_with_the_stop_button_leaves_its_two_ends(tab):
    """Stopping is a way of finishing, and it leaves what it walked.

    The Stop button breaks out of the walk, which is a different exit from
    the one a completed walk takes, and it is the one shape of interruption
    the user reaches by hand.  What he has afterwards is a start and the
    furthest point the walk got to -- two structures, and the path work is
    about them.
    """
    refs = tab
    refs['submit_ff_dd'].value = 'gfn2'
    refs['submit_scan_how'].value = 'hold'
    refs['submit_scan_whole'].value = True
    _arm(refs, '1,2', 1.53, way='to', to=8.0, steps=200)
    refs['submit_scan_run_btn'].click()
    # Pressed while it is walking, the way the user presses it: as soon as the
    # row says a second point has been reached.
    began = time.time()
    while time.time() - began < 120:
        if 'step 2 of' in _said(refs):
            break
        time.sleep(0.02)
    refs['submit_scan_run_btn'].click()
    _wait(refs)

    assert refs['editor_state'].get('scan_stop') is True
    assert refs['editor_state'].get('scan_ends')
    assert _the_path_work_is_on_screen(refs)


def _an_xtb_that_gives_up_after(monkeypatch, many):
    """Real xtb for *many* optimisations, and then the answer a real one gives.

    Measured first with a real binary: butane's central C-C walked out to 2.20
    A under GFN2, with an xtb driven into an SCF it cannot solve from the
    fourth optimisation on.  It answers ``ok`` false, ``energy`` None and a
    status naming the abort, and writes no relaxed geometry -- which is what
    the editor sees at a hard point of any real scan against its own 60-cycle
    limit.  Pinned here rather than driven, so the walk before the failure is
    still a real walk and the failure is still the shape a real one has.
    """
    from delfin.dashboard import gfn_optimize

    real = gfn_optimize.optimize_with_gfn
    counted = {'n': 0}

    def _answer(xyz_text, *args, **kwargs):
        if not kwargs.get('optimise', True):
            return real(xyz_text, *args, **kwargs)
        counted['n'] += 1
        if counted['n'] <= many:
            return real(xyz_text, *args, **kwargs)
        return {'ok': False, 'xyz': xyz_text, 'energy': None,
                'method': 'gfn2', 'seconds': 0.1, 'frames': [],
                'status': ('GFN2-xTB stopped with an error: optimizer_relax: '
                           'SCF not converged, aborting... The structure was '
                           'left as it was.')}

    monkeypatch.setattr(gfn_optimize, 'optimize_with_gfn', _answer)


@_needs_xtb
def test_a_walk_whose_xtb_stops_answering_still_leaves_its_two_ends(
        tab, monkeypatch):
    """The defect, reproduced and then fixed, in the terms he reported it.

    Three points walked, their two ends in hand, and a toolbar that showed
    ``To the saddle`` and the Climb switch and nothing else -- the editor
    saying there was nothing further to do with a structure it had just spent
    minutes making something of.  The pair is now written down in the block
    every way out of the walk goes through, so an xtb that gives up costs the
    path work nothing.
    """
    refs = tab
    _an_xtb_that_gives_up_after(monkeypatch, 3)
    _run_a_scan(refs, method='gfn2', steps=8)

    assert refs['editor_state'].get('scan_ends'), 'the walk lost its two ends'
    assert _the_path_work_is_on_screen(refs)
    assert 'It left two ends' in _said(refs)


@_needs_xtb
def test_a_walk_that_gave_up_says_so_where_the_walk_is_reported(tab,
                                                                monkeypatch):
    """The reason was written and then replaced by the verdict.

    "The scan stopped at step 4: ..." went into the row and the finishing
    block wrote the verdict over it a moment later, so what the user was left
    reading was "The scan walked 3 of 8 points ... the whole path is open" --
    a description of a path that had never been walked to the end.  The reason
    now goes into the verdict itself, which is the sentence that survives.
    """
    refs = tab
    _an_xtb_that_gives_up_after(monkeypatch, 3)
    _run_a_scan(refs, method='gfn2', steps=8)

    said = _said(refs)
    assert 'The scan walked 2 of 8 points' in said, said
    assert 'It stopped at step 3 of 8' in said, said
    assert 'SCF not converged' in said, said


@_needs_xtb
def test_a_scan_that_left_nothing_to_walk_between_says_that_instead(
        tab, monkeypatch):
    """One point priced is not two ends, and the editor says which it has.

    A walk that gets one point has its start and its end at the same geometry;
    a path finder given the same structure twice has nothing to walk.  So
    there is correctly nothing to offer -- and the verdict says so, because an
    empty toolbar is not an account of that, it is the same silence the user
    reported.
    """
    refs = tab
    _an_xtb_that_gives_up_after(monkeypatch, 2)
    _run_a_scan(refs, method='gfn2', steps=8)

    assert not refs['editor_state'].get('scan_ends')
    assert not _shown(refs['submit_saddle_btn'])
    said = _said(refs)
    assert 'It left no two ends' in said, said
    # And what to do about it, in the same clause rather than in a
    # paragraph of its own: the user is reading this after every scan.
    assert 'mark two structures by hand' in said, said


@_needs_xtb
def test_undo_back_past_a_scan_takes_the_path_work_with_it(tab):
    """A scan makes the pair, so walking back past the scan unmakes it.

    "Undo muss dann natuerlich in Teilen auch wieder Funktionen nehmen."  A
    scan is three entries -- where it came to, the highest point it crossed,
    where it started -- so the toolbar is checked at every station and not
    only at the two ends of the chain.  Inside the walk the pair is still
    about this molecule and stays; past the scan it goes, and the box's own
    value goes with it, or the press would go on meaning the walk while
    standing over a structure that never had one.
    """
    refs = tab
    _run_a_scan(refs, steps=6)
    assert _the_path_work_is_on_screen(refs)

    seen = []
    for _press in range(4):
        refs['submit_manip_undo_btn'].click()
        seen.append((bool(refs['editor_state'].get('scan_ends')),
                     _values(refs['submit_saddle_from']),
                     str(refs['submit_saddle_from'].value)))
    # Two stations inside the scan, then out the far side of it.
    assert seen[0] == (True, ['scan'], 'scan'), seen
    assert seen[1] == (True, ['scan'], 'scan'), seen
    assert seen[2] == (False, ['here'], 'here'), seen
    assert seen[3] == (False, ['here'], 'here'), seen
    assert not _shown(refs['submit_saddle_from'])
    assert not _shown(refs['submit_saddle_how'])
    # And nothing is left waiting: a wish for the pair that outlived the undo
    # would re-select itself the moment anything made a pair again.
    assert refs['editor_state'].get('saddle_start_wish') is None

    for _press in range(4):
        refs['submit_manip_redo_btn'].click()
    assert refs['editor_state'].get('scan_ends')
    assert _the_path_work_is_on_screen(refs)


def test_a_mark_from_another_molecule_is_not_offered():
    """The other start, checked the way the scan's ends are.

    A mark outlives the structure it was made on -- that is what it is for,
    the other end is loaded afterwards -- and it went on being offered
    whatever was loaded next.  An ethane marked with a water on screen put
    "the end you marked" in the box over a press that walks eight atoms into
    three: a control that is on screen and would be refused is the same
    defect as one that is missing when it would work.
    """
    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')
    part.on_submit_path_from()

    # The other end, same molecule: the start appears, which is the route the
    # user takes when he has not run a scan.
    part.coords_widget.value = _STRETCHED
    part._refresh_saddle_controls()
    assert _values(part.submit_saddle_from) == ['marked']
    assert _shown(part.submit_saddle_btn)
    part.submit_saddle_from.value = 'marked'
    assert _values(part.submit_saddle_how) == ['orca', 'neb', 'hand',
                                               'walk']
    assert _shown(part.submit_saddle_how)

    # Another molecule, the mark still held: it is not a pair for this one.
    part.coords_widget.value = _WATER
    part._refresh_saddle_controls()
    assert _values(part.submit_saddle_from) == ['here']
    assert not _shown(part.submit_saddle_from)
    assert part._marked_pair() is None
    assert part._path_ends('marked') is None
    assert 'different molecule' in ' '.join(state.get('mol_status_lines') or ())


def test_a_mark_is_not_an_undo_step_and_is_not_taken_back_by_one():
    """Marking changes no coordinate, so it is not a step to walk back.

    Everything Undo restores is a state some action left behind, and marking
    an end leaves none: it names a structure the user already has.  It is
    therefore deliberately absent from what a history entry carries -- were it
    in there, undoing a drag made *after* the mark would silently drop the
    mark as well, which is the same lie in the other direction.  The pair a
    scan leaves is a different thing: the scan computed it, so the scan owns
    it and Undo takes it back.
    """
    part, state = _an_editor()
    part.submit_ff_dd.value = _method(part, 'gfn2')
    part.on_submit_path_from()
    part.coords_widget.value = _STRETCHED
    part._refresh_saddle_controls()
    assert _values(part.submit_saddle_from) == ['marked']

    part._remember('something else')
    part._undo_structure()
    assert state.get('path_from')
    assert 'marked' in _values(part.submit_saddle_from)


def test_the_two_ends_are_written_down_where_every_exit_goes_through():
    """The one line of control flow the defect was.

    The pair used to be assigned at the end of the walk, after the loop; the
    four early exits -- xtb not answering, a push that cannot be priced, one
    with no zero to measure from -- leave through the ``finally`` instead and
    so never reached it.  Read off the source because it is a statement about
    where the line stands, which is exactly what went wrong, and a test that
    only drove a failing walk would pass again the day someone moved it back.
    """
    from editor_source import EDITOR_SOURCE as source

    scan = source.split('def on_submit_scan_run')[1].split('def _said_modes')[0]
    finally_block = scan.split('def _work():')[1].split('\n            finally:')[1]
    assert "state['scan_ends'] = (" in finally_block
    # And before the budget clamps what goes into the box, which is a
    # different question from what the walk reached.
    assert (finally_block.index("state['scan_ends'] = (")
            < finally_block.index("state['scan_walled'] = costs[walked]"))
