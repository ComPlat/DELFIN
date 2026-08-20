"""A bug report from the viewer that can be played back into a viewer.

Four defects were found in the structure editor in one afternoon of using it:
the playback stuttered, an interactive climb arrived at the wrong structure, a
scan result jumped to another geometry the moment it was grabbed, and toolbar
buttons were drawn past the right edge of the screen.  All four were invisible
to the 10,483 tests that were green while they were happening.  Nothing about
them is a wrong number -- the trajectories were right and the geometries were
right -- and one of them is still open because nobody could reproduce the
gesture.

So the report this file is about is not a complaint form.  Its whole value is
that it is *replayable*: the editor speaks a small closed vocabulary between
the page and the kernel, the journal records that vocabulary with its timings,
and a maintainer plays it back into a real editor and watches the same thing
happen.  That is what the first test here measures, and it measures it the
only way that means anything -- by recording a real gesture, writing it to
disk, reading it back and driving a second real editor with it, then comparing
the two message sequences message for message.

The rest of the file is what makes that affordable and what makes it fit on a
disk: what recording costs per message, what the ring drops when a session
outlasts it, and that a truncated report still replays from the earliest
structure it kept.
"""

import json
import statistics
import time

import numpy as np
import pytest

from delfin.dashboard import climb as _c
from delfin.dashboard import editor_journal as ej

_needs_xtb = pytest.mark.skipif(
    not _c.have_fast_gradients() and _c._gfn.find_xtb() is None,
    reason='no xtb to take gradients from')


#: Butadiene and ethene 3.35 A apart, relaxed under GFN2.  The same sixteen
#: atoms the climb tests use, so the gesture recorded here is one that has
#: been driven through this editor before and is known to do real work: the
#: drag pulls the ethene in, the release minimises, and the scan walks one of
#: the forming bonds down to 1.75 A.
_COMPLEX = """16
butadiene and ethene, relaxed apart
C      -1.511770730818    -0.078246613170     0.455478938079
C      -0.727043823107     0.952677742712     0.160901830656
C       0.726619700240     0.953646170753     0.162987634387
C       1.511848933658    -0.076870712366     0.457589167929
H      -2.584826367129     0.006765122345     0.431252533856
H      -1.117792846415    -1.040460185517     0.733672017164
H      -1.183044319988     1.897378083656    -0.109168104333
H       1.182127851949     1.899269641514    -0.104743334819
H       2.585025855997     0.009342483188     0.435050712651
H       1.118707557773    -1.040095122972     0.731436720523
C      -0.658327621243    -0.317760213959     3.689418492198
C       0.658004039895    -0.316146847946     3.686413586994
H      -1.229591522504     0.583529699796     3.544725727827
H      -1.228837918327    -1.220460679314     3.836580781023
H       1.226373947517     0.586709172705     3.537742888343
H       1.231462082883    -1.217451713317     3.830574479147
"""


def _an_editor(room):
    """A real editor over the complex, with the Submit tab's host behaviour.

    The one respect that matters here is that a write to the coordinate box
    refreshes ``current_xyz_for_copy``: without it every optimiser starts from
    the structure the session opened on rather than from the one the hand just
    made, and a replay would then be compared against a molecule nobody was
    looking at.
    """
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room.mkdir(parents=True, exist_ok=True)
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda _script: None
    state = {'editor_host': 'Submit'}
    box = widgets.Textarea(value=_COMPLEX)

    def update_view(*_a, **_k):
        raw = (box.value or '').strip()
        if not raw:
            return
        rows = [one for one in raw.split('\n') if one.strip()]
        body = rows[2:] if rows and rows[0].strip().isdigit() else rows
        state['current_xyz_for_copy'] = {
            'content': f'{len(body)}\nEdited in DELFIN viewer\n'
                       + '\n'.join(body)}
        state['manip_inflight'] = False

    part = structure_editor.build(
        ctx, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda call, *a, **k: call(*a, **k),
        update_view=update_view,
        get_smiles_charge=lambda *a, **k: None)
    box.observe(lambda _change: update_view(), names='value')
    update_view()
    return part, state, box


_SERIAL = [0]


def _command(part, verb, payload=''):
    """What the page sends when it cannot finish a gesture on its own."""
    _SERIAL[0] += 1
    part.submit_cmd_sync.value = f'{verb}:{_SERIAL[0]}:{payload}'


def _where(text):
    rows = [one for one in (text or '').split('\n') if one.strip()]
    body = rows[2:] if rows and rows[0].strip().isdigit() else rows
    return np.array([[float(x) for x in row.split()[1:4]] for row in body])


def _quiet(part, state, seconds=300.0):
    """Wait until nothing is armed and nothing is walking."""
    began = time.time()
    time.sleep(2.0)
    while time.time() - began < seconds:
        busy = (state.get('climb_run') is not None
                or state.get('optimize_run') is not None
                or state.get('scan_run')
                or state.get('gfn_settle_busy')
                or state.get('gfn_restart_armed')
                or state.get('gfn_minimise_armed')
                or state.get('gfn_settle_armed')
                or part.submit_optimize_btn.value)
        if not busy:
            break
        time.sleep(0.15)
    time.sleep(1.0)


def _a_real_gesture(part, state, box, *, drags=6, depth=0.95, pause=0.2):
    """Grab, drag, let go, scan, stop -- in the order the page sends them.

    ``ej.press`` rather than a plain assignment for every control, because
    that is the difference the journal records: a value set from Python is
    the editor moving its own switch, and only a change carrying ipywidgets'
    property lock is a hand.  A test that assigned directly would record a
    session in which the user pressed nothing.
    """
    ej.press(part.submit_ff_dd, 'gfn2')
    ej.press(part.submit_relax_btn, True)
    time.sleep(0.5)

    names = _c._elements(_COMPLEX)['symbols']
    start = _where(box.value)
    _command(part, 'gfngrab', '0,0')
    for turn in range(1, drags + 1):
        moved = start.copy()
        moved[10:16] += np.array([0.0, 0.0, -depth * turn / drags])
        part.submit_manip_sync.value = _c.xyz_document(
            names, moved, f'DELFIN drag-follow held=10,11 n={turn}')
        time.sleep(pause)
    part.submit_manip_sync.value = _c.xyz_document(
        names, _where(box.value), 'DELFIN drag-end')
    _command(part, 'gfnfree', '')
    _quiet(part, state)

    part.submit_pick_sync.value = '0,10'
    ej.press(part.submit_scan_how, 'hold')
    ej.press(part.submit_scan_way, 'to')
    ej.press(part.submit_scan_to, 1.75)
    ej.press(part.submit_scan_steps, 8)
    ej.press(part.submit_scan_whole, True)
    part.submit_scan_btn.click()
    part.submit_scan_run_btn.click()
    began = time.time()
    while state.get('scan_run') and time.time() - began < 240:
        time.sleep(0.05)
    _quiet(part, state)
    _command(part, 'gfnplay',
             f'stopped at frame 3 of run {state.get("gfn_run", 0)}')
    time.sleep(0.5)
    ej.press(part.submit_relax_btn, False)
    _quiet(part, state)


@_needs_xtb
def test_a_recorded_session_replays_as_the_same_sequence(tmp_path,
                                                         monkeypatch):
    """The whole feature, and the only test of it that proves anything.

    A real editor is driven through a real gesture -- grab, six drag-follows
    at 200 ms, drag-end, release, a minimisation, an eight-step GFN2 scan down
    to 1.75 A, and a stop at frame 3 -- the bug button is pressed, and the
    report that lands on disk is read back and played into a *second* real
    editor.  What is compared is the page-to-kernel half of both journals,
    message for message, string for string.

    Measured on this machine: the session ran 14 s and recorded 81 events, of
    which 21 are messages the page sent; the replay wrote 21 and the second
    editor's own journal recorded 21, and the two sequences are equal.  The
    structure the two sessions ended on agrees to 0.0000 A RMSD, which is the
    stronger statement: the replay did not merely repeat the messages, it
    arrived at the same molecule.

    It did not always.  The first version of this compared as equal while the
    replay never ran the scan at all -- Run scan is a *button*, a button has
    no value, and nothing about pressing one reaches a value observer.  The
    message sequences matched because a press is not a message, and the end
    geometries differed by 0.72 A.  That is why button presses are recorded as
    their own kind of event, and why this test still asserts on more than the
    sequence.

    What it asserts on is not the end geometry, and the reason is worth
    writing down because the obvious assertion is wrong.  This gesture crosses
    a barrier: the sixth step of the scan is the cycloaddition, +4.9 to about
    -21 kcal/mol in one step, and past it the surface has more than one
    product geometry within reach.  Comparing the two end structures fails
    about one run in four, always with the same value -- 0.814 A -- because it
    is a switch between two definite answers rather than drift.

    It is not the replay that is doing that.  Measured: the same scan, from a
    structure written to a file and read back byte for byte, run five times in
    five separate processes with no replay anywhere in it, walked the same
    first seven steps every time and ended at -54.8 kcal/mol four times and
    -54.9 once.  So the end geometry is not a property of the replay and this
    test may not assert it.

    What it asserts instead is everything that *is* reproducible, and all of
    it is: the run claims, which say the scan and the drag really ran; the
    whole structure history up to the barrier, which is identical to 0.000 A
    -- the six follow answers at 3.201, 3.049, 2.898, 2.748, 2.598, 2.450 A
    and the minimisation at 3.304; and the coordinate the scan was told to
    walk, which is 1.75 A on both sides by construction and 3.30 if the scan
    never ran.
    """
    pytest.importorskip('ipywidgets')
    monkeypatch.setenv('DELFIN_VIEWER_BUG_ARCHIVE', str(tmp_path / 'bugs'))

    part, state, box = _an_editor(tmp_path / 'live')
    _a_real_gesture(part, state, box)
    ej.press(part.submit_bug_btn, True)
    ej.press(part.submit_bug_note, 'the scan result jumps when I grab it')
    part.submit_bug_send.click()

    reports = sorted((tmp_path / 'bugs').iterdir())
    assert len(reports) == 1, reports
    meta, timeline = ej.load_report(reports[0])
    assert meta['description'] == 'the scan result jumps when I grab it'
    assert meta['tab'] == 'Submit'

    recorded = ej.page_messages(timeline)
    assert any(kind == 'cmd' and 'gfngrab' in payload
               for kind, _name, payload in recorded), recorded
    assert any(kind == 'manip' and 'drag-follow' in payload
               for kind, _name, payload in recorded), recorded
    assert any(kind == 'manip' and 'drag-end' in payload
               for kind, _name, payload in recorded), recorded
    assert any(kind == 'click' and name == 'submit_scan_run_btn'
               for kind, name, _payload in recorded), recorded

    # Which run numbers were claimed, and by what. Recorded rather than read
    # back off the frames, because a run that takes a number and then draws
    # nothing is the shape two of this week's defects had, and off the frames
    # alone that is indistinguishable from a number nobody ever took.
    claims = [(e['v'], e['by']) for e in timeline if e['k'] == 'run']
    assert [by for _run, by in claims].count('follow') == 1, claims
    assert 'scan' in [by for _run, by in claims], claims
    assert [run for run, _by in claims] == sorted(
        run for run, _by in claims), claims

    again, again_state, again_box = _an_editor(tmp_path / 'replay')
    wrote = ej.replay(timeline, again, pace=1.0, max_gap=2.0)
    _quiet(again, again_state)
    second = again_state['editor_journal'].timeline()
    replayed = ej.page_messages(second)

    assert wrote == recorded, (len(wrote), len(recorded))
    assert replayed == recorded, (len(replayed), len(recorded))

    # The scan really ran in the replay, and so did the drag. This is what the
    # end geometry used to stand in for, said directly and without depending
    # on any chemistry: a run number is claimed when a run begins, so a Run
    # scan that was not replayed leaves no run claimed by the scan. That is
    # exactly the case that once passed while the replay skipped the scan.
    again_claims = [(e['v'], e['by']) for e in second if e['k'] == 'run']
    assert [by for _run, by in again_claims].count('follow') == 1, again_claims
    assert 'scan' in [by for _run, by in again_claims], again_claims
    assert [by for _run, by in again_claims] == [
        by for _run, by in claims], (claims, again_claims)

    # And it was the same molecule all the way to the barrier: every structure
    # written to the box, in order, with the comment saying what produced it.
    told = ej.answers(timeline)
    barrier = next((i for i, one in enumerate(told) if one[1] == 'Scanned'),
                   len(told))
    parted = ej.first_difference(timeline, second)
    assert parted is None or parted['at'] >= barrier, parted

    # The coordinate the scan was told to walk is where it was told to walk
    # to, on both sides. Whole-molecule agreement is not asserted past this
    # point and must not be -- see the docstring.
    def _walked(text):
        rows = _where(text)
        return float(np.linalg.norm(rows[0] - rows[10]))

    assert abs(_walked(box.value) - 1.75) < 0.05, _walked(box.value)
    assert abs(_walked(again_box.value) - _walked(box.value)) < 0.05, (
        _walked(box.value), _walked(again_box.value))


def test_where_two_runs_of_the_same_session_stopped_agreeing():
    """A replay is judged on the messages and read on the answers.

    "Did the replay work" has no yes-or-no answer past a certain point in a
    session, and the first time that was asked it took three throwaway scripts
    to find out where the two runs parted.  So the reading is a function now,
    and what it reports is the sentence a maintainer wants -- "they were the
    same molecule until the sixth step of the scan" -- rather than a number.

    Two things are compared and nothing else: every structure written to the
    coordinate box with the comment saying what produced it, and every run
    number claimed with the name of what claimed it.  The status lines are
    left out although the journal keeps them, because they carry the timings
    of the run that wrote them -- "holding 2 atoms, 164 ms each" -- so two
    runs of one session never agree on them and a comparison including them
    would report a difference every time and mean nothing.
    """
    def one(comment, shift=0.0):
        return {'k': 'box', 't': 0.0,
                'v': f'2\n{comment}\nC 0.0 0.0 0.0\nH {1.0 + shift} 0.0 0.0\n'}

    walk = [one('butadiene and ethene, relaxed apart'),
            {'k': 'run', 't': 0.1, 'v': 1, 'by': 'follow'},
            {'k': 'status', 't': 0.2, 'lines': ['holding 2 atoms, 164 ms each']},
            one('Optimised in DELFIN viewer'),
            {'k': 'run', 't': 0.3, 'v': 2, 'by': 'scan'},
            one('Scanned')]

    # The same walk twice agrees, and it agrees despite the status lines
    # differing, which is the whole reason they are not compared.
    other = [dict(e) for e in walk]
    other[2] = {'k': 'status', 't': 0.2,
                'lines': ['holding 2 atoms, 402 ms each']}
    assert ej.first_difference(walk, other) is None

    # A structure that came out somewhere else is named, with how far.
    moved = [dict(e) for e in walk]
    moved[3] = one('Optimised in DELFIN viewer', shift=0.6)
    parted = ej.first_difference(walk, moved)
    assert parted['kind'] == 'box'
    assert parted['recorded'] == 'Optimised in DELFIN viewer'
    assert 0.4 < parted['apart'] < 0.5, parted
    # box, run, box, run, box -- the status line is not an answer, so
    # the second structure is the third thing compared.
    assert parted['at'] == 2, parted

    # A run nobody claimed -- which is what a button that was not replayed
    # looks like. Dropping it slides everything after it up, so what the
    # reader is told is that a scan was expected where a structure arrived.
    lost = [e for e in walk if not (e['k'] == 'run' and e['by'] == 'scan')]
    parted = ej.first_difference(walk, lost)
    assert parted['recorded'] == '2 by scan', parted
    assert parted['replayed'] == 'Scanned', parted

    # And when a side simply runs out -- the shape the real defect had, since
    # the scan was the last thing in that session -- it is named as an
    # absence rather than as a difference between two things.
    parted = ej.first_difference(walk, walk[:3])
    assert parted['kind'] == 'missing', parted
    assert parted['at'] == 2, parted
    assert 'nothing here' in parted['said'], parted

    # And a structure written under a different name is a difference even
    # when the coordinates agree: the comment says what produced it, and
    # "Scanned" arriving where "Optimised" was expected is the report.
    renamed = [dict(e) for e in walk]
    renamed[5] = one('Optimised in DELFIN viewer')
    parted = ej.first_difference(walk, renamed)
    assert parted['recorded'] == 'Scanned', parted
    assert parted['apart'] is None, parted


def test_what_a_journal_keeps_of_a_gesture(tmp_path):
    """A report has to answer what happened, not only what was pressed.

    Half of what was reported over those two days was about how a structure
    *moved* -- a stutter, a methyl appearing to act up, a jump to another
    geometry -- and none of that is visible in a list of button presses.  So
    the journal keeps the geometry as it changed hands: the structure it
    opened on, every write to the coordinate box with the comment saying what
    produced it, every window the frame channel sent, and every run number
    claimed with the name of what claimed it.

    Driven here without a server method, so what is measured is the kernel's
    own bookkeeping rather than xtb's arithmetic: the browser force field is
    chosen, an atom is picked up, six follow messages arrive and the hand
    lets go.
    """
    pytest.importorskip('ipywidgets')

    part, state, box = _an_editor(tmp_path / 'gesture')
    journal = state['editor_journal']
    names = _c._elements(_COMPLEX)['symbols']
    start = _where(box.value)

    ej.press(part.submit_relax_btn, True)
    _command(part, 'gfngrab', '4,2')
    for turn in range(1, 7):
        moved = start.copy()
        moved[10:16] += np.array([0.0, 0.0, -0.1 * turn])
        part.submit_manip_sync.value = _c.xyz_document(
            names, moved, f'DELFIN drag-follow held=10,11 n={turn}')
    part.submit_manip_sync.value = _c.xyz_document(
        names, start, 'DELFIN drag-end')
    _command(part, 'gfnfree', '')

    kinds = {}
    for event in journal.timeline():
        kinds[event['k']] = kinds.get(event['k'], 0) + 1

    # The opening is pinned outside the ring and carries the molecule.
    assert kinds.get('open') == 1
    assert journal.opened['xyz'].splitlines()[0].strip() == '16'
    assert journal.opened['widgets']['submit_ff_dd'] == 'uff'

    # Every message the page sent, verbatim, comment line included.
    follows = [e for e in journal.events
               if e['k'] == 'manip' and 'drag-follow' in e['v']]
    assert len(follows) == 6, kinds
    assert 'held=10,11' in follows[0]['v']
    assert any(e['k'] == 'manip' and 'drag-end' in e['v']
               for e in journal.events)
    assert [e['v'] for e in journal.events if e['k'] == 'cmd'] == [
        f'gfngrab:{_SERIAL[0] - 1}:4,2', f'gfnfree:{_SERIAL[0]}:']

    # The geometry as it changed hands: a structure per write to the box.
    boxes = [e for e in journal.events if e['k'] == 'box']
    assert len(boxes) >= 6, len(boxes)
    assert all(len(e['v'].splitlines()) >= 18 for e in boxes)

    # Nothing on the server ran, so nothing claimed a run number: with the
    # browser's own force field the page does the following. Which is worth
    # asserting rather than skipping -- a journal that invented a run here
    # would be inventing them everywhere. The claims are checked where they
    # really happen, in the replay test.
    assert [e for e in journal.events if e['k'] == 'run'] == []

    # And who moved each switch. Relax was pressed; the editor turned
    # Manipulate on underneath it, and those are not the same event.
    pressed = [e['n'] for e in journal.events
               if e['k'] == 'w' and e['by'] == 'page']
    by_editor = [e['n'] for e in journal.events
                 if e['k'] == 'w' and e['by'] == 'kernel']
    assert pressed == ['submit_relax_btn'], pressed
    assert 'submit_manip_btn' in by_editor, by_editor


def test_recording_costs_a_part_in_thirty_thousand_of_a_follow_step(tmp_path):
    """It has to be free enough to leave on, so this measures what it is not.

    Two paths carry traffic fast enough for the question to be real: the
    follow answers several times a second while a hand is down, and the frame
    channel carries the trajectory the player draws at up to 12,000 frames a
    second.  Both were measured with the ring already at its bound, so what is
    timed is the steady state a long session runs in and not an empty deque.

    Measured on this machine, Python 3.11:

        record() for a 1,013-byte drag-follow   2.14 us on, 0.20 us off
        record() for a 16-kB, 40-frame write    2.11 us on, 0.20 us off
        one follow message through the editor    132.5 us on, 127.6 us off
        one frame-channel write through it        23.7 us on,  20.5 us off

    So recording costs about 5 us on a follow message and 3 us on a frame
    write.  The editor's own measurement of what a GFN2 follow step costs on
    these sixteen atoms, taken in the same run, is 172 ms -- so the journal is
    about one part in thirty-five thousand of it.  On the frame channel the
    3 us is spread over the 40 frames the write carries, which is 0.08 us a
    frame: at 12,000 frames a second that is a tenth of a percent of one core.

    The assertions are deliberately loose -- twenty times the measured cost --
    because a shared machine's timings are not reproducible and a test that
    pinned them would fail for reasons that have nothing to do with this code.
    What they catch is the change that makes recording cost a millisecond.
    """
    pytest.importorskip('ipywidgets')

    names = _c._elements(_COMPLEX)['symbols']
    start = _where(_COMPLEX)
    follow = _c.xyz_document(names, start, 'DELFIN drag-follow held=10,11')
    frame = json.dumps({'run': 3, 'from': 0, 'pace': 83, 'follow': 1,
                        'frames': [[round(float(v), 4)
                                    for v in start.ravel()]] * 40})

    def per_call(journal, kind, payload, rounds=20_000):
        best = []
        for _ in range(3):
            began = time.perf_counter()
            for _ in range(rounds):
                journal.record(kind, v=payload)
            best.append((time.perf_counter() - began) / rounds)
        return min(best)

    # A small ring, so every call evicts: the steady state, not the best case.
    full = ej.Journal(max_events=2000, max_bytes=2 * 1024 * 1024)
    for _ in range(2000):
        full.record('manip', v=follow)
    off = ej.Journal(on=False)

    on_follow = per_call(full, 'manip', follow)
    off_follow = per_call(off, 'manip', follow)
    on_frame = per_call(full, 'frame', frame)
    assert full.dropped > 1000, full.dropped
    assert on_follow < 40e-6, on_follow
    assert on_frame < 40e-6, on_frame
    assert off_follow < 4e-6, off_follow

    # And through the editor's own path, where the observer chain is real.
    part, state, box = _an_editor(tmp_path / 'cost')
    journal = state['editor_journal']

    def through_the_editor(recording, rounds=200):
        journal.on = recording
        laps = []
        for turn in range(rounds):
            moved = start.copy()
            moved[10] += np.array([0.0, 0.0, -0.001 * (turn % 50)])
            text = _c.xyz_document(names, moved,
                                   f'DELFIN drag-follow held=10,11 n={turn}')
            began = time.perf_counter()
            part.submit_manip_sync.value = text
            laps.append(time.perf_counter() - began)
        return statistics.median(laps)

    through_the_editor(True, 50)                     # warm the path
    with_it = through_the_editor(True)
    without = through_the_editor(False)
    journal.on = True
    assert with_it - without < 1e-3, (with_it, without)


def test_a_session_longer_than_the_buffer_drops_its_beginning(tmp_path,
                                                              monkeypatch):
    """A long session must cost a bounded amount of memory and of disk.

    The ring is bounded twice, by a count of events and by a total of their
    bytes, and the oldest goes first.  Oldest-first is the right way round
    here and not merely the easy one: the button is pressed within seconds of
    seeing something go wrong, so the end of the buffer is the evidence and
    the beginning is a session the user has already forgotten.

    Measured for a five-minute session of continuous work at the default
    bound: 1,599 events, 1.42 MB of journal text, 0.288 MB written as
    ``journal.jsonl.gz`` and 0.324 MB for the whole report directory -- which
    is why the default bound is where it is, about two hours of that work.

    What this pins is the behaviour at the bound rather than those numbers:
    that both bounds bite, that the opening snapshot survives whatever else
    is dropped, that the report says in words how much went, and that a
    truncated journal still replays -- from the earliest structure it kept,
    because that is the one thing standing in for the opening it lost.
    """
    pytest.importorskip('ipywidgets')
    monkeypatch.setenv('DELFIN_VIEWER_BUG_ARCHIVE', str(tmp_path / 'bugs'))

    journal = ej.Journal(max_events=50, max_bytes=10 * 1024 * 1024)
    journal.opening({'submit_ff_dd': 'gfn2'}, _COMPLEX)
    for turn in range(400):
        journal.record('box', v=f'16\nstep {turn}\n' + 'C 0 0 0\n' * 16)
    told = journal.summary()
    assert told['events'] == 50
    assert told['dropped'] == 350
    assert journal.opened is not None, 'the opening was not pinned'
    assert journal.opened['xyz'] == _COMPLEX

    # The byte bound bites on its own, with the count nowhere near.
    tight = ej.Journal(max_events=1_000_000, max_bytes=20_000)
    for turn in range(400):
        tight.record('frame', v='x' * 2000)
    assert tight.bytes <= 20_000 + 2100, tight.bytes
    assert tight.dropped > 300, tight.dropped

    where = ej.write_report(journal, description='it went on a long time',
                            widgets={})
    said = (where / 'report.md').read_text(encoding='utf-8')
    assert 'longer than the buffer' in said
    assert '350 events' in said
    meta, timeline = ej.load_report(where)
    assert meta['summary']['dropped'] == 350

    # A truncated journal still knows where to start: the earliest structure
    # that survived stands in for the opening the ring dropped.
    lost = [event for event in timeline if event['k'] != 'open']
    assert ej.opening_structure(lost).splitlines()[1] == 'step 350'


def test_a_report_reads_as_a_report_and_lands_where_it_says(tmp_path,
                                                            monkeypatch):
    """The user's own sentence first, then what happened, then the sequence.

    A report nobody opens is worth as much as one nobody can replay, so the
    order is the order a person reads in: the sentence that says what to look
    for, the setup it happened under, the counts, the molecule, and only then
    the raw timeline.  The machine-readable half is a sibling file, not a
    different format for the same thing.

    Where it lands is checked here too, in the order the module resolves it:
    the environment override first, the settings key next, and the per-user
    fallback last.  It resolves without creating anything, so an interface can
    show a path on a machine where nobody has ever filed a report.
    """
    pytest.importorskip('ipywidgets')

    monkeypatch.delenv('DELFIN_VIEWER_BUG_ARCHIVE', raising=False)
    assert ej.resolve_archive_dir({}).name == 'viewer_bugs'
    assert ej.resolve_archive_dir(
        {'viewer': {'bug_archive_dir': '/tmp/somewhere'}}).as_posix() \
        == '/tmp/somewhere'
    monkeypatch.setenv('DELFIN_VIEWER_BUG_ARCHIVE', str(tmp_path / 'bugs'))
    assert ej.resolve_archive_dir({'viewer': {'bug_archive_dir': '/tmp/x'}}) \
        == tmp_path / 'bugs'
    assert not (tmp_path / 'bugs').exists(), 'resolving must not create it'

    part, _state, _box = _an_editor(tmp_path / 'read')
    ej.press(part.submit_gfn_charge, -1)
    _command(part, 'gfngrab', '7,2')
    _command(part, 'gfnfree', '')

    ej.press(part.submit_bug_btn, True)
    assert 'stays on this machine' in part.submit_bug_where.value
    assert 'bugs' in part.submit_bug_where.value
    ej.press(part.submit_bug_note, 'es zappelt beim Abspielen')
    part.submit_bug_send.click()

    where = sorted((tmp_path / 'bugs').iterdir())[0]
    assert {p.name for p in where.iterdir()} == {
        'report.md', 'report.json', 'journal.jsonl.gz',
        'structure_at_start.xyz'}

    said = (where / 'report.md').read_text(encoding='utf-8')
    head, rest = said.split('## Context', 1)
    assert 'es zappelt beim Abspielen' in head, 'the sentence is not first'
    assert rest.index('## The sequence') > rest.index(
        '## What was changed while it ran')
    assert 'gfngrab:' in rest
    assert '| tab | Submit |' in rest

    # What moved is a table of its own, ahead of the sixty controls that did
    # not: the four lines a reader came for are not worth finding in the
    # middle of a table nobody reads to the end.
    changed, standing = rest.split('## Everything else, as it stood', 1)
    changed = changed.split('## What was changed while it ran', 1)[1]
    assert 'submit_gfn_charge' in changed, changed
    assert 'submit_temperature' not in changed, changed
    assert 'submit_temperature' in standing

    meta = json.loads((where / 'report.json').read_text(encoding='utf-8'))
    assert meta['schema'] == 'delfin-viewer-bug-report/1'
    assert meta['widgets_at_send']['submit_gfn_charge'] == -1
    assert meta['atoms'] == '16'
    assert (where / 'structure_at_start.xyz').read_text().splitlines()[0] \
        == '16'

    # The line is emptied and the control closes, so the next report is not
    # filed under the last one's sentence.
    assert part.submit_bug_note.value == ''
    assert part.submit_bug_btn.value is False
    assert 'Reported.' in part.mol_status.value

    listed = ej.list_reports(tmp_path / 'bugs')
    assert [row['description'] for row in listed] == [
        'es zappelt beim Abspielen']
    assert listed[0]['tab'] == 'Submit'


def test_the_button_is_on_the_toolbar_and_out_of_the_way(tmp_path):
    """One button on the crowded row, and the line only when it is asked for.

    The rule about that row is that less is more, and it is not an aesthetic
    one: nineteen controls and 1,900 px of content were being laid out in a
    620 px row, and ten of them were painted past the right edge where nothing
    could reach them.  So the reporting control is one 44 px button until it
    is pressed, and the line to write in and Send appear underneath it and go
    away again.

    It is appended to the row rather than written into the list that builds
    it, which is what keeps it in one place of its own while the contents of
    that row are rearranged around it.

    That it can wrap is not asserted here: a row inside a row that says
    nowrap takes its whole content past the edge of the screen, and
    ``test_the_toolbar_stays_on_the_screen`` already asks that of every
    nested row in the toolbar rather than of the newest one.  It asked it of
    this one and the answer was no, which is how the group came to carry
    ``flex_flow='row wrap'``.
    """
    pytest.importorskip('ipywidgets')

    part, _state, _box = _an_editor(tmp_path / 'toolbar')
    assert part.submit_bug_group in part.submit_manip_toolbar.children
    assert part.submit_bug_group.children[0] is part.submit_bug_btn

    assert part.submit_bug_note.layout.display == 'none'
    assert part.submit_bug_send.layout.display == 'none'
    ej.press(part.submit_bug_btn, True)
    assert part.submit_bug_note.layout.display == 'flex'
    assert part.submit_bug_send.layout.display == 'flex'
    ej.press(part.submit_bug_btn, False)
    assert part.submit_bug_note.layout.display == 'none'

    # And a tab that holds the editor is handed the control with everything
    # else, under the name it has here.
    assert 'submit_bug_btn' in part.exported
    assert 'submit_bug_send' in part.exported
