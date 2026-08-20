"""What the viewer was doing, kept so that a defect can be walked through again.

Four real defects turned up in one afternoon of using the structure editor: a
stutter in the trajectory playback, an interactive climb that arrived at the
wrong structure, a scan result that jumped to a different geometry the moment
it was grabbed, and toolbar buttons drawn past the right edge of the screen.
None of the four was visible to the 10,483 tests that were green at the time.
The trajectories were right, the geometries were right, and every assertion
anyone had written held.  Only the hand on the mouse found them, and one of
them is still open because nobody could reproduce the gesture that caused it.

That is what this module is for, and it is why a bug report about the editor
is worth so much more than a sentence.  The editor is driven by a small,
closed vocabulary between the page and the kernel -- a command
``verb:serial:payload`` on one channel, an XYZ whose comment line reads
``DELFIN drag-follow held=10,11`` or ``DELFIN drag-end`` on another, an atom
pick on a third, and a widget the user changed -- and that vocabulary is
exactly what every test of the editor already drives it through.  Record the
sequence with its timings and a maintainer can play it back into a real
editor and watch the defect happen, which is the difference between "es
zappelt" and a failing test.

What is kept
------------

Everything a replay needs, and the structure as it moved:

``open``    the widget snapshot the session started from -- method, charge,
            multiplicity, solvent and its model, temperature, hand mode, the
            sliders and every toggle -- and the structure that was on screen.
``cmd``     a page-to-kernel command, verbatim: ``gfngrab:41:12,4``,
            ``gfnfree:42:``, ``gfnplay:43:stopped at frame 12 of run 4``.
``manip``   a page-to-kernel structure, verbatim, comment line included: this
            is a drag-follow, a drag-end, an undo in the browser, or the
            running force field reporting where it has got to.
``pick``    an atom picked in the picture.
``click``   a button pressed.  A button has no value, so nothing about it
            reaches a value observer -- and Run scan, Set, Centre, Reset and
            the two Undos are all buttons.  A journal without these looked
            complete and replayed a session in which the scan never ran.
``w``       a widget's value changed, to what, and by whom.  Both halves are
            worth keeping and they are not the same thing: the user pressing
            Climb to TS and the climb switching its own toggle off in the
            middle of his drag were one of this week's defects, and a journal
            that could not tell them apart would show only that the toggle
            went down.  Which one it was is read off ipywidgets' own property
            lock, which is set exactly while a change from the browser is
            being applied -- see :func:`press`.
``box``     a structure written to the coordinate box.  Each already carries
            a comment saying what produced it -- "Optimised in DELFIN
            viewer", "Past the budget: back to the last structure that was
            inside it" -- so this alone is the history of the geometry.
``frame``   a kernel-to-page write on the frame channel, verbatim.  This *is*
            how the structure moved between two boxfuls, and it is the only
            record that can show a stutter or a jump: the player draws from
            it, so whatever the player did wrong is in here.
``wall``    where the thermal budget said the hand may no longer go, and
            with what reach.
``run``     a run number claimed, and by what.  Which frames belong to which
            walk was itself the cause of two defects fixed that week, so the
            claim is recorded rather than inferred.
``status``  what the status line said, in the plain lines it was given rather
            than the rendered HTML.
``note``    a sentence from the editor about something that has no channel --
            the send itself, a refusal.

Two things a reader will look for are in here without having a kind of their
own, and deliberately so.  The structure *at a grab* is the last ``box``
before the ``gfngrab``, because the box is what the editor itself reads when
it asks where the molecule is -- a second copy would be a second answer to a
question that has one.  The structure *at a release* is the ``drag-end``
message, verbatim, which is the geometry the page sent from its mouseup
handler.  The frame number the page reported at the grab, and the walk it
counts along, ride in the ``gfngrab`` payload where the page put them.

What is deliberately not kept: no screenshot, no rendered HTML, no
force-field parameter set, no xtb scratch directory, and no second copy of a
frame -- the frame channel's own payload is stored as it went out, rounded to
four decimals, which is what the wire format already is and a thousandth of a
bond length.

Size, and what goes first
-------------------------

Recording is into memory and never touches the disk until the user presses
send.  The buffer is a ring bounded two ways -- a count of events and a total
of their bytes -- and when either is reached the *oldest* event is dropped.

Oldest-first is the right way round for this feature and not merely the easy
one: the button is pressed within seconds of seeing something go wrong, so
the end of the buffer is the part that matters and the beginning is a session
the user has already forgotten.  The ``open`` snapshot is pinned outside the
ring so it can never be dropped, and a replay of a truncated journal starts
from the earliest ``box`` event that survived.  How many events and bytes
were dropped is written into the report, in both files, so a maintainer is
never left guessing whether a replay began at the beginning.

Measured on this machine, Python 3.11, GFN2 on sixteen atoms: five minutes of
continuous work -- twenty-three rounds of grab, six drag-follows, release,
minimise, an eight-step scan and a stop -- produced 1,599 events and 1.42 MB
of journal text, which gzips to 0.288 MB, and the whole report directory came
to 0.324 MB.  The frame channel is 68% of those bytes, the coordinate box 17%
and the drag messages 12%; everything else together is 3%.  At that rate the
event bound is reached after about two hours of unbroken work and the byte
bound after about three, and nothing is dropped before then.

The frames are stored as they went out on the wire, which is already the
cheap form: the channel windows its own writes, and the coordinates in them
are rounded to four decimals before they leave -- half the JSON, and a
thousandth of a bond length.  Deltas were considered and are not worth what
they cost: gzip finds the same redundancy anyway, at 4.9x on that session,
and a delta stream cannot be read at all without replaying it from a base
frame, which is exactly the property this format is chosen for not having.

Where it lands
--------------

The same shape as the agent's bug reports next door, because that is the
convention the user asked for -- one directory per report, named for the UTC
time it was written, holding a ``report.md`` for a person and a
``report.json`` for a machine.  The path is never hard-coded; it resolves as

1. ``$DELFIN_VIEWER_BUG_ARCHIVE``,
2. ``settings["viewer"]["bug_archive_dir"]``,
3. ``~/.delfin/viewer_bugs``.

It stays on this machine.  These are the user's own structures, and there is
nothing here worth shipping anywhere: no push, no network, and the interface
says the directory it wrote to so there is never a question where the report
went.
"""

from __future__ import annotations

import gzip
import json
import os
import socket
import sys
import time
from collections import deque
from datetime import datetime, timezone
from pathlib import Path

# The agent's bug reports are next door and this follows them rather than
# inventing a second convention: the directory name, the group-readable
# fix-up and the credential redactor are theirs, and a viewer report is the
# same kind of object written from a different room.
from delfin.agent.bug_report import (
    _current_user,
    _delfin_version,
    _guarded_write,
    _make_group_readable,
    _report_dirname,
)

_FALLBACK_DIR = Path.home() / '.delfin' / 'viewer_bugs'

#: The bound.  Both are hit by the frame channel long before anything else
#: reaches them, and where they are put comes from the five-minute session
#: measured in the module docstring: about two hours of unbroken work to the
#: event bound and three to the byte bound.
MAX_EVENTS = 40_000
MAX_BYTES = 48 * 1024 * 1024

#: A single event may not be larger than this.  One write on the frame
#: channel carries at most two reads' worth of frames, so in ordinary use
#: nothing comes close; the cap exists because a pathological run -- a
#: thousand-atom chain that walked for a minute before anybody read it --
#: would otherwise put a single multi-megabyte string in the ring and evict
#: the whole gesture around it.
MAX_EVENT_BYTES = 2 * 1024 * 1024

#: The kinds the page sends the kernel.  These, and only these, are what a
#: replay writes back: everything else in a journal is the kernel's own
#: answer, and a replay that wrote those too would be reading its own
#: homework rather than driving the editor.
FROM_PAGE = ('cmd', 'manip', 'pick', 'w', 'click')

#: Not part of the gesture.  Two kinds: the channels, which are recorded as
#: messages a few lines below and would otherwise appear a second time as a
#: widget nobody pressed; and the reporting control itself, whose whole
#: content is the sentence the user typed -- which is the report's subject,
#: not a step in the session it is about.
_CHANNEL_WIDGETS = frozenset({
    'submit_cmd_sync', 'submit_manip_sync', 'submit_pick_sync',
    'submit_gfn_frame', 'submit_gfn_wall', 'submit_draw_sync',
    'submit_draw_frame', 'submit_manip_status', 'submit_internal_label',
    'submit_ff_notes',
    'submit_bug_note', 'submit_bug_btn', 'submit_bug_send',
    'submit_bug_where',
})


def resolve_archive_dir(settings=None):
    """Where viewer reports are written (see the module docstring for order).

    Does not create the directory, so a caller can show the path in the
    interface without leaving an empty folder behind on a machine where
    nobody ever files a report.
    """
    env = (os.environ.get('DELFIN_VIEWER_BUG_ARCHIVE') or '').strip()
    if env:
        return Path(env).expanduser()
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    configured = str(
        ((settings or {}).get('viewer') or {}).get('bug_archive_dir') or ''
    ).strip()
    if configured:
        return Path(configured).expanduser()
    return _FALLBACK_DIR


def _weigh(fields):
    """Roughly how many bytes this event will take as a line of JSON.

    Roughly, and on purpose: the exact answer is ``len(json.dumps(...))``,
    which is a second serialisation of every payload on a path that runs
    several times a second while a hand is being followed.  Adding up the
    lengths of the strings is what the payloads actually are -- a frame
    write and an XYZ are one long string each -- and it is a length lookup
    rather than a walk.

    Measured against the file it is estimating, over 1,500 events of drag
    messages, frame writes and status lines: 0.5% under.  That is far inside
    what the bound is for.
    """
    total = 24
    for value in fields.values():
        if type(value) is str:
            total += len(value) + 8
        elif type(value) is dict:
            total += 40 * len(value)
        else:
            total += 16
    return total


class Journal:
    """Everything the editor did, in a ring that cannot outgrow its bound.

    One of these belongs to one editor.  It is switched on by default and
    costs nothing worth measuring to leave on -- see :meth:`record` -- and
    it is written to disk only when the user presses send.
    """

    def __init__(self, *, max_events=MAX_EVENTS, max_bytes=MAX_BYTES,
                 on=True):
        self.on = bool(on)
        self.max_events = int(max_events)
        self.max_bytes = int(max_bytes)
        self.events = deque()
        self.bytes = 0
        self.dropped = 0
        self.dropped_bytes = 0
        self.opened = None
        self.began = time.time()
        self._clock = time.monotonic
        self._t0 = self._clock()

    # -- recording ---------------------------------------------------------

    def record(self, kind, **fields):
        """Put one event in the ring.

        The whole path a message takes through this module, and it is one
        dict, one append and two additions.  The numbers this was tuned
        against are in ``tests/test_the_viewer_remembers_what_you_did.py``,
        which measures the same call with recording on and off.

        Strings are kept by reference and never copied, which is why a frame
        write costs the same as a drag-follow however many frames are in it:
        the payload already exists, having just been handed to the widget.

        Called from worker threads as well as from the message handlers -- a
        status line written from inside a running xtb round arrives here on
        that round's thread -- and deliberately without a lock.  The deque's
        own append and popleft are atomic, so the *order* is safe; only the
        running byte count can drift when two threads add to it at once, and
        the effect of that is a ring that trims a little early or a little
        late.  The count is an estimate to begin with, the bound is a bound
        and not a contract, and a lock on a path that runs several times a
        second to protect an approximation would be paying for nothing.
        """
        if not self.on:
            return
        size = _weigh(fields)
        if size > MAX_EVENT_BYTES:
            # Too big to be worth the whole gesture around it. Recorded as
            # its own absence, so the report says a message was here rather
            # than showing an unexplained gap in the timeline.
            fields = {'dropped_bytes': size}
            kind = kind + '-too-large'
            size = _weigh(fields)
        event = dict(fields)
        event['t'] = round(self._clock() - self._t0, 4)
        event['k'] = kind
        event['_n'] = size
        self.events.append(event)
        self.bytes += size
        while (self.bytes > self.max_bytes
               or len(self.events) > self.max_events):
            if len(self.events) <= 1:
                break
            gone = self.events.popleft()
            self.bytes -= gone.get('_n', 0)
            self.dropped += 1
            self.dropped_bytes += gone.get('_n', 0)

    def opening(self, widgets, xyz):
        """Pin the state the session started from, outside the ring.

        Outside it because a replay has to have somewhere to start: with the
        ring full and the oldest events gone, this is the only thing left
        that says which method was chosen and what the molecule was.
        """
        self.opened = {
            't': 0.0, 'k': 'open',
            'widgets': dict(widgets or {}),
            'xyz': str(xyz or ''),
        }

    # -- the channels ------------------------------------------------------

    def watch(self, widgets):
        """Record every later change of every widget handed in.

        Handed the editor's own locals, so a control added to the toolbar
        tomorrow is recorded without anybody remembering to add it here --
        which is the failure mode this replaces, and the one that would
        leave a report silently missing the very toggle that caused the bug.
        The protocol channels are skipped: they are recorded as messages a
        few lines below, and recording them twice would put a drag-follow in
        the journal as a widget nobody pressed.

        A button is watched too, and separately, because a button has no
        value: nothing about pressing Run scan, Set, Centre or Reset reaches
        a value observer at all.  A journal without them looked complete and
        was not -- replayed, its scan never ran, and the comparison of
        message sequences still passed because a press is not a message.

        Returns the widgets whose value is watched, which is also what a
        report snapshots at send time.  Buttons are not in it: there is
        nothing about a button to snapshot.
        """
        kept = {}
        for name, thing in sorted((widgets or {}).items()):
            if not name.startswith('submit_') or name in _CHANNEL_WIDGETS:
                continue
            if hasattr(thing, 'on_click'):
                thing.on_click(self._button_pressed(name))
                continue
            if not hasattr(thing, 'observe') or not hasattr(thing, 'value'):
                continue
            kept[name] = thing
            thing.observe(self._widget_changed(name), names='value')
        return kept

    def _button_pressed(self, name):
        def note(_button):
            self.record('click', n=name)
        return note

    def _widget_changed(self, name):
        def note(change):
            if change.get('name') != 'value':
                return
            # Who moved it. ipywidgets holds a property lock over exactly the
            # window in which a value arriving from the browser is applied --
            # it is how the widget knows not to echo the change back -- and
            # observers are called inside that window. So a lock that is set
            # is a hand on the control, and a lock that is empty is the
            # editor changing its own switch.
            owner = change.get('owner')
            self.record('w', n=name, v=_plain(change.get('new')),
                        by=('page' if getattr(owner, '_property_lock', None)
                            else 'kernel'))
        return note

    def on_frame(self, change):
        """A write on the frame channel, kept verbatim."""
        if change.get('name') == 'value' and change.get('new'):
            self.record('frame', v=change['new'])

    def on_wall(self, change):
        """A write on the thermal-wall channel, kept verbatim."""
        if change.get('name') == 'value' and change.get('new'):
            self.record('wall', v=change['new'])

    def on_box(self, change):
        """A structure written to the coordinate box, kept verbatim.

        Whoever wrote it.  The comment line says what produced it, so the
        sequence of these is the history of the geometry with its own
        provenance attached and nothing has to be reconstructed from the
        events around it.
        """
        if change.get('name') == 'value' and change.get('new'):
            self.record('box', v=change['new'])

    # -- reading it back ---------------------------------------------------

    def snapshot(self, widgets):
        """What every watched widget is set to right now."""
        out = {}
        for name, thing in sorted((widgets or {}).items()):
            try:
                out[name] = _plain(thing.value)
            except Exception:
                continue
        return out

    def timeline(self):
        """The pinned opening, then the ring, in the order it happened."""
        out = [self.opened] if self.opened else []
        out.extend(self.events)
        return out

    def summary(self):
        """Counts a person reads before deciding whether to open the stream."""
        counts = {}
        for event in self.events:
            counts[event['k']] = counts.get(event['k'], 0) + 1
        span = self.events[-1]['t'] if self.events else 0.0
        return {
            'events': len(self.events),
            'bytes': self.bytes,
            'dropped': self.dropped,
            'dropped_bytes': self.dropped_bytes,
            'seconds': round(span, 2),
            'by_kind': counts,
            'max_events': self.max_events,
            'max_bytes': self.max_bytes,
        }


def _plain(value):
    """A widget value as something JSON will hold and compare.

    Tuples come back from JSON as lists, so they are made lists here and the
    two sides of a replay comparison cannot differ over a container type
    nobody chose.
    """
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (tuple, list)):
        return [_plain(one) for one in value]
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    return str(value)


# ---------------------------------------------------------------------------
# Writing one out
# ---------------------------------------------------------------------------

def opening_structure(timeline):
    """The molecule a replay of this journal has to start from.

    The pinned opening when there is one, and the earliest structure that
    survived the ring when there is not.  There are two ways to have no
    opening: a long session whose beginning was dropped, and an editor built
    over an empty coordinate box -- the Submit tab does exactly that, and
    every structure it ever shows arrives afterwards.
    """
    opening = next((e for e in timeline if e.get('k') == 'open'), None)
    text = str((opening or {}).get('xyz') or '')
    if text.strip():
        return text
    return next((str(e.get('v') or '') for e in timeline
                 if e.get('k') == 'box'), '')


def _readable_event(event):
    """One line of the timeline, for a person reading the report."""
    kind = event.get('k', '')
    when = f'{event.get("t", 0.0):9.3f}s'
    if kind == 'open':
        rows = str(event.get('xyz') or '').splitlines()
        return (f'{when}  open    {rows[0].strip() if rows else "0"} atoms, '
                f'{len(event.get("widgets") or {})} controls')
    if kind == 'cmd':
        return f'{when}  cmd     {event.get("v", "")}'
    if kind == 'w':
        hand = 'you  ' if event.get('by') == 'page' else 'editor'
        return (f'{when}  {hand}  {event.get("n", "")} = '
                f'{event.get("v", "")!r}')
    if kind in ('manip', 'box'):
        rows = str(event.get('v', '')).splitlines()
        note = rows[1].strip() if len(rows) > 1 else ''
        count = rows[0].strip() if rows else '?'
        word = 'page' if kind == 'manip' else 'box '
        return f'{when}  {word}    {count} atoms -- {note}'
    if kind == 'frame':
        raw = str(event.get('v', ''))
        try:
            said = json.loads(raw)
            frames = said.get('frames') or []
            return (f'{when}  frames  run {said.get("run")}, '
                    f'{len(frames)} from {said.get("from")}'
                    + (', final' if said.get('final') else '')
                    + (', follow' if said.get('follow') else ''))
        except ValueError:
            return f'{when}  frames  {len(raw)} bytes'
    if kind == 'wall':
        return f'{when}  wall    {str(event.get("v", ""))[:120]}'
    if kind == 'status':
        return f'{when}  said    {" / ".join(event.get("lines") or [])}'
    if kind == 'run':
        return (f'{when}  run     {event.get("v")} claimed by '
                f'{event.get("by", "?")}')
    if kind == 'click':
        return f'{when}  press   {event.get("n", "")}'
    if kind == 'pick':
        return f'{when}  pick    {event.get("v", "")}'
    if kind == 'note':
        return f'{when}  note    {event.get("v", "")}'
    return f'{when}  {kind}  {json.dumps(event, ensure_ascii=False)[:160]}'


def _render_markdown(meta, timeline, summary, *, tail=400):
    """The report a person reads: the sentence, then what happened, then the
    sequence.

    The user's own words first and unedited, because they are the only part
    nobody else can reconstruct -- everything below them is a machine's
    account of the same minutes, and it is worth reading only once the
    sentence has said what to look for.
    """
    bar = chr(92) + '|'
    lines = ['# DELFIN Viewer -- Bug Report', '']
    if meta.get('description'):
        lines += ['## What went wrong', '', meta['description'].strip(), '']
    else:
        lines += ['## What went wrong', '',
                  '_(no sentence was typed -- ask the reporter)_', '']

    lines += ['## Context', '', '| Field | Value |', '|---|---|']
    for key in ('created_at', 'user', 'host', 'tab', 'session_seconds',
                'atoms', 'delfin_version', 'python'):
        if meta.get(key) not in (None, ''):
            lines.append(f'| {key} | {str(meta[key]).replace("|", bar)} |')

    opening = next((e for e in timeline if e.get('k') == 'open'), None)
    was = (opening or {}).get('widgets') or {}
    now = meta.get('widgets_at_send') or {}
    if was or now:
        # Split, and not for tidiness. There are around sixty controls and in
        # an ordinary session four of them moved; run together, the four that
        # matter sit somewhere in the middle of a table nobody reads to the
        # end, and the one thing a reader came here for -- what was different
        # about this session -- is the hardest thing in it to find.
        moved = [name for name in sorted(set(was) | set(now))
                 if was.get(name, '') != now.get(name, '')]
        rest = [name for name in sorted(set(was) | set(now))
                if name not in moved]
        lines += ['', '## What was changed while it ran', '']
        if moved:
            lines += ['| Control | At the start | At send |', '|---|---|---|']
            for name in moved:
                lines.append(f'| {name} | `{was.get(name, "")}` | '
                             f'`{now.get(name, "")}` |')
        else:
            lines.append('Nothing: every control ended where it started.')
        lines += ['', '## Everything else, as it stood', '',
                  '| Control | Value |', '|---|---|']
        for name in rest:
            lines.append(f'| {name} | `{now.get(name, was.get(name, ""))}` |')

    lines += ['', '## What was recorded', '', '| Kind | Count |', '|---|---|']
    for kind, count in sorted((summary.get('by_kind') or {}).items()):
        lines.append(f'| {kind} | {count} |')
    lines += ['',
              f'{summary.get("events", 0)} events over '
              f'{summary.get("seconds", 0)} s, '
              f'{summary.get("bytes", 0)} bytes in memory.']
    if summary.get('dropped'):
        lines += ['',
                  f'**The session was longer than the buffer.** The oldest '
                  f'{summary["dropped"]} events ({summary["dropped_bytes"]} '
                  f'bytes) were dropped, so a replay starts part-way '
                  f'through: it begins from the earliest structure that '
                  f'survived, not from the one the session opened on. The '
                  f'bound is {summary.get("max_events")} events / '
                  f'{summary.get("max_bytes")} bytes.']

    began_on = opening_structure(timeline)
    if began_on:
        rows = began_on.splitlines()
        lines += ['', '## The structure it started from', '',
                  '```', '\n'.join(rows[:64])]
        if len(rows) > 64:
            lines.append('...')
        lines.append('```')

    shown = [e for e in timeline if e.get('k') != 'open']
    cut = shown[-tail:]
    told = ('The whole of it is `journal.jsonl.gz`, one JSON event per line, '
            'and `delfin.dashboard.editor_journal.replay` plays it back into '
            'a real editor.')
    if len(cut) < len(shown):
        told += f' The last {len(cut)} of {len(shown)} events are shown here.'
    lines += ['', '## The sequence', '', told, '', '```']
    lines += [_readable_event(e) for e in cut]
    lines += ['```', '']
    return '\n'.join(lines)


def write_report(journal, *, description='', widgets=None, tab='',
                 archive_dir=None, settings=None, extra=None):
    """Write one viewer bug report and return the directory it went to.

    Four files, and the split is the point: ``report.md`` is what a person
    reads, ``report.json`` is the header a machine reads first,
    ``journal.jsonl.gz`` is the sequence a machine replays, and
    ``structure_at_start.xyz`` is there so the molecule can be opened in
    anything without unpacking a report.

    ``report.md`` and ``report.json`` go through the agent's credential
    redactor on the way out, the way an agent report does -- they carry the
    sentence the user typed, and a sentence is the one part of this that a
    password can end up in.  The compressed stream does not: it is
    coordinates, run numbers and widget values written by the editor itself,
    and running a regex scan over several megabytes of it would cost seconds
    at the one moment the user is waiting for a button.
    """
    root = (Path(archive_dir).expanduser() if archive_dir
            else resolve_archive_dir(settings))
    user = _current_user()

    # Atomically, and with a fresh suffix on a clash: two reports in the same
    # second would otherwise share a name and the second would overwrite the
    # first. Same construction as the agent's reports next door.
    report_dir = None
    for attempt in range(12):
        name = _report_dirname(user, 'viewer')
        if attempt:
            name = f'{name}-{attempt}'
        candidate = root / name
        try:
            candidate.mkdir(parents=True, exist_ok=False)
            report_dir = candidate
            break
        except FileExistsError:
            continue
    if report_dir is None:
        import uuid
        report_dir = root / (_report_dirname(user, 'viewer') + '-'
                             + uuid.uuid4().hex)
        report_dir.mkdir(parents=True, exist_ok=True)

    timeline = journal.timeline()
    summary = journal.summary()
    start_xyz = opening_structure(timeline)
    atoms = ''
    if start_xyz:
        head = start_xyz.splitlines()[0].strip()
        atoms = head if head.isdigit() else ''

    meta = {
        'schema': 'delfin-viewer-bug-report/1',
        'created_at': datetime.now(timezone.utc).isoformat(),
        'user': user,
        'host': socket.gethostname(),
        'tab': tab,
        'description': description or '',
        'session_seconds': summary.get('seconds', 0),
        'atoms': atoms,
        'delfin_version': _delfin_version(),
        'python': sys.version.split()[0],
        'widgets_at_send': journal.snapshot(widgets or {}),
        'summary': summary,
    }
    if extra:
        meta['extra'] = extra

    # The stream. One JSON object per line, the way the agent's tool trace
    # and turn metrics are written -- greppable after a gunzip, streamable
    # while reading, and nothing has to be held whole to append to it.
    # Gzipped because the frame channel is most of it and coordinate JSON
    # compresses several-fold; the measured ratio for a real session is in
    # the test that writes one.
    body = '\n'.join(
        json.dumps({k: v for k, v in e.items() if k != '_n'},
                   ensure_ascii=False)
        for e in timeline) + '\n'
    with gzip.open(report_dir / 'journal.jsonl.gz', 'wt',
                   encoding='utf-8', compresslevel=6) as stream:
        stream.write(body)

    if start_xyz:
        (report_dir / 'structure_at_start.xyz').write_text(
            start_xyz if start_xyz.endswith('\n') else start_xyz + '\n',
            encoding='utf-8')

    redacted = _guarded_write(
        report_dir / 'report.json',
        json.dumps(meta, ensure_ascii=False, indent=2))
    redacted += _guarded_write(
        report_dir / 'report.md',
        _render_markdown(meta, timeline, summary))
    if redacted:
        try:
            (report_dir / 'REDACTIONS.txt').write_text(
                f'{redacted} credential-shaped value(s) were redacted from '
                'this report.\n', encoding='utf-8')
        except OSError:
            pass
    _make_group_readable(report_dir)
    return report_dir


def load_report(where):
    """Read a written report back as ``(meta, timeline)``.

    Takes the directory or either file in it, so a maintainer can pass
    whatever the shell completed.
    """
    path = Path(where).expanduser()
    if path.is_file():
        path = path.parent
    meta = json.loads((path / 'report.json').read_text(encoding='utf-8'))
    timeline = []
    with gzip.open(path / 'journal.jsonl.gz', 'rt', encoding='utf-8') as stream:
        for line in stream:
            line = line.strip()
            if line:
                timeline.append(json.loads(line))
    return meta, timeline


# ---------------------------------------------------------------------------
# Playing one back
# ---------------------------------------------------------------------------

def press(widget, value):
    """Set a widget the way the browser sets it, not the way the kernel does.

    Assigning ``widget.value`` from Python is the editor moving its own
    switch; a user clicking it arrives through ipywidgets' ``set_state``,
    which holds a property lock while it applies the change.  The journal
    reads that lock to tell the two apart, so anything that stands in for a
    hand -- a replay, and a headless test simulating a click -- has to hold
    the lock too, or a replayed press would come out as the editor pressing
    its own button.

    Falls back to a plain assignment on any ipywidgets that has no such lock,
    where the distinction cannot be made at all and losing it is better than
    losing the replay.

    The state is pushed to the front-end afterwards, which the real path
    deliberately does not do: a change that came *from* the browser must not
    be echoed back to it, but one that only pretends to have must, or a
    replay watched in a live dashboard would move the molecule while leaving
    the control that moved it showing its old value.  It is a no-op where
    there is no browser, which is where a replay usually runs.
    """
    lock = getattr(widget, '_lock_property', None)
    if lock is None:
        widget.value = value
        return
    with lock(value=value):
        widget.value = value
    try:
        widget.send_state('value')
    except Exception:
        pass


def page_messages(timeline, *, include_kernel=False):
    """The page-to-kernel half of a journal, in order, as comparable tuples.

    This is what a replay is judged by.  Everything else in a journal is the
    kernel's own answer -- frames, status lines, the coordinate box -- and
    those are what the replay is *for*: they come out of the editor again,
    and whether they come out the same is the question.  So the two sides
    are compared on what went in.

    A widget the *editor* moved is the kernel's answer too, and is left out
    for the same reason: a climb that converges after 39 steps one time and
    41 the next switches its own toggle off at a different moment, and a
    comparison that counted that as a difference would call every replay of
    a real calculation a failure.  Pass *include_kernel* to see them.
    """
    out = []
    for event in timeline:
        kind = event.get('k')
        if kind not in FROM_PAGE:
            continue
        if kind == 'w':
            if not include_kernel and event.get('by') != 'page':
                continue
            out.append(('w', event.get('n', ''),
                        json.dumps(event.get('v'), sort_keys=True)))
        elif kind == 'click':
            out.append(('click', event.get('n', ''), ''))
        else:
            out.append((kind, '', str(event.get('v', ''))))
    return out


def replay(timeline, editor, *, pace=1.0, max_gap=2.0, seed=True,
           sleep=None, on_event=None):
    """Drive a real editor through a recorded journal.

    *editor* is what :func:`delfin.dashboard.structure_editor.build` returns.
    The events the page sent are written back onto the same channels the page
    writes -- ``submit_cmd_sync``, ``submit_manip_sync``, ``submit_pick_sync``
    and the controls the user pressed -- and nothing else is written at all:
    the frames, the status lines, the coordinate box and every switch the
    editor moved for itself are the editor's own answers, and a replay that
    wrote those back would be reading its own homework.

    *pace* is how fast, as a multiple of the gaps that were recorded: 1.0
    replays at the speed the user worked, 0 as fast as the kernel will take
    it.  The gaps matter more than they look -- the follow is throttled, a
    scan paces its own frames, and a drag replayed instantly is not the drag
    that was recorded -- so the default is real time.  *max_gap* caps a
    single wait, so a report from a session where the user went to lunch does
    not replay the lunch.

    *seed* applies the pinned opening snapshot first, which is what puts the
    method, the charge and the toggles where they were before the first
    recorded change.

    Returns the list of messages it wrote, in the same shape
    :func:`page_messages` produces, so a caller can compare the two directly.
    """
    if sleep is None:
        sleep = time.sleep

    opening = next((e for e in timeline if e.get('k') == 'open'), None)
    if seed:
        # The molecule first. A journal whose beginning the ring dropped has
        # no opening snapshot left to trust, so the earliest structure that
        # survived stands in for it -- which is part of why the box is
        # recorded at all: without it a truncated report would replay a
        # gesture over whatever molecule happened to be loaded.
        began_on = opening_structure(timeline)
        box = getattr(editor, 'coords_widget', None)
        if began_on and box is not None:
            box.value = began_on
    if seed and opening:
        for name, value in sorted((opening.get('widgets') or {}).items()):
            thing = getattr(editor, name, None)
            if thing is None or not hasattr(thing, 'value'):
                continue
            try:
                thing.value = value
            except Exception:
                # A dropdown whose options depend on what is loaded may not
                # hold the recorded value yet. Skipped rather than raised:
                # one control that will not take its old value is not a
                # reason to abandon the replay of everything after it.
                continue

    written = []
    last = None
    for event in timeline:
        kind = event.get('k')
        if kind not in FROM_PAGE:
            continue
        now = float(event.get('t') or 0.0)
        if pace and last is not None:
            wait = (now - last) * float(pace)
            if wait > 0:
                sleep(min(wait, max_gap))
        last = now
        if kind == 'cmd':
            editor.submit_cmd_sync.value = str(event.get('v', ''))
            written.append(('cmd', '', str(event.get('v', ''))))
        elif kind == 'manip':
            editor.submit_manip_sync.value = str(event.get('v', ''))
            written.append(('manip', '', str(event.get('v', ''))))
        elif kind == 'pick':
            editor.submit_pick_sync.value = str(event.get('v', ''))
            written.append(('pick', '', str(event.get('v', ''))))
        elif kind == 'click':
            thing = getattr(editor, event.get('n', ''), None)
            if thing is None or not hasattr(thing, 'click'):
                continue
            thing.click()
            written.append(('click', event.get('n', ''), ''))
        elif kind == 'w':
            if event.get('by') != 'page':
                # The editor's own switch. It moves again by itself when the
                # replay reaches the same place, and moving it by hand here
                # would put a press in front of the thing that causes it.
                continue
            name = event.get('n', '')
            thing = getattr(editor, name, None)
            if thing is None or not hasattr(thing, 'value'):
                continue
            try:
                press(thing, event.get('v'))
            except Exception:
                continue
            written.append(('w', name,
                            json.dumps(event.get('v'), sort_keys=True)))
        if on_event is not None:
            on_event(event)
    return written


def list_reports(archive_dir=None):
    """Every report in the local archive, newest first.

    The directory names carry their UTC timestamp, so sorting them is
    sorting by time and nothing has to be opened to put the list in order.
    """
    root = (Path(archive_dir).expanduser() if archive_dir
            else resolve_archive_dir())
    if not root.is_dir():
        return []
    out = []
    for folder in root.iterdir():
        if not folder.is_dir():
            continue
        meta = {}
        header = folder / 'report.json'
        if header.is_file():
            try:
                meta = json.loads(header.read_text(encoding='utf-8'))
            except Exception:
                meta = {}
        out.append({
            'name': folder.name,
            'path': str(folder),
            'created_at': meta.get('created_at', ''),
            'user': meta.get('user', ''),
            'tab': meta.get('tab', ''),
            'description': meta.get('description', ''),
            'events': (meta.get('summary') or {}).get('events', 0),
        })
    out.sort(key=lambda row: row['name'], reverse=True)
    return out
