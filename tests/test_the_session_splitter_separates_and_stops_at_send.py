"""The session splitter drags, in both directions, and ends at Send.

The branch that built it asserted the MARKUP and that the script was
registered -- not that dragging does anything. A predicate passing is
not the gate passing, so this file runs the real handler against a DOM
and measures what it writes.

Three properties, each a thing the user reported:

* it separates dynamically -- a drag to the right widens the column and a
  drag to the LEFT narrows it again (the second commit on the branch is
  called "resizing in both directions", which is the bug this pins);
* the gap is symmetric -- it used to be the sidebar's `margin-right:
  12px`, so there were 12px on one side of the handle and none on the
  other, and the bar sat flush against the chat;
* it reaches the Send row and stops, instead of running to the bottom of
  the page like a document divider.
"""

from __future__ import annotations

import ast
import json
import pathlib
import shutil
import subprocess
import tempfile

import pytest

from delfin.dashboard import agent_sessions as AS

_NODE = shutil.which("node")


def _splitter_js() -> str:
    """The handler as it ships, pulled from the module's own constant."""
    return AS._SPLITTER_INIT_JS


_DRIVER = r"""
'use strict';
const seen = {};

function rect(o) { return () => o; }

// --- the fake page -------------------------------------------------------
const sidebar = {style: {}, className: 'delfin-sessions'};
// The shell carries the custom property the stylesheet reads.
const props = {};
const splitter = {
  style: {}, className: 'delfin-session-splitter',
  _listeners: {},
  addEventListener(t, f) { (this._listeners[t] = this._listeners[t] || []).push(f); },
  setPointerCapture() { seen.capture_taken = true; },
  releasePointerCapture() { seen.capture_released = true; },
  getBoundingClientRect: rect({top: 100, bottom: 600, left: 250, right: 258}),
};
const sendRow = {getBoundingClientRect: rect({top: 540, bottom: 580})};
const host = {};
const shell = {
  dataset: {},
  style: {
    setProperty(k, v) { props[k] = v; },
    getPropertyValue(k) { return props[k] || ''; },
  },
  getBoundingClientRect: rect({top: 100, bottom: 1400, left: 0, right: 1200}),
  querySelector(sel) {
    if (sel === '.delfin-session-splitter') return splitter;
    if (sel === '.delfin-sessions') return sidebar;
    if (sel === '.delfin-splitter-host') return host;
    if (sel === '.delfin-agent-send-row') return sendRow;
    return null;
  },
};
host.parentElement = shell;

const docL = {};
// The shell appears only after `present` flips -- the real case: the
// agent tab is not the one in front when the init script runs, so its
// DOM does not exist yet.
let present = false;
const timers = [];
global.setTimeout = (fn) => { timers.push(fn); return timers.length; };
function runTimers(n) {
  for (let i = 0; i < n && timers.length; i++) (timers.shift())();
}
global.document = {
  body: {},
  querySelectorAll: (sel) =>
    (sel === '.delfin-session-shell' && present) ? [shell] : [],
  addEventListener(t, f) { (docL[t] = docL[t] || []).push(f); },
  removeEventListener(t, f) {
    docL[t] = (docL[t] || []).filter((g) => g !== f);
  },
};
global.window = {addEventListener() { seen.resize_hooked = true; }};
global.ResizeObserver = function (fn) {
  this.observe = () => { seen.observer_attached = true; };
};
global.MutationObserver = function (fn) {
  this.observe = () => { seen.mutation_watch_attached = true; };
};

__SCRIPT__

// --- drive it ------------------------------------------------------------
function width() {
  const v = props['--delfin-sessions-w'];
  return v === undefined ? null : parseFloat(v);
}
function down() {
  // Named rather than crashing: with no retry, nothing is ever bound and
  // the array is undefined. A TypeError here tells the next reader
  // nothing, so the property says what was missing.
  const fns = splitter._listeners.pointerdown;
  if (!fns || !fns.length) {
    seen.a_handler_is_attached_at_all = false;
    console.log(JSON.stringify(seen));
    process.exit(0);
  }
  seen.a_handler_is_attached_at_all = true;
  fns.forEach((f) => f({preventDefault() {}, pointerId: 1}));
}
function drag(x) {
  down();
  (docL.pointermove || []).forEach((f) => f({clientX: x}));
  (docL.pointerup || []).forEach((f) => f({pointerId: 1}));
}

// Nothing to bind at first run, and no handler may have been attached.
seen.nothing_bound_while_absent =
  shell.dataset.delfinSplitterBound === undefined
  && (splitter._listeners.pointerdown || []).length === 0;

// The tab is opened: the shell exists now, and a later retry finds it.
present = true;
runTimers(5);
seen.bound_after_the_tab_appears =
  shell.dataset.delfinSplitterBound === '1'
  && (splitter._listeners.pointerdown || []).length === 1;

seen.bound_once = shell.dataset.delfinSplitterBound === '1';
seen.starts_unstyled = width() === null;

drag(320);
seen.drag_right_widens = width() === 320;

// Dynamic: every move writes, not just the last one before release.
down();
const seenWidths = [];
[300, 280, 260, 240].forEach((x) => {
  (docL.pointermove || []).forEach((f) => f({clientX: x}));
  seenWidths.push(width());
});
(docL.pointerup || []).forEach((f) => f({pointerId: 1}));
seen.每_move_writes = false;   // replaced below
seen.every_move_writes =
  JSON.stringify(seenWidths) === JSON.stringify([300, 280, 260, 240]);
delete seen['每_move_writes'];

drag(200);
seen.drag_left_narrows = width() === 200;
seen.both_directions = true;        // reached only if the two above ran

drag(20);
seen.clamped_at_minimum = width() === 170;

drag(5000);
seen.clamped_at_maximum = width() === 420;

// Nothing may be written as a plain inline `flex`: the stylesheet
// declares it !important, which beats an inline declaration, so such a
// write is silently discarded by the browser.
seen.no_inline_flex_written = sidebar.style.flex === undefined;

// The height was measured to the Send row's bottom, not the page bottom.
seen.height_stops_at_send = splitter.style.height === '480px';

// Listeners are taken off the document again, or every later move drags.
seen.move_listener_removed = (docL.pointermove || []).length === 0;
seen.up_listener_removed = (docL.pointerup || []).length === 0;

console.log(JSON.stringify(seen));
"""


@pytest.mark.skipif(_NODE is None, reason="node not installed")
def test_the_splitter_drags_both_ways_and_stops_at_send():
    script = _DRIVER.replace("__SCRIPT__", _splitter_js())
    folder = pathlib.Path(tempfile.mkdtemp())
    path = folder / "splitter.js"
    path.write_text(script, encoding="utf-8")
    done = subprocess.run([_NODE, str(path)], capture_output=True, text=True,
                          timeout=300)
    assert done.returncode == 0, done.stderr[-2000:]
    seen = json.loads(done.stdout.strip().splitlines()[-1])
    wrong = sorted(name for name, ok in seen.items() if not ok)
    assert not wrong, f"the splitter misbehaved: {', '.join(wrong)}"
    assert len(seen) >= 17, f"the driver checked only {len(seen)} things"


def test_the_gap_is_one_number_on_both_sides():
    """Symmetric by construction. The sidebar must carry no one-sided
    margin -- that is what put the whole gap on one side and left the
    handle flush against the chat."""
    css = AS._SIDEBAR_CSS
    i = css.index(".delfin-sessions {")
    block = css[i:css.index("}", i)]
    assert "margin: 0 12px 0 0" not in block, (
        "the one-sided margin is back; the handle will sit flush again")
    assert "margin: 0;" in block

    j = css.index(".delfin-splitter-host {")
    host = css[j:css.index("}", j)]
    assert "margin: 0 var(--delfin-splitter-gap" in host, (
        "the gap must be ONE value applied to both sides, so no second "
        "number can drift away from the first")


def test_the_height_is_not_the_whole_column():
    """It used to be `flex: 1 1 auto` inside a stretched host, which is
    the full height of the page. Both the CSS fallback and the measured
    value have to end above that."""
    css = AS._SIDEBAR_CSS
    i = css.index(".delfin-session-splitter {")
    block = css[i:css.index("}", i)]
    assert "flex: 1 1 auto" not in block
    assert "height:" in block and "100vh -" in block, (
        "a fallback height is needed for the first paint, before the "
        "measure runs")
    assert "min-height: 160px" in block

    js = AS._SPLITTER_INIT_JS
    assert ".delfin-agent-send-row" in js, "nothing measures the composer"
    assert "ResizeObserver" in js and "resize" in js, (
        "a one-shot measure is right for one layout and wrong after the "
        "first resize")


def test_the_width_is_driven_through_the_property_the_rule_reads():
    """The cascade trap, pinned.

    The sidebar's width is declared `!important`, because a widget
    stylesheet overrides a plain declaration. A stylesheet `!important`
    also beats an INLINE declaration -- so the handler as first written
    set `sidebar.style.flex` on every pointermove and the browser
    discarded all of it. The drag looked dead on the real page while a
    DOM test that models no cascade reported it working in both
    directions.

    So: the width travels through one custom property that the
    `!important` rule itself reads.
    """
    css = AS._SIDEBAR_CSS
    i = css.index(".delfin-sessions {")
    block = css[i:css.index("}", i)]
    assert "var(--delfin-sessions-w" in block, (
        "the !important rule must read the property the drag writes")
    assert "flex: 0 0 236px !important" not in block, (
        "a hardcoded basis cannot be overridden from script")

    js = AS._SPLITTER_INIT_JS
    assert "setProperty('--delfin-sessions-w'" in js
    # The ASSIGNMENT, not the mention: the comment in the handler explains
    # why this write is wrong, and a check for the bare name matched its
    # own explanation.
    import re as _re
    assert not _re.search(r"sidebar\.style\.(flex|width|maxWidth)\s*=", js), (
        "an inline width write is discarded under the !important rule")


# --- the jitter: a resize loop that fed itself ---------------------------

_OSCILLATE = r"""
'use strict';
const seen = {};
let writes = 0;
let bottom = 580;          // the send row's bottom edge
let oscillate = false;     // when true it moves by 1px each measurement
const timers = [];

const sendRow = {
  getBoundingClientRect() {
    if (oscillate) bottom = (bottom === 580) ? 581 : 580;
    return {top: 540, bottom: bottom};
  },
};
const splitter = {
  _h: '', _listeners: {},
  style: {
    set height(v) {
      writes += 1;
      splitter._h = v;
      // The feedback the real page had: writing the child's height
      // re-measures the container, which fires the observer again.
      (observers || []).forEach((cb) => cb());
    },
    get height() { return splitter._h; },
  },
  addEventListener(t, f) { (this._listeners[t] = this._listeners[t] || []).push(f); },
  setPointerCapture() {}, releasePointerCapture() {},
  getBoundingClientRect: () => ({top: 100, bottom: 600, left: 250, right: 258}),
};
const sidebar = {style: {}};
const host = {};
const props = {};
const shell = {
  dataset: {},
  style: {setProperty(k, v) { props[k] = v; },
          getPropertyValue(k) { return props[k] || ''; }},
  getBoundingClientRect: () => ({top: 100, bottom: 1400, left: 0, right: 1200}),
  querySelector(sel) {
    if (sel === '.delfin-session-splitter') return splitter;
    if (sel === '.delfin-sessions') return sidebar;
    if (sel === '.delfin-splitter-host') return host;
    if (sel === '.delfin-agent-send-row') return sendRow;
    return null;
  },
};
host.parentElement = shell;

let observers = [];
global.ResizeObserver = function (cb) { this.observe = () => { observers.push(cb); }; };
global.MutationObserver = function () { this.observe = () => {}; };
global.requestAnimationFrame = (fn) => { timers.push(fn); };
global.setTimeout = (fn) => { timers.push(fn); return timers.length; };
function flush(n) { for (let i = 0; i < n && timers.length; i++) (timers.shift())(); }

const docL = {};
let present = true;
global.document = {
  body: {},
  querySelectorAll: (s) => (s === '.delfin-session-shell' && present) ? [shell] : [],
  addEventListener(t, f) { (docL[t] = docL[t] || []).push(f); },
  removeEventListener(t, f) { docL[t] = (docL[t] || []).filter((g) => g !== f); },
};
global.window = {addEventListener() {}};

__SCRIPT__

const afterInit = writes;
seen.init_writes_once = afterInit <= 1;

// A sub-pixel layout that moves by one pixel on every measurement. The
// first version wrote on every callback and never settled.
oscillate = true;
observers.forEach((cb) => cb());
flush(500);
seen.one_pixel_jitter_settles = (writes - afterInit) === 0;
seen.not_runaway = writes < 5;

// A real change still gets through.
oscillate = false;
bottom = 700;
observers.forEach((cb) => cb());
flush(50);
seen.a_real_change_is_applied = splitter._h === '600px';

// A burst in one frame collapses to a single write.
const before = writes;
bottom = 900;
for (let i = 0; i < 50; i++) observers.forEach((cb) => cb());
flush(200);
seen.a_burst_writes_once = (writes - before) === 1;

console.log(JSON.stringify(seen));
"""


@pytest.mark.skipif(_NODE is None, reason="node not installed")
def test_the_height_fit_does_not_oscillate():
    """Reported from a real session: the dashboard window jittered so
    badly the chat could not be read.

    `fit()` wrote the height of a CHILD of `shell` while the
    ResizeObserver watched `shell` itself, so every write re-measured the
    container and fired the observer again. With sub-pixel layout the
    measurement alternates by a pixel and the loop never settles.

    Driven here rather than reasoned about: the fake observer fires on
    every write, exactly as the page did, and the measured edge moves by
    one pixel each time it is read.
    """
    script = _OSCILLATE.replace("__SCRIPT__", _splitter_js())
    folder = pathlib.Path(tempfile.mkdtemp())
    path = folder / "oscillate.js"
    path.write_text(script, encoding="utf-8")
    done = subprocess.run([_NODE, str(path)], capture_output=True, text=True,
                          timeout=300)
    assert done.returncode == 0, done.stderr[-2000:]
    seen = json.loads(done.stdout.strip().splitlines()[-1])
    wrong = sorted(name for name, ok in seen.items() if not ok)
    assert not wrong, f"the height fit still oscillates: {', '.join(wrong)}"


def test_the_fit_watches_the_row_and_not_its_own_container():
    """The loop closed because the observer watched the element that
    contains the splitter. It must watch what DECIDES the height."""
    from delfin.dashboard import agent_sessions as AS

    js = AS._SPLITTER_INIT_JS
    assert "ResizeObserver(schedule)" in js, (
        "the observer must go through the per-frame scheduler, not call "
        "fit() directly")
    assert ".observe(watched || shell)" in js, (
        "it must watch the send row, falling back to the shell only when "
        "there is no row to watch")
    assert "Math.abs(h - lastH)" in js, (
        "equality alone does not stop an alternation between two values "
        "one pixel apart")
