"""Scrolling up in one session, while another session is running.

Every open session keeps its own ``.delfin-agent-chat`` in the page; only
one of them has a viewport, which is how ``__delfinQ`` finds it. But each
of them carries the scroll tag, and the tag fires whenever THAT session's
chat is rebuilt -- which a background session does on every refresh.

The follow state was one object on ``window``, shared by all of them. A
hidden element measures ``scrollHeight`` and ``clientHeight`` as 0, so its
sync computed ``max = 0``, found the reader's offset above it, read that
as "the conversation was swapped" and turned following back on. The next
refresh of the session being read then pulled it to the bottom -- which is
the report: reading scrolled-up history, thrown to the end now and again,
with the interval set by some other session's output.

Driven rather than read: the behaviour is JavaScript in a Python string
that no import can call, so node runs the real script against two
stand-in elements.
"""

from __future__ import annotations

import ast
import json
import pathlib
import shutil
import subprocess
import tempfile

import pytest

_NODE = shutil.which("node")
_TAB_AGENT = (pathlib.Path(__file__).resolve().parents[1]
              / "delfin" / "dashboard" / "tab_agent.py")


def _chat_scroll_js() -> str:
    blocks = [
        node.value
        for node in ast.walk(ast.parse(_TAB_AGENT.read_text(encoding="utf-8")))
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and "__delfinChatScroll" in node.value
    ]
    assert len(blocks) == 1
    return blocks[0]


#: Two sessions, as the page holds them: one being read, one running behind
#: it. The hidden one reports zero for every measurement, which is what a
#: browser reports for an element whose ancestor is display:none.
_DRIVER = r"""
const docListeners = {};
let pendingScroll = [];

function makeChat(id, visible) {
  const jump = {hidden: true, textContent: ''};
  const c = {
    _top: 0,
    scrollHeight: visible ? 1000 : 0,
    clientHeight: visible ? 400 : 0,
    offsetParent: visible ? {} : null,
    classList: {contains: (k) => k === 'delfin-agent-chat'},
    getAttribute: (k) => (k === 'data-delfin-session' ? id : null),
    dataset: {delfinSession: id},
    querySelector: (s) => (s === '.delfin-chat-jump' ? jump : null),
  };
  c.closest = (s) => (s === '.delfin-agent-chat' ? c : null);
  jump.closest = c.closest;
  Object.defineProperty(c, 'scrollTop', {
    get() { return c._top; },
    set(v) {
      const max = Math.max(0, c.scrollHeight - c.clientHeight);
      const clamped = Math.max(0, Math.min(max, v));
      if (clamped === c._top) return;
      c._top = clamped;
      pendingScroll.push(c);
    },
  });
  return c;
}

const read = makeChat('session-being-read', true);
const behind = makeChat('session-running-behind', false);
const all = [read, behind];

global.window = {};
global.document = {
  querySelectorAll: () => all,
  querySelector(s) {
    if (s === '.delfin-agent-chat') return read;
    return null;
  },
  addEventListener(type, fn) {
    (docListeners[type] = docListeners[type] || []).push(fn);
  },
};
// The page's own visible-session lookup, as tab_agent defines it.
global.window.__delfinQ = function (sel) {
  const found = document.querySelectorAll(sel);
  for (const el of found) if (el.offsetParent !== null) return el;
  return found.length ? found[0] : null;
};
global.setInterval = () => 1;
global.setTimeout = () => 1;
global.clearTimeout = () => {};

__SCRIPT__

function flush() {
  const queued = pendingScroll; pendingScroll = [];
  for (const t of queued)
    for (const fn of (docListeners.scroll || [])) fn({target: t});
}
function userScroll(c, top) {
  const max = Math.max(0, c.scrollHeight - c.clientHeight);
  c._top = Math.max(0, Math.min(max, top));
  for (const fn of (docListeners.scroll || [])) fn({target: c});
}
const seen = {};
const endOf = (c) => c.scrollHeight - c.clientHeight;

// The reader scrolls up in the session they are looking at.
userScroll(read, 200);
seen['the reader is off the end'] = read.scrollTop === 200;

// The session running behind refreshes: its own scroll tag fires.
window.__delfinChatSync(behind); flush();

// ... and the session being read refreshes too, a moment later.
window.__delfinChatSync(read); flush();
seen['the reader is still where they were'] = read.scrollTop === 200;
seen['the reader was not pulled to the end'] = read.scrollTop !== endOf(read);

// Repeated background refreshes must not accumulate into a jump either.
for (let i = 0; i < 5; i++) { window.__delfinChatSync(behind); flush(); }
window.__delfinChatSync(read); flush();
seen['five background refreshes still leave it'] = read.scrollTop === 200;

// And following still works for a reader who IS at the end.
userScroll(read, endOf(read));
read.scrollHeight = 1400;
window.__delfinChatSync(read); flush();
seen['a reader at the end is still followed'] = read.scrollTop === endOf(read);

console.log(JSON.stringify(seen));
"""


@pytest.mark.skipif(_NODE is None, reason="node is not installed")
def test_a_background_session_does_not_move_the_one_being_read():
    script = _DRIVER.replace("__SCRIPT__", _chat_scroll_js())
    folder = pathlib.Path(tempfile.mkdtemp())
    path = folder / "two_sessions.js"
    path.write_text(script, encoding="utf-8")
    done = subprocess.run([_NODE, str(path)], capture_output=True, text=True,
                          timeout=300)
    assert done.returncode == 0, done.stderr[-2000:]
    seen = json.loads(done.stdout.strip().splitlines()[-1])
    wrong = sorted(name for name, ok in seen.items() if not ok)
    assert not wrong, f"with two sessions open: {', '.join(wrong)}"


def test_the_follow_state_is_not_one_object_for_every_session():
    """Read as well as driven: a single shared record is the defect itself,
    and a rewrite that keeps one would pass the driver only until the
    hidden element's measurements changed."""
    js = _chat_scroll_js()
    assert "data-delfin-session" in js or "delfinSession" in js, (
        "nothing in the script tells one session's chat from another's")


def test_every_chat_element_carries_its_session():
    """The state can only be per session if the element says which."""
    src = _TAB_AGENT.read_text(encoding="utf-8")
    assert src.count('<div class="delfin-agent-chat">') == 0, (
        "a chat element is still emitted without its session id")
