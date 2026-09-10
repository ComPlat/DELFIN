"""The chat window pulled itself to the newest output on every refresh.

A reader who scrolled up to re-read an earlier turn was thrown back to the
bottom within 150ms, so nothing above the last message could be read until
the agent had finished writing.

Following the newest output is now conditional on the reader being at the
end of the transcript, and there is a control to go back to the end.  The
behaviour lives in JavaScript inside a Python string, which no import can
call, so it is checked two ways: the emitted source is read for the shape of
the guard, and -- where node is installed -- the script itself is run against
a stand-in element and driven through the cases that matter.
"""

from __future__ import annotations

import ast
import json
import pathlib
import re
import shutil
import subprocess
import tempfile

import pytest

_NODE = shutil.which("node")
_TAB_AGENT = (pathlib.Path(__file__).resolve().parents[1]
              / "delfin" / "dashboard" / "tab_agent.py")


def _module_source() -> str:
    return _TAB_AGENT.read_text(encoding="utf-8")


def _chat_scroll_js() -> str:
    """The startup script the agent tab hands the browser.

    Taken from the parsed module rather than by string search, so the escapes
    in it are the ones the browser will see.
    """
    blocks = [
        node.value
        for node in ast.walk(ast.parse(_module_source()))
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and "__delfinChatScroll" in node.value
    ]
    assert len(blocks) == 1, (
        f"expected exactly one script to own the chat scroll, found "
        f"{len(blocks)}")
    return blocks[0]


def _poll_body() -> str:
    """The repeating callback that runs while the agent is writing."""
    js = _chat_scroll_js()
    at = js.rindex("setInterval(")
    return js[at:]


# --- what the emitted source has to say -----------------------------------


def test_the_poll_does_not_scroll_unconditionally():
    """The repeating callback used to assign the end of the content to
    scrollTop with nothing in front of it."""
    js = _chat_scroll_js()
    assert "chat.scrollTop = chat.scrollHeight" not in js, (
        "the poll still writes the end of the content straight into "
        "scrollTop, so a reader who scrolled up is pulled back")

    body = _poll_body()
    guard = body.find("S.follow")
    assert guard != -1, (
        "the poll does not consult the follow flag at all")
    write = body.find("setTop(")
    assert write != -1, "the poll no longer scrolls anything"
    assert guard < write, (
        "the poll scrolls before it checks whether the reader is at the "
        "end, so the check cannot stop it")


def test_no_render_path_forces_the_scroll():
    """Each refresh rebuilds the chat element and used to end it with a tag
    that drove the viewport to the bottom, whatever the reader was doing."""
    src = _module_source()
    assert "if(c)c.scrollTop=c.scrollHeight;" not in src, (
        "a render path still ends with an unconditional scroll to the "
        "bottom")
    assert "__delfinChatSync" in src, (
        "the render paths no longer hand the refresh to the follow-aware "
        "sync")


def test_the_jump_control_is_emitted():
    src = _module_source()
    assert 'class="delfin-chat-jump"' in src, (
        "the chat is rendered without a way back to the newest output")
    assert "window.__delfinChatToBottom(this)" in src, (
        "the jump control is rendered but is not wired to anything")
    assert ".delfin-chat-jump {" in src, (
        "the jump control has no styling of its own")
    # It floats over the transcript rather than sitting at its end, or a
    # reader far up the history could never reach it.
    style = src[src.index(".delfin-chat-jump {"):]
    style = style[:style.index("}")]
    assert "position: sticky" in style, (
        "the jump control scrolls away with the content instead of staying "
        "within reach")


def test_the_tolerance_is_named_and_real():
    """Reflow and late-loading images move the end of the content by a few
    pixels.  A tolerance of zero would drop a reader out of follow mode
    while they sat at the bottom, so the number is pinned here."""
    js = _chat_scroll_js()
    found = re.search(r"CHAT_BOTTOM_TOLERANCE_PX\s*=\s*(\d+)", js)
    assert found, "the at-the-end test names no tolerance"
    px = int(found.group(1))
    assert 40 <= px <= 80, (
        f"the tolerance is {px}px; below ~40 reflow alone ends follow mode, "
        f"above ~80 the reader is followed from well up the transcript")


def test_a_scroll_the_script_caused_is_told_apart_from_the_readers():
    js = _chat_scroll_js()
    assert "S.auto" in js, (
        "nothing marks the script's own scrollTop writes, so each of them "
        "looks like the reader scrolling away")


# --- what the script actually does ----------------------------------------


#: A stand-in for the parts of the page the script touches.  scrollTop is an
#: accessor because that is the whole difficulty: writing it queues a scroll
#: event exactly as the browser would, and the script has to tell that event
#: apart from one the reader caused.
_DRIVER = r"""
const docListeners = {};
let intervalFn = null;
let pendingScroll = [];

const jump = {hidden: true, textContent: '',
              closest: (s) => (s === '.delfin-agent-chat' ? chat : null)};
const chat = {
  scrollHeight: 1000,
  clientHeight: 400,
  _top: 0,
  classList: {contains: (c) => c === 'delfin-agent-chat'},
  querySelector: (s) => (s === '.delfin-chat-jump' ? jump : null),
  closest: (s) => (s === '.delfin-agent-chat' ? chat : null),
};
Object.defineProperty(chat, 'scrollTop', {
  get() { return chat._top; },
  set(v) {
    const max = Math.max(0, chat.scrollHeight - chat.clientHeight);
    const clamped = Math.max(0, Math.min(max, v));
    if (clamped === chat._top) return;
    chat._top = clamped;
    pendingScroll.push(chat);
  },
});
let working = null;
const sendBtn = {closest: (s) => (s === 'button' ? sendBtn : null)};

global.window = {};
global.document = {
  querySelector(s) {
    if (s === '.delfin-agent-chat') return chat;
    if (s === '.delfin-agent-working') return working;
    if (s === '.delfin-agent-send-row button') return sendBtn;
    return null;
  },
  addEventListener(type, fn) {
    (docListeners[type] = docListeners[type] || []).push(fn);
  },
};
global.setInterval = (fn) => { intervalFn = fn; return 1; };
global.setTimeout = () => 1;
global.clearTimeout = () => {};

__SCRIPT__

function flush() {
  const queued = pendingScroll; pendingScroll = [];
  for (const t of queued)
    for (const fn of (docListeners.scroll || [])) fn({target: t});
}
function userScroll(top) {
  const max = Math.max(0, chat.scrollHeight - chat.clientHeight);
  chat._top = Math.max(0, Math.min(max, top));
  for (const fn of (docListeners.scroll || [])) fn({target: chat});
}
function tick() { if (intervalFn) intervalFn(); flush(); }
function clickSend() {
  for (const fn of (docListeners.click || [])) fn({target: sendBtn});
  flush();
}
function sync() { window.__delfinChatSync(chat); flush(); }
const S = () => window.__delfinChatFollow;
const end = () => chat.scrollHeight - chat.clientHeight;
const seen = {};

// A reader at the end is carried along by new output.
sync();
seen.followed_to_the_end = chat.scrollTop === end();
seen.control_hidden_at_the_end = jump.hidden === true;

// The reader scrolls up to re-read something.
working = {working: true};
userScroll(100);
seen.scrolling_up_stops_following = S().follow === false;
seen.control_shown_when_away = jump.hidden === false;

// The agent keeps writing: the poll must leave the viewport alone.
chat.scrollHeight = 1600;
tick();
seen.poll_left_the_viewport = chat.scrollTop === 100;

// And so must a full refresh of the transcript.
chat.scrollHeight = 2200;
sync();
seen.refresh_left_the_viewport = chat.scrollTop === 100;
seen.control_reports_new_output = /New/.test(jump.textContent);

// Nearly at the end still counts as at the end; well above it does not.
userScroll(end() - 50);
seen.tolerance_holds_follow = S().follow === true;
userScroll(end() - 200);
seen.far_above_the_end_does_not = S().follow === false;

// The control returns to the end and re-enters follow mode.
window.__delfinChatToBottom(jump); flush();
seen.jump_reached_the_end = chat.scrollTop === end();
seen.jump_resumed_following = S().follow === true;
seen.jump_hid_the_control = jump.hidden === true;

// The script's own scroll must not read as the reader leaving.
chat.scrollHeight = 3000;
tick();
seen.own_scroll_kept_following = S().follow === true;
seen.poll_reached_the_new_end = chat.scrollTop === end();

// An offset past the end of what is now shown belongs to another
// conversation, so the reader is put at the end rather than mid-history.
userScroll(400);
chat.scrollHeight = 500;
sync();
seen.swapped_conversation_starts_at_the_end = S().follow === true;

// A cleared conversation drops the carried offset.
window.__delfinChatReset();
seen.reset_resumes_following = S().follow === true && S().top === 0;

// Sending is a deliberate move to the newest output, so it resumes
// following even from well up the transcript.
chat.scrollHeight = 2000;
userScroll(200);
clickSend();
seen.sending_returns_to_the_end =
  S().follow === true && chat.scrollTop === end();

console.log(JSON.stringify(seen));
"""


@pytest.mark.skipif(_NODE is None, reason="node not installed")
def test_the_script_follows_only_from_the_end():
    script = _DRIVER.replace("__SCRIPT__", _chat_scroll_js())
    folder = pathlib.Path(tempfile.mkdtemp())
    path = folder / "chat_scroll.js"
    path.write_text(script, encoding="utf-8")
    done = subprocess.run([_NODE, str(path)], capture_output=True, text=True,
                          timeout=300)
    assert done.returncode == 0, done.stderr[-2000:]
    seen = json.loads(done.stdout.strip().splitlines()[-1])
    wrong = sorted(name for name, ok in seen.items() if not ok)
    assert not wrong, f"the chat scroll misbehaved: {', '.join(wrong)}"
    assert len(seen) == 17, f"the driver checked {len(seen)} things, not 17"
