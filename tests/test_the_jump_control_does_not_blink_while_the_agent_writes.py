"""The "Newest" control blinked while the agent was writing.

Every refresh of the chat -- four a second while output streams -- replaces
the chat's content, and the control with it. It was rendered ``hidden`` and
shown a tick later by the follow script, so a reader who had scrolled up
watched it vanish with each new content and reappear right after. Its label
flickered the same way, from "Newest" to "New messages" and back.

The decision now lives on the host element, which survives a refresh, as two
marks the stylesheet reads: the button that is inserted carries no state of
its own and looks right in the frame it appears. Checked on the emitted
source, and -- where node is installed -- by driving the script through a
run of refreshes with a fresh button object each time, as the browser does.
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
    blocks = [
        node.value
        for node in ast.walk(ast.parse(_module_source()))
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and "__delfinChatScroll" in node.value
    ]
    assert len(blocks) == 1
    return blocks[0]


def _jump_markup() -> str:
    """The button as every render path emits it."""
    src = _module_source()
    at = src.index('class="delfin-chat-jump"')
    start = src.rindex("<button", 0, at)
    end = src.index("</button>", at)
    # The tag is built from adjacent string literals; collapse the quoting.
    raw = src[start:end]
    return re.sub(r"['\"]\s*\n\s*['\"]", "", raw)


def test_the_inserted_button_carries_no_state():
    markup = _jump_markup()
    assert " hidden" not in markup, (
        "the button is inserted hidden; the script shows it a tick later, "
        "and that gap is the blink on every refresh")
    assert "delfin-chat-jump-newest" in markup and \
        "delfin-chat-jump-unseen" in markup, (
            "both labels must be in the markup so neither needs writing in "
            "after the insert")


def test_the_stylesheet_reads_the_decision_off_the_host():
    src = _module_source()
    assert re.search(
        r'\.delfin-agent-chat-host:not\(\[data-delfin-follow="0"\]\)\s*'
        r'\.delfin-chat-jump', src), (
        "visibility is not decided by a mark on the host, so it cannot "
        "survive the button being replaced")
    assert '[data-delfin-unseen="1"] .delfin-chat-jump .delfin-chat-jump-unseen' \
        in src, "the label is not switched by a mark on the host"


def test_the_script_writes_nothing_to_the_button():
    """A write to the button is a write that is lost on the next refresh
    and has to be repeated after it -- which is the blink."""
    js = _chat_scroll_js()
    paint = js[js.index("function paint("):]
    paint = paint[:paint.index("\n        }")]
    assert "setAttribute('data-delfin-follow'" in paint
    assert "setAttribute('data-delfin-unseen'" in paint
    assert ".hidden" not in paint and ".textContent" not in paint, (
        "paint() still writes the button's own state")


_DRIVER = r"""
const docListeners = {};
let intervalFn = null;
let pendingScroll = [];
let buttonWrites = 0;

// A button the way the browser hands it over: a NEW object after every
// refresh, and any write to it is counted -- it would be lost with it.
function freshButton() {
  const b = {closest: (s) => (s === '.delfin-agent-chat-host' ? chat : null)};
  return new Proxy(b, {set(t, k, v) { buttonWrites++; t[k] = v; return true; }});
}
let jump = freshButton();
const chat = {
  scrollHeight: 1000,
  clientHeight: 400,
  _top: 0,
  attrs: {},
  setAttribute(k, v) { chat.attrs[k] = String(v); },
  classList: {contains: (c) => c === 'delfin-agent-chat-host'},
  querySelector: (s) => (s === '.delfin-chat-jump' ? jump : null),
  closest: (s) => (s === '.delfin-agent-chat-host' ? chat : null),
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
let working = {working: true};

global.window = {};
global.document = {
  querySelector(s) {
    if (s === '.delfin-agent-chat-host') return chat;
    if (s === '.delfin-agent-working') return working;
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
// A refresh: the content (and the button) is replaced, THEN the script
// runs. What the reader sees in between is decided by the host alone.
function refresh(grow) {
  chat.scrollHeight += grow;
  jump = freshButton();
  const between = chat.attrs['data-delfin-follow'] === '0';
  window.__delfinChatSync(chat); flush();
  return between;
}
const shown = () => chat.attrs['data-delfin-follow'] === '0';
const unseen = () => chat.attrs['data-delfin-unseen'] === '1';
const seen = {};

window.__delfinChatSync(chat); flush();
seen.hidden_at_the_end = !shown();

// The reader scrolls up; the control appears and stays through twelve
// refreshes -- also in the gap between each insert and the script's run.
userScroll(100);
seen.shown_once_away = shown();
let stayed = true, announced = true;
for (let i = 0; i < 12; i++) {
  if (!refresh(120)) stayed = false;
  if (!shown()) stayed = false;
  if (!unseen()) announced = false;
}
seen.control_never_vanished_between_refreshes = stayed;
seen.label_never_fell_back_to_newest = announced;
seen.nothing_was_written_to_the_button = buttonWrites === 0;

// Back at the end it disappears, and the label mark is cleared with it.
window.__delfinChatToBottom(jump); flush();
seen.hidden_again_at_the_end = !shown() && !unseen();
console.log(JSON.stringify(seen));
"""


@pytest.mark.skipif(_NODE is None, reason="node not installed")
def test_the_control_holds_still_through_a_run_of_refreshes():
    script = _DRIVER.replace("__SCRIPT__", _chat_scroll_js())
    folder = pathlib.Path(tempfile.mkdtemp())
    path = folder / "chat_jump.js"
    path.write_text(script, encoding="utf-8")
    done = subprocess.run([_NODE, str(path)], capture_output=True, text=True,
                          timeout=300)
    assert done.returncode == 0, done.stderr[-2000:]
    seen = json.loads(done.stdout.strip().splitlines()[-1])
    wrong = sorted(name for name, ok in seen.items() if not ok)
    assert not wrong, f"the control still blinks: {', '.join(wrong)}"
    assert len(seen) >= 6
