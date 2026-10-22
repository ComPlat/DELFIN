"""The message box returns to its resting size when a send empties it.

Input: the agent tab's auto-grow script and a message box that grew with
a multi-line message. Before: the box grew on 'input' events and was
re-measured once, 60 ms after a click in the send row. The send empties
the box from Python: the widget writes ``textarea.value`` without an
'input' event, after a kernel round trip that takes longer than 60 ms.
The box therefore kept its grown height until the next keystroke, with
Enter as with a click.
After: every write to the box's value re-measures it, whoever writes it;
an empty box drops its inline height and takes the stylesheet's 80 px;
a box that is not rendered (hidden tab) is not measured, since its
scrollHeight is 0.

Run in node over a page faked down to what the script touches, like the
viewer's player tests.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

_SRC = Path(__file__).resolve().parents[1] / "delfin/dashboard/tab_agent.py"

_PAGE_JS = r"""
const fs = require('fs');
const program = fs.readFileSync(process.argv[2], 'utf8');
const listeners = {};
let intervals = [];
const timers = [];
class HTMLTextAreaElement {
  constructor() { this._v = ''; this.style = {}; this.tagName = 'TEXTAREA';
                  this.offsetParent = {}; }
  // 20 px per line plus padding, as a browser would report it.
  get scrollHeight() { return this.offsetParent ? 24 + 20 * this._v.split('\n').length : 0; }
  closest(sel) { return sel === '.delfin-agent-input' ? {} : null; }
}
Object.defineProperty(HTMLTextAreaElement.prototype, 'value', {
  configurable: true,
  get() { return this._v; },
  set(v) { this._v = String(v); },
});
const ta = new HTMLTextAreaElement();
global.HTMLTextAreaElement = HTMLTextAreaElement;
global.window = { HTMLTextAreaElement };
global.document = {
  addEventListener(t, f) { (listeners[t] = listeners[t] || []).push(f); },
  querySelector(sel) { return sel === '.delfin-agent-input textarea' ? ta : null; },
};
global.setTimeout = (f, ms) => { timers.push(f); return timers.length; };
global.setInterval = (f, ms) => { intervals.push(f); return intervals.length; };
const fire = (t, target) => (listeners[t] || []).forEach(f => f({ target }));
const settle = () => { intervals.forEach(f => f()); while (timers.length) timers.shift()(); };
eval(program);
settle();
const out = {};
// The user types five lines: the box grows.
ta.value = 'a\nb\nc\nd\ne'; fire('input', ta); settle();
out.grown = ta.style.height;
// The send clears it from Python: a value write, no event, no keystroke.
ta.value = ''; settle();
out.after_send = ta.style.height;
// Typed again and cleared while the tab is hidden: nothing measured as 0.
ta.value = 'x\ny\nz'; fire('input', ta); settle();
ta.offsetParent = null; ta.value = 'x\ny\nz\nw'; settle();
out.hidden = ta.style.height;
process.stdout.write(JSON.stringify(out));
"""


def _grow_script() -> str:
    src = _SRC.read_text(encoding="utf-8")
    start = src.index("// --- Auto-grow the message box")
    end = src.index("})();", start) + len("})();")
    return src[start:end]


@pytest.mark.skipif(not shutil.which("node"), reason="node not installed")
def test_the_box_shrinks_when_a_send_empties_it(tmp_path):
    program = tmp_path / "grow.js"
    program.write_text(_grow_script(), encoding="utf-8")
    page = tmp_path / "page.js"
    page.write_text(_PAGE_JS, encoding="utf-8")
    done = subprocess.run(["node", str(page), str(program)],
                          capture_output=True, text=True, timeout=60)
    assert not done.returncode, done.stderr
    out = json.loads(done.stdout)
    assert out["grown"] == "124px"
    # Empty: the inline height goes, the stylesheet's 80 px applies.
    assert out["after_send"] == "", out
    # Hidden: the last real measurement stays; never 0px.
    assert out["hidden"] == "84px", out

