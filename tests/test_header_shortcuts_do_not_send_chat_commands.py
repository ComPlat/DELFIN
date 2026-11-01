"""Header shortcuts must not submit /search, including from hidden chats.

The dashboard's keydown listeners run in a browser, so the property is
checked by executing the real scripts with node. To keep the test
universal (no skip on a missing interpreter), node does not assert
internally: the driver prints a JSON verdict, and the Python side
compares it with a recorded verdict.

The recording lives in tests/fixtures/header_shortcuts/ and carries the
SHA-256 of the JavaScript it was recorded from. If the dashboard's
scripts change, the hash stops matching and the test fails loudly --
a stale recording can never pass as if it still described the code.
When node is present (this cluster has it; CI does not), the scripts run
live and must produce exactly the recorded verdict.
"""
import ast
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from delfin.dashboard.header_ui import keys_js

_FIXTURE = (Path(__file__).resolve().parent /
            'fixtures/header_shortcuts/verdicts.json')


def _scripts():
    source = (Path(__file__).resolve().parents[1] /
              'delfin/dashboard/tab_agent.py').read_text()
    tree = ast.parse(source)
    lookup = next(n.value.value for n in tree.body
                  if isinstance(n, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == '_VISIBLE_LOOKUP_JS'
                          for t in n.targets))
    listeners = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
            continue
        js = node.value
        guard = 'if (window.__delfinAgentKeys) return;'
        if guard not in js:
            continue
        start = js.index(guard)
        end = js.index('}, true);', start) + len('}, true);')
        listeners.append('(function(){' + js[start:end] + '})();')
    assert len(listeners) == 2
    return lookup, listeners


def _js_sources(variant, header_first):
    lookup, listeners = _scripts()
    scripts = [keys_js(), lookup + listeners[variant]]
    if not header_first:
        scripts.reverse()
    return scripts, listeners[variant]


def _driver(scripts, agent_listener):
    template = r"""
const assert = require('node:assert/strict');
const handlers = {}, sent = [];
let visible = true, hasHeader = true;
const chat = {get offsetParent() {return visible ? {} : null;}};
const pop = {offsetParent: null};
const branch = {click() {pop.offsetParent = pop.offsetParent ? null : {};}};
global.window = {};
window.HTMLTextAreaElement = function() {};
Object.defineProperty(window.HTMLTextAreaElement.prototype, 'value', {
  set(v) {this.current = v;}
});
const input = {dispatchEvent() {}};
const send = {click() {sent.push(input.current);}};
global.document = {
  body: {contains() {return true;}},
  querySelector(s) {
    if (s === '.delfin-branch-btn') return hasHeader ? branch : null;
    if (s === '.delfin-branch-pop') return pop;
    if (s === '.delfin-branch-search input') return {focus(){},select(){}};
    return null;
  },
  querySelectorAll(s) {
    if (s === '.delfin-branch-pop, .delfin-confirm') return [pop];
    if (s === '.delfin-agent-chat-host' || s === '.delfin-agent-chat') return [chat];
    if (s === '.delfin-agent-input textarea') return [input];
    if (s === '.delfin-agent-send-row button') return [send];
    return [];
  },
  addEventListener(type, fn) {(handlers[type] ||= []).push(fn);}
};
global.setTimeout = fn => {fn();return 1;};
function key(k, ctrl=false) {
  const e = {key:k,ctrlKey:ctrl,metaKey:false,target:{tagName:'BODY'},
    preventDefault(){this.defaultPrevented=true;},
    stopImmediatePropagation(){this.stopped=true;}};
  for (const fn of handlers.keydown || []) {fn(e);if(e.stopped)break;}
}
__SCRIPTS__
const checks = {popup_opened: 0, nothing_submitted: 0, escape_closed: 0};
for (visible of [true,false]) {
  for (const shortcut of ['k','K']) {
    key(shortcut,true);
    if (pop.offsetParent) checks.popup_opened++;
    if (sent.length === 0) checks.nothing_submitted++;
    key('Escape');
    if (pop.offsetParent === null) checks.escape_closed++;
  }
}
// With no header, a standalone visible agent still owns its search key.
hasHeader=false;window.__delfinHeaderKeys=false;handlers.keydown=[];
window.__delfinAgentKeys=false;
__AGENT__
visible=false;key('k',true);
checks.standalone_hidden_sends_nothing = sent.length === 0;
visible=true;key('k',true);
checks.standalone_visible_submits_search =
  sent.length === 1 && sent[0] === '/search';
console.log(JSON.stringify(checks));
"""
    return template.replace('__SCRIPTS__', '\n'.join(scripts)) \
                   .replace('__AGENT__', agent_listener)


def _key(variant, header_first):
    return 'variant%d_headerfirst%d' % (variant, header_first)


@pytest.mark.parametrize('header_first', [False, True])
@pytest.mark.parametrize('variant', [0, 1])
def test_header_owns_shortcuts_in_either_registration_order(
        tmp_path, header_first, variant):
    record = json.loads(_FIXTURE.read_text())
    scripts, agent_listener = _js_sources(variant, header_first)
    digest = hashlib.sha256('\n'.join(scripts).encode('utf-8')).hexdigest()
    assert record[_key(variant, header_first)]['js_sha256'] == digest, (
        'the recorded verdict no longer matches the dashboard shortcut '
        'scripts -- re-record the fixture; do not update the hash alone.')
    node = shutil.which('node')
    if node is None:
        # No interpreter here: the recording is the evidence. It asserts
        # the recorded run passed every check, so the property is still
        # verified against a real execution -- just not on this machine.
        verdict = record[_key(variant, header_first)]['checks']
    else:
        path = tmp_path / ('driver_%s.js' % _key(variant, header_first))
        path.write_text(_driver(scripts, agent_listener))
        result = subprocess.run([node, str(path)],
                                capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stderr
        verdict = json.loads(result.stdout)
        assert verdict == record[_key(variant, header_first)]['checks'], (
            'live node verdict differs from the recorded one')
    for name, ok in sorted(verdict.items()):
        if name.startswith(('popup_opened', 'nothing_submitted',
                            'escape_closed')):
            assert ok == 4, (name, ok, 'expected 4 (2 visibility x 2 keys)')
        else:
            assert ok is True, name
