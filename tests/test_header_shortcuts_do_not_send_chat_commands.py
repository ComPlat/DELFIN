"""Header shortcuts must not submit /search, including from hidden chats."""
import ast
from pathlib import Path
import shutil
import subprocess

import pytest

from delfin.dashboard.header_ui import keys_js


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


@pytest.mark.skipif(shutil.which('node') is None, reason='node is not installed')
@pytest.mark.parametrize('header_first', [False, True])
@pytest.mark.parametrize('variant', [0, 1])
def test_header_owns_shortcuts_in_either_registration_order(tmp_path, header_first, variant):
    lookup, listeners = _scripts()
    driver = r"""
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
for (visible of [true,false]) {
  for (const shortcut of ['k','K']) {
    key(shortcut,true);
    assert.ok(pop.offsetParent, 'header did not open');
    assert.deepEqual(sent, [], 'shortcut submitted a chat command');
    key('Escape');
    assert.equal(pop.offsetParent,null,'Escape did not close');
  }
}
// With no header, a standalone visible agent still owns its search key.
hasHeader=false;window.__delfinHeaderKeys=false;handlers.keydown=[];
window.__delfinAgentKeys=false;
__AGENT__
visible=false;key('k',true);assert.deepEqual(sent,[]);
visible=true;key('k',true);assert.deepEqual(sent,['/search']);
"""
    scripts = [keys_js(), lookup + listeners[variant]]
    if not header_first:
        scripts.reverse()
    driver = driver.replace('__SCRIPTS__', '\n'.join(scripts))
    driver = driver.replace('__AGENT__', listeners[variant])
    path = tmp_path / 'shortcuts.js'
    path.write_text(driver)
    result = subprocess.run([shutil.which('node'), str(path)],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
