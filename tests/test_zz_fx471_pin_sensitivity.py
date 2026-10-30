"""Reviewer probe (nacht-s7), FX471: sensitivity of the js_sha256 pin.

The pin must change when ANY script that reaches the driver changes --
keys_js(), the lookup constant, or the agent keydown listener -- and must
match the fixture verbatim at the merged tip. Not a committed guard.
"""
import ast
import hashlib
import json
from pathlib import Path

from delfin.dashboard import header_ui
from delfin.dashboard.header_ui import keys_js

FIXTURE = json.loads(
    (Path(__file__).resolve().parent / 'fixtures/header_shortcuts/verdicts.json').read_text())


def _sources():
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


def test_pins_match_the_recorded_fixture():
    lookup, listeners = _sources()
    for variant in (0, 1):
        for header_first in (False, True):
            scripts = [keys_js(), lookup + listeners[variant]]
            if not header_first:
                scripts.reverse()
            digest = hashlib.sha256(
                '\n'.join(scripts).encode('utf-8')).hexdigest()
            key = 'variant%d_headerfirst%d' % (variant, header_first)
            assert FIXTURE[key]['js_sha256'] == digest, key


def test_pin_moves_on_a_lookup_change():
    lookup, listeners = _sources()
    scripts = [keys_js(), lookup + listeners[0]]
    base = hashlib.sha256('\n'.join(scripts).encode('utf-8')).hexdigest()
    marker = 'var all = document.querySelectorAll(sel)'
    assert marker in lookup
    changed = hashlib.sha256(
        '\n'.join([keys_js(), lookup.replace(marker, 'var all2 = 1')
                   + listeners[0]]).encode('utf-8')).hexdigest()
    assert changed != base


def test_pin_moves_on_a_keys_js_change():
    lookup, listeners = _sources()
    scripts = [keys_js(), lookup + listeners[0]]
    base = hashlib.sha256('\n'.join(scripts).encode('utf-8')).hexdigest()
    original = header_ui.keys_js
    header_ui.keys_js = lambda: original().replace('Escape', 'Esc')
    try:
        scripts2 = [header_ui.keys_js(), lookup + listeners[0]]
        changed = hashlib.sha256('\n'.join(scripts2).encode('utf-8')).hexdigest()
    finally:
        header_ui.keys_js = original
    assert changed != base


def test_pin_moves_on_a_listener_change():
    lookup, listeners = _sources()
    scripts = [keys_js(), lookup + listeners[0]]
    base = hashlib.sha256('\n'.join(scripts).encode('utf-8')).hexdigest()
    tweaked = listeners[0].replace('__delfinAgentKeys', '__delfinAgentKeysX')
    changed = hashlib.sha256(
        '\n'.join([keys_js(), lookup + tweaked]).encode('utf-8')).hexdigest()
    assert changed != base
