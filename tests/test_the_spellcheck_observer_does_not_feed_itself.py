"""The spellcheck observer never writes an attribute it watches unguarded.

disable_spellcheck installs a MutationObserver on the attributes
spellcheck/autocorrect/autocapitalize/autocomplete and answered each
change by setting them again. setAttribute records a mutation even for the
value already there, so a single change fed the observer its own writes
forever: the page froze. Measured in a real browser: the Calculations tab
froze on its first open, and the DELFIN Agent tab with it (it builds the
Calculations tab first). The fix writes an attribute only when it differs.
"""
from __future__ import annotations

import re
import types

from delfin.dashboard import helpers


def _script(monkeypatch):
    seen = []
    monkeypatch.setattr(helpers, "_append_js", lambda ctx, js: seen.append(js))
    helpers.disable_spellcheck(types.SimpleNamespace(), "delfin-nospell")
    assert seen, "disable_spellcheck emitted no script"
    return seen[-1]


def test_no_watched_attribute_is_written_unguarded(monkeypatch):
    js = _script(monkeypatch)
    watched = re.search(r"attributeFilter:\s*\[([^\]]*)\]", js)
    assert watched, "the observer no longer names its attributes"
    names = re.findall(r"'([^']+)'", watched.group(1))
    assert names
    for name in names:
        assert f"el.setAttribute('{name}'" not in js, name
        assert f"put('{name}'" in js, name
    assert "el.spellcheck = false" not in js


def test_put_compares_before_it_writes(monkeypatch):
    js = _script(monkeypatch)
    assert re.search(r"if \(el\.getAttribute\(name\) !== value\) el\.setAttribute\(name, value\)", js)
