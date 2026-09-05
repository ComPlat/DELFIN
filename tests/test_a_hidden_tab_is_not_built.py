"""Hiding a tab in Settings dropped it from the bar and nothing else.

It was still built: the Python widgets, the ipywidgets models, the DOM, the
listeners, all of it, for a tab the user had said they did not want. Archive
and Office are the ones that matter, because each is another whole build of
the calculations browser -- the heaviest thing a cold start does. Between them
they carry about 1500 DOM elements and 640 JS listeners.

So a hidden tab is not built. A visible one still is, at startup, exactly as
before -- measured with nothing hidden, which is the default: 6.64 s to the
first widget, 4 long tasks / 1235 ms, 7601 elements, 2024 listeners, 14 tabs.
Unchanged.

Showing a tab again has to bring it back, and a spec whose widget is None is
dropped from the bar, so it is built at the moment it becomes visible.

Not covered here: the hidden case in a live browser. Settings live in
``~/.delfin_settings.json`` with no override, and proving it end to end would
mean editing the user's own file.
"""

import pathlib
import re

import pytest

from delfin import dashboard as dashboard_module


@pytest.fixture(scope='module')
def source():
    return pathlib.Path(dashboard_module.__file__).read_text(encoding='utf-8')


def test_a_hidden_clone_is_not_built(source):
    """The decision itself, exercised rather than read."""
    built = []

    def make(name):
        def build(_ctx):
            built.append(name)
            return f'<{name}>', {'root': name}
        return build

    def clone(tab_id, build, refs, hidden):
        # the shape of _clone_of_browser
        if tab_id in hidden:
            return None
        widget, made = build(None)
        refs.update(made or {})
        return widget

    refs_a, refs_o = {}, {}
    assert clone('archive', make('archive'), refs_a, set()) == '<archive>'
    assert clone('office', make('office'), refs_o, set()) == '<office>'
    assert built == ['archive', 'office']
    assert refs_a and refs_o, 'Settings needs its handle on a built tab'

    built.clear()
    refs_a, refs_o = {}, {}
    assert clone('archive', make('archive'), refs_a, {'archive', 'office'}) is None
    assert clone('office', make('office'), refs_o, {'archive', 'office'}) is None
    assert built == [], 'a hidden tab must cost nothing'


def test_the_hidden_set_comes_from_settings(source):
    assert "_hidden_at_start" in source
    assert "get('hidden', [])" in source
    # Read before the clones are built, or the decision comes too late.
    assert source.index('_hidden_at_start = {') < source.index(
        "tab6 = _clone_of_browser('archive'")


def test_showing_a_tab_again_builds_it(source):
    rebuild = source.split('def _rebuild_dashboard_tabs')[1].split('\n    def ')[0]
    assert "spec.get('widget') is None and callable(spec.get('build'))" in rebuild
    assert '_build_tab_now(spec)' in rebuild

    for spec_id in ("'id': 'archive'", "'id': 'office'"):
        block = source.split(spec_id)[1].split('},')[0]
        assert "'build':" in block, f'{spec_id} could not come back'


def test_a_late_tab_does_not_wipe_the_page(source):
    """Every tab's startup script goes out in one ctx.run_js, and run_js clears
    the output before it writes. A tab built after that must append."""
    builder = source.split('def _build_tab_now')[1].split('\n    def ')[0]
    # Its own docstring names run_js as the thing it must not use, so the
    # prose is dropped and only the code is looked at.
    code = builder.split('"""')[-1]
    assert '_append_js(ctx, fresh)' in code
    assert 'ctx.run_js' not in code
    assert "build = spec.pop('build', None)" in code, 'build once, not per rebuild'


def test_settings_keeps_a_live_handle(source):
    """Settings is handed the refs dicts by value at startup. If a tab is built
    later, its refs have to appear in those same dicts, or applying a new
    Archive folder silently does nothing."""
    assert 'refs6, refs_off = {}, {}' in source
    assert 'refs.update(made or {})' in source
    handed = re.search(r'tab_settings\.create_tab\((.*?)\)', source, re.S).group(1)
    assert 'archive_refs=refs6' in handed and 'office_refs=refs_off' in handed
