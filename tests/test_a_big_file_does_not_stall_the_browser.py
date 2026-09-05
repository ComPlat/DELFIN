"""A large output file must not arrive as one node of laid-out text.

The viewer used to hand the whole file to the page as markup, under
``white-space: pre-wrap``. The browser then had to line-break and lay out
every line before it could paint: on a 3.3 MB, 73 000 line ORCA output that
is around 7.5 seconds during which the tab answers nothing, which is the
"page unresponsive" prompt people were seeing on an ``S1.out``.

The text now goes to ``delfin.dashboard.text_view``, which cuts it into
blocks of whole lines that declare their height, so the browser lays out
only the blocks on screen. These tests hold the two halves of that: the file
does not travel as markup, and it does arrive at the block viewer.
"""

from __future__ import annotations

import re

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard import text_view
from delfin.dashboard.context import DashboardContext

# Distinctive enough that finding it proves this text, not some other.
NEEDLE = 'TOTAL RUN TIME MARKER 4711'


def _open_file(tmp_path, name, text):
    """Build a browser on a folder and select one file in it."""
    for folder in ('calc', 'archive', 'office'):
        (tmp_path / folder).mkdir(exist_ok=True)
    (tmp_path / 'calc' / name).write_text(text, encoding='utf-8')

    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / 'calc',
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    ctx.run_js = lambda script: sent.append(script)
    widget, refs = browser.create_tab(ctx)
    refs['calc_list_directory']()

    file_list = refs['calc_file_list']
    target = [opt for opt in file_list.options if name in str(opt)]
    assert target, f'{name} did not show up in {file_list.options}'
    value = target[0][1] if isinstance(target[0], tuple) else target[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value
    return widget, refs, sent


def _widget_values(root):
    """Every string a widget is holding, anywhere in the tab."""
    seen = set()
    out = []

    def walk(node):
        if id(node) in seen:
            return
        seen.add(id(node))
        value = getattr(node, 'value', None)
        if isinstance(value, str):
            out.append(value)
        for child in (getattr(node, 'children', None) or ()):
            walk(child)

    walk(root)
    return out


def _big_output(lines=60_000):
    """Something shaped like an ORCA .out: many short lines."""
    body = '\n'.join(
        f'  {i:6d}   -1234.567890123   0.000012345   SCF CYCLE {i % 97}'
        for i in range(lines)
    )
    return body + f'\n{NEEDLE}\n'


# Under the full-read limit, and over it, where the viewer loads windows.
SIZES = [('small.out', 2_000), ('big.out', 60_000)]


@pytest.mark.parametrize('name, lines', SIZES)
def test_the_file_is_not_sent_as_markup(tmp_path, name, lines):
    widget, _refs, _sent = _open_file(tmp_path, name, _big_output(lines))
    values = _widget_values(widget)
    assert not any(NEEDLE in v for v in values), (
        'the file is in a widget value, so it reaches the page as markup: '
        'ipywidgets carries a second copy of it and the browser lays out '
        'every line of it before it can paint'
    )
    assert any('dtv-text' in v for v in values), (
        'the content area is not the block viewer, so whatever fills it is '
        'laying out the whole file at once'
    )


@pytest.mark.parametrize('name, lines', SIZES)
def test_the_file_reaches_the_block_viewer(tmp_path, name, lines):
    _widget, _refs, sent = _open_file(tmp_path, name, _big_output(lines))
    fills = [s for s in sent if '__delfinTextView.setText(' in s]
    assert fills, 'nothing handed the text to the block viewer'
    assert any(NEEDLE in s for s in fills), (
        'the block viewer was called without the file in it'
    )


def test_the_runtime_travels_with_the_text(tmp_path):
    """A file can be opened before the page's startup scripts have run.

    If the runtime were only registered at startup, that file would silently
    not appear, so the definition ships with the first thing that needs it.
    """
    _widget, _refs, sent = _open_file(tmp_path, 'big.out', _big_output(60_000))
    fills = [s for s in sent if '__delfinTextView.setText(' in s]
    assert fills
    for script in fills:
        assert 'window.__delfinTextView = TV;' in script, (
            'setText is called without the runtime being defined in the same '
            'script, so it depends on startup order'
        )


def test_one_script_carries_the_fill_and_the_scroll(tmp_path):
    """run_js clears its output before each call.

    Two calls in a row race: the second wipes the first before the browser
    has run it. Rendering and then scrolling has to be one script.
    """
    _widget, refs, sent = _open_file(tmp_path, 'big.out', _big_output(60_000))
    sent.clear()
    refs['calc_top_btn'].click()
    fills = [s for s in sent if '__delfinTextView.setText(' in s]
    for script in fills:
        assert 'scrollTop' in script, (
            'the text went out without the scroll that follows it, so the '
            'next run_js call clears it before the browser has run it'
        )


def test_the_runtime_is_valid_javascript():
    """The Python parser is happy with a JS string that does not parse."""
    boot = text_view.text_view_bootstrap_js()
    # Placeholders must all be substituted before the script is emitted.
    assert not re.search(r'__[A-Z_]+__', boot), (
        'an unsubstituted placeholder is a syntax error in the browser'
    )
    assert boot.count('{') == boot.count('}')
    assert boot.count('(') == boot.count(')')
