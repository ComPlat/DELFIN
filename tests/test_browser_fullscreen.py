"""The file browser can take the whole window, and give it back.

A document is read at the width of the page it was written for, not at the
width left over next to a file list. The tab lifts out of the page rather
than moving its children into an overlay: everything inside already sizes
itself from the tab's height, so handing the tab the viewport is the whole
change, and no widget is reparented -- reparenting an ipywidget re-renders
it, which is how a PDF page or a 3D canvas loses its state.

Python owns the state so the button always reflects it. Escape reaches the
same button as a click, so the two can never disagree.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


@pytest.fixture
def page(tmp_path):
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    scripts: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / 'office',
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    ctx.run_js = scripts.append
    widget, refs = browser.create_tab(ctx)
    return widget, refs, scripts, ctx


def _find(root, name):
    if name in (getattr(root, '_dom_classes', ()) or ()):
        return root
    for child in getattr(root, 'children', ()) or ():
        found = _find(child, name)
        if found is not None:
            return found
    return None


def _css(root):
    out = []

    def walk(node):
        value = getattr(node, 'value', None)
        if isinstance(value, str) and '<style>' in value:
            out.append(value)
        for child in getattr(node, 'children', ()) or ():
            walk(child)

    walk(root)
    return '\n'.join(out)


def _row_containing(root, name):
    """The box that directly holds the widget carrying this class."""
    for child in getattr(root, 'children', ()) or ():
        if name in (getattr(child, '_dom_classes', ()) or ()):
            return root
        found = _row_containing(child, name)
        if found is not None:
            return found
    return None


def test_the_button_sits_in_the_row_that_is_always_shown(page):
    """Text, spreadsheet, Word, PDF and 3D: several of those hide the
    content toolbar, so the button has to live in the header row that
    stays up whatever is open."""
    widget, _refs, _scripts, _ctx = page
    assert _find(widget, 'calc-fullscreen-btn') is not None

    header = _row_containing(widget, 'calc-fullscreen-btn')
    assert header is not None
    # The header row is the one with the download and 3D controls.
    siblings = {
        cls
        for child in header.children
        for cls in (getattr(child, '_dom_classes', ()) or ())
    }
    assert 'calc-view-toggle' in siblings, (
        'the button is no longer in the header row; it would disappear '
        'whenever the content toolbar is hidden')


def test_turning_it_on_lifts_the_tab_out_of_the_page(page):
    widget, _refs, scripts, _ctx = page
    button = _find(widget, 'calc-fullscreen-btn')
    scripts.clear()

    button.value = True
    sent = '\n'.join(scripts)
    assert "classList.toggle('calc-zen', true)" in sent

    scripts.clear()
    button.value = False
    assert "classList.toggle('calc-zen', false)" in '\n'.join(scripts)


def test_the_button_says_which_way_it_goes(page):
    widget, _refs, _scripts, _ctx = page
    button = _find(widget, 'calc-fullscreen-btn')
    before = button.description
    button.value = True
    assert button.description != before
    assert 'Esc' in button.tooltip
    button.value = False
    assert button.description == before


def test_the_tab_gets_the_viewport_not_just_a_bigger_box(page):
    widget, _refs, _scripts, _ctx = page
    css = _css(widget)
    assert '.calc-tab.calc-zen' in css
    rule = css[css.index('.calc-tab.calc-zen'):][:400]
    assert 'position:fixed !important' in rule
    assert 'height:100vh !important' in rule
    # The normal height rule is calc(100vh - 145px) and would otherwise win.
    assert 'max-height:100vh !important' in rule


def test_the_viewer_measures_again_after_the_tab_resizes(page):
    """The 3D viewer sets its size in pixels from the space it measured,
    so a tab that changed size leaves it at the old one."""
    widget, _refs, scripts, _ctx = page
    _find(widget, 'calc-fullscreen-btn').value = True
    sent = '\n'.join(scripts)
    assert 'calcResizeMolViewer_' in sent
    assert "dispatchEvent(new Event('resize'))" in sent


def test_escape_presses_the_button_rather_than_undoing_the_class(page):
    """Clearing the class directly would leave the button showing that the
    tab is still in fullscreen."""
    _widget, _refs, _scripts, ctx = page
    startup = '\n'.join(ctx.init_js_parts)
    assert "e.key !== 'Escape'" in startup
    escape_block = startup[startup.index("e.key !== 'Escape'"):][:500]
    assert 'calc-fullscreen-btn' in escape_block
    assert 'btn.click()' in escape_block
    assert 'calc-zen' in escape_block
