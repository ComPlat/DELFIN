"""Fullscreen outside the Submit tab took the picture and left its controls.

In the Submit tab the whole toolbar goes fullscreen with the viewer.  In the
Calculations tab and the Archive it did not: the RMSD pair and the Fukui table
stayed behind on the page while the viewer moved into the overlay, so going
fullscreen to look closely at what those controls had produced took the
controls away.  Every part of a molecule panel is a member now.

The proportions below were measured in a real browser at 1440x900 against the
module's own CSS and its own bootstrap script -- see the class rules, not the
numbers, for what is being asserted here.
"""

import tempfile
from pathlib import Path

import pytest

from delfin.dashboard import molecule_viewer as viewer


CSS = viewer.STRUCTURE_VIEWER_FULLSCREEN_CSS


def _source(module):
    return open(module.__file__, encoding='utf-8').read()


def _find_module(root, marker, trail=()):
    """The box a tab put *marker* on, and the boxes it hangs under.

    The trail matters: the script reads a scope class by walking up from the
    module, the way ``element.closest`` does, and a tab is free to put that
    class on the module or on anything above it. The ORCA Builder puts it on
    the whole tab.
    """
    here = trail + (root,)
    if marker in (getattr(root, '_dom_classes', ()) or ()):
        return here
    for child in getattr(root, 'children', ()) or ():
        found = _find_module(child, marker, here)
        if found is not None:
            return found
    return None


@pytest.fixture(scope='module')
def built_tabs():
    """Every tab that has a molecule module, built for real."""
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import (
        tab_calculations_browser, tab_orca_builder, tab_remote_archive,
        tab_submit,
    )
    from delfin.dashboard.context import DashboardContext

    tmp = Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (tmp / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp / 'calc', archive_dir=tmp / 'archive',
                           office_dir=tmp / 'office')
    ctx.run_js = lambda _script: None
    makers = {
        'submit': tab_submit,
        'orca': tab_orca_builder,
        'calc': tab_calculations_browser,
        'remote': tab_remote_archive,
    }
    out = {}
    for kind, module in makers.items():
        widget = module.create_tab(ctx)[0]
        found = _find_module(widget, f'{kind}-structure-fs-module')
        assert found is not None, f'{kind} builds no molecule module'
        out[kind] = found
    return out


def _classes_up_the_trail(trail):
    """Every class the script would see walking up from the module."""
    seen = []
    for widget in trail:
        seen.extend(getattr(widget, '_dom_classes', ()) or ())
    return seen


def test_the_calculations_panel_travels_whole():
    from delfin.dashboard import tab_calculations_browser as calc

    source = _source(calc)
    for widget in ('calc_mol_header', 'calc_mol_view_row',
                   'calc_rmsd_controls', 'calc_fukui_panel_container'):
        assert f"{widget}.add_class('delfin-structure-fs-member')" in source, \
            f'{widget} would stay behind on the page'


def test_the_archive_panel_travels_whole():
    from delfin.dashboard import tab_remote_archive as remote

    source = _source(remote)
    for widget in ('viewer_header', 'viewer_row',
                   'remote_fukui_panel_container'):
        assert f'{widget}.add_class("delfin-structure-fs-member")' in source, \
            f'{widget} would stay behind on the page'


def test_a_panel_is_bounded_and_scrolls():
    """A filled Fukui table is taller than the screen.  Carried unbounded it
    would leave nothing of the picture, which is what fullscreen is for."""
    rule = CSS.split('.delfin-structure-fs-overlay > .delfin-structure-fs-panel {')[1]
    rule = rule.split('}')[0]
    assert 'max-height: 30vh' in rule
    assert 'overflow-y: auto' in rule
    # It may shrink -- with two panels open there is not room for both at full
    # height -- but it never grows past its share.
    assert 'flex: 0 1 auto' in rule


def test_the_picture_keeps_a_floor_under_it():
    """Measured at 1440x900 with the real CSS: the RMSD pair alone leaves the
    picture 77% of the height, a filled Fukui table alone 64%, and both
    together 46%.  Without the floor both together left it 44% and squeezed
    the RMSD controls to 55 px, which is not a control any more."""
    rule = CSS.split(
        '.delfin-structure-fs-overlay > .delfin-structure-fs-view-row {')[1]
    rule = rule.split('}')[0]
    assert 'min-height: 45vh' in rule

    # And the narrow layout gives the panel less, because a column of picture
    # and panel on a small screen cannot afford both.
    narrow = CSS.split('@media (max-width: 800px) {')[1]
    assert 'max-height: 24vh' in narrow


def test_the_orca_preview_is_untouched():
    """It has no panels and no view row -- it puts the viewer straight into the
    overlay -- so the floor and the panel rule do not reach it."""
    from delfin.dashboard import tab_orca_builder as orca

    source = _source(orca)
    assert "orca_mol_output.add_class('delfin-structure-fs-viewer')" in source
    assert 'delfin-structure-fs-view-row' not in source
    assert 'delfin-structure-fs-panel' not in source


# ---------------------------------------------------------------------------
# one fullscreen, and how a tab joins it
# ---------------------------------------------------------------------------

#: Each tab, and the prefix of the scope class its module actually carries.
#:
#: Not the kind word.  Three of these four make their own scope id and name it
#: after themselves, which makes "the prefix is the kind" look like the rule --
#: and it was written down here as one.  The ORCA Builder borrows the editor's
#: scope instead (orca_mol_module.add_class(orca_editor_scope)), so its class
#: begins with submit-scope-, and this table said orca-scope-.
#:
#: What that cost: classWithPrefix walked the ancestors, found nothing and
#: returned null, so in fullscreen the module had no scope at all and every
#: per-scope thing the overlay does on the way in was skipped -- no atom in
#: that tab could be picked in the enlarged picture, and no band drawn.  The
#: check agreed with the mistake, confidently, for as long as it existed.
TABS = {
    'submit': ('tab_submit', 'submit-scope-'),
    'orca': ('tab_orca_builder', 'submit-scope-'),
    'calc': ('tab_calculations_browser', 'calc-scope-'),
    'remote': ('tab_remote_archive', 'remote-archive-scope-'),
}


def _tab_source(name):
    import importlib

    return _source(importlib.import_module('delfin.dashboard.' + name))


def test_there_is_one_fullscreen_and_the_editor_does_not_carry_a_second():
    """The Submit tab had one of its own and the other three shared another.

    Two overlays, two stylesheets, two restore paths -- and a fix to either
    was a fix to one tab or to three, never to all four. The rescue for a
    member whose parent was replaced was in the Submit one only, and driven in
    chromium the other three lost their whole molecule panel to it.
    """
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor_js = submit_manip_bootstrap_js()
    for gone in ('enterFullscreen', 'exitFullscreen', '__delfinSubmitFullscreenBound',
                 '_submitFsByScope', 'submit-fs-overlay'):
        assert gone not in editor_js, f'the editor still carries {gone}'

    assert viewer.STRUCTURE_VIEWER_FULLSCREEN_BOOTSTRAP_JS.count(
        'function enterFullscreen') == 1

    # And the Submit tab is laid out by the shared sheet rather than by rules
    # of its own, which is what let the two drift apart.
    submit = _tab_source('tab_submit')
    assert 'structure_viewer_fullscreen_css()' in submit
    assert '.submit-fs-overlay' not in submit


def test_every_tab_declares_itself_and_the_script_knows_no_tab():
    """A builder joins by saying what it is, not by being added to a switch.

    The shared script used to name orca, calc and remote in three functions --
    which kind a module is, what its scope class is called, where its viewer
    lives. A fourth tab meant editing all three, and a TURBOMOLE builder would
    have meant editing them again.
    """
    js = viewer.STRUCTURE_VIEWER_FULLSCREEN_BOOTSTRAP_JS
    for tab in TABS:
        assert f"'{tab}'" not in js and f'"{tab}"' not in js, (
            f'the shared script should not know what {tab} is')
    assert '__delfinFsKinds' in js

    for kind, (module, prefix) in TABS.items():
        source = _tab_source(module)
        # The module box carries the word, and the tab declares that word
        # together with the prefix of the class its module really has -- which
        # is not always a prefix named after the word.  See TABS.
        assert f"add_class('{kind}-structure-fs-module')" in source \
            or f'add_class("{kind}-structure-fs-module")' in source, kind
        declaration = source.split('structure_viewer_fullscreen_kind_js(')[1]
        declaration = declaration.split(')')[0]
        assert f"'{kind}'" in declaration, kind
        assert f"'{prefix}" in declaration, (
            f'{kind} declares a prefix its module carries no class for')


def test_the_module_a_tab_builds_carries_the_scope_it_declared(built_tabs):
    """Read off the widgets rather than the source: the declaration is only
    worth anything if the class the box actually carries starts with it.

    The Submit tab's scope is not named in the Submit tab at all -- the editor
    builds it, and the tab puts it on the column. A test that read the tab's
    text for it would have found nothing and said the tab was wrong.
    """
    for kind, trail in built_tabs.items():
        prefix = TABS[kind][1]
        module = trail[-1]
        own = list(module._dom_classes)
        assert f'{kind}-structure-fs-module' in own, kind
        assert 'delfin-structure-fs-module' in own, kind
        reachable = _classes_up_the_trail(trail)
        assert any(c.startswith(prefix) for c in reachable), (
            f'{kind} declares {prefix}, and nothing from its module up to the '
            f'tab carries it: {reachable}')


def test_the_declaration_is_a_line_a_tab_can_write():
    kind_js = viewer.structure_viewer_fullscreen_kind_js(
        'turbomole', 'turbomole-scope-', ['_turbomoleViewerByScope'])
    assert '__delfinFsKinds' in kind_js
    assert '"turbomole"' in kind_js
    assert '"turbomole-scope-"' in kind_js
    assert '["_turbomoleViewerByScope"]' in kind_js
