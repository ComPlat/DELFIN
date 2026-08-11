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

from delfin.dashboard import molecule_viewer as viewer


CSS = viewer.STRUCTURE_VIEWER_FULLSCREEN_CSS


def _source(module):
    return open(module.__file__, encoding='utf-8').read()


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
