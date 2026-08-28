"""The ORCA Builder shows a structure with the editor, not with a viewer of
its own.

It used to build its own: ``window._orcaBuildViewer``, its own template, its
own labels, introduced to the editor afterwards.  That made every capability
of the editor a thing that had to work twice, and the second one kept quietly
not working -- in fullscreen no atom in this tab could be picked at all,
because the tab told the overlay a scope prefix no class of its own began
with, and a module without a scope is also what a tab with no editor looks
like.  Nothing anywhere said so.

The editor could not be used for it, because its one picture was the one
thing it wrote straight into a widget of its own instead of handing to the
host.  That is now the same channel as everything else it shows, and this
tab hands it the widget it places.

What is left of the tab's own renderer is the numbering check, whose pictures
are two structures at once and are explicitly not editable.
"""

from __future__ import annotations

import inspect
import pathlib
import tempfile

import pytest

from delfin.dashboard.context import DashboardContext

_WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n"


def _a_builder():
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import tab_orca_builder

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    sent: list[str] = []
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = sent.append
    _widget, refs = tab_orca_builder.create_tab(ctx)
    return refs, sent


def _drawn(refs):
    return ''.join(
        (one.get('data', {}) or {}).get('text/html', '') or ''
        for one in refs['orca_mol_output'].outputs)


def test_the_builder_shows_the_editors_picture():
    refs, _sent = _a_builder()
    refs['orca_coords'].value = _WATER
    html = _drawn(refs)

    assert html, 'nothing was drawn'
    # py3Dmol names the editor's viewer; the tab's own template never did.
    assert '3dmolviewer' in html, html[:200]
    # And it is registered where the editor looks for it: by scope, not by a
    # name of one tab's own.
    assert '_submitMolViewerByScope' in html
    assert refs['submit_scope_id'] in html

    # Which is what makes the toolbar live over it.
    assert not refs['submit_select_btn'].disabled
    assert not refs['submit_manip_btn'].disabled


def test_the_tab_has_no_second_renderer_for_a_structure():
    from delfin.dashboard import tab_orca_builder

    source = inspect.getsource(tab_orca_builder)
    # The in-place swap the tab used to do inside its own viewer is gone with
    # the viewer: the editor keeps the camera across a re-render itself, which
    # is what that path existed for.
    assert '_update_molecule_js' not in source
    assert '_show_molecule_in_place' not in source
    # What remains of the tab's own renderer is the comparison pictures.
    for kept in ('_overlay_viewer_html', '_numbering_check_view_html'):
        assert kept in source, kept


def test_the_editor_hands_its_picture_to_the_host():
    """The one picture the editor wrote into a widget of its own.

    Every other thing it shows goes through show_output, which is how a tab
    that places the viewer itself receives it.  This one did not, and that is
    the whole reason a second viewer had to exist.
    """
    from delfin.dashboard import structure_editor

    source = inspect.getsource(structure_editor)
    place = source.split('def _replace_mol_output_view')[1].split('\n    def ')[0]
    assert 'show_output(_build_mol_output_bundle(' in place
    assert 'mol_output.outputs = _build_mol_output_bundle' not in source
