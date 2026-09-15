"""What a viewer with nothing to show says, and where it says it.

The line was stream text, and stream text lands where console output lands:
the top left corner.  In a 560 px tall frame that one small line in the far
top corner read as a picture that had failed to arrive, not as the
invitation it was.
It is what the box is waiting for, so it goes where a picture keeps its
caption: the bottom left corner of the frame.
"""

from __future__ import annotations

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
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda js: None
    _widget, refs = tab_orca_builder.create_tab(ctx)
    return refs


def _the_notice(refs):
    """The HTML of whatever the viewer is saying, '' when it is a structure."""
    return ''.join(
        (one.get('data', {}) or {}).get('text/html', '') or ''
        for one in refs['orca_mol_output'].outputs
    )


def test_an_empty_viewer_says_it_at_the_bottom_left():
    refs = _a_builder()
    # The tab builds with the invitation already in the box.
    notice = _the_notice(refs)

    assert 'Paste XYZ coordinates to see 3D preview.' in notice, notice[:200]
    # Carried by the class the tab's sheet pins to the bottom left corner --
    # stream text put the same words at the top left instead.
    assert 'orca-mol-note' in notice, notice[:200]
    # A note, not console output: the top left is where a stream line goes.
    streams = [one for one in refs['orca_mol_output'].outputs
               if one.get('output_type') == 'stream']
    assert not streams, streams


def test_every_empty_viewer_message_is_the_same_kind_of_note():
    refs = _a_builder()
    # A SMILES is not coordinates; the box says so instead of drawing the
    # first letter of one as a one-atom model.
    refs['orca_coords'].value = 'c1ccccc1'
    notice = _the_notice(refs)

    assert 'SMILES detected.' in notice, notice[:200]
    assert 'orca-mol-note' in notice, notice[:200]

    # And the same for what is neither coordinates nor a SMILES.
    refs['orca_coords'].value = 'not a molecule at all'
    notice = _the_notice(refs)

    assert 'orca-mol-note' in notice, notice[:200]


def test_the_sheet_pins_the_note_to_the_bottom_left():
    from delfin.dashboard import tab_orca_builder
    import inspect

    source = inspect.getsource(tab_orca_builder)
    sheet = source.split('.orca-mol-note {', 1)
    assert len(sheet) == 2, 'no rule for .orca-mol-note in the tab sheet'
    rule = sheet[1].split('}', 1)[0]
    assert 'bottom' in rule, rule
    assert 'left' in rule, rule
    assert 'position: absolute' in rule, rule


def test_a_structure_still_draws_where_the_note_was():
    refs = _a_builder()
    refs['orca_coords'].value = _WATER
    drawn = _the_notice(refs)

    assert 'orca-mol-note' not in drawn, drawn[:200]
    assert '3dmolviewer' in drawn, drawn[:200]
