"""Tests for the Fukui viewer helpers in delfin.dashboard.molecule_viewer."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.dashboard import molecule_viewer as mv


_XYZ = """4
formaldehyde
C 0.00 0.00 0.00
O 1.20 0.00 0.00
H -0.50 0.90 0.00
H -0.50 -0.90 0.00
"""

_CUBE = """\
Test cube
density
    4    0.000000    0.000000    0.000000
    2    1.000000    0.000000    0.000000
    2    0.000000    1.000000    0.000000
    2    0.000000    0.000000    1.000000
    6    6.000000    0.000000    0.000000    0.000000
    8    8.000000    1.200000    0.000000    0.000000
    1    1.000000   -0.500000    0.900000    0.000000
    1    1.000000   -0.500000   -0.900000    0.000000
 1.00000E-01  2.00000E-01  3.00000E-01  4.00000E-01  5.00000E-01  6.00000E-01
 7.00000E-01  8.00000E-01
"""


def test_atom_positions_from_xyz_with_header():
    pos = mv._fukui_atom_positions_from_xyz(_XYZ)
    assert len(pos) == 4
    assert pos[0][0] == 'C'
    assert pos[0][1] == pytest.approx(0.0)
    assert pos[1][0] == 'O'
    assert pos[1][1] == pytest.approx(1.2)


def test_atom_labels_js_emits_one_per_atom():
    js = mv.fukui_atom_labels_js(_XYZ, [0.123, -0.456, 0.05, 0.05])
    assert js.count("viewer.addLabel(") == 4
    assert "+0.123" in js
    assert "-0.456" in js


def test_atom_labels_js_decimals_respected():
    js = mv.fukui_atom_labels_js(_XYZ, [0.123456, -0.456789, 0.0, 0.0], decimals=5)
    assert "+0.12346" in js  # rounded to 5 decimals
    assert "-0.45679" in js


def test_atom_labels_js_passthrough_strings():
    js = mv.fukui_atom_labels_js(_XYZ, ['1', '2', '3', '4'])
    assert "viewer.addLabel('1'" in js
    assert "viewer.addLabel('4'" in js


def test_atom_labels_js_empty_values_returns_empty_string():
    assert mv.fukui_atom_labels_js(_XYZ, None) == ''
    assert mv.fukui_atom_labels_js(_XYZ, []) == ''


def test_cube_isosurface_js_unsigned_single_lobe():
    js = mv.fukui_cube_isosurface_js(_CUBE, isoval=0.02, signed=False)
    assert js.count("addVolumetricData") == 1
    assert "#0026ff" in js


def test_cube_isosurface_js_signed_paints_both_lobes():
    js = mv.fukui_cube_isosurface_js(_CUBE, isoval=0.005, signed=True)
    assert js.count("addVolumetricData") == 2
    assert "#0026ff" in js
    assert "#b00010" in js
    # Negative lobe at -isoval
    assert "isoval:-0.005" in js
    assert "isoval:0.005" in js


def test_cube_isosurface_js_empty_input_returns_empty():
    assert mv.fukui_cube_isosurface_js('', isoval=0.02) == ''


def test_build_fukui_viewer_html_includes_xyz_and_styles():
    html = mv.build_fukui_viewer_html(_XYZ)
    assert "3Dmol-min.js" in html
    assert "addModel" in html
    # XYZ ist als JSON eingebettet, darf nicht roh stehen.
    assert "formaldehyde" in html


def test_build_fukui_viewer_html_with_labels_and_cube():
    html = mv.build_fukui_viewer_html(
        _XYZ, labels=[0.1, 0.2, 0.3, 0.4],
        cube_text=_CUBE, isoval=0.02, cube_signed=True,
    )
    assert "viewer.addLabel" in html
    assert "addVolumetricData" in html


def test_build_fukui_table_html_renders_all_rows():
    result = {
        'atoms': ['C', 'O', 'H', 'H'],
        'q_neutral': [0.1, -0.3, 0.1, 0.1],
        'q_anion':   [-0.2, -0.4, 0.05, 0.05],
        'q_cation':  [0.2, 0.5, 0.15, 0.15],
        'f_plus':    [0.3, 0.1, 0.05, 0.05],
        'f_minus':   [0.1, 0.8, 0.05, 0.05],
        'f_zero':    [0.2, 0.45, 0.05, 0.05],
        'scheme': 'mulliken',
        'orca_settings': {'functional': 'B3LYP', 'basis': 'def2-SVP'},
        'geometry_origin': 'opt',
    }
    html = mv._build_fukui_table_html(result)
    assert '<table' in html
    assert 'B3LYP' in html
    assert 'mulliken' in html
    # Exactly 4 data rows
    assert html.count('<tr>') == 4


def _write_fukui_jobdir(tmp_path: Path) -> Path:
    workdir = tmp_path / 'fukui_run'
    workdir.mkdir()
    (workdir / 'fukui_geom.xyz').write_text(_XYZ)
    (workdir / 'density_neutral.cube').write_text(_CUBE)
    (workdir / 'fukui_plus.cube').write_text(_CUBE)
    (workdir / 'fukui_result.json').write_text(json.dumps({
        'atoms': ['C', 'O', 'H', 'H'],
        'q_neutral': [0.1, -0.3, 0.1, 0.1],
        'q_anion':   [-0.2, -0.4, 0.05, 0.05],
        'q_cation':  [0.2, 0.5, 0.15, 0.15],
        'f_plus':    [0.3, 0.1, 0.05, 0.05],
        'f_minus':   [0.1, 0.8, 0.05, 0.05],
        'f_zero':    [0.2, 0.45, 0.05, 0.05],
        'scheme': 'mulliken',
        'orca_settings': {'functional': 'B3LYP'},
        'geometry_origin': 'opt',
    }))
    return workdir


def test_render_fukui_panel_returns_widget(tmp_path):
    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets
    workdir = _write_fukui_jobdir(tmp_path)
    output = widgets.Output()
    panel = mv.render_fukui_panel(workdir, output)
    assert isinstance(panel, widgets.VBox)
    # Outer wrap contains one HBox (viewer + right-hand sidebar).
    assert len(panel.children) == 1
    viewer_row = panel.children[0]
    assert isinstance(viewer_row, widgets.HBox)
    # The viewer must be an Output widget so 3Dmol's <script> bootstrap
    # actually executes (widgets.HTML would sanitise script tags away).
    viewer = viewer_row.children[0]
    assert isinstance(viewer, widgets.Output)
    sidebar = viewer_row.children[1]
    assert isinstance(sidebar, widgets.VBox)
    # Sidebar = controls VBox (4 children) + result-table HTML.
    assert len(sidebar.children) == 2
    controls = sidebar.children[0]
    assert isinstance(controls, widgets.VBox)
    assert len(controls.children) == 4


def test_render_fukui_panel_missing_json(tmp_path):
    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets
    output = widgets.Output()
    panel = mv.render_fukui_panel(tmp_path, output)
    # Returns an error HTML widget instead of raising
    assert isinstance(panel, widgets.HTML)
    assert 'Failed' in panel.value
