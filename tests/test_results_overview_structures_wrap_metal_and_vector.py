"""Results-Overview structures: wrapping, metals kept, vector PDFs.

Three behaviours the overview tab / exports must show:
- a long SMILES wraps inside its table cell instead of running off;
- the 2D drawing keeps the coordination metal (fresh draw beats the
  stale, metal-stripped SMILES.png of older runs);
- the LaTeX 3D images are vector PDFs, not rasters.
"""

import pathlib
import sys

from delfin.analysis import results_overview as ro
from delfin.dashboard import tab_results_overview as tab


def test_a_long_smiles_wraps_in_the_table_cell():
    col = [{"id": "smiles", "label": "SMILES"}]
    rows = [{"label": "c", "path": None, "values": {"smiles": "C" * 500}}]
    html = tab.rows_to_markdownish(col, rows)
    assert "overflow-wrap:anywhere" in html
    assert "word-break:break-all" in html
    # the SMILES VALUE cell must wrap; the label cell keeps nowrap on
    # purpose -- so assert on the value-cell style, not the whole body
    body = html.split("<tbody>")[1]
    value_cell = body.split("CCCC", 1)[0]
    assert "white-space:normal;overflow-wrap:anywhere" in value_cell
    assert "C" * 500 in html  # nothing cut


def test_the_2d_draw_keeps_the_coordination_metal():
    smi = "[Fe](C)(C)(n1ccccc1)"
    png = ro.smiles_png_bytes(smi)
    assert png and png[:4] == b"\x89PNG"
    svg = ro.smiles_svg_bytes(smi)
    assert svg and b"atom-0" in svg  # the metal atom glyph is drawn


def test_metal_parse_fallback_never_loses_the_metal():
    mol = ro._parse_smiles_keep_metals("[Fe](C)(C)(n1ccccc1)")
    assert mol is not None
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    assert "Fe" in syms
    # organic regression: unchanged parsing
    m2 = ro._parse_smiles_keep_metals("CCOc1ccccc1")
    assert m2 is not None and m2.GetNumAtoms() == 9


def test_the_latex_3d_images_are_vector_pdfs(tmp_path):
    import pytest

    calc = _first_metal_calc()
    if calc is None:
        pytest.skip("no metal calc under calc_dir on this host")
    for fn in (ro.initial_xyz_png_bytes, ro.optimized_3d_png_bytes):
        pdf = fn(calc, fmt="pdf")
        assert pdf[:5] == b"%PDF-", f"{fn.__name__} ignored fmt=pdf"
        assert b"/Subtype /Image" not in pdf  # no raster XObject


def _first_metal_calc():
    import os
    calc_dir = os.environ.get("DELFIN_TEST_CALC_DIR",
                              "/pfs/data6/home/ka/ka_ibcs/ka_xn0397/calc")
    if not os.path.isdir(calc_dir):
        return None
    for name in sorted(os.listdir(calc_dir))[:60]:
        path = os.path.join(calc_dir, name)
        if os.path.isdir(path) and os.path.exists(
                os.path.join(path, "initial.xyz")):
            return path
    return None
