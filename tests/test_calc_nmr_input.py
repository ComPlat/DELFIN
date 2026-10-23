import importlib.util
import sys
import types
from pathlib import Path


_ROOT = Path(__file__).resolve().parents[1] / "delfin"
_MODULE_PATH = _ROOT / "dashboard" / "tab_calculations_browser.py"

if "delfin" not in sys.modules:
    _PKG = types.ModuleType("delfin")
    _PKG.__path__ = [str(_ROOT)]
    _PKG.__version__ = "test"
    sys.modules["delfin"] = _PKG

if "delfin.dashboard" not in sys.modules:
    _DASHBOARD_PKG = types.ModuleType("delfin.dashboard")
    _DASHBOARD_PKG.__path__ = [str(_ROOT / "dashboard")]
    sys.modules["delfin.dashboard"] = _DASHBOARD_PKG

# No stub for delfin.dashboard.helpers: the browser module's "from .helpers
# import ..." must resolve against the REAL module. A leaf-module stub has no
# __file__, so any name it forgets to pre-attach fails collection here with
# "cannot import name ... (unknown location)" (that is exactly how
# js_string_literal broke this file). The real module is cheap to import —
# the browser module already pulls in ipywidgets and .constants itself.
_SPEC = importlib.util.spec_from_file_location("delfin.dashboard.tab_calculations_browser", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load calculations browser module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)


def test_build_calc_nmr_input_uses_selected_cpcm_solvent():
    inp = _MODULE.build_calc_nmr_input(
        ["C 0.0 0.0 0.0", "H 0.0 0.0 1.0"],
        pal=8,
        maxcore=2500,
        solvent="dmso",
    )

    assert "!TPSS PCSSEG-1 AUTOAUX NMR CPCM(dmso)" in inp
    assert "nprocs 8" in inp
    assert "%maxcore 2500" in inp
    assert "  C 0.0 0.0 0.0" in inp
    assert "  H 0.0 0.0 1.0" in inp
    assert "NUCLEI = ALL H {SHIFT, SSALL}" in inp


def test_build_calc_nmr_input_defaults_to_chloroform_for_blank_solvent():
    inp = _MODULE.build_calc_nmr_input(
        ["H 0.0 0.0 0.0"],
        pal=1,
        maxcore=3000,
        solvent="",
    )

    assert "CPCM(chloroform)" in inp
