import importlib.util
import json
from pathlib import Path


_MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "molecule_viewer.py"
)
_SPEC = importlib.util.spec_from_file_location("delfin_molecule_viewer", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load molecule_viewer module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


_SETTINGS_UI_SOURCE = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "tab_settings.py"
).read_text(encoding="utf-8")


def test_settings_ui_exposes_and_persists_viewer_controls():
    for widget_name in (
        "viewer_representation_dropdown",
        "viewer_atom_scale_slider",
        "viewer_bond_radius_slider",
        "viewer_multiple_bonds_toggle",
        "viewer_depth_fog_toggle",
        "viewer_ambient_occlusion_toggle",
    ):
        assert f"{widget_name} = widgets." in _SETTINGS_UI_SOURCE

    for settings_key in (
        "'representation': str(viewer_representation_dropdown.value)",
        "'atom_scale': float(viewer_atom_scale_slider.value)",
        "'bond_radius': float(viewer_bond_radius_slider.value)",
        "'multiple_bonds': bool(viewer_multiple_bonds_toggle.value)",
        "'depth_fog': bool(viewer_depth_fog_toggle.value)",
        "'ambient_occlusion': bool(viewer_ambient_occlusion_toggle.value)",
    ):
        assert settings_key in _SETTINGS_UI_SOURCE
    assert "XYZ-/ORCA-Koordinaten enthalten keine" in _SETTINGS_UI_SOURCE
    assert "nur MOL/SDF mit Bindungsordnung" in _SETTINGS_UI_SOURCE


def test_viewer_representation_styles_use_global_dimensions():
    ball_and_stick = _MODULE.build_molecule_view_style(
        "ball_and_stick",
        atom_scale=0.42,
        bond_radius=0.17,
        multiple_bonds=True,
    )
    assert set(ball_and_stick) == {"stick", "sphere"}
    assert ball_and_stick["sphere"]["scale"] == 0.42
    assert ball_and_stick["stick"]["radius"] == 0.17
    assert ball_and_stick["stick"]["singleBonds"] is False

    stick = _MODULE.build_molecule_view_style(
        "stick",
        atom_scale=0.9,
        bond_radius=0.23,
        multiple_bonds=False,
    )
    assert set(stick) == {"stick"}
    assert stick["stick"]["radius"] == 0.23
    assert stick["stick"]["singleBonds"] is True

    sphere = _MODULE.build_molecule_view_style("sphere", atom_scale=1.0)
    assert set(sphere) == {"sphere"}
    assert sphere["sphere"]["scale"] == 1.0

    line = _MODULE.build_molecule_view_style("line", bond_radius=0.22)
    assert set(line) == {"line"}
    assert line["line"]["linewidth"] == 2.0


def test_overlay_style_keeps_representation_and_replaces_element_colors():
    style = _MODULE.build_molecule_view_style(
        "ball_and_stick",
        atom_scale=0.36,
        bond_radius=0.15,
    )

    colored = json.loads(
        _MODULE.molecule_view_style_js(style, color="#1f5fff")
    )

    assert set(colored) == {"stick", "sphere"}
    assert colored["stick"]["color"] == "#1f5fff"
    assert colored["sphere"]["color"] == "#1f5fff"
    assert "colorscheme" not in colored["stick"]
    assert colored["stick"]["radius"] == 0.15
    assert colored["sphere"]["scale"] == 0.36


def test_quality_profiles_control_renderer_independently():
    low = _MODULE.VIEWER_QUALITY_PROFILES["low"]["viewer_config"]
    medium = _MODULE.VIEWER_QUALITY_PROFILES["medium"]["viewer_config"]
    high = _MODULE.VIEWER_QUALITY_PROFILES["high"]["viewer_config"]

    assert low["antialias"] is False
    assert medium["antialias"] is True
    assert "style" not in medium
    assert high["antialias"] is True
    assert "style" not in high

    assert _MODULE.build_viewer_config("high", depth_fog=True)["disableFog"] is False
    assert _MODULE.build_viewer_config("high", depth_fog=False)["disableFog"] is True
    ao = _MODULE.build_viewer_config(
        "high",
        depth_fog=True,
        ambient_occlusion=True,
    )
    assert ao["style"] == "ambientOcclusion"


def test_py3dmol_style_application_injects_renderer_config(monkeypatch):
    profile = {
        "style": {"sphere": {"colorscheme": "Jmol", "scale": 0.8}},
        "viewer_config": {
            "backgroundColor": "white",
            "antialias": False,
        },
        "viewer_config_js": '{"backgroundColor":"white","antialias":false}',
    }
    monkeypatch.setattr(_MODULE, "get_viewer_profile", lambda: profile)

    class FakeView:
        def __init__(self):
            self.startjs = (
                'viewer_UNIQUEID = $3Dmol.createViewer(el,'
                '{backgroundColor:"white"});'
            )
            self.style = None
            self.background = None

        def setStyle(self, selection, style):
            self.style = (selection, style)

        def setBackgroundColor(self, color):
            self.background = color

        def zoomTo(self):
            pass

        def center(self):
            pass

        def zoom(self, value):
            pass

        def render(self):
            pass

    view = FakeView()
    _MODULE.apply_molecule_view_style(view)

    assert '{"backgroundColor":"white","antialias":false}' in view.startjs
    assert view.style == ({}, profile["style"])
    assert view.background == "white"
