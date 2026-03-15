import importlib.util
from pathlib import Path


_MODULE_PATH = Path(__file__).resolve().parents[1] / "delfin" / "dashboard" / "energy_labels.py"
_SPEC = importlib.util.spec_from_file_location("delfin_tab_calculations_browser", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load energy_labels module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


def test_extract_comment_energy_from_goat_label():
    energy = _MODULE.extract_comment_energy("run_01 -105.779061555980 isomer")
    assert energy == -105.77906155598


def test_extract_comment_energy_from_orca_comment_path_label():
    comment = (
        "Coordinates from ORCA-job "
        "/home/qmchem_max/calc/Fabian_Co4/ox_step_1_OCCUPIER/"
        ".orca_iso_input2_309535_135210566022848_dfe6fe08/input2 "
        "E -2675.849910594636"
    )
    energy = _MODULE.extract_comment_energy(comment)
    assert energy == -2675.849910594636


def test_extract_comment_energy_from_plain_xyz_comment_line():
    energy = _MODULE.extract_comment_energy("-83.6028516932")
    assert energy == -83.6028516932


def test_format_xyz_comment_label_adds_absolute_hartree_and_kcal():
    frames = [
        ("run_01 -105.779061555980 isomer", "H 0 0 0", 1),
        ("run_02 -105.777467025622 isomer", "H 0 0 0", 1),
    ]

    rendered = _MODULE.format_xyz_comment_label(frames[1][0], frames=frames)

    assert "run_02 -105.777467025622 isomer" not in rendered
    assert "E = -105.777467025622 Hartree;" in rendered
    assert "Hartree" in rendered
    assert "kcal/mol" in rendered
    assert "ΔE" not in rendered


def test_format_xyz_comment_label_single_frame_has_only_absolute_units():
    comment = "run_01 -105.779061555980 isomer"

    rendered = _MODULE.format_xyz_comment_label(comment, frames=[(comment, "H 0 0 0", 1)])

    assert "run_01 -105.779061555980 isomer" not in rendered
    assert "E = -105.779061555980 Hartree;" in rendered
    assert "Hartree" in rendered
    assert "kcal/mol" in rendered
    assert "ΔE" not in rendered
