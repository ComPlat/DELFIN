"""Call-site switch tests: final-geometry readers go through delfin.orca_read.

The chaining readers (IMAG's displaced-geometry rounds, classic recovery's
continue-from-own-geometry, ORCA's resume path) parse the ORCA-written
``<base>.xyz`` themselves today. This test pins the switch onto
``delfin.orca_read`` with its text fallback -- behaviour-identical by the
baseline recorder (scratch/baseline.json, 207 runs) -- and that
``orca.reader=text`` restores the old path exactly.

Behaviour contract, identical pre- and post-switch:
- imag._read_xyz returns [(element, x, y, z), ...] with the same
  coordinates the sibling ORCA-written xyz carries;
- orca_recovery._continue_from_own_geometry injects the same rows into
  the next input (via orca_input.with_coordinates formatting);
- orca._read_xyz_atoms returns the same raw rows or None.
"""
from pathlib import Path
from unittest import mock

from delfin import imag, orca, orca_read, orca_recovery
from delfin.common import orca_input

FIX = Path(__file__).resolve().parent / "fixtures" / "orca611"


def _reading(value):
    return orca_read.Reading(value=value, unit="", source="text", note="test")


def _geom_rows(xyz_path: Path):
    lines = xyz_path.read_text(encoding="utf-8", errors="replace").splitlines()
    count = int(lines[0].split()[0])
    return [l.split() for l in lines[2:2 + count] if l.strip()]


def _synthetic_run(tmp_path: Path, base: str = "input") -> Path:
    """A hand-written run: .out marker + ORCA-style final-geometry xyz."""
    (tmp_path / f"{base}.out").write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
    xyz = tmp_path / f"{base}.xyz"
    xyz.write_text("3\nCoordinates from ORCA-job input E -1.0\n"
                   "H 0.0 0.0 0.0\n"
                   "O 0.0 0.0 0.97\n"
                   "H 0.0 0.76 0.58\n")
    return xyz


def test_imag_read_xyz_via_orca_read(tmp_path):
    """imag._read_xyz routes the final geometry through orca_read.read."""
    xyz = _synthetic_run(tmp_path)
    with mock.patch.object(orca_read, "read",
                           side_effect=lambda wd, bn, q, **kw:
                           _reading(_geom_rows(xyz))) as mock_read:
        atoms = imag._read_xyz(tmp_path / "input.xyz")
    assert mock_read.call_count >= 1
    assert atoms == [("H", 0.0, 0.0, 0.0),
                     ("O", 0.0, 0.0, 0.97),
                     ("H", 0.0, 0.76, 0.58)]


def test_imag_read_xyz_text_fallback_value(tmp_path):
    """Without OPI (no property files), orca_read's text path answers and
    imag._read_xyz yields the same atoms as the xyz carries."""
    xyz = _synthetic_run(tmp_path)
    atoms = imag._read_xyz(xyz)
    assert atoms == [("H", 0.0, 0.0, 0.0),
                     ("O", 0.0, 0.0, 0.97),
                     ("H", 0.0, 0.76, 0.58)]


def test_recovery_continue_from_own_geometry_via_orca_read(tmp_path):
    """Recovery's continue-from-own-geometry reads via orca_read.read."""
    xyz = _synthetic_run(tmp_path)
    items = [orca_input.Item(
        "geometry",
        ["* xyz 0 1\n", "  H  0.00000000  0.00000000  0.00000000\n",
         "  O  0.00000000  0.00000000  0.97000000\n",
         "  H  0.00000000  0.76000000  0.58000000\n", "*\n"],
        changed=False)]
    with mock.patch.object(orca_read, "read",
                           side_effect=lambda wd, bn, q, **kw:
                           _reading(_geom_rows(xyz))) as mock_read:
        orca_recovery.OrcaInputModifier._continue_from_own_geometry(
            items, tmp_path / "input.xyz")
    assert mock_read.call_count >= 1
    rows = items[0].lines[1:-1]
    assert rows == ["  H    0.00000000    0.00000000    0.00000000\n",
                    "  O    0.00000000    0.00000000    0.97000000\n",
                    "  H    0.00000000    0.76000000    0.58000000\n"]


def test_recovery_text_fallback_value(tmp_path):
    """Without OPI the text path answers; the same rows are injected."""
    _synthetic_run(tmp_path)
    items = [orca_input.Item(
        "geometry",
        ["* xyz 0 1\n", "  H  0.0 0.0 0.0\n", "  O  0.0 0.0 0.97\n",
         "  H  0.0 0.76 0.58\n", "*\n"], changed=False)]
    orca_recovery.OrcaInputModifier._continue_from_own_geometry(
        items, tmp_path / "input.xyz")
    assert [r.split()[0] for r in items[0].lines[1:-1]] == ["H", "O", "H"]
    assert "0.97000000" in items[0].lines[2]


def test_orca_read_xyz_atoms_via_orca_read(tmp_path):
    """The resume path's _read_xyz_atoms routes through orca_read.read."""
    xyz = _synthetic_run(tmp_path)
    with mock.patch.object(orca_read, "read",
                           side_effect=lambda wd, bn, q, **kw:
                           _reading(_geom_rows(xyz))) as mock_read:
        rows = orca._read_xyz_atoms(tmp_path / "input.xyz")
    assert mock_read.call_count >= 1
    assert rows is not None and len(rows) == 3


def test_orca_reader_text_restores_old_path(tmp_path):
    """orca.reader=text bypasses the new routing entirely: the original
    xyz parsing answers, no orca_read.read call for the geometry."""
    xyz = _synthetic_run(tmp_path)
    settings = {"orca": {"reader": "text"}}
    with mock.patch("delfin.orca_read.load_settings",
                    return_value=settings, create=True), \
            mock.patch.object(orca_read, "read", side_effect=AssertionError(
                "orca.reader=text must not route geometry through orca_read"
            )) as mock_read:
        atoms = imag._read_xyz(xyz)
    assert atoms == [("H", 0.0, 0.0, 0.0),
                     ("O", 0.0, 0.0, 0.97),
                     ("H", 0.0, 0.76, 0.58)]
    assert mock_read.call_count == 0


def test_orca_read_xyz_atoms_unusable_file(tmp_path):
    """A malformed xyz yields None in the resume path, as before."""
    (tmp_path / "input.out").write_text("# SYNTHETIC test fixture\n")
    (tmp_path / "input.xyz").write_text("2\ncomment\nH 0 0 0\n")
    atoms = orca._read_xyz_atoms(tmp_path / "input.xyz")
    assert atoms is None
