"""Tests for the inverse-RSS and generalised-binding features of the CO2
coordinator (delfin/co2/CO2_Coordinator6.py).

Covers: binding spec parsing (#109), Rodrigues rotation helper and scan-axis
selection (#110), and the inverse RSS control-key validation path (#108/#111).
No ORCA is executed here — main() is only tested up to the validation layer
via monkeypatched file access where needed.
"""
import os
import textwrap

import numpy as np
import pytest

from delfin.co2.CO2_Coordinator6 import (
    _is_enabled,
    _parse_binding_spec,
    rotation_about_axis,
    _rotation_axis_for_scan,
    Rz,
    METAL_SYMBOLS,
)


# ---------------------------------------------------------------------------
# Binding spec parsing (#109)
# ---------------------------------------------------------------------------

class _Atom:
    def __init__(self, symbol, position):
        self.symbol = symbol
        self.position = np.asarray(position, dtype=float)


class _Atoms:
    def __init__(self, atoms):
        self._atoms = list(atoms)

    def __len__(self):
        return len(self._atoms)

    def __iter__(self):
        return iter(self._atoms)

    def __getitem__(self, i):
        return self._atoms[i]


class _AtomsWithPositions(_Atoms):
    @property
    def positions(self):
        return np.array([a.position for a in self._atoms])


def _fake_complex():
    """Ni at origin, substrate COM on +z."""
    return _AtomsWithPositions([
        _Atom("Ni", [0.0, 0.0, 0.0]),
        _Atom("C", [0.0, 0.0, 4.0]),
        _Atom("O", [0.0, 0.0, 5.2]),
    ])


def test_binding_spec_explicit_forms():
    assert _parse_binding_spec({"binding": "atom-on:O"}) == "atom-on:O"
    assert _parse_binding_spec({"binding": " bond-on:2-3 "}) == "bond-on:2-3"
    assert _parse_binding_spec({"binding": "end-on:C"}) == "end-on:C"
    assert _parse_binding_spec({"binding": "side-on"}) == "side-on"


def test_binding_spec_legacy_mode_mapping():
    assert _parse_binding_spec({"mode": "side-on"}) == "side-on"
    assert _parse_binding_spec({"mode": "end-on"}) == "end-on"


def test_binding_spec_none_when_unset():
    assert _parse_binding_spec({}) is None
    assert _parse_binding_spec({"binding": "", "mode": ""}) is None


def test_binding_spec_binding_overrides_mode():
    assert _parse_binding_spec({"binding": "atom-on:O", "mode": "end-on"}) == "atom-on:O"


def test_is_enabled_truthy_forms():
    for v in (True, "true", "True", "yes", "1", "on"):
        assert _is_enabled(v) is True, v
    for v in (False, None, "", "no", "0", "off", "False"):
        assert _is_enabled(v) is False, v


# ---------------------------------------------------------------------------
# Rodrigues rotation helper (#110)
# ---------------------------------------------------------------------------

def test_rotation_about_axis_z_equals_rz():
    rng = np.random.default_rng(42)
    pts = rng.normal(size=(5, 3))
    R = rotation_about_axis([0.0, 0.0, 1.0], 37.0)
    assert np.allclose(pts @ R.T, pts @ Rz(37.0).T, atol=1e-10)


def test_rotation_about_axis_preserves_axis_points():
    R = rotation_about_axis([1.0, 2.0, -0.5], 123.0)
    v = np.array([2.0, 4.0, -1.0])  # parallel to the axis
    assert np.allclose(v @ R.T, v, atol=1e-10)


def test_rotation_about_axis_is_orthonormal():
    R = rotation_about_axis([1.0, 0.0, 1.0], 88.0)
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-10)
    assert abs(np.linalg.det(R) - 1.0) < 1e-10


def test_rotation_about_axis_zero_axis_raises():
    with pytest.raises(ValueError):
        rotation_about_axis([0.0, 0.0, 0.0], 10.0)


# ---------------------------------------------------------------------------
# Scan-axis selection (#110)
# ---------------------------------------------------------------------------

def test_rotation_axis_for_scan_metal_to_substrate_com():
    atoms = _fake_complex()
    axis = _rotation_axis_for_scan(atoms, [1, 2])
    expected = np.array([0.0, 0.0, 1.0])  # COM on +z, metal at origin
    assert axis is not None and np.allclose(axis, expected, atol=1e-10)


def test_rotation_axis_for_scan_degenerate_returns_none():
    atoms = _AtomsWithPositions([
        _Atom("Ni", [0.0, 0.0, 0.0]),
        _Atom("C", [0.0, 0.0, 0.0]),  # substrate COM exactly at the metal
    ])
    assert _rotation_axis_for_scan(atoms, [1]) is None


def test_rotation_axis_for_scan_is_normalised():
    atoms = _AtomsWithPositions([
        _Atom("Fe", [1.0, 1.0, 1.0]),
        _Atom("C", [4.0, 1.0, 1.0]),
    ])
    axis = _rotation_axis_for_scan(atoms, [1])
    assert axis is not None and abs(np.linalg.norm(axis) - 1.0) < 1e-10
    assert np.allclose(axis, [1.0, 0.0, 0.0], atol=1e-10)


# ---------------------------------------------------------------------------
# Inverse RSS validation in main() (#108/#111)
# ---------------------------------------------------------------------------

def test_main_rejects_scan_dissoc_without_adduct(tmp_path, monkeypatch):
    from delfin.co2 import CO2_Coordinator6 as mod
    control = tmp_path / "CONTROL.txt"
    control.write_text(textwrap.dedent("""\
        xyz=complex.xyz
        out=out.xyz
        scan_dissoc=true
        dissoc_distance=6.0
    """))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(mod, "_minimal_read_control_file",
                        lambda: mod._read_control_file(str(control))
                        if hasattr(mod, "_read_control_file")
                        else {"scan_dissoc": "true", "dissoc_distance": "6.0"})
    with pytest.raises(ValueError, match="adduct_xyz"):
        mod.main()


def test_main_rejects_scan_dissoc_with_missing_adduct_file(tmp_path, monkeypatch):
    from delfin.co2 import CO2_Coordinator6 as mod
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(mod, "_minimal_read_control_file",
                        lambda: {"scan_dissoc": "true",
                                 "adduct_xyz": str(tmp_path / "nope.xyz"),
                                 "dissoc_distance": "6.0"})
    with pytest.raises(FileNotFoundError):
        mod.main()


def test_main_rejects_invalid_dissoc_distance(tmp_path, monkeypatch):
    from delfin.co2 import CO2_Coordinator6 as mod
    adduct = tmp_path / "adduct.xyz"
    adduct.write_text("2\n\nNi 0 0 0\nC 0 0 2.0\n")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(mod, "_minimal_read_control_file",
                        lambda: {"scan_dissoc": "true",
                                 "adduct_xyz": str(adduct),
                                 "dissoc_distance": "-1.0"})
    with pytest.raises(ValueError, match="dissoc_distance"):
        mod.main()


def test_main_inverse_rss_requires_substrate_atom_index(tmp_path, monkeypatch):
    from delfin.co2 import CO2_Coordinator6 as mod
    adduct = tmp_path / "adduct.xyz"
    adduct.write_text("2\n\nNi 0 0 0\nC 0 0 2.0\n")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(mod, "_minimal_read_control_file",
                        lambda: {"scan_dissoc": "true",
                                 "adduct_xyz": str(adduct),
                                 "dissoc_distance": "6.0"})
    with pytest.raises(ValueError, match="substrate_atom_index"):
        mod.main()


def test_main_inverse_rss_rejects_bad_substrate_atom_index(tmp_path, monkeypatch):
    from delfin.co2 import CO2_Coordinator6 as mod
    adduct = tmp_path / "adduct.xyz"
    adduct.write_text("2\n\nNi 0 0 0\nC 0 0 2.0\n")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(mod, "_minimal_read_control_file",
                        lambda: {"scan_dissoc": "true",
                                 "adduct_xyz": str(adduct),
                                 "dissoc_distance": "6.0",
                                 "substrate_atom_index": "5"})
    with pytest.raises(ValueError, match="substrate_atom_index=5"):
        mod.main()


def test_main_inverse_rss_writes_scan_input_from_adduct(tmp_path, monkeypatch):
    """End-to-end (mocked ORCA): scan input must go from r0=2.0 A outward to 6.0 A."""
    from delfin.co2 import CO2_Coordinator6 as mod
    adduct = tmp_path / "adduct.xyz"
    adduct.write_text("2\n\nNi 0 0 0\nC 0 0 2.0\n")
    monkeypatch.chdir(tmp_path)

    captured = {}

    def fake_write_orca_input_and_run(atoms, xyz_path, metal_index, co2_c_index,
                                      start_distance, end_distance=1.7, steps=5, **kw):
        captured["start"] = start_distance
        captured["end"] = end_distance
        captured["steps"] = steps
        captured["metal"] = metal_index
        captured["c"] = co2_c_index
        # Simulate the expected output location so plot_scan_result is callable
        os.makedirs("relaxed_surface_scan", exist_ok=True)
        with open(os.path.join("relaxed_surface_scan", "scan.relaxscanact.dat"), "w") as f:
            f.write("# mocked\n")
        return "mocked.out"

    monkeypatch.setattr(mod, "write_orca_input_and_run", fake_write_orca_input_and_run)
    monkeypatch.setattr(mod, "plot_scan_result", lambda p: None)
    monkeypatch.setattr(mod, "_minimal_read_control_file",
                        lambda: {"scan_dissoc": "true",
                                 "adduct_xyz": str(adduct),
                                 "dissoc_distance": "6.0",
                                 "scan_steps": "20",
                                 "substrate_atom_index": "1",
                                 "PAL": "4", "maxcore": "1000",
                                 "charge": "-2", "multiplicity": "1"})
    mod.main()
    assert captured["start"] == pytest.approx(2.0)
    assert captured["end"] == pytest.approx(6.0)
    assert captured["steps"] == 20
    assert captured["metal"] == 0
    assert captured["c"] == 1
