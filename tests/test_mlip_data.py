"""MLIP training-data collection from DELFIN's own xTB/ORCA runs.

mlip_data (delfin/ml_train/mlip_data.py) reads energies, forces and stress
out of finished DELFIN runs into the extxyz format a fine-tune consumes,
recording method/basis/program per frame and converting to eV / eV/A.
Tests here write tiny SYNTHETIC fixtures in the test (never real outputs).
"""

from __future__ import annotations

import numpy as np
import pytest

from delfin.ml_train import mlip_data as md


def _write_engrad(path, n_atoms, energy_eh, grad_eh_bohr):
    """ORCA/xTB-style .engrad: n, energy Eh, 3*n gradient (Eh/bohr)."""
    lines = [
        "#",
        "# Number of atoms: must match the XYZ",
        "#",
        f"{n_atoms}",
        "#",
        "# The current total energy in Eh",
        "#",
        f"{energy_eh:.12f}",
        "#",
        "# The current gradient in Eh/bohr: Atom1X, Atom1Y, Atom1Z ...",
        "#",
    ]
    for g in np.asarray(grad_eh_bohr, dtype=float).reshape(-1):
        lines.append(f"{g:.12f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_xtb_run(dir_path, energy_eh=-10.0, grad=None):
    dir_path.mkdir(parents=True, exist_ok=True)
    # xTB writes <base>.engrad and <base>.xyz next to the run
    grad = np.array([[0.01, 0.0, 0.0],
                     [-0.01, 0.0, 0.0],
                     [0.0, 0.0, 0.0]]) if grad is None else grad
    _write_engrad(dir_path / "mol.engrad", 3, energy_eh, grad)
    (dir_path / "mol.xyz").write_text(
        "3\nxtb run\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n",
        encoding="utf-8",
    )
    (dir_path / "mol.out").write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "xtb: GFN2-xTB\nnormally converged\n",
        encoding="utf-8",
    )
    return dir_path


# ------------------------------------------------------------------ units

def test_unit_conversion_constants_are_physics():
    assert md.EV_PER_HARTREE == pytest.approx(27.211386245988)
    # Eh/bohr -> eV/A
    assert md.EV_PER_ANGSTROM_PER_HARTREE_PER_BOHR == pytest.approx(51.4220674763)


# ------------------------------------------------------------------ read_run

def test_read_xtb_run_reuses_engrad_and_xyz(tmp_path):
    run = _write_xtb_run(tmp_path / "xtb_run", energy_eh=-10.0)
    frame = md.read_run(run, "mol", program="xtb",
                        source="my xtb calc", licence="self")
    assert frame.program == "xtb"
    assert frame.method  # GFN2-xTB default for xtb
    assert frame.basis is None
    assert frame.energy_eV == pytest.approx(-10.0 * md.EV_PER_HARTREE)
    assert frame.forces.shape == (3, 3)
    # forces in eV/A, sign = negative gradient
    assert frame.forces[0, 0] == pytest.approx(-0.01 * md.EV_PER_ANGSTROM_PER_HARTREE_PER_BOHR)
    assert len(frame.symbols) == 3
    assert tuple(frame.symbols) == ("C", "H", "H")


def test_read_returns_energy_in_eV_and_forces_in_eV_per_angstrom(tmp_path):
    run = _write_xtb_run(tmp_path / "r", energy_eh=-5.0)
    frame = md.read_run(run, "mol", program="xtb",
                        source="own", licence="self")
    assert frame.energy_eV < 0
    assert frame.forces[:, :].shape == (3, 3)
    # energies in Hartree become eV (~27.2x) -- not left in Eh
    assert abs(frame.energy_eV) > abs(-5.0)


# ------------------------------------------------------------------ mixed method

def test_collect_dataset_refuses_mixed_methods(tmp_path):
    a = _write_xtb_run(tmp_path / "a", energy_eh=-10.0)
    b = _write_xtb_run(tmp_path / "b", energy_eh=-11.0)
    fa = md.read_run(a, "mol", program="xtb", method="GFN1-xTB",
                     source="own", licence="self")
    fb = md.read_run(b, "mol", program="xtb", method="GFN2-xTB",
                     source="own", licence="self")
    with pytest.raises(md.MethodMismatchError):
        md.collect_dataset([fa, fb], mixed=False)


def test_collect_dataset_allows_mixed_when_explicit(tmp_path):
    a = _write_xtb_run(tmp_path / "a", energy_eh=-10.0)
    b = _write_xtb_run(tmp_path / "b", energy_eh=-11.0)
    fa = md.read_run(a, "mol", program="xtb", method="GFN1-xTB",
                     source="own", licence="self")
    fb = md.read_run(b, "mol", program="xtb", method="GFN2-xTB",
                     source="own", licence="self")
    frames = md.collect_dataset([fa, fb], mixed=True)
    assert len(frames) == 2


# ---------------------------------------------------------------- extxyz

def test_frames_to_extxyz_records_properties_per_frame(tmp_path):
    run = _write_xtb_run(tmp_path / "run", energy_eh=-10.0)
    frame = md.read_run(run, "mol", program="xtb", method="GFN2-xTB",
                        source="own", licence="self")
    out = tmp_path / "data.extxyz"
    md.frames_to_extxyz([frame], out)
    text = out.read_text(encoding="utf-8")
    assert text.startswith("3\n")
    # Properties line names per-atom force columns
    assert "Properties=species:S:1:pos:R:3:forces:R:3" in text
    # energy recorded in eV in the comment line
    assert "energy=" in text
    assert "method=GFN2-xTB" in text
    assert "program=xtb" in text
    assert "source=own" in text
    assert "licence=self" in text


def test_frames_to_extxyz_force_columns_match_atoms(tmp_path):
    run = _write_xtb_run(tmp_path / "run", energy_eh=-10.0)
    frame = md.read_run(run, "mol", program="xtb", source="own", licence="self")
    out = tmp_path / "data.extxyz"
    md.frames_to_extxyz([frame], out)
    lines = out.read_text(encoding="utf-8").splitlines()
    # 3 atoms -> 3 atom lines, each with symbol + 3 pos + 3 force columns
    for line in lines[2:]:
        cols = line.split()
        assert len(cols) == 1 + 3 + 3


def test_frames_to_extxyz_records_stress_when_present(tmp_path):
    run = _write_xtb_run(tmp_path / "run", energy_eh=-10.0)
    frame = md.read_run(run, "mol", program="xtb", source="own", licence="self")
    frame = md.replace(frame, stress=[1.0, 2.0, 3.0, 0.0, 0.0, 0.0])
    out = tmp_path / "data.extxyz"
    md.frames_to_extxyz([frame], out)
    assert 'stress="1.0 2.0 3.0 0.0 0.0 0.0"' in out.read_text(encoding="utf-8")


# -------------------------------------------------------- licence / source

def test_csd_or_ccdc_source_is_refused():
    with pytest.raises(md.LicenceError):
        md.check_source_licence("csd 1278901", "CCDC")
    with pytest.raises(md.LicenceError):
        md.check_source_licence("Cambridge Structural Database", "unknown")
    with pytest.raises(md.LicenceError):
        md.check_source_licence("other", "BIG-IP")


def test_cod_cc0_and_own_are_allowed():
    md.check_source_licence("cod-123456", "CC0")
    md.check_source_licence("our DFT calc", "own")
    md.check_source_licence("our DFT calc", "self")


# ------------------------------------------------------------ hold-out guard

def test_guard_refuses_forbidden_hold_out_path(tmp_path):
    # synthetic forbidden path components, not real files
    for bad in ("LOCKBOX2", "stages_v2", "heldout_stage1", "pool_6000",
                "stage1_v2.txt"):
        path = tmp_path / "data" / bad / "train.xyz"
        assert md.refuses_hold_out(path) is True
    assert md.refuses_hold_out(tmp_path / "my" / "train.xyz") is False


def test_guard_refuses_provenance_naming_a_hold_out():
    assert md.refuses_hold_out("some/train.xyz",
                               provenance="from LOCKBOX2 collection") is True
    assert md.refuses_hold_out("some/train.xyz",
                               provenance="held-out part of pool_6000") is True
    assert md.refuses_hold_out("some/train.xyz", provenance="our own run") is False
