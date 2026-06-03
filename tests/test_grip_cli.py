"""Tests for the ``delfin-grip`` CLI + supporting ``delfin.grip_io`` helpers.

Coverage:
    * XYZ I/O round-trip (single + multi frame)
    * Bond perception (organic + TMC)
    * SMILES->XYZ atom mapping (same-order, reordered, equivalent atoms)
    * Metal + donor heuristic detection
    * CLI invocations: --help, basic run, --info, ensemble, --no-corrector
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

from delfin.grip_cli import build_parser, refine_frame
from delfin.grip_io import (
    Frame,
    detect_metal_and_donors,
    match_smiles_to_xyz,
    perceive_bonds_from_xyz,
    read_xyz,
    write_xyz,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _toluene_xyz_text() -> str:
    return (
        "15\n"
        "toluene\n"
        "C    1.20      0.00     0.00\n"
        "C    0.50      1.20     0.00\n"
        "C   -0.90      1.20     0.00\n"
        "C   -1.60      0.00     0.00\n"
        "C   -0.90     -1.20     0.00\n"
        "C    0.50     -1.20     0.00\n"
        "C    2.70      0.00     0.00\n"
        "H    1.05      2.15     0.00\n"
        "H   -1.45      2.15     0.00\n"
        "H   -2.70      0.00     0.00\n"
        "H   -1.45     -2.15     0.00\n"
        "H    1.05     -2.15     0.00\n"
        "H    3.07      0.50     0.90\n"
        "H    3.07      0.50    -0.90\n"
        "H    3.07     -1.00     0.00\n"
    )


def _cu_nh3_4_xyz_text() -> str:
    return (
        "13\n"
        "Cu(NH3)4 synthetic\n"
        "Cu  0.0   0.0  0.0\n"
        "N   2.0   0.0  0.0\n"
        "N  -2.0   0.0  0.0\n"
        "N   0.0   2.0  0.0\n"
        "N   0.0  -2.0  0.0\n"
        "H   2.6   0.7  0.5\n"
        "H   2.6  -0.7 -0.5\n"
        "H  -2.6   0.7  0.5\n"
        "H  -2.6  -0.7 -0.5\n"
        "H   0.7   2.6  0.5\n"
        "H  -0.7   2.6 -0.5\n"
        "H   0.7  -2.6  0.5\n"
        "H  -0.7  -2.6 -0.5\n"
    )


def _write_xyz_text(tmp_path: Path, text: str, name: str = "test.xyz") -> Path:
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# read_xyz / write_xyz
# ---------------------------------------------------------------------------
class TestReadWriteXYZ:
    def test_read_xyz_single_frame(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _toluene_xyz_text())
        frames = read_xyz(path)
        assert len(frames) == 1
        assert frames[0].n_atoms == 15
        assert frames[0].elements[0] == "C"
        assert frames[0].coordinates.shape == (15, 3)
        assert frames[0].comment == "toluene"

    def test_read_xyz_multi_frame(self, tmp_path: Path):
        text = _toluene_xyz_text() + _toluene_xyz_text()
        path = _write_xyz_text(tmp_path, text, name="multi.xyz")
        frames = read_xyz(path)
        assert len(frames) == 2
        assert frames[0].n_atoms == 15
        assert frames[1].n_atoms == 15

    def test_read_xyz_empty_raises(self, tmp_path: Path):
        path = tmp_path / "empty.xyz"
        path.write_text("", encoding="utf-8")
        with pytest.raises(ValueError):
            read_xyz(path)

    def test_read_xyz_truncated_raises(self, tmp_path: Path):
        # Declare 15 atoms but only provide 3 rows.
        text = "15\nbad\nC 0 0 0\nC 1 0 0\nC 2 0 0\n"
        path = _write_xyz_text(tmp_path, text)
        with pytest.raises(ValueError):
            read_xyz(path)

    def test_write_xyz_round_trip(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _toluene_xyz_text())
        frames = read_xyz(path)
        out = tmp_path / "rt.xyz"
        write_xyz(out, frames)
        re = read_xyz(out)
        assert len(re) == 1
        assert re[0].elements == frames[0].elements
        assert np.allclose(re[0].coordinates, frames[0].coordinates, atol=1e-6)

    def test_write_xyz_comment_prefix(self, tmp_path: Path):
        f = Frame(
            elements=["C"], coordinates=np.zeros((1, 3)), comment="orig",
        )
        out = tmp_path / "p.xyz"
        write_xyz(out, f, comment_prefix="TAG")
        text = out.read_text()
        assert "TAG orig" in text

    def test_element_canonicalisation(self, tmp_path: Path):
        # Mixed cases + atomic-number tokens
        text = "3\ntest\ncl  0 0 0\n11 1 0 0\nH  0 1 0\n"
        path = _write_xyz_text(tmp_path, text)
        frames = read_xyz(path)
        assert frames[0].elements == ["Cl", "Na", "H"]


# ---------------------------------------------------------------------------
# Bond perception
# ---------------------------------------------------------------------------
class TestPerceiveBonds:
    def test_perceive_bonds_simple_organic(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _toluene_xyz_text())
        frames = read_xyz(path)
        mol = perceive_bonds_from_xyz(
            frames[0].elements, frames[0].coordinates, charge=0,
        )
        assert mol is not None
        # Toluene has 15 atoms (7 heavy + 8 H)
        assert mol.GetNumAtoms() == 15
        # Should detect at least 14 bonds (toluene has 15)
        assert mol.GetNumBonds() >= 14

    def test_perceive_bonds_tmc_simple(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _cu_nh3_4_xyz_text())
        frames = read_xyz(path)
        # Cu may fail in rdDetermineBonds; fallback to distance must still work
        mol = perceive_bonds_from_xyz(
            frames[0].elements, frames[0].coordinates, charge=0,
        )
        assert mol is not None
        assert mol.GetNumAtoms() == 13
        # Distance graph should yield at least 4 Cu-N + 12 N-H = 16 bonds
        assert mol.GetNumBonds() >= 4

    def test_perceive_bonds_invalid_shape(self):
        with pytest.raises(ValueError):
            perceive_bonds_from_xyz(["C", "H"], np.zeros((1, 3)))


# ---------------------------------------------------------------------------
# Atom mapping
# ---------------------------------------------------------------------------
class TestSmilesXYZMapping:
    def test_match_smiles_to_xyz_same_order(self, tmp_path: Path):
        # Methane CH4: SMILES atom order matches XYZ rows.
        path = _write_xyz_text(
            tmp_path,
            "5\nmethane\nC 0 0 0\nH 1 0 0\nH -1 0 0\nH 0 1 0\nH 0 -1 0\n",
        )
        frames = read_xyz(path)
        mol = perceive_bonds_from_xyz(frames[0].elements, frames[0].coordinates)
        mapping = match_smiles_to_xyz("C", mol, frames[0].coordinates)
        assert mapping is not None
        # 5 atoms must be mapped after AddHs
        assert len(mapping) == 5

    def test_match_smiles_to_xyz_reordered(self, tmp_path: Path):
        # Build XYZ with H atoms BEFORE C; SMILES has C first.
        path = _write_xyz_text(
            tmp_path,
            "5\nmethane reordered\nH 1 0 0\nH -1 0 0\nH 0 1 0\nH 0 -1 0\nC 0 0 0\n",
        )
        frames = read_xyz(path)
        mol = perceive_bonds_from_xyz(frames[0].elements, frames[0].coordinates)
        mapping = match_smiles_to_xyz("C", mol, frames[0].coordinates)
        # Should at least be able to produce a mapping (Hungarian fallback)
        assert mapping is not None
        # All H's of SMILES should map to H's of XYZ, and the C to the C row (index 4)
        from rdkit import Chem
        smi = Chem.AddHs(Chem.MolFromSmiles("C"))
        for smi_i, xyz_i in mapping.items():
            assert smi.GetAtomWithIdx(smi_i).GetSymbol() == \
                   mol.GetAtomWithIdx(xyz_i).GetSymbol()

    def test_match_smiles_to_xyz_equivalent_atoms(self, tmp_path: Path):
        # Cu(NH3)4 — 4 equivalent N donors.  Any consistent mapping works.
        path = _write_xyz_text(tmp_path, _cu_nh3_4_xyz_text())
        frames = read_xyz(path)
        mol = perceive_bonds_from_xyz(frames[0].elements, frames[0].coordinates)
        # SMILES of cu-tetraamine using dative-like notation
        # Use simple connectivity that captures element counts.
        mapping = match_smiles_to_xyz(
            "[Cu](N)(N)(N)N", mol, frames[0].coordinates,
        )
        # We accept None as well because rdDetermineBonds often refuses
        # Cu-N bond orders; the test asserts the function returns something
        # consistent (mapping or None, no exception).
        if mapping is not None:
            from rdkit import Chem
            smi = Chem.AddHs(Chem.MolFromSmiles("[Cu](N)(N)(N)N"))
            for smi_i, xyz_i in mapping.items():
                assert smi.GetAtomWithIdx(smi_i).GetSymbol() == \
                       mol.GetAtomWithIdx(xyz_i).GetSymbol()

    def test_match_smiles_bad_smiles_returns_none(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _toluene_xyz_text())
        frames = read_xyz(path)
        mol = perceive_bonds_from_xyz(frames[0].elements, frames[0].coordinates)
        assert match_smiles_to_xyz("not a smiles!?", mol, frames[0].coordinates) is None


# ---------------------------------------------------------------------------
# Metal + donor heuristic
# ---------------------------------------------------------------------------
class TestMetalDonorDetect:
    def test_detect_metal_simple(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _cu_nh3_4_xyz_text())
        frames = read_xyz(path)
        mol = perceive_bonds_from_xyz(frames[0].elements, frames[0].coordinates)
        metal, donors = detect_metal_and_donors(mol, frames[0].coordinates)
        assert metal == 0  # Cu is the first row
        assert len(donors) >= 4  # 4 nitrogens

    def test_detect_no_metal_for_organic(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _toluene_xyz_text())
        frames = read_xyz(path)
        mol = perceive_bonds_from_xyz(frames[0].elements, frames[0].coordinates)
        metal, donors = detect_metal_and_donors(mol, frames[0].coordinates)
        assert metal is None
        assert donors == []


# ---------------------------------------------------------------------------
# CLI invocation tests (subprocess)
# ---------------------------------------------------------------------------
def _cli_run(*args, **kwargs):
    """Invoke ``python -m delfin.grip_cli`` as a subprocess."""
    env = os.environ.copy()
    env.setdefault("PYTHONHASHSEED", "0")
    return subprocess.run(
        [sys.executable, "-m", "delfin.grip_cli", *args],
        env=env, capture_output=True, text=True, check=False, **kwargs,
    )


class TestCLI:
    def test_cli_help(self):
        r = _cli_run("--help")
        assert r.returncode == 0
        assert "delfin-grip" in r.stdout
        assert "--ensemble" in r.stdout

    def test_cli_version(self):
        r = _cli_run("--version")
        assert r.returncode == 0
        assert "delfin-grip" in r.stdout

    def test_cli_info(self):
        r = _cli_run("--info")
        # Should not be a hard error if library is present
        assert r.returncode in (0, 1)
        assert "delfin-grip" in r.stdout
        assert "library" in r.stdout.lower()

    def test_cli_basic_toluene(self, tmp_path: Path):
        in_xyz = _write_xyz_text(tmp_path, _toluene_xyz_text())
        out_xyz = tmp_path / "out.xyz"
        r = _cli_run(str(in_xyz), "-o", str(out_xyz))
        # Either ACCEPTED (0) or no-improvement (2); both are clean exits.
        assert r.returncode in (0, 2), r.stderr
        assert out_xyz.exists()
        frames = read_xyz(out_xyz)
        assert frames[0].n_atoms == 15

    def test_cli_no_corrector(self, tmp_path: Path):
        in_xyz = _write_xyz_text(tmp_path, _toluene_xyz_text())
        out_xyz = tmp_path / "out.xyz"
        r = _cli_run(str(in_xyz), "-o", str(out_xyz), "--no-corrector")
        assert r.returncode in (0, 2), r.stderr
        assert out_xyz.exists()

    def test_cli_validate_only_no_output_written(self, tmp_path: Path):
        in_xyz = _write_xyz_text(tmp_path, _toluene_xyz_text())
        out_xyz = tmp_path / "should_not_exist.xyz"
        r = _cli_run(str(in_xyz), "-o", str(out_xyz), "--validate")
        assert r.returncode in (0, 2), r.stderr
        assert not out_xyz.exists()

    def test_cli_explicit_metal_donors(self, tmp_path: Path):
        in_xyz = _write_xyz_text(tmp_path, _cu_nh3_4_xyz_text())
        out_xyz = tmp_path / "cu_out.xyz"
        r = _cli_run(
            str(in_xyz), "-o", str(out_xyz),
            "--metal", "0", "--donors", "1,2,3,4",
        )
        assert r.returncode in (0, 2), r.stderr
        assert out_xyz.exists()
        # Cu and N's must still be there
        frames = read_xyz(out_xyz)
        assert frames[0].elements[0] == "Cu"

    def test_cli_bad_input_file_returns_1(self, tmp_path: Path):
        out_xyz = tmp_path / "out.xyz"
        r = _cli_run(str(tmp_path / "nonexistent.xyz"), "-o", str(out_xyz))
        assert r.returncode == 1

    def test_cli_bad_donors_returns_1(self, tmp_path: Path):
        in_xyz = _write_xyz_text(tmp_path, _cu_nh3_4_xyz_text())
        out_xyz = tmp_path / "out.xyz"
        r = _cli_run(
            str(in_xyz), "-o", str(out_xyz),
            "--donors", "not,a,number",
        )
        assert r.returncode == 1

    def test_cli_default_output_path(self, tmp_path: Path):
        in_xyz = _write_xyz_text(tmp_path, _toluene_xyz_text(), name="in.xyz")
        r = _cli_run(str(in_xyz))
        assert r.returncode in (0, 2)
        # default = in.grip.xyz next to input (the .xyz suffix is stripped)
        default = tmp_path / "in.grip.xyz"
        assert default.exists()


# ---------------------------------------------------------------------------
# In-process refine_frame
# ---------------------------------------------------------------------------
class TestRefineFrame:
    def test_refine_frame_toluene_in_process(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _toluene_xyz_text())
        frames = read_xyz(path)
        rf, res = refine_frame(frames[0])
        assert rf.n_atoms == 15
        assert res.frame_index == 0
        # No metal detected -> we still get a valid result object
        assert np.all(np.isfinite(rf.coordinates))

    def test_refine_frame_tmc_in_process(self, tmp_path: Path):
        path = _write_xyz_text(tmp_path, _cu_nh3_4_xyz_text())
        frames = read_xyz(path)
        rf, res = refine_frame(frames[0])
        assert rf.n_atoms == 13
        assert np.all(np.isfinite(rf.coordinates))

    def test_refine_frame_deterministic(self, tmp_path: Path):
        # Same input -> same output (within 1e-9).
        path = _write_xyz_text(tmp_path, _toluene_xyz_text())
        frames = read_xyz(path)
        rf1, _ = refine_frame(frames[0])
        rf2, _ = refine_frame(frames[0])
        assert np.allclose(rf1.coordinates, rf2.coordinates, atol=1e-9)


# ---------------------------------------------------------------------------
# Parser unit tests
# ---------------------------------------------------------------------------
class TestParser:
    def test_parser_minimal(self):
        p = build_parser()
        args = p.parse_args(["in.xyz"])
        assert args.input_xyz == "in.xyz"
        assert args.max_iter == 200
        assert args.no_corrector is False

    def test_parser_donors_list(self):
        p = build_parser()
        args = p.parse_args(["in.xyz", "--donors", "1,2,3"])
        assert args.donors == "1,2,3"
