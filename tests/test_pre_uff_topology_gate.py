"""Tests for Fix B: pre-UFF M-D topology gate (DELFIN_PRE_UFF_TOPOLOGY_GATE).

Universal mechanism: filter conformers whose bonded M-D distances lie outside
[0.80, 1.20] times the element-pair ideal. Default OFF; opt-in.

Failure pattern being fixed: unparam-TM HARD-fallback freezes broken
geometries that pass downstream filters because nothing checks the absolute
M-D distance against a numeric reference.
"""
from __future__ import annotations

import math
import os
import pytest

pytest.importorskip("rdkit")
from rdkit import Chem
from rdkit.Geometry import Point3D

from delfin import smiles_converter as sc


def _build_md_mol(metal_sym: str, donor_sym: str, md_dist: float):
    """Build a single-M, single-D complex at the given distance."""
    m = Chem.RWMol()
    m.AddAtom(Chem.Atom(metal_sym))
    m.AddAtom(Chem.Atom(donor_sym))
    m.AddBond(0, 1, Chem.BondType.SINGLE)
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(md_dist, 0.0, 0.0))
    m.AddConformer(conf, assignId=True)
    return m.GetMol()


def test_helper_default_off():
    """Without env-flag, _pre_uff_topology_gate_enabled returns False."""
    os.environ.pop("DELFIN_PRE_UFF_TOPOLOGY_GATE", None)
    assert sc._pre_uff_topology_gate_enabled() is False


def test_helper_env_flag_truthy(monkeypatch):
    for val in ("1", "true", "TRUE", "yes", "on"):
        monkeypatch.setenv("DELFIN_PRE_UFF_TOPOLOGY_GATE", val)
        assert sc._pre_uff_topology_gate_enabled() is True


def test_helper_env_flag_falsy(monkeypatch):
    for val in ("0", "false", "off", ""):
        monkeypatch.setenv("DELFIN_PRE_UFF_TOPOLOGY_GATE", val)
        assert sc._pre_uff_topology_gate_enabled() is False


def test_in_tolerance_accepts_ideal():
    """Exact ideal Ni-N is accepted."""
    target = float(sc._get_ml_bond_length("Ni", "N"))
    mol = _build_md_mol("Ni", "N", target)
    assert sc._md_distance_in_tolerance(mol, 0) is True


def test_in_tolerance_accepts_small_deviation():
    """+/-15% is within the default [0.80, 1.20] window."""
    target = float(sc._get_ml_bond_length("Ni", "N"))
    for factor in (0.85, 0.95, 1.05, 1.15):
        mol = _build_md_mol("Ni", "N", target * factor)
        assert sc._md_distance_in_tolerance(mol, 0) is True, (
            f"factor={factor} should be inside default tolerance"
        )


def test_in_tolerance_rejects_collapsed():
    """Ni-N at 1.40 A (~68% of ideal) is rejected."""
    mol = _build_md_mol("Ni", "N", 1.40)
    assert sc._md_distance_in_tolerance(mol, 0) is False


def test_in_tolerance_rejects_stretched():
    """Ni-N at 3.00 A (~146% of ideal) is rejected."""
    mol = _build_md_mol("Ni", "N", 3.00)
    assert sc._md_distance_in_tolerance(mol, 0) is False


def test_in_tolerance_skips_hydrogens():
    """H donors are ignored (not actually coordinating)."""
    m = Chem.RWMol()
    m.AddAtom(Chem.Atom("Ni"))
    m.AddAtom(Chem.Atom("H"))
    m.AddBond(0, 1, Chem.BondType.SINGLE)
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.20, 0.0, 0.0))  # very short M-H
    m.AddConformer(conf, assignId=True)
    mol = m.GetMol()
    # H not gated -> passes
    assert sc._md_distance_in_tolerance(mol, 0) is True


def test_in_tolerance_passes_when_no_metal():
    """Pure-organic mols always pass."""
    m = Chem.RWMol()
    m.AddAtom(Chem.Atom("C"))
    m.AddAtom(Chem.Atom("O"))
    m.AddBond(0, 1, Chem.BondType.SINGLE)
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.43, 0.0, 0.0))
    m.AddConformer(conf, assignId=True)
    mol = m.GetMol()
    assert sc._md_distance_in_tolerance(mol, 0) is True


def test_in_tolerance_custom_tolerance():
    """Caller can tighten the window."""
    target = float(sc._get_ml_bond_length("Ni", "N"))
    mol = _build_md_mol("Ni", "N", target * 1.10)
    # Default passes
    assert sc._md_distance_in_tolerance(mol, 0) is True
    # Tighter window rejects
    assert (
        sc._md_distance_in_tolerance(mol, 0, rel_low=0.95, rel_high=1.05)
        is False
    )


def test_in_tolerance_universal_across_metals():
    """Universality: same logic applies to any metal-donor pair."""
    pairs = [("Pd", "N"), ("Pt", "Cl"), ("Cu", "S"), ("Fe", "P"), ("Co", "O")]
    for m_sym, d_sym in pairs:
        target = float(sc._get_ml_bond_length(m_sym, d_sym))
        mol_ok = _build_md_mol(m_sym, d_sym, target)
        mol_short = _build_md_mol(m_sym, d_sym, target * 0.5)
        mol_long = _build_md_mol(m_sym, d_sym, target * 1.5)
        assert sc._md_distance_in_tolerance(mol_ok, 0) is True, (
            f"{m_sym}-{d_sym} ideal should pass"
        )
        assert sc._md_distance_in_tolerance(mol_short, 0) is False, (
            f"{m_sym}-{d_sym} half should fail"
        )
        assert sc._md_distance_in_tolerance(mol_long, 0) is False, (
            f"{m_sym}-{d_sym} 1.5x should fail"
        )
