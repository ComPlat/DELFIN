"""Tests for Fix A: pre-UFF M-D snap (DELFIN_PRE_UFF_MD_SNAP).

Universal mechanism: snap every bonded M-D distance to the ideal element-pair
length by rigidly translating the donor's BFS fragment along the (D - M)
direction. Default OFF.

Failure pattern being fixed: ETKDG embeds M-D bonds at ~1.5 A (organic-bond
length); for explicit-charge transition metals, OB UFF triggers an unparam-TM
HARD-fallback and freezes the geometry, preserving the broken distance.
"""
from __future__ import annotations

import math
import os
import pytest

pytest.importorskip("rdkit")
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from delfin import smiles_converter as sc


def _build_simple_complex(metal_sym: str, donor_syms, md_dist: float):
    """Build a tiny RDKit Mol with one metal at origin and one donor each.

    Uses RWMol so we can attach arbitrary metal-donor bonds without SMILES
    parser interference.
    """
    m = Chem.RWMol()
    m.AddAtom(Chem.Atom(metal_sym))  # idx 0 = metal
    donor_idxs = []
    for d in donor_syms:
        idx = m.AddAtom(Chem.Atom(d))
        donor_idxs.append(idx)
        m.AddBond(0, idx, Chem.BondType.SINGLE)
    conf = Chem.Conformer(m.GetNumAtoms())
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    # Place donors along x, y, z axes at the (broken) distance
    axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]
    for k, di in enumerate(donor_idxs):
        ax, ay, az = axes[k % len(axes)]
        conf.SetAtomPosition(di, Point3D(ax * md_dist, ay * md_dist, az * md_dist))
    m.AddConformer(conf, assignId=True)
    return m.GetMol(), donor_idxs


def test_helper_default_off():
    """Without env-flag, _pre_uff_md_snap_enabled returns False."""
    os.environ.pop("DELFIN_PRE_UFF_MD_SNAP", None)
    assert sc._pre_uff_md_snap_enabled() is False


def test_helper_env_flag_truthy_values(monkeypatch):
    """1/true/yes/on all enable the flag."""
    for val in ("1", "true", "True", "yes", "ON"):
        monkeypatch.setenv("DELFIN_PRE_UFF_MD_SNAP", val)
        assert sc._pre_uff_md_snap_enabled() is True


def test_helper_env_flag_falsy_values(monkeypatch):
    """0/false/no/empty leave the flag disabled."""
    for val in ("0", "false", "no", "off", ""):
        monkeypatch.setenv("DELFIN_PRE_UFF_MD_SNAP", val)
        assert sc._pre_uff_md_snap_enabled() is False


def test_snap_corrects_collapsed_ni_n():
    """Ni-N at 1.40 A gets snapped to the ideal ~2.05 A length."""
    mol, donor_idxs = _build_simple_complex("Ni", ["N"], md_dist=1.40)
    target = float(sc._get_ml_bond_length("Ni", "N"))
    n_snapped = sc._snap_md_distances_to_ideal(mol, 0)
    assert n_snapped == 1
    conf = mol.GetConformer(0)
    m = conf.GetAtomPosition(0)
    d = conf.GetAtomPosition(donor_idxs[0])
    new_dist = math.sqrt(
        (d.x - m.x) ** 2 + (d.y - m.y) ** 2 + (d.z - m.z) ** 2
    )
    assert abs(new_dist - target) < 0.01


def test_snap_corrects_stretched_ni_n():
    """Ni-N at 3.00 A gets snapped to the ideal length."""
    mol, donor_idxs = _build_simple_complex("Ni", ["N"], md_dist=3.00)
    target = float(sc._get_ml_bond_length("Ni", "N"))
    n_snapped = sc._snap_md_distances_to_ideal(mol, 0)
    assert n_snapped == 1
    conf = mol.GetConformer(0)
    m = conf.GetAtomPosition(0)
    d = conf.GetAtomPosition(donor_idxs[0])
    new_dist = math.sqrt(
        (d.x - m.x) ** 2 + (d.y - m.y) ** 2 + (d.z - m.z) ** 2
    )
    assert abs(new_dist - target) < 0.01


def test_snap_skips_close_to_target():
    """Distances within 10% tolerance are not modified."""
    target = float(sc._get_ml_bond_length("Ni", "N"))
    # Already ~5% off — within 10% tolerance
    mol, donor_idxs = _build_simple_complex("Ni", ["N"], md_dist=target * 1.05)
    n_snapped = sc._snap_md_distances_to_ideal(mol, 0)
    assert n_snapped == 0


def test_snap_handles_no_metal_gracefully():
    """A purely organic mol returns 0 (no metals)."""
    smi = "CCO"
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    n_snapped = sc._snap_md_distances_to_ideal(mol, 0)
    assert n_snapped == 0


def test_snap_translates_fragment_rigidly():
    """The donor's bonded fragment moves together (not just the donor atom)."""
    # Build Ni-N-C with N at 1.4 A from Ni and C at 1.4 A from N
    m = Chem.RWMol()
    m.AddAtom(Chem.Atom("Ni"))
    m.AddAtom(Chem.Atom("N"))
    m.AddAtom(Chem.Atom("C"))
    m.AddBond(0, 1, Chem.BondType.SINGLE)
    m.AddBond(1, 2, Chem.BondType.SINGLE)
    conf = Chem.Conformer(3)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.40, 0.0, 0.0))   # Ni-N collapsed
    conf.SetAtomPosition(2, Point3D(2.80, 0.0, 0.0))   # C bonded to N
    m.AddConformer(conf, assignId=True)
    mol = m.GetMol()

    n_snapped = sc._snap_md_distances_to_ideal(mol, 0)
    assert n_snapped == 1
    target = float(sc._get_ml_bond_length("Ni", "N"))
    conf = mol.GetConformer(0)
    n_pos = conf.GetAtomPosition(1)
    c_pos = conf.GetAtomPosition(2)
    # N now at target
    assert abs(n_pos.x - target) < 0.01
    # C-N bond length preserved (rigid translation)
    cn_dist = math.sqrt(
        (c_pos.x - n_pos.x) ** 2 + (c_pos.y - n_pos.y) ** 2 + (c_pos.z - n_pos.z) ** 2
    )
    assert abs(cn_dist - 1.40) < 0.01


def test_snap_does_not_cross_metals():
    """BFS fragment stops at metals (M-X-M bridge: only the donor moves,
    not the second metal)."""
    m = Chem.RWMol()
    m.AddAtom(Chem.Atom("Ni"))   # 0 = metal A
    m.AddAtom(Chem.Atom("N"))    # 1 = bridging N
    m.AddAtom(Chem.Atom("Ni"))   # 2 = metal B
    m.AddBond(0, 1, Chem.BondType.SINGLE)
    m.AddBond(1, 2, Chem.BondType.SINGLE)
    conf = Chem.Conformer(3)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.40, 0.0, 0.0))
    conf.SetAtomPosition(2, Point3D(2.80, 0.0, 0.0))
    m.AddConformer(conf, assignId=True)
    mol = m.GetMol()
    sc._snap_md_distances_to_ideal(mol, 0)
    conf = mol.GetConformer(0)
    # Second metal (Ni #2) must NOT have moved (BFS stopped at it)
    p2 = conf.GetAtomPosition(2)
    assert abs(p2.x - 2.80) < 0.01


def test_snap_skips_hydrogen_donors():
    """H donors are not snapped (they are not coordinating chemistry)."""
    m = Chem.RWMol()
    m.AddAtom(Chem.Atom("Ni"))
    m.AddAtom(Chem.Atom("H"))
    m.AddBond(0, 1, Chem.BondType.SINGLE)
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.40, 0.0, 0.0))
    m.AddConformer(conf, assignId=True)
    mol = m.GetMol()
    n_snapped = sc._snap_md_distances_to_ideal(mol, 0)
    assert n_snapped == 0


def test_snap_handles_multiple_donors():
    """All M-D pairs are snapped when each is broken."""
    mol, donor_idxs = _build_simple_complex(
        "Ni", ["N", "N", "O", "S"], md_dist=1.50
    )
    n_snapped = sc._snap_md_distances_to_ideal(mol, 0)
    assert n_snapped == 4
    conf = mol.GetConformer(0)
    m = conf.GetAtomPosition(0)
    for di, sym in zip(donor_idxs, ["N", "N", "O", "S"]):
        d = conf.GetAtomPosition(di)
        new_dist = math.sqrt(
            (d.x - m.x) ** 2 + (d.y - m.y) ** 2 + (d.z - m.z) ** 2
        )
        target = float(sc._get_ml_bond_length("Ni", sym))
        assert abs(new_dist - target) < 0.01, (
            f"donor {sym}: got {new_dist:.3f}, target {target:.3f}"
        )


def test_no_smiles_specific_paths():
    """Universal proof: snap function reads only element symbols, not SMILES.

    Build complexes with arbitrary M-D pairs across the periodic table and
    confirm snap depends only on element symbol & graph.
    """
    pairs = [("Pd", "N"), ("Pt", "Cl"), ("Cu", "S"), ("Fe", "P"), ("Co", "O")]
    for m_sym, d_sym in pairs:
        mol, donor_idxs = _build_simple_complex(m_sym, [d_sym], md_dist=1.50)
        n_snapped = sc._snap_md_distances_to_ideal(mol, 0)
        assert n_snapped == 1, f"snap failed for {m_sym}-{d_sym}"
        target = float(sc._get_ml_bond_length(m_sym, d_sym))
        conf = mol.GetConformer(0)
        m = conf.GetAtomPosition(0)
        d = conf.GetAtomPosition(donor_idxs[0])
        new_dist = math.sqrt(
            (d.x - m.x) ** 2 + (d.y - m.y) ** 2 + (d.z - m.z) ** 2
        )
        assert abs(new_dist - target) < 0.01
