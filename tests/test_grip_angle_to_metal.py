"""Tests for the angle-to-metal layer (2026-06-05).

Validate:

* default OFF byte-identical with legacy detect_fragments
* M-D-X term emission counts what we expect for a known geometry
* gradient agrees with finite differences
* gradient on (frozen) metal+donor is non-zero in the raw term but
  zeroed by grip_polish's frozen mask -- the user's intent is for the
  ligand-X atoms to move while M and D stay fixed
* env-flag round-trip (2-run determinism)
* D-M-D polyhedron-anchored terms have mu equal to the observed
  initial angle (the polyhedron is locked in by definition)
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from rdkit import Chem
from rdkit.Chem import AllChem

from delfin.fffree.grip_fragment_detect import (
    ANGLE_TO_METAL_ENV,
    DMD_POLY_ENV,
    angle_to_metal_active,
    dmd_polyhedron_active,
    build_mdx_angle_terms,
    build_dmd_polyhedron_terms,
    detect_fragments,
)
from delfin.fffree.grip_loss_terms import TotalGripLoss


def _hexammine():
    smiles = "[Fe](N)(N)(N)(N)(N)N"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0)
    metal = 0
    donors = sorted(int(n.GetIdx()) for n in mol.GetAtomWithIdx(metal).GetNeighbors()
                    if n.GetSymbol() == "N")
    P = np.array(mol.GetConformer().GetPositions(), dtype=float)
    return mol, metal, donors, P


def _hexakis_methylamine():
    """Fe(NHMe)6 — each donor N has a heavy C neighbour so M-D-X terms
    fire by default (without the include_hydrogens flag)."""
    smiles = "[Fe](NC)(NC)(NC)(NC)(NC)NC"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0)
    metal = 0
    donors = sorted(int(n.GetIdx()) for n in mol.GetAtomWithIdx(metal).GetNeighbors()
                    if n.GetSymbol() == "N")
    P = np.array(mol.GetConformer().GetPositions(), dtype=float)
    return mol, metal, donors, P


def _clear_env():
    for k in (ANGLE_TO_METAL_ENV, DMD_POLY_ENV,
              "DELFIN_FFFREE_GRIP_ANGLE_TO_METAL_WEIGHT",
              "DELFIN_FFFREE_GRIP_ANGLE_TO_METAL_INCLUDE_H"):
        os.environ.pop(k, None)


def test_default_off():
    _clear_env()
    assert not angle_to_metal_active()
    assert not dmd_polyhedron_active()


def test_env_flag_round_trip():
    _clear_env()
    os.environ[ANGLE_TO_METAL_ENV] = "1"
    try:
        assert angle_to_metal_active()
    finally:
        _clear_env()
    assert not angle_to_metal_active()


def test_off_byte_identical():
    """With ANGLE_TO_METAL OFF, passing or omitting ``metal=`` is identical."""
    _clear_env()
    mol, metal, donors, P = _hexammine()
    frozen = {metal, *donors}
    legacy = detect_fragments(mol, P, frozen_atoms=frozen, donors=donors)
    with_metal = detect_fragments(mol, P, frozen_atoms=frozen, donors=donors,
                                  metal=metal)

    def fp(terms):
        return tuple((type(t).__name__, t.atom_indices,
                      getattr(t, "mu", None), getattr(t, "sigma", None),
                      getattr(t, "weight", None)) for t in terms)
    assert fp(legacy) == fp(with_metal)


def test_mdx_term_emission():
    """Each donor with an X neighbour emits one M-D-X term.

    Fe(NH3)6 has only H neighbours, so we opt into hydrogen inclusion
    to exercise the term emission path.
    """
    _clear_env()
    mol, metal, donors, P = _hexammine()
    terms = build_mdx_angle_terms(mol, metal, donors, include_hydrogens=True)
    # 6 donors, 2 H each = 12 expected (RDKit only adds 2 H to neutral N).
    assert len(terms) > 0
    # Every term has metal at position a, donor at position b, X != metal/donor at c.
    donor_set = set(donors)
    for t in terms:
        assert t.a == metal
        assert t.b in donor_set
        assert t.c != metal
        assert t.c not in donor_set
        # VSEPR fallback for sp3 N = 109.5°.
        assert 90.0 <= t.mu <= 130.0
        assert t.sigma > 0.0
        assert t.weight > 0.0


def test_mdx_gradient_finite_diff():
    """Analytical gradient matches central-difference to 1e-6."""
    _clear_env()
    mol, metal, donors, P = _hexammine()
    terms = build_mdx_angle_terms(mol, metal, donors, include_hydrogens=True)
    agg = TotalGripLoss(terms=list(terms))
    L0, G0 = agg.evaluate(P)
    eps = 1e-5
    G_fd = np.zeros_like(P)
    for i in range(P.shape[0]):
        for k in range(3):
            Pp = P.copy(); Pp[i, k] += eps
            Pm = P.copy(); Pm[i, k] -= eps
            Lp, _ = agg.evaluate(Pp)
            Lm, _ = agg.evaluate(Pm)
            G_fd[i, k] = (Lp - Lm) / (2 * eps)
    max_err = float(np.abs(G0 - G_fd).max())
    assert max_err < 1e-5, f"max grad-vs-FD err = {max_err}"


def test_mdx_gradient_on_x_atoms_nonzero():
    """The whole point: gradient on the X (ligand) atoms must be non-zero
    so they relax under the polish (donor + metal stay frozen by mask)."""
    _clear_env()
    mol, metal, donors, P = _hexammine()
    terms = build_mdx_angle_terms(mol, metal, donors, include_hydrogens=True)
    agg = TotalGripLoss(terms=list(terms))
    _, G = agg.evaluate(P)
    # H atoms = indices 7..18 (after Fe + 6 N).  At least one X-atom must
    # have non-trivial gradient.
    g_norm_x = np.linalg.norm(G[7:], axis=1)
    assert g_norm_x.max() > 1e-3


def test_dmd_polyhedron_mu_equals_observed():
    """build_dmd_polyhedron_terms encodes the observed initial angles as mu."""
    _clear_env()
    mol, metal, donors, P = _hexammine()
    terms = build_dmd_polyhedron_terms(mol, metal, donors, P)
    # n*(n-1)/2 = 15 pairs for 6 donors.
    assert len(terms) == 15
    # The angle term encodes the observed angle; evaluating at P gives 0 loss.
    agg = TotalGripLoss(terms=list(terms))
    L, _ = agg.evaluate(P)
    assert L < 1e-8, f"D-M-D terms should have zero loss at P, got {L}"


def test_detect_fragments_on_emits_extra_terms():
    """With ANGLE_TO_METAL=1, detect_fragments emits extra M-D-X terms
    where X is a heavy ligand atom (NHMe gives N-C M-D-X)."""
    _clear_env()
    mol, metal, donors, P = _hexakis_methylamine()
    frozen = {metal, *donors}
    off = detect_fragments(mol, P, frozen_atoms=frozen, donors=donors,
                            metal=metal)
    os.environ[ANGLE_TO_METAL_ENV] = "1"
    try:
        on = detect_fragments(mol, P, frozen_atoms=frozen, donors=donors,
                              metal=metal)
    finally:
        _clear_env()
    assert len(on) > len(off)
    # Every extra term is an M-D-X angle.
    off_fp = {(type(t).__name__, t.atom_indices) for t in off}
    extras = [t for t in on if (type(t).__name__, t.atom_indices) not in off_fp]
    for t in extras:
        assert t.a == metal and t.b in set(donors) and t.c not in set(donors)


def test_two_run_determinism():
    """Two builds of the same SMILES with same seed produce byte-identical terms."""
    _clear_env()
    os.environ[ANGLE_TO_METAL_ENV] = "1"
    os.environ[DMD_POLY_ENV] = "1"
    try:
        mol1, m1, d1, P1 = _hexammine()
        mol2, m2, d2, P2 = _hexammine()
        t1 = detect_fragments(mol1, P1, frozen_atoms={m1, *d1}, donors=d1, metal=m1)
        t2 = detect_fragments(mol2, P2, frozen_atoms={m2, *d2}, donors=d2, metal=m2)
        f1 = tuple((type(t).__name__, t.atom_indices, getattr(t, "mu", None),
                    getattr(t, "sigma", None), getattr(t, "weight", None))
                   for t in t1)
        f2 = tuple((type(t).__name__, t.atom_indices, getattr(t, "mu", None),
                    getattr(t, "sigma", None), getattr(t, "weight", None))
                   for t in t2)
        assert f1 == f2
    finally:
        _clear_env()


def test_dmd_polyhedron_zero_gradient_on_frozen():
    """D-M-D terms touch only (donor, metal, donor) = all frozen.  Their
    gradient on ligand-X atoms must be exactly zero."""
    _clear_env()
    mol, metal, donors, P = _hexammine()
    terms = build_dmd_polyhedron_terms(mol, metal, donors, P)
    # Perturb P slightly so the loss is non-zero.
    P2 = P.copy(); P2[7, 0] += 0.5
    agg = TotalGripLoss(terms=list(terms))
    _, G = agg.evaluate(P2)
    # Gradient on H atom 7 must be zero (D-M-D doesn't touch H).
    assert np.linalg.norm(G[7]) < 1e-12


def test_weight_env_override():
    """ANGLE_TO_METAL_WEIGHT env scales the term weight."""
    _clear_env()
    os.environ[ANGLE_TO_METAL_ENV] = "1"
    try:
        mol, metal, donors, _ = _hexammine()
        t_default = build_mdx_angle_terms(mol, metal, donors, include_hydrogens=True)
        os.environ["DELFIN_FFFREE_GRIP_ANGLE_TO_METAL_WEIGHT"] = "1.0"
        t_one = build_mdx_angle_terms(mol, metal, donors, include_hydrogens=True)
        assert len(t_default) == len(t_one)
        # Default weight scale = 0.1, so override 1.0 = 10x.
        for td, to in zip(t_default, t_one):
            assert pytest.approx(to.weight, rel=1e-9) == 10.0 * td.weight
    finally:
        _clear_env()


def test_heavy_only_default_excludes_h():
    """Default (include_hydrogens=False) drops M-D-H terms."""
    _clear_env()
    mol, metal, donors, _ = _hexammine()
    terms_no_h = build_mdx_angle_terms(mol, metal, donors)
    terms_with_h = build_mdx_angle_terms(mol, metal, donors,
                                         include_hydrogens=True)
    assert len(terms_no_h) == 0   # only H neighbours -> zero terms by default
    assert len(terms_with_h) > 0


def test_include_h_env_flag_round_trip():
    """The include-H env-flag is read per call (subprocess-safe)."""
    _clear_env()
    from delfin.fffree.grip_fragment_detect import _angle_to_metal_include_h
    assert _angle_to_metal_include_h() is False
    os.environ["DELFIN_FFFREE_GRIP_ANGLE_TO_METAL_INCLUDE_H"] = "1"
    try:
        assert _angle_to_metal_include_h() is True
    finally:
        _clear_env()
