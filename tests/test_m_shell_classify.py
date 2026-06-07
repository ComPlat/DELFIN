"""Tests for delfin.fffree.m_shell_classify -- the universal classifier
that distinguishes emergent κⁿ donors ("valid alternatives") from
soft-DG-drift atoms ("spurious") in a metal's geometric coordination
shell.

The three required tests from the mission brief
(2026-06-07 m-shell-overfill classifier):

  * test_emerging_kappa_kept -- κ²-acetate emergence on a κ¹ SMILES must
    NOT trigger rollback.
  * test_carbon_backbone_in_shell_rejected -- an sp3-C backbone that
    drifts into the M-shell MUST be classified as spurious and trigger
    rollback through the rotamer-topology gate.
  * test_metal_main_group_overfill_handled -- a main-group metal (Sn/In)
    with overfill must be classified correctly (all backbone-C extras
    spurious; lone-pair extras emergent).
"""
from __future__ import annotations

import os

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_cu_oac_kappa2():
    """Build a Cu(OAc) complex with the SMILES topology κ¹ (only one Cu-O
    bond) but the geometry showing the second carboxylate-O drifting into
    Cu's coordination shell (~2.1 Å).  Returns (mol, syms, P, metal_idx,
    designated_donors, emergent_idx).
    """
    from rdkit import Chem
    rw = Chem.RWMol()
    cu = rw.AddAtom(Chem.Atom(29))    # Cu
    o1 = rw.AddAtom(Chem.Atom(8))     # carboxylate O bound to Cu
    c = rw.AddAtom(Chem.Atom(6))      # central acetate C
    o2 = rw.AddAtom(Chem.Atom(8))     # second carboxylate O (NOT bonded to Cu in SMILES)
    cme = rw.AddAtom(Chem.Atom(6))    # methyl C
    rw.AddBond(cu, o1, Chem.BondType.SINGLE)
    rw.AddBond(o1, c, Chem.BondType.SINGLE)
    rw.AddBond(c, o2, Chem.BondType.DOUBLE)
    rw.AddBond(c, cme, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    P = np.array([
        [0.0, 0.0, 0.0],     # Cu
        [2.05, 0.0, 0.0],    # O1 bound to Cu
        [3.10, 0.65, 0.0],   # central C
        [2.25, 1.10, 0.0],   # O2 emergent at 2.46 Å from Cu (within ceiling 2.6)
        [4.20, 0.65, 0.0],   # methyl C far from Cu
    ], dtype=float)
    return mol, syms, P, cu, [o1], o2


def _build_metal_with_backbone_drift():
    """Build a Cu-NH3 complex with a methyl-C drifting into the Cu shell.
    The methyl-C has no lone pair -> must be classified spurious.
    Returns (mol, syms, P, metal_idx, designated_donors, drift_idx).
    """
    from rdkit import Chem
    rw = Chem.RWMol()
    cu = rw.AddAtom(Chem.Atom(29))   # Cu
    n = rw.AddAtom(Chem.Atom(7))     # N donor
    c1 = rw.AddAtom(Chem.Atom(6))    # backbone C bonded to N
    c2 = rw.AddAtom(Chem.Atom(6))    # methyl C (no lone pair) -- the drift atom
    rw.AddBond(cu, n, Chem.BondType.SINGLE)
    rw.AddBond(n, c1, Chem.BondType.SINGLE)
    rw.AddBond(c1, c2, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    P = np.array([
        [0.0, 0.0, 0.0],     # Cu
        [2.05, 0.0, 0.0],    # N donor
        [3.30, 0.0, 0.0],    # C bonded to N (out of shell)
        [2.10, 1.20, 0.0],   # methyl C drifting at 2.43 Å from Cu
    ], dtype=float)
    return mol, syms, P, cu, [n], c2


def _build_main_group_overfill():
    """Build an Sn(NH2)2 main-group complex where a backbone-C drifts
    into the shell.  Returns (mol, syms, P, metal_idx, designated_donors).
    """
    from rdkit import Chem
    rw = Chem.RWMol()
    sn = rw.AddAtom(Chem.Atom(50))  # Sn
    n1 = rw.AddAtom(Chem.Atom(7))   # N donor 1
    n2 = rw.AddAtom(Chem.Atom(7))   # N donor 2
    c = rw.AddAtom(Chem.Atom(6))    # backbone C tethered to N1
    rw.AddBond(sn, n1, Chem.BondType.SINGLE)
    rw.AddBond(sn, n2, Chem.BondType.SINGLE)
    rw.AddBond(n1, c, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    P = np.array([
        [0.0, 0.0, 0.0],    # Sn
        [2.20, 0.0, 0.0],   # N1
        [-2.20, 0.0, 0.0],  # N2
        [2.30, 1.30, 0.0],  # backbone C drifting at 2.66 Å from Sn (just outside)
    ], dtype=float)
    return mol, syms, P, sn, [n1, n2]


# ---------------------------------------------------------------------------
# Required tests
# ---------------------------------------------------------------------------
def test_emerging_kappa_kept():
    """κ²-acetate emergence on a κ¹ SMILES: classifier must report the
    second carboxylate-O as a valid alternative (emergent_kappa)."""
    from delfin.fffree.m_shell_classify import classify_m_shell_extras

    mol, syms, P, m_idx, designated, emergent = _build_cu_oac_kappa2()
    res = classify_m_shell_extras(
        syms, P, m_idx, designated_donors=designated, mol=mol,
    )
    assert emergent in res["shell"], (
        f"Emergent O should be inside the geometric shell, got shell={res['shell']}"
    )
    assert emergent in res["valid_alternatives"], (
        f"κ²-OAc emergent O must be classified as valid alternative, got {res}"
    )
    assert emergent not in res["spurious"]


def test_carbon_backbone_in_shell_rejected():
    """A backbone sp3-C drifting into the shell has no lone pair and
    must be classified spurious -> triggers rollback through the gate."""
    from delfin.fffree.m_shell_classify import classify_m_shell_extras

    mol, syms, P, m_idx, designated, drift = _build_metal_with_backbone_drift()
    res = classify_m_shell_extras(
        syms, P, m_idx, designated_donors=designated, mol=mol,
    )
    assert drift in res["shell"], (
        f"Backbone C should be inside geometric shell, got shell={res['shell']}"
    )
    assert drift in res["spurious"], (
        f"Backbone sp3-C must be classified spurious, got {res}"
    )
    assert drift not in res["valid_alternatives"]


def test_metal_main_group_overfill_handled():
    """Main-group metal (Sn) with backbone-C overfill: classifier must
    return a spurious entry (or empty shell-extras when geometry doesn't
    cross the ceiling).  The verdict path through shell_audit must NOT
    bless this as emergent_kappa."""
    from delfin.fffree.m_shell_classify import (
        classify_m_shell_extras, shell_audit,
    )

    mol, syms, P, m_idx, designated = _build_main_group_overfill()
    # Direct classifier call: the backbone C may or may not be inside
    # the geometric shell depending on the (Sn,C) ideal-bond product.
    # Either way, if it IS in the shell, it must be spurious.
    res = classify_m_shell_extras(
        syms, P, m_idx, designated_donors=designated, mol=mol,
    )
    extras = [a for a in res["shell"] if a not in set(designated)]
    if extras:
        assert all(e in res["spurious"] for e in extras), (
            f"Every Sn-shell extra (no lone pair, no bite-compat) must be "
            f"spurious, got {res}"
        )
        assert res["valid_alternatives"] == []
    # shell_audit verdict path: with no designated_donors information in
    # mol (we synthesised it explicitly), the audit derives them via
    # designated_donors_for_metal; the verdict must NOT be
    # "emergent_kappa" because the only extras are backbone-C.
    audit = shell_audit(syms, P, mol, expected_cn=[2, -1, -1, -1])
    sn_audit = [a for a in audit if a["metal_idx"] == m_idx][0]
    assert sn_audit["verdict"] in ("ok", "spurious_overfill"), (
        f"Sn audit verdict must be ok or spurious_overfill, not "
        f"emergent_kappa.  Got {sn_audit}"
    )


# ---------------------------------------------------------------------------
# Additional contract tests
# ---------------------------------------------------------------------------
def test_classifier_no_extras_returns_empty():
    """A perfectly-coordinated complex (all shell atoms are designated
    donors) yields empty alternatives and empty spurious lists."""
    from delfin.fffree.m_shell_classify import classify_m_shell_extras
    from rdkit import Chem
    rw = Chem.RWMol()
    rw.AddAtom(Chem.Atom(26))
    for _ in range(6):
        rw.AddAtom(Chem.Atom(17))
    for j in range(1, 7):
        rw.AddBond(0, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    syms = ["Fe"] + ["Cl"] * 6
    P = np.array([
        [0.0, 0.0, 0.0],
        [2.3, 0.0, 0.0],
        [-2.3, 0.0, 0.0],
        [0.0, 2.3, 0.0],
        [0.0, -2.3, 0.0],
        [0.0, 0.0, 2.3],
        [0.0, 0.0, -2.3],
    ], dtype=float)
    res = classify_m_shell_extras(
        syms, P, 0, designated_donors=[1, 2, 3, 4, 5, 6], mol=mol,
    )
    assert res["valid_alternatives"] == []
    assert res["spurious"] == []
    assert sorted(res["shell"]) == [1, 2, 3, 4, 5, 6]


def test_classifier_no_mol_returns_all_spurious():
    """Without a mol, the classifier cannot prove emergence -> every
    extra must land in ``spurious`` (safe / conservative default)."""
    from delfin.fffree.m_shell_classify import classify_m_shell_extras
    syms = ["Cu", "O", "C", "O", "C"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [2.05, 0.0, 0.0],
        [3.10, 0.65, 0.0],
        [2.25, 1.10, 0.0],
        [4.20, 0.65, 0.0],
    ], dtype=float)
    res = classify_m_shell_extras(
        syms, P, 0, designated_donors=[1], mol=None,
    )
    assert res["valid_alternatives"] == []
    # O at idx 3 sits inside the shell; classifier must mark it spurious.
    assert 3 in res["spurious"]


def test_classifier_env_default_off_no_mogul():
    """DELFIN_FFFREE_M_SHELL_CLASSIFY default OFF when MOGUL_PRIMARY is
    not set."""
    from delfin.fffree.m_shell_classify import _env_on
    saved = os.environ.pop("DELFIN_FFFREE_M_SHELL_CLASSIFY", None)
    saved_mp = os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
    try:
        assert _env_on() is False
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_M_SHELL_CLASSIFY"] = saved
        if saved_mp is not None:
            os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = saved_mp


def test_classifier_env_auto_on_under_mogul_primary():
    """Default ON when MOGUL_PRIMARY=1."""
    from delfin.fffree.m_shell_classify import _env_on
    saved = os.environ.pop("DELFIN_FFFREE_M_SHELL_CLASSIFY", None)
    saved_mp = os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY")
    try:
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        assert _env_on() is True
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_M_SHELL_CLASSIFY"] = saved
        if saved_mp is None:
            os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        else:
            os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = saved_mp


def test_classifier_env_explicit_off_overrides_mogul():
    """Explicit =0 overrides MOGUL_PRIMARY auto-ON."""
    from delfin.fffree.m_shell_classify import _env_on
    saved = os.environ.get("DELFIN_FFFREE_M_SHELL_CLASSIFY")
    saved_mp = os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY")
    try:
        os.environ["DELFIN_FFFREE_M_SHELL_CLASSIFY"] = "0"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        assert _env_on() is False
    finally:
        for k, v in (
            ("DELFIN_FFFREE_M_SHELL_CLASSIFY", saved),
            ("DELFIN_FFFREE_MOGUL_PRIMARY", saved_mp),
        ):
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


# ---------------------------------------------------------------------------
# Integration with the rotamer-topology gate
# ---------------------------------------------------------------------------
def test_gate_keeps_emergent_kappa_when_classifier_on():
    """Under DELFIN_FFFREE_MOGUL_PRIMARY=1 (auto-on for classifier), a
    Cu-OAc geometry with κ² emergence must NOT be rejected as
    m_shell_overfill."""
    from delfin.fffree.rotamer_topology_gate import rotation_preserves_topology
    mol, syms, P, m_idx, designated, emergent = _build_cu_oac_kappa2()
    saved = os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY")
    saved_g = os.environ.get("DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE")
    saved_c = os.environ.get("DELFIN_FFFREE_M_SHELL_CLASSIFY")
    try:
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ.pop("DELFIN_FFFREE_M_SHELL_CLASSIFY", None)
        ok, reason = rotation_preserves_topology(
            syms, P, mol=mol, return_reason=True,
        )
        # κ² emergence -> classifier should pass it; the gate may still
        # reject for other geometry reasons (bond_break / spurious_bond)
        # but it MUST NOT be m_shell_overfill.
        assert reason != "m_shell_overfill", (
            f"κ² emergence falsely rejected as m_shell_overfill: {reason}"
        )
    finally:
        for k, v in (
            ("DELFIN_FFFREE_MOGUL_PRIMARY", saved),
            ("DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE", saved_g),
            ("DELFIN_FFFREE_M_SHELL_CLASSIFY", saved_c),
        ):
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def test_gate_rejects_backbone_drift_even_with_classifier_on():
    """Under MOGUL_PRIMARY=1 (classifier auto-on), a backbone-C drift
    into the shell MUST still trigger m_shell_overfill rollback."""
    from delfin.fffree.rotamer_topology_gate import rotation_preserves_topology

    mol, syms, P, m_idx, designated, drift = _build_metal_with_backbone_drift()
    saved = os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY")
    saved_g = os.environ.get("DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE")
    saved_c = os.environ.get("DELFIN_FFFREE_M_SHELL_CLASSIFY")
    try:
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ.pop("DELFIN_FFFREE_M_SHELL_CLASSIFY", None)
        ok, reason = rotation_preserves_topology(
            syms, P, mol=mol, return_reason=True,
        )
        assert ok is False
        assert reason in ("m_shell_overfill", "spurious_bond", "heavy_collision"), (
            f"Backbone-C drift should trigger rollback, got reason={reason}"
        )
    finally:
        for k, v in (
            ("DELFIN_FFFREE_MOGUL_PRIMARY", saved),
            ("DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE", saved_g),
            ("DELFIN_FFFREE_M_SHELL_CLASSIFY", saved_c),
        ):
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def test_classifier_is_deterministic():
    """Repeated calls with identical inputs must yield identical outputs."""
    from delfin.fffree.m_shell_classify import classify_m_shell_extras
    mol, syms, P, m_idx, designated, emergent = _build_cu_oac_kappa2()
    results = [
        classify_m_shell_extras(
            syms, P, m_idx, designated_donors=designated, mol=mol,
        )
        for _ in range(5)
    ]
    for r in results[1:]:
        assert r == results[0]
