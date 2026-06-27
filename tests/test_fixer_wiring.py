"""Unit tests for fixer dispatch wiring (F19, F25, WUXQAK) in smiles_converter.

Wiring contract (per fixer):

  * Default OFF: env-flag unset -> input results unchanged byte-for-byte.
  * Env-flag = 1: dispatcher invokes the underlying fixer; XYZ output may
    change when a real violation is present in the input frame.
  * dual_parse_done=True bypasses the dispatcher (avoids double-application
    on the inner dual-parse call).
  * mol=None safe (returns input unchanged).

The tests build minimal synthetic XYZ strings + RDKit Mol objects and call
the three ``_apply_fixer_<X>_if_enabled`` helpers directly, asserting the
default-OFF contract is bit-exact and the env-on path actually transforms
geometries when a violation exists.
"""
from __future__ import annotations

import os
from typing import List, Tuple

import pytest

# Graceful imports.
Chem = pytest.importorskip("rdkit.Chem", reason="RDKit required for tests")
AllChem = pytest.importorskip("rdkit.Chem.AllChem", reason="RDKit required")

import numpy as np  # noqa: E402

sc = pytest.importorskip(
    "delfin.smiles_converter",
    reason="delfin.smiles_converter required",
)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _mol_to_xyz(mol) -> str:
    out: List[str] = [str(mol.GetNumAtoms()), ""]
    conf = mol.GetConformer()
    for i, a in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        out.append(f"{a.GetSymbol()} {p.x:.6f} {p.y:.6f} {p.z:.6f}")
    return "\n".join(out) + "\n"


def _build_methylamine() -> Tuple[object, str]:
    """Return (mol, xyz) for embedded + UFF-relaxed CH3-NH2."""
    m = Chem.AddHs(Chem.MolFromSmiles("CN"))
    AllChem.EmbedMolecule(m, randomSeed=42)
    AllChem.UFFOptimizeMolecule(m)
    return m, _mol_to_xyz(m)


def _flatten_amine_nitrogen(mol, xyz: str) -> str:
    """Manually flatten the NH2 in methylamine to trigger F25."""
    from delfin.manta._fix_sp3_n_pyramidality import _parse_xyz, _format_xyz
    syms, pts, orig = _parse_xyz(xyz)
    n_idx = syms.index("N")
    n_pos = pts[n_idx].copy()
    for i, s in enumerate(syms):
        if s == "H" and np.linalg.norm(pts[i] - n_pos) < 1.2:
            v = pts[i] - n_pos
            v[2] = 0.0
            n = np.linalg.norm(v)
            if n > 1e-6:
                v *= 1.01 / n
                pts[i] = n_pos + v
    return _format_xyz(orig, syms, pts)


def _build_pt_methyl_with_linear_c() -> Tuple[object, str]:
    """Pt(CH3)Cl3 with one Pt-C-H deliberately forced to 180 degrees."""
    m = Chem.AddHs(Chem.MolFromSmiles("[Pt](C)(Cl)(Cl)Cl"))
    AllChem.EmbedMolecule(m, randomSeed=42)
    xyz = _mol_to_xyz(m)
    from delfin.manta._fix_wuxqak_sp3_c_linear import _parse_xyz, _format_xyz
    syms, pts, orig = _parse_xyz(xyz)
    pt_idx = syms.index("Pt")
    c_idx = syms.index("C")
    # Pick first H bonded to C
    for i, s in enumerate(syms):
        if s == "H" and np.linalg.norm(pts[i] - pts[c_idx]) < 1.2:
            v = pts[c_idx] - pts[pt_idx]
            v /= np.linalg.norm(v)
            bond = np.linalg.norm(pts[i] - pts[c_idx])
            pts[i] = pts[c_idx] + v * bond
            break
    return m, _format_xyz(orig, syms, pts)


def _scrub_fixer_env(monkeypatch):
    """Remove any pre-existing fixer env-vars so default-OFF holds."""
    for k in (
        "DELFIN_FIX_F19", "DELFIN_FIX_F19_TOL_DEG", "DELFIN_FIX_F19_CLASSES",
        "DELFIN_FIX_F25", "DELFIN_FIX_F25_THRESHOLD_DEG",
        "DELFIN_FIX_F25_TARGET_DEG", "DELFIN_FIX_F25_CLASSES",
        "DELFIN_FIX_WUXQAK", "DELFIN_FIX_WUXQAK_ANGLE_DEG",
        "DELFIN_FIX_WUXQAK_TARGET_DEG", "DELFIN_FIX_WUXQAK_CLASSES",
    ):
        monkeypatch.delenv(k, raising=False)


# --------------------------------------------------------------------------
# F19 — sp3-H tetrahedrality fixer
# --------------------------------------------------------------------------

def test_f19_default_off_bit_exact(monkeypatch):
    """DELFIN_FIX_F19 unset -> dispatcher leaves results identical."""
    _scrub_fixer_env(monkeypatch)
    mol, xyz = _build_methylamine()
    results = [(xyz, "test")]
    out = sc._apply_fixer_f19_if_enabled(mol, results, dual_parse_done=False)
    assert out == results, "F19 default-OFF must be bit-exact"


def test_f19_active_changes_xyz(monkeypatch):
    """DELFIN_FIX_F19=1 with tight tolerance -> XYZ actually changes."""
    _scrub_fixer_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_F19", "1")
    monkeypatch.setenv("DELFIN_FIX_F19_TOL_DEG", "1.0")  # force violations
    mol, xyz = _build_methylamine()
    results = [(xyz, "test")]
    out = sc._apply_fixer_f19_if_enabled(mol, results, dual_parse_done=False)
    assert len(out) == 1
    new_xyz, new_lbl = out[0]
    assert new_lbl == "test"
    assert new_xyz != xyz, "F19 must transform XYZ at tight tolerance"


def test_f19_dual_parse_bypass(monkeypatch):
    """dual_parse_done=True bypasses the dispatcher even when env=1."""
    _scrub_fixer_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_F19", "1")
    monkeypatch.setenv("DELFIN_FIX_F19_TOL_DEG", "1.0")
    mol, xyz = _build_methylamine()
    results = [(xyz, "test")]
    out = sc._apply_fixer_f19_if_enabled(mol, results, dual_parse_done=True)
    assert out == results, "dual_parse_done must short-circuit dispatcher"


# --------------------------------------------------------------------------
# F25 — sp3-N pyramidality fixer
# --------------------------------------------------------------------------

def test_f25_default_off_bit_exact(monkeypatch):
    _scrub_fixer_env(monkeypatch)
    mol, xyz = _build_methylamine()
    flat_xyz = _flatten_amine_nitrogen(mol, xyz)
    results = [(flat_xyz, "test")]
    out = sc._apply_fixer_f25_if_enabled(mol, results, dual_parse_done=False)
    assert out == results, "F25 default-OFF must be bit-exact"


def test_f25_active_changes_xyz(monkeypatch):
    _scrub_fixer_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_F25", "1")
    monkeypatch.setenv("DELFIN_FIX_F25_THRESHOLD_DEG", "340.0")
    mol, xyz = _build_methylamine()
    flat_xyz = _flatten_amine_nitrogen(mol, xyz)
    results = [(flat_xyz, "test")]
    out = sc._apply_fixer_f25_if_enabled(mol, results, dual_parse_done=False)
    assert len(out) == 1
    new_xyz, _ = out[0]
    assert new_xyz != flat_xyz, "F25 must repyramidalize flat sp3-N"


def test_f25_mol_none_safe(monkeypatch):
    _scrub_fixer_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_F25", "1")
    mol, xyz = _build_methylamine()
    results = [(xyz, "test")]
    out = sc._apply_fixer_f25_if_enabled(None, results, dual_parse_done=False)
    assert out == results, "mol=None must return input unchanged"


# --------------------------------------------------------------------------
# WUXQAK — sp3-C linear-collapse fixer
# --------------------------------------------------------------------------

def test_wuxqak_default_off_bit_exact(monkeypatch):
    _scrub_fixer_env(monkeypatch)
    mol, xyz = _build_pt_methyl_with_linear_c()
    results = [(xyz, "test")]
    out = sc._apply_fixer_wuxqak_if_enabled(mol, results, dual_parse_done=False)
    assert out == results, "WUXQAK default-OFF must be bit-exact"


def test_wuxqak_active_changes_xyz(monkeypatch):
    _scrub_fixer_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_WUXQAK", "1")
    mol, xyz = _build_pt_methyl_with_linear_c()
    results = [(xyz, "test")]
    out = sc._apply_fixer_wuxqak_if_enabled(mol, results, dual_parse_done=False)
    assert len(out) == 1
    new_xyz, _ = out[0]
    assert new_xyz != xyz, "WUXQAK must bend a 180-degree M-C-H back toward 109.5"


def test_wuxqak_no_metal_short_circuit(monkeypatch):
    """No metal in any XYZ -> dispatcher exits via fast-path unchanged."""
    _scrub_fixer_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_WUXQAK", "1")
    # Pure organic methylamine — no metal symbol
    mol, xyz = _build_methylamine()
    results = [(xyz, "test")]
    out = sc._apply_fixer_wuxqak_if_enabled(mol, results, dual_parse_done=False)
    assert out == results, "no-metal pre-check must short-circuit"


# --------------------------------------------------------------------------
# SP2N-PLANARIZE — sp2-N / nitro planarisation fixer
# --------------------------------------------------------------------------

def _build_broken_nitromethane() -> Tuple[object, str]:
    """Return (mol, xyz) for CH3-NO2 with a deliberately broken nitro:
    N pushed out of the C-O-O plane, N-O desymmetrised (1.43 / 1.12 A) --
    the AVIDAM-class defect."""
    m = Chem.AddHs(Chem.MolFromSmiles("CN(=O)=O"))
    AllChem.EmbedMolecule(m, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(m)
    from delfin.manta._fix_sp2n_planarize import (_parse_xyz, _format_xyz,
                                            detect_nitro_groups)
    xyz = _mol_to_xyz(m)
    syms, pts, orig = _parse_xyz(xyz)
    groups = detect_nitro_groups(m)
    assert groups, "test fixture must contain a detectable nitro group"
    g = groups[0]
    ni, c = g["n_idx"], g["c_idx"]
    o1, o2 = g["o_idxs"]
    # Normal of the original (good) C-O-O plane.
    nrm = np.cross(pts[o1] - pts[c], pts[o2] - pts[c])
    nrm = nrm / np.linalg.norm(nrm)
    # Push N a full 1.0 A out of plane along the plane normal (keep O fixed),
    # then desymmetrise the two N-O bonds about the displaced N.
    pts[ni] = pts[ni] + 1.0 * nrm
    pts[o1] = pts[ni] + (pts[o1] - pts[ni]) / np.linalg.norm(
        pts[o1] - pts[ni]) * 1.43                  # stretch one N-O
    pts[o2] = pts[ni] + (pts[o2] - pts[ni]) / np.linalg.norm(
        pts[o2] - pts[ni]) * 1.12                  # shorten the other
    return m, _format_xyz(orig, syms, pts)


def _nitro_metrics(mol, xyz):
    from delfin.manta._fix_sp2n_planarize import (_parse_xyz, detect_nitro_groups,
                                            _angle_deg, _oop_distance)
    syms, P, _ = _parse_xyz(xyz)
    g = detect_nitro_groups(mol)[0]
    ni, c = g["n_idx"], g["c_idx"]
    o1, o2 = g["o_idxs"]
    return (
        float(np.linalg.norm(P[ni] - P[o1])),
        float(np.linalg.norm(P[ni] - P[o2])),
        _angle_deg(P[ni], P[o1], P[o2]),
        _oop_distance(P[ni], P[o1], P[o2], P[c]),
    )


def test_sp2n_default_off_bit_exact(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_SP2N_PLANARIZE", raising=False)
    mol, xyz = _build_broken_nitromethane()
    results = [(xyz, "test")]
    out = sc._apply_fixer_sp2n_planarize_if_enabled(
        mol, results, dual_parse_done=False)
    assert out == results, "SP2N default-OFF must be bit-exact"


def test_sp2n_active_planarizes_nitro(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_SP2N_PLANARIZE", "1")
    mol, xyz = _build_broken_nitromethane()
    no1_b, no2_b, _ono_b, oop_b = _nitro_metrics(mol, xyz)
    assert oop_b > 0.3 and abs(no1_b - no2_b) > 0.15  # fixture is broken
    results = [(xyz, "test")]
    out = sc._apply_fixer_sp2n_planarize_if_enabled(
        mol, results, dual_parse_done=False)
    new_xyz, lbl = out[0]
    assert lbl == "test"
    assert new_xyz != xyz, "SP2N must transform a broken nitro"
    no1_a, no2_a, ono_a, oop_a = _nitro_metrics(mol, new_xyz)
    assert oop_a < 0.05, f"N must become planar (oop={oop_a})"
    assert abs(no1_a - no2_a) < 0.02, "N-O must symmetrise"
    assert abs(no1_a - 1.22) < 0.03 and abs(no2_a - 1.22) < 0.03
    assert abs(ono_a - 125.0) < 2.0, f"O-N-O must reach ~125 (got {ono_a})"


def test_sp2n_dual_parse_bypass(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_SP2N_PLANARIZE", "1")
    mol, xyz = _build_broken_nitromethane()
    results = [(xyz, "test")]
    out = sc._apply_fixer_sp2n_planarize_if_enabled(
        mol, results, dual_parse_done=True)
    assert out == results, "dual_parse_done must short-circuit dispatcher"


def test_sp2n_mol_none_safe(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_SP2N_PLANARIZE", "1")
    _mol, xyz = _build_broken_nitromethane()
    results = [(xyz, "test")]
    out = sc._apply_fixer_sp2n_planarize_if_enabled(
        None, results, dual_parse_done=False)
    assert out == results, "mol=None must return input unchanged"


# --------------------------------------------------------------------------
# Pipeline wiring — call-site smoke test
# --------------------------------------------------------------------------

def test_dispatchers_callable_from_module():
    """All four dispatchers exist and are call-compatible."""
    for name in (
        "_apply_fixer_f19_if_enabled",
        "_apply_fixer_f25_if_enabled",
        "_apply_fixer_sp2n_planarize_if_enabled",
        "_apply_fixer_wuxqak_if_enabled",
    ):
        fn = getattr(sc, name, None)
        assert callable(fn), f"{name} missing from smiles_converter"
        # Empty results -> empty results (no-op path)
        assert fn(None, [], False) == []
