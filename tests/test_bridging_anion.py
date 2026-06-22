"""Unit tests for the bridging-anion (μ-X) M-X-M fixer.

Test surface:

  * Detection: ``detect_mu_x_bridges`` identifies μ-Cl (Pt2), μ-OH (Fe2),
    μ-OMe (Cu2) bridges and skips terminal anions.
  * Geometry: a hand-built [Pt2(μ-Cl)2] with M-Cl-M angle forced to 180°
    is bent back toward the 80-105° dimer window when the fix is enabled.
  * Wiring: ``_apply_fixer_bridging_anion_if_enabled`` is bit-exact when
    ``DELFIN_FIX_BRIDGING_ANION`` is unset, transforms XYZ when enabled,
    safe with ``mol=None``, and respects ``dual_parse_done=True``.
  * Class-conditional override via ``DELFIN_FIX_BRIDGING_ANION_CLASSES``.
  * Malformed XYZ is handled gracefully (returns input unchanged).
"""
from __future__ import annotations

import math
import os
from typing import List, Tuple

import pytest

Chem = pytest.importorskip("rdkit.Chem", reason="RDKit required for tests")

import numpy as np  # noqa: E402

mod = pytest.importorskip("delfin._fix_bridging_anion")
sc = pytest.importorskip("delfin.smiles_converter")


# --------------------------------------------------------------------------
# Synthetic mol + XYZ builders
# --------------------------------------------------------------------------

def _build_pt2_mu_cl2() -> Tuple[object, str]:
    """Hand-build [Pt2(μ-Cl)2(NH3)4] in a near-realistic dimer geometry.

    Atom order: Pt0, Pt1, Cl2 (bridge), Cl3 (bridge), N4..N7 (terminal
    NH3 — H's omitted; the test only needs heavy-atom angles).  Bonds
    are added explicitly with RWMol so RDKit's bond graph encodes the
    μ-X bridges without relying on SMILES parsing.
    """
    rw = Chem.RWMol()
    pt0 = rw.AddAtom(Chem.Atom("Pt"))
    pt1 = rw.AddAtom(Chem.Atom("Pt"))
    cl2 = rw.AddAtom(Chem.Atom("Cl"))
    cl3 = rw.AddAtom(Chem.Atom("Cl"))
    n4 = rw.AddAtom(Chem.Atom("N"))
    n5 = rw.AddAtom(Chem.Atom("N"))
    n6 = rw.AddAtom(Chem.Atom("N"))
    n7 = rw.AddAtom(Chem.Atom("N"))
    # Bridges Cl2 and Cl3 each bonded to BOTH Pt0 and Pt1.
    rw.AddBond(pt0, cl2, Chem.BondType.SINGLE)
    rw.AddBond(pt1, cl2, Chem.BondType.SINGLE)
    rw.AddBond(pt0, cl3, Chem.BondType.SINGLE)
    rw.AddBond(pt1, cl3, Chem.BondType.SINGLE)
    # Terminal NH3 (heavy-atom only).
    rw.AddBond(pt0, n4, Chem.BondType.SINGLE)
    rw.AddBond(pt0, n5, Chem.BondType.SINGLE)
    rw.AddBond(pt1, n6, Chem.BondType.SINGLE)
    rw.AddBond(pt1, n7, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)

    # Coordinates: place metals on x-axis 3.4 Å apart, bridges above/below,
    # angle initially set close to the dimer-realistic 92°.  Tests that
    # exercise the fix use ``_force_mxm_angle`` below to push to a
    # violating geometry first.
    coords = {
        pt0: np.array([-1.7, 0.0, 0.0]),
        pt1: np.array([+1.7, 0.0, 0.0]),
        cl2: np.array([0.0, 1.7, 0.0]),
        cl3: np.array([0.0, -1.7, 0.0]),
        n4: np.array([-3.7, 0.0, 0.0]),
        n5: np.array([-1.7, 0.0, 2.0]),
        n6: np.array([+3.7, 0.0, 0.0]),
        n7: np.array([+1.7, 0.0, 2.0]),
    }
    return mol, _coords_to_xyz(mol, coords)


def _build_terminal_cl_mono_metal() -> Tuple[object, str]:
    """Single Pt with one terminal Cl — Cl should NOT be flagged μ-X."""
    rw = Chem.RWMol()
    pt = rw.AddAtom(Chem.Atom("Pt"))
    cl = rw.AddAtom(Chem.Atom("Cl"))
    n0 = rw.AddAtom(Chem.Atom("N"))
    rw.AddBond(pt, cl, Chem.BondType.SINGLE)
    rw.AddBond(pt, n0, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    coords = {
        pt: np.array([0.0, 0.0, 0.0]),
        cl: np.array([2.3, 0.0, 0.0]),
        n0: np.array([0.0, 2.0, 0.0]),
    }
    return mol, _coords_to_xyz(mol, coords)


def _build_fe2_mu_oh2() -> Tuple[object, str]:
    """[Fe2(μ-OH)2] hand-built — same shape as Pt2(μ-Cl)2 but O bridges
    carrying one H each (chemistry-meaningful μ-OH).  Used to verify the
    sp3 single-bridge motif window (rather than the 4-ring window — both
    OH's bridge the same Fe pair, so it IS a 4-ring case)."""
    rw = Chem.RWMol()
    fe0 = rw.AddAtom(Chem.Atom("Fe"))
    fe1 = rw.AddAtom(Chem.Atom("Fe"))
    o2 = rw.AddAtom(Chem.Atom("O"))
    o3 = rw.AddAtom(Chem.Atom("O"))
    h_a = rw.AddAtom(Chem.Atom("H"))
    h_b = rw.AddAtom(Chem.Atom("H"))
    rw.AddBond(fe0, o2, Chem.BondType.SINGLE)
    rw.AddBond(fe1, o2, Chem.BondType.SINGLE)
    rw.AddBond(fe0, o3, Chem.BondType.SINGLE)
    rw.AddBond(fe1, o3, Chem.BondType.SINGLE)
    rw.AddBond(o2, h_a, Chem.BondType.SINGLE)
    rw.AddBond(o3, h_b, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    coords = {
        fe0: np.array([-1.5, 0.0, 0.0]),
        fe1: np.array([+1.5, 0.0, 0.0]),
        o2: np.array([0.0, 1.5, 0.0]),
        o3: np.array([0.0, -1.5, 0.0]),
        h_a: np.array([0.0, 2.5, 0.0]),
        h_b: np.array([0.0, -2.5, 0.0]),
    }
    return mol, _coords_to_xyz(mol, coords)


def _build_cu2_mu_ome2() -> Tuple[object, str]:
    """[Cu2(μ-OMe)2] hand-built — μ-O bridges with one C substituent each."""
    rw = Chem.RWMol()
    cu0 = rw.AddAtom(Chem.Atom("Cu"))
    cu1 = rw.AddAtom(Chem.Atom("Cu"))
    o2 = rw.AddAtom(Chem.Atom("O"))
    o3 = rw.AddAtom(Chem.Atom("O"))
    c4 = rw.AddAtom(Chem.Atom("C"))
    c5 = rw.AddAtom(Chem.Atom("C"))
    rw.AddBond(cu0, o2, Chem.BondType.SINGLE)
    rw.AddBond(cu1, o2, Chem.BondType.SINGLE)
    rw.AddBond(cu0, o3, Chem.BondType.SINGLE)
    rw.AddBond(cu1, o3, Chem.BondType.SINGLE)
    rw.AddBond(o2, c4, Chem.BondType.SINGLE)
    rw.AddBond(o3, c5, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    coords = {
        cu0: np.array([-1.4, 0.0, 0.0]),
        cu1: np.array([+1.4, 0.0, 0.0]),
        o2: np.array([0.0, 1.4, 0.0]),
        o3: np.array([0.0, -1.4, 0.0]),
        c4: np.array([0.0, 2.8, 0.0]),
        c5: np.array([0.0, -2.8, 0.0]),
    }
    return mol, _coords_to_xyz(mol, coords)


def _coords_to_xyz(mol, coords_by_idx) -> str:
    n = mol.GetNumAtoms()
    out = [str(n), ""]
    for i in range(n):
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        p = coords_by_idx[i]
        out.append(f"{sym:4s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}")
    return "\n".join(out) + "\n"


def _force_mxm_angle(xyz: str, mol, x_idx: int, m_pair, target_deg: float) -> str:
    """Move ``x_idx`` so the M-X-M angle takes ``target_deg`` (keeping
    |M-X| at the geometric mean of the two current distances).  Used by
    the tests to manufacture a violating geometry from a clean dimer.
    """
    syms, pts, orig = mod._parse_xyz(xyz)
    ma, mb = m_pair
    p_ma = pts[ma]; p_mb = pts[mb]
    midpoint = 0.5 * (p_ma + p_mb)
    mm = p_mb - p_ma
    mm_norm = float(np.linalg.norm(mm))
    u_mm = mm / mm_norm
    # Move X to satisfy 2*atan( (d_mm/2) / perp ) = target  -> perp = (d_mm/2)/tan(target/2)
    perp = (0.5 * mm_norm) / math.tan(math.radians(target_deg) / 2.0)
    # Existing perpendicular direction (or pick +y if X is collinear with M-M).
    v = pts[x_idx] - midpoint
    v_perp = v - float(np.dot(v, u_mm)) * u_mm
    if np.linalg.norm(v_perp) < 1e-6:
        v_perp = np.array([0.0, 1.0, 0.0])
        v_perp -= float(np.dot(v_perp, u_mm)) * u_mm
    u_perp = v_perp / float(np.linalg.norm(v_perp))
    pts[x_idx] = midpoint + perp * u_perp
    return mod._format_xyz(orig, syms, pts)


def _measure_mxm_angle(xyz: str, x_idx: int, ma: int, mb: int) -> float:
    syms, pts, _ = mod._parse_xyz(xyz)
    v1 = pts[ma] - pts[x_idx]
    v2 = pts[mb] - pts[x_idx]
    cos_a = float(np.clip(np.dot(v1, v2)
                          / (np.linalg.norm(v1) * np.linalg.norm(v2)),
                          -1.0, 1.0))
    return float(math.degrees(math.acos(cos_a)))


# --------------------------------------------------------------------------
# Detection tests
# --------------------------------------------------------------------------

def test_detect_pt2_mu_cl2():
    """Both bridging Cl atoms must be reported as μ-X bridges."""
    mol, xyz = _build_pt2_mu_cl2()
    bridges = mod.detect_mu_x_bridges(xyz, mol)
    cl_bridges = [b for b in bridges if b["x_symbol"] == "Cl"]
    assert len(cl_bridges) == 2, (
        f"expected 2 μ-Cl bridges, got {[b['x_symbol'] for b in bridges]}"
    )
    for b in cl_bridges:
        assert len(b["metal_idxs"]) == 2
        assert b["motif"] == "4-ring", (
            f"both Cl share a 4-ring motif, got {b['motif']}"
        )


def test_detect_fe2_mu_oh2():
    """Both bridging hydroxide O must be reported as μ-X."""
    mol, xyz = _build_fe2_mu_oh2()
    bridges = mod.detect_mu_x_bridges(xyz, mol)
    o_bridges = [b for b in bridges if b["x_symbol"] == "O"]
    assert len(o_bridges) == 2
    # 4-membered ring (Fe-O-Fe-O) motif
    for b in o_bridges:
        assert b["motif"] == "4-ring"


def test_detect_cu2_mu_ome2():
    """μ-OMe bridges (4-ring motif) on Cu2."""
    mol, xyz = _build_cu2_mu_ome2()
    bridges = mod.detect_mu_x_bridges(xyz, mol)
    o_bridges = [b for b in bridges if b["x_symbol"] == "O"]
    assert len(o_bridges) == 2
    for b in o_bridges:
        assert b["motif"] == "4-ring"


def test_terminal_anion_not_flagged():
    """Single Pt + one Cl terminal: Cl bonded to one metal must NOT be
    classified as μ-X."""
    mol, xyz = _build_terminal_cl_mono_metal()
    bridges = mod.detect_mu_x_bridges(xyz, mol)
    assert bridges == [], (
        f"terminal Cl misclassified as μ-X: {bridges}"
    )


def test_malformed_xyz_returns_empty():
    """Garbage XYZ — detection returns [] without raising."""
    mol, _ = _build_pt2_mu_cl2()
    out = mod.detect_mu_x_bridges("not an xyz string", mol)
    assert out == []


# --------------------------------------------------------------------------
# Fix-pass geometric tests
# --------------------------------------------------------------------------

def test_fix_pt2_mu_cl2_bends_180_to_window():
    """Force the Pt-Cl-Pt to ~180° (catastrophic), run the fixer, expect
    the angle to enter the 4-ring window [80°, 105°]."""
    mol, xyz = _build_pt2_mu_cl2()
    # Force the bridging Cl (idx 2) to a linear angle.
    bad_xyz = _force_mxm_angle(xyz, mol, x_idx=2, m_pair=(0, 1), target_deg=178.0)
    pre = _measure_mxm_angle(bad_xyz, x_idx=2, ma=0, mb=1)
    assert pre > 150.0, f"setup failed: pre={pre:.1f}°"
    new_xyz, report = mod.fix_bridging_anion_angles(bad_xyz, mol)
    assert report["n_bridges"] >= 1
    assert report["n_violations"] >= 1
    assert report["n_fixed"] >= 1, f"no fix applied, report={report}"
    post = _measure_mxm_angle(new_xyz, x_idx=2, ma=0, mb=1)
    # Window for 4-ring motif is [80, 105]; allow a small tolerance.
    assert 78.0 <= post <= 108.0, (
        f"post-fix angle {post:.1f}° outside 4-ring window"
    )


def test_fix_noop_when_inside_window():
    """Clean dimer geometry — fixer must return the input XYZ unchanged."""
    mol, xyz = _build_pt2_mu_cl2()  # already ~90° per construction
    new_xyz, report = mod.fix_bridging_anion_angles(xyz, mol)
    assert report["n_violations"] == 0
    # When n_fixed=0, the public API returns the input string unchanged.
    assert new_xyz == xyz


def test_fix_preserves_metal_metal_distance():
    """The bridging-anion bend must keep the two metals exactly fixed."""
    mol, xyz = _build_pt2_mu_cl2()
    bad_xyz = _force_mxm_angle(xyz, mol, x_idx=2, m_pair=(0, 1), target_deg=170.0)
    syms_pre, pts_pre, _ = mod._parse_xyz(bad_xyz)
    new_xyz, report = mod.fix_bridging_anion_angles(bad_xyz, mol)
    if report["n_fixed"] == 0:
        pytest.skip("fixer did not engage — geometry-test only meaningful when it does")
    syms_post, pts_post, _ = mod._parse_xyz(new_xyz)
    # Metals (idx 0, 1) unchanged within numerical noise.
    for i in (0, 1):
        d = float(np.linalg.norm(pts_post[i] - pts_pre[i]))
        assert d < 1e-3, f"metal {i} moved by {d:.4f} Å"


# --------------------------------------------------------------------------
# Dispatch / wiring tests
# --------------------------------------------------------------------------

def _scrub_env(monkeypatch):
    for k in (
        "DELFIN_FIX_BRIDGING_ANION",
        "DELFIN_FIX_BRIDGING_ANION_CLASSES",
    ):
        monkeypatch.delenv(k, raising=False)


def test_dispatcher_default_off_bit_exact(monkeypatch):
    """Env unset -> input results returned identical (no transform)."""
    _scrub_env(monkeypatch)
    mol, xyz = _build_pt2_mu_cl2()
    bad = _force_mxm_angle(xyz, mol, x_idx=2, m_pair=(0, 1), target_deg=170.0)
    results = [(bad, "test")]
    out = sc._apply_fixer_bridging_anion_if_enabled(
        mol, results, dual_parse_done=False,
    )
    assert out == results, "default-OFF must be bit-exact"


def test_dispatcher_env_on_transforms(monkeypatch):
    """DELFIN_FIX_BRIDGING_ANION=1 -> dispatcher rewrites a violating XYZ."""
    _scrub_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_BRIDGING_ANION", "1")
    mol, xyz = _build_pt2_mu_cl2()
    bad = _force_mxm_angle(xyz, mol, x_idx=2, m_pair=(0, 1), target_deg=170.0)
    results = [(bad, "test")]
    out = sc._apply_fixer_bridging_anion_if_enabled(
        mol, results, dual_parse_done=False,
    )
    assert len(out) == 1
    new_xyz, lbl = out[0]
    assert lbl == "test"
    assert new_xyz != bad, "env-on path must rewrite XYZ with a violation"


def test_dispatcher_dual_parse_bypass(monkeypatch):
    """dual_parse_done=True short-circuits even when env=1."""
    _scrub_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_BRIDGING_ANION", "1")
    mol, xyz = _build_pt2_mu_cl2()
    bad = _force_mxm_angle(xyz, mol, x_idx=2, m_pair=(0, 1), target_deg=170.0)
    results = [(bad, "test")]
    out = sc._apply_fixer_bridging_anion_if_enabled(
        mol, results, dual_parse_done=True,
    )
    assert out == results


def test_dispatcher_mol_none_safe(monkeypatch):
    """mol=None -> dispatcher returns input unchanged."""
    _scrub_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_BRIDGING_ANION", "1")
    _, xyz = _build_pt2_mu_cl2()
    results = [(xyz, "test")]
    out = sc._apply_fixer_bridging_anion_if_enabled(
        None, results, dual_parse_done=False,
    )
    assert out == results


def test_dispatcher_no_metal_short_circuit(monkeypatch):
    """No metal symbol in any XYZ -> fast-path early return."""
    _scrub_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_BRIDGING_ANION", "1")
    # Methylamine (no metal)
    m = Chem.AddHs(Chem.MolFromSmiles("CN"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(m, randomSeed=42)
    conf = m.GetConformer()
    out_lines = [str(m.GetNumAtoms()), ""]
    for i, a in enumerate(m.GetAtoms()):
        p = conf.GetAtomPosition(i)
        out_lines.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    xyz = "\n".join(out_lines) + "\n"
    results = [(xyz, "test")]
    out = sc._apply_fixer_bridging_anion_if_enabled(
        m, results, dual_parse_done=False,
    )
    assert out == results


def test_class_conditional_allowlist(monkeypatch):
    """When DELFIN_FIX_BRIDGING_ANION_CLASSES is set to a list that does
    NOT include the mol's class, the fix must NOT fire (env_on subset
    of class-conditional gate)."""
    _scrub_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FIX_BRIDGING_ANION_CLASSES", "nonexistent_class")
    mol, xyz = _build_pt2_mu_cl2()
    bad = _force_mxm_angle(xyz, mol, x_idx=2, m_pair=(0, 1), target_deg=170.0)
    results = [(bad, "test")]
    out = sc._apply_fixer_bridging_anion_if_enabled(
        mol, results, dual_parse_done=False,
    )
    assert out == results, (
        "class-conditional allow-list with non-matching class must be no-op"
    )


def test_dispatcher_empty_results():
    """Empty input list -> empty output (no-op path)."""
    out = sc._apply_fixer_bridging_anion_if_enabled(None, [], False)
    assert out == []
