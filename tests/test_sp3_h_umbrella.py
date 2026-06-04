"""Tests for delfin/fffree/sp3_h_umbrella.py — build-time sp³ umbrella."""
from __future__ import annotations

import math
import os
import sys
import hashlib
import importlib
import io

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from delfin.fffree import sp3_h_umbrella as U  # noqa: E402

_TETRA = U._IDEAL_TETRA_DEG


def _angle(p1, vertex, p2):
    v1 = np.asarray(p1) - np.asarray(vertex)
    v2 = np.asarray(p2) - np.asarray(vertex)
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    cos = float(np.dot(v1, v2) / (n1 * n2))
    cos = max(-1.0, min(1.0, cos))
    return math.degrees(math.acos(cos))


def _make_degenerate_ch3():
    """Return (syms, coords) for a CH3 in the 90/180 broken pattern."""
    c = np.array([0.0, 0.0, 0.0])
    x = np.array([0.0, 0.0, -1.54])  # heavy C neighbour
    h0 = c + np.array([1.09, 0.0, 0.0])
    h1 = c + np.array([-1.09, 0.0, 0.0])  # 180° to h0
    h2 = c + np.array([0.0, 1.09, 0.0])    # 90° to both
    return ["C", "C", "H", "H", "H"], np.vstack([c, x, h0, h1, h2])


def _make_ideal_ch3():
    """Return (syms, coords) for a Td-perfect CH3."""
    c = np.array([0.0, 0.0, 0.0])
    x = np.array([0.0, 0.0, -1.54])
    # Build via the same primitive
    arr = U._build_methyl_umbrella(c, x, [1.09, 1.09, 1.09])
    syms = ["C", "C", "H", "H", "H"]
    coords = np.vstack([c, x, arr])
    return syms, coords


# ---------------------------------------------------------------------------
# Tests 1-3 — detection
# ---------------------------------------------------------------------------


def test_detect_degenerate_methyl_flag():
    syms, coords = _make_degenerate_ch3()
    centres = U.detect_sp3_centers(coords, syms)
    assert len(centres) == 1
    c = centres[0]
    assert c.classification == "CH3"
    assert c.center_idx == 0
    assert c.h_indices == (2, 3, 4)
    assert c.heavy_indices == (1,)
    diag = U.check_umbrella_geometry(c, coords)
    assert diag.flag_degenerate is True
    # H-H pair includes one ~180° -> max_abs_cos ≥ 0.99
    assert diag.max_hh_cos_abs > 0.99


def test_detect_ideal_methyl_pass():
    syms, coords = _make_ideal_ch3()
    centres = U.detect_sp3_centers(coords, syms)
    assert len(centres) == 1
    diag = U.check_umbrella_geometry(centres[0], coords)
    assert diag.flag_degenerate is False
    # Td angles all within 0.1° of 109.47°
    for a in diag.hh_angles_deg:
        assert abs(a - _TETRA) < 0.1


def test_no_h_no_center():
    """Centre with no H neighbours is not enumerated."""
    syms = ["C", "C"]
    coords = np.array([[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]])
    centres = U.detect_sp3_centers(coords, syms)
    assert centres == []


# ---------------------------------------------------------------------------
# Tests 4-7 — enforce
# ---------------------------------------------------------------------------


def test_enforce_methyl_recovers_tetrahedral():
    syms, coords = _make_degenerate_ch3()
    centres = U.detect_sp3_centers(coords, syms)
    new_coords = U.enforce_umbrella(coords, centres[0])
    diag = U.check_umbrella_geometry(centres[0], new_coords)
    # All H-C-H angles within (105°, 115°)
    for a in diag.hh_angles_deg:
        assert 105.0 < a < 115.0
    assert diag.flag_degenerate is False


def test_enforce_preserves_bond_lengths():
    syms, coords = _make_degenerate_ch3()
    centres = U.detect_sp3_centers(coords, syms)
    # Set individual bond lengths to distinct values 1.10, 1.11, 1.12 to check preservation
    coords[2] = coords[0] + 1.10 * (coords[2] - coords[0]) / np.linalg.norm(coords[2] - coords[0])
    coords[3] = coords[0] + 1.11 * (coords[3] - coords[0]) / np.linalg.norm(coords[3] - coords[0])
    coords[4] = coords[0] + 1.12 * (coords[4] - coords[0]) / np.linalg.norm(coords[4] - coords[0])
    new = U.enforce_umbrella(coords, centres[0], preserve_bond_lengths=True)
    for h_idx, expected in zip([2, 3, 4], [1.10, 1.11, 1.12]):
        d = float(np.linalg.norm(new[h_idx] - new[0]))
        assert abs(d - expected) < 1e-3, f"H[{h_idx}] bond {d:.4f} ≠ {expected:.2f}"


def test_enforce_does_not_move_center_or_heavy():
    syms, coords = _make_degenerate_ch3()
    centres = U.detect_sp3_centers(coords, syms)
    snap = coords.copy()
    new = U.enforce_umbrella(coords, centres[0])
    assert np.allclose(new[0], snap[0])  # centre
    assert np.allclose(new[1], snap[1])  # heavy neighbour


def test_enforce_methyl_does_not_move_unrelated_atoms():
    """Add an unrelated atom and ensure it stays put."""
    syms, coords = _make_degenerate_ch3()
    syms.append("O")
    extra = np.array([[5.0, 5.0, 5.0]])
    coords = np.vstack([coords, extra])
    centres = U.detect_sp3_centers(coords, syms)
    new = U.enforce_umbrella(coords, centres[0])
    assert np.allclose(new[-1], extra[0])


# ---------------------------------------------------------------------------
# Tests 8-11 — determinism, byte-identity
# ---------------------------------------------------------------------------


def test_determinism_same_input_same_output():
    syms, coords = _make_degenerate_ch3()
    centres = U.detect_sp3_centers(coords, syms)
    r1 = U.enforce_umbrella(coords, centres[0])
    r2 = U.enforce_umbrella(coords, centres[0])
    assert np.array_equal(r1, r2)


def test_byte_identical_default_off():
    """When DELFIN_FFFREE_SP3_H_UMBRELLA is unset, ``umbrella_active`` is False."""
    os.environ.pop("DELFIN_FFFREE_SP3_H_UMBRELLA", None)
    assert U.umbrella_active() is False


def test_env_flag_on_recognised():
    os.environ["DELFIN_FFFREE_SP3_H_UMBRELLA"] = "1"
    assert U.umbrella_active() is True
    os.environ["DELFIN_FFFREE_SP3_H_UMBRELLA"] = "true"
    assert U.umbrella_active() is True
    os.environ["DELFIN_FFFREE_SP3_H_UMBRELLA"] = "yes"
    assert U.umbrella_active() is True
    os.environ["DELFIN_FFFREE_SP3_H_UMBRELLA"] = "0"
    assert U.umbrella_active() is False
    os.environ.pop("DELFIN_FFFREE_SP3_H_UMBRELLA", None)


def test_deterministic_basis_is_function_of_axis_only():
    """Same axis, same basis — order-independent (the determinism contract)."""
    r1 = np.array([1.0, 0.0, 0.0])
    u1, v1 = U._deterministic_basis_for_axis(r1)
    u2, v2 = U._deterministic_basis_for_axis(r1.copy())
    assert np.array_equal(u1, u2)
    assert np.array_equal(v1, v2)
    # Orthogonality contract
    assert abs(float(np.dot(u1, r1))) < 1e-12
    assert abs(float(np.dot(v1, r1))) < 1e-12
    assert abs(float(np.dot(u1, v1))) < 1e-12
    # Near-z axis triggers the ŷ fallback
    rz = np.array([0.0, 0.0, 1.0])
    uz, vz = U._deterministic_basis_for_axis(rz)
    assert abs(float(np.dot(uz, rz))) < 1e-12
    assert abs(float(np.dot(vz, rz))) < 1e-12


# ---------------------------------------------------------------------------
# Tests 12-14 — CH2 / CH / NH classifications
# ---------------------------------------------------------------------------


def test_ch2_classification_and_enforce():
    """Two heavies + two H's."""
    # Glycine-like backbone fragment: N - CH2 - C
    syms = ["N", "C", "C", "H", "H"]
    c_pos = np.array([0.0, 0.0, 0.0])
    n_pos = np.array([-1.45, 0.0, 0.0])
    c2_pos = np.array([1.50, 0.0, 0.0])
    # Two H's stuck on the same axis as the heavies (broken)
    h0 = c_pos + np.array([0.0, 0.0, 1.09])
    h1 = c_pos + np.array([0.0, 0.0, -1.09])  # 180° to h0
    coords = np.vstack([n_pos, c_pos, c2_pos, h0, h1])
    centres = U.detect_sp3_centers(coords, syms)
    assert len(centres) == 1
    assert centres[0].classification == "CH2"
    diag = U.check_umbrella_geometry(centres[0], coords)
    assert diag.flag_degenerate  # 180° between the two H's
    new = U.enforce_umbrella(coords, centres[0])
    a = _angle(new[3], new[1], new[4])
    assert 95.0 < a < 125.0


def test_ch_methine_classification_and_enforce():
    """3 heavies (non-coplanar) + 1 H."""
    c_pos = np.array([0.0, 0.0, 0.0])
    # 3 heavies in a tetrahedral 3-of-4 arrangement (true tetrahedral
    # vertices); the missing 4th vertex is along +z by construction.
    a = np.array([1.50, 0.0, -0.5])
    b = np.array([-0.75, 1.3, -0.5])
    d = np.array([-0.75, -1.3, -0.5])
    h = np.array([0.0, 0.0, 1.09])
    syms = ["C", "C", "C", "C", "H"]
    coords = np.vstack([c_pos, a, b, d, h])
    centres = U.detect_sp3_centers(coords, syms)
    assert len(centres) == 1
    assert centres[0].classification == "CH"
    new = U.enforce_umbrella(coords, centres[0])
    h_new = new[4]
    direction = (h_new - c_pos) / np.linalg.norm(h_new - c_pos)
    # Out-of-plane component dominates (umbrella closure points +z because
    # the 3 heavies all have z = -0.5)
    assert direction[2] > 0.9


def test_nh2_classification():
    syms = ["N", "C", "H", "H"]
    n_pos = np.array([0.0, 0.0, 0.0])
    c_pos = np.array([0.0, 0.0, -1.45])
    h0 = n_pos + np.array([1.01, 0.0, 0.0])
    h1 = n_pos + np.array([-1.01, 0.0, 0.0])
    coords = np.vstack([n_pos, c_pos, h0, h1])
    centres = U.detect_sp3_centers(coords, syms)
    assert len(centres) == 1
    assert centres[0].classification == "NH2"


# ---------------------------------------------------------------------------
# Tests 15-18 — batch enforce_all_sp3_umbrella
# ---------------------------------------------------------------------------


def test_enforce_all_only_degenerate_skips_ideal():
    """When ``only_degenerate=True``, ideal centres are NOT rewritten."""
    syms, coords = _make_ideal_ch3()
    new, rep = U.enforce_all_sp3_umbrella(coords, syms, only_degenerate=True)
    assert rep["degenerate_flagged"] == 0
    assert rep["centres_rewritten"] == 0
    # Coords byte-identical
    assert np.array_equal(new, coords)


def test_enforce_all_rewrites_degenerate():
    syms, coords = _make_degenerate_ch3()
    new, rep = U.enforce_all_sp3_umbrella(coords, syms, only_degenerate=True)
    assert rep["degenerate_flagged"] == 1
    assert rep["centres_rewritten"] == 1
    assert rep["hs_moved"] == 3
    assert rep["max_hh_dev_after_deg"] < 1.0


def test_enforce_all_xyz_wrapper():
    """XYZ string in -> XYZ string out."""
    syms, coords = _make_degenerate_ch3()
    xyz = f"{len(syms)}\ntest\n"
    for s, c in zip(syms, coords):
        xyz += f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n"
    new_xyz, rep = U.enforce_all_sp3_umbrella_xyz(xyz)
    assert rep["centres_rewritten"] == 1
    assert new_xyz != xyz


def test_enforce_all_zero_atoms_safe():
    new, rep = U.enforce_all_sp3_umbrella(np.zeros((0, 3)), [])
    assert rep["sp3_centres_seen"] == 0
    assert new.shape == (0, 3)


# ---------------------------------------------------------------------------
# Tests 19-22 — edge cases
# ---------------------------------------------------------------------------


def test_methyl_on_chiral_center():
    """Methyl attached to a chiral C — no chirality flip from the umbrella
    rewrite because the centre itself is not moved.
    """
    # Build a stereocentre C with CH3, OH, NH2, COOH-like substituents
    c = np.array([0.0, 0.0, 0.0])
    me_c = np.array([1.54, 0.0, 0.0])      # CH3 carbon (the sp³ centre)
    o = np.array([-0.5, 1.4, 0.0])
    n = np.array([-0.5, -0.7, 1.2])
    cd = np.array([-0.5, -0.7, -1.2])
    # Degenerate methyl: 3 H's on orthogonal axes around me_c
    h0 = me_c + np.array([0.0, 1.09, 0.0])
    h1 = me_c + np.array([0.0, -1.09, 0.0])
    h2 = me_c + np.array([0.0, 0.0, 1.09])
    syms = ["C", "C", "O", "N", "C", "H", "H", "H"]
    coords = np.vstack([c, me_c, o, n, cd, h0, h1, h2])
    centres = U.detect_sp3_centers(coords, syms)
    # me_c (idx=1) is the methyl carbon
    methyl_centres = [c for c in centres if c.classification == "CH3"]
    assert len(methyl_centres) == 1
    new = U.enforce_umbrella(coords, methyl_centres[0])
    # The stereocentre c (idx=0) is not moved
    assert np.allclose(new[0], c)
    # The heavy substituents of c are not moved
    assert np.allclose(new[2], o)
    assert np.allclose(new[3], n)
    assert np.allclose(new[4], cd)


def test_missing_h_no_crash():
    """A 'CH3' centre with only 2 H's becomes 'OTHER' and is skipped."""
    syms = ["C", "C", "H", "H"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -1.54],
        [1.09, 0.0, 0.0],
        [-1.09, 0.0, 0.0],  # 180°
    ])
    # Without include_other this is CH2 if heavy_count = 2; here it's heavy=1, h=2 → "OTHER"
    centres = U.detect_sp3_centers(coords, syms, include_other=False)
    # heavy=1, h=2 -> not a canonical class -> skipped
    assert all(c.classification != "OTHER" for c in centres)


def test_partial_h_presence_ch2_only():
    """CH2 with 2 heavies + 2 H's — classified as CH2 even if H's degenerate."""
    syms = ["C", "C", "C", "H", "H"]
    c = np.array([0.0, 0.0, 0.0])
    x = np.array([-1.54, 0.0, 0.0])
    y = np.array([1.54, 0.0, 0.0])
    h0 = c + np.array([0.0, 0.0, 1.09])
    h1 = c + np.array([0.0, 0.0, -1.09])
    coords = np.vstack([c, x, y, h0, h1])
    centres = U.detect_sp3_centers(coords, syms)
    assert any(c.classification == "CH2" for c in centres)


def test_rollback_when_clash_introduced():
    """Test rollback by calling enforce_all with a pre-built degenerate
    methyl and a foreign atom that will overlap with the ideal H
    position once the umbrella is rebuilt.  We feed bond_topology
    explicitly so the foreign atom doesn't break the CH3 detection."""
    # Build the degenerate methyl
    c = np.array([0.0, 0.0, 0.0])
    x = np.array([0.0, 0.0, -1.54])
    h0 = c + np.array([1.09, 0.0, 0.0])
    h1 = c + np.array([-1.09, 0.0, 0.0])
    h2 = c + np.array([0.0, 1.09, 0.0])
    # Compute ideal landing of H[2] for this methyl
    centre = U.Sp3Center(0, "C", (2, 3, 4), (1,), "CH3")
    coords_pre = np.vstack([c, x, h0, h1, h2])
    ideal = U.enforce_umbrella(coords_pre, centre)
    # Place a heavy carbon at the IDEAL H[2] position
    foreign_c = ideal[2].copy()
    syms2 = ["C", "C", "H", "H", "H", "C"]
    coords2 = np.vstack([coords_pre, foreign_c])
    # Provide explicit bond_topology so detection sees the original CH3
    # (no edge from C[0] to the foreign C[5])
    topo = [[1, 2, 3, 4], [0], [0], [0], [0], []]
    new, rep = U.enforce_all_sp3_umbrella(
        coords2, syms2, bond_topology=topo, enable_clash_rollback=True,
    )
    # The umbrella tries to place H[2] (or H[3] / H[4]) where the foreign
    # carbon already sits -> rollback should fire on at least 1 H.
    assert rep["hs_rolled_back"] >= 1


# ---------------------------------------------------------------------------
# Test 23 — accept-if-better at heal-pathway level (smoke for sp3_h_heal)
# ---------------------------------------------------------------------------


def test_heal_accept_if_better_gate():
    from delfin.fffree.sp3_h_heal import heal_degenerate_sp3_h
    syms, coords = _make_degenerate_ch3()
    healed, rep = heal_degenerate_sp3_h(coords, mol=None, syms=syms,
                                         accept_if_better=True)
    # Original max H-H dev > 60°, healed should be near 0
    assert rep["max_hh_dev_before_deg"] > 60.0
    assert rep["max_hh_dev_after_deg"] < 1.0
    assert rep["accepted"] is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
