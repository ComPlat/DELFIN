"""Mission 3 — Tests for delfin.fffree.xrd_recall_metric.

These cover:
- isomer classifier (fac/mer/cis/trans/OH/SP)
- Kabsch heavy-atom RMSD reflexivity + translation/rotation invariance
- pool scorer with synthetic mini-reference
"""
from __future__ import annotations

import json
import math
import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

import delfin.fffree.xrd_recall_metric as xrd


# -------------------- isomer classifier --------------------
def _make_oh_fac():
    """Mer/fac heuristic example. Octahedral M with 3xN + 3xO; fac = all three
    N occupy one face (mutually cis), mer = one N pair trans."""
    # M at origin; +x, +y, +z = N (fac); -x, -y, -z = O
    syms = ["Fe", "N", "N", "N", "O", "O", "O"]
    P = np.array([
        [0, 0, 0],
        [2.0, 0, 0],
        [0, 2.0, 0],
        [0, 0, 2.0],
        [-2.0, 0, 0],
        [0, -2.0, 0],
        [0, 0, -2.0],
    ], dtype=float)
    return syms, P


def _make_oh_mer():
    """Mer: two N trans (along z) and one N perpendicular."""
    syms = ["Fe", "N", "N", "N", "O", "O", "O"]
    P = np.array([
        [0, 0, 0],
        [2.0, 0, 0],   # N1 +x
        [-2.0, 0, 0],  # N2 -x  (trans pair)
        [0, 2.0, 0],   # N3 +y
        [0, -2.0, 0],
        [0, 0, 2.0],
        [0, 0, -2.0],
    ], dtype=float)
    return syms, P


def _make_sp_cis():
    """Square-planar Pd with 2xN + 2xCl, cis."""
    syms = ["Pd", "N", "N", "Cl", "Cl"]
    P = np.array([
        [0, 0, 0],
        [2.0, 0, 0],
        [0, 2.0, 0],
        [-2.0, 0, 0],
        [0, -2.0, 0],
    ], dtype=float)
    return syms, P


def _make_sp_trans():
    """Square-planar Pd with 2xN + 2xCl, trans."""
    syms = ["Pd", "N", "N", "Cl", "Cl"]
    P = np.array([
        [0, 0, 0],
        [2.0, 0, 0],
        [-2.0, 0, 0],  # N trans
        [0, 2.0, 0],
        [0, -2.0, 0],  # Cl trans
    ], dtype=float)
    return syms, P


def test_classify_fac():
    syms, P = _make_oh_fac()
    lbl = xrd.classify_isomer(syms, P)
    assert lbl == "CN6-OH-fac", f"got {lbl}"


def test_classify_mer():
    syms, P = _make_oh_mer()
    lbl = xrd.classify_isomer(syms, P)
    assert lbl == "CN6-OH-mer", f"got {lbl}"


def test_classify_sp_cis():
    syms, P = _make_sp_cis()
    lbl = xrd.classify_isomer(syms, P)
    assert lbl == "CN4-tet-or-sp-cis", f"got {lbl}"


def test_classify_sp_trans():
    syms, P = _make_sp_trans()
    lbl = xrd.classify_isomer(syms, P)
    assert lbl == "CN4-tet-or-sp-trans", f"got {lbl}"


def test_classify_no_metal():
    syms = ["C", "C", "O", "H"]
    P = np.zeros((4, 3))
    P[1] = [1.5, 0, 0]
    P[2] = [3.0, 0, 0]
    P[3] = [3.6, 0.8, 0]
    lbl = xrd.classify_isomer(syms, P)
    assert lbl == "no-metal"


# -------------------- kabsch RMSD --------------------
def test_kabsch_rmsd_reflexive():
    syms, P = _make_oh_fac()
    r = xrd.kabsch_rmsd_heavy(syms, P, syms, P)
    assert r < 1e-6


def test_kabsch_rmsd_translation_invariant():
    syms, P = _make_oh_fac()
    Q = P + np.array([5.0, -3.0, 2.0])
    r = xrd.kabsch_rmsd_heavy(syms, P, syms, Q)
    assert r < 1e-6


def test_kabsch_rmsd_rotation_invariant():
    syms, P = _make_oh_fac()
    # 90° rotation around z
    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    Q = P @ R.T
    r = xrd.kabsch_rmsd_heavy(syms, P, syms, Q)
    assert r < 1e-6


def test_kabsch_rmsd_perturbation():
    syms, P = _make_oh_fac()
    Q = P.copy()
    Q[1] += np.array([0.1, 0, 0])  # 0.1 Å perturbation on one N
    r = xrd.kabsch_rmsd_heavy(syms, P, syms, Q)
    # 0.1 spread over 7 atoms ~= 0.038
    assert 0.0 < r < 0.1


# -------------------- refcode parser --------------------
def test_refcode_extraction():
    assert xrd._refcode_from_filename("042-ARABUR.xyz") == "ARABUR"
    assert xrd._refcode_from_filename("ARABUR.xyz") == "ARABUR"
    assert xrd._refcode_from_filename("M2-HUCPIH-foo.xyz") == "HUCPIH"
    assert xrd._refcode_from_filename("01-Fe_CO_3.xyz") is None


# -------------------- pool scorer with synthetic data --------------------
def test_score_pool_synthetic(tmp_path):
    # Build a synthetic ground truth with one refcode whose isomer is fac
    syms, P = _make_oh_fac()
    ref = {
        "version": 1,
        "records": [{
            "refcode": "FACTST",
            "ok": True,
            "symbols": list(syms),
            "positions": P.tolist(),
            "isomer_label": "CN6-OH-fac",
            "metal_sym": "Fe",
            "cn": 6,
            "polyhedron_class": "OH",
        }]
    }
    ref_path = tmp_path / "ref.json"
    ref_path.write_text(json.dumps(ref))

    # Pool with one matching fac file and one mer file
    pool = tmp_path / "pool"
    pool.mkdir()
    def _write_xyz(p, syms, P, comment="t"):
        text = f"{len(syms)}\n{comment}\n"
        for s, xyz in zip(syms, P):
            text += f"{s} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n"
        p.write_text(text)
    s1, P1 = _make_oh_fac()
    _write_xyz(pool / "001-FACTST.xyz", s1, P1)
    s2, P2 = _make_oh_mer()
    _write_xyz(pool / "002-FACTST.xyz", s2, P2)

    ref_loaded = json.loads(ref_path.read_text())
    result = xrd.score_pool(pool, ref=ref_loaded, rmsd_threshold=0.5)
    assert "per_refcode" in result
    rec = result["per_refcode"][0]
    assert rec["refcode"] == "FACTST"
    assert rec["isomer_recall"] == 1, f"expected isomer match, got {rec}"
    assert rec["rmsd_pass"] == 1, f"expected RMSD pass, got {rec}"
    assert result["isomer_recall_mean"] == 1.0
    assert result["conformer_recall_mean"] == 1.0


def test_score_pool_no_match(tmp_path):
    syms_f, Pf = _make_oh_fac()
    ref = {
        "version": 1,
        "records": [{
            "refcode": "MERTST",
            "ok": True,
            "symbols": list(syms_f),
            "positions": Pf.tolist(),  # but isomer_label says mer
            "isomer_label": "CN6-OH-mer",
        }]
    }
    pool = tmp_path / "pool"
    pool.mkdir()
    # Only fac structures in pool — no isomer match
    text = f"{len(syms_f)}\nfac\n"
    for s, xyz in zip(syms_f, Pf):
        text += f"{s} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n"
    (pool / "001-MERTST.xyz").write_text(text)
    result = xrd.score_pool(pool, ref=ref, rmsd_threshold=0.5)
    rec = result["per_refcode"][0]
    assert rec["isomer_recall"] == 0
    # RMSD might still pass since coords match — that's correct behavior (geometry
    # is identical but the LABEL doesn't match)
    assert result["isomer_recall_mean"] == 0.0


def test_score_pool_missing_refcode(tmp_path):
    syms, P = _make_oh_fac()
    ref = {
        "version": 1,
        "records": [{
            "refcode": "NOPOOL",
            "ok": True,
            "symbols": list(syms),
            "positions": P.tolist(),
            "isomer_label": "CN6-OH-fac",
        }]
    }
    pool = tmp_path / "pool_empty"
    pool.mkdir()
    result = xrd.score_pool(pool, ref=ref)
    assert result["per_refcode"][0]["n_files"] == 0
    assert result["per_refcode"][0]["isomer_recall"] == 0
    assert result["n_refcodes_missing"] == 1


# -------------------- parse_xyz robustness --------------------
def test_parse_xyz_header():
    text = "3\ncomment\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nO 0.0 1.0 0.0\n"
    syms, P = xrd.parse_xyz(text)
    assert syms == ["C", "H", "O"]
    assert P.shape == (3, 3)


def test_parse_xyz_headerless():
    text = "C 0.0 0.0 0.0\nH 1.0 0.0 0.0\nO 0.0 1.0 0.0\n"
    syms, P = xrd.parse_xyz(text)
    assert syms == ["C", "H", "O"]
    assert P.shape == (3, 3)


def test_env_flag_constant():
    assert xrd.ENV_FLAG == "DELFIN_USE_CCDC_REFERENCE_DUMP"


def test_default_rmsd_threshold():
    assert xrd.RMSD_PASS_THRESHOLD_A == pytest.approx(0.5)


# -------------------- multi-frame parser --------------------
def test_parse_xyz_frames_single():
    text = "3\ncomment\nC 0 0 0\nH 1 0 0\nO 0 1 0\n"
    frames = xrd.parse_xyz_frames(text)
    assert len(frames) == 1
    assert frames[0][0] == ["C", "H", "O"]


def test_parse_xyz_frames_multi():
    """Three-frame XYZ: ensure parse_xyz_frames returns 3 distinct frames."""
    text = (
        "2\nframe1\nC 0 0 0\nH 1 0 0\n"
        "2\nframe2\nC 0 0 0\nH 2 0 0\n"
        "2\nframe3\nN 0 0 0\nO 3 0 0\n"
    )
    frames = xrd.parse_xyz_frames(text)
    assert len(frames) == 3
    assert frames[0][1][1, 0] == pytest.approx(1.0)
    assert frames[1][1][1, 0] == pytest.approx(2.0)
    assert frames[2][0] == ["N", "O"]


def test_parse_xyz_returns_first_frame():
    text = (
        "2\nframe1\nC 0 0 0\nH 1 0 0\n"
        "2\nframe2\nC 0 0 0\nH 9 0 0\n"
    )
    syms, P = xrd.parse_xyz(text)
    assert syms == ["C", "H"]
    assert P[1, 0] == pytest.approx(1.0), \
        "parse_xyz should only return first frame"


# -------------------- CN-sphere RMSD --------------------
def test_cn_rmsd_reflexive():
    syms, P = _make_oh_fac()
    r = xrd.kabsch_rmsd_cn_sphere(syms, P, syms, P)
    assert r < 1e-6


def test_cn_rmsd_no_metal():
    syms = ["C", "C", "O"]
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    r = xrd.kabsch_rmsd_cn_sphere(syms, P, syms, P)
    # no metal → NaN
    assert math.isnan(r)


def test_cn_rmsd_translation_invariant():
    syms, P = _make_oh_fac()
    Q = P + np.array([10.0, -7.0, 4.0])
    r = xrd.kabsch_rmsd_cn_sphere(syms, P, syms, Q)
    assert r < 1e-6
