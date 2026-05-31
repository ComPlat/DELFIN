"""Tests for delfin.fffree.aromatic_snap (Phase G12)."""
import os
import sys
import importlib

import numpy as np
import pytest


def _fresh_import(flag_pt3, flag_aromsnap):
    """Re-import aromatic_snap with given env state."""
    if flag_pt3:
        os.environ["DELFIN_FFFREE_PURE_TRACK3"] = "1"
    else:
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
    if flag_aromsnap:
        os.environ["DELFIN_FFFREE_AROMSNAP"] = "1"
    else:
        os.environ.pop("DELFIN_FFFREE_AROMSNAP", None)
    for m in list(sys.modules):
        if "aromatic_snap" in m:
            sys.modules.pop(m, None)
    from delfin.fffree import aromatic_snap as AS
    return AS


def test_aromsnap_off_by_default():
    """Byte-identical default-OFF when neither PT3 nor AROMSNAP set."""
    AS = _fresh_import(False, False)
    syms = ["C"] * 6 + ["H"] * 6
    th = np.linspace(0, 2 * np.pi, 7)[:-1]
    P0 = np.array(
        [[1.4 * np.cos(t), 1.4 * np.sin(t), 0.3 * (i % 2)] for i, t in enumerate(th)]
        + [[2.5 * np.cos(t), 2.5 * np.sin(t), 0.0] for t in th]
    )
    syms_out, P_out = AS.snap_aromatic_rings(syms, P0)
    assert syms_out == syms
    assert np.allclose(P_out, P0)


def test_aromsnap_on_under_pt3():
    """PT3 alone auto-enables AROMSNAP (race-strategy: FF-free uses snap)."""
    AS = _fresh_import(True, False)
    assert AS.is_enabled() is True


def test_aromsnap_explicit_flag_works():
    """Explicit DELFIN_FFFREE_AROMSNAP=1 enables when PT3 is off."""
    AS = _fresh_import(False, True)
    assert AS.is_enabled() is True


def test_benzene_buckle_flattens():
    """Buckled benzene -> planar via SVD plane snap, C-H bonds preserved."""
    AS = _fresh_import(False, True)
    syms = ["C"] * 6 + ["H"] * 6
    th = np.linspace(0, 2 * np.pi, 7)[:-1]
    P0 = np.array(
        [[1.4 * np.cos(t), 1.4 * np.sin(t), 0.3 * (i % 2)] for i, t in enumerate(th)]
        + [[2.5 * np.cos(t), 2.5 * np.sin(t), 0.0] for t in th]
    )
    syms_out, P_out = AS.snap_aromatic_rings(syms, P0)
    # Planarity check
    pts = P_out[:6]
    centered = pts - pts.mean(axis=0)
    _, _, Vt = np.linalg.svd(centered, full_matrices=False)
    offsets = centered @ Vt[-1]
    rms = float(np.sqrt(float((offsets * offsets).mean())))
    assert rms < 1e-6, f"benzene not planar after snap: rms={rms}"
    # C-H bond lengths preserved (rigid drag of H)
    for k in range(6):
        d_before = float(np.linalg.norm(P0[k + 6] - P0[k]))
        d_after = float(np.linalg.norm(P_out[k + 6] - P_out[k]))
        assert abs(d_after - d_before) < 1e-9


def test_already_planar_no_op():
    """Already planar ring -> no movement (RMS below threshold)."""
    AS = _fresh_import(False, True)
    syms = ["C"] * 6 + ["H"] * 6
    th = np.linspace(0, 2 * np.pi, 7)[:-1]
    P0 = np.array(
        [[1.4 * np.cos(t), 1.4 * np.sin(t), 0.0] for t in th]
        + [[2.5 * np.cos(t), 2.5 * np.sin(t), 0.0] for t in th]
    )
    _, P_out = AS.snap_aromatic_rings(syms, P0)
    assert np.allclose(P_out, P0, atol=1e-9)


def test_metallacycle_skipped():
    """Ring with a metal atom is NOT snapped (universal guard)."""
    AS = _fresh_import(False, True)
    # 5-ring with one metal -> aromatic_ring_filter excludes it
    syms = ["Pd", "C", "C", "C", "C", "H", "H", "H", "H"]
    # Pentagon, one vertex is Pd
    pts = np.array(
        [[1.4 * np.cos(t), 1.4 * np.sin(t), 0.0 if i == 0 else 0.5 * (i % 2)]
         for i, t in enumerate(np.linspace(0, 2 * np.pi, 6)[:-1])]
    )
    P0 = np.vstack([pts, np.zeros((4, 3))])
    P_out = AS.snap_aromatic_rings(syms, P0)[1]
    # Pd untouched (metal frozen)
    assert np.allclose(P_out[0], P0[0])
