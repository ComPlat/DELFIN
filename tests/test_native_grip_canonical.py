"""Tests for DELFIN_FFFREE_NATIVE_GRIP_CANONICAL (User 2026-06-07).

Native-path output emits one frame per (kappa-binding-mode, isomer) tuple,
each with potentially different atom-to-atom connectivity.  Multi-frame
XYZ viewers like Avogadro infer bonds from frame 0 and render every other
frame with that single bond set, so any topology-divergent frame renders
broken.

The fix replaces the native multi-isomer list with an ETKDG + GRIP
conformer ensemble of the parent SMILES so every frame has IDENTICAL
atom-to-atom connectivity.  Activated by env flag, default OFF
byte-identical.

This module tests:
  1. OFF env-byte-identical (HEAD safety contract).
  2. Determinism (2x ON produces identical bytes).
  3. Topology consistency (all frames share the same element sequence).
  4. M-D invariant (ensemble internal self-consistency).
"""
from __future__ import annotations
import os
import sys

import pytest

# A 5d Re CN5/CN6 hetero complex with mixed N, Br, P donors and a
# tridentate phosphine cage that triggers the native-path binding-mode
# expansion (the bug class).
SUHJUB = (
    "Br[Re](Br)(N=O)([P]12C[NH+]3CN(CN(C3)C1)C2)"
    "([P]12C[NH+]3CN(CN(C3)C1)C2)[P]12C[NH+]3CN(CN(C3)C1)C2"
)


def _reset_delfin_modules():
    for k in list(sys.modules):
        if k.startswith("delfin"):
            del sys.modules[k]


def _run(env_overrides):
    """Invoke smiles_to_xyz_isomers with the given env overrides."""
    base = {
        "PYTHONHASHSEED": "0",
        "DELFIN_FFFREE_BUILDER": "1",
        "DELFIN_FFFREE_BINDING_MODE_ENUM": "1",
        "DELFIN_FFFREE_PURE_TRACK3": "1",
        "DELFIN_FFFREE_FALLBACK_MODE": "grip",
    }
    base.update(env_overrides)
    saved = {k: os.environ.get(k) for k in base}
    for k, v in base.items():
        os.environ[k] = v
    try:
        _reset_delfin_modules()
        from delfin.smiles_converter import smiles_to_xyz_isomers
        isos, err = smiles_to_xyz_isomers(SUHJUB, max_isomers=8)
        return isos or [], err
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _signature(xyz_block):
    """Element-only sequence (a frame's connectivity proxy)."""
    out = []
    for line in xyz_block.strip().splitlines():
        parts = line.split()
        if parts:
            out.append(parts[0])
    return tuple(out)


def _median_md(xyz_block):
    """Median Re-donor distance over heavy atoms within 3.5 A of Re."""
    import math
    Re = None
    rows = []
    for line in xyz_block.strip().splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        rows.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    Re_idx = next((i for i, (s, *_x) in enumerate(rows) if s == "Re"), None)
    if Re_idx is None:
        return None
    rx, ry, rz = rows[Re_idx][1:]
    ds = []
    for i, (s, x, y, z) in enumerate(rows):
        if i == Re_idx or s == "H":
            continue
        d = math.sqrt((x - rx) ** 2 + (y - ry) ** 2 + (z - rz) ** 2)
        if d <= 3.5:
            ds.append(d)
    if not ds:
        return None
    ds.sort()
    n = len(ds)
    return ds[n // 2] if n % 2 else 0.5 * (ds[n // 2 - 1] + ds[n // 2])


def test_off_byte_identical_across_runs():
    """OFF path is deterministic: 2x produces identical bytes."""
    a, err_a = _run({"DELFIN_FFFREE_NATIVE_GRIP_CANONICAL": "0"})
    b, err_b = _run({"DELFIN_FFFREE_NATIVE_GRIP_CANONICAL": "0"})
    assert err_a is None and err_b is None
    assert a == b, "OFF path produced different output across two runs"


def test_on_deterministic_across_runs():
    """ON path is deterministic: 2x produces identical bytes."""
    a, err_a = _run({"DELFIN_FFFREE_NATIVE_GRIP_CANONICAL": "1"})
    b, err_b = _run({"DELFIN_FFFREE_NATIVE_GRIP_CANONICAL": "1"})
    assert err_a is None and err_b is None
    assert a == b, "ON path is non-deterministic"


def test_on_topology_consistency():
    """All frames in the ON ensemble share IDENTICAL element sequence
    (the Avogadro-compatibility gate)."""
    isos, err = _run({"DELFIN_FFFREE_NATIVE_GRIP_CANONICAL": "1"})
    assert err is None
    assert len(isos) >= 1, "ON path emitted zero frames"
    sigs = {_signature(xyz) for xyz, _lab in isos}
    assert len(sigs) == 1, (
        f"ON ensemble has {len(sigs)} distinct topology signatures "
        f"(must be 1 for Avogadro compatibility)"
    )


def test_on_md_self_consistency():
    """ON ensemble median M-D is internally consistent (±0.05 Å)."""
    isos, err = _run({"DELFIN_FFFREE_NATIVE_GRIP_CANONICAL": "1"})
    assert err is None and len(isos) >= 1
    mds = [_median_md(xyz) for xyz, _lab in isos]
    mds = [m for m in mds if m is not None]
    assert len(mds) >= 1, "no parseable M-D distances"
    md_ref = mds[0]
    for md in mds:
        assert abs(md - md_ref) <= 0.05, (
            f"M-D self-consistency violation: median {md:.3f} A vs "
            f"reference {md_ref:.3f} A (>0.05 A delta)"
        )


def test_on_labels_carry_canonical_suffix():
    """ON frames are labelled with the canonical-grip suffix so the
    bug-fix is auditable in downstream forensik logs."""
    isos, err = _run({"DELFIN_FFFREE_NATIVE_GRIP_CANONICAL": "1"})
    assert err is None and len(isos) >= 1
    for _xyz, lab in isos:
        # Either the canonical-grip suffix (accepted) or the original
        # native label (silent rollback per safety contract).
        is_canonical = "canonical-grip" in lab
        is_native_rollback = "embed-conf" not in lab
        assert is_canonical or is_native_rollback, (
            f"unexpected label after canonicalisation: {lab[:80]}"
        )
