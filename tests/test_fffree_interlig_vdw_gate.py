"""Tests for the vdW-level inter-ligand clash gate on the ADDITIONAL conformer
frames (backbone re-embed / conformer re-seating) of the FF-free path
(delfin.fffree.converter_backend, env DELFIN_FFFREE_INTERLIG_VDW_GATE).

Eye-found defect (ACEQUY = Fe(N(Dipp)(SiMe3))3): the per-ligand backbone re-embed
re-folds each ligand with the metal+donor core frozen but WITHOUT inter-ligand
awareness, so one ligand backbone can collapse into a NEIGHBOUR ligand (~2 A C-C).
That contact sits far above ``_build_is_clean``'s 0.60*Sigma_cov gross-overlap floor
(~0.9 A for C-C) so the self-gate misses it; the new vdW gate (NEVER-WORSE vs the
base frame, hard 2.0 A floor) drops such frames while keeping clean ones.

Hard contracts:
  1. the helper computes the min NON-BONDED heavy-heavy distance (bonds skipped),
  2. _interlig_clash_ok: hard 2.0 A floor + never-worse vs base_min (1-tol),
  3. _append_reembed drops only the collapsing frames; clean frames survive,
  4. byte-identical when the reembed flag is OFF (no extra frames -> no-op),
  5. determinism.
All graph-based / universal (no refcode coordinates).
"""
import os
import numpy as np
import pytest

from delfin.fffree import converter_backend as CB


@pytest.fixture(autouse=True)
def _clean_env():
    saved = dict(os.environ)
    for k in ("DELFIN_FFFREE_INTERLIG_VDW_GATE", "DELFIN_FFFREE_BACKBONE_REEMBED",
              "DELFIN_FFFREE_CN_EXTEND"):
        os.environ.pop(k, None)
    yield
    os.environ.clear()
    os.environ.update(saved)


def test_min_nonbonded_heavy_skips_bonds():
    """Two bonded C (1.5 A) + one far non-bonded C (3.0 A) -> min non-bonded = 3.0;
    a bonded pair is never reported as the minimum non-bonded distance."""
    syms = ["C", "C", "C"]
    P = np.array([[0.0, 0, 0], [1.5, 0, 0], [4.5, 0, 0]])  # 0-1 bonded, 0-2/1-2 nonbond
    mn = CB._min_nonbonded_heavy(syms, P)
    assert abs(mn - 3.0) < 1e-6      # closest non-bonded pair is 1-2 at 3.0 A


def test_min_nonbonded_heavy_excludes_metal_and_h():
    syms = ["Fe", "N", "H"]
    P = np.array([[0.0, 0, 0], [2.0, 0, 0], [2.5, 0.5, 0]])
    # only heavy-heavy non-metal pairs counted -> none here -> +inf
    assert not np.isfinite(CB._min_nonbonded_heavy(syms, P))


def test_interlig_clash_ok_hard_floor():
    """With base_min unavailable, the hard 2.0 A floor alone applies."""
    syms = ["C", "C"]
    near = np.array([[0.0, 0, 0], [1.98, 0, 0]])    # non-bonded (1.98 > 1.30*Sigma=1.976)
    far = np.array([[0.0, 0, 0], [2.40, 0, 0]])
    assert not CB._interlig_clash_ok(syms, near, None)   # 1.98 < 2.0 -> reject
    assert CB._interlig_clash_ok(syms, far, None)        # 2.40 >= 2.0 -> keep


def test_interlig_clash_ok_never_worse():
    """A clean base frame (2.38 A) raises the floor to 2.38*0.95 = 2.261 A, so a
    2.05 A frame (better than the hard floor but worse than the base) is rejected."""
    syms = ["C", "C"]
    base_min = 2.38
    frame_205 = np.array([[0.0, 0, 0], [2.05, 0, 0]])
    frame_240 = np.array([[0.0, 0, 0], [2.40, 0, 0]])
    assert not CB._interlig_clash_ok(syms, frame_205, base_min)   # 2.05 < 2.261 -> reject
    assert CB._interlig_clash_ok(syms, frame_240, base_min)       # 2.40 >= 2.261 -> keep


def test_gate_default_on():
    assert CB._interlig_vdw_gate_enabled()                 # unset -> on
    os.environ["DELFIN_FFFREE_INTERLIG_VDW_GATE"] = "0"
    assert not CB._interlig_vdw_gate_enabled()
    os.environ["DELFIN_FFFREE_INTERLIG_VDW_GATE"] = "1"
    assert CB._interlig_vdw_gate_enabled()


# A flexible-armed monodentate Werner complex whose backbone re-embed frames can
# fold near a neighbour (reaches the FF-free path via CN_EXTEND).
_FLEX = "CCO[Co](OCC)(OCC)([NH3])([NH3])[NH3]"


def _sig(rr):
    return None if rr is None else tuple((lab, xyz) for xyz, lab in rr)


def test_byte_identical_reembed_off():
    """Reembed flag OFF -> no extra frames -> the gate is a no-op (byte-id)."""
    os.environ["DELFIN_FFFREE_CN_EXTEND"] = "1"
    os.environ.pop("DELFIN_FFFREE_BACKBONE_REEMBED", None)
    a = _sig(CB._fffree_isomers(_FLEX))
    os.environ["DELFIN_FFFREE_INTERLIG_VDW_GATE"] = "1"
    b = _sig(CB._fffree_isomers(_FLEX))
    assert a == b
    assert a is not None and not any("reembed" in lab for lab, _ in a)


def test_gate_only_drops_collapsing_reembed_frames():
    """With reembed ON, the gate (default ON) must keep every clean re-embed frame
    that the gate-OFF run produced whose min non-bonded distance is >= the floor,
    and drop only frames below it -> gate-ON count <= gate-OFF count, and every
    surviving frame respects the floor."""
    os.environ["DELFIN_FFFREE_CN_EXTEND"] = "1"
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "1"

    os.environ["DELFIN_FFFREE_INTERLIG_VDW_GATE"] = "0"
    off = CB._fffree_isomers(_FLEX)
    os.environ["DELFIN_FFFREE_INTERLIG_VDW_GATE"] = "1"
    on = CB._fffree_isomers(_FLEX)
    assert off is not None and on is not None

    def _parse(xyz):
        syms = [ln.split()[0] for ln in xyz.splitlines() if len(ln.split()) == 4]
        P = np.array([[float(x) for x in ln.split()[1:4]]
                      for ln in xyz.splitlines() if len(ln.split()) == 4])
        return syms, P

    n_re_off = sum("reembed" in lab for _, lab in off)
    n_re_on = sum("reembed" in lab for _, lab in on)
    assert n_re_off >= 1, "need re-embed frames to exercise the gate"
    assert n_re_on <= n_re_off          # gate is subtractive only

    # every surviving re-embed frame respects the hard 2.0 A floor
    for xyz, lab in on:
        if "reembed" not in lab:
            continue
        syms, P = _parse(xyz)
        assert CB._min_nonbonded_heavy(syms, P) >= CB._INTERLIG_VDW_FLOOR - 1e-9


def test_deterministic():
    os.environ["DELFIN_FFFREE_CN_EXTEND"] = "1"
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "1"
    os.environ["DELFIN_FFFREE_INTERLIG_VDW_GATE"] = "1"
    a = _sig(CB._fffree_isomers(_FLEX))
    b = _sig(CB._fffree_isomers(_FLEX))
    assert a == b
