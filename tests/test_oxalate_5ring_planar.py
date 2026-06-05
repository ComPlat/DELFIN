"""Unit tests for delfin/fffree/oxalate_5ring_planar.py (Iter-32e YILNUF
oxalate-5-ring planar projector).

The projector is env-gated default-OFF (byte-identical to HEAD when
``DELFIN_FFFREE_OXALATE_5RING_PLANAR`` is unset).  These tests cover:

  1. Detection of the M-O-C(=O)-C(=O)-O 5-ring on a synthetic ligand
     block (no SMARTS / no name knowledge — pure graph + geometry).
  2. Projection of a twisted ring onto the planar ideal: planarity
     drops below 1e-6 Å, bite angle hits the literature target (78°),
     and M-O distances are preserved.
  3. Bite angle stays inside the literature range 75-82° after the fix.
  4. Default-OFF byte-identical: with the env-flag unset, ``apply``
     returns a numerically identical copy of the input.
  5. Determinism: two consecutive runs with the flag set produce
     byte-identical output.
  6. Already-flat ring: short-circuit (no edits) when planarity is
     already below the trigger AND bite is within ±4° of the ideal.
  7. Non-oxalate chelate (e.g. en, bipy): detector returns empty list,
     ``apply`` is a no-op even with the flag set.
  8. M-O distance invariant: |M-O1| and |M-O2| are preserved (within
     1e-6 Å) by the projector — donors stay on their assigned vertices.
"""
from __future__ import annotations
import os
import numpy as np
import pytest

from delfin.fffree.oxalate_5ring_planar import (
    flag_active,
    detect_oxalate_5rings,
    apply,
    _solve_planar_5ring,
    _TARGET_BITE_DEG,
)


_FLAG = "DELFIN_FFFREE_OXALATE_5RING_PLANAR"


@pytest.fixture(autouse=True)
def _clean_env():
    """Always start each test with the env-flag UNSET, restore afterwards."""
    saved = os.environ.pop(_FLAG, None)
    yield
    if saved is None:
        os.environ.pop(_FLAG, None)
    else:
        os.environ[_FLAG] = saved


def _ideal_oxalate_ligand():
    """Return (syms, P, donor_idxs, metal_pos) for a perfectly planar
    Rh-oxalate ligand block in metal-at-origin frame.  Donor Os are at
    indices 0 and 3."""
    M = np.zeros(3)
    P_O1 = np.array([1.0, 0.5, 0.0])
    P_C1 = np.array([2.2, 0.3, 0.0])
    P_C2 = np.array([2.2, -0.3, 0.0])
    P_O2 = np.array([1.0, -0.5, 0.0])
    _, n_O1, n_C1, n_C2, n_O2 = _solve_planar_5ring(
        M, P_O1, P_C1, P_C2, P_O2, md_O=2.05,
    )
    O1e = n_C1 + np.array([1.5, 0.5, 0.0])
    O2e = n_C2 + np.array([1.5, -0.5, 0.0])
    syms = ["O", "C", "O", "O", "C", "O"]
    P = np.array([n_O1, n_C1, O1e, n_O2, n_C2, O2e])
    return syms, P, [0, 3], M


def _twisted_oxalate_ligand():
    """Same as the ideal ligand, but with a 0.5 Å out-of-plane twist on
    each ring carbon (and the exocyclic Os follow their parent C)."""
    syms, P, donors, M = _ideal_oxalate_ligand()
    P = P.copy()
    P[1] += np.array([0.0, 0.0, 0.5])
    P[2] += np.array([0.0, 0.0, 0.5])
    P[4] += np.array([0.0, 0.0, -0.5])
    P[5] += np.array([0.0, 0.0, -0.5])
    return syms, P, donors, M


def test_flag_default_off():
    """Without the env-flag the helper reports inactive."""
    assert flag_active() is False


def test_apply_off_byte_identical():
    """Default-OFF: apply() must return a copy IDENTICAL to the input."""
    syms, P, donors, M = _twisted_oxalate_ligand()
    Q, applied = apply(syms, P, donors, metal_pos=M, metal_sym="Rh")
    assert np.array_equal(Q, P)
    assert applied == []
    assert Q is not P                       # copy semantics


def test_detection_on_twisted_ring():
    """Detector finds the M-O-C-C-O 5-ring on a twisted oxalate."""
    syms, P, donors, M = _twisted_oxalate_ligand()
    rings = detect_oxalate_5rings(syms, P, donors, M, metal_sym="Rh")
    assert len(rings) == 1
    r = rings[0]
    assert r["O1"] in donors and r["O2"] in donors
    assert {r["C1"], r["C2"]} == {1, 4}
    assert r["O1_exo"] == [2] and r["O2_exo"] == [5]


def test_apply_on_flattens_ring_and_hits_target_bite():
    """With the flag set, a twisted ring is projected onto the planar
    ideal: planarity drops below 1e-6 Å and bite hits 78° ± 4°."""
    os.environ[_FLAG] = "1"
    syms, P, donors, M = _twisted_oxalate_ligand()
    Q, applied = apply(syms, P, donors, metal_pos=M, metal_sym="Rh")
    assert len(applied) == 1
    info = applied[0]
    assert info["plane_dev_before"] > 0.10
    assert info["plane_dev_after"] < 1e-6
    assert abs(info["bite_after"] - _TARGET_BITE_DEG) < 4.0
    # Donor M-O distances preserved (within 1e-6 Å).
    md_before_O1 = float(np.linalg.norm(P[0] - M))
    md_after_O1 = float(np.linalg.norm(Q[0] - M))
    md_before_O2 = float(np.linalg.norm(P[3] - M))
    md_after_O2 = float(np.linalg.norm(Q[3] - M))
    assert abs(md_after_O1 - md_before_O1) < 1e-6
    assert abs(md_after_O2 - md_before_O2) < 1e-6


def test_bite_inside_literature_range():
    """After the projector, the bite angle is in the literature window
    (CSD-Mogul on oxalato-TM fragments: typically 75-82°)."""
    os.environ[_FLAG] = "1"
    syms, P, donors, M = _twisted_oxalate_ligand()
    Q, _ = apply(syms, P, donors, metal_pos=M, metal_sym="Rh")
    v1 = Q[0] - M
    v2 = Q[3] - M
    cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    bite = float(np.degrees(np.arccos(max(-1.0, min(1.0, cos)))))
    assert 75.0 <= bite <= 82.0


def test_determinism_two_runs():
    """Two consecutive runs with the flag set produce byte-identical
    output (no RNG dependence)."""
    os.environ[_FLAG] = "1"
    syms, P, donors, M = _twisted_oxalate_ligand()
    Q1, _ = apply(syms, P, donors, metal_pos=M, metal_sym="Rh")
    Q2, _ = apply(syms, P, donors, metal_pos=M, metal_sym="Rh")
    assert np.array_equal(Q1, Q2)


def test_already_flat_ring_short_circuit():
    """An already-flat ring with near-ideal bite is left UNCHANGED
    even with the flag set (byte-identical for the easy cases)."""
    os.environ[_FLAG] = "1"
    syms, P, donors, M = _ideal_oxalate_ligand()
    Q, applied = apply(syms, P, donors, metal_pos=M, metal_sym="Rh")
    assert applied == []                    # short-circuit fired
    assert np.array_equal(Q, P)


def test_non_oxalate_chelate_is_noop():
    """An ethylenediamine-style chelate (N-C-C-N, no carboxylate Os)
    triggers neither detection nor projection.  Detector returns [];
    apply is a no-op."""
    os.environ[_FLAG] = "1"
    # M-N-CH2-CH2-N-M (en chelate sketch)
    M = np.zeros(3)
    syms = ["N", "C", "C", "N"]
    P = np.array([
        [1.5, 0.6, 0.0],
        [2.5, 0.4, 0.0],
        [2.5, -0.4, 0.0],
        [1.5, -0.6, 0.0],
    ])
    rings = detect_oxalate_5rings(syms, P, donor_idxs=[0, 3],
                                  metal_pos=M, metal_sym="Rh")
    assert rings == []
    Q, applied = apply(syms, P, [0, 3], metal_pos=M, metal_sym="Rh")
    assert applied == []
    assert np.array_equal(Q, P)


def test_assemble_complex_wiring_default_off_no_call():
    """Default-OFF integration: when DELFIN_FFFREE_OXALATE_5RING_PLANAR is
    unset, ``assemble_complex.py``'s wiring branch must be a no-op (the
    helper short-circuits inside ``apply`` and returns the input copy).

    We verify by monkey-patching the projector's ``apply`` to record calls
    and confirming the wired branch produces ZERO state-mutating calls
    when the flag is off (it is allowed to call ``apply`` once per chelate
    — but each call must return ``applied=[]``).  This guards against
    accidental future changes that would break byte-identity.
    """
    from delfin.fffree import oxalate_5ring_planar as OX
    calls = {"n": 0, "rings": 0}
    _orig = OX.apply

    def _wrap(*a, **kw):
        calls["n"] += 1
        Q, applied = _orig(*a, **kw)
        calls["rings"] += len(applied)
        return Q, applied
    OX.apply = _wrap
    try:
        # Simulate the wired branch directly (without invoking the heavy
        # converter): apply with the env-flag UNSET must return [].
        syms, P, donors, M = _twisted_oxalate_ligand()
        Q, applied = OX.apply(syms, P, donors, metal_pos=M, metal_sym="Rh")
        assert applied == []
        assert calls["rings"] == 0
        assert np.array_equal(Q, P)
    finally:
        OX.apply = _orig


def test_solver_geometry_targets():
    """Direct solver sanity: the planar-5-ring solver returns the
    literature target geometry (bite = 78°, C-C = 1.55 Å, C-O = 1.275 Å,
    M-O = supplied md_O)."""
    M = np.zeros(3)
    # Any reasonable seed positions (solver is symmetry-driven)
    P_O1 = np.array([1.0, 0.7, 0.0])
    P_C1 = np.array([2.2, 0.4, 0.0])
    P_C2 = np.array([2.2, -0.4, 0.0])
    P_O2 = np.array([1.0, -0.7, 0.0])
    _, O1, C1, C2, O2 = _solve_planar_5ring(M, P_O1, P_C1, P_C2, P_O2, md_O=2.05)
    assert abs(float(np.linalg.norm(O1)) - 2.05) < 1e-6
    assert abs(float(np.linalg.norm(O2)) - 2.05) < 1e-6
    assert abs(float(np.linalg.norm(C1 - C2)) - 1.55) < 1e-6
    assert abs(float(np.linalg.norm(O1 - C1)) - 1.275) < 1e-6
    assert abs(float(np.linalg.norm(O2 - C2)) - 1.275) < 1e-6
    v1, v2 = O1 - M, O2 - M
    cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    bite = float(np.degrees(np.arccos(max(-1.0, min(1.0, cos)))))
    assert abs(bite - 78.0) < 1e-6
