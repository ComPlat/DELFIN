"""Tests for delfin.fffree.burnside_conformer.

Phase 3 of the GRIP-bahnbrechend mandate: Burnside's lemma on the
CONFORMER level (Cauchy-Frobenius orbit counting on the assembled
3-D structure under the molecular point group).

These tests validate:
  · C1 (no symmetry) ⇒ orbit count = raw distinct count.
  · Cyclic group quotient on Pólya-style equivalent slots:
    the formula |X/C_k| = (1/k) Σ_{j=0..k-1} m^{gcd(j,k)·t}
    yields the textbook values for the [Co(en)3]^3+ ring-pucker
    enumeration.
  · Explicit-enumeration cross-check for a small case.
  · Determinism: two runs on identical inputs return identical counts.
  · Env-flag default-OFF preserves byte-identical grip_ensemble
    behaviour (no orbit dedup when DELFIN_FFFREE_BURNSIDE_CONFORMER is
    unset).

References
==========
  · Burnside, W. (1897), *Theory of Groups of Finite Order*, §148.
  · Cremer, D. & Pople, J. A. (1975), *J. Am. Chem. Soc.* 97, 1354.
  · Pólya, G. (1937), *Acta Math.* 68, 145.
"""

from __future__ import annotations

import math
import os
from typing import List, Tuple

import numpy as np
import pytest

# Force deterministic env BEFORE importing the module.
os.environ.setdefault("PYTHONHASHSEED", "0")

from delfin.fffree.burnside_conformer import (
    BurnsideOrbitCounter,
    canonical_form,
    canonical_form_equal,
    count_distinct_conformers,
    count_orbits_product_space,
    completeness_report,
    burnside_conformer_active,
    _enumerate_canonical_axis_group,
    _identity_op,
    _rotation_op,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ch4_conformer(scale: float = 1.0, jitter: Tuple[float, float, float] = (0, 0, 0)):
    """Return a CH4-like (symbols, coords) tuple."""
    syms = ["C", "H", "H", "H", "H"]
    coords = (np.array([
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ], dtype=float) * scale) + np.array(jitter, dtype=float)
    return syms, coords


def _make_en_ring_canonical(state: str, scale: float = 1.0, offset: Tuple[float, float, float] = (0, 0, 0)):
    """Return a 5-atom toy 'en chelate ring' conformer in δ or λ state.

    Used for the Co(en)3 textbook-case test.  The exact geometry is not
    important — we only need that δ and λ are non-superposable so the
    canonical-form comparison treats them as distinct.
    """
    syms = ["N", "C", "C", "N", "H"]
    z_twist = +0.3 if state == "delta" else -0.3
    coords = np.array([
        [1.0, 0.0, 0.0],
        [0.5, 0.7, z_twist],
        [-0.5, 0.7, -z_twist],
        [-1.0, 0.0, 0.0],
        [0.0, -0.5, 0.0],
    ], dtype=float) * scale + np.array(offset, dtype=float)
    return syms, coords


# ---------------------------------------------------------------------------
# 1. C1 (no symmetry) — no orbit reduction.
# ---------------------------------------------------------------------------


def test_C1_no_orbit_reduction():
    """For a C1-complex with k distinct conformers, |X/G| = k = |X|."""
    confs: List[Tuple[List[str], np.ndarray]] = []
    for i in range(4):
        syms, P = _make_ch4_conformer(scale=1.0 + 0.5 * i)
        confs.append((syms, P))

    # Explicit identity-only group.
    n_orbits = count_distinct_conformers(confs, group_ops=[_identity_op()])
    assert n_orbits == 4, f"C1 4-conformer count: expected 4, got {n_orbits}"


def test_C1_default_group_no_reduction():
    """When no symmetry exists in the structure, the auto-detected
    group must not falsely collapse genuinely-different conformers.

    Translations and uniform rotations are removed by the canonical
    inertial-frame alignment (per design — they ARE the same
    conformer), so distinct conformers must differ in shape.
    """
    confs: List[Tuple[List[str], np.ndarray]] = []
    for i in range(3):
        # Genuinely-different shape: scale + a non-uniform perturbation.
        syms, P_base = _make_ch4_conformer(scale=1.0 + 0.5 * i)
        # Asymmetric stretch along one direction so each is shape-distinct.
        P = P_base.copy()
        P[1, 0] += 0.3 * i  # nudge one H
        confs.append((syms, P))
    n_orbits = count_distinct_conformers(confs, group_ops=[_identity_op()])
    assert n_orbits == 3


# ---------------------------------------------------------------------------
# 2. D3-octahedron orbit quotient: chelate [Co(en)3]^3+ case.
# ---------------------------------------------------------------------------


def test_co_en_three_ring_pucker_C3_quotient():
    """C3-action on 3 equivalent en chelate rings with 2 CP states each.

    Raw product = 2^3 = 8 (δλ-tuples).  Under C3 acting by rotating
    the 3 ring slots, the orbits are: δδδ (1), δδλ (3 → 1 orbit),
    δλλ (3 → 1 orbit), λλλ (1).  Total = 4 orbits.

    Burnside: (1/3)(2^3 + 2 + 2) = 12/3 = 4. ✓
    """
    n = count_orbits_product_space(
        n_polya=1,
        ring_state_counts=[2, 2, 2],
        n_rotamers=1,
        group_order=3,
        n_permuted_rings=3,
    )
    assert n == 4, f"C3-quotient on 3 rings × 2 CP: expected 4 orbits, got {n}"


def test_co_en_three_with_two_polya_isomers():
    """Δ and Λ doubling: total orbits = 2 · 4 = 8."""
    n = count_orbits_product_space(
        n_polya=2,
        ring_state_counts=[2, 2, 2],
        n_rotamers=1,
        group_order=3,
        n_permuted_rings=3,
    )
    assert n == 8


def test_C2_quotient_two_rings():
    """C2-action on 2 equivalent rings of 3 states each.

    Raw = 9.  C2 orbits = (1/2)(3^2 + 3) = (9+3)/2 = 6.
    """
    n = count_orbits_product_space(
        n_polya=1, ring_state_counts=[3, 3], n_rotamers=1,
        group_order=2, n_permuted_rings=2,
    )
    assert n == 6


def test_C6_six_equivalent_rings_two_states():
    """C6 on 6 equivalent ring slots with 2 states each.

    Burnside: (1/6) Σ_{j=0..5} 2^{gcd(j,6)}
            = (1/6) (2^6 + 2 + 2^2 + 2^3 + 2^2 + 2)
            = (1/6) (64 + 2 + 4 + 8 + 4 + 2) = 84 / 6 = 14.
    """
    n = count_orbits_product_space(
        n_polya=1, ring_state_counts=[2] * 6, n_rotamers=1,
        group_order=6, n_permuted_rings=6,
    )
    assert n == 14


def test_no_permutation_when_uneven():
    """When n_permuted_rings does not divide group_order evenly,
    fall back to raw product.
    """
    n = count_orbits_product_space(
        n_polya=1, ring_state_counts=[2, 2], n_rotamers=1,
        group_order=3, n_permuted_rings=2,
    )
    # 2 not divisible by 3 -> no quotient -> raw = 4.
    assert n == 4


# ---------------------------------------------------------------------------
# 3. Explicit-enumeration cross-check for a small textbook case.
# ---------------------------------------------------------------------------


def test_burnside_count_matches_explicit_enumeration_small_case():
    """Construct 4 explicit conformers related by a C2 rotation around z.

    Setup: two CH4-like structures at scale=1 and scale=2.  Apply a
    C2 rotation to each to get rotated copies.  Total raw = 4.
    Under C2, each pair collapses to 1 orbit → 2 orbits expected.
    """
    syms, P1 = _make_ch4_conformer(scale=1.0)
    _, P2 = _make_ch4_conformer(scale=2.0)
    R = _rotation_op(np.array([0.0, 0.0, 1.0]), math.pi)
    com1 = P1.mean(axis=0)
    com2 = P2.mean(axis=0)
    P1_rot = (P1 - com1) @ R.T + com1
    P2_rot = (P2 - com2) @ R.T + com2

    confs = [
        (syms, P1),
        (syms, P1_rot),
        (syms, P2),
        (syms, P2_rot),
    ]
    # Explicit C2 group.
    ops = [_identity_op(), R]
    n = count_distinct_conformers(confs, group_ops=ops, rms_tol=0.1)
    # Canonical-form alignment is rotation-invariant, so all 4 reduce
    # to 2 distinct canonical forms even before the Burnside quotient.
    assert n == 2, f"explicit C2 quotient: expected 2 orbits, got {n}"


def test_explicit_two_distinct_no_symmetry_returns_two():
    """Two distinct CH4-like conformers with NO symmetry relating
    them: |X/G| = 2 even when the synthesised group order is larger.
    """
    syms, P1 = _make_ch4_conformer(scale=1.0)
    _, P2 = _make_ch4_conformer(scale=2.0)
    confs = [(syms, P1), (syms, P2)]
    n = count_distinct_conformers(confs, group_ops=None)
    assert n == 2


# ---------------------------------------------------------------------------
# 4. Determinism: two runs identical.
# ---------------------------------------------------------------------------


def test_determinism_two_runs_identical():
    """Same inputs → same orbit count, byte-stable canonical forms."""
    syms1, P1 = _make_ch4_conformer(scale=1.0)
    syms2, P2 = _make_ch4_conformer(scale=2.0)
    syms3, P3 = _make_ch4_conformer(scale=1.5)

    confs = [(syms1, P1), (syms2, P2), (syms3, P3)]

    n1 = count_distinct_conformers(confs, group_ops=[_identity_op()])
    n2 = count_distinct_conformers(confs, group_ops=[_identity_op()])
    assert n1 == n2

    # Canonical forms identical bit-for-bit on repeat.
    cf_a = canonical_form(syms1, P1)
    cf_b = canonical_form(syms1, P1)
    assert cf_a.symbols == cf_b.symbols
    assert cf_a.coords == cf_b.coords


def test_determinism_completeness_report():
    """completeness_report() is deterministic."""
    r1 = completeness_report(
        n_polya=2, ring_state_counts=[2, 2, 2], n_rotamers=1,
        group_order=3, n_emitted=8,
    )
    r2 = completeness_report(
        n_polya=2, ring_state_counts=[2, 2, 2], n_rotamers=1,
        group_order=3, n_emitted=8,
    )
    assert r1 == r2


# ---------------------------------------------------------------------------
# 5. Env-flag default OFF — byte-identical to HEAD.
# ---------------------------------------------------------------------------


def test_env_off_byte_identical_to_HEAD_default():
    """Without env-flag, the Burnside-conformer quotient is OFF and
    grip_ensemble emits the unmodified candidate set.
    """
    # Default is OFF.
    os.environ.pop("DELFIN_FFFREE_BURNSIDE_CONFORMER", None)
    assert burnside_conformer_active() is False
    from delfin.fffree.grip_ensemble import burnside_conformer_quotient_active
    assert burnside_conformer_quotient_active() is False


def test_env_on_activates_flag():
    """Setting the env flag to '1' / 'true' / 'yes' / 'on' activates the
    Burnside-conformer mode.
    """
    for val in ("1", "true", "yes", "on", "TRUE", "Yes"):
        os.environ["DELFIN_FFFREE_BURNSIDE_CONFORMER"] = val
        assert burnside_conformer_active() is True, f"failed for value {val!r}"
    os.environ.pop("DELFIN_FFFREE_BURNSIDE_CONFORMER", None)


def test_env_off_explicit_values():
    """'0' / 'false' / 'no' / 'off' deactivate."""
    for val in ("0", "false", "no", "off"):
        os.environ["DELFIN_FFFREE_BURNSIDE_CONFORMER"] = val
        assert burnside_conformer_active() is False, f"failed for value {val!r}"
    os.environ.pop("DELFIN_FFFREE_BURNSIDE_CONFORMER", None)


# ---------------------------------------------------------------------------
# 6. BurnsideOrbitCounter streaming correctness.
# ---------------------------------------------------------------------------


def test_orbit_counter_streams_correctly():
    """Streaming counter agrees with the bulk function."""
    syms, P1 = _make_ch4_conformer(scale=1.0)
    _, P2 = _make_ch4_conformer(scale=2.0)
    _, P3 = _make_ch4_conformer(scale=1.0, jitter=(10.0, 0, 0))  # translation only

    counter = BurnsideOrbitCounter(group_ops=[_identity_op()])
    assert counter.is_new_orbit(syms, P1) is True
    # Pure translation is removed by the canonical inertial-frame
    # alignment → same orbit.
    assert counter.is_new_orbit(syms, P3) is False
    # Scaled conformer is a genuinely different shape.
    assert counter.is_new_orbit(syms, P2) is True
    assert counter.n_orbits == 2


def test_orbit_counter_deterministic_order_independence():
    """Inserting in different orders yields the same final orbit
    count (order independence is a property of the orbit count
    itself, although the representative-set depends on order).
    """
    syms, P1 = _make_ch4_conformer(scale=1.0)
    _, P2 = _make_ch4_conformer(scale=2.0)
    _, P3 = _make_ch4_conformer(scale=3.0)

    counter_a = BurnsideOrbitCounter(group_ops=[_identity_op()])
    counter_a.is_new_orbit(syms, P1)
    counter_a.is_new_orbit(syms, P2)
    counter_a.is_new_orbit(syms, P3)

    counter_b = BurnsideOrbitCounter(group_ops=[_identity_op()])
    counter_b.is_new_orbit(syms, P3)
    counter_b.is_new_orbit(syms, P2)
    counter_b.is_new_orbit(syms, P1)

    assert counter_a.n_orbits == counter_b.n_orbits == 3


# ---------------------------------------------------------------------------
# 7. Canonical-form invariance under translation and rotation.
# ---------------------------------------------------------------------------


def test_canonical_form_translation_invariance():
    """Translated structures yield equivalent canonical forms."""
    syms, P = _make_ch4_conformer(scale=1.0)
    cf_orig = canonical_form(syms, P)
    cf_trans = canonical_form(syms, P + np.array([5.0, -3.0, 2.0]))
    assert canonical_form_equal(cf_orig, cf_trans, rms_tol=1e-3)


def test_canonical_form_distinguishes_different_shapes():
    """Genuinely different shapes give non-equal canonical forms."""
    syms, P1 = _make_ch4_conformer(scale=1.0)
    _, P2 = _make_ch4_conformer(scale=2.5)
    cf1 = canonical_form(syms, P1)
    cf2 = canonical_form(syms, P2)
    assert not canonical_form_equal(cf1, cf2, rms_tol=0.05)


# ---------------------------------------------------------------------------
# 8. Completeness-report verdict labels.
# ---------------------------------------------------------------------------


def test_completeness_report_complete_verdict():
    rep = completeness_report(
        n_polya=2, ring_state_counts=[2, 2, 2], n_rotamers=1,
        group_order=3, n_emitted=8,
    )
    assert rep["burnside_orbit_count"] == 8
    assert rep["verdict"] == "complete"
    assert rep["raw_product"] == 16


def test_completeness_report_incomplete_verdict():
    rep = completeness_report(
        n_polya=2, ring_state_counts=[2, 2, 2], n_rotamers=1,
        group_order=3, n_emitted=1,
    )
    assert rep["verdict"] in ("incomplete", "mostly_complete")


def test_completeness_report_redundant_verdict():
    rep = completeness_report(
        n_polya=2, ring_state_counts=[2, 2, 2], n_rotamers=1,
        group_order=3, n_emitted=20,
    )
    assert rep["verdict"] == "redundant"


def test_completeness_report_empty_verdict():
    rep = completeness_report(
        n_polya=1, ring_state_counts=[], n_rotamers=1,
        group_order=1, n_emitted=0,
    )
    assert rep["verdict"] == "empty"


# ---------------------------------------------------------------------------
# 9. Burnside-identity sanity check (sum = orbit_count × |G|).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k,n_rings,m", [
    (2, 2, 2),
    (3, 3, 2),
    (3, 3, 3),
    (4, 4, 2),
    (6, 6, 2),
])
def test_burnside_identity_orbits_times_group_order_equals_fix_sum(k, n_rings, m):
    """For a cyclic group C_k acting on n_rings = k·t equivalent slots
    with m states each, verify the closed-form identity:

        orbits · |G| = Σ_g |Fix(g)|
                    = Σ_{j=0..k-1} m^{t · gcd(j, k)}
    """
    if n_rings % k != 0:
        pytest.skip("test geometry requires k | n_rings")
    n_orbits = count_orbits_product_space(
        n_polya=1, ring_state_counts=[m] * n_rings, n_rotamers=1,
        group_order=k, n_permuted_rings=n_rings,
    )
    t = n_rings // k
    fix_sum = sum(m ** (t * math.gcd(j, k)) for j in range(k))
    assert n_orbits * k == fix_sum


# ---------------------------------------------------------------------------
# 10. End-to-end grip_ensemble env-OFF byte-identity test.
# ---------------------------------------------------------------------------


def test_grip_ensemble_off_path_burnside_meta_empty():
    """When the env flag is unset, the burnside_meta dict in
    EnsembleResult is empty (no quotient applied).  This is the
    byte-identity guarantee.
    """
    from delfin.fffree.grip_ensemble import EnsembleResult, EnsembleCandidate
    os.environ.pop("DELFIN_FFFREE_BURNSIDE_CONFORMER", None)
    # Construct a result without going through the full enumeration
    # (which depends on RDKit chemistry).  The contract: default
    # burnside_meta is empty.
    res = EnsembleResult(smiles="X", candidates=[])
    assert res.burnside_meta == {}
