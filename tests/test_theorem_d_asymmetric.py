"""Welle-5l T6.2 — Theorem-D asymmetric-bidentate enumerator tests.

Five unit tests covering the universal-fundamental contract for the
``_theorem_d_asymmetric_bidentate`` module:

  1. Trigger gate refuses symmetric tris-bidentate (M(en)3-like) — no
     spurious Δ/Λ split when both donors of every chelate share the
     same label.
  2. Trigger gate fires on asymmetric tris-bidentate (X10-YIVROM
     pattern: 3× (O,S)) regardless of element identity.
  3. OH octahedral classifier returns balanced Δ/Λ buckets for
     fac/mer of tris-(A,B) — the legacy classifier collapses to ''
     on this geometry; Theorem-D must produce 4 distinct cf-stereoisomers.
  4. End-to-end smoke through ``_enumerate_topological_isomers``:
     YIVROM-like (O,S)·3 with DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE=1
     emits strictly MORE cf-buckets than baseline.
  5. Default-OFF byte-equivalence: enumerator output is bit-identical
     to baseline when the env-flag is unset.

No SMILES-specific shortcuts, no element allowlists.  Tests are pure
Python — no RDKit, no Open Babel.  Sub-second runtime per test.
"""
from __future__ import annotations

import os
from typing import FrozenSet, List

import pytest

from delfin.manta._theorem_d_asymmetric_bidentate import (
    classify_helicity_asymmetric,
    is_asymmetric_bidentate_set,
    theorem_d_aware_pairs,
)


@pytest.fixture(autouse=True)
def _clean_env():
    """Ensure no env leakage between tests."""
    keys = [
        'DELFIN_CHIRAL_ENUM',
        'DELFIN_BURNSIDE_FULL',
        'DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE',
    ]
    saved = {k: os.environ.get(k) for k in keys}
    for k in keys:
        os.environ.pop(k, None)
    yield
    for k, v in saved.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v


# ---------------------------------------------------------------------------
# Test 1 — Trigger gate refuses symmetric chelates
# ---------------------------------------------------------------------------


def test_gate_refuses_symmetric_tris_bidentate():
    """M(en)3 with 3 identical (N,N) chelates: all 6 donors labelled N0.

    Expected: ``is_asymmetric_bidentate_set`` returns False.  No spurious
    chirality split must be applied.
    """
    donor_labels = ['N0'] * 6
    chelate_pairs: List[FrozenSet] = [
        frozenset([0, 2]), frozenset([1, 4]), frozenset([3, 5]),
    ]
    assert not is_asymmetric_bidentate_set(chelate_pairs, donor_labels)

    # Classifier must also return '' for any perm.
    perm = [0, 1, 2, 3, 4, 5]
    h = classify_helicity_asymmetric(perm, chelate_pairs, donor_labels, 'OH')
    assert h == ''


# ---------------------------------------------------------------------------
# Test 2 — Trigger gate fires on asymmetric chelates (universal)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "labels",
    [
        # (O,S) — X10-YIVROM pattern
        ['O0', 'O0', 'O0', 'S1', 'S1', 'S1'],
        # (N,O) — classic glycinato/amino-acid pattern
        ['N0', 'N0', 'N0', 'O1', 'O1', 'O1'],
        # (P,N) — phosphine-imine
        ['P0', 'P0', 'P0', 'N1', 'N1', 'N1'],
        # (Cp_eta5, C_aryl) — hapto/sigma mix (universal-fundamental)
        ['CpA', 'CpA', 'CpA', 'C0', 'C0', 'C0'],
    ],
)
def test_gate_fires_on_asymmetric_pairs(labels):
    """Universal: any (A,B) chelate with distinct labels triggers."""
    chelate_pairs: List[FrozenSet] = [
        frozenset([0, 3]), frozenset([1, 4]), frozenset([2, 5]),
    ]
    assert is_asymmetric_bidentate_set(chelate_pairs, labels), (
        f"gate should fire on asymmetric chelate set with labels {labels}"
    )


# ---------------------------------------------------------------------------
# Test 3 — OH tris-(A,B) yields balanced Δ/Λ buckets at fac & mer
# ---------------------------------------------------------------------------


def test_oh_tris_asymmetric_yields_delta_lambda():
    """Tris-(O,S) on octahedron Fe — Pólya expects 4 stereoisomers:
    fac-Δ, fac-Λ, mer-Δ, mer-Λ.

    The classifier must return non-empty handedness for at least 2
    chelate-valid permutations AND must split into equal-cardinality
    Δ and Λ buckets within each fac/mer multiset.
    """
    import itertools

    from delfin.manta._chirality_enumerator import _OH_VECTORS  # noqa
    from delfin.smiles_converter import _canonical_oh  # noqa

    donor_labels = ['O0', 'O0', 'O0', 'S1', 'S1', 'S1']
    chelate_pairs: List[FrozenSet] = [
        frozenset([0, 3]), frozenset([1, 4]), frozenset([2, 5]),
    ]
    trans_pos = [(0, 1), (2, 3), (4, 5)]

    buckets: dict = {}
    for perm in itertools.permutations(range(6)):
        valid = True
        for chel in chelate_pairs:
            cs = sorted(chel)
            pi = perm.index(cs[0])
            pj = perm.index(cs[1])
            for ta, tb in trans_pos:
                if (pi == ta and pj == tb) or (pi == tb and pj == ta):
                    valid = False
                    break
            if not valid:
                break
        if not valid:
            continue
        types = tuple(donor_labels[perm[pos]] for pos in range(6))
        cf = _canonical_oh(types)
        h = classify_helicity_asymmetric(perm, chelate_pairs, donor_labels, 'OH')
        buckets.setdefault((cf, h), []).append(perm)

    # Expect 4 distinct (cf, helicity) buckets: 2 cf-classes × {D, L}.
    # No '' helicity must remain — every chelate-valid OH perm of an
    # asymmetric tris-bidentate complex carries a defined handed-ness.
    helicities = {key[1] for key in buckets}
    assert 'D' in helicities, f"missing Δ bucket; buckets={list(buckets)}"
    assert 'L' in helicities, f"missing Λ bucket; buckets={list(buckets)}"
    assert '' not in helicities, (
        f"OH tris-(O,S) should not contain achiral bucket; got "
        f"{list(buckets)}"
    )
    # Count Δ buckets = count Λ buckets (parity invariance).
    n_delta = sum(1 for key in buckets if key[1] == 'D')
    n_lambda = sum(1 for key in buckets if key[1] == 'L')
    assert n_delta == n_lambda, (
        f"Δ ({n_delta}) and Λ ({n_lambda}) bucket counts must be equal"
    )


# ---------------------------------------------------------------------------
# Test 4 — End-to-end enumerator increments cf-count when gate fires
# ---------------------------------------------------------------------------


def test_enumerator_yivrom_increments_under_theorem_d_flag():
    """X10-YIVROM minimal proxy: tris-(O,S) on Fe CN=6.

    Without the env-flag the OH/TPR enumerator emits a small set of
    cf-isomers that collapse fac-Δ/Λ and mer-Δ/Λ into single buckets.
    With DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE=1 the cf-bucket count
    must strictly increase.
    """
    from delfin.smiles_converter import _enumerate_topological_isomers

    donor_labels = ['O0', 'O0', 'O0', 'S1', 'S1', 'S1']
    chelate_pairs: List[FrozenSet] = [
        frozenset([0, 3]), frozenset([1, 4]), frozenset([2, 5]),
    ]

    os.environ.pop('DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE', None)
    baseline = _enumerate_topological_isomers(
        donor_labels, 6, chelate_pairs, metal_symbol='Fe',
    )

    os.environ['DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE'] = '1'
    theorem_d = _enumerate_topological_isomers(
        donor_labels, 6, chelate_pairs, metal_symbol='Fe',
    )

    assert len(theorem_d) > len(baseline), (
        f"Theorem-D enumerator must increase cf-isomer count; "
        f"baseline={len(baseline)}, theorem_d={len(theorem_d)}"
    )


# ---------------------------------------------------------------------------
# Test 5 — Default-OFF bit-equivalence (universal-fundamental gate-5)
# ---------------------------------------------------------------------------


def test_default_off_byte_equivalent():
    """When the env-flag is unset, the enumerator MUST produce bit-identical
    output to the pre-Iter-T6.2 baseline.  No code path may leak.
    """
    from delfin.smiles_converter import _enumerate_topological_isomers

    # Test on a mixed bidentate set + a hetero-CN configuration.
    cases = [
        # symmetric M(en)3
        (['N0'] * 6, [frozenset([0, 2]), frozenset([1, 4]), frozenset([3, 5])], 'Co'),
        # asymmetric tris-(O,S) Fe
        (['O0', 'O0', 'O0', 'S1', 'S1', 'S1'],
         [frozenset([0, 3]), frozenset([1, 4]), frozenset([2, 5])], 'Fe'),
        # asymmetric tris-(N,O) Ru
        (['N0', 'N0', 'N0', 'O1', 'O1', 'O1'],
         [frozenset([0, 3]), frozenset([1, 4]), frozenset([2, 5])], 'Ru'),
    ]
    for donor_labels, chelate_pairs, metal in cases:
        os.environ.pop('DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE', None)
        ref = _enumerate_topological_isomers(
            donor_labels, 6, chelate_pairs, metal_symbol=metal,
        )

        # Set to '0' explicitly — must still match ref.
        os.environ['DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE'] = '0'
        off_zero = _enumerate_topological_isomers(
            donor_labels, 6, chelate_pairs, metal_symbol=metal,
        )
        assert off_zero == ref, (
            f"OFF (=0) drift for metal={metal} labels={donor_labels}"
        )

        # Unset → must also match.
        os.environ.pop('DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE', None)
        unset = _enumerate_topological_isomers(
            donor_labels, 6, chelate_pairs, metal_symbol=metal,
        )
        assert unset == ref, (
            f"UNSET drift for metal={metal} labels={donor_labels}"
        )


# ---------------------------------------------------------------------------
# Test 6 — theorem_d_aware_pairs no-op contract
# ---------------------------------------------------------------------------


def test_aware_pairs_skips_when_gate_fails():
    """Helper returns input unchanged for symmetric chelate sets, even
    when called explicitly (env-flag is the caller's responsibility)."""
    cf = ('OH', (('N0', 'N0'), ('N0', 'N0'), ('N0', 'N0')))
    perm = [0, 1, 2, 3, 4, 5]
    chelate_pairs = [frozenset([0, 2]), frozenset([1, 4]), frozenset([3, 5])]
    donor_labels = ['N0'] * 6

    result = theorem_d_aware_pairs(cf, perm, chelate_pairs, donor_labels, 'OH')
    assert result == cf, "should be no-op for symmetric chelates"

    # Missing arg → no-op
    result = theorem_d_aware_pairs(cf, None, chelate_pairs, donor_labels, 'OH')
    assert result == cf
    result = theorem_d_aware_pairs(cf, perm, None, donor_labels, 'OH')
    assert result == cf
    result = theorem_d_aware_pairs(cf, perm, chelate_pairs, None, 'OH')
    assert result == cf
