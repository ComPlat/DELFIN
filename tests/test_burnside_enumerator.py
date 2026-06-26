"""Pólya/Burnside enumerator correctness audit (Iter-9, 2026-05-13).

These tests pin down two completeness properties of
``delfin._burnside_groups``:

1. Every polyhedron registered in ``_TOPO_GEOMETRY_VECTORS`` (the
   enumerator's geometry dictionary in ``smiles_converter.py``) MUST
   also have a vertex set in ``_burnside_groups._GEO``.  Without it,
   :func:`get_groups` returns the identity-only fall-back and the
   ``DELFIN_BURNSIDE_FULL=1`` dedup path is a silent no-op for that
   geometry — exactly the gap the env-flag was introduced to close.
2. The constructed point-group action on vertex labels has the
   expected proper-rotation order |G+| and full-improper order
   |G_full| for every supported polyhedron, and Burnside-lemma
   orbit counts at k=1..4 colour symbols match the orbit count
   obtained by deduplicating exhaustively via
   :func:`burnside_canonical_key`.  This is a complete invariant
   test: any future change to the geometry vectors or the symmetry
   builder that drops/adds a symmetry element will be caught.

The tests are self-contained (pure Python, no SciPy/RDKit) and run
in well under a second per polyhedron.  No SMILES-specific
shortcuts are used.
"""

from __future__ import annotations

import itertools
from typing import Dict, Tuple

import pytest

from delfin._burnside_groups import (
    _GEO,
    burnside_canonical_key,
    get_groups,
    group_size,
)


# Expected effective-action group orders on vertex labels.
# Reference: standard inorganic chemistry point-group tables.
# When a polyhedron has σh-fixed in-plane vertices (LIN/TP/SQ), the
# effective action on labels is |G|/2 because σh acts as identity on
# those labels — the build_groups discovers this automatically.
EXPECTED_ORDERS: Dict[str, Tuple[int, int]] = {
    # (chiral / proper-rotation order, achiral / full order on labels)
    'LIN':   (2, 2),     # D∞h on 2 colinear vertices
    'TP':    (6, 6),     # D3h on 3 coplanar vertices (σh fixed)
    'TS':    (2, 2),     # C2v on 3 coplanar vertices
    'OH':    (24, 48),   # Oh
    'SQ':    (8, 8),     # D4h on 4 coplanar vertices (σh fixed)
    'TH':    (12, 24),   # Td
    'TBP':   (6, 12),    # D3h on 5 (axial + equatorial triangle)
    'SP':    (4, 8),     # C4v
    'PBP':   (10, 20),   # D5h
    'SAP':   (8, 16),    # D4d
    'DD':    (4, 8),     # D2d
    'SS':    (2, 2),     # C2v see-saw
    'TPR':   (6, 12),    # D3h trigonal prism
    'COH':   (3, 6),     # C3v capped octahedron
    # Iter-9 additions: CN9 + CN10/11/12.
    'TTP':   (6, 12),    # D3h tricapped trigonal prism
    'BCSAP': (8, 16),    # D4d bicapped square-antiprism
    'PAP':   (10, 20),   # D5d pentagonal antiprism
    'CPAP':  (5, 10),    # C5v capped pentagonal antiprism
    'ICOS':  (60, 120),  # Ih icosahedron
    'CUBO':  (24, 48),   # Oh cuboctahedron
    'HBP':   (12, 24),   # D6h hexagonal prism
}


def _burnside_count(group, n: int, k: int) -> int:
    """Burnside-lemma colour-orbit count: |group|^-1 * Σ k^(#cycles(g))."""
    s = 0
    for g in group:
        seen = [False] * n
        n_cycles = 0
        for i in range(n):
            if not seen[i]:
                n_cycles += 1
                j = i
                while not seen[j]:
                    seen[j] = True
                    j = g[j]
        s += k ** n_cycles
    return s // len(group)


def _enumerate_orbits_via_key(geom: str, n: int, k: int, chiral: bool = False) -> int:
    """Deduplicate exhaustively via :func:`burnside_canonical_key`."""
    symbols = [chr(ord('A') + i) for i in range(k)]
    seen = set()
    for t in itertools.product(symbols, repeat=n):
        seen.add(burnside_canonical_key(geom, t, chiral=chiral))
    return len(seen)


def _max_k_for_n(n: int) -> int:
    """Cap k so the test suite stays under ~30 s wall total.

    Higher-k cells (and full multiset cross-checks) are covered offline by
    ``scripts/polya_audit.py`` against ``results/polya_audit.csv``, so this
    only needs to be tight enough to lock the symmetry group action.
    """
    if n <= 8:
        return 4
    if n <= 9:
        return 3       # TTP n=9: 3**9=19683 fast
    # n=10..12 — Ih/Oh/D5d/C5v/D6h.  k=3 on Ih is ~90 s alone (3**12 = 531441
    # × 120 group elements); k=2 keeps each cell under 100 ms while still
    # discriminating every symmetry element (k=2 finds all bit-pattern
    # collisions among orbits).
    return 2


@pytest.mark.parametrize("geom", sorted(EXPECTED_ORDERS))
def test_geometry_registered(geom):
    """Every polyhedron the enumerator can emit MUST exist in _GEO."""
    assert geom in _GEO, (
        f"Polyhedron {geom!r} missing from _burnside_groups._GEO — "
        f"burnside_canonical_key falls back to identity for this geom, "
        f"so DELFIN_BURNSIDE_FULL=1 is a silent no-op for {geom}."
    )


@pytest.mark.parametrize("geom", sorted(EXPECTED_ORDERS))
def test_group_order(geom):
    """Effective |G+|/|G_full| on vertex labels matches point-group tables."""
    exp_p, exp_f = EXPECTED_ORDERS[geom]
    got_p, got_f = group_size(geom)
    assert got_p == exp_p, f"{geom}: |G+| got {got_p} expected {exp_p}"
    assert got_f == exp_f, f"{geom}: |G_full| got {got_f} expected {exp_f}"


# Build the (geom, k) parametrization with the per-n k cap.
_GEOM_K_CASES = [
    (geom, k)
    for geom in sorted(EXPECTED_ORDERS)
    for k in range(1, _max_k_for_n(len(_GEO[geom])) + 1)
]


@pytest.mark.parametrize("geom,k", _GEOM_K_CASES)
def test_burnside_orbit_count_achiral(geom, k):
    """Orbit count via canonical-key dedup == Burnside-lemma count (full group)."""
    n = len(_GEO[geom])
    proper, full = get_groups(geom)
    lemma = _burnside_count(full, n, k)
    via_key = _enumerate_orbits_via_key(geom, n, k, chiral=False)
    assert lemma == via_key, (
        f"{geom} k={k}: lemma={lemma} via_key={via_key} "
        f"(off-by-{lemma-via_key} — symmetry element missing or overcounted)"
    )


@pytest.mark.parametrize("geom,k", _GEOM_K_CASES)
def test_burnside_orbit_count_chiral(geom, k):
    """Same as above for the proper-rotation subgroup."""
    n = len(_GEO[geom])
    proper, full = get_groups(geom)
    lemma = _burnside_count(proper, n, k)
    via_key = _enumerate_orbits_via_key(geom, n, k, chiral=True)
    assert lemma == via_key, (
        f"{geom} k={k} (chiral): lemma={lemma} via_key={via_key}"
    )


def test_burnside_full_covers_topo_geometry_vectors():
    """Every key in _TOPO_GEOMETRY_VECTORS must have a Burnside vertex set."""
    # Import lazily to avoid pulling RDKit when _TOPO_GEOMETRY_VECTORS is
    # the only thing we need.
    from delfin.smiles_converter import _TOPO_GEOMETRY_VECTORS  # type: ignore
    missing = sorted(set(_TOPO_GEOMETRY_VECTORS) - set(_GEO))
    assert not missing, (
        f"Burnside _GEO is missing {missing} — DELFIN_BURNSIDE_FULL=1 "
        f"would be a silent no-op for these geometries."
    )


def test_burnside_canonical_key_is_orbit_invariant():
    """Applying any group element to a `types` tuple yields the same key."""
    # Spot-check a couple of polyhedra with non-trivial groups.
    cases = [
        ('OH',  ('A', 'B', 'A', 'B', 'C', 'C')),
        ('TPR', ('A', 'B', 'C', 'A', 'B', 'C')),
        ('TTP', ('A', 'B', 'C', 'A', 'B', 'C', 'D', 'D', 'D')),
        ('ICOS', ('A', 'B') * 6),
    ]
    for geom, types in cases:
        proper, full = get_groups(geom)
        n = len(types)
        base = burnside_canonical_key(geom, types, chiral=False)
        for g in full:
            permuted = tuple(types[g[i]] for i in range(n))
            assert burnside_canonical_key(geom, permuted, chiral=False) == base, (
                f"{geom}: key not invariant under group action"
            )
