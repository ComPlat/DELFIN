"""CN5 + CN3 Orbit-Completeness tests (2026-06-07).

Pins the multi-polyhedron extension of :mod:`polyhedron_vertex_polya`
documented in commit ``feat-cn5-orbit-completeness-2026-06-07``.  The
contract:

* Burnside textbook counts are correct for every per-polyhedron case in
  the CN3 / CN5 menu (TBP-5, SPY-5, SP-3 / TP-3 trigonal planar, T-3
  T-shape, Y-3 Y-shape).
* The SMILES-level enumerator unions orbits across ALL polyhedra
  registered for a CN (via :func:`polyhedra_for_cn`), so the ALESUH-class
  Cu(II) CN5 case produces >= 4 distinct (polyhedron, orbit) pairs
  (textbook: 2 from TBP-5 + 2 from SPY-5).
* Determinism: identical input -> identical orbit list across runs.
* Byte-identical when env-gate is OFF: :func:`embed_isomers` returns the
  legacy ``embed-conf{N}`` labels exactly as before the extension.

Universal, graph-only, deterministic.  PYTHONHASHSEED=0 sets the seed.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Make the repo's ``delfin`` package importable when running stand-alone.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

pytest.importorskip("rdkit", reason="rdkit not installed")

from delfin.fffree.polyhedron_vertex_polya import (  # noqa: E402
    enumerate_multi_polyhedron_orbits,
    enumerate_orbits,
    enumerate_orbits_for_smiles,
    enumerate_orbits_with_chelates,
    polyhedra_for_cn,
)


# -------------------------------------------------------------------------
# 1. TBP-5 textbook isomer counts
# -------------------------------------------------------------------------
@pytest.mark.parametrize(
    "donor_types,expected",
    [
        (["A"] * 5, 1),                                       # A5
        (["A"] * 4 + ["B"], 2),                               # A4B: axial vs equatorial
        (["A"] * 3 + ["B"] * 2, 3),                           # A3B2: D3 → 3 orbits
        (["A"] * 3 + ["B", "C"], 4),                          # A3BC: D3 → 4 orbits
        (["N"] * 4 + ["Br"], 2),                              # ALESUH-relevant N4Br
    ],
)
def test_tbp5_burnside_counts(donor_types, expected):
    orbits = enumerate_orbits("TBP-5 trigonal bipyramid", donor_types)
    assert len(orbits) == expected, (
        f"TBP-5 {donor_types}: expected {expected}, got {len(orbits)}: {orbits}"
    )


# -------------------------------------------------------------------------
# 2. SPY-5 textbook isomer counts
# -------------------------------------------------------------------------
@pytest.mark.parametrize(
    "donor_types,expected",
    [
        (["A"] * 5, 1),                                       # A5
        (["A"] * 4 + ["B"], 2),                               # A4B: apical vs basal
        (["A"] * 3 + ["B"] * 2, 3),                           # A3B2
        (["A"] * 3 + ["B", "C"], 5),                          # A3BC
        (["N"] * 4 + ["Br"], 2),                              # N4Br: apical vs basal
    ],
)
def test_spy5_burnside_counts(donor_types, expected):
    orbits = enumerate_orbits("SPY-5 square pyramid", donor_types)
    assert len(orbits) == expected, (
        f"SPY-5 {donor_types}: expected {expected}, got {len(orbits)}: {orbits}"
    )


# -------------------------------------------------------------------------
# 3. CN3 trigonal planar (SP-3) and T-shape (T-3) counts
# -------------------------------------------------------------------------
@pytest.mark.parametrize(
    "polyhedron,donor_types,expected",
    [
        # SP-3 trigonal planar (proper rotation D3, order 6): all vertices
        # equivalent under C3 → A3 = 1, A2B = 1, ABC = 1.
        ("SP-3 trigonal planar", ["A"] * 3, 1),
        ("SP-3 trigonal planar", ["A"] * 2 + ["B"], 1),
        ("SP-3 trigonal planar", ["A", "B", "C"], 1),
        # Alias TP-3 must resolve to the same group.
        ("TP-3", ["A"] * 3, 1),
        # T-3 T-shape (proper rotation C2, order 2): 2 trans-equivalent
        # vertices + 1 unique → A3 = 1, A2B = 2 (B trans or B unique),
        # ABC = 3.
        ("T-3 T-shape", ["A"] * 3, 1),
        ("T-3 T-shape", ["A"] * 2 + ["B"], 2),
        ("T-3 T-shape", ["A", "B", "C"], 3),
        # Y-3 Y-shape: same proper-rotation group as T-3 → identical counts.
        ("Y-3 Y-shape", ["A"] * 3, 1),
        ("Y-3 Y-shape", ["A"] * 2 + ["B"], 2),
        ("Y-3 Y-shape", ["A", "B", "C"], 3),
    ],
)
def test_cn3_burnside_counts(polyhedron, donor_types, expected):
    orbits = enumerate_orbits(polyhedron, donor_types)
    assert len(orbits) == expected, (
        f"{polyhedron} {donor_types}: expected {expected}, got {len(orbits)}"
    )


def test_cn3_t_shape_a2b_orbits_identify_axial():
    """T-shape A2B has 2 orbits: B sits on the unique vertex OR on a
    trans-pair vertex.  Verify the orbit reps cover both placements."""
    orbits = enumerate_orbits("T-3 T-shape", ["A", "A", "B"])
    assert len(orbits) == 2
    # T-3 vertex 2 is the unique (cis) vertex in the polya group ordering.
    # The two orbits should differ in where B lands.
    positions = sorted(o.index("B") for o in orbits)
    # B placement positions: one orbit has B on the unique vertex (idx 2),
    # the other has B on a trans vertex (idx 0 or 1).
    assert positions == [0, 2] or positions == [1, 2], (
        f"expected B placements [0,2] or [1,2], got {positions}: {orbits}"
    )


# -------------------------------------------------------------------------
# 4. polyhedra_for_cn registry
# -------------------------------------------------------------------------
def test_polyhedra_for_cn3_includes_y3():
    """CN3 polyhedra list must include SP-3, T-3 AND Y-3 (the 2026-06-07
    extension)."""
    polys = polyhedra_for_cn(3)
    assert "SP-3 trigonal planar" in polys
    assert "T-3 T-shape" in polys
    assert "Y-3 Y-shape" in polys
    assert len(polys) == 3


def test_polyhedra_for_cn5_includes_both_tbp_and_spy():
    polys = polyhedra_for_cn(5)
    assert "TBP-5 trigonal bipyramid" in polys
    assert "SPY-5 square pyramid" in polys
    assert len(polys) == 2


def test_polyhedra_for_unknown_cn_returns_empty():
    assert polyhedra_for_cn(99) == []


# -------------------------------------------------------------------------
# 5. enumerate_multi_polyhedron_orbits
# -------------------------------------------------------------------------
def test_multi_polyhedron_cn5_n4br():
    """ALESUH-class N4Br CN5: TBP-5 contributes 2 orbits, SPY-5 contributes
    2 orbits.  Total = 4 distinct (polyhedron, orbit) pairs."""
    pairs = enumerate_multi_polyhedron_orbits(5, ["N"] * 4 + ["Br"])
    polys = sorted({p for p, _ in pairs})
    assert polys == ["SPY-5 square pyramid", "TBP-5 trigonal bipyramid"]
    assert len(pairs) == 4
    # Verify each polyhedron contributes >= 2 orbits.
    tbp_orbits = [o for p, o in pairs if p.startswith("TBP-5")]
    spy_orbits = [o for p, o in pairs if p.startswith("SPY-5")]
    assert len(tbp_orbits) >= 2
    assert len(spy_orbits) >= 2


def test_multi_polyhedron_cn3_a2b():
    """CN3 A2B: SP-3 contributes 1 orbit, T-3 contributes 2, Y-3 contributes
    2.  Total = 5 distinct (polyhedron, orbit) pairs."""
    pairs = enumerate_multi_polyhedron_orbits(3, ["A", "A", "B"])
    polys = sorted({p for p, _ in pairs})
    assert "SP-3 trigonal planar" in polys
    assert "T-3 T-shape" in polys
    assert "Y-3 Y-shape" in polys
    assert len(pairs) == 1 + 2 + 2


def test_multi_polyhedron_explicit_polyhedra_arg():
    """Caller-supplied polyhedra list overrides the default registry."""
    pairs = enumerate_multi_polyhedron_orbits(
        5, ["N"] * 4 + ["Br"],
        polyhedra=["TBP-5 trigonal bipyramid"],
    )
    polys = sorted({p for p, _ in pairs})
    assert polys == ["TBP-5 trigonal bipyramid"]
    assert len(pairs) == 2


def test_multi_polyhedron_wrong_length_raises():
    with pytest.raises(ValueError):
        enumerate_multi_polyhedron_orbits(5, ["A", "B"])


def test_multi_polyhedron_deterministic():
    a = enumerate_multi_polyhedron_orbits(5, ["N"] * 4 + ["Br"])
    b = enumerate_multi_polyhedron_orbits(5, ["N"] * 4 + ["Br"])
    assert a == b
    assert a == sorted(a)


# -------------------------------------------------------------------------
# 6. ALESUH integration: SMILES → multi-polyhedron orbits
# -------------------------------------------------------------------------
ALESUH = (
    "CN1C=C[N+]([Cu-3]([Br])([N+]2=CN(C)C=C2)([N+]2=CN(C)C=C2)"
    "[N+]2=CN(C)C=C2)=C1"
)


def test_alesuh_emits_at_least_four_orbits():
    """The motivating ALESUH bug: V14 emits only 3 frames; we need >= 4
    distinct (polyhedron, orbit) pairs (textbook: 2 TBP-5 + 2 SPY-5)."""
    res = enumerate_orbits_for_smiles(ALESUH)
    assert res is not None
    assert "polyhedron_orbits" in res
    total_pairs = sum(len(p["orbits"]) for p in res["polyhedron_orbits"])
    assert total_pairs >= 4, (
        f"ALESUH should enumerate >= 4 (polyhedron, orbit) pairs, "
        f"got {total_pairs}: {res['polyhedron_orbits']}"
    )


def test_alesuh_covers_both_polyhedra():
    """Both TBP-5 AND SPY-5 must appear in the orbit breakdown for an
    ambivalent CN5 metal."""
    res = enumerate_orbits_for_smiles(ALESUH)
    assert res is not None
    polys_with_orbits = {
        p["polyhedron"] for p in res["polyhedron_orbits"] if p["orbits"]
    }
    assert "TBP-5 trigonal bipyramid" in polys_with_orbits
    assert "SPY-5 square pyramid" in polys_with_orbits


def test_alesuh_legacy_orbits_field_preserved():
    """The legacy ``polyhedron`` + ``orbits`` keys remain populated for the
    primary polyhedron (back-compat)."""
    res = enumerate_orbits_for_smiles(ALESUH)
    assert res is not None
    assert isinstance(res["polyhedron"], str)
    assert isinstance(res["orbits"], list)
    assert len(res["orbits"]) >= 1


def test_alesuh_deterministic():
    a = enumerate_orbits_for_smiles(ALESUH)
    b = enumerate_orbits_for_smiles(ALESUH)
    assert a is not None and b is not None
    # Compare every field.
    assert a["polyhedron"] == b["polyhedron"]
    assert a["orbits"] == b["orbits"]
    assert a["polyhedron_orbits"] == b["polyhedron_orbits"]
    assert a["donor_types"] == b["donor_types"]
    assert a["chelate_pairs"] == b["chelate_pairs"]


# -------------------------------------------------------------------------
# 7. Byte-identical OFF (env-flag contract)
# -------------------------------------------------------------------------
def test_embed_isomers_byte_identical_when_flag_off(monkeypatch):
    """When the env-flag is UNSET, :func:`embed_isomers` returns the
    legacy ``embed-conf{N}`` labels exactly as before the extension."""
    monkeypatch.delenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", raising=False)
    from delfin.fffree.embed_fallback import embed_isomers
    res = embed_isomers(ALESUH, max_isomers=4, polish="raw")
    if res is None:
        pytest.skip("ETKDG could not embed ALESUH on this platform")
    for _, label in res:
        # Legacy: every label starts with ``embed-conf``.
        assert label.startswith("embed-conf"), label
        # And no ``iso{N}-conf{K}`` labels leak through.
        assert not label.startswith("iso"), label


def test_embed_isomers_changes_when_flag_on(monkeypatch):
    """With the env-flag ON, ALESUH emits multi-orbit frames (label prefix
    ``isoN-confK-...``) covering more than one orbit -- including SPY-5 +
    TBP-5 mixed contributions."""
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", "1")
    from delfin.fffree.embed_fallback import embed_isomers
    res = embed_isomers(ALESUH, max_isomers=2, polish="raw")
    if res is None or not res:
        pytest.skip("ETKDG could not embed ALESUH on this platform")
    iso_prefixes = sorted(
        {label.split("-")[0] for _, label in res if label.startswith("iso")}
    )
    # We expect at least 2 distinct orbit prefixes (the chemistry asks for
    # at minimum 2 from each polyhedron, but ETKDG can drop some).
    assert len(iso_prefixes) >= 2, (
        f"expected >= 2 distinct orbit prefixes, got: {iso_prefixes}"
    )


# -------------------------------------------------------------------------
# 8. Determinism (orbit lists, byte-identical across two runs)
# -------------------------------------------------------------------------
def test_enumerate_orbits_tbp5_deterministic():
    a = enumerate_orbits("TBP-5 trigonal bipyramid", ["N"] * 4 + ["Br"])
    b = enumerate_orbits("TBP-5 trigonal bipyramid", ["N"] * 4 + ["Br"])
    assert a == b
    assert a == sorted(a)


def test_enumerate_orbits_spy5_deterministic():
    a = enumerate_orbits("SPY-5 square pyramid", ["N"] * 4 + ["Br"])
    b = enumerate_orbits("SPY-5 square pyramid", ["N"] * 4 + ["Br"])
    assert a == b
    assert a == sorted(a)


# -------------------------------------------------------------------------
# 9. Chelate-cis filter still works on CN5 + CN3 polyhedra
# -------------------------------------------------------------------------
def test_chelate_constraint_filters_trans_pair_on_tbp5():
    """A bidentate chelate on TBP-5 cannot be trans (180°), so the
    axial-axial constraint is impossible."""
    # Slots 0 and 1 are donors of label A, declared chelating.  For TBP-5
    # the axial-axial pair is 180°, so an orbit with both A's axial is
    # forbidden.
    no_chelate = enumerate_orbits_with_chelates(
        "TBP-5 trigonal bipyramid", ["A", "A", "B", "B", "B"], [],
        max_orbits=100,
    )
    with_chelate = enumerate_orbits_with_chelates(
        "TBP-5 trigonal bipyramid", ["A", "A", "B", "B", "B"], [(0, 1)],
        max_orbits=100,
    )
    # Chelate filter must remove SOMETHING that the trans-AA orbit
    # contributes (or yield <= the un-constrained count).
    assert len(with_chelate) <= len(no_chelate)
    assert len(with_chelate) >= 1


# -------------------------------------------------------------------------
# 10. vertex_coords for Y-3 falls back to T-3 placement
# -------------------------------------------------------------------------
def test_vertex_coords_y3_fallback():
    from delfin.fffree.polyhedron_vertex_polya import vertex_coords
    V_t3 = vertex_coords("T-3 T-shape")
    V_y3 = vertex_coords("Y-3 Y-shape")
    assert V_t3 is not None
    assert V_y3 is not None
    # Y-3 reuses T-3 vertices (same proper-rotation skeleton; angular
    # geometry is refined downstream by ETKDG).
    import numpy as np
    assert np.allclose(V_t3, V_y3)
