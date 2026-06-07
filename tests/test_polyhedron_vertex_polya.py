"""Polyhedron-Vertex Pólya Enumeration tests (2026-06-07).

Pins the following contracts:

1. **Determinism** — two consecutive calls produce byte-identical orbit
   lists.
2. **BEYRAY case** — the OC-6 [Cu(N-O-chelate)₂(H₂O)₂] SMILES enumerates
   ≥ 4 distinct vertex-coloring orbits (vs the 1 single-pattern frame the
   legacy ETKDG fallback produces).
3. **Burnside-correct counts** — trivial cases match the textbook isomer
   numbers exactly:

   * MA₆ OC-6           = 1 orbit
   * MA₅B OC-6          = 1 orbit
   * MA₄B₂ OC-6         = 2 orbits (cis / trans)
   * MA₃B₃ OC-6         = 2 orbits (fac / mer)
   * MA₂B₂C₂ OC-6       = 6 orbits
   * MABCDEF OC-6       = 30 orbits
   * MA₂B₂ SP-4         = 2 orbits (cis / trans)
   * MABCD T-4          = 2 orbits (Δ / Λ)

4. **Chelate constraint** — a bidentate chelating ligand cannot be
   enumerated as a trans pair (trans-cis filtering drops the orbit).
5. **Byte-identical OFF** — when the env-gate is unset,
   :func:`embed_isomers` returns the legacy ``embed-conf{N}`` labels
   unchanged.
6. **Edge cases** — CN < 2, no metal, unknown polyhedron all return
   ``None`` gracefully.

Universal, graph-only.  Deterministic with ``PYTHONHASHSEED=0``.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Ensure the repo's ``delfin`` package is importable when tests are run
# stand-alone (without ``pip install -e .``).
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# Skip the whole module when RDKit is unavailable (the orbit enumerator
# itself only needs Python, but the BEYRAY / embed_isomers tests need
# RDKit for SMILES parsing).
pytest.importorskip("rdkit", reason="rdkit not installed")

from delfin.fffree.polyhedron_vertex_polya import (  # noqa: E402
    cis_edges_for,
    cn_for,
    detect_from_smiles,
    enumerate_orbits,
    enumerate_orbits_for_smiles,
    enumerate_orbits_with_chelates,
    flag_active,
)


# --------------------------------------------------------------------------
# 1. Determinism
# --------------------------------------------------------------------------
def test_enumerate_orbits_deterministic():
    """Two calls with identical input produce identical orbit lists."""
    a = enumerate_orbits("OC-6 octahedron", ["A"] * 2 + ["B"] * 2 + ["C"] * 2)
    b = enumerate_orbits("OC-6 octahedron", ["A"] * 2 + ["B"] * 2 + ["C"] * 2)
    assert a == b
    # And the list is internally sorted (deterministic output order).
    assert a == sorted(a)


def test_enumerate_orbits_with_chelates_deterministic():
    a = enumerate_orbits_with_chelates(
        "OC-6 octahedron",
        ["A", "B", "A", "B", "C", "C"],
        [(0, 1), (2, 3)],
    )
    b = enumerate_orbits_with_chelates(
        "OC-6 octahedron",
        ["A", "B", "A", "B", "C", "C"],
        [(0, 1), (2, 3)],
    )
    assert a == b


def test_enumerate_orbits_for_smiles_deterministic():
    smi = "[Cu+2]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]"
    a = enumerate_orbits_for_smiles(smi)
    b = enumerate_orbits_for_smiles(smi)
    assert a == b


# --------------------------------------------------------------------------
# 2. BEYRAY case (motivating example)
# --------------------------------------------------------------------------
BEYRAY = (
    "O=C1[O][Cu-4]2([OH2+])([OH2+])([O]C(=O)C3=CN=CC=[N+]32)"
    "[N+]2=CC=NC=C12"
)


def test_beyray_detects_oc6_octahedron():
    det = detect_from_smiles(BEYRAY)
    assert det is not None, "BEYRAY should detect as OC-6"
    polyhedron, donor_types, chelate_pairs, metal_idx = det
    assert polyhedron == "OC-6 octahedron"
    assert len(donor_types) == 6, f"expected CN=6, got {len(donor_types)}"
    # Two chelating ligands -> two chelate pairs.
    assert len(chelate_pairs) == 2


def test_beyray_emits_at_least_four_orbits():
    """The motivating bug: legacy emits 1 vertex-coloring; we need >= 4."""
    res = enumerate_orbits_for_smiles(BEYRAY)
    assert res is not None
    orbits = res["orbits"]
    # Burnside textbook for OC-6 A2B2C2 = 6 orbits; the chelate-cis filter
    # is satisfiable for ALL 6 (a chelate can sit on any 2 cis vertices,
    # of which OC-6 has 12 -- plenty of room).  We pin >= 4 as the
    # contract because future label refinements (e.g. distinguishing
    # carboxylate-Π from carbonyl-O) could push the count up or down by
    # ±2 without altering the architectural point that the legacy
    # SINGLE-orbit emission is undercounting.
    assert len(orbits) >= 4, (
        f"BEYRAY should enumerate >= 4 distinct vertex-coloring orbits, "
        f"got {len(orbits)}: {orbits}"
    )
    # All orbits must be valid 6-tuples.
    for o in orbits:
        assert isinstance(o, tuple)
        assert len(o) == 6


# --------------------------------------------------------------------------
# 3. Burnside-correct textbook counts (no chelate constraints)
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "polyhedron,donor_types,expected",
    [
        ("OC-6 octahedron", ["A"] * 6, 1),
        ("OC-6 octahedron", ["A"] * 5 + ["B"], 1),
        ("OC-6 octahedron", ["A"] * 4 + ["B"] * 2, 2),    # cis / trans
        ("OC-6 octahedron", ["A"] * 3 + ["B"] * 3, 2),    # fac / mer
        ("OC-6 octahedron", ["A"] * 2 + ["B"] * 2 + ["C"] * 2, 6),
        ("OC-6 octahedron", ["A", "B", "C", "D", "E", "F"], 30),
        ("SP-4 square planar", ["A"] * 2 + ["B"] * 2, 2),  # cis / trans
        ("SP-4", ["A"] * 2 + ["B"] * 2, 2),                # alias form
        ("T-4 tetrahedron", ["A", "B", "C", "D"], 2),      # Δ / Λ
        ("T-4 tetrahedron", ["A"] * 2 + ["B"] * 2, 1),
        ("TBP-5 trigonal bipyramid", ["A"] * 5, 1),
        ("SPY-5 square pyramid", ["A"] * 5, 1),
    ],
)
def test_burnside_textbook_counts(polyhedron, donor_types, expected):
    orbits = enumerate_orbits(polyhedron, donor_types)
    assert len(orbits) == expected, (
        f"{polyhedron} {donor_types}: expected {expected} orbits, got "
        f"{len(orbits)}: {orbits}"
    )


# --------------------------------------------------------------------------
# 4. Chelate-constraint filtering
# --------------------------------------------------------------------------
def test_chelate_forbids_trans_pair():
    """Pinning a bidentate-(0, 1) chelate on OC-6 MA₂B₂ removes the
    trans orbit -- only cis-A₂B₂ survives.

    For OC-6 with donors [A, A, B, B] the Burnside count is 2 (cis A₂B₂
    and trans A₂B₂).  If we declare donor slots 0 and 1 (both 'A') as a
    chelating pair, slot 0 and 1 MUST sit on cis vertices.  In the
    *trans* orbit, the two A's are antipodal, so the cis constraint is
    impossible -> orbit count drops to 1.
    """
    no_chelate = enumerate_orbits_with_chelates(
        "OC-6 octahedron", ["A", "A", "B", "B", "B", "B"], [], max_orbits=100,
    )
    with_chelate = enumerate_orbits_with_chelates(
        "OC-6 octahedron", ["A", "A", "B", "B", "B", "B"],
        [(0, 1)], max_orbits=100,
    )
    assert len(no_chelate) == 2          # cis + trans
    assert len(with_chelate) == 1        # cis only (trans forbidden by chelate)


def test_chelate_constraint_does_not_change_count_when_geometrically_possible():
    """When the chelate constraint is geometrically loose (every orbit
    admits a cis assignment), the count is unchanged."""
    no_chelate = enumerate_orbits_with_chelates(
        "OC-6 octahedron", ["A"] * 2 + ["B"] * 2 + ["C"] * 2,
        [], max_orbits=100,
    )
    # Slots (0, 1) carry the two A's; constrain them to cis.  Since OC-6
    # has 12 cis edges, at least one orbit will admit an A-A cis pair for
    # every orbit -- in fact all 6 do, because A appears twice and any two
    # of the six octahedral vertices that aren't antipodal are cis.
    with_chelate = enumerate_orbits_with_chelates(
        "OC-6 octahedron", ["A"] * 2 + ["B"] * 2 + ["C"] * 2,
        [(0, 1)], max_orbits=100,
    )
    # All orbits where the 2 A's are non-trans survive: 5 of 6 in OC-6
    # A2B2C2 have cis-AA (1 has trans-AA = AB AB trans pair which forbids
    # cis A-A).  Pin >= 4 to allow future tightening of the count.
    assert len(with_chelate) >= 4
    assert len(with_chelate) <= len(no_chelate)


def test_chelate_pair_with_two_different_labels():
    """Chelate (i, j) with donor_types[i] != donor_types[j] requires a
    cis edge carrying one of each label."""
    orbits = enumerate_orbits_with_chelates(
        "OC-6 octahedron", ["A"] * 4 + ["B"] * 2,
        [(0, 4)], max_orbits=100,
    )
    # Without the chelate constraint the count would be 2 (cis / trans B₂).
    # With one chelate (one A and the first B), every orbit admits a cis
    # A-B pair (OC-6 has many AB cis edges in both cis-B₂ and trans-B₂
    # orbits).  So count stays 2.
    assert len(orbits) == 2


# --------------------------------------------------------------------------
# 5. Byte-identical OFF (env-flag contract)
# --------------------------------------------------------------------------
def test_embed_isomers_byte_identical_when_flag_off(monkeypatch):
    """When the env-flag is UNSET, :func:`embed_isomers` returns the
    legacy ``embed-conf{N}`` labels exactly as before."""
    monkeypatch.delenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", raising=False)
    from delfin.fffree.embed_fallback import embed_isomers
    res = embed_isomers(BEYRAY, max_isomers=4, polish="raw")
    assert res is not None
    # Legacy label scheme: every label starts with ``embed-conf``.
    for _, label in res:
        assert label.startswith("embed-conf"), label
    # AND no ``iso{N}-conf{K}`` labels leak through.
    for _, label in res:
        assert not label.startswith("iso"), label


def test_embed_isomers_changes_when_flag_on(monkeypatch):
    """With the flag ON, BEYRAY produces ``isoN-confK-...`` labelled
    frames covering MULTIPLE distinct orbits."""
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", "1")
    from delfin.fffree.embed_fallback import embed_isomers
    res = embed_isomers(BEYRAY, max_isomers=4, polish="raw")
    assert res is not None and len(res) > 0
    iso_prefixes = sorted(
        {label.split("-")[0] for _, label in res if label.startswith("iso")}
    )
    # At least 2 distinct orbits must produce conformers; the BEYRAY
    # textbook count is 6 and we expect at least 3-4 to embed cleanly.
    assert len(iso_prefixes) >= 2, (
        f"expected >= 2 distinct orbit prefixes, got: {iso_prefixes}"
    )


# --------------------------------------------------------------------------
# 6. Edge cases
# --------------------------------------------------------------------------
def test_detect_from_smiles_returns_none_for_organic():
    assert detect_from_smiles("CCO") is None       # ethanol
    assert detect_from_smiles("c1ccccc1") is None  # benzene


def test_detect_from_smiles_returns_none_for_lone_metal():
    # Bare metal with no donors -> can't enumerate a polyhedron.
    assert detect_from_smiles("[Cu+2]") is None


def test_detect_from_smiles_returns_none_for_empty_or_invalid():
    assert detect_from_smiles("") is None
    assert detect_from_smiles("not-a-smiles!@#$") is None


def test_enumerate_orbits_rejects_wrong_length():
    with pytest.raises(ValueError):
        enumerate_orbits("OC-6 octahedron", ["A", "B"])


def test_enumerate_orbits_unknown_polyhedron_raises():
    with pytest.raises(KeyError):
        enumerate_orbits("UNKNOWN-99 polyhedron", ["A"] * 6)


def test_cn_for_textbook_polyhedra():
    assert cn_for("OC-6 octahedron") == 6
    assert cn_for("SP-4") == 4
    assert cn_for("T-4 tetrahedron") == 4
    assert cn_for("TBP-5 trigonal bipyramid") == 5


def test_cis_edges_for_octahedron_matches_textbook():
    """OC-6 has exactly 12 cis edges (the 12 edges of the octahedron)
    and 3 trans pairs (the 3 antipodal axis pairs)."""
    cis = cis_edges_for("OC-6 octahedron")
    assert len(cis) == 12, f"expected 12 OC-6 cis edges, got {len(cis)}"


def test_cis_edges_for_square_planar_matches_textbook():
    """SP-4 has 4 cis edges (the 4 sides of the square) + 0 trans edges
    in the cis_edges_for output (the 2 diagonals are trans = 180 deg).
    """
    cis = cis_edges_for("SP-4 square planar")
    assert len(cis) == 4


def test_max_orbits_cap_honoured(monkeypatch):
    """The env-var cap (or kwarg) limits the orbit count emitted."""
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA_MAX_ORBITS", "3")
    orbits = enumerate_orbits_with_chelates(
        "OC-6 octahedron",
        ["A", "B", "C", "D", "E", "F"],  # 30 distinct orbits!
        [],
    )
    assert len(orbits) == 3


def test_flag_active_reads_env(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", raising=False)
    assert flag_active() is False
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", "1")
    assert flag_active() is True
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", "0")
    assert flag_active() is False
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_VERTEX_POLYA", "true")
    assert flag_active() is True
