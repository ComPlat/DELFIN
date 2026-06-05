"""Tests for non-f-block CN10 polyhedra (Mission A2, 2026-06-05).

Adds CN10 polyhedra BICAP-10 / CSAP-10 / SAP-10 so non-f-block CN10
complexes (e.g. Y, early-d, large-anion) can be built by the FF-free
constructor instead of falling through to the legacy UFF path.

Verifies:
  * Each CN10 polyhedron returns 10 unit-norm vertex vectors.
  * Vertex sets are well-separated (no near-coincident vertices).
  * BICAP-10 unit-vectors match the f-block module's BICAP-10 (single source
    of truth — the non-f-block path REUSES the existing f-block builder).
  * Env-OFF byte-identical: with ``DELFIN_FFFREE_CN10_POLYHEDRA`` and
    ``DELFIN_FFFREE_PURE_TRACK3`` unset, the polyhedra / decompose dispatch
    is identical to HEAD (CN10 returns empty / None for non-f-block).
  * Env-ON: ``geometries_for_cn(10, 'Y')`` returns the 3 polyhedra; the
    decompose default for Y CN10 becomes BICAP-10.
  * Determinism: two consecutive calls produce byte-identical np arrays.
  * Synthetic [Y(NH3)10] assemble produces 10 N donors at identical M-N
    distance with zero heavy-heavy <1Å duplicates (the "atom-collapse"
    pathology the legacy CN10 path was prone to).
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest


@pytest.fixture
def env_off(monkeypatch):
    """All CN10 env flags unset — legacy byte-identical path."""
    monkeypatch.delenv("DELFIN_FFFREE_CN10_POLYHEDRA", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_FBLOCK_CN8_12", raising=False)
    return None


@pytest.fixture
def env_on(monkeypatch):
    """DELFIN_FFFREE_CN10_POLYHEDRA=1 enables the new CN10 dispatch."""
    monkeypatch.setenv("DELFIN_FFFREE_CN10_POLYHEDRA", "1")
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_FBLOCK_CN8_12", raising=False)
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    return None


# ---------------------------------------------------------------------------
# 1. Vertex-set sanity (unit-norm, well-separated, 10 vertices each).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [
    "BICAP-10 bicapped square antiprism",
    "CSAP-10 capped square antiprism",
    "SAP-10 pentagonal antiprism",
])
def test_cn10_polyhedra_unit_norm_10_vertices(name):
    from delfin.fffree.polyhedra import CN10_VERTICES
    fn = CN10_VERTICES[name]
    v = fn()
    assert v.shape == (10, 3), f"{name}: expected (10,3), got {v.shape}"
    norms = np.linalg.norm(v, axis=1)
    assert float(np.max(np.abs(norms - 1.0))) < 1e-12, \
        f"{name}: vertices not unit-norm"


@pytest.mark.parametrize("name", [
    "BICAP-10 bicapped square antiprism",
    "CSAP-10 capped square antiprism",
    "SAP-10 pentagonal antiprism",
])
def test_cn10_no_coincident_vertices(name):
    from delfin.fffree.polyhedra import CN10_VERTICES
    v = CN10_VERTICES[name]()
    min_angle_deg = 180.0
    for i in range(10):
        for j in range(i + 1, 10):
            c = float(np.clip(np.dot(v[i], v[j]), -1.0, 1.0))
            ang = math.degrees(math.acos(c))
            if ang < min_angle_deg:
                min_angle_deg = ang
    # All three CN10 polyhedra should keep vertices well apart (>30°).
    assert min_angle_deg > 30.0, \
        f"{name}: min vertex-vertex angle {min_angle_deg:.2f}° too small"


def test_bicap_10_matches_fblock_module():
    """BICAP-10 in polyhedra.py must call the f-block builder (single source).

    Guards against silent divergence between the f-block CN10 default and the
    non-f-block CN10 dispatch.
    """
    from delfin.fffree.polyhedra import _bicap_10_unit
    from delfin.fffree.f_block_polyhedra import bicap_10_vertices
    a = _bicap_10_unit()
    b = bicap_10_vertices()
    assert np.array_equal(a, b), "BICAP-10 in polyhedra.py != f_block_polyhedra"


# ---------------------------------------------------------------------------
# 2. Env-OFF byte-identical (no CN10 surfaces for non-f-block when flags unset).
# ---------------------------------------------------------------------------


def test_env_off_geometries_for_cn10_empty(env_off):
    from delfin.fffree.polyhedra import geometries_for_cn
    # With ALL CN10/PT3/FBLOCK flags unset, non-f-block CN10 returns [] (HEAD).
    assert geometries_for_cn(10) == []
    assert geometries_for_cn(10, "Y") == []
    assert geometries_for_cn(10, "La") == []


def test_env_off_decompose_default_geom_cn10_none(env_off):
    from delfin.fffree.decompose import _default_geometry
    # HEAD behaviour: non-f-block CN10 -> None (caller falls to legacy).
    assert _default_geometry("Y", 10) is None
    assert _default_geometry("Sc", 10) is None


def test_env_off_legacy_cn_unchanged(env_off):
    """CN<10 geometries are byte-identical to HEAD when flags unset."""
    from delfin.fffree.polyhedra import geometries_for_cn, GEOM_BY_CN
    for cn, expected in GEOM_BY_CN.items():
        assert geometries_for_cn(cn) == expected
        # metal-aware path: same answer when no f-block / no CN10 flags set
        assert geometries_for_cn(cn, "Y") == expected


# ---------------------------------------------------------------------------
# 3. Env-ON wiring: 3 polyhedra surface for non-f-block CN10.
# ---------------------------------------------------------------------------


def test_env_on_geometries_for_cn10_yields_three(env_on):
    from delfin.fffree.polyhedra import geometries_for_cn
    got = geometries_for_cn(10, "Y")
    assert "BICAP-10 bicapped square antiprism" in got
    assert "CSAP-10 capped square antiprism" in got
    assert "SAP-10 pentagonal antiprism" in got
    assert len(got) == 3


def test_env_on_decompose_default_geom_y_cn10_bicap(env_on):
    from delfin.fffree.decompose import _default_geometry
    # Default for non-f-block CN10 = BICAP-10 (chemically most common shape).
    assert _default_geometry("Y", 10) == "BICAP-10 bicapped square antiprism"
    assert _default_geometry("Sc", 10) == "BICAP-10 bicapped square antiprism"


def test_env_on_ref_vectors_dispatches_cn10(env_on):
    """``polyhedra.ref_vectors`` must return the new vertex sets by name."""
    from delfin.fffree.polyhedra import ref_vectors
    for nm in ("BICAP-10 bicapped square antiprism",
               "CSAP-10 capped square antiprism",
               "SAP-10 pentagonal antiprism"):
        v = ref_vectors(nm)
        assert v.shape == (10, 3)


# ---------------------------------------------------------------------------
# 4. Determinism: two consecutive calls produce byte-identical arrays.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [
    "BICAP-10 bicapped square antiprism",
    "CSAP-10 capped square antiprism",
    "SAP-10 pentagonal antiprism",
])
def test_cn10_determinism_two_runs(name):
    from delfin.fffree.polyhedra import CN10_VERTICES
    a = CN10_VERTICES[name]()
    b = CN10_VERTICES[name]()
    assert np.array_equal(a, b), f"{name}: two calls disagreed"


# ---------------------------------------------------------------------------
# 5. Integration: assemble [Y(NH3)10] for all 3 polyhedra has zero
# atom-collapse and consistent M-N distance.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("geometry", [
    "BICAP-10 bicapped square antiprism",
    "CSAP-10 capped square antiprism",
    "SAP-10 pentagonal antiprism",
])
def test_cn10_assemble_y_nh3_10_no_collapse(env_on, geometry):
    """Synthetic [Y(NH3)10] build: all 10 N donors should land at the
    designed M-N distance with NO heavy-heavy <1Å duplicates."""
    import itertools
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from delfin.fffree import assemble_complex as A

    ligands = []
    for _ in range(10):
        m = Chem.AddHs(Chem.MolFromSmiles("N"))
        AllChem.EmbedMolecule(m, randomSeed=42)
        ligands.append({
            "mol": m, "donor_local_idx": 0, "donor_local_idxs": [0],
            "denticity": 1, "donor_elem": "N", "donor_elems": ["N"],
            "is_hapto": False,
        })
    config = {i: (i, 0) for i in range(10)}
    res = A.assemble_from_config("Y", geometry, config, ligands, refine=False)
    assert res is not None, f"{geometry}: assemble returned None"
    syms, P, _ = res
    P = np.asarray(P, float)
    m_pos = P[0]
    n_dists = [float(np.linalg.norm(P[i] - m_pos))
               for i, s in enumerate(syms) if s == "N"]
    assert len(n_dists) == 10, f"{geometry}: expected 10 N, got {len(n_dists)}"
    # All M-N distances within 0.05Å of each other (ideal CN10 sphere).
    assert max(n_dists) - min(n_dists) < 0.05, \
        f"{geometry}: M-N spread {max(n_dists)-min(n_dists):.3f}A too large"
    # Heavy-heavy <1Å collapse count.
    heavy = [(i, s, P[i]) for i, s in enumerate(syms) if s != "H" and i != 0]
    dups = 0
    for (i1, e1, p1), (i2, e2, p2) in itertools.combinations(heavy, 2):
        if float(np.linalg.norm(p1 - p2)) < 1.0:
            dups += 1
    assert dups == 0, f"{geometry}: {dups} heavy-heavy <1A duplicates"


# ---------------------------------------------------------------------------
# 6. f-block metal (La) at CN10 still works under env_on (f-block dispatch
# co-existence) — La must NOT be re-routed to the non-f-block injection path.
# ---------------------------------------------------------------------------


def test_env_on_la_cn10_still_handled(env_on):
    """``geometries_for_cn(10, 'La')`` under CN10_POLYHEDRA=1 returns [] (the
    f-block branch is gated by a DIFFERENT flag, FBLOCK_CN8_12).

    This codifies the separation-of-concerns: CN10_POLYHEDRA does NOT enable
    the f-block CN8-12 dispatch — that flag is independent.
    """
    from delfin.fffree.polyhedra import geometries_for_cn
    # CN10_POLYHEDRA=1 alone, La is f-block -> the CN10 injection is skipped
    # (we don't want to inject CSAP-10/SAP-10 for f-block, which has its own
    # canonical CN10 in BICAP-10).
    got = geometries_for_cn(10, "La")
    assert got == []  # f-block dispatch not enabled by this flag
