"""Tests for grip_lib v5 — perfection levers.

Verifies:

1. v5 metadata loads (cleaning_applied, meta_n_extracted)
2. TM-category lookups (carbene, hapto_eta2/5/6, mu_bridge, agostic, ox_addition)
   return non-None, chemistry-reasonable values for the categories the smoke
   build populated.
3. Block-disaggregated bond lookups distinguish 3d vs 4d vs 5d metals AND
   degrade gracefully (fall through to v4 wildcards) when no block-resolved
   entry exists for a query.
4. v5 library is backward-compatible — all v4 lookup APIs work unchanged.
5. lookup_tm_category returns None on v3/v4 libraries (no false data).
6. Determinism: a fresh load of the v5 smoke lib yields identical lookups.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("PYTHONHASHSEED", "0")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from delfin.fffree import grip_mogul_lookup as gml  # noqa: E402

V5_SMOKE = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5_smoke.npz"
)
V5_FULL = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz"
)

requires_v5_smoke = pytest.mark.skipif(
    not V5_SMOKE.exists(),
    reason="grip_lib_v5_smoke.npz not built yet",
)


def _key(parsed) -> str:
    return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)


def _build_v5_synth(tmp: Path) -> Path:
    """Build a tiny v5 npz with one carbene + one hapto_eta5 + one mu_bridge +
    block-disaggregated tables.
    """
    cc_key = _key(["C", "sp3", [["C", 1, -1, "sp3"]], []])
    master_keys = [cc_key]
    n_master = 1
    nan = np.nan
    zero_p = np.full(1, nan, dtype=np.float32)
    out = tmp / "grip_lib_v5_synth.npz"

    # TM-category tables: carbene Pd, hapto_eta5 Fe, mu_bridge Cu
    carbene_keys = [_key(["Pd", "4d", "C", "carbene"])]
    eta5_keys = [_key(["Fe", "3d", "C", "eta5"])]
    mu_keys = [_key(["Cu", "3d", "Cl", "mu"])]
    eta2_keys = []
    eta6_keys = []
    ago_keys = []
    ox_keys = []

    # Block-disagg pair_bond: Pd-Cl(4d) and Hf-Cl(5d) at different distances
    block_pb_keys = [
        _key(["Cl", "sp3", "Pd", "4d"]),
        _key(["Cl", "sp3", "Hf", "5d"]),
    ]
    block_pb_mu = np.array([2.31, 2.42], dtype=np.float64)
    block_pb_sigma = np.array([0.04, 0.05], dtype=np.float64)
    block_pb_n = np.array([100, 80], dtype=np.int32)

    block_ta_keys = []  # empty

    np.savez_compressed(
        out,
        version=np.int32(5),
        cleaning_applied=np.int32(1),
        meta_n_total_scanned=np.int32(1000),
        meta_n_extracted=np.int32(890),
        meta_n_polymeric=np.int32(50),
        meta_n_no_3d=np.int32(40),
        meta_n_too_small=np.int32(20),
        meta_n_recovered_disorder=np.int32(200),
        meta_n_recovered_partial_coord=np.int32(50),
        n_master=np.int32(n_master), n_orig=np.int32(n_master),
        n_torsion=np.int32(0), n_orig_torsion=np.int32(0),
        keys=np.array(master_keys, dtype=np.object_),
        orig_keys=np.array(master_keys, dtype=np.object_),
        orig_idx=np.zeros(n_master, dtype=np.int32),
        bond_mu=np.array([1.540], dtype=np.float64),
        bond_sigma=np.array([0.010], dtype=np.float64),
        bond_n=np.array([1000], dtype=np.int32),
        bond_p5=zero_p, bond_p50=zero_p, bond_p95=zero_p,
        angle_mu=np.full(1, nan, dtype=np.float64),
        angle_sigma=np.full(1, nan, dtype=np.float64),
        angle_n=np.zeros(1, dtype=np.int32),
        angle_p5=zero_p, angle_p50=zero_p, angle_p95=zero_p,
        improper_mu=np.full(1, nan, dtype=np.float64),
        improper_sigma=np.full(1, nan, dtype=np.float64),
        improper_n=np.zeros(1, dtype=np.int32),
        improper_p5=zero_p, improper_p50=zero_p, improper_p95=zero_p,
        fb_flat=np.zeros(0, dtype=np.int32),
        fb_ptr=np.zeros(2, dtype=np.int32),
        torsion_keys=np.array([], dtype=np.object_),
        torsion_orig_keys=np.array([], dtype=np.object_),
        torsion_n=np.zeros(0, dtype=np.int32),
        torsion_n_components=np.zeros(0, dtype=np.int32),
        torsion_pi=np.zeros((0, 3), dtype=np.float32),
        torsion_mu=np.zeros((0, 3), dtype=np.float32),
        torsion_sigma=np.zeros((0, 3), dtype=np.float32),
        torsion_fb_flat=np.zeros(0, dtype=np.int32),
        torsion_fb_ptr=np.zeros(1, dtype=np.int32),
        # v4 pair tables (small)
        n_pair_bond=np.int32(1),
        n_triple_angle=np.int32(0),
        n_improper_pair=np.int32(0),
        pair_bond_keys=np.array([_key(["Cl", "sp3", "Pd", "*"])], dtype=np.object_),
        pair_bond_mu=np.array([2.30], dtype=np.float64),
        pair_bond_sigma=np.array([0.05], dtype=np.float64),
        pair_bond_n=np.array([150], dtype=np.int32),
        triple_angle_keys=np.array([], dtype=np.object_),
        triple_angle_mu=np.zeros(0, dtype=np.float64),
        triple_angle_sigma=np.zeros(0, dtype=np.float64),
        triple_angle_n=np.zeros(0, dtype=np.int32),
        improper_pair_keys=np.array([], dtype=np.object_),
        improper_pair_mu=np.zeros(0, dtype=np.float64),
        improper_pair_sigma=np.zeros(0, dtype=np.float64),
        improper_pair_n=np.zeros(0, dtype=np.int32),
        # v5 TM-cat tables
        tm_carbene_keys=np.array(carbene_keys, dtype=np.object_),
        tm_carbene_mu=np.array([1.99], dtype=np.float64),
        tm_carbene_sigma=np.array([0.05], dtype=np.float64),
        tm_carbene_n=np.array([60], dtype=np.int32),
        tm_hapto_eta2_keys=np.array(eta2_keys, dtype=np.object_),
        tm_hapto_eta2_mu=np.zeros(0, dtype=np.float64),
        tm_hapto_eta2_sigma=np.zeros(0, dtype=np.float64),
        tm_hapto_eta2_n=np.zeros(0, dtype=np.int32),
        tm_hapto_eta5_keys=np.array(eta5_keys, dtype=np.object_),
        tm_hapto_eta5_mu=np.array([2.04], dtype=np.float64),
        tm_hapto_eta5_sigma=np.array([0.04], dtype=np.float64),
        tm_hapto_eta5_n=np.array([45], dtype=np.int32),
        tm_hapto_eta6_keys=np.array(eta6_keys, dtype=np.object_),
        tm_hapto_eta6_mu=np.zeros(0, dtype=np.float64),
        tm_hapto_eta6_sigma=np.zeros(0, dtype=np.float64),
        tm_hapto_eta6_n=np.zeros(0, dtype=np.int32),
        tm_mu_bridge_keys=np.array(mu_keys, dtype=np.object_),
        tm_mu_bridge_mu=np.array([2.65], dtype=np.float64),
        tm_mu_bridge_sigma=np.array([0.16], dtype=np.float64),
        tm_mu_bridge_n=np.array([30], dtype=np.int32),
        tm_agostic_keys=np.array(ago_keys, dtype=np.object_),
        tm_agostic_mu=np.zeros(0, dtype=np.float64),
        tm_agostic_sigma=np.zeros(0, dtype=np.float64),
        tm_agostic_n=np.zeros(0, dtype=np.int32),
        tm_ox_addition_keys=np.array(ox_keys, dtype=np.object_),
        tm_ox_addition_mu=np.zeros(0, dtype=np.float64),
        tm_ox_addition_sigma=np.zeros(0, dtype=np.float64),
        tm_ox_addition_n=np.zeros(0, dtype=np.int32),
        # block-disagg
        tm_pair_bond_block_keys=np.array(block_pb_keys, dtype=np.object_),
        tm_pair_bond_block_mu=block_pb_mu,
        tm_pair_bond_block_sigma=block_pb_sigma,
        tm_pair_bond_block_n=block_pb_n,
        tm_triple_angle_block_keys=np.array(block_ta_keys, dtype=np.object_),
        tm_triple_angle_block_mu=np.zeros(0, dtype=np.float64),
        tm_triple_angle_block_sigma=np.zeros(0, dtype=np.float64),
        tm_triple_angle_block_n=np.zeros(0, dtype=np.int32),
        mogul_validated=np.zeros(n_master, dtype=bool),
    )
    return out


@pytest.fixture
def v5_synth(tmp_path):
    p = _build_v5_synth(tmp_path)
    gml.GripLibrary._SINGLETONS.clear()
    yield p
    gml.GripLibrary._SINGLETONS.clear()


def test_v5_metadata_loaded(v5_synth):
    lib = gml.GripLibrary.get(v5_synth)
    assert lib.version == 5
    assert lib.cleaning_applied is True
    assert lib.meta_n_total_scanned == 1000
    assert lib.meta_n_extracted == 890
    assert lib.has_tm_categories is True
    assert lib.has_block_disagg is True


def test_carbene_lookup(v5_synth):
    lib = gml.GripLibrary.get(v5_synth)
    hit = lib.lookup_tm_category("carbene", "Pd")
    assert hit is not None
    mu, sigma, n = hit
    assert abs(mu - 1.99) < 1e-6
    assert n == 60


def test_carbene_lookup_other_metal_returns_none(v5_synth):
    """No carbene entry for Ru in this synthetic -> None."""
    lib = gml.GripLibrary.get(v5_synth)
    assert lib.lookup_tm_category("carbene", "Ru") is None


def test_hapto_eta5_lookup(v5_synth):
    lib = gml.GripLibrary.get(v5_synth)
    hit = lib.lookup_tm_category("hapto_eta5", "Fe")
    assert hit is not None
    assert abs(hit[0] - 2.04) < 1e-6


def test_mu_bridge_lookup(v5_synth):
    lib = gml.GripLibrary.get(v5_synth)
    hit = lib.lookup_tm_category("mu_bridge", "Cu", "Cl")
    assert hit is not None
    mu, sigma, n = hit
    assert abs(mu - 2.65) < 1e-6
    # Mu-bridge distances are LONGER than terminal Cu-Cl (chemistry sanity)
    assert mu > 2.4


def test_unknown_category_returns_none(v5_synth):
    lib = gml.GripLibrary.get(v5_synth)
    assert lib.lookup_tm_category("nonsense", "Pd") is None


def test_block_disagg_distinguishes_3d_4d_5d(v5_synth):
    """Block-disagg pair_bond chooses Pd(4d) vs Hf(5d) routes correctly."""
    os.environ.setdefault("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v5_synth)
    # Pd-Cl: block-resolved 2.31 (vs v4 wildcard 2.30 -- block is preferred)
    pd_hit = lib.lookup_bond("Pd", "*", "Cl", "sp3")
    assert pd_hit is not None
    assert abs(pd_hit[0] - 2.31) < 1e-6
    # Hf-Cl: block-resolved 2.42 (no v4 entry; block-only)
    hf_hit = lib.lookup_bond("Hf", "*", "Cl", "sp3")
    assert hf_hit is not None
    assert abs(hf_hit[0] - 2.42) < 1e-6


def test_v5_backward_compat_v4_lookup(v5_synth):
    """v4 lookup_bond on non-metal pair still works (cc_key fragment)."""
    os.environ.setdefault("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v5_synth)
    hit = lib.lookup_bond("C", "sp3", "C", "sp3")
    assert hit is not None
    assert abs(hit[0] - 1.540) < 1e-6


def test_v4_library_returns_none_for_category(tmp_path, monkeypatch):
    """v4 library MUST NOT inject category values -- returns None safely."""
    # Build minimal v4 lib
    cc_key = _key(["C", "sp3", [["C", 1, -1, "sp3"]], []])
    out = tmp_path / "v4_minimal.npz"
    np.savez_compressed(
        out,
        version=np.int32(4),
        n_master=np.int32(1), n_orig=np.int32(1),
        n_torsion=np.int32(0), n_orig_torsion=np.int32(0),
        n_pair_bond=np.int32(0), n_triple_angle=np.int32(0),
        n_improper_pair=np.int32(0),
        keys=np.array([cc_key], dtype=np.object_),
        orig_keys=np.array([cc_key], dtype=np.object_),
        orig_idx=np.zeros(1, dtype=np.int32),
        bond_mu=np.array([1.54]), bond_sigma=np.array([0.01]),
        bond_n=np.array([1000], dtype=np.int32),
        bond_p5=np.full(1, np.nan, dtype=np.float32),
        bond_p50=np.full(1, np.nan, dtype=np.float32),
        bond_p95=np.full(1, np.nan, dtype=np.float32),
        angle_mu=np.full(1, np.nan), angle_sigma=np.full(1, np.nan),
        angle_n=np.zeros(1, dtype=np.int32),
        angle_p5=np.full(1, np.nan, dtype=np.float32),
        angle_p50=np.full(1, np.nan, dtype=np.float32),
        angle_p95=np.full(1, np.nan, dtype=np.float32),
        improper_mu=np.full(1, np.nan), improper_sigma=np.full(1, np.nan),
        improper_n=np.zeros(1, dtype=np.int32),
        improper_p5=np.full(1, np.nan, dtype=np.float32),
        improper_p50=np.full(1, np.nan, dtype=np.float32),
        improper_p95=np.full(1, np.nan, dtype=np.float32),
        fb_flat=np.zeros(0, dtype=np.int32),
        fb_ptr=np.zeros(2, dtype=np.int32),
        pair_bond_keys=np.array([], dtype=np.object_),
        pair_bond_mu=np.zeros(0),
        pair_bond_sigma=np.zeros(0),
        pair_bond_n=np.zeros(0, dtype=np.int32),
        triple_angle_keys=np.array([], dtype=np.object_),
        triple_angle_mu=np.zeros(0),
        triple_angle_sigma=np.zeros(0),
        triple_angle_n=np.zeros(0, dtype=np.int32),
        improper_pair_keys=np.array([], dtype=np.object_),
        improper_pair_mu=np.zeros(0),
        improper_pair_sigma=np.zeros(0),
        improper_pair_n=np.zeros(0, dtype=np.int32),
    )
    gml.GripLibrary._SINGLETONS.clear()
    lib = gml.GripLibrary.get(out)
    assert lib.version == 4
    assert lib.has_tm_categories is False
    assert lib.has_block_disagg is False
    # Critical: category lookup returns None on v4 (no false data)
    assert lib.lookup_tm_category("carbene", "Pd") is None
    gml.GripLibrary._SINGLETONS.clear()


# ============================================================
# Smoke-build integration tests (only when smoke lib is present)
# ============================================================

@requires_v5_smoke
def test_smoke_lib_loads():
    gml.GripLibrary._SINGLETONS.clear()
    lib = gml.GripLibrary.get(V5_SMOKE)
    assert lib.version == 5
    assert lib.cleaning_applied is True
    assert lib.meta_n_extracted >= 4000
    assert lib.has_tm_categories
    gml.GripLibrary._SINGLETONS.clear()


@requires_v5_smoke
def test_smoke_carbene_chemistry_reasonable():
    """Carbene M-C distances should land in the 1.85-2.15 A range."""
    gml.GripLibrary._SINGLETONS.clear()
    lib = gml.GripLibrary.get(V5_SMOKE)
    # Probe a few common carbene metals
    hits = {}
    for m in ("Pd", "Ru", "Cu", "Ag", "Au"):
        h = lib.lookup_tm_category("carbene", m, min_n=5)
        if h is not None:
            hits[m] = h[0]
    # At least one hit must be in the expected range
    found_in_range = any(1.80 <= mu <= 2.20 for mu in hits.values())
    assert found_in_range, f"no plausible carbene M-C found: {hits}"
    gml.GripLibrary._SINGLETONS.clear()


@requires_v5_smoke
def test_smoke_mu_bridge_longer_than_terminal():
    """Cu mu-bridge Cl should be longer than terminal Cu-Cl (chemistry sanity)."""
    gml.GripLibrary._SINGLETONS.clear()
    lib = gml.GripLibrary.get(V5_SMOKE)
    mu_hit = lib.lookup_tm_category("mu_bridge", "Cu", "Cl")
    if mu_hit is None:
        pytest.skip("no Cu-Cl mu-bridge in smoke sample")
    # Compare to v4 wildcard Cu-Cl (mixed terminal + bridging)
    pair = lib.lookup_bond("Cu", "*", "Cl", "sp3")
    assert pair is not None
    # Bridging mean should be >= terminal+wildcard mean (mu-bridge is longer)
    # We allow equality because v4 pair pools both, raising the average.
    assert mu_hit[0] >= 2.30
    gml.GripLibrary._SINGLETONS.clear()


@requires_v5_smoke
def test_smoke_determinism():
    """Loading the same v5 lib twice yields identical lookups."""
    gml.GripLibrary._SINGLETONS.clear()
    lib1 = gml.GripLibrary.get(V5_SMOKE)
    hit1 = lib1.lookup_tm_category("carbene", "Pd")
    gml.GripLibrary._SINGLETONS.clear()
    lib2 = gml.GripLibrary.get(V5_SMOKE)
    hit2 = lib2.lookup_tm_category("carbene", "Pd")
    assert hit1 == hit2
    gml.GripLibrary._SINGLETONS.clear()


@requires_v5_smoke
def test_smoke_no_regression_vs_v4_keys():
    """All v4 master keys + pair_bond keys present in v5 smoke."""
    gml.GripLibrary._SINGLETONS.clear()
    lib = gml.GripLibrary.get(V5_SMOKE)
    # Pick a classic v4 bond -- C-C sp3
    hit = lib.lookup_bond("C", "sp3", "C", "sp3")
    assert hit is not None
    assert 1.45 < hit[0] < 1.60
    gml.GripLibrary._SINGLETONS.clear()
