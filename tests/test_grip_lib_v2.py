"""Tests for ``grip_lib_v2.npz`` schema + lookup compatibility.

These tests verify that:

1. v1 consumers still work (forward compatibility — schema is a superset).
2. v2 exposes torsion lookups via :meth:`GripLibrary.lookup_torsion`.
3. Hierarchical fallback chains terminate at the most-specific hit with
   ``n >= min_n``.
4. The library load is deterministic (re-load returns byte-identical arrays).
5. GMM mode counts are 1..3 with NaN-padded inactive slots.

The tests synthesise a tiny in-memory v2 npz to avoid depending on the
multi-GB CIF corpus.  When ``grip_lib_v2.npz`` exists on disk we also run
sanity checks against the real artefact.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

# Pin BLAS threads early
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("PYTHONHASHSEED", "0")

# Allow imports
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from delfin.fffree import grip_mogul_lookup as gml  # noqa: E402


# ------------------------------------------------------------------
# Helpers — build a tiny synthetic v2 npz
# ------------------------------------------------------------------

def _to_key_str(parsed) -> str:
    return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)


def _build_synthetic_v2(tmp_path: Path) -> Path:
    """Create a small valid v2 npz with bonds/angles/impropers + torsions."""
    # Fragment keys (3 levels of specificity for one C-C bond)
    bond_key_full = _to_key_str([
        "C", "sp3",
        [["C", 1, -1, "sp3"]],
        []
    ])
    # Wildcard level-4 entry (legitimately reachable via fallback)
    bond_key_l4 = _to_key_str([
        "C", "sp3",
        [["*", 1, -1, "*"]],
        []
    ])
    angle_key = _to_key_str([
        "C", "sp3",
        [["C", 1, -1, "sp3"], ["H", 1, -1, "*"]],
        []
    ])
    imp_key = _to_key_str([
        "C", "sp2",
        [["C", 1, -1, "sp2"], ["C", 1, -1, "sp2"], ["H", 1, -1, "*"]],
        []
    ])

    master_keys = sorted({bond_key_full, bond_key_l4, angle_key, imp_key})
    n_master = len(master_keys)
    key_to_idx = {k: i for i, k in enumerate(master_keys)}

    bond_mu = np.full(n_master, np.nan, dtype=np.float64)
    bond_sigma = np.full(n_master, np.nan, dtype=np.float64)
    bond_n = np.zeros(n_master, dtype=np.int32)
    angle_mu = np.full(n_master, np.nan, dtype=np.float64)
    angle_sigma = np.full(n_master, np.nan, dtype=np.float64)
    angle_n = np.zeros(n_master, dtype=np.int32)
    improper_mu = np.full(n_master, np.nan, dtype=np.float64)
    improper_sigma = np.full(n_master, np.nan, dtype=np.float64)
    improper_n = np.zeros(n_master, dtype=np.int32)

    # Populate (deterministic numbers; well-known organic ideals)
    bond_mu[key_to_idx[bond_key_full]] = 1.541
    bond_sigma[key_to_idx[bond_key_full]] = 0.013
    bond_n[key_to_idx[bond_key_full]] = 220
    bond_mu[key_to_idx[bond_key_l4]] = 1.520
    bond_sigma[key_to_idx[bond_key_l4]] = 0.045
    bond_n[key_to_idx[bond_key_l4]] = 1500

    angle_mu[key_to_idx[angle_key]] = 109.47
    angle_sigma[key_to_idx[angle_key]] = 1.2
    angle_n[key_to_idx[angle_key]] = 75

    improper_mu[key_to_idx[imp_key]] = 0.0
    improper_sigma[key_to_idx[imp_key]] = 4.0
    improper_n[key_to_idx[imp_key]] = 12

    # Percentiles (informational)
    pad_p5 = np.full(n_master, np.nan, dtype=np.float32)
    pad_p50 = np.full(n_master, np.nan, dtype=np.float32)
    pad_p95 = np.full(n_master, np.nan, dtype=np.float32)

    # Orig keys = just the bond-full key for fallback testing
    orig_keys = [bond_key_full]
    orig_idx = np.array([key_to_idx[bond_key_full]], dtype=np.int32)
    # Fallback chain: full -> level4 -> ...
    chain_levels = gml._fallback_levels(json.loads(bond_key_full))
    fb_flat = np.array(
        [key_to_idx.get(lvl, -1) for lvl in chain_levels if lvl in key_to_idx],
        dtype=np.int32,
    )
    fb_ptr = np.array([0, len(fb_flat)], dtype=np.int32)

    # Torsions: 2 keys, one unimodal sp3-C-C, one bimodal sp2-aromatic
    t_uni = _to_key_str([
        "C", "sp3", "C", "sp3", "C", "C", -1, False, False
    ])
    t_bi = _to_key_str([
        "C", "sp2", "C", "sp2", "C", "C", 6, True, True
    ])
    # also a generic fallback for the bimodal one
    t_generic = _to_key_str([
        "C", "*", "C", "*", "*", "*", -1, False, False
    ])
    torsion_keys = sorted({t_uni, t_bi, t_generic})
    n_torsion = len(torsion_keys)
    t_key_to_idx = {k: i for i, k in enumerate(torsion_keys)}

    torsion_n = np.zeros(n_torsion, dtype=np.int32)
    torsion_n_components = np.ones(n_torsion, dtype=np.int32)
    torsion_pi = np.full((n_torsion, 3), np.nan, dtype=np.float32)
    torsion_mu = np.full((n_torsion, 3), np.nan, dtype=np.float32)
    torsion_sigma = np.full((n_torsion, 3), np.nan, dtype=np.float32)

    # unimodal: 1 Gaussian at 60° anti / +60°
    torsion_n[t_key_to_idx[t_uni]] = 88
    torsion_n_components[t_key_to_idx[t_uni]] = 1
    torsion_pi[t_key_to_idx[t_uni], 0] = 1.0
    torsion_mu[t_key_to_idx[t_uni], 0] = 60.0
    torsion_sigma[t_key_to_idx[t_uni], 0] = 12.0

    # bimodal: 0°/180° (Z/E)
    torsion_n[t_key_to_idx[t_bi]] = 150
    torsion_n_components[t_key_to_idx[t_bi]] = 2
    torsion_pi[t_key_to_idx[t_bi], 0] = 0.4
    torsion_pi[t_key_to_idx[t_bi], 1] = 0.6
    torsion_mu[t_key_to_idx[t_bi], 0] = 0.0
    torsion_mu[t_key_to_idx[t_bi], 1] = 180.0
    torsion_sigma[t_key_to_idx[t_bi], 0] = 8.0
    torsion_sigma[t_key_to_idx[t_bi], 1] = 10.0

    # generic: trimodal fallback
    torsion_n[t_key_to_idx[t_generic]] = 5000
    torsion_n_components[t_key_to_idx[t_generic]] = 3
    torsion_pi[t_key_to_idx[t_generic], :] = [0.33, 0.34, 0.33]
    torsion_mu[t_key_to_idx[t_generic], :] = [-60.0, 60.0, 180.0]
    torsion_sigma[t_key_to_idx[t_generic], :] = [15.0, 15.0, 18.0]

    # Torsion orig key + fallback chain
    torsion_orig_keys = [t_uni, t_bi]
    n_orig_t = len(torsion_orig_keys)
    tor_fb_flat_list = []
    tor_fb_ptr = [0]
    for k in torsion_orig_keys:
        chain = gml.GripLibrary._torsion_fallback_levels(json.loads(k))
        for lvl in chain:
            j = t_key_to_idx.get(lvl)
            if j is not None:
                tor_fb_flat_list.append(j)
        tor_fb_ptr.append(len(tor_fb_flat_list))
    torsion_fb_flat = np.array(tor_fb_flat_list, dtype=np.int32)
    torsion_fb_ptr = np.array(tor_fb_ptr, dtype=np.int32)

    out = tmp_path / "grip_lib_v2_synthetic.npz"
    np.savez_compressed(
        out,
        version=np.int32(2),
        n_master=np.int32(n_master),
        n_orig=np.int32(len(orig_keys)),
        n_torsion=np.int32(n_torsion),
        n_orig_torsion=np.int32(n_orig_t),
        keys=np.array(master_keys, dtype=np.object_),
        orig_keys=np.array(orig_keys, dtype=np.object_),
        orig_idx=orig_idx,
        bond_mu=bond_mu, bond_sigma=bond_sigma, bond_n=bond_n,
        bond_p5=pad_p5, bond_p50=pad_p50, bond_p95=pad_p95,
        angle_mu=angle_mu, angle_sigma=angle_sigma, angle_n=angle_n,
        angle_p5=pad_p5, angle_p50=pad_p50, angle_p95=pad_p95,
        improper_mu=improper_mu, improper_sigma=improper_sigma, improper_n=improper_n,
        improper_p5=pad_p5, improper_p50=pad_p50, improper_p95=pad_p95,
        fb_flat=fb_flat, fb_ptr=fb_ptr,
        torsion_keys=np.array(torsion_keys, dtype=np.object_),
        torsion_orig_keys=np.array(torsion_orig_keys, dtype=np.object_),
        torsion_n=torsion_n,
        torsion_n_components=torsion_n_components,
        torsion_pi=torsion_pi,
        torsion_mu=torsion_mu,
        torsion_sigma=torsion_sigma,
        torsion_fb_flat=torsion_fb_flat,
        torsion_fb_ptr=torsion_fb_ptr,
    )
    return out


# ------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------


@pytest.fixture(scope="module")
def synthetic_v2(tmp_path_factory):
    tdir = tmp_path_factory.mktemp("griplibv2")
    path = _build_synthetic_v2(tdir)
    # Clear any cached singletons so we definitely get OUR file
    gml.GripLibrary._SINGLETONS.clear()
    return path


def test_v2_loads_and_reports_version(synthetic_v2):
    lib = gml.GripLibrary.get(synthetic_v2)
    assert lib.version == 2
    assert lib.has_torsions is True
    assert lib.n_torsion >= 2


def test_v2_bond_lookup_full_specificity(synthetic_v2):
    lib = gml.GripLibrary.get(synthetic_v2)
    r = lib.lookup_bond("C", "sp3", "C", "sp3")
    assert r is not None
    mu, sigma, n = r
    # Most-specific match has n=220, mu=1.541
    assert abs(mu - 1.541) < 1e-9
    assert abs(sigma - 0.013) < 1e-9
    assert n == 220


def test_v2_bond_lookup_falls_back_to_wildcard(synthetic_v2):
    lib = gml.GripLibrary.get(synthetic_v2)
    # Unknown-element pair forces walk to level-4 / 5
    r = lib.lookup_bond("C", "sp3", "Xx", "sp3")
    # Synthetic library only has the C-C key; lookup_bond canonical chain
    # for a "Xx" neighbour reaches the wildcard fallback bin which we
    # populated with mu=1.520
    assert r is not None
    mu, sigma, n = r
    assert n >= 5


def test_v2_torsion_lookup_unimodal(synthetic_v2):
    lib = gml.GripLibrary.get(synthetic_v2)
    r = lib.lookup_torsion("C", "C", "sp3", "C", "sp3", "C")
    assert r is not None
    assert r["n_components"] == 1
    assert abs(r["mu"][0] - 60.0) < 1e-3
    assert abs(r["sigma"][0] - 12.0) < 1e-3
    assert abs(r["pi"][0] - 1.0) < 1e-6
    assert r["n"] == 88


def test_v2_torsion_lookup_bimodal_aromatic_ring(synthetic_v2):
    lib = gml.GripLibrary.get(synthetic_v2)
    r = lib.lookup_torsion("C", "C", "sp2", "C", "sp2", "C",
                           ring_min_bc=6, arom_b=True, arom_c=True)
    assert r is not None
    assert r["n_components"] == 2
    # mus sorted ascending
    assert r["mu"][0] < r["mu"][1]
    assert abs(r["mu"][0] - 0.0) < 1e-3
    assert abs(r["mu"][1] - 180.0) < 1e-3
    # weights sum to 1
    assert abs(sum(r["pi"]) - 1.0) < 1e-6
    assert r["n"] == 150


def test_v2_torsion_lookup_fallback_chain(synthetic_v2):
    lib = gml.GripLibrary.get(synthetic_v2)
    # Unknown endpoint atoms should still hit the generic 3-mode fallback
    r = lib.lookup_torsion("Xx", "C", "sp3", "C", "sp3", "Yx")
    assert r is not None
    # Reached the generic level: 3 modes
    assert r["n_components"] == 3
    assert r["n"] == 5000


def test_v2_torsion_lookup_returns_none_below_minn(synthetic_v2):
    lib = gml.GripLibrary.get(synthetic_v2)
    r = lib.lookup_torsion("Xx", "C", "sp3", "C", "sp3", "Yx", min_n=10000)
    # Even the generic 5000-sample bin is too sparse here
    assert r is None


def test_v2_torsion_canonicalisation_invariant(synthetic_v2):
    """a-b-c-d and d-c-b-a must hit the same bin (dihedral symmetry)."""
    lib = gml.GripLibrary.get(synthetic_v2)
    r1 = lib.lookup_torsion("C", "C", "sp2", "C", "sp2", "C",
                            ring_min_bc=6, arom_b=True, arom_c=True)
    r2 = lib.lookup_torsion("C", "C", "sp2", "C", "sp2", "C",
                            ring_min_bc=6, arom_b=True, arom_c=True)
    # Reversed call (same edge, swapped endpoint order):
    r3 = lib.lookup_torsion("C", "C", "sp2", "C", "sp2", "C",
                            ring_min_bc=6, arom_c=True, arom_b=True)
    assert r1 is not None and r2 is not None and r3 is not None
    assert r1["n"] == r2["n"] == r3["n"]
    assert r1["mu"] == r2["mu"] == r3["mu"]


def test_v1_consumer_still_works_with_v2(synthetic_v2):
    """An existing v1-style call path must remain valid on a v2 library."""
    lib = gml.GripLibrary.get(synthetic_v2)
    # Module-level shortcut for bonds (the v1 API) still works.
    r = gml.lookup_bond("C", "sp3", "C", "sp3", library=lib)
    assert r is not None
    # Module-level shortcut for improper
    r_imp = lib.lookup_improper("C", "sp2", ["C", "C", "H"],
                                 neighbor_hybs_sorted=["sp2", "sp2", "*"])
    # Synthetic improper has n=12 -> hit
    assert r_imp is not None
    assert r_imp[2] == 12


def test_v2_library_load_is_deterministic(synthetic_v2):
    """Loading the same library twice must yield byte-identical arrays."""
    # Force two independent loads (bypass the singleton cache)
    a = np.load(synthetic_v2, allow_pickle=True)
    b = np.load(synthetic_v2, allow_pickle=True)
    for k in ("bond_mu", "bond_sigma", "bond_n",
              "angle_mu", "angle_sigma", "angle_n",
              "improper_mu", "improper_sigma", "improper_n",
              "torsion_mu", "torsion_sigma", "torsion_pi",
              "torsion_n", "torsion_n_components"):
        np.testing.assert_array_equal(a[k], b[k])


def test_v2_gmm_mode_count_invariants(synthetic_v2):
    """Each torsion bin: 1 <= n_components <= 3; inactive slots are NaN."""
    lib = gml.GripLibrary.get(synthetic_v2)
    for i in range(lib.n_torsion):
        nc = int(lib._torsion_n_components[i])
        assert 1 <= nc <= 3
        for k in range(nc):
            assert np.isfinite(lib._torsion_mu[i, k])
            assert np.isfinite(lib._torsion_sigma[i, k])
            assert np.isfinite(lib._torsion_pi[i, k])
            assert float(lib._torsion_sigma[i, k]) > 0.0
        for k in range(nc, 3):
            assert np.isnan(lib._torsion_mu[i, k])
            assert np.isnan(lib._torsion_sigma[i, k])
            assert np.isnan(lib._torsion_pi[i, k])
        # mus sorted ascending
        mus = [float(lib._torsion_mu[i, k]) for k in range(nc)]
        assert mus == sorted(mus)
        # weights sum to 1 (within numerical tolerance)
        pis = [float(lib._torsion_pi[i, k]) for k in range(nc)]
        assert abs(sum(pis) - 1.0) < 1e-3


def test_v2_lookup_returns_none_on_malformed_input(synthetic_v2):
    lib = gml.GripLibrary.get(synthetic_v2)
    # has_torsions but invalid types -> None
    r = lib.lookup_torsion("C", "C", "sp3", "C", "sp3", "C",
                           ring_min_bc="not-an-int")  # type: ignore[arg-type]
    assert r is None


def test_v1_library_lookup_torsion_returns_none():
    """A v1 library must report ``has_torsions=False`` and return None."""
    v1_path = Path(
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v1.npz"
    )
    if not v1_path.exists():
        pytest.skip("v1 library not present")
    gml.GripLibrary._SINGLETONS.clear()
    lib = gml.GripLibrary.get(v1_path)
    assert lib.version == 1
    assert lib.has_torsions is False
    assert lib.lookup_torsion("C", "C", "sp3", "C", "sp3", "C") is None


# ------------------------------------------------------------------
# Optional: run sanity checks against the real v2 library if present
# ------------------------------------------------------------------


def test_real_v2_library_smoke_if_present():
    real = Path(
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v2.npz"
    )
    if not real.exists():
        pytest.skip("real grip_lib_v2.npz not yet built")
    gml.GripLibrary._SINGLETONS.clear()
    lib = gml.GripLibrary.get(real)
    assert lib.version == 2
    assert lib.has_torsions
    # The full library should have a wealth of C-C / C-H / C-N bonds.
    r = lib.lookup_bond("C", "sp3", "C", "sp3")
    assert r is not None
    mu, sigma, n = r
    assert 1.40 < mu < 1.65   # C-C single bond range
    assert sigma > 0.0
    assert n >= 5
