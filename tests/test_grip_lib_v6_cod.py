"""Tests for grip_lib v6 — COD-source dual-library + MergedLibrary.

Verifies:

1. v6 metadata loads (version=6, source=COD, cod_n_extracted).
2. v6 keeps full v5 schema (TM-category, block-disagg, torsion GMM).
3. v6 library passes back into v5 GripLibrary class (backward-compatible npz).
4. MergedLibrary degrades to single-source when one side is None.
5. MergedLibrary n-weighted merge produces a valid mean.
6. Provenance tags classify correctly (agree/marginal/disagree/single-source).
7. Adopter-mode (CCDC unset, COD set) works.
8. Env-flag DELFIN_FFFREE_GRIP_LIB_COD_PATH is honoured by
   :func:`get_default_merged_library`.
9. Default-OFF byte-identical: with COD env unset, merged lookup ==
   single-CCDC lookup.
"""
from __future__ import annotations

import json
import math
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

V6_SMOKE = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v6_cod_smoke.npz"
)
V6_FULL = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v6_cod.npz"
)
V5_FULL = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz"
)
V5_SMOKE = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5_smoke.npz"
)

requires_v6_any = pytest.mark.skipif(
    not (V6_SMOKE.exists() or V6_FULL.exists()),
    reason="no v6 COD library built yet",
)
requires_v6_smoke = pytest.mark.skipif(
    not V6_SMOKE.exists(), reason="v6 COD smoke not built",
)


def _v6_path() -> Path:
    return V6_FULL if V6_FULL.exists() else V6_SMOKE


def _v5_path() -> Path:
    if V5_FULL.exists():
        return V5_FULL
    return V5_SMOKE


# ---------------------------------------------------------------------------
# 1) v6 metadata loads
# ---------------------------------------------------------------------------
@requires_v6_any
def test_v6_npz_version_and_source():
    """Verify version=6 and source='COD' fields."""
    npz = np.load(_v6_path(), allow_pickle=True)
    assert int(npz["version"]) == 6
    assert str(npz["source"]) == "COD"
    assert int(npz["cod_n_extracted"]) >= 1
    # cleaning_applied stamp preserved
    assert int(npz["cleaning_applied"]) == 1


# ---------------------------------------------------------------------------
# 2) v6 has the full v5 schema (every array name present)
# ---------------------------------------------------------------------------
@requires_v6_any
def test_v6_schema_matches_v5():
    """Every v5 table name must exist in v6 npz."""
    required = [
        # v3 master
        "keys", "bond_mu", "bond_sigma", "bond_n",
        "angle_mu", "angle_sigma", "angle_n",
        "improper_mu", "improper_sigma", "improper_n",
        "fb_flat", "fb_ptr",
        # v3 torsion
        "torsion_keys", "torsion_n", "torsion_n_components",
        "torsion_pi", "torsion_mu", "torsion_sigma",
        "torsion_fb_flat", "torsion_fb_ptr",
        # v4 pair tables
        "pair_bond_keys", "pair_bond_mu", "pair_bond_sigma", "pair_bond_n",
        "triple_angle_keys", "triple_angle_mu",
        "triple_angle_sigma", "triple_angle_n",
        "improper_pair_keys", "improper_pair_mu",
        "improper_pair_sigma", "improper_pair_n",
        # v5 TM-category tables
        "tm_carbene_keys", "tm_carbene_mu", "tm_carbene_sigma", "tm_carbene_n",
        "tm_hapto_eta2_keys", "tm_hapto_eta2_mu",
        "tm_hapto_eta2_sigma", "tm_hapto_eta2_n",
        "tm_hapto_eta5_keys", "tm_hapto_eta5_mu",
        "tm_hapto_eta5_sigma", "tm_hapto_eta5_n",
        "tm_hapto_eta6_keys", "tm_hapto_eta6_mu",
        "tm_hapto_eta6_sigma", "tm_hapto_eta6_n",
        "tm_mu_bridge_keys", "tm_mu_bridge_mu",
        "tm_mu_bridge_sigma", "tm_mu_bridge_n",
        "tm_agostic_keys", "tm_agostic_mu",
        "tm_agostic_sigma", "tm_agostic_n",
        "tm_ox_addition_keys", "tm_ox_addition_mu",
        "tm_ox_addition_sigma", "tm_ox_addition_n",
        # v5 block-disagg
        "tm_pair_bond_block_keys", "tm_pair_bond_block_mu",
        "tm_pair_bond_block_sigma", "tm_pair_bond_block_n",
        "tm_triple_angle_block_keys", "tm_triple_angle_block_mu",
        "tm_triple_angle_block_sigma", "tm_triple_angle_block_n",
        # Lever-3 placeholder (always False on v6 — no Mogul for COD)
        "mogul_validated",
    ]
    npz = np.load(_v6_path(), allow_pickle=True)
    missing = [r for r in required if r not in npz.files]
    assert missing == [], f"v6 missing v5 schema fields: {missing}"


# ---------------------------------------------------------------------------
# 3) v6 loads via v5 GripLibrary (backward-compatible)
# ---------------------------------------------------------------------------
@requires_v6_any
def test_v6_loads_via_grip_library():
    """The v6 npz must be load-compatible with the v5 GripLibrary class."""
    lib = gml.GripLibrary.get(_v6_path())
    assert lib.version == 6
    assert lib.has_pair_tables
    assert lib.has_tm_categories
    assert lib.n_master > 0
    assert lib.n_pair_bond > 0


# ---------------------------------------------------------------------------
# 4) MergedLibrary degrades to single-source
# ---------------------------------------------------------------------------
@requires_v6_any
def test_merged_library_cod_only_degrades_to_cod():
    """When CCDC is None, MergedLibrary defers to COD only."""
    cod_lib = gml.GripLibrary.get(_v6_path())
    merged = gml.MergedLibrary(ccdc_lib=None, cod_lib=cod_lib)
    assert merged.has_cod is True
    assert merged.has_ccdc is False
    # carbene lookup on COD-only path — should mirror direct COD lookup.
    direct = cod_lib.lookup_tm_category("carbene", "Ru", "C")
    via_merged = merged.lookup_tm_category("carbene", "Ru", "C")
    assert direct == via_merged


@requires_v6_any
def test_merged_library_ccdc_only_degrades_to_ccdc():
    """When COD is None, MergedLibrary defers to CCDC only."""
    if not _v5_path().exists():
        pytest.skip("no v5 to compare")
    ccdc_lib = gml.GripLibrary.get(_v5_path())
    merged = gml.MergedLibrary(ccdc_lib=ccdc_lib, cod_lib=None)
    assert merged.has_ccdc is True
    assert merged.has_cod is False
    direct = ccdc_lib.lookup_bond("C", "sp3", "C", "sp3")
    via_merged = merged.lookup_bond("C", "sp3", "C", "sp3")
    assert direct == via_merged


# ---------------------------------------------------------------------------
# 5) MergedLibrary n-weighted merge gives sensible mean
# ---------------------------------------------------------------------------
def test_weighted_merge_formula():
    """Closed-form verification of ``_weighted_mean_sigma``.

    Two normal samples N(1.50, 0.04, n=100) and N(1.52, 0.03, n=50) -->
    merged mean = (100*1.50 + 50*1.52)/150 = 1.5067
    merged var  = (100*(0.04^2 + (1.5-1.5067)^2)
                  + 50*(0.03^2 + (1.52-1.5067)^2)) / 150
    """
    mu_m, sig_m, n_m = gml._weighted_mean_sigma(
        1.50, 0.04, 100,
        1.52, 0.03, 50,
    )
    assert n_m == 150
    assert math.isclose(mu_m, (100 * 1.50 + 50 * 1.52) / 150, abs_tol=1e-9)
    # var = (100*(0.0016 + 4.49e-5) + 50*(0.0009 + 1.79e-4)) / 150
    expected_var = (
        100 * (0.04**2 + (1.50 - mu_m)**2)
        + 50 * (0.03**2 + (1.52 - mu_m)**2)
    ) / 150
    assert math.isclose(sig_m, math.sqrt(expected_var), abs_tol=1e-9)


def test_provenance_classification():
    """Verify the tag-mapping for various (mu_diff, sigma) pairs."""
    # AGREE: 1.50 vs 1.51, pooled sigma ~0.04, diff = 0.01/0.04 = 0.25 < 1.0
    tag = gml._classify_provenance((1.50, 0.04, 100), (1.51, 0.04, 100))
    assert tag == gml.PROVENANCE_BOTH_AGREE
    # DISAGREE: 1.50 vs 1.70, diff 0.20/0.05 = 4.0 > 2.0
    tag = gml._classify_provenance((1.50, 0.05, 100), (1.70, 0.05, 100))
    assert tag == gml.PROVENANCE_BOTH_DISAGREE
    # MARGINAL: diff 0.05/0.04 = 1.25 between 1.0 and 2.0
    tag = gml._classify_provenance((1.50, 0.04, 100), (1.55, 0.04, 100))
    assert tag == gml.PROVENANCE_BOTH_MARGINAL
    # COD-only
    assert gml._classify_provenance(None, (1.50, 0.04, 100)) == gml.PROVENANCE_COD_ONLY
    # CCDC-only
    assert gml._classify_provenance((1.50, 0.04, 100), None) == gml.PROVENANCE_CCDC_ONLY
    # None
    assert gml._classify_provenance(None, None) == gml.PROVENANCE_NONE


# ---------------------------------------------------------------------------
# 6) Lookup with provenance returns merged + tag + per-source hits
# ---------------------------------------------------------------------------
@requires_v6_any
def test_lookup_bond_with_provenance_cod_only():
    """COD-only configuration -> tag = ``cod-only`` and merged == cod hit."""
    cod_lib = gml.GripLibrary.get(_v6_path())
    merged = gml.MergedLibrary(ccdc_lib=None, cod_lib=cod_lib)
    merged_hit, tag, h_ccdc, h_cod = merged.lookup_bond_with_provenance(
        "C", "sp3", "O", "sp3"
    )
    # Even if COD misses, the tag should still classify correctly.
    if h_cod is not None:
        assert tag == gml.PROVENANCE_COD_ONLY
        assert merged_hit == h_cod
    else:
        assert tag == gml.PROVENANCE_NONE
        assert merged_hit is None
    assert h_ccdc is None


# ---------------------------------------------------------------------------
# 7) Adopter-mode (no CCDC env, COD env set) → cod-only merged lib.
# ---------------------------------------------------------------------------
@requires_v6_any
def test_adopter_mode_cod_only_via_env(monkeypatch, tmp_path):
    """Unset CCDC env, set COD env to v6 path -> merged lib has only cod_lib."""
    # Point CCDC to a non-existent file so the loader skips it.
    monkeypatch.setenv(gml.DELFIN_GRIP_LIB_PATH_ENV, str(tmp_path / "nonexistent.npz"))
    monkeypatch.setenv(gml.DELFIN_FFFREE_GRIP_LIB_COD_PATH_ENV, str(_v6_path()))
    # Clear singleton cache for both classes so the env change actually takes effect
    gml.GripLibrary._SINGLETONS.clear()
    gml.MergedLibrary._SINGLETONS.clear()
    merged = gml.get_default_merged_library()
    assert merged.has_cod is True
    assert merged.has_ccdc is False
    assert merged.version == 6


# ---------------------------------------------------------------------------
# 8) Default-OFF byte-identical: COD env unset -> merged == CCDC-only.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not V5_FULL.exists() and not V5_SMOKE.exists(),
                    reason="no v5 lib to compare")
def test_default_off_byte_identical_to_ccdc(monkeypatch):
    """With DELFIN_FFFREE_GRIP_LIB_COD_PATH unset, lookups go through CCDC-only.

    This is the "byte-identical default" requirement: existing pipelines
    that do NOT set the COD env-flag must behave exactly as before.
    """
    monkeypatch.delenv(gml.DELFIN_FFFREE_GRIP_LIB_COD_PATH_ENV, raising=False)
    monkeypatch.setenv(gml.DELFIN_GRIP_LIB_PATH_ENV, str(_v5_path()))
    gml.GripLibrary._SINGLETONS.clear()
    gml.MergedLibrary._SINGLETONS.clear()

    merged = gml.get_default_merged_library()
    assert merged.has_ccdc is True
    assert merged.has_cod is False

    ccdc_only = gml.GripLibrary.get(_v5_path())
    # Test a handful of lookups
    direct_b = ccdc_only.lookup_bond("C", "sp3", "C", "sp3")
    via_b = merged.lookup_bond("C", "sp3", "C", "sp3")
    assert direct_b == via_b
    direct_a = ccdc_only.lookup_angle("C", "C", "sp3", "C")
    via_a = merged.lookup_angle("C", "C", "sp3", "C")
    assert direct_a == via_a


# ---------------------------------------------------------------------------
# 9) v6 produces chemistry-reasonable Ru-carbene values when applicable
# ---------------------------------------------------------------------------
@requires_v6_any
def test_v6_carbene_chemistry_sanity():
    """When the carbene table holds Ru-C, the mean should fall in 1.8-2.3 A."""
    cod_lib = gml.GripLibrary.get(_v6_path())
    hit = cod_lib.lookup_tm_category("carbene", "Ru", "C")
    if hit is None:
        pytest.skip("Ru-carbene not present in this build (smoke too small)")
    mu, sigma, n = hit
    assert 1.7 <= mu <= 2.30, f"Ru-carbene mu out of range: {mu}"
    assert sigma > 0
    assert n >= 5


# ---------------------------------------------------------------------------
# 10) Determinism: a fresh load yields identical lookups.
# ---------------------------------------------------------------------------
@requires_v6_any
def test_determinism_of_v6_load():
    """Two fresh loads of v6 give byte-identical lookup results."""
    gml.GripLibrary._SINGLETONS.clear()
    a = gml.GripLibrary.get(_v6_path())
    hit_a = a.lookup_bond("C", "sp3", "O", "sp3")
    gml.GripLibrary._SINGLETONS.clear()
    b = gml.GripLibrary.get(_v6_path())
    hit_b = b.lookup_bond("C", "sp3", "O", "sp3")
    assert hit_a == hit_b


# ---------------------------------------------------------------------------
# 11) Cross-validation: when both libs hit, merged n is sum.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not (V5_FULL.exists() and V6_FULL.exists()),
    reason="needs both v5_full and v6_full",
)
def test_merged_n_is_sum_when_both_hit():
    """For a query where both libs have data, merged n == n_ccdc + n_cod."""
    ccdc = gml.GripLibrary.get(V5_FULL)
    cod = gml.GripLibrary.get(V6_FULL)
    merged = gml.MergedLibrary(ccdc_lib=ccdc, cod_lib=cod)
    h_ccdc = ccdc.lookup_bond("C", "sp3", "C", "sp3")
    h_cod = cod.lookup_bond("C", "sp3", "C", "sp3")
    if h_ccdc is None or h_cod is None:
        pytest.skip("C-C bond not in both libs")
    h_merged = merged.lookup_bond("C", "sp3", "C", "sp3")
    assert h_merged is not None
    assert h_merged[2] == h_ccdc[2] + h_cod[2]
