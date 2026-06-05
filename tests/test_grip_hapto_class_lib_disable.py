"""Tests for the class-conditional hapto/carbene library DISABLE guard.

Mission B3 (2026-06-05) — preemptive heal for the v4-lib activation
verdict regression (hapto_count +74 %, hapto_geom +74 %, cshm_mean_max
+23 % in smoke 550).  When
``DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE=1`` is set, queries that look
like metal-hapto-π / metal-carbene chemistry MUST return ``None`` at
the lookup entry, BEFORE any v3/v4/v5/v6 walk.  σ-donor queries
(M-N, M-O, M-C sp3) must still hit normally.

The tests cover:

  * default-OFF byte-identical contract (no env flag set -> guard is a
    no-op),
  * heuristic class detection on pure (Z, hyb) input (no library
    needed),
  * GripLibrary.lookup_bond / lookup_angle / lookup_improper
    suppression for hapto/carbene queries when the flag is set,
  * σ-donor queries survive the suppression,
  * v3 library + flag = legacy chain behaviour preserved (the guard
    fires BEFORE the library walk, so v3-only hits for σ-classes still
    work),
  * MergedLibrary inherits the suppression naturally (it delegates
    per-source),
  * determinism: two identical calls produce identical outputs (same
    seed / no float drift),
  * env aliases: ``"1"``, ``"true"``, ``"yes"``, ``"on"`` toggle ON;
    ``""``, ``"0"``, ``"false"``, ``"no"`` toggle OFF.
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


# ---------------------------------------------------------------------------
# Helper: build a tiny synthetic v4 library so we can exercise the lookup
# path WITHOUT the 88 MB real-CCDC artefact.
# ---------------------------------------------------------------------------
def _key(parsed) -> str:
    return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)


def _build_v4_synth(tmp: Path) -> Path:
    # σ-class baseline: C(sp3)-C(sp3) bond and C(sp3)-O(sp3) angle that the
    # legacy chain hits straight away.
    cc_key = _key(["C", "sp3", [["C", 1, -1, "sp3"]], []])
    no_angle_key = _key(["N", "sp3", [["C", 1, -1, "sp3"], ["O", 1, -1, "sp3"]], []])
    master_keys = [cc_key, no_angle_key]
    n_master = len(master_keys)

    bond_mu = np.array([1.540, np.nan], dtype=np.float64)
    bond_sigma = np.array([0.010, np.nan], dtype=np.float64)
    bond_n = np.array([1000, 0], dtype=np.int32)
    zero_p = np.full(n_master, np.nan, dtype=np.float32)
    angle_mu = np.array([np.nan, 110.0], dtype=np.float64)
    angle_sigma = np.array([np.nan, 3.0], dtype=np.float64)
    angle_n = np.array([0, 200], dtype=np.int32)
    improper_mu = np.full(n_master, np.nan, dtype=np.float64)
    improper_sigma = np.full(n_master, np.nan, dtype=np.float64)
    improper_n = np.zeros(n_master, dtype=np.int32)

    # v4 pair tables -- these are the WILDCARD pool the guard exists to
    # bypass.  Stored values are deliberately "tempting" so a test that
    # forgets the guard fails loudly.
    irc_sp2_key = _key(["C", "sp2", "Ir", "*"])     # hapto-π / carbene-like
    irc_sp3_key = _key(["C", "sp3", "Ir", "*"])     # σ-alkyl
    irn_key = _key(["Ir", "*", "N", "sp3"])         # σ-N-donor
    pair_bond_keys = [irc_sp2_key, irc_sp3_key, irn_key]
    pair_bond_mu = np.array([2.005, 2.110, 2.150], dtype=np.float64)
    pair_bond_sigma = np.array([0.030, 0.020, 0.015], dtype=np.float64)
    pair_bond_n = np.array([150, 120, 300], dtype=np.int32)

    # Triple-angle table -- C-Ir-C (centered at Ir).
    ang_circ_key = _key(["Ir", "*", "C", "C"])      # M-C-...-C wildcard pool
    # Add a "good" σ-N angle to verify it survives suppression.
    ang_sigma_key = _key(["Ir", "*", "N", "N"])     # σ-only
    triple_angle_keys = [ang_circ_key, ang_sigma_key]
    triple_angle_mu = np.array([89.5, 95.0], dtype=np.float64)
    triple_angle_sigma = np.array([3.5, 4.0], dtype=np.float64)
    triple_angle_n = np.array([80, 60], dtype=np.int32)

    # Improper at Ir w/ 3 C neighbours
    imp_key = _key(["Ir", "*", sorted(["C", "C", "C"])])
    improper_pair_keys = [imp_key]
    improper_pair_mu = np.array([0.5], dtype=np.float64)
    improper_pair_sigma = np.array([4.0], dtype=np.float64)
    improper_pair_n = np.array([60], dtype=np.int32)

    out = tmp / "grip_lib_v4_b3_synth.npz"
    np.savez_compressed(
        out,
        version=np.int32(4),
        n_master=np.int32(n_master),
        n_orig=np.int32(n_master),
        n_torsion=np.int32(0),
        n_orig_torsion=np.int32(0),
        n_pair_bond=np.int32(len(pair_bond_keys)),
        n_triple_angle=np.int32(len(triple_angle_keys)),
        n_improper_pair=np.int32(len(improper_pair_keys)),
        keys=np.array(master_keys, dtype=np.object_),
        orig_keys=np.array(master_keys, dtype=np.object_),
        orig_idx=np.arange(n_master, dtype=np.int32),
        bond_mu=bond_mu, bond_sigma=bond_sigma, bond_n=bond_n,
        bond_p5=zero_p, bond_p50=zero_p, bond_p95=zero_p,
        angle_mu=angle_mu, angle_sigma=angle_sigma, angle_n=angle_n,
        angle_p5=zero_p, angle_p50=zero_p, angle_p95=zero_p,
        improper_mu=improper_mu, improper_sigma=improper_sigma,
        improper_n=improper_n,
        improper_p5=zero_p, improper_p50=zero_p, improper_p95=zero_p,
        fb_flat=np.zeros(0, dtype=np.int32),
        fb_ptr=np.zeros(n_master + 1, dtype=np.int32),
        torsion_keys=np.array([], dtype=np.object_),
        torsion_orig_keys=np.array([], dtype=np.object_),
        torsion_n=np.zeros(0, dtype=np.int32),
        torsion_n_components=np.zeros(0, dtype=np.int32),
        torsion_pi=np.zeros((0, 3), dtype=np.float32),
        torsion_mu=np.zeros((0, 3), dtype=np.float32),
        torsion_sigma=np.zeros((0, 3), dtype=np.float32),
        torsion_fb_flat=np.zeros(0, dtype=np.int32),
        torsion_fb_ptr=np.zeros(1, dtype=np.int32),
        pair_bond_keys=np.array(pair_bond_keys, dtype=np.object_),
        pair_bond_mu=pair_bond_mu,
        pair_bond_sigma=pair_bond_sigma,
        pair_bond_n=pair_bond_n,
        triple_angle_keys=np.array(triple_angle_keys, dtype=np.object_),
        triple_angle_mu=triple_angle_mu,
        triple_angle_sigma=triple_angle_sigma,
        triple_angle_n=triple_angle_n,
        improper_pair_keys=np.array(improper_pair_keys, dtype=np.object_),
        improper_pair_mu=improper_pair_mu,
        improper_pair_sigma=improper_pair_sigma,
        improper_pair_n=improper_pair_n,
    )
    return out


@pytest.fixture
def v4_b3_synth(tmp_path):
    """Yield a fresh v4 synthetic library with the singleton cache cleared."""
    p = _build_v4_synth(tmp_path)
    gml.GripLibrary._SINGLETONS.clear()
    yield p
    gml.GripLibrary._SINGLETONS.clear()


# ---------------------------------------------------------------------------
# 1) Helper-level heuristic tests (no library, pure (Z, hyb) input)
# ---------------------------------------------------------------------------
def test_helper_bond_query_metal_pi_carbon_detected(monkeypatch):
    """Metal-C(sp2) and metal-C(sp) -> True (hapto-π / carbene)."""
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raising=False)
    f = gml._is_hapto_or_carbene_bond_query
    # Both orientations of (TM, C-sp2)
    assert f("Ir", "*", "C", "sp2") is True
    assert f("C", "sp2", "Ir", "*") is True
    # (TM, C-sp) -> carbene-like
    assert f("Pt", "*", "C", "sp") is True
    assert f("C", "sp", "Pd", "*") is True
    # Several TMs from each row
    for tm in ("Fe", "Ru", "Os", "Pd", "Pt", "Ir", "Rh", "Mo", "W", "U"):
        assert f(tm, "*", "C", "sp2") is True


def test_helper_bond_query_sigma_classes_not_suppressed(monkeypatch):
    """M-C(sp3) σ-alkyl, M-N, M-O, M-P, M-Cl -> False (KEEP lookup)."""
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raising=False)
    f = gml._is_hapto_or_carbene_bond_query
    # σ-alkyl M-CH3
    assert f("Ir", "*", "C", "sp3") is False
    assert f("C", "sp3", "Ir", "*") is False
    # σ-donors -- nitrogen, oxygen, phosphorus, sulfur, halides
    for partner in ("N", "O", "P", "S", "Cl", "Br", "F"):
        assert f("Pd", "*", partner, "sp3") is False
        assert f(partner, "sp3", "Pd", "*") is False
    # Non-metal pair (organic only) -- always False
    assert f("C", "sp2", "C", "sp2") is False
    assert f("C", "sp2", "N", "sp2") is False


def test_helper_angle_query_metal_centered_with_c_neighbour(monkeypatch):
    """M-?-C wildcard pool aggregates classes -> suppress when center is TM
    and at least one neighbour is C."""
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raising=False)
    f = gml._is_hapto_or_carbene_angle_query
    # C-Ir-C : center TM, both neighbours C -> SUPPRESS
    assert f("C", "Ir", "*", "C") is True
    # C-Pd-N : center TM, one neighbour C -> SUPPRESS (mixed pool)
    assert f("C", "Pd", "*", "N") is True
    # N-Ir-N : center TM, NO C -> KEEP
    assert f("N", "Ir", "*", "N") is False
    # O-Pd-N : center TM, NO C -> KEEP
    assert f("O", "Pd", "*", "N") is False


def test_helper_angle_query_pi_carbon_with_metal_neighbour(monkeypatch):
    """sp2-C centered angle with TM neighbour (η-C-X) -> SUPPRESS."""
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raising=False)
    f = gml._is_hapto_or_carbene_angle_query
    # Ir-C(sp2)-H (hapto-π H atom)
    assert f("Ir", "C", "sp2", "H") is True
    # H-C(sp2)-Ir reverse orientation
    assert f("H", "C", "sp2", "Ir") is True
    # C(sp3)-CH2-Ir (σ-alkyl) -> KEEP
    assert f("Ir", "C", "sp3", "H") is False
    # C(sp2)-N(sp2)-H (organic only, no metal) -> KEEP
    assert f("C", "N", "sp2", "H") is False


def test_helper_improper_metal_or_pi_centered(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raising=False)
    f = gml._is_hapto_or_carbene_improper_query
    # TM center + any C neighbour -> SUPPRESS
    assert f("Ir", "*", ["C", "C", "C"]) is True
    assert f("Pd", "*", ["N", "N", "C"]) is True
    # TM center, no C neighbour -> KEEP
    assert f("Pd", "*", ["N", "N", "O"]) is False
    # sp2-C center + TM neighbour -> SUPPRESS (η-C-X-Y improper)
    assert f("C", "sp2", ["Ir", "H", "H"]) is True
    # sp3-C center + TM neighbour -> KEEP (σ-alkyl)
    assert f("C", "sp3", ["Ir", "H", "H"]) is False


# ---------------------------------------------------------------------------
# 2) Env-flag truthiness contract
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw,expected",
    [
        # OFF cases
        ("", False),
        ("0", False),
        ("false", False),
        ("False", False),
        ("FALSE", False),
        ("no", False),
        ("No", False),
        # ON cases
        ("1", True),
        ("true", True),
        ("True", True),
        ("TRUE", True),
        ("yes", True),
        ("YES", True),
        ("on", True),
        ("On", True),
    ],
)
def test_env_flag_truthiness(monkeypatch, raw, expected):
    if raw == "":
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raising=False)
    else:
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raw)
    assert gml._hapto_lib_disable_enabled() is expected


# ---------------------------------------------------------------------------
# 3) GripLibrary lookup-level integration tests on the synthetic v4 npz
# ---------------------------------------------------------------------------
def test_default_off_byte_identical_to_pre_fix(v4_b3_synth, monkeypatch):
    """Without the env flag set, lookup_bond / lookup_angle /
    lookup_improper behave EXACTLY as before -- including the v4
    pair-table hits for metal-organic queries."""
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raising=False)
    # Activate the TM-aware fallback so the v4 pair table is queryable.
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_b3_synth)

    # v4 pair-resolved hapto/carbene bond hit (would be suppressed under
    # the new guard).
    hit_pi = lib.lookup_bond("Ir", "*", "C", "sp2")
    assert hit_pi is not None
    assert abs(hit_pi[0] - 2.005) < 1e-6
    # σ-alkyl
    hit_sigma = lib.lookup_bond("Ir", "*", "C", "sp3")
    assert hit_sigma is not None
    assert abs(hit_sigma[0] - 2.110) < 1e-6
    # σ-N (always allowed)
    hit_n = lib.lookup_bond("Ir", "*", "N", "sp3")
    assert hit_n is not None
    # σ-only angle
    hit_ang_sigma = lib.lookup_angle("N", "Ir", "*", "N")
    assert hit_ang_sigma is not None


def test_hapto_disable_suppresses_metal_pi_carbon_bond(v4_b3_synth, monkeypatch):
    """Flag ON: M-C(sp2) bond returns ``None`` BEFORE any walk."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_b3_synth)

    assert lib.lookup_bond("Ir", "*", "C", "sp2") is None
    # Reverse orientation also blocked
    assert lib.lookup_bond("C", "sp2", "Ir", "*") is None
    # Carbene C-sp also blocked
    assert lib.lookup_bond("Pt", "*", "C", "sp") is None


def test_hapto_disable_preserves_sigma_classes(v4_b3_synth, monkeypatch):
    """Flag ON: σ-classes must still hit the library."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_b3_synth)

    # σ-alkyl M-C(sp3) -- still hits (the σ-CH3 distance is real chemistry).
    hit_sigma_alkyl = lib.lookup_bond("Ir", "*", "C", "sp3")
    assert hit_sigma_alkyl is not None
    assert abs(hit_sigma_alkyl[0] - 2.110) < 1e-6
    # M-N σ-donor -- still hits
    hit_n = lib.lookup_bond("Ir", "*", "N", "sp3")
    assert hit_n is not None
    assert abs(hit_n[0] - 2.150) < 1e-6
    # Pure organic C(sp3)-C(sp3) -- still hits (legacy chain)
    hit_cc = lib.lookup_bond("C", "sp3", "C", "sp3")
    assert hit_cc is not None
    assert abs(hit_cc[0] - 1.540) < 1e-6


def test_hapto_disable_suppresses_metal_centered_angle_with_c(
    v4_b3_synth, monkeypatch,
):
    """Flag ON: C-M-C angle from the v4 wildcard pool returns ``None``."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_b3_synth)

    # C-Ir-C suppressed (TM center + any C neighbour)
    assert lib.lookup_angle("C", "Ir", "*", "C") is None
    # σ-only angle N-Ir-N still hits
    hit_n = lib.lookup_angle("N", "Ir", "*", "N")
    assert hit_n is not None


def test_hapto_disable_suppresses_metal_centered_improper(
    v4_b3_synth, monkeypatch,
):
    """Flag ON: M-centered improper with C neighbours returns ``None``."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_b3_synth)

    assert lib.lookup_improper("Ir", "*", ["C", "C", "C"]) is None


def test_module_level_lookup_honours_flag(v4_b3_synth, monkeypatch):
    """The module-level shortcut funcs (used by grip_fragment_detect) also
    honour the flag because they delegate to GripLibrary.lookup_*."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_b3_synth)
    # Pass library= so we don't depend on DEFAULT_LIB_PATH being present.
    assert gml.lookup_bond("Ir", "*", "C", "sp2", library=lib) is None
    assert gml.lookup_angle("C", "Ir", "*", "C", library=lib) is None
    # σ-class still hits.
    hit = gml.lookup_bond("C", "sp3", "C", "sp3", library=lib)
    assert hit is not None
    assert abs(hit[0] - 1.540) < 1e-6


def test_hapto_disable_determinism_two_runs(v4_b3_synth, monkeypatch):
    """Two identical calls under the same flag setting must produce
    identical outputs (no float drift, no RNG)."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_b3_synth)
    # Suppressed query: both runs None.
    r1 = lib.lookup_bond("Ir", "*", "C", "sp2")
    r2 = lib.lookup_bond("Ir", "*", "C", "sp2")
    assert r1 is None and r2 is None
    # Surviving σ-query: same tuple (mu, sigma, n) in both runs.
    r3 = lib.lookup_bond("Ir", "*", "N", "sp3")
    r4 = lib.lookup_bond("Ir", "*", "N", "sp3")
    assert r3 == r4
    assert r3 is not None


def test_merged_library_inherits_suppression(v4_b3_synth, monkeypatch):
    """When MergedLibrary delegates to a sub-lib whose lookup is guarded,
    the merged result is also ``None`` for a suppressed query."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    # Build a MergedLibrary that wraps the synthetic v4 lib in both slots.
    gml.MergedLibrary._SINGLETONS.clear()
    merged = gml.MergedLibrary.get(v4_b3_synth, v4_b3_synth)
    assert merged.lookup_bond("Ir", "*", "C", "sp2") is None
    # σ-class hits in both sub-libs and merges as a "both-agree" tuple.
    sigma_hit = merged.lookup_bond("Ir", "*", "N", "sp3")
    assert sigma_hit is not None


def test_singleton_cache_safe_under_env_toggle(v4_b3_synth, monkeypatch):
    """The guard reads the env per-call, so toggling the flag does NOT
    leave the library cache in a stale state -- the next lookup sees the
    new behaviour without a singleton flush."""
    lib = gml.GripLibrary.get(v4_b3_synth)
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")

    # Flag OFF -> v4 hit
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", raising=False)
    assert lib.lookup_bond("Ir", "*", "C", "sp2") is not None

    # Flag ON -> None (no cache flush needed)
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "1")
    assert lib.lookup_bond("Ir", "*", "C", "sp2") is None

    # Toggle OFF again -> v4 hit comes back
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_LIB_DISABLE", "0")
    assert lib.lookup_bond("Ir", "*", "C", "sp2") is not None
