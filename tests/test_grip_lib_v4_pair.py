"""Tests for v4 pair-resolved tables + TM-aware lookup fallback.

These tests verify:

1. The new env flag ``DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK`` is byte-
   identical to legacy when unset (default OFF contract).
2. With the env flag set on a v4 library, ``lookup_bond`` returns the
   correct pair-resolved value for metal-organic bonds that the legacy
   centred-fragment chain cannot satisfy.
3. ``lookup_angle`` and ``lookup_improper`` honour the v4 fallback in the
   same additive manner.
4. A v3 library (no pair tables) returns ``None`` with the flag set --
   importantly, it does NOT inject chemically wrong distances that pool
   over the centre's entire neighbour list.

We synthesise a tiny v4 npz in a tmp dir so the tests run without the full
1.4 M CCDC build.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
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


def _key(parsed) -> str:
    return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)


def _build_v4_synth(tmp: Path) -> Path:
    """Build a tiny v4 npz with one C-C fragment key + pair tables for Ir-C.
    """
    # v3-compatible fragment key (so legacy chain has one hit for C-C sp3)
    cc_key = _key(["C", "sp3", [["C", 1, -1, "sp3"]], []])
    master_keys = [cc_key]
    n_master = len(master_keys)
    bond_mu = np.array([1.540], dtype=np.float64)
    bond_sigma = np.array([0.010], dtype=np.float64)
    bond_n = np.array([1000], dtype=np.int32)
    zero_p = np.full(n_master, np.nan, dtype=np.float32)
    angle_mu = np.full(n_master, np.nan, dtype=np.float64)
    angle_sigma = np.full(n_master, np.nan, dtype=np.float64)
    angle_n = np.zeros(n_master, dtype=np.int32)
    improper_mu = np.full(n_master, np.nan, dtype=np.float64)
    improper_sigma = np.full(n_master, np.nan, dtype=np.float64)
    improper_n = np.zeros(n_master, dtype=np.int32)

    # v4 pair tables: flat sorted-pair key [Z_lo, hyb_lo, Z_hi, hyb_hi]
    # Ir-C(sp2) bond (real ~2.0 A)
    irc_sp2_key = _key(["C", "sp2", "Ir", "*"])
    # Pd-Cl bond (real ~2.3 A)
    pdcl_key = _key(["Cl", "sp3", "Pd", "*"])
    pair_bond_keys = [irc_sp2_key, pdcl_key]
    pair_bond_mu = np.array([2.005, 2.305], dtype=np.float64)
    pair_bond_sigma = np.array([0.030, 0.020], dtype=np.float64)
    pair_bond_n = np.array([150, 250], dtype=np.int32)

    # Angle C-Ir-C centered at Ir (90 deg)
    ang_circ_key = _key(["Ir", "*", "C", "C"])
    triple_angle_keys = [ang_circ_key]
    triple_angle_mu = np.array([89.5], dtype=np.float64)
    triple_angle_sigma = np.array([3.5], dtype=np.float64)
    triple_angle_n = np.array([80], dtype=np.int32)

    # Improper at Ir w/ 3 C neighbours (planar = 0 deg)
    imp_key = _key(["Ir", "*", sorted(["C", "C", "C"])])
    improper_pair_keys = [imp_key]
    improper_pair_mu = np.array([0.5], dtype=np.float64)
    improper_pair_sigma = np.array([4.0], dtype=np.float64)
    improper_pair_n = np.array([60], dtype=np.int32)

    out = tmp / "grip_lib_v4_synth.npz"
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
        orig_idx=np.zeros(n_master, dtype=np.int32),
        bond_mu=bond_mu, bond_sigma=bond_sigma, bond_n=bond_n,
        bond_p5=zero_p, bond_p50=zero_p, bond_p95=zero_p,
        angle_mu=angle_mu, angle_sigma=angle_sigma, angle_n=angle_n,
        angle_p5=zero_p, angle_p50=zero_p, angle_p95=zero_p,
        improper_mu=improper_mu, improper_sigma=improper_sigma, improper_n=improper_n,
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
def v4_synth(tmp_path):
    """Yield a clean v4 synthetic library, with the singleton cache cleared."""
    p = _build_v4_synth(tmp_path)
    gml.GripLibrary._SINGLETONS.clear()
    yield p
    gml.GripLibrary._SINGLETONS.clear()


def test_v4_pair_tables_loaded(v4_synth):
    lib = gml.GripLibrary.get(v4_synth)
    assert lib.version == 4
    assert lib.has_pair_tables is True
    assert lib.n_pair_bond == 2
    assert lib.n_triple_angle == 1
    assert lib.n_improper_pair == 1


def test_tm_fallback_default_off(v4_synth, monkeypatch):
    """With the env flag unset, lookup_bond on a metal-organic pair must
    return ``None`` (byte-identical to legacy chain)."""
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", raising=False)
    lib = gml.GripLibrary.get(v4_synth)
    # Legacy chain has no Ir-C key
    assert lib.lookup_bond("Ir", "*", "C", "sp2") is None
    # Legacy chain DOES have C(sp3)-C(sp3) -> still works
    hit = lib.lookup_bond("C", "sp3", "C", "sp3")
    assert hit is not None
    mu, sigma, n = hit
    assert abs(mu - 1.540) < 1e-6
    assert n == 1000


def test_tm_fallback_on_returns_pair_bond(v4_synth, monkeypatch):
    """With the env flag SET, lookup_bond returns the v4 pair-resolved value."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_synth)
    hit = lib.lookup_bond("Ir", "*", "C", "sp2")
    assert hit is not None
    mu, sigma, n = hit
    assert abs(mu - 2.005) < 1e-6
    assert abs(sigma - 0.030) < 1e-6
    assert n == 150
    # Reverse orientation hits the same pair key (canonical sort)
    hit_rev = lib.lookup_bond("C", "sp2", "Ir", "*")
    assert hit_rev == hit


def test_tm_fallback_hyb_wildcard(v4_synth, monkeypatch):
    """Passing a more specific hyb than stored must fall back through the
    wildcard chain inside ``_v4_lookup_bond``."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_synth)
    # The library stored Ir with hyb "*"; lookup with "sp2" must wildcard to it.
    hit = lib.lookup_bond("Ir", "sp2", "C", "sp2")
    assert hit is not None
    assert abs(hit[0] - 2.005) < 1e-6


def test_tm_fallback_min_n_filters(v4_synth, monkeypatch):
    """When min_n exceeds the stored count, the v4 fallback returns None."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_synth)
    assert lib.lookup_bond("Ir", "*", "C", "sp2", min_n=1000) is None
    assert lib.lookup_bond("Ir", "*", "C", "sp2", min_n=10) is not None


def test_tm_fallback_angle_centered_at_metal(v4_synth, monkeypatch):
    """C-Ir-C angle: centre Ir, neighbours C, C -> hits the v4 triple_angle."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_synth)
    hit = lib.lookup_angle("C", "Ir", "*", "C")
    assert hit is not None
    mu, sigma, n = hit
    assert abs(mu - 89.5) < 1e-3
    assert n == 80


def test_tm_fallback_improper_metal_centre(v4_synth, monkeypatch):
    """3-C improper at Ir (planar) -- v4 improper_pair table hit."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    lib = gml.GripLibrary.get(v4_synth)
    hit = lib.lookup_improper("Ir", "*", ["C", "C", "C"])
    assert hit is not None
    mu, _, n = hit
    assert abs(mu - 0.5) < 1e-3
    assert n == 60


def test_v3_library_returns_none_with_flag(tmp_path, monkeypatch):
    """Critical: a v3 library (no pair tables) must NOT inject wrong values
    when the TM flag is set -- it should simply return None for queries the
    legacy chain cannot satisfy.
    """
    # Build a tiny v3-shaped npz (no pair_bond_keys)
    cc_key = _key(["C", "sp3", [["C", 1, -1, "sp3"]], []])
    n_master = 1
    out = tmp_path / "grip_lib_v3_synth.npz"
    np.savez_compressed(
        out,
        version=np.int32(3),
        n_master=np.int32(n_master),
        n_orig=np.int32(n_master),
        n_torsion=np.int32(0),
        n_orig_torsion=np.int32(0),
        keys=np.array([cc_key], dtype=np.object_),
        orig_keys=np.array([cc_key], dtype=np.object_),
        orig_idx=np.zeros(n_master, dtype=np.int32),
        bond_mu=np.array([1.540], dtype=np.float64),
        bond_sigma=np.array([0.010], dtype=np.float64),
        bond_n=np.array([1000], dtype=np.int32),
        bond_p5=np.full(n_master, np.nan, dtype=np.float32),
        bond_p50=np.full(n_master, np.nan, dtype=np.float32),
        bond_p95=np.full(n_master, np.nan, dtype=np.float32),
        angle_mu=np.full(n_master, np.nan, dtype=np.float64),
        angle_sigma=np.full(n_master, np.nan, dtype=np.float64),
        angle_n=np.zeros(n_master, dtype=np.int32),
        angle_p5=np.full(n_master, np.nan, dtype=np.float32),
        angle_p50=np.full(n_master, np.nan, dtype=np.float32),
        angle_p95=np.full(n_master, np.nan, dtype=np.float32),
        improper_mu=np.full(n_master, np.nan, dtype=np.float64),
        improper_sigma=np.full(n_master, np.nan, dtype=np.float64),
        improper_n=np.zeros(n_master, dtype=np.int32),
        improper_p5=np.full(n_master, np.nan, dtype=np.float32),
        improper_p50=np.full(n_master, np.nan, dtype=np.float32),
        improper_p95=np.full(n_master, np.nan, dtype=np.float32),
        fb_flat=np.zeros(0, dtype=np.int32),
        fb_ptr=np.zeros(n_master + 1, dtype=np.int32),
    )

    monkeypatch.setenv("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", "1")
    gml.GripLibrary._SINGLETONS.clear()
    lib = gml.GripLibrary.get(out)
    assert lib.has_pair_tables is False
    # The Ir-C lookup must NOT return a pooled wrong value -- it returns None.
    assert lib.lookup_bond("Ir", "*", "C", "sp2") is None
    # Legacy chain still works on existing keys.
    hit = lib.lookup_bond("C", "sp3", "C", "sp3")
    assert hit is not None and abs(hit[0] - 1.540) < 1e-6
    gml.GripLibrary._SINGLETONS.clear()
