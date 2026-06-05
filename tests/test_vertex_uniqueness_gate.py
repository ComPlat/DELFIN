"""Unit tests for :mod:`delfin.fffree.vertex_uniqueness_gate`.

Covers the gate's default-OFF byte-identity contract, the catastrophic-collapse
detection / drop semantics, and the extended ``contains_metal`` regex that
fixes the charged-hydride misclassification responsible for the bulk of the
4.0 % voll-pool build-collapse population.
"""
from __future__ import annotations

import os

import pytest

from delfin.fffree.vertex_uniqueness_gate import (
    DEFAULT_HARD_THRESHOLD,
    contains_metal_extended,
    drop_collapsed_isomers,
    has_collapsed_pair,
    is_active,
    min_heavy_heavy_distance,
)
from delfin.smiles_converter import _METALS, contains_metal


# ---------------------------------------------------------------------------
# Helper: env scrub fixture so individual tests do not leak flags
# ---------------------------------------------------------------------------
@pytest.fixture
def clean_env(monkeypatch):
    """Strip DELFIN_FFFREE_* env vars so each test starts from the OFF path."""
    for k in list(os.environ.keys()):
        if k.startswith("DELFIN_FFFREE"):
            monkeypatch.delenv(k, raising=False)
    return monkeypatch


# ---------------------------------------------------------------------------
# Gate-activation contract
# ---------------------------------------------------------------------------
def test_gate_default_off(clean_env):
    """No env flag set -> gate inactive."""
    assert is_active() is False


def test_gate_primary_flag_activates(clean_env):
    clean_env.setenv("DELFIN_FFFREE_VERTEX_UNIQUENESS", "1")
    assert is_active() is True


def test_gate_supergate_activates(clean_env):
    """PURE_TRACK3 must also turn the gate on (universal Track-3 contract)."""
    clean_env.setenv("DELFIN_FFFREE_PURE_TRACK3", "1")
    assert is_active() is True


def test_gate_zero_value_inactive(clean_env):
    """Explicit '0' must NOT activate (typical default-off convention)."""
    clean_env.setenv("DELFIN_FFFREE_VERTEX_UNIQUENESS", "0")
    assert is_active() is False


# ---------------------------------------------------------------------------
# Collapse-detection helpers
# ---------------------------------------------------------------------------
def _xyz(coords):
    """Build a minimal XYZ string from ``[(sym, x, y, z), ...]``."""
    lines = [str(len(coords)), "test"]
    for sym, x, y, z in coords:
        lines.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines)


def test_min_heavy_heavy_distance_basic():
    xyz = _xyz([("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0), ("H", 0.0, 0.0, 1.0)])
    assert min_heavy_heavy_distance(xyz) == pytest.approx(1.5, abs=1e-9)


def test_min_heavy_heavy_distance_empty():
    assert min_heavy_heavy_distance("") == float("inf")
    assert min_heavy_heavy_distance(_xyz([("C", 0.0, 0.0, 0.0)])) == float("inf")


def test_min_heavy_heavy_distance_ignores_hydrogen():
    """H-H pair at 0.1 A must not register: the gate is HEAVY-only."""
    xyz = _xyz([
        ("C", 0.0, 0.0, 0.0),
        ("O", 1.4, 0.0, 0.0),
        ("H", 0.5, 0.0, 0.0),
        ("H", 0.5, 0.1, 0.0),
    ])
    assert min_heavy_heavy_distance(xyz) == pytest.approx(1.4, abs=1e-9)


def test_min_heavy_heavy_distance_tolerates_garbage_lines():
    xyz = "\n".join([
        "3", "header",
        "C 0 0 0",
        "garbage line here",
        "C 1.5 0 0",
        "",
        "C 3.0 0 0",
    ])
    assert min_heavy_heavy_distance(xyz) == pytest.approx(1.5, abs=1e-9)


def test_has_collapsed_pair_detects_below_threshold():
    xyz = _xyz([
        ("Fe", 0.0, 0.0, 0.0),
        ("P", 0.012, 0.0, 0.0),  # AQAVOF reference
        ("C", 1.8, 0.0, 0.0),
    ])
    assert has_collapsed_pair(xyz, threshold=0.5) is True


def test_has_collapsed_pair_clean_structure():
    xyz = _xyz([
        ("Fe", 0.0, 0.0, 0.0),
        ("P", 2.2, 0.0, 0.0),
        ("N", -2.0, 0.0, 0.0),
    ])
    assert has_collapsed_pair(xyz, threshold=0.5) is False


def test_has_collapsed_pair_threshold_edge():
    """At exactly the threshold the pair is NOT considered collapsed."""
    xyz = _xyz([
        ("C", 0.0, 0.0, 0.0),
        ("O", 0.5, 0.0, 0.0),
    ])
    assert has_collapsed_pair(xyz, threshold=0.5) is False


# ---------------------------------------------------------------------------
# drop_collapsed_isomers
# ---------------------------------------------------------------------------
COLLAPSED_XYZ = _xyz([
    ("Fe", 0.0, 0.0, 0.0),
    ("P", 0.001, 0.0, 0.0),  # 1 mA collapse
    ("C", 1.83, 0.0, 0.0),
])

CLEAN_XYZ = _xyz([
    ("Fe", 0.0, 0.0, 0.0),
    ("P", 2.22, 0.0, 0.0),
    ("C", -1.83, 0.0, 0.0),
])


def test_drop_off_path_is_identity(clean_env):
    """OFF (default) -> drop_collapsed_isomers returns the input unchanged."""
    isos = [(COLLAPSED_XYZ, "bad"), (CLEAN_XYZ, "good")]
    out = drop_collapsed_isomers(isos)
    assert out is isos  # literal identity (object pass-through)


def test_drop_on_path_removes_collapsed(clean_env):
    clean_env.setenv("DELFIN_FFFREE_VERTEX_UNIQUENESS", "1")
    isos = [(COLLAPSED_XYZ, "bad"), (CLEAN_XYZ, "good")]
    out = drop_collapsed_isomers(isos)
    assert len(out) == 1
    assert out[0][1] == "good"


def test_drop_on_path_keeps_all_clean(clean_env):
    clean_env.setenv("DELFIN_FFFREE_VERTEX_UNIQUENESS", "1")
    isos = [(CLEAN_XYZ, "a"), (CLEAN_XYZ, "b"), (CLEAN_XYZ, "c")]
    out = drop_collapsed_isomers(isos)
    assert len(out) == 3


def test_drop_on_path_empty_input(clean_env):
    clean_env.setenv("DELFIN_FFFREE_VERTEX_UNIQUENESS", "1")
    assert drop_collapsed_isomers([]) == []
    assert drop_collapsed_isomers(None) is None


def test_drop_on_path_tolerates_non_string_xyz(clean_env):
    """Defence-in-depth: malformed list items survive without raising."""
    clean_env.setenv("DELFIN_FFFREE_VERTEX_UNIQUENESS", "1")
    isos = [(None, "weird"), (CLEAN_XYZ, "clean")]
    out = drop_collapsed_isomers(isos)
    assert len(out) == 2


# ---------------------------------------------------------------------------
# contains_metal extended regex
# ---------------------------------------------------------------------------
# Five SMILES that the legacy regex misses but real voll-pool data contains.
_CHARGED_HYDRIDE_SMILES = {
    "AQAVOF":
        "CC(C)[P+]1(C(C)C)N(C)C2=CC=CC3=[N+]2[FeH2-4]1([C]#[O+])"
        "[P+](C(C)C)(C(C)C)N3C",
    "EGUKUN":
        "C[P+](C)(C1=CC=CC=C1)[IrH-3]([Cl])([Cl])([C]#[O+])"
        "[P+](C)(C)C1=CC=CC=C1",
    "HAZTIN":
        "CC(C)[P+](C(C)C)(C(C)C)[IrH-3]([NH3+])([Cl])"
        "([C]1=C(F)C(F)=NC(F)=C1F)[P+](C(C)C)(C(C)C)C(C)C",
    "RUQVAC":
        "COC1=CC=C([P+](C2=CC=C(OC)C=C2)(C2=CC=C(OC)C=C2)"
        "[RhH-3]([Cl])([Cl])([C]#[N+]C23CC4CC(CC(C4)C2)C3)"
        "[P+](C2=CC=C(OC)C=C2)(C2=CC=C(OC)C=C2)C2=CC=C(OC)C=C2)C=C1",
}


@pytest.mark.parametrize("label", sorted(_CHARGED_HYDRIDE_SMILES))
def test_contains_metal_extended_catches_charged_hydride(label):
    """Extended regex must recognise [FeH2-N] / [IrH-N] / [RhH-N] etc."""
    smi = _CHARGED_HYDRIDE_SMILES[label]
    assert contains_metal_extended(smi, _METALS) is True


def test_contains_metal_extended_negative_on_pure_organic():
    smi = "c1ccc(cc1)C(=O)O"  # benzoic acid
    assert contains_metal_extended(smi, _METALS) is False


def test_contains_metal_extended_accepts_bare_brackets():
    """``[Fe]`` (rare but legal) must still register."""
    assert contains_metal_extended("CC[Fe]CC", _METALS) is True


def test_contains_metal_extended_accepts_charge_only():
    """``[Fe+3]`` (no H) must still register (back-compat with legacy regex)."""
    assert contains_metal_extended("CCC.[Fe+3]", _METALS) is True


# ---------------------------------------------------------------------------
# Integration: contains_metal() hook (legacy regex + extended on gate-active)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("label", sorted(_CHARGED_HYDRIDE_SMILES))
def test_contains_metal_default_off_misses_hydride(clean_env, label):
    """Default OFF -> the contains_metal() public API still returns False
    for these SMILES (we preserve the legacy mis-detection to keep HEAD
    byte-identical when the gate is off)."""
    smi = _CHARGED_HYDRIDE_SMILES[label]
    assert contains_metal(smi) is False


@pytest.mark.parametrize("label", sorted(_CHARGED_HYDRIDE_SMILES))
def test_contains_metal_gate_active_catches_hydride(clean_env, label):
    """Gate ON -> contains_metal() now picks them up (routes them to the
    metal-aware build path instead of the organic-conformer pool)."""
    clean_env.setenv("DELFIN_FFFREE_VERTEX_UNIQUENESS", "1")
    smi = _CHARGED_HYDRIDE_SMILES[label]
    assert contains_metal(smi) is True


# ---------------------------------------------------------------------------
# Determinism: same input + flag -> identical output
# ---------------------------------------------------------------------------
def test_drop_deterministic_repeat(clean_env):
    """Two repeat calls produce structurally identical results (no RNG)."""
    clean_env.setenv("DELFIN_FFFREE_VERTEX_UNIQUENESS", "1")
    isos = [
        (COLLAPSED_XYZ, "bad-1"),
        (CLEAN_XYZ, "good-1"),
        (COLLAPSED_XYZ, "bad-2"),
        (CLEAN_XYZ, "good-2"),
    ]
    a = drop_collapsed_isomers(list(isos))
    b = drop_collapsed_isomers(list(isos))
    assert a == b
    assert [x[1] for x in a] == ["good-1", "good-2"]


def test_default_hard_threshold_constant():
    """Public constant exposed for downstream telemetry and tests."""
    assert DEFAULT_HARD_THRESHOLD == pytest.approx(0.5, abs=1e-9)
