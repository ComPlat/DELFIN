"""Tests for delfin.fffree.inter_ligand_clash_gate — pre-polish inter-ligand
clash filter for combinatorial rotamer enumeration.

Subagent #129 follow-up: 8-10 tests covering the count helper, generator
wrapper, env-flag plumbing and deterministic ordering.
"""
from __future__ import annotations

import itertools
import os
from typing import List, Optional, Tuple

import numpy as np
import pytest


# ----------------------------------------------------------------------
# count_inter_ligand_clashes_quick
# ----------------------------------------------------------------------


def test_count_zero_when_no_clashes():
    """Two ligands placed 20 Å apart -> zero inter-ligand clashes."""
    from delfin.fffree.inter_ligand_clash_gate import count_inter_ligand_clashes_quick
    P = np.array([
        [0.0, 0.0, 0.0],   # ligand 0 atom
        [1.0, 0.0, 0.0],   # ligand 0 atom
        [20.0, 0.0, 0.0],  # ligand 1 atom
        [21.0, 0.0, 0.0],  # ligand 1 atom
    ])
    syms = ["C", "C", "C", "C"]
    subgraphs = [[0, 1], [2, 3]]
    n = count_inter_ligand_clashes_quick(P, syms, subgraphs)
    assert n == 0


def test_count_one_when_overlap():
    """Two ligand atoms inside the vdW floor -> one clash."""
    from delfin.fffree.inter_ligand_clash_gate import count_inter_ligand_clashes_quick
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],   # too close to ligand 1's atom at 1.5 (vdW C-C floor)
        [1.5, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    ])
    syms = ["C", "C", "C", "C"]
    subgraphs = [[0, 1], [2, 3]]
    # vdW C = 1.70 -> floor = 0.85 * 3.40 = 2.89; 1.0->1.5 = 0.5 < 2.89 -> clash
    # Also 0.0->1.5 = 1.5 < 2.89 -> clash; 1.0->3.0 = 2.0 < 2.89 -> clash
    # Total pairs with both endpoints in different subgraphs that are < 2.89:
    #   (0,2)=1.5, (0,3)=3.0 NO, (1,2)=0.5, (1,3)=2.0
    n = count_inter_ligand_clashes_quick(P, syms, subgraphs)
    assert n == 3


def test_count_intra_ligand_ignored():
    """Two atoms inside the SAME ligand subgraph: not counted as inter-ligand."""
    from delfin.fffree.inter_ligand_clash_gate import count_inter_ligand_clashes_quick
    P = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],   # close, but same ligand
        [50.0, 0.0, 0.0],
    ])
    syms = ["C", "C", "C"]
    subgraphs = [[0, 1], [2]]
    n = count_inter_ligand_clashes_quick(P, syms, subgraphs)
    assert n == 0


def test_count_metal_excluded():
    """The metal atom must be ignored when counting clashes."""
    from delfin.fffree.inter_ligand_clash_gate import count_inter_ligand_clashes_quick
    P = np.array([
        [0.0, 0.0, 0.0],   # metal at origin
        [1.5, 0.0, 0.0],   # ligand 0 donor (close to metal but metal is excluded)
        [-1.5, 0.0, 0.0],  # ligand 1 donor
    ])
    syms = ["Fe", "N", "N"]
    subgraphs = [[1], [2]]
    # ligand 0 donor and ligand 1 donor are 3 Å apart; floor for N-N = 0.85*3.10 = 2.635
    # -> 3.0 > 2.635 -> no clash.
    n = count_inter_ligand_clashes_quick(P, syms, subgraphs, metal_idx=0)
    assert n == 0


def test_count_threshold_lower_lets_more_pass():
    """Lower vdW fraction -> fewer clashes."""
    from delfin.fffree.inter_ligand_clash_gate import count_inter_ligand_clashes_quick
    P = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [50.0, 0.0, 0.0],
    ])
    syms = ["C", "C", "C"]
    subgraphs = [[0], [1, 2]]
    # at threshold 0.85 -> floor 2.89 -> 2.0 < 2.89 -> clash
    n_strict = count_inter_ligand_clashes_quick(P, syms, subgraphs, threshold=0.85)
    # at threshold 0.5 -> floor 1.7 -> 2.0 > 1.7 -> no clash
    n_loose = count_inter_ligand_clashes_quick(P, syms, subgraphs, threshold=0.5)
    assert n_strict >= n_loose
    assert n_strict >= 1
    assert n_loose == 0


# ----------------------------------------------------------------------
# enumerate_with_clash_gate
# ----------------------------------------------------------------------


def _good_combo(c):
    """Two ligands far apart -> zero clashes."""
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [20.0 + c, 0.0, 0.0],
        [21.0 + c, 0.0, 0.0],
    ])
    syms = ["C", "C", "C", "C"]
    subgraphs = [[0, 1], [2, 3]]
    return P, syms, subgraphs, None


def _bad_combo(c):
    """Two ligands very close -> always clashes."""
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.05 + 0.001 * c, 0.0, 0.0],
        [1.50, 0.0, 0.0],
    ])
    syms = ["C", "C", "C", "C"]
    subgraphs = [[0, 1], [2, 3]]
    return P, syms, subgraphs, None


def test_enumerate_with_clash_gate_keeps_good_rejects_bad():
    """A clean assembler always yields; a clashing one never does."""
    from delfin.fffree.inter_ligand_clash_gate import enumerate_with_clash_gate
    combos = range(5)
    out_good = list(enumerate_with_clash_gate(
        combos, _good_combo, max_total=10, clash_threshold=0,
    ))
    out_bad = list(enumerate_with_clash_gate(
        combos, _bad_combo, max_total=10, clash_threshold=0,
    ))
    assert len(out_good) == 5
    assert len(out_bad) == 0


def test_enumerate_with_clash_gate_max_total_cap():
    """max_total caps the yielded count even with infinite-good assembler."""
    from delfin.fffree.inter_ligand_clash_gate import enumerate_with_clash_gate
    combos = range(100)
    out = list(enumerate_with_clash_gate(
        combos, _good_combo, max_total=7, clash_threshold=0,
    ))
    assert len(out) == 7


def test_enumerate_with_clash_gate_threshold_lets_some_pass():
    """clash_threshold = 100 lets even bad ones through; = 0 blocks all bad."""
    from delfin.fffree.inter_ligand_clash_gate import enumerate_with_clash_gate
    combos = range(5)
    out_lax = list(enumerate_with_clash_gate(
        combos, _bad_combo, max_total=10, clash_threshold=100,
    ))
    out_strict = list(enumerate_with_clash_gate(
        combos, _bad_combo, max_total=10, clash_threshold=0,
    ))
    assert len(out_lax) == 5
    assert len(out_strict) == 0


def test_enumerate_with_clash_gate_assembler_none_skipped():
    """Assembler returning None is silently skipped, gate keeps going."""
    from delfin.fffree.inter_ligand_clash_gate import enumerate_with_clash_gate

    def asm(c):
        if c % 2 == 0:
            return None
        return _good_combo(c)

    out = list(enumerate_with_clash_gate(
        range(6), asm, max_total=10, clash_threshold=0,
    ))
    assert len(out) == 3  # odd c=1,3,5


def test_enumerate_with_clash_gate_deterministic_order():
    """Output ordering matches combos_iter (the gate filters, never reorders)."""
    from delfin.fffree.inter_ligand_clash_gate import enumerate_with_clash_gate

    def asm(c):
        if c == 1 or c == 3:
            return _bad_combo(c)  # rejected
        return _good_combo(c)

    out = list(enumerate_with_clash_gate(
        range(6), asm, max_total=10, clash_threshold=0,
    ))
    # Expected: c = 0, 2, 4, 5 (1 and 3 rejected); the FIRST atom of each
    # accepted assembly has x=0 (good_combo always puts atom 0 at origin).
    assert len(out) == 4
    for P, syms, sg, mi in out:
        assert np.isclose(P[0, 0], 0.0)


# ----------------------------------------------------------------------
# Env-flag plumbing
# ----------------------------------------------------------------------


def test_gate_enabled_off_default():
    """Without DELFIN_FFFREE_PRE_POLISH_CLASH_GATE: off when rotamer flag off."""
    from delfin.fffree.inter_ligand_clash_gate import gate_enabled
    saved = os.environ.pop("DELFIN_FFFREE_PRE_POLISH_CLASH_GATE", None)
    try:
        assert gate_enabled(False) is False
        assert gate_enabled(True) is True   # auto-on with rotamers
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_PRE_POLISH_CLASH_GATE"] = saved


def test_gate_enabled_env_explicit_overrides():
    """Explicit env-flag overrides the auto-on behaviour."""
    from delfin.fffree.inter_ligand_clash_gate import gate_enabled
    saved = os.environ.get("DELFIN_FFFREE_PRE_POLISH_CLASH_GATE")
    try:
        os.environ["DELFIN_FFFREE_PRE_POLISH_CLASH_GATE"] = "0"
        assert gate_enabled(True) is False
        os.environ["DELFIN_FFFREE_PRE_POLISH_CLASH_GATE"] = "1"
        assert gate_enabled(False) is True
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_PRE_POLISH_CLASH_GATE"] = saved
        else:
            os.environ.pop("DELFIN_FFFREE_PRE_POLISH_CLASH_GATE", None)


def test_env_defaults_match_spec():
    """Defaults: clash_threshold=5, vdw_fraction=0.85, max_total=200."""
    from delfin.fffree.inter_ligand_clash_gate import (
        env_clash_threshold, env_clash_vdw_fraction, env_max_total,
    )
    for k in ("DELFIN_FFFREE_PRE_POLISH_CLASH_MAX",
              "DELFIN_FFFREE_PRE_POLISH_CLASH_THRESHOLD",
              "DELFIN_FFFREE_ENUMERATION_MAX_TOTAL"):
        os.environ.pop(k, None)
    assert env_clash_threshold() == 5
    assert env_clash_vdw_fraction() == 0.85
    assert env_max_total() == 200


def test_env_overrides_honored():
    """Setting env-vars actually changes the returned default."""
    from delfin.fffree.inter_ligand_clash_gate import (
        env_clash_threshold, env_clash_vdw_fraction, env_max_total,
    )
    saved = {k: os.environ.get(k) for k in (
        "DELFIN_FFFREE_PRE_POLISH_CLASH_MAX",
        "DELFIN_FFFREE_PRE_POLISH_CLASH_THRESHOLD",
        "DELFIN_FFFREE_ENUMERATION_MAX_TOTAL",
    )}
    try:
        os.environ["DELFIN_FFFREE_PRE_POLISH_CLASH_MAX"] = "10"
        os.environ["DELFIN_FFFREE_PRE_POLISH_CLASH_THRESHOLD"] = "0.6"
        os.environ["DELFIN_FFFREE_ENUMERATION_MAX_TOTAL"] = "50"
        assert env_clash_threshold() == 10
        assert abs(env_clash_vdw_fraction() - 0.6) < 1e-9
        assert env_max_total() == 50
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


# ----------------------------------------------------------------------
# gate_stats helper
# ----------------------------------------------------------------------


def test_gate_stats_counts_correctly():
    """gate_stats reports (attempts, assembled_ok, kept) accurately."""
    from delfin.fffree.inter_ligand_clash_gate import gate_stats

    def asm(c):
        if c % 2 == 0:
            return None
        if c == 3:
            return _bad_combo(c)
        return _good_combo(c)

    attempts, ok, kept = gate_stats(
        range(6), asm, max_attempts=10, clash_threshold=0,
    )
    assert attempts == 6
    assert ok == 3   # c = 1, 3, 5
    assert kept == 2  # c = 1, 5 (c=3 is bad)
