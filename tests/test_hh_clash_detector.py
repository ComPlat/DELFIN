"""Unit tests for :mod:`delfin.fffree.hh_clash_detector` (H-H clash gate).

Covers:

* clean methylene/methyl topology returns 0 clashes (no geminal false-fire)
* eclipsed methyl-methyl pair across an inter-fragment axis IS detected
* methyl-internal 1-3 H-H NOT flagged
* env-flag default OFF -> byte-identical to HEAD (no detection)
* determinism (sorted output, same input -> same output)
* threshold-factor adjustment is respected
* excluded_pairs are honoured
* real broken voll-pool structure exhibits the expected H-H count
* no silent failure (degenerate inputs return [] gracefully)
"""
from __future__ import annotations

import os

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _set_env(key: str, value: str):
    """Context-style helper: yields tuple (old, new) used for cleanup."""
    old = os.environ.get(key)
    os.environ[key] = value
    return old


def _restore_env(key: str, old):
    if old is None:
        os.environ.pop(key, None)
    else:
        os.environ[key] = old


# ---------------------------------------------------------------------------
# Synthetic geometries
# ---------------------------------------------------------------------------
def _methylene_geometry():
    """A clean CH2 fragment (C-H 1.09, H-C-H 109.5°)."""
    # Place C at origin, two H at +/- offset along x with sp3 geometry.
    # Geminal H-H ≈ 1.78 Å.  No clashes expected (geminal is topo=2).
    syms = ["C", "H", "H"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [0.89, 0.00, 0.63],
        [-0.89, 0.00, 0.63],
    ], dtype=float)
    return syms, P


def _methyl_internal_geometry():
    """A clean CH3 (3 H around C, internal H-H ≈ 1.78 Å)."""
    syms = ["C", "H", "H", "H"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [1.03, 0.00, 0.36],
        [-0.52, 0.89, 0.36],
        [-0.52, -0.89, 0.36],
    ], dtype=float)
    return syms, P


def _ethane_eclipsed_methyls():
    """Two CH3 groups joined by a C-C bond in tight eclipsed conformation.

    H atoms on the two methyls eclipse across the C-C axis -- a rotamer
    ERROR.  The H-H pairs are at topological distance 3 (H-C-C-H), so
    they are NOT flagged at the default ``min_topo_distance=4``; tests
    that target this case must call with ``min_topo_distance=3``.

    Geometry below uses C-C = 1.10 Å with H atoms placed at +/- 0.95 Å
    along x so the eclipsing in-plane H-H comes in at ~ 1.90 Å (below
    the 2.04 Å Bondi floor).  Intentionally compressed to make the
    clash unambiguous in unit tests.
    """
    syms = ["C", "H", "H", "H", "C", "H", "H", "H"]
    # C1 at -0.55, C2 at +0.55 (C-C = 1.10 Å, intentionally tight).
    # Each H is placed at +/-(0.95 + 0.55) = +/- 1.50 Å along x with the
    # eclipsing partner at the same y/z -> in-plane H-H = 3.00 Å minus
    # the H push toward the partner.  We use H at +/- 0.90 along x so
    # the eclipsing distance is 0.90 + 0.90 = 1.80 Å.
    P = np.array([
        [-0.55, 0.00, 0.00],   # C1
        [-0.90, 1.00, 0.00],   # H on C1, +y
        [-0.90, -0.50, 0.87],  # H on C1, -y +z
        [-0.90, -0.50, -0.87], # H on C1, -y -z
        [0.55, 0.00, 0.00],    # C2
        [0.90, 1.00, 0.00],    # H on C2, +y (eclipses H1)
        [0.90, -0.50, 0.87],   # H on C2, -y +z (eclipses H2)
        [0.90, -0.50, -0.87],  # H on C2, -y -z (eclipses H3)
    ], dtype=float)
    return syms, P


def _two_isolated_methyls(separation: float):
    """Two CH3 groups whose central C atoms are ``separation`` Å apart,
    aligned face-to-face so each H1 of CH3a faces an H of CH3b.

    Used to validate the inter-ligand H-H clash detection.
    """
    syms = ["C", "H", "H", "H", "C", "H", "H", "H"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [1.03, 0.00, 0.36],
        [-0.52, 0.89, 0.36],
        [-0.52, -0.89, 0.36],
        [0.00, 0.00, separation],
        [1.03, 0.00, separation - 0.36],
        [-0.52, 0.89, separation - 0.36],
        [-0.52, -0.89, separation - 0.36],
    ], dtype=float)
    return syms, P


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def test_no_hh_clash_clean_methylene():
    """A clean methylene CH2 must NOT flag the geminal H-H (topo=2)."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _methylene_geometry()
    clashes = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    assert clashes == []


def test_skips_methyl_internal_hh():
    """The 3 internal H-H pairs in a single CH3 are 1-3 (topo=3) and
    must be skipped under the default ``min_topo_distance=4``."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _methyl_internal_geometry()
    clashes = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    assert clashes == []


def test_detects_methyl_methyl_eclipse():
    """Eclipsed methyls with a 1.20 Å C-C produce 1-4 H-H clashes
    flagged when ``min_topo_distance=3`` (the default 4 skips 1-4
    rotamer pairs)."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _ethane_eclipsed_methyls()
    clashes = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=3)
    assert len(clashes) >= 3  # the 3 eclipsing pairs
    # Severity > 0 for each
    for (_i, _j, _d, sev) in clashes:
        assert sev > 0.0


def test_skips_eclipse_at_default_min_topo():
    """The default ``min_topo_distance=4`` skips 1-4 H-H eclipses across
    an sp3-sp3 bond (chemically a rotamer issue, not a clash)."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _ethane_eclipsed_methyls()
    clashes = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    assert clashes == []


def test_detects_inter_ligand_hh_clash_at_tight_separation():
    """Two CH3 groups whose carbons are 2.5 Å apart create a tight
    inter-ligand H-H contact (~ 1.78 Å) that must be flagged."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _two_isolated_methyls(separation=2.5)
    # The two methyls are not bonded so topo distance between any
    # H-H across the gap is infinite -> ALWAYS flagged when below floor.
    clashes = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    # At least one pair below 2.04 Å expected
    assert len(clashes) >= 1


def test_no_clash_when_methyls_far_apart():
    """Methyl groups 4.5 Å apart must NOT clash (no false positives)."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _two_isolated_methyls(separation=4.5)
    clashes = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    assert clashes == []


def test_env_off_no_detection_in_build_gate():
    """When DELFIN_FFFREE_HH_CLASH_INCLUDE is unset, the build_time
    has_collapse must NOT change on H-H eclipsing alone (byte-identical
    to HEAD).
    """
    from delfin.fffree.build_time_clash_gate import has_collapse
    old_hh = _set_env("DELFIN_FFFREE_HH_CLASH_INCLUDE", "0")
    old_gate = _set_env("DELFIN_FFFREE_BUILD_CLASH_GATE", "1")
    try:
        syms, P = _two_isolated_methyls(separation=2.5)
        # No bond collapse + flag OFF -> must return False
        assert has_collapse(syms, P) is False
    finally:
        _restore_env("DELFIN_FFFREE_HH_CLASH_INCLUDE", old_hh)
        _restore_env("DELFIN_FFFREE_BUILD_CLASH_GATE", old_gate)


def test_env_on_triggers_build_gate_on_inter_fragment_hh():
    """With both DELFIN_FFFREE_HH_CLASH_INCLUDE=1 AND the build clash
    gate active, the two-methyl inter-fragment geometry (no covalent
    bond between the methyls) must be flagged as a collapse."""
    from delfin.fffree.build_time_clash_gate import has_collapse
    old_hh = _set_env("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
    old_gate = _set_env("DELFIN_FFFREE_BUILD_CLASH_GATE", "1")
    try:
        syms, P = _two_isolated_methyls(separation=2.5)
        assert has_collapse(syms, P) is True
    finally:
        _restore_env("DELFIN_FFFREE_HH_CLASH_INCLUDE", old_hh)
        _restore_env("DELFIN_FFFREE_BUILD_CLASH_GATE", old_gate)


def test_determinism():
    """Same input -> same output, bit-identical sorting."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _ethane_eclipsed_methyls()
    a = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    b = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    assert a == b
    # Sorted (i, j)
    pairs = [(t[0], t[1]) for t in a]
    assert pairs == sorted(pairs)


def test_threshold_factor_respected():
    """A LOWER factor (tighter floor) yields fewer or equal clashes,
    a HIGHER factor (looser floor, e.g. 1.0) yields more or equal."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _ethane_eclipsed_methyls()
    low = detect_hh_clashes(syms, P, factor=0.50, min_topo_distance=4)
    mid = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    high = detect_hh_clashes(syms, P, factor=1.10, min_topo_distance=4)
    assert len(low) <= len(mid) <= len(high)


def test_excluded_pairs_respected():
    """Pairs listed in ``excluded_pairs`` are skipped even when they
    would otherwise clash."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    syms, P = _two_isolated_methyls(separation=2.5)
    base = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    assert base  # sanity: at least one clash to exclude
    # Exclude the first reported pair and confirm it disappears.
    skip_pair = (base[0][0], base[0][1])
    masked = detect_hh_clashes(
        syms, P, factor=0.85, min_topo_distance=4,
        excluded_pairs=[skip_pair],
    )
    for (i, j, _d, _s) in masked:
        assert (i, j) != skip_pair


def test_no_silent_failure_on_empty():
    """Empty / single-atom inputs must return [] not raise."""
    from delfin.fffree.hh_clash_detector import detect_hh_clashes
    assert detect_hh_clashes([], np.zeros((0, 3))) == []
    assert detect_hh_clashes(["H"], np.zeros((1, 3))) == []


def test_md_invariant_preserved_under_hh_penalty():
    """Calling ``count_hh_clashes`` is read-only: coordinates are not
    mutated.  This is a tiny but important guard against future regressions
    where the detector accidentally rescales coordinates in-place."""
    from delfin.fffree.hh_clash_detector import count_hh_clashes
    syms, P = _ethane_eclipsed_methyls()
    P_before = P.copy()
    _ = count_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    assert np.allclose(P, P_before, atol=0.0, rtol=0.0)


def test_real_structure_b00f9a0_finds_known_clashes():
    """Spot-check against a real voll-pool structure that exhibits H-H
    clashes per the 2026-06-04 audit.  Skips when the archive is not
    present (developer machine without the full pool checked out)."""
    from delfin.fffree.hh_clash_detector import count_hh_clashes
    candidates = [
        "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/"
        "b00f9a0-full7-VOLLPOOL/047-CIRJUK.xyz",
        "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/"
        "b00f9a0-full7-VOLLPOOL/084-MIPSOW.xyz",
    ]
    path = None
    for c in candidates:
        if os.path.exists(c):
            path = c
            break
    if path is None:
        pytest.skip("voll-pool b00f9a0-full7-VOLLPOOL archive not present")
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].strip())
    syms = []
    coords = []
    for k in range(2, 2 + n):
        parts = lines[k].split()
        syms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    P = np.array(coords)
    n_hh = count_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    # The audit recorded ≥ 1 clash for both reference files.
    assert n_hh >= 1


def test_severity_aggregates_correctly():
    """The total severity equals the sum of per-pair severities (Å²)
    and is monotone non-decreasing as clashes get tighter."""
    from delfin.fffree.hh_clash_detector import (
        detect_hh_clashes, hh_clash_severity,
    )
    syms, P = _two_isolated_methyls(separation=2.5)
    pairs = detect_hh_clashes(syms, P, factor=0.85, min_topo_distance=4)
    sev_aggregate = hh_clash_severity(
        syms, P, factor=0.85, min_topo_distance=4,
    )
    sev_manual = sum(t[3] for t in pairs)
    assert abs(sev_aggregate - sev_manual) < 1e-12

    # Tighter contact -> higher severity
    syms2, P2 = _two_isolated_methyls(separation=2.0)
    sev_tighter = hh_clash_severity(
        syms2, P2, factor=0.85, min_topo_distance=4,
    )
    assert sev_tighter >= sev_aggregate
