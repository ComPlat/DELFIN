"""Unit tests for the ASYMMETRIC permutation-invariant duplicate pruner
(``delfin.fffree.permute_dedup``).

Governing invariant (user): "nur die SCHLECHTEN aussortieren, nicht die guten".
The pruner must merge ONLY genuinely-identical frames and KEEP geometrically
distinct conformers, even when a large symmetry group could relabel one into a
sub-threshold alignment of the other (high symmetry -> MORE conservative).

Synthetic, fast, deterministic.  No heavy SMILES builds.
"""
import os

import pytest

from delfin.fffree import permute_dedup as pd


def _xyz(atoms):
    return "\n".join(f"{s} {x:.4f} {y:.4f} {z:.4f}" for (s, x, y, z) in atoms)


@pytest.fixture(autouse=True)
def _clear_env():
    saved = {k: v for k, v in os.environ.items()
             if k.startswith("DELFIN_FFFREE_PERMUTE_DEDUP")}
    for k in list(os.environ):
        if k.startswith("DELFIN_FFFREE_PERMUTE_DEDUP"):
            del os.environ[k]
    yield
    for k in list(os.environ):
        if k.startswith("DELFIN_FFFREE_PERMUTE_DEDUP"):
            del os.environ[k]
    os.environ.update(saved)


# A square-planar M with 4 identical CN donors -> high graph symmetry (the regime
# the old fixed 0.5-A threshold over-merged).
def _square_planar_cn():
    a = [("Ni", 0.0, 0.0, 0.0)]
    for (x, y) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        a.append(("C", 1.9 * x, 1.9 * y, 0.0))
        a.append(("N", 3.05 * x, 3.05 * y, 0.0))
    return a


def test_off_is_identity():
    iso = [(_xyz(_square_planar_cn()), "a"), (_xyz(_square_planar_cn()), "b")]
    os.environ.pop("DELFIN_FFFREE_PERMUTE_DEDUP", None)
    assert pd.dedup_ensemble(iso) is iso
    os.environ["DELFIN_FFFREE_PERMUTE_DEDUP"] = "0"
    assert pd.dedup_ensemble(iso) is iso


def test_true_duplicate_merges_to_one():
    # Two byte-identical frames are a TRUE duplicate -> exactly one survives.
    os.environ["DELFIN_FFFREE_PERMUTE_DEDUP"] = "1"
    frame = _xyz(_square_planar_cn())
    iso = [(frame, "a"), (frame, "b")]
    out = pd.dedup_ensemble(iso)
    assert len(out) == 1
    assert out[0][1] == "a"  # first member of the cluster kept


def test_relabeled_true_duplicate_merges():
    # Same geometry up to swapping two symmetry-equivalent CN arms (a graph
    # automorphism that just permutes labels): still a TRUE duplicate -> merge.
    os.environ["DELFIN_FFFREE_PERMUTE_DEDUP"] = "1"
    base = _square_planar_cn()
    swapped = list(base)
    # swap the +x arm (C idx1,N idx2) with the -x arm (C idx3,N idx4): identical
    # geometry, different atom labels.
    swapped[1], swapped[3] = base[3], base[1]
    swapped[2], swapped[4] = base[4], base[2]
    iso = [(_xyz(base), "a"), (_xyz(swapped), "b")]
    out = pd.dedup_ensemble(iso)
    assert len(out) == 1  # relabel-equivalent identical geometry -> merged


def test_distinct_conformer_kept_high_symmetry():
    # Two GEOMETRICALLY distinct conformers of a high-symmetry complex must be
    # KEPT.  Distort one arm clearly (bend +y arm by ~35 deg out of plane) ->
    # no symmetry relabelling can align them under the tightened threshold.
    os.environ["DELFIN_FFFREE_PERMUTE_DEDUP"] = "1"
    base = _square_planar_cn()
    distinct = list(base)
    import math
    # bend the +y donor arm out of plane
    distinct[5] = ("C", 0.0, 1.9 * math.cos(0.6), 1.9 * math.sin(0.6))
    distinct[6] = ("N", 0.0, 3.05 * math.cos(0.6), 3.05 * math.sin(0.6))
    iso = [(_xyz(base), "a"), (_xyz(distinct), "b")]
    out = pd.dedup_ensemble(iso)
    assert len(out) == 2  # distinct conformers both kept


def test_small_geometric_difference_kept():
    # A modest but real geometric difference (one arm rotated ~0.4 A) is a
    # distinct conformer, not a duplicate, under the tightened 0.25-A threshold.
    os.environ["DELFIN_FFFREE_PERMUTE_DEDUP"] = "1"
    base = _square_planar_cn()
    moved = list(base)
    moved[5] = ("C", 0.3, 1.9, 0.0)
    moved[6] = ("N", 0.45, 3.0, 0.0)
    iso = [(_xyz(base), "a"), (_xyz(moved), "b")]
    out = pd.dedup_ensemble(iso)
    assert len(out) == 2


def test_deterministic():
    os.environ["DELFIN_FFFREE_PERMUTE_DEDUP"] = "1"
    frame = _xyz(_square_planar_cn())
    iso = [(frame, "a"), (frame, "b"), (frame, "c")]
    assert pd.dedup_ensemble(iso) == pd.dedup_ensemble(iso)


def test_pair_api_true_dup_vs_distinct():
    a = _xyz(_square_planar_cn())
    assert pd.is_permutation_duplicate(a, a) is True
    base = _square_planar_cn()
    import math
    distinct = list(base)
    distinct[5] = ("C", 0.0, 1.9 * math.cos(0.6), 1.9 * math.sin(0.6))
    distinct[6] = ("N", 0.0, 3.05 * math.cos(0.6), 3.05 * math.sin(0.6))
    assert pd.is_permutation_duplicate(a, _xyz(distinct)) is False


def test_never_raises_on_garbage():
    os.environ["DELFIN_FFFREE_PERMUTE_DEDUP"] = "1"
    iso = [("garbage", "x"), ("", "y")]
    out = pd.dedup_ensemble(iso)
    assert isinstance(out, list)
