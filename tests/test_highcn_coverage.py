"""Tests for high-CN (CN7-9) coverage in the fffree builder (Option A).

Validates: ideal polyhedra (PB-7 / SQAP-8 / TTP-9) are unit-vector sets of the
right size; the polya proper-rotation groups close to the correct orders
(D5=10, D4=8, D3=6) and every element is a real geometric rotation
(distance-preserving on the vertex set); homoleptic isomer count = 1; the
decompose CN gate is env-gated (DELFIN_FFFREE_HIGHCN, default OFF = byte-identical
CN 4-6 only); and a low-denticity CN7 complex builds when the flag is on.
"""
import math
import os

import numpy as np
import pytest

from delfin.fffree import polyhedra as PH
from delfin.fffree import polya_isomer_count as PIC

_NEW = [
    ("PB-7 pentagonal bipyramid", "pentagonal_bipyramid", 7, 10),
    ("SQAP-8 square antiprism", "square_antiprism", 8, 8),
    ("TTP-9 tricapped trigonal prism", "tricapped_trigonal_prism", 9, 6),
]


@pytest.mark.parametrize("geom,key,cn,order", _NEW)
def test_ref_vectors_unit_and_count(geom, key, cn, order):
    v = PH.ref_vectors(geom)
    assert v.shape == (cn, 3)
    assert np.allclose(np.linalg.norm(v, axis=1), 1.0)


@pytest.mark.parametrize("geom,key,cn,order", _NEW)
def test_group_order(geom, key, cn, order):
    grp, n = PIC._GROUPS[key]
    assert n == cn
    assert len(grp) == order, f"{key} group order {len(grp)} != {order}"


@pytest.mark.parametrize("geom,key,cn,order", _NEW)
def test_group_elements_are_real_rotations(geom, key, cn, order):
    """Every group permutation must preserve all pairwise vertex distances
    (i.e. be an actual rotation of the polyhedron) — else isomer dedup is wrong."""
    grp, _ = PIC._GROUPS[key]
    V = PH.ref_vectors(geom)
    D = np.linalg.norm(V[:, None, :] - V[None, :, :], axis=-1)
    for p in grp:
        assert np.allclose(D[np.ix_(p, p)], D), f"{key} perm {p} not distance-preserving"


@pytest.mark.parametrize("geom,key,cn,order", _NEW)
def test_homoleptic_one_isomer(geom, key, cn, order):
    assert PIC.count_isomers(key, {"N": cn}) == 1


def test_geom_by_cn_has_high_cn():
    for cn in (7, 8, 9):
        assert cn in PH.GEOM_BY_CN and PH.GEOM_BY_CN[cn]


def _cn7_smiles():
    # [ZrF7]3- : 7 monodentate fluorides -> CN7, denticity 1 each (passes the <=3 gate)
    return "F[Zr-3](F)(F)(F)(F)(F)F"


def test_decompose_cn_gate_env_gated():
    from delfin.fffree import decompose as D
    smi = _cn7_smiles()
    os.environ.pop("DELFIN_FFFREE_HIGHCN", None)
    assert D.decompose(smi) is None, "default (HIGHCN off) must reject CN7 (byte-identical)"
    os.environ["DELFIN_FFFREE_HIGHCN"] = "1"
    try:
        d = D.decompose(smi)
        assert d is not None and d["cn"] == 7
        assert d["geometry"] == "PB-7 pentagonal bipyramid"
    finally:
        os.environ.pop("DELFIN_FFFREE_HIGHCN", None)


def test_cn7_builds_when_enabled():
    from delfin.fffree import converter_backend as CB
    os.environ["DELFIN_FFFREE_HIGHCN"] = "1"
    try:
        res = CB._fffree_isomers(_cn7_smiles(), max_isomers=4)
    finally:
        os.environ.pop("DELFIN_FFFREE_HIGHCN", None)
    assert res is not None and len(res) >= 1, "low-denticity CN7 must build under HIGHCN"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
