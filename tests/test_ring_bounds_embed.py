"""Ring-interior-angle DG bound for the metallacycle embed (DELFIN_FFFREE_RING_BOUNDS).

Eye-find QEBLOC (2026-06-24): the κ4 cage cap's rigid aromatic 5-ring embedded with
its interior angle OPEN (~115-126deg vs ideal ~108deg).  ``_tighten_ring_bounds``
adds an interior-angle 1-3 bound (108deg for 5-rings, 120deg for 6-rings) using the
bounds matrix's own 1-2 midpoints, so the ring embeds flat-and-correct.
"""
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as DG

from delfin.fffree.assemble_complex import _ring_bounds_enabled, _tighten_ring_bounds


def test_flag_off_default():
    os.environ.pop("DELFIN_FFFREE_RING_BOUNDS", None)
    assert _ring_bounds_enabled() is False
    os.environ["DELFIN_FFFREE_RING_BOUNDS"] = "1"
    assert _ring_bounds_enabled() is True
    os.environ.pop("DELFIN_FFFREE_RING_BOUNDS", None)


def _implied_angle(bm, i, j, k):
    def mid(a, b):
        lo, hi = (a, b) if a < b else (b, a)
        return 0.5 * (bm[lo][hi] + bm[hi][lo])
    d12, d23, d13 = mid(i, j), mid(j, k), mid(i, k)
    c = (d12 * d12 + d23 * d23 - d13 * d13) / (2 * d12 * d23)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def test_ring_bounds_encode_108_for_5ring():
    # imidazole: aromatic 5-ring
    m = Chem.AddHs(Chem.MolFromSmiles("c1c[nH]cn1"))
    bm = DG.GetMoleculeBoundsMatrix(m)
    ring = next(r for r in m.GetRingInfo().AtomRings() if len(r) == 5)
    # before: the embed-derived 1-3 angle may be loose
    _tighten_ring_bounds(bm, m)
    # after: every consecutive ring 1-3 triple encodes ~108 deg (regular pentagon)
    for a in range(5):
        ang = _implied_angle(bm, ring[a], ring[(a + 1) % 5], ring[(a + 2) % 5])
        assert abs(ang - 108.0) < 4.0, f"ring 1-3 angle {ang:.1f} not ~108"


def test_ring_bounds_noop_on_acyclic():
    # no ring -> nothing to tighten, must not raise / not corrupt the matrix
    m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    bm = DG.GetMoleculeBoundsMatrix(m)
    before = np.array(bm, copy=True)
    _tighten_ring_bounds(bm, m)
    assert np.allclose(np.array(bm), before)
