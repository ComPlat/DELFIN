"""The bond-valence chelate embed keeps the M-D radii AND the chemistry of the ring closure
(DELFIN_FFFREE_BV_EMBED_GEOM, 30.10.2026: lc18 FF-free chelate frames 85 % defective, AFODIJ bipy twisted
79 deg and Ru 58 deg out of the pyridine plane)."""
import numpy as np
import pytest

pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402

from delfin.manta import _ml_targets_bv as BVT  # noqa: E402
from delfin.manta import assemble_ligand_embed as LE  # noqa: E402
from delfin.manta import bv_embed_geom as G  # noqa: E402

BIPY = "c1ccc(-c2ccccn2)nc1"
R = 2.06


def _targets(n_don, r=R):
    v = [np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0])]
    return [r * v[i] for i in range(n_don)]


def _bv_ring(smi, donors, monkeypatch, geom):
    monkeypatch.setattr(BVT, "embed_at_targets", lambda: True)
    if geom:
        monkeypatch.setenv("DELFIN_FFFREE_BV_EMBED_GEOM", "1")
    else:
        monkeypatch.delenv("DELFIN_FFFREE_BV_EMBED_GEOM", raising=False)
    m = Chem.MolFromSmiles(smi)
    out = LE._embed_metallacycle(m, donors, "Ru", k=6, donor_target_pos=_targets(len(donors)))
    assert out is not None
    lsyms, coords = out
    return m, coords


def _plane_dist(P, idx, point):
    Q = P[idx]
    c = Q.mean(axis=0)
    _, _, vt = np.linalg.svd(Q - c)
    return float(abs((point - c) @ vt[2])), float(np.abs((Q - c) @ vt[2]).max())


def _bipy_geometry(m, coords):
    ri = m.GetRingInfo()
    rings = [r for r in ri.AtomRings() if len(r) == 6]
    don = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
    c2 = [b for b in m.GetBonds() if not b.IsInRing()][0]
    a, b = c2.GetBeginAtomIdx(), c2.GetEndAtomIdx()
    na = [n.GetIdx() for n in m.GetAtomWithIdx(a).GetNeighbors() if n.GetSymbol() == "N"][0]
    nb = [n.GetIdx() for n in m.GetAtomWithIdx(b).GetNeighbors() if n.GetSymbol() == "N"][0]
    res = []
    for P in coords:
        tors = abs(G._torsion_deg(P[na], P[a], P[b], P[nb]))
        mdist = max(_plane_dist(P, list(r), np.zeros(3))[0] for r in rings)
        roop = max(_plane_dist(P, list(r), np.zeros(3))[1] for r in rings)
        md = [float(np.linalg.norm(P[d])) for d in don]
        res.append((tors, mdist, roop, md))
    return res


def test_donor_kind_reads_the_bonds():
    m = Chem.AddHs(Chem.MolFromSmiles("C(C#N)c1ccncc1"))
    rw = Chem.RWMol(m)
    mi = rw.AddAtom(Chem.Atom("Ru"))
    npy = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N" and a.GetIsAromatic()][0]
    nni = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N" and not a.GetIsAromatic()][0]
    rw.AddBond(npy, mi, Chem.BondType.DATIVE)
    mh = rw.GetMol()
    assert G.donor_kind(mh, npy, mi)[0] == "sp2"
    assert G.donor_kind(mh, nni, mi)[0] == "sp"


def test_off_is_the_historic_path(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_BV_EMBED_GEOM", raising=False)
    assert not G.geom_enabled() and not G.use_mrow() and not G.use_etkdg() and not G.use_refine()


def test_bipy_flat_metal_in_plane_radii_kept(monkeypatch):
    m = Chem.MolFromSmiles(BIPY)
    don = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
    m, coords = _bv_ring(BIPY, don, monkeypatch, geom=True)
    geo = _bipy_geometry(m, coords)
    assert geo
    for tors, mdist, roop, md in geo:
        assert tors < 12.0                       # chelating bipy stays cis-planar
        assert mdist < 0.25                      # Ru in both pyridine planes
        assert roop < 0.06                       # aromatic rings flat
        assert all(abs(x - R) < 0.06 for x in md)  # bond-valence radius kept


def test_geom_beats_the_plain_bv_embed(monkeypatch):
    m = Chem.MolFromSmiles(BIPY)
    don = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
    m0, c0 = _bv_ring(BIPY, don, monkeypatch, geom=False)
    m1, c1 = _bv_ring(BIPY, don, monkeypatch, geom=True)
    g0, g1 = _bipy_geometry(m0, c0), _bipy_geometry(m1, c1)
    assert max(t for t, *_ in g1) <= max(t for t, *_ in g0) + 1e-6
    assert max(x for _, x, *_ in g1) <= max(x for _, x, *_ in g0) + 1e-6


def test_ethylenediamine_lone_pairs_toward_the_metal(monkeypatch):
    smi = "NCCN"
    don = [0, 3]
    m, coords = _bv_ring(smi, don, monkeypatch, geom=True)
    for P in coords:
        for d, c in ((0, 1), (3, 2)):
            u, v = -P[d], P[c] - P[d]
            ang = np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v)))
            assert 96.0 <= ang <= 123.0
            assert abs(np.linalg.norm(P[d]) - R) < 0.06


def test_deterministic(monkeypatch):
    m = Chem.MolFromSmiles(BIPY)
    don = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
    _, a = _bv_ring(BIPY, don, monkeypatch, geom=True)
    _, b = _bv_ring(BIPY, don, monkeypatch, geom=True)
    assert len(a) == len(b) and all(np.array_equal(x, y) for x, y in zip(a, b))
