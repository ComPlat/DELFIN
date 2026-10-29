"""Frames are emitted best first by a crystal-free realism score; nothing is dropped
(DELFIN_FRAME_REALISM_RANK, user rule 29.10.2026 "die besten, realistischsten vorne")."""
import json

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("networkx")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

from delfin.manta import frame_realism as FR  # noqa: E402

SMI = "Cl[Zn-2](Cl)(Cl)[N+]1=CC=CC=C1"


def _py_embed():
    m = Chem.AddHs(Chem.MolFromSmiles("c1ccncc1"))
    AllChem.EmbedMolecule(m, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(m)
    P = m.GetConformer().GetPositions()
    el = [a.GetSymbol() for a in m.GetAtoms()]
    n = el.index("N")
    nb = [a.GetIdx() for a in m.GetAtomWithIdx(n).GetNeighbors()]
    return el, P, n, nb


def _zn_py_frame(oop_deg):
    el, P, n, nb = _py_embed()
    u = [(P[i] - P[n]) / np.linalg.norm(P[i] - P[n]) for i in nb]
    lp = -(u[0] + u[1])
    lp /= np.linalg.norm(lp)
    nrm = np.cross(u[0], u[1])
    nrm /= np.linalg.norm(nrm)
    t = np.radians(oop_deg)
    d = np.cos(t) * lp + np.sin(t) * nrm
    zn = P[n] + 2.05 * d
    p1 = np.cross(d, nrm)
    p1 /= np.linalg.norm(p1)
    p2 = np.cross(d, p1)
    cl = []
    for k in range(3):
        ph = 2 * np.pi * k / 3
        v = np.cos(np.radians(70.5)) * d + np.sin(np.radians(70.5)) * (np.cos(ph) * p1 + np.sin(ph) * p2)
        cl.append(zn + 2.25 * v / np.linalg.norm(v))
    rows = [("Zn", zn)] + [("Cl", c) for c in cl] + list(zip(el, P))
    return "%d\nt\n" % len(rows) + "\n".join("%s %.6f %.6f %.6f" % (s, *p) for s, p in rows) + "\n"


@pytest.fixture
def table(tmp_path, monkeypatch):
    T = {"fine": {}, "coarse": {},
         "all": {q: {"n": 1000, "med": 0.0, "lo1": 0.0, "hi1": 5.0} for q in FR.CAT}}
    T["all"]["bond_len"] = {"n": 1000, "med": 1.0, "lo1": 0.9, "hi1": 1.1}
    T["all"]["ml_len"] = {"n": 1000, "med": 2.05, "lo1": 1.9, "hi1": 2.2}
    T["all"]["m_contact"] = {"n": 1000, "med": 2.0, "lo1": 1.4, "hi1": 3.0}
    T["all"]["nb_clash"] = {"n": 1000, "med": 0.9, "lo1": 0.7, "hi1": 1.0}
    p = tmp_path / "t.json"
    p.write_text(json.dumps(T))
    monkeypatch.setenv("DELFIN_FRAME_REALISM_TABLE", str(p))
    FR._TABLES.clear()
    FR._CHEM_CACHE.clear()
    return p


def test_off_is_identity(monkeypatch, table):
    monkeypatch.delenv("DELFIN_FRAME_REALISM_RANK", raising=False)
    iso = [(_zn_py_frame(60), "a"), (_zn_py_frame(0), "b")]
    assert FR.rank_isomers_realism(iso, SMI) is iso


def test_scores_see_the_metal_out_of_the_donor_plane(table):
    s = FR.score_frames([_zn_py_frame(0), _zn_py_frame(60)], SMI)
    assert s is not None and s[0] is not None and s[1] is not None
    assert s[1] > s[0]


def test_on_sorts_best_first_and_drops_nothing(monkeypatch, table):
    monkeypatch.setenv("DELFIN_FRAME_REALISM_RANK", "1")
    iso = [(_zn_py_frame(60), "bad"), ("1\nx\nZn 0 0 0\n", "unmappable"), (_zn_py_frame(0), "good"),
           (_zn_py_frame(25), "mid")]
    out = FR.rank_isomers_realism(iso, SMI)
    assert [l for _x, l in out] == ["good", "mid", "bad", "unmappable"]
    assert sorted(map(id, out)) == sorted(map(id, iso))


def test_ties_keep_the_emitted_order(monkeypatch, table):
    monkeypatch.setenv("DELFIN_FRAME_REALISM_RANK", "1")
    f = _zn_py_frame(0)
    iso = [(f, "first"), (f, "second")]
    assert [l for _x, l in FR.rank_isomers_realism(iso, SMI)] == ["first", "second"]


def test_no_table_is_a_no_op(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_FRAME_REALISM_RANK", "1")
    monkeypatch.setenv("DELFIN_FRAME_REALISM_TABLE", str(tmp_path / "missing.json"))
    FR._TABLES.clear()
    iso = [(_zn_py_frame(60), "a"), (_zn_py_frame(0), "b")]
    assert FR.rank_isomers_realism(iso, SMI) is iso


def test_mirror_visits_the_best_frame_first(monkeypatch, table):
    from delfin.manta import _mirror_enum as ME
    monkeypatch.setenv("DELFIN_FRAME_REALISM_RANK", "1")
    monkeypatch.setenv("DELFIN_MIRROR_ENUM", "1")
    monkeypatch.setenv("DELFIN_MIRROR_ONE_PER_SYSTEM", "1")
    monkeypatch.setenv("DELFIN_MIRROR_PER_ISOMER", "1")
    monkeypatch.delenv("DELFIN_MIRROR_STEREO_GATE", raising=False)
    monkeypatch.delenv("DELFIN_MIRROR_QUALITY_GATE", raising=False)
    seen = []
    real = ME.mirror_frame

    def spy(xyz):
        seen.append(xyz)
        return real(xyz)
    monkeypatch.setattr(ME, "mirror_frame", spy)
    bad, good = _zn_py_frame(60), _zn_py_frame(0)
    FR.set_context(SMI)
    try:
        ME.expand_results([(bad, "T-4-1"), (good, "T-4-1-conf2")])
    finally:
        FR.clear_context()
    assert seen and seen[0] == good
