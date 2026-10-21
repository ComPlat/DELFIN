"""Derivation of the bond-valence M-D tables (delfin.manta.ml_targets_derive, tools/derive_ml_targets.py)."""
import math
from collections import Counter

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
pytest.importorskip("networkx")

from delfin.manta import _ml_targets_bv as BV  # noqa: E402
from delfin.manta import ml_targets_derive as MD  # noqa: E402

ZN_FORMATE = "[Cl][Zn-2]([Cl])([Cl])[O]C=O"


def _zn_formate_crystal(swap):
    o1, o2 = [0.0, 0.0, 1.95], [0.0, 2.2, 2.2]          # bound O at 1.95, free O at 3.11 A
    if swap:
        o1, o2 = o2, o1
    el = ["Zn", "Cl", "Cl", "Cl", "O", "C", "O", "H"]
    X = [[0, 0, 0], [2.25, 0, 0], [-1.125, 1.95, -0.6], [-1.125, -1.95, -0.6],
         o1, [0.0, 1.05, 2.65], o2, [0.0, 1.0, 3.75]]
    return el, np.array(X, float)


@pytest.mark.parametrize("swap", [False, True])
def test_a_symmetric_donor_takes_the_bound_crystal_atom(swap):
    el, X = _zn_formate_crystal(swap)
    r = MD.record_from_smiles_crystal("TESTZN", ZN_FORMATE, el, X)
    assert r["st"] == "ok", r
    assert (r["metal"], r["ox"], r["cn"], r["nx"]) == ("Zn", 2, 4, 4)
    o = [d for d in r["donors"] if d[0] == "O"]
    assert len(o) == 1 and o[0][2] == pytest.approx(1.95, abs=1e-3)


def test_a_donor_not_bonded_in_the_crystal_refuses_the_structure():
    el, X = _zn_formate_crystal(False)
    X[4] = [0.0, 0.0, 2.9]                               # the only candidate O is now 2.9 A away
    X[5] = [0.0, 1.05, 3.6]
    X[6] = [0.0, 2.2, 3.15]
    X[7] = [0.0, 1.0, 4.7]
    r = MD.record_from_smiles_crystal("TESTZN", ZN_FORMATE, el, X)
    assert r["st"] == "map_donor_not_bonded"


def test_complex_charge_from_the_formula_unit():
    cx = Counter({"Fe": 1, "C": 30, "H": 24, "N": 6})
    fs = Counter({"Fe": 1, "C": 30, "H": 24, "N": 6, "Cl": 2, "O": 8})         # 2 ClO4-
    assert MD.complex_charge_from_formula(fs, cx) == 2
    fs2 = Counter({"Fe": 2, "C": 60, "H": 52, "N": 12, "Cl": 4, "O": 18})       # Z' = 2, + 2 H2O
    assert MD.complex_charge_from_formula(fs2, cx) == 2
    assert MD.complex_charge_from_formula(Counter(cx), cx) == 0
    assert MD.complex_charge_from_formula(Counter({"Fe": 1, "C": 30, "H": 25, "N": 6}), cx) is None


def test_consensus_needs_a_bond_valence_sum_within_the_window():
    assert MD.consensus_ok({"ox": 2, "bvs": 2.4}, 0.5)
    assert not MD.consensus_ok({"ox": 2, "bvs": 2.6}, 0.5)
    assert not MD.consensus_ok({"ox": 0, "bvs": 0.2}, 0.5)
    assert not MD.consensus_ok({"ox": 2, "bvs": None}, 0.5)


def test_the_fit_recovers_b_and_the_valence_split():
    b0, rho0, r0x, r0l = 0.20, 2.0, 1.90, 1.80
    recs = []
    k = 0
    for cn in (4, 5, 6):
        for nx in range(0, cn + 1):
            for rep in range(40):
                tot = rho0 * nx + (cn - nx)
                don = [["Cl", "t0x", round(r0x - b0 * math.log(2 * rho0 / tot), 4), None, None]] * nx
                don += [["N", "pil", round(r0l - b0 * math.log(2 / tot), 4), None, None]] * (cn - nx)
                recs.append({"id": f"S{k}", "metal": "Zn", "ox": 2, "cn": cn, "nx": nx, "donors": don})
                k += 1
    f = MD.fit_bv(recs)
    assert f["b"] == pytest.approx(b0, abs=0.026) and f["rho"] == pytest.approx(rho0, abs=0.26)
    assert f["score"] < 0.005 < f["literature_score"]


def test_the_x_count_separates_a_homoleptic_anion(tmp_path, monkeypatch):
    rows = [("N", "Zn", 2, "4/4", "Br", "t0x", 99, 2.415), ("N", "Zn", 2, "4/2", "Br", "t0x", 99, 2.36),
            ("D", "Zn", 2, 4, "Br", "t0x", 99, 2.362), ("D", "Zn", 2, 4, "N", "pil", 99, 2.06)]
    p = tmp_path / "t.tsv"
    p.write_text("\n".join("\t".join(str(x) for x in r) for r in rows) + "\n")
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setenv(BV.TABLE_ENV, str(p))
    BV._TABLE_CACHE.clear()
    BV.set_context("[Br][Zn-2]([Br])([Br])[Br]")
    assert BV.target_for_donor("Zn", "Br") == pytest.approx(2.415)
    BV.set_context("[Br][Zn-2]([Br])([N+]1=CC=CC=C1)[N+]1=CC=CC=C1")
    assert BV.target_for_donor("Zn", "Br") == pytest.approx(2.36)
    BV.clear_context()
