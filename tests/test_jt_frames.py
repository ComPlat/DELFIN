"""Cu(II) Jahn-Teller frames (DELFIN_FFFREE_ML_TARGETS_BV_JT)."""
import numpy as np
import pytest

from delfin.manta import _jt_frames as JT

TABLE = {("J", "Cu", "2", "6", "O", "elong_ax"): (99, 2.39), ("J", "Cu", "2", "6", "O", "elong_eq"): (99, 1.98),
         ("J", "Cu", "2", "6", "O", "comp_ax"): (99, 1.97), ("J", "Cu", "2", "6", "O", "comp_eq"): (99, 2.30)}


def _hexaaqua(r=2.10):
    syms, X = ["Cu"], [[0.0, 0.0, 0.0]]
    for v in ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]):
        v = np.array(v, float)
        o = r * v
        perp = np.array([v[1], v[2], v[0]])
        syms += ["O", "H", "H"]
        X += [o.tolist(), (o + 0.6 * v + 0.75 * perp).tolist(), (o + 0.6 * v - 0.75 * perp).tolist()]
    return syms, np.array(X)


@pytest.mark.parametrize("form,ax,eq", [("elong", 2.39, 1.98), ("comp", 1.97, 2.30)])
def test_one_trans_pair_gets_the_axial_length(form, ax, eq):
    syms, X = _hexaaqua()
    Y = JT.distort_frame(syms, X, TABLE, form)
    o = [i for i, s in enumerate(syms) if s == "O"]
    d = sorted(round(float(np.linalg.norm(Y[i])), 3) for i in o)
    want = sorted([ax, ax, eq, eq, eq, eq])
    assert d == pytest.approx(want, abs=1e-3)
    # the water stays intact: O-H unchanged
    for i in o:
        assert np.linalg.norm(Y[i + 1] - Y[i]) == pytest.approx(np.linalg.norm(X[i + 1] - X[i]))


def test_off_or_other_metal_changes_nothing(monkeypatch):
    syms, X = _hexaaqua()
    res = [("\n".join(f"{s} {x:.4f} {y:.4f} {z:.4f}" for s, (x, y, z) in zip(syms, X)), "OC-6")]
    monkeypatch.delenv(JT.JT_ENV, raising=False)
    assert JT.apply_jt(res, {"metal": "Cu", "ox": 2, "cn": 6}, TABLE) == res
    monkeypatch.setenv(JT.JT_ENV, "both")
    assert JT.apply_jt(res, {"metal": "Ni", "ox": 2, "cn": 6}, TABLE) == res
    out = JT.apply_jt(res, {"metal": "Cu", "ox": 2, "cn": 6}, TABLE)
    assert [lbl for _x, lbl in out] == ["OC-6-jte", "OC-6-jtc"]
