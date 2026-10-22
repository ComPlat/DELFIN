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


def test_a_chelate_donor_is_rotated_about_its_partner():
    syms, X = _hexaaqua()
    # replace the +x and +y waters by one N-C-C-N chelate (N at 2.10 A on +x and +y)

    syms = ["Cu", "N", "C", "C", "N"]
    Xn = [[0, 0, 0], [2.10, 0, 0], [2.55, 1.30, 0], [1.30, 2.55, 0], [0, 2.10, 0]]
    for v in ([-1, 0, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]):
        v = np.array(v, float)
        o = 2.10 * v
        perp = np.array([v[1], v[2], v[0]])
        syms += ["O", "H", "H"]
        Xn += [o.tolist(), (o + 0.6 * v + 0.75 * perp).tolist(), (o + 0.6 * v - 0.75 * perp).tolist()]
    X = np.array(Xn, float)
    table = dict(TABLE)
    table[("J", "Cu", "2", "6", "N", "elong_ax")] = (99, 2.50)      # makes the N-O axis the best
    table[("J", "Cu", "2", "6", "N", "elong_eq")] = (99, 2.02)
    Y = JT.distort_frame(syms, X, table, "elong")
    assert Y is not None
    d = {i: float(np.linalg.norm(Y[i])) for i in (1, 4)}
    assert sorted(round(v, 2) for v in d.values()) == [2.1, 2.5]
    # the chelate moved rigidly: N-C-C-N bond lengths unchanged
    for a, b in ((1, 2), (2, 3), (3, 4)):
        assert np.linalg.norm(Y[a] - Y[b]) == pytest.approx(np.linalg.norm(X[a] - X[b]))


def test_all_keeps_the_frame_as_built(monkeypatch):
    syms, X = _hexaaqua()
    res = [("\n".join(f"{s} {x:.4f} {y:.4f} {z:.4f}" for s, (x, y, z) in zip(syms, X)), "OC-6")]
    monkeypatch.setenv(JT.JT_ENV, "all")
    out = JT.apply_jt(res, {"metal": "Cu", "ox": 2, "cn": 6}, TABLE)
    assert [lbl for _x, lbl in out] == ["OC-6", "OC-6-jte", "OC-6-jtc"] and out[0] == res[0]


def test_square_pyramid_gets_a_long_apical_bond():
    syms, X = _hexaaqua()
    keep = [i for i in range(len(syms)) if not (16 <= i <= 18)]      # drop the -z water
    syms = [syms[i] for i in keep]
    X = X[keep]
    table = {("J", "Cu", "2", "5", "O", "sp5_ap"): (99, 2.29), ("J", "Cu", "2", "5", "O", "sp5_bas"): (99, 1.95)}
    Y = JT.distort_frame(syms, X, table, "elong")
    o = [i for i, s in enumerate(syms) if s == "O"]
    d = sorted(round(float(np.linalg.norm(Y[i])), 3) for i in o)
    assert d == pytest.approx([1.95, 1.95, 1.95, 1.95, 2.29], abs=1e-3)
    assert JT.distort_frame(syms, X, table, "comp") is None
