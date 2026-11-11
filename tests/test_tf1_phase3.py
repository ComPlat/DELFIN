"""TF1 phase 3: the mode-connects-bonds check must renumber endpoint-order
reactive pairs into the TS atom order (via atom_mapping.map_atoms) before
scoring, so a "found and verified" verdict never rests on a wrong atom index
(the reviewer's pin: renumber endpoints into TS order BEFORE any
bond-connectivity claim, else a wrong mode reads as verified).

The fixture reuses the collinear H...C-N geometry and Hessian format of the
reviewer's battery (tests/test_tf1_adversarial_tsfinder.py:54-56: H 0.5,
C 1.4, N 2.6 Angstrom; mode 0 = [H +0.5, C -0.5, N 0], $vibrational_frequencies
mode 0 = -1500). All coordinates here are SYNTHETIC fixture values, not
measurements.

Red control: the current check_imaginary_mode_connects_bonds has no
endpoint_xyz parameter, so the equivalence call below raises TypeError.
Green after phase 3 adds the optional endpoint_xyz renumbering path.
"""
import pathlib
import tempfile

from delfin.ts_finder import check_imaginary_mode_connects_bonds

H_BOHR = 0.529177  # bohr in angstrom


def _bohr(ang):
    return [v / H_BOHR for v in ang]


_TS_SYMS = ["H", "C", "N"]
# H 0.5, C 1.4, N 2.6 A on x -- synthetic, from the reviewer fixture (line 56).
_TS_ANG = [[0.5, 0.0, 0.0], [1.4, 0.0, 0.0], [2.6, 0.0, 0.0]]
# Same three atoms in a DIFFERENT (endpoint) order: C, N, H.
_END_SYMS = ["C", "N", "H"]
_END_ANG = [[1.4, 0.0, 0.0], [2.6, 0.0, 0.0], [0.5, 0.0, 0.0]]

# Reactive pair H-C given in ENDPOINT numbering (C index0, H index2 -> (2, 0));
# in TS numbering it is (0, 1).
ENDPOINT_PAIR = [(2, 0)]
TS_PAIR = [(0, 1)]


def _xyz(syms, coords_ang):
    lines = [str(len(syms)), "synthetic"]
    for s, c in zip(syms, coords_ang):
        lines.append(f"{s} {c[0]:.9f} {c[1]:.9f} {c[2]:.9f}")
    return "\n".join(lines)


TS_XYZ = _xyz(_TS_SYMS, _TS_ANG)
END_XYZ = _xyz(_END_SYMS, _END_ANG)


def _hess():
    rows = 9
    atoms = list(zip(_TS_SYMS, _TS_ANG))
    out = ["$atoms", "3"]
    for el, c in atoms:
        b = _bohr(c)
        out.append(f"{el} 1.0078 {b[0]:.9f} {b[1]:.9f} {b[2]:.9f}")
    vec = {0: [0.5, 0, 0, -0.5, 0, 0, 0, 0, 0]}
    out.append("$normal_modes")
    out.append(f"{rows} 1")
    out.append("0")
    for r in range(rows):
        out.append(f"{r} {vec[0][r]:.8f}")
    out.append("$vibrational_frequencies")
    out.append("1")
    out.append("0 -1500.000")
    d = pathlib.Path(tempfile.mkdtemp(prefix="tf1_ph3_"))
    p = d / "freq.hess"
    p.write_text("\n".join(out), encoding="utf-8")
    return str(p)


def test_renumbered_endpoint_pairs_match_the_ts_order_verdict():
    hess = _hess()
    # TS-order pair (0,1) = H-C: the mode stretches it -> verified.
    ts_order = check_imaginary_mode_connects_bonds(
        hess, TS_XYZ, TS_PAIR, mode_index=0)
    assert ts_order.ok, ts_order.problems
    # ENDPOINT-order pair (2,0) = H-C in [C,N,H] numbering must be renumbered
    # to (0,1) via map_atoms(endpoint -> TS) before scoring; it must then
    # give the SAME verdict. Current code has no endpoint_xyz -> TypeError.
    endpoint_order = check_imaginary_mode_connects_bonds(
        hess, TS_XYZ, ENDPOINT_PAIR, mode_index=0, endpoint_xyz=END_XYZ)
    assert endpoint_order.ok, endpoint_order.problems
    assert endpoint_order.ok == ts_order.ok
