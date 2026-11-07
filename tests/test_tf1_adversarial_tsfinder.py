"""TF1 (ts-finder) review battery, part 3: attack the ts_finder library.

Adversarial tests that drive the REAL functions of delfin.ts_finder (from
the builder's phase-1 commit f75fb82d6) through edge cases and both
polarities. Two groups:

1. NEVER RAISES on hostile input -- a pure check that raises on a crafted
   (but in-parser-format) input would break the expert mid-search and leak
   a traceback instead of an honest verdict. Every public check must return
   a CheckResult.
2. HONEST POLARITY -- each check must agree with a physically-built input,
   not just not-crash: a mode that stretches the forming/breaking pair
   passes, a mode that leaves that pair's length ~fixed fails (this is the
   "found but wrong mode" honesty), a mixed level of theory is rejected,
   and summarize_ts_verdict emits exactly the three honest labels.

Synthetic inputs are built in the exact format the DELFIN readers use
($atoms/$normal_modes Hessian with Bohr coordinates, xyz blocks, dicts),
never guesswork -- see delfin/imag.py displaced_geometries/_normal_mode.
"""
from __future__ import annotations

import importlib

# delfin.ts_finder imports map_atoms/ts_goat/imag; import the library and
# bind a short alias so the tests read as the product calls them.
from delfin.ts_finder import (
    CheckResult,
    VERDICT_NOT_FOUND,
    VERDICT_WRONG_MODE,
    VERDICT_VERIFIED,
    check_endpoint_mapping,
    check_irc_ends_are_minima,
    check_imaginary_mode_connects_bonds,
    check_single_level_of_theory,
    summarize_ts_verdict,
)


H = 0.529177  # bohr in angstrom, _BOHR_IN_ANGSTROM


def bohr(x_ang):
    return [v / H for v in x_ang]


def _xyz(syms, coords_ang):
    lines = [str(len(syms)), "synthetic"]
    for s, c in zip(syms, coords_ang):
        lines.append(f"{s} {c[0]:.9f} {c[1]:.9f} {c[2]:.9f}")
    return "\n".join(lines)


# A collinear 3-atom H...C-N transition-state-like geometry (Angstrom).
_TS_SYMS = ["H", "C", "N"]
_TS_ANG = [[0.5, 0.0, 0.0], [1.4, 0.0, 0.0], [2.6, 0.0, 0.0]]
TS_XYZ = _xyz(_TS_SYMS, _TS_ANG)


def _hessian(mode_rows, cols=None):
    """A Hessian text in imag.py's $atoms/$normal_modes format.

    mode_rows: list of mode index -> 3N cartesian displacement (Angstrom
    units for the vector is fine: only the relative shape matters, the code
    re-normalises by the largest element). Rows are printed in column
    blocks of <=5 like the real writer.
    """
    n_atoms = len(_TS_SYMS)
    rows = 3 * n_atoms
    cols = cols or rows
    atoms = [("H", bohr([0.5, 0.0, 0.0])), ("C", bohr([1.4, 0.0, 0.0])),
             ("N", bohr([2.6, 0.0, 0.0]))]
    out = ["$atoms", f"{n_atoms}"]
    for el, (x, y, z) in atoms:
        out.append(f"{el} 1.0078 {x:.9f} {y:.9f} {z:.9f}")
    # full zero vectors for modes not supplied
    vec = {m: [0.0] * rows for m in range(cols)}
    for m, v in mode_rows.items():
        vec[m] = list(v) + [0.0] * (rows - len(v))
    out.append("$normal_modes")
    out.append(f"{rows} {cols}")
    mode_list = list(range(cols))
    for block_start in range(0, cols, 5):
        block = mode_list[block_start:block_start + 5]
        out.append(" ".join(str(m) for m in block))
        for r in range(rows):
            out.append(f"{r} " + " ".join(f"{vec[m][r]:.8f}" for m in block))
    # The $vibrational_frequencies of the same Hessian: mode 0 the reactive
    # imaginary mode, the rest numerical-real values above the noise floor.
    out.append("$vibrational_frequencies")
    out.append(f"{cols}")
    for m in range(cols):
        freq = -1500.0 if m == 0 else 200.0 + 100.0 * m
        out.append(f"{m} {freq:.3f}")
    return "\n".join(out)


def _hess(mode_rows):
    import pathlib, tempfile
    d = pathlib.Path(tempfile.mkdtemp(prefix="tf1_adv_"))
    p = d / "freq.hess"
    p.write_text(_hessian(mode_rows), encoding="utf-8")
    return str(p)


# ---------------------------------------------------------------------------
# 1. Never raises: hostile inputs to every public check.
# ---------------------------------------------------------------------------

def test_never_raises_on_hostile_level_dict():
    # {"reactant": "ß🧪"}: a stage whose value is unicode. {"a": 1, "b": 2}
    # maps numbers -- wrong shape entirely. None and a stage-less dict are
    # the empty/None edge; a single stage is a partial barrier.
    for bad in (None, {}, {"reactant": "x"},
                {"reactant": "A", "product": "A", "ts": "A", "zpe": "A",
                 "solvent\N{COMBINING LONG SOLIDUS OVERLAY}": "ß\U0001F9EA"},
                {"a": 1, "b": 2}):
        r = check_single_level_of_theory(bad)
        assert isinstance(r, CheckResult), (bad, r)


def test_never_raises_on_hostile_endpoint_blocks():
    for a, b in ((None, "x"), ("x", None), ("", ""), ("garbage", "garbage"),
                 ("3\nno\nH a b c", _xyz(["H"], [[0, 0, 0]])),
                 ("abc\ndef", "ghi")):
        r = check_endpoint_mapping(a, b)
        assert isinstance(r, CheckResult), (a, b, r)


def test_never_raises_on_hostile_irc_counts():
    for bad in (None, {}, {"forward": 0}, {"forward": "x"}, {"f": -3},
                {"f": None}, {"f": 0, "r": None, "x": "y"}, [1, 2]):
        r = check_irc_ends_are_minima(bad)
        assert isinstance(r, CheckResult), (bad, r)


def test_never_raises_on_hostile_mode_inputs():
    # Missing hess, unparseable xyz, out-of-bounds pair, self pair.
    for args in ((_hess({0: [0.5, 0, 0, -0.5, 0, 0, 0, 0, 0]}), "garbage",
                  [(0, 1)]),
                 ("/nonexistent/freq.hess", TS_XYZ, [(0, 1)]),
                 (_hess({0: [0.5, 0, 0, -0.5, 0, 0, 0, 0, 0]}), TS_XYZ,
                  [(0, 37)]),
                 (_hess({0: [0.5, 0, 0, -0.5, 0, 0, 0, 0, 0]}), TS_XYZ,
                  [(1, 1)])):
        r = check_imaginary_mode_connects_bonds(*args)
        assert isinstance(r, CheckResult), args


# ---------------------------------------------------------------------------
# 2. Honest polarity.
# ---------------------------------------------------------------------------

def test_single_level_accepts_only_one_matching_line():
    same = {"reactant": "PBE0/def2-SVP", "ts": "PBE0/def2-SVP",
            "product": "PBE0/def2-SVP", "zpe": "PBE0/def2-SVP"}
    assert check_single_level_of_theory(same).ok
    mixed = dict(same, product="BP86/def2-SVP")
    assert not check_single_level_of_theory(mixed).ok
    assert not check_single_level_of_theory({}).ok


def test_mode_stretching_the_pair_is_verified():
    # Mode 0: H and C move in opposite x directions -> stretches H-C (pair
    # 0-1) strongly; N static. The mode moves atom 1 (C), so it acts on
    # C-N (1,2) as well -- the old "not ok" expectation for (1,2) was
    # wrong. A pair is genuinely untouched only when its two atoms move
    # rigidly together; probe that rejection branch with a rigid mode:
    # H and C both +x leaves the H-C distance constant.
    hess = _hess({0: [0.5, 0, 0, -0.5, 0, 0, 0, 0, 0]})
    r_ok = check_imaginary_mode_connects_bonds(hess, TS_XYZ, [(0, 1)])
    assert r_ok.ok, r_ok
    rigid = _hess({0: [0.5, 0, 0, 0.5, 0, 0, 0, 0, 0]})
    r_rigid = check_imaginary_mode_connects_bonds(rigid, TS_XYZ, [(0, 1)])
    assert not r_rigid.ok, r_rigid


def test_mode_moving_perpendicular_is_not_verified():
    # Mode moves only H along y (perpendicular to the H-C axis): the H-C
    # distance barely changes, so a fixed contact should NOT count.
    hess = _hess({0: [0, 0.5, 0, 0, 0, 0, 0, 0, 0]})
    r = check_imaginary_mode_connects_bonds(hess, TS_XYZ, [(0, 1)])
    assert not r.ok, r


def test_irc_counts_require_zero_imaginary_everywhere():
    assert check_irc_ends_are_minima(
        {"forward": 0, "reverse": 0}).ok
    assert not check_irc_ends_are_minima(
        {"forward": 0, "reverse": 1}).ok
    assert not check_irc_ends_are_minima({}).ok


def test_summarize_verdict_is_honest_and_stable():
    assert summarize_ts_verdict(False, True, True) == "not found"
    assert summarize_ts_verdict(True, False, True) == "found but wrong mode"
    assert summarize_ts_verdict(True, True, False) == "found but wrong mode"
    assert summarize_ts_verdict(True, True, True) == "found and verified"
    # exact stable labels, never rounded up
    assert VERDICT_NOT_FOUND == "not found"
    assert VERDICT_WRONG_MODE == "found but wrong mode"
    assert VERDICT_VERIFIED == "found and verified"
