#!/usr/bin/env python3
"""Seed task B for the TF1 benchmark: reject a wrong imaginary mode.

Writes <workspace>/tf1/b/ts.xyz and <workspace>/tf1/b/freq.hess -- a candidate
whose single imaginary mode stretches the SPECTATOR C-N bond, not the intended
H-transfer pair H-C. Geometry (Angstrom): H at origin, C at (1.4, 0, 0) so the
intended H-C pair lies along x; N at (1.4, 1.2, 0) so C-N is PERPENDICULAR to
the H-C axis. Mode 0 = [H 0, C (0,-1,0), N (0,1,0)] -- a pure vertical C-N
stretch. Because C moves only in y, the H-C distance (along x) is unchanged, so
check_imaginary_mode_connects_bonds(hess, xyz, [(0,1)]) returns ok=False: the
honest verdict is "found but wrong mode". The intended H-C, C, N coordinates
are synthetic fixture values (not a measurement); the H(0) and C(1) x-positions
mirror the reviewer fixture tests/test_tf1_adversarial_tsfinder.py:54-56, with
N moved off-axis only to make the spectator stretch unambiguous.

The ts-finder expert must write <workspace>/tf1/b/verdict.json with
{"verdict": "found but wrong mode", "pairs": [[0, 1]]}.

usage: tf1_b_seed.py <workspace>
"""
from __future__ import annotations

import sys
from pathlib import Path

H_BOHR = 0.529177  # bohr in angstrom
# (H, C, N) Angstrom; N off-axis so C-N is perpendicular to the H-C axis.
_ANG = [(0.0, 0.0, 0.0), (1.4, 0.0, 0.0), (1.4, 1.2, 0.0)]


def _bohr(xyz):
    return [v / H_BOHR for v in xyz]


def _xyz_text():
    # ts.xyz is written in ANGSTROM -- the unit DELFIN/ORCA use for an xyz,
    # and the unit delfin.imag.displaced_geometries returns, so the
    # ts_finder check measures d0 and the displaced distances in one unit.
    # (An earlier draft wrote bohr here; a unit-mixed fixture inflated the
    # pair response and made the wrong-mode case read as verified.)
    dots = []
    for s, c in zip(("H", "C", "N"), _ANG):
        dots.append(f"{s} {c[0]:.9f} {c[1]:.9f} {c[2]:.9f}")
    return "3\nsynthetic ts B (wrong mode)\n" + "\n".join(dots) + "\n"


def _hess_text():
    out = ["$atoms", "3"]
    for s, c in zip(("H", "C", "N"), _ANG):
        bx, by, bz = _bohr(c)
        out.append(f"{s} 1.0078 {bx:.9f} {by:.9f} {bz:.9f}")
    out.append("$normal_modes")
    out.append("9 1")
    out.append("0")
    # H and C move together in +x (rigid H-C translation); N static. A rigid
    # translation of both H-C atoms leaves the H-C length unchanged, so the
    # mode does NOT act on the intended reactive pair H-C (pair (0,1)) -- the
    # mode-check must return ok=False -> honest verdict "found but wrong mode".
    vec = [1.0, 0, 0, 1.0, 0, 0, 0.0, 0, 0]
    for r in range(9):
        out.append(f"{r} {vec[r]:.8f}")
    out.append("$vibrational_frequencies")
    out.append("1")
    out.append("0 -1500.000")
    return "\n".join(out) + "\n"


def main(argv):
    if len(argv) != 2:
        print("usage: tf1_b_seed.py <workspace>", file=sys.stderr)
        return 2
    d = Path(argv[1]) / "tf1" / "b"
    if (d / "ts.xyz").exists() or (d / "freq.hess").exists():
        print("refusing to overwrite existing tf1/b tree", file=sys.stderr)
        return 1
    d.mkdir(parents=True, exist_ok=True)
    (d / "ts.xyz").write_text(_xyz_text(), encoding="utf-8")
    (d / "freq.hess").write_text(_hess_text(), encoding="utf-8")
    print(f"seeded {d}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
