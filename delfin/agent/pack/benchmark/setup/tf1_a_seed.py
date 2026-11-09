#!/usr/bin/env python3
"""Seed task A for the TF1 benchmark: pick the correct transition state.

Writes <workspace>/tf1/a/ts.xyz and <workspace>/tf1/a/freq.hess -- a single
candidate TS whose IMAGINARY mode (one, at -1500 cm-1) stretches exactly the
intended H-transfer pair H(0)..C(1). The atomic coordinates are the synthetic
collinear H...C-N fixture geometry from
tests/test_tf1_adversarial_tsfinder.py:54-56 (H at 0.5, C at 1.4, N at 2.6
Angstrom along x) -- fixture values, not a measurement. The Hessian uses the
same $atoms/$normal_modes/$vibrational_frequencies syntax that delfin.imag
reads (mode 0 vector [H +1, C -1, N 0] stretches pair (0,1)).

The ts-finder expert must verify this candidate -> "found and verified",
reactive pairs [[0, 1]], and record it in <workspace>/tf1/a/verdict.json.

usage: tf1_a_seed.py <workspace>
"""
from __future__ import annotations

import sys
from pathlib import Path

H_BOHR = 0.529177  # bohr in angstrom
_X = (0.5, 1.4, 2.6)  # Angstrom, from the reviewer fixture (tsfinder test)


def _bohr(ang):
    return [v / H_BOHR for v in ang]


def _xyz_text():
    # ts.xyz in ANGSTROM (matches displaced_geometries' unit; see the
    # unit-mixing note in tf1_b_seed.py).
    dots = [f"{s} {x:.9f} 0.000000000 0.000000000"
            for s, x in zip(("H", "C", "N"), _X)]
    return f"3\nsynthetic ts A\n" + "\n".join(dots) + "\n"


def _hess_text():
    out = ["$atoms", "3"]
    for s, x in zip(("H", "C", "N"), _X):
        bx = _bohr([x])[0]
        out.append(f"{s} 1.0078 {bx:.9f} 0.000000000 0.000000000")
    out.append("$normal_modes")
    out.append("9 1")
    out.append("0")
    vec = [1.0, 0, 0, -1.0, 0, 0, 0, 0, 0]  # H +x, C -x, N static
    for r in range(9):
        out.append(f"{r} {vec[r]:.8f}")
    out.append("$vibrational_frequencies")
    out.append("1")
    out.append("0 -1500.000")
    return "\n".join(out) + "\n"


def main(argv):
    if len(argv) != 2:
        print("usage: tf1_a_seed.py <workspace>", file=sys.stderr)
        return 2
    d = Path(argv[1]) / "tf1" / "a"
    if (d / "ts.xyz").exists() or (d / "freq.hess").exists():
        print("refusing to overwrite existing tf1/a tree", file=sys.stderr)
        return 1
    d.mkdir(parents=True, exist_ok=True)
    (d / "ts.xyz").write_text(_xyz_text(), encoding="utf-8")
    (d / "freq.hess").write_text(_hess_text(), encoding="utf-8")
    print(f"seeded {d}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
