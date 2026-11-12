"""Round 9 (QS-XE1-F1): a half-aromatic ring carbon takes ONE hydrogen.

QS finding on the round-8 head (b7fa64bab): in a distorted benzene ring
where two sides are stretched to single-bond length, the carbons at the
stretched bonds carry a fractional heavy bond order (1.5 aromatic +
1.0 single = 2.5) and the C at the pushed vertex carries two plain
singles (2.0). int(round(missing)) then gave TWO hydrogens on ring
carbons, making them sp3 and destroying the pi sextet. The fractional
part of an aromatic valence sum means "at least one", never "round up
to a full CH2" for a carbon inside an aromatic ring context.

Control both ways: regular benzene stays 6 H and a plain saturated ring
(cyclohexane, all bonds single) still gets its 2 H per carbon -- the
clamp must not fire outside an aromatic ring context.
"""

from delfin.structure_edit import fix_hydrogens
from collections import Counter
import math


def _benzene_distorted() -> str:
    """Regular hexagon side 1.42 A, vertex 3 pushed 0.16 A outward.

    Sides 2-3 and 3-4 become 1.506 A (ratio 0.99 -> single order 1.0);
    the other four sides stay aromatic (1.5). Atom 3 sits between two
    single-looking bonds, atoms 2 and 4 between one aromatic and one
    single bond -- exactly the QS case.
    """
    pts = []
    for k in range(6):
        ang = math.radians(60 * k)
        r = 1.42 + (0.16 if k == 3 else 0.0)
        pts.append((r * math.cos(ang), r * math.sin(ang), 0.0))
    lines = ["C6 distorted-benzene QS-XE1-F1"]
    for x, y, z in pts:
        lines.append(f"C {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines) + "\n"


def _cyclohexane_planar() -> str:
    """Saturated ring, side 1.54 A: every carbon keeps 2 H (no clamp)."""
    pts = []
    for k in range(6):
        ang = math.radians(60 * k)
        pts.append((1.54 * math.cos(ang), 1.54 * math.sin(ang), 0.0))
    lines = ["C6 cyclohexane skeleton"]
    for x, y, z in pts:
        lines.append(f"C {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines) + "\n"


def _count_h(xyz_out: str):
    """Count element symbols in the fix_hydrogens OUTPUT xyz text."""
    rows = [ln.split() for ln in xyz_out.strip().splitlines()[2:]]
    return Counter(r[0] for r in rows)


def test_distorted_benzene_gets_six_hydrogens():
    out = fix_hydrogens(_benzene_distorted())
    assert _count_h(out["xyz"]) == Counter({"C": 6, "H": 6})


def test_regular_benzene_unchanged_six():
    pts = []
    for k in range(6):
        ang = math.radians(60 * k)
        pts.append((1.39 * math.cos(ang), 1.39 * math.sin(ang), 0.0))
    text = "C6 benzene\n" + "".join(
        f"C {x:.6f} {y:.6f} 0.000000\n" for x, y, _ in pts
    )
    out = fix_hydrogens(text)
    assert _count_h(out["xyz"]) == Counter({"C": 6, "H": 6})


def test_cyclohexane_keeps_two_h_per_carbon():
    out = fix_hydrogens(_cyclohexane_planar())
    assert _count_h(out["xyz"]) == Counter({"C": 6, "H": 12})
