"""XE1: fix_hydrogens on real heavy-atom skeletons (operator real-run FAIL).

The operator ran geometry_fix_hydrogens on real ORCA-optimised geometries
with all hydrogens removed and got chemically wrong results:

1. Methanol skeleton (C, O, C-O 1.42 A): two hydrogens were placed on O at
   the same coordinates (coincident atoms) instead of one O-H at ~0.96 A.
2. HCN skeleton (C, N, C-N 1.16 A): two hydrogens were placed on C because
   the C-N distance was counted as a double bond, instead of one C-H for a
   C#N triple bond.

Controls here (heavy-atom-only xyz, real bond lengths):
- bond order implied by the geometry (1.16 A C-N = triple, 1.42 A C-O =
  single) fixes the H count to the molecule's formula;
- no two hydrogens coincide: the minimum H-H distance must stay above a
  real X-H bond length floor (0.8 A).
"""

from delfin.structure_edit import fix_hydrogens, _parse, _dist


def _skeleton(symbols, coords):
    lines = [f"{len(symbols)}", "operator real-run skeleton"]
    for s, c in zip(symbols, coords):
        lines.append(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}")
    return "\n".join(lines) + "\n"


def _h_bonded_to(out_xyz, heavy_index, max_len=1.35):
    """Indices of hydrogens within max_len of heavy atom heavy_index."""
    g = _parse(out_xyz, "fixed")
    return [
        j
        for j, sym in enumerate(g.symbols)
        if sym == "H"
        and _dist(g.coords[j], g.coords[heavy_index]) <= max_len
    ]


def _min_hh_distance(out_xyz):
    g = _parse(out_xyz, "fixed")
    hs = [j for j, s in enumerate(g.symbols) if s == "H"]
    if len(hs) < 2:
        return None
    return min(
        _dist(g.coords[a], g.coords[b])
        for x, a in enumerate(hs)
        for b in hs[x + 1:]
    )


METHANOL = _skeleton(
    ["C", "O"],
    [(0.0, 0.0, 0.0), (1.42, 0.0, 0.0)],
)
HCN = _skeleton(
    ["C", "N"],
    [(0.0, 0.0, 0.0), (0.0, 0.0, 1.16)],
)


def test_methanol_skeleton_gets_ch3oh():
    out = fix_hydrogens(METHANOL)
    g = _parse(out["xyz"], "fixed")
    assert sorted(g.symbols) == ["C", "H", "H", "H", "H", "O"], sorted(g.symbols)
    on_c = _h_bonded_to(out["xyz"], 0)
    on_o = _h_bonded_to(out["xyz"], 1)
    assert len(on_o) == 1, f"O carries {len(on_o)} H, expected 1"
    assert len(on_c) == 3, f"C carries {len(on_c)} H, expected 3"


def test_methanol_no_coincident_hydrogens():
    out = fix_hydrogens(METHANOL)
    hh = _min_hh_distance(out["xyz"])
    assert hh is not None and hh >= 0.8, f"min H-H distance {hh} < 0.8 A"


def test_hcn_triple_bond_gets_one_h_on_c():
    out = fix_hydrogens(HCN)
    g = _parse(out["xyz"], "fixed")
    assert sorted(g.symbols) == ["C", "H", "N"], sorted(g.symbols)
    on_c = _h_bonded_to(out["xyz"], 0)
    on_n = _h_bonded_to(out["xyz"], 1)
    assert len(on_c) == 1, f"C carries {len(on_c)} H, expected 1 (C#N)"
    assert len(on_n) == 0, f"N carries {len(on_n)} H, expected 0 (triple)"


def test_hcn_no_coincident_hydrogens():
    out = fix_hydrogens(HCN)
    hh = _min_hh_distance(out["xyz"])
    assert hh is None or hh >= 0.8, f"min H-H distance {hh} < 0.8 A"
