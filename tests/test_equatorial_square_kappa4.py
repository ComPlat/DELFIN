"""κ4 equatorial-square restriction (#315): a RIGID PLANAR tetradentate (porphyrin /
phthalocyanine) must occupy 4 vertices that are COPLANAR THROUGH THE METAL — the
OC-6 equatorial square (axial sites free) or the full SP-4 square — never a folded
4-subset.  Pure geometry/combinatorics; no embed, deterministic, fast."""
import os
os.environ.setdefault("PYTHONHASHSEED", "0")
import numpy as np
import pytest
from delfin.fffree.polya_isomer_count import _equatorial_squares, enumerate_chelate_configs
from delfin.fffree import polyhedra as PLY


def _coplanar_through_metal(verts, shape):
    V = PLY.ref_vectors(shape)
    M = np.array([V[v] for v in verts], float)
    _, sv, _ = np.linalg.svd(M)            # vectors pass through origin (= metal)
    return float(sv[2] / sv[0])            # ~0 => coplanar through the metal


def test_octahedron_has_exactly_three_equatorial_squares():
    sq = _equatorial_squares("octahedron", 6)
    assert len(sq) == 3, sq
    for s in sq:
        assert len(s) == 4
        assert _coplanar_through_metal(s, "OC-6 octahedron") < 1e-6


def test_square_planar_has_one_square_all_four_vertices():
    sq = _equatorial_squares("square_planar", 4)
    assert sq == [(0, 1, 2, 3)]


@pytest.mark.parametrize("geom,n", [
    ("tetrahedron", 4),               # no square (109.47deg, not 90deg)
    ("trigonal_bipyramid", 5),        # axial pair + trigonal equator, no square
    ("square_pyramid", 5),            # basal square NOT coplanar through metal (domed)
    ("pentagonal_bipyramid", 7),
])
def test_no_equatorial_square_for_non_square_polyhedra(geom, n):
    assert _equatorial_squares(geom, n) == []


def test_porphyrin_plus_two_axial_on_octahedron_is_one_coplanar_config():
    specs = [
        {"type": "porph", "denticity": 4, "asym": False, "rigid_planar": True},
        {"type": "ax", "denticity": 1, "asym": False},
        {"type": "ax", "denticity": 1, "asym": False},
    ]
    cfg = enumerate_chelate_configs("octahedron", specs)
    assert len(cfg) == 1, f"expected single porphyrin orbit, got {len(cfg)}"
    c = cfg[0]
    porph = sorted(v for v, (li, a) in c.items() if li == 0)
    axial = sorted(v for v, (li, a) in c.items() if li != 0)
    assert len(porph) == 4 and len(axial) == 2
    # the 4 macrocycle donors must be coplanar through the metal (the equatorial square)
    assert _coplanar_through_metal(porph, "OC-6 octahedron") < 1e-6
    # the 2 axial donors must be the antipodal axis (mutually ~180deg)
    V = PLY.ref_vectors("OC-6 octahedron")
    assert float(V[axial[0]] @ V[axial[1]]) < -0.99


def test_porphyrin_alone_on_square_planar_is_one_config():
    specs = [{"type": "porph", "denticity": 4, "asym": False, "rigid_planar": True}]
    assert len(enumerate_chelate_configs("square_planar", specs)) == 1


def test_non_rigid_kappa4_uses_generic_branch_not_square():
    # a FLEXIBLE κ4 (cyclam, rigid_planar False) must NOT be square-restricted; it
    # keeps the historic over-enumerate-and-self-gate combinatorial freedom.
    specs = [
        {"type": "cyc", "denticity": 4, "asym": False, "rigid_planar": False},
        {"type": "ax", "denticity": 1, "asym": False},
        {"type": "ax", "denticity": 1, "asym": False},
    ]
    cfg = enumerate_chelate_configs("octahedron", specs)
    # generic dent>=3 branch yields >1 vertex subset (cis + trans folds), unlike the
    # single square-restricted porphyrin config.
    assert len(cfg) >= 2


def test_enumeration_is_deterministic():
    specs = [
        {"type": "porph", "denticity": 4, "asym": False, "rigid_planar": True},
        {"type": "ax", "denticity": 1, "asym": False},
        {"type": "ax", "denticity": 1, "asym": False},
    ]
    a = enumerate_chelate_configs("octahedron", specs)
    b = enumerate_chelate_configs("octahedron", specs)
    assert a == b
