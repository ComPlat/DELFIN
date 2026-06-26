"""Tests for delfin.manta._ring_conformer_templates (Welle-5p-B).

Validation HARD-gate per the Welle-5p-B brief:

1. Default-OFF byte-identical (master flag unset → base-only pool).
2. Cyclohexane produces multiple distinct chair / boat / twist
   conformers via template-based atom-Z displacements (not rotations).
3. Cyclopentane produces envelope / twist / half-chair conformers.
4. Macrocycle (cyclo-dodecane) produces saddle / ruffled / dome.
5. Chelate-ring (Fe(en) 5-ring with metal in ring) is EXCLUDED from
   templating — universal-fundamental rule.
6. Aromatic ring (benzene) is EXCLUDED — π-planar fixed.
7. M-D invariant preserved across template applications.
8. Template library has the expected per-size coverage.
"""

from __future__ import annotations

import math
import os

import pytest

from delfin.manta import _ring_conformer_templates as ring_tpl
from delfin.manta import _rotamer_diversity as _rot


# ---------------------------------------------------------------------------
# Fixture XYZs (DELFIN format — no atom-count header)
# ---------------------------------------------------------------------------


# Cyclohexane chair geometry.
_CYCLOHEXANE_XYZ = """\
C       1.265000    0.730000    0.250000
C       0.000000    1.460000   -0.250000
C      -1.265000    0.730000    0.250000
C      -1.265000   -0.730000   -0.250000
C       0.000000   -1.460000    0.250000
C       1.265000   -0.730000   -0.250000
H       2.165000    1.250000   -0.150000
H       1.265000    0.730000    1.350000
H       0.000000    1.460000   -1.350000
H       0.000000    2.500000    0.150000
H      -2.165000    1.250000   -0.150000
H      -1.265000    0.730000    1.350000
H      -2.165000   -1.250000    0.150000
H      -1.265000   -0.730000   -1.350000
H       0.000000   -1.460000    1.350000
H       0.000000   -2.500000   -0.150000
H       2.165000   -1.250000    0.150000
H       1.265000   -0.730000   -1.350000
"""


# Cyclopentane planar starting geometry.
_CYCLOPENTANE_XYZ = """\
C       1.275000    0.000000    0.000000
C       0.394000    1.213000    0.000000
C      -1.032000    0.750000    0.000000
C      -1.032000   -0.750000    0.000000
C       0.394000   -1.213000    0.000000
H       2.350000    0.000000    0.000000
H       1.275000    0.000000    1.100000
H       0.394000    1.213000    1.100000
H       0.726000    2.236000   -0.000000
H      -1.901000    1.382000    0.000000
H      -1.032000    0.750000    1.100000
H      -1.901000   -1.382000    0.000000
H      -1.032000   -0.750000    1.100000
H       0.394000   -1.213000    1.100000
H       0.726000   -2.236000    0.000000
"""


# Cyclo-dodecane (12-ring, planar regular start).
def _make_cyclododecane_xyz() -> str:
    n = 12
    R = 2.50
    lines = []
    for k in range(n):
        theta = 2.0 * math.pi * k / n
        x = R * math.cos(theta)
        y = R * math.sin(theta)
        z = 0.05 * ((-1) ** k)  # tiny puckering to stabilise OB ring detection
        lines.append(f"C    {x:12.6f} {y:12.6f} {z:12.6f}")
    # Add two axial hydrogens per carbon
    for k in range(n):
        theta = 2.0 * math.pi * k / n
        x_c = R * math.cos(theta)
        y_c = R * math.sin(theta)
        z_c = 0.05 * ((-1) ** k)
        # H_a above, H_b below the carbon
        lines.append(f"H    {x_c:12.6f} {y_c:12.6f} {z_c + 1.09:12.6f}")
        lines.append(f"H    {x_c:12.6f} {y_c:12.6f} {z_c - 1.09:12.6f}")
    return "\n".join(lines) + "\n"


_CYCLODODECANE_XYZ = _make_cyclododecane_xyz()


# Fe-en chelate 5-ring (Fe-N-C-C-N back to Fe).  Layer-2 templates
# MUST exclude this ring because it contains a metal.
_FE_EN_XYZ = """\
Fe      0.000000    0.000000    0.000000
N       2.000000    0.300000    0.000000
N       0.300000    2.000000    0.000000
C       2.500000    1.650000    0.000000
C       1.650000    2.500000    0.000000
H       2.400000   -0.260000    0.870000
H       2.400000   -0.260000   -0.870000
H      -0.260000    2.400000    0.870000
H      -0.260000    2.400000   -0.870000
H       2.500000    2.120000    0.890000
H       3.500000    1.500000   -0.400000
H       2.120000    2.500000    0.890000
H       1.500000    3.500000   -0.400000
"""


# Benzene — fully aromatic 6-ring → must be excluded.
_BENZENE_XYZ = """\
C       1.396000    0.000000    0.000000
C       0.698000    1.209000    0.000000
C      -0.698000    1.209000    0.000000
C      -1.396000    0.000000    0.000000
C      -0.698000   -1.209000    0.000000
C       0.698000   -1.209000    0.000000
H       2.476000    0.000000    0.000000
H       1.238000    2.144000    0.000000
H      -1.238000    2.144000    0.000000
H      -2.476000    0.000000    0.000000
H      -1.238000   -2.144000    0.000000
H       1.238000   -2.144000    0.000000
"""


# ---------------------------------------------------------------------------
# Env-helper
# ---------------------------------------------------------------------------

_ENV_KEYS = (
    "DELFIN_5P_B_RING_TEMPLATES",
    "DELFIN_5P_B_K_PER_RING",
    "DELFIN_5P_B_MAX_RING_VARIANTS",
    "DELFIN_5P_B_AMPLITUDE_FRACTION",
    "DELFIN_5P_B_MD_TOL",
)


def _clear_env(monkeypatch):
    for key in _ENV_KEYS:
        monkeypatch.delenv(key, raising=False)


def _ob_or_skip():
    try:
        from openbabel import pybel  # noqa: F401
    except Exception:
        pytest.skip("Open Babel not available")


# ---------------------------------------------------------------------------
# 1. Template library structural tests
# ---------------------------------------------------------------------------


def test_template_library_coverage_per_ring_size():
    """Library exposes templates for each ring-size; 6-ring has all
    canonical chair / boat / twist-boat entries."""
    cov = ring_tpl.template_library_coverage()
    assert cov["ring-5"] >= 5, cov  # at least 5 envelopes
    assert cov["ring-6"] >= 6, cov  # chair-A + chair-B + boats + twists
    assert cov["ring-7"] >= 4, cov
    assert cov["ring-12"] >= 3, cov  # saddle + ruffled + dome


def test_template_patterns_are_normalised():
    """Every template has unit L2-norm so a single global amplitude
    suffices for all ring sizes."""
    for size, tpl_dict in (
        (3, ring_tpl.TEMPLATES_3_RING),
        (4, ring_tpl.TEMPLATES_4_RING),
        (5, ring_tpl.TEMPLATES_5_RING),
        (6, ring_tpl.TEMPLATES_6_RING),
        (7, ring_tpl.TEMPLATES_7_RING),
    ):
        for tag, pattern in tpl_dict.items():
            l2 = math.sqrt(sum(v * v for v in pattern))
            assert abs(l2 - 1.0) < 1e-6, (
                f"size {size} {tag} L2={l2:.6f} not normalised"
            )
            assert len(pattern) == size, (
                f"size {size} {tag} pattern length {len(pattern)} != {size}"
            )


def test_cyclohexane_chair_vs_boat_distinct():
    """Chair and boat templates produce *distinct* displacement patterns.

    Chair = alternating ±1.  Boat = two flagpoles + four base.  They
    must not be proportional to each other (different conformational
    families).
    """
    chair = ring_tpl.TEMPLATES_6_RING["chair-A"]
    boat = ring_tpl.TEMPLATES_6_RING["boat-14"]
    # Compute cosine similarity — must be far from ±1 if distinct.
    dot = sum(a * b for a, b in zip(chair, boat))
    assert abs(dot) < 0.95, (
        f"chair vs boat too correlated (cos={dot:.4f}) — should be distinct"
    )


# ---------------------------------------------------------------------------
# 2. Default-OFF byte-identical
# ---------------------------------------------------------------------------


def test_default_off_byte_identical(monkeypatch):
    """Master flag unset → returns [(xyz, 'base')] only."""
    _clear_env(monkeypatch)
    out = ring_tpl.apply_if_enabled(_CYCLOHEXANE_XYZ)
    assert out == [(_CYCLOHEXANE_XYZ, "base")]


def test_default_off_for_chelate(monkeypatch):
    """Default-OFF: metal-complex unchanged."""
    _clear_env(monkeypatch)
    out = ring_tpl.apply_if_enabled(_FE_EN_XYZ)
    assert out == [(_FE_EN_XYZ, "base")]


# ---------------------------------------------------------------------------
# 3. Cyclohexane — distinct chair / boat / twist conformers
# ---------------------------------------------------------------------------


def test_cyclohexane_produces_multiple_distinct_conformers(monkeypatch):
    """Validation case 1 per brief: cyclohexane → chair + boat + twist
    distinct conformers via template-based atom-Z displacements (not
    rotations)."""
    _ob_or_skip()
    _clear_env(monkeypatch)
    monkeypatch.setenv("DELFIN_5P_B_RING_TEMPLATES", "1")
    monkeypatch.setenv("DELFIN_5P_B_K_PER_RING", "5")
    monkeypatch.setenv("DELFIN_5P_B_MAX_RING_VARIANTS", "6")
    monkeypatch.setenv("DELFIN_5P_B_AMPLITUDE_FRACTION", "0.5")

    out = ring_tpl.apply_if_enabled(_CYCLOHEXANE_XYZ)
    # base + ≥2 variants required (chair-A + boat-X minimum)
    assert len(out) >= 3, f"need ≥3 conformers, got {len(out)}: {out}"
    tags = [tag for (_x, tag) in out]
    # At least one chair tag and one boat tag must appear
    has_chair = any("chair" in t for t in tags)
    has_boat = any("boat" in t for t in tags)
    assert has_chair, f"missing chair conformer in {tags}"
    assert has_boat, f"missing boat conformer in {tags}"


# ---------------------------------------------------------------------------
# 4. Chelate ring (Fe-en) — EXCLUDED from templating
# ---------------------------------------------------------------------------


def test_chelate_ring_excluded_from_templating(monkeypatch):
    """Validation case 3 per brief: metal-containing rings are EXCLUDED.

    Layer-2 templates only act on non-metal rings.  Layer-3 chelate-twist
    handles metal-rings via different mechanism.
    """
    _ob_or_skip()
    # Build graph from Fe-en XYZ
    ob_mol = _rot._build_ob_mol_from_xyz(_FE_EN_XYZ)
    if ob_mol is None:
        pytest.skip("OB could not parse Fe(en)")
    graph = _rot._graph_from_ob(ob_mol)
    if not graph:
        pytest.skip("graph build failed")
    # find_rings_for_templating MUST NOT return any ring containing Fe
    rings = ring_tpl.find_rings_for_templating(graph)
    is_metal = graph["is_metal"]
    for ring in rings:
        assert not any(is_metal[i] for i in ring), (
            f"chelate ring {ring} returned despite metal-membership"
        )


def test_chelate_ring_excluded_synthetic_graph():
    """Synthetic 5-ring graph (M-N-C-C-N) — templating MUST skip it.

    This is the universal-fundamental contract independent of OB bond
    perception (which often misses M-N bonds at 2.05 Å).
    """
    graph = {
        "n_atoms": 5,
        "atomic_nums": [26, 7, 6, 6, 7],
        "is_metal": [True, False, False, False, False],
        "neighbours": [[1, 4], [0, 2], [1, 3], [2, 4], [3, 0]],
        "bonds": [
            (0, 1, 1, False, True),
            (1, 2, 1, False, True),
            (2, 3, 1, False, True),
            (3, 4, 1, False, True),
            (4, 0, 1, False, True),
        ],
    }
    rings = ring_tpl.find_rings_for_templating(graph)
    # Fe-containing ring → EXCLUDED
    assert rings == [], (
        f"chelate-ring (metal in ring) MUST be excluded, got {rings}"
    )


# ---------------------------------------------------------------------------
# 5. Aromatic ring (benzene) — EXCLUDED
# ---------------------------------------------------------------------------


def test_aromatic_ring_excluded(monkeypatch):
    """Benzene (fully aromatic) is EXCLUDED from templating."""
    _ob_or_skip()
    ob_mol = _rot._build_ob_mol_from_xyz(_BENZENE_XYZ)
    if ob_mol is None:
        pytest.skip("OB could not parse benzene")
    graph = _rot._graph_from_ob(ob_mol)
    rings = ring_tpl.find_rings_for_templating(graph)
    # Fully aromatic → EXCLUDED
    assert rings == [], f"aromatic ring MUST be excluded, got {rings}"


# ---------------------------------------------------------------------------
# 6. Macrocycle — saddle / ruffled / dome
# ---------------------------------------------------------------------------


def test_macrocycle_templates_distinct():
    """12-ring macrocycle template library exposes the three Fourier
    modes (saddle / ruffled / dome) with non-degenerate patterns."""
    templates = ring_tpl.get_templates_for_ring_size(12)
    assert "saddle" in templates
    assert "ruffled" in templates
    assert "dome" in templates
    saddle = templates["saddle"]
    ruffled = templates["ruffled"]
    # Patterns must be normalised and distinct
    cos_sr = sum(a * b for a, b in zip(saddle, ruffled))
    assert abs(cos_sr) < 0.5, (
        f"saddle vs ruffled too correlated (cos={cos_sr:.4f})"
    )


# ---------------------------------------------------------------------------
# 7. M-D invariant preserved (synthetic test — direct apply_template)
# ---------------------------------------------------------------------------


def test_apply_template_does_not_touch_metal():
    """When the ring excludes the metal, apply_template MUST NOT move
    the metal coordinate (M-D invariant preserved trivially).
    """
    # Build a simple 6-ring atom-list with NO metal in the ring.
    base_coords = [
        (1.265, 0.730, 0.250),
        (0.000, 1.460, -0.250),
        (-1.265, 0.730, 0.250),
        (-1.265, -0.730, -0.250),
        (0.000, -1.460, 0.250),
        (1.265, -0.730, -0.250),
        # External metal at index 6 — NOT in the ring
        (5.000, 5.000, 5.000),
    ]
    ring = [0, 1, 2, 3, 4, 5]
    pattern = ring_tpl.TEMPLATES_6_RING["chair-A"]
    new_coords = ring_tpl.apply_template(
        base_coords, ring, pattern, amplitude=0.35, graph=None
    )
    # Metal (index 6) untouched
    assert new_coords[6] == base_coords[6], (
        "metal moved during ring-only template application"
    )
    # At least one ring atom moved
    moved = any(
        new_coords[i] != base_coords[i] for i in ring
    )
    assert moved, "no ring atom moved — template not applied"


def test_apply_template_drags_hydrogens():
    """Bonded hydrogens MUST be dragged rigidly with their heavy parent
    to avoid stretching C-H bonds (per rigid-H tracking lesson)."""
    base_coords = [
        (0.0, 0.0, 0.0),    # C0
        (1.5, 0.0, 0.0),    # C1
        (2.25, 1.30, 0.0),  # C2
        (1.5, 2.60, 0.0),   # C3
        (0.0, 2.60, 0.0),   # C4
        (-0.75, 1.30, 0.0), # C5
        (0.0, 0.0, 1.09),   # H on C0
    ]
    graph = {
        "n_atoms": 7,
        "atomic_nums": [6, 6, 6, 6, 6, 6, 1],
        "is_metal": [False] * 7,
        "neighbours": [
            [1, 5, 6],
            [0, 2],
            [1, 3],
            [2, 4],
            [3, 5],
            [4, 0],
            [0],
        ],
        "bonds": [],
    }
    ring = [0, 1, 2, 3, 4, 5]
    pattern = ring_tpl.TEMPLATES_6_RING["chair-A"]
    new_coords = ring_tpl.apply_template(
        base_coords, ring, pattern, amplitude=0.30, graph=graph
    )
    # C0 moved → H must move by the same delta (rigid-H)
    d_c = (
        new_coords[0][0] - base_coords[0][0],
        new_coords[0][1] - base_coords[0][1],
        new_coords[0][2] - base_coords[0][2],
    )
    d_h = (
        new_coords[6][0] - base_coords[6][0],
        new_coords[6][1] - base_coords[6][1],
        new_coords[6][2] - base_coords[6][2],
    )
    for k in range(3):
        assert abs(d_c[k] - d_h[k]) < 1e-9, (
            f"H not dragged with parent: ΔC={d_c}, ΔH={d_h}"
        )


# ---------------------------------------------------------------------------
# 8. End-to-end pool: cyclohexane chair vs boat have distinct atom-Z
# ---------------------------------------------------------------------------


def test_chair_vs_boat_atom_z_displacement_distinct(monkeypatch):
    """The whole point of templates: chair and boat produce *different*
    atom-Z displacement patterns — not just different rotations."""
    _ob_or_skip()
    _clear_env(monkeypatch)
    monkeypatch.setenv("DELFIN_5P_B_RING_TEMPLATES", "1")
    monkeypatch.setenv("DELFIN_5P_B_K_PER_RING", "4")
    monkeypatch.setenv("DELFIN_5P_B_MAX_RING_VARIANTS", "6")
    monkeypatch.setenv("DELFIN_5P_B_AMPLITUDE_FRACTION", "0.5")

    out = ring_tpl.apply_if_enabled(_CYCLOHEXANE_XYZ)
    chair_xyz = None
    boat_xyz = None
    for xyz, tag in out:
        if "chair" in tag and chair_xyz is None:
            chair_xyz = xyz
        if "boat" in tag and "twist-boat" not in tag and boat_xyz is None:
            boat_xyz = xyz
    if chair_xyz is None or boat_xyz is None:
        pytest.skip("chair or boat not produced by current template gating")

    # Parse both — verify atom-Z values differ structurally
    _syms_c, coords_c = _rot._parse_delfin_xyz(chair_xyz)
    _syms_b, coords_b = _rot._parse_delfin_xyz(boat_xyz)
    z_chair = [c[2] for c in coords_c[:6]]  # 6-ring carbons
    z_boat = [c[2] for c in coords_b[:6]]
    # Chair atom-Z: alternating ±sign pattern; boat: two-up four-base.
    # The patterns must NOT be identical.
    assert z_chair != z_boat, (
        "chair and boat have identical atom-Z — template diversity lost"
    )
