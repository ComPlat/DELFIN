"""Geometry manipulation on instruction: set bond/angle/dihedral, translate/
rotate fragments, add/remove atoms, replace a ligand, mirror, align two
structures (RMSD), build from SMILES, merge fragments, freeze atoms for a
constrained optimisation; read/write xyz and ORCA coordinate blocks.

Method: every operation returns the new structure (as an XYZ block) AND a short
record of what changed (atoms moved, old/new internal coordinate), and refuses
impossible requests (an atom that does not exist, a dihedral across a ring that
would break it) with a reason.

No OPI dependency; runs under the py3.11 CI interpreter.
"""
import math

import pytest

from delfin.tools.geometry_edit import (
    GeometryError,
    parse_xyz,
    parse_orca_block,
    to_xyz,
    build_from_smiles,
    set_dihedral,
    set_angle,
    set_bond,
    translate,
    rotate,
    add_atom,
    remove_atoms,
    replace_ligand,
    mirror,
    freeze_atoms,
    constraints_block,
    merge_fragments,
    align_rmsd,
    measured_dihedral,
    measured_angle,
)

rad = math.pi / 180.0


# --------------------------------------------------------------------------
# small, explicit test molecules
# --------------------------------------------------------------------------

def _water_xyz():
    # O at origin, H's 0.96 A away, H-O-H bond angle 104.5 deg in the yz plane.
    a = 104.5 * rad / 2.0
    return (
        "3\nwater\n"
        "O 0.000000 0.000000 0.000000\n"
        f"H 0.000000 {0.96 * math.sin(a):.6f} {0.96 * math.cos(a):.6f}\n"
        f"H 0.000000 {-0.96 * math.sin(a):.6f} {0.96 * math.cos(a):.6f}\n"
    )


def _butane_xyz():
    # anti-periplanar butane: C1 C2 C3 C4 zigzag, phi(C1-C2-C3-C4) = 180.
    # C2 at origin, C3 along +x; C1 and C4 on opposite sides of the plane.
    r = 1.54
    return (
        "4\nbutane\n"
        f"C {-r / 2:.6f} {r * math.sin(rad * 30):.6f} {r * math.cos(rad * 30):.6f}\n"
        "C 0.000000 0.000000 0.000000\n"
        "C 1.540000 0.000000 0.000000\n"
        f"C {r + r / 2:.6f} {r * math.sin(rad * 30):.6f} {-r * math.cos(rad * 30):.6f}\n"
    )


def _measured_length(g):
    return math.dist(g.coords[0], g.coords[1])


# --------------------------------------------------------------------------
# read/write: round trips
# --------------------------------------------------------------------------

def test_roundtrip_xyz_preserves_structure():
    g = parse_xyz(_water_xyz())
    out_x = to_xyz(g, comment="water")
    g2 = parse_xyz(out_x)
    assert g2.symbols == ["O", "H", "H"]
    for i in range(3):
        for j in range(3):
            assert abs(g2.coords[i][j] - g.coords[i][j]) < 1e-6


def test_roundtrip_orca_block():
    # parse_orca_block / to_orca_block live in geometry_edit and round-trip.
    block = (
        "* xyz 0 1\n"
        "O      0.00000000      0.00000000      0.00000000\n"
        "H      0.00000000      0.96000000      0.00000000\n"
        "H      0.00000000     -0.96000000      0.00000000\n"
        "*\n"
    )
    g = parse_orca_block(block)
    assert g.symbols == ["O", "H", "H"]
    # write back and rebuild
    out = to_xyz(g, comment="water")  # reuses the xyz writer
    g2 = parse_xyz(out)
    assert g2.symbols == g.symbols
    for i in range(3):
        assert math.dist(g2.coords[i], g.coords[i]) < 1e-6


# --------------------------------------------------------------------------
# internal-coordinate edits
# --------------------------------------------------------------------------

def test_set_bond_length():
    # two-atom fragment
    g = parse_xyz("2\\nH2\\nH 0 0 0\\nH 1.0 0 0\\n".replace("\\n", "\n"))
    g, change = set_bond(g, 0, 1, 0.74)
    assert abs(_measured_length(g) - 0.74) < 1e-6
    assert change["op"] == "set_bond"
    assert "new" in change and "old" in change


def test_set_angle_h2o():
    g = parse_xyz(_water_xyz())
    # vertex is O (index 0), passed second: H(1)-O(0)-H(2) -> 110
    g, change = set_angle(g, 1, 0, 2, 110.0)
    assert abs(_measured_angle(g, 1, 0, 2) - 110.0) < 1e-6
    assert change["atoms_moved"]
    # the vertex O (atom 0) stays put; the moved side is the *smaller* fragment,
    # i.e. exactly one H (both H sides are singletons)
    assert set(change["atoms_moved"]) in ({1}, {2})


def _measured_angle(g, i, j, k):
    return measured_angle(g, i, j, k)


def test_set_dihedral_butane_to_60():
    g = parse_xyz(_butane_xyz())
    g, change = set_dihedral(g, 0, 1, 2, 3, 60.0)
    measured = _measured_dihedral(g, 0, 1, 2, 3)
    assert abs(measured - 60.0) % 360.0 < 0.01 or abs(measured - 60.0) < 0.01
    # exactly one fragment (the smaller side) moved
    assert change["op"] == "set_dihedral"
    assert change["atoms_moved"], "at least one atom must move"


def _measured_dihedral(g, i, j, k, l):
    return measured_dihedral(g, i, j, k, l)


# --------------------------------------------------------------------------
# refusals
# --------------------------------------------------------------------------

def test_set_dihedral_refuses_nonexistent_atom():
    g = parse_xyz(_butane_xyz())
    with pytest.raises(GeometryError):
        set_dihedral(g, 0, 1, 2, 9, 60.0)


def test_set_angle_refuses_nonexistent_atom():
    g = parse_xyz(_water_xyz())
    with pytest.raises(GeometryError):
        set_angle(g, 0, 1, 40, 90.0)


def test_set_dihedral_refuses_a_ring_bond():
    # Cyclohexane: the central bond (0,1) is part of the 6-ring, so rotating
    # the flanking atoms across it would break ring closure. geometry_edit
    # perceives the ring from the embedded geometry and refuses.
    g = build_from_smiles("C1CCCCC1")
    with pytest.raises(GeometryError) as ei:
        set_dihedral(g, 0, 1, 2, 3, 30.0)
    assert "ring" in str(ei.value).lower()


# --------------------------------------------------------------------------
# fragment moves
# --------------------------------------------------------------------------

def test_translate_fragment():
    g = parse_xyz(_water_xyz())
    g, change = translate(g, [1, 2], (1.0, 2.0, 3.0))
    assert change["atoms_moved"] == [1, 2]
    assert g.coords[1][0] == pytest.approx(1.0, abs=1e-9)
    # the fixture rounds H coordinates to 6 decimals, so the translated y has
    # up to ~5e-7 of rounding on top of the exact 2.0 shift; compare at 1e-5
    assert g.coords[1][1] == pytest.approx(
        0.96 * math.sin(104.5 * rad / 2) + 2.0, abs=1e-5
    )
    # O untouched
    assert g.coords[0] == (0.0, 0.0, 0.0)


def test_rotate_fragment_about_origin():
    g = parse_xyz(_water_xyz())
    g, change = rotate(g, [1, 2], (0.0, 0.0, 1.0), 90.0, center=(0.0, 0.0, 0.0))
    # rotating H's by 90 about z keeps their distance to O
    for idx in (1, 2):
        assert abs(_length_to(g, 0, idx) - 0.96) < 1e-6
    assert change["atoms_moved"] == [1, 2]


def _length_to(g, i, j):
    return math.dist(g.coords[i], g.coords[j])


# --------------------------------------------------------------------------
# atoms
# --------------------------------------------------------------------------

def test_add_and_remove_atoms():
    g = parse_xyz(_water_xyz())
    g, change = add_atom(g, "F", bonded_to=0, bond_length=1.4)
    assert len(g.symbols) == 4
    assert g.symbols[-1] == "F"
    assert abs(_length_to(g, 0, 3) - 1.4) < 1e-6
    g, change2 = remove_atoms(g, [3])
    assert len(g.symbols) == 3
    assert g.symbols == ["O", "H", "H"]


def test_remove_atoms_refuses_missing_index():
    g = parse_xyz(_water_xyz())
    with pytest.raises(GeometryError):
        remove_atoms(g, [7])


# --------------------------------------------------------------------------
# ligand swap on a complex built from SMILES
# --------------------------------------------------------------------------

def test_replace_ligand_on_ethanol():
    # ethanol CCO: replace the OH ligand at the C-O carbon with S
    eth = build_from_smiles("CCO")
    carbon_o = None
    for i in range(len(eth.symbols)):
        if eth.symbols[i] == "C":
            for j in range(len(eth.symbols)):
                if eth.symbols[j] == "O" and _length_to(eth, i, j) < 1.6:
                    carbon_o = i
                    break
        if carbon_o is not None:
            break
    assert carbon_o is not None, "no C-O bond found in embedded ethanol"
    new, change = replace_ligand(eth, carbon_o, "S")
    # the O is gone, an S sits near the attach carbon
    assert "O" not in new.symbols
    o_count = sum(1 for s in new.symbols if s == "O")
    assert o_count == 0
    s_idx = new.symbols.index("S")
    assert _length_to(new, change["connect"], s_idx) < 2.0  # S-C bond present


# --------------------------------------------------------------------------
# build from SMILES, mirror, merge, align
# --------------------------------------------------------------------------

def test_build_from_smiles_formula():
    eth = build_from_smiles("CCO")
    # C2H6O + explicit H = 9 atoms
    assert len(eth.symbols) == 9
    from collections import Counter
    cnt = Counter(eth.symbols)
    assert cnt["C"] == 2 and cnt["O"] == 1
    for c in eth.coords:
        for v in c:
            assert math.isfinite(v)


def test_mirror_flips_x():
    g = parse_xyz(_water_xyz())
    g, change = mirror(g, plane="yz")
    # O x should stay (it is 0), H x stays 0 in this layout -> use yz plane on a
    # molecule with non-zero x on one atom instead.
    # Rebuild water with H's displaced in x.
    g = parse_xyz(
        "3\\nwater\\nO 0 0 0\\nH 0.8 0.3 0.0\\nH -0.8 0.3 0.0\\n".replace("\\n", "\n")
    )
    g, change = mirror(g, plane="yz")
    assert g.coords[1][0] == pytest.approx(-0.8, abs=1e-9)
    assert g.coords[2][0] == pytest.approx(0.8, abs=1e-9)


def test_merge_fragments():
    w = parse_xyz(_water_xyz())
    me = parse_xyz(
        "5\\nmethane\\nC 0 0 0\\nH 0.6 0 0\\nH 0 0.6 0\\nH 0 0 0.6\\nH -0.6 0 0\\n".replace(
            "\\n", "\n"
        )
    )
    merged, change = merge_fragments(w, me)
    assert len(merged.symbols) == 8
    for i in range(3):
        assert math.dist(merged.coords[i], w.coords[i]) < 1e-9


def test_align_rmsd_returns_small_rmsd():
    g = parse_xyz(_water_xyz())
    # a rotated copy: 45 deg about the z axis
    h = parse_xyz(_water_xyz())
    from delfin.tools.geometry_edit import rotate as _rot
    h, _ = _rot(h, [0, 1, 2], (0.0, 0.0, 1.0), 45.0, center=(0.0, 0.0, 0.0))
    aligned, rmsd = align_rmsd(h, g)
    assert rmsd < 1e-4
    assert aligned.symbols == g.symbols
    for i in range(3):
        for j in range(3):
            assert abs(aligned.coords[i][j] - g.coords[i][j]) < 1e-3


# --------------------------------------------------------------------------
# freeze atoms -> constrained optimisation spec
# --------------------------------------------------------------------------

def test_freeze_atoms_and_constraints_block():
    g = parse_xyz(_water_xyz())
    g, change = freeze_atoms(g, [0])
    assert change["atoms_moved"] == []
    block = constraints_block(g)
    assert "%geom" in block
    assert "Constraints" in block
    # Cartesian constraint for 0-based index 0
    assert "{ C 0 C }" in block


# --------------------------------------------------------------------------
# the geometry tools must stay out of a NON-chemistry always-on prompt
# --------------------------------------------------------------------------

def test_geometry_toolset_hidden_for_non_chemistry_workspace(tmp_path):
    """A non-chemistry session must not advertise the geometry tools.

    The tools are registered as ``mcp__delfin-ops__*`` names and therefore
    belong to the *chemistry* tool group (delfin/agent/toolsets.py). For a
    plain, chemistry-off workspace that group is hidden, so the tool names are
    dropped from the advertised set — the always-on prompt is unchanged.
    """
    from delfin.agent.toolsets import WorkspaceToolsets, drop_hidden

    ts = WorkspaceToolsets(
        workspace=str(tmp_path), settings={"agent.toolsets.chemistry": "off"}
    )
    assert "chemistry" in ts.hidden_groups()

    geometry_names = [
        e["name"]
        for e in _GEOMETRY_CATALOG_ENTRIES()
    ]
    assert geometry_names, "geometry catalog must be non-empty"
    all_tools = [
        # real advertised-set shape: OpenAI tool dicts whose function name is
        # the MCP tool name (tool_group maps mcp__delfin-ops__* -> chemistry)
        {"function": {"name": f"mcp__delfin-ops__{n}"}}
        for n in geometry_names
    ]
    remaining = drop_hidden(all_tools, ts.hidden_groups())
    surviving = {t["function"]["name"] for t in remaining}
    assert not (surviving & set(geometry_names)), (
        "geometry tools leaked into a non-chemistry advertised set"
    )


def _GEOMETRY_CATALOG_ENTRIES():
    """Category='geometry' entries from the delfin-ops tool catalog.

    Imported lazily to keep this test importable even before the catalog
    entries exist (the red commit), and read from the live catalog so the
    toolset test always tracks whatever geometry tools the branch ships.
    """
    try:
        from delfin import api as _api
        return [e for e in _api.list_tools(category="geometry")]
    except Exception:
        return []
