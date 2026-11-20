"""Tests for delfin.manta._rotamer_diversity (Welle-5l Track-6)."""

from __future__ import annotations

import os

import pytest

from delfin.manta import _rotamer_diversity as rot


# A simple butane-like XYZ (no metals). Bond C2-C3 is a rotatable DOF.
_BUTANE_XYZ = """\
C       0.000000    0.000000    0.000000
C       1.540000    0.000000    0.000000
C       2.058000    1.452000    0.000000
C       3.598000    1.452000    0.000000
H      -0.366667    1.030400   -0.000000
H      -0.366667   -0.515200    0.892348
H      -0.366667   -0.515200   -0.892348
H       1.906667   -0.515200    0.892348
H       1.906667   -0.515200   -0.892348
H       1.691333    1.967200    0.892348
H       1.691333    1.967200   -0.892348
H       3.964667    0.937000    0.892348
H       3.964667    0.937000   -0.892348
H       3.964667    2.482400    0.000000
"""


# A toy "Co-CH3" coordination geometry, mimicking a sigma C-donor on Co.
# This exercises the metal-bond exclusion path.
_CO_METHYL_XYZ = """\
Co      0.000000    0.000000    0.000000
C       2.000000    0.000000    0.000000
H       2.380000    1.030400    0.000000
H       2.380000   -0.515200    0.892348
H       2.380000   -0.515200   -0.892348
"""


def _clear_env(monkeypatch):
    for key in (
        "DELFIN_5L_T6_ROTAMER_DIVERSITY",
        "DELFIN_5L_T6_ROTAMER_K",
        "DELFIN_5L_T6_ROTAMER_STATES",
        "DELFIN_5L_T6_ROTAMER_MAX_DOFS",
        "DELFIN_5L_T6_ROTAMER_GRID_CAP",
        "DELFIN_5L_T6_ROTAMER_MD_TOL",
    ):
        monkeypatch.delenv(key, raising=False)


def test_default_off_passthrough(monkeypatch):
    """Master flag unset → apply_if_enabled returns [xyz] byte-identical."""
    _clear_env(monkeypatch)
    out = rot.apply_if_enabled(_BUTANE_XYZ)
    assert out == [_BUTANE_XYZ]


def test_parse_and_format_roundtrip():
    symbols, coords = rot._parse_delfin_xyz(_BUTANE_XYZ)
    assert len(symbols) == 14
    assert symbols[0] == "C"
    assert pytest.approx(coords[1][0], abs=1e-9) == 1.54
    formatted = rot._format_delfin_xyz(symbols, coords)
    # round-trip stable on second parse
    symbols2, coords2 = rot._parse_delfin_xyz(formatted)
    assert symbols == symbols2
    for (a, b, c), (d, e, f) in zip(coords, coords2):
        assert pytest.approx(a, abs=1e-6) == d
        assert pytest.approx(b, abs=1e-6) == e
        assert pytest.approx(c, abs=1e-6) == f


def test_butane_has_rotatable_dofs():
    """Butane has 3 rotatable single bonds (C1-C2, C2-C3, C3-C4).

    Only the inner C-C bonds qualify because the terminal C-C bonds have an
    anchor with no heavy neighbour (anchor's other neighbours are only Hs).
    The CH3-rotors are caught by the methyl-rotor rule.
    """
    ob = rot._build_ob_mol_from_xyz(_BUTANE_XYZ)
    if ob is None:
        pytest.skip("Open Babel not installed")
    graph = rot._graph_from_ob(ob)
    dofs = rot.identify_rotamer_dofs(graph, max_dofs=8)
    # Expect at least 1 DOF (the central C2-C3 with heavy anchors)
    assert len(dofs) >= 1, "expected at least 1 rotatable DOF in butane"
    # The most bulky DOF should rotate at least 1 heavy atom + several Hs
    top = dofs[0]
    assert top["total_moved"] >= 1
    assert not top.get("is_methyl", False) or top["score"] > 0


def test_metal_bond_excluded():
    """The Co-C bond must never be flagged as a rotamer DOF."""
    ob = rot._build_ob_mol_from_xyz(_CO_METHYL_XYZ)
    if ob is None:
        pytest.skip("OB cannot parse the Co XYZ")
    graph = rot._graph_from_ob(ob)
    dofs = rot.identify_rotamer_dofs(graph, max_dofs=8)
    for dof in dofs:
        a = dof["anchor"]
        p = dof["pivot"]
        assert not graph["is_metal"][a]
        assert not graph["is_metal"][p]


# An Fe(NH2-CH2-CH2-S) chelate: 5-membered ring Fe-N-C-C-S, every
# backbone bond single, non-aromatic, ring-flagged, and without a
# metal endpoint.  The 2026-05-18 user finding: T6 rotamer diversity
# rotated this chelate backbone.
_FE_NCCS_CHELATE_XYZ = """\
Fe      0.000000    0.000000    0.000000
N       2.050000    0.000000    0.000000
C       2.770000    1.270000    0.000000
C       1.950000    2.500000    0.000000
S       0.250000    2.150000    0.000000
H       2.420000   -0.510000    0.810000
H       2.420000   -0.510000   -0.810000
H       3.420000    1.330000    0.880000
H       3.420000    1.330000   -0.880000
H       2.270000    3.080000    0.880000
H       2.270000    3.080000   -0.880000
"""


def test_ring_bonds_excluded_synthetic():
    """The ring rule alone excludes ring bonds from rotamer DOFs.

    Control pair: a 6-carbon ring whose bonds are ring-flagged yields
    zero DOFs; the identical topology with the ring flag off yields
    all 6 backbone bonds as DOFs.  Neither graph contains a metal, so
    the exclusion is attributable to the ring rule, not the metal rule.
    """
    n = 6

    def ring_graph(is_ring):
        return {
            "n_atoms": n,
            "atomic_nums": [6] * n,
            "is_metal": [False] * n,
            "neighbours": [[(i - 1) % n, (i + 1) % n] for i in range(n)],
            "bonds": [tuple(sorted((i, (i + 1) % n))) + (1, False, is_ring)
                      for i in range(n)],
        }

    assert rot.identify_rotamer_dofs(ring_graph(True), max_dofs=8) == []
    dofs = rot.identify_rotamer_dofs(ring_graph(False), max_dofs=8)
    assert len(dofs) == 6, f"expected all 6 backbone bonds, got {len(dofs)}"


def test_chelate_ring_backbone_not_rotatable():
    """The Fe(NH2-CH2-CH2-S) chelate ring yields zero rotamer DOFs.

    The ring rule in ``identify_rotamer_dofs`` is the unconditional
    successor of the Welle-5p-C hotfix that guarded this user finding;
    the hotfix's env flags and helper API are gone.  This test pins
    the finding motif end to end: the N-C, C-C and C-S backbone bonds
    are single, non-aromatic ring bonds and none has a metal endpoint,
    so the ring rule is the mechanism that excludes them (the metal
    rule would not fire on them).
    """
    ob = rot._build_ob_mol_from_xyz(_FE_NCCS_CHELATE_XYZ)
    if ob is None:
        pytest.skip("Open Babel not installed")
    graph = rot._graph_from_ob(ob)
    backbone = [
        (a, b) for a, b, order, aromatic, ring in graph["bonds"]
        if order == 1 and not aromatic and ring
        and not graph["is_metal"][a] and not graph["is_metal"][b]
    ]
    assert len(backbone) >= 3, (
        f"OB did not perceive the chelate backbone as ring bonds: {backbone}"
    )
    dofs = rot.identify_rotamer_dofs(graph, max_dofs=8)
    assert dofs == [], f"chelate ring backbone must not rotate, got {dofs}"


def test_rodrigues_360_is_identity():
    """Rotating by 2π must return the same coordinates."""
    coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.5, 1.0, 0.0)]
    rotated = rot._rodrigues_rotate(
        coords,
        axis_origin=(0.0, 0.0, 0.0),
        axis_dir=(0.0, 0.0, 1.0),
        angle_rad=2 * 3.141592653589793,
        atom_indices=[2],
    )
    assert pytest.approx(rotated[2][0], abs=1e-9) == 1.5
    assert pytest.approx(rotated[2][1], abs=1e-9) == 1.0
    assert pytest.approx(rotated[2][2], abs=1e-9) == 0.0


def test_grid_iter_cap_respected():
    """Combinatorial grid yields at most *cap* combos."""
    combos = list(rot._grid_iter(n_dofs=5, n_states=3, cap=10))
    assert len(combos) <= 10
    # First combo is identity
    assert combos[0] == (0, 0, 0, 0, 0)


def test_md_invariant_detects_breakage():
    """Moving a metal away should trigger M-D invariant failure."""
    symbols, coords = rot._parse_delfin_xyz(_CO_METHYL_XYZ)
    # OB needs to be available to test this end-to-end
    ob = rot._build_ob_mol_from_xyz(_CO_METHYL_XYZ)
    if ob is None:
        pytest.skip("OB unavailable")
    graph = rot._graph_from_ob(ob)
    # New coords: shift the C atom by +1Å (breaks Co-C distance)
    bad_coords = list(coords)
    bad_coords[1] = (3.0, 0.0, 0.0)  # was 2.0
    assert not rot._coord_bond_invariant_holds(
        graph, coords, bad_coords, tol=0.05
    )
    # Identical coords pass
    assert rot._coord_bond_invariant_holds(
        graph, coords, coords, tol=0.05
    )


def test_apply_returns_base_when_no_dofs():
    """A 1-atom XYZ has no DOFs; apply returns [xyz] verbatim."""
    xyz = "C       0.000000    0.000000    0.000000\n"
    out = rot.apply(xyz, n_per_isomer=3)
    assert out == [xyz]


def test_apply_if_enabled_obeys_flag(monkeypatch):
    _clear_env(monkeypatch)
    monkeypatch.setenv("DELFIN_5L_T6_ROTAMER_DIVERSITY", "1")
    monkeypatch.setenv("DELFIN_5L_T6_ROTAMER_K", "2")
    monkeypatch.setenv("DELFIN_5L_T6_ROTAMER_STATES", "3")
    monkeypatch.setenv("DELFIN_5L_T6_ROTAMER_MAX_DOFS", "3")
    out = rot.apply_if_enabled(_BUTANE_XYZ)
    # Always at least the base frame
    assert len(out) >= 1
    assert out[0] == _BUTANE_XYZ
    # When OB is available + DOFs found, we expect ≥1 extra
    ob = rot._build_ob_mol_from_xyz(_BUTANE_XYZ)
    if ob is not None:
        graph = rot._graph_from_ob(ob)
        dofs = rot.identify_rotamer_dofs(graph)
        if dofs:
            assert len(out) >= 2, "expected at least one rotamer when DOFs exist"
