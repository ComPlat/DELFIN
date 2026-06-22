"""Tests for the universal FF-free conformer-completeness pass
(delfin.fffree.conformer_complete).

Pure-Python / geometry-only assertions that do not require a full build:
  * identity-when-off (byte-identical default semantics),
  * DOF detection rotates the heavy-bearing subtree (bulky pendant included),
  * the topology gate rejects M-D break / spurious bond / H-in-metal-sphere,
  * Kabsch RMSD dedup,
  * ENSEMBLE-global cap + dedup,
  * determinism.
"""
from __future__ import annotations

import os

import pytest

from delfin.fffree import conformer_complete as CC


# --- helpers --------------------------------------------------------------

# A planar square Pt(NH3)2Cl2-like skeleton, headerless DELFIN XYZ.  Pt at
# origin, two N and two Cl donors; the N carry H's that the gate must keep out
# of the metal sphere.
_CISPLATIN_XYZ = """Pt   0.000000     0.000000     0.000000
N    2.050000     0.000000     0.000000
N   -2.050000     0.000000     0.000000
Cl   0.000000     2.300000     0.000000
Cl   0.000000    -2.300000     0.000000
H    2.500000     0.900000     0.000000
H    2.500000    -0.450000     0.780000
H    2.500000    -0.450000    -0.780000
H   -2.500000     0.900000     0.000000
H   -2.500000    -0.450000     0.780000
H   -2.500000    -0.450000    -0.780000
"""

# n-hexane-like chain with several rotatable HEAVY-moving bonds and no metal.
_HEXANE_XYZ = """C    0.000000     0.000000     0.000000
C    1.500000     0.000000     0.000000
C    2.100000     1.400000     0.000000
C    3.600000     1.400000     0.000000
C    4.200000     2.800000     0.000000
C    5.700000     2.800000     0.000000
H   -0.5 0.5 0.8
H   -0.5 0.5 -0.8
H   -0.5 -0.9 0.0
H    1.9 -0.6 0.8
H    1.9 -0.6 -0.8
H    1.7 2.0 0.8
H    1.7 2.0 -0.8
H    4.0 0.8 0.8
H    4.0 0.8 -0.8
H    3.8 3.4 0.8
H    3.8 3.4 -0.8
H    6.1 2.2 0.8
H    6.1 2.2 -0.8
H    6.0 3.8 0.0
"""


def _set_off(monkeypatch):
    for k in list(os.environ):
        if k.startswith("DELFIN_FFFREE_CONF_"):
            monkeypatch.delenv(k, raising=False)


def _set_on(monkeypatch, **extra):
    monkeypatch.setenv("DELFIN_FFFREE_CONF_COMPLETE", "1")
    for k, v in extra.items():
        monkeypatch.setenv(k, str(v))


# --- identity when off ----------------------------------------------------

def test_apply_to_ensemble_identity_when_off(monkeypatch):
    _set_off(monkeypatch)
    iso = [(_HEXANE_XYZ, "a"), (_HEXANE_XYZ, "b")]
    assert CC.apply_to_ensemble(iso) is iso


def test_apply_to_ensemble_identity_explicit_zero(monkeypatch):
    _set_off(monkeypatch)
    monkeypatch.setenv("DELFIN_FFFREE_CONF_COMPLETE", "0")
    iso = [(_HEXANE_XYZ, "a")]
    assert CC.apply_to_ensemble(iso) is iso


# --- DOF detection: heavy-bearing subtree --------------------------------

def test_dofs_move_heavy_atoms_only():
    ob = CC._build_ob_mol_from_xyz(_HEXANE_XYZ)
    g = CC._graph_from_ob(ob)
    if not g:
        pytest.skip("OpenBabel unavailable")
    dofs = CC.identify_complete_dofs(g, max_dofs=8)
    assert dofs, "hexane has rotatable backbone bonds"
    # every reported DOF moves at least one heavy atom beyond the rotor root
    for d in dofs:
        assert d["moved_beyond"] >= 1
        # rotating subtree contains at least one heavy atom
        assert any(g["atomic_nums"][i] > 1 for i in d["rotating"])


# --- topology gate --------------------------------------------------------

def test_gate_rejects_md_break():
    ob = CC._build_ob_mol_from_xyz(_CISPLATIN_XYZ)
    g = CC._graph_from_ob(ob)
    if not g:
        pytest.skip("OpenBabel unavailable")
    syms, base = CC._parse_delfin_xyz(_CISPLATIN_XYZ)
    bb = CC._base_bond_set(g, base)
    md = CC._md_pairs(g, base)
    rb = CC._ring_bonds(g)
    assert md, "Pt complex must have M-D pairs"
    # pull one donor 1.0 A away from the metal -> M-D break
    broken = [list(c) for c in base]
    broken[1][0] += 1.0  # move N at index 1 outward
    assert CC.topology_preserved(g, base, base, bb, md, rb) is True
    assert CC.topology_preserved(g, base, broken, bb, md, rb) is False


def test_gate_rejects_h_in_metal_sphere():
    ob = CC._build_ob_mol_from_xyz(_CISPLATIN_XYZ)
    g = CC._graph_from_ob(ob)
    if not g:
        pytest.skip("OpenBabel unavailable")
    syms, base = CC._parse_delfin_xyz(_CISPLATIN_XYZ)
    bb = CC._base_bond_set(g, base)
    md = CC._md_pairs(g, base)
    rb = CC._ring_bonds(g)
    # move an N-H (index 5) right next to the metal at the origin
    intruder = [list(c) for c in base]
    intruder[5] = [0.9, 0.0, 0.0]  # H ~0.9 A from Pt
    assert CC.topology_preserved(g, base, intruder, bb, md, rb) is False


def test_gate_rejects_spurious_bond():
    ob = CC._build_ob_mol_from_xyz(_CISPLATIN_XYZ)
    g = CC._graph_from_ob(ob)
    if not g:
        pytest.skip("OpenBabel unavailable")
    syms, base = CC._parse_delfin_xyz(_CISPLATIN_XYZ)
    bb = CC._base_bond_set(g, base)
    md = CC._md_pairs(g, base)
    rb = CC._ring_bonds(g)
    # crash the two N atoms together (forge a spurious N-N bond)
    spur = [list(c) for c in base]
    spur[2] = [base[1][0] + 0.3, base[1][1], base[1][2]]  # N idx2 onto N idx1
    # also drag idx2's H along so M-D check on that N may or may not trip;
    # the spurious-bond check alone must reject.
    assert CC.topology_preserved(g, base, spur, bb, md, rb) is False


# --- RMSD dedup -----------------------------------------------------------

def test_kabsch_rmsd_zero_for_identical():
    syms, co = CC._parse_delfin_xyz(_HEXANE_XYZ)
    assert CC._kabsch_rmsd_heavy(syms, co, co) < 1e-6


def test_kabsch_rmsd_invariant_to_rigid_rotation():
    import math
    syms, co = CC._parse_delfin_xyz(_HEXANE_XYZ)
    # rotate the whole molecule 37 deg about z -> RMSD must stay ~0
    th = math.radians(37.0)
    c, s = math.cos(th), math.sin(th)
    rot = [(x * c - y * s, x * s + y * c, z) for (x, y, z) in co]
    assert CC._kabsch_rmsd_heavy(syms, co, rot) < 1e-6


# --- ensemble cap + dedup + determinism ----------------------------------

def test_ensemble_global_cap(monkeypatch):
    _set_on(monkeypatch, DELFIN_FFFREE_CONF_MAX_PER_STRUCT=6)
    # five identical base frames; the global cap must bound total ADDED <= 6
    iso = [(_HEXANE_XYZ, f"base{i}") for i in range(5)]
    out = CC.apply_to_ensemble(iso)
    base = [l for x, l in out if "conf-" not in str(l)]
    conf = [l for x, l in out if "conf-" in str(l)]
    assert len(base) == 5  # every base frame preserved
    assert len(conf) <= 6  # global cap respected


def test_ensemble_deterministic(monkeypatch):
    _set_on(monkeypatch, DELFIN_FFFREE_CONF_MAX_PER_STRUCT=8)
    iso = [(_HEXANE_XYZ, "a")]
    out1 = CC.apply_to_ensemble(iso)
    out2 = CC.apply_to_ensemble(iso)
    assert [x for x, l in out1] == [x for x, l in out2]
    assert [l for x, l in out1] == [l for x, l in out2]


def test_ensemble_dedup_no_near_identical(monkeypatch):
    _set_on(monkeypatch, DELFIN_FFFREE_CONF_RMSD=0.5,
            DELFIN_FFFREE_CONF_MAX_PER_STRUCT=32)
    iso = [(_HEXANE_XYZ, "a")]
    out = CC.apply_to_ensemble(iso)
    coords = [CC._parse_delfin_xyz(x) for x, l in out]
    # no two emitted frames within 0.5 A heavy RMSD
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            si, ci = coords[i]
            sj, cj = coords[j]
            if len(si) == len(sj):
                assert CC._kabsch_rmsd_heavy(si, ci, cj) >= 0.5 - 1e-9


# --- metallacycle-rotation guard (DELFIN_FFFREE_METALLACYCLE_ROT_GUARD) ----
#
# Synthetic graphs (no build) exercise the graph-only metallacycle detection
# and the env-gated DOF exclusion directly.  See conformer_complete docstrings.

def _chelate_graph():
    """Ni ethylenediamine chelate: M(0)-N(1)-C(2)-C(3)-N(4) closed through M.
    Backbone bonds (1,2),(2,3),(3,4) are metallacycle-internal."""
    return {
        "n_atoms": 5,
        "atomic_nums": [28, 7, 6, 6, 7],
        "is_metal": [True, False, False, False, False],
        "neighbours": [[1, 4], [0, 2], [1, 3], [2, 4], [0, 3]],
        "bonds": [
            (0, 1, 1, False, False), (1, 2, 1, False, False),
            (2, 3, 1, False, False), (3, 4, 1, False, False),
            (0, 4, 1, False, False),
        ],
    }


def _pendant_arm_graph():
    """Monodentate pendant arm M(0)-N(1)-C(2)-C(3)-C(4): no ring through M."""
    return {
        "n_atoms": 5,
        "atomic_nums": [28, 7, 6, 6, 6],
        "is_metal": [True, False, False, False, False],
        "neighbours": [[1], [0, 2], [1, 3], [2, 4], [3]],
        "bonds": [
            (0, 1, 1, False, False), (1, 2, 1, False, False),
            (2, 3, 1, False, False), (3, 4, 1, False, False),
        ],
    }


def test_metallacycle_bonds_flags_chelate_backbone():
    g = _chelate_graph()
    mc = CC.metallacycle_bonds(g, None)
    assert mc == {(1, 2), (2, 3), (3, 4)}


def test_metallacycle_bonds_ignores_pendant_arm():
    g = _pendant_arm_graph()
    assert CC.metallacycle_bonds(g, None) == set()


def test_metallacycle_bonds_empty_without_metal():
    g = {
        "n_atoms": 3, "atomic_nums": [6, 6, 6], "is_metal": [False] * 3,
        "neighbours": [[1], [0, 2], [1]],
        "bonds": [(0, 1, 1, False, False), (1, 2, 1, False, False)],
    }
    assert CC.metallacycle_bonds(g, None) == set()


def test_guard_off_byte_identical_dof_set(monkeypatch):
    """Flag unset / '0' -> the DOF set is the base behaviour (no exclusion)."""
    g = _chelate_graph()
    monkeypatch.delenv("DELFIN_FFFREE_METALLACYCLE_ROT_GUARD", raising=False)
    dofs_unset = CC.identify_complete_dofs(g, max_dofs=8, coords=None)
    monkeypatch.setenv("DELFIN_FFFREE_METALLACYCLE_ROT_GUARD", "0")
    dofs_zero = CC.identify_complete_dofs(g, max_dofs=8, coords=None)
    keys_unset = [(d["pivot"], d["anchor"]) for d in dofs_unset]
    keys_zero = [(d["pivot"], d["anchor"]) for d in dofs_zero]
    assert keys_unset == keys_zero  # flag '0' identical to unset


def test_guard_on_excludes_metallacycle_keeps_peripheral(monkeypatch):
    """Chelate ring + a peripheral ethyl arm: ON drops the ring bonds but keeps
    the peripheral rotor (coverage preserved)."""
    # M(0)-N(1)-C(2)-C(3)-N(4) chelate, with a peripheral C5-C6 ethyl on C2.
    g = {
        "n_atoms": 7,
        "atomic_nums": [28, 7, 6, 6, 7, 6, 6],
        "is_metal": [True, False, False, False, False, False, False],
        "neighbours": [[1, 4], [0, 2], [1, 3, 5], [2, 4], [0, 3], [2, 6], [5]],
        "bonds": [
            (0, 1, 1, False, False), (1, 2, 1, False, False),
            (2, 3, 1, False, False), (3, 4, 1, False, False),
            (0, 4, 1, False, False), (2, 5, 1, False, False),
            (5, 6, 1, False, False),
        ],
    }
    mc = CC.metallacycle_bonds(g, None)
    assert (2, 3) in mc and (1, 2) in mc and (3, 4) in mc
    assert (2, 5) not in mc  # peripheral arm is NOT a metallacycle bond
    monkeypatch.setenv("DELFIN_FFFREE_METALLACYCLE_ROT_GUARD", "1")
    on = CC.identify_complete_dofs(g, max_dofs=8, coords=None)
    on_keys = {(d["pivot"], d["anchor"]) for d in on}
    # no metallacycle bond survives as a DOF
    assert not (on_keys & {(2, 3), (3, 2), (1, 2), (2, 1), (3, 4), (4, 3)})
    # the peripheral C2-C5 rotor (moving C5+C6, heavy beyond root) is kept
    assert (5, 2) in on_keys or (2, 5) in on_keys


def test_guard_deterministic(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_METALLACYCLE_ROT_GUARD", "1")
    g = _chelate_graph()
    a = CC.metallacycle_bonds(g, None)
    b = CC.metallacycle_bonds(g, None)
    assert a == b
