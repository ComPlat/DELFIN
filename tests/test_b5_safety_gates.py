"""Tests for B5 v2 safety gates (Patch P-BH-NOGO + Patch P-A5-ENV).

Two-patch safety net for the Baustein 5 PBD post-optimizer:

1. **Patch P-BH-NOGO** — μ²-B-H-M bridging-hydride precheck. When a hydrogen
   atom is bonded to both a boron atom and a metal in the same molecule
   (signature of bridging hydrides such as CANFAY Zr-CN11 borohydride),
   rigid-H tracking is forced OFF for that ``post_optimize_geometry`` call
   regardless of caller request. Dragging the H with one parent during a
   clash push otherwise breaks the bridge (4/9 B-H bonds collapsed to
   ~0.95 Å in the master_v3 pool).

2. **Patch P-A5-ENV** — Global env-flag override for Phase A.5 polyhedron
   projection (``DELFIN_B5_PHASE_A5_ENABLE``). Default 0 forces
   ``enable_symmetry=False`` regardless of caller, protecting against
   regressions if upstream defaults are ever flipped. Set the env var to 1
   to honour the caller's parameter (e.g. for forensik experiments).

Tests skip cleanly when RDKit or the B5 module are missing.
"""
from __future__ import annotations

import os
from typing import List, Sequence, Tuple

import pytest

_post_optimizer = pytest.importorskip(
    "delfin.manta._post_optimizer",
    reason="Baustein 5 post-optimizer module not yet implemented",
)
Chem = pytest.importorskip("rdkit.Chem", reason="RDKit required for tests")
_rdGeometry = pytest.importorskip("rdkit.Geometry", reason="RDKit required")

_has_bridging_bh_hydride = getattr(
    _post_optimizer, "_has_bridging_bh_hydride", None
)
post_optimize_geometry = getattr(_post_optimizer, "post_optimize_geometry", None)
if _has_bridging_bh_hydride is None or post_optimize_geometry is None:
    pytest.skip(
        "B-H precheck or post_optimize_geometry not wired",
        allow_module_level=True,
    )


Atom = Tuple[str, float, float, float]
Bond = Tuple[int, int]


def _xyz(atoms: Sequence[Atom]) -> str:
    lines = [str(len(atoms)), "test"]
    for sym, x, y, z in atoms:
        lines.append(f"{sym:<3s} {x:14.8f} {y:14.8f} {z:14.8f}")
    return "\n".join(lines) + "\n"


def _build_mol(atoms: Sequence[Atom], bonds: Sequence[Bond]):
    mol = Chem.RWMol()
    for sym, _, _, _ in atoms:
        a = Chem.Atom(sym)
        a.SetNoImplicit(True)
        mol.AddAtom(a)
    for i, j in bonds:
        mol.AddBond(i, j, Chem.BondType.SINGLE)
    conf = Chem.Conformer(len(atoms))
    for i, (_, x, y, z) in enumerate(atoms):
        conf.SetAtomPosition(
            i, _rdGeometry.Point3D(float(x), float(y), float(z))
        )
    mol.AddConformer(conf, assignId=True)
    return mol


def _parse(xyz: str) -> List[Atom]:
    lines = xyz.strip().splitlines()
    n = int(lines[0])
    out: List[Atom] = []
    for line in lines[2:2 + n]:
        parts = line.split()
        out.append(
            (parts[0], float(parts[1]), float(parts[2]), float(parts[3]))
        )
    return out


def _dist(a: Atom, b: Atom) -> float:
    return (
        (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2
    ) ** 0.5


# ============================================================================
# Fixture builders
# ============================================================================


def _zr_bh_bridge() -> Tuple[List[Atom], List[Bond]]:
    """Minimal Zr-(μ-H)-B bridging-hydride fixture (CANFAY-style).

    Zr at origin, B at +2.4 Å, bridging H at ~+1.2 Å between them, terminal H
    on B (non-bridging) at the opposite side. Bonds: Zr-H_br, B-H_br, B-H_term.
    No Zr-B direct bond (the bridge goes through the H).
    """
    atoms = [
        ("Zr", 0.0, 0.0, 0.0),
        ("B",  2.40, 0.0, 0.0),
        ("H",  1.20, 0.0, 0.0),   # bridging — bonded to both Zr and B
        ("H",  3.20, 0.0, 0.0),   # terminal B-H, no metal contact
    ]
    bonds = [(0, 2), (1, 2), (1, 3)]
    return atoms, bonds


def _isolated_bh_no_metal() -> Tuple[List[Atom], List[Bond]]:
    """BH₃ fragment with no metal in the molecule — must NOT trigger the gate."""
    atoms = [
        ("B", 0.0, 0.0, 0.0),
        ("H", 1.20, 0.0, 0.0),
        ("H", -0.60, 1.04, 0.0),
        ("H", -0.60, -1.04, 0.0),
    ]
    bonds = [(0, 1), (0, 2), (0, 3)]
    return atoms, bonds


def _pt_ch3_no_boron() -> Tuple[List[Atom], List[Bond]]:
    """Pt-CH₃: H bonded to C and C bonded to Pt — but no boron, so the
    bridging-BH detector must return False."""
    atoms = [
        ("Pt", 0.0, 0.0, 0.0),
        ("C",  2.05, 0.0, 0.0),
        ("H",  2.40, 1.00, 0.0),
        ("H",  2.40, -0.50, 0.86),
        ("H",  2.40, -0.50, -0.86),
    ]
    bonds = [(0, 1), (1, 2), (1, 3), (1, 4)]
    return atoms, bonds


def _pt_nh3_distorted() -> Tuple[List[Atom], List[Bond]]:
    """Reused from test_b5_rigid_h: Pt-NH₃ with Pt-N stretched to 2.15 Å so
    Stage 1 fires. No boron present.
    """
    atoms = [
        ("Pt", 0.0, 0.0, 0.0),
        ("N",  2.15, 0.0, 0.0),
        ("H",  2.45, 0.95, 0.0),
        ("H",  2.45, -0.475, 0.823),
        ("H",  2.45, -0.475, -0.823),
    ]
    bonds = [(0, 1), (1, 2), (1, 3), (1, 4)]
    return atoms, bonds


# ============================================================================
# 1. _has_bridging_bh_hydride direct unit tests
# ============================================================================


def test_bh_precheck_detects_bh_bridge():
    """μ²-B-H-Zr fixture → detector returns True; post_optimize_geometry
    reports ``bh_bridge_nogo=True`` and forces rigid_h OFF."""
    atoms, bonds = _zr_bh_bridge()
    mol = _build_mol(atoms, bonds)
    assert _has_bridging_bh_hydride(mol) is True, (
        "Zr-H-B bridge must be detected by the precheck"
    )
    xyz = _xyz(atoms)
    _, report = post_optimize_geometry(
        xyz, mol, class_label="sigma", max_iter=2,
        step_size=0.1, bond_tol=0.15, rigid_h=True,
    )
    assert report.get("bh_bridge_nogo") is True, (
        "post_optimize_geometry must surface bh_bridge_nogo when bridge detected"
    )


def test_bh_precheck_isolated_bh_no_metal():
    """BH₃ alone (no metal) → detector returns False. Only bridges count."""
    atoms, bonds = _isolated_bh_no_metal()
    mol = _build_mol(atoms, bonds)
    assert _has_bridging_bh_hydride(mol) is False, (
        "BH₃ without a metal must NOT trigger the bridging-hydride gate"
    )


def test_bh_precheck_no_b_atom():
    """Regular Pt-CH₃ molecule → no boron present → detector returns False."""
    atoms, bonds = _pt_ch3_no_boron()
    mol = _build_mol(atoms, bonds)
    assert _has_bridging_bh_hydride(mol) is False, (
        "Pt-CH₃ has no boron — detector must return False"
    )


def test_bh_precheck_none_input_safe():
    """Defensive: ``None`` molecule must not crash, returns False."""
    assert _has_bridging_bh_hydride(None) is False


def test_bh_precheck_b_with_only_terminal_h_no_metal_contact():
    """Boron-H with a metal in the same molecule but the metal is NOT bonded
    to the bridging H — the precheck must NOT trigger.

    Layout: Pt --- N --- B(-H) — Pt is bonded to N, B is bonded to N and H,
    but H is only on B (no M-H bond). This is a non-bridging borane (e.g.
    aminoborane ligand attached via N).
    """
    atoms = [
        ("Pt", 0.0, 0.0, 0.0),
        ("N",  2.10, 0.0, 0.0),
        ("B",  3.60, 0.0, 0.0),
        ("H",  4.80, 0.0, 0.0),  # B-H but NOT Pt-H
    ]
    bonds = [(0, 1), (1, 2), (2, 3)]
    mol = _build_mol(atoms, bonds)
    assert _has_bridging_bh_hydride(mol) is False, (
        "Terminal B-H (not μ²-bridging) must NOT trigger the gate"
    )


def test_bh_precheck_rigid_h_forced_off_for_bridge(monkeypatch):
    """Even with rigid_h=True explicitly requested, the bridging fixture
    must produce output equivalent to rigid_h=False (no H drag occurred).
    We compare to an explicit rigid_h=False run on the same fixture."""
    atoms, bonds = _zr_bh_bridge()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds)
    mol_b = _build_mol(atoms, bonds)
    # rigid_h requested but bridge present → must be auto-disabled.
    out_a, rep_a = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=5,
        step_size=0.1, bond_tol=0.15, rigid_h=True,
    )
    # Explicit rigid_h=False on the identical fixture.
    out_b, rep_b = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=5,
        step_size=0.1, bond_tol=0.15, rigid_h=False,
    )
    assert out_a == out_b, (
        "Bridge fixture with rigid_h=True must yield same XYZ as rigid_h=False"
    )
    assert rep_a.get("bh_bridge_nogo") is True
    # The rigid_h=False run is not auto-blocked — the flag is only set when
    # the auto-disable actually fires (i.e. when rigid_h was requested).
    assert rep_b.get("bh_bridge_nogo", False) is False


# ============================================================================
# 2. Phase A.5 env-gate (Patch P-A5-ENV) — default OFF, env=1 honours param
# ============================================================================


def _square_planar_4_donor_fixture() -> Tuple[List[Atom], List[Bond]]:
    """Cu(II)-like 4-donor sigma fixture that Phase A.5 would otherwise try
    to project onto an ideal symmetric polyhedron (the failure mode for
    Jahn-Teller-distorted d⁹ Cu²⁺).

    Four monodentate Cl donors arranged in a distorted square-planar shape
    (two long Cu-Cl ~ 2.40 Å along x, two shorter ~ 2.20 Å along y) — exactly
    the kind of Jahn-Teller geometry Phase A.5 would idealise.
    """
    atoms = [
        ("Cu", 0.0, 0.0, 0.0),
        ("Cl",  2.40, 0.0, 0.0),
        ("Cl", -2.40, 0.0, 0.0),
        ("Cl",  0.0,  2.20, 0.0),
        ("Cl",  0.0, -2.20, 0.0),
    ]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    return atoms, bonds


def test_phase_a5_env_default_off_ignores_enable_symmetry(monkeypatch):
    """With ``DELFIN_B5_PHASE_A5_ENABLE`` unset (default 0), passing
    ``enable_symmetry=True`` to post_optimize_geometry must have NO effect:
    output is bit-identical to enable_symmetry=False on the same input.
    """
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_ENABLE", raising=False)
    atoms, bonds = _square_planar_4_donor_fixture()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds)
    mol_b = _build_mol(atoms, bonds)
    out_true, _ = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    out_false, _ = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=False,
    )
    assert out_true == out_false, (
        "Default env (unset/0) must force enable_symmetry to no-op — "
        "enable_symmetry=True and =False outputs must match bit-exact"
    )


def test_phase_a5_env_explicit_zero_forces_off(monkeypatch):
    """Setting the env-flag to exactly "0" must behave identically to the
    unset case (force enable_symmetry → False)."""
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_ENABLE", "0")
    atoms, bonds = _square_planar_4_donor_fixture()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds)
    mol_b = _build_mol(atoms, bonds)
    out_true, _ = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    out_false, _ = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=False,
    )
    assert out_true == out_false, (
        "Env=0 must produce identical output for enable_symmetry True/False"
    )


def test_phase_a5_env_enable_respects_param(monkeypatch):
    """Setting ``DELFIN_B5_PHASE_A5_ENABLE=1`` must restore parameter
    sensitivity: ``enable_symmetry=True`` now actually invokes the Phase A.5
    polyhedron-projection helper. We instrument the helper via monkeypatch
    to count invocations — this is a precise, fixture-independent check that
    the env-gate is wired correctly in both directions.
    """
    call_log: List[bool] = []
    real_project = _post_optimizer._project_donors_to_polyhedron

    def _spy(coords, mol, nbrs, metals, md_pairs, nb_pairs, **kw):
        call_log.append(True)
        return real_project(
            coords, mol, nbrs, metals, md_pairs, nb_pairs, **kw
        )

    monkeypatch.setattr(
        _post_optimizer, "_project_donors_to_polyhedron", _spy
    )

    atoms, bonds = _square_planar_4_donor_fixture()
    xyz = _xyz(atoms)

    # env=0: enable_symmetry=True must NOT invoke the projection helper.
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_ENABLE", "0")
    mol_off = _build_mol(atoms, bonds)
    post_optimize_geometry(
        xyz, mol_off, class_label="sigma", max_iter=2,
        step_size=0.1, bond_tol=0.15, enable_symmetry=True,
    )
    n_calls_env_off = len(call_log)
    assert n_calls_env_off == 0, (
        f"Env=0 must suppress Phase A.5 entirely (got {n_calls_env_off} calls)"
    )

    # env=1: enable_symmetry=True must invoke the projection helper.
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_ENABLE", "1")
    mol_on = _build_mol(atoms, bonds)
    post_optimize_geometry(
        xyz, mol_on, class_label="sigma", max_iter=2,
        step_size=0.1, bond_tol=0.15, enable_symmetry=True,
    )
    n_calls_env_on = len(call_log)
    assert n_calls_env_on >= 1, (
        f"Env=1 + enable_symmetry=True must call Phase A.5 "
        f"(got {n_calls_env_on} calls)"
    )

    # env=1 with enable_symmetry=False must still suppress (param respected).
    mol_param_off = _build_mol(atoms, bonds)
    n_before = len(call_log)
    post_optimize_geometry(
        xyz, mol_param_off, class_label="sigma", max_iter=2,
        step_size=0.1, bond_tol=0.15, enable_symmetry=False,
    )
    assert len(call_log) == n_before, (
        "Env=1 with enable_symmetry=False must NOT invoke Phase A.5"
    )


def test_phase_a5_env_invalid_value_falls_back_to_off(monkeypatch):
    """Garbage env value ("foo") must be treated as 0 (off, fail-safe)."""
    monkeypatch.setenv("DELFIN_B5_PHASE_A5_ENABLE", "not-an-int")
    atoms, bonds = _square_planar_4_donor_fixture()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds)
    mol_b = _build_mol(atoms, bonds)
    out_true, _ = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=True,
    )
    out_false, _ = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=3,
        step_size=0.3, bond_tol=0.15, enable_symmetry=False,
    )
    assert out_true == out_false, (
        "Invalid env value must fall back to OFF — enable_symmetry no-op"
    )


# ============================================================================
# 3. Cross-cutting default-OFF contract
# ============================================================================


def test_default_off_bit_exact_no_bridge_no_env(monkeypatch):
    """Default call (no rigid_h, no enable_symmetry, env unset) must be
    bit-identical to a call with both flags explicitly False and env=0 —
    proves the two new gates are no-op on the default path.
    """
    monkeypatch.delenv("DELFIN_B5_PHASE_A5_ENABLE", raising=False)
    atoms, bonds = _pt_nh3_distorted()
    xyz = _xyz(atoms)
    mol_a = _build_mol(atoms, bonds)
    mol_b = _build_mol(atoms, bonds)
    out_default, _ = post_optimize_geometry(
        xyz, mol_a, class_label="sigma", max_iter=5,
        step_size=0.3, bond_tol=0.05,
    )
    out_explicit, _ = post_optimize_geometry(
        xyz, mol_b, class_label="sigma", max_iter=5,
        step_size=0.3, bond_tol=0.05,
        rigid_h=False, enable_symmetry=False,
    )
    assert out_default == out_explicit, (
        "Default call and explicit-False call must produce identical XYZ"
    )
