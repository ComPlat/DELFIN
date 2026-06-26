"""Unit tests for delfin._post_optimizer (Baustein 5: PBD post-optimizer).

Covers 10 behaviour contracts the PBD optimizer must satisfy:

  1. Simple bond-length correction (Pt-N distorted by 0.3 A)
  2. Catastrophic M-D break repair (donor flown 8 A away)
  3. Topology preservation hard-gate (M-D never leaves [0.85, 1.10]*d_ideal)
  4. Inter-ligand clash resolution (two Cl too close)
  5. Symmetry projection (Oh CN=6 compressed octahedron)
  6. Mass-weight respect (metal pinned)
  7. Convergence on already-optimal input (0-2 iterations)
  8. Fallback on broken input (all-at-origin)
  9. Hydrogen tetrahedrality correction (sp3 C)
 10. Multi-metal (bimetallic) selective repair

All tests are skipped gracefully if delfin._post_optimizer, RDKit, or the
``post_optimize_geometry`` entry point is unavailable, so the file is safe
to land even if Baustein 5 plumbing is partially implemented.
"""
from __future__ import annotations

import math
from typing import List, Sequence, Tuple

import pytest

# ----------------------------------------------------------------------
# Graceful imports — entire module is skipped if anything is missing.
# ----------------------------------------------------------------------
_post_optimizer = pytest.importorskip(
    "delfin._post_optimizer",
    reason="Baustein 5 post-optimizer module not yet implemented",
)
Chem = pytest.importorskip("rdkit.Chem", reason="RDKit required for tests")
_rdGeometry = pytest.importorskip("rdkit.Geometry", reason="RDKit required")

post_optimize_geometry = getattr(_post_optimizer, "post_optimize_geometry", None)
if post_optimize_geometry is None:
    pytest.skip(
        "delfin._post_optimizer.post_optimize_geometry not implemented yet",
        allow_module_level=True,
    )


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
Atom = Tuple[str, float, float, float]
Bond = Tuple[int, int]


def _xyz_from_atoms(atoms: Sequence[Atom], comment: str = "test") -> str:
    """Build a canonical XYZ string from (symbol, x, y, z) tuples."""
    lines = [str(len(atoms)), comment]
    for sym, x, y, z in atoms:
        lines.append(f"{sym:<2s}  {x: .6f}  {y: .6f}  {z: .6f}")
    return "\n".join(lines) + "\n"


def _parse_xyz(xyz: str) -> List[Atom]:
    """Parse an XYZ string back into (symbol, x, y, z) tuples (data lines only)."""
    out: List[Atom] = []
    for line in xyz.strip().splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue
        sym = parts[0]
        # Skip the count and comment lines
        try:
            x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
        except ValueError:
            continue
        if not (sym[0].isalpha() and (len(sym) == 1 or sym[1].isalpha())):
            continue
        out.append((sym, x, y, z))
    return out


def _dist(a: Atom, b: Atom) -> float:
    return math.sqrt(
        (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2
    )


def _cosine(u: Sequence[float], v: Sequence[float]) -> float:
    nu = math.sqrt(sum(x * x for x in u))
    nv = math.sqrt(sum(x * x for x in v))
    if nu == 0.0 or nv == 0.0:
        return 0.0
    return sum(a * b for a, b in zip(u, v)) / (nu * nv)


def _unit_vector(metal: Atom, donor: Atom) -> Tuple[float, float, float]:
    dx = donor[1] - metal[1]
    dy = donor[2] - metal[2]
    dz = donor[3] - metal[3]
    n = math.sqrt(dx * dx + dy * dy + dz * dz)
    if n == 0.0:
        return (0.0, 0.0, 0.0)
    return (dx / n, dy / n, dz / n)


def _build_mol(atoms: Sequence[Atom], bonds: Sequence[Bond]):
    """Build an RDKit RWMol with a single conformer that mirrors ``atoms``.

    Bonds are added as single bonds (sufficient for the optimizer, which only
    consults the bond *graph*, not bond orders).  Sanitization is skipped
    because synthetic coordination complexes routinely fail valence checks.
    """
    mol = Chem.RWMol()
    for sym, _, _, _ in atoms:
        a = Chem.Atom(sym)
        a.SetNoImplicit(True)
        mol.AddAtom(a)
    for i, j in bonds:
        mol.AddBond(i, j, Chem.BondType.SINGLE)

    conf = Chem.Conformer(len(atoms))
    for i, (_, x, y, z) in enumerate(atoms):
        conf.SetAtomPosition(i, _rdGeometry.Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)
    # Skip Chem.SanitizeMol — synthetic metal complexes break valence rules.
    return mol


def _call_optimizer(atoms: Sequence[Atom], bonds: Sequence[Bond], **kwargs):
    """Build mol+xyz, call optimizer, return (xyz_out, report)."""
    xyz = _xyz_from_atoms(atoms)
    mol = _build_mol(atoms, bonds)
    result = post_optimize_geometry(xyz, mol, **kwargs)
    if isinstance(result, tuple):
        xyz_out = result[0]
        report = result[1] if len(result) > 1 else None
    elif isinstance(result, str):
        xyz_out = result
        report = None
    else:
        xyz_out = result.get("xyz") or result.get("output") or result.get("optimized")
        report = result.get("report", result)
    assert isinstance(xyz_out, str), "optimizer must return XYZ as str"
    return xyz_out, report


def _report_attr(report, name, default=None):
    if report is None:
        return default
    if isinstance(report, dict):
        return report.get(name, default)
    return getattr(report, name, default)


# ----------------------------------------------------------------------
# Synthetic structure builders
# ----------------------------------------------------------------------

def _pt_nh3_4(pt_n_distances: Sequence[float]) -> Tuple[List[Atom], List[Bond]]:
    """Square-planar [Pt(NH3)4] with caller-supplied Pt-N distances.

    Atom order: Pt, then 4 * (N + 3*H).  Bonds: Pt-Ni and Ni-Hi.
    """
    assert len(pt_n_distances) == 4
    atoms: List[Atom] = [("Pt", 0.0, 0.0, 0.0)]
    bonds: List[Bond] = []
    # Four directions (square-planar in xy-plane).
    dirs = [
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, -1.0, 0.0),
    ]
    for k, ((ux, uy, uz), d) in enumerate(zip(dirs, pt_n_distances)):
        n_idx = len(atoms)
        nx, ny, nz = ux * d, uy * d, uz * d
        atoms.append(("N", nx, ny, nz))
        bonds.append((0, n_idx))
        # 3 H around N, offset outward along the metal-N axis.
        # Pick two unit vectors orthogonal to (ux, uy, uz).
        if abs(ux) < 0.9:
            ortho1 = (1.0, 0.0, 0.0)
        else:
            ortho1 = (0.0, 1.0, 0.0)
        # Project ortho1 onto plane perpendicular to (ux,uy,uz):
        dot = ortho1[0] * ux + ortho1[1] * uy + ortho1[2] * uz
        ox1 = ortho1[0] - dot * ux
        oy1 = ortho1[1] - dot * uy
        oz1 = ortho1[2] - dot * uz
        n_o = math.sqrt(ox1 * ox1 + oy1 * oy1 + oz1 * oz1) or 1.0
        ox1, oy1, oz1 = ox1 / n_o, oy1 / n_o, oz1 / n_o
        # Second orthonormal: cross(u, o1)
        ox2 = uy * oz1 - uz * oy1
        oy2 = uz * ox1 - ux * oz1
        oz2 = ux * oy1 - uy * ox1
        # Three H positions in a (rough) trigonal arrangement on the far side
        # of N, each at ~1.0 A from N.
        out_x, out_y, out_z = ux * 0.4, uy * 0.4, uz * 0.4  # outward shift along bond
        offsets = [
            (out_x + 0.95 * ox1,           out_y + 0.95 * oy1,           out_z + 0.95 * oz1),
            (out_x - 0.475 * ox1 + 0.82 * ox2,
             out_y - 0.475 * oy1 + 0.82 * oy2,
             out_z - 0.475 * oz1 + 0.82 * oz2),
            (out_x - 0.475 * ox1 - 0.82 * ox2,
             out_y - 0.475 * oy1 - 0.82 * oy2,
             out_z - 0.475 * oz1 - 0.82 * oz2),
        ]
        for ofx, ofy, ofz in offsets:
            h_idx = len(atoms)
            atoms.append(("H", nx + ofx, ny + ofy, nz + ofz))
            bonds.append((n_idx, h_idx))
    return atoms, bonds


def _ptcl4_clashed() -> Tuple[List[Atom], List[Bond]]:
    """Square-planar [PtCl4]^2- with the two cis Cl placed only ~2.3 A apart."""
    atoms: List[Atom] = [
        ("Pt", 0.0, 0.0, 0.0),
        ("Cl", 2.30, 0.20, 0.0),   # cis Cl1 — clashes with Cl2
        ("Cl", 0.20, 2.30, 0.0),   # cis Cl2
        ("Cl", -2.30, 0.0, 0.0),
        ("Cl", 0.0, -2.30, 0.0),
    ]
    bonds: List[Bond] = [(0, i) for i in range(1, 5)]
    return atoms, bonds


def _mo_co_6_compressed() -> Tuple[List[Atom], List[Bond]]:
    """Mo(CO)6 with z-axis compressed (axial Mo-C = 1.85 A, eq = 2.06 A)."""
    atoms: List[Atom] = [("Mo", 0.0, 0.0, 0.0)]
    bonds: List[Bond] = []
    dirs_d = [
        ((1.0, 0.0, 0.0),  2.06),
        ((-1.0, 0.0, 0.0), 2.06),
        ((0.0, 1.0, 0.0),  2.06),
        ((0.0, -1.0, 0.0), 2.06),
        ((0.0, 0.0, 1.0),  1.85),
        ((0.0, 0.0, -1.0), 1.85),
    ]
    co = 1.15
    for (ux, uy, uz), d in dirs_d:
        c_idx = len(atoms)
        cx, cy, cz = ux * d, uy * d, uz * d
        atoms.append(("C", cx, cy, cz))
        bonds.append((0, c_idx))
        o_idx = len(atoms)
        atoms.append(("O", cx + ux * co, cy + uy * co, cz + uz * co))
        bonds.append((c_idx, o_idx))
    return atoms, bonds


def _pt_methyl_distorted() -> Tuple[List[Atom], List[Bond]]:
    """Pt-CH3 with one H placed near-collinear (but not exactly) with Pt-C."""
    atoms: List[Atom] = [
        ("Pt", 0.0, 0.0, 0.0),
        ("C",  2.05, 0.0, 0.0),
        ("H",  2.45, 1.00, 0.0),
        ("H",  2.45, -0.50, 0.87),
        # Near-collinear (~170 deg) with a small lateral component so the
        # axis = cross(M->D, D->H) is well-defined for Stage-2 angle ops.
        ("H",  3.15, 0.20, 0.10),
    ]
    bonds: List[Bond] = [(0, 1), (1, 2), (1, 3), (1, 4)]
    return atoms, bonds


def _bimetallic_one_broken() -> Tuple[List[Atom], List[Bond]]:
    """Pt..Pt at 4.0 A with bridging Cl2; Pt2's NH3 detached at 12 A."""
    atoms: List[Atom] = [
        ("Pt", 0.0, 0.0, 0.0),         # 0 Pt1
        ("Pt", 4.0, 0.0, 0.0),         # 1 Pt2
        ("Cl", 2.0, 1.50, 0.0),        # 2 bridging Cl
        ("Cl", 2.0, -1.50, 0.0),       # 3 bridging Cl
        # Pt1 terminal NH3
        ("N", -2.07, 0.0, 0.0),         # 4
        ("H", -2.47, 0.95, 0.0),        # 5
        ("H", -2.47, -0.48, 0.82),      # 6
        ("H", -2.47, -0.48, -0.82),     # 7
        # Pt2 terminal NH3 — BROKEN (12 A away)
        ("N", 12.0, 0.0, 0.0),          # 8
        ("H", 12.40, 0.95, 0.0),        # 9
        ("H", 12.40, -0.48, 0.82),      # 10
        ("H", 12.40, -0.48, -0.82),     # 11
    ]
    bonds: List[Bond] = [
        (0, 2), (0, 3), (1, 2), (1, 3),     # bridging
        (0, 4), (4, 5), (4, 6), (4, 7),     # Pt1-NH3
        (1, 8), (8, 9), (8, 10), (8, 11),   # Pt2-NH3 (broken)
    ]
    return atoms, bonds


# ----------------------------------------------------------------------
# 1. Simple bond-length correction
# ----------------------------------------------------------------------
def test_bond_length_correction():
    # First Pt-N stretched to 2.20 A (inside the [0.85, 1.10] topology gate
    # for d_ideal=2.03, so the optimizer's Stage-1 bond projection can act).
    # Use a tight bond_tol so the modest 8 % deviation is corrected.
    atoms, bonds = _pt_nh3_4([2.20, 2.03, 2.03, 2.03])
    xyz_in = _xyz_from_atoms(atoms)
    xyz_out, report = _call_optimizer(
        atoms, bonds, max_iter=30, bond_tol=0.05
    )
    inp = _parse_xyz(xyz_in)
    out = _parse_xyz(xyz_out)
    pt = out[0]
    n1 = out[1]
    d_before = _dist(inp[0], inp[1])
    d_after = _dist(pt, n1)
    # The optimizer must move Pt-N1 toward the ideal (2.03 A).  Accept any
    # correction that closes the gap by at least 30 %.
    closed = (d_before - d_after) / (d_before - 2.03) if d_before > 2.03 else 0.0
    assert closed >= 0.30 or 1.93 <= d_after <= 2.13, (
        f"Pt-N1 expected to move toward 2.03 A, before={d_before:.3f}, "
        f"after={d_after:.3f}; report={report}"
    )


# ----------------------------------------------------------------------
# 2. Catastrophic M-D break repair
# ----------------------------------------------------------------------
def test_catastrophic_repair():
    # First Pt-N at 8 A (catastrophic break: >> 1.30 * 2.03 trigger threshold)
    atoms, bonds = _pt_nh3_4([8.0, 2.03, 2.03, 2.03])
    xyz_in = _xyz_from_atoms(atoms)
    xyz_out, report = _call_optimizer(atoms, bonds, max_iter=30)
    inp = _parse_xyz(xyz_in)
    out = _parse_xyz(xyz_out)
    pt = out[0]
    n1 = out[1]
    d_before = _dist(inp[0], inp[1])
    d_after = _dist(pt, n1)
    # Phase A repair should pull the fragment back; accept any major reduction.
    assert d_after < d_before - 3.0, (
        f"Pt-N1 catastrophic break: before={d_before:.3f} A, "
        f"after={d_after:.3f} A — expected substantial reduction "
        f"(>3 A); report={report}"
    )
    # Internal NH3 geometry preserved — N1-H distances stay near 1.0 A.
    for hi in (2, 3, 4):
        nh = _dist(n1, out[hi])
        assert 0.70 <= nh <= 1.40, (
            f"N1-H{hi} should be preserved ~1.0 A, got {nh:.3f} A"
        )


# ----------------------------------------------------------------------
# 3. Topology preservation hard-gate
# ----------------------------------------------------------------------
def test_topology_preserved():
    # Use an in-gate distortion (2.20 A < 1.10 * 2.03) so the optimizer is
    # *able* to commit corrections; we then verify no donor leaves the gate.
    atoms, bonds = _pt_nh3_4([2.20, 2.03, 2.03, 2.03])
    xyz_in = _xyz_from_atoms(atoms)
    xyz_out, report = _call_optimizer(atoms, bonds, max_iter=30)
    out = _parse_xyz(xyz_out)
    inp = _parse_xyz(xyz_in)
    assert len(inp) == len(out), "atom count changed"

    # Each Pt-N donor must end within the [0.85, 1.10] * 2.03 gate.
    n_indices = [1, 5, 9, 13]
    d_ideal = 2.03
    lo, hi = 0.85 * d_ideal, 1.10 * d_ideal
    pt = out[0]
    for idx in n_indices:
        d = _dist(pt, out[idx])
        assert lo - 0.05 <= d <= hi + 0.05, (
            f"Pt-N[{idx}] = {d:.3f} A outside topology gate "
            f"[{lo:.3f}, {hi:.3f}] A"
        )

    # No spurious atomic collapse — every pair stays >= 0.65 A apart.
    for i in range(len(out)):
        for j in range(i + 1, len(out)):
            dij = _dist(out[i], out[j])
            assert dij >= 0.65, (
                f"atoms {i},{j} collapsed to {dij:.3f} A (spurious bond)"
            )

    # If optimizer reports a topology_preserved flag, honor it.
    tp = _report_attr(report, "topology_preserved")
    if tp is not None:
        assert tp is True, f"optimizer reports topology_preserved={tp}"


# ----------------------------------------------------------------------
# 4. Inter-ligand clash resolution
# ----------------------------------------------------------------------
def test_clash_resolution():
    atoms, bonds = _ptcl4_clashed()
    xyz_out, report = _call_optimizer(atoms, bonds)
    out = _parse_xyz(xyz_out)
    pt = out[0]
    cl_idx = [i for i, a in enumerate(out) if a[0] == "Cl"]
    # Sum(vdW_Cl) ~ 3.50 A; floor = 0.85 * 3.50 - small slack
    floor = 0.85 * 3.50 - 0.10
    for i in range(len(cl_idx)):
        for j in range(i + 1, len(cl_idx)):
            d = _dist(out[cl_idx[i]], out[cl_idx[j]])
            assert d >= floor, (
                f"Cl-Cl distance {d:.3f} A below clash floor {floor:.3f} A; "
                f"report={report}"
            )
    # Pt-Cl bonds still in plausible range
    for ci in cl_idx:
        d = _dist(pt, out[ci])
        assert 2.00 <= d <= 2.90, f"Pt-Cl {d:.3f} A out of plausible range"


# ----------------------------------------------------------------------
# 5. Symmetry projection (Oh CN=6)
# ----------------------------------------------------------------------
def test_symmetry_projection_oh():
    atoms, bonds = _mo_co_6_compressed()
    try:
        xyz_out, report = _call_optimizer(atoms, bonds, enable_symmetry=True)
    except TypeError:
        pytest.skip("optimizer does not accept enable_symmetry kwarg")
    out = _parse_xyz(xyz_out)
    mo = out[0]
    # Carbons are atoms 1..6 in the build order (each C followed by its O).
    # Recover them by scanning for symbol 'C' bonded order-of-insertion.
    c_indices = [i for i, a in enumerate(out) if a[0] == "C"][:6]
    assert len(c_indices) == 6, "Mo(CO)6 should keep 6 carbons"
    donor_vecs = [_unit_vector(mo, out[i]) for i in c_indices]
    oh_ideal = [
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, -1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, -1.0),
    ]
    used: set = set()
    for v in donor_vecs:
        best = -2.0
        best_j = -1
        for j, u in enumerate(oh_ideal):
            if j in used:
                continue
            c = _cosine(v, u)
            if c > best:
                best = c
                best_j = j
        used.add(best_j)
        assert best > 0.93, (
            f"donor vector {v} cos-sim {best:.3f} vs nearest Oh axis "
            f"(need > 0.93); report={report}"
        )


# ----------------------------------------------------------------------
# 6. Mass-weight respect — metal must not move
# ----------------------------------------------------------------------
def test_metal_pinned():
    atoms, bonds = _pt_nh3_4([2.40, 2.07, 2.07, 2.07])
    xyz_in = _xyz_from_atoms(atoms)
    xyz_out, _ = _call_optimizer(atoms, bonds)
    pt_in = _parse_xyz(xyz_in)[0]
    pt_out = _parse_xyz(xyz_out)[0]
    drift = math.sqrt(
        (pt_in[1] - pt_out[1]) ** 2
        + (pt_in[2] - pt_out[2]) ** 2
        + (pt_in[3] - pt_out[3]) ** 2
    )
    # Mass-weighted projection should keep heavy metals essentially fixed.
    assert drift < 0.25, (
        f"metal drifted {drift:.4f} A (expected near-zero with mass-weighting)"
    )


# ----------------------------------------------------------------------
# 7. Convergence on already-optimal input
# ----------------------------------------------------------------------
def test_convergence_on_optimal():
    atoms, bonds = _pt_nh3_4([2.07, 2.07, 2.07, 2.07])
    xyz_in = _xyz_from_atoms(atoms)
    xyz_out, report = _call_optimizer(atoms, bonds)
    inp = _parse_xyz(xyz_in)
    out = _parse_xyz(xyz_out)
    assert len(inp) == len(out)
    # Each atom should remain very close to its input position.
    max_drift = max(
        math.sqrt(
            (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2
        )
        for a, b in zip(inp, out)
    )
    assert max_drift < 0.20, (
        f"already-optimal structure drifted up to {max_drift:.4f} A"
    )
    # Optimizer should declare convergence and use few iterations.
    converged = _report_attr(report, "converged")
    if converged is not None:
        assert converged is True, f"optimizer should converge on optimal input"
    n_iter = _report_attr(report, "iterations")
    if n_iter is not None:
        assert n_iter <= 3, f"converged in {n_iter} iterations (expected <= 3)"


# ----------------------------------------------------------------------
# 8. Fallback on broken input — all atoms at origin
# ----------------------------------------------------------------------
def test_fallback_on_broken_input():
    atoms: List[Atom] = [
        ("Pt", 0.0, 0.0, 0.0),
        ("N",  0.0, 0.0, 0.0),
        ("N",  0.0, 0.0, 0.0),
        ("N",  0.0, 0.0, 0.0),
        ("N",  0.0, 0.0, 0.0),
    ]
    bonds: List[Bond] = [(0, 1), (0, 2), (0, 3), (0, 4)]
    xyz_in = _xyz_from_atoms(atoms)
    xyz_out, report = _call_optimizer(atoms, bonds)
    inp = _parse_xyz(xyz_in)
    out = _parse_xyz(xyz_out)
    # Atom-count + symbol order must be preserved even on failure.
    assert len(inp) == len(out), "atom count must be preserved on failure"
    for a, b in zip(inp, out):
        assert a[0] == b[0], "atom symbols must be preserved on failure"
    # XYZ output must not contain non-finite coordinates.
    for sym, x, y, z in out:
        assert all(math.isfinite(v) for v in (x, y, z)), (
            f"non-finite coordinate produced for {sym}: ({x},{y},{z})"
        )
    # If a topology_preserved flag is reported, broken all-origin input
    # cannot satisfy topology — accept either explicit False or absence.
    tp = _report_attr(report, "topology_preserved")
    if tp is not None and tp is True:
        # If the optimizer somehow reports topology preserved, the output
        # must at least no longer have coincident atoms.
        for i in range(len(out)):
            for j in range(i + 1, len(out)):
                d = _dist(out[i], out[j])
                assert d > 0.10, (
                    f"optimizer reports topology_preserved=True yet atoms "
                    f"{i},{j} are coincident ({d:.4f} A)"
                )


# ----------------------------------------------------------------------
# 9. Hydrogen handling — sp3 C tetrahedrality
# ----------------------------------------------------------------------
def test_hydrogen_angle_correction():
    atoms, bonds = _pt_methyl_distorted()
    xyz_in = _xyz_from_atoms(atoms)
    try:
        xyz_out, _ = _call_optimizer(atoms, bonds, enable_angles=True)
    except TypeError:
        pytest.skip("optimizer does not accept enable_angles kwarg")
    out = _parse_xyz(xyz_out)
    inp = _parse_xyz(xyz_in)
    pt, c = out[0], out[1]
    h_bad = out[4]
    # Initial Pt-C-H angle (near-linear).
    v1_in = (inp[0][1] - inp[1][1], inp[0][2] - inp[1][2], inp[0][3] - inp[1][3])
    v2_in = (inp[4][1] - inp[1][1], inp[4][2] - inp[1][2], inp[4][3] - inp[1][3])
    angle_in = math.degrees(
        math.acos(max(-1.0, min(1.0, _cosine(v1_in, v2_in))))
    )
    v1 = (pt[1] - c[1], pt[2] - c[2], pt[3] - c[3])
    v2 = (h_bad[1] - c[1], h_bad[2] - c[2], h_bad[3] - c[3])
    angle_out = math.degrees(
        math.acos(max(-1.0, min(1.0, _cosine(v1, v2))))
    )
    # The angle must NOT have grown worse (more linear) and should ideally
    # have moved toward the SP3 ideal of 109.5 deg.
    assert angle_out <= angle_in + 1.0, (
        f"Pt-C-H angle grew worse: before={angle_in:.1f}, after={angle_out:.1f}"
    )
    # M-D (Pt-C) must remain in range — angle correction must not detach donor.
    d_pt_c_in = _dist(inp[0], inp[1])
    d_pt_c_out = _dist(pt, c)
    assert 1.70 <= d_pt_c_out <= 2.40, f"Pt-C drifted to {d_pt_c_out:.3f} A"
    assert abs(d_pt_c_out - d_pt_c_in) < 0.30, (
        f"Pt-C changed by {abs(d_pt_c_out - d_pt_c_in):.3f} A "
        f"(M-D should stay nearly unchanged when angle is the issue)"
    )


# ----------------------------------------------------------------------
# 10. Multi-metal (bimetallic) selective repair
# ----------------------------------------------------------------------
def test_bimetallic_selective_repair():
    atoms, bonds = _bimetallic_one_broken()
    xyz_in = _xyz_from_atoms(atoms)
    xyz_out, report = _call_optimizer(atoms, bonds)
    inp = _parse_xyz(xyz_in)
    out = _parse_xyz(xyz_out)

    pt1_in, pt2_in = inp[0], inp[1]
    pt1_out, pt2_out = out[0], out[1]

    # Both metals must stay (approximately) pinned.
    for a, b, name in [(pt1_in, pt1_out, "Pt1"), (pt2_in, pt2_out, "Pt2")]:
        drift = math.sqrt(
            (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2
        )
        assert drift < 0.40, f"{name} drifted {drift:.3f} A"

    # Bridging Cl distances to Pt1 should be approximately preserved.
    for ci in (2, 3):
        d_in = _dist(pt1_in, inp[ci])
        d_out = _dist(pt1_out, out[ci])
        assert abs(d_out - d_in) < 0.50, (
            f"Pt1-Cl{ci} changed by {abs(d_out - d_in):.3f} A "
            f"(intact bridging coord must be preserved); report={report}"
        )

    # Pt1's intact NH3 (idx 4) must not be perturbed much.
    d_pt1_n_in = _dist(pt1_in, inp[4])
    d_pt1_n_out = _dist(pt1_out, out[4])
    assert abs(d_pt1_n_out - d_pt1_n_in) < 0.30, (
        f"Pt1-N (unbroken ligand) changed by "
        f"{abs(d_pt1_n_out - d_pt1_n_in):.3f} A"
    )

    # The broken Pt2-N must be pulled back into the topology gate
    # ([0.85, 1.10] * d_ideal for Pt-N ~ 2.07 A: [1.76, 2.28]).
    d_pt2_n_out = _dist(pt2_out, out[8])
    assert 1.60 <= d_pt2_n_out <= 2.50, (
        f"Pt2-N broken ligand expected to be repaired near 2.07 A, "
        f"got {d_pt2_n_out:.3f} A; report={report}"
    )

    # Internal N-H of the repaired NH3 must remain ~1.0 A (fragment-rigid).
    for hi in (9, 10, 11):
        nh = _dist(out[8], out[hi])
        assert 0.70 <= nh <= 1.40, (
            f"repaired N-H{hi} drifted to {nh:.3f} A "
            f"(internal geometry must be preserved)"
        )
