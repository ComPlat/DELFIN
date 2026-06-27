"""Baustein 6: symmetry-operation matrices per molecular point group.

This module supplies the 3x3 operation matrices (rotations, reflections,
inversion, improper rotations Sn) of common molecular point groups used by
Baustein 6 Tier C (fragment-archetype symmetry) and Tier D (global molecular
symmetry) enforcement.

Operations are returned in a CANONICAL frame:

  * Principal Cn axis along +z.
  * For Cnv groups: one sigma_v plane is the xz-plane.
  * For Dn / Dnh: one C2' axis along +x.
  * For Td/Th/Oh: aligned with cubic xyz axes.

The caller is responsible for aligning the molecule to this frame (e.g.
via inertia tensor or RMSD-fit of an idealized template) before applying
the operations.

The implementation is pure numpy (no scipy required).  Built operations
are lazily cached on first request per point-group name.

Spec source: ``iters/BAUSTEIN6_MASTERPLAN.md`` Section 2.

No SMILES-, refcode-, or element-list shortcuts.  All operations follow
universal group-theory recipes (Rodrigues for rotations, Householder for
reflections, Sn = sigma_h @ Cn).
"""

from __future__ import annotations

from typing import Callable, Dict, List

import numpy as np

# ---------------------------------------------------------------------------
# Helper builders
# ---------------------------------------------------------------------------


def rotation_matrix_axis_angle(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """3x3 rotation matrix via Rodrigues' formula.

    Parameters
    ----------
    axis : (3,) array-like
        Rotation axis (need not be normalized).
    angle_rad : float
        Rotation angle in radians (right-hand rule).
    """
    axis = np.asarray(axis, dtype=float)
    n = float(np.linalg.norm(axis))
    if n == 0.0:
        raise ValueError("rotation_matrix_axis_angle: zero-norm axis")
    axis = axis / n
    K = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    s = float(np.sin(angle_rad))
    c = float(np.cos(angle_rad))
    R = np.eye(3) + s * K + (1.0 - c) * (K @ K)
    return R


def reflection_matrix(normal: np.ndarray) -> np.ndarray:
    """3x3 reflection (Householder) through plane with given normal."""
    n_vec = np.asarray(normal, dtype=float)
    norm = float(np.linalg.norm(n_vec))
    if norm == 0.0:
        raise ValueError("reflection_matrix: zero-norm plane normal")
    n_vec = n_vec / norm
    return np.eye(3) - 2.0 * np.outer(n_vec, n_vec)


def inversion_matrix() -> np.ndarray:
    """3x3 inversion = -I."""
    return -np.eye(3)


def sn_matrix(axis: np.ndarray, n: int) -> np.ndarray:
    """S_n improper rotation = rotation by 2 pi / n then reflection
    through plane perpendicular to axis.

    Equivalent: sigma_h @ C_n where sigma_h's normal is the axis.
    """
    if n == 0:
        raise ValueError("sn_matrix: n must be non-zero")
    rot = rotation_matrix_axis_angle(axis, 2.0 * np.pi / n)
    refl = reflection_matrix(axis)
    return refl @ rot


def sn_power_matrix(axis: np.ndarray, n: int, k: int) -> np.ndarray:
    """S_n^k = (sigma_h @ C_n)^k.

    Note S_n^k alternates between proper rotation (even k) and improper
    rotation (odd k).  Use only odd k for new operations not already
    covered by Cn powers.
    """
    rot = rotation_matrix_axis_angle(axis, 2.0 * np.pi * k / n)
    if k % 2 == 0:
        return rot
    refl = reflection_matrix(axis)
    return refl @ rot


# ---------------------------------------------------------------------------
# Canonical axes
# ---------------------------------------------------------------------------

_E = np.eye(3)
_X = np.array([1.0, 0.0, 0.0])
_Y = np.array([0.0, 1.0, 0.0])
_Z = np.array([0.0, 0.0, 1.0])


def _cn(axis: np.ndarray, n: int, k: int = 1) -> np.ndarray:
    """C_n^k = rotation by 2 pi k / n about axis."""
    return rotation_matrix_axis_angle(axis, 2.0 * np.pi * k / n)


def _cn_powers(axis: np.ndarray, n: int) -> List[np.ndarray]:
    """All proper rotations of order dividing n about axis: C_n, C_n^2, ..., C_n^(n-1).

    Identity NOT included.  Caller prepends E separately.
    """
    return [_cn(axis, n, k) for k in range(1, n)]


def _sn_odd_powers(axis: np.ndarray, n: int) -> List[np.ndarray]:
    """Odd powers of S_n that yield new improper operations.

    For S_n with n even: S_n^k for k = 1, 3, 5, ..., n - 1 are improper.
    For S_n with n odd: S_n^k for k = 1, 3, ..., 2n - 1 (mod 2n).  Full
    cycle has order 2n.  Even k => proper (already in Cn).  Odd k =>
    improper.  For n=3: S3, S3^3=sigma_h, S3^5 are improper distinct ops.
    """
    if n % 2 == 0:
        return [sn_power_matrix(axis, n, k) for k in range(1, n, 2)]
    # n odd: cycle of length 2n; odd k from 1 to 2n-1
    return [sn_power_matrix(axis, n, k) for k in range(1, 2 * n, 2)]


# ---------------------------------------------------------------------------
# Cyclic groups
# ---------------------------------------------------------------------------


def _build_C1():
    return [_E.copy()]


def _build_Cs():
    # E, sigma_h (reflection through xy-plane => normal = z)
    return [_E.copy(), reflection_matrix(_Z)]


def _build_Ci():
    return [_E.copy(), inversion_matrix()]


def _build_Cn(n: int) -> List[np.ndarray]:
    """C_n: E, C_n, C_n^2, ..., C_n^(n-1)."""
    return [_E.copy()] + _cn_powers(_Z, n)


def _build_C2():
    return _build_Cn(2)


def _build_C3():
    return _build_Cn(3)


def _build_C4():
    return _build_Cn(4)


def _build_C5():
    return _build_Cn(5)


def _build_C6():
    return _build_Cn(6)


# ---------------------------------------------------------------------------
# Cnv groups (cyclic + vertical mirrors)
# ---------------------------------------------------------------------------


def _vertical_mirror_normals(n: int) -> List[np.ndarray]:
    """Normals of n vertical mirror planes containing the z-axis.

    First plane is xz (normal = y).  Subsequent planes are rotated by
    pi / n about z, yielding n distinct planes for Cnv:
      * 2 sigma_v for C2v
      * 3 sigma_v for C3v
      * For C4v: 2 sigma_v + 2 sigma_d (4 total)
      * For C5v: 5 sigma_v
      * For C6v: 3 sigma_v + 3 sigma_d (6 total)
    """
    normals: List[np.ndarray] = []
    for k in range(n):
        angle = np.pi * k / n
        # normal rotates in xy-plane; start with +y so first plane is xz
        nx = -np.sin(angle)
        ny = np.cos(angle)
        normals.append(np.array([nx, ny, 0.0]))
    return normals


def _build_Cnv(n: int) -> List[np.ndarray]:
    ops = _build_Cn(n)
    for nrm in _vertical_mirror_normals(n):
        ops.append(reflection_matrix(nrm))
    return ops


def _build_C2v():
    return _build_Cnv(2)


def _build_C3v():
    return _build_Cnv(3)


def _build_C4v():
    return _build_Cnv(4)


def _build_C5v():
    return _build_Cnv(5)


def _build_C6v():
    return _build_Cnv(6)


# ---------------------------------------------------------------------------
# Cnh groups (cyclic + horizontal mirror)
# ---------------------------------------------------------------------------


def _build_Cnh(n: int) -> List[np.ndarray]:
    """C_nh = C_n x { E, sigma_h }.  Yields 2n operations.

    For odd n, S_n powers replace inversion-style ops; for even n,
    inversion appears as one of the S_n^k.
    """
    cn = _build_Cn(n)
    sigma_h = reflection_matrix(_Z)
    return cn + [sigma_h @ R for R in cn]


def _build_C2h():
    # E, C2(z), i, sigma_h
    return _build_Cnh(2)


def _build_C3h():
    # E, C3, C3^2, sigma_h, S3, S3^5
    return _build_Cnh(3)


# ---------------------------------------------------------------------------
# Dn groups (dihedral, proper rotations only)
# ---------------------------------------------------------------------------


def _perpendicular_C2_axes(n: int) -> List[np.ndarray]:
    """n C2 axes perpendicular to z, distributed every pi / n.

    First axis is +x.  For even n, n / 2 axes are 'C2_prime' (vertices)
    and n / 2 are 'C2_double_prime' (edge midpoints); we return all n.
    """
    axes: List[np.ndarray] = []
    for k in range(n):
        angle = np.pi * k / n
        axes.append(np.array([np.cos(angle), np.sin(angle), 0.0]))
    return axes


def _build_Dn(n: int) -> List[np.ndarray]:
    ops = _build_Cn(n)
    for ax in _perpendicular_C2_axes(n):
        ops.append(rotation_matrix_axis_angle(ax, np.pi))
    return ops


def _build_D2():
    return _build_Dn(2)


def _build_D3():
    return _build_Dn(3)


def _build_D4():
    return _build_Dn(4)


def _build_D5():
    return _build_Dn(5)


def _build_D6():
    return _build_Dn(6)


# ---------------------------------------------------------------------------
# Dnh groups (dihedral + horizontal mirror)
# ---------------------------------------------------------------------------


def _build_Dnh(n: int) -> List[np.ndarray]:
    """D_nh = D_n x { E, sigma_h }.  Yields 4n operations."""
    dn = _build_Dn(n)
    sigma_h = reflection_matrix(_Z)
    return dn + [sigma_h @ R for R in dn]


def _build_D2h():
    # 8 ops: E, C2(z), C2(x), C2(y), i, sigma(xy), sigma(xz), sigma(yz)
    return _build_Dnh(2)


def _build_D3h():
    # 12 ops: E, 2C3, 3C2, sigma_h, 2S3, 3sigma_v
    return _build_Dnh(3)


def _build_D4h():
    # 16 ops
    return _build_Dnh(4)


def _build_D5h():
    # 20 ops
    return _build_Dnh(5)


def _build_D6h():
    # 24 ops
    return _build_Dnh(6)


# ---------------------------------------------------------------------------
# Dnd groups (dihedral + diagonal mirrors, no horizontal mirror)
# ---------------------------------------------------------------------------


def _build_Dnd(n: int) -> List[np.ndarray]:
    """D_nd = D_n + S_(2n) operations + n sigma_d planes.

    Construction: take D_n proper rotations (2n ops) and combine with
    improper operations from the S_(2n) cycle.  Equivalently:
      D_nd = D_n union { S_(2n)^k @ E for odd k in [1, 2n-1] } union { n sigma_d }
    Simpler: D_nd = D_n union (S_(2n) odd-power improper ops) where the
    odd-power Sn ops automatically include the n sigma_d planes (since
    S_(2n)^n = sigma_h for even 2n => actually inversion for D3d, etc.).

    A clean explicit construction: enumerate the 4n group elements as
    D_n proper rotations + (S_(2n))^k for k = 1, 3, ..., 2n - 1 (n
    improper ops) + n sigma_d reflections.  But |D_nd| = 4n, so the
    correct count is 2n proper + 2n improper.  The 2n improper ops
    decompose as n S_(2n)^odd rotations + n sigma_d reflections.
    """
    proper = _build_Dn(n)  # 2n proper ops

    improper: List[np.ndarray] = []
    # n S_(2n) odd-power improper rotations (k = 1, 3, ..., 2n - 1)
    # S_(2n)^k = sigma_h @ C_(2n)^k = -I @ ... ; we use the helper
    for k in range(1, 2 * n, 2):
        improper.append(sn_power_matrix(_Z, 2 * n, k))

    # n sigma_d diagonal mirror planes: planes containing z-axis and
    # bisecting adjacent C2' axes.  First sigma_d normal at angle
    # pi / (2n) from +y (since first C2' is +x), then rotated by pi/n.
    for k in range(n):
        angle = np.pi / (2 * n) + np.pi * k / n
        nrm = np.array([-np.sin(angle), np.cos(angle), 0.0])
        improper.append(reflection_matrix(nrm))

    return proper + improper


def _build_D2d():
    # 8 ops: E, 2S4, C2(z), 2C2', 2sigma_d
    return _build_Dnd(2)


def _build_D3d():
    # 12 ops: E, 2C3, 3C2, i, 2S6, 3sigma_d
    return _build_Dnd(3)


def _build_D4d():
    # 16 ops
    return _build_Dnd(4)


def _build_D5d():
    # 20 ops
    return _build_Dnd(5)


# ---------------------------------------------------------------------------
# Cubic groups (Td, Th, Oh)
# ---------------------------------------------------------------------------

# Cubic frame conventions:
#   * Three C4 axes (for Oh): along +x, +y, +z.
#   * Four C3 axes (Td, Th, Oh): body diagonals (+/-1, +/-1, +/-1).
#   * Six C2' axes (Oh): edge midpoints (face-diagonal directions of
#     cube faces), e.g. (+1, +1, 0), (+1, -1, 0), etc.

_C3_AXES_CUBIC = [
    np.array([1.0, 1.0, 1.0]),
    np.array([1.0, 1.0, -1.0]),
    np.array([1.0, -1.0, 1.0]),
    np.array([-1.0, 1.0, 1.0]),
]

_C4_AXES_CUBIC = [_X, _Y, _Z]

_C2_EDGE_AXES_CUBIC = [
    np.array([1.0, 1.0, 0.0]),
    np.array([1.0, -1.0, 0.0]),
    np.array([1.0, 0.0, 1.0]),
    np.array([1.0, 0.0, -1.0]),
    np.array([0.0, 1.0, 1.0]),
    np.array([0.0, 1.0, -1.0]),
]


def _build_T_rotations() -> List[np.ndarray]:
    """Proper rotation subgroup T of order 12: E, 8 C3, 3 C2."""
    ops: List[np.ndarray] = [_E.copy()]
    # 8 C3 operations: 4 axes x 2 non-identity powers
    for ax in _C3_AXES_CUBIC:
        ops.append(_cn(ax, 3, 1))
        ops.append(_cn(ax, 3, 2))
    # 3 C2 along x, y, z
    for ax in _C4_AXES_CUBIC:
        ops.append(_cn(ax, 2, 1))
    return ops  # 12 ops


def _build_Td():
    """Td of order 24: E, 8C3, 3C2, 6S4, 6sigma_d."""
    proper = _build_T_rotations()  # 12 ops

    improper: List[np.ndarray] = []
    # 6 S4 = 2 odd powers per C4 axis (S4^1, S4^3)
    for ax in _C4_AXES_CUBIC:
        improper.append(sn_power_matrix(ax, 4, 1))
        improper.append(sn_power_matrix(ax, 4, 3))
    # 6 sigma_d: diagonal mirrors with normals along edge directions
    for nrm in _C2_EDGE_AXES_CUBIC:
        improper.append(reflection_matrix(nrm))
    return proper + improper  # 24 ops


def _build_Th():
    """Th of order 24: T x { E, i } = 12 proper + 12 improper (via inversion)."""
    proper = _build_T_rotations()  # 12 ops
    inv = inversion_matrix()
    improper = [inv @ R for R in proper]
    return proper + improper  # 24 ops


def _build_O_rotations() -> List[np.ndarray]:
    """Proper rotation subgroup O of order 24: E, 8C3, 6C4, 3C2 (= C4^2), 6C2'."""
    ops: List[np.ndarray] = [_E.copy()]
    # 8 C3 (body diagonals)
    for ax in _C3_AXES_CUBIC:
        ops.append(_cn(ax, 3, 1))
        ops.append(_cn(ax, 3, 2))
    # 6 C4 (3 axes x 2 non-trivial powers: C4, C4^3)
    for ax in _C4_AXES_CUBIC:
        ops.append(_cn(ax, 4, 1))
        ops.append(_cn(ax, 4, 3))
    # 3 C2 along principal axes (= C4^2)
    for ax in _C4_AXES_CUBIC:
        ops.append(_cn(ax, 2, 1))
    # 6 C2' along edge axes
    for ax in _C2_EDGE_AXES_CUBIC:
        ops.append(_cn(ax, 2, 1))
    return ops  # 24 ops


def _build_Oh():
    """Oh of order 48: O x { E, i }."""
    proper = _build_O_rotations()  # 24 ops
    inv = inversion_matrix()
    improper = [inv @ R for R in proper]
    return proper + improper  # 48 ops


# ---------------------------------------------------------------------------
# Linear groups (finite approximation: skip continuous rotation)
# ---------------------------------------------------------------------------


def _build_Cinfv():
    """C_infinity_v discrete approximation: E + sigma_v(xz).

    The continuous Cn(z) part is omitted -- a caller working with linear
    molecules should sample rotations explicitly or treat the molecule
    as 1D.  This provides the discrete generators needed for tier-D
    template matching.
    """
    return [_E.copy(), reflection_matrix(_Y)]


def _build_Dinfh():
    """D_infinity_h discrete approximation: E, sigma_v(xz), i, sigma_h(xy)."""
    return [
        _E.copy(),
        reflection_matrix(_Y),
        inversion_matrix(),
        reflection_matrix(_Z),
    ]


# ---------------------------------------------------------------------------
# Operation labels (human-readable)
# ---------------------------------------------------------------------------


def _labels_Cn(n: int) -> List[str]:
    labels = ["E"]
    for k in range(1, n):
        if k == 1:
            labels.append(f"C{n}")
        else:
            labels.append(f"C{n}^{k}")
    return labels


def _labels_Cnv(n: int) -> List[str]:
    return _labels_Cn(n) + [f"sigma_v{i + 1}" for i in range(n)]


def _labels_Cnh(n: int) -> List[str]:
    base = _labels_Cn(n)
    return base + ["sigma_h"] + [f"Sn_{lbl}" for lbl in base[1:]]


def _labels_Dn(n: int) -> List[str]:
    return _labels_Cn(n) + [f"C2'_{i + 1}" for i in range(n)]


def _labels_Dnh(n: int) -> List[str]:
    base = _labels_Dn(n)
    return base + ["sigma_h"] + [f"sigma_v_or_Sn_{lbl}" for lbl in base[1:]]


def _labels_Dnd(n: int) -> List[str]:
    base = _labels_Dn(n)
    improper_sn = [f"S{2 * n}^{k}" for k in range(1, 2 * n, 2)]
    sigma_d = [f"sigma_d_{i + 1}" for i in range(n)]
    return base + improper_sn + sigma_d


_OP_LABELS: Dict[str, List[str]] = {
    "C1": ["E"],
    "Cs": ["E", "sigma_h"],
    "Ci": ["E", "i"],
    "C2": _labels_Cn(2),
    "C3": _labels_Cn(3),
    "C4": _labels_Cn(4),
    "C5": _labels_Cn(5),
    "C6": _labels_Cn(6),
    "C2v": _labels_Cnv(2),
    "C3v": _labels_Cnv(3),
    "C4v": _labels_Cnv(4),
    "C5v": _labels_Cnv(5),
    "C6v": _labels_Cnv(6),
    "C2h": _labels_Cnh(2),
    "C3h": _labels_Cnh(3),
    "D2": _labels_Dn(2),
    "D3": _labels_Dn(3),
    "D4": _labels_Dn(4),
    "D5": _labels_Dn(5),
    "D6": _labels_Dn(6),
    "D2h": _labels_Dnh(2),
    "D3h": _labels_Dnh(3),
    "D4h": _labels_Dnh(4),
    "D5h": _labels_Dnh(5),
    "D6h": _labels_Dnh(6),
    "D2d": _labels_Dnd(2),
    "D3d": _labels_Dnd(3),
    "D4d": _labels_Dnd(4),
    "D5d": _labels_Dnd(5),
    "Td": (
        ["E"]
        + [f"C3_{i + 1}" for i in range(8)]
        + [f"C2_{i + 1}" for i in range(3)]
        + [f"S4_{i + 1}" for i in range(6)]
        + [f"sigma_d_{i + 1}" for i in range(6)]
    ),
    "Th": (
        ["E"]
        + [f"C3_{i + 1}" for i in range(8)]
        + [f"C2_{i + 1}" for i in range(3)]
        + ["i"]
        + [f"S6_{i + 1}" for i in range(8)]
        + [f"sigma_h_{i + 1}" for i in range(3)]
    ),
    "Oh": (
        ["E"]
        + [f"C3_{i + 1}" for i in range(8)]
        + [f"C4_{i + 1}" for i in range(6)]
        + [f"C2_{i + 1}" for i in range(3)]
        + [f"C2'_{i + 1}" for i in range(6)]
        + ["i"]
        + [f"S6_{i + 1}" for i in range(8)]
        + [f"S4_{i + 1}" for i in range(6)]
        + [f"sigma_h_{i + 1}" for i in range(3)]
        + [f"sigma_d_{i + 1}" for i in range(6)]
    ),
    "Cinfv": ["E", "sigma_v"],
    "Dinfh": ["E", "sigma_v", "i", "sigma_h"],
}


# ---------------------------------------------------------------------------
# Public dispatch
# ---------------------------------------------------------------------------


_BUILDERS: Dict[str, Callable[[], List[np.ndarray]]] = {
    "C1": _build_C1,
    "Cs": _build_Cs,
    "Ci": _build_Ci,
    "C2": _build_C2,
    "C3": _build_C3,
    "C4": _build_C4,
    "C5": _build_C5,
    "C6": _build_C6,
    "C2v": _build_C2v,
    "C3v": _build_C3v,
    "C4v": _build_C4v,
    "C5v": _build_C5v,
    "C6v": _build_C6v,
    "C2h": _build_C2h,
    "C3h": _build_C3h,
    "D2": _build_D2,
    "D3": _build_D3,
    "D4": _build_D4,
    "D5": _build_D5,
    "D6": _build_D6,
    "D2h": _build_D2h,
    "D3h": _build_D3h,
    "D4h": _build_D4h,
    "D5h": _build_D5h,
    "D6h": _build_D6h,
    "D2d": _build_D2d,
    "D3d": _build_D3d,
    "D4d": _build_D4d,
    "D5d": _build_D5d,
    "Td": _build_Td,
    "Th": _build_Th,
    "Oh": _build_Oh,
    "Cinfv": _build_Cinfv,
    "Dinfh": _build_Dinfh,
}


_PG_OPS_CACHE: Dict[str, List[np.ndarray]] = {}


def get_operations(pg_name: str) -> List[np.ndarray]:
    """Return all symmetry operations of point group as list of 3x3 matrices.

    First element is always identity (E).  Operations include proper
    rotations (det = +1) and improper operations -- reflections,
    inversion, Sn (det = -1).

    Parameters
    ----------
    pg_name : str
        Point-group Schoenflies symbol, e.g. ``"C3v"``, ``"D4h"``,
        ``"Oh"``, ``"Td"``.  Linear groups use ``"Cinfv"`` and
        ``"Dinfh"`` (discrete approximation).

    Raises
    ------
    KeyError
        If ``pg_name`` is not supported.
    """
    if pg_name in _PG_OPS_CACHE:
        return _PG_OPS_CACHE[pg_name]
    if pg_name not in _BUILDERS:
        raise KeyError(f"Point group '{pg_name}' not implemented")
    ops = _BUILDERS[pg_name]()
    _PG_OPS_CACHE[pg_name] = ops
    return ops


def get_operation_labels(pg_name: str) -> List[str]:
    """Return human-readable labels for each operation of the point group."""
    if pg_name not in _OP_LABELS:
        raise KeyError(f"Point group '{pg_name}' has no label table")
    return list(_OP_LABELS[pg_name])


def supported_point_groups() -> List[str]:
    """List of all supported point-group names."""
    return sorted(_BUILDERS.keys())


# ---------------------------------------------------------------------------
# Verification helpers
# ---------------------------------------------------------------------------


# Expected group orders for sanity check.
_EXPECTED_ORDERS: Dict[str, int] = {
    "C1": 1, "Cs": 2, "Ci": 2,
    "C2": 2, "C3": 3, "C4": 4, "C5": 5, "C6": 6,
    "C2v": 4, "C3v": 6, "C4v": 8, "C5v": 10, "C6v": 12,
    "C2h": 4, "C3h": 6,
    "D2": 4, "D3": 6, "D4": 8, "D5": 10, "D6": 12,
    "D2h": 8, "D3h": 12, "D4h": 16, "D5h": 20, "D6h": 24,
    "D2d": 8, "D3d": 12, "D4d": 16, "D5d": 20,
    "Td": 24, "Th": 24, "Oh": 48,
    "Cinfv": 2, "Dinfh": 4,
}


def _verify_pg(pg_name: str, tol: float = 1e-9) -> None:
    """Verify operations of point group satisfy basic group properties.

    Checks:
      1. Expected operation count.
      2. First operation is identity.
      3. det = +/- 1 for every operation.
      4. Operations are orthogonal: R.T @ R = I.
      5. Identity is unique (no duplicates of E).
    """
    ops = get_operations(pg_name)

    expected = _EXPECTED_ORDERS.get(pg_name)
    if expected is not None and len(ops) != expected:
        raise AssertionError(
            f"{pg_name}: expected {expected} operations, got {len(ops)}"
        )

    if not np.allclose(ops[0], _E, atol=tol):
        raise AssertionError(f"{pg_name}: first op is not identity")

    for idx, op in enumerate(ops):
        d = float(np.linalg.det(op))
        if abs(abs(d) - 1.0) > 1e-9:
            raise AssertionError(
                f"{pg_name}: op[{idx}] determinant {d} not in {{+1, -1}}"
            )
        if not np.allclose(op.T @ op, _E, atol=1e-9):
            raise AssertionError(f"{pg_name}: op[{idx}] not orthogonal")

    # Check no duplicate ops (within tolerance) -- spot-check identity only
    n_identity = sum(1 for op in ops if np.allclose(op, _E, atol=1e-9))
    if n_identity != 1:
        raise AssertionError(
            f"{pg_name}: identity appears {n_identity} times (expected 1)"
        )


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import sys

    failures: List[str] = []
    for pg in supported_point_groups():
        try:
            _verify_pg(pg)
            ops = get_operations(pg)
            n_proper = sum(1 for op in ops if np.linalg.det(op) > 0)
            n_improper = sum(1 for op in ops if np.linalg.det(op) < 0)
            print(
                f"  {pg:8s}  |G| = {len(ops):3d}   "
                f"proper = {n_proper:3d}  improper = {n_improper:3d}   OK"
            )
        except Exception as exc:  # pragma: no cover -- self-test
            failures.append(f"{pg}: {exc}")
            print(f"  {pg:8s}  FAIL: {exc}")

    # Spot checks documented in the masterplan
    spot_checks = {"C3v": 6, "D4h": 16, "Td": 24, "Th": 24, "Oh": 48}
    for pg, expected in spot_checks.items():
        ops = get_operations(pg)
        assert len(ops) == expected, f"{pg}: |G| = {len(ops)} != {expected}"

    if failures:
        print(f"\n{len(failures)} failures:")
        for msg in failures:
            print(f"  {msg}")
        sys.exit(1)

    print(f"\nAll {len(supported_point_groups())} point groups verified.")
