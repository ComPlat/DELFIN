"""delfin.fffree.conformer_symmetry — point-group order detection and
symmetry-aware cluster-representative scoring for the RMSD-dedup pipeline.

User-Insight 2026-06-02 (project_grip_symmetry_priority_rotamers_2026_06_02):

> "rotationen die höhere symmetrie erzeugen sind für uns immer sehr
>  interessant. stichwort group theory"

High-point-group-symmetry conformers (canonical XRD forms, degenerate-
orbital stabilised, cleaner spectroscopy) currently tie with C1 forms in
the RMSD-Butina cluster-representative selection.  This module adds a
universal, geometry-only point-group detector + a tunable log2-symmetry-
weighted score that the :mod:`delfin.fffree.conformer_dedup` pipeline can
opt into via the new env-flag ``DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY``
(default OFF -> byte-identical to HEAD).

Determinism:
  - Lex-sorted candidate-axis generation (atom indices ascending, pair
    indices in (i, j) with i < j).
  - No randomness, no hash-set ordering: every dict is rebuilt as a list
    in insertion order, and PYTHONHASHSEED=0 is honoured by the codebase.
  - Identical inputs -> bit-identical output.

Universal (graph/geometry only):
  - Inputs are Cartesian coordinates + atomic numbers; no SMILES, no
    metal-specific heuristics.  Works on transition-metal complexes,
    pure organics, clusters.

Backend:
  - Prefers spglib (>= 1.16, available as ``spglib==2.7.0`` in env
    ``delfin``) for crystallographic-grade detection on cell-aligned
    geometries (octahedron at xyz axes -> 48 ops).
  - For free-floating molecules whose principal axes do not align with
    the Cartesian axes (benzene at random orientation, NH3, ...), spglib
    under-counts because the vacuum-cell wrapper only sees the orthogonal-
    cell subgroup.  The custom backend is a robust fallback: candidate
    symmetry axes from atomic positions + inertia eigenvectors, then
    rotation/reflection/inversion testing with composition closure.
  - We take the MAXIMUM of the two backends.

Public API:
  - :func:`detect_point_group_order(coords, atomic_numbers, tol=0.1)`
  - :func:`score_conformer_for_representative(frame, w_sym=0.5, w_sev=1.0)`
"""
from __future__ import annotations

import os
from itertools import combinations
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np


# ----------------------------------------------------------------------
# Defaults & env helpers
# ----------------------------------------------------------------------

DEFAULT_TOL = 0.1                           # Å — position-equivalence tolerance
DEFAULT_W_SYM = 0.5                         # log2-symmetry weight
DEFAULT_W_SEV = 1.0                         # severity weight
ENV_FLAG_PRIORITY = "DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY"
ENV_W_SYM = "DELFIN_FFFREE_RMSD_DEDUP_SYM_WEIGHT"


def _env_priority_enabled(default: bool = False) -> bool:
    """``True`` when the symmetry-priority representative-selection hook is
    active.  Default OFF -> byte-identical to HEAD.
    """
    v = os.environ.get(ENV_FLAG_PRIORITY, "1" if default else "0")
    return str(v).strip() in ("1", "true", "True", "on", "yes", "YES")


def _env_w_sym(default: float = DEFAULT_W_SYM) -> float:
    try:
        return float(os.environ.get(ENV_W_SYM, str(default)))
    except (TypeError, ValueError):
        return default


# ----------------------------------------------------------------------
# Element-symbol -> atomic-number table (small subset; sufficient for the
# tests + the chemistry the dedup pipeline sees).  Used to accept
# element-symbol input transparently.
# ----------------------------------------------------------------------

_PERIODIC_TABLE = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn "
    "Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
    "In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu "
    "Am Cm Bk Cf Es Fm Md No Lr"
).split()
_SYM_TO_Z = {s: i for i, s in enumerate(_PERIODIC_TABLE)}


def _coerce_atomic_numbers(atomic_numbers: Sequence) -> List[int]:
    """Accept either integer atomic numbers (``[6, 1, 1, 1, 1]``) or element
    symbols (``["C", "H", "H", "H", "H"]``).  Returns a list of ints.

    Unknown element symbols map to 0 (the "X" dummy) which still works for
    equivalence-by-label.
    """
    out: List[int] = []
    for a in atomic_numbers:
        if isinstance(a, (int, np.integer)):
            out.append(int(a))
        else:
            s = str(a).strip()
            out.append(int(_SYM_TO_Z.get(s, 0)))
    return out


# ----------------------------------------------------------------------
# Symmetry-element primitives
# ----------------------------------------------------------------------


def _rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """Rodrigues' rotation matrix for ``theta`` radians around the unit
    vector ``axis``.
    """
    a = np.asarray(axis, dtype=float)
    a = a / max(float(np.linalg.norm(a)), 1e-12)
    c, s = float(np.cos(theta)), float(np.sin(theta))
    C = 1.0 - c
    x, y, z = float(a[0]), float(a[1]), float(a[2])
    return np.array([
        [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
        [x * y * C + z * s, c + y * y * C,     y * z * C - x * s],
        [x * z * C - y * s, y * z * C + x * s, c + z * z * C   ],
    ], dtype=float)


def _reflection_matrix(normal: np.ndarray) -> np.ndarray:
    """Householder reflection through the plane with unit ``normal``."""
    n = np.asarray(normal, dtype=float)
    n = n / max(float(np.linalg.norm(n)), 1e-12)
    return np.eye(3) - 2.0 * np.outer(n, n)


def _maps_to_same(P: np.ndarray, P_t: np.ndarray, Z: Sequence[int],
                  tol: float) -> bool:
    """Greedy nearest-neighbour matching: for every atom ``i`` in the
    transformed coordinates ``P_t``, find an unused atom ``j`` of the same
    atomic number in ``P`` with position deviation below ``tol``.

    Returns True iff a complete bijection exists.  Deterministic (scans
    j in ascending order).
    """
    n = P.shape[0]
    if P_t.shape != P.shape:
        return False
    used = [False] * n
    for i in range(n):
        zi = int(Z[i])
        best_j = -1
        best_d = float(tol)
        for j in range(n):
            if used[j] or int(Z[j]) != zi:
                continue
            d = float(np.linalg.norm(P_t[i] - P[j]))
            if d < best_d:
                best_d = d
                best_j = j
        if best_j < 0:
            return False
        used[best_j] = True
    return True


# ----------------------------------------------------------------------
# Candidate-axis enumeration (lex-sorted for determinism)
# ----------------------------------------------------------------------


def _dedup_axes(axes: Iterable[np.ndarray], precision: int = 2
                ) -> List[np.ndarray]:
    """Deduplicate candidate axis directions by rounded unit-vector key,
    treating ``+a`` and ``-a`` as identical.  Insertion-order preserved.

    Default ``precision = 2`` decimals (~0.6 degrees) merges axes that
    differ only by numerical noise from imperfect geometries -- important
    because slightly off-axis candidates lead the closure step to inflate
    the discovered group with near-redundant rotations.
    """
    seen: List[Tuple[int, ...]] = []
    out: List[np.ndarray] = []
    for a in axes:
        a = np.asarray(a, dtype=float)
        norm = float(np.linalg.norm(a))
        if norm < 1e-9:
            continue
        u = a / norm
        key = tuple(int(round(v * 10 ** precision)) for v in u)
        nkey = tuple(-k for k in key)
        if key in seen or nkey in seen:
            continue
        seen.append(key)
        out.append(u)
    return out


def _candidate_axes(P: np.ndarray, Z: Sequence[int],
                    tol: float) -> List[np.ndarray]:
    """Enumerate candidate symmetry-element axes / mirror-plane normals.

    Sources (lex-ordered):
      1. Position vectors of atoms ``i`` with ``|P[i]| > tol``.
      2. For every equivalent atom pair ``(i, j)``, ``i < j``, same Z,
         similar ``|P|``:
           - midpoint direction (axis through midpoint)
           - cross product (axis perpendicular to the pair plane)
           - difference vector (axis along the bond direction)
      3. Principal moments-of-inertia eigenvectors.
      4. Cartesian axes (x, y, z) as a permanent fallback.

    Output deduplicated up to sign and 3-decimal precision.
    """
    n = P.shape[0]
    axes: List[np.ndarray] = []

    # 1) Through each atom
    for i in range(n):
        r = float(np.linalg.norm(P[i]))
        if r > tol:
            axes.append(P[i] / r)

    # 2) Through midpoints / cross / differences of equivalent pairs
    for i, j in combinations(range(n), 2):
        if int(Z[i]) != int(Z[j]):
            continue
        ri = float(np.linalg.norm(P[i]))
        rj = float(np.linalg.norm(P[j]))
        if abs(ri - rj) > max(tol, 0.5):
            continue
        mid = (P[i] + P[j]) * 0.5
        m = float(np.linalg.norm(mid))
        if m > tol:
            axes.append(mid / m)
        cr = np.cross(P[i], P[j])
        cm = float(np.linalg.norm(cr))
        if cm > tol:
            axes.append(cr / cm)
        diff = P[i] - P[j]
        dm = float(np.linalg.norm(diff))
        if dm > tol:
            axes.append(diff / dm)

    # 3) Principal moments-of-inertia eigenvectors
    if n >= 2:
        I = np.zeros((3, 3), dtype=float)
        for i in range(n):
            x, y, z = float(P[i, 0]), float(P[i, 1]), float(P[i, 2])
            I[0, 0] += y * y + z * z
            I[1, 1] += x * x + z * z
            I[2, 2] += x * x + y * y
            I[0, 1] -= x * y
            I[0, 2] -= x * z
            I[1, 2] -= y * z
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        try:
            _, V = np.linalg.eigh(I)
            for k in range(3):
                axes.append(V[:, k])
        except np.linalg.LinAlgError:
            pass

    # 4) Cartesian fallback
    axes.append(np.array([1.0, 0.0, 0.0]))
    axes.append(np.array([0.0, 1.0, 0.0]))
    axes.append(np.array([0.0, 0.0, 1.0]))

    return _dedup_axes(axes)


# ----------------------------------------------------------------------
# Backends
# ----------------------------------------------------------------------


def _custom_point_group_order(P: np.ndarray, Z: Sequence[int],
                              tol: float = DEFAULT_TOL) -> int:
    """Custom geometric point-group order detector.

    Algorithm:
      1. Center the molecule (subtract centroid).
      2. Enumerate candidate axes via :func:`_candidate_axes`.
      3. Test inversion ``i``.
      4. For each candidate axis ``ax`` and each ``n in (2, 3, 4, 5, 6)``:
         - test Cn rotation (and add all powers if it holds).
         - test mirror plane with normal ``ax``.
         - test improper rotations Sn for ``n in (3, 4, 6)``.
      5. Close the operation set under matrix composition (up to 5 iters).
      6. Return the size of the closed group.

    Returns at least 1 (identity).  Deterministic: candidate axes are
    sorted by the dedup key, compositions iterate in insertion order.
    """
    n = int(P.shape[0])
    if n == 0:
        return 1
    Pc = P - P.mean(axis=0)
    Z_list = list(int(z) for z in Z)
    axes = _candidate_axes(Pc, Z_list, tol)

    ops: List[np.ndarray] = [np.eye(3)]

    # Use a coarse rounding key for O(1) dedup.  Two decimals + tolerance
    # zero-out is enough to merge numerically-drifted duplicates because
    # rotation/reflection matrix entries live in [-1, 1] and any unique
    # symmetry op of a point group differs from another in at least one
    # entry by O(0.1).
    def _key(M: np.ndarray) -> bytes:
        R = np.round(M, 2)
        R[np.abs(R) < 1e-3] = 0.0
        return R.tobytes()

    seen_keys = {_key(np.eye(3))}

    def has_op(op: np.ndarray) -> bool:
        return _key(op) in seen_keys

    def is_sym(op: np.ndarray) -> bool:
        return _maps_to_same(Pc, Pc @ op.T, Z_list, tol)

    def try_add(op: np.ndarray) -> bool:
        if has_op(op):
            return False
        if is_sym(op):
            ops.append(op)
            seen_keys.add(_key(op))
            return True
        return False

    # Inversion
    try_add(-np.eye(3))

    # Per-axis Cn / sigma / Sn
    for ax in axes:
        # Cn proper rotations
        for n_fold in (2, 3, 4, 5, 6):
            theta = 2.0 * np.pi / float(n_fold)
            R = _rotation_matrix(ax, theta)
            if is_sym(R):
                # add all powers C_n^k (k=1..n-1) -- group-theoretically
                # implied, no further is_sym check needed.
                Rk = np.eye(3)
                for _ in range(1, n_fold):
                    Rk = Rk @ R
                    k = _key(Rk)
                    if k not in seen_keys:
                        ops.append(Rk)
                        seen_keys.add(k)
        # Mirror plane with normal = ax
        S = _reflection_matrix(ax)
        try_add(S)
        # Improper rotations Sn = sigma_h(ax) @ Cn(ax)
        for n_fold in (3, 4, 6):
            R = _rotation_matrix(ax, 2.0 * np.pi / float(n_fold))
            sh = _reflection_matrix(ax)
            Sn = sh @ R
            if is_sym(Sn):
                # Cyclic powers -- group-theoretic, no further is_sym needed.
                M = np.eye(3)
                for _ in range(1, 2 * n_fold):
                    M = M @ Sn
                    k = _key(M)
                    if k not in seen_keys:
                        ops.append(M)
                        seen_keys.add(k)

    # Close under composition.  Group-theoretic fact: if a and b are
    # symmetries of the molecule then a @ b is too -- no further is_sym
    # verification needed.  We just need O(1) dedup on the running set.
    # ``seen_keys`` (initialised above) is reused as the rounded-bytes
    # membership set.  Bounded iterations to keep determinism predictable;
    # in practice the closed group is reached in 2-3 passes for typical
    # molecules.
    for _ in range(8):
        before = len(ops)
        # Snapshot the current op list to iterate over a stable index range.
        current = list(ops)
        for a in current:
            for b in current:
                c = a @ b
                k = _key(c)
                if k in seen_keys:
                    continue
                # Re-verify against the molecule to catch any product whose
                # numerical drift has pushed it OUT of the true symmetry
                # group (defensive; in exact arithmetic closure is automatic).
                if not is_sym(c):
                    seen_keys.add(k)  # remember the false-positive key
                    continue
                seen_keys.add(k)
                ops.append(c)
        if len(ops) == before:
            break

    return int(len(ops))


def _try_spglib_order(P: np.ndarray, Z: Sequence[int],
                      tol: float = DEFAULT_TOL) -> Optional[int]:
    """Try spglib (vacuum-cell wrapper).  Returns ``None`` if spglib is
    unavailable, the dataset is empty, or the call raises.

    Note: the vacuum-cell trick only catches symmetries whose rotations
    are compatible with the orthogonal-cell axes (e.g. octahedron aligned
    to xyz).  A molecule whose principal axes are tilted is under-counted
    by spglib; the caller takes ``max(spglib, custom)`` to recover.
    """
    try:
        import spglib  # type: ignore
    except Exception:
        return None
    try:
        if P.size == 0:
            return None
        lo = P.min(axis=0) - 10.0
        hi = P.max(axis=0) + 10.0
        cell = np.diag(hi - lo)
        frac = (P - lo) / np.maximum(hi - lo, 1e-9)
        # Map atomic numbers directly (spglib only cares about distinct
        # integer labels; element identity is unused).
        ids = [int(z) for z in Z]
        ds = spglib.get_symmetry_dataset((cell, frac, ids), symprec=float(tol))
        if ds is None:
            return None
        # spglib's dict interface is being deprecated; support both.
        try:
            rotations = ds["rotations"]
        except Exception:
            rotations = getattr(ds, "rotations", None)
        if rotations is None:
            return None
        return max(1, int(len(rotations)))
    except Exception:
        return None


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def detect_point_group_order(coords: np.ndarray,
                             atomic_numbers: Sequence,
                             tol: float = DEFAULT_TOL) -> int:
    """Return the point-group ORDER of ``(coords, atomic_numbers)``.

    Parameters
    ----------
    coords : (N, 3) array-like
        Cartesian coordinates.
    atomic_numbers : sequence of int OR element symbols
        Atomic numbers (``[6, 1, 1, 1, 1]``) or element symbols
        (``["C", "H", "H", "H", "H"]``); both forms accepted.
    tol : float, default 0.1
        Position-deviation tolerance (Å) for equivalence under a
        symmetry operation.

    Returns
    -------
    int
        Group order >= 1.  Typical values: 1 (C1), 2 (Cs/Ci/C2),
        3 (C3), 4 (C2v/C4), 6 (C3v/D3), 8 (D4/D2h), 12 (D6/T),
        24 (D6h/Td/O), 48 (Oh).

    Backend choice: ``max(spglib, custom)``.  spglib is tried first
    (deterministic, fast for cell-aligned structures); the custom
    candidate-axis algorithm covers tilted molecules where spglib
    under-counts.

    Universal: depends only on coordinates and atomic numbers; no
    SMILES, no metal lookup, no bond graph.  Deterministic.
    """
    P = np.asarray(coords, dtype=float)
    if P.ndim != 2 or P.shape[1] != 3:
        return 1
    n = int(P.shape[0])
    if n == 0:
        return 1
    Z = _coerce_atomic_numbers(atomic_numbers)
    if len(Z) != n:
        return 1
    # Center for the custom backend; spglib accepts any centering.
    Pc = P - P.mean(axis=0)
    sp = _try_spglib_order(Pc, Z, tol=tol)
    cu = _custom_point_group_order(Pc, Z, tol=tol)
    if sp is None:
        return int(cu)
    return int(max(sp, cu))


def score_conformer_for_representative(
    frame,
    w_sym: float = DEFAULT_W_SYM,
    w_sev: float = DEFAULT_W_SEV,
    tol: float = DEFAULT_TOL,
) -> float:
    """Combined symmetry / severity score for cluster-representative pick.

    Score formula (higher is better):

        score = w_sym * log2(pg_order) - w_sev * severity

    Rationale:

      - ``log2(pg_order)`` rewards group-order increases in a sub-linear
        way: C1 -> 0, C2/Cs/Ci -> 1, C3 -> ~1.58, C4 -> 2, D3 -> ~2.58,
        D4 -> 3, D6 -> ~3.58, Td/O -> ~4.58, Oh -> ~5.58.  Doubling the
        order adds one full point of bonus, keeping the trade-off with
        severity readable.
      - ``-w_sev * severity`` is the Mogul severity penalty (lower is
        better).  With the defaults ``w_sym=0.5, w_sev=1.0``, the
        symmetry bonus tops out at ~+2.8 for Oh, comparable to a one-
        unit severity drop.

    Parameters
    ----------
    frame : dict or object with attributes
        Must expose ``coords`` and ``atomic_numbers``; ``mogul_severity``
        (or ``severity``) optional, defaults to 0.0.  Accepts a dict
        with keys ``coords`` / ``atomic_numbers`` / ``mogul_severity``,
        an instance of :class:`delfin.fffree.conformer_dedup.Frame`
        (uses ``label``, ``P``, ``severity``, plus ``extra["atomic_numbers"]``
        or ``extra["symbols"]``), or any other dataclass-like with the
        expected fields.
    w_sym : float, default 0.5
        Weight on the ``log2(pg_order)`` term.
    w_sev : float, default 1.0
        Weight on the severity term.
    tol : float, default 0.1
        Position tolerance for the point-group detector.
    """
    coords = None
    atomic_numbers: Optional[Sequence] = None
    severity = 0.0

    if isinstance(frame, dict):
        coords = frame.get("coords")
        atomic_numbers = (frame.get("atomic_numbers")
                          or frame.get("symbols")
                          or frame.get("elements"))
        severity = float(frame.get("mogul_severity",
                                   frame.get("severity", 0.0)) or 0.0)
    else:
        # Object-style.  Support the Frame dataclass shape used by
        # conformer_dedup: ``P`` array + ``extra`` dict carrying symbols /
        # atomic_numbers.
        coords = getattr(frame, "coords", None)
        if coords is None:
            coords = getattr(frame, "P", None)
        atomic_numbers = getattr(frame, "atomic_numbers", None)
        if atomic_numbers is None:
            atomic_numbers = getattr(frame, "symbols", None)
        extra = getattr(frame, "extra", None)
        if atomic_numbers is None and isinstance(extra, dict):
            atomic_numbers = (extra.get("atomic_numbers")
                              or extra.get("symbols")
                              or extra.get("elements"))
        sev = getattr(frame, "mogul_severity", None)
        if sev is None:
            sev = getattr(frame, "severity", 0.0)
        try:
            severity = float(sev) if sev is not None else 0.0
        except (TypeError, ValueError):
            severity = 0.0

    if coords is None or atomic_numbers is None:
        # Without geometry, score is purely severity-driven.
        return float(-float(w_sev) * severity)
    try:
        pg = detect_point_group_order(
            np.asarray(coords, dtype=float),
            list(atomic_numbers),
            tol=tol,
        )
    except Exception:
        pg = 1
    pg = max(1, int(pg))
    sym_bonus = float(w_sym) * float(np.log2(pg))
    return float(sym_bonus - float(w_sev) * float(severity))


# ----------------------------------------------------------------------
# Cluster-representative selector (used by conformer_dedup when the env
# flag is on).  Kept here so the dedup module only has to import a
# single function.
# ----------------------------------------------------------------------


def select_symmetric_priority_representative(
    cluster_indices: Sequence[int],
    coords_list: Sequence[np.ndarray],
    severities: Sequence[float],
    atomic_numbers: Sequence,
    w_sym: Optional[float] = None,
    w_sev: float = DEFAULT_W_SEV,
    tol: float = DEFAULT_TOL,
) -> int:
    """Pick the index in ``cluster_indices`` with the highest combined
    log2-symmetry/severity score (:func:`score_conformer_for_representative`).

    Ties broken by lower severity, then by lower input index.  All-args
    interface so the caller can re-use the dedup module's normalised
    coords/severities arrays across clusters without re-instantiating
    Frame objects.
    """
    if w_sym is None:
        w_sym = _env_w_sym()
    best_idx = int(cluster_indices[0])
    best_score = -float("inf")
    Z = _coerce_atomic_numbers(atomic_numbers)
    for idx in cluster_indices:
        i = int(idx)
        try:
            pg = detect_point_group_order(coords_list[i], Z, tol=tol)
        except Exception:
            pg = 1
        sev = float(severities[i])
        score = float(w_sym) * float(np.log2(max(1, pg))) - float(w_sev) * sev
        if score > best_score:
            best_score = score
            best_idx = i
        elif abs(score - best_score) < 1e-12:
            if (sev, i) < (float(severities[best_idx]), int(best_idx)):
                best_idx = i
    return int(best_idx)
