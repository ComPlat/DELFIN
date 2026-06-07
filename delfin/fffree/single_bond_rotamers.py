"""delfin.fffree.single_bond_rotamers — deterministic single-bond rotamer
enumeration (gauche+/anti/gauche-) for fffree assembled complexes.

For each rotatable single bond (non-ring, non-metal-incident, non-double),
enumerate N dihedral states relative to the as-built geometry, build the
Cartesian product of states (capped at ``max_configs``), and yield each
distinct config as a new coordinate array.

Universal: detection is by graph topology (bond order + ring + metal flag),
not by SMILES patterns.  Deterministic: bonds sorted by atom index pair;
rotamer order fixed (identity first, then evenly spaced offsets); combinations
enumerated in ``itertools.product`` order.

Rotamer-step tightening (subagent #129 follow-up to #128):
  - default step **±60°** -> 6 states per bond (0, +60, -60, +120, -120, 180)
  - env ``DELFIN_FFFREE_ROTAMER_STEP_DEGREES`` (default 60.0) sets the step
  - default ``MAX_CONFIGS`` is **144** (cap 12² rotamer-pairs at 12 states each
    when using a 30° finer-grid step; the user-spec default 60° step yields
    6 states -> 6² = 36 cases for n=2 well under 144).

Env-gate: ``DELFIN_FFFREE_ENUMERATE_ROTAMERS=1`` (default OFF -> the
module is byte-identical to no integration when callers honor the
flag at the call-site).

Symmetry-priority rotamer pre-ranking (Phase 2, 2026-06-03):
  - ``DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS=1`` (default OFF) opt-in
    flag that re-sorts the yielded rotamer variants so configurations
    that maximise the whole-complex point-group order come FIRST.
    Implemented via :func:`rerank_rotamers_by_symmetry`, callable
    standalone for tests and used internally by
    :func:`enumerate_single_bond_rotamers` when the flag is set.
"""
from __future__ import annotations

import itertools
import os
from typing import Iterator, List, Optional, Sequence, Tuple

import numpy as np


# ------------------------------------------------------------------
# Defaults — overridable via env
# ------------------------------------------------------------------

#: Number of rotamer states per rotatable bond.  3 = gauche+/anti/gauche-.
#: When using ``_DEFAULT_STEP_DEGREES`` this is derived from the step
#: (360 / step).  Kept for callers that pass ``n_states`` explicitly.
_DEFAULT_N_STATES = 3

#: Maximum number of rotatable bonds enumerated per molecule.  Beyond this
#: the bond list is truncated (lowest atom-index pair wins) to keep the
#: combinatorial product manageable.
_DEFAULT_MAX_ROTORS = 5

#: Hard cap on the number of enumerated configurations.  144 = 12² is
#: tuned for the new finer step (#129 follow-up): with ``_DEFAULT_MAX_ROTORS``
#: fixed at 5 the cap is hit at moderate state counts and the product is
#: truncated to keep RMSD-dedup tractable.  Was 243 (= 3⁵) under the
#: previous coarse ±120° step.
_DEFAULT_MAX_CONFIGS = 144

#: Default dihedral step in degrees for the new finer-grid enumeration.
#: With ``_DEFAULT_STEP_DEGREES = 60`` the per-rotor state count is
#: ``int(round(360 / 60)) = 6`` -> 6, i.e. (0, ±60, ±120, 180).
#: Override per call via ``step_degrees`` or globally via
#: ``DELFIN_FFFREE_ROTAMER_STEP_DEGREES``.
_DEFAULT_STEP_DEGREES = 60.0


# ------------------------------------------------------------------
# Helpers — env-tunable wrappers
# ------------------------------------------------------------------


def _max_rotors() -> int:
    try:
        return max(0, int(os.environ.get(
            "DELFIN_FFFREE_ENUMERATE_ROTAMERS_MAX_ROTORS",
            str(_DEFAULT_MAX_ROTORS),
        )))
    except (TypeError, ValueError):
        return _DEFAULT_MAX_ROTORS


def _max_configs() -> int:
    try:
        return max(1, int(os.environ.get(
            "DELFIN_FFFREE_ENUMERATE_ROTAMERS_MAX_CONFIGS",
            str(_DEFAULT_MAX_CONFIGS),
        )))
    except (TypeError, ValueError):
        return _DEFAULT_MAX_CONFIGS


def _n_states() -> int:
    # When STATES is set explicitly, honor it as-is.
    raw_states = os.environ.get("DELFIN_FFFREE_ENUMERATE_ROTAMERS_STATES")
    if raw_states is not None:
        try:
            return max(1, int(raw_states))
        except (TypeError, ValueError):
            pass
    # Else derive from step-degrees (subagent #129 finer-grid default).
    step = _step_degrees()
    if step <= 0.0:
        return _DEFAULT_N_STATES
    return max(1, int(round(360.0 / step)))


def _step_degrees() -> float:
    """Per-bond dihedral step in degrees (subagent #129 finer-grid knob).

    Default = 60° -> 6 evenly-spaced states per rotor (identity + ±60° + ±120°
    + 180°).  Smaller -> denser sampling, larger product.  When the env-var
    is malformed we fall back to ``_DEFAULT_STEP_DEGREES``.
    """
    try:
        v = float(os.environ.get(
            "DELFIN_FFFREE_ROTAMER_STEP_DEGREES",
            str(_DEFAULT_STEP_DEGREES),
        ))
    except (TypeError, ValueError):
        return _DEFAULT_STEP_DEGREES
    return v if v > 0.0 else _DEFAULT_STEP_DEGREES


# ------------------------------------------------------------------
# Rotatable-bond detection
# ------------------------------------------------------------------


def _is_metal_atom(atom) -> bool:
    """Detect a metal atom by atomic number (universal, ignores formal charge).

    Mirrors delfin._bond_decollapse._is_metal but operates on an RDKit Atom.
    Includes alkali / alkaline-earth, all d-block, all f-block, and the
    p-block metals (Al/Ga/In/Tl/Sn/Pb/Bi/Po).
    """
    try:
        z = int(atom.GetAtomicNum())
    except Exception:
        return False
    # Alkali (3, 11, 19, 37, 55, 87)
    # Alk. earth (4, 12, 20, 38, 56, 88)
    # d-block (21-30, 39-48, 72-80, 104-112)
    # f-block (57-71, 89-103)
    # p-block metals (13 Al, 31 Ga, 49 In, 81 Tl, 50 Sn, 82 Pb, 83 Bi, 84 Po)
    if z in (3, 11, 19, 37, 55, 87):
        return True
    if z in (4, 12, 20, 38, 56, 88):
        return True
    if 21 <= z <= 30 or 39 <= z <= 48 or 72 <= z <= 80 or 104 <= z <= 112:
        return True
    if 57 <= z <= 71 or 89 <= z <= 103:
        return True
    if z in (13, 31, 49, 81, 50, 82, 83, 84):
        return True
    return False


def find_rotatable_bonds(mol) -> List[Tuple[int, int]]:
    """Return the sorted list of ``(a, b)`` index pairs corresponding to
    rotatable single bonds.

    A bond is "rotatable" when it satisfies ALL of the following:

    - bond order == SINGLE
    - not aromatic
    - not in a ring
    - both endpoints are heavy (non-H)
    - neither endpoint is a metal
    - both endpoints have at least one non-H neighbour besides each other
      (terminal -CH3 / -OH spins yield no distinct conformer)

    Returns the pair as ``(min(a, b), max(a, b))`` and the list is sorted
    by ``(a, b)`` for deterministic enumeration order.

    Universal: depends only on the molecular graph; no SMILES heuristics.
    """
    if mol is None:
        return []
    out: List[Tuple[int, int]] = []
    try:
        bonds = list(mol.GetBonds())
    except Exception:
        return []
    for b in bonds:
        try:
            from rdkit import Chem
            if b.GetBondType() != Chem.BondType.SINGLE:
                continue
            if b.GetIsAromatic() or b.IsInRing():
                continue
        except Exception:
            continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        if _is_metal_atom(a1) or _is_metal_atom(a2):
            continue
        # Need at least one heavy neighbour on EACH side (besides the partner)
        h1 = sum(
            1 for n in a1.GetNeighbors()
            if n.GetIdx() != a2.GetIdx() and n.GetAtomicNum() > 1
        )
        h2 = sum(
            1 for n in a2.GetNeighbors()
            if n.GetIdx() != a1.GetIdx() and n.GetAtomicNum() > 1
        )
        if h1 < 1 or h2 < 1:
            continue
        i, j = int(a1.GetIdx()), int(a2.GetIdx())
        out.append((min(i, j), max(i, j)))
    out.sort()
    return out


def _dihedral_reference_atoms(mol, b1: int, b2: int) -> Optional[Tuple[int, int]]:
    """Pick deterministic dihedral-reference atoms ``a0`` and ``a3``
    so that the dihedral ``a0-b1-b2-a3`` defines the torsion.

    Selection is the lowest-index heavy neighbour of ``b1`` other than
    ``b2`` for ``a0``, and the same for ``b2 -> a3``.  Returns None if
    either side has no qualifying neighbour (should not happen post-
    :func:`find_rotatable_bonds` filter).
    """
    try:
        a1 = sorted(
            n.GetIdx() for n in mol.GetAtomWithIdx(b1).GetNeighbors()
            if n.GetIdx() != b2 and n.GetAtomicNum() > 1
        )
        a3 = sorted(
            n.GetIdx() for n in mol.GetAtomWithIdx(b2).GetNeighbors()
            if n.GetIdx() != b1 and n.GetAtomicNum() > 1
        )
    except Exception:
        return None
    if not a1 or not a3:
        return None
    return int(a1[0]), int(a3[0])


# ------------------------------------------------------------------
# Enumeration of dihedral configurations
# ------------------------------------------------------------------


def enumerate_rotamer_configs(
    n_rotors: int,
    n_states: Optional[int] = None,
    max_configs: Optional[int] = None,
    step_degrees: Optional[float] = None,
) -> Iterator[Tuple[float, ...]]:
    """Yield each dihedral-OFFSET configuration as a tuple of length
    ``n_rotors``.  Offsets are uniformly spaced over 360° starting at 0°.

    For ``n_states = 3``: offsets are (0, +120°, -120°) -> three rotamers
    relative to the existing geometry (the "0" offset is the as-built
    conformer; the +/-120 offsets are the two gauche partners).

    For the new ``step_degrees=60`` default (subagent #129 follow-up): six
    evenly-spaced offsets per rotor (0, ±60, ±120, 180), producing 6**n_rotors
    combinations up to ``max_configs``.

    The first yielded tuple is always the all-zero offset = identity, so
    callers always get the unmodified molecule as the first variant.

    Deterministic: ``itertools.product`` over a fixed offset sequence.
    """
    mc = int(max_configs) if max_configs is not None else _max_configs()
    if n_rotors <= 0:
        yield ()
        return

    # Resolve the offset sequence -- explicit ``step_degrees`` wins, then
    # explicit ``n_states``, then env-derived defaults.
    if step_degrees is not None and float(step_degrees) > 0.0:
        sd = float(step_degrees)
        ns = max(1, int(round(360.0 / sd)))
    elif n_states is not None:
        ns = int(n_states)
    else:
        # Env-default: step-derived count.
        ns = _n_states()

    if ns == 1:
        offsets = (0.0,)
    elif ns == 2:
        offsets = (0.0, 180.0)
    elif ns == 3:
        offsets = (0.0, 120.0, -120.0)
    elif ns == 6:
        # Symmetric pairs around the as-built dihedral; identity first.
        offsets = (0.0, 60.0, -60.0, 120.0, -120.0, 180.0)
    else:
        # Generic: ns evenly spaced offsets starting at 0 (signed).
        step = 360.0 / ns
        offsets = tuple(((i if i <= ns // 2 else i - ns) * step) for i in range(ns))
    n_emitted = 0
    for combo in itertools.product(offsets, repeat=int(n_rotors)):
        if n_emitted >= mc:
            return
        n_emitted += 1
        yield combo


# ------------------------------------------------------------------
# Rotamer application
# ------------------------------------------------------------------


def _atoms_on_side(
    mol,
    src: int,
    dst: int,
    *,
    frozen_idxs: Optional[Sequence[int]] = None,
) -> List[int]:
    """BFS the molecular graph starting at ``dst`` while forbidding the
    ``src``-``dst`` bond.  Returns the list of atom indices reachable from
    ``dst`` (the rotating subtree).  If the bond is in a ring the BFS
    reaches ``src`` too; we exclude ``src`` from the result and the caller
    refuses to rotate ring bonds anyway.

    ``frozen_idxs`` (Bug #2 fix, 2026-06-07): atom indices that act as
    BFS barriers — the metal atom (and any other frozen sites) is added
    to ``seen`` before the BFS starts so the rotating subtree can NEVER
    include the metal.  This catches the metallacycle case: when a
    chelate ring is closed through the metal, RDKit's dative bonds may
    not be recognised as ring-forming, ``find_rotatable_bonds`` returns
    the bond as rotatable, and the BFS from ``dst`` reaches the metal
    through the other donor.  Without the barrier the metal gets
    rotated with the subtree and ends up off-origin in the conformer
    frame.  With the barrier the subtree truncates at the metal, the
    BFS still reaches a valid (but partial) subtree on the rotor side;
    the caller's topology gate then sees an inconsistent geometry and
    rejects the rotamer entirely — both outcomes preserve the metal-at-
    origin invariant.
    """
    n = mol.GetNumAtoms()
    seen = [False] * n
    seen[src] = True
    if frozen_idxs is not None:
        for fi in frozen_idxs:
            try:
                fii = int(fi)
            except (TypeError, ValueError):
                continue
            if 0 <= fii < n:
                seen[fii] = True
    stack = [int(dst)]
    seen[int(dst)] = True
    out: List[int] = []
    while stack:
        cur = stack.pop()
        out.append(cur)
        try:
            for nb in mol.GetAtomWithIdx(int(cur)).GetNeighbors():
                ni = int(nb.GetIdx())
                if not seen[ni]:
                    seen[ni] = True
                    stack.append(ni)
        except Exception:
            continue
    return out


def _rotation_matrix(axis: np.ndarray, theta_rad: float) -> np.ndarray:
    """Rodrigues rotation matrix around ``axis`` by ``theta_rad``."""
    a = np.asarray(axis, dtype=float)
    n = float(np.linalg.norm(a))
    if n < 1e-12:
        return np.eye(3)
    a = a / n
    c, s = float(np.cos(theta_rad)), float(np.sin(theta_rad))
    C = 1.0 - c
    x, y, z = float(a[0]), float(a[1]), float(a[2])
    return np.array([
        [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C   ],
    ], dtype=float)


def _collect_metal_indices(mol) -> List[int]:
    """Return the list of atom indices that are metals (universal —
    delegates to :func:`_is_metal_atom` for the atomic-number check).
    """
    out: List[int] = []
    try:
        for a in mol.GetAtoms():
            if _is_metal_atom(a):
                out.append(int(a.GetIdx()))
    except Exception:
        pass
    return out


def apply_rotamer_config(
    mol,
    coords: np.ndarray,
    rotors: Sequence[Tuple[int, int]],
    offsets_deg: Sequence[float],
) -> np.ndarray:
    """Return a new coordinate array obtained by rotating each subtree by
    ``offsets_deg[k]`` around the axis of ``rotors[k] = (b1, b2)``.

    Subtree = atoms reachable from ``b2`` via BFS that does NOT traverse
    the ``(b1, b2)`` bond AND DOES NOT TRAVERSE THE METAL (Bug #2 fix,
    2026-06-07).  Each rotation is applied IN ORDER over the rotors list;
    the result of one rotation is the input of the next.  Atoms not in
    any subtree are left untouched.

    Metal-as-BFS-barrier:
      RDKit dative bonds (``[Fe-]`` style) sometimes are not perceived as
      ring-forming, so a chelate-ring single bond can pass the
      :func:`find_rotatable_bonds` filter (no `IsInRing()`).  Without a
      metal barrier the BFS from ``b2`` would walk the chelate ring and
      reach the metal, then keep walking through the other donor of the
      chelate back to the rotor side — the metal would be rotated WITH
      the subtree and land off-origin.  Adding the metal indices to the
      BFS frozen set truncates the subtree at the metal and the rotation
      acts only on the genuine rotor side.  The downstream
      ``rotamer_topology_gate`` then rejects any rotation that left the
      structure inconsistent.

    Deterministic: identical inputs -> identical outputs.
    """
    P = np.asarray(coords, dtype=float).copy()
    if not rotors or not list(offsets_deg):
        return P
    metal_idxs = _collect_metal_indices(mol)
    n_rot = min(len(rotors), len(list(offsets_deg)))
    offs = list(offsets_deg)
    for k in range(n_rot):
        offset = float(offs[k])
        if abs(offset) < 1e-9:
            continue
        b1, b2 = int(rotors[k][0]), int(rotors[k][1])
        axis = P[b2] - P[b1]
        if float(np.linalg.norm(axis)) < 1e-6:
            continue
        sub = _atoms_on_side(mol, b1, b2, frozen_idxs=metal_idxs)
        if not sub:
            continue
        # Defensive: even if a metal somehow ended up in the subset
        # (e.g. the rotor itself touches the metal — should not happen
        # because find_rotatable_bonds filters metal-incident bonds, but
        # we guard against future regressions), drop it.
        if metal_idxs:
            metal_set = set(int(mi) for mi in metal_idxs)
            sub = [int(ai) for ai in sub if int(ai) not in metal_set]
            if not sub:
                continue
        R = _rotation_matrix(axis, np.radians(offset))
        pivot = P[b1].copy()
        for ai in sub:
            P[ai] = pivot + (P[ai] - pivot) @ R.T
    return P


# ------------------------------------------------------------------
# Symmetry-priority pre-ranking (Phase 2, 2026-06-03)
# ------------------------------------------------------------------


_ENV_SYMMETRY_PRIORITY = "DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS"


def _env_symmetry_priority_rotamers(default: bool = False) -> bool:
    """Return True when the symmetry-priority rotamer pre-ranking hook is
    on (env ``DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS=1``).  Default OFF
    -> byte-identical to HEAD.
    """
    v = os.environ.get(_ENV_SYMMETRY_PRIORITY, "1" if default else "0")
    return str(v).strip() in ("1", "true", "True", "on", "yes", "YES")


def rerank_rotamers_by_symmetry(
    variants: Sequence[Tuple[np.ndarray, str]],
    atomic_numbers: Sequence,
    tol: float = 0.1,
) -> List[Tuple[np.ndarray, str]]:
    """Sort rotamer variants by detected point-group order DESC.

    The pre-ranking strategy puts the most-symmetric whole-complex
    configurations first in the enumeration, so downstream consumers
    (Mogul polish, RMSD-dedup, conformer selection) preferentially see
    them as the head of the list.  Ties break by the variant's original
    input position (stable sort) -> deterministic.

    Parameters
    ----------
    variants : sequence of (P, label) tuples
        As yielded by :func:`enumerate_single_bond_rotamers`.
    atomic_numbers : sequence of int or element symbols
        Aligned with the variant coordinate arrays.
    tol : float, default 0.1
        Position tolerance for the point-group detector.

    Returns
    -------
    list of (P, label) tuples
        Re-sorted; highest detected point-group order first.

    Universal: depends only on coordinates and atomic numbers; no SMILES,
    no metal-specific heuristics.  Deterministic (stable sort).
    """
    try:
        from delfin.fffree.conformer_symmetry import detect_point_group_order
    except Exception:
        # No symmetry module -> return as-is.
        return list(variants)
    items = list(variants)
    if not items:
        return items
    scored: List[Tuple[int, int, Tuple[np.ndarray, str]]] = []
    for orig_idx, (P, lab) in enumerate(items):
        try:
            pg = detect_point_group_order(P, list(atomic_numbers), tol=tol)
        except Exception:
            pg = 1
        # Negative pg so descending sort puts max first; original index is
        # the secondary key for stable tiebreak.
        scored.append((-int(pg), int(orig_idx), (P, lab)))
    scored.sort(key=lambda t: (t[0], t[1]))
    return [t[2] for t in scored]


# ------------------------------------------------------------------
# Public API: enumerate_single_bond_rotamers
# ------------------------------------------------------------------


def enumerate_single_bond_rotamers(
    mol,
    coords: np.ndarray,
    max_rotors: Optional[int] = None,
    max_configs: Optional[int] = None,
    n_states: Optional[int] = None,
    step_degrees: Optional[float] = None,
    atomic_numbers: Optional[Sequence] = None,
    syms: Optional[Sequence[str]] = None,
) -> Iterator[Tuple[np.ndarray, str]]:
    """Yield ``(coords_variant, label)`` for each enumerated rotamer config.

    Parameters
    ----------
    mol : RDKit Mol
        Source molecule (with explicit H if you want H rotation to follow).
    coords : (N, 3) array
        Coordinates aligned with ``mol`` atom indices.
    max_rotors : int, optional
        Cap on the number of rotatable bonds considered; default from env
        ``DELFIN_FFFREE_ENUMERATE_ROTAMERS_MAX_ROTORS`` (5).
    max_configs : int, optional
        Cap on the number of yielded variants; default from env
        ``DELFIN_FFFREE_ENUMERATE_ROTAMERS_MAX_CONFIGS`` (243).
    n_states : int, optional
        Number of dihedral states per rotor; default 3 (gauche+/anti/gauche-).

    Yields
    ------
    (P_variant, label) : tuple
        ``P_variant`` is a fresh (N, 3) array; the first yielded variant is
        the original ``coords`` (all-zero offset = identity).  Labels are
        of the form ``"rot_+120+0-120"`` encoding the offsets per rotor.

    Note
    ----
    This helper does NOT require ``DELFIN_FFFREE_ENUMERATE_ROTAMERS=1`` to
    run -- the env-flag is honored by the caller (converter_backend).
    Calling this directly bypasses the gate (useful for tests).
    """
    if mol is None:
        return
    if coords is None or len(coords) == 0:
        return
    P0 = np.asarray(coords, dtype=float)
    rotors = find_rotatable_bonds(mol)
    cap = int(max_rotors) if max_rotors is not None else _max_rotors()
    rotors = rotors[: max(0, cap)]
    # Filter rotors for which we can pick reference atoms (defensive; the
    # bond filter already guarantees both sides have heavy neighbours).
    keep: List[Tuple[int, int]] = []
    for b1, b2 in rotors:
        if _dihedral_reference_atoms(mol, b1, b2) is None:
            continue
        keep.append((b1, b2))
    rotors = keep
    if not rotors:
        # No rotors -> yield the input unchanged.
        yield P0.copy(), "rot_identity"
        return

    # Optional topology-gate (universal hard-gate against rotations that
    # break SMILES topology / form spurious bonds / overfill the M-shell).
    # The gate is env-gated (default OFF -> byte-identical with HEAD; auto
    # ON under DELFIN_FFFREE_MOGUL_PRIMARY=1).
    try:
        from .rotamer_topology_gate import (
            _env_on as _topology_gate_on,
            extract_expected_bonds as _topology_expected_bonds,
            expected_metal_cn as _topology_expected_cn,
            rotation_preserves_topology as _topology_preserves,
        )
    except Exception:
        _topology_gate_on = lambda: False  # noqa: E731
        _topology_preserves = None
        _topology_expected_bonds = None
        _topology_expected_cn = None

    _gate_active = bool(_topology_gate_on()) and _topology_preserves is not None
    if _gate_active:
        # Pre-derive once per molecule -- both lookups are O(N).
        try:
            _exp_bonds = (
                _topology_expected_bonds(mol)
                if _topology_expected_bonds is not None else None
            )
        except Exception:
            _exp_bonds = None
        # Element symbols: prefer caller-supplied; else derive from mol.
        if syms is not None:
            _syms_l = [str(s) for s in syms]
        else:
            try:
                _syms_l = [
                    str(mol.GetAtomWithIdx(i).GetSymbol())
                    for i in range(mol.GetNumAtoms())
                ]
            except Exception:
                _syms_l = None
        try:
            _exp_cn = (
                _topology_expected_cn(mol, _syms_l)
                if (_topology_expected_cn is not None and _syms_l is not None)
                else None
            )
        except Exception:
            _exp_cn = None
    else:
        _exp_bonds = None
        _exp_cn = None
        _syms_l = None

    def _stream():
        identity_yielded = False
        for combo in enumerate_rotamer_configs(
            len(rotors),
            n_states=n_states,
            max_configs=max_configs,
            step_degrees=step_degrees,
        ):
            lab_parts = []
            for offset in combo:
                if abs(offset) < 1e-9:
                    lab_parts.append("0")
                elif offset > 0:
                    lab_parts.append(f"+{int(round(offset))}")
                else:
                    lab_parts.append(f"{int(round(offset))}")
            label = "rot_" + "_".join(lab_parts)
            Pv = apply_rotamer_config(mol, P0, rotors, combo)
            # Topology hard-gate: skip rotations that break the SMILES bond
            # graph, form spurious bonds, overfill the M-shell, or create
            # H-H / X-H collisions.  The IDENTITY (all-zero offset) is
            # always emitted -- it is the input geometry, must not be
            # filtered, and downstream code relies on it being first.
            is_identity = all(abs(o) < 1e-9 for o in combo)
            if _gate_active and not is_identity and _syms_l is not None:
                try:
                    keep = _topology_preserves(
                        _syms_l, Pv, mol=mol,
                        expected_bonds=_exp_bonds,
                        expected_cn=_exp_cn,
                    )
                except Exception:
                    keep = True
                if not keep:
                    continue
            if is_identity:
                identity_yielded = True
            yield Pv, label
        # Defensive: if the identity was somehow filtered (max_configs=0
        # path), still emit the input so downstream callers always see at
        # least one variant.
        if not identity_yielded:
            yield P0.copy(), "rot_identity"

    # Phase 2 hook: when the symmetry-priority env flag is on AND we have
    # atomic-number information, materialise the stream and re-rank by
    # detected whole-complex point-group order before yielding.  Default
    # OFF -> byte-identical streaming behaviour preserved.
    if _env_symmetry_priority_rotamers() and atomic_numbers is not None:
        ordered = rerank_rotamers_by_symmetry(list(_stream()), atomic_numbers)
        for item in ordered:
            yield item
    else:
        for item in _stream():
            yield item


# ------------------------------------------------------------------
# Standalone smoke check
# ------------------------------------------------------------------


if __name__ == "__main__":
    from rdkit import Chem
    from rdkit.Chem import AllChem
    for smi, name in [("CCCC", "butane"),
                      ("OCCO", "ethylene glycol"),
                      ("CCCCC", "pentane"),
                      ("C", "methane"),
                      ("c1ccccc1", "benzene (no rotors)")]:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(m, randomSeed=42)
        P = m.GetConformer().GetPositions()
        rotors = find_rotatable_bonds(m)
        n = sum(1 for _ in enumerate_single_bond_rotamers(m, P))
        print(f"{name:<32} rotors={len(rotors):>2}  variants={n:>3}")
