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


def _atoms_on_side(mol, src: int, dst: int) -> List[int]:
    """BFS the molecular graph starting at ``dst`` while forbidding the
    ``src``-``dst`` bond.  Returns the list of atom indices reachable from
    ``dst`` (the rotating subtree).  If the bond is in a ring the BFS
    reaches ``src`` too; we exclude ``src`` from the result and the caller
    refuses to rotate ring bonds anyway.
    """
    n = mol.GetNumAtoms()
    seen = [False] * n
    seen[src] = True
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


def apply_rotamer_config(
    mol,
    coords: np.ndarray,
    rotors: Sequence[Tuple[int, int]],
    offsets_deg: Sequence[float],
) -> np.ndarray:
    """Return a new coordinate array obtained by rotating each subtree by
    ``offsets_deg[k]`` around the axis of ``rotors[k] = (b1, b2)``.

    Subtree = atoms reachable from ``b2`` via BFS that does NOT traverse
    the ``(b1, b2)`` bond.  Each rotation is applied IN ORDER over the
    rotors list; the result of one rotation is the input of the next.
    Atoms not in any subtree are left untouched.

    Robust to ring bonds: if the rotor accidentally points at a ring bond
    (caller should have filtered) the BFS returns the whole ring and the
    rotation distorts the ring -- the caller is expected to honor
    :func:`find_rotatable_bonds`.

    Deterministic: identical inputs -> identical outputs.
    """
    P = np.asarray(coords, dtype=float).copy()
    if not rotors or not list(offsets_deg):
        return P
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
        sub = _atoms_on_side(mol, b1, b2)
        if not sub:
            continue
        R = _rotation_matrix(axis, np.radians(offset))
        pivot = P[b1].copy()
        for ai in sub:
            P[ai] = pivot + (P[ai] - pivot) @ R.T
    return P


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
    for combo in enumerate_rotamer_configs(
        len(rotors),
        n_states=n_states,
        max_configs=max_configs,
        step_degrees=step_degrees,
    ):
        # Build the label deterministically: e.g. rot_+120_0_-120
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
        yield Pv, label


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
