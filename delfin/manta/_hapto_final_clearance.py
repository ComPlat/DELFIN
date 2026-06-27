"""Welle-5f-F: end-of-pipeline 81f8a1f-style M-X clearance post-pass.

The 81f8a1f hapto port (commit 81f8a1f9, Apr-14) introduced
``_enforce_metal_topology`` — a radial push pass that, for every
metal, displaces non-bonded heavy atoms (or whole hapto-groups, as
rigid bodies) so they sit at least ``min_nonbonded`` Å from the
metal.  In HEAD the helper exists (``delfin/smiles_converter.py``
≈ line 18835) and is invoked exactly once, inside
``_select_best_hapto_candidate``.

Per Welle-5c-CV (2026-05-16, only real per-CV champion edge:
**+1.34 pp per-bond, +4.33 pp per-file on the 901-file hapto
intersection**) the 81f8a1f geometry edge is genuinely present
in the post-emit pipeline.  However HEAD adds several emission
paths (HD-TA σ-rotations Iter-3, Baustein-3 angle-corrector
Iter-13, Baustein-4 π-H projection Iter-14, H-VSEPR realism
Welle-5b A, optional Baustein-5 PBD optimizer) that all run
**after** ``_enforce_metal_topology`` and can re-introduce
M-X clashes (rotations move heavy atoms back inside the metal
coordination sphere; rigid-π projections snap H atoms onto a
plane that intersects the M).

This module re-applies the 81f8a1f M-X clearance push as a
**final** post-pass on every emitted ``(xyz, label)`` for hapto-
class systems, bridging the gap between the candidate-select-time
gate and the actual returned XYZ.

Universal: element-symbols + bond-graph + hapto-group detection
only.  No SMILES literals, refcodes, named-ligand patterns.
Default OFF — bit-exact when the env-flag is 0.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Mirror of ``_METAL_SET`` in smiles_converter.py.  Duplicated locally
# to avoid an import cycle (this module is loaded from the dispatch
# helper at import time inside the same module).
_METAL_SET = {
    'Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
    'Ir', 'Pt', 'Au', 'Hg', 'Al', 'Ga', 'In', 'Tl', 'Sn', 'Pb',
    'Bi', 'Po', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
}


def _parse_xyz(xyz_str: str) -> Tuple[List[str], np.ndarray, str, str]:
    """Return (symbols, Nx3 coords, header_line, blank_tail)."""
    lines = xyz_str.splitlines()
    if len(lines) < 2:
        raise ValueError("xyz too short")
    n = int(lines[0].strip())
    header = lines[1] if len(lines) > 1 else ""
    syms: List[str] = []
    rows: List[Tuple[float, float, float]] = []
    for k in range(2, 2 + n):
        parts = lines[k].split()
        if len(parts) < 4:
            raise ValueError(f"bad atom line: {lines[k]!r}")
        syms.append(parts[0])
        rows.append((float(parts[1]), float(parts[2]), float(parts[3])))
    tail = "\n".join(lines[2 + n:])
    return syms, np.asarray(rows, dtype=float), header, tail


def _format_xyz(symbols: Sequence[str], coords: np.ndarray,
                header: str, tail: str) -> str:
    """Rebuild XYZ string preserving header and any trailing data."""
    n = len(symbols)
    out: List[str] = [str(n), header]
    for sym, (x, y, z) in zip(symbols, coords):
        out.append(f"{sym:<3s} {x: .6f} {y: .6f} {z: .6f}")
    body = "\n".join(out)
    if tail:
        body += "\n" + tail
    if not body.endswith("\n"):
        body += "\n"
    return body


def _metal_indices(symbols: Sequence[str]) -> List[int]:
    return [i for i, s in enumerate(symbols) if s in _METAL_SET]


def _bonded_indices_from_mol(mol) -> Dict[int, set]:
    """Map atom index -> set of bonded neighbour indices via RDKit graph."""
    out: Dict[int, set] = {}
    if mol is None:
        return out
    for atom in mol.GetAtoms():
        out[atom.GetIdx()] = {nb.GetIdx() for nb in atom.GetNeighbors()}
    return out


def _hapto_group_map(mol, find_hapto_groups) -> Tuple[Dict[int, int],
                                                       Dict[int, List[int]]]:
    """Return (atom_idx -> group_id, group_id -> [atom indices])."""
    hapto_of: Dict[int, int] = {}
    group_atoms: Dict[int, List[int]] = {}
    if mol is None or find_hapto_groups is None:
        return hapto_of, group_atoms
    try:
        groups = find_hapto_groups(mol)
    except Exception:
        return hapto_of, group_atoms
    for gi, (_metal_idx, grp) in enumerate(groups):
        group_atoms[gi] = list(grp)
        for a in grp:
            hapto_of[a] = gi
    return hapto_of, group_atoms


def enforce_m_x_clearance_xyz(
    xyz_str: str,
    mol,
    find_hapto_groups,
    *,
    min_nonbonded: float = 2.5,
    max_passes: int = 10,
    push_gain: float = 1.1,
) -> str:
    """Apply 81f8a1f-style radial M-X clearance push to an XYZ string.

    For every metal in ``mol``, any heavy non-H, non-metal atom that
    is **not** bonded to the metal in the molecular graph but sits
    closer than ``min_nonbonded`` Å in the coordinates gets pushed
    radially outward by ``push_gain`` × intrusion.  Atoms belonging
    to a hapto group are translated rigidly with the rest of their
    group so ring geometry is preserved.

    Bit-exact passthrough if the XYZ parse fails, if there are no
    metals, or if no clash needs fixing.  Hydrogens never move
    (matches 81f8a1f's ``_enforce_metal_topology`` policy).
    """
    if mol is None:
        return xyz_str

    try:
        symbols, coords, header, tail = _parse_xyz(xyz_str)
    except Exception:
        return xyz_str

    n_mol = mol.GetNumAtoms()
    if n_mol != len(symbols):
        # Heavy/H asymmetry — XYZ may include H not in graph; refuse.
        return xyz_str

    metal_idx = _metal_indices(symbols)
    if not metal_idx:
        return xyz_str

    bonded = _bonded_indices_from_mol(mol)
    hapto_of, group_atoms = _hapto_group_map(mol, find_hapto_groups)

    coords = coords.copy()
    any_moved_at_all = False
    for _pass in range(max_passes):
        any_moved = False
        for mi in metal_idx:
            mpos = coords[mi]
            for ai in range(len(symbols)):
                if ai == mi or ai in bonded.get(mi, set()):
                    continue
                sym_ai = symbols[ai]
                if sym_ai == 'H' or sym_ai in _METAL_SET:
                    continue
                vec = coords[ai] - mpos
                d = float(np.linalg.norm(vec))
                if d >= min_nonbonded or d < 1e-8:
                    continue
                push_dir = vec / d
                push_dist = (min_nonbonded - d) * push_gain
                gi = hapto_of.get(ai)
                if gi is not None:
                    for ha in group_atoms[gi]:
                        coords[ha] = coords[ha] + push_dir * push_dist
                else:
                    coords[ai] = coords[ai] + push_dir * push_dist
                any_moved = True
                any_moved_at_all = True
        if not any_moved:
            break

    if not any_moved_at_all:
        return xyz_str

    return _format_xyz(symbols, coords, header, tail)
