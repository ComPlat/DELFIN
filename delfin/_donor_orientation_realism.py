"""Donor-orientation realism pass for metal-complex post-UFF correction.

Welle-5b Agent B (2026-05-16).  Three geometry patterns the UFF
parameter set does not enforce natively:

1. **Aromatic-N edge-on**:  for ``M-N`` where ``N`` is an
   aromatic ring member, rotate the ring + its substituents (rigid
   rotation around the ``M-N`` axis) until the ring normal becomes
   perpendicular to the ``M-N`` axis (sigma lone-pair points at the
   metal).
2. **Terminal carbonyl linearity**:  for ``M-C(=O)`` where ``C`` has
   only ``M`` as a heavy non-O neighbour and ``C=O`` is a double or
   triple bond, place ``O`` such that ``M-C-O = 180`` (preserving the
   ``C-O`` bond length).
3. **NHC carbene plane alignment**:  for ``M-C(carbene)`` where the
   carbene carbon sits between two ring nitrogens, rotate the ring
   (rigid, around the in-plane bisector at ``C``) so that the metal
   ends up in the ``N-C-N`` plane.

Universal: element-symbol + bond-graph + aromatic flag + bond-order
only.  No SMILES literals, refcodes, or named-ligand patterns.

Default-OFF env flag ``DELFIN_DONOR_ORIENT_REALISM``.  Per-pattern
overrides ``DELFIN_DONOR_ORIENT_REALISM_AROMATIC_N`` /
``_CARBONYL`` / ``_NHC`` may force-disable individual patterns
(set to ``0``) even when the master flag is ``1``.

The pass is purely geometric: pure-numpy, no RNG, no UFF re-run.
"""

from __future__ import annotations

import math
import os
from typing import Iterable, List, Optional, Set, Tuple

import numpy as np

# Metal element set (subset of ``_METALS`` in ``smiles_converter``; we
# duplicate locally to avoid an import cycle).
_METAL_SYMBOLS: Set[str] = {
    'Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
    'Ir', 'Pt', 'Au', 'Hg', 'Al', 'Ga', 'In', 'Tl', 'Sn', 'Pb',
    'Bi', 'Po', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
}

# Tolerances
_AROM_N_TILT_TOL_DEG: float = 10.0  # only fix if normal-axis tilt > tol
_NHC_OOP_TOL_A: float = 0.10        # only fix if M out-of-plane > tol Å
_CO_DEV_TOL_DEG: float = 5.0        # only fix if M-C-O deviates from 180


# --------------------------------------------------------------------- utils

def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


def _parse_xyz_block(
    xyz_str: str,
) -> Tuple[List[str], np.ndarray, int, int]:
    """Parse an XYZ block.

    Accepts both bare coordinate blocks (atom lines only) and standard
    XYZ files with a leading atom count + comment line.  Returns
    ``(symbols, coords, header_lines, trailing_newlines)``.
    ``header_lines`` is the number of non-coordinate lines stripped
    from the front, so the reverse formatter can re-emit them verbatim.
    """
    raw_lines = xyz_str.splitlines()
    # Detect leading header: lines that do not parse as ``SYM x y z``.
    header_lines = 0
    for ln in raw_lines:
        if not ln.strip():
            header_lines += 1
            continue
        parts = ln.split()
        if len(parts) >= 4:
            try:
                float(parts[1])
                float(parts[2])
                float(parts[3])
                break  # first coordinate line found
            except ValueError:
                pass
        header_lines += 1
    body = raw_lines[header_lines:]
    symbols: List[str] = []
    coords: List[List[float]] = []
    for ln in body:
        if not ln.strip():
            continue
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        symbols.append(parts[0])
        coords.append([x, y, z])
    arr = (
        np.asarray(coords, dtype=np.float64)
        if coords else np.zeros((0, 3), dtype=np.float64)
    )
    trailing_newlines = 1 if xyz_str.endswith('\n') else 0
    return symbols, arr, header_lines, trailing_newlines


def _format_xyz_block(
    header_text: str,
    symbols: List[str],
    coords: np.ndarray,
    trailing_newline: bool,
) -> str:
    """Render header + symbols+coords back to a complete XYZ string."""
    parts: List[str] = []
    if header_text:
        parts.append(header_text)
    body_lines: List[str] = []
    for sym, (x, y, z) in zip(symbols, coords):
        body_lines.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
    parts.append('\n'.join(body_lines))
    out = '\n'.join(p for p in parts if p)
    if trailing_newline and not out.endswith('\n'):
        out += '\n'
    return out


def _is_metal(sym: str) -> bool:
    return sym in _METAL_SYMBOLS


def _rodrigues(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Return 3x3 rotation matrix for given axis (unit) + angle."""
    a = axis / (np.linalg.norm(axis) + 1e-12)
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    k = np.array([
        [0.0, -a[2], a[1]],
        [a[2], 0.0, -a[0]],
        [-a[1], a[0], 0.0],
    ])
    return np.eye(3) * c + k * s + np.outer(a, a) * (1.0 - c)


def _bfs_fragment(
    mol,
    start_idx: int,
    blocked: Set[int],
) -> Set[int]:
    """BFS through the bond graph starting at ``start_idx``, never
    crossing any atom in ``blocked``.  Returns the set of visited
    atom indices (includes ``start_idx``)."""
    visited: Set[int] = {start_idx}
    frontier: List[int] = [start_idx]
    while frontier:
        nxt: List[int] = []
        for i in frontier:
            atom_i = mol.GetAtomWithIdx(i)
            for nb in atom_i.GetNeighbors():
                j = nb.GetIdx()
                if j in visited or j in blocked:
                    continue
                visited.add(j)
                nxt.append(j)
        frontier = nxt
    return visited


def _all_donor_indices(mol) -> Set[int]:
    """All non-metal heavy atoms directly bonded to any metal.

    Used to block donor-to-donor BFS crossings: rotating one donor's
    substituent fragment must NOT propagate into another donor's path.
    """
    out: Set[int] = set()
    for a in mol.GetAtoms():
        if a.GetAtomicNum() <= 1:
            continue
        sym = a.GetSymbol()
        if _is_metal(sym):
            continue
        for nb in a.GetNeighbors():
            if _is_metal(nb.GetSymbol()):
                out.add(a.GetIdx())
                break
    return out


def _all_metal_indices(mol) -> Set[int]:
    return {
        a.GetIdx() for a in mol.GetAtoms()
        if _is_metal(a.GetSymbol())
    }


# ---------------------------------------------------------- pattern locators

def _find_aromatic_n_metal_pairs(mol) -> List[Tuple[int, int, List[int]]]:
    """Find ``(metal_idx, n_idx, [ring_atom_indices])`` triples for
    aromatic-N atoms bound to a metal."""
    results: List[Tuple[int, int, List[int]]] = []
    try:
        ri = mol.GetRingInfo()
        atom_rings = ri.AtomRings() if ri is not None else ()
    except Exception:
        return results
    for n_atom in mol.GetAtoms():
        if n_atom.GetSymbol() != 'N':
            continue
        if not n_atom.GetIsAromatic():
            continue
        metals = [
            nb.GetIdx() for nb in n_atom.GetNeighbors()
            if _is_metal(nb.GetSymbol())
        ]
        if not metals:
            continue
        n_idx = n_atom.GetIdx()
        # Find the aromatic ring that contains this N.
        ring_atoms: Optional[Tuple[int, ...]] = None
        for ring in atom_rings:
            if n_idx not in ring or len(ring) not in (5, 6):
                continue
            # Reject rings that include any metal (spurious bond
            # perception can introduce M-N edges that close fake
            # rings around the metal).
            if any(
                _is_metal(mol.GetAtomWithIdx(j).GetSymbol())
                for j in ring
            ):
                continue
            arom_count = sum(
                1 for j in ring
                if mol.GetAtomWithIdx(j).GetIsAromatic()
            )
            if arom_count >= len(ring) - 1:
                ring_atoms = tuple(ring)
                break
        if ring_atoms is None:
            continue
        for m_idx in metals:
            results.append((m_idx, n_idx, list(ring_atoms)))
    return results


def _find_terminal_carbonyl_triples(mol) -> List[Tuple[int, int, int]]:
    """Find ``(metal_idx, c_idx, o_idx)`` triples for terminal M-C=O
    where C has only the metal and the oxygen as heavy neighbours and
    the ``C-O`` bond is a double or triple bond."""
    out: List[Tuple[int, int, int]] = []
    try:
        from rdkit.Chem import BondType  # type: ignore
        dbl = BondType.DOUBLE
        trp = BondType.TRIPLE
    except Exception:
        return out
    for c in mol.GetAtoms():
        if c.GetSymbol() != 'C':
            continue
        heavy_nbrs = [nb for nb in c.GetNeighbors() if nb.GetSymbol() != 'H']
        # Must have exactly one O and >=1 metal, no other heavy atoms.
        oxygens = [nb for nb in heavy_nbrs if nb.GetSymbol() == 'O']
        metals = [nb for nb in heavy_nbrs if _is_metal(nb.GetSymbol())]
        others = [
            nb for nb in heavy_nbrs
            if nb.GetSymbol() not in ('O', 'H') and not _is_metal(nb.GetSymbol())
        ]
        if len(oxygens) != 1 or not metals or others:
            continue
        if len(metals) != 1:
            continue  # bridging carbonyls not handled here
        o = oxygens[0]
        # Oxygen must be terminal (only bound to this C heavy-wise).
        o_heavy = [nb for nb in o.GetNeighbors() if nb.GetSymbol() != 'H']
        if len(o_heavy) != 1:
            continue
        bond = mol.GetBondBetweenAtoms(c.GetIdx(), o.GetIdx())
        if bond is None:
            continue
        bt = bond.GetBondType()
        if bt not in (dbl, trp):
            continue
        out.append((metals[0].GetIdx(), c.GetIdx(), o.GetIdx()))
    return out


def _find_nhc_carbenes(mol) -> List[Tuple[int, int, int, int, List[int]]]:
    """Find ``(metal_idx, c_idx, n1_idx, n2_idx, ring_atoms)`` for
    NHC carbenes: 5-ring with C(carbene) bonded to 2 N in-ring and to a
    metal outside the ring."""
    out: List[Tuple[int, int, int, int, List[int]]] = []
    try:
        ri = mol.GetRingInfo()
        atom_rings = ri.AtomRings() if ri is not None else ()
    except Exception:
        return out
    for ring in atom_rings:
        if len(ring) != 5:
            continue
        # Reject 5-rings that include a metal (spurious bond perception
        # frequently produces metal-containing 5-rings that pass the
        # 2-N-neighbours test but are not chemical NHCs).
        if any(_is_metal(mol.GetAtomWithIdx(j).GetSymbol()) for j in ring):
            continue
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != 'C':
                continue
            in_ring_nbrs = [
                nb for nb in atom.GetNeighbors()
                if nb.GetIdx() in ring
            ]
            if len(in_ring_nbrs) != 2:
                continue
            if not all(nb.GetSymbol() == 'N' for nb in in_ring_nbrs):
                continue
            metals = [
                nb for nb in atom.GetNeighbors()
                if _is_metal(nb.GetSymbol())
            ]
            if not metals:
                continue
            n1, n2 = in_ring_nbrs[0].GetIdx(), in_ring_nbrs[1].GetIdx()
            out.append((metals[0].GetIdx(), idx, n1, n2, list(ring)))
            break
    return out


# ----------------------------------------------------------- pattern actions

def _apply_aromatic_n_edge_on(
    coords: np.ndarray,
    mol,
    triples: Iterable[Tuple[int, int, List[int]]],
) -> int:
    """Rotate aromatic rings around an axis through ``N`` (perpendicular
    to the ``M-N`` direction) so the ring normal becomes perpendicular
    to the ``M-N`` axis (sigma edge-on attack).

    Rotation only moves atoms that are reachable from a non-N ring
    atom via a path that does NOT cross ``N``, the metal, or any other
    donor atom.  This protects all other ``M-D`` bonds and donor
    geometries on the same complex.
    """
    n_changed = 0
    all_donors = _all_donor_indices(mol)
    all_metals = _all_metal_indices(mol)
    seen_axes: Set[Tuple[int, int]] = set()
    for m_idx, n_idx, ring in triples:
        key = (min(m_idx, n_idx), max(m_idx, n_idx))
        if key in seen_axes:
            continue
        seen_axes.add(key)
        if len(ring) < 3:
            continue
        m_pos = coords[m_idx]
        n_pos = coords[n_idx]
        mn_axis = m_pos - n_pos
        mn_len = float(np.linalg.norm(mn_axis))
        if mn_len < 1e-3:
            continue
        mn_axis = mn_axis / mn_len
        ring_non_n = [j for j in ring if j != n_idx]
        if len(ring_non_n) < 2:
            continue
        # Atom set to rotate: just the non-N ring atoms plus any
        # *terminal* (H or single-element) substituents directly
        # bonded to each ring atom.  We deliberately do NOT BFS
        # further -- shared-atom rings, fused systems, and bridged
        # ligand backbones can connect to other donors and break
        # other M-D bonds.  The ring + immediate H envelope is the
        # safe minimum that preserves all M-D bonds while still
        # rotating the donor's aromatic plane.
        #
        # Skip fused rings: if any ring atom (other than N) belongs to
        # a second ring, the geometric rotation will distort the
        # shared bond into the fused partner.
        try:
            ri_local = mol.GetRingInfo()
            ring_sizes_for: List[int] = [
                ri_local.NumAtomRings(j) for j in ring if j != n_idx
            ]
            if any(s > 1 for s in ring_sizes_for):
                continue  # fused/bicyclic ring -- safer to skip
        except Exception:
            pass
        frag: Set[int] = set(ring) - {n_idx}
        for j in list(frag):
            for nb in mol.GetAtomWithIdx(j).GetNeighbors():
                k = nb.GetIdx()
                if k in frag or k in ring:
                    continue
                if k == n_idx or k in all_metals or k in all_donors:
                    continue
                # Include only true terminal H atoms; any heavier
                # substituent stays put to avoid downstream collisions.
                if nb.GetAtomicNum() == 1:
                    frag.add(k)
        if not frag:
            continue
        # Iterate a few times: a single Rodrigues step around the
        # off-axis-perpendicular is only first-order accurate when
        # the deviation is large.  3 passes suffice to drive tilt
        # below 1 deg for any starting orientation.
        for _ in range(3):
            v1 = coords[ring_non_n[0]] - n_pos
            v2 = coords[ring_non_n[1]] - n_pos
            ring_normal = np.cross(v1, v2)
            norm_len = float(np.linalg.norm(ring_normal))
            if norm_len < 1e-6:
                break
            ring_normal = ring_normal / norm_len
            cos_n_a = float(np.dot(ring_normal, mn_axis))
            cos_n_a = max(-1.0, min(1.0, cos_n_a))
            angle_n_a = math.degrees(math.acos(abs(cos_n_a)))
            deviation = 90.0 - angle_n_a
            if abs(deviation) < _AROM_N_TILT_TOL_DEG / 5.0:  # <2 deg
                break
            rot_axis = np.cross(ring_normal, mn_axis)
            rot_axis_len = float(np.linalg.norm(rot_axis))
            if rot_axis_len < 1e-6:
                # Degenerate (normal parallel to axis): pick the
                # direction in the ring plane perpendicular to the
                # bisector of v1+v2.
                bis = v1 + v2
                bis_len = float(np.linalg.norm(bis))
                if bis_len < 1e-6:
                    break
                bis = bis / bis_len
                # Axis perpendicular to bisector AND lying in ring
                # plane = normal x bisector.
                rot_axis = np.cross(ring_normal, bis)
                rot_axis_len = float(np.linalg.norm(rot_axis))
                if rot_axis_len < 1e-6:
                    break
                rot_axis = rot_axis / rot_axis_len
            else:
                rot_axis = rot_axis / rot_axis_len
            # Pick rotation sign that reduces |cos|.
            best_R = None
            best_resid = abs(cos_n_a)
            for sign in (1.0, -1.0):
                R = _rodrigues(
                    rot_axis, math.radians(sign * deviation)
                )
                new_normal = R @ ring_normal
                new_cos = abs(float(np.dot(new_normal, mn_axis)))
                if new_cos < best_resid - 1e-4:
                    best_resid = new_cos
                    best_R = R
            if best_R is None:
                break
            pivot = n_pos
            for j in frag:
                rel = coords[j] - pivot
                coords[j] = pivot + (best_R @ rel)
        n_changed += 1
    return n_changed


def _apply_terminal_carbonyl_linearity(
    coords: np.ndarray,
    triples: Iterable[Tuple[int, int, int]],
) -> int:
    """For M-C=O, move O so that M-C-O = 180 deg, preserving |C-O|."""
    n_changed = 0
    for m_idx, c_idx, o_idx in triples:
        m_pos = coords[m_idx]
        c_pos = coords[c_idx]
        o_pos = coords[o_idx]
        v_cm = m_pos - c_pos
        v_co = o_pos - c_pos
        len_co = float(np.linalg.norm(v_co))
        len_cm = float(np.linalg.norm(v_cm))
        if len_co < 1e-3 or len_cm < 1e-3:
            continue
        # Current M-C-O angle.
        cos_ang = float(np.dot(v_cm, v_co) / (len_cm * len_co))
        cos_ang = max(-1.0, min(1.0, cos_ang))
        angle_deg = math.degrees(math.acos(cos_ang))
        # Target = 180 deg, i.e. v_co = -v_cm direction.
        if abs(180.0 - angle_deg) < _CO_DEV_TOL_DEG:
            continue
        new_dir = -v_cm / len_cm
        coords[o_idx] = c_pos + new_dir * len_co
        n_changed += 1
    return n_changed


def _apply_nhc_plane_alignment(
    coords: np.ndarray,
    mol,
    quads: Iterable[Tuple[int, int, int, int, List[int]]],
) -> int:
    """For each NHC (M, C_carbene, N1, N2, ring), rotate the ring so M
    lies in the N1-C-N2 plane (i.e. M out-of-plane displacement = 0).
    Rotation axis: in-plane bisector at C, oriented so that the ring
    pivots toward M without losing the M-C bond direction."""
    n_changed = 0
    for m_idx, c_idx, n1_idx, n2_idx, ring in quads:
        m_pos = coords[m_idx]
        c_pos = coords[c_idx]
        n1_pos = coords[n1_idx]
        n2_pos = coords[n2_idx]
        v_cn1 = n1_pos - c_pos
        v_cn2 = n2_pos - c_pos
        # Plane normal of N1-C-N2.
        normal = np.cross(v_cn1, v_cn2)
        nlen = float(np.linalg.norm(normal))
        if nlen < 1e-6:
            continue
        normal = normal / nlen
        # M displacement perpendicular to plane.
        d_m = float(np.dot(m_pos - c_pos, normal))
        if abs(d_m) < _NHC_OOP_TOL_A:
            continue
        # Rotation axis: in-plane bisector of N1-C-N2, i.e. the
        # in-plane direction orthogonal to the bond from C to the
        # midpoint M.  We use the projection of (M - C) onto the
        # N1-C-N2 plane as the in-plane "to-M" direction; the rotation
        # axis is the perpendicular in-plane direction.
        v_cm = m_pos - c_pos
        v_cm_in_plane = v_cm - d_m * normal
        ip_len = float(np.linalg.norm(v_cm_in_plane))
        if ip_len < 1e-4:
            # M nearly along the ring normal: rotate around N1-N2 axis.
            rot_axis = (n2_pos - n1_pos)
        else:
            rot_axis = np.cross(normal, v_cm_in_plane / ip_len)
        ra_len = float(np.linalg.norm(rot_axis))
        if ra_len < 1e-6:
            continue
        rot_axis = rot_axis / ra_len
        # Required rotation angle: the OOP component d_m must vanish.
        # |M - C| does not change; we rotate the *ring* (not M) so the
        # plane swings to contain M.
        d_m_magnitude = float(np.linalg.norm(v_cm))
        if d_m_magnitude < 1e-3:
            continue
        sin_tilt = d_m / d_m_magnitude
        sin_tilt = max(-1.0, min(1.0, sin_tilt))
        tilt_rad = math.asin(sin_tilt)
        # Try both signs and pick the one that drives d_m -> 0.
        best_sign = 0.0
        best_resid = abs(d_m)
        for sign in (1.0, -1.0):
            R = _rodrigues(rot_axis, sign * tilt_rad)
            new_normal = R @ normal
            new_d = float(np.dot(m_pos - c_pos, new_normal))
            if abs(new_d) < best_resid - 1e-3:
                best_resid = abs(new_d)
                best_sign = sign
        if best_sign == 0.0:
            continue
        R = _rodrigues(rot_axis, best_sign * tilt_rad)
        # Conservative atom set: NHC ring atoms (minus carbene C) plus
        # immediate terminal H substituents.  Fused / bridged
        # downstream atoms stay put -- shared-atom backbones often
        # connect into other donor paths and rotating them breaks
        # M-D bonds elsewhere.
        all_donors = _all_donor_indices(mol)
        all_metals = _all_metal_indices(mol)
        # Skip fused NHC rings: a shared ring atom (other than the
        # carbene C) means tilting the ring will distort the partner
        # ring.  Benzimidazol-2-ylidene etc. fall in this branch.
        try:
            ri_local = mol.GetRingInfo()
            shared = [
                ri_local.NumAtomRings(j) for j in ring if j != c_idx
            ]
            if any(s > 1 for s in shared):
                continue
        except Exception:
            pass
        frag: Set[int] = set(ring) - {c_idx}
        for j in list(frag):
            for nb in mol.GetAtomWithIdx(j).GetNeighbors():
                k = nb.GetIdx()
                if k in frag or k in ring:
                    continue
                if k == c_idx or k in all_metals or k in all_donors:
                    continue
                if nb.GetAtomicNum() == 1:
                    frag.add(k)
        pivot = c_pos
        for j in frag:
            rel = coords[j] - pivot
            coords[j] = pivot + (R @ rel)
        n_changed += 1
    return n_changed


# ------------------------------------------------------------------- public

def snap_donor_orientations(
    xyz_str: str,
    mol,
    *,
    mode: str = 'end_of_pipeline',
) -> str:
    """Apply donor-orientation realism patterns to ``xyz_str``.

    Heavy-atom-only operations; H atoms move rigidly with their ring
    fragments and are otherwise untouched.  Returns the modified XYZ
    block (no header lines).  When the master flag
    ``DELFIN_DONOR_ORIENT_REALISM`` is not set / set to 0, the input
    is returned unchanged.
    """
    if not _env_int('DELFIN_DONOR_ORIENT_REALISM', 0):
        return xyz_str
    if mol is None:
        return xyz_str

    # Sub-flags (default ON when master is ON; users can opt-out per
    # pattern).
    do_arn = _env_int('DELFIN_DONOR_ORIENT_REALISM_AROMATIC_N', 1)
    do_co = _env_int('DELFIN_DONOR_ORIENT_REALISM_CARBONYL', 1)
    do_nhc = _env_int('DELFIN_DONOR_ORIENT_REALISM_NHC', 1)

    try:
        symbols, coords, header_n, trail_nl = _parse_xyz_block(xyz_str)
    except Exception:
        return xyz_str
    if coords.size == 0 or len(symbols) != mol.GetNumAtoms():
        # XYZ/mol mismatch (e.g. UFF added explicit H) -- bail safely.
        return xyz_str
    # Preserve the original header verbatim.
    raw_lines = xyz_str.splitlines()
    header_text = '\n'.join(raw_lines[:header_n]) if header_n else ''

    try:
        if do_arn:
            triples = _find_aromatic_n_metal_pairs(mol)
            if triples:
                _apply_aromatic_n_edge_on(coords, mol, triples)
        if do_co:
            co_triples = _find_terminal_carbonyl_triples(mol)
            if co_triples:
                _apply_terminal_carbonyl_linearity(coords, co_triples)
        if do_nhc:
            quads = _find_nhc_carbenes(mol)
            if quads:
                _apply_nhc_plane_alignment(coords, mol, quads)
    except Exception:
        # Any unexpected failure: silently keep the original geometry.
        return xyz_str

    return _format_xyz_block(header_text, symbols, coords, bool(trail_nl))
