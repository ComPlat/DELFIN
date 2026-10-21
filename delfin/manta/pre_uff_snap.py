"""Pre-UFF metal-donor snap, topology gate switches, metalloid clamps and d8 square-planar flattening of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _delfin_env_int,
)
from delfin.manta.ml_tables import (
    Point3D,
    RDKIT_AVAILABLE,
    _METALLOID_MD_DONORS,
    _METAL_SET,
    _get_ml_bond_length,
    _ml_bond_kind,
)

logger = get_logger("delfin.smiles_converter")


# ---------------------------------------------------------------------------
# Pre-UFF M-D snap and topology gate (universal, env-flag gated, default OFF).
#
# Failure pattern these helpers address (universal, no SMILES-specific paths):
#
#   When the metal is an unparametrized transition metal in Open Babel UFF
#   (Ni(II), Pd(II), Pt(II), Cu(II), Fe(II/III), Co(II/III), Ru(II/III),
#    Rh(III), Ir(III), Mn(II/III), etc. with explicit formal charge), OB UFF
#   triggers a HARD-fallback that freezes atom positions instead of relaxing
#   them. Whatever metal-donor distances the pre-UFF geometry already has are
#   the distances the final XYZ will have. ETKDG by itself produces M-D
#   distances in the 1.4-1.6 A range (organic-bond-like) or sometimes in the
#   2.8-3.5 A range (anti-attractive). Both extremes are catastrophic.
#
# The fix is purely geometric: snap each M-D pair to the element-pair ideal
# distance _get_ml_bond_length(M, D) before UFF runs, by translating the
# bonded-fragment attached to D along the (D - M) direction. Universal
# because it only consults element symbols and the bond graph.
# ---------------------------------------------------------------------------

def _pre_uff_md_snap_enabled() -> bool:
    """Return True iff DELFIN_PRE_UFF_MD_SNAP is set to a truthy value."""
    raw = os.environ.get("DELFIN_PRE_UFF_MD_SNAP", "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


def _pre_uff_topology_gate_enabled() -> bool:
    """Return True iff DELFIN_PRE_UFF_TOPOLOGY_GATE is set to a truthy value."""
    raw = os.environ.get("DELFIN_PRE_UFF_TOPOLOGY_GATE", "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


def _pre_uff_topology_gate_v2_enabled() -> bool:
    """Return True iff DELFIN_PRE_UFF_TOPOLOGY_GATE_V2 is set to a truthy value.

    V2 = empirically calibrated per-(M,D)-pair tolerance band. When the
    base gate is also enabled, V2 takes precedence. Default OFF.
    """
    raw = os.environ.get("DELFIN_PRE_UFF_TOPOLOGY_GATE_V2", "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


# Per-donor-element relaxed-upper-bound table (empirical, Welle-3 T5.2).
# Heavier donors with systematic ML-table under-estimate need a wider upper
# limit so legitimate CSD-realistic distances are not falsely rejected.
# Numbers from per-(M,D) p95 on smoke_500 master data + 1.5-sigma safety.
_GATE_V2_DONOR_HIGH: Dict[str, float] = {
    "As": 1.30, "Sb": 1.30, "Bi": 1.30, "Sn": 1.40,
    "Te": 1.25, "Pb": 1.30, "Si": 1.20, "Ge": 1.25,
}


# Per-donor-element relaxed-lower-bound table. Multibond M=O/M=N/M=S
# (carbene, oxo, nitrido) routinely sit at 0.7-0.8 of single-bond ideal.
_GATE_V2_DONOR_LOW_MULTIBOND: Dict[str, float] = {
    "O": 0.70, "N": 0.72, "S": 0.72, "C": 0.70, "Se": 0.72,
}


def _cn5_enum_complete_enabled() -> bool:
    """Return True iff DELFIN_CN5_ENUM_COMPLETE is set to a truthy value."""
    raw = os.environ.get("DELFIN_CN5_ENUM_COMPLETE", "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


def _bfs_ligand_fragment(
    mol,
    donor_idx: int,
    metal_indices: set,
) -> set:
    """Return the set of atom indices reachable from *donor_idx* without
    crossing any metal in *metal_indices*. Pure graph traversal.

    Used by _snap_md_distances_to_ideal to translate the donor's ligand
    fragment as a rigid body when correcting the M-D distance.
    """
    visited = {donor_idx}
    stack = [donor_idx]
    while stack:
        cur = stack.pop()
        for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in visited or ni in metal_indices:
                continue
            visited.add(ni)
            stack.append(ni)
    return visited


def _snap_md_distances_to_ideal(
    mol,
    conf_id: int,
    snap_tolerance: float = 0.10,
) -> int:
    """Snap every M-D bonded distance to the ideal length from
    ``_get_ml_bond_length(M_sym, D_sym)`` by rigidly translating the donor's
    BFS fragment along the (D - M) direction.

    Only operates on bonded M-D pairs (from the graph). Skips H donors.
    Skips pairs already within ``snap_tolerance * ideal`` of the target.
    Skips chelate donors processed once per chelate: each donor is moved
    independently, which may distort the bite, but downstream UFF/cleanup
    fixes that — the goal here is solely to escape the unparam-TM
    HARD-fallback freeze trap.

    Args:
        mol: RDKit Mol with conformer ``conf_id``.
        conf_id: Conformer index to modify in place.
        snap_tolerance: Skip pairs whose relative deviation is below this
            fraction (default 10%).

    Returns the number of M-D pairs snapped.
    """
    if not RDKIT_AVAILABLE:
        return 0
    try:
        import numpy as np
    except Exception:
        return 0
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return 0

    metal_indices = {
        a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET
    }
    if not metal_indices:
        return 0

    n_snapped = 0
    # Iterate over a sorted list so behavior is deterministic.
    for m_idx in sorted(metal_indices):
        m_atom = mol.GetAtomWithIdx(m_idx)
        m_sym = m_atom.GetSymbol()
        mp = conf.GetAtomPosition(m_idx)
        m_pos = np.array([mp.x, mp.y, mp.z], dtype=float)
        for nbr in m_atom.GetNeighbors():
            d_idx = nbr.GetIdx()
            if d_idx in metal_indices:
                continue  # M-M bonds
            if nbr.GetAtomicNum() <= 1:
                continue  # H donors not snapped
            # Welle-3 T3.3: skip bridging donors (donor bonded to >=2 metals).
            # Single-translation cannot satisfy two M-D ideals simultaneously;
            # without a skip, the second metal's pass undoes the first's
            # snap, leaving Cu-O 1.03 A in a Cu-O-Cu test (see report).
            # Leave bridging M-D for UFF distance pins to balance.
            n_metal_nbrs = sum(
                1 for nb in nbr.GetNeighbors() if nb.GetIdx() in metal_indices
            )
            if n_metal_nbrs >= 2:
                continue
            d_sym = nbr.GetSymbol()
            dp = conf.GetAtomPosition(d_idx)
            d_pos = np.array([dp.x, dp.y, dp.z], dtype=float)
            cur_d = float(np.linalg.norm(d_pos - m_pos))
            if cur_d < 1e-8:
                continue
            try:
                target_d = float(_get_ml_bond_length(m_sym, d_sym, _ml_bond_kind(mol, m_idx, d_idx) if _delfin_env_int("DELFIN_FFFREE_ME_BOND_LEN", 0) else "sigma"))
            except Exception:
                continue
            if target_d <= 0:
                continue
            rel = abs(cur_d - target_d) / target_d
            if rel < snap_tolerance:
                continue

            # BFS fragment downstream of this donor, never crossing metals.
            frag = _bfs_ligand_fragment(mol, d_idx, metal_indices)

            # Translation vector: shift fragment along (D - M) so the new
            # distance is target_d.
            unit = (d_pos - m_pos) / cur_d
            new_d_pos = m_pos + unit * target_d
            delta = new_d_pos - d_pos

            for fi in frag:
                fp = conf.GetAtomPosition(fi)
                fv = np.array([fp.x, fp.y, fp.z], dtype=float)
                new_fv = fv + delta
                conf.SetAtomPosition(
                    fi,
                    Point3D(float(new_fv[0]), float(new_fv[1]), float(new_fv[2])),
                )
            n_snapped += 1

    return n_snapped


def _clamp_metalloid_md_xyz(xyz_delfin: str, mol_template) -> str:
    """POST-UFF metalloid M-D clamp (root fix, env-gated DELFIN_FFFREE_METALLOID_MD_CLAMP).

    OB-UFF has no force-field parameters for the soft, large heavy-metalloid sigma-donors
    (Sb/As/Bi/Te/Se/Ge/Sn/Pb) and, when the seat target is short, COLLAPSES the M-metalloid bond
    below its ideal length -- QIGFOF Ag-Sb lands at 2.01 A against an ideal 2.65, wrecking realism
    even though the seat target (r_Ag+r_Sb = 2.84) was correct.  This runs AFTER UFF and rigidly
    re-snaps ONLY metalloid M-D bonds to _get_ml_bond_length by translating the donor's BFS fragment
    along the M-D axis -- the same proven mechanism as _snap_md_distances_to_ideal, but at the XYZ
    level and metalloid-only, so the N/O/P/S bonds UFF handles well are never touched.  Default OFF
    -> not called -> byte-identical.  Bare-atom-list XYZ (no header), format matches
    _flatten_sp2_atoms_xyz.  Deterministic (sorted metals).  Skips bridging donors (a single
    translation cannot satisfy two M-D ideals; leave those to the UFF pins), exactly like the pre-UFF
    snap."""
    if not RDKIT_AVAILABLE or mol_template is None or not xyz_delfin:
        return xyz_delfin
    try:
        import numpy as np
        lines = [l for l in xyz_delfin.splitlines() if l.strip()]
        if len(lines) != mol_template.GetNumAtoms():
            return xyz_delfin
        symbols: List[str] = []
        coords_list: List[List[float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return xyz_delfin
            symbols.append(parts[0])
            coords_list.append([float(parts[1]), float(parts[2]), float(parts[3])])
        coords = np.array(coords_list, dtype=float)
        metal_indices = {
            a.GetIdx() for a in mol_template.GetAtoms() if a.GetSymbol() in _METAL_SET
        }
        if not metal_indices:
            return xyz_delfin
        moved = False
        for m_idx in sorted(metal_indices):
            m_atom = mol_template.GetAtomWithIdx(m_idx)
            m_sym = m_atom.GetSymbol()
            m_pos = coords[m_idx]
            for nbr in m_atom.GetNeighbors():
                d_idx = nbr.GetIdx()
                d_sym = nbr.GetSymbol()
                if d_sym not in _METALLOID_MD_DONORS:
                    continue                       # metalloid donors ONLY -- N/O/P/S left to UFF
                if d_idx in metal_indices:
                    continue                       # M-M bonds
                if sum(1 for nb in nbr.GetNeighbors() if nb.GetIdx() in metal_indices) >= 2:
                    continue                       # bridging donor: leave to the UFF pins
                d_pos = coords[d_idx]
                cur_d = float(np.linalg.norm(d_pos - m_pos))
                if cur_d < 1e-8:
                    continue
                try:
                    # ⚠ BUG FIX 17.08.2026: here stood `_ml_bond_kind(mol, ...)`.
                    # In THIS body the parameter is called `mol_template`; a `mol` exists
                    # neither locally nor module-wide.  With ME_BOND_LEN=1 a NameError thus flew
                    # into the `except: continue` below -- and thereby skipped the
                    # METALLOID clamp ENTIRELY, although `METALLOID_MD_CLAMP` is a champion
                    # flag.  The new switch silently switched an old one off.  The same
                    # class as the `metal_idx`/`mi` error in `_manual_metal_embed`:
                    # a wrong name UNDER an `except` is invisible until one counts.
                    target_d = float(_get_ml_bond_length(m_sym, d_sym, _ml_bond_kind(mol_template, m_idx, d_idx) if _delfin_env_int("DELFIN_FFFREE_ME_BOND_LEN", 0) else "sigma"))
                except Exception:
                    continue
                if target_d <= 0 or abs(cur_d - target_d) / target_d < 0.05:
                    continue
                frag = _bfs_ligand_fragment(mol_template, d_idx, metal_indices)
                unit = (d_pos - m_pos) / cur_d
                delta = (m_pos + unit * target_d) - d_pos
                for fi in frag:
                    coords[fi] = coords[fi] + delta
                moved = True
        if not moved:
            return xyz_delfin
        out_lines = [
            f"{sym:4s} {coords[i][0]:12.6f} {coords[i][1]:12.6f} {coords[i][2]:12.6f}"
            for i, sym in enumerate(symbols)
        ]
        return "\n".join(out_lines) + "\n"
    except Exception as exc:
        logger.debug("Metalloid M-D clamp failed: %s", exc)
        return xyz_delfin


# d8 metals with a strong square-planar preference (mirror of _PREFERRED_CN4_GEOMETRY 'SQ' set).
_D8_SQ_METALS = frozenset({"Ni", "Pd", "Pt", "Au", "Rh", "Ir"})


# UNAMBIGUOUS d8 square-planar metals for the DELFIN_FFFREE_D8_SQ_ISO SP-4 imposition ONLY.  Rh and Ir
# are DROPPED here: their +1 state is d8 square-planar (Vaska/Wilkinson) but +3 is d6 OCTAHEDRAL, and a
# CN4-looking projection of a 5/6-coordinate Rh(III)/Ir(III) (e.g. a pincer HYDRIDE) has a real trans
# axis that fooled the isomer-preserving pass into forcing a square -> it collapsed distinct isomers
# (measured 2026-07-14: TIQWUN [IrH-2] pincer, n_isomers 2->1).  Pd(II)/Pt(II)/Ni(II)/Au(III) at CN4 are
# reliably SP-4, so restricting to them fixes TIQWUN while KEEPING every measured gain (all Pd/Pt).  The
# broader _D8_SQ_METALS (incl. Rh/Ir) is left untouched for the other, separate d8 paths.
_D8_SQ_ISO_METALS = frozenset({"Ni", "Pd", "Pt", "Au"})


def _flatten_d8_sq_planar_xyz(xyz_delfin: str, mol_template) -> str:
    """POST-UFF d8 square-planar flatten (root fix, env-gated DELFIN_FFFREE_D8_SQ_FLATTEN).

    A d8 CN4 centre (Pd/Pt/Ni/Au/Rh/Ir) MUST be square-planar: metal + 4 donors coplanar.  The ETKDG
    embed treats M-D like organic bonds and knows no square-planar preference, so it lifts the metal
    out of the donor plane -> SS-4 seesaw (CODSIA Pd 1.17 A out-of-plane; the chelate_oop_mer /
    poly_match SS-4->SP-4 defect, the #1 cluster).  A rigid rotation of the donors about M CANNOT fix
    this (M's out-of-plane distance is rotation-invariant); the donors must move RELATIVE to M.

    This projects each donor's M-D vector into the best-fit plane through M (v' = v - (v.n) n, rescaled
    to |v| so the M-D bond length is preserved exactly), and drags the donor's BFS fragment.  All 4
    donors then lie in one plane through M = square-planar.  For chelates the donors are backbone-
    linked, so an independent projection can strain the backbone; that strain shows up as broken
    frames and the per-system never-worse gate rejects it -- so this is SAFE to test broadly.
    Default OFF -> not called -> byte-identical.  Only fires on a genuine CN4 d8 centre.
    """
    if not RDKIT_AVAILABLE or mol_template is None or not xyz_delfin:
        return xyz_delfin
    try:
        import numpy as np
        lines = [l for l in xyz_delfin.splitlines() if l.strip()]
        if len(lines) != mol_template.GetNumAtoms():
            return xyz_delfin
        symbols: List[str] = []
        coords_list: List[List[float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return xyz_delfin
            symbols.append(parts[0])
            coords_list.append([float(parts[1]), float(parts[2]), float(parts[3])])
        coords = np.array(coords_list, dtype=float)
        metal_indices = {
            a.GetIdx() for a in mol_template.GetAtoms() if a.GetSymbol() in _METAL_SET
        }
        moved = False
        for m_idx in sorted(metal_indices):
            m_atom = mol_template.GetAtomWithIdx(m_idx)
            if m_atom.GetSymbol() not in _D8_SQ_METALS:
                continue
            donor_idxs = [nbr.GetIdx() for nbr in m_atom.GetNeighbors()
                          if nbr.GetIdx() not in metal_indices and nbr.GetSymbol() != "H"]
            if len(donor_idxs) != 4:           # square-planar flatten is defined for CN4 only
                continue
            m_pos = coords[m_idx]
            V = np.array([coords[d] - m_pos for d in donor_idxs], dtype=float)
            # best-fit plane THROUGH M (minimise sum((Di-M).n)^2): n = smallest right-singular vector.
            try:
                _u, _s, vh = np.linalg.svd(V, full_matrices=False)
            except Exception:
                continue
            n = vh[-1]
            nn = float(np.linalg.norm(n))
            if nn < 1e-9:
                continue
            n = n / nn
            # skip if already essentially planar (max out-of-plane below 0.30 A, the detector's cut)
            oop = float(np.max(np.abs(V @ n)))
            if oop <= 0.30:
                continue
            # BACKBONE-COUPLING GUARD: a chelate's donors share a backbone, so their BFS fragments
            # OVERLAP -- moving them independently would multiply-shift the shared atoms and wreck the
            # ligand.  Only flatten when the 4 donor fragments are pairwise DISJOINT (monodentate
            # donors).  Chelate d8 centres are left untouched (a rigid-body flatten is invariant; they
            # need a backbone-aware embed-level fix, not this).
            frags = [set(_bfs_ligand_fragment(mol_template, d, metal_indices)) for d in donor_idxs]
            _union = set()
            _overlap = False
            for fr in frags:
                if _union & fr:
                    _overlap = True
                    break
                _union |= fr
            if _overlap:
                continue
            # SQUARE REARRANGEMENT: place the 4 monodentate donors at ideal square-planar slots
            # (0/90/180/270 deg) in the plane through M, preserving each M-D length.  Coplanarity
            # alone leaves a tetrahedron coplanar-but-not-square; d8 CN4 must be SP-4 (90 deg spacing).
            # Runs PRE-UFF so the downstream freeze (_build_coordination_constraints_from_xyz freezes
            # monodentate donors at their positions) locks in the SQUARE, and UFF preserves it.
            vecs = [coords[d] - m_pos for d in donor_idxs]
            lens = [float(np.linalg.norm(v)) for v in vecs]
            if min(lens) < 1e-8:
                continue
            e1_raw = vecs[0] - float(np.dot(vecs[0], n)) * n
            if float(np.linalg.norm(e1_raw)) < 1e-8:
                continue
            e1 = e1_raw / float(np.linalg.norm(e1_raw))
            e2 = np.cross(n, e1)
            e2 = e2 / (float(np.linalg.norm(e2)) + 1e-15)
            import math as _math
            # BULKY-TRANS assignment: put the 2 largest donor fragments at a TRANS pair (0/180 deg,
            # 180 deg apart) and the 2 smallest at the other trans pair (90/270).  A naive angle-sort
            # can place two bulky ligands cis (90 deg) -> phenyl clash (real Pd(PPh3)2 is trans).
            # Order by fragment size (desc) and map to [0, 180, 90, 270] so big<->big and small<->small
            # are each trans.  Deterministic tie-break by donor index.
            order = sorted(range(4), key=lambda i: (-len(frags[i]), donor_idxs[i]))
            slots = [0.0, _math.pi, _math.pi / 2.0, 3.0 * _math.pi / 2.0]
            for slot_rank, di in enumerate(order):
                d = donor_idxs[di]
                phi = slots[slot_rank]
                v_ideal = lens[di] * (_math.cos(phi) * e1 + _math.sin(phi) * e2)
                delta = (m_pos + v_ideal) - coords[d]
                for fi in frags[di]:
                    coords[fi] = coords[fi] + delta
                moved = True
        if not moved:
            return xyz_delfin
        out_lines = [
            f"{sym:4s} {coords[i][0]:12.6f} {coords[i][1]:12.6f} {coords[i][2]:12.6f}"
            for i, sym in enumerate(symbols)
        ]
        return "\n".join(out_lines) + "\n"
    except Exception as exc:
        logger.debug("d8 square-planar flatten failed: %s", exc)
        return xyz_delfin


def _md_distance_in_tolerance(
    mol,
    conf_id: int,
    rel_low: float = 0.80,
    rel_high: float = 1.20,
) -> bool:
    """Return True iff every bonded M-D pair is within [rel_low, rel_high]
    times its ideal length ``_get_ml_bond_length(M, D)``.

    Heavy-atom donors only (skips H). Pure-graph: uses bonds present in the
    molecular graph. Universal — only element symbols + bond list.

    When ``DELFIN_PRE_UFF_TOPOLOGY_GATE_V2`` is enabled, the tolerance band
    is refined per-(M,D)-pair using empirical CSD-realistic ranges:
      - heavier donors (As/Sb/Sn/Te/Si/Ge/Pb/Bi) get a relaxed upper limit
        because the ML table is systematically biased low for these pairs;
      - donors involved in a multibond (double/triple, sum bond-order ≥ 1.5)
        to the metal get a relaxed lower limit because M=O / M=N / M≡N
        typically sit at 0.70-0.80 × single-bond ideal.
    The base [0.80, 1.20] band is bit-exact when V2 is OFF.
    """
    if not RDKIT_AVAILABLE:
        return True
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return True

    metal_indices = {
        a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET
    }
    if not metal_indices:
        return True

    _v2 = _pre_uff_topology_gate_v2_enabled()
    for m_idx in metal_indices:
        m_atom = mol.GetAtomWithIdx(m_idx)
        m_sym = m_atom.GetSymbol()
        mp = conf.GetAtomPosition(m_idx)
        for bond in m_atom.GetBonds():
            nbr = bond.GetOtherAtom(m_atom)
            d_idx = nbr.GetIdx()
            if d_idx in metal_indices:
                continue
            if nbr.GetAtomicNum() <= 1:
                continue
            d_sym = nbr.GetSymbol()
            dp = conf.GetAtomPosition(d_idx)
            d = math.sqrt(
                (mp.x - dp.x) ** 2 + (mp.y - dp.y) ** 2 + (mp.z - dp.z) ** 2
            )
            try:
                target_d = float(_get_ml_bond_length(m_sym, d_sym, _ml_bond_kind(mol, m_idx, d_idx) if _delfin_env_int("DELFIN_FFFREE_ME_BOND_LEN", 0) else "sigma"))
            except Exception:
                continue
            if target_d <= 0:
                continue
            lo, hi = rel_low, rel_high
            if _v2:
                try:
                    if bond.GetBondTypeAsDouble() >= 1.5:
                        lo = min(lo, _GATE_V2_DONOR_LOW_MULTIBOND.get(d_sym, lo))
                    hi = max(hi, _GATE_V2_DONOR_HIGH.get(d_sym, hi))
                except Exception:
                    pass
            rel = d / target_d
            if rel < lo or rel > hi:
                return False
    return True
