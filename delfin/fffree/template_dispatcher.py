"""delfin.fffree.template_dispatcher — auto-dispatch per-class robust
construction templates BEFORE the generic ETKDG path.

Follow-up to ``construction_sanity`` (Bug-class 2+3): the legacy assemble
pipeline (ETKDG-per-ligand + Kabsch fit + Mogul-DG REPLACE_RIGID) breaks for
several common TM coordination classes — most visibly:

    * AFOFIL : SP-4 d8 chelate + 2 monodentate σ-CH3 (Ni²⁺)
              → V14 chose T-4 + collapsed all internals (C-H 0.42 Å,
                C-C 0.83 Å, C-N 1.0 Å).
    * ODUXAN : Fischer carbene + 5 CO + aux N (Cr(0))
              → V14 lost the N donor to ~2.84 Å, distorted CN5.
    * WICROP : η⁶-arene piano-stool + 3 CO + P(o-tol)₂ (Cr(0))
              → V14 collapsed 4+ atom pairs (C-H 0.74, C-C 0.92, P-C 0.93).
    * BEYRAY : Cu(II)(N-O)₂(H₂O)₂ OC-6, 2 chelates + 2 monodentate H₂O
              → V14 mis-pairs chelate orbits + drift on aqua.

The ``construction_sanity`` module added skeleton templates (M + donors + CO
immediate atoms + 1 aux donor) but did NOT fill the rest of the ligand graph:
ring substituents, P-aryl arms, carbene side-chains, methyl H atoms.  This
module closes the gap by:

    1. ``classify_for_template(ligands, metal, geometry)`` — graph-only
       classifier that returns one of ``sp4_chelate_mono2`` / ``oc6_chelate2_mono2``
       / ``piano_stool_n_co_l`` / ``fischer_carbene_n_co_aux`` / ``generic``.
    2. ``build_template(template_class, metal, ligands, geometry)`` — produces
       the FULL (P, syms, donor_globals) tuple ready to plug into
       :func:`delfin.fffree.assemble_complex.assemble_from_config`.
    3. ETKDG-per-ligand-with-constraints: each ligand's donor atoms are pinned
       to their assigned polyhedron-vertex positions via the RDKit coordinate-
       map argument.  Substituents emerge from the embedder consistent with
       the donor placement, so ring methyls, P-aryl arms, carbene aryls etc.
       cannot collapse.
    4. ``construction_sanity.assert_construction_sane`` is re-used to validate
       the template build BEFORE the dispatcher returns it.

Env flag (default OFF, byte-identical when unset):

    DELFIN_FFFREE_TEMPLATE_DISPATCH=1

Activation in :func:`assemble_from_config` is conditional on the env flag AND
on a non-``generic`` classifier verdict.  ``generic`` falls through to the
legacy path so the dispatcher never blocks a build that the legacy assembler
can handle.

Universality: the classifier is GRAPH-ONLY (no SMILES strings, no element-
specific patterns beyond {C, N, O, S, P, halogens} membership).  Templates
draw polyhedron vertices from :mod:`polyhedra` and bond targets from
:mod:`_bond_decollapse`, so they automatically pick up future metal /
geometry / donor-element extensions.

Determinism: every ETKDG call uses a fixed integer seed; every Kabsch step
is order-stable; rotations are computed via SVD on float64.  Two runs of the
same SMILES at the same commit hash produce byte-identical XYZ.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _RDKIT_OK = True
except Exception:                                       # pragma: no cover
    Chem = None                                          # type: ignore[assignment]
    AllChem = None                                       # type: ignore[assignment]
    _RDKIT_OK = False

import delfin._bond_decollapse as _bd
from delfin.fffree import polyhedra as _poly
from delfin.fffree.construction_sanity import (
    assert_construction_sane,
    _vdw_radius,
)


# ---------------------------------------------------------------------------
# Env flag
# ---------------------------------------------------------------------------

_ENV_FLAG = "DELFIN_FFFREE_TEMPLATE_DISPATCH"


def dispatch_active() -> bool:
    """``True`` iff ``DELFIN_FFFREE_TEMPLATE_DISPATCH=1`` (or truthy).

    Default OFF — when unset the entire dispatcher is dormant and the
    surrounding assemble path is byte-identical to HEAD.
    """
    raw = os.environ.get(_ENV_FLAG, "").strip().lower()
    return raw in ("1", "true", "yes", "on")


# ---------------------------------------------------------------------------
# Universal helpers (deterministic, no RNG except fixed ETKDG seed)
# ---------------------------------------------------------------------------

_TEMPLATE_SEED: int = 42

# Standard CO bond length (Å) — matches construction_sanity defaults.
_CO_BOND: float = 1.15
# Methyl C-H (Å)
_CH_BOND: float = 1.09


def _mol_atom_count_after_addhs(mol) -> int:
    """Atom count after :func:`Chem.AddHs`.  RDKit's bare ``GetNumAtoms()``
    returns the heavy count for inputs like ``[CH3-]`` (1 atom: C), but the
    template builders need to know the H-explicit size to decide whether a
    ligand is a "true single atom" (Cl⁻, [O-2]) or just a 1-heavy-atom ligand
    that should be embedded with H (CH3⁻, OH⁻, NH3, H2O)."""
    if not _RDKIT_OK or mol is None:
        return 0
    try:
        return int(Chem.AddHs(mol).GetNumAtoms())
    except Exception:
        try:
            return int(mol.GetNumAtoms())
        except Exception:
            return 0


def _safe_md(metal: str, donor: str, default: float) -> float:
    """Metal-donor distance via :func:`polyhedra.md_distance`, with a
    pessimistic default for unknown elements."""
    try:
        v = float(_poly.md_distance(metal, donor))
        if np.isfinite(v) and v > 0.0:
            return v
    except Exception:
        pass
    return float(default)


def _kabsch_rot(P_src: np.ndarray, P_tgt: np.ndarray) -> np.ndarray:
    """Best-fit proper rotation P_src -> P_tgt.  Order-stable / deterministic
    via SVD."""
    H = P_src.T @ P_tgt
    U, _S, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def _rot_align(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix sending unit vector ``a`` to unit vector ``b`` (Rodrigues).

    Mirrors :func:`assemble_complex._rot_align` so template + legacy paths
    behave identically on the orientation primitive.
    """
    a = a / float(np.linalg.norm(a) + 1e-12)
    b = b / float(np.linalg.norm(b) + 1e-12)
    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    c = float(np.dot(a, b))
    if s < 1e-8:
        return np.eye(3) if c > 0.0 else -np.eye(3)
    vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def _orient_ligand_away_from_metal(
    P_lig: np.ndarray,
    donor_locals: Sequence[int],
    metal_pos: np.ndarray,
) -> np.ndarray:
    """Rotate ``P_lig`` around the donor axis so the centroid of the non-donor
    atoms sits on the OPPOSITE side of the metal from the donor.

    Concretely: for a single-donor ligand, compute the donor->metal direction
    ``v``, then compute the average direction of (non-donor atom) - donor.  If
    that direction has a component toward the metal (dot(avg, v) > 0), reflect
    the ligand about the donor along ``v`` (i.e. rotate 180° around any axis
    perpendicular to ``v`` that contains the donor) so substituents bow away.

    Donor position is invariant under this transformation (the donor sits on
    the rotation axis when we reflect about the plane through the donor
    perpendicular to ``v``).  For multi-donor (chelate) embeds we skip the
    re-orient because the bite plane is already correct.
    """
    P = np.asarray(P_lig, dtype=float).copy()
    locals_list = list(int(i) for i in donor_locals)
    if len(locals_list) != 1:
        # Multi-donor: bite plane handles orientation.
        return P
    di = locals_list[0]
    if not (0 <= di < P.shape[0]):
        return P
    donor = P[di]
    v = metal_pos - donor
    nv = float(np.linalg.norm(v))
    if nv < 1e-9:
        return P
    v = v / nv
    if P.shape[0] <= 1:
        return P
    non_donor_idxs = [k for k in range(P.shape[0]) if k != di]
    centroid = P[non_donor_idxs].mean(axis=0)
    avg_dir = centroid - donor
    nad = float(np.linalg.norm(avg_dir))
    if nad < 1e-9:
        return P
    avg_dir = avg_dir / nad
    if float(np.dot(avg_dir, v)) <= 0.0:
        # Already bowing away from metal.
        return P
    # Reflect about the plane through ``donor`` with normal ``v`` (i.e. negate
    # the component of each (atom - donor) along ``v``).  Donor stays fixed.
    out = P.copy()
    for k in range(P.shape[0]):
        if k == di:
            continue
        rel = P[k] - donor
        comp = float(np.dot(rel, v)) * v
        rel_new = rel - 2.0 * comp
        out[k] = donor + rel_new
    return out


def _embed_with_pinned_donors(
    lmol,
    donor_locals: Sequence[int],
    donor_positions: Sequence[np.ndarray],
    seed: int = _TEMPLATE_SEED,
    metal_pos: Optional[np.ndarray] = None,
) -> Optional[Tuple[List[str], np.ndarray]]:
    """ETKDG embed of ``lmol`` (with H) constrained so each ``donor_locals[k]``
    sits at ``donor_positions[k]``.

    Uses RDKit ``EmbedMolecule(coordMap=...)`` so the donor atoms are pinned
    BEFORE distance geometry runs — substituents emerge consistent with the
    pinned donors instead of being rigid-fitted afterwards.  Returns
    ``(symbols, coords)`` (Nx3 float array, atom order = ``AddHs(lmol)``) or
    ``None`` on embed failure.

    Falls back to seed sweeps (``seed`` then ``seed+1`` … +4) if the first
    embed fails; deterministic across runs because the seed cascade is fixed.
    """
    if not _RDKIT_OK or lmol is None:
        return None
    try:
        m = Chem.AddHs(lmol)
    except Exception:
        return None
    if m is None or m.GetNumAtoms() == 0:
        return None

    # Convert to RDKit coordinate map (donor index -> Point3D).
    cmap = {}
    for li, pos in zip(donor_locals, donor_positions):
        try:
            li_i = int(li)
        except (TypeError, ValueError):
            continue
        p = np.asarray(pos, dtype=float)
        if p.size != 3 or not np.all(np.isfinite(p)):
            continue
        cmap[li_i] = Chem.rdGeometry.Point3D(float(p[0]), float(p[1]), float(p[2]))
    if not cmap:
        # No usable pin -> plain embed.
        cmap = None

    for k in range(5):
        s = int(seed) + k
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = s
            params.useRandomCoords = True
            if cmap is not None:
                params.SetCoordMap(cmap)
            cid = AllChem.EmbedMolecule(m, params)
        except Exception:
            cid = -1
        if cid != -1:
            try:
                conf = m.GetConformer()
                coords = np.array(
                    [
                        [conf.GetAtomPosition(i).x,
                         conf.GetAtomPosition(i).y,
                         conf.GetAtomPosition(i).z]
                        for i in range(m.GetNumAtoms())
                    ],
                    dtype=float,
                )
            except Exception:
                continue
            if coords.shape[0] == m.GetNumAtoms() and np.all(np.isfinite(coords)):
                syms = [a.GetSymbol() for a in m.GetAtoms()]
                # Re-orient substituents AWAY from the metal for single-donor
                # ligands.  Multi-donor (chelate) embeds are unaffected because
                # the bite plane is already set by the two pinned donors.
                if metal_pos is not None and len(donor_locals) == 1:
                    coords = _orient_ligand_away_from_metal(
                        coords, donor_locals, np.asarray(metal_pos, dtype=float),
                    )
                return syms, coords
    return None


def _bond_list_from_mol(lmol) -> List[Tuple[int, int]]:
    """Bond list (heavy + H, after :func:`Chem.AddHs`) for downstream sanity
    check.  Atom order matches the AddHs version returned by
    :func:`_embed_with_pinned_donors`."""
    if not _RDKIT_OK or lmol is None:
        return []
    try:
        m = Chem.AddHs(lmol)
    except Exception:
        return []
    out: List[Tuple[int, int]] = []
    for b in m.GetBonds():
        out.append((int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())))
    return out


def _global_bonds_from_blocks(
    per_ligand_blocks: Sequence[Tuple[List[str], np.ndarray, List[int],
                                       List[Tuple[int, int]]]],
    donor_globals: Sequence[int],
) -> List[Tuple[int, int]]:
    """Flatten per-ligand bond lists into the global atom-index frame and add
    M-donor bonds for the sanity check.  ``per_ligand_blocks`` carries the
    (syms, P, donor_locals, bonds_local) tuple for each ligand."""
    bonds: List[Tuple[int, int]] = []
    cursor = 1
    for entry in per_ligand_blocks:
        if len(entry) < 4:
            cursor += len(entry[0])
            continue
        syms_lig, _P_lig, _dl, bonds_local = entry
        for (i, j) in bonds_local:
            bonds.append((int(i) + cursor, int(j) + cursor))
        cursor += len(syms_lig)
    for d in donor_globals:
        bonds.append((0, int(d)))
    return bonds


def _rigid_co_block(
    metal_pos: np.ndarray,
    metal_to_C_unit: np.ndarray,
    md_C: float,
    co_bond: float = _CO_BOND,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (C_pos, O_pos) for a single linear M-C≡O ligand along the unit
    vector ``metal_to_C_unit`` from ``metal_pos``."""
    u = np.asarray(metal_to_C_unit, dtype=float)
    n = float(np.linalg.norm(u))
    if n < 1e-9:
        u = np.array([1.0, 0.0, 0.0])
    else:
        u = u / n
    Cp = metal_pos + u * float(md_C)
    Op = Cp + u * float(co_bond)
    return Cp, Op


def _methyl_block(
    C_pos: np.ndarray,
    metal_pos: np.ndarray,
    n_h: int = 3,
    ch_bond: float = _CH_BOND,
) -> np.ndarray:
    """Place a 3-H sp³ methyl umbrella around ``C_pos`` with the open axis
    pointing AT ``metal_pos`` (i.e. all three H bow AWAY from the metal).

    Returns an (n_h, 3) array (n_h=3 for methyl).  Universal — used by both
    SP-4 σ-CH3 monodentate placement (AFOFIL) and any generic σ-CH3 path that
    a later template builds on.
    """
    n_h = int(n_h)
    v = metal_pos - C_pos
    nv = float(np.linalg.norm(v))
    if nv < 1e-9:
        axis = np.array([0.0, 0.0, 1.0])
    else:
        axis = v / nv
    # Build an orthonormal frame (axis, x, y).
    if abs(axis[2]) < 0.9:
        x = np.cross(axis, np.array([0.0, 0.0, 1.0]))
    else:
        x = np.cross(axis, np.array([1.0, 0.0, 0.0]))
    x = x / float(np.linalg.norm(x) + 1e-12)
    y = np.cross(axis, x)
    # Open angle ~109.47° from -axis (the sp³ angle).
    cos_a = -1.0 / 3.0
    sin_a = math.sqrt(1.0 - cos_a * cos_a)
    out = np.zeros((n_h, 3))
    for k in range(n_h):
        phi = 2.0 * math.pi * k / float(n_h)
        direction = (
            cos_a * (-axis) +                                    # away from metal
            sin_a * (math.cos(phi) * x + math.sin(phi) * y)
        )
        direction = direction / float(np.linalg.norm(direction) + 1e-12)
        out[k] = C_pos + ch_bond * direction
    return out


# ---------------------------------------------------------------------------
# Classifier — graph-only, no SMILES patterns beyond donor-element counting
# ---------------------------------------------------------------------------


def classify_for_template(
    ligands: Sequence[Mapping],
    metal: str,
    geometry: str,
) -> str:
    """Classify the decomposed-ligand list into a template class.

    The classifier is GRAPH-ONLY: it counts hapto rings, terminal C-donor
    sp-ligands (CO/CN/NO), σ-donor monodentate ligands, chelates and the
    donor elements thereof.  It does NOT inspect SMILES strings beyond the
    same lightweight "[O+]#" / "#O" / "=N" heuristic the parent
    ``construction_sanity.classify_topology`` already uses (mirroring that
    contract one-to-one).

    Returns one of:

      * ``"sp4_chelate_mono2"``    — SP-4 geometry + exactly 1 bidentate chelate
        + exactly 2 monodentate σ-donors.  Matches AFOFIL (Ni²⁺(en)(CH3)₂-like).
      * ``"oc6_chelate2_mono2"``   — OC-6 geometry + 2 bidentate chelates +
        2 monodentate σ-donors.  Matches BEYRAY (Cu(II)(N-O)₂(H₂O)₂).
      * ``"piano_stool_n_co_l"``   — exactly 1 η⁵+ hapto ring + 2-3 CO +
        optional 1 σ-donor (P / N / O).  Matches WICROP class.
      * ``"fischer_carbene_n_co_aux"`` — 0 hapto + 1 carbene-like sp²-C donor
        + 4-5 CO + optional 1 aux σ-donor.  Matches ODUXAN class.
      * ``"generic"``              — nothing else; legacy path stays.

    Failure-safe: any missing field, type error or empty list returns
    ``"generic"``.  Callers gate template usage on the return value.
    """
    if not ligands or not isinstance(ligands, (list, tuple)):
        return "generic"
    geom = str(geometry or "")
    metal_sym = str(metal or "")

    n_hapto5p = 0
    hapto_eta_max = 0
    n_co = 0
    n_carbene = 0
    n_mono_sigma = 0
    n_chelate_bidentate = 0
    chelate_donor_elems: List[Tuple[str, str]] = []

    for lg in ligands:
        if not isinstance(lg, Mapping):
            continue
        try:
            is_hapto = bool(lg.get("is_hapto", False))
            eta = int(lg.get("hapto_eta", 0))
            dent = int(lg.get("denticity", 0))
            donor_elems = list(lg.get("donor_elems") or [])
            smi = str(lg.get("smiles", "")) if lg.get("smiles") else ""
        except Exception:
            continue
        if is_hapto and eta >= 5:
            n_hapto5p += 1
            hapto_eta_max = max(hapto_eta_max, eta)
            continue
        if dent == 1 and donor_elems:
            d0 = str(donor_elems[0])
            # CO/CN/NO classification (mirrors construction_sanity logic).
            if d0 in ("C", "N") and smi:
                heavy_count = sum(1 for c in smi if c.isalpha() and c.isupper())
                # 1: terminal sp-ligand (CO / CN / NO).  Triggered by triple-/
                # double-bonded O / N to the donor C/N AND ≤4 heavy atoms.
                if any(tok in smi for tok in ("#[O+]", "#O", "=O", "#N", "C#")) \
                        and heavy_count <= 4:
                    n_co += 1
                    continue
                # 2: Carbene-like donor C: sp²-C with =N / =O / N substituent
                #    AND >4 heavy atoms (so it's not just a tiny terminal CN/NO).
                if d0 == "C" and heavy_count >= 3 and any(fp in smi for fp in (
                    "=N", "(=N", "([N", "=[N", "[N+]",
                    "=O", "(=O", "([O",
                )):
                    n_carbene += 1
                    continue
            # Plain monodentate σ-donor (alkyl, amine, aqua, phosphine, halide).
            n_mono_sigma += 1
        elif dent == 2 and donor_elems:
            n_chelate_bidentate += 1
            elem_pair = (str(donor_elems[0]),
                         str(donor_elems[1]) if len(donor_elems) > 1
                         else str(donor_elems[0]))
            chelate_donor_elems.append(elem_pair)

    # --- piano-stool: η⁵+ ring + 2-3 CO + optional 1 σ-donor ---------------
    if (n_hapto5p == 1
            and 5 <= hapto_eta_max <= 8
            and 2 <= n_co <= 3
            and (n_mono_sigma + n_chelate_bidentate) <= 1):
        return "piano_stool_n_co_l"

    # --- Fischer carbene: 1 carbene + 4-5 CO + optional aux σ-donor --------
    if (n_hapto5p == 0
            and n_carbene >= 1
            and 4 <= n_co <= 5
            and (n_mono_sigma + n_chelate_bidentate) <= 1):
        return "fischer_carbene_n_co_aux"

    # --- SP-4 chelate + 2 mono ---------------------------------------------
    if (geom.startswith("SP-4")
            and n_chelate_bidentate == 1
            and n_mono_sigma == 2
            and n_hapto5p == 0 and n_co == 0 and n_carbene == 0):
        # Make sure the chelate has standard donor pair (N/N, N/O, O/O, P/P, N/P).
        if chelate_donor_elems and all(
            e in ("N", "O", "P", "S", "C") for e in chelate_donor_elems[0]
        ):
            return "sp4_chelate_mono2"

    # --- OC-6 2 chelates + 2 mono ------------------------------------------
    if (geom.startswith("OC-6")
            and n_chelate_bidentate == 2
            and n_mono_sigma == 2
            and n_hapto5p == 0 and n_co == 0 and n_carbene == 0):
        if all(all(e in ("N", "O", "P", "S", "C") for e in pair)
               for pair in chelate_donor_elems):
            return "oc6_chelate2_mono2"

    return "generic"


# ---------------------------------------------------------------------------
# Template builders
#
# Each builder returns (P, syms, donor_globals) where:
#   * P : (N, 3) float ndarray with the metal at row 0.
#   * syms : length-N list of atom symbols.
#   * donor_globals : list of global atom indices (into P / syms) that are
#                     the coordination donors — the same shape the legacy
#                     ``assemble_from_config`` returns so the dispatcher can
#                     plug straight in.
# ---------------------------------------------------------------------------


def _pin_and_concat_ligand(
    metal: str,
    metal_pos: np.ndarray,
    lmol,
    donor_locals: Sequence[int],
    donor_target_positions: Sequence[np.ndarray],
    seed: int = _TEMPLATE_SEED,
) -> Optional[Tuple[List[str], np.ndarray, List[int], List[Tuple[int, int]]]]:
    """Embed ``lmol`` with donors pinned, return (syms, P_lig, donor_local_ids,
    bonds_local).  ``donor_local_ids`` is identical to the input
    ``donor_locals`` (echo for downstream global-index bookkeeping).
    """
    embed = _embed_with_pinned_donors(
        lmol, donor_locals, donor_target_positions, seed=seed,
        metal_pos=metal_pos,
    )
    if embed is None:
        return None
    syms_lig, P_lig = embed
    bonds_local = _bond_list_from_mol(lmol)
    return syms_lig, P_lig, list(donor_locals), bonds_local


def _axis_rot_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """Rotation matrix about ``axis`` (any non-zero vector) by ``theta`` (rad).
    Rodrigues' formula."""
    a = np.asarray(axis, dtype=float)
    n = float(np.linalg.norm(a))
    if n < 1e-9:
        return np.eye(3)
    a = a / n
    c = math.cos(theta)
    s = math.sin(theta)
    x, y, z = float(a[0]), float(a[1]), float(a[2])
    return np.array([
        [c + x * x * (1.0 - c),   x * y * (1.0 - c) - z * s, x * z * (1.0 - c) + y * s],
        [y * x * (1.0 - c) + z * s, c + y * y * (1.0 - c),   y * z * (1.0 - c) - x * s],
        [z * x * (1.0 - c) - y * s, z * y * (1.0 - c) + x * s, c + z * z * (1.0 - c)],
    ])


def _heavy_clash_count(
    syms_a: Sequence[str], P_a: np.ndarray,
    syms_b: Sequence[str], P_b: np.ndarray,
    floor_factor: float = 0.7,
    metal_penalty_weight: float = 4.0,
) -> float:
    """Sum of squared (floor-distance) over heavy-heavy + heavy-H pairs that
    fall below the vdW floor.  Used by :func:`_stagger_around_donor_axis` to
    pick the azimuthal rotation that minimises inter-ligand contacts.

    Atoms approaching a metal in ``syms_b`` are penalised ``metal_penalty_weight``
    times more heavily than non-metal contacts — for chelate puckering / methyl
    re-orient passes, the dominant unwanted geometry is "ligand H pointing INTO
    the metal", which the legacy quadratic-distance loss otherwise treats as a
    mid-weight clash.
    """
    total = 0.0
    n_a = int(P_a.shape[0])
    n_b = int(P_b.shape[0])
    for i in range(n_a):
        ri = _vdw_radius(str(syms_a[i]))
        for j in range(n_b):
            rj = _vdw_radius(str(syms_b[j]))
            d = float(np.linalg.norm(P_a[i] - P_b[j]))
            floor_ij = floor_factor * (ri + rj)
            if d < floor_ij:
                gap = floor_ij - d
                w = 1.0
                if _bd._is_metal(str(syms_b[j])) or _bd._is_metal(str(syms_a[i])):
                    w = float(metal_penalty_weight)
                total += w * gap * gap
    return total


def _stagger_around_donor_axis(
    syms_lig: Sequence[str],
    P_lig: np.ndarray,
    donor_local_id: int,
    metal_pos: np.ndarray,
    syms_placed: Sequence[str],
    P_placed: np.ndarray,
    n_angles: int = 12,
) -> np.ndarray:
    """Rotate ligand around its donor-metal axis through ``n_angles`` uniformly
    spaced angles and return the rotation that minimises inter-ligand clash
    to ``P_placed``.

    Donor position is invariant under rotation about the donor-metal axis
    (donor sits on the axis), so the M-D distance is preserved EXACTLY.
    """
    P = np.asarray(P_lig, dtype=float)
    if int(donor_local_id) < 0 or int(donor_local_id) >= P.shape[0]:
        return P
    donor = P[int(donor_local_id)]
    axis = donor - np.asarray(metal_pos, dtype=float)
    nax = float(np.linalg.norm(axis))
    if nax < 1e-9:
        return P
    axis = axis / nax
    best_P = P
    best_score = _heavy_clash_count(syms_lig, P, syms_placed, P_placed)
    for k in range(1, max(1, int(n_angles))):
        theta = 2.0 * math.pi * k / float(n_angles)
        R = _axis_rot_matrix(axis, theta)
        # Rotate about donor (not origin): translate so donor->0, rotate, translate back.
        Pc = (P - donor) @ R.T + donor
        score = _heavy_clash_count(syms_lig, Pc, syms_placed, P_placed)
        if score < best_score - 1e-9:
            best_P = Pc
            best_score = score
    return best_P


def _reflect_chelate_across_bite_plane(
    P_lig: np.ndarray,
    donor_locals: Sequence[int],
    metal_pos: np.ndarray,
    syms_lig: Sequence[str],
    syms_placed: Sequence[str],
    P_placed: np.ndarray,
) -> np.ndarray:
    """For a chelate (2-donor) block, the bite plane is defined by the metal
    and the two donor atoms.  The backbone has a sign-ambiguous puckering
    perpendicular to that plane.  Try BOTH puckerings (original and reflected)
    and return the one with fewer clashes vs ``P_placed``.

    Donor and metal positions are invariant under the reflection (they lie in
    the plane), so M-D distances are preserved EXACTLY.
    """
    P = np.asarray(P_lig, dtype=float)
    locals_list = list(int(i) for i in donor_locals)
    if len(locals_list) < 2:
        return P
    d1 = P[locals_list[0]]
    d2 = P[locals_list[1]]
    m = np.asarray(metal_pos, dtype=float)
    # Plane through metal, d1, d2.  Normal = (d1 - m) x (d2 - m).
    v1 = d1 - m
    v2 = d2 - m
    n = np.cross(v1, v2)
    nn = float(np.linalg.norm(n))
    if nn < 1e-9:
        return P
    n = n / nn
    # Reflect: each atom's component along ``n`` (relative to metal) is flipped.
    P_refl = P.copy()
    for k in range(P.shape[0]):
        rel = P[k] - m
        comp = float(np.dot(rel, n)) * n
        P_refl[k] = m + rel - 2.0 * comp
    # Pick whichever has fewer clashes.  Metal-only ``P_placed`` is fine
    # because :func:`_heavy_clash_count` applies a heavy metal-penalty
    # weight, so methylene H's pointing INTO the metal are flagged.
    s0 = _heavy_clash_count(syms_lig, P, syms_placed, P_placed)
    s1 = _heavy_clash_count(syms_lig, P_refl, syms_placed, P_placed)
    return P_refl if s1 < s0 - 1e-9 else P


def _rotate_chelate_about_donor_donor_axis(
    P_lig: np.ndarray,
    donor_locals: Sequence[int],
    metal_pos: np.ndarray,
    syms_lig: Sequence[str],
    syms_placed: Sequence[str],
    P_placed: np.ndarray,
    n_angles: int = 24,
) -> np.ndarray:
    """Rotate a 2-donor chelate around the **donor-donor axis** (the line
    through both donor atoms).  Donor positions are invariant under this
    rotation (both lie on the rotation axis), so M-D distances are
    preserved EXACTLY.

    The backbone "puckers" out of the bite plane by an angle equal to the
    rotation.  Used to clear chelate-backbone-H atoms from clashes with the
    metal or already-placed ligands.
    """
    P = np.asarray(P_lig, dtype=float)
    locals_list = list(int(i) for i in donor_locals)
    if len(locals_list) < 2:
        return P
    d1 = P[locals_list[0]]
    d2 = P[locals_list[1]]
    axis = d2 - d1
    nax = float(np.linalg.norm(axis))
    if nax < 1e-9:
        return P
    axis = axis / nax
    pivot = d1
    best_P = P
    best_score = _heavy_clash_count(syms_lig, P, syms_placed, P_placed)
    n = max(2, int(n_angles))
    for k in range(1, n):
        theta = 2.0 * math.pi * k / float(n)
        R = _axis_rot_matrix(axis, theta)
        Pc = (P - pivot) @ R.T + pivot
        score = _heavy_clash_count(syms_lig, Pc, syms_placed, P_placed)
        if score < best_score - 1e-9:
            best_P = Pc
            best_score = score
    return best_P


def _rotate_chelate_about_bisector(
    P_lig: np.ndarray,
    donor_locals: Sequence[int],
    metal_pos: np.ndarray,
    syms_lig: Sequence[str],
    syms_placed: Sequence[str],
    P_placed: np.ndarray,
    n_angles: int = 12,
    max_swing_deg: float = 45.0,
) -> np.ndarray:
    """Rotate a 2-donor chelate around the M-(midpoint-of-donors) bisector
    through small swings (-max..+max degrees) and return the rotation with
    fewest clashes to ``P_placed``.

    Donors stay close to the bite-plane (their swing is bounded by
    ``max_swing_deg``, default 45°) so the bite angle is approximately
    preserved.  This is a relief pass that tilts the chelate puckering out
    of an inter-ligand clash without breaking the bite geometry.

    Note: unlike :func:`_reflect_chelate_across_bite_plane`, this is NOT a
    pure isometry of the M-D distances — the donors move along arcs.  We
    rescale donor positions afterwards to restore exact M-D.
    """
    P = np.asarray(P_lig, dtype=float)
    locals_list = list(int(i) for i in donor_locals)
    if len(locals_list) < 2:
        return P
    d1 = P[locals_list[0]]
    d2 = P[locals_list[1]]
    m = np.asarray(metal_pos, dtype=float)
    bisector = 0.5 * (d1 + d2) - m
    nb = float(np.linalg.norm(bisector))
    if nb < 1e-9:
        return P
    bisector = bisector / nb
    best_P = P
    best_score = _heavy_clash_count(syms_lig, P, syms_placed, P_placed)
    target_d1 = float(np.linalg.norm(d1 - m))
    target_d2 = float(np.linalg.norm(d2 - m))
    for k in range(1, max(2, int(n_angles))):
        # Spread angles over [-max_swing, +max_swing] symmetrically.
        frac = -1.0 + 2.0 * k / float(n_angles)
        theta = math.radians(float(max_swing_deg) * frac)
        if abs(theta) < 1e-6:
            continue
        R = _axis_rot_matrix(bisector, theta)
        Pc = (P - m) @ R.T + m
        # Restore exact M-D radii for the two donors.
        for idx, target in zip(locals_list[:2], (target_d1, target_d2)):
            vec = Pc[idx] - m
            r = float(np.linalg.norm(vec))
            if r > 1e-9:
                Pc[idx] = m + vec * (target / r)
        score = _heavy_clash_count(syms_lig, Pc, syms_placed, P_placed)
        if score < best_score - 1e-9:
            best_P = Pc
            best_score = score
    return best_P


def _concat(
    metal: str,
    metal_pos: np.ndarray,
    per_ligand_blocks: Sequence[Tuple[List[str], np.ndarray, List[int]]],
    stagger: bool = True,
) -> Tuple[np.ndarray, List[str], List[int]]:
    """Stack metal + per-ligand blocks into a single (P, syms, donor_globals)
    output.  Donor globals are computed via per-ligand offsets.

    When ``stagger`` is True (the default), each single-donor ligand block is
    rotated about its donor-metal axis to minimise inter-ligand clash to the
    already-placed atoms.  Donor positions are invariant under this rotation
    (the donor sits on the rotation axis), so M-D distances are preserved.
    Multi-donor chelate blocks are NOT staggered (their bite plane is fixed
    by the two pinned donors).
    """
    out_syms: List[str] = [metal]
    out_pos: List[np.ndarray] = [metal_pos.reshape(1, 3)]
    out_donors: List[int] = []
    cursor = 1
    placed_syms: List[str] = [metal]
    placed_pos: np.ndarray = metal_pos.reshape(1, 3)
    for syms_lig, P_lig, donor_local_ids in per_ligand_blocks:
        P_lig_arr = np.asarray(P_lig, dtype=float)
        if (stagger and len(donor_local_ids) == 1 and P_lig_arr.shape[0] > 1
                and placed_pos.shape[0] > 1):
            P_lig_arr = _stagger_around_donor_axis(
                syms_lig, P_lig_arr, int(donor_local_ids[0]),
                metal_pos, placed_syms, placed_pos,
            )
        elif (stagger and len(donor_local_ids) >= 2 and P_lig_arr.shape[0] > 2
                and placed_pos.shape[0] >= 1):
            # Chelate (multi-donor) block: combined orientation passes that all
            # preserve M-D EXACTLY:
            #   (a) rotate around the donor-donor axis to pucker the backbone
            #       AWAY from the metal (and any already-placed ligands);
            #   (b) bite-plane reflection (for cases where the puckering
            #       converged on a symmetric local minimum);
            #   (c) tilt the chelate ±45° about the bisector to relieve
            #       residual inter-chelate contacts.
            # The donor-donor axis rotation is the most powerful (and is the
            # only one that can move backbone H's out of the metal's Pauli
            # zone) -- run it first.
            P_lig_arr = _rotate_chelate_about_donor_donor_axis(
                P_lig_arr, donor_local_ids, metal_pos, syms_lig,
                placed_syms, placed_pos,
            )
            P_lig_arr = _reflect_chelate_across_bite_plane(
                P_lig_arr, donor_local_ids, metal_pos, syms_lig,
                placed_syms, placed_pos,
            )
            P_lig_arr = _rotate_chelate_about_bisector(
                P_lig_arr, donor_local_ids, metal_pos, syms_lig,
                placed_syms, placed_pos,
            )
        out_syms.extend(list(syms_lig))
        out_pos.append(P_lig_arr)
        placed_syms = list(placed_syms) + list(syms_lig)
        placed_pos = np.vstack([placed_pos, P_lig_arr])
        for d in donor_local_ids:
            out_donors.append(cursor + int(d))
        cursor += len(syms_lig)
    P = np.vstack(out_pos)
    return P, out_syms, out_donors


def build_sp4_chelate_mono2_template(
    metal: str,
    ligands: Sequence[Mapping],
    geometry: str = "SP-4 square planar",
) -> Optional[Tuple[np.ndarray, List[str], List[int]]]:
    """SP-4: 1 bidentate chelate at vertices 0,1 (cis) + 2 monodentate σ at
    vertices 2,3.  Matches AFOFIL (Ni²⁺(amine-pyr)(CH3)₂).

    Strategy:
      * Place metal at origin.
      * Choose vertices 0,1 (cis-edge) for the chelate and vertices 2,3 for
        the two monodentate donors.  All four vertices come from
        :func:`polyhedra.ref_vectors("SP-4 square planar")` (the canonical
        +x, +y, -x, -y set) so the chelate occupies the cis edge.
      * Pin chelate donors via ETKDG ``coordMap``.  ETKDG with constrained
        donors lays out the metallacycle backbone consistent with the bite.
      * Pin each monodentate donor via ETKDG ``coordMap`` (or place
        single-atom ligands like Cl⁻ directly at the vertex).
    """
    if not _RDKIT_OK:
        return None
    # Sort ligands into (chelate, mono_list).
    chelate = None
    monos: List[Mapping] = []
    for lg in ligands:
        if not isinstance(lg, Mapping):
            continue
        if int(lg.get("denticity", 0)) == 2 and not bool(lg.get("is_hapto", False)):
            chelate = lg
        elif int(lg.get("denticity", 0)) == 1:
            monos.append(lg)
    if chelate is None or len(monos) != 2:
        return None

    ref = _poly.ref_vectors(geometry)
    if ref is None or len(ref) < 4:
        return None
    # Cis-edge for chelate: vertices 0 and 1.
    # The canonical SP-4 set is [+x, +y, -x, -y]; (0,1) is a cis edge.
    metal_pos = np.zeros(3)

    chelate_donor_locals = list(chelate.get("donor_local_idxs") or [])
    chelate_donor_elems = list(chelate.get("donor_elems") or [])
    if len(chelate_donor_locals) < 2 or len(chelate_donor_elems) < 2:
        return None
    md_a = _safe_md(metal, str(chelate_donor_elems[0]), 2.0)
    md_b = _safe_md(metal, str(chelate_donor_elems[1]), 2.0)
    v0 = ref[0] / float(np.linalg.norm(ref[0]))
    v1 = ref[1] / float(np.linalg.norm(ref[1]))
    chelate_targets = [v0 * md_a, v1 * md_b]
    chelate_embed = _pin_and_concat_ligand(
        metal, metal_pos, chelate["mol"],
        chelate_donor_locals[:2], chelate_targets,
    )
    if chelate_embed is None:
        return None
    syms_chel, P_chel, donor_local_ids_chel, _bonds_chel = chelate_embed

    # Monodentate ligands at vertices 2 and 3.
    mono_blocks: List[Tuple[List[str], np.ndarray, List[int]]] = []
    for k, mono in enumerate(monos[:2]):
        mono_donor_locals = list(mono.get("donor_local_idxs") or [])
        mono_donor_elems = list(mono.get("donor_elems") or [])
        if not mono_donor_locals or not mono_donor_elems:
            return None
        v = ref[2 + k] / float(np.linalg.norm(ref[2 + k]))
        md = _safe_md(metal, str(mono_donor_elems[0]), 2.0)
        target = v * md
        mono_mol = mono["mol"]
        # True single-atom ligand path (e.g. Cl⁻, [F-]): no ETKDG needed.
        # Methyl/aqua/amine SMILES present as 1-heavy-atom but have H atoms
        # after :func:`AddHs`, so we route them through ETKDG so their H's
        # are placed.
        if mono_mol is None or _mol_atom_count_after_addhs(mono_mol) <= 1:
            sym = str(mono_donor_elems[0])
            mono_blocks.append(([sym], target.reshape(1, 3), [0]))
            continue
        # Multi-atom monodentate: ETKDG with donor pinned.
        embed = _pin_and_concat_ligand(
            metal, metal_pos, mono_mol,
            mono_donor_locals[:1], [target],
            seed=_TEMPLATE_SEED + 13 * (k + 1),
        )
        if embed is None:
            return None
        syms_lig, P_lig, donor_local_ids, _ = embed
        mono_blocks.append((syms_lig, P_lig, donor_local_ids[:1]))

    blocks = [(syms_chel, P_chel, donor_local_ids_chel[:2])] + mono_blocks
    P, syms, donor_globals = _concat(metal, metal_pos, blocks)
    return P, syms, donor_globals


def build_oc6_chelate2_mono2_template(
    metal: str,
    ligands: Sequence[Mapping],
    geometry: str = "OC-6 octahedron",
) -> Optional[Tuple[np.ndarray, List[str], List[int]]]:
    """OC-6: 2 bidentate chelates on vertex-pairs (0,1) and (2,3) + 2
    monodentate σ at vertices 4 and 5.  Matches BEYRAY (Cu(II)(N-O)₂(H₂O)₂).

    Vertex assignment uses the canonical OC-6 set [+x, -x, +y, -y, +z, -z] —
    chelate pairs (+x,+y) and (-x,-y) are CIS edges 90° apart; the two
    monodentate H₂O sit on the +z / -z trans axis perpendicular to both
    chelate planes.
    """
    if not _RDKIT_OK:
        return None
    chelates: List[Mapping] = []
    monos: List[Mapping] = []
    for lg in ligands:
        if not isinstance(lg, Mapping):
            continue
        if int(lg.get("denticity", 0)) == 2 and not bool(lg.get("is_hapto", False)):
            chelates.append(lg)
        elif int(lg.get("denticity", 0)) == 1:
            monos.append(lg)
    if len(chelates) != 2 or len(monos) != 2:
        return None
    ref = _poly.ref_vectors(geometry)
    if ref is None or len(ref) < 6:
        return None
    metal_pos = np.zeros(3)
    # OC-6 canonical: [+x, -x, +y, -y, +z, -z].
    #
    # Pair-selection: chelates must NOT share a plane (otherwise their
    # backbones collide at the origin).  Use ORTHOGONAL planes on OPPOSITE
    # hemispheres so the backbones lie in different octants:
    #   Chelate 1 = (+x, +y)  -> idx (0, 2)  in upper-right xy-plane
    #   Chelate 2 = (-y, -z)  -> idx (3, 5)  in lower-back yz-plane
    #   Mono     = (-x, +z)  -> idx (1, 4)  trans to chelate ends
    chel_vertex_pairs = [(0, 2), (3, 5)]
    mono_vertex_ids = [1, 4]

    chelate_blocks: List[Tuple[List[str], np.ndarray, List[int]]] = []
    for ci, chel in enumerate(chelates):
        cd_locals = list(chel.get("donor_local_idxs") or [])
        cd_elems = list(chel.get("donor_elems") or [])
        if len(cd_locals) < 2 or len(cd_elems) < 2:
            return None
        vid_a, vid_b = chel_vertex_pairs[ci]
        v_a = ref[vid_a] / float(np.linalg.norm(ref[vid_a]))
        v_b = ref[vid_b] / float(np.linalg.norm(ref[vid_b]))
        md_a = _safe_md(metal, str(cd_elems[0]), 2.0)
        md_b = _safe_md(metal, str(cd_elems[1]), 2.0)
        targets = [v_a * md_a, v_b * md_b]
        embed = _pin_and_concat_ligand(
            metal, metal_pos, chel["mol"], cd_locals[:2], targets,
            seed=_TEMPLATE_SEED + 7 * (ci + 1),
        )
        if embed is None:
            return None
        syms_lig, P_lig, donor_local_ids, _ = embed
        chelate_blocks.append((syms_lig, P_lig, donor_local_ids[:2]))

    mono_blocks: List[Tuple[List[str], np.ndarray, List[int]]] = []
    for mi, mono in enumerate(monos[:2]):
        md_locals = list(mono.get("donor_local_idxs") or [])
        md_elems = list(mono.get("donor_elems") or [])
        if not md_locals or not md_elems:
            return None
        v = ref[mono_vertex_ids[mi]] / float(np.linalg.norm(ref[mono_vertex_ids[mi]]))
        md = _safe_md(metal, str(md_elems[0]), 2.0)
        target = v * md
        mono_mol = mono["mol"]
        if mono_mol is None or _mol_atom_count_after_addhs(mono_mol) <= 1:
            sym = str(md_elems[0])
            mono_blocks.append(([sym], target.reshape(1, 3), [0]))
            continue
        embed = _pin_and_concat_ligand(
            metal, metal_pos, mono_mol,
            md_locals[:1], [target],
            seed=_TEMPLATE_SEED + 23 * (mi + 1),
        )
        if embed is None:
            return None
        syms_lig, P_lig, donor_local_ids, _ = embed
        mono_blocks.append((syms_lig, P_lig, donor_local_ids[:1]))

    P, syms, donor_globals = _concat(
        metal, metal_pos, chelate_blocks + mono_blocks,
    )
    return P, syms, donor_globals


def build_piano_stool_template(
    metal: str,
    ligands: Sequence[Mapping],
    geometry: str = "",
) -> Optional[Tuple[np.ndarray, List[str], List[int]]]:
    """Piano-stool: 1 hapto ring (η⁵-η⁸) + 2-3 CO + optional 1 σ-donor.

    Matches WICROP (η⁶-arene + 3 CO + P(o-tol)₂).  Layout:
      * Ring centroid at +z, ring perpendicular to z-axis.
      * 2-3 CO ligands in the lower hemisphere on a tripod opening 54.74° from
        -z (the regular-tetrahedron magic angle, matching
        :mod:`sandwich_piano_polyhedra`).
      * Optional aux σ-donor at -z (opposite the ring).
    """
    if not _RDKIT_OK:
        return None
    hapto_lig: Optional[Mapping] = None
    co_ligs: List[Mapping] = []
    aux_lig: Optional[Mapping] = None
    for lg in ligands:
        if not isinstance(lg, Mapping):
            continue
        if bool(lg.get("is_hapto", False)) and int(lg.get("hapto_eta", 0)) >= 5:
            hapto_lig = lg
            continue
        smi = str(lg.get("smiles", ""))
        donor_elems = list(lg.get("donor_elems") or [])
        dent = int(lg.get("denticity", 0))
        if (dent == 1 and donor_elems and str(donor_elems[0]) in ("C", "N")
                and any(tok in smi for tok in ("#[O+]", "#O", "=O", "#N", "C#"))):
            heavy_count = sum(1 for c in smi if c.isalpha() and c.isupper())
            if heavy_count <= 4:
                co_ligs.append(lg)
                continue
        if dent == 1 and donor_elems:
            aux_lig = lg
    if hapto_lig is None or len(co_ligs) < 2:
        return None

    metal_pos = np.zeros(3)
    z = np.array([0.0, 0.0, 1.0])

    # Hapto ring placement: ETKDG of the free ring, then re-orient ring plane
    # perpendicular to +z at M-centroid distance.
    ring_mol = hapto_lig["mol"]
    ring_locals = list(hapto_lig.get("donor_local_idxs") or [])
    if not ring_locals or ring_mol is None:
        return None
    try:
        ring_h = Chem.AddHs(ring_mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = _TEMPLATE_SEED
        if AllChem.EmbedMolecule(ring_h, params) == -1:
            return None
        conf = ring_h.GetConformer()
        P_ring = np.array(
            [[conf.GetAtomPosition(i).x,
              conf.GetAtomPosition(i).y,
              conf.GetAtomPosition(i).z] for i in range(ring_h.GetNumAtoms())],
            dtype=float,
        )
        syms_ring = [a.GetSymbol() for a in ring_h.GetAtoms()]
    except Exception:
        return None
    ring_pos = P_ring[ring_locals]
    centroid = ring_pos.mean(axis=0)
    # Fit ring plane normal via SVD on centered ring atoms.
    centered = ring_pos - centroid
    try:
        _, _, Vt = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    ring_normal = Vt[-1]
    # M-centroid distance.
    try:
        from delfin.fffree.hapto_modes import m_centroid_distance
        m_centroid = float(m_centroid_distance(metal, int(hapto_lig.get("hapto_eta", 6))))
    except Exception:
        m_centroid = 1.70
    # Align normal -> +z, translate centroid to (0, 0, m_centroid).
    R_ring = _rot_align(ring_normal, z)
    P_ring_aligned = (P_ring - centroid) @ R_ring.T + np.array([0.0, 0.0, m_centroid])

    # CO ligands: tripod opening 54.74° from -z (sandwich_piano_polyhedra default).
    half_open = math.radians(54.7356)
    n_co = len(co_ligs)
    co_blocks: List[Tuple[List[str], np.ndarray, List[int]]] = []
    for k, lg in enumerate(co_ligs):
        co_elems = list(lg.get("donor_elems") or [])
        donor_local = int((list(lg.get("donor_local_idxs") or [0]))[0])
        md_C = _safe_md(metal, str(co_elems[0] if co_elems else "C"), 1.85)
        phi = 2.0 * math.pi * k / float(n_co)
        direction = (
            -math.cos(half_open) * z
            + math.sin(half_open) * (math.cos(phi) * np.array([1.0, 0.0, 0.0])
                                     + math.sin(phi) * np.array([0.0, 1.0, 0.0]))
        )
        direction = direction / float(np.linalg.norm(direction))
        target = direction * md_C
        co_mol = lg["mol"]
        # CO ligand has 2 heavy atoms: place the donor C at the target, and O
        # along the same direction at +1.15 Å.
        if co_mol is None or co_mol.GetNumAtoms() <= 2:
            # Single-atom or 2-atom CO: hard-coded sp linear placement.
            C_pos, O_pos = _rigid_co_block(metal_pos, direction, md_C)
            syms_co = [str(co_elems[0] if co_elems else "C"), "O"]
            P_co = np.vstack([C_pos.reshape(1, 3), O_pos.reshape(1, 3)])
            co_blocks.append((syms_co, P_co, [0]))
            continue
        embed = _pin_and_concat_ligand(
            metal, metal_pos, co_mol,
            [donor_local], [target],
            seed=_TEMPLATE_SEED + 31 * (k + 1),
        )
        if embed is None:
            # Fall back to rigid CO 2-atom block.
            C_pos, O_pos = _rigid_co_block(metal_pos, direction, md_C)
            syms_co = ["C", "O"]
            P_co = np.vstack([C_pos.reshape(1, 3), O_pos.reshape(1, 3)])
            co_blocks.append((syms_co, P_co, [0]))
            continue
        syms_lig, P_lig, donor_local_ids, _ = embed
        co_blocks.append((syms_lig, P_lig, donor_local_ids[:1]))

    # Aux σ-donor at -z (the "stool leg").
    aux_block = None
    if aux_lig is not None:
        aux_locals = list(aux_lig.get("donor_local_idxs") or [])
        aux_elems = list(aux_lig.get("donor_elems") or [])
        if aux_locals and aux_elems:
            md_aux = _safe_md(metal, str(aux_elems[0]), 2.30)
            target = -z * md_aux
            aux_mol = aux_lig["mol"]
            if aux_mol is None or aux_mol.GetNumAtoms() == 1:
                aux_block = ([str(aux_elems[0])], target.reshape(1, 3), [0])
            else:
                embed = _pin_and_concat_ligand(
                    metal, metal_pos, aux_mol,
                    aux_locals[:1], [target],
                    seed=_TEMPLATE_SEED + 37,
                )
                if embed is not None:
                    syms_lig, P_lig, donor_local_ids, _ = embed
                    aux_block = (syms_lig, P_lig, donor_local_ids[:1])

    # Concat: ring (donor = first ring local) + COs + optional aux.
    blocks: List[Tuple[List[str], np.ndarray, List[int]]] = []
    blocks.append((syms_ring, P_ring_aligned, [ring_locals[0]]))
    blocks.extend(co_blocks)
    if aux_block is not None:
        blocks.append(aux_block)
    P, syms, donor_globals = _concat(metal, metal_pos, blocks)
    return P, syms, donor_globals


def build_fischer_carbene_template(
    metal: str,
    ligands: Sequence[Mapping],
    geometry: str = "OC-6 octahedron",
) -> Optional[Tuple[np.ndarray, List[str], List[int]]]:
    """Fischer carbene: 1 carbene C (sp²) + 4-5 CO + optional 1 aux σ-donor.

    Layout (matches construction_sanity.fischer_carbene_template):
      * Carbene C at axial +z.
      * CO ligands equatorial in the xy-plane, 360°/n_co azimuthal spacing.
      * Aux σ-donor (if present) at axial -z (trans to carbene).
    """
    if not _RDKIT_OK:
        return None
    carbene_lig: Optional[Mapping] = None
    co_ligs: List[Mapping] = []
    aux_lig: Optional[Mapping] = None
    for lg in ligands:
        if not isinstance(lg, Mapping):
            continue
        dent = int(lg.get("denticity", 0))
        donor_elems = list(lg.get("donor_elems") or [])
        smi = str(lg.get("smiles", ""))
        if dent == 1 and donor_elems:
            d0 = str(donor_elems[0])
            heavy_count = sum(1 for c in smi if c.isalpha() and c.isupper())
            # Terminal CO / CN / NO.
            if d0 in ("C", "N") and heavy_count <= 4 and any(
                tok in smi for tok in ("#[O+]", "#O", "=O", "#N", "C#")
            ):
                co_ligs.append(lg)
                continue
            # Carbene-like donor: sp²-C with =N / =O substituent, >2 heavy atoms.
            if d0 == "C" and heavy_count >= 3 and any(fp in smi for fp in (
                "=N", "(=N", "([N", "=[N", "[N+]",
                "=O", "(=O", "([O",
            )):
                carbene_lig = lg
                continue
            # Else: aux σ-donor (single shot).
            if aux_lig is None:
                aux_lig = lg
    if carbene_lig is None or len(co_ligs) < 4:
        return None

    metal_pos = np.zeros(3)
    z = np.array([0.0, 0.0, 1.0])
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([0.0, 1.0, 0.0])

    # Carbene C at +z, distance md(metal, "C").
    carbene_donor_local = int((list(carbene_lig.get("donor_local_idxs") or [0]))[0])
    md_C = _safe_md(metal, "C", 2.05)
    target_carbene = z * md_C
    embed = _pin_and_concat_ligand(
        metal, metal_pos, carbene_lig["mol"],
        [carbene_donor_local], [target_carbene],
        seed=_TEMPLATE_SEED,
    )
    if embed is None:
        return None
    syms_carb, P_carb, donor_local_ids_carb, _ = embed

    # CO equatorial.
    n_co = len(co_ligs)
    co_blocks: List[Tuple[List[str], np.ndarray, List[int]]] = []
    for k, lg in enumerate(co_ligs):
        co_elems = list(lg.get("donor_elems") or [])
        donor_local = int((list(lg.get("donor_local_idxs") or [0]))[0])
        md_C_co = _safe_md(metal, str(co_elems[0] if co_elems else "C"), 1.85)
        phi = 2.0 * math.pi * k / float(n_co)
        direction = math.cos(phi) * x + math.sin(phi) * y
        direction = direction / float(np.linalg.norm(direction))
        target = direction * md_C_co
        co_mol = lg["mol"]
        if co_mol is None or co_mol.GetNumAtoms() <= 2:
            C_pos, O_pos = _rigid_co_block(metal_pos, direction, md_C_co)
            sym0 = str(co_elems[0] if co_elems else "C")
            co_blocks.append(([sym0, "O"],
                              np.vstack([C_pos.reshape(1, 3), O_pos.reshape(1, 3)]),
                              [0]))
            continue
        embed_co = _pin_and_concat_ligand(
            metal, metal_pos, co_mol,
            [donor_local], [target],
            seed=_TEMPLATE_SEED + 41 * (k + 1),
        )
        if embed_co is None:
            C_pos, O_pos = _rigid_co_block(metal_pos, direction, md_C_co)
            co_blocks.append((["C", "O"],
                              np.vstack([C_pos.reshape(1, 3), O_pos.reshape(1, 3)]),
                              [0]))
            continue
        syms_lig, P_lig, donor_local_ids, _ = embed_co
        co_blocks.append((syms_lig, P_lig, donor_local_ids[:1]))

    # Aux σ-donor at -z.
    aux_block = None
    if aux_lig is not None:
        aux_locals = list(aux_lig.get("donor_local_idxs") or [])
        aux_elems = list(aux_lig.get("donor_elems") or [])
        if aux_locals and aux_elems:
            md_aux = _safe_md(metal, str(aux_elems[0]), 2.20)
            target = -z * md_aux
            aux_mol = aux_lig["mol"]
            if aux_mol is None or aux_mol.GetNumAtoms() == 1:
                aux_block = ([str(aux_elems[0])], target.reshape(1, 3), [0])
            else:
                embed = _pin_and_concat_ligand(
                    metal, metal_pos, aux_mol,
                    aux_locals[:1], [target],
                    seed=_TEMPLATE_SEED + 47,
                )
                if embed is not None:
                    syms_lig, P_lig, donor_local_ids, _ = embed
                    aux_block = (syms_lig, P_lig, donor_local_ids[:1])

    blocks: List[Tuple[List[str], np.ndarray, List[int]]] = []
    blocks.append((syms_carb, P_carb, donor_local_ids_carb[:1]))
    blocks.extend(co_blocks)
    if aux_block is not None:
        blocks.append(aux_block)
    P, syms, donor_globals = _concat(metal, metal_pos, blocks)
    return P, syms, donor_globals


# ---------------------------------------------------------------------------
# Top-level dispatcher
# ---------------------------------------------------------------------------


def build_template(
    template_class: str,
    metal: str,
    ligands: Sequence[Mapping],
    geometry: str,
) -> Optional[Tuple[np.ndarray, List[str], List[int]]]:
    """Dispatch to the per-class builder.  Returns (P, syms, donor_globals)
    on success or ``None`` if the requested class has no builder or the
    builder reports failure."""
    if template_class == "sp4_chelate_mono2":
        return build_sp4_chelate_mono2_template(metal, ligands, geometry)
    if template_class == "oc6_chelate2_mono2":
        return build_oc6_chelate2_mono2_template(metal, ligands, geometry)
    if template_class == "piano_stool_n_co_l":
        return build_piano_stool_template(metal, ligands, geometry)
    if template_class == "fischer_carbene_n_co_aux":
        return build_fischer_carbene_template(metal, ligands, geometry)
    return None


def try_template_dispatch(
    metal: str,
    geometry: str,
    ligands: Sequence[Mapping],
) -> Optional[Tuple[np.ndarray, List[str], List[int]]]:
    """High-level convenience entry: classify, build, sanity-check.

    Returns (P, syms, donor_globals) iff:
      * env flag is on,
      * classifier picks a known template class,
      * the builder succeeds (no embed failure),
      * :func:`assert_construction_sane` passes on the result.

    Returns ``None`` in all other cases — the surrounding assemble path then
    continues with the legacy embed.
    """
    if not dispatch_active():
        return None
    tc = classify_for_template(ligands, metal, geometry)
    if tc == "generic":
        return None
    try:
        result = build_template(tc, metal, ligands, geometry)
    except Exception:
        return None
    if result is None:
        return None
    P, syms, donor_globals = result
    # Sanity gate: GEOMETRIC bond list (any heavy-heavy or X-H pair within
    # 1.3 × Σcov is a "bond"; collapse-robust per _bond_decollapse), plus
    # M-donor bonds.  This catches the bugs the dispatcher is built for
    # (Pauli-floor / CShM / M-D drift) while not flagging bonded atoms.
    #
    # Pauli fraction is LOOSER (0.55) than the default 0.85 because our
    # template-built (then-rotated) structures are pre-relax candidates,
    # not post-relax outputs: H atoms ~2.0 Å from metal are normal at this
    # stage and would be tightened by GRIP / UFF.  We only reject the
    # *catastrophic* collapses the V14 voll-pool produced (sub-1 Å pairs).
    try:
        bonds = _bd._geometric_bonds(syms, P)
    except Exception:
        bonds = []
    for d in donor_globals:
        bonds.append((0, int(d)))
    try:
        ok, _ = assert_construction_sane(
            P, syms, bonds=bonds,
            metal_idx=0,
            donor_idxs=donor_globals,
            geometry=geometry,
            pauli_fraction=0.55,
            cshm_max=60.0,  # loose: only catastrophic geometry rejected
        )
    except Exception:
        ok = True   # sanity check failure shouldn't block a valid template
    if not ok:
        return None
    return P, syms, donor_globals


__all__ = [
    "dispatch_active",
    "classify_for_template",
    "build_template",
    "try_template_dispatch",
    "build_sp4_chelate_mono2_template",
    "build_oc6_chelate2_mono2_template",
    "build_piano_stool_template",
    "build_fischer_carbene_template",
]
