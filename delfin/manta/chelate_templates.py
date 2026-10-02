"""Chelate conformer candidates (with the sigma cap override), Procrustes fragment embedding, the from-scratch topology builder and the topology template molecule of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _xyz_to_rdkit_conformer,
)
from delfin.manta.converter_flags import (
    DELFIN_CHELATE_ACCEPT_DELTA,
    DELFIN_CHELATE_CAP_30,
    DELFIN_CHELATE_CAP_60,
    DELFIN_CHELATE_CAP_90,
    DELFIN_CHELATE_N_TRIALS,
    DELFIN_CHELATE_RANK_CLASS_ALPHA,
    DELFIN_CHELATE_REJECT_DELTA,
    _CHELATE_EMBED_TIMEOUT,
    _PIPELINE_SEEDS,
    _chelate_class_donor_penalty,
    _class_conditional_flag,
    _trace_seating,
)
from delfin.manta.embed_timeout import (
    _embed_with_timeout,
)
from delfin.manta.hapto_detect import (
    _classify_complex_class,
    mol_from_smiles_rdkit,
)
from delfin.manta.isomer_labels import (
    _TOPO_GEOMETRY_VECTORS,
    _find_bridging_donors,
)
from delfin.manta.ligand_placement import (
    _align_and_orient_ligands,
    _build_multimetal_scaffold,
)
from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    RDKIT_AVAILABLE,
    _METAL_SET,
    _get_ml_bond_length,
)
from delfin.manta.single_structure import (
    smiles_to_xyz_quick,
)

logger = get_logger("delfin.smiles_converter")


# Iter-8.4a: sigma chelate-cap restoration (forward-port from 123a130).
# When the sigma class loses topology %match because chelate-cap tightening
# reduced ETKDG trial counts (HEAD: 8/15/20/40 at >90/>60/>30/>20 atom
# fragments), restore the historical wider trial counts for the sigma class
# only.  d8 / d10 metals with crowded chelates need >12 trials to find the
# correct bite-angle conformer.  Default OFF: bit-exact HEAD when env-flag
# is unset.  Class-dispatched in smiles_to_xyz_isomers entry.
_SIGMA_CHELATE_CAPS_123A: Dict[str, int] = {
    "cap_90": 8,    # >90 atoms — keep HEAD cap
    "cap_60": 15,   # >60 atoms — keep HEAD cap
    "cap_30": 40,   # >30 atoms — restore champion (HEAD: 40 already)
    "cap_20": 40,   # >20 atoms — restore champion (HEAD: 20 → 40)
}


_ITER84_SIGMA_CAPS_OVERRIDE: Optional[Dict[str, int]] = None
"""Module-global override for ``_chelate_conformer_candidates`` cap values.
Set at the entry of ``smiles_to_xyz_isomers`` when both
``DELFIN_SIGMA_PORT_123A130_ITER8=1`` AND the parent mol classifies as
'sigma'.  ``None`` (default) means use HEAD baseline caps unchanged.
Process-safe under multiprocessing pool_evaluator (each worker is a
separate process); recursive smiles_to_xyz_isomers calls within the same
mol re-set to the same value (deterministic from class)."""


def _chelate_conformer_candidates(
    mol,
    frag_atom_indices,
    donor_atom_indices,
    target_bite,
    n_trials: int = DELFIN_CHELATE_N_TRIALS,
    accept_delta: float = DELFIN_CHELATE_ACCEPT_DELTA,
    reject_delta: float = DELFIN_CHELATE_REJECT_DELTA,
    max_candidates: int = 5,
):
    """Return a deterministic list of chelate conformer coordinates
    whose donor-donor distance pattern matches the polyhedron
    vertex-pair pattern, ordered by goodness of fit (best first).

    Each entry is a dict ``{original_atom_idx: (x, y, z)}``.  The list
    never exceeds ``max_candidates``; every returned entry has
    ``delta < reject_delta``.  For bidentate chelates ``target_bite``
    is a scalar (``|vertex_i - vertex_j|``); for polydentate ligands
    it is a ``(k, k)`` pairwise distance matrix.  Returning a list
    (instead of only the best) lets the caller emit multiple ring
    puckers / backbone conformations as distinct topology isomers,
    which is essential for macrocyclic and tridentate+ ligands whose
    natural backbone space contains several clash-free realisations
    compatible with the same platonic polyhedron.
    """
    if not RDKIT_AVAILABLE:
        return []
    try:
        import numpy as _np
    except Exception:
        return []

    frag_list = sorted(frag_atom_indices)
    if len(frag_list) < 3:
        return []
    old_to_new = {old: new for new, old in enumerate(frag_list)}
    donor_new = [old_to_new[d] for d in donor_atom_indices if d in old_to_new]
    if len(donor_new) < 2:
        return []

    # --- Class-aware chelate-rank gating ----------------------------------
    # Default OFF — bit-exact when ``DELFIN_CHELATE_RANK_CLASS_AWARE`` is
    # unset.  When enabled the per-conformer composite score becomes
    #   composite = delta + alpha * (element_weight + spread_factor * pucker)
    # where ``element_weight`` is constant for a given fragment (donor
    # set fixed) and ``pucker`` is the per-conformer donor-plane RMS
    # out-of-plane deviation (Å).  Sigma class rewards puckered backbones
    # (chair / boat tridentates), hapto class is neutral on pucker.
    # Implementation note: the element weight alone cannot reorder
    # conformers from the same call (constant); the per-conformer pucker
    # is what makes the secondary score actually re-rank.
    _class_aware_enabled = False
    _class_penalty_const = 0.0
    _class_spread_factor = 0.0
    try:
        if _class_conditional_flag(
            "DELFIN_CHELATE_RANK_CLASS_AWARE", mol
        ):
            _class_aware_enabled = True
            _cls = _classify_complex_class(mol)
            _class_penalty_const = _chelate_class_donor_penalty(
                mol, donor_atom_indices, _cls,
            )
            # Per-class pucker-diversity reward.  Negative value =
            # reward (lowers the composite when the conformer is more
            # out-of-plane).  Sigma class strongly rewards pucker variety
            # because mer/fac tridentate isomers differ exactly in the
            # backbone-plane deviation.  Hapto and multi_* classes get a
            # milder reward; no_metal disables entirely.
            _class_spread_factor = {
                "sigma":       -0.50,
                "hapto":       -0.10,
                "multi_sigma": -0.30,
                "multi_hapto": -0.10,
                "no_metal":     0.00,
            }.get(_cls, 0.0)
    except Exception:
        _class_aware_enabled = False
        _class_penalty_const = 0.0
        _class_spread_factor = 0.0

    is_pairwise_matrix = not isinstance(target_bite, (int, float))
    if is_pairwise_matrix:
        try:
            target_mat = _np.asarray(target_bite, dtype=float)
        except Exception:
            return []
        if target_mat.shape != (len(donor_new), len(donor_new)):
            return []

    rw = Chem.RWMol(Chem.Mol())
    for aidx in frag_list:
        atom = mol.GetAtomWithIdx(aidx)
        new_idx = rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
        rw.GetAtomWithIdx(new_idx).SetFormalCharge(atom.GetFormalCharge())
        rw.GetAtomWithIdx(new_idx).SetNoImplicit(True)
        rw.GetAtomWithIdx(new_idx).SetNumExplicitHs(atom.GetNumExplicitHs())
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in old_to_new and bj in old_to_new:
            rw.AddBond(old_to_new[bi], old_to_new[bj], bond.GetBondType())
    try:
        Chem.SanitizeMol(rw)
    except Exception:
        try:
            rw.UpdatePropertyCache(strict=False)
        except Exception:
            pass
    frag_mol = rw.GetMol()

    # Shared pipeline seed schedule.  The first 12 entries coincide with
    # the top-level conformer sampling seeds so the chelate conformer
    # choice is cache-coherent with the rest of the pipeline; the tail
    # extends the search for heavy macrocycles and tridentate+ ligands
    # whose native-bite matching requires a wider conformational sweep.
    _SEEDS = _PIPELINE_SEEDS

    def _try_embed(m, seed):
        p = AllChem.ETKDGv3()
        p.useRandomCoords = True
        p.randomSeed = int(seed)
        try:
            return _embed_with_timeout(m, p, timeout=_CHELATE_EMBED_TIMEOUT)
        except Exception:
            return -1

    def _preopt_fragment_conformer(_mol_obj, _cid) -> None:
        """MMFF94/UFF pre-optimisation of an isolated ligand conformer.

        DISABLED pending investigation of Ir(ppy)2(acac) all-cis regression.
        """
        return

    # Collect every accepted (delta, coords_map) pair and sort by fit.
    accepted: List[Tuple[float, Dict[int, tuple]]] = []
    fallback_mol = None

    # Scale number of trials with fragment size.  Huge polydentate
    # ligands (terpyridine-NMe2, salen-biphep phosphine backbone,
    # phos-terpy) have ETKDG wall-clocks measured in seconds per seed,
    # so the default 40 trials × 6 s timeout accumulates into minutes
    # of wall-time — and with orphan thread pile-up starves the
    # subprocess long before the first isomer is written.  Caps are
    # exposed as DELFIN_CHELATE_CAP_{60,90} so regressions where the
    # correct pose lives past the 8-/15-seed window can be recovered
    # by raising the cap for that run.
    _frag_n = frag_mol.GetNumAtoms()
    # Iter-8.4a: when the sigma chelate-cap port is active (module-global
    # ``_ITER84_SIGMA_CAPS_OVERRIDE`` set by smiles_to_xyz_isomers entry),
    # use the wider 123a130 caps instead of HEAD baselines.  ``None``
    # preserves HEAD bit-exactness.
    _iter84_caps = _ITER84_SIGMA_CAPS_OVERRIDE
    if _iter84_caps is not None:
        cap_90 = _iter84_caps["cap_90"]
        cap_60 = _iter84_caps["cap_60"]
        cap_30 = _iter84_caps["cap_30"]
        cap_20 = _iter84_caps["cap_20"]
    else:
        cap_90 = DELFIN_CHELATE_CAP_90
        cap_60 = DELFIN_CHELATE_CAP_60
        cap_30 = DELFIN_CHELATE_CAP_30
        cap_20 = int(os.environ.get('DELFIN_CHELATE_CAP_20', '20'))
    if _frag_n > 90:
        n_trials = min(n_trials, cap_90)
    elif _frag_n > 60:
        n_trials = min(n_trials, cap_60)
    elif _frag_n > 30:
        n_trials = min(n_trials, cap_30)
    elif _frag_n > 20:
        # Moderate-size ligands (terpyridines ~22 atoms, salen-biphep ~25):
        # default n_trials of 40 yields diminishing returns past ~20 trials.
        # Halving here recovers ~60s per such SMILES from blocking subprocess
        # without measurable loss in conformer diversity (env override:
        # DELFIN_CHELATE_CAP_20; Iter-8.4a sigma override widens to 40).
        n_trials = min(n_trials, cap_20)

    for seed in _SEEDS[:n_trials]:
        cid = _try_embed(frag_mol, seed)
        if cid < 0:
            if fallback_mol is None:
                try:
                    rw2 = Chem.RWMol(frag_mol)
                    for a in rw2.GetAtoms():
                        if a.GetIsAromatic():
                            a.SetIsAromatic(False)
                    for b in rw2.GetBonds():
                        if (
                            b.GetIsAromatic()
                            or b.GetBondType() == Chem.BondType.AROMATIC
                        ):
                            b.SetIsAromatic(False)
                            b.SetBondType(Chem.BondType.SINGLE)
                    fallback_mol = rw2.GetMol()
                    try:
                        fallback_mol.UpdatePropertyCache(strict=False)
                    except Exception:
                        pass
                except Exception:
                    fallback_mol = None
            if fallback_mol is None:
                continue
            cid = _try_embed(fallback_mol, seed)
            if cid < 0:
                continue
            used = fallback_mol
        else:
            used = frag_mol

        # Ligand-first: relax the isolated conformer with MMFF94/UFF
        # before evaluating its donor pattern.  The brief pre-opt
        # removes bonded-term strain baked into ETKDG's initial coords
        # and gives a chemically sensible ligand shape that DFT can
        # start from after placement.  On macrocycles this is the
        # difference between a clean low-energy pucker and a ring
        # with spurious kinks.
        _preopt_fragment_conformer(used, cid)

        conf = used.GetConformer(cid)
        donor_pts = _np.array([
            [
                conf.GetAtomPosition(idx).x,
                conf.GetAtomPosition(idx).y,
                conf.GetAtomPosition(idx).z,
            ]
            for idx in donor_new
        ])
        if is_pairwise_matrix:
            diffs = donor_pts[:, None, :] - donor_pts[None, :, :]
            d_mat = _np.linalg.norm(diffs, axis=-1)
            n = len(donor_new)
            n_pairs = n * (n - 1) / 2
            ss = float(_np.triu((d_mat - target_mat) ** 2, k=1).sum())
            delta = (ss / max(n_pairs, 1.0)) ** 0.5
        else:
            p0 = donor_pts[0]
            p1 = donor_pts[1]
            d_dd = float(_np.linalg.norm(p0 - p1))
            delta = abs(d_dd - float(target_bite))
        if delta >= reject_delta:
            continue
        coords_map = {
            old: (
                conf.GetAtomPosition(new).x,
                conf.GetAtomPosition(new).y,
                conf.GetAtomPosition(new).z,
            )
            for old, new in old_to_new.items()
        }

        # Per-conformer pucker score for class-aware re-rank.  Computed
        # only when the class-aware flag is on; otherwise the secondary
        # term is 0 and the sort is bit-exact with HEAD.  Pucker = RMS
        # out-of-plane distance of donor atoms from their best-fit
        # plane (n_donors >= 3) or 0 for bidentate (no plane to fit).
        _pucker = 0.0
        if _class_aware_enabled and len(donor_new) >= 3:
            try:
                _cen = donor_pts.mean(axis=0)
                _X = donor_pts - _cen
                # SVD of centered donor matrix; smallest singular vector
                # = plane normal.  Plane RMS = smallest singular value /
                # sqrt(n).
                _U, _S, _Vt = _np.linalg.svd(_X, full_matrices=False)
                if _S.size >= 1:
                    _pucker = float(_S[-1] / max(1.0, len(donor_new) ** 0.5))
            except Exception:
                _pucker = 0.0

        accepted.append((delta, coords_map, _pucker))
        # Early-exit once we have ``max_candidates`` "good" fits
        # (delta < accept_delta).  Every extra seed after this point
        # can only replace an already-good candidate with a slightly
        # better one — not worth the 6 s per-seed wall-time on the
        # heaviest ligands where every seed costs real time.
        good = sum(1 for d, _c, _p in accepted if d < accept_delta)
        if good >= max_candidates:
            break

    if not accepted:
        return []

    # Sort by fit quality (best first), cap at max_candidates.
    if _class_aware_enabled:
        _alpha = DELFIN_CHELATE_RANK_CLASS_ALPHA
        accepted.sort(
            key=lambda item: (
                item[0]
                + _alpha * (
                    _class_penalty_const
                    + _class_spread_factor * item[2]
                )
            )
        )
    else:
        accepted.sort(key=lambda item: item[0])
    return [c for _d, c, _p in accepted[:max_candidates]]


def _best_chelate_conformer_coords(
    mol,
    frag_atom_indices,
    donor_atom_indices,
    target_bite,
    n_trials: int = DELFIN_CHELATE_N_TRIALS,
    accept_delta: float = DELFIN_CHELATE_ACCEPT_DELTA,
    rank: int = 0,
    max_candidates: int = 5,
):
    """Return the ``rank``-th best chelate conformer (``rank=0`` is the best
    fit).  When ``rank`` exceeds the number of accepted candidates, falls
    back to the best available one so callers always receive a valid
    placement when at least one conformer fits the target bite.

    ``max_candidates`` limits how many conformers the underlying search
    retains; keep it >= the highest ``rank`` the caller intends to query.
    """
    cands = _chelate_conformer_candidates(
        mol,
        frag_atom_indices,
        donor_atom_indices,
        target_bite,
        n_trials=n_trials,
        accept_delta=accept_delta,
        max_candidates=max(max_candidates, rank + 1),
    )
    if not cands:
        return None
    if rank < 0 or rank >= len(cands):
        return cands[0]
    return cands[rank]


def _embed_fragment_procrustes(
    mol,
    metal_idx: int,
    frag_atom_indices: set,
    frag_donor_indices: List[int],
    target_positions: List[Tuple[float, float, float]],
    coords: List[Tuple[float, float, float]],
    chelate_rank: int = 0,
) -> bool:
    """Embed a ligand fragment via RDKit ETKDG, then Procrustes-align donors to targets.

    Modifies *coords* in-place for atoms in *frag_atom_indices*.
    Returns True on success, False on failure (caller should use BFS fallback).
    """
    try:
        import numpy as np
    except ImportError:
        return False

    if not frag_donor_indices or not frag_atom_indices:
        return False

    # Build a sub-molecule for the fragment (non-metal atoms only)
    frag_list = sorted(frag_atom_indices)
    old_to_new = {old: new for new, old in enumerate(frag_list)}

    rw = Chem.RWMol(Chem.Mol())
    for aidx in frag_list:
        atom = mol.GetAtomWithIdx(aidx)
        new_idx = rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
        rw.GetAtomWithIdx(new_idx).SetFormalCharge(atom.GetFormalCharge())
        rw.GetAtomWithIdx(new_idx).SetNoImplicit(True)
        rw.GetAtomWithIdx(new_idx).SetNumExplicitHs(atom.GetNumExplicitHs())

    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in old_to_new and bj in old_to_new:
            rw.AddBond(old_to_new[bi], old_to_new[bj], bond.GetBondType())

    try:
        Chem.SanitizeMol(rw)
    except Exception:
        try:
            rw.UpdatePropertyCache(strict=False)
        except Exception:
            pass

    frag_mol = rw.GetMol()

    donor_new_indices = [old_to_new[d] for d in frag_donor_indices if d in old_to_new]
    if not donor_new_indices:
        return False

    # Chelate: let the ligand find its own native backbone geometry
    # whose donor-donor distances match the polyhedron vertex pairs.
    # Bidentate uses a scalar bite; tri-/tetradentate uses the full
    # pairwise distance matrix.
    if len(donor_new_indices) >= 2 and len(target_positions) >= len(donor_new_indices):
        tp = np.asarray(target_positions[:len(donor_new_indices)], dtype=float)
        if len(donor_new_indices) == 2:
            target = float(np.linalg.norm(tp[0] - tp[1]))
        else:
            # full pairwise matrix, donor order matches frag_donor_indices
            diffs = tp[:, None, :] - tp[None, :, :]
            target = np.linalg.norm(diffs, axis=-1)
        coords_map = _best_chelate_conformer_coords(
            mol, frag_atom_indices, frag_donor_indices, target,
            rank=chelate_rank,
        )
        if coords_map is not None:
            frag_coords = np.array(
                [list(coords_map[old]) for old in frag_list],
                dtype=float,
            )
        else:
            frag_coords = None
    else:
        frag_coords = None

    if frag_coords is None:
        # Single-seed fragment ETKDG fallback (monodentate, higher-denticity,
        # or bidentate chelate for which the conformer search failed).
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.randomSeed = 42
        try:
            cid = _embed_with_timeout(frag_mol, params, timeout=_CHELATE_EMBED_TIMEOUT)
        except Exception:
            cid = -1
        if cid < 0:
            try:
                rw2 = Chem.RWMol(frag_mol)
                for atom in rw2.GetAtoms():
                    if atom.GetIsAromatic():
                        atom.SetIsAromatic(False)
                for bond in rw2.GetBonds():
                    if (
                        bond.GetIsAromatic()
                        or bond.GetBondType() == Chem.BondType.AROMATIC
                    ):
                        bond.SetIsAromatic(False)
                        bond.SetBondType(Chem.BondType.SINGLE)
                frag_mol2 = rw2.GetMol()
                try:
                    frag_mol2.UpdatePropertyCache(strict=False)
                except Exception:
                    pass
                cid2 = _embed_with_timeout(
                    frag_mol2, params, timeout=_CHELATE_EMBED_TIMEOUT
                )
                if cid2 < 0:
                    return False
                frag_mol = frag_mol2
                cid = cid2
            except Exception:
                return False
        frag_conf = frag_mol.GetConformer(cid)
        frag_coords = np.array([
            [frag_conf.GetAtomPosition(i).x,
             frag_conf.GetAtomPosition(i).y,
             frag_conf.GetAtomPosition(i).z]
            for i in range(frag_mol.GetNumAtoms())
        ])

    src = frag_coords[donor_new_indices]
    tgt = np.array(target_positions[:len(donor_new_indices)])

    if len(src) != len(tgt) or len(src) == 0:
        return False

    # Procrustes alignment: translate, then rotate
    src_center = src.mean(axis=0)
    tgt_center = tgt.mean(axis=0)
    src_centered = src - src_center
    tgt_centered = tgt - tgt_center

    if len(src) >= 2:
        # SVD for optimal rotation
        H = src_centered.T @ tgt_centered
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        sign_matrix = np.diag([1, 1, 1 if d > 0 else -1])
        R = Vt.T @ sign_matrix @ U.T
    else:
        # Single donor: align the donor's lone-pair direction (anti-bisector
        # of donor -> heavy-neighbour vectors) with the donor -> metal
        # direction.  This locks the ring plane in a chemically-correct
        # orientation (LP points at M) instead of the earlier heuristic
        # that only aligned "centroid -> donor" (which left the ring
        # plane under-constrained and forced the post-build orient step
        # to fix LP alignment at clash cost).
        d_new = donor_new_indices[0]
        d_atom = frag_mol.GetAtomWithIdx(d_new)
        # ROOT FIX (DELFIN_FFFREE_SP3C_TET_SEAT=1, default OFF -> byte-identical): a monodentate sp3-C
        # donor (M-CH2-R, M-CH3) has NO lone pair -- the metal occupies the 4th tetrahedral vertex.  The
        # HEAVY-ONLY bisector below leaves only the single heavy tail (R) for an M-CH2-R, so
        # src_dir = -unit(donor->R) and aligning it with donor->metal drives R ANTI to the metal ->
        # M-C-R ~180 deg (eye: sp3c_donor_linear; a top sigma_coord defect).  Including the H neighbours
        # for an sp3-C makes the bisector the true vacant-slot direction -> M-C-X ~109 deg.  Not a
        # post-hoc bend; sets the placement.
        _sp3c_tet = (os.environ.get("DELFIN_FFFREE_SP3C_TET_SEAT", "0") == "1"
                     and d_atom.GetSymbol() == "C" and not d_atom.GetIsAromatic()
                     and not d_atom.IsInRing())        # PENDANT sp3 alkyl donor only
        nbr_idx = [
            nb.GetIdx() for nb in d_atom.GetNeighbors()
            if nb.GetAtomicNum() > 1 or _sp3c_tet
        ]
        src_dir = None
        if nbr_idx:
            d_pos_src = frag_coords[d_new]
            lp = np.zeros(3)
            for ni in nbr_idx:
                v = frag_coords[ni] - d_pos_src
                vn = np.linalg.norm(v)
                if vn > 1e-8:
                    lp += v / vn
            lp_n = np.linalg.norm(lp)
            if lp_n > 1e-6:
                # Anti-bisector: away from neighbours, toward lone pair.
                src_dir = -lp / lp_n
        if src_dir is None:
            # Atomic donor (no ring/chain neighbours): fall back to the
            # original centroid heuristic.
            frag_center = frag_coords.mean(axis=0)
            raw = src[0] - frag_center
            rn = np.linalg.norm(raw)
            src_dir = raw / rn if rn > 1e-8 else np.array([1.0, 0.0, 0.0])
        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
            _trace_seating("procrustes_single_donor elem=%s n_heavy_nbrs=%d method=%s" % (
                d_atom.GetSymbol(), len(nbr_idx),
                "anti_bisector_HEAVY_ONLY" if nbr_idx else "centroid"))
        # Target LP direction: donor -> metal (metal is at origin).
        tgt_dir = -tgt[0]
        tn = np.linalg.norm(tgt_dir)
        if tn < 1e-8:
            R = np.eye(3)
        else:
            tgt_dir = tgt_dir / tn
            v = np.cross(src_dir, tgt_dir)
            c = float(np.dot(src_dir, tgt_dir))
            if np.linalg.norm(v) < 1e-8:
                R = np.eye(3) if c > 0 else -np.eye(3)
            else:
                vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                R = np.eye(3) + vx + vx @ vx / (1 + c)

    # Transform all fragment atoms
    transformed = (frag_coords - src_center) @ R.T + tgt_center

    # Write back to coords
    for new_idx, old_idx in enumerate(frag_list):
        coords[old_idx] = tuple(transformed[new_idx])

    return True


def _build_topology_xyz_from_scratch(
    mol,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    chelate_rank: int = 0,
    inflation_schedule: Tuple[float, ...] = (0.3, 0.5, 0.7, 0.9, 1.0),
) -> Optional[str]:
    """Balloon-inflate build that needs no ETKDG template.

    Produces a DELFIN XYZ by constructing the coordination sphere
    *from scratch* — the metal sits at the origin, non-bridging donors
    are placed directly on their polyhedron-vertex × ideal-M-D
    positions, bridging donors sit on the M-M axis at the compromise
    point, and each ligand fragment is MMFF-optimised *in isolation*
    then rigidly Procrustes-aligned so that its donor atoms land on
    those vertices.  The full structure is then grown radially in
    `inflation_schedule` steps, with ligand rotations per step to
    break inter-fragment clashes.

    Returns a DELFIN-format XYZ string, or ``None`` when any phase
    fails in a way the caller cannot recover from (unknown
    geometry, fragment embedding impossible, unresolvable clash).

    This is the preferred pre-UFF builder because it guarantees
    CSD-realistic M-D distances, uses no ETKDG-template bias, and
    produces one deterministic structure per (CF, perm,
    chelate_rank) triple.  ``_build_topology_xyz_from_template``
    remains as a fallback for systems where the isolated-fragment
    ETKDG fails to converge.
    """
    if not RDKIT_AVAILABLE:
        return None
    try:
        import numpy as np
    except Exception:
        return None

    vectors = _TOPO_GEOMETRY_VECTORS.get(geometry)
    if not vectors:
        return None

    try:
        n_atoms = mol.GetNumAtoms()
        coords: List[List[float]] = [[0.0, 0.0, 0.0] for _ in range(n_atoms)]
        placed: set = set()

        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        coords[metal_idx] = [0.0, 0.0, 0.0]
        placed.add(metal_idx)

        # Phase 0a — scaffold the second metal + bridging donors if the
        # input is bimetallic.  Both helpers already return
        # {atom_idx: (x, y, z)} dicts with ideal M-M + bridge distances.
        all_metal_idxs = [
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() in _METAL_SET
        ]
        bridging = _find_bridging_donors(mol) if len(all_metal_idxs) >= 2 else []
        scaffold: Optional[Dict[int, Tuple[float, float, float]]] = None
        if len(all_metal_idxs) >= 2 and bridging:
            try:
                # Ensure metal_idx is the first entry so its position is
                # the origin the builder assumes.
                ordered_metals = [metal_idx] + [
                    m for m in all_metal_idxs if m != metal_idx
                ]
                scaffold = _build_multimetal_scaffold(
                    mol, ordered_metals, bridging
                )
            except Exception as exc:
                logger.debug(
                    "Balloon scaffold failed (falling back to mono-metal): %s",
                    exc,
                )
                scaffold = None
        if scaffold:
            for a_idx, (x, y, z) in scaffold.items():
                coords[a_idx] = [x, y, z]
                placed.add(a_idx)

        # Phase 0b — donor target positions.  Non-bridging donors of
        # ``metal_idx`` get their polyhedron-vertex × ideal-M-D;
        # already-placed bridging donors stay where the scaffold put
        # them.  Other-metal non-bridging donors stay unplaced for now
        # — they are either placed by their own per-metal build pass
        # or remain as free atoms for the downstream BFS.
        donor_target_map: Dict[int, Tuple[float, float, float]] = {}
        bridging_atom_set = {d_idx for d_idx, _ in bridging}
        for pos_idx, donor_list_idx in enumerate(perm):
            donor_atom_idx = donor_atom_indices[donor_list_idx]
            if donor_atom_idx in placed and donor_atom_idx in bridging_atom_set:
                # Bridging donor — keep scaffold position as its target.
                donor_target_map[donor_atom_idx] = tuple(coords[donor_atom_idx])
                continue
            donor_sym = mol.GetAtomWithIdx(donor_atom_idx).GetSymbol()
            bl = float(_get_ml_bond_length(metal_sym, donor_sym))
            vx, vy, vz = vectors[pos_idx]
            mag = math.sqrt(vx * vx + vy * vy + vz * vz)
            if mag > 1e-8:
                vx = vx / mag * bl
                vy = vy / mag * bl
                vz = vz / mag * bl
            donor_target_map[donor_atom_idx] = (vx, vy, vz)
            coords[donor_atom_idx] = [vx, vy, vz]
            placed.add(donor_atom_idx)

        # Phase 1 + 2 — decompose non-metal atoms into ligand fragments
        # and dock each fragment onto its donor target.  Re-uses the
        # existing ``_embed_fragment_procrustes`` which already runs
        # ``_chelate_conformer_candidates`` + MMFF-ish embedding
        # internally for polydentate ligands and Procrustes-aligns
        # the donor atoms to their targets.
        non_metal = {
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() not in _METAL_SET
        }
        adj: Dict[int, set] = {i: set() for i in non_metal}
        for bond in mol.GetBonds():
            bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in non_metal and bj in non_metal:
                adj[bi].add(bj)
                adj[bj].add(bi)

        visited_frag: set = set()
        fragments: List[set] = []
        for start in sorted(non_metal):
            if start in visited_frag:
                continue
            frag: set = set()
            stack = [start]
            while stack:
                node = stack.pop()
                if node in visited_frag:
                    continue
                visited_frag.add(node)
                frag.add(node)
                for nbr in adj.get(node, ()):
                    if nbr not in visited_frag:
                        stack.append(nbr)
            fragments.append(frag)

        any_fragment_failed = False
        for frag in fragments:
            frag_donors = [d for d in donor_atom_indices if d in frag]
            if not frag_donors:
                # Non-coordinating fragment: leave for the downstream
                # BFS-VSEPR pass to place.  Happens for pure spectator
                # fragments (counter-ion, crystallographic solvent).
                continue
            tgt_positions = [
                donor_target_map[d] for d in frag_donors
                if d in donor_target_map
            ]
            if not tgt_positions:
                continue
            try:
                ok = _embed_fragment_procrustes(
                    mol, metal_idx, frag, frag_donors, tgt_positions, coords,
                    chelate_rank=chelate_rank,
                )
            except Exception as frag_exc:
                logger.debug(
                    "Balloon fragment embed failed "
                    "(size=%d, donors=%d): %s",
                    len(frag), len(frag_donors), frag_exc,
                )
                ok = False
            if ok:
                placed.update(frag)
            else:
                any_fragment_failed = True

        # If any coordinating fragment failed its Procrustes placement
        # the balloon build cannot produce a CSD-realistic XYZ for
        # this (CF, perm) triple.  Signal None so the caller falls back
        # to the rigid-template builder, which can often rescue these
        # cases using an already-embedded full-molecule conformer.
        if any_fragment_failed:
            return None

        # Phase 2b — BFS-VSEPR to place any remaining (non-coordinating)
        # atoms.  These are usually H atoms or free fragments that the
        # Procrustes path above skipped.
        bond_len_default = 1.4
        queue = list(placed)
        while queue:
            current = queue.pop(0)
            cx, cy, cz = coords[current]
            atom = mol.GetAtomWithIdx(current)
            unplaced_nbrs = [
                n.GetIdx() for n in atom.GetNeighbors()
                if n.GetIdx() not in placed
            ]
            n_unplaced = len(unplaced_nbrs)
            for k, nbr_idx in enumerate(unplaced_nbrs):
                dx, dy, dz = 0.0, 0.0, 0.0
                for other in atom.GetNeighbors():
                    oi = other.GetIdx()
                    if oi in placed and oi != nbr_idx:
                        ox, oy, oz = coords[oi]
                        dx += cx - ox
                        dy += cy - oy
                        dz += cz - oz
                mag_base = math.sqrt(dx * dx + dy * dy + dz * dz)
                if mag_base < 1e-8:
                    dx, dy, dz = 1.0 + 0.1 * k, 0.3 * k, 0.0
                    mag_base = math.sqrt(dx * dx + dy * dy + dz * dz)
                if n_unplaced > 1 and k > 0:
                    angle = 2 * math.pi * k / n_unplaced
                    ax, ay, az = dx / mag_base, dy / mag_base, dz / mag_base
                    if abs(ax) < 0.9:
                        px, py, pz = 1.0, 0.0, 0.0
                    else:
                        px, py, pz = 0.0, 1.0, 0.0
                    dot_pa = px * ax + py * ay + pz * az
                    px -= dot_pa * ax
                    py -= dot_pa * ay
                    pz -= dot_pa * az
                    pm = math.sqrt(px * px + py * py + pz * pz)
                    if pm > 1e-8:
                        px /= pm; py /= pm; pz /= pm
                    cos_a = math.cos(angle); sin_a = math.sin(angle)
                    dx2 = dx*cos_a + (ay*dz - az*dy)*sin_a + ax*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dy2 = dy*cos_a + (az*dx - ax*dz)*sin_a + ay*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dz2 = dz*cos_a + (ax*dy - ay*dx)*sin_a + az*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dx, dy, dz = dx2, dy2, dz2
                    mag_base = math.sqrt(dx * dx + dy * dy + dz * dz)
                    if mag_base < 1e-8:
                        mag_base = 1.0
                dx = dx / mag_base * bond_len_default
                dy = dy / mag_base * bond_len_default
                dz = dz / mag_base * bond_len_default
                coords[nbr_idx] = [cx + dx, cy + dy, cz + dz]
                placed.add(nbr_idx)
                queue.append(nbr_idx)

        # Phase 3 — ligand-orientation pass to break any remaining
        # inter-fragment clashes.  ``_align_and_orient_ligands``
        # rotates each monodentate / bidentate fragment around its
        # M-donor axis (donors stay fixed on the polyhedron) and picks
        # the angle that minimises non-donor clashes.  The inflation
        # schedule is implicit: the build above already places donors
        # at full M-D ideal distance, so one orientation pass at r=1.0
        # is equivalent to the terminal step of a five-step balloon.
        try:
            _align_and_orient_ligands(
                coords, mol, metal_idx, donor_atom_indices
            )
        except Exception as orient_exc:
            logger.debug(
                "Balloon orientation optimiser failed: %s", orient_exc
            )

        if os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
            try:
                _mp = np.array(coords[metal_idx], dtype=float)
                for _d in donor_atom_indices:
                    if mol.GetAtomWithIdx(_d).GetSymbol() != "C":
                        continue
                    _dp = np.array(coords[_d], dtype=float); _mc = _mp - _dp
                    _mcn = float(np.linalg.norm(_mc))
                    if _mcn < 1e-6:
                        continue
                    _mc /= _mcn
                    for _nb in mol.GetAtomWithIdx(_d).GetNeighbors():
                        if _nb.GetAtomicNum() <= 1:
                            continue
                        _cx = np.array(coords[_nb.GetIdx()], dtype=float) - _dp
                        _cxn = float(np.linalg.norm(_cx))
                        if 1.3 < _cxn < 1.9:
                            _ang = float(np.degrees(np.arccos(
                                max(-1.0, min(1.0, float(np.dot(_mc, _cx / _cxn)))))))
                            _trace_seating("SCRATCH_POST donor=%d M-C-Xheavy=%.0f" % (_d, _ang))
            except Exception:
                pass

        # Build XYZ string.
        lines: List[str] = []
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            x, y, z = coords[i]
            lines.append(f"{atom.GetSymbol():4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(lines) + "\n"
    except Exception as exc:
        logger.debug("_build_topology_xyz_from_scratch failed: %s", exc)
        return None


def _build_topology_template_mol(smiles: str):
    """Create a topology-template RDKit Mol with a mapped 3D conformer.

    Uses the quick conversion XYZ as the geometry source and reconstructs a
    matching RDKit molecule via ``mol_from_smiles_rdkit + AddHs`` so atom
    order/length stays consistent for conformer injection.
    """
    if not RDKIT_AVAILABLE:
        return None
    try:
        xyz, err = smiles_to_xyz_quick(smiles)
        if err or not xyz:
            return None
        mol, _note = mol_from_smiles_rdkit(smiles, allow_metal=True)
        if mol is None:
            return None
        try:
            mol = Chem.AddHs(mol, addCoords=True)
        except Exception:
            pass
        conf = _xyz_to_rdkit_conformer(mol, xyz)
        if conf is None:
            return None
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        return mol
    except Exception:
        return None
