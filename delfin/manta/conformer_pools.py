"""Template conformer ranking, organic conformer pools, TFD dedup, ring puckers and d8 SP-4 variants of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
from typing import List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _xyz_to_rdkit_conformer,
)
from delfin.manta.converter_flags import (
    _PIPELINE_SEEDS,
    _delfin_env_int,
)
from delfin.manta.embed_timeout import (
    _embed_with_timeout,
)
from delfin.manta.geometry_quality import (
    _geometry_quality_score,
)
from delfin.manta.ligand_placement import (
    _snap_aromatic_rings_in_xyz,
)
from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    RDKIT_AVAILABLE,
    _METAL_SET,
)
from delfin.manta.mol_prep import (
    _prepare_mol_for_embedding,
)
from delfin.manta.openbabel_optimize import (
    _optimize_xyz_openbabel_safe,
)
from delfin.manta.topology_checks import (
    _has_atom_clash,
)

logger = get_logger("delfin.smiles_converter")


def _conf_ff_energy(mol, cid):
    """Cheap per-conformer force-field energy (MMFF, else UFF). Returns None when
    the molecule is unparameterisable (e.g. a transition metal) so callers fall
    back to geometry-only ranking. Pure-organic ligands (cyclohexane, sugars,
    macrocycle backbones) ARE parameterised -> their global-minimum pucker (the
    chair) gets a real energy and is kept."""
    try:
        from rdkit.Chem import AllChem
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=int(cid))
            if ff is not None:
                return float(ff.CalcEnergy())
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=int(cid))
        if ff is not None:
            return float(ff.CalcEnergy())
    except Exception:
        pass
    return None


def _conf_heavy_rmsd(mol, cid_a: int, cid_b: int) -> float:
    """Index-matched heavy-atom best-fit (Kabsch) RMSD between two conformers of
    the SAME molecule.  Read-only (does NOT mutate either conformer, unlike
    rdMolAlign.AlignMol), deterministic, and fast — atom indexing is identical
    across a molecule's conformers, so the index match IS the correct
    correspondence (no symmetry search needed for a basin-distinctness test).
    Returns a large value on any failure so the caller treats the pair as
    distinct (never silently merges basins)."""
    try:
        import numpy as np
        ca = mol.GetConformer(int(cid_a))
        cb = mol.GetConformer(int(cid_b))
        heavy = [at.GetIdx() for at in mol.GetAtoms() if at.GetAtomicNum() > 1]
        if len(heavy) < 3:
            heavy = [at.GetIdx() for at in mol.GetAtoms()]
        if len(heavy) < 3:
            return 1.0e9
        P = np.array([[ca.GetAtomPosition(i).x, ca.GetAtomPosition(i).y,
                       ca.GetAtomPosition(i).z] for i in heavy], float)
        Q = np.array([[cb.GetAtomPosition(i).x, cb.GetAtomPosition(i).y,
                       cb.GetAtomPosition(i).z] for i in heavy], float)
        P = P - P.mean(0)
        Q = Q - Q.mean(0)
        H = P.T @ Q
        U, _S, Vt = np.linalg.svd(H)
        d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        Pr = P @ R.T
        return float(np.sqrt(((Pr - Q) ** 2).sum() / len(heavy)))
    except Exception:
        return 1.0e9


def _rank_template_conformers(mol, *, top_k: Optional[int] = None) -> List[int]:
    """Return conformer IDs ranked deterministically by geometry quality.

    Lower ``_geometry_quality_score`` is better.  Conformers with an atom clash
    below ``min_dist=0.3`` are excluded (truly collapsed).  Ties are broken by
    ascending conformer ID so the ordering is fully reproducible.

    If ``top_k`` is given, only that many best candidates are returned.

    Completeness fix (DELFIN_FFFREE_CONF_ENERGY_RANK=1, default OFF -> byte-id):
    geometry-quality alone does NOT prefer the ENERGY global minimum, so a flexible
    ring's most important conformer (cyclohexane CHAIR) was dropped while higher-energy
    twist/half forms survived.  With the flag on, conformers are sorted by FF energy
    FIRST (then geometry score, then id), so the global-minimum conformer is always
    kept in top_k.  Falls back to geometry-only when the mol is unparameterisable
    (metal) -> byte-identical for those.

    Crystal-conformer coverage (DELFIN_FFFREE_CONF_BASIN_SPREAD=1, default OFF ->
    byte-id):  the top_k lowest-energy/best-geometry conformers are frequently
    near-duplicates of ONE basin, so the manifold never spans the conformational
    space the packing-distorted CRYSTAL conformer lives in.  With the flag on the
    selection is ADDITIVE (never-worse by construction): the original top_k frames
    are kept BYTE-IDENTICAL (same members, same order), then up to
    DELFIN_FFFREE_CONF_BASIN_EXTRA (default = top_k) ADDITIONAL energy-ordered
    RMSD-DISTINCT basins are APPENDED — each appended only if it is
    >= DELFIN_FFFREE_CONF_BASIN_RMSD (default 0.75 A) heavy-atom RMSD from every
    already-kept conformer.  The result is a strict SUPERSET of the baseline top_k
    (best-of-ensemble recall can only improve), and it is literally the WIDER
    conformer net — the higher-energy distinct minima that a crystal's packing
    forces select.  Rigid/small ligands with no further distinct basin add nothing
    (stay byte-identical).  The distinct basins are real (clash-filtered) minima,
    so this is principled widening, not spray (distinctness/validity guards stay
    clean)."""
    try:
        all_ids = [int(c.GetId()) for c in mol.GetConformers()]
    except Exception:
        return []
    if not all_ids:
        return []

    _energy_rank = os.environ.get("DELFIN_FFFREE_CONF_ENERGY_RANK", "0") == "1"
    ranked: List[Tuple[float, float, int]] = []
    for cid in all_ids:
        try:
            if _has_atom_clash(mol, cid, min_dist=0.3):
                continue
            score = float(_geometry_quality_score(mol, cid))
        except Exception:
            continue
        e = _conf_ff_energy(mol, cid) if _energy_rank else None
        # energy primary (rounded for determinism); None -> +inf so unparameterisable
        # mols keep the pure geometry-score ordering (byte-identical).
        ekey = round(e, 2) if e is not None else float("inf")
        ranked.append((ekey, score, cid))

    ranked.sort(key=lambda t: (t[0], t[1], t[2]))
    ordered = [cid for _e, _score, cid in ranked]

    _basin_spread = os.environ.get("DELFIN_FFFREE_CONF_BASIN_SPREAD", "0") == "1"
    if (_basin_spread and top_k is not None and top_k >= 1
            and len(ordered) > top_k):
        try:
            rmsd_min = float(os.environ.get("DELFIN_FFFREE_CONF_BASIN_RMSD", "0.75"))
        except ValueError:
            rmsd_min = 0.75
        try:
            n_extra = int(os.environ.get("DELFIN_FFFREE_CONF_BASIN_EXTRA", str(top_k)))
        except ValueError:
            n_extra = top_k
        # ADDITIVE: keep the baseline top_k unchanged (byte-id), then APPEND up to
        # n_extra higher-energy RMSD-distinct basins -> strict superset, wider net.
        kept: List[int] = list(ordered[:top_k])
        if n_extra > 0:
            added = 0
            for cid in ordered[top_k:]:
                if added >= n_extra:
                    break
                if all(_conf_heavy_rmsd(mol, cid, s) >= rmsd_min for s in kept):
                    kept.append(cid)
                    added += 1
        return kept

    if top_k is not None and top_k >= 0:
        ordered = ordered[:top_k]
    return ordered


def _frag_xyz_collapsed(mol, frag_list, frag_xyz, vol_min: float = 0.40) -> bool:
    """True if a fragment's TEMPLATE geometry has a planar-collapsed sp3 centre (tetra-volume < vol_min).

    Scoped version of ``_has_collapsed_sp3_centre`` operating on one fragment's numpy coords (``frag_xyz``
    rows in ``frag_list`` atom order).  Used by the ERDBEBEN isolated-fragment seating to decide whether to
    re-seat a fragment from a clean ISOLATED embed rather than the collapsed whole-complex template -- the
    decisive AQIBAE finding: the isolated ligand embeds 3D in 20/20 ETKDG seeds while the whole-complex
    ETKDG collapses one cage (a metal-context degeneracy, not a fragment-embed wall)."""
    try:
        import numpy as np
    except Exception:
        return False
    pos = {frag_list[i]: frag_xyz[i] for i in range(len(frag_list))}
    for k, ai in enumerate(frag_list):
        a = mol.GetAtomWithIdx(ai)
        if a.GetSymbol() in _METAL_SET or a.GetIsAromatic():
            continue
        if a.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        nb = [nbj.GetIdx() for nbj in a.GetNeighbors() if nbj.GetIdx() in pos]
        if len(nb) < 4 or sum(1 for j in nb if mol.GetAtomWithIdx(j).GetSymbol() != "H") < 2:
            continue
        pa = frag_xyz[k]
        near = sorted(nb, key=lambda j: float(np.sum((pos[j] - pa) ** 2)))[:4]
        p = [pos[j] for j in near]
        vol = abs(float(np.dot(np.cross(p[1] - p[0], p[2] - p[0]), p[3] - p[0]))) / 6.0
        if vol < vol_min:
            return True
    return False


def _heavy_atom_rmsd_xyz(xyz_a: str, xyz_b: str) -> float:
    """Symmetry-free heavy-atom RMSD (no alignment) between two XYZ blocks.

    Requires identical atom count and order (which holds when both XYZs
    are emitted from the same RDKit mol).  Skips hydrogens.  Returns a
    large value on any parse error so the caller conservatively treats
    the pair as distinct.
    """
    try:
        la = [l for l in xyz_a.strip().splitlines() if l.strip()]
        lb = [l for l in xyz_b.strip().splitlines() if l.strip()]
        if len(la) != len(lb) or not la:
            return 1e9
        import numpy as _np
        pa, pb = [], []
        for ra, rb in zip(la, lb):
            pa_parts = ra.split()
            pb_parts = rb.split()
            if len(pa_parts) < 4 or len(pb_parts) < 4:
                return 1e9
            if pa_parts[0] != pb_parts[0]:
                return 1e9
            if pa_parts[0].upper() == "H":
                continue
            pa.append([float(pa_parts[1]), float(pa_parts[2]), float(pa_parts[3])])
            pb.append([float(pb_parts[1]), float(pb_parts[2]), float(pb_parts[3])])
        if not pa:
            return 1e9
        A = _np.asarray(pa); B = _np.asarray(pb)
        A -= A.mean(axis=0); B -= B.mean(axis=0)
        H = A.T @ B
        U, _, Vt = _np.linalg.svd(H)
        d = _np.sign(_np.linalg.det(Vt.T @ U.T))
        D = _np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        Aaln = A @ R.T
        return float(_np.sqrt(_np.mean(_np.sum((Aaln - B) ** 2, axis=1))))
    except Exception:
        return 1e9


def _organic_mmff_ensemble(mol_h, *, max_pool: int = 8, rmsd_threshold: float = 0.5):
    """Complete low-energy conformer ENSEMBLE for an organic (DELFIN_FFFREE_CONF_ENERGY_RANK).

    Embed many ETKDGv3 seeds, MMFF-MINIMISE each so every conformer reaches its OWN
    local minimum (the chair stays a chair, the twist-boat stays a twist-boat — no
    collapse), then keep RMSD-DISTINCT basins in ASCENDING energy order: the global
    minimum is conf-1 and the OTHER genuine minima follow (twist-boat, boat, rotamers).
    This is the 'alternative to a global geometry optimisation' — the manifold contains
    the global minimum AND the populated conformers by construction.  MMFF-quality
    geometry, deterministic (fixed seed schedule, energy+id tie-break).  Returns
    [(xyz,label),...] or [] on failure (caller falls back to the legacy pool)."""
    try:
        from rdkit.Chem import AllChem
        m = Chem.Mol(mol_h)
        n = m.GetNumAtoms()
        n_emb = 250 if n <= 50 else (90 if n <= 90 else 30)
        p = AllChem.ETKDGv3()
        p.randomSeed = 42
        p.useRandomCoords = True
        p.enforceChirality = False
        p.pruneRmsThresh = -1.0
        cids = list(AllChem.EmbedMultipleConfs(m, numConfs=n_emb, params=p))
        if not cids:
            return []
        mp = AllChem.MMFFGetMoleculeProperties(m)
        ens: List[Tuple[float, int]] = []
        for cid in cids:
            try:
                ff = (AllChem.MMFFGetMoleculeForceField(m, mp, confId=int(cid))
                      if mp is not None else
                      AllChem.UFFGetMoleculeForceField(m, confId=int(cid)))
                if ff is None:
                    continue
                ff.Minimize(maxIts=1000)
                ens.append((float(ff.CalcEnergy()), int(cid)))
            except Exception:
                continue
        if not ens:
            return []
        ens.sort(key=lambda t: (round(t[0], 2), t[1]))     # energy, then id -> deterministic
        # symmetry-aware ALIGNED RMSD so equivalent conformers (ring-flip chairs,
        # rotated copies) collapse to ONE basin and only GENUINELY distinct minima
        # (chair vs twist-boat, anti vs gauche) are kept.
        from rdkit.Chem import rdMolAlign
        kept: List[int] = []
        for _e, cid in ens:
            if len(kept) >= max_pool:
                break
            dup = False
            for kcid in kept:
                try:
                    mc = Chem.Mol(m)                       # copy: don't mutate kept geometry
                    if rdMolAlign.GetBestRMS(mc, mc, prbId=int(cid), refId=int(kcid)) < rmsd_threshold:
                        dup = True
                        break
                except Exception:
                    pass
            if not dup:
                kept.append(cid)
        out: List[Tuple[str, str]] = []
        for j, cid in enumerate(kept):
            conf = m.GetConformer(int(cid))
            lines = []
            for i in range(n):
                a = m.GetAtomWithIdx(i)
                q = conf.GetAtomPosition(i)
                lines.append(f"{a.GetSymbol():4s} {q.x:12.6f} {q.y:12.6f} {q.z:12.6f}")
            out.append(("\n".join(lines) + "\n", f"conf-{j + 1}"))
        return out
    except Exception:
        return []


def _organic_conformer_pool(
    smiles: str,
    base_xyz: Optional[str],
    *,
    max_pool: int = 8,
    rmsd_threshold: float = 0.5,
    apply_uff: bool = True,
) -> List[Tuple[str, str]]:
    """Deterministic ETKDG pool for non-metal SMILES.

    Returns ``[(xyz, label), ...]`` where the first entry (if provided)
    is ``base_xyz`` with label ``"conf-1"``.  Extra conformers are
    embedded from seeded ETKDG parameters drawn from ``_PIPELINE_SEEDS``,
    UFF-polished, fused-ring-snapped, and accepted only when their
    heavy-atom RMSD exceeds ``rmsd_threshold`` against every previously
    accepted conformer.  Never raises.
    """
    if not RDKIT_AVAILABLE:
        return [(base_xyz, "")] if base_xyz else []
    pool: List[Tuple[str, str]] = []
    if base_xyz:
        pool.append((base_xyz, "conf-1"))
    try:
        mol = _prepare_mol_for_embedding(smiles, hapto_approx=False)
    except Exception:
        mol = None
    if mol is None:
        return pool
    try:
        mol_h = Chem.AddHs(mol)
    except Exception:
        mol_h = mol
    try:
        n_atoms = mol_h.GetNumAtoms()
    except Exception:
        return pool
    # Cap seed count by molecule size — huge fused-ring systems blow up
    # ETKDG wall-time; keep pool generation under a few seconds.
    if n_atoms > 120:
        seed_count = min(max_pool, 3)
    elif n_atoms > 80:
        seed_count = min(max_pool, 5)
    else:
        seed_count = max_pool
    # Completeness fix (DELFIN_FFFREE_CONF_ENERGY_RANK=1, default OFF -> byte-id):
    # the first-come RMSD-diversity selection does NOT keep the ENERGY global minimum
    # (cyclohexane's CHAIR was dropped for higher-energy twist forms).  With the flag
    # on we widen the seed pool, compute each conformer's FF energy, and accept the
    # RMSD-diverse survivors in ENERGY ORDER -> the global-minimum conformer (chair)
    # is always kept.  Default OFF keeps the exact legacy seed-order behaviour.
    _energy_rank = os.environ.get("DELFIN_FFFREE_CONF_ENERGY_RANK", "0") == "1"
    if _energy_rank:
        # Best-possible organic path: MMFF-minimised, energy-ranked, RMSD-distinct
        # basins (global minimum + all other genuine minima).  Falls through to the
        # legacy pool only if this fails.
        _ens = _organic_mmff_ensemble(mol_h, max_pool=max_pool, rmsd_threshold=rmsd_threshold)
        if _ens:
            return _ens
    if _energy_rank and n_atoms <= 80:
        seed_count = max(seed_count, min(40, len(_PIPELINE_SEEDS)))
    seeds = list(_PIPELINE_SEEDS[: max(1, seed_count)])
    candidates: List[Tuple[float, str]] = []     # (energy, xyz) when energy-ranking
    for seed in seeds:
        if not _energy_rank and len(pool) >= max_pool:
            break
        try:
            mol_try = Chem.Mol(mol_h)
            params = AllChem.ETKDGv3()
            params.randomSeed = int(seed)
            params.useRandomCoords = True
            params.enforceChirality = False
            cid = _embed_with_timeout(mol_try, params)
            if cid is None or cid < 0:
                continue
            e_key = _conf_ff_energy(mol_try, cid) if _energy_rank else None
            conf = mol_try.GetConformer(cid)
            lines = []
            for i in range(mol_try.GetNumAtoms()):
                atom = mol_try.GetAtomWithIdx(i)
                p = conf.GetAtomPosition(i)
                lines.append(
                    f"{atom.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}"
                )
            xyz = "\n".join(lines) + "\n"
            if apply_uff:
                try:
                    xyz = _optimize_xyz_openbabel_safe(xyz, mol_template=mol_h)
                except Exception:
                    pass
            try:
                xyz = _snap_aromatic_rings_in_xyz(xyz, mol_h, rms_threshold=0.05)
            except Exception:
                pass
            if not xyz or not xyz.strip():
                continue
            if _energy_rank:
                candidates.append((e_key if e_key is not None else float("inf"), xyz))
                continue
            if any(_heavy_atom_rmsd_xyz(xyz, x) < rmsd_threshold for x, _ in pool):
                continue
            pool.append((xyz, f"conf-{len(pool) + 1}"))
        except Exception as exc:
            logger.debug("Organic conformer seed %s failed: %s", seed, exc)
            continue
    if _energy_rank:
        # Discard the single-embed base anchor: it is one arbitrary ETKDG pose (often a
        # higher-energy twist) and would otherwise occupy conf-1 AND reject the true
        # global minimum as an RMSD "duplicate" (chair vs twist-boat are < threshold
        # apart).  The 40-seed candidate set is comprehensive, so rebuild the pool from
        # it in ASCENDING energy order -> the global minimum (chair) becomes conf-1 and
        # is always kept; then RMSD-diverse higher-energy minima follow.  Deterministic.
        if base_xyz:
            candidates.append((float("inf"), base_xyz))   # keep base only as a fallback
        candidates.sort(key=lambda t: t[0])
        pool = []
        for _e, xyz in candidates:
            if len(pool) >= max_pool:
                break
            if any(_heavy_atom_rmsd_xyz(xyz, x) < rmsd_threshold for x, _ in pool):
                continue
            pool.append((xyz, f"conf-{len(pool) + 1}"))
    pool = _append_ring_puckers(mol_h, pool, rmsd_threshold)
    if _delfin_env_int("DELFIN_TFD_DEDUP", 1):
        pool = _tfd_dedup_pool(mol_h, pool)
    return pool


def _tfd_dedup_pool(mol_h, pool, tfd_thr: float = 0.008):
    """Group-theoretic (symmetry-aware) dedup of the WHOLE conformer pool.

    ETKDG + rotamer sampling on a molecule with a symmetric group (a symmetric
    phosphine PCy3 / PPh3, a tert-butyl, three equivalent chelate arms, ...)
    rotates ONE bond many times and deposits the SAME rotamer over and over,
    because heavy-atom RMSD sees the relabelled-but-identical structures as
    distinct.  Torsion-Fingerprint-Deviation compares all ring + rotatable-bond
    torsions with the molecule's topological automorphisms folded in, so those
    symmetry-equivalent duplicates collapse to ONE while genuinely distinct
    conformers (a different pucker, a non-equivalent rotamer) survive — pucker
    and rotation deduped HOLISTICALLY under one symmetry-aware criterion.

    CRITICAL (user guardrail 2026-07-07): must NEVER drop a realistic,
    genuinely-distinct conformer — only TRUE duplicates.  Measured: symmetry-
    equivalent duplicates sit at TFD ~= 0.0 (TFD minimises over automorphisms),
    while genuinely distinct conformers start at TFD ~0.01-0.02 (ibuprofen 8->2
    at the 0.05 literature "same-cluster" threshold WRONGLY merged distinct
    rotamers).  So the threshold is deliberately TIGHT (0.008): it removes only
    the true relabelled duplicates and keeps every distinct realistic conformer.
    Frames whose atom order does not match ``mol_h`` are kept untouched.
    """
    if not RDKIT_AVAILABLE or len(pool) < 2:
        return pool
    try:
        from rdkit.Chem import TorsionFingerprints as _TF
    except Exception:
        return pool
    try:
        acc = Chem.Mol(mol_h)
        acc.RemoveAllConformers()
        loaded = []           # (pool_index, conf_id or None)
        for i, (xyz, _lbl) in enumerate(pool):
            conf = _xyz_to_rdkit_conformer(mol_h, xyz)
            if conf is None:
                loaded.append((i, None))
            else:
                loaded.append((i, acc.AddConformer(conf, assignId=True)))
        keep = [True] * len(pool)
        kept_ids = []
        for i, cid in loaded:
            if cid is None:
                continue          # unmatched order -> cannot compare, keep it
            dup = False
            for kid in kept_ids:
                try:
                    if _TF.GetTFDBetweenConformers(acc, [kid], [cid])[0] < tfd_thr:
                        dup = True
                        break
                except Exception:
                    pass
            if dup:
                keep[i] = False
            else:
                kept_ids.append(cid)
        deduped = [pool[i] for i in range(len(pool)) if keep[i]]
        # renumber conf-N labels sequentially, preserve non-"conf-" labels
        out = []
        for xyz, lbl in deduped:
            out.append((xyz, f"conf-{len(out) + 1}" if str(lbl).startswith("conf-") else lbl))
        return out
    except Exception:
        return pool


def _append_ring_puckers(mol_h, pool, rmsd_threshold):
    """Add explicitly-CONSTRUCTED ring-pucker conformers to an organic pool.

    ETKDG samples only the ground ring pucker (cyclohexane embeds 300/300 as the
    chair, never the twist-boat; a 5/7/8-ring never leaves its lowest pucker).
    Those higher ring basins are genuine, distinct, populated conformers a
    COMPLETE manifold must contain, but no seed count reaches them.  This pass
    constructs them for EVERY puckerable (saturated, size 5-8) ring via
    Cremer-Pople displacement + a torsion-held relax, and builds the multi-ring
    CARTESIAN PRODUCT (Cy3P: 3 cyclohexyls -> 3xchair, 2xchair+twist, ...) with a
    whole-molecule clash gate so the combined puckers stay sterically realistic.
    Deterministic, license-clean, TFD-deduped.  Byte-identical when
    ``DELFIN_RING_PUCKER=0`` or the molecule has no puckerable ring.
    """
    if not RDKIT_AVAILABLE or not pool:
        return pool
    if not _delfin_env_int("DELFIN_RING_PUCKER", 1):
        return pool
    try:
        from delfin.manta import _ring_pucker as _rpuck
        mol_base = Chem.Mol(mol_h)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol_base, params) != 0:
            return pool
        # The module already TFD-dedups internally (chair vs twist-boat differ
        # by only ~0.2 A heavy-atom RMSD, so the coarse pool RMSD gate would
        # wrongly reject the constructed puckers); append them directly, guarding
        # only against a near-exact coordinate duplicate.
        for _px, _plabel in _rpuck.generate(mol_base, budget=64):
            if any(_heavy_atom_rmsd_xyz(_px, x) < 0.05 for x, _ in pool):
                continue
            pool.append((_px, f"conf-{len(pool) + 1}"))
    except Exception as _rp_exc:
        try:
            logger.debug("ring-pucker construction skipped: %s", _rp_exc)
        except Exception:
            pass
    return pool


def _emit_ring_puckers_rp(mol, results, apply_uff, max_isomers):
    """Correct, combinatorial ring-pucker construction for the METAL path via
    ``delfin.manta._ring_pucker`` (Cremer-Pople displacement + torsion-held
    relax + whole-molecule clash gate).

    Covers BOTH ring kinds the user needs (rings with metal / rings without):
      * chelate rings that close THROUGH the metal — the metal and its
        coordinating donor atoms are FROZEN so the coordination sphere is
        preserved while the chelate backbone runs through its puckers (en
        delta/lambda, 6-ring chair/boat/twist, ...);
      * peripheral non-metal rings on the ligands (cyclohexyl, piperidinyl,
        sugar, ...), including their multi-ring CARTESIAN PRODUCT.

    Additive; every new frame still passes the final graph-topology gate.  Byte-
    identical when ``DELFIN_RING_PUCKER=0`` or no puckerable ring is present.
    """
    if not RDKIT_AVAILABLE or not results:
        return 0
    if not _delfin_env_int("DELFIN_RING_PUCKER", 1):
        return 0
    try:
        from delfin.manta import _ring_pucker as _rpuck
    except Exception:
        return 0
    base_xyz = results[0][0]
    # need an RDKit mol whose atom order matches the emitted XYZ (incl. H) so the
    # UFF relax + TFD are well-defined; try the mol as-is, then an H-added copy.
    conf = _xyz_to_rdkit_conformer(mol, base_xyz)
    work = mol
    if conf is None:
        try:
            work = Chem.AddHs(mol)
            conf = _xyz_to_rdkit_conformer(work, base_xyz)
        except Exception:
            conf = None
    if conf is None:
        return 0
    try:
        m = Chem.Mol(work)
        m.RemoveAllConformers()
        m.AddConformer(conf, assignId=True)
        frozen = set()
        for a in m.GetAtoms():
            if a.GetSymbol() in _METAL_SET:
                frozen.add(a.GetIdx())
                for nb in a.GetNeighbors():
                    frozen.add(nb.GetIdx())
        # sp3-C tetrahedral seating companion (DELFIN_FFFREE_SP3C_TET_SEAT): the metal + its donor atoms
        # are frozen above, but a monodentate sp3-C donor's HEAVY SUBSTITUENT stays FREE, so the pucker
        # relax swings it back to linear (M-C-X 180) even though the coordination sphere holds -- the
        # ring-pucker generator is the last construction path that re-linearises the seated sp3-C.  Freeze
        # each metal-bonded sp3-C donor's directly-bonded heavy substituent(s) too, locking the C-X vector
        # so the tetrahedral M-C-X survives puckering.  Element+graph only (universal); byte-identical off.
        if os.environ.get("DELFIN_FFFREE_SP3C_TET_SEAT", "0") == "1":
            _extra = set()
            for _fi in list(frozen):
                _fa = m.GetAtomWithIdx(_fi)
                # sp3 ALKYL donor only -- exclude sp2 carbene/aryl C donors (NHC, sigma-aryl): an aromatic
                # C is planar, has NO vacant tetrahedral slot, and freezing its ring neighbours distorts the
                # ring (eye's sp3c detector excludes sp2 via a pyramidality guard; mirror that here).
                if _fa.GetSymbol() != "C" or _fa.GetIsAromatic() or _fa.IsInRing():
                    continue        # PENDANT sp3 alkyl donor only (ring C: orientation fixed by scaffold)
                if not any(_nb.GetSymbol() in _METAL_SET for _nb in _fa.GetNeighbors()):
                    continue
                for _nb in _fa.GetNeighbors():
                    if _nb.GetSymbol() not in _METAL_SET and _nb.GetAtomicNum() > 1:
                        _extra.add(_nb.GetIdx())
            frozen |= _extra
        n_added = 0
        for _px, _plabel in _rpuck.generate(m, frozen=frozen, budget=48):
            if len(results) + n_added >= max_isomers:
                break
            results.append((_px, _plabel))
            n_added += 1
        return n_added
    except Exception:
        return 0


# ---------------------------------------------------------------------------
# Ring-conformer (pucker) enumeration for NON-METAL rings — the pucker
# analogue of Pólya coordination-isomer enumeration.
#
# Root cause: the generator emits one 3D structure per coordination-isomer /
# ETKDG seed; ring pucker is incidental to the random seed, and UFF cannot
# cross the ~10 kcal/mol chair<->boat barrier, so each peripheral non-metal
# ring stays frozen in whatever basin ETKDG happened to land in.  The metal
# pucker pass (``_emit_chelate_pucker_variants``) only touches rings that
# CONTAIN a metal atom; peripheral non-metal rings (e.g. the six cyclohexyls
# hanging off ZIGDOL's two As atoms) receive ZERO pucker enumeration.
#
# This sibling pass enumerates DISTINCT pucker basins (chair + boat) for the
# NON-METAL rings using the chemistry-accurate, graph-only template library
# in :mod:`delfin.manta._ring_conformer_templates` (which already excludes
# metal-chelate AND fully-aromatic rings), drives each variant cleanly into
# its target Cremer-Pople basin with a constrained geometric snap, and
# REJECTS any variant that does not reach its intended basin.  Pure additive
# (never reorders/drops existing frames); fixed deterministic ordering; no
# RNG.  Master flag DELFIN_RING_PUCKER_ENUM (default 1, validated net+ vs golden
# at 50k; set =0 to revert byte-identical to the no-op fall-through).
# ---------------------------------------------------------------------------


def _cp_basin_6ring(ring_coords) -> Tuple[str, float, float]:
    """Minimal in-delfin Cremer-Pople (1975) basin classifier for a 6-ring.

    *ring_coords* is a sequence of 6 (x, y, z) tuples IN RING ORDER.  Returns
    ``(basin, Q, theta_deg)`` where ``basin`` is ``"chair"`` / ``"boat"`` /
    ``"intermediate"`` / ``"planar"`` / ``"undefined"``.

    Replicates the q_m/phi_m formulation of Cremer & Pople, J. Am. Chem.
    Soc. 1975, 97, 1354 (the SAME math the private measurement gate uses):
    mean plane from the m=1 reference vectors, signed out-of-plane
    displacements z_j, then the (Q, theta) spherical coordinates from the
    q2 (boat/twist) and q3 (chair) puckering amplitudes.  Basin thresholds
    match the gate: chair <=> theta<=45 or >=135; boat/twist-boat <=>
    75<=theta<=105.  No import from quality_framework — standalone, ~30 LOC.
    """
    import numpy as _np
    pts = _np.asarray(ring_coords, dtype=float)
    if pts.shape[0] != 6:
        return ("undefined", 0.0, 0.0)
    R = pts - pts.mean(axis=0)
    j = _np.arange(6)
    s = _np.sin(2.0 * _np.pi * j / 6.0)
    c = _np.cos(2.0 * _np.pi * j / 6.0)
    Rp = (R * s[:, None]).sum(axis=0)
    Rpp = (R * c[:, None]).sum(axis=0)
    n = _np.cross(Rp, Rpp)
    nn = float(_np.linalg.norm(n))
    if nn < 1e-12:
        return ("undefined", 0.0, 0.0)
    n = n / nn
    z = R @ n  # signed out-of-mean-plane displacement per ring atom
    Q = float(_np.sqrt(float((z ** 2).sum())))
    # q2 (boat/twist amplitude, m=2) and q3 (chair amplitude, m=3=N/2 mode).
    c2 = _np.sqrt(2.0 / 6.0) * float((z * _np.cos(2.0 * _np.pi * 2 * j / 6.0)).sum())
    s2 = -_np.sqrt(2.0 / 6.0) * float((z * _np.sin(2.0 * _np.pi * 2 * j / 6.0)).sum())
    q2 = float(_np.sqrt(c2 * c2 + s2 * s2))
    q3 = (1.0 / _np.sqrt(6.0)) * float((z * _np.cos(_np.pi * j)).sum())  # signed
    if Q < 0.10:
        return ("planar", Q, 0.0)
    theta = math.degrees(math.atan2(q2, q3)) if (q2 or q3) else 0.0
    if theta < 0:
        theta += 360.0
    if theta > 180.0:
        theta = 360.0 - theta
    if theta <= 45.0 or theta >= 135.0:
        return ("chair", Q, theta)
    if 75.0 <= theta <= 105.0:
        return ("boat", Q, theta)
    return ("intermediate", Q, theta)


def _ring_canonical_snap_z(basin: str, size: int):
    """Return the canonical signed out-of-plane TARGET pattern (one value per
    ring atom, in ring order) that places a *size*-membered ring at the centre
    of *basin*.  Values are absolute multipliers of the snap amplitude (NOT
    L2-normalised) — the snap SETS each atom's out-of-plane component to
    ``pattern[k] * amplitude``.

    chair (6-ring): perfect D3d alternating +/- (Cremer-Pople theta -> 0/180,
                    the chair pole).
    boat  (6-ring): C2v two-flagpole pattern (atoms 0,3 up; 2,5 down; 1,4 in
                    the base plane) -> pure m=2 character, theta -> 90.
    """
    if size != 6:
        return None
    if basin == "chair":
        return [+1.0, -1.0, +1.0, -1.0, +1.0, -1.0]
    if basin == "boat":
        # Flagpoles at 0,3 (up), 2,5 (down), 1,4 in plane: cos(2*pi*2*j/6)
        # pattern = [1, -0.5, -0.5, 1, -0.5, -0.5] is the pure m=2 boat mode
        # (theta = 90).
        return [+1.0, -0.5, -0.5, +1.0, -0.5, -0.5]
    return None


def _emit_d8_sp4_variants(mol, results, apply_uff, max_isomers):
    """Additive d8-CN4 square-planar SEATING pass (eye-driven: Weddell flags "d8 CN4 built
    TETRAHEDRAL; must be square-planar" — HUPJUY Pd, VIHFAT Pt).  A d8 metal (Pt/Pd/Ni/Au/Rh/Ir)
    at CN4 is ALWAYS square-planar in the crystal; the chelate ETKDG embed can leave the donors
    TETRAHEDRAL and the SP-4 preference is silently lost (whole manifold -> ALARM).  For every
    emitted frame carrying a d8 metal, append a square-planarized variant (rigid per-ligand de-tilt
    of the donors into their common plane; ligand internals preserved).  The final graph-topology
    gate downstream drops any variant whose flattening introduced a ligand-ligand contact, so this
    is strictly ADDITIVE / never-worse.  Toggle DELFIN_D8_SP4_SEAT (default 0 -> byte-identical)."""
    import numpy as np
    if not _delfin_env_int("DELFIN_D8_SP4_SEAT", 0):
        return 0
    try:
        from delfin.manta._d8_square_planar import square_planarize_frame, _D8
    except Exception:
        return 0

    def _sig(xyz_str):
        try:
            return tuple(sorted(
                (p[0], round(float(p[1]), 2), round(float(p[2]), 2), round(float(p[3]), 2))
                for p in (ln.split() for ln in xyz_str.strip().splitlines())
                if len(p) >= 4 and p[0] not in ("H", "h")))
        except Exception:
            return tuple()

    def _parse(xyz_str):
        syms, P = [], []
        for ln in xyz_str.strip().splitlines()[2:]:
            t = ln.split()
            if len(t) >= 4:
                syms.append(t[0]); P.append([float(t[1]), float(t[2]), float(t[3])])
        return syms, np.array(P, float)

    seen = {_sig(x) for x, _ in results}
    added = 0
    for xyz, lbl in list(results):
        try:
            syms, P = _parse(xyz)
            if not any(s in _D8 for s in syms):
                continue
            sp = square_planarize_frame(syms, P)
            if sp is None:
                continue
            head = xyz.strip().splitlines()[:2]
            new_xyz = "\n".join(head + [f"{syms[i]} {sp[i][0]:.6f} {sp[i][1]:.6f} {sp[i][2]:.6f}"
                                        for i in range(len(syms))])
            sg = _sig(new_xyz)
            if sg in seen:
                continue
            seen.add(sg)
            results.append((new_xyz, (lbl or "") + " SP-4 square planar (d8 seat)"))
            added += 1
        except Exception:
            continue
    return added
