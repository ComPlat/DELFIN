"""Molecule preparation for embedding (cached), metal-donor distance rescale and robust multi-conformer embedding of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import re
from collections import OrderedDict
from typing import List, Tuple

from delfin.common.logging import get_logger
from delfin.manta.embed_timeout import (
    _embed_multiple_confs_with_timeout,
    _make_random_embed_params,
)
from delfin.manta.hapto_detect import (
    _apply_hapto_approximation,
    _find_hapto_groups,
    contains_metal,
    mol_from_smiles_rdkit,
)
from delfin.manta.metal_smiles import (
    _convert_metal_bonds_to_dative,
    _denormalize_metal_smiles,
    _fix_hapto_donor_h,
    _fix_organometallic_carbon_h,
    _hapto_h_always,
    _normalize_metal_smiles,
    _strip_h_on_metal_halogen,
)
from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    RDKIT_AVAILABLE,
    STK_AVAILABLE,
    _METAL_SET,
    _get_ml_bond_length,
    _is_metal_nitrogen_complex,
    _is_simple_organometallic,
    stk,
)

logger = get_logger("delfin.smiles_converter")


def is_smiles_string(content: str) -> bool:
    """Detect if input content is a SMILES string.

    Heuristics:
    - Single line or very few lines
    - Contains typical SMILES characters (parentheses, =, #, [, ])
    - No coordinate-like patterns (no multiple whitespace-separated numbers)
    - May contain '>' for coordination bonds in metal complexes

    Args:
        content: File content to check

    Returns:
        True if content appears to be a SMILES string
    """
    lines = [line.strip() for line in content.strip().split('\n') if line.strip()]

    if len(lines) == 0:
        return False

    # A filename (single name.ext line) is not a SMILES. The species pointer
    # "input.xyz" inside an OCCUPIER job's input.txt was read here as a
    # SMILES (the extension dot counts as the fragment separator and nearly
    # every element letter is organic), so the job demanded a SMILES
    # converter it did not need and died before reading the geometry
    # (calc 81-105_87_sub_irss, jobs 7453142/7453143).
    if len(lines) == 1 and re.match(
        r'^[A-Za-z0-9_][A-Za-z0-9._-]*\.(?:xyz|XYZ|smi|smi|mol2?|sdf|pdb|txt|cif|inp)[\s]*$',
        lines[0],
    ):
        return False

    # SMILES should be on first line (possibly with a comment on second line)
    if len(lines) > 3:
        return False

    first_line = lines[0].strip()

    # Empty or comment line
    if not first_line or first_line.startswith('#') or first_line.startswith('*'):
        return False

    # Check for coordinate patterns (element symbol followed by numbers)
    # XYZ format: "C  1.234  5.678  9.012"
    # The decimal point is not part of what makes a line coordinates: "C 0 0 0"
    # is a perfectly ordinary way to write an atom at the origin, and it was
    # read as a SMILES because it has digits in it and starts with a carbon.
    coord_pattern = re.compile(
        r'^[A-Z][a-z]?\s+[-+]?\d+(?:\.\d+)?\s+[-+]?\d+(?:\.\d+)?'
        r'\s+[-+]?\d+(?:\.\d+)?\s*$')
    if coord_pattern.match(first_line):
        return False

    # Check for XYZ header (first line is a number)
    try:
        int(first_line.split()[0])
        return False  # Looks like XYZ format (atom count)
    except (ValueError, IndexError):
        pass

    # Check for typical SMILES characters. The dot belongs in here: it is how
    # SMILES separates one molecule from another, which is what a drawing of
    # two things hands back -- and "CCO.CCO" has no bracket, no aromatic
    # letter and no ring number, so without it nothing here recognised it and
    # a two-fragment drawing could not be converted at all. Coordinates cannot
    # be caught by it: a line of them is turned away above, by the pattern for
    # an element followed by three numbers and by the atom-count check.
    smiles_chars = set('()[]=#@+-/\\>%.')
    has_smiles_chars = any(char in first_line for char in smiles_chars)

    # Check for aromatic notation (lowercase c, n, o, etc.) which is SMILES-specific
    aromatic_pattern = re.compile(r'[cnops]')
    has_aromatic = bool(aromatic_pattern.search(first_line))

    # Check for organic element symbols (both upper and lowercase for aromatic)
    # (not followed by coordinates)
    organic_pattern = re.compile(r'[CcNnOoPpSsFfClBrI]')
    has_organic = bool(organic_pattern.search(first_line))
    has_metal = contains_metal(first_line)

    # Check for numbered ring closures (like 1 in c1ccccc1)
    has_ring_numbers = bool(re.search(r'\d', first_line))

    # Simple SMILES without aromatic/lone symbols (e.g., "CC") – single token, no spaces/numbers
    simple_token = (
        len(lines) == 1
        and ' ' not in first_line
        and not any(ch.isdigit() for ch in first_line)
        and first_line.isalnum()
        and all(ch.isalpha() for ch in first_line)
    )

    # SMILES if: (special chars OR aromatic OR ring numbers OR simple token) AND has organic elements
    return (has_smiles_chars or has_aromatic or has_ring_numbers or simple_token) and (has_organic or has_metal)


# Perf (byte-identical): per-(smiles, hapto_approx) cache for the embedding-prep
# mol.  For a single metal complex this builder is invoked up to 4x with only 2
# distinct (smiles, hapto_approx) keys (η-label probe, σ-isomer enumerator,
# hapto-diversity gate, ...), and each invocation re-runs the heavy
# stk.BuildingBlock / RDKit sanitize path from scratch.  The function is a pure
# deterministic function of its two args (verified: repeated calls produce
# byte-identical Chem.Mol.ToBinary()), so we memoize the *prepared* mol and hand
# every caller a fresh Chem.Mol() copy.  The copy is bit-faithful
# (ToBinary round-trip identical) and decoupled, so callers that mutate it
# (RemoveAllConformers / AddConformer / dative edits) cannot poison the cache or
# each other.  Bounded LRU keeps memory flat across a 50k pool.
_CACHE_MISS = object()  # distinct sentinel: a None entry is a valid cached miss


_PREP_MOL_CACHE: "OrderedDict[Tuple[str, bool], object]" = OrderedDict()


_PREP_MOL_CACHE_MAX = 16


def _prepare_mol_for_embedding(smiles: str, hapto_approx: bool = False):
    """Cached wrapper over :func:`_prepare_mol_for_embedding_uncached`.

    Returns a fresh, decoupled :class:`rdkit.Chem.Mol` copy of the cached
    prepared mol (or ``None``).  Byte-identical to calling the uncached builder
    directly — the cached object is never handed out, only copies.
    """
    if not RDKIT_AVAILABLE:
        return None
    key = (smiles, bool(hapto_approx))
    cached = _PREP_MOL_CACHE.get(key, _CACHE_MISS)
    if cached is _CACHE_MISS:
        cached = _prepare_mol_for_embedding_uncached(smiles, hapto_approx)
        _PREP_MOL_CACHE[key] = cached
        if len(_PREP_MOL_CACHE) > _PREP_MOL_CACHE_MAX:
            _PREP_MOL_CACHE.popitem(last=False)
    else:
        # mark as recently used (LRU)
        _PREP_MOL_CACHE.move_to_end(key)
    if cached is None:
        return None
    return Chem.Mol(cached)


def _prepare_mol_for_embedding_uncached(smiles: str, hapto_approx: bool = False):
    """Parse SMILES and prepare an RDKit Mol for conformer embedding.

    Tries the same strategies as ``smiles_to_xyz`` /
    ``_try_multiple_strategies`` but stops before embedding so the caller
    can generate multiple conformers.  The returned Mol has **no**
    conformers.

    Returns ``None`` if the SMILES cannot be parsed by any strategy.
    """
    if not RDKIT_AVAILABLE:
        return None

    has_metal = contains_metal(smiles)
    is_metal_n = _is_metal_nitrogen_complex(smiles)

    # Build SMILES variants (same as _try_multiple_strategies)
    normalized_smiles = _normalize_metal_smiles(smiles)
    denormalized_smiles = _denormalize_metal_smiles(smiles)
    smiles_variants = [smiles]
    if normalized_smiles and normalized_smiles not in smiles_variants:
        smiles_variants.append(normalized_smiles)
    if denormalized_smiles and denormalized_smiles not in smiles_variants:
        smiles_variants.append(denormalized_smiles)

    mol = None

    # Strategy 1: stk (best for metal complexes)
    if has_metal and STK_AVAILABLE:
        for smi in smiles_variants:
            try:
                bb = stk.BuildingBlock(smi)
                mol = bb.to_rdkit_mol()
                if mol is not None:
                    break
            except Exception:
                continue

    # Strategy 2: RDKit partial-sanitize with variants
    if mol is None:
        for smi in smiles_variants:
            mol, _ = mol_from_smiles_rdkit(smi, allow_metal=has_metal)
            if mol is not None:
                break

    # Strategy 3: Unsanitized parsing (handles valence errors)
    if mol is None:
        for smi in smiles_variants:
            try:
                try:
                    p = Chem.SmilesParserParams()
                    p.sanitize = False
                    p.removeHs = False
                    p.strictParsing = False
                    mol = Chem.MolFromSmiles(smi, p)
                except Exception:
                    mol = Chem.MolFromSmiles(smi, sanitize=False)
                if mol is not None:
                    for atom in mol.GetAtoms():
                        atom.SetNoImplicit(True)
                    try:
                        mol.UpdatePropertyCache(strict=False)
                    except Exception:
                        pass
                    break
            except Exception:
                continue

    if mol is None:
        return None


    # Optional eta/hapto approximation: collapse contiguous metal-bound carbon
    # donor blocks to a single representative anchor bond per block.
    if has_metal and hapto_approx:
        try:
            hapto_groups = _find_hapto_groups(mol)
            if hapto_groups:
                mol, n_removed = _apply_hapto_approximation(mol, hapto_groups)
                logger.info(
                    "Applied experimental hapto approximation in embedding prep: "
                    "%d group(s), %d bond(s) converted to dative.",
                    len(hapto_groups), n_removed,
                )
        except Exception as e:
            logger.debug("Hapto approximation in embedding prep failed: %s", e)

    # Hydrogen handling depends on the complex type.
    # Metal-nitrogen complexes: keep N/C-metal bonds as SINGLE (required
    # for successful embedding) but convert S/O/P-metal bonds to dative
    # for correct distance bounds.
    # Other metal complexes: full dative bond conversion.
    if has_metal and not is_metal_n:
        try:
            mol = Chem.RemoveHs(mol)
            mol = _convert_metal_bonds_to_dative(mol)
            mol.UpdatePropertyCache(strict=False)
            if _is_simple_organometallic(smiles):
                mol = Chem.AddHs(mol, addCoords=False)
                mol = _strip_h_on_metal_halogen(mol)
                mol = _fix_organometallic_carbon_h(mol)
            else:
                mol = Chem.AddHs(mol, addCoords=False)
            if hapto_approx or _hapto_h_always():
                mol = _fix_hapto_donor_h(mol)
        except Exception:
            try:
                mol = Chem.AddHs(mol, addCoords=False)
                if hapto_approx or _hapto_h_always():
                    mol = _fix_hapto_donor_h(mol)
            except Exception:
                pass
    elif has_metal:
        # Metal-nitrogen: selective dative conversion for C/S/O/P only.
        # N bonds stay SINGLE so that neutral N-donor H counts are correct.
        # C bonds become dative — C donors (ppy, cyclometallated ligands) need
        # dative treatment for ETKDG to generate correct Ir-C/Ru-C distances;
        # there is no H-count issue for C since it is already at full valence.
        # Mark converted atoms NoImplicit to prevent spurious H
        # (dative bond doesn't count toward ligand valence).
        try:
            mol = _convert_metal_bonds_to_dative(mol, only_elements={'C', 'S', 'O', 'P'})
            # Mark ALL neutral atoms bonded to a metal as NoImplicit so
            # AddHs() won't add spurious H on coordinating atoms.
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in _METAL_SET:
                    continue
                has_metal_bond = any(
                    n.GetSymbol() in _METAL_SET for n in atom.GetNeighbors()
                )
                if has_metal_bond and atom.GetFormalCharge() >= 0:
                    atom.SetNoImplicit(True)
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        try:
            mol = Chem.AddHs(mol, addCoords=False)
            if hapto_approx or _hapto_h_always():
                mol = _fix_hapto_donor_h(mol)
        except Exception:
            pass
    else:
        try:
            mol = Chem.AddHs(mol, addCoords=False)
        except Exception:
            pass

    # Remove any conformers that may have been inherited from stk
    mol.RemoveAllConformers()

    return mol


def _dearomatized_embedding_copy(mol):
    """Return a temporary de-aromatized copy suitable for ETKDG fallback.

    Some charged aromatic metal complexes fail RDKit's internal kekulization
    step during conformer embedding. For those cases we generate conformers on
    a temporary copy where aromatic flags/bonds are cleared, then transfer only
    coordinates back to the original molecule.
    """
    if not RDKIT_AVAILABLE:
        return None
    try:
        tmp = Chem.Mol(mol)
        tmp.RemoveAllConformers()
        rw = Chem.RWMol(tmp)
        for atom in rw.GetAtoms():
            if atom.GetIsAromatic():
                atom.SetIsAromatic(False)
        for bond in rw.GetBonds():
            if bond.GetIsAromatic() or bond.GetBondType() == Chem.BondType.AROMATIC:
                bond.SetIsAromatic(False)
                bond.SetBondType(Chem.BondType.SINGLE)
        tmp = rw.GetMol()
        try:
            tmp.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        return tmp
    except Exception:
        return None


def _rescale_metal_donor_distances(mol, conf_id: int) -> None:
    """Scale metal positions so M-D distances match lookup-table ideals.

    ETKDG treats metals like organic atoms and places M-D bonds at
    ~1.4-1.7 Å.  Instead of moving individual donors (which breaks
    chelate/macrocyclic rings), this function scales the METAL position
    relative to the donor centroid.

    For each metal:
    1. Compute centroid of all donor positions
    2. Compute average current M-D distance and average ideal M-D
    3. Scale metal position along M→centroid vector so averages match

    This preserves the ligand geometry while correcting the coordination
    sphere radius.  Modifies the conformer in-place.
    """
    if not RDKIT_AVAILABLE:
        return
    try:
        from rdkit.Geometry import Point3D as _P3D
        conf = mol.GetConformer(conf_id)
    except Exception:
        return

    # Detect hapto metals (≥3 contiguous C donors) — skip rescaling
    # for these because moving them distorts the η-ring which breaks
    # bridging bonds to other metals.
    hapto_metal_indices: set = set()
    try:
        hapto_groups = _find_hapto_groups(mol)
        for _hm, _hc in hapto_groups:
            hapto_metal_indices.add(_hm)
    except Exception:
        pass

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        m_idx = atom.GetIdx()
        if m_idx in hapto_metal_indices:
            continue  # hapto metals keep ETKDG position
        m_sym = atom.GetSymbol()
        m_pos = conf.GetAtomPosition(m_idx)
        mx, my, mz = m_pos.x, m_pos.y, m_pos.z

        # Collect donor positions and ideal distances.
        donors: List[Tuple[float, float, float, float]] = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() <= 1 or nbr.GetSymbol() in _METAL_SET:
                continue
            dp = conf.GetAtomPosition(nbr.GetIdx())
            ideal = float(_get_ml_bond_length(m_sym, nbr.GetSymbol()))
            donors.append((dp.x, dp.y, dp.z, ideal))

        if not donors:
            continue

        # Centroid of donors.
        n = len(donors)
        cx = sum(d[0] for d in donors) / n
        cy = sum(d[1] for d in donors) / n
        cz = sum(d[2] for d in donors) / n

        # Average current distance and average ideal.
        avg_cur = sum(
            math.sqrt((d[0] - mx) ** 2 + (d[1] - my) ** 2 + (d[2] - mz) ** 2)
            for d in donors
        ) / n
        avg_ideal = sum(d[3] for d in donors) / n

        if avg_cur < 0.3 or abs(avg_ideal / avg_cur - 1.0) < 0.05:
            continue

        # Move metal: new_metal = centroid + (metal - centroid) * (ideal/current)
        scale = avg_ideal / avg_cur
        new_mx = cx + (mx - cx) * scale
        new_my = cy + (my - cy) * scale
        new_mz = cz + (mz - cz) * scale
        try:
            conf.SetAtomPosition(m_idx, _P3D(new_mx, new_my, new_mz))
        except Exception:
            pass


def _embed_multiple_confs_robust(
    mol,
    num_confs: int,
    seed: int,
) -> List[int]:
    """Embed conformers with dearomatized fallback for kekulization failures."""
    if not RDKIT_AVAILABLE or num_confs <= 0:
        return []

    def _params_for_seed(_seed: int):
        p = AllChem.ETKDGv3()
        p.useRandomCoords = True
        p.randomSeed = int(_seed)
        p.enforceChirality = False
        try:
            p.clearConfs = False
        except Exception:
            pass
        return p

    # Primary embedding on the original molecule.
    primary_ids: List[int] = []
    try:
        primary_ids = _embed_multiple_confs_with_timeout(
            mol, num_confs, _params_for_seed(seed)
        )
    except Exception as emb_exc:
        logger.debug("Primary ETKDG embedding failed (seed=%s): %s", seed, emb_exc)
    if primary_ids:
        # Rescale M-D distances from organic (~1.5 Å) to ideal values
        # BEFORE any UFF runs. This is the universal fix for metals
        # without UFF parameters.
        for cid in primary_ids:
            try:
                _rescale_metal_donor_distances(mol, cid)
            except Exception:
                pass
        return primary_ids

    # Fallback 1: embed on a temporary de-aromatized copy and transfer coords.
    tmp = _dearomatized_embedding_copy(mol)
    tmp_ids: List[int] = []
    if tmp is not None:
        try:
            tmp_ids = _embed_multiple_confs_with_timeout(
                tmp, num_confs, _params_for_seed(seed)
            )
        except Exception as emb_exc:
            logger.debug(
                "Dearomatized ETKDG embedding failed (seed=%s): %s",
                seed, emb_exc,
            )

    # Fallback 2 (on stall): permissive distance-geometry without ETKDG
    # torsion knowledge. Unlocks highly connected cage / metallacycle SMILES
    # where the ETKDG bounds matrix does not converge. Runs on the
    # dearomatized mol when available, else on the original.
    if not tmp_ids:
        base = tmp if tmp is not None else mol
        try:
            permissive = _make_random_embed_params(int(seed))
            try:
                permissive.clearConfs = False
            except Exception:
                pass
            tmp_ids = _embed_multiple_confs_with_timeout(
                base, num_confs, permissive
            )
        except Exception as emb_exc:
            logger.debug(
                "Permissive distance-geometry embedding failed (seed=%s): %s",
                seed, emb_exc,
            )
        if tmp_ids and base is mol:
            # Already on the target mol — rescale and return.
            for cid in tmp_ids:
                try:
                    _rescale_metal_donor_distances(mol, cid)
                except Exception:
                    pass
            logger.debug(
                "Embedded %d conformer(s) via permissive distance-geometry fallback (seed=%s).",
                len(tmp_ids), seed,
            )
            return list(tmp_ids)

    if not tmp_ids:
        return []

    transferred: List[int] = []
    for tid in tmp_ids:
        try:
            conf = Chem.Conformer(tmp.GetConformer(tid))
            cid = mol.AddConformer(conf, assignId=True)
            transferred.append(cid)
        except Exception:
            continue
    if transferred:
        for cid in transferred:
            try:
                _rescale_metal_donor_distances(mol, cid)
            except Exception:
                pass
        logger.debug(
            "Embedded %d conformers via dearomatized/permissive fallback (seed=%s).",
            len(transferred), seed,
        )
    return transferred
