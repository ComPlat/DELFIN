"""Metal SMILES normalisation (dative bonds, donor hydrogens, organometallic carbons) of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _delfin_env_int,
)
from delfin.manta.hapto_detect import (
    _find_hapto_groups,
    contains_metal,
)
from delfin.manta.ml_tables import (
    Chem,
    RDKIT_AVAILABLE,
    _HALOGENS,
    _METALS,
    _METAL_SET,
)

logger = get_logger("delfin.smiles_converter")


def _normalize_metal_smiles(smiles: str) -> Optional[str]:
    """Normalize neutral metal SMILES to charged form as a fallback.

    Converts common neutral metal notations to their typical charged forms:
    - Ni, Co, Fe, Cu, Zn, Mn → +2
    - Cr, Rh, Ir → +3
    - Ru → +2 (most common in coordination chemistry)
    - Neutral [N] → [N-] for coordination
    """
    # Common oxidation states for metals in coordination complexes
    metal_charges = {
        # 3d metals
        'Ni': '+2', 'Co': '+2', 'Fe': '+2', 'Cu': '+2', 'Zn': '+2',
        'Mn': '+2', 'Cr': '+3', 'V': '+3', 'Ti': '+4', 'Sc': '+3',
        # 4d metals
        'Ru': '+2', 'Rh': '+3', 'Pd': '+2', 'Ag': '+1',
        'Mo': '+4', 'Zr': '+4', 'Nb': '+5', 'Y': '+3', 'Cd': '+2',
        # 5d metals
        'Ir': '+3', 'Pt': '+2', 'Au': '+3', 'Hg': '+2',
        'Os': '+2', 'Re': '+5', 'W': '+6', 'Ta': '+5', 'Hf': '+4',
        # Lanthanides (all typically +3)
        'La': '+3', 'Ce': '+3', 'Pr': '+3', 'Nd': '+3', 'Sm': '+3',
        'Eu': '+3', 'Gd': '+3', 'Tb': '+3', 'Dy': '+3', 'Ho': '+3',
        'Er': '+3', 'Tm': '+3', 'Yb': '+3', 'Lu': '+3',
        # Actinides
        'Th': '+4', 'U': '+4', 'Np': '+4', 'Pu': '+4',
        # Main group
        'Al': '+3', 'Ga': '+3', 'In': '+3',
    }

    normalized = smiles
    found_neutral_metal = False

    for metal, charge in metal_charges.items():
        neutral_pattern = f'[{metal}]'
        if neutral_pattern in normalized:
            # Check if already has a charged version
            if f'[{metal}+' not in smiles and f'[{metal}-' not in smiles:
                normalized = normalized.replace(neutral_pattern, f'[{metal}{charge}]')
                found_neutral_metal = True

    # Also convert neutral [N] to [N-] for coordination
    if found_neutral_metal and '[N]' in normalized and '[N-]' not in smiles:
        normalized = normalized.replace('[N]', '[N-]')

    # Normalize CO ligand variants: add charges to balance valence WITHOUT
    # changing atom order. User notation O#[C][M] means O≡C-M (C bonded
    # to metal, O terminal). We must preserve this order — just add [O+]
    # and [C-] to make valences legal for RDKit.  The carbon may already
    # carry any formal-charge label ([C], [C+], [C-]); those are all
    # equivalent SMILES variants of the same coordinated CO motif and
    # must canonicalise to the same zwitterion so downstream processing
    # is invariant to the user's charge bookkeeping.
    if contains_metal(normalized):
        # CO with O terminal, C bonded to metal: O#[C{+,-,}][M] -> [O+]#[C-][M]
        normalized = re.sub(
            r'(?<![\[\w])O#\[C[+-]?\](?=[\[\d])', '[O+]#[C-]', normalized
        )
        # CO with C bonded to metal, O terminal: [M][C{+,-,}]#O -> [C-]#[O+]
        normalized = re.sub(r'\[C[+-]?\]#O(?![a-zA-Z])', '[C-]#[O+]', normalized)
        # Nitrosyl terminal: [N{+,-,}]=O adjacent to metal -> [N+]=O
        normalized = re.sub(r'\[N[+-]?\]=O(?![a-zA-Z])', '[N+]=O', normalized)

    # Metal formal-charge canonicalisation.  Writing the metal as [Fe-2],
    # [Fe-3], [Fe-5] etc. is user-side bookkeeping of the neutral overall
    # charge and should not alter the coordination topology.  Replace any
    # non-canonical metal charge with the charge from ``metal_charges``
    # when that table has an entry for the element; this keeps enumerator
    # inputs invariant across SMILES variants of the same complex.
    #
    # Welle-3 T1.1 (2026-05-15): env-flag to preserve user-explicit POSITIVE
    # charges ([Pt+4], [Fe+3], [Cu+1], [Ni+3], etc.).  Negative metal charges
    # are always treated as CCDC-style bookkeeping (true negative oxidation
    # states are extremely rare in TM coord chem and would not appear in the
    # smiles_master pool).
    #
    # Welle-5g Step-0b (2026-05-17) — DEFAULT-FLIP 0 -> 1 per 5f-N +5g-V
    # default-flip bisect: PRESERVE_METAL_CHARGE is the best-single-flag
    # measured (+2204 NET on the 5-archive bisect, sigma +1331 / hapto +863).
    # Verdict re-confirmed by 5g-V DEFAULT-FLIP verification (math reproduced
    # byte-for-byte from JSON; no per-class negatives at single-flag level).
    # Disable via DELFIN_PRESERVE_METAL_CHARGE=0 to restore legacy stomping.
    #
    # Welle-5j Step-B (2026-05-17) — charge-magnitude-conditional protection.
    # Per Welle-5i Agent E apples-to-apples recal-battery, the Step-0b default
    # flip introduced real same-detector regressions on F3_angle (+7.01pp),
    # F3_bond (+2.50pp), vdw clashes (+2.37pp) and lig_pct_realistic (-3.95pp).
    # Root-cause hypothesis: for SMILES with extreme metal formal charges,
    # pinning them through to OB-UFF causes the force-field to pick different
    # local minima for some ligands (OB UFF unparametrized-TM-cation bug —
    # see memory/feedback_ob_uff_unparametrized_tm.md), tilting internal
    # organic angles past detector cutoffs even though topology is preserved.
    #
    # DELFIN_PRESERVE_CHARGE_MAGNITUDE_THRESHOLD (int, default 0 = no
    # threshold = current behaviour preserved byte-for-byte) defines a
    # magnitude cutoff N: when N >= 1, any metal whose preserved positive
    # formal charge satisfies abs(q) > N falls back to canonicalisation
    # (replaced with the table-standard charge from ``metal_charges``).  The
    # check is universal — magnitude-based, no per-element allowlist — and
    # only affects the OB-UFF / RDKit-parse normalisation step.  Downstream
    # chemistry, user-facing labels and the original input SMILES (which
    # callers preserve when ``_normalize_metal_smiles`` returns None /
    # unchanged) are untouched.
    _preserve_pos = _delfin_env_int("DELFIN_PRESERVE_METAL_CHARGE", 1) == 1
    _mag_threshold = _delfin_env_int(
        "DELFIN_PRESERVE_CHARGE_MAGNITUDE_THRESHOLD", 0
    )
    for metal, charge in metal_charges.items():
        if _preserve_pos:
            # Only canonicalise negative bookkeeping charges ([Fe-2], [Fe-5]);
            # positive charges ([Fe+3], [Pt+4]) are user chemistry and kept
            # — unless DELFIN_PRESERVE_CHARGE_MAGNITUDE_THRESHOLD is set, in
            # which case positives whose magnitude exceeds the threshold are
            # also canonicalised (see Welle-5j Step-B comment block above).
            pattern = rf'\[{metal}-\d*\]'
        else:
            pattern = rf'\[{metal}[+-]\d*\]'
        if re.search(pattern, normalized):
            normalized = re.sub(pattern, f'[{metal}{charge}]', normalized)
        # Magnitude-threshold canonicalisation of preserved positive charges.
        # Active only when threshold >= 1 AND positive-preservation is on
        # (threshold without preservation is a no-op because the all-charges
        # pattern above already canonicalised every positive metal charge).
        if _preserve_pos and _mag_threshold >= 1:
            # Match [Metal+N] with N >= 1 and canonicalise when N > threshold.
            pos_pattern = rf'\[{metal}\+(\d+)\]'
            def _replace_if_extreme(match: 're.Match[str]') -> str:
                try:
                    q = int(match.group(1))
                except (TypeError, ValueError):
                    return match.group(0)
                if q > _mag_threshold:
                    return f'[{metal}{charge}]'
                return match.group(0)
            normalized = re.sub(pos_pattern, _replace_if_extreme, normalized)

    return normalized if normalized != smiles else None


def _denormalize_metal_smiles(smiles: str) -> Optional[str]:
    """Convert charged metal SMILES back to neutral form as a fallback.

    Handles various charge notations including stereochemistry markers for all metals.
    """
    # Check if we have any charged notation
    if not re.search(r'\[[A-Z][a-z]?(?:@+|@@)?[+-]\d*\]', smiles) and '[N-]' not in smiles:
        return None

    denorm = smiles

    # Convert charged metals to neutral (handle stereochemistry)
    # Pattern matches: [Metal@+2], [Metal@@+3], [Metal+2], [Metal+], [Metal-], etc.
    for metal in _METALS:
        # Pattern: [Metal with optional stereochemistry and charge]
        pattern = rf'\[{metal}(?:@+|@@)?[+-]\d*\]'
        denorm = re.sub(pattern, f'[{metal}]', denorm)

    # Convert [N-] back to [N]
    denorm = denorm.replace('[N-]', '[N]')

    return denorm if denorm != smiles else None


def _strip_h_on_coordinated_p(mol):
    """Remove hydrogens on P/As atoms that are coordinated to a metal.

    When a P-Metal bond is not converted to dative, RDKit may assign
    implicit H to P (P default valence = 3 or 5).  Tertiary phosphine
    ligands (PR₃) coordinated to metals should never carry P-H bonds.
    This strips H from any P or As that is bonded to a metal AND already
    has ≥ 3 non-H, non-metal neighbours (i.e. a full set of R groups).
    """
    if not RDKIT_AVAILABLE:
        return mol
    rwmol = Chem.RWMol(mol)
    to_remove = []
    for atom in rwmol.GetAtoms():
        if atom.GetSymbol() not in ('P', 'As'):
            continue
        # Check if bonded to a metal
        has_metal = any(n.GetSymbol() in _METAL_SET for n in atom.GetNeighbors())
        if not has_metal:
            continue
        # Count non-H, non-metal neighbours (the "R" groups)
        r_count = sum(
            1 for n in atom.GetNeighbors()
            if n.GetSymbol() != 'H' and n.GetSymbol() not in _METAL_SET
        )
        if r_count < 3:
            continue
        # Remove all H bonded to this P/As
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() == 'H':
                to_remove.append(nbr.GetIdx())
    for idx in sorted(set(to_remove), reverse=True):
        rwmol.RemoveAtom(idx)
    return rwmol.GetMol()


def _strip_h_on_metal_halogen(mol):
    """Remove hydrogens attached to metals/halogens (keep H on carbon)."""
    if not RDKIT_AVAILABLE:
        return mol
    rwmol = Chem.RWMol(mol)
    to_remove = []
    for atom in rwmol.GetAtoms():
        if atom.GetSymbol() == 'H':
            nbrs = atom.GetNeighbors()
            if not nbrs:
                continue
            nbr = nbrs[0]
            sym = nbr.GetSymbol()
            if sym in _METAL_SET or sym in _HALOGENS:
                to_remove.append(atom.GetIdx())
    for idx in sorted(to_remove, reverse=True):
        rwmol.RemoveAtom(idx)
    return rwmol.GetMol()


def _fix_organometallic_carbon_h(mol):
    """Ensure carbon bonded to a metal has a reasonable number of H atoms."""
    if not RDKIT_AVAILABLE:
        return mol
    rwmol = Chem.RWMol(mol)
    to_remove = []
    for atom in rwmol.GetAtoms():
        if atom.GetSymbol() != 'C':
            continue
        # Check if carbon is bonded to a metal
        nbrs = atom.GetNeighbors()
        if not any(n.GetSymbol() in _METAL_SET for n in nbrs):
            continue
        # Count heavy (non-H) neighbors
        heavy_nbrs = [n for n in nbrs if n.GetSymbol() != 'H']
        desired_h = max(0, 4 - len(heavy_nbrs))
        h_nbrs = [n for n in nbrs if n.GetSymbol() == 'H']
        if len(h_nbrs) > desired_h:
            # Remove extra H atoms (deterministically by index)
            remove_count = len(h_nbrs) - desired_h
            for h in sorted(h_nbrs, key=lambda a: a.GetIdx())[:remove_count]:
                to_remove.append(h.GetIdx())
    for idx in sorted(set(to_remove), reverse=True):
        rwmol.RemoveAtom(idx)
    return rwmol.GetMol()


def _hapto_h_always() -> bool:
    """Should the hapto donor-H repair run outside the experimental hapto mode?

    WHY THIS EXISTS (user, 2026-07-31, on KEHWEA).  The user opened a built frame and saw
    the hydrogens missing on the eta2 alkene coordinating the Pd.  They are missing in the
    INPUT: the SMILES writes the pi interaction as a sigma bond,

        ... [C+]1=[C+](C6H5) -> [Pd-4] ...

    and a bracketed carbon with no H spec carries ZERO hydrogens, so every eta2 alkene CH,
    eta5 Cp CH and eta6 arene CH loses its proton before construction ever begins.  The
    builder is not at fault; it builds exactly what it is handed, and no seating can rescue
    a molecule that is short an atom.

    Measured over the 1000-system pool with the hapto discriminator (in a hapto block
    SEVERAL CONTIGUOUS carbons bond the SAME metal, in a sigma bond exactly one -- which is
    what separates a real eta2 alkene from an NHC carbene carbon or a sigma-aryl ipso
    carbon, both correctly H-free): 52 systems, 322 hydrogens missing.

    The repair itself already existed and was correct -- it was simply never reached,
    because it only ran under the experimental ``hapto_approx`` mode.  Same shape as the
    other holes found this week: not a missing mechanism, a mechanism out of scope.

    Default OFF -> byte-identical.
    """
    return _delfin_env_int("DELFIN_FFFREE_HAPTO_H", 0) == 1


def _fix_hapto_donor_h(mol):
    """Adjust H on hapto donor carbons so each C has at most 4 bonding partners.

    Simple rule: desired_h = max(0, 4 - number_of_non_H_neighbors).
    Metal neighbors ARE counted (they occupy a coordination site).

    Scope is already exactly right and must stay that way: it acts only inside blocks of
    >= 2 contiguous metal-bound carbons, so a carbene or a sigma-aryl is never touched.
    For the hapto cases the arithmetic lands on the correct answer -- an eta-bound CH has
    two heavy ring/chain neighbours plus the metal, so 4 - 3 = 1 hydrogen.
    """
    if not RDKIT_AVAILABLE or mol is None:
        return mol
    try:
        hapto_groups = _find_hapto_groups(mol)
    except Exception:
        return mol
    if not hapto_groups:
        return mol

    rwmol = Chem.RWMol(mol)
    to_remove = []
    to_add: List[int] = []
    for _metal_idx, grp in hapto_groups:
        if len(grp) < 2:
            continue
        for c_idx in grp:
            atom = rwmol.GetAtomWithIdx(c_idx)
            if atom.GetSymbol() != 'C':
                continue
            nbrs = list(atom.GetNeighbors())
            h_nbrs = [n for n in nbrs if n.GetSymbol() == 'H']
            non_h_count = sum(1 for n in nbrs if n.GetSymbol() != 'H')
            desired_h = max(0, 4 - non_h_count)

            if len(h_nbrs) > desired_h:
                remove_count = len(h_nbrs) - desired_h
                for h in sorted(h_nbrs, key=lambda a: a.GetIdx())[:remove_count]:
                    to_remove.append(h.GetIdx())
            elif len(h_nbrs) < desired_h:
                to_add.extend([c_idx] * (desired_h - len(h_nbrs)))

    for idx in sorted(set(to_remove), reverse=True):
        rwmol.RemoveAtom(idx)
    for c_idx in to_add:
        if c_idx < 0 or c_idx >= rwmol.GetNumAtoms():
            continue
        h_idx = rwmol.AddAtom(Chem.Atom('H'))
        rwmol.AddBond(c_idx, h_idx, Chem.BondType.SINGLE)
    out = rwmol.GetMol()
    try:
        out.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return out


def _convert_metal_bonds_to_dative(mol, only_elements=None):
    """Convert single bonds from NEUTRAL atoms to metals to dative bonds.

    RDKit counts metal coordination bonds towards the valence of ligand atoms,
    but dative/coordinative bonds should not count towards ligand valence.
    By converting SINGLE bonds to DATIVE bonds, RDKit correctly calculates
    implicit hydrogens on the ligand atoms.

    IMPORTANT: Only NEUTRAL atoms get their bonds converted to dative.
    - [N] bound to metal → dative bond → H atoms calculated normally
    - [N+] bound to metal → remains covalent → no extra H atoms

    Based on RDKit Cookbook: https://www.rdkit.org/docs/Cookbook.html

    Args:
        mol: RDKit Mol object (will be modified)
        only_elements: If set, only convert bonds from these elements
            (e.g. ``{'S', 'O', 'P'}``).  None means convert all.

    Returns:
        RDKit Mol object with dative bonds to metals
    """
    if not RDKIT_AVAILABLE:
        return mol

    rwmol = Chem.RWMol(mol)
    # Preserve explicit H annotations from the input SMILES for non-donor
    # atoms (e.g., neutral ring N-H). RDKit often cannot reconstruct these
    # from implicit valence on partially-sanitized metal-complex mols.
    orig_explicit_h: Dict[int, int] = {
        a.GetIdx(): int(a.GetNumExplicitHs()) for a in rwmol.GetAtoms()
    }
    orig_no_implicit: Dict[int, bool] = {
        a.GetIdx(): bool(a.GetNoImplicit()) for a in rwmol.GetAtoms()
    }

    # Track dative donors that already exist in the input graph.
    # Needed for hapto approximation where some M-C contacts are pre-marked
    # as dative before this function runs.
    dative_donor_indices = set()
    for bond in rwmol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DATIVE:
            continue
        b = bond.GetBeginAtom()
        e = bond.GetEndAtom()
        sb, se = b.GetSymbol(), e.GetSymbol()
        if sb in _METAL_SET and se not in _METAL_SET:
            dative_donor_indices.add(e.GetIdx())
        elif se in _METAL_SET and sb not in _METAL_SET:
            dative_donor_indices.add(b.GetIdx())

    # Find all single bonds between metals and NEUTRAL non-metals
    bonds_to_convert = []
    for bond in rwmol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        s1, s2 = a1.GetSymbol(), a2.GetSymbol()

        # Check if one is metal and one is not
        is_metal_1 = s1 in _METAL_SET
        is_metal_2 = s2 in _METAL_SET

        if is_metal_1 != is_metal_2:  # XOR - exactly one is metal
            # Determine which atom is the ligand (non-metal)
            ligand_atom = a2 if is_metal_1 else a1

            # Convert if the ligand atom is neutral or POSITIVELY charged.
            # Positive donors ([P+], [N+], [N@@+]) coordinate via lone pairs
            # (dative). Negative donors ([O-], [N-]) have genuinely covalent
            # bonds to the metal and must stay as-is.
            if ligand_atom.GetFormalCharge() >= 0:
                # If only_elements is set, skip elements not in the set
                if only_elements and ligand_atom.GetSymbol() not in only_elements:
                    continue
                bonds_to_convert.append((
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    is_metal_1  # True if atom1 is the metal
                ))

    if not bonds_to_convert and not dative_donor_indices:
        return mol

    # Convert bonds to dative (ligand -> metal direction)
    for idx1, idx2, atom1_is_metal in bonds_to_convert:
        rwmol.RemoveBond(idx1, idx2)
        if atom1_is_metal:
            # atom2 is ligand, atom1 is metal: ligand -> metal
            rwmol.AddBond(idx2, idx1, Chem.BondType.DATIVE)
            dative_donor_indices.add(idx2)
        else:
            # atom1 is ligand, atom2 is metal: ligand -> metal
            rwmol.AddBond(idx1, idx2, Chem.BondType.DATIVE)
            dative_donor_indices.add(idx1)

    result_mol = rwmol.GetMol()

    # Map hapto donor carbons to eta group size so H handling can be tuned:
    # eta>=4 (Cp-like) should not gain implicit H; eta3 should still be able
    # to carry H where valence allows (e.g., allyl/propenyl motifs).
    hapto_group_size_by_atom: Dict[int, int] = {}
    try:
        for _midx, _grp in _find_hapto_groups(result_mol):
            _sz = len(_grp)
            for _aidx in _grp:
                _prev = hapto_group_size_by_atom.get(_aidx, 0)
                if _sz > _prev:
                    hapto_group_size_by_atom[_aidx] = _sz
    except Exception:
        hapto_group_size_by_atom = {}

    # Reset atom properties to allow recalculation of implicit hydrogens.
    # Dative donors: set NoImplicit=True so AddHs() won't add spurious H
    # (trust the H count from the original SMILES).
    # Non-coordinating neutral atoms: recalculate H normally.
    for atom in result_mol.GetAtoms():
        if atom.GetFormalCharge() >= 0 and atom.GetSymbol() not in _METAL_SET:
            if atom.GetIdx() in dative_donor_indices:
                hapto_sz = hapto_group_size_by_atom.get(atom.GetIdx(), 0)
                if atom.GetSymbol() == 'C' and hapto_sz == 3:
                    atom.SetNoImplicit(False)
                    atom.SetNumExplicitHs(0)
                else:
                    atom.SetNoImplicit(True)
            else:
                atom.SetNoImplicit(False)
                # Preserve explicit H on atoms where RDKit's implicit H
                # calculation can't recover them (e.g., tetracoordinate B in
                # pyrazolylborate/scorpionate ligands: B standard valence = 3
                # but actual valence = 4 with the B-H bond).
                orig_h = orig_explicit_h.get(atom.GetIdx(), 0)
                if orig_h > 0:
                    atom.SetNumExplicitHs(orig_h)
                    if orig_no_implicit.get(atom.GetIdx(), False):
                        atom.SetNoImplicit(True)
                elif atom.GetSymbol() == 'B' and atom.GetNumExplicitHs() > 0:
                    atom.SetNoImplicit(True)
                else:
                    atom.SetNumExplicitHs(0)

    logger.info(f"Converted {len(bonds_to_convert)} neutral ligand-metal bonds to dative bonds")

    return result_mol
