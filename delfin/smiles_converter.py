"""SMILES to XYZ conversion using RDKit with metal complex support.

Note: RDKit generates rough starting geometries that are NOT chemically accurate,
especially for metal complexes. These coordinates should ALWAYS be optimized with
GOAT/xTB before running ORCA calculations.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available - SMILES conversion will not work. Install with: pip install rdkit")

# Try to import stk
try:
    import stk
    STK_AVAILABLE = True
except ImportError:
    STK_AVAILABLE = False


_METALS = [
    'Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
    'Ir', 'Pt', 'Au', 'Hg', 'Al', 'Ga', 'In', 'Tl', 'Sn', 'Pb',
    'Bi', 'Po', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
]

_METAL_SET = set(_METALS)
_ORGANOMETALLIC_METALS = {'Li', 'Na', 'K', 'Mg', 'Zn', 'Al'}
_HALOGENS = {'F', 'Cl', 'Br', 'I'}


def _is_simple_organometallic(smiles: str) -> bool:
    """Heuristic: simple organometallics like C[Mg]Br where adding Hs is wrong."""
    # Look for bracketed metals commonly used in organometallic reagents
    for m in _ORGANOMETALLIC_METALS:
        if f'[{m}]' in smiles:
            # If a halogen is present, it's likely a reagent like Grignard
            if any(h in smiles for h in _HALOGENS):
                return True
            # Direct carbon-metal pattern (very rough)
            if f'C[{m}]' in smiles or f'[{m}]C' in smiles:
                return True
    return False


def _prefer_no_sanitize(smiles: str) -> bool:
    """Heuristic: avoid initial sanitize for neutral metal + [N] SMILES."""
    # Check for any neutral metal pattern
    has_neutral_metal = any(f'[{m}]' in smiles for m in _METALS)
    has_neutral_n = '[N]' in smiles

    if has_neutral_metal and has_neutral_n:
        # Only prefer no-sanitize if no charged forms are present
        has_charged = '[N-]' in smiles or any(f'[{m}+' in smiles or f'[{m}-' in smiles for m in _METALS)
        if not has_charged:
            return True
    return False


def _is_metal_nitrogen_complex(smiles: str) -> bool:
    """Check if SMILES represents a metal-nitrogen coordination complex.

    Matches both neutral ([Metal]) and charged ([Metal+2], [Metal-]) notations,
    as well as stereochemical variants ([Metal@@+2], etc.) for all metals.
    """
    # Build pattern for all metals (neutral or charged, with optional stereochemistry)
    metal_pattern = '|'.join(re.escape(m) for m in _METALS)
    has_metal = bool(re.search(rf'\[(?:{metal_pattern})(?:@+|@@)?(?:[+-]\d*)?\]', smiles))

    # Check for coordinating nitrogen (neutral [N] or negative [N-] or ring-closing N)
    has_coord_n = '[N]' in smiles or '[N-]' in smiles or bool(re.search(r'N\d+=', smiles))

    return has_metal and has_coord_n


def _try_multiple_strategies(smiles: str, output_path: Optional[str] = None) -> Tuple[Optional[str], Optional[str]]:
    """Try multiple parsing strategies for metal complexes.

    IMPORTANT: Does NOT convert bonds to dative - both neutral and charged
    SMILES notations should produce the SAME result (no extra H atoms).

    This function tries different approaches in order of preference:
    1. stk (if available) - often handles complex coordination well
    2. RDKit with original SMILES
    3. RDKit with normalized (charged) SMILES
    4. RDKit with denormalized (neutral) SMILES
    5. RDKit unsanitized fallback (no H addition)

    Returns:
        Tuple of (xyz_content, error_message)
    """
    errors = []
    normalized = _normalize_metal_smiles(smiles)
    denormalized = _denormalize_metal_smiles(smiles)

    # Build list of SMILES variants to try (original + alternatives)
    smiles_variants = [smiles]
    if normalized and normalized not in smiles_variants:
        smiles_variants.append(normalized)
    if denormalized and denormalized not in smiles_variants:
        smiles_variants.append(denormalized)

    # Strategy 1: Try stk first (good for metal complexes)
    # stk handles H atoms correctly based on the SMILES notation
    if STK_AVAILABLE:
        for smi in smiles_variants:
            try:
                bb = stk.BuildingBlock(smi)
                mol = bb.to_rdkit_mol()
                if mol is not None:
                    # Embed if needed
                    if mol.GetNumConformers() == 0:
                        params = AllChem.ETKDGv3()
                        params.randomSeed = 42
                        params.useRandomCoords = True
                        result = AllChem.EmbedMolecule(mol, params)
                        if result != 0:
                            continue

                    xyz_content = _mol_to_xyz(mol)
                    if output_path:
                        Path(output_path).write_text(xyz_content, encoding='utf-8')
                        logger.info(f"Converted SMILES to XYZ using stk: {output_path}")
                    return xyz_content, None
            except Exception as e:
                errors.append(f"stk({smi[:30]}...): {e}")

    # Strategy 2: RDKit with partial sanitize (no dative conversion)
    for smi in smiles_variants:
        try:
            mol, note = mol_from_smiles_rdkit(smi, allow_metal=True)
            if mol is not None:
                # Try to add H without modifying bond types
                try:
                    mol = Chem.AddHs(mol, addCoords=True)
                except Exception:
                    pass  # Continue without H if it fails

                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                params.useRandomCoords = True
                result = AllChem.EmbedMolecule(mol, params)
                if result == 0:
                    xyz_content = _mol_to_xyz(mol)
                    if output_path:
                        Path(output_path).write_text(xyz_content, encoding='utf-8')
                        logger.info(f"Converted SMILES to XYZ using RDKit: {output_path}")
                    return xyz_content, None
        except Exception as e:
            errors.append(f"RDKit({smi[:30]}...): {e}")

    # Strategy 3: Unsanitized fallback (preserves SMILES as-is)
    for smi in smiles_variants:
        xyz_content, err = _smiles_to_xyz_unsanitized_fallback(smi)
        if xyz_content:
            if output_path:
                Path(output_path).write_text(xyz_content, encoding='utf-8')
                logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
            return xyz_content, None
        if err:
            errors.append(f"unsanitized({smi[:30]}...): {err}")

    # Strategy 4: Ultra-permissive (bypass all valence checks)
    for smi in smiles_variants:
        xyz_content, err = _smiles_to_xyz_no_valence_check(smi)
        if xyz_content:
            if output_path:
                Path(output_path).write_text(xyz_content, encoding='utf-8')
                logger.info(f"Converted SMILES to XYZ using no-valence-check fallback: {output_path}")
            return xyz_content, None
        if err:
            errors.append(f"no-valence({smi[:30]}...): {err}")

    return None, f"All strategies failed: {'; '.join(errors)}"


def _smiles_to_xyz_no_valence_check(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Ultra-permissive fallback: bypass all valence checks and geometry rules.

    This is for problematic metal complexes where RDKit's valence rules
    don't apply (unusual coordination numbers, complex ring systems, etc.).
    No H atoms are added - only what's in the SMILES.
    Uses minimal geometry knowledge to handle exotic structures.
    """
    if not RDKIT_AVAILABLE:
        return None, "RDKit is not installed"

    try:
        # Parse without sanitization
        try:
            params = Chem.SmilesParserParams()
            params.sanitize = False
            params.removeHs = False
            params.strictParsing = False
            mol = Chem.MolFromSmiles(smiles, params)
        except Exception:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)

        if mol is None:
            return None, "Failed to parse SMILES"

        # Disable all implicit H handling to avoid valence checks
        for atom in mol.GetAtoms():
            atom.SetNoImplicit(True)

        # Use minimal geometry requirements - this handles exotic metal complexes
        # that don't follow standard geometry rules
        embed_params = AllChem.EmbedParameters()
        embed_params.randomSeed = 42
        embed_params.useRandomCoords = True
        embed_params.ignoreSmoothingFailures = True
        embed_params.useExpTorsionAnglePrefs = False  # Don't enforce torsion angles
        embed_params.useBasicKnowledge = False  # Don't use standard geometry knowledge
        embed_params.enforceChirality = False

        result = AllChem.EmbedMolecule(mol, embed_params)

        if result != 0:
            # Try with different random seeds
            for seed in [123, 456, 789, 1000]:
                embed_params.randomSeed = seed
                result = AllChem.EmbedMolecule(mol, embed_params)
                if result == 0:
                    break

        if result != 0:
            return None, "Failed to generate 3D coordinates"

        # Try to add H atoms after embedding (coordinates will be placed)
        try:
            # Reset noImplicit for non-metal atoms to allow H calculation
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    atom.SetNoImplicit(False)
            mol.UpdatePropertyCache(strict=False)
            mol = Chem.AddHs(mol, addCoords=True)
        except Exception:
            # If H addition fails, continue with what we have
            pass

        xyz_content = _mol_to_xyz(mol)
        return xyz_content, None

    except Exception as e:
        return None, f"No-valence-check error: {e}"


def _normalize_metal_smiles(smiles: str) -> Optional[str]:
    """Normalize neutral metal SMILES to charged form as a fallback.

    Converts common neutral metal notations to their typical charged forms:
    - Ni, Co, Fe → +2
    - Ru, Rh, Ir → +3
    - Neutral [N] → [N-] for coordination
    """
    # Common oxidation states for metals in coordination complexes
    metal_charges = {
        'Ni': '+2', 'Co': '+2', 'Fe': '+2', 'Cu': '+2', 'Zn': '+2',
        'Mn': '+2', 'Cr': '+3', 'Ru': '+2', 'Rh': '+3', 'Ir': '+3',
        'Pd': '+2', 'Pt': '+2', 'Au': '+3', 'Ag': '+1',
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


def _convert_metal_bonds_to_dative(mol):
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

    Returns:
        RDKit Mol object with dative bonds to metals
    """
    if not RDKIT_AVAILABLE:
        return mol

    rwmol = Chem.RWMol(mol)

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

            # Only convert if the ligand atom is NEUTRAL (no formal charge)
            # Charged atoms like [N+] have covalent bonds that should stay covalent
            if ligand_atom.GetFormalCharge() == 0:
                bonds_to_convert.append((
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    is_metal_1  # True if atom1 is the metal
                ))

    if not bonds_to_convert:
        return mol

    # Convert bonds to dative (ligand -> metal direction)
    for idx1, idx2, atom1_is_metal in bonds_to_convert:
        rwmol.RemoveBond(idx1, idx2)
        if atom1_is_metal:
            # atom2 is ligand, atom1 is metal: ligand -> metal
            rwmol.AddBond(idx2, idx1, Chem.BondType.DATIVE)
        else:
            # atom1 is ligand, atom2 is metal: ligand -> metal
            rwmol.AddBond(idx1, idx2, Chem.BondType.DATIVE)

    result_mol = rwmol.GetMol()

    # Reset atom properties to allow recalculation of implicit hydrogens
    # Only for NEUTRAL atoms - charged atoms should keep their state
    for atom in result_mol.GetAtoms():
        if atom.GetFormalCharge() == 0 and atom.GetSymbol() not in _METAL_SET:
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)

    logger.info(f"Converted {len(bonds_to_convert)} neutral ligand-metal bonds to dative bonds")

    return result_mol


def contains_metal(smiles: str) -> bool:
    """Return True if SMILES likely contains a metal atom."""
    for metal in _METALS:
        if re.search(rf'\[{metal}[+\-\d\]@]', smiles, re.IGNORECASE):
            return True
        if re.search(rf'\[{metal}\]', smiles, re.IGNORECASE):
            return True
    return False


def mol_from_smiles_rdkit(smiles: str, allow_metal: bool = False):
    """Create RDKit Mol, with relaxed sanitizing for metal complexes."""
    try:
        if not _prefer_no_sanitize(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return mol, None
        if not allow_metal:
            return None, "Failed to parse SMILES string"
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None, "Failed to parse SMILES (no sanitize)"
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                    ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ),
            )
        except Exception:
            pass
        # Update property cache to compute implicit H counts for AddHs
        # Avoid property-cache updates for neutral Ni/Co+[N] SMILES to prevent valence errors
        if not _prefer_no_sanitize(smiles):
            try:
                mol.UpdatePropertyCache(strict=False)
            except Exception:
                pass
        return mol, "partial sanitize"
    except Exception as e:
        return None, f"RDKit error: {e}"

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

    # SMILES should be on first line (possibly with a comment on second line)
    if len(lines) > 3:
        return False

    first_line = lines[0].strip()

    # Empty or comment line
    if not first_line or first_line.startswith('#') or first_line.startswith('*'):
        return False

    # Check for coordinate patterns (element symbol followed by numbers)
    # XYZ format: "C  1.234  5.678  9.012"
    coord_pattern = re.compile(r'^[A-Z][a-z]?\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+')
    if coord_pattern.match(first_line):
        return False

    # Check for XYZ header (first line is a number)
    try:
        int(first_line.split()[0])
        return False  # Looks like XYZ format (atom count)
    except (ValueError, IndexError):
        pass

    # Check for typical SMILES characters
    smiles_chars = set('()[]=#@+-/\\>%')
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


def _prepare_mol_for_embedding(smiles: str):
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

    # Hydrogen handling depends on the complex type.
    # Metal-nitrogen complexes: match _try_multiple_strategies (no dative
    # conversion, simple AddHs).
    # Other metal complexes: dative bond conversion first.
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
        except Exception:
            try:
                mol = Chem.AddHs(mol, addCoords=False)
            except Exception:
                pass
    elif has_metal:
        # Metal-nitrogen: try AddHs but tolerate failure
        try:
            mol = Chem.AddHs(mol, addCoords=False)
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


def _angle_class(pos_metal, pos_a, pos_b) -> str:
    """Classify the angle A-Metal-B as 'cis' (<120 deg) or 'trans'."""
    v1 = (pos_a.x - pos_metal.x, pos_a.y - pos_metal.y, pos_a.z - pos_metal.z)
    v2 = (pos_b.x - pos_metal.x, pos_b.y - pos_metal.y, pos_b.z - pos_metal.z)
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    if mag1 < 1e-8 or mag2 < 1e-8:
        return 'cis'
    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    angle_deg = math.degrees(math.acos(cos_angle))
    return 'cis' if angle_deg < 120 else 'trans'


def _compute_coordination_fingerprint(mol, conf_id: int) -> tuple:
    """Compute a hashable fingerprint describing the coordination geometry.

    Uses **all** pairwise angles between coordinating atoms (not just
    same-element pairs).  This captures bidentate ligand flips where
    donor atoms of different elements swap positions.

    For each metal centre:
    1. Find all coordinating atoms (neighbours of the metal).
    2. For every pair, compute the angle through the metal.
    3. Classify as cis/trans, label with the sorted element pair.
    4. Return sorted tuple of ((elemA, elemB), angle_class) entries.
    """
    conf = mol.GetConformer(conf_id)
    fp_parts: List[tuple] = []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue

        metal_pos = conf.GetAtomPosition(atom.GetIdx())

        # Collect coordinating atoms as (element, idx) sorted by element
        coord_atoms = sorted(
            ((nbr.GetSymbol(), nbr.GetIdx()) for nbr in atom.GetNeighbors()),
            key=lambda x: (x[0], x[1]),
        )

        pair_classes: List[tuple] = []
        for i in range(len(coord_atoms)):
            for j in range(i + 1, len(coord_atoms)):
                sym_a, idx_a = coord_atoms[i]
                sym_b, idx_b = coord_atoms[j]
                cls = _angle_class(
                    metal_pos,
                    conf.GetAtomPosition(idx_a),
                    conf.GetAtomPosition(idx_b),
                )
                pair_key = tuple(sorted((sym_a, sym_b)))
                pair_classes.append((pair_key, cls))

        fp_parts.append(tuple(sorted(pair_classes)))

    return tuple(sorted(fp_parts))


def _classify_isomer_label(fingerprint: tuple, mol) -> str:
    """Translate a coordination fingerprint to a human-readable label.

    Extracts same-element pair data from the full pairwise fingerprint,
    then applies the rules:
    - 3 identical donors, all cis -> **fac**
    - 3 identical donors, 2 cis + 1 trans -> **mer**
    - 2 identical donors, all consistent cis -> **cis**
    - 2 identical donors, all consistent trans -> **trans**
    - Otherwise -> empty string (caller assigns "Isomer N")
    """
    labels_found: List[str] = []
    for metal_fp in fingerprint:
        # Collect same-element pair angles grouped by element
        same_elem: Dict[str, List[str]] = {}
        for pair_key, cls in metal_fp:
            if pair_key[0] == pair_key[1]:
                same_elem.setdefault(pair_key[0], []).append(cls)

        for sym in sorted(same_elem):
            classes = same_elem[sym]
            n_cis = classes.count('cis')
            n_trans = classes.count('trans')
            total = n_cis + n_trans
            # 3 identical donors -> 3 pairs
            if total == 3:
                if n_cis == 3:
                    return 'fac'
                if n_cis == 2 and n_trans == 1:
                    return 'mer'
            # 2 identical donors -> 1 pair
            if total == 1:
                if n_cis == 1:
                    labels_found.append('cis')
                elif n_trans == 1:
                    labels_found.append('trans')
    # Only return cis/trans if all donor groups agree
    if labels_found and len(set(labels_found)) == 1:
        return labels_found[0]
    return ''


def smiles_to_xyz_isomers(
    smiles: str,
    num_confs: int = 50,
    max_isomers: int = 10,
) -> Tuple[List[Tuple[str, str]], Optional[str]]:
    """Generate distinct coordination isomers for a SMILES string.

    For non-metal molecules a single geometry is returned (delegates to
    ``smiles_to_xyz``).  For metal complexes multiple conformers are
    embedded, their coordination fingerprints are computed, and one
    representative per unique fingerprint is returned.

    Returns ``([(xyz_string, label), ...], error)``.
    """
    if not RDKIT_AVAILABLE:
        return [], "RDKit is not installed"

    # Non-metal molecules: single geometry
    if not contains_metal(smiles):
        xyz, err = smiles_to_xyz(smiles)
        if err:
            return [], err
        return [(xyz, '')], None

    # Prepare molecule for embedding
    mol = _prepare_mol_for_embedding(smiles)
    if mol is None:
        # Fall back to single-conformer conversion
        xyz, err = smiles_to_xyz(smiles)
        if err:
            return [], err
        return [(xyz, '')], None

    # Embed multiple conformers with random seeds
    try:
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.randomSeed = -1  # true random for diversity
        params.enforceChirality = False
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    except Exception as e:
        logger.warning("Multi-conformer embedding failed: %s", e)
        xyz, err = smiles_to_xyz(smiles)
        if err:
            return [], err
        return [(xyz, '')], None

    if not conf_ids:
        xyz, err = smiles_to_xyz(smiles)
        if err:
            return [], err
        return [(xyz, '')], None

    # Classify each conformer by label, deduplicate.
    # Known labels (fac/mer/cis/trans): keep one representative per label.
    # Unknown geometry: keep one representative per raw fingerprint.
    seen_labels: Dict[str, int] = {}       # label -> conf_id
    seen_unknown_fps: Dict[tuple, int] = {}  # fingerprint -> conf_id
    for cid in conf_ids:
        try:
            fp = _compute_coordination_fingerprint(mol, cid)
        except Exception:
            continue
        label = _classify_isomer_label(fp, mol)
        if label:
            if label not in seen_labels:
                seen_labels[label] = cid
        else:
            if fp not in seen_unknown_fps:
                seen_unknown_fps[fp] = cid
        if len(seen_labels) + len(seen_unknown_fps) >= max_isomers:
            break

    if not seen_labels and not seen_unknown_fps:
        xyz, err = smiles_to_xyz(smiles)
        if err:
            return [], err
        return [(xyz, '')], None

    # Build results: known labels first; only include unknowns when no
    # known labels were found (unknowns alongside fac/mer/cis/trans are
    # geometry artifacts from RDKit's rough embedding).
    results: List[Tuple[str, str]] = []
    for label, cid in seen_labels.items():
        try:
            xyz = _mol_to_xyz_conformer(mol, cid)
        except Exception:
            continue
        results.append((xyz, label))

    if not results:
        unknown_counter = 0
        for fp, cid in seen_unknown_fps.items():
            unknown_counter += 1
            try:
                xyz = _mol_to_xyz_conformer(mol, cid)
            except Exception:
                continue
            results.append((xyz, f'Isomer {unknown_counter}'))

    if not results:
        xyz, err = smiles_to_xyz(smiles)
        if err:
            return [], err
        return [(xyz, '')], None

    # Single isomer: fall back to smiles_to_xyz() which uses a
    # deterministic seed and produces more realistic geometry.
    if len(results) == 1:
        xyz, err = smiles_to_xyz(smiles)
        if err:
            return results, None  # keep the multi-conf result as fallback
        return [(xyz, results[0][1])], None

    return results, None


def smiles_to_xyz(smiles: str, output_path: Optional[str] = None) -> Tuple[Optional[str], Optional[str]]:
    """Convert SMILES string to XYZ coordinates using RDKit.

    Uses RDKit's ETKDG (Experimental Torsion Knowledge Distance Geometry) method
    for generating reasonable 3D conformers. This is suitable for initial geometries
    that will be further optimized with GOAT/xTB/ORCA.

    Supports coordination bonds (>) in SMILES for metal complexes.

    Args:
        smiles: SMILES string to convert
        output_path: Optional path to write XYZ file

    Returns:
        Tuple of (xyz_content, error_message)
        - xyz_content: XYZ format string if successful, None on error
        - error_message: Error description if failed, None on success
    """
    if not RDKIT_AVAILABLE:
        error = "RDKit is not installed. Install with: pip install rdkit"
        logger.error(error)
        return None, error

    try:
        has_metal = contains_metal(smiles)
        method = None

        # For metal-nitrogen coordination complexes (both neutral and charged notation),
        # use the multi-strategy approach that tries multiple parsing methods
        if _is_metal_nitrogen_complex(smiles):
            logger.info("Detected metal-nitrogen complex, using multi-strategy approach")
            return _try_multiple_strategies(smiles, output_path)

        # Parse SMILES - try stk first for metal complexes, then RDKit
        mol = None
        normalized_smiles = _normalize_metal_smiles(smiles)
        if has_metal and STK_AVAILABLE:
            try:
                bb = stk.BuildingBlock(smiles)
                mol = bb.to_rdkit_mol()
                method = "stk"
            except Exception as e:
                logger.info("stk conversion failed, falling back to RDKit: %s", e)
                mol = None
                method = None

        if mol is None and normalized_smiles:
            mol, rdkit_note = mol_from_smiles_rdkit(normalized_smiles, allow_metal=True)
            if mol is not None:
                method = "RDKit (normalized metal SMILES)"

        if mol is None:
            mol, rdkit_note = mol_from_smiles_rdkit(smiles, allow_metal=has_metal)
            if mol is None:
                # Try normalized charged form for neutral metal SMILES
                if normalized_smiles:
                    mol2, rdkit_note2 = mol_from_smiles_rdkit(normalized_smiles, allow_metal=True)
                    if mol2 is not None:
                        mol = mol2
                        method = "RDKit (normalized metal SMILES)"
                        rdkit_note = None
                    else:
                        rdkit_note = rdkit_note2 or rdkit_note

                if rdkit_note and ("Explicit valence" in rdkit_note or "kekulize" in rdkit_note):
                    legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                    if legacy_err is None and legacy_xyz:
                        if output_path:
                            Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                            logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                        return legacy_xyz, None
                error = f"Failed to parse SMILES string: {rdkit_note}"
                logger.error(error)
                return None, error
            method = "RDKit" if rdkit_note is None else f"RDKit ({rdkit_note})"

        # For metal complexes: convert bonds to dative and recalculate hydrogens
        # This fixes the issue where metal coordination bonds are counted towards ligand valence
        if has_metal:
            try:
                # Remove existing explicit H atoms first (stk may have added incorrect ones)
                mol = Chem.RemoveHs(mol)

                # Convert single bonds to metals to dative bonds
                mol = _convert_metal_bonds_to_dative(mol)
                mol.UpdatePropertyCache(strict=False)

                # For simple organometallics (e.g., Grignard), add Hs but strip
                # any hydrogens attached to metals/halogens.
                if _is_simple_organometallic(smiles):
                    mol = Chem.AddHs(mol, addCoords=True)
                    mol = _strip_h_on_metal_halogen(mol)
                    mol = _fix_organometallic_carbon_h(mol)
                else:
                    # Now add hydrogens with correct valence calculation
                    mol = Chem.AddHs(mol, addCoords=True)
            except Exception as e:
                logger.warning(f"Could not fix metal coordination hydrogens: {e}")
                # Fallback: just try to add Hs
                try:
                    if _is_simple_organometallic(smiles):
                        mol = Chem.AddHs(mol, addCoords=True)
                        mol = _strip_h_on_metal_halogen(mol)
                        mol = _fix_organometallic_carbon_h(mol)
                    else:
                        mol = Chem.AddHs(mol, addCoords=True)
                except Exception as e2:
                    if "Explicit valence" in str(e2):
                        legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                        if legacy_err is None and legacy_xyz:
                            if output_path:
                                Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                                logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                            return legacy_xyz, None
                    pass
        else:
            # Non-metal molecules: just add hydrogens normally
            try:
                mol = Chem.AddHs(mol, addCoords=True)
            except Exception as e:
                if "Explicit valence" in str(e):
                    legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                    if legacy_err is None and legacy_xyz:
                        if output_path:
                            Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                            logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                        return legacy_xyz, None
                pass

        # Generate 3D coordinates using ETKDG method
        # ETKDGv3 is the latest version with improved accuracy
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  # For reproducibility

        result = AllChem.EmbedMolecule(mol, params)

        if result != 0:
            # Try with random coordinates if ETKDG fails
            logger.warning("ETKDG embedding failed, trying random coordinates")
            result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        if result != 0:
            # Fallback to legacy embedding logic (used in older dashboards)
            logger.warning("ETKDG embedding failed, trying legacy SMILES embedding fallback")
            legacy_xyz, legacy_err = _smiles_to_xyz_legacy(smiles, has_metal=has_metal)
            if legacy_err is None and legacy_xyz:
                if output_path:
                    Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                    logger.info(f"Converted SMILES to XYZ using legacy fallback: {output_path}")
                return legacy_xyz, None
            error = f"Failed to generate 3D coordinates for SMILES: {smiles}"
            if legacy_err:
                error = f"{error} (legacy fallback failed: {legacy_err})"
            logger.error(error)
            return None, error

        # Optimize geometry with UFF force field for better starting structure
        # Note: UFF often fails for metal complexes, which is OK - the raw ETKDG
        # geometry is good enough for further optimization with xTB/GOAT/ORCA
        if not has_metal:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                logger.debug("UFF optimization successful")
            except Exception as e:
                logger.info(f"UFF optimization skipped (common for metal complexes): {e}")

        # Convert to XYZ format
        xyz_content = _mol_to_xyz(mol)

        # Check if molecule contains metals and warn about geometry quality
        has_metals = any(atom.GetAtomicNum() in range(21, 31) or  # 3d metals: Sc-Zn
                        atom.GetAtomicNum() in range(39, 49) or  # 4d metals: Y-Cd
                        atom.GetAtomicNum() in range(57, 81)     # Lanthanides + 5d metals
                        for atom in mol.GetAtoms())

        if has_metals:
            logger.warning(
                "SMILES contains metal atoms - generated geometry may be unrealistic! "
                "Coordination geometries from RDKit are rough approximations. "
                "STRONGLY RECOMMENDED: Use GOAT or manual geometry optimization before ORCA calculations."
            )

        if output_path:
            Path(output_path).write_text(xyz_content, encoding='utf-8')
            logger.info(f"Converted SMILES to XYZ using {method}: {output_path}")

        return xyz_content, None

    except Exception as e:
        msg = str(e)
        if "kekulize" in msg or "Explicit valence" in msg:
            legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
            if legacy_err is None and legacy_xyz:
                if output_path:
                    Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                    logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                return legacy_xyz, None
        error = f"Error converting SMILES to XYZ: {e}"
        logger.error(error, exc_info=True)
        return None, error


def _smiles_to_xyz_legacy(smiles: str, has_metal: bool) -> Tuple[Optional[str], Optional[str]]:
    """Legacy embedding fallback (matches older dashboard behavior)."""
    if not RDKIT_AVAILABLE:
        return None, "RDKit is not installed"

    stk_error = None
    mol = None

    # Try stk first for metal complexes
    if has_metal and STK_AVAILABLE:
        try:
            bb = stk.BuildingBlock(smiles)
            mol = bb.to_rdkit_mol()
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
        except Exception as e:
            stk_error = str(e)
            mol = None

    # RDKit fallback
    if mol is None:
        mol, rdkit_note = mol_from_smiles_rdkit(smiles, allow_metal=has_metal)
        if mol is None:
            if stk_error:
                return None, f"RDKit: {rdkit_note}; stk: {stk_error}"
            return None, rdkit_note

        try:
            mol = Chem.AddHs(mol)
        except Exception:
            pass

        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useRandomCoords = True
            try:
                params.maxAttempts = 100
            except Exception:
                pass

            try:
                result = AllChem.EmbedMolecule(mol, params, maxAttempts=50)
            except TypeError:
                result = AllChem.EmbedMolecule(mol, params)

            if result == -1:
                params2 = AllChem.ETKDGv3()
                params2.randomSeed = 42
                params2.useRandomCoords = True
                try:
                    params2.maxAttempts = 200
                except Exception:
                    pass
                try:
                    result = AllChem.EmbedMolecule(mol, params2, maxAttempts=200)
                except TypeError:
                    result = AllChem.EmbedMolecule(mol, params2)
                if result == -1:
                    return None, "Legacy embed failed to generate 3D structure"
        except Exception as e:
            return None, f"Legacy RDKit error: {e}"

    try:
        xyz_content = _mol_to_xyz(mol)
        return xyz_content, None
    except Exception as e:
        return None, f"Legacy coordinate error: {e}"


def _smiles_to_xyz_unsanitized_fallback(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Last-resort fallback: embed without sanitization (handles valence/kekulize errors)."""
    if not RDKIT_AVAILABLE:
        return None, "RDKit is not installed"

    try:
        normalized = _normalize_metal_smiles(smiles)
        smiles_try = normalized or smiles
        try:
            params = Chem.SmilesParserParams()
            params.sanitize = False
            params.removeHs = False
            params.strictParsing = False
            mol = Chem.MolFromSmiles(smiles_try, params)
        except Exception:
            mol = Chem.MolFromSmiles(smiles_try, sanitize=False)
        if mol is None:
            return None, "Failed to parse SMILES (no sanitize)"
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass

        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.useRandomCoords = True
        try:
            params.maxAttempts = 200
        except Exception:
            pass

        try:
            result = AllChem.EmbedMolecule(mol, params, maxAttempts=200)
        except TypeError:
            result = AllChem.EmbedMolecule(mol, params)

        if result != 0:
            try:
                result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
            except TypeError:
                result = AllChem.EmbedMolecule(mol, useRandomCoords=True)

        if result != 0:
            return None, "Failed to generate 3D coordinates (unsanitized)"

        # Try to add hydrogens after embedding unless this is a neutral Ni/Co+[N] case
        if not _prefer_no_sanitize(smiles):
            try:
                mol_h = Chem.AddHs(mol)
                # Re-embed to place H coordinates
                try:
                    result_h = AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=42)
                except TypeError:
                    result_h = AllChem.EmbedMolecule(mol_h, useRandomCoords=True)
                if result_h == 0:
                    mol = mol_h
            except Exception:
                pass

        xyz_content = _mol_to_xyz(mol)
        return xyz_content, None
    except Exception as e:
        # Extra-permissive fallback for explicit valence errors
        if "Explicit valence" in str(e):
            try:
                normalized = _normalize_metal_smiles(smiles)
                smiles_try = normalized or smiles
                mol = Chem.MolFromSmiles(smiles_try, sanitize=False)
                if mol is None:
                    return None, "Failed to parse SMILES (no sanitize)"
                # Disable implicit H handling to avoid valence checks
                for atom in mol.GetAtoms():
                    atom.SetNoImplicit(True)
                    atom.SetNumExplicitHs(0)
                try:
                    result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
                except TypeError:
                    result = AllChem.EmbedMolecule(mol, useRandomCoords=True)
                if result == 0:
                    xyz_content = _mol_to_xyz(mol)
                    return xyz_content, None
            except Exception:
                pass
        return None, f"Unsanitized fallback error: {e}"


def _mol_to_xyz(mol) -> str:
    """Convert RDKit molecule to DELFIN coordinate format.

    DELFIN expects coordinates without XYZ header (no atom count, no comment line).
    Just element symbols followed by x, y, z coordinates.

    Args:
        mol: RDKit molecule with 3D coordinates

    Returns:
        Coordinate string in DELFIN format (no header)
    """
    # Strip spurious H on coordinated P/As (affects legacy/unsanitized paths)
    if any(a.GetSymbol() in _METAL_SET for a in mol.GetAtoms()):
        mol = _strip_h_on_coordinated_p(mol)

    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()

    # Atom coordinates only (no XYZ header for DELFIN)
    lines = []
    for i in range(num_atoms):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        lines.append(f"{symbol:4s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")

    return '\n'.join(lines) + '\n'


def _mol_to_xyz_conformer(mol, conf_id: int) -> str:
    """Convert a specific conformer of an RDKit molecule to DELFIN coordinate format.

    Like ``_mol_to_xyz`` but uses *conf_id* instead of the first conformer.
    """
    conf = mol.GetConformer(conf_id)
    num_atoms = mol.GetNumAtoms()
    lines = []
    for i in range(num_atoms):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        lines.append(f"{symbol:4s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")
    return '\n'.join(lines) + '\n'


def convert_input_if_smiles(input_path: Path) -> Tuple[bool, Optional[str]]:
    """Check if input file contains SMILES and convert if needed.

    This function is called by the pipeline to automatically handle SMILES input.

    Args:
        input_path: Path to input file to check

    Returns:
        Tuple of (was_converted, error_message)
        - was_converted: True if file was SMILES and was converted
        - error_message: Error description if conversion failed, None if success or not SMILES
    """
    if not input_path.exists():
        return False, f"Input file does not exist: {input_path}"

    try:
        content = input_path.read_text(encoding='utf-8', errors='ignore')
    except Exception as e:
        return False, f"Could not read input file: {e}"

    if not is_smiles_string(content):
        return False, None

    # Extract SMILES (first non-comment line)
    smiles = None
    for line in content.split('\n'):
        line = line.strip()
        if line and not line.startswith('#') and not line.startswith('*'):
            smiles = line
            break

    if not smiles:
        return False, "No SMILES string found in input file"

    logger.info(f"Detected SMILES string in {input_path.name}: {smiles}")

    # Convert SMILES to XYZ
    xyz_content, error = smiles_to_xyz(smiles)

    if error:
        return False, error

    # Write XYZ content back to input file (replacing SMILES)
    try:
        input_path.write_text(xyz_content, encoding='utf-8')
        logger.info(f"Converted SMILES to XYZ coordinates in {input_path}")
        return True, None
    except Exception as e:
        return False, f"Could not write converted coordinates: {e}"
