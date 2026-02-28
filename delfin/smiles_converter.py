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

# Try to import Open Babel (has full UFF parameters for transition metals)
try:
    from openbabel import pybel
    OPENBABEL_AVAILABLE = True
    try:
        # Reduce Open Babel console noise in notebook/voila runs.
        # Keep only error-level output (suppress repeated warning spam like
        # "Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders").
        pybel.ob.obErrorLog.SetOutputLevel(pybel.ob.obError)
    except Exception:
        pass
except ImportError:
    OPENBABEL_AVAILABLE = False

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

# Atomic numbers of metal elements — used to filter metal-containing SSSR rings
# in the roundtrip check so that OB-perceived M-L rings (from short distances)
# do not skew the comparison against the SMILES organic-ring count.
try:
    if RDKIT_AVAILABLE:
        _pt = Chem.GetPeriodicTable()
        _METAL_ATOMICNUMS: frozenset = frozenset(
            _pt.GetAtomicNumber(sym) for sym in _METALS
            if _pt.GetAtomicNumber(sym) > 0
        )
    else:
        raise RuntimeError("rdkit unavailable")
except Exception:
    # Hardcoded fallback covering all metals in _METALS list
    _METAL_ATOMICNUMS = frozenset([
        3, 11, 19, 37, 55,          # Li Na K Rb Cs
        4, 12, 20, 38, 56,          # Be Mg Ca Sr Ba
        13, 31, 49, 81,             # Al Ga In Tl
        50, 82, 83, 84,             # Sn Pb Bi Po
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30,   # Sc-Zn
        39, 40, 41, 42, 43, 44, 45, 46, 47, 48,   # Y-Cd
        72, 73, 74, 75, 76, 77, 78, 79, 80,       # Hf-Hg
        57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,  # La-Lu
        89, 90, 91, 92, 93, 94,     # Ac Th Pa U Np Pu
    ])


# ---------------------------------------------------------------------------
# Covalent radii (Pyykkö 2009, single-bond values in Å)
# Used as fallback for M-L bond length estimation when no specific value
# is available in _METAL_LIGAND_BOND_LENGTHS.
# ---------------------------------------------------------------------------
_COVALENT_RADII: Dict[str, float] = {
    # 3d transition metals
    'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39,
    'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
    # 4d transition metals
    'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47,
    'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
    # 5d transition metals
    'La': 2.07, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62, 'Re': 1.51,
    'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
    # Lanthanides
    'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
    'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92,
    'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87,
    # Actinides
    'Ac': 2.15, 'Th': 2.06, 'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
    # Main-group metals
    'Li': 1.28, 'Na': 1.66, 'K': 2.03, 'Rb': 2.20, 'Cs': 2.44,
    'Be': 0.96, 'Mg': 1.41, 'Ca': 1.76, 'Sr': 1.95, 'Ba': 2.15,
    'Al': 1.21, 'Ga': 1.22, 'In': 1.42, 'Tl': 1.45,
    'Sn': 1.39, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40,
    # Non-metals (common donor atoms and halogens)
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Se': 1.20, 'Br': 1.20,
    'Te': 1.38, 'I': 1.39,
}

# ---------------------------------------------------------------------------
# Typical M-L bond lengths (Å) from crystallographic databases.
# Keyed as (metal_symbol, donor_element) → distance.
# ---------------------------------------------------------------------------
_METAL_LIGAND_BOND_LENGTHS: Dict[Tuple[str, str], float] = {
    # Iridium
    ('Ir', 'C'): 2.02, ('Ir', 'N'): 2.15, ('Ir', 'O'): 2.12,
    ('Ir', 'P'): 2.30, ('Ir', 'Cl'): 2.38, ('Ir', 'S'): 2.38,
    # Ruthenium
    ('Ru', 'C'): 2.00, ('Ru', 'N'): 2.10, ('Ru', 'O'): 2.08,
    ('Ru', 'P'): 2.30, ('Ru', 'Cl'): 2.40, ('Ru', 'S'): 2.35,
    # Rhodium
    ('Rh', 'C'): 2.00, ('Rh', 'N'): 2.10, ('Rh', 'O'): 2.05,
    ('Rh', 'P'): 2.28, ('Rh', 'Cl'): 2.38,
    # Platinum
    ('Pt', 'C'): 2.00, ('Pt', 'N'): 2.05, ('Pt', 'O'): 2.02,
    ('Pt', 'P'): 2.25, ('Pt', 'Cl'): 2.30, ('Pt', 'S'): 2.30,
    ('Pt', 'Br'): 2.43,
    # Palladium
    ('Pd', 'C'): 2.00, ('Pd', 'N'): 2.05, ('Pd', 'O'): 2.02,
    ('Pd', 'P'): 2.28, ('Pd', 'Cl'): 2.30, ('Pd', 'S'): 2.30,
    # Gold
    ('Au', 'C'): 2.00, ('Au', 'N'): 2.05, ('Au', 'P'): 2.28,
    ('Au', 'Cl'): 2.28, ('Au', 'S'): 2.30,
    # Iron
    ('Fe', 'C'): 1.91, ('Fe', 'N'): 1.97, ('Fe', 'O'): 2.00,
    ('Fe', 'P'): 2.20, ('Fe', 'Cl'): 2.28, ('Fe', 'S'): 2.30,
    # Cobalt
    ('Co', 'C'): 1.90, ('Co', 'N'): 1.95, ('Co', 'O'): 1.95,
    ('Co', 'P'): 2.18, ('Co', 'Cl'): 2.25,
    # Nickel
    ('Ni', 'C'): 1.88, ('Ni', 'N'): 1.90, ('Ni', 'O'): 1.88,
    ('Ni', 'P'): 2.15, ('Ni', 'Cl'): 2.20,
    # Copper
    ('Cu', 'C'): 1.95, ('Cu', 'N'): 2.00, ('Cu', 'O'): 1.97,
    ('Cu', 'P'): 2.20, ('Cu', 'Cl'): 2.25, ('Cu', 'S'): 2.30,
    # Zinc
    ('Zn', 'N'): 2.05, ('Zn', 'O'): 2.10, ('Zn', 'S'): 2.30,
    ('Zn', 'Cl'): 2.20,
    # Manganese
    ('Mn', 'N'): 2.05, ('Mn', 'O'): 2.10, ('Mn', 'Cl'): 2.35,
    # Chromium
    ('Cr', 'N'): 2.05, ('Cr', 'O'): 2.00, ('Cr', 'Cl'): 2.30,
    # Titanium
    ('Ti', 'N'): 2.10, ('Ti', 'O'): 1.95, ('Ti', 'Cl'): 2.30,
    # Zirconium
    ('Zr', 'N'): 2.20, ('Zr', 'O'): 2.10, ('Zr', 'Cl'): 2.45,
    # Molybdenum
    ('Mo', 'N'): 2.15, ('Mo', 'O'): 2.05, ('Mo', 'Cl'): 2.40,
    ('Mo', 'S'): 2.40,
    # Tungsten
    ('W', 'N'): 2.15, ('W', 'O'): 2.05, ('W', 'Cl'): 2.40,
    # Rhenium
    ('Re', 'N'): 2.15, ('Re', 'O'): 2.05, ('Re', 'Cl'): 2.40,
    # Osmium
    ('Os', 'N'): 2.10, ('Os', 'O'): 2.05, ('Os', 'Cl'): 2.35,
    ('Os', 'P'): 2.30,
    # Silver
    ('Ag', 'N'): 2.15, ('Ag', 'O'): 2.20, ('Ag', 'P'): 2.35,
    ('Ag', 'S'): 2.40, ('Ag', 'Cl'): 2.30,
    # Lanthanides (representative: La, Ce, Gd, Lu)
    ('La', 'N'): 2.55, ('La', 'O'): 2.45, ('La', 'Cl'): 2.75,
    ('Ce', 'N'): 2.50, ('Ce', 'O'): 2.40, ('Ce', 'Cl'): 2.70,
    ('Gd', 'N'): 2.45, ('Gd', 'O'): 2.35, ('Gd', 'Cl'): 2.65,
    ('Lu', 'N'): 2.30, ('Lu', 'O'): 2.20, ('Lu', 'Cl'): 2.50,
    # Actinides
    ('U', 'N'): 2.55, ('U', 'O'): 2.30, ('U', 'Cl'): 2.65,
    ('Th', 'N'): 2.55, ('Th', 'O'): 2.35, ('Th', 'Cl'): 2.70,
}


def _get_ml_bond_length(metal_symbol: str, donor_symbol: str) -> float:
    """Return estimated M-L bond length in Å.

    Lookup order:
    1. Specific (metal, donor) pair in _METAL_LIGAND_BOND_LENGTHS
    2. Sum of covalent radii + 0.5 Å (coordination bond correction)
    3. Default 2.0 Å
    """
    key = (metal_symbol, donor_symbol)
    if key in _METAL_LIGAND_BOND_LENGTHS:
        return _METAL_LIGAND_BOND_LENGTHS[key]
    r_m = _COVALENT_RADII.get(metal_symbol)
    r_d = _COVALENT_RADII.get(donor_symbol)
    if r_m is not None and r_d is not None:
        return r_m + r_d + 0.5
    return 2.0


# ---------------------------------------------------------------------------
# Preferred CN=4 geometry per metal (d8 → square-planar, others → tetrahedral)
# ---------------------------------------------------------------------------
_PREFERRED_CN4_GEOMETRY: Dict[str, str] = {
    # d8 metals (strong square-planar preference)
    'Ni': 'SQ', 'Pd': 'SQ', 'Pt': 'SQ', 'Au': 'SQ',
    'Rh': 'SQ', 'Ir': 'SQ',
    # Tetrahedral preference
    'Zn': 'TH', 'Cd': 'TH', 'Hg': 'TH',
    'Co': 'TH', 'Fe': 'TH', 'Mn': 'TH', 'Cu': 'TH',
    'Ti': 'TH', 'Zr': 'TH', 'Hf': 'TH',
    'V': 'TH', 'Cr': 'TH',
    'Al': 'TH', 'Ga': 'TH', 'In': 'TH',
}


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

    Returns True only for [N] (neutral), [N-] (anionic), or ring-closing N
    donors (pyridyl-type, e.g. N1=, N%10=).  [N+] is intentionally excluded:
    complexes like Ir(ppy)3 use [N+] notation and must go through the full
    dative-bond path so that ETKDG can generate fac/mer conformers correctly.

    Matches both neutral ([Metal]) and charged ([Metal+2], [Metal-]) metal
    notations, as well as stereochemical variants ([Metal@@+2], etc.).
    """
    metal_pattern = '|'.join(re.escape(m) for m in _METALS)
    has_metal = bool(re.search(rf'\[(?:{metal_pattern})(?:@+|@@)?(?:[+-]\d*)?\]', smiles))

    # [N+] intentionally excluded — routes Ir(ppy)3 etc. to full-dative path.
    has_coord_n = (
        '[N]' in smiles
        or '[N-]' in smiles
        or bool(re.search(r'N\d+=', smiles))
        or bool(re.search(r'N%\d+=', smiles))
    )

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

    # Strategy 5: Open Babel / Avogadro-like 3D generation.
    if OPENBABEL_AVAILABLE:
        xyz_blocks, ob_err = _openbabel_generate_conformer_xyz(smiles, num_confs=1)
        if xyz_blocks:
            xyz_content = xyz_blocks[0]
            if output_path:
                Path(output_path).write_text(xyz_content, encoding='utf-8')
                logger.info("Converted SMILES to XYZ using Open Babel: %s", output_path)
            return xyz_content, None
        if ob_err:
            errors.append(f"openbabel({smiles[:30]}...): {ob_err}")

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


def _manual_metal_embed(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Last-resort fallback: manually construct coordinates for metal complexes.

    Places the metal at the origin and coordinating atoms at idealized
    positions (octahedral for 6-coord, tetrahedral for 4-coord, etc.).
    Remaining ligand atoms are placed along bond directions using a simple
    BFS traversal.  The geometry is rough but usable as a GOAT/xTB input.

    For 4-coordinate metals, both tetrahedral and square-planar geometries are
    tried, OB UFF is applied to each, and the better-scoring geometry is kept.
    """
    if not RDKIT_AVAILABLE:
        return None, "RDKit is not installed"

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None, "Failed to parse SMILES"

        # Suppress valence checks
        for atom in mol.GetAtoms():
            atom.SetNoImplicit(True)

        # Try to add H for non-metal atoms
        try:
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    atom.SetNoImplicit(False)
            mol.UpdatePropertyCache(strict=False)
            mol = Chem.AddHs(mol)
        except Exception:
            pass

        # Find metal centers
        metal_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
        if not metal_indices:
            return None, "No metal found"

        # Idealized coordination vectors for common coordination numbers
        _COORD_VECTORS_BASE = {
            2: [(1, 0, 0), (-1, 0, 0)],
            3: [(1, 0, 0), (-0.5, 0.866, 0), (-0.5, -0.866, 0)],
            4: [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)],  # tetrahedral
            5: [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -0.5, 0.866), (0, -0.5, -0.866)],
            6: [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)],
            7: [(0, 0, 1), (0, 0, -1), (1, 0, 0), (0.309, 0.951, 0),
                (-0.809, 0.588, 0), (-0.809, -0.588, 0), (0.309, -0.951, 0)],
            8: [(1.414, 0, 0.8), (0, 1.414, 0.8), (-1.414, 0, 0.8), (0, -1.414, 0.8),
                (1, 1, -0.8), (-1, 1, -0.8), (-1, -1, -0.8), (1, -1, -0.8)],
        }

        n_atoms = mol.GetNumAtoms()

        def _build_xyz_for_4coord_vecs(override_4coord):
            """Build a raw XYZ string using given vectors for 4-coord metals."""
            local_coords = [(0.0, 0.0, 0.0)] * n_atoms
            local_placed: set = set()

            for mi in metal_indices:
                local_coords[mi] = (0.0, 0.0, 0.0)
                local_placed.add(mi)

                neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(mi).GetNeighbors()]
                n_coord = len(neighbors)
                if n_coord == 4:
                    vectors = override_4coord
                else:
                    vectors = _COORD_VECTORS_BASE.get(n_coord, _COORD_VECTORS_BASE.get(6, []))

                metal_sym = mol.GetAtomWithIdx(mi).GetSymbol()
                for i, nbr_idx in enumerate(neighbors):
                    donor_sym = mol.GetAtomWithIdx(nbr_idx).GetSymbol()
                    bond_len = _get_ml_bond_length(metal_sym, donor_sym)
                    if i < len(vectors):
                        vx, vy, vz = vectors[i]
                        mag = math.sqrt(vx**2 + vy**2 + vz**2)
                        if mag > 1e-8:
                            vx = vx/mag * bond_len
                            vy = vy/mag * bond_len
                            vz = vz/mag * bond_len
                    else:
                        angle = 2 * math.pi * i / n_coord
                        vx = bond_len * math.cos(angle)
                        vy = bond_len * math.sin(angle)
                        vz = 0.0
                    local_coords[nbr_idx] = (vx, vy, vz)
                    local_placed.add(nbr_idx)

            # BFS to place remaining atoms
            bond_len_default = 1.4
            queue = list(local_placed)
            while queue:
                current = queue.pop(0)
                cx, cy, cz = local_coords[current]
                atom = mol.GetAtomWithIdx(current)
                unplaced_nbrs = [
                    n.GetIdx() for n in atom.GetNeighbors()
                    if n.GetIdx() not in local_placed
                ]
                for k, nbr_idx in enumerate(unplaced_nbrs):
                    dx, dy, dz = 0.0, 0.0, 0.0
                    for other in atom.GetNeighbors():
                        oi = other.GetIdx()
                        if oi in local_placed and oi != nbr_idx:
                            ox, oy, oz = local_coords[oi]
                            dx += cx - ox
                            dy += cy - oy
                            dz += cz - oz
                    mag = math.sqrt(dx**2 + dy**2 + dz**2)
                    if mag < 1e-8:
                        dx, dy, dz = 1.0 + 0.1 * k, 0.3 * k, 0.0
                        mag = math.sqrt(dx**2 + dy**2 + dz**2)
                    dx = dx/mag * bond_len_default
                    dy = dy/mag * bond_len_default
                    dz = dz/mag * bond_len_default
                    local_coords[nbr_idx] = (cx + dx, cy + dy, cz + dz)
                    local_placed.add(nbr_idx)
                    queue.append(nbr_idx)

            lines = []
            for i in range(n_atoms):
                atom = mol.GetAtomWithIdx(i)
                x, y, z = local_coords[i]
                lines.append(f"{atom.GetSymbol():4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            return '\n'.join(lines) + '\n'

        # Check whether any metal is 4-coordinate (warrants dual geometry trial)
        has_4coord = any(
            len([nbr.GetIdx() for nbr in mol.GetAtomWithIdx(mi).GetNeighbors()]) == 4
            for mi in metal_indices
        )

        if has_4coord:
            # Try both tetrahedral and square-planar vectors
            tetra_vecs = [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)]
            sq_vecs = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]

            # Determine preferred geometry from metal identity
            metal_syms = [mol.GetAtomWithIdx(mi).GetSymbol() for mi in metal_indices]
            pref_geom = _PREFERRED_CN4_GEOMETRY.get(metal_syms[0], 'SQ')

            xyz_tetra = _build_xyz_for_4coord_vecs(tetra_vecs)
            xyz_sq = _build_xyz_for_4coord_vecs(sq_vecs)

            # OB UFF refinement for both
            xyz_tetra = _optimize_xyz_openbabel(xyz_tetra)
            xyz_sq = _optimize_xyz_openbabel(xyz_sq)

            # Score via temporary RDKit conformers
            try:
                mol_tmp = Chem.RWMol(mol)
                mol_tmp.RemoveAllConformers()
                conf_t = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_tetra)
                conf_s = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_sq)
                if conf_t is not None and conf_s is not None:
                    cid_t = mol_tmp.AddConformer(conf_t, assignId=True)
                    cid_s = mol_tmp.AddConformer(conf_s, assignId=True)
                    score_t = _geometry_quality_score(mol_tmp.GetMol(), cid_t)
                    score_s = _geometry_quality_score(mol_tmp.GetMol(), cid_s)
                    xyz_str = xyz_tetra if score_t <= score_s else xyz_sq
                else:
                    # Use metal-preferred geometry as fallback
                    xyz_str = xyz_tetra if pref_geom == 'TH' else xyz_sq
            except Exception:
                xyz_str = xyz_tetra if pref_geom == 'TH' else xyz_sq
        else:
            # Standard path for non-4-coord metals (use default vectors)
            xyz_str = _build_xyz_for_4coord_vecs(
                _COORD_VECTORS_BASE.get(4, [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)])
            )
            xyz_str = _optimize_xyz_openbabel(xyz_str)

        return xyz_str, None

    except Exception as e:
        return None, f"Manual metal embed error: {e}"


# ---------------------------------------------------------------------------
# Open Babel helper functions (Avogadro-equivalent 3D generation pipeline)
# ---------------------------------------------------------------------------

def _normalize_conversion_backend(backend: Optional[str]) -> str:
    """Normalize conversion backend token to ``'rdkit'`` or ``'avogadro'``."""
    if backend is None:
        return "rdkit"
    token = str(backend).strip().lower()
    aliases = {
        "rdkit": "rdkit",
        "default": "rdkit",
        "avogadro": "avogadro",
        "openbabel": "avogadro",
        "obabel": "avogadro",
        "ob": "avogadro",
    }
    return aliases.get(token, "rdkit")


def _openbabel_smiles_variants(smiles: str) -> List[str]:
    """Return unique SMILES variants used for robust Open Babel parsing."""
    variants = [smiles]
    normalized = _normalize_metal_smiles(smiles)
    denormalized = _denormalize_metal_smiles(smiles)
    if normalized and normalized not in variants:
        variants.append(normalized)
    if denormalized and denormalized not in variants:
        variants.append(denormalized)
    return variants


def _pick_openbabel_forcefield(smiles: str) -> str:
    """Pick an Open Babel force field close to Avogadro defaults.

    UFF has full parameters for transition metals; MMFF94 is preferred for
    purely organic systems.
    """
    if not OPENBABEL_AVAILABLE:
        return "uff"
    if contains_metal(smiles):
        return "uff"
    if "mmff94" in pybel._forcefields:
        return "mmff94"
    return "uff"


def _obmol_to_delfin_xyz(ob_mol) -> str:
    """Convert an Open Babel OBMol to DELFIN XYZ lines (no header)."""
    lines = []
    for ob_atom in pybel.ob.OBMolAtomIter(ob_mol):
        symbol = pybel.ob.GetSymbol(ob_atom.GetAtomicNum())
        x, y, z = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()
        lines.append(f"{symbol:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"


def _openbabel_generate_conformer_xyz(
    smiles: str,
    *,
    num_confs: int = 200,
    forcefield: Optional[str] = None,
    rotor_steps: int = 120,
    localopt_steps: int = 250,
    optimize: bool = True,
) -> Tuple[List[str], Optional[str]]:
    """Generate 3D conformers using Open Babel in an Avogadro-like workflow.

    Pipeline (mirrors Avogadro's ``gen3d`` quality):
    1. Fragment-based 3D initialization via ``make3D`` (uses CSD crystal data)
    2. ``WeightedRotorSearch`` — generates a pool of ``num_confs`` conformers
    3. Per-conformer ``ConjugateGradients`` local optimization

    Tries multiple SMILES variants (original / normalized / denormalized) for
    robustness with metal complexes.  Duplicates are removed via text-hash.

    Returns:
        (xyz_blocks, error) — ``xyz_blocks`` is a list of DELFIN-format XYZ
        strings (one per unique conformer); ``error`` is set only when *all*
        variants failed.
    """
    if not OPENBABEL_AVAILABLE:
        return [], "Open Babel is not installed"

    target = max(1, int(num_confs))
    errors: List[str] = []
    chosen_ff = (forcefield or "").strip().lower() or _pick_openbabel_forcefield(smiles)

    for smi in _openbabel_smiles_variants(smiles):
        try:
            ob_mol = pybel.readstring("smi", smi)
        except Exception as exc:
            errors.append(f"parse({smi[:30]}...): {exc}")
            continue

        if ob_mol.OBMol.NumAtoms() <= 0:
            errors.append(f"parse({smi[:30]}...): empty molecule")
            continue

        ff_name = chosen_ff if chosen_ff in pybel._forcefields else "uff"
        if ff_name not in pybel._forcefields:
            errors.append("no usable Open Babel forcefield")
            continue

        make3d_steps = max(50, min(200, localopt_steps // 2)) if optimize else 0
        try:
            ob_mol.make3D(forcefield=ff_name, steps=make3d_steps)
        except Exception as exc:
            errors.append(f"make3D({ff_name}): {exc}")
            # Fallback once to UFF if initial forcefield failed.
            if ff_name != "uff" and "uff" in pybel._forcefields:
                try:
                    ff_name = "uff"
                    fallback_steps = max(50, min(200, localopt_steps // 2)) if optimize else 1
                    ob_mol.make3D(forcefield=ff_name, steps=fallback_steps)
                except Exception as exc2:
                    errors.append(f"make3D(uff): {exc2}")
                    continue
            else:
                continue

        if not optimize:
            xyz_text = _obmol_to_delfin_xyz(ob_mol.OBMol)
            if xyz_text.strip():
                return [xyz_text], None
            errors.append(f"conformers({smi[:30]}...): empty geometry")
            continue

        ff = pybel._forcefields.get(ff_name)
        if ff is None:
            errors.append(f"forcefield({ff_name}): unavailable")
            continue

        try:
            if not ff.Setup(ob_mol.OBMol):
                errors.append(f"forcefield({ff_name}): setup failed")
                continue
        except Exception as exc:
            errors.append(f"forcefield({ff_name}): {exc}")
            continue

        try:
            ff.WeightedRotorSearch(target, max(25, int(rotor_steps)))
            ff.GetConformers(ob_mol.OBMol)
        except Exception as exc:
            logger.debug("Open Babel conformer search failed, using initial 3D geometry: %s", exc)

        num_ob_confs = int(ob_mol.OBMol.NumConformers() or 0)
        if num_ob_confs <= 0:
            num_ob_confs = 1

        xyz_blocks: List[str] = []
        seen: set = set()
        for conf_idx in range(min(target, num_ob_confs)):
            try:
                ob_mol.OBMol.SetConformer(conf_idx)
            except Exception:
                pass

            try:
                if ff.Setup(ob_mol.OBMol):
                    ff.ConjugateGradients(max(50, int(localopt_steps)))
                    ff.GetCoordinates(ob_mol.OBMol)
            except Exception:
                # Keep current coordinates if local optimization fails.
                pass

            xyz_text = _obmol_to_delfin_xyz(ob_mol.OBMol)
            key = "\n".join(line.strip() for line in xyz_text.splitlines() if line.strip())
            if key in seen:
                continue
            seen.add(key)
            xyz_blocks.append(xyz_text)

        if xyz_blocks:
            return xyz_blocks, None

        errors.append(f"conformers({smi[:30]}...): none generated")

    if errors:
        return [], f"Open Babel conversion failed: {'; '.join(errors)}"
    return [], "Open Babel conversion failed"


def _xyz_to_rdkit_conformer(mol, xyz_delfin: str):
    """Build an RDKit Conformer from DELFIN XYZ lines if atom order matches.

    Validates that the number of atom lines equals ``mol.GetNumAtoms()`` and
    that element symbols appear in the same order.  Returns ``None`` if
    validation fails or any coordinate cannot be parsed.
    """
    if not RDKIT_AVAILABLE:
        return None

    lines = [ln.strip() for ln in xyz_delfin.splitlines() if ln.strip()]
    if len(lines) != mol.GetNumAtoms():
        return None

    conf = Chem.Conformer(mol.GetNumAtoms())
    for idx, line in enumerate(lines):
        parts = line.split()
        if len(parts) < 4:
            return None
        atom = mol.GetAtomWithIdx(idx)
        if parts[0] != atom.GetSymbol():
            return None
        try:
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
        except ValueError:
            return None
        conf.SetAtomPosition(idx, (x, y, z))
    return conf


def _inject_openbabel_conformers_into_mol(mol, xyz_blocks: List[str]) -> List[int]:
    """Populate *mol* with conformers parsed from Open Babel XYZ blocks.

    Clears all existing RDKit conformers and injects the OB-generated
    conformers.  Skips any xyz block whose atom count or symbol order does
    not match *mol* (atom-ordering mismatch between OB and RDKit).

    Returns the list of assigned conformer IDs.
    """
    if not xyz_blocks:
        return []
    mol.RemoveAllConformers()
    conf_ids: List[int] = []
    for xyz_text in xyz_blocks:
        conf = _xyz_to_rdkit_conformer(mol, xyz_text)
        if conf is None:
            continue
        conf_ids.append(int(mol.AddConformer(conf, assignId=True)))
    return conf_ids


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
                # If only_elements is set, skip elements not in the set
                if only_elements and ligand_atom.GetSymbol() not in only_elements:
                    continue
                bonds_to_convert.append((
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    is_metal_1  # True if atom1 is the metal
                ))

    if not bonds_to_convert:
        return mol

    # Track which atoms are dative donors (coordinating atoms whose H
    # count should be preserved from the SMILES, not recalculated).
    dative_donor_indices = set()

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

    # Reset atom properties to allow recalculation of implicit hydrogens.
    # Dative donors: set NoImplicit=True so AddHs() won't add spurious H
    # (trust the H count from the original SMILES).
    # Non-coordinating neutral atoms: recalculate H normally.
    for atom in result_mol.GetAtoms():
        if atom.GetFormalCharge() == 0 and atom.GetSymbol() not in _METAL_SET:
            if atom.GetIdx() in dative_donor_indices:
                atom.SetNoImplicit(True)
            else:
                atom.SetNoImplicit(False)
                # Preserve explicit H on atoms where RDKit's implicit H
                # calculation can't recover them (e.g., tetracoordinate B in
                # pyrazolylborate/scorpionate ligands: B standard valence = 3
                # but actual valence = 4 with the B-H bond).
                if atom.GetSymbol() == 'B' and atom.GetNumExplicitHs() > 0:
                    atom.SetNoImplicit(True)
                else:
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
        except Exception:
            try:
                mol = Chem.AddHs(mol, addCoords=False)
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
                if has_metal_bond and atom.GetFormalCharge() == 0:
                    atom.SetNoImplicit(True)
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
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


def _xyz_to_canonical_smiles(xyz_delfin: str) -> Optional[str]:
    """Convert DELFIN-format XYZ to a canonical SMILES via OpenBabel.

    Metal–ligand bonds perceived by OB from interatomic distances are
    intentionally removed before canonicalisation so that the comparison
    reflects only the organic/ligand connectivity (bond perception for
    M–L is unreliable and distance-dependent).

    Returns the canonical SMILES string, or ``None`` on failure.
    """
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return None
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        ob_mol = pybel.readstring('xyz', std_xyz)

        # Remove bonds involving metal atoms — OB perceives them from
        # distance alone which is geometry-dependent and unreliable.
        metal_ob_idxs = {a.OBAtom.GetIdx() for a in ob_mol.atoms
                         if a.atomicnum in _METAL_ATOMICNUMS}
        if metal_ob_idxs:
            raw = ob_mol.OBMol
            bonds_to_del = []
            for bond in pybel.ob.OBMolBondIter(raw):
                bi = bond.GetBeginAtomIdx()
                ei = bond.GetEndAtomIdx()
                if bi in metal_ob_idxs or ei in metal_ob_idxs:
                    bonds_to_del.append(bond)
            for bond in bonds_to_del:
                raw.DeleteBond(bond)

        # Canonical SMILES (fragment-aware — disconnected parts separated by '.')
        conv = pybel.ob.OBConversion()
        conv.SetOutFormat('can')
        raw_smi = conv.WriteString(ob_mol.OBMol).strip()
        if not raw_smi:
            return None
        # Normalise: sort dot-separated fragments for stable comparison
        parts = sorted(raw_smi.split('.'))
        return '.'.join(parts)
    except Exception:
        return None


def _roundtrip_ring_count_ok(xyz_delfin: str, original_smiles: str, tolerance: int = 2) -> bool:
    """Validate a 3D structure by round-tripping coordinates → SMILES via OpenBabel.

    Only organic rings (not containing metal atoms) are compared. Metal chelate
    and coordination rings are excluded because OB bond perception from XYZ is
    unreliable for metal-ligand bonds (e.g. Ir-N ~2.1 Å), which would cause
    false rejections for complexes like fac/mer-Ir(ppy)3.
    """
    if not OPENBABEL_AVAILABLE:
        return True
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return True
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        rt_mol = pybel.readstring('xyz', std_xyz)
        # Count ORGANIC-ONLY rings from OB: exclude any ring that contains a
        # metal atom (OB may perceive M-L bonds from short atom-atom distances,
        # adding chelate rings to the SSSR that are absent in the SMILES graph).
        try:
            _metal_ob_idxs = {a.idx for a in rt_mol.atoms
                              if a.atomicnum in _METAL_ATOMICNUMS}
            rt_rings = sum(
                1 for ring in rt_mol.sssr
                if not any(i in _metal_ob_idxs for i in ring._path)
            )
        except Exception:
            # OB failed to kekulize or ring._path unavailable — be permissive.
            # Using total SSSR here would include chelate rings and cause false
            # rejections for aromatic metal complexes (e.g. Ir(ppy)3).
            return True
        # Organic-only ring count from original SMILES.
        # OB is used first (same engine as XYZ parsing → consistent ring perception).
        # RDKit SMILES parsing often fails for metal-cyclometallated SMILES where
        # all ring-closure bonds pass through the metal atom (e.g. Ir(ppy)3),
        # causing RDKit to report 0 organic rings and rejecting all conformers.
        orig_rings = None
        try:
            orig_mol_ob = pybel.readstring('smi', original_smiles)
            _metal_smi_idxs = {a.idx for a in orig_mol_ob.atoms
                               if a.atomicnum in _METAL_ATOMICNUMS}
            orig_rings = sum(
                1 for ring in orig_mol_ob.sssr
                if not any(i in _metal_smi_idxs for i in ring._path)
            )
        except Exception:
            pass
        if orig_rings is None and RDKIT_AVAILABLE:
            try:
                orig_mol_rd = Chem.MolFromSmiles(original_smiles, sanitize=False)
                if orig_mol_rd is not None:
                    try:
                        orig_mol_rd.UpdatePropertyCache(strict=False)
                    except Exception:
                        pass
                    metal_indices = {
                        a.GetIdx() for a in orig_mol_rd.GetAtoms()
                        if a.GetSymbol() in _METAL_SET
                    }
                    ring_info = orig_mol_rd.GetRingInfo()
                    orig_rings = sum(
                        1 for ring in ring_info.AtomRings()
                        if not any(idx in metal_indices for idx in ring)
                    )
            except Exception:
                pass
        if orig_rings is None:
            return True
        return abs(rt_rings - orig_rings) <= tolerance
    except Exception:
        return True


def _no_spurious_bonds(xyz_delfin: str, original_smiles: str) -> bool:
    """Return True if OB-perceived XYZ contains no clearly spurious homodiatomic bonds.

    Only checks a short whitelist of homodiatomic element pairs that are very
    unlikely to appear in typical coordination-chemistry SMILES but can be
    falsely perceived by OpenBabel from close atom distances in bad geometries:

        O-O  (peroxide artifact)
        F-F, Cl-Cl, Br-Br, I-I  (halogen-halogen artifact)

    N-N, P-P, S-S are intentionally NOT checked: these appear in legitimate
    ligands (hydrazine, phosphine dimers, disulfides) and also produce false
    positives for metal complexes where heteroatom donors from different
    ligands end up close in ETKDG-generated geometries.

    Returns True (permissive) on any error or if required libraries are absent.
    """
    # Homodiatomic pairs to check (frozenset of atomic number, atomic number)
    _CHECKED_HOMODIATOMIC = {
        frozenset([8, 8]),    # O-O
        frozenset([9, 9]),    # F-F
        frozenset([17, 17]),  # Cl-Cl
        frozenset([35, 35]),  # Br-Br
        frozenset([53, 53]),  # I-I
    }

    if not OPENBABEL_AVAILABLE or not RDKIT_AVAILABLE:
        return True
    try:
        # Homodiatomic pairs present in the original SMILES (e.g. real peroxides)
        orig_mol = Chem.MolFromSmiles(original_smiles, sanitize=False)
        if orig_mol is None:
            return True
        try:
            orig_mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        orig_homo: set = set()
        for bond in orig_mol.GetBonds():
            n1 = bond.GetBeginAtom().GetAtomicNum()
            n2 = bond.GetEndAtom().GetAtomicNum()
            pair = frozenset([n1, n2])
            if pair in _CHECKED_HOMODIATOMIC:
                orig_homo.add(pair)

        # Homodiatomic pairs OB perceives in XYZ
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return True
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        rt_mol = pybel.readstring('xyz', std_xyz)
        try:
            from openbabel import openbabel as _ob
            for bond in _ob.OBMolBondIter(rt_mol.OBMol):
                n1 = bond.GetBeginAtom().GetAtomicNum()
                n2 = bond.GetEndAtom().GetAtomicNum()
                pair = frozenset([n1, n2])
                if pair in _CHECKED_HOMODIATOMIC and pair not in orig_homo:
                    logger.debug("Spurious homodiatomic bond in XYZ: %s-%s", n1, n2)
                    return False
        except Exception:
            return True

        return True
    except Exception:
        return True


def _has_atom_clash(mol, conf_id: int, min_dist: float = 0.5) -> bool:
    """Return True if any pair of non-bonded atoms is closer than *min_dist* Å.

    Checks only heavy atoms (non-H) for efficiency.  Bonded atom pairs
    are excluded since their distance is governed by bond length.
    """
    conf = mol.GetConformer(conf_id)
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    bonded = set()
    for bond in mol.GetBonds():
        bonded.add((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        bonded.add((bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()))
    for i in range(len(heavy)):
        pi = conf.GetAtomPosition(heavy[i])
        for j in range(i + 1, len(heavy)):
            if (heavy[i], heavy[j]) in bonded:
                continue
            pj = conf.GetAtomPosition(heavy[j])
            dx = pi.x - pj.x
            dy = pi.y - pj.y
            dz = pi.z - pj.z
            if dx*dx + dy*dy + dz*dz < min_dist * min_dist:
                return True
    return False


def _ml_distance_range(metal_symbol: str, donor_symbol: str) -> Tuple[float, float]:
    """Return (min_dist, max_dist) in Å for a metal-donor pair.

    Uses covalent radii sum with scaling factors:
    - min_dist = 0.7 × (r_M + r_D)
    - max_dist = 1.8 × (r_M + r_D)
    Falls back to (1.4, 3.5) if radii are unknown.
    """
    r_m = _COVALENT_RADII.get(metal_symbol)
    r_d = _COVALENT_RADII.get(donor_symbol)
    if r_m is not None and r_d is not None:
        r_sum = r_m + r_d
        return (0.7 * r_sum, 1.8 * r_sum)
    return (1.4, 3.5)


def _has_bad_geometry(mol, conf_id: int) -> bool:
    """Return True if the conformer has unrealistic metal-ligand geometry.

    Checks:
    1. Metal-ligand bond lengths must be within metal-specific distance range
       (derived from covalent radii, fallback to 1.4-3.5 Å)
    2. L-M-L angles: chelate bite angles (donors in the same chelate ring)
       may be as small as 40°; all other L-M-L pairs must be >=50°.
       This correctly allows 5- and 6-membered chelate ring bite angles
       (~55-70°) while still rejecting collapsed non-chelate geometries.
    """
    conf = mol.GetConformer(conf_id)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        metal_pos = conf.GetAtomPosition(metal_idx)
        neighbors = list(atom.GetNeighbors())
        if not neighbors:
            continue

        metal_sym = atom.GetSymbol()
        nbr_indices = [nb.GetIdx() for nb in neighbors]
        coord_positions = []
        for nbr in neighbors:
            nbr_pos = conf.GetAtomPosition(nbr.GetIdx())
            dx = nbr_pos.x - metal_pos.x
            dy = nbr_pos.y - metal_pos.y
            dz = nbr_pos.z - metal_pos.z
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            min_d, max_d = _ml_distance_range(metal_sym, nbr.GetSymbol())
            if dist < min_d or dist > max_d:
                return True
            coord_positions.append(nbr_pos)

        # Pre-compute chelate pairs (donors connected through non-metal path)
        chelate = _chelate_pairs(mol, metal_idx, nbr_indices)
        chelate_set = {frozenset(p) for p in chelate}
        n = len(coord_positions)
        n_pairs = n * (n - 1) // 2

        all_angles = []
        for i in range(n):
            for j in range(i + 1, n):
                pa, pb = coord_positions[i], coord_positions[j]
                v1 = (pa.x - metal_pos.x, pa.y - metal_pos.y, pa.z - metal_pos.z)
                v2 = (pb.x - metal_pos.x, pb.y - metal_pos.y, pb.z - metal_pos.z)
                dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
                mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
                mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                if mag1 < 1e-8 or mag2 < 1e-8:
                    return True
                cos_a = max(-1.0, min(1.0, dot / (mag1 * mag2)))
                angle = math.degrees(math.acos(cos_a))
                all_angles.append((angle, frozenset([nbr_indices[i], nbr_indices[j]])))
                is_chelate_pair = frozenset([nbr_indices[i], nbr_indices[j]]) in chelate_set
                min_angle = 40.0 if is_chelate_pair else 50.0
                if angle < min_angle:
                    return True

        # NOTE: The old macrocyclic tetrahedral-rejection check (CN=4, all
        # pairs chelate, max angle < 115°) was removed.  It incorrectly
        # rejected valid tetrahedral structures (ideal 109.5°).  With the
        # BFS path-length cutoff in _chelate_pairs (max_path=4), macrocyclic
        # opposite-donor pairs are no longer marked as chelate, so the
        # "all pairs chelate" condition rarely triggers anyway.  The
        # _geometry_quality_score function already handles tetrahedral vs
        # square-planar ranking correctly.
    return False


def _geometry_quality_score(mol, conf_id: int) -> float:
    """Score how regular the metal coordination geometry is (lower = better).

    Measures:
    - Spread of M-L bond lengths (std-dev, ideally 0 for identical ligands)
    - Angular regularity:
      - 4-coordinate centers: best of tetrahedral (109.5 deg) OR
        square-planar (90/180 deg) targets
      - Other coordinations: deviation from nearest ideal (90/180 deg)
    """
    conf = mol.GetConformer(conf_id)
    total_penalty = 0.0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_pos = conf.GetAtomPosition(atom.GetIdx())
        neighbors = list(atom.GetNeighbors())
        if not neighbors:
            continue

        # Bond length uniformity
        dists = []
        coord_positions = []
        for nbr in neighbors:
            nbr_pos = conf.GetAtomPosition(nbr.GetIdx())
            dx = nbr_pos.x - metal_pos.x
            dy = nbr_pos.y - metal_pos.y
            dz = nbr_pos.z - metal_pos.z
            dists.append(math.sqrt(dx*dx + dy*dy + dz*dz))
            coord_positions.append(nbr_pos)

        if dists:
            mean_d = sum(dists) / len(dists)
            std_d = math.sqrt(sum((d - mean_d)**2 for d in dists) / len(dists))
            total_penalty += std_d * 10  # weight bond-length spread

        # Angle regularity. For 4-coordinate centers, allow both tetrahedral
        # and square-planar patterns and keep the better one.
        angles = []
        for i in range(len(coord_positions)):
            for j in range(i + 1, len(coord_positions)):
                pa, pb = coord_positions[i], coord_positions[j]
                v1 = (pa.x - metal_pos.x, pa.y - metal_pos.y, pa.z - metal_pos.z)
                v2 = (pb.x - metal_pos.x, pb.y - metal_pos.y, pb.z - metal_pos.z)
                dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
                mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
                mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                if mag1 < 1e-8 or mag2 < 1e-8:
                    continue
                cos_a = max(-1.0, min(1.0, dot / (mag1 * mag2)))
                angles.append(math.degrees(math.acos(cos_a)))

        if not angles:
            continue

        if len(neighbors) == 2 and len(angles) == 1:
            # CN=2 (linear): ideal angle = 180°
            total_penalty += abs(angles[0] - 180.0)
        elif len(neighbors) == 3 and len(angles) == 3:
            # CN=3: best of trigonal-planar (all 120°) vs T-shaped (90°,90°,180°)
            tp_pen = sum(abs(a - 120.0) for a in angles)
            ts_targets = sorted([90.0, 90.0, 180.0])
            ts_pen = sum(
                abs(a - t) for a, t in zip(sorted(angles), ts_targets)
            )
            total_penalty += min(tp_pen, ts_pen)
        elif len(neighbors) == 4 and len(angles) == 6:
            tetra_pen = sum(abs(a - 109.5) for a in angles)
            square_targets = [90.0, 90.0, 90.0, 90.0, 180.0, 180.0]
            square_pen = sum(
                abs(a - t) for a, t in zip(sorted(angles), square_targets)
            )
            total_penalty += min(tetra_pen, square_pen)
        elif len(neighbors) == 5:
            # CN=5: score against best of TBP or SP ideal angles.
            # TBP (D3h): 1×180° (ax-ax), 6×90° (ax-eq), 3×120° (eq-eq)
            tbp_targets = sorted([90, 90, 90, 90, 90, 90, 120, 120, 120, 180])
            # SP (C4v): 4×90° (cis-basal), 2×180° (trans-basal), 4×~100° (apical-basal)
            sp_targets = sorted([90, 90, 90, 90, 100, 100, 100, 100, 180, 180])
            sorted_a = sorted(angles)
            tbp_pen = sum(abs(a - t) for a, t in zip(sorted_a, tbp_targets))
            sp_pen = sum(abs(a - t) for a, t in zip(sorted_a, sp_targets))
            total_penalty += min(tbp_pen, sp_pen)
        elif len(neighbors) == 7:
            # CN=7 (PBP, D5h): penalize each angle against nearest of 72°/90°/144°/180°
            ideal_7 = [72.0, 90.0, 144.0, 180.0]
            for a in angles:
                total_penalty += min(abs(a - t) for t in ideal_7)
        elif len(neighbors) == 8:
            # CN=8: score against best of SAP (D4d) or DD (D2d)
            # SAP: 8×~73° + 8×~118° + 4×~143° + 4×~70° + 4×~180° (28 pairs)
            # Simplified: penalize each angle against nearest ideal
            ideal_sap = [52.4, 73.1, 118.5, 143.1, 180.0]
            ideal_dd = [62.2, 73.7, 117.4, 143.6, 180.0]
            sap_pen = sum(min(abs(a - t) for t in ideal_sap) for a in angles)
            dd_pen = sum(min(abs(a - t) for t in ideal_dd) for a in angles)
            total_penalty += min(sap_pen, dd_pen)
        else:
            # General: penalize distance from nearest ideal (90 or 180)
            for a in angles:
                dev_90 = abs(a - 90)
                dev_180 = abs(a - 180)
                total_penalty += min(dev_90, dev_180)

    # Chelate-ring planarity bonus: aromatic rings that coordinate to a metal
    # should be flat.  Penalize RMSD of ring atoms from the best-fit plane.
    try:
        Chem.FastFindRings(mol)
        ring_info = mol.GetRingInfo()
        metal_atom_indices = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET}
        if metal_atom_indices and ring_info and ring_info.NumRings() > 0:
            for ring in ring_info.AtomRings():
                # Check whether ≥2 ring atoms are direct neighbours of a metal
                ring_set = set(ring)
                n_coord_in_ring = sum(
                    1 for ridx in ring_set
                    if any(nbr.GetIdx() in metal_atom_indices
                           for nbr in mol.GetAtomWithIdx(ridx).GetNeighbors())
                )
                if n_coord_in_ring < 2:
                    continue  # Not a chelate ring

                # Collect ring-atom 3D positions
                positions = []
                for ridx in ring:
                    pos = conf.GetAtomPosition(ridx)
                    positions.append((pos.x, pos.y, pos.z))

                if len(positions) < 3:
                    continue

                # Fit a plane via centroid + SVD-like normal estimation using
                # cross products of consecutive edge vectors (rotation-invariant).
                cx = sum(p[0] for p in positions) / len(positions)
                cy = sum(p[1] for p in positions) / len(positions)
                cz = sum(p[2] for p in positions) / len(positions)
                vecs = [(p[0]-cx, p[1]-cy, p[2]-cz) for p in positions]

                # Accumulate a rough normal via cross products of consecutive vectors
                nx, ny, nz = 0.0, 0.0, 0.0
                n_v = len(vecs)
                for vi in range(n_v):
                    a_v = vecs[vi]
                    b_v = vecs[(vi + 1) % n_v]
                    nx += a_v[1]*b_v[2] - a_v[2]*b_v[1]
                    ny += a_v[2]*b_v[0] - a_v[0]*b_v[2]
                    nz += a_v[0]*b_v[1] - a_v[1]*b_v[0]
                n_mag = math.sqrt(nx*nx + ny*ny + nz*nz)
                if n_mag < 1e-8:
                    continue
                nx /= n_mag
                ny /= n_mag
                nz /= n_mag

                # RMSD of atoms from the plane
                deviations = [abs(v[0]*nx + v[1]*ny + v[2]*nz) for v in vecs]
                planarity_rmsd = math.sqrt(
                    sum(d*d for d in deviations) / len(deviations)
                )
                total_penalty += planarity_rmsd * 5.0
    except Exception:
        pass  # Ring-planarity scoring is optional; never block normal scoring

    return total_penalty


def _segment_distance_sq(p1, p2, p3, p4) -> float:
    """Squared minimum distance between 3D line segments P1-P2 and P3-P4.

    Each point is a tuple ``(x, y, z)``.  Uses the parametric approach:
    closest points on segment ``P1 + s*(P2-P1)`` and ``P3 + t*(P4-P3)``
    with ``s, t`` clamped to ``[0, 1]``.
    """
    d1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    d2 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])
    r = (p1[0] - p3[0], p1[1] - p3[1], p1[2] - p3[2])

    a = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]
    e = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]
    f = d2[0] * r[0] + d2[1] * r[1] + d2[2] * r[2]

    EPS = 1e-12
    if a <= EPS and e <= EPS:
        return r[0] ** 2 + r[1] ** 2 + r[2] ** 2

    if a <= EPS:
        t = max(0.0, min(1.0, f / e))
        v = (r[0] - t * d2[0], r[1] - t * d2[1], r[2] - t * d2[2])
        return v[0] ** 2 + v[1] ** 2 + v[2] ** 2

    c = d1[0] * r[0] + d1[1] * r[1] + d1[2] * r[2]
    if e <= EPS:
        s = max(0.0, min(1.0, -c / a))
        v = (r[0] + s * d1[0], r[1] + s * d1[1], r[2] + s * d1[2])
        return v[0] ** 2 + v[1] ** 2 + v[2] ** 2

    b = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2]
    denom = a * e - b * b

    if abs(denom) > EPS:
        s = max(0.0, min(1.0, (b * f - c * e) / denom))
    else:
        s = 0.0

    t = (b * s + f) / e
    if t < 0.0:
        t = 0.0
        s = max(0.0, min(1.0, -c / a))
    elif t > 1.0:
        t = 1.0
        s = max(0.0, min(1.0, (b - c) / a))

    v = (r[0] + s * d1[0] - t * d2[0],
         r[1] + s * d1[1] - t * d2[1],
         r[2] + s * d1[2] - t * d2[2])
    return v[0] ** 2 + v[1] ** 2 + v[2] ** 2


def _has_ligand_intertwining(mol, conf_id: int, threshold: float = 0.3) -> bool:
    """Return True if ligands are physically intertwined in the conformer.

    Decomposes the molecule into ligand fragments by removing metal centers,
    then checks if any heavy-atom bond from one fragment comes closer than
    *threshold* Å to any heavy-atom bond from another fragment (segment-to-
    segment distance).  This catches unphysical conformers where ligand arms
    pass through each other.
    """
    conf = mol.GetConformer(conf_id)
    metal_idxs = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET}
    if not metal_idxs:
        return False

    # Build adjacency graph for non-metal atoms (bonds to metals removed)
    non_metal = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() not in _METAL_SET}
    adj: Dict[int, set] = {i: set() for i in non_metal}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in metal_idxs or j in metal_idxs:
            continue
        adj[i].add(j)
        adj[j].add(i)

    # Connected components = ligand fragments
    visited: set = set()
    fragments: List[set] = []
    for start in sorted(non_metal):
        if start in visited:
            continue
        frag: set = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            frag.add(node)
            for nbr in adj.get(node, ()):
                if nbr not in visited:
                    stack.append(nbr)
        fragments.append(frag)

    if len(fragments) < 2:
        return False

    # Map atom index -> fragment index
    atom_frag: Dict[int, int] = {}
    for fi, frag in enumerate(fragments):
        for aidx in frag:
            atom_frag[aidx] = fi

    # Collect heavy-atom bonds per fragment (skip bonds involving H)
    frag_bonds: Dict[int, List[Tuple[int, int]]] = {}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in metal_idxs or j in metal_idxs:
            continue
        if mol.GetAtomWithIdx(i).GetAtomicNum() <= 1:
            continue
        if mol.GetAtomWithIdx(j).GetAtomicNum() <= 1:
            continue
        fi = atom_frag.get(i)
        if fi is None:
            continue
        frag_bonds.setdefault(fi, []).append((i, j))

    # Check inter-fragment bond-segment distances
    thresh_sq = threshold * threshold
    fkeys = sorted(frag_bonds)
    for a_idx in range(len(fkeys)):
        for b_idx in range(a_idx + 1, len(fkeys)):
            for i1, j1 in frag_bonds[fkeys[a_idx]]:
                p1 = conf.GetAtomPosition(i1)
                p2 = conf.GetAtomPosition(j1)
                pa, pb = (p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z)
                for i2, j2 in frag_bonds[fkeys[b_idx]]:
                    p3 = conf.GetAtomPosition(i2)
                    p4 = conf.GetAtomPosition(j2)
                    pc, pd = (p3.x, p3.y, p3.z), (p4.x, p4.y, p4.z)
                    if _segment_distance_sq(pa, pb, pc, pd) < thresh_sq:
                        return True
    return False


def _conformer_rmsd(mol, cid_a: int, cid_b: int) -> float:
    """Heavy-atom RMSD between two conformers after Kabsch alignment.

    Uses RDKit's ``AlignMol`` on a temporary copy so the original molecule
    coordinates are not modified.  Falls back to a sorted distance-matrix
    comparison (rotation-invariant) if alignment fails.
    """
    try:
        from rdkit.Chem import rdMolAlign
        heavy_map = [(i, i) for i, a in enumerate(mol.GetAtoms())
                     if a.GetAtomicNum() > 1]
        if not heavy_map:
            return float('inf')
        # AlignMol modifies probe coordinates → work on a copy
        probe = Chem.RWMol(mol)
        return rdMolAlign.AlignMol(probe, mol, prbCid=cid_a, refCid=cid_b,
                                   atomMap=heavy_map)
    except Exception:
        pass

    # Fallback: sorted distance-matrix comparison (rotation-invariant)
    try:
        conf_a = mol.GetConformer(cid_a)
        conf_b = mol.GetConformer(cid_b)
        heavy = [i for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() > 1]
        if len(heavy) < 2:
            return float('inf')
        dists_a: List[float] = []
        dists_b: List[float] = []
        for i in range(len(heavy)):
            pa = conf_a.GetAtomPosition(heavy[i])
            pb = conf_b.GetAtomPosition(heavy[i])
            for j in range(i + 1, len(heavy)):
                qa = conf_a.GetAtomPosition(heavy[j])
                qb = conf_b.GetAtomPosition(heavy[j])
                dists_a.append(math.sqrt(
                    (pa.x - qa.x) ** 2 + (pa.y - qa.y) ** 2 + (pa.z - qa.z) ** 2))
                dists_b.append(math.sqrt(
                    (pb.x - qb.x) ** 2 + (pb.y - qb.y) ** 2 + (pb.z - qb.z) ** 2))
        dists_a.sort()
        dists_b.sort()
        return math.sqrt(
            sum((a - b) ** 2 for a, b in zip(dists_a, dists_b)) / len(dists_a))
    except Exception:
        return float('inf')


def _angle_class(pos_metal, pos_a, pos_b) -> str:
    """Classify the angle A-Metal-B as 'cis' (<135 deg) or 'trans'.

    Uses 135 deg as threshold (midpoint between ideal octahedral 90 and
    180 deg) for robust classification despite geometric noise from
    distance-geometry embedding.
    """
    v1 = (pos_a.x - pos_metal.x, pos_a.y - pos_metal.y, pos_a.z - pos_metal.z)
    v2 = (pos_b.x - pos_metal.x, pos_b.y - pos_metal.y, pos_b.z - pos_metal.z)
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    if mag1 < 1e-8 or mag2 < 1e-8:
        return 'cis'
    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    angle_deg = math.degrees(math.acos(cos_angle))
    return 'cis' if angle_deg < 135 else 'trans'


def _donor_type_map(mol) -> Dict[int, tuple]:
    """Compute a donor type for each atom based on its chemical environment.

    Uses Morgan fingerprint bits at radius 2 to distinguish chemically
    inequivalent atoms of the same element (e.g. N atoms in a tridentate
    terpyridine vs. a bidentate phenylpyridine).

    Returns:
        Dict mapping atom index to ``(element_symbol, env_hash)`` where
        ``env_hash`` is a frozenset of Morgan bits unique to that
        chemical environment.
    """
    result: Dict[int, tuple] = {}
    if not RDKIT_AVAILABLE:
        return result

    # Compute Morgan fingerprint with bit info.
    # Ensure ring info is initialized (required by Morgan FP; may be missing
    # on partially-sanitized metal-complex mols).
    atom_bits: Dict[int, set] = {}
    try:
        Chem.FastFindRings(mol)
    except Exception:
        pass
    try:
        from rdkit.Chem import rdFingerprintGenerator
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
        ao = rdFingerprintGenerator.AdditionalOutput()
        ao.AllocateBitInfoMap()
        gen.GetFingerprint(mol, additionalOutput=ao)
        for bit_id, entries in ao.GetBitInfoMap().items():
            for atom_idx, _radius in entries:
                if atom_idx not in atom_bits:
                    atom_bits[atom_idx] = set()
                atom_bits[atom_idx].add(bit_id)
    except Exception:
        # Fallback to legacy API
        try:
            from rdkit.Chem import rdMolDescriptors
            fp_info: Dict[int, list] = {}
            rdMolDescriptors.GetMorganFingerprint(mol, 2, bitInfo=fp_info)
            for bit_id, entries in fp_info.items():
                for atom_idx, _radius in entries:
                    if atom_idx not in atom_bits:
                        atom_bits[atom_idx] = set()
                    atom_bits[atom_idx].add(bit_id)
        except Exception:
            pass

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        env_hash = frozenset(atom_bits.get(idx, set()))
        result[idx] = (atom.GetSymbol(), env_hash)

    return result


def _compute_coordination_fingerprint(mol, conf_id: int, dtype_map: Optional[Dict[int, tuple]] = None) -> tuple:
    """Compute a hashable fingerprint describing the coordination geometry.

    Hybrid approach for robustness:
    1. **Trans pairs** (top floor(N/2) angles): which element pairs sit
       across the metal.  Robust against noise for mixed-element donors.
    2. **Same-element cis/trans pattern**: for each pair of identical-
       element donors, classify as cis (<135 deg) or trans.  This adds
       sensitivity for all-same-element complexes (e.g. Fe with 6 N)
       where trans-pair elements alone cannot distinguish isomers.
    3. **Detailed trans pairs**: uses Morgan-invariant-enriched donor types
       to distinguish same-element donors in different chemical environments
       (e.g. N in terpyridine vs. N in phenylpyridine of an MABC complex).

    Args:
        mol: RDKit molecule with conformers
        conf_id: Conformer ID to analyze
        dtype_map: Pre-computed donor type map from ``_donor_type_map(mol)``.
            If None, computed on the fly.
    """
    conf = mol.GetConformer(conf_id)
    fp_parts: List[tuple] = []

    # Use pre-computed donor types or compute on the fly
    if dtype_map is None:
        dtype_map = _donor_type_map(mol)

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue

        metal_pos = conf.GetAtomPosition(atom.GetIdx())

        coord_atoms = sorted(
            ((nbr.GetSymbol(), nbr.GetIdx()) for nbr in atom.GetNeighbors()),
            key=lambda x: (x[0], x[1]),
        )
        n_coord = len(coord_atoms)
        if n_coord < 2:
            continue

        # Compute all pairwise angles
        # entries: (angle, symA, symB, idxA, idxB)
        angle_pairs: List[Tuple[float, str, str, int, int]] = []
        for i in range(n_coord):
            for j in range(i + 1, n_coord):
                sym_a, idx_a = coord_atoms[i]
                sym_b, idx_b = coord_atoms[j]
                v1 = (conf.GetAtomPosition(idx_a).x - metal_pos.x,
                      conf.GetAtomPosition(idx_a).y - metal_pos.y,
                      conf.GetAtomPosition(idx_a).z - metal_pos.z)
                v2 = (conf.GetAtomPosition(idx_b).x - metal_pos.x,
                      conf.GetAtomPosition(idx_b).y - metal_pos.y,
                      conf.GetAtomPosition(idx_b).z - metal_pos.z)
                dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
                mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
                mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                if mag1 < 1e-8 or mag2 < 1e-8:
                    continue
                cos_a = max(-1.0, min(1.0, dot / (mag1 * mag2)))
                angle = math.degrees(math.acos(cos_a))
                angle_pairs.append((angle, sym_a, sym_b, idx_a, idx_b))

        # Part 1: trans pair elements via disjoint pairing that maximizes
        # total trans-angle sum (more robust than taking top N/2 raw angles).
        n_trans = n_coord // 2
        trans_pairs_raw: List[Tuple[str, str]] = []
        # Also track detailed trans pairs using donor types
        detailed_trans_raw: List[Tuple[tuple, tuple]] = []
        if n_coord % 2 == 0 and n_coord <= 8 and angle_pairs:
            idx_by_symbol = {idx: sym for sym, idx in coord_atoms}
            angle_map: Dict[Tuple[int, int], float] = {}
            for angle, _sa, _sb, ia, ib in angle_pairs:
                key = (ia, ib) if ia < ib else (ib, ia)
                angle_map[key] = angle

            donor_indices = tuple(sorted(idx_by_symbol.keys()))
            memo: Dict[Tuple[int, ...], Tuple[float, List[Tuple[int, int]]]] = {}

            def _best_pairing(rem: Tuple[int, ...]) -> Tuple[float, List[Tuple[int, int]]]:
                if len(rem) < 2:
                    return 0.0, []
                if rem in memo:
                    return memo[rem]

                first = rem[0]
                best_score = -1.0
                best_pairs: List[Tuple[int, int]] = []
                for k in range(1, len(rem)):
                    second = rem[k]
                    key = (first, second) if first < second else (second, first)
                    angle_val = angle_map.get(key, 0.0)
                    rest = rem[1:k] + rem[k + 1:]
                    sub_score, sub_pairs = _best_pairing(rest)
                    score = angle_val + sub_score
                    if score > best_score:
                        best_score = score
                        best_pairs = [(first, second)] + sub_pairs

                memo[rem] = (best_score, best_pairs)
                return memo[rem]

            _score, pair_idx = _best_pairing(donor_indices)
            for ia, ib in pair_idx[:n_trans]:
                trans_pairs_raw.append((idx_by_symbol[ia], idx_by_symbol[ib]))
                detailed_trans_raw.append((dtype_map.get(ia, (idx_by_symbol[ia],)),
                                           dtype_map.get(ib, (idx_by_symbol[ib],))))
        else:
            angle_pairs.sort(key=lambda x: -x[0])  # descending
            for _angle, sa, sb, ia, ib in angle_pairs[:n_trans]:
                trans_pairs_raw.append((sa, sb))
                detailed_trans_raw.append((dtype_map.get(ia, (sa,)),
                                           dtype_map.get(ib, (sb,))))

        trans_pairs = tuple(sorted(
            tuple(sorted((sa, sb))) for sa, sb in trans_pairs_raw
        ))
        detailed_trans = tuple(sorted(
            tuple(sorted((da, db))) for da, db in detailed_trans_raw
        ))

        # Part 2: pairwise cis/trans pattern for ALL same-element donor
        # pairs (not just when all donors share one element).  This adds
        # sensitivity for mixed complexes like MA2B2C2 where two different
        # "cis-A" arrangements produce identical trans-pair signatures but
        # differ in which same-element pairs are cis vs trans.
        same_elem_pattern: List[tuple] = []
        for angle, sa, sb, _ia, _ib in angle_pairs:
            if sa == sb:  # same element pair
                cls = 'cis' if angle < 135 else 'trans'
                same_elem_pattern.append((sa, cls))
        same_elem_pattern.sort()
        fp_parts.append((trans_pairs, tuple(same_elem_pattern), detailed_trans))

    return tuple(sorted(fp_parts))


def _classify_isomer_label(fingerprint: tuple, mol) -> str:
    """Translate a coordination fingerprint to a human-readable label.

    The fingerprint is ``((trans_pairs, same_elem_pattern, detailed_trans), ...)``
    per metal center.  Uses trans-pair analysis to classify all common
    octahedral and 4-coordinate patterns:

    6-coordinate (octahedral):
    - **MA6** / **MA5B**: single isomer → ``''``
    - **MA4B2**: cis / trans (minority pair B trans or not)
    - **MA3B3**: fac / mer
    - **MA2B2C2**: all-cis / all-trans / X-trans (X = element
      whose pair is trans while the other two are cis)
    - **MA3B2C**: fac/mer for A + cis/trans for B

    4-coordinate:
    - **MA2B2**: cis / trans
    """
    # Count coordinating atoms per element from the mol
    n_coord = 0
    elem_counts: Dict[str, int] = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        for nbr in atom.GetNeighbors():
            sym = nbr.GetSymbol()
            elem_counts[sym] = elem_counts.get(sym, 0) + 1
            n_coord += 1

    count_signature = sorted(elem_counts.values())

    for metal_fp in fingerprint:
        trans_pairs, same_elem_pattern = metal_fp[0], metal_fp[1]

        # Count same-element trans pairs from trans_pairs
        same_trans: Dict[str, int] = {}
        for pair in trans_pairs:
            if pair[0] == pair[1]:
                same_trans[pair[0]] = same_trans.get(pair[0], 0) + 1

        # Supplement from same_elem_pattern ONLY when trans_pairs contain
        # exclusively same-element pairs (= all donors share one element).
        # In this case trans_pairs cannot distinguish isomers by element
        # and the cis/trans angle classification is the only differentiator.
        # For mixed-element complexes, trans_pairs already carry the needed
        # information; adding same_elem_pattern would double-count.
        has_hetero_trans = any(p[0] != p[1] for p in trans_pairs)
        if not has_hetero_trans:
            for sym, cls in same_elem_pattern:
                if cls == 'trans':
                    if sym not in same_trans:
                        same_trans[sym] = 0
                    same_trans[sym] += 1

        # --- 6-coordinate patterns ---
        if n_coord == 6:
            # MA6 or MA5B with homogeneous element set: usually a single isomer.
            # Exception: when all donors share one element symbol but have two
            # distinct Morgan-hash types (3+3), we can still classify fac/mer
            # using the detailed_trans information.  This handles tris-bidentate
            # complexes like Fe(citrate)3 where e.g. all donors are O but split
            # into alkoxo-O and carboxylato-O.
            if len(elem_counts) == 1 or count_signature == [1, 5]:
                if count_signature != [1, 5] and len(metal_fp) > 2:
                    detailed = metal_fp[2]  # tuple of (type_i, type_j) trans pairs
                    # Collect all donor types present
                    all_types = set()
                    for pair in detailed:
                        all_types.add(pair[0])
                        all_types.add(pair[1])
                    if len(all_types) == 2:
                        # Two distinct Morgan types, 3+3 split: fac or mer
                        n_same_type_trans = sum(1 for p in detailed if p[0] == p[1])
                        if n_same_type_trans == 0:
                            return 'fac'
                        return 'mer'
                return ''

            # MA4B2: cis/trans based on minority element (count==2)
            if count_signature == [2, 4]:
                minority = [s for s, c in elem_counts.items() if c == 2][0]
                if same_trans.get(minority, 0) >= 1:
                    return 'trans'
                return 'cis'

            # MA3B3: fac/mer
            if count_signature == [3, 3]:
                for sym, count in elem_counts.items():
                    if count == 3:
                        if same_trans.get(sym, 0) == 0:
                            return 'fac'
                        return 'mer'

            # MA2B2C2: all-cis / all-trans / X-trans
            if count_signature == [2, 2, 2]:
                elems_with_2 = [s for s, c in elem_counts.items() if c == 2]
                trans_elems = [s for s in elems_with_2 if same_trans.get(s, 0) >= 1]
                n_trans = len(trans_elems)
                if n_trans == 0:
                    return 'all-cis'
                if n_trans == 3:
                    return 'all-trans'
                if n_trans == 1:
                    return f'{trans_elems[0]}-trans'
                # 2 trans pairs: label by the one that is cis
                cis_elems = [s for s in elems_with_2 if s not in trans_elems]
                if len(cis_elems) == 1:
                    return f'{cis_elems[0]}-cis'
                return ''

            # MA3B2C: fac/mer for A (count==3) + cis/trans for B (count==2)
            if count_signature == [1, 2, 3]:
                a_sym = [s for s, c in elem_counts.items() if c == 3][0]
                b_sym = [s for s, c in elem_counts.items() if c == 2][0]
                a_label = 'mer' if same_trans.get(a_sym, 0) >= 1 else 'fac'
                b_label = 'trans' if same_trans.get(b_sym, 0) >= 1 else 'cis'
                return f'{a_label}-{b_label}'

            # MA4BC: the majority element has 4, the other two have 1 each.
            # Classify by which pairs are trans.
            if count_signature == [1, 1, 4]:
                majority = [s for s, c in elem_counts.items() if c == 4][0]
                minorities = sorted(s for s, c in elem_counts.items() if c == 1)
                # Check if the two minorities are trans to each other
                minority_trans = any(
                    (sorted(p) == sorted(minorities)) for p in trans_pairs
                )
                n_maj_trans = same_trans.get(majority, 0)
                if minority_trans:
                    return f'{"/".join(minorities)}-trans'
                if n_maj_trans >= 2:
                    return f'{majority}-all-trans'
                return f'{"/".join(minorities)}-cis'

        # --- 2-coordinate patterns ---
        if n_coord == 2:
            return 'linear'

        # --- 3-coordinate patterns ---
        if n_coord == 3:
            if len(elem_counts) == 1:
                return ''  # MA3: single isomer (trigonal-planar)
            # MA2B: cis/trans only for T-shaped geometry (where trans pair exists)
            if count_signature == [1, 2] and trans_pairs:
                minority = [s for s, c in elem_counts.items() if c == 1][0]
                is_trans = any(minority in p for p in trans_pairs)
                if is_trans:
                    return 'trans'
                return 'cis'
            return ''

        # --- 4-coordinate patterns ---
        if n_coord == 4:
            if count_signature == [2, 2]:
                minority = sorted(elem_counts.keys())[0]
                if same_trans.get(minority, 0) >= 1:
                    return 'trans'
                return 'cis'

        # --- 5-coordinate patterns ---
        if n_coord == 5:
            if len(elem_counts) == 1:
                return ''  # MA5: single isomer
            if count_signature == [1, 4]:
                minority = [s for s, c in elem_counts.items() if c == 1][0]
                # Axial donors appear in a trans-pair; equatorial don't
                is_axial = any(minority in p for p in trans_pairs)
                return 'axial' if is_axial else 'equatorial'
            if count_signature == [2, 3]:
                minority = [s for s, c in elem_counts.items() if c == 2][0]
                n_trans = same_trans.get(minority, 0)
                return 'diaxial' if n_trans >= 1 else 'ax-eq'

        # --- 7-coordinate patterns ---
        if n_coord == 7:
            if len(elem_counts) == 1:
                return ''
            minority = min(elem_counts, key=elem_counts.get)
            c = elem_counts[minority]
            if c == 1:
                return 'axial' if any(minority in p for p in trans_pairs) else 'equatorial'
            if c == 2:
                n_trans = same_trans.get(minority, 0)
                return 'ax-ax' if n_trans >= 1 else 'ax-eq'

        # --- 8-coordinate patterns ---
        if n_coord == 8:
            if len(elem_counts) == 1:
                return ''  # MA8: single isomer
            # For mixed-element CN=8: count trans-pair types for differentiation
            n_total_same_trans = sum(same_trans.values())
            if count_signature == [4, 4]:
                # MA4B4: classify by number of same-element trans pairs
                if n_total_same_trans == 0:
                    return 'all-cis'
                if n_total_same_trans >= 4:
                    return 'all-trans'
                return f'{n_total_same_trans}-trans'
            return ''

    return ''


# ---------------------------------------------------------------------------
# Topological isomer enumerator (Feature 1)
# ---------------------------------------------------------------------------

def _chelate_pairs(mol, metal_idx: int, donor_indices: List[int],
                    max_path: int = 4) -> List[FrozenSet]:
    """Find donor pairs connected through a short non-metal path (chelate constraints).

    BFS from each donor atom to every other donor, blocking the metal.
    A successful path with length <= *max_path* bonds means the two donors
    form a chelate ring small enough to enforce a *cis* constraint (up to a
    ``max_path + 1``-membered chelate ring counting the metal).

    Longer paths (e.g. opposite donors of a 14-membered macrocycle) are
    **not** marked as chelate because they can adopt trans arrangements.

    Args:
        max_path: Maximum number of bonds in the non-metal path between two
            donors to count as a chelate pair.  Default 4 corresponds to a
            6-membered chelate ring (donor–4 bridge atoms–donor + metal).

    Returns a list of frozensets {donor_i_idx, donor_j_idx}.
    """
    pairs: List[FrozenSet] = []
    n = len(donor_indices)
    for i in range(n):
        for j in range(i + 1, n):
            start = donor_indices[i]
            target = donor_indices[j]
            # BFS with distance tracking, blocking the metal atom
            visited = {metal_idx, start}
            queue: List[Tuple[int, int]] = [(start, 0)]  # (atom_idx, distance)
            found_dist = -1
            while queue and found_dist < 0:
                current, dist = queue.pop(0)
                if dist >= max_path:
                    # Cannot reach target within max_path from here
                    continue
                for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni == target:
                        found_dist = dist + 1
                        break
                    if ni not in visited:
                        visited.add(ni)
                        queue.append((ni, dist + 1))
            if 0 < found_dist <= max_path:
                pairs.append(frozenset([donor_indices[i], donor_indices[j]]))
    return pairs


def _canonical_oh(types: tuple) -> tuple:
    """Canonical form for octahedral (3 trans pairs: 0-1, 2-3, 4-5)."""
    pairs = tuple(sorted([
        tuple(sorted([types[0], types[1]])),
        tuple(sorted([types[2], types[3]])),
        tuple(sorted([types[4], types[5]])),
    ]))
    return ('OH', pairs)


def _canonical_sq(types: tuple) -> tuple:
    """Canonical form for square-planar (2 trans pairs: 0-1, 2-3)."""
    pairs = tuple(sorted([
        tuple(sorted([types[0], types[1]])),
        tuple(sorted([types[2], types[3]])),
    ]))
    return ('SQ', pairs)


def _canonical_tbp(types: tuple) -> tuple:
    """Canonical form for trigonal-bipyramidal (axial: 0,1; equatorial: 2,3,4)."""
    axial = tuple(sorted([types[0], types[1]]))
    equatorial = tuple(sorted([types[2], types[3], types[4]]))
    return ('TBP', axial, equatorial)


def _canonical_sp(types: tuple) -> tuple:
    """Canonical form for square-pyramidal (apical: 0; basal: 1,2,3,4; trans pairs 1-3, 2-4)."""
    basal_pairs = tuple(sorted([
        tuple(sorted([types[1], types[3]])),
        tuple(sorted([types[2], types[4]])),
    ]))
    return ('SP', types[0], basal_pairs)


def _canonical_pbp(types: tuple) -> tuple:
    """Canonical form for pentagonal-bipyramidal (axial: 0,1; equatorial: 2-6)."""
    axial = tuple(sorted([types[0], types[1]]))
    equatorial = tuple(sorted([types[2], types[3], types[4], types[5], types[6]]))
    return ('PBP', axial, equatorial)


def _canonical_th(types: tuple) -> tuple:
    """Canonical form for tetrahedral (all 4 sites equivalent, no trans pairs)."""
    return ('TH', tuple(sorted(types)))


def _canonical_lin(types: tuple) -> tuple:
    """Canonical form for linear (2 sites, 1 trans pair: 0-1)."""
    return ('LIN', tuple(sorted(types)))


def _canonical_tp(types: tuple) -> tuple:
    """Canonical form for trigonal-planar (3 equivalent sites, no trans pairs)."""
    return ('TP', tuple(sorted(types)))


def _canonical_ts(types: tuple) -> tuple:
    """Canonical form for T-shaped (trans pair: 0-1; unique: 2)."""
    trans = tuple(sorted([types[0], types[1]]))
    return ('TS', trans, types[2])


def _canonical_sap(types: tuple) -> tuple:
    """Canonical form for square-antiprismatic (D4d).

    Positions 0-3 = top square, 4-7 = bottom square (staggered by 45°).
    Trans pairs: 0-6, 1-7, 2-4, 3-5 (opposite vertices through center).
    """
    trans_pairs = tuple(sorted([
        tuple(sorted([types[0], types[6]])),
        tuple(sorted([types[1], types[7]])),
        tuple(sorted([types[2], types[4]])),
        tuple(sorted([types[3], types[5]])),
    ]))
    return ('SAP', trans_pairs)


def _canonical_dd(types: tuple) -> tuple:
    """Canonical form for dodecahedral (D2d, = two interpenetrating trapezoids).

    Positions 0-3 = A sites (larger trapezoid), 4-7 = B sites (smaller).
    Trans pairs: 0-2, 1-3 (A-site), 4-6, 5-7 (B-site).
    """
    a_pairs = tuple(sorted([
        tuple(sorted([types[0], types[2]])),
        tuple(sorted([types[1], types[3]])),
    ]))
    b_pairs = tuple(sorted([
        tuple(sorted([types[4], types[6]])),
        tuple(sorted([types[5], types[7]])),
    ]))
    return ('DD', a_pairs, b_pairs)


# Trans position pairs (indices into the geometry vector list) for each geometry
_TOPO_TRANS_POSITIONS: Dict[str, List[Tuple[int, int]]] = {
    'LIN': [(0, 1)],
    'TP':  [],
    'TS':  [(0, 1)],
    'OH':  [(0, 1), (2, 3), (4, 5)],
    'SQ':  [(0, 1), (2, 3)],
    'TH':  [],               # tetrahedral has no trans pairs
    'TBP': [(0, 1)],          # only axial-axial is strictly 180°
    'SP':  [(1, 3), (2, 4)],  # trans basal pairs
    'PBP': [(0, 1)],          # axial-axial
    'SAP': [(0, 6), (1, 7), (2, 4), (3, 5)],  # square-antiprismatic
    'DD':  [(0, 2), (1, 3), (4, 6), (5, 7)],  # dodecahedral
}

_TOPO_CANONICAL_FNS = {
    'LIN': _canonical_lin,
    'TP':  _canonical_tp,
    'TS':  _canonical_ts,
    'OH':  _canonical_oh,
    'SQ':  _canonical_sq,
    'TH':  _canonical_th,
    'TBP': _canonical_tbp,
    'SP':  _canonical_sp,
    'PBP': _canonical_pbp,
    'SAP': _canonical_sap,
    'DD':  _canonical_dd,
}

# Idealized coordination vectors per geometry (bond length ~2 Å)
_TOPO_GEOMETRY_VECTORS: Dict[str, List[Tuple[float, float, float]]] = {
    'LIN': [(2, 0, 0), (-2, 0, 0)],
    'TP':  [(2, 0, 0), (-1, 1.732, 0), (-1, -1.732, 0)],
    'TS':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0)],
    'OH':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)],
    'SQ':  [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0)],
    'TH':  [(1.155, 1.155, 1.155), (-1.155, -1.155, 1.155),
            (-1.155, 1.155, -1.155), (1.155, -1.155, -1.155)],
    'TBP': [(0, 0, 2), (0, 0, -2), (2, 0, 0), (-1, 1.732, 0), (-1, -1.732, 0)],
    'SP':  [(0, 0, 2.2), (2, 0, 0.4), (0, 2, 0.4), (-2, 0, 0.4), (0, -2, 0.4)],
    'PBP': [(0, 0, 2), (0, 0, -2), (2, 0, 0), (0.618, 1.902, 0),
            (-1.618, 1.176, 0), (-1.618, -1.176, 0), (0.618, -1.902, 0)],
    # Square-antiprismatic (D4d): top square at z=+0.8, bottom at z=-0.8, rotated 45°
    'SAP': [(1.414, 0, 0.8), (0, 1.414, 0.8), (-1.414, 0, 0.8), (0, -1.414, 0.8),
            (1, 1, -0.8), (-1, 1, -0.8), (-1, -1, -0.8), (1, -1, -0.8)],
    # Dodecahedral (D2d): two interpenetrating trapezoids
    'DD':  [(1.414, 0, 1), (0, 1.414, -1), (-1.414, 0, 1), (0, -1.414, -1),
            (0.9, 0.9, 0), (-0.9, 0.9, 0), (-0.9, -0.9, 0), (0.9, -0.9, 0)],
}


def _enumerate_topological_isomers(
    donor_labels: List[str],
    n_coord: int,
    chelate_pairs: List[FrozenSet],
    metal_symbol: str = '',
) -> List[Tuple[tuple, List[int]]]:
    """Return all unique (canonical_form, permutation) pairs.

    ``perm[position] = donor_list_index`` — permutations where chelate-
    constrained donor pairs are never placed in trans positions.

    Args:
        donor_labels: Element symbol per donor atom (index = position in
            donor_indices list passed by the caller).
        n_coord: Coordination number (2–8).
        chelate_pairs: frozensets of donor-list indices that must stay cis.
        metal_symbol: Metal element symbol (used for CN=4 geometry preference).

    Returns:
        List of (canonical_form_tuple, perm_list) for every unique isomer.
    """
    import itertools

    if n_coord == 2:
        geometries = ['LIN']
    elif n_coord == 3:
        geometries = ['TP', 'TS']
    elif n_coord == 4:
        # Metal-aware geometry ordering: preferred geometry first
        pref = _PREFERRED_CN4_GEOMETRY.get(metal_symbol, 'SQ')
        other = 'TH' if pref == 'SQ' else 'SQ'
        geometries = [pref, other]
    elif n_coord == 5:
        geometries = ['TBP', 'SP']
    elif n_coord == 6:
        geometries = ['OH']
    elif n_coord == 7:
        geometries = ['PBP']
    elif n_coord == 8:
        geometries = ['SAP', 'DD']
    else:
        return []

    results: List[Tuple[tuple, List[int]]] = []
    seen_canonical: set = set()

    for geom_name in geometries:
        canonical_fn = _TOPO_CANONICAL_FNS[geom_name]
        trans_pos = _TOPO_TRANS_POSITIONS[geom_name]

        for perm in itertools.permutations(range(n_coord)):
            # Validate chelate constraints: paired donors must not sit trans
            valid = True
            for chelate in chelate_pairs:
                chelate_list = sorted(chelate)
                # Find the geometry positions of these two donors.
                # chelate_list contains donor-list indices (0..n_coord-1),
                # which are always present in perm (a permutation of the
                # same range), so .index() will always succeed.
                pos_i = perm.index(chelate_list[0])
                pos_j = perm.index(chelate_list[1])
                for ta, tb in trans_pos:
                    if (pos_i == ta and pos_j == tb) or (pos_i == tb and pos_j == ta):
                        valid = False
                        break
                if not valid:
                    break
            if not valid:
                continue

            # Build canonical form from donor element symbols at each position
            types = tuple(donor_labels[perm[pos]] for pos in range(n_coord))
            cf = canonical_fn(types)
            if cf not in seen_canonical:
                seen_canonical.add(cf)
                results.append((cf, list(perm)))

    return results


def _embed_fragment_procrustes(
    mol,
    metal_idx: int,
    frag_atom_indices: set,
    frag_donor_indices: List[int],
    target_positions: List[Tuple[float, float, float]],
    coords: List[Tuple[float, float, float]],
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

    # Embed the fragment
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.randomSeed = 42
    cid = AllChem.EmbedMolecule(frag_mol, params)
    if cid < 0:
        return False

    # Extract fragment coordinates
    frag_conf = frag_mol.GetConformer(cid)
    frag_coords = np.array([
        [frag_conf.GetAtomPosition(i).x,
         frag_conf.GetAtomPosition(i).y,
         frag_conf.GetAtomPosition(i).z]
        for i in range(frag_mol.GetNumAtoms())
    ])

    # Source: donor positions in the fragment embedding
    donor_new_indices = [old_to_new[d] for d in frag_donor_indices if d in old_to_new]
    if not donor_new_indices:
        return False

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
        # Single donor: align fragment-donor→COM vector to target direction
        frag_center = frag_coords.mean(axis=0)
        src_dir = src[0] - frag_center
        tgt_dir = tgt[0] - np.zeros(3)  # target from origin (metal)
        src_norm = np.linalg.norm(src_dir)
        tgt_norm = np.linalg.norm(tgt_dir)
        if src_norm < 1e-8 or tgt_norm < 1e-8:
            R = np.eye(3)
        else:
            src_dir /= src_norm
            tgt_dir /= tgt_norm
            v = np.cross(src_dir, tgt_dir)
            c = np.dot(src_dir, tgt_dir)
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


def _build_topology_xyz(
    mol,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    apply_uff: bool,
) -> Optional[str]:
    """Build a DELFIN XYZ for one topological arrangement.

    Places the metal at the origin, donor atoms at idealized geometry
    vectors, then attempts RDKit fragment embedding with Procrustes
    alignment for each ligand fragment.  Falls back to BFS placement
    if fragment embedding fails.  Optionally applies OB UFF refinement.

    Args:
        mol: RDKit mol (with H atoms, as from ``_prepare_mol_for_embedding``).
        metal_idx: Index of the metal atom in *mol*.
        donor_atom_indices: List of donor atom indices (length == n_coord).
        perm: ``perm[position] = index into donor_atom_indices``.
        geometry: Key in ``_TOPO_GEOMETRY_VECTORS`` ('OH', 'SQ', …).
        apply_uff: Whether to run OB UFF optimization after placement.

    Returns:
        DELFIN-format XYZ string, or None on failure.
    """
    try:
        vectors = _TOPO_GEOMETRY_VECTORS[geometry]
        n_atoms = mol.GetNumAtoms()
        coords: List[Tuple[float, float, float]] = [(0.0, 0.0, 0.0)] * n_atoms
        placed: set = set()

        # Metal at origin
        coords[metal_idx] = (0.0, 0.0, 0.0)
        placed.add(metal_idx)

        # Donors at geometry positions
        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        donor_target_map: Dict[int, Tuple[float, float, float]] = {}
        for pos_idx, donor_list_idx in enumerate(perm):
            donor_atom_idx = donor_atom_indices[donor_list_idx]
            donor_sym = mol.GetAtomWithIdx(donor_atom_idx).GetSymbol()
            bl = _get_ml_bond_length(metal_sym, donor_sym)
            vx, vy, vz = vectors[pos_idx]
            mag = math.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
            if mag > 1e-8:
                vx = vx / mag * bl
                vy = vy / mag * bl
                vz = vz / mag * bl
            coords[donor_atom_idx] = (vx, vy, vz)
            donor_target_map[donor_atom_idx] = (vx, vy, vz)
            placed.add(donor_atom_idx)

        # --- Fragment embedding approach ---
        # Decompose non-metal atoms into ligand fragments (connected components
        # after removing metal bonds), then embed each fragment with RDKit and
        # Procrustes-align the donor atoms to their target positions.
        frag_embed_placed: set = set()
        try:
            non_metal = {a.GetIdx() for a in mol.GetAtoms()
                         if a.GetSymbol() not in _METAL_SET}
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

            for frag in fragments:
                # Find donors in this fragment
                frag_donors = [d for d in donor_atom_indices if d in frag]
                if not frag_donors:
                    continue  # non-coordinating fragment, BFS will handle it
                tgt_positions = [donor_target_map[d] for d in frag_donors
                                 if d in donor_target_map]
                if not tgt_positions:
                    continue
                ok = _embed_fragment_procrustes(
                    mol, metal_idx, frag, frag_donors, tgt_positions, coords
                )
                if ok:
                    frag_embed_placed.update(frag)
                    placed.update(frag)
        except Exception as _frag_exc:
            logger.debug("Fragment embedding failed, using BFS fallback: %s", _frag_exc)

        # BFS fallback for any atoms not placed by fragment embedding
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
                mag_base = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                if mag_base < 1e-8:
                    dx, dy, dz = 1.0 + 0.1 * k, 0.3 * k, 0.0
                    mag_base = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                # Rotate around the base direction for each additional neighbour
                # so they fan out instead of all pointing the same way.
                if n_unplaced > 1 and k > 0:
                    angle = 2 * math.pi * k / n_unplaced
                    # Build a perpendicular vector to (dx,dy,dz)
                    ax, ay, az = dx / mag_base, dy / mag_base, dz / mag_base
                    if abs(ax) < 0.9:
                        px, py, pz = 1.0, 0.0, 0.0
                    else:
                        px, py, pz = 0.0, 1.0, 0.0
                    # Gram-Schmidt: subtract projection onto ax
                    dot_pa = px * ax + py * ay + pz * az
                    px -= dot_pa * ax; py -= dot_pa * ay; pz -= dot_pa * az
                    pm = math.sqrt(px ** 2 + py ** 2 + pz ** 2)
                    if pm > 1e-8:
                        px /= pm; py /= pm; pz /= pm
                    # Rodrigues rotation of (dx,dy,dz) by angle around (ax,ay,az)
                    cos_a = math.cos(angle); sin_a = math.sin(angle)
                    dx2 = dx*cos_a + (ay*dz - az*dy)*sin_a + ax*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dy2 = dy*cos_a + (az*dx - ax*dz)*sin_a + ay*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dz2 = dz*cos_a + (ax*dy - ay*dx)*sin_a + az*(ax*dx+ay*dy+az*dz)*(1-cos_a)
                    dx, dy, dz = dx2, dy2, dz2
                    mag_base = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                    if mag_base < 1e-8:
                        mag_base = 1.0
                dx = dx / mag_base * bond_len_default
                dy = dy / mag_base * bond_len_default
                dz = dz / mag_base * bond_len_default
                coords[nbr_idx] = (cx + dx, cy + dy, cz + dz)
                placed.add(nbr_idx)
                queue.append(nbr_idx)

        # Build XYZ string
        lines = []
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            x, y, z = coords[i]
            lines.append(f"{atom.GetSymbol():4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        xyz = '\n'.join(lines) + '\n'

        if apply_uff:
            try:
                # Free UFF optimization — do NOT constrain the metal or
                # M-L distances here.  Constraints prevent proper ligand
                # backbone relaxation and cause roundtrip connectivity
                # failures (structures get dropped by the topology filter).
                # The constraint API in _optimize_xyz_openbabel is available
                # for callers that need it, but the topology builder relies
                # on unconstrained relaxation for best results.
                xyz = _optimize_xyz_openbabel(xyz)
            except Exception as uff_exc:
                # Keep the generated topology geometry when UFF fails.
                # Dropping the isomer here can hide valid alternatives.
                logger.debug("Topology UFF optimization failed, keeping unoptimized XYZ: %s", uff_exc)

        return xyz
    except Exception as e:
        logger.debug("_build_topology_xyz failed: %s", e)
        return None


def _generate_topological_isomers(
    mol,
    smiles: str,
    apply_uff: bool = True,
    max_isomers: int = 50,
) -> List[Tuple[str, str]]:
    """Guarantee-complete isomer enumeration via topological permutation.

    For each metal center: enumerates all unique canonical arrangements of
    donor atoms (respecting chelate cis-constraints), builds an idealized
    3D structure for each, runs OB UFF, and labels the result.

    Returns [(xyz_string, label), …].
    """
    results: List[Tuple[str, str]] = []
    dtype_map = _donor_type_map(mol)

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        donor_indices = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        n_coord = len(donor_indices)

        if n_coord < 2 or n_coord > 8:
            continue

        # Use element symbols as donor labels for canonical form comparison
        donor_labels = [
            dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset()))[0]
            for d in donor_indices
        ]

        chelate_ps = _chelate_pairs(mol, metal_idx, donor_indices)

        # Skip macrocyclic complexes: when ALL donor pairs are chelate-connected
        # (complete chelate graph), every donor is part of one big ring.
        # The BFS in _build_topology_xyz cannot respect ring-closure constraints
        # and always produces broken structures for macrocycles.
        max_pairs = n_coord * (n_coord - 1) // 2
        if len(chelate_ps) >= max_pairs:
            logger.debug(
                "Skipping topo enumerator for macrocyclic complex (CN=%d, "
                "chelate_pairs=%d/%d)", n_coord, len(chelate_ps), max_pairs
            )
            continue

        # Convert chelate pairs from atom indices to donor-list indices.
        # _chelate_pairs returns frozensets of atom indices (e.g. {1, 12}),
        # but _enumerate_topological_isomers works with donor-list indices
        # (0..n_coord-1).  Without this mapping the constraint check always
        # hits ValueError → silently skipped → all constraints ignored.
        atom_to_listidx = {atom_idx: li for li, atom_idx in enumerate(donor_indices)}
        chelate_list_pairs: List[FrozenSet] = []
        for cp in chelate_ps:
            pair = sorted(cp)
            if pair[0] in atom_to_listidx and pair[1] in atom_to_listidx:
                chelate_list_pairs.append(frozenset([
                    atom_to_listidx[pair[0]], atom_to_listidx[pair[1]]
                ]))

        isomers = _enumerate_topological_isomers(
            donor_labels, n_coord, chelate_list_pairs,
            metal_symbol=atom.GetSymbol(),
        )

        for canonical_form, perm in isomers:
            if len(results) >= max_isomers:
                return results
            geom_name = canonical_form[0]  # 'OH', 'SQ', 'TH', 'LIN', 'TP', 'TS', …
            try:
                xyz = _build_topology_xyz(
                    mol, metal_idx, donor_indices, perm, geom_name, apply_uff
                )
                if xyz is None:
                    continue

                # Inject as temp conformer to run quality checks + fingerprint
                mol_tmp = Chem.RWMol(mol)
                mol_tmp.RemoveAllConformers()
                conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz)
                if conf is None:
                    # Atom count/order mismatch: XYZ is unreliable, skip it.
                    continue
                cid = mol_tmp.AddConformer(conf, assignId=True)
                try:
                    # Relaxed clash check for topo-generated structures:
                    # use min_dist=0.3 instead of default 0.5 since the
                    # topological enumerator guarantees correct topology
                    # and minor steric overlaps are expected before UFF.
                    # Skip _has_bad_geometry entirely — the topo enumerator
                    # places donors at ideal geometry positions, so angle-
                    # based rejection would discard valid structures
                    # (especially TH with 109.5° angles).
                    if _has_atom_clash(mol_tmp.GetMol(), cid, min_dist=0.3):
                        continue
                except Exception:
                    pass
                fp = _compute_coordination_fingerprint(
                    mol_tmp.GetMol(), cid, dtype_map=dtype_map
                )
                label = _classify_isomer_label(fp, mol_tmp.GetMol())

                # For CNs with multiple geometry types (e.g. CN=4: SQ + TH),
                # prefix the geometry name on secondary geometries to
                # distinguish them from the primary geometry's isomers.
                # Primary geometry labels stay unprefixed so they dedup
                # correctly against sampling results.
                _PRIMARY_GEOM_BASE = {
                    2: 'LIN', 3: 'TP', 5: 'TBP',
                    6: 'OH', 7: 'PBP', 8: 'SAP',
                }
                # CN=4: use metal-specific preference (d8 → SQ, else → TH)
                _PRIMARY_GEOM = dict(_PRIMARY_GEOM_BASE)
                _PRIMARY_GEOM[4] = _PREFERRED_CN4_GEOMETRY.get(
                    atom.GetSymbol(), 'SQ'
                )
                _GEOM_PRETTY = {
                    'LIN': 'linear', 'TP': 'trigonal-planar', 'TS': 'T-shaped',
                    'SQ': 'square-planar', 'TH': 'tetrahedral',
                    'TBP': 'trigonal-bipyramidal', 'SP': 'square-pyramidal',
                    'OH': 'octahedral', 'PBP': 'pentagonal-bipyramidal',
                    'SAP': 'square-antiprismatic', 'DD': 'dodecahedral',
                }
                primary = _PRIMARY_GEOM.get(n_coord)
                if primary and geom_name != primary:
                    geom_pretty = _GEOM_PRETTY.get(geom_name, geom_name)
                    if label:
                        label = f'{geom_pretty} {label}'
                    else:
                        label = geom_pretty

                results.append((xyz, label))
            except Exception as e:
                logger.debug("Topology isomer build failed (%s): %s", geom_name, e)
                continue

    return results


# ---------------------------------------------------------------------------
# Linkage isomers (Feature 2)
# ---------------------------------------------------------------------------

def _find_linkage_alternatives(mol) -> List[Tuple[int, int, int, str]]:
    """Detect ambidentate ligands and return alternative coordination modes.

    Supported patterns:
      - NO2⁻ (N-donor → O-donor, label 'nitrito')
      - SCN⁻ (S-donor → N-donor, label 'isothiocyanato-N')
      - CN⁻  (C-donor → N-donor, label 'isocyano')

    Returns list of (metal_idx, current_donor_idx, alt_donor_idx, label).
    """
    alternatives: List[Tuple[int, int, int, str]] = []
    metal_neighbor_sets: Dict[int, set] = {}

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        bonded = {nbr.GetIdx() for nbr in atom.GetNeighbors()}
        metal_neighbor_sets[metal_idx] = bonded

        for donor in atom.GetNeighbors():
            donor_idx = donor.GetIdx()
            donor_sym = donor.GetSymbol()

            # --- NO2 → nitrito (N-donor to O-donor) ---
            if donor_sym == 'N':
                o_nbrs = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O' and n.GetIdx() not in bonded
                ]
                if len(o_nbrs) >= 2:
                    alt_o = o_nbrs[0]
                    alternatives.append((metal_idx, donor_idx, alt_o.GetIdx(), 'nitrito'))

            # --- SCN → isothiocyanato-N (S-donor to N-donor) ---
            if donor_sym == 'S':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for n_nbr in c_nbr.GetNeighbors():
                            if (n_nbr.GetSymbol() == 'N'
                                    and n_nbr.GetIdx() != donor_idx
                                    and n_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, n_nbr.GetIdx(), 'isothiocyanato-N')
                                )

            # --- CN → isocyano (C-donor to N-donor via triple bond) ---
            if donor_sym == 'C':
                for n_nbr in donor.GetNeighbors():
                    if (n_nbr.GetSymbol() == 'N'
                            and n_nbr.GetIdx() not in bonded):
                        bond = mol.GetBondBetweenAtoms(donor_idx, n_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 2.5:
                            alternatives.append(
                                (metal_idx, donor_idx, n_nbr.GetIdx(), 'isocyano')
                            )

            # --- NCS (N-donor → S-donor via C): thiocyanato-S ---
            if donor_sym == 'N':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for s_nbr in c_nbr.GetNeighbors():
                            if (s_nbr.GetSymbol() == 'S'
                                    and s_nbr.GetIdx() != donor_idx
                                    and s_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, s_nbr.GetIdx(), 'thiocyanato-S')
                                )

            # --- Carboxylate (O1 → O2): alternative carboxylate oxygen ---
            if donor_sym == 'O':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C':
                        o_nbrs = [
                            n for n in c_nbr.GetNeighbors()
                            if n.GetSymbol() == 'O'
                            and n.GetIdx() != donor_idx
                            and n.GetIdx() not in bonded
                        ]
                        if o_nbrs:
                            alternatives.append(
                                (metal_idx, donor_idx, o_nbrs[0].GetIdx(), 'carboxylato-alt')
                            )

            # --- Sulfoxide (S-donor → O-donor via double bond) ---
            if donor_sym == 'S':
                o_nbrs = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O'
                    and n.GetIdx() not in bonded
                    and mol.GetBondBetweenAtoms(donor_idx, n.GetIdx()) is not None
                    and mol.GetBondBetweenAtoms(donor_idx, n.GetIdx()).GetBondTypeAsDouble() >= 1.5
                ]
                if o_nbrs:
                    alternatives.append(
                        (metal_idx, donor_idx, o_nbrs[0].GetIdx(), 'sulfoxide-O')
                    )

            # --- Pyrazole / Imidazole / Triazole (N1 → N2 in same ring) ---
            if donor_sym == 'N':
                try:
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if donor_idx in ring:
                            other_n = [
                                idx for idx in ring
                                if mol.GetAtomWithIdx(idx).GetSymbol() == 'N'
                                and idx != donor_idx
                                and idx not in bonded
                            ]
                            for alt_n in other_n:
                                alternatives.append(
                                    (metal_idx, donor_idx, alt_n, 'N-isomer')
                                )
                except Exception:
                    pass

    return alternatives


def _rewire_linkage(mol, metal_idx: int, old_donor: int, new_donor: int) -> Optional[object]:
    """Return a new mol with the M–old_donor bond replaced by M–new_donor.

    Returns None if the rewiring would duplicate an existing bond or fails.
    """
    try:
        rw = Chem.RWMol(mol)
        # Abort if new_donor is already bonded to metal
        if rw.GetBondBetweenAtoms(metal_idx, new_donor) is not None:
            return None
        rw.RemoveBond(metal_idx, old_donor)
        rw.AddBond(metal_idx, new_donor, Chem.BondType.SINGLE)
        rw.UpdatePropertyCache(strict=False)
        return rw.GetMol()
    except Exception as e:
        logger.debug("_rewire_linkage failed: %s", e)
        return None


def _generate_linkage_isomers(
    mol,
    smiles: str,
    apply_uff: bool = True,
) -> List[Tuple[str, str]]:
    """Build linkage isomers by rewiring ambidentate ligands and generating XYZ.

    For each alternative coordination mode, rewires the metal bond and builds
    an idealized topology structure (OB UFF optimized).

    Returns [(xyz_string, label), …] where label is e.g. 'nitrito'.
    """
    results: List[Tuple[str, str]] = []
    alternatives = _find_linkage_alternatives(mol)

    for metal_idx, old_donor, new_donor, type_label in alternatives:
        alt_mol = _rewire_linkage(mol, metal_idx, old_donor, new_donor)
        if alt_mol is None:
            continue
        try:
            # Find metal in alt_mol (same index)
            metal_atom = alt_mol.GetAtomWithIdx(metal_idx)
            donor_indices = [nbr.GetIdx() for nbr in metal_atom.GetNeighbors()]
            n_coord = len(donor_indices)
            geom_map = {2: 'LIN', 3: 'TP', 4: 'SQ', 5: 'TBP', 6: 'OH', 7: 'PBP'}
            if n_coord not in geom_map:
                continue
            geom = geom_map[n_coord]
            perm = list(range(n_coord))  # identity: donor i → position i
            xyz = _build_topology_xyz(
                alt_mol, metal_idx, donor_indices, perm, geom, apply_uff
            )
            if xyz is not None:
                results.append((xyz, type_label))
        except Exception as e:
            logger.debug("Linkage isomer (%s) failed: %s", type_label, e)

    return results


# ---------------------------------------------------------------------------
# Alternative binding-site exploration (Feature 5)
# ---------------------------------------------------------------------------

def _is_viable_donor(mol, atom_idx: int, metal_bonded: set) -> bool:
    """Check if an atom is a viable donor (has lone pairs available for coordination).

    Rejects atoms that are:
    - Double-bonded to carbon without additional lone pairs (C=O carbonyl)
    - Already bonded to a metal
    - Terminal hydrogen
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    sym = atom.GetSymbol()

    # Already bonded to a metal → not an *alternative* donor
    if any(nbr.GetSymbol() in _METAL_SET for nbr in atom.GetNeighbors()):
        return False

    # Standard heteroatom donors: N, O, S, P
    if sym in ('N', 'O', 'S', 'P'):
        # Count double bonds to carbon — carbonyl O (C=O) has only one neighbor
        # and that bond is a double bond, so it has no remaining lone pair for sigma donation
        if sym == 'O':
            nbrs = list(atom.GetNeighbors())
            if len(nbrs) == 1 and nbrs[0].GetSymbol() == 'C':
                bond = mol.GetBondBetweenAtoms(atom_idx, nbrs[0].GetIdx())
                if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                    return False  # carbonyl oxygen — not a viable donor
        return True

    # Carbon donors: cyclometalated aromatic C (ppy-type) or NHC carbens
    if sym == 'C':
        # Aromatic C that could be cyclometalated
        if atom.GetIsAromatic():
            return True
        # NHC-type carbene: C bonded to ≥2 N atoms
        n_nitrogen_nbrs = sum(
            1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'N'
        )
        if n_nitrogen_nbrs >= 2:
            return True

    return False


def _find_bridging_donors(mol) -> List[Tuple[int, List[int]]]:
    """Detect bridging donor atoms bound to ≥2 metal centers.

    Returns list of ``(donor_idx, [metal_idx_1, metal_idx_2, …])``.
    """
    bridging: List[Tuple[int, List[int]]] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METAL_SET:
            continue
        metal_nbrs = [
            nbr.GetIdx() for nbr in atom.GetNeighbors()
            if nbr.GetSymbol() in _METAL_SET
        ]
        if len(metal_nbrs) >= 2:
            bridging.append((atom.GetIdx(), metal_nbrs))
    return bridging


def _ligand_fragments(mol, metal_idx: int) -> List[set]:
    """Decompose non-metal atoms into ligand fragments (connected components after removing metal).

    Bridging donors (atoms bound to ≥2 metals) are included in ALL
    fragments they connect to, so each metal center's fragment list is
    complete.
    """
    non_metal = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() not in _METAL_SET}
    adj: Dict[int, set] = {i: set() for i in non_metal}
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in non_metal and bj in non_metal:
            adj[bi].add(bj)
            adj[bj].add(bi)
    visited: set = set()
    fragments: List[set] = []
    for start in sorted(non_metal):
        if start in visited:
            continue
        frag: set = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            frag.add(node)
            for nbr in adj.get(node, ()):
                if nbr not in visited:
                    stack.append(nbr)
        fragments.append(frag)

    # Bridging donors: ensure they appear in every fragment whose metal
    # they coordinate to.  This prevents incorrect fragment splitting
    # for μ-Cl, μ-OR, μ-oxo bridged complexes.
    bridging = _find_bridging_donors(mol)
    if bridging:
        for donor_idx, metal_list in bridging:
            for frag in fragments:
                # If any atom in frag is a neighbour of donor_idx (via adj),
                # the donor should be in this fragment.
                if donor_idx not in frag and frag & adj.get(donor_idx, set()):
                    frag.add(donor_idx)

    return fragments


def _generate_alternative_binding_modes(
    mol,
    smiles: str,
    apply_uff: bool = True,
    max_alternatives: int = 20,
) -> List[Tuple[str, str]]:
    """Generate isomers by swapping donor atoms with alternative binding sites.

    Only considers viable donors (atoms with available lone pairs) on the
    SAME ligand fragment as the current donor.  This avoids generating
    nonsensical structures from swapping donors across unrelated ligands
    or using C=O carbonyl oxygens as coordination donors.

    Returns [(xyz_string, label), ...].
    """
    results: List[Tuple[str, str]] = []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        bonded = {nbr.GetIdx() for nbr in atom.GetNeighbors()}
        current_donors = list(bonded)

        # Build ligand fragments to ensure we only swap within the same ligand
        fragments = _ligand_fragments(mol, metal_idx)
        atom_to_frag: Dict[int, int] = {}
        for fi, frag in enumerate(fragments):
            for aidx in frag:
                atom_to_frag[aidx] = fi

        # For each current donor, find viable alternative donors on the same fragment
        for current_d in current_donors:
            current_frag = atom_to_frag.get(current_d, -1)
            if current_frag < 0:
                continue
            current_sym = mol.GetAtomWithIdx(current_d).GetSymbol()

            for alt_d in fragments[current_frag]:
                if alt_d == current_d or alt_d in bonded:
                    continue
                if not _is_viable_donor(mol, alt_d, bonded):
                    continue
                if len(results) >= max_alternatives:
                    return results

                alt_sym = mol.GetAtomWithIdx(alt_d).GetSymbol()
                alt_mol = _rewire_linkage(mol, metal_idx, current_d, alt_d)
                if alt_mol is None:
                    continue
                try:
                    new_metal = alt_mol.GetAtomWithIdx(metal_idx)
                    donor_indices = [nbr.GetIdx() for nbr in new_metal.GetNeighbors()]
                    n_coord = len(donor_indices)
                    geom_map = {2: 'LIN', 3: 'TP', 4: 'SQ', 5: 'TBP', 6: 'OH', 7: 'PBP'}
                    if n_coord not in geom_map:
                        continue
                    geom = geom_map[n_coord]
                    perm = list(range(n_coord))
                    xyz = _build_topology_xyz(
                        alt_mol, metal_idx, donor_indices, perm, geom, apply_uff
                    )
                    if xyz is not None:
                        if current_sym == alt_sym:
                            label = f'alt-{alt_sym}-isomer'
                        else:
                            label = f'alt-bind-{alt_sym}'
                        results.append((xyz, label))
                except Exception as e:
                    logger.debug("Alternative binding mode failed: %s", e)
                    continue

    return results


def smiles_to_xyz_isomers(
    smiles: str,
    num_confs: int = 200,
    max_isomers: int = 50,
    apply_uff: bool = True,
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
        xyz, err = smiles_to_xyz(smiles, apply_uff=apply_uff)
        if err:
            return [], err
        return [(xyz, '')], None

    # Prepare molecule for embedding
    mol = _prepare_mol_for_embedding(smiles)
    if mol is None:
        # Fall back to single-conformer conversion
        xyz, err = smiles_to_xyz(smiles, apply_uff=apply_uff)
        if err:
            return [], err
        return [(xyz, '')], None

    # For metal complexes: prepend OB conformers to the pool so that
    # Avogadro-quality geometries are always considered during isomer search.
    has_metal = contains_metal(smiles)
    conf_ids: List[int] = []
    if has_metal and OPENBABEL_AVAILABLE:
        try:
            # Run OB conformer search in 3 independent restarts.
            # WeightedRotorSearch is non-deterministic (global PRNG); separate
            # calls from fresh OBMol instances give complementary pools.
            _n_ob_restarts = 3
            _per_ob = max(10, int(num_confs) // _n_ob_restarts)
            ob_xyz_blocks: List[str] = []
            _ob_seen: set = set()
            ob_error: Optional[str] = None
            for _restart in range(_n_ob_restarts):
                _blocks, _err = _openbabel_generate_conformer_xyz(
                    smiles, num_confs=_per_ob
                )
                if _blocks:
                    for _b in _blocks:
                        _key = "\n".join(
                            l.strip() for l in _b.splitlines() if l.strip()
                        )
                        if _key not in _ob_seen:
                            _ob_seen.add(_key)
                            ob_xyz_blocks.append(_b)
                elif _err and not ob_xyz_blocks:
                    ob_error = _err
            if ob_xyz_blocks:
                ob_ids = _inject_openbabel_conformers_into_mol(mol, ob_xyz_blocks)
                conf_ids.extend(ob_ids)
                logger.debug(
                    "OB conformers injected for isomer search: %d (3 restarts)",
                    len(ob_ids),
                )
            elif ob_error:
                logger.debug("OB conformer generation: %s", ob_error)
        except Exception as ob_exc:
            logger.debug("OB conformer generation exception: %s", ob_exc)
            mol.RemoveAllConformers()
            conf_ids = []

    # Embed multiple conformers with deterministic seed schedule.
    # This improves reproducibility while still sampling diverse geometries.
    try:
        # Fixed seeds keep results reproducible across runs.
        # Multiple diverse seeds improve sampling of coordination isomers.
        seeds = [31, 42, 7, 97, 13, 61, 83, 127, 211, 307, 401, 503]
        n_rounds = len(seeds)
        per_round = max(1, int(math.ceil(num_confs / n_rounds)))
        for seed in seeds:
            params = AllChem.ETKDGv3()
            params.useRandomCoords = True
            params.randomSeed = seed
            params.enforceChirality = False
            try:
                params.clearConfs = False
            except Exception:
                pass
            conf_ids.extend(
                list(AllChem.EmbedMultipleConfs(mol, numConfs=per_round, params=params))
            )
    except Exception as e:
        logger.warning("Multi-conformer embedding failed: %s", e)
        # Do not abort here: keep already injected OB conformers if available.
        # If none exist, continue with an empty pool so topo/linkage
        # enumeration can still add valid isomers.
        if not conf_ids:
            try:
                mol.RemoveAllConformers()
            except Exception:
                pass
            conf_ids = []

    if not conf_ids:
        logger.debug(
            "No conformers generated by OB/ETKDG; continuing with fallback + topo/linkage enumeration."
        )

    # Classify each conformer, skip obvious artifacts, then deduplicate by
    # full coordination fingerprint. This keeps distinct coordination
    # arrangements even when they share a coarse textual label.
    # Pre-compute donor types once (Morgan-based) to avoid repeated calls.
    dtype_map = _donor_type_map(mol)
    fp_label_pairs: List[Tuple[tuple, str, int, float]] = []
    for cid in conf_ids:
        try:
            # Hard-reject only truly collapsed structures (< 0.3 Å clash
            # or M-L < 1.2 Å); borderline issues become score penalties.
            if _has_atom_clash(mol, cid, min_dist=0.3):
                continue
            penalty = 0.0
            if _has_atom_clash(mol, cid):
                penalty += 500.0
            if _has_bad_geometry(mol, cid):
                penalty += 300.0
            if _has_ligand_intertwining(mol, cid):
                penalty += 200.0
            fp = _compute_coordination_fingerprint(mol, cid, dtype_map=dtype_map)
            score = _geometry_quality_score(mol, cid) + penalty
        except Exception:
            continue
        label = _classify_isomer_label(fp, mol)
        fp_label_pairs.append((fp, label, cid, score))

    # Second pass: deduplicate, keeping the best-scoring conformer per group.
    # fp -> (label, conf_id, score)
    seen_fps: Dict[tuple, Tuple[str, int, float]] = {}
    for fp, label, cid, score in fp_label_pairs:
        if fp not in seen_fps or score < seen_fps[fp][2]:
            seen_fps[fp] = (label, cid, score)
        if len(seen_fps) >= max_isomers:
            break

    # RMSD-based dedup: remove geometrically identical conformers that
    # slipped through fingerprint-based dedup (e.g. borderline angle
    # classifications producing different fingerprints for the same isomer).
    if len(seen_fps) > 1:
        def _base_label(lbl: str) -> str:
            if not lbl:
                return ''
            return re.sub(r'-\d+$', '', str(lbl))

        fps_list = list(seen_fps.keys())
        removed: set = set()
        for i in range(len(fps_list)):
            if i in removed:
                continue
            _li, cid_i, si = seen_fps[fps_list[i]]
            base_i = _base_label(_li)
            for j in range(i + 1, len(fps_list)):
                if j in removed:
                    continue
                _lj, cid_j, sj = seen_fps[fps_list[j]]
                base_j = _base_label(_lj)
                # Only collapse near-duplicates within the same label class.
                # Different labels (e.g. fac vs mer) are kept even at low RMSD.
                if base_i != base_j:
                    continue
                rmsd = _conformer_rmsd(mol, cid_i, cid_j)
                if rmsd < 1.5:
                    # Keep the conformer with the better geometry score
                    if si <= sj:
                        removed.add(j)
                    else:
                        removed.add(i)
                        break
        if removed:
            logger.debug("RMSD dedup removed %d duplicate(s)", len(removed))
            seen_fps = {fps_list[i]: seen_fps[fps_list[i]]
                        for i in range(len(fps_list)) if i not in removed}

    # Build results
    results: List[Tuple[str, str]] = []
    unknown_counter = 0
    if not seen_fps:
        # Sampling failed — get a single fallback geometry and still allow
        # the topological enumerator / linkage detector to augment it below.
        _fb_xyz, _fb_err = smiles_to_xyz(smiles, apply_uff=apply_uff)
        if _fb_err:
            return [], _fb_err
        results = [(_fb_xyz, '')]
    else:
        # Number duplicate labels (e.g. multiple "mer" with different fingerprints)
        label_counts: Dict[str, int] = {}
        for fp in seen_fps:
            lbl = seen_fps[fp][0] or ''
            label_counts[lbl] = label_counts.get(lbl, 0) + 1
        label_seen: Dict[str, int] = {}
        for fp, (label, cid, _score) in seen_fps.items():
            if not label:
                unknown_counter += 1
                display = f'Isomer {unknown_counter}'
            elif label_counts[label] > 1:
                label_seen[label] = label_seen.get(label, 0) + 1
                display = f'{label}-{label_seen[label]}'
            else:
                display = label
            try:
                xyz = _mol_to_xyz_conformer(mol, cid)
                if apply_uff:
                    try:
                        xyz = _optimize_xyz_openbabel(xyz)
                    except Exception as uff_exc:
                        # Preserve the isomer if UFF cannot optimize this geometry.
                        logger.debug(
                            "UFF optimization failed for conformer %s (%s), keeping unoptimized XYZ.",
                            cid, uff_exc,
                        )
            except Exception:
                continue
            # Topology check: organic ring count from OB XYZ must match original
            # SMILES (charge-insensitive: [N+]/[Fe-2] == [N]/[Fe] topologically).
            if not _roundtrip_ring_count_ok(xyz, smiles):
                logger.debug("Skipping conformer %d: topology mismatch", cid)
                continue
            # Bond check: reject conformers where OB perceives spurious bonds
            # between non-metal atoms (e.g. O-O, N-N) absent in original SMILES.
            if not _no_spurious_bonds(xyz, smiles):
                logger.debug("Skipping conformer %d: spurious bonds", cid)
                continue
            results.append((xyz, display))

        if not results:
            _fb_xyz, _fb_err = smiles_to_xyz(smiles, apply_uff=apply_uff)
            if _fb_err:
                return [], _fb_err
            results = [(_fb_xyz, '')]

    # --- Topological enumerator: guarantee completeness ---
    # Runs after sampling-based dedup; adds any isomers not found by sampling.
    # Dedup is now fingerprint-based (not label-based) so that distinct
    # topo-isomers with the same label (e.g. two different "mer" arrangements)
    # are not incorrectly suppressed.
    if has_metal:
        try:
            # Collect fingerprints from sampling results
            existing_fps: set = set()
            dtype_map_topo = _donor_type_map(mol)
            for existing_xyz, _existing_display in results:
                try:
                    mol_tmp = Chem.RWMol(mol)
                    mol_tmp.RemoveAllConformers()
                    conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), existing_xyz)
                    if conf is not None:
                        cid = mol_tmp.AddConformer(conf, assignId=True)
                        fp = _compute_coordination_fingerprint(
                            mol_tmp.GetMol(), cid, dtype_map=dtype_map_topo
                        )
                        existing_fps.add(fp)
                except Exception:
                    pass

            existing_displays = {display for _, display in results}

            topo_results = _generate_topological_isomers(
                mol, smiles, apply_uff=apply_uff, max_isomers=max_isomers
            )

            for topo_xyz, topo_label in topo_results:
                # Compute fingerprint of topo structure for dedup
                topo_fp = None
                try:
                    mol_tmp = Chem.RWMol(mol)
                    mol_tmp.RemoveAllConformers()
                    conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), topo_xyz)
                    if conf is not None:
                        cid = mol_tmp.AddConformer(conf, assignId=True)
                        topo_fp = _compute_coordination_fingerprint(
                            mol_tmp.GetMol(), cid, dtype_map=dtype_map_topo
                        )
                except Exception:
                    pass

                # Skip if this fingerprint was already seen
                if topo_fp is not None and topo_fp in existing_fps:
                    continue

                norm = topo_label or ''
                if not norm:
                    unknown_counter += 1
                    display = f'Isomer {unknown_counter}'
                else:
                    # Number duplicate labels
                    if norm in existing_displays:
                        suffix = 2
                        while f'{norm}-{suffix}' in existing_displays:
                            suffix += 1
                        display = f'{norm}-{suffix}'
                    else:
                        display = norm
                existing_displays.add(display)
                if topo_fp is not None:
                    existing_fps.add(topo_fp)
                results.append((topo_xyz, display))
        except Exception as _topo_exc:
            logger.debug("Topological isomer generation failed: %s", _topo_exc)

    # --- Linkage isomers ---
    if has_metal:
        try:
            link_results = _generate_linkage_isomers(mol, smiles, apply_uff=apply_uff)
            results.extend(link_results)
        except Exception as _link_exc:
            logger.debug("Linkage isomer generation failed: %s", _link_exc)

    # --- Alternative binding-site exploration ---
    if has_metal:
        try:
            existing_displays = {display for _, display in results}
            existing_base = {re.sub(r'-\d+$', '', d) for d in existing_displays}
            alt_results = _generate_alternative_binding_modes(
                mol, smiles, apply_uff=apply_uff
            )
            for alt_xyz, alt_label in alt_results:
                if alt_label not in existing_base:
                    results.append((alt_xyz, alt_label))
                    existing_base.add(alt_label)
        except Exception as _alt_exc:
            logger.debug("Alternative binding mode generation failed: %s", _alt_exc)

    # --- Topology validation: majority-vote SMILES filter ---
    # Convert every result XYZ back to canonical SMILES (metal bonds stripped)
    # and keep only structures whose topology matches the majority.
    # This catches fragmented or topology-corrupted structures from UFF /
    # topo-embedding / linkage permutation without relying on the input SMILES.
    if len(results) > 1 and OPENBABEL_AVAILABLE:
        try:
            smi_per_result: List[Optional[str]] = []
            for xyz, _label in results:
                smi_per_result.append(_xyz_to_canonical_smiles(xyz))

            # Count occurrences — only non-None SMILES participate in voting
            from collections import Counter
            valid_smiles = [s for s in smi_per_result if s is not None]
            if valid_smiles:
                counts = Counter(valid_smiles)
                majority_smi, majority_count = counts.most_common(1)[0]

                # quick XYZ as tie-breaker: if it matches a minority SMILES,
                # use that instead (quick usually preserves topology well)
                quick_xyz, _qerr = smiles_to_xyz_quick(smiles)
                if quick_xyz:
                    quick_smi = _xyz_to_canonical_smiles(quick_xyz)
                    if quick_smi and quick_smi in counts and quick_smi != majority_smi:
                        # quick disagrees with majority — trust quick only if
                        # the majority is very slim (≤50 % of valid results)
                        if majority_count <= len(valid_smiles) / 2:
                            majority_smi = quick_smi
                            logger.debug(
                                "Topology vote: quick SMILES overrides slim "
                                "majority (%d/%d)", majority_count, len(valid_smiles)
                            )

                filtered: List[Tuple[str, str]] = []
                n_dropped = 0
                for (xyz, label), smi in zip(results, smi_per_result):
                    if smi is None:
                        # Could not convert → keep (benefit of the doubt)
                        filtered.append((xyz, label))
                    elif smi == majority_smi:
                        filtered.append((xyz, label))
                    else:
                        n_dropped += 1
                        logger.debug(
                            "Topology filter dropped '%s': SMILES mismatch "
                            "('%s' vs majority '%s')",
                            label, smi[:80], majority_smi[:80],
                        )
                if filtered:
                    results = filtered
                if n_dropped:
                    logger.info(
                        "Topology filter: kept %d, dropped %d structure(s) "
                        "with inconsistent connectivity",
                        len(results), n_dropped,
                    )
        except Exception as _topo_filt_exc:
            logger.debug("Topology filter failed, keeping all results: %s", _topo_filt_exc)

    return results, None


def smiles_to_xyz_quick(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Fast single-conformer SMILES → XYZ (DELFIN format, no header).

    Uses the same strategy as the pre-Avogadro smiles_to_xyz: single
    ETKDGv3 embedding per strategy (no OB restarts, no multi-seed loop).
    For metal complexes delegates to _try_multiple_strategies which tries
    stk → RDKit → unsanitized → no-valence-check → OB in order and
    returns as soon as one succeeds.  Returns ``(xyz, error)``.
    """
    if not RDKIT_AVAILABLE:
        return None, "RDKit not available"

    has_metal = contains_metal(smiles)

    if has_metal:
        return _try_multiple_strategies(smiles)

    # Non-metal: single ETKDG embedding
    mol, note = mol_from_smiles_rdkit(smiles, allow_metal=False)
    if mol is None:
        return None, f"Could not parse SMILES: {note}"
    try:
        mol = Chem.AddHs(mol)
    except Exception:
        pass
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    result = AllChem.EmbedMolecule(mol, params)
    if result != 0:
        result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
    if result != 0:
        return None, f"Quick embedding failed for SMILES: {smiles[:60]}"
    return _mol_to_xyz(mol), None


def smiles_to_xyz(
    smiles: str,
    output_path: Optional[str] = None,
    apply_uff: bool = True,
) -> Tuple[Optional[str], Optional[str]]:
    """Convert SMILES string to XYZ coordinates using RDKit.

    Uses RDKit's ETKDG (Experimental Torsion Knowledge Distance Geometry) method
    for generating reasonable 3D conformers. This is suitable for initial geometries
    that will be further optimized with GOAT/xTB/ORCA.

    Supports coordination bonds (>) in SMILES for metal complexes.

    Args:
        smiles: SMILES string to convert
        output_path: Optional path to write XYZ file
        apply_uff: If True, apply UFF refinement (RDKit/Open Babel) where available.

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

        # Parse SMILES - try stk first for metal complexes, then RDKit.
        # IMPORTANT: Prefer the original SMILES before charge-normalized
        # variants to avoid introducing artificial oxidation states
        # (e.g., neutral Cu complexes rewritten as Cu+2).
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

        if mol is None:
            mol, rdkit_note = mol_from_smiles_rdkit(smiles, allow_metal=has_metal)
            method = "RDKit" if (mol is not None and rdkit_note is None) else (
                f"RDKit ({rdkit_note})" if mol is not None else None
            )
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

                # Try denormalized (neutral) SMILES as fallback
                if mol is None:
                    denormalized_smiles = _denormalize_metal_smiles(smiles)
                    if denormalized_smiles:
                        mol3, rdkit_note3 = mol_from_smiles_rdkit(denormalized_smiles, allow_metal=True)
                        if mol3 is not None:
                            mol = mol3
                            method = "RDKit (denormalized)"
                            rdkit_note = None

                if rdkit_note and ("Explicit valence" in rdkit_note or "kekulize" in rdkit_note):
                    legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                    if legacy_err is None and legacy_xyz:
                        if output_path:
                            Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                            logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                        return legacy_xyz, None
                # Last resort: try multi-strategy approach for metal complexes
                if has_metal:
                    logger.info("Trying multi-strategy fallback for unparseable metal SMILES")
                    return _try_multiple_strategies(smiles, output_path)
                error = f"Failed to parse SMILES string: {rdkit_note}"
                logger.error(error)
                return None, error

        # For metal complexes: convert bonds to dative and recalculate hydrogens
        # This fixes the issue where metal coordination bonds are counted towards ligand valence
        # Each step is independent so failure of one (e.g., RemoveHs on tetracoordinate B)
        # doesn't block dative conversion or AddHs.
        if has_metal:
            # Step 1: Remove existing explicit H atoms (stk may have added incorrect ones)
            try:
                mol = Chem.RemoveHs(mol)
            except Exception as e:
                logger.debug(f"RemoveHs skipped (non-standard valence, e.g. tetracoordinate B): {e}")

            # Step 2: Convert single bonds to metals to dative bonds
            try:
                mol = _convert_metal_bonds_to_dative(mol)
                mol.UpdatePropertyCache(strict=False)
            except Exception as e:
                logger.debug(f"Dative conversion skipped: {e}")

            # Step 3: Add hydrogens with correct valence calculation
            try:
                if _is_simple_organometallic(smiles):
                    mol = Chem.AddHs(mol, addCoords=True)
                    mol = _strip_h_on_metal_halogen(mol)
                    mol = _fix_organometallic_carbon_h(mol)
                else:
                    mol = Chem.AddHs(mol, addCoords=True)
            except Exception as e:
                logger.warning(f"Could not add hydrogens to metal complex: {e}")
                # Fallback for valence errors: try unsanitized path
                if "Explicit valence" in str(e):
                    legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                    if legacy_err is None and legacy_xyz:
                        if output_path:
                            Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                            logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                        return legacy_xyz, None

            # Step 4: Fix atoms with non-standard valence (e.g., tetracoordinate B
            # in pyrazolylborate/scorpionate ligands) to prevent embedding failures.
            # B with 4 bonds is chemically B- (borate) - set charge so RDKit
            # accepts valence 4 during embedding.
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'B' and atom.GetDegree() >= 4:
                    atom.SetNoImplicit(True)
                    if atom.GetFormalCharge() == 0:
                        atom.SetFormalCharge(-1)
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

        # Generate 3D coordinates using a hybrid OB+RDKit conformer pool.
        # For metal complexes: OB WeightedRotorSearch in 3 independent restarts
        # (Avogadro-quality, non-deterministic diversity) + RDKit ETKDG with
        # 12 diverse fixed seeds → up to ~500 conformers total.
        # Best geometry is selected by _geometry_quality_score.
        result = -1
        if has_metal:
            all_conf_ids: List[int] = []

            # --- OB conformers (Avogadro-equivalent pipeline) ---
            ob_injection_ok = False
            if OPENBABEL_AVAILABLE:
                try:
                    _n_ob_r = 3
                    _per_ob_r = max(10, 200 // _n_ob_r)
                    ob_xyz_blocks = []
                    _ob_seen_s: set = set()
                    ob_err: Optional[str] = None
                    for _ri in range(_n_ob_r):
                        _bl, _er = _openbabel_generate_conformer_xyz(
                            smiles, num_confs=_per_ob_r
                        )
                        if _bl:
                            for _b in _bl:
                                _k = "\n".join(
                                    l.strip() for l in _b.splitlines() if l.strip()
                                )
                                if _k not in _ob_seen_s:
                                    _ob_seen_s.add(_k)
                                    ob_xyz_blocks.append(_b)
                        elif _er and not ob_xyz_blocks:
                            ob_err = _er
                    if ob_xyz_blocks:
                        ob_ids = _inject_openbabel_conformers_into_mol(mol, ob_xyz_blocks)
                        if ob_ids:
                            all_conf_ids.extend(ob_ids)
                            ob_injection_ok = True
                            logger.debug(
                                "OB conformer injection: %d conformers (3 restarts) for %s",
                                len(ob_ids), smiles[:40],
                            )
                        else:
                            logger.debug("OB atom-order mismatch; using RDKit-only pool")
                            mol.RemoveAllConformers()
                    elif ob_err:
                        logger.debug("OB conformer generation: %s", ob_err)
                except Exception as ob_exc:
                    logger.debug("OB conformer generation exception: %s", ob_exc)
                    mol.RemoveAllConformers()

            # --- RDKit ETKDG with 12 diverse fixed seeds (~17 confs/seed) ---
            seeds = [31, 42, 7, 97, 13, 61, 83, 127, 211, 307, 401, 503]
            per_seed = max(1, 200 // len(seeds))
            try:
                for seed in seeds:
                    params_multi = AllChem.ETKDGv3()
                    params_multi.randomSeed = seed
                    params_multi.useRandomCoords = True
                    params_multi.enforceChirality = False
                    try:
                        params_multi.clearConfs = False  # append to OB conformers
                    except Exception:
                        pass
                    new_ids = list(AllChem.EmbedMultipleConfs(
                        mol, numConfs=per_seed, params=params_multi
                    ))
                    all_conf_ids.extend(new_ids)
            except Exception:
                pass

            if all_conf_ids:
                best_conf = None
                best_score = float("inf")
                for cid in all_conf_ids:
                    try:
                        if _has_atom_clash(mol, cid):
                            continue
                        if _has_bad_geometry(mol, cid):
                            continue
                        score = _geometry_quality_score(mol, cid)
                        if score < best_score:
                            best_score = score
                            best_conf = cid
                    except Exception:
                        continue

                if best_conf is None:
                    best_conf = all_conf_ids[0]

                # Keep only the selected conformer so downstream conversion
                # can continue using the default conformer accessor.
                selected = Chem.Conformer(mol.GetConformer(int(best_conf)))
                mol.RemoveAllConformers()
                mol.AddConformer(selected, assignId=True)
                result = 0

        if result != 0:
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
            # Try multi-strategy approach before giving up
            if has_metal:
                logger.info("Trying multi-strategy fallback after embedding failure")
                multi_xyz, multi_err = _try_multiple_strategies(smiles, output_path)
                if multi_xyz:
                    return multi_xyz, None
                # Last resort: manual coordinate construction for metal complexes
                logger.info("Trying manual metal coordinate construction")
                manual_xyz, manual_err = _manual_metal_embed(smiles)
                if manual_xyz:
                    if output_path:
                        Path(output_path).write_text(manual_xyz, encoding='utf-8')
                        logger.info(f"Converted SMILES to XYZ using manual metal embed: {output_path}")
                    logger.warning(
                        "Used manual coordinate construction - geometry is very rough! "
                        "GOAT/xTB optimization is essential before any calculations."
                    )
                    return manual_xyz, None
            error = f"Failed to generate 3D coordinates for SMILES: {smiles}"
            if legacy_err:
                error = f"{error} (legacy fallback failed: {legacy_err})"
            logger.error(error)
            return None, error

        # Optional UFF refinement for better starting structures
        if apply_uff and not has_metal:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                logger.debug("RDKit UFF optimization successful")
            except Exception as e:
                logger.info(f"RDKit UFF optimization skipped: {e}")

        # Convert to XYZ format
        xyz_content = _mol_to_xyz(mol)

        # For metal complexes, Open Babel UFF has full metal parameters.
        if apply_uff and has_metal:
            xyz_content = _optimize_xyz_openbabel(xyz_content)

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
        # Last resort: multi-strategy fallback for metal complexes
        if has_metal:
            logger.info("Trying multi-strategy fallback after exception: %s", e)
            multi_xyz, multi_err = _try_multiple_strategies(smiles, output_path)
            if multi_xyz:
                return multi_xyz, None
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


def _optimize_xyz_openbabel(
    xyz_delfin: str,
    steps: int = 500,
    constraints: Optional[Dict] = None,
) -> str:
    """Optimize a DELFIN-format XYZ string using Open Babel's UFF force field.

    Open Babel's UFF implementation includes full parameters for transition
    metals (Fe, Co, Ni, Ti, etc.) that RDKit lacks.  This makes it ideal as
    a post-ETKDG refinement step for metal complexes.

    Args:
        xyz_delfin: Coordinate block in DELFIN format (``symbol x y z`` per line,
            no atom-count header).
        steps: Maximum number of conjugate-gradient steps.
        constraints: Optional dict with keys:
            - ``'fix_atoms'``: list of 0-based atom indices to freeze in place
            - ``'distances'``: list of ``(idx_a, idx_b, target_dist)``
            - ``'angles'``: list of ``(idx_a, idx_b, idx_c, target_angle_deg)``

    Returns:
        Optimized DELFIN-format XYZ string, or the original string if
        optimization fails for any reason.
    """
    if not OPENBABEL_AVAILABLE:
        return xyz_delfin

    try:
        # Parse DELFIN lines (skip blank lines)
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        n_atoms = len(lines)
        if n_atoms == 0:
            return xyz_delfin

        # Build standard XYZ string (atom count + comment + coordinates)
        std_xyz = f"{n_atoms}\n\n"
        for line in lines:
            parts = line.split()
            # DELFIN format: "symbol x y z"
            std_xyz += f"{parts[0]}  {parts[1]}  {parts[2]}  {parts[3]}\n"

        # Read into Open Babel
        ob_mol = pybel.readstring("xyz", std_xyz)

        # Run UFF conjugate-gradient optimization
        ff = pybel._forcefields["uff"]
        if not ff.Setup(ob_mol.OBMol):
            logger.debug("Open Babel UFF setup failed, returning unoptimized geometry")
            return xyz_delfin

        # Apply constraints (if provided) to preserve coordination geometry
        if constraints:
            try:
                ob_constraints = pybel.ob.OBFFConstraints()
                # Fix atom positions (1-based indices for OB)
                for idx in constraints.get('fix_atoms', []):
                    ob_constraints.AddAtomConstraint(idx + 1)
                # Distance constraints
                for idx_a, idx_b, target in constraints.get('distances', []):
                    ob_constraints.AddDistanceConstraint(
                        idx_a + 1, idx_b + 1, target
                    )
                # Angle constraints
                for idx_a, idx_b, idx_c, target in constraints.get('angles', []):
                    ob_constraints.AddAngleConstraint(
                        idx_a + 1, idx_b + 1, idx_c + 1, target
                    )
                ff.SetConstraints(ob_constraints)
            except Exception as e:
                logger.debug("OB UFF constraint setup failed, running unconstrained: %s", e)

        ff.ConjugateGradients(steps)
        ff.GetCoordinates(ob_mol.OBMol)

        # Convert back to DELFIN format
        opt_lines = []
        for ob_atom in pybel.ob.OBMolAtomIter(ob_mol.OBMol):
            symbol = pybel.ob.GetSymbol(ob_atom.GetAtomicNum())
            x, y, z = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()
            opt_lines.append(f"{symbol:4s} {x:12.6f} {y:12.6f} {z:12.6f}")

        logger.debug("Open Babel UFF optimization completed (%d steps max)", steps)
        return '\n'.join(opt_lines) + '\n'

    except Exception as e:
        logger.debug("Open Babel UFF optimization failed: %s", e)
        return xyz_delfin


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
