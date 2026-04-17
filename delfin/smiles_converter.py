"""SMILES to XYZ conversion using RDKit with metal complex support.

Note: RDKit generates rough starting geometries that are NOT chemically accurate,
especially for metal complexes. These coordinates should ALWAYS be optimized with
GOAT/xTB before running ORCA calculations.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import concurrent.futures
import math
import os
import re
import threading
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)
_HAPTO_QUICK_PREVIEW_CACHE: Dict[str, List[Tuple[str, str]]] = {}


@dataclass
class _HybridHaptoFragment:
    """Metal-free fragment plus atom mappings for hybrid hapto assembly."""

    atom_indices: List[int]
    donor_atom_indices: List[int]
    metal_neighbor_indices: List[int]
    bridging_donor_indices: List[int]
    hapto_group_ids: List[int]
    anchor_atom_indices: List[int]
    original_to_fragment: Dict[int, int]
    fragment_to_original: Dict[int, int]
    fragment_mol: object
    use_scaffold_only: bool = False


@dataclass
class _HybridHaptoDecomposition:
    """Fragmented view of a metal complex with hapto metadata preserved."""

    metal_indices: List[int]
    donor_atom_indices: List[int]
    hapto_groups: List[Tuple[int, List[int]]]
    fragments: List[_HybridHaptoFragment]


@dataclass
class _PrimaryOrganometalModule:
    """Primary hapto-metal module with eligible non-hapto donor fragments."""

    metal_idx: int
    hapto_group_ids: List[int]
    donor_atom_indices: List[int]
    correlated_fragment_indices: List[int]
    terminal_fragment_indices: List[int]

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D
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
# Period-aware offset for M-L bond length fallback (Å).
# 5d metals are shorter than expected from covalent radii due to
# relativistic contraction; main-group metals are ionic with larger gaps.
# ---------------------------------------------------------------------------
_METAL_ROW_OFFSET: Dict[str, float] = {}
for _m in ('Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'):
    _METAL_ROW_OFFSET[_m] = 0.40
for _m in ('Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'):
    _METAL_ROW_OFFSET[_m] = 0.45
for _m in ('La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'):
    _METAL_ROW_OFFSET[_m] = 0.30
for _m in ('Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
           'Er', 'Tm', 'Yb', 'Lu'):
    _METAL_ROW_OFFSET[_m] = 0.45
for _m in ('Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'):
    _METAL_ROW_OFFSET[_m] = 0.40
for _m in ('Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba',
           'Al', 'Ga', 'In', 'Tl', 'Sn', 'Pb', 'Bi', 'Po'):
    _METAL_ROW_OFFSET[_m] = 0.50

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
    # Yttrium
    ('Y', 'C'): 2.55, ('Y', 'N'): 2.35, ('Y', 'O'): 2.25, ('Y', 'Cl'): 2.60,
    # Scandium
    ('Sc', 'C'): 2.30, ('Sc', 'N'): 2.15, ('Sc', 'O'): 2.05, ('Sc', 'Cl'): 2.40,
    # Chromium (carbon)
    ('Cr', 'C'): 2.05,
    # Vanadium
    ('V', 'C'): 2.15, ('V', 'N'): 2.10, ('V', 'O'): 2.00, ('V', 'Cl'): 2.35,
    # Manganese (carbon)
    ('Mn', 'C'): 2.10,
    # Titanium (carbon)
    ('Ti', 'C'): 2.25,
    # Zirconium (carbon)
    ('Zr', 'C'): 2.40,
    # Hafnium
    ('Hf', 'C'): 2.35, ('Hf', 'N'): 2.20, ('Hf', 'O'): 2.10, ('Hf', 'Cl'): 2.40,
    # Ruthenium (S)
    ('Ru', 'S'): 2.35,
    # Cadmium (d10, larger ionic radius than first-row TMs)
    ('Cd', 'C'): 2.25, ('Cd', 'N'): 2.33, ('Cd', 'O'): 2.30,
    ('Cd', 'P'): 2.55, ('Cd', 'Cl'): 2.55, ('Cd', 'S'): 2.55,
    ('Cd', 'Br'): 2.68, ('Cd', 'F'): 2.22,
    # Mercury (d10, even larger)
    ('Hg', 'C'): 2.12, ('Hg', 'N'): 2.20, ('Hg', 'O'): 2.25,
    ('Hg', 'P'): 2.48, ('Hg', 'Cl'): 2.50, ('Hg', 'S'): 2.45,
    ('Hg', 'Br'): 2.62,
    # Alkaline earth carboxylate/aqua coordination
    ('Ca', 'N'): 2.55, ('Ca', 'O'): 2.40, ('Ca', 'Cl'): 2.75,
    ('Mg', 'N'): 2.20, ('Mg', 'O'): 2.05, ('Mg', 'Cl'): 2.45,
    ('Sr', 'N'): 2.65, ('Sr', 'O'): 2.55, ('Sr', 'Cl'): 2.85,
    ('Ba', 'N'): 2.85, ('Ba', 'O'): 2.75, ('Ba', 'Cl'): 3.05,
    # Zn completions
    ('Zn', 'C'): 2.05, ('Zn', 'P'): 2.35, ('Zn', 'Br'): 2.35, ('Zn', 'F'): 1.95,
    # First-row TM: P, S, Br, F completions
    ('Sc', 'P'): 2.55, ('Sc', 'S'): 2.50, ('Sc', 'Br'): 2.65, ('Sc', 'F'): 2.00,
    ('Ti', 'P'): 2.50, ('Ti', 'S'): 2.40, ('Ti', 'Br'): 2.55, ('Ti', 'F'): 1.90,
    ('V', 'P'): 2.35, ('V', 'S'): 2.35, ('V', 'Br'): 2.50, ('V', 'F'): 1.85,
    ('Cr', 'P'): 2.30, ('Cr', 'S'): 2.35, ('Cr', 'Br'): 2.50, ('Cr', 'F'): 1.90,
    ('Mn', 'P'): 2.25, ('Mn', 'S'): 2.35, ('Mn', 'Br'): 2.50, ('Mn', 'F'): 1.90,
    ('Fe', 'Br'): 2.40, ('Fe', 'F'): 1.90,
    ('Co', 'P'): 2.18, ('Co', 'S'): 2.25, ('Co', 'Br'): 2.38, ('Co', 'F'): 1.88,
    ('Ni', 'S'): 2.20, ('Ni', 'Br'): 2.35, ('Ni', 'F'): 1.85,
    ('Cu', 'Br'): 2.40, ('Cu', 'F'): 1.90,
    # Second-row TM completions
    ('Zr', 'P'): 2.60, ('Zr', 'S'): 2.55, ('Zr', 'Br'): 2.60, ('Zr', 'F'): 2.05,
    ('Nb', 'N'): 2.15, ('Nb', 'O'): 2.05, ('Nb', 'Cl'): 2.40,
    ('Nb', 'P'): 2.50, ('Nb', 'S'): 2.50,
    ('Mo', 'P'): 2.45, ('Mo', 'Br'): 2.55, ('Mo', 'F'): 2.00,
    ('Tc', 'N'): 2.10, ('Tc', 'O'): 2.05, ('Tc', 'Cl'): 2.35,
    ('Ru', 'Br'): 2.50, ('Ru', 'F'): 2.00,
    ('Rh', 'S'): 2.35, ('Rh', 'Br'): 2.50, ('Rh', 'F'): 2.00,
    ('Pd', 'Br'): 2.45, ('Pd', 'F'): 2.00,
    ('Ag', 'Br'): 2.50, ('Ag', 'F'): 2.15,
    # Third-row TM completions
    ('Hf', 'P'): 2.55, ('Hf', 'S'): 2.50, ('Hf', 'Br'): 2.55, ('Hf', 'F'): 2.00,
    ('Ta', 'N'): 2.10, ('Ta', 'O'): 2.00, ('Ta', 'Cl'): 2.35,
    ('Ta', 'P'): 2.50, ('Ta', 'S'): 2.45,
    ('W', 'P'): 2.45, ('W', 'S'): 2.40, ('W', 'Br'): 2.55, ('W', 'F'): 2.00,
    ('Re', 'P'): 2.40, ('Re', 'S'): 2.40, ('Re', 'Br'): 2.50, ('Re', 'F'): 2.00,
    ('Os', 'S'): 2.35, ('Os', 'Br'): 2.50, ('Os', 'F'): 2.00,
    ('Ir', 'Br'): 2.50, ('Ir', 'F'): 2.00,
    ('Pt', 'F'): 2.00,
    ('Au', 'Br'): 2.40, ('Au', 'F'): 2.00, ('Au', 'O'): 2.10,
    # Selenium/Tellurium donors
    ('Fe', 'Se'): 2.40, ('Ru', 'Se'): 2.45, ('Pd', 'Se'): 2.40,
    ('Pt', 'Se'): 2.40, ('Cu', 'Se'): 2.40, ('Ni', 'Se'): 2.30,
    # Alkali (hard Lewis acids)
    ('Li', 'N'): 2.10, ('Li', 'O'): 1.95, ('Li', 'Cl'): 2.35,
    ('Na', 'N'): 2.45, ('Na', 'O'): 2.35, ('Na', 'Cl'): 2.70,
    ('K', 'N'): 2.85, ('K', 'O'): 2.75, ('K', 'Cl'): 3.10,
}

# ---------------------------------------------------------------------------
# Metal-Metal bond lengths (Å) from CSD averages.
# Keyed as frozenset({sym1, sym2}) → distance to handle both orderings.
# ---------------------------------------------------------------------------
_METAL_METAL_BOND_LENGTHS: Dict[frozenset, float] = {
    frozenset({'Cr'}): 1.83,  # Cr-Cr quintuple bond
    frozenset({'Mo'}): 2.10,  # Mo-Mo quadruple bond
    frozenset({'W'}): 2.20,
    frozenset({'Re'}): 2.24,
    frozenset({'Ru'}): 2.28,
    frozenset({'Rh'}): 2.40,
    frozenset({'Ir'}): 2.55,
    frozenset({'Pd'}): 2.58,
    frozenset({'Pt'}): 2.60,
    frozenset({'Cu'}): 2.55,
    frozenset({'Ag'}): 2.90,
    frozenset({'Au'}): 2.70,  # aurophilic
    frozenset({'Fe'}): 2.50,
    frozenset({'Co'}): 2.50,
    frozenset({'Ni'}): 2.45,
    frozenset({'Mn'}): 2.60,
    frozenset({'Ti'}): 2.75,
    frozenset({'Zr'}): 2.90,
    frozenset({'Zn'}): 2.90,
    frozenset({'Cd'}): 3.10,
    # Heterometallic (common pairs)
    frozenset({'Cu', 'Fe'}): 2.55,
    frozenset({'Mo', 'Cu'}): 2.65,
    frozenset({'Rh', 'Ir'}): 2.50,
    frozenset({'Pt', 'Pd'}): 2.60,
}

# ---------------------------------------------------------------------------
# Crystallographic Metal-Centroid distances (Å) for hapto coordination.
# Keyed as (metal_symbol, eta) → centroid distance.
# Values from CSD averages / literature data.
# ---------------------------------------------------------------------------
_HAPTO_CENTROID_DISTANCES: Dict[Tuple[str, int], float] = {
    # Iron
    ('Fe', 5): 1.65, ('Fe', 6): 1.55, ('Fe', 3): 1.90, ('Fe', 2): 2.05,
    ('Fe', 4): 1.75,
    # Ruthenium
    ('Ru', 5): 1.82, ('Ru', 6): 1.70, ('Ru', 3): 2.00, ('Ru', 2): 2.10,
    # Osmium
    ('Os', 5): 1.84, ('Os', 6): 1.72,
    # Chromium
    ('Cr', 5): 1.80, ('Cr', 6): 1.61, ('Cr', 3): 1.95, ('Cr', 2): 2.10,
    # Manganese
    ('Mn', 5): 1.78, ('Mn', 6): 1.65,
    # Vanadium
    ('V', 5): 1.90, ('V', 6): 1.75,
    # Titanium
    ('Ti', 5): 2.04, ('Ti', 6): 1.90, ('Ti', 3): 2.15,
    # Zirconium
    ('Zr', 5): 2.20, ('Zr', 6): 2.05, ('Zr', 3): 2.30,
    # Hafnium
    ('Hf', 5): 2.18, ('Hf', 6): 2.03,
    # Yttrium
    ('Y', 5): 2.35, ('Y', 6): 2.30, ('Y', 3): 2.45,
    # Scandium
    ('Sc', 5): 2.15, ('Sc', 6): 2.00,
    # Cobalt
    ('Co', 5): 1.66, ('Co', 3): 1.88, ('Co', 2): 2.02,
    # Rhodium
    ('Rh', 5): 1.85, ('Rh', 6): 1.75, ('Rh', 3): 1.95, ('Rh', 2): 2.05,
    # Iridium
    ('Ir', 5): 1.85, ('Ir', 6): 1.75, ('Ir', 3): 1.95, ('Ir', 2): 2.05,
    # Nickel
    ('Ni', 5): 1.74, ('Ni', 3): 1.88, ('Ni', 2): 1.98,
    # Palladium
    ('Pd', 3): 2.00, ('Pd', 2): 2.10,
    # Platinum
    ('Pt', 2): 2.02, ('Pt', 3): 2.00,
    # Molybdenum
    ('Mo', 5): 1.95, ('Mo', 6): 1.78, ('Mo', 3): 2.05,
    # Tungsten
    ('W', 5): 1.95, ('W', 6): 1.78,
    # Rhenium
    ('Re', 5): 1.90,
    # Lanthanides
    ('La', 5): 2.55, ('Ce', 5): 2.50, ('Nd', 5): 2.45, ('Sm', 5): 2.40,
    ('Gd', 5): 2.38, ('Lu', 5): 2.25,
    # Actinides
    ('U', 5): 2.45, ('U', 8): 1.92, ('Th', 5): 2.50, ('Th', 8): 2.00,
}


def _target_mc_dist(metal_sym: str, eta: int) -> float:
    """Return target metal-centroid distance for a hapto group.

    Lookup order:
    1. Exact match in _HAPTO_CENTROID_DISTANCES
    2. Geometric estimate from M-C bond length and ring radius
    3. Conservative fallback based on covalent radius
    """
    key = (metal_sym, eta)
    if key in _HAPTO_CENTROID_DISTANCES:
        return _HAPTO_CENTROID_DISTANCES[key]

    # Geometric estimate: M-Centroid = sqrt(M_C^2 - R_ring^2)
    try:
        base_mc = float(_get_ml_bond_length(metal_sym, 'C'))
    except Exception:
        base_mc = 2.20
    if not math.isfinite(base_mc):
        base_mc = 2.20

    if eta >= 3:
        ring_radius = 1.40 / (2.0 * math.sin(math.pi / eta))
        mc_sq = base_mc * base_mc - ring_radius * ring_radius
        if mc_sq > 0:
            dist = math.sqrt(mc_sq)
        else:
            dist = 0.80 * base_mc
    else:
        # eta2: offset from midpoint of C=C bond (~1.38 A / 2 = 0.69 A)
        mc_sq = base_mc * base_mc - 0.69 * 0.69
        if mc_sq > 0:
            dist = math.sqrt(mc_sq)
        else:
            dist = 0.85 * base_mc

    # Eta-specific clamps
    if eta >= 5:
        dist = max(dist, 1.50)
    elif eta == 4:
        dist = max(dist, 1.45)
    elif eta == 3:
        dist = max(dist, 1.40)
    else:
        dist = max(dist, 1.35)

    return dist


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
        offset = _METAL_ROW_OFFSET.get(metal_symbol, 0.5)
        return r_m + r_d + offset
    return 2.0


def _secondary_donor_target_length(mol, metal_symbol: str, donor_idx: int) -> float:
    """Return an environment-adjusted secondary-metal donor distance.

    The base value comes from the generic metal-ligand table and is then nudged
    toward common experimental average distances for planar N/O donors, anionic
    donors, and carbon donors in conjugated coordination environments.
    """
    if mol is None:
        return float(_get_ml_bond_length(metal_symbol, 'C'))

    atom = mol.GetAtomWithIdx(donor_idx)
    donor_sym = atom.GetSymbol()
    target = float(_get_ml_bond_length(metal_symbol, donor_sym))
    formal_charge = int(atom.GetFormalCharge())
    aromatic = bool(atom.GetIsAromatic())
    try:
        hyb = atom.GetHybridization()
    except Exception:
        hyb = None
    planar_like = aromatic or hyb in {
        Chem.rdchem.HybridizationType.SP,
        Chem.rdchem.HybridizationType.SP2,
    }
    conjugated = any(bond.GetIsConjugated() for bond in atom.GetBonds())
    multiple_bonded = any(bond.GetBondTypeAsDouble() >= 1.5 for bond in atom.GetBonds())

    if donor_sym == 'N':
        if planar_like:
            target -= 0.02
        if formal_charge < 0:
            target -= 0.04
        elif formal_charge > 0 and not aromatic:
            target += 0.03
    elif donor_sym == 'O':
        if planar_like or conjugated or multiple_bonded:
            target -= 0.03
        if formal_charge < 0:
            target -= 0.07
        elif formal_charge > 0 and not conjugated:
            target += 0.02
    elif donor_sym == 'S':
        if planar_like or conjugated:
            target -= 0.02
        if formal_charge < 0:
            target -= 0.05
    elif donor_sym == 'P':
        if formal_charge < 0:
            target -= 0.03
    elif donor_sym == 'C':
        if planar_like:
            target -= 0.01
        if formal_charge < 0:
            target -= 0.05
        elif formal_charge > 0:
            target += 0.03

    base = float(_get_ml_bond_length(metal_symbol, donor_sym))
    return float(min(max(target, base - 0.12), base + 0.08))


def _secondary_donor_fit_weight(mol, metal_symbol: str, donor_idx: int) -> float:
    """Return a relative fit weight for secondary-metal donors.

    Experimental averages are typically tighter for heterodonors in planar
    conjugated fragments than for carbon donors in bridged organometallic arms.
    """
    if mol is None:
        return 1.0

    atom = mol.GetAtomWithIdx(donor_idx)
    donor_sym = atom.GetSymbol()
    if donor_sym == 'C':
        weight = 0.85
    elif donor_sym == 'N':
        weight = 1.30
    elif donor_sym == 'O':
        weight = 1.35
    elif donor_sym in {'S', 'P'}:
        weight = 1.15
    else:
        weight = 1.00

    formal_charge = int(atom.GetFormalCharge())
    aromatic = bool(atom.GetIsAromatic())
    try:
        hyb = atom.GetHybridization()
    except Exception:
        hyb = None
    planar_like = aromatic or hyb in {
        Chem.rdchem.HybridizationType.SP,
        Chem.rdchem.HybridizationType.SP2,
    }
    if planar_like and donor_sym in {'N', 'O', 'S'}:
        weight += 0.20
    if donor_sym == 'O' and formal_charge < 0:
        weight += 0.20
    if donor_sym == 'C' and formal_charge < 0:
        weight += 0.10
    return float(weight)


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

# ---------------------------------------------------------------------------
# Preferred CN=5 geometry per metal
# ---------------------------------------------------------------------------
_PREFERRED_CN5_GEOMETRY: Dict[str, str] = {
    # d8 metals: square-pyramidal (SP) preferred
    'Ni': 'SP', 'Pd': 'SP', 'Pt': 'SP',
    # d6 low-spin: SP preferred
    'Co': 'SP', 'Rh': 'SP', 'Ir': 'SP',
    # d0–d5, d7, d10: trigonal-bipyramidal (TBP) preferred
    'Fe': 'TBP', 'Mn': 'TBP', 'Cu': 'TBP',
    'Zn': 'TBP', 'Cd': 'TBP',
    'Ti': 'TBP', 'V': 'TBP', 'Cr': 'TBP',
    'Mo': 'TBP', 'W': 'TBP',
}

# ---------------------------------------------------------------------------
# Preferred CN=6 geometry per metal
# ---------------------------------------------------------------------------
_PREFERRED_CN6_GEOMETRY: Dict[str, str] = {
    # Almost all TMs: octahedral strongly preferred.
    # Trigonal-prismatic only for d0 with specific dithiolene ligands.
    # Default is OH for everything not listed here.
    'Mo': 'OH', 'W': 'OH', 'Re': 'OH',
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


# ---------------------------------------------------------------------------
# Timeout-guarded embedding
# ---------------------------------------------------------------------------
_EMBED_TIMEOUT: float = float(os.environ.get("DELFIN_EMBED_TIMEOUT", "10"))
"""Per-call timeout (seconds) for RDKit EmbedMolecule / stk embedding.

Large metal complexes with highly connected ring systems can cause ETKDG
distance-bounds calculation to hang indefinitely.  This timeout prevents
the converter from blocking forever.  Override via env var if needed.
"""

_OB_ROTOR_TIMEOUT: float = float(os.environ.get("DELFIN_OB_ROTOR_TIMEOUT", "15"))
"""Per-call timeout (seconds) for Open Babel rotor search."""


def _embed_with_timeout(mol, params=None, timeout: Optional[float] = None):
    """Run AllChem.EmbedMolecule with a timeout guard.

    Returns the embed result (0 on success, -1 on failure/timeout).
    The *mol* is modified in-place on success (same as EmbedMolecule).
    """
    if timeout is None:
        timeout = _EMBED_TIMEOUT

    # Fast path: skip timeout machinery for very small molecules
    # BUT always use timeout for highly connected molecules (many ring bonds)
    # which can cause ETKDG to hang even at low atom counts (e.g. borane cages).
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()
    highly_connected = n_bonds > 2 * n_atoms  # cage/cluster topology
    if n_atoms < 40 and not highly_connected:
        if params is not None:
            return AllChem.EmbedMolecule(mol, params)
        return AllChem.EmbedMolecule(mol)

    result = [-1]
    exc_holder = [None]

    def _do_embed():
        try:
            if params is not None:
                result[0] = AllChem.EmbedMolecule(mol, params)
            else:
                result[0] = AllChem.EmbedMolecule(mol)
        except Exception as e:
            exc_holder[0] = e

    t = threading.Thread(target=_do_embed, daemon=True)
    t.start()
    t.join(timeout=timeout)

    if t.is_alive():
        logger.debug(
            "EmbedMolecule timed out after %.1fs for mol with %d atoms",
            timeout, mol.GetNumAtoms(),
        )
        # Thread cannot be killed; it will finish eventually in the background.
        # Return failure so the caller moves on to the next strategy.
        return -1

    if exc_holder[0] is not None:
        raise exc_holder[0]

    return result[0]


def _make_random_embed_params(seed: int = 42):
    """Create minimal embed params that skip ETKDG torsion/knowledge checks.

    These are much faster for complex metal ring systems where the ETKDG
    distance bounds matrix computation can hang.
    """
    params = AllChem.EmbedParameters()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.useBasicKnowledge = False
    params.useExpTorsionAnglePrefs = False
    params.enforceChirality = False
    return params


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
    # stk handles H atoms correctly based on the SMILES notation.
    # stk.BuildingBlock internally runs ETKDG embedding which can hang
    # on complex metal ring systems → guard with timeout.
    if STK_AVAILABLE:
        for smi in smiles_variants:
            try:
                _stk_result = [None]

                def _stk_build(s=smi):
                    try:
                        bb = stk.BuildingBlock(s)
                        _stk_result[0] = bb.to_rdkit_mol()
                    except Exception:
                        _stk_result[0] = None

                _st = threading.Thread(target=_stk_build, daemon=True)
                _st.start()
                _st.join(timeout=_EMBED_TIMEOUT)
                if _st.is_alive():
                    errors.append(f"stk({smi[:30]}...): timed out")
                    continue
                mol = _stk_result[0]
                if mol is not None:
                    # Embed if needed
                    if mol.GetNumConformers() == 0:
                        params = AllChem.ETKDGv3()
                        params.randomSeed = 42
                        params.useRandomCoords = True
                        result = _embed_with_timeout(mol, params)
                        if result != 0:
                            # Fallback: permissive embed (no ETKDG knowledge)
                            result = _embed_with_timeout(
                                mol, _make_random_embed_params(42))
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
                result = _embed_with_timeout(mol, params)
                if result != 0:
                    # Fallback: permissive embed (no ETKDG knowledge)
                    result = _embed_with_timeout(
                        mol, _make_random_embed_params(42))
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

        result = _embed_with_timeout(mol, embed_params)

        if result != 0:
            # Try with different random seeds
            for seed in [123, 456, 789, 1000]:
                embed_params.randomSeed = seed
                result = _embed_with_timeout(mol, embed_params)
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

            # Offset each metal so multi-metal complexes don't overlap
            _metal_offset = 4.0  # Angstrom separation between metal centers
            for metal_rank, mi in enumerate(metal_indices):
                mx = _metal_offset * metal_rank
                local_coords[mi] = (mx, 0.0, 0.0)
                local_placed.add(mi)

                neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(mi).GetNeighbors()]
                n_coord = len(neighbors)
                if n_coord == 4:
                    vectors = override_4coord
                else:
                    vectors = _COORD_VECTORS_BASE.get(n_coord, _COORD_VECTORS_BASE.get(6, []))

                metal_sym = mol.GetAtomWithIdx(mi).GetSymbol()
                # For high-CN systems (>8), use Fibonacci sphere for even distribution
                if n_coord > 8:
                    golden_ratio = (1 + math.sqrt(5)) / 2
                    vectors = []
                    for fi in range(n_coord):
                        theta = math.acos(1 - 2 * (fi + 0.5) / n_coord)
                        phi = 2 * math.pi * fi / golden_ratio
                        vectors.append((math.sin(theta) * math.cos(phi),
                                        math.sin(theta) * math.sin(phi),
                                        math.cos(theta)))
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
                    local_coords[nbr_idx] = (mx + vx, vy, vz)
                    local_placed.add(nbr_idx)

            # BFS to place remaining atoms
            bond_len_default = 1.4
            _bfs_counter = [0]  # mutable counter for unique directions
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
                        # Use Fibonacci sphere direction for unique placement
                        _bfs_counter[0] += 1
                        golden = (1 + math.sqrt(5)) / 2
                        theta = math.acos(1 - 2 * (_bfs_counter[0] % 50 + 0.5) / 50)
                        phi = 2 * math.pi * _bfs_counter[0] / golden
                        dx = math.sin(theta) * math.cos(phi)
                        dy = math.sin(theta) * math.sin(phi)
                        dz = math.cos(theta)
                        mag = 1.0
                    dx = dx/mag * bond_len_default
                    dy = dy/mag * bond_len_default
                    dz = dz/mag * bond_len_default
                    # Clash check: if target position is too close to any placed atom, nudge
                    nx, ny, nz = cx + dx, cy + dy, cz + dz
                    for _attempt in range(5):
                        too_close = False
                        for pi in local_placed:
                            px, py, pz = local_coords[pi]
                            dd = math.sqrt((nx-px)**2 + (ny-py)**2 + (nz-pz)**2)
                            if dd < 0.8:
                                too_close = True
                                break
                        if not too_close:
                            break
                        # Rotate direction by 60 degrees around z
                        cos60, sin60 = 0.5, 0.866
                        dx, dy = dx*cos60 - dy*sin60, dx*sin60 + dy*cos60
                        nx, ny, nz = cx + dx, cy + dy, cz + dz
                    local_coords[nbr_idx] = (nx, ny, nz)
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

            # OB UFF refinement for both — through safe wrapper for
            # universal coordination constraints.
            xyz_tetra = _optimize_xyz_openbabel_safe(xyz_tetra, mol_template=mol)
            xyz_sq = _optimize_xyz_openbabel_safe(xyz_sq, mol_template=mol)

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
            xyz_str = _optimize_xyz_openbabel_safe(xyz_str, mol_template=mol)

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
    deterministic: bool = True,
) -> Tuple[List[str], Optional[str]]:
    """Generate 3D conformers using Open Babel in an Avogadro-like workflow.

    Pipeline (mirrors Avogadro's ``gen3d`` quality):
    1. Fragment-based 3D initialization via ``make3D`` (uses CSD crystal data)
    2. Rotor search for conformer pool generation:
       - deterministic mode: ``SystematicRotorSearch``
       - non-deterministic mode: ``WeightedRotorSearch``
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

        # Run rotor search.  Skip for large molecules where
        # SystematicRotorSearch is combinatorially explosive and holds the
        # GIL, making thread-based timeouts ineffective.
        n_heavy = sum(1 for a in pybel.ob.OBMolAtomIter(ob_mol.OBMol)
                       if a.GetAtomicNum() > 1)
        if n_heavy > 50:
            logger.debug(
                "Skipping OB rotor search for large molecule (%d heavy atoms), "
                "using initial make3D geometry", n_heavy,
            )
        else:
            try:
                if deterministic:
                    ff.SystematicRotorSearch(target)
                else:
                    ff.WeightedRotorSearch(target, max(25, int(rotor_steps)))
                ff.GetConformers(ob_mol.OBMol)
            except Exception as exc_inner:
                logger.debug("Open Babel conformer search failed: %s", exc_inner)

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


def _xyz_to_rdkit_conformer_via_ob_mapping(mol, xyz_delfin: str):
    """Build an RDKit conformer from XYZ by graph-based OB↔RDKit atom mapping.

    This is a fallback for Open Babel conformers where atom order in XYZ does
    not match the RDKit molecule order.  It reconstructs an OB-perceived graph
    from XYZ, then maps atoms to *mol* by element/topology.
    """
    if not (RDKIT_AVAILABLE and OPENBABEL_AVAILABLE):
        return None

    lines = [ln.strip() for ln in xyz_delfin.splitlines() if ln.strip()]
    n_atoms = mol.GetNumAtoms()
    if len(lines) < n_atoms:
        return None

    # Parse XYZ symbols once; coordinates for the mapped conformer are taken
    # from the OB-derived RDKit molecule (which keeps a consistent atom order
    # with its own graph representation).
    xyz_symbols: List[str] = []
    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            return None
        xyz_symbols.append(parts[0])

    rd_symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]
    xyz_counts = Counter(xyz_symbols)
    rd_counts = Counter(rd_symbols)
    if any(xyz_counts.get(sym, 0) < cnt for sym, cnt in rd_counts.items()):
        return None

    # Build standard XYZ and let Open Babel perceive connectivity.
    try:
        std_xyz = f"{n_atoms}\n\n" + "\n".join(
            f"{ln.split()[0]}  {ln.split()[1]}  {ln.split()[2]}  {ln.split()[3]}"
            for ln in lines
        ) + "\n"
        ob_py = pybel.readstring("xyz", std_xyz)
        conv = pybel.ob.OBConversion()
        if not conv.SetOutFormat("mol"):
            return None
        mol_block = conv.WriteString(ob_py.OBMol)
        if not mol_block:
            return None
        ob_rd = Chem.MolFromMolBlock(mol_block, sanitize=False, removeHs=False)
        if ob_rd is None or ob_rd.GetNumAtoms() < n_atoms:
            return None
        try:
            ob_rd.UpdatePropertyCache(strict=False)
        except Exception:
            pass
    except Exception:
        return None

    ob_symbols = [ob_rd.GetAtomWithIdx(i).GetSymbol() for i in range(ob_rd.GetNumAtoms())]
    ob_counts = Counter(ob_symbols)
    if any(ob_counts.get(sym, 0) < cnt for sym, cnt in rd_counts.items()):
        return None

    def _build_topology_submol(src_mol, keep_indices: List[int]):
        """Create a single-bond topology-only submol; return (submol, old_order)."""
        rw = Chem.RWMol()
        old_order: List[int] = []
        old_to_new: Dict[int, int] = {}
        for old_idx in keep_indices:
            at = src_mol.GetAtomWithIdx(old_idx)
            nat = Chem.Atom(int(at.GetAtomicNum()))
            nat.SetFormalCharge(0)
            nat.SetNoImplicit(True)
            old_to_new[old_idx] = rw.AddAtom(nat)
            old_order.append(old_idx)
        keep_set = set(keep_indices)
        for bond in src_mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            if i in keep_set and j in keep_set:
                ni = old_to_new[i]
                nj = old_to_new[j]
                if rw.GetBondBetweenAtoms(ni, nj) is None:
                    rw.AddBond(ni, nj, Chem.BondType.SINGLE)
        sub = rw.GetMol()
        try:
            sub.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        return sub, old_order

    # 1) Match non-metal heavy-atom skeleton (most robust anchor).
    rd_core = [
        i for i in range(n_atoms)
        if mol.GetAtomWithIdx(i).GetAtomicNum() > 1
        and mol.GetAtomWithIdx(i).GetSymbol() not in _METAL_SET
    ]
    ob_core = [
        i for i in range(ob_rd.GetNumAtoms())
        if ob_rd.GetAtomWithIdx(i).GetAtomicNum() > 1
        and ob_rd.GetAtomWithIdx(i).GetSymbol() not in _METAL_SET
    ]
    if len(rd_core) != len(ob_core):
        return None
    if Counter(mol.GetAtomWithIdx(i).GetSymbol() for i in rd_core) != Counter(
        ob_rd.GetAtomWithIdx(i).GetSymbol() for i in ob_core
    ):
        return None

    rd_sub, rd_order = _build_topology_submol(mol, rd_core)
    ob_sub, ob_order = _build_topology_submol(ob_rd, ob_core)

    matches = ()
    if rd_sub.GetNumAtoms() > 0:
        try:
            matches = ob_sub.GetSubstructMatches(
                rd_sub, uniquify=False, useChirality=False, maxMatches=256
            )
        except Exception:
            matches = ()
        if not matches:
            return None

    # Pick the match with best heavy-neighbor-degree consistency.
    def _core_deg(src_mol, idx: int) -> int:
        atom = src_mol.GetAtomWithIdx(idx)
        return sum(
            1
            for nb in atom.GetNeighbors()
            if nb.GetAtomicNum() > 1 and nb.GetSymbol() not in _METAL_SET
        )

    best_match = ()
    if matches:
        best_score = float("inf")
        for cand in matches:
            score = 0.0
            for q_new, t_new in enumerate(cand):
                rd_old = rd_order[q_new]
                ob_old = ob_order[t_new]
                score += abs(_core_deg(mol, rd_old) - _core_deg(ob_rd, ob_old))
            if score < best_score:
                best_score = score
                best_match = cand
        if not best_match:
            best_match = matches[0]

    mapping: Dict[int, int] = {}
    if rd_sub.GetNumAtoms() > 0:
        for q_new, t_new in enumerate(best_match):
            mapping[rd_order[q_new]] = ob_order[t_new]

    used_ob = set(mapping.values())

    # 2) Match metal atoms by element symbol.
    metal_symbols = sorted({
        mol.GetAtomWithIdx(i).GetSymbol()
        for i in range(n_atoms)
        if mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET
    })
    for sym in metal_symbols:
        rd_m = sorted(
            i for i in range(n_atoms)
            if mol.GetAtomWithIdx(i).GetSymbol() == sym
        )
        ob_m = sorted(
            i for i in range(ob_rd.GetNumAtoms())
            if ob_rd.GetAtomWithIdx(i).GetSymbol() == sym and i not in used_ob
        )
        if len(ob_m) < len(rd_m):
            return None
        for ri, oi in zip(rd_m, ob_m[:len(rd_m)]):
            mapping[ri] = oi
            used_ob.add(oi)

    # 3) Map H atoms via the mapped heavy-atom neighbour if possible.
    ob_conf = ob_rd.GetConformer()

    def _dist_sq(i: int, j: int) -> float:
        pi = ob_conf.GetAtomPosition(i)
        pj = ob_conf.GetAtomPosition(j)
        dx = pi.x - pj.x
        dy = pi.y - pj.y
        dz = pi.z - pj.z
        return dx * dx + dy * dy + dz * dz

    rd_h = sorted(
        i for i in range(n_atoms)
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 1
    )
    ob_h_all = {
        i for i in range(ob_rd.GetNumAtoms())
        if ob_rd.GetAtomWithIdx(i).GetAtomicNum() == 1 and i not in used_ob
    }
    for h_idx in rd_h:
        if h_idx in mapping:
            continue
        h_atom = mol.GetAtomWithIdx(h_idx)
        anchor_rd = None
        for nb in h_atom.GetNeighbors():
            if nb.GetAtomicNum() > 1:
                anchor_rd = nb.GetIdx()
                break

        candidates: List[int] = []
        if anchor_rd is not None and anchor_rd in mapping:
            anchor_ob = mapping[anchor_rd]
            anchor_atom_ob = ob_rd.GetAtomWithIdx(anchor_ob)
            candidates = [
                nb.GetIdx()
                for nb in anchor_atom_ob.GetNeighbors()
                if nb.GetAtomicNum() == 1 and nb.GetIdx() in ob_h_all
            ]
            if not candidates:
                candidates = sorted(ob_h_all, key=lambda j: _dist_sq(anchor_ob, j))
        else:
            candidates = sorted(ob_h_all)

        if not candidates:
            return None
        chosen = candidates[0]
        mapping[h_idx] = chosen
        used_ob.add(chosen)
        ob_h_all.discard(chosen)

    # 4) Map any remaining atoms by symbol (rare fallback).
    remaining_rd = [i for i in range(n_atoms) if i not in mapping]
    remaining_ob = [i for i in range(ob_rd.GetNumAtoms()) if i not in used_ob]
    if remaining_rd:
        by_sym_rd: Dict[str, List[int]] = {}
        by_sym_ob: Dict[str, List[int]] = {}
        for i in remaining_rd:
            by_sym_rd.setdefault(mol.GetAtomWithIdx(i).GetSymbol(), []).append(i)
        for i in remaining_ob:
            by_sym_ob.setdefault(ob_rd.GetAtomWithIdx(i).GetSymbol(), []).append(i)
        if set(by_sym_rd) != set(by_sym_ob):
            return None
        for sym in sorted(by_sym_rd):
            r_list = sorted(by_sym_rd[sym])
            o_list = sorted(by_sym_ob[sym])
            if len(o_list) < len(r_list):
                return None
            for ri, oi in zip(r_list, o_list[:len(r_list)]):
                mapping[ri] = oi
                used_ob.add(oi)

    if len(mapping) != n_atoms:
        return None

    conf = Chem.Conformer(n_atoms)
    for rd_idx in range(n_atoms):
        ob_idx = mapping.get(rd_idx)
        if ob_idx is None:
            return None
        if mol.GetAtomWithIdx(rd_idx).GetSymbol() != ob_rd.GetAtomWithIdx(ob_idx).GetSymbol():
            return None
        pos = ob_conf.GetAtomPosition(ob_idx)
        conf.SetAtomPosition(rd_idx, (pos.x, pos.y, pos.z))
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
    n_mapped_fallback = 0
    for xyz_text in xyz_blocks:
        conf = _xyz_to_rdkit_conformer(mol, xyz_text)
        used_fallback = False
        if conf is None:
            conf = _xyz_to_rdkit_conformer_via_ob_mapping(mol, xyz_text)
            used_fallback = conf is not None
        if conf is None:
            continue
        conf_ids.append(int(mol.AddConformer(conf, assignId=True)))
        if used_fallback:
            n_mapped_fallback += 1
    if n_mapped_fallback:
        logger.debug(
            "OB conformer mapping fallback succeeded for %d block(s).",
            n_mapped_fallback,
        )
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


def _fix_hapto_donor_h(mol):
    """Adjust H on hapto donor carbons so each C has at most 4 bonding partners.

    Simple rule: desired_h = max(0, 4 - number_of_non_H_neighbors).
    Metal neighbors ARE counted (they occupy a coordination site).
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


def contains_metal(smiles: str) -> bool:
    """Return True if SMILES likely contains a metal atom."""
    for metal in _METALS:
        if re.search(rf'\[{metal}[+\-\d\]@]', smiles, re.IGNORECASE):
            return True
        if re.search(rf'\[{metal}\]', smiles, re.IGNORECASE):
            return True
    return False


def _hapto_approx_enabled(flag: Optional[bool] = None) -> bool:
    """Return True when the experimental hapto approximation is enabled."""
    if flag is not None:
        return bool(flag)
    env = os.environ.get("DELFIN_HAPTO_APPROX", "").strip().lower()
    return env in {"1", "true", "yes", "on"}


def _find_hapto_groups(mol) -> List[Tuple[int, List[int]]]:
    """Detect likely eta/hapto groups as contiguous metal-bound carbon sets."""
    if not RDKIT_AVAILABLE or mol is None:
        return []

    groups: List[Tuple[int, List[int]]] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        all_neighbors = list(atom.GetNeighbors())
        c_neighbors = [n.GetIdx() for n in all_neighbors if n.GetAtomicNum() == 6]
        if len(c_neighbors) < 2:
            continue
        # Skip cluster compounds (borane/carborane cages): if >50% of
        # metal neighbours are B atoms, the C atoms are cage vertices,
        # not part of a cyclopentadienyl-type hapto ligand.
        b_count = sum(1 for n in all_neighbors if n.GetAtomicNum() == 5)
        if b_count > len(all_neighbors) / 2:
            continue

        c_set = set(c_neighbors)
        seen: set = set()
        for start in c_neighbors:
            if start in seen:
                continue
            comp: List[int] = []
            stack = [start]
            seen.add(start)
            while stack:
                cur = stack.pop()
                comp.append(cur)
                cur_atom = mol.GetAtomWithIdx(cur)
                for nbr in cur_atom.GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni not in c_set or ni in seen:
                        continue
                    seen.add(ni)
                    stack.append(ni)

            if len(comp) < 2:
                continue
            # Avoid false positives for ordinary C,C chelation: classify as hapto
            # when the contiguous donor block is ring-like or has >=3 atoms.
            ring_like = any(mol.GetAtomWithIdx(i).IsInRing() for i in comp)
            if len(comp) >= 3 or ring_like:
                groups.append((metal_idx, sorted(comp)))

    return groups


def _probe_hapto_groups_from_smiles(smiles: str) -> List[Tuple[int, List[int]]]:
    """Parse SMILES quickly and report likely hapto groups."""
    if not RDKIT_AVAILABLE or not contains_metal(smiles):
        return []

    mol, _note = mol_from_smiles_rdkit(smiles, allow_metal=True)
    if mol is None:
        try:
            p = Chem.SmilesParserParams()
            p.sanitize = False
            p.removeHs = False
            p.strictParsing = False
            mol = Chem.MolFromSmiles(smiles, p)
        except Exception:
            try:
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
            except Exception:
                mol = None
    return _find_hapto_groups(mol)


def _hapto_failfast_error(hapto_groups: List[Tuple[int, List[int]]]) -> str:
    """Build a concise user-facing error for unsupported hapto coordination."""
    if not hapto_groups:
        return (
            "Hapto (eta) coordination detected. "
            "Enable DELFIN_HAPTO_APPROX=1 for experimental approximation mode."
        )
    max_group = max(len(g[1]) for g in hapto_groups)
    return (
        "Hapto (eta) coordination detected "
        f"({len(hapto_groups)} group(s), max eta~{max_group}). "
        "Standard conversion does not support this reliably. "
        "Set DELFIN_HAPTO_APPROX=1 to enable experimental approximation mode."
    )


def _select_multihapto_anchors(
    mol,
    metal_idx: int,
    metal_groups: List[List[int]],
) -> List[int]:
    """Choose one anchor per hapto group, maximally spread around metal_idx.

    For metals with 2+ groups: pick anchors from opposite ends of the graph
    to avoid two anchors being direct neighbours.
    """
    if len(metal_groups) <= 1:
        # Single group: just pick best anchor
        grp = metal_groups[0]
        grp_set = set(grp)

        def _akey(idx: int) -> Tuple[int, int, int, int]:
            a = mol.GetAtomWithIdx(idx)
            aromatic = 1 if a.GetIsAromatic() else 0
            in_ring = 1 if a.IsInRing() else 0
            c_in_grp = sum(
                1 for n in a.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() in grp_set
            )
            return (aromatic, in_ring, c_in_grp, -idx)
        return [max(grp, key=_akey)]

    # For each group, compute a "position score" = average graph distance
    # from the metal through the group members.  Pick anchors that are
    # maximally spread by assigning them greedily.
    anchors: List[int] = []
    used_anchors: set = set()
    for grp in metal_groups:
        grp_set = set(grp)

        def _anchor_score(idx: int) -> Tuple[float, int, int, int, int]:
            a = mol.GetAtomWithIdx(idx)
            aromatic = 1 if a.GetIsAromatic() else 0
            in_ring = 1 if a.IsInRing() else 0
            c_in_grp = sum(
                1 for n in a.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() in grp_set
            )
            # Penalty for being adjacent to an already-chosen anchor
            adj_penalty = 0
            for n in a.GetNeighbors():
                if n.GetIdx() in used_anchors:
                    adj_penalty = -2
                    break
            return (adj_penalty, aromatic, in_ring, c_in_grp, -idx)

        anchor = max(grp, key=_anchor_score)
        anchors.append(anchor)
        used_anchors.add(anchor)

    return anchors


def _apply_hapto_approximation(
    mol,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
):
    """Experimental eta->anchor approximation for embedding/optimization.

    For each contiguous metal-bound carbon block, keep one representative
    metal-carbon bond as a normal bond and convert the other M-C contacts to
    dative bonds (C->M). For metals with multiple hapto groups, anchors are
    chosen to be maximally spread to improve embedding convergence.
    """
    if not RDKIT_AVAILABLE or mol is None:
        return mol, 0

    groups = hapto_groups if hapto_groups is not None else _find_hapto_groups(mol)
    if not groups:
        return mol, 0

    # Group hapto groups by metal for coordinated anchor selection
    by_metal: Dict[int, List[List[int]]] = {}
    by_metal_order: Dict[int, List[int]] = {}  # track group indices
    for gi, (metal_idx, grp) in enumerate(groups):
        if len(grp) < 2:
            continue
        by_metal.setdefault(metal_idx, []).append(grp)
        by_metal_order.setdefault(metal_idx, []).append(gi)

    # Compute coordinated anchors per metal
    anchor_map: Dict[int, int] = {}  # group_index -> anchor_atom_idx
    for metal_idx, metal_groups in by_metal.items():
        anchors = _select_multihapto_anchors(mol, metal_idx, metal_groups)
        for i, anchor in enumerate(anchors):
            gi = by_metal_order[metal_idx][i]
            anchor_map[gi] = anchor

    rw = Chem.RWMol(mol)
    converted = 0
    for gi, (metal_idx, grp) in enumerate(groups):
        if len(grp) < 2:
            continue
        anchor = anchor_map.get(gi)
        if anchor is None:
            grp_set = set(grp)

            def _anchor_key(idx: int) -> Tuple[int, int, int, int]:
                a = rw.GetAtomWithIdx(idx)
                aromatic = 1 if a.GetIsAromatic() else 0
                in_ring = 1 if a.IsInRing() else 0
                c_in_grp = sum(
                    1 for n in a.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() in grp_set
                )
                return (aromatic, in_ring, c_in_grp, -idx)
            anchor = max(grp, key=_anchor_key)

        for c_idx in grp:
            c_atom = rw.GetAtomWithIdx(c_idx)
            c_atom.SetNoImplicit(False)
            if c_idx == anchor:
                continue
            bond = rw.GetBondBetweenAtoms(metal_idx, c_idx)
            if bond is None:
                continue
            if bond.GetBondType() == Chem.BondType.DATIVE:
                continue
            rw.RemoveBond(metal_idx, c_idx)
            rw.AddBond(c_idx, metal_idx, Chem.BondType.DATIVE)
            converted += 1

    out = rw.GetMol()
    try:
        out.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return out, converted


def _copy_fragment_atom(atom, atom_map_num: int):
    """Create a detached copy of an atom with its valence/H metadata."""
    new_atom = Chem.Atom(atom.GetAtomicNum())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetNumExplicitHs(atom.GetNumExplicitHs())
    new_atom.SetNoImplicit(atom.GetNoImplicit())
    new_atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons())
    new_atom.SetIsAromatic(atom.GetIsAromatic())
    new_atom.SetChiralTag(atom.GetChiralTag())
    new_atom.SetAtomMapNum(atom_map_num)
    try:
        new_atom.SetHybridization(atom.GetHybridization())
    except Exception:
        pass
    try:
        new_atom.SetIsotope(atom.GetIsotope())
    except Exception:
        pass
    return new_atom


def _extract_fragment_mol(mol, atom_indices: List[int]):
    """Copy a connected metal-free fragment out of ``mol``."""
    if not RDKIT_AVAILABLE or mol is None or not atom_indices:
        return None, {}, {}

    atom_set = set(atom_indices)
    rw = Chem.RWMol()
    original_to_fragment: Dict[int, int] = {}
    fragment_to_original: Dict[int, int] = {}

    for old_idx in atom_indices:
        old_atom = mol.GetAtomWithIdx(old_idx)
        new_idx = rw.AddAtom(_copy_fragment_atom(old_atom, old_idx + 1))
        original_to_fragment[old_idx] = new_idx
        fragment_to_original[new_idx] = old_idx

    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx not in atom_set or end_idx not in atom_set:
            continue
        new_begin = original_to_fragment[begin_idx]
        new_end = original_to_fragment[end_idx]
        rw.AddBond(new_begin, new_end, bond.GetBondType())
        new_bond = rw.GetBondBetweenAtoms(new_begin, new_end)
        if new_bond is None:
            continue
        new_bond.SetIsAromatic(bond.GetIsAromatic())
        new_bond.SetIsConjugated(bond.GetIsConjugated())
        try:
            new_bond.SetBondDir(bond.GetBondDir())
            new_bond.SetStereo(bond.GetStereo())
            stereo_atoms = list(bond.GetStereoAtoms())
            if len(stereo_atoms) == 2:
                new_bond.SetStereoAtoms(
                    original_to_fragment[stereo_atoms[0]],
                    original_to_fragment[stereo_atoms[1]],
                )
        except Exception:
            pass

    frag = rw.GetMol()
    try:
        frag.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    try:
        Chem.FastFindRings(frag)
    except Exception:
        pass
    return frag, original_to_fragment, fragment_to_original


def _choose_fragment_anchor_atoms(
    mol,
    frag_atoms: List[int],
    donor_indices: set,
    hapto_groups: List[Tuple[int, List[int]]],
    hapto_atom_to_group: Dict[int, int],
) -> Tuple[List[int], List[int], List[int]]:
    """Choose chemically meaningful anchors for rigid fragment alignment."""
    frag_set = set(frag_atoms)
    frag_donors = sorted(idx for idx in frag_atoms if idx in donor_indices)
    frag_group_ids = sorted({hapto_atom_to_group[idx] for idx in frag_atoms if idx in hapto_atom_to_group})

    anchors: List[int] = []
    for gid in frag_group_ids:
        for atom_idx in hapto_groups[gid][1]:
            if atom_idx in frag_set and atom_idx not in anchors:
                anchors.append(atom_idx)

    for atom_idx in frag_donors:
        if atom_idx not in anchors:
            anchors.append(atom_idx)

    seed_atoms = list(anchors) if anchors else list(frag_donors) if frag_donors else list(frag_atoms)
    for atom_idx in seed_atoms:
        for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in frag_set and nbr_idx not in anchors and nbr.GetAtomicNum() > 1:
                anchors.append(nbr_idx)
                if len(anchors) >= 6:
                    break
        if len(anchors) >= 6:
            break

    if len(anchors) < 3:
        for atom_idx in frag_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom_idx in anchors:
                continue
            anchors.append(atom_idx)
            if len(anchors) >= 6:
                break

    return frag_donors, frag_group_ids, anchors


def _decompose_hapto_complex(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> Optional[_HybridHaptoDecomposition]:
    """Split a metal complex into ligand fragments while keeping hapto metadata."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return None

    metal_indices = sorted(
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in _METAL_SET
    )
    if not metal_indices:
        return None

    metal_set = set(metal_indices)
    donor_indices: set = set()
    bonds_to_remove: List[Tuple[int, int]] = []
    hapto_atom_to_group: Dict[int, int] = {}
    for group_idx, (_metal_idx, group_atoms) in enumerate(hapto_groups):
        for atom_idx in group_atoms:
            hapto_atom_to_group[atom_idx] = group_idx

    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in metal_set:
            donor_indices.add(end_idx)
            bonds_to_remove.append((begin_idx, end_idx))
        elif end_idx in metal_set:
            donor_indices.add(begin_idx)
            bonds_to_remove.append((begin_idx, end_idx))

    if not bonds_to_remove:
        return None

    rw = Chem.RWMol(mol)
    for begin_idx, end_idx in bonds_to_remove:
        if rw.GetBondBetweenAtoms(begin_idx, end_idx) is not None:
            rw.RemoveBond(begin_idx, end_idx)
    for metal_idx in sorted(metal_indices, reverse=True):
        rw.RemoveAtom(metal_idx)

    metal_free = rw.GetMol()
    try:
        metal_free.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    old_to_new: Dict[int, int] = {}
    new_to_old: Dict[int, int] = {}
    removed_count = 0
    for old_idx in range(mol.GetNumAtoms()):
        if old_idx in metal_set:
            removed_count += 1
            continue
        new_idx = old_idx - removed_count
        old_to_new[old_idx] = new_idx
        new_to_old[new_idx] = old_idx

    try:
        frags = list(Chem.GetMolFrags(metal_free, asMols=False, sanitizeFrags=False))
    except Exception:
        frags = []

    fragments: List[_HybridHaptoFragment] = []
    for frag_new in frags:
        frag_atoms = sorted(new_to_old[idx] for idx in frag_new if idx in new_to_old)
        if not frag_atoms:
            continue
        fragment_mol, original_to_fragment, fragment_to_original = _extract_fragment_mol(mol, frag_atoms)
        if fragment_mol is None:
            continue
        frag_donors, frag_group_ids, anchors = _choose_fragment_anchor_atoms(
            mol,
            frag_atoms,
            donor_indices,
            hapto_groups,
            hapto_atom_to_group,
        )
        metal_neighbors = sorted({
            nbr.GetIdx()
            for atom_idx in frag_donors
            for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
            if nbr.GetSymbol() in _METAL_SET
        })
        bridging_donors = sorted(
            atom_idx
            for atom_idx in frag_donors
            if sum(
                1
                for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
                if nbr.GetSymbol() in _METAL_SET
            ) >= 2
        )
        use_scaffold_only = len(metal_neighbors) >= 2 or bool(bridging_donors)
        fragments.append(
            _HybridHaptoFragment(
                atom_indices=frag_atoms,
                donor_atom_indices=frag_donors,
                metal_neighbor_indices=metal_neighbors,
                bridging_donor_indices=bridging_donors,
                hapto_group_ids=frag_group_ids,
                anchor_atom_indices=anchors,
                original_to_fragment=original_to_fragment,
                fragment_to_original=fragment_to_original,
                fragment_mol=fragment_mol,
                use_scaffold_only=use_scaffold_only,
            )
        )

    if not fragments:
        return None

    fragments.sort(
        key=lambda frag: (
            frag.use_scaffold_only,
            -len(frag.hapto_group_ids),
            -len(frag.anchor_atom_indices),
            -len(frag.donor_atom_indices),
            -len(frag.atom_indices),
        )
    )
    return _HybridHaptoDecomposition(
        metal_indices=metal_indices,
        donor_atom_indices=sorted(donor_indices),
        hapto_groups=hapto_groups,
        fragments=fragments,
    )


def _embed_hybrid_fragment(fragment_mol):
    """Embed a detached ligand fragment while preserving its local chemistry."""
    if not RDKIT_AVAILABLE or fragment_mol is None:
        return None

    work = Chem.Mol(fragment_mol)
    work.RemoveAllConformers()

    for atom in work.GetAtoms():
        if atom.GetDegree() < 3 and atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        if (
            atom.GetSymbol() == 'O'
            and atom.GetFormalCharge() == -1
            and any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())
        ):
            atom.SetFormalCharge(0)

    try:
        work.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    seeds = [42, 1337, 7]
    embedded = False
    for seed in seeds:
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.useRandomCoords = True
            params.enforceChirality = False
            result = AllChem.EmbedMolecule(work, params)
        except Exception:
            result = -1
        if result == 0:
            embedded = True
            break

    if not embedded:
        try:
            result = AllChem.EmbedMolecule(work, useRandomCoords=True, randomSeed=42)
        except Exception:
            result = -1
        embedded = result == 0

    if not embedded:
        return None

    try:
        AllChem.UFFOptimizeMolecule(work, maxIters=200)
    except Exception:
        try:
            AllChem.MMFFOptimizeMolecule(work, maxIters=200)
        except Exception:
            pass

    return work


def _detect_primary_organometal_module(
    decomposition: Optional[_HybridHaptoDecomposition],
    hapto_groups: List[Tuple[int, List[int]]],
) -> Optional[_PrimaryOrganometalModule]:
    """Detect a primary hapto center with buildable non-hapto donor blocks."""
    if decomposition is None or not hapto_groups:
        return None

    hapto_metals = sorted({metal_idx for metal_idx, _grp in hapto_groups})
    if len(hapto_metals) != 1:
        return None

    metal_idx = hapto_metals[0]
    hapto_group_ids = [
        group_idx for group_idx, (group_metal_idx, _group_atoms) in enumerate(hapto_groups)
        if group_metal_idx == metal_idx
    ]
    if not hapto_group_ids:
        return None

    correlated_fragment_indices: List[int] = []
    terminal_fragment_indices: List[int] = []
    donor_atom_indices: List[int] = []
    for frag_idx, fragment in enumerate(decomposition.fragments):
        if metal_idx not in fragment.metal_neighbor_indices:
            continue
        if any(group_id in hapto_group_ids for group_id in fragment.hapto_group_ids):
            continue
        if fragment.use_scaffold_only:
            continue
        if any(neighbor_idx != metal_idx for neighbor_idx in fragment.metal_neighbor_indices):
            continue
        frag_donors = [
            donor_idx for donor_idx in fragment.donor_atom_indices
            if donor_idx not in donor_atom_indices
        ]
        if not frag_donors:
            continue
        donor_atom_indices.extend(frag_donors)
        if len(fragment.donor_atom_indices) >= 2:
            correlated_fragment_indices.append(frag_idx)
        elif len(fragment.donor_atom_indices) == 1:
            terminal_fragment_indices.append(frag_idx)

    if not correlated_fragment_indices and not terminal_fragment_indices:
        return None

    return _PrimaryOrganometalModule(
        metal_idx=metal_idx,
        hapto_group_ids=hapto_group_ids,
        donor_atom_indices=sorted(donor_atom_indices),
        correlated_fragment_indices=correlated_fragment_indices,
        terminal_fragment_indices=terminal_fragment_indices,
    )


def _primary_organometal_module_quality_ok(
    mol,
    decomposition: Optional[_HybridHaptoDecomposition],
    module: Optional[_PrimaryOrganometalModule],
) -> bool:
    """Cheap quality gate for the primary-metal organometal module path."""
    if not RDKIT_AVAILABLE or mol is None or decomposition is None or module is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    metal_idx = int(module.metal_idx)
    metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)

    donor_set = set(module.donor_atom_indices)
    for donor_idx in sorted(donor_set):
        donor_pos = np.array(conf.GetAtomPosition(donor_idx), dtype=float)
        target_len = float(_get_ml_bond_length(metal_sym, mol.GetAtomWithIdx(donor_idx).GetSymbol()))
        dist = float(np.linalg.norm(donor_pos - metal_pos))
        max_err = 0.70 if mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'C' else 0.60
        if abs(dist - target_len) > max_err:
            return False

    inspected_atoms: set = set()
    for frag_idx in module.correlated_fragment_indices + module.terminal_fragment_indices:
        if not (0 <= int(frag_idx) < len(decomposition.fragments)):
            continue
        fragment = decomposition.fragments[int(frag_idx)]
        for atom_idx in fragment.atom_indices:
            if atom_idx in donor_set or atom_idx in inspected_atoms:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            inspected_atoms.add(atom_idx)
            atom_pos = np.array(conf.GetAtomPosition(atom_idx), dtype=float)
            dist = float(np.linalg.norm(atom_pos - metal_pos))
            min_allowed = max(
                1.65,
                0.90 * float(_get_ml_bond_length(metal_sym, atom.GetSymbol())),
            )
            if atom.GetIsAromatic():
                min_allowed = max(min_allowed, 1.85)
            elif atom.GetSymbol() in {'C', 'N', 'O'}:
                min_allowed = max(min_allowed, 1.75)
            if dist < min_allowed - 0.08:
                return False

    return True


def _append_hapto_preview_xyz(
    mol,
    preview_store: Optional[List[Tuple[str, str]]],
    *,
    label: str,
) -> None:
    """Append the current XYZ of a preview candidate if available."""
    if preview_store is None or mol is None:
        return
    try:
        preview_xyz = _mol_to_xyz(mol)
    except Exception:
        return
    preview_entry = (preview_xyz, label)
    if preview_entry not in preview_store:
        preview_store.append(preview_entry)


def _store_hapto_preview_candidate(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
    preview_store: Optional[List[Tuple[str, str]]],
    *,
    label: str,
) -> None:
    """Store a preview XYZ for a hapto scaffold candidate if available."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups or preview_store is None:
        return
    try:
        preview_mol = Chem.Mol(mol)
        if not _build_hapto_scaffold(preview_mol, hapto_groups):
            return
        try:
            _correct_hapto_geometry(preview_mol, 0, hapto_groups)
        except Exception:
            pass
        try:
            _propagate_non_hapto_atoms(
                preview_mol,
                0,
                hapto_groups,
                extra_fixed_indices=_hapto_primary_donor_indices(preview_mol, hapto_groups),
            )
        except Exception:
            pass
        _append_hapto_preview_xyz(
            preview_mol,
            preview_store,
            label=label,
        )
    except Exception:
        pass


def _rotation_matrix_from_vectors(vec_a, vec_b):
    """Return a rotation matrix that maps ``vec_a`` onto ``vec_b``."""
    import numpy as np

    a = np.asarray(vec_a, dtype=float)
    b = np.asarray(vec_b, dtype=float)
    norm_a = float(np.linalg.norm(a))
    norm_b = float(np.linalg.norm(b))
    if norm_a < 1e-12 or norm_b < 1e-12:
        return np.eye(3)

    a /= norm_a
    b /= norm_b
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = float(np.linalg.norm(v))

    if s < 1e-12:
        if c > 0.0:
            return np.eye(3)
        axis = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.8:
            axis = np.array([0.0, 1.0, 0.0])
        v = np.cross(a, axis)
        v /= max(float(np.linalg.norm(v)), 1e-12)
        vx = np.array([
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ])
        return np.eye(3) + 2.0 * vx @ vx

    vx = np.array([
        [0.0, -v[2], v[1]],
        [v[2], 0.0, -v[0]],
        [-v[1], v[0], 0.0],
    ])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def _align_hybrid_fragment_onto_scaffold(
    scaffold_mol,
    fragment: _HybridHaptoFragment,
    embedded_fragment,
) -> bool:
    """Rigidly align an embedded fragment onto scaffold anchor coordinates."""
    if not RDKIT_AVAILABLE or scaffold_mol is None or embedded_fragment is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        scaffold_conf = scaffold_mol.GetConformer()
        fragment_conf = embedded_fragment.GetConformer()
    except Exception:
        return False

    anchor_pairs: List[Tuple[int, int]] = []
    for original_idx in fragment.anchor_atom_indices:
        fragment_idx = fragment.original_to_fragment.get(original_idx)
        if fragment_idx is None:
            continue
        scaffold_pos = scaffold_conf.GetAtomPosition(original_idx)
        target_vec = np.array([scaffold_pos.x, scaffold_pos.y, scaffold_pos.z], dtype=float)
        if not np.all(np.isfinite(target_vec)):
            continue
        anchor_pairs.append((fragment_idx, original_idx))

    if not anchor_pairs:
        return False

    source = []
    target = []
    for fragment_idx, original_idx in anchor_pairs:
        frag_pos = fragment_conf.GetAtomPosition(fragment_idx)
        source.append([frag_pos.x, frag_pos.y, frag_pos.z])
        scaf_pos = scaffold_conf.GetAtomPosition(original_idx)
        target.append([scaf_pos.x, scaf_pos.y, scaf_pos.z])

    src = np.asarray(source, dtype=float)
    tgt = np.asarray(target, dtype=float)
    all_coords = np.asarray(
        [
            [
                fragment_conf.GetAtomPosition(i).x,
                fragment_conf.GetAtomPosition(i).y,
                fragment_conf.GetAtomPosition(i).z,
            ]
            for i in range(embedded_fragment.GetNumAtoms())
        ],
        dtype=float,
    )

    if len(anchor_pairs) == 1:
        delta = tgt[0] - src[0]
        aligned = all_coords + delta
    elif len(anchor_pairs) == 2:
        rot = _rotation_matrix_from_vectors(src[1] - src[0], tgt[1] - tgt[0])
        aligned = (all_coords - src[0]) @ rot.T + tgt[0]
    else:
        src_centroid = src.mean(axis=0)
        tgt_centroid = tgt.mean(axis=0)
        covariance = (src - src_centroid).T @ (tgt - tgt_centroid)
        try:
            u, _s, vt = np.linalg.svd(covariance)
        except Exception:
            return False
        rot = vt.T @ u.T
        if float(np.linalg.det(rot)) < 0.0:
            vt[-1, :] *= -1.0
            rot = vt.T @ u.T
        aligned = (all_coords - src_centroid) @ rot.T + tgt_centroid

    # Keep scaffold-defining anchors exact and let the subsequent hapto-only
    # local relaxation repair nearby internal strain.
    for pair_idx, (fragment_idx, _original_idx) in enumerate(anchor_pairs):
        aligned[fragment_idx] = tgt[pair_idx]

    for fragment_idx, original_idx in fragment.fragment_to_original.items():
        scaffold_conf.SetAtomPosition(
            original_idx,
            Point3D(
                float(aligned[fragment_idx, 0]),
                float(aligned[fragment_idx, 1]),
                float(aligned[fragment_idx, 2]),
            ),
        )
    return True


def _align_hybrid_fragment_to_targets(
    scaffold_mol,
    fragment: _HybridHaptoFragment,
    embedded_fragment,
    target_atom_indices: List[int],
    target_positions: Dict[int, object],
    reference_center=None,
    metal_idx: Optional[int] = None,
    exact_target_count: int = 0,
) -> bool:
    """Rigidly align a fragment to explicit target coordinates."""
    if not RDKIT_AVAILABLE or scaffold_mol is None or embedded_fragment is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        scaffold_conf = scaffold_mol.GetConformer()
        fragment_conf = embedded_fragment.GetConformer()
    except Exception:
        return False

    anchor_pairs: List[Tuple[int, int]] = []
    seen_original: set = set()
    for original_idx in target_atom_indices:
        if original_idx in seen_original:
            continue
        seen_original.add(original_idx)
        fragment_idx = fragment.original_to_fragment.get(original_idx)
        if fragment_idx is None:
            continue
        target_value = target_positions.get(original_idx)
        if target_value is None:
            continue
        target_vec = np.asarray(target_value, dtype=float)
        if target_vec.shape != (3,) or not np.all(np.isfinite(target_vec)):
            continue
        anchor_pairs.append((fragment_idx, original_idx))

    if not anchor_pairs:
        return False

    src = np.asarray(
        [
            [
                fragment_conf.GetAtomPosition(fragment_idx).x,
                fragment_conf.GetAtomPosition(fragment_idx).y,
                fragment_conf.GetAtomPosition(fragment_idx).z,
            ]
            for fragment_idx, _original_idx in anchor_pairs
        ],
        dtype=float,
    )
    tgt = np.asarray(
        [np.asarray(target_positions[original_idx], dtype=float) for _fragment_idx, original_idx in anchor_pairs],
        dtype=float,
    )
    all_coords = np.asarray(
        [
            [
                fragment_conf.GetAtomPosition(i).x,
                fragment_conf.GetAtomPosition(i).y,
                fragment_conf.GetAtomPosition(i).z,
            ]
            for i in range(embedded_fragment.GetNumAtoms())
        ],
        dtype=float,
    )
    src_coord_map = {
        original_idx: all_coords[fragment_idx]
        for fragment_idx, original_idx in fragment.fragment_to_original.items()
    }
    src_planar_frame = _fragment_planar_donor_frame(
        scaffold_mol,
        fragment.atom_indices,
        target_atom_indices,
        src_coord_map,
    )

    if len(anchor_pairs) == 1:
        if reference_center is not None:
            ref = np.asarray(reference_center, dtype=float)
            frag_center = all_coords.mean(axis=0)
            if src_planar_frame is not None:
                src_dir = np.asarray(src_planar_frame["outward"], dtype=float)
                tgt_dir = ref - tgt[0]
            else:
                src_dir = frag_center - src[0]
                tgt_dir = tgt[0] - ref
            rot = _rotation_matrix_from_vectors(src_dir, tgt_dir)
            aligned = (all_coords - src[0]) @ rot.T + tgt[0]
        else:
            aligned = all_coords + (tgt[0] - src[0])
    elif len(anchor_pairs) == 2:
        rot = _rotation_matrix_from_vectors(src[1] - src[0], tgt[1] - tgt[0])
        src_mid = 0.5 * (src[0] + src[1])
        tgt_mid = 0.5 * (tgt[0] + tgt[1])
        aligned = (all_coords - src_mid) @ rot.T + tgt_mid
    else:
        src_centroid = src.mean(axis=0)
        tgt_centroid = tgt.mean(axis=0)
        covariance = (src - src_centroid).T @ (tgt - tgt_centroid)
        try:
            u, _s, vt = np.linalg.svd(covariance)
        except Exception:
            return False
        rot = vt.T @ u.T
        if float(np.linalg.det(rot)) < 0.0:
            vt[-1, :] *= -1.0
            rot = vt.T @ u.T
        aligned = (all_coords - src_centroid) @ rot.T + tgt_centroid

    if reference_center is not None and metal_idx is not None and len(anchor_pairs) in (1, 2):
        ref = np.asarray(reference_center, dtype=float)

        def _alignment_score(candidate_coords) -> float:
            donor_err = 0.0
            coord_map = {}
            donor_original_indices = []
            for pair_idx, (fragment_idx, original_idx) in enumerate(anchor_pairs):
                donor_original_indices.append(original_idx)
                donor_err += float(np.sum((candidate_coords[fragment_idx] - tgt[pair_idx]) ** 2))
            for fragment_idx, original_idx in fragment.fragment_to_original.items():
                coord_map[original_idx] = np.asarray(candidate_coords[fragment_idx], dtype=float)
            return (
                25.0 * donor_err
                + _secondary_non_donor_contact_penalty(
                    scaffold_mol,
                    metal_idx,
                    ref,
                    fragment.atom_indices,
                    donor_original_indices,
                    coord_map,
                )
                + _fragment_donor_selectivity_penalty(
                    scaffold_mol,
                    metal_idx,
                    ref,
                    fragment.atom_indices,
                    donor_original_indices,
                    coord_map,
                )
                + _fragment_environment_clash_penalty(
                    scaffold_mol,
                    fragment.atom_indices,
                    coord_map,
                    ignore_atom_indices={metal_idx},
                    primary_metal_idx=metal_idx,
                )
                + _planar_fragment_metal_coplanarity_penalty(
                    scaffold_mol,
                    fragment.atom_indices,
                    donor_original_indices,
                    ref,
                    coord_map,
                )
                + _planar_fragment_donor_approach_penalty(
                    scaffold_mol,
                    fragment.atom_indices,
                    donor_original_indices,
                    ref,
                    coord_map,
                )
                + _fragment_bite_direction_penalty(
                    scaffold_mol,
                    fragment.atom_indices,
                    donor_original_indices,
                    coord_map,
                    ref,
                )
            )

        aligned_coord_map = {
            original_idx: np.asarray(aligned[fragment_idx], dtype=float)
            for fragment_idx, original_idx in fragment.fragment_to_original.items()
        }
        aligned_planar_frame = _fragment_planar_donor_frame(
            scaffold_mol,
            fragment.atom_indices,
            [original_idx for _fragment_idx, original_idx in anchor_pairs],
            aligned_coord_map,
        )

        if len(anchor_pairs) == 1:
            pivot = tgt[0]
            axis = ref - tgt[0]
        else:
            pivot = tgt[0]
            axis = tgt[1] - tgt[0]
        axis_candidates = []
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm > 1e-12:
            axis_candidates.append(axis / axis_norm)
        if aligned_planar_frame is not None:
            normal = np.asarray(aligned_planar_frame["normal"], dtype=float)
            normal_norm = float(np.linalg.norm(normal))
            if normal_norm > 1e-12:
                axis_candidates.append(normal / normal_norm)
            outward = np.asarray(aligned_planar_frame["outward"], dtype=float)
            outward_norm = float(np.linalg.norm(outward))
            if outward_norm > 1e-12:
                axis_candidates.append(outward / outward_norm)
            if axis_candidates:
                cross_axis = np.cross(axis_candidates[0], normal if normal_norm > 1e-12 else outward)
                cross_norm = float(np.linalg.norm(cross_axis))
                if cross_norm > 1e-12:
                    axis_candidates.append(cross_axis / cross_norm)

        if axis_candidates:
            best_aligned = aligned.copy()
            best_score = _alignment_score(best_aligned)
            step = 5 if aligned_planar_frame is not None else 15
            for axis in axis_candidates:
                for angle_deg in range(-180, 181, step):
                    if angle_deg == 0:
                        continue
                    rot = _axis_angle_rotation_matrix(axis, math.radians(float(angle_deg)))
                    trial = (aligned - pivot) @ rot.T + pivot
                    trial_score = _alignment_score(trial)
                    if trial_score + 1e-9 < best_score:
                        best_score = trial_score
                        best_aligned = trial
            aligned = best_aligned

    for pair_idx, (fragment_idx, _original_idx) in enumerate(anchor_pairs[:max(0, exact_target_count)]):
        aligned[fragment_idx] = tgt[pair_idx]

    for fragment_idx, original_idx in fragment.fragment_to_original.items():
        scaffold_conf.SetAtomPosition(
            original_idx,
            Point3D(
                float(aligned[fragment_idx, 0]),
                float(aligned[fragment_idx, 1]),
                float(aligned[fragment_idx, 2]),
            ),
        )
    return True


def _hapto_primary_donor_indices(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> set:
    """Return non-metal donors bound directly to hapto metal centers."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return set()

    hapto_metals = {metal_idx for metal_idx, _grp in hapto_groups}
    donors: set = set()
    for metal_idx in hapto_metals:
        for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors():
            if nbr.GetSymbol() not in _METAL_SET:
                donors.add(nbr.GetIdx())
    return donors


def _secondary_metal_geometry_codes(n_coord: int, metal_symbol: str) -> List[str]:
    """Return geometry candidates for explicit secondary-metal placement."""
    if n_coord <= 1:
        return []
    if n_coord == 2:
        return ['LIN']
    if n_coord == 3:
        return ['TP', 'TS']
    if n_coord == 4:
        pref = _PREFERRED_CN4_GEOMETRY.get(metal_symbol, 'SQ')
        other = 'TH' if pref == 'SQ' else 'SQ'
        return [pref, other]
    if n_coord == 5:
        return ['TBP', 'SP']
    if n_coord == 6:
        return ['OH']
    if n_coord == 7:
        return ['PBP']
    if n_coord == 8:
        return ['SAP', 'DD']
    return []


def _fit_secondary_metal_model_to_targets(
    model_positions,
    target_positions,
    fit_indices: List[int],
):
    """Fit a full coordination model onto explicit donor targets."""
    import numpy as np

    model = np.asarray(model_positions, dtype=float)
    targets = np.asarray(target_positions, dtype=float)
    if model.ndim != 2 or targets.ndim != 2 or model.shape != targets.shape:
        return None
    if not fit_indices:
        return None

    model_fit = model[fit_indices]
    target_fit = targets[fit_indices]

    if len(fit_indices) == 1:
        delta = target_fit[0] - model_fit[0]
        transformed = model + delta
        metal_pos = delta
        return transformed, metal_pos

    if len(fit_indices) == 2:
        rot = _rotation_matrix_from_vectors(
            model_fit[1] - model_fit[0],
            target_fit[1] - target_fit[0],
        )
        model_mid = 0.5 * (model_fit[0] + model_fit[1])
        target_mid = 0.5 * (target_fit[0] + target_fit[1])
        transformed = (model - model_mid) @ rot.T + target_mid
        metal_pos = (-model_mid) @ rot.T + target_mid
        return transformed, metal_pos

    model_centroid = model_fit.mean(axis=0)
    target_centroid = target_fit.mean(axis=0)
    covariance = (model_fit - model_centroid).T @ (target_fit - target_centroid)
    try:
        u, _s, vt = np.linalg.svd(covariance)
    except Exception:
        return None
    rot = vt.T @ u.T
    if float(np.linalg.det(rot)) < 0.0:
        vt[-1, :] *= -1.0
        rot = vt.T @ u.T
    transformed = (model - model_centroid) @ rot.T + target_centroid
    metal_pos = (-model_centroid) @ rot.T + target_centroid
    return transformed, metal_pos


def _axis_angle_rotation_matrix(axis, angle_rad: float):
    """Return the rotation matrix for a rotation around ``axis``."""
    import numpy as np

    axis_arr = np.asarray(axis, dtype=float)
    norm = float(np.linalg.norm(axis_arr))
    if norm < 1e-12:
        return np.eye(3)
    axis_arr /= norm
    ux, uy, uz = axis_arr
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    one_c = 1.0 - c
    return np.array([
        [c + ux * ux * one_c, ux * uy * one_c - uz * s, ux * uz * one_c + uy * s],
        [uy * ux * one_c + uz * s, c + uy * uy * one_c, uy * uz * one_c - ux * s],
        [uz * ux * one_c - uy * s, uz * uy * one_c + ux * s, c + uz * uz * one_c],
    ], dtype=float)


def _project_displacements_to_rigid_body(
    positions: Dict[int, object],
    displacement_map: Dict[int, object],
    atom_indices: List[int],
) -> Dict[int, object]:
    """Approximate a set of atomic displacements by one rigid-body motion."""
    import numpy as np

    ordered = [atom_idx for atom_idx in atom_indices if atom_idx in positions]
    if len(ordered) < 2:
        return {
            atom_idx: np.asarray(
                displacement_map.get(atom_idx, np.zeros(3, dtype=float)),
                dtype=float,
            )
            for atom_idx in ordered
        }

    pts = np.asarray([positions[atom_idx] for atom_idx in ordered], dtype=float)
    disp = np.asarray(
        [displacement_map.get(atom_idx, np.zeros(3, dtype=float)) for atom_idx in ordered],
        dtype=float,
    )
    centroid = pts.mean(axis=0)

    a = np.zeros((3 * len(ordered), 6), dtype=float)
    b = disp.reshape(-1)
    ident = np.eye(3, dtype=float)
    for row_idx, atom_idx in enumerate(ordered):
        rel = np.asarray(positions[atom_idx], dtype=float) - centroid
        cross_block = np.array(
            [
                [0.0, rel[2], -rel[1]],
                [-rel[2], 0.0, rel[0]],
                [rel[1], -rel[0], 0.0],
            ],
            dtype=float,
        )
        start = 3 * row_idx
        a[start:start + 3, 0:3] = ident
        a[start:start + 3, 3:6] = cross_block

    try:
        sol, *_rest = np.linalg.lstsq(a, b, rcond=None)
    except Exception:
        return {
            atom_idx: np.asarray(
                displacement_map.get(atom_idx, np.zeros(3, dtype=float)),
                dtype=float,
            )
            for atom_idx in ordered
        }

    trans = np.asarray(sol[:3], dtype=float)
    omega = np.asarray(sol[3:6], dtype=float)
    projected: Dict[int, object] = {}
    for atom_idx in ordered:
        rel = np.asarray(positions[atom_idx], dtype=float) - centroid
        projected[atom_idx] = trans + np.cross(omega, rel)
    return projected


def _secondary_non_donor_contact_penalty(
    mol,
    metal_idx: int,
    candidate_metal_pos,
    atom_indices: List[int],
    donor_indices: List[int],
    coord_map: Dict[int, object],
) -> float:
    """Penalty for non-donor atoms collapsing onto a secondary metal center."""
    import numpy as np

    if mol is None or not atom_indices:
        return 0.0

    metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    donor_set = set(donor_indices)
    atom_set = set(atom_indices)
    alpha_atoms: set = set()
    beta_atoms: set = set()

    for donor_idx in donor_indices:
        atom = mol.GetAtomWithIdx(donor_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in atom_set and nbr_idx not in donor_set and nbr.GetAtomicNum() > 1:
                alpha_atoms.add(nbr_idx)
    for atom_idx in alpha_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in atom_set and nbr_idx not in donor_set and nbr_idx not in alpha_atoms and nbr.GetAtomicNum() > 1:
                beta_atoms.add(nbr_idx)

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    def _min_allowed(atom_idx: int, atom_sym: str) -> float:
        ml_ref = float(_get_ml_bond_length(metal_sym, atom_sym))
        min_allowed = max(1.70, 0.92 * ml_ref)
        if atom_idx in alpha_atoms:
            min_allowed = max(min_allowed, 2.05)
        elif atom_idx in beta_atoms:
            min_allowed = max(min_allowed, 1.85)
        atom = mol.GetAtomWithIdx(atom_idx)
        if _is_planar_like(atom):
            if atom_idx in alpha_atoms:
                min_allowed = max(min_allowed, 2.18)
            elif atom_idx in beta_atoms:
                min_allowed = max(min_allowed, 1.98)
            else:
                min_allowed = max(min_allowed, 1.82)
        return min_allowed

    metal_pos = np.asarray(candidate_metal_pos, dtype=float)
    penalty = 0.0
    for atom_idx in atom_indices:
        if atom_idx in donor_set:
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            continue
        coord = coord_map.get(atom_idx)
        if coord is None:
            continue
        atom_pos = np.asarray(coord, dtype=float)
        dist = float(np.linalg.norm(atom_pos - metal_pos))
        min_allowed = _min_allowed(atom_idx, atom.GetSymbol())
        if dist < min_allowed:
            gap = min_allowed - dist
            if _is_planar_like(atom):
                weight = 40.0 if atom_idx in alpha_atoms else 26.0
            else:
                weight = 28.0 if atom_idx in alpha_atoms else 18.0
            penalty += weight * gap * gap
        if dist < 1.55:
            penalty += 120.0 * (1.55 - dist) ** 2
    return penalty


def _secondary_fragment_donor_selectivity_penalty(
    mol,
    metal_idx: int,
    candidate_metal_pos,
    donor_indices: List[int],
    donor_to_fragment: Dict[int, _HybridHaptoFragment],
    coord_map: Dict[int, object],
) -> float:
    """Penalize non-donor atoms approaching a secondary metal as closely as the true donor."""
    import numpy as np

    if mol is None or not donor_indices or not coord_map:
        return 0.0

    metal_pos = np.asarray(candidate_metal_pos, dtype=float)
    donor_set = set(donor_indices)
    penalty = 0.0

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    def _is_carbonyl_like_carbon(atom_idx: int) -> bool:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            return False
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() < 1.5:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetSymbol() == 'O':
                return True
        return False

    for donor_idx in donor_indices:
        donor_pos = coord_map.get(donor_idx)
        fragment = donor_to_fragment.get(donor_idx)
        if donor_pos is None or fragment is None:
            continue
        donor_atom = mol.GetAtomWithIdx(donor_idx)
        donor_sym = donor_atom.GetSymbol()
        donor_dist = float(np.linalg.norm(np.asarray(donor_pos, dtype=float) - metal_pos))

        for atom_idx in fragment.atom_indices:
            if atom_idx in donor_set or atom_idx == metal_idx:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            atom_pos = coord_map.get(atom_idx)
            if atom_pos is None:
                continue
            atom_dist = float(np.linalg.norm(np.asarray(atom_pos, dtype=float) - metal_pos))
            try:
                topo_sep = len(Chem.GetShortestPath(mol, donor_idx, atom_idx)) - 1
            except Exception:
                topo_sep = 99
            if topo_sep < 1 or topo_sep > 3:
                continue

            planar_like = _is_planar_like(atom)
            delta = 0.14
            weight = 10.0
            if topo_sep == 1:
                delta = 0.28 if donor_sym == 'O' else 0.32
                weight = 28.0 if donor_sym == 'O' else 34.0
            elif topo_sep == 2:
                delta = 0.20 if donor_sym == 'O' else 0.24
                weight = 18.0 if donor_sym == 'O' else 24.0
            elif topo_sep == 3:
                delta = 0.10 if donor_sym == 'O' else 0.12
                weight = 10.0 if donor_sym == 'O' else 12.0

            if planar_like:
                delta += 0.08
                weight += 10.0
            if donor_sym == 'N' and atom.GetSymbol() == 'C':
                delta += 0.06
                weight += 8.0
            if donor_sym == 'O' and atom.GetSymbol() == 'C':
                delta += 0.05
                weight += 6.0
                if _is_carbonyl_like_carbon(atom_idx):
                    delta += 0.16
                    weight += 18.0

            preferred_min = donor_dist + delta
            if atom_dist < preferred_min:
                gap = preferred_min - atom_dist
                penalty += weight * gap * gap
            if atom_dist < donor_dist + 0.02:
                gap = donor_dist + 0.02 - atom_dist
                penalty += (34.0 + weight) * gap * gap

    return penalty


def _fragment_donor_selectivity_penalty(
    mol,
    metal_idx: int,
    candidate_metal_pos,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    coord_map: Dict[int, object],
) -> float:
    """Penalize non-donor atoms in one donor fragment approaching as closely as the true donor."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not donor_indices or not coord_map:
        return 0.0

    metal_pos = np.asarray(candidate_metal_pos, dtype=float)
    donor_set = set(donor_indices)
    fragment_set = set(fragment_atom_indices)
    penalty = 0.0

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    def _is_carbonyl_like_carbon(atom_idx: int) -> bool:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            return False
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() < 1.5:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetSymbol() == 'O':
                return True
        return False

    for donor_idx in donor_indices:
        donor_pos = coord_map.get(donor_idx)
        if donor_pos is None:
            continue
        donor_atom = mol.GetAtomWithIdx(donor_idx)
        donor_sym = donor_atom.GetSymbol()
        donor_dist = float(np.linalg.norm(np.asarray(donor_pos, dtype=float) - metal_pos))
        for atom_idx in fragment_atom_indices:
            if atom_idx not in fragment_set or atom_idx in donor_set or atom_idx == metal_idx:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            atom_pos = coord_map.get(atom_idx)
            if atom_pos is None:
                continue
            atom_dist = float(np.linalg.norm(np.asarray(atom_pos, dtype=float) - metal_pos))
            try:
                topo_sep = len(Chem.GetShortestPath(mol, donor_idx, atom_idx)) - 1
            except Exception:
                topo_sep = 99
            if topo_sep < 1 or topo_sep > 3:
                continue

            planar_like = _is_planar_like(atom)
            delta = 0.14
            weight = 10.0
            if topo_sep == 1:
                delta = 0.28 if donor_sym == 'O' else 0.32
                weight = 28.0 if donor_sym == 'O' else 34.0
            elif topo_sep == 2:
                delta = 0.20 if donor_sym == 'O' else 0.24
                weight = 18.0 if donor_sym == 'O' else 24.0
            elif topo_sep == 3:
                delta = 0.10 if donor_sym == 'O' else 0.12
                weight = 10.0 if donor_sym == 'O' else 12.0

            if planar_like:
                delta += 0.08
                weight += 10.0
            if donor_sym == 'N' and atom.GetSymbol() == 'C':
                delta += 0.06
                weight += 8.0
            if donor_sym == 'O' and atom.GetSymbol() == 'C':
                delta += 0.05
                weight += 6.0
                if _is_carbonyl_like_carbon(atom_idx):
                    delta += 0.16
                    weight += 18.0

            preferred_min = donor_dist + delta
            if atom_dist < preferred_min:
                gap = preferred_min - atom_dist
                penalty += weight * gap * gap
            if atom_dist < donor_dist + 0.02:
                gap = donor_dist + 0.02 - atom_dist
                penalty += (34.0 + weight) * gap * gap

    return penalty


def _planar_fragment_metal_coplanarity_penalty(
    mol,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    metal_pos,
    coord_map: Dict[int, object],
) -> float:
    """Penalty when a planar donor fragment is not approximately coplanar with the metal."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not donor_indices or not coord_map:
        return 0.0

    donor_set = set(donor_indices)

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    planar_atoms = tuple(
        atom_idx for atom_idx in fragment_atom_indices
        if _is_planar_like(mol.GetAtomWithIdx(atom_idx))
    )
    donor_overlap = [
        donor_idx for donor_idx in donor_indices
        if donor_idx in set(planar_atoms)
    ]
    if len(planar_atoms) < 4 or not donor_overlap:
        return 0.0

    pts = []
    for atom_idx in planar_atoms:
        arr = coord_map.get(atom_idx)
        if arr is None:
            return 0.0
        pts.append(np.asarray(arr, dtype=float))
    pts_arr = np.asarray(pts, dtype=float)
    centroid = pts_arr.mean(axis=0)
    try:
        _u, _s, vt = np.linalg.svd(pts_arr - centroid)
    except Exception:
        return 0.0
    normal = vt[-1]
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm < 1e-12:
        return 0.0
    normal = normal / normal_norm

    weight = 42.0 + 12.0 * len(donor_overlap)
    if any(mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'N' for donor_idx in donor_overlap):
        weight += 14.0
    if any(mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'O' for donor_idx in donor_overlap):
        weight += 14.0
    if any(mol.GetAtomWithIdx(donor_idx).GetIsAromatic() for donor_idx in donor_overlap):
        weight += 12.0

    plane_offset = float(np.dot(np.asarray(metal_pos, dtype=float) - centroid, normal))
    penalty = weight * plane_offset * plane_offset
    if abs(plane_offset) > 0.28:
        penalty += 75.0 * (abs(plane_offset) - 0.28) ** 2
    return penalty


def _fragment_planar_donor_frame(
    mol,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    coord_map: Dict[int, object],
):
    """Return a simple planar donor frame for rigid early docking."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not donor_indices or not coord_map:
        return None

    fragment_set = set(fragment_atom_indices)

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    planar_atoms = [
        atom_idx for atom_idx in fragment_atom_indices
        if _is_planar_like(mol.GetAtomWithIdx(atom_idx))
    ]
    donor_overlap = [
        donor_idx for donor_idx in donor_indices
        if donor_idx in fragment_set and donor_idx in planar_atoms
    ]
    if len(planar_atoms) < 4 or not donor_overlap:
        return None

    pts = []
    for atom_idx in planar_atoms:
        arr = coord_map.get(atom_idx)
        if arr is None:
            return None
        pts.append(np.asarray(arr, dtype=float))
    pts_arr = np.asarray(pts, dtype=float)
    centroid = pts_arr.mean(axis=0)
    try:
        _u, _s, vt = np.linalg.svd(pts_arr - centroid)
    except Exception:
        return None
    normal = vt[-1]
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm < 1e-12:
        return None
    normal = normal / normal_norm

    donor_pts = np.asarray([coord_map[atom_idx] for atom_idx in donor_overlap], dtype=float)
    donor_centroid = donor_pts.mean(axis=0)

    body_pts = []
    for atom_idx in planar_atoms:
        if atom_idx in donor_overlap:
            continue
        body_pts.append(np.asarray(coord_map[atom_idx], dtype=float))
    if body_pts:
        body_centroid = np.asarray(body_pts, dtype=float).mean(axis=0)
    else:
        body_centroid = centroid

    outward = donor_centroid - body_centroid
    outward = outward - float(np.dot(outward, normal)) * normal
    outward_norm = float(np.linalg.norm(outward))
    if outward_norm < 1e-12:
        outward = donor_centroid - centroid
        outward = outward - float(np.dot(outward, normal)) * normal
        outward_norm = float(np.linalg.norm(outward))
    if outward_norm < 1e-12:
        return None
    outward = outward / outward_norm

    return {
        "planar_atoms": tuple(planar_atoms),
        "donor_overlap": tuple(donor_overlap),
        "centroid": centroid,
        "donor_centroid": donor_centroid,
        "normal": normal,
        "outward": outward,
    }


def _planar_fragment_donor_approach_penalty(
    mol,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    metal_pos,
    coord_map: Dict[int, object],
) -> float:
    """Penalty when a planar donor fragment presents adjacent pi atoms instead of the donor side."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not donor_indices or not coord_map:
        return 0.0

    fragment_set = set(fragment_atom_indices)

    def _is_planar_like(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    planar_atoms = tuple(
        atom_idx for atom_idx in fragment_atom_indices
        if _is_planar_like(mol.GetAtomWithIdx(atom_idx))
    )
    donor_overlap = [
        donor_idx for donor_idx in donor_indices
        if donor_idx in fragment_set and donor_idx in planar_atoms
    ]
    if len(planar_atoms) < 4 or not donor_overlap:
        return 0.0

    pts = []
    for atom_idx in planar_atoms:
        arr = coord_map.get(atom_idx)
        if arr is None:
            return 0.0
        pts.append(np.asarray(arr, dtype=float))
    pts_arr = np.asarray(pts, dtype=float)
    centroid = pts_arr.mean(axis=0)
    try:
        _u, _s, vt = np.linalg.svd(pts_arr - centroid)
    except Exception:
        return 0.0
    normal = vt[-1]
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm < 1e-12:
        return 0.0
    normal = normal / normal_norm

    penalty = 0.0
    metal_vec_ref = np.asarray(metal_pos, dtype=float)
    for donor_idx in donor_overlap:
        donor_atom = mol.GetAtomWithIdx(donor_idx)
        donor_pos = coord_map.get(donor_idx)
        if donor_pos is None:
            continue
        donor_pos = np.asarray(donor_pos, dtype=float)

        outward_terms = []
        for nbr in donor_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in fragment_set:
                continue
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            if not _is_planar_like(nbr_atom):
                continue
            nbr_pos = coord_map.get(nbr_idx)
            if nbr_pos is None:
                continue
            vec = donor_pos - np.asarray(nbr_pos, dtype=float)
            vec = vec - float(np.dot(vec, normal)) * normal
            vec_norm = float(np.linalg.norm(vec))
            if vec_norm < 1e-8:
                continue
            outward_terms.append(vec / vec_norm)

        if outward_terms:
            outward_vec = np.sum(np.asarray(outward_terms, dtype=float), axis=0)
        else:
            outward_vec = donor_pos - centroid
            outward_vec = outward_vec - float(np.dot(outward_vec, normal)) * normal
        outward_norm = float(np.linalg.norm(outward_vec))
        if outward_norm < 1e-8:
            continue
        outward_unit = outward_vec / outward_norm

        metal_vec = metal_vec_ref - donor_pos
        metal_in_plane = metal_vec - float(np.dot(metal_vec, normal)) * normal
        metal_in_plane_norm = float(np.linalg.norm(metal_in_plane))
        if metal_in_plane_norm < 1e-8:
            continue
        metal_unit = metal_in_plane / metal_in_plane_norm

        cosang = float(np.dot(outward_unit, metal_unit))
        donor_sym = donor_atom.GetSymbol()
        min_cos = 0.18
        weight = 10.0
        if donor_sym == 'N':
            min_cos = 0.34 if donor_atom.GetIsAromatic() else 0.26
            weight = 18.0 if donor_atom.GetIsAromatic() else 14.0
        elif donor_sym == 'O':
            min_cos = 0.28
            weight = 16.0
            if any(bond.GetIsConjugated() for bond in donor_atom.GetBonds()):
                min_cos += 0.06
                weight += 6.0
        elif donor_sym == 'S':
            min_cos = 0.22
            weight = 12.0

        if cosang < min_cos:
            penalty += weight * (min_cos - cosang) ** 2
        if cosang < -0.05:
            penalty += (weight + 18.0) * (abs(cosang) + 0.05) ** 2

    return penalty


def _fragment_environment_clash_penalty(
    mol,
    fragment_atom_indices: List[int],
    coord_map: Dict[int, object],
    ignore_atom_indices: Optional[set] = None,
    primary_metal_idx: Optional[int] = None,
) -> float:
    """Penalty for a candidate rigid fragment pose colliding with the fixed environment."""
    import numpy as np

    if mol is None or not fragment_atom_indices or not coord_map:
        return 0.0

    try:
        conf = mol.GetConformer()
    except Exception:
        return 0.0

    ignore = set(ignore_atom_indices or set())
    fragment_set = set(fragment_atom_indices)
    penalty = 0.0

    for atom_idx in fragment_atom_indices:
        if atom_idx in ignore:
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        sym_i = atom.GetSymbol()
        if atom.GetAtomicNum() <= 1 or sym_i in _METAL_SET:
            continue
        pos_i = coord_map.get(atom_idx)
        if pos_i is None:
            continue
        pos_i = np.asarray(pos_i, dtype=float)
        if pos_i.shape != (3,) or not np.all(np.isfinite(pos_i)):
            continue
        for other_idx in range(mol.GetNumAtoms()):
            if other_idx in fragment_set or other_idx in ignore:
                continue
            other = mol.GetAtomWithIdx(other_idx)
            sym_j = other.GetSymbol()
            if other.GetAtomicNum() <= 1 or sym_j in _METAL_SET:
                continue
            if mol.GetBondBetweenAtoms(atom_idx, other_idx) is not None:
                continue
            other_pos = conf.GetAtomPosition(other_idx)
            pos_j = np.array([other_pos.x, other_pos.y, other_pos.z], dtype=float)
            dist = float(np.linalg.norm(pos_i - pos_j))
            other_metal_bound = any(
                nbr.GetSymbol() in _METAL_SET and nbr.GetIdx() != primary_metal_idx
                for nbr in other.GetNeighbors()
            )
            min_dist = max(
                1.76,
                _COVALENT_RADII.get(sym_i, 0.76) + _COVALENT_RADII.get(sym_j, 0.76) + 0.34,
            )
            if sym_i in {'N', 'O', 'S', 'P'} or sym_j in {'N', 'O', 'S', 'P'}:
                min_dist = max(min_dist, 1.90)
            if other_metal_bound:
                min_dist = max(min_dist, 2.12)
            if dist < min_dist:
                weight = 18.0 if other_metal_bound else 10.0
                penalty += weight * (min_dist - dist) ** 2
            if dist < 1.55:
                weight = 120.0 if other_metal_bound else 70.0
                penalty += weight * (1.55 - dist) ** 2

    return penalty


def _fragment_bite_direction_penalty(
    mol,
    fragment_atom_indices: List[int],
    donor_indices: List[int],
    coord_map: Dict[int, object],
    metal_pos,
) -> float:
    """Penalty when a multidentate rigid fragment presents its back side to the metal."""
    import numpy as np

    if mol is None or len(donor_indices) < 2 or not coord_map:
        return 0.0

    donor_coords = []
    body_coords = []
    fallback_coords = []
    for atom_idx in fragment_atom_indices:
        coord = coord_map.get(atom_idx)
        if coord is None:
            continue
        arr = np.asarray(coord, dtype=float)
        if arr.shape != (3,) or not np.all(np.isfinite(arr)):
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            continue
        fallback_coords.append(arr)
        if atom_idx in donor_indices:
            donor_coords.append(arr)
        else:
            body_coords.append(arr)

    if len(donor_coords) < 2:
        return 0.0
    if not body_coords:
        body_coords = fallback_coords
    if not body_coords:
        return 0.0

    donor_centroid = np.mean(np.asarray(donor_coords, dtype=float), axis=0)
    body_centroid = np.mean(np.asarray(body_coords, dtype=float), axis=0)
    pocket_vec = donor_centroid - body_centroid
    metal_vec = np.asarray(metal_pos, dtype=float) - donor_centroid
    pocket_norm = float(np.linalg.norm(pocket_vec))
    metal_norm = float(np.linalg.norm(metal_vec))
    if pocket_norm < 1e-8 or metal_norm < 1e-8:
        return 0.0

    pocket_unit = pocket_vec / pocket_norm
    metal_unit = metal_vec / metal_norm
    cosang = float(np.dot(pocket_unit, metal_unit))
    if cosang >= 0.55:
        return 0.0
    gap = 0.55 - cosang
    return 6.0 * gap * gap


def _fit_secondary_metal_geometry_from_donors(
    donor_positions,
    donor_symbols: List[str],
    metal_symbol: str,
    constrained_indices: Optional[List[int]] = None,
    bite_distance_constraints: Optional[List[Tuple[int, int, float]]] = None,
    donor_target_lengths: Optional[List[float]] = None,
    donor_fit_weights: Optional[List[float]] = None,
):
    """Fit an ideal coordination geometry while optionally honoring anchor donors."""
    fits = _enumerate_secondary_metal_geometry_fits(
        donor_positions,
        donor_symbols,
        metal_symbol,
        constrained_indices=constrained_indices,
        bite_distance_constraints=bite_distance_constraints,
        donor_target_lengths=donor_target_lengths,
        donor_fit_weights=donor_fit_weights,
        max_candidates=1,
    )
    if not fits:
        return None
    metal_pos, geom_code, transformed, _score = fits[0]
    return metal_pos, geom_code, transformed


def _enumerate_secondary_metal_geometry_fits(
    donor_positions,
    donor_symbols: List[str],
    metal_symbol: str,
    constrained_indices: Optional[List[int]] = None,
    bite_distance_constraints: Optional[List[Tuple[int, int, float]]] = None,
    donor_target_lengths: Optional[List[float]] = None,
    donor_fit_weights: Optional[List[float]] = None,
    max_candidates: int = 4,
) -> List[Tuple[object, str, object, float]]:
    """Return ranked geometry fits for a secondary metal donor set."""
    try:
        import itertools
        import numpy as np
    except ImportError:
        return []

    donors = np.asarray(donor_positions, dtype=float)
    if donors.ndim != 2 or donors.shape[0] != len(donor_symbols) or donors.shape[0] < 2:
        return []

    n_coord = donors.shape[0]
    target_lengths = [
        float(donor_target_lengths[i]) if donor_target_lengths and i < len(donor_target_lengths)
        else float(_get_ml_bond_length(metal_symbol, donor_symbols[i]))
        for i in range(n_coord)
    ]
    fit_weights = np.asarray(
        [
            float(donor_fit_weights[i]) if donor_fit_weights and i < len(donor_fit_weights) else 1.0
            for i in range(n_coord)
        ],
        dtype=float,
    )
    geometry_codes = _secondary_metal_geometry_codes(n_coord, metal_symbol)
    if not geometry_codes:
        return []

    fit_indices = sorted({
        int(idx) for idx in (constrained_indices or [])
        if 0 <= int(idx) < n_coord
    })
    if len(fit_indices) < 2:
        fit_indices = list(range(n_coord))

    nonpreferred_cn4_penalty = 0.0
    if n_coord == 4:
        preferred_cn4 = _PREFERRED_CN4_GEOMETRY.get(metal_symbol)
        if preferred_cn4:
            if metal_symbol in {'Pt', 'Pd', 'Rh', 'Ir', 'Au'}:
                nonpreferred_cn4_penalty = 0.24
            elif metal_symbol == 'Ni':
                nonpreferred_cn4_penalty = 0.12
            else:
                nonpreferred_cn4_penalty = 0.16
        else:
            preferred_cn4 = None
    else:
        preferred_cn4 = None

    ranked_candidates: List[Tuple[float, float, float, float, object, str, object]] = []
    for geom_code in geometry_codes:
        vectors = _TOPO_GEOMETRY_VECTORS.get(geom_code)
        if not vectors or len(vectors) != n_coord:
            continue
        normed_vectors = []
        for vec in vectors:
            arr = np.asarray(vec, dtype=float)
            norm = float(np.linalg.norm(arr))
            if norm < 1e-12:
                break
            normed_vectors.append(arr / norm)
        if len(normed_vectors) != n_coord:
            continue

        if n_coord <= 6:
            permutations = itertools.permutations(range(n_coord))
        else:
            identity = tuple(range(n_coord))
            permutations = [identity, tuple(reversed(identity))]

        for perm in permutations:
            model = np.asarray(
                [
                    normed_vectors[perm[i]] * target_lengths[i]
                    for i in range(n_coord)
                ],
                dtype=float,
            )
            fitted = _fit_secondary_metal_model_to_targets(model, donors, fit_indices)
            if fitted is None:
                continue
            transformed, metal_pos = fitted
            fit_resid = np.sum((transformed[fit_indices] - donors[fit_indices]) ** 2, axis=1)
            fit_w = fit_weights[fit_indices]
            rmsd_fit = float(np.sqrt(np.sum(fit_w * fit_resid) / max(np.sum(fit_w), 1e-12)))
            all_resid = np.sum((transformed - donors) ** 2, axis=1)
            rmsd_all = float(np.sqrt(np.sum(fit_weights * all_resid) / max(np.sum(fit_weights), 1e-12)))
            bite_rmsd = 0.0
            if bite_distance_constraints:
                bite_errors = []
                for idx_i, idx_j, expected_dist in bite_distance_constraints:
                    if idx_i < 0 or idx_j < 0 or idx_i >= n_coord or idx_j >= n_coord:
                        continue
                    actual_dist = float(np.linalg.norm(transformed[idx_i] - transformed[idx_j]))
                    bite_errors.append((actual_dist - expected_dist) ** 2)
                if bite_errors:
                    bite_rmsd = float(np.sqrt(np.mean(bite_errors)))
                score = rmsd_fit + 0.10 * rmsd_all + 0.55 * bite_rmsd
                if preferred_cn4 is not None and geom_code != preferred_cn4:
                    score += nonpreferred_cn4_penalty
                ranked_candidates.append(
                    (
                        float(score),
                        float(rmsd_fit),
                        float(rmsd_all),
                        float(bite_rmsd),
                        np.asarray(metal_pos, dtype=float),
                        geom_code,
                        np.asarray(transformed, dtype=float),
                    )
                )

    if not ranked_candidates:
        return []

    ranked_candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))

    selected: List[Tuple[object, str, object, float]] = []

    def _is_distinct_geometry_candidate(
        transformed,
        geom_code: str,
        existing: List[Tuple[object, str, object, float]],
    ) -> bool:
        transformed_arr = np.asarray(transformed, dtype=float)
        for _pos, prev_geom, prev_transformed, _score in existing:
            if geom_code != prev_geom:
                continue
            prev_arr = np.asarray(prev_transformed, dtype=float)
            if prev_arr.shape != transformed_arr.shape:
                continue
            deltas = np.linalg.norm(prev_arr - transformed_arr, axis=1)
            if float(np.max(deltas)) < 0.22:
                return False
        return True

    first = ranked_candidates[0]
    selected.append((first[4], first[5], first[6], first[0]))
    seen_geometries = {first[5]}

    for candidate in ranked_candidates[1:]:
        if len(selected) >= max(1, int(max_candidates)):
            break
        geom_code = candidate[5]
        if geom_code in seen_geometries:
            continue
        if not _is_distinct_geometry_candidate(candidate[6], geom_code, selected):
            continue
        selected.append((candidate[4], geom_code, candidate[6], candidate[0]))
        seen_geometries.add(geom_code)

    for candidate in ranked_candidates[1:]:
        if len(selected) >= max(1, int(max_candidates)):
            break
        geom_code = candidate[5]
        if not _is_distinct_geometry_candidate(candidate[6], geom_code, selected):
            continue
        selected.append((candidate[4], geom_code, candidate[6], candidate[0]))

    return selected[:max(1, int(max_candidates))]


def _fit_secondary_metal_position_from_donors(
    donor_positions,
    donor_symbols: List[str],
    metal_symbol: str,
    donor_target_lengths: Optional[List[float]] = None,
    donor_fit_weights: Optional[List[float]] = None,
):
    """Fit a secondary metal center to already-placed donor atoms."""
    fits = _enumerate_secondary_metal_position_fits(
        donor_positions,
        donor_symbols,
        metal_symbol,
        donor_target_lengths=donor_target_lengths,
        donor_fit_weights=donor_fit_weights,
        max_candidates=1,
    )
    if not fits:
        return None
    metal_pos, geom_code, _score = fits[0]
    return metal_pos, geom_code


def _enumerate_secondary_metal_position_fits(
    donor_positions,
    donor_symbols: List[str],
    metal_symbol: str,
    donor_target_lengths: Optional[List[float]] = None,
    donor_fit_weights: Optional[List[float]] = None,
    max_candidates: int = 4,
) -> List[Tuple[object, str, float]]:
    """Return ranked metal-position fits for already placed donors."""
    try:
        import itertools
        import numpy as np
    except ImportError:
        return []

    donors = np.asarray(donor_positions, dtype=float)
    if donors.ndim != 2 or donors.shape[0] != len(donor_symbols) or donors.shape[0] < 2:
        return []

    n_coord = donors.shape[0]
    target_lengths = [
        float(donor_target_lengths[i]) if donor_target_lengths and i < len(donor_target_lengths)
        else float(_get_ml_bond_length(metal_symbol, donor_symbols[i]))
        for i in range(n_coord)
    ]
    fit_weights = np.asarray(
        [
            float(donor_fit_weights[i]) if donor_fit_weights and i < len(donor_fit_weights) else 1.0
            for i in range(n_coord)
        ],
        dtype=float,
    )
    geometry_codes = _secondary_metal_geometry_codes(n_coord, metal_symbol)
    if not geometry_codes:
        return []

    ranked_candidates: List[Tuple[float, float, object, str]] = []
    for geom_code in geometry_codes:
        vectors = _TOPO_GEOMETRY_VECTORS.get(geom_code)
        if not vectors or len(vectors) != n_coord:
            continue
        normed_vectors = []
        for vec in vectors:
            arr = np.asarray(vec, dtype=float)
            norm = float(np.linalg.norm(arr))
            if norm < 1e-12:
                break
            normed_vectors.append(arr / norm)
        if len(normed_vectors) != n_coord:
            continue

        if n_coord <= 6:
            permutations = itertools.permutations(range(n_coord))
        else:
            identity = tuple(range(n_coord))
            permutations = [identity, tuple(reversed(identity))]

        for perm in permutations:
            model = np.asarray(
                [
                    normed_vectors[perm[i]] * target_lengths[i]
                    for i in range(n_coord)
                ],
                dtype=float,
            )
            model_centroid = model.mean(axis=0)
            donor_centroid = donors.mean(axis=0)
            covariance = (model - model_centroid).T @ (donors - donor_centroid)
            try:
                u, _s, vt = np.linalg.svd(covariance)
            except Exception:
                continue
            rot = vt.T @ u.T
            if float(np.linalg.det(rot)) < 0.0:
                vt[-1, :] *= -1.0
                rot = vt.T @ u.T
            fitted = (model - model_centroid) @ rot.T + donor_centroid
            resid = np.sum((fitted - donors) ** 2, axis=1)
            rmsd = float(np.sqrt(np.sum(fit_weights * resid) / max(np.sum(fit_weights), 1e-12)))
            metal_pos = donor_centroid - model_centroid @ rot.T
            score = rmsd
            ranked_candidates.append(
                (float(score), float(rmsd), np.asarray(metal_pos, dtype=float), geom_code)
            )

    if not ranked_candidates:
        return []

    ranked_candidates.sort(key=lambda item: (item[0], item[1]))

    selected: List[Tuple[object, str, float]] = []

    def _is_distinct_position_candidate(
        metal_pos,
        geom_code: str,
        existing: List[Tuple[object, str, float]],
    ) -> bool:
        metal_arr = np.asarray(metal_pos, dtype=float)
        for prev_pos, prev_geom, _score in existing:
            if geom_code != prev_geom:
                continue
            prev_arr = np.asarray(prev_pos, dtype=float)
            if float(np.linalg.norm(prev_arr - metal_arr)) < 0.20:
                return False
        return True

    first = ranked_candidates[0]
    selected.append((first[2], first[3], first[0]))
    seen_geometries = {first[3]}

    for candidate in ranked_candidates[1:]:
        if len(selected) >= max(1, int(max_candidates)):
            break
        geom_code = candidate[3]
        if geom_code in seen_geometries:
            continue
        if not _is_distinct_position_candidate(candidate[2], geom_code, selected):
            continue
        selected.append((candidate[2], geom_code, candidate[0]))
        seen_geometries.add(geom_code)

    for candidate in ranked_candidates[1:]:
        if len(selected) >= max(1, int(max_candidates)):
            break
        geom_code = candidate[3]
        if not _is_distinct_position_candidate(candidate[2], geom_code, selected):
            continue
        selected.append((candidate[2], geom_code, candidate[0]))

    return selected[:max(1, int(max_candidates))]


def _secondary_module_fragment_is_movable(
    fragment: _HybridHaptoFragment,
    metal_idx: int,
) -> bool:
    """Return whether a fragment may be moved as part of one secondary-metal module."""
    return (
        metal_idx in fragment.metal_neighbor_indices
        and len(fragment.metal_neighbor_indices) == 1
        and not fragment.bridging_donor_indices
        and not fragment.hapto_group_ids
    )


def _find_scaffold_secondary_branch(
    mol,
    fragment: _HybridHaptoFragment,
    hapto_groups: List[Tuple[int, List[int]]],
    donor_idx: int,
) -> Optional[Tuple[int, set]]:
    """Return a rotatable donor branch off a fixed hapto/bridge scaffold, if present."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None:
        return None
    if donor_idx not in fragment.atom_indices or donor_idx in fragment.bridging_donor_indices:
        return None

    fragment_atoms = set(fragment.atom_indices)
    core_seed_atoms: set = set(fragment.bridging_donor_indices)
    for group_id in fragment.hapto_group_ids:
        if 0 <= int(group_id) < len(hapto_groups):
            _metal_idx, group_atoms = hapto_groups[int(group_id)]
            core_seed_atoms.update(atom_idx for atom_idx in group_atoms if atom_idx in fragment_atoms)
    for other_donor in fragment.donor_atom_indices:
        if other_donor != donor_idx and other_donor in fragment_atoms:
            core_seed_atoms.add(other_donor)
    if donor_idx in core_seed_atoms:
        return None

    movable_pool = fragment_atoms - core_seed_atoms
    if donor_idx not in movable_pool:
        return None

    branch_atoms: set = set()
    queue = [donor_idx]
    while queue:
        atom_idx = queue.pop()
        if atom_idx in branch_atoms or atom_idx not in movable_pool:
            continue
        branch_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in movable_pool and nbr_idx not in branch_atoms:
                queue.append(nbr_idx)

    if donor_idx not in branch_atoms or len(branch_atoms) < 2:
        return None

    attachment_atoms: set = set()
    for atom_idx in branch_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_seed_atoms:
                attachment_atoms.add(nbr_idx)
    interface_atoms = {
        atom_idx
        for atom_idx in branch_atoms
        if any(nbr.GetIdx() in core_seed_atoms for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors())
    }
    if len(interface_atoms) == 1:
        pivot_idx = next(iter(interface_atoms))
        branch_atoms = set(branch_atoms)
        branch_atoms.discard(pivot_idx)
    else:
        if len(attachment_atoms) != 1:
            return None
        pivot_idx = next(iter(attachment_atoms))

    if donor_idx == pivot_idx or donor_idx not in branch_atoms:
        return None
    return pivot_idx, branch_atoms


def _reorient_scaffold_secondary_branch(
    mol,
    fragment: _HybridHaptoFragment,
    hapto_groups: List[Tuple[int, List[int]]],
    metal_idx: int,
    donor_idx: int,
    donor_target,
) -> bool:
    """Rotate a branch off a fixed hapto/bridge core toward a secondary-metal target."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    branch_info = _find_scaffold_secondary_branch(
        mol,
        fragment,
        hapto_groups,
        donor_idx,
    )
    if branch_info is None:
        return False
    pivot_idx, branch_atoms = branch_info

    pivot_pos = np.array(conf.GetAtomPosition(pivot_idx), dtype=float)
    current_coords = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in branch_atoms
    }
    donor_target_vec = np.asarray(donor_target, dtype=float)
    if donor_target_vec.shape != (3,) or not np.all(np.isfinite(donor_target_vec)):
        return False

    current_donor_vec = current_coords[donor_idx] - pivot_pos
    target_donor_vec = donor_target_vec - pivot_pos
    if float(np.linalg.norm(current_donor_vec)) < 1e-12 or float(np.linalg.norm(target_donor_vec)) < 1e-12:
        return False

    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
    base_rot = _rotation_matrix_from_vectors(current_donor_vec, target_donor_vec)
    rotated_coords = {
        atom_idx: (coord - pivot_pos) @ base_rot.T + pivot_pos
        for atom_idx, coord in current_coords.items()
    }
    branch_frame = _fragment_planar_donor_frame(
        mol,
        sorted(branch_atoms),
        [donor_idx],
        rotated_coords,
    )
    if branch_frame is not None:
        outward_src = np.asarray(branch_frame["outward"], dtype=float)
        metal_dir = metal_pos - donor_target_vec
        metal_dir = metal_dir - float(np.dot(metal_dir, branch_frame["normal"])) * np.asarray(branch_frame["normal"], dtype=float)
        if float(np.linalg.norm(metal_dir)) > 1e-8:
            outward_tgt = metal_dir / float(np.linalg.norm(metal_dir))
            donor_center = rotated_coords[donor_idx]
            branch_rot = _rotation_matrix_from_vectors(outward_src, outward_tgt)
            trial_coords = {
                atom_idx: (coord - donor_center) @ branch_rot.T + donor_center
                for atom_idx, coord in rotated_coords.items()
            }
            if donor_idx in trial_coords:
                donor_shift = donor_target_vec - trial_coords[donor_idx]
                trial_coords = {
                    atom_idx: coord + donor_shift
                    for atom_idx, coord in trial_coords.items()
                }
            rotated_coords = trial_coords

    def _score(coord_map: Dict[int, object]) -> Tuple[float, float]:
        donor_err = float(np.linalg.norm(coord_map[donor_idx] - donor_target_vec))
        penalty = _secondary_non_donor_contact_penalty(
            mol,
            metal_idx,
            metal_pos,
            sorted(branch_atoms),
            [donor_idx],
            coord_map,
        )
        penalty += _fragment_donor_selectivity_penalty(
            mol,
            metal_idx,
            metal_pos,
            sorted(branch_atoms),
            [donor_idx],
            coord_map,
        )
        penalty += _planar_fragment_metal_coplanarity_penalty(
            mol,
            sorted(branch_atoms),
            [donor_idx],
            metal_pos,
            coord_map,
        )
        penalty += _planar_fragment_donor_approach_penalty(
            mol,
            sorted(branch_atoms),
            [donor_idx],
            metal_pos,
            coord_map,
        )
        return donor_err, 30.0 * donor_err * donor_err + penalty

    current_err, current_score = _score(current_coords)
    best_coords = rotated_coords
    best_err, best_score = _score(rotated_coords)

    axis = rotated_coords[donor_idx] - pivot_pos
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm > 1e-12:
        axis /= axis_norm
        step = 5 if branch_frame is not None else 15
        for angle_deg in range(-180, 181, step):
            if angle_deg == 0:
                continue
            rot = _axis_angle_rotation_matrix(axis, math.radians(float(angle_deg)))
            trial_coords = {
                atom_idx: (coord - pivot_pos) @ rot.T + pivot_pos
                for atom_idx, coord in rotated_coords.items()
            }
            trial_err, trial_score = _score(trial_coords)
            if (
                trial_score + 1e-9 < best_score
                or (
                    abs(trial_score - best_score) < 1e-9
                    and trial_err + 1e-6 < best_err
                )
            ):
                best_coords = trial_coords
                best_err = trial_err
                best_score = trial_score

    if best_score + 1e-9 >= current_score and best_err + 1e-6 >= current_err:
        return False

    for atom_idx, coord in best_coords.items():
        conf.SetAtomPosition(atom_idx, Point3D(float(coord[0]), float(coord[1]), float(coord[2])))
    return True


def _optimize_secondary_fragment_pose(
    mol,
    fragment: _HybridHaptoFragment,
    metal_idx: int,
    donor_target_map: Dict[int, object],
) -> bool:
    """Improve a rigid donor fragment around one anchored donor without distorting it."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None or not donor_target_map:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
    donor_indices = [
        donor_idx for donor_idx in fragment.donor_atom_indices
        if donor_idx in donor_target_map
    ]
    if len(donor_indices) < 2:
        return False

    base_coords = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in fragment.atom_indices
    }
    metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    target_lengths = {
        donor_idx: float(_get_ml_bond_length(metal_sym, mol.GetAtomWithIdx(donor_idx).GetSymbol()))
        for donor_idx in donor_indices
    }

    donor_errors = {
        donor_idx: abs(float(np.linalg.norm(base_coords[donor_idx] - metal_pos)) - target_lengths[donor_idx])
        for donor_idx in donor_indices
    }
    pivot_idx = min(donor_indices, key=lambda donor_idx: donor_errors[donor_idx])
    pivot_pos = np.array(base_coords[pivot_idx], dtype=float)
    current_errors = donor_errors.copy()

    def _score(coord_map: Dict[int, object]) -> float:
        score = 0.0
        for donor_idx in donor_indices:
            donor_pos = coord_map[donor_idx]
            target_len = target_lengths[donor_idx]
            err = float(np.linalg.norm(donor_pos - metal_pos)) - target_len
            weight = 0.3 if donor_idx == pivot_idx else 3.0
            score += weight * err * err
        score += _secondary_non_donor_contact_penalty(
            mol,
            metal_idx,
            metal_pos,
            fragment.atom_indices,
            donor_indices,
            coord_map,
        )
        score += _fragment_donor_selectivity_penalty(
            mol,
            metal_idx,
            metal_pos,
            fragment.atom_indices,
            donor_indices,
            coord_map,
        )
        score += _fragment_environment_clash_penalty(
            mol,
            fragment.atom_indices,
            coord_map,
            ignore_atom_indices={metal_idx},
            primary_metal_idx=metal_idx,
        )
        score += _planar_fragment_metal_coplanarity_penalty(
            mol,
            fragment.atom_indices,
            donor_indices,
            metal_pos,
            coord_map,
        )
        score += _planar_fragment_donor_approach_penalty(
            mol,
            fragment.atom_indices,
            donor_indices,
            metal_pos,
            coord_map,
        )
        score += _fragment_bite_direction_penalty(
            mol,
            fragment.atom_indices,
            donor_indices,
            coord_map,
            metal_pos,
        )
        return score

    best_coords = {idx: vec.copy() for idx, vec in base_coords.items()}
    best_score = _score(best_coords)
    best_errors = current_errors.copy()
    donor_axis = metal_pos - pivot_pos
    donor_axis_norm = float(np.linalg.norm(donor_axis))
    if donor_axis_norm < 1e-12:
        return False
    donor_axis /= donor_axis_norm

    rng = np.random.default_rng(
        int((metal_idx + 1) * 1000 + min(fragment.atom_indices) + len(fragment.atom_indices))
    )
    axes = [donor_axis]
    for donor_idx in donor_indices:
        if donor_idx == pivot_idx:
            continue
        pair_axis = np.asarray(base_coords[donor_idx] - pivot_pos, dtype=float)
        pair_norm = float(np.linalg.norm(pair_axis))
        if pair_norm < 1e-12:
            continue
        axes.insert(0, pair_axis / pair_norm)
        break
    for _ in range(192):
        axis = rng.normal(size=3)
        norm = float(np.linalg.norm(axis))
        if norm < 1e-12:
            continue
        axes.append(axis / norm)

    improved = False
    for axis in axes:
        for angle_rad in (
            math.radians(v) for v in (-180, -150, -120, -90, -75, -60, -45, -30, -20, -10, 10, 20, 30, 45, 60, 75, 90, 120, 150, 180)
        ):
            rot = _axis_angle_rotation_matrix(axis, angle_rad)
            trial_coords = {}
            for atom_idx, vec in base_coords.items():
                if atom_idx == pivot_idx:
                    trial_coords[atom_idx] = vec.copy()
                    continue
                trial_coords[atom_idx] = (vec - pivot_pos) @ rot.T + pivot_pos
            trial_score = _score(trial_coords)
            trial_errors = {
                donor_idx: abs(float(np.linalg.norm(trial_coords[donor_idx] - metal_pos)) - target_lengths[donor_idx])
                for donor_idx in donor_indices
            }
            trial_max_nonpivot = max(
                (err for donor_idx, err in trial_errors.items() if donor_idx != pivot_idx),
                default=0.0,
            )
            best_max_nonpivot = max(
                (err for donor_idx, err in best_errors.items() if donor_idx != pivot_idx),
                default=0.0,
            )
            if (
                trial_score + 1e-9 < best_score
                and trial_max_nonpivot <= best_max_nonpivot + 1e-4
            ):
                best_score = trial_score
                best_coords = trial_coords
                best_errors = trial_errors
                improved = True

    if not improved:
        return False

    for atom_idx, vec in best_coords.items():
        conf.SetAtomPosition(
            atom_idx,
            Point3D(float(vec[0]), float(vec[1]), float(vec[2])),
        )
    return True


def _secondary_fragment_prefers_rigid_pose(
    mol,
    fragment: _HybridHaptoFragment,
    targeted_donor_indices: List[int],
) -> bool:
    """Return whether a CN=4 secondary-metal fragment should keep its rigid aligned pose."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None or not targeted_donor_indices:
        return False

    targeted_donor_set = set(targeted_donor_indices)
    if len(targeted_donor_set) >= 2:
        return True

    planar_like_atoms = 0
    for atom_idx in fragment.atom_indices:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            continue
        if atom.GetIsAromatic():
            planar_like_atoms += 1
            continue
        try:
            hyb = atom.GetHybridization()
        except Exception:
            continue
        if hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }:
            planar_like_atoms += 1

    return planar_like_atoms >= 4 and bool(targeted_donor_set & set(fragment.donor_atom_indices))


_HAPTO_SECONDARY_CN4_GEOMETRIES = {'SQ', 'TH', 'TET'}
_HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES = {'SQ', 'TET'}


def _relieve_secondary_oo_chelate_contacts(
    mol,
    fragment: _HybridHaptoFragment,
    metal_idx: int,
    donor_target_map: Dict[int, object],
) -> bool:
    """Rotate a rigid O,O chelate around the O-O axis to relieve non-donor contacts."""
    if not RDKIT_AVAILABLE or mol is None or fragment is None or not donor_target_map:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    donor_indices = [
        donor_idx for donor_idx in fragment.donor_atom_indices
        if donor_idx in donor_target_map and mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'O'
    ]
    if len(donor_indices) < 2:
        return False

    axis_start = np.array(conf.GetAtomPosition(donor_indices[0]), dtype=float)
    axis_end = np.array(conf.GetAtomPosition(donor_indices[1]), dtype=float)
    axis = axis_end - axis_start
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-12:
        return False
    axis = axis / axis_norm

    rotatable_atoms = [
        atom_idx for atom_idx in fragment.atom_indices
        if atom_idx not in donor_indices
        and mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
        and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
    ]
    if not rotatable_atoms:
        return False

    base_coords = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in fragment.atom_indices
    }
    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)

    def _score(coord_map: Dict[int, object]) -> float:
        return (
            _secondary_non_donor_contact_penalty(
                mol,
                metal_idx,
                metal_pos,
                fragment.atom_indices,
                donor_indices,
                coord_map,
            )
            + _fragment_donor_selectivity_penalty(
                mol,
                metal_idx,
                metal_pos,
                fragment.atom_indices,
                donor_indices,
                coord_map,
            )
            + _planar_fragment_metal_coplanarity_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                metal_pos,
                coord_map,
            )
            + _planar_fragment_donor_approach_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                metal_pos,
                coord_map,
            )
            + _fragment_environment_clash_penalty(
                mol,
                fragment.atom_indices,
                coord_map,
                ignore_atom_indices={metal_idx},
                primary_metal_idx=metal_idx,
            )
            + _fragment_bite_direction_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                coord_map,
                metal_pos,
            )
            + _planar_fragment_metal_coplanarity_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                metal_pos,
                coord_map,
            )
            + _planar_fragment_donor_approach_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                metal_pos,
                coord_map,
            )
        )

    best_coords = {atom_idx: vec.copy() for atom_idx, vec in base_coords.items()}
    best_score = _score(best_coords)
    improved = False
    pivot = 0.5 * (axis_start + axis_end)

    for angle_deg in (-30, -25, -20, -15, -10, -5, 5, 10, 15, 20, 25, 30):
        rot = _axis_angle_rotation_matrix(axis, math.radians(float(angle_deg)))
        trial_coords = {atom_idx: vec.copy() for atom_idx, vec in base_coords.items()}
        for atom_idx in rotatable_atoms:
            trial_coords[atom_idx] = (base_coords[atom_idx] - pivot) @ rot.T + pivot
        trial_score = _score(trial_coords)
        if trial_score + 1e-9 < best_score:
            best_score = trial_score
            best_coords = trial_coords
            improved = True

    if not improved:
        return False

    for atom_idx, vec in best_coords.items():
        conf.SetAtomPosition(atom_idx, Point3D(float(vec[0]), float(vec[1]), float(vec[2])))
    return True


def _restore_secondary_rigid_fragment_geometry(
    mol,
    fragment: _HybridHaptoFragment,
    metal_idx: int,
    donor_target_map: Dict[int, object],
    reference_coords: Dict[int, object],
) -> bool:
    """Restore a rigid secondary fragment against its saved pre-optimizer geometry."""
    if (
        not RDKIT_AVAILABLE
        or mol is None
        or fragment is None
        or not donor_target_map
        or not reference_coords
    ):
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    donor_indices = [
        donor_idx
        for donor_idx in fragment.donor_atom_indices
        if donor_idx in donor_target_map and donor_idx in reference_coords
    ]
    if len(donor_indices) < 2:
        return False

    atom_indices = [
        atom_idx for atom_idx in fragment.atom_indices
        if atom_idx in reference_coords
    ]
    if len(atom_indices) < 3:
        return False

    src = np.asarray([reference_coords[donor_idx] for donor_idx in donor_indices], dtype=float)
    tgt = np.asarray(
        [np.array(conf.GetAtomPosition(donor_idx), dtype=float) for donor_idx in donor_indices],
        dtype=float,
    )
    all_ref = np.asarray([reference_coords[atom_idx] for atom_idx in atom_indices], dtype=float)

    if len(donor_indices) == 2:
        rot = _rotation_matrix_from_vectors(src[1] - src[0], tgt[1] - tgt[0])
        src_mid = 0.5 * (src[0] + src[1])
        tgt_mid = 0.5 * (tgt[0] + tgt[1])
        aligned = (all_ref - src_mid) @ rot.T + tgt_mid
    else:
        src_centroid = src.mean(axis=0)
        tgt_centroid = tgt.mean(axis=0)
        covariance = (src - src_centroid).T @ (tgt - tgt_centroid)
        try:
            u, _s, vt = np.linalg.svd(covariance)
        except Exception:
            return False
        rot = vt.T @ u.T
        if float(np.linalg.det(rot)) < 0.0:
            vt[-1, :] *= -1.0
            rot = vt.T @ u.T
        aligned = (all_ref - src_centroid) @ rot.T + tgt_centroid

    donor_set = set(donor_indices)
    current_coords = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in atom_indices
    }
    rigid_coords = {
        atom_idx: aligned[idx].copy()
        for idx, atom_idx in enumerate(atom_indices)
    }
    for row_idx, donor_idx in enumerate(donor_indices):
        rigid_coords[donor_idx] = tgt[row_idx].copy()

    heavy_atoms = [
        atom_idx for atom_idx in atom_indices
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
        and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
    ]
    if len(heavy_atoms) < 3:
        return False

    ref_pair_targets: List[Tuple[int, int, float]] = []
    for i in range(len(heavy_atoms)):
        for j in range(i + 1, len(heavy_atoms)):
            atom_i = heavy_atoms[i]
            atom_j = heavy_atoms[j]
            ref_pair_targets.append(
                (
                    atom_i,
                    atom_j,
                    float(
                        np.linalg.norm(
                            np.asarray(reference_coords[atom_i], dtype=float)
                            - np.asarray(reference_coords[atom_j], dtype=float)
                        )
                    ),
                )
            )

    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
    planar_atoms: List[int] = []
    for atom_idx in heavy_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic():
            planar_atoms.append(atom_idx)
            continue
        try:
            hyb = atom.GetHybridization()
        except Exception:
            continue
        if hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }:
            planar_atoms.append(atom_idx)

    def _external_score(coord_map: Dict[int, object]) -> float:
        return (
            _secondary_non_donor_contact_penalty(
                mol,
                metal_idx,
                metal_pos,
                fragment.atom_indices,
                donor_indices,
                coord_map,
            )
            + _fragment_donor_selectivity_penalty(
                mol,
                metal_idx,
                metal_pos,
                fragment.atom_indices,
                donor_indices,
                coord_map,
            )
            + _fragment_environment_clash_penalty(
                mol,
                fragment.atom_indices,
                coord_map,
                ignore_atom_indices={metal_idx},
                primary_metal_idx=metal_idx,
            )
            + _fragment_bite_direction_penalty(
                mol,
                fragment.atom_indices,
                donor_indices,
                coord_map,
                metal_pos,
            )
        )

    def _shape_penalty(coord_map: Dict[int, object]) -> float:
        penalty = 0.0
        for atom_i, atom_j, target_len in ref_pair_targets:
            vec_i = np.asarray(coord_map.get(atom_i), dtype=float)
            vec_j = np.asarray(coord_map.get(atom_j), dtype=float)
            if vec_i.shape != (3,) or vec_j.shape != (3,):
                continue
            penalty += (float(np.linalg.norm(vec_i - vec_j)) - target_len) ** 2
        return penalty

    def _planarity_penalty(coord_map: Dict[int, object]) -> float:
        if len(planar_atoms) < 4:
            return 0.0
        pts = np.asarray([coord_map[atom_idx] for atom_idx in planar_atoms], dtype=float)
        if pts.shape[0] < 4:
            return 0.0
        centroid = pts.mean(axis=0)
        try:
            _u, _s, vt = np.linalg.svd(pts - centroid)
        except Exception:
            return 0.0
        normal = vt[-1]
        normal_norm = float(np.linalg.norm(normal))
        if normal_norm < 1e-12:
            return 0.0
        normal = normal / normal_norm
        rms = float(np.sqrt(np.mean(((pts - centroid) @ normal) ** 2)))
        return rms * rms

    def _total_score(coord_map: Dict[int, object]) -> float:
        return (
            _external_score(coord_map)
            + 6.0 * _shape_penalty(coord_map)
            + 28.0 * _planarity_penalty(coord_map)
        )

    current_score = _total_score(current_coords)
    best_coords = {atom_idx: vec.copy() for atom_idx, vec in rigid_coords.items()}
    best_score = _total_score(best_coords)

    if len(donor_indices) == 2:
        axis = tgt[1] - tgt[0]
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm > 1e-12:
            axis = axis / axis_norm
            pivot = 0.5 * (tgt[0] + tgt[1])
            rotatable_atoms = [atom_idx for atom_idx in atom_indices if atom_idx not in donor_set]
            for angle_deg in range(-180, 181, 15):
                if angle_deg == 0:
                    continue
                rot = _axis_angle_rotation_matrix(axis, math.radians(float(angle_deg)))
                trial_coords = {
                    atom_idx: vec.copy() for atom_idx, vec in rigid_coords.items()
                }
                for atom_idx in rotatable_atoms:
                    trial_coords[atom_idx] = (rigid_coords[atom_idx] - pivot) @ rot.T + pivot
                trial_score = _total_score(trial_coords)
                if trial_score + 1e-9 < best_score:
                    best_score = trial_score
                    best_coords = trial_coords

    if best_score + 1e-9 >= current_score:
        return False

    for row_idx, donor_idx in enumerate(donor_indices):
        best_coords[donor_idx] = tgt[row_idx].copy()
    for atom_idx, vec in best_coords.items():
        conf.SetAtomPosition(atom_idx, Point3D(float(vec[0]), float(vec[1]), float(vec[2])))
    return True


def _optimize_secondary_metal_module_local(
    mol,
    fragments: List[_HybridHaptoFragment],
    hapto_groups: List[Tuple[int, List[int]]],
    metal_idx: int,
    donor_indices: List[int],
    donor_to_fragment: Dict[int, _HybridHaptoFragment],
    geometry_code: Optional[str] = None,
    geometry_metal_pos=None,
    geometry_target_map: Optional[Dict[int, object]] = None,
) -> set:
    """Locally relax one hapto-coupled secondary-metal module with the hapto core fixed."""
    if not RDKIT_AVAILABLE or mol is None or not donor_indices:
        return set()
    try:
        import numpy as np
    except ImportError:
        return set()

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return set()

    metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    donor_set = set(donor_indices)
    module_atoms: set = {metal_idx} | donor_set
    movable_atoms: set = {metal_idx}
    donor_pair_targets: List[Tuple[int, int, float]] = []
    rigid_pair_targets: List[Tuple[int, int, float, float]] = []
    seen_donor_pairs: set = set()
    rigid_units: List[Tuple[int, ...]] = []
    seen_units: set = set()
    metal_repulsion_targets: List[Tuple[int, float]] = []
    atom_to_fragment_idx: Dict[int, int] = {}

    hapto_metals = {idx for idx, _grp in hapto_groups}
    hapto_atoms = {atom_idx for _idx, grp in hapto_groups for atom_idx in grp}
    for frag_idx, fragment in enumerate(fragments):
        for atom_idx in fragment.atom_indices:
            atom_to_fragment_idx[atom_idx] = frag_idx
    for atom in mol.GetAtoms():
        other_idx = atom.GetIdx()
        if other_idx == metal_idx or atom.GetSymbol() not in _METAL_SET:
            continue
        min_dist = max(
            3.05,
            _COVALENT_RADII.get(metal_sym, 1.35) + _COVALENT_RADII.get(atom.GetSymbol(), 1.35) + 0.45,
        )
        metal_repulsion_targets.append((other_idx, float(min_dist)))

    for donor_idx in donor_indices:
        fragment = donor_to_fragment.get(donor_idx)
        if fragment is None:
            continue
        module_atoms.update(fragment.atom_indices)
        fragment_donors = [
            idx for idx in fragment.donor_atom_indices
            if idx in donor_set
        ]
        for i in range(len(fragment_donors)):
            for j in range(i + 1, len(fragment_donors)):
                pair = tuple(sorted((fragment_donors[i], fragment_donors[j])))
                if pair in seen_donor_pairs:
                    continue
                pos_i = np.array(conf.GetAtomPosition(pair[0]), dtype=float)
                pos_j = np.array(conf.GetAtomPosition(pair[1]), dtype=float)
                donor_pair_targets.append((pair[0], pair[1], float(np.linalg.norm(pos_i - pos_j))))
                seen_donor_pairs.add(pair)

        if _secondary_module_fragment_is_movable(fragment, metal_idx):
            movable_atoms.update(fragment.atom_indices)
            unit = tuple(
                sorted(
                    atom_idx for atom_idx in fragment.atom_indices
                    if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
                    and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
                )
            )
            if len(unit) >= 3 and unit not in seen_units:
                rigid_units.append(unit)
                seen_units.add(unit)
            continue
        if fragment.use_scaffold_only:
            branch_info = _find_scaffold_secondary_branch(
                mol,
                fragment,
                hapto_groups,
                donor_idx,
            )
            if branch_info is not None:
                pivot_idx, branch_atoms = branch_info
                module_atoms.add(pivot_idx)
                module_atoms.update(branch_atoms)
                movable_atoms.update(branch_atoms)
                unit = tuple(
                    sorted(
                        atom_idx for atom_idx in ({pivot_idx} | set(branch_atoms))
                        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
                        and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
                    )
                )
                if len(unit) >= 3 and unit not in seen_units:
                    rigid_units.append(unit)
                    seen_units.add(unit)

    movable_atoms -= hapto_atoms
    movable_atoms -= (hapto_metals - {metal_idx})
    movable_atoms.add(metal_idx)
    movable_atoms &= module_atoms
    if len(movable_atoms) <= 1:
        return set()

    bonded_pairs: set = set()
    bond_targets: List[Tuple[int, int, float]] = []
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx not in module_atoms or end_idx not in module_atoms:
            continue
        bonded_pairs.add((min(begin_idx, end_idx), max(begin_idx, end_idx)))
        if begin_idx in movable_atoms or end_idx in movable_atoms:
            bond_targets.append((begin_idx, end_idx, _hybrid_bond_target_length(bond)))

    def _donor_target_profile(donor_idx: int) -> Tuple[float, float]:
        atom = mol.GetAtomWithIdx(donor_idx)
        donor_sym = atom.GetSymbol()
        target_len = float(_secondary_donor_target_length(mol, metal_sym, donor_idx))
        weight = 2.0 * float(_secondary_donor_fit_weight(mol, metal_sym, donor_idx))
        formal_charge = int(atom.GetFormalCharge())
        aromatic = bool(atom.GetIsAromatic())
        try:
            hyb = atom.GetHybridization()
        except Exception:
            hyb = None
        planar_like = aromatic or hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

        if donor_sym == 'N':
            if planar_like:
                weight += 0.35
            if formal_charge > 0:
                target_len += 0.05
                weight -= 0.45
            elif formal_charge < 0:
                target_len -= 0.04
                weight += 0.20
        elif donor_sym == 'O':
            if any(b.GetBondTypeAsDouble() >= 1.5 for b in atom.GetBonds()):
                weight += 0.20
            if any(b.GetIsConjugated() for b in atom.GetBonds()):
                weight += 0.15
            if formal_charge < 0:
                target_len -= 0.08
                weight += 0.70
            elif formal_charge > 0:
                target_len += 0.04
                weight -= 0.20
        elif donor_sym == 'C':
            if formal_charge < 0:
                target_len -= 0.06
                weight += 0.40
            elif formal_charge > 0:
                target_len += 0.05
                weight -= 0.30

        return target_len, max(weight, 0.85)

    ml_targets = [
        (
            donor_idx,
            *_donor_target_profile(donor_idx),
        )
        for donor_idx in donor_indices
    ]
    cn4_geometry_pair_targets: List[Tuple[int, int, float, float]] = []
    cn4_square_planar_normal = None
    if (
        geometry_code in _HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES
        and len(donor_indices) == 4
        and geometry_metal_pos is not None
        and geometry_target_map
    ):
        pair_weight = 0.55 if geometry_code == 'SQ' else 0.42
        for i in range(len(donor_indices)):
            pos_i = geometry_target_map.get(donor_indices[i])
            if pos_i is None:
                cn4_geometry_pair_targets = []
                break
            pos_i = np.asarray(pos_i, dtype=float)
            if pos_i.shape != (3,) or not np.all(np.isfinite(pos_i)):
                cn4_geometry_pair_targets = []
                break
            for j in range(i + 1, len(donor_indices)):
                pos_j = geometry_target_map.get(donor_indices[j])
                if pos_j is None:
                    cn4_geometry_pair_targets = []
                    break
                pos_j = np.asarray(pos_j, dtype=float)
                if pos_j.shape != (3,) or not np.all(np.isfinite(pos_j)):
                    cn4_geometry_pair_targets = []
                    break
                target_dist = float(np.linalg.norm(pos_i - pos_j))
                cn4_geometry_pair_targets.append(
                    (donor_indices[i], donor_indices[j], target_dist, pair_weight)
                )
            if not cn4_geometry_pair_targets and i > 0:
                break
        if geometry_code == 'SQ' and len(cn4_geometry_pair_targets) >= 2:
            target_vecs = []
            for donor_idx in donor_indices:
                target_pos = geometry_target_map.get(donor_idx)
                if target_pos is None:
                    target_vecs = []
                    break
                target_arr = np.asarray(target_pos, dtype=float) - np.asarray(geometry_metal_pos, dtype=float)
                if target_arr.shape != (3,) or not np.all(np.isfinite(target_arr)):
                    target_vecs = []
                    break
                target_vecs.append(target_arr)
            best_norm = 0.0
            for i in range(len(target_vecs)):
                for j in range(i + 1, len(target_vecs)):
                    normal = np.cross(target_vecs[i], target_vecs[j])
                    norm = float(np.linalg.norm(normal))
                    if norm > best_norm:
                        best_norm = norm
                        cn4_square_planar_normal = normal / norm

    initial_positions = {
        atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
        for atom_idx in movable_atoms
    }
    initial_positions.update(
        {
            atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
            for unit in rigid_units
            for atom_idx in unit
            if atom_idx not in initial_positions
        }
    )

    for unit in rigid_units:
        unit_atoms = list(unit)
        unit_planar = False
        planar_like_count = 0
        for atom_idx in unit_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            if atom.GetIsAromatic():
                planar_like_count += 1
                continue
            try:
                hyb = atom.GetHybridization()
            except Exception:
                continue
            if hyb in {
                Chem.rdchem.HybridizationType.SP,
                Chem.rdchem.HybridizationType.SP2,
            }:
                planar_like_count += 1
        unit_planar = planar_like_count >= 4

        if len(unit_atoms) <= 6 or (unit_planar and len(unit_atoms) <= 14):
            pair_iter = [
                (unit_atoms[i], unit_atoms[j])
                for i in range(len(unit_atoms))
                for j in range(i + 1, len(unit_atoms))
            ]
        else:
            pair_iter = []
            unit_set = set(unit_atoms)
            for atom_idx in unit_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in unit_set and nbr_idx > atom_idx:
                        pair_iter.append((atom_idx, nbr_idx))
                    if nbr_idx not in unit_set:
                        continue
                    for nbr2 in nbr.GetNeighbors():
                        nbr2_idx = nbr2.GetIdx()
                        if nbr2_idx in unit_set and nbr2_idx > atom_idx and nbr2_idx != atom_idx:
                            pair = (atom_idx, nbr2_idx)
                            if pair not in pair_iter:
                                pair_iter.append(pair)
        seen_pairs = set()
        for atom_i, atom_j in pair_iter:
            pair = (min(atom_i, atom_j), max(atom_i, atom_j))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            target = float(np.linalg.norm(initial_positions[pair[0]] - initial_positions[pair[1]]))
            bond = mol.GetBondBetweenAtoms(pair[0], pair[1])
            if bond is not None:
                weight = 1.60
            else:
                weight = 0.85 if unit_planar else 0.55
            rigid_pair_targets.append((pair[0], pair[1], target, weight))

    def _is_planar_like_atom(atom) -> bool:
        if atom is None or atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            return False
        if atom.GetIsAromatic():
            return True
        try:
            hyb = atom.GetHybridization()
        except Exception:
            return False
        return hyb in {
            Chem.rdchem.HybridizationType.SP,
            Chem.rdchem.HybridizationType.SP2,
        }

    planar_units: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = []
    for unit in rigid_units:
        unit_set = set(unit)
        movable_overlap = unit_set & movable_atoms
        if len(movable_overlap) < 2:
            continue
        donor_overlap = [
            donor_idx for donor_idx in donor_indices
            if donor_idx in unit_set and _is_planar_like_atom(mol.GetAtomWithIdx(donor_idx))
        ]
        if not donor_overlap:
            continue
        planar_atoms = tuple(
            atom_idx for atom_idx in unit
            if _is_planar_like_atom(mol.GetAtomWithIdx(atom_idx))
        )
        if len(planar_atoms) < 4:
            continue
        pts = np.asarray([initial_positions[atom_idx] for atom_idx in planar_atoms], dtype=float)
        centroid = pts.mean(axis=0)
        try:
            _u, _s, vt = np.linalg.svd(pts - centroid)
        except Exception:
            continue
        normal = vt[-1]
        rms_offset = float(np.sqrt(np.mean(((pts - centroid) @ normal) ** 2)))
        if rms_offset <= 0.18:
            planar_units.append((tuple(sorted(unit)), planar_atoms))

    metal_coplanar_units: List[Tuple[Tuple[int, ...], float]] = []
    seen_metal_coplanar_units: set = set()
    for unit in rigid_units:
        unit_set = set(unit)
        donor_overlap = [
            donor_idx for donor_idx in donor_indices
            if donor_idx in unit_set and _is_planar_like_atom(mol.GetAtomWithIdx(donor_idx))
        ]
        if not donor_overlap:
            continue
        planar_atoms = tuple(
            atom_idx for atom_idx in unit
            if _is_planar_like_atom(mol.GetAtomWithIdx(atom_idx))
        )
        if len(planar_atoms) < 4:
            continue
        planar_key = tuple(sorted(planar_atoms))
        if planar_key in seen_metal_coplanar_units:
            continue
        seen_metal_coplanar_units.add(planar_key)
        weight = 34.0 + 10.0 * len(donor_overlap)
        if any(mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'N' for donor_idx in donor_overlap):
            weight += 12.0
        if any(mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'O' for donor_idx in donor_overlap):
            weight += 12.0
        if any(mol.GetAtomWithIdx(donor_idx).GetIsAromatic() for donor_idx in donor_overlap):
            weight += 10.0
        metal_coplanar_units.append((planar_key, weight))

    def _gp(atom_idx: int):
        pos = conf.GetAtomPosition(atom_idx)
        return np.array([pos.x, pos.y, pos.z], dtype=float)

    def _sp(atom_idx: int, vec):
        conf.SetAtomPosition(atom_idx, Point3D(float(vec[0]), float(vec[1]), float(vec[2])))

    def _score() -> float:
        metal_pos = _gp(metal_idx)
        coord_map = {
            atom_idx: _gp(atom_idx)
            for atom_idx in module_atoms
            if atom_idx != metal_idx
        }
        score = 0.0
        for begin_idx, end_idx, target_len in bond_targets:
            dist = float(np.linalg.norm(_gp(begin_idx) - _gp(end_idx)))
            score += 0.45 * (dist - target_len) ** 2
        for donor_idx, target_len, weight in ml_targets:
            dist = float(np.linalg.norm(_gp(donor_idx) - metal_pos))
            score += weight * (dist - target_len) ** 2
        for donor_i, donor_j, target_len in donor_pair_targets:
            dist = float(np.linalg.norm(_gp(donor_i) - _gp(donor_j)))
            score += 0.20 * (dist - target_len) ** 2
        for donor_i, donor_j, target_len, weight in cn4_geometry_pair_targets:
            dist = float(np.linalg.norm(_gp(donor_i) - _gp(donor_j)))
            score += weight * (dist - target_len) ** 2
        if cn4_square_planar_normal is not None:
            normal = np.asarray(cn4_square_planar_normal, dtype=float)
            for donor_idx in donor_indices:
                offset = float(np.dot(_gp(donor_idx) - metal_pos, normal))
                score += 0.35 * offset * offset
        for planar_atoms, weight in metal_coplanar_units:
            pts = []
            for atom_idx in planar_atoms:
                arr = coord_map.get(atom_idx)
                if arr is None:
                    arr = _gp(atom_idx)
                pts.append(np.asarray(arr, dtype=float))
            if len(pts) < 4:
                continue
            pts_arr = np.asarray(pts, dtype=float)
            centroid = pts_arr.mean(axis=0)
            try:
                _u, _s, vt = np.linalg.svd(pts_arr - centroid)
            except Exception:
                continue
            normal = vt[-1]
            normal_norm = float(np.linalg.norm(normal))
            if normal_norm < 1e-12:
                continue
            normal = normal / normal_norm
            plane_offset = float(np.dot(metal_pos - centroid, normal))
            score += weight * plane_offset * plane_offset
            if abs(plane_offset) > 0.28:
                score += 75.0 * (abs(plane_offset) - 0.28) ** 2
        seen_fragment_planes: set = set()
        for donor_idx in donor_indices:
            fragment = donor_to_fragment.get(donor_idx)
            if fragment is None:
                continue
            fragment_key = tuple(sorted(fragment.atom_indices))
            if fragment_key in seen_fragment_planes:
                continue
            seen_fragment_planes.add(fragment_key)
            fragment_donors = [
                idx for idx in fragment.donor_atom_indices
                if idx in donor_set
            ]
            if not fragment_donors:
                continue
            score += _planar_fragment_donor_approach_penalty(
                mol,
                fragment.atom_indices,
                fragment_donors,
                metal_pos,
                coord_map,
            )
        for atom_i, atom_j, target_len, weight in rigid_pair_targets:
            dist = float(np.linalg.norm(_gp(atom_i) - _gp(atom_j)))
            score += weight * (dist - target_len) ** 2
        for other_metal_idx, min_dist in metal_repulsion_targets:
            dist = float(np.linalg.norm(_gp(other_metal_idx) - metal_pos))
            if dist < min_dist:
                score += 16.0 * (min_dist - dist) ** 2
        module_atom_list = sorted(module_atoms - {metal_idx})
        for idx_i, atom_i in enumerate(module_atom_list):
            frag_i = atom_to_fragment_idx.get(atom_i, -1)
            sym_i = mol.GetAtomWithIdx(atom_i).GetSymbol()
            if sym_i == 'H' or sym_i in _METAL_SET:
                continue
            pos_i = _gp(atom_i)
            for atom_j in module_atom_list[idx_i + 1:]:
                if (min(atom_i, atom_j), max(atom_i, atom_j)) in bonded_pairs:
                    continue
                frag_j = atom_to_fragment_idx.get(atom_j, -1)
                if frag_i == frag_j:
                    continue
                sym_j = mol.GetAtomWithIdx(atom_j).GetSymbol()
                if sym_j == 'H' or sym_j in _METAL_SET:
                    continue
                pos_j = _gp(atom_j)
                dist = float(np.linalg.norm(pos_i - pos_j))
                min_dist = max(
                    1.78,
                    _COVALENT_RADII.get(sym_i, 0.76) + _COVALENT_RADII.get(sym_j, 0.76) + 0.36,
                )
                if sym_i in {'N', 'O', 'S', 'P'} or sym_j in {'N', 'O', 'S', 'P'}:
                    min_dist = max(min_dist, 1.90)
                if dist < min_dist:
                    score += 12.0 * (min_dist - dist) ** 2
        score += _secondary_non_donor_contact_penalty(
            mol,
            metal_idx,
            metal_pos,
            sorted(module_atoms - {metal_idx}),
            donor_indices,
            coord_map,
        )
        score += _secondary_fragment_donor_selectivity_penalty(
            mol,
            metal_idx,
            metal_pos,
            donor_indices,
            donor_to_fragment,
            coord_map,
        )
        return score

    def _cn4_geometry_deviation() -> float:
        if len(donor_indices) != 4 or geometry_code not in _HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES:
            return 0.0
        deviation = 0.0
        metal_pos = _gp(metal_idx)
        for donor_i, donor_j, target_len, _weight in cn4_geometry_pair_targets:
            dist = float(np.linalg.norm(_gp(donor_i) - _gp(donor_j)))
            deviation += (dist - target_len) ** 2
        if cn4_square_planar_normal is not None:
            normal = np.asarray(cn4_square_planar_normal, dtype=float)
            for donor_idx in donor_indices:
                offset = float(np.dot(_gp(donor_idx) - metal_pos, normal))
                deviation += 0.5 * offset * offset
        return deviation

    best_positions = {atom_idx: vec.copy() for atom_idx, vec in initial_positions.items()}
    start_score = _score()
    best_score = start_score
    start_geometry_deviation = _cn4_geometry_deviation()
    best_geometry_deviation = start_geometry_deviation
    rng = np.random.default_rng(int(1009 * (metal_idx + 1) + len(module_atoms)))

    alpha_atoms: set = set()
    beta_atoms: set = set()
    for donor_idx in donor_indices:
        atom = mol.GetAtomWithIdx(donor_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in module_atoms and nbr_idx not in donor_set and nbr.GetAtomicNum() > 1:
                alpha_atoms.add(nbr_idx)
    for atom_idx in list(alpha_atoms):
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if (
                nbr_idx in module_atoms
                and nbr_idx not in donor_set
                and nbr_idx not in alpha_atoms
                and nbr.GetAtomicNum() > 1
            ):
                beta_atoms.add(nbr_idx)

    for _pass in range(120):
        displacements: Dict[int, np.ndarray] = {
            atom_idx: np.zeros(3, dtype=float) for atom_idx in movable_atoms
        }
        active = False

        for begin_idx, end_idx, target_len in bond_targets:
            begin_pos = _gp(begin_idx)
            end_pos = _gp(end_idx)
            diff = end_pos - begin_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.01:
                continue
            active = True
            strength = 0.42 if ((begin_idx in movable_atoms) ^ (end_idx in movable_atoms)) else 0.24
            move = strength * delta * unit
            if begin_idx in movable_atoms and end_idx in movable_atoms:
                displacements[begin_idx] = displacements[begin_idx] + 0.5 * move
                displacements[end_idx] = displacements[end_idx] - 0.5 * move
            elif begin_idx in movable_atoms:
                displacements[begin_idx] = displacements[begin_idx] + move
            elif end_idx in movable_atoms:
                displacements[end_idx] = displacements[end_idx] - move

        metal_pos = _gp(metal_idx)
        for donor_idx, target_len, _weight in ml_targets:
            donor_pos = _gp(donor_idx)
            diff = donor_pos - metal_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.01:
                continue
            active = True
            strength = 0.40 if donor_idx not in movable_atoms else 0.22
            move = strength * delta * unit
            if metal_idx in movable_atoms and donor_idx in movable_atoms:
                displacements[metal_idx] = displacements[metal_idx] + 0.55 * move
                displacements[donor_idx] = displacements[donor_idx] - 0.45 * move
            elif metal_idx in movable_atoms:
                displacements[metal_idx] = displacements[metal_idx] + move
            elif donor_idx in movable_atoms:
                displacements[donor_idx] = displacements[donor_idx] - move

        for donor_i, donor_j, target_len in donor_pair_targets:
            pos_i = _gp(donor_i)
            pos_j = _gp(donor_j)
            diff = pos_j - pos_i
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.015:
                continue
            active = True
            move = 0.12 * delta * unit
            if donor_i in movable_atoms and donor_j in movable_atoms:
                displacements[donor_i] = displacements[donor_i] + 0.5 * move
                displacements[donor_j] = displacements[donor_j] - 0.5 * move
            elif donor_i in movable_atoms:
                displacements[donor_i] = displacements[donor_i] + move
            elif donor_j in movable_atoms:
                displacements[donor_j] = displacements[donor_j] - move

        for donor_i, donor_j, target_len, weight in cn4_geometry_pair_targets:
            pos_i = _gp(donor_i)
            pos_j = _gp(donor_j)
            diff = pos_j - pos_i
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.01:
                continue
            active = True
            move = min(0.18, 0.08 + 0.10 * weight) * delta * unit
            if donor_i in movable_atoms and donor_j in movable_atoms:
                displacements[donor_i] = displacements[donor_i] + 0.5 * move
                displacements[donor_j] = displacements[donor_j] - 0.5 * move
            elif donor_i in movable_atoms:
                displacements[donor_i] = displacements[donor_i] + move
            elif donor_j in movable_atoms:
                displacements[donor_j] = displacements[donor_j] - move

        if cn4_square_planar_normal is not None:
            normal = np.asarray(cn4_square_planar_normal, dtype=float)
            metal_pos = _gp(metal_idx)
            for donor_idx in donor_indices:
                donor_pos = _gp(donor_idx)
                offset = float(np.dot(donor_pos - metal_pos, normal))
                if abs(offset) < 0.01:
                    continue
                active = True
                move = 0.16 * offset * normal
                if donor_idx in movable_atoms and metal_idx in movable_atoms:
                    displacements[donor_idx] = displacements[donor_idx] - 0.55 * move
                    displacements[metal_idx] = displacements[metal_idx] + 0.45 * move
                elif donor_idx in movable_atoms:
                    displacements[donor_idx] = displacements[donor_idx] - move
                elif metal_idx in movable_atoms:
                    displacements[metal_idx] = displacements[metal_idx] + move

        for atom_i, atom_j, target_len, weight in rigid_pair_targets:
            pos_i = _gp(atom_i)
            pos_j = _gp(atom_j)
            diff = pos_j - pos_i
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target_len
            if abs(delta) < 0.008:
                continue
            active = True
            move = min(0.18, 0.10 + 0.04 * weight) * delta * unit
            if atom_i in movable_atoms and atom_j in movable_atoms:
                displacements[atom_i] = displacements[atom_i] + 0.5 * move
                displacements[atom_j] = displacements[atom_j] - 0.5 * move
            elif atom_i in movable_atoms:
                displacements[atom_i] = displacements[atom_i] + move
            elif atom_j in movable_atoms:
                displacements[atom_j] = displacements[atom_j] - move

        metal_pos = _gp(metal_idx)
        for other_metal_idx, min_dist in metal_repulsion_targets:
            other_pos = _gp(other_metal_idx)
            diff = metal_pos - other_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            if dist >= min_dist:
                continue
            active = True
            unit = diff / max(dist, 1e-12)
            gap = min_dist - dist
            displacements[metal_idx] = displacements[metal_idx] + 0.70 * gap * unit

        for atom_idx in sorted(module_atoms - donor_set - {metal_idx}):
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            atom_pos = _gp(atom_idx)
            diff = atom_pos - metal_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            min_allowed = max(1.70, 0.92 * float(_get_ml_bond_length(metal_sym, atom.GetSymbol())))
            if atom_idx in alpha_atoms:
                if _is_planar_like_atom(atom):
                    min_allowed = max(min_allowed, 2.18)
                else:
                    min_allowed = max(min_allowed, 2.05)
            elif atom_idx in beta_atoms:
                if _is_planar_like_atom(atom):
                    min_allowed = max(min_allowed, 1.96)
                else:
                    min_allowed = max(min_allowed, 1.85)
            if dist >= min_allowed:
                continue
            active = True
            gap = min_allowed - dist
            move = 0.60 * gap * unit
            if metal_idx in movable_atoms and atom_idx in movable_atoms:
                displacements[metal_idx] = displacements[metal_idx] - 0.65 * move
                displacements[atom_idx] = displacements[atom_idx] + 0.35 * move
            elif metal_idx in movable_atoms:
                displacements[metal_idx] = displacements[metal_idx] - move
            elif atom_idx in movable_atoms:
                displacements[atom_idx] = displacements[atom_idx] + move

        for atom_i in sorted(module_atoms):
            for atom_j in sorted(module_atoms):
                if atom_j <= atom_i:
                    continue
                if (min(atom_i, atom_j), max(atom_i, atom_j)) in bonded_pairs:
                    continue
                if atom_i not in movable_atoms and atom_j not in movable_atoms:
                    continue
                pos_i = _gp(atom_i)
                pos_j = _gp(atom_j)
                diff = pos_j - pos_i
                dist = float(np.linalg.norm(diff))
                if dist < 1e-8:
                    diff = rng.standard_normal(3)
                    dist = float(np.linalg.norm(diff))
                unit = diff / max(dist, 1e-12)
                sym_i = mol.GetAtomWithIdx(atom_i).GetSymbol()
                sym_j = mol.GetAtomWithIdx(atom_j).GetSymbol()
                if sym_i in _METAL_SET or sym_j in _METAL_SET:
                    min_dist = 2.0
                elif sym_i == 'H' and sym_j == 'H':
                    min_dist = 1.5
                elif sym_i == 'H' or sym_j == 'H':
                    min_dist = 1.0
                else:
                    min_dist = 1.18
                    frag_i = atom_to_fragment_idx.get(atom_i, -1)
                    frag_j = atom_to_fragment_idx.get(atom_j, -1)
                    if frag_i != frag_j:
                        min_dist = max(
                            min_dist,
                            _COVALENT_RADII.get(sym_i, 0.76) + _COVALENT_RADII.get(sym_j, 0.76) + 0.36,
                            1.78,
                        )
                        if sym_i in {'N', 'O', 'S', 'P'} or sym_j in {'N', 'O', 'S', 'P'}:
                            min_dist = max(min_dist, 1.90)
                if dist >= min_dist:
                    continue
                active = True
                gap = min_dist - dist
                move_scale = 0.38 if atom_to_fragment_idx.get(atom_i, -1) != atom_to_fragment_idx.get(atom_j, -1) else 0.18
                move = move_scale * gap * unit
                if atom_i in movable_atoms and atom_j in movable_atoms:
                    displacements[atom_i] = displacements[atom_i] - 0.5 * move
                    displacements[atom_j] = displacements[atom_j] + 0.5 * move
                elif atom_i in movable_atoms:
                    displacements[atom_i] = displacements[atom_i] - move
                elif atom_j in movable_atoms:
                    displacements[atom_j] = displacements[atom_j] + move

        if not active:
            break

        capped_displacements: Dict[int, np.ndarray] = {}
        for atom_idx, disp in displacements.items():
            norm = float(np.linalg.norm(disp))
            if norm < 1e-8:
                continue
            if norm > 0.22:
                disp = disp * (0.22 / norm)
            capped_displacements[atom_idx] = disp

        for unit, planar_atoms in planar_units:
            if not any(atom_idx in capped_displacements for atom_idx in unit):
                continue
            unit_positions = {
                atom_idx: _gp(atom_idx)
                for atom_idx in unit
            }
            unit_displacements = {
                atom_idx: capped_displacements.get(atom_idx, np.zeros(3, dtype=float))
                for atom_idx in unit
            }
            rigid_projected = _project_displacements_to_rigid_body(
                unit_positions,
                unit_displacements,
                list(unit),
            )
            for atom_idx in unit:
                if atom_idx in capped_displacements and atom_idx in rigid_projected:
                    capped_displacements[atom_idx] = np.asarray(rigid_projected[atom_idx], dtype=float)

        for atom_idx, disp in capped_displacements.items():
            _sp(atom_idx, _gp(atom_idx) + disp)

        trial_score = _score()
        if trial_score + 1e-9 < best_score:
            best_score = trial_score
            best_geometry_deviation = _cn4_geometry_deviation()
            best_positions = {
                atom_idx: _gp(atom_idx).copy()
                for atom_idx in movable_atoms
            }

    if best_score + 1e-6 >= start_score:
        for atom_idx, vec in initial_positions.items():
            _sp(atom_idx, vec)
        return set()

    if len(donor_indices) == 4 and geometry_code in _HAPTO_SECONDARY_CN4_OPTIMIZER_GEOMETRIES:
        allowed_geometry_deviation = max(
            start_geometry_deviation + 0.03,
            1.35 * start_geometry_deviation + 1e-6,
        )
        if best_geometry_deviation > allowed_geometry_deviation:
            for atom_idx, vec in initial_positions.items():
                _sp(atom_idx, vec)
            logger.info(
                "Hybrid hapto rejected local secondary metal module optimization for %s%d: CN=4 %s geometry deviation %.3f -> %.3f",
                metal_sym,
                metal_idx,
                geometry_code,
                start_geometry_deviation,
                best_geometry_deviation,
            )
            return set()

    for atom_idx, vec in best_positions.items():
        _sp(atom_idx, vec)

    logger.info(
        "Hybrid hapto locally optimized secondary metal module %s%d: %.3f -> %.3f",
        metal_sym,
        metal_idx,
        start_score,
        best_score,
    )
    return set(movable_atoms)


def _assemble_secondary_metal_coordination_modules(
    mol,
    decomposition: _HybridHaptoDecomposition,
    hapto_groups: List[Tuple[int, List[int]]],
    fit_variant_plan: Optional[Dict[int, int]] = None,
) -> Tuple[set, set, set]:
    """Build secondary metal coordination spheres from donor fragments."""
    if not RDKIT_AVAILABLE or mol is None or decomposition is None or not hapto_groups:
        return set(), set(), set()
    try:
        import numpy as np
    except ImportError:
        return set(), set(), set()

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return set(), set(), set()

    hapto_metals = {metal_idx for metal_idx, _grp in hapto_groups}
    placed_secondary: set = set()
    module_atoms: set = set()
    relaxable_atoms: set = set()

    donor_to_fragment: Dict[int, _HybridHaptoFragment] = {}
    embedded_fragment_cache: Dict[Tuple[int, ...], object] = {}
    for fragment in decomposition.fragments:
        for donor_idx in fragment.donor_atom_indices:
            donor_to_fragment[donor_idx] = fragment

    for atom in mol.GetAtoms():
        metal_idx = atom.GetIdx()
        metal_sym = atom.GetSymbol()
        if metal_sym not in _METAL_SET or metal_idx in hapto_metals:
            continue
        variant_rank = max(0, int((fit_variant_plan or {}).get(metal_idx, 0)))

        donor_indices = [
            nbr.GetIdx()
            for nbr in atom.GetNeighbors()
            if nbr.GetSymbol() not in _METAL_SET and nbr.GetAtomicNum() > 1
        ]
        if len(donor_indices) == 0:
            continue
        if len(donor_indices) == 1:
            # Single-donor fallback: place metal along donor→COM direction
            donor_idx = donor_indices[0]
            donor_pos = conf.GetAtomPosition(donor_idx)
            dp = np.array([donor_pos.x, donor_pos.y, donor_pos.z])
            all_pos = []
            for a2 in mol.GetAtoms():
                if a2.GetSymbol() not in _METAL_SET:
                    p2 = conf.GetAtomPosition(a2.GetIdx())
                    all_pos.append(np.array([p2.x, p2.y, p2.z]))
            com = np.mean(all_pos, axis=0) if all_pos else dp
            direction = dp - com
            norm = np.linalg.norm(direction)
            if norm < 1e-6:
                direction = np.array([1.0, 0.0, 0.0])
                norm = 1.0
            direction = direction / norm
            try:
                bond_len = _secondary_donor_target_length(mol, metal_sym, donor_idx)
            except Exception:
                bond_len = 2.1
            metal_pos = dp + direction * bond_len
            conf.SetAtomPosition(metal_idx, Point3D(*metal_pos))
            placed_secondary.add(metal_idx)
            module_atoms.add(metal_idx)
            module_atoms.add(donor_idx)
            continue

        donor_positions = []
        donor_symbols = []
        donor_target_lengths = []
        donor_fit_weights = []
        heterodonor_count = 0
        carbon_donor_count = 0
        for donor_idx in donor_indices:
            pos = conf.GetAtomPosition(donor_idx)
            donor_positions.append(np.array([pos.x, pos.y, pos.z], dtype=float))
            donor_sym = mol.GetAtomWithIdx(donor_idx).GetSymbol()
            donor_symbols.append(donor_sym)
            donor_target_lengths.append(_secondary_donor_target_length(mol, metal_sym, donor_idx))
            donor_fit_weights.append(_secondary_donor_fit_weight(mol, metal_sym, donor_idx))
            if donor_sym == 'C':
                carbon_donor_count += 1
            elif donor_sym in {'N', 'O', 'S', 'P'}:
                heterodonor_count += 1

        movable_fragments: List[_HybridHaptoFragment] = []
        seen_fragments: set = set()
        anchored_slots: List[int] = []
        supplemental_fit_slots: List[int] = []
        bite_distance_constraints: List[Tuple[int, int, float]] = []
        for slot_idx, donor_idx in enumerate(donor_indices):
            fragment = donor_to_fragment.get(donor_idx)
            donor_sym = mol.GetAtomWithIdx(donor_idx).GetSymbol()
            if fragment is None:
                anchored_slots.append(slot_idx)
                continue
            if heterodonor_count and carbon_donor_count and donor_sym == 'C':
                if fragment.bridging_donor_indices or fragment.use_scaffold_only:
                    donor_fit_weights[slot_idx] *= 0.72
                elif len(fragment.hapto_group_ids) > 0:
                    donor_fit_weights[slot_idx] *= 0.82
            if donor_sym in {'N', 'O', 'S', 'P'}:
                if len(fragment.donor_atom_indices) >= 2:
                    donor_fit_weights[slot_idx] *= 1.12
                if any(
                    mol.GetAtomWithIdx(atom_idx).GetIsAromatic()
                    or mol.GetAtomWithIdx(atom_idx).GetHybridization() in {
                        Chem.rdchem.HybridizationType.SP,
                        Chem.rdchem.HybridizationType.SP2,
                    }
                    for atom_idx in fragment.atom_indices
                    if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1
                    and mol.GetAtomWithIdx(atom_idx).GetSymbol() not in _METAL_SET
                ):
                    donor_fit_weights[slot_idx] *= 1.08
            if _secondary_module_fragment_is_movable(fragment, metal_idx):
                frag_key = tuple(fragment.atom_indices)
                if frag_key not in seen_fragments:
                    movable_fragments.append(fragment)
                    seen_fragments.add(frag_key)
                supplemental_fit_slots.append(slot_idx)
                continue
            if (
                fragment.use_scaffold_only
                and _find_scaffold_secondary_branch(
                    mol,
                    fragment,
                    hapto_groups,
                    donor_idx,
                ) is not None
            ):
                supplemental_fit_slots.append(slot_idx)
                continue
            anchored_slots.append(slot_idx)

        fit_slot_indices: Optional[List[int]]
        if anchored_slots:
            fit_slot_indices = list(anchored_slots)
            for slot_idx in supplemental_fit_slots:
                if slot_idx not in fit_slot_indices:
                    fit_slot_indices.append(slot_idx)
                if len(fit_slot_indices) >= 2:
                    break
        else:
            fit_slot_indices = None

        has_multidentate_movable_fragment = any(
            len([
                donor_idx
                for donor_idx in fragment.donor_atom_indices
                if donor_idx in donor_indices
            ]) >= 2
            for fragment in movable_fragments
        )
        has_secondary_bridge_branch = any(
            fragment.use_scaffold_only and metal_idx in fragment.metal_neighbor_indices
            for fragment in decomposition.fragments
        )
        if (
            len(donor_indices) == 4
            and has_multidentate_movable_fragment
            and has_secondary_bridge_branch
        ):
            fit_slot_indices = None

        for fragment in movable_fragments:
            frag_key = tuple(fragment.atom_indices)
            embedded_fragment = embedded_fragment_cache.get(frag_key)
            if embedded_fragment is None:
                embedded_fragment = _embed_hybrid_fragment(fragment.fragment_mol)
                if embedded_fragment is not None:
                    embedded_fragment_cache[frag_key] = embedded_fragment
            if embedded_fragment is None:
                continue
            try:
                frag_conf = embedded_fragment.GetConformer()
            except Exception:
                continue
            frag_donors = [
                donor_idx for donor_idx in fragment.donor_atom_indices
                if donor_idx in donor_indices
            ]
            if len(frag_donors) < 2:
                continue
            donor_slots = {
                donor_idx: donor_indices.index(donor_idx)
                for donor_idx in frag_donors
            }
            for i in range(len(frag_donors)):
                for j in range(i + 1, len(frag_donors)):
                    donor_i = frag_donors[i]
                    donor_j = frag_donors[j]
                    frag_i = fragment.original_to_fragment.get(donor_i)
                    frag_j = fragment.original_to_fragment.get(donor_j)
                    if frag_i is None or frag_j is None:
                        continue
                    pos_i = frag_conf.GetAtomPosition(frag_i)
                    pos_j = frag_conf.GetAtomPosition(frag_j)
                    expected_dist = math.sqrt(
                        (pos_i.x - pos_j.x) ** 2
                        + (pos_i.y - pos_j.y) ** 2
                        + (pos_i.z - pos_j.z) ** 2
                    )
                    bite_distance_constraints.append(
                        (donor_slots[donor_i], donor_slots[donor_j], float(expected_dist))
                    )

        fit_candidates = _enumerate_secondary_metal_geometry_fits(
            donor_positions,
            donor_symbols,
            metal_sym,
            constrained_indices=fit_slot_indices,
            bite_distance_constraints=bite_distance_constraints,
            donor_target_lengths=donor_target_lengths,
            donor_fit_weights=donor_fit_weights,
            max_candidates=max(variant_rank + 1, 4),
        )
        fit = fit_candidates[min(variant_rank, len(fit_candidates) - 1)] if fit_candidates else None
        if fit is None:
            fallback_candidates = _enumerate_secondary_metal_position_fits(
                donor_positions,
                donor_symbols,
                metal_sym,
                donor_target_lengths=donor_target_lengths,
                donor_fit_weights=donor_fit_weights,
                max_candidates=max(variant_rank + 1, 4),
            )
            fallback = fallback_candidates[min(variant_rank, len(fallback_candidates) - 1)] if fallback_candidates else None
            if fallback is None:
                continue
            metal_pos, geom_code, _fit_score = fallback
            conf.SetAtomPosition(
                metal_idx,
                Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
            )
            placed_secondary.add(metal_idx)
            module_atoms.add(metal_idx)
            logger.info(
                "Hybrid hapto built secondary metal %s%d via fallback %s fit to %d donor(s)",
                metal_sym,
                metal_idx,
                geom_code,
                len(donor_indices),
            )
            continue

        metal_pos, geom_code, target_positions, _fit_score = fit
        donor_target_map: Dict[int, np.ndarray] = {
            donor_idx: np.asarray(target_positions[slot_idx], dtype=float)
            for slot_idx, donor_idx in enumerate(donor_indices)
        }
        conf.SetAtomPosition(
            metal_idx,
            Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
        )

        reoriented_scaffold_branches = 0
        for fragment in decomposition.fragments:
            if not fragment.use_scaffold_only or metal_idx not in fragment.metal_neighbor_indices:
                continue
            for donor_idx in fragment.donor_atom_indices:
                if donor_idx not in donor_target_map:
                    continue
                branch_info = _find_scaffold_secondary_branch(
                    mol,
                    fragment,
                    hapto_groups,
                    donor_idx,
                )
                if branch_info is None:
                    continue
                try:
                    if _reorient_scaffold_secondary_branch(
                        mol,
                        fragment,
                        hapto_groups,
                        metal_idx,
                        donor_idx,
                        donor_target_map[donor_idx],
                    ):
                        reoriented_scaffold_branches += 1
                        pivot_idx, branch_atoms = branch_info
                        module_atoms.add(pivot_idx)
                        module_atoms.update(branch_atoms)
                        relaxable_atoms.update(branch_atoms)
                except Exception:
                    continue

        if reoriented_scaffold_branches:
            donor_positions = []
            for donor_idx in donor_indices:
                pos = conf.GetAtomPosition(donor_idx)
                donor_positions.append(np.array([pos.x, pos.y, pos.z], dtype=float))
            refit_candidates = _enumerate_secondary_metal_geometry_fits(
                donor_positions,
                donor_symbols,
                metal_sym,
                constrained_indices=fit_slot_indices,
                bite_distance_constraints=bite_distance_constraints,
                donor_target_lengths=donor_target_lengths,
                donor_fit_weights=donor_fit_weights,
                max_candidates=max(variant_rank + 1, 4),
            )
            refit = refit_candidates[min(variant_rank, len(refit_candidates) - 1)] if refit_candidates else None
            if refit is not None:
                metal_pos, geom_code, target_positions, _fit_score = refit
                donor_target_map = {
                    donor_idx: np.asarray(target_positions[slot_idx], dtype=float)
                    for slot_idx, donor_idx in enumerate(donor_indices)
                }
                conf.SetAtomPosition(
                    metal_idx,
                    Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
                )

        moved_fragments = 0
        has_multidentate_movable_fragment = False
        rigid_reference_fragments: List[Tuple[_HybridHaptoFragment, Dict[int, np.ndarray]]] = []
        oo_relief_fragments: List[_HybridHaptoFragment] = []
        for fragment in movable_fragments:
            frag_key = tuple(fragment.atom_indices)
            embedded_fragment = embedded_fragment_cache.get(frag_key)
            if embedded_fragment is None:
                embedded_fragment = _embed_hybrid_fragment(fragment.fragment_mol)
                if embedded_fragment is not None:
                    embedded_fragment_cache[frag_key] = embedded_fragment
            if embedded_fragment is None:
                continue
            frag_target_indices = [
                donor_idx for donor_idx in fragment.donor_atom_indices
                if donor_idx in donor_target_map
            ]
            if not frag_target_indices:
                continue
            if len(frag_target_indices) >= 2:
                has_multidentate_movable_fragment = True
            if _align_hybrid_fragment_to_targets(
                mol,
                fragment,
                embedded_fragment,
                frag_target_indices,
                donor_target_map,
                reference_center=metal_pos,
                metal_idx=metal_idx,
                exact_target_count=0,
            ):
                moved_fragments += 1
                module_atoms.update(fragment.atom_indices)
                relaxable_atoms.update(fragment.atom_indices)
                keep_rigid_pose = (
                    len(donor_indices) == 4
                    and geom_code in _HAPTO_SECONDARY_CN4_GEOMETRIES
                    and has_secondary_bridge_branch
                    and _secondary_fragment_prefers_rigid_pose(
                        mol,
                        fragment,
                        frag_target_indices,
                    )
                )
                if not keep_rigid_pose:
                    try:
                        _optimize_secondary_fragment_pose(
                            mol,
                            fragment,
                            metal_idx,
                        donor_target_map,
                    )
                    except Exception:
                        pass
                else:
                    rigid_reference_fragments.append(
                        (
                            fragment,
                            {
                                atom_idx: np.array(conf.GetAtomPosition(atom_idx), dtype=float)
                                for atom_idx in fragment.atom_indices
                            },
                        )
                    )
                if keep_rigid_pose and sum(
                    1 for donor_idx in frag_target_indices
                    if mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'O'
                ) >= 2:
                    oo_relief_fragments.append(fragment)

        allow_post_fragment_refit = not (
            len(donor_indices) == 4
            and geom_code in {'SQ', 'TET'}
            and has_multidentate_movable_fragment
            and has_secondary_bridge_branch
        )

        if moved_fragments and allow_post_fragment_refit:
            def _module_error(candidate_metal_pos) -> float:
                err = 0.0
                metal_arr = np.asarray(candidate_metal_pos, dtype=float)
                coord_map: Dict[int, np.ndarray] = {}
                for fragment in decomposition.fragments:
                    if metal_idx not in fragment.metal_neighbor_indices:
                        continue
                    for atom_idx in fragment.atom_indices:
                        if atom_idx in coord_map:
                            continue
                        pos = conf.GetAtomPosition(atom_idx)
                        coord_map[atom_idx] = np.array([pos.x, pos.y, pos.z], dtype=float)
                for donor_idx in donor_indices:
                    donor_pos = np.array(conf.GetAtomPosition(donor_idx), dtype=float)
                    target_len = float(_secondary_donor_target_length(mol, metal_sym, donor_idx))
                    weight = 2.0 * float(_secondary_donor_fit_weight(mol, metal_sym, donor_idx))
                    err += weight * (float(np.linalg.norm(donor_pos - metal_arr)) - target_len) ** 2
                seen_fragments: set = set()
                for donor_idx in donor_indices:
                    fragment = donor_to_fragment.get(donor_idx)
                    if fragment is None:
                        continue
                    frag_key = tuple(fragment.atom_indices)
                    if frag_key in seen_fragments:
                        continue
                    seen_fragments.add(frag_key)
                    fragment_donors = [
                        idx for idx in fragment.donor_atom_indices
                        if idx in set(donor_indices)
                    ]
                    if not fragment_donors:
                        continue
                    err += _secondary_non_donor_contact_penalty(
                        mol,
                        metal_idx,
                        metal_arr,
                        fragment.atom_indices,
                        fragment_donors,
                        coord_map,
                    )
                    err += _planar_fragment_metal_coplanarity_penalty(
                        mol,
                        fragment.atom_indices,
                        fragment_donors,
                        metal_arr,
                        coord_map,
                    )
                    err += _planar_fragment_donor_approach_penalty(
                        mol,
                        fragment.atom_indices,
                        fragment_donors,
                        metal_arr,
                        coord_map,
                    )
                    if len(fragment_donors) >= 2:
                        err += _fragment_bite_direction_penalty(
                            mol,
                            fragment.atom_indices,
                            fragment_donors,
                            coord_map,
                            metal_arr,
                        )
                return err

            current_metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
            best_metal_pos = current_metal_pos
            best_metal_err = _module_error(current_metal_pos)
            donor_positions = []
            for donor_idx in donor_indices:
                pos = conf.GetAtomPosition(donor_idx)
                donor_positions.append(np.array([pos.x, pos.y, pos.z], dtype=float))
            refit_candidates = _enumerate_secondary_metal_geometry_fits(
                donor_positions,
                donor_symbols,
                metal_sym,
                constrained_indices=fit_slot_indices,
                bite_distance_constraints=bite_distance_constraints,
                donor_target_lengths=donor_target_lengths,
                donor_fit_weights=donor_fit_weights,
                max_candidates=max(variant_rank + 1, 4),
            )
            refit = refit_candidates[min(variant_rank, len(refit_candidates) - 1)] if refit_candidates else None
            if refit is not None:
                refit_metal_pos, refit_geom_code, _target_positions, _refit_score = refit
                refit_err = _module_error(refit_metal_pos)
                if refit_err + 1e-9 < best_metal_err:
                    best_metal_pos = np.asarray(refit_metal_pos, dtype=float)
                    best_metal_err = refit_err
                    geom_code = refit_geom_code
            metal_pos = best_metal_pos
        elif moved_fragments and not allow_post_fragment_refit:
            logger.info(
                "Hybrid hapto kept pre-fragment-fit metal position for %s%d: skipped post-fragment refit for bridged CN=4 module",
                metal_sym,
                metal_idx,
            )

        conf.SetAtomPosition(
            metal_idx,
            Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
        )
        optimized_atoms: set = set()
        try:
            optimized_atoms = _optimize_secondary_metal_module_local(
                mol,
                decomposition.fragments,
                hapto_groups,
                metal_idx,
                donor_indices,
                donor_to_fragment,
                geometry_code=geom_code,
                geometry_metal_pos=metal_pos,
                geometry_target_map=donor_target_map,
            )
        except Exception as e:
            logger.debug("Hybrid hapto local secondary-metal optimization skipped: %s", e)
            optimized_atoms = set()
        if optimized_atoms:
            module_atoms.update(optimized_atoms)
            relaxable_atoms.update(optimized_atoms - {metal_idx})
            pos = conf.GetAtomPosition(metal_idx)
            metal_pos = np.array([pos.x, pos.y, pos.z], dtype=float)
        restored_rigid_count = 0
        for fragment, reference_coords in rigid_reference_fragments:
            try:
                if _restore_secondary_rigid_fragment_geometry(
                    mol,
                    fragment,
                    metal_idx,
                    donor_target_map,
                    reference_coords,
                ):
                    restored_rigid_count += 1
            except Exception:
                continue
        oo_relief_count = 0
        for fragment in oo_relief_fragments:
            try:
                if _relieve_secondary_oo_chelate_contacts(
                    mol,
                    fragment,
                    metal_idx,
                    donor_target_map,
                ):
                    oo_relief_count += 1
            except Exception:
                continue
        placed_secondary.add(metal_idx)
        module_atoms.add(metal_idx)
        module_atoms.update(donor_indices)
        logger.info(
            "Hybrid hapto built secondary metal module %s%d via %s fit to %d donor(s); moved %d fragment(s), reoriented %d scaffold branch(es), locally optimized %d atom(s), restored %d rigid fragment(s), relieved %d O,O chelate(s)",
            metal_sym,
            metal_idx,
            geom_code,
            len(donor_indices),
            moved_fragments,
            reoriented_scaffold_branches,
            len(optimized_atoms),
            restored_rigid_count,
            oo_relief_count,
        )

    return placed_secondary, module_atoms, relaxable_atoms


def _build_primary_organometal_hapto_module(
    mol,
    decomposition: Optional[_HybridHaptoDecomposition],
    hapto_groups: List[Tuple[int, List[int]]],
    module: Optional[_PrimaryOrganometalModule],
    preview_only: bool = False,
    preview_store: Optional[List[Tuple[str, str]]] = None,
):
    """Build a single-metal hapto organometal module before falling back to the legacy path."""
    if (
        not RDKIT_AVAILABLE
        or mol is None
        or decomposition is None
        or module is None
        or not hapto_groups
    ):
        return None
    try:
        import numpy as np
    except ImportError:
        return None

    scaffold_mol = Chem.Mol(mol)
    if not _build_hapto_scaffold(scaffold_mol, hapto_groups):
        return None

    try:
        _correct_hapto_geometry(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    try:
        conf = scaffold_mol.GetConformer(0)
    except Exception:
        return None

    metal_idx = int(module.metal_idx)
    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
    donor_target_map: Dict[int, np.ndarray] = {}
    for donor_idx in module.donor_atom_indices:
        pos = conf.GetAtomPosition(donor_idx)
        donor_target_map[donor_idx] = np.array([pos.x, pos.y, pos.z], dtype=float)

    trusted_atoms: set = {
        atom_idx
        for group_idx in module.hapto_group_ids
        if 0 <= int(group_idx) < len(hapto_groups)
        for atom_idx in hapto_groups[int(group_idx)][1]
    }
    aligned_fragments = 0
    aligned_correlated = 0

    target_fragment_indices = set(module.correlated_fragment_indices) | set(module.terminal_fragment_indices)
    for frag_idx, fragment in enumerate(decomposition.fragments):
        if frag_idx not in target_fragment_indices:
            if any(group_id in module.hapto_group_ids for group_id in fragment.hapto_group_ids):
                trusted_atoms.update(fragment.atom_indices)
            continue

        embedded_fragment = _embed_hybrid_fragment(fragment.fragment_mol)
        if embedded_fragment is None:
            if frag_idx in module.correlated_fragment_indices:
                return None
            continue

        frag_target_indices = [
            donor_idx for donor_idx in fragment.donor_atom_indices
            if donor_idx in donor_target_map
        ]
        if not frag_target_indices:
            if frag_idx in module.correlated_fragment_indices:
                return None
            continue

        if not _align_hybrid_fragment_to_targets(
            scaffold_mol,
            fragment,
            embedded_fragment,
            frag_target_indices,
            donor_target_map,
            reference_center=metal_pos,
            metal_idx=metal_idx,
            exact_target_count=0,
        ):
            if frag_idx in module.correlated_fragment_indices:
                return None
            continue

        aligned_fragments += 1
        if frag_idx in module.correlated_fragment_indices:
            aligned_correlated += 1
        trusted_atoms.update(fragment.atom_indices)

        if (
            len(frag_target_indices) >= 2
            and not _secondary_fragment_prefers_rigid_pose(
                scaffold_mol,
                fragment,
                frag_target_indices,
            )
        ):
            try:
                _optimize_secondary_fragment_pose(
                    scaffold_mol,
                    fragment,
                    metal_idx,
                    donor_target_map,
                )
            except Exception:
                pass

    if aligned_correlated < len(module.correlated_fragment_indices):
        return None
    if aligned_fragments == 0:
        return None

    try:
        _propagate_non_hapto_atoms(
            scaffold_mol,
            0,
            hapto_groups,
            extra_fixed_indices=trusted_atoms | _hapto_primary_donor_indices(scaffold_mol, hapto_groups),
        )
    except Exception:
        pass

    if preview_store is not None:
        try:
            preview_xyz = _mol_to_xyz(scaffold_mol)
            preview_entry = (preview_xyz, 'quick-primary-preview')
            if preview_entry not in preview_store:
                preview_store.append(preview_entry)
        except Exception:
            pass

    if (
        not preview_only
        and not _primary_organometal_module_quality_ok(
            scaffold_mol,
            decomposition,
            module,
        )
    ):
        return None

    logger.info(
        "Hybrid hapto primary organometal module builder succeeded for %s%d with %d correlated fragment(s)%s",
        mol.GetAtomWithIdx(metal_idx).GetSymbol(),
        metal_idx,
        len(module.correlated_fragment_indices),
        " [preview]" if preview_only else "",
    )
    return scaffold_mol


def _place_secondary_metals_in_hapto_fragments(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> set:
    """Place non-hapto metals from their donor clouds in hapto-coupled systems."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return set()
    try:
        import numpy as np
    except ImportError:
        return set()

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return set()

    hapto_metals = {metal_idx for metal_idx, _grp in hapto_groups}
    placed_secondary: set = set()

    for atom in mol.GetAtoms():
        metal_idx = atom.GetIdx()
        metal_sym = atom.GetSymbol()
        if metal_sym not in _METAL_SET or metal_idx in hapto_metals:
            continue

        donor_indices = [
            nbr.GetIdx()
            for nbr in atom.GetNeighbors()
            if nbr.GetSymbol() not in _METAL_SET and nbr.GetAtomicNum() > 1
        ]
        if len(donor_indices) < 2:
            continue

        donor_positions = []
        donor_symbols = []
        donor_target_lengths = []
        donor_fit_weights = []
        skip = False
        for donor_idx in donor_indices:
            pos = conf.GetAtomPosition(donor_idx)
            arr = np.array([pos.x, pos.y, pos.z], dtype=float)
            if not np.all(np.isfinite(arr)):
                skip = True
                break
            donor_positions.append(arr)
            donor_symbols.append(mol.GetAtomWithIdx(donor_idx).GetSymbol())
            donor_target_lengths.append(_secondary_donor_target_length(mol, metal_sym, donor_idx))
            donor_fit_weights.append(_secondary_donor_fit_weight(mol, metal_sym, donor_idx))
        if skip:
            continue

        fit = _fit_secondary_metal_position_from_donors(
            donor_positions,
            donor_symbols,
            metal_sym,
            donor_target_lengths=donor_target_lengths,
            donor_fit_weights=donor_fit_weights,
        )
        if fit is None:
            continue
        metal_pos, geom_code = fit
        conf.SetAtomPosition(
            metal_idx,
            Point3D(float(metal_pos[0]), float(metal_pos[1]), float(metal_pos[2])),
        )
        placed_secondary.add(metal_idx)
        logger.info(
            "Hybrid hapto placed secondary metal %s%d via %s fit to %d donor(s)",
            metal_sym,
            metal_idx,
            geom_code,
            len(donor_indices),
        )

    return placed_secondary


def _hybrid_bond_target_length(bond) -> float:
    """Approximate covalent bond target length for hybrid hapto relaxation."""
    begin_atom = bond.GetBeginAtom()
    end_atom = bond.GetEndAtom()
    begin_sym = begin_atom.GetSymbol()
    end_sym = end_atom.GetSymbol()

    if begin_sym in _METAL_SET or end_sym in _METAL_SET:
        metal_sym = begin_sym if begin_sym in _METAL_SET else end_sym
        donor_sym = end_sym if begin_sym in _METAL_SET else begin_sym
        return float(_get_ml_bond_length(metal_sym, donor_sym))

    pair = frozenset([begin_sym, end_sym])
    bond_type = bond.GetBondType()

    if 'H' in pair:
        other = end_sym if begin_sym == 'H' else begin_sym
        return {'C': 1.08, 'N': 1.01, 'O': 0.96, 'Si': 1.48, 'B': 1.19}.get(other, 1.08)

    if bond.GetIsAromatic() or bond_type == Chem.BondType.AROMATIC:
        aromatic_map = {
            frozenset(['C', 'C']): 1.39,
            frozenset(['C', 'N']): 1.35,
            frozenset(['C', 'O']): 1.36,
        }
        if pair in aromatic_map:
            return aromatic_map[pair]
        return max(_COVALENT_RADII.get(begin_sym, 0.76) + _COVALENT_RADII.get(end_sym, 0.76) - 0.12, 0.95)

    if bond_type == Chem.BondType.TRIPLE:
        triple_map = {
            frozenset(['C', 'C']): 1.20,
            frozenset(['C', 'N']): 1.16,
            frozenset(['N', 'N']): 1.10,
        }
        if pair in triple_map:
            return triple_map[pair]
        return max(_COVALENT_RADII.get(begin_sym, 0.76) + _COVALENT_RADII.get(end_sym, 0.76) - 0.22, 0.90)

    if bond_type == Chem.BondType.DOUBLE:
        double_map = {
            frozenset(['C', 'C']): 1.34,
            frozenset(['C', 'N']): 1.29,
            frozenset(['C', 'O']): 1.23,
            frozenset(['N', 'O']): 1.22,
            frozenset(['N', 'N']): 1.25,
        }
        if pair in double_map:
            return double_map[pair]
        return max(_COVALENT_RADII.get(begin_sym, 0.76) + _COVALENT_RADII.get(end_sym, 0.76) - 0.12, 0.95)

    single_map = {
        frozenset(['Si', 'C']): 1.87,
        frozenset(['C', 'C']): 1.50,
        frozenset(['C', 'N']): 1.47,
        frozenset(['C', 'O']): 1.43,
        frozenset(['C', 'Si']): 1.87,
        frozenset(['C', 'S']): 1.82,
        frozenset(['C', 'P']): 1.84,
        frozenset(['N', 'N']): 1.45,
        frozenset(['N', 'O']): 1.40,
        frozenset(['Si', 'Si']): 2.34,
    }
    if pair in single_map:
        return single_map[pair]
    return _COVALENT_RADII.get(begin_sym, 0.76) + _COVALENT_RADII.get(end_sym, 0.76)


def _refine_hybrid_hapto_complex(
    mol,
    decomposition: _HybridHaptoDecomposition,
    hapto_groups: List[Tuple[int, List[int]]],
    trusted_atom_indices: Optional[set] = None,
    secondary_metal_indices: Optional[set] = None,
    apply_final_hapto_correction: bool = True,
) -> bool:
    """Local hapto-only relaxation around a fixed coordination scaffold."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        import numpy as np
    except ImportError:
        return False
    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    metals = set(decomposition.metal_indices)
    if secondary_metal_indices:
        metals = metals | set(secondary_metal_indices)
    hapto_atoms = {idx for _metal_idx, grp in hapto_groups for idx in grp}
    primary_hapto_donors = _hapto_primary_donor_indices(mol, hapto_groups)
    fixed = metals | hapto_atoms | primary_hapto_donors
    if trusted_atom_indices:
        fixed = fixed | set(trusted_atom_indices)
    movable = {idx for idx in range(mol.GetNumAtoms()) if idx not in fixed}
    if not movable:
        return True

    def _gp(idx):
        pos = conf.GetAtomPosition(idx)
        return np.array([pos.x, pos.y, pos.z], dtype=float)

    def _sp(idx, arr):
        conf.SetAtomPosition(idx, Point3D(float(arr[0]), float(arr[1]), float(arr[2])))

    bond_targets = []
    bonded_pairs: set = set()
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        bonded_pairs.add((begin_idx, end_idx))
        bonded_pairs.add((end_idx, begin_idx))
        if begin_idx not in movable and end_idx not in movable:
            continue
        bond_targets.append((begin_idx, end_idx, _hybrid_bond_target_length(bond)))

    rng = np.random.default_rng(42)
    for _pass in range(80):
        displacements: Dict[int, np.ndarray] = {idx: np.zeros(3, dtype=float) for idx in movable}
        max_err = 0.0
        for begin_idx, end_idx, target in bond_targets:
            begin_pos = _gp(begin_idx)
            end_pos = _gp(end_idx)
            diff = end_pos - begin_pos
            dist = float(np.linalg.norm(diff))
            if dist < 1e-8:
                diff = rng.standard_normal(3)
                dist = float(np.linalg.norm(diff))
            unit = diff / max(dist, 1e-12)
            delta = dist - target
            if abs(delta) < 0.01:
                continue
            max_err = max(max_err, abs(delta))
            strength = 0.35 if ((begin_idx in fixed) ^ (end_idx in fixed)) else 0.22
            move = strength * delta * unit
            if begin_idx in movable and end_idx in movable:
                displacements[begin_idx] = displacements[begin_idx] + 0.5 * move
                displacements[end_idx] = displacements[end_idx] - 0.5 * move
            elif begin_idx in movable:
                displacements[begin_idx] = displacements[begin_idx] + move
            elif end_idx in movable:
                displacements[end_idx] = displacements[end_idx] - move

        if max_err < 0.04:
            break
        for atom_idx, disp in displacements.items():
            if float(np.linalg.norm(disp)) < 1e-6:
                continue
            _sp(atom_idx, _gp(atom_idx) + disp)

    for _pass in range(24):
        any_push = False
        for i in range(mol.GetNumAtoms()):
            for j in range(i + 1, mol.GetNumAtoms()):
                if (i, j) in bonded_pairs:
                    continue
                if i not in movable and j not in movable:
                    continue
                pi = _gp(i)
                pj = _gp(j)
                dist = float(np.linalg.norm(pi - pj))
                sym_i = mol.GetAtomWithIdx(i).GetSymbol()
                sym_j = mol.GetAtomWithIdx(j).GetSymbol()

                if sym_i in _METAL_SET or sym_j in _METAL_SET:
                    min_dist = 2.0
                elif sym_i == 'H' and sym_j == 'H':
                    min_dist = 1.5
                elif sym_i == 'H' or sym_j == 'H':
                    min_dist = 1.0
                else:
                    min_dist = 1.18

                if dist >= min_dist:
                    continue
                if dist < 1e-8:
                    direction = rng.standard_normal(3)
                    direction /= max(float(np.linalg.norm(direction)), 1e-12)
                else:
                    direction = (pj - pi) / dist

                gap = min_dist - dist
                if i in movable and j in movable:
                    _sp(i, pi - 0.5 * gap * direction)
                    _sp(j, pj + 0.5 * gap * direction)
                elif i in movable:
                    _sp(i, pi - gap * direction)
                elif j in movable:
                    _sp(j, pj + gap * direction)
                any_push = True
        if not any_push:
            break

    try:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
        if ff is not None:
            for atom_idx in sorted(fixed):
                ff.AddFixedPoint(atom_idx)
            ff.Minimize(maxIts=250)
    except Exception:
        pass

    if apply_final_hapto_correction:
        try:
            _correct_hapto_geometry(mol, 0, hapto_groups)
        except Exception:
            pass
    return True


def _enforce_donor_pi_coplanarity(
    mol,
    conf_id: int = 0,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
    max_rotation_deg: float = 75.0,
) -> bool:
    """Rigidly rotate planar donor fragments so the metal lies in the pi-plane.

    For each metal, identifies planar (aromatic / sp2-conjugated) fragments
    that bear donor atoms, then rotates each fragment as a rigid body so that
    the metal atom ends up in (or very near) the fragment's best-fit plane.

    * Bidentate+ donors  -> rotation around the donor-donor axis.
    * Monodentate donor   -> rotation around the single donor atom.

    Donor-metal bond lengths are preserved exactly because the rotation axis
    passes through the donor atom(s).
    """
    if not RDKIT_AVAILABLE or mol is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return False

    hapto_atom_set: set = set()
    hapto_metal_set: set = set()
    if hapto_groups:
        for _mi, grp in hapto_groups:
            hapto_atom_set.update(grp)
            hapto_metal_set.add(_mi)

    max_rot = np.radians(max_rotation_deg)

    def _gp(idx):
        p = conf.GetAtomPosition(idx)
        return np.array([p.x, p.y, p.z], dtype=float)

    def _sp(idx, arr):
        conf.SetAtomPosition(
            idx, Point3D(float(arr[0]), float(arr[1]), float(arr[2])))

    def _rod(v, k, angle):
        c, s = np.cos(angle), np.sin(angle)
        return v * c + np.cross(k, v) * s + k * float(np.dot(k, v)) * (1.0 - c)

    def _is_planar_atom(a):
        if a.GetIsAromatic():
            return True
        try:
            hyb = a.GetHybridization()
        except Exception:
            return False
        if hyb in {Chem.rdchem.HybridizationType.SP,
                    Chem.rdchem.HybridizationType.SP2}:
            return True
        for bond in a.GetBonds():
            if bond.GetBondTypeAsDouble() >= 1.5 or bond.GetIsConjugated():
                return True
        return False

    changed = False

    for atom in mol.GetAtoms():
        metal_idx = atom.GetIdx()
        if atom.GetSymbol() not in _METAL_SET:
            continue
        # Skip hapto metals: piano-stool geometry should not be coplanarised
        if metal_idx in hapto_metal_set:
            continue

        # Non-hapto, non-metal, heavy-atom neighbours = candidate donors
        donors = []
        for nbr in atom.GetNeighbors():
            ni = nbr.GetIdx()
            if ni in hapto_atom_set or nbr.GetSymbol() in _METAL_SET or nbr.GetAtomicNum() <= 1:
                continue
            donors.append(ni)
        if not donors:
            continue

        m_pos = _gp(metal_idx)
        processed: set = set()

        for start_donor in donors:
            if start_donor in processed:
                continue

            if not _is_planar_atom(mol.GetAtomWithIdx(start_donor)):
                processed.add(start_donor)
                continue

            # --- BFS: find connected planar core ---
            planar_core: set = {start_donor}
            bfs_q = [start_donor]
            while bfs_q:
                cur = bfs_q.pop(0)
                for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in planar_core or ni in hapto_atom_set:
                        continue
                    if nbr.GetSymbol() in _METAL_SET or nbr.GetAtomicNum() <= 1:
                        continue
                    if _is_planar_atom(nbr):
                        planar_core.add(ni)
                        bfs_q.append(ni)

            if len(planar_core) < 3:
                processed.add(start_donor)
                continue

            frag_donors = [d for d in donors if d in planar_core]
            processed.update(frag_donors)
            frag_donor_set = set(frag_donors)

            # --- Best-fit plane (SVD) ---
            core_list = sorted(planar_core)
            pts = np.array([_gp(i) for i in core_list])
            centroid = pts.mean(axis=0)
            q = pts - centroid
            try:
                _u, _s, vh = np.linalg.svd(q, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            normal = vh[-1].copy()
            n_len = float(np.linalg.norm(normal))
            if n_len < 1e-12:
                continue
            normal /= n_len

            # --- Rigid body: planar core + direct substituents + their H ---
            rigid_body: set = set(planar_core)
            for ci in list(planar_core):
                for nbr in mol.GetAtomWithIdx(ci).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in rigid_body or ni in hapto_atom_set or nbr.GetSymbol() in _METAL_SET:
                        continue
                    rigid_body.add(ni)
                    if nbr.GetAtomicNum() > 1:
                        for nbr2 in nbr.GetNeighbors():
                            ni2 = nbr2.GetIdx()
                            if ni2 not in rigid_body and nbr2.GetAtomicNum() == 1:
                                rigid_body.add(ni2)

            # ==========================================================
            # Bidentate+: rotate around the donor-donor axis
            # ==========================================================
            if len(frag_donors) >= 2:
                d1, d2 = frag_donors[0], frag_donors[1]
                d1_pos, d2_pos = _gp(d1), _gp(d2)
                axis = d2_pos - d1_pos
                axis_len = float(np.linalg.norm(axis))
                if axis_len < 1e-8:
                    continue
                u_ax = axis / axis_len

                # Gram-Schmidt: make normal exactly perpendicular to axis
                normal = normal - float(np.dot(normal, u_ax)) * u_ax
                n_len = float(np.linalg.norm(normal))
                if n_len < 1e-8:
                    continue
                normal /= n_len

                # a*cos(t) + b*sin(t) = 0  =>  t = atan2(-a, b)
                v_M = m_pos - d1_pos
                a_coeff = float(np.dot(v_M, normal))
                b_dir = np.cross(u_ax, normal)
                b_coeff = float(np.dot(v_M, b_dir))

                if abs(a_coeff) < 0.05:
                    continue

                theta = float(np.arctan2(-a_coeff, b_coeff))

                # Outward-side check: metal should face donors, not ring body
                body_atoms = [i for i in planar_core if i not in frag_donor_set]
                if body_atoms:
                    body_ctr = np.mean([_gp(i) for i in body_atoms], axis=0)
                    body_rot = d1_pos + _rod(body_ctr - d1_pos, u_ax, theta)
                    d_ctr = (d1_pos + d2_pos) / 2.0
                    outward = d_ctr - body_rot
                    if float(np.dot(m_pos - d_ctr, outward)) < 0:
                        theta += np.pi

                if abs(theta) > max_rot:
                    theta = float(np.sign(theta)) * max_rot

                for ai in rigid_body:
                    p = _gp(ai) - d1_pos
                    _sp(ai, d1_pos + _rod(p, u_ax, theta))
                changed = True

            # ==========================================================
            # Monodentate: rotate around the single donor
            # ==========================================================
            else:
                d_idx = frag_donors[0]
                d_pos = _gp(d_idx)

                v = m_pos - d_pos
                v_len = float(np.linalg.norm(v))
                if v_len < 1e-8:
                    continue

                h = float(np.dot(v, normal))
                if abs(h) < 0.05:
                    continue

                # Desired normal: component of normal perpendicular to v
                n_proj = normal - (float(np.dot(normal, v)) / (v_len ** 2)) * v
                n_proj_len = float(np.linalg.norm(n_proj))
                if n_proj_len < 1e-8:
                    continue
                desired = n_proj / n_proj_len

                rot_axis = np.cross(normal, desired)
                ra_len = float(np.linalg.norm(rot_axis))
                if ra_len < 1e-10:
                    continue
                rot_axis /= ra_len
                cos_a = float(np.clip(np.dot(normal, desired), -1.0, 1.0))
                theta = float(np.arccos(cos_a))

                if theta > max_rot:
                    theta = max_rot

                for ai in rigid_body:
                    p = _gp(ai) - d_pos
                    _sp(ai, d_pos + _rod(p, rot_axis, theta))
                changed = True

    return changed


def _build_multimetal_hapto_sequential(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
    decomposition: '_HybridHaptoDecomposition',
    secondary_variant_plan: Optional[Dict[int, int]] = None,
):
    """Sequential multi-metal builder with rigid-body clash-free placement.

    Strategy:
      1. Build hapto scaffold (ferrocene / Cp rings + hapto metal)
      2. Embed & align bridge fragments (Procrustes + UFF relax)
      3. Place secondary metal(s) from donor positions
      4. Place remaining fragments as rigid bodies, rotating around
         metal-donor axis to avoid clashes with already-placed atoms
      5. Whole-fragment M-L correction (rigid shift, no bond breaking)
      6. Final topology-preserving clash resolution
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups or decomposition is None:
        return None
    try:
        import numpy as np
        from rdkit.Geometry import Point3D
    except ImportError:
        return None

    # ---- helpers ----
    def _gp(conf, idx):
        p = conf.GetAtomPosition(idx)
        return np.array([p.x, p.y, p.z], dtype=float)

    def _sp(conf, idx, pos):
        conf.SetAtomPosition(idx, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))

    def _check_clash(conf, new_positions, placed_set, bonded_pairs, threshold=1.5):
        """Check if new_positions clash with placed atoms.
        Returns worst clash distance (0 = no clash, smaller = worse)."""
        worst = 999.0
        for new_idx, new_pos in new_positions.items():
            for placed_idx in placed_set:
                if placed_idx == new_idx:
                    continue
                pair = (min(new_idx, placed_idx), max(new_idx, placed_idx))
                if pair in bonded_pairs:
                    continue
                placed_pos = _gp(conf, placed_idx)
                d = float(np.linalg.norm(new_pos - placed_pos))
                if d < worst:
                    worst = d
        return worst

    def _rodrigues_rotate(points, axis, angle):
        """Rotate points around axis by angle (radians) using Rodrigues."""
        axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        return (points * cos_a
                + np.cross(axis, points) * sin_a
                + axis[np.newaxis, :] * (points @ axis[:, np.newaxis]) * (1 - cos_a))

    # ---- classify atoms ----
    hapto_metal_set = {mi for mi, _g in hapto_groups}
    all_hapto_carbons: set = set()
    hapto_carbons_by_metal: Dict[int, set] = {}
    for mi, grp in hapto_groups:
        all_hapto_carbons.update(grp)
        hapto_carbons_by_metal.setdefault(mi, set()).update(grp)

    all_metals: set = set()
    non_hapto_metals: set = set()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METAL_SET:
            all_metals.add(atom.GetIdx())
            if atom.GetIdx() not in hapto_metal_set:
                non_hapto_metals.add(atom.GetIdx())

    if not non_hapto_metals and len(hapto_metal_set) < 2:
        return None  # not a multi-metal system

    # ==== STEP 1: Build hapto scaffold (Fe + Cp rings) ====
    scaffold_mol = Chem.Mol(mol)
    if not _build_hapto_scaffold(scaffold_mol, hapto_groups):
        return None
    try:
        _correct_hapto_geometry(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    try:
        conf = scaffold_mol.GetConformer(0)
    except Exception:
        return None

    # Only trust hapto metals initially — secondary metals are NOT placed yet
    trusted: set = set(hapto_metal_set) | set(all_hapto_carbons)
    logger.info(
        "Sequential multi-metal: scaffold built, %d hapto atoms + %d hapto metals trusted",
        len(all_hapto_carbons), len(hapto_metal_set),
    )

    # Helper: Rodrigues rotation (used in Steps 2 and 4)
    def _rotate_fragment_around_axis(orig_positions, pivot, axis, angle):
        """Rotate a dict {idx: pos} around pivot+axis by angle radians."""
        ax = axis / max(float(np.linalg.norm(axis)), 1e-12)
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        result = {}
        for idx, pos in orig_positions.items():
            v = pos - pivot
            v_rot = (v * cos_a
                     + np.cross(ax, v) * sin_a
                     + ax * float(np.dot(ax, v)) * (1 - cos_a))
            result[idx] = pivot + v_rot
        return result

    # ==== STEP 2: Embed & align all ligand fragments ====
    # Strategy: ETKDG embed (gives correct ring geometry) → Procrustes align
    # to scaffold anchors (with anchor snap) → UFF relaxation with anchors
    # fixed to repair any bond stretching caused by the snap.
    for fragment in decomposition.fragments:
        # Pure hapto-ring fragments: trust scaffold geometry already
        if all(ai in all_hapto_carbons for ai in fragment.atom_indices):
            trusted.update(fragment.atom_indices)
            continue

        # Check if fragment has anchors in trusted set
        has_anchor = any(ai in trusted for ai in fragment.atom_indices)
        if not has_anchor:
            continue  # defer to Step 4

        # Embed fragment with ETKDG (correct ring geometry)
        embedded = _embed_hybrid_fragment(fragment.fragment_mol)
        if embedded is None:
            continue

        # Procrustes align to scaffold (with anchor snap for junction correctness)
        if not _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded):
            continue

        # UFF relaxation: fix anchor atoms, let non-anchors relax to fix
        # any bond stretching caused by the anchor snap.
        try:
            # Build a temporary molecule from the fragment with aligned coords
            frag_mol = Chem.Mol(fragment.fragment_mol)
            frag_mol.RemoveAllConformers()
            frag_conf_obj = Chem.Conformer(frag_mol.GetNumAtoms())
            for frag_idx, orig_idx in fragment.fragment_to_original.items():
                p = conf.GetAtomPosition(orig_idx)
                frag_conf_obj.SetAtomPosition(frag_idx, Point3D(p.x, p.y, p.z))
            frag_mol.AddConformer(frag_conf_obj, assignId=True)

            ff = AllChem.UFFGetMoleculeForceField(frag_mol)
            if ff is not None:
                # Fix anchor atoms (Cp carbons that must stay at scaffold positions)
                n_fixed = 0
                for orig_idx in fragment.anchor_atom_indices:
                    frag_idx = fragment.original_to_fragment.get(orig_idx)
                    if frag_idx is not None and orig_idx in trusted:
                        ff.AddFixedPoint(frag_idx)
                        n_fixed += 1
                e_before = ff.CalcEnergy()
                ff.Minimize(maxIts=500)
                e_after = ff.CalcEnergy()
                logger.debug(
                    "UFF relax: %d fixed, E %.1f -> %.1f",
                    n_fixed, e_before, e_after,
                )

                # Copy ALL relaxed coordinates back to scaffold
                # (anchors were fixed so their positions didn't change)
                relaxed_conf = frag_mol.GetConformer()
                for frag_idx, orig_idx in fragment.fragment_to_original.items():
                    p = relaxed_conf.GetAtomPosition(frag_idx)
                    conf.SetAtomPosition(orig_idx, Point3D(p.x, p.y, p.z))
        except Exception as e:
            logger.debug("UFF relaxation failed for fragment: %s", e)

        # ---- Clash-aware torsion for bridge fragments ----
        # If this fragment contains hapto carbons, the non-hapto substituents
        # (e.g., pyridine ring) may clash with the OTHER Cp ring of the
        # metallocene.  Fix by rotating the non-hapto substructure around
        # the bond connecting it to the hapto ring.
        frag_hapto = [ai for ai in fragment.atom_indices if ai in all_hapto_carbons]
        frag_non_hapto_heavy = [ai for ai in fragment.atom_indices
                                if ai not in all_hapto_carbons
                                and ai not in all_metals
                                and mol.GetAtomWithIdx(ai).GetAtomicNum() > 1]
        if frag_hapto and frag_non_hapto_heavy:
            # Find the OTHER Cp ring's atoms (for clash checking)
            other_cp_atoms: set = set()
            for hm in hapto_metal_set:
                for nbr in mol.GetAtomWithIdx(hm).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in all_hapto_carbons and ni not in frag_hapto:
                        other_cp_atoms.add(ni)

            # Find the bond connecting hapto → non-hapto HEAVY atom
            bridge_bond = None
            for hi in frag_hapto:
                for nbr in mol.GetAtomWithIdx(hi).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in frag_non_hapto_heavy:
                        bridge_bond = (hi, ni)
                        break
                if bridge_bond is not None:
                    break

            if bridge_bond is not None and other_cp_atoms:
                pivot_idx, next_idx = bridge_bond
                pivot_pos = _gp(conf, pivot_idx)
                next_pos = _gp(conf, next_idx)
                rot_axis = next_pos - pivot_pos
                rot_axis_len = float(np.linalg.norm(rot_axis))

                if rot_axis_len > 1e-6:
                    # Check all clash-relevant atoms (other Cp + all trusted
                    # not in this fragment)
                    clash_targets = (other_cp_atoms | trusted) - set(fragment.atom_indices)
                    # Collect atoms to rotate: all non-hapto atoms (heavy + H)
                    frag_non_hapto_all = [ai for ai in fragment.atom_indices
                                          if ai not in all_hapto_carbons
                                          and ai not in all_metals]
                    rot_atoms = list(frag_non_hapto_all)
                    orig_positions = {ai: _gp(conf, ai) for ai in rot_atoms}

                    def _score_rotation(angle):
                        if abs(angle) < 1e-9:
                            positions = orig_positions
                        else:
                            positions = _rotate_fragment_around_axis(
                                orig_positions, pivot_pos, rot_axis, angle)
                        min_d = 999.0
                        for ai, pos in positions.items():
                            for ti in clash_targets:
                                d = float(np.linalg.norm(pos - _gp(conf, ti)))
                                if d < min_d:
                                    min_d = d
                        return min_d

                    best_angle = 0.0
                    best_min_d = _score_rotation(0.0)
                    for step in range(1, 36):
                        angle = step * (2 * np.pi / 36)
                        min_d = _score_rotation(angle)
                        if min_d > best_min_d:
                            best_min_d = min_d
                            best_angle = angle

                    if best_angle != 0.0:
                        rotated = _rotate_fragment_around_axis(
                            orig_positions, pivot_pos, rot_axis, best_angle)
                        for ai, pos in rotated.items():
                            _sp(conf, ai, pos)
                        logger.info(
                            "Step 2 torsion: rotated %d non-hapto atoms by %.0f° "
                            "around %s%d-%s%d, min_d=%.2f",
                            len(rot_atoms), np.degrees(best_angle),
                            mol.GetAtomWithIdx(pivot_idx).GetSymbol(), pivot_idx,
                            mol.GetAtomWithIdx(next_idx).GetSymbol(), next_idx,
                            best_min_d,
                        )

        trusted.update(fragment.atom_indices)
        logger.info(
            "Sequential: embed+align+relax fragment (%d atoms, bridge=%s)",
            len(fragment.atom_indices), fragment.use_scaffold_only,
        )

    # Mark all atoms reachable from trusted set (excluding secondary metals)
    from collections import deque
    bfs_visited: set = set(trusted)
    bfs_q: deque = deque(trusted)
    while bfs_q:
        curr = bfs_q.popleft()
        for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in bfs_visited:
                continue
            if ni in non_hapto_metals:
                continue  # don't cross into secondary metals yet
            bfs_visited.add(ni)
            bfs_q.append(ni)
    trusted = bfs_visited
    logger.info("Sequential Step 2: %d atoms trusted after fragment alignment", len(trusted))

    # ==== STEP 3: Place secondary (non-hapto) metals ====
    for metal_idx in sorted(non_hapto_metals):
        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        variant_rank = max(0, int((secondary_variant_plan or {}).get(metal_idx, 0)))
        # Collect donors that are already placed (in trusted set)
        donor_indices = []
        for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in all_metals:
                continue
            if nbr.GetAtomicNum() <= 1:
                continue
            donor_indices.append(ni)

        placed_donors = [d for d in donor_indices if d in trusted]
        unplaced_donors = [d for d in donor_indices if d not in trusted]

        if len(placed_donors) < 1:
            logger.debug(
                "Sequential: skipping %s%d (no placed donors yet)", metal_sym, metal_idx,
            )
            continue

        # Get positions and symbols of placed donors
        donor_positions = []
        donor_symbols = []
        target_lengths = []
        for d in placed_donors:
            donor_positions.append(_gp(conf, d))
            d_sym = mol.GetAtomWithIdx(d).GetSymbol()
            donor_symbols.append(d_sym)
            target_lengths.append(_get_ml_bond_length(metal_sym, d_sym))

        if len(placed_donors) >= 2:
            # For bimetallic complexes, the LIN/centroid fit from 2 donors
            # often places the metal BETWEEN close donors (inside the
            # ferrocene cage). Strategy:
            #   1. Try SVD fit first
            #   2. Check inter-metal distance; if too close to hapto metals,
            #      use outward placement from hapto center through donor midpoint
            hapto_center = np.mean([_gp(conf, mi) for mi in hapto_metal_set], axis=0)
            avg_len = np.mean(target_lengths)

            fit_candidates = _enumerate_secondary_metal_position_fits(
                donor_positions, donor_symbols, metal_sym,
                donor_target_lengths=target_lengths,
                max_candidates=max(variant_rank + 1, 4),
            )
            fit_result = (
                fit_candidates[min(variant_rank, len(fit_candidates) - 1)]
                if fit_candidates else None
            )
            geom_code = "SVD"
            if fit_result is not None:
                metal_pos, geom_code, _fit_score = fit_result
            else:
                metal_pos = np.mean(donor_positions, axis=0)

            # Check if fitted position is too close to any hapto metal
            too_close = False
            for hm in hapto_metal_set:
                hm_pos = _gp(conf, hm)
                if float(np.linalg.norm(metal_pos - hm_pos)) < 2.5:
                    too_close = True
                    break

            if too_close:
                # Intersection-circle search: for 2+ donors, the ideal
                # metal position lies on the intersection circle of the
                # donor spheres. Search on this circle first, then fall
                # back to multi-sphere search.
                best_score = -1e9
                best_cand = metal_pos
                geom_code = "SEARCH"

                # Build candidate list from intersection circle + sphere
                candidates: list = []

                if len(donor_positions) >= 2:
                    # Compute intersection circle of first two donors
                    c1 = np.array(donor_positions[0])
                    c2 = np.array(donor_positions[1])
                    r1 = target_lengths[0]
                    r2 = target_lengths[1]
                    axis_vec = c2 - c1
                    d_donors = float(np.linalg.norm(axis_vec))
                    if d_donors > 1e-6 and d_donors < r1 + r2:
                        axis_unit = axis_vec / d_donors
                        # Distance from c1 to intersection plane
                        h = (r1 ** 2 - r2 ** 2 + d_donors ** 2) / (2 * d_donors)
                        circle_r_sq = r1 ** 2 - h ** 2
                        if circle_r_sq > 0:
                            circle_r = np.sqrt(circle_r_sq)
                            circle_center = c1 + h * axis_unit
                            # Build orthonormal basis in the plane
                            arb = np.array([1.0, 0.0, 0.0])
                            if abs(float(np.dot(axis_unit, arb))) > 0.9:
                                arb = np.array([0.0, 1.0, 0.0])
                            u = np.cross(axis_unit, arb)
                            u = u / float(np.linalg.norm(u))
                            v = np.cross(axis_unit, u)
                            # Sample circle at 5° resolution
                            for deg in range(0, 360, 5):
                                angle = np.radians(deg)
                                cand = (circle_center
                                        + circle_r * (np.cos(angle) * u + np.sin(angle) * v))
                                candidates.append(cand)

                # Also sample spheres around each donor (10° resolution)
                for center_pos, center_len in zip(donor_positions, target_lengths):
                    center = np.array(center_pos)
                    for theta_deg in range(0, 180, 15):
                        for phi_deg in range(0, 360, 15):
                            theta = np.radians(theta_deg)
                            phi = np.radians(phi_deg)
                            d_vec = np.array([
                                np.sin(theta) * np.cos(phi),
                                np.sin(theta) * np.sin(phi),
                                np.cos(theta),
                            ])
                            candidates.append(center + d_vec * center_len)

                scored_candidates: List[Tuple[float, float, object]] = []
                if fit_result is not None:
                    min_clash_fit = min(
                        (float(np.linalg.norm(np.asarray(metal_pos, dtype=float) - _gp(conf, ti)))
                         for ti in trusted if ti not in all_metals),
                        default=999.0,
                    )
                    min_hm_fit = min(
                        float(np.linalg.norm(np.asarray(metal_pos, dtype=float) - _gp(conf, hm)))
                        for hm in hapto_metal_set
                    )
                    ml_err_fit = sum(
                        (float(np.linalg.norm(np.asarray(metal_pos, dtype=float) - np.array(dp2))) - tl2) ** 2
                        for dp2, tl2 in zip(donor_positions, target_lengths)
                    )
                    scored_candidates.append(
                        (
                            -ml_err_fit * 2.0 + min(min_clash_fit, 3.0) * 0.3 + min(min_hm_fit, 4.0) * 0.12,
                            min_clash_fit,
                            np.asarray(metal_pos, dtype=float),
                        )
                    )

                for clash_thresh in (1.8, 1.0):
                    for cand in candidates:
                        min_hm = min(
                            float(np.linalg.norm(cand - _gp(conf, hm)))
                            for hm in hapto_metal_set
                        )
                        if min_hm < 2.0:
                            continue
                        ml_err = sum(
                            (float(np.linalg.norm(cand - np.array(dp2))) - tl2) ** 2
                            for dp2, tl2 in zip(donor_positions, target_lengths)
                        )
                        min_clash = min(
                            (float(np.linalg.norm(cand - _gp(conf, ti)))
                             for ti in trusted if ti not in all_metals),
                            default=999.0,
                        )
                        if min_clash < clash_thresh:
                            continue
                        score = -ml_err * 2.0 + min(min_clash, 3.0) * 0.3
                        score += min(min_hm, 4.0) * 0.12
                        scored_candidates.append((score, min_clash, np.asarray(cand, dtype=float)))
                    if scored_candidates:
                        break

                if scored_candidates:
                    scored_candidates.sort(key=lambda item: (-item[0], -item[1]))
                    unique_candidates: List[Tuple[float, float, object]] = []
                    for score, min_clash, cand in scored_candidates:
                        if any(
                            float(np.linalg.norm(np.asarray(cand, dtype=float) - np.asarray(prev_cand, dtype=float))) < 0.28
                            for _score_prev, _clash_prev, prev_cand in unique_candidates
                        ):
                            continue
                        unique_candidates.append((score, min_clash, cand))
                    chosen_idx = min(variant_rank, len(unique_candidates) - 1)
                    best_score, best_min_clash, best_cand = unique_candidates[chosen_idx]
                    metal_pos = np.asarray(best_cand, dtype=float)
                    logger.info(
                        "SEARCH: selected ranked position %d/%d, score=%.2f, "
                        "min_clash=%.2f, ml_dists=%s",
                        chosen_idx + 1,
                        len(unique_candidates),
                        best_score,
                        best_min_clash,
                        [f"{float(np.linalg.norm(best_cand - np.array(dp))):.2f}"
                         for dp in donor_positions],
                    )

            _sp(conf, metal_idx, metal_pos)
            trusted.add(metal_idx)
            logger.info(
                "Sequential: placed %s%d via %s fit from %d donors, "
                "dist to hapto=%.2f",
                metal_sym, metal_idx, geom_code, len(placed_donors),
                float(np.linalg.norm(metal_pos - hapto_center)),
            )
        elif len(placed_donors) == 1:
            # Single placed donor: place metal along donor→outward direction
            d_pos = np.array(donor_positions[0])
            # Direction: away from hapto metal(s)
            hapto_center = np.mean([_gp(conf, mi) for mi in hapto_metal_set], axis=0)
            direction = d_pos - hapto_center
            norm = np.linalg.norm(direction)
            if norm < 1e-6:
                direction = np.array([1.0, 0.0, 0.0])
            else:
                direction = direction / norm
            _sp(conf, metal_idx, d_pos + direction * target_lengths[0])
            trusted.add(metal_idx)
            logger.info(
                "Sequential: placed %s%d from single donor %s%d at %.2f A",
                metal_sym, metal_idx,
                donor_symbols[0], placed_donors[0], target_lengths[0],
            )

        # Now place unplaced donors around the newly placed metal
        if metal_idx in trusted and unplaced_donors:
            metal_pos = _gp(conf, metal_idx)
            # Compute used directions (from metal to placed donors)
            used_dirs = []
            for d in placed_donors:
                vec = _gp(conf, d) - metal_pos
                n = np.linalg.norm(vec)
                if n > 1e-6:
                    used_dirs.append(vec / n)

            for ud in unplaced_donors:
                ud_sym = mol.GetAtomWithIdx(ud).GetSymbol()
                bl = _get_ml_bond_length(metal_sym, ud_sym)
                # Find direction maximally separated from ALL used directions.
                # For each candidate direction on a sphere, compute the
                # minimum angle to any used direction. Pick the candidate
                # with the largest minimum angle.
                if not used_dirs:
                    new_dir = np.array([0.0, 0.0, 1.0])
                else:
                    best_dir = np.array([0.0, 0.0, 1.0])
                    best_min_angle = -1.0
                    for t_deg in range(0, 180, 15):
                        for p_deg in range(0, 360, 15):
                            t = np.radians(t_deg)
                            p = np.radians(p_deg)
                            cand = np.array([
                                np.sin(t) * np.cos(p),
                                np.sin(t) * np.sin(p),
                                np.cos(t),
                            ])
                            # Minimum angle to any used direction
                            min_ang = min(
                                np.arccos(np.clip(float(np.dot(cand, ud)), -1, 1))
                                for ud in used_dirs
                            )
                            if min_ang > best_min_angle:
                                best_min_angle = min_ang
                                best_dir = cand
                    new_dir = best_dir

                _sp(conf, ud, metal_pos + new_dir * bl)
                used_dirs.append(new_dir)
                trusted.add(ud)

    # ==== Build bonded-pairs set for clash detection ====
    bonded_pairs: set = set()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonded_pairs.add((min(a, b), max(a, b)))
    # Also treat hapto-metal ↔ hapto-carbon as bonded
    for mi, grp in hapto_groups:
        for ci in grp:
            bonded_pairs.add((min(mi, ci), max(mi, ci)))
    # Also 1-3 pairs (atoms 2 bonds apart) — these are NOT clashes
    pairs_13: set = set()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        for nbr_a in mol.GetAtomWithIdx(a).GetNeighbors():
            ni = nbr_a.GetIdx()
            if ni != b:
                pairs_13.add((min(ni, b), max(ni, b)))
        for nbr_b in mol.GetAtomWithIdx(b).GetNeighbors():
            ni = nbr_b.GetIdx()
            if ni != a:
                pairs_13.add((min(ni, a), max(ni, a)))
    bonded_pairs |= pairs_13

    def _min_dist_to_placed(positions_dict, placed_set, skip_set=frozenset()):
        """Return minimum non-bonded distance between new atoms and placed atoms."""
        worst = 999.0
        for new_idx, new_pos in positions_dict.items():
            for placed_idx in placed_set:
                if placed_idx in skip_set:
                    continue
                if placed_idx == new_idx:
                    continue
                pair = (min(new_idx, placed_idx), max(new_idx, placed_idx))
                if pair in bonded_pairs:
                    continue
                placed_pos = _gp(conf, placed_idx)
                d = float(np.linalg.norm(new_pos - placed_pos))
                if d < worst:
                    worst = d
        return worst

    # ==== STEP 4: Embed & place remaining fragments WITH CLASH-AWARE ROTATION ====
    # For fragments connected to secondary metals whose scaffold positions are
    # NOT yet set, we cannot use _align_hybrid_fragment_onto_scaffold (it reads
    # uninitialised scaffold coordinates). Instead:
    #   1. Embed with ETKDG (correct internal geometry)
    #   2. Find donor atoms connecting to the secondary metal
    #   3. Translate fragment so donor centroid sits at correct M-L distance
    #   4. Rotate 360° around metal→donor axis, pick clash-free orientation
    for fragment in decomposition.fragments:
        if all(ai in trusted for ai in fragment.atom_indices):
            continue
        # Must have at least one anchor already placed
        has_anchor = any(ai in trusted for ai in fragment.atom_indices)
        if not has_anchor:
            # Check if this fragment has a metal neighbor that's trusted
            for ai in fragment.atom_indices:
                for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                    if nbr.GetIdx() in trusted and nbr.GetIdx() in all_metals:
                        has_anchor = True
                        break
                if has_anchor:
                    break
        if not has_anchor:
            continue

        embedded = _embed_hybrid_fragment(fragment.fragment_mol)
        if embedded is None:
            continue

        # Find donor atoms in this fragment that connect to a trusted metal
        donor_to_metal: Dict[int, int] = {}  # orig_donor_idx → metal_idx
        for ai in fragment.atom_indices:
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in trusted and ni in all_metals and ni not in fragment.atom_indices:
                    donor_to_metal[ai] = ni

        # Also check anchor atoms already in trusted set
        anchors_in_trusted = [ai for ai in fragment.anchor_atom_indices if ai in trusted]

        # If we have anchor atoms in the scaffold (from Step 2), use standard alignment
        if anchors_in_trusted and not donor_to_metal:
            if _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded):
                # UFF relaxation with fixed anchors
                try:
                    frag_mol = Chem.Mol(fragment.fragment_mol)
                    frag_mol.RemoveAllConformers()
                    frag_conf_obj = Chem.Conformer(frag_mol.GetNumAtoms())
                    for frag_idx, orig_idx in fragment.fragment_to_original.items():
                        p = conf.GetAtomPosition(orig_idx)
                        frag_conf_obj.SetAtomPosition(frag_idx, Point3D(p.x, p.y, p.z))
                    frag_mol.AddConformer(frag_conf_obj, assignId=True)
                    ff = AllChem.UFFGetMoleculeForceField(frag_mol)
                    if ff is not None:
                        for orig_idx in fragment.anchor_atom_indices:
                            frag_idx = fragment.original_to_fragment.get(orig_idx)
                            if frag_idx is not None and orig_idx in trusted:
                                ff.AddFixedPoint(frag_idx)
                        ff.Minimize(maxIts=500)
                        relaxed_conf = frag_mol.GetConformer()
                        for frag_idx, orig_idx in fragment.fragment_to_original.items():
                            if orig_idx in trusted and orig_idx in all_hapto_carbons:
                                continue
                            p = relaxed_conf.GetAtomPosition(frag_idx)
                            conf.SetAtomPosition(orig_idx, Point3D(p.x, p.y, p.z))
                except Exception:
                    pass
                trusted.update(fragment.atom_indices)
                logger.info(
                    "Sequential Step 4: scaffold-aligned fragment (%d atoms)",
                    len(fragment.atom_indices),
                )
                continue

        # ---- Placement for metal-connected fragments ----
        # The donor atoms (e.g. O11, O15) were placed at M-L distance from the
        # metal in Step 3. Use those positions as alignment targets: Procrustes
        # align the ETKDG-embedded fragment so its donors match the Step 3
        # positions. This preserves internal fragment geometry while placing
        # donors at correct M-L distances.
        if not donor_to_metal:
            # Fallback: try standard alignment anyway
            if _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded):
                trusted.update(fragment.atom_indices)
            continue

        donor_indices_list = list(donor_to_metal.keys())
        metal_idx_for_frag = list(donor_to_metal.values())[0]
        metal_pos = _gp(conf, metal_idx_for_frag)
        metal_sym = mol.GetAtomWithIdx(metal_idx_for_frag).GetSymbol()

        # Target positions for donor atoms = their current scaffold positions
        # (set in Step 3's unplaced-donor VSEPR placement)
        target_positions = {d: _gp(conf, d) for d in donor_indices_list}
        for d in donor_indices_list:
            tp = target_positions[d]
            td = float(np.linalg.norm(tp - metal_pos))
            logger.info("Step 4 target: %s%d at %.2f from %s%d, pos=%s",
                        mol.GetAtomWithIdx(d).GetSymbol(), d, td,
                        metal_sym, metal_idx_for_frag, tp)

        # Procrustes alignment of ETKDG embedding to donor targets
        aligned = _align_hybrid_fragment_to_targets(
            scaffold_mol, fragment, embedded,
            target_atom_indices=donor_indices_list,
            target_positions=target_positions,
            reference_center=metal_pos,
            metal_idx=metal_idx_for_frag,
        )
        if not aligned:
            # Fallback: try standard alignment
            logger.info("Step 4: Procrustes alignment FAILED, using scaffold alignment")
            _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded)

        # ---- Snap donors to their target positions ----
        # Procrustes alignment may leave donors at wrong M-L distances
        # because the fragment's donor-donor distance differs from the
        # chelate bite distance. Snap each donor atom to its target
        # position, then UFF-relax the body with donors fixed.
        for d_idx in donor_indices_list:
            tp = target_positions[d_idx]
            _sp(conf, d_idx, np.asarray(tp, dtype=float))
            d_dist = float(np.linalg.norm(np.asarray(tp) - metal_pos))
            logger.info("Step 4 snap: %s%d → %.2f from %s%d",
                        mol.GetAtomWithIdx(d_idx).GetSymbol(), d_idx,
                        d_dist, metal_sym, metal_idx_for_frag)

        # UFF relaxation: fix donor atoms at correct positions, relax
        # body to accommodate the snap (repairs stretched bonds)
        try:
            frag_mol = Chem.Mol(fragment.fragment_mol)
            frag_mol.RemoveAllConformers()
            frag_conf_obj = Chem.Conformer(frag_mol.GetNumAtoms())
            for frag_idx, orig_idx in fragment.fragment_to_original.items():
                p = conf.GetAtomPosition(orig_idx)
                frag_conf_obj.SetAtomPosition(frag_idx, Point3D(p.x, p.y, p.z))
            frag_mol.AddConformer(frag_conf_obj, assignId=True)
            ff = AllChem.UFFGetMoleculeForceField(frag_mol)
            if ff is not None:
                for d_idx in donor_indices_list:
                    frag_idx = fragment.original_to_fragment.get(d_idx)
                    if frag_idx is not None:
                        ff.AddFixedPoint(frag_idx)
                ff.Minimize(maxIts=500)
                relaxed_conf = frag_mol.GetConformer()
                for frag_idx, orig_idx in fragment.fragment_to_original.items():
                    p = relaxed_conf.GetAtomPosition(frag_idx)
                    conf.SetAtomPosition(orig_idx, Point3D(p.x, p.y, p.z))
        except Exception:
            pass

        # ---- Clash-aware rotation around metal→donor axis ----
        # Rotate the ENTIRE fragment as a rigid body (including donors)
        # to avoid breaking donor-to-body bonds.
        all_frag_atoms = list(fragment.atom_indices)
        non_frag_trusted = trusted - set(all_frag_atoms)
        if all_frag_atoms and non_frag_trusted:
            # Rotation axis: metal → donor centroid
            donor_centroid = np.mean([_gp(conf, d) for d in donor_indices_list], axis=0)
            axis = donor_centroid - metal_pos
            if float(np.linalg.norm(axis)) < 1e-6:
                axis = np.array([0.0, 0.0, 1.0])
            pivot = metal_pos
            orig_positions = {ai: _gp(conf, ai) for ai in all_frag_atoms}

            best_angle = 0.0
            best_min_d = _min_dist_to_placed(orig_positions, non_frag_trusted)

            for step in range(1, 36):
                angle = step * (2 * np.pi / 36)
                rotated = _rotate_fragment_around_axis(orig_positions, pivot, axis, angle)
                min_d = _min_dist_to_placed(rotated, non_frag_trusted)
                if min_d > best_min_d:
                    best_min_d = min_d
                    best_angle = angle

            if best_angle != 0.0:
                rotated = _rotate_fragment_around_axis(orig_positions, pivot, axis, best_angle)
                for ai, pos in rotated.items():
                    _sp(conf, ai, pos)

            logger.info(
                "Sequential Step 4: placed fragment (%d atoms) on %s%d, "
                "rot=%.0f°, min_clash=%.2f",
                len(fragment.atom_indices), metal_sym, metal_idx_for_frag,
                np.degrees(best_angle), best_min_d,
            )

        trusted.update(fragment.atom_indices)

    # BFS-propagate any truly remaining atoms (H atoms, small substituents only)
    remaining = set(range(mol.GetNumAtoms())) - trusted
    if remaining:
        try:
            _propagate_non_hapto_atoms(
                scaffold_mol, 0, hapto_groups,
                extra_fixed_indices=trusted,
            )
        except Exception as e:
            logger.debug("Sequential Step 4 BFS fallback: %s", e)
    trusted = set(range(mol.GetNumAtoms()))

    # ==== STEP 5: Build rigid body cache (skip M-L re-fit to preserve layout) ====
    # The metal was placed in Step 3, and fragments in Step 4 — both with
    # clash-aware placement. Re-fitting the metal position now would risk
    # moving it into fragment atoms. M-L distances may not be perfect but
    # topology is preserved.
    stop_set = all_metals | all_hapto_carbons
    frag_body_cache: Dict[int, Set[int]] = {}

    for metal_idx in sorted(non_hapto_metals):
        donor_indices = [
            nbr.GetIdx() for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
            if nbr.GetIdx() not in all_metals and nbr.GetAtomicNum() > 1
        ]
        for d_idx in donor_indices:
            if d_idx in all_hapto_carbons or d_idx in frag_body_cache:
                continue
            body = {d_idx}
            bfs_q_rb = [d_idx]
            while bfs_q_rb:
                curr = bfs_q_rb.pop()
                for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in body or ni in stop_set:
                        continue
                    body.add(ni)
                    bfs_q_rb.append(ni)
            frag_body_cache[d_idx] = body

        logger.info(
            "Sequential Step 5: %s%d kept, dists=%s",
            mol.GetAtomWithIdx(metal_idx).GetSymbol(), metal_idx,
            [f"{float(np.linalg.norm(_gp(conf, d) - _gp(conf, metal_idx))):.2f}" for d in donor_indices],
        )

    # ==== STEP 6: Final fixes ====
    conf = scaffold_mol.GetConformer(0)
    n_atoms = scaffold_mol.GetNumAtoms()

    # 6a: Clash-resolution rotations DISABLED — they compound and create
    # new clashes worse than the original. Step 4 clash-aware rotation
    # handles initial placement; post-hoc rotation just moves problems around.

    # 6b: Fix H atoms with bad parent distances or clashes
    for _h_pass in range(3):
        any_fixed = False
        for i in range(n_atoms):
            if scaffold_mol.GetAtomWithIdx(i).GetAtomicNum() != 1:
                continue
            nbrs = [n.GetIdx() for n in scaffold_mol.GetAtomWithIdx(i).GetNeighbors()]
            if not nbrs:
                continue
            parent_idx = nbrs[0]
            parent_pos = _gp(conf, parent_idx)
            h_pos = _gp(conf, i)
            parent_dist = float(np.linalg.norm(h_pos - parent_pos))

            needs_fix = parent_dist < 0.8 or parent_dist > 1.5
            if not needs_fix:
                for j in range(n_atoms):
                    if j == i or j == parent_idx:
                        continue
                    d = float(np.linalg.norm(h_pos - _gp(conf, j)))
                    if d < (1.0 if scaffold_mol.GetAtomWithIdx(j).GetAtomicNum() > 1 else 1.2):
                        needs_fix = True
                        break
            if not needs_fix:
                continue
            any_fixed = True

            # Place H opposite to parent's other neighbors
            parent_nbrs = [n.GetIdx() for n in scaffold_mol.GetAtomWithIdx(parent_idx).GetNeighbors()
                           if n.GetIdx() != i]
            if parent_nbrs:
                avg_dir = np.zeros(3)
                for pn in parent_nbrs:
                    v = _gp(conf, pn) - parent_pos
                    vn = float(np.linalg.norm(v))
                    if vn > 1e-6:
                        avg_dir += v / vn
                avg_len = float(np.linalg.norm(avg_dir))
                h_dir = -avg_dir / avg_len if avg_len > 1e-6 else np.array([0.0, 0.0, 1.0])
            else:
                h_dir = np.array([0.0, 0.0, 1.0])

            # Try multiple rotations, pick best
            best_pos = parent_pos + h_dir * 1.09
            best_min_d = 0.0
            for j2 in range(n_atoms):
                if j2 == i or j2 == parent_idx:
                    continue
                d2 = float(np.linalg.norm(best_pos - _gp(conf, j2)))
                if best_min_d == 0.0 or d2 < best_min_d:
                    best_min_d = d2

            if best_min_d < 1.0:
                perp = np.cross(h_dir, np.array([1.0, 0.0, 0.0]))
                if float(np.linalg.norm(perp)) < 1e-6:
                    perp = np.cross(h_dir, np.array([0.0, 1.0, 0.0]))
                perp = perp / float(np.linalg.norm(perp))
                for angle in [1.047, 2.094, 3.142, 4.189, 5.236]:
                    cos_a, sin_a = np.cos(angle), np.sin(angle)
                    rot_vec = (h_dir * cos_a
                               + np.cross(perp, h_dir) * sin_a
                               + perp * float(np.dot(perp, h_dir)) * (1 - cos_a))
                    rot_vec = rot_vec / max(float(np.linalg.norm(rot_vec)), 1e-12)
                    cand = parent_pos + rot_vec * 1.09
                    min_d_cand = min(
                        (float(np.linalg.norm(cand - _gp(conf, j3)))
                         for j3 in range(n_atoms)
                         if j3 != i and j3 != parent_idx),
                        default=999.0,
                    )
                    if min_d_cand > best_min_d:
                        best_min_d = min_d_cand
                        best_pos = cand
                        if best_min_d >= 1.0:
                            break

            _sp(conf, i, best_pos)

        if not any_fixed:
            break

    # 6c: Skip pi-coplanarity for multi-metal — it moves atoms and
    # can re-introduce clashes that Step 6a resolved.

    logger.info(
        "Sequential multi-metal: %d/%d atoms placed",
        len(trusted), mol.GetNumAtoms(),
    )
    return scaffold_mol


def _secondary_metal_variant_plans(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
    max_alternates_per_metal: int = 2,
) -> List[Tuple[str, Dict[int, int]]]:
    """Return ranked secondary-metal variant plans for hybrid hapto building."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return [("hybrid", {})]

    hapto_metals = {metal_idx for metal_idx, _grp in hapto_groups}
    plans: List[Tuple[str, Dict[int, int]]] = [("hybrid", {})]
    for atom in mol.GetAtoms():
        metal_idx = atom.GetIdx()
        if atom.GetSymbol() not in _METAL_SET or metal_idx in hapto_metals:
            continue
        donor_count = sum(
            1
            for nbr in atom.GetNeighbors()
            if nbr.GetAtomicNum() > 1 and nbr.GetSymbol() not in _METAL_SET
        )
        if donor_count < 2:
            continue
        for rank in range(1, max(1, int(max_alternates_per_metal)) + 1):
            plans.append(
                (f"hybrid-secondary-{atom.GetSymbol()}{metal_idx}-alt{rank}", {metal_idx: rank})
            )
    return plans


def _build_hybrid_hapto_complex(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
    preview_store: Optional[List[Tuple[str, str]]] = None,
    secondary_variant_plan: Optional[Dict[int, int]] = None,
):
    """Assemble a hapto complex from a scaffold plus rigid ligand fragments."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return None

    decomposition = _decompose_hapto_complex(mol, hapto_groups)
    if decomposition is None:
        return None

    # Multi-metal: use sequential builder
    n_metals = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET)
    if n_metals >= 2 and decomposition is not None:
        try:
            seq_result = _build_multimetal_hapto_sequential(
                mol,
                hapto_groups,
                decomposition,
                secondary_variant_plan=secondary_variant_plan,
            )
        except Exception as e:
            logger.debug("Sequential multi-metal builder failed: %s", e)
            seq_result = None
        if seq_result is not None:
            return seq_result
        logger.info("Sequential multi-metal builder fell back to standard path")

    primary_module = _detect_primary_organometal_module(decomposition, hapto_groups)
    if primary_module is not None:
        try:
            primary_candidate = _build_primary_organometal_hapto_module(
                mol,
                decomposition,
                hapto_groups,
                primary_module,
                preview_store=preview_store,
            )
        except Exception as e:
            logger.debug("Hybrid hapto primary organometal builder skipped: %s", e)
            primary_candidate = None
        if primary_candidate is not None:
            return primary_candidate
        logger.info(
            "Hybrid hapto primary organometal module builder fell back to standard scaffold assembly for %s%d",
            mol.GetAtomWithIdx(primary_module.metal_idx).GetSymbol(),
            primary_module.metal_idx,
        )
        if preview_store is not None and not preview_store:
            _store_hapto_preview_candidate(
                mol,
                hapto_groups,
                preview_store,
                label='quick-hapto-preview',
            )
    elif preview_store is not None:
        _store_hapto_preview_candidate(
            mol,
            hapto_groups,
            preview_store,
            label='quick-hapto-preview',
        )

    scaffold_mol = Chem.Mol(mol)
    if not _build_hapto_scaffold(scaffold_mol, hapto_groups):
        return None

    aligned_fragments = 0
    scaffold_only_fragments = 0
    deferred_secondary_fragments = 0
    trusted_atoms: set = set()
    hapto_metal_indices = {metal_idx for metal_idx, _grp in hapto_groups}
    for fragment in decomposition.fragments:
        if fragment.use_scaffold_only:
            scaffold_only_fragments += 1
            trusted_atoms.update(fragment.atom_indices)
            logger.debug(
                "Hybrid hapto: keeping scaffold geometry for fragment "
                "(atoms=%d, metals=%d, bridging_donors=%d)",
                len(fragment.atom_indices),
                len(fragment.metal_neighbor_indices),
                len(fragment.bridging_donor_indices),
            )
            continue
        if any(metal_idx not in hapto_metal_indices for metal_idx in fragment.metal_neighbor_indices):
            deferred_secondary_fragments += 1
            logger.debug(
                "Hybrid hapto: deferring fragment with secondary-metal coordination "
                "(atoms=%d, metals=%s) to secondary module assembly",
                len(fragment.atom_indices),
                sorted(fragment.metal_neighbor_indices),
            )
            continue
        embedded_fragment = _embed_hybrid_fragment(fragment.fragment_mol)
        if embedded_fragment is None:
            continue
        if _align_hybrid_fragment_onto_scaffold(scaffold_mol, fragment, embedded_fragment):
            aligned_fragments += 1
            trusted_atoms.update(fragment.atom_indices)

    if aligned_fragments == 0 and scaffold_only_fragments == 0:
        return None

    try:
        _correct_hapto_geometry(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    secondary_metals, secondary_module_atoms, relaxable_secondary_atoms = _assemble_secondary_metal_coordination_modules(
        scaffold_mol,
        decomposition,
        hapto_groups,
        fit_variant_plan=secondary_variant_plan,
    )
    trusted_atoms |= secondary_module_atoms

    if not secondary_metals:
        secondary_metals = _place_secondary_metals_in_hapto_fragments(
            scaffold_mol,
            hapto_groups,
        )

    if secondary_metals:
        logger.debug(
            "Hybrid hapto built %d explicit secondary metal module(s)",
            len(secondary_metals),
        )
    else:
        logger.debug("Hybrid hapto used donor-cloud fallback for secondary metals")

    if aligned_fragments + scaffold_only_fragments < len(decomposition.fragments):
        try:
            _propagate_non_hapto_atoms(
                scaffold_mol,
                0,
                hapto_groups,
                extra_fixed_indices=trusted_atoms | _hapto_primary_donor_indices(scaffold_mol, hapto_groups) | secondary_metals,
            )
        except Exception as e:
            logger.debug("Hybrid hapto non-hapto propagation skipped: %s", e)
        if not secondary_metals:
            secondary_metals |= _place_secondary_metals_in_hapto_fragments(
                scaffold_mol,
                hapto_groups,
            )

    try:
        _enforce_donor_pi_coplanarity(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    try:
        _refine_hybrid_hapto_complex(
            scaffold_mol,
            decomposition,
            hapto_groups,
            trusted_atom_indices=trusted_atoms,
            secondary_metal_indices=secondary_metals,
            apply_final_hapto_correction=not bool(secondary_metals),
        )
    except Exception as e:
        logger.debug("Hybrid hapto local refinement skipped: %s", e)

    try:
        _enforce_donor_pi_coplanarity(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    if not secondary_metals:
        try:
            _correct_hapto_geometry(scaffold_mol, 0, hapto_groups)
        except Exception:
            pass

    # Final clash resolution: _enforce_donor_pi_coplanarity may have moved
    # atoms (especially H) into metal clash zones.
    try:
        _final_clash_resolution(scaffold_mol, 0, hapto_groups)
    except Exception:
        pass

    logger.info(
        "Hybrid hapto assembly aligned %d fragment(s), scaffold-kept %d, deferred-secondary %d/%d",
        aligned_fragments,
        scaffold_only_fragments,
        deferred_secondary_fragments,
        len(decomposition.fragments),
    )

    # Bond-length repair for severely distorted non-metal bonds, then quality gate
    try:
        import numpy as np
        conf = scaffold_mol.GetConformer(0)
        _expected = {
            frozenset(['C', 'C']): 1.45, frozenset(['C', 'N']): 1.40,
            frozenset(['C', 'O']): 1.35, frozenset(['C', 'S']): 1.78,
            frozenset(['N', 'N']): 1.40, frozenset(['N', 'O']): 1.35,
            frozenset(['O', 'O']): 1.45, frozenset(['C', 'Si']): 1.87,
        }
        hapto_frozen = set()
        for _mi, catoms in hapto_groups:
            hapto_frozen.add(_mi)
            hapto_frozen.update(catoms)
        for _ai in range(scaffold_mol.GetNumAtoms()):
            if scaffold_mol.GetAtomWithIdx(_ai).GetSymbol() in _METAL_SET:
                hapto_frozen.add(_ai)

        def _gp2(i):
            p = conf.GetAtomPosition(i)
            return np.array([p.x, p.y, p.z])
        def _sp2(i, pos):
            conf.SetAtomPosition(i, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))

        # Build H groups
        h_of2: Dict[int, List[int]] = {}
        for _ai in range(scaffold_mol.GetNumAtoms()):
            a = scaffold_mol.GetAtomWithIdx(_ai)
            if a.GetAtomicNum() == 1:
                for nbr in a.GetNeighbors():
                    if nbr.GetAtomicNum() > 1:
                        h_of2.setdefault(nbr.GetIdx(), []).append(_ai)
                        break

        # Iterative bond-length repair (60 iterations)
        for _rep_pass in range(60):
            max_err = 0.0
            for bond in scaffold_mol.GetBonds():
                bi = bond.GetBeginAtomIdx()
                bj = bond.GetEndAtomIdx()
                si = scaffold_mol.GetAtomWithIdx(bi).GetSymbol()
                sj = scaffold_mol.GetAtomWithIdx(bj).GetSymbol()
                if si in _METAL_SET or sj in _METAL_SET:
                    continue
                pi, pj = _gp2(bi), _gp2(bj)
                d = float(np.linalg.norm(pj - pi))
                if si == 'H' or sj == 'H':
                    expected = 1.08
                else:
                    expected = _expected.get(frozenset([si, sj]), 1.50)
                err = abs(d - expected)
                if err < 0.05:
                    continue
                if err > max_err:
                    max_err = err
                if d < 1e-8:
                    continue
                direction = (pj - pi) / d
                correction = 0.3 * (d - expected)
                bi_frozen = bi in hapto_frozen
                bj_frozen = bj in hapto_frozen
                if bi_frozen and bj_frozen:
                    continue
                if bi_frozen:
                    delta = -correction * direction
                    _sp2(bj, pj + delta)
                    for hi in h_of2.get(bj, []):
                        if hi not in hapto_frozen:
                            _sp2(hi, _gp2(hi) + delta)
                elif bj_frozen:
                    delta = correction * direction
                    _sp2(bi, pi + delta)
                    for hi in h_of2.get(bi, []):
                        if hi not in hapto_frozen:
                            _sp2(hi, _gp2(hi) + delta)
                else:
                    _sp2(bi, pi + 0.5 * correction * direction)
                    _sp2(bj, pj - 0.5 * correction * direction)
            if max_err < 0.1:
                break

        # Re-center secondary metals after bond repair shifted donors
        _fix_secondary_metal_distances(scaffold_mol, 0, hapto_groups)

        # Re-run clash resolution after bond repair + metal recentering
        _final_clash_resolution(scaffold_mol, 0, hapto_groups)

        # Quality gate after repair
        for bond in scaffold_mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            si = scaffold_mol.GetAtomWithIdx(bi).GetSymbol()
            sj = scaffold_mol.GetAtomWithIdx(bj).GetSymbol()
            if si in _METAL_SET or sj in _METAL_SET:
                continue
            pi, pj = _gp2(bi), _gp2(bj)
            d = float(np.linalg.norm(pj - pi))
            if si == 'H' or sj == 'H':
                expected = 1.08
            else:
                expected = _expected.get(frozenset([si, sj]), 1.50)
            ratio = d / expected
            if ratio < 0.50 or ratio > 3.00:
                logger.info(
                    "Hybrid hapto rejected: bond %s[%d]-%s[%d] = %.3f A (expected ~%.2f)",
                    si, bi, sj, bj, d, expected,
                )
                return None
    except Exception:
        pass

    return scaffold_mol


def _final_clash_resolution(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
) -> None:
    """Push apart non-bonded atoms that are too close after all geometry steps.

    Metals and hapto ring atoms are frozen; everything else can move.
    When a heavy atom is pushed, its bonded H atoms move rigidly with it.
    Runs after _enforce_donor_pi_coplanarity which may introduce clashes.
    """
    import numpy as np

    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return

    def _gp(i):
        p = conf.GetAtomPosition(i)
        return np.array([p.x, p.y, p.z])

    def _sp(i, pos):
        conf.SetAtomPosition(i, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))

    def _move_with_h(ai, displacement):
        """Move atom ai and its bonded H atoms by displacement."""
        _sp(ai, _gp(ai) + displacement)
        for hi in h_of.get(ai, []):
            if hi not in frozen:
                _sp(hi, _gp(hi) + displacement)

    frozen = set()
    for mi, catoms in hapto_groups:
        frozen.add(mi)
        frozen.update(catoms)
    for ai in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET:
            frozen.add(ai)

    bonded_pairs = set()
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonded_pairs.add((bi, bj))
        bonded_pairs.add((bj, bi))

    # Build H→parent and parent→H maps
    h_of: Dict[int, List[int]] = {}
    h_parent: Dict[int, int] = {}
    for ai in range(mol.GetNumAtoms()):
        a = mol.GetAtomWithIdx(ai)
        if a.GetAtomicNum() == 1:
            for nbr in a.GetNeighbors():
                if nbr.GetAtomicNum() > 1:
                    h_parent[ai] = nbr.GetIdx()
                    h_of.setdefault(nbr.GetIdx(), []).append(ai)
                    break

    # Only check heavy atoms (H moves with parent)
    heavy = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    n_atoms = mol.GetNumAtoms()
    rng = np.random.default_rng(99)

    for _pass in range(30):
        any_push = False
        # Heavy-heavy clashes
        for ii, i in enumerate(heavy):
            for j in heavy[ii + 1:]:
                if (i, j) in bonded_pairs:
                    continue
                if i in frozen and j in frozen:
                    continue
                pi, pj = _gp(i), _gp(j)
                d = float(np.linalg.norm(pi - pj))
                si = mol.GetAtomWithIdx(i).GetSymbol()
                sj = mol.GetAtomWithIdx(j).GetSymbol()
                if si in _METAL_SET or sj in _METAL_SET:
                    min_d = 2.0
                else:
                    min_d = 1.2
                if d >= min_d:
                    continue
                if d < 1e-8:
                    direction = rng.standard_normal(3)
                    direction /= max(float(np.linalg.norm(direction)), 1e-12)
                else:
                    direction = (pj - pi) / d
                gap = min_d - d
                if i in frozen:
                    _move_with_h(j, gap * direction)
                elif j in frozen:
                    _move_with_h(i, -gap * direction)
                else:
                    _move_with_h(i, -0.5 * gap * direction)
                    _move_with_h(j, 0.5 * gap * direction)
                any_push = True

        # H-metal clashes: move the parent heavy atom (+ all its H) so
        # the offending H clears the metal without breaking C-H bonds.
        for i in range(n_atoms):
            a = mol.GetAtomWithIdx(i)
            if a.GetAtomicNum() != 1 or i in frozen:
                continue
            parent = h_parent.get(i)
            if parent is not None and parent in frozen:
                continue  # parent is frozen, can't fix
            for j in frozen:
                if (i, j) in bonded_pairs:
                    continue
                if mol.GetAtomWithIdx(j).GetSymbol() not in _METAL_SET:
                    continue
                pi, pj = _gp(i), _gp(j)
                d = float(np.linalg.norm(pi - pj))
                if d >= 2.0:
                    continue
                if d < 1e-8:
                    direction = rng.standard_normal(3)
                    direction /= max(float(np.linalg.norm(direction)), 1e-12)
                else:
                    direction = (pi - pj) / d
                displacement = (2.0 - d) * direction
                if parent is not None and parent not in frozen:
                    _move_with_h(parent, displacement)
                else:
                    _sp(i, _gp(i) + displacement)
                any_push = True

        if not any_push:
            break


def _fix_secondary_metal_distances(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
) -> None:
    """Fix non-hapto metal positions by moving them into their donor cloud.

    ETKDG doesn't understand metal-ligand bonds, so metal positions are random.
    For each non-hapto metal, compute the optimal position given its donors'
    current positions and target M-L distances, then move the metal there.
    This avoids moving donor atoms (which may share fragments).
    """
    import numpy as np

    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return

    def _gp(i):
        p = conf.GetAtomPosition(i)
        return np.array([p.x, p.y, p.z])

    def _sp(i, pos):
        conf.SetAtomPosition(i, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))

    hapto_metals = {mi for mi, _ in hapto_groups}

    for atom in mol.GetAtoms():
        mi = atom.GetIdx()
        msym = atom.GetSymbol()
        if msym not in _METAL_SET or mi in hapto_metals:
            continue

        donors = []
        targets = []
        for nbr in atom.GetNeighbors():
            ni = nbr.GetIdx()
            if nbr.GetSymbol() in _METAL_SET or nbr.GetAtomicNum() <= 1:
                continue
            dsym = nbr.GetSymbol()
            donors.append(ni)
            targets.append(float(_get_ml_bond_length(msym, dsym)))

        if len(donors) < 2:
            continue

        # Iteratively move metal toward the position that minimizes
        # sum of (actual_dist - target_dist)^2
        mpos = _gp(mi)
        for _it in range(80):
            grad = np.zeros(3)
            for di, target in zip(donors, targets):
                dpos = _gp(di)
                vec = mpos - dpos
                dist = float(np.linalg.norm(vec))
                if dist < 1e-8:
                    continue
                grad += 2.0 * (dist - target) * vec / dist
            step = -0.15 * grad
            if float(np.linalg.norm(step)) < 0.001:
                break
            mpos = mpos + step

        _sp(mi, mpos)

        # If any donor is still >30% off target, move it toward the metal.
        # Use a limited BFS (max 6 atoms) to avoid dragging shared fragments.
        for di, target in zip(donors, targets):
            dpos = _gp(di)
            dist = float(np.linalg.norm(dpos - mpos))
            if dist < 1e-8 or abs(dist - target) / target < 0.30:
                continue
            direction = (dpos - mpos) / dist
            delta = direction * (target - dist)
            # Limited BFS: donor + up to 5 non-metal, non-hapto neighbors
            local_frag = {di}
            bfs_q = [di]
            while bfs_q and len(local_frag) < 7:
                cur = bfs_q.pop(0)
                for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                    ni2 = nbr.GetIdx()
                    if ni2 in local_frag or ni2 in hapto_metals:
                        continue
                    if nbr.GetSymbol() in _METAL_SET:
                        continue
                    # Don't cross into another donor's territory
                    if ni2 in set(donors) and ni2 != di:
                        continue
                    local_frag.add(ni2)
                    bfs_q.append(ni2)
            for ai in local_frag:
                _sp(ai, _gp(ai) + delta)


def _find_ansa_bridges(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> List[Tuple[int, int, int]]:
    """Find bridging atoms (Si, C, Ge, etc.) connecting two hapto groups on the same metal.

    Returns list of (bridge_atom_idx, group_index_1, group_index_2).
    """
    if not RDKIT_AVAILABLE or mol is None or len(hapto_groups) < 2:
        return []

    # Build map: atom_idx -> list of (group_index, metal_idx)
    atom_to_groups: Dict[int, List[Tuple[int, int]]] = {}
    for gi, (metal_idx, grp) in enumerate(hapto_groups):
        for a_idx in grp:
            atom_to_groups.setdefault(a_idx, []).append((gi, metal_idx))

    bridges: List[Tuple[int, int, int]] = []
    seen_bridges: set = set()
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in _METAL_SET or sym == 'H':
            continue
        a_idx = atom.GetIdx()
        if a_idx in atom_to_groups:
            continue  # skip hapto ring atoms themselves

        # Check if this atom connects to atoms from 2+ different groups on the same metal
        connected_groups: Dict[int, set] = {}  # metal_idx -> set of group indices
        for nbr in atom.GetNeighbors():
            ni = nbr.GetIdx()
            if ni in atom_to_groups:
                for gi, mi in atom_to_groups[ni]:
                    connected_groups.setdefault(mi, set()).add(gi)

        for _mi, group_set in connected_groups.items():
            if len(group_set) >= 2:
                glist = sorted(group_set)
                for i in range(len(glist)):
                    for j in range(i + 1, len(glist)):
                        key = (a_idx, glist[i], glist[j])
                        if key not in seen_bridges:
                            seen_bridges.add(key)
                            bridges.append(key)

    return bridges



def _apply_hapto_centroid_bias(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
    min_centroid_sep: float = 2.5,
) -> bool:
    """Conservative eta-group regularization without aggressive expansion.

    Applies only rigid group translations (no internal ring distortion):
    1) enforce a minimum metal-centroid distance for each eta-group,
    2) apply a light anti-collapse repulsion only when groups are too close.
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return False

    def _centroid(indices: List[int]) -> Tuple[float, float, float]:
        pts = [conf.GetAtomPosition(i) for i in indices]
        return (
            sum(p.x for p in pts) / len(pts),
            sum(p.y for p in pts) / len(pts),
            sum(p.z for p in pts) / len(pts),
        )

    def _norm(vx: float, vy: float, vz: float) -> float:
        return math.sqrt(vx * vx + vy * vy + vz * vz)

    def _translate(indices: List[int], dx: float, dy: float, dz: float) -> None:
        for idx in indices:
            p = conf.GetAtomPosition(idx)
            conf.SetAtomPosition(idx, Point3D(p.x + dx, p.y + dy, p.z + dz))

    def _group_min_distance(a: List[int], b: List[int]) -> float:
        out = float("inf")
        for ai in a:
            pa = conf.GetAtomPosition(ai)
            for bi in b:
                pb = conf.GetAtomPosition(bi)
                d = _norm(pa.x - pb.x, pa.y - pb.y, pa.z - pb.z)
                if d < out:
                    out = d
        return out if math.isfinite(out) else 999.0

    # Normalize/validate groups once.
    group_entries: List[Dict[str, object]] = []
    for metal_idx, grp in hapto_groups:
        if metal_idx < 0 or metal_idx >= mol.GetNumAtoms():
            continue
        grp_valid = sorted(set(i for i in grp if 0 <= i < mol.GetNumAtoms()))
        if len(grp_valid) < 2:
            continue
        group_entries.append({
            "metal": metal_idx,
            "atoms": grp_valid,
        })
    if not group_entries:
        return False

    by_metal: Dict[int, List[int]] = {}
    for gi, ge in enumerate(group_entries):
        by_metal.setdefault(int(ge["metal"]), []).append(gi)

    moved = False

    # Step 1: enforce minimal metal-centroid distance by rigidly pushing groups.
    for metal_idx, gidxs in by_metal.items():
        mpos = conf.GetAtomPosition(metal_idx)
        msym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        for gi in gidxs:
            atoms = list(group_entries[gi]["atoms"])
            cx, cy, cz = _centroid(atoms)
            vx, vy, vz = (cx - mpos.x, cy - mpos.y, cz - mpos.z)
            dist = _norm(vx, vy, vz)
            min_mc = _target_mc_dist(msym, len(atoms))

            if dist < min_mc:
                if dist < 1e-8:
                    ux, uy, uz = (1.0, 0.0, 0.0)
                else:
                    ux, uy, uz = (vx / dist, vy / dist, vz / dist)
                shift = float(min_mc - dist)
                _translate(atoms, ux * shift, uy * shift, uz * shift)
                moved = True

    # Detect ansa bridges for separation limits
    ansa_bridged_pairs: set = set()
    try:
        bridges = _find_ansa_bridges(mol, hapto_groups)
        for _bridge_atom, g1, g2 in bridges:
            ansa_bridged_pairs.add((min(g1, g2), max(g1, g2)))
    except Exception:
        pass

    # Step 2: light anti-collapse guard: only separate groups that are too close.
    safe_min_cc = max(1.85, min(2.20, float(min_centroid_sep)))
    for metal_idx, gidxs in by_metal.items():
        if len(gidxs) < 2:
            continue
        for _iter in range(4):
            changed = False
            centroids: Dict[int, Tuple[float, float, float]] = {}
            for gi in gidxs:
                atoms = list(group_entries[gi]["atoms"])
                centroids[gi] = _centroid(atoms)

            for ii in range(len(gidxs)):
                gi = gidxs[ii]
                ai = set(group_entries[gi]["atoms"])
                ci = centroids[gi]
                for jj in range(ii + 1, len(gidxs)):
                    gj = gidxs[jj]
                    aj = set(group_entries[gj]["atoms"])
                    if ai.intersection(aj):
                        continue

                    cj = centroids[gj]
                    dx, dy, dz = (cj[0] - ci[0], cj[1] - ci[1], cj[2] - ci[2])
                    dcc = _norm(dx, dy, dz)
                    eta_i = len(group_entries[gi]["atoms"])
                    eta_j = len(group_entries[gj]["atoms"])
                    pair_min_cc = safe_min_cc
                    if eta_i >= 4 and eta_j >= 4:
                        pair_min_cc = max(pair_min_cc, 3.00)
                    elif eta_i >= 4 or eta_j >= 4:
                        pair_min_cc = max(pair_min_cc, 2.60)
                    else:  # eta3/eta3
                        pair_min_cc = max(pair_min_cc, 2.30)

                    # Ansa-bridged pairs: limit max centroid separation
                    # Si-C bond ~1.88 A -> max separation ~ 2*1.88 + ring_radii
                    pair_key = (min(gi, gj), max(gi, gj))
                    max_sep = float('inf')
                    if pair_key in ansa_bridged_pairs:
                        max_sep = 5.0  # conservative limit for Si-bridged systems

                    min_cc = _group_min_distance(
                        list(group_entries[gi]["atoms"]),
                        list(group_entries[gj]["atoms"]),
                    )
                    if min_cc >= pair_min_cc:
                        continue

                    if dcc < 1e-8:
                        ux, uy, uz = (1.0, 0.0, 0.0)
                    else:
                        ux, uy, uz = (dx / dcc, dy / dcc, dz / dcc)

                    shift = min(0.45, 0.5 * float(pair_min_cc - min_cc))
                    # Limit shift for ansa-bridged pairs
                    if pair_key in ansa_bridged_pairs and dcc + 2 * shift > max_sep:
                        shift = max(0.0, 0.5 * (max_sep - dcc))
                    _translate(list(group_entries[gi]["atoms"]), -ux * shift, -uy * shift, -uz * shift)
                    _translate(list(group_entries[gj]["atoms"]), ux * shift, uy * shift, uz * shift)
                    changed = True
                    moved = True

            if not changed:
                break

    return moved




def _enforce_hapto_ring_planarity(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
    target_cc: float = 1.40,
    n_iter: int = 5,
) -> bool:
    """Place hapto ring atoms as regular polygon with correct C-C distances.

    For cyclic groups (eta>=3, all atoms have >=2 intra-group neighbors):
      constructs a regular polygon in the plane perpendicular to the
      metal-centroid axis with C-C = target_cc.
    For chain groups: iteratively regularizes C-C distances.
    Substituent atoms (H + direct non-ring neighbors) are rigidly
    translated along with their parent ring atom.
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        import numpy as np
    except ImportError:
        return False
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return False

    moved_substituents: set = set()
    changed = False
    for metal_idx, grp in hapto_groups:
        atoms = sorted(set(i for i in grp if 0 <= i < mol.GetNumAtoms()))
        eta = len(atoms)
        if eta < 3:
            continue
        if metal_idx < 0 or metal_idx >= mol.GetNumAtoms():
            continue

        # Current positions
        pts = np.array([
            [conf.GetAtomPosition(i).x,
             conf.GetAtomPosition(i).y,
             conf.GetAtomPosition(i).z]
            for i in atoms
        ], dtype=float)
        mpos = np.array([
            conf.GetAtomPosition(metal_idx).x,
            conf.GetAtomPosition(metal_idx).y,
            conf.GetAtomPosition(metal_idx).z,
        ])
        centroid = pts.mean(axis=0)

        # Normal = metal -> centroid direction
        mc_vec = centroid - mpos
        mc_dist = float(np.linalg.norm(mc_vec))
        if mc_dist < 1e-8:
            q = pts - centroid
            try:
                _u, _s, vh = np.linalg.svd(q, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            normal = vh[-1].copy()
        else:
            normal = mc_vec / mc_dist
        n_norm = float(np.linalg.norm(normal))
        if n_norm < 1e-12:
            continue
        normal = normal / n_norm

        # In-plane orthonormal basis
        if abs(normal[0]) < 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        else:
            ref = np.array([0.0, 1.0, 0.0])
        u = np.cross(normal, ref)
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        v /= np.linalg.norm(v)

        # Build intra-group adjacency
        atom_set = set(atoms)
        idx_map = {a: i for i, a in enumerate(atoms)}
        adj = [[] for _ in range(eta)]
        for i_local, atom_idx in enumerate(atoms):
            a = mol.GetAtomWithIdx(atom_idx)
            for nbr in a.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set and ni != atom_idx:
                    j_local = idx_map[ni]
                    if j_local not in adj[i_local]:
                        adj[i_local].append(j_local)

        is_cyclic = all(len(adj[i]) >= 2 for i in range(eta))

        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        target_mc = _target_mc_dist(metal_sym, eta)

        if is_cyclic:
            # --- Regular polygon placement for cyclic rings ---
            # Find ring traversal order via graph walk
            order = []
            visited_r = [False] * eta
            cur, prev = 0, -1
            for _ in range(eta):
                order.append(cur)
                visited_r[cur] = True
                nxt = -1
                for nb in adj[cur]:
                    if nb != prev and not visited_r[nb]:
                        nxt = nb
                        break
                if nxt == -1:
                    break
                prev, cur = cur, nxt

            if len(order) != eta:
                # Fallback: angular ordering
                rel = pts - centroid
                angles = np.arctan2(rel @ v, rel @ u)
                order = list(np.argsort(angles))

            # Circumradius for regular polygon with side = target_cc
            R = target_cc / (2.0 * np.sin(np.pi / eta))

            # Starting angle: preserve first atom's angular position
            rel0 = pts[order[0]] - centroid
            start_angle = float(np.arctan2(rel0 @ v, rel0 @ u))

            new_centroid = mpos + normal * target_mc

            pts_new = np.zeros_like(pts)
            for k, idx_in_order in enumerate(order):
                angle = start_angle + 2.0 * np.pi * k / eta
                pts_new[idx_in_order] = (
                    new_centroid
                    + R * np.cos(angle) * u
                    + R * np.sin(angle) * v
                )
        else:
            # --- Chain/non-cyclic: project + iterative regularization ---
            q = pts - centroid
            proj = q - np.outer(q @ normal, normal)
            pts_new = proj + centroid

            edges = []
            for i_local in range(eta):
                for j_local in adj[i_local]:
                    if j_local > i_local:
                        edges.append((i_local, j_local))

            if edges:
                for _it in range(20):
                    for i_l, j_l in edges:
                        pi = pts_new[i_l]
                        pj = pts_new[j_l]
                        vec = pj - pi
                        dist = float(np.linalg.norm(vec))
                        if dist < 1e-8:
                            vec = np.random.default_rng(42).standard_normal(3)
                            dist = float(np.linalg.norm(vec))
                        correction = 0.5 * (dist - target_cc) / dist
                        pts_new[i_l] = pi + correction * vec
                        pts_new[j_l] = pj - correction * vec

                ctr2 = pts_new.mean(axis=0)
                q2 = pts_new - ctr2
                proj2 = q2 - np.outer(q2 @ normal, normal)
                pts_new = proj2 + ctr2

            # Adjust M-centroid distance
            centroid_new = pts_new.mean(axis=0)
            mc_v = centroid_new - mpos
            mc_d = float(np.linalg.norm(mc_v))
            if mc_d > 1e-8:
                shift = (target_mc / mc_d - 1.0) * mc_v
                pts_new = pts_new + shift

        # Write back ring atom coordinates
        new_centroid_final = pts_new.mean(axis=0)
        for i_local, atom_idx in enumerate(atoms):
            conf.SetAtomPosition(
                atom_idx,
                Point3D(float(pts_new[i_local, 0]),
                        float(pts_new[i_local, 1]),
                        float(pts_new[i_local, 2])),
            )

        # Place substituents: correct bond distance + radial outward direction
        all_hapto_atoms = set()
        for _mi, _gi in hapto_groups:
            all_hapto_atoms.update(_gi)

        # Detect bridge atoms (bonded to ring atoms in 2+ different groups)
        bridge_atoms: set = set()
        for pot in range(mol.GetNumAtoms()):
            if pot in all_hapto_atoms:
                continue
            pa = mol.GetAtomWithIdx(pot)
            grps_connected: set = set()
            for pnbr in pa.GetNeighbors():
                pni = pnbr.GetIdx()
                for g_idx, (_gm, g_atoms) in enumerate(hapto_groups):
                    if pni in set(g_atoms):
                        grps_connected.add(g_idx)
            if len(grps_connected) >= 2:
                bridge_atoms.add(pot)

        for i_local, atom_idx in enumerate(atoms):
            ring_pos = pts_new[i_local]
            outward = ring_pos - new_centroid_final
            outward_len = float(np.linalg.norm(outward))
            if outward_len > 1e-8:
                outward_unit = outward / outward_len
            else:
                outward_unit = u  # fallback

            a = mol.GetAtomWithIdx(atom_idx)
            subs = []
            for nbr in a.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set or ni == metal_idx or ni in moved_substituents:
                    continue
                if ni in all_hapto_atoms or ni in bridge_atoms:
                    continue
                subs.append(nbr)

            if not subs:
                continue

            # For single substituent: place radially outward
            # For multiple: spread in a fan around the outward direction
            for s_k, nbr in enumerate(subs):
                ni = nbr.GetIdx()
                nbr_sym = nbr.GetSymbol()
                if nbr_sym == 'H':
                    bond_d = 1.08
                elif nbr_sym == 'Si':
                    bond_d = 1.87
                else:
                    bond_d = 1.50

                if len(subs) == 1:
                    # Single sub: radially outward, tilted away from metal
                    direction = outward_unit + 0.3 * normal
                    d_norm = float(np.linalg.norm(direction))
                    if d_norm > 1e-8:
                        direction = direction / d_norm
                    else:
                        direction = outward_unit
                else:
                    # Multiple subs: spread around outward direction
                    tilt_angle = (s_k - (len(subs) - 1) / 2.0) * 1.2
                    direction = (
                        outward_unit * np.cos(tilt_angle)
                        + normal * np.sin(tilt_angle)
                    )
                    d_norm = float(np.linalg.norm(direction))
                    if d_norm > 1e-8:
                        direction = direction / d_norm

                sub_pos = ring_pos + bond_d * direction
                conf.SetAtomPosition(
                    ni,
                    Point3D(float(sub_pos[0]),
                            float(sub_pos[1]),
                            float(sub_pos[2])),
                )
                moved_substituents.add(ni)
        changed = True

    # --- Post-processing: ansa bridge constraint + bridge atom placement ---
    if changed:
        try:
            import numpy as np
            conf = mol.GetConformer(conf_id)
        except Exception:
            return changed

        all_hapto_set = set()
        group_of_atom: dict = {}
        for g_idx, (_mi, _gi) in enumerate(hapto_groups):
            all_hapto_set.update(_gi)
            for ai in _gi:
                group_of_atom[ai] = g_idx

        # Find bridge atoms connecting different hapto groups
        bridges = []
        bridge_atom_set: set = set()
        for pot in range(mol.GetNumAtoms()):
            if pot in all_hapto_set:
                continue
            pa = mol.GetAtomWithIdx(pot)
            connections = []
            for pnbr in pa.GetNeighbors():
                pni = pnbr.GetIdx()
                if pni in group_of_atom:
                    connections.append((pni, group_of_atom[pni]))
            groups_seen = set(gi for _, gi in connections)
            if len(groups_seen) >= 2:
                bridges.append((pot, connections))
                bridge_atom_set.add(pot)

        # Collect movable atoms per group (ring + substituents, excluding bridges)
        group_atom_sets: dict = {}
        for g_idx, (_mi, _gi) in enumerate(hapto_groups):
            g_set = set(_gi)
            for ai in _gi:
                a = mol.GetAtomWithIdx(ai)
                for nbr in a.GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni not in all_hapto_set and ni != _mi and ni not in bridge_atom_set:
                        g_set.add(ni)
            group_atom_sets[g_idx] = g_set

        def _get_pos(idx):
            p = conf.GetAtomPosition(idx)
            return np.array([p.x, p.y, p.z])

        def _set_pos(idx, arr):
            conf.SetAtomPosition(idx, Point3D(float(arr[0]),
                                               float(arr[1]),
                                               float(arr[2])))

        def _rodrigues_rotate(points, center, axis, angle):
            k = axis / np.linalg.norm(axis)
            cos_a, sin_a = np.cos(angle), np.sin(angle)
            result = []
            for p in points:
                v = p - center
                v_rot = (v * cos_a + np.cross(k, v) * sin_a
                         + k * np.dot(k, v) * (1.0 - cos_a))
                result.append(center + v_rot)
            return result

        def _orient_bridge_atom_inward(ring_atoms, all_atoms, c_bridge,
                                        centroid_self, centroid_other, normal):
            """Rotate ring around its normal so c_bridge faces toward other ring."""
            toward = centroid_other - centroid_self
            toward_proj = toward - np.dot(toward, normal) * normal
            tp_len = float(np.linalg.norm(toward_proj))
            if tp_len < 0.1:
                return  # can't determine direction
            toward_proj /= tp_len

            v_bridge = _get_pos(c_bridge) - centroid_self
            v_proj = v_bridge - np.dot(v_bridge, normal) * normal
            vp_len = float(np.linalg.norm(v_proj))
            if vp_len < 1e-8:
                return
            v_proj /= vp_len

            cos_r = float(np.clip(np.dot(v_proj, toward_proj), -1.0, 1.0))
            sin_r = float(np.dot(np.cross(v_proj, toward_proj), normal))
            rot_angle = np.arctan2(sin_r, cos_r)
            if abs(rot_angle) < 0.01:
                return

            pts = [_get_pos(a) for a in all_atoms]
            pts_rot = _rodrigues_rotate(pts, centroid_self, normal, rot_angle)
            for ai, new_p in zip(all_atoms, pts_rot):
                _set_pos(ai, new_p)

        # Process each ansa bridge using analytical tilt angle
        for bridge_idx, connections in bridges:
            bridge_sym = mol.GetAtomWithIdx(bridge_idx).GetSymbol()
            if bridge_sym == 'Si':
                bridge_bond = 1.87
                bridge_angle_deg = 93.0
            elif bridge_sym == 'C':
                bridge_bond = 1.54
                bridge_angle_deg = 109.5
            else:
                bridge_bond = 1.80
                bridge_angle_deg = 100.0
            d_target = 2.0 * bridge_bond * np.sin(np.radians(bridge_angle_deg) / 2.0)

            by_group: dict = {}
            for ring_atom, g_idx in connections:
                by_group.setdefault(g_idx, []).append(ring_atom)
            group_indices = list(by_group.keys())
            if len(group_indices) < 2:
                continue
            gi, gj = group_indices[0], group_indices[1]
            c_a = by_group[gi][0]
            c_b = by_group[gj][0]
            metal_idx = hapto_groups[gi][0]
            metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
            grp_i = hapto_groups[gi][1]
            grp_j = hapto_groups[gj][1]
            atoms_i = list(group_atom_sets.get(gi, set(grp_i)))
            atoms_j = list(group_atom_sets.get(gj, set(grp_j)))

            eta_i = len(grp_i)
            eta_j = len(grp_j)
            r_i = _target_mc_dist(metal_sym, eta_i)
            r_j = _target_mc_dist(metal_sym, eta_j)
            R_i = target_cc / (2.0 * np.sin(np.pi / max(eta_i, 3)))
            R_j = target_cc / (2.0 * np.sin(np.pi / max(eta_j, 3)))
            r = (r_i + r_j) / 2.0
            R = (R_i + R_j) / 2.0

            # Analytical tilt: psi = arctan(R/r) + arcsin(d_target/(2*sqrt(r²+R²)))
            hyp = np.sqrt(r ** 2 + R ** 2)
            ratio = d_target / (2.0 * hyp)
            if abs(ratio) > 1.0:
                continue  # bridge too long for ring geometry
            psi_target = np.arctan2(R, r) + np.arcsin(ratio)
            alpha_target = 2.0 * psi_target

            pos_m = _get_pos(metal_idx)
            cent_i = np.mean([_get_pos(a) for a in grp_i], axis=0)
            cent_j = np.mean([_get_pos(a) for a in grp_j], axis=0)
            n_i = cent_i - pos_m
            n_j = cent_j - pos_m
            ni_len = float(np.linalg.norm(n_i))
            nj_len = float(np.linalg.norm(n_j))
            if ni_len < 1e-8 or nj_len < 1e-8:
                continue
            n_i /= ni_len
            n_j /= nj_len

            cos_alpha_cur = float(np.dot(n_i, n_j))
            alpha_current = np.arccos(np.clip(cos_alpha_cur, -1.0, 1.0))

            delta_tilt = (alpha_current - alpha_target) / 2.0
            if delta_tilt < 0.01:
                pass  # skip tilt, go straight to orient + place bridge
            else:
                # Tilt axis: perpendicular to plane of M, cent_i, cent_j
                rot_axis = np.cross(n_i, n_j)
                ra_len = float(np.linalg.norm(rot_axis))
                if ra_len < 1e-8:
                    if abs(n_i[0]) < 0.9:
                        rot_axis = np.cross(n_i, np.array([1.0, 0.0, 0.0]))
                    else:
                        rot_axis = np.cross(n_i, np.array([0.0, 1.0, 0.0]))
                    ra_len = float(np.linalg.norm(rot_axis))
                rot_axis /= ra_len

                # Check direction: +delta should increase n_i·n_j (reduce angle)
                test_ci = _rodrigues_rotate([cent_i], pos_m, rot_axis, 0.01)[0]
                test_ni = test_ci - pos_m
                test_ni /= np.linalg.norm(test_ni)
                if float(np.dot(test_ni, n_j)) < cos_alpha_cur:
                    rot_axis = -rot_axis

                # Apply single tilt to ring_i (+delta_tilt)
                pts_i = [_get_pos(a) for a in atoms_i]
                pts_i_rot = _rodrigues_rotate(pts_i, pos_m, rot_axis, +delta_tilt)
                for ai, new_p in zip(atoms_i, pts_i_rot):
                    _set_pos(ai, new_p)

                # Apply single tilt to ring_j (-delta_tilt)
                pts_j = [_get_pos(a) for a in atoms_j]
                pts_j_rot = _rodrigues_rotate(pts_j, pos_m, rot_axis, -delta_tilt)
                for aj, new_p in zip(atoms_j, pts_j_rot):
                    _set_pos(aj, new_p)

            # Orient bridge atoms to face each other (now that rings are tilted)
            cent_i = np.mean([_get_pos(a) for a in grp_i], axis=0)
            cent_j = np.mean([_get_pos(a) for a in grp_j], axis=0)
            n_i = cent_i - pos_m
            ni_len = float(np.linalg.norm(n_i))
            n_j = cent_j - pos_m
            nj_len = float(np.linalg.norm(n_j))
            if ni_len > 1e-8:
                n_i /= ni_len
            if nj_len > 1e-8:
                n_j /= nj_len
            _orient_bridge_atom_inward(
                grp_i, atoms_i, c_a, cent_i, cent_j, n_i)
            _orient_bridge_atom_inward(
                grp_j, atoms_j, c_b, cent_j, cent_i, n_j)

            # --- Place bridge atom with correct geometry ---
            pos_a = _get_pos(c_a)
            pos_b = _get_pos(c_b)
            mid = (pos_a + pos_b) / 2.0
            d_ab = float(np.linalg.norm(pos_a - pos_b))
            half_d = d_ab / 2.0
            h_sq = bridge_bond ** 2 - half_d ** 2
            h = float(np.sqrt(max(h_sq, 0.01)))

            away = mid - _get_pos(metal_idx)
            a_norm = float(np.linalg.norm(away))
            if a_norm > 1e-8:
                away = away / a_norm
            else:
                away = np.array([0.0, 1.0, 0.0])
            bridge_pos = mid + h * away
            _set_pos(bridge_idx, bridge_pos)

            # Place bridge substituents in local frame
            pa = mol.GetAtomWithIdx(bridge_idx)
            ca_cb = pos_b - pos_a
            x_ax = ca_cb / max(d_ab, 1e-8)
            z_ax = away
            y_ax = np.cross(z_ax, x_ax)
            y_norm = float(np.linalg.norm(y_ax))
            if y_norm > 1e-8:
                y_ax = y_ax / y_norm
            else:
                y_ax = np.array([0.0, 0.0, 1.0])

            sub_count = 0
            for pnbr in pa.GetNeighbors():
                pni = pnbr.GetIdx()
                if pni in all_hapto_set:
                    continue
                if pnbr.GetSymbol() == 'C':
                    sub_bond = 1.87
                elif pnbr.GetSymbol() == 'H':
                    sub_bond = 1.48
                else:
                    sub_bond = 1.50
                sub_angle = (sub_count - 0.5) * 2.1
                direction = (np.cos(sub_angle) * y_ax
                             + np.sin(sub_angle) * z_ax)
                d_norm = float(np.linalg.norm(direction))
                if d_norm > 1e-8:
                    direction /= d_norm
                sub_pos = bridge_pos + sub_bond * direction
                _set_pos(pni, sub_pos)
                sub_count += 1

    return changed


def _correct_hapto_geometry(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
    target_cc: float = 1.40,
) -> bool:
    """Topology-preserving hapto geometry correction for ETKDG embeddings.

    Uses rigid-body transforms on connected fragments to fix M-centroid
    distances and ring orientations while preserving all non-metal bond
    lengths and angles from the ETKDG embedding.

    Steps:
      1. Identify fragment for each hapto group (BFS from ring atoms,
         excluding metal and other groups' ring atoms).
      2. For groups on the same metal, enforce angular separation
         (opposite sides for 2 groups, ~120° for 3, etc.).
      3. Rigid-body translate each fragment so M-centroid = target distance.
      4. Flatten ring atoms onto plane perpendicular to M-centroid axis (SVD).
      5. Gentle C-C regularization via spring iterations.
      6. For ansa-bridged systems, adjust tilt and place bridge atoms.
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        import numpy as np
    except ImportError:
        return False
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return False

    def _gp(i):
        p = conf.GetAtomPosition(i)
        return np.array([p.x, p.y, p.z])

    def _sp(i, arr):
        conf.SetAtomPosition(i, Point3D(float(arr[0]), float(arr[1]),
                                         float(arr[2])))

    # -- Collect all hapto atoms and bridge atoms --
    all_hapto = set()
    group_of = {}
    for gi, (mi, catoms) in enumerate(hapto_groups):
        for a in catoms:
            all_hapto.add(a)
            group_of[a] = gi

    metals = set(mi for mi, _ in hapto_groups)
    # Also include ALL metal atoms (not just hapto metals) so that
    # BFS fragments don't cross into other metals' coordination spheres.
    for ai in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET:
            metals.add(ai)

    bridge_atoms: set = set()
    for ai in range(mol.GetNumAtoms()):
        if ai in all_hapto or ai in metals:
            continue
        a = mol.GetAtomWithIdx(ai)
        groups_seen = set()
        for nbr in a.GetNeighbors():
            ni = nbr.GetIdx()
            if ni in group_of:
                groups_seen.add(group_of[ni])
        if len(groups_seen) >= 2:
            bridge_atoms.add(ai)

    # -- Build exclusion set: atoms in rings containing sigma donors to
    #    other metals.  These must NOT be pulled along when rotating
    #    a hapto fragment, as they belong to another metal's coordination
    #    sphere (e.g. pyridine ring bonded to Pt in an Fe/Pt complex). --
    sigma_ring_exclude: set = set()
    ri = mol.GetRingInfo()
    try:
        all_rings = ri.AtomRings()
    except Exception:
        all_rings = []
    for m_idx in metals:
        for nbr in mol.GetAtomWithIdx(m_idx).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in all_hapto or ni in metals:
                continue
            # ni is a sigma donor to metal m_idx
            # Exclude all ring atoms of rings containing ni
            for ring in all_rings:
                if ni in ring:
                    sigma_ring_exclude.update(ring)

    # -- Find fragment for each group (BFS, excluding metal/other groups/bridges) --
    fragments: List[set] = []
    for gi, (mi, catoms) in enumerate(hapto_groups):
        frag = set(catoms)
        queue = list(catoms)
        while queue:
            cur = queue.pop(0)
            for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in frag or ni in metals or ni in bridge_atoms:
                    continue
                if ni in all_hapto and group_of.get(ni) != gi:
                    continue
                # Exclude sigma donors to other metals and their ring atoms
                if ni in sigma_ring_exclude:
                    continue
                bonded_to_other = False
                for nn in mol.GetAtomWithIdx(ni).GetNeighbors():
                    nni = nn.GetIdx()
                    if nni in metals and nni != mi:
                        bonded_to_other = True
                        break
                if bonded_to_other:
                    continue
                frag.add(ni)
                queue.append(ni)
        fragments.append(frag)

    # -- Group hapto entries by metal --
    by_metal: Dict[int, List[int]] = {}
    for gi, (mi, _) in enumerate(hapto_groups):
        by_metal.setdefault(mi, []).append(gi)

    changed = False

    # =====================================================================
    # STEP 1: Place hapto groups at analytically ideal positions on the
    # coordination sphere.  Uses rigid-body rotation around the metal
    # centre (preserves ALL internal distances within each fragment).
    #
    # 2 groups without ansa bridge → sandwich (~175°)
    # 2 groups with ansa bridge   → bent (~130°)
    # 3+ groups                   → evenly distributed in a plane
    # =====================================================================

    def _rodrigues(v, k, angle):
        """Rodrigues rotation of vector *v* around unit axis *k*."""
        ca, sa = np.cos(angle), np.sin(angle)
        return v * ca + np.cross(k, v) * sa + k * np.dot(k, v) * (1.0 - ca)

    def _rotate_fragment(frag_indices, centre, cur_dir, tgt_dir):
        """Rigid-body rotation of fragment so *cur_dir* maps to *tgt_dir*."""
        cd = np.asarray(cur_dir, dtype=float)
        td = np.asarray(tgt_dir, dtype=float)
        cos_a = float(np.clip(np.dot(cd, td), -1.0, 1.0))
        if cos_a > 0.9999:
            return  # already aligned
        cross = np.cross(cd, td)
        cl = float(np.linalg.norm(cross))
        if cl < 1e-8:
            # Anti-parallel – pick any perpendicular axis
            perp = np.array([1, 0, 0]) if abs(cd[0]) < 0.9 else np.array([0, 1, 0])
            cross = np.cross(cd, perp)
            cl = float(np.linalg.norm(cross))
        k = cross / cl
        rot_angle = np.arccos(cos_a)
        for ai in frag_indices:
            v = _gp(ai) - centre
            _sp(ai, centre + _rodrigues(v, k, rot_angle))

    def _groups_bridged(gi, gj):
        """True if groups *gi* and *gj* share an ansa bridge atom."""
        for bi in bridge_atoms:
            gs = set()
            for nbr in mol.GetAtomWithIdx(bi).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in group_of:
                    gs.add(group_of[ni])
            if gi in gs and gj in gs:
                return True
        return False

    for mi, gidxs in by_metal.items():
        if len(gidxs) < 2:
            continue
        mpos = _gp(mi)
        n_grp = len(gidxs)

        # Current centroid unit-directions from metal
        cur_dirs: Dict[int, np.ndarray] = {}
        for gi in gidxs:
            c = np.mean([_gp(a) for a in hapto_groups[gi][1]], axis=0)
            v = c - mpos
            d = float(np.linalg.norm(v))
            cur_dirs[gi] = v / d if d > 1e-8 else np.array([1.0, 0, 0])

        # --- Compute ideal directions on sphere ---
        ideal_dirs: Dict[int, np.ndarray] = {}

        if n_grp == 2:
            gi, gj = gidxs
            bridged = _groups_bridged(gi, gj)
            target_angle = np.radians(130.0 if bridged else 175.0)

            d_i = cur_dirs[gi]
            ideal_dirs[gi] = d_i  # keep group i fixed

            # Rotation axis: perpendicular to d_i in the (d_i, d_j) plane
            cross = np.cross(d_i, cur_dirs[gj])
            cl = float(np.linalg.norm(cross))
            if cl < 1e-8:
                perp = np.array([1, 0, 0]) if abs(d_i[0]) < 0.9 else np.array([0, 1, 0])
                cross = np.cross(d_i, perp)
                cl = float(np.linalg.norm(cross))
            rot_ax = cross / cl
            ideal_dirs[gj] = _rodrigues(d_i, rot_ax, target_angle)

        elif n_grp >= 3:
            # Best-fit plane through current centroid directions
            pts = np.array([cur_dirs[gi] for gi in gidxs])
            try:
                _, _, vh = np.linalg.svd(
                    pts - pts.mean(axis=0), full_matrices=False)
                px = vh[0] / np.linalg.norm(vh[0])
                py = vh[1] / np.linalg.norm(vh[1])
            except np.linalg.LinAlgError:
                px = np.array([1, 0, 0])
                py = np.array([0, 1, 0])

            # Starting angle from first group's projection
            start = np.arctan2(
                float(np.dot(cur_dirs[gidxs[0]], py)),
                float(np.dot(cur_dirs[gidxs[0]], px)),
            )
            for k, gi in enumerate(gidxs):
                a = start + 2.0 * np.pi * k / n_grp
                d = np.cos(a) * px + np.sin(a) * py
                dl = float(np.linalg.norm(d))
                ideal_dirs[gi] = d / dl if dl > 1e-8 else px

        # --- Rotate each fragment to its ideal direction ---
        for gi in gidxs:
            if gi not in ideal_dirs:
                continue
            _rotate_fragment(fragments[gi], mpos, cur_dirs[gi], ideal_dirs[gi])
            changed = True

    # =====================================================================
    # STEP 2: Rigid-body translate fragments to target M-centroid distance
    # =====================================================================
    for gi, (mi, catoms) in enumerate(hapto_groups):
        mpos = _gp(mi)
        centroid = np.mean([_gp(a) for a in catoms], axis=0)
        mc_vec = centroid - mpos
        mc_dist = float(np.linalg.norm(mc_vec))
        if mc_dist < 1e-8:
            mc_dir = np.array([1.0, 0.0, 0.0])
            mc_dist = 1e-8
        else:
            mc_dir = mc_vec / mc_dist

        metal_sym = mol.GetAtomWithIdx(mi).GetSymbol()
        target_mc = _target_mc_dist(metal_sym, len(catoms))

        shift = mc_dir * (target_mc - mc_dist)
        if abs(target_mc - mc_dist) > 0.01:
            for ai in fragments[gi]:
                _sp(ai, _gp(ai) + shift)
            changed = True

    # =====================================================================
    # STEP 3: Rotate ring plane to be perpendicular to M-centroid axis
    # Uses rigid-body rotation (preserves ALL internal distances).
    # =====================================================================
    for gi, (mi, catoms) in enumerate(hapto_groups):
        atoms = sorted(set(catoms))
        eta = len(atoms)
        if eta < 3:
            continue
        mpos = _gp(mi)
        pts = np.array([_gp(a) for a in atoms])
        centroid = pts.mean(axis=0)
        mc_vec = centroid - mpos
        mc_dist = float(np.linalg.norm(mc_vec))
        if mc_dist < 1e-8:
            continue
        target_normal = mc_vec / mc_dist

        # Find current ring plane normal via SVD
        q = pts - centroid
        try:
            _u, _s, vh = np.linalg.svd(q, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        ring_normal = vh[-1].copy()
        # Ensure ring_normal points same way as target_normal
        if np.dot(ring_normal, target_normal) < 0:
            ring_normal = -ring_normal

        # Compute rotation from ring_normal to target_normal
        cos_a = float(np.clip(np.dot(ring_normal, target_normal), -1, 1))
        if cos_a > 0.9999:
            continue  # already aligned
        rot_axis = np.cross(ring_normal, target_normal)
        ra_len = float(np.linalg.norm(rot_axis))
        if ra_len < 1e-10:
            continue
        rot_axis /= ra_len
        angle = np.arccos(cos_a)
        cos_r = np.cos(angle)
        sin_r = np.sin(angle)
        k = rot_axis

        # Apply Rodrigues rotation to entire fragment (centered at centroid)
        for ai in fragments[gi]:
            v = _gp(ai) - centroid
            v_rot = (v * cos_r + np.cross(k, v) * sin_r
                     + k * np.dot(k, v) * (1.0 - cos_r))
            _sp(ai, centroid + v_rot)
        changed = True

    # =====================================================================
    # STEP 4: Ring shape correction.
    # For cyclic rings: reconstruct as regular polygon (optimal geometry
    # for Cp/arene hapto groups).  For non-cyclic: gentle spring correction.
    # Only ring atoms + their direct H substituents are moved.
    # Polygon starting angle optimized to minimize clashes with environment.
    # =====================================================================
    # Determine max number of groups on a single metal (for adaptive behavior)
    max_groups_per_metal = max(
        (len(gidxs) for gidxs in by_metal.values()), default=1
    )

    for gi, (mi, catoms) in enumerate(hapto_groups):
        atoms = sorted(set(catoms))
        eta = len(atoms)
        if eta < 3:
            continue
        atom_set = set(atoms)

        # Build edges (intra-ring bonds)
        edges = []
        for ai in atoms:
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set and ni > ai:
                    edges.append((ai, ni))

        if not edges:
            continue

        # Collect H substituents for each ring atom
        h_subs: Dict[int, List[int]] = {}
        for ai in atoms:
            hs = []
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if nbr.GetSymbol() == 'H' and ni not in all_hapto:
                    hs.append(ni)
            h_subs[ai] = hs

        # Reconstruct cyclic rings as regular polygons.
        # Choose starting angle to minimize clashes with non-ring atoms.
        mpos = _gp(mi)
        ring_pts = np.array([_gp(a) for a in atoms])
        ring_centroid = ring_pts.mean(axis=0)
        mc_vec = ring_centroid - mpos
        mc_norm = float(np.linalg.norm(mc_vec))
        if mc_norm < 1e-8:
            continue
        plane_normal = mc_vec / mc_norm

        # In-plane basis
        if abs(plane_normal[0]) < 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        else:
            ref = np.array([0.0, 1.0, 0.0])
        u = np.cross(plane_normal, ref)
        u /= np.linalg.norm(u)
        v = np.cross(plane_normal, u)
        v /= np.linalg.norm(v)

        # Build adjacency for ring traversal
        adj_r = [[] for _ in range(eta)]
        idx_map = {a: i for i, a in enumerate(atoms)}
        for ai in atoms:
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set and ni != ai:
                    jl = idx_map[ni]
                    il = idx_map[ai]
                    if jl not in adj_r[il]:
                        adj_r[il].append(jl)

        is_cyclic = all(len(adj_r[i]) >= 2 for i in range(eta))

        if not is_cyclic:
            # Non-cyclic groups (rare): gentle spring-based C-C correction
            for _it in range(30):
                max_err = 0.0
                for ai, aj in edges:
                    pi = _gp(ai)
                    pj = _gp(aj)
                    vec = pj - pi
                    d = float(np.linalg.norm(vec))
                    if d < 1e-8:
                        continue
                    err = d - target_cc
                    max_err = max(max_err, abs(err))
                    c = 0.15 * err / d * vec
                    c = c - np.dot(c, plane_normal) * plane_normal
                    _sp(ai, pi + c)
                    _sp(aj, pj - c)
                if max_err < 0.1:
                    break
            changed = True
            continue

        # -- Cyclic rings: always reconstruct as regular polygon --
        # For each ring atom, find its "branch" — all fragment atoms
        # reachable exclusively through that ring atom.  These move
        # rigidly with the ring atom (preserves substituent bond lengths).

        # Build branch for each ring atom (BFS into fragment, not
        # crossing other ring atoms or metal or bridges)
        branches: Dict[int, List[int]] = {a: [] for a in atoms}
        for ai in atoms:
            visited_br: set = set(atoms) | metals | bridge_atoms | sigma_ring_exclude
            queue_br = []
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if ni not in visited_br:
                    queue_br.append(ni)
                    visited_br.add(ni)
            while queue_br:
                cur_br = queue_br.pop(0)
                branches[ai].append(cur_br)
                for nbr in mol.GetAtomWithIdx(cur_br).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni not in visited_br:
                        visited_br.add(ni)
                        queue_br.append(ni)

        # Ring traversal order
        order = []
        visited_r = [False] * eta
        cur, prev = 0, -1
        for _ in range(eta):
            order.append(cur)
            visited_r[cur] = True
            nxt = -1
            for nb in adj_r[cur]:
                if nb != prev and not visited_r[nb]:
                    nxt = nb
                    break
            if nxt == -1:
                break
            prev, cur = cur, nxt

        if len(order) != eta:
            rel = ring_pts - ring_centroid
            angles_f = np.arctan2(rel @ v, rel @ u)
            order = list(np.argsort(angles_f))

        # Circumradius for regular polygon
        R = target_cc / (2.0 * np.sin(np.pi / eta))

        # Collect atom indices that belong to THIS group's branches
        # (these move with the ring, so exclude from clash scoring)
        own_branch_atoms = set(atoms)
        for ai in atoms:
            own_branch_atoms.update(branches[ai])

        # Collect non-ring, non-metal, non-own-branch heavy atom positions
        # for clash scoring
        env_atoms = []
        for eidx in range(mol.GetNumAtoms()):
            if eidx in own_branch_atoms or eidx in metals or eidx in bridge_atoms:
                continue
            if mol.GetAtomWithIdx(eidx).GetSymbol() == 'H':
                continue
            env_atoms.append(_gp(eidx))

        # Try multiple starting angles, pick the one with fewest clashes
        best_pts = None
        best_clash_score = float('inf')
        n_angles = 12
        for a_try in range(n_angles):
            start_angle = 2.0 * np.pi * a_try / n_angles
            pts_try = np.zeros((eta, 3))
            for k, idx_in_order in enumerate(order):
                angle = start_angle + 2.0 * np.pi * k / eta
                pts_try[idx_in_order] = (
                    ring_centroid
                    + R * np.cos(angle) * u
                    + R * np.sin(angle) * v
                )
            # Score: sum of 1/d for close contacts with environment
            clash_score = 0.0
            for ep in env_atoms:
                for pi in pts_try:
                    d = float(np.linalg.norm(pi - ep))
                    if d < 1.5:
                        clash_score += (1.5 - d) ** 2
            if clash_score < best_clash_score:
                best_clash_score = clash_score
                best_pts = pts_try

        if best_pts is None:
            continue

        # Apply polygon positions, moving entire branches with each ring atom
        for i_loc, ai in enumerate(atoms):
            delta = best_pts[i_loc] - _gp(ai)
            if float(np.linalg.norm(delta)) < 0.001:
                continue
            _sp(ai, best_pts[i_loc])
            for bi in branches[ai]:
                _sp(bi, _gp(bi) + delta)

        # Orient branches outward: if a branch root atom is closer to the
        # ring centroid than its ring parent, reflect it across the parent
        # in the outward direction (away from centroid).
        for ai in atoms:
            if not branches[ai]:
                continue
            ring_pos = _gp(ai)
            outward = ring_pos - ring_centroid
            outward_len = float(np.linalg.norm(outward))
            if outward_len < 1e-8:
                continue
            outward_unit = outward / outward_len
            # Check first branch atom (root of substituent)
            root = branches[ai][0]
            root_pos = _gp(root)
            branch_dir = root_pos - ring_pos
            # If branch points inward (toward centroid), reflect it
            if float(np.dot(branch_dir, outward_unit)) < 0:
                # Reflect all branch atoms across the plane perpendicular
                # to outward through the ring atom
                for bi in branches[ai]:
                    bp = _gp(bi)
                    v = bp - ring_pos
                    proj = float(np.dot(v, outward_unit))
                    # Reflect across plane perpendicular to outward
                    _sp(bi, bp - 2.0 * proj * outward_unit)

        changed = True

    # =====================================================================
    # STEP 5: Ansa bridge handling
    # =====================================================================
    if bridge_atoms and changed:
        for bridge_idx in bridge_atoms:
            ba = mol.GetAtomWithIdx(bridge_idx)
            bridge_sym = ba.GetSymbol()
            if bridge_sym == 'Si':
                bridge_bond = 1.87
            elif bridge_sym == 'C':
                bridge_bond = 1.54
            else:
                bridge_bond = 1.80

            # Find ring atoms this bridge connects
            ring_nbrs = []
            for nbr in ba.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in all_hapto:
                    ring_nbrs.append((ni, group_of[ni]))
            if len(ring_nbrs) < 2:
                continue

            # Place bridge at midpoint of connected ring atoms, pushed outward
            ring_pos = [_gp(rn) for rn, _ in ring_nbrs]
            mid = np.mean(ring_pos, axis=0)

            # Find the metal for outward direction
            mi_bridge = hapto_groups[ring_nbrs[0][1]][0]
            mpos = _gp(mi_bridge)
            away = mid - mpos
            aw_len = float(np.linalg.norm(away))
            if aw_len > 1e-8:
                away /= aw_len
            else:
                away = np.array([0.0, 1.0, 0.0])

            # Distance between the two ring attachment points
            if len(ring_pos) >= 2:
                d_ab = float(np.linalg.norm(ring_pos[0] - ring_pos[1]))
                half_d = d_ab / 2.0
                h_sq = bridge_bond ** 2 - half_d ** 2
                h = float(np.sqrt(max(h_sq, 0.01)))
            else:
                h = bridge_bond * 0.7

            bridge_pos = mid + h * away

            # Ensure bridge atom is at least 2.5A from metal
            bm_vec = bridge_pos - mpos
            bm_dist = float(np.linalg.norm(bm_vec))
            min_bridge_metal = 2.5
            if bm_dist < min_bridge_metal and bm_dist > 1e-8:
                bridge_pos = mpos + bm_vec / bm_dist * min_bridge_metal
            elif bm_dist < 1e-8:
                bridge_pos = mpos + away * min_bridge_metal

            _sp(bridge_idx, bridge_pos)
            changed = True

            # Place non-hapto substituents on the bridge atom
            sub_count = 0
            ca_cb = ring_pos[1] - ring_pos[0] if len(ring_pos) >= 2 else np.array([1, 0, 0])
            ca_cb_len = float(np.linalg.norm(ca_cb))
            if ca_cb_len > 1e-8:
                x_ax = ca_cb / ca_cb_len
            else:
                x_ax = np.array([1.0, 0.0, 0.0])
            z_ax = away
            y_ax = np.cross(z_ax, x_ax)
            y_norm = float(np.linalg.norm(y_ax))
            if y_norm > 1e-8:
                y_ax /= y_norm
            else:
                y_ax = np.array([0.0, 0.0, 1.0])

            for nbr in ba.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in all_hapto:
                    continue
                nsym = nbr.GetSymbol()
                if nsym == 'C':
                    sub_bond = 1.87
                elif nsym == 'H':
                    sub_bond = 1.48
                else:
                    sub_bond = 1.50
                sub_angle = (sub_count - 0.5) * 2.1
                direction = (np.cos(sub_angle) * y_ax
                             + np.sin(sub_angle) * z_ax)
                d_n = float(np.linalg.norm(direction))
                if d_n > 1e-8:
                    direction /= d_n
                _sp(ni, bridge_pos + sub_bond * direction)
                sub_count += 1

    # =====================================================================
    # STEP 6: Fix H orientations on ring carbons (radially outward + tilt)
    # =====================================================================
    for gi, (mi, catoms) in enumerate(hapto_groups):
        atoms = sorted(set(catoms))
        mpos = _gp(mi)
        centroid = np.mean([_gp(a) for a in atoms], axis=0)
        mc_vec = centroid - mpos
        mc_dist = float(np.linalg.norm(mc_vec))
        if mc_dist < 1e-8:
            continue
        normal = mc_vec / mc_dist

        for ai in atoms:
            ring_pos = _gp(ai)
            outward = ring_pos - centroid
            out_len = float(np.linalg.norm(outward))
            if out_len > 1e-8:
                outward_unit = outward / out_len
            else:
                continue

            a = mol.GetAtomWithIdx(ai)
            h_nbrs = []
            for nbr in a.GetNeighbors():
                ni = nbr.GetIdx()
                if nbr.GetSymbol() == 'H' and ni not in all_hapto:
                    h_nbrs.append(ni)

            if not h_nbrs:
                continue

            # Place H atoms radially outward, tilted away from metal
            direction = outward_unit + 0.3 * normal
            d_n = float(np.linalg.norm(direction))
            if d_n > 1e-8:
                direction /= d_n
            for hi in h_nbrs:
                _sp(hi, ring_pos + 1.08 * direction)
                changed = True

    # =====================================================================
    # STEP 7: Clash resolution — push apart any non-bonded atoms that
    # are too close.  Hapto ring atoms are frozen; only non-ring atoms
    # are pushed outward from the metal or from each other.
    # Also enforce minimum distance from metals for non-bonded atoms.
    # =====================================================================
    if changed:
        # Determine minimum allowed distances
        frozen = set()
        for _mi, catoms in hapto_groups:
            frozen.update(catoms)
            frozen.add(_mi)

        # Build set of atom pairs that are directly bonded (should not be
        # pushed apart — bond length is the embedding's responsibility)
        bonded_pairs: set = set()
        for bond in mol.GetBonds():
            ai, aj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bonded_pairs.add((min(ai, aj), max(ai, aj)))

        # Build "rigid group" for each non-frozen heavy atom: the atom
        # itself plus all its directly bonded H atoms that are also not
        # frozen.  When a heavy atom is pushed, its H group moves with it.
        rigid_group: Dict[int, List[int]] = {}
        h_parent: Dict[int, int] = {}  # H atom → parent heavy atom
        for ai in range(mol.GetNumAtoms()):
            if ai in frozen:
                continue
            a = mol.GetAtomWithIdx(ai)
            if a.GetSymbol() == 'H':
                continue
            hs = []
            for nbr in a.GetNeighbors():
                ni = nbr.GetIdx()
                if nbr.GetSymbol() == 'H' and ni not in frozen:
                    hs.append(ni)
                    h_parent[ni] = ai
            rigid_group[ai] = hs

        def _push_atom(idx, displacement):
            """Push atom + its bonded H group by displacement vector.
            If idx is an H atom, redirect to its parent heavy atom
            (unless the parent is frozen)."""
            if idx in h_parent:
                parent = h_parent[idx]
                if parent not in frozen:
                    idx = parent
                # else: push H directly (parent frozen, can't move it)
            if idx in frozen:
                return  # safety: never move frozen atoms
            _sp(idx, _gp(idx) + displacement)
            if idx in rigid_group:
                for hi in rigid_group[idx]:
                    _sp(hi, _gp(hi) + displacement)

        n_atoms = mol.GetNumAtoms()
        for _pass in range(40):
            any_push = False
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    i_frozen = i in frozen
                    j_frozen = j in frozen
                    if i_frozen and j_frozen:
                        continue
                    if (i, j) in bonded_pairs:
                        continue

                    pi = _gp(i)
                    pj = _gp(j)
                    d = float(np.linalg.norm(pi - pj))

                    si = mol.GetAtomWithIdx(i).GetSymbol()
                    sj = mol.GetAtomWithIdx(j).GetSymbol()

                    i_is_metal = i in metals
                    j_is_metal = j in metals
                    if i_is_metal or j_is_metal:
                        other_sym = sj if i_is_metal else si
                        min_d = 2.5 if other_sym == 'H' else 2.2
                    elif si == 'H' and sj == 'H':
                        min_d = 1.5
                    elif si == 'H' or sj == 'H':
                        min_d = 1.0
                    else:
                        min_d = 1.2

                    if d >= min_d:
                        continue

                    if d < 1e-8:
                        direction = np.random.default_rng(
                            42 + _pass).standard_normal(3)
                        direction /= np.linalg.norm(direction)
                    else:
                        direction = (pj - pi) / d

                    gap = min_d - d
                    if i_frozen:
                        _push_atom(j, gap * direction)
                    elif j_frozen:
                        _push_atom(i, -gap * direction)
                    else:
                        _push_atom(i, -0.5 * gap * direction)
                        _push_atom(j, 0.5 * gap * direction)
                    any_push = True

            if not any_push:
                break

        # Final bond repair: fix any bonds distorted by the clash pushes.
        # For each bond, if the distance deviates >20% from the expected
        # bond length, move the lighter/non-frozen atom to correct it.
        for bond in mol.GetBonds():
            ai = bond.GetBeginAtomIdx()
            aj = bond.GetEndAtomIdx()
            si = mol.GetAtomWithIdx(ai).GetSymbol()
            sj = mol.GetAtomWithIdx(aj).GetSymbol()
            if si in _METAL_SET or sj in _METAL_SET:
                continue  # skip metal bonds
            pi = _gp(ai)
            pj = _gp(aj)
            d = float(np.linalg.norm(pi - pj))
            # Expected bond lengths
            if si == 'H' or sj == 'H':
                expected = 1.08
            elif si == 'Si' or sj == 'Si':
                expected = 1.87
            else:
                expected = 1.50
            if d < 1e-8 or abs(d - expected) / expected < 0.20:
                continue
            # Move the non-frozen (or lighter) atom
            direction = (pj - pi) / d
            correction = expected - d
            ai_frozen = ai in frozen
            aj_frozen = aj in frozen
            if ai_frozen and not aj_frozen:
                _sp(aj, pj + correction * direction)
            elif aj_frozen and not ai_frozen:
                _sp(ai, pi - correction * direction)
            elif not ai_frozen and not aj_frozen:
                # Move H toward its parent, or split equally
                if si == 'H':
                    _sp(ai, pi - correction * direction)
                elif sj == 'H':
                    _sp(aj, pj + correction * direction)
                else:
                    _sp(ai, pi - 0.5 * correction * direction)
                    _sp(aj, pj + 0.5 * correction * direction)

    return changed


# ---------------------------------------------------------------------------
# BFS propagation for non-hapto atoms after Phase 1 correction
# ---------------------------------------------------------------------------

def _propagate_non_hapto_atoms(
    mol,
    conf_id,
    hapto_groups,
    extra_fixed_indices: Optional[set] = None,
):
    """Re-place non-hapto atoms using BFS from the fixed hapto scaffold.

    After Phase 1 corrects hapto ring geometry, atoms NOT reachable from
    hapto rings (without crossing metals) may still be in bad ETKDG positions.
    These are typically sigma ligands and their substituents.  This function
    identifies them and re-places them using VSEPR-based local geometry,
    followed by bond-length relaxation and clash resolution.
    """
    from collections import deque
    import numpy as np

    conf = mol.GetConformer(conf_id)
    n_atoms = mol.GetNumAtoms()

    def _gp(i):
        p = conf.GetAtomPosition(i)
        return np.array([p.x, p.y, p.z])

    def _sp(i, pos):
        conf.SetAtomPosition(i, Point3D(
            float(pos[0]), float(pos[1]), float(pos[2])))

    # Identify metals and hapto ring atoms
    metals = set()
    hapto_ring_atoms = set()
    for mi, catoms in hapto_groups:
        metals.add(mi)
        for a in catoms:
            hapto_ring_atoms.add(a)

    # Compute hapto-reachable set: BFS from ring atoms, not crossing metals.
    # These atoms were already corrected by branch-aware rotation in Phase 1.
    hapto_reachable = set(hapto_ring_atoms)
    bfs_q = deque(list(hapto_ring_atoms))
    while bfs_q:
        cur = bfs_q.popleft()
        for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in hapto_reachable and ni not in metals:
                hapto_reachable.add(ni)
                bfs_q.append(ni)

    # Fixed = metals + hapto_reachable (already in correct positions)
    fixed = metals | hapto_reachable
    if extra_fixed_indices:
        fixed = fixed | set(extra_fixed_indices)

    # Find atoms that need re-placement
    needs_placement = set()
    for ai in range(n_atoms):
        if ai not in fixed:
            needs_placement.add(ai)

    if not needs_placement:
        return  # all atoms already in good positions

    # Bond length helper
    def _bl(sym_a, sym_b):
        if sym_a == 'H' or sym_b == 'H':
            other = sym_b if sym_a == 'H' else sym_a
            return {'C': 1.08, 'N': 1.01, 'O': 0.96, 'Si': 1.48,
                    'S': 1.34, 'P': 1.42}.get(other, 1.08)
        pair = frozenset([sym_a, sym_b])
        _bl_map = {
            frozenset(['C', 'C']): 1.50, frozenset(['C', 'N']): 1.47,
            frozenset(['C', 'O']): 1.43, frozenset(['C', 'Si']): 1.87,
            frozenset(['C', 'S']): 1.82, frozenset(['C', 'P']): 1.84,
            frozenset(['N', 'N']): 1.45, frozenset(['N', 'O']): 1.40,
            frozenset(['Si', 'Si']): 2.34,
        }
        if pair in _bl_map:
            return _bl_map[pair]
        r1 = _COVALENT_RADII.get(sym_a, 0.76)
        r2 = _COVALENT_RADII.get(sym_b, 0.76)
        return r1 + r2

    # BFS from fixed atoms to place non-hapto atoms
    placed = set(fixed)
    queue = deque()
    in_queue = set()

    for ai in fixed:
        for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in placed and ni not in in_queue:
                queue.append((ni, ai))
                in_queue.add(ni)

    while queue:
        ai, parent_idx = queue.popleft()
        if ai in placed:
            continue

        parent_pos = _gp(parent_idx)
        child_sym = mol.GetAtomWithIdx(ai).GetSymbol()
        parent_sym = mol.GetAtomWithIdx(parent_idx).GetSymbol()

        if parent_sym in _METAL_SET:
            bond_len = float(_get_ml_bond_length(parent_sym, child_sym))
        elif child_sym in _METAL_SET:
            bond_len = float(_get_ml_bond_length(child_sym, parent_sym))
        else:
            bond_len = _bl(parent_sym, child_sym)

        # VSEPR: direction away from already-placed neighbors of parent
        parent_atom = mol.GetAtomWithIdx(parent_idx)
        used_dirs = []
        for pnbr in parent_atom.GetNeighbors():
            pni = pnbr.GetIdx()
            if pni in placed and pni != ai:
                d = _gp(pni) - parent_pos
                d_len = float(np.linalg.norm(d))
                if d_len > 1e-8:
                    used_dirs.append(d / d_len)

        if len(used_dirs) == 0:
            direction = np.array([0.0, 0.0, 1.0])
        elif len(used_dirs) == 1:
            d0 = used_dirs[0]
            ref = (np.array([1.0, 0.0, 0.0]) if abs(d0[0]) < 0.9
                   else np.array([0.0, 1.0, 0.0]))
            perp = np.cross(d0, ref)
            perp = perp / max(float(np.linalg.norm(perp)), 1e-12)
            # ~109.5° from existing bond (tetrahedral)
            direction = (-d0 * np.cos(np.radians(70.5))
                         + perp * np.sin(np.radians(70.5)))
        elif len(used_dirs) == 2:
            avg = (used_dirs[0] + used_dirs[1]) / 2.0
            avg_len = float(np.linalg.norm(avg))
            if avg_len > 1e-8:
                direction = -avg / avg_len
            else:
                perp = np.cross(used_dirs[0], used_dirs[1])
                p_len = float(np.linalg.norm(perp))
                direction = (perp / p_len if p_len > 1e-8
                             else np.array([0.0, 0.0, 1.0]))
        else:
            avg = sum(used_dirs) / len(used_dirs)
            avg_len = float(np.linalg.norm(avg))
            direction = (-avg / avg_len if avg_len > 1e-8
                         else np.array([0.0, 0.0, 1.0]))

        d_norm = float(np.linalg.norm(direction))
        if d_norm > 1e-8:
            direction = direction / d_norm

        _sp(ai, parent_pos + bond_len * direction)
        placed.add(ai)

        for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in placed and ni not in in_queue:
                queue.append((ni, ai))
                in_queue.add(ni)

    # Bond-length relaxation for non-metal bonds
    bond_targets = []
    for bond in mol.GetBonds():
        bi = bond.GetBeginAtomIdx()
        bj = bond.GetEndAtomIdx()
        si = mol.GetAtomWithIdx(bi).GetSymbol()
        sj = mol.GetAtomWithIdx(bj).GetSymbol()
        if si in _METAL_SET or sj in _METAL_SET:
            continue
        # Only relax bonds involving at least one re-placed atom
        if bi not in needs_placement and bj not in needs_placement:
            continue
        bond_targets.append((bi, bj, _bl(si, sj)))

    rng = np.random.default_rng(42)
    for _relax_pass in range(60):
        max_err = 0.0
        forces = {}
        for bi, bj, target in bond_targets:
            diff = _gp(bj) - _gp(bi)
            d = float(np.linalg.norm(diff))
            err = abs(d - target)
            if err < 0.01:
                continue
            if d < 1e-8:
                diff = rng.standard_normal(3)
                d = float(np.linalg.norm(diff))
            unit = diff / d
            force = 0.3 * (d - target) * unit
            forces.setdefault(bi, np.zeros(3))
            forces.setdefault(bj, np.zeros(3))
            forces[bi] = forces[bi] + force
            forces[bj] = forces[bj] - force
            if err > max_err:
                max_err = err
        if max_err < 0.05:
            break
        for ai, f in forces.items():
            if ai in metals:
                continue
            w = 0.1 if ai in hapto_ring_atoms else 1.0
            _sp(ai, _gp(ai) + f * w)

    # Clash resolution for re-placed atoms
    bonded_pairs = set()
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonded_pairs.add((bi, bj))
        bonded_pairs.add((bj, bi))

    for _pass in range(20):
        any_push = False
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if (i, j) in bonded_pairs:
                    continue
                # At least one must be a re-placed atom
                if i not in needs_placement and j not in needs_placement:
                    continue
                pi = _gp(i)
                pj = _gp(j)
                d = float(np.linalg.norm(pi - pj))

                si = mol.GetAtomWithIdx(i).GetSymbol()
                sj = mol.GetAtomWithIdx(j).GetSymbol()

                if si in _METAL_SET or sj in _METAL_SET:
                    min_d = 2.0
                elif si == 'H' and sj == 'H':
                    min_d = 1.5
                elif si == 'H' or sj == 'H':
                    min_d = 1.0
                else:
                    min_d = 1.2

                if d >= min_d:
                    continue
                if d < 1e-8:
                    direction = rng.standard_normal(3)
                    direction /= max(float(np.linalg.norm(direction)), 1e-12)
                else:
                    direction = (pj - pi) / d

                gap = min_d - d
                i_fixed = i in fixed
                j_fixed = j in fixed
                if i_fixed:
                    _sp(j, pj + gap * direction)
                elif j_fixed:
                    _sp(i, pi - gap * direction)
                else:
                    _sp(i, pi - 0.5 * gap * direction)
                    _sp(j, pj + 0.5 * gap * direction)
                any_push = True
        if not any_push:
            break


# ---------------------------------------------------------------------------
# Sphere-based constructive hapto scaffold builder
# ---------------------------------------------------------------------------

def _build_hapto_scaffold(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> bool:
    """Construct hapto complex geometry using sphere-based ligand placement.

    Instead of relying on RDKit's distance geometry (ETKDG), which has no
    knowledge of eta-coordination, this builds the geometry analytically:

    1. Place each metal at a fixed position
    2. Distribute hapto-group centroids on a sphere (r = M-centroid distance)
    3. Construct each hapto ring as a regular polygon perpendicular to M->centroid
    4. Place sigma-bonded ligands in remaining coordination directions
    5. Place substituents with correct bond lengths and angles
    6. BFS-propagate remaining atoms with local VSEPR geometry rules

    The mol gets a new conformer with the constructed coordinates.
    Returns True on success.
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    n_atoms = mol.GetNumAtoms()
    coords = np.full((n_atoms, 3), np.nan)
    placed: set = set()

    # ---- Classify atoms ----
    by_metal: Dict[int, List[List[int]]] = {}
    all_hapto_atoms: set = set()
    atom_to_hapto_group: Dict[int, int] = {}  # atom_idx -> group index
    for g_idx, (metal_idx, grp) in enumerate(hapto_groups):
        by_metal.setdefault(metal_idx, []).append(grp)
        all_hapto_atoms.update(grp)
        for a in grp:
            atom_to_hapto_group[a] = g_idx
    metal_indices = sorted(by_metal.keys())

    sigma_donors: Dict[int, List[int]] = {}
    for mi in metal_indices:
        donors = []
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in all_hapto_atoms and nbr.GetSymbol() not in _METAL_SET:
                donors.append(ni)
        sigma_donors[mi] = donors

    # Detect ansa bridges (multiple bridges per group pair supported)
    ansa_bridge_map: Dict[Tuple[int, int], List[int]] = {}
    all_bridge_atoms: set = set()
    try:
        bridges = _find_ansa_bridges(mol, hapto_groups)
        for bridge_atom, gi, gj in bridges:
            key = (min(gi, gj), max(gi, gj))
            ansa_bridge_map.setdefault(key, []).append(bridge_atom)
            all_bridge_atoms.add(bridge_atom)
    except Exception:
        pass

    # ---- Helper functions ----
    def _ortho_basis(normal):
        n = normal / max(float(np.linalg.norm(normal)), 1e-12)
        ref = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        u = np.cross(n, ref)
        u = u / max(float(np.linalg.norm(u)), 1e-12)
        v = np.cross(n, u)
        v = v / max(float(np.linalg.norm(v)), 1e-12)
        return u, v

    def _sphere_dirs(n):
        if n <= 0:
            return []
        if n == 1:
            return [np.array([0.0, 0.0, 1.0])]
        if n == 2:
            return [np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0])]
        dirs = []
        golden = (1.0 + np.sqrt(5.0)) / 2.0
        for k in range(n):
            theta = np.arccos(1.0 - 2.0 * (k + 0.5) / n)
            phi = 2.0 * np.pi * k / golden
            dirs.append(np.array([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(theta),
            ]))
        return dirs

    def _ring_traversal(grp):
        atom_set = set(grp)
        idx_map = {a: i for i, a in enumerate(grp)}
        eta = len(grp)
        adj = [[] for _ in range(eta)]
        for i_loc, atom_idx in enumerate(grp):
            for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set and ni != atom_idx:
                    j_loc = idx_map[ni]
                    if j_loc not in adj[i_loc]:
                        adj[i_loc].append(j_loc)
        is_cyclic = all(len(adj[i]) >= 2 for i in range(eta))
        if is_cyclic and eta >= 3:
            order = []
            visited = [False] * eta
            cur, prev = 0, -1
            for _ in range(eta):
                order.append(cur)
                visited[cur] = True
                nxt = -1
                for nb in adj[cur]:
                    if nb != prev and not visited[nb]:
                        nxt = nb
                        break
                if nxt == -1:
                    break
                prev, cur = cur, nxt
            if len(order) == eta:
                return [grp[i] for i in order]
        return grp

    def _build_chain_order(grp, mol_ref):
        """Order a non-cyclic hapto group as a chain following molecular topology.

        For disconnected hapto groups (e.g., two eta2 pairs bridged by
        non-hapto atoms), determines ordering by shortest path through
        the full molecular graph.
        """
        grp_set = set(grp)
        adj_intra: Dict[int, List[int]] = {a: [] for a in grp}
        for a in grp:
            for nb in mol_ref.GetAtomWithIdx(a).GetNeighbors():
                if nb.GetIdx() in grp_set:
                    adj_intra[a].append(nb.GetIdx())

        # Find connected components within the group
        comp_id: Dict[int, int] = {}
        components: List[List[int]] = []
        for a in grp:
            if a in comp_id:
                continue
            comp: List[int] = []
            stack = [a]
            while stack:
                cur = stack.pop()
                if cur in comp_id:
                    continue
                comp_id[cur] = len(components)
                comp.append(cur)
                for nb in adj_intra[cur]:
                    if nb not in comp_id:
                        stack.append(nb)
            components.append(comp)

        if len(components) == 1:
            # Single connected component: simple chain walk
            start = grp[0]
            for a in grp:
                if len(adj_intra[a]) <= 1:
                    start = a
                    break
            chain = [start]
            visited = {start}
            while len(chain) < len(grp):
                cur = chain[-1]
                nxt = None
                for nb in adj_intra[cur]:
                    if nb not in visited:
                        nxt = nb
                        break
                if nxt is None:
                    break
                chain.append(nxt)
                visited.add(nxt)
            return chain

        # Multiple disconnected components: order them by BFS shortest
        # path through non-hapto atoms in the molecular graph
        def _bfs_dist(src, tgt_set):
            """BFS from src to any atom in tgt_set, through non-metal atoms."""
            from collections import deque
            q = deque([(src, 0)])
            seen = {src}
            while q:
                cur, d = q.popleft()
                if cur in tgt_set:
                    return d
                for nb in mol_ref.GetAtomWithIdx(cur).GetNeighbors():
                    ni = nb.GetIdx()
                    if ni not in seen and nb.GetSymbol() not in _METAL_SET:
                        seen.add(ni)
                        q.append((ni, d + 1))
            return 999

        # Order components by shortest path between consecutive pairs
        ordered_comps = [components[0]]
        remaining = list(range(1, len(components)))
        while remaining:
            last_comp = ordered_comps[-1]
            best_ci = remaining[0]
            best_d = 999
            for ci in remaining:
                # Find min BFS distance from any atom in last_comp to comp[ci]
                comp_set = set(components[ci])
                for a in last_comp:
                    d = _bfs_dist(a, comp_set)
                    if d < best_d:
                        best_d = d
                        best_ci = ci
            ordered_comps.append(components[best_ci])
            remaining.remove(best_ci)

        # Within each component, order as chain from endpoint
        chain = []
        for comp in ordered_comps:
            if len(comp) == 1:
                chain.append(comp[0])
                continue
            start = comp[0]
            for a in comp:
                if len(adj_intra[a]) <= 1:
                    start = a
                    break
            sub = [start]
            visited = {start}
            while len(sub) < len(comp):
                cur = sub[-1]
                nxt = None
                for nb in adj_intra[cur]:
                    if nb not in visited:
                        nxt = nb
                        break
                if nxt is None:
                    break
                sub.append(nxt)
                visited.add(nxt)
            # Orient: if previous chain end is closer to sub[-1], reverse
            if chain:
                prev_end = chain[-1]
                d_fwd = _bfs_dist(prev_end, {sub[0]})
                d_rev = _bfs_dist(prev_end, {sub[-1]})
                if d_rev < d_fwd:
                    sub = sub[::-1]
            chain.extend(sub)
        return chain

    def _bond_len(sym1, sym2):
        pair = frozenset([sym1, sym2])
        if 'H' in pair:
            other = sym2 if sym1 == 'H' else sym1
            return {'C': 1.08, 'N': 1.01, 'O': 0.96, 'Si': 1.48, 'B': 1.19}.get(other, 1.08)
        if pair == frozenset(['Si', 'C']):
            return 1.87
        if pair == frozenset(['C', 'C']):
            return 1.50
        if pair == frozenset(['C', 'N']):
            return 1.47
        if pair == frozenset(['C', 'O']):
            return 1.43
        r1 = _COVALENT_RADII.get(sym1, 0.76)
        r2 = _COVALENT_RADII.get(sym2, 0.76)
        return r1 + r2

    def _rodrigues(points, center, axis, angle):
        k = axis / max(float(np.linalg.norm(axis)), 1e-12)
        c, s = np.cos(angle), np.sin(angle)
        out = []
        for p in points:
            v = p - center
            v_rot = v * c + np.cross(k, v) * s + k * float(np.dot(k, v)) * (1.0 - c)
            out.append(center + v_rot)
        return out

    # ---- Place metals ----
    if len(metal_indices) == 1:
        coords[metal_indices[0]] = [0.0, 0.0, 0.0]
        placed.add(metal_indices[0])
    else:
        for i, mi in enumerate(metal_indices):
            if i == 0:
                coords[mi] = [0.0, 0.0, 0.0]
            else:
                prev_mi = metal_indices[i - 1]
                si = mol.GetAtomWithIdx(mi).GetSymbol()
                sj = mol.GetAtomWithIdx(prev_mi).GetSymbol()
                bond = mol.GetBondBetweenAtoms(mi, prev_mi)
                d = (_COVALENT_RADII.get(si, 1.5) + _COVALENT_RADII.get(sj, 1.5) + 0.1
                     if bond is not None else 3.5)
                coords[mi] = coords[prev_mi] + np.array([d, 0.0, 0.0])
            placed.add(mi)

    # ---- Build global group list ----
    group_list: List[Tuple[int, List[int]]] = []
    for mi in metal_indices:
        for grp in by_metal[mi]:
            group_list.append((mi, grp))

    # ---- For each metal: build coordination sphere ----
    for mi in metal_indices:
        m_pos = coords[mi].copy()
        m_sym = mol.GetAtomWithIdx(mi).GetSymbol()
        groups = by_metal[mi]
        n_hapto = len(groups)
        n_sigma = len(sigma_donors.get(mi, []))
        n_total = n_hapto + n_sigma
        if n_total == 0:
            continue

        # Map local group indices -> global
        local_to_global = []
        for g_global, (gm, _) in enumerate(group_list):
            if gm == mi:
                local_to_global.append(g_global)

        # Detect ansa constraints for this metal
        ansa_angle = None
        ansa_local_i = None
        ansa_local_j = None
        ansa_bridge_indices: List[int] = []
        for (gi_g, gj_g), br_list in ansa_bridge_map.items():
            try:
                li = local_to_global.index(gi_g)
                lj = local_to_global.index(gj_g)
            except ValueError:
                continue
            ansa_local_i, ansa_local_j = li, lj
            ansa_bridge_indices = list(br_list)
            g_i, g_j = groups[li], groups[lj]
            eta_i, eta_j = len(g_i), len(g_j)
            r_i = _target_mc_dist(m_sym, eta_i)
            r_j = _target_mc_dist(m_sym, eta_j)
            R_i = 1.40 / (2.0 * np.sin(np.pi / max(eta_i, 3)))
            R_j = 1.40 / (2.0 * np.sin(np.pi / max(eta_j, 3)))
            r_avg = (r_i + r_j) / 2.0
            R_avg = (R_i + R_j) / 2.0
            # Use first bridge for angle computation
            br_sym = mol.GetAtomWithIdx(br_list[0]).GetSymbol()
            br_bond = {'Si': 1.87, 'C': 1.54, 'Ge': 1.94}.get(br_sym, 1.80)
            br_angle_deg = {'Si': 93.0, 'C': 109.5, 'Ge': 95.0}.get(br_sym, 100.0)
            d_target = 2.0 * br_bond * np.sin(np.radians(br_angle_deg) / 2.0)
            hyp = np.sqrt(r_avg ** 2 + R_avg ** 2)
            ratio = d_target / (2.0 * hyp)
            if abs(ratio) <= 1.0:
                psi = np.arctan2(R_avg, r_avg) + np.arcsin(ratio)
                ansa_angle = 2.0 * psi
            else:
                ansa_angle = np.radians(140.0)
            break

        # Generate direction vectors
        hapto_dirs: List[Optional[np.ndarray]] = [None] * n_hapto

        if ansa_angle is not None and ansa_local_i is not None:
            half = ansa_angle / 2.0
            d1 = np.array([0.0, np.sin(half), np.cos(half)])
            d2 = np.array([0.0, -np.sin(half), np.cos(half)])
            hapto_dirs[ansa_local_i] = d1 / float(np.linalg.norm(d1))
            hapto_dirs[ansa_local_j] = d2 / float(np.linalg.norm(d2))

        # Linked groups: hapto pairs connected through short non-hapto paths
        # (e.g., COD eta2+eta2 bridged by CH2-CH2).
        # Detect and mark for combined rectangular placement.
        linked_pairs: List[Tuple[int, int]] = []  # (local_i, local_j)
        if n_hapto >= 2:
            from collections import deque
            linked_used: set = set()
            for li in range(n_hapto):
                if li in linked_used or hapto_dirs[li] is not None:
                    continue
                for lj in range(li + 1, n_hapto):
                    if lj in linked_used or hapto_dirs[lj] is not None:
                        continue
                    grp_i_set = set(groups[li])
                    grp_j_set = set(groups[lj])
                    linked = False
                    for start in groups[li]:
                        q = deque([(start, 0)])
                        seen = {start}
                        while q:
                            cur, d_bfs = q.popleft()
                            if cur in grp_j_set:
                                linked = True
                                break
                            if d_bfs >= 4:
                                continue
                            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                                ni = nb.GetIdx()
                                if ni not in seen and nb.GetSymbol() not in _METAL_SET:
                                    seen.add(ni)
                                    q.append((ni, d_bfs + 1))
                        if linked:
                            break
                    if linked:
                        linked_pairs.append((li, lj))
                        linked_used.add(li)
                        linked_used.add(lj)
                        # Assign same direction (will be placed as rectangle)
                        hapto_dirs[li] = np.array([0.0, 0.0, 1.0])
                        hapto_dirs[lj] = np.array([0.0, 0.0, 1.0])
                        break

        # ---- Assign ideal directions for unassigned hapto groups ----
        unassigned_hapto = [i for i in range(n_hapto) if hapto_dirs[i] is None]
        assigned_dirs = [hapto_dirs[i] for i in range(n_hapto)
                         if hapto_dirs[i] is not None]

        # In multi-metal systems, compute anti-M2 direction so hapto
        # groups point away from the other metal(s).
        _anti_other_metal = None
        if len(metal_indices) > 1:
            other_positions = [coords[mj] for mj in metal_indices if mj != mi
                               and mj in placed]
            if other_positions:
                avg_other = np.mean(other_positions, axis=0)
                away = m_pos - avg_other
                away_len = float(np.linalg.norm(away))
                if away_len > 1e-8:
                    _anti_other_metal = away / away_len

        if unassigned_hapto:
            if not assigned_dirs:
                # All hapto groups unassigned: use analytical placement
                if len(unassigned_hapto) == 1:
                    hapto_dirs[unassigned_hapto[0]] = (_anti_other_metal
                        if _anti_other_metal is not None
                        else np.array([0.0, 0.0, 1.0]))
                elif len(unassigned_hapto) == 2:
                    # Sandwich: ~175° apart (nearly anti-parallel).
                    # In multi-metal: orient perpendicular to M-M axis.
                    if _anti_other_metal is not None:
                        u_s, _ = _ortho_basis(_anti_other_metal)
                        hapto_dirs[unassigned_hapto[0]] = u_s
                        a175 = np.radians(175.0)
                        hapto_dirs[unassigned_hapto[1]] = (
                            np.cos(a175) * u_s + np.sin(a175)
                            * _anti_other_metal)
                    else:
                        hapto_dirs[unassigned_hapto[0]] = np.array([0.0, 0.0, 1.0])
                        a175 = np.radians(175.0)
                        hapto_dirs[unassigned_hapto[1]] = np.array(
                            [0.0, np.sin(a175), np.cos(a175)])
                else:
                    # Evenly distributed in a plane (120° for 3, 90° for 4, …)
                    for k, ui in enumerate(unassigned_hapto):
                        a = 2.0 * np.pi * k / len(unassigned_hapto)
                        hapto_dirs[ui] = np.array(
                            [np.cos(a), np.sin(a), 0.0])
            else:
                # Some already assigned (ansa/linked): place rest maximally
                # separated from all existing directions.
                used = list(assigned_dirs)
                candidates = _sphere_dirs(60)
                for ui in unassigned_hapto:
                    best_dir = None
                    best_max_dot = 1.0
                    for cd in candidates:
                        max_dot = max(float(np.dot(cd, ud)) for ud in used)
                        if max_dot < best_max_dot:
                            best_max_dot = max_dot
                            best_dir = cd
                    hapto_dirs[ui] = (best_dir if best_dir is not None
                                      else np.array([0.0, 0.0, 1.0]))
                    used.append(hapto_dirs[ui])

        # Fill sigma donor directions, maximally separated from hapto dirs.
        # Piano-stool special case: 1 hapto + n sigma → place sigma donors
        # analytically on a cone at ~125° from the hapto axis.
        sigma_dir_list = []
        if n_sigma > 0:
            all_used = [hd for hd in hapto_dirs if hd is not None]
            if n_hapto == 1 and len(all_used) == 1:
                hapto_ax = all_used[0] / max(float(np.linalg.norm(all_used[0])), 1e-12)
                cone_angle = np.radians(125.0)
                u_cone, v_cone = _ortho_basis(hapto_ax)
                for k in range(n_sigma):
                    phi = 2.0 * np.pi * k / n_sigma
                    d = (np.cos(cone_angle) * hapto_ax
                         + np.sin(cone_angle) * (np.cos(phi) * u_cone + np.sin(phi) * v_cone))
                    sigma_dir_list.append(d / max(float(np.linalg.norm(d)), 1e-12))
            else:
                candidates = _sphere_dirs(max(n_sigma + len(all_used) + 8, 60))
                candidates.sort(key=lambda cd: max(
                    float(np.dot(cd, ud)) for ud in all_used) if all_used else 0.0)
                for cd in candidates:
                    if len(sigma_dir_list) >= n_sigma:
                        break
                    ok = True
                    for ud in all_used + sigma_dir_list:
                        if float(np.dot(cd, ud)) > 0.50:
                            ok = False
                            break
                    if ok:
                        sigma_dir_list.append(cd)
            while len(sigma_dir_list) < n_sigma:
                sigma_dir_list.append(np.array([1.0, 0.0, 0.0]))

        for i in range(n_hapto):
            if hapto_dirs[i] is None:
                hapto_dirs[i] = np.array([0.0, 0.0, 1.0])

        # ---- Place linked pairs as rectangles first ----
        linked_placed_locals: set = set()
        for li, lj in linked_pairs:
            direction = hapto_dirs[li]
            grp_i, grp_j = groups[li], groups[lj]
            mc_i = _target_mc_dist(m_sym, len(grp_i))
            mc_j = _target_mc_dist(m_sym, len(grp_j))
            mc_avg = (mc_i + mc_j) / 2.0
            combined_centroid = m_pos + direction * mc_avg
            u_lnk, v_lnk = _ortho_basis(direction)

            # Find which atoms bridge the two groups, to orient the rectangle
            # Chain order: build chain for combined atoms through mol graph
            all_linked = list(grp_i) + list(grp_j)
            chain_i = _build_chain_order(grp_i, mol)
            chain_j = _build_chain_order(grp_j, mol)

            # Orient chains: find which ends face each other (BFS shortest)
            from collections import deque
            def _bfs_d(src, tgt):
                q = deque([(src, 0)])
                seen = {src}
                while q:
                    cur, dd = q.popleft()
                    if cur == tgt:
                        return dd
                    if dd >= 6:
                        return 99
                    for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                        ni = nb.GetIdx()
                        if ni not in seen and nb.GetSymbol() not in _METAL_SET:
                            seen.add(ni)
                            q.append((ni, dd + 1))
                return 99

            # Find the pair of endpoints with shortest BFS distance
            best_di, best_dj = 0, 0
            best_bfs = 99
            for di_end in [0, -1]:
                for dj_end in [0, -1]:
                    bd = _bfs_d(chain_i[di_end], chain_j[dj_end])
                    if bd < best_bfs:
                        best_bfs = bd
                        best_di, best_dj = di_end, dj_end
            # Orient chains so facing ends are the last atoms
            if best_di == 0:
                chain_i = chain_i[::-1]
            if best_dj == 0:
                chain_j = chain_j[::-1]

            # Place as rectangle: two parallel rows separated by bridge gap
            # Bridge gap ≈ n_bridge_atoms * 1.54A with 109.5° angles
            bridge_gap = max(best_bfs * 1.00, 2.2)  # empirical
            half_gap = bridge_gap / 2.0
            cc = 1.40  # C-C distance within each pair

            for k, atom_idx in enumerate(chain_i):
                offset_v = (k - (len(chain_i) - 1) / 2.0) * cc
                coords[atom_idx] = combined_centroid + half_gap * u_lnk + offset_v * v_lnk
                placed.add(atom_idx)
            for k, atom_idx in enumerate(chain_j):
                offset_v = (k - (len(chain_j) - 1) / 2.0) * cc
                coords[atom_idx] = combined_centroid - half_gap * u_lnk + offset_v * v_lnk
                placed.add(atom_idx)
            linked_placed_locals.add(li)
            linked_placed_locals.add(lj)

        # ---- Build hapto rings/chains on sphere ----
        for g_local, grp in enumerate(groups):
            if g_local in linked_placed_locals:
                continue  # already placed as linked rectangle
            direction = hapto_dirs[g_local]
            eta = len(grp)
            mc_dist = _target_mc_dist(m_sym, eta)
            centroid = m_pos + direction * mc_dist

            if eta >= 3:
                R_ring = 1.40 / (2.0 * np.sin(np.pi / eta))
            else:
                R_ring = 0.69

            u, v = _ortho_basis(direction)
            ordered = _ring_traversal(grp)

            # Check if group is cyclic (all atoms have >=2 intra-group bonds)
            grp_set = set(grp)
            is_cyclic = True
            for a_idx in grp:
                intra = sum(1 for nb in mol.GetAtomWithIdx(a_idx).GetNeighbors()
                            if nb.GetIdx() in grp_set)
                if intra < 2:
                    is_cyclic = False
                    break

            if eta == 2:
                coords[ordered[0]] = centroid + R_ring * u
                coords[ordered[1]] = centroid - R_ring * u
            elif is_cyclic:
                # Regular polygon (closed ring)
                for k, atom_idx in enumerate(ordered):
                    angle = 2.0 * np.pi * k / eta
                    coords[atom_idx] = centroid + R_ring * (
                        np.cos(angle) * u + np.sin(angle) * v
                    )
            else:
                # Non-cyclic group: check for disconnected components
                # (e.g., COD eta4 = two eta2 pairs bridged by CH2-CH2)
                chain = _build_chain_order(grp, mol)
                # Find connected components within the hapto group
                comp_list = []
                cur_comp = [chain[0]]
                for k in range(1, len(chain)):
                    # Check if chain[k] bonds to chain[k-1] within grp
                    bonded = False
                    for nb in mol.GetAtomWithIdx(chain[k]).GetNeighbors():
                        if nb.GetIdx() == chain[k - 1]:
                            bonded = True
                            break
                    if bonded:
                        cur_comp.append(chain[k])
                    else:
                        comp_list.append(cur_comp)
                        cur_comp = [chain[k]]
                comp_list.append(cur_comp)

                if len(comp_list) >= 2:
                    # Disconnected components: place as parallel pairs
                    # on a rectangle to keep both gaps bridgeable
                    n_comp = len(comp_list)
                    for ci, comp in enumerate(comp_list):
                        # Angle offset for this component on the circle
                        comp_center_angle = 2.0 * np.pi * ci / n_comp
                        cc_half = 1.40 * (len(comp) - 1) / 2.0
                        for k, atom_idx in enumerate(comp):
                            # Spread atoms within component along v axis
                            offset = (k - (len(comp) - 1) / 2.0) * 1.40
                            pos = centroid + (
                                R_ring * np.cos(comp_center_angle) * u
                                + offset * v
                            )
                            coords[atom_idx] = pos
                else:
                    # Single connected chain: place along arc
                    step_angle = 2.0 * np.arcsin(
                        min(1.40 / (2.0 * R_ring), 1.0))
                    total_span = step_angle * (eta - 1)
                    start_angle = -total_span / 2.0
                    for k, atom_idx in enumerate(chain):
                        angle = start_angle + k * step_angle
                        coords[atom_idx] = centroid + R_ring * (
                            np.cos(angle) * u + np.sin(angle) * v
                        )
            for atom_idx in grp:
                placed.add(atom_idx)

        # ---- Orient ansa-bridged rings so bridge carbons face each other ----
        if ansa_bridge_indices and ansa_local_i is not None:
            grp_i = groups[ansa_local_i]
            grp_j = groups[ansa_local_j]
            # Collect all neighbors of ALL bridge atoms
            br_nbrs: set = set()
            for _bi in ansa_bridge_indices:
                for _bn in mol.GetAtomWithIdx(_bi).GetNeighbors():
                    br_nbrs.add(_bn.GetIdx())
            for grp_src, dir_src, grp_tgt_local in [
                (grp_i, hapto_dirs[ansa_local_i], ansa_local_j),
                (grp_j, hapto_dirs[ansa_local_j], ansa_local_i),
            ]:
                c_src = next((a for a in grp_src if a in br_nbrs), None)
                if c_src is None:
                    continue
                cent_src = m_pos + dir_src * _target_mc_dist(m_sym, len(grp_src))
                dir_tgt = hapto_dirs[grp_tgt_local]
                cent_tgt = m_pos + dir_tgt * _target_mc_dist(
                    m_sym, len(groups[grp_tgt_local])
                )
                toward = cent_tgt - cent_src
                toward_proj = toward - float(np.dot(toward, dir_src)) * dir_src
                tp_len = float(np.linalg.norm(toward_proj))
                if tp_len < 1e-8:
                    continue
                toward_proj = toward_proj / tp_len
                vc = coords[c_src] - cent_src
                vc_proj = vc - float(np.dot(vc, dir_src)) * dir_src
                vp_len = float(np.linalg.norm(vc_proj))
                if vp_len < 1e-8:
                    continue
                vc_proj = vc_proj / vp_len
                cos_r = float(np.clip(np.dot(vc_proj, toward_proj), -1.0, 1.0))
                sin_r = float(np.dot(np.cross(vc_proj, toward_proj), dir_src))
                rot_angle = float(np.arctan2(sin_r, cos_r))
                if abs(rot_angle) > 0.01:
                    pts = [coords[a].copy() for a in grp_src]
                    pts_rot = _rodrigues(pts, cent_src, dir_src, rot_angle)
                    for a_idx, new_p in zip(grp_src, pts_rot):
                        coords[a_idx] = new_p

        # ---- Place sigma donors ----
        for s_idx, donor_idx in enumerate(sigma_donors.get(mi, [])):
            if donor_idx in placed:
                continue
            direction = (sigma_dir_list[s_idx] if s_idx < len(sigma_dir_list)
                         else np.array([1.0, 0.0, 0.0]))
            donor_sym = mol.GetAtomWithIdx(donor_idx).GetSymbol()
            bl = _get_ml_bond_length(m_sym, donor_sym)
            coords[donor_idx] = m_pos + direction * bl
            placed.add(donor_idx)

    # ---- Place substituents on hapto ring atoms ----
    for metal_idx, grp in hapto_groups:
        ring_set = set(grp)
        m_pos = coords[metal_idx]
        centroid = np.mean([coords[a] for a in grp], axis=0)
        mc_dir = centroid - m_pos
        mc_len = float(np.linalg.norm(mc_dir))
        normal = mc_dir / mc_len if mc_len > 1e-8 else np.array([0.0, 0.0, 1.0])

        for atom_idx in grp:
            ring_pos = coords[atom_idx]
            outward = ring_pos - centroid
            outward_len = float(np.linalg.norm(outward))
            outward_unit = (outward / outward_len if outward_len > 1e-8
                            else np.array([1.0, 0.0, 0.0]))
            subs = []
            for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in ring_set or ni == metal_idx or ni in placed:
                    continue
                if ni in all_bridge_atoms:
                    continue
                subs.append(nbr)

            for s_k, nbr in enumerate(subs):
                ni = nbr.GetIdx()
                bl = _bond_len(
                    mol.GetAtomWithIdx(atom_idx).GetSymbol(), nbr.GetSymbol()
                )
                if len(subs) == 1:
                    direction = outward_unit + 0.3 * normal
                else:
                    tilt = (s_k - (len(subs) - 1) / 2.0) * 1.2
                    direction = outward_unit * np.cos(tilt) + normal * np.sin(tilt)
                d_norm = float(np.linalg.norm(direction))
                if d_norm > 1e-8:
                    direction = direction / d_norm
                coords[ni] = ring_pos + bl * direction
                placed.add(ni)

    # ---- Place ALL ansa bridge atoms ----
    for (gi, gj), bridge_list in ansa_bridge_map.items():
        grp_i = hapto_groups[gi][1]
        grp_j = hapto_groups[gj][1]
        metal_idx = hapto_groups[gi][0]

        for bridge_idx in bridge_list:
            if bridge_idx in placed:
                continue
            bridge_atom = mol.GetAtomWithIdx(bridge_idx)
            bridge_sym = bridge_atom.GetSymbol()
            br_nbrs = set(n.GetIdx() for n in bridge_atom.GetNeighbors())

            c_i = next((a for a in grp_i if a in br_nbrs), None)
            c_j = next((a for a in grp_j if a in br_nbrs), None)
            if c_i is not None and c_j is not None:
                pos_a, pos_b = coords[c_i], coords[c_j]
                mid = (pos_a + pos_b) / 2.0
                bl = _bond_len(bridge_sym, 'C')
                d_ab = float(np.linalg.norm(pos_a - pos_b))
                h_sq = bl ** 2 - (d_ab / 2.0) ** 2
                h = float(np.sqrt(max(h_sq, 0.01)))
                away = mid - coords[metal_idx]
                a_norm = float(np.linalg.norm(away))
                away = away / a_norm if a_norm > 1e-8 else np.array([0.0, 1.0, 0.0])
                coords[bridge_idx] = mid + h * away
                placed.add(bridge_idx)

                # Place bridge substituents using robust local frame
                ca_cb = pos_b - pos_a
                x_ax = ca_cb / max(d_ab, 1e-8)
                z_ax = away
                y_ax = np.cross(z_ax, x_ax)
                y_norm = float(np.linalg.norm(y_ax))
                if y_norm < 1e-8:
                    # Degenerate: z_ax parallel to x_ax, pick perpendicular
                    ref = (np.array([1.0, 0.0, 0.0]) if abs(z_ax[0]) < 0.9
                           else np.array([0.0, 1.0, 0.0]))
                    y_ax = np.cross(z_ax, ref)
                    y_norm = float(np.linalg.norm(y_ax))
                y_ax = y_ax / max(y_norm, 1e-12)

                sub_count = 0
                for pnbr in bridge_atom.GetNeighbors():
                    pni = pnbr.GetIdx()
                    if pni in placed or pni in all_hapto_atoms:
                        continue
                    sub_bl = _bond_len(bridge_sym, pnbr.GetSymbol())
                    sub_angle = (sub_count - 0.5) * 2.1
                    direction = np.cos(sub_angle) * y_ax + np.sin(sub_angle) * z_ax
                    d_norm = float(np.linalg.norm(direction))
                    if d_norm > 1e-8:
                        direction = direction / d_norm
                    coords[pni] = coords[bridge_idx] + sub_bl * direction
                    placed.add(pni)
                    sub_count += 1
                sub_count += 1

    # ---- Pre-compute ring info for inline ring placement during BFS ----
    try:
        Chem.FastFindRings(mol)
        _ri = mol.GetRingInfo()
        _all_rings = list(_ri.AtomRings())
    except Exception:
        _all_rings = []
    _organic_rings = []
    for ring in _all_rings:
        ring_set = set(ring)
        if ring_set <= all_hapto_atoms:
            continue
        if any(mol.GetAtomWithIdx(a).GetSymbol() in _METAL_SET for a in ring):
            continue
        _organic_rings.append(ring)
    _atom_to_rings: Dict[int, List[int]] = {}
    for ri_idx, ring in enumerate(_organic_rings):
        for a in ring:
            _atom_to_rings.setdefault(a, []).append(ri_idx)
    _ring_done: set = set()

    def _place_ring_inline(ring_tuple, anchor_idx):
        ring_set = set(ring_tuple)
        rsize = len(ring_tuple)
        cc_ring = 1.40 if rsize <= 6 else 1.50
        radius = cc_ring / (2.0 * np.sin(np.pi / max(rsize, 3)))

        anchor_pos = coords[anchor_idx]
        parent_of_anchor = None
        for nbr in mol.GetAtomWithIdx(anchor_idx).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in placed and ni not in ring_set:
                parent_of_anchor = ni
                break
        if parent_of_anchor is not None:
            bond_dir = anchor_pos - coords[parent_of_anchor]
            bd_len = float(np.linalg.norm(bond_dir))
            if bd_len > 1e-8:
                bond_dir = bond_dir / bd_len
            else:
                bond_dir = np.array([1.0, 0.0, 0.0])
        else:
            bond_dir = np.array([1.0, 0.0, 0.0])

        u_r, v_r = _ortho_basis(bond_dir)
        ring_center = anchor_pos + bond_dir * radius

        ring_adj = {a: [] for a in ring_tuple}
        for a in ring_tuple:
            for nbr in mol.GetAtomWithIdx(a).GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    ring_adj[a].append(nbr.GetIdx())
        chain = [anchor_idx]
        visited_r = {anchor_idx}
        while len(chain) < rsize:
            cur = chain[-1]
            nxt = None
            for nb in ring_adj[cur]:
                if nb not in visited_r:
                    nxt = nb
                    break
            if nxt is None:
                break
            chain.append(nxt)
            visited_r.add(nxt)
        if len(chain) != rsize:
            return False

        anchor_offset = np.arctan2(
            float(np.dot(anchor_pos - ring_center, v_r)),
            float(np.dot(anchor_pos - ring_center, u_r)),
        )
        for k, atom_idx in enumerate(chain):
            angle = anchor_offset + 2.0 * np.pi * k / rsize
            coords[atom_idx] = ring_center + radius * (
                np.cos(angle) * u_r + np.sin(angle) * v_r)
            placed.add(atom_idx)
        return True

    # ---- BFS propagate remaining atoms with VSEPR local geometry ----
    for _iteration in range(n_atoms * 3):
        progress = False
        for atom in mol.GetAtoms():
            ai = atom.GetIdx()
            if ai in placed:
                continue
            parent = None
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in placed:
                    parent = nbr
                    break
            if parent is None:
                continue

            # If this atom belongs to an unplaced organic ring, place
            # the entire ring analytically instead of via BFS.
            ring_handled = False
            for _ri_idx in _atom_to_rings.get(ai, []):
                if _ri_idx in _ring_done:
                    continue
                _ring_t = _organic_rings[_ri_idx]
                _ring_s = set(_ring_t)
                _ring_anchors = [a for a in _ring_t if a in placed]
                if _ring_anchors:
                    if _place_ring_inline(_ring_t, _ring_anchors[0]):
                        _ring_done.add(_ri_idx)
                        ring_handled = True
                        progress = True
                        break
            if ring_handled:
                continue

            pi = parent.GetIdx()
            parent_pos = coords[pi]
            child_sym = atom.GetSymbol()
            parent_sym = parent.GetSymbol()

            if parent_sym in _METAL_SET:
                bl = _get_ml_bond_length(parent_sym, child_sym)
            elif child_sym in _METAL_SET:
                bl = _get_ml_bond_length(child_sym, parent_sym)
            else:
                bl = _bond_len(parent_sym, child_sym)

            # VSEPR: direction away from already-placed neighbors of parent
            used_dirs = []
            for pnbr in parent.GetNeighbors():
                if pnbr.GetIdx() in placed and pnbr.GetIdx() != ai:
                    d = coords[pnbr.GetIdx()] - parent_pos
                    d_len = float(np.linalg.norm(d))
                    if d_len > 1e-8:
                        used_dirs.append(d / d_len)

            # Metal avoidance: if child is NOT bonded to any metal,
            # add metal directions as repulsive so BFS pushes away.
            _child_bonded_metals = set()
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() in _METAL_SET:
                    _child_bonded_metals.add(nbr.GetIdx())
            for mi in metal_indices:
                if mi in _child_bonded_metals:
                    continue
                m_to_p = parent_pos - coords[mi]
                m_to_p_len = float(np.linalg.norm(m_to_p))
                if m_to_p_len < 4.0 and m_to_p_len > 1e-8:
                    used_dirs.append(-m_to_p / m_to_p_len)

            if len(used_dirs) == 0:
                direction = np.array([0.0, 0.0, 1.0])
            elif len(used_dirs) == 1:
                d0 = used_dirs[0]
                ref = (np.array([1.0, 0.0, 0.0]) if abs(d0[0]) < 0.9
                       else np.array([0.0, 1.0, 0.0]))
                perp = np.cross(d0, ref)
                perp = perp / max(float(np.linalg.norm(perp)), 1e-12)
                direction = (-d0 * np.cos(np.radians(70.5))
                             + perp * np.sin(np.radians(70.5)))
            elif len(used_dirs) == 2:
                avg = (used_dirs[0] + used_dirs[1]) / 2.0
                avg_len = float(np.linalg.norm(avg))
                if avg_len > 1e-8:
                    direction = -avg / avg_len
                else:
                    perp = np.cross(used_dirs[0], used_dirs[1])
                    p_len = float(np.linalg.norm(perp))
                    direction = (perp / p_len if p_len > 1e-8
                                 else np.array([0.0, 0.0, 1.0]))
            else:
                avg = sum(used_dirs) / len(used_dirs)
                avg_len = float(np.linalg.norm(avg))
                direction = (-avg / avg_len if avg_len > 1e-8
                             else np.array([0.0, 0.0, 1.0]))

            d_norm = float(np.linalg.norm(direction))
            if d_norm > 1e-8:
                direction = direction / d_norm

            coords[ai] = parent_pos + bl * direction
            placed.add(ai)
            progress = True

        if not progress:
            break

    # Handle orphan atoms
    rng = np.random.default_rng(42)
    for ai in range(n_atoms):
        if ai not in placed:
            coords[ai] = rng.standard_normal(3) * 3.0

    # ---- Bond-length relaxation (fix topology) ----
    # BFS only uses one parent per atom; bonds to other placed atoms may be
    # stretched.  Iteratively correct all bond lengths toward their targets.
    bond_targets: List[Tuple[int, int, float]] = []
    for bond in mol.GetBonds():
        bi = bond.GetBeginAtomIdx()
        bj = bond.GetEndAtomIdx()
        si = mol.GetAtomWithIdx(bi).GetSymbol()
        sj = mol.GetAtomWithIdx(bj).GetSymbol()
        if si in _METAL_SET or sj in _METAL_SET:
            continue  # skip metal bonds - already placed correctly
        target_bl = _bond_len(si, sj)
        bond_targets.append((bi, bj, target_bl))

    # Spring-based relaxation: accumulate all bond forces, then apply
    for _relax_pass in range(120):
        forces = np.zeros_like(coords)
        max_err = 0.0
        for bi, bj, target_bl in bond_targets:
            diff = coords[bj] - coords[bi]
            d = float(np.linalg.norm(diff))
            err = abs(d - target_bl)
            if err < 0.01:
                continue
            if d < 1e-8:
                diff = rng.standard_normal(3)
                d = float(np.linalg.norm(diff))
            unit = diff / d
            # Spring force proportional to displacement
            force = 0.3 * (d - target_bl) * unit
            forces[bi] += force
            forces[bj] -= force
            if err > max_err:
                max_err = err
        if max_err < 0.05:
            break
        # Apply forces — hapto ring atoms and metals are frozen
        for ai in range(n_atoms):
            if mol.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET:
                continue
            if ai in all_hapto_atoms:
                continue
            coords[ai] = coords[ai] + forces[ai]

    # ---- Clash resolution (push apart, preserve intra-ring geometry) ----
    # Build metal bond set for metal proximity check
    metal_bonded: set = set()
    for mi in metal_indices:
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            metal_bonded.add((mi, nbr.GetIdx()))
            metal_bonded.add((nbr.GetIdx(), mi))

    rng_clash = np.random.default_rng(123)
    for _pass in range(30):
        moved = False
        for i in range(n_atoms):
            is_metal_i = mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET
            for j in range(i + 1, n_atoms):
                is_metal_j = mol.GetAtomWithIdx(j).GetSymbol() in _METAL_SET
                if is_metal_i and is_metal_j:
                    continue
                diff = coords[j] - coords[i]
                d = float(np.linalg.norm(diff))
                # Metal proximity: non-bonded atoms too close to metal
                if (is_metal_i or is_metal_j) and (i, j) not in metal_bonded:
                    min_ml = 2.8  # minimum non-bonded M-X distance
                    if d < min_ml:
                        if d < 1e-8:
                            direction = rng_clash.standard_normal(3)
                            direction /= max(np.linalg.norm(direction), 1e-12)
                            push = 0.5 * direction
                        else:
                            push = 0.8 * (min_ml - d) / d * diff
                        non_metal = j if is_metal_i else i
                        metal_atom = i if is_metal_i else j
                        sign = 1.0 if non_metal == j else -1.0
                        # Allow push if atom is not hapto, OR if it's hapto
                        # but not bonded to THIS metal (bimetallic case).
                        is_own_hapto = (non_metal in all_hapto_atoms
                                        and (metal_atom, non_metal) in metal_bonded)
                        if not is_own_hapto:
                            coords[non_metal] = coords[non_metal] + sign * push
                            moved = True
                    continue
                if is_metal_i or is_metal_j:
                    continue
                if d < 0.8:
                    if d < 1e-8:
                        direction = rng_clash.standard_normal(3)
                        direction /= max(np.linalg.norm(direction), 1e-12)
                        push = 0.4 * direction
                    else:
                        push = 0.5 * (0.8 - d) / d * diff
                    gi = atom_to_hapto_group.get(i, -1)
                    gj = atom_to_hapto_group.get(j, -2)
                    same_group = (gi == gj and gi >= 0)
                    if not same_group:
                        can_move_i = i not in all_hapto_atoms
                        can_move_j = j not in all_hapto_atoms
                        if can_move_i:
                            coords[i] = coords[i] - push * 0.5
                            moved = True
                        if can_move_j:
                            coords[j] = coords[j] + push * 0.5
                            moved = True
        if not moved:
            break

    # ---- Rigid-body separation for multi-metal systems ----
    # If hapto atoms of one metal intrude into another metal's coordination
    # sphere, translate the entire hapto fragment (metal + all its bonded
    # atoms) as a rigid body to increase M-M distance.
    if len(metal_indices) > 1:
        min_nonbonded_ml = 2.8
        for _sep_pass in range(5):
            any_moved = False
            for mi in metal_indices:
                mi_hapto = set()
                for grp in by_metal.get(mi, []):
                    mi_hapto.update(grp)
                if not mi_hapto:
                    continue
                for mj in metal_indices:
                    if mj == mi:
                        continue
                    mj_pos = coords[mj]
                    worst_intrusion = 0.0
                    push_dir = np.zeros(3)
                    for ha in mi_hapto:
                        d = float(np.linalg.norm(coords[ha] - mj_pos))
                        if d < min_nonbonded_ml and (mj, ha) not in metal_bonded:
                            intrusion = min_nonbonded_ml - d
                            if intrusion > worst_intrusion:
                                worst_intrusion = intrusion
                                v = coords[ha] - mj_pos
                                vl = float(np.linalg.norm(v))
                                push_dir = v / vl if vl > 1e-8 else np.array([1, 0, 0])
                    if worst_intrusion > 0.05:
                        shift = push_dir * worst_intrusion * 1.2
                        mi_frag = {mi} | mi_hapto
                        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
                            mi_frag.add(nbr.GetIdx())
                        for ai in mi_frag:
                            coords[ai] = coords[ai] + shift
                        any_moved = True
            if not any_moved:
                break

    # ---- Write conformer to mol ----
    conf = Chem.Conformer(n_atoms)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, Point3D(
            float(coords[i, 0]), float(coords[i, 1]), float(coords[i, 2])))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return True



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


def _prepare_mol_for_embedding(smiles: str, hapto_approx: bool = False):
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
            if hapto_approx:
                mol = _fix_hapto_donor_h(mol)
        except Exception:
            try:
                mol = Chem.AddHs(mol, addCoords=False)
                if hapto_approx:
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
            if hapto_approx:
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
        primary_ids = list(
            AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=_params_for_seed(seed))
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

    # Fallback: embed on a temporary de-aromatized copy and transfer coords.
    tmp = _dearomatized_embedding_copy(mol)
    if tmp is None:
        return []
    try:
        tmp_ids = list(
            AllChem.EmbedMultipleConfs(tmp, numConfs=num_confs, params=_params_for_seed(seed))
        )
    except Exception as emb_exc:
        logger.debug("Fallback ETKDG embedding failed (seed=%s): %s", seed, emb_exc)
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
            "Embedded %d conformers via dearomatized fallback (seed=%s).",
            len(transferred), seed,
        )
    return transferred


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
        # For multi-metal complexes with bridging O atoms, O-O perception
        # by OB is expected (two oxo/hydroxo on neighbouring metals) — skip.
        metal_count = sum(1 for a in orig_mol.GetAtoms() if a.GetSymbol() in _METAL_SET)
        if metal_count >= 2:
            o_on_metal = 0
            for a in orig_mol.GetAtoms():
                if a.GetAtomicNum() == 8 and any(
                    n.GetSymbol() in _METAL_SET for n in a.GetNeighbors()
                ):
                    o_on_metal += 1
            if o_on_metal >= 2:
                orig_homo.add(frozenset([8, 8]))

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


def _organic_graph_signature(
    symbol_by_idx: Dict[int, str],
    adj: Dict[int, set],
    wl_rounds: int = 3,
) -> "frozenset | None":
    """Return a WL-like multiset signature for disconnected organic graphs.

    The signature is atom-order independent and captures connectivity pattern
    per component much stronger than plain element-count fragment signatures.
    Bond order is intentionally ignored (distance-perception ambiguity in XYZ);
    connectivity preservation is the primary target.
    """
    try:
        comp_mult: Dict[tuple, int] = {}
        visited: set = set()
        for start in sorted(adj.keys()):
            if start in visited:
                continue
            comp: List[int] = []
            stack = [start]
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                comp.append(node)
                stack.extend(adj.get(node, set()) - visited)

            labels: Dict[int, str] = {
                i: str(symbol_by_idx.get(i, '?')) for i in comp
            }
            for _ in range(max(1, int(wl_rounds))):
                new_labels: Dict[int, str] = {}
                for i in comp:
                    neigh = sorted(
                        labels.get(j, '?') for j in adj.get(i, set()) if j in labels
                    )
                    new_labels[i] = labels[i] + "|" + ",".join(neigh)
                labels = new_labels

            comp_sig = (
                len(comp),
                tuple(sorted(symbol_by_idx.get(i, '?') for i in comp)),
                tuple(sorted(len(adj.get(i, set())) for i in comp)),
                tuple(sorted(labels[i] for i in comp)),
            )
            comp_mult[comp_sig] = comp_mult.get(comp_sig, 0) + 1

        return frozenset(comp_mult.items())
    except Exception:
        return None


def _component_stats_from_adj(adj: Dict[int, set]) -> Optional[Tuple[int, int, int]]:
    """Return ``(n_components, largest_component_size, n_nodes)`` for adjacency map."""
    try:
        if not adj:
            return 0, 0, 0
        visited: set = set()
        n_comp = 0
        largest = 0
        for start in adj.keys():
            if start in visited:
                continue
            n_comp += 1
            stack = [start]
            size = 0
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                size += 1
                stack.extend(adj.get(node, set()) - visited)
            if size > largest:
                largest = size
        return n_comp, largest, len(adj)
    except Exception:
        return None


def _heavy_graph_edges_smiles(
    smiles: str,
) -> Optional[Tuple[Dict[int, str], set]]:
    """Return non-metal heavy symbols and bond-order-neutral edge set from SMILES."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass

        heavy_symbols: Dict[int, str] = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() <= 1:
                continue
            if atom.GetSymbol() in _METAL_SET:
                continue
            heavy_symbols[atom.GetIdx()] = atom.GetSymbol()

        edges: set = set()
        for bond in mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            if bi not in heavy_symbols or bj not in heavy_symbols:
                continue
            edges.add(tuple(sorted((int(bi), int(bj)))))
        return heavy_symbols, edges
    except Exception:
        return None


def _heavy_graph_edges_xyz(
    xyz_delfin: str,
    expected_heavy_symbols: Dict[int, str],
) -> Optional[set]:
    """Return bond-order-neutral non-metal heavy edge set perceived from XYZ."""
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return None

        for idx, sym in expected_heavy_symbols.items():
            if idx >= len(lines):
                return None
            parts = lines[idx].split()
            if not parts or parts[0] != sym:
                return None

        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        ob_mol = pybel.readstring('xyz', std_xyz).OBMol
        try:
            from openbabel import openbabel as _ob
        except ImportError:
            return None

        heavy_idx_set = set(int(i) for i in expected_heavy_symbols.keys())
        edges: set = set()
        for bond in _ob.OBMolBondIter(ob_mol):
            i1 = int(bond.GetBeginAtomIdx()) - 1
            i2 = int(bond.GetEndAtomIdx()) - 1
            if i1 not in heavy_idx_set or i2 not in heavy_idx_set:
                continue
            edges.add(tuple(sorted((i1, i2))))
        return edges
    except Exception:
        return None


def _heavy_graph_exact_match_ok(
    xyz_delfin: str,
    original_smiles: str,
) -> bool:
    """Return True if non-metal heavy connectivity matches exactly, ignoring bond order."""
    try:
        ref = _heavy_graph_edges_smiles(original_smiles)
        if ref is None:
            return True
        heavy_symbols, ref_edges = ref
        xyz_edges = _heavy_graph_edges_xyz(xyz_delfin, heavy_symbols)
        if xyz_edges is None:
            return True
        if xyz_edges == ref_edges:
            return True

        missing = sorted(ref_edges - xyz_edges)
        added = sorted(xyz_edges - ref_edges)
        if missing:
            logger.debug(
                "Heavy-graph mismatch: missing %d bond(s), first=%s",
                len(missing),
                missing[0],
            )
        if added:
            logger.debug(
                "Heavy-graph mismatch: added %d bond(s), first=%s",
                len(added),
                added[0],
            )
        return False
    except Exception:
        return True


def _shortest_path_length_excluding(
    adj: Dict[int, set],
    start: int,
    goal: int,
    blocked: int,
    max_depth: int = 6,
) -> Optional[int]:
    """Return shortest path length from ``start`` to ``goal`` while skipping ``blocked``."""
    from collections import deque
    if start == goal:
        return 0
    visited = {blocked, start}
    queue = deque([(start, 0)])
    while queue:
        node, depth = queue.popleft()
        if depth >= max_depth:
            continue
        for nbr in adj.get(node, set()):
            if nbr in visited:
                continue
            if nbr == goal:
                return depth + 1
            visited.add(nbr)
            queue.append((nbr, depth + 1))
    return None


def _cycle_size_signature_from_adj(
    adj: Dict[int, set],
    max_cycle_size: int = 8,
) -> Dict[int, Tuple[int, ...]]:
    """Return per-atom cycle-size signatures derived from the adjacency graph."""
    cycle_sizes: Dict[int, set] = {idx: set() for idx in adj}
    for node, neighbors in adj.items():
        nbrs = sorted(neighbors)
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                path_len = _shortest_path_length_excluding(
                    adj,
                    nbrs[i],
                    nbrs[j],
                    blocked=node,
                    max_depth=max(2, int(max_cycle_size) - 2),
                )
                if path_len is None:
                    continue
                cycle_size = path_len + 2
                if 3 <= cycle_size <= max_cycle_size:
                    cycle_sizes[node].add(int(cycle_size))
    return {
        idx: tuple(sorted(sizes))
        for idx, sizes in cycle_sizes.items()
    }


def _heavy_local_signature_multiset(
    symbol_by_idx: Dict[int, str],
    adj: Dict[int, set],
    wl_rounds: int = 2,
) -> Counter:
    """Return a multiset of local heavy-atom graph signatures."""
    cycle_sig = _cycle_size_signature_from_adj(adj)
    labels: Dict[int, str] = {
        idx: f"{symbol_by_idx.get(idx, '?')}|d{len(adj.get(idx, set()))}|c{','.join(map(str, cycle_sig.get(idx, ()) ))}"
        for idx in adj
    }
    for _ in range(max(1, int(wl_rounds))):
        new_labels: Dict[int, str] = {}
        for idx in adj:
            neigh = sorted(labels.get(nbr, '?') for nbr in adj.get(idx, set()))
            new_labels[idx] = labels[idx] + "|" + ";".join(neigh)
        labels = new_labels
    signatures = Counter()
    for idx in adj:
        signatures[
            (
                symbol_by_idx.get(idx, '?'),
                len(adj.get(idx, set())),
                cycle_sig.get(idx, ()),
                labels.get(idx, '?'),
            )
        ] += 1
    return signatures


def _heavy_local_signature_multiset_smiles(
    smiles: str,
) -> Optional[Counter]:
    """Return local heavy-atom environment multiset from SMILES."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        keep = {
            atom.GetIdx(): atom.GetSymbol()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() > 1 and atom.GetSymbol() not in _METAL_SET
        }
        adj: Dict[int, set] = {idx: set() for idx in keep}
        for bond in mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            if bi in adj and bj in adj:
                adj[bi].add(bj)
                adj[bj].add(bi)
        return _heavy_local_signature_multiset(keep, adj)
    except Exception:
        return None


def _heavy_local_signature_multiset_xyz(
    xyz_delfin: str,
    expected_heavy_symbols: Dict[int, str],
) -> Optional[Counter]:
    """Return local heavy-atom environment multiset from XYZ via OB perception."""
    xyz_edges = _heavy_graph_edges_xyz(xyz_delfin, expected_heavy_symbols)
    if xyz_edges is None:
        return None
    adj: Dict[int, set] = {idx: set() for idx in expected_heavy_symbols}
    for i1, i2 in xyz_edges:
        if i1 in adj and i2 in adj:
            adj[i1].add(i2)
            adj[i2].add(i1)
    return _heavy_local_signature_multiset(expected_heavy_symbols, adj)


def _heavy_local_signature_match_ok(
    xyz_delfin: str,
    original_smiles: str,
) -> bool:
    """Return True when local heavy-atom graph environments still match."""
    try:
        ref = _heavy_graph_edges_smiles(original_smiles)
        if ref is None:
            return True
        heavy_symbols, _ref_edges = ref
        sig_smiles = _heavy_local_signature_multiset_smiles(original_smiles)
        sig_xyz = _heavy_local_signature_multiset_xyz(xyz_delfin, heavy_symbols)
        if sig_smiles is None or sig_xyz is None:
            return True
        if sig_smiles == sig_xyz:
            return True
        logger.debug(
            "Heavy local-signature mismatch: smiles=%d xyz=%d",
            sum(sig_smiles.values()),
            sum(sig_xyz.values()),
        )
        return False
    except Exception:
        return True


def _nonmetal_fragment_ids(mol) -> Dict[int, int]:
    """Return connected-component ids for heavy non-metal atoms."""
    if not RDKIT_AVAILABLE or mol is None:
        return {}
    adj: Dict[int, set] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            continue
        adj[atom.GetIdx()] = set()
    for bond in mol.GetBonds():
        bi = bond.GetBeginAtomIdx()
        bj = bond.GetEndAtomIdx()
        if bi in adj and bj in adj:
            adj[bi].add(bj)
            adj[bj].add(bi)

    frag_id_by_atom: Dict[int, int] = {}
    next_frag_id = 0
    for start in sorted(adj):
        if start in frag_id_by_atom:
            continue
        stack = [start]
        while stack:
            node = stack.pop()
            if node in frag_id_by_atom:
                continue
            frag_id_by_atom[node] = next_frag_id
            stack.extend(adj.get(node, set()) - set(frag_id_by_atom))
        next_frag_id += 1
    return frag_id_by_atom


def _metal_aware_coordination_ok(
    mol,
    conf_id: int,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
) -> bool:
    """Return True when metal-centered donor fragments match the intended model."""
    if not RDKIT_AVAILABLE or mol is None:
        return True
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return True

    frag_id_by_atom = _nonmetal_fragment_ids(mol)
    by_metal: Dict[int, List[List[int]]] = {}
    hapto_atom_to_metal: Dict[int, int] = {}
    for metal_idx, group_atoms in (hapto_groups or []):
        by_metal.setdefault(int(metal_idx), []).append(list(group_atoms))
        for atom_idx in group_atoms:
            hapto_atom_to_metal[int(atom_idx)] = int(metal_idx)

    def _dist(i: int, j: int) -> float:
        pi = conf.GetAtomPosition(i)
        pj = conf.GetAtomPosition(j)
        return math.sqrt(
            (pi.x - pj.x) ** 2 + (pi.y - pj.y) ** 2 + (pi.z - pj.z) ** 2
        )

    for metal in mol.GetAtoms():
        if metal.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = metal.GetIdx()
        metal_sym = metal.GetSymbol()
        metal_hapto_groups = by_metal.get(metal_idx, [])
        metal_hapto_atoms = {
            atom_idx for group_atoms in metal_hapto_groups for atom_idx in group_atoms
        }

        # Hapto consistency: each expected group should retain a plausible
        # metal-centroid distance after candidate selection.
        for group_atoms in metal_hapto_groups:
            if not group_atoms:
                continue
            pts = [conf.GetAtomPosition(atom_idx) for atom_idx in group_atoms]
            centroid = (
                sum(p.x for p in pts) / len(pts),
                sum(p.y for p in pts) / len(pts),
                sum(p.z for p in pts) / len(pts),
            )
            mpos = conf.GetAtomPosition(metal_idx)
            mc_dist = math.sqrt(
                (centroid[0] - mpos.x) ** 2
                + (centroid[1] - mpos.y) ** 2
                + (centroid[2] - mpos.z) ** 2
            )
            target_mc = _target_mc_dist(metal_sym, len(group_atoms))
            if abs(mc_dist - target_mc) > 0.70:
                return False

        # Fragment-aware donor consistency for non-hapto donors.
        expected_frag_donors: Dict[int, Counter] = {}
        expected_donor_atom_indices: set = set()
        for nbr in metal.GetNeighbors():
            donor_idx = nbr.GetIdx()
            if nbr.GetAtomicNum() <= 1 or nbr.GetSymbol() in _METAL_SET:
                continue
            if donor_idx in metal_hapto_atoms:
                continue
            expected_donor_atom_indices.add(donor_idx)
            frag_id = frag_id_by_atom.get(donor_idx)
            if frag_id is None:
                continue
            expected_frag_donors.setdefault(frag_id, Counter())[nbr.GetSymbol()] += 1

        if expected_frag_donors:
            observed_frag_donors: Dict[int, Counter] = {}
            for atom in mol.GetAtoms():
                atom_idx = atom.GetIdx()
                if atom_idx == metal_idx:
                    continue
                if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                    continue
                if atom_idx in metal_hapto_atoms:
                    continue
                frag_id = frag_id_by_atom.get(atom_idx)
                if frag_id is None:
                    continue
                dist = _dist(metal_idx, atom_idx)
                target = _get_ml_bond_length(metal_sym, atom.GetSymbol())
                # Count atoms that genuinely occupy the coordination shell.
                if dist <= target + 0.45:
                    observed_frag_donors.setdefault(frag_id, Counter())[atom.GetSymbol()] += 1

            for frag_id, expected_counter in expected_frag_donors.items():
                observed_counter = observed_frag_donors.get(frag_id, Counter())
                for donor_sym, expected_count in expected_counter.items():
                    if observed_counter.get(donor_sym, 0) < expected_count:
                        return False

            # For secondary/non-hapto metals, reject foreign fragments that
            # enter the first coordination shell unexpectedly.
            if not metal_hapto_groups:
                expected_frag_ids = set(expected_frag_donors.keys())
                for atom in mol.GetAtoms():
                    atom_idx = atom.GetIdx()
                    if atom_idx == metal_idx:
                        continue
                    if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                        continue
                    frag_id = frag_id_by_atom.get(atom_idx)
                    if frag_id is None or frag_id in expected_frag_ids:
                        continue
                    dist = _dist(metal_idx, atom_idx)
                    intrusion_limit = _get_ml_bond_length(metal_sym, atom.GetSymbol()) + 0.15
                    if dist <= intrusion_limit:
                        return False

    return True


def _organic_fragment_signature(smiles: str) -> "frozenset | None":
    """Return an order-independent connectivity signature for organic fragments."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass

        keep = [
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() > 1 and a.GetSymbol() not in _METAL_SET
        ]
        symbol_by_idx: Dict[int, str] = {
            idx: mol.GetAtomWithIdx(idx).GetSymbol() for idx in keep
        }
        adj: Dict[int, set] = {idx: set() for idx in keep}
        for bond in mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            if bi in adj and bj in adj:
                adj[bi].add(bj)
                adj[bj].add(bi)
        return _organic_graph_signature(symbol_by_idx, adj)
    except Exception:
        return None


def _heavy_component_stats_smiles(smiles: str) -> Optional[Tuple[int, int, int]]:
    """Return heavy-atom component stats for SMILES, including metal atoms."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        keep = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        adj: Dict[int, set] = {idx: set() for idx in keep}
        for bond in mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            if bi in adj and bj in adj:
                adj[bi].add(bj)
                adj[bj].add(bi)
        return _component_stats_from_adj(adj)
    except Exception:
        return None


def _organic_fragment_signature_xyz(xyz_delfin: str) -> "frozenset | None":
    """Return organic connectivity signature for DELFIN XYZ (via OB perception)."""
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return None
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        ob_mol = pybel.readstring('xyz', std_xyz).OBMol
        try:
            from openbabel import openbabel as _ob
        except ImportError:
            return None

        metal_ob_idx = {a.GetIdx() for a in _ob.OBMolAtomIter(ob_mol)
                        if a.GetAtomicNum() in _METAL_ATOMICNUMS}

        # Build heavy-atom adjacency excluding metals.
        n_atoms = ob_mol.NumAtoms()
        adj: dict = {i: set() for i in range(1, n_atoms + 1)
                     if ob_mol.GetAtom(i).GetAtomicNum() not in _METAL_ATOMICNUMS
                     and ob_mol.GetAtom(i).GetAtomicNum() not in (0, 1)}
        symbol_by_idx: Dict[int, str] = {
            i: (
                Chem.GetPeriodicTable().GetElementSymbol(ob_mol.GetAtom(i).GetAtomicNum())
                if RDKIT_AVAILABLE else str(ob_mol.GetAtom(i).GetAtomicNum())
            )
            for i in adj
        }

        for bond in _ob.OBMolBondIter(ob_mol):
            i1 = bond.GetBeginAtomIdx()
            i2 = bond.GetEndAtomIdx()
            if i1 in adj and i2 in adj:
                adj[i1].add(i2)
                adj[i2].add(i1)

        return _organic_graph_signature(symbol_by_idx, adj)
    except Exception:
        return None


def _heavy_component_stats_xyz(xyz_delfin: str) -> Optional[Tuple[int, int, int]]:
    """Return heavy-atom component stats for XYZ via Open Babel, incl. metals."""
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return None
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        ob_mol = pybel.readstring('xyz', std_xyz).OBMol
        try:
            from openbabel import openbabel as _ob
        except ImportError:
            return None

        n_atoms = ob_mol.NumAtoms()
        adj: Dict[int, set] = {
            i: set()
            for i in range(1, n_atoms + 1)
            if ob_mol.GetAtom(i).GetAtomicNum() not in (0, 1)
        }
        for bond in _ob.OBMolBondIter(ob_mol):
            i1 = bond.GetBeginAtomIdx()
            i2 = bond.GetEndAtomIdx()
            if i1 in adj and i2 in adj:
                adj[i1].add(i2)
                adj[i2].add(i1)
        return _component_stats_from_adj(adj)
    except Exception:
        return None


def _global_heavy_connectivity_ok(
    xyz_delfin: str,
    original_smiles: str,
    max_extra_components: int = 1,
    min_largest_frac: float = 0.70,
) -> bool:
    """Reject severely fragmented heavy-atom graphs vs original SMILES.

    Complements ``_fragment_topology_ok``: the organic-only signature can miss
    cases where metal-ligand connectivity collapses while organic fragments stay
    internally intact.
    """
    orig_stats = _heavy_component_stats_smiles(original_smiles)
    xyz_stats = _heavy_component_stats_xyz(xyz_delfin)
    if orig_stats is None or xyz_stats is None:
        return True

    o_comp, o_largest, o_total = orig_stats
    x_comp, x_largest, _x_total = xyz_stats

    if x_comp > o_comp + int(max_extra_components):
        logger.debug(
            "Heavy-graph fragmentation: SMILES has %d component(s), XYZ has %d",
            o_comp, x_comp,
        )
        return False

    # If original is mostly one connected heavy graph, demand that XYZ keeps
    # a comparably large main component.
    if o_total > 0 and o_largest >= max(3, int(math.ceil(0.80 * o_total))):
        min_keep = max(2, int(math.floor(float(min_largest_frac) * float(o_largest))))
        if x_largest < min_keep:
            logger.debug(
                "Heavy-graph largest component collapsed: SMILES=%d, XYZ=%d (min=%d)",
                o_largest, x_largest, min_keep,
            )
            return False

    return True


def _verify_topology_from_graph(
    xyz_delfin: str,
    mol_template,
) -> bool:
    """Graph-based topology verification — no OB perception, no roundtrip.

    Checks the XYZ coordinates directly against the template molecular
    graph (from SMILES).  This is the SINGLE authoritative check for
    structural integrity.

    Three rules:
    1. Every bond in the template graph must have a reasonable distance
       in the XYZ:
       - M-D bonds: within [0.75, 1.60] × ``_get_ml_bond_length``
       - M-M bonds: within [0.70, 1.80] × ``_METAL_METAL_BOND_LENGTHS``
       - Covalent bonds (non-metal): < 2.2 Å
    2. No two heavy atoms closer than 0.7 Å (collapsed structure)
    3. No H-H closer than 0.4 Å (collapsed hydrogens)
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return True
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        n_atoms = mol_template.GetNumAtoms()
        if len(lines) != n_atoms:
            # Atom count mismatch → try AddHs fallback
            try:
                mol_h = Chem.AddHs(mol_template)
                if len(lines) == mol_h.GetNumAtoms():
                    mol_template = mol_h
                    n_atoms = mol_h.GetNumAtoms()
                else:
                    return True  # can't validate → permissive
            except Exception:
                return True

        coords: List[Tuple[float, float, float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return True
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

        # Rule 1: Every bond must have a reasonable distance.
        # For bridging donors: collect ALL M-D distances, pass if at
        # least one metal is within normal range (a bridging atom can't
        # be at ideal distance from ALL metals simultaneously).
        bridging_donor_bonds: Dict[int, List[Tuple[int, str, float, float]]] = {}

        for bond in mol_template.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            i1, i2 = a1.GetIdx(), a2.GetIdx()
            s1, s2 = a1.GetSymbol(), a2.GetSymbol()
            dx = coords[i1][0] - coords[i2][0]
            dy = coords[i1][1] - coords[i2][1]
            dz = coords[i1][2] - coords[i2][2]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)

            is_metal_1 = s1 in _METAL_SET
            is_metal_2 = s2 in _METAL_SET

            if is_metal_1 and is_metal_2:
                mm_key = frozenset({s1, s2})
                ideal = _METAL_METAL_BOND_LENGTHS.get(mm_key)
                if ideal is None:
                    r1 = _COVALENT_RADII.get(s1)
                    r2 = _COVALENT_RADII.get(s2)
                    ideal = (r1 + r2 + 0.3) if r1 and r2 else 2.5
                if d < 0.70 * ideal or d > 1.80 * ideal:
                    return False
            elif is_metal_1 or is_metal_2:
                m_sym = s1 if is_metal_1 else s2
                d_sym = s2 if is_metal_1 else s1
                d_atom = a2 if is_metal_1 else a1
                m_idx = i1 if is_metal_1 else i2
                if a1.GetAtomicNum() <= 1 or a2.GetAtomicNum() <= 1:
                    continue
                ideal = float(_get_ml_bond_length(m_sym, d_sym))
                if ideal <= 0:
                    continue
                n_metal_nbrs = sum(
                    1 for nbr in d_atom.GetNeighbors()
                    if nbr.GetSymbol() in _METAL_SET
                )
                if n_metal_nbrs >= 2:
                    # Bridging: defer to collective check below
                    bridging_donor_bonds.setdefault(
                        d_atom.GetIdx(), []
                    ).append((m_idx, m_sym, d, ideal))
                else:
                    # Terminal donor: strict check
                    if d < 0.75 * ideal or d > 1.60 * ideal:
                        return False
            else:
                if a1.GetAtomicNum() <= 1 or a2.GetAtomicNum() <= 1:
                    if d > 1.8:
                        return False
                else:
                    # Covalent bonds near metal centers (e.g. Cp ring C-C)
                    # can stretch slightly due to coordination effects.
                    if d > 2.4:
                        return False

        # Bridging donors: ALL M-D bonds must be within [0.60, 2.00].
        # A bridging atom IS bonded to all its metals per the SMILES —
        # if any bond is broken (>2.0× ideal), the topology is lost.
        for d_idx, metal_entries in bridging_donor_bonds.items():
            for _m_idx, _m_sym, d, ideal in metal_entries:
                ratio = d / ideal
                if ratio < 0.60 or ratio > 2.00:
                    return False

        # Rule 2+3: No collapsed heavy atoms or hydrogens.
        heavy_indices = [
            i for i in range(n_atoms)
            if mol_template.GetAtomWithIdx(i).GetAtomicNum() > 1
        ]
        for i in range(len(heavy_indices)):
            xi, yi, zi = coords[heavy_indices[i]]
            for j in range(i + 1, min(i + 50, len(heavy_indices))):
                xj, yj, zj = coords[heavy_indices[j]]
                dsq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
                if dsq < 0.49:  # 0.7²
                    return False

        # Rule 4: Pi-ring planarity — reject if any ring with sp2 atoms
        # has severe out-of-plane distortion (>0.5 Å).
        try:
            ring_info = mol_template.GetRingInfo()
            if ring_info is not None:
                for ring in ring_info.AtomRings():
                    if len(ring) < 5 or len(ring) > 7:
                        continue
                    # Check if ring has sp2 character (aromatic or conjugated)
                    n_sp2 = sum(
                        1 for ri in ring
                        if mol_template.GetAtomWithIdx(ri).GetIsAromatic()
                        or mol_template.GetAtomWithIdx(ri).GetHybridization()
                        in (Chem.rdchem.HybridizationType.SP2,
                            Chem.rdchem.HybridizationType.SP)
                    )
                    if n_sp2 < len(ring) * 0.6:
                        continue  # not a pi ring
                    # Compute planarity deviation via SVD
                    try:
                        import numpy as _np
                        pts = _np.array([coords[ri] for ri in ring])
                        centered = pts - pts.mean(axis=0)
                        _u, _s, vh = _np.linalg.svd(centered, full_matrices=False)
                        deviations = _np.abs(centered @ vh[-1])
                        if deviations.max() > 0.5:
                            return False
                    except Exception:
                        pass
        except Exception:
            pass

        # Rule 5: Inter-ligand proximity — non-bonded heavy atoms from
        # DIFFERENT ligand fragments must not be closer than 1.2 Å.
        try:
            metal_idxs = {
                a.GetIdx() for a in mol_template.GetAtoms()
                if a.GetSymbol() in _METAL_SET
            }
            non_metal = {
                a.GetIdx() for a in mol_template.GetAtoms()
                if a.GetSymbol() not in _METAL_SET and a.GetAtomicNum() > 1
            }
            adj_nm: Dict[int, set] = {i: set() for i in non_metal}
            for bond in mol_template.GetBonds():
                bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if bi in non_metal and bj in non_metal:
                    adj_nm[bi].add(bj)
                    adj_nm[bj].add(bi)
            visited: set = set()
            frag_id: Dict[int, int] = {}
            fid = 0
            for start in sorted(non_metal):
                if start in visited:
                    continue
                stack = [start]
                while stack:
                    node = stack.pop()
                    if node in visited:
                        continue
                    visited.add(node)
                    frag_id[node] = fid
                    for nb in adj_nm.get(node, ()):
                        if nb not in visited:
                            stack.append(nb)
                fid += 1
            if fid >= 2:
                # Check inter-fragment heavy-atom distances
                nm_list = sorted(non_metal)
                for i in range(len(nm_list)):
                    fi = frag_id.get(nm_list[i], -1)
                    xi, yi, zi = coords[nm_list[i]]
                    for j in range(i + 1, min(i + 30, len(nm_list))):
                        fj = frag_id.get(nm_list[j], -1)
                        if fi == fj:
                            continue
                        xj, yj, zj = coords[nm_list[j]]
                        dsq = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
                        if dsq < 1.44:  # 1.2² = 1.44
                            return False
        except Exception:
            pass

        return True
    except Exception:
        return False


def _fragment_topology_relaxed_fallback_ok(
    xyz_delfin: str,
    original_smiles: str,
) -> bool:
    """Allow only mildly fragmented fallback candidates.

    Used when strict fragment-topology checks reject all sampling conformers:
    keep only structures where the organic fragment signature is unchanged and
    heavy-atom connectivity degradation is limited.
    """
    try:
        orig_sig = _organic_fragment_signature(original_smiles)
        xyz_sig = _organic_fragment_signature_xyz(xyz_delfin)
        if orig_sig is None or xyz_sig is None or orig_sig != xyz_sig:
            return False
        return _global_heavy_connectivity_ok(
            xyz_delfin,
            original_smiles,
            max_extra_components=3,
            min_largest_frac=0.85,
        )
    except Exception:
        return False


def _metal_donor_distances_realistic(
    xyz_delfin: str,
    mol_template=None,
    min_abs_ml: float = 1.70,
    max_frac: float = 1.50,
) -> bool:
    """Reject structures with unphysical metal-donor distances.

    Two checks per metal atom in the XYZ:

    1. No heavy atom (any element, coordinating or not) closer than
       ``min_abs_ml`` to any metal — catches collapsed alt-binding
       artifacts where OB perception created false short bonds.
    2. When ``mol_template`` is available, each bonded non-metal donor
       must stay below ``max_frac`` × ``_get_ml_bond_length`` to flag
       broken connectivity.
    """
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        atoms: List[Tuple[str, Tuple[float, float, float]]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return True
            atoms.append(
                (parts[0], (float(parts[1]), float(parts[2]), float(parts[3])))
            )
        metals = [
            (sym, pos) for sym, pos in atoms if sym in _METAL_SET
        ]
        if not metals:
            return True

        min_abs_sq = float(min_abs_ml) * float(min_abs_ml)
        for m_sym, m_pos in metals:
            for a_sym, a_pos in atoms:
                if a_sym == m_sym and a_pos == m_pos:
                    continue
                if a_sym == 'H':
                    continue
                dx = m_pos[0] - a_pos[0]
                dy = m_pos[1] - a_pos[1]
                dz = m_pos[2] - a_pos[2]
                dsq = dx * dx + dy * dy + dz * dz
                if 0.0 < dsq < min_abs_sq:
                    return False

        if mol_template is not None and RDKIT_AVAILABLE:
            if len(atoms) != mol_template.GetNumAtoms():
                return True
            coords = [pos for _sym, pos in atoms]
            for atom in mol_template.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    continue
                m_idx = atom.GetIdx()
                m_sym_t = atom.GetSymbol()
                for nbr in atom.GetNeighbors():
                    if nbr.GetSymbol() in _METAL_SET:
                        continue
                    if nbr.GetAtomicNum() <= 1:
                        continue
                    n_idx = nbr.GetIdx()
                    dx = coords[m_idx][0] - coords[n_idx][0]
                    dy = coords[m_idx][1] - coords[n_idx][1]
                    dz = coords[m_idx][2] - coords[n_idx][2]
                    d = math.sqrt(dx * dx + dy * dy + dz * dz)
                    ideal = float(_get_ml_bond_length(m_sym_t, nbr.GetSymbol()))
                    if ideal <= 0:
                        continue
                    if d > float(max_frac) * ideal:
                        return False
        return True
    except Exception:
        return False


def _verify_metal_connectivity(
    xyz_delfin: str,
    mol_template,
    max_donor_frac: float = 1.60,
    min_donor_frac: float = 0.75,
    min_nonddonor_frac: float = 0.80,
) -> bool:
    """Verify that the XYZ preserves the metal-donor connectivity from the template.

    Fundamental principles:
    1. Every defined donor must be within ``max_donor_frac × ideal`` of
       its metal (not drifted away).
    2. Every defined donor must be farther than ``min_donor_frac × ideal``
       from its metal (not collapsed onto it).
    3. Bridging donors (bonded to 2+ metals) only need to satisfy the
       distance window for *at least one* of their metals — they cannot
       be at ideal distance from all metals simultaneously.
    4. Non-bonded atoms must not have collapsed onto a metal
       (< min_nondonor_frac × ideal).
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return True
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if len(lines) != mol_template.GetNumAtoms():
            return True
        coords: List[Tuple[float, float, float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return True
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

        # Pre-compute which donor atoms are bridging (bonded to 2+ metals)
        # and collect all (metal_idx, donor_idx) pairs for deferred bridging check.
        bridging_pairs: Dict[int, List[Tuple[int, str, float]]] = {}
        for atom in mol_template.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() <= 1 or nbr.GetSymbol() in _METAL_SET:
                    continue
                n_metal_nbrs = sum(
                    1 for nn in nbr.GetNeighbors()
                    if nn.GetSymbol() in _METAL_SET
                )
                if n_metal_nbrs >= 2:
                    d_idx = nbr.GetIdx()
                    m_idx = atom.GetIdx()
                    m_sym = atom.GetSymbol()
                    d_sym = nbr.GetSymbol()
                    mx, my, mz = coords[m_idx]
                    dx, dy, dz = coords[d_idx]
                    d = math.sqrt((mx - dx) ** 2 + (my - dy) ** 2 + (mz - dz) ** 2)
                    ideal = float(_get_ml_bond_length(m_sym, d_sym))
                    if ideal > 0:
                        bridging_pairs.setdefault(d_idx, []).append((m_idx, m_sym, d / ideal))

        for atom in mol_template.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            m_idx = atom.GetIdx()
            m_sym = atom.GetSymbol()
            mx, my, mz = coords[m_idx]
            bonded_set = {nbr.GetIdx() for nbr in atom.GetNeighbors()}

            for d_idx in bonded_set:
                d_atom = mol_template.GetAtomWithIdx(d_idx)
                if d_atom.GetAtomicNum() <= 1:
                    continue
                d_sym = d_atom.GetSymbol()
                if d_sym in _METAL_SET:
                    continue
                dx, dy, dz = coords[d_idx]
                d = math.sqrt((mx - dx) ** 2 + (my - dy) ** 2 + (mz - dz) ** 2)
                ideal = float(_get_ml_bond_length(m_sym, d_sym))
                if ideal <= 0:
                    continue
                ratio = d / ideal

                # Bridging donors: defer to collective check below
                if d_idx in bridging_pairs:
                    continue

                # Terminal donor: must be within [min_donor_frac, max_donor_frac] × ideal
                if ratio > max_donor_frac or ratio < min_donor_frac:
                    return False

            # Non-donor collapse check — only for atoms ≥3 bonds away from
            # ANY metal.  Atoms within 2 bonds of any metal center are
            # naturally close in macrocyclic/chelating/bimetallic systems
            # (e.g. porphyrin C_alpha, bridging-ligand fragments).
            near_any_metal: set = set()
            for m_atom2 in mol_template.GetAtoms():
                if m_atom2.GetSymbol() not in _METAL_SET:
                    continue
                near_any_metal.add(m_atom2.GetIdx())
                for nbr1 in m_atom2.GetNeighbors():
                    near_any_metal.add(nbr1.GetIdx())
                    for nbr2 in nbr1.GetNeighbors():
                        near_any_metal.add(nbr2.GetIdx())
            for other in mol_template.GetAtoms():
                o_idx = other.GetIdx()
                if o_idx == m_idx or o_idx in near_any_metal:
                    continue
                if other.GetAtomicNum() <= 1:
                    continue
                if other.GetSymbol() in _METAL_SET:
                    continue
                ox, oy, oz = coords[o_idx]
                d = math.sqrt((mx - ox) ** 2 + (my - oy) ** 2 + (mz - oz) ** 2)
                ideal = float(_get_ml_bond_length(m_sym, other.GetSymbol()))
                if ideal <= 0:
                    continue
                if d < min_nonddonor_frac * ideal:
                    return False

        # Bridging donor check: each bridging donor must satisfy the distance
        # window for AT LEAST ONE of its metal partners.
        for d_idx, metal_entries in bridging_pairs.items():
            any_ok = False
            for _m_idx, _m_sym, ratio in metal_entries:
                if min_donor_frac <= ratio <= max_donor_frac:
                    any_ok = True
                    break
            if not any_ok:
                return False

        return True
    except Exception:
        return False


def _xyz_passes_final_geometry_checks(
    xyz_delfin: str,
    mol_template,
    skip_angle_check: bool = False,
) -> bool:
    """Final geometry sanity check for accepted XYZ outputs.

    Intended for relaxed-fallback candidates: keep only structures that still
    pass hard geometric plausibility checks when mapped back to the template
    molecular graph.

    When ``skip_angle_check`` is True, only covalent bond distortion is
    checked (not L-M-L angles). This is correct for topology-built
    structures whose donor positions are set by ideal-polyhedron vectors
    — their angles may not satisfy the sampling-oriented thresholds in
    ``_has_bad_geometry``.
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return True
    try:
        # Ensure mol has explicit H so atom count matches XYZ (which
        # includes H atoms from the topology builder).
        mol_tmp = Chem.RWMol(mol_template)
        mol_tmp.RemoveAllConformers()
        conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_delfin)
        if conf is None:
            # XYZ has H but mol doesn't → add explicit H and retry.
            try:
                mol_h = Chem.AddHs(mol_template)
                mol_h = Chem.RWMol(mol_h)
                mol_h.RemoveAllConformers()
                conf = _xyz_to_rdkit_conformer(mol_h.GetMol(), xyz_delfin)
                if conf is None:
                    return True  # Can't validate → permissive
                cid = mol_h.AddConformer(conf, assignId=True)
                if _has_severe_covalent_distortion(mol_h.GetMol(), cid):
                    return False
                if _has_bad_geometry(mol_h.GetMol(), cid):
                    return False
                return True
            except Exception:
                return True  # Can't validate → permissive
        cid = mol_tmp.AddConformer(conf, assignId=True)
        if _has_severe_covalent_distortion(mol_tmp.GetMol(), cid):
            return False
        if not skip_angle_check and _has_bad_geometry(mol_tmp.GetMol(), cid):
            return False
        return True
    except Exception:
        return False


def _fragment_topology_ok(xyz_delfin: str, original_smiles: str) -> bool:
    """Return True if topology is consistent with the original SMILES."""
    orig_sig = _organic_fragment_signature(original_smiles)
    xyz_sig = _organic_fragment_signature_xyz(xyz_delfin)
    if orig_sig is not None and xyz_sig is not None and orig_sig != xyz_sig:
        def _total_frags(sig):
            return sum(mult for _, mult in sig)

        logger.debug(
            "Organic graph mismatch: SMILES has %d fragment(s), XYZ has %d",
            _total_frags(orig_sig), _total_frags(xyz_sig),
        )
        return False

    if not _global_heavy_connectivity_ok(xyz_delfin, original_smiles):
        return False

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


def _has_unphysical_metal_nonbonded_contact(
    mol,
    conf_id: int,
    min_abs_dist: float = 1.10,
    rel_ml_scale: float = 0.55,
) -> bool:
    """Return True for unrealistically short non-bonded metal-heavy contacts.

    Some generated conformers place a non-coordinating heavy atom extremely
    close to a metal center (e.g. <1.0 A) while no bond exists in the graph.
    These geometries are physically implausible and should be rejected.
    """
    if not RDKIT_AVAILABLE:
        return False
    try:
        conf = mol.GetConformer(conf_id)
        for m_atom in mol.GetAtoms():
            if m_atom.GetSymbol() not in _METAL_SET:
                continue
            m_idx = m_atom.GetIdx()
            m_pos = conf.GetAtomPosition(m_idx)
            m_sym = m_atom.GetSymbol()
            for atom in mol.GetAtoms():
                a_idx = atom.GetIdx()
                if a_idx == m_idx:
                    continue
                if atom.GetAtomicNum() <= 1:
                    continue
                if atom.GetSymbol() in _METAL_SET:
                    continue
                if mol.GetBondBetweenAtoms(m_idx, a_idx) is not None:
                    continue

                a_pos = conf.GetAtomPosition(a_idx)
                d = math.sqrt(
                    (m_pos.x - a_pos.x) ** 2
                    + (m_pos.y - a_pos.y) ** 2
                    + (m_pos.z - a_pos.z) ** 2
                )
                expected_ml = _get_ml_bond_length(m_sym, atom.GetSymbol())
                min_d = max(float(min_abs_dist), float(rel_ml_scale) * float(expected_ml))
                if d < min_d:
                    return True
    except Exception:
        return False
    return False


def _has_unphysical_oco_geometry(
    mol,
    conf_id: int,
    min_oco_angle: float = 100.0,
    max_oco_angle: float = 145.0,
    min_co_dist: float = 1.05,
    max_co_dist: float = 1.45,
) -> bool:
    """Return True if a carboxyl-like O-C-O unit is unrealistically distorted.

    For carbon atoms bound to exactly two oxygen atoms (typical carboxyl/carbonyl
    motif), the O-C-O angle should be trigonal-planar-like (~120 deg), not
    linear like free CO2. Distances are also checked against broad covalent
    ranges to catch collapsed or stretched C-O bonds.
    """
    if not RDKIT_AVAILABLE:
        return False
    try:
        conf = mol.GetConformer(conf_id)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetSymbol() in _METAL_SET:
                continue

            o_neighbors = [
                n for n in atom.GetNeighbors()
                if n.GetAtomicNum() == 8 and n.GetSymbol() not in _METAL_SET
            ]
            if len(o_neighbors) != 2:
                continue

            c_idx = atom.GetIdx()
            c_pos = conf.GetAtomPosition(c_idx)
            o1_idx = o_neighbors[0].GetIdx()
            o2_idx = o_neighbors[1].GetIdx()
            o1_pos = conf.GetAtomPosition(o1_idx)
            o2_pos = conf.GetAtomPosition(o2_idx)

            d1 = math.sqrt(
                (c_pos.x - o1_pos.x) ** 2
                + (c_pos.y - o1_pos.y) ** 2
                + (c_pos.z - o1_pos.z) ** 2
            )
            d2 = math.sqrt(
                (c_pos.x - o2_pos.x) ** 2
                + (c_pos.y - o2_pos.y) ** 2
                + (c_pos.z - o2_pos.z) ** 2
            )
            if d1 < min_co_dist or d1 > max_co_dist or d2 < min_co_dist or d2 > max_co_dist:
                return True

            v1 = (o1_pos.x - c_pos.x, o1_pos.y - c_pos.y, o1_pos.z - c_pos.z)
            v2 = (o2_pos.x - c_pos.x, o2_pos.y - c_pos.y, o2_pos.z - c_pos.z)
            m1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
            m2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
            if m1 < 1e-10 or m2 < 1e-10:
                return True
            cos_a = max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (m1 * m2)))
            angle = math.degrees(math.acos(cos_a))
            if angle < min_oco_angle or angle > max_oco_angle:
                return True
    except Exception:
        return False
    return False


def _has_pi_ring_nonplanarity(
    mol,
    conf_id: int,
    max_ring_rms: float = 0.35,
) -> bool:
    """Return True if unsaturated 5-7 membered rings are strongly non-planar."""
    if not RDKIT_AVAILABLE:
        return False
    try:
        import numpy as np
    except Exception:
        return False
    try:
        conf = mol.GetConformer(conf_id)
        try:
            Chem.GetSymmSSSR(mol)
        except Exception:
            pass
        ri = mol.GetRingInfo()
        if ri is None:
            return False

        for ring in ri.AtomRings():
            if len(ring) < 5 or len(ring) > 7:
                continue
            if any(
                mol.GetAtomWithIdx(i).GetAtomicNum() <= 1
                or mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET
                for i in ring
            ):
                continue

            unsat = 0
            valid_edges = 0
            for i in range(len(ring)):
                a = ring[i]
                b = ring[(i + 1) % len(ring)]
                bond = mol.GetBondBetweenAtoms(a, b)
                if bond is None:
                    continue
                valid_edges += 1
                if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
                    unsat += 1
            if valid_edges < max(3, len(ring) - 1):
                continue
            if unsat < 2:
                continue

            pts = np.array([
                [
                    conf.GetAtomPosition(a).x,
                    conf.GetAtomPosition(a).y,
                    conf.GetAtomPosition(a).z,
                ]
                for a in ring
            ], dtype=float)
            ctr = pts.mean(axis=0)
            q = pts - ctr
            _u, _s, vh = np.linalg.svd(q, full_matrices=False)
            n = vh[-1]
            n_norm = float(np.linalg.norm(n))
            if n_norm < 1e-12:
                continue
            n = n / n_norm
            d = np.abs(q @ n)
            rms = float(np.sqrt(np.mean(d * d)))
            if rms > max_ring_rms:
                return True
    except Exception:
        return False
    return False


def _has_severe_covalent_distortion(
    mol,
    conf_id: int,
    max_abs_bond: float = 2.4,
    max_covalent_scale: float = 1.8,
) -> bool:
    """Return True if non-metal covalent bonds are unrealistically stretched."""
    if not RDKIT_AVAILABLE:
        return False
    try:
        conf = mol.GetConformer(conf_id)
        pt = Chem.GetPeriodicTable()
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetSymbol() in _METAL_SET or a2.GetSymbol() in _METAL_SET:
                continue

            p1 = conf.GetAtomPosition(a1.GetIdx())
            p2 = conf.GetAtomPosition(a2.GetIdx())
            d = math.sqrt(
                (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2
            )
            if d > max_abs_bond:
                return True

            try:
                rc1 = float(pt.GetRcovalent(a1.GetAtomicNum()))
                rc2 = float(pt.GetRcovalent(a2.GetAtomicNum()))
                if rc1 > 0 and rc2 > 0 and d > max_covalent_scale * (rc1 + rc2):
                    return True
            except Exception:
                pass
    except Exception:
        return False
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


def _cn4_geometry_penalties(
    metal_pos,
    coord_positions: List[object],
    angles: List[float],
) -> Tuple[float, float]:
    """Return ``(tetra_pen, square_pen)`` for a 4-coordinate center."""
    tetra_pen = sum(abs(a - 109.5) for a in angles)
    square_targets = [90.0, 90.0, 90.0, 90.0, 180.0, 180.0]
    square_angle_pen = sum(
        abs(a - t) for a, t in zip(sorted(angles), square_targets)
    )

    # Square-planar candidates should also be planar. Add a donor-plane
    # penalty to distinguish flattened tetrahedral-like solutions.
    square_planarity_pen = 0.0
    try:
        import numpy as np

        pts = np.array(
            [[p.x, p.y, p.z] for p in coord_positions],
            dtype=float,
        )
        if pts.shape == (4, 3):
            centroid = pts.mean(axis=0)
            q = pts - centroid
            _u, _s, vh = np.linalg.svd(q, full_matrices=False)
            normal = vh[-1]
            n_norm = float(np.linalg.norm(normal))
            if n_norm > 1e-12:
                normal = normal / n_norm
                donor_dev = np.abs(q @ normal)
                rms_donor = float(np.sqrt(np.mean(donor_dev * donor_dev)))
                mp = np.array([metal_pos.x, metal_pos.y, metal_pos.z], dtype=float)
                metal_dev = abs(float(np.dot(mp - centroid, normal)))
                square_planarity_pen = 35.0 * rms_donor + 45.0 * metal_dev
    except Exception:
        pass

    return tetra_pen, (square_angle_pen + square_planarity_pen)


def _preferred_cn4_geometry_score(
    metal_symbol: str,
    donor_symbols: List[str],
    metal_pos,
    coord_positions: List[object],
    angles: List[float],
) -> float:
    """Return a metal-aware geometry score for 4-coordinate centers."""
    tetra_pen, square_pen = _cn4_geometry_penalties(
        metal_pos,
        coord_positions,
        angles,
    )

    if metal_symbol in {'Pt', 'Pd'}:
        return min(square_pen, tetra_pen + 40.0)
    if metal_symbol == 'Au':
        return min(square_pen, tetra_pen + 28.0)
    if metal_symbol == 'Ni':
        n_n = donor_symbols.count('N')
        n_o = donor_symbols.count('O')
        n_p = donor_symbols.count('P')
        n_c = donor_symbols.count('C')
        # d8 Ni(II) with cyclometalated / N-rich donor sets strongly tends
        # toward square-planar arrangements.
        if n_c >= 1 or (n_n + n_o) >= 3 or n_p >= 2:
            return min(square_pen, tetra_pen + 22.0)
        return min(square_pen, tetra_pen + 10.0)
    if metal_symbol in {'Zn', 'Cd', 'Hg'}:
        return min(tetra_pen, square_pen + 22.0)
    if metal_symbol in {'Cu', 'Ag'}:
        return min(tetra_pen, square_pen + 10.0)

    return min(tetra_pen, square_pen)


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
    if _has_unphysical_metal_nonbonded_contact(mol, conf_id):
        return True
    if _has_unphysical_oco_geometry(mol, conf_id):
        return True
    if _has_pi_ring_nonplanarity(mol, conf_id):
        return True
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
        donor_symbols: List[str] = []
        for nbr in neighbors:
            nbr_pos = conf.GetAtomPosition(nbr.GetIdx())
            dx = nbr_pos.x - metal_pos.x
            dy = nbr_pos.y - metal_pos.y
            dz = nbr_pos.z - metal_pos.z
            dists.append(math.sqrt(dx*dx + dy*dy + dz*dz))
            coord_positions.append(nbr_pos)
            donor_symbols.append(nbr.GetSymbol())

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
            total_penalty += _preferred_cn4_geometry_score(
                atom.GetSymbol(),
                donor_symbols,
                metal_pos,
                coord_positions,
                angles,
            )
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


def _hapto_candidate_topology_ok(
    xyz_delfin: str,
    original_smiles: str,
    mol=None,
    conf_id: int = 0,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
) -> bool:
    """Return True when a hapto candidate preserves the input topology."""
    try:
        if mol is not None and not _metal_aware_coordination_ok(mol, conf_id, hapto_groups):
            return False
        if not _heavy_graph_exact_match_ok(xyz_delfin, original_smiles):
            return False
        if not _heavy_local_signature_match_ok(xyz_delfin, original_smiles):
            return False
        if not _roundtrip_ring_count_ok(xyz_delfin, original_smiles):
            return False
        if not _no_spurious_bonds(xyz_delfin, original_smiles):
            return False
        if not _fragment_topology_ok(xyz_delfin, original_smiles):
            if not _fragment_topology_relaxed_fallback_ok(xyz_delfin, original_smiles):
                return False
        return True
    except Exception:
        # If topology perception is unavailable (for example without OB),
        # keep the candidate and let graph/geometry checks decide.
        return True


def _hapto_geometry_quality_score(
    mol,
    conf_id: int = 0,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
) -> float:
    """Return a hapto-specific geometry penalty (lower = better)."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return 0.0
    try:
        import numpy as np
        conf = mol.GetConformer(conf_id)
    except Exception:
        return 0.0

    def _metal_atomic_number(sym: str) -> int:
        try:
            if RDKIT_AVAILABLE:
                return int(Chem.GetPeriodicTable().GetAtomicNumber(sym))
        except Exception:
            pass
        return 0

    def _is_f_block_like(sym: str) -> bool:
        z = _metal_atomic_number(sym)
        return (
            57 <= z <= 71
            or 89 <= z <= 103
            or sym in {'Y', 'Sc'}
        )

    def _hapto_weight_profile(metal_sym: str, eta: int, n_groups_for_metal: int) -> Dict[str, float]:
        profile = {
            'mc': 180.0,
            'radial': 30.0,
            'plane': 120.0,
            'axis': 40.0,
            'lateral': 35.0,
            'min_pair_angle': 90.0,
            'pair_sep': 0.45,
            'equal_eta_mc_std': 40.0,
        }

        if eta >= 5:
            profile.update({
                'mc': 220.0,
                'radial': 36.0,
                'plane': 165.0,
                'axis': 72.0,
                'lateral': 55.0,
                'min_pair_angle': 118.0,
                'pair_sep': 0.80,
                'equal_eta_mc_std': 55.0,
            })
        elif eta == 4:
            profile.update({
                'mc': 195.0,
                'radial': 32.0,
                'plane': 135.0,
                'axis': 52.0,
                'lateral': 42.0,
                'min_pair_angle': 102.0,
                'pair_sep': 0.60,
                'equal_eta_mc_std': 45.0,
            })
        elif eta == 3:
            profile.update({
                'mc': 150.0,
                'radial': 24.0,
                'plane': 82.0,
                'axis': 24.0,
                'lateral': 20.0,
                'min_pair_angle': 74.0,
                'pair_sep': 0.28,
                'equal_eta_mc_std': 28.0,
            })
        else:
            profile.update({
                'mc': 135.0,
                'radial': 18.0,
                'plane': 55.0,
                'axis': 14.0,
                'lateral': 12.0,
                'min_pair_angle': 60.0,
                'pair_sep': 0.18,
                'equal_eta_mc_std': 18.0,
            })

        if _is_f_block_like(metal_sym):
            profile['mc'] *= 0.80
            profile['radial'] *= 0.75
            profile['plane'] *= 0.45
            profile['axis'] *= 0.35
            profile['lateral'] *= 0.35
            profile['min_pair_angle'] = min(profile['min_pair_angle'], 72.0 if eta >= 5 else 58.0)
            profile['pair_sep'] *= 0.40
            profile['equal_eta_mc_std'] *= 0.70
        elif metal_sym in {'Fe', 'Co', 'Ni', 'Ru', 'Rh', 'Ir', 'Os', 'Mo', 'W', 'Re'} and eta >= 5:
            profile['plane'] *= 1.10
            profile['axis'] *= 1.20
            profile['lateral'] *= 1.15

        if n_groups_for_metal >= 2 and eta >= 4 and not _is_f_block_like(metal_sym):
            profile['axis'] *= 1.10
            profile['lateral'] *= 1.10
            profile['pair_sep'] *= 1.15

        return profile

    penalty = 0.0
    groups_by_metal: Dict[int, List[Tuple[int, List[int]]]] = {}
    for metal_idx, group_atoms in hapto_groups:
        groups_by_metal.setdefault(int(metal_idx), []).append((int(metal_idx), list(group_atoms)))

    centroid_records: Dict[int, List[Tuple[int, np.ndarray, np.ndarray, float, Dict[str, float]]]] = {}
    for metal_idx, group_atoms in hapto_groups:
        if not group_atoms:
            continue
        try:
            metal_idx = int(metal_idx)
            metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
            profile = _hapto_weight_profile(metal_sym, len(group_atoms), len(groups_by_metal.get(metal_idx, [])))
            metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)
            pts = np.array([conf.GetAtomPosition(int(atom_idx)) for atom_idx in group_atoms], dtype=float)
        except Exception:
            continue
        if pts.ndim != 2 or pts.shape[0] < 3:
            continue

        centroid = pts.mean(axis=0)
        mc_dist = float(np.linalg.norm(metal_pos - centroid))
        target_mc = float(_target_mc_dist(metal_sym, len(group_atoms)))
        penalty += profile['mc'] * (mc_dist - target_mc) ** 2

        mc_atom_dists = np.linalg.norm(pts - metal_pos, axis=1)
        if mc_atom_dists.size:
            penalty += profile['radial'] * float(np.std(mc_atom_dists) ** 2)

        q = pts - centroid
        try:
            _u, _s, vh = np.linalg.svd(q, full_matrices=False)
        except Exception:
            continue
        normal = vh[-1]
        n_norm = float(np.linalg.norm(normal))
        if n_norm < 1e-12:
            continue
        normal = normal / n_norm
        plane_dev = np.abs(q @ normal)
        plane_rms = float(np.sqrt(np.mean(plane_dev * plane_dev)))
        penalty += profile['plane'] * plane_rms * plane_rms

        axis = metal_pos - centroid
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm > 1e-12:
            axis = axis / axis_norm
            cosang = abs(float(np.dot(normal, axis)))
            penalty += profile['axis'] * (1.0 - cosang) ** 2
            lateral_offset = float(np.linalg.norm((metal_pos - centroid) - np.dot(metal_pos - centroid, normal) * normal))
            penalty += profile['lateral'] * lateral_offset * lateral_offset
            centroid_records.setdefault(metal_idx, []).append(
                (len(group_atoms), centroid, axis, mc_dist, profile)
            )

    for metal_idx, records in centroid_records.items():
        if len(records) < 2:
            continue
        for i in range(len(records)):
            eta_i, centroid_i, axis_i, mc_i, profile_i = records[i]
            for j in range(i + 1, len(records)):
                eta_j, centroid_j, axis_j, mc_j, profile_j = records[j]
                cosang = max(-1.0, min(1.0, float(np.dot(axis_i, axis_j))))
                angle = math.degrees(math.acos(cosang))
                min_pair_angle = min(profile_i['min_pair_angle'], profile_j['min_pair_angle'])
                if angle < min_pair_angle:
                    gap = (min_pair_angle - angle) / max(min_pair_angle, 1.0)
                    penalty += 220.0 * min(profile_i['pair_sep'], profile_j['pair_sep']) * gap * gap
                if eta_i == eta_j:
                    penalty += min(profile_i['equal_eta_mc_std'], profile_j['equal_eta_mc_std']) * (mc_i - mc_j) ** 2
    return penalty


def _hapto_candidate_quality_score(
    mol,
    conf_id: int = 0,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
) -> float:
    """Return a penalty score for a hapto candidate (lower = better)."""
    if not RDKIT_AVAILABLE or mol is None:
        return float("inf")

    score = 0.0
    try:
        score += float(_geometry_quality_score(mol, conf_id))
    except Exception:
        score += 1.0e6

    try:
        score += float(_hapto_geometry_quality_score(mol, conf_id, hapto_groups))
    except Exception:
        score += 1.0e5

    try:
        if _has_severe_covalent_distortion(mol, conf_id):
            score += 1.0e6
    except Exception:
        score += 1.0e5

    try:
        if _has_bad_geometry(mol, conf_id):
            score += 2.5e5
    except Exception:
        score += 1.0e5

    try:
        if _has_atom_clash(mol, conf_id, min_dist=0.80):
            score += 2.5e5
    except Exception:
        pass

    try:
        if _has_unphysical_metal_nonbonded_contact(mol, conf_id):
            score += 2.0e5
    except Exception:
        pass

    try:
        if _has_ligand_intertwining(mol, conf_id):
            score += 8.0e4
    except Exception:
        pass

    return score


def _hapto_mol_from_xyz_template(mol_template, xyz_delfin: str):
    """Map DELFIN XYZ coordinates back onto a template molecule."""
    if not RDKIT_AVAILABLE or mol_template is None:
        return None
    try:
        mol_tmp = Chem.Mol(mol_template)
        mol_tmp.RemoveAllConformers()
        conf = _xyz_to_rdkit_conformer(mol_tmp, xyz_delfin)
        if conf is None:
            return None
        mol_tmp.AddConformer(conf, assignId=True)
        return mol_tmp
    except Exception:
        return None


def _refine_hapto_candidate_with_rdkit_uff(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
):
    """Run local RDKit UFF while freezing the hapto scaffold."""
    if not RDKIT_AVAILABLE or mol is None:
        return None
    try:
        work = Chem.Mol(mol)
        if work.GetNumConformers() == 0:
            return None
        ff = AllChem.UFFGetMoleculeForceField(work, confId=0)
        if ff is None:
            return None

        fixed: set = set()
        for atom in work.GetAtoms():
            if atom.GetSymbol() in _METAL_SET:
                fixed.add(atom.GetIdx())
        for metal_idx, catoms in hapto_groups:
            fixed.add(metal_idx)
            fixed.update(catoms)
        for atom_idx in sorted(fixed):
            ff.AddFixedPoint(int(atom_idx))

        ff.Minimize(maxIts=600)
        try:
            _enforce_donor_pi_coplanarity(work, 0, hapto_groups)
        except Exception:
            pass
        try:
            _fix_secondary_metal_distances(work, 0, hapto_groups)
        except Exception:
            pass
        try:
            _final_clash_resolution(work, 0, hapto_groups)
        except Exception:
            pass
        return work
    except Exception:
        return None


def _select_best_hapto_candidate(
    smiles: str,
    hapto_groups: List[Tuple[int, List[int]]],
    candidates: List[Tuple[str, object]],
    *,
    apply_uff: bool,
):
    """Select the best topology-preserving hapto candidate."""
    accepted: List[Tuple[float, object, str]] = []
    relaxed: List[Tuple[float, object, str]] = []

    for label, candidate_mol in candidates:
        if candidate_mol is None or candidate_mol.GetNumConformers() == 0:
            continue

        try:
            xyz_candidate = _mol_to_xyz(candidate_mol)
        except Exception:
            continue

        topo_ok = _hapto_candidate_topology_ok(
            xyz_candidate,
            smiles,
            mol=candidate_mol,
            conf_id=0,
            hapto_groups=hapto_groups,
        )
        base_score = _hapto_candidate_quality_score(candidate_mol, 0, hapto_groups)
        target_bucket = accepted if topo_ok else relaxed
        target_bucket.append((base_score, candidate_mol, label))

        if not apply_uff:
            continue

        refined_trials: List[Tuple[str, object]] = []
        rdkit_refined = _refine_hapto_candidate_with_rdkit_uff(candidate_mol, hapto_groups)
        if rdkit_refined is not None:
            refined_trials.append((f"{label}+rdkit-uff", rdkit_refined))

        try:
            xyz_refined = _optimize_xyz_openbabel_safe(
                xyz_candidate,
                mol_template=candidate_mol,
                smiles=smiles,
                steps=750,
                apply_template_constraints=True,
            )
        except Exception:
            xyz_refined = xyz_candidate
        if xyz_refined and xyz_refined != xyz_candidate:
            ob_refined = _hapto_mol_from_xyz_template(candidate_mol, xyz_refined)
            if ob_refined is not None:
                refined_trials.append((f"{label}+ob-uff", ob_refined))

        for refined_label, refined_mol in refined_trials:
            try:
                xyz_refined = _mol_to_xyz(refined_mol)
            except Exception:
                continue
            if not _hapto_candidate_topology_ok(
                xyz_refined,
                smiles,
                mol=refined_mol,
                conf_id=0,
                hapto_groups=hapto_groups,
            ):
                continue
            refined_score = _hapto_candidate_quality_score(refined_mol, 0, hapto_groups)
            if refined_score + 1e-6 < base_score:
                accepted.append((refined_score, refined_mol, refined_label))

    pool = accepted if accepted else relaxed
    if not pool:
        return None
    pool.sort(key=lambda item: (item[0], item[2]))
    best_score, best_mol, best_label = pool[0]
    logger.info(
        "Selected hapto candidate %s (score=%.2f, strict_topology=%s, pool=%d)",
        best_label,
        best_score,
        bool(accepted),
        len(pool),
    )
    _enforce_metal_topology(best_mol, conf_id=0)
    return best_mol


def _enforce_metal_topology(mol, conf_id: int = 0, min_nonbonded: float = 2.5):
    """Push non-bonded atoms away from metals to preserve SMILES topology.

    Operates as a final post-processing step: for every metal, any atom
    NOT bonded to it in the molecular graph but closer than *min_nonbonded*
    gets radially pushed outward.  Hapto ring atoms are moved as a rigid
    group to preserve ring geometry.
    """
    if not RDKIT_AVAILABLE or mol is None:
        return
    try:
        import numpy as np
        conf = mol.GetConformer(conf_id)
    except Exception:
        return

    def _gp(i):
        p = conf.GetAtomPosition(i)
        return np.array([p.x, p.y, p.z])

    def _sp(i, arr):
        conf.SetAtomPosition(i, Point3D(float(arr[0]), float(arr[1]),
                                         float(arr[2])))

    n = mol.GetNumAtoms()
    metal_indices = [i for i in range(n)
                     if mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET]
    if not metal_indices:
        return

    metal_bonded: Dict[int, set] = {}
    for mi in metal_indices:
        bonded = set()
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            bonded.add(nbr.GetIdx())
        metal_bonded[mi] = bonded

    hapto_groups = _find_hapto_groups(mol)
    hapto_of: Dict[int, int] = {}
    for gi, (_, grp) in enumerate(hapto_groups):
        for a in grp:
            hapto_of[a] = gi
    group_atoms: Dict[int, List[int]] = {}
    for gi, (_, grp) in enumerate(hapto_groups):
        group_atoms[gi] = list(grp)

    for _pass in range(10):
        any_moved = False
        for mi in metal_indices:
            mpos = _gp(mi)
            bonded = metal_bonded[mi]
            for ai in range(n):
                if ai == mi or ai in bonded:
                    continue
                if mol.GetAtomWithIdx(ai).GetSymbol() == 'H':
                    continue
                if mol.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET:
                    continue
                apos = _gp(ai)
                d = float(np.linalg.norm(apos - mpos))
                if d >= min_nonbonded or d < 1e-8:
                    continue
                push_dir = (apos - mpos) / d
                push_dist = (min_nonbonded - d) * 1.1
                gi = hapto_of.get(ai)
                if gi is not None:
                    for ha in group_atoms[gi]:
                        _sp(ha, _gp(ha) + push_dir * push_dist)
                else:
                    _sp(ai, apos + push_dir * push_dist)
                any_moved = True
        if not any_moved:
            break
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

        # --- 9-coordinate patterns ---
        if n_coord == 9:
            if len(elem_counts) == 1:
                return ''
            return ''

    return ''


# ---------------------------------------------------------------------------
# Topological isomer enumerator (Feature 1)
# ---------------------------------------------------------------------------

def _chelate_pairs(mol, metal_idx: int, donor_indices: List[int],
                    max_path: int = 5) -> List[FrozenSet]:
    """Find donor pairs connected through a short non-metal path (chelate constraints).

    BFS from each donor atom to every other donor, blocking the metal.
    A successful path with length <= *max_path* bonds means the two donors
    form a chelate ring small enough to enforce a *cis* constraint (up to a
    ``max_path + 1``-membered chelate ring counting the metal).

    Longer paths (e.g. opposite donors of a 14-membered macrocycle) are
    **not** marked as chelate because they can adopt trans arrangements.

    Args:
        max_path: Maximum number of bonds in the non-metal path between two
            donors to count as a chelate pair.  Default 5 corresponds to a
            7-membered chelate ring (donor–5 bridge atoms–donor + metal),
            capturing common diphosphine (dppp) and diamine chelates.

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


def _canonical_ss(types: tuple) -> tuple:
    """Canonical form for see-saw / C2v (CN=4).

    Positions 0,1 = axial (trans pair), 2,3 = equatorial (cis pair).
    """
    axial = tuple(sorted([types[0], types[1]]))
    equatorial = tuple(sorted([types[2], types[3]]))
    return ('SS', axial, equatorial)


def _canonical_tpr(types: tuple) -> tuple:
    """Canonical form for trigonal-prismatic (D3h, CN=6).

    Positions 0-2 = top triangle, 3-5 = bottom triangle (eclipsed).
    No 180° trans pairs exist in a trigonal prism.
    """
    top = tuple(sorted(types[:3]))
    bottom = tuple(sorted(types[3:6]))
    combined = tuple(sorted([top, bottom]))
    return ('TPR', combined)


def _canonical_coh(types: tuple) -> tuple:
    """Canonical form for capped octahedral (C3v, CN=7).

    Positions 0-5 = octahedral base (3 trans pairs: 0-1, 2-3, 4-5),
    position 6 = capping atom above one triangular face.
    """
    base_pairs = tuple(sorted([
        tuple(sorted([types[0], types[1]])),
        tuple(sorted([types[2], types[3]])),
        tuple(sorted([types[4], types[5]])),
    ]))
    return ('COH', base_pairs, types[6])


def _canonical_ttp(types: tuple) -> tuple:
    """Canonical form for tricapped trigonal-prismatic (D3h, CN=9).

    Positions 0-5 = prism vertices (top 0-2, bottom 3-5 eclipsed),
    positions 6-8 = equatorial caps.
    """
    prism = tuple(sorted(types[:6]))
    caps = tuple(sorted(types[6:9]))
    return ('TTP', prism, caps)


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
    'SS':  [(0, 1)],          # see-saw: axial pair only
    'TPR': [],                # trigonal-prismatic: no 180° pairs
    'COH': [(0, 1), (2, 3), (4, 5)],  # capped-oh: octahedral base trans pairs
    'TTP': [],                # tricapped trigonal prism: no 180° pairs
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
    'SS':  _canonical_ss,
    'TPR': _canonical_tpr,
    'COH': _canonical_coh,
    'TTP': _canonical_ttp,
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
    # See-saw / C2v (CN=4): axial pair along z, equatorial pair in xz-plane
    'SS':  [(0, 0, 2), (0, 0, -2), (2, 0, 0.4), (-2, 0, 0.4)],
    # Trigonal-prismatic (D3h, CN=6): eclipsed top/bottom triangles
    'TPR': [(2, 0, 1), (-1, 1.732, 1), (-1, -1.732, 1),
            (2, 0, -1), (-1, 1.732, -1), (-1, -1.732, -1)],
    # Capped octahedral (C3v, CN=7): octahedral base + capping 7th
    'COH': [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2),
            (1.155, 1.155, 1.155)],
    # Tricapped trigonal prism (D3h, CN=9): 6 prism + 3 equatorial caps
    'TTP': [(1.633, 0, 1.155), (-0.816, 1.414, 1.155), (-0.816, -1.414, 1.155),
            (1.633, 0, -1.155), (-0.816, 1.414, -1.155), (-0.816, -1.414, -1.155),
            (1.0, 1.732, 0), (-2.0, 0, 0), (1.0, -1.732, 0)],
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
        geometries = [pref, other, 'SS']
    elif n_coord == 5:
        pref5 = _PREFERRED_CN5_GEOMETRY.get(metal_symbol, 'TBP')
        other5 = 'SP' if pref5 == 'TBP' else 'TBP'
        geometries = [pref5, other5]
    elif n_coord == 6:
        pref6 = _PREFERRED_CN6_GEOMETRY.get(metal_symbol, 'OH')
        other6 = 'TPR' if pref6 == 'OH' else 'OH'
        geometries = [pref6, other6]
    elif n_coord == 7:
        geometries = ['PBP', 'COH']
    elif n_coord == 8:
        geometries = ['SAP', 'DD']
    elif n_coord == 9:
        geometries = ['TTP']
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
    try:
        cid = AllChem.EmbedMolecule(frag_mol, params)
    except Exception:
        cid = -1
    if cid < 0:
        # Fallback for problematic aromatic/charged fragments where ETKDG
        # fails in kekulization preprocessing. We only need a geometric
        # scaffold here, so a temporarily dearomatized embedding is fine.
        try:
            rw2 = Chem.RWMol(frag_mol)
            for atom in rw2.GetAtoms():
                if atom.GetIsAromatic():
                    atom.SetIsAromatic(False)
            for bond in rw2.GetBonds():
                if bond.GetIsAromatic() or bond.GetBondType() == Chem.BondType.AROMATIC:
                    bond.SetIsAromatic(False)
                    bond.SetBondType(Chem.BondType.SINGLE)
            frag_mol2 = rw2.GetMol()
            try:
                frag_mol2.UpdatePropertyCache(strict=False)
            except Exception:
                pass
            cid2 = AllChem.EmbedMolecule(frag_mol2, params)
            if cid2 < 0:
                return False
            frag_mol = frag_mol2
            cid = cid2
        except Exception:
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


def _build_multimetal_scaffold(
    mol,
    metal_indices: List[int],
    bridging_donors: List[Tuple[int, List[int]]],
) -> Optional[Dict[int, Tuple[float, float, float]]]:
    """Place metals + bridging donors as a rigid scaffold.

    Strategy:
    1. Metal1 at origin
    2. Each bridging donor at ideal M1-D distance along geometry vectors
    3. Metal2 placed so M2-bridge distances match ideals
    4. Returns {atom_idx: (x,y,z)} for all scaffold atoms

    This ensures the M-bridge-M core has correct geometry BEFORE
    peripheral donors are placed.
    """
    if len(metal_indices) < 2 or not bridging_donors:
        return None

    try:
        import numpy as np

        coords: Dict[int, Tuple[float, float, float]] = {}
        m1 = metal_indices[0]
        m2 = metal_indices[1]
        m1_sym = mol.GetAtomWithIdx(m1).GetSymbol()
        m2_sym = mol.GetAtomWithIdx(m2).GetSymbol()

        # Place M1 at origin.
        coords[m1] = (0.0, 0.0, 0.0)

        # Identify bridging donor atoms shared between M1 and M2.
        shared_bridges = []
        for d_idx, m_list in bridging_donors:
            if m1 in m_list and m2 in m_list:
                shared_bridges.append(d_idx)

        if not shared_bridges:
            return None

        # Place first bridging donor on +x axis at ideal M1-D distance.
        d0 = shared_bridges[0]
        d0_sym = mol.GetAtomWithIdx(d0).GetSymbol()
        d_m1 = float(_get_ml_bond_length(m1_sym, d0_sym))
        d_m2 = float(_get_ml_bond_length(m2_sym, d0_sym))
        coords[d0] = (d_m1, 0.0, 0.0)

        # Place M2: on the M1-D0-M2 plane.
        # M2 is at distance d_m2 from D0.
        # M1-D0 = d_m1 along x. M2 must be at d_m2 from D0.
        # Choose M1-D0-M2 angle ~130° (typical for μ-O bridge).
        bridge_angle_rad = math.radians(130.0)
        m2_x = d_m1 + d_m2 * math.cos(math.pi - bridge_angle_rad)
        m2_y = d_m2 * math.sin(math.pi - bridge_angle_rad)
        coords[m2] = (m2_x, m2_y, 0.0)

        # Place additional bridging donors (if any) between the two metals.
        for k, dk in enumerate(shared_bridges[1:], 1):
            dk_sym = mol.GetAtomWithIdx(dk).GetSymbol()
            dk_m1 = float(_get_ml_bond_length(m1_sym, dk_sym))
            dk_m2 = float(_get_ml_bond_length(m2_sym, dk_sym))
            # Place below the M1-M2 plane (alternating above/below)
            angle_offset = math.radians(130.0)
            sign = -1.0 if k % 2 == 1 else 1.0
            mid_x = (coords[m1][0] + coords[m2][0]) / 2.0
            mid_y = (coords[m1][1] + coords[m2][1]) / 2.0
            coords[dk] = (mid_x, mid_y, sign * dk_m1 * 0.8)

        return coords
    except Exception as exc:
        logger.debug("_build_multimetal_scaffold failed: %s", exc)
        return None


def _build_topology_xyz(
    mol,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    apply_uff: bool,
    conf_id: Optional[int] = None,
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
        conf_id: Optional template conformer ID. If None, the rigid-fragment
            builder auto-selects the best-scoring conformer.  Callers iterating
            over multiple alternative binding modes can use this to avoid
            depending on a single (possibly poor) template.

    Returns:
        DELFIN-format XYZ string, or None on failure.
    """
    try:
        # If a template conformer is available, place whole ligand fragments
        # as rigid bodies onto the target donor geometry. This preserves
        # intraligand structure much better than de-novo fragment embedding.
        if mol.GetNumConformers() > 0:
            xyz_from_template = _build_topology_xyz_from_template(
                mol, metal_idx, donor_atom_indices, perm, geometry, apply_uff,
                conf_id=conf_id,
            )
            if xyz_from_template is not None:
                return xyz_from_template

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
                try:
                    ok = _embed_fragment_procrustes(
                        mol, metal_idx, frag, frag_donors, tgt_positions, coords
                    )
                except Exception as frag_exc:
                    logger.debug(
                        "Fragment embedding failed for one fragment (size=%d, donors=%d): %s",
                        len(frag), len(frag_donors), frag_exc,
                    )
                    continue
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
                # Coordination-preserving UFF: M-D distances and L-M-L
                # angles are pinned to the idealized polyhedron so UFF
                # cannot distort the octahedral/PBP/etc. cage even when
                # it lacks force-field parameters for the metal.
                coord_constraints = None
                try:
                    coord_constraints = _build_coordination_constraints_from_xyz(
                        mol, xyz,
                    )
                except Exception as cexc:
                    logger.debug(
                        "Coordination constraint build failed, falling back to template constraints: %s",
                        cexc,
                    )
                xyz = _optimize_xyz_openbabel_safe(
                    xyz,
                    mol_template=mol,
                    coord_constraints=coord_constraints,
                )
            except Exception as uff_exc:
                # Keep the generated topology geometry when UFF fails.
                # Dropping the isomer here can hide valid alternatives.
                logger.debug("Topology UFF optimization failed, keeping unoptimized XYZ: %s", uff_exc)

        return xyz
    except Exception as e:
        logger.debug("_build_topology_xyz failed: %s", e)
        return None


def _rank_template_conformers(mol, *, top_k: Optional[int] = None) -> List[int]:
    """Return conformer IDs ranked deterministically by geometry quality.

    Lower ``_geometry_quality_score`` is better.  Conformers with an atom clash
    below ``min_dist=0.3`` are excluded (truly collapsed).  Ties are broken by
    ascending conformer ID so the ordering is fully reproducible.

    If ``top_k`` is given, only that many best candidates are returned.
    """
    try:
        all_ids = [int(c.GetId()) for c in mol.GetConformers()]
    except Exception:
        return []
    if not all_ids:
        return []

    ranked: List[Tuple[float, int]] = []
    for cid in all_ids:
        try:
            if _has_atom_clash(mol, cid, min_dist=0.3):
                continue
            score = float(_geometry_quality_score(mol, cid))
        except Exception:
            continue
        ranked.append((score, cid))

    ranked.sort(key=lambda pair: (pair[0], pair[1]))
    ordered = [cid for _score, cid in ranked]
    if top_k is not None and top_k >= 0:
        ordered = ordered[:top_k]
    return ordered


def _build_topology_xyz_from_template(
    mol,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    apply_uff: bool,
    conf_id: Optional[int] = None,
) -> Optional[str]:
    """Rigid-fragment topology builder using an existing template conformer.

    The ligand fragments are transformed as rigid bodies so their internal
    geometry remains close to the template. This is especially useful for
    aromatic/charged chelating ligands where ETKDG fragment embedding can fail.

    If ``conf_id`` is None, the best-scored conformer (per
    :func:`_rank_template_conformers`) is used.  Pass a specific ID when the
    caller wants to iterate over several candidates (e.g. to try multiple
    templates for the same alternative binding mode).
    """
    if not RDKIT_AVAILABLE:
        return None
    if mol.GetNumConformers() == 0:
        return None

    try:
        import numpy as np
    except Exception:
        return None

    if conf_id is None:
        ranked = _rank_template_conformers(mol, top_k=1)
        if not ranked:
            return None
        conf_id = ranked[0]

    try:
        conf = mol.GetConformer(int(conf_id))
        vectors = _TOPO_GEOMETRY_VECTORS[geometry]
        n_atoms = mol.GetNumAtoms()

        # Original coordinates translated so the metal sits at origin.
        mpos = conf.GetAtomPosition(metal_idx)
        orig = np.zeros((n_atoms, 3), dtype=float)
        for i in range(n_atoms):
            p = conf.GetAtomPosition(i)
            orig[i, 0] = p.x - mpos.x
            orig[i, 1] = p.y - mpos.y
            orig[i, 2] = p.z - mpos.z

        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        donor_target_map: Dict[int, np.ndarray] = {}
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
            donor_target_map[donor_atom_idx] = np.array([vx, vy, vz], dtype=float)

        # Build non-metal fragments.
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

        coords = np.array(orig, copy=True)
        coords[metal_idx] = np.array([0.0, 0.0, 0.0], dtype=float)
        placed = {metal_idx}

        def _rotation_from_vectors(v_from: np.ndarray, v_to: np.ndarray) -> np.ndarray:
            nf = np.linalg.norm(v_from)
            nt = np.linalg.norm(v_to)
            if nf < 1e-10 or nt < 1e-10:
                return np.eye(3)
            a = v_from / nf
            b = v_to / nt
            v = np.cross(a, b)
            s = np.linalg.norm(v)
            c = float(np.clip(np.dot(a, b), -1.0, 1.0))
            if s < 1e-10:
                if c > 0:
                    return np.eye(3)
                # 180° rotation around any axis perpendicular to a
                axis = np.array([1.0, 0.0, 0.0])
                if abs(a[0]) > 0.9:
                    axis = np.array([0.0, 1.0, 0.0])
                axis = axis - np.dot(axis, a) * a
                axis = axis / max(np.linalg.norm(axis), 1e-10)
                K = np.array([
                    [0, -axis[2], axis[1]],
                    [axis[2], 0, -axis[0]],
                    [-axis[1], axis[0], 0],
                ])
                return np.eye(3) + 2.0 * (K @ K)
            K = np.array([
                [0, -v[2], v[1]],
                [v[2], 0, -v[0]],
                [-v[1], v[0], 0],
            ])
            return np.eye(3) + K + K @ K * ((1.0 - c) / (s ** 2))

        def _metal_proximity_penalty(
            xyz_frag: np.ndarray,
            frag_atoms: List[int],
            donor_atoms: List[int],
        ) -> float:
            """Penalty for non-donor heavy atoms placed too close to the metal."""
            donor_set = set(donor_atoms)
            pen = 0.0
            for li, atom_idx in enumerate(frag_atoms):
                if atom_idx in donor_set:
                    continue
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() <= 1:
                    continue
                if atom.GetSymbol() in _METAL_SET:
                    continue
                d = float(np.linalg.norm(xyz_frag[li]))
                sym = atom.GetSymbol()
                # Keep non-donor atoms clearly outside the first coordination shell.
                ml_ref = float(_get_ml_bond_length(metal_sym, sym))
                min_allowed = max(1.15, 0.65 * ml_ref)
                if d < min_allowed:
                    dd = (min_allowed - d)
                    pen += dd * dd
                if d < 1.0:
                    pen += 5.0
            return pen

        for frag in fragments:
            frag_list = sorted(frag)
            frag_donors = [d for d in donor_atom_indices if d in frag]
            if not frag_donors:
                continue

            frag_xyz = orig[frag_list, :]
            donor_local = [frag_list.index(d) for d in frag_donors]
            src = frag_xyz[donor_local, :]
            tgt = np.array([donor_target_map[d] for d in frag_donors], dtype=float)
            if len(src) != len(tgt) or len(src) == 0:
                continue

            if len(src) >= 2:
                src_center = src.mean(axis=0)
                tgt_center = tgt.mean(axis=0)
                X = src - src_center
                Y = tgt - tgt_center
                H = X.T @ Y
                U, S, Vt = np.linalg.svd(H)
                R = Vt.T @ U.T
                if np.linalg.det(R) < 0:
                    Vt[-1, :] *= -1
                    R = Vt.T @ U.T
                transformed = (frag_xyz - src_center) @ R.T + tgt_center

                # Bidentate/multidentate fragments can be mirrored around the
                # donor-donor axis, yielding two plausible orientations with
                # identical donor placement. Pick the orientation that keeps
                # non-donor atoms farther from the metal center.
                try:
                    if len(frag_donors) >= 2:
                        d0 = donor_target_map[frag_donors[0]]
                        d1 = donor_target_map[frag_donors[1]]
                        axis = d1 - d0
                        axis_norm = float(np.linalg.norm(axis))
                        if axis_norm > 1e-10:
                            u = axis / axis_norm
                            pivot = 0.5 * (d0 + d1)
                            v = transformed - pivot
                            # 180° rotation around donor-donor axis:
                            # v' = -v + 2*(u·v)*u
                            v_rot = -v + 2.0 * np.outer(v @ u, u)
                            transformed_flip = v_rot + pivot

                            p0 = _metal_proximity_penalty(transformed, frag_list, frag_donors)
                            p1 = _metal_proximity_penalty(transformed_flip, frag_list, frag_donors)
                            if p1 < p0:
                                transformed = transformed_flip
                except Exception:
                    pass
            else:
                src_d = src[0]
                tgt_d = tgt[0]
                src_com = frag_xyz.mean(axis=0)
                # Keep the fragment extending away from the metal.
                v_from = src_com - src_d
                v_to = tgt_d
                R = _rotation_from_vectors(v_from, v_to)
                transformed = (frag_xyz - src_d) @ R.T + tgt_d

            for li, atom_idx in enumerate(frag_list):
                coords[atom_idx] = transformed[li]
                placed.add(atom_idx)

        # Any atom not touched by fragment placement keeps template-relative coords.
        for i in range(n_atoms):
            if i not in placed:
                coords[i] = orig[i]

        lines = []
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            x, y, z = coords[i]
            lines.append(f"{atom.GetSymbol():4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}")
        xyz = '\n'.join(lines) + '\n'

        # Conservative UFF: only keep UFF result if it preserves topology.
        # For metals without UFF parameters, UFF can BREAK the structure.
        # The pre-UFF Procrustes geometry has correct M-D distances and
        # is often better than what UFF produces.
        if apply_uff:
            xyz_pre_uff = xyz
            try:
                coord_constraints = None
                try:
                    coord_constraints = _build_coordination_constraints_from_xyz(
                        mol, xyz,
                    )
                except Exception:
                    pass
                xyz_uff = _optimize_xyz_openbabel_safe(
                    xyz,
                    mol_template=mol,
                    coord_constraints=coord_constraints,
                )
                # Keep UFF only if topology survives.
                if _verify_topology_from_graph(xyz_uff, mol):
                    xyz = xyz_uff
                else:
                    logger.debug(
                        "UFF broke topology in template builder — keeping pre-UFF geometry"
                    )
                    xyz = xyz_pre_uff
            except Exception as uff_exc:
                logger.debug(
                    "Template-topology UFF failed, keeping pre-UFF: %s", uff_exc
                )
        return xyz
    except Exception as e:
        logger.debug("_build_topology_xyz_from_template failed: %s", e)
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

    def _passes_chelate_distance_feasibility(
        _mol,
        _metal_idx: int,
        _donor_indices: List[int],
        _perm: List[int],
        _geom_name: str,
        _chelate_atom_pairs: List[FrozenSet],
        abs_tol: float = 0.7,
        rel_tol: float = 0.35,
    ) -> bool:
        """Reject geometrically impossible chelate placements.

        For each chelate pair, compare donor-donor distance in the template
        conformer to the idealized target distance implied by geometry+perm.
        If they differ too much, this arrangement is likely non-physical for
        the ligand bite and tends to collapse into unrealistic structures.
        """
        try:
            if _mol.GetNumConformers() == 0:
                return True
            conf = _mol.GetConformer(0)
            vectors = _TOPO_GEOMETRY_VECTORS.get(_geom_name)
            if not vectors:
                return True

            m_sym = _mol.GetAtomWithIdx(_metal_idx).GetSymbol()
            target_by_donor: Dict[int, Tuple[float, float, float]] = {}
            for pos_idx, donor_list_idx in enumerate(_perm):
                d_atom = _donor_indices[donor_list_idx]
                d_sym = _mol.GetAtomWithIdx(d_atom).GetSymbol()
                bl = _get_ml_bond_length(m_sym, d_sym)
                vx, vy, vz = vectors[pos_idx]
                mag = math.sqrt(vx * vx + vy * vy + vz * vz)
                if mag > 1e-8:
                    vx = vx / mag * bl
                    vy = vy / mag * bl
                    vz = vz / mag * bl
                target_by_donor[d_atom] = (vx, vy, vz)

            for cp in _chelate_atom_pairs:
                pair = sorted(cp)
                if len(pair) != 2:
                    continue
                a, b = pair
                if a not in target_by_donor or b not in target_by_donor:
                    continue

                pa = conf.GetAtomPosition(a)
                pb = conf.GetAtomPosition(b)
                d_src = math.sqrt(
                    (pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2
                )
                ta = target_by_donor[a]
                tb = target_by_donor[b]
                d_tgt = math.sqrt(
                    (ta[0] - tb[0]) ** 2 + (ta[1] - tb[1]) ** 2 + (ta[2] - tb[2]) ** 2
                )

                tol = max(abs_tol, rel_tol * max(d_src, 1e-8))
                if abs(d_tgt - d_src) > tol:
                    return False
            return True
        except Exception:
            return True

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        donor_indices = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        n_coord = len(donor_indices)

        if n_coord < 2 or n_coord > 9:
            continue

        # Use donor environment classes (element + Morgan environment) instead
        # of plain element symbols. This preserves chemically meaningful
        # distinctions like aqua-O vs carboxylate-O, which are essential for
        # CN=7/8 trans-pattern completeness.
        donor_keys = [
            dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset()))
            for d in donor_indices
        ]
        uniq_keys = sorted(
            set(donor_keys),
            key=lambda k: (k[0], tuple(sorted(k[1]))),
        )
        key_to_class = {k: i for i, k in enumerate(uniq_keys)}
        donor_labels = [f"{k[0]}{key_to_class[k]}" for k in donor_keys]

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

        # Chelate-distance feasibility is a useful guard, but can over-prune
        # higher-coordination systems (notably CN=7) when idealized vectors and
        # template distances differ systematically. If it rejects everything,
        # fall back to the unfiltered topological set.
        feasible_isomers: List[Tuple[tuple, List[int]]] = []
        for canonical_form, perm in isomers:
            geom_name = canonical_form[0]
            if _passes_chelate_distance_feasibility(
                mol, metal_idx, donor_indices, perm, geom_name, chelate_ps
            ):
                feasible_isomers.append((canonical_form, perm))
        if not feasible_isomers and isomers:
            logger.debug(
                "Chelate-distance feasibility rejected all %d topo isomer(s) "
                "for CN=%d; using unfiltered set.",
                len(isomers), n_coord,
            )
            feasible_isomers = isomers

        # Pre-compute ranked template conformers once per metal centre so each
        # permutation can retry against several templates when the default
        # (best-scored) one produces no viable XYZ.
        topo_template_cids = _rank_template_conformers(mol, top_k=3) or [None]

        _PRIMARY_GEOM_BASE = {
            2: 'LIN', 3: 'TP', 5: 'TBP',
            6: 'OH', 7: 'PBP', 8: 'SAP', 9: 'TTP',
        }
        _PRIMARY_GEOM = dict(_PRIMARY_GEOM_BASE)
        _PRIMARY_GEOM[4] = _PREFERRED_CN4_GEOMETRY.get(
            atom.GetSymbol(), 'SQ'
        )
        _GEOM_PRETTY = {
            'LIN': 'linear', 'TP': 'trigonal-planar', 'TS': 'T-shaped',
            'SQ': 'square-planar', 'TH': 'tetrahedral', 'SS': 'see-saw',
            'TBP': 'trigonal-bipyramidal', 'SP': 'square-pyramidal',
            'OH': 'octahedral', 'TPR': 'trigonal-prismatic',
            'PBP': 'pentagonal-bipyramidal', 'COH': 'capped-octahedral',
            'SAP': 'square-antiprismatic', 'DD': 'dodecahedral',
            'TTP': 'tricapped-trigonal-prismatic',
        }
        _primary_geom = _PRIMARY_GEOM.get(n_coord)

        def _build_one_topo(args):
            cf, pm = args
            gn = cf[0]
            try:
                xyz = None
                for _tc in topo_template_cids:
                    xyz = _build_topology_xyz(
                        mol, metal_idx, donor_indices, pm, gn,
                        apply_uff, conf_id=_tc,
                    )
                    if xyz is not None:
                        break
                if xyz is None:
                    return None
                mt = Chem.RWMol(mol)
                mt.RemoveAllConformers()
                c = _xyz_to_rdkit_conformer(mt.GetMol(), xyz)
                if c is None:
                    return None
                ci = mt.AddConformer(c, assignId=True)
                try:
                    if _has_atom_clash(mt.GetMol(), ci, min_dist=0.3):
                        return None
                    if _has_unphysical_metal_nonbonded_contact(mt.GetMol(), ci):
                        return None
                    if _has_unphysical_oco_geometry(mt.GetMol(), ci):
                        return None
                    if _has_pi_ring_nonplanarity(mt.GetMol(), ci):
                        return None
                    if _has_severe_covalent_distortion(mt.GetMol(), ci):
                        return None
                except Exception:
                    pass
                fp = _compute_coordination_fingerprint(
                    mt.GetMol(), ci, dtype_map=dtype_map
                )
                lbl = _classify_isomer_label(fp, mt.GetMol())
                if _primary_geom and gn != _primary_geom:
                    gp = _GEOM_PRETTY.get(gn, gn)
                    lbl = f'{gp} {lbl}' if lbl else gp
                return (xyz, lbl)
            except Exception as exc:
                logger.debug("Topo isomer build failed (%s): %s", gn, exc)
                return None

        # Build 3D for each permutation. OB UFF holds the GIL so
        # ProcessPoolExecutor is needed for real parallelism.
        # Strategy: build pre-UFF Procrustes XYZ in main process (fast),
        # batch all UFF calls through ProcessPool, then quality-check.
        _pre_uff_batch: List[Tuple[tuple, List[int], str, str, Optional[Dict]]] = []
        for cf, pm in feasible_isomers:
            if len(_pre_uff_batch) + len(results) >= max_isomers * 3:
                break
            gn = cf[0]
            try:
                xyz = None
                for _tc in topo_template_cids:
                    xyz = _build_topology_xyz(
                        mol, metal_idx, donor_indices, pm, gn,
                        False, conf_id=_tc,
                    )
                    if xyz is not None:
                        break
                if xyz is None:
                    continue
                # No pre-UFF clash gate — Procrustes can have transient
                # overlaps that UFF resolves. Graph check after UFF.
                # Build constraints for this permutation.
                coord_c = None
                if apply_uff:
                    try:
                        coord_c = _build_coordination_constraints_from_xyz(
                            mol, xyz,
                        )
                    except Exception:
                        pass
                _pre_uff_batch.append((cf, pm, gn, xyz, coord_c))
            except Exception as exc:
                logger.debug("Topo pre-UFF build failed (%s): %s", cf[0], exc)

        # Batch UFF via ProcessPool (OB holds GIL → threads don't help).
        if apply_uff and _pre_uff_batch:
            _n_uff_workers = min(len(_pre_uff_batch), os.cpu_count() or 4, 32)
            _uff_inputs = [
                (xyz, 500, cstr) for _cf, _pm, _gn, xyz, cstr in _pre_uff_batch
            ]
            try:
                if _n_uff_workers > 1 and len(_uff_inputs) > 2:
                    with concurrent.futures.ProcessPoolExecutor(
                        max_workers=_n_uff_workers
                    ) as _pp:
                        _uff_results = list(_pp.map(
                            _optimize_xyz_openbabel,
                            [inp[0] for inp in _uff_inputs],
                            [inp[1] for inp in _uff_inputs],
                            [inp[2] for inp in _uff_inputs],
                        ))
                else:
                    _uff_results = [
                        _optimize_xyz_openbabel(inp[0], inp[1], inp[2])
                        for inp in _uff_inputs
                    ]
            except Exception as _ppe:
                logger.debug("ProcessPool UFF failed, falling back to sequential: %s", _ppe)
                _uff_results = [
                    _optimize_xyz_openbabel(inp[0], inp[1], inp[2])
                    for inp in _uff_inputs
                ]
            for idx, (cf, pm, gn, xyz_pre, _cstr) in enumerate(_pre_uff_batch):
                xyz_opt = _uff_results[idx] if idx < len(_uff_results) else xyz_pre
                if not xyz_opt:
                    xyz_opt = xyz_pre
                # Conservative UFF: if UFF broke topology, keep pre-UFF.
                if not _verify_topology_from_graph(xyz_opt, mol):
                    xyz_opt = xyz_pre
                _pre_uff_batch[idx] = (cf, pm, gn, xyz_opt, _cstr)

        # Post-UFF: graph-based topology check (replaces the 5 legacy
        # checks that were too aggressive for topo-generated structures).
        for cf, pm, gn, xyz, _cstr in _pre_uff_batch:
            if len(results) >= max_isomers:
                break
            try:
                if not _verify_topology_from_graph(xyz, mol):
                    continue
                mt = Chem.RWMol(mol)
                mt.RemoveAllConformers()
                c = _xyz_to_rdkit_conformer(mt.GetMol(), xyz)
                if c is None:
                    continue
                ci = mt.AddConformer(c, assignId=True)
                fp = _compute_coordination_fingerprint(
                    mt.GetMol(), ci, dtype_map=dtype_map
                )
                lbl = _classify_isomer_label(fp, mt.GetMol())
                if _primary_geom and gn != _primary_geom:
                    gp = _GEOM_PRETTY.get(gn, gn)
                    lbl = f'{gp} {lbl}' if lbl else gp
                results.append((xyz, lbl))
            except Exception as exc:
                logger.debug("Topo post-UFF check failed (%s): %s", gn, exc)
                continue

    # --- Multinuclear coupled enumeration for 2-metal clusters ---
    # Detect metals connected through bridging donors and enumerate the
    # Cartesian product of their per-metal isomers to capture arrangements
    # that per-metal enumeration misses.
    try:
        bridging = _find_bridging_donors(mol)
        if bridging and len(results) < max_isomers:
            # Build metal cluster graph via bridging donors.
            metal_indices = [
                a.GetIdx() for a in mol.GetAtoms()
                if a.GetSymbol() in _METAL_SET
            ]
            if len(metal_indices) == 2:
                m1, m2 = metal_indices
                d1 = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(m1).GetNeighbors()]
                d2 = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(m2).GetNeighbors()]
                n1, n2 = len(d1), len(d2)
                if 2 <= n1 <= 9 and 2 <= n2 <= 9:
                    # Per-metal isomers.
                    dk1 = [dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset())) for d in d1]
                    dk2 = [dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset())) for d in d2]
                    uk1 = sorted(set(dk1), key=lambda k: (k[0], tuple(sorted(k[1]))))
                    uk2 = sorted(set(dk2), key=lambda k: (k[0], tuple(sorted(k[1]))))
                    kc1 = {k: i for i, k in enumerate(uk1)}
                    kc2 = {k: i for i, k in enumerate(uk2)}
                    dl1 = [f"{k[0]}{kc1[k]}" for k in dk1]
                    dl2 = [f"{k[0]}{kc2[k]}" for k in dk2]
                    cp1 = _chelate_pairs(mol, m1, d1)
                    cp2 = _chelate_pairs(mol, m2, d2)
                    al1 = {atom_idx: li for li, atom_idx in enumerate(d1)}
                    al2 = {atom_idx: li for li, atom_idx in enumerate(d2)}
                    clp1 = [frozenset([al1[sorted(cp)[0]], al1[sorted(cp)[1]]]) for cp in cp1
                            if sorted(cp)[0] in al1 and sorted(cp)[1] in al1]
                    clp2 = [frozenset([al2[sorted(cp)[0]], al2[sorted(cp)[1]]]) for cp in cp2
                            if sorted(cp)[0] in al2 and sorted(cp)[1] in al2]
                    ms1 = mol.GetAtomWithIdx(m1).GetSymbol()
                    ms2 = mol.GetAtomWithIdx(m2).GetSymbol()
                    iso1 = _enumerate_topological_isomers(dl1, n1, clp1, metal_symbol=ms1)
                    iso2 = _enumerate_topological_isomers(dl2, n2, clp2, metal_symbol=ms2)

                    # Scaffold-first approach: build M-bridge-M core first,
                    # then place non-bridging donors around each metal.
                    # Use the ETKDG template as scaffold base (it has correct
                    # M-bridge-M topology from SMILES).
                    import itertools as _it
                    max_combos = max(1, max_isomers - len(results))
                    scaffold = _build_multimetal_scaffold(
                        mol, metal_indices, bridging
                    )
                    topo_template_cids = _rank_template_conformers(mol, top_k=3) or [None]
                    combo_count = 0
                    for (cf1, pm1), (cf2, pm2) in _it.product(iso1, iso2):
                        if combo_count >= max_combos:
                            break
                        gn1 = cf1[0]
                        gn2 = cf2[0]
                        try:
                            # Build each metal's donors on the template, keeping
                            # the template's OTHER metal + bridge positions.
                            for _tpl_cid in topo_template_cids:
                                xyz1 = _build_topology_xyz(
                                    mol, m1, d1, pm1, gn1, False, conf_id=_tpl_cid
                                )
                                if xyz1 is None:
                                    continue
                                mol_tmp = Chem.RWMol(mol)
                                mol_tmp.RemoveAllConformers()
                                conf_tmp = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz1)
                                if conf_tmp is None:
                                    continue
                                cid_tmp = mol_tmp.AddConformer(conf_tmp, assignId=True)
                                # Rescale M-D for Metal1's arrangement.
                                _rescale_metal_donor_distances(mol_tmp, cid_tmp)
                                xyz_combined = _build_topology_xyz(
                                    mol_tmp.GetMol(), m2, d2, pm2, gn2, False,
                                    conf_id=cid_tmp,
                                )
                                if xyz_combined is not None:
                                    # Rescale M-D again for Metal2.
                                    mol_tmp2 = Chem.RWMol(mol)
                                    mol_tmp2.RemoveAllConformers()
                                    conf_c = _xyz_to_rdkit_conformer(mol_tmp2.GetMol(), xyz_combined)
                                    if conf_c is not None:
                                        cid_c = mol_tmp2.AddConformer(conf_c, assignId=True)
                                        _rescale_metal_donor_distances(mol_tmp2, cid_c)
                                        xyz_combined = _mol_to_xyz_conformer(mol_tmp2, cid_c)
                                    # Apply UFF with constraints.
                                    if apply_uff:
                                        xyz_combined = _optimize_xyz_openbabel_safe(
                                            xyz_combined, mol_template=mol
                                        )
                                    break
                            if xyz_combined is not None:
                                # Graph-based topology check.
                                if not _verify_topology_from_graph(xyz_combined, mol):
                                    continue
                                fp = _compute_coordination_fingerprint(
                                    Chem.RWMol(mol), 0, dtype_map=dtype_map
                                ) if False else None
                                label = f"multi-{gn1}/{gn2}"
                                results.append((xyz_combined, label))
                                combo_count += 1
                        except Exception as _cexc:
                            logger.debug("Multinuclear combo build failed: %s", _cexc)
                            continue
    except Exception as _mn_exc:
        logger.debug("Multinuclear enumeration failed: %s", _mn_exc)

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

            # --- Sulfite / Sulfonate (S-donor → O-donor, S with ≥3 O) ---
            if donor_sym == 'S':
                o_nbrs_all = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O' and n.GetIdx() not in bonded
                ]
                if len(o_nbrs_all) >= 2:
                    for alt_o in o_nbrs_all:
                        alternatives.append(
                            (metal_idx, donor_idx, alt_o.GetIdx(), 'sulfito-O')
                        )

            # --- Nitrosyl (N=O, N-donor → O-donor) ---
            if donor_sym == 'N':
                for o_nbr in donor.GetNeighbors():
                    if (o_nbr.GetSymbol() == 'O'
                            and o_nbr.GetIdx() not in bonded
                            and len(list(o_nbr.GetNeighbors())) == 1):
                        bond = mol.GetBondBetweenAtoms(donor_idx, o_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                            alternatives.append(
                                (metal_idx, donor_idx, o_nbr.GetIdx(), 'isonitrosyl-O')
                            )

            # --- Nitrosyl O-bound → N-bound ---
            if donor_sym == 'O':
                for n_nbr in donor.GetNeighbors():
                    if (n_nbr.GetSymbol() == 'N'
                            and n_nbr.GetIdx() not in bonded
                            and len(list(donor.GetNeighbors())) == 1):
                        bond = mol.GetBondBetweenAtoms(donor_idx, n_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                            alternatives.append(
                                (metal_idx, donor_idx, n_nbr.GetIdx(), 'nitrosyl-N')
                            )

            # --- Selenocyanate SeCN (Se-donor → N-donor via C) ---
            if donor_sym == 'Se':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for n_nbr in c_nbr.GetNeighbors():
                            if (n_nbr.GetSymbol() == 'N'
                                    and n_nbr.GetIdx() != donor_idx
                                    and n_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, n_nbr.GetIdx(), 'isoselenocyanato-N')
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
    max_template_tries: int = 8,
) -> List[Tuple[str, str]]:
    """Build linkage isomers by rewiring ambidentate ligands and generating XYZ.

    For each alternative coordination mode, rewires the metal bond and builds
    an idealized topology structure (OB UFF optimized).  The builder is tried
    against the top ``max_template_tries`` template conformers (ranked by
    :func:`_rank_template_conformers`) so one bad conformer does not silently
    discard valid linkage isomers.

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

            candidate_cids = _rank_template_conformers(
                alt_mol, top_k=max_template_tries
            )
            if not candidate_cids:
                candidate_cids = [None]

            for cid in candidate_cids:
                # Build WITHOUT UFF first (fast), check topology, THEN UFF.
                xyz = _build_topology_xyz(
                    alt_mol, metal_idx, donor_indices, perm, geom, False,
                    conf_id=cid,
                )
                if xyz is None:
                    continue
                if not _fragment_topology_ok(xyz, smiles):
                    logger.debug(
                        "linkage %s via template cid=%s: "
                        "fragment topology mismatch, trying next template",
                        type_label, cid,
                    )
                    continue
                if apply_uff:
                    xyz = _optimize_xyz_openbabel_safe(
                        xyz, mol_template=alt_mol
                    )
                results.append((xyz, type_label))
                break
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
    max_template_tries: int = 8,
) -> List[Tuple[str, str]]:
    """Generate isomers by swapping donor atoms with alternative binding sites.

    Only considers viable donors (atoms with available lone pairs) on the
    SAME ligand fragment as the current donor.  This avoids generating
    nonsensical structures from swapping donors across unrelated ligands
    or using C=O carbonyl oxygens as coordination donors.

    For each rewired candidate the builder is tried against the top
    ``max_template_tries`` template conformers (ranked by geometry quality)
    and the first result that passes ``_fragment_topology_ok`` is kept.  This
    decouples constitutional-isomer generation from a single, possibly poor,
    template conformer (see Codex findings for the reasoning).

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

                    if current_sym == alt_sym:
                        label = f'alt-{alt_sym}-isomer'
                    else:
                        label = f'alt-bind-{alt_sym}'

                    # Multi-template trial: iterate over the best-scored
                    # conformer IDs and keep the first XYZ that validates.
                    candidate_cids = _rank_template_conformers(
                        alt_mol, top_k=max_template_tries
                    )
                    if not candidate_cids:
                        # No usable template conformer — try de-novo once.
                        candidate_cids = [None]

                    accepted = False
                    for cid in candidate_cids:
                        # Build WITHOUT UFF first, check topology, THEN UFF.
                        xyz = _build_topology_xyz(
                            alt_mol, metal_idx, donor_indices, perm,
                            geom, False, conf_id=cid,
                        )
                        if xyz is None:
                            continue
                        if not _fragment_topology_ok(xyz, smiles):
                            logger.debug(
                                "alt-mode %s via template cid=%s: "
                                "fragment topology mismatch, trying next template",
                                label, cid,
                            )
                            continue
                        if apply_uff:
                            xyz = _optimize_xyz_openbabel_safe(
                                xyz, mol_template=alt_mol
                            )
                        results.append((xyz, label))
                        accepted = True
                        break

                    if not accepted:
                        logger.debug(
                            "alt-mode %s: no template conformer produced a "
                            "topology-consistent XYZ (tried %d)",
                            label, len(candidate_cids),
                        )
                except Exception as e:
                    logger.debug("Alternative binding mode failed: %s", e)
                    continue

    return results


def _enumerate_hapto_sigma_isomers(
    smiles: str,
    base_xyz: str,
    apply_uff: bool = True,
    max_isomers: int = 20,
) -> List[Tuple[str, str]]:
    """Enumerate σ-donor permutations for hapto complexes.

    The η-ring positions are kept fixed from ``base_xyz``.  Only the
    σ-donor ligand fragments are permuted via rigid-body Procrustes
    swaps, giving constitutionally distinct hapto isomers that differ
    in which σ-donor occupies which coordination slot around the
    η-ring(s).
    """
    if not RDKIT_AVAILABLE:
        return []
    try:
        import numpy as np
        import itertools as _it

        mol = _prepare_mol_for_embedding(smiles, hapto_approx=True)
        if mol is None:
            return []

        # Inject base_xyz as conformer.
        conf = _xyz_to_rdkit_conformer(mol, base_xyz)
        if conf is None:
            return []
        mol.RemoveAllConformers()
        cid = mol.AddConformer(conf, assignId=True)

        hapto_groups = _find_hapto_groups(mol)
        if not hapto_groups:
            return []

        hapto_atoms: set = set()
        for _midx, members in hapto_groups:
            hapto_atoms.update(members)

        results: List[Tuple[str, str]] = []
        dtype_map = _donor_type_map(mol)

        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            metal_idx = atom.GetIdx()
            all_donors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
            sigma_donors = [d for d in all_donors if d not in hapto_atoms]

            if len(sigma_donors) < 2:
                continue

            # Donor labels for σ-donors.
            donor_keys = [
                dtype_map.get(d, (mol.GetAtomWithIdx(d).GetSymbol(), frozenset()))
                for d in sigma_donors
            ]
            # If all σ-donors are equivalent → only 1 arrangement.
            if len(set(donor_keys)) <= 1:
                continue

            # Current σ-donor positions from conformer.
            conf_obj = mol.GetConformer(cid)
            sigma_positions = []
            for d in sigma_donors:
                p = conf_obj.GetAtomPosition(d)
                sigma_positions.append(np.array([p.x, p.y, p.z]))

            # Build ligand fragments per σ-donor (BFS excluding metal + η).
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

            donor_frag_atoms: Dict[int, set] = {}
            for d in sigma_donors:
                frag: set = set()
                stack = [d]
                visited: set = set()
                while stack:
                    node = stack.pop()
                    if node in visited or node in hapto_atoms:
                        continue
                    # Don't cross into OTHER σ-donor territories.
                    if node != d and node in sigma_donors:
                        continue
                    visited.add(node)
                    frag.add(node)
                    for nbr in adj.get(node, ()):
                        if nbr not in visited:
                            stack.append(nbr)
                donor_frag_atoms[d] = frag

            # Generate unique permutations of σ-donor labels.
            original_label_tuple = tuple(donor_keys)
            seen_labels: set = {original_label_tuple}

            for perm in _it.permutations(range(len(sigma_donors))):
                perm_labels = tuple(donor_keys[p] for p in perm)
                if perm_labels in seen_labels:
                    continue
                seen_labels.add(perm_labels)
                if len(results) >= max_isomers:
                    break

                # Build swapped XYZ: move fragment of donor perm[i] to
                # the position of donor i via rigid-body alignment.
                new_coords = np.zeros((mol.GetNumAtoms(), 3))
                for ai in range(mol.GetNumAtoms()):
                    p = conf_obj.GetAtomPosition(ai)
                    new_coords[ai] = [p.x, p.y, p.z]

                swap_ok = True
                for slot_idx in range(len(sigma_donors)):
                    src_donor = sigma_donors[perm[slot_idx]]
                    tgt_donor = sigma_donors[slot_idx]
                    if src_donor == tgt_donor:
                        continue
                    src_frag = sorted(donor_frag_atoms[src_donor])
                    tgt_pos = sigma_positions[slot_idx]
                    src_pos = sigma_positions[perm[slot_idx]]

                    # Translate fragment so donor lands at target position.
                    delta = tgt_pos - src_pos
                    orig_frag_coords = np.array([
                        [conf_obj.GetAtomPosition(ai).x,
                         conf_obj.GetAtomPosition(ai).y,
                         conf_obj.GetAtomPosition(ai).z]
                        for ai in src_frag
                    ])
                    for fi, ai in enumerate(src_frag):
                        new_coords[ai] = orig_frag_coords[fi] + delta

                # Write XYZ.
                lines = []
                for ai in range(mol.GetNumAtoms()):
                    sym = mol.GetAtomWithIdx(ai).GetSymbol()
                    x, y, z = new_coords[ai]
                    lines.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
                xyz_new = '\n'.join(lines) + '\n'

                if apply_uff:
                    try:
                        xyz_new = _optimize_xyz_openbabel_safe(
                            xyz_new, mol_template=mol
                        )
                    except Exception:
                        pass

                # Quality gate.
                if not _metal_donor_distances_realistic(xyz_new, mol):
                    continue

                try:
                    mol_chk = Chem.RWMol(mol)
                    mol_chk.RemoveAllConformers()
                    conf_chk = _xyz_to_rdkit_conformer(mol_chk.GetMol(), xyz_new)
                    if conf_chk is None:
                        continue
                    cid_chk = mol_chk.AddConformer(conf_chk, assignId=True)
                    if _has_atom_clash(mol_chk.GetMol(), cid_chk, min_dist=0.3):
                        continue
                    fp = _compute_coordination_fingerprint(
                        mol_chk.GetMol(), cid_chk, dtype_map=dtype_map
                    )
                    label = _classify_isomer_label(fp, mol_chk.GetMol())
                    if not label:
                        label = f'hapto-sigma-{len(results) + 1}'
                    results.append((xyz_new, label))
                except Exception:
                    continue

    except Exception as exc:
        logger.debug("Hapto sigma isomer enumeration failed: %s", exc)
    return results


def smiles_to_xyz_isomers(
    smiles: str,
    num_confs: int = 200,
    max_isomers: int = 50,
    apply_uff: bool = True,
    collapse_label_variants: bool = True,
    include_binding_mode_isomers: bool = True,
    deterministic: bool = True,
    hapto_approx: Optional[bool] = None,
) -> Tuple[List[Tuple[str, str]], Optional[str]]:
    """Generate distinct coordination isomers for a SMILES string.

    For non-metal molecules a single geometry is returned (delegates to
    ``smiles_to_xyz``).  For metal complexes multiple conformers are
    embedded, their coordination fingerprints are computed, and one
    representative per unique fingerprint is returned.

    Returns ``([(xyz_string, label), ...], error)``.

    When ``collapse_label_variants`` is ``True`` (default), additional cleanup
    merges numbered variants with the same base label (e.g. ``trans-1`` and
    ``trans-2``). Set it to ``False`` to keep these variants for workflows
    that prefer broader structural diversity.

    When ``deterministic`` is ``True`` (default), the metal-isomer path avoids
    Open Babel conformer injection (its ``make3D``/rotor pipeline is not fully
    reproducible across runs) and relies on seeded RDKit embedding only.
    """
    if not RDKIT_AVAILABLE:
        return [], "RDKit is not installed"

    has_metal = contains_metal(smiles)
    hapto_mode = _hapto_approx_enabled(hapto_approx)
    mol = None
    hapto_groups: List[Tuple[int, List[int]]] = []
    if has_metal:
        hapto_groups = _probe_hapto_groups_from_smiles(smiles)
        if hapto_groups and not hapto_mode:
            return [], _hapto_failfast_error(hapto_groups)
        if hapto_groups and hapto_mode:
            xyz, err = smiles_to_xyz(
                smiles,
                apply_uff=apply_uff,
                hapto_approx=True,
            )
            if err:
                return [], err
            results_hapto: List[Tuple[str, str]] = [(xyz, 'hapto-approx')]
            # Enumerate σ-donor permutations around fixed η-positions.
            try:
                hapto_sigma_iso = _enumerate_hapto_sigma_isomers(
                    smiles, xyz, apply_uff=apply_uff
                )
                if hapto_sigma_iso:
                    results_hapto.extend(hapto_sigma_iso)
            except Exception as _hse:
                logger.debug("Hapto sigma enumeration failed: %s", _hse)
            # For mixed hapto+regular multi-metal: DON'T early-return.
            # Fall through to multi-metal sampling augmentation so the
            # non-hapto metal's coordination isomers are explored.
            _n_metals_hapto = 0
            try:
                _mol_h = _prepare_mol_for_embedding(smiles, hapto_approx=True)
                if _mol_h:
                    _n_metals_hapto = sum(
                        1 for a in _mol_h.GetAtoms()
                        if a.GetSymbol() in _METAL_SET
                    )
            except Exception:
                pass
            if _n_metals_hapto <= 1:
                return results_hapto, None
            # Multi-metal hapto: the hapto builder already produced the
            # best possible Cp geometry (perfectly planar rings). ETKDG
            # sampling would produce conformers with broken Cp rings.
            # Return the hapto results directly.
            return results_hapto, None

    # Non-metal molecules: single geometry
    if not has_metal:
        xyz, err = smiles_to_xyz(smiles, apply_uff=apply_uff, hapto_approx=hapto_mode)
        if err:
            return [], err
        return [(xyz, '')], None

    # Prepare molecule for embedding (skip if already set by hapto multi-metal path)
    if mol is None:
        mol = _prepare_mol_for_embedding(smiles, hapto_approx=hapto_mode)
    if mol is None:
        # Fall back to single-conformer conversion
        xyz, err = smiles_to_xyz(smiles, apply_uff=apply_uff, hapto_approx=hapto_mode)
        if err:
            return [], err
        return [(xyz, '')], None

    # For metal complexes: prepend OB conformers to the pool so that
    # Avogadro-quality geometries are always considered during isomer search.
    conf_ids: List[int] = []
    if has_metal and OPENBABEL_AVAILABLE and not deterministic:
        try:
            # Non-deterministic enrichment path: augment RDKit pool with OB
            # conformers to increase diversity.
            _n_ob_restarts = 3
            _per_ob = max(10, int(num_confs) // _n_ob_restarts)
            ob_xyz_blocks: List[str] = []
            _ob_seen: set = set()
            ob_error: Optional[str] = None
            for _restart in range(_n_ob_restarts):
                _blocks, _err = _openbabel_generate_conformer_xyz(
                    smiles, num_confs=_per_ob, deterministic=False
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
                    "OB conformers injected for isomer search: %d",
                    len(ob_ids),
                )
            elif ob_error:
                logger.debug("OB conformer generation: %s", ob_error)
        except Exception as ob_exc:
            logger.debug("OB conformer generation exception: %s", ob_exc)
            mol.RemoveAllConformers()
            conf_ids = []

    # Embed multiple conformers with deterministic seed schedule.
    # Seeds are independent → parallelize with ThreadPoolExecutor.
    try:
        seeds = [31, 42, 7, 97, 13, 61, 83, 127, 211, 307, 401, 503]
        n_rounds = len(seeds)
        per_round = max(1, int(math.ceil(num_confs / n_rounds)))
        _n_workers = min(len(seeds), os.cpu_count() or 4)
        with concurrent.futures.ThreadPoolExecutor(max_workers=_n_workers) as _pool:
            _futs = [
                _pool.submit(_embed_multiple_confs_robust, mol, per_round, s)
                for s in seeds
            ]
            for _fut in concurrent.futures.as_completed(_futs):
                try:
                    conf_ids.extend(_fut.result())
                except Exception:
                    pass
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

    def _classify_one_conf(_cid, _relax):
        try:
            if _has_atom_clash(mol, _cid, min_dist=0.3):
                return None
            try:
                xyz_check = _mol_to_xyz_conformer(mol, _cid)
                if not _metal_donor_distances_realistic(xyz_check, mol):
                    return None
            except Exception:
                return None
            if _has_pi_ring_nonplanarity(mol, _cid):
                return None
            penalty = 0.0
            if _has_unphysical_metal_nonbonded_contact(mol, _cid):
                if not _relax:
                    return None
                penalty += 350.0
            if _has_unphysical_oco_geometry(mol, _cid):
                if not _relax:
                    return None
                penalty += 250.0
            if _has_atom_clash(mol, _cid):
                penalty += 500.0
            if _has_bad_geometry(mol, _cid):
                penalty += 300.0
            if _has_ligand_intertwining(mol, _cid):
                penalty += 200.0
            fp = _compute_coordination_fingerprint(mol, _cid, dtype_map=dtype_map)
            score = _geometry_quality_score(mol, _cid) + penalty
            label = _classify_isomer_label(fp, mol)
            return (fp, label, _cid, score)
        except Exception:
            return None

    def _collect_fp_label_pairs(relax_hard_chem_filters: bool = False) -> List[Tuple[tuple, str, int, float]]:
        _n_cw = min(len(conf_ids), os.cpu_count() or 4) if conf_ids else 1
        if _n_cw > 1 and len(conf_ids) > 4:
            with concurrent.futures.ThreadPoolExecutor(max_workers=_n_cw) as _cp:
                _futs = [
                    _cp.submit(_classify_one_conf, c, relax_hard_chem_filters)
                    for c in conf_ids
                ]
                return [
                    r for r in (f.result() for f in _futs) if r is not None
                ]
        return [
            r for r in (_classify_one_conf(c, relax_hard_chem_filters) for c in conf_ids)
            if r is not None
        ]

    fp_label_pairs: List[Tuple[tuple, str, int, float]] = _collect_fp_label_pairs(
        relax_hard_chem_filters=False
    )
    if not fp_label_pairs and conf_ids and has_metal:
        logger.debug(
            "Strict conformer filters rejected all sampled structures; "
            "retrying with relaxed OCO/nonbonded penalties."
        )
        fp_label_pairs = _collect_fp_label_pairs(relax_hard_chem_filters=True)

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
                if rmsd < 2.5:
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
        _fb_xyz, _fb_err = smiles_to_xyz(smiles, apply_uff=apply_uff, hapto_approx=hapto_mode)
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
        relaxed_fragment_results: List[Tuple[str, str]] = []
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
                        xyz = _optimize_xyz_openbabel_safe(
                            xyz,
                            mol_template=mol,
                            smiles=smiles,
                            apply_template_constraints=True,
                        )
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
            # Fragment topology check: organic ligand fragments must match the
            # original SMILES (catches broken/fused ligands while preserving
            # fac/mer isomers which have identical fragment sets).
            # For multi-metal complexes, skip strict fragment check — the
            # complex bridging topology causes frequent false-positive
            # mismatches in OB bond perception.
            _n_metals_in_mol = sum(
                1 for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET
            )
            if not _fragment_topology_ok(xyz, smiles):
                if _n_metals_in_mol >= 2:
                    # Accept multi-metal conformers with relaxed topology
                    logger.debug(
                        "Accepting conformer %d despite fragment mismatch "
                        "(multi-metal complex)", cid)
                else:
                    logger.debug("Skipping conformer %d: fragment topology mismatch", cid)
                    if (
                        _fragment_topology_relaxed_fallback_ok(xyz, smiles)
                        and _xyz_passes_final_geometry_checks(xyz, mol)
                    ):
                        relaxed_fragment_results.append((xyz, display))
                    continue
            results.append((xyz, display))

        if not results:
            if relaxed_fragment_results:
                logger.debug(
                    "Using %d conformer(s) despite fragment-topology mismatch "
                    "after stricter checks.",
                    len(relaxed_fragment_results),
                )
                results = relaxed_fragment_results
            else:
                _fb_xyz, _fb_err = smiles_to_xyz(smiles, apply_uff=apply_uff, hapto_approx=hapto_mode)
                if _fb_err:
                    return [], _fb_err
                results = [(_fb_xyz, '')]

    # --- Topological enumerator: guarantee completeness ---
    # For MONO-metallic: use topology enumerator (guaranteed complete).
    # For MULTI-metallic: the topo builder destroys bridge ligand geometry
    # by moving metals independently. Instead, rely on ETKDG sampling
    # (which preserves ligand topology) + M-D rescaling + graph check.
    try:
        _n_metals_total = sum(
            1 for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET
        ) if RDKIT_AVAILABLE else 0
    except Exception:
        _n_metals_total = 0
    if has_metal and _n_metals_total <= 1:
        try:
            # Collect fingerprints from sampling results
            existing_fps: set = set()
            topo_mol = mol
            dtype_map_topo = _donor_type_map(topo_mol)
            for existing_xyz, _existing_display in results:
                try:
                    mol_tmp = Chem.RWMol(topo_mol)
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
            existing_xyz_keys = {
                "\n".join(l.strip() for l in xyz.splitlines() if l.strip())
                for xyz, _d in results
            }

            topo_results = _generate_topological_isomers(
                topo_mol, smiles, apply_uff=apply_uff, max_isomers=max_isomers
            )

            for topo_xyz, topo_label in topo_results:
                # Compute fingerprint of topo structure for dedup
                topo_fp = None
                try:
                    mol_tmp = Chem.RWMol(topo_mol)
                    mol_tmp.RemoveAllConformers()
                    conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), topo_xyz)
                    if conf is not None:
                        cid = mol_tmp.AddConformer(conf, assignId=True)
                        topo_fp = _compute_coordination_fingerprint(
                            mol_tmp.GetMol(), cid, dtype_map=dtype_map_topo
                        )
                except Exception:
                    pass

                # Always let topo isomers through; the final geometry-based
                # dedup selects the best-scoring candidate per base label.
                # Sampling conformers are UFF-distorted, so topo versions
                # (donors pinned to ideal polyhedron) typically win.
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
                # Graph-based topology verification: checks every bond in
                # the template graph against actual XYZ distances.  No OB
                # perception, no roundtrip — works for mono AND multi-metal.
                if not _verify_topology_from_graph(topo_xyz, topo_mol):
                    logger.debug("Skipping topo isomer %s: graph topology check failed", display)
                    continue
                topo_key = "\n".join(l.strip() for l in topo_xyz.splitlines() if l.strip())
                if topo_key in existing_xyz_keys:
                    logger.debug("Skipping topo isomer %s: duplicate XYZ", display)
                    continue
                existing_displays.add(display)
                existing_xyz_keys.add(topo_key)
                if topo_fp is not None:
                    existing_fps.add(topo_fp)
                results.append((topo_xyz, display))
        except Exception as _topo_exc:
            logger.debug("Topological isomer generation failed: %s", _topo_exc)

    # --- Multi-metal sampling augmentation ---
    # For multi-metallic systems: augment with extra ETKDG conformers.
    # EXCEPTION: hapto complexes — ETKDG can't build Cp rings correctly
    # (they come out non-planar). The hapto builder already produced
    # perfect Cp geometry, so skip ETKDG augmentation for hapto systems.
    _skip_mm_augmentation = bool(hapto_groups and hapto_mode)
    if has_metal and _n_metals_total >= 2 and len(results) < max_isomers and not _skip_mm_augmentation:
        try:
            _extra_seeds = [137, 251, 353, 461, 571, 683, 797, 911, 1013, 1117]
            _extra_ids: List[int] = []
            _n_extra = min(len(_extra_seeds), os.cpu_count() or 4)
            with concurrent.futures.ThreadPoolExecutor(max_workers=_n_extra) as _xp:
                _xfuts = [
                    _xp.submit(_embed_multiple_confs_robust, mol, 3, s)
                    for s in _extra_seeds
                ]
                for _xf in concurrent.futures.as_completed(_xfuts):
                    try:
                        _extra_ids.extend(_xf.result())
                    except Exception:
                        pass
            logger.debug("Multi-metal augmentation: %d extra conformers", len(_extra_ids))

            # Classify + dedup with existing results
            _existing_fps_mm: set = set()
            for _xyz, _d in results:
                try:
                    _mt = Chem.RWMol(mol)
                    _mt.RemoveAllConformers()
                    _c = _xyz_to_rdkit_conformer(_mt.GetMol(), _xyz)
                    if _c is not None:
                        _ci = _mt.AddConformer(_c, assignId=True)
                        _fp = _compute_coordination_fingerprint(
                            _mt.GetMol(), _ci, dtype_map=dtype_map
                        )
                        _existing_fps_mm.add(_fp)
                except Exception:
                    pass

            for _xid in _extra_ids:
                if len(results) >= max_isomers:
                    break
                try:
                    _xyz = _mol_to_xyz_conformer(mol, _xid)
                    _xyz = _optimize_xyz_openbabel_safe(_xyz, mol_template=mol)
                    if not _verify_topology_from_graph(_xyz, mol):
                        continue
                    _mt = Chem.RWMol(mol)
                    _mt.RemoveAllConformers()
                    _c = _xyz_to_rdkit_conformer(_mt.GetMol(), _xyz)
                    if _c is None:
                        continue
                    _ci = _mt.AddConformer(_c, assignId=True)
                    _fp = _compute_coordination_fingerprint(
                        _mt.GetMol(), _ci, dtype_map=dtype_map
                    )
                    if _fp in _existing_fps_mm:
                        continue
                    _existing_fps_mm.add(_fp)
                    _lbl = _classify_isomer_label(_fp, _mt.GetMol())
                    if not _lbl:
                        unknown_counter += 1
                        _lbl = f'Isomer {unknown_counter}'
                    results.append((_xyz, _lbl))
                except Exception:
                    continue
        except Exception as _mm_exc:
            logger.debug("Multi-metal augmentation failed: %s", _mm_exc)

    # --- Linkage isomers ---
    if has_metal and include_binding_mode_isomers:
        try:
            link_results = _generate_linkage_isomers(mol, smiles, apply_uff=apply_uff)
            for _lxyz, _llabel in link_results:
                if not _metal_donor_distances_realistic(_lxyz, mol):
                    logger.debug("Skipping linkage isomer %s: unphysical M-D distance", _llabel)
                    continue
                if _fragment_topology_ok(_lxyz, smiles):
                    results.append((_lxyz, _llabel))
                else:
                    logger.debug("Skipping linkage isomer %s: fragment topology mismatch", _llabel)
        except Exception as _link_exc:
            logger.debug("Linkage isomer generation failed: %s", _link_exc)

    # --- Alternative binding-site exploration ---
    if has_metal and include_binding_mode_isomers:
        try:
            existing_displays = {display for _, display in results}
            existing_base = {re.sub(r'-\d+$', '', d) for d in existing_displays}
            alt_results = _generate_alternative_binding_modes(
                mol, smiles, apply_uff=apply_uff
            )
            for alt_xyz, alt_label in alt_results:
                if alt_label not in existing_base:
                    if not _metal_donor_distances_realistic(alt_xyz, mol):
                        logger.debug(
                            "Skipping alt-binding isomer %s: unphysical M-D distance",
                            alt_label,
                        )
                        continue
                    if not _fragment_topology_ok(alt_xyz, smiles):
                        logger.debug("Skipping alt-binding isomer %s: fragment topology mismatch", alt_label)
                        continue
                    results.append((alt_xyz, alt_label))
                    existing_base.add(alt_label)
        except Exception as _alt_exc:
            logger.debug("Alternative binding mode generation failed: %s", _alt_exc)

    # --- Final label dedup (safety net) ---
    # Collapse entries that share the same base label (e.g. "trans-1" and
    # "trans-2" that slipped through fingerprint/RMSD dedup).  Pick the
    # geometrically best candidate per base label (sampling conformers
    # from UFF-distorted ETKDG are usually worse than the topology-
    # enumerator output which pins donors to the idealized polyhedron).
    # "Isomer N" labels use spaces not dashes, so they are never collapsed.
    if collapse_label_variants and results and has_metal:
        try:
            score_mol = mol

            def _score_xyz(_xyz: str) -> float:
                if score_mol is None:
                    return float("inf")
                try:
                    tmp = Chem.RWMol(score_mol)
                    tmp.RemoveAllConformers()
                    conf = _xyz_to_rdkit_conformer(tmp.GetMol(), _xyz)
                    if conf is None:
                        return float("inf")
                    cid = tmp.AddConformer(conf, assignId=True)
                    return float(_geometry_quality_score(tmp.GetMol(), cid))
                except Exception:
                    return float("inf")

            # Pick best-score candidate per base label, preserve first
            # appearance order for the final result list.
            base_best: Dict[str, Tuple[int, float]] = {}
            scored_cache: List[float] = []
            for idx, (xyz, lbl) in enumerate(results):
                base = re.sub(r'-\d+$', '', lbl) if lbl else ''
                score = _score_xyz(xyz)
                scored_cache.append(score)
                if base not in base_best or score < base_best[base][1]:
                    base_best[base] = (idx, score)

            keep_indices = {idx for idx, _s in base_best.values()}
            new_results: List[Tuple[str, str]] = []
            for idx, (xyz, lbl) in enumerate(results):
                if idx not in keep_indices:
                    logger.debug(
                        "Final dedup: dropping %r (base=%r, score=%.3f) — better kept",
                        lbl,
                        re.sub(r'-\d+$', '', lbl) if lbl else '',
                        scored_cache[idx],
                    )
                    continue
                new_results.append((xyz, re.sub(r'-\d+$', '', lbl) if lbl else lbl))
            results = new_results
        except Exception as _dedup_exc:
            logger.debug("Geometry-based dedup failed, falling back to first-keep: %s", _dedup_exc)
            _seen_base: Dict[str, int] = {}
            _keep: List[bool] = [True] * len(results)
            for _idx, (_, _lbl) in enumerate(results):
                _base = re.sub(r'-\d+$', '', _lbl) if _lbl else ''
                if _base in _seen_base:
                    _keep[_idx] = False
                else:
                    _seen_base[_base] = _idx
            results = [
                (xyz, re.sub(r'-\d+$', '', lbl))
                for (xyz, lbl), keep in zip(results, _keep) if keep
            ]

    # --- Final output gate: graph-based topology verification ---
    # Every output structure must preserve the bond topology from the
    # input SMILES.  Uses _verify_topology_from_graph which checks
    # every bond distance against the template graph — no OB perception.
    if has_metal and results:
        verified: List[Tuple[str, str]] = []
        for xyz, lbl in results:
            if _verify_topology_from_graph(xyz, mol):
                verified.append((xyz, lbl))
            else:
                logger.info(
                    "Output gate rejected isomer %r: graph topology violated", lbl
                )
        if verified:
            results = verified

    return results, None


def smiles_to_xyz_quick_hapto_previews(
    smiles: str,
    hapto_approx: Optional[bool] = None,
) -> List[Tuple[str, str]]:
    """Return extra quick-convert preview structures for hapto-specific builders."""
    _ = _hapto_approx_enabled(hapto_approx)
    return list(_HAPTO_QUICK_PREVIEW_CACHE.get(smiles, []))


def smiles_to_xyz_quick(
    smiles: str,
    hapto_approx: Optional[bool] = None,
) -> Tuple[Optional[str], Optional[str]]:
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
    hapto_mode = _hapto_approx_enabled(hapto_approx)

    if has_metal:
        hapto_groups = _probe_hapto_groups_from_smiles(smiles)
        if hapto_groups:
            if not hapto_mode:
                return None, _hapto_failfast_error(hapto_groups)
            # Experimental eta/hapto fallback: reuse the main converter with
            # approximation enabled and no UFF to keep quick-mode behavior.
            return smiles_to_xyz(smiles, apply_uff=False, hapto_approx=True)
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
    hapto_approx: Optional[bool] = None,
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

    import numpy as np  # local import: keeps module import cheap

    try:
        has_metal = contains_metal(smiles)
        hapto_mode = _hapto_approx_enabled(hapto_approx)
        hapto_groups = _probe_hapto_groups_from_smiles(smiles) if has_metal else []
        _HAPTO_QUICK_PREVIEW_CACHE[smiles] = []
        legacy_hapto_xyz: Optional[str] = None
        if hapto_groups and not hapto_mode:
            error = _hapto_failfast_error(hapto_groups)
            logger.warning(error)
            return None, error
        if has_metal and hapto_mode and hapto_groups:
            # For single-hapto systems, try legacy path as backup.
            # For multi-hapto (>=2 groups), skip the slow legacy path.
            if len(hapto_groups) <= 1:
                legacy_xyz, legacy_err = _try_multiple_strategies(smiles)
                if legacy_err is None and legacy_xyz:
                    legacy_hapto_xyz = legacy_xyz
                    logger.info(
                        "Hapto detected: collected legacy multi-strategy candidate "
                        "for universal candidate selection."
                    )
        method = None

        # Cluster compounds (borane/carborane cages, etc.): extremely dense
        # ring-closure topologies cause ETKDG to hang.  Route directly to
        # the manual metal embed fallback which places atoms geometrically.
        if has_metal and not hapto_groups:
            _cluster_mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if _cluster_mol is not None:
                _nb = _cluster_mol.GetNumBonds()
                _na = _cluster_mol.GetNumAtoms()
                if _na > 0 and _nb > 2 * _na:
                    logger.info(
                        "Detected cluster topology (%d bonds / %d atoms = %.1f), "
                        "using manual metal embed", _nb, _na, _nb / _na)
                    xyz_content, manual_err = _manual_metal_embed(smiles)
                    if xyz_content:
                        if apply_uff:
                            xyz_content = _optimize_xyz_openbabel_safe(
                                xyz_content, mol_template=mol
                            )
                        if output_path:
                            Path(output_path).write_text(xyz_content, encoding='utf-8')
                        return xyz_content, None

        # For metal-nitrogen coordination complexes (both neutral and charged notation),
        # use the multi-strategy approach that tries multiple parsing methods
        if _is_metal_nitrogen_complex(smiles) and not hapto_groups:
            logger.info("Detected metal-nitrogen complex, using multi-strategy approach")
            return _try_multiple_strategies(smiles, output_path)

        # Parse SMILES - try stk first for metal complexes, then RDKit.
        # IMPORTANT: Prefer the original SMILES before charge-normalized
        # variants to avoid introducing artificial oxidation states
        # (e.g., neutral Cu complexes rewritten as Cu+2).
        mol = None
        normalized_smiles = _normalize_metal_smiles(smiles)
        if has_metal and STK_AVAILABLE and len(hapto_groups) <= 1:
            # Skip STK for multi-hapto systems (slow and usually fails).
            # Guard with timeout: stk internally runs ETKDG which can hang.
            _stk_mol_holder = [None]

            def _stk_parse():
                try:
                    bb = stk.BuildingBlock(smiles)
                    _stk_mol_holder[0] = bb.to_rdkit_mol()
                except Exception:
                    _stk_mol_holder[0] = None

            _stk_t = threading.Thread(target=_stk_parse, daemon=True)
            _stk_t.start()
            _stk_t.join(timeout=_EMBED_TIMEOUT)
            if _stk_t.is_alive():
                logger.info("stk conversion timed out after %.1fs, falling back to RDKit", _EMBED_TIMEOUT)
                mol = None
                method = None
            elif _stk_mol_holder[0] is not None:
                mol = _stk_mol_holder[0]
                method = "stk"
            else:
                logger.info("stk conversion failed, falling back to RDKit")
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
                    # For hapto complexes: try unsanitized mol (no valence check)
                    # so we can still go through the hapto correction path.
                    if hapto_groups and hapto_mode:
                        try:
                            mol_unsan = Chem.MolFromSmiles(smiles, sanitize=False)
                            if mol_unsan is not None:
                                mol_unsan.UpdatePropertyCache(strict=False)
                                mol = mol_unsan
                                method = "RDKit (unsanitized for hapto)"
                                rdkit_note = None
                                logger.info(
                                    "Using unsanitized mol for hapto path "
                                    "(valence error bypassed)")
                        except Exception:
                            pass  # fall through to legacy

                    if mol is None:
                        legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                        if legacy_err is None and legacy_xyz:
                            if hapto_groups and hapto_mode:
                                # Save as fallback but don't return yet —
                                # try hapto path first via unsanitized mol
                                legacy_hapto_xyz = legacy_xyz
                            else:
                                if output_path:
                                    Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                                return legacy_xyz, None
                # Last resort: try multi-strategy approach for metal complexes
                if has_metal and mol is None:
                    if hapto_groups and hapto_mode:
                        error = (
                            "Failed to parse hapto complex for experimental approximation. "
                            "Try simplifying the SMILES around eta-coordination."
                        )
                        logger.error(error)
                        return None, error
                    logger.info("Trying multi-strategy fallback for unparseable metal SMILES")
                    return _try_multiple_strategies(smiles, output_path)
                error = f"Failed to parse SMILES string: {rdkit_note}"
                logger.error(error)
                return None, error

        if has_metal and hapto_mode:
            try:
                parsed_hapto = _find_hapto_groups(mol)
                if parsed_hapto:
                    mol, n_removed = _apply_hapto_approximation(mol, parsed_hapto)
                    logger.info(
                        "Applied experimental hapto approximation: "
                        "%d group(s), %d bond(s) converted to dative.",
                        len(parsed_hapto), n_removed,
                    )
            except Exception as e:
                logger.debug("Hapto approximation failed, continuing with original graph: %s", e)

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
                if hapto_mode:
                    mol = _fix_hapto_donor_h(mol)
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

        # Hapto path: ETKDG for topology preservation + sphere correction
        # for correct hapto geometry.  Sphere scaffold only as last fallback.
        if has_metal and hapto_mode and hapto_groups:
            # Re-identify hapto groups on the current mol (atom indices may have
            # shifted due to RemoveHs/AddHs/dative bond conversion above).
            hapto_groups = _find_hapto_groups(mol)
            if not hapto_groups:
                hapto_groups = _probe_hapto_groups_from_smiles(smiles)
            preview_candidates: List[Tuple[str, str]] = []
            hapto_candidate_mols: List[Tuple[str, object]] = []

            # Primary hybrid path: analytical hapto scaffold plus rigidly
            # aligned ligand fragments. This avoids global ETKDG on the full
            # metal graph where multi-hapto systems are most brittle.
            try:
                hybrid_variant_plans = _secondary_metal_variant_plans(mol, hapto_groups)
                seen_hybrid_xyz: set = set()
                for plan_idx, (hybrid_label, variant_plan) in enumerate(hybrid_variant_plans):
                    hybrid_mol = _build_hybrid_hapto_complex(
                        mol,
                        hapto_groups,
                        preview_store=preview_candidates if plan_idx == 0 else None,
                        secondary_variant_plan=variant_plan,
                    )
                    if hybrid_mol is None:
                        continue
                    try:
                        hybrid_xyz_key = "\n".join(
                            line.strip()
                            for line in _mol_to_xyz(hybrid_mol).splitlines()
                            if line.strip()
                        )
                    except Exception:
                        hybrid_xyz_key = ""
                    if hybrid_xyz_key and hybrid_xyz_key in seen_hybrid_xyz:
                        continue
                    if hybrid_xyz_key:
                        seen_hybrid_xyz.add(hybrid_xyz_key)
                    hapto_candidate_mols.append((hybrid_label, hybrid_mol))
                if seen_hybrid_xyz:
                    logger.info(
                        "Hybrid hapto scaffold/fragment builder produced %d ranked candidate(s)",
                        len(seen_hybrid_xyz),
                    )
            except Exception as e:
                logger.debug("Hybrid hapto scaffold/fragment builder failed: %s", e)
            _HAPTO_QUICK_PREVIEW_CACHE[smiles] = list(preview_candidates)

            # Fix O- atoms with double bonds that cause valence errors
            # during ETKDG embedding (e.g. C=[O-] in metal chelates).
            # Temporarily neutralize for embedding; formal charges don't
            # affect distance geometry.
            _fixed_o_atoms = []
            for _oa in mol.GetAtoms():
                if (_oa.GetSymbol() == 'O'
                        and _oa.GetFormalCharge() == -1
                        and any(b.GetBondType() == Chem.BondType.DOUBLE
                                for b in _oa.GetBonds())):
                    _oa.SetFormalCharge(0)
                    _fixed_o_atoms.append(_oa.GetIdx())
            if _fixed_o_atoms:
                try:
                    mol.UpdatePropertyCache(strict=False)
                except Exception:
                    pass

            # SECONDARY: Two-phase ETKDG embedding
            # Phase 1: ETKDG + analytical hapto correction → fix hapto atoms
            # Phase 2: candidate selection via topology/geometry filter.
            etkdg_seeds = [42, 1337, 2027, 7, 97, 13, 271, 911]
            for seed in etkdg_seeds:
                try:
                    mol_try = Chem.Mol(mol)
                    params = AllChem.ETKDGv3()
                    params.randomSeed = seed
                    params.useRandomCoords = True
                    params.enforceChirality = False
                    r = _embed_with_timeout(mol_try, params)
                    if r != 0:
                        r = _embed_with_timeout(
                            mol_try,
                            _make_random_embed_params(seed),
                        )
                    if r != 0:
                        continue
                    cid = int(r)

                    # Phase 1: Analytical hapto correction
                    try:
                        _correct_hapto_geometry(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    # Phase 1.5: Fix sigma metal-ligand distances.
                    # ETKDG treats M-L bonds like organic bonds (~1.5Å)
                    # which is too short.  Scale entire ligand fragments
                    # rigidly to correct M-L distances.
                    try:
                        from collections import deque as _deque
                        _conf = mol_try.GetConformer(cid)
                        hapto_atom_set = set()
                        metal_set_idx = set()
                        for _mi, catoms in hapto_groups:
                            metal_set_idx.add(_mi)
                            for _ca in catoms:
                                hapto_atom_set.add(_ca)
                        for _ai in range(mol_try.GetNumAtoms()):
                            if mol_try.GetAtomWithIdx(_ai).GetSymbol() in _METAL_SET:
                                metal_set_idx.add(_ai)

                        for _ai in metal_set_idx:
                            _atom = mol_try.GetAtomWithIdx(_ai)
                            m_sym = _atom.GetSymbol()
                            mp = _conf.GetAtomPosition(_ai)
                            m_pos = np.array([mp.x, mp.y, mp.z])
                            for _nbr in _atom.GetNeighbors():
                                _ni = _nbr.GetIdx()
                                if _ni in hapto_atom_set:
                                    continue  # handled by Phase 1
                                if _ni in metal_set_idx:
                                    continue  # M-M bonds
                                l_sym = _nbr.GetSymbol()
                                lp = _conf.GetAtomPosition(_ni)
                                l_pos = np.array([lp.x, lp.y, lp.z])
                                cur_d = float(np.linalg.norm(l_pos - m_pos))
                                if cur_d < 1e-8:
                                    continue
                                target_d = float(
                                    _get_ml_bond_length(m_sym, l_sym))
                                if abs(cur_d - target_d) / target_d < 0.15:
                                    continue
                                # BFS to find entire fragment attached to
                                # this donor (not crossing metals)
                                frag = set()
                                q = _deque([_ni])
                                frag.add(_ni)
                                while q:
                                    cur = q.popleft()
                                    for fn in mol_try.GetAtomWithIdx(cur).GetNeighbors():
                                        fi = fn.GetIdx()
                                        if (fi not in frag
                                                and fi not in metal_set_idx
                                                and fi not in hapto_atom_set):
                                            frag.add(fi)
                                            q.append(fi)
                                # Rigid translation of entire fragment
                                delta = (target_d - cur_d) / cur_d * (
                                    l_pos - m_pos)
                                for fi in frag:
                                    fp = _conf.GetAtomPosition(fi)
                                    fv = np.array([fp.x, fp.y, fp.z])
                                    new_fv = fv + delta
                                    _conf.SetAtomPosition(
                                        fi, Point3D(float(new_fv[0]),
                                                    float(new_fv[1]),
                                                    float(new_fv[2])))
                    except Exception:
                        pass

                    # Phase 1.6: BFS propagation for non-hapto atoms
                    try:
                        _propagate_non_hapto_atoms(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    try:
                        _enforce_donor_pi_coplanarity(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    # Fix secondary metal positions AFTER coplanarity
                    try:
                        _fix_secondary_metal_distances(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    try:
                        _final_clash_resolution(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    hapto_candidate_mols.append((f"etkdg-seed-{seed}", mol_try))
                except Exception:
                    continue

            # FALLBACK: Sphere-based constructive scaffold builder
            try:
                mol_scaffold = Chem.Mol(mol)
                if _build_hapto_scaffold(mol_scaffold, hapto_groups):
                    try:
                        _enforce_donor_pi_coplanarity(mol_scaffold, 0, hapto_groups)
                    except Exception:
                        pass
                    try:
                        _final_clash_resolution(mol_scaffold, 0, hapto_groups)
                    except Exception:
                        pass
                    hapto_candidate_mols.append(("scaffold", mol_scaffold))
            except Exception as e:
                logger.debug("Sphere scaffold failed: %s", e)

            best_mol = _select_best_hapto_candidate(
                smiles,
                hapto_groups,
                hapto_candidate_mols,
                apply_uff=apply_uff,
            )

            if best_mol is None:
                if legacy_hapto_xyz:
                    if output_path:
                        Path(output_path).write_text(legacy_hapto_xyz, encoding='utf-8')
                    return legacy_hapto_xyz, None
                return None, "Hapto-approx embedding failed"

            mol = best_mol
            xyz_content = _mol_to_xyz(mol)
            if output_path:
                Path(output_path).write_text(xyz_content, encoding='utf-8')
            return xyz_content, None

        # Generate 3D coordinates using a hybrid OB+RDKit conformer pool.
        # For metal complexes: OB WeightedRotorSearch in 3 independent restarts
        # (non-deterministic diversity) + RDKit ETKDG with
        # 12 diverse fixed seeds → up to ~500 conformers total.
        # Best geometry is selected by _geometry_quality_score.
        #
        # For large molecules (>50 heavy atoms): skip the expensive OB/ETKDG
        # multi-conformer pipeline.  OB's SystematicRotorSearch is
        # combinatorially explosive and holds the GIL (no thread-based timeout).
        # ETKDG can also hang on complex metal ring systems.  These are initial
        # geometries destined for GOAT/xTB anyway.
        _n_heavy_atoms = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
        _skip_expensive_pipeline = _n_heavy_atoms > 50
        result = -1
        if has_metal and not _skip_expensive_pipeline:
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
            # Use a thread-pool with a global timeout to prevent hangs on
            # complex metal ring systems where ETKDG can stall.
            seeds = [31, 42, 7, 97, 13, 61, 83, 127, 211, 307, 401, 503]
            per_seed = max(1, 200 // len(seeds))
            _etkdg_deadline = _EMBED_TIMEOUT * 3  # total budget for all seeds
            _etkdg_start = __import__('time').monotonic()
            try:
                for seed in seeds:
                    if __import__('time').monotonic() - _etkdg_start > _etkdg_deadline:
                        logger.debug("ETKDG multi-seed budget exhausted after %.1fs", _etkdg_deadline)
                        break
                    params_multi = AllChem.ETKDGv3()
                    params_multi.randomSeed = seed
                    params_multi.useRandomCoords = True
                    params_multi.enforceChirality = False
                    try:
                        params_multi.clearConfs = False  # append to OB conformers
                    except Exception:
                        pass
                    # Run with timeout to avoid hanging on single seed
                    _multi_ids = [None]
                    def _do_multi(m=mol, n=per_seed, p=params_multi):
                        _multi_ids[0] = list(AllChem.EmbedMultipleConfs(m, numConfs=n, params=p))
                    _mt = threading.Thread(target=_do_multi, daemon=True)
                    _mt.start()
                    _mt.join(timeout=_EMBED_TIMEOUT)
                    if _mt.is_alive():
                        logger.debug("EmbedMultipleConfs timed out for seed %d", seed)
                        continue
                    if _multi_ids[0]:
                        all_conf_ids.extend(_multi_ids[0])
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
            if _skip_expensive_pipeline:
                # For large molecules: go straight to permissive embed.
                # ETKDG with standard knowledge can hang on complex metal
                # ring systems (GIL-holding C code, thread timeout ineffective).
                logger.info(
                    "Large molecule (%d heavy atoms): using permissive embed",
                    _n_heavy_atoms,
                )
                result = AllChem.EmbedMolecule(mol, _make_random_embed_params(42))
            else:
                params = AllChem.ETKDGv3()
                params.randomSeed = 42  # For reproducibility

                result = _embed_with_timeout(mol, params)

                if result != 0:
                    # Try with random coordinates if ETKDG fails
                    logger.warning("ETKDG embedding failed, trying random coordinates")
                    result = _embed_with_timeout(mol, AllChem.ETKDG())

                if result != 0:
                    # Last resort: permissive embed without ETKDG knowledge
                    logger.warning("ETKDG timed out or failed, trying permissive embed")
                    result = AllChem.EmbedMolecule(mol, _make_random_embed_params(42))

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

        # For metal complexes: UFF with universal coordination constraints.
        if apply_uff and has_metal:
            xyz_content = _optimize_xyz_openbabel_safe(
                xyz_content, mol_template=mol
            )

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

        # Save explicit H counts and original NoImplicit state, then zero out
        # H and set NoImplicit=True to prevent valence errors during embedding
        # (e.g. [CH+] with 3 bonds = valence 4 > permitted 3 for C+).
        saved_explicit_h = {}
        originally_no_implicit = set()
        for atom in mol.GetAtoms():
            if atom.GetNoImplicit():
                originally_no_implicit.add(atom.GetIdx())
            eh = atom.GetNumExplicitHs()
            if eh > 0:
                saved_explicit_h[atom.GetIdx()] = eh
                atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(True)

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
            result = _embed_with_timeout(mol, params)
        except (TypeError, Exception):
            result = -1

        if result != 0:
            # Relaxed embedding for complex metal geometries (ferrocene,
            # multi-center Pd complexes) where ETKDG distance bounds fail.
            relaxed = AllChem.EmbedParameters()
            relaxed.randomSeed = 42
            relaxed.useRandomCoords = True
            relaxed.useExpTorsionAnglePrefs = False
            relaxed.useBasicKnowledge = False
            relaxed.enforceChirality = False
            try:
                relaxed.ignoreSmoothingFailures = True
            except Exception:
                pass
            result = _embed_with_timeout(mol, relaxed)

        if result != 0:
            return None, "Failed to generate 3D coordinates (unsanitized)"

        # Restore explicit H and add hydrogens with coordinates
        if not _prefer_no_sanitize(smiles):
            try:
                # Restore saved explicit H counts
                for idx, eh in saved_explicit_h.items():
                    mol.GetAtomWithIdx(idx).SetNumExplicitHs(eh)
                # Allow implicit H calculation only for organic subset atoms
                # (originally NoImplicit=False, e.g. bare C in phenyl rings).
                # Bracket atoms like [C], [CH+], [C@@] keep NoImplicit=True
                # to respect the SMILES-specified H count.
                for atom in mol.GetAtoms():
                    if (atom.GetIdx() not in originally_no_implicit
                            and atom.GetSymbol() not in _METAL_SET):
                        atom.SetNoImplicit(False)
                        atom.SetNumExplicitHs(0)
                mol.UpdatePropertyCache(strict=False)
                mol = Chem.AddHs(mol, addCoords=True)
                mol = _fix_zero_coord_hydrogens(mol)
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


def _fix_zero_coord_hydrogens(mol):
    """Fix H atoms stuck at (0,0,0) after AddHs(addCoords=True).

    When RDKit cannot compute coordinates for some H atoms, it places them
    at the origin.  This function detects those H atoms and places them at
    a reasonable position (~1.09 A from the parent atom) using the parent's
    existing neighbors to determine the correct direction.
    """
    if mol.GetNumConformers() == 0:
        return mol
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H':
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        if abs(pos.x) > 0.01 or abs(pos.y) > 0.01 or abs(pos.z) > 0.01:
            continue
        # This H is at (0,0,0) - fix it
        parent = atom.GetNeighbors()[0]
        p_pos = conf.GetAtomPosition(parent.GetIdx())
        # Compute direction away from other neighbors
        import numpy as np
        p = np.array([p_pos.x, p_pos.y, p_pos.z])
        neighbor_vecs = []
        for nbr in parent.GetNeighbors():
            if nbr.GetIdx() == atom.GetIdx():
                continue
            n_pos = conf.GetAtomPosition(nbr.GetIdx())
            v = np.array([n_pos.x, n_pos.y, n_pos.z]) - p
            norm = np.linalg.norm(v)
            if norm > 0.01:
                neighbor_vecs.append(v / norm)
        if neighbor_vecs:
            # Place H opposite to the average neighbor direction
            avg_dir = np.mean(neighbor_vecs, axis=0)
            norm = np.linalg.norm(avg_dir)
            if norm > 0.01:
                h_dir = -avg_dir / norm
            else:
                h_dir = np.array([0.0, 0.0, 1.0])
        else:
            h_dir = np.array([0.0, 0.0, 1.0])
        h_pos = p + h_dir * 1.09  # C-H bond length
        from rdkit.Geometry import Point3D
        conf.SetAtomPosition(atom.GetIdx(), Point3D(*h_pos.tolist()))
    return mol


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


def _build_uff_constraints_from_template(
    mol_template,
    xyz_delfin: Optional[str] = None,
) -> Optional[Dict]:
    """Build conservative OB-UFF constraints from the RDKit template graph.

    The goal is not to freeze the structure, but to keep fragile motifs
    chemically reasonable during UFF relaxation:
    - carboxyl-like O-C-O groups near trigonal-planar geometry
    - aromatic ring torsions near planarity
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return None

    constraints: Dict[str, list] = {
        "fix_atoms": [],
        "distances": [],
        "angles": [],
        "torsions": [],
    }
    seen_dist: set = set()
    seen_angle: set = set()
    seen_tors: set = set()

    # Optional source coordinates (same atom ordering as mol_template) to pick
    # the closest planar torsion target (0 or 180 deg).
    coords: Optional[List[Tuple[float, float, float]]] = None
    if xyz_delfin:
        try:
            lines = [l for l in xyz_delfin.splitlines() if l.strip()]
            if len(lines) == mol_template.GetNumAtoms():
                parsed: List[Tuple[float, float, float]] = []
                for line in lines:
                    parts = line.split()
                    if len(parts) < 4:
                        parsed = []
                        break
                    parsed.append((float(parts[1]), float(parts[2]), float(parts[3])))
                if len(parsed) == mol_template.GetNumAtoms():
                    coords = parsed
        except Exception:
            coords = None

    def _nearest_planar_target(a: int, b: int, c: int, d: int) -> float:
        if coords is None:
            return 0.0
        try:
            p1 = coords[a]
            p2 = coords[b]
            p3 = coords[c]
            p4 = coords[d]
            b1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
            b2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
            b3 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])

            n1 = (
                b1[1] * b2[2] - b1[2] * b2[1],
                b1[2] * b2[0] - b1[0] * b2[2],
                b1[0] * b2[1] - b1[1] * b2[0],
            )
            n2 = (
                b2[1] * b3[2] - b2[2] * b3[1],
                b2[2] * b3[0] - b2[0] * b3[2],
                b2[0] * b3[1] - b2[1] * b3[0],
            )
            n1m = math.sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2])
            n2m = math.sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2])
            if n1m < 1e-10 or n2m < 1e-10:
                return 0.0
            n1u = (n1[0] / n1m, n1[1] / n1m, n1[2] / n1m)
            n2u = (n2[0] / n2m, n2[1] / n2m, n2[2] / n2m)
            dot = max(-1.0, min(1.0, n1u[0] * n2u[0] + n1u[1] * n2u[1] + n1u[2] * n2u[2]))
            dihedral = math.degrees(math.acos(dot))
            # Planar torsions are near 0 or 180.
            return 180.0 if abs(dihedral - 180.0) < abs(dihedral - 0.0) else 0.0
        except Exception:
            return 0.0

    try:
        # Carboxyl/carbonyl-like O-C-O units: keep broad C-O distances and ~120 deg angle.
        for atom in mol_template.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetSymbol() in _METAL_SET:
                continue
            o_neighbors = [
                n for n in atom.GetNeighbors()
                if n.GetAtomicNum() == 8 and n.GetSymbol() not in _METAL_SET
            ]
            if len(o_neighbors) != 2:
                continue

            c_idx = atom.GetIdx()
            o1_idx, o2_idx = sorted([o_neighbors[0].GetIdx(), o_neighbors[1].GetIdx()])

            angle_key = (o1_idx, c_idx, o2_idx)
            if angle_key not in seen_angle:
                seen_angle.add(angle_key)
                constraints["angles"].append((o1_idx, c_idx, o2_idx, 120.0))

            for o_idx in (o1_idx, o2_idx):
                dist_key = tuple(sorted((c_idx, o_idx)))
                if dist_key in seen_dist:
                    continue
                seen_dist.add(dist_key)
                bond = mol_template.GetBondBetweenAtoms(c_idx, o_idx)
                target = 1.27
                if bond is not None:
                    btype = bond.GetBondType()
                    if btype == Chem.BondType.DOUBLE:
                        target = 1.22
                    elif btype == Chem.BondType.SINGLE:
                        target = 1.31
                    elif btype == Chem.BondType.AROMATIC:
                        target = 1.27
                constraints["distances"].append((c_idx, o_idx, target))

            # Keep carboxyl group coplanar with its carbon backbone when possible.
            x_neighbors = sorted(
                n.GetIdx()
                for n in atom.GetNeighbors()
                if n.GetIdx() not in (o1_idx, o2_idx)
                and n.GetAtomicNum() > 1
                and n.GetSymbol() not in _METAL_SET
            )
            for x_idx in x_neighbors:
                x_atom = mol_template.GetAtomWithIdx(x_idx)
                anchor_candidates = sorted(
                    n.GetIdx()
                    for n in x_atom.GetNeighbors()
                    if n.GetIdx() != c_idx
                    and n.GetAtomicNum() > 1
                    and n.GetSymbol() not in _METAL_SET
                )
                if not anchor_candidates:
                    continue
                a_idx = anchor_candidates[0]
                for o_idx in (o1_idx, o2_idx):
                    key = (a_idx, x_idx, c_idx, o_idx)
                    rev = (o_idx, c_idx, x_idx, a_idx)
                    tors_key = key if key <= rev else rev
                    if tors_key in seen_tors:
                        continue
                    seen_tors.add(tors_key)
                    constraints["torsions"].append(
                        (a_idx, x_idx, c_idx, o_idx, _nearest_planar_target(a_idx, x_idx, c_idx, o_idx))
                    )

        # Aromatic ring planarity via torsion constraints on consecutive ring quartets.
        try:
            Chem.GetSymmSSSR(mol_template)
        except Exception:
            pass
        ring_info = mol_template.GetRingInfo()
        if ring_info is not None:
            for ring in ring_info.AtomRings():
                if len(ring) < 5 or len(ring) > 7:
                    continue
                if not all(mol_template.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    continue
                if any(mol_template.GetAtomWithIdx(i).GetSymbol() in _METAL_SET for i in ring):
                    continue

                n = len(ring)
                for i in range(n):
                    a = ring[(i - 1) % n]
                    b = ring[i]
                    c = ring[(i + 1) % n]
                    d = ring[(i + 2) % n]
                    if len({a, b, c, d}) < 4:
                        continue
                    key = (a, b, c, d)
                    rev = (d, c, b, a)
                    tors_key = key if key <= rev else rev
                    if tors_key in seen_tors:
                        continue
                    seen_tors.add(tors_key)
                    constraints["torsions"].append((a, b, c, d, _nearest_planar_target(a, b, c, d)))

        # Metal-metal bond distance constraints.
        for bond in mol_template.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetSymbol() not in _METAL_SET or a2.GetSymbol() not in _METAL_SET:
                continue
            i1, i2 = a1.GetIdx(), a2.GetIdx()
            dist_key = tuple(sorted((i1, i2)))
            if dist_key in seen_dist:
                continue
            seen_dist.add(dist_key)
            mm_key = frozenset({a1.GetSymbol(), a2.GetSymbol()})
            target = _METAL_METAL_BOND_LENGTHS.get(mm_key)
            if target is None:
                r1 = _COVALENT_RADII.get(a1.GetSymbol())
                r2 = _COVALENT_RADII.get(a2.GetSymbol())
                if r1 is not None and r2 is not None:
                    target = r1 + r2 + 0.3
                else:
                    target = 2.5
            constraints["distances"].append((i1, i2, float(target)))
    except Exception:
        return None

    if not (
        constraints["fix_atoms"]
        or constraints["distances"]
        or constraints["angles"]
        or constraints["torsions"]
    ):
        return None
    return constraints


def _build_coordination_constraints_from_xyz(
    mol_template,
    xyz_delfin: str,
) -> Optional[Dict]:
    """Auto-detect metal coordination from template graph and pin it during UFF.

    Unlike ``_build_coordination_uff_constraints`` (which needs explicit
    perm/geometry arguments from the topology enumerator), this function
    derives constraints directly from the template bond graph.  This
    makes it usable for ANY UFF call — sampling, linkage, alt-binding,
    hapto, etc.

    Only M-D DISTANCES are constrained (from the lookup table).  L-M-L
    angles are left free so UFF can optimize toward ideal polyhedral
    geometry.  This combines Principle 2 (distances from chemistry) with
    Principle 3 (UFF improves angles toward ideal).
    """
    base = _build_uff_constraints_from_template(mol_template, xyz_delfin=xyz_delfin)
    if base is None:
        base = {"fix_atoms": [], "distances": [], "angles": [], "torsions": []}

    if not RDKIT_AVAILABLE or mol_template is None or not xyz_delfin:
        return base if any(base.values()) else None

    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if len(lines) != mol_template.GetNumAtoms():
            return base if any(base.values()) else None
        coords: List[Tuple[float, float, float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return base if any(base.values()) else None
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

        seen_dist = {tuple(sorted((a, b))) for a, b, _t in base["distances"]}
        seen_angle = set()
        for a, b, c, _t in base["angles"]:
            seen_angle.add((a, b, c))
            seen_angle.add((c, b, a))

        for atom in mol_template.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            m_idx = atom.GetIdx()
            m_sym = atom.GetSymbol()

            donor_indices: List[int] = []
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() <= 1:
                    continue
                donor_indices.append(nbr.GetIdx())

            # M-D distance constraints from lookup table.
            for d_idx in donor_indices:
                key = tuple(sorted((m_idx, d_idx)))
                if key in seen_dist:
                    continue
                seen_dist.add(key)
                d_sym = mol_template.GetAtomWithIdx(d_idx).GetSymbol()
                if d_sym in _METAL_SET:
                    continue
                bl = float(_get_ml_bond_length(m_sym, d_sym))
                base["distances"].append((m_idx, d_idx, bl))

            # L-M-L angles are NOT constrained here — UFF should be free
            # to optimize angles toward the ideal polyhedron.

        # Hapto protection: fix hapto metal atoms + all Cp ring atoms
        # during UFF. Without metal FF parameters, UFF treats the metal
        # as empty space, allowing ring atoms to drift out of plane.
        # Fixing these atoms preserves the hapto builder's geometry.
        try:
            hapto_groups = _find_hapto_groups(mol_template)
            fixed = set()
            for _hm, members in hapto_groups:
                fixed.add(_hm)  # fix the hapto metal
                for ci in members:
                    fixed.add(ci)  # fix all ring atoms
            for fi in sorted(fixed):
                if fi not in base["fix_atoms"]:
                    base["fix_atoms"].append(fi)
        except Exception:
            pass

    except Exception as exc:
        logger.debug("_build_coordination_constraints_from_xyz failed: %s", exc)

    return base if any(base.values()) else None


def _build_coordination_uff_constraints(
    mol_template,
    metal_idx: int,
    donor_atom_indices: List[int],
    perm: List[int],
    geometry: str,
    xyz_delfin: Optional[str] = None,
) -> Optional[Dict]:
    """UFF constraints that preserve an idealized coordination polyhedron.

    Extends :func:`_build_uff_constraints_from_template` with hard
    metal-donor distance pins (from :func:`_get_ml_bond_length`) and
    donor-metal-donor angle targets derived from
    ``_TOPO_GEOMETRY_VECTORS[geometry]`` combined with ``perm``.

    These constraints keep octahedral/PBP/etc. geometry intact during
    UFF relaxation even when UFF lacks parameters for the metal.
    """
    base = _build_uff_constraints_from_template(
        mol_template, xyz_delfin=xyz_delfin
    )
    if base is None:
        base = {
            "fix_atoms": [],
            "distances": [],
            "angles": [],
            "torsions": [],
        }

    vectors = _TOPO_GEOMETRY_VECTORS.get(geometry)
    if not vectors or metal_idx is None:
        return base if any(base.values()) else None

    try:
        metal_sym = mol_template.GetAtomWithIdx(metal_idx).GetSymbol()
    except Exception:
        return base if any(base.values()) else None

    n_pos = len(vectors)
    if n_pos != len(perm) or n_pos != len(donor_atom_indices):
        return base if any(base.values()) else None

    # Unit vectors for each geometry position.
    unit_vecs: List[Tuple[float, float, float]] = []
    for vx, vy, vz in vectors:
        m = math.sqrt(vx * vx + vy * vy + vz * vz)
        if m < 1e-8:
            unit_vecs.append((0.0, 0.0, 0.0))
        else:
            unit_vecs.append((vx / m, vy / m, vz / m))

    # Track existing distance keys so we don't duplicate OCO pins.
    seen_dist = {tuple(sorted((a, b))) for a, b, _t in base["distances"]}

    # M-D distance pins.
    donor_pos_map: Dict[int, int] = {}
    for pos_idx, donor_list_idx in enumerate(perm):
        if donor_list_idx < 0 or donor_list_idx >= len(donor_atom_indices):
            continue
        donor_atom_idx = donor_atom_indices[donor_list_idx]
        donor_pos_map[donor_atom_idx] = pos_idx
        try:
            donor_sym = mol_template.GetAtomWithIdx(donor_atom_idx).GetSymbol()
        except Exception:
            continue
        bl = float(_get_ml_bond_length(metal_sym, donor_sym))
        key = tuple(sorted((metal_idx, donor_atom_idx)))
        if key in seen_dist:
            continue
        seen_dist.add(key)
        base["distances"].append((metal_idx, donor_atom_idx, bl))

    # L-M-L angle targets from ideal geometry vectors.
    seen_angle = {
        (a, b, c) for a, b, c, _t in base["angles"]
    }
    seen_angle.update(
        (c, b, a) for a, b, c, _t in base["angles"]
    )
    donors_with_pos = [
        d for d in donor_atom_indices if d in donor_pos_map
    ]
    for i in range(len(donors_with_pos)):
        for j in range(i + 1, len(donors_with_pos)):
            d_i = donors_with_pos[i]
            d_j = donors_with_pos[j]
            pi = donor_pos_map[d_i]
            pj = donor_pos_map[d_j]
            vi = unit_vecs[pi]
            vj = unit_vecs[pj]
            dot = vi[0] * vj[0] + vi[1] * vj[1] + vi[2] * vj[2]
            dot = max(-1.0, min(1.0, dot))
            ideal = math.degrees(math.acos(dot))
            key = (d_i, metal_idx, d_j)
            if key in seen_angle:
                continue
            seen_angle.add(key)
            seen_angle.add((d_j, metal_idx, d_i))
            base["angles"].append((d_i, metal_idx, d_j, ideal))

    if not (
        base["fix_atoms"]
        or base["distances"]
        or base["angles"]
        or base["torsions"]
    ):
        return None
    return base


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
            - ``'torsions'``: list of ``(idx_a, idx_b, idx_c, idx_d, target_deg)``

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
                # Torsion constraints
                for idx_a, idx_b, idx_c, idx_d, target in constraints.get('torsions', []):
                    ob_constraints.AddTorsionConstraint(
                        idx_a + 1, idx_b + 1, idx_c + 1, idx_d + 1, target
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


def _optimize_xyz_openbabel_safe(
    xyz_delfin: str,
    mol_template=None,
    smiles: Optional[str] = None,
    steps: int = 500,
    apply_template_constraints: bool = False,
    coord_constraints: Optional[Dict] = None,
) -> str:
    """Run OB-UFF, but keep original XYZ if optimization breaks topology.

    OB force fields operate on perceived connectivity from XYZ and can
    occasionally stretch/break covalent bonds for charged metal complexes.
    This wrapper accepts the optimized geometry only if it passes structural
    sanity checks against the original molecular graph.

    When ``coord_constraints`` is provided, those constraints are used
    directly (topology-enumerator path with M-D + L-M-L pinning).
    Otherwise, when ``apply_template_constraints`` is True and a
    ``mol_template`` is available, mild OCO/aromatic planarity constraints
    are passed to OB-UFF.
    """
    constraints = coord_constraints
    if constraints is None and mol_template is not None:
        # Universal principle: EVERY UFF call on a metal complex gets
        # coordination constraints (M-D distances + L-M-L angles).
        # This prevents UFF from inventing distances for metals it has
        # no parameters for (Sc, Cd, lanthanides, etc.).
        _has_metal_in_template = any(
            a.GetSymbol() in _METAL_SET for a in mol_template.GetAtoms()
        ) if RDKIT_AVAILABLE else False
        if _has_metal_in_template:
            try:
                constraints = _build_coordination_constraints_from_xyz(
                    mol_template, xyz_delfin
                )
            except Exception as exc:
                logger.debug("Coordination constraint generation failed: %s", exc)
                constraints = None
        elif apply_template_constraints:
            try:
                constraints = _build_uff_constraints_from_template(
                    mol_template, xyz_delfin=xyz_delfin
                )
            except Exception as exc:
                logger.debug("Template constraint generation failed: %s", exc)
                constraints = None

    xyz_opt = _optimize_xyz_openbabel(xyz_delfin, steps=steps, constraints=constraints)
    if not xyz_opt or xyz_opt == xyz_delfin:
        return xyz_delfin

    # Fundamental check: UFF must not change the metal-donor connectivity.
    # If a donor drifted away or a non-donor collapsed onto the metal,
    # discard the UFF result and keep the original XYZ.
    if mol_template is not None and RDKIT_AVAILABLE:
        if not _verify_metal_connectivity(xyz_opt, mol_template):
            logger.debug("Discarding UFF geometry: metal connectivity changed.")
            return xyz_delfin

    if mol_template is not None and RDKIT_AVAILABLE:
        try:
            # Build baseline quality from the pre-UFF geometry so we can keep
            # constrained UFF structures that are improved even if still not
            # fully "good" by strict thresholds.
            orig_bad = False
            orig_score = float("inf")
            mol_orig = Chem.RWMol(mol_template)
            mol_orig.RemoveAllConformers()
            conf_orig = _xyz_to_rdkit_conformer(mol_orig.GetMol(), xyz_delfin)
            if conf_orig is not None:
                cid_orig = mol_orig.AddConformer(conf_orig, assignId=True)
                orig_bad = _has_bad_geometry(mol_orig.GetMol(), cid_orig)
                try:
                    orig_score = _geometry_quality_score(mol_orig.GetMol(), cid_orig)
                except Exception:
                    orig_score = float("inf")

            mol_tmp = Chem.RWMol(mol_template)
            mol_tmp.RemoveAllConformers()
            conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_opt)
            if conf is None:
                logger.debug("Discarding UFF geometry: atom mapping failed.")
                return xyz_delfin
            cid = mol_tmp.AddConformer(conf, assignId=True)
            if _has_severe_covalent_distortion(mol_tmp.GetMol(), cid):
                logger.debug("Discarding UFF geometry: severe covalent distortion.")
                return xyz_delfin
            opt_bad = _has_bad_geometry(mol_tmp.GetMol(), cid)
            if opt_bad:
                # Keep constrained UFF if it improves relative to the original
                # geometry; otherwise fall back to the unoptimized structure.
                try:
                    opt_score = _geometry_quality_score(mol_tmp.GetMol(), cid)
                except Exception:
                    opt_score = float("inf")
                if (not orig_bad) or (opt_score >= orig_score - 1e-6):
                    logger.debug(
                        "Discarding UFF geometry: no bad-geometry improvement "
                        "(orig_bad=%s orig_score=%.3f opt_score=%.3f).",
                        orig_bad, orig_score, opt_score,
                    )
                    return xyz_delfin
        except Exception as e:
            logger.debug("Discarding UFF geometry: post-check failed (%s).", e)
            return xyz_delfin

    # For constrained, template-guided UFF in the isomer path we defer
    # topology/spurious checks to the caller (which applies the same checks
    # afterwards) to avoid double-rejecting all optimized geometries.
    if smiles and not apply_template_constraints and coord_constraints is None:
        try:
            if not _roundtrip_ring_count_ok(xyz_opt, smiles):
                logger.debug("Discarding UFF geometry: ring-count mismatch.")
                return xyz_delfin
            if not _no_spurious_bonds(xyz_opt, smiles):
                logger.debug("Discarding UFF geometry: spurious bonds.")
                return xyz_delfin
            if not _fragment_topology_ok(xyz_opt, smiles):
                logger.debug("Discarding UFF geometry: fragment topology mismatch.")
                return xyz_delfin
        except Exception:
            pass

    return xyz_opt


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


def _smiles_to_architector_input(smiles: str) -> Optional[dict]:
    """Decompose a metal-complex SMILES into an Architector input dict."""
    if not RDKIT_AVAILABLE:
        return None

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None

    metal_indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METALS:
            metal_indices.append(atom.GetIdx())
    if not metal_indices:
        return None

    metal_idx = metal_indices[0]
    metal_atom = mol.GetAtomWithIdx(metal_idx)
    metal_symbol = metal_atom.GetSymbol()
    metal_charge = metal_atom.GetFormalCharge()

    coord_atom_to_metal = {}
    bonds_to_remove = []
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in metal_indices:
            coord_atom_to_metal[b] = a
            bonds_to_remove.append((a, b))
        elif b in metal_indices:
            coord_atom_to_metal[a] = b
            bonds_to_remove.append((a, b))

    edit_mol = Chem.RWMol(mol)
    orig_props = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        orig_props[idx] = {
            'formal_charge': atom.GetFormalCharge(),
            'explicit_h': atom.GetNumExplicitHs(),
            'no_implicit': atom.GetNoImplicit(),
            'rad_e': atom.GetNumRadicalElectrons(),
        }

    for a, b in bonds_to_remove:
        edit_mol.RemoveBond(a, b)

    metal_set = set(metal_indices)
    idx_map = {}
    removed = 0
    for old_idx in range(mol.GetNumAtoms()):
        if old_idx in metal_set:
            removed += 1
        else:
            idx_map[old_idx] = old_idx - removed

    for mi in sorted(metal_indices, reverse=True):
        edit_mol.RemoveAtom(mi)

    mol_no_metal = edit_mol.GetMol()

    for old_idx, props in orig_props.items():
        if old_idx in metal_set:
            continue
        new_idx = idx_map.get(old_idx)
        if new_idx is None:
            continue
        try:
            atom = mol_no_metal.GetAtomWithIdx(new_idx)
            atom.SetNumRadicalElectrons(props['rad_e'])
            if old_idx in coord_atom_to_metal:
                typical_valence = {
                    'C': 4, 'N': 3, 'O': 2, 'S': 2, 'Se': 2,
                    'P': 3, 'As': 3, 'Si': 4, 'B': 3, 'Te': 2,
                    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
                }
                sym = atom.GetSymbol()
                typ_val = typical_valence.get(sym, 2)
                bo_sum = 0.0
                for bond in mol.GetBonds():
                    ba, bb = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    if old_idx in (ba, bb):
                        neighbor = bb if ba == old_idx else ba
                        if neighbor not in metal_set:
                            bo_sum += bond.GetBondTypeAsDouble()
                saturated = (bo_sum + props['explicit_h']) >= typ_val

                if props['formal_charge'] != 0:
                    atom.SetFormalCharge(props['formal_charge'])
                    atom.SetNumExplicitHs(0)
                    atom.SetNoImplicit(True)
                elif saturated:
                    atom.SetFormalCharge(0)
                else:
                    atom.SetFormalCharge(-1)
                    atom.SetNumExplicitHs(0)
                    atom.SetNoImplicit(True)
            else:
                atom.SetNumExplicitHs(props['explicit_h'])
                atom.SetNoImplicit(props['no_implicit'])
        except Exception:
            pass

    try:
        mol_no_metal.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    frag_atom_lists = Chem.GetMolFrags(mol_no_metal, asMols=False)
    frag_mols = Chem.GetMolFrags(mol_no_metal, asMols=True, sanitizeFrags=False)
    new_to_old = {v: k for k, v in idx_map.items()}

    ligands = []
    for frag_atoms, frag_mol in zip(frag_atom_lists, frag_mols):
        coord_positions = []
        for pos_in_frag, new_idx in enumerate(frag_atoms):
            old_idx = new_to_old.get(new_idx)
            if old_idx is not None and old_idx in coord_atom_to_metal:
                coord_positions.append(pos_in_frag)

        for atom in frag_mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)

        mapped_smiles = Chem.MolToSmiles(frag_mol, canonical=True)
        mapped_mol = Chem.MolFromSmiles(mapped_smiles, sanitize=False)
        if mapped_mol is None:
            continue

        frag_to_canon = {}
        for atom in mapped_mol.GetAtoms():
            orig_frag_idx = atom.GetAtomMapNum() - 1
            frag_to_canon[orig_frag_idx] = atom.GetIdx()

        canon_coord = []
        for pos in coord_positions:
            mapped_idx = frag_to_canon.get(pos)
            if mapped_idx is not None:
                canon_coord.append(mapped_idx)

        if not canon_coord:
            continue

        for atom in mapped_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        canon_smiles = Chem.MolToSmiles(mapped_mol, canonical=False)
        ligands.append({'smiles': canon_smiles, 'coordList': sorted(canon_coord)})

    if not ligands:
        return None

    total_charge = metal_charge
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in metal_set:
            total_charge += atom.GetFormalCharge()

    return {
        'core': {'metal': metal_symbol},
        'ligands': ligands,
        'parameters': {'full_charge': total_charge},
    }


def smiles_to_xyz_architector(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Convert a metal-complex SMILES to XYZ using Architector."""
    if not contains_metal(smiles):
        return None, 'Architector requires a metal-containing SMILES string'

    try:
        import importlib.util
        if importlib.util.find_spec('architector') is None:
            return None, 'architector is not installed'

        from architector.complex_construction import build_complex_driver
        from architector.io_process_input import inparse
        from architector import io_ptable
    except Exception as exc:
        return None, f'Could not import architector: {exc}'

    try:
        input_dict = _smiles_to_architector_input(smiles)
        if input_dict is None:
            return None, 'Could not decompose SMILES into metal + ligands for Architector'

        input_dict = inparse(input_dict)
        results = build_complex_driver(input_dict)

        real_keys = [k for k in results if '_init_only' not in k]
        ligands = input_dict.get('ligands', [])
        max_dent = max((len(l.get('coordList', [])) for l in ligands), default=0)
        if not real_keys and max_dent >= 2:
            for larger in (True, False):
                scaled = io_ptable.map_metal_radii(
                    inparse(input_dict), larger=larger,
                )
                extra = build_complex_driver(scaled)
                suffix = '_larger_scaled' if larger else '_smaller_scaled'
                for key, value in extra.items():
                    results[key + suffix] = value
                if any('_init_only' not in key for key in extra):
                    break

        ranked = []
        for key, mol in results.items():
            if '_init_only' in key:
                continue
            atoms = mol.get('ase_atoms')
            if atoms is None:
                continue
            energy = mol.get('energy', None)
            ranked.append((float(energy) if energy is not None else float('inf'), atoms, key))

        if not ranked:
            return None, 'Architector returned no valid 3D structures'

        ranked.sort(key=lambda item: item[0])
        atoms = ranked[0][1]
        xyz_lines = []
        for atom, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            xyz_lines.append(f'{atom}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}')
        return '\n'.join(xyz_lines), None
    except Exception as exc:
        return None, f'Architector conversion failed: {exc}'
