"""Element sets, covalent radii, metal-ligand / metal-metal / metal-centroid bond-length tables, preferred coordination polyhedra per metal and the polyhedron catalogue of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
import re
import threading
from typing import Dict, FrozenSet, List, Optional, Tuple

from delfin.common.logging import get_logger

logger = get_logger("delfin.smiles_converter")


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


# stk is optional, but only in one sense: NOT INSTALLED. Any other import failure (a broken
# install, or a shared-library conflict such as an older libstdc++ loaded first from another
# program's directory on LD_LIBRARY_PATH) is raised, never swallowed: a silent fallback would make
# the construction depend on the order in which a process happened to import its modules.
try:
    import stk
    STK_AVAILABLE = True
except ModuleNotFoundError as _stk_err:
    if _stk_err.name != "stk":
        raise
    stk = None
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
    # Less common but legitimate coordinating main-group atoms
    # (boryl, silyl, germyl, arsine, stiboranyl, etc.).  Cordero 2008.
    'B': 0.84, 'Si': 1.11, 'Ge': 1.20, 'As': 1.19, 'Sb': 1.39,
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
    # Iridium — CSD averages for Ir(III) octahedral coordination.
    # Ir-C and Ir-N are tuned to the cyclometallated regime (Ir(ppy)2
    # type, 2.01/2.04 Å) because pure σ-amine Ir-N is very rare in
    # output targets; Ir-O is the acac / β-diketonate range.
    ('Ir', 'C'): 2.02, ('Ir', 'N'): 2.05, ('Ir', 'O'): 2.05,
    ('Ir', 'P'): 2.30, ('Ir', 'Cl'): 2.38, ('Ir', 'S'): 2.38,
    # Ruthenium — CSD averages for Ru(II) polypyridyl / cyclometallated
    # complexes: Ru-C(Ph-pyridyl) ~2.02, Ru-N(bpy/ppy) ~2.05, Ru-O(acac) ~2.05.
    ('Ru', 'C'): 2.02, ('Ru', 'N'): 2.06, ('Ru', 'O'): 2.06,
    ('Ru', 'P'): 2.30, ('Ru', 'Cl'): 2.40, ('Ru', 'S'): 2.35,
    # Rhodium
    # Rhodium — CSD averages (cyclometallated Rh(ppy)2 type: Rh-C 1.99,
    # Rh-N(ppy trans) 2.04; reference [Rh(ppy)2Cl(CH3CN)]).
    ('Rh', 'C'): 2.00, ('Rh', 'N'): 2.05, ('Rh', 'O'): 2.05,
    ('Rh', 'P'): 2.28, ('Rh', 'Cl'): 2.38,
    # Platinum — cyclometallated Pt(II) Pt-C ~2.01, Pt-N ~2.02.
    ('Pt', 'C'): 2.01, ('Pt', 'N'): 2.03, ('Pt', 'O'): 2.02,
    ('Pt', 'P'): 2.25, ('Pt', 'Cl'): 2.30, ('Pt', 'S'): 2.30,
    ('Pt', 'Br'): 2.43,
    # Palladium — cyclometallated Pd(II) Pd-C ~1.99, Pd-N ~2.03 (Δ ~0.04 Å
    # vs. longer values seen in non-cyclometallated amines).
    ('Pd', 'C'): 1.99, ('Pd', 'N'): 2.03, ('Pd', 'O'): 2.02,
    ('Pd', 'P'): 2.28, ('Pd', 'Cl'): 2.30, ('Pd', 'S'): 2.30,
    # Gold
    ('Au', 'C'): 2.00, ('Au', 'N'): 2.05, ('Au', 'P'): 2.28,
    ('Au', 'Cl'): 2.28, ('Au', 'S'): 2.30,
    # Iron — Fe-N(py/imine) 2.15 HS / 1.97 LS; 2.05 is the CSD average
    # across coordination environments (low-spin carbene / high-spin py).
    ('Fe', 'C'): 1.95, ('Fe', 'N'): 2.05, ('Fe', 'O'): 2.00,
    ('Fe', 'P'): 2.20, ('Fe', 'Cl'): 2.28, ('Fe', 'S'): 2.30,
    # Cobalt — Co(III) octahedral Co-N 1.97 LS / Co(II) 2.15 HS;
    # 2.00 is a reasonable single-value average.
    ('Co', 'C'): 1.95, ('Co', 'N'): 2.00, ('Co', 'O'): 1.95,
    ('Co', 'P'): 2.18, ('Co', 'Cl'): 2.25,
    # Nickel — Ni(II)-N(pyridyl/amine) ~2.05 octahedral, 1.90 square-planar
    # low-spin; 2.00 covers both regimes.
    ('Ni', 'C'): 1.95, ('Ni', 'N'): 2.00, ('Ni', 'O'): 1.95,
    ('Ni', 'P'): 2.18, ('Ni', 'Cl'): 2.25,
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
    ('Co', 'S'): 2.25, ('Co', 'Br'): 2.38, ('Co', 'F'): 1.88,
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
    # Metal-Iodide — CSD averages.  The covalent-radius fallback
    # systematically overshoots M-I because the Pyykkö 2009 radii
    # under-estimate the iodide contraction in ionic bonds, leaving
    # real terminal M-I pairs at the fallback's 3.1-3.4 Å instead of
    # the experimental 2.6-2.9 Å.  Explicit entries below fix that.
    ('Sc', 'I'): 2.90, ('Ti', 'I'): 2.79, ('V',  'I'): 2.76,
    ('Cr', 'I'): 2.70, ('Mn', 'I'): 2.70, ('Fe', 'I'): 2.57,
    ('Co', 'I'): 2.56, ('Ni', 'I'): 2.53, ('Cu', 'I'): 2.57,
    ('Zn', 'I'): 2.64,
    ('Y',  'I'): 2.95, ('Zr', 'I'): 2.85, ('Nb', 'I'): 2.85,
    ('Mo', 'I'): 2.75, ('Tc', 'I'): 2.75, ('Ru', 'I'): 2.70,
    ('Rh', 'I'): 2.68, ('Pd', 'I'): 2.62, ('Ag', 'I'): 2.80,
    ('Cd', 'I'): 2.75,
    ('Hf', 'I'): 2.85, ('Ta', 'I'): 2.80, ('W',  'I'): 2.80,
    ('Re', 'I'): 2.75, ('Os', 'I'): 2.72, ('Ir', 'I'): 2.67,
    ('Pt', 'I'): 2.60, ('Au', 'I'): 2.58, ('Hg', 'I'): 2.76,
    ('La', 'I'): 3.10, ('Ce', 'I'): 3.05, ('Gd', 'I'): 3.00,
    ('Lu', 'I'): 2.90,
    ('Al', 'I'): 2.50, ('Ga', 'I'): 2.55, ('In', 'I'): 2.70,
    ('Sn', 'I'): 2.75, ('Pb', 'I'): 2.90, ('Bi', 'I'): 2.85,
    ('U',  'I'): 3.05, ('Th', 'I'): 3.05,
    ('Li', 'I'): 2.55, ('Na', 'I'): 2.85, ('K',  'I'): 3.20,
    # Metal–Silyl bonds.  The covalent-radius fallback gives ~2.75-3.02 Å
    # which is 0.4-0.6 Å too long versus CSD averages for σ-silyl ligands
    # (M-SiR3 is shorter than the additive-radii estimate because Si
    # acts as a soft σ-donor with significant π-backbonding).  Without
    # these entries, σ-Si donors on hapto metals end up far outside the
    # bonding range during topology check (observed 3.08-3.45 Å).
    ('Sc', 'Si'): 2.55, ('Ti', 'Si'): 2.50, ('V',  'Si'): 2.45,
    ('Cr', 'Si'): 2.42, ('Mn', 'Si'): 2.40, ('Fe', 'Si'): 2.35,
    ('Co', 'Si'): 2.30, ('Ni', 'Si'): 2.25, ('Cu', 'Si'): 2.30,
    ('Zn', 'Si'): 2.40,
    ('Y',  'Si'): 2.85, ('Zr', 'Si'): 2.75,
    ('Nb', 'Si'): 2.65, ('Mo', 'Si'): 2.55, ('Tc', 'Si'): 2.50,
    ('Ru', 'Si'): 2.40, ('Rh', 'Si'): 2.30, ('Pd', 'Si'): 2.30,
    ('Ag', 'Si'): 2.50, ('Cd', 'Si'): 2.55,
    ('Hf', 'Si'): 2.75, ('Ta', 'Si'): 2.65, ('W',  'Si'): 2.55,
    ('Re', 'Si'): 2.50, ('Os', 'Si'): 2.40, ('Ir', 'Si'): 2.40,
    ('Pt', 'Si'): 2.35, ('Au', 'Si'): 2.40, ('Hg', 'Si'): 2.50,
    # -------------------------------------------------------------------
    # Comprehensive fill-in for all remaining (M, L) pairs of practical
    # interest (L ∈ {C, N, O, P, S, F, Cl, Br, Se}).  Values are CSD /
    # Cordero 2008 averages with minor adjustments for common
    # coordination environments; a missing entry here means the
    # covalent-radius fallback is already acceptable (< 5 % deviation).
    # -------------------------------------------------------------------
    # Alkali completions
    ('Li', 'C'): 2.15, ('Li', 'P'): 2.50, ('Li', 'S'): 2.40,
    ('Li', 'F'): 1.80, ('Li', 'Br'): 2.60, ('Li', 'Se'): 2.55,
    ('Na', 'C'): 2.70, ('Na', 'P'): 2.85, ('Na', 'S'): 2.75,
    ('Na', 'F'): 2.20, ('Na', 'Br'): 2.95, ('Na', 'Se'): 2.90,
    ('K',  'C'): 3.10, ('K',  'P'): 3.25, ('K',  'S'): 3.15,
    ('K',  'F'): 2.55, ('K',  'Br'): 3.30, ('K',  'Se'): 3.25,
    ('Rb', 'N'): 2.95, ('Rb', 'O'): 2.85, ('Rb', 'Cl'): 3.20,
    ('Rb', 'Br'): 3.40, ('Rb', 'I'): 3.55, ('Rb', 'F'): 2.80,
    ('Cs', 'N'): 3.10, ('Cs', 'O'): 3.00, ('Cs', 'Cl'): 3.40,
    ('Cs', 'Br'): 3.60, ('Cs', 'I'): 3.75, ('Cs', 'F'): 3.00,
    # Alkaline earth completions
    ('Be', 'N'): 1.75, ('Be', 'O'): 1.65, ('Be', 'Cl'): 1.85,
    ('Be', 'Br'): 2.00, ('Be', 'I'): 2.20, ('Be', 'F'): 1.55,
    ('Be', 'C'): 1.70, ('Be', 'P'): 2.05, ('Be', 'S'): 1.95,
    ('Mg', 'C'): 2.15, ('Mg', 'P'): 2.45, ('Mg', 'S'): 2.40,
    ('Mg', 'F'): 1.90, ('Mg', 'Br'): 2.55, ('Mg', 'I'): 2.75,
    ('Ca', 'C'): 2.65, ('Ca', 'P'): 2.85, ('Ca', 'S'): 2.75,
    ('Ca', 'F'): 2.25, ('Ca', 'Br'): 2.90, ('Ca', 'I'): 3.10,
    ('Sr', 'C'): 2.80, ('Sr', 'P'): 3.00, ('Sr', 'S'): 2.90,
    ('Sr', 'F'): 2.40, ('Sr', 'Br'): 3.05, ('Sr', 'I'): 3.25,
    ('Ba', 'C'): 3.00, ('Ba', 'P'): 3.20, ('Ba', 'S'): 3.10,
    ('Ba', 'F'): 2.60, ('Ba', 'Br'): 3.25, ('Ba', 'I'): 3.45,
    # Y completions
    ('Y',  'P'): 2.75, ('Y',  'S'): 2.65, ('Y',  'F'): 2.05,
    ('Y',  'Br'): 2.75, ('Y',  'Se'): 2.80,
    # Nb / Tc extra (missing C, F, Br)
    ('Nb', 'C'): 2.20, ('Nb', 'F'): 1.95, ('Nb', 'Br'): 2.55,
    ('Mo', 'C'): 2.10,
    ('Tc', 'C'): 2.10, ('Tc', 'P'): 2.40, ('Tc', 'S'): 2.40,
    ('Tc', 'F'): 1.95, ('Tc', 'Br'): 2.55,
    # Ag extra
    ('Ag', 'C'): 2.10,
    # Ta / W / Re / Os extra
    ('Ta', 'C'): 2.15, ('Ta', 'F'): 1.90, ('Ta', 'Br'): 2.55,
    ('W',  'C'): 2.10,
    ('Re', 'C'): 2.10,
    ('Os', 'C'): 2.10,
    # Hg completions
    ('Hg', 'F'): 2.05,
    # Lanthanide completions (Pr / Nd / Pm / Sm / Eu / Tb / Dy / Ho / Er / Tm / Yb)
    ('Pr', 'N'): 2.50, ('Pr', 'O'): 2.40, ('Pr', 'Cl'): 2.70,
    ('Pr', 'I'): 3.05,
    ('Nd', 'N'): 2.48, ('Nd', 'O'): 2.38, ('Nd', 'Cl'): 2.68,
    ('Nd', 'I'): 3.02,
    ('Pm', 'N'): 2.46, ('Pm', 'O'): 2.36, ('Pm', 'Cl'): 2.66,
    ('Pm', 'I'): 3.00,
    ('Sm', 'N'): 2.45, ('Sm', 'O'): 2.35, ('Sm', 'Cl'): 2.65,
    ('Sm', 'I'): 2.98,
    ('Eu', 'N'): 2.44, ('Eu', 'O'): 2.34, ('Eu', 'Cl'): 2.64,
    ('Eu', 'I'): 2.96,
    ('Tb', 'N'): 2.42, ('Tb', 'O'): 2.32, ('Tb', 'Cl'): 2.62,
    ('Tb', 'I'): 2.94,
    ('Dy', 'N'): 2.40, ('Dy', 'O'): 2.30, ('Dy', 'Cl'): 2.60,
    ('Dy', 'I'): 2.92,
    ('Ho', 'N'): 2.38, ('Ho', 'O'): 2.28, ('Ho', 'Cl'): 2.58,
    ('Ho', 'I'): 2.90,
    ('Er', 'N'): 2.36, ('Er', 'O'): 2.26, ('Er', 'Cl'): 2.56,
    ('Er', 'I'): 2.88,
    ('Tm', 'N'): 2.34, ('Tm', 'O'): 2.24, ('Tm', 'Cl'): 2.54,
    ('Tm', 'I'): 2.86,
    ('Yb', 'N'): 2.32, ('Yb', 'O'): 2.22, ('Yb', 'Cl'): 2.52,
    ('Yb', 'I'): 2.88,
    # Main-group p-block metals (Al / Ga / In / Tl / Sn / Pb)
    ('Al', 'C'): 2.00, ('Al', 'N'): 2.00, ('Al', 'O'): 1.85,
    ('Al', 'P'): 2.35, ('Al', 'S'): 2.30, ('Al', 'F'): 1.75,
    ('Al', 'Cl'): 2.30, ('Al', 'Br'): 2.35,
    ('Ga', 'C'): 2.00, ('Ga', 'N'): 2.05, ('Ga', 'O'): 1.95,
    ('Ga', 'P'): 2.40, ('Ga', 'S'): 2.35, ('Ga', 'F'): 1.85,
    ('Ga', 'Cl'): 2.30, ('Ga', 'Br'): 2.40,
    ('In', 'C'): 2.20, ('In', 'N'): 2.25, ('In', 'O'): 2.15,
    ('In', 'P'): 2.55, ('In', 'S'): 2.55, ('In', 'F'): 2.00,
    ('In', 'Cl'): 2.45, ('In', 'Br'): 2.55,
    ('Tl', 'N'): 2.55, ('Tl', 'O'): 2.40, ('Tl', 'Cl'): 2.55,
    ('Tl', 'Br'): 2.70, ('Tl', 'I'): 2.85,
    ('Sn', 'C'): 2.15, ('Sn', 'N'): 2.25, ('Sn', 'O'): 2.10,
    ('Sn', 'P'): 2.55, ('Sn', 'S'): 2.45, ('Sn', 'F'): 1.95,
    ('Sn', 'Cl'): 2.40, ('Sn', 'Br'): 2.55,
    ('Pb', 'C'): 2.30, ('Pb', 'N'): 2.45, ('Pb', 'O'): 2.30,
    ('Pb', 'P'): 2.70, ('Pb', 'S'): 2.60, ('Pb', 'F'): 2.10,
    ('Pb', 'Cl'): 2.55, ('Pb', 'Br'): 2.70,
    # Bi
    ('Bi', 'C'): 2.20, ('Bi', 'N'): 2.40, ('Bi', 'O'): 2.25,
    ('Bi', 'P'): 2.65, ('Bi', 'S'): 2.55, ('Bi', 'F'): 2.00,
    ('Bi', 'Cl'): 2.50, ('Bi', 'Br'): 2.65, ('Bi', 'Se'): 2.65,
    # Selenium donor completions for remaining 3d/4d metals
    ('Sc', 'Se'): 2.65, ('Ti', 'Se'): 2.50, ('V',  'Se'): 2.45,
    ('Cr', 'Se'): 2.45, ('Mn', 'Se'): 2.45, ('Co', 'Se'): 2.35,
    ('Zn', 'Se'): 2.40,
    ('Zr', 'Se'): 2.65, ('Nb', 'Se'): 2.55, ('Mo', 'Se'): 2.50,
    ('Rh', 'Se'): 2.40, ('Cd', 'Se'): 2.65,
    ('Hf', 'Se'): 2.65, ('Ta', 'Se'): 2.55, ('W',  'Se'): 2.50,
    ('Re', 'Se'): 2.45, ('Os', 'Se'): 2.40, ('Ir', 'Se'): 2.40,
    ('Au', 'Se'): 2.45, ('Hg', 'Se'): 2.60,
    # Actinide completions (Ac, Pa, Np, Pu)
    ('Ac', 'N'): 2.60, ('Ac', 'O'): 2.40, ('Ac', 'Cl'): 2.75,
    ('Pa', 'N'): 2.50, ('Pa', 'O'): 2.30, ('Pa', 'Cl'): 2.65,
    ('Np', 'N'): 2.50, ('Np', 'O'): 2.25, ('Np', 'Cl'): 2.60,
    ('Pu', 'N'): 2.50, ('Pu', 'O'): 2.25, ('Pu', 'Cl'): 2.60,
    # -------------------------------------------------------------------
    # Welle-3 T3.1 (2026-05-15) — M-H (terminal hydride) entries.
    # Pure covalent-radius fallback gives 2.0-3.3 A which is >0.3-1.5 A
    # too long for crystallographically known terminal M-H (typically
    # 1.5-1.85 A).  Without these entries the topology check rejects
    # legitimate M-H sigma donors as "outside the bonding range" and the
    # downstream `_snap_md_distances_to_ideal` translates the H far from
    # the metal.  Universal — only element symbols.  Gated by env-flag
    # ``DELFIN_MH_TABLE_FALLBACK`` (default 0 → bit-exact HEAD).
    # Activation moved to module-tail (see ``_apply_mh_table_fallback``)
    # so the env-flag is read after init, not on import.
    # Values: 3d TM ~1.55-1.70, 4d TM ~1.60-1.75, 5d TM ~1.60-1.80
    # (relativistic contraction), main-group ~1.45-1.85, lanthanide
    # ~2.05-2.20.  Sources: CSD averages, Cordero 2008.
}


_METAL_HYDRIDE_BOND_LENGTHS: Dict[Tuple[str, str], float] = {
    # 3d TMs
    ('Sc','H'): 1.85, ('Ti','H'): 1.75, ('V','H'): 1.70,
    ('Cr','H'): 1.65, ('Mn','H'): 1.62, ('Fe','H'): 1.55,
    ('Co','H'): 1.53, ('Ni','H'): 1.50, ('Cu','H'): 1.55,
    ('Zn','H'): 1.60,
    # 4d TMs
    ('Y','H'): 1.95, ('Zr','H'): 1.85, ('Nb','H'): 1.75,
    ('Mo','H'): 1.70, ('Tc','H'): 1.68, ('Ru','H'): 1.60,
    ('Rh','H'): 1.58, ('Pd','H'): 1.55, ('Ag','H'): 1.65,
    ('Cd','H'): 1.70,
    # 5d TMs
    ('La','H'): 2.05, ('Hf','H'): 1.85, ('Ta','H'): 1.75,
    ('W','H'): 1.70, ('Re','H'): 1.68, ('Os','H'): 1.62,
    ('Ir','H'): 1.58, ('Pt','H'): 1.55, ('Au','H'): 1.62,
    ('Hg','H'): 1.70,
    # Lanthanides
    ('Ce','H'): 2.05, ('Pr','H'): 2.05, ('Nd','H'): 2.05,
    ('Pm','H'): 2.00, ('Sm','H'): 2.00, ('Eu','H'): 2.00,
    ('Gd','H'): 1.98, ('Tb','H'): 1.95, ('Dy','H'): 1.95,
    ('Ho','H'): 1.92, ('Er','H'): 1.92, ('Tm','H'): 1.90,
    ('Yb','H'): 1.90, ('Lu','H'): 1.88,
    # Actinides
    ('Ac','H'): 2.10, ('Th','H'): 2.00, ('Pa','H'): 1.95,
    ('U','H'): 1.95, ('Np','H'): 1.92, ('Pu','H'): 1.90,
    # Main-group metals
    ('Li','H'): 1.60, ('Na','H'): 1.85, ('K','H'): 2.20,
    ('Rb','H'): 2.30, ('Cs','H'): 2.45,
    ('Be','H'): 1.35, ('Mg','H'): 1.70, ('Ca','H'): 2.00,
    ('Sr','H'): 2.15, ('Ba','H'): 2.30,
    ('Al','H'): 1.55, ('Ga','H'): 1.55, ('In','H'): 1.70,
    ('Tl','H'): 1.85, ('Sn','H'): 1.70, ('Pb','H'): 1.80,
    ('Bi','H'): 1.85,
}


def _apply_mh_table_fallback() -> None:
    """Merge M-H entries into the master table.

    Default OFF for Welle-3 T3.1 (env-flag-gated).  Set
    ``DELFIN_MH_TABLE_FALLBACK=1`` to enable; default ``0`` keeps the
    legacy covalent-radius fallback (bit-exact pre-patch behaviour).
    Skips pairs already populated so explicit table entries always win.
    """
    import os as _os
    if int(_os.environ.get("DELFIN_MH_TABLE_FALLBACK", "0") or "0") <= 0:
        return
    for key, val in _METAL_HYDRIDE_BOND_LENGTHS.items():
        if key not in _METAL_LIGAND_BOND_LENGTHS:
            _METAL_LIGAND_BOND_LENGTHS[key] = val


_apply_mh_table_fallback()


# ---------------------------------------------------------------------------
# Welle-3 T3.3 (2026-05-15) — 84 new M-L pairs (As/Sb/Te/B/Se/La extras).
# Gated by DELFIN_NEW_ML_PAIRS (default 0).  Welle-5g Step-0 targeted revert
# per 5f-A bisect: a3edabe always-on inclusion attributed ~70% of -4064 NET
# (hapto -767 full-pool).  Default-OFF restores 9b1f541-equivalent behaviour.
# ---------------------------------------------------------------------------
_METAL_LIGAND_BOND_LENGTHS_T33: Dict[Tuple[str, str], float] = {
    # M-As donors (arsine, CSD 2.30-2.65 A)
    ('Sc','As'): 2.65, ('Ti','As'): 2.55, ('V','As'): 2.50,
    ('Cr','As'): 2.45, ('Mn','As'): 2.45, ('Fe','As'): 2.40,
    ('Co','As'): 2.35, ('Ni','As'): 2.30, ('Cu','As'): 2.35,
    ('Zn','As'): 2.45, ('Y','As'): 2.85, ('Zr','As'): 2.75,
    ('Nb','As'): 2.65, ('Mo','As'): 2.55, ('Tc','As'): 2.50,
    ('Ru','As'): 2.45, ('Rh','As'): 2.40, ('Pd','As'): 2.40,
    ('Ag','As'): 2.55, ('Cd','As'): 2.60, ('Hf','As'): 2.70,
    ('Ta','As'): 2.60, ('W','As'): 2.55, ('Re','As'): 2.50,
    ('Os','As'): 2.45, ('Ir','As'): 2.42, ('Pt','As'): 2.40,
    ('Au','As'): 2.45,
    # M-Sb donors (stibine, CSD 2.50-2.70 A)
    ('Cr','Sb'): 2.65, ('Mn','Sb'): 2.65, ('Co','Sb'): 2.55,
    ('Ni','Sb'): 2.50, ('Cu','Sb'): 2.55, ('Mo','Sb'): 2.65,
    ('Ru','Sb'): 2.60, ('Rh','Sb'): 2.55, ('Pd','Sb'): 2.55,
    ('Ag','Sb'): 2.65, ('W','Sb'): 2.65, ('Re','Sb'): 2.60,
    ('Os','Sb'): 2.55, ('Pt','Sb'): 2.55, ('Au','Sb'): 2.55,
    # M-Te donors (CSD 2.50-2.85 A)
    ('Sc','Te'): 2.85, ('V','Te'): 2.65, ('Cr','Te'): 2.70,
    ('Mn','Te'): 2.65, ('Co','Te'): 2.55, ('Ni','Te'): 2.50,
    ('Cu','Te'): 2.55, ('Zn','Te'): 2.65, ('Nb','Te'): 2.70,
    ('Mo','Te'): 2.65, ('Tc','Te'): 2.60, ('Ru','Te'): 2.60,
    ('Rh','Te'): 2.60, ('Pd','Te'): 2.65, ('Ag','Te'): 2.75,
    ('Cd','Te'): 2.75, ('W','Te'): 2.65, ('Re','Te'): 2.60,
    ('Ir','Te'): 2.55, ('Pt','Te'): 2.55, ('Au','Te'): 2.55,
    ('Hg','Te'): 2.70,
    # M-B donors (boryl, M-BR2/3; pi-backbond shortens vs covalent sum)
    ('Sc','B'): 2.40, ('Fe','B'): 2.00, ('Co','B'): 1.95,
    ('Ni','B'): 1.95, ('Cu','B'): 2.00, ('Ru','B'): 2.05,
    ('Rh','B'): 2.00, ('Pd','B'): 2.05, ('Cd','B'): 2.30,
    ('Os','B'): 2.05, ('Ir','B'): 2.00, ('Pt','B'): 2.05,
    ('Au','B'): 2.05,
    # M-Se/M-La extras
    ('Ag','Se'): 2.60, ('Tc','Se'): 2.50,
    ('La','Br'): 2.95, ('La','C'): 2.65, ('La','S'): 2.85, ('La','Se'): 2.95,
}


def _apply_new_ml_pairs_t33() -> None:
    """Merge Welle-3 T3.3 84 new M-L pairs into the master table.

    Default OFF (Welle-5g Step-0 targeted revert per 5f-A bisect).  Set
    ``DELFIN_NEW_ML_PAIRS=1`` to enable; default ``0`` restores the
    9b1f541-equivalent behaviour (bit-exact pre-a3edabe T3.3).  Skips
    pairs already populated so explicit table entries always win.
    """
    import os as _os
    if int(_os.environ.get("DELFIN_NEW_ML_PAIRS", "0") or "0") <= 0:
        return
    for key, val in _METAL_LIGAND_BOND_LENGTHS_T33.items():
        if key not in _METAL_LIGAND_BOND_LENGTHS:
            _METAL_LIGAND_BOND_LENGTHS[key] = val


_apply_new_ml_pairs_t33()


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


# Heavy metalloid sigma-donor elements (stibine Sb, arsine As, bismuthine Bi,
# telluro/seleno-ether Te/Se, heavy tetrel Ge/Sn/Pb).  Mirrors
# decompose._METALLOID_DONORS.  Used ONLY by the DELFIN_FFFREE_METALLOID_MD_LEN
# root fix below (default OFF -> byte-identical).
_METALLOID_MD_DONORS = frozenset({"Sb", "As", "Bi", "Te", "Se", "Ge", "Sn", "Pb"})


# Thread-local override for the DELFIN_FFFREE_METALLOID_MD_BOTH dual-distance pass:
# when its ``.value`` is True the metalloid M-D length is forced to the SHORT
# (covalent-sum) target even if the env flag is off, so smiles_to_xyz_isomers can
# build a SECOND, short-distance manifold and append it to the offset one.  Default
# unset -> byte-identical.  (Same threading.local pattern as _CONF_COMPLETE_ACTIVE.)
_METALLOID_MD_FORCE_SHORT = threading.local()


_ML_ME_MIN_N = 50          # minimum sample size per pair -- see _ml_me_band


_ML_ME_CACHE = {}          # path -> {(metal, donor): mu}


def _ml_me_band(metal_symbol: str, donor_symbol: str):
    """Terminal M=E MULTIPLE-BOND length, if a calibrated table is supplied.

    WHAT FOR (measured 16.08.2026).  `ml_len_realized` misses 192 of 953 systems, and the
    breakdown by element pair is extremely sharp:
        O-Re  8.86x | N-Re 7.88x | C-Ir 6.89x | O-V, O-W, O-Rh: ONLY failures
        Fe-N  0.25x | Mn-N 0.27x  (protected -- there is no multiple bond there)
    They are, throughout, terminal OXO and NITRIDO bonds on high-valent metals.  Cause:
    DELFIN writes them as a SINGLE BOND with a charged end atom (`[O-][Re]`), the bond
    order is NOT in the RDKit graph -- and `_get_ml_bond_length` looks up solely by
    (metal, donor element).  A terminal Re=O thereby gets exactly the length of a
    bridging Re-O.

    ⚠ WHY NO UNIFORM FACTOR.  From five textbook cases ~0.80 looked plausible.
    Measured on the crystals (76 pairs with n >= 50) the shortening is PAIR-SPECIFIC:
        U=O 0.743 (n=1100) | Re=O 0.813 | W=O 0.850 | V=O 0.890 | Mo=O 0.895 (n=4014)
        the MAJORITY 0.98-1.02 -- NO shortening
        Os=C 1.102, W=C 1.055 -- LENGTHENING
    Only 10 of 76 pairs shorten at all.  A blanket 0.80 factor produced 0.4-0.5 A errors at
    S/Se/C where today there is none.  Hence: the VALUE is looked up, not computed.

    ⚠ LICENSE.  The values are CCDC-derived and remain PRIVATE -- this file contains NO
    number from them.  The path comes via `DELFIN_ML_ME_BANDS`; if it is missing or the
    file is missing, this function returns None and the behaviour is unchanged.  The same
    pattern as for the torsion bands.
    """
    path = os.environ.get("DELFIN_ML_ME_BANDS", "")
    if not path:
        return None
    tbl = _ML_ME_CACHE.get(path)
    if tbl is None:
        tbl = {}
        try:
            with open(path) as fh:
                for line in fh:
                    if not line or line[0] == "#":
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 6 or f[0] == "metal":
                        continue
                    try:
                        # The threshold is enforced HERE once more, not only at
                        # generation time: a hand-edited table must not smuggle in a
                        # row with a tiny sample.
                        if int(f[5]) < _ML_ME_MIN_N:
                            continue
                        tbl[(f[0], f[1])] = float(f[2])
                    except Exception:
                        continue
        except Exception:
            tbl = {}
        _ML_ME_CACHE[path] = tbl
    return tbl.get((metal_symbol, donor_symbol))


def _ml_bond_kind(mol, metal_idx: int, donor_idx: int) -> str:
    """"me" for a TERMINAL multiple bond (oxo/nitrido/carbyne), otherwise "sigma".

    THE CRITERION IS STRUCTURAL, NOT FROM THE SMILES.  DELFIN writes terminal oxo as
    `[O-][Re]` -- a SINGLE BOND with a charged end atom --, so the multiple-bond order
    is not in the graph at all.  That is exactly why the eye, too, detects them
    structurally: a donor whose ONLY heavy neighbour is the metal can only be terminal, and
    a terminal O/N/C donor on a metal is chemically a multiple bond.

    Hydrogen does not count as a heavy neighbour; an OH or NH2 is thus correctly NOT an
    oxo.  Only the elements for which terminal M=E chemistry exists at all.

    Default OFF -> byte-identical (the switch is checked at the caller).
    """
    try:
        a = mol.GetAtomWithIdx(donor_idx)
        if a.GetSymbol() not in ("O", "N", "C"):
            return "sigma"
        if a.GetTotalNumHs() > 0:
            return "sigma"                       # OH / NH2 is not an oxo/imido
        heavy = [n.GetIdx() for n in a.GetNeighbors() if n.GetSymbol() != "H"]
        return "me" if heavy == [metal_idx] else "sigma"
    except Exception:
        return "sigma"


# M-C median from 40 000 CCDC crystals, per metal -- the TARGET VALUE, not estimated.
# Included only where |builder - crystal| >= 0.10 A AND n >= 30 measured bonds.
# Collected with `MANTA2/harness/mc_laenge_kristallmedian.py` (donor set of the EYE on the
# CRYSTAL = the coordination sphere; no second bond perception).
# Format: metal -> (median in A, n, previous value)  -- n and the old value stand next to
# it deliberately, so that every row can be checked without looking anything up.
_MC_LEN_CRYSTAL_SOURCE = {
    # carbonyl metals: all sat on the smooth placeholder 2.10 (Cr 2.05)
    "Mn": (1.809, 2137, 2.10), "Os": (1.908, 2834, 2.10), "Re": (1.925, 2453, 2.10),
    "Cr": (1.886, 2078, 2.05), "Tc": (1.918,   31, 2.10),
    # early transition metals: were too SHORT
    "Ti": (2.403,  875, 2.25), "Zr": (2.548,  751, 2.40), "Nb": (2.408,   97, 2.20),
    "Hf": (2.552,   99, 2.35), "Ta": (2.439,   94, 2.15), "Sc": (2.537,   55, 2.30),
    "V":  (2.281,  199, 2.15), "Y":  (2.670,  223, 2.55),
    # lanthanoids: were too LONG
    "La": (2.850,   33, 3.13), "Ce": (2.815,   59, 3.25), "Lu": (2.642,   86, 3.08),
    # isolated case
    "Ge": (1.964,  308, 2.46),
}


_MC_LEN_CRYSTAL_MEDIAN = {_m: _v[0] for _m, _v in _MC_LEN_CRYSTAL_SOURCE.items()}


def _get_ml_bond_length(metal_symbol: str, donor_symbol: str,
                        kind: str = "sigma") -> float:
    """Return estimated M-L bond length in Å.

    ``kind`` (16.08.2026): ``"sigma"`` = ordinary dative/covalent M-L bond (default,
    behaviour unchanged); ``"me"`` = TERMINAL multiple bond (oxo/nitrido/carbyne).  The
    structural criterion for "me" belongs to the CALLER, not here: a donor whose only
    heavy neighbour is the metal -- exactly how the eye detects it too, because the bond
    order is not in the graph.

    Lookup order:
    0. kind == "me" and a calibrated table supplied -> its value (see _ml_me_band)
    1. Specific (metal, donor) pair in _METAL_LIGAND_BOND_LENGTHS
    2. Sum of covalent radii + 0.5 Å (coordination bond correction)
    3. Default 2.0 Å
    """
    if kind == "me":
        _me = _ml_me_band(metal_symbol, donor_symbol)
        if _me is not None:
            return float(_me)
    # ===== M-C FROM THE CRYSTALS INSTEAD OF FROM PLACEHOLDERS (24.08.2026) ============
    # Measured on 40 000 CCDC crystals (14 566 with M-C), median per metal against the
    # value assumed here -- `harness/mc_laenge_kristallmedian.py`.  Over 32 metals the
    # median of the deviation is +0.033 A, so at its core the table is right.
    # It is wrong for THREE COHERENT FAMILIES, each in one direction:
    #   * the carbonyl metals all sit on the smooth placeholder 2.10 and are thereby
    #     0.16-0.29 A TOO LONG        (Mn n=2137, Os 2834, Re 2453, Cr 2078, Tc 31)
    #   * the early transition metals are 0.12-0.29 A TOO SHORT, where M-C really is
    #     long                         (Ti 875, Zr 751, Nb 97, Hf 99, Ta 94, Sc, V, Y)
    #   * the lanthanoids are 0.28-0.44 A TOO LONG              (La 33, Ce 59, Lu 86)
    #   * Ge is the isolated case with +0.496
    # Correctly calibrated are, of all things, the most frequent: Fe (n=11 428) -0.024,
    # Pt -0.000, Cu -0.001, Rh -0.003, Ir +0.007, Pd/Au -0.014.
    #
    # WHY THIS MATTERS: the eye's donor model accepts C as a donor only up to 1.25 x the
    # ideal bond.  An Os-C that is 0.19 A too long slips above that, no longer counts as
    # a donor, the metal CN drops -- and because the CN gate is EXACT, NO frame passes
    # the topology check any more.  Measured on `ACAHUH`: crystal CN 4, build CN 0 in
    # all twelve frames, `topo_correct_frame = false`, although the frames are clean.
    #
    # ⚠ ONLY WHERE THE EVIDENCE HOLDS: included are metals with |delta| >= 0.10 A AND
    # n >= 30 measured bonds.  The value is the MEDIAN, not the covalent sum and not a
    # blanket-deleted offset -- a blanket intervention would be the same error as the
    # one it fixes.
    # Default OFF -> byte-identical.  ONE read site.
    # ── THE TABLE CAN BE RESTRICTED TO A LIST OF METALS ──────────────────────────
    # OCCASION (27.08.2026, evening).  `mclen6k` has been recovered and answers its
    # pre-registration SPLIT: `topo_correct_frame` NET **+40** (46 won against 6
    # lost) -- the chain M-C length -> donor cut -> CN -> topology gate is thereby
    # PROVEN, H0 refuted.  But the run does not land: 2 hard `ccdc_isomer_lost`
    # (ECOZUV, WIKBOJ, fingerprint IDENTICAL in both arms, hence REAL losses),
    # 5 polyhedra, 7 isomers_lost.
    #
    # 🔑 AND THE LOSSES ARE NOT EVERYWHERE.  Counted per metal:
    #        WINNERS 46    Re 22 · Mn 14 · Os 6 · Tc 3 · Ti 1
    #        LOSERS  6     Os  3 · Cr  1 · Mn 1 · Zr 1
    #    Re stands at 22 : 0, Mn 14 : 1 -- Os on the other hand 6 : 3, and Os has the
    #    STRONGEST crystal evidence (n=2834).  Prior knowledge and A/B contradict there.
    #
    # ⚠️ THIS SELECTION IS POST HOC.  It comes from the result of the very run it is
    #    supposed to improve -- that is selection on the test set and would be
    #    worthless as a finding.  The switch therefore exists ONLY so that a NEW,
    #    pre-registered run can confirm or refute it.
    #    Until that one is through, "Re+Mn is better" is a HYPOTHESIS.
    #
    # Default: empty -> ALL metals of the table -> byte-identical to before.
    if (donor_symbol in ("C", "Si")
            and os.environ.get("DELFIN_FFFREE_MC_LEN_CRYSTAL", "0") == "1"):
        _mc = _MC_LEN_CRYSTAL_MEDIAN.get(metal_symbol)
        _only = os.environ.get("DELFIN_FFFREE_MC_LEN_METALS", "").strip()
        if _mc is not None and _only:
            _erlaubt = {_s.strip() for _s in _only.split(",") if _s.strip()}
            if metal_symbol not in _erlaubt:
                _mc = None                 # not in the list -> old path
        if _mc is not None:
            return float(_mc)
    key = (metal_symbol, donor_symbol)
    if key in _METAL_LIGAND_BOND_LENGTHS:
        return _METAL_LIGAND_BOND_LENGTHS[key]
    r_m = _COVALENT_RADII.get(metal_symbol)
    r_d = _COVALENT_RADII.get(donor_symbol)
    if r_m is not None and r_d is not None:
        # Heavy metalloid sigma-donors (Sb/As/Bi/Te/Se/Ge/Sn/Pb): the generic
        # +_METAL_ROW_OFFSET coordination-bond correction is calibrated for HARD
        # N/O/P/S donors, where the covalent-radius sum UNDER-estimates the M-L
        # bond.  For these SOFT, already-large-radius donors it OVER-shoots: the
        # seated AND UFF-pinned M-D bond lands ~0.4 Å too long (e.g. Cu-Sb
        # 1.32+1.39+0.40 = 3.11 vs real ~2.6), so the donor reads as DETACHED
        # (eye: mdbreak "M-D too_long" + smiles_topology on ~70% of heavy_donor
        # frames; silent on the CCDC crystals -> a real construction defect).
        # When DELFIN_FFFREE_METALLOID_MD_LEN=1, drop the row offset for metalloid
        # donors so the target is the plain covalent sum (~2.7 Å == polyhedra.
        # md_distance), the physical M-metalloid distance.  Env-gated, default
        # byte-identical; a construction-chain ROOT fix -- it sets the seat / UFF-pin
        # / snap / tolerance target BEFORE emit, so the raw geometry is built
        # correct (NOT a post-hoc geometry repair).
        if donor_symbol in _METALLOID_MD_DONORS and (
                getattr(_METALLOID_MD_FORCE_SHORT, "value", False)
                or os.environ.get("DELFIN_FFFREE_METALLOID_MD_LEN", "0") == "1"):
            return r_m + r_d
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
# CN=5 chemistry-aware classifier (universal: element + formal charge ->
# d-count + ligand-field signals).  Used when DELFIN_CN5_GEOM_AWARE=1.
# Default OFF for bit-exact HEAD behaviour.
# ---------------------------------------------------------------------------

# Periodic-table group number per transition-metal centre.  d_count = group -
# oxidation_state, clamped to [0, 10].  Non-d-block returns None.
_METAL_GROUP_NUMBER: Dict[str, int] = {
    'Sc': 3, 'Y': 3, 'La': 3, 'Ac': 3,
    'Ti': 4, 'Zr': 4, 'Hf': 4,
    'V': 5,  'Nb': 5, 'Ta': 5,
    'Cr': 6, 'Mo': 6, 'W': 6,
    'Mn': 7, 'Tc': 7, 'Re': 7,
    'Fe': 8, 'Ru': 8, 'Os': 8,
    'Co': 9, 'Rh': 9, 'Ir': 9,
    'Ni': 10, 'Pd': 10, 'Pt': 10,
    'Cu': 11, 'Ag': 11, 'Au': 11,
    'Zn': 12, 'Cd': 12, 'Hg': 12,
}


# pi-acceptor donor elements (phosphines/arsines/stibines).  Carbonyl C and
# NHC C are detected graph-based in the rich helper, not via label.
_PI_ACCEPTOR_DONOR_ELEMS: Tuple[str, ...] = ('P', 'As', 'Sb')


def _cn5_d_electron_count(
    metal_symbol: str,
    formal_charge: int,
) -> Optional[int]:
    """Return d-electron count for ``metal_symbol`` at ox-state ``formal_charge``.

    Returns ``None`` for non-d-block (main-group, lanthanide/actinide, unknown).
    Uses ``d = group - oxidation_state`` clamped to ``[0, 10]``.
    """
    grp = _METAL_GROUP_NUMBER.get(metal_symbol)
    if grp is None:
        return None
    try:
        ox = int(formal_charge)
    except Exception:
        return None
    d = grp - ox
    if d < 0:
        return 0
    if d > 10:
        return 10
    return d


def _cn5_count_pi_acceptor_donors(donor_labels: List[str]) -> int:
    """Count pi-acceptor donors among ``donor_labels`` (P/As/Sb only)."""
    n = 0
    for lbl in donor_labels:
        sym = "".join(ch for ch in lbl if ch.isalpha())
        if sym in _PI_ACCEPTOR_DONOR_ELEMS:
            n += 1
    return n


def _cn5_has_tridentate_chelate(
    chelate_pairs: List[FrozenSet],
    n_coord: int,
) -> bool:
    """Detect chelating ligand with denticity >= 3 via union-find on
    ``chelate_pairs`` (donor-list-index pairs).  Tridentate-+ chelates
    force three donors into one plane -> SP preference."""
    if n_coord != 5 or not chelate_pairs:
        return False
    parent = list(range(n_coord))

    def _find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def _union(a: int, b: int) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[ra] = rb

    for cp in chelate_pairs:
        items = list(cp)
        if len(items) != 2:
            continue
        a, b = int(items[0]), int(items[1])
        if 0 <= a < n_coord and 0 <= b < n_coord:
            _union(a, b)

    from collections import Counter
    sizes = Counter(_find(i) for i in range(n_coord))
    return any(s >= 3 for s in sizes.values())


def _classify_cn5_geometry_from_labels(
    metal_symbol: str,
    formal_charge: int,
    donor_labels: List[str],
    chelate_pairs: Optional[List[FrozenSet]] = None,
) -> str:
    """Universal CN=5 polyhedron classifier (label-only facade).

    Returns ``'SP'`` (C4v square-pyramidal) or ``'TBP'`` (D3h
    trigonal-bipyramidal).  Decision tree:

    1. d-electron count = group - formal_charge:
       - d8  -> SP   (Ni(II), Pd(II), Pt(II), Co(I), Rh(I), Ir(I), Au(III))
       - d4/d6 LS -> SP   (Cr(II), Mn(III), Fe(II) LS, Co(III), Ru(II),
                          Rh(III), Ir(III))
       - d0/d1/d2/d10 -> TBP
       - d3/d5/d7/d9  -> TBP (default; weak preference)
    2. Override #1 - >= 2 pi-acceptor donors (P/As/Sb) -> SP.
    3. Override #2 - tridentate (or larger) chelate present -> SP.
    4. Safe fallback for non-d-block / unknown elements: ``'TBP'``.

    Pure-scalar API; no RDKit needed.  Rich graph-based variant:
    :func:`_classify_cn5_geometry`.
    """
    d = _cn5_d_electron_count(metal_symbol, formal_charge)
    if d is None:
        pref = 'TBP'
    elif d == 8:
        pref = 'SP'
    elif d in (4, 6):
        pref = 'SP'
    elif d in (0, 1, 2, 10):
        pref = 'TBP'
    else:
        pref = 'TBP'

    if _cn5_count_pi_acceptor_donors(donor_labels) >= 2 and pref == 'TBP':
        pref = 'SP'

    chelate_list = list(chelate_pairs) if chelate_pairs else []
    if _cn5_has_tridentate_chelate(chelate_list, len(donor_labels)) and pref == 'TBP':
        pref = 'SP'

    return pref


def _classify_cn5_geometry(
    mol,
    metal_idx: int,
    donor_idxs: List[int],
) -> str:
    """Universal CN=5 polyhedron classifier (rich graph-based version).

    Same return contract as :func:`_classify_cn5_geometry_from_labels`
    (``'SP'`` or ``'TBP'``) but uses the RDKit molecule directly so it
    can also detect:

      * Carbonyl donors (C double-bonded to terminal O) -> pi-acceptor.
      * NHC donors (sp2 aromatic C with >= 2 N neighbours) -> pi-acceptor.
      * Sterically bulky donors (>= 3 heavy non-metal neighbours).

    Args:
        mol: RDKit ``Mol`` containing the metal centre and its donors.
        metal_idx: 0-based atom index of the metal centre.
        donor_idxs: list of 0-based atom indices for the five donor atoms.

    Returns:
        ``'SP'`` or ``'TBP'``.  ``'TBP'`` on any parse error (safe fallback).
    """
    try:
        if mol is None or metal_idx is None or donor_idxs is None:
            return 'TBP'
        if len(donor_idxs) != 5:
            return 'TBP'

        metal_atom = mol.GetAtomWithIdx(int(metal_idx))
        metal_symbol = metal_atom.GetSymbol()
        formal_charge = int(metal_atom.GetFormalCharge() or 0)

        donor_labels: List[str] = []
        n_pi_rich = 0
        n_bulky = 0
        for li, di in enumerate(donor_idxs):
            try:
                datom = mol.GetAtomWithIdx(int(di))
            except Exception:
                donor_labels.append(f"X{li}")
                continue
            sym = datom.GetSymbol()
            donor_labels.append(f"{sym}{li}")
            if sym == 'C':
                is_carbonyl = False
                for nb in datom.GetNeighbors():
                    if nb.GetIdx() == int(metal_idx):
                        continue
                    if nb.GetSymbol() != 'O':
                        continue
                    bnd = mol.GetBondBetweenAtoms(datom.GetIdx(), nb.GetIdx())
                    if bnd is not None and bnd.GetBondTypeAsDouble() >= 1.5:
                        heavy_nbs = [
                            x for x in nb.GetNeighbors()
                            if x.GetSymbol() != 'H'
                        ]
                        if len(heavy_nbs) == 1:
                            is_carbonyl = True
                            break
                if is_carbonyl:
                    n_pi_rich += 1
                else:
                    n_nbs = [
                        x for x in datom.GetNeighbors()
                        if x.GetSymbol() == 'N' and x.GetIdx() != int(metal_idx)
                    ]
                    if len(n_nbs) >= 2 and datom.GetIsAromatic():
                        n_pi_rich += 1
            elif sym in _PI_ACCEPTOR_DONOR_ELEMS:
                n_pi_rich += 1
            heavy_nbs_d = [
                x for x in datom.GetNeighbors()
                if x.GetSymbol() != 'H' and x.GetIdx() != int(metal_idx)
            ]
            if len(heavy_nbs_d) >= 3:
                n_bulky += 1

        chelate_pairs: List[FrozenSet] = []
        try:
            ring_info = mol.GetRingInfo()
            for i in range(len(donor_idxs)):
                for j in range(i + 1, len(donor_idxs)):
                    ai = int(donor_idxs[i])
                    aj = int(donor_idxs[j])
                    for ring in ring_info.AtomRings():
                        if (
                            ai in ring
                            and aj in ring
                            and int(metal_idx) in ring
                            and len(ring) <= 7
                        ):
                            chelate_pairs.append(frozenset([i, j]))
                            break
        except Exception:
            pass

        pref = _classify_cn5_geometry_from_labels(
            metal_symbol=metal_symbol,
            formal_charge=formal_charge,
            donor_labels=donor_labels,
            chelate_pairs=chelate_pairs,
        )

        if n_pi_rich >= 2 and pref == 'TBP':
            pref = 'SP'

        if n_bulky >= 4 and pref == 'TBP':
            pref = 'SP'

        return pref
    except Exception:
        return 'TBP'


# ---------------------------------------------------------------------------
# Preferred CN=6 geometry per metal
# ---------------------------------------------------------------------------
_PREFERRED_CN6_GEOMETRY: Dict[str, str] = {
    # Almost all TMs: octahedral strongly preferred.
    # Trigonal-prismatic only for d0 with specific dithiolene ligands.
    # Default is OH for everything not listed here.
    'Mo': 'OH', 'W': 'OH', 'Re': 'OH',
}


# ---------------------------------------------------------------------------
# Iter-2 Subagent A - complete polyhedra per CN (Polya enumeration).
# Used when DELFIN_ALL_POLYHEDRA=1 (default ON) for call sites that
# previously picked ONE preferred polyhedron per metal.  Stable-first order
# keeps legacy preferred polyhedron leading; new entries appended LAST.
# Bit-exact HEAD when env-flag is 0.
# ---------------------------------------------------------------------------
_ALL_CN3_POLYHEDRA: Tuple[str, ...] = ('TP', 'TS')          # trig-planar, T-shape


_ALL_CN4_POLYHEDRA: Tuple[str, ...] = ('TH', 'SQ', 'SS')    # tet, sq-planar, see-saw


_ALL_CN5_POLYHEDRA: Tuple[str, ...] = ('TBP', 'SP')         # trig-bipyr, sq-pyr


_ALL_CN6_POLYHEDRA: Tuple[str, ...] = ('OH', 'TPR')         # octahedral, trig-prism


_ALL_CN7_POLYHEDRA: Tuple[str, ...] = ('PBP', 'COH')        # pent-bipyr, capped-oct


_ALL_CN8_POLYHEDRA: Tuple[str, ...] = ('SAP', 'DD')         # sq-antiprism, dodecahedral


_ALL_POLYHEDRA_BY_CN: Dict[int, Tuple[str, ...]] = {
    3: _ALL_CN3_POLYHEDRA,
    4: _ALL_CN4_POLYHEDRA,
    5: _ALL_CN5_POLYHEDRA,
    6: _ALL_CN6_POLYHEDRA,
    7: _ALL_CN7_POLYHEDRA,
    8: _ALL_CN8_POLYHEDRA,
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
