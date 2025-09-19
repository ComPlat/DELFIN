import re
import os
import logging
from typing import List, Tuple, Optional

# ------------------------------------------------------------------------------------
# Transition metals (IUPAC blocks)
# ------------------------------------------------------------------------------------
_TM_LIST: List[str] = [
    'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
    'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
    'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Rf',
    'Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn'
]

_TM_D3  = {'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn'}
_TM_D45 = {'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
           'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',
           'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn'}

# Regex: extract an element symbol from a token like "Fe", "Fe(1)", "Fe1"
_ELEM_FROM_TOKEN = re.compile(r'^([A-Za-z]{1,2})(?![A-Za-z])')


# ------------------------------------------------------------------------------------
# Metal detection
# ------------------------------------------------------------------------------------
def search_transition_metals(inputfile: str) -> List[str]:
    """
    Return a list of unique transition-metal symbols present in a coordinate-like file.
    Tolerant to formats (plain XYZ, ORCA '* xyz' block, raw lines).
    """
    try:
        with open(inputfile, 'r', errors="ignore") as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"File '{inputfile}' not found.")
        return []

    found: List[str] = []
    seen = set()
    for line in lines:
        ls = line.strip()
        if not ls or ls.startswith('*'):
            continue
        parts = ls.split()
        if not parts:
            continue
        m = _ELEM_FROM_TOKEN.match(parts[0])
        if not m:
            continue
        sym = m.group(1)
        # Normalize case, e.g., 'fe' -> 'Fe'
        sym = sym[0].upper() + (sym[1].lower() if len(sym) > 1 else "")
        if sym in _TM_LIST and sym not in seen:
            seen.add(sym)
            found.append(sym)
    return found


def classify_tm_presence(found_metals: List[str]) -> str:
    """
    Classify the TM set:
      - 'none'  : no TMs
      - '3d'    : only 3d TMs present (subset of _TM_D3)
      - '4d5d'  : only 4d/5d TMs present (subset of _TM_D45)
      - 'mixed' : both 3d and 4d/5d present
    """
    if not found_metals:
        return 'none'
    has_3d  = any(m in _TM_D3  for m in found_metals)
    has_d45 = any(m in _TM_D45 for m in found_metals)
    if has_3d and has_d45:
        return 'mixed'
    if has_d45:
        return '4d5d'
    return '3d'


# ------------------------------------------------------------------------------------
# Relativity policy and basis selection
# d3 metals → non-rel; any 4d/5d (or mixed) → relativistic
# ------------------------------------------------------------------------------------
def _rel_method_token(config: dict) -> str:
    """
    Map CONTROL 'relativity' to the ORCA method token (printed on the '!' line).
    Supports: none|zora|x2c|dkh|dkh2. Returns '' if none/unknown.
    """
    rel = str(config.get("relativity", "none")).strip().lower()
    return {"zora": "ZORA", "x2c": "X2C", "dkh": "DKH", "dkh2": "DKH"}.get(rel, "")


def _should_use_rel(found_metals: List[str]) -> bool:
    """
    Relativity if classification is '4d5d' or 'mixed'.
    (Uses _TM_D3 actively via classify_tm_presence to satisfy policy.)
    """
    cls = classify_tm_presence(found_metals)
    return cls in ('4d5d', 'mixed')


def select_rel_and_aux(found_metals: List[str], config: dict) -> Tuple[str, str, bool]:
    """
    Decide relativity token and aux-JK set following the policy:
      - only 3d → non-rel ⇒ rel_token='', aux=aux_jk
      - any 4d/5d (or mixed) → rel ⇒ rel_token from CONTROL, aux=aux_jk_rel
    Returns (rel_token, aux_jk_value, use_rel_bool).
    """
    if _should_use_rel(found_metals):
        return _rel_method_token(config), str(config.get("aux_jk_rel", "")).strip(), True
    return "", str(config.get("aux_jk", "")).strip(), False


def set_main_basisset(found_metals: List[str], config: dict) -> Tuple[str, Optional[str]]:
    """
    Choose orbital bases for the '!' line and the metal override (inline NewGTO):
      - any 4d/5d (or mixed) → use main_basisset_rel / metal_basisset_rel
      - only 3d or none      → use main_basisset / metal_basisset
    There is no 'basisset_org' fallback anymore.
    Returns (main_basis, metal_basis or None). If no TM is present → metal_basis=None.
    """
    use_rel = _should_use_rel(found_metals)
    if use_rel:
        main  = str(config.get("main_basisset_rel", "")).strip() or "ZORA-def2-SVP"
        metal = str(config.get("metal_basisset_rel", "")).strip() or None
    else:
        main  = str(config.get("main_basisset", "")).strip() or "def2-SVP"
        metal = str(config.get("metal_basisset", "")).strip() or None

    if not found_metals:
        metal = None  # no TM → do not attach any per-atom overrides

    # Explicit debug helps to see why Fe (3d) stays non-rel:
    cls = classify_tm_presence(found_metals)
    logging.debug(f"[utils] TM class={cls} → use_rel={use_rel}; main='{main}'; metal='{metal}'")
    return main, metal


# ------------------------------------------------------------------------------------
# Electron counting (from CONTROL + geometry file)
# ------------------------------------------------------------------------------------
_ATOM_ELECTRONS = {
    "H": 1, "He": 2,
    "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
    "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98,
    "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
    "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113,
    "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
}


def _parse_control_for_input_file(control_file_path: str) -> str:
    """
    Find 'input_file=...' in CONTROL; fallback to 'input.txt' in the same folder.
    """
    control_dir = os.path.dirname(control_file_path)
    candidate = os.path.join(control_dir, "input.txt")
    try:
        with open(control_file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.strip().lower().startswith("input_file"):
                    eq = line.split("=", 1)
                    if len(eq) == 2:
                        val = eq[1].strip()
                        if val:
                            return os.path.join(control_dir, val)
    except Exception as e:
        logging.info(f"CONTROL parse fallback (input_file): {e}")
    return candidate


def _looks_like_xyz(lines: List[str]) -> bool:
    """Return True if first line is an int (Natoms) and there are >= 2 lines."""
    if len(lines) < 2:
        return False
    try:
        int(str(lines[0]).strip().split()[0])
        return True
    except Exception:
        return False


def calculate_total_electrons_txt(control_file_path: str) -> Optional[Tuple[int, int]]:
    """
    Read the geometry file referenced by CONTROL ('input_file=...') if present,
    otherwise 'input.txt' next to CONTROL. Count total electrons by summing
    atomic numbers of element tokens in the coordinate block.
    Returns (total_electrons, guess_multiplicity_1_or_2), or None on error.
    """
    input_file_path = _parse_control_for_input_file(control_file_path)
    if not os.path.exists(input_file_path):
        logging.error(f"Input file '{input_file_path}' not found.")
        return None

    try:
        with open(input_file_path, 'r', encoding='utf-8', errors="ignore") as input_file:
            lines = input_file.readlines()
    except Exception as e:
        logging.error(f"Error reading the file '{input_file_path}': {e}")
        return None

    # If XYZ: skip first two lines (natoms + comment)
    coord_lines = lines[2:] if _looks_like_xyz(lines) else lines[:]

    total_electrons_txt = 0
    for line in coord_lines:
        ls = line.strip()
        if not ls or ls == '*':
            continue
        parts = ls.split()
        if not parts:
            continue
        m = _ELEM_FROM_TOKEN.match(parts[0])
        if not m:
            continue
        sym = m.group(1)
        sym = sym[0].upper() + (sym[1].lower() if len(sym) > 1 else "")
        if sym in _ATOM_ELECTRONS:
            total_electrons_txt += _ATOM_ELECTRONS[sym]
        else:
            logging.warning(f"Unknown element '{sym}' in line: {line.strip()}")

    multiplicity_guess = 1 if (total_electrons_txt % 2 == 0) else 2
    return total_electrons_txt, multiplicity_guess
