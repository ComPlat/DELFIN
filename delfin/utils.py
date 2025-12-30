import re
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path

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

logger = get_logger(__name__)


# ------------------------------------------------------------------------------------
# String utilities
# ------------------------------------------------------------------------------------

def normalize_str(value: Any) -> str:
    """Normalize a value to lowercase trimmed string.

    Args:
        value: Any value to normalize

    Returns:
        Lowercase trimmed string representation
    """
    return str(value).strip().lower()


# ------------------------------------------------------------------------------------
# Git version tracking
# ------------------------------------------------------------------------------------
def get_git_commit_info() -> Optional[str]:
    """Get current git commit hash and status for reproducibility tracking.

    Returns:
        String like "0558f4b" or "0558f4b-dirty" if uncommitted changes,
        or None if not in a git repository.
    """
    try:
        # Get the directory where DELFIN is installed
        delfin_dir = Path(__file__).parent.parent

        # Get commit hash (short)
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=delfin_dir,
            capture_output=True,
            text=True,
            timeout=2
        )

        if result.returncode != 0:
            return None

        commit_hash = result.stdout.strip()

        # Check if there are uncommitted changes
        result = subprocess.run(
            ["git", "diff-index", "--quiet", "HEAD", "--"],
            cwd=delfin_dir,
            timeout=2
        )

        if result.returncode != 0:
            commit_hash += "-dirty"

        return commit_hash

    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        return None


# ------------------------------------------------------------------------------------
# Metal detection
# ------------------------------------------------------------------------------------
def search_transition_metals(inputfile: str) -> List[str]:
    """Detect transition metals in molecular geometry file.

    Scans coordinate file for transition metal elements (Sc-Cn)
    and returns unique symbols found. Tolerant to various formats
    including XYZ files and ORCA input blocks.

    Args:
        inputfile: Path to coordinate file (XYZ, ORCA input, etc.)

    Returns:
        List of unique transition metal symbols (e.g., ['Fe', 'Co'])
    """
    try:
        with open(inputfile, 'r', errors="ignore") as f:
            lines = f.readlines()
    except FileNotFoundError:
        logger.error(f"File '{inputfile}' not found.")
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
    """Classify transition metal composition for computational policy.

    Determines computational approach based on d-orbital period:
    - 3d metals: Sc-Zn (typically no relativistic effects needed)
    - 4d/5d metals: Y-Cd, Hf-Hg (require relativistic corrections)

    Args:
        found_metals: List of transition metal symbols

    Returns:
        Classification: 'none', '3d', '4d5d', or 'mixed'
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
def _rel_method_token(config: Dict[str, Any]) -> str:
    """Convert relativity setting to ORCA method keyword.

    Maps configuration relativity option to appropriate ORCA token
    for the method line (!-line).

    Args:
        config: Configuration dictionary with 'relativity' key

    Returns:
        ORCA relativity token ('ZORA', 'X2C', 'DKH', or '')
    """
    rel = str(config.get("relativity", "none")).strip().lower()
    return {"zora": "ZORA", "x2c": "X2C", "dkh": "DKH", "dkh2": "DKH"}.get(rel, "")


def _should_use_rel(found_metals: List[str]) -> bool:
    """Determine if relativistic corrections are needed.

    Policy: Apply relativistic methods for 4d/5d metals or mixed systems.
    Pure 3d metal systems use non-relativistic methods.

    Args:
        found_metals: List of transition metal symbols

    Returns:
        True if relativistic corrections should be applied
    """
    cls = classify_tm_presence(found_metals)
    return cls in ('4d5d', 'mixed')


def select_rel_and_aux(found_metals: List[str], config: Dict[str, Any]) -> Tuple[str, str, bool]:
    """Select relativity method and auxiliary basis sets based on metal composition.

    Implements intelligent basis set policy:
    - 3d metals only: non-relativistic with standard auxiliary sets
    - 4d/5d metals: relativistic methods with appropriate auxiliary sets

    Args:
        found_metals: List of transition metal symbols found
        config: Configuration dictionary with basis set settings

    Returns:
        Tuple of (relativity_token, auxiliary_basis, use_relativistic_flag)
    """
    if _should_use_rel(found_metals):
        return _rel_method_token(config), str(config.get("aux_jk_rel", "")).strip(), True
    return "", str(config.get("aux_jk", "")).strip(), False


def set_main_basisset(found_metals: List[str], config: Dict[str, Any]) -> Tuple[str, Optional[str]]:
    """Select appropriate basis sets based on transition metal composition.

    Implements basis set selection policy:
    - 3d metals: standard basis sets (def2-SVP family)
    - 4d/5d metals: relativistic basis sets (ZORA-def2-SVP family)
    - No metals: main basis only, no per-atom overrides

    Args:
        found_metals: List of transition metal symbols detected
        config: Configuration with basis set preferences

    Returns:
        Tuple of (main_basisset, metal_basisset_or_None)
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
    logger.debug(f"[utils] TM class={cls} → use_rel={use_rel}; main='{main}'; metal='{metal}'")
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


def _parse_control_for_input_file(control_file_path: str) -> Path:
    """Find 'input_file=...' in CONTROL; fallback to 'input.txt' in the same folder."""
    control_path = resolve_path(control_file_path)
    control_dir = control_path.parent
    candidate = control_dir / "input.txt"
    try:
        with control_path.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.strip().lower().startswith("input_file"):
                    eq = line.split("=", 1)
                    if len(eq) == 2:
                        val = eq[1].strip()
                        if val:
                            return resolve_path(control_dir / val)
    except Exception as e:
        logger.info(f"CONTROL parse fallback (input_file): {e}")
    return resolve_path(candidate)


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
    """Calculate total electron count and suggest spin multiplicity.

    Reads molecular geometry from CONTROL-referenced file and sums
    atomic numbers to determine total electron count. Reads charge from
    CONTROL.txt and calculates actual electron count after charge correction.
    Suggests initial spin multiplicity based on even/odd electron number
    AFTER charge correction.

    Args:
        control_file_path: Path to CONTROL.txt file

    Returns:
        Tuple of (total_electrons_neutral, suggested_multiplicity) or None on error
        Note: total_electrons_neutral is the electron count BEFORE charge correction
    """
    # Read charge from CONTROL.txt
    control_path = resolve_path(control_file_path)
    charge = 0
    try:
        with control_path.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.strip().lower().startswith("charge"):
                    eq = line.split("=", 1)
                    if len(eq) == 2:
                        charge_str = eq[1].strip()
                        try:
                            charge = int(charge_str)
                        except ValueError:
                            logger.warning(f"Could not parse charge '{charge_str}', using 0")
                            charge = 0
                        break
    except Exception as e:
        logger.warning(f"Could not read charge from CONTROL.txt: {e}, using 0")

    input_file_path = _parse_control_for_input_file(control_file_path)
    if not input_file_path.exists():
        logger.error(f"Input file '{input_file_path}' not found.")
        return None

    try:
        with input_file_path.open('r', encoding='utf-8', errors="ignore") as input_file:
            file_content = input_file.read()
            lines = file_content.splitlines(keepends=True)
    except Exception as e:
        logger.error(f"Error reading the file '{input_file_path}': {e}")
        return None

    # Check if input is SMILES (simple heuristic: single line = SMILES)
    non_empty_lines = [l for l in lines if l.strip()]
    if len(non_empty_lines) == 1:
        try:
            from delfin.smiles_converter import smiles_to_xyz
            logger.info(f"Single line detected in {input_file_path.name}, treating as SMILES")
            xyz_content, error = smiles_to_xyz(file_content.strip())
            if error:
                logger.warning(f"SMILES conversion failed: {error}")
            elif xyz_content:
                lines = xyz_content.splitlines(keepends=True)
                logger.debug("Successfully converted SMILES to XYZ for electron counting")
        except ImportError:
            logger.warning("SMILES converter not available, cannot convert single-line input")
        except Exception as e:
            logger.warning(f"SMILES conversion failed: {e}")

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
            logger.warning(f"Unknown element '{sym}' in line: {line.strip()}")

    # Calculate actual electron count after charge correction
    total_electrons_actual = total_electrons_txt - charge

    # Determine multiplicity based on ACTUAL electron count (after charge correction)
    multiplicity_guess = 1 if (total_electrons_actual % 2 == 0) else 2

    logger.debug(f"Electron count: neutral={total_electrons_txt}, charge={charge}, "
                 f"actual={total_electrons_actual}, multiplicity={multiplicity_guess}")

    return total_electrons_txt, multiplicity_guess
