# OCCUPIER.py

import math
import os
import re
import threading
import time
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from delfin.common.logging import get_logger, add_file_handler
from delfin.common.paths import resolve_path
from delfin.global_manager import get_global_manager

from .config import OCCUPIER_parser, read_control_file
from .occupier_sequences import (
    infer_species_delta,
    resolve_sequences_for_delta,
    append_sequence_overrides,
    remove_existing_sequence_blocks,
)
from .utils import (
    set_main_basisset,
    search_transition_metals,
    select_rel_and_aux,
)
from .reporting import generate_summary_report_OCCUPIER
from .orca import run_orca
from .xyz_io import (
    build_qmmm_block,
    split_qmmm_sections,
    _ensure_qmmm_implicit_model,
    normalize_xyz_body,
)
from .parallel_classic_manually import (
    _WorkflowManager,
    WorkflowJob,
    _update_pal_block,
    determine_effective_slots,
    normalize_parallel_token,
)
from .process_checker import check_and_warn_competing_processes

logger = get_logger(__name__)

# Global lock to prevent race conditions when multiple FoBs read input0.xyz simultaneously
_geometry_file_lock = threading.RLock()

# Safe file reader with retries to handle transient FS races
def _read_file_lines_with_retries(
    path: Path,
    *,
    attempts: int = 5,
    delay: float = 0.5,
    label: str = "file",
    min_lines: int = 1,
):
    """Read lines from a file with retry/backoff to avoid partial reads."""
    last_exc: Optional[Exception] = None
    for attempt in range(1, attempts + 1):
        try:
            with path.open('r') as handle:
                lines = handle.readlines()
            if len(lines) < min_lines:
                raise ValueError(f"{label} '{path}' is too short ({len(lines)} lines)")
            return lines
        except Exception as exc:
            last_exc = exc
            if attempt < attempts:
                wait = delay * attempt
                logger.warning(
                    "Read attempt %d/%d for %s '%s' failed: %s; retrying in %.1fs",
                    attempt,
                    attempts,
                    label,
                    path,
                    exc,
                    wait,
                )
                time.sleep(wait)
            else:
                logger.error(
                    "Read failed for %s '%s' after %d attempts: %s",
                    label,
                    path,
                    attempts,
                    exc,
                )
                raise last_exc

# ======================== Helper functions (extracted for flat architecture) ========================

def _stem(i: int, base: str = "input") -> str:
    """Generate filename stem for FoB index i."""
    return base if i == 1 else f"{base}{i}"


def _resolve_primary_source(raw_from, fallback: int) -> int:
    """Resolve the primary source index from raw_from (can be list/tuple/int)."""
    if isinstance(raw_from, (list, tuple, set)):
        for candidate in raw_from:
            try:
                parsed = int(str(candidate).strip())
            except (TypeError, ValueError):
                continue
            else:
                if parsed >= 0:
                    return parsed
        return fallback if fallback >= 0 else 0
    try:
        parsed = int(str(raw_from).strip())
    except (TypeError, ValueError):
        return fallback if fallback >= 0 else 0
    return parsed if parsed >= 0 else (fallback if fallback >= 0 else 0)


def _parse_dependency_indices(raw_from):
    """Parse dependency indices from various formats (list, string, int)."""
    deps: set[int] = set()
    if raw_from in (None, "", 0):
        return deps

    tokens = []
    if isinstance(raw_from, (list, tuple, set)):
        tokens = list(raw_from)
    else:
        text = str(raw_from)
        tokens = [tok for tok in re.split(r"[;,|]", text) if tok.strip()]

    if not tokens:
        tokens = [raw_from]

    for token in tokens:
        try:
            parsed = int(str(token).strip())
        except (TypeError, ValueError):
            continue
        if parsed > 0:
            deps.add(parsed)
    return deps


def _parse_geometry_wait(value: Any, default: Optional[float] = None) -> Optional[float]:
    """Convert CONTROL option to seconds; <=0 or 'inf' disables timeout."""
    try:
        from delfin.utils import normalize_str
        text = normalize_str(value)
    except Exception:
        text = ""

    if not text:
        return default
    if text in {"inf", "+inf", "infinite", "none", "disable", "disabled"}:
        return None

    try:
        parsed = float(text)
    except (TypeError, ValueError):
        return default

    if parsed <= 0:
        return None
    return parsed


def _wait_for_geometry_source(
    folder: Path,
    source_idx: int,
    *,
    label: str,
    fob_idx: int,
    timeout: Optional[float] = 900.0,
    failure_check: Optional[Callable[[], bool]] = None,
    completion_check: Optional[Callable[[], bool]] = None,
    poll_interval: float = 2.0,
) -> Path:
    """Block until the source geometry file exists or timeout expires."""
    if source_idx <= 1:
        candidate_name = "input.xyz"
    else:
        candidate_name = f"input{source_idx}.xyz"
    candidate = folder / candidate_name
    start = time.time()
    last_log = 0.0

    while True:
        if completion_check:
            try:
                completed = completion_check()
            except Exception:
                completed = False
        else:
            completed = True

        if candidate.exists() and completed:
            break

        elapsed = time.time() - start
        if not candidate.exists():
            if timeout is not None and elapsed >= timeout:
                raise FileNotFoundError(
                    f"[{label}] FoB {fob_idx}: geometry '{candidate_name}' not available after {int(elapsed)}s"
                )
        if failure_check and failure_check() and not candidate.exists():
            raise FileNotFoundError(
                f"[{label}] FoB {fob_idx}: source FoB {source_idx} failed and no geometry '{candidate_name}' was produced"
            )
        if elapsed - last_log >= 30.0:
            if candidate.exists() and not completed:
                msg = f"{candidate_name} (upstream ORCA still running)"
            else:
                msg = candidate_name
            logger.debug(
                "[%s] FoB %d waiting for %s (%.0fs elapsed)",
                label,
                fob_idx,
                msg,
                elapsed,
            )
            last_log = elapsed
        time.sleep(poll_interval)
    return candidate


# Covalent radii cache
COVALENT_RADII_CACHE: Dict[str, Optional[Dict[str, float]]] = {}

COVALENT_RADII_FALLBACK = {
    "H": 0.31, "He": 0.28,
    "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.25, "Co": 1.26, "Ni": 1.21, "Cu": 1.38, "Zn": 1.31,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
}


def load_covalent_radii(source="pyykko2009"):
    """Load single-bond covalent radii (Å) for H–Og using 'mendeleev'."""
    key = str(source).lower()
    if key in COVALENT_RADII_CACHE:
        return COVALENT_RADII_CACHE[key]
    try:
        from mendeleev import element
    except Exception:
        COVALENT_RADII_CACHE[key] = None
        return None

    attr = {"pyykko2009": "covalent_radius_pyykko", "cordero2008": "covalent_radius_cordero"}.get(
        key, "covalent_radius_pyykko"
    )
    radii = {}
    for Z in range(1, 119):
        el = element(Z)
        r = getattr(el, attr, None)
        if r is None:
            alt = "covalent_radius_cordero" if attr == "covalent_radius_pyykko" else "covalent_radius_pyykko"
            r = getattr(el, alt, None)
        if r is not None:
            radii[el.symbol] = float(r) / 100.0
        else:
            logger.warning("Missing covalent radius for %s (Z=%d); will use fallback.", el.symbol, Z)
    COVALENT_RADII_CACHE[key] = radii
    return radii


def _elem_from_label(label: str) -> str:
    """Extract the chemical symbol from 'Fe(1)' or 'Fe'."""
    m = re.match(r"([A-Za-z]{1,2})", label.strip())
    return m.group(1) if m else label.strip()


def _dist(a, b):
    """Euclidean distance for atom dicts with x,y,z."""
    return math.sqrt((a['x'] - b['x']) ** 2 + (a['y'] - b['y']) ** 2 + (a['z'] - b['z']) ** 2)


def _parse_xyz_atoms_with_indices(xyz_lines):
    """Parse atoms from a cleaned XYZ block (no count/comment lines)."""
    atoms = []
    for idx, line in enumerate(xyz_lines):
        ls = line.strip()
        if not ls or ls == '*':
            break
        parts = ls.split()
        if len(parts) < 4:
            continue
        raw_label = parts[0]
        elem = _elem_from_label(raw_label)
        try:
            x, y, z = map(float, parts[1:4])
        except ValueError:
            continue
        atoms.append({"line_idx": idx, "raw": raw_label, "elem": elem, "x": x, "y": y, "z": z})
    return atoms


def _first_coordination_sphere_indices(atoms, metal_indices, scale, radii_map):
    """Return a set of indices (in 'atoms') belonging to the first sphere of any metal."""
    def _rcov(sym: str) -> float:
        if radii_map and sym in radii_map:
            return float(radii_map[sym])
        return float(COVALENT_RADII_FALLBACK.get(sym, 1.20))

    first = set()
    for im in metal_indices:
        m = atoms[im]
        r_m = _rcov(m['elem'])
        for i, a in enumerate(atoms):
            if i == im:
                continue
            r_a = _rcov(a['elem'])
            cutoff = scale * (r_m + r_a)
            if _dist(m, a) <= cutoff:
                first.add(i)
    return first


def read_and_modify_file_OCCUPIER(from_index, output_file_path, charge, multiplicity,
                                  solvent, found_metals, metal_basisset, main_basisset,
                                  config, additions):
    """Build the ORCA input with per-atom NewGTO for metals and first coordination sphere."""
    xyz_file = "input.xyz" if from_index == 1 else f"input{from_index}.xyz"

    # Thread-safe file reading to prevent race conditions when multiple FoBs
    # access input0.xyz simultaneously
    with _geometry_file_lock:
        xyz_path = resolve_path(xyz_file)
        if not xyz_path.exists():
            if from_index == 1:
                alt_path = resolve_path("input0.xyz")
                if alt_path.exists():
                    logger.warning("Primary geometry '%s' missing; falling back to '%s'.", xyz_path, alt_path)
                    xyz_path = alt_path
                else:
                    error_msg = f"XYZ input file '{xyz_path}' not found (alt: '{alt_path}')"
                    logger.error(error_msg)
                    raise FileNotFoundError(error_msg)
            else:
                error_msg = f"XYZ input file '{xyz_path}' not found (from_index={from_index})"
                logger.error(error_msg)
                raise FileNotFoundError(error_msg)

        lines = _read_file_lines_with_retries(
            xyz_path,
            attempts=5,
            delay=0.5,
            label=f"XYZ geometry (from_index={from_index})",
            min_lines=2,
        )

    cleaned_xyz = normalize_xyz_body(lines[2:])  # skip count/comment lines
    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(cleaned_xyz, xyz_path)
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    atoms = _parse_xyz_atoms_with_indices(geom_lines)
    if not atoms:
        error_msg = f"No atoms parsed from XYZ file '{xyz_file}'"
        logger.error(error_msg)
        raise ValueError(error_msg)

    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes', 'true', '1', 'on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

    if enable_first:
        if sphere_scale_raw:
            sphere_scale = float(sphere_scale_raw)
            radii_all = None
        else:
            radii_all = load_covalent_radii(source=str(config.get('covalent_radii_source', 'pyykko2009')))
            sphere_scale = 1.20
    else:
        radii_all = None
        sphere_scale = float(sphere_scale_raw or 1.20)

    # Relativity token + AUX-JK token follow 3d/4d5 policy
    rel_token, aux_jk_token, _use_rel = select_rel_and_aux(found_metals or [], config)

    # Optional implicit solvent token
    implicit = ""
    if config.get('implicit_solvation_model') and solvent:
        implicit = f"{config['implicit_solvation_model']}({solvent})"
    elif config.get('implicit_solvation_model'):
        implicit = config['implicit_solvation_model']

    # Initial guess (trim accidental trailing text)
    initial_guess = (str(config.get('initial_guess', '')).split() or [''])[0]

    # Whether to add FREQ
    freq_flag = "FREQ" if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes' else ""

    # Build the '!' line
    tokens = ["!"]
    if qmmm_range:
        qmmm_method = str(config.get('qmmm_option', 'QM/XTB')).strip()
        tokens.append(qmmm_method)
    tokens.extend([
        str(config['functional']),
        rel_token,
        str(main_basisset),
        str(config.get('disp_corr', '')),
        str(config.get('ri_jkx', '')),
        aux_jk_token,
        implicit,
        str(config.get('geom_opt_OCCUPIER', config.get('geom_opt', ''))),
    ])
    if freq_flag:
        tokens.append(freq_flag)
    tokens.append(initial_guess)
    bang = " ".join(t for t in tokens if t).replace("  ", " ").strip()

    modified_lines = []
    if additions and "moinp" in additions.lower():
        bang += " MORead"
    modified_lines.append(bang + "\n")

    # Parallel / SCF controls
    modified_lines.append(f"%pal nprocs {config['PAL']} end\n")
    modified_lines.append(f"%scf maxiter {config['maxiter_occupier']} end\n")
    modified_lines.append(f"%maxcore {config['maxcore']}\n")

    # Extra blocks (e.g., %moinp, %scf BrokenSym, etc.)
    if additions and additions.strip():
        modified_lines.append(f"{additions.strip()}\n")

    # Add %freq block with temperature if FREQ is enabled
    if freq_flag:
        from .xyz_io import _build_freq_block
        freq_block = _build_freq_block(config)
        modified_lines.append(freq_block)

    # Start XYZ
    modified_lines.extend(build_qmmm_block(qmmm_range))
    modified_lines.append(f"* xyz {charge} {multiplicity}\n")

    # Metals found by symbol
    metal_syms = {m.strip().capitalize() for m in (found_metals or [])}
    metal_indices = [i for i, a in enumerate(atoms) if a['elem'].capitalize() in metal_syms]

    # First sphere indices (optional)
    first_sphere = set()
    if enable_first and metal_indices and metal_basisset:
        first_sphere = _first_coordination_sphere_indices(atoms, metal_indices, sphere_scale, radii_all)

    # Logging
    def _fmt_idx(i): return f"{i}({atoms[i]['elem']})"
    logger.info("Metals: %s | 1st sphere: %s",
                 ", ".join(map(_fmt_idx, metal_indices)) if metal_indices else "-",
                 ", ".join(map(_fmt_idx, sorted(first_sphere))) if first_sphere else "-")

    metal_line_set = {atoms[i]['line_idx'] for i in metal_indices}
    sphere_line_set = {atoms[i]['line_idx'] for i in first_sphere}

    # Write coordinates, appending per-atom NewGTO where needed
    for idx, line in enumerate(geom_lines):
        ls = line.strip()
        if not ls:
            continue
        apply_metal_basis = False
        if metal_basisset:
            if idx in metal_line_set:
                apply_metal_basis = True
            elif enable_first and idx in sphere_line_set:
                apply_metal_basis = True

        if apply_metal_basis:
            line = line.rstrip() + f'   NewGTO "{metal_basisset}" end'

        modified_lines.append(line if line.endswith("\n") else line + "\n")

    if not geom_lines or geom_lines[-1].strip() != '*':
        modified_lines.append("*\n")

    with open(output_file_path, 'w') as file:
        file.writelines(modified_lines)

    logger.info(f"Input file from '{xyz_file}' modified and saved as '{output_file_path}'")


# ======================== End of extracted helper functions ========================


def run_OCCUPIER():

    print("""
                      *******************
                      *     OCCUPIER    *
                      *******************
    """)

    global_log = os.environ.get("DELFIN_GLOBAL_LOG")
    if global_log:
        add_file_handler(global_log)
    occ_log_path = Path("occupier.log")
    add_file_handler(occ_log_path)
    logger.info("Occupier log attached at %s", occ_log_path.resolve())

    # Bootstrap global manager if running with scheduler-driven core allocation
    # Skip if running as subprocess (DELFIN_SUBPROCESS=1) to save startup time
    is_subprocess = str(os.environ.get("DELFIN_SUBPROCESS", "0")).lower() in ("1", "true", "yes", "on")
    if not is_subprocess:
        from delfin.global_manager import bootstrap_global_manager_from_env
        bootstrap_global_manager_from_env()

    # --------------------------- helpers ---------------------------

    def calculate_total_electrons(control_file_path):
        """Compute the total number of electrons from the XYZ referenced in CONTROL.txt."""
        atom_electrons = {
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

        control_path = resolve_path(control_file_path)
        config = read_control_file(str(control_path))
        input_file_entry = config.get("input_file")
        if not input_file_entry:
            logger.error("'input_file' not in CONTROL.txt.")
            return

        input_file_path = resolve_path(control_path.parent / input_file_entry)
        if not input_file_path.exists():
            logger.error(f"Input file '{input_file_path}' not found.")
            return

        total_electrons = 0
        try:
            with input_file_path.open('r') as input_file:
                lines = input_file.readlines()
                for line in lines[2:]:
                    parts = line.split()
                    if not parts:
                        continue
                    element = parts[0]
                    if element == "$":
                        continue
                    if element in atom_electrons:
                        total_electrons += atom_electrons[element]
                    else:
                        logger.warning(f"Unknown element '{element}' in line: {line.strip()}")
            return total_electrons
        except Exception as e:
            logger.error(f"Error reading file '{input_file_path}': {e}")

    def clean_xyz_block(lines):
        """Return coordinate lines without stray markers or empty rows."""
        return normalize_xyz_body(lines)

    # --------- covalent radii (all elements via mendeleev, with fallback) ---------

    def load_covalent_radii(source="pyykko2009"):
        """
        Load single-bond covalent radii (Å) for H–Og using 'mendeleev'.
        Returns None if the package is unavailable; a small fallback map is used then.
        """
        key = str(source).lower()
        if key in COVALENT_RADII_CACHE:
            return COVALENT_RADII_CACHE[key]
        try:
            from mendeleev import element
        except Exception as e:
            #logger.warning("Could not import 'mendeleev' (%s). Falling back to internal radii.", e)
            COVALENT_RADII_CACHE[key] = None
            return None

        attr = {"pyykko2009": "covalent_radius_pyykko", "cordero2008": "covalent_radius_cordero"}.get(
            key, "covalent_radius_pyykko"
        )
        radii = {}
        for Z in range(1, 119):
            el = element(Z)
            r = getattr(el, attr, None)
            if r is None:
                alt = "covalent_radius_cordero" if attr == "covalent_radius_pyykko" else "covalent_radius_pyykko"
                r = getattr(el, alt, None)
            if r is not None:
                # mendeleev reports radii in pm → convert to Å for downstream distance checks
                radii[el.symbol] = float(r) / 100.0
            else:
                logger.warning("Missing covalent radius for %s (Z=%d); will use fallback.", el.symbol, Z)
        COVALENT_RADII_CACHE[key] = radii
        return radii

    COVALENT_RADII_FALLBACK = {
        "H": 0.31, "He": 0.28,
        "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
        "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
        "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
        "Fe": 1.25, "Co": 1.26, "Ni": 1.21, "Cu": 1.38, "Zn": 1.31,
        "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
        "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
        "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    }

    COVALENT_RADII_CACHE: Dict[str, Optional[Dict[str, float]]] = {}


    def _elem_from_label(label: str) -> str:
        """Extract the chemical symbol from 'Fe(1)' or 'Fe'."""
        m = re.match(r"([A-Za-z]{1,2})", label.strip())
        return m.group(1) if m else label.strip()

    def _dist(a, b):
        """Euclidean distance for atom dicts with x,y,z."""
        return math.sqrt((a['x'] - b['x']) ** 2 + (a['y'] - b['y']) ** 2 + (a['z'] - b['z']) ** 2)

    def _parse_xyz_atoms_with_indices(xyz_lines):
        """
        Parse atoms from a cleaned XYZ block (no count/comment lines).
        Returns a list of dicts including the line index for later per-atom basis tagging.
        """
        atoms = []
        for idx, line in enumerate(xyz_lines):
            ls = line.strip()
            if not ls or ls == '*':
                break
            parts = ls.split()
            if len(parts) < 4:
                continue
            raw_label = parts[0]
            elem = _elem_from_label(raw_label)
            try:
                x, y, z = map(float, parts[1:4])
            except ValueError:
                continue
            atoms.append({"line_idx": idx, "raw": raw_label, "elem": elem, "x": x, "y": y, "z": z})
        return atoms

    def _first_coordination_sphere_indices(atoms, metal_indices, scale, radii_map):
        """
        Return a set of indices (in 'atoms') belonging to the first sphere of any metal:
            d(M-X) <= scale * (r_cov(M) + r_cov(X))
        """
        def _rcov(sym: str) -> float:
            if radii_map and sym in radii_map:
                return float(radii_map[sym])
            return float(COVALENT_RADII_FALLBACK.get(sym, 1.20))

        first = set()
        for im in metal_indices:
            m = atoms[im]
            r_m = _rcov(m['elem'])
            for i, a in enumerate(atoms):
                if i == im:
                    continue
                r_a = _rcov(a['elem'])
                cutoff = scale * (r_m + r_a)
                if _dist(m, a) <= cutoff:
                    first.add(i)
        return first

    # ---------------- ORCA input writer (per-atom NewGTO) ----------------

    def read_and_modify_file_OCCUPIER(from_index, output_file_path, charge, multiplicity,
                                      solvent, found_metals, metal_basisset, main_basisset,
                                      config, additions):
        """
        Build the ORCA input with:
          - global method line '!' using main_basisset and aux-JK (via select_rel_and_aux),
          - per-atom NewGTO "metal_basisset" for metals (always),
          - optional per-atom NewGTO for the first coordination sphere.
        """
        xyz_file = "input.xyz" if from_index == 1 else f"input{from_index}.xyz"

        # Thread-safe file reading to prevent race conditions when multiple FoBs
        # access input0.xyz simultaneously
        with _geometry_file_lock:
            xyz_path = resolve_path(xyz_file)
            if not xyz_path.exists():
                if from_index == 1:
                    # Some workflows maintain only input0.xyz; fall back gracefully.
                    alt_path = resolve_path("input0.xyz")
                    if alt_path.exists():
                        logger.warning(
                            "Primary geometry '%s' missing; falling back to '%s'.",
                            xyz_path,
                            alt_path,
                        )
                        xyz_path = alt_path
                    else:
                        error_msg = f"XYZ input file '{xyz_path}' not found (alt: '{alt_path}')"
                        logger.error(error_msg)
                        raise FileNotFoundError(error_msg)
                else:
                    error_msg = f"XYZ input file '{xyz_path}' not found (from_index={from_index})"
                    logger.error(error_msg)
                    raise FileNotFoundError(error_msg)

            lines = _read_file_lines_with_retries(
                xyz_path,
                attempts=5,
                delay=0.5,
                label=f"XYZ geometry (from_index={from_index})",
                min_lines=2,
            )

        cleaned_xyz = clean_xyz_block(lines[2:])  # skip count/comment lines
        geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(cleaned_xyz, xyz_path)
        _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
        atoms = _parse_xyz_atoms_with_indices(geom_lines)
        if not atoms:
            error_msg = f"No atoms parsed from XYZ file '{xyz_file}'"
            logger.error(error_msg)
            raise ValueError(error_msg)

        enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes', 'true', '1', 'on')
        sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

        if enable_first:
            if sphere_scale_raw:
                sphere_scale = float(sphere_scale_raw)
                radii_all = None
            else:
                radii_all = load_covalent_radii(source=str(config.get('covalent_radii_source', 'pyykko2009')))
                sphere_scale = 1.20
        else:
            radii_all = None
            sphere_scale = float(sphere_scale_raw or 1.20)

        # Relativity token + AUX-JK token follow 3d/4d5 policy
        rel_token, aux_jk_token, _use_rel = select_rel_and_aux(found_metals or [], config)

        # Optional implicit solvent token
        implicit = ""
        if config.get('implicit_solvation_model') and solvent:
            implicit = f"{config['implicit_solvation_model']}({solvent})"
        elif config.get('implicit_solvation_model'):
            implicit = config['implicit_solvation_model']

        # Initial guess (trim accidental trailing text)
        initial_guess = (str(config.get('initial_guess', '')).split() or [''])[0]

        # Whether to add FREQ
        freq_flag = "FREQ" if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes' else ""

        # Build the '!' line
        tokens = ["!"]
        if qmmm_range:
            qmmm_method = str(config.get('qmmm_option', 'QM/XTB')).strip()
            tokens.append(qmmm_method)
        tokens.extend([
            str(config['functional']),
            rel_token,                       # '' or ZORA/X2C/DKH
            str(main_basisset),
            str(config.get('disp_corr', '')),
            str(config.get('ri_jkx', '')),
            aux_jk_token,                    # '' or def2/J | SARC/J
            implicit,
            str(config.get('geom_opt_OCCUPIER', config.get('geom_opt', ''))),
        ])
        if freq_flag:
            tokens.append(freq_flag)
        tokens.append(initial_guess)
        bang = " ".join(t for t in tokens if t).replace("  ", " ").strip()

        modified_lines = []
        if additions and "moinp" in additions.lower():
            bang += " MORead"
        modified_lines.append(bang + "\n")

        # Parallel / SCF controls
        modified_lines.append(f"%pal nprocs {config['PAL']} end\n")
        modified_lines.append(f"%scf maxiter {config['maxiter_occupier']} end\n")
        modified_lines.append(f"%maxcore {config['maxcore']}\n")

        # Extra blocks (e.g., %moinp, %scf BrokenSym, etc.)
        if additions and additions.strip():
            modified_lines.append(f"{additions.strip()}\n")

        # Add %freq block with temperature if FREQ is enabled
        if freq_flag:
            from .xyz_io import _build_freq_block
            freq_block = _build_freq_block(config)
            modified_lines.append(freq_block)

        # Start XYZ
        modified_lines.extend(build_qmmm_block(qmmm_range))
        modified_lines.append(f"* xyz {charge} {multiplicity}\n")

        # Metals found by symbol
        metal_syms = {m.strip().capitalize() for m in (found_metals or [])}
        metal_indices = [i for i, a in enumerate(atoms) if a['elem'].capitalize() in metal_syms]

        # First sphere indices (optional)
        first_sphere = set()
        if enable_first and metal_indices and metal_basisset:
            first_sphere = _first_coordination_sphere_indices(atoms, metal_indices, sphere_scale, radii_all)

        # Logging
        def _fmt_idx(i): return f"{i}({atoms[i]['elem']})"
        logger.info("Metals: %s | 1st sphere: %s",
                     ", ".join(map(_fmt_idx, metal_indices)) if metal_indices else "-",
                     ", ".join(map(_fmt_idx, sorted(first_sphere))) if first_sphere else "-")

        metal_line_set = {atoms[i]['line_idx'] for i in metal_indices}
        sphere_line_set = {atoms[i]['line_idx'] for i in first_sphere}

        # Write coordinates, appending per-atom NewGTO where needed
        for idx, line in enumerate(geom_lines):
            ls = line.strip()
            if not ls:
                continue
            apply_metal_basis = False
            if metal_basisset:
                if idx in metal_line_set:
                    apply_metal_basis = True
                elif enable_first and idx in sphere_line_set:
                    apply_metal_basis = True

            if apply_metal_basis:
                line = line.rstrip() + f'   NewGTO "{metal_basisset}" end'

            modified_lines.append(line if line.endswith("\n") else line + "\n")

        if not geom_lines or geom_lines[-1].strip() != '*':
            modified_lines.append("*\n")

        with open(output_file_path, 'w') as file:
            file.writelines(modified_lines)

        logger.info(f"Input file from '{xyz_file}' modified and saved as '{output_file_path}'")

    # ------------------------ output parsers ------------------------

    def find_FSPE(filename):
        """Scan ORCA output bottom-up for 'FINAL SINGLE POINT ENERGY'."""
        search_text = "FINAL SINGLE POINT ENERGY   "
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError:
            logger.warning(f"File {filename} not found. Skipping.")
            return None

        for line in reversed(lines):
            idx = line.find(search_text)
            if idx != -1:
                energy_str = line[idx + len(search_text):].strip().split()[0]
                try:
                    return float(energy_str)
                except ValueError:
                    logger.error(f"Could not convert '{energy_str}' to float.")
                    continue
        logger.warning("No FINAL SINGLE POINT ENERGY found in the file.")
        return None

    def find_G(filename):
        """Scan ORCA output for 'Final Gibbs free energy ...'."""
        search_text = "Final Gibbs free energy         ... "
        try:
            with open(filename, 'r') as file:
                for line in file:
                    idx = line.find(search_text)
                    if idx != -1:
                        energy_str = line[idx + len(search_text):].strip().split()[0]
                        try:
                            return float(energy_str)
                        except ValueError:
                            logger.error(f"Could not convert '{energy_str}' to float.")
                            return None
        except FileNotFoundError:
            logger.warning(f"File {filename} not found. Skipping.")
        return None

    # --------------------------- main flow -------------------------

    try:
        start_time = time.time()
        control_file_path = resolve_path("CONTROL.txt")
        config = OCCUPIER_parser(str(control_file_path))

        input_file_entry = config.get("input_file")
        if not input_file_entry:
            logger.error("Missing 'input_file' in the configuration.")
            return
        input_file_path = resolve_path(control_file_path.parent / input_file_entry)
        if not input_file_path.exists():
            logger.error(f"Input file '{input_file_path}' does not exist.")
            return

        input_file = str(input_file_path)

        solvent = config.get("solvent", None)
        total_electrons = calculate_total_electrons(str(control_file_path))
        charge = int(config.get("charge", 0))
        total_electrons = (total_electrons or 0) - charge
        is_even = (total_electrons % 2 == 0)

        # Detect metals for policy & per-atom tagging
        metals = search_transition_metals(input_file)

        # Select orbital bases based on 3d/4d5 policy (from utils)
        main_basisset, metal_basisset = set_main_basisset(metals, config)
        # If no metals present, do not apply per-atom overrides
        if not metals:
            metal_basisset = None

        species_delta = infer_species_delta(Path.cwd())
        seq_bundle = resolve_sequences_for_delta(config, species_delta)

        method_token = str(config.get("OCCUPIER_method", "auto") or "auto").strip().lower()
        if method_token == "auto" and seq_bundle:
            try:
                remove_existing_sequence_blocks(control_file_path, force=True)
                append_sequence_overrides(control_file_path, seq_bundle)
            except Exception as exc:  # noqa: BLE001
                logger.debug("[occupier] Failed to refresh auto sequence overrides: %s", exc)

        # Choose which sequence to run
        seq_key = "even_seq" if is_even else "odd_seq"
        sequence = seq_bundle.get(seq_key) or config.get(seq_key, [])
        if not sequence:
            logger.error(f"No sequence found under '{seq_key}' in CONTROL.txt.")
            return

        # Helper to build filename stems: 1 -> "input", n>=2 -> "input{n}"
        def _stem(i: int, base: str = "input") -> str:
            return base if i == 1 else f"{base}{i}"

        def _resolve_primary_source(raw_from, fallback: int) -> int:
            if isinstance(raw_from, (list, tuple, set)):
                for candidate in raw_from:
                    try:
                        parsed = int(str(candidate).strip())
                    except (TypeError, ValueError):
                        continue
                    else:
                        if parsed >= 0:
                            return parsed
                return fallback if fallback >= 0 else 0
            try:
                parsed = int(str(raw_from).strip())
            except (TypeError, ValueError):
                return fallback if fallback >= 0 else 0
            return parsed if parsed >= 0 else (fallback if fallback >= 0 else 0)

        def _parse_dependency_indices(raw_from):
            deps: set[int] = set()
            if raw_from in (None, "", 0):
                return deps

            tokens = []
            if isinstance(raw_from, (list, tuple, set)):
                tokens = list(raw_from)
            else:
                text = str(raw_from)
                tokens = [tok for tok in re.split(r"[;,|]", text) if tok.strip()]

            if not tokens:
                tokens = [raw_from]

            for token in tokens:
                try:
                    parsed = int(str(token).strip())
                except (TypeError, ValueError):
                    continue
                if parsed > 0:
                    deps.add(parsed)
            return deps

        def _resolve_pal_jobs_value(config_dict) -> int:
            raw = config_dict.get('pal_jobs')
            try:
                parsed = int(str(raw).strip()) if raw not in (None, "") else 0
            except (TypeError, ValueError):
                parsed = 0
            if parsed <= 0:
                try:
                    total = int(str(config_dict.get('PAL', 1)).strip())
                except (TypeError, ValueError):
                    total = 1
                parsed = max(1, min(4, max(1, total // 2)))
            return parsed

        # Energy extractor: FSPE unless frequency run is requested
        use_gibbs = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes'
        finder = find_G if use_gibbs else find_FSPE

        # Parallel OCCUPIER execution
        resolved_pal_jobs = _resolve_pal_jobs_value(config)

        def _sequence_parallel_width(seq) -> int:
            if not seq:
                return 0
            known = {int(entry["index"]) for entry in seq if "index" in entry}
            deps_map = {
                idx: {
                    dep for dep in _parse_dependency_indices(entry.get("from", idx - 1))
                    if dep in known
                }
                for entry, idx in ((entry, int(entry["index"])) for entry in seq if "index" in entry)
            }

            completed: set[int] = set()
            remaining = set(known)
            max_width = 0
            guard = 0

            while remaining and guard <= len(known) * 2:
                ready = {idx for idx in remaining if deps_map[idx] <= completed}
                if not ready:
                    break
                max_width = max(max_width, len(ready))
                completed.update(ready)
                remaining -= ready
                guard += 1

            if remaining:
                return max_width or 1

            return max(max_width, 1)

        sequence_width = _sequence_parallel_width(sequence)
        parallel_setting = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
        orca_parallel_strategy = str(config.get('orca_parallel_strategy', 'auto')).strip().lower()

        # Interpret parallel strategy:
        # - 'auto': Use MPI if available → FoBs can run in parallel
        # - 'threads': Force OpenMP-only (no MPI) → FoBs can STILL run in parallel!
        # - 'serial': Force truly sequential execution → only 1 FoB at a time
        strategy_forces_serial = orca_parallel_strategy == 'serial'

        parallel_allowed = (
            not strategy_forces_serial
            and (
                parallel_setting == 'enable'
                or (parallel_setting == 'auto' and sequence_width > 1)
            )
        )
        effective_pal_jobs = max(1, min(resolved_pal_jobs, sequence_width)) if parallel_allowed else 1

        if strategy_forces_serial and parallel_setting != 'disable':
            logger.info(
                "[occupier] ORCA parallel strategy=%s → forcing sequential OCCUPIER scheduling",
                orca_parallel_strategy,
            )
        elif orca_parallel_strategy == 'threads':
            logger.info(
                "[occupier] ORCA parallel strategy=threads → using OpenMP (no MPI), "
                "but FoBs can still run in parallel"
            )

        logger.info(
            "Scheduling %d OCCUPIER FoBs (parallel=%s, pal_jobs=%s → resolved=%d, width=%d)",
            len(sequence),
            parallel_setting,
            config.get('pal_jobs'),
            effective_pal_jobs,
            sequence_width,
        )

        known_indices = {int(entry["index"]) for entry in sequence if "index" in entry}
        force_sequential = parallel_setting == "disable"
        previous_idx: Optional[int] = None
        dependencies_map = {
            int(entry["index"]): {
                f"occ_{dep}"
                for dep in _parse_dependency_indices(entry.get("from", entry["index"] - 1))
                if dep in known_indices
            }
            for entry in sequence
        }
        if force_sequential:
            for entry in sequence:
                idx = int(entry["index"])
                if previous_idx is not None and idx != previous_idx:
                    dependencies_map[idx].add(f"occ_{previous_idx}")
                previous_idx = idx

        # Always use the global job manager (even in subprocess mode)
        global_mgr = get_global_manager()
        manager_initialized = False
        try:
            manager_initialized = global_mgr.is_initialized()
        except Exception:  # noqa: BLE001
            manager_initialized = False

        if not manager_initialized:
            # Ensure a global pool exists so parallel scheduling can proceed
            global_mgr.initialize(config)

        manager = _WorkflowManager(config, label="occupier_core", max_jobs_override=effective_pal_jobs)

        if effective_pal_jobs <= 1 and manager.pool.max_concurrent_jobs != 1:
            manager.pool.max_concurrent_jobs = 1
            manager.max_jobs = 1
            manager._sync_parallel_flag()

        fspe_results: dict[int, float | None] = {}
        results_lock = threading.Lock()

        OK = "ORCA TERMINATED NORMALLY"
        recalc = str(os.environ.get("DELFIN_RECALC", "0")).lower() in ("1", "true", "yes", "on")
        geometry_wait_timeout = _parse_geometry_wait(
            config.get("occupier_geometry_wait_s"), default=None
        )

        def _has_ok_marker(path: str) -> bool:
            candidate = resolve_path(path)
            if not candidate.exists():
                return False
            try:
                # Check file size to avoid reading incomplete files
                if candidate.stat().st_size < 100:  # Files with OK marker should be larger
                    return False
                with candidate.open("r", encoding="utf-8", errors="replace") as _f:
                    content = _f.read()
                    return OK in content
            except Exception as e:
                logger.debug("[recalc] could not check %s (%s) -> will run", path, e)
                return False

        freq_enabled = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes'
        pass_wf_enabled = str(config.get('pass_wavefunction', 'no')).strip().lower() in ('yes', 'true', '1', 'on', 'y')
        raw_apm = config.get("approximate_spin_projection_APMethod")
        apm = str(raw_apm).strip() if raw_apm not in (None, "", 0, "0") else ""

        def make_work(entry_dict: dict) -> callable:
            idx = int(entry_dict["index"])
            multiplicity = entry_dict["m"]
            bs = entry_dict.get("BS", "")
            raw_from = entry_dict.get("from", idx - 1)
            src_idx = _resolve_primary_source(raw_from, idx - 1)

            stem = _stem(idx)
            inp = f"{stem}.inp"
            out = f"output{'' if idx == 1 else idx}.out"

            def _work(cores: int) -> None:
                if recalc and _has_ok_marker(out):
                    logger.info("[recalc] Skipping ORCA for %s; found '%s'.", out, OK)
                    parsed_val = finder(out)
                    with results_lock:
                        fspe_results[idx] = parsed_val
                    return

                parts: list[str] = []

                if pass_wf_enabled:
                    gbw_candidate = "input.gbw" if src_idx == 1 else f"input{src_idx}.gbw"
                    gbw_path = resolve_path(gbw_candidate)
                    if gbw_path.exists():
                        parts.append(f'%moinp "{gbw_path}"')
                    else:
                        logger.info(
                            "No GBW found for from=%s (%s) – starting with standard guess.",
                            src_idx,
                            gbw_path,
                        )

                if bs:
                    scf_lines = [f"  BrokenSym {bs}"]
                    if not freq_enabled and apm:
                        scf_lines.append(f"  APMethod {apm}")
                    parts.append("%scf\n" + "\n".join(scf_lines) + "\nend")

                additions = "\n".join(parts)

                workdir = Path.cwd()
                label = workdir.name or "occupier_core"
                def _source_failed() -> bool:
                    effective_idx = src_idx
                    if effective_idx <= 1:
                        return False
                    with results_lock:
                        return (
                            effective_idx in fspe_results
                            and fspe_results[effective_idx] is None
                        )

                def _source_completed() -> bool:
                    effective_idx = src_idx
                    if effective_idx <= 1:
                        return True
                    with results_lock:
                        return effective_idx in fspe_results

                _wait_for_geometry_source(
                    workdir,
                    src_idx,
                    label=label,
                    fob_idx=idx,
                    timeout=geometry_wait_timeout,
                    failure_check=_source_failed,
                    completion_check=_source_completed,
                )

                read_and_modify_file_OCCUPIER(
                    src_idx,
                    inp,
                    charge,
                    multiplicity,
                    solvent,
                    metals,
                    metal_basisset,
                    main_basisset,
                    config,
                    additions,
                )

                if not Path(inp).exists():
                    raise RuntimeError(f"Failed to create OCCUPIER input '{inp}'")

                _update_pal_block(inp, cores)

                # Second check right before execution (race condition protection)
                if recalc and _has_ok_marker(out):
                    logger.info("[recalc] Skipping ORCA for %s; completed by another process.", out)
                    parsed_val = finder(out)
                    with results_lock:
                        fspe_results[idx] = parsed_val
                    return

                # Check for competing ORCA processes that might block execution
                check_and_warn_competing_processes(inp)

                # Run ORCA and check if it succeeded
                orca_success = run_orca(inp, out)
                if not orca_success:
                    error_msg = f"ORCA failed for {inp}. Check {out} for errors."
                    logger.error(error_msg)
                    # Mark this FoB as failed by storing None
                    with results_lock:
                        fspe_results[idx] = None
                    raise RuntimeError(error_msg)

                parsed_val = finder(out)
                with results_lock:
                    fspe_results[idx] = parsed_val

            return _work

        try:
            for entry in sequence:
                idx = int(entry["index"])
                job = WorkflowJob(
                    job_id=f"occ_{idx}",
                    work=make_work(entry),
                    description=f"FoB index {idx}",
                    dependencies=dependencies_map.get(idx, set()),
                )
                cores_bounds = manager.derive_core_bounds(hint=job.description)
                job.cores_min, job.cores_optimal, job.cores_max = cores_bounds
                manager.add_job(job)

            dynamic_slots = determine_effective_slots(
                manager.total_cores,
                manager._jobs.values(),
                effective_pal_jobs,
                sequence_width,
            )
            if dynamic_slots != manager.pool.max_concurrent_jobs:
                logger.info(
                    "[occupier_core] Updating parallel slots to %d (width=%d, requested=%d)",
                    dynamic_slots,
                    sequence_width,
                    effective_pal_jobs,
                )
                manager.pool.max_concurrent_jobs = dynamic_slots
                manager.max_jobs = dynamic_slots
                manager._sync_parallel_flag()

            manager.run()
        finally:
            if manager is not None:
                manager.shutdown()

        fspe_values = [fspe_results.get(int(entry["index"])) for entry in sequence]

        duration = time.time() - start_time
        generate_summary_report_OCCUPIER(duration, fspe_values, is_even, charge, solvent, config, main_basisset, sequence)
        logger.info("Summary report generated and saved as 'OCCUPIER.txt'")

    except Exception as e:
        logger.error(f"An error occurred in the OCCUPIER workflow: {e}")
