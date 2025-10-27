import json
import logging
import math
import re
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from delfin.common.orca_blocks import OrcaInputBuilder, collect_output_blocks, resolve_maxiter

# import canonical selection helpers from utils
from .utils import set_main_basisset, select_rel_and_aux

# Cache QMMM splits once detected so subsequent geometries (e.g. ORCA outputs without '$')
# continue to receive the QM/XTB setup.
_QMMM_CACHE: Dict[Tuple[str, ...], Tuple[int, int]] = {}
_QMMM_CACHE_LOCK = threading.Lock()
_QMMM_CACHE_BASENAME = ".qmmm_cache.json"


def _cache_key(signature: Tuple[str, ...]) -> str:
    """Convert a geometry signature into a stable string key."""
    return "|".join(signature)


def _normalise_source_path(source_path: Optional[Path]) -> Optional[Path]:
    """Resolve the provided source path if possible."""
    if source_path is None:
        return None
    try:
        return source_path.resolve()
    except Exception:
        return source_path


def _candidate_cache_paths(source_path: Optional[Path]) -> List[Path]:
    """
    Return cache file locations to probe, starting from the geometry directory
    and walking up the parent chain.
    """
    if source_path is None:
        return []
    base = source_path if source_path.is_dir() else source_path.parent
    base = _normalise_source_path(base)
    if base is None:
        return []
    candidates: List[Path] = []
    seen: Set[Path] = set()
    current = base
    for _ in range(8):  # guard against deep/recursive paths
        if current in seen:
            break
        seen.add(current)
        candidates.append(current / _QMMM_CACHE_BASENAME)
        parent = current.parent
        if parent == current:
            break
        current = parent
    return candidates


def _load_qmmm_signature_from_disk(signature: Tuple[str, ...], source_path: Optional[Path]) -> Optional[Tuple[int, int]]:
    """Attempt to load a persisted QM/XTB split for the given signature."""
    key = _cache_key(signature)
    for cache_file in _candidate_cache_paths(source_path):
        try:
            with cache_file.open("r", encoding="utf-8") as fh:
                data = json.load(fh)
        except FileNotFoundError:
            continue
        except Exception:
            continue
        raw = data.get(key)
        if isinstance(raw, dict):
            start = raw.get("start")
            end = raw.get("end")
            if start is not None and end is not None:
                try:
                    return int(start), int(end)
                except (TypeError, ValueError):
                    continue
        if isinstance(raw, list) and len(raw) == 2:
            try:
                return int(raw[0]), int(raw[1])
            except (TypeError, ValueError):
                continue
        # fallback: use directory-level default mapping if present
        default_raw = data.get("__default__")
        if isinstance(default_raw, list) and len(default_raw) == 2:
            try:
                return int(default_raw[0]), int(default_raw[1])
            except (TypeError, ValueError):
                continue
        if isinstance(default_raw, dict):
            start = default_raw.get("start")
            end = default_raw.get("end")
            if start is not None and end is not None:
                try:
                    return int(start), int(end)
                except (TypeError, ValueError):
                    continue
    return None


def _persist_qmmm_signature(signature: Tuple[str, ...], qmmm_range: Tuple[int, int], source_path: Optional[Path]) -> None:
    """Store the detected QM/XTB split alongside the geometry and its parent directory."""
    if source_path is None:
        return
    key = _cache_key(signature)
    payload = [int(qmmm_range[0]), int(qmmm_range[1])]
    base = source_path if source_path.is_dir() else source_path.parent
    base = _normalise_source_path(base)
    if base is None:
        return
    targets = {base}
    parent = base.parent
    if parent != base:
        targets.add(parent)
    for target in targets:
        cache_file = target / _QMMM_CACHE_BASENAME
        try:
            with _QMMM_CACHE_LOCK:
                try:
                    with cache_file.open("r", encoding="utf-8") as fh:
                        data = json.load(fh)
                except FileNotFoundError:
                    data = {}
                except Exception:
                    data = {}
                data[key] = payload
                data.setdefault("__default__", payload)
                with cache_file.open("w", encoding="utf-8") as fh:
                    json.dump(data, fh)
        except Exception:
            continue


def _geometry_signature(lines: List[str]) -> Optional[Tuple[str, ...]]:
    """
    Build a signature for a geometry based solely on element ordering.
    This stays constant across optimizations, allowing us to reuse cached
    QM/XTB splits even if coordinates change.
    """
    elements: List[str] = []
    for ln in lines:
        stripped = ln.strip()
        if not stripped or stripped in {"*", "$"}:
            continue
        parts = stripped.split()
        if not parts:
            continue
        elements.append(_elem_from_label(parts[0]))
    if not elements:
        return None
    return tuple(elements)

# -------------------------------------------------------------------------
# generic IO helpers (unchanged API)
# -------------------------------------------------------------------------
def write_to_file(lines: List[str], output_file_path: str) -> None:
    with open(output_file_path, 'w') as file:
        for line in lines:
            file.write(line + '\n')
    logging.info(f"Lines written to '{output_file_path}'")

def modify_file(file_path: str) -> None:
    with open(file_path, 'r') as file:
        lines = file.readlines()
    modified_lines = lines[2:-3]
    with open(file_path, 'w') as file:
        for line in modified_lines:
            file.write(line)
    logging.info(f"File '{file_path}' modified")

def extract_orbital_energies(file_path: str) -> List[str]:
    with open(file_path, 'r') as file:
        lines = file.readlines()
    start_marker = "ORBITAL ENERGIES"
    end_marker = "MULLIKEN POPULATION ANALYSIS"
    extracting = False
    extracted_lines = []
    for line in lines:
        if start_marker in line:
            extracting = True
            continue
        if extracting and end_marker in line:
            extracting = False
            break
        if extracting:
            extracted_lines.append(line.strip())
    logging.info(f"Extracted orbital energies from '{file_path}'")
    return extracted_lines

def find_homo_index_and_extract_lines(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    last_occurrence_index = -1
    for i, line in enumerate(lines):
        if "2.0000" in line:
            last_occurrence_index = i
    if last_occurrence_index == -1:
        logging.warning("No occurrence of '2.0000' found in the file.")
        return -1, []
    else:
        homo_index = last_occurrence_index - 1
        logging.info(f"HOMO index found at {homo_index}")
        start_index = max(homo_index - 10, 0)
        end_index = min(homo_index + 11, len(lines))
        mos_for_print = [line.strip() for line in lines[start_index:end_index]]
        return homo_index, mos_for_print

def modify_file2(target_file, header, footer):
    with open(target_file, "r") as f:
        content = f.read()
    with open(target_file, "w") as f:
        f.write(header + content + footer)

# -------------------------------------------------------------------------
# geometry + basis helpers
# -------------------------------------------------------------------------

# Compact fallback radii (Å); a full table can be loaded via 'mendeleev' if available
_COVALENT_RADII_FALLBACK: Dict[str, float] = {
    "H": 0.31, "He": 0.28,
    "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.25, "Co": 1.26, "Ni": 1.21, "Cu": 1.38, "Zn": 1.31,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
}

_RADII_CACHE: Dict[str, Optional[Dict[str, float]]] = {}


def _load_covalent_radii(source: str = "pyykko2009") -> Optional[Dict[str, float]]:
    """Return {symbol: radius Å} using 'mendeleev' if available; otherwise None."""
    key = str(source).lower()
    if key in _RADII_CACHE:
        return _RADII_CACHE[key]
    try:
        from mendeleev import element
    except Exception as e:
        #logging.info("mendeleev not available (%s) – using fallback radii.", e)
        _RADII_CACHE[key] = None
        return None
    attr = {
        "pyykko2009": "covalent_radius_pyykko",
        "cordero2008": "covalent_radius_cordero",
    }.get(key, "covalent_radius_pyykko")
    radii: Dict[str, float] = {}
    for Z in range(1, 119):
        el = element(Z)
        r = getattr(el, attr, None)
        if r is None:
            alt = "covalent_radius_cordero" if attr == "covalent_radius_pyykko" else "covalent_radius_pyykko"
            r = getattr(el, alt, None)
        if r is not None:
            # mendeleev radii are reported in picometers; convert to Å to stay consistent
            radii[el.symbol] = float(r) / 100.0
    _RADII_CACHE[key] = radii
    return radii

def _elem_from_label(label: str) -> str:
    """Extract chemical symbol from lines like 'Fe(1)' or 'Fe'."""
    m = re.match(r"([A-Za-z]{1,2})", label.strip())
    return m.group(1) if m else label.strip()

def _dist(a, b):
    return math.sqrt((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2)

def _parse_xyz_atoms(xyz_lines: List[str]):
    """
    Parse atoms from a list of coordinate lines (no count/comment).
    Returns a list of dicts with coords, element and original line index.
    """
    atoms = []
    for idx, line in enumerate(xyz_lines):
        ls = line.strip()
        if not ls or ls == '*':
            break
        parts = ls.split()
        if len(parts) < 4:
            continue
        raw = parts[0]
        elem = _elem_from_label(raw)
        try:
            x, y, z = map(float, parts[1:4])
        except ValueError:
            continue
        atoms.append({"line_idx": idx, "elem": elem, "x": x, "y": y, "z": z})
    return atoms


def split_qmmm_sections(
    coord_lines: List[str],
    source_path: Optional[Path] = None,
) -> Tuple[List[str], Optional[Tuple[int, int]], bool]:
    """
    Split coordinate lines into QM/XTB sections using a line that contains only '$' as separator.

    Returns a tuple ``(all coordinate lines without separators, qmmm_range, explicit)``, where
    ``qmmm_range`` is the ``(start, end)`` tuple for QMAtoms (or ``None`` if no split exists) and
    ``explicit`` is ``True`` when the separator was present in the current file (``False`` when the
    split is restored from cache).
    """
    qm_lines: List[str] = []
    mm_lines: List[str] = []
    separator_seen = False

    for raw in coord_lines:
        stripped = raw.strip()
        if not stripped or stripped == "*":
            continue
        if stripped == "$":
            if not separator_seen:
                separator_seen = True
            continue

        line = raw if raw.endswith("\n") else raw + "\n"
        if separator_seen:
            mm_lines.append(line)
        else:
            qm_lines.append(line)

    combined = qm_lines + mm_lines

    signature = _geometry_signature(combined)

    if separator_seen and qm_lines:
        qmmm_range = (0, len(qm_lines) - 1)
        if signature:
            with _QMMM_CACHE_LOCK:
                _QMMM_CACHE[signature] = qmmm_range
            _persist_qmmm_signature(signature, qmmm_range, source_path)
        return combined, qmmm_range, True

    if signature:
        cached: Optional[Tuple[int, int]] = None
        with _QMMM_CACHE_LOCK:
            cached = _QMMM_CACHE.get(signature)
        if cached:
            return combined, cached, False

        disk_cached = _load_qmmm_signature_from_disk(signature, source_path)
        if disk_cached:
            with _QMMM_CACHE_LOCK:
                _QMMM_CACHE[signature] = disk_cached
            return combined, disk_cached, False

    return combined, None, False


def build_qmmm_block(qmmm_range: Optional[Tuple[int, int]]) -> List[str]:
    """
    Build the %QMMM block for ORCA when a QM/XTB split is present.
    """
    if not qmmm_range:
        return []

    start, end = qmmm_range
    if end < start:
        return []

    return [
        "%QMMM\n",
        f"  QMATOMS {{{start}:{end}}}\n",
        "END\n",
        "END\n",
    ]

def _rcov(sym: str, radii_map: Optional[Dict[str, float]]) -> float:
    if radii_map and sym in radii_map:
        return float(radii_map[sym])
    return float(_COVALENT_RADII_FALLBACK.get(sym, 1.20))

def _first_sphere_indices(atoms, metal_indices, scale, radii_map):
    """Return a set of indices of atoms that belong to the first sphere of any metal."""
    first = set()
    for im in metal_indices:
        m = atoms[im]
        r_m = _rcov(m["elem"], radii_map)
        for i, a in enumerate(atoms):
            if i == im:
                continue
            r_a = _rcov(a["elem"], radii_map)
            cutoff = scale * (r_m + r_a)
            if _dist(m, a) <= cutoff:
                first.add(i)
    return first

def _implicit_token(config, solvent):
    """Build implicit solvent token."""
    mdl = str(config.get('implicit_solvation_model','') or '').strip()
    if not mdl:
        return ""
    return f"{mdl}({solvent})" if solvent else mdl


def _ensure_qmmm_implicit_model(
    config: Dict[str, Any],
    qmmm_range: Optional[Tuple[int, int]],
    qmmm_explicit: bool = False,
) -> None:
    """
    Ensure the implicit solvation model is compatible with QM/MM separators ('$').
    If CPCM is requested while QM/MM is active, automatically switch to ALPB
    and emit a warning; other models are left untouched.
    """
    if not qmmm_range:
        return
    raw_model = config.get("implicit_solvation_model")
    if raw_model is None:
        return
    model = str(raw_model).strip()
    if not model:
        return
    if not qmmm_explicit:
        # Split restored from cache – respect the user-provided model.
        logging.debug(
            "QM/MM split inferred from cache; keeping implicit solvation model '%s'.",
            model,
        )
        return
    if model.upper() == "CPCM":
        logging.warning('Detected "$" separator: CPCM is incompatible with QM/MM. Switching implicit_solvation_model to ALPB.')
        config["implicit_solvation_model"] = "ALPB"

def _build_freq_block(config):
    """
    Build %freq block with temperature configuration.
    Always returns %freq block when called, using temperature from config.
    """
    temperature = config.get('temperature', '298.15')
    return f"%freq\n  Temp {temperature}\nend\n"

def _build_bang_line(config, rel_token, main_basis, aux_jk, implicit,
                     include_freq=False, geom_key="geom_opt", qmmm_method: Optional[str] = None):
    """
    Construct the ORCA '!' line according to new CONTROL keys.
    include_freq=True adds the 'FREQ' keyword.
    geom_key selects which geometry token to use (e.g., 'geom_opt' or 'geom_opt_OCCUPIER').
    """
    ri_jkx = str(config.get("ri_jkx", "")).strip()
    disp   = str(config.get("disp_corr", "")).strip()
    geom   = str(config.get(geom_key, "")).strip()
    initg = (str(config.get("initial_guess", "")).split() or [""])[0]

    tokens = ["!"]
    if qmmm_method:
        tokens.append(qmmm_method)
    tokens.append(str(config["functional"]).strip())
    if rel_token:
        tokens.append(rel_token)
    tokens.append(str(main_basis).strip())
    if disp:
        tokens.append(disp)
    if ri_jkx:
        tokens.append(ri_jkx)
    if aux_jk:
        tokens.append(aux_jk)
    if implicit:
        tokens.append(implicit)
    if geom:
        tokens.append(geom)
    if include_freq:
        tokens.append("FREQ")
    tokens.append(initg)

    # normalize spacing
    return " ".join(t for t in tokens if t).replace("  ", " ").strip()

def _apply_per_atom_newgto(geom_lines: List[str], found_metals: List[str],
                           metal_basisset: Optional[str], config, radii_map):
    """
    Append per-atom 'NewGTO "metal_basisset" end' to
    - all metal atoms (always when metal_basisset provided),
    - atoms in first coordination sphere when enabled in CONTROL.
    """
    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes','true','1','on')

    if not metal_basisset:
        return geom_lines[:]  # nothing to do

    if not found_metals and not enable_first:
        return geom_lines[:]

    atoms = _parse_xyz_atoms(geom_lines)
    if not atoms:
        return geom_lines[:]

    # metal indices by symbol
    metal_syms = {m.strip().capitalize() for m in (found_metals or [])}
    metal_indices = [i for i, a in enumerate(atoms) if a["elem"].capitalize() in metal_syms]

    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()
    if sphere_scale_raw:
        scale = float(sphere_scale_raw)
    else:
        scale = 1.20
    first = _first_sphere_indices(atoms, metal_indices, scale, radii_map) if (enable_first and metal_indices) else set()

    metal_line_set = {atoms[i]['line_idx'] for i in metal_indices}
    first_line_set = {atoms[i]['line_idx'] for i in first}

    out = []
    for idx, line in enumerate(geom_lines):
        ls = line.strip()
        if not ls or ls == "*":
            out.append(line if line.endswith("\n") else line + "\n")
            continue
        if idx in metal_line_set or idx in first_line_set:
            line = line.rstrip() + f'   NewGTO "{metal_basisset}" end'
        out.append(line if line.endswith("\n") else line + "\n")
    return out

# -------------------------------------------------------------------------
# main writers (updated to new CONTROL keys + per-atom basis via utils)
# -------------------------------------------------------------------------

def read_and_modify_file(input_file_path, output_file_path, charge, multiplicity, solvent,
                         found_metals, metal_basisset, main_basisset, config, additions):
    """
    Build a generic ORCA input from an existing coordinate file (plain XYZ-like block).
    Applies: new '!' line (with ri_jkx/aux_jk/relativity via utils), optional print blocks,
    and per-atom NewGTO for metals (+ optional first sphere).
    """
    input_path = Path(input_file_path)
    with input_path.open('r') as file:
        coord_lines = [ln for ln in file.readlines() if ln.strip() and ln.strip() != "*"]

    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(coord_lines, input_path)
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = "QM/XTB" if qmmm_range else None

    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes','true','1','on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

    load_radii = enable_first and not sphere_scale_raw
    radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009")) if load_radii else None

    # decide main/metal bases per d3 vs. d4/5 policy; allow explicit overrides
    auto_main, auto_metal = set_main_basisset(found_metals, config)
    main  = main_basisset  or auto_main
    metal = metal_basisset or auto_metal

    # relativity & aux-JK selection (only active for 4d/5d per utils)
    rel_token, aux_jk, _ = select_rel_and_aux(found_metals, config)
    implicit = _implicit_token(config, solvent)

    # include FREQ only if frequency_calculation_OCCUPIER=yes
    include_freq = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes'
    bang = _build_bang_line(
        config,
        rel_token,
        main,
        aux_jk,
        implicit,
        include_freq=include_freq,
        geom_key="geom_opt",
        qmmm_method=qmmm_token,
    )

    output_blocks = collect_output_blocks(config)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(config['maxcore'], config['PAL'], resolve_maxiter(config))
    builder.add_additions(additions)
    if include_freq:
        builder.add_block(_build_freq_block(config))
    builder.add_blocks(output_blocks)

    lines = builder.lines

    # geometry
    lines.extend(build_qmmm_block(qmmm_range))
    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = _apply_per_atom_newgto(geom_lines, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    with open(output_file_path, 'w') as file:
        file.writelines(lines)

def read_and_modify_file_1(input_file_path, output_file_path, charge, multiplicity, solvent,
                         found_metals, metal_basisset, main_basisset, config, additions):
    """
    Build a generic ORCA input from an existing coordinate file (plain XYZ-like block).
    Applies: new '!' line (with ri_jkx/aux_jk/relativity via utils), optional print blocks,
    and per-atom NewGTO for metals (+ optional first sphere).

    NOTE: FREQ is always included on the '!' line.
    """
    input_path = Path(input_file_path)
    with input_path.open('r') as file:
        coord_lines = [ln for ln in file.readlines() if ln.strip() and ln.strip() != "*"]

    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(coord_lines, input_path)
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = "QM/XTB" if qmmm_range else None

    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes','true','1','on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

    load_radii = enable_first and not sphere_scale_raw
    radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009")) if load_radii else None

    # decide main/metal bases per d3 vs. d4/5 policy; allow explicit overrides
    auto_main, auto_metal = set_main_basisset(found_metals, config)
    main  = main_basisset  or auto_main
    metal = metal_basisset or auto_metal

    # relativity & aux-JK selection (only active for 4d/5d per utils)
    rel_token, aux_jk, _ = select_rel_and_aux(found_metals, config)
    implicit = _implicit_token(config, solvent)

    # ALWAYS include FREQ
    include_freq = True
    bang = _build_bang_line(
        config,
        rel_token,
        main,
        aux_jk,
        implicit,
        include_freq=include_freq,
        geom_key="geom_opt",
        qmmm_method=qmmm_token,
    )

    # Fallback guard: ensure 'FREQ' really present (in case _build_bang_line ignores the flag)
    if "FREQ" not in bang.upper():
        if bang.endswith("\n"):
            bang = bang.rstrip("\n") + " FREQ\n"
        else:
            bang = bang + " FREQ"

    output_blocks = collect_output_blocks(config)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(config['maxcore'], config['PAL'], resolve_maxiter(config))
    builder.add_additions(additions)
    if include_freq:
        builder.add_block(_build_freq_block(config))
    builder.add_blocks(output_blocks)

    lines = builder.lines

    # geometry
    lines.extend(build_qmmm_block(qmmm_range))
    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = _apply_per_atom_newgto(geom_lines, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    with open(output_file_path, 'w') as file:
        file.writelines(lines)


def read_xyz_and_create_input2(xyz_file_path: str, output_file_path: str, charge: int, multiplicity: int, solvent: str,
                               found_metals: List[str], config: Dict[str, Any], main_basisset: str, metal_basisset: Optional[str], additions: str) -> None:
    """
    TDDFT single-point builder (no freq). Uses new CONTROL keys and per-atom basis tagging.
    """
    xyz_path = Path(xyz_file_path)
    try:
        with xyz_path.open('r') as file:
            xyz_lines = file.readlines()[2:]  # skip natoms + comment
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(xyz_lines, xyz_path)
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = "QM/XTB" if qmmm_range else None

    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes','true','1','on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

    load_radii = enable_first and not sphere_scale_raw
    radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009")) if load_radii else None

    # TDDFT block
    triplet_flag = str(config.get("triplet_flag", "FALSE")).upper()
    tddft_block = [
        f"%TDDFT  NROOTS  {config['NROOTS']}\n",
        f"        Triplets   {triplet_flag}\n",
        f"        DOSOC     {config['DOSOC']}\n",
        f"        TDA       {config['TDA']}\n",
        f"        DONTO       {config['DONTO']}\n",
        "END\n"
    ]

    # bases (auto per utils, with optional explicit overrides)
    auto_main, auto_metal = set_main_basisset(found_metals, config)
    main  = main_basisset  or auto_main
    metal = metal_basisset or auto_metal

    # relativity token + aux-JK
    rel_token, aux_jk, _ = select_rel_and_aux(found_metals, config)
    implicit = _implicit_token(config, solvent)

    # method line (no freq)
    bang = _build_bang_line(
        config,
        rel_token,
        main,
        aux_jk,
        implicit,
        include_freq=False,
        geom_key="",
        qmmm_method=qmmm_token,
    )

    input_lines: List[str] = []
    input_lines.append(bang + "\n")
    input_lines.extend(tddft_block)
    input_lines.append(f"%maxcore {config['maxcore']}\n%pal nprocs {config['PAL']} end\n")
    maxiter_val = resolve_maxiter(config)
    if maxiter_val is not None:
        input_lines.append(f"%scf maxiter {maxiter_val} end\n")
    if additions and additions.strip():
        input_lines.append(f"{additions.strip()}\n")

    # geometry
    input_lines.extend(build_qmmm_block(qmmm_range))
    input_lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = _apply_per_atom_newgto(geom_lines, found_metals, metal, config, radii_all)
    input_lines.extend(geom)
    input_lines.append("*\n")

    try:
        with open(output_file_path, 'w') as file:
            file.writelines(input_lines)
        logging.info(f"XYZ file '{xyz_file_path}' processed and saved as '{output_file_path}'")
    except Exception as e:
        logging.error(f"Error writing '{output_file_path}': {e}")

def read_xyz_and_create_input2_2(xyz_file_path, output_file_path, charge, multiplicity,
                                 solvent, found_metals, config, main_basisset, metal_basisset, additions):
    """
    TDDFT + OPT builder (no freq). Uses new CONTROL keys and per-atom basis tagging.
    """
    xyz_path = Path(xyz_file_path)
    try:
        with xyz_path.open('r') as file:
            xyz_lines = file.readlines()[2:]
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(xyz_lines, xyz_path)
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = "QM/XTB" if qmmm_range else None

    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes','true','1','on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

    load_radii = enable_first and not sphere_scale_raw
    radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009")) if load_radii else None

    # TDDFT block
    tddft_block = [
        f"%TDDFT  NROOTS  {config['NROOTS']}\n",
        f"        DOSOC   {config['DOSOC']}\n",
        f"        TDA     {config['TDA']}\n",
        "END\n"
    ]

    # bases
    auto_main, auto_metal = set_main_basisset(found_metals, config)
    main  = main_basisset  or auto_main
    metal = metal_basisset or auto_metal

    # relativity + aux-JK
    rel_token, aux_jk, _ = select_rel_and_aux(found_metals, config)
    implicit = _implicit_token(config, solvent)

    # method line with OPT (no freq)
    bang = _build_bang_line(
        config,
        rel_token,
        main,
        aux_jk,
        implicit,
        include_freq=False,
        geom_key="geom_opt",
        qmmm_method=qmmm_token,
    )

    input_lines: List[str] = []
    input_lines.append(bang + "\n")
    input_lines.extend(tddft_block)
    input_lines.append(f"%maxcore {config['maxcore']}\n%pal nprocs {config['PAL']} end\n")
    maxiter_val = resolve_maxiter(config)
    if maxiter_val is not None:
        input_lines.append(f"%scf maxiter {maxiter_val} end\n")
    if additions and additions.strip():
        input_lines.append(f"{additions.strip()}\n")

    # geometry
    input_lines.extend(build_qmmm_block(qmmm_range))
    input_lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = _apply_per_atom_newgto(geom_lines, found_metals, metal, config, radii_all)
    input_lines.extend(geom)
    input_lines.append("*\n")

    try:
        with open(output_file_path, 'w') as file:
            file.writelines(input_lines)
        logging.info(f"XYZ file '{xyz_file_path}' processed and saved as '{output_file_path}'")
    except Exception as e:
        logging.error(f"Error writing '{output_file_path}': {e}")

def read_xyz_and_create_input3(xyz_file_path: str, output_file_path: str, charge: int, multiplicity: int,
                               solvent: str, found_metals: List[str], metal_basisset: Optional[str], main_basisset: str, config: Dict[str, Any], additions: str) -> None:
    """
    Frequency job builder (adds FREQ). Uses new CONTROL keys and per-atom basis tagging.
    """
    xyz_path = Path(xyz_file_path)
    try:
        with xyz_path.open('r') as file:
            xyz_lines = file.readlines()[2:]
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(xyz_lines, xyz_path)
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = "QM/XTB" if qmmm_range else None

    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes','true','1','on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

    load_radii = enable_first and not sphere_scale_raw
    radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009")) if load_radii else None

    # bases
    auto_main, auto_metal = set_main_basisset(found_metals, config)
    main  = main_basisset  or auto_main
    metal = metal_basisset or auto_metal

    # relativity + aux-JK
    rel_token, aux_jk, _ = select_rel_and_aux(found_metals, config)
    implicit = _implicit_token(config, solvent)

    # method line with FREQ
    include_freq = True
    bang = _build_bang_line(
        config,
        rel_token,
        main,
        aux_jk,
        implicit,
        include_freq=include_freq,
        geom_key="geom_opt",
        qmmm_method=qmmm_token,
    )

    output_blocks = collect_output_blocks(config)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(config['maxcore'], config['PAL'], resolve_maxiter(config))
    builder.add_additions(additions)
    if include_freq:
        builder.add_block(_build_freq_block(config))
    builder.add_blocks(output_blocks)

    lines = builder.lines

    lines.extend(build_qmmm_block(qmmm_range))
    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = _apply_per_atom_newgto(geom_lines, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    with open(output_file_path, 'w') as file:
        file.writelines(lines)
    logging.info(f"XYZ file '{xyz_file_path}' processed and saved as '{output_file_path}'")

def read_xyz_and_create_input4(xyz_file_path: str, output_file_path: str, charge: int, multiplicity: int,
                               solvent: str, found_metals: List[str], metal_basisset: Optional[str], main_basisset: str, config: Dict[str, Any], additions: str) -> None:
    """
    E00 / selected-root TDDFT builder (with %TDDFT IROOT/FOLLOWIROOT).
    Uses new CONTROL keys and per-atom basis tagging.
    """
    xyz_path = Path(xyz_file_path)
    try:
        with xyz_path.open('r') as file:
            xyz_lines = file.readlines()[2:]
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(xyz_lines, xyz_path)
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = "QM/XTB" if qmmm_range else None

    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes','true','1','on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

    load_radii = enable_first and not sphere_scale_raw
    radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009")) if load_radii else None

    # TDDFT block (root following)
    tddft_block = [
        "%TDDFT\n",
        f" NROOTS {config['NROOTS']}\n",
        f" IROOT {config['IROOT']}\n",
        f" FOLLOWIROOT   {config['FOLLOWIROOT']}\n",
        "end\n"
    ]

    # bases
    auto_main, auto_metal = set_main_basisset(found_metals, config)
    main  = main_basisset  or auto_main
    metal = metal_basisset or auto_metal

    # relativity + aux-JK
    rel_token, aux_jk, _ = select_rel_and_aux(found_metals, config)
    implicit = _implicit_token(config, solvent)

    # method line (with FREQ to mirror previous behavior)
    bang = _build_bang_line(
        config,
        rel_token,
        main,
        aux_jk,
        implicit,
        include_freq=True,
        geom_key="geom_opt",
        qmmm_method=qmmm_token,
    )

    lines: List[str] = []
    lines.append(bang + "\n")
    lines.extend(tddft_block)
    # special mcore for E00 flows (kept from your original)
    lines.append(f"%maxcore {config['mcore_E00']}\n%pal nprocs {config['PAL']} end\n")
    maxiter_val = resolve_maxiter(config)
    if maxiter_val is not None:
        lines.append(f"%scf maxiter {maxiter_val} end\n")
    if additions and additions.strip():
        lines.append(f"{additions.strip()}\n")

    lines.extend(build_qmmm_block(qmmm_range))
    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = _apply_per_atom_newgto(geom_lines, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    with open(output_file_path, 'w') as file:
        file.writelines(lines)
    logging.info(f"XYZ file '{xyz_file_path}' processed and saved as '{output_file_path}'")


def _create_s1_deltascf_input(xyz_file_path: str, output_file_path: str, charge: int, multiplicity: int,
                               solvent: str, found_metals: List[str], metal_basisset: Optional[str], main_basisset: str,
                               config: Dict[str, Any], additions: str) -> None:
    xyz_path = Path(xyz_file_path)
    try:
        with xyz_path.open('r') as file:
            xyz_lines = file.readlines()[2:]
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(xyz_lines, xyz_path)
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = "QM/XTB" if qmmm_range else None

    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes', 'true', '1', 'on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()

    load_radii = enable_first and not sphere_scale_raw
    radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009")) if load_radii else None

    auto_main, auto_metal = set_main_basisset(found_metals, config)
    main = main_basisset or auto_main
    metal = metal_basisset or auto_metal

    rel_token, aux_jk, _ = select_rel_and_aux(found_metals, config)
    implicit = _implicit_token(config, solvent)

    functional = str(config.get("functional", "PBE0")).strip()
    disp = str(config.get("disp_corr", "")).strip()
    ri_jkx = str(config.get("ri_jkx", "")).strip()
    geom_token = str(config.get("geom_opt", "")).strip()
    init_guess = (str(config.get("initial_guess", "")).split() or [""])[0]

    tokens = ["!"]
    if qmmm_token:
        tokens.append(qmmm_token)
    tokens.extend([functional, "UKS"])
    if rel_token:
        tokens.append(rel_token)
    tokens.append(str(main).strip())
    if disp:
        tokens.append(disp)
    if ri_jkx:
        tokens.append(ri_jkx)
    if aux_jk:
        tokens.append(aux_jk)
    if implicit:
        tokens.append(implicit)
    tokens.append("deltaSCF")
    if geom_token:
        tokens.append(geom_token)
    tokens.append("FREQ")
    if init_guess:
        tokens.append(init_guess)

    bang = " ".join(token for token in tokens if token)

    builder = OrcaInputBuilder(bang)
    lines = builder.lines

    maxiter_val = resolve_maxiter(config)
    scf_lines = ["%scf\n"]
    if maxiter_val is not None:
        scf_lines.append(f" maxiter {maxiter_val}\n")

    def _maybe_value(key: str, default: str) -> Optional[str]:
        raw = config.get(key, default)
        if raw is None:
            return None
        text = str(raw).strip()
        return text or None

    def _maybe_bool(key: str, default: str) -> Optional[str]:
        value = _maybe_value(key, default)
        if value is None:
            return None
        lowered = value.lower()
        if lowered in {"true", "false"}:
            return lowered
        return value

    scf_params = (
        ("DOMOM", _maybe_bool('deltaSCF_DOMOM', 'true')),
        ("pmom", _maybe_bool('deltaSCF_PMOM', 'true')),
        ("keepinitialref", _maybe_bool('deltaSCF_keepinitialref', 'true')),
        ("alphaconf", "0,1"),
        ("betaconf", "0"),
        ("SOSCFHESSUP", _maybe_value('deltaSCF_SOSCFHESSUP', 'LBFGS')),
    )

    for key, value in scf_params:
        if value is not None:
            scf_lines.append(f" {key} {value}\n")

    scf_lines.append("end\n")
    lines.append("".join(scf_lines))

    builder.add_resources(config['maxcore'], config['PAL'], None)
    builder.add_additions(additions)
    builder.add_blocks(collect_output_blocks(config))

    lines.extend(build_qmmm_block(qmmm_range))
    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = _apply_per_atom_newgto(geom_lines, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    try:
        with open(output_file_path, 'w') as file:
            file.writelines(lines)
        logging.info(f"XYZ file '{xyz_file_path}' processed and saved as '{output_file_path}'")
    except Exception as exc:  # noqa: BLE001
        logging.error(f"Error writing '{output_file_path}': {exc}")


def create_s1_optimization_input(xyz_file_path: str, output_file_path: str, charge: int, multiplicity: int,
                                 solvent: str, found_metals: List[str], metal_basisset: Optional[str],
                                 main_basisset: str, config: Dict[str, Any], additions: str) -> None:
    method_value = config.get('s1_opt', config.get('S1_opt', 'TDDFT'))
    method_raw = str(method_value).strip().lower()
    if method_raw in {'deltascf', 'delta scf', 'dscf'}:
        _create_s1_deltascf_input(
            xyz_file_path,
            output_file_path,
            charge,
            multiplicity,
            solvent,
            found_metals,
            metal_basisset,
            main_basisset,
            config,
            additions,
        )
        return

    # Fallback to existing TDDFT workflow
    read_xyz_and_create_input4(
        xyz_file_path,
        output_file_path,
        charge,
        multiplicity,
        solvent,
        found_metals,
        metal_basisset,
        main_basisset,
        config,
        additions,
    )
