import json
import logging
import math
import re
import threading
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

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
                with cache_file.open("w", encoding="utf-8") as fh:
                    json.dump(data, fh)
        except Exception:
            continue


def normalize_xyz_body(lines: Iterable[str]) -> List[str]:
    """Normalize XYZ coordinate lines by stripping blanks and trailing markers."""
    normalized: List[str] = []
    for raw in lines:
        if raw is None:
            continue
        stripped = str(raw).strip()
        if not stripped or stripped == "*":
            continue
        normalized.append(stripped + ("\n" if not stripped.endswith("\n") else ""))
    return normalized


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
    normalized = list(normalize_xyz_body(coord_lines))

    qm_lines: List[str] = []
    mm_lines: List[str] = []
    separator_seen = False

    for raw in normalized:
        stripped = raw.strip()
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
        f"  QMAtoms {{{start}:{end}}} end\n",
        "  Charge_Medium 0\n",
        "  Mult_Medium   1\n",
        "end\n",
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
    If CPCM is requested while QM/XTB is active, automatically switch to ALPB
    and emit a warning; other models are left untouched.
    For other QM/MM methods (QM/PBEH-3C, etc.), CPCM is kept as-is.
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

    # Only switch to ALPB for QM/XTB (xTB requires ALPB for implicit solvation)
    # Other methods like QM/PBEH-3C support CPCM directly
    qmmm_method = str(config.get('qmmm_option', 'QM/XTB')).strip().upper()
    if model.upper() == "CPCM" and qmmm_method == "QM/XTB":
        logging.warning('Detected "$" separator with QM/XTB: CPCM is incompatible with xTB. Switching implicit_solvation_model to ALPB.')
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
    include_freq=True adds the frequency keyword (FREQ or NUMFREQ based on freq_type).
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
        freq_type = str(config.get("freq_type", "FREQ")).strip().upper()
        if freq_type not in ["FREQ", "NUMFREQ"]:
            freq_type = "FREQ"
        tokens.append(freq_type)
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
    qmmm_token = str(config.get('qmmm_option', 'QM/XTB')).strip() if qmmm_range else None

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

    output_blocks = collect_output_blocks(config, allow=True)
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
    qmmm_token = str(config.get('qmmm_option', 'QM/XTB')).strip() if qmmm_range else None

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

    # Fallback guard: ensure freq keyword really present (in case _build_bang_line ignores the flag)
    freq_type = str(config.get("freq_type", "FREQ")).strip().upper()
    if freq_type not in ["FREQ", "NUMFREQ"]:
        freq_type = "FREQ"
    if freq_type not in bang.upper():
        if bang.endswith("\n"):
            bang = bang.rstrip("\n") + f" {freq_type}\n"
        else:
            bang = bang + f" {freq_type}"

    output_blocks = collect_output_blocks(config, allow=True)
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
    qmmm_token = str(config.get('qmmm_option', 'QM/XTB')).strip() if qmmm_range else None

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

    output_blocks = collect_output_blocks(config, allow=True)
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
