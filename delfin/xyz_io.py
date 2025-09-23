import logging, re, math
from typing import List, Dict, Optional, Any

from delfin.common.orca_blocks import OrcaInputBuilder, collect_output_blocks

# import canonical selection helpers from utils
from .utils import set_main_basisset, select_rel_and_aux

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
            radii[el.symbol] = float(r)
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

def _build_freq_block(config):
    """
    Build %freq block with temperature configuration.
    Always returns %freq block when called, using temperature from config.
    """
    temperature = config.get('temperature', '298.15')
    return f"%freq\n  Temp {temperature}\nend\n"

def _build_bang_line(config, rel_token, main_basis, aux_jk, implicit, include_freq=False, geom_key="geom_opt"):
    """
    Construct the ORCA '!' line according to new CONTROL keys.
    include_freq=True adds the 'FREQ' keyword.
    geom_key selects which geometry token to use (e.g., 'geom_opt' or 'geom_opt_OCCUPIER').
    """
    ri_jkx = str(config.get("ri_jkx", "")).strip()
    disp   = str(config.get("disp_corr", "")).strip()
    geom   = str(config.get(geom_key, "")).strip()
    initg = (str(config.get("initial_guess", "")).split() or [""])[0]

    tokens = ["!", str(config["functional"]).strip()]
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
    with open(input_file_path, 'r') as file:
        coord_lines = [ln for ln in file.readlines() if ln.strip() and ln.strip() != "*"]

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
    bang = _build_bang_line(config, rel_token, main, aux_jk, implicit,
                            include_freq=include_freq, geom_key="geom_opt")

    output_blocks = collect_output_blocks(config)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(config['maxcore'], config['PAL'])
    builder.add_additions(additions)
    if include_freq:
        builder.add_block(_build_freq_block(config))
    builder.add_blocks(output_blocks)

    lines = builder.lines

    # geometry
    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = [ln if ln.endswith("\n") else ln + "\n" for ln in coord_lines]
    geom = _apply_per_atom_newgto(geom, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    with open(output_file_path, 'w') as file:
        file.writelines(lines)
    logging.info(f"Input file '{input_file_path}' modified and saved as '{output_file_path}'")

def read_and_modify_file_1(input_file_path, output_file_path, charge, multiplicity, solvent,
                         found_metals, metal_basisset, main_basisset, config, additions):
    """
    Build a generic ORCA input from an existing coordinate file (plain XYZ-like block).
    Applies: new '!' line (with ri_jkx/aux_jk/relativity via utils), optional print blocks,
    and per-atom NewGTO for metals (+ optional first sphere).

    NOTE: FREQ is always included on the '!' line.
    """
    with open(input_file_path, 'r') as file:
        coord_lines = [ln for ln in file.readlines() if ln.strip() and ln.strip() != "*"]

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
        config, rel_token, main, aux_jk, implicit,
        include_freq=include_freq, geom_key="geom_opt"
    )

    # Fallback guard: ensure 'FREQ' really present (in case _build_bang_line ignores the flag)
    if "FREQ" not in bang.upper():
        if bang.endswith("\n"):
            bang = bang.rstrip("\n") + " FREQ\n"
        else:
            bang = bang + " FREQ"

    output_blocks = collect_output_blocks(config)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(config['maxcore'], config['PAL'])
    builder.add_additions(additions)
    if include_freq:
        builder.add_block(_build_freq_block(config))
    builder.add_blocks(output_blocks)

    lines = builder.lines

    # geometry
    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = [ln if ln.endswith("\n") else ln + "\n" for ln in coord_lines]
    geom = _apply_per_atom_newgto(geom, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    with open(output_file_path, 'w') as file:
        file.writelines(lines)
    logging.info(f"Input file '{input_file_path}' modified and saved as '{output_file_path}'")


def read_xyz_and_create_input2(xyz_file_path: str, output_file_path: str, charge: int, multiplicity: int, solvent: str,
                               found_metals: List[str], config: Dict[str, Any], main_basisset: str, metal_basisset: Optional[str], additions: str) -> None:
    """
    TDDFT single-point builder (no freq). Uses new CONTROL keys and per-atom basis tagging.
    """
    try:
        with open(xyz_file_path, 'r') as file:
            xyz_lines = file.readlines()[2:]  # skip natoms + comment
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

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
    bang = _build_bang_line(config, rel_token, main, aux_jk, implicit,
                            include_freq=False, geom_key="")

    input_lines: List[str] = []
    input_lines.append(bang + "\n")
    input_lines.extend(tddft_block)
    input_lines.append(f"%maxcore {config['maxcore']}\n%pal nprocs {config['PAL']} end\n")
    if additions and additions.strip():
        input_lines.append(f"{additions.strip()}\n")

    # geometry
    input_lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = [ln if ln.endswith("\n") else ln + "\n" for ln in xyz_lines]
    geom = _apply_per_atom_newgto(geom, found_metals, metal, config, radii_all)
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
    try:
        with open(xyz_file_path, 'r') as file:
            xyz_lines = file.readlines()[2:]
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

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
    bang = _build_bang_line(config, rel_token, main, aux_jk, implicit,
                            include_freq=False, geom_key="geom_opt")

    input_lines: List[str] = []
    input_lines.append(bang + "\n")
    input_lines.extend(tddft_block)
    input_lines.append(f"%maxcore {config['maxcore']}\n%pal nprocs {config['PAL']} end\n")
    if additions and additions.strip():
        input_lines.append(f"{additions.strip()}\n")

    # geometry
    input_lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = [ln if ln.endswith("\n") else ln + "\n" for ln in xyz_lines]
    geom = _apply_per_atom_newgto(geom, found_metals, metal, config, radii_all)
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
    try:
        with open(xyz_file_path, 'r') as file:
            xyz_lines = file.readlines()[2:]
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

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
    bang = _build_bang_line(config, rel_token, main, aux_jk, implicit,
                            include_freq=include_freq, geom_key="geom_opt")

    output_blocks = collect_output_blocks(config)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(config['maxcore'], config['PAL'])
    builder.add_additions(additions)
    if include_freq:
        builder.add_block(_build_freq_block(config))
    builder.add_blocks(output_blocks)

    lines = builder.lines

    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = [ln if ln.endswith("\n") else ln + "\n" for ln in xyz_lines]
    geom = _apply_per_atom_newgto(geom, found_metals, metal, config, radii_all)
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
    try:
        with open(xyz_file_path, 'r') as file:
            xyz_lines = file.readlines()[2:]
    except FileNotFoundError:
        logging.error(f"File not found: {xyz_file_path}")
        return

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
    bang = _build_bang_line(config, rel_token, main, aux_jk, implicit,
                            include_freq=True, geom_key="geom_opt")

    lines: List[str] = []
    lines.append(bang + "\n")
    lines.extend(tddft_block)
    # special mcore for E00 flows (kept from your original)
    lines.append(f"%maxcore {config['mcore_E00']}\n%pal nprocs {config['PAL']} end\n")
    if additions and additions.strip():
        lines.append(f"{additions.strip()}\n")

    lines.append(f"* xyz {charge} {multiplicity}\n")
    geom = [ln if ln.endswith("\n") else ln + "\n" for ln in xyz_lines]
    geom = _apply_per_atom_newgto(geom, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    with open(output_file_path, 'w') as file:
        file.writelines(lines)
    logging.info(f"XYZ file '{xyz_file_path}' processed and saved as '{output_file_path}'")
