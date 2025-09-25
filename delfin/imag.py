import os, re, sys, shutil, logging, subprocess, math
from pathlib import Path

from delfin.common.orca_blocks import OrcaInputBuilder, collect_output_blocks

from .utils import search_transition_metals, set_main_basisset, select_rel_and_aux
from .orca import run_orca_IMAG

OK_MARKER = "ORCA TERMINATED NORMALLY"

# ------------------------ small helpers ------------------------

# compact fallback radii (Å); good enough for 1st-sphere detection
_COVALENT_RADII_FALLBACK = {
    "H": 0.31, "He": 0.28,
    "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.25, "Co": 1.26, "Ni": 1.21, "Cu": 1.38, "Zn": 1.31,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
}

def _elem_from_label(label: str) -> str:
    m = re.match(r"([A-Za-z]{1,2})", label.strip())
    return m.group(1) if m else label.strip()

def _dist(a, b):
    return math.sqrt((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2)

def _parse_xyz_atoms(xyz_lines):
    """Parse atom lines (stop at '*' or blank). Returns list with coords and original line index."""
    atoms = []
    for idx, line in enumerate(xyz_lines):
        ls = line.strip()
        if not ls or ls == '*':
            break
        parts = ls.split()
        if len(parts) < 4:
            continue
        elem = _elem_from_label(parts[0])
        try:
            x, y, z = map(float, parts[1:4])
        except ValueError:
            continue
        atoms.append({"line_idx": idx, "elem": elem, "x": x, "y": y, "z": z})
    return atoms

def _rcov(sym: str):
    return float(_COVALENT_RADII_FALLBACK.get(sym, 1.20))

def _first_sphere_indices(atoms, metal_indices, scale):
    """Return indices of atoms in the first sphere of any metal using covalent radii rule."""
    first = set()
    for im in metal_indices:
        m = atoms[im]
        r_m = _rcov(m["elem"])
        for i, a in enumerate(atoms):
            if i == im:
                continue
            r_a = _rcov(a["elem"])
            cutoff = scale * (r_m + r_a)
            if _dist(m, a) <= cutoff:
                first.add(i)
    return first

def _build_bang_line_IMAG(config, rel_token, main_basisset, aux_jk, implicit):
    """
    Construct the '!' line for IMAG iterations:
      functional [REL] main_basis [disp] [ri_jkx] [aux_jk] [implicit] [geom_opt] FREQ initial_guess
    """
    ri_jkx = str(config.get("ri_jkx", "")).strip()
    disp   = str(config.get("disp_corr", "")).strip()
    geom   = str(config.get("geom_opt", "OPT")).strip()
    init_tokens = str(config.get("initial_guess", "PModel")).split()
    initg  = init_tokens[0] if init_tokens else "PModel"

    tokens = ["!", str(config["functional"])]
    if rel_token:
        tokens.append(rel_token)           # ZORA / X2C / DKH or ''
    tokens.append(str(main_basisset))
    if disp:
        tokens.append(disp)
    if ri_jkx:
        tokens.append(ri_jkx)
    if aux_jk:
        tokens.append(aux_jk)              # def2/J or SARC/J
    if implicit:
        tokens.append(implicit)
    if geom:
        tokens.append(geom)
    tokens.append("FREQ")
    tokens.append(initg)
    return " ".join(t for t in tokens if t).replace("  ", " ").strip()

def _has_ok_marker(path: str | Path) -> bool:
    candidate = Path(path)
    if not candidate.exists():
        return False
    try:
        with candidate.open("r", errors="ignore") as f:
            return OK_MARKER in f.read()
    except Exception:
        return False

def search_imaginary_mode2(log_file):
    """Return the (most negative) imaginary freq value in cm**-1, or None if none is present."""
    try:
        with open(log_file, 'r', errors="ignore") as file:
            for line in file:
                if "***imaginary mode***" in line:
                    m = re.search(r'(-?\d+(?:\.\d+)?)\s+cm\*\*-1', line)
                    if m:
                        freq_value = float(m.group(1))
                        logging.info(f"Imaginary mode found: {freq_value} cm**-1 in {log_file}")
                        return freq_value
        logging.info(f"No imaginary mode found in {log_file}.")
        return None
    except FileNotFoundError:
        logging.error(f"Log file '{log_file}' not found.")
        sys.exit(1)

def _imag_resolved(out_path: str | Path, threshold: float) -> bool:
    candidate = Path(out_path)
    if not candidate.exists():
        return False
    freq = search_imaginary_mode2(str(candidate))
    return (freq is None) or (freq >= threshold)

def _find_last_ok_iteration(folder: str | Path):
    best_i, best_path = None, None
    folder_path = Path(folder)
    for entry in folder_path.iterdir():
        if not entry.is_file():
            continue
        m = re.fullmatch(r"output_(\d+)\.out", entry.name)
        if not m:
            continue
        i = int(m.group(1))
        if _has_ok_marker(entry) and (best_i is None or i > best_i):
            best_i, best_path = i, entry
    return best_i, best_path

def run_plotvib(iteration):
    try:
        plotvib_cmd = f"orca_pltvib input_{iteration}.hess 6"
        subprocess.run(plotvib_cmd, shell=True, check=True)
        logging.info(f"plotvib run successful for 'input_{iteration}.hess'")
        new_structure_file = Path(f"input_{iteration}.hess.v006.xyz")
        if not new_structure_file.exists():
            logging.error(f"Expected structure file '{new_structure_file}' not found after plotvib.")
            sys.exit(1)
        return str(new_structure_file)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running plotvib: {e}")
        sys.exit(1)

def run_plotvib0(input_file):
    try:
        plotvib_cmd = f"orca_pltvib {input_file}.hess 6"
        subprocess.run(plotvib_cmd, shell=True, check=True)
        logging.info(f"plotvib run successful for '{input_file}.hess'")
        new_structure_file = Path(f"{input_file}.hess.v006.xyz")
        if not new_structure_file.exists():
            logging.error(f"Expected structure file '{new_structure_file}' not found after plotvib.")
            sys.exit(1)
        return str(new_structure_file)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running plotvib: {e}")
        sys.exit(1)

# ----------------------- main writer (updated) -----------------------

def read_and_modify_xyz_IMAG(input_file_path, output_file_path, charge, multiplicity,
                             solvent, metals, config, main_basisset, metal_basisset, additions):
    """
    Build IMAG iteration input:
      - policy-based "!" line with ri_jkx / aux_jk and ZORA/X2C/DKH enabled for 4d/5d metals
      - main/metal basis sets via ``utils.set_main_basisset`` (arguments act as overrides)
      - per-atom ``NewGTO`` tags for metals plus optional first coordination sphere entries
      - FREQ job configuration and optional %output blocks just like the original workflow
    """
    try:
        with open(input_file_path, 'r') as file:
            xyz_lines = file.readlines()
        coord_block = [ln for ln in xyz_lines if ln.strip()]  # keep raw non-empty lines
    except Exception as e:
        logging.error(f"Error reading '{input_file_path}': {e}")
        sys.exit(1)

    # Determine metals from the current structure
    found_metals_local = search_transition_metals(input_file_path)

    # Select base rates according to 3d/4d5 policy (arguments serve only as fallback)
    main_sel, metal_sel = set_main_basisset(found_metals_local, config)
    main_eff  = main_basisset  or main_sel
    metal_eff = metal_basisset or metal_sel

    # Relativity/AUX policy (3d → non-rel + aux_jk; 4d/5d → rel + aux_jk_rel)
    rel_token, aux_jk_token, _ = select_rel_and_aux(found_metals_local, config)

    # implicit solvation
    implicit = ""
    model = str(config.get('implicit_solvation_model', '')).strip()
    if model:
        implicit = f"{model}({solvent})" if solvent else model

    # '!' line for IMAG (always FREQ)
    bang = _build_bang_line_IMAG(config, rel_token, main_eff, aux_jk_token, implicit)
    if additions and "moinp" in additions.lower():
        bang += " MORead"

    # 1. Coordination sphere
    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes', 'true', '1', 'on')
    sphere_scale_raw = str(config.get('first_coordination_sphere_scale', '')).strip()
    scale = float(sphere_scale_raw) if sphere_scale_raw else 1.20

    # only output the first contiguous coordinate block
    coord_only = []
    for ln in coord_block:
        ls = ln.strip()
        if not ls or ls == '*':
            break
        coord_only.append(ln)

    atoms = _parse_xyz_atoms(coord_only)
    metal_syms = {m.strip().capitalize() for m in (found_metals_local or [])}
    metal_indices = [i for i, a in enumerate(atoms) if a["elem"].capitalize() in metal_syms]
    first_indices = _first_sphere_indices(atoms, metal_indices, scale) if (enable_first and metal_indices and metal_eff) else set()
    metal_line_set = {atoms[i]['line_idx'] for i in metal_indices}
    first_line_set = {atoms[i]['line_idx'] for i in first_indices}

    # Assemble input
    output_blocks = collect_output_blocks(config)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(config['maxcore'], config['PAL'])
    builder.add_additions(additions)

    # Add %freq block with temperature (IMAG always uses FREQ)
    from .xyz_io import _build_freq_block
    freq_block = _build_freq_block(config)
    builder.add_block(freq_block)
    builder.add_blocks(output_blocks)

    out = builder.lines
    out.append(f"* xyz {charge} {multiplicity}\n")

    # Write coordinates + optional NewGTO tags
    for idx, ln in enumerate(coord_only):
        parts = ln.split()
        if len(parts) < 4:
            continue
        elem, x, y, z = parts[0], parts[1], parts[2], parts[3]
        line = f"{elem} {x} {y} {z}"
        apply_basis = metal_eff and (idx in metal_line_set or idx in first_line_set)
        if apply_basis:
            line += f'   NewGTO "{metal_eff}" end'
        out.append(line + "\n")

    out.append("*\n")

    try:
        with open(output_file_path, 'w') as f:
            f.writelines(out)
        logging.info(f"Input file '{output_file_path}' created successfully.")
    except Exception as e:
        logging.error(f"Error writing '{output_file_path}': {e}")
        sys.exit(1)

# ----------------------- rest of the logic -----------------------

def extract_structure(input_file, iteration):
    try:
        with open(input_file, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        logging.info(f"Process file: {input_file}")
        star_indices = [i for i, line in enumerate(lines) if '*' in line]
        logging.info(f"Number of found '*' lines: {len(star_indices)}")
        if len(star_indices) != 20:
            raise ValueError(f"File '{input_file}' does not contain exactly 20 '*' lines.")
        start_index = star_indices[4] - 1
        end_index = star_indices[5]
        extracted_lines = lines[start_index:end_index]
        output_file = f"input_{iteration}_structure_5.xyz"
        with open(output_file, 'w', encoding='utf-8') as file:
            file.writelines(extracted_lines)
        logging.info(f"Structure extracted to '{output_file}'")
    except Exception as e:
        logging.error(f"Error extracting structure: {e}")
        sys.exit(1)

def extract_structure0(input_file):
    input_file2 = f"{input_file}.hess.v006.xyz"
    try:
        with open(input_file2, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        logging.info(f"Process file: {input_file2}")
        star_indices = [i for i, line in enumerate(lines) if '*' in line]
        logging.info(f"Number of found '*' lines: {len(star_indices)}")
        if len(star_indices) != 20:
            raise ValueError(f"File '{input_file2}' does not contain exactly 20 '*' lines.")
        start_index = star_indices[4] - 1
        end_index = star_indices[5]
        extracted_lines = lines[start_index:end_index]
        output_file = f"{input_file}_structure_5.xyz"
        with open(output_file, 'w', encoding='utf-8') as file:
            file.writelines(extracted_lines)
        logging.info(f"Structure extracted to '{output_file}'")
    except Exception as e:
        logging.error(f"Error extracting structure: {e}")
        sys.exit(1)

def run_IMAG(input_file, hess_file, charge, multiplicity, solvent, metals, config,
             main_basisset, metal_basisset, additions):
    if str(config.get("IMAG", "no")).lower() != "yes":
        return
    print("""
                      *******************
                      *       IMAG      *
                      *******************
    """)
    allow_raw = float(config.get('allow_imaginary_freq', 0))
    threshold = allow_raw if allow_raw <= 0 else -allow_raw
    recalc = str(os.environ.get("DELFIN_RECALC", "0")).lower() in ("1", "true", "yes", "on")

    first_freq = search_imaginary_mode2(input_file)
    if first_freq is None or first_freq >= threshold:
        logging.info(f"Imag within tolerance (freq={first_freq}, thr={threshold}) -> no IMAG.")
        return

    imag_folder = Path(f"{hess_file}_IMAG")
    if recalc and imag_folder.is_dir():
        last_i, last_out = _find_last_ok_iteration(imag_folder)
        if last_i is not None and last_out and _imag_resolved(last_out, threshold):
            logging.info(f"[recalc] IMAG complete at iter {last_i}; skip.")
            final_log = Path(last_out)
            final_xyz = imag_folder / f"input_{last_i}.xyz"
            parent = imag_folder.resolve().parent
            input_path = Path(input_file)
            dest_log = input_path if input_path.is_absolute() else parent / input_path
            dest_xyz = parent / f"{hess_file}.xyz"
            try: shutil.copy2(final_log, dest_log)
            except Exception as e: logging.warning(f"copy log: {e}")
            if final_xyz.exists():
                try: shutil.copy2(final_xyz, dest_xyz)
                except Exception as e: logging.warning(f"copy xyz: {e}")
            return
        else:
            try: shutil.rmtree(imag_folder)
            except Exception: pass

    run_plotvib0(hess_file)
    extract_structure0(hess_file)
    seed_xyz = f"{hess_file}_structure_5.xyz"
    imag_folder.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copy2(seed_xyz, imag_folder / seed_xyz)
        logging.info(f"Seed structure copied to {imag_folder / seed_xyz}")
    except Exception as e:
        logging.error(f"Seed copy failed: {e}")
        return

    cwd = Path.cwd()
    os.chdir(imag_folder)
    iteration = 1
    try:
        while True:
            current_input_file = seed_xyz if iteration == 1 else f"input_{iteration - 1}_structure_5.xyz"
            output_file = f"input_{iteration}.inp"
            if iteration == 1:
                candidates = [Path("..") / f"{hess_file}.gbw", Path("..") / "input.gbw"]
                gbw_path = next((c for c in candidates if c.exists()), None)
            else:
                gbw_prev = f"input_{iteration - 1}.gbw"
                gbw_path = Path(gbw_prev) if Path(gbw_prev).exists() else None
            parts = []
            if gbw_path: parts.append(f'%moinp "{gbw_path}"')
            if additions and additions.strip(): parts.append(additions.strip())
            additions_eff = "\n".join(parts)

            read_and_modify_xyz_IMAG(current_input_file, output_file, charge, multiplicity,
                                     solvent, metals, config, main_basisset, metal_basisset, additions_eff)
            run_orca_IMAG(output_file, iteration)

            log_file = f"output_{iteration}.out"
            if not _has_ok_marker(log_file):
                logging.warning(f"Iteration {iteration}: ORCA not normal; stop.")
                break

            freq = search_imaginary_mode2(log_file)
            if (freq is None) or (freq >= threshold):
                logging.info(f"Iteration {iteration}: Imag resolved (freq={freq}, thr={threshold}).")
                break

            new_structure_file = run_plotvib(iteration)
            extract_structure(new_structure_file, iteration)
            iteration += 1
    finally:
        os.chdir(cwd)

    final_log_file = imag_folder / f"output_{iteration}.out"
    final_xyz_file = imag_folder / f"input_{iteration}.xyz"
    destination_folder = Path.cwd()
    input_path = Path(input_file)
    destination_log = input_path if input_path.is_absolute() else destination_folder / input_path
    destination_structure = destination_folder / f"{hess_file}.xyz"
    if final_log_file.exists():
        try: shutil.copy2(final_log_file, destination_log); print(f"Log file 'output_{iteration}.out' copied back as '{destination_log.name}'.")
        except Exception as e: logging.warning(f"copy back log: {e}")
    else:
        print(f"ERROR: Log file 'output_{iteration}.out' not found.")
    if final_xyz_file.exists():
        try: shutil.copy2(final_xyz_file, destination_structure); print(f"Structure file 'input_{iteration}.xyz' copied back as '{hess_file}.xyz'.")
        except Exception as e: logging.warning(f"copy back xyz: {e}")
    else:
        print(f"ERROR: Structure file 'input_{iteration}.xyz' not found.")
