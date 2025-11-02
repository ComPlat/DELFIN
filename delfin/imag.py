import os, re, sys, shutil, logging, subprocess, math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

from delfin.common.orca_blocks import OrcaInputBuilder, collect_output_blocks, resolve_maxiter

from .utils import search_transition_metals, set_main_basisset, select_rel_and_aux
from .energies import find_electronic_energy
from .orca import run_orca_IMAG, run_orca
from .xyz_io import (
    split_qmmm_sections,
    _ensure_qmmm_implicit_model,
    build_qmmm_block,
    _apply_per_atom_newgto,
    _load_covalent_radii,
)

OK_MARKER = "ORCA TERMINATED NORMALLY"
ENERGY_IMPROVEMENT_TOL = 1e-6

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

def _strip_xyz_header(lines):
    """Remove leading atom-count/comment lines from XYZ-like fragments."""
    if not lines:
        return lines
    working = list(lines)
    first = working[0].strip()
    try:
        int(first)
        working = working[1:]
        if working:
            head = working[0].strip().split()
            if (len(head) < 4) or (not head[0]) or (not head[0][0].isalpha()):
                working = working[1:]
        return working
    except ValueError:
        return lines

def _trim_xyz_columns(lines):
    """Keep only element + XYZ coordinates (preserve separators and QMMM markers)."""
    out = []
    for raw in lines:
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.isdigit():
            break
        if stripped in {"*", "$"}:
            out.append(stripped + "\n")
            continue
        parts = stripped.split()
        if len(parts) >= 4 and parts[0] and parts[0][0].isalpha():
            tail = []
            if "NewGTO" in parts:
                idx = parts.index("NewGTO")
                tail = parts[idx:]
            base = parts[:4]
            rebuilt = " ".join(base + tail)
            out.append(rebuilt + "\n")
        else:
            out.append(stripped + "\n")
    return out

_COORD_LINE_RE = re.compile(
    r"^(?P<lead>\s*)(?P<elem>\S+)(?P<sp1>\s+)(?P<x>\S+)(?P<sp2>\s+)"
    r"(?P<y>\S+)(?P<sp3>\s+)(?P<z>\S+)(?P<rest>.*)$"
)


def _load_inp_template(template_path: Path | str):
    path = Path(template_path)
    try:
        with path.open("r", encoding="utf-8") as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        logging.warning(f"IMAG template input '{path}' not found; falling back to generated input.")
        return None
    except Exception as exc:
        logging.warning(f"Failed to read IMAG template '{path}': {exc}; falling back to generated input.")
        return None

    geom_start = None
    for idx, line in enumerate(lines):
        if line.strip().lower().startswith("* xyz"):
            geom_start = idx + 1
            break
    if geom_start is None:
        logging.warning(f"Template '{path}' missing '* xyz' section; fallback to generated input.")
        return None

    geom_end = geom_start
    while geom_end < len(lines) and lines[geom_end].strip() != "*":
        geom_end += 1
    if geom_end >= len(lines):
        logging.warning(f"Template '{path}' missing terminating '*' for geometry; fallback to generated input.")
        return None

    coord_count = sum(
        1
        for idx in range(geom_start, geom_end)
        if _COORD_LINE_RE.match(lines[idx].rstrip("\n"))
    )
    return {
        "path": path,
        "lines": lines,
        "geom_start": geom_start,
        "geom_end": geom_end,
        "coord_count": coord_count,
    }


def _extract_resources_from_input(input_path: Path | str) -> tuple[int | None, int | None]:
    path = Path(input_path)
    if not path.exists():
        return None, None
    pal_val = None
    maxcore_val = None
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                stripped = line.strip().lower()
                if stripped.startswith("%pal") and "nprocs" in stripped:
                    parts = stripped.replace("=", " ").split()
                    for idx, token in enumerate(parts):
                        if token == "nprocs" and idx + 1 < len(parts):
                            try:
                                pal_val = int(parts[idx + 1])
                            except ValueError:
                                pal_val = None
                            break
                elif stripped.startswith("%maxcore"):
                    parts = stripped.split()
                    if len(parts) >= 2:
                        try:
                            maxcore_val = int(parts[1])
                        except ValueError:
                            maxcore_val = None
                if pal_val is not None and maxcore_val is not None:
                    break
    except Exception:
        pal_val = pal_val
        maxcore_val = maxcore_val
    return pal_val, maxcore_val


def _extract_xyz_coordinates(xyz_path: Path | str) -> list[tuple[str, str, str, str]]:
    path = Path(xyz_path)
    try:
        with path.open("r", encoding="utf-8") as fh:
            raw_lines = [ln for ln in fh.readlines() if ln.strip()]
    except Exception as exc:
        logging.error(f"Error reading geometry '{path}': {exc}")
        return []

    lines = _strip_xyz_header(raw_lines)
    coords: list[tuple[str, str, str, str]] = []
    for raw in lines:
        stripped = raw.strip()
        if not stripped or stripped == "*" or stripped.startswith("*"):
            break
        parts = stripped.split()
        if len(parts) < 4:
            continue
        coords.append((parts[0], parts[1], parts[2], parts[3]))
    return coords


def _sanitize_template_lines(lines: list[str]) -> list[str]:
    sanitized: list[str] = []
    for line in lines:
        stripped = line.strip()
        lower = stripped.lower()
        if lower.startswith("%moinp"):
            continue
        if stripped.startswith("!"):
            tokens = stripped.split()
            filtered = [tok for tok in tokens if tok.lower() != "moread"]
            if not filtered:
                continue
            rebuilt = " ".join(filtered)
            sanitized.append(rebuilt + ("\n" if not rebuilt.endswith("\n") else ""))
        else:
            sanitized.append(line if line.endswith("\n") else line + "\n")
    return sanitized


def _write_input_from_template(
    template_ctx,
    coords,
    output_path: Path | str,
    additions_text: str,
    geom_source_path: Path | str,
    config,
    main_basisset,
    metal_basisset,
    pal_override,
    maxcore_override,
) -> bool:
    if not template_ctx:
        return False
    expected = template_ctx["coord_count"]
    if expected and len(coords) != expected:
        logging.warning(
            "Template coordinate count mismatch (expected %d, got %d); falling back to generated input.",
            expected,
            len(coords),
        )
        return False

    lines = _sanitize_template_lines(list(template_ctx["lines"]))

    pal_val = int(pal_override) if pal_override is not None else int(config["PAL"])
    maxcore_val = int(maxcore_override) if maxcore_override is not None else int(config["maxcore"])

    for idx, line in enumerate(lines):
        lower = line.strip().lower()
        if lower.startswith("%pal"):
            lines[idx] = f"%pal nprocs {pal_val} end\n"
        elif lower.startswith("%maxcore"):
            lines[idx] = f"%maxcore {maxcore_val}\n"

    geom_start = None
    for idx, line in enumerate(lines):
        if line.strip().lower().startswith("* xyz"):
            geom_start = idx + 1
            break
    if geom_start is None:
        logging.warning("Sanitized template lost '* xyz' marker; fallback to generated input.")
        return False

    geom_end = geom_start
    while geom_end < len(lines) and lines[geom_end].strip() != "*":
        geom_end += 1
    if geom_end >= len(lines):
        logging.warning("Sanitized template missing terminal '*' marker; fallback to generated input.")
        return False
    coord_idx = 0

    raw_geom_lines = []
    for elem, x, y, z in coords:
        raw_geom_lines.append(f"{elem} {x} {y} {z}\n")

    found_metals_local = search_transition_metals(str(geom_source_path))
    main_sel, metal_sel = set_main_basisset(found_metals_local, config)
    metal_eff = metal_basisset or metal_sel
    enable_first = str(config.get("first_coordination_sphere_metal_basisset", "no")).lower() in (
        "yes",
        "true",
        "1",
        "on",
    )
    sphere_scale_raw = str(config.get("first_coordination_sphere_scale", "")).strip()
    radii_map = (
        _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009"))
        if (enable_first and not sphere_scale_raw)
        else None
    )

    geom_with_basis = _apply_per_atom_newgto(
        raw_geom_lines,
        found_metals_local,
        metal_eff,
        config,
        radii_map,
    )

    for idx in range(geom_start, geom_end):
        match = _COORD_LINE_RE.match(lines[idx].rstrip("\n"))
        if not match:
            continue
        if coord_idx >= len(geom_with_basis):
            logging.warning("Insufficient coordinates to populate IMAG template; fallback required.")
            return False
        lines[idx] = geom_with_basis[coord_idx]
        coord_idx += 1

    if coord_idx != len(geom_with_basis):
        logging.warning(
            "Template consumed %d coordinates but geometry list contained %d; falling back to generated input.",
            coord_idx,
            len(geom_with_basis),
        )
        return False

    additions_text = additions_text.strip()
    if additions_text:
        additions_lines = [ln if ln.endswith("\n") else ln + "\n" for ln in additions_text.splitlines()]
        insertion_index = geom_end
        lines = lines[:insertion_index] + additions_lines + lines[insertion_index:]

    try:
        with Path(output_path).open("w", encoding="utf-8") as fh:
            fh.writelines(lines)
        return True
    except Exception as exc:
        logging.error(f"Failed to write IMAG input '{output_path}' from template: {exc}")
        return False


def _normalize_additions_payload(additions) -> str:
    if not additions:
        return ""
    if isinstance(additions, str):
        candidate = additions.strip()
    if isinstance(additions, dict):
        chunks = []
        for key, value in additions.items():
            key_str = str(key).strip()
            val_str = str(value).strip()
            if key_str and val_str:
                chunks.append(f"{key_str}={val_str}")
        candidate = "\n".join(chunks)
    else:
        candidate = str(additions).strip()

    if not candidate:
        return ""

    filtered_lines = []
    for line in candidate.splitlines():
        if "%moinp" in line.lower():
            continue
        clean = line.strip()
        if clean:
            filtered_lines.append(clean)
    return "\n".join(filtered_lines)

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

def _build_bang_line_IMAG(
    config,
    rel_token,
    main_basisset,
    aux_jk,
    implicit,
    qmmm_method=None,
    *,
    include_freq: bool = True,
    geom_override: str | None = None,
):
    """
    Construct the '!' line for IMAG iterations:
      functional [REL] main_basis [disp] [ri_jkx] [aux_jk] [implicit] [geom_opt] FREQ initial_guess
    """
    ri_jkx = str(config.get("ri_jkx", "")).strip()
    disp   = str(config.get("disp_corr", "")).strip()
    geom   = geom_override if geom_override is not None else str(config.get("geom_opt", "OPT")).strip()
    init_tokens = str(config.get("initial_guess", "PModel")).split()
    initg  = init_tokens[0] if init_tokens else "PModel"

    tokens = ["!"]
    if qmmm_method:
        tokens.append(qmmm_method)
    tokens.append(str(config["functional"]))
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
    if include_freq:
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

def collect_imaginary_modes(log_file: str) -> List[Tuple[int, float]]:
    """Return list of (mode_index, frequency) for imaginary modes sorted by frequency."""
    modes: List[Tuple[int, float]] = []
    try:
        with open(log_file, "r", errors="ignore") as fh:
            for line in fh:
                if "***imaginary mode***" not in line:
                    continue
                m = re.search(r"^\s*(\d+):\s*([-+]?\d+(?:\.\d+)?)\s+cm\*\*-1", line)
                if m:
                    try:
                        idx = int(m.group(1))
                        freq = float(m.group(2))
                        modes.append((idx, freq))
                    except ValueError:
                        continue
    except FileNotFoundError:
        logging.error(f"Log file '{log_file}' not found when collecting imaginary modes.")
        return []
    return sorted(modes, key=lambda item: item[1])

@dataclass
class ImagModeCandidate:
    mode_index: int
    frequency: float
    energy: float
    geometry_path: Path
    sp_output: Path

def run_plotvib_mode(
    hess_file: Path,
    mode_index: int,
    *,
    workdir: Path,
    amplitude: float | None = None,
) -> Path:
    """Run orca_pltvib for a specific mode and return generated XYZ path."""
    hess_file = Path(hess_file).resolve()
    if not hess_file.exists():
        logging.error(f"Hessian file '{hess_file}' not found for mode {mode_index}.")
        raise FileNotFoundError(hess_file)

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    cmd = ["orca_pltvib", str(hess_file), str(mode_index)]
    if amplitude is not None and not math.isclose(amplitude, 1.0):
        cmd.append(f"{amplitude:.6f}")

    try:
        subprocess.run(cmd, cwd=str(workdir), check=True)
    except subprocess.CalledProcessError as exc:
        logging.error(f"orca_pltvib failed for '{hess_file}' mode {mode_index}: {exc}")
        raise

    xyz_name = f"{hess_file.name}.v{mode_index:03d}.xyz"
    xyz_path = workdir / xyz_name
    if not xyz_path.exists():
        logging.error(f"Expected vib geometry '{xyz_path}' not created by orca_pltvib.")
        raise FileNotFoundError(xyz_path)
    return xyz_path

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

def run_plotvib(iteration, workdir: str | Path | None = None):
    try:
        base_dir = Path(workdir) if workdir is not None else Path.cwd()
        hess_path = base_dir / f"input_{iteration}.hess"
        new_structure = run_plotvib_mode(hess_path, 6, workdir=base_dir)
        logging.info(f"plotvib run successful for 'input_{iteration}.hess'")
        return str(new_structure)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running plotvib: {e}")
        sys.exit(1)

def run_plotvib0(input_file, workdir: str | Path | None = None):
    try:
        base_dir = Path(workdir) if workdir is not None else Path.cwd()
        hess_path = Path(input_file).with_suffix(".hess")
        if not hess_path.is_absolute():
            hess_path = (Path.cwd() / hess_path).resolve()
        new_structure = run_plotvib_mode(hess_path, 6, workdir=base_dir)
        logging.info(f"plotvib run successful for '{input_file}.hess'")
        return str(new_structure)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running plotvib: {e}")
        sys.exit(1)

# ----------------------- main writer (updated) -----------------------

def read_and_modify_xyz_IMAG(
    input_file_path,
    output_file_path,
    charge,
    multiplicity,
    solvent,
    metals,
    config,
    main_basisset,
    metal_basisset,
    additions,
    pal_override=None,
    maxcore_override=None,
    include_freq: bool = True,
):
    """Construct an ORCA input for an IMAG iteration.

    The builder now mirrors the generic OCCUPIER writers, including support for
    QM/XTB splits and cached %QMMM blocks, while still allowing per-iteration
    overrides supplied via ``additions``.
    """
    try:
        with open(input_file_path, "r") as file:
            coord_lines = [ln for ln in file.readlines() if ln.strip()]
        coord_lines = _strip_xyz_header(coord_lines)
        coord_lines = _trim_xyz_columns(coord_lines)
    except Exception as e:
        logging.error(f"Error reading '{input_file_path}': {e}")
        sys.exit(1)

    # Determine metals from the current structure
    found_metals_local = search_transition_metals(input_file_path)

    # Resolve basis settings
    main_sel, metal_sel = set_main_basisset(found_metals_local, config)
    main_eff = main_basisset or main_sel
    metal_eff = metal_basisset or metal_sel

    # Relativity/AUX policy (3d → non-rel + aux_jk; 4d/5d → rel + aux_jk_rel)
    rel_token, aux_jk_token, _ = select_rel_and_aux(found_metals_local, config)

    # implicit solvation
    implicit = ""
    model = str(config.get("implicit_solvation_model", "")).strip()
    if model:
        implicit = f"{model}({solvent})" if solvent else model

    # QM/XTB partition handling
    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(coord_lines, Path(input_file_path))
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = "QM/XTB" if qmmm_range else None

    # Load radii for optional first coordination sphere tagging when needed
    enable_first = str(config.get("first_coordination_sphere_metal_basisset", "no")).lower() in (
        "yes",
        "true",
        "1",
        "on",
    )
    sphere_scale_raw = str(config.get("first_coordination_sphere_scale", "")).strip()
    radii_map = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009")) if (enable_first and not sphere_scale_raw) else None

    # '!' line for IMAG (always FREQ)
    geom_override = "SP" if not include_freq else None
    bang = _build_bang_line_IMAG(
        config,
        rel_token,
        main_eff,
        aux_jk_token,
        implicit,
        qmmm_token,
        include_freq=include_freq,
        geom_override=geom_override,
    )

    output_blocks = collect_output_blocks(config)
    builder = OrcaInputBuilder(bang)
    pal_val = int(pal_override) if pal_override is not None else int(config["PAL"])
    maxcore_val = int(maxcore_override) if maxcore_override is not None else int(config["maxcore"])
    builder.add_resources(maxcore_val, pal_val, resolve_maxiter(config))
    builder.add_additions(additions)

    if include_freq:
        # Add %freq block with temperature (IMAG always uses FREQ)
        from .xyz_io import _build_freq_block

        freq_block = _build_freq_block(config)
        builder.add_block(freq_block)
    builder.add_blocks(output_blocks)

    lines = builder.lines
    lines.extend(build_qmmm_block(qmmm_range))
    lines.append(f"* xyz {charge} {multiplicity}\n")

    geom = _apply_per_atom_newgto(geom_lines, found_metals_local, metal_eff, config, radii_map)
    lines.extend(geom)
    if not lines or not lines[-1].strip() == "*":
        lines.append("*\n")

    try:
        with open(output_file_path, "w") as f:
            f.writelines(lines)
        logging.info(f"Input file '{output_file_path}' created successfully.")
    except Exception as e:
        logging.error(f"Error writing '{output_file_path}': {e}")
        sys.exit(1)

# ----------------------- rest of the logic -----------------------

def extract_structure(input_file, iteration):
    try:
        source_path = Path(input_file).resolve()
        with source_path.open('r', encoding='utf-8') as file:
            lines = file.readlines()
        logging.info(f"Process file: {input_file}")
        star_indices = [i for i, line in enumerate(lines) if '*' in line]
        logging.info(f"Number of found '*' lines: {len(star_indices)}")
        if len(star_indices) != 20:
            raise ValueError(f"File '{input_file}' does not contain exactly 20 '*' lines.")
        start_index = star_indices[4] - 1
        end_index = star_indices[5]
        extracted_lines = lines[start_index:end_index]
        output_path = source_path.parent / f"input_{iteration}_structure_5.xyz"
        with output_path.open('w', encoding='utf-8') as file:
            file.writelines(extracted_lines)
        logging.info(f"Structure extracted to '{output_path}'")
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

def run_IMAG(
    input_file,
    hess_file,
    charge,
    multiplicity,
    solvent,
    metals,
    config,
    main_basisset,
    metal_basisset,
    additions,
    step_name="initial",
    source_input=None,
    pal_override=None,
    maxcore_override=None,
):
    """Run IMAG elimination for a given calculation step.

    Args:
        step_name: Name of the calculation step (e.g., "initial", "red_step_1", "ox_step_2")
    """
    if str(config.get("IMAG", "no")).lower() != "yes":
        return

    # Check IMAG_scope setting
    imag_scope = str(config.get("IMAG_scope", "initial")).lower()
    if imag_scope == "initial" and step_name != "initial":
        logging.info(f"Skipping IMAG for '{step_name}' (IMAG_scope=initial)")
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
    original_out = Path(input_file)
    template_ctx = _load_inp_template(source_input) if source_input else None

    if source_input:
        src_pal, src_maxcore = _extract_resources_from_input(source_input)
        if src_pal is not None:
            pal_override = src_pal
        if src_maxcore is not None:
            maxcore_override = src_maxcore

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

    imag_folder.mkdir(parents=True, exist_ok=True)
    def _write_imag_input(
        geometry_source: Path,
        destination_path: Path,
        include_freq: bool,
        additions_payload: str,
    ) -> str:
        """Create ORCA input from geometry, using template when available."""
        coords = _extract_xyz_coordinates(str(geometry_source))
        used_template = False
        if include_freq and template_ctx and coords:
            used_template = _write_input_from_template(
                template_ctx,
                coords,
                destination_path,
                additions_payload,
                geometry_source,
                config,
                main_basisset,
                metal_basisset,
                pal_override,
                maxcore_override,
            )

        if not used_template:
            read_and_modify_xyz_IMAG(
                str(geometry_source),
                str(destination_path),
                charge,
                multiplicity,
                solvent,
                metals,
                config,
                main_basisset,
                metal_basisset,
                additions_payload,
                pal_override=pal_override,
                maxcore_override=maxcore_override,
                include_freq=include_freq,
            )
        return additions_payload

    if original_out.exists():
        try:
            shutil.copy2(original_out, imag_folder / f"{step_name}_0.out")
        except Exception as exc:
            logging.warning(f"Could not archive original output '{original_out}': {exc}")
    if template_ctx and template_ctx.get("path") and template_ctx["path"].exists():
        try:
            shutil.copy2(template_ctx["path"], imag_folder / f"{step_name}_0.inp")
        except Exception as exc:
            logging.debug(f"Failed to archive original input '{template_ctx['path']}': {exc}")

    additions_base = _normalize_additions_payload(additions)

    best_energy = find_electronic_energy(str(original_out))
    if best_energy is None:
        logging.warning("Could not extract FINAL SINGLE POINT ENERGY from original output; aborting IMAG.")
        return

    displacement_scale_raw = config.get("IMAG_displacement_scale", 1.0)
    try:
        displacement_scale = float(displacement_scale_raw)
        if displacement_scale <= 0:
            raise ValueError
    except (TypeError, ValueError):
        logging.warning("IMAG_displacement_scale must be a positive number; defaulting to 1.0.")
        displacement_scale = 1.0

    energy_tol_raw = config.get("IMAG_energy_tol", ENERGY_IMPROVEMENT_TOL)
    try:
        energy_tol = float(energy_tol_raw)
        if energy_tol <= 0:
            raise ValueError
    except (TypeError, ValueError):
        logging.warning("IMAG_energy_tol must be a positive number; defaulting to %.1e.", ENERGY_IMPROVEMENT_TOL)
        energy_tol = ENERGY_IMPROVEMENT_TOL

    current_log_path = original_out
    current_hess_path = Path(hess_file).with_suffix(".hess")
    if not current_hess_path.is_absolute():
        current_hess_path = (Path.cwd() / current_hess_path).resolve()
    if not current_hess_path.exists():
        logging.error(f"Hessian file '{current_hess_path}' not found; cannot run IMAG.")
        return

    iteration = 0
    last_success_iteration = 0
    additions_eff = additions_base

    while True:
        modes = [
            (idx, freq)
            for idx, freq in collect_imaginary_modes(str(current_log_path))
            if freq < threshold
        ]
        if not modes:
            logging.info("All imaginary frequencies within threshold; stopping IMAG loop.")
            break

        best_candidate: ImagModeCandidate | None = None

        for mode_index, freq in modes:
            try:
                vib_xyz_path = run_plotvib_mode(
                    current_hess_path,
                    mode_index,
                    workdir=imag_folder,
                    amplitude=displacement_scale,
                )
            except Exception:
                logging.debug(f"Mode {mode_index}: failed to generate displacement geometry.", exc_info=True)
                continue

            sp_label = f"iter{iteration}_mode{mode_index:03d}"
            sp_input_path = imag_folder / f"{sp_label}_sp.inp"
            sp_output_path = imag_folder / f"{sp_label}_sp.out"

            for tmp in (sp_input_path, sp_output_path):
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except OSError:
                        pass

            _write_imag_input(vib_xyz_path, sp_input_path, include_freq=False, additions_payload=additions_eff)
            sp_success = run_orca(
                str(sp_input_path),
                str(sp_output_path),
                working_dir=imag_folder,
            )
            if not sp_success:
                logging.warning(f"Mode {mode_index}: single-point calculation failed; skipping.")
                continue

            trial_energy = find_electronic_energy(str(sp_output_path))
            if trial_energy is None:
                logging.warning(f"Mode {mode_index}: could not extract single-point energy; skipping.")
                continue

            if trial_energy < best_energy - energy_tol:
                if best_candidate is None or trial_energy < best_candidate.energy - energy_tol:
                    best_candidate = ImagModeCandidate(
                        mode_index=mode_index,
                        frequency=freq,
                        energy=trial_energy,
                        geometry_path=vib_xyz_path,
                        sp_output=sp_output_path,
                    )

        if best_candidate is None:
            logging.info("No imaginary-mode displacement lowered the energy; stopping IMAG.")
            break

        iteration += 1
        freq_input_path = imag_folder / f"input_{iteration}.inp"
        freq_output_path = imag_folder / f"output_{iteration}.out"
        for tmp in (freq_input_path, freq_output_path):
            if tmp.exists():
                try:
                    tmp.unlink()
                except OSError:
                    pass

        additions_current = additions_eff
        _write_imag_input(
            best_candidate.geometry_path,
            freq_input_path,
            include_freq=True,
            additions_payload=additions_current,
        )
        success = run_orca_IMAG(str(freq_input_path), iteration, working_dir=imag_folder)
        if not success:
            geometry_mismatch = False
            log_file = imag_folder / f"output_{iteration}.out"
            if log_file.exists():
                try:
                    with log_file.open("r", errors="ignore") as fh:
                        if "Input geometry does not match current geometry" in fh.read():
                            geometry_mismatch = True
                except Exception as exc:
                    logging.debug(f"Iteration {iteration}: failed to inspect IMAG log: {exc}")

            if geometry_mismatch and additions_current:
                logging.warning(
                    f"Iteration {iteration}: geometry mismatch detected; retrying without supplemental additions."
                )
                additions_current = ""
                additions_eff = ""
                additions_base = ""
                _write_imag_input(
                    best_candidate.geometry_path,
                    freq_input_path,
                    include_freq=True,
                    additions_payload=additions_current,
                )
                success = run_orca_IMAG(str(freq_input_path), iteration, working_dir=imag_folder)

            if not success:
                logging.warning(f"Iteration {iteration}: frequency calculation failed; stopping IMAG loop.")
                break

        current_log_path = imag_folder / f"output_{iteration}.out"
        new_energy = find_electronic_energy(str(current_log_path))
        if new_energy is None:
            logging.warning(
                f"Iteration {iteration}: no FINAL SINGLE POINT ENERGY found; keeping single-point estimate."
            )
            new_energy = best_candidate.energy
        best_energy = new_energy
        last_success_iteration = iteration

        current_hess_path = (imag_folder / f"input_{iteration}.hess").resolve()
        if not current_hess_path.exists():
            logging.warning(f"Iteration {iteration}: Hessian file missing; terminating IMAG loop.")
            break

    if last_success_iteration == 0:
        logging.info("IMAG finished without improvements; original files remain unchanged.")
        return

    final_log_file = imag_folder / f"output_{last_success_iteration}.out"
    final_xyz_file = imag_folder / f"input_{last_success_iteration}.xyz"
    input_path = Path(input_file)
    destination_folder = input_path.parent if input_path.is_absolute() else Path.cwd()
    destination_log = input_path if input_path.is_absolute() else destination_folder / input_path
    destination_structure = destination_folder / f"{hess_file}.xyz"
    if final_log_file.exists():
        try:
            shutil.copy2(final_log_file, destination_log)
            print(
                f"Log file 'output_{last_success_iteration}.out' copied back as '{destination_log.name}'."
            )
        except Exception as e:
            logging.warning(f"copy back log: {e}")
    else:
        print(f"ERROR: Log file 'output_{last_success_iteration}.out' not found.")
    if final_xyz_file.exists():
        try:
            shutil.copy2(final_xyz_file, destination_structure)
            print(
                f"Structure file 'input_{last_success_iteration}.xyz' copied back as '{hess_file}.xyz'."
            )
        except Exception as e:
            logging.warning(f"copy back xyz: {e}")
    else:
        print(f"ERROR: Structure file 'input_{last_success_iteration}.xyz' not found.")
