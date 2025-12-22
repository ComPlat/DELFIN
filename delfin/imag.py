import logging
import math
import os
import re
import shutil
import subprocess
import sys
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional

from delfin.common.orca_blocks import OrcaInputBuilder, collect_output_blocks, resolve_maxiter
from delfin.dynamic_pool import PoolJob, JobPriority
from delfin.global_manager import get_global_manager

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


def _coerce_int(value, fallback: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return fallback

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
    def _is_coord_line(parts: list[str]) -> bool:
        """Return True if parts look like 'Elem x y z' with numeric coords."""
        if len(parts) < 4 or not parts[0] or not parts[0][0].isalpha():
            return False
        try:
            float(parts[1])
            float(parts[2])
            float(parts[3])
            return True
        except Exception:
            return False

    try:
        int(first)
        working = working[1:]
        if working:
            head = working[0].strip().split()
            if not _is_coord_line(head):
                working = working[1:]
        return working
    except ValueError:
        lower = first.lower()
        if lower.startswith("coordinates from orca-job") or lower.startswith("coordinates from orca job"):
            return working[1:]
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
    *,
    geom_override: str | None = None,
    include_freq: bool = True,
) -> bool:
    if not template_ctx:
        return False
    expected = template_ctx["coord_count"]
    if expected:
        actual_len = len(coords)
        if actual_len < expected:
            logging.warning(
                "Template coordinate count mismatch (expected %d, got %d); falling back to generated input.",
                expected,
                actual_len,
            )
            return False
        if actual_len > expected:
            logging.warning(
                "Template coordinate count mismatch (expected %d, got %d); trimming coordinates for IMAG template.",
                expected,
                actual_len,
            )
            coords = coords[:expected]

    lines = _sanitize_template_lines(list(template_ctx["lines"]))

    pal_val = int(pal_override) if pal_override is not None else int(config["PAL"])
    maxcore_val = int(maxcore_override) if maxcore_override is not None else int(config["maxcore"])

    for idx, line in enumerate(lines):
        lower = line.strip().lower()
        if lower.startswith("%pal"):
            lines[idx] = f"%pal nprocs {pal_val} end\n"
        elif lower.startswith("%maxcore"):
            lines[idx] = f"%maxcore {maxcore_val}\n"
        elif line.strip().startswith("!"):
            geom_tokens = str(config.get("geom_opt", "")).split()
            override_tokens = [tok for tok in (geom_override.split() if geom_override else []) if tok]
            existing_tokens = line.strip().split()
            tokens = list(existing_tokens)
            if geom_override is not None:
                tokens = [tok for tok in tokens if tok not in geom_tokens]
                if override_tokens:
                    tokens.extend(override_tokens)
            if include_freq and "FREQ" not in {tok.upper() for tok in tokens}:
                tokens.append("FREQ")
            if not include_freq:
                tokens = [tok for tok in tokens if tok.upper() != "FREQ"]
            lines[idx] = " ".join(tokens)
            if not lines[idx].endswith("\n"):
                lines[idx] += "\n"

    if not include_freq:
        cleaned: list[str] = []
        skipping_freq = False
        for line in lines:
            stripped = line.strip().lower()
            if stripped.startswith("%freq"):
                skipping_freq = True
                continue
            if skipping_freq:
                if stripped == "end":
                    skipping_freq = False
                continue
            cleaned.append(line)
        lines = cleaned

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
        existing_normalized = {line.strip().lower() for line in lines}
        additions_lines = [ln for ln in additions_lines if ln.strip().lower() not in existing_normalized]
        if additions_lines:
            insertion_index = geom_start - 1 if geom_start > 0 else 0
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
    direction: str
    optimized_geometry: Path | None = None


@dataclass
class _IMAGCandidateJob:
    label: str
    mode_index: int
    frequency: float
    geometry_path: Path
    sp_input_path: Path
    sp_output_path: Path
    direction: str


@dataclass
class _ImagSinglePointResult:
    job: _IMAGCandidateJob
    success: bool = False
    energy: Optional[float] = None
    optimized_geometry: Optional[Path] = None
    message: Optional[str] = None


def _locate_relaxed_xyz(sp_input_path: Path, sp_output_path: Path) -> Optional[Path]:
    candidates = [
        Path(sp_input_path).with_suffix(".xyz"),
        Path(sp_output_path).with_suffix(".xyz"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _run_sp_candidates_sequential(
    jobs: List[_IMAGCandidateJob],
) -> List[_ImagSinglePointResult]:
    results: List[_ImagSinglePointResult] = []
    recalc = str(os.environ.get("DELFIN_RECALC", "0")).lower() in ("1", "true", "yes", "on")

    for job in jobs:
        result = _ImagSinglePointResult(job=job)

        # RECALC mode: Check if output already exists and is valid
        if recalc and job.sp_output_path.exists():
            try:
                if job.sp_output_path.stat().st_size >= 100:
                    with job.sp_output_path.open("r", encoding="utf-8", errors="replace") as f:
                        content = f.read()
                        if OK_MARKER in content:
                            energy = find_electronic_energy(str(job.sp_output_path))
                            if energy is not None:
                                logging.info(f"[recalc] IMAG SP mode {job.mode_index} ({job.direction}): skip (energy={energy})")
                                result.success = True
                                result.energy = energy
                                result.optimized_geometry = _locate_relaxed_xyz(job.sp_input_path, job.sp_output_path)
                                results.append(result)
                                continue
            except Exception as e:
                logging.debug(f"[recalc] Failed to check existing IMAG SP output for mode {job.mode_index}: {e}")

        sp_success = run_orca(
            str(job.sp_input_path),
            str(job.sp_output_path),
            working_dir=job.sp_output_path.parent,
        )
        if not sp_success:
            msg = f"Mode {job.mode_index}: single-point calculation failed; skipping."
            logging.warning(msg)
            result.message = msg
        else:
            energy = find_electronic_energy(str(job.sp_output_path))
            if energy is None:
                msg = f"Mode {job.mode_index}: could not extract single-point energy; skipping."
                logging.warning(msg)
                result.message = msg
            else:
                result.success = True
                result.energy = energy
        result.optimized_geometry = _locate_relaxed_xyz(job.sp_input_path, job.sp_output_path)
        results.append(result)
    return results


def _run_sp_candidates_parallel(
    jobs: List[_IMAGCandidateJob],
    pool,
    *,
    core_allocation: int,
    maxcore_value: int,
) -> List[_ImagSinglePointResult]:
    handles: List[tuple[_ImagSinglePointResult, threading.Event]] = []
    recalc = str(os.environ.get("DELFIN_RECALC", "0")).lower() in ("1", "true", "yes", "on")

    for job in jobs:
        result = _ImagSinglePointResult(job=job)
        event = threading.Event()

        # RECALC mode: Check if output already exists and is valid before submitting to pool
        if recalc and job.sp_output_path.exists():
            try:
                if job.sp_output_path.stat().st_size >= 100:
                    with job.sp_output_path.open("r", encoding="utf-8", errors="replace") as f:
                        content = f.read()
                        if OK_MARKER in content:
                            energy = find_electronic_energy(str(job.sp_output_path))
                            if energy is not None:
                                logging.info(f"[recalc] IMAG SP mode {job.mode_index} ({job.direction}): skip (energy={energy})")
                                result.success = True
                                result.energy = energy
                                result.optimized_geometry = _locate_relaxed_xyz(job.sp_input_path, job.sp_output_path)
                                event.set()
                                handles.append((result, event))
                                continue
            except Exception as e:
                logging.debug(f"[recalc] Failed to check existing IMAG SP output for mode {job.mode_index}: {e}")

        handles.append((result, event))

        def runner(*_args, cur_job=job, cur_result=result, cur_event=event, **kwargs) -> None:
            try:
                sp_success = run_orca(
                    str(cur_job.sp_input_path),
                    str(cur_job.sp_output_path),
                    working_dir=cur_job.sp_output_path.parent,
                )
                if not sp_success:
                    msg = f"Mode {cur_job.mode_index}: single-point calculation failed; skipping."
                    logging.warning(msg)
                    cur_result.message = msg
                    return

                energy = find_electronic_energy(str(cur_job.sp_output_path))
                if energy is None:
                    msg = f"Mode {cur_job.mode_index}: could not extract single-point energy; skipping."
                    logging.warning(msg)
                    cur_result.message = msg
                    return

                cur_result.success = True
                cur_result.energy = energy
            except Exception as exc:  # noqa: BLE001
                logging.error("Mode %s: parallel single-point run crashed: %s", cur_job.mode_index, exc, exc_info=True)
                cur_result.message = str(exc)
            finally:
                if cur_result.optimized_geometry is None:
                    cur_result.optimized_geometry = _locate_relaxed_xyz(cur_job.sp_input_path, cur_job.sp_output_path)
                cur_event.set()

        job_id = f"IMAG_SP::{job.label}"
        cores_requested = max(1, core_allocation)
        mem_mb = max(256, maxcore_value) * cores_requested
        pool_job = PoolJob(
            job_id=job_id,
            cores_min=cores_requested,
            cores_optimal=cores_requested,
            cores_max=cores_requested,
            memory_mb=mem_mb,
            priority=JobPriority.NORMAL,
            execute_func=runner,
            args=(),
            kwargs={},
            estimated_duration=1800.0,
            working_dir=job.sp_output_path.parent,
        )
        pool_job.suppress_pool_logs = True
        pool.submit_job(pool_job)

    results: List[_ImagSinglePointResult] = []
    for result, event in handles:
        event.wait()
        if result.optimized_geometry is None:
            result.optimized_geometry = _locate_relaxed_xyz(
                result.job.sp_input_path,
                result.job.sp_output_path,
            )
        results.append(result)
    return results


def _execute_sp_candidates(
    jobs: List[_IMAGCandidateJob],
    *,
    pool,
    core_allocation: int,
    maxcore_value: int,
) -> List[_ImagSinglePointResult]:
    if not jobs:
        return []
    use_pool = False
    parent_job_id = None
    effective_core_allocation = core_allocation

    if pool is not None:
        # Check if we're running within a pool job (potential parent for borrowing)
        from delfin.dynamic_pool import get_current_job_id
        parent_job_id = get_current_job_id()

        try:
            status = pool.get_status()
        except Exception:
            status = None

        if status:
            total = status.get("total_cores", 0)
            allocated = status.get("allocated_cores", 0)
            queued = status.get("queued_jobs", 0)
            running = status.get("running_jobs", 0)
            available = total - allocated if total else 0
            cores_needed = len(jobs) * core_allocation

            # If we have a parent job, IMAG child jobs can borrow cores from parent
            # This makes parallel execution much more viable
            if parent_job_id is not None:
                # If nothing is free, run sequential to avoid deadlock on the parent's pool slot
                if available <= 0:
                    logging.debug(
                        f"[IMAG] Falling back to sequential (parent={parent_job_id}): "
                        f"0/{total} cores available for {len(jobs)} SP jobs"
                    )
                else:
                    # Parent exists - child jobs can borrow cores. Cap to available so we start instead of hanging.
                    effective_core_allocation = max(1, min(core_allocation, available))
                    if running < pool.max_concurrent_jobs:
                        use_pool = True
                        logging.debug(
                            f"[IMAG] Using pool with parent-child borrowing for {len(jobs)} SP jobs: "
                            f"parent_job={parent_job_id}, {available}/{total} pool cores available, "
                            f"alloc_per_job={effective_core_allocation}"
                        )
                    else:
                        logging.debug(
                            f"[IMAG] Falling back to sequential (parent={parent_job_id}): "
                            f"{available}/{total} cores, {running}/{pool.max_concurrent_jobs} running"
                        )
            else:
                # No parent job - use conservative approach to avoid deadlock
                # Require significant headroom (at least 25% of total cores)
                min_available = max(cores_needed, total // 4) if total else cores_needed

                if available >= min_available and running < pool.max_concurrent_jobs:
                    use_pool = True
                    logging.debug(
                        f"[IMAG] Using pool for {len(jobs)} SP calculations: {available}/{total} cores available, "
                        f"{cores_needed} cores needed, {running}/{pool.max_concurrent_jobs} jobs running"
                    )
                else:
                    logging.debug(
                        f"[IMAG] Falling back to sequential execution to avoid deadlock: "
                        f"{available}/{total} cores available, {cores_needed} cores needed for {len(jobs)} jobs, "
                        f"{running}/{pool.max_concurrent_jobs} jobs running, {queued} jobs queued"
                    )
        else:
            use_pool = True

    if not use_pool:
        return _run_sp_candidates_sequential(jobs)
    return _run_sp_candidates_parallel(
        jobs,
        pool,
        core_allocation=effective_core_allocation,
        maxcore_value=maxcore_value,
    )


def _inject_moinp_block(input_path: Path, gbw_path: Path) -> None:
    """Inject %moinp block and ensure MOREAD keyword for the supplied GBW file."""
    input_path = Path(input_path)
    gbw_path = Path(gbw_path)
    try:
        lines = input_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except FileNotFoundError:
        logging.error(f"Input file '{input_path}' missing when adding %moinp.")
        return

    rel_gbw = os.path.relpath(gbw_path, start=input_path.parent).replace("\\", "/")
    moinp_line = f'%moinp "{rel_gbw}"'

    if not any(line.strip().startswith("%moinp") for line in lines):
        insert_idx = 0
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("%maxcore"):
                insert_idx = idx
                break
            if stripped.startswith("%") and insert_idx == 0:
                insert_idx = idx
        lines.insert(insert_idx, moinp_line)
    else:
        lines = [moinp_line if line.strip().startswith("%moinp") else line for line in lines]

    for idx, line in enumerate(lines):
        if line.strip().startswith("!"):
            # Add MOREAD if not present (required when using MOINP with .gbw)
            if "MOREAD" not in line:
                lines[idx] = line.rstrip() + " MOREAD\n"
            break

    input_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _extract_structure_section(lines: list[str], structure_number: int) -> list[str]:
    """Return the coordinate block for the requested structure index (1-based)."""
    star_indices = [i for i, line in enumerate(lines) if "*" in line]
    if len(star_indices) < structure_number + 1:
        raise ValueError(
            f"Displaced geometry contains only {len(star_indices) - 1} structures; requested {structure_number}."
        )
    start_index = star_indices[structure_number - 1] - 1
    end_index = star_indices[structure_number]
    return lines[start_index:end_index]


def _extract_structure_xyz(source_xyz: Path, destination: Path, structure_number: int) -> Path:
    """Extract the specified structure fragment into destination."""
    source_xyz = Path(source_xyz)
    try:
        with source_xyz.open("r", encoding="utf-8") as fh:
            lines = fh.readlines()
    except FileNotFoundError as exc:
        logging.error(f"Unable to read displaced geometry '{source_xyz}': {exc}")
        raise

    extracted = _extract_structure_section(lines, structure_number)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as out:
        out.writelines(extracted)
    return destination


def run_plotvib_mode(
    hess_file: Path,
    mode_index: int,
    *,
    workdir: Path,
    amplitude: float | None = None,
) -> dict[str, Path]:
    """Run orca_pltvib for a specific mode and return {pos, neg} structure_5/15 geometries."""
    hess_file = Path(hess_file).resolve()
    if not hess_file.exists():
        logging.error(f"Hessian file '{hess_file}' not found for mode {mode_index}.")
        raise FileNotFoundError(hess_file)

    workdir_path = Path(workdir)
    workdir_path.mkdir(parents=True, exist_ok=True)

    cmd = ["orca_pltvib", str(hess_file), str(mode_index)]
    if amplitude is not None and not math.isclose(amplitude, 1.0):
        cmd.append(f"{amplitude:.6f}")

    try:
        subprocess.run(cmd, cwd=str(workdir_path), check=True)
    except subprocess.CalledProcessError as exc:
        logging.error(f"orca_pltvib failed for '{hess_file}' mode {mode_index}: {exc}")
        raise

    xyz_name = f"{hess_file.name}.v{mode_index:03d}.xyz"
    displaced_xyz = workdir_path / xyz_name

    if not displaced_xyz.exists():
        alt_path = hess_file.parent / xyz_name
        if alt_path.exists():
            displaced_xyz = alt_path
        else:
            logging.error(f"Expected vib geometry '{xyz_name}' not created by orca_pltvib.")
            raise FileNotFoundError(displaced_xyz)

    try:
        pos_path = _extract_structure_xyz(displaced_xyz, workdir_path / f"{xyz_name}.structure_5.xyz", 5)
        neg_path = _extract_structure_xyz(displaced_xyz, workdir_path / f"{xyz_name}.structure_15.xyz", 15)
    except Exception as exc:
        logging.error(f"Failed to extract displaced structures from '{displaced_xyz}': {exc}")
        raise
    return {"pos": pos_path, "neg": neg_path}

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
        structures = run_plotvib_mode(hess_path, 6, workdir=base_dir)
        new_structure = structures["pos"]
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
        structures = run_plotvib_mode(hess_path, 6, workdir=base_dir)
        new_structure = structures["pos"]
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
    geom_override: str | None = None,
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

    # QM/MM partition handling
    geom_lines, qmmm_range, qmmm_explicit = split_qmmm_sections(coord_lines, Path(input_file_path))
    _ensure_qmmm_implicit_model(config, qmmm_range, qmmm_explicit)
    qmmm_token = str(config.get('qmmm_option', 'QM/XTB')).strip() if qmmm_range else None

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
    *,
    manager=None,
):
    """Run IMAG elimination for a given calculation step.

    Args:
        step_name: Name of the calculation step (e.g., "initial", "red_step_1", "ox_step_2")
        manager: Optional _WorkflowManager to submit IMAG jobs to for parallel execution
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

    cwd_snapshot = Path.cwd().resolve(strict=False)

    def _resolve_path(candidate: str | os.PathLike[str],
                      *,
                      preferred_dir: Path | None = None,
                      require_exists: bool = False) -> Path:
        path_obj = Path(candidate)
        if path_obj.is_absolute():
            return path_obj

        search_dirs: list[Path] = []
        if preferred_dir is not None:
            preferred_dir = Path(preferred_dir).resolve(strict=False)
            search_dirs.append(preferred_dir)
            for parent in preferred_dir.parents:
                if parent not in search_dirs:
                    search_dirs.append(parent)

        if cwd_snapshot not in search_dirs:
            search_dirs.append(cwd_snapshot)
        for parent in cwd_snapshot.parents:
            if parent not in search_dirs:
                search_dirs.append(parent)

        fallback_candidate: Path | None = None
        for base in search_dirs:
            candidate_abs = (base / path_obj).resolve(strict=False)
            if candidate_abs.exists():
                return candidate_abs
            if fallback_candidate is None:
                fallback_candidate = candidate_abs

        if not require_exists and fallback_candidate is not None:
            return fallback_candidate
        return fallback_candidate if fallback_candidate is not None else (cwd_snapshot / path_obj).resolve(strict=False)

    source_input_path: Path | None = None
    preferred_dir: Path | None = None
    if source_input:
        source_input_path = _resolve_path(source_input, require_exists=False)
        preferred_dir = source_input_path.parent

    input_path = _resolve_path(input_file, preferred_dir=preferred_dir, require_exists=True)
    first_freq = search_imaginary_mode2(str(input_path))
    if first_freq is None or first_freq >= threshold:
        logging.info(f"Imag within tolerance (freq={first_freq}, thr={threshold}) -> no IMAG.")
        return

    imag_folder = input_path.parent / f"{hess_file}_IMAG"
    original_out = input_path
    template_ctx = _load_inp_template(source_input_path) if source_input_path else None

    if source_input_path:
        src_pal, src_maxcore = _extract_resources_from_input(source_input_path)
        if pal_override is None and src_pal is not None:
            pal_override = src_pal
        if maxcore_override is None and src_maxcore is not None:
            maxcore_override = src_maxcore

    # RECALC mode: Check if IMAG was already successfully completed
    if recalc and imag_folder.is_dir():
        last_i, last_out = _find_last_ok_iteration(imag_folder)
        if last_i is not None and last_out and _imag_resolved(last_out, threshold):
            logging.info(f"[recalc] IMAG complete at iter {last_i}; skip.")
            final_log = Path(last_out)
            final_xyz = imag_folder / f"input_{last_i}.xyz"
            parent = imag_folder.resolve().parent
            dest_log = input_path
            dest_xyz = parent / f"{hess_file}.xyz"
            try: shutil.copy2(final_log, dest_log)
            except Exception as e: logging.warning(f"copy log: {e}")
            if final_xyz.exists():
                try: shutil.copy2(final_xyz, dest_xyz)
                except Exception as e: logging.warning(f"copy xyz: {e}")
            return
        # If IMAG folder exists but wasn't successful, check if we have partial work
        # (don't delete it like before - this allows proper resume functionality)
        if last_i is not None:
            logging.info(f"[recalc] IMAG incomplete at iter {last_i}; resuming from existing work.")
        else:
            # Check if any single-point calculations exist (iter0_mode*_sp.out)
            has_sp_jobs = False
            try:
                for entry in imag_folder.iterdir():
                    if entry.is_file() and entry.name.startswith("iter") and "_sp.out" in entry.name:
                        has_sp_jobs = True
                        break
            except Exception:
                pass

            if has_sp_jobs:
                logging.info(f"[recalc] IMAG folder contains existing single-point calculations; resuming.")
            else:
                logging.info(f"[recalc] IMAG folder exists but no work found; cleaning and restarting.")
                try: shutil.rmtree(imag_folder)
                except Exception: pass

    imag_folder.mkdir(parents=True, exist_ok=True)
    def _write_imag_input(
        geometry_source: Path,
        destination_path: Path,
        include_freq: bool,
        additions_payload: str,
        *,
        geom_override: str | None = None,
        pal_override_local: int | None = None,
        maxcore_override_local: int | None = None,
    ) -> str:
        """Create ORCA input from geometry, using template when available."""
        coords = _extract_xyz_coordinates(str(geometry_source))
        used_template = False
        if template_ctx and coords:
            used_template = _write_input_from_template(
                template_ctx,
                coords,
                destination_path,
                additions_payload,
                geometry_source,
                config,
                main_basisset,
                metal_basisset,
                pal_override_local if pal_override_local is not None else pal_override,
                maxcore_override_local if maxcore_override_local is not None else maxcore_override,
                geom_override=geom_override,
                include_freq=include_freq,
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
                pal_override=pal_override_local if pal_override_local is not None else pal_override,
                maxcore_override=maxcore_override_local if maxcore_override_local is not None else maxcore_override,
                geom_override=geom_override,
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

    sp_window_raw = config.get(
        "IMAG_sp_energy_window",
        config.get("IMAG_energy_tol", ENERGY_IMPROVEMENT_TOL),
    )
    try:
        sp_energy_window = float(sp_window_raw)
        if sp_energy_window <= 0:
            raise ValueError
    except (TypeError, ValueError):
        logging.warning(
            "IMAG_sp_energy_window must be a positive number; defaulting to %.1e.",
            ENERGY_IMPROVEMENT_TOL,
        )
        sp_energy_window = ENERGY_IMPROVEMENT_TOL

    optimize_candidates = str(config.get("IMAG_optimize_candidates", "no")).strip().lower() in (
        "yes",
        "true",
        "1",
        "on",
    )

    current_log_path = original_out
    current_hess_path = Path(hess_file).with_suffix(".hess")
    if not current_hess_path.is_absolute():
        current_hess_path = (Path.cwd() / current_hess_path).resolve()
    if not current_hess_path.exists():
        logging.error(f"Hessian file '{current_hess_path}' not found; cannot run IMAG.")
        return

    maxcore_effective = _coerce_int(
        maxcore_override,
        _coerce_int(config.get("maxcore"), 1000),
    )
    shared_pool = None
    try:
        manager = get_global_manager()
        manager.ensure_initialized(config)
        if manager.is_initialized():
            shared_pool = manager.get_pool()
    except Exception:
        logging.debug("Global job manager unavailable for IMAG single-point pool runs; falling back to sequential.", exc_info=True)
        shared_pool = None

    resolved_pal: Optional[int]
    if pal_override is not None:
        try:
            resolved_pal = max(1, int(pal_override))
        except (TypeError, ValueError):
            resolved_pal = max(1, _coerce_int(config.get("PAL"), 1))
    else:
        core_share = _coerce_int(config.get("PAL"), 1)
        if shared_pool is not None:
            try:
                total_cores = max(1, getattr(shared_pool, "total_cores", core_share))
                concurrent = max(1, getattr(shared_pool, "max_concurrent_jobs", 1))
                base_share = max(1, total_cores // concurrent)
                status = shared_pool.get_status()
                available = status.get("total_cores", total_cores) - status.get("allocated_cores", 0)
                if available > 0:
                    base_share = min(base_share, available)
                core_share = max(1, min(base_share, total_cores))
            except Exception:
                logging.debug(
                    "Failed to derive core share from global pool; falling back to CONTROL PAL.",
                    exc_info=True,
                )
        resolved_pal = max(1, core_share)

    pal_override = resolved_pal

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

        candidate_jobs: List[_IMAGCandidateJob] = []
        candidate_geom_override = None if optimize_candidates else ""

        for mode_index, freq in modes:
            try:
                structures = run_plotvib_mode(
                    current_hess_path,
                    mode_index,
                    workdir=imag_folder,
                    amplitude=abs(displacement_scale),
                )
            except Exception:
                logging.debug(f"Mode {mode_index}: failed to generate displacement geometries.", exc_info=True)
                continue

            for direction_label, vib_xyz_path in structures.items():
                candidate_geom = imag_folder / f"iter{iteration}_mode{mode_index:03d}_{direction_label}.xyz"
                try:
                    shutil.copy2(vib_xyz_path, candidate_geom)
                except OSError:
                    candidate_geom = vib_xyz_path

                sp_label = f"iter{iteration}_mode{mode_index:03d}_{direction_label}"
                sp_input_path = imag_folder / f"{sp_label}_sp.inp"
                sp_output_path = imag_folder / f"{sp_label}_sp.out"

                # RECALC mode: Check if this single-point calculation already exists and is valid
                # Do this BEFORE deleting files!
                sp_exists_and_valid = False
                if recalc and sp_output_path.exists():
                    try:
                        if sp_output_path.stat().st_size >= 100:
                            with sp_output_path.open("r", encoding="utf-8", errors="replace") as f:
                                content = f.read()
                                if OK_MARKER in content:
                                    energy = find_electronic_energy(str(sp_output_path))
                                    if energy is not None:
                                        sp_exists_and_valid = True
                                        logging.info(f"[recalc] IMAG SP mode {mode_index} ({direction_label}): reusing existing result (energy={energy})")
                    except Exception as e:
                        logging.debug(f"[recalc] Failed to check existing IMAG SP output for mode {mode_index}: {e}")

                # Only delete and recreate if job doesn't exist or is invalid
                if not sp_exists_and_valid:
                    for tmp in (sp_input_path, sp_output_path):
                        if tmp.exists():
                            try:
                                tmp.unlink()
                            except OSError:
                                pass

                    _write_imag_input(
                        candidate_geom,
                        sp_input_path,
                        include_freq=False,
                        additions_payload=additions_eff,
                        geom_override=candidate_geom_override,
                        pal_override_local=pal_override,
                        maxcore_override_local=maxcore_effective,
                    )

                candidate_jobs.append(
                    _IMAGCandidateJob(
                        label=sp_label,
                        mode_index=mode_index,
                        frequency=freq,
                        geometry_path=candidate_geom,
                        sp_input_path=sp_input_path,
                        sp_output_path=sp_output_path,
                        direction=direction_label,
                    )
                )

        if not candidate_jobs:
            logging.info("No imaginary-mode displacements available; stopping IMAG.")
            break

        sp_results = _execute_sp_candidates(
            candidate_jobs,
            pool=shared_pool,
            core_allocation=pal_override,
            maxcore_value=maxcore_effective,
        )

        best_candidate: ImagModeCandidate | None = None
        for result in sp_results:
            if not result.success or result.energy is None:
                continue
            if result.energy < best_energy - sp_energy_window:
                if best_candidate is None or result.energy < best_candidate.energy - sp_energy_window:
                    best_candidate = ImagModeCandidate(
                        mode_index=result.job.mode_index,
                        frequency=result.job.frequency,
                        energy=result.energy,
                        geometry_path=result.job.geometry_path,
                        sp_output=result.job.sp_output_path,
                        direction=result.job.direction,
                        optimized_geometry=result.optimized_geometry,
                    )

        if best_candidate is None:
            logging.info("No imaginary-mode displacement lowered the energy; stopping IMAG.")
            break

        iteration += 1
        freq_input_path = imag_folder / f"input_{iteration}.inp"
        freq_output_path = imag_folder / f"output_{iteration}.out"

        # RECALC mode: Check if frequency calculation was already completed successfully
        freq_calc_needed = True
        additions_current = additions_eff
        gbw_candidate = best_candidate.sp_output.with_suffix(".gbw")

        if recalc and freq_output_path.exists():
            try:
                if freq_output_path.stat().st_size >= 100:
                    with freq_output_path.open("r", encoding="utf-8", errors="replace") as f:
                        content = f.read()
                        if OK_MARKER in content:
                            # Check if this iteration resolved the imaginary frequency
                            check_freq = search_imaginary_mode2(str(freq_output_path))
                            if check_freq is None or check_freq >= threshold:
                                logging.info(f"[recalc] IMAG iteration {iteration} already complete and successful; using existing result.")
                                freq_calc_needed = False
                                current_log_path = freq_output_path
                                current_hess_path = freq_output_path.with_suffix(".hess")
                                if not current_hess_path.exists():
                                    current_hess_path = imag_folder / f"input_{iteration}.hess"
                                best_energy = find_electronic_energy(str(freq_output_path)) or best_energy
                                last_success_iteration = iteration
            except Exception as e:
                logging.debug(f"[recalc] Failed to check existing IMAG freq output for iteration {iteration}: {e}")

        if freq_calc_needed:
            # Clean up old files before creating new input
            for tmp in (freq_input_path, freq_output_path):
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except OSError:
                        pass

            freq_geometry_source = (
                best_candidate.optimized_geometry
                if best_candidate.optimized_geometry and Path(best_candidate.optimized_geometry).exists()
                else best_candidate.geometry_path
            )
            _write_imag_input(
                freq_geometry_source,
                freq_input_path,
                include_freq=True,
                additions_payload=additions_current,
                geom_override=None,
                pal_override_local=pal_override,
                maxcore_override_local=maxcore_effective,
            )
            if gbw_candidate.exists():
                _inject_moinp_block(freq_input_path, gbw_candidate)
            else:
                logging.warning(
                    f"Iteration {iteration}: Expected GBW '{gbw_candidate}' missing; proceeding without MOREAD."
                )
            success = run_orca_IMAG(str(freq_input_path), iteration, working_dir=imag_folder)
        else:
            # Skip running ORCA but mark as success since result already exists
            success = True
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

            if geometry_mismatch and additions_current and freq_calc_needed:
                logging.warning(
                    f"Iteration {iteration}: geometry mismatch detected; retrying without supplemental additions."
                )
                additions_current = ""
                additions_eff = ""
                additions_base = ""
                freq_geometry_source = (
                    best_candidate.optimized_geometry
                    if best_candidate.optimized_geometry and Path(best_candidate.optimized_geometry).exists()
                    else best_candidate.geometry_path
                )
                _write_imag_input(
                    freq_geometry_source,
                    freq_input_path,
                    include_freq=True,
                    additions_payload=additions_current,
                    geom_override=None,
                    pal_override_local=pal_override,
                    maxcore_override_local=maxcore_effective,
                )
                if gbw_candidate.exists():
                    _inject_moinp_block(freq_input_path, gbw_candidate)
                success = run_orca_IMAG(str(freq_input_path), iteration, working_dir=imag_folder)

            if not success:
                logging.warning(f"Iteration {iteration}: frequency calculation failed; stopping IMAG loop.")
                break

        # Update paths and energies based on completed (or skipped) calculation
        if not freq_calc_needed:
            # Already set in the recalc check above
            pass
        else:
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
    destination_folder = input_path.parent
    destination_log = input_path
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
            # Also propagate the IMAG-refined geometry back to OCCUPIER-propagated files
            propagated_name = f"input_{step_name}_OCCUPIER.xyz"
            propagated_path = destination_folder / propagated_name
            if propagated_path.exists():
                try:
                    shutil.copy2(destination_structure, propagated_path)
                    logging.info(
                        "[IMAG] Propagated refined geometry to %s",
                        propagated_path,
                    )
                except Exception as copy_exc:  # noqa: BLE001
                    logging.warning(
                        "[IMAG] Failed to propagate refined geometry to %s: %s",
                        propagated_path,
                        copy_exc,
                    )
        except Exception as e:
            logging.warning(f"copy back xyz: {e}")
    else:
        print(f"ERROR: Structure file 'input_{last_success_iteration}.xyz' not found.")
