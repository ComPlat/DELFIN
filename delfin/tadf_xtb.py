from __future__ import annotations

import argparse
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import threading
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path
from delfin.dynamic_pool import JobPriority, PoolJob
from delfin.global_manager import get_global_manager
from delfin.orca import run_orca
from delfin.qm_runtime import get_xtb4stda_runtime_status, resolve_tool, run_tool

logger = get_logger(__name__)

_STATE_LINE_RE = re.compile(
    r"^\s*(\d+)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)"
)
_LABEL_SAFE_RE = re.compile(r"[^A-Za-z0-9._-]+")
_XTB_TOTAL_ENERGY_RE = re.compile(
    r"TOTAL ENERGY\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+Eh"
)
_ALLOWED_FOSC_MIN = 1.0e-3
_EV_PER_HARTREE = 27.211386245988
_HC_EV_NM = 1239.8419843320026


@dataclass(frozen=True)
class WorkflowEntry:
    label: str
    smiles: Optional[str] = None
    xyz_text: Optional[str] = None
    xyz_path: Optional[str] = None


@dataclass
class TadfXtbResult:
    label: str
    smiles: Optional[str]
    charge: int
    multiplicity: int
    xtb_method: str
    excited_method: str
    workdir: str
    initial_xyz: str
    s0_xyz: str
    crest_xyz: Optional[str]
    goat_xyz: Optional[str]
    crest_output: Optional[str]
    goat_output: Optional[str]
    s0_opt_output: Optional[str]
    s0_xtb_energy_eh: Optional[float]
    s0_sp_output: Optional[str]
    t1_opt_xyz: Optional[str]
    t1_vertical_xtb_output: Optional[str]
    t1_opt_output: Optional[str]
    t1_xtb_energy_eh: Optional[float]
    t1_vertical_xtb_ev: Optional[float]
    t1_adiabatic_ev: Optional[float]
    t1_relaxation_ev: Optional[float]
    stokes_shift_est_ev: Optional[float]
    s1_relaxed_est_ev: Optional[float]
    t1_relaxed_ev: Optional[float]
    singlet_output: str
    triplet_output: str
    tda_singlet: Optional[str]
    tda_triplet: Optional[str]
    s1_state: Optional[int]
    s1_ev: Optional[float]
    s1_nm: Optional[float]
    lambda_abs_nm: Optional[float]
    lambda_em_s1_nm: Optional[float]
    lambda_pl_est_nm: Optional[float]
    s1_f_osc: Optional[float]
    first_allowed_singlet_state: Optional[int]
    first_allowed_singlet_ev: Optional[float]
    first_allowed_singlet_nm: Optional[float]
    first_allowed_singlet_f_osc: Optional[float]
    brightest_singlet_state: Optional[int]
    brightest_singlet_ev: Optional[float]
    brightest_singlet_nm: Optional[float]
    brightest_singlet_f_osc: Optional[float]
    t1_ev: Optional[float]
    t1_nm: Optional[float]
    delta_est_ev: Optional[float]
    delta_est_relaxed_est_ev: Optional[float]


@dataclass(frozen=True)
class PreflightCheck:
    name: str
    ok: bool
    detail: str


def _safe_label(label: str, fallback: str) -> str:
    text = _LABEL_SAFE_RE.sub("_", (label or "").strip()).strip("._-")
    return text or fallback


def _rdkit_available() -> bool:
    return importlib.util.find_spec("rdkit") is not None


def _collect_preflight_checks(
    *,
    need_rdkit: bool,
    need_xtb: bool,
    need_crest: bool,
    need_orca: bool,
    excited_method: str,
) -> list[PreflightCheck]:
    checks: list[PreflightCheck] = []

    if need_rdkit:
        if _rdkit_available():
            checks.append(PreflightCheck("rdkit", True, sys.executable))
        else:
            checks.append(
                PreflightCheck(
                    "rdkit",
                    False,
                    f"missing in {sys.executable}; install project Python dependencies",
                )
            )

    tool_names = ["xtb4stda", excited_method]
    if need_xtb:
        tool_names.insert(0, "xtb")
    if need_crest:
        tool_names.append("crest")
    if need_orca:
        tool_names.insert(0, "orca")

    for tool_name in tool_names:
        resolved = resolve_tool(tool_name)
        if resolved is None:
            checks.append(PreflightCheck(tool_name, False, "not found"))
            continue
        checks.append(
            PreflightCheck(
                tool_name,
                True,
                f"{resolved.path} [{resolved.source}]",
            )
        )

    runtime_root, missing_runtime = get_xtb4stda_runtime_status()
    if missing_runtime:
        checks.append(
            PreflightCheck(
                "xtb4stda-data",
                False,
                f"{runtime_root} missing {', '.join(missing_runtime)}",
            )
        )
    else:
        checks.append(PreflightCheck("xtb4stda-data", True, str(runtime_root)))

    return checks


def _print_preflight_checks(checks: list[PreflightCheck]) -> None:
    for check in checks:
        status = "OK" if check.ok else "MISSING"
        print(f"{check.name}: {status} - {check.detail}")


def _write_smiles_inputs(smiles: str, label: str, workdir: Path) -> tuple[Path, Path]:
    from delfin.smiles_converter import smiles_to_xyz_quick

    xyz_body, err = smiles_to_xyz_quick(smiles)
    if err or not xyz_body:
        raise RuntimeError(err or "SMILES to 3D conversion failed")

    coord_lines = [line.rstrip() for line in xyz_body.splitlines() if line.strip()]
    if not coord_lines:
        raise RuntimeError("SMILES conversion returned no coordinates")

    start_path = workdir / "start.txt"
    start_path.write_text("\n".join(coord_lines) + "\n", encoding="utf-8")

    xyz_path = workdir / f"{label}.xyz"
    xyz_path.write_text(
        f"{len(coord_lines)}\n{label}\n" + "\n".join(coord_lines) + "\n",
        encoding="utf-8",
    )
    return start_path, xyz_path


def _write_xyz_inputs(input_xyz: str, label: str, workdir: Path) -> tuple[Path, Path]:
    source_xyz = resolve_path(input_xyz)
    if not source_xyz.is_file():
        raise RuntimeError(f"XYZ file not found: {source_xyz}")

    coord_block = _xyz_coord_block(source_xyz)
    coord_lines = [line.rstrip() for line in coord_block.splitlines() if line.strip()]
    if not coord_lines:
        raise RuntimeError(f"XYZ file contains no coordinates: {source_xyz}")

    start_path = workdir / "start.txt"
    start_path.write_text("\n".join(coord_lines) + "\n", encoding="utf-8")

    xyz_path = workdir / f"{label}.xyz"
    source_lines = source_xyz.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(source_lines) < 2:
        raise RuntimeError(f"XYZ file is incomplete: {source_xyz}")
    xyz_path.write_text("\n".join(source_lines).rstrip() + "\n", encoding="utf-8")
    return start_path, xyz_path


def _write_xyz_text_inputs(xyz_text: str, label: str, workdir: Path) -> tuple[Path, Path]:
    raw_lines = [line.rstrip() for line in str(xyz_text).splitlines()]
    content_lines = [line for line in raw_lines if line.strip()]
    if not content_lines:
        raise RuntimeError("XYZ input text is empty")

    coord_lines: list[str]
    xyz_lines: list[str]
    try:
        atom_count = int(content_lines[0].strip())
    except ValueError:
        coord_lines = content_lines
        xyz_lines = [str(len(coord_lines)), label, *coord_lines]
    else:
        if atom_count < 1:
            raise RuntimeError("XYZ atom count must be positive")
        if len(content_lines) < atom_count + 2:
            raise RuntimeError("XYZ input text is incomplete")
        xyz_lines = content_lines[: atom_count + 2]
        coord_lines = xyz_lines[2:]

    start_path = workdir / "start.txt"
    start_path.write_text("\n".join(coord_lines) + "\n", encoding="utf-8")

    xyz_path = workdir / f"{label}.xyz"
    xyz_path.write_text("\n".join(xyz_lines).rstrip() + "\n", encoding="utf-8")
    return start_path, xyz_path


def _xyz_coord_block(xyz_path: Path) -> str:
    lines = xyz_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 3:
        raise RuntimeError(f"XYZ file is incomplete: {xyz_path}")
    coord_lines = [line.rstrip() for line in lines[2:] if line.strip()]
    if not coord_lines:
        raise RuntimeError(f"XYZ file contains no coordinates: {xyz_path}")
    return "\n".join(coord_lines) + "\n"


_XTB_MAX_THREADS = 8  # xTB scales poorly beyond ~8 threads; higher values trigger Fortran I/O bugs


def _tool_env(cores: int) -> dict[str, str]:
    return {
        "OMP_NUM_THREADS": str(min(max(1, int(cores)), _XTB_MAX_THREADS)),
        "MKL_NUM_THREADS": "1",
    }


def _resolve_global_limits(requested_cores: int, requested_maxcore: int) -> tuple[int, int]:
    manager = get_global_manager()
    try:
        if manager.is_initialized():
            return manager.resolve_job_resources(
                requested_cores=requested_cores,
                requested_maxcore=requested_maxcore,
            )
    except Exception:
        logger.debug("Could not resolve global limits for tadf_xtb", exc_info=True)
    return max(1, int(requested_cores)), max(256, int(requested_maxcore))


def _parse_xtb_total_energy(output_path: Path) -> Optional[float]:
    matches = _XTB_TOTAL_ENERGY_RE.findall(
        output_path.read_text(encoding="utf-8", errors="ignore")
    )
    if not matches:
        return None
    try:
        return float(matches[-1])
    except ValueError:
        return None


def _run_xtb_opt(
    input_xyz: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    cores: int,
    stage_name: str,
    output_name: str,
) -> tuple[Path, Path, Optional[float]]:
    stage_dir = workdir / stage_name
    stage_dir.mkdir(parents=True, exist_ok=True)

    local_input = stage_dir / input_xyz.name
    shutil.copyfile(input_xyz, local_input)

    output_path = stage_dir / "xtb_opt.out"
    uhf = max(0, int(multiplicity) - 1)
    with output_path.open("w", encoding="utf-8") as log_file:
        run_tool(
            "xtb",
            [
                local_input.name,
                "--opt",
                "--chrg",
                str(charge),
                "--uhf",
                str(uhf),
            ],
            cwd=stage_dir,
            env=_tool_env(cores),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
            track_process=True,
        )

    optimized_xyz = stage_dir / "xtbopt.xyz"
    if not optimized_xyz.exists():
        raise RuntimeError(f"xTB optimization did not produce {optimized_xyz.name}")

    final_xyz = workdir / output_name
    shutil.copyfile(optimized_xyz, final_xyz)
    return output_path, final_xyz, _parse_xtb_total_energy(output_path)


def _run_xtb_singlepoint(
    input_xyz: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    cores: int,
    stage_name: str,
) -> tuple[Path, Optional[float]]:
    stage_dir = workdir / stage_name
    stage_dir.mkdir(parents=True, exist_ok=True)

    local_input = stage_dir / input_xyz.name
    shutil.copyfile(input_xyz, local_input)

    output_path = stage_dir / "xtb_sp.out"
    uhf = max(0, int(multiplicity) - 1)
    with output_path.open("w", encoding="utf-8") as log_file:
        run_tool(
            "xtb",
            [
                local_input.name,
                "--chrg",
                str(charge),
                "--uhf",
                str(uhf),
            ],
            cwd=stage_dir,
            env=_tool_env(cores),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
            track_process=True,
        )

    return output_path, _parse_xtb_total_energy(output_path)


def _run_crest(
    input_xyz: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    cores: int,
) -> tuple[Path, Path]:
    stage_dir = workdir / "crest"
    stage_dir.mkdir(parents=True, exist_ok=True)

    local_input = stage_dir / "input.xyz"
    shutil.copyfile(input_xyz, local_input)

    output_path = stage_dir / "CREST.out"
    uhf = max(0, int(multiplicity) - 1)
    with output_path.open("w", encoding="utf-8") as log_file:
        run_tool(
            "crest",
            [
                local_input.name,
                "--chrg",
                str(charge),
                "--uhf",
                str(uhf),
            ],
            cwd=stage_dir,
            env=_tool_env(cores),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
            track_process=True,
        )

    crest_best = stage_dir / "crest_best.xyz"
    if not crest_best.exists():
        raise RuntimeError(f"CREST did not produce {crest_best.name}")

    final_xyz = workdir / "crest_best.xyz"
    shutil.copyfile(crest_best, final_xyz)
    return output_path, final_xyz


def _write_goat_input(
    input_xyz: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    xtb_method: str,
    cores: int,
    maxcore: int,
) -> Path:
    coord_block = _xyz_coord_block(input_xyz)
    input_path = workdir / "XTB_GOAT.inp"
    input_path.write_text(
        f"!{xtb_method} GOAT\n"
        f"%pal nprocs {max(1, int(cores))} end\n"
        f"%maxcore {max(1, int(maxcore))}\n"
        f"*xyz {charge} {multiplicity}\n"
        f"{coord_block}"
        "*\n",
        encoding="utf-8",
    )
    return input_path


def _run_goat(
    input_xyz: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    xtb_method: str,
    cores: int,
    maxcore: int,
) -> tuple[Path, Path]:
    stage_dir = workdir / "goat"
    stage_dir.mkdir(parents=True, exist_ok=True)

    local_input = stage_dir / input_xyz.name
    shutil.copyfile(input_xyz, local_input)

    input_path = _write_goat_input(
        local_input,
        stage_dir,
        charge=charge,
        multiplicity=multiplicity,
        xtb_method=xtb_method,
        cores=cores,
        maxcore=maxcore,
    )
    output_path = stage_dir / "output_XTB_GOAT.out"
    ok = run_orca(str(input_path), str(output_path), working_dir=stage_dir)
    if not ok:
        raise RuntimeError(f"GOAT failed; see {output_path}")
    goat_xyz = stage_dir / "XTB_GOAT.globalminimum.xyz"
    if not goat_xyz.exists():
        raise RuntimeError(f"GOAT did not produce {goat_xyz.name}")
    final_xyz = workdir / "goat_globalminimum.xyz"
    shutil.copyfile(goat_xyz, final_xyz)
    return output_path, final_xyz


def _run_xtb4stda(
    geometry_xyz: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    cores: int,
) -> Path:
    uhf = max(0, int(multiplicity) - 1)
    output_path = workdir / "xtb4stda.out"
    with output_path.open("w", encoding="utf-8") as log_file:
        run_tool(
            "xtb4stda",
            [geometry_xyz.name, "-chrg", str(charge), "-uhf", str(uhf)],
            cwd=workdir,
            env=_tool_env(cores),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
            track_process=True,
        )
    wfn_path = workdir / "wfn.xtb"
    if not wfn_path.exists():
        raise RuntimeError(f"xtb4stda did not produce {wfn_path.name}")
    return output_path


def _run_excited_method(
    workdir: Path,
    *,
    method: str,
    triplet: bool,
    energy_window: float,
    cores: int,
) -> tuple[Path, Optional[Path]]:
    output_stem = f"{method}_{'triplet' if triplet else 'singlet'}"
    output_path = workdir / f"{output_stem}.out"
    args = ["-xtb", "-e", str(energy_window)]
    if triplet:
        args.insert(1, "-t")
    with output_path.open("w", encoding="utf-8") as log_file:
        run_tool(
            method,
            args,
            cwd=workdir,
            env=_tool_env(cores),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
            track_process=True,
        )

    tda_path = workdir / "tda.dat"
    renamed = workdir / ("tda_triplet.dat" if triplet else "tda_singlet.dat")
    if tda_path.exists():
        tda_path.replace(renamed)
        return output_path, renamed
    return output_path, None


def _parse_stda_states(output_path: Path) -> list[dict[str, float | int]]:
    states: list[dict[str, float | int]] = []
    in_table = False
    for raw_line in output_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.rstrip()
        if "state    eV" in line and "fL" in line:
            in_table = True
            continue
        if not in_table:
            continue
        if not line.strip():
            if states:
                break
            continue
        match = _STATE_LINE_RE.match(line)
        if not match:
            if states:
                break
            continue
        states.append(
            {
                "state": int(match.group(1)),
                "ev": float(match.group(2)),
                "nm": float(match.group(3)),
                "f_osc": float(match.group(4)),
                "rv_corr": float(match.group(5)),
            }
        )
    return states


def _first_allowed_singlet(
    states: list[dict[str, float | int]],
    *,
    min_f_osc: float = _ALLOWED_FOSC_MIN,
) -> Optional[dict[str, float | int]]:
    for state in states:
        if float(state["f_osc"]) >= min_f_osc:
            return state
    return None


def _brightest_singlet(states: list[dict[str, float | int]]) -> Optional[dict[str, float | int]]:
    if not states:
        return None
    return max(states, key=lambda state: float(state["f_osc"]))


def _energy_to_wavelength_nm(energy_ev: Optional[float]) -> Optional[float]:
    if energy_ev is None or energy_ev <= 0.0:
        return None
    return _HC_EV_NM / energy_ev


def run_single_tadf_xtb(
    entry: WorkflowEntry,
    *,
    charge: int,
    multiplicity: int,
    xtb_method: str,
    excited_method: str,
    energy_window: float,
    cores: int,
    maxcore: int,
    workdir: Path,
    use_crest: bool,
    use_goat: bool,
    run_t1_opt: bool,
    t1_multiplicity: int,
    optimize_s0: bool = True,
) -> TadfXtbResult:
    workdir.mkdir(parents=True, exist_ok=True)
    if entry.xyz_path:
        _start_path, initial_xyz = _write_xyz_inputs(entry.xyz_path, entry.label, workdir)
    elif entry.xyz_text:
        _start_path, initial_xyz = _write_xyz_text_inputs(entry.xyz_text, entry.label, workdir)
    elif entry.smiles:
        _start_path, initial_xyz = _write_smiles_inputs(entry.smiles, entry.label, workdir)
    else:
        raise RuntimeError("workflow entry has neither SMILES nor XYZ input")

    current_xyz = initial_xyz
    crest_xyz: Optional[Path] = None
    goat_xyz: Optional[Path] = None
    crest_output: Optional[Path] = None
    goat_output: Optional[Path] = None
    s0_opt_output: Optional[Path] = None
    s0_xtb_energy_eh: Optional[float] = None
    s0_sp_output: Optional[Path] = None
    t1_opt_xyz: Optional[Path] = None
    t1_vertical_xtb_output: Optional[Path] = None
    t1_opt_output: Optional[Path] = None
    t1_xtb_energy_eh: Optional[float] = None
    t1_vertical_xtb_ev: Optional[float] = None
    t1_adiabatic_ev: Optional[float] = None
    t1_relaxation_ev: Optional[float] = None
    stokes_shift_est_ev: Optional[float] = None
    s1_relaxed_est_ev: Optional[float] = None
    lambda_pl_est_nm: Optional[float] = None
    t1_relaxed_ev: Optional[float] = None
    delta_est_relaxed_est_ev: Optional[float] = None

    if use_crest:
        crest_output, crest_xyz = _run_crest(
            current_xyz,
            workdir,
            charge=charge,
            multiplicity=multiplicity,
            cores=cores,
        )
        current_xyz = crest_xyz

    if use_goat:
        goat_output, goat_xyz = _run_goat(
            current_xyz,
            workdir,
            charge=charge,
            multiplicity=multiplicity,
            xtb_method=xtb_method,
            cores=cores,
            maxcore=maxcore,
        )
        current_xyz = goat_xyz
    elif optimize_s0:
        s0_opt_output, current_xyz, s0_xtb_energy_eh = _run_xtb_opt(
            current_xyz,
            workdir,
            charge=charge,
            multiplicity=multiplicity,
            cores=cores,
            stage_name="s0_xtb_opt",
            output_name="s0_xtbopt.xyz",
        )

    _run_xtb4stda(
        current_xyz,
        workdir,
        charge=charge,
        multiplicity=multiplicity,
        cores=cores,
    )
    singlet_out, tda_singlet = _run_excited_method(
        workdir,
        method=excited_method,
        triplet=False,
        energy_window=energy_window,
        cores=cores,
    )
    triplet_out, tda_triplet = _run_excited_method(
        workdir,
        method=excited_method,
        triplet=True,
        energy_window=energy_window,
        cores=cores,
    )

    singlets = _parse_stda_states(singlet_out)
    triplets = _parse_stda_states(triplet_out)
    s1 = singlets[0] if singlets else None
    first_allowed = _first_allowed_singlet(singlets)
    brightest = _brightest_singlet(singlets)
    t1 = triplets[0] if triplets else None

    if run_t1_opt:
        if s0_xtb_energy_eh is None:
            s0_sp_output, s0_xtb_energy_eh = _run_xtb_singlepoint(
                current_xyz,
                workdir,
                charge=charge,
                multiplicity=multiplicity,
                cores=cores,
                stage_name="s0_xtb_sp",
            )
        t1_vertical_xtb_output, t1_vertical_xtb_energy_eh = _run_xtb_singlepoint(
            current_xyz,
            workdir,
            charge=charge,
            multiplicity=t1_multiplicity,
            cores=cores,
            stage_name="t1_xtb_vertical_sp",
        )
        if s0_xtb_energy_eh is not None and t1_vertical_xtb_energy_eh is not None:
            t1_vertical_xtb_ev = (t1_vertical_xtb_energy_eh - s0_xtb_energy_eh) * _EV_PER_HARTREE
        t1_opt_output, t1_opt_xyz, t1_xtb_energy_eh = _run_xtb_opt(
            current_xyz,
            workdir,
            charge=charge,
            multiplicity=t1_multiplicity,
            cores=cores,
            stage_name="t1_xtb_opt",
            output_name="t1_xtbopt.xyz",
        )
        if s0_xtb_energy_eh is not None and t1_xtb_energy_eh is not None:
            t1_adiabatic_ev = (t1_xtb_energy_eh - s0_xtb_energy_eh) * _EV_PER_HARTREE
            t1_relaxed_ev = t1_adiabatic_ev

    if t1_vertical_xtb_ev is not None and t1_adiabatic_ev is not None:
        t1_relaxation_ev = t1_vertical_xtb_ev - t1_adiabatic_ev
        stokes_shift_est_ev = t1_relaxation_ev

    if s1 is not None and stokes_shift_est_ev is not None:
        s1_relaxed_est_ev = float(s1["ev"]) - stokes_shift_est_ev
        lambda_pl_est_nm = _energy_to_wavelength_nm(s1_relaxed_est_ev)

    if s1_relaxed_est_ev is not None and t1_adiabatic_ev is not None:
        delta_est_relaxed_est_ev = s1_relaxed_est_ev - t1_adiabatic_ev

    return TadfXtbResult(
        label=entry.label,
        smiles=entry.smiles,
        charge=charge,
        multiplicity=multiplicity,
        xtb_method=xtb_method,
        excited_method=excited_method,
        workdir=str(workdir),
        initial_xyz=str(initial_xyz),
        s0_xyz=str(current_xyz),
        crest_xyz=str(crest_xyz) if crest_xyz else None,
        goat_xyz=str(goat_xyz) if goat_xyz else None,
        crest_output=str(crest_output) if crest_output else None,
        goat_output=str(goat_output) if goat_output else None,
        s0_opt_output=str(s0_opt_output) if s0_opt_output else None,
        s0_xtb_energy_eh=s0_xtb_energy_eh,
        s0_sp_output=str(s0_sp_output) if s0_sp_output else None,
        t1_opt_xyz=str(t1_opt_xyz) if t1_opt_xyz else None,
        t1_vertical_xtb_output=str(t1_vertical_xtb_output) if t1_vertical_xtb_output else None,
        t1_opt_output=str(t1_opt_output) if t1_opt_output else None,
        t1_xtb_energy_eh=t1_xtb_energy_eh,
        t1_vertical_xtb_ev=t1_vertical_xtb_ev,
        t1_adiabatic_ev=t1_adiabatic_ev,
        t1_relaxation_ev=t1_relaxation_ev,
        stokes_shift_est_ev=stokes_shift_est_ev,
        s1_relaxed_est_ev=s1_relaxed_est_ev,
        t1_relaxed_ev=t1_relaxed_ev,
        singlet_output=str(singlet_out),
        triplet_output=str(triplet_out),
        tda_singlet=str(tda_singlet) if tda_singlet else None,
        tda_triplet=str(tda_triplet) if tda_triplet else None,
        s1_state=int(s1["state"]) if s1 else None,
        s1_ev=float(s1["ev"]) if s1 else None,
        s1_nm=float(s1["nm"]) if s1 else None,
        lambda_abs_nm=float(s1["nm"]) if s1 else None,
        lambda_em_s1_nm=lambda_pl_est_nm,
        lambda_pl_est_nm=lambda_pl_est_nm,
        s1_f_osc=float(s1["f_osc"]) if s1 else None,
        first_allowed_singlet_state=int(first_allowed["state"]) if first_allowed else None,
        first_allowed_singlet_ev=float(first_allowed["ev"]) if first_allowed else None,
        first_allowed_singlet_nm=float(first_allowed["nm"]) if first_allowed else None,
        first_allowed_singlet_f_osc=float(first_allowed["f_osc"]) if first_allowed else None,
        brightest_singlet_state=int(brightest["state"]) if brightest else None,
        brightest_singlet_ev=float(brightest["ev"]) if brightest else None,
        brightest_singlet_nm=float(brightest["nm"]) if brightest else None,
        brightest_singlet_f_osc=float(brightest["f_osc"]) if brightest else None,
        t1_ev=float(t1["ev"]) if t1 else None,
        t1_nm=float(t1["nm"]) if t1 else None,
        delta_est_ev=(float(s1["ev"]) - float(t1["ev"])) if s1 and t1 else None,
        delta_est_relaxed_est_ev=delta_est_relaxed_est_ev,
    )


def _load_entries(
    smiles_values: list[str],
    smiles_file: Optional[str],
    xyz_values: list[str],
    xyz_file_values: list[str],
    label: Optional[str],
) -> list[WorkflowEntry]:
    entries: list[WorkflowEntry] = []
    for idx, smiles in enumerate(smiles_values, start=1):
        fallback = f"mol_{idx:03d}"
        entries.append(WorkflowEntry(label=_safe_label(label or fallback, fallback), smiles=smiles.strip()))

    single_xyz_label = (
        len(xyz_values) + len(xyz_file_values) == 1
        and not smiles_values
        and not smiles_file
    )
    for idx, xyz_value in enumerate(xyz_values, start=1):
        fallback = f"xyz_{idx:03d}"
        requested = label if single_xyz_label else fallback
        entries.append(
            WorkflowEntry(
                label=_safe_label(requested, fallback),
                xyz_text=xyz_value,
            )
        )

    for idx, xyz_value in enumerate(xyz_file_values, start=1):
        source_path = resolve_path(xyz_value)
        fallback = _safe_label(source_path.stem, f"xyz_{idx:03d}")
        requested = label if single_xyz_label else fallback
        entries.append(
            WorkflowEntry(
                label=_safe_label(requested, fallback),
                xyz_path=str(source_path),
            )
        )

    if smiles_file:
        path = resolve_path(smiles_file)
        line_index = 0
        for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            line_index += 1
            if "\t" in line:
                raw_label, raw_smiles = line.split("\t", 1)
            else:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    raw_label, raw_smiles = parts
                else:
                    raw_label, raw_smiles = f"mol_{line_index:03d}", parts[0]
            fallback = f"mol_{line_index:03d}"
            entries.append(WorkflowEntry(label=_safe_label(raw_label, fallback), smiles=raw_smiles.strip()))

    unique_entries: list[WorkflowEntry] = []
    seen: dict[str, int] = {}
    for entry in entries:
        count = seen.get(entry.label, 0) + 1
        seen[entry.label] = count
        if count == 1:
            unique_entries.append(entry)
            continue
        unique_entries.append(
            WorkflowEntry(
                label=f"{entry.label}_{count:02d}",
                smiles=entry.smiles,
                xyz_text=entry.xyz_text,
                xyz_path=entry.xyz_path,
            )
        )
    return unique_entries


def _print_result(result: TadfXtbResult) -> None:
    def _fmt(value: Optional[float], digits: int) -> str:
        if value is None:
            return "n/a"
        return f"{value:.{digits}f}"

    parts = [
        f"method={result.excited_method}",
        f"S1(v)={_fmt(result.s1_ev, 3)} eV ({_fmt(result.s1_nm, 1)} nm, f={_fmt(result.s1_f_osc, 4)})",
        f"T1={_fmt(result.t1_ev, 3)} eV ({_fmt(result.t1_nm, 1)} nm)",
        f"lambda_abs(S0->S1)={_fmt(result.lambda_abs_nm, 1)} nm",
        f"lambda_PL_est={_fmt(result.lambda_pl_est_nm, 1)} nm",
        f"DeltaEST_vert={_fmt(result.delta_est_ev, 3)} eV",
    ]
    if result.first_allowed_singlet_ev is not None:
        parts.append(
            "first_allowed="
            f"S{result.first_allowed_singlet_state} "
            f"{_fmt(result.first_allowed_singlet_ev, 3)} eV "
            f"({_fmt(result.first_allowed_singlet_nm, 1)} nm, "
            f"f={_fmt(result.first_allowed_singlet_f_osc, 4)})"
        )
    if result.brightest_singlet_ev is not None:
        parts.append(
            "brightest="
            f"S{result.brightest_singlet_state} "
            f"{_fmt(result.brightest_singlet_ev, 3)} eV "
            f"({_fmt(result.brightest_singlet_nm, 1)} nm, "
            f"f={_fmt(result.brightest_singlet_f_osc, 4)})"
        )
    if result.t1_adiabatic_ev is not None:
        parts.append(f"E_T1_adiab={_fmt(result.t1_adiabatic_ev, 3)} eV")
    if result.t1_vertical_xtb_ev is not None:
        parts.append(f"E_T1_vert_xTB={_fmt(result.t1_vertical_xtb_ev, 3)} eV")
    if result.t1_relaxation_ev is not None:
        parts.append(f"DeltaE_relax(T1)={_fmt(result.t1_relaxation_ev, 3)} eV")
    if result.stokes_shift_est_ev is not None:
        parts.append(f"Stokes_est={_fmt(result.stokes_shift_est_ev, 3)} eV")
    if result.s1_relaxed_est_ev is not None:
        parts.append(f"E_S1_rel_est={_fmt(result.s1_relaxed_est_ev, 3)} eV")
    if result.delta_est_relaxed_est_ev is not None:
        parts.append(f"DeltaEST_rel_est={_fmt(result.delta_est_relaxed_est_ev, 3)} eV")
    parts.append(f"workdir={result.workdir}")
    print(f"{result.label}: " + ", ".join(parts))


def run_cli(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin tadf_xtb",
        description=(
            "SMILES/XYZ -> optional CREST -> GOAT or xTB S0 optimization "
            "-> xtb4stda -> stda/std2(singlet/triplet) TADF screening workflow."
        ),
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run dependency preflight for the workflow and exit.",
    )
    parser.add_argument("--smiles", action="append", default=[], help="SMILES string. Repeatable.")
    parser.add_argument("--smiles-file", help="Optional file with one SMILES per line, optionally prefixed by label.")
    parser.add_argument("--xyz", action="append", default=[], help="Inline XYZ text. Repeatable.")
    parser.add_argument("--xyz-file", action="append", default=[], help="Input XYZ file path. Repeatable.")
    parser.add_argument("--label", help="Optional label for a single direct input (--smiles, --xyz, or --xyz-file).")
    parser.add_argument("--workdir", default="tadf_xtb_runs", help="Base output directory.")
    parser.add_argument("--charge", type=int, default=0, help="Molecular charge.")
    parser.add_argument("--multiplicity", type=int, default=1, help="Ground-state multiplicity.")
    parser.add_argument("--xtb-method", default="XTB2", help="ORCA xTB keyword for optional GOAT, e.g. XTB2.")
    parser.add_argument(
        "--excited-method",
        choices=("stda", "std2"),
        default="stda",
        help="Excited-state engine: stda = sTDA-xTB, std2 = sTD-DFT-xTB.",
    )
    parser.add_argument("--crest", action="store_true", help="Run CREST conformer search before the excited-state step.")
    parser.add_argument("--goat", action="store_true", help="Use ORCA GOAT instead of plain xTB S0 optimization.")
    parser.add_argument("--t1-opt", action="store_true", help="Run an additional xTB triplet optimization after stda.")
    parser.add_argument("--t1-multiplicity", type=int, default=3, help="Multiplicity for optional T1 xTB optimization.")
    parser.add_argument("--energy-window", type=float, default=10.0, help="Excitation window in eV for stda.")
    parser.add_argument("--pal", type=int, default=4, help="Cores per molecule workflow.")
    parser.add_argument("--parallel-jobs", type=int, default=1, help="Number of molecule workflows to run in parallel.")
    parser.add_argument("--maxcore", type=int, default=1000, help="Memory per core in MB for pool scheduling.")
    parser.add_argument("--json-out", help="Optional JSON summary path.")
    args = parser.parse_args(argv)

    if args.multiplicity < 1:
        parser.error("--multiplicity must be >= 1")
    if args.t1_multiplicity < 1:
        parser.error("--t1-multiplicity must be >= 1")

    # Load user settings and apply runtime environment (qm_tools PATH, ORCA, etc.)
    try:
        from delfin.user_settings import load_settings
        from delfin.runtime_setup import apply_runtime_environment
        _settings = load_settings()
        _runtime = _settings.get("runtime", {}) or {}
        apply_runtime_environment(
            qm_tools_root=_runtime.get("qm_tools_root", ""),
            orca_base=_runtime.get("orca_base", ""),
            csp_tools_root=_runtime.get("csp_tools_root", ""),
            tool_binaries=_runtime.get("tool_binaries", {}) or {},
        )
    except Exception:
        pass  # non-fatal: tools may still be in PATH

    preflight_checks = _collect_preflight_checks(
        need_rdkit=bool(args.smiles or args.smiles_file) or not bool(args.xyz or args.xyz_file),
        need_xtb=(not args.goat) or args.t1_opt or args.crest,
        need_crest=args.crest,
        need_orca=args.goat,
        excited_method=args.excited_method,
    )
    if args.check:
        _print_preflight_checks(preflight_checks)
        return 0 if all(check.ok for check in preflight_checks) else 1

    failed_checks = [check for check in preflight_checks if not check.ok]
    if failed_checks:
        _print_preflight_checks(failed_checks)
        print("tadf_xtb preflight failed. Run `delfin tadf_xtb --check` after fixing the missing items.")
        return 1

    entries = _load_entries(args.smiles, args.smiles_file, args.xyz, args.xyz_file, args.label)
    if not entries:
        parser.error("provide --smiles, --smiles-file, --xyz, and/or --xyz-file")

    base_workdir = resolve_path(args.workdir)
    base_workdir.mkdir(parents=True, exist_ok=True)

    results: list[TadfXtbResult] = []
    errors: dict[str, str] = {}

    if len(entries) == 1 or args.parallel_jobs <= 1:
        for entry in entries:
            try:
                cores_limit, maxcore_limit = _resolve_global_limits(args.pal, args.maxcore)
                result = run_single_tadf_xtb(
                    entry,
                    charge=args.charge,
                    multiplicity=args.multiplicity,
                    xtb_method=args.xtb_method,
                    excited_method=args.excited_method,
                    energy_window=args.energy_window,
                    cores=cores_limit,
                    maxcore=maxcore_limit,
                    workdir=base_workdir / entry.label,
                    use_crest=args.crest,
                    use_goat=args.goat,
                    run_t1_opt=args.t1_opt,
                    t1_multiplicity=args.t1_multiplicity,
                )
                results.append(result)
                _print_result(result)
            except Exception as exc:  # noqa: BLE001
                errors[entry.label] = str(exc)
                logger.error("[%s] tadf_xtb failed: %s", entry.label, exc, exc_info=True)
    else:
        manager = get_global_manager()
        manager.ensure_initialized(
            {
                "PAL": args.pal,
                "maxcore": args.maxcore,
                "pal_jobs": args.parallel_jobs,
                "parallel_workflows": "enable",
            }
        )
        pool = manager.get_pool()
        lock = threading.Lock()
        order = {entry.label: idx for idx, entry in enumerate(entries)}

        for entry in entries:
            workdir = base_workdir / entry.label

            def runner(*_args, _entry=entry, _workdir=workdir, **kwargs) -> None:
                allocated = kwargs.get("cores", args.pal)
                try:
                    cores_limit, maxcore_limit = _resolve_global_limits(int(allocated), args.maxcore)
                    result = run_single_tadf_xtb(
                        _entry,
                        charge=args.charge,
                        multiplicity=args.multiplicity,
                        xtb_method=args.xtb_method,
                        excited_method=args.excited_method,
                        energy_window=args.energy_window,
                        cores=cores_limit,
                        maxcore=maxcore_limit,
                        workdir=_workdir,
                        use_crest=args.crest,
                        use_goat=args.goat,
                        run_t1_opt=args.t1_opt,
                        t1_multiplicity=args.t1_multiplicity,
                    )
                    with lock:
                        results.append(result)
                except Exception as exc:  # noqa: BLE001
                    with lock:
                        errors[_entry.label] = str(exc)
                    raise

            pool_job = PoolJob(
                job_id=f"TADF_XTB_{entry.label}",
                cores_min=1,
                cores_optimal=max(1, int(args.pal)),
                cores_max=max(1, int(args.pal)),
                memory_mb=max(256, int(args.maxcore) * max(1, int(args.pal))),
                priority=JobPriority.NORMAL,
                execute_func=runner,
                args=(),
                kwargs={},
                estimated_duration=1800.0,
                working_dir=workdir,
            )
            pool.submit_job(pool_job)

        pool.wait_for_completion()
        results.sort(key=lambda item: order.get(item.label, 10**9))
        for result in results:
            _print_result(result)

    payload = {
        "results": [asdict(item) for item in results],
        "errors": errors,
    }
    if args.json_out:
        resolve_path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    if errors:
        print(json.dumps(payload, indent=2))
        return 1
    return 0
