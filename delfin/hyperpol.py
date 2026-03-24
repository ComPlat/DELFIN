from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import threading
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Optional

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path
from delfin.dynamic_pool import JobPriority, PoolJob
from delfin.global_manager import get_global_manager
from delfin.parser import calculate_beta_zzz_aligned
from delfin.qm_runtime import get_xtb4stda_runtime_status, resolve_tool, run_tool
from delfin.tadf_xtb import (
    WorkflowEntry,
    _load_entries,
    _rdkit_available,
    _run_crest,
    _run_xtb_opt,
    _run_xtb_singlepoint,
    _safe_label,
    _tool_env,
    _write_smiles_inputs,
    _write_xyz_inputs,
    _write_xyz_text_inputs,
)

logger = get_logger(__name__)

_HC_EV_NM = 1239.8419843320026
_AU_TO_ESU = 8.6393e-33
_AU_TO_ESU_30 = 8.6393e-3
_DEBYE_PER_AU = 2.541746473
_FLOAT_RE = r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?"
_HYPERPOL_HEADER_RE = re.compile(
    r"^\s*SHG first hyperpolarizability\s*\(\s*(.*?)\s*;\s*(.*?)\s*,\s*(.*?)\s*\)\s*$"
)
_TENSOR_ROW_RE = re.compile(
    rf"^\s*([xyz]{{2}})\.\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s*$"
)
_BETA_VEC_RE = re.compile(rf"^\s*beta_VEC_([XYZ])\s+({_FLOAT_RE})\s*$")
_BETA_SQ_RE = re.compile(
    rf"^\s*beta\*\*2_ZZZ\s+({_FLOAT_RE})\s+beta\*\*2_XZZ\s+({_FLOAT_RE})\s*$"
)
_BETA_HRS_RE = re.compile(rf"^\s*beta_HRS\s+({_FLOAT_RE})\s+DR\s+({_FLOAT_RE})\s*$")
_XTB_DIPOLE_RE = re.compile(
    rf"^\s*full:\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s*$",
    re.MULTILINE,
)


@dataclass(frozen=True)
class PreflightCheck:
    name: str
    ok: bool
    detail: str


@dataclass
class HyperpolResult:
    label: str
    smiles: Optional[str]
    charge: int
    multiplicity: int
    engine: str
    preopt: str
    energy_window_ev: float
    requested_wavelengths_nm: list[float]
    workdir: str
    initial_xyz: str
    selected_xyz: str
    preopt_output: Optional[str]
    dipole_output: str
    dipole_x_au: Optional[float]
    dipole_y_au: Optional[float]
    dipole_z_au: Optional[float]
    dipole_total_debye: Optional[float]
    xtb4stda_output: str
    response_output: str
    wavelength_file: str
    beta_hrs_file: str
    beta_tensor_file: str
    response_points: list[dict[str, Any]]


def _print_preflight_checks(checks: list[PreflightCheck]) -> None:
    for check in checks:
        status = "OK" if check.ok else "MISSING"
        print(f"{check.name}: {status} - {check.detail}")


def _collect_preflight_checks(
    *,
    need_rdkit: bool,
    preopt: str,
    engine: str,
) -> list[PreflightCheck]:
    checks: list[PreflightCheck] = []

    if need_rdkit:
        if _rdkit_available():
            checks.append(PreflightCheck("rdkit", True, "available"))
        else:
            checks.append(PreflightCheck("rdkit", False, "missing in current Python"))

    tool_names = ["xtb4stda", engine]
    if preopt in {"xtb", "crest"}:
        tool_names.insert(0, "xtb")
    if preopt == "crest":
        tool_names.append("crest")

    for tool_name in tool_names:
        resolved = resolve_tool(tool_name)
        if resolved is None:
            checks.append(PreflightCheck(tool_name, False, "not found"))
        else:
            checks.append(
                PreflightCheck(tool_name, True, f"{resolved.path} [{resolved.source}]")
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


def _write_inputs(entry: WorkflowEntry, workdir: Path) -> Path:
    if entry.xyz_path:
        _start_path, initial_xyz = _write_xyz_inputs(entry.xyz_path, entry.label, workdir)
    elif entry.xyz_text:
        _start_path, initial_xyz = _write_xyz_text_inputs(entry.xyz_text, entry.label, workdir)
    elif entry.smiles:
        _start_path, initial_xyz = _write_smiles_inputs(entry.smiles, entry.label, workdir)
    else:
        raise RuntimeError("workflow entry has neither SMILES nor XYZ input")
    return initial_xyz


def _write_wavelength_file(workdir: Path, wavelengths_nm: list[float]) -> Path:
    path = workdir / "wavelength"
    entries: list[str] = []
    for value in wavelengths_nm:
        if math.isinf(value):
            entries.append("inf")
        else:
            entries.append(f"{value:.8f}")
    contents = "\n".join(entries) + "\n"
    path.write_text(contents, encoding="utf-8")
    return path


def _parse_xtb_dipole(output_path: Path) -> tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    text = output_path.read_text(encoding="utf-8", errors="ignore")
    match = _XTB_DIPOLE_RE.search(text)
    if not match:
        return None, None, None, None
    dipole_x_debye = float(match.group(1))
    dipole_y_debye = float(match.group(2))
    dipole_z_debye = float(match.group(3))
    dipole_total_debye = float(match.group(4))
    return (
        dipole_x_debye / _DEBYE_PER_AU,
        dipole_y_debye / _DEBYE_PER_AU,
        dipole_z_debye / _DEBYE_PER_AU,
        dipole_total_debye,
    )


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
    if not (workdir / "wfn.xtb").exists():
        raise RuntimeError("xtb4stda did not produce wfn.xtb")
    return output_path


def _run_response(
    workdir: Path,
    *,
    engine: str,
    energy_window_ev: float,
    wavelength_count: int,
    cores: int,
) -> tuple[Path, Path, Path]:
    output_path = workdir / f"{engine}_hyperpol.out"
    with output_path.open("w", encoding="utf-8") as log_file:
        # std2 returns non-zero exit codes (e.g. 191) even on success;
        # rely on output file checks below instead of check=True.
        run_tool(
            engine,
            ["-xtb", "-e", str(energy_window_ev), "-resp", str(wavelength_count)],
            cwd=workdir,
            env=_tool_env(cores),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=False,
            track_process=True,
        )

    beta_hrs_path = workdir / "beta_HRS"
    beta_tensor_path = workdir / "beta_tensor"
    if not beta_hrs_path.exists():
        output_text = output_path.read_text(encoding="utf-8", errors="ignore")
        if "STOP no CSF" in output_text:
            raise RuntimeError(
                f"{engine} found no CSFs in the requested response window; increase --energy-window"
            )
        raise RuntimeError(f"{engine} did not produce beta_HRS")
    if not beta_tensor_path.exists():
        raise RuntimeError(f"{engine} did not produce beta_tensor")
    return output_path, beta_hrs_path, beta_tensor_path


def _parse_header_wavelength_nm(token: str) -> Optional[float]:
    text = str(token or "").strip()
    if not text:
        return None
    lowered = text.lower()
    if "inf" in lowered:
        return math.inf
    try:
        value = float(text)
    except ValueError:
        return None
    if value <= 0.0:
        return math.inf
    if value < 0.0:
        return None
    return value


def _parse_beta_hrs_file(beta_hrs_path: Path) -> list[dict[str, Any]]:
    points: list[dict[str, Any]] = []
    for raw_line in beta_hrs_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            energy_ev = float(parts[0])
            beta_hrs_au = float(parts[1])
        except ValueError:
            continue
        wavelength_nm = math.inf if abs(energy_ev) < 1.0e-12 else (_HC_EV_NM / energy_ev)
        points.append(
            {
                "energy_ev": energy_ev,
                "wavelength_nm": wavelength_nm,
                "is_static": abs(energy_ev) < 1.0e-12,
                "beta_hrs_au": beta_hrs_au,
                "beta_hrs_esu": beta_hrs_au * _AU_TO_ESU,
                "beta_hrs_esu_30": beta_hrs_au * _AU_TO_ESU_30,
            }
        )
    return points


def _parse_response_blocks(output_path: Path) -> list[dict[str, Any]]:
    lines = output_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    blocks: list[dict[str, Any]] = []
    idx = 0
    while idx < len(lines):
        header_match = _HYPERPOL_HEADER_RE.match(lines[idx])
        if not header_match:
            idx += 1
            continue

        block: dict[str, Any] = {
            "header": lines[idx].strip(),
            "response_wavelength_nm": _parse_header_wavelength_nm(header_match.group(2)),
            "tensor_au": {},
        }
        idx += 1
        while idx < len(lines):
            line = lines[idx]
            if _HYPERPOL_HEADER_RE.match(line):
                break

            tensor_match = _TENSOR_ROW_RE.match(line)
            if tensor_match:
                row = tensor_match.group(1)
                values = [float(tensor_match.group(pos)) for pos in (2, 3, 4)]
                for col, value in zip(("x", "y", "z"), values):
                    block["tensor_au"][f"{row}{col}"] = value
                idx += 1
                continue

            vec_match = _BETA_VEC_RE.match(line)
            if vec_match:
                axis = vec_match.group(1).lower()
                value_au = float(vec_match.group(2))
                block[f"beta_vec_{axis}_au"] = value_au
                block[f"beta_vec_{axis}_esu"] = value_au * _AU_TO_ESU
                block[f"beta_vec_{axis}_esu_30"] = value_au * _AU_TO_ESU_30
                idx += 1
                continue

            sq_match = _BETA_SQ_RE.match(line)
            if sq_match:
                block["beta_sq_zzz_au2"] = float(sq_match.group(1))
                block["beta_sq_xzz_au2"] = float(sq_match.group(2))
                idx += 1
                continue

            hrs_match = _BETA_HRS_RE.match(line)
            if hrs_match:
                block["beta_hrs_au"] = float(hrs_match.group(1))
                block["beta_hrs_esu"] = block["beta_hrs_au"] * _AU_TO_ESU
                block["beta_hrs_esu_30"] = block["beta_hrs_au"] * _AU_TO_ESU_30
                block["dr"] = float(hrs_match.group(2))
                idx += 1
                break

            idx += 1

        if len(block["tensor_au"]) == 27 or "beta_hrs_au" in block:
            blocks.append(block)

    return blocks


def _merge_response_points(
    beta_hrs_points: list[dict[str, Any]],
    blocks: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    merged: list[dict[str, Any]] = []
    count = max(len(beta_hrs_points), len(blocks))
    for idx in range(count):
        point: dict[str, Any] = {}
        if idx < len(beta_hrs_points):
            point.update(beta_hrs_points[idx])
        if idx < len(blocks):
            point.update(blocks[idx])
        if "wavelength_nm" not in point and "response_wavelength_nm" in point:
            point["wavelength_nm"] = point["response_wavelength_nm"]
        wavelength_nm = point.get("wavelength_nm")
        if wavelength_nm == math.inf:
            point["is_static"] = True
        merged.append(point)
    deduped: list[dict[str, Any]] = []
    for point in merged:
        if deduped:
            previous = deduped[-1]
            same_wavelength = previous.get("wavelength_nm") == point.get("wavelength_nm")
            same_beta = previous.get("beta_hrs_au") == point.get("beta_hrs_au")
            same_tensor = previous.get("tensor_au") == point.get("tensor_au")
            if same_wavelength and same_beta and same_tensor:
                continue
        deduped.append(point)
    return deduped


def _annotate_tensor_metrics(
    response_points: list[dict[str, Any]],
    *,
    dipole_x_au: Optional[float],
    dipole_y_au: Optional[float],
    dipole_z_au: Optional[float],
) -> list[dict[str, Any]]:
    if dipole_x_au is None or dipole_y_au is None or dipole_z_au is None:
        return response_points

    for point in response_points:
        tensor = point.get("tensor_au")
        if not isinstance(tensor, dict) or len(tensor) < 27:
            continue
        beta_zzz_au = float(tensor.get("zzz", 0.0))
        beta_zzz_aligned_au, kleinman_score, kleinman_applied = calculate_beta_zzz_aligned(
            tensor,
            dipole_x_au,
            dipole_y_au,
            dipole_z_au,
            kleinman_mode="off",
        )
        point["beta_zzz_au"] = beta_zzz_au
        point["beta_zzz_esu"] = beta_zzz_au * _AU_TO_ESU
        point["beta_zzz_esu_30"] = beta_zzz_au * _AU_TO_ESU_30
        point["beta_zzz_aligned_au"] = beta_zzz_aligned_au
        point["beta_zzz_aligned_esu"] = beta_zzz_aligned_au * _AU_TO_ESU
        point["beta_zzz_aligned_esu_30"] = beta_zzz_aligned_au * _AU_TO_ESU_30
        point["dipole_x_au"] = dipole_x_au
        point["dipole_y_au"] = dipole_y_au
        point["dipole_z_au"] = dipole_z_au
        point["kleinman_score_raw_rotation"] = kleinman_score
        point["kleinman_applied_raw_rotation"] = kleinman_applied
    return response_points


def _normalize_response_points(
    response_points: list[dict[str, Any]],
    *,
    wavelengths_nm: list[float],
) -> list[dict[str, Any]]:
    static_only = len(wavelengths_nm) == 1 and math.isinf(wavelengths_nm[0])
    if not static_only:
        return response_points

    static_points = [point for point in response_points if point.get("is_static")]
    if len(static_points) <= 1:
        return response_points

    first_static = static_points[0]
    if any(point != first_static for point in static_points[1:]):
        logger.warning(
            "response backend produced %d static points for static-only mode; keeping the first one",
            len(static_points),
        )

    normalized = [first_static]
    normalized.extend(point for point in response_points if not point.get("is_static"))
    return normalized


def _extract_response_points(
    output_path: Path,
    beta_hrs_path: Path,
) -> list[dict[str, Any]]:
    return _merge_response_points(
        _parse_beta_hrs_file(beta_hrs_path),
        _parse_response_blocks(output_path),
    )


def run_single_hyperpol_workflow(
    entry: WorkflowEntry,
    *,
    charge: int,
    multiplicity: int,
    engine: str,
    preopt: str,
    wavelengths_nm: list[float],
    energy_window_ev: float,
    cores: int,
    workdir: Path,
) -> HyperpolResult:
    workdir.mkdir(parents=True, exist_ok=True)
    initial_xyz = _write_inputs(entry, workdir)
    current_xyz = initial_xyz
    preopt_output: Optional[Path] = None
    dipole_output: Optional[Path] = None
    dipole_x_au: Optional[float] = None
    dipole_y_au: Optional[float] = None
    dipole_z_au: Optional[float] = None
    dipole_total_debye: Optional[float] = None

    if preopt == "xtb":
        preopt_output, current_xyz, _xtb_energy = _run_xtb_opt(
            current_xyz,
            workdir,
            charge=charge,
            multiplicity=multiplicity,
            cores=cores,
            stage_name="xtb_opt",
            output_name="xtbopt.xyz",
        )
    elif preopt == "crest":
        preopt_output, current_xyz = _run_crest(
            current_xyz,
            workdir,
            charge=charge,
            multiplicity=multiplicity,
            cores=cores,
        )

    dipole_output, _xtb_sp_energy = _run_xtb_singlepoint(
        current_xyz,
        workdir,
        charge=charge,
        multiplicity=multiplicity,
        cores=cores,
        stage_name="dipole_xtb_sp",
    )
    dipole_x_au, dipole_y_au, dipole_z_au, dipole_total_debye = _parse_xtb_dipole(dipole_output)

    wavelength_file = _write_wavelength_file(workdir, wavelengths_nm)
    xtb4stda_output = _run_xtb4stda(
        current_xyz,
        workdir,
        charge=charge,
        multiplicity=multiplicity,
        cores=cores,
    )
    response_output, beta_hrs_file, beta_tensor_file = _run_response(
        workdir,
        engine=engine,
        energy_window_ev=energy_window_ev,
        wavelength_count=len(wavelengths_nm),
        cores=cores,
    )
    response_points = _annotate_tensor_metrics(
        _extract_response_points(response_output, beta_hrs_file),
        dipole_x_au=dipole_x_au,
        dipole_y_au=dipole_y_au,
        dipole_z_au=dipole_z_au,
    )
    response_points = _normalize_response_points(
        response_points,
        wavelengths_nm=wavelengths_nm,
    )

    return HyperpolResult(
        label=entry.label,
        smiles=entry.smiles,
        charge=charge,
        multiplicity=multiplicity,
        engine=engine,
        preopt=preopt,
        energy_window_ev=energy_window_ev,
        requested_wavelengths_nm=list(wavelengths_nm),
        workdir=str(workdir),
        initial_xyz=str(initial_xyz),
        selected_xyz=str(current_xyz),
        preopt_output=str(preopt_output) if preopt_output else None,
        dipole_output=str(dipole_output),
        dipole_x_au=dipole_x_au,
        dipole_y_au=dipole_y_au,
        dipole_z_au=dipole_z_au,
        dipole_total_debye=dipole_total_debye,
        xtb4stda_output=str(xtb4stda_output),
        response_output=str(response_output),
        wavelength_file=str(wavelength_file),
        beta_hrs_file=str(beta_hrs_file),
        beta_tensor_file=str(beta_tensor_file),
        response_points=response_points,
    )


def _load_wavelengths(
    wavelength_values: list[float],
    wavelength_file: Optional[str],
    *,
    static_only: bool,
) -> list[float]:
    if static_only:
        if wavelength_values or wavelength_file:
            raise ValueError("--static-only cannot be combined with --wavelength or --wavelength-file")
        return [math.inf]
    values = [float(value) for value in wavelength_values]
    if wavelength_file:
        path = resolve_path(wavelength_file)
        for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            values.append(float(line))
    cleaned: list[float] = []
    for value in values:
        if value <= 0.0:
            raise ValueError("response wavelengths must be positive in nm")
        cleaned.append(value)
    return cleaned or [1064.0]


def _print_result(result: HyperpolResult) -> None:
    def _fmt_beta(label: str, point: dict[str, Any], au_key: str) -> str:
        value_au = point.get(au_key)
        if value_au is None:
            return f"{label}=n/a"
        esu_key = au_key.replace("_au", "_esu")
        esu_30_key = au_key.replace("_au", "_esu_30")
        value_esu = point.get(esu_key)
        value_esu_30 = point.get(esu_30_key)
        return (
            f"{label}={float(value_au):.3f} au "
            f"({float(value_esu):.3e} esu; {float(value_esu_30):.3f} x10^-30 esu)"
        )

    parts = [f"engine={result.engine}", f"preopt={result.preopt}"]
    static_point = next(
        (item for item in result.response_points if item.get("is_static")),
        None,
    )
    if static_point is not None:
        parts.append(_fmt_beta("beta_HRS(static)", static_point, "beta_hrs_au"))
        if static_point.get("beta_zzz_au") is not None:
            parts.append(_fmt_beta("beta_zzz(static)", static_point, "beta_zzz_au"))
        if static_point.get("beta_zzz_aligned_au") is not None:
            parts.append(_fmt_beta("beta_zzz_aligned(static)", static_point, "beta_zzz_aligned_au"))
    for point in result.response_points:
        if point.get("is_static"):
            continue
        wavelength_nm = point.get("wavelength_nm")
        if wavelength_nm in (None, math.inf):
            continue
        dr = point.get("dr")
        segment = _fmt_beta(f"beta_HRS({float(wavelength_nm):.1f} nm)", point, "beta_hrs_au")
        if dr is not None:
            segment += f" (DR={float(dr):.3f})"
        parts.append(segment)
        if point.get("beta_zzz_au") is not None:
            parts.append(_fmt_beta(f"beta_zzz({float(wavelength_nm):.1f} nm)", point, "beta_zzz_au"))
        if point.get("beta_zzz_aligned_au") is not None:
            parts.append(_fmt_beta(f"beta_zzz_aligned({float(wavelength_nm):.1f} nm)", point, "beta_zzz_aligned_au"))
    parts.append(f"workdir={result.workdir}")
    print(f"{result.label}: " + ", ".join(parts))


def run_cli(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin hyperpol_xtb",
        description=(
            "SMILES/XYZ -> optional xTB/CREST preoptimization -> xtb4stda -> "
            "std2/stda nonlinear-response SHG first hyperpolarizability workflow."
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
    parser.add_argument("--workdir", default="hyperpol_xtb_runs", help="Base output directory.")
    parser.add_argument("--charge", type=int, default=0, help="Molecular charge.")
    parser.add_argument("--multiplicity", type=int, default=1, help="Ground-state multiplicity.")
    parser.add_argument(
        "--engine",
        choices=("std2", "stda"),
        default="std2",
        help="Response engine: std2 = sTD-DFT-xTB (paper default), stda = sTDA-xTB.",
    )
    parser.add_argument(
        "--preopt",
        choices=("xtb", "crest", "none"),
        default="xtb",
        help="Geometry preparation before xtb4stda.",
    )
    parser.add_argument(
        "--wavelength",
        action="append",
        type=float,
        default=[],
        help="Fundamental SHG wavelength in nm. Repeatable. Static (infinity) is emitted automatically by std2/stda -resp.",
    )
    parser.add_argument(
        "--wavelength-file",
        help="Optional file with one positive wavelength in nm per line.",
    )
    parser.add_argument(
        "--static-only",
        action="store_true",
        help="Run only the static limit (infinite wavelength, 0 eV).",
    )
    parser.add_argument(
        "--energy-window",
        type=float,
        default=15.0,
        help="Excitation window in eV for the response calculation. Paper-like default: 15 eV.",
    )
    parser.add_argument("--pal", type=int, default=4, help="Cores per molecule workflow.")
    parser.add_argument("--parallel-jobs", type=int, default=1, help="Number of molecule workflows to run in parallel.")
    parser.add_argument("--maxcore", type=int, default=1000, help="Memory per core in MB for pool scheduling.")
    parser.add_argument("--json-out", help="Optional JSON summary path.")
    args = parser.parse_args(argv)

    try:
        wavelengths_nm = _load_wavelengths(
            args.wavelength,
            args.wavelength_file,
            static_only=args.static_only,
        )
    except ValueError as exc:
        parser.error(str(exc))

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
        preopt=args.preopt,
        engine=args.engine,
    )
    if args.check:
        _print_preflight_checks(preflight_checks)
        return 0 if all(check.ok for check in preflight_checks) else 1

    failed_checks = [check for check in preflight_checks if not check.ok]
    if failed_checks:
        _print_preflight_checks(failed_checks)
        print("hyperpol_xtb preflight failed. Run `delfin hyperpol_xtb --check` after fixing the missing items.")
        return 1

    entries = _load_entries(args.smiles, args.smiles_file, args.xyz, args.xyz_file, args.label)
    if not entries:
        parser.error("provide --smiles, --smiles-file, --xyz, and/or --xyz-file")

    base_workdir = resolve_path(args.workdir)
    base_workdir.mkdir(parents=True, exist_ok=True)

    results: list[HyperpolResult] = []
    errors: dict[str, str] = {}

    if len(entries) == 1 or args.parallel_jobs <= 1:
        for entry in entries:
            try:
                result = run_single_hyperpol_workflow(
                    entry,
                    charge=args.charge,
                    multiplicity=args.multiplicity,
                    engine=args.engine,
                    preopt=args.preopt,
                    wavelengths_nm=wavelengths_nm,
                    energy_window_ev=args.energy_window,
                    cores=args.pal,
                    workdir=base_workdir / _safe_label(entry.label, entry.label),
                )
                results.append(result)
                _print_result(result)
            except Exception as exc:  # noqa: BLE001
                errors[entry.label] = str(exc)
                logger.error("[%s] hyperpol_xtb failed: %s", entry.label, exc, exc_info=True)
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
            workdir = base_workdir / _safe_label(entry.label, entry.label)

            def runner(*_args, _entry=entry, _workdir=workdir, **kwargs) -> None:
                allocated = kwargs.get("cores", args.pal)
                try:
                    result = run_single_hyperpol_workflow(
                        _entry,
                        charge=args.charge,
                        multiplicity=args.multiplicity,
                        engine=args.engine,
                        preopt=args.preopt,
                        wavelengths_nm=wavelengths_nm,
                        energy_window_ev=args.energy_window,
                        cores=max(1, int(allocated)),
                        workdir=_workdir,
                    )
                    with lock:
                        results.append(result)
                except Exception as exc:  # noqa: BLE001
                    with lock:
                        errors[_entry.label] = str(exc)
                    raise

            pool_job = PoolJob(
                job_id=f"HYPERPOL_XTB_{entry.label}",
                cores_min=1,
                cores_optimal=max(1, int(args.pal)),
                cores_max=max(1, int(args.pal)),
                memory_mb=max(256, int(args.maxcore) * max(1, int(args.pal))),
                priority=JobPriority.NORMAL,
                execute_func=runner,
                args=(),
                kwargs={},
                estimated_duration=3600.0,
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
