"""Helpers for running single-file XYZ workflows launched from the browser UI."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

from delfin.common.paths import resolve_path
from delfin.dashboard.input_processing import smiles_to_xyz_quick
from delfin.ensemble_nmr import (
    build_anmrrc_text,
    build_censo_anmr_rc,
    build_orca_reference_input,
    normalize_cpcm_solvent_name,
    write_text_file,
    xyz_body_to_coord_text,
)
from delfin.global_manager import get_global_manager
from delfin.hyperpol import (
    WorkflowEntry as HyperpolWorkflowEntry,
    _load_wavelengths as _hyperpol_load_wavelengths,
    _print_result as _print_hyperpol_result,
    run_single_hyperpol_workflow,
)
from delfin.nmr_spectrum import parse_nmr_orca
from delfin.qm_runtime import resolve_tool
from delfin.runtime_setup import run_analysis_tools_installer
from delfin.tadf_xtb import (
    WorkflowEntry as TadfWorkflowEntry,
    _print_result as _print_tadf_result,
    run_single_tadf_xtb,
)


def _resolve_active_job_limits(*, requested_cores: int, requested_maxcore: int) -> tuple[int, int]:
    manager = get_global_manager()
    try:
        if manager.is_initialized():
            return manager.resolve_job_resources(
                requested_cores=requested_cores,
                requested_maxcore=requested_maxcore,
            )
    except Exception:
        pass
    return max(1, int(requested_cores)), max(256, int(requested_maxcore))


def _write_json(path: str | None, payload: dict) -> None:
    if not path:
        return
    target = resolve_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def _default_summary_paths(workdir: Path) -> tuple[Path, Path]:
    stem = workdir.name or "workflow"
    return workdir / f"{stem}.txt", workdir / f"{stem}.json"


def _fmt(value: object, digits: int = 3) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        if math.isinf(value):
            return "inf"
        return f"{value:.{digits}f}"
    return str(value)


def _hyperpol_summary_payload(result) -> dict:
    static_point = next(
        (item for item in result.response_points if item.get("is_static")),
        None,
    )
    return {
        "label": result.label,
        "engine": result.engine,
        "preopt": result.preopt,
        "charge": result.charge,
        "multiplicity": result.multiplicity,
        "dipole_total_debye": result.dipole_total_debye,
        "beta_hrs_static_au": None if static_point is None else static_point.get("beta_hrs_au"),
        "beta_hrs_static_esu": None if static_point is None else static_point.get("beta_hrs_esu"),
        "beta_hrs_static_esu_30": None if static_point is None else static_point.get("beta_hrs_esu_30"),
        "beta_zzz_static_au": None if static_point is None else static_point.get("beta_zzz_au"),
        "beta_zzz_aligned_static_au": None if static_point is None else static_point.get("beta_zzz_aligned_au"),
        "selected_xyz": result.selected_xyz,
        "beta_hrs_file": result.beta_hrs_file,
        "beta_tensor_file": result.beta_tensor_file,
        "response_output": result.response_output,
        "response_points": result.response_points,
    }


def _hyperpol_summary_text(result) -> str:
    summary = _hyperpol_summary_payload(result)
    lines = [
        f"label: {summary['label']}",
        f"engine: {summary['engine']}",
        f"preopt: {summary['preopt']}",
        f"charge: {summary['charge']}",
        f"multiplicity: {summary['multiplicity']}",
        f"dipole_total_debye: {_fmt(summary['dipole_total_debye'], 4)}",
        f"beta_HRS_static_au: {_fmt(summary['beta_hrs_static_au'], 6)}",
        f"beta_HRS_static_esu: {_fmt(summary['beta_hrs_static_esu'], 6)}",
        f"beta_HRS_static_x10^-30_esu: {_fmt(summary['beta_hrs_static_esu_30'], 6)}",
        f"beta_zzz_static_au: {_fmt(summary['beta_zzz_static_au'], 6)}",
        f"beta_zzz_aligned_static_au: {_fmt(summary['beta_zzz_aligned_static_au'], 6)}",
        f"selected_xyz: {summary['selected_xyz']}",
        f"beta_hrs_file: {summary['beta_hrs_file']}",
        f"beta_tensor_file: {summary['beta_tensor_file']}",
        f"response_output: {summary['response_output']}",
    ]
    return "\n".join(lines)


def _tadf_summary_payload(result) -> dict:
    return {
        "label": result.label,
        "excited_method": result.excited_method,
        "charge": result.charge,
        "multiplicity": result.multiplicity,
        "s1_state": result.s1_state,
        "s1_ev": result.s1_ev,
        "s1_nm": result.s1_nm,
        "s1_f_osc": result.s1_f_osc,
        "first_allowed_singlet_state": result.first_allowed_singlet_state,
        "first_allowed_singlet_ev": result.first_allowed_singlet_ev,
        "first_allowed_singlet_nm": result.first_allowed_singlet_nm,
        "first_allowed_singlet_f_osc": result.first_allowed_singlet_f_osc,
        "brightest_singlet_state": result.brightest_singlet_state,
        "brightest_singlet_ev": result.brightest_singlet_ev,
        "brightest_singlet_nm": result.brightest_singlet_nm,
        "brightest_singlet_f_osc": result.brightest_singlet_f_osc,
        "t1_ev": result.t1_ev,
        "t1_nm": result.t1_nm,
        "delta_est_vert_ev": result.delta_est_ev,
        "lambda_abs_nm": result.lambda_abs_nm,
        "selected_xyz": result.s0_xyz,
        "singlet_output": result.singlet_output,
        "triplet_output": result.triplet_output,
    }


def _tadf_summary_text(result) -> str:
    summary = _tadf_summary_payload(result)
    lines = [
        f"label: {summary['label']}",
        f"excited_method: {summary['excited_method']}",
        f"charge: {summary['charge']}",
        f"multiplicity: {summary['multiplicity']}",
        f"S1_state: {_fmt(summary['s1_state'])}",
        f"S1_eV: {_fmt(summary['s1_ev'], 6)}",
        f"S1_nm: {_fmt(summary['s1_nm'], 3)}",
        f"S1_f_osc: {_fmt(summary['s1_f_osc'], 6)}",
        f"first_allowed_singlet_state: {_fmt(summary['first_allowed_singlet_state'])}",
        f"first_allowed_singlet_eV: {_fmt(summary['first_allowed_singlet_ev'], 6)}",
        f"first_allowed_singlet_nm: {_fmt(summary['first_allowed_singlet_nm'], 3)}",
        f"first_allowed_singlet_f_osc: {_fmt(summary['first_allowed_singlet_f_osc'], 6)}",
        f"brightest_singlet_state: {_fmt(summary['brightest_singlet_state'])}",
        f"brightest_singlet_eV: {_fmt(summary['brightest_singlet_ev'], 6)}",
        f"brightest_singlet_nm: {_fmt(summary['brightest_singlet_nm'], 3)}",
        f"brightest_singlet_f_osc: {_fmt(summary['brightest_singlet_f_osc'], 6)}",
        f"T1_eV: {_fmt(summary['t1_ev'], 6)}",
        f"T1_nm: {_fmt(summary['t1_nm'], 3)}",
        f"DeltaEST_vert_eV: {_fmt(summary['delta_est_vert_ev'], 6)}",
        f"lambda_abs_nm: {_fmt(summary['lambda_abs_nm'], 3)}",
        f"selected_xyz: {summary['selected_xyz']}",
        f"singlet_output: {summary['singlet_output']}",
        f"triplet_output: {summary['triplet_output']}",
    ]
    return "\n".join(lines)


def _hyperpol_formatted_report(result) -> str:
    """Generate DELFIN-style formatted hyperpol_xtb report (same as pipeline)."""
    def _fmt_beta(label, point, au_key):
        v_au = point.get(au_key)
        if v_au is None:
            return f"  {label:<36s} = n/a"
        v_esu = point.get(au_key.replace("_au", "_esu"))
        v_esu_30 = point.get(au_key.replace("_au", "_esu_30"))
        return (
            f"  {label:<36s} = {float(v_au):>12.3f} au"
            f"    ({float(v_esu):.3e} esu; {float(v_esu_30):.3f} x10^-30 esu)"
        )

    sep = "-" * 72
    lines = [
        sep,
        "DELFIN \u2014 xTB Hyperpolarizability (sTD-DFT-xTB)",
        sep,
        "",
        "Calculation settings:",
        f"  Label          = {result.label}",
        f"  Engine         = {result.engine}",
        f"  Preopt         = {result.preopt}",
        f"  Charge         = {result.charge}",
        f"  Multiplicity   = {result.multiplicity}",
        f"  Energy window  = {result.energy_window_ev} eV",
        "",
        sep,
        "Dipole moment:",
        f"  mu_x           = {result.dipole_x_au:.6f} au" if result.dipole_x_au is not None else "  mu_x           = n/a",
        f"  mu_y           = {result.dipole_y_au:.6f} au" if result.dipole_y_au is not None else "  mu_y           = n/a",
        f"  mu_z           = {result.dipole_z_au:.6f} au" if result.dipole_z_au is not None else "  mu_z           = n/a",
        f"  |mu|           = {result.dipole_total_debye:.4f} Debye" if result.dipole_total_debye is not None else "  |mu|           = n/a",
        "",
    ]
    for point in result.response_points:
        wl = point.get("wavelength_nm")
        is_static = point.get("is_static", False)
        tag = "static" if is_static else f"{float(wl):.1f} nm"
        lines.append(sep)
        lines.append(f"First hyperpolarizability ({tag}):")
        lines.append("")
        lines.append(_fmt_beta(f"beta_HRS({tag})", point, "beta_hrs_au"))
        if point.get("beta_zzz_au") is not None:
            lines.append(_fmt_beta(f"beta_zzz({tag})", point, "beta_zzz_au"))
        if point.get("beta_zzz_aligned_au") is not None:
            lines.append(_fmt_beta(f"|beta_zzz| aligned to dipole({tag})", point, "beta_zzz_aligned_au"))
        for axis in ("x", "y", "z"):
            key = f"beta_vec_{axis}_au"
            if point.get(key) is not None:
                lines.append(_fmt_beta(f"beta_vec_{axis}({tag})", point, key))
        if point.get("dr") is not None:
            lines.append(f"  {'Depolarization ratio (DR)':<36s} = {float(point['dr']):>12.3f}")
        if point.get("beta_sq_zzz_au2") is not None:
            lines.append(f"  {'<beta_zzz^2>':<36s} = {float(point['beta_sq_zzz_au2']):>12.3f} au^2")
        if point.get("beta_sq_xzz_au2") is not None:
            lines.append(f"  {'<beta_xzz^2>':<36s} = {float(point['beta_sq_xzz_au2']):>12.3f} au^2")
        if point.get("kleinman_score_raw_rotation") is not None:
            lines.append(f"  {'Kleinman symmetry score':<36s} = {float(point['kleinman_score_raw_rotation']):>12.6f}")
        lines.append("")
    lines.append(sep)
    lines.append(f"Working directory: {result.workdir}")
    lines.append(sep)
    return "\n".join(lines)


def _tadf_formatted_report(result, preopt: str = "none", energy_window: float = 10.0, run_t1_opt: bool = True) -> str:
    """Generate DELFIN-style formatted tadf_xtb report (same as pipeline)."""
    sep = "-" * 72
    lines = [
        sep,
        "DELFIN \u2014 xTB TADF Screening (sTD-DFT-xTB)",
        sep,
        "",
        "Calculation settings:",
        f"  Label             = {result.label}",
        f"  Excited method    = {result.excited_method}",
        f"  xTB method        = {result.xtb_method}",
        f"  Charge            = {result.charge}",
        f"  Multiplicity      = {result.multiplicity}",
        f"  Preopt            = {preopt}",
        f"  Energy window     = {energy_window} eV",
        f"  T1 optimization   = {'yes' if run_t1_opt else 'no'}",
        "",
        sep,
        "Vertical excited states (at S0 geometry):",
        f"  S1                = {_fmt(result.s1_ev, 4):>10s} eV    ({_fmt(result.s1_nm, 1):>8s} nm,  f = {_fmt(result.s1_f_osc, 5)})",
        f"  T1                = {_fmt(result.t1_ev, 4):>10s} eV    ({_fmt(result.t1_nm, 1):>8s} nm)",
        f"  Delta_E(S-T) vert = {_fmt(result.delta_est_ev, 4):>10s} eV",
        "",
    ]

    if result.first_allowed_singlet_ev is not None:
        lines.append(
            f"  First allowed     = S{result.first_allowed_singlet_state:<3d}"
            f" {_fmt(result.first_allowed_singlet_ev, 4):>10s} eV"
            f"    ({_fmt(result.first_allowed_singlet_nm, 1):>8s} nm,"
            f"  f = {_fmt(result.first_allowed_singlet_f_osc, 5)})"
        )
    if result.brightest_singlet_ev is not None:
        lines.append(
            f"  Brightest singlet = S{result.brightest_singlet_state:<3d}"
            f" {_fmt(result.brightest_singlet_ev, 4):>10s} eV"
            f"    ({_fmt(result.brightest_singlet_nm, 1):>8s} nm,"
            f"  f = {_fmt(result.brightest_singlet_f_osc, 5)})"
        )
    if result.first_allowed_singlet_ev is not None or result.brightest_singlet_ev is not None:
        lines.append("")

    lines.append(sep)
    lines.append("Total energies:")
    if result.s0_xtb_energy_eh is not None:
        lines.append(f"  E(S0, xTB)        = {_fmt(result.s0_xtb_energy_eh, 8):>16s} Eh")
    if result.t1_xtb_energy_eh is not None:
        lines.append(f"  E(T1, xTB)        = {_fmt(result.t1_xtb_energy_eh, 8):>16s} Eh")
    lines.append("")

    lines.append(sep)
    lines.append("Adiabatic & relaxed energetics:")
    if result.t1_vertical_xtb_ev is not None:
        lines.append(f"  E(T1, vertical)   = {_fmt(result.t1_vertical_xtb_ev, 4):>10s} eV")
    if result.t1_adiabatic_ev is not None:
        lines.append(f"  E(T1, adiabatic)  = {_fmt(result.t1_adiabatic_ev, 4):>10s} eV")
    if result.t1_relaxed_ev is not None:
        lines.append(f"  E(T1, relaxed)    = {_fmt(result.t1_relaxed_ev, 4):>10s} eV")
    if result.t1_relaxation_ev is not None:
        lines.append(f"  T1 relaxation     = {_fmt(result.t1_relaxation_ev, 4):>10s} eV")
    if result.s1_relaxed_est_ev is not None:
        lines.append(f"  E(S1, rel. est.)  = {_fmt(result.s1_relaxed_est_ev, 4):>10s} eV")
    if result.delta_est_relaxed_est_ev is not None:
        lines.append(f"  Delta_E(S-T) rel. = {_fmt(result.delta_est_relaxed_est_ev, 4):>10s} eV")
    lines.append("")

    lines.append(sep)
    lines.append("Photophysical properties:")
    lines.append(f"  lambda_abs(S0->S1)  = {_fmt(result.lambda_abs_nm, 1):>8s} nm")
    if result.lambda_em_s1_nm is not None:
        lines.append(f"  lambda_em(S1)       = {_fmt(result.lambda_em_s1_nm, 1):>8s} nm")
    lines.append(f"  lambda_PL(est.)     = {_fmt(result.lambda_pl_est_nm, 1):>8s} nm")
    if result.stokes_shift_est_ev is not None:
        lines.append(f"  Stokes shift (est.) = {_fmt(result.stokes_shift_est_ev, 4):>10s} eV")
    lines.append("")

    lines.append(sep)
    lines.append(f"Working directory: {result.workdir}")
    lines.append(sep)
    return "\n".join(lines)


def _write_success_summaries(
    workdir: Path,
    *,
    payload: dict,
    summary_text: str,
    json_out: str | None,
) -> None:
    txt_path, json_path = _default_summary_paths(workdir)
    _write_text(txt_path, summary_text)
    _write_json(str(json_path), payload)
    if json_out:
        resolved_json = resolve_path(json_out)
        if resolved_json != json_path:
            _write_json(str(resolved_json), payload)


def _write_error_summaries(
    workdir: Path,
    *,
    label: str,
    error_message: str,
    json_out: str | None,
) -> None:
    txt_path, json_path = _default_summary_paths(workdir)
    payload = {"result": None, "errors": {label: error_message}}
    _write_text(txt_path, f"label: {label}\nstatus: failed\nerror: {error_message}")
    _write_json(str(json_path), payload)
    if json_out:
        resolved_json = resolve_path(json_out)
        if resolved_json != json_path:
            _write_json(str(resolved_json), payload)


def _resolved_path_or_empty(tool_name: str) -> str:
    resolved = resolve_tool(tool_name)
    return resolved.path if resolved is not None else ""


_AUTO_INSTALLABLE_ANALYSIS_TOOLS = frozenset({"censo", "anmr", "c2anmr", "nmrplot"})


def _minimum_supported_censo_version() -> tuple[int, int, int]:
    return (3, 0, 0) if sys.version_info >= (3, 12) else (2, 1, 4)


def _censo_requires_upgrade(censo_path: str) -> bool:
    version = _detect_censo_version(censo_path)
    if version is None:
        return False
    return version < _minimum_supported_censo_version()


def _censo_refinement_threshold(censo_path: str) -> float:
    version = _detect_censo_version(censo_path)
    if version is not None and version >= (3, 0, 0):
        return 0.001
    return 0.0


def _auto_install_analysis_tools_enabled() -> bool:
    value = str(os.environ.get("DELFIN_AUTO_INSTALL_ANALYSIS_TOOLS", "1")).strip().lower()
    return value not in {"0", "false", "no", "off"}


def _maybe_auto_install_analysis_tools(missing_tools: list[str], *, workdir: Path) -> list[str]:
    supported_missing = [
        tool_name for tool_name in missing_tools if tool_name in _AUTO_INSTALLABLE_ANALYSIS_TOOLS
    ]
    if not supported_missing or not _auto_install_analysis_tools_enabled():
        return []

    installer_env = {
        "INSTALL_MORFEUS": "0",
        "INSTALL_CCLIB": "0",
        "INSTALL_NGLVIEW": "0",
        "INSTALL_PACKMOL": "0",
        "INSTALL_MULTIWFN": "0",
        "CENSO_PREFER_LATEST": "1",
        "INSTALL_CENSO": "1" if any(
            tool in supported_missing for tool in ("censo", "c2anmr", "nmrplot")
        ) else "0",
        "INSTALL_ANMR": "1" if "anmr" in supported_missing else "0",
        "FORCE_REINSTALL": "1" if "censo" in supported_missing else "0",
    }
    target, result = run_analysis_tools_installer(extra_env=installer_env)
    log_lines = [
        "Auto-install missing analysis tools",
        f"Requested: {', '.join(supported_missing)}",
        f"Target: {target}",
        "",
        result.stdout or "(no installer output)",
    ]
    write_text_file(workdir / "analysis_tools_auto_install.log", "\n".join(log_lines))
    if result.returncode != 0:
        raise RuntimeError(
            "Automatic installation of missing analysis tools failed. "
            f"See {workdir / 'analysis_tools_auto_install.log'} for details."
        )

    active_bin = str(Path(sys.executable).resolve().parent)
    current_path = os.environ.get("PATH", "")
    path_parts = current_path.split(os.pathsep) if current_path else []
    if active_bin not in path_parts:
        os.environ["PATH"] = active_bin if not current_path else f"{active_bin}{os.pathsep}{current_path}"

    installed_now: list[str] = []
    for tool_name in supported_missing:
        if _resolved_path_or_empty(tool_name):
            installed_now.append(tool_name)
    return installed_now


def _require_tools(tool_names: list[str], *, workdir: Path | None = None) -> dict[str, str]:
    resolved: dict[str, str] = {}
    missing: list[str] = []
    for tool_name in tool_names:
        path = _resolved_path_or_empty(tool_name)
        if tool_name == "censo" and path and _censo_requires_upgrade(path):
            missing.append(tool_name)
        elif path:
            resolved[tool_name] = path
        else:
            missing.append(tool_name)
    if missing and workdir is not None:
        _maybe_auto_install_analysis_tools(missing, workdir=workdir)
        resolved.clear()
        missing = []
        for tool_name in tool_names:
            path = _resolved_path_or_empty(tool_name)
            if tool_name == "censo" and path and _censo_requires_upgrade(path):
                missing.append(tool_name)
            elif path:
                resolved[tool_name] = path
            else:
                missing.append(tool_name)
    if missing:
        raise RuntimeError(
            "Missing required tools for CENSO/ANMR workflow: "
            + ", ".join(missing)
            + ". Configure them in Settings or ensure they are on PATH. "
            + "DELFIN auto-installs supported analysis tools by default; unsupported tools such as crest, xtb, and orca must already exist."
        )
    return resolved


def _run_logged(
    cmd: list[str],
    *,
    cwd: Path,
    log_path: Path,
    env: dict[str, str] | None = None,
    input_text: str | None = None,
) -> None:
    env_map = os.environ.copy()
    if env:
        env_map.update({str(k): str(v) for k, v in env.items() if str(v).strip()})
    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env_map,
        capture_output=True,
        text=True,
        input=input_text,
        check=False,
    )
    log_text = (
        "$ " + " ".join(shlex.quote(part) for part in cmd) + "\n\n"
        + result.stdout
        + ("\n" if result.stdout and not result.stdout.endswith("\n") else "")
        + result.stderr
    )
    write_text_file(log_path, log_text)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(shlex.quote(part) for part in cmd)}"
        )


def _resolve_orca_reference_result_path(*, workdir: Path, stem: str) -> Path:
    """Return the best ORCA result file for downstream parsing."""
    out_path = workdir / f"{stem}.out"
    if out_path.is_file():
        return out_path

    log_path = workdir / f"{stem}.log"
    if log_path.is_file():
        log_text = log_path.read_text(encoding="utf-8", errors="replace")
        if "ORCA TERMINATED NORMALLY" in log_text:
            return log_path

    raise RuntimeError(f"ORCA reference run finished without {stem}.out")


def _extract_xyz_body_from_orca_input(inp_path: Path) -> str:
    """Extract the `* xyz` geometry block from an ORCA input file."""
    lines = inp_path.read_text(encoding="utf-8", errors="replace").splitlines()
    in_xyz = False
    xyz_lines: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not in_xyz:
            if stripped.lower().startswith("* xyz"):
                in_xyz = True
            continue
        if stripped == "*":
            break
        if stripped:
            xyz_lines.append(stripped)
    return "\n".join(xyz_lines)


def _ensure_censo_nmr_coords(workdir: Path) -> None:
    """Backfill missing `coord` files for CENSO NMR conformers from ORCA inputs."""
    nmr_root = workdir / "4_NMR"
    if not nmr_root.is_dir():
        return

    for conf_dir in sorted(path for path in nmr_root.iterdir() if path.is_dir()):
        coord_path = conf_dir / "coord"
        if coord_path.is_file():
            continue
        inp_path = conf_dir / "nmr" / "nmr.inp"
        if not inp_path.is_file():
            continue
        xyz_body = _extract_xyz_body_from_orca_input(inp_path)
        if not xyz_body.strip():
            continue
        write_text_file(coord_path, xyz_body_to_coord_text(xyz_body))


def _ensure_anmr_coord_inputs(*, workdir: Path, anmr_dir: Path) -> None:
    """Copy CENSO NMR `coord` files into the ANMR workspace layout."""
    nmr_root = workdir / "4_NMR"
    if not nmr_root.is_dir() or not anmr_dir.is_dir():
        return

    first_coord_text = ""
    for conf_dir in sorted(path for path in anmr_dir.iterdir() if path.is_dir() and path.name.startswith("CONF")):
        source_coord = nmr_root / conf_dir.name / "coord"
        if not source_coord.is_file():
            continue
        coord_text = source_coord.read_text(encoding="utf-8", errors="replace")
        if not coord_text.strip():
            continue
        if not first_coord_text:
            first_coord_text = coord_text
        target_coord = conf_dir / "coord"
        if not target_coord.is_file():
            write_text_file(target_coord, coord_text)

    root_coord = anmr_dir / "coord"
    if first_coord_text and not root_coord.is_file():
        write_text_file(root_coord, first_coord_text)


def _helper_launch_command(tool_path: str, extra_args: list[str] | None = None) -> list[str]:
    """Build a launcher command that respects shell-vs-Python helper scripts."""
    extra_args = list(extra_args or [])
    path = Path(tool_path).expanduser().resolve()
    try:
        header = path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return [str(path), *extra_args]

    exec_match = re.search(r'^exec\s+"([^"]+)"\s+"([^"]+)"\s+"\$@"', header, re.MULTILINE)
    if exec_match:
        wrapper_interpreter = exec_match.group(1)
        target_path = Path(exec_match.group(2)).expanduser()
        if target_path.is_file():
            return _helper_launch_command(str(target_path), extra_args=extra_args)
        return [wrapper_interpreter, str(target_path), *extra_args]

    first_line = header.splitlines()[0].strip() if header.splitlines() else ""
    shebang = first_line[2:].strip().lower() if first_line.startswith("#!") else ""
    if "bash" in shebang or shebang.endswith("/sh") or " sh" in shebang:
        return ["bash", str(path), *extra_args]
    if "python" in shebang:
        return [sys.executable, str(path), *extra_args]
    return [str(path), *extra_args]


def _parse_censo_version_text(text: str) -> tuple[int, int, int] | None:
    match = re.search(r"\b(\d+)\.(\d+)\.(\d+)\b", str(text or ""))
    if not match:
        return None
    return tuple(int(part) for part in match.groups())


def _detect_censo_version(censo_path: str) -> tuple[int, int, int] | None:
    try:
        result = subprocess.run(
            [censo_path, "-v"],
            capture_output=True,
            text=True,
            check=False,
        )
    except Exception:
        return None
    return _parse_censo_version_text((result.stdout or "") + "\n" + (result.stderr or ""))


def _build_censo_command(
    *,
    censo_path: str,
    ensemble_name: str,
    charge: int,
    multiplicity: int,
    rc_name: str,
    pal: int,
    solvent: str,
) -> list[str]:
    version = _detect_censo_version(censo_path)
    solvent_name = normalize_cpcm_solvent_name(solvent)
    base_cmd = [
        censo_path,
        "-i",
        ensemble_name,
        "-c",
        str(int(charge)),
        "-u",
        str(max(0, int(multiplicity) - 1)),
        "--inprc",
        rc_name,
        "--maxcores",
        str(max(1, int(pal))),
    ]
    if version and version[0] >= 3:
        base_cmd[1:1] = ["--prescreening", "--screening", "--optimization", "--nmr"]
        base_cmd.extend(["--omp-min", "1"])
        if solvent_name != "gas":
            base_cmd.extend(["--solvent", solvent_name])
        else:
            base_cmd.append("--gas-phase")
        return base_cmd

    base_cmd.extend(["-O", "1"])
    if solvent_name != "gas":
        base_cmd.extend(["-s", solvent_name])
    else:
        base_cmd.append("--gas-phase")
    return base_cmd


def _run_logged_shell(
    command: str,
    *,
    cwd: Path,
    log_path: Path,
    env: dict[str, str] | None = None,
) -> None:
    env_map = os.environ.copy()
    if env:
        env_map.update({str(k): str(v) for k, v in env.items() if str(v).strip()})
    result = subprocess.run(
        ["bash", "-lc", command],
        cwd=str(cwd),
        env=env_map,
        capture_output=True,
        text=True,
        check=False,
    )
    log_text = "$ bash -lc " + shlex.quote(command) + "\n\n" + result.stdout + result.stderr
    write_text_file(log_path, log_text)
    if result.returncode != 0:
        raise RuntimeError(f"Shell command failed ({result.returncode}): {command}")


def _guess_plot_window(solvent: str) -> tuple[float, float]:
    if str(solvent).strip().lower() == "h2o":
        return 0.0, 12.0
    return 0.0, 10.5


def _render_anmr_png(anmr_data_path: Path, output_png: Path, *, title: str) -> Path | None:
    """Render a simple headless PNG spectrum from `anmr.dat`."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    x_vals: list[float] = []
    y_vals: list[float] = []
    for line in anmr_data_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            x_vals.append(float(parts[0]))
            y_vals.append(float(parts[1]))
        except ValueError:
            continue
    if not x_vals:
        return None

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(x_vals, y_vals, color="black", linewidth=1.0)
    ax.fill_between(x_vals, y_vals, color="black", alpha=0.08)
    ax.set_xlabel("Chemical Shift (ppm)")
    ax.set_ylabel("Intensity")
    ax.set_title(title)
    ax.set_xlim(max(x_vals), min(x_vals))
    ax.margins(x=0.01, y=0.05)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_png if output_png.is_file() else None


def _build_censo_anmr_summary(payload: dict[str, object]) -> str:
    lines = [
        f"label: {payload.get('label', 'censo_anmr')}",
        "status: ok",
        f"source_xyz: {payload.get('source_xyz', '')}",
        f"workdir: {payload.get('workdir', '')}",
        f"solvent: {payload.get('solvent', '')}",
        f"charge: {payload.get('charge', 0)}",
        f"multiplicity: {payload.get('multiplicity', 1)}",
        f"resonance_frequency_mhz: {payload.get('resonance_frequency_mhz', 400.0)}",
        f"files: {json.dumps(payload.get('files', {}), indent=2)}",
    ]
    return "\n".join(lines)


def _truthy_env(name: str) -> bool:
    value = str(os.environ.get(name, "")).strip().lower()
    return value in {"1", "true", "yes", "on"}


def _file_contains(path: Path, needle: str) -> bool:
    if not path.is_file():
        return False
    try:
        return needle in path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False


def _censo_resume_requested(cli_flag: bool) -> bool:
    return bool(cli_flag) or _truthy_env("DELFIN_CENSO_NMR_RESUME")


def _crest_outputs_complete(workdir: Path) -> bool:
    required = [
        workdir / "crest_conformers.xyz",
        workdir / "anmr_nucinfo",
        workdir / "anmr_rotamer",
    ]
    return all(path.is_file() and path.stat().st_size > 0 for path in required)


def _censo_outputs_complete(workdir: Path) -> bool:
    return _file_contains(workdir / "censo.out", "CENSO all done!")


def _c2anmr_outputs_complete(workdir: Path) -> bool:
    anmr_dir = workdir / "anmr"
    return anmr_dir.is_dir() and (anmr_dir / "anmr_enso").is_file()


def _tms_reference_complete(workdir: Path, anmr_dir: Path) -> bool:
    return (
        (anmr_dir / ".anmrrc").is_file()
        and (
            _file_contains(workdir / "tms_reference.log", "ORCA TERMINATED NORMALLY")
            or (workdir / "tms_reference.out").is_file()
        )
    )


def _anmr_outputs_complete(anmr_dir: Path) -> bool:
    return (
        (anmr_dir / "anmr.dat").is_file()
        and _file_contains(anmr_dir / "anmr.out", "All done.")
    )


def _existing_anmr_plot(anmr_dir: Path) -> str:
    for name in ("anmr_spectrum.png", "anmr_spectrum.pdf", "anmr_spectrum.svg"):
        candidate = anmr_dir / name
        if candidate.is_file():
            return str(candidate)
    return ""


def _run_censo_anmr(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m delfin.dashboard.browser_workflows censo_anmr",
        description="Run a CREST + CENSO + ANMR ensemble NMR workflow in a dedicated workdir.",
    )
    parser.add_argument("--xyz-file", required=True, help="Source XYZ file.")
    parser.add_argument("--label", required=True, help="Workflow label.")
    parser.add_argument("--workdir", required=True, help="Exact work directory to write into.")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--solvent", default="chcl3")
    parser.add_argument("--pal", type=int, default=8)
    parser.add_argument("--maxcore", type=int, default=3000)
    parser.add_argument("--mhz", type=float, default=400.0)
    parser.add_argument("--resume", action="store_true", help="Reuse completed workflow stages in the existing workdir.")
    parser.add_argument("--json-out", help="Optional JSON summary file.")
    args = parser.parse_args(argv)

    workdir = resolve_path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    label = str(args.label).strip() or Path(args.xyz_file).stem or "censo_anmr"
    source_xyz = resolve_path(args.xyz_file)
    if not source_xyz.is_file():
        _write_error_summaries(
            workdir,
            label=label,
            error_message=f"XYZ file not found: {source_xyz}",
            json_out=args.json_out,
        )
        return 1

    try:
        tools = _require_tools(
            ["crest", "censo", "c2anmr", "anmr", "orca", "xtb"],
            workdir=workdir,
        )
        nmrplot_path = _resolved_path_or_empty("nmrplot")
        resume_mode = _censo_resume_requested(args.resume)

        xyz_target = workdir / source_xyz.name
        if source_xyz.resolve() != xyz_target.resolve():
            shutil.copy2(source_xyz, xyz_target)

        xyz_lines = xyz_target.read_text(encoding="utf-8", errors="replace").splitlines()
        xyz_body = "\n".join(line for line in xyz_lines[2:] if line.strip())
        if not xyz_body.strip():
            raise RuntimeError(f"XYZ file contains no coordinates: {xyz_target}")

        coord_path = write_text_file(workdir / "coord", xyz_body_to_coord_text(xyz_body))

        ensemble_path = workdir / "crest_conformers.xyz"
        if not (resume_mode and _crest_outputs_complete(workdir)):
            crest_cmd = [tools["crest"], str(coord_path)]
            if str(args.solvent).strip().lower() != "gas":
                crest_cmd.extend(["-g", str(args.solvent).strip().lower()])
            crest_cmd.extend(["-gfn2", "-T", str(max(1, int(args.pal))), "-nmr"])
            _run_logged(crest_cmd, cwd=workdir, log_path=workdir / "crest.out")
        if not ensemble_path.is_file():
            raise RuntimeError("CREST finished without writing crest_conformers.xyz")

        rc_path = write_text_file(
            workdir / "workflow.censo2rc",
            build_censo_anmr_rc(
                solvent=str(args.solvent).strip().lower(),
                resonance_frequency=float(args.mhz),
                active_nuclei=("h",),
                orca_path=tools["orca"],
                xtb_path=tools["xtb"],
                refinement_threshold=_censo_refinement_threshold(tools["censo"]),
            ),
        )

        censo_cmd = _build_censo_command(
            censo_path=tools["censo"],
            ensemble_name=str(ensemble_path.name),
            charge=int(args.charge),
            multiplicity=int(args.multiplicity),
            rc_name=str(rc_path.name),
            pal=max(1, int(args.pal)),
            solvent=str(args.solvent).strip().lower(),
        )
        if not (resume_mode and _censo_outputs_complete(workdir)):
            _run_logged(censo_cmd, cwd=workdir, log_path=workdir / "censo.out")
        _ensure_censo_nmr_coords(workdir)

        anmr_dir = workdir / "anmr"
        if not (resume_mode and _c2anmr_outputs_complete(workdir)):
            _run_logged(
                _helper_launch_command(tools["c2anmr"]),
                cwd=workdir,
                log_path=workdir / "c2anmr.out",
            )
        if not anmr_dir.is_dir():
            raise RuntimeError("c2anmr finished without creating the anmr/ folder")
        _ensure_anmr_coord_inputs(workdir=workdir, anmr_dir=anmr_dir)

        if not (resume_mode and _tms_reference_complete(workdir, anmr_dir)):
            tms_xyz, _num_atoms, _method, tms_error = smiles_to_xyz_quick("C[Si](C)(C)C")
            if tms_error or not tms_xyz:
                raise RuntimeError(f"Failed to generate TMS reference geometry: {tms_error or 'unknown error'}")
            tms_body = "\n".join(line for line in tms_xyz.splitlines() if line.strip())
            tms_inp = write_text_file(
                workdir / "tms_reference.inp",
                build_orca_reference_input(
                    tms_body,
                    solvent=str(args.solvent).strip().lower(),
                    pal=max(1, int(args.pal)),
                    maxcore=max(100, int(args.maxcore)),
                ),
            )
            _run_logged([tools["orca"], str(tms_inp.name)], cwd=workdir, log_path=workdir / "tms_reference.log")
        tms_out = _resolve_orca_reference_result_path(workdir=workdir, stem="tms_reference")

        tms_result = parse_nmr_orca(tms_out)
        h_vals = [item.isotropic_ppm for item in tms_result.shieldings if item.element == "H"]
        c_vals = [item.isotropic_ppm for item in tms_result.shieldings if item.element == "C"]
        if not h_vals or not c_vals:
            raise RuntimeError("Could not extract 1H/13C shieldings from the ORCA TMS reference run")

        write_text_file(
            anmr_dir / ".anmrrc",
            build_anmrrc_text(
                solvent=str(args.solvent).strip().lower(),
                resonance_frequency=float(args.mhz),
                shielding_ref_h=sum(h_vals) / len(h_vals),
                shielding_ref_c=sum(c_vals) / len(c_vals),
            ),
        )

        if not (resume_mode and _anmr_outputs_complete(anmr_dir)):
            anmr_cmd = f"ulimit -s unlimited && {shlex.quote(tools['anmr'])} -plain -mf {float(args.mhz):.1f}"
            _run_logged_shell(anmr_cmd, cwd=anmr_dir, log_path=anmr_dir / "anmr.out")

        plot_path = _existing_anmr_plot(anmr_dir) if resume_mode else ""
        if nmrplot_path and not plot_path:
            ppm_min, ppm_max = _guess_plot_window(str(args.solvent).strip().lower())
            plot_cmd = [
                *_helper_launch_command(nmrplot_path),
                "-i",
                "anmr.dat",
                "-start",
                str(ppm_min),
                "-end",
                str(ppm_max),
                "-o",
                "anmr_spectrum",
            ]
            _run_logged(
                plot_cmd,
                cwd=anmr_dir,
                log_path=anmr_dir / "nmrplot.out",
                env={"MPLBACKEND": "Agg"},
                input_text="n\n",
            )
            pdf_candidate = anmr_dir / "anmr_spectrum.pdf"
            png_candidate = anmr_dir / "anmr_spectrum.png"
            if pdf_candidate.is_file():
                plot_path = str(pdf_candidate)
            elif png_candidate.is_file():
                plot_path = str(png_candidate)

        rendered_png = _render_anmr_png(
            anmr_dir / "anmr.dat",
            anmr_dir / "anmr_spectrum.png",
            title=f"ANMR Spectrum: {label}",
        )
        if rendered_png is not None:
            plot_path = str(rendered_png)

        payload = {
            "label": label,
            "source_xyz": str(source_xyz),
            "workdir": str(workdir),
            "solvent": str(args.solvent).strip().lower(),
            "charge": int(args.charge),
            "multiplicity": int(args.multiplicity),
            "resonance_frequency_mhz": float(args.mhz),
            "files": {
                "crest_output": str(workdir / "crest.out"),
                "ensemble_xyz": str(ensemble_path),
                "censo_rc": str(rc_path),
                "censo_log": str(workdir / "censo.out"),
                "anmrrc": str(anmr_dir / ".anmrrc"),
                "anmr_output": str(anmr_dir / "anmr.out"),
                "anmr_data": str(anmr_dir / "anmr.dat"),
                "nmrplot_output": plot_path,
                "nmrplot_png": plot_path,
                "tms_reference_out": str(tms_out),
            },
            "errors": {},
        }
        _write_success_summaries(
            workdir,
            payload=payload,
            summary_text=_build_censo_anmr_summary(payload),
            json_out=args.json_out,
        )
        print(json.dumps(payload, indent=2))
        return 0
    except Exception as exc:  # noqa: BLE001
        _write_error_summaries(
            workdir,
            label=label,
            error_message=str(exc),
            json_out=args.json_out,
        )
        print(json.dumps({"result": None, "errors": {label: str(exc)}}, indent=2))
        return 1


def _run_hyperpol_xtb(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m delfin.dashboard.browser_workflows hyperpol_xtb",
        description="Run a single hyperpol_xtb XYZ workflow directly inside the chosen workdir.",
    )
    parser.add_argument("--xyz-file", required=True, help="XYZ file for the workflow.")
    parser.add_argument("--label", required=True, help="Workflow label.")
    parser.add_argument("--workdir", required=True, help="Exact work directory to write into.")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--engine", choices=("std2", "stda"), default="std2")
    parser.add_argument("--preopt", choices=("xtb", "crest", "none"), default="none")
    parser.add_argument("--static-only", action="store_true")
    parser.add_argument("--energy-window", type=float, default=15.0)
    parser.add_argument("--bfw", action="store_true", help="Pass -BFW to stda/std2.")
    parser.add_argument("--pal", type=int, default=4)
    parser.add_argument("--maxcore", type=int, default=1000)
    parser.add_argument("--json-out", help="Optional JSON summary file.")
    args = parser.parse_args(argv)

    label = str(args.label).strip() or Path(args.xyz_file).stem or "hyperpol_xtb"
    payload: dict[str, object]
    try:
        workdir = resolve_path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        wavelengths_nm = _hyperpol_load_wavelengths([], None, static_only=args.static_only)
        entry = HyperpolWorkflowEntry(label=label, xyz_path=str(resolve_path(args.xyz_file)))
        resolved_cores, resolved_maxcore = _resolve_active_job_limits(
            requested_cores=args.pal,
            requested_maxcore=getattr(args, "maxcore", 1000),
        )
        result = run_single_hyperpol_workflow(
            entry,
            charge=args.charge,
            multiplicity=args.multiplicity,
            engine=args.engine,
            preopt=args.preopt,
            wavelengths_nm=wavelengths_nm,
            energy_window_ev=args.energy_window,
            cores=resolved_cores,
            maxcore=resolved_maxcore,
            workdir=workdir,
            use_bfw=bool(args.bfw),
        )
        payload = {
            "summary": _hyperpol_summary_payload(result),
            "result": asdict(result),
            "errors": {},
        }
        _write_success_summaries(
            workdir,
            payload=payload,
            summary_text=_hyperpol_summary_text(result),
            json_out=args.json_out,
        )
        # Write DELFIN-style formatted report (same as pipeline)
        _write_text(workdir / "hyperpol_xtb_summary.txt", _hyperpol_formatted_report(result))
        # Write structured JSON (same as pipeline)
        _write_json(str(workdir / "hyperpol_xtb_summary.json"), {
            "hyperpol_xtb": {
                "settings": {
                    "engine": result.engine,
                    "preopt": result.preopt,
                    "charge": result.charge,
                    "multiplicity": result.multiplicity,
                    "energy_window_ev": result.energy_window_ev,
                    "requested_wavelengths_nm": result.requested_wavelengths_nm,
                },
                "dipole_moment": {
                    "x_au": result.dipole_x_au,
                    "y_au": result.dipole_y_au,
                    "z_au": result.dipole_z_au,
                    "total_debye": result.dipole_total_debye,
                },
                "response_points": result.response_points,
                "files": {
                    "initial_xyz": result.initial_xyz,
                    "selected_xyz": result.selected_xyz,
                    "xtb4stda_output": result.xtb4stda_output,
                    "response_output": result.response_output,
                    "beta_hrs_file": result.beta_hrs_file,
                    "beta_tensor_file": result.beta_tensor_file,
                },
            }
        })
        _print_hyperpol_result(result)
        return 0
    except Exception as exc:  # noqa: BLE001
        workdir = resolve_path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        payload = {"result": None, "errors": {label: str(exc)}}
        _write_error_summaries(
            workdir,
            label=label,
            error_message=str(exc),
            json_out=args.json_out,
        )
        print(json.dumps(payload, indent=2))
        return 1


def _run_tadf_xtb(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m delfin.dashboard.browser_workflows tadf_xtb",
        description="Run a single tadf_xtb XYZ workflow directly inside the chosen workdir.",
    )
    parser.add_argument("--xyz-file", required=True, help="XYZ file for the workflow.")
    parser.add_argument("--label", required=True, help="Workflow label.")
    parser.add_argument("--workdir", required=True, help="Exact work directory to write into.")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--xtb-method", default="XTB2")
    parser.add_argument("--preopt", choices=("none", "xtb", "crest", "goat"), default="none")
    parser.add_argument("--excited-method", choices=("stda", "std2"), default="stda")
    parser.add_argument("--energy-window", type=float, default=10.0)
    parser.add_argument("--crest", action="store_true")
    parser.add_argument("--goat", action="store_true")
    parser.add_argument("--t1-opt", action="store_true")
    parser.add_argument("--t1-multiplicity", type=int, default=3)
    parser.add_argument("--bfw", action="store_true", help="Pass -BFW to stda/std2.")
    parser.add_argument("--pal", type=int, default=4)
    parser.add_argument("--maxcore", type=int, default=1000)
    parser.add_argument("--json-out", help="Optional JSON summary file.")
    args = parser.parse_args(argv)

    label = str(args.label).strip() or Path(args.xyz_file).stem or "tadf_xtb"
    payload: dict[str, object]
    try:
        workdir = resolve_path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        entry = TadfWorkflowEntry(label=label, xyz_path=str(resolve_path(args.xyz_file)))
        use_crest = bool(args.crest or args.preopt == "crest")
        use_goat = bool(args.goat or args.preopt == "goat")
        resolved_cores, resolved_maxcore = _resolve_active_job_limits(
            requested_cores=args.pal,
            requested_maxcore=args.maxcore,
        )
        result = run_single_tadf_xtb(
            entry,
            charge=args.charge,
            multiplicity=args.multiplicity,
            xtb_method=args.xtb_method,
            excited_method=args.excited_method,
            energy_window=args.energy_window,
            cores=resolved_cores,
            maxcore=resolved_maxcore,
            workdir=workdir,
            use_crest=use_crest,
            use_goat=use_goat,
            run_t1_opt=bool(args.t1_opt),
            t1_multiplicity=args.t1_multiplicity,
            optimize_s0=(args.preopt == "xtb"),
            use_bfw=bool(args.bfw),
        )
        payload = {
            "summary": _tadf_summary_payload(result),
            "result": asdict(result),
            "errors": {},
        }
        _write_success_summaries(
            workdir,
            payload=payload,
            summary_text=_tadf_summary_text(result),
            json_out=args.json_out,
        )
        # Write DELFIN-style formatted report (same as pipeline)
        _preopt = "crest" if args.crest else ("goat" if args.goat else "none")
        _write_text(
            workdir / "tadf_xtb_summary.txt",
            _tadf_formatted_report(
                result,
                preopt=_preopt,
                energy_window=args.energy_window,
                run_t1_opt=bool(args.t1_opt),
            ),
        )
        # Write structured JSON (same as pipeline)
        def _opt_f(v, d=6):
            return round(float(v), d) if v is not None else None

        _write_json(str(workdir / "tadf_xtb_summary.json"), {
            "tadf_xtb": {
                "settings": {
                    "label": result.label,
                    "excited_method": result.excited_method,
                    "xtb_method": result.xtb_method,
                    "charge": result.charge,
                    "multiplicity": result.multiplicity,
                    "energy_window_ev": args.energy_window,
                    "preopt": _preopt,
                    "run_t1_opt": bool(args.t1_opt),
                },
                "excited_states": {
                    "s1_ev": _opt_f(result.s1_ev, 4),
                    "s1_nm": _opt_f(result.s1_nm, 1),
                    "s1_f_osc": _opt_f(result.s1_f_osc, 5),
                    "t1_ev": _opt_f(result.t1_ev, 4),
                    "t1_nm": _opt_f(result.t1_nm, 1),
                    "delta_est_vertical_ev": _opt_f(result.delta_est_ev, 4),
                    "first_allowed_singlet": {
                        "state": result.first_allowed_singlet_state,
                        "ev": _opt_f(result.first_allowed_singlet_ev, 4),
                        "nm": _opt_f(result.first_allowed_singlet_nm, 1),
                        "f_osc": _opt_f(result.first_allowed_singlet_f_osc, 5),
                    } if result.first_allowed_singlet_ev is not None else None,
                    "brightest_singlet": {
                        "state": result.brightest_singlet_state,
                        "ev": _opt_f(result.brightest_singlet_ev, 4),
                        "nm": _opt_f(result.brightest_singlet_nm, 1),
                        "f_osc": _opt_f(result.brightest_singlet_f_osc, 5),
                    } if result.brightest_singlet_ev is not None else None,
                },
                "energetics": {
                    "s0_xtb_energy_eh": _opt_f(result.s0_xtb_energy_eh, 8),
                    "t1_vertical_xtb_ev": _opt_f(result.t1_vertical_xtb_ev, 4),
                    "t1_adiabatic_ev": _opt_f(result.t1_adiabatic_ev, 4),
                    "t1_relaxed_ev": _opt_f(result.t1_relaxed_ev, 4),
                    "t1_relaxation_ev": _opt_f(result.t1_relaxation_ev, 4),
                    "s1_relaxed_est_ev": _opt_f(result.s1_relaxed_est_ev, 4),
                },
                "photophysics": {
                    "lambda_abs_nm": _opt_f(result.lambda_abs_nm, 1),
                    "lambda_em_s1_nm": _opt_f(result.lambda_em_s1_nm, 1),
                    "lambda_pl_est_nm": _opt_f(result.lambda_pl_est_nm, 1),
                    "stokes_shift_est_ev": _opt_f(result.stokes_shift_est_ev, 4),
                    "delta_est_relaxed_est_ev": _opt_f(result.delta_est_relaxed_est_ev, 4),
                },
                "files": {
                    "workdir": str(result.workdir),
                    "s0_xyz": result.s0_xyz,
                    "t1_opt_xyz": result.t1_opt_xyz,
                    "crest_xyz": result.crest_xyz,
                    "goat_xyz": result.goat_xyz,
                },
            }
        })
        _print_tadf_result(result)
        return 0
    except Exception as exc:  # noqa: BLE001
        workdir = resolve_path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        payload = {"result": None, "errors": {label: str(exc)}}
        _write_error_summaries(
            workdir,
            label=label,
            error_message=str(exc),
            json_out=args.json_out,
        )
        print(json.dumps(payload, indent=2))
        return 1


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("Usage: python -m delfin.dashboard.browser_workflows <hyperpol_xtb|tadf_xtb|censo_anmr> [...]")
        return 1

    command = args[0]
    if command == "hyperpol_xtb":
        return _run_hyperpol_xtb(args[1:])
    if command == "tadf_xtb":
        return _run_tadf_xtb(args[1:])
    if command == "censo_anmr":
        return _run_censo_anmr(args[1:])

    print(f"Unknown browser workflow: {command}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
