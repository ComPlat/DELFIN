"""Helpers for running single-file XYZ workflows launched from the browser UI."""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict
from pathlib import Path

from delfin.common.paths import resolve_path
from delfin.hyperpol import (
    WorkflowEntry as HyperpolWorkflowEntry,
    _load_wavelengths as _hyperpol_load_wavelengths,
    _print_result as _print_hyperpol_result,
    run_single_hyperpol_workflow,
)
from delfin.tadf_xtb import (
    WorkflowEntry as TadfWorkflowEntry,
    _print_result as _print_tadf_result,
    run_single_tadf_xtb,
)


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
        result = run_single_hyperpol_workflow(
            entry,
            charge=args.charge,
            multiplicity=args.multiplicity,
            engine=args.engine,
            preopt=args.preopt,
            wavelengths_nm=wavelengths_nm,
            energy_window_ev=args.energy_window,
            cores=args.pal,
            workdir=workdir,
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
    parser.add_argument("--excited-method", choices=("stda", "std2"), default="stda")
    parser.add_argument("--energy-window", type=float, default=10.0)
    parser.add_argument("--crest", action="store_true")
    parser.add_argument("--goat", action="store_true")
    parser.add_argument("--t1-opt", action="store_true")
    parser.add_argument("--t1-multiplicity", type=int, default=3)
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
        result = run_single_tadf_xtb(
            entry,
            charge=args.charge,
            multiplicity=args.multiplicity,
            xtb_method=args.xtb_method,
            excited_method=args.excited_method,
            energy_window=args.energy_window,
            cores=args.pal,
            maxcore=args.maxcore,
            workdir=workdir,
            use_crest=bool(args.crest),
            use_goat=bool(args.goat),
            run_t1_opt=bool(args.t1_opt),
            t1_multiplicity=args.t1_multiplicity,
            optimize_s0=False,
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
        print("Usage: python -m delfin.dashboard.browser_workflows <hyperpol_xtb|tadf_xtb> [...]")
        return 1

    command = args[0]
    if command == "hyperpol_xtb":
        return _run_hyperpol_xtb(args[1:])
    if command == "tadf_xtb":
        return _run_tadf_xtb(args[1:])

    print(f"Unknown browser workflow: {command}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
