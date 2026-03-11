from __future__ import annotations

import argparse
import json
import re
import subprocess
import threading
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path
from delfin.dynamic_pool import JobPriority, PoolJob
from delfin.global_manager import get_global_manager
from delfin.orca import run_orca
from delfin.qm_runtime import run_tool
from delfin.smiles_converter import smiles_to_xyz_quick

logger = get_logger(__name__)

_STATE_LINE_RE = re.compile(
    r"^\s*(\d+)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)"
)
_LABEL_SAFE_RE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass(frozen=True)
class WorkflowEntry:
    label: str
    smiles: str


@dataclass
class TadfXtbResult:
    label: str
    smiles: str
    charge: int
    multiplicity: int
    xtb_method: str
    workdir: str
    goat_xyz: str
    singlet_output: str
    triplet_output: str
    tda_singlet: Optional[str]
    tda_triplet: Optional[str]
    s1_ev: Optional[float]
    s1_nm: Optional[float]
    s1_f_osc: Optional[float]
    t1_ev: Optional[float]
    t1_nm: Optional[float]
    delta_est_ev: Optional[float]


def _safe_label(label: str, fallback: str) -> str:
    text = _LABEL_SAFE_RE.sub("_", (label or "").strip()).strip("._-")
    return text or fallback


def _write_smiles_inputs(smiles: str, label: str, workdir: Path) -> tuple[Path, Path]:
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


def _write_goat_input(
    start_path: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    xtb_method: str,
    cores: int,
    maxcore: int,
    skip_initial_opt: bool,
) -> Path:
    coord_block = start_path.read_text(encoding="utf-8")
    goat_block = ""
    if skip_initial_opt:
        goat_block = "%goat\n  SKIPINITIALOPT TRUE\nend\n"
    input_path = workdir / "XTB_GOAT.inp"
    input_path.write_text(
        f"!{xtb_method} GOAT\n"
        f"%pal nprocs {max(1, int(cores))} end\n"
        f"%maxcore {max(1, int(maxcore))}\n"
        f"{goat_block}"
        f"*xyz {charge} {multiplicity}\n"
        f"{coord_block}"
        "*\n",
        encoding="utf-8",
    )
    return input_path


def _run_goat(
    start_path: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    xtb_method: str,
    cores: int,
    maxcore: int,
    skip_initial_opt: bool,
) -> Path:
    input_path = _write_goat_input(
        start_path,
        workdir,
        charge=charge,
        multiplicity=multiplicity,
        xtb_method=xtb_method,
        cores=cores,
        maxcore=maxcore,
        skip_initial_opt=skip_initial_opt,
    )
    output_path = workdir / "output_XTB_GOAT.out"
    ok = run_orca(str(input_path), str(output_path), working_dir=workdir)
    if not ok:
        raise RuntimeError(f"GOAT failed; see {output_path}")
    goat_xyz = workdir / "XTB_GOAT.globalminimum.xyz"
    if not goat_xyz.exists():
        raise RuntimeError(f"GOAT did not produce {goat_xyz.name}")
    return goat_xyz


def _run_xtb4stda(
    goat_xyz: Path,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
) -> Path:
    uhf = max(0, int(multiplicity) - 1)
    output_path = workdir / "xtb4stda.out"
    with output_path.open("w", encoding="utf-8") as log_file:
        run_tool(
            "xtb4stda",
            [goat_xyz.name, "-chrg", str(charge), "-uhf", str(uhf)],
            cwd=workdir,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
            track_process=True,
        )
    wfn_path = workdir / "wfn.xtb"
    if not wfn_path.exists():
        raise RuntimeError(f"xtb4stda did not produce {wfn_path.name}")
    return output_path


def _run_stda(workdir: Path, *, triplet: bool, energy_window: float) -> tuple[Path, Optional[Path]]:
    output_path = workdir / ("stda_triplet.out" if triplet else "stda_singlet.out")
    args = ["-xtb", "-e", str(energy_window)]
    if triplet:
        args.insert(1, "-t")
    with output_path.open("w", encoding="utf-8") as log_file:
        run_tool(
            "stda",
            args,
            cwd=workdir,
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


def run_single_tadf_xtb(
    entry: WorkflowEntry,
    *,
    charge: int,
    multiplicity: int,
    xtb_method: str,
    energy_window: float,
    cores: int,
    maxcore: int,
    goat_skip_initial_opt: bool,
    workdir: Path,
) -> TadfXtbResult:
    workdir.mkdir(parents=True, exist_ok=True)
    start_path, _ = _write_smiles_inputs(entry.smiles, entry.label, workdir)
    goat_xyz = _run_goat(
        start_path,
        workdir,
        charge=charge,
        multiplicity=multiplicity,
        xtb_method=xtb_method,
        cores=cores,
        maxcore=maxcore,
        skip_initial_opt=goat_skip_initial_opt,
    )
    _run_xtb4stda(
        goat_xyz,
        workdir,
        charge=charge,
        multiplicity=multiplicity,
    )
    singlet_out, tda_singlet = _run_stda(workdir, triplet=False, energy_window=energy_window)
    triplet_out, tda_triplet = _run_stda(workdir, triplet=True, energy_window=energy_window)

    singlets = _parse_stda_states(singlet_out)
    triplets = _parse_stda_states(triplet_out)
    s1 = singlets[0] if singlets else None
    t1 = triplets[0] if triplets else None

    return TadfXtbResult(
        label=entry.label,
        smiles=entry.smiles,
        charge=charge,
        multiplicity=multiplicity,
        xtb_method=xtb_method,
        workdir=str(workdir),
        goat_xyz=str(goat_xyz),
        singlet_output=str(singlet_out),
        triplet_output=str(triplet_out),
        tda_singlet=str(tda_singlet) if tda_singlet else None,
        tda_triplet=str(tda_triplet) if tda_triplet else None,
        s1_ev=float(s1["ev"]) if s1 else None,
        s1_nm=float(s1["nm"]) if s1 else None,
        s1_f_osc=float(s1["f_osc"]) if s1 else None,
        t1_ev=float(t1["ev"]) if t1 else None,
        t1_nm=float(t1["nm"]) if t1 else None,
        delta_est_ev=(float(s1["ev"]) - float(t1["ev"])) if s1 and t1 else None,
    )


def _load_entries(smiles_values: list[str], smiles_file: Optional[str], label: Optional[str]) -> list[WorkflowEntry]:
    entries: list[WorkflowEntry] = []
    for idx, smiles in enumerate(smiles_values, start=1):
        fallback = f"mol_{idx:03d}"
        entries.append(WorkflowEntry(label=_safe_label(label or fallback, fallback), smiles=smiles.strip()))

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
        unique_entries.append(WorkflowEntry(label=f"{entry.label}_{count:02d}", smiles=entry.smiles))
    return unique_entries


def _print_result(result: TadfXtbResult) -> None:
    def _fmt(value: Optional[float], digits: int) -> str:
        if value is None:
            return "n/a"
        return f"{value:.{digits}f}"

    print(
        f"{result.label}: "
        f"S1={_fmt(result.s1_ev, 3)} eV, "
        f"T1={_fmt(result.t1_ev, 3)} eV, "
        f"DeltaEST={_fmt(result.delta_est_ev, 3)} eV, "
        f"f_osc={_fmt(result.s1_f_osc, 4)}, "
        f"workdir={result.workdir}"
    )


def run_cli(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin tadf_xtb",
        description="SMILES -> GOAT -> xtb4stda -> stda(singlet/triplet) TADF screening workflow.",
    )
    parser.add_argument("--smiles", action="append", default=[], help="SMILES string. Repeatable.")
    parser.add_argument("--smiles-file", help="Optional file with one SMILES per line, optionally prefixed by label.")
    parser.add_argument("--label", help="Optional label for a single --smiles input.")
    parser.add_argument("--workdir", default="tadf_xtb_runs", help="Base output directory.")
    parser.add_argument("--charge", type=int, default=0, help="Molecular charge.")
    parser.add_argument("--multiplicity", type=int, default=1, help="Ground-state multiplicity.")
    parser.add_argument("--xtb-method", default="XTB2", help="ORCA xTB keyword for GOAT, e.g. XTB2.")
    parser.add_argument("--energy-window", type=float, default=10.0, help="Excitation window in eV for stda.")
    parser.add_argument("--pal", type=int, default=4, help="Cores per molecule workflow.")
    parser.add_argument("--parallel-jobs", type=int, default=1, help="Number of molecule workflows to run in parallel.")
    parser.add_argument("--maxcore", type=int, default=1000, help="Memory per core in MB for pool scheduling.")
    parser.add_argument(
        "--goat-skip-initial-opt",
        action="store_true",
        help="Pass SKIPINITIALOPT TRUE to ORCA GOAT to speed up exploratory runs.",
    )
    parser.add_argument("--json-out", help="Optional JSON summary path.")
    args = parser.parse_args(argv)

    entries = _load_entries(args.smiles, args.smiles_file, args.label)
    if not entries:
        parser.error("provide --smiles and/or --smiles-file")

    base_workdir = resolve_path(args.workdir)
    base_workdir.mkdir(parents=True, exist_ok=True)

    results: list[TadfXtbResult] = []
    errors: dict[str, str] = {}

    if len(entries) == 1 or args.parallel_jobs <= 1:
        for entry in entries:
            try:
                result = run_single_tadf_xtb(
                    entry,
                    charge=args.charge,
                    multiplicity=args.multiplicity,
                    xtb_method=args.xtb_method,
                    energy_window=args.energy_window,
                    cores=args.pal,
                    maxcore=args.maxcore,
                    goat_skip_initial_opt=args.goat_skip_initial_opt,
                    workdir=base_workdir / entry.label,
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
                    result = run_single_tadf_xtb(
                        _entry,
                        charge=args.charge,
                        multiplicity=args.multiplicity,
                        xtb_method=args.xtb_method,
                        energy_window=args.energy_window,
                        cores=max(1, int(allocated)),
                        maxcore=args.maxcore,
                        goat_skip_initial_opt=args.goat_skip_initial_opt,
                        workdir=_workdir,
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
