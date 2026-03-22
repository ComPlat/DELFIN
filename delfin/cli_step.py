"""CLI entry point for ``delfin-step`` — run a single tool step.

Usage::

    delfin-step xtb_opt --geometry mol.xyz --charge 0 --cores 4
    delfin-step orca_sp --geometry mol.xyz --charge 0 --method B3LYP --basis def2-SVP
    delfin-step smiles_to_xyz --smiles "CCO"
    delfin-step --list              # show all available steps

Output is a one-line JSON summary printed to stdout.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="delfin-step",
        description="Run a single DELFIN tool step.",
    )
    parser.add_argument(
        "step_name",
        nargs="?",
        help="Registered step name (e.g. xtb_opt, orca_sp, smiles_to_xyz)",
    )
    parser.add_argument(
        "--list", "-l",
        action="store_true",
        dest="list_steps",
        help="List all available steps and exit",
    )
    parser.add_argument(
        "--geometry", "-g",
        type=str,
        default=None,
        help="Path to input XYZ file",
    )
    parser.add_argument(
        "--cores", "-j",
        type=str,
        default="1",
        help='CPU cores (default: 1). Use "auto" to detect from cluster/system.',
    )
    parser.add_argument(
        "--work-dir", "-d",
        type=str,
        default=None,
        help="Working directory (default: auto-created temp dir)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Full JSON output instead of summary line",
    )
    parser.add_argument(
        "--slurm",
        action="store_true",
        help="Submit as SLURM job instead of running locally",
    )
    return parser


def _parse_extra_kwargs(remaining: list[str]) -> dict:
    """Parse --key=value or --key value pairs into a dict.

    Tries to coerce values to int, then float, then leaves as str.
    Handles --flag (no value) as True.
    """
    kwargs = {}
    i = 0
    while i < len(remaining):
        arg = remaining[i]
        if not arg.startswith("--"):
            i += 1
            continue

        key = arg[2:].replace("-", "_")

        if "=" in key:
            key, val_str = key.split("=", 1)
        elif i + 1 < len(remaining) and not remaining[i + 1].startswith("--"):
            val_str = remaining[i + 1]
            i += 1
        else:
            # --flag with no value
            kwargs[key] = True
            i += 1
            continue

        # Type coercion
        kwargs[key] = _coerce(val_str)
        i += 1

    return kwargs


def _coerce(val: str):
    """Try int → float → str."""
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    # Boolean strings
    if val.lower() in ("true", "yes"):
        return True
    if val.lower() in ("false", "no"):
        return False
    return val


def _result_to_dict(result) -> dict:
    """Convert StepResult to a JSON-serializable dict."""
    return {
        "step_name": result.step_name,
        "status": result.status.value,
        "ok": result.ok,
        "geometry": str(result.geometry) if result.geometry else None,
        "output_file": str(result.output_file) if result.output_file else None,
        "work_dir": str(result.work_dir) if result.work_dir else None,
        "data": result.data,
        "artifacts": {k: str(v) for k, v in result.artifacts.items()},
        "error": result.error,
        "elapsed_seconds": round(result.elapsed_seconds, 2),
    }


def _resolve_cores_arg(cores_str: str) -> int:
    """Resolve --cores argument: integer or 'auto'."""
    if cores_str.lower() == "auto":
        from delfin.cluster_utils import detect_cluster_environment
        info = detect_cluster_environment()
        detected = info.get("cpus_available") or 1
        print(f"Auto-detected {detected} cores ({info.get('scheduler', 'system')})")
        return detected
    return int(cores_str)


def _submit_step_slurm(args, cores: int, extra: dict) -> int:
    """Submit a single step as a SLURM job."""
    import shutil
    import subprocess

    sbatch = shutil.which("sbatch")
    if not sbatch:
        print("Error: sbatch not found. Is SLURM installed?", file=sys.stderr)
        return 1

    # Build the command to run inside SLURM
    cmd_parts = ["delfin-step", args.step_name, f"--cores={cores}"]
    if args.geometry:
        cmd_parts.append(f"--geometry={args.geometry}")
    if args.work_dir:
        cmd_parts.append(f"--work-dir={args.work_dir}")
    for k, v in extra.items():
        cmd_parts.append(f"--{k}={v}")

    step_cmd = " ".join(cmd_parts)
    job_name = f"delfin-{args.step_name}"

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cores}
#SBATCH --time=24:00:00
#SBATCH --output={job_name}_%j.out
#SBATCH --error={job_name}_%j.err

echo "DELFIN step: {args.step_name}"
echo "Cores: {cores}"
echo "Start: $(date)"

{step_cmd}

echo "End: $(date)"
"""
    script_path = Path(args.work_dir or ".") / f"{job_name}_slurm.sh"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(script)

    result = subprocess.run([sbatch, str(script_path)], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"SLURM job submitted: {result.stdout.strip()}")
        print(f"Script: {script_path}")
        return 0
    else:
        print(f"SLURM submission failed: {result.stderr}", file=sys.stderr)
        return 1


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args, remaining = parser.parse_known_args(argv)

    # --list: show available steps
    if args.list_steps:
        from delfin.tools import list_steps
        steps = list_steps()
        if not steps:
            print("No steps registered.")
            return 0
        max_name = max(len(n) for n in steps)
        for name, adapter in sorted(steps.items()):
            geom = "→ XYZ" if adapter.produces_geometry else "      "
            print(f"  {name:<{max_name}}  {geom}  {adapter.description}")
        return 0

    if not args.step_name:
        parser.print_help()
        return 1

    # Parse extra kwargs from remaining args
    extra = _parse_extra_kwargs(remaining)

    # Resolve cores (supports "auto" for cluster detection)
    cores = _resolve_cores_arg(args.cores)

    # SLURM submission
    if args.slurm:
        return _submit_step_slurm(args, cores, extra)

    from delfin.tools import run_step

    try:
        result = run_step(
            args.step_name,
            geometry=args.geometry,
            cores=cores,
            work_dir=Path(args.work_dir) if args.work_dir else None,
            **extra,
        )
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    if args.json_output:
        print(json.dumps(_result_to_dict(result), indent=2))
    else:
        status = "OK" if result.ok else "FAILED"
        elapsed = f"{result.elapsed_seconds:.1f}s"
        parts = [f"[{status}] {result.step_name} ({elapsed})"]
        if result.geometry:
            parts.append(f"geometry={result.geometry}")
        if result.data:
            for k, v in result.data.items():
                if isinstance(v, float):
                    parts.append(f"{k}={v:.6f}")
                else:
                    parts.append(f"{k}={v}")
        if result.error:
            parts.append(f"error={result.error}")
        print("  ".join(parts))

    return 0 if result.ok else 1


if __name__ == "__main__":
    sys.exit(main())
