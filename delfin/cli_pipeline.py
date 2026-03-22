"""CLI entry point for ``delfin-pipeline`` — run a pipeline from YAML.

Usage::

    delfin-pipeline workflow.yaml --cores 8
    delfin-pipeline workflow.yaml --cores 8 --geometry input.xyz
    delfin-pipeline workflow.yaml --cores 8 --param charge=2 --param method=B3LYP
    delfin-pipeline workflow.yaml --cores 8 --scheduled   # use DynamicCorePool
    delfin-pipeline workflow.yaml --cores 8 --slurm       # submit via SLURM

YAML format::

    name: opt_freq
    defaults:
      charge: 0
      method: B3LYP
      basis: def2-SVP
    steps:
      - step: smiles_to_xyz
        smiles: "CCO"
      - step: xtb_opt
      - step: orca_opt
        ri: RIJCOSX
        aux_basis: def2/J
      - step: orca_freq

Parametric template (with {placeholders})::

    name: redox_template
    template: true
    defaults:
      method: "{method}"
      basis: "{basis}"
    steps:
      - step: xtb_opt
        charge: "{charge}"
      - step: orca_opt
        charge: "{charge}"
    branches:
      oxidation:
        - step: orca_opt
          charge: "{charge}+1"
          mult: "{mult_ox}"
      reduction:
        - step: orca_opt
          charge: "{charge}-1"
          mult: "{mult_red}"
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="delfin-pipeline",
        description="Run a DELFIN pipeline from a YAML definition.",
    )
    parser.add_argument(
        "yaml_file",
        type=str,
        help="Path to pipeline YAML file",
    )
    parser.add_argument(
        "--cores", "-j",
        type=int,
        default=1,
        help="Number of CPU cores (default: 1)",
    )
    parser.add_argument(
        "--geometry", "-g",
        type=str,
        default=None,
        help="Initial input geometry (XYZ file)",
    )
    parser.add_argument(
        "--work-dir", "-d",
        type=str,
        default=None,
        help="Root working directory (default: current dir)",
    )
    parser.add_argument(
        "--param", "-p",
        action="append",
        default=[],
        help="Template parameter as key=value (repeatable)",
    )
    parser.add_argument(
        "--scheduled",
        action="store_true",
        help="Use DynamicCorePool scheduler instead of simple sequential execution",
    )
    parser.add_argument(
        "--slurm",
        action="store_true",
        help="Submit as SLURM job (requires SLURM environment)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Output results as JSON",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse and validate the YAML without executing",
    )
    return parser


def _parse_params(param_list: List[str]) -> Dict[str, Any]:
    """Parse key=value param strings into a dict with type coercion."""
    params = {}
    for item in param_list:
        if "=" not in item:
            print(f"Warning: ignoring param without '=': {item}", file=sys.stderr)
            continue
        key, val = item.split("=", 1)
        try:
            params[key] = int(val)
        except ValueError:
            try:
                params[key] = float(val)
            except ValueError:
                params[key] = val
    return params


def _load_yaml(path: str) -> Dict[str, Any]:
    """Load and validate a pipeline YAML file."""
    import yaml

    yaml_path = Path(path)
    if not yaml_path.is_file():
        raise FileNotFoundError(f"Pipeline file not found: {yaml_path}")

    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict):
        raise ValueError(f"Pipeline YAML must be a mapping, got {type(data).__name__}")

    if "name" not in data:
        raise ValueError("Pipeline YAML must have a 'name' field")

    if "steps" not in data:
        raise ValueError("Pipeline YAML must have a 'steps' list")

    return data


def _build_pipeline_from_yaml(data: Dict[str, Any]):
    """Build a Pipeline or PipelineTemplate from parsed YAML."""
    is_template = data.get("template", False)
    name = data["name"]
    defaults = data.get("defaults", {})
    steps = data["steps"]
    branches = data.get("branches", {})

    if is_template:
        from delfin.tools.pipeline import PipelineTemplate
        pipe = PipelineTemplate(name, defaults=defaults)
    else:
        from delfin.tools.pipeline import Pipeline
        pipe = Pipeline(name, defaults=defaults)

    for spec in steps:
        if not isinstance(spec, dict) or "step" not in spec:
            raise ValueError(f"Each step must be a dict with 'step' key, got: {spec}")
        step_name = spec.pop("step")
        label = spec.pop("label", "")
        pipe.add(step_name, label=label, **spec)

    for bname, bsteps in branches.items():
        branch = pipe.branch(bname)
        for spec in bsteps:
            if not isinstance(spec, dict) or "step" not in spec:
                raise ValueError(f"Branch step must be a dict with 'step' key, got: {spec}")
            step_name = spec.pop("step")
            label = spec.pop("label", "")
            branch.add(step_name, label=label, **spec)

    return pipe


def _result_to_dict(result) -> dict:
    """Convert PipelineResult to JSON-serializable dict."""
    d = {
        "name": result.name,
        "ok": result.ok,
        "steps": [],
        "branches": {},
    }
    for r in result.results:
        d["steps"].append({
            "step_name": r.step_name,
            "status": r.status.value,
            "ok": r.ok,
            "geometry": str(r.geometry) if r.geometry else None,
            "data": r.data,
            "error": r.error,
            "elapsed_seconds": round(r.elapsed_seconds, 2),
        })
    for bname, br in result.branch_results.items():
        d["branches"][bname] = _result_to_dict(br)
    return d


def _submit_slurm(yaml_path: str, args) -> int:
    """Submit the pipeline as a SLURM job."""
    import shutil
    import subprocess

    sbatch = shutil.which("sbatch")
    if not sbatch:
        print("Error: sbatch not found. Is SLURM installed?", file=sys.stderr)
        return 1

    # Build the delfin-pipeline command to run inside the SLURM job
    cmd_parts = ["delfin-pipeline", yaml_path, f"--cores={args.cores}"]
    if args.geometry:
        cmd_parts.append(f"--geometry={args.geometry}")
    if args.work_dir:
        cmd_parts.append(f"--work-dir={args.work_dir}")
    for p in args.param:
        cmd_parts.append(f"--param={p}")

    pipeline_cmd = " ".join(cmd_parts)

    # Load YAML to get pipeline name for SLURM job name
    data = _load_yaml(yaml_path)
    job_name = data.get("name", "delfin-pipeline")

    # SLURM directives
    slurm_config = data.get("slurm", {})
    partition = slurm_config.get("partition", "")
    time_limit = slurm_config.get("time", "24:00:00")
    mem = slurm_config.get("mem", "4G")
    nodes = slurm_config.get("nodes", 1)

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.cores}
#SBATCH --mem={mem}
#SBATCH --time={time_limit}
#SBATCH --nodes={nodes}
"""
    if partition:
        script += f"#SBATCH --partition={partition}\n"

    script += f"""
#SBATCH --output={job_name}_%j.out
#SBATCH --error={job_name}_%j.err

echo "DELFIN pipeline: {job_name}"
echo "Cores: {args.cores}"
echo "Start: $(date)"

{pipeline_cmd}

echo "End: $(date)"
"""

    # Write script
    script_path = Path(args.work_dir or ".") / f"{job_name}_slurm.sh"
    script_path.write_text(script)

    # Submit
    result = subprocess.run(
        [sbatch, str(script_path)],
        capture_output=True, text=True,
    )

    if result.returncode == 0:
        print(f"SLURM job submitted: {result.stdout.strip()}")
        print(f"Script: {script_path}")
        return 0
    else:
        print(f"SLURM submission failed: {result.stderr}", file=sys.stderr)
        return 1


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Load YAML
    try:
        data = _load_yaml(args.yaml_file)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Build pipeline
    try:
        pipe = _build_pipeline_from_yaml(data)
    except (ValueError, KeyError) as e:
        print(f"Error building pipeline: {e}", file=sys.stderr)
        return 1

    # Dry run
    if args.dry_run:
        print(f"Pipeline: {pipe}")
        is_tpl = data.get("template", False)
        n_steps = len(data["steps"])
        n_branches = len(data.get("branches", {}))
        print(f"  Type: {'template' if is_tpl else 'concrete'}")
        print(f"  Steps: {n_steps}")
        print(f"  Branches: {n_branches}")
        if data.get("defaults"):
            print(f"  Defaults: {data['defaults']}")
        print("  Dry run — no execution.")
        return 0

    # SLURM submission
    if args.slurm:
        return _submit_slurm(args.yaml_file, args)

    # Parse template params
    params = _parse_params(args.param)

    work_dir = Path(args.work_dir) if args.work_dir else None

    # Execute
    from delfin.tools.pipeline import PipelineTemplate
    is_template = isinstance(pipe, PipelineTemplate)

    try:
        if is_template:
            if args.scheduled:
                concrete = pipe.build(**params)
                result = concrete.run_scheduled(
                    cores=args.cores,
                    geometry=args.geometry,
                    work_dir=work_dir,
                )
            else:
                result = pipe.run(
                    cores=args.cores,
                    geometry=args.geometry,
                    work_dir=work_dir,
                    **params,
                )
        else:
            if args.scheduled:
                result = pipe.run_scheduled(
                    cores=args.cores,
                    geometry=args.geometry,
                    work_dir=work_dir,
                )
            else:
                result = pipe.run(
                    cores=args.cores,
                    geometry=args.geometry,
                    work_dir=work_dir,
                )
    except Exception as e:
        print(f"Pipeline execution error: {e}", file=sys.stderr)
        return 1

    # Output
    if args.json_output:
        print(json.dumps(_result_to_dict(result), indent=2))
    else:
        print(result.summary())

    return 0 if result.ok else 1


if __name__ == "__main__":
    sys.exit(main())
