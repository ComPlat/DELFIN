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
        type=str,
        default="1",
        help='CPU cores (default: 1). Use "auto" to detect from cluster/system.',
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


def _build_condition(expr: str) -> Any:
    """Build a condition callable from a YAML expression string.

    Supported expressions:
        - ``last.data.key == value``
        - ``last.data.key > value``
        - ``last.data.key < value``
        - ``last.data.key >= value``
        - ``last.data.key != value``
        - ``last.ok``
        - ``not last.ok``
    """
    import re

    expr = expr.strip()

    # "not last.ok"
    if expr == "not last.ok":
        return lambda results, last: last is not None and not last.ok

    # "last.ok"
    if expr == "last.ok":
        return lambda results, last: last is not None and last.ok

    # "last.data.KEY OP VALUE"
    m = re.match(r"last\.data\.(\w+)\s*(==|!=|>=|<=|>|<)\s*(.+)", expr)
    if m:
        key, op, raw_val = m.group(1), m.group(2), m.group(3).strip()
        # Parse value
        try:
            val = int(raw_val)
        except ValueError:
            try:
                val = float(raw_val)
            except ValueError:
                val = raw_val.strip("'\"")

        ops = {
            "==": lambda a, b: a == b,
            "!=": lambda a, b: a != b,
            ">": lambda a, b: a > b,
            "<": lambda a, b: a < b,
            ">=": lambda a, b: a >= b,
            "<=": lambda a, b: a <= b,
        }
        op_fn = ops[op]
        return lambda results, last, _k=key, _fn=op_fn, _v=val: (
            last is not None and _fn(last.data.get(_k), _v)
        )

    raise ValueError(f"Cannot parse condition expression: {expr!r}")


def _build_until(expr: str) -> Any:
    """Build a loop-until callable from a YAML expression string.

    Supported: ``result.data.KEY OP VALUE``
    """
    import re

    m = re.match(r"result\.data\.(\w+)\s*(==|!=|>=|<=|>|<)\s*(.+)", expr.strip())
    if m:
        key, op, raw_val = m.group(1), m.group(2), m.group(3).strip()
        try:
            val = int(raw_val)
        except ValueError:
            try:
                val = float(raw_val)
            except ValueError:
                val = raw_val.strip("'\"")
        ops = {
            "==": lambda a, b: a == b,
            "!=": lambda a, b: a != b,
            ">": lambda a, b: a > b,
            "<": lambda a, b: a < b,
            ">=": lambda a, b: a >= b,
            "<=": lambda a, b: a <= b,
        }
        op_fn = ops[op]
        return lambda result, iteration, _k=key, _fn=op_fn, _v=val: (
            _fn(result.data.get(_k), _v)
        )

    raise ValueError(f"Cannot parse until expression: {expr!r}")


def _add_step_to_pipeline(pipe, spec: dict) -> None:
    """Add a step (possibly with flow control) to a Pipeline.

    YAML step types:
        - Normal: ``{step: xtb_opt, charge: 0}``
        - Conditional: ``{step: imag_fix, type: if, condition: "last.data.n_imaginary > 0"}``
        - Loop: ``{step: imag_fix, type: loop, until: "result.data.n_imaginary == 0", max_iter: 5}``
        - Retry: ``{step: orca_sp, type: retry, max_attempts: 3}``
        - Checkpoint: ``{type: checkpoint}``
        - Compute: ``{type: compute, label: "calc_energy", module: "my_module", function: "calc"}``
    """
    spec = dict(spec)  # don't mutate original
    step_type = spec.pop("type", "normal")

    if step_type == "checkpoint":
        label = spec.pop("label", "checkpoint")
        pipe.add_checkpoint(label=label)
        return

    if step_type == "compute":
        label = spec.pop("label", "compute")
        module_name = spec.pop("module", None)
        func_name = spec.pop("function", None)
        if not module_name or not func_name:
            raise ValueError("compute step requires 'module' and 'function'")
        import importlib
        mod = importlib.import_module(module_name)
        fn = getattr(mod, func_name)
        pipe.add_compute(fn, label=label)
        return

    if "step" not in spec:
        raise ValueError(f"Step must have a 'step' key, got: {spec}")

    step_name = spec.pop("step")
    label = spec.pop("label", "")

    if step_type == "if":
        condition_expr = spec.pop("condition", None)
        if not condition_expr:
            raise ValueError("'if' step requires 'condition'")
        cond_fn = _build_condition(condition_expr)
        pipe.add_if(cond_fn, step_name, label=label, **spec)

    elif step_type == "loop":
        until_expr = spec.pop("until", None)
        max_iter = spec.pop("max_iter", 10)
        if not until_expr:
            raise ValueError("'loop' step requires 'until'")
        until_fn = _build_until(until_expr)
        pipe.add_loop(step_name, until=until_fn, max_iter=max_iter, label=label, **spec)

    elif step_type == "retry":
        max_attempts = spec.pop("max_attempts", 3)
        delay = spec.pop("delay", 0.0)
        pipe.add_retry(step_name, max_attempts=max_attempts, delay=delay, label=label, **spec)

    elif step_type == "normal" or step_type == "step":
        pipe.add(step_name, label=label, **spec)

    else:
        raise ValueError(f"Unknown step type: {step_type!r}")


def _build_pipeline_from_yaml(data: Dict[str, Any]):
    """Build a Pipeline or PipelineTemplate from parsed YAML.

    Supports flow control types in YAML::

        steps:
          - step: xtb_opt
            charge: 0
          - step: orca_freq
            charge: 0
          - step: imag_fix
            type: if
            condition: "last.data.n_imaginary > 0"
            charge: 0
          - step: orca_sp
            type: retry
            max_attempts: 3
            charge: 0
          - type: checkpoint
          - step: orca_opt
            type: loop
            until: "result.data.energy_change < 0.0001"
            max_iter: 10
            charge: 0
    """
    is_template = data.get("template", False)
    name = data["name"]
    defaults = data.get("defaults", {})
    steps = data["steps"]
    branches = data.get("branches", {})

    if is_template:
        from delfin.tools.pipeline import PipelineTemplate
        pipe = PipelineTemplate(name, defaults=defaults)
        # Templates only support basic add() — no flow control
        for spec in steps:
            if not isinstance(spec, dict) or "step" not in spec:
                raise ValueError(f"Each step must be a dict with 'step' key, got: {spec}")
            spec = dict(spec)
            step_name = spec.pop("step")
            spec.pop("type", None)  # ignore type for templates
            label = spec.pop("label", "")
            pipe.add(step_name, label=label, **spec)
    else:
        from delfin.tools.pipeline import Pipeline
        pipe = Pipeline(name, defaults=defaults)
        for spec in steps:
            if not isinstance(spec, dict):
                raise ValueError(f"Each step must be a dict, got: {spec}")
            _add_step_to_pipeline(pipe, spec)

    for bname, bsteps in branches.items():
        branch = pipe.branch(bname)
        for spec in bsteps:
            if not isinstance(spec, dict):
                raise ValueError(f"Branch step must be a dict, got: {spec}")
            if is_template:
                if "step" not in spec:
                    raise ValueError(f"Branch step must have 'step' key, got: {spec}")
                spec = dict(spec)
                step_name = spec.pop("step")
                spec.pop("type", None)
                label = spec.pop("label", "")
                branch.add(step_name, label=label, **spec)
            else:
                _add_step_to_pipeline(branch, spec)

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


def _resolve_cores_arg(cores_str: str) -> int:
    """Resolve --cores argument: integer or 'auto'."""
    if cores_str.lower() == "auto":
        from delfin.cluster_utils import detect_cluster_environment
        info = detect_cluster_environment()
        detected = info.get("cpus_available") or 1
        print(f"Auto-detected {detected} cores ({info.get('scheduler', 'system')})")
        return detected
    return int(cores_str)


def _submit_slurm(yaml_path: str, args) -> int:
    """Submit the pipeline as a SLURM job."""
    import shutil
    import subprocess

    sbatch = shutil.which("sbatch")
    if not sbatch:
        print("Error: sbatch not found. Is SLURM installed?", file=sys.stderr)
        return 1

    cores = _resolve_cores_arg(args.cores)

    # Inside SLURM job, use "auto" to detect allocated resources
    cmd_parts = ["delfin-pipeline", yaml_path, "--cores=auto"]
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
#SBATCH --cpus-per-task={cores}
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
echo "Cores: {cores} (auto-detect inside job)"
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
