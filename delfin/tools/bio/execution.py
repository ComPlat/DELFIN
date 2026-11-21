"""Explicit DAG bindings onto DELFIN's existing run_step execution path."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import importlib.util
import importlib.metadata
import inspect
import json
from pathlib import Path
import subprocess

from delfin.tools.bio.artifacts import child_path, file_hash, read_json, validate_artifact, write_json


@dataclass(frozen=True)
class Binding:
    adapter: str
    inputs: dict[str, str]
    outputs: dict[str, str]


BINDINGS = {
    "ligand_parameterize": Binding("bio_ligand_parameterize", {"ligand": "ligand_file"}, {"parameters": "ligand_parameters"}),
    "parameterize_complex": Binding("bio_openmm_complex_parameterize", {"complex": "structure_file", "parameters": "ligand_parameters_file"}, {"system": "simulation_system"}),
    "mdanalysis": Binding("bio_mdanalysis", {"trajectory": "trajectory_file", "system": "system_file"}, {"report": "trajectory_report"}),
    "openmm": Binding("bio_openmm_run", {"system": "system_file"}, {"trajectory": "trajectory", "checkpoint": "checkpoint"}),
    "parameterize": Binding("bio_openmm_parameterize", {"complex": "structure_file"}, {"system": "simulation_system"}),
    "prolif": Binding("prolif_analyze", {"poses": "poses_file", "receptor": "receptor_file"}, {"interactions": "interaction_report"}),
    "esmfold": Binding("esmfold_predict", {"request": "request_file"}, {"models": "prediction_bundle"}),
    "receptor_prepare": Binding("bio_receptor_prepare", {"structure": "structure_file"}, {"receptor": "prepared_receptor"}),
    "assemble_complex": Binding("bio_assemble_complex", {"receptor": "receptor_file", "poses": "poses_file"}, {"complex": "complex_structure"}),
    "assay_join": Binding("bio_assay_join", {"properties": "properties_file", "assays": "assays_file"}, {"dataset": "medchem_dataset"}),
    "ligand_select": Binding("bio_ligand_select", {"ligands": "ligands_file"}, {"ligand": "ligand"}),
    "pose_select": Binding("bio_ligand_select", {"poses": "ligands_file"}, {"ligand": "ligand"}),
    "state_select": Binding("bio_ligand_select", {"states": "ligands_file"}, {"ligand": "ligand"}),
    "box_from_ligand": Binding("bio_box_from_ligand", {"ligand": "ligand_file"}, {"site": "binding_box"}),
    "posebusters": Binding("posebusters_validate", {"poses": "poses_file", "reference": "reference_file", "receptor": "receptor_file"}, {"report": "pose_validation"}),
    "meeko_ligand": Binding("meeko_ligand_prepare", {"ligand": "input_file"}, {"ligand": "ligand_pdbqt"}),
    "meeko_receptor": Binding("meeko_receptor_prepare", {"receptor": "input_file"}, {"receptor": "receptor_pdbqt"}),
    "meeko_export": Binding("meeko_pose_export", {"poses": "input_file"}, {"poses": "pose_set"}),
    "plip": Binding("plip_analyze", {"complex": "input_file"}, {"interactions": "interaction_report"}),
    "vina": Binding("vina_dock", {"receptor": "receptor_file", "ligand": "ligand_file", "site": "site_file"}, {"poses": "poses_pdbqt"}),
    "gnina": Binding("gnina_dock", {"receptor": "receptor_file", "ligand": "ligand_file", "site": "site_file"}, {"poses": "pose_set"}),
    "smina": Binding("smina_dock", {"receptor": "receptor_file", "ligand": "ligand_file", "site": "site_file"}, {"poses": "pose_set"}),
    "rdkit_prepare": Binding("ligand_prepare", {"ligand": "ligand_file"}, {"ligand": "ligand"}),
    "ligand_states": Binding("bio_ligand_states", {"ligand": "ligand_file"}, {"states": "ligand_state_set"}),
    "medchem_properties": Binding("bio_ligand_properties", {"ligands": "ligands_file"}, {"properties": "property_table"}),
    "select_model": Binding("bio_model_select", {"models": "bundle_file"}, {"structure": "macromolecule"}),
    "structure_select": Binding("bio_structure_select", {"structure": "structure_file"}, {"structure": "macromolecule"}),
    "ligand_depict": Binding("bio_ligand_depict", {"ligands": "ligands_file"}, {"image": "svg"}),
    "visualize": Binding("bio_visualize", {"structure": "structure_file"}, {"view": "html_view", "manifest": "view_manifest"}),
    "visualize_poses": Binding("bio_visualize", {"structure": "structure_file", "poses": "poses_file", "site": "site_file"}, {"view": "html_view", "manifest": "view_manifest"}),
}


def _implementation_hash(adapter) -> str:
    from delfin import qm_runtime
    from delfin.tools import _base, _runner, _products
    paths = sorted(Path(__file__).parent.glob("*.py"))
    paths.append(Path(inspect.getfile(type(adapter))))
    paths.extend(Path(inspect.getfile(module)) for module in (qm_runtime, _base, _runner, _products))
    return hashlib.sha256("".join(file_hash(p) for p in paths).encode()).hexdigest()


def preflight_plan(plan: dict, *, provision: bool = False) -> dict:
    from delfin.tools.bio.workflows import validate_plan, _INDEX
    from delfin.tools._registry import get
    structural = validate_plan(plan)
    if not structural["valid"]:
        return {**structural, "inputs": {}}
    errors, inputs, steps, pending = [], {}, [], []
    for name, item in plan["inputs"].items():
        try:
            if not isinstance(item.get("path"), str) or not item["path"]:
                raise ValueError("input path is required")
            inputs[f"inputs.{name}"] = validate_artifact(Path(item["path"]), item["type"])
        except (ValueError, OSError, KeyError, ImportError) as exc:
            errors.append(f"inputs.{name}: {exc}")
    for node in plan["nodes"]:
        binding = BINDINGS.get(node["backend"])
        if binding is None:
            errors.append(f"{node['id']}: execution binding not implemented for {node['backend']}")
            continue
        adapter = get(binding.adapter)
        if adapter is None:
            errors.append(f"{node['id']}: registered adapter {binding.adapter} unavailable")
            continue
        spec = _INDEX[node["backend"]]
        if set(binding.inputs) != set(dict(spec.inputs)) or set(binding.outputs) != set(dict(spec.outputs)):
            errors.append(f"{node['id']}: binding does not match backend ports")
            continue
        contract = adapter.contract()
        if not set(binding.outputs.values()) <= set(contract.produces):
            errors.append(f"{node['id']}: adapter does not declare all bound outputs")
        settings = node.get("settings", {})
        declared = {p.name for p in contract.params}
        unknown = set(settings) - declared
        if unknown or set(settings) & set(binding.inputs.values()):
            errors.append(f"{node['id']}: unknown or input-overriding settings {sorted(unknown or set(settings) & set(binding.inputs.values()))}")
        kwargs = dict(settings)
        for port, parameter in binding.inputs.items():
            ref = node["inputs"][port]
            kwargs[parameter] = inputs[ref]["path"] if ref in inputs else "/planned-output/" + ref
        dependencies, settings_inputs = {}, {}
        try:
            adapter.validate_params(**kwargs)
            for parameter in contract.params:
                if parameter.type == "path" and parameter.name in settings and parameter.name != "python":
                    settings_inputs[parameter.name] = validate_artifact(Path(settings[parameter.name]), "file")
            for requirement in contract.requires_python:
                if importlib.util.find_spec(requirement) is None:
                    errors.append(f"{node['id']}: Python dependency unavailable: {requirement}")
                else:
                    distributions = importlib.metadata.packages_distributions().get(requirement.split(".")[0], [])
                    dependencies[requirement] = {name: importlib.metadata.version(name) for name in distributions}
            pending.append((adapter, kwargs, node['id'], contract.requires_binaries))
        except (ValueError, ImportError, ModuleNotFoundError, OSError, subprocess.TimeoutExpired) as exc:
            runtime = None
            errors.append(f"{node['id']}: {exc}")
        steps.append({"id": node["id"], "adapter": binding.adapter,
                      "runtime": None,
                      "python_dependencies": dependencies, "settings_inputs": settings_inputs,
                      "implementation_sha256": _implementation_hash(adapter)})
    # Do not install anything until every node's inputs/settings/binding passed.
    if not errors:
        from delfin import qm_runtime
        by_id = {step['id']: step for step in steps}
        for adapter, kwargs, identity, binaries in pending:
            try:
                runtime = adapter.runtime_preflight(_allow_install=provision, **kwargs) if hasattr(adapter, 'runtime_preflight') else None
                by_id[identity]['runtime'] = runtime
                for binary in binaries:
                    if not (runtime and runtime.get('executable')) and qm_runtime.resolve_tool(binary) is None:
                        errors.append(f'{identity}: executable unavailable: {binary}')
            except (ValueError, ImportError, OSError, subprocess.TimeoutExpired) as exc:
                errors.append(f'{identity}: {exc}')
    return {"valid": True, "execution_ready": not errors, "execution_performed": False,
            "errors": [], "blockers": errors, "inputs": inputs, "steps": steps,
            "note": "Local preflight only; scientific applicability must be established separately"}


def execute_plan(plan: dict, work_dir: Path, *, cores: int = 1, resume: bool = False, cancelled=None) -> dict:
    from delfin.tools import run_step
    from delfin.tools.bio.workflows import _INDEX
    if type(cores) is not int or cores < 1:
        raise ValueError("cores must be a positive integer")
    if cancelled is not None and cancelled():
        return {'status': 'cancelled', 'execution_performed': False, 'steps': []}
    root = Path(work_dir).resolve()
    if root.exists() and (not resume or not (root / 'run.json').is_file()):
        raise ValueError('existing run directory requires resume=True and its run.json')
    preflight = preflight_plan(plan, provision=True)
    if not preflight["execution_ready"]:
        return {"status": "blocked", "execution_performed": False, "preflight": preflight}
    signature = hashlib.sha256(json.dumps({"plan": plan, "inputs": preflight["inputs"],
                                         "implementations": preflight["steps"], "cores": cores}, sort_keys=True, allow_nan=False).encode()).hexdigest()
    root = Path(work_dir).resolve()
    manifest_path = root / "run.json"
    if root.exists():
        if not resume or not manifest_path.is_file():
            raise ValueError("existing run directory requires resume=True and its run.json")
        manifest = read_json(manifest_path)
        if manifest.get("signature") != signature:
            raise ValueError("resume refused: plan, inputs, cores or implementation changed")
    else:
        root.mkdir(parents=True)
        manifest = {"schema_version": 1, "signature": signature, "plan": plan,
                    "inputs": preflight["inputs"], "steps": [], "status": "running", "execution_performed": False}
    available = {key: Path(value["path"]) for key, value in preflight["inputs"].items()}
    completed = {row["id"]: row for row in manifest["steps"] if row["status"] == "success"}
    for node in plan["nodes"]:
        name = node["id"]
        if cancelled is not None and cancelled():
            manifest["status"] = "cancelled"
            write_json(manifest_path, manifest)
            return manifest
        binding = BINDINGS[node["backend"]]
        if name in completed:
            row = completed[name]
            for relative, checksum in row["files"].items():
                if file_hash(child_path(root, relative)) != checksum:
                    raise ValueError(f"resume refused: output changed: {relative}")
            for port, artifact in row["outputs"].items():
                available[f"{name}.{port}"] = child_path(root, artifact["path"])
            continue
        attempt = 1 + sum(row["id"] == name for row in manifest["steps"])
        directory = root / f"{name}_{attempt}"
        kwargs = dict(node.get("settings", {}))
        for port, parameter in binding.inputs.items():
            kwargs[parameter] = available[node["inputs"][port]]
        row = {"id": name, "adapter": binding.adapter, "attempt": attempt, "status": "running"}
        manifest["steps"].append(row)
        manifest["status"] = "running"
        write_json(manifest_path, manifest)
        try:
            result = run_step(binding.adapter, work_dir=directory, cores=cores, **kwargs)
            manifest["execution_performed"] = True
            if not result.ok:
                raise ValueError(f"{result.error_kind.value}: {result.error}")
            outputs = {}
            for port, kind in _INDEX[node["backend"]].outputs:
                actual = Path(result.artifacts[binding.outputs[port]]).resolve()
                if not actual.is_relative_to(directory.resolve()):
                    raise ValueError("producer output is outside its step directory")
                artifact = validate_artifact(actual, kind)
                artifact["path"] = str(actual.relative_to(root))
                outputs[port] = artifact
                available[f"{name}.{port}"] = actual
            files = {str(path.relative_to(root)): file_hash(path)
                     for path in directory.rglob("*") if path.is_file()}
            row.update(status="success", outputs=outputs, files=files, data=result.data)
        except Exception as exc:
            row.update(status="failed", error=str(exc))
            manifest["status"] = "failed"
            write_json(manifest_path, manifest)
            return manifest
        write_json(manifest_path, manifest)
    manifest["status"] = "success"
    write_json(manifest_path, manifest)
    return manifest
