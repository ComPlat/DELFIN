"""Separate agent entry points for planning, input preparation and execution."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec, DataKeySpec
from delfin.tools._types import StepStatus, ErrorKind
from delfin.tools.bio.workflows import backend_catalog, validate_plan, workflow_template
from delfin.tools.bio.prediction_inputs import EXPORTERS, build_prediction_input


class BioWorkflowPlan(StepAdapter):
    name = "bio_workflow_plan"
    description = "Validate a typed biomolecular planning graph or generate a template; report missing adapters, never execute planned tools"
    category = "meta"
    produces_geometry = False
    produces = ("json_report",)
    required_artifacts = produces
    example_kwargs = {"template": "docking"}
    params = (ParamSpec("plan", "dict"), ParamSpec("template", "str",
              description="processing, docking, medchem, md, or predictor ID; exactly one of plan/template"))
    data_keys = (DataKeySpec("valid", "bool"), DataKeySpec("execution_ready", "bool"), DataKeySpec("blockers", "list"))

    def validate_params(self, **kwargs):
        if (kwargs.get("plan") is None) == (kwargs.get("template") is None):
            raise ValueError("exactly one of plan or template is required")
        if kwargs.get("template") is not None:
            workflow_template(kwargs["template"])

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start = time.monotonic()
        try:
            self.validate_params(**kwargs)
            plan = kwargs.get("plan") if kwargs.get("template") is None else workflow_template(kwargs["template"])
            report = validate_plan(plan)
            report["plan"] = plan
            path = Path(work_dir) / "bio_workflow_plan.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
            return self._make_result(self.name, StepStatus.SUCCESS if report["valid"] else StepStatus.FAILED,
                                     work_dir, start, data=report, artifacts={"json_report": path},
                                     error=None if report["valid"] else "; ".join(report["errors"]),
                                     error_kind=ErrorKind.NONE if report["valid"] else ErrorKind.UNSUITABLE_INPUT)
        except (ValueError, TypeError) as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc), error_kind=ErrorKind.MISSING_PARAM)


class BioBackendCatalog(StepAdapter):
    name = "bio_backend_catalog"
    description = "Describe biomolecular backends, typed ports, prerequisites and integration status; optional PATH probe never runs or installs engines"
    category = "meta"
    produces_geometry = False
    produces = ("json_report",)
    required_artifacts = produces
    example_kwargs = {"task": "prediction", "probe": False}
    params = (ParamSpec("task", "str"), ParamSpec("probe", "bool", default=False))
    data_keys = (DataKeySpec("backends", "list"),)

    def validate_params(self, **kwargs):
        if type(kwargs.get("probe", False)) is not bool:
            raise ValueError("probe must be boolean")
        backend_catalog(task=kwargs.get("task"))

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start = time.monotonic()
        try:
            self.validate_params(**kwargs)
            report = backend_catalog(task=kwargs.get("task"), probe=kwargs.get("probe", False))
            path = Path(work_dir) / "bio_backends.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(report, indent=2), encoding="utf-8")
            return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start,
                                     data=report, artifacts={"json_report": path})
        except (ValueError, TypeError) as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc), error_kind=ErrorKind.MISSING_PARAM)


class BioPredictionInput(StepAdapter):
    name = "bio_prediction_input"
    description = "Convert explicit protein/DNA/RNA/SMILES entities to a native prediction input and provenance manifest; no prediction execution"
    category = "structure"
    produces_geometry = False
    produces = ("prediction_input", "json_report")
    required_artifacts = produces
    example_kwargs = {"backend": "alphafold3", "request": {
        "schema_version": 1, "name": "example", "seeds": [1], "msa_policy": "single_sequence",
        "entities": [{"id": "A", "kind": "protein", "sequence": "MKWVTFISLL"}]}}
    requires_python = ("rdkit",)
    params = (ParamSpec("backend", "str", required=True, enum=EXPORTERS),
              ParamSpec("request", "dict", required=True,
                        description="Version 1 request with name, typed entities, seeds and explicit msa_policy; schema from bio_backends"))
    data_keys = (DataKeySpec("backend", "str"), DataKeySpec("artifact_type", "str"),
                 DataKeySpec("execution_ready", "bool"))

    def validate_params(self, **kwargs):
        if not kwargs.get("backend") or kwargs.get("request") is None:
            raise ValueError("backend and request are required")
        build_prediction_input(kwargs["backend"], kwargs["request"])

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start = time.monotonic()
        try:
            result = build_prediction_input(kwargs.get("backend"), kwargs.get("request"))
            content = result.pop("content")
            directory = Path(work_dir)
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / result["filename"]
            manifest = directory / "prediction_input_manifest.json"
            if path.exists() or manifest.exists():
                raise ValueError("prediction input outputs already exist; use a new step directory")
            result["input_sha256"] = hashlib.sha256(content.encode()).hexdigest()
            path.write_text(content, encoding="utf-8")
            manifest.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
            return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start,
                                     data=result, artifacts={"prediction_input": path, "json_report": manifest})
        except (ValueError, TypeError) as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc), error_kind=ErrorKind.MISSING_PARAM)


register(BioWorkflowPlan())
register(BioBackendCatalog())
register(BioPredictionInput())


class BioWorkflowExecution(StepAdapter):
    category = "meta"
    produces_geometry = False
    produces = ("json_report",)
    required_artifacts = produces
    params = (ParamSpec("plan", "dict", required=True), ParamSpec("resume", "bool", default=False))
    data_keys = (DataKeySpec("status", "str"), DataKeySpec("execution_performed", "bool"))
    example_kwargs = {"plan": {"schema_version": 1, "inputs": {"ligand": {"type": "ligand", "path": "ligand.sdf"}},
                               "nodes": [{"id": "prepare", "backend": "rdkit_prepare", "inputs": {"ligand": "inputs.ligand"}}]}}

    def __init__(self, preflight=False):
        self.preflight = preflight
        self.name = "bio_workflow_preflight" if preflight else "bio_workflow_execute"
        self.description = ("Check actual input files, adapter bindings and local environments; optional version probes run but no calculation"
                            if preflight else "Execute a bound biomolecular DAG through DELFIN adapters with explicit edges, failure stop and checksum-verified resume")

    def validate_params(self, **kwargs):
        if kwargs.get("plan") is None:
            raise ValueError("plan is required")
        result = validate_plan(kwargs["plan"])
        if not result["valid"]:
            raise ValueError("; ".join(result["errors"]))
        if type(kwargs.get("resume", False)) is not bool:
            raise ValueError("resume must be boolean")

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        from delfin.tools.bio.execution import execute_plan, preflight_plan
        start = time.monotonic()
        try:
            self.validate_params(**kwargs)
            if self.preflight:
                data = preflight_plan(kwargs["plan"])
                ok = data["execution_ready"]
            else:
                data = execute_plan(kwargs["plan"], Path(work_dir) / "workflow", cores=cores, resume=kwargs.get("resume", False))
                ok = data["status"] == "success"
            path = Path(work_dir) / "workflow_report.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")
            return self._make_result(self.name, StepStatus.SUCCESS if ok else StepStatus.FAILED, work_dir, start,
                                     data=data, artifacts={"json_report": path},
                                     error=None if ok else "Workflow blocked or failed; see workflow_report.json",
                                     error_kind=ErrorKind.NONE if ok else ErrorKind.TOOL_FAILED)
        except (ValueError, OSError) as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc), error_kind=ErrorKind.UNSUITABLE_INPUT)


register(BioWorkflowExecution(preflight=True))
register(BioWorkflowExecution())
