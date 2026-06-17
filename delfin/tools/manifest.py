"""Machine-readable manifest of the platform — the binding contract.

Emits a single JSON document describing everything an external consumer (the
agent system, a UI, another tool) needs to bind to DELFIN's tools without any
Python coupling:

* **capabilities** — every building block's full contract (params, ports,
  requirements)
* **applications** — every registered workflow's input/output contract
* **keys** — the central well-known-key vocabulary (``functional`` / ``solvent``
  / … with their allowed values), so producers know which values are valid
* **schemas** — JSON Schemas for a PipelineSpec and an Application, so producers
  can author and validate workflow definitions as data

Pure data; no third-party dependency, fully testable.
"""

from __future__ import annotations

from typing import Any, Dict

from delfin.tools._application import (
    _outputspec_to_dict,
    _paramspec_to_dict,
    list_applications,
)
from delfin.tools._keys import list_keys
from delfin.tools._registry import list_steps

MANIFEST_VERSION = "1"


def _datakey_to_dict(d) -> Dict[str, Any]:
    out: Dict[str, Any] = {"name": d.name, "type": d.type}
    if d.unit:
        out["unit"] = d.unit
    if d.description:
        out["description"] = d.description
    return out


def contract_to_dict(c) -> Dict[str, Any]:
    """Serialize a :class:`StepContract` to a plain dict for the manifest."""
    return {
        "name": c.name,
        "description": c.description,
        "category": c.category,
        "produces_geometry": c.produces_geometry,
        "params": [_paramspec_to_dict(p) for p in c.params],
        "consumes": sorted(c.consumes),
        "produces": sorted(c.produces),
        "data_keys": [_datakey_to_dict(d) for d in c.data_keys],
        "requires_binaries": sorted(c.requires_binaries),
        "requires_python": sorted(c.requires_python),
    }


def _keyspec_to_dict(ks) -> Dict[str, Any]:
    out: Dict[str, Any] = {"name": ks.name, "type": ks.type}
    if ks.description:
        out["description"] = ks.description
    if ks.default is not None:
        out["default"] = ks.default
    if ks.enum is not None:
        out["enum"] = list(ks.enum)
    if ks.unit:
        out["unit"] = ks.unit
    if ks.enum_source:
        out["enum_source"] = ks.enum_source
    return out


def _application_summary(app) -> Dict[str, Any]:
    return {
        "name": app.name,
        "description": app.description,
        "category": app.category,
        "version": app.version,
        "inputs": [_paramspec_to_dict(p) for p in app.inputs],
        "outputs": [_outputspec_to_dict(o) for o in app.outputs],
    }


# --- JSON Schemas (draft 2020-12) ----------------------------------------


def pipeline_spec_schema() -> Dict[str, Any]:
    """JSON Schema for a serialized pipeline (the ``from_dict`` input)."""
    step = {
        "type": "object",
        "properties": {
            "step": {"type": "string", "description": "Registered capability name"},
            "type": {
                "type": "string",
                "enum": ["normal", "step", "if", "loop", "retry", "checkpoint",
                         "compute", "reactive", "transform", "fan_out", "map",
                         "sub_pipeline"],
                "default": "normal",
            },
            "label": {"type": "string"},
            "kwargs": {"type": "object", "description": "Adapter parameters"},
            "condition": {"type": "string", "description": "DSL for type=if"},
            "until": {"type": "string", "description": "DSL for type=loop"},
            "max_iter": {"type": "integer"},
            "max_attempts": {"type": "integer"},
            "ref": {"type": "string", "description": "Named-callable ref"},
        },
        "additionalProperties": True,
    }
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "DelfinPipelineSpec",
        "type": "object",
        "required": ["name", "steps"],
        "properties": {
            "name": {"type": "string"},
            "kind": {"type": "string", "enum": ["pipeline", "template"]},
            "template": {"type": "boolean"},
            "defaults": {"type": "object"},
            "steps": {"type": "array", "items": step},
            "branches": {
                "type": "object",
                "additionalProperties": {"type": "array", "items": step},
            },
        },
    }


def _param_schema() -> Dict[str, Any]:
    return {
        "type": "object",
        "required": ["name"],
        "properties": {
            "name": {"type": "string"},
            "type": {"type": "string",
                     "enum": ["str", "int", "float", "bool", "list", "dict", "path"]},
            "required": {"type": "boolean"},
            "default": {},
            "enum": {"type": "array"},
            "description": {"type": "string"},
            "unit": {"type": "string"},
        },
    }


def application_schema() -> Dict[str, Any]:
    """JSON Schema for a serialized application (``Application.to_dict``)."""
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "DelfinApplication",
        "type": "object",
        "required": ["name", "spec"],
        "properties": {
            "name": {"type": "string"},
            "description": {"type": "string"},
            "category": {"type": "string"},
            "delfin_application": {"type": "string", "description": "schema version"},
            "inputs": {"type": "array", "items": _param_schema()},
            "outputs": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": ["name", "step", "key"],
                    "properties": {
                        "name": {"type": "string"},
                        "step": {"type": "string"},
                        "key": {"type": "string"},
                        "branch": {"type": "string"},
                        "type": {"type": "string"},
                        "unit": {"type": "string"},
                        "description": {"type": "string"},
                    },
                },
            },
            "spec": {"$ref": "#/$defs/pipeline_spec"},
        },
        "$defs": {"pipeline_spec": pipeline_spec_schema()},
    }


# --- the build guide (read this first) -----------------------------------


def build_guide() -> Dict[str, Any]:
    """A short, authoritative recipe for assembling/validating/running/extending
    a pipeline — the rails a model (even a weak one) should follow."""
    return {
        "overview": (
            "Build a workflow as a pipeline spec (JSON matching schemas."
            "pipeline_spec), validate it, then run it. Specify only the truly "
            "required params — the framework fills declared defaults and wires "
            "artifacts between steps. If no building block does what you need, "
            "build one and register it. Study failures via diagnostics and iterate."
        ),
        "rules": [
            "Specify only required params; declared defaults (method, basis, "
            "maxcore, …) and auto-wiring (e.g. gbw→moread, hessian→hess_file) "
            "are filled for you — call resolve_spec to see the full resolved spec.",
            "Prefer values allowed by a param's enum / the key vocabulary (see "
            "keys / describe_capability); validate flags others with the allowed set.",
            "Validate before running; each diagnostic names the concrete fix.",
            "License-restricted engines (ORCA, Turbomole) are never auto-installed.",
        ],
        "workflow": [
            {"step": "discover", "do": "Find building blocks and allowed values",
             "tools": ["get_manifest", "list_capabilities", "describe_capability",
                       "catalog", "list_keys", "describe_key", "compatible_successors"]},
            {"step": "assemble", "do": "Write a pipeline spec (schemas.pipeline_spec)",
             "tools": ["pipeline_spec schema", "resolve_spec"]},
            {"step": "validate", "do": "Statically check it; apply the suggested fixes",
             "tools": ["validate_spec"]},
            {"step": "build-missing", "do": "No block fits? scaffold + integrate one",
             "tools": ["new_capability_template", "register_module"]},
            {"step": "persist", "do": "Save it so it appears in the Pipelines tab",
             "tools": ["save_application"]},
            {"step": "run", "do": "Execute locally or on SLURM; results land in ~/calc",
             "tools": ["run_application", "submit_application", "run_status"]},
            {"step": "diagnose", "do": "Study what failed and iterate",
             "tools": ["run_diagnostics", "run_metrics"]},
        ],
    }


# --- the manifest ---------------------------------------------------------


def build_manifest() -> Dict[str, Any]:
    """Assemble the full platform manifest as a plain dict."""
    capabilities = [contract_to_dict(a.contract())
                    for _, a in sorted(list_steps().items())]
    applications = [_application_summary(app)
                    for _, app in sorted(list_applications().items())]
    keys = [_keyspec_to_dict(ks) for _, ks in sorted(list_keys().items())]
    return {
        "delfin_tools_manifest": MANIFEST_VERSION,
        "guide": build_guide(),
        "capabilities": capabilities,
        "applications": applications,
        "keys": keys,
        "schemas": {
            "pipeline_spec": pipeline_spec_schema(),
            "application": application_schema(),
        },
    }


def manifest_json(*, indent: int = 2) -> str:
    """The manifest as a JSON string."""
    import json
    return json.dumps(build_manifest(), indent=indent, default=str)


__all__ = [
    "MANIFEST_VERSION",
    "contract_to_dict",
    "pipeline_spec_schema",
    "application_schema",
    "build_guide",
    "build_manifest",
    "manifest_json",
]
