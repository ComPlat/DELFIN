"""User-defined workflows: save, load, version, list, export.

The agent-side store behind "build my own workflow": a checked PipelineSpec
(``delfin.tools._serialize`` shape, ``delfin_pipeline_spec: "1"``) is saved as
a versioned JSON file in the user's DELFIN data directory, can be listed and
loaded back (round-trip identical), re-validated without running anything,
and exported as a contribution template the user can PR into
``delfin/workflows/contrib/``.

Path resolution: the store root is resolved PER CALL — never an import-time
constant — so a redirected home (``state_paths.USER_STATE_RESOLVERS``) moves
it, and no module constant lands in the real-home guard tables
(``tests/test_a_constant_into_the_real_home_is_accounted_for.py``).
"""

from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

_SPEC_MARKER = "delfin_pipeline_spec"

_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_-]{0,63}$")


def default_root() -> Path:
    """The user's workflow store root, resolved per call.

    Honours ``DELFIN_USER_WORKFLOWS_DIR`` (a documented override, used by
    tests and by users who keep DELFIN state elsewhere); otherwise
    ``~/.delfin/workflows``.
    """
    override = os.environ.get("DELFIN_USER_WORKFLOWS_DIR")
    if override:
        return Path(override)
    return Path.home() / ".delfin" / "workflows"


def _check_name(name: str) -> str:
    if not isinstance(name, str) or not _NAME_RE.match(name):
        raise ValueError(
            "workflow name must be [A-Za-z][A-Za-z0-9_-]{0,63} "
            f"(a file name is derived from it), got: {name!r}")
    return name


def _check_spec(spec: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(spec, dict) or "name" not in spec or "steps" not in spec:
        raise ValueError(
            "spec must be a mapping with 'name' and 'steps' "
            "(the delfin_pipeline_spec shape)")
    return spec


def _dir(root: Optional[Path], name: str) -> Path:
    base = Path(root) if root is not None else default_root()
    return base / _check_name(name)


def _write_record(dir_path: Path, spec: Dict[str, Any], version: int) -> Path:
    record = {
        "name": spec.get("name"),
        "version": version,
        "saved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "spec": spec,
    }
    path = dir_path / f"v{version:03d}.json"
    path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    return path


def save(spec: Dict[str, Any], *, name: Optional[str] = None,
         root: Optional[Path] = None) -> Path:
    """Save a checked spec as the next version of *name*.

    Returns the file path written. The spec itself is stored unchanged —
    validation is the caller's job (:func:`check_spec` / the workflow
    builder); saving an unchecked spec is allowed but marked in
    :func:`list_workflows`.
    """
    _check_spec(spec)
    name = _check_name(name or spec.get("name", ""))
    dir_path = _dir(root, name)
    dir_path.mkdir(parents=True, exist_ok=True)
    version = _next_version(dir_path)
    return _write_record(dir_path, spec, version)


def _versions(dir_path: Path) -> List[int]:
    if not dir_path.is_dir():
        return []
    out = []
    for p in dir_path.glob("v*.json"):
        try:
            out.append(int(p.stem[1:]))
        except ValueError:
            continue
    return sorted(out)


def _next_version(dir_path: Path) -> int:
    known = _versions(dir_path)
    return (known[-1] + 1) if known else 1


def load_path(name: str, *, version: Optional[int] = None,
              root: Optional[Path] = None) -> Path:
    """File path of a saved workflow version (newest if *version* is None)."""
    dir_path = _dir(root, name)
    if version is None:
        versions = _versions(dir_path)
        if not versions:
            raise KeyError(name)
        version = versions[-1]
    path = dir_path / f"v{version:03d}.json"
    if not path.is_file():
        raise KeyError(f"{name} version {version}")
    return path


def load(name: str, *, version: Optional[int] = None,
         root: Optional[Path] = None) -> Dict[str, Any]:
    """Load a saved workflow's spec — round-trip identical to what was saved."""
    record = json.loads(load_path(name, version=version, root=root)
                        .read_text(encoding="utf-8"))
    return record["spec"]


def list_versions(name: str, *, root: Optional[Path] = None) -> List[Dict[str, Any]]:
    """Saved versions of one workflow, oldest first."""
    dir_path = _dir(root, name)
    out = []
    for v in _versions(dir_path):
        rec = json.loads((dir_path / f"v{v:03d}.json")
                         .read_text(encoding="utf-8"))
        out.append({"version": v, "saved_at": rec.get("saved_at", ""),
                    "description": rec["spec"].get("description", "")})
    return out


def list_workflows(*, root: Optional[Path] = None) -> List[Dict[str, Any]]:
    """All saved workflows: name, newest version, description, checked state."""
    base = Path(root) if root is not None else default_root()
    out = []
    if not base.is_dir():
        return out
    for dir_path in sorted(base.iterdir()):
        versions = _versions(dir_path)
        if not versions:
            continue
        spec = json.loads((dir_path / f"v{versions[-1]:03d}.json")
                          .read_text(encoding="utf-8"))["spec"]
        out.append({
            "name": dir_path.name,
            "version": versions[-1],
            "n_versions": len(versions),
            "description": spec.get("description", ""),
            "steps": [s.get("step", "?") for s in spec.get("steps", [])],
        })
    return out


# ---------------------------------------------------------------------------
# Extending: start from a named template or a saved workflow, mutate, return
# a NEW spec. The original is never touched.
# ---------------------------------------------------------------------------


def _source_dict(name: str, *, root: Optional[Path] = None) -> Dict[str, Any]:
    """Spec dict of a saved workflow OR a named built-in template."""
    from delfin.tools import templates

    template = getattr(templates, name, None)
    if template is not None and hasattr(template, "to_dict"):
        return template.to_dict()
    return load(name, root=root)


def _apply_edits(spec: Dict[str, Any], *, add: Optional[List[Dict[str, Any]]] = None,
                 remove: Optional[List[str]] = None,
                 replace: Optional[Dict[str, Dict[str, Any]]] = None,
                 overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    steps = list(spec.get("steps", []))
    if remove:
        labels = set(remove)
        kept = [s for s in steps if s.get("label", s.get("step")) not in labels]
        if len(kept) == len(steps):
            raise KeyError(f"no step matched remove={sorted(labels)}")
        steps = kept
    if replace:
        replaced = False
        for i, s in enumerate(steps):
            key = s.get("label", s.get("step"))
            if key in replace:
                steps[i] = dict(s, kwargs=dict(s.get("kwargs", {}),
                                               **replace[key]["kwargs"]))
                replaced = True
        if not replaced:
            raise KeyError(f"no step matched replace={sorted(replace)}")
    if overrides:
        spec = dict(spec)
        spec["defaults"] = dict(spec.get("defaults", {}), **overrides)
    spec = dict(spec)
    spec["steps"] = steps + [dict(a) for a in (add or [])]
    return spec


def extend_template(name: str, *, add: Optional[List[Dict[str, Any]]] = None,
                    remove: Optional[List[str]] = None,
                    replace: Optional[Dict[str, Dict[str, Any]]] = None,
                    overrides: Optional[Dict[str, Any]] = None,
                    new_name: Optional[str] = None,
                    root: Optional[Path] = None) -> Dict[str, Any]:
    """Build a NEW spec from a saved workflow or named template.

    ``add`` steps go last; ``remove``/``replace`` match a step by its label
    (falling back to the step name); ``overrides`` extend the spec's
    ``defaults``. The source template or saved file is read-only here:
    a fresh dict is returned, the original stays untouched.
    """
    spec = dict(_source_dict(name, root=root))
    new_name = _check_name(new_name or spec.get("name", "") or "workflow")
    spec["name"] = new_name
    spec.pop("delfin_pipeline_spec", None)
    return _apply_edits(spec, add=add, remove=remove, replace=replace,
                        overrides=overrides)


# ---------------------------------------------------------------------------
# Checking (read-only) and re-running: validation through the platform.
# ---------------------------------------------------------------------------


def check_spec(spec: Dict[str, Any], *, inputs: Optional[Dict[str, Any]] = None,
               geometry: bool = False) -> Dict[str, Any]:
    """The mandatory pre-save check chain, in order.

    resolve_spec -> validate_spec -> scientific_lint -> dry run (validate with
    geometry flag). Returns ``{"ok", "resolved", "report", "lint", "dry_run"}``;
    ``ok`` is False when ANY check fails. Nothing is executed — a dry run is
    ``Pipeline.validate``, never ``Pipeline.run``.
    """
    from delfin.tools import platform

    resolved = platform.resolve_spec(spec)
    report = platform.validate_spec(spec, inputs=inputs, geometry=geometry)
    lint = platform.scientific_lint(spec)
    errors = list(report.errors) if report is not None else []
    if not lint_ok(lint):
        errors = errors + [f"[lint] {f['step']}: {f['message']}" for f in lint
                           if f.get("level") == "error"]
    return {"ok": report is not None and report.ok and lint_ok(lint),
            "resolved": resolved, "report": report, "lint": lint,
            "errors": errors, "dry_run": report}


def lint_ok(lint: List[Dict[str, Any]]) -> bool:
    """A lint result passes when no finding has level ``error``."""
    return not any(f.get("level") == "error" for f in lint)


def validate_saved(name: str, *, version: Optional[int] = None,
                   root: Optional[Path] = None, geometry: bool = False):
    """Re-validate a saved workflow through the platform (no execution)."""
    return check_spec(load(name, version=version, root=root),
                      geometry=geometry)["report"]


def load_spec_from_export(data: Dict[str, Any]) -> Dict[str, Any]:
    """The spec embedded in an exported contribution template.

    Runs the pipeline round-trip (``_serialize.from_dict``) on the embedded
    spec so an export that cannot load back fails loudly here.
    """
    from delfin.tools._serialize import from_dict

    if not isinstance(data, dict) or "spec" not in data:
        raise ValueError("contribution template must be an envelope with 'spec'")
    from_dict(data["spec"])
    return data["spec"]


def export_contrib(name: str, *, out_dir: Union[str, Path],
                   root: Optional[Path] = None,
                   author: str = "") -> Path:
    """Export a saved workflow as a contribution template file.

    The file is the spec wrapped in a JSON envelope with provenance; it can
    be PR-ed into ``delfin/workflows/contrib/``. Re-loading it through
    :func:`load_spec_from_export` must reproduce the saved spec.
    """
    spec = load(name, root=root)
    envelope = {
        "delfin_contrib_workflow": "1",
        "name": spec.get("name", name),
        "author": author,
        "exported_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "spec": spec,
    }
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{_check_name(name)}.json"
    path.write_text(json.dumps(envelope, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    return path


__all__ = [
    "default_root", "save", "load", "load_path", "list_workflows",
    "list_versions", "extend_template", "check_spec", "lint_ok",
    "validate_saved", "export_contrib", "load_spec_from_export",
]
