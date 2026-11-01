"""Scaffold a new DELFIN step adapter: skeleton, contract test, registry entry.

``delfin-agent new-adapter <name> --kind <program|parser|analysis>`` calls
into this module.  Everything it writes is English source that passes the
contract tests in ``tests/test_tool_contracts.py`` unchanged: a
:class:`~delfin.tools._base.StepAdapter` subclass with a declared category,
declared params and ports, ``validate_params`` and an ``execute`` that
returns a :class:`~delfin.tools._types.StepResult`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

# Same rule the register_module writer enforces (delfin/tools/platform.py:122):
# the name becomes a file name and is imported immediately.
_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

_KINDS = {
    "program": (
        "program",
        "Wrap an external program as a pipeline step",
        "executable",
        "Run {name} on the geometry",
    ),
    "parser": (
        "parser",
        "Parse a tool's output into structured data",
        "",
        "Parse {name} output files",
    ),
    "analysis": (
        "analysis",
        "Analyse results downstream of a calculation",
        "",
        "Analyse {name} results",
    ),
}


@dataclass
class ScaffoldPlan:
    """What ``build_scaffold`` intends to write, and where."""

    module_name: str
    kind: str
    class_name: str
    adapter_name: str
    files: Dict[Path, str] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def _class_name(name: str) -> str:
    return "".join(part.capitalize() for part in name.split("_")) + "Adapter"


def _adapter_name(name: str, kind: str) -> str:
    return name if name.endswith("_" + kind) else f"{name}_{kind}"


def _kind_name(name: str) -> str:
    for kind in _KINDS:
        if name.endswith("_" + kind):
            return kind
    return "program"


def _skeleton(kind: str, class_name: str, adapter_name: str, description: str) -> str:
    params_doc = {
        "program": 'ParamSpec("geometry", "path", description="Input geometry file"),',
        "parser": 'ParamSpec("output", "path", required=True,\n                    description="Output file to parse"),',
        "analysis": 'ParamSpec("results", "path", required=True,\n                     description="Results file to analyse"),',
    }[kind]
    body = (
        "from __future__ import annotations\n"
        "\n"
        "import time\n"
        "from pathlib import Path\n"
        "from typing import Any, Optional\n"
        "\n"
        "from delfin.tools._base import StepAdapter\n"
        "from delfin.tools._registry import register\n"
        "from delfin.tools._spec import ParamSpec\n"
        "from delfin.tools._types import StepResult, StepStatus\n"
        "\n"
        "\n"
        "class CLASS_NAME(StepAdapter):\n"
        '    name = "ADAPTER_NAME"\n'
        '    description = "DESCRIPTION"\n'
        "    produces_geometry = False\n"
        '    category = "KIND"\n'
        "    params = (\n"
        "PARAMS_DOC\n"
        "    )\n"
        "\n"
        "    def validate_params(self, **kwargs: Any) -> None:\n"
        "        for spec in self.params:\n"
        "            if spec.required and spec.name not in kwargs:\n"
        "                raise ValueError(\n"
        "                    \"'{}' parameter is required\".format(spec.name))\n"
        "\n"
        "    def execute(\n"
        "        self,\n"
        "        work_dir: Path,\n"
        "        *,\n"
        "        geometry: Optional[Path] = None,\n"
        "        cores: int = 1,\n"
        "        **kwargs: Any,\n"
        "    ) -> StepResult:\n"
        "        start = time.monotonic()\n"
        "        raise NotImplementedError(\n"
        "            \"TODO: implement ADAPTER_NAME.execute()\"\n"
        "        )\n"
        "\n"
        "\n"
        "register(CLASS_NAME())\n"
    )
    body = (
        body.replace("CLASS_NAME", class_name)
        .replace("ADAPTER_NAME", adapter_name)
        .replace("DESCRIPTION", description)
        .replace("KIND", kind)
        .replace("PARAMS_DOC", params_doc)
    )
    header = (
        '"""' + adapter_name + " — " + description + ".\n"
        "\n"
        "Scaffolded by ``delfin-agent new-adapter``; fill in the TODO in\n"
        "``execute``.\n"
        '"""\n'
    )
    return header + body


def build_scaffold(
    name: str,
    kind: str = "program",
    *,
    directory=None,
    registry_dir=None,
) -> ScaffoldPlan:
    """Plan (and return) the files a new adapter needs; the caller writes them.

    Writes: (1) the adapter skeleton, (2) a contract test in the
    ``tests/test_tool_contracts.py`` style, (3) the docs stub entry.  The
    registry entry is not a file — the skeleton ends with ``register(...)``,
    so the entry exists the moment the module is imported (the wave-20 OR1
    defect was adapters whose module was never imported at all).
    """
    plan = ScaffoldPlan(module_name=name, kind=kind, class_name="", adapter_name="")

    if not isinstance(name, str) or not _NAME_RE.match(name):
        plan.errors.append(
            f"invalid adapter name {name!r}: it becomes a file name and is "
            "imported immediately, so it must be a plain Python identifier"
        )
        return plan

    if kind not in _KINDS:
        plan.errors.append(
            f"unknown kind {kind!r}: choose one of {sorted(_KINDS)}"
        )
        return plan

    category, description, _requires, _docs = _KINDS[kind]
    class_name = _class_name(name)
    adapter_name = _adapter_name(name, kind)
    plan.class_name = class_name
    plan.adapter_name = adapter_name

    base = Path(directory) if directory else None
    if base is None:
        plan.errors.append(
            "a target directory is required — pass directory= (the repo root "
            "or a user adapter directory)"
        )
        return plan

    # Destination layout: inside a REPO, the skeleton lives in
    # delfin/tools/adapters/ where production discovery walks it; inside a
    # USER adapter directory (flat *.py, exactly what _load_user_adapters
    # reads — delfin/tools/_registry.py) the skeleton must sit at the top
    # level, or the loader never sees it.  Detection: a REPO is a tree
    # containing delfin/tools/adapters/ (a checkout); an empty directory
    # defaults to the repo layout, matching the documented CLI usage
    # ("new-adapter ... --dir <repo-root>").  Everything else — notably a
    # directory named for adapters, e.g. $DELFIN_ADAPTERS_DIR — is a user
    # adapter directory and gets the flat layout.
    is_repo = (base / "delfin" / "tools" / "adapters").is_dir()
    is_repo = is_repo or (
        not any(base.iterdir()) and not base.name.startswith("adapter")
    )
    if is_repo:
        skeleton_path = base / "delfin" / "tools" / "adapters" / f"{name}.py"
        test_path = base / "tests" / f"test_adapter_{name}.py"
        docs_path = base / "docs" / "extending_delfin.md"
        skeleton_link = f"delfin/tools/adapters/{name}.py"
    else:
        skeleton_path = base / f"{name}.py"
        test_path = base / f"test_adapter_{name}.py"
        docs_path = base / "ADAPTERS.md"
        skeleton_link = f"{name}.py (user adapter directory)"

    plan.files[skeleton_path] = _skeleton(
        kind, class_name, adapter_name, description
    )

    test_body = f'''"""Contract test for the scaffolded ``{{adapter_name}}`` adapter."""

from __future__ import annotations

from delfin.tools._registry import list_steps


def test_{{name}}_is_registered():
    names = set(list_steps())
    assert "{{adapter_name}}" in names


def test_{{name}}_contract_builds():
    adapter = list_steps()["{{adapter_name}}"]
    contract = adapter.contract()
    assert contract.category == "{{kind}}"
    assert contract.params, "declared params must not be empty"
'''
    plan.files[test_path] = test_body

    docs_body = f"""## {adapter_name}

{description}.

- Kind: `{kind}`
- Skeleton: `{skeleton_link}`
- Registry: registered on import (no hand edit needed)

TODO: describe parameters, inputs and outputs.
"""
    plan.files[docs_path] = docs_body

    return plan


def write_scaffold(
    name: str,
    kind: str = "program",
    *,
    directory=None,
    registry_dir=None,
) -> ScaffoldPlan:
    """``build_scaffold`` plus actually writing the files to disk.

    A file that already exists elsewhere in the plan (e.g. the shared
    docs stub when scaffolding several adapters into one tree) is only
    *appended to*, never overwritten; a pre-existing per-adapter file is
    refused.
    """
    plan = build_scaffold(name, kind, directory=directory, registry_dir=registry_dir)
    if not plan.ok:
        return plan
    for path, content in plan.files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            text = path.read_text(encoding="utf-8")
            if name in text:
                plan.errors.append(
                    f"{path} already documents {name} — refusing to overwrite"
                )
                continue
            path.write_text(text + "\n" + content, encoding="utf-8")
            continue
        path.write_text(content, encoding="utf-8")
    return plan


# ---------------------------------------------------------------------------
# New-workflow scaffold — descriptor + steps + acceptance test
# ---------------------------------------------------------------------------

@dataclass
class WorkflowScaffoldPlan:
    """What ``build_workflow_scaffold`` writes for a new workflow, and where."""

    name: str
    kind: str
    spec: Dict[str, Any]
    module_path: Path
    test_path: Path
    errors: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def _workflow_skeleton(name: str, kind: str, spec: Dict[str, Any]) -> str:
    import json

    body = '"""{} — a DELFIN workflow (scaffolded)."""\n\n'.format(name)
    body += (
        "from delfin.tools.workflows import WorkflowDescriptor\n"
        "from delfin.tools._application import register_application\n\n"
        "# The ONE contract: a WorkflowDescriptor extends Application (no\n"
        "# parallel type).  It describes the workflow; the steps are the\n"
        "# building blocks that run in order (a pipeline) or the jobs a\n"
        "# classic CONTROL.txt workflow writes.  Fill in purpose/entry and\n"
        "# extend the steps with building blocks from delfin.tools.\n"
        "workflow = WorkflowDescriptor(\n"
        "    name={name!r},\n"
        "    kind={kind!r},\n"
        "    purpose='TODO: one sentence on what this workflow is for',\n"
        "    entry=\"run_application('{name}')\",\n"
        "    steps={steps_repr},\n"
        "    spec={spec_repr},\n"
        ")\n"
        "register_application(workflow)\n"
    ).format(name=name, kind=kind,
             steps_repr=json.dumps(spec.get("steps", []), indent=4, default=str),
             spec_repr=json.dumps(spec, indent=4, default=str))
    return body


def build_workflow_scaffold(
    name: str,
    *,
    kind: str = "pipeline",
    directory=None,
    register: bool = False,
    steps: List[Dict[str, Any]] = (),
) -> WorkflowScaffoldPlan:
    """Scaffold a NEW workflow in the style of the existing adapter scaffold.

    Plans a :class:`WorkflowScaffoldPlan` whose ``spec`` is the validated
    PipelineSpec the workflow runs, then writes (1) a ``<name>.py`` module in
    ``delfin/tools/adapters/`` (repo layout) or at the given directory that
    registers a :class:`~delfin.tools.workflows.WorkflowDescriptor`, and
    (2) an acceptance test asserting the workflow is listed and its contract
    is intact.

    The spec is checked through the SAME existing checks a real workflow
    must pass (:func:`delfin.tools.user_workflows.check_spec`), so a scaffold
    that would not run is refused at creation time rather than after a run.
    """
    from pathlib import Path

    errors: List[str] = []
    if not isinstance(name, str) or not _NAME_RE.match(name):
        errors.append(
            f"invalid workflow name {name!r}: it becomes a module and file "
            "name, so it must be a plain Python identifier"
        )
    kk = kind if kind in ("classic-control", "pipeline", "user") else None
    if kk is None:
        errors.append(
            "unknown kind: choose classic-control | pipeline | user"
        )

    spec: Dict[str, Any] = {
        "delfin_pipeline_spec": "1",
        "name": name,
        "kind": "pipeline",
        "defaults": {},
        "steps": [dict(s) for s in steps],
    }

    # The same checks a real workflow must pass.
    if not errors:
        from delfin.tools import user_workflows as uw

        try:
            report = uw.check_spec(spec)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"workflow spec rejected before writing: {exc}")
            report = None
        if report is not None and not report["ok"]:
            errors.extend(f"[spec] {e}" for e in report["errors"])

    base = Path(directory) if directory else Path.cwd()
    if not errors and not base.is_dir():
        errors.append(f"target directory does not exist: {base}")

    module_repr = base / f"{name}.py"
    test_repr = base / f"test_workflow_{name}.py"
    plan = WorkflowScaffoldPlan(
        name=name, kind=kk or kind, spec=spec,
        module_path=module_repr, test_path=test_repr, errors=errors,
    )
    if not plan.ok:
        return plan

    base.mkdir(parents=True, exist_ok=True)
    if module_repr.exists():
        plan.errors.append(f"{module_repr} already exists — refusing to overwrite")
        return plan
    module_repr.write_text(_workflow_skeleton(name, kk, spec), encoding="utf-8")

    test_body = f'''"""Acceptance test for the scaffolded ``{{name}}`` workflow."""

from __future__ import annotations

from delfin.tools.workflows import workflow_for, register_classic_workflows


def test_{{name}}_is_registered_as_a_pipeline_workflow():
    "A saved umbrella workflow still lists {name} once the scaffold is imported."
    from delfin.tools._application import list_applications
    apps = list_applications()
    assert "{name}" in apps, "scaffolded workflow is not registered"


def test_{{name}}_spec_passes_the_existing_checks():
    """The scaffolded spec must pass the checks a real workflow passes."""
    from delfin.tools import user_workflows as uw
    app = list_applications()["{name}"]
    report = uw.check_spec(app.spec)
    assert report["ok"], report["errors"]
'''
    plan.test_path.write_text(test_body, encoding="utf-8")
    return plan


def write_workflow_scaffold(
    name: str,
    *,
    kind: str = "pipeline",
    directory=None,
    register: bool = False,
    steps: List[Dict[str, Any]] = (),
) -> WorkflowScaffoldPlan:
    """Alias of :func:`build_workflow_scaffold` — the function writes the
    files itself, matching how the test drives it and the CLI usage."""
    return build_workflow_scaffold(
        name, kind=kind, directory=directory, register=register, steps=steps,
    )
