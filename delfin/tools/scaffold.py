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
from typing import Dict, List

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
    return name if name.endswith(f"_{kind}") else f"{name}_{kind}"


def _kind_name(name: str) -> str:
    for kind in _KINDS:
        if name.endswith(f"_{kind}"):
            return kind
    return "program"


def _skeleton(kind: str, class_name: str, adapter_name: str, description: str) -> str:
    params_doc = {
        "program": 'ParamSpec("geometry", "path", description="Input geometry file"),',
        "parser": 'ParamSpec("output", "path", required=True,\n                    description="Output file to parse"),',
        "analysis": 'ParamSpec("results", "path", required=True,\n                     description="Results file to analyse"),',
    }[kind]
    return f'''"""{{adapter_name}} — {{description}}.

Scaffolded by ``delfin-agent new-adapter``; fill in the TODO in ``execute``.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec
from delfin.tools._types import StepResult, StepStatus


class {{class_name}}(StepAdapter):
    name = "{{adapter_name}}"
    description = "{{description}}"
    produces_geometry = False
    category = "{{kind}}"
    params = (
        {{params_doc}}
    )

    def validate_params(self, **kwargs: Any) -> None:
        # required params are enforced by the contract; add custom checks here
        pass

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()
        raise NotImplementedError(
            "TODO: implement {{adapter_name}}.execute()"
        )
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
        )


register({{class_name}}())
'''.replace("{adapter_name}", adapter_name).replace("{description}", description).replace("{class_name}", class_name).replace("{kind}", kind).replace("{params_doc}", params_doc)


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

    plan.files[base / "delfin" / "tools" / "adapters" / f"{name}.py"] = _skeleton(
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
    plan.files[base / "tests" / f"test_adapter_{name}.py"] = test_body

    docs_body = f"""## {adapter_name}

{description}.

- Kind: `{kind}`
- Skeleton: `delfin/tools/adapters/{name}.py`
- Contract test: `tests/test_adapter_{name}.py`
- Registry: registered on import (no hand edit needed)

TODO: describe parameters, inputs and outputs.
"""
    plan.files[base / "docs" / "extending_delfin.md"] = docs_body

    return plan


def write_scaffold(
    name: str,
    kind: str = "program",
    *,
    directory=None,
    registry_dir=None,
) -> ScaffoldPlan:
    """``build_scaffold`` plus actually writing the files to disk."""
    plan = build_scaffold(name, kind, directory=directory, registry_dir=registry_dir)
    if not plan.ok:
        return plan
    for path, content in plan.files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            plan.errors.append(f"{path} already exists — refusing to overwrite")
            continue
        path.write_text(content, encoding="utf-8")
    return plan
