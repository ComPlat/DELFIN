"""Every packaged benchmark task is still offered the tools it is graded on.

For each task the workspace is resolved the way benchmark_runner does it
(a fixture folder by task class / mode, else the checkout, a DELFIN
tree), the fixture is copied into tmp, what the task's setup script
writes is seeded (chem/<molecule>/start.xyz), and the tool-group
decision is taken there. A required expected signal that names tools
is unsatisfiable when EVERY tool it names belongs to a hidden group.

Measured when written: 126 tasks; 28 name a group or DELFIN-only tool.
Red on the first content sniff: the science fixture's xtb outputs say
"xtb output" in lower case, so the chemistry group stayed hidden and two
science tasks graded on the ops parsers (with search_calcs, DELFIN-only,
as the only other option) could not pass.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pytest

from delfin.agent import benchmark as B
from delfin.agent.api_client import (
    _DELFIN_ONLY_TOOL_NAMES,
    _DOC_TOOLS_OPENAI,
    _MCP_READONLY_TOOL_BASES,
    _OFFICE_TOOL_NAMES,
)
from delfin.agent.benchmark_runner import workspace_for
from delfin.agent.toolsets import GROUPS, WorkspaceToolsets

_REPO = Path(__file__).resolve().parents[1]

_NATIVE = {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}
_GROUP_OF: dict[str, str | None] = {n: None for n in _NATIVE}
for _n in _MCP_READONLY_TOOL_BASES - _NATIVE:
    _GROUP_OF[_n] = "chemistry"
_GROUP_OF["show_molecule"] = "chemistry"
for _n in _OFFICE_TOOL_NAMES:
    _GROUP_OF[_n] = "office"
# Not a tool group, but equally absent outside a DELFIN tree: the doc
# and calc index tools are stripped there (api_client, stream_message).
for _n in _DELFIN_ONLY_TOOL_NAMES:
    _GROUP_OF[_n] = "delfin_only"
_NAME_RE = re.compile(r"\b(" + "|".join(
    sorted(map(re.escape, _GROUP_OF), key=len, reverse=True)) + r")\b")

_SEEDED = ("acetaminophen", "chloronitrobenzene", "au_cyanide",
           "cr_carbonyl")


def _toolsets_for(task, tmp_path: Path) -> WorkspaceToolsets:
    auto = {g: "auto" for g in GROUPS}
    ws = workspace_for(_REPO, mode=task.mode, task_class=task.task_class)
    if ws is None:
        return WorkspaceToolsets(_REPO, is_delfin_workspace=True,
                                 settings=auto)
    dst = tmp_path / "ws"
    shutil.copytree(ws, dst)
    if (task.task_class or "") == "chemistry":
        for mol in _SEEDED:
            (dst / "chem" / mol).mkdir(parents=True, exist_ok=True)
            (dst / "chem" / mol / "start.xyz").write_text("1\n\nH 0 0 0\n")
    return WorkspaceToolsets(dst, is_delfin_workspace=False, settings=auto)


_TASKS = [t for t in B.load_tasks()
          if any(_NAME_RE.search(s.pattern or "")
                 and any(_GROUP_OF[n] for n in _NAME_RE.findall(s.pattern))
                 for s in t.expected_signals)]


def test_the_check_sees_the_tasks_it_is_about():
    assert len(_TASKS) >= 10, len(_TASKS)


@pytest.mark.parametrize("task", _TASKS, ids=[t.id for t in _TASKS])
def test_the_task_is_offered_a_tool_its_signal_names(task, tmp_path):
    ts = _toolsets_for(task, tmp_path)
    role = "office_agent" if task.mode == "office" else ""
    hidden = set(ts.hidden_groups(role=role))
    if not ts.is_delfin_workspace:
        hidden.add("delfin_only")
    for sig in task.expected_signals:
        if getattr(sig, "optional", False):
            continue
        named = set(_NAME_RE.findall(sig.pattern or ""))
        if named and all(_GROUP_OF[n] in hidden for n in named):
            pytest.fail(f"{task.id}: every tool of {sorted(named)} is in a "
                        f"hidden group {sorted(hidden)} ({ts.decided()})")
