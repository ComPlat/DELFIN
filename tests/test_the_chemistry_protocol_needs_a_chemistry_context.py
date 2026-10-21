"""The chemistry protocol is injected only where search_docs is offered.

Input: a PromptLoader with ``is_delfin_workspace`` set as the engine sets
it, and a task text. Output: the ``--- Chemistry Protocol ---`` block.

The block orders "Use search_docs (ORCA manual, xTB docs) BEFORE
answering". In a foreign workspace ``search_docs`` is stripped from the
tool list (``_DELFIN_ONLY_TOOL_NAMES``), and the keyword classifier
reads "method", "functional", "scan", "spin" and "charge" as chemistry,
so an ordinary coding task in a TypeScript repository was told to call
a tool it does not have. Measured on main: 4 of 4 foreign probe repos
carried the block for "Refactor this method into a pure function and
update the functional component".
"""
from __future__ import annotations

import pytest

from delfin.agent.prompt_loader import PromptLoader

CODING = ("Refactor this method into a pure function and update the "
          "functional component")
CHEMISTRY = "Which DFT functional and basis set should the ORCA optimization use?"


def _loader(tmp_path, delfin_ws):
    loader = PromptLoader()
    loader.workspace_root = tmp_path
    loader.is_delfin_workspace = delfin_ws
    return loader


@pytest.mark.parametrize("task", [
    CODING,
    "Add a spin loader to the React page",
    "Make the functional component handle the scan button",
    CHEMISTRY,
])
def test_a_foreign_workspace_never_gets_the_protocol(tmp_path, task):
    assert _loader(tmp_path, False)._load_chemistry_reminder(task) == ""


def test_the_built_foreign_prompt_has_no_protocol(tmp_path):
    prompt = _loader(tmp_path, False).build_system_prompt(
        role_id="solo_agent", mode_id="solo", task_text=CODING)
    assert "Chemistry Protocol" not in prompt
    assert "search_docs (ORCA manual" not in prompt


def test_a_delfin_workspace_keeps_it_for_a_chemistry_task(tmp_path):
    """Control: the protocol still fires where its tool exists."""
    text = _loader(tmp_path, True)._load_chemistry_reminder(CHEMISTRY)
    assert "search_docs" in text
