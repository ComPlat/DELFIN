"""Two anti-drift guards that only ever fired on one backend.

``_maybe_pin_project_dir`` records the directory of the first project
write so every later prompt can re-pin the agent there.
``_note_stray_write`` says so when a write lands somewhere else. Both
read the tool name off a tool_use event and match it against one set —
which held only the OpenAI-compatible spellings, so on the CLI backend
neither had ever fired.

The mirror image of the same fault was three lines further down: the
auto-verification hook matched only {"Edit", "Write"} and so could never
fire on the KIT and Ollama backends. That one turned out to have no
caller at all and is gone; this one is real, and now reads both.
"""

from __future__ import annotations

import pytest

from delfin.agent.engine import AgentEngine


@pytest.mark.parametrize("name", [
    "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
    "Write", "Edit", "MultiEdit", "NotebookEdit",
])
def test_every_spelling_of_a_write_counts_as_one(name):
    assert name in AgentEngine._MUTATE_TOOLS_FOR_PIN, name


@pytest.mark.parametrize("name", ["read_file", "Read", "grep_file", "Grep",
                                  "bash", "Bash", "search_docs"])
def test_reading_is_not_a_write(name):
    assert name not in AgentEngine._MUTATE_TOOLS_FOR_PIN, name


def test_the_pin_survives_the_mcp_namespace(tmp_path):
    """Both guards strip the prefix before matching, and the KIT backend
    is where the namespaced names come from."""
    from unittest.mock import MagicMock, patch
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                          model="kit.glm-5.3", mode="solo")
    eng._maybe_pin_project_dir(
        "mcp__kit-coding__write_file",
        '{"path": "/work/project/app.py"}')
    assert eng._project_dir == "/work/project", eng._project_dir


def test_the_cli_spelling_pins_the_same_directory(tmp_path):
    from unittest.mock import MagicMock, patch
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = AgentEngine(repo_dir=tmp_path, backend="cli", mode="solo")
    eng._maybe_pin_project_dir("Write", '{"file_path": "/work/project/app.py"}')
    assert eng._project_dir == "/work/project", eng._project_dir
