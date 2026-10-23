"""Told not to build in the source tree, and not told where it is.

Asked to read ``delfin/agent/cli.py`` for orientation from a working
directory two levels inside the checkout, a live agent wrote "the path is
wrong, the file must be elsewhere" and then spent eighteen tool calls
looking for it, ``find /`` among them. The prompt had told it the DELFIN
source tree is a different directory and must not be built into — true,
and useless without an address.

Not naming the root never kept the agent out. It only made getting in
expensive.
"""

from pathlib import Path

import pytest

from delfin.agent.prompt_loader import PromptLoader

_REPO = Path(__file__).resolve().parents[1]


def _block(workspace: Path | None) -> str:
    loader = PromptLoader()
    loader.workspace_root = workspace
    return loader._build_session_env_block()


def test_a_working_directory_inside_the_repo_is_told_the_root():
    ws = _REPO / "tests" / "fixtures" / "user_project_workspace"
    block = _block(ws)
    assert str(ws) in block
    assert "project root" in block
    assert str(_REPO) in block
    # And the address arrives with the rule it belongs to, not instead
    # of it: naming the tree as a workspace once had the agent building
    # the user's project inside DELFIN's own checkout.
    assert "never build into it" in block
    assert "do not build the user's project inside it" in block


def test_the_root_line_names_what_a_repo_relative_path_means():
    ws = _REPO / "tests" / "fixtures" / "user_project_workspace"
    block = _block(ws)
    assert "delfin/agent/cli.py" in block
    assert "not\nunder your working directory" in block.replace("\n", "\n") \
        or "not under your working directory" in block.replace("\n", " ")


def test_a_workspace_that_is_the_root_says_nothing_extra():
    """No second address when there is only one directory in play."""
    block = _block(None)
    assert "project root" not in block
    assert "cwd:" in block


def test_a_workspace_outside_any_repo_degrades_quietly(tmp_path):
    block = _block(tmp_path)
    assert str(tmp_path) in block
    assert "project root" not in block


def test_the_block_never_raises(tmp_path):
    loader = PromptLoader()
    for ws in (None, tmp_path, tmp_path / "does-not-exist",
               _REPO / "tests" / "fixtures" / "user_project_workspace"):
        loader.workspace_root = ws
        assert isinstance(loader._build_session_env_block(), str)
