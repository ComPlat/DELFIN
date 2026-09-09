"""Both paths were named and the model still could not put them together.

The session-environment block already told the agent two true things: its
working directory, and the project root above it that a repo-relative
path resolves under. What it never said is that ONE repo-relative path --
the working directory's own -- means here.

That is the whole content of the second-largest group in the denial log.
Measured over ~/.delfin/audit.log on 2026-09-09: 27 refusals naming a
path outside the workspace roots and 9 more `cwd is not a directory`, and
every one of them was the working directory's repo-relative spelling,
because the person who wrote the request was standing in the root:

    Bau mir in tests/fixtures/user_project_workspace/ ein Skript …

The agent IS that directory. Ten of the thirteen went further and passed
the project root itself as an absolute cwd. Neither is the agent getting
lost -- it is standing in the right place and describing it from where
the request was written, which is what a person does too.

Two things changed. The refusals name the mistake and the path that works
(api_client._hint_if_it_names_the_workspace), and the environment block
names the path before anyone gets it wrong. This file is about the
second: a rule is a rule when it is in the BUILT prompt, so it is checked
there and not in the source that produces it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.prompt_loader import PromptLoader


@pytest.fixture
def root(tmp_path):
    """A git repo with a project directory inside it."""
    import subprocess
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg" / "core.py").write_text("x = 1\n", encoding="utf-8")
    ws = tmp_path / "projects" / "bookmarks"
    ws.mkdir(parents=True)
    return tmp_path, ws


def _block(repo: Path, workspace: Path | None) -> str:
    pl = PromptLoader(repo)
    pl.workspace_root = workspace
    return pl._build_session_env_block()


def test_the_workspace_is_named_the_way_the_request_names_it(root):
    repo, ws = root
    block = _block(repo, ws)
    # The exact string a request written from the root would use.
    assert "projects/bookmarks/" in block
    assert "IS your working directory" in block


def test_it_says_what_to_do_instead(root):
    """Naming the ambiguity without naming the fix is how the project-root
    line itself once cost eighteen tool calls."""
    repo, ws = root
    block = _block(repo, ws)
    assert "relative paths" in block
    assert "cwd" in block


def test_the_project_root_line_is_still_there(root):
    """The new line is an exception to that one, not a replacement: a
    repo-relative path that is NOT the workspace still resolves above."""
    repo, ws = root
    block = _block(repo, ws)
    assert str(repo) in block
    assert "never build into it" in block


def test_nothing_is_said_when_the_workspace_is_the_root(root):
    """No ambiguity, no sentence. A line that fires for everyone is a line
    everyone learns to skip."""
    repo, _ = root
    block = _block(repo, repo)
    assert "IS your working directory" not in block
    assert "…except" not in block


def test_nothing_is_said_when_the_workspace_is_outside_the_repo(tmp_path):
    """Then there is no repo-relative spelling of it to confuse anyone
    with, and claiming one would be false."""
    import subprocess
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    block = _block(repo, outside)
    assert "IS your working directory" not in block


def test_the_block_survives_a_workspace_that_is_not_a_repo(tmp_path):
    """Never raises: the environment block is built for every session."""
    ws = tmp_path / "plain"
    ws.mkdir()
    _block(tmp_path, ws)
