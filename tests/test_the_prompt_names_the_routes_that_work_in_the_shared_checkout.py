"""What the agent is told to run works in the checkout DELFIN runs from.

The shared checkout stays on main and refuses a branch move
(shared_checkout.refusal). Every suggestion in the prompts, skills and
gate messages to make a branch with ``git switch -c`` / ``git checkout -b``
therefore names the route that works there -- ``enter_worktree`` (and
``git branch <name> origin/main`` without a switch, or the dashboard's Test
branch) -- or is gone. Checked in the BUILT solo and dashboard prompts, not
only in the source files: a rule is a rule when it is in the prompt.

The approved shape ``git switch -c pr/<name> origin/main`` stays approved
inside a worktree and is refused in the shared checkout; its doc says so.
"""

from __future__ import annotations

import inspect
import os
import re
import subprocess
from pathlib import Path

import pytest

from delfin.agent import held_branch, read_only_shapes
from delfin.agent import shared_checkout as SC

PACK = Path(read_only_shapes.__file__).resolve().parent / "pack"
_MOVE = re.compile(r"git (switch -c|checkout -b)")
_ROUTE = re.compile(r"enter_worktree|Test branch|git branch <")
_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _unrouted(text: str, window: int = 400) -> list[str]:
    """Branch-making suggestions with no working route within *window*."""
    flat = re.sub(r"\s+", " ", text)
    bad = []
    for m in _MOVE.finditer(flat):
        near = flat[max(0, m.start() - window):m.end() + window]
        if not _ROUTE.search(near):
            bad.append(near[window - 80:window + 120])
    return bad


@pytest.mark.parametrize("role,mode", [("solo_agent", "solo"),
                                       ("dashboard_agent", "dashboard")])
def test_the_built_prompt_routes_every_branch_suggestion(role, mode):
    from delfin.agent.prompt_loader import PromptLoader
    built = PromptLoader().build_system_prompt(
        role_id=role, mode_id=mode, task_text="open a pull request for DELFIN",
        session_key=f"routes-{role}")
    assert _unrouted(built) == []
    if role == "solo_agent":
        assert 'enter_worktree(base_ref="origin/main")' in built
        assert not re.search(r"git switch -c\s+<user>/<feature>", built)


@pytest.mark.parametrize("name", ["skills/new-branch.md", "skills/git-errors.md",
                                  "skills/git-merge-conflicts.md",
                                  "shared/git_workflow_addendum.md",
                                  "agents/solo_agent.md"])
def test_the_pack_text_routes_every_branch_suggestion(name):
    assert _unrouted((PACK / name).read_text(encoding="utf-8")) == []


def test_the_gate_messages_route_every_branch_suggestion():
    from delfin.agent import api_client
    src = inspect.getsource(api_client._DocToolExecutor._run_permission_gate) \
        if hasattr(api_client, "_DocToolExecutor") else \
        inspect.getsource(api_client)
    assert _unrouted(src) == []


def test_held_branch_names_the_commits_without_switching():
    src = inspect.getsource(held_branch.detach_at_branch)
    assert "git switch -c" not in src and "git branch <name>" in src


def test_the_shape_doc_says_where_it_is_refused():
    doc = read_only_shapes.__doc__ or ""
    assert "shared checkout" in doc and "worktree" in doc


@pytest.fixture
def shared(tmp_path, monkeypatch):
    d = tmp_path / "shared"
    d.mkdir()
    for args in (("init", "-q", "-b", "main"),):
        subprocess.run(["git", *args], cwd=d, env=_ENV, check=True)
    (d / "f").write_text("x")
    subprocess.run(["git", "-c", "user.email=t@t", "-c", "user.name=t", "add", "f"],
                   cwd=d, env=_ENV, check=True)
    subprocess.run(["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit",
                    "-qm", "f"], cwd=d, env=_ENV, check=True)
    monkeypatch.setenv(SC.ENV, str(d))
    return d


def test_the_pr_shape_is_refused_in_the_shared_checkout_only(shared, tmp_path):
    cmd = "git switch -c pr/fix-1 origin/main"
    assert SC.refusal(cmd, shared)
    wt = tmp_path / "wt"
    subprocess.run(["git", "worktree", "add", "-q", "--detach", str(wt), "main"],
                   cwd=shared, env=_ENV, check=True)
    assert SC.refusal(cmd, wt) == ""
