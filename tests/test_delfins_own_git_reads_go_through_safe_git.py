"""DELFIN's own git reads on the workspace go through safe_git.

Every module listed here runs git on the session's workspace from outside
the sandbox. A plain ``subprocess`` call with ``"git"`` as the program runs
whatever the repository's config and attributes name; ``safe_git.run``
does not (tests/test_git_runs_no_program_the_repository_names.py). This
guard finds a plain call in these modules by reading their syntax tree.
Further modules move onto safe_git in later changes and join the list.
"""
import ast
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1] / "delfin" / "agent"
_MODULES = ("launch_guard", "prompt_loader", "project_introspect",
            "evidence_freshness", "held_branch", "worktree", "egress_proxy")
_RUNNERS = {"run", "check_output", "check_call", "call", "Popen"}


def _plain_git_calls(path: Path) -> list[int]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    lines = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not node.args:
            continue
        f = node.func
        name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")
        if name not in _RUNNERS:
            continue
        first = node.args[0]
        if isinstance(first, (ast.List, ast.Tuple)) and first.elts:
            head = first.elts[0]
            if isinstance(head, ast.Constant) and head.value == "git":
                lines.append(node.lineno)
    return lines


@pytest.mark.parametrize("module", _MODULES)
def test_the_module_runs_no_plain_git(module):
    found = _plain_git_calls(_ROOT / f"{module}.py")
    assert not found, (
        f"delfin/agent/{module}.py runs git directly at lines {found}; "
        "use delfin.agent.safe_git.run")


def test_the_guard_sees_a_plain_call(tmp_path):
    p = tmp_path / "m.py"
    p.write_text("import subprocess\nsubprocess.run(['git', 'status'])\n")
    assert _plain_git_calls(p) == [2]
