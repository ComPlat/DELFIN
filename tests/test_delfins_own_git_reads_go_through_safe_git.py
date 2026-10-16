"""DELFIN's own git commands go through safe_git.

DELFIN's code runs git from outside the sandbox on repositories it does
not control: the session's workspace, a worktree the agent worked in, the
DELFIN checkout itself when that is the workspace. A plain ``subprocess``
call with ``"git"`` as the program runs whatever the repository's config
and attributes name; ``safe_git.run`` does not
(tests/test_git_runs_no_program_the_repository_names.py).

The guard reads the syntax tree of every module under ``delfin/`` and
reports each git command line written as a literal (a list or tuple whose
first element is ``"git"``, or a ``"git ..."`` string handed to a
subprocess runner) outside the reviewed exceptions in ``ALLOWED``. A
literal that is only compared with (``w[:2] == ["git", "apply"]``,
``name in ("git", "uv")``) is a pattern, not a command, and is not
reported. A command line assembled from a variable program name is not
seen; that is the remaining blind spot.
"""
import ast
import functools
import warnings
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_ROOT = _REPO / "delfin"
_RUNNERS = {"run", "check_output", "check_call", "call", "Popen",
            "getoutput", "getstatusoutput", "system",
            "create_subprocess_exec", "create_subprocess_shell"}

#: ``path:function`` -> why that git command line does not go through
#: safe_git. Reviewed; a new entry needs the same scrutiny.
ALLOWED = {
    "delfin/agent/safe_git.py:_base":
        "the runner itself",
    "delfin/agent/git_config_watch.py:_program_entries":
        "'git config --file <temporary copy> --list' parses bytes; no "
        "repository configuration is read and no program runs",
    "delfin/agent/api_client.py:_github_push_plan":
        "the push the agent asked for, run by _run_github_push outside the "
        "cage: hooks and fsmonitor off, the transport (ssh options or the gh "
        "credential) fixed on the command line, refused when the checkout "
        "configures its own transport; the user's credentials must work",
    "delfin/agent/api_client.py:_git_value_runs":
        "a word list handed to the bash gate's analyser; nothing runs",
    "delfin/agent/doctor.py:_check_push":
        "ls-remote from an empty temporary directory with discovery stopped "
        "above it (no repository configuration applies) and in its own "
        "session so ssh cannot prompt on the terminal; and 'config --get "
        "credential.helper', which reads configuration and runs nothing -- "
        "through safe_git the override would hide the checkout's own "
        "helper, the one this row reports",
    "delfin/agent/pack/benchmark/setup/a_checkout_whose_main_is_held.py:git":
        "benchmark setup: builds the fixture repository it creates, before "
        "any model command runs in it",
    "delfin/agent/pack/benchmark/setup/a_fix_claim_repo.py:_git":
        "benchmark setup: builds its own fixture repository",
    "delfin/agent/pack/benchmark/setup/a_fix_claim_repo.py:_head_test_on_base_code":
        "benchmark setup: builds its own fixture repository",
    "delfin/agent/pack/benchmark/setup/a_repo_whose_test_was_already_broken.py:git":
        "benchmark setup: builds its own fixture repository",
    "delfin/agent/pack/benchmark/setup/generalist_large_repo.py:build":
        "benchmark setup: builds its own fixture repository",
    "delfin/agent/pack/benchmark/setup/generalist_merge_conflict.py:build":
        "benchmark setup: builds its own fixture repository",
    "delfin/agent/pack/benchmark/setup/generalist_repo.py:git":
        "benchmark setup: builds its own fixture repository",
}


def _is_git_head(node) -> bool:
    if not isinstance(node, (ast.List, ast.Tuple)) or not node.elts:
        return False
    head = node.elts[0]
    return isinstance(head, ast.Constant) and head.value == "git"


def git_command_lines(path: Path) -> list[tuple[int, str]]:
    """``(line, enclosing function)`` of each git command line in *path*."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        warnings.simplefilter("ignore", DeprecationWarning)
        tree = ast.parse(path.read_text(encoding="utf-8"))
    parent = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parent[child] = node

    def function_of(node) -> str:
        while node in parent:
            node = parent[node]
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                return node.name
        return "<module>"

    found = []
    for node in ast.walk(tree):
        if _is_git_head(node):
            if isinstance(parent.get(node), ast.Compare):
                continue
            found.append((node.lineno, function_of(node)))
        elif isinstance(node, ast.Call) and node.args:
            f = node.func
            name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")
            first = node.args[0]
            if (name in _RUNNERS and isinstance(first, ast.Constant)
                    and isinstance(first.value, str)
                    and first.value.lstrip().startswith("git ")):
                found.append((node.lineno, function_of(node)))
    return sorted(found)


@functools.lru_cache(maxsize=1)
def _all_sites() -> dict[str, list[int]]:
    sites: dict[str, list[int]] = {}
    for path in sorted(_ROOT.rglob("*.py")):
        rel = path.relative_to(_REPO).as_posix()
        for line, func in git_command_lines(path):
            sites.setdefault(f"{rel}:{func}", []).append(line)
    return sites


def test_every_git_command_line_in_delfin_goes_through_safe_git():
    plain = {k: v for k, v in _all_sites().items() if k not in ALLOWED}
    assert not plain, (
        "git run directly (use delfin.agent.safe_git.run, or add a reviewed "
        "entry to ALLOWED with the reason):\n"
        + "\n".join(f"  {k} lines {v}" for k, v in sorted(plain.items())))


def test_every_exception_still_names_a_git_command_line():
    stale = sorted(set(ALLOWED) - set(_all_sites()))
    assert not stale, f"ALLOWED entries with no git command line left: {stale}"


@pytest.mark.parametrize("src,expected", [
    ("import subprocess\nsubprocess.run(['git', 'status'])\n", [(2, "<module>")]),
    ("def f():\n    argv = ['git', 'log']\n    return argv\n", [(2, "f")]),
    ("import os\nos.system('git status')\n", [(2, "<module>")]),
    ("ok = w[:2] == ['git', 'apply']\n", []),
    ("ok = name in ('git', 'uv')\n", []),
    ("from delfin.agent import safe_git\nsafe_git.run(['status'], '.')\n", []),
])
def test_the_guard_tells_a_command_from_a_pattern(tmp_path, src, expected):
    p = tmp_path / "m.py"
    p.write_text(src)
    assert git_command_lines(p) == expected
