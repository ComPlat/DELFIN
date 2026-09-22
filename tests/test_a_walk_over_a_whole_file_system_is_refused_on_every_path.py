"""A walk over a whole file system is refused -- on every path, the same way.

Rooted at a home directory or a mount point, find/du/tree/rg, grep -r and
ls -R cost one metadata request per file on a shared parallel file system;
a site's operations team reads that as abuse (it once did, for a periodic
du over a home). The rule existed only in the CLI backend's approval
runner. The API and terminal sessions -- most of the work -- walked past
it: `find ~ -name '*.py'` ran three times in one session, and walks over
the archive ran unasked (measured in the audit log, 2026-09-22).

Both paths now call one function. Every command below is asked of both,
and both must give the same answer.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import sandbox as S


@pytest.fixture
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    (h / "proj" / "calc").mkdir(parents=True)
    (h / "notes.txt").write_text("x\n")
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: h))
    return h


def _walks(home):
    return [f"find {home} -name '*.py'", "find ~ -name '*.out'",
            f"du -sh {home}", "du -sh ~", f"tree {home}", f"rg TODO {home}",
            f"grep -rn TODO {home}", f"grep -R x ~", f"ls -R {home}",
            f"ls -laR {home}", "find / -name passwd", f"nice -n 10 du -s {home}",
            f"cd {home}/proj && du -sh {home} | tail -1"]


def _fine(home):
    return [f"find {home}/proj -name '*.py'", f"du -sh {home}/proj/calc",
            f"grep -rn TODO {home}/proj", f"ls -la {home}",
            f"grep -n x {home}/notes.txt", f"rg TODO {home}/proj",
            f"tree {home}/proj"]


def test_every_walk_is_refused_by_the_shared_rule(home):
    for cmd in _walks(home):
        assert S.tree_walk_refusal(cmd), cmd


def test_a_subdirectory_or_a_non_recursive_read_is_not_a_walk(home):
    for cmd in _fine(home):
        assert S.tree_walk_refusal(cmd) is None, cmd


def test_the_cli_path_and_the_api_gate_give_the_same_answer(home):
    perms = A.KitToolPermissions(workspace=str(home / "proj"),
                                 mode="bypassPermissions")
    for cmd in _walks(home):
        cli = S.is_allowed(cmd)
        gate = A._doc_executor._run_permission_gate(
            "bash", {"command": cmd}, perms)
        # The CLI runner may refuse earlier, on its own allow-list (nice is
        # not on it); what matters is that neither path lets a walk run.
        assert not cli.allowed, (cmd, cli)
        assert gate and "walk" in gate, (cmd, gate)


def test_even_bypass_mode_does_not_open_it(home):
    perms = A.KitToolPermissions(workspace=str(home / "proj"),
                                 mode="bypassPermissions")
    err = A._doc_executor._run_permission_gate(
        "bash", {"command": "du -sh ~"}, perms)
    assert err and "delfin.quota.home_usage" in err


def test_relative_paths_are_read_from_where_the_command_runs(home):
    # From the project, "." is the project: not a walk.
    assert S.tree_walk_refusal("find . -name '*.py'", cwd=home / "proj") is None
    assert S.tree_walk_refusal("grep -rn TODO", cwd=home / "proj") is None
    # From the home, it is.
    assert S.tree_walk_refusal("find . -name '*.py'", cwd=home)
    # And a cd inside the command moves it.
    assert S.tree_walk_refusal("cd ~ && find . -name x", cwd=home / "proj")
    assert S.tree_walk_refusal("cd .. && cd .. && du -sh .",
                               cwd=home / "proj" / "calc")


def test_the_gate_reads_relative_paths_from_the_workspace(home):
    perms = A.KitToolPermissions(workspace=str(home / "proj"),
                                 mode="bypassPermissions")
    gate = A._doc_executor._run_permission_gate
    assert gate("bash", {"command": "find . -name '*.py'"}, perms) is None
    err = gate("bash", {"command": "cd ~ && find . -name x"}, perms)
    assert err and "walk" in err
