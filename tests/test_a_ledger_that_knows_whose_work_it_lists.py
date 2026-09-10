"""The changes report lists what THIS workspace did, and the worktree
tools end their own sequence without an error.

Both found by driving the tools the way a model does, through the
executor, in a throwaway repository:

  list_changes_made in a fresh probe workspace listed a benchmark's
  commands -- `python3 pipeline.py`, `rm _verify_boltzmann.py` -- as
  "everything recorded under <probe>". The audit record of a command
  carried a cwd only when the model had passed one, so nearly every
  record had none, and a record without an absolute cwd passes the
  workspace filter. Records now name their workspace, and the report
  judges by it exactly.

  enter -> merge -> exit, the order the descriptions suggest, ended in
  {"error": "worktree path missing"}: a clean merge removes the
  worktree itself and said nothing about it. The merge now says so, and
  exit answers a path that is already gone with ok/removed=false.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from delfin.agent import audit_log
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------

def _rec(**kw):
    """A command record in the shape the executor writes: the command sits
    at the top level, which is where the report reads it."""
    base = {"tool": "bash", "decision": "ok", "command": kw.pop("command", "ls")}
    base.update(kw)
    return base


def test_a_record_from_another_workspace_is_not_listed(tmp_path):
    log = tmp_path / "audit.log"
    mine = tmp_path / "mine"; other = tmp_path / "other"
    mine.mkdir(); other.mkdir()
    audit_log.append(_rec(command="python3 pipeline.py", workspace=str(other), cwd=str(other)), log_path=log)
    audit_log.append(_rec(command="ls -la", workspace=str(mine), cwd=str(mine)), log_path=log)
    rep = audit_log.build_changes_report(None, log_path=log, workspace=str(mine))
    cmds = [c["command"] for c in rep["commands"]]
    assert cmds == ["ls -la"], cmds


def test_a_tagged_record_is_judged_by_its_tag_even_with_a_relative_cwd(tmp_path):
    """The old filter let a relative cwd through; the tag decides now."""
    log = tmp_path / "audit.log"
    mine = tmp_path / "mine"; other = tmp_path / "other"
    mine.mkdir(); other.mkdir()
    audit_log.append(_rec(command="rm x", workspace=str(other), cwd="."), log_path=log)
    rep = audit_log.build_changes_report(None, log_path=log, workspace=str(mine))
    assert rep["commands"] == []


def test_an_untagged_old_record_still_falls_back_to_the_path_test(tmp_path):
    log = tmp_path / "audit.log"
    mine = tmp_path / "mine"; other = tmp_path / "other"
    mine.mkdir(); other.mkdir()
    audit_log.append(_rec(command="old absolute elsewhere", cwd=str(other)), log_path=log)
    audit_log.append(_rec(command="old relative", cwd="."), log_path=log)
    rep = audit_log.build_changes_report(None, log_path=log, workspace=str(mine))
    cmds = [c["command"] for c in rep["commands"]]
    assert cmds == ["old relative"], cmds


def test_the_same_workspace_spelled_differently_still_matches(tmp_path):
    log = tmp_path / "audit.log"
    mine = tmp_path / "mine"; mine.mkdir()
    audit_log.append(_rec(command="pwd", workspace=str(mine) + "/."), log_path=log)
    rep = audit_log.build_changes_report(None, log_path=log, workspace=str(mine))
    assert [c["command"] for c in rep["commands"]] == ["pwd"]


# ---------------------------------------------------------------------------
# What the executor records
# ---------------------------------------------------------------------------

def test_a_command_run_through_the_executor_records_workspace_and_absolute_cwd(tmp_path, monkeypatch):
    log = tmp_path / "audit.log"
    monkeypatch.setattr(audit_log, "_default_log_path", lambda: log)
    ws = tmp_path / "ws"; ws.mkdir()
    ex = _DocToolExecutor()
    out = ex.execute("bash", {"command": "pwd"}, KitToolPermissions(workspace=str(ws)))
    assert json.loads(out).get("exit_code") == 0
    recs = [json.loads(l) for l in log.read_text().splitlines() if l.strip()]
    bash = [r for r in recs if r.get("tool") == "bash"]
    assert bash, recs
    rec = bash[-1]
    assert Path(rec["workspace"]).resolve() == ws.resolve()
    assert Path(rec["cwd"]).is_absolute() and Path(rec["cwd"]).resolve() == ws.resolve()


def test_a_probe_workspace_does_not_inherit_another_workspaces_commands(tmp_path, monkeypatch):
    """The end-to-end shape of the finding: two workspaces, one report each."""
    log = tmp_path / "audit.log"
    monkeypatch.setattr(audit_log, "_default_log_path", lambda: log)
    a = tmp_path / "a"; b = tmp_path / "b"; a.mkdir(); b.mkdir()
    ex = _DocToolExecutor()
    ex.execute("bash", {"command": "echo from-a"}, KitToolPermissions(workspace=str(a)))
    ex.execute("bash", {"command": "echo from-b"}, KitToolPermissions(workspace=str(b)))
    report_a = ex.execute("list_changes_made", {}, KitToolPermissions(workspace=str(a)))
    assert "echo from-a" in report_a
    assert "echo from-b" not in report_a


# ---------------------------------------------------------------------------
# The worktree sequence
# ---------------------------------------------------------------------------

@pytest.fixture
def repo(tmp_path):
    r = tmp_path / "repo"; r.mkdir()
    subprocess.run(["git", "init", "-q", "-b", "main", str(r)], check=True)
    subprocess.run(["git", "-C", str(r), "config", "user.email", "t@t"], check=True)
    subprocess.run(["git", "-C", str(r), "config", "user.name", "t"], check=True)
    (r / "README.md").write_text("hello\n")
    subprocess.run(["git", "-C", str(r), "add", "."], check=True)
    subprocess.run(["git", "-C", str(r), "commit", "-q", "-m", "init"], check=True)
    return r


def test_enter_merge_exit_ends_without_an_error(repo):
    ex = _DocToolExecutor()
    perms = KitToolPermissions(workspace=str(repo))
    entered = json.loads(ex.execute("enter_worktree", {"repo_dir": str(repo), "branch_prefix": "probe"}, perms))
    wt = entered["path"]
    (Path(wt) / "note.txt").write_text("from the worktree\n")

    merged = json.loads(ex.execute("worktree_merge", {"path": wt}, perms))
    assert merged["status"] == "ok" and merged["applied"] is True
    assert merged["worktree_removed"] is True
    assert "removed" in merged["message"]
    assert not Path(wt).exists()
    assert (repo / "note.txt").read_text() == "from the worktree\n"

    left = json.loads(ex.execute("exit_worktree", {"path": wt}, perms))
    assert "error" not in left, left
    assert left["status"] == "ok" and left["removed"] is False
    assert "nothing to tear down" in left["note"]


def test_exit_on_a_live_worktree_reports_that_it_removed_it(repo):
    ex = _DocToolExecutor()
    perms = KitToolPermissions(workspace=str(repo))
    entered = json.loads(ex.execute("enter_worktree", {"repo_dir": str(repo), "branch_prefix": "probe"}, perms))
    wt = entered["path"]
    left = json.loads(ex.execute("exit_worktree", {"path": wt}, perms))
    assert left["status"] == "ok" and left["removed"] is True
    assert not Path(wt).exists()


def test_the_descriptions_say_who_removes_the_worktree():
    from delfin.agent import api_client as a
    src = open(a.__file__).read()
    assert "apply removes the worktree and its branch" in src
    assert "(worktree_removed; no exit needed)" in src


# ---------------------------------------------------------------------------
# A notebook is cells, not one JSON line
# ---------------------------------------------------------------------------

def test_find_references_reads_a_notebook_by_cell(tmp_path):
    """Driven: a match in run.ipynb came back as line 1, col 67 with a
    preview of raw JSON. True, and useless to anyone reading the notebook."""
    from delfin.agent.code_nav import find_references

    nb = {"cells": [
        {"cell_type": "markdown", "source": ["# Energies\n"]},
        {"cell_type": "code", "source": ["E = -76.4\n", "print(E)\n"], "outputs": []},
        {"cell_type": "code", "source": "E * 2\n", "outputs": []},
    ], "nbformat": 4, "nbformat_minor": 5, "metadata": {}}
    (tmp_path / "run.ipynb").write_text(json.dumps(nb))
    out = find_references(tmp_path, "E", language="any")
    hits = [(m["line"], m["preview"]) for m in out["matches"] if m["path"] == "run.ipynb"]
    assert (1, "cell 1: E = -76.4") in hits
    assert (2, "cell 1: print(E)") in hits
    assert (1, "cell 2: E * 2") in hits
    assert all(not p.startswith("{") for _, p in hits), "raw JSON leaked into a preview"


def test_a_broken_notebook_falls_back_to_plain_lines(tmp_path):
    from delfin.agent.code_nav import find_references

    (tmp_path / "broken.ipynb").write_text("not json but mentions E here\n")
    out = find_references(tmp_path, "E", language="any")
    assert any(m["path"] == "broken.ipynb" and m["line"] == 1 for m in out["matches"])
