"""A sub-agent's background job was orphaned by its own worktree teardown.

``bash_background`` registers a job under the CALLER's workspace. For an
isolated sub-agent that workspace IS the temporary git worktree, so the
job record was written to ``<worktree>/.delfin/bash_jobs.json`` -- and
``.delfin/`` is git-ignored, so the tree looked unchanged and was removed
with ``git worktree remove --force``, taking the record with it.

The child process survives that: it is started in its own session, so it
kept running to its 24-hour cap with its working directory deleted
underneath it. Its completion event was destroyed with the registry file,
the workspace sweep skips a directory that no longer exists (by design),
and ``bash_status(job_id)`` answered ``unknown job_id`` -- indistinguishable
from a typo, while a live process burned the node unreaped and uncounted
against the background-job cap.

The teardown now refuses while the tree still owns running jobs and says
so, with the ids, in the sub-agent's report. Re-registering the record
against the parent workspace would have rescued the bookkeeping and still
left a running process whose cwd had been deleted.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from pathlib import Path

import pytest

from delfin.agent import bash_jobs as bj
from delfin.agent import subagents as sa
from delfin.agent import worktree as wt


@pytest.fixture(autouse=True)
def _iso(monkeypatch, tmp_path):
    monkeypatch.setattr(sa, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(sa, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(sa, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")
    monkeypatch.setattr(bj, "_INDEX_PATH", tmp_path / "jobs_index.json")


def _register(workspace: Path, job_id: str, *, pid: int,
              finished_at=None) -> None:
    bj._persist_job_start(str(workspace), {
        "job_id": job_id,
        "pid": pid,
        "proc_start_ticks": bj._proc_start_ticks(pid),
        "command": "xtb --opt mol.xyz",
        "description": "",
        "cwd": str(workspace),
        "workspace": str(workspace),
        "stdout_path": "", "stderr_path": "",
        "started_at": time.time(),
        "timeout_s": 3600,
        "exit_code": None,
        "finished_at": finished_at,
        "acknowledged": False,
    })


@pytest.fixture
def live_child():
    proc = subprocess.Popen(["sleep", "30"], stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    try:
        yield proc
    finally:
        proc.kill()
        proc.wait()


# ---------------------------------------------------------------------------
# Who still owns the directory
# ---------------------------------------------------------------------------

def test_a_live_job_is_reported_as_holding_its_workspace(tmp_path, live_child):
    ws = tmp_path / "wt"
    ws.mkdir()
    _register(ws, "ab12", pid=live_child.pid)
    held = bj.running_jobs_for_workspace(ws)
    assert [j["job_id"] for j in held] == ["ab12"]
    assert held[0]["pid"] == live_child.pid


def test_a_finished_job_holds_nothing(tmp_path, live_child):
    ws = tmp_path / "wt"
    ws.mkdir()
    _register(ws, "ab12", pid=live_child.pid, finished_at=time.time())
    assert bj.running_jobs_for_workspace(ws) == []


def test_a_dead_pid_holds_nothing(tmp_path):
    ws = tmp_path / "wt"
    ws.mkdir()
    proc = subprocess.Popen(["true"])
    proc.wait()
    _register(ws, "ab12", pid=proc.pid)
    assert bj.running_jobs_for_workspace(ws) == []


def test_an_unreadable_registry_reports_nothing(tmp_path):
    ws = tmp_path / "wt"
    (ws / ".delfin").mkdir(parents=True)
    (ws / ".delfin" / "bash_jobs.json").write_text("{not json",
                                                   encoding="utf-8")
    assert bj.running_jobs_for_workspace(ws) == []


# ---------------------------------------------------------------------------
# The teardown refuses, and says why
# ---------------------------------------------------------------------------

def _fake_info(path: Path):
    class _Info:
        branch = "agent/deadbeef"
        repo_dir = path.parent
        base_ref = "HEAD"
        cleaned_up = False
        had_changes = False
        final_path = None
    _Info.path = path
    return _Info()


def test_the_tree_is_kept_while_a_job_still_runs(tmp_path, live_child):
    tree = tmp_path / "delfin-wt-1"
    tree.mkdir()
    _register(tree, "ab12", pid=live_child.pid)

    held = sa._jobs_holding_worktree(tree)
    assert held
    summary = sa._kept_for_running_jobs(_fake_info(tree), held)
    assert summary["cleaned_up"] is False
    assert summary["running_jobs"] == ["ab12"]
    assert "NOT removed" in summary["warning"]
    assert "ab12" in summary["warning"]
    assert "bash_kill" in summary["warning"]
    assert str(tree) in summary["warning"]


def test_the_job_registry_survives_the_run(tmp_path, live_child):
    """The record is the only map back to the pid, the output files and
    the completion event."""
    tree = tmp_path / "delfin-wt-2"
    tree.mkdir()
    _register(tree, "ab12", pid=live_child.pid)
    registry = tree / ".delfin" / "bash_jobs.json"
    assert registry.is_file()

    # What the teardown decision does when jobs are still running: nothing.
    assert sa._jobs_holding_worktree(tree)
    assert registry.is_file()
    assert json.loads(registry.read_text())["jobs"]["ab12"]["pid"] \
        == live_child.pid


def test_a_tree_nobody_is_using_is_still_removed(tmp_path):
    """The refusal must not turn every isolated run into a leaked folder."""
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.email", "a@b.c"],
                   check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "t"],
                   check=True)
    (repo / "f.txt").write_text("x", encoding="utf-8")
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "init"],
                   check=True)

    info = wt.enter_worktree(repo, parent=tmp_path)
    assert sa._jobs_holding_worktree(info.path) == []
    wt.exit_worktree(info, keep_if_changed=True)
    assert not info.path.exists()


def test_the_kept_tree_is_named_in_the_stored_report(tmp_path, live_child):
    """A BACKGROUND run answers with an id, not a payload, so the account
    of the tree has to reach the parent through the session store."""
    tree = tmp_path / "delfin-wt-3"
    tree.mkdir()
    _register(tree, "ab12", pid=live_child.pid)
    summary = sa._kept_for_running_jobs(
        _fake_info(tree), sa._jobs_holding_worktree(tree))

    sa._save_subagent_session(
        "bg1", subagent_type="general-purpose", description="run xtb",
        messages=[{"role": "assistant", "content": "started the optimisation"}],
        interactions=[], worktree=summary)
    collected = sa.get_subagent_result("bg1")
    assert "NOT removed" in collected["worktree"]["warning"]


def _git_repo(path: Path) -> Path:
    path.mkdir()
    subprocess.run(["git", "init", "-q", str(path)], check=True)
    subprocess.run(["git", "-C", str(path), "config", "user.email", "a@b.c"],
                   check=True)
    subprocess.run(["git", "-C", str(path), "config", "user.name", "t"],
                   check=True)
    (path / "f.txt").write_text("x", encoding="utf-8")
    subprocess.run(["git", "-C", str(path), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(path), "commit", "-qm", "init"],
                   check=True)
    return path


def test_a_whole_isolated_run_keeps_the_tree_its_job_lives_in(
        tmp_path, live_child, monkeypatch):
    """End to end: the sub-agent starts a job in its worktree, the run ends,
    and the tree with the job's registry is still there."""
    from delfin.agent.api_client import StreamEvent, create_client

    repo = _git_repo(tmp_path / "repo")
    parent = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(repo))

    def _fake_stream(**kwargs):
        # The sub-agent's workspace is named in its own system prompt.
        line = [ln for ln in str(kwargs.get("system", "")).splitlines()
                if ln.startswith("Workspace: ")][0]
        _register(Path(line.split("Workspace: ", 1)[1].strip()), "ab12",
                  pid=live_child.pid)
        yield StreamEvent(type="text_delta", text="started the optimisation")

    monkeypatch.setattr(parent, "stream_message", _fake_stream)
    res = sa.run_subagent(
        subagent_type="general-purpose", description="run xtb",
        prompt="optimise it", parent_client=parent,
        parent_perms=parent._permissions, isolation="worktree")

    tree = Path(res.workspace)
    try:
        assert tree != repo, "the run was not isolated, so this proves nothing"
        assert tree.is_dir(), "the tree was removed while its job was running"
        assert (tree / ".delfin" / "bash_jobs.json").is_file()
        assert res.worktree["cleaned_up"] is False
        assert "ab12" in res.worktree["warning"]
    finally:
        # The tree is kept ON PURPOSE here, so the test has to clear it.
        shutil.rmtree(tree, ignore_errors=True)


def test_an_isolated_run_without_jobs_still_cleans_up(tmp_path, monkeypatch):
    from delfin.agent.api_client import StreamEvent, create_client

    repo = _git_repo(tmp_path / "repo")
    parent = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(repo))

    def _fake_stream(**kwargs):
        yield StreamEvent(type="text_delta", text="looked around")

    monkeypatch.setattr(parent, "stream_message", _fake_stream)
    res = sa.run_subagent(
        subagent_type="general-purpose", description="look",
        prompt="look", parent_client=parent,
        parent_perms=parent._permissions, isolation="worktree")
    assert res.worktree.get("cleaned_up") is True
    assert not Path(res.workspace).exists()


def test_the_source_registers_jobs_under_the_workspace(tmp_path):
    """Where the orphaning starts: the job's home IS the sub-agent's tree."""
    from delfin.agent import api_client as A
    source = Path(A.__file__).read_text(encoding="utf-8")
    block = source[source.index("_bj.get_registry().start("):]
    assert "workspace=str(perms.workspace)" in block[:600]
