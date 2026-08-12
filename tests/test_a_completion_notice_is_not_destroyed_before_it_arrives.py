"""The one notice a twelve-hour calculation gets was destroyed three ways.

ACKNOWLEDGED BEFORE DELIVERED. ``drain_finished_events`` marked each event
acknowledged and wrote the file BEFORE the caller had done anything with the
returned list, and nothing ever un-acknowledged. A turn that died between
the drain and the prompt it was being built into took the only completion
notice of a multi-hour run with it -- permanently, and without a trace.

LOST UPDATES BETWEEN PROCESSES. The registry lock was a ``threading.Lock``:
in-process only. Two front ends sharing the per-user locator index raced
each other's read-modify-write, and the loser's entry vanished -- which is
how a live job becomes unaddressable, with ``bash_status`` answering
"unknown job_id" about a process that is running right now.

A DELETED WORKSPACE TOOK THE EVENTS WITH IT. The "jobs holding this
worktree" check existed in the sub-agent teardown only; the ``exit_worktree``
tool and the post-merge cleanup still force-removed the tree with no check
at all -- deleting the working directory of a live process along with the
registry file that was the only record of it.

The shared-node cost is the same in all three: the calculation keeps its
cores either way, and losing the record is what makes it unkillable and
uncounted. A job nobody can name is a job nobody can stop, on a machine
whose CPU belongs to everyone using it.

Real short-lived subprocesses; no SLURM, no LLM.
"""

from __future__ import annotations

import json
import multiprocessing
import subprocess
import time
from pathlib import Path

import pytest

from delfin.agent import bash_jobs as BJ
from delfin.agent import worktree as WT


@pytest.fixture(autouse=True)
def _isolated(tmp_path, monkeypatch):
    monkeypatch.setattr(BJ, "_INDEX_PATH", tmp_path / "index.json")
    BJ._REGISTRY._jobs.clear()
    yield
    BJ._REGISTRY._jobs.clear()


@pytest.fixture
def ws(tmp_path) -> Path:
    d = tmp_path / "ws"
    d.mkdir()
    return d


def _finished_record(ws: Path, job_id: str, **over) -> dict:
    rec = {
        "job_id": job_id,
        "pid": 999999999,
        "proc_start_ticks": None,
        "command": "orca big.inp",
        "description": "twelve-hour run",
        "cwd": str(ws),
        "workspace": str(ws),
        "stdout_path": str(ws / f"{job_id}.out"),
        "stderr_path": str(ws / f"{job_id}.err"),
        "started_at": time.time() - 120,
        "timeout_s": 3600,
        "deadline_at": time.time() + 3600,
        "cores": 1,
        "exit_code": 0,
        "finished_at": time.time() - 1,
        "acknowledged": False,
    }
    rec.update(over)
    BJ._persist_job_start(str(ws), rec)
    BJ._note_job_workspace(job_id, str(ws))
    return rec


def _ids(events) -> list[str]:
    return sorted(e["job_id"] for e in events)


# ---------------------------------------------------------------------------
# Two-phase acknowledgement
# ---------------------------------------------------------------------------

def test_a_drain_that_is_never_confirmed_delivers_the_event_again(ws, monkeypatch):
    """The turn died before the prompt was built. The notice must survive it."""
    _finished_record(ws, "longrun")

    assert _ids(BJ.drain_finished_events(ws)) == ["longrun"]
    # ...the turn dies here. Nothing confirmed. After the claim expires:
    monkeypatch.setattr(BJ, "_EVENT_CLAIM_GRACE_S", 0.0)

    assert _ids(BJ.drain_finished_events(ws)) == ["longrun"]


def test_a_confirmed_event_is_never_delivered_twice(ws, monkeypatch):
    _finished_record(ws, "longrun")

    events = BJ.drain_finished_events(ws)
    assert BJ.confirm_finished_events(events) == 1

    monkeypatch.setattr(BJ, "_EVENT_CLAIM_GRACE_S", 0.0)
    assert BJ.drain_finished_events(ws) == []


def test_a_confirmation_survives_a_restart(ws, monkeypatch):
    _finished_record(ws, "longrun")
    BJ.confirm_finished_events(BJ.drain_finished_events(ws))

    BJ._REGISTRY._jobs.clear()                      # the restart
    monkeypatch.setattr(BJ, "_EVENT_CLAIM_GRACE_S", 0.0)

    assert BJ.drain_finished_events(ws) == []


def test_a_claim_holds_off_a_concurrent_drain(ws):
    """Two front ends must not both announce the same completion."""
    _finished_record(ws, "longrun")

    assert _ids(BJ.drain_finished_events(ws)) == ["longrun"]
    assert BJ.drain_finished_events(ws) == []       # inside the grace


def test_the_confirmation_carries_the_workspace_the_job_belongs_to(
        ws, tmp_path):
    """A job belongs to the workspace it was STARTED in, not the current one."""
    far = tmp_path / "far"
    far.mkdir()
    _finished_record(far, "elsewhere")

    events = BJ.drain_all_finished_events(ws)

    assert [e["workspace"] for e in events] == [str(far)]
    assert BJ.confirm_finished_events(events) == 1


def test_confirming_nothing_is_harmless():
    assert BJ.confirm_finished_events([]) == 0
    assert BJ.confirm_finished_events([{"job_id": "x"}]) == 0


# ---------------------------------------------------------------------------
# The lock has to reach across processes
# ---------------------------------------------------------------------------

def _writer(index_path: str, workspace: str, first: int, count: int) -> None:
    from delfin.agent import bash_jobs as inner
    inner._INDEX_PATH = Path(index_path)
    for n in range(first, first + count):
        inner._note_job_workspace(f"job{n:04d}", workspace)


def test_two_processes_writing_the_index_do_not_lose_each_others_entries(
        ws, tmp_path):
    """A lost index entry is a job the agent can no longer name."""
    index = tmp_path / "index.json"
    ctx = multiprocessing.get_context("spawn")
    procs = [ctx.Process(target=_writer, args=(str(index), str(ws), base, 25))
             for base in (0, 1000)]
    for p in procs:
        p.start()
    for p in procs:
        p.join(timeout=60)

    entries = json.loads(index.read_text())["jobs"]

    assert len(entries) == 50, "each process's writes must all survive"


def test_the_registry_read_modify_write_is_locked_across_processes():
    """The lock must be an OS-level one; a threading lock orders one process."""
    source = Path(BJ.__file__).read_text(encoding="utf-8")
    body = source[source.index("def _locked("):]
    body = body[:body.index("\n@dataclass")]
    assert "flock" in body, (
        "two front ends share these files; only an OS lock orders them")


# ---------------------------------------------------------------------------
# A teardown must not delete a live job's tree
# ---------------------------------------------------------------------------

def _git(repo: Path, *args: str) -> None:
    subprocess.run(["git", *args], cwd=str(repo), check=True,
                   capture_output=True, text=True)


@pytest.fixture
def repo(tmp_path) -> Path:
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q")
    _git(r, "config", "user.email", "t@example.invalid")
    _git(r, "config", "user.name", "t")
    (r / "a.txt").write_text("hello\n")
    _git(r, "add", "a.txt")
    _git(r, "commit", "-qm", "initial")
    return r


def _hold(worktree: Path, pid: int) -> None:
    BJ._persist_job_start(str(worktree), {
        "job_id": "holder01",
        "pid": pid,
        "proc_start_ticks": BJ._proc_start_ticks(pid),
        "command": "orca big.inp",
        "description": "",
        "cwd": str(worktree),
        "workspace": str(worktree),
        "stdout_path": "", "stderr_path": "",
        "started_at": time.time(),
        "timeout_s": 3600,
        "deadline_at": time.time() + 3600,
        "cores": 1,
        "exit_code": None, "finished_at": None, "acknowledged": False,
    })
    BJ._note_job_workspace("holder01", str(worktree))


def test_exit_worktree_keeps_a_tree_that_still_runs_a_job(repo, tmp_path):
    """The tool path force-removed it: the child survives ``setsid``, so it
    kept running with its working directory deleted underneath it."""
    info = WT.enter_worktree(repo, parent=tmp_path / "wts")
    child = subprocess.Popen(["sleep", "30"], start_new_session=True)
    try:
        _hold(info.path, child.pid)

        WT.exit_worktree(info, keep_if_changed=True)

        assert info.held_by_jobs, "the guard must say what is holding the tree"
        assert info.cleaned_up is False
        assert info.path.is_dir(), "the live job's working directory survives"
        assert (info.path / ".delfin" / "bash_jobs.json").is_file()
    finally:
        child.kill()
        child.wait()


def test_exit_worktree_still_removes_a_tree_nobody_is_using(repo, tmp_path):
    info = WT.enter_worktree(repo, parent=tmp_path / "wts")

    WT.exit_worktree(info, keep_if_changed=True)

    assert info.held_by_jobs == []
    assert info.cleaned_up is True
    assert not info.path.exists()


def test_the_post_merge_cleanup_keeps_a_tree_that_still_runs_a_job(
        repo, tmp_path):
    """A merge does not stop the processes the tree is running."""
    info = WT.enter_worktree(repo, parent=tmp_path / "wts")
    (info.path / "b.txt").write_text("new\n")
    child = subprocess.Popen(["sleep", "30"], start_new_session=True)
    try:
        _hold(info.path, child.pid)

        result = WT.merge_worktree(info, cleanup=True)

        assert result.applied is True, "the merge itself must still happen"
        assert (repo / "b.txt").is_file()
        assert info.path.is_dir(), "but the live job keeps its directory"
        assert "still running" in result.message
        assert "holder01" in result.message
    finally:
        child.kill()
        child.wait()


def test_the_post_merge_cleanup_still_removes_an_idle_tree(repo, tmp_path):
    info = WT.enter_worktree(repo, parent=tmp_path / "wts")
    (info.path / "b.txt").write_text("new\n")

    result = WT.merge_worktree(info, cleanup=True)

    assert result.applied is True
    assert not info.path.exists()


def test_all_three_teardown_paths_ask_the_same_question():
    """One implementation, so they cannot drift apart on which of them checks."""
    from delfin.agent import subagents as SA

    assert SA._jobs_holding_worktree.__module__ == "delfin.agent.subagents"
    src = Path(SA.__file__).read_text(encoding="utf-8")
    body = src[src.index("def _jobs_holding_worktree"):]
    body = body[:body.index("def _kept_for_running_jobs")]
    assert "jobs_holding_worktree" in body and "from .worktree" in body

    wt_src = Path(WT.__file__).read_text(encoding="utf-8")
    for fn in ("def exit_worktree", "def _cleanup_merged_worktree"):
        body = wt_src[wt_src.index(fn):]
        body = body[:body.index("\n\ndef ")] if "\n\ndef " in body else body
        assert "jobs_holding_worktree" in body, f"{fn} does not check"


# ---------------------------------------------------------------------------
# The engine confirms only what actually reached the prompt
# ---------------------------------------------------------------------------

def _engine(tmp_path):
    from unittest.mock import MagicMock, patch

    from delfin.agent import engine as E
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                             model="kit.qwen3.5-397b-A17b", mode="solo")


def test_the_engine_confirms_the_events_it_put_in_the_prompt(
        tmp_path, monkeypatch):
    confirmed: list = []
    monkeypatch.setattr(
        "delfin.agent.bash_jobs.confirm_finished_events",
        lambda events: confirmed.extend(e["job_id"] for e in events))
    monkeypatch.setattr(
        "delfin.agent.bash_jobs.drain_all_finished_events",
        lambda ws: [{"job_id": "4711", "workspace": str(ws),
                     "command": "sbatch opt.slurm", "exit_code": 0,
                     "runtime_s": 3600.0, "stdout_tail": "done",
                     "stderr_tail": ""}])
    monkeypatch.setattr("delfin.agent.job_monitor.check_agent_jobs",
                        lambda ws: [])

    block = _engine(tmp_path)._build_finished_jobs_block()

    assert "4711" in block
    assert confirmed == ["4711"], (
        "the notice exists now, so the claim may become permanent")


def test_a_turn_that_never_built_the_block_confirms_nothing(
        tmp_path, monkeypatch):
    """Nothing was delivered, so nothing may be retired."""
    confirmed: list = []
    monkeypatch.setattr(
        "delfin.agent.bash_jobs.confirm_finished_events",
        lambda events: confirmed.extend(events))
    monkeypatch.setattr("delfin.agent.bash_jobs.drain_all_finished_events",
                        lambda ws: [])
    monkeypatch.setattr("delfin.agent.job_monitor.check_agent_jobs",
                        lambda ws: [])

    assert _engine(tmp_path)._build_finished_jobs_block() == ""
    assert confirmed == []


def test_a_wrapper_shell_that_left_children_says_so_in_the_prompt(
        tmp_path, monkeypatch):
    """'ok, 2s' for a twelve-hour cluster job is the failure this prevents."""
    monkeypatch.setattr(
        "delfin.agent.bash_jobs.drain_all_finished_events",
        lambda ws: [{"job_id": "4711", "workspace": str(ws),
                     "command": "sbatch opt.slurm", "exit_code": 0,
                     "runtime_s": 2.0, "stdout_tail": "", "stderr_tail": "",
                     "children_running": True,
                     "watched_slurm_jobs": ["99123"]}])
    monkeypatch.setattr("delfin.agent.job_monitor.check_agent_jobs",
                        lambda ws: [])

    block = _engine(tmp_path)._build_finished_jobs_block()

    assert "children still running" in block
    assert "99123" in block and "now watched" in block


def test_a_job_killed_at_its_cap_is_not_reported_as_a_clean_exit(
        tmp_path, monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.bash_jobs.drain_all_finished_events",
        lambda ws: [{"job_id": "4711", "workspace": str(ws),
                     "command": "orca big.inp", "exit_code": None,
                     "runtime_s": 86400.0, "stdout_tail": "",
                     "stderr_tail": "", "timed_out": True}])
    monkeypatch.setattr("delfin.agent.job_monitor.check_agent_jobs",
                        lambda ws: [])

    assert "wall-clock cap" in _engine(tmp_path)._build_finished_jobs_block()


def test_an_unreachable_queue_reaches_the_prompt_as_a_degradation(
        tmp_path, monkeypatch):
    monkeypatch.setattr("delfin.agent.bash_jobs.drain_all_finished_events",
                        lambda ws: [])
    monkeypatch.setattr(
        "delfin.agent.job_monitor.check_agent_jobs",
        lambda ws: [{"job_id": "555", "kind": "slurm", "state": "UNAVAILABLE",
                     "description": "twelve-hour opt", "ok": False,
                     "exit_code": None, "signatures": [],
                     "degraded": "the scheduler could not be asked about "
                                 "this job — its state is unknown, not "
                                 "unchanged"}])

    block = _engine(tmp_path)._build_finished_jobs_block()

    assert "unknown, not unchanged" in block
