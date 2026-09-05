"""The bookkeeping died with the process; the job did not.

Every background command is started with ``os.setsid``, so the child
deliberately survives the agent process. Three consequences of that were
never followed through, and each one costs a shared node directly.

THE WALL-CLOCK CAP DIED WITH ITS THREAD. The 24 h cap was enforced by one
daemon thread inside the starting process. Re-attach read ``timeout_s``
into the recovered view and never looked at it again, so after any restart
every surviving job had NO wall-clock bound at all.

THE PRUNE THEN DELETED THE LIVE JOB. ``_prune_old_records`` dropped
records older than ~7 days with no aliveness check, on the premise its own
docstring stated -- "the 24 h hard timeout guarantees any such job is long
over" -- which the point above falsifies. The result was a live 32-core
process that nobody can name: ``bash_status`` and ``bash_kill`` both answer
"unknown job_id", which reads as a typo rather than as a calculation
holding a node.

THE CAP COUNTED THE WRONG THING, IN THE WRONG PLACE. It counted an
in-memory list, so a restart reset it to zero while N jobs ran and two
front ends gave two registries and twice the cap. And it counted
PROCESSES: eight ``orca -np 64`` passed a cap of eight, which is 512 cores
of somebody else's queue, from a module whose stated premise is that this
is other people's CPU.

Plus the completion signal itself: a background job's "finished" was the
wrapper shell's exit, so ``bash_background("sbatch run.sh")`` reported
"ok, 2s" for a twelve-hour cluster job.

No SLURM here -- ``sbatch`` output is written by hand and the registrar is
injected. Real subprocesses are short-lived sleeps.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import time
from pathlib import Path

import pytest

from delfin.agent import bash_jobs as BJ


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


def _record(ws: Path, job_id: str, **over) -> dict:
    rec = {
        "job_id": job_id,
        "pid": 999999999,                 # definitely dead
        "proc_start_ticks": None,
        "command": "orca big.inp -np 64",
        "description": "",
        "cwd": str(ws),
        "workspace": str(ws),
        "stdout_path": str(ws / f"{job_id}.out"),
        "stderr_path": str(ws / f"{job_id}.err"),
        "started_at": time.time() - 60,
        "timeout_s": 3600,
        "cores": 1,
        "exit_code": None,
        "finished_at": None,
        "acknowledged": False,
    }
    rec.update(over)
    return rec


def _register(ws: Path, rec: dict) -> None:
    BJ._persist_job_start(str(ws), rec)
    BJ._note_job_workspace(rec["job_id"], str(ws))


def _read(ws: Path, job_id: str) -> dict:
    data = json.loads((ws / ".delfin" / "bash_jobs.json").read_text())
    return data["jobs"][job_id]


# ---------------------------------------------------------------------------
# The prune must not delete a running job's record
# ---------------------------------------------------------------------------

def test_a_live_job_keeps_its_record_past_the_prune_horizon(ws):
    """The record IS the only way to name, watch or kill that process."""
    child = subprocess.Popen(["sleep", "30"], start_new_session=True)
    try:
        _register(ws, _record(
            ws, "liveold1", pid=child.pid,
            proc_start_ticks=BJ._proc_start_ticks(child.pid),
            started_at=time.time() - 9 * 24 * 3600,
            timeout_s=BJ._DEFAULT_BG_TIMEOUT_S,
            deadline_at=time.time() + 3600))

        jobs = BJ._load_registry_file(ws)["jobs"]
        pruned = BJ._prune_old_records(jobs, time.time())

        assert pruned is False
        assert "liveold1" in jobs
        assert BJ.get_registry().get("liveold1", ws) is not None
    finally:
        child.kill()
        child.wait()


def test_a_dead_job_older_than_the_horizon_is_still_pruned(ws):
    """The prune still does its job — it just checks first."""
    _register(ws, _record(ws, "deadold1",
                          started_at=time.time() - 9 * 24 * 3600))

    jobs = BJ._load_registry_file(ws)["jobs"]

    assert BJ._prune_old_records(jobs, time.time()) is True
    assert jobs == {}


# ---------------------------------------------------------------------------
# The wall-clock cap is re-armed after a restart
# ---------------------------------------------------------------------------

def test_the_deadline_is_absolute_so_it_survives_the_process(ws):
    job = BJ.get_registry().start("sleep 0.1", cwd=str(ws), workspace=ws,
                                  timeout_s=120)
    rec = _read(ws, job.job_id)

    assert rec["deadline_at"] == pytest.approx(rec["started_at"] + 120, abs=2)


def test_a_reattached_job_past_its_deadline_is_ended_not_left_running(ws):
    """After a restart nothing else will ever stop it."""
    child = subprocess.Popen(["sleep", "60"], start_new_session=True)
    try:
        _register(ws, _record(
            ws, "overdue1", pid=child.pid,
            proc_start_ticks=BJ._proc_start_ticks(child.pid),
            started_at=time.time() - 500, timeout_s=100,
            deadline_at=time.time() - 400))

        job = BJ.get_registry().get("overdue1", ws)
        assert job.poll() is not None, "an overdue job must not report running"

        deadline = time.monotonic() + 5
        while child.poll() is None and time.monotonic() < deadline:
            time.sleep(0.05)
        assert child.poll() is not None, "the process itself must be stopped"
        assert _read(ws, "overdue1")["timed_out"] is True
    finally:
        if child.poll() is None:
            child.kill()
        child.wait()


def test_a_reattached_job_inside_its_deadline_still_reports_running(ws):
    child = subprocess.Popen(["sleep", "30"], start_new_session=True)
    try:
        _register(ws, _record(
            ws, "intime1", pid=child.pid,
            proc_start_ticks=BJ._proc_start_ticks(child.pid),
            started_at=time.time() - 10, timeout_s=3600,
            deadline_at=time.time() + 3590))

        assert BJ.get_registry().get("intime1", ws).poll() is None
        assert child.poll() is None
    finally:
        child.kill()
        child.wait()


def test_a_record_without_a_deadline_falls_back_to_the_hard_cap(ws):
    """Records written before deadlines existed must still get a bound."""
    rec = _record(ws, "legacy1", started_at=1000.0, timeout_s=7200)
    rec.pop("deadline_at", None)

    assert BJ._record_deadline(rec) == 1000.0 + 7200


def test_the_fallback_deadline_can_never_exceed_the_hard_cap(ws):
    rec = _record(ws, "legacy2", started_at=1000.0,
                  timeout_s=10 * BJ._DEFAULT_BG_TIMEOUT_S)

    assert BJ._record_deadline(rec) == 1000.0 + BJ._DEFAULT_BG_TIMEOUT_S


# ---------------------------------------------------------------------------
# The cap counts live pids across registries, and it counts cores
# ---------------------------------------------------------------------------

def test_a_restart_does_not_reset_the_cap_to_zero(ws):
    """The in-memory list is exactly what a restart loses; the job is not."""
    child = subprocess.Popen(["sleep", "30"], start_new_session=True)
    try:
        _register(ws, _record(
            ws, "survivor", pid=child.pid,
            proc_start_ticks=BJ._proc_start_ticks(child.pid),
            deadline_at=time.time() + 3600))
        BJ._REGISTRY._jobs.clear()              # the restart

        assert BJ.get_registry().count_running() == 1
    finally:
        child.kill()
        child.wait()


def test_a_second_front_end_sees_the_first_ones_jobs(ws, tmp_path):
    """Two dashboards used to mean two registries and twice the cap."""
    other = tmp_path / "other_ws"
    other.mkdir()
    kids = [subprocess.Popen(["sleep", "30"], start_new_session=True)
            for _ in range(2)]
    try:
        for i, (folder, kid) in enumerate(zip((ws, other), kids)):
            _register(folder, _record(
                folder, f"front{i}", pid=kid.pid,
                proc_start_ticks=BJ._proc_start_ticks(kid.pid),
                deadline_at=time.time() + 3600))

        assert BJ.get_registry().count_running() == 2
    finally:
        for kid in kids:
            kid.kill()
            kid.wait()


@pytest.mark.parametrize("command,cores", [
    ("orca big.inp -np 64", 64),
    ("mpirun -np 8 xtb", 8),
    ("xtb --threads 12 mol.xyz", 12),
    ("crest -T 24 in.xyz", 24),
    ("sbatch --cpus-per-task=32 run.sh", 32),
    ("OMP_NUM_THREADS=16 python run.py", 16),
    ("orca job.inp PAL8", 8),
    ("pytest -q", 1),
    ("head -c 4096 big.log", 1),
    ("bash -c 'echo hi'", 1),
])
def test_the_requested_width_is_read_from_the_command(command, cores):
    assert BJ._requested_cores(command) == cores


def test_eight_wide_jobs_do_not_pass_a_cap_of_eight_processes(ws, monkeypatch):
    """512 cores of somebody else's queue used to be eight jobs, i.e. fine."""
    monkeypatch.setattr(BJ, "_core_budget", lambda: 64)
    monkeypatch.setattr(BJ, "_max_background_jobs", lambda: 8)
    child = subprocess.Popen(["sleep", "30"], start_new_session=True)
    try:
        _register(ws, _record(
            ws, "wide1", pid=child.pid, cores=64,
            proc_start_ticks=BJ._proc_start_ticks(child.pid),
            deadline_at=time.time() + 3600))

        with pytest.raises(ValueError) as excinfo:
            BJ.get_registry().start("orca second.inp -np 64", cwd=str(ws),
                                    workspace=ws)

        message = str(excinfo.value)
        assert "64 core(s)" in message
        assert "budget is 64" in message
        assert "bash_status" in message, (
            "the refusal must name the action that clears it")
    finally:
        child.kill()
        child.wait()


def test_a_narrow_job_still_starts_beside_a_wide_one(ws, monkeypatch):
    """The cap must not become a reason to stop working."""
    monkeypatch.setattr(BJ, "_core_budget", lambda: 64)
    monkeypatch.setattr(BJ, "_max_background_jobs", lambda: 8)
    child = subprocess.Popen(["sleep", "30"], start_new_session=True)
    try:
        _register(ws, _record(
            ws, "wide2", pid=child.pid, cores=60,
            proc_start_ticks=BJ._proc_start_ticks(child.pid),
            deadline_at=time.time() + 3600))

        job = BJ.get_registry().start("echo ok", cwd=str(ws), workspace=ws)
        assert job.job_id
    finally:
        child.kill()
        child.wait()


def test_one_command_wider_than_the_whole_grant_is_refused_outright(
        ws, monkeypatch):
    monkeypatch.setattr(BJ, "_core_budget", lambda: 8)

    with pytest.raises(ValueError) as excinfo:
        BJ.get_registry().start("orca huge.inp -np 512", cwd=str(ws),
                                workspace=ws)

    assert "was granted 8" in str(excinfo.value)


def test_the_started_job_records_the_width_it_asked_for(ws, monkeypatch):
    monkeypatch.setattr(BJ, "_core_budget", lambda: 64)
    job = BJ.get_registry().start("xtb --threads 4 m.xyz", cwd=str(ws),
                                  workspace=ws)

    assert _read(ws, job.job_id)["cores"] == 4


# ---------------------------------------------------------------------------
# "Finished" must mean the work, not the wrapper shell
# ---------------------------------------------------------------------------

def test_a_shell_that_left_children_behind_is_not_reported_as_done(ws):
    """The wrapper exits; what it started keeps the node."""
    job = BJ.get_registry().start(
        "sleep 20 & echo started", cwd=str(ws), workspace=ws)
    deadline = time.monotonic() + 10
    while job.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    assert job.poll() is not None, "the wrapper shell should have exited"

    status = job.status_dict()

    try:
        assert status["children_running"] is True
        assert "children still running" in status["note"]
    finally:
        try:
            os.killpg(os.getpgid(job.proc.pid), signal.SIGKILL)
        except Exception:
            pass


def test_a_shell_that_left_children_behind_keeps_its_slot(ws, monkeypatch):
    """Releasing the slot is the other half of the harm: the cores are still
    spent, so the cap must still count them."""
    monkeypatch.setattr(BJ, "_core_budget", lambda: 64)
    job = BJ.get_registry().start(
        "sleep 20 & echo started", cwd=str(ws), workspace=ws)
    deadline = time.monotonic() + 10
    while job.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    assert job.poll() is not None
    BJ._REGISTRY._jobs.clear()          # only the persisted record remains

    try:
        assert [j["job_id"] for j in BJ.live_jobs()] == [job.job_id]
        assert BJ.get_registry().count_running() == 1
    finally:
        try:
            os.killpg(os.getpgid(job.proc.pid), signal.SIGKILL)
        except Exception:
            pass


def test_a_worktree_running_only_orphaned_children_is_still_held(ws, monkeypatch):
    job = BJ.get_registry().start(
        "sleep 20 & echo started", cwd=str(ws), workspace=ws)
    deadline = time.monotonic() + 10
    while job.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)

    try:
        held = BJ.running_jobs_for_workspace(ws)
        assert [h["job_id"] for h in held] == [job.job_id]
    finally:
        try:
            os.killpg(os.getpgid(job.proc.pid), signal.SIGKILL)
        except Exception:
            pass


def test_an_empty_group_releases_the_slot(ws, monkeypatch):
    """The cap must not leak slots — that stops the agent working."""
    monkeypatch.setattr(BJ, "_core_budget", lambda: 64)
    job = BJ.get_registry().start("echo done", cwd=str(ws), workspace=ws)
    deadline = time.monotonic() + 10
    while job.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    BJ._REGISTRY._jobs.clear()

    assert BJ.live_jobs() == []
    assert BJ.running_jobs_for_workspace(ws) == []


def test_a_job_that_really_finished_says_nothing_about_children(ws):
    job = BJ.get_registry().start("echo done", cwd=str(ws), workspace=ws)
    deadline = time.monotonic() + 10
    while job.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)

    status = job.status_dict()

    assert status["running"] is False
    assert status.get("children_running") is not True
    assert "note" not in status


def test_a_submitted_cluster_job_is_watched_even_from_the_background_path(
        ws, monkeypatch):
    """The auto-watch was wired into the FOREGROUND bash path only, so the
    one submission nobody reads was the one nobody watched."""
    registered: list[tuple] = []
    monkeypatch.setattr(
        "delfin.agent.job_monitor.register_agent_job",
        lambda workspace, jid, description="": registered.append(
            (str(workspace), jid)))

    rec = _record(ws, "submitter", finished_at=time.time(), exit_code=0)
    Path(rec["stdout_path"]).write_text(
        "preparing\nSubmitted batch job 991155\n")
    _register(ws, rec)

    events = BJ.drain_finished_events(ws)

    assert registered == [(str(ws), "991155")]
    assert events[0]["watched_slurm_jobs"] == ["991155"]


def test_the_registered_watch_is_not_repeated_on_every_drain(ws, monkeypatch):
    calls: list = []
    monkeypatch.setattr(
        "delfin.agent.job_monitor.register_agent_job",
        lambda workspace, jid, description="": calls.append(jid))
    rec = _record(ws, "submitter2", finished_at=time.time(), exit_code=0)
    Path(rec["stdout_path"]).write_text("Submitted batch job 42\n")
    _register(ws, rec)

    BJ.drain_finished_events(ws)
    BJ.confirm_finished_events([{"workspace": str(ws), "job_id": "submitter2"}])
    BJ.drain_finished_events(ws)

    assert calls == ["42"]
