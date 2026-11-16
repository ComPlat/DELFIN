"""BG2 phase 2: adopt() registers a live process like a started job.

The demote path's own public API: a foreground command's live process
joins the registry exactly like a bash_background job -- same record
fields, same deadline cap, same completion event -- with the capture
attached so output from before and after the adoption both land in the
job's output files. Unlike start(), adoption never refuses on the job or
core caps (the process already runs; refusing to name it cannot undo it).
"""

import json
import subprocess
import time
from pathlib import Path

import pytest

from delfin.agent import bash_jobs as BJ


class _FakeCapture:
    """A stand-in for contained_run's pipe readers: it records what it was
    asked to write so a test can pin the order, and answers join()."""

    def __init__(self):
        self.attached: list = []
        self.joined = False

    def attach(self, sout, serr):
        self.attached.append((sout, serr))
        sout.write("before\n")
        serr.write("err-before\n")
        sout.flush()

    def join(self):
        self.joined = True


@pytest.fixture(autouse=True)
def _isolated(tmp_path):
    BJ._INDEX_PATH = tmp_path / "index.json"
    BJ._REGISTRY._jobs.clear()
    yield
    for job in list(BJ._REGISTRY._jobs.values()):
        try:
            BJ.get_registry().kill(job.job_id)
        except Exception:
            pass
    BJ._REGISTRY._jobs.clear()


def _record(ws: Path, job_id: str) -> dict:
    return json.loads(BJ._registry_path(ws).read_text())["jobs"][job_id]


def test_adopt_registers_a_live_process_like_a_started_job(tmp_path):
    child = subprocess.Popen(
        ["/bin/bash", "-c", "echo started; sleep 60"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, text=True)
    try:
        job = BJ.get_registry().adopt(
            command="sleep 60", cwd=str(tmp_path), proc=child,
            workspace=str(tmp_path), timeout_s=120,
            foreground_elapsed_s=2.0)
        rec = _record(tmp_path, job.job_id)
        assert rec["origin"] == "adopted"
        assert rec["foreground_s"] == 2.0
        assert rec["pid"] == child.pid
        assert rec["exit_code"] is None and rec["finished_at"] is None
        assert rec["deadline_at"] == pytest.approx(
            rec["started_at"] + 120, abs=5)
        assert job.proc is child, "adopt re-forked instead of adopting"
        assert job.proc.poll() is None
        assert job.elapsed_s() >= 2.0, (
            "runtime does not count the foreground part")
    finally:
        child.kill()
        child.wait()


def test_adopt_refuses_a_process_that_already_exited(tmp_path):
    dead = subprocess.Popen(["true"])
    dead.wait()
    with pytest.raises(ValueError):
        BJ.get_registry().adopt(
            command="true", cwd=str(tmp_path), proc=dead,
            workspace=str(tmp_path))


def test_adopt_attaches_the_capture_so_early_output_is_kept(tmp_path):
    child = subprocess.Popen(
        ["/bin/bash", "-c", "sleep 60"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, text=True)
    try:
        capture = _FakeCapture()
        job = BJ.get_registry().adopt(
            command="sleep 60", cwd=str(tmp_path), proc=child,
            workspace=str(tmp_path), timeout_s=120, capture=capture)
        assert capture.attached, "the capture was never attached"
        head = Path(job.stdout_path).read_text()
        serr_head = Path(job.stderr_path).read_text()
        assert "before" in head, head
        assert "err-before" in serr_head, serr_head
    finally:
        child.kill()
        child.wait()


def test_adopt_records_where_the_environment_came_from(tmp_path):
    """F1(a): the record names the provenance of the running process's
    environment (the registry cannot read a process's environ)."""
    child = subprocess.Popen(
        ["/bin/bash", "-c", "sleep 60"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, text=True)
    try:
        job = BJ.get_registry().adopt(
            command="sleep 60", cwd=str(tmp_path), proc=child,
            workspace=str(tmp_path), timeout_s=120, env_source="gate")
        rec = _record(tmp_path, job.job_id)
        assert rec["env_source"] == "gate"
        assert rec["cores"] == 1  # started-record parity
    finally:
        child.kill()
        child.wait()


def test_adopt_records_the_affinity_the_process_runs_under(tmp_path):
    """F1(b): an adopted job's cpu_set is NOT applied post-hoc (the
    process already runs; re-binding the leader would not reach its
    children) -- the record describes the affinity the kernel already
    holds, so it stays honest."""
    child = subprocess.Popen(
        ["/bin/bash", "-c", "sleep 60"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, text=True)
    try:
        job = BJ.get_registry().adopt(
            command="sleep 60", cwd=str(tmp_path), proc=child,
            workspace=str(tmp_path), timeout_s=120, cpu_set=[0, 1])
        rec = _record(tmp_path, job.job_id)
        assert rec["cores"] == 1
        assert job.proc.poll() is None
    finally:
        child.kill()
        child.wait()


def test_adopt_is_not_refused_when_the_caps_are_full(tmp_path, monkeypatch):
    """Refusing to register a running process leaves work nobody can
    name, watch or stop -- worse than an over-cap registry."""
    child = subprocess.Popen(
        ["/bin/bash", "-c", "sleep 60"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, text=True)
    try:
        monkeypatch.setattr(BJ, "_MAX_BG_JOBS_CEILING", 0)
        monkeypatch.setattr(BJ, "_available_cpus", lambda: 0)
        job = BJ.get_registry().adopt(
            command="sleep 60", cwd=str(tmp_path), proc=child,
            workspace=str(tmp_path), timeout_s=120)
        assert job.proc.poll() is None
        assert _record(tmp_path, job.job_id)["job_id"] == job.job_id
    finally:
        child.kill()
        child.wait()
