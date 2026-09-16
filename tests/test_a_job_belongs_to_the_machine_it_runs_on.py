"""A background job's record is judged on the machine that runs it.

Login nodes share the home directory and with it every job registry. A pid
from node A names nothing on node B, or a stranger's process with the same
number (audit 2026-09-16). Node B used to report A's running job as
finished, and past its deadline it would signal whatever local process
carried that pid. After a reboot of the same machine, the pid is equally
a stranger.
"""
import json
import subprocess
import time

import pytest

from delfin.agent import bash_jobs as BJ


@pytest.fixture
def stranger():
    """A live local process that happens to carry the recorded pid."""
    proc = subprocess.Popen(["sleep", "60"], start_new_session=True)
    yield proc
    proc.kill()
    proc.wait()


def _rec(pid, **over):
    rec = {"job_id": "bg-remote", "pid": pid, "proc_start_ticks": None,
           "command": "orca big.inp", "cwd": "/x", "workspace": "/x",
           "stdout_path": "", "stderr_path": "", "started_at": time.time() - 60,
           "timeout_s": 3600, "finished_at": None, "exit_code": None,
           "acknowledged": False, "host": "another-login-node", "boot_id": "b-1"}
    rec.update(over)
    return rec


def _write(ws, rec):
    (ws / ".delfin").mkdir(parents=True, exist_ok=True)
    (ws / ".delfin" / "bash_jobs.json").write_text(json.dumps({"jobs": {rec["job_id"]: rec}}))


def test_a_new_record_names_its_machine():
    here = BJ.this_machine()
    assert here["host"]
    rec = _rec(1, **here)
    assert not BJ.record_runs_elsewhere(rec)


def test_another_nodes_running_job_is_not_announced_finished(tmp_path):
    _write(tmp_path, _rec(999999999))
    assert BJ.drain_finished_events(str(tmp_path)) == []
    data = json.loads((tmp_path / ".delfin" / "bash_jobs.json").read_text())
    assert data["jobs"]["bg-remote"]["finished_at"] is None


def test_a_pid_from_another_node_is_never_signalled(tmp_path, stranger):
    rec = _rec(stranger.pid, started_at=time.time() - 3 * 86400, timeout_s=60)
    assert BJ._enforce_deadline(rec, time.time()) is False
    _write(tmp_path, rec)
    BJ.drain_finished_events(str(tmp_path))
    job = BJ._reattach("bg-remote", str(tmp_path))
    assert job is not None and job.poll() is not None     # past its cap
    time.sleep(0.2)
    assert stranger.poll() is None, "a stranger's process was signalled"


def test_ending_another_nodes_job_from_here_is_refused(tmp_path, stranger, monkeypatch):
    _write(tmp_path, _rec(stranger.pid))
    reg = BJ.get_registry()
    monkeypatch.setattr(reg, "get", lambda jid, ws=None: BJ._reattach(jid, str(tmp_path)))
    ok, message = reg.kill("bg-remote")
    assert ok is False and "another-login-node" in message
    time.sleep(0.2)
    assert stranger.poll() is None


def test_another_nodes_job_takes_no_core_here_and_keeps_its_folder(stranger):
    rec = _rec(stranger.pid)
    assert BJ._record_holds_the_node(rec, None, time.time()) is False
    assert BJ._record_alive(rec) is True          # its worktree is not torn down


def test_a_job_from_before_a_reboot_is_gone(stranger):
    rec = _rec(stranger.pid, host=BJ.this_machine()["host"], boot_id="an-earlier-boot")
    if not BJ.this_machine()["boot_id"]:
        pytest.skip("no boot id on this host")
    assert BJ._record_alive(rec) is False
    assert BJ._record_holds_the_node(rec, None, time.time()) is False
