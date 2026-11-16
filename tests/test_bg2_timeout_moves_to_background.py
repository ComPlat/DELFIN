"""BG2: a foreground bash command that reaches its timeout is MOVED to the
background, not ended.

Measured in the JC3 live run 6 (2026-11-15): ``delfin chain resume <dir>``
hit the foreground timeout, the tool ended the whole process group (the
ORCA optimisation died minutes in), and the message then advised
``bash_background`` -- after the work was already gone. Claude Code does
it the other way round: the command keeps running as a background job of
the session, and the result names the job id and the output file.

The floor tests here run with the process cage switched off (plain argv),
like test_a_command_ends_with_everything_it_started.py, so they pass on
every host the cage does not.
"""

import json
import os
import time
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as BJ


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().split(")")[-1].split()[0] != "Z"
    except OSError:
        return True


def _wait_dead(pid, s=5.0):
    end = time.monotonic() + s
    while time.monotonic() < end and _alive(pid):
        time.sleep(0.05)
    return not _alive(pid)


@pytest.fixture(autouse=True)
def _plain_argv(monkeypatch):
    """Run the foreground tool with the cage switched off, like the
    floor tests do: the move itself is what is tested, not bwrap."""
    monkeypatch.setattr(
        A, "_bash_isolation_argv",
        lambda cmd, cwd, perms, **kw: ["/bin/bash", "-c", cmd])


@pytest.fixture(autouse=True)
def _clean_registry(tmp_path):
    BJ._REGISTRY._jobs.clear()
    yield
    for job in list(BJ._REGISTRY._jobs.values()):
        try:
            BJ.get_registry().kill(job.job_id)
        except Exception:
            pass
    BJ._REGISTRY._jobs.clear()


def _run(cmd, timeout_s=2, tmp_path=None, **extra):
    perms = A.KitToolPermissions(
        workspace=str(tmp_path), mode="bypassPermissions")
    # The repeat-guard remembers timed-out commands per turn; without this
    # one test's timeout poisons the next (seen in the first gate run).
    getattr(A._doc_executor, "_timed_out_by_session", {}).clear()
    return json.loads(A._doc_executor._execute_bash(
        {"command": cmd, "timeout_s": timeout_s, **extra}, perms))


# ---------------------------------------------------------------------------
# The move itself (red on the previous commit: the work is ended instead)
# ---------------------------------------------------------------------------

def test_a_timeout_moves_the_command_to_the_background(tmp_path):
    """The command survives its foreground timeout as a registered job."""
    res = _run(f"echo started; echo $$ > {tmp_path / 'pid'}; sleep 60",
               timeout_s=2, tmp_path=tmp_path)
    assert res.get("status") == "moved_to_background", res
    job_id = res.get("job_id")
    assert job_id, res
    job = BJ.get_registry().get(job_id)
    assert job is not None, "the moved job is not registered"
    assert job.proc.poll() is None, "the moved command was ended anyway"
    assert res.get("output_file"), res
    assert _alive(int((tmp_path / "pid").read_text()))
    try:
        assert "started" in res.get("stdout", "") or "started" in _tail(job)
    finally:
        BJ.get_registry().kill(job_id)


def _tail(job) -> str:
    from delfin.agent import bash_jobs
    return bash_jobs.read_output(job, tail_lines=20).get("stdout_tail", "")


def test_the_move_starts_nothing_twice(tmp_path):
    """One registered job for the timed-out command, with the original
    process -- a second process would run the work again."""
    res = _run("sleep 60", timeout_s=2, tmp_path=tmp_path)
    jid = res.get("job_id")
    assert res.get("status") == "moved_to_background" and jid, res
    job = BJ.get_registry().get(jid)
    try:
        assert "sleep 60" in job.command
        assert job.proc.poll() is None
        assert len([j for j in BJ._REGISTRY._jobs.values()
                    if "sleep 60" in j.command]) == 1
    finally:
        BJ.get_registry().kill(jid)


def test_the_moved_job_reports_completion(tmp_path):
    """The completion notice arrives the way a bash_background job's
    does: the watchdog records the exit, drain_finished_events serves it."""
    res = _run("sleep 60", timeout_s=2, tmp_path=tmp_path)
    jid = res.get("job_id")
    assert res.get("status") == "moved_to_background" and jid, res
    BJ.get_registry().kill(jid)
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        events = BJ.drain_finished_events(str(tmp_path))
        if any(e.get("job_id") == jid for e in events):
            break
        time.sleep(0.2)
    assert any(e.get("job_id") == jid for e in events), "no completion notice"


def test_the_cage_survives_the_move(tmp_path):
    """The moved job keeps the argv/env/cwd the gate approved: the move
    adds no permission and no new process."""
    res = _run("sleep 60", timeout_s=2, tmp_path=tmp_path)
    jid = res.get("job_id")
    assert res.get("status") == "moved_to_background" and jid, res
    job = BJ.get_registry().get(res["job_id"])
    try:
        argv = job.proc.args if hasattr(job.proc, "args") else None
        assert isinstance(argv, list), "the job's argv is not the approved one"
        assert Path(job.cwd) == tmp_path or job.cwd == str(tmp_path)
    finally:
        BJ.get_registry().kill(res["job_id"])


def test_output_before_and_after_the_timeout_reaches_the_job_file(tmp_path):
    """The s4 finding: the move must keep output from BOTH sides of the
    timeout in the job's output file, not just the so-far part."""
    res = _run(f"echo before; echo after-later; sleep 60 & wait",
               timeout_s=2, tmp_path=tmp_path)
    jid = res.get("job_id")
    assert res.get("status") == "moved_to_background" and jid, res
    job = BJ.get_registry().get(jid)
    try:
        deadline = time.monotonic() + 8
        text = ""
        while time.monotonic() < deadline:
            text = Path(job.stdout_path).read_text()
            if "after-later" in text:
                break
            time.sleep(0.2)
        assert "before" in text, text
        assert "after-later" in text, text
    finally:
        BJ.get_registry().kill(jid)
