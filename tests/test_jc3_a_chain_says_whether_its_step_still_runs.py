"""A chain says whether its open step is still being run.

Input: a chain whose step is ``submitted`` and whose process is alive, has
ended, ran on another machine (with or without a SLURM job), or was
written by an older DELFIN without a runner record.

Before: the state said ``submitted`` and nothing else. In JC2's live run 4
(2026-11-15) the restarted agent found the ORCA temporary files of the
killed optimisation and told the user the step was "actively running";
it never relayed the resume offer.
Now each submission records its runner (pid, start fingerprint, host,
PID namespace, heartbeat file). A pid is asked only where it names the
same process (same machine, same PID namespace, matching start): in live
run 5 (2026-11-15) the agent's sandbox called the runner "pid 3", and
from outside pid 3 is a kernel thread, so the killed step read "still
running". Elsewhere the heartbeat the running step touches decides; a
SLURM job is no sign (the agent's own job outlives its chains). The resume offer and ``delfin chain status`` say "still
running" / "not running, interrupted" / "cannot be checked from here",
the offer says leftover files are no sign of life, and ``delfin chain
resume`` refuses a chain whose step verifiably still runs.
"""
from __future__ import annotations

import os
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from delfin.agent import proc_identity
from delfin.tools import chain_liveness as cl
from delfin.tools import chain_state as cs


@pytest.fixture()
def home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    monkeypatch.setattr("pathlib.Path.home", lambda: home)
    monkeypatch.delenv("DELFIN_STATE", raising=False)
    import delfin.agent.session_store as ss
    monkeypatch.setattr(ss, "_SESSIONS_DIR", home / ".delfin" / "agent_sessions")
    return home


def _ended_runner() -> dict:
    p = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(0.5)"])
    start = proc_identity.process_start(p.pid)
    p.wait()
    return {"pid": p.pid, "start": start, "host": socket.gethostname(),
            "pid_ns": cl._pid_namespace()}


def _chain(wd: Path, runner) -> cs.ChainState:
    chain = cs.ChainState.create(name="demo", work_dir=wd)
    chain.add_step("xtb_opt", "xtb_opt", params={"charge": 0})
    chain.add_step("orca_opt", "orca_opt", params={"charge": 0})
    chain.steps[0].mark_submitted("local", runner=runner)
    chain.save()
    return chain


def test_this_process_is_running():
    step = cs.ChainStep(name="s", tool="t")
    step.mark_submitted("local", runner=cl.this_runner())
    assert cl.runner_liveness(step)[0] == cl.RUNNING


def test_an_ended_process_is_gone():
    step = cs.ChainStep(name="s", tool="t")
    step.mark_submitted("local", runner=_ended_runner())
    alive, why = cl.runner_liveness(step)
    assert alive == cl.GONE and "interrupted" in why


def test_an_old_state_without_a_runner_is_unknown():
    step = cs.ChainStep.from_dict({"name": "s", "tool": "t", "state": "submitted"})
    assert cl.runner_liveness(step)[0] == cl.UNKNOWN


def test_another_host_without_a_job_is_unknown():
    step = cs.ChainStep(name="s", tool="t")
    step.mark_submitted("local", runner={"pid": 1, "start": "x",
                                         "host": "some-other-node"})
    alive, why = cl.runner_liveness(step)
    assert alive == cl.UNKNOWN and "some-other-node" in why


def test_a_pid_from_another_pid_namespace_is_not_asked():
    """Live run 5: the sandbox's "pid 3" is a root kernel thread outside."""
    step = cs.ChainStep(name="s", tool="t")
    step.mark_submitted("local", runner={
        "pid": 1, "start": "221626440", "host": socket.gethostname(),
        "pid_ns": "pid:[1]"})
    assert cl.runner_liveness(step)[0] == cl.UNKNOWN


def test_another_users_pid_is_never_running():
    """kill(pid, 0) answers EPERM for a pid we cannot signal; that is not
    our runner, whatever ``proc_identity.alive`` makes of it."""
    step = cs.ChainStep(name="s", tool="t")
    step.mark_submitted("local", runner={
        "pid": 1, "start": "221626440", "host": socket.gethostname(),
        "pid_ns": cl._pid_namespace()})
    assert cl.runner_liveness(step)[0] != cl.RUNNING


def _with_heartbeat(tmp_path, age_s):
    beat = tmp_path / cl.HEARTBEAT_NAME
    beat.touch()
    t = time.time() - age_s
    os.utime(beat, (t, t))
    step = cs.ChainStep(name="s", tool="t")
    step.mark_submitted("local", runner={
        "pid": 3, "start": "x", "host": socket.gethostname(),
        "pid_ns": "pid:[1]", "heartbeat": str(beat)})
    return step


def test_a_fresh_heartbeat_is_running(tmp_path):
    assert cl.runner_liveness(_with_heartbeat(tmp_path, 5))[0] == cl.RUNNING


def test_a_stale_heartbeat_is_gone(tmp_path):
    alive, why = cl.runner_liveness(
        _with_heartbeat(tmp_path, cl.STALE_SECONDS + 60))
    assert alive == cl.GONE and "no sign of life" in why


def test_the_heartbeat_beats_and_stops(tmp_path, monkeypatch):
    monkeypatch.setattr(cl, "BEAT_SECONDS", 0.05)
    beat = tmp_path / "hb"
    with cl.beating(str(beat)):
        first = beat.stat().st_mtime_ns
        old = time.time() - 600
        os.utime(beat, (old, old))
        time.sleep(0.3)
        assert beat.stat().st_mtime_ns >= first
    assert not [t for t in threading.enumerate() if t.name == "chain-heartbeat"]


def test_finished_and_pending_steps_are_not_asked():
    step = cs.ChainStep(name="s", tool="t")
    assert cl.runner_liveness(step) == ("", "")


def test_the_runner_is_persisted_and_cleared_on_retry(tmp_path):
    chain = _chain(tmp_path / "c", cl.this_runner())
    again = cs.ChainState.load(chain.state_path)
    assert again.steps[0].runner["pid"] == os.getpid()
    step = again.steps[0]
    step.mark_failed("killed")
    step.reset_for_retry()
    assert step.runner is None


def test_a_pipeline_run_records_its_runner(home, monkeypatch):
    def fake_run_step(step_name, geometry=None, work_dir=None, cores=1, **kw):
        from delfin.tools.pipeline import StepResult, StepStatus
        out = Path(work_dir) / f"{step_name}.out"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("fake run output\n")
        return StepResult(step_name=step_name, status=StepStatus.SUCCESS,
                          work_dir=Path(work_dir))

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)
    from delfin.tools.pipeline import Pipeline

    wd = home / "runs" / "a"
    geom = wd.parent / "start.xyz"
    geom.parent.mkdir(parents=True, exist_ok=True)
    geom.write_text("2\n\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", charge=0)
    pipe.run(work_dir=wd, geometry=geom, cores=1)
    step = cs.ChainState.load(wd / "chain_state.json").steps[0]
    assert step.runner and step.runner["pid"] == os.getpid()
    assert step.runner["host"] == socket.gethostname()
    assert Path(step.runner["heartbeat"]).parent == wd
    assert Path(step.runner["heartbeat"]).exists()


def test_the_offer_says_an_interrupted_step_is_not_running(home):
    from delfin.agent import session_store as ss

    ws = home / "ws"
    _chain(ws / "chain1", _ended_runner())
    text = ss._render_chain_offer("sess", workspace=ws)
    assert "do not mean" in text and "temporary files" in text
    assert "not running" in text and "interrupted" in text
    assert "resume with: delfin chain resume" in text


def test_the_offer_does_not_offer_a_resume_of_a_running_chain(home):
    from delfin.agent import session_store as ss

    ws = home / "ws"
    _chain(ws / "chain1", cl.this_runner())
    text = ss._render_chain_offer("sess", workspace=ws)
    assert "still running" in text and "do NOT resume" in text
    assert "resume with:" not in text


def test_status_names_the_process(tmp_path, capsys):
    from delfin.cli_chain import main

    _chain(tmp_path / "c", _ended_runner())
    assert main(["status", str(tmp_path / "c")]) == 0
    assert "process: not running" in capsys.readouterr().out


def test_resume_refuses_a_chain_that_still_runs(tmp_path, monkeypatch, capsys):
    from delfin.cli_chain import main

    chain = _chain(tmp_path / "c", cl.this_runner())
    before = chain.state_path.read_bytes()
    ran = []
    monkeypatch.setattr("delfin.tools.pipeline.Pipeline.run",
                        lambda self, **kw: ran.append(kw))
    assert main(["resume", str(tmp_path / "c")]) == 1
    assert not ran
    assert chain.state_path.read_bytes() == before
    assert "NOT resumed" in capsys.readouterr().out


def test_a_run_that_raised_in_a_living_process_is_gone(home, monkeypatch):
    """The process lives on (a service, the tool server, a kernel) after
    the step's run raised: its pid answers, yet the step is not running —
    otherwise the chain would read "still running" for ever and every
    resume would be refused."""
    def raising_run_step(step_name, **kw):
        raise RuntimeError("simulated failure inside the step")

    monkeypatch.setattr("delfin.tools._runner.run_step", raising_run_step)
    from delfin.tools.pipeline import Pipeline

    wd = home / "runs" / "b"
    geom = wd.parent / "start.xyz"
    geom.parent.mkdir(parents=True, exist_ok=True)
    geom.write_text("2\n\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", charge=0)
    with pytest.raises(RuntimeError):
        pipe.run(work_dir=wd, geometry=geom, cores=1)
    step = cs.ChainState.load(wd / "chain_state.json").steps[0]
    assert step.state == "submitted"          # JC1: a crash leaves it so
    alive, why = cl.runner_liveness(step)
    assert alive == cl.GONE and "RuntimeError" in why
