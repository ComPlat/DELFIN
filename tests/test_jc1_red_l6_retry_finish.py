"""JC1 live-run-#3 finding L6: a retried step's success/failure cannot be
recorded because the retry path never marks the step submitted/running
after reset_for_retry — mark_finished/mark_failed are then refused from
'pending' and an 11-minute converged run is thrown away (live evidence).

Also: next_attempt_dir numbers the SECOND attempt 'attempt_1' because the
first attempt's numbered dir (01_<name>) is not counted.

Controls (red on d537996a1):
- a retried step that SUCCEEDS must end finished with evidence;
- a retried step that FAILS must end failed (not the mark_failed refused
  crash), the failure recorded, attempt 2 the last ledger entry;
- attempt numbering: the new dir for a step whose 1st attempt was
  01_<name> must be attempt_2.
"""

from __future__ import annotations

from pathlib import Path

import pytest


def _fail_result(step_name, work_dir):
    class R:
        ok = False
        label = step_name

        def __init__(self, wd):
            self.work_dir = str(wd)
            self.error = "killed"

    return R(work_dir)


def _ok_result(step_name, work_dir):
    class R:
        ok = True
        label = step_name
        geometry = "C 0 0 0\nH 1.09 0 0 0\n"
        artifacts = None

        def __init__(self, wd):
            self.work_dir = str(wd)

    return R(work_dir)


def _chain_with_killed_attempt(tmp_path):
    from delfin.tools import chain_state

    wd = tmp_path / "chain_wd"
    wd.mkdir()
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    chain.add_step("orca_opt", "orca_opt",
                   params={"charge": 0, "bang": "PBE0 TightOpt"})
    step = chain.steps[0]
    # The killed first attempt: numbered dir like a real run wrote it,
    # state submitted (the job was still thought alive when the process
    # crashed) — the resume classifies it failed and resubmits.
    step.mark_submitted(job_id="8812345")
    step.work_dir = str(wd / "01_orca_opt")
    chain.save()
    return chain, wd


def test_retried_step_success_is_recorded(tmp_path, monkeypatch):
    chain, wd = _chain_with_killed_attempt(tmp_path)
    from delfin.tools import chain_state
    from delfin.tools.pipeline import Pipeline

    seq = {"n": 0}

    def fake_run_step(step_name, geometry=None, work_dir=None, **kw):
        seq["n"] += 1
        return _ok_result(step_name, work_dir)

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)
    monkeypatch.setattr("delfin.scheduler_client.query", lambda *a, **k: "")

    pipe = Pipeline(name="demo")
    pipe.add("orca_opt", label="orca_opt", charge=0,
             bang="PBE0 TightOpt")
    # One resume: the submitted step's job is gone -> classified failed
    # -> resubmitted (attempt 2, new dir) -> SUCCEEDS.  The success must
    # be recorded (live run #3 threw an 11-minute converged run away
    # because mark_finished was refused from 'pending').
    res2 = pipe.run(work_dir=wd, resume_chain=True)
    assert res2.ok

    step = chain_state.ChainState.load(wd / "chain_state.json").steps[0]
    assert step.state == "finished", (
        f"expected finished, got {step.state!r} "
        f"(mark_finished refused from pending?)"
    )
    assert step.attempt == 2
    assert [a["attempt"] for a in step.attempts] == [1]
    assert Path(step.work_dir) == Path(
        wd / "orca_opt" / "attempt_2")
    assert res2.ok  # the retried attempt's success is in the result


def test_retried_step_failure_is_recorded_not_crashed(tmp_path, monkeypatch):
    chain, wd = _chain_with_killed_attempt(tmp_path)
    from delfin.tools import chain_state
    from delfin.tools.pipeline import Pipeline

    def fake_run_step(step_name, geometry=None, work_dir=None, **kw):
        return _fail_result(step_name, work_dir)

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)
    monkeypatch.setattr("delfin.scheduler_client.query", lambda *a, **k: "")

    pipe = Pipeline(name="demo")
    pipe.add("orca_opt", label="orca_opt", charge=0,
             bang="PBE0 TightOpt")
    # One resume: job gone -> resubmit attempt 2 -> fails again.  The
    # failure must be recorded (mark_failed refused from 'pending' was
    # the live crash), not swallowed.
    res = pipe.run(work_dir=wd, resume_chain=True)

    step = chain_state.ChainState.load(wd / "chain_state.json").steps[0]
    assert step.state == "failed", (
        f"expected failed, got {step.state!r} "
        f"(mark_failed refused from pending?)"
    )
    assert step.attempt == 2
    assert [a["attempt"] for a in step.attempts] == [1]
    assert step.error  # the failure is recorded, not lost
    assert not res.ok


def test_second_attempt_dir_is_attempt_2(tmp_path):
    from delfin.tools import chain_state

    wd = tmp_path / "chain_wd"
    wd.mkdir()
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    chain.add_step("orca_opt", "orca_opt")
    step = chain.steps[0]
    step.mark_submitted(job_id="8812345")
    step.work_dir = str(wd / "01_orca_opt")  # numbered 1st attempt
    chain.save()
    assert chain_state.next_attempt_dir(step, wd).name == "attempt_2"
