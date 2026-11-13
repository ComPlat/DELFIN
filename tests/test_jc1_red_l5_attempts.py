"""JC1 live-run finding L5: attempt bookkeeping after a resume of a chain
with one killed attempt.  Expected (operator):
- attempts[] lists the killed first attempt (its dir, end = killed/job gone)
- the new attempt is 2
- the recorded work_dir is the dir actually used (a NEW attempt dir, the
  old calculation is never overwritten).
"""

from __future__ import annotations

from pathlib import Path


def _chain_with_killed_attempt(tmp_path: Path):
    """Step 1 finished + verified; step 2 was submitted, its job died,
    no proof in its output -> the resume must resubmit it."""
    from delfin.tools import chain_state

    wd = tmp_path / "chain_wd"
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    s0 = wd / "00_step_0"
    s0.mkdir()
    (s0 / "step_0.out").write_text("ORCA TERMINATED NORMALLY\n",
                                   encoding="utf-8")
    (s0 / "step_0.xyz").write_text("2\n\nC 0.0 0.0 0.0\nH 1.09 0.0 0.0\n",
                                   encoding="utf-8")
    chain.add_step("step_0", "xtb_opt")
    chain.steps[0].mark_submitted(job_id="local")
    chain.steps[0].work_dir = str(s0)
    chain.steps[0].mark_running()
    chain.steps[0].mark_finished(evidence={
        "normal_termination": True, "scf_converged": True,
        "tool_reported_success": True,
    })

    s1 = wd / "01_orca_opt"
    s1.mkdir()
    (s1 / "orca_opt.out").write_text("killed mid-run\n", encoding="utf-8")
    chain.add_step("orca_opt", "orca_opt", params={"charge": 0})
    chain.steps[1].mark_submitted(job_id="8716164")
    chain.steps[1].work_dir = str(s1)
    chain.save()
    return chain, wd


def _fail_result(step_name, work_dir):
    from delfin.tools._types import StepResult, StepStatus
    return StepResult(step_name=step_name, status=StepStatus.FAILED,
                      error="job gone", work_dir=work_dir)


def test_resume_records_killed_attempt_and_new_attempt(tmp_path, monkeypatch):
    chain, wd = _chain_with_killed_attempt(tmp_path)
    from delfin.tools import chain_state
    from delfin.tools.pipeline import Pipeline

    received = {}

    def fake_run_step(step_name, geometry=None, work_dir=None, **kw):
        received["work_dir"] = work_dir
        return _fail_result(step_name, work_dir)

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)
    # Pin the scheduler answer: the job is GONE (the throttle under a
    # loaded gate process must not decide this test's branch).
    monkeypatch.setattr("delfin.scheduler_client.query",
                        lambda *a, **k: "")

    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", label="step_0", charge=0)
    pipe.add("orca_opt", label="orca_opt", charge=0)
    pipe.run(work_dir=wd, resume_chain=True)

    state = chain_state.ChainState.load(wd / "chain_state.json")
    step = state.step("orca_opt")
    # The killed first attempt is listed with its dir and its end.
    assert step.attempts, "the killed attempt must be recorded"
    first = step.attempts[0]
    assert first["work_dir"].endswith("01_orca_opt")
    assert first["attempt"] == 1
    assert first["end"]
    # The new attempt is 2.
    assert step.attempt == 2
    # The recorded work_dir is the dir actually used — a NEW dir, not the
    # killed attempt's (never overwrite a calculation).
    assert received["work_dir"] is not None
    assert Path(step.work_dir) == Path(received["work_dir"])
    assert Path(step.work_dir) != wd / "orca_opt" / "01_orca_opt"
    # The old dir still exists, untouched (kept, not deleted).
    assert (wd / "01_orca_opt" / "orca_opt.out").is_file()
