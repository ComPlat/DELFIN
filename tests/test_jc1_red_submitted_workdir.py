"""F-D (reviewer suggestion, 13.11.): the attempt dir must be recorded at
SUBMISSION time, not after the run returns.

Today pipeline.py sets cstep.work_dir only after run_step comes back; a
process crash while the job runs leaves a persisted `submitted` step with
no work_dir, so a resume finds no output and marks the step failed +
resubmits it although the run's output exists.  Control: this file, the
persistence test red on a6cedb6d5 (work_dir None on disk), green after
the fix records the dir and saves at mark_submitted time.
"""

from pathlib import Path

import pytest

from delfin.tools import chain_state


def _fake_step(monkeypatch, tmp_path, crash=True):
    """A Pipeline.run whose run_step records the submission then crashes."""
    from delfin.tools.pipeline import Pipeline

    def fake_run_step(step_name, work_dir=None, **kwargs):
        assert work_dir is not None
        Path(work_dir).mkdir(parents=True, exist_ok=True)
        if crash:
            raise RuntimeError("simulated node death before the run returned")
        return None

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)
    pipe = Pipeline(name="t")
    pipe.add("xtb_opt")
    pipe.run(work_dir=tmp_path, resume_chain=True)


def test_submitted_step_records_work_dir_before_the_run_returns(
    monkeypatch, tmp_path
):
    chain_wd = tmp_path / "chain"
    chain_wd.mkdir()
    with pytest.raises(RuntimeError):
        _fake_step(monkeypatch, chain_wd, crash=True)
    state = chain_state.ChainState.load(chain_wd / "chain_state.json")
    assert len(state.steps) == 1
    step = state.steps[0]
    # The submission is on disk with the attempt dir it runs in.
    assert step.state == "submitted"
    assert step.work_dir is not None, (
        "a submitted step must carry its work dir at submission time, "
        "so a crash mid-run still leaves the output findable on resume"
    )


def test_submitted_work_dir_makes_resume_find_the_output(monkeypatch, tmp_path):
    """With the dir recorded at submission, a vanished job whose output
    exists (normal termination marker written by the run) is finished,
    not failed+resubmitted."""
    from delfin.tools import chain_state as cs

    chain_wd = tmp_path / "chain"
    chain_wd.mkdir()
    try:
        _fake_step(monkeypatch, chain_wd, crash=True)
    except RuntimeError:
        pass
    state = chain_state.ChainState.load(chain_wd / "chain_state.json")
    step = state.steps[0]
    # The job actually completed and wrote its output before the crash.
    out = Path(step.work_dir) / "xtb_opt.out"
    out.write_text("normal termination\n")
    plan = cs.build_resume_plan(state, job_state=lambda j: None)
    assert plan.get("wait") is False
    assert plan.get("complete") is True, (
        f"a submitted step with recorded work_dir and real output proof "
        f"must be finished on resume, plan was: {plan}"
    )
