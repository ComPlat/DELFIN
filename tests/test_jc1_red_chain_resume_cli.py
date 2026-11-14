"""JC1 phase 3: the explicit ``delfin chain resume`` path, driven through
its public entry point (:func:`delfin.cli_chain.main`) with a fake
``run_step`` via monkeypatch (operator rule: probe tests go through the
gate, never a python -c smoke check).

Cases:
1. resume of a pending chain actually runs the pending step through
   ``run_step`` — and only because the resume command was given;
2. resume of a chain whose step finished is a no-op (no rerun of
   finished work — JC1 target 1);
3. resume without a state file is a clean error.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def _chain_with_state(tmp_path: Path, tool: str = "xtb_opt"):
    from delfin.tools import chain_state

    wd = tmp_path / "chain_wd"
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    chain.add_step("step_0", tool)
    chain.save()
    return chain, wd


def _patched_runner(monkeypatch, calls, tool="xtb_opt"):
    """Fake run_step as a module attribute of delfin.tools._runner."""
    import delfin.tools._runner as runner_mod
    from delfin.tools._types import StepResult, StepStatus

    def fake_run_step(step_name, *, geometry=None, cores=1, work_dir=None,
                      **kw):
        calls.append(step_name)
        wd = Path(work_dir)
        wd.mkdir(parents=True, exist_ok=True)
        (wd / f"{step_name}.out").write_text(
            "ORCA TERMINATED NORMALLY\nSCF CONVERGED\n"
        )
        return StepResult(
            step_name=step_name, status=StepStatus.SUCCESS,
            work_dir=wd, data={},
        )

    monkeypatch.setattr(runner_mod, "run_step", fake_run_step)


def test_resume_runs_the_pending_step_only_when_asked(tmp_path, monkeypatch):
    chain, wd = _chain_with_state(tmp_path)
    calls: list = []
    _patched_runner(monkeypatch, calls)

    from delfin.cli_chain import main as chain_main

    rc = chain_main(["resume", str(wd)])
    assert rc == 0
    assert calls == [chain.steps[0].tool], (
        "resume must actually run the pending step through run_step, "
        f"and only after the explicit command (got {calls})"
    )
    from delfin.tools import chain_state

    reloaded = chain_state.ChainState.load(wd / "chain_state.json")
    assert reloaded.steps[0].state == "finished"


def test_resume_of_a_finished_chain_is_a_no_op(tmp_path, monkeypatch):
    from delfin.tools import chain_state

    chain, wd = _chain_with_state(tmp_path)
    # Mark the step finished with full ORCA-marker evidence and a real
    # output so re-verification can prove it.
    step = chain.steps[0]
    step.mark_submitted(job_id="1")
    step.mark_running()
    step_dir = wd / "00_step_0"
    step_dir.mkdir(parents=True)
    (step_dir / "step_0.out").write_text(
        "ORCA TERMINATED NORMALLY\nSCF CONVERGED\n"
    )
    step.mark_finished(evidence={
        "normal_termination": True, "scf_converged": True,
    })
    chain.save()

    calls: list = []
    _patched_runner(monkeypatch, calls)
    from delfin.cli_chain import main as chain_main

    rc = chain_main(["resume", str(wd)])
    assert rc == 0
    assert calls == [], (
        "a verified finished step must not be rerun on resume"
    )


def test_resume_without_state_file_is_a_clean_error(tmp_path, capsys):
    from delfin.cli_chain import main as chain_main

    rc = chain_main(["resume", str(tmp_path / "nowhere")])
    assert rc != 0
    assert "no chain state file" in capsys.readouterr().out
