"""JC1 L8/L9 (live rerun #5) + reviewer note on the L7 re-review.

L8: params are persisted only when a step STARTS, so a chain killed
before its later steps started left them with params {} and a resume
could never finish them ("'charge' parameter is required" on step 3).
Every step's params (+ inputs hash) must be persisted when the chain is
created / Pipeline.run first writes the state, not at submission.

L9: a failure in a later step rolled the WHOLE state back — a step that
finished with evidence during the same resume lost its fact (operator
corrected the L2 wording: "leave the state as found" is right for a
crash BEFORE anything new is proven; restore per step, not per resume).

Reviewer note (L7 follow-up): the persisted params must hold only JSON
types — _prev_artifacts carries pathlib.Path values into the params
assignment.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def _fake(step_status_for):
    """Return a run_step fake whose per-step outcome comes from
    ``step_status_for(step_name)`` -> "ok" | "fail" | "raise"."""
    from delfin.tools._types import StepResult, StepStatus

    def fake_run_step(step_name, geometry=None, work_dir=None, **kw):
        Path(work_dir).mkdir(parents=True, exist_ok=True)
        out = Path(work_dir) / "calc.out"
        out.write_text("ORCA TERMINATED NORMALLY\n", encoding="utf-8")
        what = step_status_for(step_name)
        if what == "raise":
            raise ValueError("'charge' parameter is required")
        if what == "fail":
            return StepResult(step_name=step_name, status=StepStatus.FAILED,
                              work_dir=Path(work_dir), error="killed")
        return StepResult(step_name=step_name, status=StepStatus.SUCCESS,
                          work_dir=Path(work_dir), output_file=out,
                          geometry=out)
    return fake_run_step


def _three_step_pipe(monkeypatch, wd, statuses):
    from delfin.tools import chain_state
    from delfin.tools.pipeline import Pipeline

    monkeypatch.setattr("delfin.tools._runner.run_step", _fake(statuses))
    monkeypatch.setattr("delfin.scheduler_client.query", lambda *a, **k: "")
    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", label="xtb_opt", charge=0)
    pipe.add("orca_opt", label="orca_opt", charge=0, bang="PBE0 TightOpt")
    pipe.add("orca_sp", label="orca_sp", charge=0, basis="def2-SVP")
    pipe.run(work_dir=wd, resume_chain=True)
    return chain_state.ChainState.load(wd / "chain_state.json")


def test_l8_never_started_step_has_params_and_hash(tmp_path, monkeypatch):
    """L8: a step that never started still carries its params + inputs
    hash from chain creation, so a resume can run it."""
    wd = tmp_path / "chain_wd"
    wd.mkdir()
    chain = _three_step_pipe(monkeypatch, wd,
                             lambda n: "fail" if n == "orca_opt" else "ok")
    sp = next(s for s in chain.steps if s.name == "orca_sp")
    assert sp.params, "never-started step must carry its params (L8)"
    assert sp.params.get("charge") == 0
    assert sp.params.get("basis") == "def2-SVP"
    assert sp.inputs_hash, "never-started step must carry its inputs hash (L8)"


def test_l9_verified_finish_survives_a_later_failure(tmp_path, monkeypatch):
    """L9: step 2 finishes with evidence, step 3 raises — the finish is
    a fact and must survive; the next resume skips step 2."""
    from delfin.tools import chain_state
    from delfin.tools.pipeline import Pipeline

    wd = tmp_path / "chain_wd"
    wd.mkdir()
    chain = _three_step_pipe(monkeypatch, wd,
                             lambda n: "fail" if n == "orca_opt" else "ok")

    def fake_run_step(step_name, geometry=None, work_dir=None, **kw):
        from delfin.tools._types import StepResult, StepStatus

        Path(work_dir).mkdir(parents=True, exist_ok=True)
        if step_name == "orca_sp":
            raise ValueError("'charge' parameter is required")
        out = Path(work_dir) / "calc.out"
        out.write_text("ORCA TERMINATED NORMALLY\n", encoding="utf-8")
        return StepResult(step_name=step_name, status=StepStatus.SUCCESS,
                          work_dir=Path(work_dir), output_file=out,
                          geometry=out)

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)
    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", label="xtb_opt", charge=0)
    pipe.add("orca_opt", label="orca_opt", charge=0, bang="PBE0 TightOpt")
    pipe.add("orca_sp", label="orca_sp", charge=0, basis="def2-SVP")
    with pytest.raises(Exception):
        pipe.run(work_dir=wd, resume_chain=True)

    after = chain_state.ChainState.load(wd / "chain_state.json")
    opt = next(s for s in after.steps if s.name == "orca_opt")
    # The verified finish survived the later failure (L9).
    assert opt.state == "finished"
    assert opt.evidence, "the finish's evidence must survive the later failure"


def test_reviewer_params_hold_only_json_types(tmp_path, monkeypatch):
    """Reviewer note: the persisted params hold no pathlib.Path values."""
    from delfin.tools import chain_state

    wd = tmp_path / "chain_wd"
    wd.mkdir()
    chain = _three_step_pipe(monkeypatch, wd, lambda n: "ok")
    for step in chain.steps:
        for value in step.params.values():
            assert not isinstance(value, Path), (
                f"{step.name}: param {value!r} is a Path — only JSON "
                "types belong in the persisted params"
            )


def test_json_safe_helper_converts_paths():
    from delfin.tools import chain_state

    out = chain_state.json_safe(
        {"a": Path("/x"), "b": [Path("/y"), 1], "c": {"d": Path("/z")}})
    assert out == {"a": "/x", "b": ["/y", 1], "c": {"d": "/z"}}
    assert json.dumps(out)
