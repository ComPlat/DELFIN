"""JC1 live-run finding L4: a resumed chain must hand the geometry from the
skipped step's verified output to the next step, exactly as an uninterrupted
run would.  Live run #2 failed with "geometry is required" on step 2 because
the skip path (pipeline.py chain_skip continue) never rehydrated
current_geom — an uninterrupted run sets it from result.geometry.

The control drives the real geometry kwarg that run_step receives.
"""

from __future__ import annotations

from pathlib import Path


def _chain_with_finished_step(tmp_path: Path, tool: str = "xtb_opt"):
    """A chain whose step 1 is finished with a real output dir: the final
    xyz the chaining readers use, next to the output file."""
    from delfin.tools import chain_state

    wd = tmp_path / "chain_wd"
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    step_dir = wd / "00_step_0"
    step_dir.mkdir()
    (step_dir / "step_0.out").write_text(
        "ORCA TERMINATED NORMALLY\n", encoding="utf-8"
    )
    (step_dir / "step_0.xyz").write_text(
        "2\n\nC 0.0 0.0 0.0\nH 1.09 0.0 0.0\n", encoding="utf-8"
    )
    chain.add_step("step_0", tool)
    chain.steps[0].mark_submitted(job_id="local")
    chain.steps[0].work_dir = str(step_dir)
    chain.steps[0].mark_running()
    chain.steps[0].mark_finished(evidence={
        "normal_termination": True,
        "scf_converged": True,
        "tool_reported_success": True,
    })
    chain.save()
    return chain, wd


def _ok_result(step_name, work_dir):
    from delfin.tools._types import StepResult, StepStatus
    return StepResult(step_name=step_name, status=StepStatus.SUCCESS,
                      geometry=None, work_dir=work_dir)


def test_resume_hands_geometry_from_skipped_step(tmp_path, monkeypatch):
    """L4: the NEXT step's run_step call carries the geometry read from
    the finished step's .xyz."""
    chain, wd = _chain_with_finished_step(tmp_path)
    from delfin.tools import chain_state
    from delfin.tools.pipeline import Pipeline

    chain.add_step("step_1", "xtb_opt", params={"charge": 0})
    chain.save()

    received = {}

    def fake_run_step(step_name, geometry=None, work_dir=None, **kw):
        received["geometry"] = geometry
        received["work_dir"] = work_dir
        return _ok_result(step_name, work_dir)

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)

    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", label="step_0", charge=0)
    pipe.add("xtb_opt", label="step_1", charge=0)
    res = pipe.run(work_dir=wd, resume_chain=True)
    assert res.results, "the pending step must run"
    assert received["geometry"] is not None, (
        "the skipped step's output geometry must feed the next step"
    )
    # It is the real file: the xyz next to the verified output.
    assert Path(received["geometry"]).name == "step_0.xyz"


def test_resume_without_xyz_does_not_crash(tmp_path, monkeypatch):
    """A finished step whose xyz is missing does not crash the resume —
    the next step gets no geometry from the skip path (an explicit
    geometry kwarg on the step still wins)."""
    chain, wd = _chain_with_finished_step(tmp_path)
    (wd / "00_step_0" / "step_0.xyz").unlink()
    from delfin.tools import chain_state
    from delfin.tools.pipeline import Pipeline

    received = {}

    def fake_run_step(step_name, geometry=None, work_dir=None, **kw):
        received["geometry"] = geometry
        return _ok_result(step_name, work_dir)

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)

    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", label="step_0", charge=0)
    pipe.add("xtb_opt", label="step_1", charge=0)
    res = pipe.run(work_dir=wd, resume_chain=True)
    assert res.results
    assert received["geometry"] is None


def test_explicit_geometry_wins_over_skipped_output(tmp_path, monkeypatch):
    """A step with an explicit geometry kwarg is not overridden by the
    skip path's rehydrated geometry."""
    chain, wd = _chain_with_finished_step(tmp_path)
    from delfin.tools import chain_state
    from delfin.tools.pipeline import Pipeline

    explicit = wd / "explicit.xyz"
    explicit.write_text("2\n\nC 0.0 0.0 0.0\nH 1.09 0.0 0.0\n",
                        encoding="utf-8")

    received = {}

    def fake_run_step(step_name, geometry=None, work_dir=None, **kw):
        received["geometry"] = geometry
        return _ok_result(step_name, work_dir)

    monkeypatch.setattr("delfin.tools._runner.run_step", fake_run_step)

    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", label="step_0", charge=0)
    pipe.add("xtb_opt", label="step_1", charge=0, geometry=str(explicit))
    res = pipe.run(work_dir=wd, resume_chain=True)
    assert res.results
    assert Path(received["geometry"]) == explicit
