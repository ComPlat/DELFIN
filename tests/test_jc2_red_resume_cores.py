"""The plan's builder requirement (live-run plan, Job B step 5):
``delfin chain resume`` defaults to --cores 1, so a chain started with
N cores resumes with 1 — a silently changed PAL.  The offer must carry
the chain's cores (``--cores N``) or the resume must take them from the
state; a test pins that a resume never silently changes PAL.

Control: RED on 56f8c176c — nothing records the cores a chain was
submitted with, the offer carries no --cores, and the resume default is
1 regardless of the chain.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture()
def home(tmp_path, monkeypatch):
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path / "home")
    import delfin.agent.session_store as ss
    monkeypatch.setattr(
        ss, "_SESSIONS_DIR", tmp_path / "home" / ".delfin" / "agent_sessions")
    return tmp_path / "home"


def _fake_run_step(step_name, geometry=None, work_dir=None, cores=1, **kwargs):
    from delfin.tools.pipeline import StepResult, StepStatus

    out = Path(work_dir) / f"{step_name}.out"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("fake run output\n")
    return StepResult(step_name=step_name, status=StepStatus.SUCCESS,
                      work_dir=Path(work_dir))


def _run_chain(wd, monkeypatch, cores):
    monkeypatch.setattr("delfin.tools._runner.run_step", _fake_run_step)
    from delfin.tools.pipeline import Pipeline

    geom = Path(wd) / "start.xyz"
    geom.parent.mkdir(parents=True, exist_ok=True)
    geom.write_text("2\n\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", charge=0)
    pipe.run(work_dir=wd, geometry=geom, cores=cores)
    return wd


def test_chain_records_its_submitted_cores(home, monkeypatch):
    from delfin.tools import chain_state

    wd = _run_chain(home / "runs" / "a", monkeypatch, cores=4)
    chain = chain_state.ChainState.load(wd / "chain_state.json")
    recorded = {s.cores for s in chain.steps}
    assert recorded == {4}, (
        f"the chain must record the cores each step was submitted with "
        f"(got {recorded}) so a resume never silently changes PAL")


def test_resume_takes_cores_from_the_state(monkeypatch, tmp_path):
    """`delfin chain resume` with NO --cores flag uses the chain's
    recorded cores, never a silent 1."""

    from delfin.tools import chain_state
    from delfin.cli_chain import main

    wd = tmp_path / "chain"
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    chain.add_step("xtb_opt", "xtb_opt", params={"charge": 0})
    chain.steps[0].cores = 4
    chain.save()

    captured = {}

    def fake_run(self, *, cores=1, **kwargs):
        captured["cores"] = cores
        from delfin.tools.pipeline import PipelineResult
        return PipelineResult(name="demo", results=[], branch_results={})

    monkeypatch.setattr("delfin.tools.pipeline.Pipeline.run", fake_run)
    rc = main(["resume", str(wd)])
    assert rc == 0
    assert captured["cores"] == 4, (
        "a resume without --cores must use the chain's recorded cores — "
        "PAL must never change silently")


def test_explicit_cores_flag_wins(monkeypatch, tmp_path):
    from delfin.tools import chain_state
    from delfin.cli_chain import main

    wd = tmp_path / "chain"
    chain = chain_state.ChainState.create(name="demo", work_dir=wd)
    chain.add_step("xtb_opt", "xtb_opt", params={"charge": 0})
    chain.steps[0].cores = 4
    chain.save()

    captured = {}

    def fake_run(self, *, cores=1, **kwargs):
        captured["cores"] = cores
        from delfin.tools.pipeline import PipelineResult
        return PipelineResult(name="demo", results=[], branch_results={})

    monkeypatch.setattr("delfin.tools.pipeline.Pipeline.run", fake_run)
    rc = main(["resume", str(wd), "--cores", "2"])
    assert rc == 0
    assert captured["cores"] == 2, "an explicit flag wins over the state"
