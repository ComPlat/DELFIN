"""Red controls for the operator's live-chain-run findings L1/L2/L3.

Code under test: agent/s5-jc1b @ 51076d69f (the snapshot the operator
ran).  All four cases must be RED before the fix (L1, L2, L3, NOTE).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.tools import chain_state as cs
from delfin.tools import pipeline as pipeline_mod
from delfin.tools.pipeline import Pipeline


def _fake_run_step(monkeypatch, calls):
    """A run_step stand-in that RECORDS the kwargs it receives (unlike the
    fakes that ignore kwargs — that is why L1 passed the unit tests)."""

    def fake(step_name, geometry=None, cores=None, work_dir=None, **kwargs):
        calls.append({
            "step_name": step_name,
            "geometry": str(geometry) if geometry else None,
            "work_dir": str(work_dir) if work_dir else None,
            "kwargs": kwargs,
        })
        from delfin.tools._types import StepResult, StepStatus

        return StepResult(step_name=step_name, status=StepStatus.SUCCESS)

    monkeypatch.setattr(pipeline_mod, "run_step", fake, raising=False)
    monkeypatch.setattr("delfin.tools._runner.run_step", fake, raising=False)
    return fake


def _run_two_step_chain(monkeypatch, tmp_path, calls):
    _fake_run_step(monkeypatch, calls)
    wd = tmp_path / "chain"
    pipe = Pipeline(name="live")
    pipe.add("xtb_opt", charge=0, mult=1)
    pipe.add("orca_sp", charge=0, mult=1, method="PBE0", basis="def2-SVP")
    pipe.run(work_dir=wd, cores=2)
    state = json.loads((wd / "chain_state.json").read_text())
    return wd, state


def test_l1_step_params_persisted_and_validate_against_the_real_adapter(
    monkeypatch, tmp_path
):
    """L1: the persisted state must carry the step's parameters so the
    resume can rebuild a REAL run (charge/method/basis, ...), and the
    persisted params must satisfy the real adapter's validate_params."""
    wd, state = _run_two_step_chain(
        monkeypatch, tmp_path, []
    )
    orca = [s for s in state["steps"] if s["tool"] == "orca_sp"][0]
    assert orca.get("params"), "step params not persisted in chain state"
    assert orca["params"]["method"] == "PBE0"

    from delfin.tools._registry import get as get_adapter

    adapter = get_adapter("orca_sp")
    assert adapter is not None
    # The REAL validation the resume rebuild must satisfy:
    adapter.validate_params(**orca["params"])


def test_l2_failed_resume_leaves_state_untouched(monkeypatch, tmp_path):
    """L2: a resume that raises must leave chain_state.json EXACTLY as it
    found it (the crashed run's submitted/attempt facts survive)."""
    wd, state = _run_two_step_chain(monkeypatch, tmp_path, [])
    before = (wd / "chain_state.json").read_bytes()
    loaded = cs.ChainState.load(wd / "chain_state.json")
    loaded.steps[0].state = "submitted"
    loaded.steps[0].job_id = "8716164"
    loaded.steps[0].work_dir = str(wd / "01_orca_opt")
    loaded.save()
    before = (wd / "chain_state.json").read_bytes()

    # Build the CLI resume path and let Pipeline.run raise the operator's
    # exact error before step 2 executes.
    from delfin.tools import pipeline as pl

    def boom(self, *a, **kw):
        raise ValueError("'charge' parameter is required")

    monkeypatch.setattr(pl.Pipeline, "run", boom)
    from delfin import cli_chain

    rc = cli_chain.main(["resume", str(wd)])
    assert rc == 1, "resume must exit 1, not crash"
    after = (wd / "chain_state.json").read_bytes()
    assert after == before, "failed resume mutated the persisted state"


def test_l3_inputs_hash_is_stable_and_nonempty(monkeypatch, tmp_path):
    """L3: inputs_hash must be computed from the step's actual inputs and
    be stable for identical inputs."""
    wd, state = _run_two_step_chain(monkeypatch, tmp_path, [])
    hashes = {s["name"]: s.get("inputs_hash") for s in state["steps"]}
    assert all(hashes.values()), f"inputs_hash empty: {hashes}"

    # Same inputs -> same hash (stability)
    wd2 = tmp_path / "chain2"
    pipe = Pipeline(name="live2")
    pipe.add("xtb_opt", charge=0, mult=1)
    pipe.add("orca_sp", charge=0, mult=1, method="PBE0", basis="def2-SVP")
    pipe.run(work_dir=wd2, cores=2)
    state2 = json.loads((wd2 / "chain_state.json").read_text())
    hashes2 = {s["name"]: s.get("inputs_hash") for s in state2["steps"]}
    assert hashes == hashes2


def test_note_resume_crash_is_one_clear_line_not_a_traceback(
    monkeypatch, tmp_path, capsys
):
    """NOTE: `delfin chain resume` must exit non-zero with one clear line
    when the resume raises."""
    wd, state = _run_two_step_chain(monkeypatch, tmp_path, [])
    from delfin.tools import pipeline as pl

    def boom(self, *a, **kw):
        raise ValueError("'charge' parameter is required")

    monkeypatch.setattr(pl.Pipeline, "run", boom)
    from delfin import cli_chain

    rc = cli_chain.main(["resume", str(wd)])
    out = capsys.readouterr().out
    assert rc == 1
    assert "'charge' parameter is required" in out
    assert "Traceback" not in out
