"""JC2 reviewer adversarial tests, part 1: persisted params stay clean.

JC1's live run #6 showed orca_sp's persisted params holding pipeline-internal
wiring: ``_prev_artifacts`` (pipeline.py:1428) and the absolute ``moread``
path autowire_kwargs injects (_wiring.py).  Per the JC2 brief, persisted
params must be the USER's step parameters only — what ``Pipeline.add(...)``
was given plus the charge/mult defaults the pipeline applies — and
``inputs_hash`` must be identical for the same inputs regardless of machine
or directory.

Control: this file, red on main (026eba8f0) because the submission path
still persists merged/autowired kwargs (pipeline.py:1523).
"""

from pathlib import Path

import pytest

from delfin.tools import chain_state


def _fake_run_step(step_name, geometry=None, work_dir=None, **kwargs):
    """A stand-in step that writes a non-empty output like a real run."""
    from delfin.tools.pipeline import StepResult, StepStatus

    out = Path(work_dir) / f"{step_name}.out"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("fake run output\n")
    return StepResult(
        step_name=step_name,
        status=StepStatus.SUCCESS,
        # Like a real ORCA adapter: the run produces a .gbw the next
        # orca_* step consumes via the gbw→moread wiring rule.
        artifacts={"gbw": out.with_suffix(".gbw")},
    )


def _run_two_step_chain(monkeypatch, wd):
    """A real Pipeline.run over two steps in *wd*; returns the chain state."""
    monkeypatch.setattr("delfin.tools._runner.run_step", _fake_run_step)
    from delfin.tools.pipeline import Pipeline

    pipe = Pipeline(name="jc2-params")
    pipe.add("xtb_opt", charge=0)
    # An orca step downstream of another orca step is where the gbw→moread
    # wiring fires (_wiring.DEFAULT_RULES) and _prev_artifacts is injected.
    pipe.add("orca_sp", charge=0, method="B3LYP", basis="def2-SVP")
    pipe.run(work_dir=wd, resume_chain=True)
    return chain_state.ChainState.load(wd / "chain_state.json")


def test_persisted_params_hold_only_user_params(monkeypatch, tmp_path):
    state = _run_two_step_chain(monkeypatch, tmp_path / "a")
    sp = state.step("orca_sp")
    internal = {
        k for k in sp.params
        if k.startswith("_") or k in ("moread", "gbw", "hess_file")
    }
    assert not internal, (
        f"persisted params must be the user's step parameters only; got "
        f"pipeline-internal wiring keys {sorted(internal)} (full params: "
        f"{sorted(sp.params)})"
    )


def test_inputs_hash_is_directory_independent(monkeypatch, tmp_path):
    state_a = _run_two_step_chain(monkeypatch, tmp_path / "a")
    state_b = _run_two_step_chain(monkeypatch, tmp_path / "b")
    for name in ("xtb_opt", "orca_sp"):
        ha = state_a.step(name).inputs_hash
        hb = state_b.step(name).inputs_hash
        assert ha == hb, (
            f"inputs_hash of step {name!r} differs between two directories "
            f"({ha} vs {hb}): it must hash (tool, user params, geometry "
            f"content) only, never absolute paths of one machine"
        )


def test_load_of_jc1_state_drops_internal_keys_with_one_log_line(
    monkeypatch, tmp_path, caplog
):
    """A chain_state.json written by JC1 (with the internal keys) still
    loads and resumes; the internal keys are dropped on load, one log
    line."""
    wd = tmp_path / "legacy"
    wd.mkdir()
    state = {
        "name": "legacy-chain",
        "work_dir": str(wd),
        "steps": [
            {
                "name": "orca_sp",
                "tool": "orca_sp",
                "state": "pending",
                "params": {
                    "charge": 0,
                    "method": "B3LYP",
                    "basis": "def2-SVP",
                    "moread": "/abs/path/from/another/machine/ref.gbw",
                    "_prev_artifacts": {"gbw": "/abs/path/ref.gbw"},
                },
                "inputs_hash": "legacy-hash",
            }
        ],
    }
    (wd / "chain_state.json").write_text(
        __import__("json").dumps(state), encoding="utf-8"
    )
    import logging

    with caplog.at_level(logging.INFO):
        loaded = chain_state.ChainState.load(wd / "chain_state.json")
    sp = loaded.step("orca_sp")
    assert "_prev_artifacts" not in sp.params, (
        "internal keys must be dropped on load so a legacy state resumes "
        "with clean params"
    )
    assert "moread" not in sp.params
    assert sp.params.get("charge") == 0, "user params survive the load"
    assert any(
        "internal" in r.message.lower() or "_prev_artifacts" in r.message
        for r in caplog.records
    ), "dropping internal keys must leave one log line"
