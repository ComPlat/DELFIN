"""JC1 phase-3 red control: bounded resubmission + Pipeline wiring + CLI
(builder nacht-s5).

Pins JC1 targets 3-5. Red on 1ba4744e9, where the resume decider exists but
Pipeline.run neither writes nor reads the chain state, no bounded
resubmission through the recalc-style attempt dirs exists, and the
``delfin chain`` subcommand is unregistered.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def _module_available(name: str) -> bool:
    try:
        __import__(name)
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# 3. bounded resubmission
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _module_available("delfin.tools.chain_state"),
    reason="chain_state not present",
)
class TestBoundedResubmission:
    def _chain(self, tmp_path: Path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = chain_state.ChainState.create(name="demo", work_dir=wd)
        for i in range(2):
            chain.add_step(f"step_{i}", "orca_opt")
        return chain

    def test_failed_step_resubmitted_in_new_attempt_dir(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        old_dir = wd / "step_0_attempt_1"
        old_dir.mkdir()
        (old_dir / "step_0.out").write_text("CARTESIAN COORDINATES\n")
        chain.steps[0].work_dir = str(old_dir)
        chain.steps[0].mark_submitted(job_id="1")
        chain.steps[0].mark_running()
        chain.steps[0].mark_failed("SCF did not converge", "convergence")
        chain.save()

        plan = chain_state.build_resume_plan(chain, work_dir=wd)
        assert plan["next_step"] is chain.steps[0]
        new_dir = chain_state.next_attempt_dir(chain.steps[0], wd)
        assert new_dir != Path(chain.steps[0].work_dir)
        assert not new_dir.exists() or list(new_dir.iterdir()) == []
        # The old attempt is kept: no overwrite.
        assert (old_dir / "step_0.out").exists()

    def test_bounded_attempts_stop_the_chain_with_reason(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        step = chain.steps[0]
        step.mark_submitted(job_id="1")
        step.mark_running()
        step.mark_failed("missing parameter", "missing_param")
        chain.save()

        decision = chain_state.resubmission_decision(step, max_attempts=3)
        # A non-recoverable cause (per _recovery.suggest_recovery) stops the
        # chain — it is not retried forever.
        assert decision.action == "stop"
        assert decision.reason

    def test_recoverable_cause_allows_a_bounded_retry(self, tmp_path):
        from delfin.tools import chain_state

        chain = self._chain(tmp_path)
        step = chain.steps[0]
        step.mark_submitted(job_id="1")
        step.mark_running()
        step.mark_failed("SCF did not converge", "convergence")

        decision = chain_state.resubmission_decision(step, max_attempts=3)
        assert decision.action == "retry"

    def test_retry_budget_exhausted_stops(self, tmp_path):
        from delfin.tools import chain_state

        chain = self._chain(tmp_path)
        step = chain.steps[0]
        step.mark_submitted(job_id="1")
        step.mark_running()
        step.mark_failed("SCF did not converge", "convergence")
        step.attempt = 3

        decision = chain_state.resubmission_decision(step, max_attempts=3)
        assert decision.action == "stop"
        assert "attempt" in decision.reason.lower()

    def test_pal_maxcore_time_never_changed(self, tmp_path):
        from delfin.tools import chain_state

        chain = self._chain(tmp_path)
        step = chain.steps[0]
        step.mark_submitted(job_id="1")
        step.mark_running()
        step.mark_failed("SCF did not converge", "convergence")

        decision = chain_state.resubmission_decision(step, max_attempts=3)
        assert decision.kwargs_overrides == {} or all(
            k not in ("PAL", "maxcore", "time") for k in decision.kwargs_overrides
        )


# ---------------------------------------------------------------------------
# 4. Pipeline.run writes and reads the chain state (one mechanism)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _module_available("delfin.tools.pipeline"),
    reason="pipeline not present",
)
class TestPipelineChainWiring:
    def _pipeline(self):
        from delfin.tools.pipeline import Pipeline

        pipe = Pipeline(name="jc1-chain")
        pipe.add("xtb_opt", charge=0)
        pipe.add("orca_sp", charge=0)
        return pipe

    def test_run_writes_chain_state(self, tmp_path, monkeypatch):
        from delfin.tools import chain_state
        from delfin.tools._types import StepResult, StepStatus

        def fake_run_step(step_name, *, geometry=None, cores=1, work_dir=None, **kw):
            wd = Path(work_dir)
            wd.mkdir(parents=True, exist_ok=True)
            return StepResult(
                step_name=step_name, status=StepStatus.SUCCESS,
                work_dir=wd, data={"energy_Eh": -1.0},
            )

        monkeypatch.setattr(
            "delfin.tools._runner.run_step", fake_run_step
        )
        pipe = self._pipeline()
        pipe.run(cores=1, work_dir=tmp_path)
        state_path = tmp_path / "chain_state.json"
        assert state_path.is_file(), (
            "Pipeline.run must persist the chain state next to its work dir"
        )
        chain = chain_state.ChainState.load(state_path)
        assert [s.state for s in chain.steps] == ["finished", "finished"]
        assert all(
            s.evidence is not None for s in chain.steps
        )

    def test_resume_from_chain_state_skips_verified_finished(
        self, tmp_path, monkeypatch
    ):
        from delfin.tools import chain_state
        from delfin.tools._types import StepResult, StepStatus

        calls = []

        def fake_run_step(step_name, *, geometry=None, cores=1, work_dir=None, **kw):
            calls.append(step_name)
            wd = Path(work_dir)
            wd.mkdir(parents=True, exist_ok=True)
            (wd / f"{step_name}.out").write_text(
                "ORCA TERMINATED NORMALLY\nSCF CONVERGED\n"
            )
            return StepResult(
                step_name=step_name, status=StepStatus.SUCCESS,
                work_dir=wd, data={"energy_Eh": -1.0},
            )

        monkeypatch.setattr(
            "delfin.tools._runner.run_step", fake_run_step
        )
        pipe = self._pipeline()
        pipe.run(cores=1, work_dir=tmp_path)
        assert calls == ["xtb_opt", "orca_sp"]

        # Resume: step 1's output still proves it finished -> skipped.
        calls.clear()
        pipe.run(cores=1, work_dir=tmp_path, resume_chain=True)
        assert calls == [], (
            "a verified finished step must not be rerun on resume"
        )

    def test_resume_reruns_a_step_whose_output_lost_the_proof(
        self, tmp_path, monkeypatch
    ):
        from delfin.tools import chain_state
        from delfin.tools._types import StepResult, StepStatus

        def fake_run_step(step_name, *, geometry=None, cores=1, work_dir=None, **kw):
            wd = Path(work_dir)
            wd.mkdir(parents=True, exist_ok=True)
            (wd / f"{step_name}.out").write_text(
                "ORCA TERMINATED NORMALLY\nSCF CONVERGED\n"
            )
            return StepResult(
                step_name=step_name, status=StepStatus.SUCCESS,
                work_dir=wd, data={"energy_Eh": -1.0},
            )

        monkeypatch.setattr(
            "delfin.tools._runner.run_step", fake_run_step
        )
        pipe = self._pipeline()
        pipe.run(cores=1, work_dir=tmp_path)

        # Corrupt step 1's output (proof gone — e.g. node crash mid-write).
        out = tmp_path / "01_xtb_opt" / "xtb_opt.out"
        out.write_text("CARTESIAN COORDINATES\n")

        calls = []

        def fake_run_step2(step_name, *, geometry=None, cores=1, work_dir=None, **kw):
            calls.append(step_name)
            wd = Path(work_dir)
            wd.mkdir(parents=True, exist_ok=True)
            (wd / f"{step_name}.out").write_text(
                "ORCA TERMINATED NORMALLY\nSCF CONVERGED\n"
            )
            return StepResult(
                step_name=step_name, status=StepStatus.SUCCESS,
                work_dir=wd, data={"energy_Eh": -2.0},
            )

        monkeypatch.setattr(
            "delfin.tools._runner.run_step", fake_run_step2
        )
        pipe.run(cores=1, work_dir=tmp_path, resume_chain=True)
        assert "xtb_opt" in calls, (
            "a step whose output lost the proof must be rerun, not skipped"
        )


# ---------------------------------------------------------------------------
# 5. entry point: delfin chain status|resume
# ---------------------------------------------------------------------------


class TestChainCli:
    def test_cli_chain_module_exists(self):
        try:
            import delfin.cli_chain  # noqa: F401
        except ImportError as exc:
            pytest.fail(f"delfin.cli_chain missing: {exc}")

    def test_delfin_chain_status_lists_steps(self, tmp_path, capsys):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = chain_state.ChainState.create(name="demo", work_dir=wd)
        chain.add_step("step_0", "xtb_opt")
        chain.steps[0].mark_submitted(job_id="42")
        chain.save()

        from delfin.cli_chain import main as chain_main

        rc = chain_main(["status", str(wd)])
        assert rc == 0
        out = capsys.readouterr().out
        assert "step_0" in out
        assert "submitted" in out

    def test_delfin_chain_resume_is_explicit_not_automatic(self, tmp_path):
        from delfin.tools import chain_state
        from delfin.tools import pipeline as pipeline_mod
        from delfin.tools._types import StepResult, StepStatus

        wd = tmp_path / "chain_wd"
        chain = chain_state.ChainState.create(name="demo", work_dir=wd)
        chain.add_step("step_0", "xtb_opt")

        calls = []

        def fake_run_step(step_name, *, geometry=None, cores=1, work_dir=None, **kw):
            calls.append(step_name)
            wd2 = Path(work_dir)
            wd2.mkdir(parents=True, exist_ok=True)
            (wd2 / f"{step_name}.out").write_text(
                "ORCA TERMINATED NORMALLY\nSCF CONVERGED\n"
            )
            return StepResult(
                step_name=step_name, status=StepStatus.SUCCESS,
                work_dir=wd2, data={},
            )

        import delfin.cli_chain as cli_chain

        import delfin.tools._runner as runner_mod

        original = runner_mod.run_step
        runner_mod.run_step = fake_run_step
        try:
            rc = cli_chain.main(["resume", str(wd)])
        finally:
            runner_mod.run_step = original
        assert rc == 0
        assert calls == ["step_0"], (
            "resume must actually run the pending step — and only after the "
            "explicit resume command (never on its own)"
        )

    def test_status_of_a_chain_with_no_state_file_is_a_clean_error(
        self, tmp_path, capsys
    ):
        from delfin.cli_chain import main as chain_main

        rc = chain_main(["status", str(tmp_path / "nowhere")])
        assert rc != 0
