"""JC1 phase-2 red control: resume decides from facts, not from the state
file (builder nacht-s5).

Pins JC1 target 2. Red on 65f5c8eaf, where ChainState exists but no
resume-decision logic does.
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


def test_resume_decider_module_surface_exists():
    from delfin.tools import chain_state

    for name in ("verify_finished_step", "classify_running_step"):
        assert hasattr(chain_state, name), (
            f"chain_state.{name} missing — resume must decide from facts "
            "(JC1 target 2)"
        )


@pytest.mark.skipif(
    not _module_available("delfin.tools.chain_state"),
    reason="chain_state not present",
)
class TestResumeDecidesFromFacts:
    """A step marked finished is re-verified from its ORCA output before it
    is skipped; a step marked running is checked against the scheduler."""

    def _chain(self, tmp_path: Path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = chain_state.ChainState.create(name="demo", work_dir=wd)
        for i in range(3):
            chain.add_step(f"step_{i}", "orca_opt")
        return chain

    # -- finished re-verification ------------------------------------

    def _mark_finished(self, chain, work_dir: Path):
        chain.steps[0].work_dir = str(work_dir)
        chain.steps[0].mark_submitted(job_id="1")
        chain.steps[0].mark_running()
        chain.steps[0].mark_finished(
            evidence={"normal_termination": True, "scf_converged": True}
        )
        chain.save()

    def test_finished_step_with_good_output_is_verified(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        step_dir = wd / "step_0"
        step_dir.mkdir()
        (step_dir / "step_0.out").write_text(
            "ORCA TERMINATED NORMALLY\nSCF CONVERGED\n"
        )
        self._mark_finished(chain, step_dir)
        verdict = chain_state.verify_finished_step(chain.steps[0], wd)
        assert verdict.verified is True

    def test_finished_step_whose_output_is_missing_is_downgraded(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        empty = wd / "step_0"
        empty.mkdir()
        self._mark_finished(chain, empty)
        verdict = chain_state.verify_finished_step(chain.steps[0], wd)
        assert verdict.verified is False
        assert verdict.downgraded_to == "pending"

    def test_finished_step_without_scf_convergence_is_downgraded(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        step_dir = wd / "step_0"
        step_dir.mkdir()
        (step_dir / "step_0.out").write_text(
            "ORCA TERMINATED NORMALLY\nSCF not converged\n"
        )
        self._mark_finished(chain, step_dir)
        verdict = chain_state.verify_finished_step(chain.steps[0], wd)
        assert verdict.verified is False
        assert verdict.downgraded_to == "pending"

    def test_downgrade_is_persisted_on_save(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        empty = wd / "step_0"
        empty.mkdir()
        self._mark_finished(chain, empty)
        verdict = chain_state.verify_finished_step(chain.steps[0], wd)
        assert verdict.downgraded_to == "pending"
        chain.save()
        reloaded = chain_state.ChainState.load(chain.state_path)
        assert reloaded.steps[0].state == "pending"
        assert reloaded.steps[0].evidence is None

    def test_downgraded_step_yields_resume_plan_skip_false(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        empty = wd / "step_0"
        empty.mkdir()
        self._mark_finished(chain, empty)
        chain_state.verify_finished_step(chain.steps[0], wd)
        plan = chain_state.build_resume_plan(chain)
        assert plan["steps_to_skip"] == []
        assert plan["next_step"] == chain.steps[0]

    # -- running step vs scheduler ------------------------------------

    def test_running_step_still_queued_or_running_waits(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        chain.steps[0].mark_submitted(job_id="4242")
        chain.steps[0].mark_running()
        chain.save()

        verdict = chain_state.classify_running_step(
            chain.steps[0], job_state=lambda jid: "RUNNING"
        )
        assert verdict.action == "wait"

    def test_running_step_gone_from_scheduler_is_classified(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        chain.steps[0].mark_submitted(job_id="4242")
        chain.steps[0].mark_running()
        chain.save()

        verdict = chain_state.classify_running_step(
            chain.steps[0], job_state=lambda jid: None
        )
        # Job gone: the end is classified from the step's own output via
        # the recovery vocabulary, not guessed.
        assert verdict.action == "classify"
        assert verdict.error_kind is not None
        assert verdict.error_kind != "none"

    def test_running_step_gone_with_dead_output_is_failed(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        chain.steps[0].mark_submitted(job_id="4242")
        chain.steps[0].mark_running()
        chain.save()

        step_dir = wd / "step_0"
        step_dir.mkdir()
        (step_dir / "step_0.out").write_text(
            "CARTESIAN COORDINATES\n...some partial output...\n"
        )
        verdict = chain_state.classify_running_step(chain.steps[0])
        assert verdict.action == "mark_failed"
        assert chain.steps[0].state == "failed"

    def test_throttled_scheduler_answer_is_not_a_fact(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        chain.steps[0].mark_submitted(job_id="4242")
        chain.steps[0].mark_running()
        chain.save()

        from delfin.scheduler_client import THROTTLED

        verdict = chain_state.classify_running_step(
            chain.steps[0], job_state=lambda jid: THROTTLED
        )
        assert verdict.action == "wait"

    # -- resume plan ----------------------------------------------------

    def test_resume_plan_skips_only_verified_finished(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        # step_0 finished with real output, step_1 running (gone), step_2 pending
        step0 = wd / "step_0"
        step0.mkdir()
        (step0 / "step_0.out").write_text(
            "ORCA TERMINATED NORMALLY\nSCF CONVERGED\n"
        )
        chain.steps[0].work_dir = str(step0)
        chain.steps[0].mark_submitted(job_id="1")
        chain.steps[0].mark_running()
        chain.steps[0].mark_finished(
            evidence={"normal_termination": True, "scf_converged": True}
        )
        chain.steps[1].mark_submitted(job_id="2")
        chain.steps[1].mark_running()
        chain.save()

        plan = chain_state.build_resume_plan(
            chain, job_state=lambda jid: None, work_dir=wd
        )
        assert plan["steps_to_skip"] == [chain.steps[0]]
        assert plan["next_step"] is chain.steps[1]
        assert plan["wait"] is False

    def test_resume_plan_waits_on_a_live_job(self, tmp_path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = self._chain(tmp_path)
        chain.steps[0].mark_submitted(job_id="2")
        chain.steps[0].mark_running()
        chain.save()
        plan = chain_state.build_resume_plan(
            chain, job_state=lambda jid: "RUNNING", work_dir=wd
        )
        assert plan["wait"] is True
        assert plan["next_step"] is None
