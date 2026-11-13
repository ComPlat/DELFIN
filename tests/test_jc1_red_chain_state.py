"""JC1 phase-1 red control: the chain state machine (builder nacht-s5).

Pins the behaviour JC1 target 1 demands. Red on main 9524f9dde, where
delfin/tools/chain_state.py does not exist. Committed BEFORE the
implementation; turned green by the same commit's implementation.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest


def _module_available(name: str) -> bool:
    try:
        __import__(name)
        return True
    except ImportError:
        return False


def test_chain_state_module_exists():
    try:
        import delfin.tools.chain_state  # noqa: F401
    except ImportError as exc:
        pytest.fail(f"delfin.tools.chain_state missing: {exc}")


@pytest.mark.skipif(
    not _module_available("delfin.tools.chain_state"),
    reason="chain_state not implemented yet",
)
class TestChainStateMachine:
    def _chain(self, tmp_path: Path):
        from delfin.tools import chain_state

        wd = tmp_path / "chain_wd"
        chain = chain_state.ChainState.create(name="demo", work_dir=wd)
        for i in range(3):
            chain.add_step(f"step_{i}", "xtb_opt")
        return chain

    def test_steps_start_pending_and_transition(self, tmp_path):
        chain = self._chain(tmp_path)
        assert [s.state for s in chain.steps] == ["pending"] * 3
        step = chain.steps[0]
        step.mark_submitted(job_id="12345")
        step.mark_running()
        step.mark_finished(
            evidence={"normal_termination": True, "scf_converged": True}
        )
        chain.save()
        reloaded = chain_state.ChainState.load(chain.state_path)
        assert reloaded.steps[0].state == "finished"
        assert reloaded.steps[0].job_id == "12345"
        assert reloaded.steps[0].evidence["normal_termination"] is True

    def test_invalid_transition_refused(self, tmp_path):
        chain = self._chain(tmp_path)
        with pytest.raises(Exception):
            chain.steps[0].mark_finished(
                evidence={"normal_termination": True, "scf_converged": True}
            )

    def test_finished_requires_evidence_never_a_guess(self, tmp_path):
        chain = self._chain(tmp_path)
        chain.steps[0].mark_submitted(job_id="1")
        chain.steps[0].mark_running()
        with pytest.raises(Exception):
            chain.steps[0].mark_finished(evidence=None)
        with pytest.raises(Exception):
            chain.steps[0].mark_finished(evidence={"normal_termination": True})
        with pytest.raises(Exception):
            chain.steps[0].mark_finished(evidence={"scf_converged": True})

    def test_atomic_save_survives_crash_before_rename(self, tmp_path):
        chain = self._chain(tmp_path)
        chain.steps[0].mark_submitted(job_id="7")
        chain.save()
        good = json.loads(chain.state_path.read_text())

        class CrashBeforeRename:
            def replace(self, src, dst):
                raise RuntimeError("simulated crash before rename")

        original_replace = os.replace
        try:
            os.replace = CrashBeforeRename().replace
            with pytest.raises(RuntimeError):
                chain.save()
        finally:
            os.replace = original_replace

        survived = json.loads(chain.state_path.read_text())
        assert survived == good

    def test_no_temp_files_after_clean_save(self, tmp_path):
        chain = self._chain(tmp_path)
        chain.save()
        leftovers = [
            p for p in chain.state_path.parent.iterdir()
            if p.name != chain.state_path.name
        ]
        assert leftovers == []

    def test_load_of_torn_file_never_reports_finished(self, tmp_path):
        chain = self._chain(tmp_path)
        chain.steps[0].mark_submitted(job_id="9")
        chain.save()
        # Tear it: truncate mid-write, as a crash would.
        raw = chain.state_path.read_text()
        chain.state_path.write_text(raw[: len(raw) // 2])
        reloaded = chain_state.ChainState.load(chain.state_path)
        for step in reloaded.steps:
            assert step.state != "finished"

    def test_state_file_next_to_work_dir(self, tmp_path):
        chain = self._chain(tmp_path)
        assert chain.state_path.parent == chain.work_dir

    def test_inputs_hash_recorded(self, tmp_path):
        chain = self._chain(tmp_path)
        h = chain.steps[0].inputs_hash
        assert isinstance(h, str) and h
