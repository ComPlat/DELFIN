"""Adversarial tests for JC1 phases 1-2 (reviewer nacht-s7).

Attacks the resume-decides-from-facts mechanism where it can silently do
the wrong thing.  All cases drive the real delfin.tools.chain_state API
on real on-disk state files -- no shadow copies, no .gate copies.

Findings under test (each red against agent/s5-jc1b @ 1ba4744e9):
  F-A  a ``submitted`` step is invisible to build_resume_plan: the plan
       walks past it and would run the NEXT step while the first job is
       still alive in the scheduler.
  F-B  a finished NON-ORCA step (xtb_opt, the live chain's step 1 per
       target 6) can never re-verify: the output carries neither an ORCA
       normal-termination marker nor "SCF CONVERGED", so a verified
       finished step is downgraded to pending and rerun -- exactly the
       rerun of finished work target 1 forbids.
  F-C  mark_finished accepts falsy evidence: only key PRESENCE is
       checked, so {"normal_termination": False, "scf_converged": False}
       marks a step finished with proof of the opposite.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.tools.chain_state import (
    ChainState,
    ChainStateError,
    build_resume_plan,
    verify_finished_step,
)


# ---------------------------------------------------------------------------
# Helpers: real on-disk chains
# ---------------------------------------------------------------------------

def _chain(tmp_path: Path, n_steps: int = 2) -> ChainState:
    chain = ChainState.create(name="test-chain", work_dir=tmp_path / "chain")
    for i in range(1, n_steps + 1):
        chain.add_step(f"step{i}", "orca_sp")
    chain.save()
    return chain


def _orca_out(dir_path: Path, name: str = "out") -> Path:
    """A minimal real-shaped ORCA output: normal termination + SCF."""
    p = dir_path / f"{name}.out"
    p.write_text(
        "*****************\n* ORCA lord *\n*****************\n"
        "SCF CONVERGED AFTER  12 ITERATIONS\n"
        "TOTAL ENERGY      -1.00000000 Eh\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )
    return p


# ---------------------------------------------------------------------------
# F-A: submitted steps must block the plan, not be walked past
# ---------------------------------------------------------------------------

def test_submitted_step_alive_blocks_the_plan(tmp_path):
    """A step in `submitted` whose job is alive must yield wait=False but
    next_step=None (the caller must watch the job, not start the next
    step).  Walking past a submitted step risks two jobs on one chain."""
    chain = _chain(tmp_path)
    s1, s2 = chain.steps
    s1.mark_submitted(job_id="12345")
    s2.mark_submitted(job_id="12346")
    chain.save()
    plan = build_resume_plan(chain, job_state=lambda jid: "R")
    assert plan.get("wait") is True or plan.get("next_step") is None, (
        f"plan walked past live submitted steps: {plan}"
    )


def test_submitted_step_gone_without_output_is_not_skipped(tmp_path):
    """A submitted step whose job vanished without output must NOT be
    skipped: the plan must name it as the next step (resubmit) rather
    than continuing to the next step."""
    chain = _chain(tmp_path)
    s1, s2 = chain.steps
    s1.mark_submitted(job_id="12345")
    chain.save()
    plan = build_resume_plan(chain, job_state=lambda jid: None)
    assert plan.get("next_step") is s1, (
        f"plan skipped a submitted step whose job vanished: {plan}"
    )


# ---------------------------------------------------------------------------
# F-B: non-ORCA finished steps (xtb) must re-verify, not rerun
# ---------------------------------------------------------------------------

def test_finished_xtb_step_without_orca_markers_stays_skipped(tmp_path):
    """Target 6's chain starts with xtb_opt.  A finished xtb step's
    output has no ORCA termination marker and no 'SCF CONVERGED' line;
    re-verification must still recognise it as finished (its own tool's
    success evidence), not downgrade it and rerun finished work."""
    chain = _chain(tmp_path)
    s1, _ = chain.steps
    s1.tool = "xtb_opt"
    s1.state = "submitted"
    s1.mark_running()
    s1.mark_finished(evidence={"normal_termination": True,
                               "scf_converged": True})
    wd = tmp_path / "chain" / "step1"
    wd.mkdir()
    # An xtb-style output: no ORCA markers at all.
    (wd / "xtb.out").write_text("xtb4stda finished\nTOTAL ENERGY -5.0 Eh\n",
                                encoding="utf-8")
    chain.save()
    verdict = verify_finished_step(s1, chain_work_dir=tmp_path / "chain")
    assert verdict.verified, (
        f"finished xtb step downgraded to {s1.state}: {verdict.reason}"
    )
    assert s1.state == "finished"


# ---------------------------------------------------------------------------
# F-C: falsy evidence must not mark a step finished
# ---------------------------------------------------------------------------

def test_mark_finished_rejects_falsy_evidence(tmp_path):
    chain = _chain(tmp_path)
    s1 = chain.steps[0]
    s1.mark_submitted(job_id="1")
    s1.mark_running()
    with pytest.raises(ChainStateError):
        s1.mark_finished(evidence={"normal_termination": False,
                                   "scf_converged": False})
    assert s1.state != "finished"


def test_mark_finished_rejects_missing_truthiness_per_key(tmp_path):
    """Each required evidence key must be truthy -- an empty/false value
    is a guess, not proof."""
    chain = _chain(tmp_path)
    s1 = chain.steps[0]
    s1.mark_submitted(job_id="1")
    s1.mark_running()
    for bad in (
        {"normal_termination": True, "scf_converged": False},
        {"normal_termination": "", "scf_converged": True},
        {"normal_termination": None, "scf_converged": True},
    ):
        with pytest.raises(ChainStateError):
            s1.mark_finished(evidence=bad)
    assert s1.state != "finished"
