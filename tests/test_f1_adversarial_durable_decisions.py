"""F1 adversarial battery vs. durable.decide_recovery (builder hash 6d1789304).

Four attacks on the decision engine, all through the public surface with
records produced by the REAL writers:

A. A checker that RAISES must not abort recovery: an unreadable side
   effect is "unknown" -> MUST_ASK, not a crashed restart.
B. A journal record that stores side_effect=False for a tool the
   classifier calls side-effecting (corrupted/torn record) must still
   block an automatic re-run: the classification, not the stored flag,
   is the safety authority for OPEN steps.
C. A stale checkpoint (turn 2) above a journal that already finished
   turn 5 must not send the resume point backwards.
D. would_block_rerun's own contract (docstring): a FINISHED step must
   also be blocked -- the guard is the last line of defence if the loop
   re-dispatches one.

Each case runs red or green against the merged builder hash; reds are
findings for nacht-s1.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_store, step_journal
from delfin.agent.durable import RecoveryVerdict, decide_recovery


@pytest.fixture
def store_root(tmp_path, monkeypatch):
    root = tmp_path / "store"
    root.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", root)
    return root


def test_A_checker_raising_does_not_abort_recovery(store_root, tmp_path):
    """A broken checker degrades to MUST_ASK; recovery still returns a plan."""
    session = "adv-a"
    step_journal.start_step(session, turn=3, step=1, tool="bash",
                            arguments=f"touch {tmp_path / 'x'}", turn_key="t3")
    session_store.save_turn_checkpoint(session, {"task": "x", "turn": 3})

    def broken(_step):
        raise RuntimeError("checker blew up")

    plan = decide_recovery(session, root=store_root, checkers={"bash": broken})
    decision = plan.decision_for_open_step(step=1, turn=3)
    assert decision is not None
    assert decision.verdict is RecoveryVerdict.MUST_ASK
    assert decision.auto_rerun is False


def test_B_stored_side_effect_false_cannot_unlock_a_side_effecting_rerun(
        store_root, tmp_path):
    """A corrupted record's side_effect=False must not allow an auto re-run
    of a tool the classifier knows is side-effecting."""
    session = "adv-b"
    target = tmp_path / "f1-side-effect.sentinel"
    step_journal.start_step(session, turn=7, step=1, tool="bash",
                            arguments=f"touch {target}", turn_key="t7")
    session_store.save_turn_checkpoint(session, {"task": "x", "turn": 7})
    assert not target.exists()

    # Corrupt the record the way a torn/edited journal could: claim read-only.
    open_steps, _ = step_journal.open_steps(session, root=store_root)
    assert len(open_steps) == 1
    # Rewrite the journal line with side_effect explicitly False.
    jp = step_journal.journal_path(session, root=store_root)
    lines = jp.read_text(encoding="utf-8").splitlines(keepends=True)
    fixed = []
    for line in lines:
        rec = line[:-1]
        import json as _json
        obj = _json.loads(rec)
        if obj.get("kind") == "start" and obj.get("tool") == "bash":
            obj["side_effect"] = False
        fixed.append(_json.dumps(obj) + "\n")
    jp.write_text("".join(fixed), encoding="utf-8")

    plan = decide_recovery(session, root=store_root)
    decision = plan.decision_for_open_step(step=1, turn=7)
    assert decision is not None
    assert decision.side_effecting is True
    assert decision.auto_rerun is False
    assert decision.verdict is not RecoveryVerdict.RE_RUN_OK


def test_C_stale_checkpoint_does_not_move_resume_backwards(store_root, tmp_path):
    """A checkpoint from turn 2 above a journal that finished turn 5 must not
    make the plan resume turn 2 (redoing finished work)."""
    session = "adv-c"
    session_store.save_turn_checkpoint(session, {"task": "old", "turn": 2})
    rec = step_journal.auto_start(session, tool="bash",
                                  arguments="git commit -am done", turn=5,
                                  turn_key="t5")
    step_journal.auto_finish(session, rec, ok=True, result="done")

    plan = decide_recovery(session, root=store_root)
    finished = [d for d in plan.decisions if d.verdict is RecoveryVerdict.FINISHED]
    assert finished, "the finished turn-5 step must be in the plan"
    assert plan.resume_turn >= 5, (
        f"resume point went backwards: checkpoint turn 2 vs journal turn 5 "
        f"(resume_turn={plan.resume_turn})"
    )


def test_D_would_block_rerun_blocks_a_finished_step(store_root, tmp_path):
    """would_block_rerun's docstring: a finished step must also be blocked."""
    session = "adv-d"
    rec = step_journal.auto_start(session, tool="bash",
                                  arguments="git push origin HEAD", turn=4,
                                  turn_key="t4")
    step_journal.auto_finish(session, rec, ok=True, result="pushed")

    plan = decide_recovery(session, root=store_root)
    blocked = plan.would_block_rerun(session, tool="bash",
                                     step=rec["step"], turn=rec["turn"])
    assert blocked is True, "a FINISHED step must not be re-dispatchable"
