"""SW2a phase 3 red control: triage dialog actions + pause/resume probes.

Reviewer round-3 scope (nacht-s16): triage Approve/Deny actions that go
through SW1's own triage machinery (the tab must never forge a verdict --
triage_answer re-derives it, triage.py:112-116), and the fail-closed
pause/resume capability probe in BOTH directions (QS ruling nacht-s11).
"""
from pathlib import Path

import pytest

from delfin.dashboard import tab_swarm


TRIAGE_ITEM = {
    "id": "t1",
    "member": "m1",
    "kind": "session_message",
    "text": "please run the test suite",
}


def test_triage_decide_approve_delegates_to_sw1_triage_answer(tmp_path: Path, monkeypatch) -> None:
    """Approve = SW1 triage_answer; SW1 re-derives and refuses escalations."""
    called = {}

    def fake_answer(item, verdict, *, answer):
        called["item"] = item
        called["answer"] = answer
        return {"answered": True, "item_id": item["id"], "logged": "x"}

    monkeypatch.setattr(tab_swarm._triage, "triage_answer", fake_answer)
    out = tab_swarm.triage_decide(TRIAGE_ITEM, "approve", reason="ok")
    assert out["answered"] is True
    assert called["item"]["id"] == "t1"
    assert called["answer"] == "ok"


def test_triage_decide_never_forges_a_verdict(monkeypatch) -> None:
    """The tab must NOT pass its own verdict dict into triage_answer."""
    def fake_answer(item, verdict, *, answer):
        # SW1's real contract: verdict is ignored and re-derived. A tab
        # that pre-computes a verdict dict here would be forging one.
        assert verdict is None
        return {"answered": True, "item_id": item["id"], "logged": "x"}

    monkeypatch.setattr(tab_swarm._triage, "triage_answer", fake_answer)
    tab_swarm.triage_decide(TRIAGE_ITEM, "approve", reason="ok")


def test_triage_decide_surfaces_sw1_refusal(monkeypatch) -> None:
    """An escalation is refused BY SW1 (TriageError), surfaced, not swallowed."""
    from delfin.agent.swarm.triage import TriageError

    def refuse(item, verdict, *, answer):
        raise TriageError("item 't1' is 'escalate' -- goes to the user")

    monkeypatch.setattr(tab_swarm._triage, "triage_answer", refuse)
    out = tab_swarm.triage_decide(TRIAGE_ITEM, "approve", reason="ok")
    assert out["answered"] is False
    assert "escalate" in out["error"]


def test_triage_decide_deny_is_fail_closed() -> None:
    """SW1 ships no deny verb: deny marks the item left for the user."""
    out = tab_swarm.triage_decide(TRIAGE_ITEM, "deny", reason="no")
    assert out["answered"] is False
    assert out["left_for_user"] is True


# ---- pause/resume capability probe: BOTH directions (QS nacht-s11) ----


def test_pause_resume_absent_when_sw1_lacks_the_verbs() -> None:
    """The real SW1 surface today: no pause/resume -> no verbs offered."""
    actions = tab_swarm.available_actions({"phase": "building"})
    assert "pause_all" not in actions
    assert "resume_all" not in actions


def test_pause_resume_present_once_sw1_ships_the_verbs(monkeypatch) -> None:
    """When SW1 grows pause_all/resume_all, the tab exposes them."""
    import delfin.agent.swarm as _swarm_pkg
    monkeypatch.setattr(_swarm_pkg, "pause_all", lambda: None, raising=False)
    monkeypatch.setattr(_swarm_pkg, "resume_all", lambda: None, raising=False)
    actions = tab_swarm.available_actions({"phase": "building"})
    assert "pause_all" in actions
    assert "resume_all" in actions
