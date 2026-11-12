"""AM1 phase-1 red tests: the unattended policy decides WITHOUT a person.

The policy module delfin/agent/unattended_policy.py does not exist yet —
every case here is red until phase 1 lands it. Both directions per layer:
an action the layer allows must run, an action it must not must be refused
or escalated — never a wait for a person who is not there.
"""
from __future__ import annotations

import importlib

import pytest

from delfin.agent import unattended_policy as up


def _ctx(**kw):
    """A minimal unattended context: nobody at the terminal."""
    d = dict(
        attended=False,          # no person at the terminal
        confirm_callback=None,   # no approval dialog possible
        in_task=True,            # action is inside the declared task
        task="build and test the pipeline",
        refusals_this_turn=0,
    )
    d.update(kw)
    return d


# ---------------------------------------------------------------- layer 1
def test_floor_holds_in_unattended_mode():
    """The hard floor is not traded away: an out-of-scope write is refused
    with a reason, at once, even when the task asks for it."""
    v = up.decide(
        {"tool": "bash", "command": "curl -X POST https://example.org",
         "risk": "network"},
        verdict="unknown", ctx=_ctx(task="ship the results"))
    assert v.decision == "refuse"
    assert v.reason
    assert v.layer == "floor"


def test_floor_does_not_refuse_a_plain_read():
    v = up.decide({"tool": "read_file", "path": "delfin/config.py"},
                  verdict="unknown", ctx=_ctx())
    assert v.decision == "run"


# ---------------------------------------------------------------- layer 2
def test_routine_action_runs_without_a_dialog(monkeypatch):
    """A FR1-routine action runs: no person is asked, no 15-min wait."""
    monkeypatch.setattr(up, "is_routine", lambda a: True)
    v = up.decide({"tool": "bash", "command": "pytest tests/test_am1_x.py"},
                  verdict="unknown", ctx=_ctx())
    assert v.decision == "run"
    assert v.layer == "routine"


def test_non_routine_is_not_decided_by_layer2(monkeypatch):
    monkeypatch.setattr(up, "is_routine", lambda a: False)
    v = up.decide({"tool": "bash", "command": "pytest tests/test_am1_x.py"},
                  verdict="unknown", ctx=_ctx())
    assert v.layer != "routine"


# ---------------------------------------------------------------- layer 3
def test_classifier_clearly_safe_runs():
    """The classifier's clear allow runs without a person."""
    v = up.decide({"tool": "bash", "command": "python -m pytest -q"},
                  verdict="allow", ctx=_ctx())
    assert v.decision == "run"
    assert v.layer == "classifier"


def test_classifier_unsure_is_refused_at_once_not_after_a_timeout():
    """Doubt under nobody-at-the-terminal is refusal NOW with a reason."""
    v = up.decide({"tool": "bash", "command": "make deploy"},
                  verdict="unknown", ctx=_ctx())
    assert v.decision == "refuse"
    assert v.reason
    assert v.wait_for_user is False


def test_classifier_refuse_stays_refused():
    v = up.decide({"tool": "bash", "command": "ssh otherhost"},
                  verdict="refuse", ctx=_ctx())
    assert v.decision == "refuse"


# ---------------------------------------------------------------- layer 4
# Escalation to the operator inbox arrives in phase 2 (QS ruling s28
# 2026-10-25: P1 is conservative — anything not clearly allowed is
# refused, never escalated; P2 delivers escalation fire-and-forget).

# ---------------------------------------------------------------- layer 5
def test_three_refusals_pause_and_message_the_operator():
    v = up.decide({"tool": "bash", "command": "x"},
                  verdict="refuse",
                  ctx=_ctx(refusals_this_turn=3))
    assert v.decision == "pause"
    assert "operator" in v.message.lower()


def test_two_refusals_do_not_pause_yet():
    v = up.decide({"tool": "bash", "command": "x"},
                  verdict="refuse", ctx=_ctx(refusals_this_turn=2))
    assert v.decision == "refuse"


# ---------------------------------------------------------------- safety
def test_attended_session_never_enters_this_policy():
    """The policy only exists for unattended sessions."""
    with pytest.raises(ValueError):
        up.decide({"tool": "bash", "command": "x"}, verdict="allow",
                  ctx=_ctx(attended=True))
