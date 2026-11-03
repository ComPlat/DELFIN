"""AM1 phase-2 red tests: the unattended 3-refusal lock pauses, not waits.

Driven through the REAL call site (``_classifier_second_opinion``, the
function the gate calls) with no confirm callback: nobody at the
terminal. On the unpatched code the third refusal's text tells the agent
to wait for the user — a wait for nobody (s29 lost two hours this way).
After the px patch it says: paused, queue ONE operator message
fire-and-forget, end the turn, never retry. The wake-shaped path with no
user input is asserted by the turn ending (the refusal text returns; no
dialog opens; nothing waits).
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A


class _NoDialog:
    """There is no person: no dialog object, no callback at all."""
    def ask(self, name, args, preview):  # pragma: no cover - must not run
        raise AssertionError("no dialog may open in unattended mode")


@pytest.fixture
def perms(tmp_path, monkeypatch):
    monkeypatch.setattr(
        A, "_classifier_settings",
        lambda: {"model": "m", "rules": [], "broad": False})
    ws = tmp_path / "wt"
    ws.mkdir()
    p = A.KitToolPermissions(workspace=ws.resolve(), mode="bypassPermissions")
    p.classifier_runner = lambda prompt, model: (
        "DECISION: refuse\nREASON: not the user's task")
    p.classifier_messages = [{"role": "user",
                              "content": "build and test the pipeline"}]
    p._dialog = _NoDialog()
    p.confirm_callback = None          # nobody at the terminal
    return p


def _second_opinion(perms, cmd):
    return A._DocToolExecutor()._classifier_second_opinion(
        perms, "bash", {"command": cmd}, "bypassPermissions")


def _refuse_three(perms, cmds):
    out = []
    for cmd in cmds:
        out.append(_second_opinion(perms, cmd))
    return out


def test_unattended_lock_pauses_with_operator_message_not_a_wait(perms):
    first = _second_opinion(perms, "cmd-one --flag")
    assert first and "refused" in first
    second = _second_opinion(perms, "cmd-two --flag")
    assert second and "refused" in second
    third = _second_opinion(perms, "cmd-three --flag")
    assert third and "refused" in third
    # The fourth call hits the LOCK: unattended -> pause text, no wait.
    locked = _second_opinion(perms, "cmd-four --flag")
    assert locked and "PAUSED" in locked
    assert "operator" in locked
    assert "Do not wait for an answer" in locked
    assert "retry" in locked


def test_attended_lock_keeps_the_user_wait_text(perms):
    """The lock text for an ATTENDED session is unchanged (both
    directions: the patch must not touch the attended path)."""
    perms.confirm_callback = lambda name, args, preview: False
    _second_opinion(perms, "cmd-one --flag")
    _second_opinion(perms, "cmd-two --flag")
    _second_opinion(perms, "cmd-three --flag")
    locked = _second_opinion(perms, "cmd-four --flag")
    assert locked and "three actions were refused in this turn" in locked
    assert "until the user has answered" in locked
    assert "PAUSED" not in locked


def test_pause_decision_is_fire_and_forget_shape(perms):
    """The pause text sends no dialog and carries no fused ask: the call
    RETURNS (the turn can end with nobody at the keyboard)."""
    _second_opinion(perms, "cmd-one --flag")
    _second_opinion(perms, "cmd-two --flag")
    _second_opinion(perms, "cmd-three --flag")
    locked = _second_opinion(perms, "cmd-four --flag")
    # A plain string return — the caller reports it and the turn ends;
    # nothing here blocks on input (the calls above already proved it:
    # they returned instead of hanging).
    assert isinstance(locked, str)
