"""Pressing Enter must not erase the user's side of the conversation.

The raw input box is deliberately cleared on submit.  Until the submitted
text was copied into the transcript, a slow first turn left only a cache
notice and a spinner on screen; there was no way to verify which task was
actually running.  The same happened at the idle prompt after every answer.
"""

from __future__ import annotations

import io

from delfin.agent import repl as R


class _Tty(io.StringIO):
    def isatty(self):
        return True


class _Engine:
    session_id = "session-1"
    messages = []
    kit_permissions = None


def _quiet_lifecycle(agent):
    """Keep this loop test away from process signals and user history."""
    for name in (
        "_install_sigint", "_load_history", "_install_completer",
        "_announce_presence", "_save_history", "_withdraw_presence",
        "_restore_sigint", "_restore_sigwinch",
    ):
        setattr(agent, name, lambda: None)


def test_the_first_prompt_is_visible_before_its_turn_starts():
    err = _Tty()
    agent = R.TerminalAgent(
        _Engine(), out=io.StringIO(), err=err,
        read_line=lambda _prompt: "/exit",
        opts=R.ReplOptions(color="never"),
    )
    _quiet_lifecycle(agent)
    seen = []

    def _turn(prompt):
        seen.append((prompt, err.getvalue()))
        return R.TurnResult()

    agent.turn = _turn
    assert agent.run("meine Nachricht bleibt sichtbar") == 0
    assert seen and seen[0][0] == "meine Nachricht bleibt sichtbar"
    assert "» meine Nachricht bleibt sichtbar" in seen[0][1], (
        "the permanent copy must exist before a cache wait can begin")


def test_a_message_queued_during_a_turn_is_shown_with_its_text():
    from delfin.agent import repl_keys as rk

    err = _Tty()
    agent = R.TerminalAgent(
        _Engine(), out=io.StringIO(), err=err,
        opts=R.ReplOptions(color="never"),
    )
    agent._on_key(
        rk.KeyEvent(rk.SUBMIT, text="und auf main dann bringen"),
        rk.KeyDecoder(),
    )
    shown = err.getvalue()
    assert "» und auf main dann bringen" in shown
    assert "queued" in shown
    assert agent.queued == ["und auf main dann bringen"]


def test_pasted_control_sequences_cannot_repaint_the_transcript():
    err = _Tty()
    agent = R.TerminalAgent(
        _Engine(), out=io.StringIO(), err=err,
        opts=R.ReplOptions(color="never"),
    )
    agent._show_user_input("vorher\x1b[2Jnachher")
    assert "vorhernachher" in err.getvalue()
    assert "\x1b[2J" not in err.getvalue()
