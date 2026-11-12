"""A turn that ends on an announced step says that nothing runs.

Bug report davh 414e5a17 (2026-10-20): "he claims he is solving it now, and
then there is no feedback at all -- no progress bar, no sign. The user cannot
tell whether the agent works or the model hangs." The turns had ended on
"Ich speichere zuerst den Feature-Entscheid des Users, dann implementiere
ich:" and "Ich lese die Imports (Z.1-50) ..., um meine Edits anzusetzen:";
the indicator went out -- correctly, the turn was over -- and nothing said so.

Driven through the real AgentEngine.stream_response. Also pinned: the tool
loop sends an announcing answer back up to three times per turn, the second
and third only after tool activity (tests/test_an_announced_action_is_taken.py
pins the wiring).
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


def _client(text):
    class _C:
        def stream_message(self, *_a, **_k):
            yield StreamEvent(type="text_delta", text=text)
            yield StreamEvent(type="message_delta", stop_reason="end_turn",
                              input_tokens=1, output_tokens=1)
    return _C()


def _run(tmp_path, text, open_tasks=0):
    with patch("delfin.agent.engine.create_client",
               return_value=_client(text)):
        eng = E.AgentEngine(repo_dir=str(tmp_path), backend="api",
                            provider="kit", model="kit.deepseek-v4-flash",
                            mode="solo")
    eng._open_task_count = lambda: open_tasks
    shown = []
    eng.stream_response("bau den Plot ein", on_token=shown.append)
    return eng, "".join(shown)


def test_an_announcing_end_without_a_queued_follow_up_says_nothing_runs(
        tmp_path):
    # The line reference ("Z.1-50") hid this ending from both announcement
    # detectors, so the engine queues no follow-up; the colon close is
    # still an announcement, and the user is told nothing runs.
    eng, shown = _run(tmp_path, "Ich lese die Imports (Z.1–50) und den "
                      "_build-Bereich mit main-VBox (Z.160–215), um meine "
                      "Edits anzusetzen:")
    assert not eng.pending_turn_continuation
    assert "Nothing is running now" in shown and "weiter" in shown


@pytest.mark.parametrize("text", [
    "Ich speichere zuerst den Feature-Entscheid des Users, dann "
    "implementiere ich:",
    "Let me read the panel constructor first.",
])
def test_a_queued_follow_up_is_not_called_nothing(tmp_path, text):
    # The engine queues a follow-up the surface runs (up to three in a row);
    # the surface says "nothing runs" when it does not run it.
    eng, shown = _run(tmp_path, text)
    assert eng.pending_turn_continuation
    assert "Nothing is running now" not in shown


@pytest.mark.parametrize("text", [
    "Fertig: Plot eingebaut, Tests grün.",
    "Soll ich den Plot so einbauen?",
    "Ergebnis:\n\n- Breite 10\n- Overlay eingebaut",
])
def test_a_finished_answer_or_a_question_says_nothing_extra(tmp_path, text):
    _, shown = _run(tmp_path, text)
    assert "Nothing is running now" not in shown


# -- the terminal, when its follow-up cap is spent -------------------------

class _Transcript:
    theme = None

    def __init__(self):
        self.lines = []

    def chrome(self, line):
        self.lines.append(line)


def _terminal(done):
    from types import SimpleNamespace
    from delfin.agent.repl import TerminalAgent
    eng = SimpleNamespace(pending_turn_continuation="System note: ...",
                          clear_turn_continuation=lambda: None)
    t = SimpleNamespace(engine=eng, transcript=_Transcript(),
                        _announce_followups=done,
                        _ANNOUNCE_FOLLOWUP_CAP=3,
                        _presence_key=lambda: "")
    t._say_nothing_runs = lambda: TerminalAgent._say_nothing_runs(t)
    return t, TerminalAgent._continue_after_announcement(t)


def test_the_terminal_runs_the_follow_up_below_the_cap():
    t, note = _terminal(done=0)
    assert note and t.transcript.lines == []


def test_the_terminal_says_nothing_runs_when_the_cap_is_spent():
    t, note = _terminal(done=3)
    assert note == ""
    assert any("Nothing is running now" in ln for ln in t.transcript.lines)
