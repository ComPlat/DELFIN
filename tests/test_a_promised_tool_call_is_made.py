"""A turn that PROMISES a tool call in prose and makes none is corrected.

From the swarm harness (main d44702503, 216 runs, mostly kit.
deepseek-v4-flash): the final text claims the action is happening NOW
("I am executing the write_file tool call now", "I am now emitting the
actual tool call block", "I'll apply the multi_edit now ...", "I'll run
git status/log ...") and the turn ends with no tool call. PR #600's DSML
recovery covers leaked tool-call MARKUP; prose promises fall through and
the turn ends "done" with the task list untouched.

The existing announced-action trigger (test_an_announced_action_is_taken)
does not catch them either: it examines only the LAST sentence -- shapes 1
and 2 put the promise in the first sentence and declarative filler last --
and its vocabulary lacks the present-continuous "I am executing/emitting"
forms.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import german as G


# The four verbatim field shapes.
SHAPE1 = ("I am executing the write_file tool call now -- creating "
          "`tools/adapters/orca_scan_opi.py` with its complete content. "
          "This is the actual tool invocation being executed.")
SHAPE2 = ("I am now emitting the actual tool call block. No prose -- "
          "just the invocation.")
SHAPE3 = ("I'll apply the `multi_edit` now with the exact baseline from "
          "lines 355-360 ...")
SHAPE4 = ("I'll run git status/log and the two greps: "
          "<system_warning>THIS IS THE ACTION BLOCK -- executing the "
          "batched tools.</system_warning>")


class TestWhatCountsAsAPromisedToolCall:
    @pytest.mark.parametrize("text", [
        SHAPE1,
        SHAPE2,
        SHAPE3,
        SHAPE4,
        # Bare present-continuous forms.
        "I am writing the file now.",
        "I'm applying the patch now.",
        # German: the matcher is the one place German is allowed.
        "Ich führe jetzt den write_file Aufruf aus.",
        "Ich schreibe jetzt die Datei.",
        "Ich erstelle nun den Test.",
    ])
    def test_a_present_tense_promise_matches(self, text):
        assert A._promises_an_unmade_tool_call(text), text

    @pytest.mark.parametrize("text", [
        "",
        # A finished report: the action is in the PAST, not promised.
        "I executed the write_file call and the file exists.",
        "I wrote the file; tests pass.",
        "Ich habe die Datei geschrieben.",
        # A question waits on the user.
        "Soll ich die Datei jetzt schreiben?",
        "Should I apply the edit now?",
        # Plain finished answers without any promise (the control).
        "All three checks pass; nothing else to do.",
        "The energy of the complex is -42.7 kcal/mol at B3LYP/def2-SVP.",
    ])
    def test_a_report_or_question_or_plain_answer_does_not(self, text):
        assert not A._promises_an_unmade_tool_call(text), text

    def test_the_vocabulary_lives_in_the_shared_module(self):
        assert hasattr(G, "PROMISES_TOOL_CALL_RE")
        assert G.PROMISES_TOOL_CALL_RE.search(
            "I am executing the write_file tool call now")


class TestTheCorrectiveTurn:
    """The bounded corrective re-request: one or two corrective user
    turns with the short wording from the brief, then the turn ends with
    a clear (non-done) stop reason so a supervisor resumes it. Never
    after a tool ran in that response."""

    def test_the_corrective_message_wording(self):
        assert A._PROMISED_CONT_TEXT == (
            "No tool call arrived. Make the call now; if the file is "
            "large, write it in parts.")

    def test_budget_and_stop_reason_are_bounded_and_clear(self):
        import inspect
        src = inspect.getsource(A)
        assert "_PROMISED_CONT_MAX = 2" in src
        assert "_promised_count += 1" in src
        assert "_promised_no_tool_call" in src
        assert '_promises_an_unmade_tool_call("".join(_text_chunks))' in src
        # The trigger must stand down once a tool call ran in this
        # response (no double execution).
        assert "not _tool_calls" in src or "not _tool_calls and" in src

    def test_end_state_is_judged_open_or_unknown(self):
        import inspect
        src = inspect.getsource(A)
        # Gated on the task list: open or unknown, never "none".
        assert '_promised_end_state.get("state", "none")' in src


class TestNoDoubleExecution:
    def test_a_tool_round_ignores_the_trigger(self):
        """When tool calls ran in this response the recovery must not
        fire: the matcher is text-only, so the no-double-execution
        guarantee lives in the wiring guard, asserted in
        TestTheCorrectiveTurn and verified behaviourally against the
        fake endpoint (see test_a_promise_repeated_over_a_cut_call)."""
        assert A._promises_an_unmade_tool_call(SHAPE3)
