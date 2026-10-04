"""R1 / phase 2 — turn continuation: detector tests.

Contract (package R1, finding 1): when a turn ends with announced-but-not-
done intent ("Let me …", "Next I will …", "I'll now …") AND open work remains,
the session should be nudged to continue — at most once. This file tests the
pure detector in :mod:`delfin.agent.turn_continuation`, NOT the engine wiring
(that is a .gate patch for the operator).

The control: turn_continuation.py does not exist on a2ed0773, so every
import here is red before the module ships.
"""

from delfin.agent.turn_continuation import (
    announced_action,
    should_continue,
    CONTINUATION_LIMIT,
)


class TestAnnouncedActionPositive:
    """Answers that say the agent will do more work in this session."""

    def test_let_me_future_action(self):
        assert announced_action("Let me write the red test first.") is True

    def test_next_i_will(self):
        assert announced_action("Next I will fix the failing test.") is True

    def test_i_will_run_and_report(self):
        assert announced_action("I will run the suite and report.") is True

    def test_ill_now_create(self):
        assert announced_action("I'll now create the plan file.") is True

    def test_im_going_to_implement(self):
        assert announced_action("I'm going to implement the detector.") is True

    def test_let_me_grep_first(self):
        assert announced_action("Let me grep the file first.") is True


class TestAnnouncedActionNegative:
    """Answers that hand back to the user, wait, or finish — never a nudge."""

    def test_let_me_know_is_a_handback(self):
        assert announced_action("Let me know if you want any changes.") is False

    def test_let_me_ask_is_a_question(self):
        assert announced_action("Let me ask you a question.") is False

    def test_i_will_wait(self):
        assert announced_action("I will wait for your answer.") is False

    def test_ill_be_here(self):
        assert announced_action("I'll be here if you need me.") is False

    def test_im_done(self):
        assert announced_action("I'm done here.") is False

    def test_i_will_stop(self):
        assert announced_action("I will stop here.") is False


class TestAnnouncedActionEdges:
    """Boundary and casing."""

    def test_empty_is_false(self):
        assert announced_action("") is False

    def test_whitespace_is_false(self):
        assert announced_action("   \n") is False

    def test_case_insensitive(self):
        assert announced_action("I'LL fix it.") is True

    def test_multiline(self):
        text = "The suite is green.\nNext I will commit the fix."
        assert announced_action(text) is True
class TestShouldContinue:
    """The gate: announced intent AND open work must BOTH hold."""

    def test_announced_and_open_is_true(self):
        ok, note = should_continue("Let me fix the test.", 2)
        assert ok is True
        assert isinstance(note, str) and note.strip()

    def test_no_announcement_never_continues(self):
        ok, _ = should_continue("I'm done here.", 2)
        assert ok is False

    def test_zero_open_tasks_never_continues(self):
        ok, _ = should_continue("Let me fix the test.", 0)
        assert ok is False

    def test_negative_open_is_treated_as_zero(self):
        ok, _ = should_continue("I will run the suite.", -3)
        assert ok is False

    def test_note_present_only_when_continuing(self):
        ok, note = should_continue("Let me fix the test.", 1)
        assert ok and note
        ok2, note2 = should_continue("I'm done.", 1)
        assert ok2 is False and note2 == ""


class TestContinuationLimit:
    """At most ONE follow-up turn per wake cycle — the loop guard."""

    def test_limit_is_one(self):
        assert CONTINUATION_LIMIT == 1
class TestAdversarialBoundary:
    """The cases a strict reviewer probes: partial words and near-hand-backs."""

    def test_third_person_will_is_not_signal(self):
        assert announced_action("The suite will run again.") is False

    def test_the_word_ill_is_not_a_signal(self):
        assert announced_action("Ill never write that file.") is False

    def test_ill_as_noun_phrase(self):
        assert announced_action("I ill not go further.") is False

    def test_let_me_stuck_inside_word(self):
        assert announced_action("theylet me down.") is False

    def test_handback_plus_announced_word_stays_closed(self):
        assert announced_action("I'll be here. Let me know if needed.") is False

    def test_question_to_user_is_not_signal(self):
        assert announced_action("Should I continue, or stop here?") is False

    def test_wait_on_you_is_not_signal(self):
        assert announced_action("Waiting for your reply before I act.") is False

    def test_im_going_to_with_space(self):
        assert announced_action("I am going to run the full suite.") is True

    def test_ammended_announcement_with_handback(self):
        assert announced_action("I will wait for the cluster.") is False


class TestShouldContinueNote:
    def test_note_names_the_two_conditions(self):
        ok, note = should_continue("Let me fix it.", 1)
        assert ok and "announcing" in note and "not empty" in note

    def test_returns_tuple_shape_always(self):
        assert isinstance(should_continue("x", 0), tuple)
        assert isinstance(should_continue("x", 5), tuple)
class TestTypographicApostrophe:
    """U+2019 right single quote — review finding: real transcripts use it."""

    def test_curly_apostrophe_ill(self):
        # "I\u2019ll fix the failing test now." — RED before the fix.
        assert announced_action("I\u2019ll fix the failing test now.") is True

    def test_curly_apostrophe_im(self):
        # "I\u2019m going to run the suite." — RED before the fix.
        assert announced_action("I\u2019m going to run the suite.") is True

    def test_curly_apostrophe_let_me(self):
        assert announced_action("Let me check the value first.") is True

    def test_curly_apostrophe_handback_still_suppresses(self):
        # "I\u2019ll wait for your reply." must stay a hand-back.
        assert announced_action("I\u2019ll wait for your reply.") is False


class TestReviewerDecisions:
    """Cases the reviewer flagged; decisions documented here."""

    def test_i_will_follow_up_later_is_not_an_action(self):
        # Back-channel, not an in-session action — accepted as a hand-back.
        assert announced_action("I will follow up later.") is False

    def test_signal_inside_quoted_instruction_is_accepted(self):
        # ACCEPTED tradeoff (documented): stripping code blocks / quoted
        # spans before scanning is fragile and heavy; a quoted "I will …"
        # inside an answer is rare, while the detector's job is to catch a
        # REAL first-person announcement. Pin True so the behaviour is
        # explicit rather than accidental.
        assert announced_action(
            'You asked me to "I will book the room".') is True

    def test_signal_inside_code_comment_is_accepted(self):
        # ACCEPTED tradeoff: same as above.
        assert announced_action("# let me check the value\n") is True
