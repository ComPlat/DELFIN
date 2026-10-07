"""R1 / phase 2 — adversarial tests (s12 reviewer), green against shipped code.

Aggressor pass over the builder's turn-continuation detector
(delfin/agent/turn_continuation.py, commit 5eb9e446 + wiring probe 7cfbe6f3).
Every case here is one where the detector correctly refuses a nudge or catches
it, i.e. where the fix WITHSTOOD the attack. Green adversarial tests pin that
behaviour so a future regression is caught.

Three red adversarial findings (curly U+2019 apostrophe false-negative; signal
inside quote/code false-positive; "I will follow up later" false-positive) were
SENT to the builder and are deliberately NOT here — they are findings, not
committed tests. These cases are the complement.

Commits: 5eb9e446 (detector), 7cfbe6f3 (wiring red control).
"""
from delfin.agent.turn_continuation import announced_action, should_continue


class TestDetectorWithstandsCurlyApostrophes:
    """An ASCII-signal sentence beside a typographic apostrophe is still
    detected through the signal word that carries no apostrophe."""

    def test_let_me_beside_curly_apostrophe(self):
        assert announced_action("Let me fix this. I\u2019ll be right back.") is True


class TestHandbackExclusionsHeld:
    """Hand-back phrases in the builder's list close the turn even when a
    signal word also appears nearby."""

    def test_get_back_handback(self):
        assert announced_action("Let me get back to you on that fix.") is False

    def test_handback_then_announce_is_still_closed(self):
        # "i'll be here" is a hand-back; a later "let me know" stays closed
        assert announced_action("I'll be here. Let me know if needed.") is False

    def test_wait_backchannel_is_closed(self):
        assert announced_action("Waiting for your reply before I act.") is False


class TestDetectorRefusesNonActs:
    """Statements that are not an in-session agent action must never nudge."""

    def test_do_nothing_until_condition(self):
        assert announced_action("Will do nothing until the cluster returns.") is False

    def test_third_person_will(self):
        assert announced_action("The suite will run again.") is False

    def test_question_to_user(self):
        assert announced_action("Should I continue, or stop here?") is False


class TestShouldContinueGuard:
    """The gate: announce AND open work, with a sane note when both hold."""

    # The task count no longer gates the nudge (see
    # tests/test_an_announced_push_is_done_or_explained.py for the two
    # field reports that changed it). What these pin now is that no
    # count -- missing, zero or negative -- raises, and that the
    # announcement alone decides.
    def test_open_tasks_none_still_nudges(self):
        ok, _ = should_continue("Let me fix it.", None)
        assert ok is True

    def test_open_tasks_zero_still_nudges(self):
        ok, _ = should_continue("Let me fix it.", 0)
        assert ok is True

    def test_negative_open_tasks_still_nudges(self):
        ok, _ = should_continue("I will run the suite.", -3)
        assert ok is True

    def test_no_announcement_is_no_nudge_at_any_count(self):
        """The hand-back table carries the whole weight now."""
        for count in (None, 0, -3, 7):
            assert should_continue("Let me know what you want next.",
                                   count)[0] is False, count

    def test_the_note_names_both_conditions(self):
        ok, note = should_continue("Let me fix it.", 1)
        assert ok
        assert "announcing" in note and "not empty" in note
