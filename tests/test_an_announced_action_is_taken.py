"""A turn that ends by announcing an action is sent back to take it.

From the field (2026-10-08). One session wrote

    Ich korrigiere die dokumentierten Fehler und die Inkonsistenz. Zuerst
    lese ich die betroffenen Zeilen, um sauber zu editieren:

and ended its turn there. Another spent seven consecutive turns with no
tool call at all, each one announcing. The existing auto-continue could
not fire for either: it is gated on OPEN TASKS, and neither session used
the task tools.

An announcement is a stronger signal than an open task -- the model has
named the action itself. So a second trigger: the last sentence announces
an action, the turn is not waiting on a question or a watched job, and it
has not already fired this turn.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A


_THE_REAL_ONE = ("Ich korrigiere die dokumentierten Fehler und die "
                 "Inkonsistenz. Zuerst lese ich die betroffenen Zeilen, "
                 "um sauber zu editieren:")


class TestWhatCountsAsAnAnnouncement:
    @pytest.mark.parametrize("text", [
        _THE_REAL_ONE,
        "Zuerst lese ich die betroffenen Zeilen.",
        "Ich prüfe jetzt die Datei.",
        "Jetzt schreibe ich den CSV-Writer.",
        "Als Nächstes erstelle ich den Test.",
        "Dann korrigiere ich die Zuordnung:",
        "Ich werde jetzt die Datei öffnen.",
        "Let me read the affected lines.",
        "I'll start by reading the file.",
        "Now I will run the tests:",
        "First, I read the blueprint.",
        "Done with the review.\n\nNext, I'll fix the IBM attribution.",
    ])
    def test_an_announcement_at_the_end(self, text):
        assert A._announces_an_action(text), text

    @pytest.mark.parametrize("text", [
        "",
        "Soll ich git push origin main jetzt ausführen?",
        "Let me read the file first — is that what you want?",
        # Narration in the middle of a finished answer is not a promise.
        "Zuerst lese ich die Zeilen. Danach habe ich sie korrigiert und "
        "der Test ist grün.",
        "I'll read the file. Done: three lines changed, tests pass.",
        "Die Zuordnung ist korrigiert; alle 12 Tests laufen.",
        "The review found two errors; both are fixed.",
        "Ergebnis: 31 Belege, Summe 4.512,30 EUR.",
    ])
    def test_a_finished_answer_or_a_question_is_not(self, text):
        assert not A._announces_an_action(text), text

    def test_trailing_markdown_does_not_hide_it(self):
        assert A._announces_an_action("**Zuerst lese ich die Datei:**")
        assert not A._announces_an_action("**Soll ich das tun?**")


class TestItReachesTheTurnLoop:
    def test_the_trigger_is_wired_once_per_turn(self):
        """One source check, for the wiring a behavioural test cannot
        reach without a full fake of the streaming endpoint: the
        predicate decides, it fires at most once per turn, it counts
        against the shared cap, and it never fires over a question --
        which the predicate itself refuses, so a question is not sent
        back in and asked twice (report 20260915-132613)."""
        import inspect

        src = inspect.getsource(A)
        assert '_announces_an_action("".join(_text_chunks))' in src
        assert "not _announced_once" in src
        assert "_announced_once = True" in src
        assert "_auto_cont_count < _AUTO_CONT_CAP\n" in src or \
            "_auto_cont_count < _AUTO_CONT_CAP" in src
        assert "announced an action → taking it" in src

    def test_the_vocabulary_lives_in_the_shared_module(self):
        from delfin.agent import german as G
        assert hasattr(G, "ANNOUNCES_ACTION_RE")
        assert G.ANNOUNCES_ACTION_RE.search("Zuerst lese ich die Zeilen:")
