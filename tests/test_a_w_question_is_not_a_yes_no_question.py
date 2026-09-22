"""Yes/No buttons appeared under a greeting, and there was no question.

Reported from the dashboard on 2026-09-22. The session had done nothing
but say hello:

    Hallo! Ich bin Session C, bereit im Worktree … Was soll ich tun?

Under it: two buttons, Yes and No. The user pressed Yes, and the agent
answered that "yes" alone gave it no instruction — which is true, and is
the point. The user's question was the right one: why is there a Yes/No
on a greeting at all?

``detect_question`` decides what to draw. Its yes/no list holds the
German "soll ich" and the English "should i", meant for

    Soll ich die Tests laufen lassen?      Shall I run the tests?

Both substrings sit inside a W-question just as happily:

    Was soll ich tun?                      What should I do?

and with a question mark somewhere in the tail that was enough. A
W-question is OPEN by construction — "was", "wie", "welche", "warum" ask
for a thing, and a thing is not yes or no. So is a question offering
alternatives: "Soll ich A oder B?" has two answers and neither is "yes".

The buttons are not cosmetic. They put a word into the conversation that
the user did not mean, and the next turn is spent recovering from it.
"""

from __future__ import annotations

import pytest

from delfin.dashboard.tab_agent import detect_question as _detect


GREETING = ("Hallo! Ich bin Session C, bereit im Worktree "
            "delfin-wt-591e51b8 (Branch session/591e51b8, Status: clean). "
            "Was soll ich tun?")


# -- the report ------------------------------------------------------------

def test_the_greeting_that_started_this_is_open_not_yes_no():
    got = _detect(GREETING)
    assert got and got["type"] == "open", got


@pytest.mark.parametrize("text", [
    "Was soll ich tun?",
    "Wie soll ich vorgehen?",
    "Welche Datei soll ich zuerst lesen?",
    "Warum soll ich das ändern?",
    "Wann soll ich den Gate starten?",
    "Wo soll ich anfangen?",
    "Wer soll das entscheiden?",
    "What should I do?",
    "Which test should I run first?",
    "How should I proceed?",
])
def test_a_w_question_is_open(text):
    got = _detect("Ich bin bereit. " + text)
    assert got and got["type"] == "open", f"{text!r} -> {got}"


@pytest.mark.parametrize("text", [
    "Soll ich A oder B nehmen?",
    "Should I use the cage or the gate?",
])
def test_a_question_offering_alternatives_is_open(text):
    """Two answers, and neither of them is "yes"."""
    got = _detect("Ich habe zwei Wege gefunden. " + text)
    assert got and got["type"] == "open", f"{text!r} -> {got}"


# -- what must keep working ------------------------------------------------

@pytest.mark.parametrize("text", [
    "Soll ich die Tests laufen lassen?",
    "Shall I run the tests?",
    "Should I commit this?",
    "Möchtest du dass ich den Branch pushe?",
    "Proceed?",
])
def test_a_real_yes_no_question_still_gets_its_buttons(text):
    got = _detect("Der Gate ist grün. " + text)
    assert got and got["type"] == "yesno", f"{text!r} -> {got}"


def test_a_statement_asks_nothing():
    assert _detect("Der Gate ist grün, 17374 Tests sind durch.") is None


def test_numbered_options_are_still_a_choice():
    got = _detect("Ich sehe zwei Wege. Welchen nimmst du?\n"
                  "1) Die Wurzel reparieren\n2) Das Symptom abfangen")
    assert got and got["type"] == "numbered", got
    assert len(got["options"]) == 2


def test_it_never_raises_on_anything():
    for junk in ("", "?", "x" * 5000, "Was soll ich tun? " * 200):
        try:
            _detect(junk)
        except Exception as exc:                        # noqa: BLE001
            pytest.fail(f"{junk[:40]!r} raised {type(exc).__name__}: {exc}")


def test_the_tab_uses_the_one_at_module_level():
    """Pulled out of the builder so it can be driven at all. A second
    copy inside would drift from this one, and the drawing would follow
    the copy nobody tests."""
    import inspect

    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    assert src.count("def detect_question(") == 1
    assert "def _detect_question(" not in src, (
        "the nested copy is still there; the tab would keep calling it")
