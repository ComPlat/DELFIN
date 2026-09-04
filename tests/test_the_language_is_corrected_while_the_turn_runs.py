"""The wrong language is said mid-turn, because the end is too late.

Field report 20260904-083658, filed with the question "why still English
and German?". Its roles, in order:

    user · thinking · assistant · tool · tool · thinking · tool · tool
         · thinking · tool · tool · system

One assistant message — the opening sentence before the first tool call —
in German, to an English question. The end-of-turn guard corrects the
ANSWER, and this turn had not produced one yet. So the only text the user
had seen was the one piece of text no guard could reach.

The scanner was never the problem: fed those two strings it returns "en".
What was missing is a place to say so while the turn is still running.
The rail exists — `push_run_note`, drained between tool rounds — and this
is that rail carrying the language.
"""

from __future__ import annotations

from delfin.agent import verify_guard as vg
from delfin.agent.engine import AgentEngine


# Verbatim from the report, minus the trailing quote mark the user's
# paste carried.
ASKED_EN = ("Extract the hyperpolarizability tensor from all folders in "
            "`/pfs/data6/home/ka/ka_ibcs/ka_ew7404/archive` and calculate "
            "the Hyper-Rayleigh scattering hyperpolarizability BHRS from it.")
SAID_DE = ("Ich extrahiere die Hyperpolarisierbarkeit-Tensoren aus allen "
           "Ordnern im Archive und berechne das. Zuerst verschaffe ich mir "
           "einen Ueberblick ueber die Struktur der Ordner.")
SAID_EN = ("I will extract the hyperpolarizability tensors from every "
           "folder in the archive and compute the value. First let me get "
           "an overview of how the folders are laid out.")


class _Client:
    """Just the one method the nudge uses."""

    def __init__(self):
        self.notes: list[str] = []

    def push_run_note(self, text):
        self.notes.append(text)


def _engine(asked=ASKED_EN):
    engine = AgentEngine.__new__(AgentEngine)      # no backend, no network
    engine.client = _Client()
    engine._last_user_message = asked
    engine._language_note_sent = False
    return engine


def test_the_scanner_had_always_seen_it():
    """Stated first, because it is what makes the rest a wiring defect
    rather than a detection one."""
    assert vg.scan_for_language_mismatch(SAID_DE, ASKED_EN) == "en"


def test_a_german_opening_to_an_english_question_is_told_mid_turn():
    engine = _engine()
    engine._nudge_language_if_wrong([SAID_DE])
    assert len(engine.client.notes) == 1
    assert "English" in engine.client.notes[0]
    # It is advice to the model, not a verdict shown to the user.
    assert not engine.client.notes[0].startswith("[Verify]")


def test_an_answer_in_the_right_language_is_left_alone():
    engine = _engine()
    engine._nudge_language_if_wrong([SAID_EN])
    assert engine.client.notes == []


def test_german_in_german_out_is_left_alone():
    """The rule is the user's language, not English."""
    engine = _engine(asked=(
        "Extrahiere den Hyperpolarisierbarkeits-Tensor aus allen Ordnern "
        "im Archiv und berechne daraus den Wert fuer jede Rechnung."))
    engine._nudge_language_if_wrong([SAID_DE])
    assert engine.client.notes == []


def test_it_is_said_once_and_not_on_every_round():
    """A model told the same thing at every tool call has been given no
    chance to act on the first telling."""
    engine = _engine()
    for _round in range(4):
        engine._nudge_language_if_wrong([SAID_DE])
    assert len(engine.client.notes) == 1


def test_a_correction_turn_does_not_demand_english_of_a_german_session():
    """A `[Verify]` body is the harness talking, always in English.
    Reading one back as the user's language would flip every German
    session to English at the first correction."""
    engine = _engine(asked="[Verify] The user wrote their message in German")
    engine._nudge_language_if_wrong([SAID_DE])
    assert engine.client.notes == []


def test_an_opening_too_short_to_judge_says_nothing():
    engine = _engine()
    engine._nudge_language_if_wrong(["Ok."])
    assert engine.client.notes == []


def test_a_client_that_cannot_take_notes_does_not_cost_the_turn():
    """It is advice, not a gate."""
    engine = _engine()
    engine.client = object()                        # no push_run_note
    engine._nudge_language_if_wrong([SAID_DE])      # must not raise


def test_an_engine_that_never_set_the_flag_still_speaks_up():
    """The once-per-turn latch must not be able to fail CLOSED.

    An engine reaching the nudge without the attribute — an older session
    restored, a subagent built by a path that skips the turn setup —
    would otherwise be permanently silent, which is the failure mode that
    looks exactly like the defect being fixed.
    """
    engine = _engine()
    del engine._language_note_sent
    engine._nudge_language_if_wrong([SAID_DE])
    assert len(engine.client.notes) == 1


def test_clearing_the_latch_lets_the_next_turn_say_it_again():
    """What the turn setup does at the top of every turn. A session that
    said it once and never again would leave the second wrong answer
    standing."""
    engine = _engine()
    engine._nudge_language_if_wrong([SAID_DE])
    engine._language_note_sent = False              # as stream_response does
    engine._nudge_language_if_wrong([SAID_DE])
    assert len(engine.client.notes) == 2
