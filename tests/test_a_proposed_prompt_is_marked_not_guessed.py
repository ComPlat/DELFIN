"""A proposed message is read from an explicit marker, not from prose.

Input: the model's answer text. Output: the message to place in the
dashboard input box, or None.

    PROMPT: Hallo, wie geht es dir?

Marker semantics:
  - anchored at line start; an occurrence inside a sentence is ignored
  - the last marker in an answer wins
  - an empty or whitespace-only remainder yields None

Precedence: a marked prompt overrides the open-task list, since it is a
statement by the agent rather than an inference from stored state.

Alternative considered and rejected: extracting the proposal from prose.
That requires selecting one sentence of an answer as the offer. Two
defects with the same shape were measured in this repository on
2026-09-22 (substring match instead of question type; file mode instead
of path reachability), and the failure mode here writes into an input
field the user is about to send.

Presentation: the placeholder is coloured while a proposal is pending and
grey otherwise; Tab copies it into the field without sending.
"""

from __future__ import annotations

import pytest

from delfin.dashboard.tab_agent import (input_placeholder, proposed_prompt)


HINT = "Message the agent... (Enter = send)"


# -- reading the marker -----------------------------------------------------

def test_the_marked_line_is_the_proposal():
    assert proposed_prompt(
        "Klar! Hier wäre ein Vorschlag:\n"
        "PROMPT: Hallo, wie geht es dir?") == "Hallo, wie geht es dir?"


def test_the_rest_of_the_answer_is_not_taken():
    got = proposed_prompt("Ich habe den Gate gefahren, 17813 Tests sind "
                          "durch.\nPROMPT: bring die Kette auf main\n"
                          "Sag Bescheid, wenn etwas fehlt.")
    assert got == "bring die Kette auf main"


def test_an_answer_without_the_marker_proposes_nothing():
    assert proposed_prompt("Hier wäre ein Vorschlag: Hallo, wie geht es "
                           "dir? Das könntest du senden.") is None


def test_the_last_marker_wins():
    """A long answer may reconsider; the final word is the offer."""
    assert proposed_prompt("PROMPT: erst dies\nPROMPT: dann doch das") == (
        "dann doch das")


def test_a_marker_without_text_is_no_proposal():
    assert proposed_prompt("PROMPT:   ") is None
    assert proposed_prompt("PROMPT:") is None


def test_a_marker_inside_a_sentence_is_not_a_marker():
    """It has to open the line, or every answer mentioning the word
    starts putting things in the user's box."""
    assert proposed_prompt(
        "Den Marker PROMPT: schreibt das Modell an den Zeilenanfang.") is None


def test_it_never_raises():
    for junk in ("", None, "PROMPT:" * 500, "x" * 20000):
        proposed_prompt(junk or "")


# -- who wins ---------------------------------------------------------------

def test_a_marked_prompt_outranks_the_open_tasks():
    """A task is what is open; a marked prompt is what the agent just
    decided to offer."""
    got = input_placeholder(["run the full suite"], HINT,
                            proposed="Hallo, wie geht es dir?")
    assert got == "Hallo, wie geht es dir?"


def test_without_a_marker_the_tasks_still_fill_the_line():
    assert input_placeholder(["run the full suite"], HINT) == (
        "run the full suite")


def test_neither_brings_the_hint_back():
    assert input_placeholder([], HINT, proposed=None) == HINT
    assert input_placeholder([], HINT, proposed="  ") == HINT


# -- green, and only when somebody put it there -----------------------------

def test_the_box_is_marked_green_only_while_a_proposal_waits():
    import inspect

    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    assert "delfin-agent-input-proposed" in src, (
        "the class that colours the box has to exist")
    css = src.split("<style>")[1].split("</style>")[0]
    assert "delfin-agent-input-proposed" in css, "and be styled"
    assert "remove_class" in src, (
        "it has to come off again, or the box stays green over a hint")
    # The colour is bound to the placeholder being VISIBLE, not merely to
    # a proposal being pending: one typed character hides the placeholder,
    # and the frame must stop being green around text the user wrote.
    assert "placeholder-shown::placeholder" in css, (
        "typing would leave the colour behind on the user's own text")
    rule = css.split(".delfin-agent-input-proposed")[1].split("}")[0]
    assert "border" not in rule and "background" not in rule, (
        "the text carries the difference; restyling the whole field for a "
        "line that vanishes on the first keystroke is louder than the "
        "information in it")


def test_tab_still_only_fills_and_never_sends():
    import inspect

    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    branch = src.split("e.key === 'Tab'")[1].split("if (e.key ===")[0]
    assert "placeholder" in branch
    assert "send-row" not in branch


# -- and the model is told, in the prompt it actually reads -----------------

def test_the_marker_is_in_the_role_prompt():
    """A rule is a rule when it is IN the built prompt."""
    import pathlib

    from delfin.agent import prompt_loader as PL

    pack = pathlib.Path(PL.__file__).parent / "pack" / "agents"
    text = (pack / "solo_agent.md").read_text(encoding="utf-8")
    assert "PROMPT:" in text, (
        "the model cannot use a marker nobody told it about")
