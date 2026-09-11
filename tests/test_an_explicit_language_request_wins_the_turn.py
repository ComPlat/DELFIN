""""Antworte auf Englisch" is an instruction, not a change of language.

Steering drive, 2026-09-11: a German session, a mid-run steer "nur mod1
und mod2, und antworte auf Englisch", an English answer -- and the
language guard rewrote it into German with "Sie haben recht, die Antwort
muss in Deutsch sein". The session pin is right for a session; a request
in the turn's own text is narrower and more recent, and wins that turn.
"""
import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import verify_guard as vg
from delfin.agent.api_client import StreamEvent


@pytest.mark.parametrize("text,want", [
    ("Änderung: nur mod1 und mod2, und antworte auf Englisch.", "en"),
    ("Bitte antworte auf Englisch: was gibt f1 zurück?", "en"),
    ("Kannst du das bitte auf Deutsch erklären?", "de"),
    ("Please answer in English.", "en"),
    ("Summarize the file in German, please.", "de"),
    ("Auf Englisch bitte.", "en"),
    ("In German please", "de"),
])
def test_an_instruction_names_the_language(text, want):
    assert vg.explicit_language_request(text) == want


@pytest.mark.parametrize("text", [
    "Die Doku ist auf Deutsch geschrieben.",
    "The manual is in English and the code is in Python.",
    "Was macht f1 im Ordner mod1?",
    "",
])
def test_a_statement_is_not_a_request(text):
    assert vg.explicit_language_request(text) == ""


def test_the_last_request_wins():
    assert vg.explicit_language_request(
        "Antworte auf Deutsch. Nein, doch bitte antworte auf Englisch.") == "en"


def test_the_language_correction_says_who_asks():
    text = vg.language_mismatch_feedback("de")
    assert "Automatic check, not a message from the user" in text
    assert "Do not apologise" in text


# --- through the engine ---------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


GERMAN_OPENER = "Kannst du mir bitte die Module in diesem Ordner erklären?"
GERMAN_ANSWER = "Gern, hier ist die Erklärung der Module in diesem Ordner."
ENGLISH_ANSWER = ("f1 returns the value of x unchanged, so nothing else "
                  "happens there and the result is the input itself.")


def _client(replies, hook=None):
    fake = MagicMock()
    fake._observed_files_session = set()
    calls = {"n": 0}

    def _stream(*a, **k):
        i = calls["n"]
        calls["n"] += 1
        if hook is not None:
            hook(i)
        yield StreamEvent(type="text_delta",
                          text=replies[min(i, len(replies) - 1)])
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    return fake


def _engine(agent_tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


def test_a_request_in_the_message_keeps_the_english_answer(agent_tree):
    client = _client([GERMAN_ANSWER, ENGLISH_ANSWER])
    engine = _engine(agent_tree, client)
    engine.stream_response(GERMAN_OPENER)
    assert engine._session_language == "de"
    out = engine.stream_response(
        "Bitte antworte auf Englisch: was gibt f1 zurück?")
    assert client.stream_message.call_count == 2, "no correction turn"
    assert "[Verify]" not in out and "returns" in out
    assert engine._session_language == "de", "the pin does not move"


def test_without_a_request_the_pin_still_corrects(agent_tree):
    client = _client([GERMAN_ANSWER, ENGLISH_ANSWER, GERMAN_ANSWER])
    engine = _engine(agent_tree, client)
    engine.stream_response(GERMAN_OPENER)
    engine.stream_response("Und was gibt die Funktion f1 genau zurück?")
    assert client.stream_message.call_count == 3, "one correction turn"


def test_a_request_in_a_mid_run_steer_counts(agent_tree):
    holder = {}

    def hook(i):
        if i == 1:
            holder["engine"].steer("Nur mod1, und antworte auf Englisch.")

    client = _client([GERMAN_ANSWER, ENGLISH_ANSWER], hook=hook)
    engine = _engine(agent_tree, client)
    holder["engine"] = engine
    engine.stream_response(GERMAN_OPENER)
    out = engine.stream_response("Lies die Module und erkläre jede Funktion.")
    assert client.stream_message.call_count == 2, "no correction turn"
    assert "returns" in out


def test_steers_are_per_turn(agent_tree):
    client = _client([GERMAN_ANSWER])
    engine = _engine(agent_tree, client)
    engine.steer("antworte auf Englisch")
    assert engine._turn_steers == ["antworte auf Englisch"]
    engine.stream_response(GERMAN_OPENER)
    assert engine._turn_steers == []
