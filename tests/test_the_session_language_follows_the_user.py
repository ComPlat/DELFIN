"""The answer language follows the user's latest message that says one.

Input: the messages of one session, in order. Field report 2026-10-22
(KIT cluster, DeepSeek flash): the session's first long message was an
English task, the user then asked "was kannst du?" in German, and the
answer came in German AND English. The session language had been pinned
by the first message and "a later message in another language does not
change it", so the model held a rule saying English and a question in
German, and wrote both.

The user decided (2026-10-22): follow the language of the user's latest
message that says one. So:
- a message the user typed that the detector can place (three function
  words of one language, twice the other's) sets the language, the first
  one and every later one;
- a message that cannot be placed ("ok", "hi", "ja mach das") keeps it;
- a turn the system starts (a finished job, a wake-up, a message from
  another session: HarnessText) keeps it, so a long unattended run does
  not drift;
- the line in front of the model says this, and that one answer is in one
  language.
"""
from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.harness_rows import HarnessText


@pytest.fixture
def engine(tmp_path):
    from delfin.agent.engine import AgentEngine
    modes = tmp_path / "pack_lite" / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (tmp_path / "pack_lite" / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    fake = MagicMock()
    fake._observed_files_session = set()
    with patch("delfin.agent.engine.create_client", return_value=fake):
        return AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick",
                           pack_dir=tmp_path)


def test_the_latest_message_that_says_one_sets_the_language(engine):
    engine._note_session_language(
        "Which conformation of butane is more stable, anti or gauche?")
    assert engine._session_language == "en"
    engine._note_session_language("was kannst du?")
    assert engine._session_language == "de"
    engine._note_session_language("what can you do?")
    assert engine._session_language == "en"


@pytest.mark.parametrize("short", ["ok", "hi", "ja mach das", "ok danke",
                                   "Thanks, go on"])
def test_a_message_that_cannot_be_placed_keeps_it(engine, short):
    engine._note_session_language("Kannst du mir das bitte erklären?")
    engine._note_session_language(short)
    assert engine._session_language == "de"


def test_a_turn_the_system_starts_keeps_it(engine):
    engine._note_session_language("Kannst du mir das bitte erklären?")
    engine._note_session_language(HarnessText(
        "[watch — a system event, not the user] A job you were watching "
        "has finished: the shell job is done and the output is in the file"))
    assert engine._session_language == "de"


def test_the_line_in_front_of_the_model_says_it_follows_the_user(engine):
    engine._note_session_language("Kannst du mir das bitte erklären?")
    block = engine._session_language_block()
    assert "German" in block
    assert "does not change it" not in block
    assert "latest message" in block
    assert "one language" in block


def test_the_system_prompt_says_the_same(engine):
    from delfin.agent.prompt_loader import PromptLoader
    anchor = PromptLoader._build_critical_anchor("solo_agent", "de")
    assert "set by its first message" not in anchor
    assert "latest message" in anchor
