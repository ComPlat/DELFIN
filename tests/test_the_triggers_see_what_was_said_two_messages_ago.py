"""The conversation reached the trigger matcher through nobody.

``build_system_prompt`` takes ``conversation_text`` and hands it down
four layers to ``_detect_active_modules``, whose docstring names the
exact case it is for -- "die Tabelle" said once, then "ja, mach das".
Every caller passed ``task_text`` alone. The parameter was threaded end
to end, tested at the bottom, and fed by no one: the mechanism existed
and one hop in the middle dropped it.

What that cost: a follow-up with no trigger words of its own loaded a
prompt with the modules for the subject stripped out. The sticky union
hides it after the first triggering message of a session and does
nothing for the case that starts on a follow-up -- a resumed session,
or a user who names the subject once and then works.

Bounded on purpose: the last few USER turns, capped in characters. The
transcript is unbounded, and a pasted logfile deciding which modules
load is a different failure of the same kind.
"""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# solo mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def _engine(agent_tree):
    from delfin.agent.engine import AgentEngine
    client = MagicMock()
    client.model = "kit.qwen3.5-397b-A17b"
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api", provider="kit",
                          model="kit.qwen3.5-397b-A17b", mode="solo",
                          pack_dir=agent_tree)
    return eng


# ---------------------------------------------------------------------------
# What the engine offers the matcher
# ---------------------------------------------------------------------------

def test_the_last_user_turns_are_offered(agent_tree):
    eng = _engine(agent_tree)
    eng.messages = [
        {"role": "user", "content": "schau dir die Geometrie an"},
        {"role": "assistant", "content": "ok"},
        {"role": "user", "content": "ja, mach das"},
    ]
    text = eng._recent_user_text()
    assert "Geometrie" in text
    assert "ja, mach das" in text


def test_the_assistants_own_words_are_not_offered(agent_tree):
    """The modules shape how the agent phrases things, so matching on
    that phrasing would close a loop with nothing outside it."""
    eng = _engine(agent_tree)
    eng.messages = [
        {"role": "user", "content": "hallo"},
        {"role": "assistant", "content": "Ich starte eine Frequenzrechnung."},
    ]
    assert "Frequenzrechnung" not in eng._recent_user_text()


def test_the_history_offered_is_bounded(agent_tree):
    eng = _engine(agent_tree)
    eng.messages = [{"role": "user", "content": "x" * 5000}
                    for _ in range(40)]
    text = eng._recent_user_text()
    assert len(text) <= eng._TRIGGER_HISTORY_CHARS + 40


def test_a_pasted_logfile_cannot_crowd_out_the_current_line(agent_tree):
    """The most recent turn is the one that must survive the cap."""
    eng = _engine(agent_tree)
    eng.messages = [
        {"role": "user", "content": "noise " * 5000},
        {"role": "user", "content": "und jetzt die Geometrie"},
    ]
    assert "Geometrie" in eng._recent_user_text()


def test_an_image_message_contributes_its_text_and_does_not_raise(agent_tree):
    eng = _engine(agent_tree)
    eng.messages = [{"role": "user", "content": [
        {"type": "text", "text": "was steht in der Tabelle?"},
        {"type": "image", "source": {"data": "..."}},
    ]}]
    assert "Tabelle" in eng._recent_user_text()


def test_an_empty_history_is_an_empty_string(agent_tree):
    eng = _engine(agent_tree)
    eng.messages = []
    assert eng._recent_user_text() == ""


# ---------------------------------------------------------------------------
# It arrives where it was always meant to
# ---------------------------------------------------------------------------

def test_the_prompt_builder_is_given_the_conversation(agent_tree):
    eng = _engine(agent_tree)
    eng.messages = [
        {"role": "user", "content": "schau dir die Geometrie an"},
        {"role": "assistant", "content": "ok"},
    ]
    seen = {}

    def _spy(**kw):
        seen.update(kw)
        return "PROMPT"

    eng.loader.build_system_prompt = _spy
    eng._build_current_system_prompt(task_text="ja, mach das")
    assert "Geometrie" in seen.get("conversation_text", "")


def test_a_follow_up_keeps_the_module_its_subject_needs(agent_tree):
    """End to end through the real matcher, on a fresh session key so the
    sticky union cannot be what makes it pass."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader.__new__(PromptLoader)
    loader._prompt_state = {}

    without = loader._detect_active_modules(
        "ja, mach das", mode_id="solo", session_key="a")
    with_history = loader._detect_active_modules(
        "ja, mach das", mode_id="solo", session_key="b",
        conversation_text="schau dir bitte die Geometrie an")

    assert "chemistry" not in without
    assert "chemistry" in with_history


def test_the_conversation_cannot_activate_everything(agent_tree):
    """Still a trigger match, not a switch that turns the stripping off."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader.__new__(PromptLoader)
    loader._prompt_state = {}
    active = loader._detect_active_modules(
        "ja, mach das", mode_id="solo", session_key="c",
        conversation_text="schau dir bitte die Geometrie an")
    assert active != set(loader._MODULE_TRIGGERS)
