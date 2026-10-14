"""Per-model tuning for the two KIT flagships, each knob switchable.

Field data (shared bug archive, 33 sessions, 2026-09-16 .. 10-14) and the
benchmark store (950 DeepSeek attempts, 273 GLM attempts):

* kit.deepseek-v4-flash decodes at about 6 s per 100 output tokens
  (request p50, archive) and spends 81% of its turn time generating. Every
  one of its 38 archived tool-call messages carried prose beside the call,
  p50 593 characters ("Let me read ... and then look at ..."), against
  259 of 712 on kit.glm-5.3 at p50 160.

Each knob is a ModelProfile field, so ``agent.model_overrides`` switches
it per model without editing Python, and a benchmark arm that does so is
stamped with the override in its profile name.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import model_profiles as mp
from delfin.agent.prompt_loader import PromptLoader


@pytest.fixture
def settings(tmp_path, monkeypatch):
    """A settings file this test owns; the user's is never read."""
    path = tmp_path / "settings.json"
    path.write_text(json.dumps({"agent": {}}))

    def _write(overrides):
        path.write_text(json.dumps({"agent": {"model_overrides": overrides}}))
        mp._OVERRIDES_CACHE = None

    monkeypatch.setattr("delfin.user_settings.get_settings_path",
                        lambda *a, **k: path)
    mp._OVERRIDES_CACHE = None
    yield _write
    mp._OVERRIDES_CACHE = None


_DEEPSEEK = "kit.deepseek-v4-flash"
_GLM = "kit.glm-5.3-flash"


def _solo_sections(model: str):
    return PromptLoader().compose_sections(
        role_id="solo_agent", mode_id="solo", task_text="hallo",
        session_key=f"s-{model}", model=model)


def _solo_prompt(model: str) -> str:
    return "\n\n".join(s.content for s in _solo_sections(model))


# --- terse tool rounds ------------------------------------------------------

def test_deepseek_is_told_to_keep_tool_rounds_terse(settings):
    assert mp.get_profile(_DEEPSEEK).terse_tool_rounds is True
    note = mp.MODEL_NOTES["terse_tool_rounds"]
    assert note in _solo_prompt(_DEEPSEEK)


def test_glm_is_not(settings):
    """GLM-flash decodes ~7x faster (0.9 s per 100 tokens) and narrates in
    36% of its tool rounds; the note would cost it prompt and buy little."""
    assert mp.get_profile(_GLM).terse_tool_rounds is False
    assert mp.MODEL_NOTES["terse_tool_rounds"] not in _solo_prompt(_GLM)


def test_the_note_can_be_switched_off_per_model(settings):
    settings({"kit.deepseek": {"terse_tool_rounds": False}})
    assert mp.get_profile(_DEEPSEEK).terse_tool_rounds is False
    assert mp.MODEL_NOTES["terse_tool_rounds"] not in _solo_prompt(_DEEPSEEK)


def test_the_note_sits_in_the_stable_head(settings):
    """It is fixed for the model, so it belongs before every per-turn
    section: a note after the volatile material would be re-sent cold
    whenever the git status line moved."""
    sections = _solo_sections(_DEEPSEEK)
    names = [s.name for s in sections]
    assert "model_notes" in names
    note = sections[names.index("model_notes")]
    assert note.layer == PromptLoader.LAYER_STABLE
    layers = [s.layer for s in sections]
    assert layers == sorted(layers)


# --- concise final answers --------------------------------------------------

def test_deepseek_is_told_to_answer_concisely(settings):
    """Benchmark, 950 DeepSeek attempts: within one task, answers longer
    than that task's median passed 0.79 (n=406) and shorter ones 0.75
    (n=399) -- no gain -- at a duration p50 of 55 s against 40 s."""
    assert mp.get_profile(_DEEPSEEK).concise_answers is True
    assert mp.MODEL_NOTES["concise_answers"] in _solo_prompt(_DEEPSEEK)


def test_glm_answers_are_left_as_they_are(settings):
    """GLM's passing answers are already short (p50 710 chars against
    DeepSeek's 2146), and its within-task longer answers passed less
    often, not more -- nothing for the note to buy."""
    assert mp.get_profile(_GLM).concise_answers is False
    assert mp.MODEL_NOTES["concise_answers"] not in _solo_prompt(_GLM)


def test_each_note_switches_alone(settings):
    settings({"kit.deepseek": {"concise_answers": False}})
    text = _solo_prompt(_DEEPSEEK)
    assert mp.MODEL_NOTES["concise_answers"] not in text
    assert mp.MODEL_NOTES["terse_tool_rounds"] in text
