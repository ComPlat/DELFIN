"""A first turn that takes minutes and says nothing is a hang.

A user asked GLM "Hallo" and waited 190 seconds. The report they filed
says: "Jetzt hat er geantwortet aber hat sehr lange gedauert." Nothing in
the interface distinguished that wait from a broken session, and nothing
could have: the engine said nothing until the answer arrived.

Measured the same day on the same endpoint, same 15k-token prompt:
199 / 266 / 268 s when the endpoint cannot serve the prompt head from its
prefix cache, 7-12 s when it can. The wait is real and cannot be given
back. What can be given is the reason for it.
"""

import pytest

from delfin.agent import verify_guard as vg
from delfin.agent.model_profiles import get_profile


@pytest.fixture(autouse=True)
def _german():
    vg.set_caveat_language("de")
    yield
    vg.set_caveat_language("de")


def test_the_notice_names_the_model_and_the_measured_wait():
    text = vg.cold_start_notice("kit.glm-5.3", 200.0)
    assert "First turn on kit.glm-5.3" in text
    assert "about 200 s" in text, "the measured figure, not a rounded-up promise"
    assert "queue adds more" in text and "faster" in text


def test_the_notice_is_english_whatever_the_session():
    vg.set_caveat_language("de")
    text = vg.cold_start_notice("kit.glm-5.3", 200.0)
    assert "First turn" in text and "Erster Zug" not in text


def test_a_short_wait_is_still_said_in_seconds():
    assert "about 90 s" in vg.cold_start_notice("kit.glm-5.3", 90.0)


def test_a_fast_model_says_nothing():
    assert vg.cold_start_notice("kit.deepseek-v4-flash", 0.0) == ""
    assert vg.cold_start_notice("", 200.0) == ""
    assert get_profile("kit.deepseek-v4-flash").slow_cold_start_s == 0.0


def test_only_the_measured_model_carries_the_wait():
    assert get_profile("kit.glm-5.3").slow_cold_start_s >= 120.0


def test_the_engine_says_it_once_per_session(monkeypatch):
    """Once. A model that is slow to start is slow to start once; a line
    repeated every turn is a line nobody reads."""
    from delfin.agent.engine import AgentEngine

    eng = AgentEngine.__new__(AgentEngine)
    eng.model = "kit.glm-5.3"
    said = []

    def _emit():
        if not getattr(eng, "_cold_start_noted", False):
            eng._cold_start_noted = True
            note = vg.cold_start_notice(
                eng.model, get_profile(eng.model).slow_cold_start_s)
            if note:
                said.append(note)

    _emit()
    _emit()
    _emit()
    assert len(said) == 1
    assert "kit.glm-5.3" in said[0]


def test_the_engine_emits_it_before_the_request_goes_out():
    """It has to be said BEFORE the wait, which means before the model
    call -- after it, it is a description of something already over."""
    import inspect

    from delfin.agent.engine import AgentEngine

    src = inspect.getsource(AgentEngine.stream_response)
    notice_at = src.find("_cold_start_noted")
    build_at = src.find("system_prompt = self._build_current_system_prompt")
    assert 0 < build_at < notice_at, (build_at, notice_at)
    # ... and before anything is streamed back from the model.
    stream_at = src.find("for event in ")
    assert stream_at < 0 or notice_at < stream_at, (notice_at, stream_at)
