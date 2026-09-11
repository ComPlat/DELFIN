"""Field report 2026-09-11 (kit.glm-5.3, dashboard mode, effort high):
"glm führt nicht zu ergebnissen, er denkt ewig lang nach". The turn log
shows the two halves. First the endpoint: 408 s to the first token, then
601 s with nothing. Then ours: "stream went silent for 121 s (stall
budget 120 s)" -- GLM's profile says 420 s, because between two rounds it
thinks, and the dashboard read the budget for the empty model name. The
engine never had a ``model`` attribute; every ``getattr(engine, "model",
"")`` in the dashboard fell back to "". And nothing said that effort
"high" on that model is 989 output tokens for a 224-character answer.
"""

from __future__ import annotations

import pathlib
from unittest.mock import MagicMock, patch

_TAB = (pathlib.Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_agent.py")


def _engine(tmp_path, model="kit.glm-5.3", mode="dashboard"):
    from delfin.agent import engine as E
    with patch("delfin.agent.engine.create_client",
               return_value=MagicMock(model=model)):
        return E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                             model=model, mode=mode)


def test_the_engine_says_which_model_it_talks_to(tmp_path):
    eng = _engine(tmp_path)
    assert eng.model == "kit.glm-5.3"
    from delfin.agent.model_profiles import get_profile
    assert get_profile(eng.model).stale_kill_after_s == 420.0
    assert get_profile("").stale_kill_after_s == 120.0, "what the dashboard used to read"


def test_the_resolver_never_hands_back_an_empty_name_by_accident(tmp_path):
    from delfin.dashboard.tab_agent import _engine_model_name
    eng = _engine(tmp_path)
    assert _engine_model_name(eng, "kit.deepseek-v4-flash") == "kit.glm-5.3"
    assert _engine_model_name(None, "kit.deepseek-v4-flash") == "kit.deepseek-v4-flash"

    class _Bare:
        client = MagicMock(model="kit.glm-5.3")
    assert _engine_model_name(_Bare(), "") == "kit.glm-5.3"

    class _Nothing:
        pass
    assert _engine_model_name(_Nothing(), "fallback") == "fallback"


def test_every_dashboard_read_of_the_model_goes_through_the_resolver():
    text = _TAB.read_text(encoding="utf-8")
    body = text[text.index("def _arm_stale_watcher("):][:14000]
    assert '_model = _engine_model_name(_engine, model_dropdown.value or "")' in body
    assert 'getattr(_engine, "model", "")' not in body
    # the two metrics rows read it the same way
    assert text.count("_engine_model_name(") >= 4


def test_a_reasoning_model_above_its_profile_is_told_what_it_costs():
    from delfin.dashboard.tab_agent import _effort_above_profile_note
    note = _effort_above_profile_note("kit.glm-5.3", "high")
    assert note.startswith("Effort **high** on kit.glm-5.3: its profile recommends **low**")
    assert "/effort low" in note and "your choice" in note
    assert _effort_above_profile_note("kit.glm-5.3", "low") == ""
    assert _effort_above_profile_note("kit.glm-5.3", "") == ""


def test_a_model_that_does_not_reason_is_not_lectured():
    from delfin.dashboard.tab_agent import _effort_above_profile_note
    assert _effort_above_profile_note("kit.deepseek-v4-flash", "xhigh") == ""


def test_the_note_is_posted_when_the_engine_is_built():
    text = _TAB.read_text(encoding="utf-8")
    i = text.index("_probe_endpoint_in_background(engine, provider, model)")
    assert "_effort_above_profile_note(" in text[i:i + 500]
