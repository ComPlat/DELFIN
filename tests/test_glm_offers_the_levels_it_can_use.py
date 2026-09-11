"""Four field reports (2026-09-07 .. 2026-09-11) ran kit.glm-5.3 at
effort "high"; one of them spent 989 output tokens on a 224-character
reply, and the benchmark arms at medium and low did not differ. Above
the default the level buys hidden reasoning and minutes, not answers.
The profile now states the highest level the model can use, the
dashboard offers only those, a saved choice above the ceiling is
brought down and said so, and the request never carries more.
"""

from __future__ import annotations

import pathlib

from delfin.agent.model_profiles import (
    clamp_effort, effort_choices, get_profile, _COERCE,
)

_TAB = (pathlib.Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_agent.py")


def test_glm_stops_at_medium_and_deepseek_has_no_ceiling():
    assert get_profile("kit.glm-5.3").max_effort == "medium"
    assert get_profile("kit.deepseek-v4-flash").max_effort == ""
    assert effort_choices("kit.glm-5.3") == ("low", "medium")
    assert effort_choices("kit.deepseek-v4-flash") == ("low", "medium", "high", "xhigh")
    assert "max_effort" in _COERCE, "a user may set the ceiling from the settings file"


def test_a_level_above_the_ceiling_is_brought_down():
    assert clamp_effort("kit.glm-5.3", "high") == "medium"
    assert clamp_effort("kit.glm-5.3", "xhigh") == "medium"
    assert clamp_effort("kit.glm-5.3", "low") == "low"
    assert clamp_effort("kit.deepseek-v4-flash", "xhigh") == "xhigh"
    assert clamp_effort("kit.glm-5.3", "") == "" and clamp_effort("kit.glm-5.3", "odd") == "odd"


def test_the_dashboard_resolver_clamps_a_saved_choice():
    from delfin.dashboard.tab_agent import _effort_for_model, _effort_options_for_model
    assert _effort_for_model("kit.glm-5.3", "high") == "medium"
    assert _effort_for_model("kit.glm-5.3", "") == "low"
    assert _effort_for_model("kit.deepseek-v4-flash", "high") == "high"
    opts = _effort_options_for_model("kit.glm-5.3")
    assert [v for _, v in opts] == ["low", "medium"]
    assert any("profile default" in label and v == "low" for label, v in opts)
    assert [v for _, v in _effort_options_for_model("kit.deepseek-v4-flash")] == ["low", "medium", "high", "xhigh"]


def test_the_client_never_sends_a_level_above_the_ceiling():
    import inspect
    from delfin.agent import api_client
    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    i = src.index('_eff_setting = getattr(self, "effort", "")')
    assert "_clamp_eff(self.model, _eff_setting)" in src[i:i + 900]


def test_the_control_and_the_slash_command_follow_the_model():
    text = _TAB.read_text(encoding="utf-8")
    assert text.count("effort_dropdown.options = _effort_options_for_model(") >= 2, (
        "at build time and on every model change")
    i = text.index('if cmd.startswith("/effort "):')
    body = text[i:i + 1200]
    assert "is not on offer for" in body
    assert 'valid = {"low", "medium", "high", "xhigh"}' not in body
