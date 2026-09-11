"""The Effort control did nothing for a KIT model on the API backend.

``create_client`` takes ``effort`` and hands it to the CLI backend, which
passes ``--effort`` to a subprocess. The OpenAI-compatible client — the
one the dashboard uses — never received it, so the dropdown moved and
nothing changed. A field bug report shows a session at effort "high"
answering the word "Hallo" in 190 seconds; the setting was not reaching
the model either way.

Probed against the KIT endpoint on 2026-09-07, one trivial question:

    kit.glm-5.3            no setting  25 completion tokens
    kit.glm-5.3            high         8
    kit.glm-5.3            low          3
    kit.deepseek-v4-flash  no setting   2
    kit.deepseek-v4-flash  low         11

Hidden reasoning is most of what a GLM turn spends, and DeepSeek is not a
reasoning line at all — it moves the wrong way. So the parameter is sent
by declared capability, never by default.
"""

import pytest

from delfin.agent.api_client import _reasoning_effort_param
from delfin.agent.model_capabilities import resolve


def _caps(model: str):
    return resolve("kit", model, "", allow_live=False)


def test_a_reasoning_model_gets_the_level():
    caps = _caps("kit.glm-5.3")
    assert caps.is_reasoning
    assert _reasoning_effort_param("low", caps, "kit") == "low"
    assert _reasoning_effort_param("medium", caps, "kit") == "medium"
    assert _reasoning_effort_param("high", caps, "kit") == "high"


def test_the_two_levels_above_high_are_still_high():
    """DELFIN offers five levels, the parameter takes three, and there is
    nothing above high to ask for."""
    caps = _caps("kit.glm-5.3")
    assert _reasoning_effort_param("xhigh", caps, "kit") == "high"
    assert _reasoning_effort_param("max", caps, "kit") == "high"


def test_a_model_that_does_not_reason_is_sent_nothing():
    """DeepSeek spends 2 completion tokens on the question and 11 with the
    parameter set — asking it to reason costs and buys nothing."""
    caps = _caps("kit.deepseek-v4-flash")
    assert not caps.is_reasoning
    assert _reasoning_effort_param("low", caps, "kit") == ""
    assert _reasoning_effort_param("high", caps, "kit") == ""


def test_an_unset_effort_sends_nothing():
    """A session that never touched the control must read exactly what it
    read before."""
    caps = _caps("kit.glm-5.3")
    for value in ("", None, "   ", "nonsense"):
        assert _reasoning_effort_param(value, caps, "kit") == ""


def test_ollama_is_never_sent_the_parameter():
    """It 400s on it — the branch above this one says so already."""
    caps = _caps("kit.glm-5.3")
    assert _reasoning_effort_param("low", caps, "ollama") == ""


def test_missing_capabilities_send_nothing():
    assert _reasoning_effort_param("low", None, "kit") == ""

    class _Broken:
        @property
        def is_reasoning(self):
            raise RuntimeError("no")

    assert _reasoning_effort_param("low", _Broken(), "kit") == ""


def test_the_client_keeps_the_session_effort():
    """It was dropped one call after create_client received it."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.__init__)
    assert "self.effort = (effort or" in src
    create = inspect.getsource(api_client.create_client)
    assert create.count("effort=effort") >= 2


def test_the_parameter_rides_beside_max_tokens_not_instead_of_it():
    """vLLM takes max_tokens; max_completion_tokens is the Azure spelling,
    and the branch that uses it is the name-gated one above."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    i = src.index('kwargs["max_tokens"] = max_tokens')
    tail = src[i:i + 1900]
    assert "_reasoning_effort_param(" in tail
    # The assignment, not the word: the comment beside it names the Azure
    # spelling in order to say why this branch does not use it.
    assert 'kwargs["max_completion_tokens"]' not in tail


def test_the_profile_decides_when_the_user_did_not():
    """effort_default was a field on ModelProfile, documented and tested,
    and read by nothing in the product — there was no path for it to take
    until the parameter existed."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    i = src.index('kwargs["max_tokens"] = max_tokens')
    tail = src[i:i + 1800]
    assert 'if not _eff_setting:' in tail
    assert '.effort_default' in tail
    # The session's own choice still wins over the profile.
    assert tail.index('getattr(self, "effort"') < tail.index(".effort_default")


def test_a_non_reasoning_model_is_unaffected_by_its_profile_default():
    """DeepSeek's profile says medium; sending it would cost and buy
    nothing, so the capability gate has the last word."""
    from delfin.agent.model_profiles import get_profile

    assert get_profile("kit.deepseek-v4-flash").effort_default == "medium"
    caps = _caps("kit.deepseek-v4-flash")
    assert _reasoning_effort_param("medium", caps, "kit") == ""


# ---------------------------------------------------------------------------
# The dashboard's effort follows the model unless the user chose one
# ---------------------------------------------------------------------------
#
# The dropdown starts at "medium" and its value went to the engine for
# every model, so GLM's profile default of "low" -- the one knob that
# shortens its hidden reasoning -- never applied from the dashboard, and
# a "hallo" took minutes (2026-09-11).

def test_the_dashboard_starts_from_the_models_profile_default():
    from delfin.dashboard.tab_agent import _effort_for_model
    from delfin.agent.model_profiles import get_profile
    assert get_profile("kit.glm-5.3").effort_default == "low"
    assert _effort_for_model("kit.glm-5.3", "") == "low"


def test_a_saved_choice_wins_over_the_profile():
    from delfin.dashboard.tab_agent import _effort_for_model
    assert _effort_for_model("kit.glm-5.3", "high") == "high"
    assert _effort_for_model("kit.glm-5.3", "nonsense") == "low"


def test_a_model_without_a_profile_default_falls_back_to_medium():
    from delfin.dashboard.tab_agent import _effort_for_model
    assert _effort_for_model("", "") in ("low", "medium", "high", "xhigh")
    assert _effort_for_model("no-such-model-xyz", "") in ("low", "medium", "high", "xhigh")


def test_the_model_change_follows_the_profile_and_a_sync_is_not_a_choice():
    from pathlib import Path
    text = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    i = text.index("def _on_model_change(change):")
    body = text[i:i + 3000]
    assert "_effort_for_model(change[\"new\"], \"\")" in body
    j = text.index("def _on_effort_change(change):")
    assert "_controls_sync_internal" in text[j:j + 400]
