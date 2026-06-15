"""Tests for provider-agnostic model routing.

The framework contract: ANY model from ANY provider plugs in via the
tier table / settings — no hardcoded single model, opt-in, and never
route into an unavailable or known-broken model.
"""

from __future__ import annotations

import pytest

from delfin.agent.model_routing import (
    RouteDecision,
    is_known_broken,
    route_model,
    tier_model,
    _PROVIDER_TIERS,
)


_ON = {"agent": {"routing": {"enabled": True}}}
_OFF = {"agent": {"routing": {"enabled": False}}}


# ---------------------------------------------------------------------------
# Opt-in: disabled = exactly today's behaviour
# ---------------------------------------------------------------------------

def test_disabled_returns_user_model_untouched():
    d = route_model(provider="kit", user_model="azure.gpt-5.4",
                    complexity="simple", settings=_OFF)
    assert d.model == "azure.gpt-5.4"
    assert d.routed is False


def test_missing_settings_default_to_disabled():
    d = route_model(provider="kit", user_model="M", complexity="complex",
                    settings={})
    assert d.model == "M" and d.routed is False


# ---------------------------------------------------------------------------
# Provider-agnostic tiers: every provider routes within its own namespace
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("provider", sorted(_PROVIDER_TIERS))
def test_simple_routes_to_cheap_tier_per_provider(provider):
    d = route_model(provider=provider, user_model="user-choice",
                    complexity="simple", settings=_ON)
    assert d.routed and d.tier == "cheap"
    assert d.model == _PROVIDER_TIERS[provider]["cheap"]


@pytest.mark.parametrize("provider", sorted(_PROVIDER_TIERS))
def test_complex_routes_to_strong_tier_per_provider(provider):
    d = route_model(provider=provider, user_model="user-choice",
                    complexity="complex", settings=_ON)
    assert d.routed and d.tier == "strong"
    assert d.model == _PROVIDER_TIERS[provider]["strong"]


def test_moderate_keeps_user_model():
    d = route_model(provider="kit", user_model="M", complexity="moderate",
                    settings=_ON)
    assert d.model == "M" and d.routed is False


def test_unknown_provider_falls_back_to_user_model():
    d = route_model(provider="someday-new-provider", user_model="M",
                    complexity="simple", settings=_ON)
    assert d.model == "M" and d.routed is False


# ---------------------------------------------------------------------------
# Settings overrides: global and per-provider — any model plugs in
# ---------------------------------------------------------------------------

def test_global_override_wins_over_builtin():
    s = {"agent": {"routing": {"enabled": True, "cheap_model": "my.tiny"}}}
    d = route_model(provider="kit", user_model="M", complexity="simple",
                    settings=s)
    assert d.model == "my.tiny"


def test_per_provider_override_wins_over_global():
    s = {"agent": {"routing": {
        "enabled": True,
        "cheap_model": "global.cheap",
        "providers": {"ollama": {"cheap_model": "phi:mini"}},
    }}}
    assert tier_model("ollama", "cheap", s) == "phi:mini"
    assert tier_model("kit", "cheap", s) == "global.cheap"


# ---------------------------------------------------------------------------
# Safety: never route into broken / unavailable models
# ---------------------------------------------------------------------------

def test_never_routes_into_known_broken_model(monkeypatch, tmp_path):
    import delfin.agent.model_routing as mr
    # _KNOWN_BROKEN is empty now (azure.gpt-5.5 works again) — use the runtime
    # broken mechanism (mark_broken), which is the durable detector.
    monkeypatch.setattr(mr, "_RUNTIME_BROKEN_PATH", tmp_path / "broken.json")
    mr.mark_broken("azure.broken-model")
    assert is_known_broken("azure.broken-model") is True
    s = {"agent": {"routing": {"enabled": True,
                                "strong_model": "azure.broken-model"}}}
    d = route_model(provider="kit", user_model="M", complexity="complex",
                    settings=s)
    assert d.model == "M" and d.routed is False
    assert "known-broken" in d.reason


def test_unavailable_candidate_falls_back_to_user():
    d = route_model(provider="kit", user_model="M", complexity="simple",
                    settings=_ON, available_models=["M", "other"])
    assert d.model == "M" and d.routed is False


def test_candidate_in_live_list_is_used():
    d = route_model(provider="kit", user_model="M", complexity="simple",
                    settings=_ON,
                    available_models=["M", "azure.gpt-5-nano"])
    assert d.model == "azure.gpt-5-nano" and d.routed


def test_same_as_user_model_is_not_reported_as_routed():
    d = route_model(provider="kit", user_model="azure.gpt-5-nano",
                    complexity="simple", settings=_ON)
    assert d.model == "azure.gpt-5-nano" and d.routed is False


# ---------------------------------------------------------------------------
# Settings default exists (opt-in wiring present)
# ---------------------------------------------------------------------------

def test_routing_settings_default_present_and_disabled():
    from delfin.user_settings import DEFAULT_SETTINGS
    routing = DEFAULT_SETTINGS["agent"]["routing"]
    assert routing["enabled"] is False
    assert routing["strong_model"] == ""
    assert routing["cheap_model"] == ""
