"""Tests for the per-model behavioural profile registry."""

from __future__ import annotations

import pytest

from delfin.agent.model_profiles import (
    ModelProfile, STRONG_DEFAULT, WEAK_DEFAULT,
    get_profile, register_profile, list_profiles,
)


# ---------------------------------------------------------------------------
# Tier defaults
# ---------------------------------------------------------------------------


def test_strong_default_is_full_surface():
    p = STRONG_DEFAULT
    assert p.compact_prompt is False
    assert p.core_tools_only is False
    assert p.max_tool_rounds >= 30
    assert p.effort_default in ("low", "medium", "high", "xhigh")


def test_weak_default_clamps_everything():
    p = WEAK_DEFAULT
    assert p.compact_prompt is True
    assert p.core_tools_only is True
    assert p.effort_default == "low"
    assert p.max_tool_rounds <= 15
    assert p.stale_kill_after_s <= 90


# ---------------------------------------------------------------------------
# Exact + prefix matching
# ---------------------------------------------------------------------------


def test_kit_qwen3_5_397b_is_registered_as_strong_tool_model():
    """The primary KIT target — must NOT use core-tools-only or
    compact-prompt; Qwen3.5 routes the full surface cleanly."""
    p = get_profile("kit.qwen3.5-397b-A17b")
    assert p.compact_prompt is False
    assert p.core_tools_only is False
    assert p.effort_default == "medium"
    assert "Qwen 3.5 397B" in p.notes


def test_kit_gpt_oss_120b_handled_via_prefix():
    """``kit.gpt-oss-*`` should resolve via the prefix table even for
    a future variant we haven't pinned by exact name."""
    p = get_profile("kit.gpt-oss-120b")
    assert p.compact_prompt is False
    p_future = get_profile("kit.gpt-oss-200b")
    assert p_future.compact_prompt is False
    assert p_future.notes == p.notes


def test_azure_gpt5_variants_all_get_low_effort_default():
    """All Azure GPT-5 family members must default to effort=low so
    they don't burn budget on hidden reasoning."""
    for model in ("azure.gpt-5", "azure.gpt-5.4", "azure.gpt-5-mini",
                  "azure.gpt-5-nano", "azure.gpt-5.1"):
        p = get_profile(model)
        assert p.effort_default == "low", (
            f"{model} effort should be low, got {p.effort_default}"
        )


def test_haiku_gets_low_effort():
    p = get_profile("haiku")
    assert p.effort_default == "low"


def test_opus_keeps_medium_effort():
    p = get_profile("opus")
    assert p.effort_default == "medium"


# ---------------------------------------------------------------------------
# Fallback tier resolution
# ---------------------------------------------------------------------------


def test_unknown_strong_model_falls_back_to_strong_default():
    p = get_profile("some-future-frontier-model")
    assert p is STRONG_DEFAULT or p == STRONG_DEFAULT


def test_unknown_weak_model_falls_back_to_weak_default():
    """A 7B-tier model the profile registry doesn't list explicitly
    should still get the weak tier via PromptLoader._is_weak_model."""
    p = get_profile("qwen2.5-coder:7b")
    assert p.compact_prompt is True
    assert p.core_tools_only is True


def test_empty_model_string_returns_strong_default():
    assert get_profile("") is STRONG_DEFAULT


# ---------------------------------------------------------------------------
# Runtime registration
# ---------------------------------------------------------------------------


def test_register_profile_adds_lookup():
    custom = ModelProfile(
        compact_prompt=True,
        core_tools_only=True,
        effort_default="low",
        notes="test profile",
    )
    register_profile("test:experimental", custom)
    out = get_profile("test:experimental")
    assert out.notes == "test profile"
    assert out.compact_prompt is True


def test_list_profiles_includes_kit_qwen():
    names = {n for n, _ in list_profiles()}
    assert "kit.qwen3.5-397b-A17b" in names
    assert "kit.gemma4-31b-it" in names
    assert "azure.gpt-5.4" in names
    assert "opus" in names


# ---------------------------------------------------------------------------
# Capability-aware tier fallback (get_profile(model, caps))
# ---------------------------------------------------------------------------


def _caps(**kw):
    from delfin.agent.model_capabilities import ModelCapabilities
    base = dict(model="m", provider="ollama", context_window=131_072,
                supports_tools=True)
    base.update(kw)
    return ModelCapabilities(**base)


def test_get_profile_caps_none_is_unchanged_legacy_behaviour():
    # An 8B-by-name model is weak via the name heuristic regardless of caps.
    assert get_profile("llama3.1:8b").core_tools_only is True
    # An unknown strong-looking name stays strong.
    assert get_profile("mystery-model-x").core_tools_only is False


def test_get_profile_caps_no_tools_forces_weak():
    p = get_profile("mystery-model-x", _caps(supports_tools=False))
    assert p is WEAK_DEFAULT


def test_get_profile_caps_large_window_does_not_promote_small_model():
    # 8B model with a 131k (capped) window must STILL be weak — window is not
    # a strength proxy. Regression guard for the core-tools surface.
    p = get_profile("llama3.1:8b", _caps(context_window=131_072))
    assert p.core_tools_only is True


def test_get_profile_caps_tiny_window_forces_weak():
    p = get_profile("mystery-model-x", _caps(context_window=4_096))
    assert p is WEAK_DEFAULT
