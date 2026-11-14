"""Characterization tests for the model-capability layer's decisions.

This file pins the values the layer decides *today*, so a later change that
moves code without changing behaviour can be proven neutral: the same tests
must pass before and after. It states the current answer, not a preference —
nothing here should be edited to make a refactor pass; a difference is a
behaviour change and has to be argued on its own.

Pinned, for every profile the registry exposes and for the layer as a whole:

  - the weak/strong decision for each registered model
    (``model_capabilities._is_weak`` — the single heuristic that decides it),
  - the tier flags a weak profile carries (compact prompt, core-tools-only),
  - the per-profile auto-continuation cap (``ModelProfile.max_tool_rounds``),
  - the per-role output-token cap that a request actually carries
    (``engine._ROLE_MAX_TOKENS`` / ``engine._DEFAULT_MAX_TOKENS``),
  - the one KIT endpoint URL the layer owns (``model_capabilities._KIT_BASE_URL``).

The output-token cap pinned here is the one that reaches a request. The
``ModelCapabilities.max_output_tokens`` field is not it: nothing reads that
field, so it decides nothing and is not part of what behaviour depends on.
"""

from __future__ import annotations

import pytest

from delfin.agent import engine
from delfin.agent import model_capabilities as mc
from delfin.agent.model_profiles import list_profiles

# The KIT Toolbox endpoint the capability layer owns, in one place.
_EXPECTED_KIT_BASE_URL = "https://ki-toolbox.scc.kit.edu/api/v1"

# Per registered model: (is_weak, max_tool_rounds, compact_prompt, core_tools_only)
# Recorded from the shipped registry as it stands. A change to any cell is a
# behaviour change, not a refactor.
_PINNED_PROFILES = {
    "azure.gpt-5":           (False, 100, False, False),
    "azure.gpt-5-mini":      (False, 100, False, False),
    "azure.gpt-5-nano":      (False, 100, False, False),
    "azure.gpt-5.1":         (False, 100, False, False),
    "azure.gpt-5.4":         (False, 100, False, False),
    "haiku":                 (False, 100, False, False),
    "kit.deepseek-v4-flash": (False, 100, False, False),
    "kit.gemma4-31b-it":     (False,  10, False, False),
    "kit.glm-5.3":           (False, 100, False, False),
    "kit.gpt-oss-120b":      (False,  15, False, True),
    "kit.qwen3.5-397b-A17b": (False, 100, False, False),
    "opus":                  (False, 100, False, False),
    "sonnet":                (False, 100, False, False),
}

# Per role: the output-token cap a request carries. Pinned so the layer cannot
# silently change how long a turn may be.
_PINNED_ROLE_MAX_TOKENS = {
    "builder_agent": 32768,
    "chief_agent": 8192,
    "critic_agent": 16384,
    "dashboard_agent": 16384,
    "research_agent": 8192,
    "reviewer_agent": 16384,
    "runtime_agent": 8192,
    "session_manager": 16384,
    "solo_agent": 32768,
    "test_agent": 8192,
}

_PINNED_DEFAULT_MAX_TOKENS = 8192


def test_every_pinned_model_is_still_registered():
    """Every shipped profile is reachable from the registry.

    Not an exact-set check: another test module may register a scratch
    profile at runtime, and that must not make this characterization fail.
    What this pins is that the shipped defaults are all still there.
    """
    names = {name for name, _ in list_profiles()}
    missing = sorted(set(_PINNED_PROFILES) - names)
    assert not missing, f"profiles dropped from the registry: {missing}"


@pytest.mark.parametrize("model", sorted(_PINNED_PROFILES))
def test_each_registered_profile_keeps_its_pinned_decisions(model):
    weak, rounds, compact, core_only = _PINNED_PROFILES[model]
    profile = dict(list_profiles())[model]
    # The weak/strong decision (the name heuristic the layer uses).
    assert mc._is_weak(model) is weak
    # What a weak tier flags: the compact prompt and the core-tool surface.
    assert profile.compact_prompt is compact
    assert profile.core_tools_only is core_only
    # The per-profile auto-continuation cap.
    assert profile.max_tool_rounds == rounds


def test_the_per_role_output_token_caps_are_pinned():
    """The cap a request carries per role must not drift."""
    assert engine._ROLE_MAX_TOKENS == _PINNED_ROLE_MAX_TOKENS


def test_the_default_output_token_cap_is_pinned():
    assert engine._DEFAULT_MAX_TOKENS == _PINNED_DEFAULT_MAX_TOKENS


def test_the_kit_endpoint_url_is_one_pinned_value():
    """The layer owns one KIT endpoint URL; pin it so all users agree."""
    assert mc._KIT_BASE_URL == _EXPECTED_KIT_BASE_URL


def test_a_cache_entry_with_the_removed_field_still_loads(tmp_path, monkeypatch):
    """A capability cache written before the field was removed carries
    ``max_output_tokens``; it loads instead of being dropped by a TypeError."""
    import dataclasses
    import json
    import time

    rec = dataclasses.asdict(mc.ModelCapabilities(
        model="m-old", provider="ollama", context_window=4096,
        source="live", discovered_at=time.time()))
    rec["max_output_tokens"] = 8192
    path = tmp_path / "caps.json"
    path.write_text(json.dumps({"_v": mc._CACHE_VERSION,
                                "entries": {"ollama:m-old": rec}}))
    monkeypatch.setattr(mc, "_CACHE_PATH", path)
    monkeypatch.setattr(mc, "_CACHE", {})
    mc._load_disk_cache()
    assert mc._CACHE["ollama:m-old"].context_window == 4096
