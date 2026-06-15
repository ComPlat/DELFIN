"""Integration tests: the model-profile registry now drives the
existing code paths (compact prompt / core-tool filter / stale-kill).
"""

from __future__ import annotations

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Profile → compact prompt
# ---------------------------------------------------------------------------


def test_prompt_loader_reads_profile_for_compact_decision():
    """Manually setting profile.compact_prompt should make the loader
    use the compact path even for a model name that wouldn't trip the
    legacy weak heuristic."""
    from delfin.agent.prompt_loader import PromptLoader
    from delfin.agent.model_profiles import (
        ModelProfile, register_profile,
    )
    register_profile(
        "custom:strong-but-prefers-compact",
        ModelProfile(compact_prompt=True, core_tools_only=False,
                     notes="test"),
    )
    loader = PromptLoader()
    text = loader.load_role_prompt("solo_agent")
    # The fake-name model isn't matched by _is_weak_model regex
    assert loader._is_weak_model("custom:strong-but-prefers-compact") is False
    # …but the profile flag overrides
    slim = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo",
        model="custom:strong-but-prefers-compact",
    )
    normal = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo", model="opus",
    )
    assert len(slim) < len(normal), (
        f"profile.compact_prompt=True must shrink the prompt: "
        f"slim={len(slim)} normal={len(normal)}"
    )


def test_kit_qwen3_5_keeps_full_prompt_per_profile():
    """Primary KIT target — profile says compact=False, prompt stays full."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    text = loader.load_role_prompt("solo_agent")
    out_strong = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo", model="opus",
    )
    out_qwen = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo",
        model="kit.qwen3.5-397b-A17b",
    )
    # Same size — profile registers it as a strong model
    assert len(out_qwen) == len(out_strong)


# ---------------------------------------------------------------------------
# Profile → core-tool filter (in api_client)
# ---------------------------------------------------------------------------


def test_api_client_reads_profile_for_core_tools_filter():
    """The advertised_tools filter must consult profile.core_tools_only
    instead of just the model-name regex. The filter now also passes the
    resolved capabilities so the weak/strong split can use real facts
    (no native tools / tiny window), not only the name. Check via static
    source scan since we can't construct an OpenAIClient without keys."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "agent" / "api_client.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("Weak-model core-tool filter")
    assert idx > 0
    snippet = text[idx: idx + 2500]
    assert "from .model_profiles import get_profile" in snippet
    assert "_get_profile(self.model, _caps).core_tools_only" in snippet
    assert "_WEAK_MODEL_CORE_TOOLS" in snippet


def test_kit_qwen3_5_profile_says_full_tools():
    from delfin.agent.model_profiles import get_profile
    p = get_profile("kit.qwen3.5-397b-A17b")
    assert p.core_tools_only is False


def test_weak_model_profile_says_core_only():
    """Auto-detected weak model gets the WEAK_DEFAULT which has
    core_tools_only=True."""
    from delfin.agent.model_profiles import get_profile
    p = get_profile("qwen2.5-coder:7b")
    assert p.core_tools_only is True


# ---------------------------------------------------------------------------
# Profile → stale-kill threshold (in tab_agent._arm_stale_watcher)
# ---------------------------------------------------------------------------


def test_arm_stale_watcher_reads_profile_kill_threshold():
    """The watcher must check profile.stale_kill_after_s for dashboard
    mode so per-model kill timing (180s for gemma vs 90s for qwen) is
    honoured. Static source scan."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("Stale-kill threshold")
    assert idx > 0
    snippet = text[idx: idx + 2500]
    assert "from delfin.agent.model_profiles import get_profile" in snippet
    assert "stale_kill_after_s" in snippet
    # The user setting still wins, profile is the second priority
    assert "_user_kill > 0" in snippet
    assert "_profile_kill > 0" in snippet


def test_gemma_profile_has_longer_kill_threshold():
    """Gemma on KIT is slower than other models — profile gives it
    more room (180s vs 120s) so we don't kill legitimate slow runs."""
    from delfin.agent.model_profiles import get_profile
    p = get_profile("kit.gemma4-31b-it")
    assert p.stale_kill_after_s >= 150


def test_qwen3_5_profile_has_aggressive_kill_threshold():
    """Qwen 3.5 is responsive — kill faster on hang to save budget."""
    from delfin.agent.model_profiles import get_profile
    p = get_profile("kit.qwen3.5-397b-A17b")
    assert p.stale_kill_after_s <= 120


# ---------------------------------------------------------------------------
# Profile lookup priority — user-setting always wins
# ---------------------------------------------------------------------------


def test_user_setting_overrides_profile_compact(tmp_path, monkeypatch):
    """``agent.compact_prompt: true`` in user settings must force
    compact mode even for a model whose profile says False."""
    from delfin.agent.prompt_loader import PromptLoader
    from delfin import user_settings as _us
    # Spy on load_settings via monkeypatch
    monkeypatch.setattr(_us, "load_settings",
                        lambda: {"agent": {"compact_prompt": True}})
    loader = PromptLoader()
    text = loader.load_role_prompt("solo_agent")
    out_opus = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo", model="opus",
    )
    # Opus profile says compact=False but the user setting forces it
    monkeypatch.setattr(_us, "load_settings",
                        lambda: {"agent": {"compact_prompt": False}})
    out_opus_normal = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo", model="opus",
    )
    assert len(out_opus) < len(out_opus_normal)
