"""Per-model behavioural profiles.

Different models have different sweet-spots for the same framework.
Azure GPT-5.4 is a reasoning model that silently consumes its budget
unless ``reasoning_effort`` is set; KIT's Qwen3.5-397B MoE handles
tool calls perfectly with no special handling; tiny Ollama models
need a much smaller prompt + tool surface. This module lets each
model carry the knobs it needs in one place instead of scattering
``if model.startswith(...)`` across the codebase.

Lookup rules:
1. Exact model-name match wins (``kit.qwen3.5-397b-A17b``).
2. Otherwise a longest-prefix match across registered profiles
   (``azure.gpt-5`` matches ``azure.gpt-5-mini``).
3. Fallback to a tiered default chosen by ``PromptLoader._is_weak_model``:
   weak → ``WEAK_DEFAULT``, otherwise → ``STRONG_DEFAULT``.

Profiles are intentionally compact dataclasses — add fields when
something useful turns up, not preemptively.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace


@dataclass(frozen=True)
class ModelProfile:
    """Per-model behavioural knobs."""

    # Slim the system prompt + collapse prose paragraphs?
    compact_prompt: bool = False

    # Strip the tool schema down to the 15-tool weak-model surface?
    core_tools_only: bool = False

    # Default reasoning_effort for dashboard / solo when user hasn't
    # explicitly picked one. "low" | "medium" | "high" | "xhigh".
    effort_default: str = "medium"

    # Hard cap on auto-continuation rounds in the dashboard's
    # ACTION-execute loop. Weak models hallucinate longer chains.
    max_tool_rounds: int = 50

    # Tool-result truncation cap in KB. Weak models choke on 5KB.
    tool_result_cap_kb: int = 5

    # If True, the agent's output MUST start ACTION lines with the
    # ``ACTION:`` prefix. If False, bare /cmd lines are accepted too
    # (helps weak models who drop the prefix).
    strict_action_prefix: bool = False

    # Default cooperative-stop threshold in seconds (dashboard mode).
    # Reasoning models silently consume tokens for minutes; chat
    # models should respond quickly so we kill earlier.
    stale_kill_after_s: float = 120.0

    # Free-form notes — useful in /agents stats / /model output and
    # for the human reading this file.
    notes: str = ""


# Default tiers. Anything that doesn't have an explicit profile lands
# on one of these.
STRONG_DEFAULT = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="medium",
    max_tool_rounds=50,
    tool_result_cap_kb=5,
    strict_action_prefix=False,
    stale_kill_after_s=120.0,
    notes="Default for strong cloud / large MoE models.",
)

WEAK_DEFAULT = ModelProfile(
    compact_prompt=True,
    core_tools_only=True,
    effort_default="low",
    max_tool_rounds=10,
    tool_result_cap_kb=3,
    strict_action_prefix=False,
    stale_kill_after_s=60.0,
    notes="Default for 7-13B local models (gemma/llama/qwen/phi/mistral).",
)


# --- Concrete profiles ------------------------------------------------------

# kit.qwen3.5-397b-A17b — currently the strongest tool-calling-capable
# model on KIT Toolbox for our scope. 397B MoE with 17B active params,
# excellent agentic tool routing, no silent-reasoning footguns.
_QWEN35_397B = ModelProfile(
    compact_prompt=False,        # full 7.5k slim prompt is fine
    core_tools_only=False,       # handles the 45-tool surface cleanly
    effort_default="medium",     # MoE doesn't need high; medium is sharp
    max_tool_rounds=20,          # rarely needs more than 6-8 in practice
    tool_result_cap_kb=5,
    strict_action_prefix=False,  # accept fault-tolerant form anyway
    stale_kill_after_s=90.0,     # responsive — fail fast if hung
    notes=(
        "KIT Qwen 3.5 397B MoE (17B active). Strong agentic tool use, "
        "no reasoning-effort dance needed. Default choice for the "
        "DELFIN agent on KIT Toolbox."
    ),
)

# kit.gpt-oss-120b — OpenAI's open-weight model. Trained for tool use,
# good but a notch below Qwen3.5-397b on context handling.
# 2026-05-20: core_tools_only=True after P2 27-task baseline showed
# tool-call explosion on chemistry tasks (dash_chemistry_basis_diff
# spent $2.51 / 24 tools / 77s where gemma takes $0.12 / 1 tool / 6s).
# Hypothesis: cap to 15 core tools removes the over-tooling without
# losing the actual chemistry knowledge.
_GPT_OSS_120B = ModelProfile(
    compact_prompt=False,
    core_tools_only=True,
    effort_default="medium",
    max_tool_rounds=15,
    tool_result_cap_kb=5,
    strict_action_prefix=False,
    stale_kill_after_s=120.0,
    notes="KIT GPT-OSS 120B — core_tools_only on (2026-05-20 P2 iter).",
)

# kit.gemma4-31b-it — 31B dense, decent but slower than the MoEs.
_GEMMA4_31B = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="low",
    max_tool_rounds=10,
    tool_result_cap_kb=4,
    strict_action_prefix=False,
    stale_kill_after_s=180.0,    # gemma is slow on KIT — give it room
    notes="KIT Gemma4-31B-IT — dense 31B, slower; low effort by default.",
)

# Azure GPT-5.4 — reasoning model. Without low/medium reasoning_effort
# it spends minutes thinking before emitting tokens. Our budget right-
# sizing in dashboard mode already handles that; this profile pins it.
_AZURE_GPT5 = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="low",         # critical — see api_client reasoning detector
    max_tool_rounds=20,
    tool_result_cap_kb=5,
    strict_action_prefix=False,
    stale_kill_after_s=180.0,    # reasoning legitimately takes time
    notes=(
        "Azure GPT-5.x — reasoning model. Must use reasoning_effort=low "
        "for dashboard or it silently consumes the budget. The api_client "
        "reasoning detector handles this; this profile pins effort default."
    ),
)

# Sonnet — strong all-rounder, no special handling.
_SONNET = ModelProfile(
    compact_prompt=False,
    core_tools_only=False,
    effort_default="medium",
    notes="Sonnet — strong default, no quirks.",
)


# Registry — exact match first, longest-prefix second.
_PROFILES: dict[str, ModelProfile] = {
    # KIT Toolbox
    "kit.qwen3.5-397b-A17b": _QWEN35_397B,
    "kit.gpt-oss-120b": _GPT_OSS_120B,
    "kit.gemma4-31b-it": _GEMMA4_31B,
    # Azure (via KIT or OpenAI)
    "azure.gpt-5.4": _AZURE_GPT5,
    "azure.gpt-5.1": _AZURE_GPT5,
    "azure.gpt-5": _AZURE_GPT5,
    "azure.gpt-5-mini": _AZURE_GPT5,
    "azure.gpt-5-nano": _AZURE_GPT5,
    # Frontier-tier (top capability, no special routing needed)
    "sonnet": _SONNET,
    "opus": replace(_SONNET,
                    notes="Opus — top tier, no quirks."),
    "haiku": replace(_SONNET, effort_default="low",
                     notes="Haiku — fast/cheap, low effort."),
}


# Prefixes that map all matching model names onto the same profile.
# Useful when a provider ships many variants (azure.gpt-5-*).
_PREFIX_PROFILES: tuple[tuple[str, ModelProfile], ...] = (
    ("azure.gpt-5", _AZURE_GPT5),
    ("kit.gpt-oss", _GPT_OSS_120B),
)


def get_profile(model: str) -> ModelProfile:
    """Return the profile for ``model``, falling back to a tier default."""
    if not model:
        return STRONG_DEFAULT
    if model in _PROFILES:
        return _PROFILES[model]
    # Longest-prefix match across the prefix table.
    best: tuple[int, ModelProfile] | None = None
    for prefix, prof in _PREFIX_PROFILES:
        if model.startswith(prefix):
            length = len(prefix)
            if best is None or length > best[0]:
                best = (length, prof)
    if best is not None:
        return best[1]
    # Tier fallback via PromptLoader weak-model heuristic.
    try:
        from .prompt_loader import PromptLoader
        if PromptLoader()._is_weak_model(model):
            return WEAK_DEFAULT
    except Exception:
        pass
    return STRONG_DEFAULT


def register_profile(model: str, profile: ModelProfile) -> None:
    """Register an exact-name profile at runtime. Mostly useful for
    tests and per-session experiments."""
    _PROFILES[model] = profile


def list_profiles() -> list[tuple[str, ModelProfile]]:
    """Return all registered exact-name profiles, sorted by name."""
    return sorted(_PROFILES.items(), key=lambda kv: kv[0])


__all__ = [
    "ModelProfile",
    "STRONG_DEFAULT",
    "WEAK_DEFAULT",
    "get_profile",
    "register_profile",
    "list_profiles",
]
