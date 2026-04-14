"""Shared + provider-specific learning profiles for the DELFIN agent.

Profiles combine a repo-wide DELFIN baseline (``shared``) with
provider-specific overlays (Claude, OpenAI/Codex, KIT Toolbox). The
merged result drives adaptive routing and thinking budget selection,
while outcome updates stay scoped to the active provider.

Safety:
- Only modifies ``delfin/agent/learned_profiles.json``
- All values are bounded (no runaway drift)
- Rate-limited: max 1 update per cycle
- Fully reversible: reset_profile() restores provider defaults
"""

from __future__ import annotations

import json
import time
from copy import deepcopy
from pathlib import Path
from typing import Any

from delfin.agent.outcome_tracker import CycleOutcome


# Store in the repo so learned profiles are version-controlled and shared
# across machines (pushed with git).  Outcome history stays local (~/.delfin/).
_DEFAULT_PATH = Path(__file__).resolve().parent / "learned_profiles.json"

# Exponential moving average smoothing factor (0.1 = slow, 0.3 = fast)
_EMA_ALPHA = 0.15

# Bounded ranges for self-optimization values
_BOUNDS = {
    "thinking_budget_mult": (0.0, 3.0),
    "success_rate": (0.0, 1.0),
    "avg_cost_per_task": (0.0, 100.0),
}

_PROVIDER_ALIASES = {
    "anthropic": "claude",
    "claude": "claude",
    "codex": "openai",
    "openai": "openai",
    "kit": "kit",
}

_SHARED_PROFILE_KEY = "shared"
_NO_DIFF = object()

# Default profile for new providers
_DEFAULT_PROFILE: dict[str, Any] = {
    "preferred_model": {},
    "thinking_budget_mult": 1.0,
    "avg_cost_per_task": {},
    "success_rate": {},
    "task_performance": {},
    "common_failures": [],
    "total_cycles": 0,
    "updated_at": "",
}


def canonical_provider_name(provider: str) -> str:
    """Normalize provider aliases to the profile storage key."""
    normalized = (provider or "").strip().lower()
    if not normalized:
        return "claude"
    return _PROVIDER_ALIASES.get(normalized, normalized)


def _read_all(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _write_all(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def _clamp(value: float, key: str) -> float:
    """Clamp value to bounded range."""
    lo, hi = _BOUNDS.get(key, (0.0, 100.0))
    return max(lo, min(hi, value))


def _ema(old: float, new: float, alpha: float = _EMA_ALPHA) -> float:
    """Exponential moving average update."""
    return old * (1 - alpha) + new * alpha


def _merge_profile_layers(base: Any, overlay: Any) -> Any:
    """Recursively merge shared and provider-specific profile layers."""
    if isinstance(base, dict) and isinstance(overlay, dict):
        merged = deepcopy(base)
        for key, value in overlay.items():
            if key in merged:
                merged[key] = _merge_profile_layers(merged[key], value)
            else:
                merged[key] = deepcopy(value)
        return merged

    if isinstance(base, list) and isinstance(overlay, list):
        merged = deepcopy(base)
        for item in overlay:
            if item not in merged:
                merged.append(deepcopy(item))
        return merged

    return deepcopy(overlay)


def _profile_overlay(base: Any, merged: Any) -> Any:
    """Return only the provider-specific delta relative to shared defaults."""
    if isinstance(base, dict) and isinstance(merged, dict):
        diff: dict[str, Any] = {}
        for key, value in merged.items():
            if key in base:
                child = _profile_overlay(base[key], value)
                if child is not _NO_DIFF:
                    diff[key] = child
            else:
                diff[key] = deepcopy(value)
        return diff if diff else _NO_DIFF

    if isinstance(base, list) and isinstance(merged, list):
        extras = [deepcopy(item) for item in merged if item not in base]
        return extras if extras else _NO_DIFF

    if base == merged:
        return _NO_DIFF

    return deepcopy(merged)


def _default_profile() -> dict[str, Any]:
    return deepcopy(_DEFAULT_PROFILE)


def _load_shared_profile_from_data(data: dict[str, Any]) -> dict[str, Any]:
    shared = data.get(_SHARED_PROFILE_KEY, {})
    return _merge_profile_layers(_default_profile(), shared)


def _profile_has_context(profile: dict[str, Any]) -> bool:
    if profile.get("total_cycles", 0) > 0:
        return True

    for key in (
        "common_failures",
        "denied_patterns",
        "anti_patterns",
        "communication",
        "tool_usage",
        "domain",
        "playbooks",
        "codebase_map",
        "task_performance",
        "error_patterns",
    ):
        if profile.get(key):
            return True

    return False


def load_shared_profile(path: Path | None = None) -> dict[str, Any]:
    """Load the shared DELFIN baseline profile."""
    data = _read_all(path or _DEFAULT_PATH)
    return _load_shared_profile_from_data(data)


def load_provider_profile(
    provider: str,
    path: Path | None = None,
) -> dict[str, Any]:
    """Load the merged shared + provider-specific profile."""
    provider = canonical_provider_name(provider)
    data = _read_all(path or _DEFAULT_PATH)
    shared = _load_shared_profile_from_data(data)
    provider_profile = data.get(provider, {})
    return _merge_profile_layers(shared, provider_profile)


def load_task_profile(
    provider: str,
    task_class: str,
    path: Path | None = None,
) -> dict[str, Any]:
    """Load the merged task-specific profile for a task class."""
    profile = load_provider_profile(provider, path)
    return profile.get("task_performance", {}).get(task_class or "", {})


def save_provider_profile(
    provider: str,
    profile: dict[str, Any],
    path: Path | None = None,
) -> None:
    """Save only the provider-specific overlay, preserving shared rules."""
    provider = canonical_provider_name(provider)
    p = path or _DEFAULT_PATH
    data = _read_all(p)
    shared = _load_shared_profile_from_data(data)
    profile_to_save = deepcopy(profile)
    profile_to_save["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    overlay = _profile_overlay(shared, profile_to_save)
    data[provider] = overlay if isinstance(overlay, dict) else {}
    _write_all(p, data)


def reset_profile(
    provider: str,
    path: Path | None = None,
) -> None:
    """Reset a provider profile to defaults."""
    save_provider_profile(provider, _default_profile(), path)


def update_from_outcome(
    provider: str,
    outcome: CycleOutcome,
    path: Path | None = None,
) -> dict[str, str]:
    """Update provider profile based on a cycle outcome.

    Returns a dict of changes made (for transparency logging).
    """
    provider = canonical_provider_name(provider)
    profile = load_provider_profile(provider, path)
    changes: dict[str, str] = {}

    profile["total_cycles"] = profile.get("total_cycles", 0) + 1

    # 1. Update success rate per mode
    mode = outcome.mode
    if mode and outcome.verdict:
        old_rate = profile.get("success_rate", {}).get(mode, 0.5)
        new_val = 1.0 if outcome.verdict == "PASS" else 0.0
        new_rate = _clamp(_ema(old_rate, new_val), "success_rate")
        profile.setdefault("success_rate", {})[mode] = round(new_rate, 3)
        if abs(new_rate - old_rate) > 0.01:
            changes["success_rate"] = (
                f"{mode}: {old_rate:.0%} → {new_rate:.0%}"
            )

    # 2. Update average cost per task class
    task_class = outcome.task_class or "general"
    if outcome.cost_usd > 0:
        old_avg = profile.get("avg_cost_per_task", {}).get(
            task_class, outcome.cost_usd
        )
        new_avg = _clamp(
            _ema(old_avg, outcome.cost_usd), "avg_cost_per_task"
        )
        profile.setdefault("avg_cost_per_task", {})[task_class] = round(
            new_avg, 3
        )
        if abs(new_avg - old_avg) > 0.01:
            changes["avg_cost"] = (
                f"{task_class}: ${old_avg:.2f} → ${new_avg:.2f}"
            )

    # 3. Update task-class success rate for adaptive routing/budgeting
    if outcome.task_class and outcome.verdict:
        task_perf = profile.setdefault("task_performance", {}).setdefault(
            outcome.task_class,
            {},
        )
        old_task_rate = task_perf.get("success_rate", 0.5)
        new_task_val = 1.0 if outcome.verdict == "PASS" else 0.0
        new_task_rate = _clamp(
            _ema(old_task_rate, new_task_val), "success_rate"
        )
        task_perf["success_rate"] = round(new_task_rate, 3)
        if abs(new_task_rate - old_task_rate) > 0.01:
            changes["task_success_rate"] = (
                f"{outcome.task_class}: {old_task_rate:.0%} → "
                f"{new_task_rate:.0%}"
            )

    # 4. Track failure patterns
    if outcome.verdict == "FAIL" and outcome.error_type:
        failures = profile.setdefault("common_failures", [])
        if outcome.error_type not in failures:
            failures.append(outcome.error_type)
            # Keep only last 10 failure types
            profile["common_failures"] = failures[-10:]
            changes["new_failure"] = outcome.error_type

    # 5. Track denied command patterns
    if outcome.denied_commands:
        denied = profile.setdefault("denied_patterns", [])
        for cmd in outcome.denied_commands:
            if cmd not in denied:
                denied.append(cmd)
        profile["denied_patterns"] = denied[-20:]

    save_provider_profile(provider, profile, path)
    return changes


def format_profile_context(provider: str, path: Path | None = None) -> str:
    """Format provider profile as context for the system prompt."""
    provider = canonical_provider_name(provider)
    data = _read_all(path or _DEFAULT_PATH)
    shared = _load_shared_profile_from_data(data)
    profile = _merge_profile_layers(shared, data.get(provider, {}))
    provider_overlay = data.get(provider, {})
    total = profile.get("total_cycles", 0)
    if not _profile_has_context(shared) and not _profile_has_context(profile):
        return ""

    parts = []

    shared_failures = shared.get("common_failures", [])
    if shared_failures:
        parts.append(
            f"Shared DELFIN failure patterns: {', '.join(shared_failures[-5:])}"
        )

    shared_rules = shared.get("tool_usage", {}).get("rules", [])
    if shared_rules:
        parts.append(f"Shared DELFIN tool rules: {', '.join(shared_rules[:5])}")

    parts.append(f"Provider: {provider} ({total} cycles)")

    # Success rates
    rates = profile.get("success_rate", {})
    if rates:
        rate_strs = [f"{k}: {v:.0%}" for k, v in sorted(rates.items())]
        parts.append(f"Success rates: {', '.join(rate_strs)}")

    # Common failures
    failures = provider_overlay.get("common_failures", [])
    if failures:
        parts.append(
            f"Provider-specific failure patterns: {', '.join(failures[-5:])}"
        )

    # Denied patterns
    denied = provider_overlay.get("denied_patterns", [])
    if denied:
        parts.append(
            f"Previously denied commands: {', '.join(denied[-5:])}"
        )

    return "\n".join(parts)
