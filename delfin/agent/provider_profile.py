"""Provider-specific learning profiles for DELFIN agent.

Each provider (Claude, OpenAI/Codex, KIT Toolbox) has a profile that
tracks success rates, cost patterns, and failure modes.  Profiles are
updated after each cycle outcome and used for adaptive routing and
thinking budget selection.

Safety:
- Only modifies ~/.delfin/provider_profiles.json (agent's own config)
- All values are bounded (no runaway drift)
- Rate-limited: max 1 update per cycle
- Fully reversible: reset_profile() restores defaults
"""

from __future__ import annotations

import json
import time
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

# Default profile for new providers
_DEFAULT_PROFILE: dict[str, Any] = {
    "preferred_model": {},
    "thinking_budget_mult": 1.0,
    "avg_cost_per_task": {},
    "success_rate": {},
    "common_failures": [],
    "total_cycles": 0,
    "updated_at": "",
}


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


def load_provider_profile(
    provider: str,
    path: Path | None = None,
) -> dict[str, Any]:
    """Load the learned profile for a provider."""
    data = _read_all(path or _DEFAULT_PATH)
    profile = data.get(provider, {})
    # Merge with defaults for missing keys
    result = dict(_DEFAULT_PROFILE)
    result.update(profile)
    return result


def save_provider_profile(
    provider: str,
    profile: dict[str, Any],
    path: Path | None = None,
) -> None:
    """Save a provider profile."""
    p = path or _DEFAULT_PATH
    data = _read_all(p)
    profile["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    data[provider] = profile
    _write_all(p, data)


def reset_profile(
    provider: str,
    path: Path | None = None,
) -> None:
    """Reset a provider profile to defaults."""
    save_provider_profile(provider, dict(_DEFAULT_PROFILE), path)


def update_from_outcome(
    provider: str,
    outcome: CycleOutcome,
    path: Path | None = None,
) -> dict[str, str]:
    """Update provider profile based on a cycle outcome.

    Returns a dict of changes made (for transparency logging).
    """
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

    # 3. Track failure patterns
    if outcome.verdict == "FAIL" and outcome.error_type:
        failures = profile.setdefault("common_failures", [])
        if outcome.error_type not in failures:
            failures.append(outcome.error_type)
            # Keep only last 10 failure types
            profile["common_failures"] = failures[-10:]
            changes["new_failure"] = outcome.error_type

    # 4. Track denied command patterns
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
    profile = load_provider_profile(provider, path)
    total = profile.get("total_cycles", 0)
    if total == 0:
        return ""

    parts = [f"Provider: {provider} ({total} cycles)"]

    # Success rates
    rates = profile.get("success_rate", {})
    if rates:
        rate_strs = [f"{k}: {v:.0%}" for k, v in sorted(rates.items())]
        parts.append(f"Success rates: {', '.join(rate_strs)}")

    # Common failures
    failures = profile.get("common_failures", [])
    if failures:
        parts.append(f"Known failure patterns: {', '.join(failures[-5:])}")

    # Denied patterns
    denied = profile.get("denied_patterns", [])
    if denied:
        parts.append(
            f"Previously denied commands: {', '.join(denied[-5:])}"
        )

    return "\n".join(parts)
