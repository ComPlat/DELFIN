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
import re
import time
from copy import deepcopy
from pathlib import Path
from typing import Any

from delfin.agent.outcome_tracker import CycleOutcome


# Store in the repo so learned profiles are version-controlled and shared
# across machines (pushed with git).  Outcome history stays local (~/.delfin/).
_DEFAULT_PATH = Path(__file__).resolve().parent / "learned_profiles.json"
_PLAYBOOKS_PATH = Path(__file__).resolve().parent / "profile_playbooks.json"
_LOCAL_STATE_PATH = Path.home() / ".delfin" / "provider_profile_state.json"
_LOCAL_ONLY_KEYS = {"next_steps", "denied_patterns"}

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

_PLAYBOOK_STOPWORDS = {"editing", "changes", "change"}

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


def _resolve_local_state_path(
    profile_path: Path | None = None,
    local_state_path: Path | None = None,
) -> Path:
    if local_state_path is not None:
        return local_state_path
    if profile_path is not None and profile_path != _DEFAULT_PATH:
        return profile_path.with_name(f"{profile_path.stem}.local.json")
    return _LOCAL_STATE_PATH


def _read_local_state(path: Path | None = None) -> dict[str, Any]:
    return _read_all(path or _LOCAL_STATE_PATH)


def _write_local_state(data: dict[str, Any], path: Path | None = None) -> None:
    _write_all(path or _LOCAL_STATE_PATH, data)


def _split_local_overlay(profile: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    repo_profile = deepcopy(profile)
    local_profile: dict[str, Any] = {}
    for key in list(repo_profile.keys()):
        if key in _LOCAL_ONLY_KEYS:
            local_profile[key] = repo_profile.pop(key)
    return repo_profile, local_profile


def _merge_local_overlay(provider: str, profile: dict[str, Any], local_state_path: Path | None = None) -> dict[str, Any]:
    local_state = _read_local_state(local_state_path)
    local_overlay = local_state.get(provider, {})
    if not isinstance(local_overlay, dict):
        return profile
    return _merge_profile_layers(profile, local_overlay)


def _load_playbook_data(path: Path | None = None) -> dict[str, Any]:
    return _read_all(path or _PLAYBOOKS_PATH)


def _load_playbook_catalog(
    provider: str,
    path: Path | None = None,
    profile_data: dict[str, Any] | None = None,
) -> dict[str, Any]:
    provider = canonical_provider_name(provider)
    data = profile_data or _read_all(_DEFAULT_PATH)
    shared = data.get(_SHARED_PROFILE_KEY, {})
    provider_block = data.get(provider, {})
    legacy_shared = {
        "playbooks": shared.get("playbooks", {}),
        "codebase_map": shared.get("codebase_map", {}),
    }
    legacy_provider = {
        "playbooks": provider_block.get("playbooks", {}),
        "codebase_map": provider_block.get("codebase_map", {}),
    }
    if (
        path is None
        and profile_data is not None
        and (legacy_shared["playbooks"] or legacy_provider["playbooks"])
    ):
        return _merge_profile_layers(legacy_shared, legacy_provider)

    catalog = _load_playbook_data(path)
    if catalog:
        shared = catalog.get("shared", {})
        provider_block = catalog.get(provider, {})
        return _merge_profile_layers(shared, provider_block)

    # Backward-compatible fallback: load legacy embedded playbooks from profile data.
    return _merge_profile_layers(legacy_shared, legacy_provider)


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
        "next_steps",
        "playbooks",
        "codebase_map",
        "task_performance",
        "error_patterns",
    ):
        if profile.get(key):
            return True

    return False


def _tokenize_identifier(value: str) -> set[str]:
    """Split identifiers like ``dashboard/tab_agent.py`` into search tokens."""
    tokens = {
        token
        for token in re.split(r"[^a-z0-9]+", (value or "").lower())
        if token and token not in _PLAYBOOK_STOPWORDS
    }
    return tokens


def _playbook_aliases(playbook_key: str, module_name: str) -> list[str]:
    """Return canonical aliases that may appear in a task description."""
    aliases = {
        playbook_key.lower(),
        playbook_key.lower().removesuffix("_editing"),
        playbook_key.lower().removesuffix("_changes"),
        module_name.lower(),
    }
    module_path = Path(module_name)
    aliases.add(module_path.name.lower())
    aliases.add(module_path.stem.lower())
    if module_name and not module_name.startswith("delfin/"):
        aliases.add(f"delfin/{module_name.lower()}")
    return [alias for alias in aliases if alias]


def _find_relevant_playbook(
    profile: dict[str, Any],
    task_text: str,
) -> tuple[str, dict[str, Any], str] | None:
    """Resolve the best-matching playbook for a task description."""
    lower_task = (task_text or "").lower()
    if not lower_task.strip():
        return None

    playbooks = profile.get("playbooks", {})
    if not isinstance(playbooks, dict) or not playbooks:
        return None

    modules = profile.get("codebase_map", {}).get("modules", {})
    module_names = list(modules.keys()) if isinstance(modules, dict) else []
    best: tuple[int, int, str, dict[str, Any], str] | None = None

    for playbook_key, playbook in playbooks.items():
        playbook_tokens = _tokenize_identifier(playbook_key)
        matched_module = ""
        matched_score = -1
        for module_name in module_names:
            overlap = len(
                playbook_tokens & _tokenize_identifier(module_name)
            )
            if overlap > matched_score:
                matched_score = overlap
                matched_module = module_name

        for alias in _playbook_aliases(playbook_key, matched_module):
            if "/" in alias or "." in alias:
                position = lower_task.find(alias)
                if position < 0:
                    continue
            else:
                match = re.search(
                    rf"(?<![a-z0-9_]){re.escape(alias)}(?![a-z0-9_])",
                    lower_task,
                )
                if not match:
                    continue
                position = match.start()

            candidate = (
                position,
                -len(alias),
                playbook_key,
                playbook,
                matched_module,
            )
            if best is None or candidate < best:
                best = candidate

    if best is None:
        return None

    _, _, key, playbook, module_name = best
    return key, playbook, module_name


def load_shared_profile(path: Path | None = None) -> dict[str, Any]:
    """Load the shared DELFIN baseline profile."""
    data = _read_all(path or _DEFAULT_PATH)
    return _load_shared_profile_from_data(data)


def load_provider_profile(
    provider: str,
    path: Path | None = None,
    local_state_path: Path | None = None,
) -> dict[str, Any]:
    """Load the merged shared + provider-specific profile."""
    provider = canonical_provider_name(provider)
    data = _read_all(path or _DEFAULT_PATH)
    shared = _load_shared_profile_from_data(data)
    provider_profile = data.get(provider, {})
    merged = _merge_profile_layers(shared, provider_profile)
    resolved_local = _resolve_local_state_path(path, local_state_path)
    return _merge_local_overlay(provider, merged, resolved_local)


def load_task_profile(
    provider: str,
    task_class: str,
    path: Path | None = None,
    local_state_path: Path | None = None,
) -> dict[str, Any]:
    """Load the merged task-specific profile for a task class."""
    profile = load_provider_profile(provider, path, local_state_path)
    return profile.get("task_performance", {}).get(task_class or "", {})


def _format_next_step_entry(step: Any) -> str:
    """Render a persisted next-step entry for prompt context."""
    if isinstance(step, str):
        return step.strip()
    if not isinstance(step, dict):
        return ""

    task = str(step.get("task", "")).strip()
    status = str(step.get("status", "")).strip()
    why = str(step.get("why", "")).strip()

    parts: list[str] = []
    if status and task:
        parts.append(f"[{status}] {task}")
    elif task:
        parts.append(task)
    elif status:
        parts.append(f"[{status}]")

    if why:
        if parts:
            parts.append(f"({why})")
        else:
            parts.append(why)

    return " ".join(parts).strip()


def save_provider_profile(
    provider: str,
    profile: dict[str, Any],
    path: Path | None = None,
    local_state_path: Path | None = None,
) -> None:
    """Save only the provider-specific overlay, preserving shared rules."""
    provider = canonical_provider_name(provider)
    p = path or _DEFAULT_PATH
    resolved_local = _resolve_local_state_path(p, local_state_path)
    data = _read_all(p)
    shared = _load_shared_profile_from_data(data)
    profile_to_save = deepcopy(profile)
    profile_to_save["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    repo_profile, local_profile = _split_local_overlay(profile_to_save)
    overlay = _profile_overlay(shared, repo_profile)
    data[provider] = overlay if isinstance(overlay, dict) else {}
    _write_all(p, data)

    local_state = _read_local_state(resolved_local)
    if local_profile:
        local_state[provider] = local_profile
    else:
        local_state.pop(provider, None)
    _write_local_state(local_state, resolved_local)


def reset_profile(
    provider: str,
    path: Path | None = None,
    local_state_path: Path | None = None,
) -> None:
    """Reset a provider profile to defaults."""
    save_provider_profile(provider, _default_profile(), path, local_state_path)


def update_from_outcome(
    provider: str,
    outcome: CycleOutcome,
    path: Path | None = None,
    local_state_path: Path | None = None,
) -> dict[str, str]:
    """Update provider profile based on a cycle outcome.

    Returns a dict of changes made (for transparency logging).
    """
    provider = canonical_provider_name(provider)
    profile = load_provider_profile(provider, path, local_state_path)
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

    save_provider_profile(provider, profile, path, local_state_path)
    return changes


def format_profile_context(
    provider: str,
    path: Path | None = None,
    *,
    mode_id: str = "",
) -> str:
    """Format provider profile as context for the system prompt."""
    provider = canonical_provider_name(provider)
    data = _read_all(path or _DEFAULT_PATH)
    shared = _load_shared_profile_from_data(data)
    resolved_local = _resolve_local_state_path(path, None)
    profile = _merge_local_overlay(
        provider,
        _merge_profile_layers(shared, data.get(provider, {})),
        resolved_local,
    )
    provider_overlay = data.get(provider, {})
    total = profile.get("total_cycles", 0)
    if not _profile_has_context(shared) and not _profile_has_context(profile):
        return ""

    parts = []

    shared_failures = shared.get("common_failures", [])
    if shared_failures:
        parts.append(
            f"Shared failures: {', '.join(shared_failures[-2:])}"
        )

    shared_rules = shared.get("tool_usage", {}).get("rules", [])
    if shared_rules:
        parts.append(f"Shared tool rules: {', '.join(shared_rules[:2])}")

    parts.append(f"Provider: {provider} ({total} cycles)")

    rates = profile.get("success_rate", {})
    if rates:
        if mode_id and mode_id in rates:
            parts.append(f"Mode success: {mode_id}={rates[mode_id]:.0%}")
        else:
            top_rates = sorted(rates.items(), key=lambda item: item[1], reverse=True)[:2]
            parts.append(
                "Best mode success: "
                + ", ".join(f"{k}={v:.0%}" for k, v in top_rates)
            )

    failures = provider_overlay.get("common_failures", [])
    if failures:
        parts.append(
            f"Provider failures: {', '.join(failures[-2:])}"
        )

    communication_rules = (
        provider_overlay.get("communication", {}).get("rules", [])
    )
    if communication_rules:
        parts.append(
            "Communication: " + ", ".join(str(rule) for rule in communication_rules[:2])
        )

    provider_tool_rules = (
        provider_overlay.get("tool_usage", {}).get("rules", [])
    )
    if provider_tool_rules:
        parts.append(
            "Tool rules: " + ", ".join(str(rule) for rule in provider_tool_rules[:2])
        )

    denied = profile.get("denied_patterns", [])
    if denied:
        parts.append(
            f"Denied: {', '.join(denied[-2:])}"
        )

    next_steps = profile.get("next_steps", [])
    if isinstance(next_steps, list):
        rendered_steps = [
            rendered
            for rendered in (_format_next_step_entry(step) for step in next_steps[:1])
            if rendered
        ]
        if rendered_steps:
            parts.append(f"Next: {'; '.join(rendered_steps)}")

    return "\n".join(parts)


def format_relevant_playbook_context(
    provider: str,
    task_text: str,
    path: Path | None = None,
    *,
    playbooks_path: Path | None = None,
) -> str:
    """Format the single most relevant playbook for the current task."""
    provider = canonical_provider_name(provider)
    profile_data = _read_all(path or _DEFAULT_PATH)
    catalog = _load_playbook_catalog(provider, playbooks_path, profile_data)
    resolved = _find_relevant_playbook(catalog, task_text)
    if not resolved:
        return ""

    playbook_key, playbook, module_name = resolved
    tests = (
        catalog.get("codebase_map", {})
        .get("test_mapping", {})
        .get(module_name, [])
    )

    lines = ["Relevant Playbook"]
    if module_name:
        lines.append(f"Target module: {module_name}")
    description = playbook.get("description", "")
    if description:
        lines.append(f"Focus: {description}")

    steps = playbook.get("steps", [])
    if steps:
        lines.append("Steps:")
        lines.extend(str(step) for step in steps)

    invariants = playbook.get("key_invariants", [])
    if invariants:
        lines.append("Key invariants:")
        lines.extend(f"- {item}" for item in invariants)

    if tests:
        lines.append(f"Related tests: {', '.join(str(test) for test in tests)}")

    lines.append(f"Playbook key: {playbook_key}")
    return "\n".join(lines)
