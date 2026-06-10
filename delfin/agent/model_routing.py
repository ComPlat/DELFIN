"""Provider-agnostic model routing — the right model for each turn.

Fundamental framework requirement: ANY model from ANY provider (Anthropic,
OpenAI, KIT Toolbox, Ollama/vLLM, …) must plug in without code changes.
Nothing here hardcodes a single "the" model; routing resolves a **tier**
(strong / cheap) **within the active provider**, from three layers:

1. ``settings["agent"]["routing"]`` — explicit user/team override
   (``strong_model`` / ``cheap_model``, either global or per-provider).
2. Built-in per-provider tier defaults (curated, updated with the
   provider landscape — a *table*, not scattered ``if model == …``).
3. The user's currently selected model (always the final fallback, so
   routing degrades gracefully to today's behaviour).

A model is only routed *to* if it is actually **available** (live model
list when supplied, else assumed) and not on the known-broken list
(e.g. ``azure.gpt-5.5`` 401s on the KIT endpoint).  Routing is **opt-in**
via ``agent.routing.enabled`` — off by default, zero behaviour change.

Per-model *behaviour* (reasoning_effort, tool surface, round caps) stays
in ``model_profiles.py``; this module only decides *which* model runs.
"""

from __future__ import annotations

from dataclasses import dataclass


# Built-in per-provider tier defaults.  Adding a provider = adding a row.
# These are FALLBACKS — settings["agent"]["routing"] overrides them.
_PROVIDER_TIERS: dict[str, dict[str, str]] = {
    "claude": {"strong": "sonnet",         "cheap": "haiku"},
    "openai": {"strong": "gpt-5.4",        "cheap": "gpt-5.4-mini"},
    "kit":    {"strong": "azure.gpt-5.4",  "cheap": "azure.gpt-5-nano"},
    "ollama": {"strong": "qwen3-coder:32b", "cheap": "qwen2.5-coder:7b"},
}

# Models advertised by an endpoint but known to fail (kept short; the
# actionable-401 handler is the runtime detector — this list prevents
# routing INTO a known failure).
_KNOWN_BROKEN: frozenset[str] = frozenset({
    "azure.gpt-5.5",     # listed on the KIT endpoint but 401s (2026-06)
})


@dataclass(frozen=True)
class RouteDecision:
    model: str
    tier: str        # "strong" | "cheap" | "user"
    reason: str      # short human-readable why

    @property
    def routed(self) -> bool:
        return self.tier != "user"


def _routing_settings(settings: dict | None) -> dict:
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    return ((settings or {}).get("agent") or {}).get("routing") or {}


def tier_model(provider: str, tier: str, settings: dict | None = None) -> str:
    """Resolve the model id for ``tier`` within ``provider``.

    Settings override first (per-provider dict or flat global), then the
    built-in provider table.  Empty string when nothing is configured for
    this provider (caller falls back to the user's model).
    """
    cfg = _routing_settings(settings)
    key = f"{tier}_model"
    # Per-provider override: routing.providers.<provider>.strong_model
    per_provider = ((cfg.get("providers") or {}).get(provider) or {})
    if per_provider.get(key):
        return str(per_provider[key])
    # Global override (applies when it belongs to this provider's namespace
    # or the user explicitly set it — we trust explicit config).
    if cfg.get(key):
        return str(cfg[key])
    return _PROVIDER_TIERS.get(provider, {}).get(tier, "")


_RUNTIME_BROKEN_PATH = None  # default: ~/.delfin/broken_models.json


def _runtime_broken_path():
    from pathlib import Path as _P
    return _RUNTIME_BROKEN_PATH or (_P.home() / ".delfin" / "broken_models.json")


def runtime_broken() -> set[str]:
    """Models observed failing (401/unavailable) at runtime — learned, so
    nobody routes/picks them again until the file is cleared."""
    import json as _json
    try:
        p = _runtime_broken_path()
        return set(_json.loads(p.read_text(encoding="utf-8")))
    except Exception:
        return set()


def mark_broken(model: str) -> None:
    """Record a model as broken (called by the 401/unavailable handler)."""
    import json as _json
    m = (model or "").strip()
    if not m:
        return
    try:
        p = _runtime_broken_path()
        cur = runtime_broken()
        if m in cur:
            return
        cur.add(m)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(_json.dumps(sorted(cur)), encoding="utf-8")
    except Exception:
        pass


def is_known_broken(model: str) -> bool:
    m = model or ""
    return m in _KNOWN_BROKEN or m in runtime_broken()


def route_model(
    *,
    provider: str,
    user_model: str,
    complexity: str,
    settings: dict | None = None,
    available_models: list[str] | None = None,
) -> RouteDecision:
    """Pick the model for this turn.

    - routing disabled → the user's model, untouched.
    - complexity "simple"  → provider's cheap tier.
    - complexity "complex" → provider's strong tier.
    - "moderate"           → the user's model (their choice is the default).
    - Candidate not available / known-broken → fall back to the user's
      model (never route INTO a failure).
    """
    cfg = _routing_settings(settings)
    if not cfg.get("enabled"):
        return RouteDecision(user_model, "user", "routing disabled")

    if complexity == "simple":
        tier = "cheap"
    elif complexity == "complex":
        tier = "strong"
    else:
        return RouteDecision(user_model, "user", "moderate → user's choice")

    candidate = tier_model(provider, tier, settings)
    if not candidate or candidate == user_model:
        return RouteDecision(user_model, "user", "no distinct tier model")
    if is_known_broken(candidate):
        return RouteDecision(user_model, "user",
                             f"{candidate} known-broken — kept user model")
    if available_models is not None and candidate not in available_models:
        return RouteDecision(user_model, "user",
                             f"{candidate} not in live model list")
    # Never route the user DOWN from an explicitly chosen strong model to
    # cheap unless they enabled it — enabled IS the consent; proceed.
    return RouteDecision(candidate, tier, f"{complexity} task → {tier} tier")
