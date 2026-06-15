"""Per-model capability resolution — the single source of truth for a
model's *real* context window and what it can actually do.

Fundamental framework requirement (mirrors ``model_routing``): ANY model
from ANY provider (Anthropic, OpenAI, KIT Toolbox, Ollama/vLLM/LM-Studio,
…) must run at its FULL potential without code changes. Two failures block
that today and this module fixes both:

1. The engine hard-codes a 100k context window for every model. Small local
   models overflow; large models (128k/200k/1M) are throttled to 100k.
2. Ollama's OpenAI-compatible endpoint silently truncates to ``num_ctx``
   (2-4k by default) unless we pass ``options.num_ctx`` — so even a
   128k-capable model only sees a few thousand tokens.

``resolve(provider, model, base_url)`` returns a :class:`ModelCapabilities`
with the real window plus tool / vision / reasoning flags, discovered LIVE
where possible (Ollama ``/api/show``, OpenAI-compatible ``/v1/models``) and
falling back to a curated static table then name heuristics. Results are
cached (process + ``~/.delfin`` overlay, 24h TTL).

CRITICAL Ollama invariant: ``context_window == num_ctx_override``. Ollama
only attends to ``num_ctx`` tokens, so the compaction budget MUST equal what
we actually told it to use — never the model's theoretical maximum.

Network calls use ``urllib.request`` (no new dependency), short timeouts,
and never raise: any failure degrades to static/heuristic. Per-model
*behaviour* knobs (effort, round caps, prompt size) stay in
``model_profiles``; this module only resolves *facts* about a model.
"""

from __future__ import annotations

import json
import time
import urllib.request
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Tunables
# ---------------------------------------------------------------------------

# Default ceiling for Ollama ``num_ctx``. Exposing a model's full window
# (e.g. 128k) on a local GPU multiplies VRAM use and can OOM a shared box,
# so we cap by default and let users raise it per deployment via
# ``settings["agent"]["ollama"]["num_ctx"]``.
DEFAULT_NUM_CTX_CAP = 32_768

# Conservative fallback windows when nothing better is known.
_WEAK_FALLBACK_WINDOW = 8_192
_STRONG_FALLBACK_WINDOW = 32_768

# Engine's historical default — used as the absolute last resort so behaviour
# never regresses below today's.
_LEGACY_WINDOW = 100_000

_CACHE_TTL_S = 24 * 3600
_CACHE_PATH = Path.home() / ".delfin" / "model_capabilities_cache.json"

# KIT Toolbox hosts many models; only a few are worth driving the agent
# (strong agentic tool use, large window). Selecting a weak one should warn
# and point at the best. ``KIT_BEST_MODEL`` is the curated default; the
# recommended set is the small pool that "lohnt sich" for agent work.
KIT_BEST_MODEL = "kit.qwen3.5-397b-A17b"
_KIT_RECOMMENDED: frozenset[str] = frozenset({
    "kit.qwen3.5-397b-A17b",   # best: strong agentic tool routing, no footguns
    "kit.gpt-oss-120b",        # solid, tool-trained — a notch below qwen3.5
})
# Below this window a KIT model is not worth the agent regardless of tools.
_KIT_MIN_WINDOW = 32_000

# Network timeouts (seconds). Tags is a cheap liveness probe; show parses
# metadata; models is the OpenAI-compatible list.
_TIMEOUT_TAGS = 1.5
_TIMEOUT_SHOW = 4.0
_TIMEOUT_MODELS = 3.0


# ---------------------------------------------------------------------------
# Capability record
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ModelCapabilities:
    """Resolved facts about one (provider, model) pair."""

    model: str
    provider: str                       # "claude"|"openai"|"kit"|"ollama"|…
    context_window: int                 # real usable token window
    max_output_tokens: int = 8_192
    supports_tools: bool = True
    supports_vision: bool = False
    is_reasoning: bool = False          # o-series / gpt-5 / deepseek-r1 / qwq …
    thinking_tagged: bool = False       # emits <think>…</think> in text channel
    recommended_effort: str = "medium"  # "low"|"medium"|"high"
    num_ctx_override: int | None = None  # value for options.num_ctx (ollama only)
    source: str = "static"              # "live"|"static"|"heuristic"|"fallback"
    discovered_at: float = 0.0


# ---------------------------------------------------------------------------
# Curated static table (fallback when live discovery is unreachable)
# ---------------------------------------------------------------------------

# Partial specs keyed by exact model name. Anything omitted takes the
# dataclass default. Windows are deliberately conservative — they drive the
# compaction threshold, so under-shooting just compacts a little earlier.
_STATIC: dict[str, dict[str, Any]] = {
    # Anthropic
    "sonnet": {"context_window": 200_000, "supports_vision": True},
    "opus": {"context_window": 200_000, "supports_vision": True},
    "haiku": {"context_window": 200_000, "supports_vision": True,
              "recommended_effort": "low"},
    # OpenAI / Azure GPT-5 (reasoning family)
    "gpt-5.4": {"context_window": 256_000, "is_reasoning": True,
                "supports_vision": True, "recommended_effort": "low"},
    "gpt-5": {"context_window": 256_000, "is_reasoning": True,
              "supports_vision": True, "recommended_effort": "low"},
    "o3": {"context_window": 200_000, "is_reasoning": True,
           "recommended_effort": "medium"},
    "o4-mini": {"context_window": 200_000, "is_reasoning": True,
                "supports_vision": True, "recommended_effort": "medium"},
    "gpt-4.1": {"context_window": 200_000, "supports_vision": True},
    # KIT Toolbox (served via vLLM — windows reflect the deployment config)
    "kit.qwen3.5-397b-A17b": {"context_window": 128_000},
    "kit.gpt-oss-120b": {"context_window": 128_000},
    "kit.gemma4-31b-it": {"context_window": 128_000, "recommended_effort": "low"},
    # Common Ollama tags (theoretical max; the ollama branch re-caps to num_ctx)
    "qwen3-coder:32b": {"context_window": 256_000},
    "qwen2.5-coder:7b": {"context_window": 131_072},
    "qwen2.5-coder:14b": {"context_window": 131_072},
    "qwen2.5-coder:32b": {"context_window": 131_072},
    "llama3.3:70b": {"context_window": 131_072, "supports_vision": False},
    "llama3.1:8b": {"context_window": 131_072},
    "deepseek-r1:7b": {"context_window": 131_072, "is_reasoning": True,
                       "thinking_tagged": True},
    "deepseek-r1:32b": {"context_window": 131_072, "is_reasoning": True,
                        "thinking_tagged": True},
    "qwq:32b": {"context_window": 131_072, "is_reasoning": True,
                "thinking_tagged": True},
    "llava:7b": {"context_window": 32_768, "supports_vision": True},
}

# Longest-prefix fallbacks (provider many-variant families).
_STATIC_PREFIX: tuple[tuple[str, dict[str, Any]], ...] = (
    ("azure.gpt-5", {"context_window": 256_000, "is_reasoning": True,
                     "supports_vision": True, "recommended_effort": "low"}),
    ("kit.gpt-oss", {"context_window": 128_000}),
    ("gpt-5", {"context_window": 256_000, "is_reasoning": True,
               "supports_vision": True, "recommended_effort": "low"}),
    ("deepseek-r1", {"context_window": 131_072, "is_reasoning": True,
                     "thinking_tagged": True}),
    ("qwq", {"context_window": 131_072, "is_reasoning": True,
             "thinking_tagged": True}),
    ("qwen3-coder", {"context_window": 256_000}),
    ("qwen2.5-coder", {"context_window": 131_072}),
    ("llama3", {"context_window": 131_072}),
    ("llava", {"context_window": 32_768, "supports_vision": True}),
)


def _static_spec(model: str) -> dict[str, Any] | None:
    """Return the static partial spec for ``model`` (exact then prefix)."""
    if not model:
        return None
    if model in _STATIC:
        return _STATIC[model]
    base = model.split(".", 1)[-1] if model.startswith(("azure.", "kit.")) else model
    if base in _STATIC:
        return _STATIC[base]
    best: tuple[int, dict[str, Any]] | None = None
    for prefix, spec in _STATIC_PREFIX:
        if model.startswith(prefix) or base.startswith(prefix):
            if best is None or len(prefix) > best[0]:
                best = (len(prefix), spec)
    return best[1] if best else None


# ---------------------------------------------------------------------------
# Heuristics (last resort when neither live nor static knows the model)
# ---------------------------------------------------------------------------

import re as _re

_REASONING_NAME = _re.compile(
    r"deepseek-?r1|qwq|:thinking|-thinking|^o\d|gpt-5", _re.IGNORECASE
)
_THINKING_TAG_NAME = _re.compile(
    r"deepseek-?r1|qwq|:thinking|-thinking|qwen3", _re.IGNORECASE
)


def _is_weak(model: str) -> bool:
    """Reuse the PromptLoader weak-model name heuristic (single source)."""
    try:
        from .prompt_loader import PromptLoader
        return PromptLoader()._is_weak_model(model)
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


def _configured_num_ctx_cap() -> int:
    """User-tunable Ollama num_ctx ceiling (``agent.ollama.num_ctx``)."""
    try:
        from delfin.user_settings import load_settings
        val = (((load_settings() or {}).get("agent") or {})
               .get("ollama") or {}).get("num_ctx")
        if val and int(val) > 0:
            return int(val)
    except Exception:
        pass
    return DEFAULT_NUM_CTX_CAP


# ---------------------------------------------------------------------------
# Live discovery
# ---------------------------------------------------------------------------


def _ollama_root(base_url: str) -> str:
    """Strip the OpenAI-compat ``/v1`` suffix to reach Ollama's native API."""
    root = (base_url or "").rstrip("/")
    if root.endswith("/v1"):
        root = root[: -len("/v1")]
    return root or "http://localhost:11434"


def _http_get_json(url: str, timeout: float) -> Any:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310
        return json.loads(resp.read().decode("utf-8"))


def _http_post_json(url: str, payload: dict, timeout: float) -> Any:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=data,
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310
        return json.loads(resp.read().decode("utf-8"))


def _ollama_installed(root: str) -> list[str] | None:
    """Installed model tags via ``/api/tags`` (None if unreachable)."""
    try:
        data = _http_get_json(f"{root}/api/tags", _TIMEOUT_TAGS)
        return [m.get("name", "") for m in (data.get("models") or [])]
    except Exception:
        return None


def _discover_ollama(base_url: str, model: str) -> dict[str, Any] | None:
    """Query Ollama ``/api/show`` for context_length + capabilities.

    Returns a partial spec (context_window/supports_*/is_reasoning/…) or
    None when the endpoint is unreachable / the model isn't pulled.
    """
    root = _ollama_root(base_url)
    try:
        data = _http_post_json(f"{root}/api/show", {"model": model}, _TIMEOUT_SHOW)
    except Exception:
        return None
    if not isinstance(data, dict):
        return None

    spec: dict[str, Any] = {}

    # context length — model_info carries a "<arch>.context_length" key.
    info = data.get("model_info") or {}
    ctx = 0
    for key, val in info.items():
        if key.endswith(".context_length"):
            try:
                ctx = int(val)
            except (TypeError, ValueError):
                ctx = 0
            break
    if ctx > 0:
        spec["context_window"] = ctx

    # capabilities array → tools / vision / thinking
    caps = data.get("capabilities") or []
    if isinstance(caps, list) and caps:
        spec["supports_tools"] = "tools" in caps
        spec["supports_vision"] = "vision" in caps
        if "thinking" in caps:
            spec["is_reasoning"] = True
            spec["thinking_tagged"] = True

    return spec or None


def _models_url(base_url: str) -> str:
    url = (base_url or "").rstrip("/")
    return (url if url.endswith("/v1") else url + "/v1") + "/models"


def _fetch_openai_models(base_url: str) -> list[dict] | None:
    """Raw ``/v1/models`` entries (None if unreachable)."""
    try:
        data = _http_get_json(_models_url(base_url), _TIMEOUT_MODELS)
        entries = data.get("data") or []
        return [e for e in entries if isinstance(e, dict)]
    except Exception:
        return None


def _model_matches(entry_id: str, model: str) -> bool:
    """Lenient match between a /v1/models id and our model name.

    Handles DELFIN's provider prefixes (``kit.``/``azure.``) which the
    served id does not carry, plus ``:tag`` suffixes.
    """
    a = (entry_id or "").strip().lower()
    b = (model or "").strip().lower()
    if not a or not b:
        return False
    b_base = b.split(".", 1)[-1] if b.startswith(("kit.", "azure.")) else b
    a0, b0 = a.split(":")[0], b_base.split(":")[0]
    return a == b_base or a0 == b0 or b_base in a or a in b_base


def _discover_openai_models(base_url: str, model: str) -> dict[str, Any] | None:
    """Refine the context window from an OpenAI-compatible ``/v1/models``.

    vLLM (which serves KIT Toolbox) reports ``max_model_len`` per entry;
    some backends use ``context_length``/``max_context_length``. When the
    matched entry carries one, it is the server's TRUE window — authoritative
    over the static table. Returns a partial spec or None.
    """
    entries = _fetch_openai_models(base_url)
    if not entries:
        return None

    def _win(e: dict) -> int:
        for k in ("max_model_len", "context_length", "max_context_length",
                  "context_window"):
            try:
                v = int(e.get(k) or 0)
            except (TypeError, ValueError):
                v = 0
            if v > 0:
                return v
        return 0

    matched = [e for e in entries if _model_matches(e.get("id", ""), model)]
    pool = matched or (entries if len(entries) == 1 else [])
    for e in pool:
        w = _win(e)
        if w > 0:
            return {"context_window": w}
    return None


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

_CACHE: dict[str, ModelCapabilities] = {}


def _cache_key(provider: str, model: str, base_url: str) -> str:
    return f"{provider}\x1f{model}\x1f{base_url}"


def _load_disk_cache() -> None:
    try:
        raw = json.loads(_CACHE_PATH.read_text(encoding="utf-8"))
    except Exception:
        return
    now = time.time()
    for key, rec in (raw or {}).items():
        try:
            if now - float(rec.get("discovered_at", 0)) > _CACHE_TTL_S:
                continue
            _CACHE.setdefault(key, ModelCapabilities(**rec))
        except Exception:
            continue


def _save_disk_cache() -> None:
    try:
        _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        out = {k: v.__dict__ for k, v in _CACHE.items() if v.source == "live"}
        _CACHE_PATH.write_text(json.dumps(out), encoding="utf-8")
    except Exception:
        pass


_disk_loaded = False


def clear_cache() -> None:
    """Drop all cached capabilities (process + best-effort disk)."""
    _CACHE.clear()
    try:
        _CACHE_PATH.unlink()
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


def _apply_ollama_window(spec: dict[str, Any]) -> dict[str, Any]:
    """Re-cap an Ollama spec so context_window == num_ctx_override.

    Ollama only attends to ``num_ctx`` tokens, so the compaction budget
    (context_window) must equal what we actually send as ``options.num_ctx``
    — never the model's theoretical maximum.
    """
    model_max = int(spec.get("context_window") or _STRONG_FALLBACK_WINDOW)
    cap = _configured_num_ctx_cap()
    num_ctx = max(1, min(model_max, cap))
    spec = dict(spec)
    spec["num_ctx_override"] = num_ctx
    spec["context_window"] = num_ctx
    return spec


def resolve(
    provider: str,
    model: str,
    base_url: str = "",
    *,
    allow_live: bool = True,
) -> ModelCapabilities:
    """Resolve capabilities for ``(provider, model)``.

    Order: cache → live discovery → static table → name heuristic. Never
    raises; degrades to a safe fallback so the engine always gets a window.
    """
    global _disk_loaded
    provider = (provider or "").strip().lower()
    model = (model or "").strip()

    if not _disk_loaded:
        _disk_loaded = True
        _load_disk_cache()

    key = _cache_key(provider, model, base_url)
    cached = _CACHE.get(key)
    if cached is not None:
        return cached

    spec: dict[str, Any] = {}
    source = "fallback"

    # 1. Live discovery -----------------------------------------------------
    # Ollama exposes per-model context_length + capabilities via /api/show.
    # KIT Toolbox (and any vLLM/OpenAI-compatible backend) reports the served
    # window as max_model_len via /v1/models — the server's TRUE limit, so we
    # use it for KIT too. Probes are cached (process + ~/.delfin) so this is a
    # one-time cost per (provider, model, endpoint).
    if allow_live and base_url:
        if provider == "ollama":
            live = _discover_ollama(base_url, model)
        else:
            live = _discover_openai_models(base_url, model)
        if live:
            spec = live
            source = "live"

    # 2. Static table (fill gaps the live probe didn't cover) ---------------
    static = _static_spec(model)
    if static:
        merged = dict(static)
        merged.update(spec)          # live values win over static
        spec = merged
        if source != "live":
            source = "static"

    # 3. Heuristic window when still unknown --------------------------------
    if "context_window" not in spec:
        spec["context_window"] = (
            _WEAK_FALLBACK_WINDOW if _is_weak(model) else _STRONG_FALLBACK_WINDOW
        )
        if source == "fallback":
            source = "heuristic"

    # Heuristic reasoning/thinking flags when neither live nor static set them.
    if "is_reasoning" not in spec and _REASONING_NAME.search(model):
        spec["is_reasoning"] = True
    if "thinking_tagged" not in spec and _THINKING_TAG_NAME.search(model):
        spec["thinking_tagged"] = True

    # Ollama: re-cap window to the num_ctx we will actually send.
    if provider == "ollama":
        spec = _apply_ollama_window(spec)

    caps = ModelCapabilities(
        model=model,
        provider=provider,
        context_window=int(spec.get("context_window") or _LEGACY_WINDOW),
        max_output_tokens=int(spec.get("max_output_tokens") or 8_192),
        supports_tools=bool(spec.get("supports_tools", True)),
        supports_vision=bool(spec.get("supports_vision", False)),
        is_reasoning=bool(spec.get("is_reasoning", False)),
        thinking_tagged=bool(spec.get("thinking_tagged", False)),
        recommended_effort=str(spec.get("recommended_effort", "medium")),
        num_ctx_override=spec.get("num_ctx_override"),
        source=source,
        discovered_at=time.time(),
    )

    _CACHE[key] = caps
    if source == "live":
        _save_disk_cache()
    return caps


def register_static(model: str, spec: dict[str, Any]) -> None:
    """Register/override an exact-name static spec (tests / experiments)."""
    _STATIC[model] = dict(spec)


# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------


def kit_recommendation(model: str, caps: "ModelCapabilities") -> str:
    """Warn-and-redirect string for a KIT model not worth the agent.

    Empty string when the model is fine. Only the strong, tool-capable,
    large-window KIT models "lohnen sich"; everything else gets a warning
    pointing at :data:`KIT_BEST_MODEL`.
    """
    if model in _KIT_RECOMMENDED:
        return ""
    if caps.supports_tools and caps.context_window >= _KIT_MIN_WINDOW \
            and not _is_weak(model):
        return ""
    reasons: list[str] = []
    if not caps.supports_tools:
        reasons.append("kein Tool-Support")
    if _is_weak(model):
        reasons.append("schwaches/kleines Modell")
    if caps.context_window < _KIT_MIN_WINDOW:
        reasons.append(f"kleines Kontextfenster ({caps.context_window})")
    why = ", ".join(reasons) or "für Agent-Aufgaben nicht empfohlen"
    return (
        f"KIT-Modell `{model}` lohnt sich für Agent-Aufgaben kaum ({why}). "
        f"Bestes KIT-Modell: `{KIT_BEST_MODEL}`."
    )


def preflight(
    provider: str,
    model: str,
    base_url: str = "",
    *,
    needs_tools: bool = True,
) -> tuple[bool, str]:
    """Cheap readiness check before a turn. Returns ``(ok, message)``.

    Contract:
      - ``(False, msg)`` → hard block (endpoint down, model not installed,
        or no tool support when the run needs tools). ``msg`` is actionable.
      - ``(True, msg)`` with non-empty ``msg`` → soft warning: the run can
        proceed but the user is steered to a better model.
      - ``(True, "")`` → all clear.
    """
    provider = (provider or "").strip().lower()
    model = (model or "").strip()

    if provider == "ollama":
        root = _ollama_root(base_url)
        installed = _ollama_installed(root)
        if installed is None:
            return False, (
                f"Ollama unter {root} nicht erreichbar — läuft `ollama serve`? "
                f"(Endpoint via OLLAMA_HOST überschreibbar.)"
            )
        # Ollama tags include the ":latest" suffix; match leniently.
        def _match(tag: str) -> bool:
            return tag == model or tag.split(":")[0] == model.split(":")[0]
        if model and not any(_match(t) for t in installed):
            avail = ", ".join(sorted(installed)[:8]) or "(keine)"
            return False, (
                f"Modell `{model}` ist in Ollama nicht installiert — "
                f"`ollama pull {model}`. Verfügbar: {avail}"
            )

    caps = resolve(provider, model, base_url)

    # Hard block: agent work with a model that can't call tools at all.
    if needs_tools and not caps.supports_tools:
        return False, (
            f"Modell `{model}` hat keine native Tool-Unterstützung und kann "
            f"daher keine Agent-Werkzeuge (bash/edit/read …) ausführen. Wähle "
            f"ein tool-fähiges Modell (z. B. `qwen2.5-coder`, `llama3.3`, "
            f"`{KIT_BEST_MODEL}`) oder nutze dieses Modell im reinen Chat-Modus."
        )

    # Soft warning: a KIT model that works but isn't worth the agent.
    if provider == "kit":
        warn = kit_recommendation(model, caps)
        if warn:
            return True, warn

    return True, ""


__all__ = [
    "ModelCapabilities",
    "resolve",
    "preflight",
    "kit_recommendation",
    "register_static",
    "clear_cache",
    "DEFAULT_NUM_CTX_CAP",
    "KIT_BEST_MODEL",
]
