"""Tests for the model-capability + context-adaptivity layer.

No live Ollama / KIT endpoint is required: ``urllib.request.urlopen`` is
monkeypatched with canned ``/api/show``, ``/api/tags`` and ``/v1/models``
payloads. Covers the core "full potential" guarantees:

  - Ollama context window is discovered and num_ctx is capped (and raisable).
  - context_window == num_ctx_override for Ollama (compaction matches reality).
  - KIT/vLLM window is discovered live from /v1/models max_model_len.
  - Tool / vision / reasoning flags resolve from capabilities.
  - Graceful static/heuristic fallback when the network is unreachable.
  - Preflight surfaces actionable messages (down / not pulled / no tools /
    weak-KIT recommendation).
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import model_capabilities as mc


# ---------------------------------------------------------------------------
# Fake HTTP plumbing
# ---------------------------------------------------------------------------


class _FakeResp:
    def __init__(self, payload):
        self._b = json.dumps(payload).encode("utf-8")

    def read(self):
        return self._b

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _urlopen_router(routes):
    """Return a fake urlopen dispatching by URL substring.

    ``routes`` maps a URL fragment → payload dict (or an Exception class to
    raise). Unmatched URLs raise to mimic an unreachable endpoint.
    """
    def _fake(req, timeout=None):
        url = getattr(req, "full_url", None) or str(req)
        for frag, payload in routes.items():
            if frag in url:
                if isinstance(payload, type) and issubclass(payload, Exception):
                    raise payload("boom")
                return _FakeResp(payload)
        raise OSError(f"no route for {url}")
    return _fake


@pytest.fixture(autouse=True)
def _isolate(monkeypatch, tmp_path):
    """Isolate cache + pin the num_ctx cap so tests are deterministic."""
    monkeypatch.setattr(mc, "_CACHE_PATH", tmp_path / "caps_cache.json")
    monkeypatch.setattr(mc, "_disk_loaded", True)   # skip disk load
    mc._CACHE.clear()
    monkeypatch.setattr(mc, "_configured_num_ctx_cap", lambda: 32_768)
    yield
    mc._CACHE.clear()


_OLLAMA_BASE = "http://localhost:11434/v1"
_KIT_BASE = "https://ki-toolbox.scc.kit.edu/api/v1"


def _show_payload(ctx_len, caps):
    return {"model_info": {"qwen3.context_length": ctx_len},
            "capabilities": list(caps)}


# ---------------------------------------------------------------------------
# Ollama live discovery
# ---------------------------------------------------------------------------


def test_ollama_live_discovery_window_and_caps(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({"/api/show": _show_payload(32_768, ["tools", "vision"])}),
    )
    caps = mc.resolve("ollama", "qwen3:8b", _OLLAMA_BASE)
    assert caps.source == "live"
    assert caps.context_window == 32_768
    assert caps.num_ctx_override == 32_768          # == context_window (invariant)
    assert caps.supports_tools is True
    assert caps.supports_vision is True


def test_ollama_num_ctx_is_capped_below_model_max(monkeypatch):
    # Model advertises 131k but the default cap (32k) must clamp num_ctx,
    # and the compaction window must equal that cap — never the 131k max.
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({"/api/show": _show_payload(131_072, ["tools"])}),
    )
    caps = mc.resolve("ollama", "qwen2.5-coder:32b", _OLLAMA_BASE)
    assert caps.num_ctx_override == 32_768
    assert caps.context_window == 32_768            # invariant: not 131_072


def test_ollama_num_ctx_cap_is_raisable(monkeypatch):
    monkeypatch.setattr(mc, "_configured_num_ctx_cap", lambda: 131_072)
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({"/api/show": _show_payload(131_072, ["tools"])}),
    )
    caps = mc.resolve("ollama", "qwen2.5-coder:32b", _OLLAMA_BASE)
    assert caps.num_ctx_override == 131_072
    assert caps.context_window == 131_072


def test_ollama_no_tools_capability(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({"/api/show": _show_payload(8_192, ["completion"])}),
    )
    caps = mc.resolve("ollama", "tinyllama:1b", _OLLAMA_BASE)
    assert caps.supports_tools is False


def test_ollama_thinking_capability(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({"/api/show": _show_payload(65_536,
                                                    ["tools", "thinking"])}),
    )
    caps = mc.resolve("ollama", "qwq:32b", _OLLAMA_BASE)
    assert caps.is_reasoning is True
    assert caps.thinking_tagged is True


# ---------------------------------------------------------------------------
# Fallback when the network is unreachable
# ---------------------------------------------------------------------------


def test_static_fallback_on_network_error(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/api/show": OSError}),
    )
    caps = mc.resolve("ollama", "llama3.1:8b", _OLLAMA_BASE)
    assert caps.source in {"static", "heuristic"}
    assert caps.context_window > 0
    assert caps.num_ctx_override == 32_768          # still re-capped for ollama


def test_heuristic_reasoning_flags_without_network():
    # No base_url → no live probe; deepseek-r1 known reasoning/thinking.
    caps = mc.resolve("ollama", "deepseek-r1:7b", "")
    assert caps.is_reasoning is True
    assert caps.thinking_tagged is True


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------


def test_cache_hit_survives_later_network_failure(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({"/api/show": _show_payload(16_384, ["tools"])}),
    )
    first = mc.resolve("ollama", "phi3:mini", _OLLAMA_BASE)
    assert first.source == "live" and first.context_window == 16_384
    # Now the endpoint dies — the cached value must still be returned.
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/api/show": OSError}),
    )
    second = mc.resolve("ollama", "phi3:mini", _OLLAMA_BASE)
    assert second.context_window == 16_384 and second.source == "live"


# ---------------------------------------------------------------------------
# KIT / vLLM live window from /v1/models max_model_len
# ---------------------------------------------------------------------------


def test_kit_live_window_from_v1_models(monkeypatch):
    payload = {"data": [
        {"id": "qwen3.5-397b-A17b", "object": "model", "max_model_len": 262_144},
        {"id": "gemma4-31b-it", "object": "model", "max_model_len": 131_072},
    ]}
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/v1/models": payload}),
    )
    caps = mc.resolve("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE)
    assert caps.source == "live"
    assert caps.context_window == 262_144           # the served max, not static 128k
    assert caps.num_ctx_override is None            # cloud honours context server-side


def test_kit_static_fallback_when_models_unreachable(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/v1/models": OSError}),
    )
    caps = mc.resolve("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE)
    assert caps.context_window == 128_000           # curated static value
    assert caps.source == "static"


def test_kit_live_window_requires_auth_header(monkeypatch):
    """KIT /v1/models rejects an unauthenticated probe (401) and falls back to
    static; with the provider key threaded through, the Bearer header is sent
    and the server's TRUE window is read live."""
    payload = {"data": [
        {"id": "qwen3.5-397b-A17b", "object": "model", "max_model_len": 262_144},
    ]}
    seen = {}

    def _fake(req, timeout=None):
        auth = req.headers.get("Authorization")
        seen["auth"] = auth
        if "/v1/models" not in (getattr(req, "full_url", "") or str(req)):
            raise OSError("no route")
        if not auth:
            raise OSError("401 Unauthorized")       # mimic KIT's rejection
        return _FakeResp(payload)

    monkeypatch.setattr(mc.urllib.request, "urlopen", _fake)

    # No key → unauthenticated probe is rejected → static fallback.
    caps0 = mc.resolve("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE)
    assert caps0.source == "static"
    assert caps0.context_window == 128_000

    # With the key → Bearer header sent, live window discovered (distinct cache
    # entry, so the earlier key-less miss is not served).
    caps1 = mc.resolve("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE, api_key="SECRET")
    assert seen["auth"] == "Bearer SECRET"
    assert caps1.source == "live"
    assert caps1.context_window == 262_144


def test_kit_vision_capability_parsed_live(monkeypatch):
    """KIT exposes a per-model vision flag under info.meta.capabilities — read it
    so a multimodal model (qwen3.5-397b, gemma4) is recognised as vision-capable
    instead of falling back to the name heuristic (which would say no)."""
    def _payload(vis):
        return {"data": [{"id": "kit.m", "info": {"id": "kit.m", "meta": {
            "description": "… Context length ≈ 256K",
            "capabilities": {"vision": vis, "file_context": True}}}}]}
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/v1/models": _payload(True)}))
    assert mc.resolve("kit", "kit.m", _KIT_BASE, api_key="K").supports_vision is True
    mc._CACHE.clear()
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/v1/models": _payload(False)}))
    assert mc.resolve("kit", "kit.m", _KIT_BASE, api_key="K").supports_vision is False


def test_nonchat_reason_filters_modality_models():
    """Embedding / reranker / speech / image models (listed by KIT alongside
    chat models, but they 400 on the chat endpoint) get a reason; real chat
    models — including the working standard-* aliases — return None."""
    nonchat = mc.nonchat_reason
    assert nonchat("kit.qwen3-embedding-8b") == "embedding model"
    assert nonchat("kit.qwen3-reranker-8b") == "reranker"
    assert nonchat("kit.whisper-large-v3") == "speech-to-text"
    assert nonchat("kit.voxtral-4b-tts-2603") == "speech model"
    assert nonchat("kit.flux.2-dev") == "image generation"
    for chat in ("azure.gpt-5.4", "kit.qwen3.5-397b-A17b", "kit.gpt-oss-120b",
                 "kit.gemma4-31b-it", "kit.minimax-m2.7-229b",
                 "kit.mistral-small-4-119b-a8b", "azure.o3",
                 "standard-extern", "standard-local", ""):
        assert nonchat(chat) is None, f"{chat!r} wrongly flagged non-chat"


def test_kit_window_from_open_webui_description(monkeypatch):
    """KIT Toolbox is an Open-WebUI proxy: no numeric max_model_len — the
    window lives as prose in info.meta.description ("Context length ≈ 256K").
    Parse it as the live signal (the real KIT /v1/models shape)."""
    payload = {"data": [
        {"id": "kit.qwen3.5-397b-A17b",
         "info": {"id": "kit.qwen3.5-397b-A17b",
                  "meta": {"description":
                           "MoE multimodal model • Context length ≈ 256K\n"
                           "Host: KIT • Origin: Alibaba"}}},
    ]}
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/v1/models": payload}),
    )
    caps = mc.resolve("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE, api_key="K")
    assert caps.source == "live"
    assert caps.context_window == 256 * 1024        # 262144, parsed from prose


# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------


def test_preflight_ollama_unreachable(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/api/tags": OSError}),
    )
    ok, msg = mc.preflight("ollama", "qwen3-coder:32b", _OLLAMA_BASE)
    assert ok is False
    assert "not reachable" in msg


def test_preflight_model_not_installed(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({"/api/tags": {"models": [{"name": "llama3.1:8b"}]}}),
    )
    ok, msg = mc.preflight("ollama", "qwen3-coder:32b", _OLLAMA_BASE)
    assert ok is False
    assert "pull" in msg


def test_preflight_no_tool_support_blocks(monkeypatch):
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({
            "/api/tags": {"models": [{"name": "tinyllama:1b"}]},
            "/api/show": _show_payload(4_096, ["completion"]),
        }),
    )
    ok, msg = mc.preflight("ollama", "tinyllama:1b", _OLLAMA_BASE)
    assert ok is False
    assert "no native tool support" in msg


def test_preflight_kit_weak_model_warns_and_recommends(monkeypatch):
    # No /v1/models route → static/heuristic; a weak-by-name KIT model.
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({}),
    )
    ok, msg = mc.preflight("kit", "kit.gemma2-2b-it", _KIT_BASE)
    assert ok is True                                # soft warning, not a block
    assert mc.KIT_BEST_MODEL in msg


def test_preflight_kit_best_model_clean(monkeypatch):
    payload = {"data": [{"id": "qwen3.5-397b-A17b", "max_model_len": 262_144}]}
    monkeypatch.setattr(
        mc.urllib.request, "urlopen", _urlopen_router({"/v1/models": payload}),
    )
    ok, msg = mc.preflight("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE)
    assert ok is True
    assert msg == ""


# ---------------------------------------------------------------------------
# Disk-cache schema versioning (bug 140214/172931: flush pre-v2 live entries
# that wrongly cached supports_vision=False for vision-capable KIT models)
# ---------------------------------------------------------------------------


def test_disk_cache_version_envelope_roundtrips(monkeypatch, tmp_path):
    p = tmp_path / "caps.json"
    monkeypatch.setattr(mc, "_CACHE_PATH", p)
    mc._CACHE.clear()
    key = mc._cache_key("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE, True)
    mc._CACHE[key] = mc.ModelCapabilities(
        model="kit.qwen3.5-397b-A17b", provider="kit",
        context_window=262_144, supports_vision=True, source="live",
        discovered_at=time.time(),
    )
    mc._save_disk_cache()
    raw = json.loads(p.read_text())
    assert raw["_v"] == mc._CACHE_VERSION and "entries" in raw
    mc._CACHE.clear()
    mc._load_disk_cache()
    assert key in mc._CACHE and mc._CACHE[key].supports_vision is True


def test_disk_cache_old_flat_format_is_discarded(monkeypatch, tmp_path):
    # Pre-v2 flat file (no "_v" envelope) holding a stale vision=False live
    # entry — must be dropped wholesale so a fresh live resolve can fix vision.
    p = tmp_path / "caps.json"
    monkeypatch.setattr(mc, "_CACHE_PATH", p)
    key = mc._cache_key("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE, True)
    stale = {"model": "kit.qwen3.5-397b-A17b", "provider": "kit",
             "context_window": 262_144, "supports_vision": False,
             "source": "live", "discovered_at": time.time()}
    p.write_text(json.dumps({key: stale}))      # old flat shape, no "_v"
    mc._CACHE.clear()
    mc._load_disk_cache()
    assert key not in mc._CACHE


# ---------------------------------------------------------------------------
# KIT-served model ids (observed on the live /v1/models listing, 2026-07-29):
# every chat-capable id must resolve to a curated static entry — never fall
# through to the silent heuristic default that mis-gates tools/windows.
# ---------------------------------------------------------------------------

_KIT_SERVED_CHAT_IDS = (
    "azure.gpt-5",
    "azure.gpt-5-mini",
    "azure.gpt-5-nano",
    "azure.gpt-5.1",
    "azure.gpt-5.4",
    "azure.gpt-5.5",
    "azure.gpt-5.6-luna",
    "azure.gpt-5.6-sol",
    "azure.gpt-5.6-terra",
    "google.claude-fable-5",
    "google.claude-haiku-4.5",
    "google.claude-opus-4.8",
    "google.claude-sonnet-4.6",
    "google.claude-sonnet-5",
    "google.gemini-2.5-flash",
    "google.gemini-2.5-flash-lite",
    "google.gemini-2.5-pro",
    "google.gemini-3.1-flash-lite",
    "google.gemini-3.5-flash",
    "kit.gemma4-31b-it",
    "kit.gpt-oss-120b",
    "kit.minimax-m2.7-229b",
    "kit.mistral-small-4-119b-a8b",
    "kit.qwen3.5-397b-A17b",
    "standard-extern",
    "standard-local",
)

_KIT_SERVED_NONCHAT_IDS = (
    "kit.flux.2-dev",
    "kit.qwen3-embedding-8b",
    "kit.qwen3-reranker-8b",
    "kit.voxtral-4b-tts-2603",
    "kit.whisper-large-v3",
)


@pytest.mark.parametrize("model", _KIT_SERVED_CHAT_IDS)
def test_kit_served_chat_id_resolves_statically(model):
    # base_url="" → no live probe; the static registry alone must know the id.
    caps = mc.resolve("kit", model, "")
    assert caps.source == "static", f"{model} fell through to {caps.source}"
    assert caps.context_window >= 32_000
    assert caps.supports_tools is True


@pytest.mark.parametrize("model", _KIT_SERVED_NONCHAT_IDS)
def test_kit_served_nonchat_id_is_flagged(model):
    assert mc.nonchat_reason(model) is not None


def test_kit_served_chat_ids_pass_preflight_without_warning():
    # The strong KIT-hosted chat models must not trip the weak-KIT warning.
    for model in ("kit.qwen3.5-397b-A17b", "kit.gpt-oss-120b",
                  "kit.minimax-m2.7-229b", "kit.mistral-small-4-119b-a8b"):
        ok, msg = mc.preflight("kit", model, "")
        assert ok is True
        assert msg == "", f"unexpected warning for {model}: {msg}"


# ---------------------------------------------------------------------------
# Resolution bugs (each failed before the corresponding registry fix)
# ---------------------------------------------------------------------------


def test_google_gateway_prefix_is_stripped():
    # Served ids carry a google. prefix; before the fix only kit./azure. were
    # stripped, so these fell to the heuristic (32k, no vision).
    caps = mc.resolve("kit", "google.gemini-2.5-pro", "")
    assert caps.source == "static"
    assert caps.supports_vision is True
    assert caps.context_window > 100_000


def test_static_lookup_is_case_insensitive():
    # The served id has a mixed-case active-parameter suffix (…-A17b); a
    # lowercased config value must still hit the curated exact entry.
    caps = mc.resolve("kit", "kit.qwen3.5-397b-a17b", "")
    assert caps.source == "static"
    assert caps.context_window == 128_000           # exact entry, not a prefix


def test_gpt_oss_is_reasoning_for_token_floor():
    # gpt-oss is a reasoning line: without is_reasoning the max_tokens floor
    # never applies and small budgets yield empty replies mid-think.
    assert mc.resolve("kit", "kit.gpt-oss-120b", "").is_reasoning is True
    assert mc.resolve("openai", "gpt-oss-20b", "").is_reasoning is True


def test_deepseek_r1_static_has_no_native_tools():
    # R1 line: reasoning without native tool calling. The static entry must
    # say so — advertising tools to it makes the backend reject the request.
    caps = mc.resolve("ollama", "deepseek-r1:32b", "")
    assert caps.supports_tools is False
    assert caps.is_reasoning is True


def test_deepseek_r1_live_tools_capability_overrides_static(monkeypatch):
    # A live probe that reports tool support wins over the conservative
    # static claim (newer builds may add native tool calling).
    monkeypatch.setattr(
        mc.urllib.request, "urlopen",
        _urlopen_router({"/api/show": _show_payload(65_536,
                                                    ["tools", "thinking"])}),
    )
    caps = mc.resolve("ollama", "deepseek-r1:32b", _OLLAMA_BASE)
    assert caps.source == "live"
    assert caps.supports_tools is True


def test_qwen3_coder_is_not_thinking_tagged():
    # The coder line is non-thinking; the qwen3 name heuristic must not mark
    # it thinking-tagged (static entry pins the flags explicitly).
    caps = mc.resolve("ollama", "qwen3-coder:32b", "")
    assert caps.thinking_tagged is False
    assert caps.is_reasoning is False


# ---------------------------------------------------------------------------
# New family entries
# ---------------------------------------------------------------------------


def test_kimi_k2_family():
    caps = mc.resolve("openai", "kimi-k2:1t", "")
    assert caps.source == "static"
    assert caps.supports_tools is True
    assert caps.context_window == 131_072
    thinking = mc.resolve("openai", "kimi-k2-thinking", "")
    assert thinking.is_reasoning is True
    assert thinking.thinking_tagged is True
    assert thinking.context_window == 262_144


def test_deepseek_v3_family_keeps_tools():
    caps = mc.resolve("openai", "deepseek-v3.1", "")
    assert caps.source == "static"
    assert caps.supports_tools is True
    assert caps.is_reasoning is False


def test_glm_family_windows():
    assert mc.resolve("openai", "glm-4.6", "").context_window == 200_000
    assert mc.resolve("openai", "glm-4.5-air", "").context_window == 131_072


def test_llama_hyphenated_and_v4_ids():
    assert mc.resolve("openai", "llama-3.3-70b-instruct",
                      "").context_window == 131_072
    assert mc.resolve("openai", "llama-4-scout-17b-16e-instruct",
                      "").context_window == 131_072
    assert mc.resolve("openai", "llama4:16x17b", "").context_window == 131_072


def test_mistral_family_specific_beats_broad():
    # Longest-prefix: mistral-small (128k line) must not be shadowed by the
    # broad 32k "mistral" floor.
    assert mc.resolve("openai", "mistral-small-3.2", "").context_window \
        == 131_072
    assert mc.resolve("openai", "mixtral:8x7b", "").context_window == 32_768
    assert mc.resolve("openai", "mistral:7b", "").context_window == 32_768


def test_minimax_m2_reasoning_flags():
    caps = mc.resolve("kit", "kit.minimax-m2.7-229b", "")
    assert caps.is_reasoning is True
    assert caps.thinking_tagged is True
    generic = mc.resolve("openai", "minimax-m2", "")
    assert generic.is_reasoning is True and generic.thinking_tagged is True


def test_gemma_generation_windows_differ():
    assert mc.resolve("openai", "gemma3-27b-it", "").context_window == 131_072
    assert mc.resolve("openai", "gemma2-9b-it", "").context_window == 8_192


def test_standard_aliases_resolve_conservatively():
    for alias in ("standard-extern", "standard-local"):
        caps = mc.resolve("kit", alias, "")
        assert caps.source == "static"
        assert caps.context_window == 32_768
        assert caps.supports_tools is True


# ---------------------------------------------------------------------------
# Provenance note + unchanged fallback
# ---------------------------------------------------------------------------


def test_note_field_carries_registry_provenance():
    caps = mc.resolve("kit", "kit.qwen3.5-397b-A17b", "")
    assert "endpoint-verified" in caps.note
    family = mc.resolve("openai", "kimi-k2:1t", "")
    assert "family-knowledge" in family.note


def test_unknown_model_fallback_unchanged():
    caps = mc.resolve("kit", "totally-unknown-model-z", "")
    assert caps.source == "heuristic"
    assert caps.context_window == 32_768            # strong fallback window
    assert caps.supports_tools is True              # default stays permissive
    assert caps.note == ""


def test_disk_cache_version_mismatch_is_discarded(monkeypatch, tmp_path):
    p = tmp_path / "caps.json"
    monkeypatch.setattr(mc, "_CACHE_PATH", p)
    key = mc._cache_key("kit", "kit.qwen3.5-397b-A17b", _KIT_BASE, True)
    rec = {"model": "kit.qwen3.5-397b-A17b", "provider": "kit",
           "context_window": 262_144, "supports_vision": True,
           "source": "live", "discovered_at": time.time()}
    p.write_text(json.dumps({"_v": mc._CACHE_VERSION - 1, "entries": {key: rec}}))
    mc._CACHE.clear()
    mc._load_disk_cache()
    assert key not in mc._CACHE
