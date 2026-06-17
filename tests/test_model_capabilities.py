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
    assert "erreichbar" in msg


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
    assert "Tool-Unterst" in msg


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
