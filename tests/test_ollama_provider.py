"""Tests for the Ollama (and OpenAI-compatible local server) provider.

The Ollama agent reuses the same OpenAIClient agentic-while-loop +
KIT-style sandbox/hooks/failure-log stack as the cloud providers —
only the base_url differs and no api key is required.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# create_client builds an OpenAIClient pointed at the right endpoint
# ---------------------------------------------------------------------------


def test_create_client_ollama_default_localhost(monkeypatch):
    """No env var → localhost:11434/v1 (the canonical Ollama default)."""
    monkeypatch.delenv("OLLAMA_HOST", raising=False)
    monkeypatch.delenv("OLLAMA_BASE_URL", raising=False)
    from delfin.agent.api_client import create_client, OpenAIClient
    client = create_client(backend="api", provider="ollama", model="qwen3-coder:32b")
    assert isinstance(client, OpenAIClient)
    base = str(client.client.base_url)
    assert "localhost:11434" in base
    assert base.rstrip("/").endswith("/v1")


def test_create_client_ollama_honours_OLLAMA_HOST(monkeypatch):
    monkeypatch.setenv("OLLAMA_HOST", "http://gpu-server:11434")
    from delfin.agent.api_client import create_client
    client = create_client(backend="api", provider="ollama", model="llama3.3:70b")
    assert "gpu-server" in str(client.client.base_url)


def test_create_client_ollama_honours_OLLAMA_BASE_URL_as_fallback(monkeypatch):
    monkeypatch.delenv("OLLAMA_HOST", raising=False)
    monkeypatch.setenv("OLLAMA_BASE_URL", "http://lmstudio.local:1234/v1")
    from delfin.agent.api_client import create_client
    client = create_client(backend="api", provider="ollama", model="qwen2.5-coder:14b")
    assert "lmstudio.local" in str(client.client.base_url)


def test_create_client_ollama_appends_v1_when_missing(monkeypatch):
    monkeypatch.setenv("OLLAMA_HOST", "http://localhost:8080")
    from delfin.agent.api_client import create_client
    client = create_client(backend="api", provider="ollama", model="x")
    assert str(client.client.base_url).rstrip("/").endswith("/v1")


def test_create_client_ollama_explicit_host_via_api_key_arg(monkeypatch):
    """``api_key='http://...'`` is interpreted as a host override —
    handy for the headless CLI's ``--api-key`` channel."""
    monkeypatch.delenv("OLLAMA_HOST", raising=False)
    monkeypatch.delenv("OLLAMA_BASE_URL", raising=False)
    from delfin.agent.api_client import create_client
    client = create_client(
        backend="api", provider="ollama",
        api_key="http://remote-gpu:11434",
        model="qwen3-coder:32b",
    )
    assert "remote-gpu" in str(client.client.base_url)


def test_create_client_ollama_attaches_kit_perms_sandbox(monkeypatch, tmp_path):
    """Local-model agent gets the same sandbox / allow-list / deny-list
    as the KIT-Toolbox cloud provider so a malicious prompt can't escape
    just because we're running locally."""
    monkeypatch.delenv("OLLAMA_HOST", raising=False)
    from delfin.agent.api_client import create_client, KitToolPermissions
    client = create_client(
        backend="api", provider="ollama",
        model="qwen3-coder:32b", cwd=str(tmp_path),
    )
    assert client._permissions is not None
    assert isinstance(client._permissions, KitToolPermissions)
    assert client._permissions.workspace == tmp_path.resolve()


def test_create_client_ollama_no_api_key_required(monkeypatch):
    """Local Ollama doesn't need auth; the openai-python client must
    still construct cleanly without OPENAI_API_KEY set."""
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.delenv("OLLAMA_HOST", raising=False)
    from delfin.agent.api_client import create_client
    # Should not raise
    client = create_client(backend="api", provider="ollama", model="x")
    assert client is not None


# ---------------------------------------------------------------------------
# Provider dropdown wiring
# ---------------------------------------------------------------------------


def test_ollama_model_fallback_present_in_dashboard():
    """The Ollama fallback list must exist in the dashboard so the
    dropdown isn't empty when /api/tags can't be reached."""
    src = Path(__file__).resolve().parent.parent / "delfin" / "dashboard" / "tab_agent.py"
    text = src.read_text(encoding="utf-8")
    # Locate the _PROVIDER_MODELS_FALLBACK block + verify ollama is in it
    idx = text.find('_PROVIDER_MODELS_FALLBACK = {')
    assert idx > 0
    block = text[idx: idx + 4000]
    assert '"ollama":' in block
    assert "qwen3-coder" in block or "qwen2.5-coder" in block


def test_ollama_pricing_treated_as_free():
    """Local-model models must not be charged at cloud prices in the
    /cost / /usage output."""
    src = Path(__file__).resolve().parent.parent / "delfin" / "dashboard" / "tab_agent.py"
    text = src.read_text(encoding="utf-8")
    # The free-pricing branch must exist
    assert 'provider_dropdown.value == "ollama"' in text
    # Free pricing means input=0, output=0
    free_block_idx = text.find('provider_dropdown.value == "ollama"')
    snippet = text[free_block_idx: free_block_idx + 500]
    assert '"input": 0.0' in snippet
    assert '"output": 0.0' in snippet


# ---------------------------------------------------------------------------
# num_ctx is sent ONLY for Ollama (the "full potential" fix)
# ---------------------------------------------------------------------------


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Delta:
    content = "ok"
    tool_calls = None


class _Choice:
    delta = _Delta()
    finish_reason = "stop"


class _Chunk:
    usage = _Usage()
    choices = [_Choice()]


class _FakeStream:
    def __iter__(self):
        return iter([_Chunk()])

    def close(self):
        pass


def _capture_create(client):
    """Replace the model's create() with a capturing stub; return the dict
    that will hold the kwargs of the last call."""
    captured: dict = {}

    def _fake_create(**kwargs):
        captured.clear()
        captured.update(kwargs)
        return _FakeStream()

    client.client.chat.completions.create = _fake_create
    return captured


def _fixed_caps(monkeypatch, *, provider, num_ctx):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(
        model="m", provider=provider, context_window=num_ctx or 200_000,
        supports_tools=True, num_ctx_override=num_ctx, source="static",
    )
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)


def test_stream_message_sends_num_ctx_for_ollama(monkeypatch, tmp_path):
    from delfin.agent.api_client import create_client
    _fixed_caps(monkeypatch, provider="ollama", num_ctx=4096)
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))
    captured = _capture_create(client)
    list(client.stream_message("sys", [{"role": "user", "content": "hi"}],
                               max_tokens=100))
    assert "extra_body" in captured
    assert captured["extra_body"]["options"]["num_ctx"] == 4096
    # Ollama must NOT receive reasoning_effort / max_completion_tokens.
    assert "reasoning_effort" not in captured
    assert "max_tokens" in captured


def test_stream_message_no_num_ctx_for_kit(monkeypatch, tmp_path):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "dummy-key")
    from delfin.agent.api_client import create_client
    _fixed_caps(monkeypatch, provider="kit", num_ctx=None)
    client = create_client(backend="api", provider="kit",
                           model="kit.qwen3.5-397b-A17b", cwd=str(tmp_path))
    captured = _capture_create(client)
    list(client.stream_message("sys", [{"role": "user", "content": "hi"}],
                               max_tokens=100))
    assert "extra_body" not in captured
