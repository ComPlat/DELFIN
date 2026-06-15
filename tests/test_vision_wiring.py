"""Vision wiring: image pixels reach a vision-capable model as image_url
content; non-vision models / Claude get a text fallback."""

from __future__ import annotations

import pytest

from delfin.agent import model_capabilities as mc
from delfin.agent.engine import AgentEngine

# 1x1 transparent PNG
_PNG = bytes.fromhex(
    "89504e470d0a1a0a0000000d49484452000000010000000108060000001f15c489"
    "0000000d49444154789c626001000000050001a5f645400000000049454e44ae426082"
)


def _img(tmp_path):
    p = tmp_path / "x.png"
    p.write_bytes(_PNG)
    return str(p)


def _engine(tmp_path, monkeypatch, *, supports_vision):
    caps = mc.ModelCapabilities(
        model="qwen3-vl:4b", provider="ollama", context_window=32768,
        supports_tools=True, supports_vision=supports_vision,
        num_ctx_override=32768)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    eng = AgentEngine(repo_dir=str(tmp_path), backend="api",
                      provider="ollama", model="qwen3-vl:4b")
    eng._active_capabilities = caps
    return eng


def test_vision_model_gets_image_url_content(tmp_path, monkeypatch):
    eng = _engine(tmp_path, monkeypatch, supports_vision=True)
    msg = eng._build_user_message("what is this?", [_img(tmp_path)])
    assert isinstance(msg["content"], list)
    types = [p.get("type") for p in msg["content"]]
    assert "image_url" in types and "text" in types
    url = next(p["image_url"]["url"] for p in msg["content"]
               if p.get("type") == "image_url")
    assert url.startswith("data:image/png;base64,")


def test_non_vision_model_text_fallback(tmp_path, monkeypatch):
    eng = _engine(tmp_path, monkeypatch, supports_vision=False)
    msg = eng._build_user_message("hi", [_img(tmp_path)])
    assert msg["content"] == "hi"


def test_no_images_plain_text(tmp_path, monkeypatch):
    eng = _engine(tmp_path, monkeypatch, supports_vision=True)
    assert eng._build_user_message("hi", None)["content"] == "hi"


def test_claude_provider_text_fallback(tmp_path, monkeypatch):
    # Claude uses a different image format → must not get OpenAI image_url.
    eng = _engine(tmp_path, monkeypatch, supports_vision=True)
    eng.provider = "claude"
    msg = eng._build_user_message("hi", [_img(tmp_path)])
    assert msg["content"] == "hi"
