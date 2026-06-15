"""Smoke tests for the multimodal image-input layer.

Covers MIME validation, size cap, vision-model detection, and the
plain-text fallback when the active model isn't vision-capable.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.image_input import (
    ImageAttachment, ImageError, load_image,
    model_supports_vision, to_openai_content,
)


# 1x1 transparent PNG — just enough for the file-format checker.
_PNG_BYTES = (
    b"\x89PNG\r\n\x1a\n"
    b"\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x06\x00\x00\x00"
    b"\x1f\x15\xc4\x89\x00\x00\x00\rIDATx\x9cc\xfc\xcf\xc0\x00\x00\x00\x03"
    b"\x00\x01<>\xa6\xeb\x00\x00\x00\x00IEND\xaeB`\x82"
)


def test_load_image_accepts_png(tmp_path):
    p = tmp_path / "tiny.png"
    p.write_bytes(_PNG_BYTES)
    img = load_image(p)
    assert isinstance(img, ImageAttachment)
    assert img.mime == "image/png"
    assert img.size_bytes() == len(_PNG_BYTES)
    assert img.source_path == p


def test_load_image_rejects_missing_file(tmp_path):
    with pytest.raises(ImageError):
        load_image(tmp_path / "nope.png")


def test_load_image_rejects_unsupported_mime(tmp_path):
    p = tmp_path / "data.txt"
    p.write_text("not an image")
    with pytest.raises(ImageError) as exc:
        load_image(p)
    assert "unsupported MIME" in str(exc.value)


def test_load_image_rejects_oversized(tmp_path, monkeypatch):
    """Cap is 8 MB raw; we monkey-patch the limit so the test stays fast."""
    import delfin.agent.image_input as mod
    monkeypatch.setattr(mod, "_MAX_BYTES", 16)
    p = tmp_path / "big.png"
    p.write_bytes(_PNG_BYTES)  # 67 bytes > 16
    with pytest.raises(ImageError) as exc:
        load_image(p)
    assert "too large" in str(exc.value)


def test_data_uri_round_trip(tmp_path):
    p = tmp_path / "a.png"
    p.write_bytes(_PNG_BYTES)
    img = load_image(p)
    uri = img.data_uri()
    assert uri.startswith("data:image/png;base64,")


@pytest.mark.parametrize(
    "model, expected",
    [
        ("gpt-5.4", True),
        ("gpt-4.1-mini", True),
        ("o4-mini", True),
        ("claude-sonnet-4-6", True),
        ("claude-opus-4-7", True),
        ("haiku", False),  # bare "haiku" without "claude-" prefix
        ("gemini-2-pro", True),
        ("text-only-model", False),
        ("", False),
    ],
)
def test_model_supports_vision(model, expected):
    assert model_supports_vision(model) is expected


def test_to_openai_content_text_only_when_no_images():
    assert to_openai_content("hi") == "hi"


def test_to_openai_content_builds_multimodal_parts_for_vision_model(tmp_path):
    p = tmp_path / "x.png"
    p.write_bytes(_PNG_BYTES)
    img = load_image(p)
    content = to_openai_content(
        "what do you see?", images=[img], model="gpt-5.4"
    )
    assert isinstance(content, list)
    assert content[0]["type"] == "text"
    assert content[0]["text"] == "what do you see?"
    assert content[1]["type"] == "image_url"
    assert content[1]["image_url"]["url"].startswith("data:image/png;base64,")


def test_to_openai_content_falls_back_to_text_for_non_vision_model(tmp_path):
    p = tmp_path / "x.png"
    p.write_bytes(_PNG_BYTES)
    img = load_image(p)
    content = to_openai_content(
        "what do you see?", images=[img], model="text-only-model"
    )
    assert isinstance(content, str)
    assert "image attachments suppressed" in content
    assert "x.png" in content


def test_to_openai_content_force_bypasses_capability_check(tmp_path):
    p = tmp_path / "x.png"
    p.write_bytes(_PNG_BYTES)
    img = load_image(p)
    content = to_openai_content(
        "look", images=[img], model="text-only-model", force=True
    )
    assert isinstance(content, list)


def test_to_openai_content_rejects_too_many_images(tmp_path):
    p = tmp_path / "x.png"
    p.write_bytes(_PNG_BYTES)
    img = load_image(p)
    with pytest.raises(ImageError) as exc:
        to_openai_content("hi", images=[img] * 11, model="gpt-5.4")
    assert "too many images" in str(exc.value)


# ---------------------------------------------------------------------------
# Capability-aware vision detection (model_supports_vision(model, caps))
# ---------------------------------------------------------------------------


def _vcaps(supports_vision: bool):
    from delfin.agent.model_capabilities import ModelCapabilities
    return ModelCapabilities(model="m", provider="ollama",
                             context_window=32_768,
                             supports_vision=supports_vision)


def test_vision_via_caps_overrides_name_allowlist():
    # A local model not on the static allow-list, but discovered vision-capable.
    assert model_supports_vision("some-local-vlm:7b", _vcaps(True)) is True


def test_vision_caps_false_blocks_even_known_name():
    assert model_supports_vision("llava:7b", _vcaps(False)) is False


def test_vision_caps_none_falls_back_to_name():
    assert model_supports_vision("llava:7b") is True
    assert model_supports_vision("llama3.1:8b") is False
