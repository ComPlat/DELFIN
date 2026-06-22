"""Multimodal image input helpers.

Some KIT-Toolbox-hosted models accept image content. This module
provides a thin layer that:

  - validates an image (size cap, MIME type)
  - encodes it as base64 with a data: URI
  - converts to OpenAI's chat-completions ``image_url`` content shape
  - falls back gracefully (returns a text-only message) when the
    model is known to NOT support vision

Model capability detection is best-effort: a small allow-list lives
in ``_VISION_CAPABLE_PATTERNS``. Unknown models get an OPT-IN
behaviour — they receive the image only when the caller passes
``force=True``.

Hard caps:

  - 8 MB raw, 12 MB base64
  - PNG / JPEG / WebP / GIF
  - max 10 images per message (the API will reject more anyway)
"""

from __future__ import annotations

import base64
import mimetypes
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

_MAX_BYTES = 8 * 1024 * 1024
_ALLOWED_MIMES = {"image/png", "image/jpeg", "image/webp", "image/gif"}
_MAX_PER_MESSAGE = 10

# Patterns matched against model name (case-insensitive). Models we
# know support vision today; conservative on purpose.
_VISION_CAPABLE_PATTERNS = (
    r"gpt-4\.?\d*",            # GPT-4 family (vision enabled in most variants)
    r"gpt-5",                  # GPT-5 family
    r"o4-mini",
    r"claude-(sonnet|opus|haiku)",
    r"gemini-1\.5",
    r"gemini-2",
    r"llava",
    r"qwen2-vl",
)


@dataclass
class ImageAttachment:
    bytes_: bytes
    mime: str
    source_path: Path | None = None

    def data_uri(self) -> str:
        b64 = base64.b64encode(self.bytes_).decode("ascii")
        return f"data:{self.mime};base64,{b64}"

    def size_bytes(self) -> int:
        return len(self.bytes_)


class ImageError(ValueError):
    pass


def model_supports_vision(model: str) -> bool:
    if not model:
        return False
    m = model.lower()
    return any(re.search(p, m) for p in _VISION_CAPABLE_PATTERNS)


def load_image(path: Path | str) -> ImageAttachment:
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise ImageError(f"image not found: {p}")
    mime, _ = mimetypes.guess_type(str(p))
    if mime is None or mime not in _ALLOWED_MIMES:
        raise ImageError(
            f"unsupported MIME type {mime!r} for {p.name} — "
            f"allowed: {sorted(_ALLOWED_MIMES)}"
        )
    raw = p.read_bytes()
    if len(raw) > _MAX_BYTES:
        raise ImageError(
            f"image too large ({len(raw):,} > {_MAX_BYTES:,} B)"
        )
    return ImageAttachment(bytes_=raw, mime=mime, source_path=p)


def to_openai_content(
    text: str,
    images: list[ImageAttachment] | None = None,
    *,
    model: str = "",
    force: bool = False,
) -> Any:
    """Build a chat-completions ``content`` value for a multimodal turn.

    Returns either a plain string (no images / unsupported model) or a
    list of content parts in OpenAI's mixed-content format.
    """
    images = images or []
    if not images:
        return text
    if len(images) > _MAX_PER_MESSAGE:
        raise ImageError(
            f"too many images ({len(images)} > {_MAX_PER_MESSAGE})"
        )
    if not (force or model_supports_vision(model)):
        # Fall back to text-only with a note so the agent knows.
        names = ", ".join(
            (img.source_path.name if img.source_path else "<inline>")
            for img in images
        )
        return (
            f"{text}\n\n[image attachments suppressed — model "
            f"{model!r} not detected as vision-capable: {names}]"
        )
    parts: list[dict] = [{"type": "text", "text": text}] if text else []
    for img in images:
        parts.append({
            "type": "image_url",
            "image_url": {"url": img.data_uri()},
        })
    return parts


__all__ = [
    "ImageAttachment",
    "ImageError",
    "model_supports_vision",
    "load_image",
    "to_openai_content",
]
