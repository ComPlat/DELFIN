"""Credentials masked in what the agent service shows and keeps: a
session's events (the worker), and what needs the user (the attention
feed, which reads records other DELFIN processes wrote).

Light on purpose: the service process imports it, never the worker's
engine side.
"""

from __future__ import annotations

import time
from typing import Any, Callable


class Redactor:
    """Masks credentials in what a session shows and keeps.

    Two passes, as ``output_guard.scrub_secrets`` does: the exact values of
    the credentials DELFIN holds (a KIT key has no recognisable shape), then
    everything shaped like a credential. The known values are read once a
    minute, not once per string.
    """

    _REFRESH_S = 60.0

    def __init__(self):
        self._known: list[str] = []
        self._loaded = -1e9
        try:
            from ..output_guard import _redact_secrets
            self._shapes: Callable[[str, list], str] | None = _redact_secrets
        except Exception:
            self._shapes = None

    def _known_values(self) -> list[str]:
        now = time.monotonic()
        if now - self._loaded > self._REFRESH_S:
            self._loaded = now
            try:
                from ..web_tools import _known_secret_values
                self._known = sorted(set(_known_secret_values()), key=len,
                                     reverse=True)
            except Exception:
                self._known = []
        return self._known

    def text(self, text: str) -> str:
        for secret in self._known_values():
            if secret in text:
                text = text.replace(secret, "[redacted:credential]")
        if self._shapes is not None:
            try:
                text = self._shapes(text, [])
            except Exception:
                pass
        return text

    def value(self, value: Any) -> Any:
        if isinstance(value, str):
            return self.text(value)
        if isinstance(value, dict):
            return {k: self.value(v) for k, v in value.items()}
        if isinstance(value, list):
            return [self.value(v) for v in value]
        return value

    def event(self, event: dict, identifiers: frozenset = frozenset()) -> dict:
        return {k: (v if k in identifiers else self.value(v))
                for k, v in event.items()}
