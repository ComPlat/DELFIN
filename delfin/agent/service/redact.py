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
        self._lookup_failed = False
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
                # A successful refresh ends the degraded window -- the
                # store recovered, masking resumes normally.
                self._lookup_failed = False
            except Exception as exc:  # noqa: BLE001 - the redactor keeps
                # serving, but the degraded mask list is VISIBLE, never
                # silent: with [] a shapeless held credential (a KIT key
                # has none) would pass the shape pass untouched, so
                # masking is WITHHELD for this refresh window instead.
                import logging
                logging.getLogger(__name__).warning(
                    "held-credential lookup failed (%s); masking is "
                    "withheld until the next refresh", type(exc).__name__)
                self._known = []
                self._lookup_failed = True
        return self._known

    def text(self, text: str) -> str:
        # The lookup must run BEFORE the degraded-window check: if the
        # flag were read first, the one call that hits the failure would
        # still return the text unmasked (operator node run, 20.11).
        self._known_values()
        if self._lookup_failed:
            return "[redaction failed -- text withheld]"
        for secret in self._known:
            if secret in text:
                text = text.replace(secret, "[redacted:credential]")
        if self._shapes is not None:
            try:
                text = self._shapes(text, [])
            except Exception as exc:  # noqa: BLE001 - fail closed, the
                # same failure mode as output_guard's callers.
                import logging
                logging.getLogger(__name__).warning(
                    "secret-shape pass failed (%s)", type(exc).__name__)
                return "[redaction failed -- text withheld]"
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
