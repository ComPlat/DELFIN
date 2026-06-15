"""Sanitise corrupted model output (harmony tool-channel leaks + glitch tokens).

Observed in production with azure.gpt-5.x served through the KIT
OpenAI-compatible endpoint: the model signals tool calls in its native
"harmony" channel syntax (``to=<tool> <json>``) which the endpoint passes
through as **text** instead of structured ``tool_calls``, and the channel's
special tokens decode into low-frequency multilingual "glitch" tokens
(e.g. ``手机天天中彩票``, ``출장샵``, ``ացին``).  The user sees Chinese / Korean /
Armenian garbage interleaved with a leaked ``to=search_docs {…}`` fragment.

This module repairs that text: it strips the leaked tool-channel fragments
(reporting which tools the model *intended* to call) and removes runs of
scripts that never legitimately appear in DELFIN's German/English chemistry
output.  It is a pure, dependency-free function — safe to run on every
response (a no-op on clean text).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field


# Leaked harmony tool-channel fragment: ``to=<tool> <junk> {json}``.
# The junk between the tool name and the JSON is the corrupted special
# tokens (``json_schema``, glitch chars).  DOTALL so the JSON can wrap.
_LEAKED_TOOL = re.compile(
    r"to=\s*([A-Za-z_][\w\-]*)[^\{\n]*?(\{.*?\})",
    re.DOTALL,
)

# Bare leftover ``to=<tool>`` with no JSON (defensive — strip the marker).
_BARE_TOOL_MARKER = re.compile(r"to=\s*[A-Za-z_][\w\-]*")

# Harmony special-token leftovers that decode as literal words.
_HARMONY_TOKENS = re.compile(r"\b(?:json_schema|json|constrain)\b(?=[^\sA-Za-z])")

# Reasoning-tag models (deepseek-r1, qwq, qwen3-thinking, …) emit their chain
# of thought as <think>…</think> in the visible text channel. Strip the whole
# block; also strip a dangling unterminated <think> (streaming can cut off
# before </think>) so no reasoning leaks into the user-visible answer.
_THINK_BLOCK = re.compile(r"<think>.*?</think>", re.DOTALL | re.IGNORECASE)
_THINK_DANGLING = re.compile(r"<think>.*$", re.DOTALL | re.IGNORECASE)

# Scripts that do not occur in DELFIN's de/en chemistry output — a run of
# these is glitch-token corruption, not content.  (CJK, kana, Hangul,
# Armenian, Cyrillic, Malayalam, Georgian, Thai, Devanagari, fullwidth.)
_GLITCH = re.compile(
    "["
    "　-〿぀-ヿ㐀-䶿一-鿿豈-﫿"
    "가-힯ᄀ-ᇿ"
    "Ѐ-ӿ"
    "԰-֏"
    "ഀ-ൿ"
    "Ⴀ-ჿ"
    "฀-๿"
    "ऀ-ॿ"
    "＀-￯"
    "]+"
)


@dataclass
class SanitizeResult:
    text: str
    leaked_tools: list[str] = field(default_factory=list)
    glitch_chars: int = 0
    think_stripped: bool = False

    @property
    def changed(self) -> bool:
        return (bool(self.leaked_tools) or self.glitch_chars > 0
                or self.think_stripped)


def sanitize_agent_text(text: str) -> SanitizeResult:
    """Return a cleaned copy of ``text`` plus what was repaired.

    - Strips leaked ``to=<tool> {json}`` fragments (the model's intended
      tool calls that leaked into the text channel) and reports the tool
      names in ``leaked_tools``.
    - Removes runs of non-DELFIN scripts (glitch tokens), counting how many
      characters were dropped in ``glitch_chars``.

    Clean text passes through unchanged (``changed`` is False).
    """
    if not text:
        return SanitizeResult(text=text or "")

    leaked: list[str] = []
    for m in _LEAKED_TOOL.finditer(text):
        name = m.group(1)
        if name not in leaked:
            leaked.append(name)

    think_stripped = bool(_THINK_BLOCK.search(text) or _THINK_DANGLING.search(text))
    cleaned = _THINK_BLOCK.sub(" ", text)
    cleaned = _THINK_DANGLING.sub(" ", cleaned)
    cleaned = _LEAKED_TOOL.sub(" ", cleaned)
    cleaned = _BARE_TOOL_MARKER.sub(" ", cleaned)
    cleaned = _HARMONY_TOKENS.sub(" ", cleaned)

    glitch_chars = sum(len(m.group(0)) for m in _GLITCH.finditer(cleaned))
    cleaned = _GLITCH.sub("", cleaned)

    # Collapse the whitespace/newlines left behind by the removals.
    cleaned = re.sub(r"[ \t]{2,}", " ", cleaned)
    cleaned = re.sub(r"\n{3,}", "\n\n", cleaned)
    cleaned = cleaned.strip()

    return SanitizeResult(
        text=cleaned,
        leaked_tools=leaked,
        glitch_chars=glitch_chars,
        think_stripped=think_stripped,
    )


def parse_leaked_tool_calls(text: str) -> list[dict]:
    """Best-effort recovery of the JSON args from leaked tool fragments.

    Returns ``[{"name": str, "arguments": dict}]`` for each parseable
    fragment.  Used only for diagnostics / future re-dispatch — NOT executed
    automatically (parsing tool calls out of corrupted text is unsafe).
    """
    out: list[dict] = []
    for m in _LEAKED_TOOL.finditer(text):
        name = m.group(1)
        try:
            args = json.loads(m.group(2))
        except Exception:
            args = {}
        out.append({"name": name, "arguments": args})
    return out
