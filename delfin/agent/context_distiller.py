"""Cheap-model context distillation for expensive inference calls.

When the estimated input token count exceeds a threshold, uses a fast/cheap
model (Haiku) to compress the system prompt before sending it to the
expensive model (Opus/Sonnet).  Layer 0 (role identity + critical rules)
is never distilled — it is prepended verbatim.

This is opt-in and only triggers when the context is large enough that
the distillation cost (~$0.04) is outweighed by the savings on the
expensive model (~$0.45 for 30k fewer Opus tokens).
"""

from __future__ import annotations

import re
from typing import Any

# Rough chars-per-token estimate.  Anthropic models average ~3.5-4 chars/token.
_CHARS_PER_TOKEN = 4

# Only distill when estimated input tokens exceed this threshold.
_TOKEN_THRESHOLD = 40_000

# Target output size for the distilled context.
_TARGET_CHARS = 12_000  # ~3000 tokens

# Layer 0 marker: everything before this marker is kept verbatim.
_LAYER0_END_MARKER = "--- Project Context ---"
_LAYER0_END_MARKERS = [
    "--- Project Context ---",
    "--- Repo Map ---",
    "--- Provider Profile ---",
    "--- Task Briefing ---",
]

_DISTILL_PROMPT = """\
Compress the following system context for a coding agent. The agent will \
use this context to complete a software engineering task.

KEEP (verbatim where possible):
- Role identity and critical rules
- Specific file paths, function names, line numbers
- Acceptance criteria and success metrics
- Error patterns and denied command lists
- Concrete numbers (percentages, thresholds, costs)

REMOVE:
- Verbose explanations and boilerplate
- Redundant examples
- Generic advice the model already knows
- Repeated information

Output a compressed version under {target_tokens} tokens. \
Preserve structure (headers, bullet points).

CONTEXT TO COMPRESS:
{context}"""


def estimate_tokens(text: str) -> int:
    """Rough token count estimate from character count."""
    return len(text) // _CHARS_PER_TOKEN


def _split_layer0(system_prompt: str) -> tuple[str, str]:
    """Split system prompt into layer 0 (kept verbatim) and the rest.

    Layer 0 ends at the first occurrence of a known section marker.
    """
    best_pos = len(system_prompt)
    for marker in _LAYER0_END_MARKERS:
        pos = system_prompt.find(marker)
        if pos != -1 and pos < best_pos:
            best_pos = pos
    if best_pos == len(system_prompt):
        # No marker found — treat first 30% as layer 0
        best_pos = len(system_prompt) * 3 // 10
    return system_prompt[:best_pos], system_prompt[best_pos:]


class ContextDistiller:
    """Compress large context using a cheap model before expensive inference."""

    def __init__(
        self,
        *,
        token_threshold: int = _TOKEN_THRESHOLD,
        enabled: bool = False,
    ):
        self.token_threshold = token_threshold
        self.enabled = enabled
        self._client: Any = None

    def should_distill(
        self,
        system_prompt: str,
        messages: list[dict[str, str]],
    ) -> bool:
        """Check whether distillation would be beneficial."""
        if not self.enabled:
            return False
        total_chars = len(system_prompt) + sum(
            len(m.get("content", "")) for m in messages
        )
        return (total_chars // _CHARS_PER_TOKEN) > self.token_threshold

    def distill(
        self,
        system_prompt: str,
        task_text: str = "",
    ) -> str:
        """Compress the system prompt using a cheap model.

        Returns the distilled system prompt with layer 0 prepended verbatim.
        Falls back to extractive compression if the API call fails.
        """
        layer0, rest = _split_layer0(system_prompt)
        if not rest.strip():
            return system_prompt

        # Try API-based distillation first
        try:
            distilled = self._api_distill(rest)
            if distilled and len(distilled) < len(rest):
                return layer0 + distilled
        except Exception:
            pass

        # Fallback: extractive compression
        return layer0 + self._extractive_compress(rest)

    def _api_distill(self, context: str) -> str:
        """Call a cheap model to compress the context."""
        if self._client is None:
            try:
                import anthropic
                self._client = anthropic.Anthropic()
            except Exception:
                return ""

        target_tokens = _TARGET_CHARS // _CHARS_PER_TOKEN
        prompt = _DISTILL_PROMPT.format(
            target_tokens=target_tokens,
            context=context[:80_000],  # cap input to avoid excessive cost
        )
        try:
            response = self._client.messages.create(
                model="claude-haiku-4-5-20251001",
                max_tokens=target_tokens + 500,
                messages=[{"role": "user", "content": prompt}],
            )
            return response.content[0].text if response.content else ""
        except Exception:
            return ""

    @staticmethod
    def _extractive_compress(text: str, target: int = _TARGET_CHARS) -> str:
        """Extract structured lines (headers, bullets, paths) from text.

        No API call — pure heuristic extraction.  Used as fallback.
        """
        lines = text.splitlines()
        selected: list[str] = []
        total = 0
        for line in lines:
            stripped = line.strip()
            if not stripped:
                continue
            is_structural = (
                stripped.startswith("#")
                or stripped.startswith("- ")
                or stripped.startswith("* ")
                or re.match(r"^\d+\.\s", stripped)
                or stripped.startswith("**")
                or ":" in stripped[:40]
                or "/" in stripped  # likely a file path
                or stripped.startswith("<critical>")
                or stripped.startswith("</critical>")
            )
            if is_structural:
                selected.append(line)
                total += len(line)
                if total >= target:
                    break
        return "\n".join(selected)
