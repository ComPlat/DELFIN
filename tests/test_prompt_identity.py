"""The agent prompts must not give the model a fixed 'Claude' identity.

DELFIN runs on whatever backend the user picked (qwen/KIT, GPT, Ollama, Claude).
A subagent is a fresh instance of the SAME model — not "a fresh Claude". Bug
20260625 (ka_ew7404, kit.qwen3.5-397b): the prompt said subagents "spawn a fresh
Claude", so the qwen agent told the user its subagents were "Claude instances"
("Warum redet er von claude es ist doch qwen"). These guards keep the prompts
model-neutral. Legitimate references (the "Claude CLI" backend option, the
CLAUDE.md filename, ~/.claude path conventions) are NOT identity claims and are
deliberately not matched.
"""

from __future__ import annotations

from pathlib import Path

import pytest

_PACK = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"

# Phrasings that ascribe a Claude identity to the agent or its subagents.
_IDENTITY_LEAKS = (
    "fresh Claude",
    "a Claude ",
    "Claude instance",
    "Claude with its own",
    "isolated Claude",
    "Claude-Instanz",
)


@pytest.mark.parametrize("md", sorted(_PACK.rglob("*.md")), ids=lambda p: p.name)
def test_prompt_does_not_assign_a_claude_identity(md):
    text = md.read_text(encoding="utf-8")
    for phrase in _IDENTITY_LEAKS:
        assert phrase not in text, (
            f"{md.relative_to(_PACK)} describes the agent/subagent as Claude "
            f"({phrase!r}) — leaks a wrong identity on non-Claude backends."
        )
