"""Token budget guards for dashboard / solo agent prompts.

These tests fail if a future edit blows the role-prompt size past a
fixed budget. Token estimation uses ``len(text) / 4`` (a tight
upper-bound vs. real BPE tokenization).
"""
from __future__ import annotations

from pathlib import Path

import pytest


_REPO = Path(__file__).resolve().parents[1]
_PROMPT_DIR = _REPO / "delfin" / "agent" / "pack" / "agents"


def _estimate_tokens(text: str) -> int:
    """Fast ~upper-bound token count (1 token ≈ 4 chars in English/German)."""
    return (len(text) + 3) // 4


@pytest.mark.parametrize(
    "filename, max_tokens",
    [
        ("dashboard_agent.md", 3000),
        ("solo_agent.md", 4000),
    ],
)
def test_role_prompt_within_token_budget(filename, max_tokens):
    """Role prompts must stay below their per-role token budget."""
    path = _PROMPT_DIR / filename
    assert path.exists(), f"missing prompt file: {path}"
    text = path.read_text()
    actual = _estimate_tokens(text)
    assert actual <= max_tokens, (
        f"{filename}: {actual} tokens (>{max_tokens} budget). "
        f"Trim before extending."
    )


def test_dashboard_prompt_keeps_essential_sections():
    """Slimmed prompt must still cover the irreducible safety + tool basics."""
    text = (_PROMPT_DIR / "dashboard_agent.md").read_text()
    must_have = [
        "ACTION:",                              # how commands work
        "Safety rules",                         # safety policy
        "agent_workspace",                      # where Write/Bash are allowed
        "Mutating MCP-ops calls",               # explicit allow_mutate gate
        "Background tasks — anti-stall rule",   # learnt-the-hard-way rule
        "Live state",                           # tells agent to read injected state
    ]
    missing = [k for k in must_have if k not in text]
    assert not missing, f"essential sections dropped: {missing}"
