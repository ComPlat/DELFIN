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


# Budgets reflect the post-P4 prompts (verify-enforcement, the 8-pattern
# playbook, ORCA-builder grounding docs were all deliberate additions
# after the original S2 slim-down) with ~8% headroom, so any further
# creep fails immediately. A deliberate "prompt diet" — trimming these
# back down WITH benchmark validation — is on the backlog; shrink the
# budgets again when it lands.
@pytest.mark.parametrize(
    "filename, max_tokens",
    [
        ("dashboard_agent.md", 7500),
        # Raised 13800 -> 14000 for the workspace/CLI orchestration guidance
        # (launch-dir = workspace). Content is intentional; keep growth in check.
        ("solo_agent.md", 14000),
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
    """The prompt must keep the irreducible safety + grounding sections.

    The list pins TODAY's intentional design (dashboard mode is
    guide+UI only since 3a2c802; ORCA claims must be grounded in the
    manual since P4) — update it consciously when the design changes,
    never to silence a failure.
    """
    text = (_PROMPT_DIR / "dashboard_agent.md").read_text()
    must_have = [
        "ACTION:",                       # how commands work
        "Safety rules",                  # safety policy
        "Hard scope limits",             # guide+UI-only contract (3a2c802)
        "Ground every ORCA",             # manual-grounding rule (P4)
        "Tools you may NOT use",         # forbidden tool surface
        "agent_workspace",               # spelled out as NOT available
    ]
    missing = [k for k in must_have if k not in text]
    assert not missing, f"essential sections dropped: {missing}"
