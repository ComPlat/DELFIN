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


# Budgets are set just above the CURRENT prompt sizes, so any regrowth
# fails immediately. The prompt diet removed duplicated rules (the same
# contract stated in a role prompt AND a shared addendum), worked examples
# that restated a rule already given, historical incident narrative, and
# blocks that re-listed what a tool schema already declares. No behavioral
# contract was dropped — extending one is fine, but pay for it by trimming
# elsewhere rather than by raising the number.
@pytest.mark.parametrize(
    "filename, max_tokens",
    [
        # 7950 -> 7150: dropped the two worked ACTION examples (they restate
        # the plan-first / verify-after-set rules stated above them), the
        # ORCA counter-example lists, and the duplicated tab-set + command-
        # discovery blocks.
        ("dashboard_agent.md", 7150),
        # 14200 -> 10600: dropped the worked-example dialogs and the
        # "how these compound" walk-through, folded the three separate
        # workspace-location statements into one, compressed the sandbox
        # internals, and removed tool-signature listings that duplicate the
        # tool schemas.
        ("solo_agent.md", 10600),
        # Written lean from the start: the shared addenda carry the general
        # contracts, so this prompt only states what is specific to working
        # on someone's real records. Raised as the mode's surface grew —
        # series work, record-addressed edits, PDF assembly, remembered
        # folder conventions, the working folder — each replacing something
        # the model would otherwise improvise, so each has to be named. The
        # prose was tightened three times on the way; what is left is
        # contract, not explanation.
        ("office_agent.md", 1600),
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


# The file budgets above guard the markdown. This one guards what the model
# actually receives: the CACHEABLE HEAD of the composed prompt (role prompt +
# shared addenda + project context). It is the part that is re-sent verbatim
# on every request of a session, so it is where prompt cost is decided.
@pytest.mark.parametrize(
    "role_id, mode_id, route, max_stable_tokens",
    [
        ("solo_agent", "solo", ["solo_agent"], 11400),
        ("dashboard_agent", "dashboard", ["dashboard_agent"], 10400),
    ],
)
def test_composed_stable_head_within_budget(
        monkeypatch, role_id, mode_id, route, max_stable_tokens):
    from delfin import user_settings
    from delfin.agent.prompt_loader import PromptLoader

    # Pin the lazy-module setting so the budget measures the prompt, not the
    # machine's local configuration.
    monkeypatch.setattr(
        user_settings, "load_settings",
        lambda *a, **k: {"agent": {"slim_prompt": True}})

    report = PromptLoader().prompt_size_report(
        role_id=role_id, mode_id=mode_id, route=route,
        task_text="fix the failing test in foo.py",
        session_key="budget-1")
    actual = report["stable_tokens"]
    assert actual <= max_stable_tokens, (
        f"{role_id}: cacheable head is {actual} tokens "
        f"(>{max_stable_tokens} budget). Trim before extending."
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
