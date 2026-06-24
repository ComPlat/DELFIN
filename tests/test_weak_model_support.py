"""Tests for the weak-model audit pass: compact prompt + /plan
approve/reject fallback."""

from __future__ import annotations

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Weak-model auto-detection
# ---------------------------------------------------------------------------


@pytest.fixture
def loader():
    from delfin.agent.prompt_loader import PromptLoader
    return PromptLoader()


@pytest.mark.parametrize("model, expected_weak", [
    # Strong cloud / API models — never compact
    ("opus",                            False),
    ("sonnet",                          False),
    ("haiku",                           False),
    ("gpt-5.4",                         False),
    ("azure.gpt-5.4",                   False),
    ("azure.gpt-5-mini",                False),
    # Strong local / MoE models — handle the slim prompt fine
    ("kit.gpt-oss-120b",                False),
    ("kit.qwen3.5-397b-A17b",           False),
    ("qwen3-coder:32b",                 False),
    ("llama3.3:70b",                    False),
    ("gemma:27b",                       False),
    ("mistral-small:22b",               False),   # >14b via generic size
    ("deepseek-r1:671b",                False),   # huge MoE
    # Weak local models — compact prompt activates
    ("qwen2.5-coder:7b",                True),
    ("qwen2.5-coder:14b",               True),
    ("llama3.2:8b",                     True),
    ("llama3.2:3b",                     True),
    ("gemma4:7b",                       True),
    ("gemma2:2b",                       True),
    ("phi-3.5",                         True),
    ("mistral-7b",                      True),
    ("codellama:7b",                    True),
    # Generic size signal catches tags the family regex missed (the bug):
    ("qwen3-vl:4b",                     True),    # was wrongly classified strong
    ("gemma3:2b",                       True),
    ("qwen3:0.6b",                      True),
    ("deepseek-r1:7b",                  True),
])
def test_is_weak_model_classification(loader, model, expected_weak):
    assert loader._is_weak_model(model) is expected_weak, (
        f"{model}: expected weak={expected_weak}, "
        f"got {loader._is_weak_model(model)}"
    )


def test_is_weak_model_empty_or_unknown(loader):
    assert loader._is_weak_model("") is False
    assert loader._is_weak_model("totally-unknown-model") is False


# ---------------------------------------------------------------------------
# Compact prompt actually shrinks for weak models
# ---------------------------------------------------------------------------


def test_compact_prompt_shrinks_solo_for_weak_model(loader):
    """For a weak model the slim prompt should be materially smaller
    than for a strong one — paragraphs collapsed to first-sentence."""
    text = loader.load_role_prompt("solo_agent")
    strong = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo", model="opus",
    )
    weak = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo", model="qwen2.5-coder:7b",
    )
    assert len(weak) < len(strong), (
        f"compact mode didn't shrink: weak={len(weak)} strong={len(strong)}"
    )
    # Target: at least 15% smaller than the already-slim prompt
    savings_ratio = (len(strong) - len(weak)) / len(strong)
    assert savings_ratio > 0.15, (
        f"compact savings too small: {savings_ratio:.1%} "
        f"(expected >15%)"
    )


def test_compact_prompt_preserves_section_headers(loader):
    """Compaction must keep structural anchors — section headers,
    list items, code blocks — so the model can navigate the prompt."""
    text = loader.load_role_prompt("solo_agent")
    compact = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo", model="qwen2.5-coder:7b",
    )
    # Core sections must survive. Memory guidance now lives once in the
    # universal memory_addendum.md (injected for every role), not in the role
    # prompt — so it's no longer expected here.
    for needle in (
        "## Strategies for approaching tasks",
        "## Subagents",
        "## Context management",
    ):
        assert needle in compact, f"compact lost section: {needle}"


def test_memory_guidance_present_once_in_full_prompt(loader):
    """The memory section was duplicated (role prompt + universal addendum);
    it's now consolidated into memory_addendum.md so the assembled prompt
    carries it exactly once, with the key guidance intact."""
    sp = loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo", mode_description="Code",
        route=["solo_agent"], role_index=0, prior_outputs=None,
        memory_context="", task_text="x", session_key="m", live_state=None,
        model="kit.qwen3.5-397b-A17b", permission_mode="ask_all")
    assert sp.count("## Memory") == 1                 # no duplication
    assert "remember" in sp and "durable" in sp        # proactive-save guidance
    assert "still holds" in sp                          # verify-before-acting bit
    assert "~/.claude/" not in sp                       # the wrong path is gone


def test_compact_prose_keeps_code_blocks_verbatim(loader):
    """Code blocks are examples — high signal for weak models."""
    sample = (
        "# Title\n\n"
        "Some long prose. This sentence should be dropped.\n"
        "More wordy explanation that should be collapsed.\n\n"
        "```python\n"
        "x = 1\n"
        "y = 2\n"
        "```\n\n"
        "Another paragraph. Trailing text.\n"
    )
    out = loader._compact_prose(sample)
    assert "```python" in out and "x = 1" in out and "y = 2" in out
    # First sentence of each prose paragraph survives
    assert "Some long prose." in out
    assert "Another paragraph." in out


def test_compact_prose_keeps_bullets_and_tables(loader):
    sample = (
        "## Section\n\n"
        "Lead-in prose. Detail to drop.\n\n"
        "- Bullet one\n"
        "- Bullet two\n\n"
        "| col1 | col2 |\n"
        "|---|---|\n"
        "| a | b |\n"
    )
    out = loader._compact_prose(sample)
    assert "- Bullet one" in out
    assert "- Bullet two" in out
    assert "| col1 | col2 |" in out


# ---------------------------------------------------------------------------
# /plan approve | /plan reject fallback
# ---------------------------------------------------------------------------


def test_plan_approve_fallback_present_in_dashboard():
    """The handler must exist + be wired into _SLASH_COMMANDS so the
    palette + autocomplete pick it up."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    # Handler
    idx = text.find('if cmd in ("/plan approve", "/plan reject"):')
    assert idx > 0, "fallback handler missing"
    snippet = text[idx: idx + 2500]
    assert "_plan_approval_event" in snippet
    # Approve/reject is delegated to the shared _finalize_plan_decision
    # helper (same code path as the 'Plan akzeptieren' button) so the two
    # can't drift — the handler picks the branch via approved=(cmd==...).
    assert "_finalize_plan_decision(" in snippet
    assert 'approved=(cmd == "/plan approve")' in snippet
    assert "Plan approved via /plan approve fallback" in snippet
    assert "Plan rejected via /plan reject fallback" in snippet
    # Slash registry entry so the palette shows it
    assert '"/plan approve"' in text
    assert '"/plan reject"' in text


def test_plan_approve_without_pending_plan_is_friendly():
    """Calling /plan approve when nothing is awaiting must NOT crash
    — it should print a helpful message pointing at /mode plan."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("No plan is currently awaiting approval")
    assert idx > 0
    # And the surrounding code does an early-return so the rest of
    # the handler doesn't fire with None values
    block = text[max(0, idx-400): idx + 400]
    assert "if ev is None or result is None:" in block


# ---------------------------------------------------------------------------
# Threading: build_system_prompt accepts model + plumbs to stripper
# ---------------------------------------------------------------------------


def test_build_system_prompt_threads_model_param(loader):
    """The new ``model`` kwarg on build_system_prompt must reach the
    compact-prompt decision — without it, weak-model auto-detection
    is dead code."""
    text_w = loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo",
        task_text="hi", model="qwen2.5-coder:7b",
    )
    text_s = loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo",
        task_text="hi", model="opus",
    )
    # Weak model gets a meaningfully smaller solo prompt
    assert len(text_w) < len(text_s)
