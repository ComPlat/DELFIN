"""Tests for the live-regression fixes after the user reported
'Öffne Calculations' going silent in dashboard mode."""

from __future__ import annotations

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Reasoning-model detection (GPT-5 family now flagged)
# ---------------------------------------------------------------------------


def test_is_reasoning_detector_matches_gpt5_family():
    """Azure GPT-5.x are reasoning models that REQUIRE
    ``reasoning_effort`` to be set; otherwise they consume the budget
    on internal reasoning and emit zero text tokens (the live
    regression we just hit)."""
    from delfin.agent.api_client import OpenAIClient

    import re as _re

    def _is_reasoning(model_id: str) -> bool:
        base = (
            model_id.split(".", 1)[-1]
            if model_id.startswith(("azure.", "kit.")) else model_id
        )
        # ``o\d`` — o1/o3/o4-mini etc. ``opus`` must NOT match (it
        # starts with 'o' but isn't a reasoning model in the OpenAI
        # sense — and it goes through APIClient anyway).
        return bool(
            _re.match(r"^o\d", base) or base.startswith("gpt-5")
        )

    # GPT-5 family
    assert _is_reasoning("azure.gpt-5") is True
    assert _is_reasoning("azure.gpt-5.4") is True
    assert _is_reasoning("azure.gpt-5-mini") is True
    assert _is_reasoning("azure.gpt-5-nano") is True
    # o-series stays detected
    assert _is_reasoning("azure.o3") is True
    assert _is_reasoning("azure.o4-mini") is True
    assert _is_reasoning("o4-mini") is True
    # Non-reasoning models stay non-reasoning
    assert _is_reasoning("azure.gpt-4.1") is False
    assert _is_reasoning("azure.gpt-4.1-mini") is False
    assert _is_reasoning("opus") is False
    assert _is_reasoning("sonnet") is False
    assert _is_reasoning("haiku") is False


def test_is_reasoning_baked_into_source():
    """Regression guard against an accidental revert."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "agent" / "api_client.py")
    text = p.read_text(encoding="utf-8")
    # Tightened detector: r"^o\d" prevents Anthropic "opus" from
    # being mis-classified as a reasoning model (would still go
    # through APIClient, but a future refactor could route it here
    # and the false positive would manifest).
    assert r'_re_reason.match(r"^o\d"' in text
    assert '_base.startswith("gpt-5")' in text
    # Comment explaining WHY (so a future cleanup doesn't drop one)
    assert "GPT-5 family" in text
    assert "silently consume the budget" in text


# ---------------------------------------------------------------------------
# Dashboard mode auto-clamps high/xhigh effort to low
# ---------------------------------------------------------------------------


def test_dashboard_thinking_budget_keeps_reasoning_at_low():
    """The base thinking-budget for dashboard_agent must keep the
    derived ``reasoning_effort`` at 'low' across every effort dropdown
    setting (low / medium / high / xhigh). Otherwise Azure GPT-5.x
    silently burns 60-120 s on hidden reasoning for tab-open-style
    dispatch requests.

    Mapping (from OpenAIClient.stream_message):
      <  16 k → reasoning_effort='low'
      < 64 k → 'medium'
      ≥ 64 k → 'high'

    With base=4 k and the four effort multipliers (0.5/1.0/2.0/3.0)
    the largest possible budget is 12 k — still < 16 k → low.
    """
    from delfin.agent.engine import _ROLE_THINKING_BUDGETS
    base = _ROLE_THINKING_BUDGETS["dashboard_agent"]
    assert base * 3.0 < 16000, (
        f"dashboard_agent base budget {base} × xhigh-mult (3.0) = "
        f"{base * 3.0} would push reasoning_effort above 'low' — "
        f"Azure GPT-5 would hang again."
    )


def test_dashboard_mode_clamps_high_effort_to_low():
    """Dashboard tab-open should never need deep reasoning. High/xhigh
    against a reasoning model takes 2+ minutes for what should be a
    100 ms routing decision — the auto-clamp prevents the user from
    accidentally paying that cost."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    # The clamp block sits next to the effort-multiplier table
    idx = text.find("Dashboard mode auto-clamps effort")
    assert idx > 0, "missing dashboard effort auto-clamp block"
    snippet = text[idx: idx + 1200]
    assert 'mode_dropdown.value == "dashboard"' in snippet
    assert '_effective_effort in ("high", "xhigh")' in snippet
    assert '_effective_effort = "low"' in snippet


# ---------------------------------------------------------------------------
# /tab fuzzy-match (typos like "claultions" land on Calculations)
# ---------------------------------------------------------------------------


def test_tab_fuzzy_match_threshold_catches_real_typos():
    """The threshold (0.6 SequenceMatcher ratio) must catch the typos
    users actually make — 1-2 char transpositions/drops on tab names."""
    import difflib

    # Mirror the alias keys the production /tab handler exposes
    candidates = [
        "submit", "recalc", "jobs", "orca", "calc",
        "calcs", "calculation", "calculations",
        "berechnung", "berechnungen", "archive", "archiv",
        "literature", "literatur", "agent", "settings",
        "einstellungen", "job",
    ]

    def _best(typo: str) -> tuple[float, str]:
        scored = sorted(
            (
                (difflib.SequenceMatcher(None, typo.lower(), c.lower()).ratio(), c)
                for c in candidates
            ),
            reverse=True,
        )
        return scored[0]

    # The live regression: "Öffne Calculations" → agent emits /tab claultions
    ratio, target = _best("claultions")
    assert target in ("calculations", "calcs", "calculation")
    assert ratio >= 0.6
    # Other realistic typos
    for typo, expected in [
        ("Calcultaions", "calculations"),
        ("clacuations",  "calculations"),
        ("submmit",      "submit"),
        ("joobs",        "jobs"),
        ("setings",      "settings"),
    ]:
        ratio, target = _best(typo)
        assert ratio >= 0.6, f"{typo} → {ratio:.2f} (below threshold)"
        assert target == expected or expected in target, (
            f"{typo} → {target} (expected {expected})"
        )


def test_tab_fuzzy_match_block_present_in_source():
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find('# Fuzzy match: typos like "claultions"')
    assert idx > 0, "fuzzy-match block missing in /tab handler"
    snippet = text[idx: idx + 1200]
    assert "SequenceMatcher" in snippet
    assert ">= 0.6" in snippet
    assert "Fuzzy-matched" in snippet
