"""Built-in subagent slash commands (/explore /review /plan /delegate).

A direct user lever to delegate to a sub-agent: each expands to a directive
that steers the agent to call the subagent tool with the right preset.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.slash_commands import builtin_subagent_command


def test_explore_maps_to_explore_preset():
    out = builtin_subagent_command("explore", "find all call-sites of foo()")
    assert "subagent_type='explore'" in out
    assert "find all call-sites of foo()" in out


def test_review_maps_to_code_reviewer():
    out = builtin_subagent_command("review", "audit the diff in module x")
    assert "subagent_type='code-reviewer'" in out


def test_plan_maps_to_plan_preset():
    out = builtin_subagent_command("plan", "design the migration")
    assert "subagent_type='plan'" in out


def test_delegate_names_a_specialist_and_labels_it_a_suggestion():
    """CHANGED ON PURPOSE (wave 43 s1) — because the TEST asserted the
    defect, not a contract.

    It read ``test_delegate_defaults_to_general_purpose`` and asserted
    ``general-purpose`` for the text "run the benchmark and summarise".
    The text names work one specialist's own description claims, and the
    assertion held for EVERY possible text: /delegate was a constant, which
    is exactly what the labelled set measured
    (tests/fixtures/specialist_routing_w43s1.yaml: 28 of 30 requests
    reached the same row, accuracy 0.067 on 41fbf18e2). Pinning a constant
    is pinning the thing this package exists to remove.

    What is pinned now is the contract, not one scorer's opinion about a
    debatable row: the directive names a registered preset, it says the
    choice was a SUGGESTION rather than the user's own naming, and the
    user's text survives into the briefing. Which specialist a debatable
    request deserves is measured by the matrix, not asserted here.
    """
    out = builtin_subagent_command("delegate", "run the benchmark and summarise")
    assert "subagent_type='" in out
    assert "suggest" in out.lower(), (
        "a specialist the user did not name must arrive labelled a "
        "suggestion, or a wrong pick is indistinguishable from a wish")
    assert "run the benchmark and summarise" in out


def test_a_request_with_one_clear_owner_is_not_left_with_the_catch_all():
    """The falsifiable half: a request whose cue names exactly one
    specialist must not be answered with the preset whose own description
    says it is for when no narrower one fits.

    FALSIFIED on the base (which answers general-purpose to everything)
    and falsified by any scorer that ranks the catch-all first.
    """
    out = builtin_subagent_command(
        "delegate",
        "Find all callers of parse_control and tell me where it is defined.")
    assert "subagent_type='explore'" in out, out[:200]
    assert "suggest" in out.lower()


def test_delegate_falls_back_to_general_purpose_and_says_why():
    """The default still exists — it is no longer silent."""
    out = builtin_subagent_command("delegate", "Mach das bitte schnell.")
    assert "subagent_type='general-purpose'" in out
    assert "no specialist" in out.lower() or "unsure" in out.lower(), out[:250]


def test_delegate_with_explicit_type():
    out = builtin_subagent_command("delegate", "explore where config is parsed")
    assert "subagent_type='explore'" in out
    assert "where config is parsed" in out


def test_leading_slash_is_tolerated():
    assert builtin_subagent_command("/explore", "x y z") is not None


def test_unknown_command_returns_none():
    assert builtin_subagent_command("frobnicate", "x") is None


def test_empty_task_asks_for_one():
    out = builtin_subagent_command("explore", "")
    assert out is not None
    assert "explore" in out.lower()


def test_mentions_background_and_result_collection():
    out = builtin_subagent_command("delegate", "long independent research task")
    assert "background=true" in out
    assert "subagent_result" in out


def test_dashboard_wires_builtin_subagent_commands():
    """The dashboard dispatch must call builtin_subagent_command and the
    autocomplete must list the four commands."""
    src = (Path(__file__).resolve().parent.parent
           / "delfin" / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    assert "builtin_subagent_command(_name, _rest)" in src
    for cmd in ("/explore", "/review", "/plan", "/delegate"):
        assert f'"{cmd}"' in src
