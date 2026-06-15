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


def test_delegate_defaults_to_general_purpose():
    out = builtin_subagent_command("delegate", "run the benchmark and summarise")
    assert "subagent_type='general-purpose'" in out


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
