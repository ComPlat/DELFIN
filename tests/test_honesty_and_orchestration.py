"""Honesty addendum (universal grounding) + subagent auto-isolation."""

from __future__ import annotations

import pytest


# --- Honesty addendum -------------------------------------------------------


def test_honesty_addendum_file_exists():
    from pathlib import Path
    p = (Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"
         / "shared" / "honesty_addendum.md")
    body = p.read_text(encoding="utf-8")
    assert "Verify before you claim" in body
    assert "Never invent" in body


def test_honesty_addendum_injected_into_prompt():
    from delfin.agent.prompt_loader import PromptLoader
    p = PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo", task_text="fix a bug")
    assert "Honesty & grounding" in p or "Verify before you claim" in p


def test_honesty_addendum_present_for_dashboard_role():
    from delfin.agent.prompt_loader import PromptLoader
    p = PromptLoader().build_system_prompt(
        role_id="dashboard_agent", mode_id="dashboard", task_text="set BP86")
    assert "Verify before you claim" in p


# --- Subagent writer detection + auto-isolation -----------------------------


def test_is_writer_preset():
    from delfin.agent import subagents as sa
    assert sa.is_writer_preset("general-purpose") is True
    assert sa.is_writer_preset("explore") is False
    assert sa.is_writer_preset("plan") is False
    assert sa.is_writer_preset("code-reviewer") is False
    assert sa.is_writer_preset("does-not-exist") is False


def _tc(i, sa_type, isolation=None):
    args = {"subagent_type": sa_type, "description": "d", "prompt": "p" * 25}
    if isolation is not None:
        args["isolation"] = isolation
    import json
    return {"id": str(i), "function": {"name": "subagent",
            "arguments": json.dumps(args)}}


def test_fanout_auto_isolates_parallel_writers(monkeypatch):
    from delfin.agent import api_client as A
    captured = []
    monkeypatch.setattr(A._doc_executor, "execute",
                        lambda name, args, permissions=None: captured.append(dict(args)) or "{}")
    futures, ex = A._fan_out_subagents(
        [_tc(1, "general-purpose"), _tc(2, "explore")], permissions=None)
    for f in futures.values():
        f.result()
    ex.shutdown(wait=True)
    by_type = {c["subagent_type"]: c for c in captured}
    assert by_type["general-purpose"].get("isolation") == "worktree"
    assert by_type["explore"].get("isolation") in (None, "")


def test_fanout_single_subagent_not_isolated(monkeypatch):
    from delfin.agent import api_client as A
    # <2 subagents → no fan-out, no auto-isolation
    futures, ex = A._fan_out_subagents([_tc(1, "general-purpose")], permissions=None)
    assert futures == {} and ex is None


def test_fanout_respects_explicit_isolation(monkeypatch):
    from delfin.agent import api_client as A
    captured = []
    monkeypatch.setattr(A._doc_executor, "execute",
                        lambda name, args, permissions=None: captured.append(dict(args)) or "{}")
    futures, ex = A._fan_out_subagents(
        [_tc(1, "general-purpose", isolation="worktree"),
         _tc(2, "general-purpose")], permissions=None)
    for f in futures.values():
        f.result()
    ex.shutdown(wait=True)
    assert all(c.get("isolation") == "worktree" for c in captured)
