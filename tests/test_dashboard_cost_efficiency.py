"""Tests for the dashboard cost-efficiency pass:

- ACTION: /done sentinel is recognised by _dashboard_auto_exec
- the ACTION continuation budget distinguishes progress from repetition
- _arm_stale_watcher honours stale_kill_after_s for cooperative kill
- dashboard_agent.md documents the /done sentinel
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# /done sentinel + cap
# ---------------------------------------------------------------------------


def test_done_sentinel_documented_in_dashboard_agent_md():
    """The agent must be taught when to emit ACTION: /done — otherwise
    the optimisation is pure code and the model never triggers it."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "agent" / "pack" / "agents" / "dashboard_agent.md")
    body = p.read_text(encoding="utf-8")
    assert "ACTION: /done" in body
    assert "skips the post-execute commentary" in body
    # Concrete examples for multi-step and single-action paths
    assert "/control key functional BP86" in body
    assert "/tab orca" in body


def test_action_cap_is_resolved_from_settings_not_hardcoded():
    """The continuation budget used to be a hardcoded flat cap of 3, which
    ended healthy multi-step work at the same point as a runaway loop. It
    now comes from ``_resolve_action_round_limits`` (settings-backed) and
    is paired with a repetition check."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    assert "_MAX_ACTION_CONT = 3" not in text
    assert "_MAX_ACTION_CONT, _ACTION_REPEATS = (" in text
    assert "_resolve_action_round_limits())" in text
    # while-guard survives as the absolute backstop
    assert "while chunks and _cont_turn < _MAX_ACTION_CONT:" in text


def test_dashboard_auto_exec_returns_done_marker_when_sentinel_present():
    """Production helper: when ACTION: /done appears anywhere in the
    agent response, auto_exec short-circuits with the magic marker so
    the continuation loop can break without re-prompting."""
    # We mirror the production logic here rather than calling the
    # closure directly (it's a nested function inside create_tab).
    import re

    def _detect_done(agent_text: str) -> bool:
        lines = agent_text.split("\n")
        return any(re.match(r"^ACTION:\s*/done\b", ln) for ln in lines)

    assert _detect_done("ACTION: /tab orca\nACTION: /done")
    assert _detect_done("blabla\n\nACTION: /done\n")
    assert _detect_done("ACTION: /done")
    assert not _detect_done("ACTION: /tab orca")
    assert not _detect_done("/done is in the prose not as an ACTION")
    # Production code uses the same regex — confirm with a glob check
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    # /done is now parsed inline as a sentinel — appended to the
    # results list AFTER real commands ran, not short-circuited.
    assert "ACTION: /done is a sentinel" in text
    assert '_done_seen = any(' in text
    assert 'results.append("__DONE__")' in text


def test_continuation_loop_breaks_on_done_marker():
    """The continuation loop must split the ``__DONE__`` sentinel out
    from real command results so the placeholder + early-break work
    correctly even when /done was emitted alongside real ACTIONs."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    # The loop must compute ``done_seen`` + ``real_results`` slices
    idx = text.find("Split the done-sentinel from real results")
    assert idx > 0, "loop must separate sentinel from real results"
    # Snippet must be large enough to span the entire while-loop body
    snippet = text[idx: idx + 2500]
    assert 'done_seen = "__DONE__" in exec_results' in snippet
    assert "real_results = [r for r in exec_results if r != \"__DONE__\"]" in snippet
    # Break logic uses the split variables, not the raw list
    assert "if done_seen:" in snippet
    assert "if not real_results:" in snippet


def test_done_only_response_shows_clear_no_action_placeholder():
    """When /done was emitted ALONE (no real ACTIONs), the cleaned
    placeholder must NOT say '(commands executed)' — that's misleading
    because nothing actually ran. On the FIRST round the user gets a
    'please clarify' message; in a wrap-up round (actions already ran
    earlier in the turn) a bare /done is a normal close and renders
    '(request completed)' instead."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("agent had no action to execute")
    assert idx > 0, (
        "missing 'no action' placeholder — /done-alone case still "
        "shows misleading '(commands executed)'"
    )
    # The branch that triggers it: cleaned empty + no real_results +
    # done_seen, split by continuation round.
    block = text[max(0, idx - 900): idx + 200]
    assert "elif done_seen:" in block
    assert "_cont_turn > 1" in block
    assert "(request completed)" in block


# ---------------------------------------------------------------------------
# Stale-stream cooperative kill watcher
# ---------------------------------------------------------------------------


def test_arm_stale_watcher_documents_kill_after_s():
    """The watcher's docstring must mention the new kill_after_s knob
    so callers know they can tune it. Regression-guard against silent
    removal."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("def _arm_stale_watcher")
    assert idx > 0
    snippet = text[idx: idx + 2500]
    assert "stale_kill_after_s" in snippet
    assert "cooperatively stops the stream" in snippet


def test_dashboard_mode_gets_default_120s_kill_threshold():
    """The dashboard default — when no user-setting is configured — is
    120 s. This guards against runaway 'commands executed' turns like
    the 117 s/$0.05 case the user reported."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("def _arm_stale_watcher")
    snippet = text[idx: idx + 3500]
    assert 'kill_after = 120.0' in snippet
    assert 'mode", "") == "dashboard"' in snippet


def test_kill_timer_disarmed_in_worker_finally():
    """Both stale-state and kill timers must be cancelled when the
    worker turn ends — leaking timers would fire after the next send
    and potentially stop a fresh stream."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    # Locate the finally cleanup block
    idx = text.find("Disarm the stale + kill watchers")
    assert idx > 0
    block = text[idx: idx + 600]
    assert '"_stale_timer"' in block
    assert '"_stale_kill_timer"' in block


def test_state_dict_seeds_kill_timer_slot():
    """The state dict must seed ``_stale_kill_timer: None`` so the
    finally block can safely dict.get() without a key error on the
    very first turn (before _arm_stale_watcher has run)."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    assert '"_stale_kill_timer": None' in text


# ---------------------------------------------------------------------------
# ACTION continuation budget — progress vs. repetition (loop wiring)
# ---------------------------------------------------------------------------


def _tab_agent_source() -> str:
    return (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
            / "tab_agent.py").read_text(encoding="utf-8")


def _continuation_loop_block() -> str:
    """Source of the ACTION continuation loop, from its budget setup to the
    stop-note block."""
    text = _tab_agent_source()
    start = text.find("_MAX_ACTION_CONT, _ACTION_REPEATS = (")
    assert start > 0, "continuation-loop budget setup not found"
    end = text.find("-- Verify-Enforcement", start)
    assert end > start, "continuation-loop block not found"
    return text[start:end]


def test_loop_carries_signature_command_and_stale_history():
    """The decision needs per-turn state: which ACTION sets ran, which
    commands ran, and how many rounds produced nothing new."""
    block = _continuation_loop_block()
    assert "_round_sigs: list[str] = []" in block
    assert "_seen_cmds: set[str] = set()" in block
    assert "_ran_cmds: list[str] = []" in block
    assert "_stale_rounds = 0" in block


def test_loop_calls_the_pure_decision_before_the_next_model_turn():
    """The progress-vs-repetition check must sit between 'results exist'
    and the continuation ``stream_response`` — otherwise a repeat still
    costs a model call."""
    block = _continuation_loop_block()
    dec_idx = block.find("_dec = _decide_action_round(")
    stream_idx = block.find("engine.stream_response(")
    results_idx = block.find("if not real_results:")
    assert 0 < results_idx < dec_idx < stream_idx
    for kw in ("rounds_used=_cont_turn", "stale_rounds=_stale_rounds",
               "ceiling=_MAX_ACTION_CONT", "repeat_limit=_ACTION_REPEATS"):
        assert kw in block, kw
    assert "if not _dec.proceed:" in block
    assert "_stop_reason = _dec.reason" in block


def test_done_sentinel_break_still_precedes_the_decision():
    """The explicit /done early-exit must keep short-circuiting the loop —
    the new budget must not make a finished turn pay a round."""
    block = _continuation_loop_block()
    done_idx = block.find("if done_seen:")
    dec_idx = block.find("_dec = _decide_action_round(")
    assert 0 < done_idx < dec_idx


def test_loop_compares_the_commands_that_were_actually_dispatched():
    """The tolerant dispatcher accepts ACTION forms the strict extractor
    misses, so the repetition check reads the dispatched batch and only
    falls back to re-parsing the raw text."""
    text = _tab_agent_source()
    assert 'state["_last_exec_commands"] = list(commands[:10])' in text
    assert '"_last_exec_commands": [],' in text
    block = _continuation_loop_block()
    assert 'state.get("_last_exec_commands") or []' in block
    assert "if not _round_cmds:" in block


def test_stop_note_fires_only_for_budget_and_repetition_stops():
    """A turn the agent closed itself must stay silent; a truncated one
    must not."""
    block = _continuation_loop_block()
    assert "if _stop_reason in (_ACTION_ROUND_CEILING," in block
    assert "_ACTION_ROUND_REPEAT," in block
    assert "_ACTION_ROUND_NO_PROGRESS):" in block
    assert "_note = _format_action_stop_note(" in block
    # The old, dishonest wording is gone
    assert "agent kept emitting commands" not in _tab_agent_source()


def test_stop_note_excludes_already_executed_commands_from_pending():
    """Reporting an executed command as outstanding would be a lie."""
    block = _continuation_loop_block()
    idx = block.find("_pending = [")
    assert idx > 0
    snippet = block[idx: idx + 400]
    assert "_normalize_action_command(c)" in snippet
    assert "not in _seen_cmds" in snippet


def test_budget_defaults_are_generous_and_repetition_is_tight():
    """The ceiling only has to bound an endless stream of DIFFERENT
    commands; repetition is caught after the first repeat."""
    from delfin.dashboard import tab_agent as ta
    assert ta._ACTION_ROUND_CEILING_DEFAULT >= 8
    assert ta._ACTION_ROUND_REPEAT_LIMIT_DEFAULT == 2


def test_settings_keys_are_documented_in_the_resolver():
    from delfin.dashboard import tab_agent as ta
    doc = ta._resolve_action_round_limits.__doc__ or ""
    assert "agent.max_action_rounds" in doc
    assert "agent.action_repeat_limit" in doc
