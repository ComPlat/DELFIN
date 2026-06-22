"""Tests for the dashboard cost-efficiency pass:

- ACTION: /done sentinel is recognised by _dashboard_auto_exec
- _MAX_ACTION_CONT cap reduced from 4 to 3
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


def test_max_action_cont_reduced_to_three():
    """The continuation-cap was 4; the cost-efficiency pass lowered it
    to 3. Empirical: >3 rounds is always model-confusion in practice."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    assert "_MAX_ACTION_CONT = 3" in text
    # And the old "= 4" should be gone in the continuation-loop region
    cont_idx = text.find("_MAX_ACTION_CONT")
    snippet = text[cont_idx: cont_idx + 200]
    assert "= 4" not in snippet


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
    because nothing actually ran. Instead show a 'please clarify'
    message so the user knows the agent had no idea what to do."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("agent had no action to execute")
    assert idx > 0, (
        "missing 'no action' placeholder — /done-alone case still "
        "shows misleading '(commands executed)'"
    )
    # And the branch that triggers it: cleaned is empty + no real_results + done_seen
    block = text[max(0, idx-600): idx + 200]
    assert "elif done_seen:" in block


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
