"""A greeting must not be answered with a pipeline report.

Field report 20260805-085306 (office, kit.qwen3.5-397b): the user typed a
greeting and got back "--- Cycle complete ---", the pipeline line, a
self-optimization note and "Continue chatting or /reset for new task." --
four lines of machinery around a one-sentence answer.

The discriminator is what the turn DID, not which mode it ran in. A
single-role office session and a multi-role code route both report their
cycle when a cycle actually ran, and both stay quiet when nothing did.

Suppressed is the chat noise only: the cycle event, the persistent cycle
memory and the self-optimization measurement are still recorded, so the
Cycle Inspector and the cost history keep working.
"""

from __future__ import annotations

import pathlib
import re

from delfin.dashboard.tab_agent import _turn_warrants_a_cycle_report

_SOURCE = pathlib.Path(
    __import__("delfin.dashboard.tab_agent", fromlist=["x"]).__file__
).read_text(encoding="utf-8")


def test_a_turn_without_a_tool_call_gets_no_ceremony():
    assert _turn_warrants_a_cycle_report(0) is False


def test_a_turn_that_did_work_still_reports():
    assert _turn_warrants_a_cycle_report(1) is True
    assert _turn_warrants_a_cycle_report(37) is True


def test_the_rule_does_not_look_at_the_mode():
    """It must read the same in office and in every other mode -- a rule
    keyed on the route would have to be remembered for each new mode."""
    import inspect

    src = inspect.getsource(_turn_warrants_a_cycle_report)
    body = src.split('"""')[-1]
    for leak in ("mode", "route", "office", "dashboard", "solo"):
        assert leak not in body, (
            f'the rule inspects {leak!r}; it must decide on work done alone')


# ---------------------------------------------------------------------------
# What must stay wired
# ---------------------------------------------------------------------------

def test_all_three_chat_lines_are_gated():
    """Cycle verdict, self-optimization note and the follow-up hint."""
    assert re.search(r"if _report_cycle:\s*\n\s*_append_system_message\(\s*\n\s*f\"--- Cycle complete",
                     _SOURCE), "the cycle verdict line is no longer gated"
    assert "if _opt_changes and _report_cycle:" in _SOURCE, (
        "the self-optimization note is no longer gated")
    assert re.search(r"if _report_cycle:\s*\n\s*_append_system_message\(\s*\n\s*f\"Continue chatting",
                     _SOURCE), "the follow-up hint is no longer gated"


def test_the_measurements_are_not_gated():
    """Suppressing the chat line must not suppress the recording behind it:
    the Cycle Inspector and the cost history depend on these firing."""
    decided = _SOURCE.index("_report_cycle = _turn_warrants_a_cycle_report")
    tail = _SOURCE[decided:decided + 3000]
    for call, what in (
        ('_record_cycle_event("cycle", "Cycle complete")', "the cycle event"),
        ("_check_acceptance_gate(engine)", "the acceptance verdict"),
        ("engine.record_cycle_outcome(", "the self-optimization measurement"),
    ):
        idx = tail.find(call)
        assert idx != -1, f"{what} vanished from the cycle-completion block"
        line_start = tail.rfind("\n", 0, idx) + 1
        assert not tail[line_start:idx].strip().startswith("if _report_cycle"), (
            f"{what} is now gated on the chat line; it must always run")


def test_the_follow_up_state_is_set_regardless():
    """The chat hint is cosmetic; the state change it describes is not."""
    decided = _SOURCE.index("_report_cycle = _turn_warrants_a_cycle_report")
    tail = _SOURCE[decided:decided + 4000]
    state_idx = tail.index('state["_follow_up"] = True')
    hint_idx = tail.index('f"Continue chatting or /reset for new task."')
    assert state_idx < hint_idx, (
        "the follow-up state must be set before -- and independently of -- "
        "the hint that announces it")


def test_a_user_triggered_completion_still_answers_the_user():
    """Gate approval and the manual handoff button are direct responses to
    something the user just clicked or typed; those stay unconditional."""
    for marker in ('f"--- Gate approved: {prev_label}',
                   'f"--- Manual handoff: {prev_role}'):
        assert marker in _SOURCE, f"{marker} disappeared"


# ---------------------------------------------------------------------------
# Per-turn telemetry belongs on the status row, not in the conversation
# ---------------------------------------------------------------------------

def test_the_status_row_carries_the_last_turn_cost():
    """Taking the chat lines out must not take the number away."""
    from delfin.dashboard.tab_agent import _render_status

    html = _render_status("office", "api", "office_agent", 0, 1,
                          77_109, 1_254, 0.164, provider="kit",
                          last_turn_cost_usd=0.1643)
    assert "0.164" in html, "the last turn's cost is not on the status row"
    assert "last turn" in html


def test_a_turn_that_cost_nothing_adds_nothing_to_the_status_row():
    from delfin.dashboard.tab_agent import _render_status

    html = _render_status("office", "api", "office_agent", 0, 1,
                          0, 0, 0.0, provider="kit", last_turn_cost_usd=0.0)
    assert "last turn" not in html


def test_a_tiny_turn_cost_is_not_rounded_to_zero():
    """Three decimals would print $0.000 for a cheap turn and read as free."""
    from delfin.dashboard.tab_agent import _render_status

    html = _render_status("office", "api", "office_agent", 0, 1,
                          900, 20, 0.0002, provider="kit",
                          last_turn_cost_usd=0.0002)
    assert "$0.000 " not in html
    assert "0.0002" in html


def test_neither_per_turn_line_is_posted_into_the_chat_any_more():
    assert "⏱ turn" not in _SOURCE, (
        'the clock line is back in the conversation')
    assert 'f"{role_label}: {_role_in:,} in / {_role_out:,} out"' not in _SOURCE, (
        'the per-role usage line is back in the conversation')


def test_cost_milestones_still_interrupt():
    """Running past a whole dollar is worth saying out loud; a turn is not."""
    assert 'f"Cost milestone: ${engine.cost_usd:.2f}' in _SOURCE


# ---------------------------------------------------------------------------
# No gap between the role label and the first line
# ---------------------------------------------------------------------------

def test_a_reply_opening_with_a_newline_does_not_start_with_a_break():
    """Every newline becomes a <br>, so a leading one reads as a blank line
    between the agent's name and its first sentence."""
    from delfin.dashboard.tab_agent import _md_to_html

    html = _md_to_html("\n\nIch lese alle Buchungsdaten aus.")
    assert not html.startswith("<br>"), html[:40]
    assert html.startswith("Ich lese")


def test_a_trailing_newline_does_not_pad_the_bottom():
    from delfin.dashboard.tab_agent import _md_to_html

    assert not _md_to_html("Fertig.\n\n").endswith("<br>")


def test_content_inside_the_reply_keeps_its_line_breaks():
    """Only the edges are trimmed -- paragraphs must survive."""
    from delfin.dashboard.tab_agent import _md_to_html

    assert "<br>" in _md_to_html("Zeile eins\nZeile zwei")


def test_the_streaming_view_trims_the_same_way():
    assert 'strip(chr(13) + chr(10))' in _SOURCE, (
        'the streaming branch renders raw content again, so the gap comes '
        'back for the whole time a reply is being written')


def test_the_role_label_sits_tight_against_its_line():
    """Both levers, not just the margin: default leading on an 11px label
    reopens the gap on its own."""
    block = _SOURCE[_SOURCE.index(".delfin-chat-role {"):]
    block = block[:block.index("}")]
    assert "margin-bottom: 1px" in block, block
    assert "line-height: 1.1" in block, block


# ---------------------------------------------------------------------------
# The pipeline line is chat, not a widget
# ---------------------------------------------------------------------------

def test_the_pipeline_line_is_gated_like_the_rest():
    """I first read _update_pipeline_display as a widget refresh and
    protected it from the gate. It posts a chat message, so a greeting was
    still answered with "Pipeline: [ok] Office Agent"."""
    decided = _SOURCE.index("_report_cycle = _turn_warrants_a_cycle_report")
    tail = _SOURCE[decided:decided + 1200]
    idx = tail.index("_update_pipeline_display(engine)")
    indent = len(tail[tail.rfind("\n", 0, idx) + 1:idx])
    gate = tail.rindex("if _report_cycle:", 0, idx)
    gate_indent = len(tail[tail.rfind("\n", 0, gate) + 1:gate])
    assert indent > gate_indent, (
        "the pipeline line is posted outside the gate again")


def test_a_route_of_one_is_not_announced_as_a_pipeline():
    """Single-role modes list the one agent that just answered."""
    assert "if len(steps) < 2:" in _SOURCE, (
        'a one-step route is announced as a pipeline again')


def test_the_length_rule_is_not_a_list_of_mode_names():
    """Named exclusions have to be extended for every new single-role mode;
    the length check is the rule behind them."""
    block = _SOURCE[_SOURCE.index("def _update_pipeline_display"):]
    block = block[:block.index("_append_system_message")]
    assert "office" not in block, (
        'the pipeline line excludes office by name instead of by shape')
