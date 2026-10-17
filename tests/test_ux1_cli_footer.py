"""UX1 phase 3: the CLI turn footer is told by turn_summary.

The footer a person sees after a turn in the CLI is the summary line:

    Surfaced after 31s · done 11:11 · 1 shell still running

built by ``delfin.agent.turn_summary.format_turn_summary`` — not the
old token counter. render_turn_footer stays the public entry point; it
maps its stats onto the summary's inputs (outcome, ended_at,
running_jobs) and keeps the token figures as a dimmed suffix when the
event carries them.
"""

import datetime

from delfin.agent.cli_stream import render_turn_footer
from delfin.agent.turn_summary import format_turn_summary


def test_footer_is_the_summary_line():
    stats = {"duration_s": 31, "ended_at": "2026-10-18T11:11:30",
             "outcome": "ok", "running_jobs": 1, "is_tty": False}
    line = render_turn_footer(stats)
    assert line == "Surfaced after 31s · done 11:11 · 1 shell still running"


def test_footer_without_time_or_jobs():
    line = render_turn_footer(
        {"duration_s": 31, "outcome": "ok", "is_tty": False})
    assert line == "Surfaced after 31s"


def test_footer_failure_never_gets_the_success_verb():
    line = render_turn_footer(
        {"duration_s": 12, "outcome": "stopped", "is_tty": False})
    assert line == "Stopped after 12s"


def test_footer_keeps_the_token_suffix_dimmed():
    stats = {"duration_s": 31, "outcome": "ok", "input_tokens": 10,
             "output_tokens": 20, "tool_calls": 2, "is_tty": False}
    line = render_turn_footer(stats)
    assert line.startswith("Surfaced after 31s · ")
    assert "in=10" in line and "out=20" in line and "tools=2" in line


def test_footer_matches_turn_summary_exactly():
    stats = {"duration_s": 133, "ended_at": "2026-10-18T11:11:30",
             "outcome": "ok", "running_jobs": 2, "is_tty": False}
    assert render_turn_footer(stats) == format_turn_summary(
        duration_s=133,
        ended_at=datetime.datetime.fromisoformat("2026-10-18T11:11:30"),
        outcome="ok", running_jobs=2)


def test_streamrenderer_error_turn_marks_the_outcome():
    # A turn that ended through an error event never renders "Surfaced":
    # the renderer maps what it saw onto the outcome it reports.
    from delfin.agent.cli_stream import StreamRenderer
    events = [
        {"type": "text", "text": "partial answer"},
        {"type": "error", "text": "boom"},
        {"type": "turn_end", "duration_s": 12},
    ]
    lines = list(StreamRenderer(iter(events), is_tty=False))
    assert lines[-1] == "Stopped after 12s"
