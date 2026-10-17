"""UX1 phase 3: the CLI turn footer is told by turn_summary.

At turn end the CLI renders two lines: the UX1 summary

    Surfaced after 31s · done 11:11 · 1 shell still running

from ``delfin.agent.turn_summary.format_turn_summary`` (via
``render_turn_summary_line``), then the legacy stats footer
("→ in=100 · out=20 · tools=3 · 12.3s", bold head on a TTY) unchanged
from before UX1 — its format is a pinned contract.
"""

import datetime

from delfin.agent.cli_stream import (StreamRenderer, render_turn_footer,
                                     render_turn_summary_line)
from delfin.agent.turn_summary import format_turn_summary


def test_summary_line_is_the_formatter_output():
    stats = {"duration_s": 31, "ended_at": "2026-10-18T11:11:30",
             "outcome": "ok", "running_jobs": 1}
    line = render_turn_summary_line(stats)
    assert line == "Surfaced after 31s · done 11:11 · 1 shell still running"


def test_summary_line_without_time_or_jobs():
    line = render_turn_summary_line({"duration_s": 31, "outcome": "ok"})
    assert line == "Surfaced after 31s"


def test_summary_line_failure_never_gets_the_success_verb():
    line = render_turn_summary_line({"duration_s": 12, "outcome": "stopped"})
    assert line == "Stopped after 12s"


def test_summary_line_is_plain_even_on_tty():
    # The summary is the person-facing line; no ANSI escapes in it, and
    # a pipe sees exactly the same characters.
    stats = {"duration_s": 31, "outcome": "ok", "is_tty": True}
    assert render_turn_summary_line(stats) == "Surfaced after 31s"
    assert "\x1b" not in render_turn_summary_line(stats)


def test_footer_is_the_legacy_stats_line():
    # The legacy footer keeps its exact pinned format (and its bold
    # head on a TTY — pinned by tests/test_cli_stream.py).
    stats = {"duration_s": 31, "outcome": "ok", "input_tokens": 10,
             "output_tokens": 20, "tool_calls": 2, "is_tty": False}
    assert render_turn_footer(stats) == \
        "→ in=10 · out=20 · tools=2 · 31.0s"
    stats["is_tty"] = True
    assert render_turn_footer(stats) == \
        "\x1b[1m→ in=10 · out=20 · tools=2 · 31.0s\x1b[0m"


def test_streamrenderer_error_turn_marks_the_outcome():
    # A turn that ended through an error event never renders "Surfaced":
    # the renderer maps what it saw onto the outcome it reports.
    events = [
        {"type": "text", "text": "partial answer"},
        {"type": "error", "text": "boom"},
        {"type": "turn_end", "duration_s": 12},
    ]
    lines = list(StreamRenderer(iter(events), is_tty=False))
    assert lines[-1] == "→ 12.0s"
    assert lines[-2] == "Stopped after 12s"
