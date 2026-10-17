"""UX1 review-fixup tests (s37 findings on f05be3a4).

(1) inf duration must not crash. (2) the footer in a pipe stays
ANSI-free. (3) an error line collapses whitespace and strips escapes.
"""

import datetime

from delfin.agent.cli_stream import StreamRenderer, render_turn_footer
from delfin.agent.turn_summary import format_duration, format_turn_summary


def test_inf_duration_renders_nothing_not_crash():
    assert format_duration(float("inf")) == ""
    assert format_duration(float("nan")) == ""
    line = format_turn_summary(duration_s=float("inf"), ended_at=None,
                               outcome="ok")
    assert "after" not in line


def test_footer_pipe_has_no_ansi():
    stats = {"duration_s": 1.0, "outcome": "ok", "input_tokens": 1,
             "output_tokens": 2, "tool_calls": 1, "is_tty": False}
    line = render_turn_footer(stats)
    assert "\x1b" not in line
    assert line == "Surfaced after 1s · in=1 · out=2 · tools=1"


def test_footer_tty_keeps_dimmed_suffix():
    stats = {"duration_s": 1.0, "outcome": "ok", "input_tokens": 1,
             "tool_calls": 1, "is_tty": True}
    line = render_turn_footer(stats)
    assert line == "Surfaced after 1s · \x1b[2min=1 · tools=1\x1b[0m"


def test_error_line_collapses_whitespace_and_strips_ansi():
    events = [
        {"type": "error", "text": "boom\x1b[31m!\nsecond  line"},
        {"type": "turn_end", "duration_s": 12},
    ]
    lines = list(StreamRenderer(iter(events), is_tty=False))
    assert lines[0] == "⚠ boom! second line"
    assert "\x1b" not in lines[0]
    assert lines[-1] == "Stopped after 12s"


def test_error_turn_still_reports_stopped_after_sanitising():
    # The outcome rule survives sanitisation: an error event maps the
    # turn_end footer onto "Stopped" even when its text was mangled.
    events = [
        {"type": "error", "text": "\x1b[1mfail"},
        {"type": "turn_end", "duration_s": 3},
    ]
    lines = list(StreamRenderer(iter(events), is_tty=False))
    assert lines[-1] == "Stopped after 3s"


def test_ended_at_datetime_still_renders_done_hhmm():
    # Phase-2 pinned behaviour must survive the fixup: a datetime
    # ended_at renders "done HH:MM" exactly as before.
    dt = datetime.datetime(2026, 10, 18, 11, 11, 30)
    assert format_turn_summary(duration_s=31, ended_at=dt,
                               outcome="ok") == \
        "Surfaced after 31s · done 11:11"
