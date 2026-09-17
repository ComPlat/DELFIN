"""Tests for delfin.agent.cli_stream — the rendering only.

Events are fabricated here in the exact shapes cli.py::_run_once
emits ({"type": "text"}, {"type": "tool_use"}) plus this module's
own tool_result / tick / turn_end. No engine is imported.
"""

from __future__ import annotations

from delfin.agent import cli_stream
from delfin.agent.cli_stream import (
    SPINNER_FRAMES,
    StreamRenderer,
    render_tool_call,
    render_turn_footer,
)


# --- render_tool_call -------------------------------------------------

def test_render_tool_call_tty_bold_with_args():
    line = render_tool_call(
        "read_file", {"path": "delfin/agent/cli.py", "limit": 40},
        {"is_tty": True, "elapsed_s": 1.2})
    assert line.startswith("\x1b[1mread_file")
    assert "path=delfin/agent/cli.py" in line
    assert "limit=40" in line
    assert "(1.2s)" in line
    assert "\x1b[0m" in line  # reset present


def test_render_tool_call_spinner_frame_in_tty():
    line = render_tool_call("bash", {"command": "ls"},
                            {"is_tty": True, "spinner": 3, "elapsed_s": 0.0})
    assert line.startswith(f"{SPINNER_FRAMES[3]} ")


def test_render_tool_call_not_a_tty_is_plain_and_greppable():
    line = render_tool_call("grep_file", {"pattern": "stream"},
                            {"is_tty": False, "elapsed_s": 2.5})
    assert line == "· grep_file pattern=stream (2.5s)"
    assert "\x1b" not in line


def test_render_tool_call_truncates_long_values():
    long_path = "x" * 100
    line = render_tool_call("read_file", {"path": long_path},
                            {"is_tty": False})
    assert "…" in line
    assert "x" * 100 not in line


def test_render_tool_call_without_elapsed():
    assert render_tool_call("ls", None, {"is_tty": False}) == "· ls"


def test_render_tool_call_flattens_newlines_in_values():
    line = render_tool_call("write_file", {"content": "a\nb"},
                            {"is_tty": False})
    assert "\n" not in line


# --- render_turn_footer -----------------------------------------------

def test_render_turn_footer_full_stats():
    line = render_turn_footer({"input_tokens": 100, "output_tokens": 20,
                               "tool_calls": 3, "duration_s": 12.34,
                               "is_tty": False})
    assert line == "→ in=100 · out=20 · tools=3 · 12.3s"


def test_render_turn_footer_tty_is_bold():
    line = render_turn_footer({"input_tokens": 1, "output_tokens": 2,
                               "tool_calls": 0, "duration_s": 1.0,
                               "is_tty": True})
    assert line.startswith("\x1b[1m→")


def test_render_turn_footer_missing_fields_omitted_not_zero():
    line = render_turn_footer({"is_tty": False})
    assert line == "→ turn ended"
    assert "in=" not in line and "tools=" not in line


# --- StreamRenderer ---------------------------------------------------

def _events():
    return [
        {"type": "text", "text": "Reading the file first."},
        {"type": "tool_use", "name": "read_file",
         "input": {"path": "delfin/agent/cli.py"}},
        {"type": "tool_result", "name": "read_file", "elapsed_s": 1.5},
        {"type": "tool_use", "name": "grep_file",
         "input": {"pattern": "stream"}},
        {"type": "tool_result", "name": "grep_file", "elapsed_s": 0.5},
        {"type": "text", "text": "Found it at cli.py:288."},
        {"type": "turn_end", "input_tokens": 512, "output_tokens": 64,
         "duration_s": 9.0},
    ]


def test_renderer_tty_yields_text_calls_and_bold_footer():
    clock = iter([0.0, 0.0, 10.0, 10.0, 19.0])  # only monotonic reads used
    lines = list(StreamRenderer(
        _events(), is_tty=True, now=lambda: next(clock, 99.0)))
    assert lines[0] == "Reading the file first."
    # tool_use line: spinner frame, bold name, arg summary
    assert lines[1].startswith(f"{SPINNER_FRAMES[0]} \x1b[1mread_file")
    assert "path=delfin/agent/cli.py" in lines[1]
    # closing line of the call carries the elapsed time
    assert "(1.5s)" in lines[2]
    assert "grep_file pattern=stream" in lines[3]
    assert "(0.5s)" in lines[4]
    assert lines[5] == "Found it at cli.py:288."
    assert lines[6].startswith("\x1b[1m→ in=512 · out=64 · tools=2 · 9.0s")


def test_renderer_not_a_tty_no_spinner_no_ansi_one_line_per_event():
    lines = list(StreamRenderer(_events(), is_tty=False))
    assert len(lines) == 7  # one per non-tick event
    for line in lines:
        assert "\x1b" not in line
    assert not any(l.startswith("⠋") for l in lines)
    assert lines[1] == "· read_file path=delfin/agent/cli.py"
    assert lines[2] == "· read_file (1.5s)"
    assert lines[6] == "→ in=512 · out=64 · tools=2 · 9.0s"


def test_renderer_counts_tool_calls_for_footer_fallback():
    events = _events()
    events[-1] = {"type": "turn_end", "duration_s": 5.0}
    lines = list(StreamRenderer(events, is_tty=False))
    assert lines[-1] == "→ tools=2 · 5.0s"


def test_renderer_tick_spins_only_while_a_call_is_open():
    events = [
        {"type": "tool_use", "name": "bash", "input": {"command": "pytest"}},
        {"type": "tick"},
        {"type": "tick"},
        {"type": "tool_result", "name": "bash", "elapsed_s": 3.0},
        {"type": "tick"},          # nothing open: must render nothing
        {"type": "turn_end", "duration_s": 4.0},
    ]
    clock = iter([0.0, 1.0, 2.0, 3.0, 4.0])
    lines = list(StreamRenderer(events, is_tty=True,
                                now=lambda: next(clock, 99.0)))
    spinner_lines = [l for l in lines if l.startswith(SPINNER_FRAMES[1])]
    # tick 1 -> frame[1]; tick 2 -> frame[2]; tick 3 (idle) -> nothing
    assert lines[1].startswith(f"{SPINNER_FRAMES[1]} \x1b[1mbash")
    assert lines[2].startswith(f"{SPINNER_FRAMES[2]} \x1b[1mbash")
    assert "(3.0s)" in lines[3]
    assert lines[4].startswith("\x1b[1m→ tools=1 · 4.0s")
    assert len(lines) == 5
    assert spinner_lines  # a spinner that is not a lie: frames did advance


def test_renderer_tty_tick_hidden_when_not_tty():
    events = [
        {"type": "tool_use", "name": "bash", "input": {}},
        {"type": "tick"},
        {"type": "tool_result", "name": "bash", "elapsed_s": 1.0},
        {"type": "turn_end", "duration_s": 1.0},
    ]
    lines = list(StreamRenderer(events, is_tty=False,
                                now=lambda: 0.0))
    assert lines == ["· bash", "· bash (1.0s)", "→ tools=1 · 1.0s"]


def test_renderer_empty_stream_yields_nothing():
    assert list(StreamRenderer([], is_tty=False)) == []


def test_renderer_default_is_tty_follows_stdout(monkeypatch):
    class FakeOut:
        def isatty(self):
            return True
    monkeypatch.setattr(cli_stream.sys, "stdout", FakeOut())
    r = StreamRenderer([{"type": "text", "text": "hi"}])
    assert r._is_tty is True


def test_no_engine_import():
    import sys
    mod = sys.modules["delfin.agent.cli_stream"]
    assert "delfin.agent.engine" not in getattr(mod, "__dict__", {})
