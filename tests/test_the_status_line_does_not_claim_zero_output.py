"""The status line does not claim a turn produced zero output tokens.

Measured on supervised runs (21./22.9.2026, kit.glm-5.3 via the
OpenAI-compatible endpoint): while a turn was running the status line
showed ``↑<n> ↓0`` the whole time, although ``turn_metrics`` recorded
30-500 output tokens for the same requests afterwards.

The number is genuinely unavailable DURING the turn on that endpoint:
usage arrives only with the last stream chunk (``include_usage``), the
engine adds it to ``token_usage`` on the final ``message_delta``, and
the status line reads the counters before that. A hard ``↓0`` is
therefore a false statement about the present, not a missing value.

Judged here:

  a mid-turn estimate, marked as one
                        when the engine reports no output tokens yet but
                        the pump has seen streamed text, the line shows
                        a marked estimate (``↓~<n>``), never ``↓0``

  nothing seen yet is honest about it
                        before the first token, ``↓0`` would again be a
                        false exact number; the line says ``↓…`` instead

  a real count still wins
                        an engine that DOES report output tokens
                        mid-turn (Anthropic-style message_delta) keeps
                        its exact ``↓<n>`` -- the estimate is a fallback,
                        not a replacement
"""

from __future__ import annotations

import time

from delfin.agent.repl import RenderItem
from delfin.agent.repl import TerminalAgent
from delfin.agent.repl import ReplOptions

from pathlib import Path


class _Engine:
    """Cumulative counters, exactly what get_status() reports."""

    session_id = "statusline-0002"
    token_usage = {"input": 0, "output": 0}
    last_turn_stop_reason = ""

    def __init__(self, *, output_tokens: int):
        self._output = output_tokens
        self.messages = []

    def get_status(self):
        return {"input_tokens": 10, "output_tokens": self._output}


def _agent(tmp_path, *, output_tokens):
    import io

    engine = _Engine(output_tokens=output_tokens)
    out, err = io.StringIO(), io.StringIO()
    out.isatty = lambda: False
    err.isatty = lambda: False
    agent = TerminalAgent(
        engine, opts=ReplOptions(cwd=Path(tmp_path), max_tokens=0),
        out=out, err=err, read_line=lambda _p="": "exit")
    agent._turn_active.set()
    agent._turn_t0 = time.monotonic()
    agent._turn_base = (0, 0, 0.0)
    return agent


def _status_of(agent) -> str:
    line = agent._status_line()
    plain = line
    # The theme may wrap the line in ANSI; the arrow segment must be
    # findable regardless.
    for tag in ("↓",):
        assert tag in plain or "\u2193" in plain, (
            f"no output-token segment in status line: {line!r}")
    return line


def test_streamed_text_shows_a_marked_estimate_not_zero(tmp_path):
    """Endpoint reports no output tokens mid-turn (usage comes with the
    last chunk), but the turn has streamed text: the line must not
    claim ``↓0``."""
    agent = _agent(tmp_path, output_tokens=0)
    agent._count_streamed(RenderItem(
        "text", text="word " * 80))          # 400 chars -> ~100 tokens
    line = _status_of(agent)
    assert "↓~" in line, (
        f"a turn with streamed text must show a MARKED estimate, got: "
        f"{line!r}")
    assert "↓0" not in line, (
        f"↓0 is a false exact claim while tokens are streaming: {line!r}")


def test_nothing_seen_yet_shows_an_ellipsis_not_zero(tmp_path):
    """Before the first token nothing is known, not even for an
    estimate: ``↓…`` instead of a false ``↓0``."""
    agent = _agent(tmp_path, output_tokens=0)
    line = _status_of(agent)
    assert "↓…" in line, (
        f"before the first token the output count is UNKNOWN, not zero: "
        f"{line!r}")


def test_a_real_mid_turn_count_stays_exact(tmp_path):
    """An engine that reports output tokens during the turn keeps the
    exact number; the estimate never overrides provider truth."""
    agent = _agent(tmp_path, output_tokens=123)
    agent._count_streamed(RenderItem(
        "text", text="word " * 80))
    line = _status_of(agent)
    assert "↓123" in line, (
        f"provider-reported output tokens must win over the estimate: "
        f"{line!r}")
