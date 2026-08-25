"""Harness speech and answer text shared one channel, with no way to tell.

``on_token`` carried both the model's answer and the harness talking about
itself — retry banners, the stop notice, the cost ceiling, a blocking
hook. Every consumer that treated a notice as answer text drew a wrong
conclusion from it; the engine's own comment records a benchmark run whose
entire recorded output was three retry banners, scored as a model that
answered badly and written into the file baselines compare against.

A terminal has the same problem one layer further out, and worse: it
prints. Machinery speech lands in the middle of a sentence, and
``delfin-agent -p "..." > answer.txt`` captures the banner as if it were
the answer.

The second half of this file is the mirror image on the tool side. A
blocked call and a successful one were indistinguishable to any UI,
because ``on_tool_result`` hands over a 2000-character head slice and no
verdict — even though the engine computes the verdict two lines later for
its own trace.
"""

from __future__ import annotations

import io
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent import repl, repl_render as rr
from delfin.agent.api_client import StreamEvent


@pytest.fixture
def eng(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")


def _stream(*events):
    def _gen(**kwargs):
        for event in events:
            yield event
        yield StreamEvent(type="message_delta", input_tokens=5,
                          output_tokens=7, cost_usd=0.0)
    return _gen


# ---------------------------------------------------------------------------
# The seam
# ---------------------------------------------------------------------------

def test_a_notice_and_the_answer_arrive_on_different_channels(eng):
    eng.client.stream_message = _stream(
        StreamEvent(type="text_delta", text="part one "),
        StreamEvent(type="notice", text="retrying after 429"),
        StreamEvent(type="text_delta", text="part two"),
    )
    answer: list[str] = []
    notices: list[str] = []
    eng.stream_response(
        "go", on_token=answer.append, on_notice=notices.append)

    assert "".join(answer) == "part one part two"
    assert notices == ["retrying after 429"]


def test_without_the_new_channel_nothing_changes_for_old_callers(eng):
    """The non-breaking guarantee.

    The dashboard, cli._run_once, the benchmark and the scheduler all pass
    no on_notice. They have to keep seeing exactly what they saw.
    """
    eng.client.stream_message = _stream(
        StreamEvent(type="text_delta", text="answer"),
        StreamEvent(type="notice", text="retrying after 429"),
    )
    seen: list[str] = []
    eng.stream_response("go", on_token=seen.append)
    assert "retrying after 429" in seen


def test_a_notice_is_still_not_counted_as_the_answer(eng):
    """It was never appended to `chunks`, and must stay that way."""
    eng.client.stream_message = _stream(
        StreamEvent(type="notice", text="retrying after 429"),
        StreamEvent(type="text_delta", text="the real answer"),
    )
    out = eng.stream_response("go", on_notice=lambda _t: None)
    assert out.strip() == "the real answer"


# ---------------------------------------------------------------------------
# The tool side
# ---------------------------------------------------------------------------

def test_a_blocked_tool_is_not_reported_as_a_successful_one(eng):
    eng.client.stream_message = _stream(
        StreamEvent(type="tool_result", tool_name="write_file",
                    tool_output='{"error": "path is on the deny-list"}'),
        StreamEvent(type="text_delta", text="ok"),
    )
    meta: list[tuple[str, dict]] = []
    eng.stream_response("go", on_tool_result_meta=lambda n, m: meta.append((n, m)))

    assert meta, "no meta was emitted for the tool result"
    name, payload = meta[0]
    assert name == "write_file"
    assert payload["ok"] is False
    assert "deny-list" in payload["error"]


def test_a_result_that_worked_says_so(eng):
    eng.client.stream_message = _stream(
        StreamEvent(type="tool_result", tool_name="read_file",
                    tool_output="line one\nline two\n"),
        StreamEvent(type="text_delta", text="ok"),
    )
    meta: list[dict] = []
    eng.stream_response("go", on_tool_result_meta=lambda n, m: meta.append(m))
    assert meta and meta[0]["ok"] is True
    assert meta[0]["error"] == ""


def test_a_truncated_result_carries_its_true_size(eng):
    eng.client.stream_message = _stream(
        StreamEvent(type="tool_result", tool_name="read_file",
                    tool_output="head slice only",
                    output_truncated=True, output_chars=100_000),
        StreamEvent(type="text_delta", text="ok"),
    )
    meta: list[dict] = []
    eng.stream_response("go", on_tool_result_meta=lambda n, m: meta.append(m))
    assert meta[0]["truncated"] is True
    assert meta[0]["chars"] == 100_000


def test_the_old_tool_callback_is_untouched(eng):
    eng.client.stream_message = _stream(
        StreamEvent(type="tool_result", tool_name="read_file",
                    tool_output="body"),
        StreamEvent(type="text_delta", text="ok"),
    )
    seen: list[tuple[str, str]] = []
    eng.stream_response("go", on_tool_result=lambda n, o: seen.append((n, o)))
    assert seen == [("read_file", "body")]


# ---------------------------------------------------------------------------
# What the terminal does with it
# ---------------------------------------------------------------------------

def test_the_banner_never_reaches_stdout(eng):
    """The reason the seam exists at all, asserted at the streams."""
    eng.client.stream_message = _stream(
        StreamEvent(type="text_delta", text="part one "),
        StreamEvent(type="notice", text="retrying after 429"),
        StreamEvent(type="text_delta", text="part two"),
    )
    out, err = io.StringIO(), io.StringIO()
    transcript = repl.Transcript(out, err, theme=rr.Theme(enabled=False))
    result = repl.run_turn(eng, "go", sink=transcript.render)
    transcript.finish()

    assert out.getvalue() == "part one part two\n"
    assert "retrying after 429" in err.getvalue()
    assert "retrying" not in out.getvalue()
    assert result.text == "part one part two"


def test_a_blocked_call_is_drawn_as_blocked(eng):
    eng.client.stream_message = _stream(
        StreamEvent(type="tool_result", tool_name="write_file",
                    tool_output='{"error": "path is on the deny-list"}'),
        StreamEvent(type="text_delta", text="ok"),
    )
    out, err = io.StringIO(), io.StringIO()
    transcript = repl.Transcript(out, err, theme=rr.Theme(enabled=False))
    repl.run_turn(eng, "go", sink=transcript.render)

    assert "blocked" in err.getvalue()
    assert "blocked" not in out.getvalue()


def test_one_tool_result_draws_one_line(eng):
    """Found by running it, not by reading it.

    A tool result arrives on TWO callbacks — the head slice on
    on_tool_result, the verdict and the true size on on_tool_result_meta —
    and rendering each as it came drew every call twice, the second time
    with a body it did not have (`9+ lines` then `0 lines`). Two stubs in
    two unit tests each looked right on its own.
    """
    eng.client.stream_message = _stream(
        StreamEvent(type="tool_result", tool_name="read_file",
                    tool_output="a\nb\nc"),
        StreamEvent(type="tool_result", tool_name="grep_file",
                    tool_output="hit"),
        StreamEvent(type="text_delta", text="ok"),
    )
    out, err = io.StringIO(), io.StringIO()
    transcript = repl.Transcript(out, err, theme=rr.Theme(enabled=False))
    repl.run_turn(eng, "go", sink=transcript.render)

    results = [ln for ln in err.getvalue().splitlines() if "⎿" in ln]
    assert len(results) == 2, f"expected one line per result, got {results}"
    # And with both halves merged the count is exact, not a floor.
    assert "+" not in results[0]


def test_a_result_with_no_verdict_is_still_drawn():
    """The merge must not be able to swallow a line.

    If a path ever reports a result without the meta callback, holding the
    text for a verdict that never comes would lose it entirely.
    """
    seen: list[repl.RenderItem] = []

    class _NoMeta:
        token_usage = {"input": 0, "output": 0}

        def stream_response(self, **kw):
            kw["on_tool_result"]("read_file", "body")
            return "done"

    repl.run_turn(_NoMeta(), "go", sink=seen.append)
    results = [i for i in seen if i.kind == "tool_result"]
    assert len(results) == 1
    assert results[0].text == "body"


def test_a_remote_tool_says_which_server_it_ran_on():
    """mcp__delfin-docs__read_file is unreadable; the server still matters.

    It is the difference between a call DELFIN gated and one that ran
    somewhere else, outside the workspace sandbox and outside the bash
    deny-list.
    """
    assert rr.short_tool_name(
        "mcp__delfin-docs__read_file") == "delfin-docs:read_file"
    assert rr.short_tool_name("read_file") == "read_file"
    line = rr.tool_headline("mcp__delfin-docs__read_file",
                            {"path": "calc.py"}, theme=rr.Theme(enabled=False))
    assert "delfin-docs:read_file" in line
    assert "mcp__" not in line
    assert "calc.py" in line, "the argument lookup must survive the rename"


def test_the_interactive_turn_reports_what_the_headless_turn_reports():
    """One JSON contract, not two.

    cli._run_once and repl.run_turn are separate implementations on
    purpose — changing _run_once's signature would break the injected
    doubles in the benchmark and scheduler tests. Separate implementations
    need this assertion instead.
    """
    from delfin.agent import cli as agent_cli

    class _Stub:
        token_usage = {"input": 0, "output": 0}

        def stream_response(self, **kwargs):
            return "x"

    headless = set(agent_cli._run_once(_Stub(), "hi"))
    interactive = set(repl.TurnResult().to_dict())
    assert headless == interactive == set(repl.TURN_KEYS)
