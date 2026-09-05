"""There are three ways a turn goes wrong, and the report named one.

A turn can DIE (an exception, and the record carries an error), it can
never START (the backend stayed silent and no token ever arrived), or it
can STALL (something arrived, late). They are told apart in the log by
design, and the roll-up a person reads counted crashes and never-started
turns without printing either -- so a week of dying turns and a week of
silent endpoints both rendered as a healthy week with a slightly odd
average.

The averages had the matching problem one step further on: a mean over
three first-token samples out of ninety turns is printed identically to a
mean over ninety, and the difference is whether the number means
anything.
"""

from __future__ import annotations

import pytest

from delfin.agent import turn_metrics as tm


def _recorded(module) -> list[dict]:
    """Every turn entry the log holds, oldest first."""
    import json
    out: list[dict] = []
    for fp in sorted(module._DIR.glob("*.jsonl")):
        for line in fp.read_text(encoding="utf-8").splitlines():
            try:
                out.append(json.loads(line))
            except Exception:
                continue
    return out


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setattr(tm, "_DIR", tmp_path / "turns")
    return tmp_path


# ---------------------------------------------------------------------------
# The roll-up can express all three
# ---------------------------------------------------------------------------

def test_a_silent_backend_is_counted_and_is_not_a_crash(home):
    """No error, no token, wall clock spent."""
    tm.record("s", model="m", total_ms=30_000, ttft_ms=None,
              output_chars=0, tool_calls=0, output_tokens=0)
    s = tm.aggregate_turn_stats()
    assert s["never_started"] == 1
    assert s["crashes"] == 0


def test_a_turn_that_died_is_a_crash_and_not_a_silent_backend(home):
    """The same empty record, and an error in it. Without the
    distinction, every crash that took more than twenty seconds to fail
    was absorbed into the stall count."""
    tm.record("s", model="m", total_ms=30_000, ttft_ms=None,
              error="APIError: upstream refused", output_chars=0)
    s = tm.aggregate_turn_stats()
    assert s["crashes"] == 1
    assert s["never_started"] == 0


def test_the_healthy_turn_is_neither(home):
    tm.record("s", model="m", total_ms=4000, ttft_ms=800,
              output_chars=500, output_tokens=120)
    s = tm.aggregate_turn_stats()
    assert (s["crashes"], s["never_started"], s["stalls"]) == (0, 0, 0)


def test_the_ttft_average_says_how_many_turns_it_covers(home):
    """Three answered turns among many is a different statement from
    three answered turns among three."""
    tm.record("s", model="m", total_ms=1000, ttft_ms=500, output_chars=10)
    for _ in range(4):
        tm.record("s", model="m", total_ms=30_000, ttft_ms=None,
                  output_chars=0, output_tokens=0)
    s = tm.aggregate_turn_stats()
    assert s["turns"] == 5
    assert s["ttft_sample"] == 1
    assert s["never_started"] == 4


# ---------------------------------------------------------------------------
# ... and the report a person reads says so
# ---------------------------------------------------------------------------

def test_the_health_section_names_every_kind(home, monkeypatch):
    from delfin.agent import eval_loop

    monkeypatch.setattr(
        "delfin.agent.turn_metrics.aggregate_turn_stats",
        lambda *a, **k: {"turns": 90, "avg_ttft_ms": 800, "p90_ttft_ms": 1200,
                         "stalls": 2, "crashes": 7, "stopped_count": 1,
                         "never_started": 11, "ttft_sample": 3})
    text = "\n".join(eval_loop._turn_health_lines())
    assert "crashes: 7" in text
    assert "never started: 11" in text
    assert "over 3 of 90 turns" in text


def test_a_quiet_week_still_renders(home, monkeypatch):
    """The half that keeps the section readable: nothing wrong, nothing
    shouted."""
    from delfin.agent import eval_loop

    monkeypatch.setattr(
        "delfin.agent.turn_metrics.aggregate_turn_stats",
        lambda *a, **k: {"turns": 40, "avg_ttft_ms": 700, "p90_ttft_ms": 900,
                         "stalls": 0, "crashes": 0, "stopped_count": 0,
                         "never_started": 0, "ttft_sample": 40})
    text = "\n".join(eval_loop._turn_health_lines())
    assert "crashes: 0" in text and "never started: 0" in text
    assert "over 40 of 40 turns" in text


# ---------------------------------------------------------------------------
# A turn that only thought still started
# ---------------------------------------------------------------------------

def test_a_reasoning_token_counts_as_the_backend_starting(tmp_path,
                                                          monkeypatch):
    """A turn that reasoned for two minutes and then answered recorded no
    first-token stamp at all, which is the record a backend that said
    nothing writes -- and on an endpoint that reports no usage there was
    nothing left to tell them apart."""
    import textwrap
    from unittest.mock import MagicMock, patch
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# solo mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))

    monkeypatch.setattr(tm, "_DIR", tmp_path / "turns")

    fake = MagicMock()
    fake._observed_files_session = set()
    fake.stream_message = MagicMock(side_effect=lambda *a, **k: iter([
        StreamEvent(type="thinking_delta", text="weighing the options…"),
        StreamEvent(type="text_delta", text="Die Antwort ist 42."),
        StreamEvent(type="message_delta", output_tokens=7, cost_usd=0.0),
    ]))
    with patch("delfin.agent.engine.create_client", return_value=fake):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick",
                             pack_dir=tmp_path)
    engine.stream_response("frag was")

    # Read back what the LOG holds -- that record is the only evidence
    # left after the turn, and it is what the health report reads.
    entries = _recorded(tm)
    assert entries, "the turn was not logged at all"
    assert entries[-1].get("ttft_ms") is not None, \
        "a reasoning token is the backend starting"
    assert not tm.is_never_started(entries[-1])


def test_a_turn_that_only_reasoned_is_not_a_silent_backend(tmp_path,
                                                           monkeypatch):
    """The case that actually isolates it. With a text delta anywhere in
    the turn the stamp comes from there whatever thinking does, so a test
    carrying both proves nothing about thinking -- mine did, and the
    mutation check said so. Here the model reasons and stops: on an
    endpoint that reports no usage, that record was indistinguishable
    from one written by a backend that never spoke."""
    import textwrap
    from unittest.mock import MagicMock, patch
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# solo mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    monkeypatch.setattr(tm, "_DIR", tmp_path / "turns")

    fake = MagicMock()
    fake._observed_files_session = set()
    fake.stream_message = MagicMock(side_effect=lambda *a, **k: iter([
        StreamEvent(type="thinking_delta", text="weighing the options…"),
        StreamEvent(type="thinking_delta", text="still weighing…"),
        # No usage reported, which is the endpoint this matters on.
        StreamEvent(type="message_delta", output_tokens=0, cost_usd=0.0),
    ]))
    with patch("delfin.agent.engine.create_client", return_value=fake):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick",
                             pack_dir=tmp_path)
    engine.stream_response("frag was")

    entries = _recorded(tm)
    assert entries, "the turn was not logged at all"
    assert entries[-1].get("ttft_ms") is not None, entries[-1]
    assert not tm.is_never_started(entries[-1]), entries[-1]
