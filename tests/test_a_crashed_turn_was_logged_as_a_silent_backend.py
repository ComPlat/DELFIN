"""The turn log could not tell a silent backend from a turn that crashed.

``turn_metrics.record`` has taken an ``error`` argument since it was
written. The engine's only call site never passed one -- while holding the
exception three lines above it, in an ``except Exception as _turn_exc``
that even classifies it (``context_overflow`` / ``transient_api``) for the
provider profiles. The field existed, the value existed, and they never
met.

What that cost: ``is_never_started`` -- "the backend produced nothing for
20+ seconds" -- keys on ``ttft_ms is None`` plus a silent rest of the
record. A turn that raised before any event (a connect or read timeout, an
auth rejection, a 400) writes exactly that shape, so the log called it a
backend stall. After the fact the log is the only evidence, and it could
not say which of the two had happened.

Two more discriminators were missing:

* Time-to-first-TOKEN is stamped only by ``text_delta``. A stream that
  delivered ``message_start`` and then went quiet, and a stream that never
  delivered a byte, both record ``ttft_ms=None`` -- one is a live
  connection to a silent model, the other is a dead transport, and nothing
  in the record separated them. ``first_event_ms`` is stamped by ANY
  event, so ``None`` now means the transport delivered nothing at all.
* The client retries a transient failure three times (1.5 + 3 + 6 s of
  sleep plus request latency). Nothing in the record said so, and a turn
  that fought for twenty seconds looked exactly like a quiet one.

Provenance is the rule the module already followed and these fields keep:
a MISSING key is a record written before the field existed and must never
be judged; present-and-``None`` is the marker of an instrumented turn.
"""

from __future__ import annotations

import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import turn_metrics as tm
from delfin.agent.api_client import StreamEvent


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(tm, "_DIR", tmp_path / ".delfin" / "turn_metrics")
    return tmp_path


def _silent(**kw) -> dict:
    """The record shape a turn writes when no token ever arrived."""
    base = {"total_ms": 90_000, "ttft_ms": None, "first_event_ms": None,
            "output_chars": 0, "tool_calls": 0, "output_tokens": 0,
            "stopped": False, "error": "", "retries": None}
    base.update(kw)
    return base


# ---------------------------------------------------------------------------
# 1. A recorded error is what tells the two apart
# ---------------------------------------------------------------------------

def test_a_silent_backend_is_still_the_stall_it_always_was():
    assert tm.is_never_started(_silent()) is True
    assert tm.is_crashed(_silent()) is False


def test_a_turn_that_raised_is_a_crash_and_not_a_silent_backend():
    crashed = _silent(error="transient_api: Connection reset by peer")
    assert tm.is_crashed(crashed) is True
    assert tm.is_never_started(crashed) is False


def test_a_crash_no_longer_counts_towards_the_backend_stall_rate():
    assert tm.is_stall(_silent(error="APIConnectionError: dns failure")) is False
    assert tm.is_stall(_silent()) is True


def test_a_record_written_before_the_error_field_is_not_judged_as_crashed():
    legacy = {"total_ms": 90_000, "ttft_ms": None, "output_chars": 0,
              "tool_calls": 0}
    assert "error" not in legacy
    assert tm.is_crashed(legacy) is False


def test_the_roll_up_counts_crashes_instead_of_losing_them(home):
    tm.record("s", model="m", total_ms=90_000, ttft_ms=None, output_chars=0)
    tm.record("s", model="m", total_ms=90_000, ttft_ms=None, output_chars=0,
              error="transient_api: read timeout")
    s = tm.aggregate_turn_stats()
    assert s["turns"] == 2
    assert s["stalls"] == 1          # only the genuinely silent one
    assert s["crashes"] == 1         # the other is visible, not gone


def test_the_summary_line_names_the_crash_rather_than_blaming_the_backend(home):
    tm.record("s", model="m", total_ms=90_000, ttft_ms=None, output_chars=0,
              error="transient_api: read timeout")
    out = tm.format_summary(tm.read("s"))
    assert "crashed" in out
    assert "transient_api" in out
    assert "backend-stall" not in out


# ---------------------------------------------------------------------------
# 2. First EVENT, not first token: was the transport alive at all?
# ---------------------------------------------------------------------------

def test_a_dead_transport_is_told_apart_from_a_model_that_stayed_silent():
    dead = _silent(first_event_ms=None)
    alive = _silent(first_event_ms=40)
    assert tm.silence_kind(dead) == "transport"
    assert tm.silence_kind(alive) == "model"


def test_a_record_from_before_the_first_event_field_is_never_judged():
    legacy = {"total_ms": 90_000, "ttft_ms": None, "output_chars": 0,
              "tool_calls": 0, "output_tokens": 0, "error": ""}
    assert "first_event_ms" not in legacy
    assert tm.is_never_started(legacy) is True     # unchanged verdict
    assert tm.silence_kind(legacy) == "unclassified"


def test_an_answered_turn_has_no_silence_to_classify():
    answered = _silent(ttft_ms=800, first_event_ms=40, output_chars=20,
                       total_ms=2_000)
    assert tm.silence_kind(answered) == ""


def test_the_summary_says_which_kind_of_silence_it_was(home):
    tm.record("s", model="m", total_ms=90_000, ttft_ms=None,
              first_event_ms=None, output_chars=0)
    assert "no transport" in tm.format_summary(tm.read("s"))
    tm.record("s2", model="m", total_ms=90_000, ttft_ms=None,
              first_event_ms=40, output_chars=0)
    assert "model silent" in tm.format_summary(tm.read("s2"))


def test_the_field_is_written_as_none_rather_than_omitted(home):
    """Present-and-None is the marker of an instrumented turn."""
    tm.record("s", model="m", total_ms=10)
    entry = tm.read("s")[0]
    assert "first_event_ms" in entry and entry["first_event_ms"] is None
    assert "retries" in entry and entry["retries"] is None


# ---------------------------------------------------------------------------
# 3. Retries: a turn that fought is not a quiet turn
# ---------------------------------------------------------------------------

def test_a_retried_turn_no_longer_reads_like_a_quiet_one(home):
    tm.record("s", model="m", total_ms=90_000, ttft_ms=None, output_chars=0,
              retries=3, error="transient_api: read timeout")
    entry = tm.read("s")[0]
    assert entry["retries"] == 3
    assert "retries=3" in tm.format_summary([entry])


def test_the_recorded_error_does_not_carry_a_credential_to_the_disk(home):
    """The first free-form field in this log; a bug report bundles it."""
    secret = "sk-ant-api03-THISLOOKSLIKEAKEY0123456789abcdefghijklmnop"
    tm.record("s", model="m", total_ms=100,
              error=f"AuthenticationError: api_key={secret} rejected")
    entry = tm.read("s")[0]
    assert secret not in entry["error"]
    assert "AuthenticationError" in entry["error"]
    assert tm.is_crashed(entry) is True


def test_a_backend_that_reports_no_retries_says_none_not_zero(home):
    """``None`` is "this backend does not report retries", 0 is a count."""
    tm.record("s", model="m", total_ms=100, ttft_ms=50)
    tm.record("s", model="m", total_ms=100, ttft_ms=50, retries=0)
    a, b = tm.read("s")
    assert a["retries"] is None
    assert b["retries"] == 0


# ---------------------------------------------------------------------------
# The engine end: it had all three and passed none
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def _engine(agent_tree, stream, **client_attrs):
    from delfin.agent.engine import AgentEngine
    client = MagicMock()
    client.model = "kit.qwen3.5-397b-A17b"
    for k, v in client_attrs.items():
        setattr(client, k, v)
    client.stream_message = MagicMock(side_effect=stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api", provider="kit",
                          mode="quick", pack_dir=agent_tree)
    eng.client = client
    return eng


def _last(eng) -> dict:
    entries = tm.read(eng.trace_session())
    assert entries, "the turn should have been recorded"
    return entries[-1]


def test_the_engine_records_the_error_it_already_held(agent_tree, home):
    def boom(system, messages, **kw):
        raise RuntimeError("backend went away")
        yield                                   # pragma: no cover

    eng = _engine(agent_tree, boom)
    with pytest.raises(RuntimeError):
        eng.stream_response("hi")
    entry = _last(eng)
    assert "backend went away" in entry["error"]
    assert tm.is_crashed(entry) is True
    assert tm.is_never_started(entry) is False


class ReadTimeout(Exception):
    """A transport timeout, classified by class name like the real ones."""


def test_a_crashed_turn_reports_the_classification_the_engine_made(
        agent_tree, home):
    """The engine already classifies the exception for the profiles."""
    def boom(system, messages, **kw):
        raise ReadTimeout("timed out")
        yield                                   # pragma: no cover

    eng = _engine(agent_tree, boom)
    with pytest.raises(ReadTimeout):
        eng.stream_response("hi")
    assert _last(eng)["error"].startswith("transient_api")


def test_a_turn_that_never_reached_the_transport_stamps_no_first_event(
        agent_tree, home):
    def boom(system, messages, **kw):
        raise RuntimeError("connection refused")
        yield                                   # pragma: no cover

    eng = _engine(agent_tree, boom)
    with pytest.raises(RuntimeError):
        eng.stream_response("hi")
    assert _last(eng)["first_event_ms"] is None


def test_a_live_connection_to_a_silent_model_stamps_the_first_event(
        agent_tree, home):
    """The model thought, and then said nothing to the user.

    A thinking delta IS the model producing output — it is billed as
    output and it is the moment the wait ended. Treating it as no token
    put a turn that reasoned for three minutes and a turn that never
    heard from the endpoint into the same bucket, which is the whole
    distinction this file exists for.
    """
    def quiet(system, messages, **kw):
        yield StreamEvent(type="message_start", input_tokens=900)
        yield StreamEvent(type="thinking_delta", text="hmm")

    eng = _engine(agent_tree, quiet)
    eng.stream_response("hi")
    entry = _last(eng)
    assert entry["ttft_ms"] is not None         # thinking is production
    assert entry["first_event_ms"] is not None  # the stream was alive
    # ... and it still answered nothing, which is a different fault.
    assert tm.silence_kind(entry) == "model"


def test_a_backend_that_produced_nothing_at_all_stamps_no_ttft(
        agent_tree, home):
    """The other side of the same distinction: the envelope arrived and
    the model never began. That is the turn with no first token."""
    def quiet(system, messages, **kw):
        yield StreamEvent(type="message_start", input_tokens=900)

    eng = _engine(agent_tree, quiet)
    eng.stream_response("hi")
    entry = _last(eng)
    assert entry["ttft_ms"] is None
    assert entry["first_event_ms"] is not None
    assert tm.silence_kind(entry) == "model"


def test_the_engine_carries_the_retry_count_a_backend_reports(
        agent_tree, home):
    def quiet(system, messages, **kw):
        yield StreamEvent(type="message_start", input_tokens=1)

    eng = _engine(agent_tree, quiet, last_turn_retries=3)
    eng.stream_response("hi")
    assert _last(eng)["retries"] == 3


def test_a_backend_without_a_retry_counter_records_none(agent_tree, home):
    def quiet(system, messages, **kw):
        yield StreamEvent(type="message_start", input_tokens=1)

    eng = _engine(agent_tree, quiet)
    del eng.client.last_turn_retries             # a client that has no counter
    eng.stream_response("hi")
    assert _last(eng)["retries"] is None
