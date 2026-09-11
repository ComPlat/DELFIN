"""What the dashboard's controls say must be what the model gets.

Three of them did not. The Effort selector persisted the choice and
stopped: the engine kept the level it was built with, so a session that
started at "low" answered at "low" whatever the control said. The Model
selector moved during a running turn while the engine, and with it the
model, stayed put -- silently, and without even saving the choice. And
every switch that did rebuild the engine set the reference to None and
nothing else: the chat panel still showed every message, the new engine
had none of them, and the next question was answered without the
context on screen.

The wait is the fourth thing. In solo mode the first-token watchdog was
never armed, and the spinner said "Thinking..." for as long as a queued
endpoint liked; measured on 2026-09-11, kit.glm-5.3 sent nothing for a
one-word question with a one-line system prompt for over ten minutes,
while kit.deepseek-v4-flash answered the same in 0.6-4 s. The spinner
now counts the wait, names the queue past the profile's cold start, and
the budget applies in every mode.
"""

from __future__ import annotations

import pathlib
from unittest.mock import MagicMock, patch

import pytest

_TAB = (pathlib.Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_agent.py")


def _source() -> str:
    return _TAB.read_text(encoding="utf-8")


def _body(name: str, text: str, span: int = 3500) -> str:
    i = text.index(f"def {name}(")
    return text[i:i + span]


# ---------------------------------------------------------------------------
# Effort reaches the engine that exists
# ---------------------------------------------------------------------------

def _engine(tmp_path, client):
    from delfin.agent import engine as E
    with patch("delfin.agent.engine.create_client", return_value=client):
        return E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                             model="kit.glm-5.3", mode="solo", effort="low")


def test_the_engine_hands_a_new_level_to_a_client_that_reads_it_per_request(tmp_path):
    client = MagicMock()
    client.EFFORT_PER_REQUEST = True
    eng = _engine(tmp_path, client)
    assert eng.effort == "low"
    assert eng.set_effort("High") is True
    assert eng.effort == "high"
    assert client.effort == "high", "the client sends what it holds; it held the old level"


def test_a_backend_that_fixes_the_level_at_start_says_so(tmp_path):
    client = MagicMock()
    client.EFFORT_PER_REQUEST = False
    eng = _engine(tmp_path, client)
    assert eng.set_effort("high") is False, (
        "a CLI subprocess took --effort when it started; claiming the change "
        "applied would be the lie this file is about")
    assert eng.effort == "high"


def test_the_api_clients_read_the_level_per_request_and_the_cli_does_not():
    from delfin.agent.api_client import APIClient, CLIClient, OpenAIClient
    assert OpenAIClient.EFFORT_PER_REQUEST is True
    assert APIClient.EFFORT_PER_REQUEST is True
    assert CLIClient.EFFORT_PER_REQUEST is False


def test_the_openai_client_sends_the_level_it_holds_now():
    """The request builder reads ``self.effort`` at call time, which is
    what makes a live change take effect on the next round."""
    import inspect
    from delfin.agent import api_client
    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    assert '_eff_setting = getattr(self, "effort", "")' in src


def test_the_effort_control_hands_the_level_to_the_engine():
    body = _body("_on_effort_change", _source())
    assert "engine.set_effort(level)" in body
    assert 'if state["streaming"] or state.get("_controls_sync_internal")' not in body, (
        "a change during a running turn was dropped on the floor, not even saved")
    assert "_engine_rebuild_pending" in body, (
        "a backend that cannot take the level live must apply it after the turn")
    assert "_drop_engine()" in body


# ---------------------------------------------------------------------------
# A model switch during a turn is honoured, later, and said so
# ---------------------------------------------------------------------------

def test_a_model_switch_during_a_turn_is_queued_not_swallowed():
    body = _body("_on_model_change", _source())
    assert '        if state["streaming"]:\n            return' not in body
    assert '_engine_rebuild_pending"] = "model"' in body
    assert "when the running turn" in body


def test_the_turn_end_rebuilds_for_every_pending_switch():
    text = _source()
    i = text.index('_pending_other = state.pop("_engine_rebuild_pending", "")')
    window = text[i - 200:i + 200]
    assert '_pending_perm = state.pop("_perm_rebuild_pending", "")' in window
    assert "_drop_engine()" in window


def test_the_slash_commands_leave_the_announcement_to_the_observer():
    """Two messages for one switch read as two switches."""
    text = _source()
    i = text.index('if cmd.startswith("/model "):')
    assert 'Model switched to {name}.' not in text[i:i + 600]
    j = text.index('if cmd.startswith("/effort "):')
    assert 'Effort set to {level}.' not in text[j:j + 600]


# ---------------------------------------------------------------------------
# The conversation survives an engine rebuild
# ---------------------------------------------------------------------------

def test_the_dropped_engines_conversation_is_carried_to_the_next():
    text = _source()
    body = _body("_drop_engine", text, 2000)
    assert "engine.export_state()" in body
    assert 'state["_engine_carry_over"] = carry' in body
    assert "engine.client.kill()" in body
    ensure = _body("_ensure_engine", text, 12000)
    assert 'carry = state.pop("_engine_carry_over", None)' in ensure
    assert "engine.restore_state({**carry" in ensure


@pytest.mark.parametrize("handler", ["_on_model_change", "_on_provider_change",
                                     "_on_perm_change"])
def test_every_switch_that_rebuilds_goes_through_the_carrying_drop(handler):
    body = _body(handler, _source(), 5000)
    assert "_drop_engine()" in body
    assert 'state["engine"] = None' not in body, (
        f"{handler} still drops the engine bare, losing its conversation")


def test_export_then_restore_carries_the_messages_between_two_engines(tmp_path):
    """The mechanism the dashboard relies on, exercised on real engines."""
    client = MagicMock()
    old = _engine(tmp_path, client)
    old.messages.extend([
        {"role": "user", "content": "which basis set did the archive use?"},
        {"role": "assistant", "content": "def2-TZVP, in all three folders."},
    ])
    carry = dict(old.export_state())
    assert carry["engine_messages"]
    new = _engine(tmp_path, MagicMock())
    assert not new.messages
    new.restore_state({**carry, "mode": "solo"})
    assert [m["content"] for m in new.messages[-2:]] == [
        "which basis set did the archive use?",
        "def2-TZVP, in all three folders.",
    ]


# ---------------------------------------------------------------------------
# The wait is counted, and bounded in every mode
# ---------------------------------------------------------------------------

def test_the_label_counts_the_wait_and_names_the_queue():
    from delfin.dashboard.tab_agent import _waiting_label
    early = _waiting_label("kit.glm-5.3", 25, first=True, slow_cold_start_s=200)
    assert early == "Waiting for kit.glm-5.3 · 0:25 without a token"
    cold = _waiting_label("kit.glm-5.3", 187, first=True, slow_cold_start_s=200)
    assert "3:07 without a token" in cold and "usually takes ~200 s" in cold
    queued = _waiting_label("kit.glm-5.3", 640, first=True, slow_cold_start_s=200)
    assert "10:40" in queued and "queueing" in queued and "pick another model" in queued
    later = _waiting_label("kit.deepseek-v4-flash", 65, first=False)
    assert later == "Waiting for kit.deepseek-v4-flash · 1:05 since the last token"


def test_a_model_without_a_cold_start_figure_is_not_told_one():
    from delfin.dashboard.tab_agent import _waiting_label
    assert "usually" not in _waiting_label("x", 120, first=True, slow_cold_start_s=0)
    assert "queueing" in _waiting_label("x", 301, first=True, slow_cold_start_s=0)


def test_the_first_token_budget_is_armed_in_every_mode():
    body = _body("_arm_stale_watcher", _source(), 14000)
    assert "kill_after = 0  # disabled (solo mode)" in body
    assert "_due = [b for b in (kill_after, first_token_kill) if b > 0]" in body, (
        "solo mode armed no timer at all, so the first token had no budget")
    assert "if budget <= 0:\n                    return" in body, (
        "the mid-stream kill stays off where it was off")
    assert "3.0 * _slow_cold" in body
    assert "kill_after * 4.0" not in body


def test_the_first_token_budget_is_ten_minutes_for_glm_not_twenty_eight():
    from delfin.agent.model_profiles import get_profile
    p = get_profile("kit.glm-5.3")
    assert p.stale_kill_after_s == 420.0 and p.slow_cold_start_s == 200.0
    assert max(600.0, 3.0 * p.slow_cold_start_s) == 600.0
    assert max(600.0, p.stale_kill_after_s * 4.0) == 1680.0, "the old formula"


def test_the_spinner_ticks_while_the_provider_is_silent():
    body = _body("_arm_stale_watcher", _source(), 14000)
    assert "def _tick_wait():" in body
    assert "_waiting_label(" in body
    assert "_threading.Timer(10.0, _tick_wait)" in body
    text = _source()
    assert text.count("                    _mark_output()\n") == 3, "every output site stamps the first token"
    assert '_turn_first_output_monotonic"] = now' in text


def test_the_status_row_says_how_long_the_last_turn_took():
    from delfin.dashboard.tab_agent import _render_status, _turn_timing_text
    text = _turn_timing_text(303.4, 252.0, 3)
    assert text == "4:12 to first token · 5:03 turn · 3 tool calls"
    assert _turn_timing_text(12.0, -1.0, 0) == "0:12 turn"
    assert _turn_timing_text(70.0, 3.2, 1) == "0:03 to first token · 1:10 turn · 1 tool call"
    html = _render_status("solo", "api", "solo_agent", 0, 0, 100, 20, 0.01,
                          provider="kit", model="kit.glm-5.3",
                          last_turn_cost_usd=0.004, last_turn_timing=text)
    assert "4:12 to first token" in html and "5:03 turn" in html


def test_the_turn_end_writes_the_timing_the_status_row_reads():
    text = _source()
    assert 'state["_last_turn_timing"] = _turn_timing_text(' in text
    assert 'last_turn_timing=str(state.get("_last_turn_timing") or "")' in text


# ---------------------------------------------------------------------------
# A retry after ten silent minutes says what it is
# ---------------------------------------------------------------------------

def test_a_zero_byte_cut_is_named_and_a_dropped_connection_is_not_confused_with_it():
    from delfin.agent.api_client import _retry_notice, _ZERO_BYTE_CUT_S
    cut = _retry_notice("InternalServerError", 600.3, False, 1, 3, 2.0)
    assert "No byte from the endpoint in 600 s" in cut
    assert "cut the request before the model started" in cut
    assert "joins its queue again" in cut and "1/3" in cut
    dropped = _retry_notice("APIConnectionError", 600.3, True, 2, 3, 3.0)
    assert "connection lost mid-answer" in dropped and "2/3" in dropped
    quick = _retry_notice("RateLimitError", 4.0, False, 1, 3, 2.0)
    assert quick == "\n⏳ Transient API error (RateLimitError); retrying 1/3 in 2s…\n"
    assert _ZERO_BYTE_CUT_S == 300.0


def test_the_stream_loop_stamps_the_round_and_uses_the_notice():
    import inspect
    from delfin.agent import api_client
    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    assert "_round_t0 = time.monotonic()" in src
    assert "_retry_notice(" in src


def test_the_watchdog_does_not_advise_a_longer_budget_where_it_cannot_help():
    body = _body("_arm_stale_watcher", _source(), 16000)
    assert "The KIT gateway cuts a request the model has not " in body
    assert "started after 600 s, so a longer budget cannot " in body
    assert 'provider_dropdown.value or "") == "kit"' in body
    assert "raise " in body, "the other providers keep the setting advice"
