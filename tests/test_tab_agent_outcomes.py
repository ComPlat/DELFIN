"""Focused tests for solo outcome tracking in the dashboard agent tab."""

from unittest.mock import MagicMock


def test_record_solo_turn_outcome_persists_completed_turn():
    from delfin.dashboard.tab_agent import _record_solo_turn_outcome

    engine = MagicMock()
    engine.mode = "solo"
    engine.record_cycle_outcome.return_value = {"success_rate": "solo: 93% → 94%"}

    changes = _record_solo_turn_outcome(
        engine,
        user_task="Fix prompt_loader.py",
        response_text="Patched and tested.",
        state={"_denied_commands": ["Bash(git commit)"]},
        start_time=12.5,
    )

    assert changes == {"success_rate": "solo: 93% → 94%"}
    engine.record_cycle_outcome.assert_called_once_with(
        verdict="PASS",
        user_task="Fix prompt_loader.py",
        denied_commands=["Bash(git commit)"],
        start_time=12.5,
    )


def test_record_solo_turn_outcome_skips_non_solo_or_empty_turns():
    from delfin.dashboard.tab_agent import _record_solo_turn_outcome

    engine = MagicMock()
    engine.mode = "dashboard"

    assert _record_solo_turn_outcome(
        engine,
        user_task="Inspect dashboard state",
        response_text="Done.",
        state={},
        start_time=1.0,
    ) == {}
    engine.record_cycle_outcome.assert_not_called()

    engine.mode = "solo"
    assert _record_solo_turn_outcome(
        engine,
        user_task="Inspect dashboard state",
        response_text="",
        state={},
        start_time=1.0,
    ) == {}
    engine.record_cycle_outcome.assert_not_called()
