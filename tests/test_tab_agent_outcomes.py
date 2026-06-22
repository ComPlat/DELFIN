"""Focused tests for conversational-turn outcome tracking in the agent tab.

Outcome tracking used to fire only in solo mode. Today it tracks every
conversational mode (solo / dashboard / quick / …) so the Activity tab
shows real PASS/FAIL/PARTIAL across the user's actual usage. The legacy
``_record_solo_turn_outcome`` symbol is preserved as an alias.
"""

from unittest.mock import MagicMock


def test_record_turn_outcome_persists_pass_when_clean_completion():
    """No denied commands, normal response → PASS verdict, persisted."""
    from delfin.dashboard.tab_agent import _record_turn_outcome

    engine = MagicMock()
    engine.mode = "solo"
    engine._stop_requested = False
    engine.record_cycle_outcome.return_value = {
        "success_rate": "solo: 93% → 94%",
    }

    changes = _record_turn_outcome(
        engine,
        user_task="Fix prompt_loader.py",
        response_text="Patched and tested.",
        state={"_denied_commands": []},
        start_time=12.5,
    )

    assert changes == {"success_rate": "solo: 93% → 94%"}
    engine.record_cycle_outcome.assert_called_once_with(
        verdict="PASS",
        user_task="Fix prompt_loader.py",
        denied_commands=[],
        error_type=None,
        start_time=12.5,
    )


def test_record_turn_outcome_partial_when_denied_commands():
    """Denied tool calls without surrender → PARTIAL, still persisted."""
    from delfin.dashboard.tab_agent import _record_turn_outcome

    engine = MagicMock()
    engine.mode = "solo"
    engine._stop_requested = False
    engine.record_cycle_outcome.return_value = {
        "success_rate": "solo: 93% → 92%",
    }

    _record_turn_outcome(
        engine,
        user_task="Fix prompt_loader.py",
        response_text="Patched and tested.",
        state={"_denied_commands": ["Bash(git commit)"]},
        start_time=12.5,
    )

    engine.record_cycle_outcome.assert_called_once_with(
        verdict="PARTIAL",
        user_task="Fix prompt_loader.py",
        denied_commands=["Bash(git commit)"],
        error_type=None,
        start_time=12.5,
    )


def test_record_turn_outcome_tracks_dashboard_mode_too():
    """Dashboard mode is no longer skipped — every mode is tracked."""
    from delfin.dashboard.tab_agent import _record_turn_outcome

    engine = MagicMock()
    engine.mode = "dashboard"
    engine._stop_requested = False
    engine.record_cycle_outcome.return_value = {"ok": True}

    _record_turn_outcome(
        engine,
        user_task="Inspect dashboard state",
        response_text="Done.",
        state={},
        start_time=1.0,
    )
    engine.record_cycle_outcome.assert_called_once()


def test_record_turn_outcome_skips_empty_task_or_response():
    """Either empty task or empty response → skip, no record_cycle_outcome."""
    from delfin.dashboard.tab_agent import _record_turn_outcome

    engine = MagicMock()
    engine.mode = "solo"

    assert _record_turn_outcome(
        engine,
        user_task="",
        response_text="Done.",
        state={},
        start_time=1.0,
    ) == {}
    engine.record_cycle_outcome.assert_not_called()

    assert _record_turn_outcome(
        engine,
        user_task="Inspect dashboard state",
        response_text="",
        state={},
        start_time=1.0,
    ) == {}
    engine.record_cycle_outcome.assert_not_called()


def test_record_solo_turn_outcome_alias_still_works():
    """Backwards-compat alias points at the new function."""
    from delfin.dashboard.tab_agent import (
        _record_solo_turn_outcome,
        _record_turn_outcome,
    )
    assert _record_solo_turn_outcome is _record_turn_outcome
