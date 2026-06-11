"""Focused tests for turn outcome tracking in the dashboard agent tab.

Contract under test (current behavior of ``_record_turn_outcome``):
- every conversational mode is tracked (solo, dashboard, quick, ...)
  so the Activity tab shows real PASS/FAIL/PARTIAL across modes;
- empty turns are never recorded;
- verdict heuristics: clean turn -> PASS, denied commands without
  surrender -> PARTIAL, user stop with minimal output -> FAIL.
"""

from unittest.mock import MagicMock


def _engine(mode="solo"):
    engine = MagicMock()
    engine.mode = mode
    engine._stop_requested = False  # real engines expose a plain bool
    return engine


def test_record_turn_outcome_persists_clean_pass():
    from delfin.dashboard.tab_agent import _record_solo_turn_outcome

    engine = _engine()
    engine.record_cycle_outcome.return_value = {"success_rate": "solo: 93% → 94%"}

    changes = _record_solo_turn_outcome(
        engine,
        user_task="Fix prompt_loader.py",
        response_text="Patched and tested.",
        state={},
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


def test_denied_commands_without_surrender_are_partial():
    from delfin.dashboard.tab_agent import _record_solo_turn_outcome

    engine = _engine()
    engine.record_cycle_outcome.return_value = {}
    _record_solo_turn_outcome(
        engine,
        user_task="Fix prompt_loader.py",
        response_text="Patched and tested.",
        state={"_denied_commands": ["Bash(git commit)"]},
        start_time=12.5,
    )
    kwargs = engine.record_cycle_outcome.call_args.kwargs
    assert kwargs["verdict"] == "PARTIAL"
    assert kwargs["denied_commands"] == ["Bash(git commit)"]


def test_record_turn_outcome_tracks_every_mode_but_skips_empty():
    from delfin.dashboard.tab_agent import _record_solo_turn_outcome

    # Dashboard mode IS tracked — the Activity tab covers all modes.
    engine = _engine(mode="dashboard")
    engine.record_cycle_outcome.return_value = {"success_rate": "dashboard: 90%"}
    changes = _record_solo_turn_outcome(
        engine,
        user_task="Inspect dashboard state",
        response_text="Done.",
        state={},
        start_time=1.0,
    )
    assert changes == {"success_rate": "dashboard: 90%"}
    engine.record_cycle_outcome.assert_called_once()

    # Empty turns never produce an outcome row.
    engine2 = _engine()
    assert _record_solo_turn_outcome(
        engine2,
        user_task="Inspect dashboard state",
        response_text="",
        state={},
        start_time=1.0,
    ) == {}
    engine2.record_cycle_outcome.assert_not_called()


def test_stop_with_minimal_output_records_fail():
    from delfin.dashboard.tab_agent import _record_solo_turn_outcome

    engine = _engine()
    engine._stop_requested = True
    engine.record_cycle_outcome.return_value = {}
    _record_solo_turn_outcome(
        engine,
        user_task="long refactor task",
        response_text="ok",  # <80 chars after a user stop
        state={},
        start_time=0.0,
    )
    assert engine.record_cycle_outcome.call_args.kwargs["verdict"] == "FAIL"
