"""The number that ended a dashboard turn had no field, and was the lower one.

Two limits bound a turn. ``agent.max_tool_rounds`` caps the tool-call
rounds and has had a field in the settings tab for a while.
``agent.max_action_rounds`` caps the CONTINUATION loop that carries out
the commands the agent proposes -- and it defaulted to 12 while the other
was commonly set to 500, so it was the limit that bit, every time, and it
had no field at all.

A field report calls the per-agent step budget too low and says the value
is not reliably saved. Both halves follow from that: the saved value was
the wrong number, and the right number was unreachable. Every settings
writer in the product merges rather than replaces the ``agent`` section,
so nothing was dropping it -- it was simply never written.

The ceiling is also not what stops a degenerate loop:
``action_repeat_limit`` ends the turn on the first command re-issued in
the same turn, so the ceiling only ever binds work where every round
brings new commands. Raising it from 12 to 40 therefore costs no loop
protection.
"""

from __future__ import annotations

import inspect

from delfin.dashboard import tab_agent as T
from delfin.dashboard import tab_settings as TS
from delfin.user_settings import DEFAULT_SETTINGS


def _agent_defaults() -> dict:
    return (DEFAULT_SETTINGS.get("agent") or {}) if isinstance(
        DEFAULT_SETTINGS, dict) else {}


# ---------------------------------------------------------------------------
# The two numbers agree between the module and the settings file
# ---------------------------------------------------------------------------

def test_the_ceiling_default_is_one_number_in_two_places():
    """The module constant and the shipped setting are the same answer;
    two copies that drift would make the field show a value the loop does
    not use."""
    assert _agent_defaults()["max_action_rounds"] == (
        T._ACTION_ROUND_CEILING_DEFAULT)


def test_the_ceiling_leaves_room_for_honest_multi_step_work():
    """Not a magic number: the repeat guard is the loop protection, so this
    one only bounds rounds that each bring new commands."""
    assert T._ACTION_ROUND_CEILING_DEFAULT >= 40


def test_the_repeat_guard_is_what_stops_a_loop_and_is_unchanged():
    assert T._ACTION_ROUND_REPEAT_LIMIT_DEFAULT == 2
    assert _agent_defaults()["action_repeat_limit"] == 2


def test_a_resolved_limit_follows_the_setting(monkeypatch):
    import delfin.user_settings as US

    monkeypatch.setattr(
        US, "load_settings",
        lambda *a, **k: {"agent": {"max_action_rounds": 7,
                                   "action_repeat_limit": 3}})
    ceiling, repeat = T._resolve_action_round_limits()
    assert (ceiling, repeat) == (7, 3)


def test_zero_still_means_no_ceiling(monkeypatch):
    """A documented power-user choice: the repeat guard remains the stop."""
    import delfin.user_settings as US

    monkeypatch.setattr(US, "load_settings",
                        lambda *a, **k: {"agent": {"max_action_rounds": 0}})
    ceiling, _ = T._resolve_action_round_limits()
    assert ceiling > 1000


def test_a_corrupt_value_falls_back_to_the_default(monkeypatch):
    import delfin.user_settings as US

    monkeypatch.setattr(
        US, "load_settings",
        lambda *a, **k: {"agent": {"max_action_rounds": "many"}})
    ceiling, _ = T._resolve_action_round_limits()
    assert ceiling == T._ACTION_ROUND_CEILING_DEFAULT


# ---------------------------------------------------------------------------
# The number the stop note names can be reached
# ---------------------------------------------------------------------------

def test_the_stop_note_names_the_setting():
    note = T._format_action_stop_note(
        T._ACTION_ROUND_CEILING, ["git status"], 40)
    assert "max_action_rounds" in note
    assert "40" in note


def test_the_settings_tab_has_a_field_for_it():
    src = inspect.getsource(TS)
    assert "max_action_rounds_input" in src, (
        "the limit that ends a turn first must be reachable from the UI")


def test_the_field_is_loaded_saved_and_laid_out():
    src = inspect.getsource(TS)
    # loaded back, or the field shows a default over a value already chosen
    assert "max_action_rounds_input.value = int(" in src
    # written by the button next to it
    assert "['agent']['max_action_rounds'] = int(" in src
    # and reachable on the page
    assert "widgets.HBox([max_action_rounds_input" in src


def test_both_limits_are_written_by_the_same_button():
    """They sit in one box; a save that persisted only one of them is how
    the report's "not reliably saved" would survive this change."""
    src = inspect.getsource(TS)
    i = src.index("def _on_agentopt_save")
    body = src[i:src.index("agentopt_save_btn.on_click", i)]
    assert "max_tool_rounds" in body
    assert "max_action_rounds" in body


def test_the_two_fields_are_told_apart_in_the_text():
    """They are both "rounds per turn"; the one that stops the turn first
    has to say so or the reader raises the wrong one."""
    src = inspect.getsource(TS)
    i = src.index("max_action_rounds_input, max_action_rounds_hint")
    following = src[i:i + 1200]
    assert "max_action_rounds" in following
    assert "first" in following.lower()
