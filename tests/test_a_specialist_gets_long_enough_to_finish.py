"""The caps were sized for `explore`, and the specialists cannot finish.

A delegated subagent had 40 tool calls and 900 s. That fits `explore`,
which uses about ten calls. It does not fit a specialist: `verifier`
exists to run the covering tests and then the same tests on the base as a
control, and the full suite takes 2520 s on this installation. A cap
below the thing a preset exists to do makes the preset useless and bills
the run for a report nobody can use.

120 calls and 3600 s. Caps, not targets -- the cost circuit-breaker and
the consecutive-failure abort are unchanged, and a delegate that finishes
in four calls still costs four.

The figures live in four places that have to move together: the
constants, the module docstring pinned to them, the shipped default, and
the limits block in the role prompt (pinned by
tests/test_subagent_fanout.py, because a stale budget in the prompt once
made the model decline to delegate at all). A fifth copy was in the
dashboard's help text, still reading 300 s from a default two changes
ago; it is generated now from the same function the tool description
uses.

And a file that still carried an OLD default follows the new one once,
because otherwise a raise is inert for every installation that has ever
opened the settings tab.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from delfin.agent import subagents as SA
from delfin.user_settings import DEFAULT_SETTINGS, load_settings

_SHIPPED = (DEFAULT_SETTINGS.get("agent") or {}).get("subagents") or {}
_PACK = Path(SA.__file__).resolve().parent / "pack" / "agents"


# ---------------------------------------------------------------------------
# One answer in every place that states it
# ---------------------------------------------------------------------------

def test_the_constants_and_the_shipped_default_agree():
    assert _SHIPPED["max_tool_calls"] == SA._MAX_TOOL_CALLS
    assert _SHIPPED["max_wall_s"] == int(SA._MAX_WALL_S)
    assert _SHIPPED["max_output_tokens"] == SA._MAX_OUTPUT_TOKENS


def test_the_module_docstring_states_them():
    doc = SA.__doc__ or ""
    assert f"max {SA._MAX_TOOL_CALLS} tool calls" in doc
    assert f"max {int(SA._MAX_WALL_S)} seconds" in doc


def test_the_role_prompt_states_them():
    text = (_PACK / "solo_agent.md").read_text(encoding="utf-8")
    i = text.index("**Backend limits per subagent run**")
    block = text[i:i + 260]
    assert f"{SA._MAX_TOOL_CALLS} tool calls" in block
    assert f"{int(SA._MAX_WALL_S)} s wall-clock" in block


def test_the_dashboard_help_is_generated_not_typed():
    """It read 300 s -- a default from two changes ago. The function that
    generates it is the one the tool description already uses, so there is
    one answer rather than a third copy to go stale."""
    import inspect

    from delfin.dashboard import tab_settings as TS
    src = inspect.getsource(TS)
    i = src.index("Caps for each delegated subagent")
    arm = src[i:i + 500]
    assert "_subagent_caps_now()" in arm
    assert not re.search(r"\b300\b", arm), "a figure is typed here again"
    assert not re.search(r"Defaults: \d", arm)


def test_the_generated_phrase_reports_what_is_in_force(monkeypatch):
    from delfin.agent import api_client as AC

    monkeypatch.setattr(SA, "_subagent_limits", lambda: {
        "max_tool_calls": 7, "max_wall_s": 11, "max_output_tokens": 13000})
    phrase = AC._subagent_caps_phrase()
    assert "7 tool calls" in phrase
    assert "11s" in phrase
    assert "13k" in phrase


# ---------------------------------------------------------------------------
# The numbers themselves
# ---------------------------------------------------------------------------

def test_the_wall_clock_fits_the_job_the_verifier_is_for():
    """It runs the suite and then a control run on the base. A cap under
    one suite run cannot be met by a preset whose method requires two."""
    assert SA._MAX_WALL_S >= 3600


def test_the_call_budget_is_not_explore_sized():
    assert SA._MAX_TOOL_CALLS >= 120


def test_the_other_stops_are_untouched():
    """The caps are a ceiling, not the loop protection: a run still ends
    on cost and on repeated failure."""
    assert SA._MAX_OUTPUT_TOKENS == 16000
    limits = SA._subagent_limits()
    assert set(limits) >= {"max_tool_calls", "max_wall_s",
                           "max_output_tokens"}


def test_a_caller_can_still_ask_for_less(monkeypatch):
    """A per-call argument beats the default in both directions; the raise
    must not become a floor."""
    monkeypatch.setattr(SA, "_subagent_limits", lambda: {
        "max_tool_calls": 2, "max_wall_s": 3, "max_output_tokens": 4000})
    got = SA._subagent_limits()
    assert got["max_wall_s"] == 3


# ---------------------------------------------------------------------------
# A file carrying an old default follows the new one, once
# ---------------------------------------------------------------------------

def _settings(tmp_path, subs) -> dict:
    p = tmp_path / "delfin_settings.json"
    p.write_text(json.dumps({"agent": {"subagents": subs}}), encoding="utf-8")
    return load_settings(p)["agent"]["subagents"]


@pytest.mark.parametrize("old_wall", [300, 900])
def test_an_old_shipped_default_follows_the_new_one(tmp_path, old_wall):
    """300 and 900 were the wall-clock default at different times; a file
    holding either is holding a default, not a decision."""
    got = _settings(tmp_path, {"max_wall_s": old_wall, "max_tool_calls": 40})
    assert got["max_wall_s"] == int(SA._MAX_WALL_S)
    assert got["max_tool_calls"] == SA._MAX_TOOL_CALLS


@pytest.mark.parametrize("subs", [
    {"max_wall_s": 1800, "max_tool_calls": 75},
    {"max_wall_s": 120, "max_tool_calls": 5},
    {"max_wall_s": 7200, "max_tool_calls": 200},
])
def test_a_number_somebody_chose_is_left_alone(tmp_path, subs):
    """The point of a setting is that it is respected. Overwriting one
    would make the dashboard lie about what is in force -- including a
    value LOWER than the default, which is a deliberate choice too."""
    assert _settings(tmp_path, dict(subs)) == pytest.approx(subs) or \
        _settings(tmp_path, dict(subs))["max_wall_s"] == subs["max_wall_s"]
    got = _settings(tmp_path, dict(subs))
    assert got["max_wall_s"] == subs["max_wall_s"]
    assert got["max_tool_calls"] == subs["max_tool_calls"]


def test_it_runs_once_and_a_later_choice_survives(tmp_path):
    p = tmp_path / "delfin_settings.json"
    p.write_text(json.dumps({"agent": {"subagents": {"max_wall_s": 900}}}),
                 encoding="utf-8")
    load_settings(p)                       # migrates
    data = json.loads(p.read_text(encoding="utf-8"))
    assert data["agent"]["subagents"]["max_wall_s"] == int(SA._MAX_WALL_S)
    # the user sets the old number back on purpose
    data["agent"]["subagents"]["max_wall_s"] = 900
    p.write_text(json.dumps(data), encoding="utf-8")
    assert load_settings(p)["agent"]["subagents"]["max_wall_s"] == 900


def test_the_migration_is_recorded_under_an_honest_name(tmp_path):
    """Not under the security updates: raising a cap is not a security
    fix, and the record must not say it was."""
    p = tmp_path / "delfin_settings.json"
    p.write_text(json.dumps({"agent": {"subagents": {"max_wall_s": 900}}}),
                 encoding="utf-8")
    load_settings(p)
    data = json.loads(p.read_text(encoding="utf-8"))
    assert data["default_updates_applied"]
    assert all("subagent" in k for k in data["default_updates_applied"])


def test_it_says_so_rather_than_changing_the_file_quietly(tmp_path):
    p = tmp_path / "delfin_settings.json"
    p.write_text(json.dumps({"agent": {"subagents": {"max_wall_s": 900}}}),
                 encoding="utf-8")
    load_settings(p)
    notices = json.loads(p.read_text(encoding="utf-8")).get("security_notices")
    assert notices
    joined = " ".join(notices)
    assert "900" in joined and str(int(SA._MAX_WALL_S)) in joined
    assert "will not be changed again" in joined


def test_an_absent_key_needs_no_migration(tmp_path):
    """The default already governs; touching the file would be noise."""
    got = _settings(tmp_path, {})
    assert got["max_wall_s"] == int(SA._MAX_WALL_S)


def test_a_corrupt_value_is_not_migrated(tmp_path):
    got = _settings(tmp_path, {"max_wall_s": "soon", "max_tool_calls": 40})
    assert got["max_wall_s"] == "soon"
    assert got["max_tool_calls"] == SA._MAX_TOOL_CALLS
