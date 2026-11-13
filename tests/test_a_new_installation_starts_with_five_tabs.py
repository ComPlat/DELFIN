"""A new installation starts with five tabs; every other one is a switch away.

Input: the dashboard's tab specs and the person's ``ui.tabs`` settings
(``order``: the arrangement saved from Settings > Dashboard Tabs,
``hidden``: the tabs switched off there). Output:
tab_defaults.resolve_tab_preferences -> (order, hidden), which the
dashboard's tab bar and nothing else reads.

Decision (2026-11-13): the defaults -- what a person sees who never
arranged the tabs, a new installation first of all -- are Submit Job,
Job Status, DELFIN Agent, Home, Calculations (and Settings, which is
always there). Every other tab stays available and is switched on in
Settings. A person who arranged the tabs keeps exactly that arrangement:
a saved order lists every tab of its day, and its hidden list decides
those. A tab added in a later version follows the default. An existing
installation that never arranged its tabs gets the new default too:
nothing distinguishes it from a new one (``order`` is empty in both; the
installer itself writes the settings file).
"""
from __future__ import annotations

import inspect

from delfin.dashboard import tab_defaults

IDS = ["submit", "recalc", "orca_builder", "chemdarwin", "agent", "agent_activity",
       "job_status", "calculations", "archive", "home", "office", "literature",
       "remote_archive", "tools", "swarm", "reactions", "settings"]


def _specs():
    specs = [{"id": tab_id, "default_order": 10 * (n + 1), "default_visible": True,
              "fixed": tab_id == "settings"} for n, tab_id in enumerate(IDS)]
    return tab_defaults.apply_first_tabs(specs)


def _shown(prefs):
    order, hidden = tab_defaults.resolve_tab_preferences(_specs(), prefs)
    return [t for t in order if t not in hidden]


def test_a_person_who_never_arranged_the_tabs_sees_the_five():
    assert _shown({}) == ["submit", "job_status", "agent", "home", "calculations"]
    assert _shown({"order": [], "hidden": []}) == list(tab_defaults.FIRST_TABS)


def test_a_saved_arrangement_is_kept_whole():
    saved = {"order": ["calculations", "submit", "recalc", "archive", "swarm", "agent",
                       "job_status", "home", "tools"],
             "hidden": ["home"]}
    assert _shown(saved) == ["calculations", "submit", "recalc", "archive", "swarm",
                             "agent", "job_status", "tools"]


def test_a_tab_new_since_the_saved_arrangement_follows_the_default():
    saved = {"order": ["submit", "recalc"], "hidden": []}
    shown = _shown(saved)
    assert shown[:2] == ["submit", "recalc"]
    assert "literature" not in shown          # newer than the arrangement: default off
    assert "agent" in shown                   # newer, and one of the five


def test_settings_is_never_one_of_the_switchable_tabs():
    specs = _specs()
    settings = next(s for s in specs if s["id"] == "settings")
    assert settings["default_visible"] is True
    order, _hidden = tab_defaults.resolve_tab_preferences(specs, {})
    assert "settings" not in order


def test_the_dashboard_applies_the_defaults_to_every_tab_it_has():
    from delfin import dashboard
    source = inspect.getsource(dashboard)
    apply_at = source.index("apply_first_tabs(tab_specs)")
    assert source.index("tab_specs.extend(registered_tab_specs(ctx))") < apply_at
    assert apply_at < source.index("ctx.tab_specs = tab_specs")
    assert "resolve_tab_preferences(" in source
