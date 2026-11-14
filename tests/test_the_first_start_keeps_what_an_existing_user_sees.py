"""The first start decides once: five tabs for a new installation, the old
view for anyone who used DELFIN before.

Input: the person's settings (``ui.tabs``, ``ui.first_start``) and the
traces of earlier use -- a saved agent session, a job in the local job
list, an entry in the calculations or archive folder.

Output: tab_defaults.first_start writes, once, a saved arrangement and the
marker ``ui.first_start``:

* no arrangement, no trace -> a new installation: the five FIRST_TABS
  visible, every other tab saved as hidden (switched on in Settings);
* no arrangement, any trace -> an existing user: the tab list and the
  visibility they saw before 2026-11-13 (legacy defaults), so nothing
  changes for them;
* an arrangement already saved -> left alone, only the marker.

From then on the saved arrangement decides; a marker already there means
nothing is written again. Rejected as a trace: the settings file itself
(the installer writes it before the dashboard ever starts).
"""
from __future__ import annotations

import inspect
import json

import pytest

from delfin import user_settings
from delfin.dashboard import tab_defaults

IDS = ["submit", "recalc", "orca_builder", "agent", "job_status", "calculations",
       "archive", "home", "office", "literature", "tools", "settings"]
#: The view before 2026-11-13: everything on except Home.
OLD_HIDDEN = {"home"}


def _specs():
    specs = [{"id": tab_id, "default_order": 10 * (n + 1),
              "default_visible": tab_id not in OLD_HIDDEN,
              "fixed": tab_id == "settings"} for n, tab_id in enumerate(IDS)]
    return tab_defaults.apply_first_tabs(specs)


@pytest.fixture
def place(tmp_path):
    paths = {"settings_path": tmp_path / ".delfin_settings.json",
             "calc_dir": tmp_path / "calc", "archive_dir": tmp_path / "archive",
             "sessions_dir": tmp_path / ".delfin" / "agent_sessions",
             "jobs_file": tmp_path / ".delfin_jobs.json"}
    paths["calc_dir"].mkdir()
    paths["archive_dir"].mkdir()
    return paths


def _first_start(place, specs=None):
    return tab_defaults.first_start(specs or _specs(), **place)


def _saved(place):
    return user_settings.load_settings(place["settings_path"])["ui"]


def _shown(place):
    ui = _saved(place)
    order, hidden = tab_defaults.resolve_tab_preferences(_specs(), ui["tabs"])
    return [t for t in order if t not in hidden]


def test_a_new_installation_saves_the_five(place):
    assert _first_start(place) == "new"
    assert _shown(place) == list(tab_defaults.FIRST_TABS)
    ui = _saved(place)
    assert ui["first_start"]["tabs"] == "new"
    assert set(ui["tabs"]["hidden"]) == set(IDS) - set(tab_defaults.FIRST_TABS) - {"settings"}


def test_a_settings_file_alone_is_no_trace_of_use(place):
    # The installer writes the settings file before the first start.
    user_settings.save_settings({"runtime": {"orca_base": "/opt/orca"}},
                                place["settings_path"])
    assert _first_start(place) == "new"


@pytest.mark.parametrize("trace", ["session", "job", "calc", "archive"])
def test_an_existing_user_keeps_the_old_view(place, trace):
    if trace == "session":
        place["sessions_dir"].mkdir(parents=True)
        (place["sessions_dir"] / "abc.json").write_text("{}")
    elif trace == "job":
        place["jobs_file"].write_text(json.dumps({"next_job_id": 1002,
                                                  "jobs": [{"job_id": 1001}]}))
    else:
        (place[f"{trace}_dir"] / "Fe_complex").mkdir()
    assert _first_start(place) == "existing"
    assert _shown(place) == [t for t in IDS if t not in OLD_HIDDEN and t != "settings"]
    assert _saved(place)["first_start"]["tabs"] == "existing"
    assert trace in _saved(place)["first_start"]["traces"][0]


def test_an_empty_job_list_or_hidden_files_are_no_trace(place):
    place["jobs_file"].write_text(json.dumps({"next_job_id": 1001, "jobs": []}))
    (place["calc_dir"] / ".DS_Store").write_text("")
    place["sessions_dir"].mkdir(parents=True)
    assert _first_start(place) == "new"


def test_a_saved_arrangement_is_left_alone(place):
    user_settings.save_settings({"ui": {"tabs": {"order": ["submit", "recalc"],
                                                  "hidden": ["recalc"]}}},
                                place["settings_path"])
    assert _first_start(place) == "arranged"
    ui = _saved(place)
    assert ui["tabs"] == {"order": ["submit", "recalc"], "hidden": ["recalc"]}
    assert ui["first_start"]["tabs"] == "arranged"


def test_the_marker_makes_it_happen_once(place):
    assert _first_start(place) == "new"
    (place["calc_dir"] / "later_job").mkdir()       # used since: still decided
    assert _first_start(place) is None
    assert _shown(place) == list(tab_defaults.FIRST_TABS)


def test_the_dashboard_decides_before_settings_reads_the_arrangement():
    from delfin import dashboard
    source = inspect.getsource(dashboard)
    seed = source.index("first_start(tab_specs")
    assert source.index("apply_first_tabs(tab_specs)") < seed
    assert seed < source.index("tab_settings.create_tab(ctx")
