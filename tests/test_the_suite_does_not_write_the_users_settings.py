"""The test suite does not write the user's settings file.

On 2026-09-15 a test run wrote ``~/.delfin_settings.json``: save_settings
merges every missing default into what it writes, so a newly added default
(agent.git_role = contributor) appeared in the maintainer's own settings,
where it would have refused their pushes to main. The shared isolation
fixture covered ``~/.delfin`` and nothing beside it.
"""

from __future__ import annotations

import pathlib

from delfin import user_settings as us


def _real_settings_file() -> pathlib.Path:
    return pathlib.Path.home() / us.SETTINGS_FILE_NAME


def test_saving_settings_in_a_test_leaves_the_real_file_alone():
    real = _real_settings_file()
    before = real.stat().st_mtime_ns if real.exists() else None

    written = us.get_settings_path()
    us.save_settings({"agent": {"effort": "low"}})

    assert written != real
    assert written.exists()
    after = real.stat().st_mtime_ns if real.exists() else None
    assert before == after


def test_an_explicit_path_is_still_used(tmp_path):
    target = tmp_path / "settings.json"
    us.save_settings({"agent": {"effort": "low"}}, target)
    assert us.load_settings(target)["agent"]["effort"] == "low"
