"""Per-model knobs a user can set.

The numbers in ``model_profiles`` are measurements, and a measurement is a
claim about one deployment on one day. The KIT roster changed twice this
year and the same model name can be re-pointed at other hardware
underneath, so the person in front of the endpoint is sometimes the only
one who can see that a number has gone wrong. ``agent.model_overrides``
lets them fix it in the settings file.
"""

import json
from dataclasses import fields

import pytest

from delfin.agent import model_profiles as mp
from delfin.agent.model_profiles import ModelProfile


@pytest.fixture
def settings(tmp_path, monkeypatch):
    """Point the loader at a settings file this test owns."""
    path = tmp_path / "settings.json"

    def _write(overrides):
        path.write_text(json.dumps({"agent": {"model_overrides": overrides}}))
        mp._OVERRIDES_CACHE = None

    monkeypatch.setattr("delfin.user_settings.get_settings_path",
                        lambda *a, **k: path)
    mp._OVERRIDES_CACHE = None
    yield _write
    mp._OVERRIDES_CACHE = None


def test_an_exact_name_is_tuned(settings):
    settings({"kit.glm-5.3": {"stale_kill_after_s": 600,
                              "effort_default": "low"}})
    p = mp.get_profile("kit.glm-5.3")
    assert p.stale_kill_after_s == 600.0
    assert p.effort_default == "low"
    # Untouched knobs keep the measured value.
    assert p.max_tool_rounds == 20


def test_a_prefix_tunes_the_family(settings):
    settings({"kit.deepseek": {"max_tool_rounds": 30}})
    assert mp.get_profile("kit.deepseek-v4-flash").max_tool_rounds == 30
    assert mp.get_profile("kit.glm-5.3").max_tool_rounds == 20


def test_the_longer_key_wins(settings):
    settings({"kit.": {"max_tool_rounds": 5},
              "kit.glm": {"max_tool_rounds": 7}})
    assert mp.get_profile("kit.glm-5.3").max_tool_rounds == 7
    assert mp.get_profile("kit.deepseek-v4-flash").max_tool_rounds == 5


def test_an_exact_name_beats_a_prefix(settings):
    settings({"kit.glm": {"max_tool_rounds": 7},
              "kit.glm-5.3": {"max_tool_rounds": 9}})
    assert mp.get_profile("kit.glm-5.3").max_tool_rounds == 9


def test_a_model_with_no_registry_entry_can_still_be_tuned(settings):
    """The reason the key exists: a model DELFIN has never heard of."""
    settings({"kit.something-new": {"stale_kill_after_s": 900,
                                    "core_tools_only": True}})
    p = mp.get_profile("kit.something-new")
    assert p.stale_kill_after_s == 900.0
    assert p.core_tools_only is True


def test_an_override_says_so_in_the_notes(settings):
    """/model and the agent stats print notes. A number that differs from
    the file has to be visible, or the next person debugging it reads the
    source and believes it."""
    settings({"kit.glm-5.3": {"stale_kill_after_s": 600}})
    notes = mp.get_profile("kit.glm-5.3").notes
    assert "user override" in notes
    assert "stale_kill_after_s=600.0" in notes


def test_a_typo_does_not_take_the_agent_down(settings):
    settings({"kit.glm-5.3": {"stale_kil_after_s": 600,
                              "max_tool_rounds": "not a number",
                              "effort_default": "low"}})
    p = mp.get_profile("kit.glm-5.3")
    assert p.stale_kill_after_s == 420.0      # unknown key ignored
    assert p.max_tool_rounds == 20            # unusable value ignored
    assert p.effort_default == "low"          # the good one still applies


def test_broken_settings_are_not_fatal(tmp_path, monkeypatch):
    path = tmp_path / "settings.json"
    path.write_text("{not json")
    monkeypatch.setattr("delfin.user_settings.get_settings_path",
                        lambda *a, **k: path)
    mp._OVERRIDES_CACHE = None
    assert mp.get_profile("kit.glm-5.3").stale_kill_after_s == 420.0
    mp._OVERRIDES_CACHE = None


def test_no_settings_file_means_the_measured_values(tmp_path, monkeypatch):
    monkeypatch.setattr("delfin.user_settings.get_settings_path",
                        lambda *a, **k: tmp_path / "absent.json")
    mp._OVERRIDES_CACHE = None
    p = mp.get_profile("kit.glm-5.3")
    assert p.stale_kill_after_s == 420.0
    assert "user override" not in p.notes
    mp._OVERRIDES_CACHE = None


def test_reading_the_overrides_never_rewrites_the_settings_file(settings):
    """``load_settings`` fills in defaults and writes the file back. This
    runs while a prompt is being built; it must only read."""
    from delfin.user_settings import get_settings_path
    settings({"kit.glm-5.3": {"stale_kill_after_s": 600}})
    path = get_settings_path()
    before = path.read_bytes()
    mp.get_profile("kit.glm-5.3")
    assert path.read_bytes() == before


def test_every_knob_has_a_declared_type():
    """A new ModelProfile field must be a deliberate decision about whether
    a user may set it -- not silently un-settable because nobody updated
    the coercion table."""
    assert set(mp._COERCE) == {f.name for f in fields(ModelProfile)}
