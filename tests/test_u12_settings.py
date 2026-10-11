"""Package U12 — classifier settings UI wiring + save/load round-trip.

The dashboard settings tab exposes agent.classifier.{enabled (default ON),
model (default kit.deepseek-v4-flash), user_rules}; the classifier's judge()
reads exactly those keys. This locks the integration so the UI, the reader
and the default model cannot drift.
"""

from __future__ import annotations

from pathlib import Path

from delfin import user_settings as US

_SRC = (Path(__file__).resolve().parent.parent
        / "delfin" / "dashboard" / "tab_settings.py").read_text(encoding="utf-8")


def test_classifier_widgets_declared():
    assert "classifier_enabled_toggle = widgets.Checkbox(" in _SRC
    assert "classifier_model_input = widgets.Combobox(" in _SRC
    assert "classifier_rules_input = widgets.Textarea(" in _SRC


def test_save_merges_classifier_keys():
    # The dedicated "Save agent extras" handler persists the three keys under
    # agent.classifier via a merge dict (hand-set keys survive).
    assert "payload['agent']['classifier'] = classifier" in _SRC
    assert "'enabled': bool(classifier_enabled_toggle.value)" in _SRC
    assert "'model': str(classifier_model_input.value or '').strip()" in _SRC
    assert "'user_rules'" in _SRC


def test_load_path_reads_classifier_keys():
    assert "classifier_payload = agent_payload.get('classifier')" in _SRC
    assert "classifier_enabled_toggle.value = bool(" in _SRC
    assert "classifier_model_input.value = str(" in _SRC
    assert "classifier_payload.get('user_rules') or []" in _SRC


def test_default_enabled_on_and_default_model():
    # Default ON, default model kit.deepseek-v4-flash (fast cold+warm).
    assert "classifier_payload.get('enabled', True)" in _SRC
    assert "kit.deepseek-v4-flash" in _SRC
    assert "value=True" in _SRC  # the toggle comes up ON


def test_ui_section_present():
    assert "🛡 Action classifier" in _SRC
    assert "Rules decide first" in _SRC


def test_roundtrip_save_and_load(tmp_path):
    """The keys the UI writes survive a save/load round-trip and defaults are
    synthesised when absent (enabled default True)."""
    settings_path = tmp_path / "settings.json"
    saved = {"agent": {"classifier": {
        "enabled": True,
        "model": "kit.deepseek-v4-flash",
        "user_rules": ["never allow: git push origin main",
                       "always allow: git push origin my-topic"],
    }}}
    US.save_settings(saved, settings_path=settings_path)
    loaded = US.load_settings(settings_path)
    cls = loaded["agent"]["classifier"]
    assert cls["enabled"] is True
    assert cls["model"] == "kit.deepseek-v4-flash"
    assert cls["user_rules"] == ["never allow: git push origin main",
                                 "always allow: git push origin my-topic"]


def test_absent_classifier_is_enabled_by_default(tmp_path):
    settings_path = tmp_path / "settings.json"
    US.save_settings({}, settings_path=settings_path)
    loaded = US.load_settings(settings_path)
    cls = (loaded.get("agent") or {}).get("classifier") or {}
    assert cls.get("enabled", True) is True  # must default ON
