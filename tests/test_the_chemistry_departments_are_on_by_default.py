"""The two chemistry department switches are ON in the shipped defaults.

Input: ``user_settings.DEFAULT_SETTINGS`` and the preset files under
``delfin/agent/pack/agents``. Output: the preset registry composed by
``subagents._build_preset_registry()``. Semantics: a preset declaring
``enabled_setting: agent.departments.<name>`` enters the registry only
while that key is truthy, so the switch decides whether the model is
offered the type at all. Both keys now default True, which puts five
shipped specialists in the registry of a default installation; either key
set False still removes its own presets, so the switch keeps switching
off.

Every settings read here is pinned to a default installation
(``_merge_missing_defaults({}, DEFAULT_SETTINGS)``) instead of the host's
own settings file, so the file measures the product default on any
machine, including one whose personal file turns a department off.
"""

import delfin.agent.subagents as SA
import delfin.user_settings as us

REACTIVITY_KEY = "agent.departments.reactivity"
STRUCTURE_KEY = "agent.departments.structure"

REACTIVITY_PRESETS = ("ts-finder", "irc-checker", "reactivity-department")
STRUCTURE_PRESETS = ("metal-smiles-builder", "xyz-editor")


def _default_installation_settings():
    """The effective settings of a person with no overrides."""
    return us._merge_missing_defaults({}, us.DEFAULT_SETTINGS)


def _switch_off(cfg, name):
    """Turn one department off in ``cfg``, whatever the defaults hold.

    The off path is the never-worse half and must be measurable even
    where the key is absent, so this writes the key rather than
    assuming the shipped block is there.
    """
    cfg.setdefault("agent", {}).setdefault("departments", {})[name] = False


def _pin_settings(monkeypatch, cfg):
    """Make every settings read in this process return ``cfg``."""
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: cfg)


def test_both_department_keys_ship_in_the_defaults():
    departments = (us.DEFAULT_SETTINGS.get("agent") or {}).get("departments")
    assert isinstance(departments, dict), (
        "DEFAULT_SETTINGS['agent']['departments'] missing: a department "
        "switch that is absent reads as off")
    for key in ("reactivity", "structure"):
        assert key in departments, f"departments.{key} not in DEFAULT_SETTINGS"
        assert departments[key], f"departments.{key} ships off"


def test_the_switches_read_as_on_for_a_default_installation(monkeypatch):
    _pin_settings(monkeypatch, _default_installation_settings())
    assert SA._setting_on(REACTIVITY_KEY) is True
    assert SA._setting_on(STRUCTURE_KEY) is True


def test_the_five_department_presets_are_in_the_registry(monkeypatch):
    """The point of the default: the types exist, so the model can be
    offered them and a call does not hit the unknown-type refusal."""
    _pin_settings(monkeypatch, _default_installation_settings())
    registry = SA._build_preset_registry()
    missing = [name for name in REACTIVITY_PRESETS + STRUCTURE_PRESETS
               if name not in registry]
    assert missing == [], (
        f"not offered with the shipped defaults: {missing}; "
        f"registry has {sorted(registry)}")


def test_a_switch_turned_off_still_removes_its_presets(monkeypatch):
    """The never-worse half. A switch that can no longer be turned off is
    a regression, so the off path is measured too: reactivity off removes
    its three presets while structure, still on, keeps its two."""
    cfg = _default_installation_settings()
    _switch_off(cfg, "reactivity")
    _pin_settings(monkeypatch, cfg)
    assert SA._setting_on(REACTIVITY_KEY) is False
    registry = SA._build_preset_registry()
    for name in REACTIVITY_PRESETS:
        assert name not in registry, (
            f"{name} offered although its department switch is off")
    for name in STRUCTURE_PRESETS:
        assert name in registry, (
            f"{name} lost although its own switch is on: the off path "
            f"must remove one department, not the registry")


def test_the_other_switch_turned_off_removes_its_presets(monkeypatch):
    """The mirror case, so neither key is proven only through the other."""
    cfg = _default_installation_settings()
    _switch_off(cfg, "structure")
    _pin_settings(monkeypatch, cfg)
    assert SA._setting_on(STRUCTURE_KEY) is False
    registry = SA._build_preset_registry()
    for name in STRUCTURE_PRESETS:
        assert name not in registry, (
            f"{name} offered although its department switch is off")
    for name in REACTIVITY_PRESETS:
        assert name in registry, f"{name} lost although its switch is on"
