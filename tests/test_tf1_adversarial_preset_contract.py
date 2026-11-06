"""TF1 (ts-finder) review battery, part 2: preset loading contract.

The ts-finder expert is a department preset behind a switch. Before the
builder's ts-finder_subagent.md lands, this pins the REAL loading contract
that every *_subagent.md preset must satisfy, verified against current main
(04bf9703b, 2026-11-06):

  1. validate_preset() does NOT check 'department' or 'enabled_setting' --
     a preset whose switch is absent, off, or mis-pathed can still pass
     validation and be silently never offered. So a test on the TF1 spec
     must check BOTH validate_preset() AND the enabled_setting declaration.
     (Probe name is deliberately NOT 'ts-finder' now that the real pack
     preset of that name exists -- a shadow ts-finder would correctly trip
     the pack-override rule. Any other kebab-case name probes the same rule.)

  2. _setting_on(dotted_key) resolves an agent user-settings path (defaults
     off for a missing path); a department preset only enters the registry
     when its enabled_setting path is truthy in user settings
     (delfin/agent/subagents.py:1493-1509, called at 1690).
"""
from delfin.agent.subagents import validate_preset, _setting_on


def _preset(name, enabled=None):
    """Synthetic preset in the parser's format (frontmatter + body)."""
    lines = [
        "---",
        "name: " + name,
        "description: Finds and checks transition states.",
        "mode: default",
    ]
    if enabled is not None:
        lines.append("enabled_setting: " + enabled)
    lines += [
        "---",
        "# " + name,
        "Instructions for finding and verifying transition states.",
    ]
    return "\n".join(lines)


def test_validate_preset_accepts_a_plain_generic_preset():
    # A preset with a kebab-case name, description and body passes even
    # without a department / enabled_setting switch. Name is NOT 'ts-finder':
    # the real pack preset exists now, so a shadow name would collide.
    assert validate_preset(_preset("calc-verifier")) == []


def test_validate_preset_ignores_enabled_setting_and_department():
    # validate_preset is not the gate that guarantees the expert is offered:
    # a preset with NO switch passes validation (so a missing / mis-pathed
    # agent.departments.reactivity would slip through a validation-only test).
    assert validate_preset(_preset("calc-verifier")) == []
    assert validate_preset(_preset("calc-verifier", enabled="wrong.path")) == []


def test_setting_on_defaults_off_for_a_missing_path(monkeypatch):
    # The switch defaults OFF: an empty settings map never offers the
    # department. This is the shape the builder must satisfy: the frontmatter
    # must spell agent.departments.reactivity, and it must be truthy only when
    # the user has turned the department on.
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    assert _setting_on("agent.departments.reactivity") is False


def test_setting_on_is_true_when_the_department_switch_is_on(monkeypatch):
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda: {"agent": {"departments": {"reactivity": True}}},
    )
    assert _setting_on("agent.departments.reactivity") is True
