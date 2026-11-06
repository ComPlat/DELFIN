"""Phase-1 red controls for TF1.

These tests pin the package skeleton the task requires: the three packaged
presets of the reactivity department and the pure-check library module
``delfin.ts_finder``. They are written against the public packaging surface,
not against .gate copies.

Preset files and their bodies are absent before this package lands, and
``delfin.ts_finder`` does not import -- so all tests below are red on the
preceding commit, green once the package skeleton exists. The numeric red
control is recorded in the phase-1 commit message.
"""
from __future__ import annotations

import importlib
from pathlib import Path

from delfin.agent import subagents as SA


REPO = Path(__file__).resolve().parents[1]
AGENTS = REPO / "delfin" / "agent" / "pack" / "agents"

TS_FINDER = AGENTS / "ts-finder_subagent.md"
IRC_CHECKER = AGENTS / "irc-checker_subagent.md"
REACTIVITY_DEPT = AGENTS / "reactivity-department_subagent.md"

_PRESET_PATHS = (TS_FINDER, IRC_CHECKER, REACTIVITY_DEPT)


def _frontmatter(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8")
    parts = text.split("---", 2)
    assert len(parts) >= 3, f"{path.name} lacks YAML frontmatter"
    data: dict[str, str] = {}
    for raw in parts[1].splitlines():
        if ":" not in raw:
            continue
        key, value = raw.split(":", 1)
        data[key.strip()] = value.strip()
    return data


def test_tf1_packaged_preset_files_exist():
    missing = [p.name for p in _PRESET_PATHS if not p.is_file()]
    assert missing == [], f"missing packaged preset files: {missing}"


def test_tf1_preset_files_validate():
    """Each preset passes validate_preset except the one rule that cannot
    apply to a shipped pack file: its own name *is* an existing pack id, so
    the uniqueness rule flags it against itself. The framework's own oracle
    (test_sa0_validate_preset.py::test_every_existing_preset_passes_validation)
    filters exactly that self-collision; we do the same. Every other rule --
    kebab-case name, body present, tools known, no protected write_globs --
    must pass with zero problems.
    """
    for path in _PRESET_PATHS:
        assert path.is_file(), f"missing preset: {path.name}"
        problems = SA.validate_preset(str(path))
        # A shipped preset name IS an existing id; filter the self-collision
        # the way the framework's oracle does: the problem names the preset's
        # own name AND says override/exists. Anything else stays.
        own = path.stem.replace("_subagent", "")
        remaining = [p for p in problems
                     if not (own in p and ("override" in p or "exists" in p))]
        assert remaining == [], (
            f"{path.name}: validate_preset problems beyond the expected "
            f"self-collision: {problems!r}")


def test_tf1_preset_metadata_declares_the_reactivity_department():
    ts_meta = _frontmatter(TS_FINDER)
    irc_meta = _frontmatter(IRC_CHECKER)
    dept_meta = _frontmatter(REACTIVITY_DEPT)
    assert ts_meta["department"] == "reactivity"
    assert irc_meta["department"] == "reactivity"
    assert dept_meta["department"] == "reactivity"
    assert ts_meta["enabled_setting"] == "agent.departments.reactivity"
    assert irc_meta["enabled_setting"] == "agent.departments.reactivity"
    assert dept_meta["enabled_setting"] == "agent.departments.reactivity"


def test_tf1_department_is_off_when_the_switch_is_off(monkeypatch):
    """agent.departments.reactivity defaults to OFF: the department presets
    are not offered to the model. This is the framework's own behaviour for
    any preset with an enabled_setting (subagents._build_preset_registry) --
    the presets only register when the switch is on."""
    monkeypatch.setattr(SA, "_setting_on", lambda key: False)
    SA.reload_subagent_presets()
    try:
        for name in ("ts-finder", "irc-checker", "reactivity-department"):
            assert name not in SA.SUBAGENT_PRESETS, (
                f"{name} offered with the department switch OFF")
    finally:
        SA.reload_subagent_presets()


def test_tf1_department_registers_when_the_switch_is_on(monkeypatch):
    """With the switch on, the three presets (and only they) register; the
    model can then be offered a ts-finder / irc-checker / manager type."""
    monkeypatch.setattr(SA, "_setting_on", lambda key: True)
    SA.reload_subagent_presets()
    try:
        for name in ("ts-finder", "irc-checker", "reactivity-department"):
            assert name in SA.SUBAGENT_PRESETS, (
                f"{name} not registered with the switch ON")
    finally:
        SA.reload_subagent_presets()


def test_ts_finder_module_exports_phase1_checks():
    mod = importlib.import_module("delfin.ts_finder")
    required = {
        "check_endpoint_mapping",
        "check_single_level_of_theory",
        "check_imaginary_mode_connects_bonds",
        "check_irc_endpoint_match",
        "summarize_ts_verdict",
    }
    missing = sorted(name for name in required if not hasattr(mod, name))
    assert missing == [], f"missing public API: {missing}"


def test_ts_finder_module_uses_honest_verdict_labels():
    mod = importlib.import_module("delfin.ts_finder")
    assert mod.summarize_ts_verdict(found=False, mode_ok=False,
                                    irc_ok=False) == "not found"
    assert mod.summarize_ts_verdict(found=True, mode_ok=False,
                                    irc_ok=False) == "found but wrong mode"
    assert mod.summarize_ts_verdict(found=True, mode_ok=True,
                                    irc_ok=True) == "found and verified"
