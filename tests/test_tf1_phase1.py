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


# ---------------------------------------------------------------------------
# Reviewer findings F1-F3 (tests/test_tf1_adversarial_tsfinder.py battery).
# Each pins a defect the reviewer's battery caught; all three were red on
# f75fb82d6 and are green on the fix.
# ---------------------------------------------------------------------------


def test_single_level_accepts_a_fully_matching_barrier():
    """F1: check_single_level_of_theory must return ok=True when every
    stage carries the SAME method line. On f75fb82d6 this raised TypeError
    ('list' object is not callable) from `key=list(levels)` in the sort,
    on every call -- an honest barrier never got a green verdict."""
    mod = importlib.import_module("delfin.ts_finder")
    levels = {
        "reactant": "PBE0-D3BJ/def2-SVP//:",
        "ts": "PBE0-D3BJ/def2-SVP//:",
        "product": "PBE0-D3BJ/def2-SVP//:",
        "zpe": "PBE0-D3BJ/def2-SVP//:",
    }
    result = mod.check_single_level_of_theory(levels)
    assert result.ok is True, f"matching barrier rejected: {result.problems}"
    assert result.problems == []


def test_single_level_reports_a_differing_stage():
    """F1 still rejects a barrier where one stage differs -- the check is
    meaningful, not vacuously green."""
    mod = importlib.import_module("delfin.ts_finder")
    levels = {
        "reactant": "PBE0-D3BJ/def2-SVP//:",
        "ts": "PBE0-D3BJ/def2-TZVP//:",  # different basis
        "product": "PBE0-D3BJ/def2-SVP//:",
        "zpe": "PBE0-D3BJ/def2-SVP//:",
    }
    result = mod.check_single_level_of_theory(levels)
    assert result.ok is False
    assert any("ts" in p for p in result.problems)


def test_single_level_shapes_come_back_as_a_check_result():
    """F1 hostile shape: a non-dict (including a list) must return a
    CheckResult (ok=False, reason), never raise. The reviewer's battery
    passed a hostile dict/list that previously crashed the check."""
    mod = importlib.import_module("delfin.ts_finder")
    for bad in ([("reactant", "PBE0/def2-SVP")], [], None, "PBE0/def2-SVP"):
        result = mod.check_single_level_of_theory(bad)
        assert result.ok is False
        assert isinstance(result.problems, list) and result.problems


def test_single_level_complains_about_a_missing_stage():
    """F1: a stage omitted from the barrier dict is reported, because a
    barrier must not be silent about an endpoint it did not run."""
    mod = importlib.import_module("delfin.ts_finder")
    levels = {"reactant": "PBE0-D3BJ/def2-SVP//:",
              "ts": "PBE0-D3BJ/def2-SVP//:",
              "product": "PBE0-D3BJ/def2-SVP//:"}  # zpe absent
    result = mod.check_single_level_of_theory(levels)
    assert result.ok is False
    assert any("zpe" in p for p in result.problems)


def test_irc_ends_are_minima_bad_shape_is_a_check_result():
    """F2: check_irc_ends_are_minima([1, 2]) raised AttributeError on
    f75fb82d6 ('list' object has no attribute 'items'). A wrong-shaped
    input must come back as CheckResult(False, reason), not a traceback --
    the pure-check layer must survive any one bad shape."""
    mod = importlib.import_module("delfin.ts_finder")
    result = mod.check_irc_ends_are_minima([1, 2])
    assert result.ok is False
    assert isinstance(result.problems, list) and result.problems


def test_irc_ends_are_minima_reports_a_saddle_end():
    """F2 keeps its job: an end that is not a minimum (>=1 imaginary mode)
    is reported, so the mechanism cannot be claimed ending at a saddle."""
    mod = importlib.import_module("delfin.ts_finder")
    result = mod.check_irc_ends_are_minima({"forward": 0, "reverse": 1})
    assert result.ok is False
    assert any("reverse" in p for p in result.problems)


def _write_synthetic_hess(tmp_path, bond_angstrom=1.11):
    """A minimal ORCA-style .hess for a two-atom H-C system whose only
    mode (-1500 cm-1) stretches the H-C bond by ~0.5 A at 0.25 A max shift
    (well over the 0.095 A response threshold). The $normal_modes vector is
    an asymmetric H(+x)/C(-x) stretch; $atoms coordinates in Bohr."""
    bohr = 0.529177
    xc_bohr = bond_angstrom / bohr
    text = (
        "$atoms\n"
        "2\n"
        f"H  1.008   0.00000000   0.00000000   0.00000000\n"
        f"C 12.011   {xc_bohr:.8f}   0.00000000   0.00000000\n"
        "$normal_modes\n"
        "6 1\n"
        "0\n"
        "0     1.00000000   0.00000000   0.00000000\n"
        "1     0.00000000   0.00000000   0.00000000\n"
        "2     0.00000000   0.00000000   0.00000000\n"
        "3    -1.00000000   0.00000000   0.00000000\n"
        "4     0.00000000   0.00000000   0.00000000\n"
        "5     0.00000000   0.00000000   0.00000000\n"
        "$vibrational_frequencies\n"
        "1\n"
        "0  -1500.000000\n"
    )
    path = tmp_path / "ts.hess"
    path.write_text(text, encoding="utf-8")
    return path


def test_imag_mode_that_acts_on_the_pair_says_verified(tmp_path):
    """F3: a mode that acts on all expected bond pairs must return ok=True
    (CheckResult(True, [])). On f75fb82d6 the success message was appended
    to problems and the result was ok=FALSE even for a correct mode -- the
    verifier could never say "found and verified". This is the reviewer's
    synthetic-Hessian case: mode 0 = -1500 cm-1 stretching the H-C pair."""
    mod = importlib.import_module("delfin.ts_finder")
    hess = _write_synthetic_hess(tmp_path)
    ts_xyz = (
        "2\n"
        "comment\n"
        "H     0.00000   0.00000   0.00000\n"
        "C     1.11000   0.00000   0.00000\n"
    )
    result = mod.check_imaginary_mode_connects_bonds(
        str(hess), ts_xyz, expected_pairs=[(0, 1)], mode_index=0)
    assert result.ok is True, f"acting mode rejected: {result.problems}"
    assert result.problems == []


def test_imag_mode_without_a_hessian_is_a_check_result(tmp_path):
    """F3: a missing Hessian must come back as CheckResult(False, reason),
    not crash -- the check is pure, bad input is a verdict not a traceback."""
    mod = importlib.import_module("delfin.ts_finder")
    ts_xyz = (
        "2\n"
        "comment\n"
        "H     0.00000   0.00000   0.00000\n"
        "C     1.11000   0.00000   0.00000\n"
    )
    result = mod.check_imaginary_mode_connects_bonds(
        str(tmp_path / "absent.hess"), ts_xyz,
        expected_pairs=[(0, 1)], mode_index=0)
    assert result.ok is False
    assert isinstance(result.problems, list) and result.problems

