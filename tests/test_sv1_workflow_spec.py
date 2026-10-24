"""SV1 phase 2b: the solvation workflow as a WB2 pipeline spec.

The agent-side path from a request to a checked, savable workflow:
resolve_spec -> validate_spec -> scientific_lint -> dry run -> cost
(delfin/agent/workflow_builder.run_checks, WB2's chain). This file
builds the SV1 solvation workflow AS A SPEC (the shape the user can
save with delfin.tools.user_workflows.save and PR as a contrib
template) and runs the FULL check chain over it.

Reviewer angle R2 (nacht-s19): the spec must carry the
level-of-theory separation (xTB structure preparation vs the SOLVATOR
level) and charge/spin inheritance — the checks and the lint must see
them.
"""

import pytest


def _solva_spec() -> dict:
    """The SV1 solvation workflow as a serialized pipeline spec.

    solvent_xyz and the solute geometry are user-supplied at run time —
    the spec carries them as `{placeholders}` (the WB2 template
    convention; _validate.py:112 treats them as resolve-at-runtime).
    """
    return {
        "delfin_pipeline_spec": "1",
        "name": "explicit_solvation",
        "steps": [
            {"step": "packmol_solvate",
             "label": "packmol_box",
             "kwargs": {"solvent_xyz": "{solvent_xyz}",
                        "n_solvent": 50, "box_side": 20.0,
                        "tolerance": 2.0}},
            {"step": "xtb_opt",
             "label": "xtb_relax",
             "kwargs": {"charge": 0, "mult": 1}},
            {"step": "xtb_solvator",
             "label": "explicit_shell",
             "kwargs": {"solvent": "water", "n_solvent": 12,
                        "method": "GFN2-xTB", "charge": 0, "mult": 1}},
        ],
    }


def test_resolve_fills_the_declared_defaults():
    from delfin.tools.platform import resolve_spec

    resolved = resolve_spec(_solva_spec())
    steps = {s["step"]: s for s in resolved["steps"]}
    # packmol defaults the adapter declares
    assert steps["packmol_solvate"]["resolved_kwargs"]["n_solvent"] == 50
    assert steps["packmol_solvate"]["resolved_kwargs"]["box_side"] == 20.0
    # the solvator's declared xTB-only method survives resolve
    assert steps["xtb_solvator"]["resolved_kwargs"]["method"] == "GFN2-xTB"


def test_the_full_check_chain_passes():
    from delfin.agent.workflow_builder import run_checks

    report = run_checks(_solva_spec(), system_size=4)
    assert report["ok"] is True, report["errors"]
    assert report["order"] == ["resolve", "validate", "lint",
                               "dry_run", "cost"]
    # the cost block is an ESTIMATE and each entry says so
    assert report["cost"], report
    # the lint must have seen the real steps
    linted_steps = {f["step"] for f in report["lint"]}
    assert linted_steps <= {"packmol_solvate", "xtb_relax",
                            "explicit_shell"}


def test_the_lint_refuses_a_dft_method_in_solvator():
    """The xTB-only limit must be visible to the WB2 chain: a DFT method
    in the SOLVATOR step is a physics error, not a style preference."""
    from delfin.agent.workflow_builder import run_checks

    spec = _solva_spec()
    spec["steps"][2]["kwargs"]["method"] = "B3LYP"
    report = run_checks(spec, system_size=4)
    lint_rules = {(f["rule"], f["level"]) for f in report["lint"]
                  if f["step"] in ("explicit_shell", "xtb_solvator")}
    assert ("xtb-solvator-method-not-xtb", "error") in lint_rules, report["lint"]
    assert report["ok"] is False


def test_the_lint_flags_implicit_and_explicit_double_counting():
    """A DFT solvent model (CPCM) alongside the explicit SOLVATOR shell
    double-counts solvation — the chain must flag it."""
    from delfin.agent.workflow_builder import run_checks

    spec = _solva_spec()
    spec["steps"].append(
        {"step": "orca_sp", "label": "dft_energy",
         "kwargs": {"method": "B3LYP", "basis": "def2-TZVP",
                    "solvent": "water", "solvent_model": "CPCM"}})
    report = run_checks(spec, system_size=4)
    rules = {f["rule"] for f in report["lint"]}
    assert "explicit-and-implicit-solvation-double-count" in rules, report["lint"]
    assert report["ok"] is False


# ---------------------------------------------------------------------------
# The rule functions themselves (green now — the hook into
# delfin/tools/_science.py is a protected file and goes to the operator
# as .gate/sv1_sciencelint_hook.patch; the run_checks tests above stay
# red until that hook is built).
# ---------------------------------------------------------------------------

def test_rule_module_flags_dft_method_directly():
    from delfin.tools.adapters.xtb_solvator_realism_lint import (
        lint_solvator_step)

    findings = lint_solvator_step(
        "explicit_shell", "xtb_solvator", {"method": "B3LYP"})
    assert len(findings) == 1
    assert findings[0].level == "error"
    assert findings[0].rule == "xtb-solvator-method-not-xtb"


def test_rule_module_accepts_xtb_methods_and_other_steps():
    from delfin.tools.adapters.xtb_solvator_realism_lint import (
        lint_solvator_step)

    assert lint_solvator_step(
        "explicit_shell", "xtb_solvator", {"method": "GFN2-xTB"}) == []
    assert lint_solvator_step(
        "explicit_shell", "xtb_solvator",
        {"method": "{solvator_method}"}) == []
    # a DFT method on a NON-solvator step is not this rule's business
    assert lint_solvator_step(
        "dft_energy", "orca_sp", {"method": "B3LYP"}) == []


def test_rule_module_flags_double_counting_directly():
    from delfin.tools.adapters.xtb_solvator_realism_lint import (
        lint_solvation_chain)

    steps = [
        {"step": "packmol_solvate", "label": "packmol_box",
         "kwargs": {"solvent_xyz": "{solvent_xyz}"}},
        {"step": "xtb_solvator", "label": "explicit_shell",
         "kwargs": {"solvent": "water"}},
        {"step": "orca_sp", "label": "dft_energy",
         "kwargs": {"solvent": "water", "solvent_model": "CPCM"}},
    ]
    findings = lint_solvation_chain(steps)
    assert len(findings) == 1
    assert findings[0].rule == "explicit-and-implicit-solvation-double-count"
    # F2 (nacht-s19): xTB continuum models double-count too
    for model in ("ALPB", "GBSA", "SMD", "CPCM"):
        kw = dict(steps[2]["kwargs"], solvent_model=model)
        assert lint_solvation_chain(steps[:2] + [
            {"step": "orca_sp", "label": "dft_energy", "kwargs": kw}]), model
    # no double-count without a shell step before the implicit model
    # (a chain that only optimises and runs the DFT single point is the
    # ordinary implicit-solvation case, not a double count)
    no_shell = [
        {"step": "xtb_opt", "label": "relax", "kwargs": {}},
        {"step": "orca_sp", "label": "dft_energy",
         "kwargs": {"solvent": "water", "solvent_model": "CPCM"}},
    ]
    assert lint_solvation_chain(no_shell) == []


def test_shell_step_itself_never_fires_the_double_count_rule():
    """F3 (nacht-s19): the SOLVATOR shell's own implicit component comes
    from the `solvent` key (ORCA: explicit shell + ALPB background) —
    that is deliberate ORCA design, not a double count. The rule keys on
    `solvent_model` only and excludes shell steps outright."""
    from delfin.tools.adapters.xtb_solvator_realism_lint import (
        lint_solvation_chain)

    steps = [
        {"step": "packmol_solvate", "label": "packmol_box",
         "kwargs": {"solvent_xyz": "{solvent_xyz}"}},
        {"step": "xtb_solvator", "label": "explicit_shell",
         "kwargs": {"solvent": "water", "method": "GFN2-xTB",
                    "solvent_model": "ALPB"}},
    ]
    assert lint_solvation_chain(steps) == []
