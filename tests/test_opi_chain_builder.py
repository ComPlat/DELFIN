"""The OPI chain is offered as a first-class workflow: a request like
"r2SCAN-3c opt+freq, then DLPNO-CCSD(T) single point" must draft a valid
pipeline (the ``opi_opt_freq_dlpno`` template), and every OPI ORCA step
must appear with a reason in the builder's output.

Red control: on the previous commit the request fell through to
``classic_opt_freq`` (no opi_opt_freq_dlpno template) and the OPI steps had
no catalog reasons.
"""

from __future__ import annotations

from delfin.agent import workflow_builder as wb

REQUEST = ("Run r2SCAN-3c opt+freq, then DLPNO-CCSD(T) single point, "
           "on ethanol (SMILES CCO).")


def test_request_maps_to_the_opi_chain_template():
    result = wb.build_workflow(REQUEST)
    assert result["source"] == "template:opi_opt_freq_dlpno"
    steps = [s["step"] for s in result["spec"]["steps"]]
    assert steps == ["smiles_to_xyz", "orca_opt_opi", "orca_freq_opi",
                     "orca_dlpno_opi"]


def test_opi_chain_spec_passes_all_checks():
    result = wb.build_workflow(REQUEST)
    assert result["checks"]["ok"], result["checks"]["errors"]
    order = result["checks"]["order"]
    assert order == ["resolve", "validate", "lint", "dry_run", "cost"]
    assert result["cost"], "a cost estimate is mandatory before submit"


def test_opi_steps_carry_catalog_reasons():
    result = wb.build_workflow(REQUEST)
    reasons = {e["step"]: e["reason"] for e in result["steps_reasons"]}
    for step in ("orca_opt_opi", "orca_freq_opi", "orca_dlpno_opi"):
        assert reasons.get(step, "").strip(), f"no reason for {step}"


def test_the_wft_step_does_not_inherit_the_dft_method():
    result = wb.build_workflow(REQUEST)
    for s in result["spec"]["steps"]:
        if s["step"] == "orca_dlpno_opi":
            assert s["kwargs"].get("method", "") == "", (
                "DLPNO step must not inherit the DFT functional broadcast")
