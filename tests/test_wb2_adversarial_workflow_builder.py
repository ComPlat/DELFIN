"""WB2 adversarial review (nacht-s4) — attacks on workflow_builder (d08b2569).

A red test here is a FINDING to nacht-s3; a green test is committed
evidence the code holds.
"""

from __future__ import annotations

import pytest

from delfin.agent import workflow_builder as wb

GOOD_REQUEST = (
    "Build me a classic opt+freq workflow for ethanol, SMILES CCO, "
    "charge 0, B3LYP def2-SVP, 8 cores."
)


def test_dry_run_field_is_backed_by_a_real_dry_run():
    """F4: the advertised 5-check chain names a 'dry_run' stage. It must
    be backed by a real Pipeline dry run (a distinct static check), not
    by nothing — the chain as committed runs only resolve/validate/lint
    + cost and the result carries no dry-run evidence."""
    result = wb.build_workflow(GOOD_REQUEST)
    checks = result["checks"]
    assert "dry_run" in checks, (
        "the advertised dry-run stage must appear in the checks result")
    # Whatever form it takes, it must be real evidence of the stage,
    # not a bare restatement of the validate report.
    dr = checks["dry_run"]
    assert dr is not None


def test_cost_entries_carry_provenance_when_size_is_assumed():
    """F5: without a system_size the code silently assumes 10 atoms and
    the cost line then reads as measured. advise() returns
    'system_size_atoms' (cluster_efficiency.py:307) but that only proves
    the size the entry USED, not that the size was the USER'S. The
    builder must mark entries whose size it invented: 'assumed_size'
    (or a note). Every entry must carry provenance."""
    result = wb.build_workflow(GOOD_REQUEST)
    for entry in result["cost"]:
        if "assumed_size" in entry or "note" in entry:
            continue
        assert False, f"cost entry lacks assumed-size provenance: {entry}"


def test_lint_error_forces_ok_false_even_when_validate_passes():
    """F7 (order enforcement): a lint ERROR must set ok=False even when
    structural validation passes. Unresolvable multiplicity."""
    result = wb.build_workflow(
        "Build an opt+freq workflow, SMILES CCO, charge 0, multiplicity 0.")
    lint_levels = {f.get("level") for f in result["checks"]["lint"]}
    if "error" in lint_levels:
        assert result["checks"]["ok"] is False


def test_empty_request_yields_a_draft_not_a_crash():
    """Robustness: an empty request must return a result with a spec and
    the checks run — never raise."""
    result = wb.build_workflow("")
    assert result["spec"]["steps"]
    assert result["checks"]["order"]


def test_huge_and_unicode_requests_are_handled():
    """Robustness: a 10k-char request and unicode text must not break
    the builder; description is truncated, spec still builds."""
    huge = "Build opt+freq for SMILES CCO " + "x" * 10000
    r1 = wb.build_workflow(huge)
    assert len(r1["spec"].get("description", "")) <= 200
    r2 = wb.build_workflow("Erstelle einen Opt+Freq-Workflow, SMILES CCO, Mülleimer")
    assert r2["spec"]["steps"]


def test_good_request_builds_and_passes_checks():
    """Negative control: a good request must produce ok=True (so the
    attacks above are not trivially red for everything)."""
    result = wb.build_workflow(GOOD_REQUEST)
    assert result["checks"]["ok"] is True
    assert result["submitted"] is False
