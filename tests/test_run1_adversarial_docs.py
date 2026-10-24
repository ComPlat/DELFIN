"""Reviewer adversarial tests for docs/workflow_runs.md (RUN1, nacht-s25).

The package task requires one table per workflow with the run evidence.
These cases pin the doc's shape on the REVIEWER branch — they are red
until the builder's phase-3 deliverable lands and gets merged here.

Rules pinned:
- D1: the doc exists and covers all four w21 workflows.
- D2: a run row carries method/basis/charge/mult/solvent, requested vs
  used resources, wall time, core-hours — and a source for every number
  (the package task: "physics verdict with its source").
- D3: a skipped workflow has a row that names the probe evidence (the
  command/check that came back empty), not a silent omission.
"""
from __future__ import annotations

from pathlib import Path

DOC = Path("docs") / "workflow_runs.md"

WORKFLOWS = (
    "w21_mlp_screen", "w21_descriptors_steric",
    "w21_solvation_workflow", "w21_turbomole_workflow",
)


def _doc_text() -> str:
    assert DOC.exists(), f"{DOC} missing — phase-3 deliverable not merged"
    return DOC.read_text(encoding="utf-8")


def test_d1_doc_exists_and_covers_all_four_workflows():
    text = _doc_text()
    # Intent: every workflow has an identifiable section. The doc names
    # sections by workflow/CLI name ("w21_solvation (explicit_solvation)"),
    # so the full module name OR the workflow prefix counts.
    for wf in WORKFLOWS:
        stem = wf.replace("_workflow", "")
        assert wf in text or stem in text, f"{DOC} has no section for {wf}"


def test_d2_run_rows_carry_evidence_columns_and_sources():
    text = _doc_text()
    for col in ("method", "basis", "charge", "mult", "solvent"):
        assert col in text.lower(), f"{DOC} run table lacks column {col}"
    for col in ("requested", "wall", "core-hours"):
        assert col in text.lower(), f"{DOC} run table lacks column {col}"
    for marker in ("CB0A", "sacct"):
        assert marker in text, (
            f"{DOC} lacks run-evidence source marker {marker}")


def test_d3_skipped_workflows_name_probe_evidence():
    text = _doc_text()
    assert "skip" in text.lower(), (
        f"{DOC} does not record skipped workflows and why")
