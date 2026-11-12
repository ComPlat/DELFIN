"""Review finding s25 on 477cdc3c (non-blocking): GlueGrader._claim_layer
fails the verdict whenever ANY claim exists — an ANSWERED claim whose
value has a source in the evidence also fails. result_claims.check
separates them: claim-no-source is level "error", a flagged-but-
sourced value is level "warn" (result_claims.py:170-189).

Red controls (before the fix):
- answered sourced number keeps the rubric verdict;
- a number with NO source still fails (direction check: the layer must
  keep firing).
"""
import types

from delfin.agent.bench_holdout_glue import GlueGrader


def _rubric_pass_result(answer, evidence_texts):
    # _rubric_verdict reads the run_pilot fields of the result; give it
    # a passing shape and put the answer/evidence on the attributes the
    # grader reads (result.answer, result.evidence).
    steps = [types.SimpleNamespace(text=t) for t in evidence_texts]
    return types.SimpleNamespace(
        success=True, quality_0_100=80,
        answer=answer, evidence=steps, unmeasured=False,
    )


def test_answered_sourced_number_keeps_verdict():
    g = GlueGrader()
    result = _rubric_pass_result(
        "PBE0/def2-TZVP energy is -15.3 Eh as reported in the output.",
        ["final single point energy -15.30 Eh, SCF converged"],
    )
    success, quality = g.grade(object(), result)
    assert success is True
    assert quality > 0


def test_number_without_source_still_fails():
    g = GlueGrader()
    result = _rubric_pass_result(
        "The energy is -15.3 Eh (invented).",
        ["some unrelated output text without the number"],
    )
    success, quality = g.grade(object(), result)
    assert success is False
    assert quality == 0
