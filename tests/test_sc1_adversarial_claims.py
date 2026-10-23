"""Adversarial tests for SC1: delfin/agent/result_claims.check.

Red control: every test here must FAIL on the current tree (no
result_claims module yet). They measure the DELIVERABLE, not existing
code. Written by the reviewer (nacht-s4) before the builder's code
landed; nobody may weaken them to get green.

Focus: numbers stated in an answer must trace to evidence (CB0A-shaped
records), missing provenance flags, estimates are exempt, and the
function never raises.
"""

from __future__ import annotations

import pytest


def _ev(ref="benzene.out", value=-76.34, unit="Eh", method="B3LYP",
        basis="def2-SVP", solvent=None):
    rec = {"kind": "calc", "ref": ref, "value": value, "unit": unit,
           "method": method, "basis": basis}
    if solvent is not None:
        rec["solvent"] = solvent
    return rec


def test_module_and_check_exist():
    from delfin.agent import result_claims
    assert callable(result_claims.check)


def test_made_up_energy_is_flagged():
    from delfin.agent import result_claims
    answer = ("The electronic energy is -76.3421 Eh at B3LYP/def2-SVP.")
    findings = result_claims.check(answer, [])
    assert findings, "a number with no evidence must be flagged"
    assert any("76.3421" in str(f) or "76.3421" in str(getattr(f, "what", ""))
               for f in findings), findings


def test_estimate_marker_is_not_flagged():
    from delfin.agent import result_claims
    answer = ("Rough estimate: the energy is about -76.3 Eh "
              "(order-of-magnitude estimate, not from a run).")
    findings = result_claims.check(answer, [])
    assert not any("estimate" not in str(f).lower() and
                   "76.3" in str(f) for f in findings), findings


def test_energy_without_unit_is_flagged():
    from delfin.agent import result_claims
    answer = "The energy converged to -76.3421, so the job is fine."
    findings = result_claims.check(answer, [])
    assert any("unit" in str(f).lower() for f in findings), findings


def test_energy_without_method_basis_is_flagged():
    from delfin.agent import result_claims
    answer = "The final energy is -76.3421 Eh."
    findings = result_claims.check(answer, [])
    assert any(("method" in str(f).lower() or "basis" in str(f).lower())
               for f in findings), findings


def test_number_in_path_is_not_a_claim():
    from delfin.agent import result_claims
    answer = "See benzene_76_34.out for details."
    findings = result_claims.check(answer, [])
    assert not any("-76.34" in str(f) or "76_34" in str(f)
                   for f in findings), findings


def test_never_raises_on_weird_input():
    from delfin.agent import result_claims
    for text, ev in ((None, []), ("", None), ("x" * 200000, None),
                     ("１.２３４ Eh", "not-a-list"),
                     (object(), object())):
        try:
            result_claims.check(text, ev)
        except Exception as exc:  # noqa: BLE001
            pytest.fail(f"check({text!r:.30}, ...) raised {exc!r}")
